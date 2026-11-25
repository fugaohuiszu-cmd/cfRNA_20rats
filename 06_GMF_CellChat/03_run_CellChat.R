#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(reticulate)
  library(CellChat)
  library(dplyr)
  library(readr)
  library(future)
})

set.seed(2025)
use_python(Sys.getenv("CELLFREEGMF_PYTHON", unset = "python"), required = TRUE)

work_dir <- file.path(".", "03_results", "06_GMF_CellChat")
gmf_dir <- file.path(work_dir, "output", "gmf")
out_dir <- file.path(work_dir, "output", "cellchat")
table_dir <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

npz_paths <- c(
  N = file.path(gmf_dir, "N", "gnmf_results.npz"),
  OP = file.path(gmf_dir, "OP", "gnmf_results.npz")
)
missing_npz <- npz_paths[!file.exists(npz_paths)]
if (length(missing_npz) > 0) {
  stop("Missing GMF result: ", paste(missing_npz, collapse = "; "))
}

cell_levels <- c(
  "Neutrophils",
  "T/NK cells",
  "Myeloid cells",
  "B cells",
  "Pro-B cells",
  "Plasma cells",
  "Erythroid cells",
  "Osteogenic cells",
  "Endothelial cells",
  "Mast cells",
  "pDCs"
)

np <- import("numpy")

read_npz <- function(path) {
  z <- np$load(path, allow_pickle = TRUE)
  available <- py_to_r(z$files)
  required <- c("Result_cfRNA_cell", "cfRNA_name", "cell_id", "cell_type")
  missing <- setdiff(required, available)
  if (length(missing) > 0) {
    stop("Missing NPZ fields in ", path, ": ", paste(missing, collapse = ", "))
  }
  list(
    matrix = py_to_r(z$f[["Result_cfRNA_cell"]]),
    genes = as.character(py_to_r(z$f[["cfRNA_name"]])),
    cell_ids = as.character(py_to_r(z$f[["cell_id"]])),
    cell_types = as.character(py_to_r(z$f[["cell_type"]]))
  )
}

gmf <- lapply(npz_paths, read_npz)

if (!identical(gmf$N$genes, gmf$OP$genes)) {
  stop("N and OP gene orders are not identical")
}
if (!identical(gmf$N$cell_ids, gmf$OP$cell_ids)) {
  stop("N and OP cell orders are not identical")
}
if (!identical(gmf$N$cell_types, gmf$OP$cell_types)) {
  stop("N and OP cell-type annotations are not identical")
}
if (anyDuplicated(gmf$N$genes) > 0) {
  stop("Duplicated human gene symbols in GMF output")
}
if (!setequal(unique(gmf$N$cell_types), cell_levels)) {
  stop("Cell types in GMF output do not match the predefined 11 cell types")
}

prepare_group <- function(x, group_name) {
  mat <- as.matrix(x$matrix)
  if (nrow(mat) != length(x$genes) || ncol(mat) != length(x$cell_ids)) {
    stop(group_name, " matrix dimensions do not match gene and cell metadata")
  }
  if (any(!is.finite(mat))) {
    stop(group_name, " matrix contains non-finite values")
  }
  if (any(mat < 0)) {
    stop(group_name, " matrix contains negative values")
  }

  rownames(mat) <- x$genes
  colnames(mat) <- paste0(x$cell_ids, "_", group_name)
  meta <- data.frame(
    cell_id = colnames(mat),
    cell_type = factor(x$cell_types, levels = cell_levels),
    group = group_name,
    samples = factor(group_name, levels = c("N", "OP")),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$cell_id
  list(matrix = mat, meta = meta)
}

prepared <- list(
  N = prepare_group(gmf$N, "N"),
  OP = prepare_group(gmf$OP, "OP")
)

cell_counts <- as.data.frame(table(prepared$N$meta$cell_type), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("cell_type", "cell_count")
write_csv(cell_counts, file.path(table_dir, "cell_counts.csv"))

run_cellchat <- function(mat, meta, group_name) {
  database <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
  database_genes <- extractGene(list(
    interaction = database$interaction,
    complex = CellChatDB.human$complex,
    cofactor = CellChatDB.human$cofactor,
    geneInfo = CellChatDB.human$geneInfo
  ))
  keep_genes <- intersect(rownames(mat), database_genes)
  if (length(keep_genes) == 0) {
    stop(group_name, " has no genes overlapping CellChatDB.human")
  }

  input_summary <- data.frame(
    group = group_name,
    input_genes = nrow(mat),
    CellChatDB_genes = length(keep_genes),
    cells = ncol(mat),
    cell_types = nlevels(meta$cell_type),
    minimum_cells_per_type = min(table(meta$cell_type)),
    stringsAsFactors = FALSE
  )
  write_csv(input_summary, file.path(table_dir, paste0("input_summary_", group_name, ".csv")))

  mat <- mat[keep_genes, , drop = FALSE]
  cellchat <- createCellChat(object = mat, meta = meta, group.by = "cell_type")
  cellchat@DB <- database
  cellchat <- subsetData(cellchat, features = rownames(mat))
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  netAnalysis_computeCentrality(cellchat, slot.name = "netP")
}

extract_edges <- function(cellchat, group_name) {
  weight <- as.data.frame(as.table(cellchat@net$weight), stringsAsFactors = FALSE)
  count <- as.data.frame(as.table(cellchat@net$count), stringsAsFactors = FALSE)
  colnames(weight) <- c("source", "target", "interaction_weight")
  colnames(count) <- c("source", "target", "interaction_count")
  left_join(weight, count, by = c("source", "target")) |>
    mutate(group = group_name, .before = 1) |>
    arrange(desc(interaction_weight))
}

options(future.globals.maxSize = 128 * 1024^3)
future::plan(future::sequential)

cellchat_N <- run_cellchat(prepared$N$matrix, prepared$N$meta, "N")
saveRDS(cellchat_N, file.path(object_dir, "cellchat_N.rds"), compress = "xz")

cellchat_OP <- run_cellchat(prepared$OP$matrix, prepared$OP$meta, "OP")
saveRDS(cellchat_OP, file.path(object_dir, "cellchat_OP.rds"), compress = "xz")

cellchat_merged <- mergeCellChat(
  list(N = cellchat_N, OP = cellchat_OP),
  cell.prefix = TRUE,
  add.names = c("N", "OP")
)
saveRDS(cellchat_merged, file.path(object_dir, "cellchat_merged.rds"), compress = "xz")

communication_N <- subsetCommunication(cellchat_N) |>
  mutate(group = "N", .before = 1)
communication_OP <- subsetCommunication(cellchat_OP) |>
  mutate(group = "OP", .before = 1)
write_csv(communication_N, file.path(table_dir, "communications_N.csv"))
write_csv(communication_OP, file.path(table_dir, "communications_OP.csv"))

edges_N <- extract_edges(cellchat_N, "N")
edges_OP <- extract_edges(cellchat_OP, "OP")
write_csv(edges_N, file.path(table_dir, "network_edges_N.csv"))
write_csv(edges_OP, file.path(table_dir, "network_edges_OP.csv"))

differential_edges <- full_join(
  edges_N |>
    select(source, target, N_weight = interaction_weight, N_count = interaction_count),
  edges_OP |>
    select(source, target, OP_weight = interaction_weight, OP_count = interaction_count),
  by = c("source", "target")
) |>
  mutate(
    across(c(N_weight, N_count, OP_weight, OP_count), ~ coalesce(.x, 0)),
    OP_minus_N_weight = OP_weight - N_weight,
    OP_minus_N_count = OP_count - N_count
  ) |>
  arrange(desc(abs(OP_minus_N_weight)))
write_csv(differential_edges, file.path(table_dir, "differential_network_edges.csv"))

network_summary <- bind_rows(
  edges_N |>
    summarise(
      group = "N",
      total_weight = sum(interaction_weight),
      total_count = sum(interaction_count),
      nonzero_edges = sum(interaction_weight > 0)
    ),
  edges_OP |>
    summarise(
      group = "OP",
      total_weight = sum(interaction_weight),
      total_count = sum(interaction_count),
      nonzero_edges = sum(interaction_weight > 0)
    )
) |>
  mutate(
    relative_weight = total_weight / total_weight[group == "N"],
    relative_count = total_count / total_count[group == "N"]
  )
write_csv(network_summary, file.path(table_dir, "group_interaction_values.csv"))

lr_keys <- c("source", "target", "interaction_name", "interaction_name_2", "pathway_name", "ligand", "receptor")
lr_N <- communication_N |>
  select(all_of(lr_keys), N_probability = prob, N_pvalue = pval)
lr_OP <- communication_OP |>
  select(all_of(lr_keys), OP_probability = prob, OP_pvalue = pval)
differential_lr <- full_join(lr_N, lr_OP, by = lr_keys) |>
  mutate(
    N_probability = coalesce(N_probability, 0),
    OP_probability = coalesce(OP_probability, 0),
    OP_minus_N_probability = OP_probability - N_probability
  ) |>
  arrange(desc(abs(OP_minus_N_probability)))
write_csv(differential_lr, file.path(table_dir, "differential_ligand_receptor.csv"))

package_versions <- data.frame(
  package = c("R", "CellChat", "reticulate", "dplyr", "readr"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("CellChat")),
    as.character(packageVersion("reticulate")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("readr"))
  )
)
write_csv(package_versions, file.path(table_dir, "package_versions.csv"))
