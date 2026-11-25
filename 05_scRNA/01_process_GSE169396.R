# 用途：GSE169396单细胞参考预处理。
# 输入：GSE169396_RAW.tar。
# 输出：Seurat对象、细胞注释和质控表。

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

set.seed(2025)

project <- Sys.getenv("CFRNA_GMF_ROOT", unset = ".")
geo_dir <- file.path(project, "GEO_singlecell_raw_downloads_20260525",
                     "11_human_femoral_head_scRNA_GSE169396")
work_dir <- file.path(project, "CellFreeGMF_work",
                      "09_Set7_GSE169396_human_femoral_head_GMF")
out_dir <- file.path(work_dir, "GSE169396_human_femoral_head_reference")
raw10x_dir <- file.path(out_dir, "10x")
dir.create(raw10x_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

raw_tar <- file.path(geo_dir, "raw_downloads", "GSE169396_RAW.tar")
if (!file.exists(raw_tar)) {
  stop("Missing GSE169396_RAW.tar: ", raw_tar)
}
untar(raw_tar, exdir = raw10x_dir)

sample_info <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4"),
  gsm = c("GSM5201883", "GSM5201884", "GSM5201885", "GSM5201886"),
  sex = c("female", "female", "male", "male"),
  age = c(61, 45, 66, 31),
  tissue = "femoral_head",
  stringsAsFactors = FALSE
)

read_one_sample <- function(sample_row) {
  prefix <- paste0(sample_row$gsm, "_", sample_row$sample_id)
  barcodes <- file.path(raw10x_dir, paste0(prefix, "_barcodes.tsv.gz"))
  features <- file.path(raw10x_dir, paste0(prefix, "_features.tsv.gz"))
  matrix <- file.path(raw10x_dir, paste0(prefix, "_matrix.mtx.gz"))
  if (!all(file.exists(c(barcodes, features, matrix)))) {
    stop("Missing one or more 10x files for ", prefix)
  }

  counts <- ReadMtx(
    mtx = matrix,
    cells = barcodes,
    features = features,
    feature.column = 2,
    unique.features = TRUE
  )
  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE169396_FHTC",
    min.cells = 0,
    min.features = 0
  )
  obj$sample_id <- sample_row$sample_id
  obj$gsm <- sample_row$gsm
  obj$sex <- sample_row$sex
  obj$age <- sample_row$age
  obj$tissue <- sample_row$tissue
  obj$barcode_raw <- colnames(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  low <- as.numeric(quantile(obj$nFeature_RNA, probs = 0.02, na.rm = TRUE))
  high <- as.numeric(quantile(obj$nFeature_RNA, probs = 0.98, na.rm = TRUE))
  obj$qc_keep_original_policy <- obj$nFeature_RNA > low &
    obj$nFeature_RNA < high &
    obj$percent.mt <= 10
  obj$qc_low_nFeature_2pct <- low
  obj$qc_high_nFeature_98pct <- high

  obj <- subset(obj, subset = qc_keep_original_policy)
  obj <- RenameCells(obj, new.names = paste(sample_row$sample_id, obj$barcode_raw, sep = "_"))
  obj
}

objs <- lapply(seq_len(nrow(sample_info)), function(i) read_one_sample(sample_info[i, ]))
names(objs) <- sample_info$sample_id
obj <- Reduce(function(x, y) merge(x, y), objs)
obj <- JoinLayers(obj)

qc_by_sample <- obj@meta.data %>%
  count(sample_id, sex, age, tissue, name = "n_cells_after_qc") %>%
  left_join(
    bind_rows(lapply(objs, function(x) {
      data.frame(
        sample_id = unique(x$sample_id),
        qc_low_nFeature_2pct = unique(x$qc_low_nFeature_2pct),
        qc_high_nFeature_98pct = unique(x$qc_high_nFeature_98pct),
        median_nCount_RNA = median(x$nCount_RNA),
        median_nFeature_RNA = median(x$nFeature_RNA),
        median_percent_mt = median(x$percent.mt)
      )
    })),
    by = "sample_id"
  )
write.csv(qc_by_sample, file.path(out_dir, "tables", "GSE169396_QC_by_sample.csv"), row.names = FALSE)

qc_summary <- data.frame(
  item = c(
    "input_source",
    "qc_policy",
    "n_cells_after_qc",
    "n_genes",
    "median_nCount_RNA",
    "median_nFeature_RNA",
    "median_percent_mt"
  ),
  value = c(
    "GEO 10x matrix files from GSE169396_RAW.tar",
    "Original-paper-style QC: remove cells in bottom/top 2 percent of nFeature_RNA within each sample and percent.mt > 10",
    ncol(obj),
    nrow(obj),
    median(obj$nCount_RNA),
    median(obj$nFeature_RNA),
    median(obj$percent.mt)
  )
)
write.csv(qc_summary, file.path(out_dir, "tables", "GSE169396_QC_summary.csv"), row.names = FALSE)

write.csv(
  obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id, sample_id, gsm, sex, age, tissue, barcode_raw,
           nCount_RNA, nFeature_RNA, percent.mt,
           qc_low_nFeature_2pct, qc_high_nFeature_98pct),
  file.path(out_dir, "tables", "GSE169396_cell_QC_metrics.csv"),
  row.names = FALSE
)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)

res_grid <- seq(0.1, 2.0, by = 0.1)
grid_rows <- list()
for (res in res_grid) {
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  col <- paste0("RNA_snn_res.", res)
  grid_rows[[as.character(res)]] <- data.frame(
    resolution = res,
    n_clusters = length(unique(obj@meta.data[[col]]))
  )
}
resolution_grid <- bind_rows(grid_rows)
resolution_grid$distance_to_16 <- abs(resolution_grid$n_clusters - 16)
resolution_grid$distance_to_0p5 <- abs(resolution_grid$resolution - 0.5)
exact_16 <- resolution_grid %>% filter(n_clusters == 16)
if (nrow(exact_16) > 0) {
  chosen <- exact_16 %>% arrange(distance_to_0p5) %>% slice(1)
  resolution_policy <- "auto_exact_16_clusters_then_closest_to_0.5"
} else {
  chosen <- resolution_grid %>% arrange(distance_to_16, distance_to_0p5) %>% slice(1)
  resolution_policy <- "auto_nearest_to_original_16_clusters_no_exact_16_found"
}
chosen_res <- chosen$resolution[[1]]
chosen_col <- paste0("RNA_snn_res.", chosen_res)
obj$cluster_reconstructed <- as.character(obj@meta.data[[chosen_col]])
write.csv(
  resolution_grid,
  file.path(out_dir, "tables", "GSE169396_resolution_grid.csv"),
  row.names = FALSE
)

obj <- RunTSNE(obj, dims = 1:20, reduction = "pca", seed.use = 2025, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca", seed.use = 2025, verbose = FALSE)

marker_sets <- list(
  Granulocytes = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "CEACAM8", "MPO", "ELANE"),
  T_cells = c("CD3D", "CD3E", "CD3G", "TRAC", "CD2", "IL7R", "CD8A", "CD8B"),
  Monocytes = c("CD14", "LYZ", "LST1", "FCN1", "S100A8", "S100A9", "MS4A7", "CTSS"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD74", "HLA-DRA"),
  Red_blood_cells = c("HBA1", "HBA2", "HBB", "HBD", "GYPA", "ALAS2"),
  Osteoblastic_lineage_cells = c("COL1A1", "COL1A2", "SPARC", "BGLAP", "RUNX2", "ALPL", "IBSP", "SPP1", "PRSS23", "MXRA8"),
  Endothelial_cells = c("PECAM1", "VWF", "KDR", "CLDN5", "EMCN", "RAMP2", "ESAM"),
  Endothelial_progenitor_cells = c("KIT", "PROM1", "CD34", "FLT1", "KDR", "ENG"),
  Plasmacytoid_dendritic_cells = c("GZMB", "IL3RA", "CLEC4C", "LILRA4", "TCF4", "IRF7")
)

marker_df <- bind_rows(lapply(names(marker_sets), function(nm) {
  data.frame(cell_type_marker_set = nm, gene = marker_sets[[nm]])
}))
write.csv(marker_df, file.path(out_dir, "tables", "GSE169396_marker_sets_used.csv"), row.names = FALSE)

Idents(obj) <- "cluster_reconstructed"
norm_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
clusters <- sort(unique(obj$cluster_reconstructed))

score_rows <- list()
for (cl in clusters) {
  cells <- WhichCells(obj, idents = cl)
  for (set_name in names(marker_sets)) {
    genes <- intersect(marker_sets[[set_name]], rownames(norm_mat))
    score <- if (length(genes) > 0) {
      mean(Matrix::colMeans(norm_mat[genes, cells, drop = FALSE]))
    } else {
      NA_real_
    }
    score_rows[[length(score_rows) + 1]] <- data.frame(
      cluster = cl,
      cell_type_marker_set = set_name,
      n_markers_present = length(genes),
      mean_log_normalized_expr = score
    )
  }
}
score_df <- bind_rows(score_rows)
score_wide <- tidyr::pivot_wider(
  score_df,
  id_cols = cluster,
  names_from = cell_type_marker_set,
  values_from = mean_log_normalized_expr
)
score_mat <- as.matrix(score_wide[, -1, drop = FALSE])
rownames(score_mat) <- score_wide$cluster
score_z <- apply(score_mat, 2, function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    rep(0, length(x))
  } else {
    as.numeric(scale(x))
  }
})
score_z <- as.data.frame(score_z)
score_z$cluster <- rownames(score_mat)
score_z <- score_z[, c("cluster", setdiff(colnames(score_z), "cluster"))]

type_cols <- setdiff(colnames(score_z), "cluster")
score_numeric <- as.matrix(score_z[, type_cols, drop = FALSE])
best_idx <- max.col(score_numeric, ties.method = "first")
cluster_anno <- data.frame(
  cluster = score_z$cluster,
  FHTC_major_type_reconstructed = type_cols[best_idx],
  marker_score_z = score_numeric[cbind(seq_len(nrow(score_numeric)), best_idx)],
  stringsAsFactors = FALSE
)

cluster_sizes <- obj@meta.data %>%
  count(cluster_reconstructed, name = "n_cells") %>%
  rename(cluster = cluster_reconstructed)
cluster_anno <- cluster_sizes %>%
  left_join(cluster_anno, by = "cluster") %>%
  arrange(as.numeric(cluster)) %>%
  mutate(reconstructed_cluster_id = sprintf("C%02d", row_number()))

obj$FHTC_major_type_reconstructed <- cluster_anno$FHTC_major_type_reconstructed[
  match(obj$cluster_reconstructed, cluster_anno$cluster)
]
obj$reconstructed_cluster_id <- cluster_anno$reconstructed_cluster_id[
  match(obj$cluster_reconstructed, cluster_anno$cluster)
]

tsne <- Embeddings(obj, "tsne")
umap <- Embeddings(obj, "umap")
cell_anno <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    barcode = barcode_raw,
    tsne_1 = tsne[cell_id, 1],
    tsne_2 = tsne[cell_id, 2],
    umap_1 = umap[cell_id, 1],
    umap_2 = umap[cell_id, 2],
    annotation_status = "reconstructed marker-based annotation following original broad FHTC cell types"
  ) %>%
  select(
    cell_id,
    sample_id,
    gsm,
    barcode,
    sex,
    age,
    tissue,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    cluster_reconstructed,
    reconstructed_cluster_id,
    FHTC_major_type_reconstructed,
    annotation_status,
    tsne_1,
    tsne_2,
    umap_1,
    umap_2
  )

write.csv(cluster_anno, file.path(out_dir, "tables", "GSE169396_cluster_to_FHTC_reconstructed_annotation.csv"), row.names = FALSE)
write.csv(score_df, file.path(out_dir, "tables", "GSE169396_cluster_marker_scores_long.csv"), row.names = FALSE)
write.csv(score_z, file.path(out_dir, "tables", "GSE169396_cluster_marker_scores_z.csv"), row.names = FALSE)
write.csv(cell_anno, file.path(out_dir, "tables", "GSE169396_FHTC_reconstructed_cell_annotation.csv"), row.names = FALSE)

present_marker_genes <- intersect(unique(unlist(marker_sets)), rownames(obj))

png(file.path(out_dir, "figures", "GSE169396_tSNE_clusters.png"), width = 2200, height = 1800, res = 180)
print(DimPlot(obj, reduction = "tsne", group.by = "reconstructed_cluster_id", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed clusters"))
dev.off()

png(file.path(out_dir, "figures", "GSE169396_tSNE_FHTC_major_types.png"), width = 2400, height = 1800, res = 180)
print(DimPlot(obj, reduction = "tsne", group.by = "FHTC_major_type_reconstructed", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed FHTC major cell types"))
dev.off()

png(file.path(out_dir, "figures", "GSE169396_UMAP_FHTC_major_types.png"), width = 2400, height = 1800, res = 180)
print(DimPlot(obj, reduction = "umap", group.by = "FHTC_major_type_reconstructed", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed FHTC major cell types"))
dev.off()

if (length(present_marker_genes) > 0) {
  png(file.path(out_dir, "figures", "GSE169396_FHTC_marker_dotplot.png"), width = 3600, height = 1800, res = 180)
  print(DotPlot(obj, features = present_marker_genes, group.by = "FHTC_major_type_reconstructed") +
          RotatedAxis() +
          ggtitle("GSE169396 marker expression by reconstructed FHTC type"))
  dev.off()
}

saveRDS(obj, file.path(out_dir, "GSE169396_human_femoral_head_reconstructed_seurat.rds"))

run_summary <- data.frame(
  item = c(
    "chosen_resolution",
    "chosen_n_clusters",
    "resolution_policy",
    "n_cells_final",
    "n_genes",
    "n_FHTC_major_types",
    "sample_count",
    "annotation_status"
  ),
  value = c(
    chosen_res,
    chosen$n_clusters[[1]],
    resolution_policy,
    ncol(obj),
    nrow(obj),
    length(unique(obj$FHTC_major_type_reconstructed)),
    length(unique(obj$sample_id)),
    "reconstructed marker-based annotation; GEO does not provide barcode-level original annotation"
  )
)
write.csv(run_summary, file.path(out_dir, "tables", "GSE169396_reconstructed_reference_summary.csv"), row.names = FALSE)

if (tolower(Sys.getenv("SKIP_FULL_MARKERS", unset = "false")) %in% c("1", "true", "yes")) {
  cat("GSE169396 human femoral head reconstruction finished before full FindAllMarkers.\n")
  cat("Full Wilcoxon marker table skipped because SKIP_FULL_MARKERS=true.\n")
  cat("Output:", out_dir, "\n")
  print(run_summary)
  quit(save = "no", status = 0)
}

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  verbose = FALSE
)
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

tsne <- Embeddings(obj, "tsne")
umap <- Embeddings(obj, "umap")
cell_anno <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    barcode = barcode_raw,
    tsne_1 = tsne[cell_id, 1],
    tsne_2 = tsne[cell_id, 2],
    umap_1 = umap[cell_id, 1],
    umap_2 = umap[cell_id, 2],
    annotation_status = "reconstructed marker-based annotation following original broad FHTC cell types"
  ) %>%
  select(
    cell_id,
    sample_id,
    gsm,
    barcode,
    sex,
    age,
    tissue,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    cluster_reconstructed,
    reconstructed_cluster_id,
    FHTC_major_type_reconstructed,
    annotation_status,
    tsne_1,
    tsne_2,
    umap_1,
    umap_2
  )

write.csv(cluster_anno, file.path(out_dir, "tables", "GSE169396_cluster_to_FHTC_reconstructed_annotation.csv"), row.names = FALSE)
write.csv(score_df, file.path(out_dir, "tables", "GSE169396_cluster_marker_scores_long.csv"), row.names = FALSE)
write.csv(score_z, file.path(out_dir, "tables", "GSE169396_cluster_marker_scores_z.csv"), row.names = FALSE)
write.csv(top50, file.path(out_dir, "tables", "GSE169396_cluster_top50_markers.csv"), row.names = FALSE)
write.csv(cell_anno, file.path(out_dir, "tables", "GSE169396_FHTC_reconstructed_cell_annotation.csv"), row.names = FALSE)

present_marker_genes <- intersect(unique(unlist(marker_sets)), rownames(obj))

png(file.path(out_dir, "figures", "GSE169396_tSNE_clusters.png"), width = 2200, height = 1800, res = 180)
print(DimPlot(obj, reduction = "tsne", group.by = "reconstructed_cluster_id", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed clusters"))
dev.off()

png(file.path(out_dir, "figures", "GSE169396_tSNE_FHTC_major_types.png"), width = 2400, height = 1800, res = 180)
print(DimPlot(obj, reduction = "tsne", group.by = "FHTC_major_type_reconstructed", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed FHTC major cell types"))
dev.off()

png(file.path(out_dir, "figures", "GSE169396_UMAP_FHTC_major_types.png"), width = 2400, height = 1800, res = 180)
print(DimPlot(obj, reduction = "umap", group.by = "FHTC_major_type_reconstructed", label = TRUE, repel = TRUE) +
        ggtitle("GSE169396 reconstructed FHTC major cell types"))
dev.off()

if (length(present_marker_genes) > 0) {
  png(file.path(out_dir, "figures", "GSE169396_FHTC_marker_dotplot.png"), width = 3600, height = 1800, res = 180)
  print(DotPlot(obj, features = present_marker_genes, group.by = "FHTC_major_type_reconstructed") +
          RotatedAxis() +
          ggtitle("GSE169396 marker expression by reconstructed FHTC type"))
  dev.off()
}

saveRDS(obj, file.path(out_dir, "GSE169396_human_femoral_head_reconstructed_seurat.rds"))

run_summary <- data.frame(
  item = c(
    "chosen_resolution",
    "chosen_n_clusters",
    "resolution_policy",
    "n_cells_final",
    "n_genes",
    "n_FHTC_major_types",
    "sample_count",
    "annotation_status"
  ),
  value = c(
    chosen_res,
    chosen$n_clusters[[1]],
    resolution_policy,
    ncol(obj),
    nrow(obj),
    length(unique(obj$FHTC_major_type_reconstructed)),
    length(unique(obj$sample_id)),
    "reconstructed marker-based annotation; GEO does not provide barcode-level original annotation"
  )
)
write.csv(run_summary, file.path(out_dir, "tables", "GSE169396_reconstructed_reference_summary.csv"), row.names = FALSE)

cat("GSE169396 human femoral head reconstruction finished.\n")
cat("Output:", out_dir, "\n")
print(run_summary)
