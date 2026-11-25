#!/usr/bin/env Rscript
# 用途：11类参考细胞注释及绘图。
# 输入：参考Seurat对象和人工整理注释；输出：11类细胞注释、UMAP和标记基因图。

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

project <- Sys.getenv("CFRNA_GMF_ROOT", unset = ".")
ref_dir <- file.path(
  project,
  "CellFreeGMF_work",
  "09_Set7_GSE169396_human_femoral_head_GMF",
  "GSE169396_human_femoral_head_reference"
)

rds_path <- file.path(ref_dir, "GSE169396_human_femoral_head_reconstructed_seurat.rds")
curated_dir <- file.path(ref_dir, "curated_annotation")
fig_dir <- file.path(ref_dir, "figures", "curated_res0p3", "proB_split_11_major")
dir.create(curated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cell_path <- file.path(curated_dir, "GSE169396_res0p3_curated_cell_metadata_for_GMF.csv")
cluster_path <- file.path(curated_dir, "GSE169396_res0p3_curated_cluster_annotation_for_GMF.csv")

if (!file.exists(rds_path)) stop("Missing Seurat RDS: ", rds_path)
if (!file.exists(cell_path)) stop("Missing curated cell metadata: ", cell_path)
if (!file.exists(cluster_path)) stop("Missing curated cluster annotation: ", cluster_path)

obj <- readRDS(rds_path)
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
cell_meta <- read.csv(cell_path, stringsAsFactors = FALSE, check.names = FALSE)
cluster_anno <- read.csv(cluster_path, stringsAsFactors = FALSE, check.names = FALSE)

required_cell_cols <- c("cell_id", "raw_cluster", "fine_annotation", "GMF_major_celltype")
required_cluster_cols <- c("raw_cluster", "fine_annotation", "GMF_major_celltype", "confidence", "top_marker_evidence")
missing_cell_cols <- setdiff(required_cell_cols, colnames(cell_meta))
missing_cluster_cols <- setdiff(required_cluster_cols, colnames(cluster_anno))
if (length(missing_cell_cols) > 0) {
  stop("Curated cell metadata missing columns: ", paste(missing_cell_cols, collapse = ", "))
}
if (length(missing_cluster_cols) > 0) {
  stop("Curated cluster annotation missing columns: ", paste(missing_cluster_cols, collapse = ", "))
}
if (!"umap" %in% Reductions(obj)) stop("Seurat object has no UMAP reduction")
if (!all(colnames(obj) %in% cell_meta$cell_id)) {
  missing_cells <- setdiff(colnames(obj), cell_meta$cell_id)
  stop("Curated cell metadata missing Seurat cells, n=", length(missing_cells))
}

cell_meta <- cell_meta[match(colnames(obj), cell_meta$cell_id), ]
if (!identical(cell_meta$cell_id, colnames(obj))) {
  stop("Cell metadata order does not match Seurat cells")
}

fine_to_11_major <- c(
  "Mature neutrophils / granulocytes" = "Granulocytes / neutrophils",
  "Immature neutrophils / promyelocyte-like granulocytes" = "Granulocytes / neutrophils",
  "Activated inflammatory neutrophils / granulocytes" = "Granulocytes / neutrophils",
  "IFN-stimulated neutrophils / granulocytes" = "Granulocytes / neutrophils",
  "CD4-like / CCR7+ T cells" = "T cells / cytotoxic lymphocytes",
  "Cytotoxic T/NK-like cells" = "T cells / cytotoxic lymphocytes",
  "CD1C+ monocyte/DC-like myeloid cells" = "Monocyte/macrophage/DC-like myeloid",
  "C1QC+ tissue-resident-like macrophages" = "Monocyte/macrophage/DC-like myeloid",
  "RAG1/RAG2+ pre-B / immature B cells" = "B-lineage cells",
  "DNTT+/VPREB1+ pro-B / early B-lineage cells" = "Pro-B / early B-lineage cells",
  "Plasma cells" = "Plasma cells",
  "Erythroid precursors / erythroblasts" = "Erythroid cells",
  "Mature erythrocytes / RBC" = "Erythroid cells",
  "Osteoblastic lineage / skeletal stromal-chondro-osteogenic cells" = "Osteoblastic lineage cells",
  "Endothelial cells" = "Endothelial cells",
  "Mast cells" = "Mast cells",
  "Plasmacytoid dendritic cells / pDCs" = "pDCs"
)

major_11_levels <- c(
  "Granulocytes / neutrophils",
  "T cells / cytotoxic lymphocytes",
  "Monocyte/macrophage/DC-like myeloid",
  "B-lineage cells",
  "Pro-B / early B-lineage cells",
  "Plasma cells",
  "Erythroid cells",
  "Osteoblastic lineage cells",
  "Endothelial cells",
  "Mast cells",
  "pDCs"
)

celltype_colors <- c(
  "Neutrophils" = "#EE7B20",
  "T/NK cells" = "#8963A7",
  "Myeloid cells" = "#2171A3",
  "B cells" = "#6AA5C7",
  "Pro-B cells" = "#4B4FA3",
  "Plasma cells" = "#F9D6C2",
  "Erythroid cells" = "#C92928",
  "Osteogenic cells" = "#329838",
  "Endothelial cells" = "#88564A",
  "Mast cells" = "#F2B166",
  "pDCs" = "#BFF9D6"
)

celltype_short_labels <- c(
  "Granulocytes / neutrophils" = "Neutrophils",
  "T cells / cytotoxic lymphocytes" = "T/NK cells",
  "Monocyte/macrophage/DC-like myeloid" = "Myeloid cells",
  "B-lineage cells" = "B cells",
  "Pro-B / early B-lineage cells" = "Pro-B cells",
  "Plasma cells" = "Plasma cells",
  "Erythroid cells" = "Erythroid cells",
  "Osteoblastic lineage cells" = "Osteogenic cells",
  "Endothelial cells" = "Endothelial cells",
  "Mast cells" = "Mast cells",
  "pDCs" = "pDCs"
)
short_11_levels <- unname(celltype_short_labels[major_11_levels])
celltype_colors <- celltype_colors[short_11_levels]

missing_fine_cell <- setdiff(unique(cell_meta$fine_annotation), names(fine_to_11_major))
missing_fine_cluster <- setdiff(unique(cluster_anno$fine_annotation), names(fine_to_11_major))
if (length(missing_fine_cell) > 0 || length(missing_fine_cluster) > 0) {
  stop(
    "Fine annotation labels do not match mapping. Cell extra: ",
    paste(missing_fine_cell, collapse = "; "),
    " Cluster extra: ",
    paste(missing_fine_cluster, collapse = "; ")
  )
}

cell_meta$GMF_major_celltype_11_proB_split <- unname(fine_to_11_major[cell_meta$fine_annotation])
cluster_anno$GMF_major_celltype_11_proB_split <- unname(fine_to_11_major[cluster_anno$fine_annotation])
obj$GMF_major_celltype_11_proB_split <- factor(
  cell_meta$GMF_major_celltype_11_proB_split,
  levels = major_11_levels
)
obj$fine_annotation_res0p3 <- factor(
  cell_meta$fine_annotation,
  levels = cluster_anno$fine_annotation[order(as.integer(cluster_anno$raw_cluster))]
)

cell_meta_out <- file.path(curated_dir, "GSE169396_res0p3_curated_cell_metadata_for_GMF_proB_split_11_major.csv")
cluster_anno_out <- file.path(curated_dir, "GSE169396_res0p3_curated_cluster_annotation_for_GMF_proB_split_11_major.csv")
cell_counts_out <- file.path(curated_dir, "GSE169396_res0p3_proB_split_11_major_cell_counts.csv")
write.csv(cell_meta, cell_meta_out, row.names = FALSE)
write.csv(cluster_anno, cluster_anno_out, row.names = FALSE)
write.csv(
  data.frame(
    GMF_major_celltype_11_proB_split = names(table(obj$GMF_major_celltype_11_proB_split)),
    n_cells = as.integer(table(obj$GMF_major_celltype_11_proB_split)),
    stringsAsFactors = FALSE
  ),
  cell_counts_out,
  row.names = FALSE
)

marker_panel <- list(
  "Granulocytes / neutrophils" = c("LTF", "MMP8", "LCN2", "CEACAM8"),
  "T cells / cytotoxic lymphocytes" = c("CD3D", "TRAC", "IL7R", "NKG7"),
  "Monocyte/macrophage/DC-like myeloid" = c("LYZ", "CD1C", "C1QA", "C1QC"),
  "B-lineage cells" = c("RAG1", "RAG2", "TCL1A", "CD79A"),
  "Pro-B / early B-lineage cells" = c("DNTT", "VPREB1", "IGLL1", "PTCRA"),
  "Plasma cells" = c("MZB1", "XBP1", "JCHAIN", "TNFRSF17"),
  "Erythroid cells" = c("HBA1", "HBB", "ALAS2", "AHSP"),
  "Osteoblastic lineage cells" = c("COL1A1", "SPARC", "IBSP", "BGLAP"),
  "Endothelial cells" = c("PECAM1", "VWF", "CLDN5", "EMCN"),
  "Mast cells" = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  "pDCs" = c("LILRA4", "CLEC4C", "IL3RA", "TCF4")
)

marker_table <- do.call(
  rbind,
  lapply(names(marker_panel), function(cell_type) {
    data.frame(
      GMF_major_celltype_11_proB_split = cell_type,
      selected_rank = seq_along(marker_panel[[cell_type]]),
      marker = marker_panel[[cell_type]],
      stringsAsFactors = FALSE
    )
  })
)
missing_markers <- setdiff(marker_table$marker, rownames(obj))
if (length(missing_markers) > 0) {
  stop("Selected markers absent from Seurat object: ", paste(missing_markers, collapse = ", "))
}

norm_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
raw_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
plot_features <- unique(marker_table$marker)

metric_rows <- list()
for (cell_type in major_11_levels) {
  cells <- colnames(obj)[obj$GMF_major_celltype_11_proB_split == cell_type]
  if (length(cells) == 0) stop("No cells found for 11-major cell type: ", cell_type)
  for (gene in plot_features) {
    expr <- as.numeric(norm_mat[gene, cells, drop = TRUE])
    raw <- as.numeric(raw_mat[gene, cells, drop = TRUE])
    metric_rows[[length(metric_rows) + 1]] <- data.frame(
      cell_type = cell_type,
      gene_symbol = gene,
      n_cells = length(cells),
      pct_expressing = mean(raw > 0) * 100,
      avg_log_expression = mean(expr),
      stringsAsFactors = FALSE
    )
  }
}

dotplot_values <- do.call(rbind, metric_rows)
dotplot_values$avg_expression_scaled_0_1 <- ave(
  dotplot_values$avg_log_expression,
  dotplot_values$gene_symbol,
  FUN = function(x) {
    if (max(x, na.rm = TRUE) > min(x, na.rm = TRUE)) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    } else {
      rep(0, length(x))
    }
  }
)
dotplot_values$cell_type <- factor(dotplot_values$cell_type, levels = rev(major_11_levels))
dotplot_values$cell_type_display <- factor(
  unname(celltype_short_labels[as.character(dotplot_values$cell_type)]),
  levels = rev(short_11_levels)
)
dotplot_values$gene_symbol <- factor(dotplot_values$gene_symbol, levels = plot_features)

selected_marker_out <- file.path(curated_dir, "GSE169396_res0p3_proB_split_11_major_selected_marker_table.csv")
dotplot_values_out <- file.path(curated_dir, "GSE169396_res0p3_proB_split_11_major_marker_dotplot_values_set1_style.csv")
write.csv(marker_table, selected_marker_out, row.names = FALSE)
write.csv(dotplot_values, dotplot_values_out, row.names = FALSE)

umap_png <- file.path(fig_dir, "GSE169396_res0p3_curated_UMAP_GMF_major_celltype_proB_split_11.png")
umap_pdf <- file.path(fig_dir, "GSE169396_res0p3_curated_UMAP_GMF_major_celltype_proB_split_11.pdf")
obj$GMF_major_celltype_plot <- factor(
  unname(celltype_short_labels[as.character(obj$GMF_major_celltype_11_proB_split)]),
  levels = short_11_levels
)
p_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "GMF_major_celltype_plot",
  label = FALSE,
  raster = FALSE,
  cols = celltype_colors
) +
  ggtitle("GSE169396 scRNA cell type overview") +
  labs(color = "Cell type", x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    axis.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.35, "cm")
  )
png(umap_png, width = 3000, height = 2100, res = 180)
print(p_umap)
dev.off()
ggsave(umap_pdf, p_umap, width = 3000 / 180, height = 2100 / 180, bg = "white")

dotplot_png <- file.path(fig_dir, "GSE169396_res0p3_proB_split_11_major_marker_dotplot_set1_style.png")
dotplot_pdf <- file.path(fig_dir, "GSE169396_res0p3_proB_split_11_major_marker_dotplot_set1_style.pdf")
p_set1 <- ggplot(dotplot_values, aes(x = gene_symbol, y = cell_type_display)) +
  geom_point(
    aes(size = pct_expressing, fill = avg_expression_scaled_0_1),
    shape = 21,
    color = "grey70",
    stroke = 0.25,
    alpha = 0.98
  ) +
  scale_size(
    name = "Fraction of cells\nin group (%)",
    range = c(0.15, 3.1),
    breaks = c(20, 40, 60, 80),
    limits = c(0, 100)
  ) +
  scale_fill_gradientn(
    name = "Mean expression\nin group",
    colors = c("#FFFFFF", "#FDE0DD", "#F76868", "#CB181D"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1.0)
  ) +
  labs(title = "Expression of marker genes in cell types", x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      order = 1,
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(fill = "white", color = "grey45", alpha = 1)
    ),
    fill = guide_colorbar(
      order = 2,
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth = unit(0.9, "cm"),
      barheight = unit(0.18, "cm")
    )
  ) +
  theme_classic(base_size = 6) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8, margin = margin(b = 1)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.2, color = "black"),
    axis.text.y = element_text(size = 5.6, color = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.title = element_text(size = 5.6, hjust = 0.5),
    legend.text = element_text(size = 5.4),
    legend.spacing.x = unit(0.12, "cm"),
    legend.key.height = unit(0.12, "cm"),
    legend.key.width = unit(0.18, "cm"),
    legend.box.margin = margin(0, 0, -3, 0),
    plot.margin = margin(1, 3, 2, 3)
  )
ggsave(dotplot_pdf, p_set1, width = 9.8, height = 3.35, bg = "white")
ggsave(dotplot_png, p_set1, width = 9.8, height = 3.35, dpi = 700, bg = "white")

seurat_dotplot_png <- file.path(fig_dir, "GSE169396_res0p3_proB_split_11_major_marker_dotplot_seurat.png")
png(seurat_dotplot_png, width = 4400, height = 2100, res = 180)
print(
  DotPlot(obj, features = plot_features, group.by = "GMF_major_celltype_11_proB_split") +
    RotatedAxis() +
    ggtitle("GSE169396 selected markers by 11 GMF major cell types") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 9)
    )
)
dev.off()

cat("Wrote 11-major cell metadata: ", cell_meta_out, "\n", sep = "")
cat("Wrote 11-major cluster annotation: ", cluster_anno_out, "\n", sep = "")
cat("Wrote 11-major cell counts: ", cell_counts_out, "\n", sep = "")
cat("Wrote selected marker table: ", selected_marker_out, "\n", sep = "")
cat("Wrote dotplot values: ", dotplot_values_out, "\n", sep = "")
cat("Wrote UMAP PNG: ", umap_png, "\n", sep = "")
cat("Wrote UMAP PDF: ", umap_pdf, "\n", sep = "")
cat("Wrote Set1-style dotplot PNG: ", dotplot_png, "\n", sep = "")
cat("Wrote Set1-style dotplot PDF: ", dotplot_pdf, "\n", sep = "")
cat("Wrote Seurat dotplot: ", seurat_dotplot_png, "\n", sep = "")
cat("11-major cell counts:\n")
print(table(obj$GMF_major_celltype_11_proB_split))
