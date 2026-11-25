#!/usr/bin/env Rscript
# 用途：绘制本研究使用的单细胞参考。
# 输入：cell_metadata_11.csv、配色和marker_expression.csv；输出：UMAP、样本细胞组成和标记基因图。

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) > 0) normalizePath(args[[1]]) else normalizePath(".")

metadata_path <- file.path(
  project_dir,
  "03_results",
  "05_cell_source",
  "matrices",
  "cell_metadata_11.csv"
)
color_path <- file.path(
  project_dir,
  "01_input",
  "reference",
  "GSE169396_celltype_colors.csv"
)
output_dir <- file.path(
  project_dir,
  "03_results",
  "05_cell_source",
  "scrna"
)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
marker_path <- file.path(table_dir, "marker_expression.csv")

required_files <- c(metadata_path, color_path, marker_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = "; "))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
colors <- read.csv(color_path, stringsAsFactors = FALSE, check.names = FALSE)
marker_values <- read.csv(marker_path, stringsAsFactors = FALSE, check.names = FALSE)

required_metadata <- c("cell_id", "sample_id", "umap_1", "umap_2", "cell_type")
required_colors <- c("order", "cell_type_short", "umap_color")
required_markers <- c(
  "cell_type_short",
  "gene_symbol",
  "pct_expressing",
  "avg_expression_scaled_0_1"
)
if (length(setdiff(required_metadata, colnames(metadata))) > 0) {
  stop("Cell metadata columns do not match the expected structure")
}
if (length(setdiff(required_colors, colnames(colors))) > 0) {
  stop("Cell-type color columns do not match the expected structure")
}
if (length(setdiff(required_markers, colnames(marker_values))) > 0) {
  stop("Marker-expression columns do not match the expected structure")
}
if (anyDuplicated(metadata$cell_id) > 0) {
  stop("Duplicated cell IDs in cell metadata")
}
if (any(!is.finite(metadata$umap_1)) || any(!is.finite(metadata$umap_2))) {
  stop("Non-finite UMAP coordinates in cell metadata")
}

colors <- colors[order(colors$order), ]
cell_levels <- colors$cell_type_short
cell_colors <- setNames(colors$umap_color, colors$cell_type_short)
if (!setequal(unique(metadata$cell_type), cell_levels)) {
  stop("Cell types in metadata do not match the 11-cell-type color map")
}
if (!all(unique(marker_values$cell_type_short) %in% cell_levels)) {
  stop("Marker table contains cell types outside the 11-cell-type reference")
}

metadata$cell_type <- factor(metadata$cell_type, levels = cell_levels)

p_umap <- ggplot(metadata, aes(umap_1, umap_2, color = cell_type)) +
  geom_point(size = 0.18, alpha = 0.95) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  coord_equal() +
  labs(
    title = "GSE169396 scRNA cell type overview",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.38, "cm"),
    legend.key.width = unit(0.38, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(figure_dir, "UMAP.png"), p_umap, width = 12, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "UMAP.pdf"), p_umap, width = 12, height = 8.5, bg = "white")

sample_levels <- sort(unique(metadata$sample_id))
composition <- as.data.frame(
  table(
    sample_id = factor(metadata$sample_id, levels = sample_levels),
    cell_type = factor(metadata$cell_type, levels = cell_levels)
  ),
  stringsAsFactors = FALSE
)
colnames(composition)[3] <- "n_cells"
sample_totals <- aggregate(n_cells ~ sample_id, composition, sum)
colnames(sample_totals)[2] <- "total_cells"
composition <- left_join(composition, sample_totals, by = "sample_id") |>
  mutate(
    proportion = n_cells / total_cells,
    percent = proportion * 100,
    percent_label = ifelse(percent >= 7, sprintf("%.1f%%", percent), ""),
    sample_id = factor(sample_id, levels = sample_levels),
    cell_type = factor(cell_type, levels = cell_levels)
  )

write.csv(composition, file.path(table_dir, "cell_composition.csv"), row.names = FALSE)

sample_labels <- setNames(
  paste0(sample_totals$sample_id, "\nn=", sample_totals$total_cells),
  sample_totals$sample_id
)

p_composition <- ggplot(composition, aes(sample_id, percent, fill = cell_type)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = percent_label),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 3.4
  ) +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_x_discrete(labels = sample_labels) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "GSE169396 scRNA cell-type composition by sample",
    x = "Sample",
    y = "Cell-type proportion (%)",
    fill = "Cell type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(figure_dir, "sample_composition.png"), p_composition, width = 10, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "sample_composition.pdf"), p_composition, width = 10, height = 6.5, bg = "white")

gene_levels <- unique(marker_values$gene_symbol)
marker_values$cell_type_short <- factor(marker_values$cell_type_short, levels = rev(cell_levels))
marker_values$gene_symbol <- factor(marker_values$gene_symbol, levels = gene_levels)

p_dotplot <- ggplot(marker_values, aes(gene_symbol, cell_type_short)) +
  geom_point(
    aes(size = pct_expressing, fill = avg_expression_scaled_0_1),
    shape = 21,
    color = "grey70",
    stroke = 0.25
  ) +
  scale_size(
    name = "Fraction of cells\nin group (%)",
    range = c(0.15, 3.2),
    limits = c(0, 100),
    breaks = c(20, 40, 60, 80)
  ) +
  scale_fill_gradientn(
    name = "Mean expression\nin group",
    colors = c("#FFFFFF", "#FDE0DD", "#F76868", "#CB181D"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(title = "Expression of marker genes in cell types", x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      order = 1,
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(fill = "white", color = "grey45")
    ),
    fill = guide_colorbar(
      order = 2,
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth = grid::unit(1.0, "cm"),
      barheight = grid::unit(0.2, "cm")
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45),
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left"
  )

ggsave(file.path(figure_dir, "marker_dotplot.png"), p_dotplot, width = 13.5, height = 4.2, dpi = 500, bg = "white")
ggsave(file.path(figure_dir, "marker_dotplot.pdf"), p_dotplot, width = 13.5, height = 4.2, bg = "white")

cat("Cells:", nrow(metadata), "\n")
cat("Samples:", paste(sample_levels, collapse = ", "), "\n")
cat("Cell types:", length(cell_levels), "\n")
cat("Output:", output_dir, "\n")
