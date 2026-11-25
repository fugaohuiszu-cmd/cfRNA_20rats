suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_dir <- "."
result_dir <- file.path(project_dir, "03_results", "04_source_tracing", "female_reference")
figure_dir <- file.path(result_dir, "figures")
reference <- read.delim(file.path(result_dir, "female_reference_FPKM.tsv"),
                        row.names = 1, check.names = FALSE)
assignment <- read.csv(file.path(result_dir, "tissue_assignment.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(reference) == 44, !anyNA(reference),
          !anyDuplicated(assignment$Geneid),
          all(assignment$Geneid %in% rownames(reference)))

tissues <- c("BoneMarrow", "Cartilage",
             setdiff(colnames(reference), c("BoneMarrow", "Cartilage")))
groups <- c("Bone marrow", "Cartilage", "Other tissues")
group_colors <- c("Bone marrow" = "#D73027", "Cartilage" = "#FC8D59",
                  "Other tissues" = "#D9D9D9")
color_scale <- colorRamp2(c(-2.5, 0, 2.5), c("#3F7FBF", "#FFFFFF", "#B9413E"))
plot_names <- c(Upregulated = "Figure_S4_upregulated_heatmap",
                Downregulated = "Figure_S5_downregulated_heatmap")
expected_counts <- c(Upregulated = 554L, Downregulated = 1227L)
column_group <- factor(c("Bone marrow", "Cartilage", rep("Other tissues", 42)),
                       levels = groups)
dir.create(figure_dir, showWarnings = FALSE)

for (direction in names(plot_names)) {
  selected <- assignment[assignment[["Expression change"]] == direction, ]
  stopifnot(nrow(selected) == expected_counts[[direction]])
  expression <- as.matrix(reference[selected$Geneid, tissues])
  log_expression <- log2(expression + 1)
  stopifnot(all(apply(log_expression, 1, sd) > 0))
  row_z <- t(scale(t(log_expression)))
  row_group <- ifelse(selected[["Highest-expression tissue"]] == "BoneMarrow", "Bone marrow",
                      ifelse(selected[["Highest-expression tissue"]] == "Cartilage", "Cartilage", "Other tissues"))
  row_group <- factor(row_group, levels = groups)
  row_order <- order(row_group, selected[["Highest-expression tissue"]], selected$pvalue)
  row_z <- row_z[row_order, , drop = FALSE]
  row_group <- row_group[row_order]

  top_annotation <- HeatmapAnnotation(
    Tissue = column_group, col = list(Tissue = group_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Tissue", at = groups)
  )
  left_annotation <- rowAnnotation(
    `Highest-expression tissue` = row_group,
    col = list(`Highest-expression tissue` = group_colors),
    show_annotation_name = FALSE, show_legend = FALSE, width = unit(3, "mm")
  )
  heatmap <- Heatmap(
    row_z, name = "Row Z-score", col = color_scale,
    top_annotation = top_annotation, left_annotation = left_annotation,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = row_group, row_title = NULL, row_gap = unit(1.5, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 60, column_names_gp = gpar(fontsize = 7),
    column_title = paste(direction, "plasma cf-mRNAs across 44 female rat tissues"),
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    heatmap_legend_param = list(title = "Row Z-score", at = c(-2.5, 0, 2.5)),
    use_raster = FALSE
  )

  pdf(file.path(figure_dir, paste0(plot_names[[direction]], ".pdf")),
      width = 15, height = 13, useDingbats = FALSE)
  draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(dev.off())
  png(file.path(figure_dir, paste0(plot_names[[direction]], ".png")),
      width = 4500, height = 3900, res = 300)
  draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(dev.off())
  cat(direction, ":", nrow(row_z), "genes x", ncol(row_z), "tissues\n")
}
cat("Figures:", figure_dir, "\n")
