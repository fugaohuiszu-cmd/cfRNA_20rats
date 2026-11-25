#!/usr/bin/env Rscript
# 用途：绘制细胞通讯结果。
# 输入：已保存的CellChat对象和结果表；输出：通讯网络、通路和配体受体图。

rm(list = ls())

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

work_dir <- "/home/fugaohui/rat-siglecell-cfRNA/CellFreeGMF_work/20sample_OP_vs_N_pvalue_GMF"
cellchat_dir <- file.path(work_dir, "output", "cellchat")
table_dir <- file.path(cellchat_dir, "tables")
object_dir <- file.path(cellchat_dir, "objects")
figure_dir <- file.path(cellchat_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

object_paths <- c(
  N = file.path(object_dir, "cellchat_N.rds"),
  OP = file.path(object_dir, "cellchat_OP.rds"),
  merged = file.path(object_dir, "cellchat_merged.rds")
)
missing_objects <- object_paths[!file.exists(object_paths)]
if (length(missing_objects) > 0) {
  stop("Missing CellChat objects: ", paste(missing_objects, collapse = "; "))
}

cellchat_N <- readRDS(object_paths[["N"]])
cellchat_OP <- readRDS(object_paths[["OP"]])
cellchat_merged <- readRDS(object_paths[["merged"]])

cell_colors <- c(
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

cell_levels <- levels(cellchat_N@idents)
if (!identical(cell_levels, levels(cellchat_OP@idents))) {
  stop("N and OP CellChat objects have different cell-type levels")
}
if (!all(cell_levels %in% names(cell_colors))) {
  stop("Missing colors for: ", paste(setdiff(cell_levels, names(cell_colors)), collapse = ", "))
}
color_use <- cell_colors[cell_levels]

save_circle <- function(cellchat, group_name, group_size, edge_max) {
  plot_one <- function() {
    netVisual_circle(
      cellchat@net$weight,
      vertex.weight = group_size,
      vertex.weight.max = max(c(table(cellchat_N@idents), table(cellchat_OP@idents))),
      vertex.size.max = 18,
      edge.weight.max = edge_max,
      weight.scale = TRUE,
      label.edge = FALSE,
      color.use = color_use,
      vertex.label.cex = 0.72,
      title.name = NULL
    )
    title(main = paste("Interaction for", group_name, "samples"), line = 0.5, cex.main = 1.25, font.main = 2)
  }

  png(file.path(figure_dir, paste0(group_name, "_network.png")), width = 2400, height = 2100, res = 260)
  plot_one()
  dev.off()
  pdf(file.path(figure_dir, paste0(group_name, "_network.pdf")), width = 9.0, height = 7.8)
  plot_one()
  dev.off()
}

edge_max <- max(cellchat_N@net$weight, cellchat_OP@net$weight, na.rm = TRUE)
save_circle(cellchat_N, "N", as.numeric(table(cellchat_N@idents)), edge_max)
save_circle(cellchat_OP, "OP", as.numeric(table(cellchat_OP@idents)), edge_max)

group_values <- read_csv(file.path(table_dir, "group_interaction_values.csv"), show_col_types = FALSE) |>
  mutate(group = factor(group, levels = c("N", "OP")))

p_strength <- ggplot(group_values, aes(group, relative_weight, fill = group)) +
  geom_col(width = 0.54, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", relative_weight)), vjust = -0.45, size = 4.0, fontface = "bold") +
  scale_fill_manual(values = c(N = "#22B3B8", OP = "#EE6A65")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Group interaction comparison", x = NULL, y = "Relative interaction strength") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )
ggsave(file.path(figure_dir, "interaction_strength.png"), p_strength, width = 4.0, height = 4.5, dpi = 500, bg = "white")
ggsave(file.path(figure_dir, "interaction_strength.pdf"), p_strength, width = 4.0, height = 4.5, bg = "white")

png(file.path(figure_dir, "differential_network.png"), width = 2400, height = 2200, res = 260)
netVisual_diffInteraction(
  cellchat_merged,
  comparison = c(1, 2),
  weight.scale = TRUE,
  measure = "weight",
  color.use = color_use,
  title.name = "Differential interaction strength"
)
dev.off()
pdf(file.path(figure_dir, "differential_network.pdf"), width = 9.0, height = 8.2)
netVisual_diffInteraction(
  cellchat_merged,
  comparison = c(1, 2),
  weight.scale = TRUE,
  measure = "weight",
  color.use = color_use,
  title.name = "Differential interaction strength"
)
dev.off()

png(file.path(figure_dir, "differential_heatmap.png"), width = 2400, height = 2050, res = 260)
print(netVisual_heatmap(
  cellchat_merged,
  comparison = c(1, 2),
  measure = "weight",
  color.use = color_use,
  title.name = "Differential interaction strength"
))
dev.off()
pdf(file.path(figure_dir, "differential_heatmap.pdf"), width = 9.0, height = 7.7)
print(netVisual_heatmap(
  cellchat_merged,
  comparison = c(1, 2),
  measure = "weight",
  color.use = color_use,
  title.name = "Differential interaction strength"
))
dev.off()

source_cell <- "Osteogenic cells"
if (!(source_cell %in% levels(cellchat_merged@idents$joint))) {
  stop("Osteogenic cells are absent from the merged CellChat object")
}
target_use <- seq_along(levels(cellchat_merged@idents$joint))

p_increased <- netVisual_bubble(
  cellchat_merged,
  sources.use = source_cell,
  targets.use = target_use,
  comparison = c(1, 2),
  max.dataset = 2,
  title.name = "Increased signaling in OP: Osteogenic cells",
  angle.x = 45,
  remove.isolate = TRUE
)
p_decreased <- netVisual_bubble(
  cellchat_merged,
  sources.use = source_cell,
  targets.use = target_use,
  comparison = c(1, 2),
  max.dataset = 1,
  title.name = "Decreased signaling in OP: Osteogenic cells",
  angle.x = 45,
  remove.isolate = TRUE
)

p_increased <- p_increased +
  theme(plot.margin = margin(t = 10, r = 25, b = 100, l = 100))
p_decreased <- p_decreased +
  theme(plot.margin = margin(t = 10, r = 25, b = 100, l = 100))

ggsave(file.path(figure_dir, "osteogenic_LR_increased.png"), p_increased, width = 12.5, height = 8.5, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "osteogenic_LR_increased.pdf"), p_increased, width = 12.5, height = 8.5, bg = "white")
ggsave(file.path(figure_dir, "osteogenic_LR_decreased.png"), p_decreased, width = 12.5, height = 8.5, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "osteogenic_LR_decreased.pdf"), p_decreased, width = 12.5, height = 8.5, bg = "white")

message("CellChat figures completed: ", figure_dir)
