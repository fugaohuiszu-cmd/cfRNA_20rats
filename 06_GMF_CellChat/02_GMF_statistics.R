# 用途：比较两组细胞来源评分。
# 输入：GMF_scores_all.csv和细胞类型配色。
# 输出：组间统计表和细胞来源图。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

set.seed(2025)

project_dir <- "F:/DE-analysis/final/cfRNA_project"
score_file <- file.path(
  project_dir,
  "03_results",
  "05_cell_source",
  "matrices",
  "GMF_scores_all.csv"
)
color_file <- file.path(
  project_dir,
  "01_input",
  "reference",
  "GSE169396_celltype_colors.csv"
)
figure_dir <- file.path(project_dir, "03_results", "05_cell_source", "figures")
table_dir <- file.path(project_dir, "03_results", "05_cell_source", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(score_file)) stop("Missing GMF score file: ", score_file)
if (!file.exists(color_file)) stop("Missing cell-type color file: ", color_file)

scores <- read_csv(score_file, show_col_types = FALSE)
colors <- read_csv(color_file, show_col_types = FALSE) |>
  arrange(order)

required_score_columns <- c(
  "sample_id", "group", "cell_type", "n_reference_cells", "GMF_score"
)
missing_score_columns <- setdiff(required_score_columns, names(scores))
if (length(missing_score_columns) > 0) {
  stop("GMF score file missing columns: ", paste(missing_score_columns, collapse = ", "))
}

cell_order <- colors$cell_type_short
if (!setequal(unique(scores$cell_type), cell_order)) {
  stop("GMF score cell types do not match the 11-cell-type color table")
}
if (!identical(sort(unique(scores$group)), c("N", "OP"))) {
  stop("GMF score groups must be N and OP")
}
if (n_distinct(scores$sample_id) != 20) {
  stop("GMF score table must contain 20 samples")
}

cliff_delta <- function(op, n) {
  comparison <- outer(op, n, "-")
  (sum(comparison > 0) - sum(comparison < 0)) / (length(op) * length(n))
}

statistics <- bind_rows(lapply(cell_order, function(cell_type_name) {
  current <- filter(scores, cell_type == cell_type_name)
  n_values <- current$GMF_score[current$group == "N"]
  op_values <- current$GMF_score[current$group == "OP"]
  if (length(n_values) != 10 || length(op_values) != 10) {
    stop("Each cell type must contain 10 N and 10 OP scores: ", cell_type_name)
  }
  data.frame(
    cell_type = cell_type_name,
    n_reference_cells = unique(current$n_reference_cells),
    N_n = length(n_values),
    OP_n = length(op_values),
    N_mean = mean(n_values),
    OP_mean = mean(op_values),
    OP_minus_N = mean(op_values) - mean(n_values),
    N_median = median(n_values),
    OP_median = median(op_values),
    wilcoxon_p = wilcox.test(op_values, n_values, exact = FALSE)$p.value,
    cliff_delta = cliff_delta(op_values, n_values),
    stringsAsFactors = FALSE
  )
})) |>
  mutate(BH_FDR = p.adjust(wilcoxon_p, method = "BH")) |>
  arrange(match(cell_type, cell_order))

format_p <- function(p) {
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.3f", p)
}

scores <- scores |>
  mutate(
    group = factor(group, levels = c("N", "OP")),
    cell_type = factor(cell_type, levels = cell_order)
  )

plot_labels <- statistics |>
  mutate(
    cell_type = factor(cell_type, levels = cell_order),
    label = paste0("P = ", vapply(wilcoxon_p, format_p, character(1)))
  ) |>
  left_join(
    scores |>
      group_by(cell_type) |>
      summarise(y = max(GMF_score), .groups = "drop"),
    by = "cell_type"
  )

group_colors <- c(N = "#4C78A8", OP = "#E45756")

p <- ggplot(scores, aes(x = group, y = GMF_score, color = group)) +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.45
  ) +
  geom_jitter(width = 0.14, size = 1.5, alpha = 0.82) +
  geom_text(
    data = plot_labels,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    color = "#B2182B",
    fontface = "bold",
    size = 3.1,
    vjust = -0.45
  ) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "CellFreeGMF cell source scores",
    x = "Group",
    y = "GMF score"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "#F0F0F0", color = "#BDBDBD"),
    strip.text = element_text(face = "bold", size = 9.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "GMF_scores_boxplot.pdf"),
  p,
  width = 14,
  height = 9,
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "GMF_scores_boxplot.png"),
  p,
  width = 14,
  height = 9,
  dpi = 300,
  bg = "white"
)

write_csv(scores |> mutate(group = as.character(group), cell_type = as.character(cell_type)),
          file.path(table_dir, "GMF_sample_scores.csv"))
write_csv(statistics, file.path(table_dir, "GMF_group_statistics.csv"))

print(statistics)
