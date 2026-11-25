# 用途：绘制各类cfRNA差异数量
# 输入：DE_summary.csv
# 输出：差异数量表及图

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
result_dir <- file.path(project_dir, "03_results", "02_DE")
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")
summary_file <- file.path(table_dir, "DE_summary.csv")

if (!file.exists(summary_file)) {
  stop("Missing DE summary table: ", summary_file)
}

summary_df <- read.csv(summary_file, check.names = FALSE)
required_columns <- c("RNA_type", "upregulated", "downregulated")
missing_columns <- setdiff(required_columns, colnames(summary_df))
if (length(missing_columns) > 0L) {
  stop("DE summary table is missing columns: ", paste(missing_columns, collapse = ", "))
}

rna_levels <- c("mRNA", "lncRNA", "miRNA", "piRNA", "snRNA", "snoRNA")
plot_source <- summary_df %>%
  filter(RNA_type %in% rna_levels) %>%
  select(RNA_type, upregulated, downregulated)

if (!setequal(plot_source$RNA_type, rna_levels) || anyDuplicated(plot_source$RNA_type)) {
  stop("DE summary table does not contain one row for each of the six RNA types.")
}

plot_df <- bind_rows(
  plot_source %>%
    transmute(RNA_type, direction = "Upregulated in OP", count = upregulated),
  plot_source %>%
    transmute(RNA_type, direction = "Downregulated in OP", count = downregulated)
) %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Downregulated in OP", "Upregulated in OP")
    ),
    RNA_type = factor(RNA_type, levels = rev(rna_levels)),
    signed_count = ifelse(direction == "Downregulated in OP", -count, count)
  )

if (anyNA(plot_df$count) || any(plot_df$count < 0)) {
  stop("Differential cfRNA counts must be non-negative and complete.")
}

write.csv(
  plot_df %>%
    mutate(RNA_type = as.character(RNA_type)) %>%
    select(RNA_type, direction, count),
  file.path(table_dir, "DE_counts_by_RNA_type.csv"),
  row.names = FALSE
)

max_count <- max(plot_df$count)
label_offset <- max(8, ceiling(max_count * 0.025))
x_limit <- ceiling((max_count + label_offset * 4) / 50) * 50

p <- ggplot(plot_df, aes(x = signed_count, y = RNA_type, fill = direction)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.55) +
  geom_text(
    aes(
      x = signed_count + ifelse(signed_count < 0, -label_offset, label_offset),
      label = count,
      hjust = ifelse(signed_count < 0, 1, 0)
    ),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Downregulated in OP" = "#46559F",
      "Upregulated in OP" = "#870F45"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = scales::breaks_pretty(n = 6),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Differential cfRNAs by RNA type",
    subtitle = "Adjusted P value < 0.05 and |log2FC| > 0.5",
    x = "Number of differential cfRNAs",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10.5),
    legend.position = "top",
    legend.text = element_text(size = 11),
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave(
  file.path(figure_dir, "DE_counts_by_RNA_type.pdf"),
  p, width = 8.2, height = 5.4, device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "DE_counts_by_RNA_type.png"),
  p, width = 8.2, height = 5.4, dpi = 300, bg = "white"
)

cat("Differential cfRNA count plot completed.\n")
print(plot_df %>% select(RNA_type, direction, count) %>% arrange(RNA_type, direction))
