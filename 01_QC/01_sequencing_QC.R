# 用途：测序质控与六类cfRNA组成
# 输入：测序读数统计和原始计数矩阵
# 输出：质控、RNA组成表及图

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
input_dir <- file.path(project_dir, "01_input")
output_dir <- file.path(project_dir, "03_results", "01_QC")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

total_file <- file.path(input_dir, "Total_reads.txt")
clean_file <- file.path(input_dir, "Trimmed_reads.txt")
unmapped_file <- file.path(input_dir, "Unmapped_reads.txt")
matrix_file <- file.path(input_dir, "all_RNA_final_merged.csv")

required_files <- c(total_file, clean_file, unmapped_file, matrix_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- c(paste0("N", 1:10), paste0("OP", 1:10))
sample_group <- data.frame(
  sample_id = sample_ids,
  group = factor(c(rep("N", 10), rep("OP", 10)), levels = c("N", "OP"))
)

read_count_file <- function(path, value_name) {
  x <- read.table(
    path,
    header = FALSE,
    sep = "",
    quote = "",
    comment.char = "",
    colClasses = c("numeric", "character")
  )
  if (ncol(x) != 2L) {
    stop("Expected exactly two columns in: ", path)
  }

  ids <- sub("_.*$", "", basename(x[[2]]))
  if (anyDuplicated(ids)) {
    stop("Duplicated sample IDs in: ", path)
  }
  if (!setequal(ids, sample_ids)) {
    stop(
      "Sample IDs do not match N1-N10 and OP1-OP10 in: ", path,
      "\nObserved: ", paste(ids, collapse = ", ")
    )
  }

  out <- data.frame(sample_id = ids, value = as.numeric(x[[1]]))
  names(out)[2] <- value_name
  out[match(sample_ids, out$sample_id), , drop = FALSE]
}

total_reads <- read_count_file(total_file, "total_reads")
clean_reads <- read_count_file(clean_file, "clean_reads")
unmapped_reads <- read_count_file(unmapped_file, "genome_unmapped_reads")

qc <- sample_group %>%
  left_join(total_reads, by = "sample_id") %>%
  left_join(clean_reads, by = "sample_id") %>%
  left_join(unmapped_reads, by = "sample_id") %>%
  mutate(
    low_quality_reads = total_reads - clean_reads,
    genome_mapped_reads = clean_reads - genome_unmapped_reads,
    clean_ratio = clean_reads / total_reads,
    mapping_rate = genome_mapped_reads / clean_reads,
    genome_mapped_proportion = genome_mapped_reads / total_reads,
    genome_unmapped_proportion = genome_unmapped_reads / total_reads,
    low_quality_proportion = low_quality_reads / total_reads
  )

if (anyNA(qc)) {
  stop("Missing values were produced while joining QC inputs.")
}
if (any(qc$total_reads <= 0) || any(qc$clean_reads <= 0)) {
  stop("Total and clean read counts must be positive for every sample.")
}
if (any(qc$low_quality_reads < 0) || any(qc$genome_mapped_reads < 0)) {
  stop("Negative read categories detected; QC inputs are inconsistent.")
}
if (any(abs(
  qc$genome_mapped_proportion + qc$genome_unmapped_proportion +
    qc$low_quality_proportion - 1
) > 1e-10)) {
  stop("Read composition does not sum to 100% for every sample.")
}

rna <- read.csv(matrix_file, check.names = FALSE)
required_matrix_columns <- c("RNA_type", "geneid", sample_ids)
missing_matrix_columns <- setdiff(required_matrix_columns, colnames(rna))
if (length(missing_matrix_columns) > 0L) {
  stop(
    "Missing required columns in all_RNA_final_merged.csv: ",
    paste(missing_matrix_columns, collapse = ", ")
  )
}

six_types <- c("mRNA", "lncRNA", "snoRNA", "snRNA", "miRNA", "piRNA")
rna6 <- rna[rna$RNA_type %in% six_types, c("RNA_type", sample_ids), drop = FALSE]
if (!setequal(unique(rna6$RNA_type), six_types)) {
  stop("Not all six expected RNA types are present in the expression matrix.")
}

count_matrix <- as.matrix(rna6[, sample_ids, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix) || any(count_matrix < 0)) {
  stop("The six-type count matrix contains missing or negative values.")
}

type_counts <- rowsum(count_matrix, group = rna6$RNA_type, reorder = FALSE)
type_counts <- type_counts[six_types, sample_ids, drop = FALSE]

composition <- as.data.frame(type_counts) %>%
  mutate(RNA_type = rownames(type_counts), .before = 1) %>%
  pivot_longer(
    cols = all_of(sample_ids),
    names_to = "sample_id",
    values_to = "RNA_type_count"
  ) %>%
  left_join(sample_group, by = "sample_id") %>%
  group_by(sample_id, group) %>%
  mutate(
    six_type_total_count = sum(RNA_type_count),
    relative_proportion = RNA_type_count / six_type_total_count,
    relative_percent = relative_proportion * 100
  ) %>%
  ungroup() %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_ids),
    RNA_type = factor(RNA_type, levels = six_types)
  )

if (any(!is.finite(composition$relative_proportion))) {
  stop("Non-finite six-type relative proportions were produced.")
}

qc_group_summary <- qc %>%
  group_by(group) %>%
  summarise(
    n = n(),
    total_reads_mean = mean(total_reads),
    total_reads_sd = sd(total_reads),
    total_reads_median = median(total_reads),
    clean_reads_mean = mean(clean_reads),
    clean_reads_sd = sd(clean_reads),
    clean_reads_median = median(clean_reads),
    clean_ratio_mean_percent = mean(clean_ratio) * 100,
    clean_ratio_sd_percent = sd(clean_ratio) * 100,
    clean_ratio_median_percent = median(clean_ratio) * 100,
    mapping_rate_mean_percent = mean(mapping_rate) * 100,
    mapping_rate_sd_percent = sd(mapping_rate) * 100,
    mapping_rate_median_percent = median(mapping_rate) * 100,
    .groups = "drop"
  )

composition_group_summary <- composition %>%
  group_by(group, RNA_type) %>%
  summarise(
    n = n(),
    relative_percent_mean = mean(relative_percent),
    relative_percent_sd = sd(relative_percent),
    relative_percent_median = median(relative_percent),
    relative_percent_min = min(relative_percent),
    relative_percent_max = max(relative_percent),
    .groups = "drop"
  )

qc_export <- qc %>%
  mutate(
    group = as.character(group),
    clean_ratio_percent = clean_ratio * 100,
    mapping_rate_percent = mapping_rate * 100,
    genome_mapped_percent_of_total = genome_mapped_proportion * 100,
    genome_unmapped_percent_of_total = genome_unmapped_proportion * 100,
    low_quality_percent_of_total = low_quality_proportion * 100
  ) %>%
  select(
    sample_id, group, total_reads, clean_reads, low_quality_reads,
    genome_mapped_reads, genome_unmapped_reads,
    clean_ratio_percent, mapping_rate_percent,
    genome_mapped_percent_of_total, genome_unmapped_percent_of_total,
    low_quality_percent_of_total
  )

composition_export <- composition %>%
  transmute(
    sample_id = as.character(sample_id),
    group = as.character(group),
    RNA_type = as.character(RNA_type),
    RNA_type_count,
    six_type_total_count,
    relative_percent
  )

write.csv(
  qc_export,
  file.path(table_dir, "sequencing_QC_per_sample.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  composition_export,
  file.path(table_dir, "six_cfRNA_composition_per_sample.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  qc_group_summary,
  file.path(table_dir, "sequencing_QC_group_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  composition_group_summary,
  file.path(table_dir, "six_cfRNA_composition_group_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

group_colors <- c(N = "#D73027", OP = "#4575B4")
group_shapes <- c(N = 16, OP = 17)
read_colors <- c(
  "Genome-mapped" = "#A11D1F",
  "Genome-unmapped" = "#F4B7A0",
  "Low-quality" = "#85B8CC"
)
rna_colors <- c(
  mRNA = "#8EC9EB",
  lncRNA = "#3B8CCB",
  snoRNA = "#9A72C7",
  snRNA = "#D0BBDB",
  miRNA = "#FCA3A3",
  piRNA = "#ED5F5F"
)

box_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(size = 18, face = "plain", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, colour = "#4D4D4D"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

composition_theme <- theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, colour = "#333333"),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#BDBDBD"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 12, 8, 10)
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

read_composition <- qc %>%
  select(
    sample_id, group,
    `Genome-mapped` = genome_mapped_proportion,
    `Genome-unmapped` = genome_unmapped_proportion,
    `Low-quality` = low_quality_proportion
  ) %>%
  pivot_longer(
    cols = all_of(names(read_colors)),
    names_to = "read_category",
    values_to = "proportion"
  ) %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_ids),
    read_category = factor(read_category, levels = names(read_colors))
  )

p_read_composition <- ggplot(
  read_composition,
  aes(x = sample_id, y = proportion * 100, fill = read_category)
) +
  geom_col(width = 0.86, colour = "white", linewidth = 0.25) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = read_colors, name = "Read category") +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = label_number(suffix = "%"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Read quality and genome-alignment composition",
    x = NULL,
    y = "Proportion of total reads"
  ) +
  composition_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

save_plot(p_read_composition, "Fig1_read_quality_composition", 10.5, 5.6)

read_count_long <- qc %>%
  select(sample_id, group, `Total Reads` = total_reads, `Clean Reads` = clean_reads) %>%
  pivot_longer(
    cols = c(`Total Reads`, `Clean Reads`),
    names_to = "metric",
    values_to = "reads"
  ) %>%
  mutate(metric = factor(metric, levels = c("Total Reads", "Clean Reads")))

p_read_counts <- ggplot(read_count_long, aes(x = group, y = reads, colour = group)) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(shape = group), size = 2.5, width = 0.2, stroke = 0.5) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title = "N vs OP",
    x = NULL,
    y = "Read Count"
  ) +
  box_theme

save_plot(p_read_counts, "Fig2_raw_and_clean_reads", 8, 4)

rate_long <- qc %>%
  select(sample_id, group, `Clean ratio` = clean_ratio, `Mapping rate` = mapping_rate) %>%
  pivot_longer(
    cols = c(`Clean ratio`, `Mapping rate`),
    names_to = "metric",
    values_to = "rate"
  ) %>%
  mutate(metric = factor(metric, levels = c("Clean ratio", "Mapping rate")))

p_rates <- ggplot(rate_long, aes(x = group, y = rate * 100, colour = group)) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(shape = group), size = 2.5, width = 0.2, stroke = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("N", "OP")),
    size = 3
  ) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.22))
  ) +
  labs(
    title = "N vs OP",
    x = NULL,
    y = "Rate (%)"
  ) +
  box_theme

save_plot(p_rates, "Fig3_clean_ratio_and_mapping_rate", 8, 4)

p_rna_composition <- ggplot(
  composition,
  aes(x = sample_id, y = relative_percent, fill = RNA_type)
) +
  geom_col(width = 0.86, colour = "white", linewidth = 0.25) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = rna_colors, name = "RNA type", drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = label_number(suffix = "%"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Relative composition of six cfRNA types",
    x = NULL,
    y = "Relative abundance among six types"
  ) +
  composition_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

save_plot(p_rna_composition, "Fig4_six_cfRNA_relative_composition", 10.5, 5.6)

rna_ratio_order <- c("miRNA", "piRNA", "mRNA", "lncRNA", "snRNA", "snoRNA")
composition_box <- composition %>%
  mutate(RNA_type = factor(as.character(RNA_type), levels = rna_ratio_order))

p_rna_ratio <- ggplot(
  composition_box,
  aes(x = group, y = relative_percent, colour = group)
) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(shape = group), size = 2.5, width = 0.2, stroke = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("N", "OP")),
    size = 3
  ) +
  facet_wrap(~RNA_type, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "Relative abundance (%)") +
  box_theme +
  theme(plot.title = element_blank())

save_plot(p_rna_ratio, "Fig5_six_cfRNA_relative_abundance_by_group", 12, 6)

cat("QC analysis completed successfully.\n")
cat("Samples: ", nrow(qc), " (N=10, OP=10)\n", sep = "")
cat("Six-type features: ", nrow(rna6), "\n", sep = "")
cat(
  "Clean ratio range: ",
  sprintf("%.2f%%-%.2f%%", min(qc$clean_ratio) * 100, max(qc$clean_ratio) * 100),
  "\n",
  sep = ""
)
cat(
  "Mapping rate range: ",
  sprintf("%.2f%%-%.2f%%", min(qc$mapping_rate) * 100, max(qc$mapping_rate) * 100),
  "\n",
  sep = ""
)
cat("Figures: ", normalizePath(figure_dir, winslash = "/"), "\n", sep = "")
cat("Tables: ", normalizePath(table_dir, winslash = "/"), "\n", sep = "")
