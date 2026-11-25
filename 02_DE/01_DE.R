options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- "."
input_file <- file.path(project_dir, "01_input", "all_RNA_final_merged.csv")
result_dir <- file.path(project_dir, "03_results", "02_DE")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
rds_dir <- file.path(result_dir, "rds")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- c(paste0("N", 1:10), paste0("OP", 1:10))
six_types <- c("mRNA", "lncRNA", "snoRNA", "snRNA", "miRNA", "piRNA")
fdr_cutoff <- 0.05
lfc_cutoff <- 0.5
min_count_samples <- 1L

if (!file.exists(input_file)) {
  stop("Missing input matrix: ", input_file)
}

dat <- read.csv(input_file, check.names = FALSE)
required_columns <- c("RNA_type", "geneid", sample_ids)
missing_columns <- setdiff(required_columns, colnames(dat))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

dat <- dat[dat$RNA_type %in% six_types, required_columns, drop = FALSE]
if (!setequal(unique(dat$RNA_type), six_types)) {
  stop("Not all six expected RNA types are present after filtering.")
}
if (anyNA(dat$geneid) || any(dat$geneid == "") || anyDuplicated(dat$geneid)) {
  stop("geneid contains missing, empty, or duplicated values.")
}

counts <- as.matrix(dat[, sample_ids, drop = FALSE])
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) {
  stop("Count matrix must contain non-negative integer values without missing data.")
}
rownames(counts) <- dat$geneid

keep <- rowSums(counts >= 1) >= min_count_samples
if (sum(keep) == 0L) {
  stop("No features remain after the expression filter.")
}
counts_keep <- counts[keep, , drop = FALSE]
type_map <- setNames(dat$RNA_type, dat$geneid)

coldata <- data.frame(
  sample_id = sample_ids,
  group = factor(c(rep("N", 10), rep("OP", 10)), levels = c("N", "OP")),
  row.names = sample_ids
)
if (!identical(colnames(counts_keep), rownames(coldata))) {
  stop("Count matrix columns do not match metadata rows.")
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_keep,
  colData = coldata,
  design = ~ group
)
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
res <- results(dds, contrast = c("group", "OP", "N"), alpha = fdr_cutoff)

res_df <- as.data.frame(res)
res_df$feature_id <- rownames(res_df)
res_df$RNA_type <- unname(type_map[res_df$feature_id])
res_df$contrast <- "OP_vs_N"
res_df$case_group <- "OP"
res_df$reference_group <- "N"
res_df$n_OP <- sum(coldata$group == "OP")
res_df$n_N <- sum(coldata$group == "N")
res_df$n_features_before_filter <- nrow(counts)
res_df$n_features_after_filter <- nrow(counts_keep)
res_df$status <- "Not significant"

is_sig <- !is.na(res_df$padj) &
  res_df$padj < fdr_cutoff &
  is.finite(res_df$log2FoldChange) &
  abs(res_df$log2FoldChange) > lfc_cutoff
res_df$status[is_sig & res_df$log2FoldChange > lfc_cutoff] <- "Upregulated"
res_df$status[is_sig & res_df$log2FoldChange < -lfc_cutoff] <- "Downregulated"

res_df <- res_df[, c(
  "contrast", "feature_id", "RNA_type", "baseMean", "log2FoldChange",
  "lfcSE", "stat", "pvalue", "padj", "status", "case_group",
  "reference_group", "n_OP", "n_N", "n_features_before_filter",
  "n_features_after_filter"
)]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), , drop = FALSE]

sig_df <- res_df[res_df$status %in% c("Upregulated", "Downregulated"), , drop = FALSE]
summary_by_type <- res_df %>%
  group_by(RNA_type) %>%
  summarise(
    total_features = n(),
    features_with_padj = sum(!is.na(padj)),
    upregulated = sum(status == "Upregulated"),
    downregulated = sum(status == "Downregulated"),
    not_significant = sum(status == "Not significant"),
    .groups = "drop"
  )
summary_all <- res_df %>%
  summarise(
    RNA_type = "All",
    total_features = n(),
    features_with_padj = sum(!is.na(padj)),
    upregulated = sum(status == "Upregulated"),
    downregulated = sum(status == "Downregulated"),
    not_significant = sum(status == "Not significant")
  )
summary_df <- bind_rows(summary_all, summary_by_type)

write.csv(res_df, file.path(table_dir, "DE_all.csv"), row.names = FALSE, quote = FALSE)
write.csv(sig_df, file.path(table_dir, "DE_sig.csv"), row.names = FALSE, quote = FALSE)
write.csv(summary_df, file.path(table_dir, "DE_summary.csv"), row.names = FALSE, quote = FALSE)
saveRDS(list(dds = dds, result = res_df, metadata = coldata), file.path(rds_dir, "DE_objects.rds"))

plot_df <- res_df[
  is.finite(res_df$log2FoldChange) & !is.na(res_df$padj),
  ,
  drop = FALSE
]
positive_fdr <- plot_df$padj[is.finite(plot_df$padj) & plot_df$padj > 0]
if (length(positive_fdr) == 0L) {
  stop("No positive adjusted P values are available for the volcano plot.")
}
plot_df$plot_fdr <- plot_df$padj
plot_df$plot_fdr[plot_df$plot_fdr <= 0 | !is.finite(plot_df$plot_fdr)] <- min(positive_fdr) / 10
plot_df$neg_log10_fdr <- -log10(plot_df$plot_fdr)
plot_df$status <- factor(
  plot_df$status,
  levels = c("Upregulated", "Not significant", "Downregulated")
)

status_counts <- table(plot_df$status)
legend_labels <- c(
  "Upregulated" = sprintf("Upregulated (%d)", status_counts[["Upregulated"]]),
  "Not significant" = sprintf("Not Significant (%d)", status_counts[["Not significant"]]),
  "Downregulated" = sprintf("Downregulated (%d)", status_counts[["Downregulated"]])
)
volcano_colors <- c(
  "Upregulated" = "#870F45",
  "Not significant" = "#D2DAE2",
  "Downregulated" = "#46559F"
)
xmax <- ceiling(max(abs(plot_df$log2FoldChange), lfc_cutoff * 1.5) * 1.08 * 2) / 2
ymax <- ceiling(max(plot_df$neg_log10_fdr, -log10(fdr_cutoff) * 1.5) * 1.08)

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_fdr)) +
  geom_point(
    data = plot_df[plot_df$status == "Not significant", , drop = FALSE],
    aes(colour = status), shape = 19, alpha = 0.50, size = 1
  ) +
  geom_point(
    data = plot_df[plot_df$status == "Upregulated", , drop = FALSE],
    aes(colour = status), shape = 19, alpha = 0.90, size = 1.5
  ) +
  geom_point(
    data = plot_df[plot_df$status == "Downregulated", , drop = FALSE],
    aes(colour = status), shape = 19, alpha = 0.90, size = 1.5
  ) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = 4, colour = "black", linewidth = 0.4) +
  scale_colour_manual(values = volcano_colors, labels = legend_labels, drop = FALSE) +
  coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
  labs(
    title = "OP vs N",
    subtitle = "Adjusted P value < 0.05 and |log2FC| > 0.5",
    x = "log2(Fold change, OP vs N)",
    y = "-log10(Adjusted P value)",
    colour = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = grid::unit(0.15, "cm"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

ggsave(file.path(figure_dir, "DE_volcano.pdf"), p_volcano, width = 7, height = 6, device = cairo_pdf)
ggsave(file.path(figure_dir, "DE_volcano.png"), p_volcano, width = 7, height = 6, dpi = 300, bg = "white")
