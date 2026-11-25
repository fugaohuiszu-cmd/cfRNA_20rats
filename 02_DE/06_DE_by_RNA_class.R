options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(DESeq2))

project_dir <- "."
count_file <- file.path(project_dir, "01_input", "all_RNA_final_merged.csv")
table_dir <- file.path(project_dir, "03_results", "02_DE", "tables")
result_dir <- commandArgs(trailingOnly = TRUE)
stopifnot(length(result_dir) == 1L, dir.exists(result_dir))
outputs <- file.path(result_dir, paste0("DE_class_specific_", c("comparison", "summary", "five_features"), ".csv"))
stopifnot(!any(file.exists(outputs)))

classes <- c("mRNA", "lncRNA", "miRNA", "piRNA", "snRNA", "snoRNA")
sample_ids <- c(paste0("N", 1:10), paste0("OP", 1:10))
samples <- data.frame(group = factor(rep(c("N", "OP"), each = 10), levels = c("N", "OP")), row.names = sample_ids)
panel <- c("ENSRNOG00000000047", "ENSRNOG00000036688", "ENSRNOG00000054973", "piR-rno-2281228", "piR-rno-2096851")
fdr <- 0.05
lfc <- 0.5

raw <- read.csv(count_file, check.names = FALSE)
raw <- raw[raw$RNA_type %in% classes, c("RNA_type", "geneid", sample_ids)]
stopifnot(nrow(raw) == 44038L, !anyNA(raw), all(nzchar(raw$geneid)), !anyDuplicated(raw$geneid))
counts <- as.matrix(raw[, sample_ids])
stopifnot(all(counts >= 0), all(counts == floor(counts)), identical(colnames(counts), rownames(samples)))
rownames(counts) <- raw$geneid
keep <- rowSums(counts >= 1L) >= 1L
strata <- setNames(lapply(classes, function(k) counts[keep & raw$RNA_type == k, , drop = FALSE]), classes)

pooled <- read.csv(file.path(table_dir, "DE_all.csv"), check.names = FALSE)
stopifnot(setequal(pooled$feature_id, raw$geneid[keep]), !anyDuplicated(pooled$feature_id), all(pooled$contrast == "OP_vs_N"))
key <- match(pooled$feature_id, raw$geneid)
stopifnot(all(raw$RNA_type[key] == pooled$RNA_type), all(panel %in% pooled$feature_id))

fit_class <- function(m, samples, fdr) {
  dds <- DESeqDataSetFromMatrix(countData = m, colData = samples, design = ~ group)
  dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  de <- as.data.frame(results(dds, contrast = c("group", "OP", "N"), alpha = fdr))
  data.frame(Geneid = rownames(de), class_baseMean = de$baseMean,
             class_log2FC = de$log2FoldChange, class_P = de$pvalue, class_adjusted_P = de$padj)
}
fits <- lapply(strata, fit_class, samples = samples, fdr = fdr)
significance <- function(fc, p) {
  result <- rep("Not significant", length(p))
  pass <- is.finite(fc) & !is.na(p) & p < fdr
  result[pass & fc > lfc] <- "Upregulated"
  result[pass & fc < -lfc] <- "Downregulated"
  result
}
stopifnot(identical(significance(pooled$log2FoldChange, pooled$padj), pooled$status))
comparison <- do.call(rbind, lapply(classes, function(k) {
  fit <- fits[[k]]
  p <- pooled[match(fit$Geneid, pooled$feature_id), ]
  stopifnot(identical(p$feature_id, rownames(strata[[k]])), all(p$RNA_type == k))
  data.frame(RNA_type = k, Geneid = fit$Geneid, pooled_log2FC = p$log2FoldChange,
             pooled_P = p$pvalue, pooled_adjusted_P = p$padj, pooled_status = p$status,
             fit[, c("class_baseMean", "class_log2FC", "class_P", "class_adjusted_P")],
             class_status = significance(fit$class_log2FC, fit$class_adjusted_P),
             stringsAsFactors = FALSE)
}))
summarize_class <- function(x, k) {
  a <- x$pooled_status; b <- x$class_status
  sig_a <- a != "Not significant"; sig_b <- b != "Not significant"
  data.frame(RNA_type = k, features = nrow(x), pooled_nonmissing_adjusted_P = sum(!is.na(x$pooled_adjusted_P)),
             class_nonmissing_adjusted_P = sum(!is.na(x$class_adjusted_P)),
             pooled_up = sum(a == "Upregulated"), pooled_down = sum(a == "Downregulated"),
             class_up = sum(b == "Upregulated"), class_down = sum(b == "Downregulated"),
             significant_in_both_same_direction = sum(sig_a & sig_b & a == b),
             significant_in_both_opposite_directions = sum(sig_a & sig_b & a != b),
             significant_only_in_pooled = sum(sig_a & !sig_b), significant_only_in_class = sum(!sig_a & sig_b))
}
summary_table <- do.call(rbind, lapply(classes, function(k) summarize_class(comparison[comparison$RNA_type == k, ], k)))
summary_table <- rbind(summary_table, summarize_class(comparison, "All"))
orig <- read.csv(file.path(table_dir, "DE_summary.csv"))
idx <- match(summary_table$RNA_type, orig$RNA_type)
stopifnot(identical(summary_table$pooled_up, orig$upregulated[idx]), identical(summary_table$pooled_down, orig$downregulated[idx]),
          all(summary_table$features == orig$total_features[idx]), all(summary_table$pooled_nonmissing_adjusted_P == orig$features_with_padj[idx]))
panel_table <- comparison[match(panel, comparison$Geneid), ]
stopifnot(nrow(panel_table) == 5L, !anyNA(panel_table$Geneid), sum(keep) == nrow(comparison))
for (i in seq_along(outputs)) {
  write.csv(list(comparison, summary_table, panel_table)[[i]], outputs[[i]], row.names = FALSE, na = "NA")
}
print(summary_table, row.names = FALSE)
print(panel_table[, c("Geneid","pooled_log2FC","pooled_adjusted_P","class_log2FC","class_adjusted_P","class_status")], row.names = FALSE)
