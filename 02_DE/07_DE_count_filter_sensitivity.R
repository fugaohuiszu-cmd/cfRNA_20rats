options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(DESeq2))

project_dir <- "."
table_dir <- file.path(project_dir, "03_results", "02_DE", "tables")
out_dir <- commandArgs(trailingOnly = TRUE)
stopifnot(length(out_dir) == 1L, dir.exists(out_dir))

classes <- c("mRNA", "lncRNA", "miRNA", "piRNA", "snRNA", "snoRNA")
animals <- c(paste0("N", 1:10), paste0("OP", 1:10))
min_count <- 10L
min_animals <- 5L
alpha <- 0.05
min_lfc <- 0.5
panel <- c("ENSRNOG00000000047", "ENSRNOG00000036688", "ENSRNOG00000054973", "piR-rno-2281228", "piR-rno-2096851")
files <- file.path(out_dir, paste0("DE_count_filter_", c("comparison", "summary", "five_features"), ".csv"))

raw <- read.csv(file.path(project_dir, "01_input", "all_RNA_final_merged.csv"), check.names = FALSE)
raw <- raw[raw$RNA_type %in% classes, c("RNA_type", "geneid", animals)]
stopifnot(nrow(raw) == 44038L, !anyNA(raw), !anyDuplicated(raw$geneid), all(nzchar(raw$geneid)))
m <- as.matrix(raw[, animals])
stopifnot(all(m >= 0), all(m == floor(m)))
rownames(m) <- raw$geneid
metadata <- data.frame(group = factor(rep(c("N","OP"), each = 10), levels = c("N","OP")), row.names = animals)
stopifnot(identical(colnames(m), rownames(metadata)))
original_keep <- rowSums(m >= 1L) >= 1L
n_detected <- rowSums(m >= min_count)
strict_keep <- n_detected >= min_animals
original <- read.csv(file.path(table_dir, "DE_all.csv"), check.names = FALSE)
stopifnot(setequal(original$feature_id, raw$geneid[original_keep]), !anyDuplicated(original$feature_id),
          all(original$contrast == "OP_vs_N"), all(original$RNA_type == raw$RNA_type[match(original$feature_id, raw$geneid)]))

dds <- DESeqDataSetFromMatrix(countData = m[strict_keep, , drop = FALSE], colData = metadata, design = ~ group)
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("group", "OP", "N"), alpha = alpha))
stopifnot(identical(rownames(res), rownames(m)[strict_keep]))

status <- function(keep, p, q, fc) {
  z <- rep("Not significant", length(keep))
  z[keep & is.na(q)] <- "No adjusted P"
  z[keep & is.na(p)] <- "No P value"
  z[!keep] <- "Below count filter"
  pass <- keep & !is.na(q) & q < alpha & is.finite(fc)
  z[pass & fc > min_lfc] <- "Upregulated"
  z[pass & fc < -min_lfc] <- "Downregulated"
  z
}
old_idx <- match(raw$geneid, original$feature_id)
new_idx <- match(raw$geneid, rownames(res))
x <- data.frame(Geneid = raw$geneid, RNA_type = raw$RNA_type, rats_with_at_least_10_reads = as.integer(n_detected),
                original_filter_pass = original_keep, strict_filter_pass = strict_keep,
                original_log2FC = original$log2FoldChange[old_idx], original_P = original$pvalue[old_idx],
                original_adjusted_P = original$padj[old_idx], strict_log2FC = res$log2FoldChange[new_idx],
                strict_P = res$pvalue[new_idx], strict_adjusted_P = res$padj[new_idx])
x$original_status <- with(x, status(original_filter_pass, original_P, original_adjusted_P, original_log2FC))
x$strict_status <- with(x, status(strict_filter_pass, strict_P, strict_adjusted_P, strict_log2FC))
summarise <- function(d, cls) {
  a <- d$original_status
  b <- d$strict_status
  sig <- c("Upregulated", "Downregulated")
  old <- a %in% sig
  new <- b %in% sig
  data.frame(RNA_type = cls, raw_features = nrow(d), original_kept = sum(d$original_filter_pass),
             strict_kept = sum(d$strict_filter_pass), strict_with_adjusted_P = sum(!is.na(d$strict_adjusted_P)),
             original_up = sum(a == "Upregulated"), original_down = sum(a == "Downregulated"),
             strict_up = sum(b == "Upregulated"), strict_down = sum(b == "Downregulated"),
             significant_both_same_direction = sum(old & new & a == b),
             significant_both_opposite_direction = sum(old & new & a != b),
             original_significant_below_strict_filter = sum(old & !d$strict_filter_pass),
             original_significant_retained_no_adjusted_P = sum(old & d$strict_filter_pass & is.na(d$strict_adjusted_P)),
             original_significant_retained_not_significant = sum(old & d$strict_filter_pass & !is.na(d$strict_adjusted_P) & !new),
             significant_only_after_strict_filter = sum(!old & new))
}
totals <- do.call(rbind, lapply(classes, function(k) summarise(x[x$RNA_type == k, ], k)))
totals <- rbind(totals, summarise(x, "All"))
baseline <- read.csv(file.path(table_dir, "DE_summary.csv"))
idx <- match(totals$RNA_type, baseline$RNA_type)
stopifnot(all(totals$original_up == baseline$upregulated[idx]), all(totals$original_down == baseline$downregulated[idx]),
          all(totals$original_kept == baseline$total_features[idx]),
          with(totals, all(original_up + original_down == significant_both_same_direction + significant_both_opposite_direction +
            original_significant_below_strict_filter + original_significant_retained_no_adjusted_P + original_significant_retained_not_significant)),
          with(totals, all(strict_up + strict_down == significant_both_same_direction + significant_both_opposite_direction + significant_only_after_strict_filter)))
panel_results <- x[match(panel, x$Geneid), ]
stopifnot(identical(panel_results$Geneid, panel), !any(file.exists(files)))
write.csv(x, files[1], row.names = FALSE, na = "NA")
write.csv(totals, files[2], row.names = FALSE, na = "NA")
write.csv(panel_results, files[3], row.names = FALSE, na = "NA")
print(totals, row.names = FALSE)
print(panel_results[, c("Geneid", "rats_with_at_least_10_reads", "original_log2FC", "strict_log2FC", "strict_adjusted_P", "strict_status")], row.names = FALSE)
