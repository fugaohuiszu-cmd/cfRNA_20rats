project_dir <- "F:/DE-analysis/final/cfRNA_project"
result_dir <- file.path(project_dir, "03_results/04_source_tracing/female_reference")
output_dir <- file.path(result_dir, "background_comparison")
p_cutoff <- 0.05
log2fc_cutoff <- 0.5

de <- read.csv(file.path(project_dir, "03_results/02_DE/tables/DE_all.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
reference <- read.delim(file.path(result_dir, "female_reference_FPKM.tsv"),
                        row.names = 1, check.names = FALSE)
assignment <- read.csv(file.path(result_dir, "tissue_assignment.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(reference) == 44, !anyDuplicated(de$feature_id))
mrna <- de[de$contrast == "OP_vs_N" & de$RNA_type == "mRNA", ]
tested <- mrna[is.finite(mrna$pvalue) & is.finite(mrna$log2FoldChange), ]
matched <- tested[tested$feature_id %in% rownames(reference), ]
expression <- as.matrix(reference[matched$feature_id, ])
stopifnot(!anyNA(expression), all(expression >= 0))
peak <- apply(expression, 1, max)
peak_count <- rowSums(expression == peak)
eligible <- peak > 0 & peak_count == 1L
tissue <- rep(NA_character_, nrow(expression))
tissue[peak > 0] <- apply(expression[peak > 0, ], 1, function(x) {
  paste(colnames(expression)[x == max(x)], collapse = "; ")
})
background <- matched[eligible, ]
background$tissue <- tissue[eligible]
stopifnot(nrow(background) == 16719L, all(assignment$Geneid %in% background$feature_id))

background_reference <- data.frame(
  "Geneid" = matched$feature_id,
  "Highest-expression tissue" = tissue,
  "Highest FPKM" = peak,
  "log2FoldChange" = matched$log2FoldChange,
  "P value" = matched$pvalue,
  check.names = FALSE
)
directions <- c(Upregulated = 1, Downregulated = -1)
results <- list()
selected_counts <- integer(2)
names(selected_counts) <- names(directions)

for (direction in names(directions)) {
  selected <- background$pvalue < p_cutoff &
    background$log2FoldChange * directions[[direction]] > log2fc_cutoff
  selected_counts[[direction]] <- sum(selected)
  existing <- assignment[assignment[["Expression change"]] == direction, ]
  stopifnot(setequal(background$feature_id[selected], existing$Geneid))
  stopifnot(identical(background$tissue[match(existing$Geneid, background$feature_id)],
                      existing[["Highest-expression tissue"]]))

  for (target in colnames(reference)) {
    in_tissue <- background$tissue == target
    observed <- sum(selected & in_tissue)
    tissue_total <- sum(in_tissue)
    selected_total <- sum(selected)
    background_total <- nrow(background)
    contingency <- matrix(c(
      observed, selected_total - observed,
      tissue_total - observed,
      background_total - selected_total - tissue_total + observed
    ), nrow = 2, byrow = TRUE)
    stopifnot(all(contingency >= 0), sum(contingency) == background_total)
    test <- fisher.test(contingency, alternative = "greater")
    expected <- selected_total * tissue_total / background_total
    results[[length(results) + 1L]] <- data.frame(
      "Expression change" = direction, "Tissue" = target,
      "Observed genes" = observed, "Background genes" = tissue_total,
      "Expected genes" = expected, "Fold enrichment" = observed / expected,
      "P value" = test$p.value, check.names = FALSE
    )
  }
}
results <- do.call(rbind, results)
stopifnot(nrow(results) == 88L)
results[["BH-adjusted P value"]] <- p.adjust(results[["P value"]], method = "BH")
results <- results[order(match(results[["Expression change"]], names(directions)),
                         results[["BH-adjusted P value"]], results[["P value"]]), ]

counts <- data.frame(
  "Stage" = c("Detected mRNAs", "mRNAs with valid differential-expression results",
               "Matched to female reference", "All-zero reference expression",
               "Tied highest reference expression", "Enrichment background",
               "Upregulated mRNAs in background", "Downregulated mRNAs in background"),
  "Genes (n)" = c(nrow(mrna), nrow(tested), nrow(matched), sum(peak == 0),
                   sum(peak > 0 & peak_count > 1), nrow(background), unname(selected_counts)),
  check.names = FALSE
)
stopifnot(!dir.exists(output_dir))
dir.create(output_dir)
write.csv(results, file.path(output_dir, "tissue_enrichment.csv"), row.names = FALSE)
write.csv(background_reference, file.path(output_dir, "background_reference.csv"),
          row.names = FALSE, na = "")
write.csv(counts, file.path(output_dir, "analysis_counts.csv"), row.names = FALSE)
print(results[results$Tissue %in% c("BoneMarrow", "Cartilage", "Liver", "WholeBlood"), ],
      row.names = FALSE, digits = 4)
cat("BH-adjusted P < 0.05:\n")
print(results[results[["BH-adjusted P value"]] < 0.05, ], row.names = FALSE, digits = 4)
cat("Results:", output_dir, "\n")
