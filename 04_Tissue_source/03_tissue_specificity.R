project_dir <- "."
result_dir <- file.path(project_dir, "03_results/04_source_tracing/female_reference")
output_dir <- file.path(result_dir, "tissue_specificity")

reference <- read.delim(file.path(result_dir, "female_reference_FPKM.tsv"),
                        row.names = 1, check.names = FALSE)
assignment <- read.csv(file.path(result_dir, "tissue_assignment.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
unassigned <- read.csv(file.path(result_dir, "unassigned_expression.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
de <- read.csv(file.path(project_dir, "03_results/02_DE/tables/DE_all.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
genes <- rbind(assignment[, c("Geneid", "Gene symbol", "Expression change")],
               unassigned[, c("Geneid", "Gene symbol", "Expression change")])
stopifnot(ncol(reference) == 44, nrow(genes) == 1784,
          !anyDuplicated(genes$Geneid), !anyDuplicated(de$feature_id),
          all(genes$Geneid %in% rownames(reference)),
          all(genes$Geneid %in% de$feature_id))

expression <- as.matrix(reference[genes$Geneid, ])
stopifnot(!anyNA(expression), all(expression >= 0))
log_expression <- log2(expression + 1)
peak <- apply(log_expression, 1, max)
positive <- peak > 0
stopifnot(setequal(genes$Geneid[!positive],
                   c("ENSRNOG00000057953", "ENSRNOG00000003227")))

tau <- rep(NA_real_, nrow(expression))
tau[positive] <- rowSums(1 - sweep(log_expression[positive, ], 1,
                                  peak[positive], "/")) / (ncol(expression) - 1)
stopifnot(all(tau[positive] >= 0 & tau[positive] <= 1))
top_tissue <- rep(NA_character_, nrow(expression))
top_tissue[positive] <- apply(expression[positive, ], 1, function(x) {
  paste(colnames(expression)[x == max(x)], collapse = "; ")
})
de_rows <- de[match(genes$Geneid, de$feature_id), ]
stopifnot(all(de_rows$contrast == "OP_vs_N"), all(de_rows$RNA_type == "mRNA"),
          all(de_rows$pvalue < 0.05), all(abs(de_rows$log2FoldChange) > 0.5))

results <- genes
results[["Highest-expression tissue"]] <- top_tissue
results[["Highest FPKM"]] <- apply(expression, 1, max)
results[["Tau"]] <- tau
results[["log2FoldChange"]] <- de_rows$log2FoldChange
results[["P value"]] <- de_rows$pvalue
results <- results[order(match(results[["Expression change"]],
                               c("Upregulated", "Downregulated")), results[["P value"]]), ]
bone <- results[which(results[["Highest-expression tissue"]] %in% c("BoneMarrow", "Cartilage")), ]
stopifnot(nrow(bone) == 68)
old_rows <- match(assignment$Geneid, results$Geneid)
stopifnot(identical(results[["Highest-expression tissue"]][old_rows],
                    assignment[["Highest-expression tissue"]]))

summary_rows <- list()
for (direction in c("Upregulated", "Downregulated")) {
  current <- results[results[["Expression change"]] == direction & !is.na(results$Tau), ]
  subsets <- list(
    "All evaluated mRNAs" = current,
    "Bone marrow" = current[which(current[["Highest-expression tissue"]] == "BoneMarrow"), ],
    "Cartilage" = current[which(current[["Highest-expression tissue"]] == "Cartilage"), ]
  )
  for (label in names(subsets)) {
    values <- subsets[[label]]$Tau
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      "Expression change" = direction, "Gene set" = label, "Genes (n)" = length(values),
      "Median Tau" = median(values), "Minimum Tau" = min(values),
      "Maximum Tau" = max(values), check.names = FALSE
    )
  }
}
summary_table <- do.call(rbind, summary_rows)
stopifnot(!dir.exists(output_dir))
dir.create(output_dir)
write.csv(results, file.path(output_dir, "tissue_specificity.csv"), row.names = FALSE, na = "")
write.csv(bone, file.path(output_dir, "bone_gene_specificity.csv"), row.names = FALSE, na = "")
write.csv(summary_table, file.path(output_dir, "tissue_specificity_summary.csv"), row.names = FALSE)
print(summary_table, row.names = FALSE, digits = 4)
cat("Results:", output_dir, "\n")
