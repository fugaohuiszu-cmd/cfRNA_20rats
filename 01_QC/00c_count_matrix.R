args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
call_dir <- args[[1]]
output_dir <- args[[2]]
pairs_file <- args[[3]]
sample_ids <- c(paste0("N", 1:10), paste0("OP", 1:10))
rna_classes <- c("mRNA", "lncRNA", "miRNA", "piRNA", "snRNA", "snoRNA")

pairs <- read.csv(pairs_file, stringsAsFactors = FALSE)
stopifnot("SampleID_run1" %in% names(pairs))
stopifnot(nrow(pairs) == 20L, !anyDuplicated(pairs$clinicalID),
          setequal(pairs$clinicalID, sample_ids),
          identical(pairs$clinicalID, pairs$clinicalID2))
pairs <- pairs[match(sample_ids, pairs$clinicalID), ]
sample_names <- paste0(sample_ids, "_", sub("^.*_", "", pairs$SampleID_run1))

read_small_rna <- function(rna_class) {
  paths <- file.path(call_dir, rna_class, paste0(sample_names, "_", rna_class, ".txt"))
  stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0))
  tables <- lapply(paths, function(path) read.table(path, stringsAsFactors = FALSE))
  stopifnot(all(vapply(tables, function(x) ncol(x) == 2L && !anyDuplicated(x[[2]]), logical(1))))
  gene_ids <- sort(unique(unlist(lapply(tables, function(x) x[[2]]))))
  counts <- matrix(0, nrow = length(gene_ids), ncol = length(sample_ids),
                   dimnames = list(gene_ids, sample_ids))
  for (i in seq_along(tables)) {
    counts[match(tables[[i]][[2]], gene_ids), i] <- tables[[i]][[1]]
  }
  counts
}

read_feature_counts <- function(rna_class) {
  path <- file.path(call_dir, rna_class, paste0(rna_class, "_S1_transcript.txt"))
  counts <- read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(names(counts)[1:6], c("Geneid", "Chr", "Start", "End", "Strand", "Length")),
            !anyDuplicated(counts$Geneid))
  sample_columns <- sub("_map2GRCr8\\.sam$", "", basename(names(counts)[-(1:6)]))
  stopifnot(!anyDuplicated(sample_columns), setequal(sample_columns, sample_names))
  matrix <- as.matrix(counts[, match(sample_names, sample_columns) + 6L, drop = FALSE])
  dimnames(matrix) <- list(counts$Geneid, sample_ids)
  matrix
}

matrices <- setNames(lapply(rna_classes, function(rna_class) {
  if (rna_class %in% c("miRNA", "piRNA")) {
    read_small_rna(rna_class)
  } else {
    read_feature_counts(rna_class)
  }
}), rna_classes)
stopifnot(all(vapply(matrices, function(x) all(is.finite(x) & x >= 0 & x == floor(x)), logical(1))))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (rna_class in rna_classes) {
  write.csv(matrices[[rna_class]], file.path(output_dir, paste0(rna_class, ".csv")))
}
combined <- do.call(rbind, lapply(rna_classes, function(rna_class) {
  data.frame(RNA_type = rna_class, geneid = rownames(matrices[[rna_class]]),
             matrices[[rna_class]], row.names = NULL, check.names = FALSE)
}))
stopifnot(!anyDuplicated(combined$geneid))
write.csv(combined, file.path(output_dir, "all_RNA_final_merged.csv"))
cat(nrow(combined), "features;", length(sample_ids), "samples\n")
