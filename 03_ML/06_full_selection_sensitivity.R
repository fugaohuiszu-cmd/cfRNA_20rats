options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(DESeq2)
  library(glmnet)
  library(Boruta)
  library(ggplot2)
})

project_dir <- "."
output_dir <- file.path(project_dir, "03_results", "03_ML", "03_selection_sensitivity")
raw_file <- file.path(project_dir, "03_results/02_DE/rds/DE_objects.rds")
selection_file <- file.path(project_dir, "03_results/03_ML/02_strict_model/01_feature_selection/rds/strict_feature_selection.rds")
seed <- 20260717L
boruta_seed_offset <- 100000L
stability_iterations <- 100L
within_subset_cutoff <- 0.10
frequency_cutoffs <- c(0.15, 0.20, 0.25)
fdr_cutoff <- 0.05
lfc_cutoff <- 0.50
lasso_cv_folds <- 4L
boruta_pvalue <- 0.001
boruta_max_runs <- 300L
boruta_num_trees <- 300L
boruta_threads <- 4L
workers <- 8L

obj <- readRDS(raw_file)
selection <- readRDS(selection_file)
stopifnot(all(c("dds", "metadata", "result") %in% names(obj)),
          all(c("train_ids", "test_ids", "final_features") %in% names(selection)))
raw_counts <- counts(obj$dds, normalized = FALSE)
metadata <- obj$metadata
train_ids <- selection$train_ids
test_ids <- selection$test_ids
original_features <- selection$final_features$feature_id
stopifnot(identical(rownames(metadata), colnames(raw_counts)), length(train_ids) == 12L,
          length(test_ids) == 8L, !any(train_ids %in% test_ids), length(original_features) == 5L)
stopifnot(identical(levels(metadata$group), c("N", "OP")),
          identical(as.integer(table(metadata[train_ids, "group"])), c(6L, 6L)))
type_map <- setNames(obj$result$RNA_type, obj$result$feature_id)
stopifnot(!anyDuplicated(names(type_map)), all(rownames(raw_counts) %in% names(type_map)))
n_ids <- train_ids[metadata[train_ids, "group"] == "N"]
op_ids <- train_ids[metadata[train_ids, "group"] == "OP"]
sampling <- expand.grid(omitted_N = n_ids, omitted_OP = op_ids, stringsAsFactors = FALSE)
sampling$subsample <- seq_len(nrow(sampling))
sampling$included_samples <- vapply(seq_len(nrow(sampling)), function(i) {
  paste(setdiff(train_ids, c(sampling$omitted_N[i], sampling$omitted_OP[i])), collapse = ";")
}, character(1))
sampling <- sampling[, c("subsample", "omitted_N", "omitted_OP", "included_samples")]
stopifnot(nrow(sampling) == 36L, !anyDuplicated(sampling$included_samples))
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir, recursive = TRUE)
dir.create(file.path(output_dir, "progress"))
dir.create(file.path(output_dir, "subsamples"))
write.csv(sampling, file.path(output_dir, "subsamples.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = c(train_ids, test_ids), dataset = c(rep("Training", 12), rep("Internal hold-out", 8))),
          file.path(output_dir, "sample_split.csv"), row.names = FALSE)

make_foldid <- function(group, nfolds) {
  if (min(table(group)) < nfolds) stop("Too few animals per group for LASSO cross-validation.")
  foldid <- integer(length(group))
  for (g in levels(group)) {
    idx <- which(group == g)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  foldid
}

run_subsample <- function(i) {
  stop_file <- file.path(output_dir, "STOP.txt")
  if (file.exists(stop_file)) stop("Another worker reported an error.")
  withCallingHandlers({
    start <- proc.time()[[3]]
    ids <- setdiff(train_ids, c(sampling$omitted_N[i], sampling$omitted_OP[i]))
    meta <- metadata[ids, , drop = FALSE]
    stopifnot(length(ids) == 10L, !any(ids %in% test_ids), identical(as.integer(table(meta$group)), c(5L, 5L)))
    progress_file <- file.path(output_dir, "progress", sprintf("subsample_%02d.csv", i))
    write_progress <- function(k, lasso_n, boruta_n, final_n) {
      write.csv(data.frame(subsample = i, completed_runs = k, total_runs = stability_iterations,
                           candidates = ncol(x), LASSO_features = lasso_n, Boruta_features = boruta_n,
                           shared_features = final_n, elapsed_seconds = proc.time()[[3]] - start),
                progress_file, row.names = FALSE)
    }

    x_counts <- raw_counts[, ids, drop = FALSE]
    x_counts <- x_counts[rowSums(x_counts >= 1L) >= 1L, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = x_counts, colData = meta, design = ~ group)
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    de <- as.data.frame(results(dds, contrast = c("group", "OP", "N"), alpha = fdr_cutoff))
    de$Geneid <- rownames(de)
    candidate <- de[!is.na(de$padj) & de$padj < fdr_cutoff & is.finite(de$log2FoldChange) & abs(de$log2FoldChange) > lfc_cutoff, , drop = FALSE]
    candidate <- candidate[order(candidate$padj, -abs(candidate$log2FoldChange)), , drop = FALSE]
    if (nrow(candidate) < 2L) stop("Fewer than two DE candidates in subsample ", i)
    vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
    x <- t(assay(vst)[candidate$Geneid, ids, drop = FALSE])
    stopifnot(all(is.finite(x)))
    x <- x[, apply(x, 2, sd) > 1e-8, drop = FALSE]
    candidate <- candidate[match(colnames(x), candidate$Geneid), , drop = FALSE]
    if (ncol(x) < 2L) stop("Fewer than two variable candidates in subsample ", i)
    genes <- colnames(x)
    lasso_hits <- setNames(integer(length(genes)), genes)
    boruta_hits <- setNames(integer(length(genes)), genes)
    lasso_abs_coef <- setNames(numeric(length(genes)), genes)
    runs <- vector("list", stability_iterations)
    warning_counts <- c(small_class = 0L, grouped_FALSE = 0L, no_tentative = 0L)
    write_progress(0L, 0L, 0L, 0L)

    for (k in seq_len(stability_iterations)) {
      if (file.exists(stop_file)) stop("Another worker reported an error.")
      runs[[k]] <- withCallingHandlers({
        set.seed(seed + k)
        foldid <- make_foldid(meta$group, lasso_cv_folds)
        lasso <- cv.glmnet(x = x, y = as.integer(meta$group == "OP"), family = "binomial",
                          alpha = 1, type.measure = "deviance", nfolds = lasso_cv_folds,
                          foldid = foldid, standardize = TRUE)
        stopifnot(all(is.finite(lasso$cvm)), is.finite(lasso$lambda.1se))
        coefs <- as.matrix(coef(lasso, s = "lambda.1se"))
        coefs <- coefs[rownames(coefs) != "(Intercept)", 1]
        selected_lasso <- names(coefs)[abs(coefs) > 0]
        lasso_hits[selected_lasso] <- lasso_hits[selected_lasso] + 1L
        lasso_abs_coef[names(coefs)] <- lasso_abs_coef[names(coefs)] + abs(coefs)

        set.seed(seed + boruta_seed_offset + k)
        boruta <- Boruta(x = as.data.frame(x, check.names = FALSE), y = meta$group,
                         getImp = getImpRfGini, pValue = boruta_pvalue, maxRuns = boruta_max_runs,
                         doTrace = 0, num.trees = boruta_num_trees, num.threads = boruta_threads)
        boruta <- TentativeRoughFix(boruta)
        bstats <- attStats(boruta)
        if (any(!is.finite(bstats$meanImp))) stop("Non-finite Boruta importance in subsample ", i, ", run ", k)
        selected_boruta <- rownames(bstats)[bstats$decision == "Confirmed"]
        boruta_hits[selected_boruta] <- boruta_hits[selected_boruta] + 1L
        list(run = k, LASSO = selected_lasso, Boruta = selected_boruta,
             lambda_1se = lasso$lambda.1se, foldid = setNames(foldid, ids))
      }, warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("one multinomial or binomial class has fewer than 8", msg, fixed = TRUE)) {
          warning_counts["small_class"] <<- warning_counts["small_class"] + 1L
        } else if (grepl("Option grouped=FALSE enforced in cv.glmnet", msg, fixed = TRUE)) {
          warning_counts["grouped_FALSE"] <<- warning_counts["grouped_FALSE"] + 1L
        } else if (grepl("There are no Tentative attributes! Returning original object.", msg, fixed = TRUE)) {
          warning_counts["no_tentative"] <<- warning_counts["no_tentative"] + 1L
        } else {
          stop("Subsample ", i, ", run ", k, ": ", msg, call. = FALSE)
        }
        invokeRestart("muffleWarning")
      })
      if (k %% 10L == 0L) {
        l <- lasso_hits / k >= within_subset_cutoff
        b <- boruta_hits / k >= within_subset_cutoff
        write_progress(k, sum(l), sum(b), sum(l & b))
      }
    }

    records <- data.frame(subsample = i, Geneid = genes, RNA_type = unname(type_map[genes]),
                          log2FoldChange = candidate$log2FoldChange, padj = candidate$padj,
                          LASSO_count = as.integer(lasso_hits[genes]), Boruta_count = as.integer(boruta_hits[genes]),
                          LASSO_frequency = as.numeric(lasso_hits[genes]) / stability_iterations,
                          Boruta_frequency = as.numeric(boruta_hits[genes]) / stability_iterations,
                          LASSO_mean_absolute_coefficient = as.numeric(lasso_abs_coef[genes]) / stability_iterations)
    records$LASSO_selected <- records$LASSO_frequency >= within_subset_cutoff
    records$Boruta_selected <- records$Boruta_frequency >= within_subset_cutoff
    records$final_selected <- records$LASSO_selected & records$Boruta_selected
    stopifnot(!anyNA(records), all(records$LASSO_count <= 100L), all(records$Boruta_count <= 100L))
    result <- list(subsample = i, samples = ids, records = records, runs = runs,
                   warnings = warning_counts, elapsed_seconds = proc.time()[[3]] - start)
    saveRDS(result, file.path(output_dir, "subsamples", sprintf("subsample_%02d.rds", i)))
    write.csv(records, file.path(output_dir, "subsamples", sprintf("subsample_%02d_features.csv", i)), row.names = FALSE)
    result
  }, error = function(e) {
    if (!file.exists(stop_file)) writeLines(paste("Subsample", i, conditionMessage(e)), stop_file, useBytes = TRUE)
  })
}

run_all <- function() {
  cl <- parallel::makePSOCKcluster(workers, outfile = file.path(output_dir, "workers.log"))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("raw_counts", "metadata", "train_ids", "test_ids", "type_map", "sampling",
                               "output_dir", "seed", "boruta_seed_offset", "stability_iterations", "within_subset_cutoff",
                               "fdr_cutoff", "lfc_cutoff", "lasso_cv_folds", "boruta_pvalue", "boruta_max_runs",
                               "boruta_num_trees", "boruta_threads", "make_foldid", "run_subsample"), envir = .GlobalEnv)
  parallel::clusterEvalQ(cl, {
    options(stringsAsFactors = FALSE, warn = 1)
    suppressPackageStartupMessages({library(DESeq2); library(glmnet); library(Boruta)})
    NULL
  })
  parallel::parLapplyLB(cl, seq_len(nrow(sampling)), run_subsample)
}

results <- run_all()
stopifnot(length(results) == 36L, all(vapply(results, function(x) length(x$runs), integer(1)) == 100L))
records <- do.call(rbind, lapply(results, `[[`, "records"))
all_ids <- unique(c(original_features, records$Geneid))
count_ids <- function(ids) as.integer(table(factor(ids, levels = all_ids)))
frequency <- data.frame(Geneid = all_ids, RNA_type = unname(type_map[all_ids]),
                        DE_count = count_ids(records$Geneid),
                        LASSO_count = count_ids(records$Geneid[records$LASSO_selected]),
                        Boruta_count = count_ids(records$Geneid[records$Boruta_selected]),
                        final_count = count_ids(records$Geneid[records$final_selected]))
frequency$DE_frequency <- frequency$DE_count / 36
frequency$LASSO_frequency <- frequency$LASSO_count / 36
frequency$Boruta_frequency <- frequency$Boruta_count / 36
frequency$final_frequency <- frequency$final_count / 36
frequency$original_model_feature <- frequency$Geneid %in% original_features
frequency <- frequency[order(-as.integer(frequency$original_model_feature), -frequency$final_frequency,
                             -frequency$LASSO_frequency - frequency$Boruta_frequency, frequency$Geneid), ]
stopifnot(all(frequency$final_count <= pmin(frequency$LASSO_count, frequency$Boruta_count)),
          all(frequency$LASSO_count <= frequency$DE_count), all(frequency$Boruta_count <= frequency$DE_count),
          all(frequency$DE_count <= 36L), !anyNA(frequency))
threshold_summary <- do.call(rbind, lapply(frequency_cutoffs, function(cutoff) {
  keep <- frequency$final_frequency >= cutoff
  data.frame(frequency_threshold = cutoff, minimum_selections = ceiling(cutoff * 36),
             LASSO_features = sum(frequency$LASSO_frequency >= cutoff),
             Boruta_features = sum(frequency$Boruta_frequency >= cutoff),
             shared_features = sum(keep), original_features_retained = sum(keep & frequency$original_model_feature))
}))
threshold_features <- do.call(rbind, lapply(frequency_cutoffs, function(cutoff) {
  keep <- frequency$final_frequency >= cutoff
  data.frame(frequency_threshold = rep(cutoff, sum(keep)), frequency[keep, , drop = FALSE])
}))
original_frequency <- frequency[match(original_features, frequency$Geneid), ]
progress <- do.call(rbind, lapply(results, function(x) data.frame(subsample = x$subsample, runs = length(x$runs),
  candidates = nrow(x$records), LASSO_selected = sum(x$records$LASSO_selected),
  Boruta_selected = sum(x$records$Boruta_selected), shared_features = sum(x$records$final_selected),
  elapsed_seconds = x$elapsed_seconds)))
warning_summary <- do.call(rbind, lapply(results, function(x) data.frame(subsample = x$subsample,
  small_class = x$warnings["small_class"], grouped_FALSE = x$warnings["grouped_FALSE"], no_tentative = x$warnings["no_tentative"])))
write.csv(records, file.path(output_dir, "subsample_feature_selection.csv"), row.names = FALSE)
write.csv(frequency, file.path(output_dir, "feature_selection_frequency.csv"), row.names = FALSE)
write.csv(original_frequency, file.path(output_dir, "original_five_feature_frequency.csv"), row.names = FALSE)
write.csv(threshold_summary, file.path(output_dir, "threshold_sensitivity.csv"), row.names = FALSE)
write.csv(threshold_features, file.path(output_dir, "features_by_threshold.csv"), row.names = FALSE)
write.csv(progress, file.path(output_dir, "subsample_results.csv"), row.names = FALSE)
write.csv(warning_summary, file.path(output_dir, "warnings.csv"), row.names = FALSE)
parameters <- data.frame(parameter = c("seed", "subsamples", "animals_per_group", "runs_per_subsample", "within_subset_cutoff",
  "DE_padj", "DE_abs_log2FC", "LASSO_folds", "Boruta_pvalue", "Boruta_maxRuns", "Boruta_num_trees", "Boruta_threads", "parallel_workers"),
  value = c(seed, 36, 5, stability_iterations, within_subset_cutoff, fdr_cutoff, lfc_cutoff, lasso_cv_folds,
            boruta_pvalue, boruta_max_runs, boruta_num_trees, boruta_threads, workers))
write.csv(parameters, file.path(output_dir, "parameters.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(list(frequency = frequency, original_features = original_frequency, thresholds = threshold_summary,
                         threshold_features = threshold_features, subsamples = sampling, progress = progress, parameters = parameters),
                    file.path(output_dir, "results.json"), pretty = TRUE, dataframe = "rows", digits = 12)

plot_data <- original_frequency
plot_data$Geneid <- factor(plot_data$Geneid, levels = rev(original_features))
p <- ggplot(plot_data, aes(x = Geneid, y = final_frequency)) + geom_col(width = 0.65, fill = "#3B6D9A") +
  coord_flip() + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::label_percent()) +
  labs(x = NULL, y = "Final selection frequency across 36 subsets") + theme_classic(base_size = 11)
ggsave(file.path(output_dir, "feature_selection_frequency.pdf"), p, width = 8, height = 4)
ggsave(file.path(output_dir, "feature_selection_frequency.png"), p, width = 8, height = 4, dpi = 300)
