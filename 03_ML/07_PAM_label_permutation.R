options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(DESeq2)
  library(glmnet)
  library(Boruta)
  library(caret)
  library(pamr)
})

project_dir <- "."
result_dir <- file.path(project_dir, "03_results", "03_ML", "02_models")
output_dir <- file.path(project_dir, "03_results/03_ML/03_label_permutation")
seed <- 20260717L
stability_iterations <- 100L
selection_cutoff <- 0.10
fdr_cutoff <- 0.05
lfc_cutoff <- 0.50
lasso_cv_folds <- 4L
boruta_seed_offset <- 100000L
boruta_pvalue <- 0.001
boruta_max_runs <- 300L
boruta_num_trees <- 300L
boruta_threads <- 4L
workers <- 24L

parameters <- read.csv(file.path(result_dir, "best_parameters.csv"))
pam_threshold <- as.numeric(parameters$Value[parameters$Model == "PAM" & parameters$Parameter == "threshold"])
stopifnot(length(pam_threshold) == 1L, is.finite(pam_threshold))
de <- readRDS(file.path(project_dir, "03_results/02_DE/rds/DE_objects.rds"))
selection <- readRDS(file.path(project_dir,
  "03_results/03_ML/02_strict_model/01_feature_selection/rds/strict_feature_selection.rds"))
stopifnot(all(c("dds", "metadata", "result") %in% names(de)),
          all(c("train_ids", "test_ids", "final_features") %in% names(selection)))
all_counts <- counts(de$dds, normalized = FALSE)
metadata <- de$metadata
train_ids <- selection$train_ids
test_ids <- selection$test_ids
type_map <- setNames(de$result$RNA_type, de$result$feature_id)
stopifnot(length(train_ids) == 12L, length(test_ids) == 8L,
          !any(train_ids %in% test_ids), identical(rownames(metadata), colnames(all_counts)),
          identical(levels(metadata$group), c("N", "OP")),
          identical(as.integer(table(metadata[train_ids, "group"])), c(6L, 6L)),
          identical(as.integer(table(metadata[test_ids, "group"])), c(4L, 4L)))
holdout_labels <- ifelse(metadata[test_ids, "group"] == "OP", "OVX", "Sham")
arrangements <- combn(seq_along(train_ids), 6L)
observed_positions <- which(metadata[train_ids, "group"] == "OP")
observed_id <- which(apply(arrangements, 2L, function(x) identical(as.integer(x), observed_positions)))
stopifnot(ncol(arrangements) == 924L, length(observed_id) == 1L, !dir.exists(output_dir))
dir.create(file.path(output_dir, "permutations"), recursive = TRUE)
write.csv(data.frame(permutation = seq_len(ncol(arrangements)),
  observed = seq_len(ncol(arrangements)) == observed_id,
  OVX_training_samples = apply(arrangements, 2L, function(x) paste(train_ids[x], collapse = ";"))),
  file.path(output_dir, "label_arrangements.csv"), row.names = FALSE)

make_foldid <- function(group, nfolds) {
  stopifnot(min(table(group)) >= nfolds)
  foldid <- integer(length(group))
  for (g in levels(group)) {
    idx <- which(group == g)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  foldid
}

initialize_lasso_check <- function() {
  stopifnot(as.character(packageVersion("glmnet")) == "4.1.10")
  invisible(trace("cv.glmnet.raw", where = asNamespace("glmnet"),
    exit = quote({
      if (exists("outlist", inherits = FALSE)) assign(".lasso_fold_paths", outlist, envir = .GlobalEnv)
    }), print = FALSE))
}

fit_lasso <- function(x, y, foldid) {
  convergence_messages <- character()
  assign(".lasso_fold_paths", NULL, envir = .GlobalEnv)
  fit <- withCallingHandlers(
    cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance",
               nfolds = lasso_cv_folds, foldid = foldid, standardize = TRUE),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (startsWith(msg, "from glmnet C++ code (error code -") &&
          grepl("Convergence for [0-9]+th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned$", msg)) {
        convergence_messages <<- c(convergence_messages, msg)
        invokeRestart("muffleWarning")
      }
    })
  stopifnot(all(is.finite(fit$cvm)), is.finite(fit$lambda.1se))
  if (length(convergence_messages) > 0L) {
    paths <- get(".lasso_fold_paths", envir = .GlobalEnv)
    stopifnot(fit$glmnet.fit$jerr == 0L, length(paths) == lasso_cv_folds)
    affected <- which(vapply(paths, function(p) p$jerr < 0L, logical(1)))
    stopifnot(length(affected) == length(convergence_messages))
    for (fold in affected) {
      path <- paths[[fold]]
      stopifnot(path$jerr > -10000L, length(path$lambda) >= 2L,
                all(is.finite(path$lambda)), all(is.finite(fit$lambda)),
                min(fit$lambda) >= min(path$lambda))
    }
  }
  fit
}

holdout_auc <- function(prob) {
  positive <- holdout_labels == "OVX"
  differences <- outer(prob[positive], prob[!positive], "-")
  mean((differences > 0) + 0.5 * (differences == 0))
}

run_permutation <- function(id) {
  withCallingHandlers({
    meta <- metadata[train_ids, , drop = FALSE]
    labels <- rep("N", 12L)
    labels[arrangements[, id]] <- "OP"
    meta$group <- factor(labels, levels = c("N", "OP"))
    train_counts <- all_counts[, train_ids, drop = FALSE]
    train_counts <- train_counts[rowSums(train_counts >= 1) >= 1L, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = train_counts, colData = meta, design = ~ group)
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("group", "OP", "N"), alpha = fdr_cutoff))
    res$Geneid <- rownames(res)
    candidates <- res[!is.na(res$padj) & res$padj < fdr_cutoff &
      is.finite(res$log2FoldChange) & abs(res$log2FoldChange) > lfc_cutoff, , drop = FALSE]
    candidates <- candidates[order(candidates$padj, -abs(candidates$log2FoldChange)), , drop = FALSE]
    n_candidates <- nrow(candidates)
    selected <- character()
    selection_records <- data.frame()
    selection_runs <- list()
    if (n_candidates >= 2L) {
      vst_train <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
      vst_matrix <- assay(vst_train)
      x <- t(vst_matrix[candidates$Geneid, train_ids, drop = FALSE])
      stopifnot(all(is.finite(x)))
      x <- x[, apply(x, 2L, sd) > 1e-8, drop = FALSE]
      candidates <- candidates[match(colnames(x), candidates$Geneid), , drop = FALSE]
      if (ncol(x) >= 2L) {
        genes <- colnames(x)
        lasso_hits <- boruta_hits <- setNames(integer(length(genes)), genes)
        selection_runs <- vector("list", stability_iterations)
        for (k in seq_len(stability_iterations)) {
          set.seed(seed + k)
          foldid <- make_foldid(meta$group, lasso_cv_folds)
          lasso <- fit_lasso(x, as.integer(meta$group == "OP"), foldid)
          coefs <- as.matrix(coef(lasso, s = "lambda.1se"))
          coefs <- coefs[rownames(coefs) != "(Intercept)", 1]
          lasso_selected <- names(coefs)[abs(coefs) > 0]
          lasso_hits[lasso_selected] <- lasso_hits[lasso_selected] + 1L
          set.seed(seed + boruta_seed_offset + k)
          boruta <- Boruta(as.data.frame(x, check.names = FALSE), meta$group, getImp = getImpRfGini,
            pValue = boruta_pvalue, maxRuns = boruta_max_runs, doTrace = 0,
            num.trees = boruta_num_trees, num.threads = boruta_threads)
          bstats <- attStats(TentativeRoughFix(boruta))
          stopifnot(all(is.finite(bstats$meanImp)))
          boruta_selected <- rownames(bstats)[bstats$decision == "Confirmed"]
          boruta_hits[boruta_selected] <- boruta_hits[boruta_selected] + 1L
          selection_runs[[k]] <- list(run = k, LASSO = lasso_selected, Boruta = boruta_selected,
            lambda_1se = lasso$lambda.1se, foldid = setNames(foldid, train_ids))
        }
        selection_records <- data.frame(Geneid = genes, RNA_type = unname(type_map[genes]),
          log2FoldChange = candidates$log2FoldChange, padj = candidates$padj,
          LASSO_count = as.integer(lasso_hits), Boruta_count = as.integer(boruta_hits))
        selection_records$final_selected <- lasso_hits / stability_iterations >= selection_cutoff &
          boruta_hits / stability_iterations >= selection_cutoff
        ordering <- order(-as.integer(selection_records$final_selected),
          -(selection_records$LASSO_count + selection_records$Boruta_count), selection_records$padj)
        selection_records <- selection_records[ordering, , drop = FALSE]
        selected <- selection_records$Geneid[selection_records$final_selected]
      }
    }
    probability <- setNames(rep(0.5, length(test_ids)), test_ids)
    if (length(selected) >= 2L) {
      test_dds <- DESeqDataSetFromMatrix(countData = round(all_counts[rownames(train_counts), test_ids]),
        colData = metadata[test_ids, , drop = FALSE], design = ~ group)
      test_dds <- estimateSizeFactors(test_dds, geoMeans = exp(rowMeans(log(train_counts))))
      dispersionFunction(test_dds) <- dispersionFunction(dds)
      test_vst <- varianceStabilizingTransformation(test_dds, blind = FALSE)
      x_train <- as.data.frame(t(vst_matrix[selected, train_ids, drop = FALSE]))
      x_test <- as.data.frame(t(assay(test_vst)[selected, test_ids, drop = FALSE]))
      colnames(x_train) <- colnames(x_test) <- make.unique(make.names(selected))
      scaling <- caret::preProcess(x_train, method = c("center", "scale"))
      x_train <- predict(scaling, x_train)
      x_test <- predict(scaling, x_test)
      y <- factor(ifelse(labels == "OP", "OVX", "Sham"), levels = c("OVX", "Sham"))
      set.seed(seed + 8L)
      fit <- pamr::pamr.train(list(x = t(x_train), y = y), threshold = pam_threshold)
      p <- pamr::pamr.predict(fit, t(x_test), threshold = pam_threshold, type = "posterior")
      probability[] <- p[, "OVX"]
    }
    stopifnot(all(is.finite(probability)), all(probability >= 0 & probability <= 1))
    summary <- data.frame(permutation = id, observed = id == observed_id,
      baseline = length(selected) < 2L, DE_candidates = n_candidates,
      selected_features = length(selected), holdout_AUC = holdout_auc(probability))
    predictions <- data.frame(permutation = id, Sample_ID = test_ids,
      Observed = holdout_labels, Probability = probability)
    saveRDS(list(summary = summary, training_labels = setNames(labels, train_ids),
      selection = selection_records, selection_runs = selection_runs,
      selected_features = selected, predictions = predictions),
      file.path(output_dir, "permutations", sprintf("permutation_%03d.rds", id)))
    list(summary = summary, predictions = predictions)
  }, warning = function(w) {
    expected <- c("one multinomial or binomial class has fewer than 8",
                  "Option grouped=FALSE enforced in cv.glmnet",
                  "There are no Tentative attributes! Returning original object.")
    if (!any(vapply(expected, grepl, logical(1), x = conditionMessage(w), fixed = TRUE))) stop(w)
    invokeRestart("muffleWarning")
  })
}

run_all <- function() {
  cluster <- parallel::makePSOCKcluster(workers, outfile = file.path(output_dir, "workers.log"))
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterExport(cluster, c("all_counts", "metadata", "train_ids", "test_ids", "type_map",
    "holdout_labels", "arrangements", "observed_id", "output_dir", "seed", "stability_iterations",
    "selection_cutoff", "fdr_cutoff", "lfc_cutoff", "lasso_cv_folds", "boruta_seed_offset",
    "boruta_pvalue", "boruta_max_runs", "boruta_num_trees", "boruta_threads", "pam_threshold",
    "make_foldid", "initialize_lasso_check", "fit_lasso", "holdout_auc", "run_permutation"))
  parallel::clusterEvalQ(cluster, {
    options(stringsAsFactors = FALSE, warn = 1)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    suppressPackageStartupMessages({library(DESeq2); library(glmnet); library(Boruta); library(caret); library(pamr)})
    initialize_lasso_check()
    NULL
  })
  parallel::parLapplyLB(cluster, seq_len(ncol(arrangements)), run_permutation, chunk.size = 1L)
}

runs <- run_all()
summary <- do.call(rbind, lapply(runs, `[[`, "summary"))
predictions <- do.call(rbind, lapply(runs, `[[`, "predictions"))
stopifnot(identical(summary$permutation, seq_len(924L)), nrow(predictions) == 924L * 8L,
          all(is.finite(summary$holdout_AUC)))
observed_auc <- summary$holdout_AUC[observed_id]
at_least_observed <- sum(summary$holdout_AUC >= observed_auc)
result <- data.frame(PAM_shrinkage_parameter = pam_threshold, Internal_holdout_AUC = observed_auc,
  Arrangements_at_least_observed = at_least_observed, Total_arrangements = nrow(summary),
  Unadjusted_P = at_least_observed / nrow(summary))
write.csv(summary, file.path(output_dir, "permutation_details.csv"), row.names = FALSE)
write.csv(predictions, file.path(output_dir, "permutation_predictions.csv"), row.names = FALSE)
write.csv(result, file.path(output_dir, "permutation_summary.csv"), row.names = FALSE)
print(result)
