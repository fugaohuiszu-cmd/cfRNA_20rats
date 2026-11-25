# 用途：训练集差异分析及稳定特征筛选
# 输入：DE_objects.rds中的原始计数
# 输出：样本划分、特征表及strict_feature_selection.rds

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(glmnet)
  library(Boruta)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
input_rds <- file.path(project_dir, "03_results", "02_DE", "rds", "DE_objects.rds")
result_dir <- file.path(project_dir, "03_results", "03_ML", "02_strict_model", "01_feature_selection")
table_dir <- file.path(result_dir, "tables")
rds_dir <- file.path(result_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260717L
boruta_seed_offset <- 100000L
train_fraction <- 0.60
stability_iterations <- 100L
stability_sample_mode <- "full_training"
lasso_stable_cutoff <- 0.10
boruta_stable_cutoff <- 0.10
fdr_cutoff <- 0.05
lfc_cutoff <- 0.50
lasso_cv_folds <- 4L
boruta_pvalue <- 0.001
boruta_max_runs <- 300L
boruta_num_trees <- 300L
boruta_threads <- 4L

obj <- readRDS(input_rds)
required_objects <- c("dds", "result", "metadata")
missing_objects <- setdiff(required_objects, names(obj))
if (length(missing_objects) > 0L) {
  stop("DE object is missing: ", paste(missing_objects, collapse = ", "))
}

metadata <- obj$metadata
raw_counts <- counts(obj$dds, normalized = FALSE)
if (!identical(rownames(metadata), colnames(raw_counts))) {
  stop("Metadata rows do not match count-matrix columns.")
}
metadata$group <- factor(metadata$group, levels = c("N", "OP"))
if (!identical(as.integer(table(metadata$group)), c(10L, 10L))) {
  stop("Expected 10 N and 10 OP samples.")
}

set.seed(seed)
train_ids_unordered <- unlist(lapply(levels(metadata$group), function(group_name) {
  group_ids <- rownames(metadata)[metadata$group == group_name]
  sample(group_ids, size = floor(length(group_ids) * train_fraction), replace = FALSE)
}), use.names = FALSE)
train_ids <- rownames(metadata)[rownames(metadata) %in% train_ids_unordered]
test_ids <- setdiff(rownames(metadata), train_ids)

split_table <- data.frame(
  sample_id = rownames(metadata),
  group = as.character(metadata$group),
  dataset = ifelse(rownames(metadata) %in% train_ids, "Training", "Validation")
)
if (!identical(as.integer(table(split_table$group, split_table$dataset)[, "Training"]), c(6L, 6L))) {
  stop("Training split must contain 6 N and 6 OP samples.")
}
if (!identical(as.integer(table(split_table$group, split_table$dataset)[, "Validation"]), c(4L, 4L))) {
  stop("Validation split must contain 4 N and 4 OP samples.")
}

metadata_train <- metadata[train_ids, , drop = FALSE]
metadata_train$group <- droplevels(metadata_train$group)
counts_train <- raw_counts[, train_ids, drop = FALSE]
counts_train <- counts_train[rowSums(counts_train >= 1) >= 1L, , drop = FALSE]

dds_train <- DESeqDataSetFromMatrix(
  countData = counts_train,
  colData = metadata_train,
  design = ~ group
)
dds_train <- DESeq(dds_train, fitType = "local", quiet = TRUE)
de_train <- as.data.frame(results(
  dds_train,
  contrast = c("group", "OP", "N"),
  alpha = fdr_cutoff
))
de_train$feature_id <- rownames(de_train)

type_map <- setNames(obj$result$RNA_type, obj$result$feature_id)
de_train$RNA_type <- unname(type_map[de_train$feature_id])
if (anyNA(de_train$RNA_type)) {
  stop("RNA type is missing for one or more training DE features.")
}

de_train <- de_train[, c(
  "feature_id", "RNA_type", "baseMean", "log2FoldChange",
  "lfcSE", "stat", "pvalue", "padj"
)]
candidate_table <- de_train %>%
  filter(
    !is.na(.data$padj),
    .data$padj < fdr_cutoff,
    is.finite(.data$log2FoldChange),
    abs(.data$log2FoldChange) > lfc_cutoff
  ) %>%
  arrange(.data$padj, desc(abs(.data$log2FoldChange)))

candidate_features <- unique(candidate_table$feature_id)
if (length(candidate_features) < 2L) {
  stop("Fewer than two training-set DE candidates were found.")
}

vst_train <- varianceStabilizingTransformation(dds_train, blind = TRUE, fitType = "local")
vst_matrix_train <- assay(vst_train)
x_all <- t(vst_matrix_train[candidate_features, train_ids, drop = FALSE])
variable_features <- apply(x_all, 2, sd, na.rm = TRUE) > 1e-8
x_all <- x_all[, variable_features, drop = FALSE]
candidate_features <- colnames(x_all)
candidate_table <- candidate_table %>% filter(.data$feature_id %in% candidate_features)
if (ncol(x_all) < 2L) {
  stop("Fewer than two variable training candidates remain after VST.")
}

group_train <- metadata_train$group
names(group_train) <- rownames(metadata_train)

make_subsample <- function(group, fraction) {
  ids_by_group <- split(names(group), group)
  sampled <- lapply(ids_by_group, function(ids) {
    n_take <- round(length(ids) * fraction)
    if (n_take < 2L || n_take >= length(ids)) {
      stop("Invalid class-specific subsample size: ", n_take)
    }
    sample(ids, n_take, replace = FALSE)
  })
  unname(unlist(sampled, use.names = FALSE))
}

make_foldid <- function(group, nfolds) {
  if (min(table(group)) < nfolds) {
    stop("Too few samples per class for LASSO cross-validation.")
  }
  foldid <- integer(length(group))
  for (group_name in levels(group)) {
    idx <- which(group == group_name)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  foldid
}

extract_lasso <- function(x, y, foldid) {
  fit <- suppressWarnings(cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    type.measure = "deviance",
    nfolds = length(unique(foldid)),
    foldid = foldid,
    standardize = TRUE
  ))
  coefficient_matrix <- as.matrix(coef(fit, s = "lambda.1se"))
  coefficient_vector <- coefficient_matrix[rownames(coefficient_matrix) != "(Intercept)", 1]
  names(coefficient_vector) <- rownames(coefficient_matrix)[rownames(coefficient_matrix) != "(Intercept)"]
  list(
    selected = names(coefficient_vector)[abs(coefficient_vector) > 0],
    coefficients = coefficient_vector
  )
}

feature_names <- colnames(x_all)
lasso_hits <- setNames(integer(length(feature_names)), feature_names)
boruta_hits <- setNames(integer(length(feature_names)), feature_names)
lasso_abs_coef_sum <- setNames(numeric(length(feature_names)), feature_names)

iteration_summary <- data.frame(
  iteration = seq_len(stability_iterations),
  n_samples = integer(stability_iterations),
  n_lasso_selected = integer(stability_iterations),
  n_boruta_confirmed = integer(stability_iterations),
  boruta_elapsed_seconds = numeric(stability_iterations)
)

for (iteration in seq_len(stability_iterations)) {
  set.seed(seed + iteration)
  sample_ids <- names(group_train)
  x_sub <- x_all[sample_ids, , drop = FALSE]
  group_sub <- droplevels(group_train[sample_ids])
  y_sub <- ifelse(group_sub == "OP", 1L, 0L)

  foldid <- make_foldid(group_sub, lasso_cv_folds)
  lasso_result <- extract_lasso(x_sub, y_sub, foldid)
  lasso_selected <- lasso_result$selected
  lasso_hits[lasso_selected] <- lasso_hits[lasso_selected] + 1L
  lasso_abs_coef_sum[names(lasso_result$coefficients)] <-
    lasso_abs_coef_sum[names(lasso_result$coefficients)] + abs(lasso_result$coefficients)

  set.seed(seed + boruta_seed_offset + iteration)
  start_time <- proc.time()[[3]]
  getimp_warning_count <- 0L
  boruta_fit <- withCallingHandlers(
    Boruta(
      x = as.data.frame(x_sub, check.names = FALSE),
      y = group_sub,
      getImp = getImpRfGini,
      pValue = boruta_pvalue,
      maxRuns = boruta_max_runs,
      doTrace = 0,
      num.trees = boruta_num_trees,
      num.threads = boruta_threads
    ),
    warning = function(w) {
      if (grepl("getImp result contains", conditionMessage(w), fixed = TRUE)) {
        getimp_warning_count <<- getimp_warning_count + 1L
      }
    }
  )
  if (getimp_warning_count > 0L) {
    stop("Boruta getImp produced non-finite values in iteration ", iteration, ".")
  }
  boruta_fit <- TentativeRoughFix(boruta_fit)
  boruta_stats <- attStats(boruta_fit)
  boruta_confirmed <- rownames(boruta_stats)[boruta_stats$decision == "Confirmed"]
  boruta_hits[boruta_confirmed] <- boruta_hits[boruta_confirmed] + 1L

  iteration_summary$n_samples[iteration] <- length(sample_ids)
  iteration_summary$n_lasso_selected[iteration] <- length(lasso_selected)
  iteration_summary$n_boruta_confirmed[iteration] <- length(boruta_confirmed)
  iteration_summary$boruta_elapsed_seconds[iteration] <- proc.time()[[3]] - start_time

  if (iteration %% 10L == 0L) {
    message("Completed stability iteration ", iteration, "/", stability_iterations)
  }
}

stability_table <- candidate_table %>%
  select(all_of(c(
    "feature_id", "RNA_type", "baseMean", "log2FoldChange",
    "lfcSE", "stat", "pvalue", "padj"
  ))) %>%
  mutate(
    lasso_selection_count = as.integer(lasso_hits[.data$feature_id]),
    lasso_selection_frequency = as.numeric(lasso_hits[.data$feature_id]) / stability_iterations,
    lasso_mean_absolute_coefficient = as.numeric(lasso_abs_coef_sum[.data$feature_id]) / stability_iterations,
    boruta_confirmed_count = as.integer(boruta_hits[.data$feature_id]),
    boruta_confirmed_frequency = as.numeric(boruta_hits[.data$feature_id]) / stability_iterations,
    stable_lasso = .data$lasso_selection_frequency >= lasso_stable_cutoff,
    stable_boruta = .data$boruta_confirmed_frequency >= boruta_stable_cutoff,
    final_feature = .data$stable_lasso & .data$stable_boruta
  ) %>%
  arrange(
    desc(.data$final_feature),
    desc(.data$lasso_selection_frequency + .data$boruta_confirmed_frequency),
    .data$padj
  )

lasso_stable <- stability_table %>% filter(.data$stable_lasso)
boruta_stable <- stability_table %>% filter(.data$stable_boruta)
final_features <- stability_table %>% filter(.data$final_feature)

write.csv(split_table, file.path(table_dir, "sample_split.csv"), row.names = FALSE, quote = FALSE)
write.csv(de_train, file.path(table_dir, "training_DE_all.csv"), row.names = FALSE, quote = FALSE)
write.csv(candidate_table, file.path(table_dir, "training_DE_candidates.csv"), row.names = FALSE, quote = FALSE)
write.csv(stability_table, file.path(table_dir, "feature_stability.csv"), row.names = FALSE, quote = FALSE)
write.csv(lasso_stable, file.path(table_dir, "LASSO_stable_features.csv"), row.names = FALSE, quote = FALSE)
write.csv(boruta_stable, file.path(table_dir, "Boruta_stable_features.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_features, file.path(table_dir, "final_features.csv"), row.names = FALSE, quote = FALSE)
write.csv(iteration_summary, file.path(table_dir, "stability_iterations.csv"), row.names = FALSE, quote = FALSE)

parameters <- data.frame(
  parameter = c(
    "seed", "train_fraction", "stability_iterations", "stability_sample_mode",
    "DE_padj_cutoff", "DE_abs_log2FC_cutoff", "LASSO_frequency_cutoff",
    "Boruta_frequency_cutoff", "LASSO_folds", "Boruta_pvalue",
    "Boruta_max_runs", "Boruta_num_trees"
  ),
  value = c(
    seed, train_fraction, stability_iterations, stability_sample_mode,
    fdr_cutoff, lfc_cutoff, lasso_stable_cutoff,
    boruta_stable_cutoff, lasso_cv_folds, boruta_pvalue,
    boruta_max_runs, boruta_num_trees
  )
)
write.csv(parameters, file.path(table_dir, "parameters.csv"), row.names = FALSE, quote = FALSE)

saveRDS(list(
  parameters = parameters,
  train_ids = train_ids,
  test_ids = test_ids,
  metadata_train = metadata_train,
  dds_train = dds_train,
  vst_matrix_train = vst_matrix_train,
  training_de = de_train,
  candidates = candidate_table,
  stability = stability_table,
  final_features = final_features
), file.path(rds_dir, "strict_feature_selection.rds"))

cat("Training samples: ", length(train_ids), "\n", sep = "")
cat("Validation samples: ", length(test_ids), "\n", sep = "")
cat("Training DE candidates: ", nrow(candidate_table), "\n", sep = "")
cat("Stable LASSO: ", nrow(lasso_stable), "\n", sep = "")
cat("Stable Boruta: ", nrow(boruta_stable), "\n", sep = "")
cat("Final intersection: ", nrow(final_features), "\n", sep = "")

if (nrow(final_features) < 2L) {
  stop("Fewer than two strict-training final features were retained.")
}
