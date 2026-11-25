# 用途：固定特征的内部交叉验证及模型训练
# 输入：strict_feature_selection.rds
# 输出：交叉验证指标及Internal_CV_models.rds

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(foreach)
  library(ggplot2)
  library(pROC)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
input_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "01_feature_selection", "rds", "strict_feature_selection.rds"
)
result_dir <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model", "02_internal_CV"
)
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")
rds_dir <- file.path(result_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260717L
cv_folds <- 4L
cv_repeats <- 50L
expected_feature_count <- 5L
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- 2L
n_cores <- max(1L, min(8L, detected_cores - 2L))

model_methods <- c(
  GLMNETLasso = "glmnet",
  GLMNETRidge = "glmnet",
  SVMLIN = "svmLinear2",
  SVMRAD = "svmRadial",
  RF = "ranger",
  GBM = "gbm",
  LDA = "lda",
  PAM = "pam",
  XGB = "native_xgboost",
  LGBM = "native_lightgbm"
)

required_packages <- c(
  "caret", "pROC", "ggplot2", "dplyr", "foreach", "doParallel",
  "glmnet", "e1071", "kernlab", "ranger", "gbm", "MASS", "pamr",
  "xgboost", "lightgbm", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

metric_values <- function(obs, prob, threshold = 0.5) {
  obs_chr <- as.character(obs)
  prob <- as.numeric(prob)
  if (length(obs_chr) != length(prob) || any(!is.finite(prob))) {
    stop("Invalid observations or prediction probabilities.")
  }
  if (!all(obs_chr %in% c("N_like", "OP_like"))) {
    stop("Unexpected class labels in predictions.")
  }
  if (length(unique(obs_chr)) != 2L) {
    stop("Each assessment fold must contain both classes.")
  }

  roc_obj <- pROC::roc(
    response = factor(obs_chr, levels = c("N_like", "OP_like")),
    predictor = prob,
    levels = c("N_like", "OP_like"),
    direction = "<",
    quiet = TRUE
  )
  pred_op <- prob >= threshold
  true_op <- obs_chr == "OP_like"
  tp <- sum(pred_op & true_op)
  fp <- sum(pred_op & !true_op)
  fn <- sum(!pred_op & true_op)
  tn <- sum(!pred_op & !true_op)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  accuracy <- (tp + tn) / length(obs_chr)
  precision <- if ((tp + fp) == 0L) 0 else tp / (tp + fp)
  f1 <- if ((precision + sensitivity) == 0) {
    0
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }

  c(
    AUC = as.numeric(pROC::auc(roc_obj)),
    Spec = specificity,
    Sens = sensitivity,
    F1 = f1,
    Accuracy = accuracy
  )
}

choose_threshold <- function(obs, prob) {
  roc_obj <- pROC::roc(
    response = factor(as.character(obs), levels = c("N_like", "OP_like")),
    predictor = as.numeric(prob),
    levels = c("N_like", "OP_like"),
    direction = "<",
    quiet = TRUE
  )
  threshold <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )
  threshold <- as.numeric(unlist(threshold, use.names = FALSE))[1]
  if (!is.finite(threshold)) stop("Failed to determine a CV threshold.")
  pmin(pmax(threshold, 0), 1)
}

scale_pair <- function(x_train, x_test) {
  center <- colMeans(x_train)
  scale_value <- apply(x_train, 2L, stats::sd)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2L, center, "-"), 2L, scale_value, "/"),
    test = sweep(sweep(x_test, 2L, center, "-"), 2L, scale_value, "/"),
    center = center,
    scale = scale_value
  )
}

fit_native_model <- function(model_name, x, y, parameter_row, model_seed) {
  y_binary <- ifelse(as.character(y) == "OP_like", 1, 0)

  if (model_name == "XGB") {
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x), label = y_binary)
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = as.integer(parameter_row$max_depth),
      eta = as.numeric(parameter_row$eta),
      min_child_weight = as.numeric(parameter_row$min_child_weight),
      subsample = as.numeric(parameter_row$subsample),
      colsample_bytree = as.numeric(parameter_row$colsample_bytree),
      gamma = as.numeric(parameter_row$gamma),
      nthread = 1L,
      seed = as.integer(model_seed)
    )
    return(xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = as.integer(parameter_row$nrounds),
      verbose = 0
    ))
  }

  if (model_name == "LGBM") {
    dtrain <- lightgbm::lgb.Dataset(data = as.matrix(x), label = y_binary)
    params <- list(
      objective = "binary",
      metric = "binary_logloss",
      boosting = "gbdt",
      max_depth = as.integer(parameter_row$max_depth),
      num_leaves = as.integer(parameter_row$num_leaves),
      learning_rate = as.numeric(parameter_row$learning_rate),
      feature_fraction = as.numeric(parameter_row$feature_fraction),
      min_data_in_leaf = as.integer(parameter_row$min_data_in_leaf),
      lambda_l2 = as.numeric(parameter_row$lambda_l2),
      num_threads = 1L,
      seed = as.integer(model_seed),
      deterministic = TRUE,
      force_col_wise = TRUE,
      verbosity = -1
    )
    return(lightgbm::lgb.train(
      params = params,
      data = dtrain,
      nrounds = as.integer(parameter_row$nrounds),
      verbose = -1
    ))
  }

  if (model_name == "CATBOOST") {
    pool <- catboost::catboost.load_pool(data = as.matrix(x), label = y_binary)
    params <- list(
      loss_function = "Logloss",
      eval_metric = "AUC",
      iterations = as.integer(parameter_row$iterations),
      depth = as.integer(parameter_row$depth),
      learning_rate = as.numeric(parameter_row$learning_rate),
      l2_leaf_reg = as.numeric(parameter_row$l2_leaf_reg),
      random_strength = as.numeric(parameter_row$random_strength),
      random_seed = as.integer(model_seed),
      thread_count = 1L,
      verbose = FALSE,
      allow_writing_files = FALSE
    )
    return(catboost::catboost.train(pool, NULL, params = params))
  }

  stop("Unknown native model: ", model_name)
}

predict_native_model <- function(model_name, model, x) {
  if (model_name == "XGB") {
    return(as.numeric(predict(
      model,
      xgboost::xgb.DMatrix(data = as.matrix(x))
    )))
  }
  if (model_name == "LGBM") {
    return(as.numeric(predict(model, as.matrix(x))))
  }
  if (model_name == "CATBOOST") {
    pool <- catboost::catboost.load_pool(data = as.matrix(x))
    return(as.numeric(catboost::catboost.predict(
      model,
      pool,
      prediction_type = "Probability"
    )))
  }
  stop("Unknown native model: ", model_name)
}

native_grid <- function(model_name) {
  if (model_name == "XGB") {
    return(expand.grid(
      nrounds = c(50L, 100L, 200L),
      max_depth = c(1L, 2L),
      eta = c(0.01, 0.05, 0.10),
      min_child_weight = c(1, 3),
      subsample = 0.8,
      colsample_bytree = 0.8,
      gamma = 0,
      KEEP.OUT.ATTRS = FALSE
    ))
  }
  if (model_name == "LGBM") {
    grid <- expand.grid(
      nrounds = c(50L, 100L, 200L),
      max_depth = c(1L, 2L),
      learning_rate = c(0.01, 0.05, 0.10),
      num_leaves = c(2L, 3L),
      feature_fraction = 0.8,
      min_data_in_leaf = c(1L, 3L),
      lambda_l2 = c(0, 1),
      KEEP.OUT.ATTRS = FALSE
    )
    return(grid[grid$num_leaves <= 2^grid$max_depth, , drop = FALSE])
  }
  if (model_name == "CATBOOST") {
    return(expand.grid(
      iterations = c(50L, 100L, 200L),
      depth = c(1L, 2L),
      learning_rate = c(0.01, 0.05, 0.10),
      l2_leaf_reg = c(1, 3),
      random_strength = c(0, 1),
      KEEP.OUT.ATTRS = FALSE
    ))
  }
  stop("Unknown native model: ", model_name)
}

tune_native_model <- function(model_name, x, y, fold_index, model_position) {
  grid <- native_grid(model_name)
  package_name <- switch(
    model_name,
    XGB = "xgboost",
    LGBM = "lightgbm",
    CATBOOST = "catboost"
  )

  tuning <- foreach::foreach(
    grid_row = seq_len(nrow(grid)),
    .combine = rbind,
    .inorder = TRUE,
    .packages = c(package_name, "pROC"),
    .export = c(
      "scale_pair", "fit_native_model", "predict_native_model",
      "metric_values", "seed"
    )
  ) %dopar% {
    parameter_row <- grid[grid_row, , drop = FALSE]
    fold_scores <- lapply(seq_along(fold_index), function(fold_number) {
      analysis_rows <- fold_index[[fold_number]]
      assessment_rows <- setdiff(seq_len(nrow(x)), analysis_rows)
      scaled <- scale_pair(
        x[analysis_rows, , drop = FALSE],
        x[assessment_rows, , drop = FALSE]
      )
      fit <- fit_native_model(
        model_name,
        scaled$train,
        y[analysis_rows],
        parameter_row,
        seed + model_position * 1000000L + grid_row * 1000L + fold_number
      )
      prob <- predict_native_model(model_name, fit, scaled$test)
      metric_values(y[assessment_rows], prob, threshold = 0.5)
    })
    fold_scores <- do.call(rbind, fold_scores)
    cbind(
      parameter_row,
      Mean_AUC = mean(fold_scores[, "AUC"]),
      Mean_Spec = mean(fold_scores[, "Spec"]),
      Mean_Sens = mean(fold_scores[, "Sens"])
    )
  }

  if (nrow(tuning) != nrow(grid) || any(!is.finite(tuning$Mean_AUC))) {
    stop(model_name, " native tuning did not return complete finite results.")
  }
  best_index <- order(
    -tuning$Mean_AUC,
    -tuning$Mean_Sens,
    -tuning$Mean_Spec
  )[1]
  best_tune <- tuning[best_index, colnames(grid), drop = FALSE]

  cv_predictions <- foreach::foreach(
    fold_number = seq_along(fold_index),
    .combine = rbind,
    .inorder = TRUE,
    .packages = package_name,
    .export = c(
      "scale_pair", "fit_native_model", "predict_native_model", "seed"
    )
  ) %dopar% {
    analysis_rows <- fold_index[[fold_number]]
    assessment_rows <- setdiff(seq_len(nrow(x)), analysis_rows)
    scaled <- scale_pair(
      x[analysis_rows, , drop = FALSE],
      x[assessment_rows, , drop = FALSE]
    )
    fit <- fit_native_model(
      model_name,
      scaled$train,
      y[analysis_rows],
      best_tune,
      seed + model_position * 1000000L + 500000L + fold_number
    )
    prob <- predict_native_model(model_name, fit, scaled$test)
    data.frame(
      rowIndex = assessment_rows,
      obs = as.character(y[assessment_rows]),
      prob = prob,
      Resample = names(fold_index)[fold_number],
      stringsAsFactors = FALSE
    )
  }

  full_scaled <- scale_pair(x, x)
  final_model <- fit_native_model(
    model_name,
    full_scaled$train,
    y,
    best_tune,
    seed + model_position * 1000000L + 900000L
  )

  list(
    best_tune = best_tune,
    tuning = tuning,
    cv_predictions = cv_predictions,
    final_model = final_model,
    preprocessing = list(center = full_scaled$center, scale = full_scaled$scale)
  )
}

caret_grid <- function(model_name, x) {
  if (model_name %in% c("GLMNETElastic", "GLMNETLasso", "GLMNETRidge")) {
    alpha_value <- switch(
      model_name,
      GLMNETElastic = 0.5,
      GLMNETLasso = 1,
      GLMNETRidge = 0
    )
    return(expand.grid(
      alpha = alpha_value,
      lambda = 10^seq(-4, 1, length.out = 30)
    ))
  }
  if (model_name == "SVMLIN") {
    return(expand.grid(cost = 10^seq(-2, 2, length.out = 9)))
  }
  if (model_name == "SVMRAD") {
    sigma_estimate <- kernlab::sigest(scale(x), frac = 1)
    sigma_center <- stats::median(as.numeric(sigma_estimate))
    if (!is.finite(sigma_center) || sigma_center <= 0) {
      stop("Unable to estimate a valid radial-SVM sigma value.")
    }
    return(expand.grid(
      sigma = sigma_center * c(0.1, 0.3, 1, 3, 10),
      C = 10^seq(-2, 2, length.out = 9)
    ))
  }
  if (model_name == "RF") {
    return(expand.grid(
      mtry = seq_len(ncol(x)),
      splitrule = "gini",
      min.node.size = c(1L, 2L, 3L)
    ))
  }
  if (model_name == "GBM") {
    return(expand.grid(
      n.trees = c(50L, 100L, 200L),
      interaction.depth = c(1L, 2L),
      shrinkage = c(0.01, 0.05, 0.10),
      n.minobsinnode = c(1L, 2L)
    ))
  }
  if (model_name == "PAM") {
    return(expand.grid(threshold = seq(0, 5, length.out = 20)))
  }
  NULL
}

train_caret_model <- function(model_name, x, y, train_control, model_position) {
  method <- model_methods[[model_name]]
  grid <- caret_grid(model_name, x)
  args <- list(
    x = as.data.frame(x, check.names = FALSE),
    y = y,
    method = method,
    metric = "ROC",
    trControl = train_control,
    preProcess = c("center", "scale")
  )
  if (model_name %in% c("LR", "LDA")) {
    args$tuneLength <- 1L
  } else {
    args$tuneGrid <- grid
  }
  if (model_name == "RF") {
    args$importance <- "impurity"
    args$num.trees <- 500L
    args$num.threads <- 1L
  }
  if (model_name == "GBM") {
    args$distribution <- "bernoulli"
    args$bag.fraction <- 0.7
    args$verbose <- FALSE
  }

  set.seed(seed + model_position)
  fit <- do.call(caret::train, args)
  if (is.null(fit$pred) || nrow(fit$pred) == 0L) {
    stop(model_name, " did not return saved CV predictions.")
  }
  if (!"OP_like" %in% colnames(fit$pred)) {
    stop(model_name, " predictions do not contain OP_like probabilities.")
  }
  cv_predictions <- data.frame(
    rowIndex = fit$pred$rowIndex,
    obs = as.character(fit$pred$obs),
    prob = as.numeric(fit$pred$OP_like),
    Resample = fit$pred$Resample,
    stringsAsFactors = FALSE
  )
  list(
    best_tune = fit$bestTune,
    tuning = fit$results,
    cv_predictions = cv_predictions,
    final_model = fit
  )
}

collect_resample_metrics <- function(model_name, predictions, threshold) {
  pieces <- split(predictions, predictions$Resample)
  rows <- lapply(names(pieces), function(resample_name) {
    values <- metric_values(
      pieces[[resample_name]]$obs,
      pieces[[resample_name]]$prob,
      threshold
    )
    data.frame(
      Model = model_name,
      Resample = resample_name,
      Threshold = threshold,
      AUC = values[["AUC"]],
      Spec = values[["Spec"]],
      Sens = values[["Sens"]],
      F1 = values[["F1"]],
      Accuracy = values[["Accuracy"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

metrics_to_long <- function(metrics) {
  metric_names <- c("AUC", "Spec", "Sens", "F1", "Accuracy")
  do.call(rbind, lapply(metric_names, function(metric_name) {
    data.frame(
      Model = metrics$Model,
      Resample = metrics$Resample,
      Threshold = metrics$Threshold,
      Metric = metric_name,
      Value = metrics[[metric_name]],
      stringsAsFactors = FALSE
    )
  }))
}

summarize_metrics <- function(metrics_long) {
  pieces <- split(
    metrics_long,
    paste(metrics_long$Model, metrics_long$Metric, sep = "\r")
  )
  do.call(rbind, lapply(pieces, function(x) {
    values <- x$Value
    if (any(!is.finite(values))) {
      stop("Non-finite internal-CV metric values were found.")
    }
    n <- length(values)
    mean_value <- mean(values)
    sd_value <- stats::sd(values)
    se_value <- sd_value / sqrt(n)
    data.frame(
      Model = x$Model[1],
      Metric = x$Metric[1],
      N_Resamples = n,
      Mean = mean_value,
      SD = sd_value,
      CI_Low = max(0, mean_value - 1.96 * se_value),
      CI_High = min(1, mean_value + 1.96 * se_value),
      stringsAsFactors = FALSE
    )
  }))
}

make_plot <- function(metrics_long, summary_table) {
  metric_order <- c("AUC", "Spec", "Sens", "F1", "Accuracy")
  auc_order <- summary_table %>%
    filter(.data$Metric == "AUC") %>%
    arrange(desc(.data$Mean)) %>%
    pull(.data$Model)
  y_levels <- rev(auc_order)
  base_palette <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#0072B2", "#D55E00", "#CC79A7"
  )
  palette <- setNames(base_palette[seq_along(model_methods)], names(model_methods))

  metrics_long$Metric <- factor(metrics_long$Metric, levels = metric_order)
  metrics_long$Model <- factor(metrics_long$Model, levels = y_levels)
  summary_table$Metric <- factor(summary_table$Metric, levels = metric_order)
  summary_table$Model <- factor(summary_table$Model, levels = y_levels)
  summary_table$Label <- sprintf(
    "%.3f \u00b1 %.3f\n95%% CI %.3f-%.3f",
    summary_table$Mean,
    summary_table$SD,
    summary_table$CI_Low,
    summary_table$CI_High
  )

  ggplot(metrics_long, aes(x = .data$Value, y = .data$Model)) +
    geom_jitter(
      aes(colour = .data$Model),
      width = 0,
      height = 0.14,
      alpha = 0.16,
      size = 0.48
    ) +
    geom_boxplot(
      aes(fill = .data$Model, colour = .data$Model),
      width = 0.42,
      outlier.shape = NA,
      alpha = 0.26,
      linewidth = 0.45
    ) +
    geom_text(
      data = summary_table,
      aes(x = 0.018, y = .data$Model, label = .data$Label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = 1.95,
      lineheight = 0.82,
      colour = "black"
    ) +
    facet_grid(cols = vars(.data$Metric)) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(
      title = paste0(
        "The internal validation results of ",
        length(model_methods),
        " machine learning (ML) models."
      ),
      x = NULL,
      y = NULL
    ) +
    guides(colour = "none", fill = "none") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(colour = "black", size = 8.5),
      axis.text.y = element_text(colour = "black", size = 8.8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "#8A8A8A", fill = NA, linewidth = 0.45),
      strip.background = element_rect(fill = "#F1F1F1", colour = "#8A8A8A", linewidth = 0.45),
      strip.text = element_text(face = "bold", colour = "black", size = 10),
      plot.margin = grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")
    )
}

main <- function() {
  object <- readRDS(input_rds)
  required_objects <- c(
    "train_ids", "test_ids", "metadata_train", "vst_matrix_train",
    "final_features"
  )
  missing_objects <- setdiff(required_objects, names(object))
  if (length(missing_objects) > 0L) {
    stop("Strict feature object is missing: ", paste(missing_objects, collapse = ", "))
  }

  feature_ids <- object$final_features$feature_id
  if (length(feature_ids) != expected_feature_count) {
    stop(
      "Expected ", expected_feature_count,
      " final features, but found ", length(feature_ids), "."
    )
  }
  if (!all(feature_ids %in% rownames(object$vst_matrix_train))) {
    stop("Final features are missing from the training VST matrix.")
  }

  metadata_train <- object$metadata_train
  x_train <- t(object$vst_matrix_train[feature_ids, , drop = FALSE])
  if (!identical(rownames(x_train), rownames(metadata_train))) {
    stop("Training expression rows do not match training metadata rows.")
  }
  feature_map <- data.frame(
    feature_id = feature_ids,
    model_feature = make.unique(make.names(feature_ids)),
    stringsAsFactors = FALSE
  )
  colnames(x_train) <- feature_map$model_feature
  y_train <- factor(
    ifelse(metadata_train$group == "OP", "OP_like", "N_like"),
    levels = c("OP_like", "N_like")
  )
  names(y_train) <- rownames(x_train)
  if (!identical(as.integer(table(y_train)), c(6L, 6L))) {
    stop("Training labels must contain 6 OP_like and 6 N_like samples.")
  }

  set.seed(seed)
  fold_index <- caret::createMultiFolds(
    y_train,
    k = cv_folds,
    times = cv_repeats
  )
  if (length(fold_index) != cv_folds * cv_repeats) {
    stop("Unexpected number of repeated-CV folds.")
  }
  fold_index_out <- lapply(
    fold_index,
    function(index) setdiff(seq_along(y_train), index)
  )
  if (any(vapply(
    fold_index_out,
    function(index) length(unique(y_train[index])) != 2L,
    logical(1)
  ))) {
    stop("At least one assessment fold does not contain both classes.")
  }

  set.seed(seed + 1L)
  caret_seeds <- vector("list", length(fold_index) + 1L)
  for (i in seq_along(fold_index)) {
    caret_seeds[[i]] <- sample.int(10000000L, 200L)
  }
  caret_seeds[[length(caret_seeds)]] <- sample.int(10000000L, 1L)
  train_control <- caret::trainControl(
    method = "repeatedcv",
    number = cv_folds,
    repeats = cv_repeats,
    index = fold_index,
    indexOut = fold_index_out,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = TRUE,
    seeds = caret_seeds
  )

  cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cluster)
  on.exit({
    parallel::stopCluster(cluster)
    foreach::registerDoSEQ()
  }, add = TRUE)
  cat("Parallel workers: ", n_cores, "\n", sep = "")

  model_results <- vector("list", length(model_methods))
  names(model_results) <- names(model_methods)
  all_predictions <- vector("list", length(model_methods))
  all_metrics <- vector("list", length(model_methods))
  best_parameter_rows <- list()

  for (model_position in seq_along(model_methods)) {
    model_name <- names(model_methods)[model_position]
    cat(
      "Training ", model_name, " (", model_position, "/",
      length(model_methods), ")\n",
      sep = ""
    )

    if (model_name %in% c("XGB", "LGBM", "CATBOOST")) {
      result <- tune_native_model(
        model_name,
        x_train,
        y_train,
        fold_index,
        model_position
      )
    } else {
      result <- train_caret_model(
        model_name,
        x_train,
        y_train,
        train_control,
        model_position
      )
    }

    predictions <- result$cv_predictions
    expected_predictions <- length(y_train) * cv_repeats
    if (nrow(predictions) != expected_predictions) {
      stop(
        model_name, " returned ", nrow(predictions),
        " CV predictions; expected ", expected_predictions, "."
      )
    }
    if (length(unique(predictions$Resample)) != cv_folds * cv_repeats) {
      stop(model_name, " returned an unexpected number of resamples.")
    }
    threshold <- choose_threshold(predictions$obs, predictions$prob)
    predictions$Model <- model_name
    predictions$sample_id <- rownames(x_train)[predictions$rowIndex]
    predictions$Threshold <- threshold

    model_results[[model_name]] <- c(result, list(threshold = threshold))
    all_predictions[[model_name]] <- predictions
    all_metrics[[model_name]] <- collect_resample_metrics(
      model_name,
      predictions,
      threshold
    )

    best_tune <- result$best_tune
    best_parameter_rows[[model_name]] <- data.frame(
      Model = model_name,
      Parameter = colnames(best_tune),
      Value = vapply(best_tune, function(x) as.character(x[1]), character(1)),
      stringsAsFactors = FALSE
    )
  }

  predictions <- bind_rows(all_predictions)
  resample_metrics <- bind_rows(all_metrics)
  metrics_long <- metrics_to_long(resample_metrics)
  summary_table <- summarize_metrics(metrics_long)

  if (length(unique(resample_metrics$Model)) != length(model_methods)) {
    stop("Not all requested models returned complete results.")
  }
  expected_metric_rows <- length(model_methods) * cv_folds * cv_repeats
  if (nrow(resample_metrics) != expected_metric_rows) {
    stop("Unexpected number of resample-metric rows.")
  }

  fold_table <- bind_rows(lapply(names(fold_index_out), function(resample_name) {
    assessment_rows <- fold_index_out[[resample_name]]
    data.frame(
      Resample = resample_name,
      sample_id = rownames(x_train)[assessment_rows],
      group = as.character(y_train[assessment_rows]),
      stringsAsFactors = FALSE
    )
  }))
  best_parameters <- bind_rows(best_parameter_rows)

  write.csv(
    feature_map,
    file.path(table_dir, "feature_mapping.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    fold_table,
    file.path(table_dir, "CV_fold_assignments.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    predictions,
    file.path(table_dir, "Internal_CV_predictions.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    resample_metrics,
    file.path(table_dir, "Internal_CV_resample_metrics.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    summary_table,
    file.path(table_dir, "Internal_CV_summary.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    best_parameters,
    file.path(table_dir, "best_parameters.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(list(
    parameters = list(
      seed = seed,
      cv_folds = cv_folds,
      cv_repeats = cv_repeats,
      n_cores = n_cores,
      train_ids = object$train_ids,
      validation_ids = object$test_ids
    ),
    feature_map = feature_map,
    fold_index = fold_index,
    model_results = model_results,
    predictions = predictions,
    resample_metrics = resample_metrics,
    summary = summary_table
  ), file.path(rds_dir, "Internal_CV_models.rds"))

  plot_object <- make_plot(metrics_long, summary_table)
  pdf_file <- file.path(figure_dir, "Internal_CV_metrics.pdf")
  png_file <- file.path(figure_dir, "Internal_CV_metrics.png")
  ggsave(
    pdf_file,
    plot_object,
    width = 17,
    height = 7.2,
    device = grDevices::cairo_pdf
  )
  ggsave(
    png_file,
    plot_object,
    width = 17,
    height = 7.2,
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )

  cat("Models completed: ", length(model_methods), "\n", sep = "")
  cat("Resamples per model: ", cv_folds * cv_repeats, "\n", sep = "")
  cat("Figure: ", png_file, "\n", sep = "")
}

main()
