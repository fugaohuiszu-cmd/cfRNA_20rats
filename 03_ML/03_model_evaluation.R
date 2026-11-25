options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

project_dir <- "."
result_dir <- file.path(project_dir, "03_results", "03_ML", "02_models")
feature_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "01_feature_selection", "rds", "strict_feature_selection.rds"
)
cv_rds <- file.path(result_dir, "Internal_CV_models.rds")
de_rds <- file.path(
  project_dir, "03_results", "02_DE", "rds", "DE_objects.rds"
)
figure_dir <- file.path(result_dir, "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260717L
set_labels <- c("Training set" = "Training set", "Validation set" = "Internal hold-out set")
bootstrap_n <- 1000L
model_order <- c(
  "RF", "XGB", "PAM", "LGBM", "SVMRAD",
  "GLMNETRidge", "LDA", "SVMLIN", "GBM", "GLMNETLasso"
)
model_colors <- c(
  RF = "#E64B35",
  XGB = "#4DBBD5",
  PAM = "#3C5488",
  LGBM = "#F39B7F",
  SVMRAD = "#8491B4",
  GLMNETRidge = "#91D1C2",
  LDA = "#7E6148",
  SVMLIN = "#B09C85",
  GBM = "#0072B2",
  GLMNETLasso = "#D55E00"
)
radar_metrics <- c("Sensitivity", "Specificity", "Accuracy", "F1", "Precision", "NPV")

classification_metrics <- function(obs, probability, threshold) {
  obs <- as.character(obs)
  probability <- as.numeric(probability)
  true_op <- obs == "OVX"
  pred_op <- probability >= threshold
  tp <- sum(true_op & pred_op)
  fp <- sum(!true_op & pred_op)
  fn <- sum(true_op & !pred_op)
  tn <- sum(!true_op & !pred_op)
  stopifnot(tp + fn > 0L, tn + fp > 0L)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  accuracy <- (tp + tn) / length(obs)
  precision <- if ((tp + fp) == 0L) 0 else tp / (tp + fp)
  npv <- if ((tn + fn) == 0L) 0 else tn / (tn + fn)
  f1 <- if ((precision + sensitivity) == 0) {
    0
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }
  roc_object <- pROC::roc(
    response = factor(obs, levels = c("Sham", "OVX")),
    predictor = probability,
    levels = c("Sham", "OVX"),
    direction = "<",
    quiet = TRUE
  )
  c(
    AUC = as.numeric(pROC::auc(roc_object)),
    Brier = mean((as.integer(true_op) - probability)^2),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    Precision = precision,
    NPV = npv,
    F1 = f1,
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn
  )
}

bootstrap_metrics <- function(obs, probability, threshold, n_boot, bootstrap_seed) {
  set.seed(bootstrap_seed)
  obs <- as.character(obs)
  class_indices <- split(seq_along(obs), factor(obs, levels = c("Sham", "OVX")))
  values <- matrix(
    NA_real_,
    nrow = n_boot,
    ncol = 8L,
    dimnames = list(
      NULL,
      c("AUC", "Brier", "Sensitivity", "Specificity", "Accuracy", "Precision", "NPV", "F1")
    )
  )
  for (iteration in seq_len(n_boot)) {
    index <- unlist(lapply(class_indices, function(x) sample(x, length(x), replace = TRUE)))
    one <- classification_metrics(obs[index], probability[index], threshold)
    values[iteration, ] <- one[colnames(values)]
  }
  data.frame(
    Metric = colnames(values),
    CI_low = apply(values, 2L, stats::quantile, probs = 0.025),
    CI_high = apply(values, 2L, stats::quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
}

predict_locked_model <- function(model_name, model_result, new_data) {
  fit <- model_result$final_model
  if (inherits(fit, "train")) {
    ovx_class <- fit$levels[[1L]]
    stopifnot(identical(rownames(fit$trainingData), train_ids),
              identical(as.character(fit$trainingData$.outcome) == ovx_class, training_labels == "OVX"))
    probability <- predict(
      fit,
      newdata = as.data.frame(new_data, check.names = FALSE),
      type = "prob"
    )[, ovx_class]
    return(as.numeric(probability))
  }
  center <- model_result$preprocessing$center[colnames(new_data)]
  scale_value <- model_result$preprocessing$scale[colnames(new_data)]
  if (any(!is.finite(center)) || any(!is.finite(scale_value)) || any(scale_value == 0)) {
    stop("Invalid locked preprocessing parameters for model: ", model_name)
  }
  scaled_data <- sweep(sweep(as.matrix(new_data), 2L, center, "-"), 2L, scale_value, "/")
  if (model_name == "XGB") {
    return(as.numeric(predict(
      model_result$final_model,
      xgboost::xgb.DMatrix(data = scaled_data)
    )))
  }
  if (model_name == "LGBM") {
    return(as.numeric(predict(model_result$final_model, scaled_data)))
  }
  stop("Unsupported model object: ", model_name)
}

build_roc_data <- function(prediction_data) {
  curves <- list()
  statistics <- list()
  for (set_name in c("Training set", "Validation set")) {
    for (model_name in model_order) {
      one <- prediction_data[
        prediction_data$Set == set_name & prediction_data$Model == model_name,
        , drop = FALSE
      ]
      roc_object <- pROC::roc(
        factor(one$Observed, levels = c("Sham", "OVX")),
        one$Probability,
        levels = c("Sham", "OVX"),
        direction = "<",
        quiet = TRUE
      )
      auc_ci <- as.numeric(pROC::ci.auc(roc_object, conf.level = 0.95))
      stopifnot(all(is.finite(auc_ci)))
      curve <- data.frame(
        FPR = 1 - as.numeric(roc_object$specificities),
        TPR = as.numeric(roc_object$sensitivities),
        Model = model_name,
        Set = set_name,
        stringsAsFactors = FALSE
      )
      curves[[paste(set_name, model_name)]] <- curve
      statistics[[paste(set_name, model_name)]] <- data.frame(
        Set = set_name,
        Model = model_name,
        AUC = as.numeric(pROC::auc(roc_object)),
        AUC_CI_low = auc_ci[1],
        AUC_CI_high = auc_ci[3],
        stringsAsFactors = FALSE
      )
    }
  }
  list(curve = bind_rows(curves), statistics = bind_rows(statistics))
}

build_calibration_data <- function(prediction_data) {
  pieces <- split(
    prediction_data,
    interaction(prediction_data$Set, prediction_data$Model, drop = TRUE)
  )
  bind_rows(lapply(pieces, function(one) {
    one <- one[order(one$Probability), , drop = FALSE]
    n_bins <- min(4L, floor(nrow(one) / 2L))
    one$Bin <- ceiling(seq_len(nrow(one)) * n_bins / nrow(one))
    one |>
      group_by(Set, Model, Bin) |>
      summarise(
        Predicted_probability = mean(Probability),
        Observed_frequency = mean(Observed == "OVX"),
        N = n(),
        .groups = "drop"
      )
  }))
}

build_dca_data <- function(prediction_data) {
  thresholds <- seq(0.01, 0.99, by = 0.01)
  pieces <- list()
  for (set_name in c("Training set", "Validation set")) {
    set_data <- prediction_data[prediction_data$Set == set_name, , drop = FALSE]
    prevalence <- mean(set_data$Observed == "OVX")
    for (model_name in model_order) {
      one <- set_data[set_data$Model == model_name, , drop = FALSE]
      y <- as.integer(one$Observed == "OVX")
      rows <- lapply(thresholds, function(threshold) {
        high_risk <- one$Probability >= threshold
        tp_rate <- sum(high_risk & y == 1L) / length(y)
        fp_rate <- sum(high_risk & y == 0L) / length(y)
        net_benefit <- tp_rate - fp_rate * threshold / (1 - threshold)
        data.frame(
          Set = set_name,
          Model = model_name,
          Threshold = threshold,
          Standardized_net_benefit = net_benefit / prevalence,
          stringsAsFactors = FALSE
        )
      })
      pieces[[paste(set_name, model_name)]] <- bind_rows(rows)
    }
    treat_all <- prevalence - (1 - prevalence) * thresholds / (1 - thresholds)
    pieces[[paste(set_name, "All")]] <- data.frame(
      Set = set_name,
      Model = "All",
      Threshold = thresholds,
      Standardized_net_benefit = treat_all / prevalence,
      stringsAsFactors = FALSE
    )
    pieces[[paste(set_name, "None")]] <- data.frame(
      Set = set_name,
      Model = "None",
      Threshold = thresholds,
      Standardized_net_benefit = 0,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(pieces)
}

make_radar_plot <- function(performance_data, set_name) {
  long_data <- performance_data |>
    filter(Set == set_name) |>
    select(Model, all_of(radar_metrics)) |>
    pivot_longer(all_of(radar_metrics), names_to = "Metric", values_to = "Value") |>
    mutate(
      Model = factor(Model, levels = model_order),
      Metric = factor(Metric, levels = radar_metrics)
    )
  axis_data <- data.frame(
    Metric = factor(radar_metrics, levels = radar_metrics),
    Angle = pi / 2 - 2 * pi * (seq_along(radar_metrics) - 1) / length(radar_metrics)
  )
  plot_data <- long_data |>
    left_join(axis_data, by = "Metric") |>
    mutate(X = Value * cos(Angle), Y = Value * sin(Angle)) |>
    group_by(Model) |>
    group_modify(~ bind_rows(.x, .x[1, , drop = FALSE])) |>
    ungroup()
  inner_ring_data <- expand.grid(
    Radius = c(0.25, 0.50, 0.75),
    Angle = seq(0, 2 * pi, length.out = 361)
  ) |>
    mutate(X = Radius * cos(Angle), Y = Radius * sin(Angle))
  outer_ring_data <- data.frame(
    Angle = seq(0, 2 * pi, length.out = 361)
  ) |>
    mutate(X = cos(Angle), Y = sin(Angle))
  spoke_data <- axis_data |>
    mutate(X = 0, Y = 0, Xend = cos(Angle), Yend = sin(Angle))
  label_data <- axis_data |>
    mutate(
      X = 1.12 * cos(Angle),
      Y = 1.12 * sin(Angle),
      Hjust = case_when(X > 0.15 ~ 0, X < -0.15 ~ 1, TRUE ~ 0.5),
      Vjust = case_when(Y > 0.15 ~ 0, Y < -0.15 ~ 1, TRUE ~ 0.5)
    )
  scale_label_data <- data.frame(
    X = c(0, 0.25, 0.50, 0.75, 1.00),
    Y = 0.035,
    Label = c("0", "0.25", "0.5", "0.75", "1")
  )
  ggplot(plot_data, aes(X, Y, group = Model, colour = Model)) +
    geom_path(
      data = inner_ring_data,
      aes(X, Y, group = Radius),
      inherit.aes = FALSE,
      colour = "grey72",
      linewidth = 0.65,
      linetype = "dashed"
    ) +
    geom_path(
      data = outer_ring_data,
      aes(X, Y),
      inherit.aes = FALSE,
      colour = "grey70",
      linewidth = 1.05
    ) +
    geom_segment(
      data = spoke_data,
      aes(x = X, y = Y, xend = Xend, yend = Yend),
      inherit.aes = FALSE,
      colour = "grey76",
      linewidth = 0.60,
      linetype = "dashed"
    ) +
    geom_path(linewidth = 1.35, alpha = 0.95, na.rm = TRUE) +
    geom_point(size = 2.8, alpha = 0.98, na.rm = TRUE) +
    geom_text(
      data = label_data,
      aes(X, Y, label = Metric, hjust = Hjust, vjust = Vjust),
      inherit.aes = FALSE,
      size = 4.1
    ) +
    geom_text(
      data = scale_label_data,
      aes(X, Y, label = Label),
      inherit.aes = FALSE,
      size = 3.2,
      colour = "grey45",
      hjust = 0
    ) +
    scale_colour_manual(values = model_colors, breaks = model_order) +
    coord_equal(xlim = c(-1.30, 1.30), ylim = c(-1.25, 1.25), clip = "off") +
    labs(title = unname(set_labels[set_name]), colour = NULL) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.text = element_text(size = 10),
      plot.margin = margin(14, 24, 14, 24)
    )
}

save_plot <- function(plot_object, file_name, width, height) {
  ggsave(
    file.path(figure_dir, paste0(file_name, ".pdf")),
    plot_object,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  )
  ggsave(
    file.path(figure_dir, paste0(file_name, ".png")),
    plot_object,
    width = width,
    height = height,
    dpi = 320,
    device = ragg::agg_png,
    bg = "white"
  )
}

feature_object <- readRDS(feature_rds)
cv_object <- readRDS(cv_rds)
de_object <- readRDS(de_rds)

if (!setequal(names(cv_object$model_results), model_order)) {
  stop("The locked model list does not match the expected 10-model design.")
}
if (!identical(cv_object$parameters$train_ids, feature_object$train_ids)) {
  stop("Training sample IDs differ between feature-selection and model objects.")
}
if (!identical(cv_object$parameters$validation_ids, feature_object$test_ids)) {
  stop("Validation sample IDs differ between feature-selection and model objects.")
}

feature_map <- cv_object$feature_map
feature_ids <- feature_map$feature_id
train_ids <- feature_object$train_ids
validation_ids <- feature_object$test_ids
if (!all(feature_ids %in% rownames(feature_object$vst_matrix_train))) {
  stop("One or more locked features are absent from the training VST matrix.")
}

training_data <- t(feature_object$vst_matrix_train[feature_ids, train_ids, drop = FALSE])
colnames(training_data) <- feature_map$model_feature

training_counts <- counts(feature_object$dds_train, normalized = FALSE)
all_counts <- counts(de_object$dds, normalized = FALSE)
if (!all(rownames(training_counts) %in% rownames(all_counts))) {
  stop("Training genes are not fully represented in the complete count matrix.")
}
validation_counts <- all_counts[
  rownames(training_counts), validation_ids, drop = FALSE
]
geometric_means <- exp(rowMeans(log(training_counts)))
validation_metadata <- as.data.frame(colData(de_object$dds))[
  validation_ids, , drop = FALSE
]
validation_metadata$group <- factor(validation_metadata$group, levels = c("N", "OP"))
validation_dds <- DESeqDataSetFromMatrix(
  countData = round(validation_counts),
  colData = validation_metadata,
  design = ~ group
)
validation_dds <- estimateSizeFactors(validation_dds, geoMeans = geometric_means)
dispersionFunction(validation_dds) <- dispersionFunction(feature_object$dds_train)
validation_vst <- varianceStabilizingTransformation(validation_dds, blind = FALSE)
if (!all(feature_ids %in% rownames(validation_vst))) {
  stop("One or more locked features are absent from the validation VST matrix.")
}
validation_data <- t(assay(validation_vst)[feature_ids, validation_ids, drop = FALSE])
colnames(validation_data) <- feature_map$model_feature

groups <- c(N = "Sham", OP = "OVX")
group_by_id <- setNames(unname(groups[as.character(de_object$metadata$group)]), rownames(de_object$metadata))
training_labels <- unname(group_by_id[train_ids])
validation_labels <- unname(group_by_id[validation_ids])
stopifnot(!anyNA(training_labels), !anyNA(validation_labels))

prediction_rows <- list()
for (model_name in model_order) {
  model_result <- cv_object$model_results[[model_name]]
  threshold <- as.numeric(model_result$threshold)
  if (!is.finite(threshold) || threshold < 0 || threshold > 1) {
    stop("Invalid locked threshold for model: ", model_name)
  }
  training_probability <- predict_locked_model(model_name, model_result, training_data)
  validation_probability <- predict_locked_model(model_name, model_result, validation_data)
  prediction_rows[[paste(model_name, "training")]] <- data.frame(
    Model = model_name,
    Set = "Training set",
    Sample_ID = train_ids,
    Observed = training_labels,
    Probability = training_probability,
    Threshold = threshold,
    Predicted = ifelse(training_probability >= threshold, "OVX", "Sham"),
    stringsAsFactors = FALSE
  )
  prediction_rows[[paste(model_name, "validation")]] <- data.frame(
    Model = model_name,
    Set = "Validation set",
    Sample_ID = validation_ids,
    Observed = validation_labels,
    Probability = validation_probability,
    Threshold = threshold,
    Predicted = ifelse(validation_probability >= threshold, "OVX", "Sham"),
    stringsAsFactors = FALSE
  )
}
prediction_data <- bind_rows(prediction_rows)
if (any(!is.finite(prediction_data$Probability))) {
  stop("Non-finite prediction probabilities were generated.")
}
if (any(prediction_data$Probability < 0 | prediction_data$Probability > 1)) {
  stop("Prediction probabilities outside [0, 1] were generated.")
}

performance_rows <- list()
ci_rows <- list()
counter <- 0L
for (set_name in c("Training set", "Validation set")) {
  for (model_name in model_order) {
    counter <- counter + 1L
    one <- prediction_data[
      prediction_data$Set == set_name & prediction_data$Model == model_name,
      , drop = FALSE
    ]
    metrics <- classification_metrics(one$Observed, one$Probability, unique(one$Threshold))
    bootstrap_ci <- bootstrap_metrics(
      one$Observed,
      one$Probability,
      unique(one$Threshold),
      bootstrap_n,
      seed + counter
    )
    if (set_name == "Validation set") {
      counts <- list(
        Sensitivity = c(metrics[["TP"]], metrics[["TP"]] + metrics[["FN"]]),
        Specificity = c(metrics[["TN"]], metrics[["TN"]] + metrics[["FP"]]),
        Accuracy = c(metrics[["TP"]] + metrics[["TN"]], nrow(one))
      )
      for (metric_name in names(counts)) {
        interval <- stats::binom.test(counts[[metric_name]][1], counts[[metric_name]][2])$conf.int
        bootstrap_ci[bootstrap_ci$Metric == metric_name, c("CI_low", "CI_high")] <-
          as.list(as.numeric(interval))
      }
    }
    performance_rows[[paste(set_name, model_name)]] <- data.frame(
      Set = set_name,
      Model = model_name,
      Threshold = unique(one$Threshold),
      N = nrow(one),
      AUC = metrics[["AUC"]],
      Brier = metrics[["Brier"]],
      Sensitivity = metrics[["Sensitivity"]],
      Specificity = metrics[["Specificity"]],
      Accuracy = metrics[["Accuracy"]],
      Precision = metrics[["Precision"]],
      NPV = metrics[["NPV"]],
      F1 = metrics[["F1"]],
      TP = metrics[["TP"]],
      FP = metrics[["FP"]],
      FN = metrics[["FN"]],
      TN = metrics[["TN"]],
      stringsAsFactors = FALSE
    )
    ci_rows[[paste(set_name, model_name)]] <- bootstrap_ci |>
      mutate(Set = set_name, Model = model_name)
  }
}
performance_data <- bind_rows(performance_rows)
confidence_intervals <- bind_rows(ci_rows) |>
  select(Set, Model, Metric, CI_low, CI_high)
performance_with_ci <- performance_data |>
  pivot_longer(
    cols = c(AUC, Brier, Sensitivity, Specificity, Accuracy, Precision, NPV, F1),
    names_to = "Metric",
    values_to = "Estimate"
  ) |>
  left_join(confidence_intervals, by = c("Set", "Model", "Metric")) |>
  arrange(factor(Set, levels = c("Training set", "Validation set")), factor(Model, levels = model_order))

roc_data <- build_roc_data(prediction_data)
calibration_data <- build_calibration_data(prediction_data)
dca_data <- build_dca_data(prediction_data)

plot_roc <- function(set_name) {
  statistics <- roc_data$statistics |>
    filter(Set == set_name) |>
    mutate(Display_order = match(Model, model_order)) |>
    arrange(desc(AUC), Display_order) |>
    mutate(
      Legend = ifelse(
        is.finite(AUC_CI_low) & is.finite(AUC_CI_high),
        sprintf("%s  AUC: %.4f (%.4f-%.4f)", Model, AUC, AUC_CI_low, AUC_CI_high),
        sprintf("%s  AUC: %.4f", Model, AUC)
      )
    )
  plot_data <- roc_data$curve |>
    filter(Set == set_name) |>
    left_join(statistics[, c("Model", "Legend")], by = "Model") |>
    mutate(Legend = factor(Legend, levels = statistics$Legend))
  legend_colors <- setNames(model_colors[statistics$Model], statistics$Legend)
  ggplot(plot_data, aes(FPR, TPR, colour = Legend, group = Legend)) +
    geom_path(linewidth = 0.95) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
    scale_colour_manual(values = legend_colors) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    labs(title = unname(set_labels[set_name]), x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

plot_calibration <- function(set_name) {
  brier_ci <- confidence_intervals |>
    filter(Set == set_name, Metric == "Brier")
  brier_labels <- performance_data |>
    filter(Set == set_name) |>
    left_join(brier_ci[, c("Model", "CI_low", "CI_high")], by = "Model") |>
    arrange(Brier) |>
    mutate(Legend = sprintf("%s  Brier: %.3f (%.3f-%.3f)", Model, Brier, CI_low, CI_high))
  plot_data <- calibration_data |>
    filter(Set == set_name) |>
    left_join(brier_labels[, c("Model", "Legend")], by = "Model") |>
    mutate(Legend = factor(Legend, levels = brier_labels$Legend))
  rug_data <- prediction_data |>
    filter(Set == set_name) |>
    left_join(brier_labels[, c("Model", "Legend")], by = "Model") |>
    mutate(Legend = factor(Legend, levels = brier_labels$Legend))
  legend_colors <- setNames(model_colors[brier_labels$Model], brier_labels$Legend)
  ggplot(plot_data, aes(Predicted_probability, Observed_frequency, colour = Legend, group = Legend)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 1.00) +
    geom_rug(
      data = rug_data,
      aes(x = Probability, y = NULL, colour = Legend),
      inherit.aes = FALSE,
      sides = "b",
      linewidth = 0.35,
      alpha = 0.28
    ) +
    scale_colour_manual(values = legend_colors) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = unname(set_labels[set_name]), x = "Predicted probability", y = "Observed frequency", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

draw_dca <- function(set_name) {
  plot_data <- dca_data |>
    filter(Set == set_name)
  dca_colors <- c(model_colors, All = "grey65", None = "black")
  dca_linetypes <- c(setNames(rep(1, length(model_order)), model_order), All = 1, None = 1)
  graphics::par(mar = c(8.0, 5.0, 3.5, 1.0) + 0.1, xpd = FALSE)
  graphics::plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "",
    ylab = "Standardized Net Benefit",
    main = paste0("Decision curve analysis - ", unname(set_labels[set_name])),
    xaxt = "n",
    bty = "n",
    cex.main = 1.35,
    cex.lab = 1.15,
    cex.axis = 1.0
  )
  graphics::grid(nx = 5, ny = 5, col = "grey90", lty = 1)
  graphics::axis(1, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)))
  graphics::mtext("High Risk Threshold", side = 1, line = 2.1, cex = 1.15)
  graphics::axis(
    1,
    at = c(0.01, 0.20, 0.40, 0.60, 0.80, 0.99),
    labels = c("1:100", "1:4", "2:3", "3:2", "4:1", "100:1"),
    line = 3.7,
    tick = TRUE
  )
  graphics::mtext("Cost:Benefit Ratio", side = 1, line = 5.9, cex = 1.15)
  for (model_name in c(model_order, "All", "None")) {
    one <- plot_data[plot_data$Model == model_name, , drop = FALSE]
    graphics::lines(
      one$Threshold,
      one$Standardized_net_benefit,
      col = dca_colors[[model_name]],
      lty = dca_linetypes[[model_name]],
      lwd = ifelse(model_name %in% c("All", "None"), 1.2, 1.6)
    )
  }
  graphics::box()
  graphics::legend(
    "topright",
    legend = c(model_order, "All", "None"),
    col = dca_colors[c(model_order, "All", "None")],
    lty = dca_linetypes[c(model_order, "All", "None")],
    lwd = c(rep(1.6, length(model_order)), 1.2, 1.2),
    bg = "white",
    box.col = "black",
    cex = 0.78,
    y.intersp = 0.85
  )
}

save_dca <- function(set_name, file_name, width = 10, height = 8) {
  grDevices::cairo_pdf(
    file.path(figure_dir, paste0(file_name, ".pdf")),
    width = width,
    height = height
  )
  draw_dca(set_name)
  grDevices::dev.off()
  ragg::agg_png(
    file.path(figure_dir, paste0(file_name, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 320,
    background = "white"
  )
  draw_dca(set_name)
  grDevices::dev.off()
}

roc_training <- plot_roc("Training set")
roc_validation <- plot_roc("Validation set")
calibration_training <- plot_calibration("Training set")
calibration_validation <- plot_calibration("Validation set")
radar_training <- make_radar_plot(performance_data, "Training set")
radar_validation <- make_radar_plot(performance_data, "Validation set")

save_plot(roc_training, "ROC_training", 10, 8)
save_plot(roc_validation, "ROC_validation", 10, 8)
save_plot(calibration_training, "Calibration_training", 10, 8)
save_plot(calibration_validation, "Calibration_validation", 10, 8)
save_dca("Training set", "DCA_training")
save_dca("Validation set", "DCA_validation")
save_plot(radar_training, "Radar_training", 10, 8)
save_plot(radar_validation, "Radar_validation", 10, 8)

write.csv(
  prediction_data |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "Model_predictions.csv"),
  row.names = FALSE
)
write.csv(
  performance_data |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "Model_performance.csv"),
  row.names = FALSE
)
write.csv(
  performance_with_ci |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "Model_performance_CI.csv"),
  row.names = FALSE
)
write.csv(
  roc_data$statistics |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "ROC_summary.csv"),
  row.names = FALSE
)
write.csv(
  calibration_data |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "Calibration_data.csv"),
  row.names = FALSE
)
write.csv(
  dca_data |> mutate(Set = unname(set_labels[Set])),
  file.path(result_dir, "DCA_data.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Sample_ID = c(train_ids, validation_ids),
    Set = c(rep("Training set", length(train_ids)), rep("Internal hold-out set", length(validation_ids))),
    Group = c(training_labels, validation_labels),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "Sample_split.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    parameters = list(
      seed = seed,
      bootstrap_n = bootstrap_n,
      models = model_order,
      train_ids = train_ids,
      validation_ids = validation_ids,
      thresholds = setNames(
        vapply(cv_object$model_results, function(x) as.numeric(x$threshold), numeric(1)),
        names(cv_object$model_results)
      ),
      validation_normalization = "Training geometric means and training dispersion function"
    ),
    feature_map = feature_map,
    validation_size_factors = sizeFactors(validation_dds),
    predictions = prediction_data,
    performance = performance_data,
    performance_ci = performance_with_ci,
    roc = roc_data$statistics,
    calibration = calibration_data,
    dca = dca_data
  ),
  file.path(result_dir, "Model_evaluation.rds")
)
