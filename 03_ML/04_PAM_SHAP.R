options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(caret)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(kernelshap)
  library(shapviz)
})

project_dir <- "."
analysis_dir <- file.path(project_dir, "03_results", "03_ML", "02_models")
feature_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "01_feature_selection", "rds", "strict_feature_selection.rds"
)
cv_rds <- file.path(analysis_dir, "Internal_CV_models.rds")
de_rds <- file.path(
  project_dir, "03_results", "02_DE", "rds", "DE_objects.rds"
)
evaluation_rds <- file.path(analysis_dir, "Model_evaluation.rds")
result_dir <- file.path(analysis_dir, "SHAP")
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
rds_dir <- file.path(result_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260717L
set.seed(seed)

feature_object <- readRDS(feature_rds)
cv_object <- readRDS(cv_rds)
de_object <- readRDS(de_rds)
evaluation_object <- readRDS(evaluation_rds)

if (!identical(cv_object$parameters$train_ids, feature_object$train_ids)) {
  stop("Training sample IDs differ between feature-selection and model objects.")
}
if (!identical(cv_object$parameters$validation_ids, feature_object$test_ids)) {
  stop("Validation sample IDs differ between feature-selection and model objects.")
}

pam_result <- cv_object$model_results[["PAM"]]
pam_model <- pam_result$final_model
if (!inherits(pam_model, "train") || !identical(pam_model$method, "pam")) {
  stop("The locked PAM model is not a caret PAM model.")
}

feature_map <- cv_object$feature_map
stopifnot(all(c("feature_id", "model_feature") %in% names(feature_map)))
if (anyDuplicated(feature_map$feature_id) || anyDuplicated(feature_map$model_feature)) {
  stop("Feature map contains duplicated feature names.")
}

feature_ids <- feature_map$feature_id
model_features <- feature_map$model_feature
train_ids <- feature_object$train_ids
validation_ids <- feature_object$test_ids

if (!all(feature_ids %in% rownames(feature_object$vst_matrix_train))) {
  stop("One or more locked features are absent from the training VST matrix.")
}

training_data <- t(
  feature_object$vst_matrix_train[feature_ids, train_ids, drop = FALSE]
)
colnames(training_data) <- model_features

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
validation_data <- t(
  assay(validation_vst)[feature_ids, validation_ids, drop = FALSE]
)
colnames(validation_data) <- model_features

groups <- c(N = "Sham", OP = "OVX")
group_by_id <- setNames(unname(groups[as.character(de_object$metadata$group)]), rownames(de_object$metadata))
training_labels <- unname(group_by_id[train_ids])
validation_labels <- unname(group_by_id[validation_ids])
stopifnot(!anyNA(training_labels), !anyNA(validation_labels))
ovx_class <- pam_model$levels[[1L]]
sham_class <- pam_model$levels[[2L]]
stopifnot(identical(rownames(pam_model$trainingData), train_ids),
          identical(as.character(pam_model$trainingData$.outcome) == ovx_class, training_labels == "OVX"))

predict_pam_probability <- function(object, X, ...) {
  probability <- predict(
    object,
    newdata = as.data.frame(X, check.names = FALSE),
    type = "prob"
  )[, ovx_class]
  as.numeric(probability)
}

training_probability <- predict_pam_probability(pam_model, training_data)
validation_probability <- predict_pam_probability(pam_model, validation_data)

locked <- subset(evaluation_object$predictions, Model == "PAM")
sample_ids <- c(train_ids, validation_ids)
stopifnot(!anyDuplicated(locked$Sample_ID), setequal(locked$Sample_ID, sample_ids))
expected <- locked$Probability[match(sample_ids, locked$Sample_ID)]
stopifnot(max(abs(expected - c(training_probability, validation_probability))) < 1e-12)

kernel_object <- kernelshap::kernelshap(
  object = pam_model,
  X = as.data.frame(validation_data, check.names = FALSE),
  bg_X = as.data.frame(training_data, check.names = FALSE),
  pred_fun = predict_pam_probability,
  feature_names = model_features,
  exact = TRUE,
  verbose = FALSE,
  seed = seed
)

shap_matrix <- as.matrix(kernel_object$S)
rownames(shap_matrix) <- rownames(validation_data)
colnames(shap_matrix) <- colnames(validation_data)
baseline <- as.numeric(kernel_object$baseline)
reconstructed_probability <- baseline + rowSums(shap_matrix)
local_accuracy_error <- max(abs(validation_probability - reconstructed_probability))
if (!is.finite(local_accuracy_error) || local_accuracy_error > 1e-8) {
  stop("SHAP local accuracy check failed: ", local_accuracy_error)
}

preprocess_mean <- pam_model$preProcess$mean[model_features]
preprocess_sd <- pam_model$preProcess$std[model_features]
if (any(!is.finite(preprocess_mean)) || any(!is.finite(preprocess_sd)) ||
    any(preprocess_sd <= 0)) {
  stop("Invalid locked PAM preprocessing parameters.")
}
scaled_validation <- sweep(
  sweep(validation_data, 2L, preprocess_mean, "-"),
  2L, preprocess_sd, "/"
)

prediction_table <- data.frame(
  Sample_ID = validation_ids,
  Observed = validation_labels,
  PAM_OVX_probability = validation_probability,
  PAM_predicted_class = ifelse(
    validation_probability >= as.numeric(pam_result$threshold),
    "OVX", "Sham"
  ),
  stringsAsFactors = FALSE
)

shap_long <- data.frame(
  Sample_ID = rep(rownames(shap_matrix), times = ncol(shap_matrix)),
  model_feature = rep(colnames(shap_matrix), each = nrow(shap_matrix)),
  SHAP_value = as.vector(shap_matrix),
  VST_expression = as.vector(validation_data),
  scaled_expression = as.vector(scaled_validation),
  stringsAsFactors = FALSE
) |>
  left_join(feature_map, by = "model_feature") |>
  left_join(prediction_table, by = "Sample_ID") |>
  select(
    Sample_ID, Observed, PAM_predicted_class, PAM_OVX_probability,
    feature_id, model_feature, VST_expression, scaled_expression, SHAP_value
  )

importance_table <- shap_long |>
  group_by(feature_id, model_feature) |>
  summarise(
    mean_abs_SHAP = mean(abs(SHAP_value)),
    mean_SHAP = mean(SHAP_value),
    expression_SHAP_spearman = cor(scaled_expression, SHAP_value, method = "spearman"),
    .groups = "drop"
  ) |>
  arrange(desc(mean_abs_SHAP), feature_id) |>
  mutate(SHAP_importance_rank = row_number())

centroid_matrix <- pam_model$finalModel$centroids
if (!all(model_features %in% rownames(centroid_matrix)) ||
    !all(c(ovx_class, sham_class) %in% colnames(centroid_matrix))) {
  stop("PAM centroid matrix does not match the locked features and classes.")
}
centroid_table <- data.frame(
  model_feature = model_features,
  PAM_centroid_OVX = centroid_matrix[model_features, ovx_class],
  PAM_centroid_Sham = centroid_matrix[model_features, sham_class],
  stringsAsFactors = FALSE
) |>
  left_join(feature_map, by = "model_feature") |>
  mutate(
    PAM_centroid_difference_OVX_minus_Sham =
      PAM_centroid_OVX - PAM_centroid_Sham
  ) |>
  select(
    feature_id, model_feature, PAM_centroid_OVX, PAM_centroid_Sham,
    PAM_centroid_difference_OVX_minus_Sham
  )

pick_representative <- function(prediction_data, observed_class, sample_id) {
  candidates <- prediction_data |>
    filter(
      Observed == observed_class,
      PAM_predicted_class == observed_class
    )
  if (nrow(candidates) == 0L) {
    stop("No correctly classified internal hold-out sample for ", observed_class)
  }
  target_probability <- median(candidates$PAM_OVX_probability)
  selected <- candidates |>
    mutate(distance_to_group_median = abs(PAM_OVX_probability - target_probability)) |>
    filter(Sample_ID == sample_id)
  stopifnot(nrow(selected) == 1L)
  selected
}

representative_op <- pick_representative(prediction_table, "OVX", "OP4")
representative_n <- pick_representative(prediction_table, "Sham", "N1")
representative_samples <- bind_rows(representative_op, representative_n)

display_order <- importance_table$feature_id
plot_long <- shap_long |>
  mutate(feature_id = factor(feature_id, levels = rev(display_order)))

theme_shap <- theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(color = "black", face = "bold", size = 10),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, color = "#333333"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line = element_blank(),
    legend.title = element_text(face = "bold")
  )

importance_plot_data <- importance_table |>
  mutate(feature_id = factor(feature_id, levels = rev(display_order)))

p_importance <- ggplot(
  importance_plot_data,
  aes(x = mean_abs_SHAP, y = feature_id)
) +
  geom_col(fill = "#9C2F5C", width = 0.66) +
  labs(
    title = "PAM SHAP feature importance",
    subtitle = "Internal hold-out set",
    x = "Mean |SHAP value|",
    y = NULL
  ) +
  theme_shap

p_summary <- ggplot(
  plot_long,
  aes(x = SHAP_value, y = feature_id, colour = scaled_expression)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.45
  ) +
  ggbeeswarm::geom_quasirandom(
    groupOnX = FALSE,
    orientation = "y",
    width = 0.22,
    alpha = 0.90,
    size = 2.6
  ) +
  scale_colour_gradient2(
    low = "#3C5488",
    mid = "#F2F2F2",
    high = "#E64B35",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +
  labs(
    title = "PAM SHAP summary",
    subtitle = "Positive SHAP values increase OVX probability",
    x = "SHAP value for OVX probability",
    y = NULL
  ) +
  theme_shap +
  theme(legend.position = "right")

display_shap_matrix <- shap_matrix
display_validation_data <- validation_data
colnames(display_shap_matrix) <- feature_map$feature_id[
  match(colnames(display_shap_matrix), feature_map$model_feature)
]
colnames(display_validation_data) <- feature_map$feature_id[
  match(colnames(display_validation_data), feature_map$model_feature)
]
shapviz_object <- shapviz::shapviz(
  display_shap_matrix,
  X = as.data.frame(display_validation_data, check.names = FALSE),
  baseline = baseline
)

make_waterfall <- function(sample_id, observed_class, probability) {
  row_index <- match(sample_id, rownames(display_shap_matrix))
  if (is.na(row_index)) {
    stop("Waterfall sample is absent from the SHAP matrix: ", sample_id)
  }
  shapviz::sv_waterfall(
    shapviz_object,
    row_id = row_index,
    max_display = ncol(display_shap_matrix),
    fill_colors = c("#E64B35", "#3C5488"),
    format_shap = function(x) sprintf("%.3f", x),
    format_feat = function(x) sprintf("%.3g", x)
  ) +
    labs(
      title = paste0(
        "PAM SHAP waterfall: ",
        observed_class
      ),
      subtitle = paste0(
        sample_id,
        "; predicted OVX probability = ",
        sprintf("%.3f", probability)
      ),
      x = "SHAP contribution to OVX probability"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

p_waterfall_op <- make_waterfall(
  representative_op$Sample_ID,
  representative_op$Observed,
  representative_op$PAM_OVX_probability
)
p_waterfall_n <- make_waterfall(
  representative_n$Sample_ID,
  representative_n$Observed,
  representative_n$PAM_OVX_probability
)

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

write.csv(
  prediction_table,
  file.path(table_dir, "PAM_validation_predictions.csv"),
  row.names = FALSE
)
write.csv(
  shap_long,
  file.path(table_dir, "PAM_SHAP_values.csv"),
  row.names = FALSE
)
write.csv(
  importance_table,
  file.path(table_dir, "PAM_SHAP_importance.csv"),
  row.names = FALSE
)
write.csv(
  centroid_table,
  file.path(table_dir, "PAM_centroids.csv"),
  row.names = FALSE
)
write.csv(
  representative_samples,
  file.path(table_dir, "PAM_representative_samples.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    parameters = list(
      seed = seed,
      model = "PAM",
      SHAP_method = "exact Kernel SHAP",
      background_set = "Training set",
      explained_set = "Validation set",
      n_background = nrow(training_data),
      n_explained = nrow(validation_data),
      n_features = ncol(validation_data),
      baseline = baseline,
      local_accuracy_max_error = local_accuracy_error,
      package_versions = vapply(
        c("kernelshap", "shapviz", "pamr", "caret", "ggplot2"),
        function(package_name) as.character(packageVersion(package_name)),
        character(1)
      )
    ),
    feature_map = feature_map,
    predictions = prediction_table,
    SHAP_values = shap_long,
    SHAP_importance = importance_table,
    PAM_centroids = centroid_table,
    representative_samples = representative_samples,
    kernelshap_object = kernel_object,
    shapviz_object = shapviz_object
  ),
  file.path(rds_dir, "PAM_SHAP_objects.rds")
)

save_plot(p_importance, "PAM_SHAP_importance", 7.0, 5.2)
save_plot(p_summary, "PAM_SHAP_summary", 8.0, 5.8)
save_plot(p_waterfall_op, "PAM_SHAP_waterfall_OP", 8.2, 5.8)
save_plot(p_waterfall_n, "PAM_SHAP_waterfall_N", 8.2, 5.8)
