# 用途：解释PAM模型的验证集预测
# 输入：特征、模型、DESeq2及评价对象
# 输出：SHAP结果表、对象及图

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

project_dir <- "F:/DE-analysis/final/cfRNA_project"
feature_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "01_feature_selection", "rds", "strict_feature_selection.rds"
)
cv_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "02_internal_CV", "rds", "Internal_CV_models.rds"
)
de_rds <- file.path(
  project_dir, "03_results", "02_DE", "rds", "DE_objects.rds"
)
evaluation_rds <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "03_model_evaluation", "rds", "Model_evaluation.rds"
)
result_dir <- file.path(
  project_dir, "03_results", "03_ML", "02_strict_model",
  "04_PAM_interpretation"
)
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")
rds_dir <- file.path(result_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260717L
set.seed(seed)

required_packages <- c(
  "DESeq2", "caret", "dplyr", "ggplot2", "ggbeeswarm",
  "kernelshap", "shapviz", "pamr", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

required_files <- c(feature_rds, cv_rds, de_rds, evaluation_rds)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

feature_object <- readRDS(feature_rds)
cv_object <- readRDS(cv_rds)
de_object <- readRDS(de_rds)
evaluation_object <- readRDS(evaluation_rds)

required_cv_names <- c("parameters", "feature_map", "model_results")
missing_cv_names <- setdiff(required_cv_names, names(cv_object))
if (length(missing_cv_names) > 0L) {
  stop("Internal CV object is missing: ", paste(missing_cv_names, collapse = ", "))
}
if (!"PAM" %in% names(cv_object$model_results)) {
  stop("PAM is absent from the locked model object.")
}
if (!identical(cv_object$parameters$train_ids, feature_object$train_ids)) {
  stop("Training sample IDs differ between feature-selection and model objects.")
}
if (!identical(cv_object$parameters$validation_ids, feature_object$test_ids)) {
  stop("Validation sample IDs differ between feature-selection and model objects.")
}

pam_result <- cv_object$model_results[["PAM"]]
required_pam_names <- c("final_model", "threshold")
missing_pam_names <- setdiff(required_pam_names, names(pam_result))
if (length(missing_pam_names) > 0L) {
  stop("PAM result is missing: ", paste(missing_pam_names, collapse = ", "))
}
pam_model <- pam_result$final_model
if (!inherits(pam_model, "train") || !identical(pam_model$method, "pam")) {
  stop("The locked PAM model is not a caret PAM model.")
}

feature_map <- cv_object$feature_map
required_map_columns <- c("feature_id", "model_feature")
missing_map_columns <- setdiff(required_map_columns, colnames(feature_map))
if (length(missing_map_columns) > 0L) {
  stop("Feature map is missing: ", paste(missing_map_columns, collapse = ", "))
}
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

all_metadata <- as.data.frame(colData(de_object$dds))
training_group <- as.character(all_metadata[train_ids, "group"])
validation_group <- as.character(all_metadata[validation_ids, "group"])
training_labels <- ifelse(training_group == "OP", "OP_like", "N_like")
validation_labels <- ifelse(validation_group == "OP", "OP_like", "N_like")

predict_pam_probability <- function(object, X, ...) {
  probability <- predict(
    object,
    newdata = as.data.frame(X, check.names = FALSE),
    type = "prob"
  )[, "OP_like"]
  as.numeric(probability)
}

training_probability <- predict_pam_probability(pam_model, training_data)
validation_probability <- predict_pam_probability(pam_model, validation_data)

locked_predictions <- evaluation_object$predictions |>
  filter(Model == "PAM") |>
  arrange(factor(Set, levels = c("Training set", "Validation set")), Sample_ID)
current_predictions <- bind_rows(
  data.frame(
    Set = "Training set",
    Sample_ID = train_ids,
    Observed = training_labels,
    Probability = training_probability
  ),
  data.frame(
    Set = "Validation set",
    Sample_ID = validation_ids,
    Observed = validation_labels,
    Probability = validation_probability
  )
) |>
  arrange(factor(Set, levels = c("Training set", "Validation set")), Sample_ID)

if (!identical(locked_predictions$Sample_ID, current_predictions$Sample_ID) ||
    max(abs(locked_predictions$Probability - current_predictions$Probability)) > 1e-12) {
  stop("PAM predictions do not reproduce the locked model-evaluation results.")
}

cat("[PAM SHAP] Background: Training set,", nrow(training_data), "samples\n")
cat("[PAM SHAP] Explained set: Validation set,", nrow(validation_data), "samples\n")
cat("[PAM SHAP] Features:", ncol(validation_data), "\n")
cat("[PAM SHAP] Method: exact Kernel SHAP\n")

kernel_object <- kernelshap::kernelshap(
  object = pam_model,
  X = as.data.frame(validation_data, check.names = FALSE),
  bg_X = as.data.frame(training_data, check.names = FALSE),
  pred_fun = predict_pam_probability,
  feature_names = model_features,
  exact = TRUE,
  verbose = TRUE,
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
  PAM_OP_like_probability = validation_probability,
  PAM_predicted_class = ifelse(
    validation_probability >= as.numeric(pam_result$threshold),
    "OP_like", "N_like"
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
    Sample_ID, Observed, PAM_predicted_class, PAM_OP_like_probability,
    feature_id, model_feature, VST_expression, scaled_expression, SHAP_value
  )

importance_table <- shap_long |>
  group_by(feature_id, model_feature) |>
  summarise(
    mean_abs_SHAP = mean(abs(SHAP_value)),
    mean_SHAP = mean(SHAP_value),
    expression_SHAP_spearman = suppressWarnings(
      cor(scaled_expression, SHAP_value, method = "spearman")
    ),
    .groups = "drop"
  ) |>
  arrange(desc(mean_abs_SHAP), feature_id) |>
  mutate(SHAP_importance_rank = row_number())

centroid_matrix <- pam_model$finalModel$centroids
if (!all(model_features %in% rownames(centroid_matrix)) ||
    !all(c("OP_like", "N_like") %in% colnames(centroid_matrix))) {
  stop("PAM centroid matrix does not match the locked features and classes.")
}
centroid_table <- data.frame(
  model_feature = model_features,
  PAM_centroid_OP_like = centroid_matrix[model_features, "OP_like"],
  PAM_centroid_N_like = centroid_matrix[model_features, "N_like"],
  stringsAsFactors = FALSE
) |>
  left_join(feature_map, by = "model_feature") |>
  mutate(
    PAM_centroid_difference_OP_minus_N =
      PAM_centroid_OP_like - PAM_centroid_N_like
  ) |>
  select(
    feature_id, model_feature, PAM_centroid_OP_like, PAM_centroid_N_like,
    PAM_centroid_difference_OP_minus_N
  )

pick_representative <- function(prediction_data, observed_class) {
  candidates <- prediction_data |>
    filter(
      Observed == observed_class,
      PAM_predicted_class == observed_class
    )
  if (nrow(candidates) == 0L) {
    stop("No correctly classified validation sample for: ", observed_class)
  }
  target_probability <- median(candidates$PAM_OP_like_probability)
  candidates |>
    mutate(distance_to_group_median = abs(PAM_OP_like_probability - target_probability)) |>
    arrange(distance_to_group_median, Sample_ID) |>
    slice(1L)
}

representative_op <- pick_representative(prediction_table, "OP_like")
representative_n <- pick_representative(prediction_table, "N_like")
representative_samples <- bind_rows(
  mutate(representative_op, Explanation = "Representative OP-like validation sample"),
  mutate(representative_n, Explanation = "Representative N-like validation sample")
)

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
    subtitle = "Validation set",
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
    subtitle = "Positive SHAP values increase OP-like probability",
    x = "SHAP value for OP-like probability",
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
        gsub("_", "-", observed_class, fixed = TRUE)
      ),
      subtitle = paste0(
        sample_id,
        "; predicted OP-like probability = ",
        sprintf("%.3f", probability)
      ),
      x = "SHAP contribution to OP-like probability"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

p_waterfall_op <- make_waterfall(
  representative_op$Sample_ID,
  representative_op$Observed,
  representative_op$PAM_OP_like_probability
)
p_waterfall_n <- make_waterfall(
  representative_n$Sample_ID,
  representative_n$Observed,
  representative_n$PAM_OP_like_probability
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

cat("[PAM SHAP] Baseline:", sprintf("%.8f", baseline), "\n")
cat("[PAM SHAP] Local accuracy max error:", format(local_accuracy_error, scientific = TRUE), "\n")
cat("[PAM SHAP] Representative samples:\n")
print(representative_samples)
cat("[PAM SHAP] Importance ranking:\n")
print(importance_table)
cat("[PAM SHAP] Output:", normalizePath(result_dir, winslash = "/"), "\n")
