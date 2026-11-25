options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260722)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(ggbeeswarm)
  library(shapviz)
  library(pROC)
})

project_dir <- "."
analysis_dir <- file.path(project_dir, "03_results", "03_ML", "02_models")
set_labels <- c("Training set" = "Training set", "Validation set" = "Internal hold-out set")
strict_dir <- file.path(project_dir, "03_results", "03_ML", "02_strict_model")
output_dir <- file.path(analysis_dir, "figures")

input_files <- c(
  feature = file.path(
    strict_dir, "01_feature_selection", "rds", "strict_feature_selection.rds"
  ),
  internal_cv = file.path(analysis_dir, "Internal_CV_models.rds"),
  evaluation = file.path(analysis_dir, "Model_evaluation.rds"),
  shap = file.path(analysis_dir, "SHAP", "rds", "PAM_SHAP_objects.rds"),
  de = file.path(project_dir, "03_results", "02_DE", "rds", "DE_objects.rds")
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

module_dirs <- setNames(
  file.path(
    output_dir,
    c(
      "01_feature_selection",
      "02_internal_validation",
      "03_model_evaluation",
      "04_PAM_interpretation"
    )
  ),
  c("feature", "internal_cv", "evaluation", "shap")
)
invisible(lapply(module_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

feature_object <- readRDS(input_files[["feature"]])
cv_object <- readRDS(input_files[["internal_cv"]])
evaluation_object <- readRDS(input_files[["evaluation"]])
shap_object <- readRDS(input_files[["shap"]])
de_object <- readRDS(input_files[["de"]])

stopifnot(
  all(c("parameters", "train_ids", "test_ids", "training_de", "stability", "final_features") %in% names(feature_object)),
  all(c("resample_metrics", "summary") %in% names(cv_object)),
  all(c("predictions", "performance", "roc", "calibration", "dca") %in% names(evaluation_object)),
  all(c("SHAP_values", "SHAP_importance", "representative_samples", "shapviz_object") %in% names(shap_object)),
  all(c("dds", "metadata") %in% names(de_object))
)

final_features <- feature_object$final_features
if (nrow(final_features) != 5L) {
  stop("The locked feature panel does not contain five features.")
}

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

if (!setequal(unique(cv_object$resample_metrics$Model), model_order)) {
  stop("Internal-CV models do not match the locked 10-model design.")
}
if (!setequal(unique(evaluation_object$predictions$Model), model_order)) {
  stop("Evaluation models do not match the locked 10-model design.")
}

parameter_value <- function(parameter_name) {
  value <- feature_object$parameters$value[
    feature_object$parameters$parameter == parameter_name
  ]
  if (length(value) != 1L) {
    stop("Parameter is not uniquely defined: ", parameter_name)
  }
  as.numeric(value)
}

fdr_cutoff <- parameter_value("DE_padj_cutoff")
lfc_cutoff <- parameter_value("DE_abs_log2FC_cutoff")
lasso_cutoff <- parameter_value("LASSO_frequency_cutoff")
boruta_cutoff <- parameter_value("Boruta_frequency_cutoff")

save_plot_pair <- function(plot_object, file_stem, width, height, dpi = 320) {
  ggsave(
    paste0(file_stem, ".pdf"),
    plot_object,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  )
  ggsave(
    paste0(file_stem, ".png"),
    plot_object,
    width = width,
    height = height,
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )
}

theme_article <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.85),
    axis.title = element_text(face = "bold", colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    axis.ticks = element_line(colour = "black", linewidth = 0.45),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey25", size = 10),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.background = element_blank()
  )

theme_model <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(face = "bold", colour = "black", size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65)
  )

training_de <- feature_object$training_de |>
  mutate(
    direction = case_when(
      !is.na(padj) & padj < fdr_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
      !is.na(padj) & padj < fdr_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    minus_log10_padj = -log10(pmax(coalesce(padj, 1), .Machine$double.xmin)),
    final_feature = feature_id %in% final_features$feature_id
  )

training_de$direction <- factor(
  training_de$direction,
  levels = c("Upregulated", "Not significant", "Downregulated")
)

direction_counts <- table(training_de$direction)
direction_labels <- c(
  "Upregulated" = paste0("Upregulated (", direction_counts[["Upregulated"]], ")"),
  "Not significant" = paste0("Not significant (", direction_counts[["Not significant"]], ")"),
  "Downregulated" = paste0("Downregulated (", direction_counts[["Downregulated"]], ")")
)
volcano_colors <- c(
  "Upregulated" = "#B2182B",
  "Not significant" = "#D9D9D9",
  "Downregulated" = "#2166AC"
)

x_abs_limit <- ceiling(max(abs(training_de$log2FoldChange), na.rm = TRUE) * 10) / 10
y_upper <- ceiling(max(training_de$minus_log10_padj, na.rm = TRUE) * 10) / 10

p_volcano <- ggplot(training_de, aes(log2FoldChange, minus_log10_padj)) +
  geom_point(
    data = filter(training_de, direction == "Not significant"),
    aes(colour = direction),
    size = 0.42,
    alpha = 0.16,
    stroke = 0
  ) +
  geom_point(
    data = filter(training_de, direction != "Not significant"),
    aes(colour = direction),
    position = position_jitter(width = 0.018, height = 0.018, seed = 20260722),
    size = 1.85,
    alpha = 0.88,
    stroke = 0
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = 2,
    linewidth = 0.42,
    colour = "black"
  ) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = 4,
    linewidth = 0.42,
    colour = "black"
  ) +
  geom_point(
    data = filter(training_de, final_feature),
    shape = 21,
    size = 2.8,
    stroke = 0.65,
    colour = "black",
    fill = "white"
  ) +
  geom_text_repel(
    data = filter(training_de, final_feature),
    aes(label = feature_id),
    size = 3.0,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.colour = "grey45",
    segment.size = 0.3,
    seed = 20260722,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = volcano_colors, labels = direction_labels, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_cartesian(
    xlim = c(-x_abs_limit, x_abs_limit),
    ylim = c(0, y_upper * 1.05),
    clip = "off"
  ) +
  labs(
    title = "Differential cfRNAs before feature selection",
    subtitle = "Training set: OVX vs Sham",
    x = "log2(Fold change, OVX vs Sham)",
    y = expression(-log[10](adjusted~P~value)),
    colour = NULL
  ) +
  theme_article

save_plot_pair(
  p_volcano,
  file.path(module_dirs[["feature"]], "01_DE_volcano"),
  8.6,
  6.2
)

stability_data <- feature_object$stability |>
  mutate(
    status = case_when(
      stable_lasso & stable_boruta ~ "Final panel",
      stable_lasso ~ "LASSO only",
      stable_boruta ~ "Boruta only",
      TRUE ~ "Not selected"
    ),
    status = factor(
      status,
      levels = c("Not selected", "LASSO only", "Boruta only", "Final panel")
    )
  )

status_counts <- table(stability_data$status)
status_labels <- setNames(
  paste0(names(status_counts), " (", as.integer(status_counts), ")"),
  names(status_counts)
)
stability_colors <- c(
  "Not selected" = "#CFCFCF",
  "LASSO only" = "#74ADD1",
  "Boruta only" = "#F4A582",
  "Final panel" = "#D73027"
)

boruta_axis_max <- max(
  0.35,
  ceiling(max(stability_data$boruta_confirmed_frequency) * 20) / 20 + 0.05
)

set.seed(20260722)
stability_data <- stability_data |>
  mutate(
    plot_x = pmin(
      1.045,
      pmax(0, lasso_selection_frequency + runif(n(), -0.006, 0.006))
    ),
    plot_y = pmin(
      boruta_axis_max - 0.002,
      pmax(0, boruta_confirmed_frequency + runif(n(), -0.004, 0.004))
    )
  )

p_stability <- ggplot(
  stability_data,
  aes(plot_x, plot_y)
) +
  annotate(
    "rect",
    xmin = lasso_cutoff,
    xmax = Inf,
    ymin = boruta_cutoff,
    ymax = Inf,
    fill = "#D73027",
    alpha = 0.04
  ) +
  geom_hline(
    yintercept = boruta_cutoff,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_vline(
    xintercept = lasso_cutoff,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_point(
    aes(fill = status),
    shape = 21,
    size = 3.2,
    colour = "black",
    stroke = 0.35,
    alpha = 0.96
  ) +
  geom_text_repel(
    data = filter(stability_data, final_feature),
    aes(label = feature_id),
    size = 2.75,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.25,
    segment.colour = "grey45",
    segment.size = 0.25,
    seed = 20260722
  ) +
  scale_fill_manual(values = stability_colors, labels = status_labels, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, boruta_axis_max),
    breaks = seq(0, boruta_axis_max, 0.05),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Feature selection stability",
    x = "LASSO selection frequency",
    y = "Boruta confirmation frequency",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "#EAEAEA", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.85),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.margin = margin(8, 24, 8, 8)
  )

save_plot_pair(
  p_stability,
  file.path(module_dirs[["feature"]], "02_feature_stability"),
  7.4,
  6.6
)

frequency_order <- final_features |>
  arrange(desc(lasso_selection_frequency), desc(boruta_confirmed_frequency)) |>
  pull(feature_id)

frequency_data <- final_features |>
  mutate(
    direction = if_else(log2FoldChange > 0, "Upregulated", "Downregulated"),
    feature_id = factor(feature_id, levels = rev(frequency_order))
  )

direction_colors <- c(
  "Upregulated" = "#B2182B",
  "Downregulated" = "#2166AC"
)

make_frequency_plot <- function(value_column, cutoff, title, y_title, y_max) {
  ggplot(
    frequency_data,
    aes(x = feature_id, y = .data[[value_column]], fill = direction)
  ) +
    geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
    geom_hline(
      yintercept = cutoff,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "black"
    ) +
    geom_text(
      aes(label = scales::percent(.data[[value_column]], accuracy = 1)),
      hjust = -0.14,
      size = 3.8
    ) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = direction_colors) +
    scale_y_continuous(
      limits = c(0, y_max),
      labels = scales::label_percent(accuracy = 1),
      expand = c(0, 0)
    ) +
    labs(title = title, x = NULL, y = y_title, fill = NULL) +
    theme_article +
    theme(
      legend.position = "top",
      plot.margin = margin(8, 34, 8, 8),
      axis.text.y = element_text(size = 9.5)
    )
}

p_lasso <- make_frequency_plot(
  "lasso_selection_frequency",
  lasso_cutoff,
  "LASSO selection frequency",
  "Selection frequency",
  1.10
)
p_boruta <- make_frequency_plot(
  "boruta_confirmed_frequency",
  boruta_cutoff,
  "Boruta confirmation frequency",
  "Confirmation frequency",
  max(0.36, max(frequency_data$boruta_confirmed_frequency) + 0.06)
)

save_plot_pair(
  p_lasso,
  file.path(module_dirs[["feature"]], "03_LASSO_frequency"),
  8.0,
  4.8
)
save_plot_pair(
  p_boruta,
  file.path(module_dirs[["feature"]], "04_Boruta_frequency"),
  8.0,
  4.8
)

metadata <- as.data.frame(de_object$metadata)
if (!identical(rownames(metadata), colnames(de_object$dds))) {
  stop("Metadata rows do not match the DESeq2 object columns.")
}
if (!all(final_features$feature_id %in% rownames(de_object$dds))) {
  stop("Selected features are not all present in the DESeq2 object.")
}

vst_all <- varianceStabilizingTransformation(de_object$dds, blind = FALSE)
heatmap_matrix <- assay(vst_all)[final_features$feature_id, , drop = FALSE]
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix[heatmap_matrix > 2.5] <- 2.5
heatmap_matrix[heatmap_matrix < -2.5] <- -2.5

sample_number <- as.integer(gsub("[^0-9]", "", rownames(metadata)))
sample_order <- rownames(metadata)[
  order(factor(metadata$group, levels = c("N", "OP")), sample_number)
]
heatmap_matrix <- heatmap_matrix[, sample_order, drop = FALSE]
stopifnot(!anyDuplicated(c(feature_object$train_ids, feature_object$test_ids)),
          setequal(sample_order, c(feature_object$train_ids, feature_object$test_ids)))

heatmap_annotation <- data.frame(
  Group = factor(metadata[sample_order, "group"], levels = c("N", "OP"), labels = c("Sham", "OVX")),
  Dataset = factor(
    ifelse(sample_order %in% feature_object$train_ids, "Training", "Validation"),
    levels = c("Training", "Validation"), labels = c("Training", "Internal hold-out")
  ),
  row.names = sample_order
)

heatmap_col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#FFFFFF", "#B2182B")
)
heatmap_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Group = heatmap_annotation$Group,
  Dataset = heatmap_annotation$Dataset,
  col = list(
    Group = c(Sham = "#66C2A5", OVX = "#FC8D62"),
    Dataset = c(Training = "#4C72B0", "Internal hold-out" = "#DD8452")
  ),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = grid::unit(3, "mm")
)

heatmap_plot <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-score",
  col = heatmap_col_fun,
  top_annotation = heatmap_top_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = heatmap_annotation$Group,
  cluster_column_slices = FALSE,
  column_gap = grid::unit(3, "mm"),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 9),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(title = "Z-score", at = c(-2, -1, 0, 1, 2)),
  rect_gp = grid::gpar(col = NA)
)

heatmap_pdf <- file.path(module_dirs[["feature"]], "05_selected_feature_heatmap.pdf")
grDevices::pdf(heatmap_pdf, width = 11.2, height = 5.3, useDingbats = FALSE)
ComplexHeatmap::draw(
  heatmap_plot,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  column_title = "Selected cfRNA features",
  column_title_gp = grid::gpar(fontsize = 15, fontface = "bold")
)
grDevices::dev.off()

heatmap_png <- file.path(module_dirs[["feature"]], "05_selected_feature_heatmap.png")
ragg::agg_png(
  heatmap_png,
  width = 11.2,
  height = 5.3,
  units = "in",
  res = 400,
  background = "white"
)
ComplexHeatmap::draw(
  heatmap_plot,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  column_title = "Selected cfRNA features",
  column_title_gp = grid::gpar(fontsize = 15, fontface = "bold")
)
grDevices::dev.off()

cv_long <- cv_object$resample_metrics |>
  select(Model, Resample, Threshold, AUC, Spec, Sens, F1, Accuracy) |>
  pivot_longer(
    cols = c(AUC, Spec, Sens, F1, Accuracy),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Model = factor(Model, levels = rev(model_order)),
    Metric = factor(Metric, levels = c("AUC", "Spec", "Sens", "F1", "Accuracy"))
  )

if (any(!is.finite(cv_long$Value))) {
  stop("Internal-CV metrics contain non-finite values.")
}

cv_summary <- cv_long |>
  group_by(Model, Metric) |>
  summarise(
    N = n(),
    Mean = mean(Value),
    SD = sd(Value),
    CI_low = pmax(0, Mean - 1.96 * SD / sqrt(N)),
    CI_high = pmin(1, Mean + 1.96 * SD / sqrt(N)),
    .groups = "drop"
  ) |>
  mutate(
    Label = sprintf(
      "%.3f +/- %.3f\n95%% CI (%.3f-%.3f)",
      Mean,
      SD,
      CI_low,
      CI_high
    )
  )

p_internal_cv <- ggplot(cv_long, aes(Value, Model)) +
  geom_jitter(
    aes(colour = Model),
    width = 0,
    height = 0.14,
    alpha = 0.16,
    size = 0.48
  ) +
  geom_boxplot(
    aes(fill = Model, colour = Model),
    width = 0.42,
    outlier.shape = NA,
    alpha = 0.26,
    linewidth = 0.45,
    median.colour = "black",
    median.linewidth = 0.55
  ) +
  geom_text(
    data = cv_summary,
    aes(x = 0.018, y = Model, label = Label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.0,
    lineheight = 0.84,
    colour = "black"
  ) +
  facet_grid(cols = vars(Metric)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_colour_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Internal validation of 10 machine-learning models",
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

save_plot_pair(
  p_internal_cv,
  file.path(module_dirs[["internal_cv"]], "01_internal_CV_metrics"),
  17.0,
  6.4
)

prediction_data <- evaluation_object$predictions |>
  mutate(
    Observed = factor(metadata[Sample_ID, "group"], levels = c("N", "OP"), labels = c("Sham", "OVX")),
    Model = factor(Model, levels = model_order),
    Set = factor(Set, levels = c("Training set", "Validation set"))
  )
stopifnot(!anyNA(prediction_data$Observed))
performance_data <- evaluation_object$performance |>
  mutate(
    Model = factor(Model, levels = model_order),
    Set = factor(Set, levels = c("Training set", "Validation set"))
  )

roc_curve_data <- bind_rows(lapply(
  split(prediction_data, interaction(prediction_data$Set, prediction_data$Model, drop = TRUE)),
  function(one) {
    roc_result <- pROC::roc(
      response = one$Observed,
      predictor = one$Probability,
      levels = c("Sham", "OVX"),
      direction = "<",
      quiet = TRUE
    )
    data.frame(
      Set = as.character(one$Set[[1L]]),
      Model = as.character(one$Model[[1L]]),
      FPR = 1 - as.numeric(roc_result$specificities),
      TPR = as.numeric(roc_result$sensitivities)
    )
  }
))

plot_roc <- function(set_name) {
  roc_statistics <- evaluation_object$roc |>
    filter(Set == set_name) |>
    mutate(model_index = match(Model, model_order)) |>
    arrange(desc(AUC), model_index) |>
    mutate(Legend = sprintf("%s  AUC: %.4f", Model, AUC))

  plot_data <- roc_curve_data |>
    filter(Set == set_name) |>
    left_join(select(roc_statistics, Model, Legend), by = "Model") |>
    mutate(Legend = factor(Legend, levels = roc_statistics$Legend))

  legend_colors <- setNames(model_colors[roc_statistics$Model], roc_statistics$Legend)

  ggplot(plot_data, aes(FPR, TPR, colour = Legend, group = Legend)) +
    geom_path(linewidth = 0.95) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
    scale_colour_manual(values = legend_colors) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    labs(
      title = unname(set_labels[set_name]),
      x = "1 - Specificity",
      y = "Sensitivity",
      colour = NULL
    ) +
    theme_model
}

plot_calibration <- function(set_name) {
  brier_labels <- performance_data |>
    filter(Set == set_name) |>
    arrange(Brier) |>
    mutate(Legend = sprintf("%s  Brier: %.3f", Model, Brier))

  plot_data <- evaluation_object$calibration |>
    filter(Set == set_name) |>
    left_join(select(brier_labels, Model, Legend), by = "Model") |>
    mutate(Legend = factor(Legend, levels = brier_labels$Legend))

  rug_data <- prediction_data |>
    filter(Set == set_name) |>
    left_join(select(brier_labels, Model, Legend), by = "Model") |>
    mutate(Legend = factor(Legend, levels = brier_labels$Legend))

  legend_colors <- setNames(model_colors[as.character(brier_labels$Model)], brier_labels$Legend)

  ggplot(
    plot_data,
    aes(Predicted_probability, Observed_frequency, colour = Legend, group = Legend)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 1.0) +
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
    labs(
      title = unname(set_labels[set_name]),
      x = "Predicted probability",
      y = "Observed frequency",
      colour = NULL
    ) +
    theme_model
}

draw_dca <- function(set_name) {
  plot_data <- evaluation_object$dca |>
    filter(Set == set_name)
  dca_colors <- c(model_colors, All = "grey65", None = "black")
  line_order <- c(model_order, "All", "None")

  graphics::par(mar = c(8.0, 5.0, 3.5, 1.0) + 0.1, xpd = FALSE, family = "sans")
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
  graphics::axis(
    1,
    at = seq(0, 1, 0.2),
    labels = sprintf("%.1f", seq(0, 1, 0.2))
  )
  graphics::mtext("High Risk Threshold", side = 1, line = 2.1, cex = 1.15)
  graphics::axis(
    1,
    at = c(0.01, 0.20, 0.40, 0.60, 0.80, 0.99),
    labels = c("1:100", "1:4", "2:3", "3:2", "4:1", "100:1"),
    line = 3.7,
    tick = TRUE
  )
  graphics::mtext("Cost:Benefit Ratio", side = 1, line = 5.9, cex = 1.15)

  for (model_name in line_order) {
    one <- plot_data[plot_data$Model == model_name, , drop = FALSE]
    graphics::lines(
      one$Threshold,
      one$Standardized_net_benefit,
      col = dca_colors[[model_name]],
      lty = 1,
      lwd = ifelse(model_name %in% c("All", "None"), 1.2, 1.6)
    )
  }

  graphics::box()
  graphics::legend(
    "topright",
    legend = line_order,
    col = dca_colors[line_order],
    lty = 1,
    lwd = c(rep(1.6, length(model_order)), 1.2, 1.2),
    bg = "white",
    box.col = "black",
    cex = 0.78,
    y.intersp = 0.85
  )
}

save_dca_pair <- function(set_name, file_stem, width = 10, height = 8) {
  grDevices::cairo_pdf(paste0(file_stem, ".pdf"), width = width, height = height)
  draw_dca(set_name)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(file_stem, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 320,
    background = "white"
  )
  draw_dca(set_name)
  grDevices::dev.off()
}

radar_metrics <- c("Sensitivity", "Specificity", "Accuracy", "F1", "Precision", "NPV")

make_radar_plot <- function(set_name) {
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

  outer_ring_data <- data.frame(Angle = seq(0, 2 * pi, length.out = 361)) |>
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

evaluation_plots <- list(
  "01_ROC_training" = plot_roc("Training set"),
  "02_ROC_validation" = plot_roc("Validation set"),
  "03_calibration_training" = plot_calibration("Training set"),
  "04_calibration_validation" = plot_calibration("Validation set"),
  "07_metrics_training" = make_radar_plot("Training set"),
  "08_metrics_validation" = make_radar_plot("Validation set")
)

for (plot_name in names(evaluation_plots)) {
  save_plot_pair(
    evaluation_plots[[plot_name]],
    file.path(module_dirs[["evaluation"]], plot_name),
    10,
    8
  )
}

save_dca_pair(
  "Training set",
  file.path(module_dirs[["evaluation"]], "05_DCA_training")
)
save_dca_pair(
  "Validation set",
  file.path(module_dirs[["evaluation"]], "06_DCA_validation")
)

common_shap_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 9.5),
    axis.title = element_text(face = "bold", colour = "black", size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65)
  )

importance_order <- shap_object$SHAP_importance |>
  arrange(desc(mean_abs_SHAP)) |>
  pull(feature_id)

importance_data <- shap_object$SHAP_importance |>
  mutate(feature_id = factor(feature_id, levels = rev(importance_order)))

p_shap_importance <- ggplot(importance_data, aes(mean_abs_SHAP, feature_id)) +
  geom_col(fill = "#E64B35", width = 0.66) +
  geom_text(
    aes(label = sprintf("%.3f", mean_abs_SHAP)),
    hjust = -0.15,
    size = 3.2
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "PAM SHAP feature importance",
    subtitle = "Internal hold-out set",
    x = "mean(|SHAP value|)",
    y = NULL
  ) +
  common_shap_theme +
  theme(plot.margin = margin(8, 28, 8, 8))

shap_long <- shap_object$SHAP_values |>
  mutate(feature_id = factor(feature_id, levels = rev(importance_order)))

expression_range <- range(shap_long$scaled_expression, na.rm = TRUE)

p_shap_summary <- ggplot(
  shap_long,
  aes(SHAP_value, feature_id, colour = scaled_expression)
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
    alpha = 0.88,
    size = 2.6
  ) +
  scale_colour_gradient2(
    low = "#3C5488",
    mid = "grey90",
    high = "#E64B35",
    midpoint = 0,
    breaks = expression_range,
    labels = c("Low", "High"),
    name = "Expression\nlevel"
  ) +
  labs(
    title = "PAM SHAP summary",
    subtitle = "Internal hold-out set",
    x = "SHAP value (impact toward OVX prediction)",
    y = NULL
  ) +
  common_shap_theme +
  theme(legend.position = "right")

make_shap_waterfall <- function(observed_class, sample_id) {
  representative <- prediction_data |>
    filter(Model == "PAM", Set == "Validation set", Observed == observed_class,
           Sample_ID == sample_id)
  if (nrow(representative) != 1L) {
    stop("Representative SHAP sample is not unique for: ", observed_class)
  }

  sample_id <- representative$Sample_ID[[1L]]
  row_id <- match(sample_id, rownames(shap_object$shapviz_object$S))
  if (is.na(row_id)) {
    stop("Representative sample is absent from the shapviz object: ", sample_id)
  }

  shapviz::sv_waterfall(
    shap_object$shapviz_object,
    row_id = row_id,
    max_display = ncol(shap_object$shapviz_object$S),
    fill_colors = c("#E64B35", "#3C5488"),
    format_shap = function(x) sprintf("%.3f", x),
    format_feat = function(x) format(signif(x, 4), trim = TRUE),
    show_connection = TRUE,
    show_annotation = TRUE,
    annotation_size = 3.2
  ) +
    labs(
      title = paste0("PAM SHAP waterfall - ", sample_id),
      subtitle = sprintf(
        "%s internal hold-out sample; OVX probability = %.3f",
        observed_class, representative$Probability[[1L]]
      ),
      x = "Cumulative SHAP value",
      y = NULL
    ) +
    common_shap_theme +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_waterfall_op <- make_shap_waterfall("OVX", "OP4")
p_waterfall_n <- make_shap_waterfall("Sham", "N1")

save_plot_pair(
  p_shap_importance,
  file.path(module_dirs[["shap"]], "01_SHAP_importance"),
  7.0,
  5.2
)
save_plot_pair(
  p_shap_summary,
  file.path(module_dirs[["shap"]], "02_SHAP_summary"),
  8.0,
  5.8
)
save_plot_pair(
  p_waterfall_op,
  file.path(module_dirs[["shap"]], "03_SHAP_waterfall_OP"),
  8.2,
  5.8
)
save_plot_pair(
  p_waterfall_n,
  file.path(module_dirs[["shap"]], "04_SHAP_waterfall_N"),
  8.2,
  5.8
)

