# 用途：聚类分辨率比较。
# 输入：参考Seurat对象。
# 输出：分辨率比较、标记基因及细胞簇评估表。

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(presto)
})

set.seed(2025)

project <- Sys.getenv("CFRNA_GMF_ROOT", unset = ".")
work_dir <- file.path(project, "CellFreeGMF_work",
                      "09_Set7_GSE169396_human_femoral_head_GMF")
ref_dir <- file.path(work_dir, "GSE169396_human_femoral_head_reference")
out_dir <- file.path(ref_dir, "resolution_assessment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "markers"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(ref_dir, "GSE169396_human_femoral_head_reconstructed_seurat.rds")
if (!file.exists(rds_path)) {
  stop("Missing reconstructed Seurat RDS: ", rds_path)
}

cat("presto_installed=", requireNamespace("presto", quietly = TRUE), "\n")
cat("presto_version=", as.character(packageVersion("presto")), "\n")

obj <- readRDS(rds_path)
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
if (length(res_cols) == 0) {
  stop("No RNA_snn_res.* columns found in Seurat metadata.")
}

parse_res <- function(x) as.numeric(sub("^RNA_snn_res\\.", "", x))
res_tbl <- data.frame(
  resolution = vapply(res_cols, parse_res, numeric(1)),
  cluster_col = res_cols,
  stringsAsFactors = FALSE
) %>%
  arrange(resolution)

marker_sets <- list(
  Granulocytes = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "CEACAM8", "MPO", "ELANE"),
  T_cells = c("CD3D", "CD3E", "CD3G", "TRAC", "CD2", "IL7R", "CD8A", "CD8B"),
  Monocytes = c("CD14", "LYZ", "LST1", "FCN1", "S100A8", "S100A9", "MS4A7", "CTSS"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD74", "HLA-DRA"),
  Red_blood_cells = c("HBA1", "HBA2", "HBB", "HBD", "GYPA", "ALAS2"),
  Osteoblastic_lineage_cells = c("COL1A1", "COL1A2", "SPARC", "BGLAP", "RUNX2", "ALPL", "IBSP", "SPP1", "PRSS23", "MXRA8"),
  Endothelial_cells = c("PECAM1", "VWF", "KDR", "CLDN5", "EMCN", "RAMP2", "ESAM"),
  Endothelial_progenitor_cells = c("KIT", "PROM1", "CD34", "FLT1", "KDR", "ENG"),
  Plasmacytoid_dendritic_cells = c("GZMB", "IL3RA", "CLEC4C", "LILRA4", "TCF4", "IRF7")
)

norm_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

score_resolution <- function(resolution, cluster_col) {
  clusters <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
  cluster_sizes <- table(as.character(obj@meta.data[[cluster_col]]))

  rows <- list()
  for (cl in clusters) {
    cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[cluster_col]]) == cl]
    for (set_name in names(marker_sets)) {
      genes <- intersect(marker_sets[[set_name]], rownames(norm_mat))
      score <- if (length(genes) > 0) {
        mean(Matrix::colMeans(norm_mat[genes, cells, drop = FALSE]))
      } else {
        NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        resolution = resolution,
        cluster = cl,
        cell_type_marker_set = set_name,
        n_markers_present = length(genes),
        mean_log_normalized_expr = score,
        stringsAsFactors = FALSE
      )
    }
  }
  long <- bind_rows(rows)
  wide <- pivot_wider(
    long,
    id_cols = c(resolution, cluster),
    names_from = cell_type_marker_set,
    values_from = mean_log_normalized_expr
  )
  score_mat <- as.matrix(wide[, names(marker_sets), drop = FALSE])
  rownames(score_mat) <- wide$cluster
  score_z <- apply(score_mat, 2, function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
      rep(0, length(x))
    } else {
      as.numeric(scale(x))
    }
  })
  score_z <- as.data.frame(score_z)
  score_z$resolution <- resolution
  score_z$cluster <- rownames(score_mat)
  type_cols <- names(marker_sets)
  score_numeric <- as.matrix(score_z[, type_cols, drop = FALSE])
  ord <- t(apply(score_numeric, 1, order, decreasing = TRUE))
  best_idx <- ord[, 1]
  second_idx <- ord[, 2]

  cluster_anno <- data.frame(
    resolution = resolution,
    cluster = score_z$cluster,
    n_cells = as.integer(cluster_sizes[score_z$cluster]),
    best_marker_type = type_cols[best_idx],
    best_score_z = score_numeric[cbind(seq_len(nrow(score_numeric)), best_idx)],
    second_marker_type = type_cols[second_idx],
    second_score_z = score_numeric[cbind(seq_len(nrow(score_numeric)), second_idx)],
    marker_margin = score_numeric[cbind(seq_len(nrow(score_numeric)), best_idx)] -
      score_numeric[cbind(seq_len(nrow(score_numeric)), second_idx)],
    stringsAsFactors = FALSE
  )
  list(long = long, score_z = score_z, cluster_anno = cluster_anno)
}

all_long <- list()
all_score_z <- list()
all_cluster_anno <- list()
for (i in seq_len(nrow(res_tbl))) {
  cat("Scoring resolution", res_tbl$resolution[[i]], "\n")
  x <- score_resolution(res_tbl$resolution[[i]], res_tbl$cluster_col[[i]])
  all_long[[i]] <- x$long
  all_score_z[[i]] <- x$score_z
  all_cluster_anno[[i]] <- x$cluster_anno
}
marker_scores_long <- bind_rows(all_long)
marker_scores_z <- bind_rows(all_score_z)
cluster_anno_all <- bind_rows(all_cluster_anno)

write.csv(
  marker_scores_long,
  file.path(out_dir, "GSE169396_resolution_marker_scores_long.csv"),
  row.names = FALSE
)
write.csv(
  marker_scores_z,
  file.path(out_dir, "GSE169396_resolution_marker_scores_z.csv"),
  row.names = FALSE
)
write.csv(
  cluster_anno_all,
  file.path(out_dir, "GSE169396_resolution_cluster_marker_type_calls.csv"),
  row.names = FALSE
)

size_summary <- cluster_anno_all %>%
  group_by(resolution) %>%
  summarise(
    n_clusters = n(),
    min_cluster_size = min(n_cells),
    q10_cluster_size = as.numeric(quantile(n_cells, 0.10)),
    median_cluster_size = median(n_cells),
    max_cluster_size = max(n_cells),
    n_clusters_lt_50 = sum(n_cells < 50),
    n_clusters_lt_100 = sum(n_cells < 100),
    n_clusters_lt_200 = sum(n_cells < 200),
    n_best_marker_types = n_distinct(best_marker_type),
    n_duplicate_marker_type_clusters = n_clusters - n_best_marker_types,
    median_marker_margin = median(marker_margin, na.rm = TRUE),
    n_ambiguous_margin_lt_0p25 = sum(marker_margin < 0.25, na.rm = TRUE),
    n_ambiguous_margin_lt_0p5 = sum(marker_margin < 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    distance_to_original_16 = abs(n_clusters - 16),
    overcluster_flag = n_clusters >= 24 | n_clusters_lt_100 > 0,
    undercluster_flag = n_clusters < 9 | n_best_marker_types < 8
  ) %>%
  arrange(resolution)

write.csv(
  size_summary,
  file.path(out_dir, "GSE169396_resolution_size_marker_summary.csv"),
  row.names = FALSE
)

candidate_res <- size_summary %>%
  filter(
    n_clusters >= 13,
    n_clusters <= 24,
    n_clusters_lt_100 == 0
  ) %>%
  arrange(distance_to_original_16, n_duplicate_marker_type_clusters, desc(median_marker_margin)) %>%
  pull(resolution)

if (length(candidate_res) == 0) {
  candidate_res <- size_summary %>%
    arrange(distance_to_original_16, n_duplicate_marker_type_clusters, desc(median_marker_margin)) %>%
    slice_head(n = 4) %>%
    pull(resolution)
}
candidate_res <- unique(candidate_res)
cat("Candidate resolutions for marker testing:", paste(candidate_res, collapse = ", "), "\n")

marker_quality <- list()
for (res in candidate_res) {
  col <- res_tbl$cluster_col[match(res, res_tbl$resolution)]
  label <- paste0("res", gsub("\\.", "p", as.character(res)))
  cat("Running FindAllMarkers for", label, "\n")
  obj$resolution_assessment_cluster <- as.character(obj@meta.data[[col]])
  Idents(obj) <- "resolution_assessment_cluster"

  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = TRUE
  ) %>%
    arrange(cluster, p_val_adj, desc(avg_log2FC))

  top50 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
    ungroup()

  write.csv(
    markers,
    file.path(out_dir, "markers", paste0("GSE169396_", label, "_all_positive_markers_presto.csv")),
    row.names = FALSE
  )
  write.csv(
    top50,
    file.path(out_dir, "markers", paste0("GSE169396_", label, "_top50_markers.csv")),
    row.names = FALSE
  )

  top_per_cluster <- markers %>%
    group_by(cluster) %>%
    slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
    ungroup()
  top10 <- markers %>%
    group_by(cluster) %>%
    slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
    ungroup()

  marker_quality[[label]] <- data.frame(
    resolution = res,
    n_clusters = length(unique(obj$resolution_assessment_cluster)),
    n_marker_rows = nrow(markers),
    n_clusters_without_markers = length(unique(obj$resolution_assessment_cluster)) -
      length(unique(markers$cluster)),
    median_top_marker_avg_log2FC = median(top_per_cluster$avg_log2FC, na.rm = TRUE),
    median_top_marker_pct_diff = median(top_per_cluster$pct.1 - top_per_cluster$pct.2, na.rm = TRUE),
    median_top10_avg_log2FC = median(top10$avg_log2FC, na.rm = TRUE),
    median_top10_pct_diff = median(top10$pct.1 - top10$pct.2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

marker_quality <- bind_rows(marker_quality)
write.csv(
  marker_quality,
  file.path(out_dir, "GSE169396_resolution_marker_quality_summary.csv"),
  row.names = FALSE
)

combined_summary <- size_summary %>%
  left_join(marker_quality, by = c("resolution", "n_clusters")) %>%
  mutate(
    recommendation_rank = rank(
      distance_to_original_16 * 2 +
        n_clusters_lt_100 * 5 +
        n_ambiguous_margin_lt_0p5 * 0.5 +
        n_duplicate_marker_type_clusters * 0.5 -
        coalesce(median_top10_pct_diff, 0),
      ties.method = "min"
    )
  ) %>%
  arrange(recommendation_rank, resolution)

write.csv(
  combined_summary,
  file.path(out_dir, "GSE169396_resolution_combined_recommendation_table.csv"),
  row.names = FALSE
)

png(file.path(out_dir, "figures", "GSE169396_resolution_cluster_count.png"), width = 2200, height = 1500, res = 180)
print(
  ggplot(size_summary, aes(x = resolution, y = n_clusters)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 16, linetype = "dashed", color = "red") +
    theme_bw() +
    labs(x = "Resolution", y = "Number of clusters", title = "GSE169396 cluster count across resolutions")
)
dev.off()

png(file.path(out_dir, "figures", "GSE169396_resolution_marker_margin.png"), width = 2200, height = 1500, res = 180)
print(
  ggplot(cluster_anno_all, aes(x = factor(resolution), y = marker_margin)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Resolution", y = "Top1 - Top2 canonical marker score z", title = "Canonical marker separability across resolutions")
)
dev.off()

cat("Resolution assessment finished.\n")
cat("Output:", out_dir, "\n")
print(combined_summary)
