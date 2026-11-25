# 用途：单细胞簇标记基因分析。
# 输入：参考Seurat对象。
# 输出：各簇标记基因及汇总表。

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(presto)
})

set.seed(2025)

project <- Sys.getenv("CFRNA_GMF_ROOT", unset = ".")
work_dir <- file.path(project, "CellFreeGMF_work",
                      "09_Set7_GSE169396_human_femoral_head_GMF")
ref_dir <- file.path(work_dir, "GSE169396_human_femoral_head_reference")
table_dir <- file.path(ref_dir, "tables")

rds_path <- file.path(ref_dir, "GSE169396_human_femoral_head_reconstructed_seurat.rds")
if (!file.exists(rds_path)) {
  stop("Missing reconstructed Seurat RDS: ", rds_path)
}

cat("presto_installed=", requireNamespace("presto", quietly = TRUE), "\n")
cat("presto_version=", as.character(packageVersion("presto")), "\n")

obj <- readRDS(rds_path)
if (!"cluster_reconstructed" %in% colnames(obj@meta.data)) {
  stop("Missing cluster_reconstructed in Seurat metadata.")
}

Idents(obj) <- "cluster_reconstructed"

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  verbose = TRUE
)

markers <- markers %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

write.csv(
  markers,
  file.path(table_dir, "GSE169396_cluster_all_positive_markers_presto.csv"),
  row.names = FALSE
)
write.csv(
  top50,
  file.path(table_dir, "GSE169396_cluster_top50_markers.csv"),
  row.names = FALSE
)

summary <- data.frame(
  item = c(
    "marker_method",
    "test_use",
    "only_pos",
    "min_pct",
    "logfc_threshold",
    "n_clusters",
    "n_marker_rows",
    "presto_version"
  ),
  value = c(
    "Seurat FindAllMarkers accelerated by presto",
    "wilcox",
    TRUE,
    0.1,
    0.25,
    length(unique(obj$cluster_reconstructed)),
    nrow(markers),
    as.character(packageVersion("presto"))
  )
)
write.csv(
  summary,
  file.path(table_dir, "GSE169396_FindAllMarkers_presto_summary.csv"),
  row.names = FALSE
)

cat("FindAllMarkers finished.\n")
print(summary)
