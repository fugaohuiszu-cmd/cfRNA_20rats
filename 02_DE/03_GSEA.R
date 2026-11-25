options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
options(clusterProfiler.download.method = "libcurl")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(AnnotationDbi)
  library(dplyr)
})

set.seed(20260714)

project_dir <- "."
result_dir <- file.path(project_dir, "03_results", "02_DE")
input_file <- file.path(result_dir, "tables", "DE_all.csv")
table_dir <- file.path(result_dir, "tables")
rds_dir <- file.path(result_dir, "rds")

if (!file.exists(input_file)) {
  stop("Missing DESeq2 result table: ", input_file)
}

de_df <- read.csv(input_file, check.names = FALSE)
required_columns <- c("feature_id", "RNA_type", "log2FoldChange", "stat", "pvalue", "padj")
missing_columns <- setdiff(required_columns, colnames(de_df))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

rank_df <- de_df %>%
  mutate(
    feature_id = as.character(feature_id),
    ENSEMBL = sub("\\.\\d+$", "", feature_id),
    rank_stat = as.numeric(stat)
  ) %>%
  filter(grepl("^ENSRNOG", ENSEMBL), is.finite(rank_stat))

if (nrow(rank_df) == 0L) {
  stop("No finite Ensembl-ranked features are available for GSEA.")
}

ensembl_keys <- unique(rank_df$ENSEMBL)
entrez_map <- AnnotationDbi::mapIds(
  x = org.Rn.eg.db,
  keys = ensembl_keys,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rank_df$ENTREZID <- unname(entrez_map[rank_df$ENSEMBL])

mapped_df <- rank_df %>%
  filter(!is.na(ENTREZID), ENTREZID != "") %>%
  arrange(desc(abs(rank_stat)), ENTREZID) %>%
  distinct(ENTREZID, .keep_all = TRUE)

if (nrow(mapped_df) < 100L) {
  stop("Too few unique Entrez genes remain after mapping: ", nrow(mapped_df))
}

gene_list <- mapped_df$rank_stat
names(gene_list) <- mapped_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

write.csv(rank_df, file.path(table_dir, "GSEA_gene_map.csv"), row.names = FALSE)
write.csv(mapped_df, file.path(table_dir, "GSEA_ranked_genes.csv"), row.names = FALSE)

gsea_go <- clusterProfiler::gseGO(
  geneList = gene_list,
  ont = "ALL",
  OrgDb = org.Rn.eg.db,
  keyType = "ENTREZID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE,
  by = "fgsea"
)

gsea_kegg <- clusterProfiler::gseKEGG(
  geneList = gene_list,
  organism = "rno",
  keyType = "ncbi-geneid",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  use_internal_data = FALSE,
  seed = TRUE,
  by = "fgsea"
)

go_df <- as.data.frame(gsea_go)
kegg_df <- as.data.frame(gsea_kegg)
if (!"ONTOLOGY" %in% colnames(go_df)) {
  stop("GO GSEA result does not contain the expected ONTOLOGY column.")
}

go_df <- go_df %>%
  mutate(
    contrast = "OP_vs_N",
    direction = ifelse(NES >= 0, "Up_in_OP", "Down_in_OP")
  ) %>%
  relocate(contrast, direction, ONTOLOGY, .before = 1)
kegg_df <- kegg_df %>%
  mutate(
    contrast = "OP_vs_N",
    direction = ifelse(NES >= 0, "Up_in_OP", "Down_in_OP")
  ) %>%
  relocate(contrast, direction, .before = 1)

go_p05 <- go_df %>% filter(!is.na(pvalue), pvalue < 0.05) %>% arrange(pvalue, desc(abs(NES)))
kegg_p05 <- kegg_df %>% filter(!is.na(pvalue), pvalue < 0.05) %>% arrange(pvalue, desc(abs(NES)))
go_fdr05 <- go_df %>% filter(!is.na(p.adjust), p.adjust < 0.05) %>% arrange(p.adjust, desc(abs(NES)))
kegg_fdr05 <- kegg_df %>% filter(!is.na(p.adjust), p.adjust < 0.05) %>% arrange(p.adjust, desc(abs(NES)))

write.csv(go_df, file.path(table_dir, "GSEA_GO.csv"), row.names = FALSE)
write.csv(kegg_df, file.path(table_dir, "GSEA_KEGG.csv"), row.names = FALSE)
write.csv(go_p05, file.path(table_dir, "GSEA_GO_p05.csv"), row.names = FALSE)
write.csv(kegg_p05, file.path(table_dir, "GSEA_KEGG_p05.csv"), row.names = FALSE)
write.csv(go_fdr05, file.path(table_dir, "GSEA_GO_FDR05.csv"), row.names = FALSE)
write.csv(kegg_fdr05, file.path(table_dir, "GSEA_KEGG_FDR05.csv"), row.names = FALSE)

summary_df <- data.frame(
  item = c(
    "DE_features_with_finite_Ensembl_stat",
    "unique_Entrez_genes_for_GSEA",
    "GO_terms_total", "GO_terms_p05", "GO_terms_FDR05",
    "KEGG_terms_total", "KEGG_terms_p05", "KEGG_terms_FDR05"
  ),
  value = c(
    nrow(rank_df), nrow(mapped_df),
    nrow(go_df), nrow(go_p05), nrow(go_fdr05),
    nrow(kegg_df), nrow(kegg_p05), nrow(kegg_fdr05)
  )
)
write.csv(summary_df, file.path(table_dir, "GSEA_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    gene_list = gene_list,
    gene_map = rank_df,
    ranked_genes = mapped_df,
    gsea_go = gsea_go,
    gsea_kegg = gsea_kegg
  ),
  file.path(rds_dir, "GSEA_objects.rds")
)
