# 用途：绘制选定通路的富集曲线
# 输入：GSEA对象和富集汇总绘图数据
# 输出：通路富集曲线

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(patchwork)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
result_dir <- file.path(project_dir, "03_results", "02_DE")
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")
single_curve_dir <- file.path(figure_dir, "GSEA_curves")
go_curve_dir <- file.path(single_curve_dir, "GO")
kegg_curve_dir <- file.path(single_curve_dir, "KEGG")
gsea_rds_file <- file.path(result_dir, "rds", "GSEA_objects.rds")
go_plot_file <- file.path(table_dir, "GSEA_GO_plot.csv")
kegg_plot_file <- file.path(table_dir, "GSEA_KEGG_plot.csv")

required_files <- c(gsea_rds_file, go_plot_file, kegg_plot_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing GSEA files: ", paste(missing_files, collapse = ", "))
}

dir.create(go_curve_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kegg_curve_dir, recursive = TRUE, showWarnings = FALSE)

gsea_objects <- readRDS(gsea_rds_file)
if (!is.list(gsea_objects) ||
    !all(c("gsea_go", "gsea_kegg", "gene_list") %in% names(gsea_objects))) {
  stop("GSEA_objects.rds does not contain gsea_go, gsea_kegg, and gene_list.")
}

go_plot_df <- read.csv(go_plot_file, check.names = FALSE)
kegg_plot_df <- read.csv(kegg_plot_file, check.names = FALSE)
required_columns <- c("module", "module_rank", "ID", "Description", "NES", "p.adjust")
if (!all(required_columns %in% colnames(go_plot_df)) ||
    !all(required_columns %in% colnames(kegg_plot_df))) {
  stop("GSEA plot tables do not contain the expected module selection columns.")
}

go_module_levels <- c("Splicing", "Cytoskeleton", "ECM", "Bone", "Calcium")
kegg_module_levels <- c("Mineral", "Hormone", "Inflammation", "Immune", "Signaling")

select_representative_terms <- function(plot_df, module_levels, source_name) {
  selected <- plot_df %>%
    filter(module_rank == 1L, module %in% module_levels) %>%
    mutate(module = factor(module, levels = module_levels)) %>%
    arrange(module) %>%
    mutate(module = as.character(module), source = source_name)

  if (nrow(selected) != length(module_levels) ||
      !setequal(selected$module, module_levels) ||
      anyDuplicated(selected$module)) {
    stop(source_name, " does not contain one rank-1 term for each module.")
  }

  selected
}

go_terms <- select_representative_terms(go_plot_df, go_module_levels, "GO")
kegg_terms <- select_representative_terms(kegg_plot_df, kegg_module_levels, "KEGG")
curve_terms <- bind_rows(go_terms, kegg_terms) %>%
  select(source, module, ID, Description, NES, p.adjust)

write.csv(
  curve_terms,
  file.path(table_dir, "GSEA_curve_terms.csv"),
  row.names = FALSE
)

go_colors <- c(
  Splicing = "#C94C3C",
  Cytoskeleton = "#4C78A8",
  ECM = "#E69F00",
  Bone = "#2E8B57",
  Calcium = "#7A5195"
)
kegg_colors <- c(
  Mineral = "#4C78A8",
  Hormone = "#CC79A7",
  Inflammation = "#D95F02",
  Immune = "#2E8B57",
  Signaling = "#7A5195"
)

make_single_curve <- function(gsea_result, term_row, module_color, base_size = 9) {
  plot_title <- sprintf(
    "%s | %s\nNES = %.2f; FDR = %.2g",
    term_row$module,
    term_row$Description,
    term_row$NES,
    term_row$p.adjust
  )

  plot_parts <- enrichplot::gseaplot2(
    gsea_result,
    geneSetID = term_row$ID,
    title = plot_title,
    color = module_color,
    base_size = base_size,
    rel_heights = c(1.5, 0.5, 1),
    subplots = 1:3,
    pvalue_table = FALSE
  )

  patchwork::wrap_plots(
    plotlist = plot_parts,
    ncol = 1,
    heights = c(1.5, 0.5, 1)
  )
}

save_single_curves <- function(gsea_result, term_df, module_colors,
                               output_dir, source_name) {
  output_files <- character(nrow(term_df))

  for (i in seq_len(nrow(term_df))) {
    term_row <- term_df[i, , drop = FALSE]
    curve_plot <- make_single_curve(
      gsea_result,
      term_row,
      unname(module_colors[term_row$module]),
      base_size = 11
    )
    file_stem <- sprintf("%s_%02d_%s", source_name, i, term_row$module)

    ggsave(
      file.path(output_dir, paste0(file_stem, ".pdf")),
      curve_plot, width = 7.2, height = 5.6, device = cairo_pdf
    )
    ggsave(
      file.path(output_dir, paste0(file_stem, ".png")),
      curve_plot, width = 7.2, height = 5.6, dpi = 300, bg = "white"
    )
    output_files[i] <- file_stem
  }

  output_files
}

make_curve_collection <- function(gsea_result, term_df, module_colors, figure_title) {
  single_plots <- lapply(seq_len(nrow(term_df)), function(i) {
    term_row <- term_df[i, , drop = FALSE]
    make_single_curve(
      gsea_result,
      term_row,
      unname(module_colors[term_row$module])
    )
  })

  patchwork::wrap_plots(plotlist = single_plots, ncol = 2) +
    patchwork::plot_annotation(
      title = figure_title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
}

p_go <- make_curve_collection(
  gsea_objects$gsea_go,
  go_terms,
  go_colors,
  "Representative GO GSEA enrichment curves"
)
p_kegg <- make_curve_collection(
  gsea_objects$gsea_kegg,
  kegg_terms,
  kegg_colors,
  "Representative KEGG GSEA enrichment curves"
)

go_single_files <- save_single_curves(
  gsea_objects$gsea_go,
  go_terms,
  go_colors,
  go_curve_dir,
  "GO"
)
kegg_single_files <- save_single_curves(
  gsea_objects$gsea_kegg,
  kegg_terms,
  kegg_colors,
  kegg_curve_dir,
  "KEGG"
)

ggsave(
  file.path(figure_dir, "GSEA_GO_curves.pdf"),
  p_go, width = 15, height = 14, device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "GSEA_GO_curves.png"),
  p_go, width = 15, height = 14, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "GSEA_KEGG_curves.pdf"),
  p_kegg, width = 15, height = 14, device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "GSEA_KEGG_curves.png"),
  p_kegg, width = 15, height = 14, dpi = 300, bg = "white"
)

cat("Representative GSEA enrichment curves completed.\n")
cat("Individual GO curves:", paste(go_single_files, collapse = ", "), "\n")
cat("Individual KEGG curves:", paste(kegg_single_files, collapse = ", "), "\n")
print(curve_terms)
