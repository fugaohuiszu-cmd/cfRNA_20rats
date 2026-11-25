# 用途：绘制GO与KEGG富集汇总图
# 输入：GSEA_GO.csv和GSEA_KEGG.csv
# 输出：富集汇总图及绘图数据

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_dir <- "F:/DE-analysis/final/cfRNA_project"
result_dir <- file.path(project_dir, "03_results", "02_DE")
table_dir <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")

go_file <- file.path(table_dir, "GSEA_GO.csv")
kegg_file <- file.path(table_dir, "GSEA_KEGG.csv")
required_files <- c(go_file, kegg_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing GSEA result files: ", paste(missing_files, collapse = ", "))
}

fdr_cutoff <- 0.05
top_n_per_module <- 5L

go_modules <- data.frame(
  module = c("Splicing", "Cytoskeleton", "ECM", "Bone", "Calcium"),
  pattern = c(
    "splice|splicing|spliceosome|snrnp|pre-mrna|mrna .*site|rna splicing",
    "cytoskeleton|actin|focal adhesion|intermediate filament|microtubule|cytokinesis|muscle contraction|muscle cells",
    "extracellular matrix|collagen|ecm-receptor|matrix disassembly|matrix organization|fibril|glycosaminoglycan|chondroitin|dermatan",
    "bone|ossification|osteoclast|osteoblast|chondrocyte|skeletal system|skeletal.*morphogenesis|cartilage",
    "calcium|mineral absorption|mineralization|ion transport|store-operated calcium"
  ),
  stringsAsFactors = FALSE
)

kegg_modules <- data.frame(
  module = c("Mineral", "Hormone", "Inflammation", "Immune", "Signaling"),
  pattern = c(
    "mineral absorption|calcium signaling|calcium reabsorption|osteoclast differentiation|parathyroid hormone|aldosterone-regulated sodium reabsorption",
    "thyroid hormone|renin secretion|insulin|growth hormone|estrogen|progesterone|steroid hormone",
    "chemokine|cytokine|nf-kappa|tnf|jak-stat|inflammatory|infection|hepatitis|influenza|chagas|yersinia|cytomegalovirus|herpes|virus|viral",
    "t cell|b cell|th17|natural killer|leukemia|platelet activation|complement|antigen|allograft|autoimmune|hematopoietic",
    "ras signaling|rap1 signaling|pi3k-akt|mapk|sphingolipid signaling|focal adhesion|foxo|notch|phospholipase|cell adhesion molecules|ecm-receptor"
  ),
  stringsAsFactors = FALSE
)

go_df <- read.csv(go_file, check.names = FALSE)
kegg_df <- read.csv(kegg_file, check.names = FALSE)
required_columns <- c(
  "direction", "ID", "Description", "NES", "pvalue", "p.adjust",
  "setSize", "core_enrichment"
)
if (!all(c(required_columns, "ONTOLOGY") %in% colnames(go_df))) {
  stop("GO table does not contain the expected columns.")
}
if (!all(required_columns %in% colnames(kegg_df))) {
  stop("KEGG table does not contain the expected columns.")
}

wrap_text <- function(x, width = 44) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

count_core_genes <- function(x) {
  vapply(
    x,
    function(z) {
      if (is.na(z) || !nzchar(z)) {
        return(0L)
      }
      length(strsplit(z, "/", fixed = TRUE)[[1]])
    },
    integer(1)
  )
}

select_module_terms <- function(gsea_df, module_def, source_name) {
  module_levels <- module_def$module
  significant_df <- gsea_df %>%
    filter(
      is.finite(NES),
      !is.na(p.adjust),
      p.adjust > 0,
      p.adjust < fdr_cutoff
    )

  candidates <- bind_rows(lapply(seq_len(nrow(module_def)), function(i) {
    significant_df %>%
      filter(grepl(module_def$pattern[i], tolower(Description), perl = TRUE)) %>%
      mutate(module = module_def$module[i])
  }))

  if (nrow(candidates) == 0L) {
    stop("No FDR-significant ", source_name, " terms matched the module definitions.")
  }

  candidates <- candidates %>%
    mutate(module_order = match(module, module_levels)) %>%
    arrange(ID, Description, module_order, p.adjust, pvalue, desc(abs(NES))) %>%
    group_by(ID, Description) %>%
    slice_head(n = 1L) %>%
    ungroup()

  missing_modules <- setdiff(module_levels, unique(candidates$module))
  if (length(missing_modules) > 0L) {
    stop(
      source_name,
      " has no FDR-significant terms in modules: ",
      paste(missing_modules, collapse = ", ")
    )
  }

  selected <- candidates %>%
    mutate(module = factor(module, levels = module_levels)) %>%
    arrange(module, p.adjust, pvalue, desc(abs(NES))) %>%
    group_by(module) %>%
    slice_head(n = top_n_per_module) %>%
    mutate(module_rank = row_number()) %>%
    ungroup() %>%
    mutate(
      module = as.character(module),
      source = source_name,
      comparison = "OP vs N",
      plot_logFDR = -log10(p.adjust),
      core_gene_count = count_core_genes(core_enrichment),
      leading_edge_percent = 100 * core_gene_count / setSize,
      term_key = paste(module, ID, sep = "__"),
      term_label = wrap_text(Description)
    )

  if (any(selected$setSize <= 0) ||
      any(selected$core_gene_count > selected$setSize) ||
      any(!is.finite(selected$leading_edge_percent))) {
    stop("Invalid leading-edge gene ratio in selected ", source_name, " terms.")
  }

  selected
}

go_plot_df <- select_module_terms(go_df, go_modules, "GO")
kegg_plot_df <- select_module_terms(kegg_df, kegg_modules, "KEGG")

prepare_plot_factors <- function(plot_df, module_levels) {
  term_levels <- plot_df %>%
    mutate(module = factor(module, levels = module_levels)) %>%
    arrange(module, desc(module_rank)) %>%
    pull(term_key)

  plot_df %>%
    mutate(
      module = factor(module, levels = module_levels),
      term_key = factor(term_key, levels = term_levels)
    )
}

go_plot_df <- prepare_plot_factors(go_plot_df, go_modules$module)
kegg_plot_df <- prepare_plot_factors(kegg_plot_df, kegg_modules$module)

write.csv(go_plot_df, file.path(table_dir, "GSEA_GO_plot.csv"), row.names = FALSE)
write.csv(kegg_plot_df, file.path(table_dir, "GSEA_KEGG_plot.csv"), row.names = FALSE)

all_nes <- c(go_plot_df$NES, kegg_plot_df$NES)
nes_limit <- max(2, ceiling(max(abs(all_nes), na.rm = TRUE) * 10) / 10)
all_edge_percent <- c(go_plot_df$leading_edge_percent, kegg_plot_df$leading_edge_percent)
edge_limit <- ceiling(max(all_edge_percent, na.rm = TRUE) / 5) * 5

make_module_plot <- function(plot_df, plot_title) {
  term_labels <- setNames(plot_df$term_label, as.character(plot_df$term_key))

  ggplot(plot_df, aes(x = NES, y = term_key)) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      color = "#737373",
      linewidth = 0.45
    ) +
    geom_point(
      aes(size = plot_logFDR, fill = leading_edge_percent),
      shape = 21,
      color = "black",
      stroke = 0.35
    ) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = term_labels) +
    scale_size_continuous(
      name = expression(-log[10](adjusted~P~value)),
      range = c(3.5, 9),
      breaks = scales::breaks_pretty(n = 4)
    ) +
    scale_fill_gradientn(
      colors = c("#DCEAF4", "#52A6A4", "#E9B44C", "#C94C3C"),
      limits = c(0, edge_limit),
      oob = scales::squish,
      name = "Leading-edge\ngenes (%)"
    ) +
    scale_x_continuous(
      limits = c(-nes_limit, nes_limit),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    labs(
      title = plot_title,
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.6),
      panel.grid.minor = element_blank(),
      strip.background.y = element_rect(fill = "white", color = "black", linewidth = 0.9),
      strip.text.y.right = element_text(
        angle = -90,
        face = "bold",
        size = 10,
        margin = margin(4, 4, 4, 4)
      ),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 10.5, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.spacing.y = grid::unit(0.08, "lines")
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(
        order = 2,
        override.aes = list(fill = "grey70", color = "black")
      )
    )
}

p_go <- make_module_plot(go_plot_df, "GO GSEA: OP vs N")
p_kegg <- make_module_plot(kegg_plot_df, "KEGG GSEA: OP vs N")

ggsave(file.path(figure_dir, "GSEA_GO.pdf"), p_go, width = 10, height = 10.2, device = cairo_pdf)
ggsave(
  file.path(figure_dir, "GSEA_GO.png"), p_go,
  width = 10, height = 10.2, dpi = 300, bg = "white"
)
ggsave(file.path(figure_dir, "GSEA_KEGG.pdf"), p_kegg, width = 10, height = 9.8, device = cairo_pdf)
ggsave(
  file.path(figure_dir, "GSEA_KEGG.png"), p_kegg,
  width = 10, height = 9.8, dpi = 300, bg = "white"
)

cat("Module-based GSEA plots completed.\n")
cat("GO plotted terms by module:\n")
print(go_plot_df %>% count(module, name = "n_terms"))
cat("KEGG plotted terms by module:\n")
print(kegg_plot_df %>% count(module, name = "n_terms"))
cat("Criterion: adjusted P value < 0.05; ranked by adjusted P value.\n")
