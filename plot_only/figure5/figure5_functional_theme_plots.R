PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else file.path(PROJECT_ROOT, "3-<U+7ED3><U+679C><U+4E0E><U+5206><U+6790>", "3-<U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4>", "Figure5_ctniche8")
recolor_only <- length(args) >= 2 && identical(args[2], "recolor_only")
out <- file.path(root, "27-5C_voom_limma929_function_and_cell_class")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(x) read.delim(x, check.names = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, n) write.table(x, file.path(out, n), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
save_plot <- function(p, s, w, h) {
  if (recolor_only && !grepl("recolored$", s)) {
    return(invisible(NULL))
  }
  ggsave(file.path(out, paste0(s, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}

enr <- file.path(root, "8-5C_strict929_pathway_enrichment")
go <- read_tsv(file.path(enr, "5C_54_strict929_GO_BP_simplified.tsv"))
react <- read_tsv(file.path(enr, "5C_55_strict929_Reactome_all.tsv"))
mods <- read_tsv(file.path(root, "9-5C_core_signature_modules", "5C_65_gene_function_module_long.tsv"))
cells <- read_tsv(file.path(root, "26-5C_voom_limma929_singlecell_marker_distribution", "26_02_shared929_all_subtype_marker_matches.tsv"))

theme_fun <- function(x) {
  z <- tolower(x)
  case_when(
    str_detect(z, "extracellular|collagen|matrix|elastic|proteoglycan") ~ "ECM / collagen remodeling",
    str_detect(z, "immune|inflamm|leukocyte|complement|interleukin|neutrophil|bacteri|cytokine") ~ "Immune / inflammatory response",
    str_detect(z, "chemotaxis|taxis|migration|locomotion|adhesion|integrin") ~ "Migration / adhesion",
    str_detect(z, "vascular|vessel|angiogen|endothelial") ~ "Vascular remodeling",
    str_detect(z, "development|morphogenesis|growth|wound|regeneration|ossification|bone|skeletal|igf") ~ "Growth / tissue repair",
    str_detect(z, "stress|hypoxi|oxidative|apopto") ~ "Stress response", TRUE ~ "Other"
  )
}
theme_cols <- c("ECM / collagen remodeling" = "#B2182B", "Immune / inflammatory response" = "#E08214", "Migration / adhesion" = "#7B3294", "Vascular remodeling" = "#008837", "Growth / tissue repair" = "#2166AC", "Stress response" = "#4D4D4D", "Other" = "grey70")

# GO BP: representative non-redundant terms grouped into biological themes.
gop <- go %>%
  mutate(theme = theme_fun(Description)) %>%
  filter(theme != "Other", p.adjust < .05) %>%
  group_by(theme) %>%
  arrange(p.adjust, desc(FoldEnrichment), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label = str_wrap(Description, 34), label = factor(label, levels = rev(unique(label))))
write_tsv(gop, "27_01_GO_BP_representative_by_theme.tsv")
p1 <- ggplot(gop, aes(FoldEnrichment, label, size = Count, colour = -log10(pmax(p.adjust, 1e-300)))) +
  geom_point(alpha = .9) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = "#74ADD1", high = "#B2182B") +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(x = "Fold enrichment", y = NULL, size = "Genes", colour = expression(-log[10](FDR)), title = "Functional landscape of the conserved 929-gene ctniche_8 program", subtitle = "GO Biological Process; representative terms grouped into biological themes") +
  theme_bw(base_size = 10) +
  theme(strip.text.y = element_text(angle = 0, face = "bold"), panel.grid.minor = element_blank())
save_plot(p1, "27_Fig1_voom929_GO_BP_functional_themes", 10, 10)
p1_recolored <- p1 + scale_colour_gradient(low = "#F4D9A6", high = "#8E1B3A")
save_plot(p1_recolored, "27_Fig1_voom929_GO_BP_functional_themes_recolored", 10, 10)

# Reactome: core pathways, keeping major ECM and immune/remodeling programs visible.
rp <- react %>%
  mutate(theme = theme_fun(Description)) %>%
  filter(p.adjust < .05) %>%
  group_by(theme) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(theme, p.adjust) %>%
  mutate(label = str_wrap(Description, 38), label = factor(label, levels = rev(unique(label))))
write_tsv(rp, "27_02_Reactome_representative_by_theme.tsv")
p2 <- ggplot(rp, aes(FoldEnrichment, label, size = Count, colour = -log10(pmax(p.adjust, 1e-300)))) +
  geom_point(alpha = .9) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = "#74ADD1", high = "#B2182B") +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(x = "Fold enrichment", y = NULL, size = "Genes", colour = expression(-log[10](FDR)), title = "Core pathways of the conserved 929-gene ctniche_8 program", subtitle = "Reactome enrichment; representative pathways by biological theme") +
  theme_bw(base_size = 10) +
  theme(strip.text.y = element_text(angle = 0, face = "bold"), panel.grid.minor = element_blank())
save_plot(p2, "27_Fig2_voom929_Reactome_functional_themes", 10, 9)
p2_recolored <- p2 + scale_colour_gradient(low = "#F4D9A6", high = "#8E1B3A")
save_plot(p2_recolored, "27_Fig2_voom929_Reactome_functional_themes_recolored", 10, 9)

# All-marker overlap: one gene may map to multiple modules and multiple cell classes.
module_names <- c(Adhesion_migration = "Adhesion / migration", ECM_collagen = "ECM / collagen", Growth_factor_signaling = "Growth-factor signaling", Immune_myeloid = "Immune / myeloid", Stress_hypoxia = "Stress / hypoxia", Vascular_endothelial = "Vascular / endothelial")
mc <- mods %>%
  filter(module %in% names(module_names)) %>%
  select(gene, module) %>%
  distinct() %>%
  inner_join(cells %>% select(gene, cell_class) %>% distinct(), by = "gene", relationship = "many-to-many") %>%
  distinct(gene, module, cell_class)
heat <- mc %>%
  count(module, cell_class, name = "n_genes") %>%
  group_by(module) %>%
  mutate(module_total = n_distinct(mc$gene[mc$module == first(module)]), fraction_of_module = n_genes / module_total) %>%
  ungroup() %>%
  mutate(module = recode(module, !!!module_names), module = factor(module, levels = unname(module_names)), cell_class = factor(cell_class, levels = c("Fibroblast", "Myeloid", "Endothelial", "Cardiomyocyte", "Lymphoid / mast", "Other")))
write_tsv(mc, "27_03_gene_function_module_cell_class_long.tsv")
write_tsv(heat, "27_04_function_module_by_cell_class_summary.tsv")
p3 <- ggplot(heat, aes(module, cell_class, fill = n_genes)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = n_genes), size = 4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B") +
  labs(x = NULL, y = NULL, fill = "Genes", title = "Cell-class distribution of conserved ctniche_8 functional modules", subtitle = "All significant subtype-marker overlaps; genes may contribute to more than one cell class") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p3, "27_Fig3_voom929_function_module_by_cell_class_counts", 9, 5.4)

# Within-module percentages improve comparison across modules of different sizes.
p4 <- ggplot(heat, aes(module, cell_class, fill = fraction_of_module * 100)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", fraction_of_module * 100)), size = 3.8) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC") +
  labs(x = NULL, y = NULL, fill = "Module genes (%)", title = "Relative cell-class coverage of conserved ctniche_8 functional modules", subtitle = "Denominator: 929 genes assigned to each functional module; multi-cell-class marker membership is retained") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p4, "27_Fig4_voom929_function_module_cell_class_fraction", 9, 5.4)

# Functional modules within the major fibroblast and myeloid subtypes.
sub <- mods %>%
  filter(module %in% names(module_names)) %>%
  select(gene, module) %>%
  distinct() %>%
  inner_join(cells %>% filter(cell_class %in% c("Fibroblast", "Myeloid")) %>% select(gene, subtype, cell_class) %>% distinct(), by = "gene", relationship = "many-to-many") %>%
  distinct(gene, module, subtype, cell_class) %>%
  count(cell_class, subtype, module, name = "n_genes") %>%
  mutate(module = recode(module, !!!module_names))
write_tsv(sub, "27_05_fibroblast_myeloid_subtype_function_summary.tsv")
p5 <- ggplot(sub, aes(module, subtype, fill = n_genes)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = n_genes), size = 3.5) +
  facet_grid(cell_class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#FFF7EC", high = "#B2182B") +
  labs(x = NULL, y = NULL, fill = "Genes", title = "Functional programs across fibroblast and myeloid subtypes", subtitle = "Based on all marker overlaps rather than exclusive primary-subtype assignment") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0, face = "bold"))
save_plot(p5, "27_Fig5_voom929_fibroblast_myeloid_function_heatmap", 9, 6.2)

writeLines(c("voom-limma 929 functional interpretation figures", "Fig1 and Fig2 summarize complete GO BP and Reactome enrichment.", "Fig3-Fig5 connect pre-defined enriched function modules to cell classes using all significant subtype-marker overlaps.", "A gene may map to multiple functions, subtypes and cell classes; therefore counts and percentages are not mutually exclusive.", "Primary subtype assignment is not used in the main function-to-cell-class panels.", "KEGG enrichment was empty in the existing analysis and was not plotted."), file.path(out, "27_06_results_and_figure_notes.txt"))
message("Completed: ", out)
