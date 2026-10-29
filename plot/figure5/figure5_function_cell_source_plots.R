PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else file.path(PROJECT_ROOT, "3-<U+7ED3><U+679C><U+4E0E><U+5206><U+6790>", "3-<U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4>", "Figure5_ctniche8")
input <- file.path(root, "20-5C_region_corrected_gene_integration")
out <- file.path(root, "22-5C_signature_function_and_cell_source")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(x) read.delim(x, check.names = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, n) write.table(x, file.path(out, n), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
save_plot <- function(p, stem, w, h) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}

dat <- read_tsv(file.path(input, "20_03_original929_region_corrected_integrated_evidence.tsv"))
fun0 <- read_tsv(file.path(input, "20_05_function_module_by_corrected_tier.tsv"))
cell0 <- read_tsv(file.path(input, "20_06_cell_source_by_corrected_tier.tsv"))

tier_labels <- c(
  Tier_A_region_corrected_ultra_strict = "Tier A: LOO strict (472)",
  Tier_B_region_corrected_three_FDR = "Tier B: three-disease FDR (287)",
  Tier_C_region_corrected_direction_stable = "Tier C: direction stable (170)"
)
module_labels <- c(
  Adhesion_migration = "Adhesion / migration",
  ECM_collagen = "ECM / collagen",
  Growth_factor_signaling = "Growth-factor signaling",
  Immune_myeloid = "Immune / myeloid",
  Stress_hypoxia = "Stress / hypoxia",
  Vascular_endothelial = "Vascular / endothelial"
)
module_cols <- c(
  "Adhesion / migration" = "#7A5195", "ECM / collagen" = "#BC5090",
  "Growth-factor signaling" = "#EF5675", "Immune / myeloid" = "#FFA600",
  "Stress / hypoxia" = "#4C78A8", "Vascular / endothelial" = "#59A14F"
)
cell_cols <- c(
  Fibroblast = "#D95F5F", Myeloid = "#4C78A8", Endothelial = "#59A14F",
  Cardiomyocyte = "#B279A2", Other = "#F2A65A", No_marker_assignment = "grey78"
)

# Function-module comparison. A gene may occur in more than one module, so percentages do not sum to 100%.
fun <- fun0 %>%
  filter(function_modules != "Unassigned_by_selected_rules") %>%
  mutate(
    tier = factor(recode(corrected_evidence_tier, !!!tier_labels), levels = unname(tier_labels)),
    module = recode(function_modules, !!!module_labels),
    module = factor(module, levels = rev(unname(module_labels)))
  )
write_tsv(fun, "22_01_function_module_by_evidence_tier.tsv")

p1 <- ggplot(fun, aes(fraction_within_tier * 100, module, colour = tier, group = tier)) +
  geom_line(aes(group = module), colour = "grey88", linewidth = .45) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c(
    "Tier A: LOO strict (472)" = "#2E8B57",
    "Tier B: three-disease FDR (287)" = "#4C78A8",
    "Tier C: direction stable (170)" = "#F2A65A"
  )) +
  labs(
    x = "Genes in each evidence tier (%)", y = NULL, colour = NULL,
    title = "Functional programs retained after region correction",
    subtitle = "Multi-label module assignment; percentages are calculated against all genes in each tier"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")
save_plot(p1, "22_Fig1_function_modules_by_evidence_tier", 8.6, 5.4)

# Tier-A subtype sources, ordered so fibroblasts and myeloid populations stay together.
subtype_order <- c(
  "FN1+ Fib", "HIF3A+ Fib", "THBS4+ Fib",
  "CPM+ MP", "FOLR2+ MP", "Monocyte", "SPP1+ MP",
  "Arterial_endo", "Capillary_endo", "Endocardial_endo", "Lymphatic_endo", "Venous_endo",
  "Normal_vCM", "Stress_vCM", "Adipocyte", "CD4+T_cells", "CD8+T_cells", "Mast", "NK",
  "Neuronal", "Pericyte", "SMC"
)
tierA_cell <- cell0 %>%
  filter(
    corrected_evidence_tier == "Tier_A_region_corrected_ultra_strict",
    primary_subtype_plot != "No_marker_assignment"
  ) %>%
  mutate(primary_subtype_plot = factor(primary_subtype_plot, levels = rev(subtype_order))) %>%
  arrange(cell_source_group, desc(n_genes))
write_tsv(tierA_cell, "22_02_TierA472_primary_subtype_source.tsv")

p2 <- ggplot(tierA_cell, aes(n_genes, primary_subtype_plot, fill = cell_source_group)) +
  geom_col(width = .72) +
  geom_text(aes(label = n_genes), hjust = -.18, size = 3.2) +
  scale_fill_manual(values = cell_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, .15))) +
  labs(
    x = "Tier-A genes assigned as primary subtype markers", y = NULL, fill = "Cell class",
    title = "Single-cell subtype sources of the strict spatial signature",
    subtitle = "Primary marker assignment; 268 genes without a selected subtype marker are not plotted"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
save_plot(p2, "22_Fig2_TierA472_cell_subtype_sources", 8.3, 8.2)

# Gene-level multi-label expansion for the Tier-A function-by-cell-class map.
tierA <- dat %>% filter(corrected_evidence_tier == "Tier_A_region_corrected_ultra_strict")
long <- tierA %>%
  select(gene, cell_source_group, function_modules) %>%
  separate_rows(function_modules, sep = ";") %>%
  filter(
    !is.na(function_modules), function_modules != "", function_modules != "Unassigned_by_selected_rules",
    !is.na(cell_source_group), cell_source_group != "No_marker_assignment"
  ) %>%
  distinct(gene, cell_source_group, function_modules) %>%
  count(cell_source_group, function_modules, name = "n_genes") %>%
  mutate(
    module = recode(function_modules, !!!module_labels),
    module = factor(module, levels = unname(module_labels)),
    cell_source_group = factor(cell_source_group,
      levels = c("Fibroblast", "Myeloid", "Endothelial", "Cardiomyocyte", "Other")
    )
  )
write_tsv(long, "22_03_TierA472_function_by_cell_class.tsv")

p3 <- ggplot(long, aes(module, cell_source_group, fill = n_genes)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = n_genes), size = 4) +
  scale_fill_gradient(low = "#F3F6F8", high = "#B2182B") +
  labs(
    x = NULL, y = NULL, fill = "Genes",
    title = "Functional programs and their principal cellular sources",
    subtitle = "Tier-A genes only; both function and cell-source assignments are required"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
save_plot(p3, "22_Fig3_TierA472_function_by_cell_class", 8.5, 5.0)

# Cell-class composition across evidence tiers; unassigned genes shown explicitly.
cell_tier <- cell0 %>%
  group_by(corrected_evidence_tier, corrected_tier_order, cell_source_group) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  group_by(corrected_evidence_tier) %>%
  mutate(fraction = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  mutate(
    tier = factor(recode(corrected_evidence_tier, !!!tier_labels), levels = unname(tier_labels)),
    cell_source_group = factor(cell_source_group,
      levels = c("Fibroblast", "Myeloid", "Endothelial", "Cardiomyocyte", "Other", "No_marker_assignment")
    )
  )
write_tsv(cell_tier, "22_04_cell_class_composition_by_evidence_tier.tsv")

p4 <- ggplot(cell_tier, aes(tier, fraction * 100, fill = cell_source_group)) +
  geom_col(width = .7) +
  scale_fill_manual(values = cell_cols, labels = c(No_marker_assignment = "No marker assignment")) +
  labs(
    x = NULL, y = "Genes in evidence tier (%)", fill = "Cell class",
    title = "Cell-class attribution across spatial evidence tiers",
    subtitle = "Single primary subtype assignment per gene"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "right")
save_plot(p4, "22_Fig4_cell_class_by_evidence_tier", 8.4, 5.5)

notes <- c(
  "Figure5C Local22: function and cellular-source figures",
  "Primary biological set: the 472 genes in Tier A, which belong to the original 929 and pass the region-corrected edgeR and LOO criteria.",
  "Function modules are multi-label rule-based annotations; one gene may contribute to several modules, so module percentages do not sum to 100%.",
  "Cell-source panels use the single-cell subtype marker table and one primary-subtype assignment per gene.",
  "No_marker_assignment means that a spatial gene was not a selected subtype marker; it does not mean that the gene has no cellular origin.",
  "New edgeR-only genes cannot be assigned fairly from the old929 annotation table alone; they should be remapped with the same single-cell marker rules before being included in these biological panels."
)
writeLines(notes, file.path(out, "22_05_figure_interpretation_notes.txt"), useBytes = TRUE)
message("Completed: ", out)
