PROJECT_ROOT <- Sys.getenv('PROJECT_ROOT', unset = '.')

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', 'Figure5_ctniche8')
in19 <- file.path(root, "19-5C_region_corrected_pseudobulk")
in20 <- file.path(root, "20-5C_region_corrected_gene_integration")
out <- file.path(root, "21-5C_region_corrected_signature_figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(x) read.delim(x, check.names = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
save_plot <- function(p, stem, width, height) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
         dpi = 400, bg = "white", device = grDevices::png, type = "cairo")
}

old <- read_tsv(file.path(in20, "20_03_original929_region_corrected_integrated_evidence.tsv"))
new874 <- read_tsv(file.path(in20, "20_08_de_novo874_shared_annotated.tsv"))
new510 <- read_tsv(file.path(in20, "20_09_de_novo510_ultra_strict_annotated.tsv"))

# Only remove canonical ribosomal structural genes found in the new statistical set.
# RPS6KA1 encodes a kinase and is deliberately retained.
canonical_ribosomal <- c("RPL13", "RPS11", "RPL28", "RPL3")
new870 <- new874 %>% filter(!gene %in% canonical_ribosomal)
strict_filtered <- new510 %>% filter(!gene %in% canonical_ribosomal)

old_genes <- unique(old$gene)
new_genes <- unique(new874$gene)
new870_genes <- unique(new870$gene)
strict_genes <- unique(new510$gene)
overlap_genes <- intersect(old_genes, new_genes)

membership <- tibble(
  gene = sort(unique(c(old_genes, new_genes))),
  voom_limma_929 = gene %in% old_genes,
  edgeR_874 = gene %in% new_genes,
  cross_method_759 = gene %in% overlap_genes,
  edgeR_LOO_strict_510 = gene %in% strict_genes,
  recommended_set = case_when(
    gene %in% strict_genes & gene %in% overlap_genes ~ "Cross-method + LOO strict",
    gene %in% strict_genes ~ "edgeR-new + LOO strict",
    gene %in% overlap_genes ~ "Cross-method stable",
    gene %in% old_genes ~ "voom-limma only",
    TRUE ~ "edgeR only"
  )
)

write_tsv(new870, "21_01_de_novo870_biological_signature.tsv")
write_tsv(strict_filtered, "21_02_de_novo_LOO_strict_technical_filtered.tsv")
write_tsv(membership, "21_03_signature_set_membership.tsv")
write_tsv(tibble(gene = sort(overlap_genes)), "21_04_cross_method_stable759_genes.tsv")
write_tsv(tibble(gene = sort(setdiff(old_genes, new_genes))), "21_05_old929_only170_genes.tsv")
write_tsv(tibble(gene = sort(setdiff(new_genes, old_genes))), "21_06_new874_only115_genes.tsv")

cols <- c(
  "voom-limma only" = "#9E9E9E",
  "Cross-method stable" = "#4C78A8",
  "edgeR only" = "#F2A65A",
  "Cross-method + LOO strict" = "#2E8B57",
  "edgeR-new + LOO strict" = "#C44E52"
)

# Fig1: nested evidence sets, emphasizing analytical roles rather than replacement.
set_counts <- tibble(
  set = factor(c("Original voom-limma", "Region-corrected edgeR", "Cross-method overlap", "LOO strict core"),
               levels = rev(c("Original voom-limma", "Region-corrected edgeR", "Cross-method overlap", "LOO strict core"))),
  n = c(length(old_genes), length(new_genes), length(overlap_genes), length(strict_genes)),
  label = c("929\noriginal discovery", "874\nregion-corrected replication", "759\ncross-method stable", "510\nLOO strict")
)
p1 <- ggplot(set_counts, aes(n, set, fill = set)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = label), hjust = -0.08, size = 4.1, lineheight = .95) +
  scale_fill_manual(values = c("Original voom-limma" = "#8DA0CB", "Region-corrected edgeR" = "#FC8D62",
                               "Cross-method overlap" = "#66C2A5", "LOO strict core" = "#2E8B57"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .28))) +
  labs(x = "Number of genes", y = NULL, title = "ctniche_8 signature evidence sets",
       subtitle = "AMI-IZ, ICM-FZ and HCM; sample-by-niche paired pseudobulk") +
  theme_classic(base_size = 12) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_plot(p1, "21_Fig1_signature_evidence_sets", 8.2, 4.8)

# Fig2: exact old/new partition.
partition <- membership %>% count(recommended_set, name = "n") %>%
  mutate(recommended_set = factor(recommended_set,
    levels = c("voom-limma only", "Cross-method stable", "Cross-method + LOO strict", "edgeR-new + LOO strict", "edgeR only")))
p2 <- ggplot(partition, aes(recommended_set, n, fill = recommended_set)) +
  geom_col(width = .72) + geom_text(aes(label = n), vjust = -.35, fontface = "bold") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = NULL, y = "Number of genes", title = "Cross-model membership and LOO stability",
       subtitle = "The strict core contains both replicated and edgeR-new genes") +
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p2, "21_Fig2_cross_model_membership", 8.4, 5.2)

# Fig3: disease-specific effect sizes of the edgeR shared set.
hm <- new874 %>%
  transmute(gene, AMI_IZ = logFC__AMI_IZ, ICM_FZ = logFC__ICM_FZ, HCM = logFC__HCM,
            strict = gene %in% strict_genes) %>%
  pivot_longer(c(AMI_IZ, ICM_FZ, HCM), names_to = "disease", values_to = "logFC") %>%
  mutate(disease = factor(disease, c("AMI_IZ", "ICM_FZ", "HCM"), c("AMI-IZ", "ICM-FZ", "HCM")),
         strict = factor(strict, c(FALSE, TRUE), c("Other shared genes", "LOO strict core")))
p3 <- ggplot(hm, aes(logFC, colour = disease, fill = disease)) +
  geom_density(alpha = .13, linewidth = .8) +
  geom_vline(xintercept = .1, linetype = 2, colour = "grey45") +
  facet_wrap(~strict, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("AMI-IZ" = "#49c2d9", "ICM-FZ" = "#f59c76", "HCM" = "#fded95")) +
  scale_fill_manual(values = c("AMI-IZ" = "#49c2d9", "ICM-FZ" = "#f59c76", "HCM" = "#fded95")) +
  labs(x = "Region-corrected pseudobulk logFC", y = "Density", colour = NULL, fill = NULL,
       title = "Effect-size distributions across disease states",
       subtitle = "ctniche_8 versus other niches within the matched disease region") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
save_plot(p3, "21_Fig3_disease_logFC_distributions", 8.2, 6.4)

# Fig4: leave-one-sample-out FDR support in each disease.
loo <- new874 %>%
  transmute(gene, strict = gene %in% strict_genes,
            `AMI-IZ` = fraction_LOO_FDR05__AMI_IZ,
            `ICM-FZ` = fraction_LOO_FDR05__ICM_FZ,
            HCM = fraction_LOO_FDR05__HCM) %>%
  pivot_longer(c(`AMI-IZ`, `ICM-FZ`, HCM), names_to = "disease", values_to = "fraction") %>%
  mutate(strict = factor(strict, c(FALSE, TRUE), c("Other shared genes", "LOO strict core")),
         disease = factor(disease, c("AMI-IZ", "ICM-FZ", "HCM")))
p4 <- ggplot(loo, aes(disease, fraction, fill = disease)) +
  geom_violin(scale = "width", trim = TRUE, alpha = .75, colour = "white") +
  geom_boxplot(width = .14, outlier.shape = NA, fill = "white", linewidth = .35) +
  geom_hline(yintercept = .8, linetype = 2, colour = "#B2182B") +
  facet_wrap(~strict, nrow = 1) +
  scale_fill_manual(values = c("AMI-IZ" = "#49c2d9", "ICM-FZ" = "#f59c76", "HCM" = "#fded95"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = NULL, y = "Fraction of LOO runs with FDR < 0.05", title = "Leave-one-sample-out significance stability",
       subtitle = "Dashed line indicates the 80% stability requirement") +
  theme_classic(base_size = 12)
save_plot(p4, "21_Fig4_LOO_FDR_stability", 8.5, 5.2)

# Fig5: concordance between disease effect sizes.
pairdat <- new874 %>% mutate(class = ifelse(gene %in% strict_genes, "LOO strict core", "Other shared genes"))
p5 <- ggplot(pairdat, aes(logFC__AMI_IZ, logFC__ICM_FZ, colour = class, size = -log10(pmax(FDR__HCM, 1e-300)))) +
  geom_point(alpha = .62) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  scale_colour_manual(values = c("Other shared genes" = "grey72", "LOO strict core" = "#2E8B57")) +
  scale_size_continuous(range = c(.5, 3.2), name = expression(-log[10]("HCM FDR"))) +
  labs(x = "AMI-IZ logFC", y = "ICM-FZ logFC", colour = NULL,
       title = "Cross-disease concordance of the region-corrected signature",
       subtitle = "Point size represents independent HCM significance") +
  theme_classic(base_size = 12) + theme(legend.position = "right")
save_plot(p5, "21_Fig5_cross_disease_effect_concordance", 7.3, 6.2)

summary_lines <- c(
  "Figure5C Local21: region-corrected signature figures",
  paste0("Original voom-limma discovery set: ", length(old_genes)),
  paste0("Region-corrected edgeR statistical set: ", length(new_genes)),
  paste0("Refined technical-filtered edgeR set: ", length(new870_genes)),
  paste0("Cross-method overlap: ", length(overlap_genes)),
  paste0("edgeR LOO strict set: ", length(strict_genes)),
  paste0("Technical-filtered LOO strict set: ", nrow(strict_filtered)),
  paste0("Removed canonical ribosomal genes: ", paste(intersect(canonical_ribosomal, new_genes), collapse = ", ")),
  "RPS6KA1 retained: it encodes a kinase rather than a canonical ribosomal structural protein.",
  "The 929, 874, 759 and 510 sets represent discovery, model replication, cross-model stability and LOO robustness, respectively."
)
writeLines(summary_lines, file.path(out, "21_07_results_and_figure_notes.txt"), useBytes = TRUE)
message("Completed: ", out)