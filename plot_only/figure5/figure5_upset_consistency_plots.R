PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else file.path(PROJECT_ROOT, "3-<U+7ED3><U+679C><U+4E0E><U+5206><U+6790>", "3-<U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4>", "Figure5_ctniche8")
in19 <- file.path(root, "19-5C_region_corrected_pseudobulk")
in20 <- file.path(root, "20-5C_region_corrected_gene_integration")
out <- file.path(root, "23-5C_cross_disease_DE_consistency")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(x) read.delim(x, check.names = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, n) write.table(x, file.path(out, n), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
save_plot <- function(p, stem, w, h) ggsave(file.path(out, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)

all <- read_tsv(file.path(in19, "C18_10_de_novo_all_gene_cross_disease_evidence.tsv")) %>%
  filter(tested_all_three %in% TRUE)
strict <- read_tsv(file.path(in20, "20_09_de_novo510_ultra_strict_annotated.tsv"))
old <- read_tsv(file.path(in20, "20_03_original929_region_corrected_integrated_evidence.tsv"))

make_membership <- function(direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (direction == "up") {
    a <- all$logFC__AMI_IZ > .1 & all$FDR__AMI_IZ < .05
    i <- all$logFC__ICM_FZ > .1 & all$FDR__ICM_FZ < .05
    h <- all$logFC__HCM > .1 & all$FDR__HCM < .05
  } else {
    a <- all$logFC__AMI_IZ < -.1 & all$FDR__AMI_IZ < .05
    i <- all$logFC__ICM_FZ < -.1 & all$FDR__ICM_FZ < .05
    h <- all$logFC__HCM < -.1 & all$FDR__HCM < .05
  }
  tibble(gene = all$gene, AMI_IZ = a, ICM_FZ = i, HCM = h) %>%
    mutate(combination = case_when(
      AMI_IZ & ICM_FZ & HCM ~ "AMI-IZ & ICM-FZ & HCM",
      AMI_IZ & ICM_FZ ~ "AMI-IZ & ICM-FZ",
      AMI_IZ & HCM ~ "AMI-IZ & HCM",
      ICM_FZ & HCM ~ "ICM-FZ & HCM",
      AMI_IZ ~ "AMI-IZ only", ICM_FZ ~ "ICM-FZ only", HCM ~ "HCM only",
      TRUE ~ "None"
    ))
}

up <- make_membership("up")
down <- make_membership("down")
write_tsv(up %>% filter(combination != "None"), "23_01_upregulated_gene_set_membership.tsv")
write_tsv(down %>% filter(combination != "None"), "23_02_downregulated_gene_set_membership.tsv")

disease_cols <- c("AMI-IZ" = "#49c2d9", "ICM-FZ" = "#f59c76", "HCM" = "#fded95")
comb_cols <- c(
  "AMI-IZ only" = "#49c2d9", "ICM-FZ only" = "#f59c76", "HCM only" = "#fded95",
  "AMI-IZ & ICM-FZ" = "#7DB7B2", "AMI-IZ & HCM" = "#74B6A2",
  "ICM-FZ & HCM" = "#D5A96E", "AMI-IZ & ICM-FZ & HCM" = "#B2182B"
)

make_upset <- function(mem, title, subtitle, stem) {
  cnt <- mem %>%
    filter(combination != "None") %>%
    count(combination, name = "n") %>%
    arrange(desc(n))
  cnt$combination <- factor(cnt$combination, levels = cnt$combination)
  keys <- as.character(cnt$combination)
  mat <- tidyr::expand_grid(combination = keys, disease = c("AMI-IZ", "ICM-FZ", "HCM")) %>%
    mutate(
      active = case_when(
        disease == "AMI-IZ" ~ grepl("AMI-IZ", combination),
        disease == "ICM-FZ" ~ grepl("ICM-FZ", combination),
        TRUE ~ grepl("HCM", combination)
      ),
      combination = factor(combination, levels = keys),
      disease = factor(disease, levels = rev(c("AMI-IZ", "ICM-FZ", "HCM")))
    )
  setn <- tibble(
    disease = c("AMI-IZ", "ICM-FZ", "HCM"),
    n = c(sum(mem$AMI_IZ), sum(mem$ICM_FZ), sum(mem$HCM))
  ) %>%
    mutate(disease = factor(disease, levels = c("AMI-IZ", "ICM-FZ", "HCM")))
  pbar <- ggplot(cnt, aes(combination, n, fill = combination)) +
    geom_col(width = .68) +
    geom_text(aes(label = n), vjust = -.35, fontface = "bold", size = 3.7) +
    scale_fill_manual(values = comb_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(x = NULL, y = "Intersection size", title = title, subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  pmat <- ggplot(mat, aes(combination, disease)) +
    geom_line(
      data = mat %>% filter(active) %>% group_by(combination), aes(group = combination),
      colour = "grey35", linewidth = .75
    ) +
    geom_point(aes(fill = active), shape = 21, size = 3.7, colour = "grey35", stroke = .25) +
    scale_fill_manual(values = c(`FALSE` = "grey88", `TRUE` = "grey20"), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    )
  pset <- ggplot(setn, aes(n, disease, fill = disease)) +
    geom_col(width = .58) +
    geom_text(aes(x = n / 2, label = n), colour = "grey15", fontface = "bold", size = 3.4) +
    scale_x_reverse(expand = expansion(mult = c(.04, .2))) +
    scale_fill_manual(values = disease_cols, guide = "none") +
    labs(x = "Set size", y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  full <- (plot_spacer() + pbar + plot_layout(widths = c(1.45, 5))) /
    (pset + pmat + plot_layout(widths = c(1.45, 5))) + plot_layout(heights = c(3.4, 1.35))
  save_plot(full, stem, 9.4, 6.8)
  cnt
}

up_cnt <- make_upset(
  up,
  "Conserved ctniche_8 upregulated genes across disease states",
  "Paired sample-by-niche edgeR pseudobulk; common tested universe = 13,580 genes\nlogFC > 0.1 and FDR < 0.05",
  "23_Fig1_upregulated_DEG_UpSet"
)
down_cnt <- make_upset(
  down,
  "Conserved ctniche_8 downregulated genes across disease states",
  "Supplementary analysis; common tested universe = 13,580 genes\nlogFC < -0.1 and FDR < 0.05",
  "23_FigS1_downregulated_DEG_UpSet"
)
write_tsv(up_cnt, "23_03_upregulated_intersection_counts.tsv")
write_tsv(down_cnt, "23_04_downregulated_intersection_counts.tsv")

# Evidence decomposition of the 874 shared-upregulated genes.
shared <- up %>%
  filter(AMI_IZ, ICM_FZ, HCM) %>%
  pull(gene)
old_genes <- unique(old$gene)
strict_genes <- unique(strict$gene)
evidence <- tibble(
  comparison = factor(c("Cross-model support", "Cross-model support", "LOO robustness", "LOO robustness"),
    levels = c("Cross-model support", "LOO robustness")
  ),
  class = c("Also in original voom-limma set", "edgeR-only", "LOO strict", "Not LOO-FDR80 strict"),
  n = c(
    sum(shared %in% old_genes), sum(!shared %in% old_genes),
    sum(shared %in% strict_genes), sum(!shared %in% strict_genes)
  )
)
write_tsv(evidence, "23_05_shared874_evidence_decomposition.tsv")
p2 <- ggplot(evidence, aes(comparison, n, fill = class)) +
  geom_col(width = .62) +
  geom_text(aes(label = n), position = position_stack(vjust = .5), size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "Also in original voom-limma set" = "#4C78A8", "edgeR-only" = "#F2A65A",
    "LOO strict" = "#2E8B57", "Not LOO-FDR80 strict" = "grey78"
  )) +
  labs(
    x = NULL, y = "Genes within the shared 874", fill = NULL,
    title = "Independent support for the 874 shared-upregulated genes"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
save_plot(p2, "23_Fig2_shared874_evidence_decomposition", 7.4, 5.3)

# Representative strict genes: prioritize conservative rank and show raw disease-specific logFC.
top <- old %>%
  filter(gene %in% strict_genes) %>%
  arrange(conservative_rank) %>%
  slice_head(n = 30) %>%
  select(gene, conservative_rank, logFC__AMI_IZ, logFC__ICM_FZ, logFC__HCM)
top$gene <- factor(top$gene, levels = rev(top$gene))
heat <- top %>%
  pivot_longer(starts_with("logFC__"), names_to = "disease", values_to = "logFC") %>%
  mutate(
    disease = recode(disease, logFC__AMI_IZ = "AMI-IZ", logFC__ICM_FZ = "ICM-FZ", logFC__HCM = "HCM"),
    disease = factor(disease, levels = c("AMI-IZ", "ICM-FZ", "HCM"))
  )
write_tsv(heat, "23_06_top30_LOO_strict_gene_logFC.tsv")
p3 <- ggplot(heat, aes(disease, gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = .65) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0.8) +
  labs(
    x = NULL, y = NULL, fill = "logFC",
    title = "Representative LOO-robust conserved ctniche_8 genes",
    subtitle = "Top 30 genes ranked by conservative cross-disease evidence"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))
save_plot(p3, "23_Fig3_top30_LOO_strict_logFC_heatmap", 6.2, 8.2)

summary <- c(
  "Cross-disease ctniche_8 DEG consistency figures",
  paste0("Common tested universe: ", nrow(all)),
  paste0("Upregulated sets: AMI-IZ=", sum(up$AMI_IZ), ", ICM-FZ=", sum(up$ICM_FZ), ", HCM=", sum(up$HCM)),
  paste0("Shared upregulated genes: ", length(shared)),
  paste0("Shared upregulated and original voom-limma supported: ", sum(shared %in% old_genes)),
  paste0("Shared upregulated and LOO strict: ", sum(shared %in% strict_genes)),
  paste0("Shared downregulated genes: ", sum(down$AMI_IZ & down$ICM_FZ & down$HCM)),
  "Main figure recommendation: Fig1 UpSet plus Fig2 evidence decomposition; Fig3 is a representative-gene panel.",
  "Upregulated and downregulated genes are displayed separately."
)
writeLines(summary, file.path(out, "23_07_results_summary.txt"), useBytes = TRUE)
message("Completed: ", out)
