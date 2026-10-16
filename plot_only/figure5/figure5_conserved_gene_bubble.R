suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

infile <- "5B_conserved_ctniche8_signature.tsv"
stopifnot(file.exists(infile))
x <- read.delim(infile, check.names = FALSE)

# <U+53EA><U+5C55><U+793A><U+6700><U+7A33><U+5B9A><U+3001><U+6700><U+6709><U+533A><U+5206><U+5EA6><U+7684><U+6838><U+5FC3><U+57FA><U+56E0><U+FF0C><U+907F><U+514D><U+628A>1788<U+4E2A><U+57FA><U+56E0><U+5168><U+90E8><U+585E><U+5165><U+4E3B><U+56FE><U+3002>
# <U+6392><U+5E8F><U+4F18><U+5148><U+8003><U+8651><U+4E09><U+79CD><U+75BE><U+75C5><U+4E2D><U+6700><U+5F31><U+7684><U+6548><U+5E94><U+91CF><U+FF0C><U+800C><U+4E0D><U+662F><U+5E73><U+5747><U+6548><U+5E94><U+91CF><U+FF0C><U+907F><U+514D><U+7531>HCM<U+5355><U+72EC><U+4E3B><U+5BFC><U+3002>
N_SHOW <- 30
core <- x %>%
  filter(
    logFC_AMI_IZ > 0,
    logFC_ICM_FZ > 0,
    logFC_HCM > 0,
    patient_positive_fraction >= 0.80,
    !grepl("^(MT-|RPL|RPS|HBA|HBB)", gene)
  ) %>%
  mutate(
    min_logFC_plot = pmin(logFC_AMI_IZ, logFC_ICM_FZ, logFC_HCM),
    max_FDR_plot = pmax(FDR_AMI_IZ, FDR_ICM_FZ, FDR_HCM)
  ) %>%
  arrange(desc(min_logFC_plot), max_FDR_plot) %>%
  slice_head(n = N_SHOW)

if (nrow(core) == 0) stop("<U+6CA1><U+6709><U+6EE1><U+8DB3><U+4F5C><U+56FE><U+6761><U+4EF6><U+7684><U+4FDD><U+5B88><U+57FA><U+56E0>")

long <- core %>%
  select(
    gene,
    logFC_AMI_IZ, logFC_ICM_FZ, logFC_HCM,
    FDR_AMI_IZ, FDR_ICM_FZ, FDR_HCM
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("measure", "Disease"),
    names_pattern = "(logFC|FDR)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(
    Disease = factor(Disease,
      levels = c("AMI_IZ", "ICM_FZ", "HCM"),
      labels = c("AMI-IZ", "ICM-FZ", "HCM")
    ),
    minus_log10_FDR = pmin(-log10(pmax(FDR, 1e-300)), 20)
  )

gene_order <- core %>%
  arrange(min_logFC_plot) %>%
  pull(gene)
long$gene <- factor(long$gene, levels = gene_order)
core$gene <- factor(core$gene, levels = gene_order)

p1 <- ggplot(long, aes(Disease, gene)) +
  geom_point(aes(size = minus_log10_FDR, colour = logFC)) +
  scale_colour_gradient2(
    low = "#4575B4", mid = "white", high = "#C51B3A",
    midpoint = 0, name = "log2FC"
  ) +
  scale_size_continuous(range = c(1.5, 7), name = expression(-log[10](FDR))) +
  labs(
    x = NULL, y = NULL,
    title = "Conserved ctniche_8-associated genes across disease states"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

p2 <- ggplot(core, aes(patient_positive_fraction, gene)) +
  geom_segment(aes(
    x = 0.5, xend = patient_positive_fraction,
    y = gene, yend = gene
  ), colour = "grey75") +
  geom_point(aes(colour = min_logFC_plot), size = 3) +
  scale_x_continuous(
    limits = c(0.5, 1.02), breaks = c(0.5, 0.75, 1.0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradient(
    low = "#FDAE61", high = "#A50026",
    name = "Minimum\nlog2FC"
  ) +
  labs(x = "Patient pairs with positive change", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p <- p1 + p2 + plot_layout(widths = c(2.4, 1.2), guides = "collect")

ggsave("5B_ctniche8_conserved_genes_bubble.pdf", p, width = 10.5, height = 9)
ggsave("5B_ctniche8_conserved_genes_bubble.png", p,
  width = 10.5, height = 9,
  dpi = 300, bg = "white"
)

write.table(
  core %>% arrange(desc(min_logFC_plot)),
  "5B_ctniche8_conserved_genes_plotted.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Plotted ", nrow(core), " conserved genes.")
