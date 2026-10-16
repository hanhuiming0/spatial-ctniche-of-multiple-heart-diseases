PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(
    PROJECT_ROOT, "3-<U+7ED3><U+679C><U+4E0E><U+5206><U+6790>", "3-<U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4>", "2-ctniche", "5-<U+529F><U+80FD>", "enrichG_up1.RDS"
  )
}
output_file <- if (length(args) >= 2) args[2] else "figure2_ctniche_go_bp.pdf"

enrichment <- readRDS(input_rds)
plot_data <- as.data.frame(enrichment) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    GeneRatioValue = vapply(strsplit(as.character(GeneRatio), "/", fixed = TRUE), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    }, numeric(1)),
    Cluster = factor(Cluster, levels = paste0("ctniche_", 1:12)),
    Description = factor(Description, levels = rev(unique(Description)))
  )

p <- ggplot(plot_data, aes(Cluster, Description)) +
  geom_point(aes(size = GeneRatioValue, color = -log10(p.adjust))) +
  scale_color_gradient(low = "#78B7C5", high = "#B5403A") +
  labs(x = NULL, y = NULL, size = "Gene ratio", color = "-log10(FDR)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())

ggsave(output_file, p, width = 7.4, height = 6, device = cairo_pdf)
