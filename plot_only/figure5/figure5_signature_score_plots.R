PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else file.path(PROJECT_ROOT, "results", "figure5", "strict81_bulk")
output_dir <- if (length(args) >= 2) args[2] else file.path(PROJECT_ROOT, "outputs", "figure5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) read.delim(file.path(input_dir, name), check.names = FALSE)
scores <- read_tsv("36_03_all_sample_scores.tsv")
qc <- read_tsv("36_05_dataset_QC_and_gene_coverage.tsv")
heatmap_data <- read_tsv("36_10_top40_heatmap_long.tsv")

group_colors <- c(Control = "#d0e2c0", AMI = "#49c2d9", ICM = "#f59c76", HCM = "#fded95")
scores <- scores %>% mutate(display_group = if_else(status == "Control", "Control", disease))
score_long <- scores %>%
  pivot_longer(c(rank_score, z_score), names_to = "method", values_to = "score") %>%
  mutate(method = recode(method, rank_score = "Within-sample rank score", z_score = "Gene-wise z score"))

p1 <- ggplot(score_long, aes(display_group, score, fill = display_group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "grey35") +
  geom_point(position = position_jitter(width = 0.10, seed = 42), size = 1.3, alpha = 0.75) +
  facet_grid(method ~ dataset, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Signature score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(output_dir, "figure5_signature_scores_rank_and_z.pdf"), p1, width = 11, height = 8, device = cairo_pdf)

count_data <- bind_rows(
  transmute(qc, dataset, class = "Significant up", n = significant_up),
  transmute(qc, dataset, class = "Significant down", n = significant_down)
)
p2 <- ggplot(count_data, aes(dataset, n, fill = class)) +
  geom_col(position = "dodge", width = 0.68) +
  geom_text(aes(label = n), position = position_dodge(0.68), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Significant up" = "#B5403A", "Significant down" = "#3568A8")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of signature genes", fill = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(output_dir, "figure5_signature_de_counts.pdf"), p2, width = 8, height = 6.5, device = cairo_pdf)

limit <- max(abs(quantile(heatmap_data$logFC, c(0.02, 0.98), na.rm = TRUE)))
p3 <- ggplot(heatmap_data, aes(dataset, gene, fill = pmax(pmin(logFC, limit), -limit))) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(
    data = filter(heatmap_data, class %in% c("Significant_Up", "Significant_Down"), !is.na(logFC)),
    shape = 21, size = 1.5, stroke = 0.3, fill = "white", color = "black"
  ) +
  scale_fill_gradient2(low = "#3568A8", mid = "#F5F3EC", high = "#B5403A", midpoint = 0, na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "logFC") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "italic"))
ggsave(file.path(output_dir, "figure5_top40_direction_heatmap.pdf"), p3, width = 9, height = 10, device = cairo_pdf)
