options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(scales)
})

# MISTy <U+4E0D><U+540C><U+7A7A><U+95F4><U+5C3A><U+5EA6><U+7684><U+7EC6><U+80DE><U+4E9A><U+578B><U+7A7A><U+95F4><U+4F9D><U+8D56><U+7ED3><U+679C><U+6C47><U+603B>
# <U+8BF4><U+660E><U+FF1A>importance <U+8868><U+793A> predictor <U+5BF9> target <U+7684><U+9884><U+6D4B><U+8D21><U+732E><U+FF0C><U+4E0D><U+4EE3><U+8868><U+76F8><U+5173><U+65B9><U+5411><U+6216><U+7269><U+7406><U+63A5><U+89E6><U+3002>
# para.2/para.5 <U+4E3A><U+9AD8><U+65AF><U+52A0><U+6743><U+7A7A><U+95F4><U+89C6><U+56FE><U+FF0C><U+4E8C><U+8005><U+4E0D><U+662F><U+4E92><U+65A5><U+7684><U+8DDD><U+79BB><U+73AF><U+3002>

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) stop("<U+8BF7><U+4F7F><U+7528> Rscript <U+8FD0><U+884C><U+672C><U+811A><U+672C><U+FF0C><U+4EE5><U+4FBF><U+81EA><U+52A8><U+5B9A><U+4F4D><U+8F93><U+5165><U+76EE><U+5F55><U+3002>")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- dirname(script_file)
input_dir <- file.path(root, "SubType_MISTy2")
project_root <- normalizePath(file.path(root, "../../../../.."), winslash = "/", mustWork = TRUE)
meta_file <- file.path(project_root, "2-<U+6570><U+636E>", "Sample_Info_2025.8.20.txt")
out_dir <- file.path(root, "SubType_MISTy2_<U+4E0D><U+540C><U+8DDD><U+79BB><U+6C47><U+603B><U+56FE>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(input_dir), file.exists(meta_file))

group_levels <- c(
  "Normal", "AMI-IZ", "AMI-BZ", "AMI-RZ",
  "ICM-FZ", "ICM-IZ", "ICM-BZ", "ICM-RZ", "HCM"
)
main_groups <- c("Normal", "AMI-IZ", "ICM-FZ", "HCM")
view_levels <- c("intra", "para.2", "para.5")
view_labels <- c(intra = "Same spot", para.2 = "l=2", para.5 = "l=5")
group_cols <- c(
  "Normal" = "#6B7280", "AMI-IZ" = "#D55E00",
  "ICM-FZ" = "#0072B2", "HCM" = "#8E5AA9"
)

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", colour = "grey65"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"), legend.position = "right"
    )
}

meta <- read.delim(meta_file, check.names = FALSE) %>%
  filter(Dataset_Type == "Visium") %>%
  mutate(Group = case_when(
    Sample_Type == "Normal" ~ "Normal",
    Sample_Type == "HCM" ~ "HCM",
    TRUE ~ paste(Sample_Type, Disease_State, sep = "-")
  )) %>%
  select(Sample_ID, Patient_ID, Sample_Type, Disease_State, Group)

sample_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
samples <- basename(sample_dirs)
sample_map <- tibble(Sample_ID = samples, sample_dir = sample_dirs) %>%
  left_join(meta, by = "Sample_ID")
if (anyNA(sample_map$Group)) stop("<U+4EE5><U+4E0B><U+6837><U+672C><U+672A><U+5339><U+914D><U+6837><U+672C><U+4FE1><U+606F><U+FF1A>", paste(sample_map$Sample_ID[is.na(sample_map$Group)], collapse = ", "))
sample_map$Group <- factor(sample_map$Group, levels = group_levels)
write_tsv(sample_map %>% select(-sample_dir), file.path(out_dir, "00_<U+6837><U+672C><U+5206><U+7EC4><U+8868>.tsv"))

read_performance <- function(path, sample) {
  read_table(file.path(path, "performance.txt"), show_col_types = FALSE) %>%
    mutate(
      Sample_ID = sample, Delta_R2 = multi.R2 - intra.R2,
      Delta_RMSE = intra.RMSE - multi.RMSE
    )
}

read_coefficients <- function(path, sample) {
  read_table(file.path(path, "coefficients.txt"), show_col_types = FALSE) %>%
    mutate(Sample_ID = sample)
}

read_importance_file <- function(f, sample) {
  nm <- basename(f)
  hit <- str_match(nm, "^importances_(.+)_(intra|para\\.2|para\\.5)\\.txt$")
  if (anyNA(hit)) {
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE) %>%
    transmute(
      Sample_ID = sample, Predictor = target,
      Target = hit[, 2], View = hit[, 3], Importance = imp
    )
}

performance <- map2_dfr(sample_map$sample_dir, sample_map$Sample_ID, read_performance) %>%
  left_join(sample_map %>% select(Sample_ID, Group), by = "Sample_ID")
coefficients <- map2_dfr(sample_map$sample_dir, sample_map$Sample_ID, read_coefficients) %>%
  left_join(sample_map %>% select(Sample_ID, Group), by = "Sample_ID")
importance <- map2_dfr(sample_map$sample_dir, sample_map$Sample_ID, function(path, sample) {
  fs <- list.files(path, pattern = "^importances_.*_(intra|para\\.2|para\\.5)\\.txt$", full.names = TRUE)
  map_dfr(fs, read_importance_file, sample = sample)
}) %>% left_join(sample_map %>% select(Sample_ID, Group), by = "Sample_ID")

importance$View <- factor(importance$View, levels = view_levels)
write_tsv(performance, file.path(out_dir, "01_performance_<U+9010><U+6837><U+672C>.tsv"))
write_tsv(coefficients, file.path(out_dir, "02_view_coefficients_<U+9010><U+6837><U+672C>.tsv"))
write_tsv(importance, file.path(out_dir, "03_importance_<U+9010><U+6837><U+672C><U+957F><U+8868>.tsv"))

# 1. <U+7A7A><U+95F4><U+89C6><U+56FE><U+76F8><U+5BF9><U+540C> spot <U+6A21><U+578B><U+5E26><U+6765><U+7684> R2 <U+589E><U+76CA> ----------------------------------
perf_main <- performance %>%
  filter(as.character(Group) %in% main_groups) %>%
  mutate(Group = factor(as.character(Group), levels = main_groups))
p_gain <- ggplot(perf_main, aes(Group, Delta_R2, colour = Group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.45) +
  geom_jitter(width = 0.14, alpha = 0.45, size = 0.9) +
  facet_wrap(~target, scales = "free_y", ncol = 6) +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(
    title = "Predictive gain after adding spatial-neighborhood views",
    x = NULL, y = expression(Delta * R^2 ~ "(multi - intra, percentage points)")
  ) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(out_dir, "Figure1_<U+7A7A><U+95F4><U+6A21><U+578B>_DeltaR2_<U+5168><U+90E8><U+9776><U+7EC6><U+80DE>.pdf"), p_gain, width = 14, height = 11, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure1_<U+7A7A><U+95F4><U+6A21><U+578B>_DeltaR2_<U+5168><U+90E8><U+9776><U+7EC6><U+80DE>.png"), p_gain, width = 14, height = 11, dpi = 300)

gain_summary <- perf_main %>%
  group_by(Group, target) %>%
  summarise(
    n = n(), median_Delta_R2 = median(Delta_R2, na.rm = TRUE),
    Q1 = quantile(Delta_R2, .25, na.rm = TRUE), Q3 = quantile(Delta_R2, .75, na.rm = TRUE), .groups = "drop"
  )
write_tsv(gain_summary, file.path(out_dir, "Figure1_<U+7A7A><U+95F4><U+6A21><U+578B>_DeltaR2_<U+6C47><U+603B>.tsv"))

# 2. <U+6BCF><U+4E2A><U+4E3B><U+8981><U+75BE><U+75C5><U+7EC4><U+3001><U+6BCF><U+4E2A><U+7A7A><U+95F4><U+5C3A><U+5EA6><U+7684><U+5B8C><U+6574> interaction heatmap -------------------
agg <- importance %>%
  filter(as.character(Group) %in% main_groups) %>%
  group_by(Group, View, Predictor, Target) %>%
  summarise(
    n = sum(!is.na(Importance)), median_importance = median(Importance, na.rm = TRUE),
    mean_importance = mean(Importance, na.rm = TRUE), .groups = "drop"
  )
write_tsv(agg, file.path(out_dir, "Figure2_<U+5404><U+7EC4><U+5404><U+5C3A><U+5EA6>_importance<U+6C47><U+603B>.tsv"))

cell_order <- sort(unique(c(importance$Predictor, importance$Target)))
heat_limit <- quantile(agg$median_importance, .99, na.rm = TRUE)
plot_heat <- function(dat, ttl, divergent = FALSE, lim = NULL) {
  dat <- dat %>% mutate(Predictor = factor(Predictor, levels = cell_order), Target = factor(Target, levels = rev(cell_order)))
  p <- ggplot(dat, aes(Predictor, Target, fill = value)) +
    geom_tile() +
    coord_equal() +
    labs(title = ttl, x = "Predictor cell subtype", y = "Target cell subtype", fill = NULL) +
    theme_pub(8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  if (divergent) {
    p + scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-lim, lim), oob = squish
    )
  } else {
    p + scale_fill_gradient(low = "white", high = "#8B1A1A", limits = c(0, lim), oob = squish)
  }
}

pdf(file.path(out_dir, "Figure2_<U+4E3B><U+8981><U+75BE><U+75C5><U+7EC4>_para2_para5<U+5B8C><U+6574><U+4E92><U+4F5C><U+70ED><U+56FE>.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in main_groups) {
  for (v in c("para.2", "para.5")) {
    dd <- agg %>%
      filter(as.character(Group) == g, View == v) %>%
      rename(value = median_importance)
    print(plot_heat(dd, paste0(g, " | ", view_labels[[v]], " | median importance"), lim = heat_limit))
  }
}
dev.off()

# 3. l=5 <U+4E0E> l=2 <U+7684><U+76F8><U+5BF9><U+5C3A><U+5EA6><U+504F><U+597D> ---------------------------------------------
scale_diff <- agg %>%
  filter(View %in% c("para.2", "para.5")) %>%
  select(Group, View, Predictor, Target, median_importance) %>%
  pivot_wider(names_from = View, values_from = median_importance) %>%
  mutate(Diff_l5_minus_l2 = para.5 - para.2)
write_tsv(scale_diff, file.path(out_dir, "Figure3_l5<U+51CF>l2_<U+5C3A><U+5EA6><U+504F><U+597D>.tsv"))
diff_lim <- quantile(abs(scale_diff$Diff_l5_minus_l2), .99, na.rm = TRUE)
pdf(file.path(out_dir, "Figure3_l5<U+51CF>l2_<U+5C3A><U+5EA6><U+504F><U+597D><U+70ED><U+56FE>.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in main_groups) {
  dd <- scale_diff %>%
    filter(as.character(Group) == g) %>%
    rename(value = Diff_l5_minus_l2)
  print(plot_heat(dd, paste0(g, " | relative scale preference (l=5 minus l=2)"), TRUE, diff_lim))
}
dev.off()

# 4. <U+91CD><U+70B9><U+6210><U+7EA4><U+7EF4><U+7EC6><U+80DE><U+2014><U+9AD3><U+7CFB><U+7EC6><U+80DE><U+65B9><U+5411><U+6027><U+5173><U+7CFB><U+7684><U+4E09><U+89C6><U+56FE><U+5206><U+5E03> --------------------------
key_pairs <- tribble(
  ~Predictor, ~Target,
  "SPP1MP", "FN1Fib", "FN1Fib", "SPP1MP",
  "FOLR2MP", "FN1Fib", "FN1Fib", "FOLR2MP",
  "CPMMP", "FN1Fib", "FN1Fib", "CPMMP"
) %>% mutate(Pair = paste(Predictor, "<U+2192>", Target))
key_pair_levels <- key_pairs$Pair

# <U+5728><U+7B5B><U+9009><U+91CD><U+70B9><U+7EC6><U+80DE><U+5BF9><U+4E4B><U+524D><U+8FDB><U+884C><U+6807><U+51C6><U+5316><U+FF1A><U+6BCF><U+4E2A><U+6837><U+672C> <U+00D7> target <U+00D7> view<U+4E2D><U+FF0C>
# <U+5168><U+90E8>predictor<U+7684><U+76F8><U+5BF9>importance<U+4E4B><U+548C><U+4E3A>2<U+3002><U+8BE5><U+503C><U+8868><U+793A><U+76F8><U+5BF9><U+8D21><U+732E><U+FF0C><U+4E0D><U+662F><U+539F><U+59CB>importance<U+3002>
importance_relative <- importance %>%
  group_by(Sample_ID, Target, View) %>%
  mutate(Importance_relative = ifelse(
    sum(Importance, na.rm = TRUE) > 0,
    Importance / sum(Importance, na.rm = TRUE) * 2,
    NA_real_
  )) %>%
  ungroup()

key_dat <- importance_relative %>%
  inner_join(key_pairs, by = c("Predictor", "Target")) %>%
  filter(as.character(Group) %in% main_groups) %>%
  mutate(
    Group = factor(as.character(Group), levels = main_groups),
    Pair = factor(Pair, levels = key_pair_levels),
    View_label = factor(view_labels[as.character(View)], levels = unname(view_labels))
  )
p_key <- ggplot(key_dat, aes(View_label, Importance_relative, group = Sample_ID, colour = Group)) +
  geom_line(alpha = .22, linewidth = .35) +
  geom_point(alpha = .38, size = .8) +
  stat_summary(aes(group = Group), fun = median, geom = "line", linewidth = 1.05) +
  stat_summary(aes(group = Group), fun = median, geom = "point", size = 2) +
  facet_grid(Pair ~ Group, scales = "free_y") +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(
    title = "Scale-dependent fibroblast-myeloid spatial dependencies",
    subtitle = "Relative importance normalized within each sample, target and view (sum across predictors = 2)",
    x = NULL, y = "Relative MISTy importance (sum = 2)"
  ) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "Figure4_<U+91CD><U+70B9>Fib_Myeloid<U+7EC6><U+80DE><U+5BF9>_<U+76F8><U+5BF9><U+6807><U+51C6><U+5316>importance.pdf"), p_key, width = 8.5, height = 11, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure4_<U+91CD><U+70B9>Fib_Myeloid<U+7EC6><U+80DE><U+5BF9>_<U+76F8><U+5BF9><U+6807><U+51C6><U+5316>importance.png"), p_key, width = 8.5, height = 11, dpi = 300)
write_tsv(key_dat, file.path(out_dir, "Figure4_<U+91CD><U+70B9>Fib_Myeloid<U+7EC6><U+80DE><U+5BF9>_<U+76F8><U+5BF9><U+6807><U+51C6><U+5316>importance.tsv"))

# l=2 vs l=5 <U+914D><U+5BF9> Wilcoxon<U+FF1B><U+68C0><U+9A8C><U+7684><U+662F><U+5C3A><U+5EA6><U+5DEE><U+5F02><U+FF0C><U+4E0D><U+89E3><U+91CA><U+4E3A><U+6B63>/<U+8D1F><U+5171><U+5B9A><U+4F4D><U+3002>
scale_tests <- importance %>%
  filter(View %in% c("para.2", "para.5"), as.character(Group) %in% main_groups) %>%
  select(Sample_ID, Group, Predictor, Target, View, Importance) %>%
  pivot_wider(names_from = View, values_from = Importance) %>%
  group_by(Group, Predictor, Target) %>%
  summarise(
    n = sum(complete.cases(para.2, para.5)),
    median_l2 = median(para.2, na.rm = TRUE), median_l5 = median(para.5, na.rm = TRUE),
    median_difference = median(para.5 - para.2, na.rm = TRUE),
    p_value = if (n >= 3 && any((para.5 - para.2) != 0, na.rm = TRUE)) {
      suppressWarnings(wilcox.test(para.5, para.2, paired = TRUE, exact = FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(FDR = p.adjust(p_value, "BH")) %>%
  ungroup()
write_tsv(scale_tests, file.path(out_dir, "Figure4_l2_vs_l5_<U+9010><U+7EC6><U+80DE><U+5BF9><U+914D><U+5BF9><U+68C0><U+9A8C>.tsv"))

# 5. <U+5404><U+75BE><U+75C5><U+6838><U+5FC3><U+72B6><U+6001><U+76F8><U+5BF9> Normal <U+7684><U+4E2D><U+4F4D> importance <U+5DEE><U+503C> ------------------------
contrast <- agg %>%
  filter(View %in% c("para.2", "para.5")) %>%
  select(Group, View, Predictor, Target, median_importance) %>%
  pivot_wider(names_from = Group, values_from = median_importance) %>%
  pivot_longer(cols = all_of(c("AMI-IZ", "ICM-FZ", "HCM")), names_to = "Disease", values_to = "Disease_median") %>%
  mutate(Difference_vs_Normal = Disease_median - Normal)
write_tsv(contrast, file.path(out_dir, "Figure5_<U+75BE><U+75C5><U+6838><U+5FC3><U+72B6><U+6001>_vs_Normal_importance<U+5DEE><U+503C>.tsv"))
contrast_lim <- quantile(abs(contrast$Difference_vs_Normal), .99, na.rm = TRUE)
pdf(file.path(out_dir, "Figure5_<U+75BE><U+75C5><U+6838><U+5FC3><U+72B6><U+6001>_vs_Normal_<U+5DEE><U+5F02><U+70ED><U+56FE>.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in c("AMI-IZ", "ICM-FZ", "HCM")) {
  for (v in c("para.2", "para.5")) {
    dd <- contrast %>%
      filter(Disease == g, View == v) %>%
      rename(value = Difference_vs_Normal)
    print(plot_heat(dd, paste0(g, " vs Normal | ", view_labels[[v]]), TRUE, contrast_lim))
  }
}
dev.off()

# 6. <U+7A33><U+5B9A><U+8FB9><U+8868><U+FF1A>importance <U+8F83><U+5F3A><U+4E14><U+5728><U+591A><U+4E2A><U+6837><U+672C><U+4E2D><U+91CD><U+590D><U+51FA><U+73B0> --------------------------
# <U+6BCF><U+4E2A><U+6837><U+672C>/<U+5C3A><U+5EA6><U+5185><U+90E8><U+4F7F><U+7528><U+5176><U+5168><U+90E8><U+8FB9><U+7684><U+524D> 10% <U+4F5C><U+4E3A><U+201C><U+8F83><U+5F3A><U+201D><U+FF0C><U+907F><U+514D><U+4EFB><U+610F><U+56FA><U+5B9A> importance <U+9608><U+503C><U+3002>
stable <- importance %>%
  filter(View %in% c("para.2", "para.5"), as.character(Group) %in% main_groups) %>%
  group_by(Sample_ID, View) %>%
  mutate(sample_q90 = quantile(Importance, .90, na.rm = TRUE), Strong = Importance >= sample_q90) %>%
  ungroup() %>%
  group_by(Group, View, Predictor, Target) %>%
  summarise(
    n_samples = n(), prevalence = mean(Strong, na.rm = TRUE),
    median_importance = median(Importance, na.rm = TRUE), .groups = "drop"
  ) %>%
  arrange(Group, View, desc(prevalence), desc(median_importance))
write_tsv(stable, file.path(out_dir, "Figure6_<U+7A33><U+5B9A><U+7A7A><U+95F4><U+4F9D><U+8D56><U+8FB9><U+8868>.tsv"))

# <U+8F93><U+51FA><U+5206><U+6790><U+8BF4><U+660E><U+548C><U+8D28><U+91CF><U+6838><U+5BF9>
qc <- c(
  paste0("<U+6837><U+672C><U+6570>: ", n_distinct(importance$Sample_ID)),
  paste0("<U+7EC6><U+80DE><U+4E9A><U+578B><U+6570>: ", length(cell_order)),
  paste0("importance<U+8BB0><U+5F55><U+6570>: ", nrow(importance)),
  paste0("<U+4E3B><U+8981><U+6BD4><U+8F83><U+7EC4>: ", paste(main_groups, collapse = ", ")),
  "<U+89E3><U+91CA><U+8FB9><U+754C>: importance<U+4E3A><U+65B9><U+5411><U+6027><U+9884><U+6D4B><U+8D21><U+732E><U+FF0C><U+4E0D><U+63D0><U+4F9B><U+6B63><U+8D1F><U+76F8><U+5173><U+65B9><U+5411><U+FF0C><U+4E5F><U+4E0D><U+80FD><U+5355><U+72EC><U+8BC1><U+660E><U+7269><U+7406><U+63A5><U+89E6><U+3002>",
  "<U+5C3A><U+5EA6><U+8FB9><U+754C>: para.2<U+548C>para.5<U+662F><U+91CD><U+53E0><U+7684><U+9AD8><U+65AF><U+52A0><U+6743><U+90BB><U+57DF><U+FF1B>l<U+503C><U+9700><U+7ED3><U+5408><U+8F93><U+5165>geometry<U+5355><U+4F4D><U+624D><U+80FD><U+6362><U+7B97><U+4E3A><U+5B9E><U+9645><U+5FAE><U+7C73><U+3002>",
  "<U+7EDF><U+8BA1><U+5355><U+4F4D>: <U+6837><U+672C><U+3002>l=2<U+4E0E>l=5<U+4F7F><U+7528><U+6837><U+672C><U+5185><U+914D><U+5BF9>Wilcoxon<U+68C0><U+9A8C><U+5E76><U+5728><U+75BE><U+75C5><U+7EC4><U+5185>BH<U+6821><U+6B63><U+3002>",
  "<U+7A33><U+5B9A><U+8FB9>: <U+6BCF><U+4E2A><U+6837><U+672C><U+6BCF><U+4E2A><U+5C3A><U+5EA6><U+5185><U+90E8>importance<U+524D>10%<U+FF0C><U+518D><U+8BA1><U+7B97><U+8DE8><U+6837><U+672C><U+51FA><U+73B0><U+6BD4><U+4F8B><U+3002>"
)
writeLines(qc, file.path(out_dir, "README_<U+7ED3><U+679C><U+89E3><U+91CA><U+4E0E><U+6587><U+4EF6><U+8BF4><U+660E>.txt"), useBytes = TRUE)

message("<U+5B8C><U+6210><U+3002><U+7ED3><U+679C><U+76EE><U+5F55><U+FF1A>", out_dir)
