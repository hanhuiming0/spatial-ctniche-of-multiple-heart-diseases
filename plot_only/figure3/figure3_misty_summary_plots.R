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

# MISTy 不同空间尺度的细胞亚型空间依赖结果汇总
# 说明：importance 表示 predictor 对 target 的预测贡献，不代表相关方向或物理接触。
# para.2/para.5 为高斯加权空间视图，二者不是互斥的距离环。

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) stop("请使用 Rscript 运行本脚本，以便自动定位输入目录。")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- dirname(script_file)
input_dir <- file.path(root, "SubType_MISTy2")
project_root <- normalizePath(file.path(root, "../../../../.."), winslash = "/", mustWork = TRUE)
meta_file <- file.path(project_root, "2-数据", "Sample_Info_2025.8.20.txt")
out_dir <- file.path(root, "SubType_MISTy2_不同距离汇总图")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(input_dir), file.exists(meta_file))

group_levels <- c("Normal", "AMI-IZ", "AMI-BZ", "AMI-RZ",
                  "ICM-FZ", "ICM-IZ", "ICM-BZ", "ICM-RZ", "HCM")
main_groups <- c("Normal", "AMI-IZ", "ICM-FZ", "HCM")
view_levels <- c("intra", "para.2", "para.5")
view_labels <- c(intra = "Same spot", para.2 = "l=2", para.5 = "l=5")
group_cols <- c("Normal" = "#6B7280", "AMI-IZ" = "#D55E00",
                "ICM-FZ" = "#0072B2", "HCM" = "#8E5AA9")

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25),
          strip.background = element_rect(fill = "grey95", colour = "grey65"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"), legend.position = "right")
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
if (anyNA(sample_map$Group)) stop("以下样本未匹配样本信息：", paste(sample_map$Sample_ID[is.na(sample_map$Group)], collapse = ", "))
sample_map$Group <- factor(sample_map$Group, levels = group_levels)
write_tsv(sample_map %>% select(-sample_dir), file.path(out_dir, "00_样本分组表.tsv"))

read_performance <- function(path, sample) {
  read_table(file.path(path, "performance.txt"), show_col_types = FALSE) %>%
    mutate(Sample_ID = sample, Delta_R2 = multi.R2 - intra.R2,
           Delta_RMSE = intra.RMSE - multi.RMSE)
}

read_coefficients <- function(path, sample) {
  read_table(file.path(path, "coefficients.txt"), show_col_types = FALSE) %>%
    mutate(Sample_ID = sample)
}

read_importance_file <- function(f, sample) {
  nm <- basename(f)
  hit <- str_match(nm, "^importances_(.+)_(intra|para\\.2|para\\.5)\\.txt$")
  if (anyNA(hit)) return(NULL)
  read_csv(f, show_col_types = FALSE) %>%
    transmute(Sample_ID = sample, Predictor = target,
              Target = hit[, 2], View = hit[, 3], Importance = imp)
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
write_tsv(performance, file.path(out_dir, "01_performance_逐样本.tsv"))
write_tsv(coefficients, file.path(out_dir, "02_view_coefficients_逐样本.tsv"))
write_tsv(importance, file.path(out_dir, "03_importance_逐样本长表.tsv"))

# 1. 空间视图相对同 spot 模型带来的 R2 增益 ----------------------------------
perf_main <- performance %>% filter(as.character(Group) %in% main_groups) %>%
  mutate(Group = factor(as.character(Group), levels = main_groups))
p_gain <- ggplot(perf_main, aes(Group, Delta_R2, colour = Group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.45) +
  geom_jitter(width = 0.14, alpha = 0.45, size = 0.9) +
  facet_wrap(~target, scales = "free_y", ncol = 6) +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(title = "Predictive gain after adding spatial-neighborhood views",
       x = NULL, y = expression(Delta*R^2~"(multi - intra, percentage points)")) +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(out_dir, "Figure1_空间模型_DeltaR2_全部靶细胞.pdf"), p_gain, width = 14, height = 11, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure1_空间模型_DeltaR2_全部靶细胞.png"), p_gain, width = 14, height = 11, dpi = 300)

gain_summary <- perf_main %>% group_by(Group, target) %>%
  summarise(n = n(), median_Delta_R2 = median(Delta_R2, na.rm = TRUE),
            Q1 = quantile(Delta_R2, .25, na.rm = TRUE), Q3 = quantile(Delta_R2, .75, na.rm = TRUE), .groups = "drop")
write_tsv(gain_summary, file.path(out_dir, "Figure1_空间模型_DeltaR2_汇总.tsv"))

# 2. 每个主要疾病组、每个空间尺度的完整 interaction heatmap -------------------
agg <- importance %>% filter(as.character(Group) %in% main_groups) %>%
  group_by(Group, View, Predictor, Target) %>%
  summarise(n = sum(!is.na(Importance)), median_importance = median(Importance, na.rm = TRUE),
            mean_importance = mean(Importance, na.rm = TRUE), .groups = "drop")
write_tsv(agg, file.path(out_dir, "Figure2_各组各尺度_importance汇总.tsv"))

cell_order <- sort(unique(c(importance$Predictor, importance$Target)))
heat_limit <- quantile(agg$median_importance, .99, na.rm = TRUE)
plot_heat <- function(dat, ttl, divergent = FALSE, lim = NULL) {
  dat <- dat %>% mutate(Predictor = factor(Predictor, levels = cell_order), Target = factor(Target, levels = rev(cell_order)))
  p <- ggplot(dat, aes(Predictor, Target, fill = value)) + geom_tile() + coord_equal() +
    labs(title = ttl, x = "Predictor cell subtype", y = "Target cell subtype", fill = NULL) +
    theme_pub(8) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  if (divergent) p + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                                           limits = c(-lim, lim), oob = squish)
  else p + scale_fill_gradient(low = "white", high = "#8B1A1A", limits = c(0, lim), oob = squish)
}

pdf(file.path(out_dir, "Figure2_主要疾病组_para2_para5完整互作热图.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in main_groups) for (v in c("para.2", "para.5")) {
  dd <- agg %>% filter(as.character(Group) == g, View == v) %>% rename(value = median_importance)
  print(plot_heat(dd, paste0(g, " | ", view_labels[[v]], " | median importance"), lim = heat_limit))
}
dev.off()

# 3. l=5 与 l=2 的相对尺度偏好 ---------------------------------------------
scale_diff <- agg %>% filter(View %in% c("para.2", "para.5")) %>%
  select(Group, View, Predictor, Target, median_importance) %>%
  pivot_wider(names_from = View, values_from = median_importance) %>%
  mutate(Diff_l5_minus_l2 = para.5 - para.2)
write_tsv(scale_diff, file.path(out_dir, "Figure3_l5减l2_尺度偏好.tsv"))
diff_lim <- quantile(abs(scale_diff$Diff_l5_minus_l2), .99, na.rm = TRUE)
pdf(file.path(out_dir, "Figure3_l5减l2_尺度偏好热图.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in main_groups) {
  dd <- scale_diff %>% filter(as.character(Group) == g) %>% rename(value = Diff_l5_minus_l2)
  print(plot_heat(dd, paste0(g, " | relative scale preference (l=5 minus l=2)"), TRUE, diff_lim))
}
dev.off()

# 4. 重点成纤维细胞—髓系细胞方向性关系的三视图分布 --------------------------
key_pairs <- tribble(
  ~Predictor, ~Target,
  "SPP1MP", "FN1Fib", "FN1Fib", "SPP1MP",
  "FOLR2MP", "FN1Fib", "FN1Fib", "FOLR2MP",
  "CPMMP", "FN1Fib", "FN1Fib", "CPMMP"
) %>% mutate(Pair = paste(Predictor, "→", Target))
key_pair_levels <- key_pairs$Pair

# 在筛选重点细胞对之前进行标准化：每个样本 × target × view中，
# 全部predictor的相对importance之和为2。该值表示相对贡献，不是原始importance。
importance_relative <- importance %>%
  group_by(Sample_ID, Target, View) %>%
  mutate(Importance_relative = ifelse(
    sum(Importance, na.rm = TRUE) > 0,
    Importance / sum(Importance, na.rm = TRUE) * 2,
    NA_real_
  )) %>%
  ungroup()

key_dat <- importance_relative %>% inner_join(key_pairs, by = c("Predictor", "Target")) %>%
  filter(as.character(Group) %in% main_groups) %>%
  mutate(Group = factor(as.character(Group), levels = main_groups),
         Pair = factor(Pair, levels = key_pair_levels),
         View_label = factor(view_labels[as.character(View)], levels = unname(view_labels)))
p_key <- ggplot(key_dat, aes(View_label, Importance_relative, group = Sample_ID, colour = Group)) +
  geom_line(alpha = .22, linewidth = .35) + geom_point(alpha = .38, size = .8) +
  stat_summary(aes(group = Group), fun = median, geom = "line", linewidth = 1.05) +
  stat_summary(aes(group = Group), fun = median, geom = "point", size = 2) +
  facet_grid(Pair ~ Group, scales = "free_y") + scale_colour_manual(values = group_cols, guide = "none") +
  labs(title = "Scale-dependent fibroblast-myeloid spatial dependencies",
       subtitle = "Relative importance normalized within each sample, target and view (sum across predictors = 2)",
       x = NULL, y = "Relative MISTy importance (sum = 2)") + theme_pub(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "Figure4_重点Fib_Myeloid细胞对_相对标准化importance.pdf"), p_key, width = 8.5, height = 11, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure4_重点Fib_Myeloid细胞对_相对标准化importance.png"), p_key, width = 8.5, height = 11, dpi = 300)
write_tsv(key_dat, file.path(out_dir, "Figure4_重点Fib_Myeloid细胞对_相对标准化importance.tsv"))

# l=2 vs l=5 配对 Wilcoxon；检验的是尺度差异，不解释为正/负共定位。
scale_tests <- importance %>% filter(View %in% c("para.2", "para.5"), as.character(Group) %in% main_groups) %>%
  select(Sample_ID, Group, Predictor, Target, View, Importance) %>%
  pivot_wider(names_from = View, values_from = Importance) %>%
  group_by(Group, Predictor, Target) %>%
  summarise(n = sum(complete.cases(para.2, para.5)),
            median_l2 = median(para.2, na.rm = TRUE), median_l5 = median(para.5, na.rm = TRUE),
            median_difference = median(para.5 - para.2, na.rm = TRUE),
            p_value = if (n >= 3 && any((para.5 - para.2) != 0, na.rm = TRUE))
              suppressWarnings(wilcox.test(para.5, para.2, paired = TRUE, exact = FALSE)$p.value) else NA_real_,
            .groups = "drop") %>% group_by(Group) %>% mutate(FDR = p.adjust(p_value, "BH")) %>% ungroup()
write_tsv(scale_tests, file.path(out_dir, "Figure4_l2_vs_l5_逐细胞对配对检验.tsv"))

# 5. 各疾病核心状态相对 Normal 的中位 importance 差值 ------------------------
contrast <- agg %>% filter(View %in% c("para.2", "para.5")) %>%
  select(Group, View, Predictor, Target, median_importance) %>%
  pivot_wider(names_from = Group, values_from = median_importance) %>%
  pivot_longer(cols = all_of(c("AMI-IZ", "ICM-FZ", "HCM")), names_to = "Disease", values_to = "Disease_median") %>%
  mutate(Difference_vs_Normal = Disease_median - Normal)
write_tsv(contrast, file.path(out_dir, "Figure5_疾病核心状态_vs_Normal_importance差值.tsv"))
contrast_lim <- quantile(abs(contrast$Difference_vs_Normal), .99, na.rm = TRUE)
pdf(file.path(out_dir, "Figure5_疾病核心状态_vs_Normal_差异热图.pdf"), width = 9.5, height = 8.5, onefile = TRUE)
for (g in c("AMI-IZ", "ICM-FZ", "HCM")) for (v in c("para.2", "para.5")) {
  dd <- contrast %>% filter(Disease == g, View == v) %>% rename(value = Difference_vs_Normal)
  print(plot_heat(dd, paste0(g, " vs Normal | ", view_labels[[v]]), TRUE, contrast_lim))
}
dev.off()

# 6. 稳定边表：importance 较强且在多个样本中重复出现 --------------------------
# 每个样本/尺度内部使用其全部边的前 10% 作为“较强”，避免任意固定 importance 阈值。
stable <- importance %>% filter(View %in% c("para.2", "para.5"), as.character(Group) %in% main_groups) %>%
  group_by(Sample_ID, View) %>% mutate(sample_q90 = quantile(Importance, .90, na.rm = TRUE), Strong = Importance >= sample_q90) %>%
  ungroup() %>% group_by(Group, View, Predictor, Target) %>%
  summarise(n_samples = n(), prevalence = mean(Strong, na.rm = TRUE),
            median_importance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Group, View, desc(prevalence), desc(median_importance))
write_tsv(stable, file.path(out_dir, "Figure6_稳定空间依赖边表.tsv"))

# 输出分析说明和质量核对
qc <- c(
  paste0("样本数: ", n_distinct(importance$Sample_ID)),
  paste0("细胞亚型数: ", length(cell_order)),
  paste0("importance记录数: ", nrow(importance)),
  paste0("主要比较组: ", paste(main_groups, collapse = ", ")),
  "解释边界: importance为方向性预测贡献，不提供正负相关方向，也不能单独证明物理接触。",
  "尺度边界: para.2和para.5是重叠的高斯加权邻域；l值需结合输入geometry单位才能换算为实际微米。",
  "统计单位: 样本。l=2与l=5使用样本内配对Wilcoxon检验并在疾病组内BH校正。",
  "稳定边: 每个样本每个尺度内部importance前10%，再计算跨样本出现比例。"
)
writeLines(qc, file.path(out_dir, "README_结果解释与文件说明.txt"), useBytes = TRUE)

message("完成。结果目录：", out_dir)