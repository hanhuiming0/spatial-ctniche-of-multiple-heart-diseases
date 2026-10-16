ST2_ROOT <- Sys.getenv('ST2_ROOT', unset = '/path/to/ST2')

# ctniche 信息
cluster_info <- readRDS(file.path(ST2_ROOT, 'analysis', 'ST', 'Niche', 'cluster_info', 'cluster_info_split5.RDS'))
cluster_info35 <- cluster_info[,c('spot_id','k_35')]
colnames(cluster_info35) <- c("row_id","ct_niche")
cluster_info35$ct_niche <- paste0("ctniche_", cluster_info35$ct_niche)

norm_weights <- readRDS(file.path(ST2_ROOT, 'analysis', 'ST', 'RCTD', 'cell_type_proportions.RDS')) # spot-SubType矩阵
norm_weights <- as.data.frame(as.matrix(norm_weights))

# 计算每个ctniche的细胞类型比例，绘制扇形图
library(ggplot2)
library(dplyr)
library(tidyr)

# 合并ctniche信息和细胞类型比例
data <- merge(cluster_info35, norm_weights, by.x = "row_id", by.y = "row.names")
data_long <- data %>%
    pivot_longer(cols = -c(row_id, ct_niche), names_to = "cell_type", values_to = "proportion")

# 计算每个ctniche的细胞类型比例总和
data_long[1:10,]
{
  data_long$sample <- sub('_.*','',data_long$row_id)

  Sample_Info <- read.delim(file.path(ST2_ROOT, 'data', 'sample_info', 'Sample_Info_2025.8.20.txt'))
  Sample_Info <- Sample_Info[which(Sample_Info$Dataset_Type == 'Visium'),]
  niche_cols = c('#A7C5DD','#F9C7A1','#9DCB9C','#CACACA','#C8B7E2',
                 '#C6ACA4','#F5B6E3','#EA9292','#E5EAA5','#A7E0E7',
                 '#AE6D73','#01549A')

  Ischaemic <- Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'AMI' & Sample_Info$Disease_State == 'IZ')]
  Fibrotic <- Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'ICM' & Sample_Info$Disease_State == 'FZ')]
  HCM <- Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'HCM')]
  Border <- Sample_Info$Sample_ID[which(Sample_Info$Disease_State %in% 'BZ')]
  Remote <- Sample_Info$Sample_ID[which(Sample_Info$Disease_State %in% 'RZ')]

  data_long <- data_long[-which(data_long$sample == 'ACH0010'),]
  data_long$ST_Type <- 'Normal'
  data_long$ST_Type[which(data_long$sample %in% Ischaemic)] <- 'Ischaemic'
  data_long$ST_Type[which(data_long$sample %in% Fibrotic)] <- 'Fibrotic'
  data_long$ST_Type[which(data_long$sample %in% HCM)] <- 'HCM'
  data_long$ST_Type[which(data_long$sample %in% Border)] <- 'Border'
  data_long$ST_Type[which(data_long$sample %in% Remote)] <- 'Remote'
}
data_summary <- data_long %>%
    group_by(ct_niche, cell_type) %>%
    summarise(total_proportion = sum(proportion), .groups = 'drop')
{
  data_summary.Ischaemic <- data_long[which(data_long$ST_Type == 'Ischaemic'),] %>%
    group_by(ct_niche, cell_type) %>%
    summarise(total_proportion = sum(proportion), .groups = 'drop')

  data_summary.Fibrotic <- data_long[which(data_long$ST_Type == 'Fibrotic'),] %>%
    group_by(ct_niche, cell_type) %>%
    summarise(total_proportion = sum(proportion), .groups = 'drop')

  data_summary.HCM <- data_long[which(data_long$ST_Type == 'HCM'),] %>%
    group_by(ct_niche, cell_type) %>%
    summarise(total_proportion = sum(proportion), .groups = 'drop')

}
# 每个ctniche的细胞类型比例总和为1
data_summary <- data_summary %>%
    group_by(ct_niche) %>%
    mutate(total_proportion = total_proportion / sum(total_proportion)) %>%
    ungroup()
{
  data_summary.Ischaemic <- data_summary.Ischaemic %>%
    group_by(ct_niche) %>%
    mutate(total_proportion = total_proportion / sum(total_proportion)) %>%
    ungroup()

  data_summary.Fibrotic <- data_summary.Fibrotic %>%
    group_by(ct_niche) %>%
    mutate(total_proportion = total_proportion / sum(total_proportion)) %>%
    ungroup()

  data_summary.HCM <- data_summary.HCM %>%
    group_by(ct_niche) %>%
    mutate(total_proportion = total_proportion / sum(total_proportion)) %>%
    ungroup()

    }
# 绘制扇形图
SubType_Color <- c(
  '#d985b4','#f8caaa',
  '#87cabf','#ed7c6f','#bab6d5',
  '#7eaac9','#f3ae63','#acd16a','#f5c9dd','#d8d8d8',
  '#dac8e0','#f9ecab','#c7e2c0','#ce5d15',
  '#e47faf','#e1edc6','#f3afaa',
  '#AB3FFD',"#E278C1",
  "#AFC6E7","#8D554F","#FF7D0E"
)
names(SubType_Color) <- c(
  'Normal_vCM','Stress_vCM',
  'THBS4+ Fib', 'FN1+ Fib', 'HIF3A Fib',
  'Arterial_endo','Capillary_endo','Venous_endo','Lymphatic_endo','Endocardial_endo',
  'CPM+ MP','FOLR2+ MP','Monocyte','SPP1+ MP',
  'CD8+T_cells','CD4+T_cells','NK',
  'Pericyte','SMC',
  'Neuronal','Mast','Adipocyte'
)
data_summary$ct_niche <- factor(data_summary$ct_niche, levels = paste0("ctniche_", 1:12))
ggplot(data_summary, aes(x = "", y = total_proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = SubType_Color) +
  facet_wrap(~ ct_niche) +
  theme_void() +
  labs(title = "Cell Type Proportions in Each ct_niche") +
  theme(legend.position = "right")
ggsave(file.path(ST2_ROOT, 'Niche', 'ctniche_celltype_proportion_pie_chart.pdf'),
width = 12, height = 8 )

# 20260624
data_summary.Ischaemic$ct_niche <- factor(data_summary.Ischaemic$ct_niche, levels = paste0("ctniche_", 1:12))
# ggplot(data_summary.Ischaemic, aes(x = "", y = total_proportion, fill = cell_type)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y", start = 0) +
#   scale_fill_manual(values = SubType_Color) +
#   facet_wrap(~ ct_niche) +
#   theme_void() +
#   labs(title = "Cell Type Proportions in Each ct_niche") +
#   theme(legend.position = "right")
ggplot(data_summary.Ischaemic, aes(x = 2, y = total_proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = SubType_Color) +
  facet_wrap(~ ct_niche) +
  theme_void() +
  labs(title = "Ischaemic Cell Type Proportions in Each ct_niche") +
  theme(legend.position = "right")
ggsave(file.path(ST2_ROOT, 'Niche', 'AMI.ctniche_celltype_proportion.pdf'),
       width = 6.8, height = 3.8 )

data_summary.Fibrotic$ct_niche <- factor(data_summary.Fibrotic$ct_niche, levels = paste0("ctniche_", 1:12))
ggplot(data_summary.Fibrotic, aes(x = 2, y = total_proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = SubType_Color) +
  facet_wrap(~ ct_niche) +
  theme_void() +
  labs(title = "Fibrotic Cell Type Proportions in Each ct_niche") +
  theme(legend.position = "right")
ggsave(file.path(ST2_ROOT, 'Niche', 'ICM.ctniche_celltype_proportion.pdf'),
       width = 6.8, height = 3.8 )

data_summary.HCM$ct_niche <- factor(data_summary.HCM$ct_niche, levels = paste0("ctniche_", 1:12))
ggplot(data_summary.HCM, aes(x = 2, y = total_proportion, fill = cell_type)) +
         geom_bar(stat = "identity", width = 1) +
         coord_polar("y", start = 0) +
         xlim(0.5, 2.5) +
         scale_fill_manual(values = SubType_Color) +
         facet_wrap(~ ct_niche) +
         theme_void() +
         labs(title = "HCM Cell Type Proportions in Each ct_niche") +
         theme(legend.position = "right")
ggsave(file.path(ST2_ROOT, 'Niche', 'HCM.ctniche_celltype_proportion.pdf'),
       width = 6.8, height = 3.8 )

print(data_summary.HCM[which(data_summary.HCM$ct_niche == 'ctniche_8'), ], n = 22)