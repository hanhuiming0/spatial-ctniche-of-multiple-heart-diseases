PROJECT_ROOT <- Sys.getenv('PROJECT_ROOT', unset = '.')

setwd(file.path(PROJECT_ROOT, '2-数据'))
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)

cluster_info <- readRDS("cluster_info_split5.RDS")
cluster_info <- as.data.frame(cluster_info)
rownames(cluster_info) <- cluster_info$spot_id
cluster_info_35 <- as.data.frame(cluster_info[,'k_35',drop=F])
cluster_info_35$Sample <- sub('_.*','',rownames(cluster_info_35))
cluster_info_35$k_35 = paste0('ctniche_', cluster_info_35$k_35)
cluster_info_35$k_35 <- factor(cluster_info_35$k_35, levels = paste0('ctniche_', 1:12))
unique(cluster_info_35$Sample)

Sample_Info <- read.delim("Sample_Info_2025.8.20.txt")
Sample_Info <- Sample_Info[which(Sample_Info$Dataset_Type == 'Visium'),]

Cluster_Sample <- as.data.frame(table(cluster_info_35$k_35,cluster_info_35$Sample))
colnames(Cluster_Sample) <- c('ctniche','sample','spot')

sample_rank <- c(Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'Normal' & Sample_Info$Sampling_Tissue == 'LV')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'Normal' & Sample_Info$Sampling_Tissue == 'S')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'AMI' & Sample_Info$Disease_State == 'RZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'AMI' & Sample_Info$Disease_State == 'BZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'ICM' & Sample_Info$Disease_State == 'RZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'ICM' & Sample_Info$Disease_State == 'BZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'AMI' & Sample_Info$Disease_State == 'IZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'ICM' & Sample_Info$Disease_State == 'IZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'ICM' & Sample_Info$Disease_State == 'FZ')],
                 Sample_Info$Sample_ID[which(Sample_Info$Sample_Type == 'HCM')])
Cluster_Sample$sample <- factor(Cluster_Sample$sample, levels = sample_rank)

ctniche_color = c('#A7C5DD','#F9C7A1','#9DCB9C','#CACACA','#C8B7E2',
                              '#C6ACA4','#F5B6E3','#EA9292','#E5EAA5','#A7E0E7',
                              '#AE6D73','#01549A')
p_bar <- ggplot(Cluster_Sample, aes(x = sample, y = spot, fill = ctniche)) +
  geom_bar(stat = "identity", width = 0.95, position = "fill",
           colour = "black", alpha = 0.8, linewidth = 0.1) +
  scale_fill_manual(values = ctniche_color) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.75),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.y = element_text(vjust = .5),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
  labs(x = "Sample", y = "Proportion") +
  scale_y_continuous(expand = c(0,0))

annotation_df <- data.frame(
  Sample = sample_rank,
  SampleType = factor(Sample_Info$Sample_Type[match(sample_rank,Sample_Info$Sample_ID)],
                      levels = c('Normal','AMI','ICM','HCM')),
  TissueState = Sample_Info$Disease_State[match(sample_rank,Sample_Info$Sample_ID)],
  Sampling_Tissue = Sample_Info$Sampling_Tissue[match(sample_rank,Sample_Info$Sample_ID)]
)
anno_colors <- list(
  SampleType = c(Normal="#d0e2c0", HCM="#fded95", AMI="#49c2d9", ICM="#f59c76"),
  TissueState = c(Normal="#d0e2c0", HCM="#fded95", IZ='#F0F0F0',
                  BZ='#E3998A', RZ='#AA0A1E', FZ='#9E619A'),
  Sampling_Tissue = c(S='#FF9556', LV='#B9E4E1')
)
anno_colors <- c(Normal="#d0e2c0", HCM="#fded95", AMI="#49c2d9", ICM="#f59c76",
                 IZ='#F0F0F0', BZ='#E3998A', RZ='#AA0A1E', FZ='#9E619A',S='#FF9556', LV='#B9E4E1')
annotation_df$Sample <- factor(annotation_df$Sample, levels = sample_rank)
p_anno <- ggplot() +
  geom_tile(data = annotation_df,
            aes(x = Sample, y = 1, fill = SampleType),
            height = 0.95) +
  geom_tile(data = annotation_df,
            aes(x = Sample, y = 2, fill = TissueState),
            height = 0.95) +
  geom_tile(data = annotation_df,
            aes(x = Sample, y = 3, fill = Sampling_Tissue),
            height = 0.95) +
  scale_fill_manual(values = anno_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")) +
  ylim(0.5, 3.5)

p_final <- p_anno / p_bar +  plot_layout(heights = c(1, 5))

pdf(file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', '2-ctniche', '3-ctniche占比变化', 'sample-ctniche.pdf'),
    width = 10,height = 4.5)
print(p_final)
dev.off()