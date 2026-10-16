ST2_ROOT <- Sys.getenv("ST2_ROOT", unset = "/path/to/ST2")

setwd(file.path(ST2_ROOT, "CellChat", "Data"))
library(CellChat)

Disease.Color <- c(Normal = "#d0e2c0", HCM = "#fded95", AMI = "#49c2d9", ICM = "#f59c76")

cellchat.Sample_Type <- readRDS("cellchat.Sample_Type_2026.1.7.rds")
cellchat.Sample_Type
cellchat_merged <- mergeCellChat(cellchat.Sample_Type[c(1, 2)], add.names = names(cellchat.Sample_Type)[c(1, 2)])
gg1 <- netVisual_heatmap(cellchat_merged, measure = "count")
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
pdf(file.path(ST2_ROOT, "CellChat", "Figure2", "AMI vs Normal netVisual_heatmap.pdf"), width = 7.5, height = 4.2)
gg1 + gg2
dev.off()

cellchat_merged <- mergeCellChat(cellchat.Sample_Type[c(1, 3)], add.names = names(cellchat.Sample_Type)[c(1, 3)])
gg1 <- netVisual_heatmap(cellchat_merged, measure = "count")
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
pdf(file.path(ST2_ROOT, "CellChat", "Figure2", "ICM vs Normal netVisual_heatmap.pdf"), width = 7.5, height = 4.2)
gg1 + gg2
dev.off()

cellchat_merged <- mergeCellChat(cellchat.Sample_Type[c(1, 4)], add.names = names(cellchat.Sample_Type)[c(1, 4)])
gg1 <- netVisual_heatmap(cellchat_merged, measure = "count")
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
pdf(file.path(ST2_ROOT, "CellChat", "Figure2", "HCM vs Normal netVisual_heatmap.pdf"), width = 7.5, height = 4.2)
gg1 + gg2
dev.off()

CellChat_Sample_Type_Comparison <- readRDS("CellChat_Sample_Type_Comparison_2026.1.7.rds")
CellChat_Sample_Type_Comparison
