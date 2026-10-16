EXTERNAL_DATA_ROOT <- Sys.getenv('EXTERNAL_DATA_ROOT', unset = '.')
LEGACY_RESULTS_ROOT <- Sys.getenv('LEGACY_RESULTS_ROOT', unset = '.')

PROJECT_ROOT <- Sys.getenv('PROJECT_ROOT', unset = '.')

# ACH007

#- 单样本测试结果,ACH007,'FN1','SDC2'-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
sample <- 'ACH007'
sample <- 'ACH0016'
ACH007 <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
ACH007 <- NormalizeData(ACH007) %>% FindVariableFeatures() %>% ScaleData()

#- 单基因表达----------
pdf(file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'ACH0016.gene.FeaturePlot_2026.6.4.pdf'),width = 15, height = 6.3)
SpatialFeaturePlot(ACH007,features = c('ITGA4','SPP1'),pt.size.factor = 2) +
  plot_layout(
    ncol = 2, nrow = 1,
    guides = "collect",
    heights = c(1,1,1),  # 每行高度相等
    widths = c(1,1,1)    # 每列宽度相等
  )
dev.off()

SpatialFeaturePlot(ACH007,features = c('VEGFA'),pt.size.factor = 2)
SpatialFeaturePlot(ACH007,features = c('PSAP'),pt.size.factor = 2)
SpatialFeaturePlot(ACH007,features = c('LRP1'),pt.size.factor = 2)

pdf(file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'ACH0016.PSAP-LRP1.pdf'),width = 15, height = 6.3)
SpatialFeaturePlot(ACH007,features = c('PSAP','LRP1'),pt.size.factor = 2) +
  plot_layout(
    ncol = 2, nrow = 1,
    guides = "collect",
    heights = c(1,1,1),  # 每行高度相等
    widths = c(1,1,1)    # 每列宽度相等
  )
dev.off()

pdf(file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'ACH0016.FN1-PLAUR.pdf'),width = 15, height = 6.3)
SpatialFeaturePlot(ACH007,features = c('FN1','PLAUR'),pt.size.factor = 2) +
  plot_layout(
    ncol = 2, nrow = 1,
    guides = "collect",
    heights = c(1,1,1),  # 每行高度相等
    widths = c(1,1,1)    # 每列宽度相等
  )
dev.off()

PSAP-LRP1#- ggplot两基因表达----------
knn = 7
pt.size=1
alpha.min=0.1
max.cut=0.95# 共定位强度的截断值
# 获得位置信息
spatial_coord <- data.frame(GetTissueCoordinates(ACH007))
colnames(spatial_coord) <- c('imagerow','imagecol','barcodeID')
location = spatial_coord[,c('imagerow','imagecol')]
topn=floor(0.2*dim(location)[1])# 高表达位点阈值：前20%的位点定义为高表达
# 所有基因
GenePairs = c('FN1','SDC2')
expr = GetAssayData(ACH007, assay = "Spatial", layer = "data")
expr = t(expr[GenePairs,])
ncell <- dim(expr)[1]

# 计算每个spot的7个最近邻位点的索引（nn2返回nn.idx=近邻索引矩阵）
nnmatrix <- RANN::nn2(location,k=knn)$nn.idx# 高效计算 k 近邻
countsum <- Matrix::colSums(expr)
####normalize
expr <- Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))

# 计算每个spot的近邻位点中目标基因的最大表达量（排除自身，取2:knn）
gene1 <- expr[,GenePairs[1]]
gene2 <- expr[,GenePairs[2]]
pairsexp <- rbind(gene1,gene2)
neighexp <- apply(nnmatrix,1,function(x){apply(pairsexp[,x[2:knn]],1,max)})
# 计算共定位强度（两个基因的表达量×近邻最大表达量的最大值）
LRadd <- pmax(pairsexp[1,]*neighexp[2,],pairsexp[2,]*neighexp[1,])
LRadd_max <- quantile(LRadd,probs=max.cut)
LRadd[LRadd>LRadd_max] <- LRadd_max

if(sum(gene1>0)>topn){n1<-order(gene1,sample(ncell,ncell),decreasing=T)[1:topn]}else{n1<-which(gene1>0)} # Top20高表达基因spot
if(sum(gene2>0)>topn){n2<-order(gene2,sample(ncell,ncell),decreasing=T)[1:topn]}else{n2<-which(gene2>0)} #
expcol <- rep(0,ncell)
expcol[n1] <- 1
expcol[n2] <- 2
expcol[intersect(n1,n2)] <- 3

tmp <- data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))
tmpLRadd <- data.frame(x=location[,1],y=location[,2],Genes=LRadd)
alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min) + alpha.min

location = spatial_coord[,c('imagerow','imagecol')]
p1 <- ggplot(tmp,aes(x = x,y = y,col = Exp)) +
  geom_point(size=pt.size) +
  scale_color_manual(values=c("gray","#c7e2c0","#7eaac9","#ed7c6f"),
                     labels=c("Both_low","FN1_high","SDC2_high","Both_high")) +
  ggtitle(paste0(GenePairs,collapse="_")) + xlab("") + ylab("") +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

p2 <- ggplot(tmpLRadd,aes(x = x,y = y,col = Genes)) +
  geom_point(size=pt.size,alpha=alpha)+
  scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="white",mid="white",high="#ce5d15")+
  xlab("")+
  ylab("")+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  labs(color="colocalization") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

pdf(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '2-单细胞', '5-基因', 'ACH007.colocalization.ggplot2.pdf'),width = 10,height = 4)
p1 + p2 & scale_x_reverse() & coord_flip()
dev.off()

#- 单样本测试结果,ACH004-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
sample <- 'ACH004'
ACH004 <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
ACH004 <- NormalizeData(ACH004) %>% FindVariableFeatures() %>% ScaleData()

#- 单基因表达----------
pdf(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '2-单细胞', '5-基因', 'ACH004.gene.FeaturePlot.pdf'),width = 15, height = 6.3)
SpatialFeaturePlot(ACH004,features = c('FN1','SDC2'),pt.size.factor = 2) +
  plot_layout(
    ncol = 2, nrow = 1,
    guides = "collect",
    heights = c(1,1,1),  # 每行高度相等
    widths = c(1,1,1)    # 每列宽度相等
  )
dev.off()

#- ggplot两基因表达----------
knn = 7
pt.size=1
alpha.min=0.1
max.cut=0.95# 共定位强度的截断值
# 获得位置信息
spatial_coord <- data.frame(GetTissueCoordinates(ACH004))
colnames(spatial_coord) <- c('imagerow','imagecol','barcodeID')
location = spatial_coord[,c('imagerow','imagecol')]
topn=floor(0.2*dim(location)[1])# 高表达位点阈值：前20%的位点定义为高表达
# 所有基因
GenePairs = c('FN1','SDC2')
expr = GetAssayData(ACH004, assay = "Spatial", layer = "data")
expr = t(expr[GenePairs,])
ncell <- dim(expr)[1]

# 计算每个spot的7个最近邻位点的索引（nn2返回nn.idx=近邻索引矩阵）
nnmatrix <- RANN::nn2(location,k=knn)$nn.idx# 高效计算 k 近邻
countsum <- Matrix::colSums(expr)
####normalize
expr <- Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))

# 计算每个spot的近邻位点中目标基因的最大表达量（排除自身，取2:knn）
gene1 <- expr[,GenePairs[1]]
gene2 <- expr[,GenePairs[2]]
pairsexp <- rbind(gene1,gene2)
neighexp <- apply(nnmatrix,1,function(x){apply(pairsexp[,x[2:knn]],1,max)})
# 计算共定位强度（两个基因的表达量×近邻最大表达量的最大值）
LRadd <- pmax(pairsexp[1,]*neighexp[2,],pairsexp[2,]*neighexp[1,])
LRadd_max <- quantile(LRadd,probs=max.cut)
LRadd[LRadd>LRadd_max] <- LRadd_max

if(sum(gene1>0)>topn){n1<-order(gene1,sample(ncell,ncell),decreasing=T)[1:topn]}else{n1<-which(gene1>0)} # Top20高表达基因spot
if(sum(gene2>0)>topn){n2<-order(gene2,sample(ncell,ncell),decreasing=T)[1:topn]}else{n2<-which(gene2>0)} #
expcol.normal <- rep(0,ncell)
expcol.normal[n1] <- 1
expcol.normal[n2] <- 2
expcol.normal[intersect(n1,n2)] <- 3

tmp.normal <- data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol.normal))
tmpLRadd.normal<- data.frame(x=location[,1],y=location[,2],Genes=LRadd)
alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min) + alpha.min

location = spatial_coord[,c('imagerow','imagecol')]
p1 <- ggplot(tmp.normal,aes(x = x,y = y,col = Exp)) +
  geom_point(size=pt.size) +
  scale_color_manual(values=c("gray","#c7e2c0","#7eaac9","#ed7c6f"),
                     labels=c("Both_low","FN1_high","SDC2_high","Both_high")) +
  ggtitle(paste0(GenePairs,collapse="_")) + xlab("") + ylab("") +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

p2 <- ggplot(tmpLRadd.normal,aes(x = x,y = y,col = Genes)) +
  geom_point(size=pt.size,alpha=alpha)+
  scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="white",mid="white",high="#ce5d15")+
  xlab("")+
  ylab("")+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  labs(color="colocalization") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

pdf(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '2-单细胞', '5-基因', 'ACH004.colocalization.ggplot2.pdf'),width = 10,height = 4)
p1 + p2 & scale_x_reverse() & coord_flip()
dev.off()

#- for循环样本测试结果-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
Samples <- dir()

knn = 7
pt.size=1
alpha.min=0.1
max.cut=0.95
GenePairs = c('FN1','SDC2')

for(sample in Samples){
  sample.seurat <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
  sample.seurat <- NormalizeData(sample.seurat) %>% FindVariableFeatures() %>% ScaleData()

  p <- SpatialFeaturePlot(sample.seurat,features = c('FN1','SDC2'),pt.size.factor = 2) +
    plot_layout(
      ncol = 2, nrow = 1,
      guides = "collect",
      heights = c(1,1,1),  # 每行高度相等
      widths = c(1,1,1)    # 每列宽度相等
    )
  pdf(paste0(paste0(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '2-单细胞', '5-基因', 'FeaturePlot'), '/'),sample,'.gene.FeaturePlot.pdf'),width = 15, height = 6.3)
  print(p)
  dev.off()

  # 获得位置信息
  spatial_coord <- data.frame(GetTissueCoordinates(sample.seurat))
  colnames(spatial_coord) <- c('imagerow','imagecol','barcodeID')
  location = spatial_coord[,c('imagerow','imagecol')]
  topn=floor(0.2*dim(location)[1])# 高表达位点阈值：前20%的位点定义为高表达
  # 所有基因

  expr = GetAssayData(sample.seurat, assay = "Spatial", layer = "data")
  expr = t(expr[GenePairs,])
  ncell <- dim(expr)[1]

  # 计算每个spot的7个最近邻位点的索引（nn2返回nn.idx=近邻索引矩阵）
  nnmatrix <- RANN::nn2(location,k=knn)$nn.idx# 高效计算 k 近邻
  countsum <- Matrix::colSums(expr)
  ####normalize
  expr <- Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))

  # 计算每个spot的近邻位点中目标基因的最大表达量（排除自身，取2:knn）
  gene1 <- expr[,GenePairs[1]]
  gene2 <- expr[,GenePairs[2]]
  pairsexp <- rbind(gene1,gene2)
  neighexp <- apply(nnmatrix,1,function(x){apply(pairsexp[,x[2:knn]],1,max)})
  # 计算共定位强度（两个基因的表达量×近邻最大表达量的最大值）
  LRadd <- pmax(pairsexp[1,]*neighexp[2,],pairsexp[2,]*neighexp[1,])
  LRadd_max <- quantile(LRadd,probs=max.cut)
  LRadd[LRadd>LRadd_max] <- LRadd_max

  if(sum(gene1>0)>topn){n1<-order(gene1,sample(ncell,ncell),decreasing=T)[1:topn]}else{n1<-which(gene1>0)} # Top20高表达基因spot
  if(sum(gene2>0)>topn){n2<-order(gene2,sample(ncell,ncell),decreasing=T)[1:topn]}else{n2<-which(gene2>0)} #
  expcol.normal <- rep(0,ncell)
  expcol.normal[n1] <- 1
  expcol.normal[n2] <- 2
  expcol.normal[intersect(n1,n2)] <- 3

  tmp.normal <- data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol.normal))
  tmpLRadd.normal<- data.frame(x=location[,1],y=location[,2],Genes=LRadd)
  alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min) + alpha.min

  location = spatial_coord[,c('imagerow','imagecol')]
  p1 <- ggplot(tmp.normal,aes(x = x,y = y,col = Exp)) +
    geom_point(size=pt.size) +
    scale_color_manual(values=c("gray","#c7e2c0","#7eaac9","#ed7c6f"),
                       labels=c("Both_low","FN1_high","SDC2_high","Both_high")) +
    ggtitle(paste0(GenePairs,collapse="_")) + xlab("") + ylab("") +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank()) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())

  p2 <- ggplot(tmpLRadd.normal,aes(x = x,y = y,col = Genes)) +
    geom_point(size=pt.size,alpha=alpha)+
    scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="white",mid="white",high="#ce5d15")+
    xlab("")+
    ylab("")+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    labs(color="colocalization") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())

  p3 <- p1 + p2 & scale_x_reverse() & coord_flip()
  pdf(paste0(paste0(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '2-单细胞', '5-基因', 'FeaturePlot'), '/'),sample,'.colocalization.ggplot2.pdf'),width = 10,height = 4)
  print(p3)
  dev.off()
}

#- 单样本测试结果,ACH0013,FN1 − (ITGAV+ITGB1)-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
sample <- 'ACH0013'
ACH0013 <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
ACH0013 <- NormalizeData(ACH0013) %>% FindVariableFeatures() %>% ScaleData()

#- 单基因表达----------
pdf(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'ACH0013.gene.FeaturePlot.pdf'),width = 15, height = 6.3)
SpatialFeaturePlot(ACH0013,features = c('FN1','ITGAV','ITGB1'),pt.size.factor = 2) +
  plot_layout(
    ncol = 3, nrow = 1,
    guides = "collect",
    heights = c(1,1,1),  # 每行高度相等
    widths = c(1,1,1)    # 每列宽度相等
  )
dev.off()

#- for循环样本测试结果-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
Samples <- dir()

GenePairs = c('FN1','ITGAV','ITGB1')

for(sample in Samples){
  sample.seurat <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
  sample.seurat <- NormalizeData(sample.seurat) %>% FindVariableFeatures() %>% ScaleData()

  p <- SpatialFeaturePlot(sample.seurat,features = c('FN1','ITGAV','ITGB1'),pt.size.factor = 2) +
    plot_layout(
      ncol = 3, nrow = 1,
      guides = "collect",
      heights = c(1,1,1),  # 每行高度相等
      widths = c(1,1,1)    # 每列宽度相等
    )
  pdf(paste0(paste0(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'FN1-ITGAV-ITGB1'), '/'),sample,'.gene.FeaturePlot.pdf'),width = 15, height = 6.3)
  print(p)
  dev.off()
}

for(sample in Samples){
  sample.seurat <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
  sample.seurat <- NormalizeData(sample.seurat) %>% FindVariableFeatures() %>% ScaleData()

  p <- SpatialFeaturePlot(sample.seurat,features = c('THBS1','THBS2','CD36'),pt.size.factor = 2) +
    plot_layout(
      ncol = 3, nrow = 1,
      guides = "collect",
      heights = c(1,1,1),  # 每行高度相等
      widths = c(1,1,1)    # 每列宽度相等
    )
  pdf(paste0(paste0(file.path(LEGACY_RESULTS_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位', 'THBS1-THBS2-CD36'), '/'),sample,'.gene.FeaturePlot.pdf'),width = 15, height = 6.3)
  print(p)
  dev.off()
}

#- MIF − CD74 PDGFC − PDGFRA LAMA2 − CD44 for循环样本测试结果-------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

knn = 7
pt.size=1
alpha.min=0.1
max.cut=0.95

GenePairs = c('MIF','CD74')
GenePairs = c('PDGFC','PDGFRA')
GenePairs = c('LAMA2','CD44')

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '2.Zenodo6580069', 'Data'))
Samples <- c('ACH006','ACH0011')

setwd(file.path(EXTERNAL_DATA_ROOT, 'HF2', 'Hs', '1.figshare_5777948', 'Data'))
Samples <- c('HCM1225D')

for(sample in Samples){
  sample.seurat <- Load10X_Spatial(data.dir = sample, assay = "Spatial", slice = sample)
  sample.seurat <- NormalizeData(sample.seurat) %>% FindVariableFeatures() %>% ScaleData()

  # 获得位置信息
  spatial_coord <- data.frame(GetTissueCoordinates(sample.seurat))
  colnames(spatial_coord) <- c('imagerow','imagecol','barcodeID')
  location = spatial_coord[,c('imagerow','imagecol')]
  topn=floor(0.2*dim(location)[1])# 高表达位点阈值：前20%的位点定义为高表达
  # 所有基因

  expr = GetAssayData(sample.seurat, assay = "Spatial", layer = "data")
  expr = t(expr[GenePairs,])
  ncell <- dim(expr)[1]

  # 计算每个spot的7个最近邻位点的索引（nn2返回nn.idx=近邻索引矩阵）
  nnmatrix <- RANN::nn2(location,k=knn)$nn.idx# 高效计算 k 近邻
  countsum <- Matrix::colSums(expr)
  ####normalize
  expr <- Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))

  # 计算每个spot的近邻位点中目标基因的最大表达量（排除自身，取2:knn）
  gene1 <- expr[,GenePairs[1]]
  gene2 <- expr[,GenePairs[2]]
  pairsexp <- rbind(gene1,gene2)
  neighexp <- apply(nnmatrix,1,function(x){apply(pairsexp[,x[2:knn]],1,max)})
  # 计算共定位强度（两个基因的表达量×近邻最大表达量的最大值）
  LRadd <- pmax(pairsexp[1,]*neighexp[2,],pairsexp[2,]*neighexp[1,])
  LRadd_max <- quantile(LRadd,probs=max.cut)
  LRadd[LRadd>LRadd_max] <- LRadd_max

  if(sum(gene1>0)>topn){n1<-order(gene1,sample(ncell,ncell),decreasing=T)[1:topn]}else{n1<-which(gene1>0)} # Top20高表达基因spot
  if(sum(gene2>0)>topn){n2<-order(gene2,sample(ncell,ncell),decreasing=T)[1:topn]}else{n2<-which(gene2>0)} #
  expcol.normal <- rep(0,ncell)
  expcol.normal[n1] <- 1
  expcol.normal[n2] <- 2
  expcol.normal[intersect(n1,n2)] <- 3

  tmp.normal <- data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol.normal))
  tmpLRadd.normal<- data.frame(x=location[,1],y=location[,2],Genes=LRadd)
  alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min) + alpha.min

  location = spatial_coord[,c('imagerow','imagecol')]
  p1 <- ggplot(tmp.normal,aes(x = x,y = y,col = Exp)) +
    geom_point(size=pt.size) +
    scale_color_manual(values=c("gray","#c7e2c0","#7eaac9","#ed7c6f"),
                       # labels=c("Both_low","MIF_high","CD74_high","Both_high")
                       # labels=c("Both_low","PDGFC_high","PDGFRA_high","Both_high")
                       labels=c("Both_low","LAMA2_high","CD44_high","Both_high")
                       ) +
    ggtitle(paste0(GenePairs,collapse="_")) + xlab("") + ylab("") +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank()) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())

  p2 <- ggplot(tmpLRadd.normal,aes(x = x,y = y,col = Genes)) +
    geom_point(size=pt.size,alpha=alpha)+
    scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="white",mid="white",high="#ce5d15")+
    xlab("")+
    ylab("")+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    labs(color="colocalization") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())

  p3 <- p1 + p2 & scale_x_reverse() & coord_flip()
  pdf(paste0(paste0(file.path(PROJECT_ROOT, '3-结果与分析', '3-空间转录组', '8-受配体互作', '3-基因共定位'), '/'),sample,'.colocalization.ggplot2.pdf'),width = 10,height = 4)
  print(p3)
  dev.off()
}
