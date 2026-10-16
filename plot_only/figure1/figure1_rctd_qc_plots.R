#!/usr/bin/env Rscript
# Reviewer analysis 02: RCTD reference coverage, weight QC, and spatial marker concordance.
# Usage:

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<5) stop("Required: weights.rds single_cell.rds spatial_seurat.rds subtype_markers.tsv output_dir")
WEIGHT_FILE<-args[1];SC_FILE<-args[2];ST_FILE<-args[3];MARKER_FILE<-args[4];OUT<-args[5]
TOP_MARKERS<-15L;MIN_SPOTS<-30L
pkgs<-c("Matrix","SeuratObject","ggplot2")
miss<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop("Missing packages: ",paste(miss,collapse=", "))
suppressPackageStartupMessages({library(Matrix);library(SeuratObject);library(ggplot2)})
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
wtab<-function(x,f) write.table(x,file.path(OUT,f),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
canon<-function(x)tolower(gsub("[^A-Za-z0-9]","",x))
find_col<-function(d,candidates){hit<-candidates[candidates%in%names(d)];if(!length(hit))stop("Missing metadata column; tried: ",paste(candidates,collapse=", "));hit[1]}

weights<-readRDS(WEIGHT_FILE);weights<-as.matrix(weights);storage.mode(weights)<-"double"
if(is.null(rownames(weights))||is.null(colnames(weights)))stop("Weights need spot row names and subtype column names")
rs<-rowSums(weights);qc_spot<-data.frame(spot=rownames(weights),weight_sum=rs,n_negative=rowSums(weights<0),
  n_na=rowSums(!is.finite(weights)),max_weight=apply(weights,1,max,na.rm=TRUE),n_nonzero=rowSums(weights>0))
# Sample names are assigned later from the spatial object; preserve raw sums before normalization.
good<-is.finite(rs)&rs>0&qc_spot$n_negative==0&qc_spot$n_na==0
weights_norm<-weights[good,,drop=FALSE]/rs[good]
wtab(qc_spot,"R02_01_RCTD_spot_weight_QC.tsv")
qc_summary<-data.frame(n_spots=nrow(weights),n_valid=sum(good),n_invalid=sum(!good),n_negative=sum(qc_spot$n_negative>0),
  n_with_na=sum(qc_spot$n_na>0),weight_sum_min=min(rs,na.rm=TRUE),weight_sum_median=median(rs,na.rm=TRUE),weight_sum_max=max(rs,na.rm=TRUE))
wtab(qc_summary,"R02_02_RCTD_weight_QC_summary.tsv")

# Reference coverage: numbers of cells and donors supporting every subtype.
sc<-readRDS(SC_FILE);meta_sc<-sc@meta.data
sub_col<-find_col(meta_sc,c("SubType","subtype","Cell_SubType","cell_subtype"))
donor_col<-find_col(meta_sc,c("orig.ident","Sample_ID","sample","sample_id"))
disease_col<-find_col(meta_sc,c("Sample_Type","disease","Disease","group"))
ref_n<-aggregate(rep(1,nrow(meta_sc)),list(subtype=as.character(meta_sc[[sub_col]])),sum);names(ref_n)[2]<-"n_cells"
ref_donor<-aggregate(as.character(meta_sc[[donor_col]]),list(subtype=as.character(meta_sc[[sub_col]])),function(x)length(unique(x)));names(ref_donor)[2]<-"n_donors"
ref_qc<-merge(ref_n,ref_donor,by="subtype",all=TRUE)
dis_tab<-as.data.frame.matrix(table(meta_sc[[sub_col]],meta_sc[[disease_col]]));dis_tab$subtype<-rownames(dis_tab);rownames(dis_tab)<-NULL
ref_qc<-merge(ref_qc,dis_tab,by="subtype",all.x=TRUE);ref_qc<-ref_qc[order(ref_qc$n_cells),]
wtab(ref_qc,"R02_03_single_cell_reference_subtype_coverage.tsv")

# Marker selection. This is a concordance check, not fully independent validation if these markers trained RCTD.
mk<-read.delim(MARKER_FILE,check.names=FALSE)
gene_col<-find_col(mk,c("gene","Gene","features","feature"));cluster_col<-find_col(mk,c("cluster","SubType","subtype","cell_type"))
padj_col<-intersect(c("p_val_adj","p_adj","FDR","padj"),names(mk));fc_col<-intersect(c("avg_log2FC","avg_logFC","logFC"),names(mk))
if(length(padj_col))mk<-mk[is.na(mk[[padj_col[1]]])|mk[[padj_col[1]]]<0.05,,drop=FALSE]
if(length(fc_col))mk<-mk[is.na(mk[[fc_col[1]]])|mk[[fc_col[1]]]>0,,drop=FALSE]
mk$gene_symbol<-as.character(mk[[gene_col]]);mk$subtype_raw<-as.character(mk[[cluster_col]])
mk<-mk[!grepl("^(MT-|RPL|RPS)",mk$gene_symbol,ignore.case=TRUE),,drop=FALSE]
if(length(fc_col))mk<-mk[order(mk$subtype_raw,-mk[[fc_col[1]]]),] else mk<-mk[order(mk$subtype_raw),]
mk_top<-do.call(rbind,lapply(split(mk,mk$subtype_raw),function(d)head(d[!duplicated(d$gene_symbol),],TOP_MARKERS)))
marker_out_cols<-unique(na.omit(c("subtype_raw","gene_symbol",padj_col[1],fc_col[1])))
wtab(mk_top[,marker_out_cols,drop=FALSE],"R02_04_markers_used_for_concordance.tsv")

# Spatial counts and sample annotation.
st<-readRDS(ST_FILE);meta_st<-st@meta.data
st_sample_col<-find_col(meta_st,c("orig.ident","Sample_ID","sample","sample_id"))
assay_name<-if("Spatial"%in%names(st@assays))"Spatial" else names(st@assays)[1]
counts<-tryCatch(LayerData(st[[assay_name]],layer="counts"),error=function(e)GetAssayData(st,assay=assay_name,slot="counts"))
common_spots<-Reduce(intersect,list(colnames(counts),rownames(meta_st),rownames(weights_norm)))
counts<-counts[,common_spots,drop=FALSE];weights_norm<-weights_norm[common_spots,,drop=FALSE];meta_st<-meta_st[common_spots,,drop=FALSE]
spot_sample<-as.character(meta_st[[st_sample_col]])

# Map marker subtype labels to RCTD column names after punctuation-insensitive canonicalization.
map<-data.frame(marker_subtype=unique(mk_top$subtype_raw),marker_key=canon(unique(mk_top$subtype_raw)))
map$rctd_subtype<-colnames(weights_norm)[match(map$marker_key,canon(colnames(weights_norm)))]
wtab(map,"R02_05_marker_to_RCTD_subtype_mapping.tsv")
valid_map<-map[!is.na(map$rctd_subtype),,drop=FALSE]
all_marker_genes<-unique(mk_top$gene_symbol[mk_top$subtype_raw%in%valid_map$marker_subtype])
all_marker_genes<-intersect(all_marker_genes,rownames(counts))
lib<-Matrix::colSums(counts);lib[lib<=0]<-NA_real_
lognorm<-log1p(t(t(counts[all_marker_genes,,drop=FALSE])/lib)*1e4)

cors<-list();scores_long<-list();j<-0L
for(i in seq_len(nrow(valid_map))){
  ms<-valid_map$marker_subtype[i];rsu<-valid_map$rctd_subtype[i]
  gs<-intersect(mk_top$gene_symbol[mk_top$subtype_raw==ms],rownames(lognorm));if(length(gs)<3)next
  score<-Matrix::colMeans(lognorm[gs,,drop=FALSE]);abund<-weights_norm[,rsu]
  for(s in unique(spot_sample)){
    idx<-which(spot_sample==s&is.finite(score)&is.finite(abund));if(length(idx)<MIN_SPOTS)next
    ct<-suppressWarnings(cor.test(score[idx],abund[idx],method="spearman",exact=FALSE));j<-j+1L
    cors[[j]]<-data.frame(sample=s,marker_subtype=ms,rctd_subtype=rsu,n_spots=length(idx),n_markers=length(gs),rho=unname(ct$estimate),p_value=ct$p.value)
  }
  scores_long[[i]]<-data.frame(spot=common_spots,sample=spot_sample,marker_subtype=ms,rctd_subtype=rsu,marker_score=score,rctd_abundance=abund)
}
cor_tab<-do.call(rbind,cors);if(is.null(cor_tab))stop("No marker-subtype correlations could be computed")
cor_tab$fdr<-p.adjust(cor_tab$p_value,"BH");wtab(cor_tab,"R02_06_sample_level_marker_RCTD_correlations.tsv")
cor_summary<-do.call(rbind,lapply(split(cor_tab,cor_tab$rctd_subtype),function(d)data.frame(rctd_subtype=d$rctd_subtype[1],n_samples=nrow(d),
  median_rho=median(d$rho,na.rm=TRUE),positive_fraction=mean(d$rho>0,na.rm=TRUE),significant_positive_fraction=mean(d$rho>0&d$fdr<0.05,na.rm=TRUE))))
wtab(cor_summary,"R02_07_marker_RCTD_correlation_summary.tsv")

p<-ggplot(cor_tab,aes(x=reorder(rctd_subtype,rho,FUN=median),y=rho))+geom_hline(yintercept=0,lty=2,color="grey60")+
  geom_boxplot(outlier.shape=NA,fill="#D9E8F5")+geom_jitter(width=.15,size=1,alpha=.55)+coord_flip()+
  labs(x=NULL,y="Within-sample Spearman rho",title="RCTD abundance and spatial marker-score concordance")+theme_bw()
ggsave(file.path(OUT,"R02_Fig1_marker_RCTD_concordance.pdf"),p,width=7,height=max(5,.25*nrow(cor_summary)))

writeLines(c("RCTD weights were checked and row-normalized only after raw-sum QC.",
  "Marker concordance is not fully independent if the same markers were used to construct the RCTD reference.",
  "Singlet/doublet/reject fractions require the original per-sample RCTD result objects and cannot be reconstructed from the weight matrix."),
  file.path(OUT,"R02_LIMITATIONS.txt"))
saveRDS(list(weight_qc=qc_summary,reference_qc=ref_qc,marker_correlations=cor_tab,marker_summary=cor_summary),file.path(OUT,"R02_all_results.rds"))
message("Reviewer02 completed: ",OUT)