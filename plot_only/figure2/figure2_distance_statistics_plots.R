#!/usr/bin/env Rscript
# Reviewer analysis 03: sample-level statistics for distance gradients.
# Usage:

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<4)stop("Required: ctniche8_sample_csv fibrotic_sample_csv infarct_sample_csv output_dir [threshold_root] [continuous_score_file]")
FILES<-c(ctniche8=args[1],fibrotic=args[2],infarct=args[3]);OUT<-args[4]
THRESHOLD_ROOT<-if(length(args)>=5)args[5] else "";CONTINUOUS_FILE<-if(length(args)>=6)args[6] else ""
set.seed(20260928);N_BOOT<-1000L
if(!requireNamespace("ggplot2",quietly=TRUE))stop("Missing ggplot2")
library(ggplot2);dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
wtab<-function(x,f)write.table(x,file.path(OUT,f),sep="\t",quote=FALSE,row.names=FALSE,na="NA")

read_distance<-function(file,target){
  d<-read.csv(file,check.names=FALSE);req<-c("ring","cell_type","proportion","sample");if(!all(req%in%names(d)))stop(file," lacks: ",paste(setdiff(req,names(d)),collapse=", "))
  if(!"group"%in%names(d))d$group<-if(target=="infarct")"AMI-IZ" else "Unknown"
  if(!"n_spots"%in%names(d))d$n_spots<-NA_real_
  lev<-if(target=="ctniche8")c("0","1","2","3","4","5",">5") else c("Core","radius1","radius2","radius3","radius4","radius5","Outside")
  d$distance_index<-match(as.character(d$ring),lev)-1L;d$target<-target
  d$proportion<-as.numeric(d$proportion);d$n_spots<-as.numeric(d$n_spots)
  d
}
all<-do.call(rbind,lapply(names(FILES),function(n)read_distance(FILES[[n]],n)))
wtab(all,"R03_01_all_sample_level_distance_data.tsv")

# Bin QC and sample-bootstrap 95% confidence intervals. Samples, not spots, are resampled.
key<-interaction(all$target,all$group,all$cell_type,all$ring,drop=TRUE,lex.order=TRUE)
boot_one<-function(d){
  d<-d[is.finite(d$proportion),,drop=FALSE];ns<-length(unique(d$sample))
  if(!nrow(d))return(data.frame(n_samples=0,total_spots=0,mean=NA,median=NA,ci_low=NA,ci_high=NA))
  vals<-d$proportion
  b<-if(ns>=2)replicate(N_BOOT,mean(sample(vals,length(vals),replace=TRUE),na.rm=TRUE)) else rep(mean(vals),N_BOOT)
  data.frame(n_samples=ns,total_spots=sum(d$n_spots,na.rm=TRUE),median_spots_per_sample=median(d$n_spots,na.rm=TRUE),
    mean=mean(d$proportion),median=median(d$proportion),ci_low=quantile(b,.025,na.rm=TRUE),ci_high=quantile(b,.975,na.rm=TRUE))
}
bs<-do.call(rbind,lapply(split(all,key),function(d)cbind(d[1,c("target","group","cell_type","ring","distance_index")],boot_one(d))))
rownames(bs)<-NULL;wtab(bs,"R03_02_distance_bin_bootstrap_95CI_and_QC.tsv")

# Per-sample monotonic trend, followed by sample-level group test against zero.
trend_one<-function(d){
  d<-d[is.finite(d$proportion)&is.finite(d$distance_index),];if(nrow(d)<4||length(unique(d$distance_index))<4)return(NULL)
  ct<-suppressWarnings(cor.test(d$distance_index,d$proportion,method="spearman",exact=FALSE))
  data.frame(target=d$target[1],group=d$group[1],cell_type=d$cell_type[1],sample=d$sample[1],n_bins=nrow(d),rho=unname(ct$estimate),p_value=ct$p.value)
}
tr_key<-interaction(all$target,all$group,all$cell_type,all$sample,drop=TRUE)
tr<-do.call(rbind,lapply(split(all,tr_key),trend_one));if(is.null(tr))stop("No sample trends")
tr$fdr_within_all_tests<-p.adjust(tr$p_value,"BH");wtab(tr,"R03_03_per_sample_spearman_trends.tsv")
g_key<-interaction(tr$target,tr$group,tr$cell_type,drop=TRUE)
gt<-do.call(rbind,lapply(split(tr,g_key),function(d){
  wt<-if(nrow(d)>=3)wilcox.test(d$rho,mu=0,exact=FALSE) else NULL
  data.frame(target=d$target[1],group=d$group[1],cell_type=d$cell_type[1],n_samples=nrow(d),median_rho=median(d$rho),mean_rho=mean(d$rho),
    positive_fraction=mean(d$rho>0),wilcoxon_p=if(is.null(wt))NA else wt$p.value)
}));gt$fdr<-ave(gt$wilcoxon_p,interaction(gt$target,gt$group),FUN=function(p)p.adjust(p,"BH"));wtab(gt,"R03_04_group_level_trend_tests.tsv")

# Optional mixed-effects linear trend on logit proportions.
if(requireNamespace("lme4",quietly=TRUE)){
  mix<-do.call(rbind,lapply(split(all,interaction(all$target,all$group,all$cell_type,drop=TRUE)),function(d){
    d<-d[is.finite(d$proportion)&is.finite(d$distance_index),];if(length(unique(d$sample))<3||nrow(d)<12)return(NULL)
    eps<-1e-5;d$logit_p<-qlogis(pmin(pmax(d$proportion,eps),1-eps))
    fit<-try(lme4::lmer(logit_p~distance_index+(1|sample),data=d,REML=FALSE),silent=TRUE);if(inherits(fit,"try-error"))return(NULL)
    cf<-coef(summary(fit))["distance_index",]
    data.frame(target=d$target[1],group=d$group[1],cell_type=d$cell_type[1],n_samples=length(unique(d$sample)),estimate=cf[1],se=cf[2],t_value=cf[3])
  }));if(!is.null(mix)){mix$p_value<-2*pnorm(-abs(mix$t_value));mix$fdr<-ave(mix$p_value,interaction(mix$target,mix$group),FUN=function(p)p.adjust(p,"BH"));wtab(mix,"R03_05_mixed_effect_trend_models.tsv")}
}else writeLines("Package lme4 not installed; mixed-effects model skipped.",file.path(OUT,"R03_05_MIXED_MODEL_NOT_RUN.txt"))

# Main plot contains selected fibroblast/myeloid subtypes; all subtypes remain in the tables.
priority<-c("FN1+ Fib","HIF3A+ Fib","THBS4+ Fib","FOLR2+ MP","SPP1+ MP","CPM+ MP","Monocyte")
show_types<-intersect(priority,unique(bs$cell_type));if(!length(show_types))show_types<-head(unique(bs$cell_type),8)
plotdat<-bs[bs$cell_type%in%show_types,];plotdat$ring_label<-factor(plotdat$ring,levels=c("0","Core","1","radius1","2","radius2","3","radius3","4","radius4","5","radius5",">5","Outside"))
p<-ggplot(plotdat,aes(distance_index,mean,color=cell_type,fill=cell_type,group=cell_type))+geom_ribbon(aes(ymin=ci_low,ymax=ci_high),alpha=.12,color=NA)+
  geom_line(linewidth=.75)+geom_point(size=1.5)+facet_grid(target~group,scales="free_x")+
  labs(x="Distance bin (not physical mm unless separately calibrated)",y="Mean cell-type proportion",color=NULL,fill=NULL)+theme_bw()+theme(legend.position="bottom")
ggsave(file.path(OUT,"R03_Fig1_key_subtype_distance_trends_95CI.pdf"),p,width=12,height=7)

# Optional threshold sensitivity: recursively process all sample-level distance CSVs.
if(nzchar(THRESHOLD_ROOT)&&dir.exists(THRESHOLD_ROOT)){
  fs<-list.files(THRESHOLD_ROOT,pattern="sample_level.*distance.*\\.csv$",recursive=TRUE,full.names=TRUE)
  sens<-list();k<-0L
  for(f in fs){
    d<-try(read.csv(f,check.names=FALSE),silent=TRUE);if(inherits(d,"try-error")||!all(c("ring","cell_type","proportion","sample")%in%names(d)))next
    thr<-sub(".*阈值([0-9.]+).*","\\1",f);if(identical(thr,f))thr<-NA_character_
    d$distance_index<-match(d$ring,c("Core","radius1","radius2","radius3","radius4","radius5","Outside"))-1L
    d$threshold<-thr;d$source_file<-f;if(!"group"%in%names(d))d$group<-"AMI-IZ"
    z<-do.call(rbind,lapply(split(d,interaction(d$group,d$cell_type,d$sample,drop=TRUE)),function(q){q$target<-"threshold_sensitivity";trend_one(q)}))
    if(!is.null(z)){z$threshold<-thr;z$source_file<-f;k<-k+1L;sens[[k]]<-z}
  }
  if(length(sens))wtab(do.call(rbind,sens),"R03_06_threshold_sensitivity_sample_trends.tsv")
}

# Optional continuous Fibrotic Score file: sample, cell_type, proportion, fibrotic_score.
if(nzchar(CONTINUOUS_FILE)&&file.exists(CONTINUOUS_FILE)){
  d<-read.delim(CONTINUOUS_FILE,check.names=FALSE);req<-c("sample","cell_type","proportion","fibrotic_score")
  if(!all(req%in%names(d)))stop("Continuous score file lacks: ",paste(setdiff(req,names(d)),collapse=", "))
  cr<-do.call(rbind,lapply(split(d,interaction(d$sample,d$cell_type,drop=TRUE)),function(q){
    q<-q[is.finite(q$proportion)&is.finite(q$fibrotic_score),];if(nrow(q)<20)return(NULL);ct<-cor.test(q$fibrotic_score,q$proportion,method="spearman",exact=FALSE)
    data.frame(sample=q$sample[1],cell_type=q$cell_type[1],n_spots=nrow(q),rho=unname(ct$estimate),p_value=ct$p.value)
  }));cr$fdr<-p.adjust(cr$p_value,"BH");wtab(cr,"R03_07_continuous_fibrotic_score_correlations.tsv")
}
writeLines(c("All inference treats the sample as the independent unit.","Bootstrap intervals resample samples, not spots.",
  "ctniche_8 rings are spot-spacing bins; fibrosis/infarct rings are coordinate-distance bins unless physical calibration is supplied.",
  "Missing sample-by-bin combinations are not converted to zero."),file.path(OUT,"R03_METHOD_NOTES.txt"))
message("Reviewer03 completed: ",OUT)