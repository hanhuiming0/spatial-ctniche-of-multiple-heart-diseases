options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
out<-Sys.getenv('MISTY_FIG_DIR');stopifnot(dir.exists(out))
f<-Sys.getenv('MISTY_DATA_FILE');stopifnot(file.exists(f))
d<-read.delim(f,check.names=FALSE)
kp<-data.frame(Predictor=c('SPP1MP','FN1Fib','FOLR2MP','FN1Fib','CPMMP','FN1Fib'),Target=c('FN1Fib','SPP1MP','FN1Fib','FOLR2MP','FN1Fib','CPMMP'))
kp$Pair<-paste(kp$Predictor,'\u2192',kp$Target)
groups<-c('Normal','AMI-IZ','ICM-FZ','HCM')
vl<-c(intra='Same spot',para.2='l=2',para.5='l=5')
dd<-inner_join(d,kp,by=c('Predictor','Target')) %>% filter(Group%in%groups,View%in%names(vl)) %>% mutate(Group=factor(Group,levels=groups),Pair=factor(Pair,levels=kp$Pair),View_label=factor(vl[as.character(View)],levels=unname(vl)))
stopifnot(n_distinct(dd$Pair)==6,n_distinct(dd$Group)==4)
cols<-c(Normal='#6B7280','AMI-IZ'='#D55E00','ICM-FZ'='#0072B2',HCM='#8E5AA9')
p<-ggplot(dd,aes(View_label,Importance,group=Sample_ID,colour=Group))+
 geom_line(alpha=.22,linewidth=.35)+geom_point(alpha=.38,size=.8)+
 stat_summary(aes(group=Group),fun=median,geom='line',linewidth=1.05)+
 stat_summary(aes(group=Group),fun=median,geom='point',size=2)+
 facet_grid(Pair~Group,scales='free_y')+scale_colour_manual(values=cols,guide='none')+
 labs(title='Scale-dependent fibroblast-myeloid spatial dependencies',x=NULL,y='MISTy importance')+
 theme_bw(base_size=8)+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=.25),strip.background=element_rect(fill='grey95',colour='grey65'),axis.text.x=element_text(angle=35,hjust=1),plot.title=element_text(face='bold'))
ggsave(file.path(out,'Figure4_Fib_Myeloid_revised_height75pct.pdf'),p,width=8.5,height=8.25,device=cairo_pdf)
cat('Saved vector PDF: 8.5 x 8.25 inches; raw Importance; ',nrow(dd),' observations.\n')