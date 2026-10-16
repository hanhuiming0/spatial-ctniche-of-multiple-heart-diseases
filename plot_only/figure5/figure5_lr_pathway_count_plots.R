options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(CellChat);library(ggplot2)})
fig5<-Sys.getenv('CTNICHE_FIG5');stopifnot(dir.exists(fig5))
input<-file.path(fig5,'33-5C_shared929_LR_CellChatDB_CellPhoneDB')
out<-file.path(fig5,'34-5C_shared929_LR_pathway_counts');dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd<-function(f)read.delim(f,check.names=FALSE,quote='',comment.char='')
wr<-function(x,f)write.table(x,file.path(out,f),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
data(CellChatDB.human)
inter<-CellChatDB.human$interaction;cx<-CellChatDB.human$complex
canon<-function(x){x<-sub('^COMPLEX:','',toupper(trimws(x)));z<-unlist(strsplit(x,'[_+;|]'));paste(sort(unique(z[!is.na(z)&nzchar(z)])),collapse=';')}
resolve<-function(x){if(x%in%rownames(cx)){v<-as.character(unlist(cx[x,,drop=FALSE]));v<-v[!is.na(v)&nzchar(v)];canon(paste(v,collapse=';'))}else canon(x)}
annot<-data.frame(reference_interaction=rownames(inter),ligand=inter$ligand,receptor=inter$receptor,pathway=inter$pathway_name)
annot$ligand_key<-vapply(annot$ligand,resolve,character(1));annot$receptor_key<-vapply(annot$receptor,resolve,character(1))
annot$key<-paste(annot$ligand_key,annot$receptor_key,sep=' -> ')
annot<-annot[!is.na(annot$pathway)&nzchar(annot$pathway),]
wr(annot,'34_01_CellChat_pathway_reference.tsv')
anno_pairs<-unique(annot[c('key','pathway')])
allstats<-list();summaries<-list()
for(db in c('CellChatDB','CellPhoneDB')){
 d<-rd(file.path(input,db,paste0(db,'_02_pairs_at_least_one_end.tsv')))
 d$pair_id<-paste(d$source_genesymbol,d$target_genesymbol,sep=' -> ')
 d$key<-paste(vapply(d$source_genesymbol,canon,character(1)),vapply(d$target_genesymbol,canon,character(1)),sep=' -> ')
 stopifnot(!anyDuplicated(d$key))
 z<-merge(d,anno_pairs,by='key',all.x=TRUE,sort=FALSE)
 z$pathway[is.na(z$pathway)]<-'Unmapped'
 z$pathway_source<-paste0('CellChat ',packageVersion('CellChat'),' human DB; exact directed subunit-set match')
 z<-z[!duplicated(z[c('key','pathway')]),]
 wr(z,paste0(db,'_01_pair_pathway_annotations.tsv'))
 wr(z[z$pathway=='Unmapped',],paste0(db,'_02_unmapped_pairs.tsv'))
 subsets<-list(At_least_one_end=rep(TRUE,nrow(z)),Both_ends_any_subunit=z$ligand_has_929&z$receptor_has_929,Both_ends_complete=z$both_ends_complete_in_929)
 ss<-lapply(names(subsets),function(k){zz<-z[subsets[[k]],];t<-as.data.frame(table(zz$pathway));names(t)<-c('pathway','n_pairs');t$subset<-k;t$database<-db;t})
 ss<-do.call(rbind,ss);ss$pathway<-as.character(ss$pathway);wr(ss,paste0(db,'_03_pathway_pair_counts.tsv'))
 summaries[[db]]<-ss
 allstats[[db]]<-do.call(rbind,lapply(names(subsets),function(k){zz<-z[subsets[[k]],];data.frame(database=db,subset=k,total_pairs=length(unique(zz$key)),mapped_pairs=length(unique(zz$key[zz$pathway!='Unmapped'])),unmapped_pairs=length(unique(zz$key[zz$pathway=='Unmapped'])),n_named_pathways=length(setdiff(unique(zz$pathway),'Unmapped')))}))
 # Figure 1: all candidate pairs. Top 25 pathways plus Unmapped; complete data in TSV.
 a<-ss[ss$subset=='At_least_one_end',];a<-a[order(-a$n_pairs,a$pathway),];keep<-head(a$pathway[a$pathway!='Unmapped'],25)
 a<-a[a$pathway%in%c(keep,'Unmapped'),];a$pathway<-factor(a$pathway,levels=rev(a$pathway));a$mapped<-ifelse(a$pathway=='Unmapped','Unmapped','Mapped')
 p<-ggplot(a,aes(pathway,n_pairs,fill=mapped))+geom_col(width=.72)+geom_text(aes(label=n_pairs),hjust=-.15,size=3.3)+coord_flip()+scale_fill_manual(values=c(Mapped=if(db=='CellChatDB')'#457B9D' else '#C47B51',Unmapped='grey70'),guide='none')+scale_y_continuous(expand=expansion(mult=c(0,.15)))+theme_bw(base_size=11)+labs(x=NULL,y='Number of unique ligand-receptor pairs',title=paste(db,'- candidate pairs by signaling pathway'),subtitle='At least one end hits shared929 | top 25 mapped pathways + Unmapped',caption='Pathway labels: CellChat human database; exact directed complex matching.')
 ggsave(file.path(out,paste0(db,'_Fig1_pathways_at_least_one_end.pdf')),p,width=9,height=max(5,2.5+.24*nrow(a)))
 # Figure 2: both ends; mutually exclusive partial vs complete coverage.
 b<-z[z$ligand_has_929&z$receptor_has_929,];b$coverage<-ifelse(b$both_ends_complete_in_929,'All subunits in shared929','Partial subunit coverage')
 t<-as.data.frame(table(b$pathway,b$coverage));names(t)<-c('pathway','coverage','n_pairs');t<-t[t$n_pairs>0,]
 tot<-aggregate(n_pairs~pathway,t,sum);tot<-tot[order(tot$n_pairs,tot$pathway),];t$pathway<-factor(t$pathway,levels=tot$pathway)
 p<-ggplot(t,aes(pathway,n_pairs,fill=coverage))+geom_col(width=.72)+coord_flip()+scale_fill_manual(values=c('All subunits in shared929'='#457B9D','Partial subunit coverage'='#E9AE83'))+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(x=NULL,y='Number of unique ligand-receptor pairs',fill=NULL,title=paste(db,'- pathways with both ends in shared929'),subtitle='Complete and partial subunit coverage shown separately',caption='Pathway labels: CellChat human database; Unmapped = no exact reference match.')
 ggsave(file.path(out,paste0(db,'_Fig2_pathways_both_ends.pdf')),p,width=9,height=max(5,2.5+.28*nrow(tot)))
}
wr(do.call(rbind,allstats),'34_02_mapping_and_pathway_summary.tsv')
wr(do.call(rbind,summaries),'34_03_all_pathway_pair_counts.tsv')
writeLines(c('Input interactions are exactly the existing shared929 CellChatDB and CellPhoneDB selections (folder33).',
 'Pathway reference is installed CellChat human DB; see version in sessionInfo. CellPhoneDB results use transferred CellChat pathway labels, not native CellPhoneDB annotations.',
 'Ligand and receptor direction is retained; complete complex subunit sets must match. No ligand-only or reverse-direction fallback.',
 'One pair may map to multiple pathways: pathway counts are not necessarily additive to the number of unique pairs.',
 'Unmapped does not mean biologically unclassified; it means no exact match to this reference version.',
 'Figure1 shows top25 named pathways plus Unmapped. All pathways are in TSV. Figure2 shows every pathway among both-end hits.',
 'The counts describe catalog membership, NOT pathway activation, significance, or inferred communication strength.',
 'Existing figures and source files were not modified.'),file.path(out,'README.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(do.call(rbind,allstats))