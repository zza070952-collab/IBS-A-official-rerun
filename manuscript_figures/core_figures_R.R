# Original manuscript graphics from frozen native source tables; no inference.
suppressPackageStartupMessages({library(grid);library(ggplot2);library(gridExtra);library(data.table);library(jsonlite)})
cf_path<-function(ctx,n)file.path(ctx$out,'02_Main_Figures/source_data/core',n)
cf_read<-function(ctx,n)as.data.frame(data.table::fread(cf_path(ctx,n),na.strings=c('NA','NaN'),check.names=FALSE))
cf_cols<-c(HV='#747B83',IBSD='#168B8A',IBSC='#B76D9A')
cf_contrasts<-c('Duo_IBSD_vs_HV','Duo_IBSC_vs_HV','Jej_IBSD_vs_HV','Jej_IBSC_vs_HV')
cf_short<-c('Duo IBS-D','Duo IBS-C','Jej IBS-D','Jej IBS-C')
cf_theme<-function(ctx)theme_classic(base_size=8,base_family=ctx$font_family)+theme(axis.line=element_line(linewidth=.32,colour='#27333A'),axis.ticks=element_line(linewidth=.32),axis.ticks.length=unit(1.2,'mm'),axis.text=element_text(size=8,colour='#27333A'),axis.title=element_text(size=8),legend.title=element_text(size=8),legend.text=element_text(size=8),plot.title=element_text(size=9,face='bold',hjust=0),plot.tag=element_text(size=11,face='bold'),plot.margin=margin(6,5,5,5),strip.text=element_text(size=8,face='bold'),strip.background=element_blank(),legend.position='bottom',legend.margin=margin(0,0,0,0))
cf_tag<-function(p,letter,title) p+labs(tag=letter,title=title)
cf_place<-function(g,x,y,w,h)grobTree(g,vp=viewport(x=x,y=y,width=w,height=h,just=c('left','bottom'),clip='off'))
cf_label<-function(ctx,text,x,y,size=8,bold=FALSE,col='#27333A',just='centre')textGrob(text,x=x,y=y,just=just,gp=gpar(fontfamily=ctx$font_family,fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))

build_figure1<-function(ctx){
 r<-cf_read(ctx,'FIGURE1_RESOURCES.tsv');value<-function(q)format(as.numeric(r$value[r$quantity==q]),big.mark=',',scientific=FALSE,trim=TRUE)
 groups<-cf_read(ctx,'BULK_GROUP_COUNTS_NATIVE.tsv');dots<-cf_read(ctx,'SPATIAL_SPOT_FLOW_NATIVE.tsv')
 primary_stage<-dots$stage[grepl('primary',dots$stage,ignore.case=TRUE)]
 # The existing primary cohort is read from final frozen domain records.
 ds<-as.data.frame(fread(file.path(ctx$out,'02_Main_Figures/source_data/spatial/DONOR_SUMMARY.tsv')))
 genes<-cf_read(ctx,'A6_ORIGINAL_RANK.tsv')$gene_symbol
 wet<-as.data.frame(fread(file.path(ctx$out,'source_data/wetlab/QPCR_SAMPLE_DERIVED_VALUES.tsv')))
 n_batch<-length(unique(wet$Batch_ID))
 if(n_batch==0)n_batch<-length(unique(wet$Biological_Batch))
 # Column names are checked against the author-derived data, not guessed values.
 if(n_batch==0){b<-grep('batch',names(wet),ignore.case=TRUE,value=TRUE);stopifnot(length(b)>0);n_batch<-length(unique(wet[[b[1]]]))}
 colors<-c('#24578A','#168B8A','#B78B48');fill<-c('#EEF4F9','#EAF5F3','#FAF4E9');xs<-c(.18,.50,.82)
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)))
 rect<-function(x,y,w,h,fc,lc='#CCD4D8')roundrectGrob(x=x,y=y,width=w,height=h,r=unit(2,'mm'),gp=gpar(fill=fc,col=lc,lwd=.65))
 arrow<-function(x0,y0,x1,y1)segmentsGrob(x0=x0,y0=y0,x1=x1,y1=y1,gp=gpar(col='#72818A',lwd=.8),arrow=grid::arrow(length=unit(1.4,'mm'),type='closed'))
 titles<-c('Human genetic association','Patient intestinal RNA-seq','Adult spatial reference')
 details<-c(paste0('GCST90016564\n',value('cases_max'),' cases\n',value('controls_max'),' controls'),paste0('GSE166869\n',value('main_tissue_libraries'),' tissue libraries\nDuodenum and jejunum'),paste0(nrow(ds),' non-IBS surgical donors\n',format(sum(ds$n_spots),big.mark=','),' retained spots\nProximal small intestine'))
 for(i in 1:3){
  x<-xs[i];gs<-c(gs,list(rect(x,.806,.29,.315,fill[i],colors[i]),cf_label(ctx,LETTERS[i],x-.125,.945,11,TRUE,colors[i],just='left'),cf_label(ctx,titles[i],x,.915,8.5,TRUE,colors[i]),cf_label(ctx,details[i],x,.746,8.5)))
  if(i==1){t<-seq(0,1,length.out=80);yy<-.825+.060*t;xx<-.04*sin(2*pi*t);gs<-c(gs,list(linesGrob(x+xx,yy,gp=gpar(col=colors[i],lwd=1.3)),linesGrob(x-xx,yy,gp=gpar(col=colors[i],lwd=1.3))));for(j in seq(4,75,by=10))gs<-c(gs,list(segmentsGrob(x0=x+xx[j],x1=x-xx[j],y0=yy[j],y1=yy[j],gp=gpar(col=colors[i],lwd=.7))))}
  if(i==2){tt<-seq(0,1,length.out=120);gs<-c(gs,list(linesGrob(x=x+.065*sin(6*pi*tt),y=.884-.064*tt,gp=gpar(col=colors[i],lwd=2.2,lineend='round'))))}
  if(i==3){a<-expand.grid(xx=0:4,yy=0:3);gs<-c(gs,list(pointsGrob(x=x+(a$xx-2)*.018+(a$yy%%2)*.009,y=.834+a$yy*.013,default.units='npc',pch=21,size=unit(1.2,'mm'),gp=gpar(col=colors[i],fill='#E8D6B7',lwd=.6))))}
  gs<-c(gs,list(arrow(x,.64,x,.603),rect(x,.487,.29,.223,'white',colors[i])))
 }
 subheads<-c('GWAS QC\nMAGMA gene-level analysis','Count QC\nTissue-stratified DESeq2','Author-logic preprocessing\nSpatial reference analysis')
 tags<-c('GRCh37 / European LD\nGene-body association','GSEA / fixed-QC sensitivity\nIBS-D and IBS-C vs HV','BayesSpace domains\nSpatialDM global and local')
 for(i in 1:3)gs<-c(gs,list(cf_label(ctx,subheads[i],xs[i],.535,8.7,TRUE,colors[i]),cf_label(ctx,tags[i],xs[i],.442,8)))
 for(x in xs)gs<-c(gs,list(segmentsGrob(x0=x,x1=x,y0=.375,y1=.340,gp=gpar(col='#72818A',lwd=.8))))
 gs<-c(gs,list(segmentsGrob(x0=xs[1],x1=xs[3],y0=.340,y1=.340,gp=gpar(col='#72818A',lwd=.8)),arrow(.5,.340,.5,.31),rect(.5,.247,.93,.117,'#F4F6F7'),cf_label(ctx,'Evidence integration and candidate prioritization',.5,.281,9.2,TRUE),cf_label(ctx,'Unified HGNC identities | Genetic association + patient expression + spatial context',.5,.252,8)))
 for(i in 1:6)gs<-c(gs,list(rect(.116+(i-1)*.153,.216,.135,.030,'white','#D1D8DC'),cf_label(ctx,genes[i],.116+(i-1)*.153,.216,8)))
 gs<-c(gs,list(arrow(.5,.187,.5,.165),rect(.5,.099,.80,.118,'#F1EDF6','#8064A2'),cf_label(ctx,'Targeted transcript-level follow-up',.5,.135,9,TRUE,'#71568F'),cf_label(ctx,paste0('SH-SY5Y | 2 NCAM1 siRNAs | 20 nM | 48 h | n = ',n_batch,' biological batches/group'),.5,.101,8),cf_label(ctx,'RT-qPCR: NCAM1 with RPLP0 and HPRT1 reference transcripts',.5,.073,8),cf_label(ctx,'Spatial reference context; not an IBS case-control spatial comparison',.5,.022,8)))
 do.call(grobTree,gs)
}

build_figure2<-function(ctx){
 s<-cf_read(ctx,'FIG2A_GWAS_POINTS.tsv.gz');m<-cf_read(ctx,'FIG2C_MAGMA_POINTS.tsv');q<-cf_read(ctx,'FIG2B_GWAS_QQ.tsv');c<-cf_read(ctx,'GENOME_AXIS_CENTERS.tsv');d<-cf_read(ctx,'FIG2D_GENE_DISPLAY.tsv');a6<-cf_read(ctx,'A6_ORIGINAL_RANK.tsv');facts<-fromJSON(cf_path(ctx,'CORE_FACTS.json'))
 th<-cf_theme(ctx);s$chrcol<-factor(as.integer(s$CHR)%%2)
 p1<-ggplot(s,aes(plot_x_Mb,plot_y_minus_log10_P,colour=chrcol))+geom_point(size=.45,stroke=0)+scale_colour_manual(values=c('#8E969C','#24578A'),guide='none')+geom_hline(yintercept=-log10(5e-8),linetype=2,linewidth=.35,colour='#B78B48')+scale_x_continuous(breaks=c$center_Mb,labels=c$CHR,guide=guide_axis(n.dodge=2),expand=expansion(mult=.01))+labs(x='Chromosome (GRCh37)',y=expression(-log[10](italic(P))))+th+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
 p1<-cf_tag(p1,'A','IBS genome-wide association')
 q$x<--log10(q$expected_P);q$y<--log10(q$observed_P)
 p2<-ggplot(q,aes(x,y))+geom_abline(slope=1,intercept=0,linewidth=.35,colour='#BCC4C9')+geom_point(size=.6,colour='#24578A')+labs(x=expression(Expected~-log[10](italic(P))),y=expression(Observed~-log[10](italic(P))))+annotate('text',x=0,y=max(q$y)*.98,label=paste0('lambda GC = ',format(facts$GWAS_lambda_GC,digits=4)),hjust=0,vjust=1,size=8/2.845,family=ctx$font_family)+th
 p2<-cf_tag(p2,'B','GWAS quantiles')
 m$chrcol<-factor(as.integer(m$CHR)%%2);m$is_A6<-m$hgnc_id%in%a6$identity_key
 lower<-min(m$plot_x_Mb[m$CHR==6])+(28477797-min(m$START[m$CHR==6]))/1e6
 # Exact interval position uses the same chromosome offset as the plotted rows.
 offset6<-m$plot_x_Mb[which(m$CHR==6)[1]]-(m$START[which(m$CHR==6)[1]]+m$STOP[which(m$CHR==6)[1]])/2/1e6
 interval<-data.frame(xmin=offset6+28477797/1e6,xmax=offset6+33448354/1e6)
 threshold<-min(m$P[m$BH_FDR<.05]) # used only for descriptive annotation below, not a recomputed cutoff
 p3<-ggplot(m,aes(plot_x_Mb,plot_y_minus_log10_P))+geom_rect(data=interval,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill='#E7D4AA',alpha=.5)+geom_point(aes(colour=chrcol),size=.4,stroke=0)+scale_colour_manual(values=c('#9AA0A6','#526775'),guide='none')+geom_point(data=m[m$BH_FDR<.05,],colour='#B78B48',size=.8)+geom_point(data=m[m$is_A6,],shape=21,size=1.6,stroke=.5,fill='#8064A2',colour='#46344F')+geom_hline(yintercept=-log10(.05/nrow(m)),linetype=2,colour='#24578A',linewidth=.35)+scale_x_continuous(breaks=c$center_Mb,labels=c$CHR,guide=guide_axis(n.dodge=2),expand=expansion(mult=.01))+labs(x='Chromosome (GRCh37)',y=expression(-log[10](italic(P))))+th+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
 p3<-p3+annotate('text',x=offset6+31,y=max(m$plot_y_minus_log10_P)*.98,label='MHC',angle=90,vjust=-.3,hjust=1,size=8/2.845,family=ctx$font_family,colour='#9A7537')
 p3<-cf_tag(p3,'C',paste0('MAGMA: ',format(nrow(m),big.mark=','),' tested genes'))
 d$gene_label<-ifelse(d$gene_symbol==''|is.na(d$gene_symbol),d$reference_symbol,d$gene_symbol);d$gene_label<-factor(d$gene_label,levels=rev(unique(d$gene_label)));d$A6<-ifelse(d$is_A6,'Frozen A candidate','Other top-12 association')
 p4<-ggplot(d,aes(-log10(P),gene_label,colour=A6))+geom_segment(aes(x=0,xend=-log10(P),yend=gene_label),linewidth=.35,colour='#D4DADD')+geom_point(size=2)+scale_colour_manual(values=c('Frozen A candidate'='#8064A2','Other top-12 association'='#B78B48'),name=NULL)+labs(x=expression(-log[10](italic(P))),y=NULL)+th+theme(legend.position='bottom',legend.key.size=unit(2,'mm'),legend.text=element_text(size=8))
 p4<-cf_tag(p4,'D','Gene associations and six candidates')
 arrangeGrob(p1,p2,p3,p4,ncol=2,widths=c(1.15,1),heights=c(1,1.1))
}

cf_pca<-function(ctx,tissue,tag=NULL){
 d<-cf_read(ctx,paste0(tissue,'_PCA_SAMPLE_QC.tsv'));v<-cf_read(ctx,paste0(tissue,'_PCA_VARIANCE.tsv'))
 xx<-100*v$variance_fraction[v$PC=='PC1'];yy<-100*v$variance_fraction[v$PC=='PC2']
 d$group<-factor(d$group,levels=c('HV','IBSD','IBSC'));p<-ggplot(d,aes(PC1,PC2,colour=group,shape=sex))+geom_point(size=1.8,stroke=.3)+scale_colour_manual(values=cf_cols,labels=c('HV','IBS-D','IBS-C'),name=NULL)+scale_shape_manual(values=c(F=16,M=17,'NA'=15),na.value=15,name='Sex')+labs(x=sprintf('PC1 (%.1f%%)',xx),y=sprintf('PC2 (%.1f%%)',yy),title=paste0(tissue,': n = ',nrow(d),' libraries'))+cf_theme(ctx)+theme(legend.key.size=unit(2,'mm'),legend.box='horizontal')
 if(!is.null(tag))p<-p+labs(tag=tag);p
}
cf_volcano<-function(ctx,contrast,cohort='main',tag=NULL,title=NULL){
 d<-cf_read(ctx,paste0('VOLCANO_',cohort,'_',contrast,'.tsv.gz'));full_n<-nrow(d)
 d$finite<-is.finite(d$log2FoldChange)&is.finite(d$padj)&d$padj>0
 a6<-cf_read(ctx,'A6_ORIGINAL_RANK.tsv');gr<-a6[[paste0(contrast,'_gene_id')]]
 plot<-d[d$finite,];plot$signif<-plot$padj<.05;plot$A6<-plot$gene_id%in%gr
 p<-ggplot(plot,aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=signif),size=.45,stroke=0,alpha=.7)+scale_colour_manual(values=c('FALSE'='#BDC3C7','TRUE'='#168B8A'),guide='none')+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.35,colour='#949EA5')+geom_vline(xintercept=0,linewidth=.3,colour='#C8CED2')+geom_point(data=plot[plot$A6,],shape=21,size=2,stroke=.6,colour='#634976',fill='white')+labs(x='log2 fold change vs HV',y=expression(-log[10](BH-FDR)))+cf_theme(ctx)
 labs<-plot[plot$A6,];labs$label<-a6$gene_symbol[match(labs$gene_id,gr)]
 p<-p+ggrepel::geom_text_repel(data=labs,aes(label=label),size=8/2.845,family=ctx$font_family,colour='#634976',min.segment.length=0,box.padding=.2,max.overlaps=Inf,seed=20261108,segment.size=.3)
 if(is.null(title))title<-paste0(cf_short[match(contrast,cf_contrasts)],': ',sum(plot$signif),' FDR < 0.05')
 p<-p+labs(title=title,subtitle=paste0('Native FDR unavailable for ',sum(is.na(d$padj)),' / ',full_n,' features'))+theme(plot.subtitle=element_text(size=8))
 if(!is.null(tag))p<-p+labs(tag=tag);p
}
cf_gsea<-function(ctx,tag=NULL){
 d<-cf_read(ctx,'FIG3D_GSEA_DISPLAY.tsv');d$contrast<-factor(d$contrast,levels=cf_contrasts,labels=cf_short)
 d$pathway_label<-gsub('^(HALLMARK_|REACTOME_)','',d$pathway);d$pathway_label<-gsub('_',' ',d$pathway_label);d$pathway_label<-vapply(d$pathway_label,function(x)paste(strwrap(x,width=30),collapse='\n'),'')
 d$pathway_label<-factor(d$pathway_label,levels=rev(unique(d$pathway_label)));d$display_q<-pmin(-log10(d$padj),12);d$support<-d$padj<.05
 cap<-ceiling(max(abs(d$NES),na.rm=TRUE))
 p<-ggplot(d,aes(contrast,pathway_label))+geom_tile(fill='#F5F6F6',colour='white')+geom_point(aes(size=display_q,fill=NES,shape=support),colour='#525B61',stroke=.3)+scale_shape_manual(values=c('TRUE'=21,'FALSE'=22),guide='none')+scale_size_continuous(range=c(.5,3),breaks=c(3,8,12),name='-log10 FDR\n(cap 12)')+scale_fill_gradient2(low='#24578A',mid='white',high='#B78B48',midpoint=0,limits=c(-cap,cap),breaks=c(-cap,0,cap),name='NES')+guides(fill=guide_colourbar(order=1,barwidth=unit(20,'mm'),barheight=unit(2,'mm'),title.position='top'),size=guide_legend(order=2,title.position='top',nrow=1))+labs(x=NULL,y=NULL,title='Cross-contrast pathway patterns')+cf_theme(ctx)+theme(axis.text.x=element_text(angle=45,hjust=1),axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',legend.box='horizontal',legend.key.size=unit(2,'mm'),legend.text=element_text(size=8),legend.title=element_text(size=8),plot.title.position='plot')
 if(!is.null(tag))p<-p+labs(tag=tag);p
}
cf_qc<-function(ctx,tag=NULL){
 d<-cf_read(ctx,'DEG_COUNTS_ALL_QC.tsv');d$contrast<-factor(d$contrast,levels=cf_contrasts,labels=cf_short);d$cohort<-factor(d$cohort,levels=c('main','qc_low_750k','qc_high_1250k'),labels=c('Primary (63 libraries)','750k (66 libraries)','1250k (59 libraries)'))
 p<-ggplot(d,aes(contrast,cohort))+geom_tile(aes(fill=log10(1+FDR05)),colour='white',linewidth=1)+geom_text(aes(label=format(FDR05,big.mark=',',trim=TRUE)),family=ctx$font_family,size=8/2.845)+scale_fill_gradient(low='#F2F6F5',high='#8ABCB7',guide='none')+labs(x=NULL,y=NULL,title='FDR-significant feature counts under fixed QC')+cf_theme(ctx)+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='none')
 if(!is.null(tag))p<-p+labs(tag=tag);p
}
build_figure3<-function(ctx){
 pa<-arrangeGrob(cf_pca(ctx,'Duo','A')+guides(shape='none'),cf_pca(ctx,'Jej',' ')+guides(colour='none'),ncol=2)
 counts<-cf_read(ctx,'DEG_COUNTS_ALL_QC.tsv');counts<-counts[counts$cohort=='main',];dd<-rbind(transform(counts,direction='Up',count=up),transform(counts,direction='Down',count=down));dd$contrast<-factor(dd$contrast,levels=cf_contrasts,labels=cf_short)
 pb<-ggplot(dd,aes(contrast,count,fill=direction))+geom_col(position='stack',width=.65)+scale_fill_manual(values=c(Up='#168B8A',Down='#24578A'),name=NULL)+labs(x=NULL,y='Features',title='Primary differential expression',tag='B')+cf_theme(ctx)+theme(axis.text.x=element_text(angle=45,hjust=1),legend.key.size=unit(2,'mm'))
 pc<-cf_volcano(ctx,'Jej_IBSD_vs_HV','main','C','Jejunum IBS-D vs HV')
 pd<-cf_gsea(ctx,'D');pe<-cf_qc(ctx,'E')
 # Each existing-data panel is distinct; no workflow or mechanism diagrams.
 top<-arrangeGrob(pa,pb,ncol=2,widths=c(1.8,1));middle<-arrangeGrob(pc,pd,ncol=2,widths=c(1,1.1))
 arrangeGrob(top,middle,pe,ncol=1,heights=c(.30,.50,.20))
}

build_core_supp<-function(ctx){
 m<-cf_read(ctx,'FIG2C_MAGMA_POINTS.tsv');q<-data.frame(expected=-log10((seq_len(nrow(m))-.5)/nrow(m)),observed=-log10(sort(m$P)))
 pp<-ggplot(q,aes(expected,observed))+geom_abline(slope=1,intercept=0,colour='#AAB4BC',linewidth=.35)+geom_point(size=.6,colour='#24578A')+labs(x='Expected -log10 P',y='Observed -log10 P',title='All tested MAGMA gene quantiles',tag='A')+cf_theme(ctx)
 qc<-fromJSON(cf_path(ctx,'CORE_FACTS.json'))$GWAS_counts;dd<-data.frame(stage=c('Raw rows','Unique valid autosomal','Exact EUR LD match'),n=as.numeric(c(qc$raw_rows,qc$valid_unique_autosomal_snps,qc$eligible_exact_ld_matches)))
 pp2<-ggplot(dd,aes(n/1e6,factor(stage,levels=rev(stage))))+geom_col(fill='#24578A',width=.6)+geom_text(aes(label=format(n,big.mark=',',trim=TRUE)),hjust=1.04,size=8/2.845,family=ctx$font_family,colour='white')+labs(x='Variants (millions)',y=NULL,title='Recorded GWAS QC counts',tag='B')+cf_theme(ctx)
 S1<-arrangeGrob(pp,pp2,ncol=1)
 vols<-lapply(seq_along(cf_contrasts),function(i)cf_volcano(ctx,cf_contrasts[i],tag=LETTERS[i]))
 S2_1<-arrangeGrob(cf_pca(ctx,'Duo','A'),cf_pca(ctx,'Jej','B'),ncol=1)
 S2_2<-arrangeGrob(grobs=vols,ncol=2)
 S3_1<-arrangeGrob(grobs=lapply(seq_along(cf_contrasts),function(i)cf_volcano(ctx,cf_contrasts[i],'qc_low_750k',LETTERS[i],paste0(cf_short[i],' | QC 750k'))),ncol=2)
 S3_2<-arrangeGrob(grobs=lapply(seq_along(cf_contrasts),function(i)cf_volcano(ctx,cf_contrasts[i],'qc_high_1250k',LETTERS[i],paste0(cf_short[i],' | QC 1250k'))),ncol=2)
 result<-list(S1=S1,S2_1=S2_1,S2_2=S2_2,S3_1=S3_1,S3_2=S3_2)
 for(tissue in c('Duo','Jej')){
  path<-cf_path(ctx,paste0(tissue,'_VST_SAMPLE_CORRELATION.tsv'))
  if(file.exists(path)){
   d<-cf_read(ctx,paste0(tissue,'_VST_SAMPLE_CORRELATION.tsv'));ids<-unique(d$row_library)
   d$row_library<-factor(d$row_library,levels=rev(ids));d$column_library<-factor(d$column_library,levels=ids)
   p<-ggplot(d,aes(column_library,row_library,fill=Pearson_r))+geom_tile()+scale_fill_gradient(low='#EFF5F5',high='#24578A',name='Pearson r')+coord_equal()+labs(x=NULL,y=NULL,title=paste(tissue,'sample correlation: archived blind VST'))+cf_theme(ctx)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8),axis.text.y=element_text(size=8),axis.line=element_blank(),axis.ticks=element_blank(),plot.title.position='plot')
   result[[if(tissue=='Duo')'S2_3' else 'S2_4']]<-ggplotGrob(p)
  }
 }
 all<-cf_read(ctx,'GSEA_ALL_COHORTS.tsv.gz');chosen<-cf_read(ctx,'FIG3D_GSEA_DISPLAY.tsv');keys<-unique(paste(chosen$collection,chosen$pathway,sep='::'));d<-all[paste(all$collection,all$pathway,sep='::')%in%keys,]
 d$pathway_label<-gsub('^(HALLMARK_|REACTOME_)','',d$pathway);d$pathway_label<-gsub('_',' ',d$pathway_label);d$pathway_label<-vapply(d$pathway_label,function(x)paste(strwrap(x,width=28),collapse='\n'),'')
 d$contrast<-factor(d$contrast,levels=cf_contrasts,labels=c('Duo D','Duo C','Jej D','Jej C'));d$cohort<-factor(d$cohort,levels=c('main','qc_low_750k','qc_high_1250k'),labels=c('Primary','QC 750k','QC 1250k'));d$qdisplay<-pmin(-log10(d$padj),12);cap<-ceiling(max(abs(d$NES),na.rm=TRUE))
 p<-ggplot(d,aes(contrast,pathway_label))+geom_tile(fill='#F3F4F4',colour='white')+geom_point(aes(size=qdisplay,fill=NES),shape=21,stroke=.3,colour='#505B62',na.rm=TRUE)+facet_grid(.~cohort)+scale_size_continuous(range=c(.5,3),breaks=c(3,8,12),name='-log10 FDR (cap 12)')+scale_fill_gradient2(low='#24578A',mid='white',high='#B78B48',limits=c(-cap,cap),name='NES')+labs(x=NULL,y=NULL,title='Identical primary-selected pathways across overlapping QC cohorts')+cf_theme(ctx)+theme(axis.text.x=element_text(angle=45,hjust=1),axis.line=element_blank(),axis.ticks=element_blank(),plot.title.position='plot',legend.position='bottom')
 result[['S3_3']]<-ggplotGrob(p)
 result[order(match(names(result),c('S1','S2_1','S2_2','S2_3','S2_4','S3_1','S3_2','S3_3')))]
}
