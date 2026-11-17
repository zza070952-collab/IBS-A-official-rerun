# Frozen spatial result visualisation only. No inference/model-fitting calls.
# Public entry points return grobs and never open a graphics device or export.
suppressPackageStartupMessages({library(grid);library(data.table)})

sf_source <- function(ctx) if(!is.null(ctx$source_dir)) ctx$source_dir else file.path(ctx$out,'02_Main_Figures/source_data/spatial')
sf_read <- function(ctx,name) data.table::fread(file.path(sf_source(ctx),name),na.strings=c('NA','NaN'),data.table=FALSE,check.names=FALSE)
sf_font <- function(ctx) if(is.null(ctx$font_family)) 'Arial' else ctx$font_family
sf_col <- function(ctx) list(spatial='#B78B48',dark='#27333A',muted='#737B82',line='#D9DEDF',purple='#8064A2',patient='#168B8A')
sf_text <- function(ctx,label,x=.02,y=.98,size=8,bold=FALSE,just=c('left','top'),col='#27333A',rot=0) textGrob(label,x=x,y=y,just=just,rot=rot,gp=gpar(fontfamily=sf_font(ctx),fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
sf_place <- function(g,x,y,w,h) grobTree(g,vp=viewport(x=x,y=y,width=w,height=h,just=c('left','bottom'),clip='off'))
sf_panel <- function(ctx,g,title,letter=NULL) {
 if(is.null(letter)) return(grobTree(sf_text(ctx,title,.01,.99,8.5,TRUE),sf_place(g,0,0,1,.89)))
 grobTree(sf_text(ctx,letter,.01,.99,11,TRUE),sf_text(ctx,title,if(nchar(letter)>1).09 else .055,.99,9,TRUE),sf_place(g,0,0,1,.89))
}
sf_empty <- function(ctx,text) grobTree(rectGrob(gp=gpar(fill='white',col=NA)),sf_text(ctx,text,.5,.5,8,just='centre'))
sf_domain_palette <- function(n) rep(c('#4E79A7','#59A14F','#E6AB48','#B07AA1','#76B7B2','#E15759','#9C755F','#BAB0AC'),length.out=n)
sf_color_values <- function(values,maximum=NULL) {
 if(is.null(maximum)) maximum<-max(values[is.finite(values)],na.rm=TRUE)
 if(!is.finite(maximum)||maximum<=0) maximum<-1
 pal<-grDevices::colorRampPalette(c('#F4F0E8','#DEBF88','#A7782D','#573B16'))(101)
 cols<-rep('#E3E5E5',length(values));ok<-is.finite(values)
 cols[ok]<-pal[1+pmin(100,pmax(0,round(100*values[ok]/maximum)))];cols
}
sf_image <- function(ctx,donor,domains=NULL,expression=NULL,maximum=NULL,histology=FALSE) {
 summary<-sf_read(ctx,'DONOR_SUMMARY.tsv');r<-summary[summary$donor==donor,,drop=FALSE]
 stopifnot(nrow(r)==1)
 image<-png::readPNG(file.path(sf_source(ctx),r$image_filename))
 hp<-dim(image)[1];wp<-dim(image)[2]
 content<-list(rectGrob(gp=gpar(fill='white',col=NA)))
 if(histology) content<-c(content,list(rasterGrob(image,width=1,height=1,interpolate=FALSE)))
 if(!is.null(domains)) {
  z<-domains[domains$donor==donor,,drop=FALSE];num<-as.integer(z$domain);cc<-sf_domain_palette(max(num))[num]
  content<-c(content,list(pointsGrob(default.units='npc',x=as.numeric(z$hires_x)/wp,y=1-as.numeric(z$hires_y)/hp,pch=16,size=unit(.32,'mm'),gp=gpar(col=cc))))
 }
 if(!is.null(expression)) {
  z<-expression[expression$donor==donor,,drop=FALSE];values<-as.numeric(z$log1p_CP10000)
  content<-c(content,list(pointsGrob(default.units='npc',x=as.numeric(z$hires_x)/wp,y=1-as.numeric(z$hires_y)/hp,pch=16,size=unit(.33,'mm'),gp=gpar(col=sf_color_values(values,maximum)))))
 }
 # Calibrated image-coordinate bar only. No unverified pixel-to-micron claim.
 bar_width<-5000*as.numeric(r$tissue_hires_scalef)/wp
 content<-c(content,list(rectGrob(x=.19,y=.068,width=.28,height=.105,gp=gpar(col=NA,fill='white')),segmentsGrob(x0=.06,x1=.06+bar_width,y0=.10,y1=.10,gp=gpar(col='#1D252A',lwd=1.0)),sf_text(ctx,'5,000 px',.06,.075,8,just=c('left','top'))))
 # Null-unit layout enforces the true image aspect ratio within every panel.
 inner<-do.call(grobTree,c(content,list(vp=viewport(layout.pos.row=1,layout.pos.col=1))))
 image_grob<-grobTree(inner,vp=viewport(layout=grid.layout(1,1,widths=unit(wp,'null'),heights=unit(hp,'null'),respect=TRUE)))
 if(!is.null(domains)) {
  q<-max(as.integer(domains$domain[domains$donor==donor]));x<-.10+(seq_len(q)-.5)*.80/q
  legend_y<-rep(.028,length(x))
  legend<-grobTree(pointsGrob(default.units='npc',x=x-.020,y=legend_y,pch=15,size=unit(1.1,'mm'),gp=gpar(col=sf_domain_palette(q))),textGrob(seq_len(q),x=x+.016,y=legend_y,gp=gpar(fontfamily=sf_font(ctx),fontsize=8)))
  return(grobTree(sf_place(image_grob,0,.08,1,.92),legend))
 }
 grobTree(sf_place(image_grob,0,.08,1,.92))
}
sf_gradient_legend <- function(ctx,label,max_value,x=.55,y=.03,w=.30,h=.03,size=8) {
 colors<-sf_color_values(seq(0,max_value,length.out=100),max_value)
 grobTree(rectGrob(x=seq(x,x+w,length.out=100),y=y,width=w/100,height=h,gp=gpar(col=NA,fill=colors)),sf_text(ctx,label,x,y+h*1.6,size,just=c('left','bottom')),sf_text(ctx,'0',x,y-h,size,just=c('left','top')),sf_text(ctx,format(max_value,digits=2),x+w,y-h,size,just=c('right','top')))
}
sf_marker_heatmap <- function(ctx,fontsize=8) {
 d<-sf_read(ctx,'FIXED_MARKER_NATIVE_EFFECTS.tsv');donors<-c('P1','P5','P8','P9','P10','P11','P12','P13')
 markers<-c('EPCAM','VIL1','ALPI','APOA4','DEFA5','DEFA6','MUC2','PTPRC','CD74','COL1A1','DCN','PECAM1','ACTA2','TAGLN')
 keys<-unique(d[c('donor','domain')]);keys<-keys[order(match(keys$donor,donors),as.integer(keys$domain)),];keys$key<-paste(keys$donor,keys$domain,sep=':')
 d$key<-paste(d$donor,d$domain,sep=':');x<-.16+(match(d$key,keys$key)-.5)*.82/nrow(keys);y<-.15+(length(markers)-match(d$marker,markers)+.5)*.76/length(markers)
 vals<-as.numeric(d$native_summary_logFC);cap<-max(1,min(4,max(abs(vals[is.finite(vals)]),na.rm=TRUE)))
 pal<-grDevices::colorRampPalette(c('#24578A','#FFFFFF','#A7782D'))(201)
 color<-rep('#DFE2E3',length(vals));ok<-is.finite(vals)
 color[ok]<-pal[1+pmax(0,pmin(200,round(100+100*vals[ok]/cap)))]
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),rectGrob(x=x,y=y,width=.82/nrow(keys),height=.76/length(markers),gp=gpar(col='white',lwd=.4,fill=color)),sf_text(ctx,'Original scran summary.logFC (log2 expression); grey = unavailable',.16,.99,fontsize))
 for(i in seq_along(markers))gs<-c(gs,list(sf_text(ctx,markers[i],.15,.15+(length(markers)-i+.5)*.76/length(markers),fontsize,just=c('right','centre'))))
 for(donor in donors) {
  ids<-which(keys$donor==donor);cx<-.16+(mean(ids)-.5)*.82/nrow(keys)
  gs<-c(gs,list(sf_text(ctx,donor,cx,.10,fontsize,just=c('centre','top'))))
  if(min(ids)>1) gs<-c(gs,list(segmentsGrob(x0=.16+(min(ids)-1)*.82/nrow(keys),x1=.16+(min(ids)-1)*.82/nrow(keys),y0=.15,y1=.91,gp=gpar(col='white',lwd=1.5))))
 }
 lx<-seq(.63,.92,length.out=201)
 gs<-c(gs,list(rectGrob(x=lx,y=.025,width=.29/201,height=.025,gp=gpar(fill=pal,col=NA)),sf_text(ctx,paste0('-',format(cap,digits=2)),.62,.025,8,just=c('right','centre')),sf_text(ctx,'0',.775,.025,8,just='centre'),sf_text(ctx,paste0('+',format(cap,digits=2)),.93,.025,8,just=c('left','centre'))))
 do.call(grobTree,gs)
}
sf_candidate_dot <- function(ctx) {
 d<-sf_read(ctx,'A6_DONOR_DETECTION_NATIVE.tsv');genes<-c('NCAM1','PCLO','SKIC2','C4B','VWA7','RBM33');donors<-c('P1','P5','P8','P9','P10','P11','P12','P13')
 values<-as.numeric(d$mean_log_expression);fraction<-as.numeric(d$detection_fraction);x<-.15+(match(d$donor,donors)-.5)*.62/8;y<-.14+(6-match(d$canonical_symbol,genes)+.5)*.78/6
 radius<-3.6*sqrt(pmax(fraction,0))
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),circleGrob(x=x,y=y,r=unit(radius,'mm'),gp=gpar(col=ifelse(fraction>0,'#7C6A50',NA),lwd=.4,fill=ifelse(fraction>0,sf_color_values(values),NA))))
 for(i in seq_along(genes))gs<-c(gs,list(sf_text(ctx,genes[i],.14,.14+(6-i+.5)*.78/6,8,just=c('right','centre'))))
 for(i in seq_along(donors))gs<-c(gs,list(sf_text(ctx,donors[i],.15+(i-.5)*.62/8,.10,8,just=c('centre','top'))))
 gs<-c(gs,list(sf_text(ctx,'Positive spots',.81,.91,8),sf_text(ctx,'1%   10%   30%',.80,.57,8)))
 for(i in seq_along(c(.01,.10,.30)))gs<-c(gs,list(circleGrob(x=.825+(i-1)*.07,y=.72,r=unit(3.6*sqrt(c(.01,.10,.30)[i]),'mm'),gp=gpar(col='#7C6A50',lwd=.4,fill='#D9C195'))))
 gs<-c(gs,list(sf_gradient_legend(ctx,'Mean ln(1+CP10k)',max(values,na.rm=TRUE),.81,.28,.16,.045,8)))
 do.call(grobTree,gs)
}
build_figure4 <- function(ctx) {
 donors<-if(is.null(ctx$representative_donors))c('P8','P9') else ctx$representative_donors
 stopifnot(length(donors)==2)
 dom<-sf_read(ctx,'ALL_DOMAIN_SPOTS.tsv.gz');exp<-sf_read(ctx,'A6_SPOT_EXPRESSION.tsv.gz');gg<-c('NCAM1','SKIC2','VWA7');maxexp<-max(as.numeric(exp$log1p_CP10000[exp$gene_symbol%in%gg&exp$donor%in%donors]),na.rm=TRUE)
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),sf_text(ctx,'A',.01,.995,11,TRUE),sf_text(ctx,'Tissue images and frozen spatial domains',.055,.995,9,TRUE))
 for(i in seq_along(donors)) {
  d<-donors[i];x<-.01+(i-1)*.50
  gs<-c(gs,list(sf_place(sf_panel(ctx,sf_image(ctx,d,histology=TRUE),paste(d,'histology')),x,.755,.235,.205),sf_place(sf_panel(ctx,sf_image(ctx,d,domains=dom),paste(d,'domains')),x+.25,.755,.235,.205)))
 }
 gs<-c(gs,list(sf_place(sf_panel(ctx,sf_marker_heatmap(ctx),'Marker programs from existing domain statistics','B'),.01,.47,.98,.275),sf_place(sf_panel(ctx,sf_candidate_dot(ctx),'Six transcripts across all eight reference donors','C'),.01,.295,.98,.17),sf_text(ctx,'D',.01,.285,11,TRUE),sf_text(ctx,'Actual spot expression in the same selected donors',.055,.285,9,TRUE)))
 for(i in seq_along(donors))for(j in seq_along(gg)) {
  d<-donors[i];gene<-gg[j];sub<-exp[exp$gene_symbol==gene,,drop=FALSE]
  gs<-c(gs,list(sf_place(sf_panel(ctx,sf_image(ctx,d,expression=sub,maximum=maxexp),paste(d,gene)),.01+(j-1)*.325,.050+(2-i)*.110,.315,.100)))
 }
 gs<-c(gs,list(sf_gradient_legend(ctx,'Shared expression scale: ln(1+CP10k)',maxexp,.56,.023,.36,.006,8)))
 do.call(grobTree,gs)
}

sf_counts <- function(ctx) {
 d<-sf_read(ctx,'DONOR_SUMMARY.tsv');d<-d[match(c('P1','P5','P8','P9','P10','P11','P12','P13'),d$donor),]
 n<-nrow(d);xx<-.10+(seq_len(n)-.5)*.86/n;maximum<-max(as.numeric(d$global_tested))*1.18
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)))
 for(i in seq_len(n)) {
  a<-as.numeric(d$global_tested[i]);b<-as.numeric(d$global_significant[i]);w<-.055
  gs<-c(gs,list(rectGrob(x=xx[i],y=.15+a/maximum*.68/2,width=w,height=a/maximum*.68,gp=gpar(col=NA,fill='#D8DCDD')),rectGrob(x=xx[i],y=.15+b/maximum*.68/2,width=w,height=b/maximum*.68,gp=gpar(col=NA,fill='#B78B48')),sf_text(ctx,paste0(b,'/',a),xx[i],.18+a/maximum*.68,8,just=c('centre','bottom'),rot=55),sf_text(ctx,d$donor[i],xx[i],.11,8,just=c('centre','top'))))
 }
 gs<-c(gs,list(sf_text(ctx,'BH-FDR <0.05 / tested LR pairs',.08,.98,8),segmentsGrob(x0=.07,x1=.98,y0=.15,y1=.15,gp=gpar(col='#7E8589',lwd=.6))))
 do.call(grobTree,gs)
}
sf_lr_matrix <- function(ctx,data,pairs,title_note='',fontsize=8) {
 donors<-c('P1','P5','P8','P9','P10','P11','P12','P13');nr<-length(pairs)
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),sf_text(ctx,title_note,.34,.995,fontsize))
 for(i in seq_along(pairs))for(j in seq_along(donors)) {
  r<-data[data$pair_id==pairs[i]&data$donor==donors[j],,drop=FALSE]
  tested<-nrow(r)==1&&r$test_status=='TESTED';v<-if(tested)as.numeric(r$padj) else NA_real_
  value<-if(is.finite(v)&&v>0)min(3,-log10(v)) else if(tested&&v==0)3 else NA_real_
  x<-.34+(j-.5)*.64/8;y<-.31+(nr-i+.5)*.59/nr
  fill<-if(!tested)'#E2E4E4' else if(is.finite(v)&&v<.05)sf_color_values(value,3) else '#FAF8F3'
  gs<-c(gs,list(rectGrob(x=x,y=y,width=.64/8,height=.59/nr,gp=gpar(fill=fill,col='white',lwd=.45))))
  if(!tested)gs<-c(gs,list(segmentsGrob(x0=x-.64/8*.25,x1=x+.64/8*.25,y0=y-.59/nr*.25,y1=y+.59/nr*.25,gp=gpar(col='#92999C',lwd=.6))))
  else if(v<.05)gs<-c(gs,list(pointsGrob(default.units='npc',x=x,y=y,pch=if(is.finite(as.numeric(r$pvalue))&&as.numeric(r$pvalue)==0)17 else 16,size=unit(1.0,'mm'),gp=gpar(col=if(value>=2)'white' else '#292825'))))
 }
 for(i in seq_along(pairs))gs<-c(gs,list(sf_text(ctx,gsub('_',' | ',pairs[i],fixed=TRUE),.33,.31+(nr-i+.5)*.59/nr,fontsize,just=c('right','centre'))))
 for(j in seq_along(donors))gs<-c(gs,list(sf_text(ctx,donors[j],.34+(j-.5)*.64/8,.255,fontsize,just=c('centre','top'))))
 gs<-c(gs,list(sf_text(ctx,'Pale: tested, not significant; /: untested',.34,unit(4.7,'mm'),8,just=c('left','bottom')),sf_text(ctx,'Dot: BH support; triangle: empirical P=0',.34,unit(1.3,'mm'),8,just=c('left','bottom'))))
 lx<-seq(.06,.28,length.out=100)
 gs<-c(gs,list(rectGrob(x=lx,y=.10,width=.22/100,height=.015,gp=gpar(col=NA,fill=sf_color_values(seq(0,3,length.out=100),3))),sf_text(ctx,'-log10(q), cap 3',.06,.135,8,just=c('left','bottom')),textGrob(c('0','1.5','3'),x=c(.06,.17,.28),y=rep(.015,3),just=c('centre','bottom'),gp=gpar(fontfamily=sf_font(ctx),fontsize=8))))
 do.call(grobTree,gs)
}
sf_local_map <- function(ctx) {
 d<-sf_read(ctx,'REPRESENTATIVE_ACTUAL_LOCAL.tsv');stopifnot(length(unique(d$donor))==1,length(unique(d$pair_id))==1)
 donor<-unique(d$donor);pair<-unique(d$pair_id);sm<-sf_read(ctx,'DONOR_SUMMARY.tsv');s<-sm[sm$donor==donor,]
 image<-png::readPNG(file.path(sf_source(ctx),s$image_filename));hp<-dim(image)[1];wp<-dim(image)[2]
 x<-as.numeric(d$pixel_x)*s$tissue_hires_scalef/wp;y<-1-as.numeric(d$pixel_y)*s$tissue_hires_scalef/hp
 yes<-tolower(as.character(d$significant))=='true'
 g<-grobTree(rectGrob(gp=gpar(col=NA,fill='white')),pointsGrob(default.units='npc',x=x,y=y,pch=16,size=unit(.38,'mm'),gp=gpar(col=ifelse(yes,'#8F5821','#DCDDDB'))),vp=viewport(layout.pos.row=1,layout.pos.col=1))
 body<-grobTree(g,vp=viewport(layout=grid.layout(1,1,widths=unit(wp,'null'),heights=unit(hp,'null'),respect=TRUE)))
 grobTree(sf_text(ctx,paste(donor,pair),.02,.98,8),sf_place(body,0,.11,1,.76),sf_text(ctx,paste0(sum(yes),' / ',nrow(d),' spots: native local BH <0.05'),.02,.08,8))
}
build_figure5 <- function(ctx) {
 top<-sf_read(ctx,'GLOBAL_TOP10_MATRIX.tsv');rank<-sf_read(ctx,'GLOBAL_PAIR_PLOT_RANKING.tsv');pairs<-head(rank$pair_id,10)
 a6<-sf_read(ctx,'A6_EXACT_PAIR_DETAIL_NATIVE.tsv');a6<-a6[a6$canonical_symbol=='NCAM1',,drop=FALSE]
 d<-data.frame(donor=a6$donor,pair_id=a6$pair_id,test_status=ifelse(a6$candidate_global_test_status=='TESTED','TESTED','UNTESTED'),pvalue=as.numeric(a6$global_pvalue),padj=as.numeric(a6$global_padj))
 nc<-c('NCAM1_NCAM1','NCAM1_NCAM2','NCAM1_FGFR1','NCAM1_L1CAM')
 grobTree(rectGrob(gp=gpar(fill='white',col=NA)),sf_place(sf_panel(ctx,sf_counts(ctx),'Global test denominators','A'),.01,.735,.49,.255),sf_place(sf_panel(ctx,sf_local_map(ctx),'An actually tested local pattern','C'),.52,.735,.47,.255),sf_place(sf_panel(ctx,sf_lr_matrix(ctx,top,pairs,'Existing permutation BH-FDR; colour capped at -log10(q)=3'),'Recurrent reference LR patterns','B'),.01,.285,.98,.44),sf_place(sf_panel(ctx,sf_lr_matrix(ctx,d,nc,'Strict NCAM1-NCAM1 global support: P10 and P13 only'),'NCAM1 exact-pair scope across all eight donors','D'),.01,.015,.98,.255))
}

sf_ari_ess <- function(ctx) {
 sm<-sf_read(ctx,'DONOR_SUMMARY.tsv');ch<-sf_read(ctx,'ALL_CHAIN_DIAGNOSTICS_NATIVE.tsv');donors<-c('P1','P5','P8','P9','P10','P11','P12','P13')
 ess_values<-as.numeric(ch$effective_size_pseudologlik)
 stopifnot(length(ess_values)==16,all(is.finite(ess_values)),all(ess_values>0),all(ch$stage%in%c('formal','stability')))
 ess_lo<-floor(min(log10(ess_values)));ess_hi<-ceiling(max(log10(ess_values)))
 if(ess_hi==ess_lo)ess_hi<-ess_lo+1
 ess_x<-function(value).66+.28*(log10(value)-ess_lo)/(ess_hi-ess_lo)
 gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),sf_text(ctx,'Independent-seed ARI and native pseudologlik ESS',.02,.99,9,TRUE))
 gs<-c(gs,list(sf_text(ctx,'ARI',.17,.945,8,TRUE),sf_text(ctx,'ESS: log10 x-axis; labels show actual values',.66,.945,8,TRUE)))
 for(i in seq_along(donors)) {
  donor<-donors[i];y<-.16+(8-i+.5)*.74/8;r<-sm[sm$donor==donor,];ari<-as.numeric(r$independent_seed_ARI)
  gs<-c(gs,list(sf_text(ctx,donor,.12,y,8,just=c('right','centre')),pointsGrob(default.units='npc',x=.17+.32*ari,y=y,pch=16,size=unit(1.6,'mm'),gp=gpar(col='#B78B48')),sf_text(ctx,format(ari,digits=4),.51,y,8,just=c('left','centre'))))
  cc<-ch[ch$donor==donor,]
  for(k in seq_len(nrow(cc))) {
   value<-as.numeric(cc$effective_size_pseudologlik[k]);xx<-ess_x(value);yy<-y+(if(cc$stage[k]=='formal').014 else -.014)
   color<-if(cc$stage[k]=='formal')'#24578A' else '#737B82'
   label<-format(signif(value,4),scientific=FALSE,trim=TRUE)
   tx<-if(xx>.85)xx-.008 else xx+.008;just<-if(xx>.85)c('right','centre') else c('left','centre')
   gs<-c(gs,list(pointsGrob(default.units='npc',x=xx,y=yy,pch=if(cc$stage[k]=='formal')16 else 1,size=unit(1.5,'mm'),gp=gpar(col=color)),sf_text(ctx,label,tx,yy,8,just=just,col=color)))
  }
 }
 ax<-c(.17,.33,.49);ay<-rep(.145,length(ax))
 ex<-ess_x(10^(ess_lo:ess_hi));ey<-rep(.145,length(ex))
 gs<-c(gs,list(segmentsGrob(x0=.17,x1=.49,y0=.145,y1=.145,gp=gpar(col='#A5ABAD',lwd=.6)),segmentsGrob(x0=ax,x1=ax,y0=ay,y1=ay-.009,gp=gpar(col='#A5ABAD',lwd=.6)),textGrob(c('0','0.5','1'),x=ax,y=rep(.12,length(ax)),just=c('centre','top'),gp=gpar(fontfamily=sf_font(ctx),fontsize=8)),segmentsGrob(x0=.66,x1=.94,y0=.145,y1=.145,gp=gpar(col='#A5ABAD',lwd=.6)),segmentsGrob(x0=ex,x1=ex,y0=ey,y1=ey-.009,gp=gpar(col='#A5ABAD',lwd=.6)),textGrob(paste0('10^',ess_lo:ess_hi),x=ex,y=rep(.12,length(ex)),just=c('centre','top'),gp=gpar(fontfamily=sf_font(ctx),fontsize=8)),sf_text(ctx,'Formal: solid  /  stability: open',.64,.055,8),sf_text(ctx,'Labels: 4 significant digits; full values in table',.64,.025,8)))
 do.call(grobTree,gs)
}
build_spatial_supp <- function(ctx) {
 donors<-c('P1','P5','P8','P9','P10','P11','P12','P13');dom<-sf_read(ctx,'ALL_DOMAIN_SPOTS.tsv.gz');result<-list()
 for(page in 1:2) {
  dd<-donors[((page-1)*4+1):(page*4)];gs<-list(rectGrob(gp=gpar(fill='white',col=NA)),sf_text(ctx,paste('Supplementary Figure S4',page,'All reference donor histology and frozen domains'),.02,.995,10,TRUE))
  for(i in seq_along(dd)) {
   y<-.015+(4-i)*.237
   gs<-c(gs,list(sf_place(sf_panel(ctx,sf_image(ctx,dd[i],histology=TRUE),paste(dd[i],'histology')),.01,y,.48,.225),sf_place(sf_panel(ctx,sf_image(ctx,dd[i],domains=dom),paste(dd[i],'domains')),.51,y,.48,.225)))
  }
  result[[paste0('S4_',page)]]<-do.call(grobTree,gs)
 }
 result[['S5']]<-sf_ari_ess(ctx)
 result[['S6']]<-sf_panel(ctx,sf_marker_heatmap(ctx,8),'Fixed marker panel across all 44 original domains; post hoc interpretation','S6')
 global<-sf_read(ctx,'ALL_GLOBAL_LR_NATIVE.tsv.gz');elig<-sf_read(ctx,'ALL_LR_ELIGIBILITY_NATIVE.tsv.gz');rank<-sf_read(ctx,'GLOBAL_PAIR_PLOT_RANKING.tsv');pairs<-head(rank$pair_id,30);matrix<-list()
 for(pair in pairs)for(donor in donors) {
  r<-global[global$pair_id==pair&global$donor==donor,,drop=FALSE]
  matrix[[length(matrix)+1]]<-data.frame(donor=donor,pair_id=pair,test_status=if(nrow(r))'TESTED' else 'UNTESTED',pvalue=if(nrow(r))as.numeric(r$pvalue[1]) else NA_real_,padj=if(nrow(r))as.numeric(r$padj[1]) else NA_real_)
 }
 result[['S7']]<-sf_panel(ctx,sf_lr_matrix(ctx,do.call(rbind,matrix),pairs,'All native values/statuses; descriptive top30 by donor support',8),'Additional reference LR patterns','S7')
 result
}
