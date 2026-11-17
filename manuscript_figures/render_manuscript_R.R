suppressPackageStartupMessages({library(grid);library(jsonlite);library(ragg);library(svglite)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=2);out<-normalizePath(args[1],mustWork=TRUE);phase<-args[2]
selected<-if(length(args)>=3)strsplit(args[3],',',fixed=TRUE)[[1]] else as.character(1:6)
source(file.path(out,'code/core_figures_R.R'),encoding='UTF-8')
source(file.path(out,'code/spatial_figures_R.R'),encoding='UTF-8')
source(file.path(out,'code/build_figure6_and_s8.R'),encoding='UTF-8')
ctx<-list(out=out,root='/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105',font_family='Arial',representative_donors=c('P8','P9'),source_dir=file.path(out,'02_Main_Figures/source_data/spatial'),figure6_source_dir=file.path(out,'02_Main_Figures/source_data/figure6'))
options(warn=1)
# Use the selected R raster device for font metrics, without an implicit Rplots.pdf.
ragg::agg_capture(width=180,height=220,units='mm',res=120,background='white')
metrics_device<-grDevices::dev.cur()
for(p in c('02_Main_Figures/preview','02_Main_Figures/editable','02_Main_Figures/submission','03_Supplementary/figures/preview','03_Supplementary/figures/editable','logs'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
status<-list(stage=phase,R=R.version.string,backend='R/grid and ggplot2; existing results only',font_requested=ctx$font_family,font_match=systemfonts::match_font(ctx$font_family),biological_models_run=FALSE,figures=list(),supplementary=list())
draw<-function(g){grid.newpage();grid.draw(g)}
record_file<-function(path)list(path=path,bytes=file.info(path)$size)
save_grob<-function(g,name,width,height,main=TRUE,dpi=600){
 prefix<-if(main)'02_Main_Figures' else '03_Supplementary/figures'
 files<-list();png<-file.path(out,prefix,'preview',paste0(name,'.png'))
 ragg::agg_png(png,width=width,height=height,units='mm',res=200,background='white');tryCatch(draw(g),finally=dev.off());files$preview<-record_file(png)
 if(phase=='final'){
  pdf<-file.path(out,prefix,'editable',paste0(name,'.pdf'));grDevices::cairo_pdf(pdf,width=width/25.4,height=height/25.4,family=ctx$font_family,onefile=FALSE);tryCatch(draw(g),finally=dev.off());files$pdf<-record_file(pdf)
  svg<-file.path(out,prefix,'editable',paste0(name,'.svg'));svglite::svglite(svg,width=width/25.4,height=height/25.4,bg='white',standalone=TRUE);tryCatch(draw(g),finally=dev.off());files$svg<-record_file(svg)
  if(main){
   tif<-file.path(out,prefix,'submission',paste0(name,'.tiff'));ragg::agg_tiff(tif,width=width,height=height,units='mm',res=dpi,background='white',bitsize=8,compression='lzw');tryCatch(draw(g),finally=dev.off());files$tiff<-record_file(tif)
   qa<-file.path(out,'02_Main_Figures/preview',paste0(name,'_QA_200percent.png'));ragg::agg_png(qa,width=width,height=height,units='mm',res=400,background='white');tryCatch(draw(g),finally=dev.off());files$qa_200_percent<-record_file(qa)
  }
 }
 list(width_mm=width,height_mm=height,dpi=if(main)dpi else NA,files=files,render_completed=TRUE,visual_review='PENDING_ROOT_ACTUAL_VIEW',backend='R')
}
sizes<-list(c(180,165),c(180,175),c(180,220),c(180,220),c(180,190),c(180,180))
for(i in selected[selected!='supp']){
 n<-as.integer(i);fun<-get(paste0('build_figure',n));g<-fun(ctx);size<-sizes[[n]]
 status$figures[[paste0('Figure_',i)]]<-save_grob(g,paste0('Figure_',i),size[1],size[2],TRUE,if(n==1)1200 else 600)
 if(n==1){
  pp<-file.path(out,'02_Main_Figures/preview/Figure_1_Embedded.png')
  ragg::agg_png(pp,width=150,height=130,units='mm',res=300,background='white');tryCatch(draw(g),finally=dev.off())
  status$figures[[paste0('Figure_',i)]]$internal_review_layout<-list(path=pp,width_mm=150,height_mm=130,font_size_preserved=TRUE)
 }
 write_json(status,file.path(out,'logs',paste0('FIGURE_RENDER_',toupper(phase),'_STATUS.json')),auto_unbox=TRUE,pretty=TRUE,na='null')
 cat(paste0('RENDERED_FIGURE_',i,'\n'));flush.console()
}
 if(length(selected)==6||identical(selected,'supp')){
 supplementary<-c(build_core_supp(ctx),build_spatial_supp(ctx),list(S8=build_s8(ctx)))
 for(name in names(supplementary)){
  status$supplementary[[name]]<-save_grob(supplementary[[name]],name,180,if(name%in%c('S4_1','S4_2','S6','S7'))220 else 190,FALSE)
  cat(paste0('RENDERED_SUPPLEMENT_',name,'\n'));flush.console()
 }
}
write_json(status,file.path(out,'logs',paste0('FIGURE_RENDER_',toupper(phase),'_STATUS.json')),auto_unbox=TRUE,pretty=TRUE,na='null')
writeLines(capture.output(sessionInfo()),file.path(out,'logs',paste0('R_MANUSCRIPT_',phase,'_SESSIONINFO.txt')))
cat('MANUSCRIPT_GRAPHICS_STAGE_FINISHED\n')
if(metrics_device%in%grDevices::dev.list())grDevices::dev.off(metrics_device)
