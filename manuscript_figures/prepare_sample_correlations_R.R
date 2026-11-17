# Descriptive QC from archived VST objects; no transformation or model re-fit.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);root<-args[1];out<-args[2];dest<-file.path(out,'02_Main_Figures/source_data/core');rows<-list()
for(tissue in c('Duo','Jej')){
 p<-file.path(root,'results/bulk/main',paste0(tissue,'_VST_blind_QC.rds'))
 if(!file.exists(p)){rows[[tissue]]<-list(tissue=tissue,status='SOURCE_OBJECT_NOT_AVAILABLE',source=p);next}
 before<-digest::digest(file=p,algo='sha256')
 suppressPackageStartupMessages({requireNamespace('DESeq2');requireNamespace('SummarizedExperiment')})
 obj<-readRDS(p);mat<-if(is.matrix(obj))obj else SummarizedExperiment::assay(obj)
 stopifnot(is.matrix(mat),all(is.finite(mat)),ncol(mat)>1)
 r<-stats::cor(mat,method='pearson',use='everything')
 d<-as.data.frame(as.table(r),stringsAsFactors=FALSE);names(d)<-c('row_library','column_library','Pearson_r');d$tissue<-tissue;d$n_existing_VST_features<-nrow(mat)
 write.table(d,file.path(dest,paste0(tissue,'_VST_SAMPLE_CORRELATION.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
 after<-digest::digest(file=p,algo='sha256');stopifnot(identical(before,after))
 rows[[tissue]]<-list(tissue=tissue,status='DESCRIPTIVE_EXISTING_VST_CORRELATION_PASS',source=p,sha256_before=before,sha256_after=after,n_existing_features=nrow(mat),n_samples=ncol(mat),new_model_fit=FALSE,new_PCA=FALSE,new_VST=FALSE)
}
write_json(rows,file.path(out,'logs/EXISTING_VST_CORRELATION_QC.json'),pretty=TRUE,auto_unbox=TRUE)
cat('EXISTING_VST_DESCRIPTIVE_CORRELATION_FINISHED\n')
