#!/usr/bin/env Rscript
# IBS clean rerun: new DESeq2 models and full preranked fgsea analyses.
# Only raw-derived files in this NEW project and immutable reference sets are read.
# The parent run_logged.py records the OS command; this file times actual R calls.

options(warn=1, stringsAsFactors=FALSE, mc.cores=1)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

bulk_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz="UTC")
bulk_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
  if(length(missing)) stop("Missing packages in the independent project library: ", paste(missing, collapse=", "))
}
bulk_sha <- function(path) {
  if(!file.exists(path)) stop("Missing hash input: ", path)
  digest::digest(file=path, algo="sha256", serialize=FALSE)
}
bulk_json <- function(path, obj) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(obj, tmp, pretty=TRUE, auto_unbox=TRUE, na="null", null="null", digits=17)
  if(!file.rename(tmp,path)) stop("Cannot atomically write ", path)
}
bulk_tsv <- function(x,path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  data.table::fwrite(x, path, sep="\t", quote=FALSE, na="NA")
}
bulk_bool <- function(x) as.character(x) %in% c("True","TRUE","true","1")
bulk_parse <- function(args) {
  out <- list(mode="run", tag="main", selection="main_analysis_pass")
  i <- 1L
  while(i <= length(args)) {
    if(!startsWith(args[i],"--") || i == length(args)) stop("All arguments require --name value: ", args[i])
    out[[gsub("-","_",substring(args[i],3L))]] <- args[i+1L]
    i <- i+2L
  }
  out
}
bulk_root <- function(root) {
  root <- normalizePath(root, mustWork=TRUE)
  if(basename(root) != "IBS_A_OFFICIAL_RERUN_20261105") stop("Unauthorized project root")
  root
}
bulk_input_hashes <- function(paths) {
  paths <- unique(paths)
  stats::setNames(vapply(paths, bulk_sha, character(1)), paths)
}

bulk_context <- function(root,out,run_id,input_hashes) {
  ctx <- new.env(parent=emptyenv())
  ctx$root <- root; ctx$out <- out; ctx$run_id <- run_id
  ctx$input_hashes <- input_hashes
  ctx$trace <- file.path(out,"BULK_MODEL_EXECUTION.tsv")
  ctx$versions <- paste(c(paste0("R=", getRversion()),
      vapply(c("DESeq2","fgsea","data.table","jsonlite","msigdbr","digest"),
      function(p) paste0(p,"=",utils::packageVersion(p)), character(1))),collapse=";")
  ctx
}

bulk_timed <- function(ctx,stage,tissue="NA",contrast="NA",collection="NA",command,fun,is_fit=TRUE) {
  start <- bulk_now(); timer <- proc.time()[["elapsed"]]
  exit <- 0L; message <- ""; result <- NULL; error <- NULL
  tryCatch({result <- fun()}, error=function(e) {
    exit <<- 1L; message <<- conditionMessage(e); error <<- e
  })
  finish <- bulk_now(); elapsed <- proc.time()[["elapsed"]]-timer
  entry <- data.frame(module=stage,tissue=tissue,contrast=contrast,collection=collection,
    source_command=command,start_iso8601=start,finish_iso8601=finish,
    duration_seconds=elapsed,exit_code=exit,hostname=Sys.info()[["nodename"]],
    software_version=ctx$versions,run_id=ctx$run_id,is_new_fit=is_fit,
    reused_only_immutable_raw_inputs=TRUE,
    input_hashes=as.character(jsonlite::toJSON(as.list(ctx$input_hashes),auto_unbox=TRUE)),
    notes=message,stringsAsFactors=FALSE)
  data.table::fwrite(entry,ctx$trace,sep="\t",quote=TRUE,append=file.exists(ctx$trace),
                     col.names=!file.exists(ctx$trace),na="NA")
  if(!is.null(error)) stop(error)
  result
}

bulk_prepare_reference <- function(root) {
  bulk_require(c("data.table","jsonlite","digest","msigdbr"))
  data.table::setDTthreads(1L)
  folder <- file.path(root,"data/reference/bulk")
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  manifest_path <- file.path(folder,"GENE_SETS_MANIFEST.json")
  if(file.exists(manifest_path)) {
    old <- jsonlite::read_json(manifest_path,simplifyVector=TRUE)
    for(k in names(old$files_sha256)) {
      if(bulk_sha(file.path(folder,k)) != old$files_sha256[[k]]) stop("Gene-set reference hash mismatch: ",k)
    }
    return(manifest_path)
  }
  start <- bulk_now()
  formals_names <- names(formals(msigdbr::msigdbr))
  call_sets <- function(category,subcategory=NULL) {
    if("collection" %in% formals_names) {
      a <- list(species="Homo sapiens",collection=category)
      if(!is.null(subcategory)) a$subcollection <- subcategory
    } else {
      a <- list(species="Homo sapiens",category=category)
      if(!is.null(subcategory)) a$subcategory <- subcategory
    }
    as.data.frame(do.call(msigdbr::msigdbr,a))
  }
  sets <- list(Hallmark=call_sets("H"),Reactome=call_sets("C2","CP:REACTOME"))
  files <- character(); snapshots <- list(); counts <- list()
  for(collection in names(sets)) {
    tab <- sets[[collection]]
    if(!all(c("gs_name","gene_symbol") %in% names(tab)) || !nrow(tab)) stop("Invalid msigdbr schema for ",collection)
    for(col in names(tab)) if(is.list(tab[[col]])) tab[[col]] <- vapply(tab[[col]],paste,collapse=";",FUN.VALUE=character(1))
    source_name <- paste0(collection,"_msigdbr_source.tsv")
    bulk_tsv(tab,file.path(folder,source_name))
    pathways <- split(as.character(tab$gene_symbol),as.character(tab$gs_name))
    pathways <- lapply(pathways,function(x) sort(unique(x[!is.na(x) & nzchar(x)])))
    pathways <- pathways[order(names(pathways))]
    gmt_name <- paste0(collection,"_Homo_sapiens.gmt")
    lines <- vapply(names(pathways),function(n) paste(c(n,"MSigDB_via_msigdbr",pathways[[n]]),collapse="\t"),character(1))
    writeLines(lines,file.path(folder,gmt_name),useBytes=TRUE)
    snapshots[[collection]] <- if("db_version" %in% names(tab)) unique(tab$db_version) else
      paste0("msigdbr packaged snapshot; package ",utils::packageVersion("msigdbr"),"; database version not exposed as field")
    counts[[collection]] <- length(pathways)
    files <- c(files,source_name,gmt_name)
  }
  hashes <- stats::setNames(vapply(file.path(folder,files),bulk_sha,character(1)),files)
  bulk_json(manifest_path,list(status="REFERENCE_FROZEN",created_start_utc=start,created_finish_utc=bulk_now(),
    source="New extraction from version-locked msigdbr package; no previous-project GMT used",
    source_urls=c("https://igordot.github.io/msigdbr/","https://www.gsea-msigdb.org/gsea/msigdb/"),
    package_version=as.character(utils::packageVersion("msigdbr")),species="Homo sapiens",identifier_namespace="gene_symbol",
    database_snapshot=snapshots,pathway_counts=counts,files_sha256=as.list(hashes)))
  manifest_path
}

bulk_read_reference <- function(root) {
  manifest <- bulk_prepare_reference(root)
  folder <- dirname(manifest)
  sets <- lapply(c("Hallmark","Reactome"),function(collection) {
    fgsea::gmtPathways(file.path(folder,paste0(collection,"_Homo_sapiens.gmt")))
  })
  names(sets) <- c("Hallmark","Reactome")
  list(manifest=manifest,sets=sets,files=c(manifest,
    file.path(folder,c("Hallmark_Homo_sapiens.gmt","Reactome_Homo_sapiens.gmt"))))
}

bulk_load_inputs <- function(root,gate,selection) {
  gate <- normalizePath(gate,mustWork=TRUE)
  if(!startsWith(gate,paste0(root,"/")) && !startsWith(gate,paste0(root,"\\"))) stop("Global gate must be in the new project")
  g <- jsonlite::read_json(gate,simplifyVector=TRUE)
  if(!identical(g$status,"PASS")) stop("Global Phase0 has not PASSed; no statistical fit may start")
  paths <- c(status=file.path(root,"data/metadata/bulk/BULK_PHASE0_STATUS.json"),
    protocol=file.path(root,"data/metadata/bulk/QC_PROTOCOL_FROZEN.json"),
    samples=file.path(root,"data/metadata/bulk/SAMPLE_MANIFEST.tsv"),
    matrix=file.path(root,"data/processed/bulk/raw_non_rrna_counts.tsv.gz"),
    annotation=file.path(root,"data/processed/bulk/gene_annotation.tsv"))
  s <- jsonlite::read_json(paths[["status"]],simplifyVector=TRUE)
  if(!identical(s$status,"PASS")) stop("bulk Phase0 is not PASS")
  expected <- c(protocol=s$protocol_sha256,samples=s$manifest_sha256,matrix=s$matrix_sha256,annotation=s$annotation_sha256)
  for(k in names(expected)) if(bulk_sha(paths[[k]]) != expected[[k]]) stop("Frozen Phase0 hash mismatch: ",k)
  protocol <- jsonlite::read_json(paths[["protocol"]],simplifyVector=TRUE)
  if(!grepl("smallest retained group",protocol$low_expression_filter,fixed=TRUE)) stop("Unexpected frozen low-count filter")
  if(!grepl("pooled",protocol$gsea$FDR,ignore.case=TRUE)) stop("Protocol must state primary pooled GSEA FDR")
  samples <- as.data.frame(data.table::fread(paths[["samples"]],sep="\t",colClasses="character",na.strings=NULL))
  if(!selection %in% names(samples)) stop("Unknown frozen sample-selection column: ",selection)
  selected <- bulk_bool(samples[[selection]])
  if(!any(selected) || any(samples$mapping_status[selected] != "UNIQUE_EXACT_GEO_TITLE")) stop("Invalid or unmapped sample selected")
  if(anyDuplicated(samples$library_id)) stop("Duplicate count-library ID")
  counts <- as.data.frame(data.table::fread(cmd=paste("gzip -dc",shQuote(paths[["matrix"]])),sep="\t",check.names=FALSE))
  if(!all(c("gene_id","gene_symbol","gene_type") %in% names(counts))) stop("Processed raw count schema missing annotation fields")
  if(anyDuplicated(counts$gene_id)) stop("Duplicate Ensembl count features")
  required <- samples$library_id[selected]
  if(!all(required %in% names(counts))) stop("Selected count columns missing")
  x <- as.matrix(counts[,required,drop=FALSE])
  if(any(!is.finite(x)) || any(x < 0) || any(x != floor(x)) || any(x > .Machine$integer.max)) stop("Noninteger/invalid raw counts")
  storage.mode(x) <- "integer"; rownames(x) <- counts$gene_id
  anno <- as.data.frame(data.table::fread(paths[["annotation"]],sep="\t",colClasses="character",na.strings=NULL))
  if(anyDuplicated(anno$gene_id)) stop("Duplicate annotation genes")
  if(!setequal(counts$gene_id,anno$gene_id[!bulk_bool(anno$rrna_excluded)])) stop("Non-rRNA annotation/matrix mismatch")
  meta <- samples[selected,,drop=FALSE]
  meta <- meta[match(colnames(x),meta$library_id),,drop=FALSE]
  rownames(meta) <- meta$library_id
  if(!identical(meta$library_id,colnames(x))) stop("Count metadata order mismatch")
  list(paths=paths,gate=gate,protocol=protocol,phase0=s,samples=meta,all_samples=samples,
       counts=x,annotation=anno,matrix_annotation=counts[,c("gene_id","gene_symbol","gene_type")])
}

bulk_qc_pca <- function(dds,meta,out,tissue,ctx) {
  vsd <- bulk_timed(ctx,"VST_PCA",tissue,command=
    "varianceStabilizingTransformation(dds, blind=TRUE); prcomp(t(top500_variable_genes)); report_only_no_sample_deletion",fun=function() {
      v <- DESeq2::varianceStabilizingTransformation(dds,blind=TRUE)
      x <- SummarizedExperiment::assay(v)
      vars <- matrixStats::rowVars(x)
      selected <- head(order(-vars,rownames(x)),min(500L,nrow(x)))
      p <- stats::prcomp(t(x[selected,,drop=FALSE]),center=TRUE,scale.=FALSE)
      pcs <- as.data.frame(p$x[,seq_len(min(10L,ncol(p$x))),drop=FALSE])
      pcs$library_id <- rownames(pcs)
      pcs <- merge(meta,pcs,by="library_id",sort=FALSE)
      pc_dist <- sqrt(pcs$PC1^2+pcs$PC2^2)
      limit <- stats::median(pc_dist)+4*stats::mad(pc_dist)
      pcs$PCA_PC1_PC2_distance <- pc_dist
      pcs$PCA_report_only_4MAD_flag <- pc_dist > limit
      pcs$PCA_report_only_4MAD_cutoff <- limit
      bulk_tsv(pcs,file.path(out,paste0(tissue,"_PCA_SAMPLE_QC.tsv")))
      bulk_tsv(data.frame(PC=paste0("PC",seq_along(p$sdev)),variance_fraction=p$sdev^2/sum(p$sdev^2)),
               file.path(out,paste0(tissue,"_PCA_VARIANCE.tsv")))
      saveRDS(v,file.path(out,paste0(tissue,"_VST_blind_QC.rds")))
      v
    },is_fit=FALSE)
  invisible(vsd)
}

bulk_fit_tissue <- function(input, tissue, out, ctx) {
  meta <- input$samples[input$samples$tissue == tissue,,drop=FALSE]
  meta$group <- factor(meta$group,levels=c("HV","IBSC","IBSD"))
  group_counts <- table(meta$group)
  if(any(group_counts < 2L)) stop(tissue,": each group requires at least two libraries")
  valid_subject <- !is.na(meta$subject_id) & !meta$subject_id %in% c("NA","")
  if(anyDuplicated(meta$subject_id[valid_subject])) stop(tissue,": verified within-tissue repeats require an explicit model decision")
  meta$sex <- factor(meta$sex,levels=c("F","M"))
  design <- stats::as.formula("~group")
  design_reason <- "Sex unavailable, invariant or confounded: use simplest full-rank group model"
  if(!anyNA(meta$sex) && nlevels(droplevels(meta$sex)) == 2L) {
    candidate <- stats::model.matrix(~sex+group,data=meta)
    if(qr(candidate)$rank == ncol(candidate)) {
      design <- stats::as.formula("~sex+group")
      design_reason <- "Complete explicitly recorded sex; sex+group design full rank"
    }
  }
  mm <- stats::model.matrix(design,data=meta)
  if(qr(mm)$rank != ncol(mm) || nrow(mm) <= ncol(mm)) stop("Nonidentifiable tissue model")
  bulk_tsv(data.frame(library_id=meta$library_id,mm,check.names=FALSE),file.path(out,paste0(tissue,"_MODEL_MATRIX.tsv")))
  bulk_tsv(as.data.frame(table(group=meta$group,sex=meta$sex,useNA="ifany")),file.path(out,paste0(tissue,"_GROUP_SEX_CONTINGENCY.tsv")))
  bulk_tsv(data.frame(library_id=meta$library_id,group=meta$group,batch=meta$batch),file.path(out,paste0(tissue,"_BATCH_AUDIT.tsv")))
  x_all <- input$counts[,meta$library_id,drop=FALSE]
  qc <- data.frame(library_id=meta$library_id,gsm=meta$gsm,tissue=tissue,group=meta$group,sex=meta$sex,
    non_rrna_library_counts=colSums(x_all),detected_non_rrna_genes=colSums(x_all>0),
    subject_id=meta$subject_id,subject_id_basis=meta$subject_id_basis,batch=meta$batch)
  bulk_tsv(qc,file.path(out,paste0(tissue,"_PREMODEL_LIBRARY_QC.tsv")))
  min_group <- min(as.integer(group_counts))
  keep <- rowSums(x_all >= 10L) >= min_group
  filter_tab <- data.frame(gene_id=rownames(x_all),n_libraries_count_ge10=rowSums(x_all >= 10L),
    required_libraries=min_group,passed_low_expression_filter=keep,all_zero=rowSums(x_all)==0)
  bulk_tsv(filter_tab,file.path(out,paste0(tissue,"_GENE_FILTER.tsv")))
  if(sum(keep) < 100L) stop("Insufficient genes for a reliable RNA-seq model")
  x <- x_all[keep,,drop=FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=x,colData=meta,design=design)
  # QC VST/PCA is blind and precedes formal fitting. Flags do not delete samples.
  bulk_qc_pca(dds,meta,out,tissue,ctx)
  dds <- bulk_timed(ctx,"DESEQ2_NEW_FIT",tissue,command=
    "DESeq(dds,test='Wald',fitType='parametric',sfType='ratio',betaPrior=FALSE,minReplicatesForReplace=7,parallel=FALSE)",fun=function() {
      DESeq2::DESeq(dds,test="Wald",fitType="parametric",sfType="ratio",betaPrior=FALSE,
                   minReplicatesForReplace=7,parallel=FALSE,quiet=FALSE)
    })
  saveRDS(dds,file.path(out,paste0(tissue,"_DESeq2_NEW_FIT.rds")))
  fitted_meta <- as.data.frame(SummarizedExperiment::colData(dds))
  fitted_meta$library_id <- rownames(fitted_meta)
  fitted_meta$size_factor <- DESeq2::sizeFactors(dds)
  bulk_tsv(fitted_meta,file.path(out,paste0(tissue,"_FITTED_SAMPLE_METADATA.tsv")))
  gene_qc <- as.data.frame(S4Vectors::mcols(dds))
  gene_qc$gene_id <- rownames(gene_qc)
  for(k in names(gene_qc)) if(is.list(gene_qc[[k]])) gene_qc[[k]] <- vapply(gene_qc[[k]],paste,collapse=";",FUN.VALUE=character(1))
  bulk_tsv(gene_qc,file.path(out,paste0(tissue,"_MODEL_GENE_DIAGNOSTICS.tsv")))
  conv <- S4Vectors::mcols(dds)$betaConv
  dispersion <- DESeq2::dispersions(dds)
  sf <- DESeq2::sizeFactors(dds)
  if(is.null(conv) || anyNA(conv) || any(!conv) || any(!is.finite(dispersion)) || any(!is.finite(sf)) || any(sf <= 0))
    stop("New model has nonconverged genes or invalid dispersion/size factors; diagnostics saved, do not mark complete")
  cooks <- SummarizedExperiment::assays(dds)[["cooks"]]
  if(!is.null(cooks)) {
    f_cut <- stats::qf(.99,ncol(mm),nrow(mm)-ncol(mm))
    bulk_tsv(data.frame(library_id=colnames(cooks),max_cooks=apply(cooks,2,max,na.rm=TRUE),
      n_genes_cooks_above_reference_F99=colSums(cooks>f_cut,na.rm=TRUE),reference_F99=f_cut),
      file.path(out,paste0(tissue,"_COOKS_SAMPLE_DIAGNOSTICS.tsv")))
  }
  model_info <- list(tissue=tissue,n_samples=ncol(dds),n_genes_fitted=nrow(dds),
    group_counts=as.list(group_counts),design=paste(deparse(design),collapse=""),design_reason=design_reason,
    matrix_rank=qr(mm)$rank,matrix_columns=ncol(mm),prefilter_min_count=10,prefilter_min_libraries=min_group,
    beta_converged=sum(conv),beta_not_converged=sum(!conv),all_count_assay_integer=is.integer(DESeq2::counts(dds)),
    batch_status=if(all(meta$batch %in% c("NA",""))) "NO_VERIFIED_TECHNICAL_BATCH_AVAILABLE" else "BATCH_FIELDS_RECORDED_NOT_MODELED_DISCLOSE_CONFOUNDING",
    subject_pairing=if(all(!valid_subject)) "NO_VERIFIED_SUBJECT_ID_NO_PAIRING_INFERENCE" else "EXPLICIT_IDS_ONLY_TISSUES_MODELED_SEPARATELY",
    fitted_object=file.path(out,paste0(tissue,"_DESeq2_NEW_FIT.rds")),fitted_object_sha256=bulk_sha(file.path(out,paste0(tissue,"_DESeq2_NEW_FIT.rds"))))
  bulk_json(file.path(out,paste0(tissue,"_MODEL_AUDIT.json")),model_info)
  list(dds=dds,filter=filter_tab,info=model_info)
}

bulk_result <- function(input,fit,tissue,group,out,ctx) {
  label <- paste0(tissue,"_",group,"_vs_HV")
  res <- bulk_timed(ctx,"DESEQ2_RESULTS",tissue,label,command=paste0(
    "results(dds,contrast=c('group','",group,"','HV'),alpha=0.05,pAdjustMethod='BH',independentFiltering=TRUE,lfcThreshold=0)"),fun=function() {
      DESeq2::results(fit$dds,contrast=c("group",group,"HV"),alpha=.05,
        pAdjustMethod="BH",independentFiltering=TRUE,lfcThreshold=0)
    },is_fit=FALSE)
  model <- as.data.frame(res)
  model$gene_id <- rownames(model)
  tab <- input$annotation
  tab$gene_id_stable <- sub("\\.[0-9]+$","",tab$gene_id)
  idx <- match(tab$gene_id,model$gene_id)
  for(col in c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")) tab[[col]] <- model[[col]][idx]
  tab$tissue <- tissue; tab$contrast <- label
  tab$P_NA_reason <- "NONE"; tab$padj_NA_reason <- "NONE"
  rrna <- bulk_bool(tab$rrna_excluded)
  ft <- fit$filter[match(tab$gene_id,fit$filter$gene_id),,drop=FALSE]
  low <- !rrna & (is.na(ft$passed_low_expression_filter) | !ft$passed_low_expression_filter)
  zero <- !rrna & !is.na(ft$all_zero) & ft$all_zero
  tab$P_NA_reason[rrna] <- "RRNA_EXCLUDED_BEFORE_MODEL"
  tab$P_NA_reason[low] <- "LOW_EXPRESSION_PREFILTER_NOT_FITTED"
  tab$P_NA_reason[zero] <- "ALL_ZERO_COUNTS_NOT_FITTED"
  fitted <- !is.na(idx)
  pna <- fitted & is.na(tab$pvalue)
  mc <- S4Vectors::mcols(fit$dds)
  maxcooks <- if("maxCooks" %in% colnames(mc)) mc$maxCooks[idx] else rep(NA_real_,nrow(tab))
  mm <- stats::model.matrix(DESeq2::design(fit$dds),data=as.data.frame(SummarizedExperiment::colData(fit$dds)))
  cut <- stats::qf(.99,ncol(mm),nrow(mm)-ncol(mm))
  tab$P_NA_reason[pna] <- "DESEQ2_PVALUE_NA_REQUIRES_DIAGNOSTIC_REVIEW"
  outlier <- pna & is.finite(maxcooks) & maxcooks>cut
  tab$P_NA_reason[outlier] <- "COOKS_DISTANCE_OUTLIER_DEFAULT_FILTER"
  tab$padj_NA_reason[is.na(tab$padj)] <- tab$P_NA_reason[is.na(tab$padj)]
  tab$padj_NA_reason[is.na(tab$padj) & is.finite(tab$pvalue)] <- "INDEPENDENT_FILTERING_LOW_BASEMEAN"
  tab$DEG_FDR05 <- !is.na(tab$padj) & tab$padj<.05
  tab$DEG_FDR05_ABSLFC1 <- tab$DEG_FDR05 & !is.na(tab$log2FoldChange) & abs(tab$log2FoldChange)>=1
  tab$direction <- ifelse(is.na(tab$log2FoldChange),"NA",ifelse(tab$log2FoldChange>0,"UP",ifelse(tab$log2FoldChange<0,"DOWN","ZERO")))
  order_cols <- c("gene_symbol","gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  tab <- tab[,c(order_cols,setdiff(names(tab),order_cols)),drop=FALSE]
  bulk_tsv(tab,file.path(out,paste0(label,".tsv")))
  metadata <- S4Vectors::metadata(res)
  bulk_json(file.path(out,paste0(label,"_RESULTS_METADATA.json")),list(contrast=label,
    direction=paste(group,"relative to HV"),alpha=.05,lfcThreshold=0,
    independent_filter_threshold=metadata$filterThreshold,
    independent_filter_theta=metadata$filterTheta,
    pvalue_na_reasons=as.list(table(tab$P_NA_reason[is.na(tab$pvalue)])),
    padj_na_reasons=as.list(table(tab$padj_NA_reason[is.na(tab$padj)])),
    result_descriptions=as.data.frame(S4Vectors::mcols(res))))
  list(label=label,table=tab)
}

bulk_gsea <- function(result,refs,out,ctx,seed) {
  tab <- result$table; label <- result$label
  eligible <- is.finite(tab$stat) & is.finite(tab$pvalue) & !is.na(tab$gene_symbol) & nzchar(tab$gene_symbol) & tab$gene_symbol!="NA"
  ranking <- tab
  ranking$GSEA_selected <- FALSE
  ranking$GSEA_exclusion_reason <- ifelse(eligible,"DUPLICATE_SYMBOL_NONREPRESENTATIVE",
    ifelse(is.na(tab$gene_symbol) | !nzchar(tab$gene_symbol) | tab$gene_symbol=="NA","MISSING_GENE_SYMBOL","NONFINITE_WALD_OR_PVALUE"))
  candidates <- which(eligible)
  candidates <- candidates[order(-tab$baseMean[candidates],tab$gene_id[candidates],na.last=TRUE)]
  chosen <- candidates[!duplicated(tab$gene_symbol[candidates])]
  ranking$GSEA_selected[chosen] <- TRUE
  ranking$GSEA_exclusion_reason[chosen] <- "SELECTED_HIGHEST_BASEMEAN_THEN_GENE_ID"
  ranks <- stats::setNames(tab$stat[chosen],tab$gene_symbol[chosen])
  ranks <- ranks[order(-ranks,names(ranks))]
  if(length(ranks)<100L || anyDuplicated(names(ranks))) stop("Invalid full gene ranking")
  bulk_tsv(ranking,file.path(out,paste0(label,"_GSEA_RANK_SELECTION_ALL_FEATURES.tsv")))
  bulk_tsv(data.frame(gene_symbol=names(ranks),wald_stat=unname(ranks),rank=seq_along(ranks)),file.path(out,paste0(label,"_GSEA_WALD_RANKS.tsv")))
  tables <- list()
  for(i in seq_along(refs$sets)) {
    collection <- names(refs$sets)[i]
    pathways <- refs$sets[[i]]
    overlap <- vapply(pathways,function(x) length(intersect(x,names(ranks))),integer(1))
    selected <- overlap>=15L & overlap<=500L
    if(!any(selected)) stop("No eligible pathways for ",collection)
    this_seed <- as.integer(seed+i-1L)
    raw <- bulk_timed(ctx,"FGSEA_MULTILEVEL_NEW",substr(label,1L,3L),label,collection,
      command=paste0("set.seed(",this_seed,"); fgseaMultilevel(pathways=eligible_",collection,
        ",stats=all_finite_Wald_ranks,minSize=15,maxSize=500,eps=0,nPermSimple=1000,nproc=1)"),fun=function() {
        set.seed(this_seed)
        fgsea::fgseaMultilevel(pathways=pathways[selected],stats=ranks,minSize=15,maxSize=500,
          eps=0,nPermSimple=1000,nproc=1,scoreType="std")
      })
    saveRDS(raw,file.path(out,paste0(label,"_",collection,"_FGSEA_NEW.rds")))
    raw <- as.data.frame(raw)
    if(anyDuplicated(raw$pathway) || !setequal(raw$pathway,names(pathways)[selected])) stop("fgsea output does not cover all eligible pathways")
    raw$leadingEdge <- vapply(raw$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
    # The package padj is within-collection. The primary pooled BH is applied
    # once both collections are available, without looking at target genes.
    raw$padj_collection <- raw$padj
    raw$padj <- NULL
    full <- data.frame(collection=collection,pathway=names(pathways),
      gene_set_size_original=vapply(pathways,length,integer(1)),rank_overlap_size=as.integer(overlap),
      tested=selected,stringsAsFactors=FALSE)
    mi <- match(full$pathway,raw$pathway)
    for(k in setdiff(names(raw),"pathway")) full[[k]] <- raw[[k]][mi]
    full$NA_reason <- ifelse(!selected,"UNTESTED_RANK_OVERLAP_OUTSIDE_15_500",
      ifelse(is.na(full$pval),"FGSEA_P_NA_ESTIMATION_IMBALANCE_OR_FAILURE_CHECK_LOG","NONE"))
    full$seed <- this_seed; full$contrast <- label
    tables[[collection]] <- full
  }
  all <- data.table::rbindlist(tables,use.names=TRUE,fill=TRUE)
  # Failed estimates stay NA and still count in the tested family, rather than
  # making the adjustment more permissive by silently reducing its size.
  all$padj <- stats::p.adjust(all$pval,method="BH",n=sum(all$tested))
  all$padj_pooled <- all$padj
  bulk_tsv(all,file.path(out,paste0(label,"_GSEA_ALL.tsv")))
  bulk_json(file.path(out,paste0(label,"_GSEA_AUDIT.json")),list(
    status="ALL_ELIGIBLE_PATHWAYS_ATTEMPTED",contrast=label,n_ranked_symbols=length(ranks),
    n_rank_statistic_ties=sum(duplicated(unname(ranks))),rank_direction="positive Wald=IBS subgroup higher than HV",
    ranking_selection="Finite Wald and finite pvalue; no significance threshold; highest baseMean duplicate-symbol representative",
    primary_FDR="BH across both tested collections within this contrast",supplementary_FDR="within-collection BH",
    n_pathways_declared=nrow(all),n_pathways_tested=sum(all$tested),
    n_pathways_estimation_NA=sum(all$tested & is.na(all$pval)),n_FDR05=sum(all$padj<.05,na.rm=TRUE),
    nPermSimple_interpretation="Initial simple sampling in the adaptive multilevel procedure; not fixed final 1000 permutations per pathway"))
  as.data.frame(all)
}

run_bulk_analysis <- function(root,gate,selection="main_analysis_pass",tag="main") {
  bulk_require(c("DESeq2","fgsea","data.table","jsonlite","digest","msigdbr","matrixStats","BiocParallel"))
  data.table::setDTthreads(1L)
  BiocParallel::register(BiocParallel::SerialParam())
  root <- bulk_root(root)
  if(!grepl("^[A-Za-z0-9_-]+$",tag)) stop("Invalid output tag")
  out <- if(tag=="main") file.path(root,"results/bulk/main") else file.path(root,"results/bulk/sensitivity",tag)
  if(dir.exists(out) && length(list.files(out,all.files=FALSE))) stop(
    "Output directory already contains this NEW project's execution artifacts. Do not silently overwrite or reuse fits: ",out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  started <- bulk_now(); run_id <- paste0("bulk_",tag,"_",format(Sys.time(),"%Y%m%dT%H%M%S",tz="UTC"),"_",Sys.getpid())
  status_path <- file.path(out,"BULK_STATUS.json")
  status <- list(status="RUNNING",tag=tag,selection_column=selection,run_id=run_id,start_utc=started,
                 clean_rerun=TRUE,old_models_or_results_read=FALSE)
  bulk_json(status_path,status)
  tryCatch({
    input <- bulk_load_inputs(root,gate,selection)
    refs <- bulk_read_reference(root)
    hashes <- bulk_input_hashes(c(unname(input$paths),input$gate,refs$files))
    ctx <- bulk_context(root,out,run_id,hashes)
    writeLines(capture.output(sessionInfo()),file.path(out,"R_SESSIONINFO_START.txt"))
    bulk_json(file.path(out,"INPUT_SHA256.json"),as.list(hashes))
    bulk_tsv(input$samples,file.path(out,"ANALYSIS_SAMPLE_MANIFEST.tsv"))
    fits <- list(); summary <- list(); contrasts <- list(); index <- 0L
    for(tissue in c("Duo","Jej")) {
      fit <- bulk_fit_tissue(input,tissue,out,ctx)
      fits[[tissue]] <- fit$info
      for(group in c("IBSD","IBSC")) {
        index <- index+1L
        result <- bulk_result(input,fit,tissue,group,out,ctx)
        gsea <- bulk_gsea(result,refs,out,ctx,as.integer(input$protocol$gsea$seed_base)+index*10L)
        tab <- result$table
        summary[[result$label]] <- data.frame(contrast=result$label,tissue=tissue,
          n_samples=fit$info$n_samples,n_genes_fitted=fit$info$n_genes_fitted,
          n_feature_rows=nrow(tab),n_finite_p=sum(is.finite(tab$pvalue)),n_padj_NA=sum(is.na(tab$padj)),
          n_DEG_FDR05=sum(tab$DEG_FDR05),n_DEG_FDR05_ABSLFC1=sum(tab$DEG_FDR05_ABSLFC1),
          n_GSEA_tested=sum(gsea$tested),n_GSEA_primary_pooled_FDR05=sum(gsea$padj<.05,na.rm=TRUE),
          n_GSEA_estimation_NA=sum(gsea$tested & is.na(gsea$pval)))
        contrasts[[result$label]] <- list(DESEQ2=file.path(out,paste0(result$label,".tsv")),
          GSEA=file.path(out,paste0(result$label,"_GSEA_ALL.tsv")),
          RANKS=file.path(out,paste0(result$label,"_GSEA_WALD_RANKS.tsv")))
      }
      rm(fit); gc(verbose=FALSE)
    }
    bulk_tsv(data.table::rbindlist(summary),file.path(out,"BULK_CONTRAST_SUMMARY.tsv"))
    writeLines(capture.output(sessionInfo()),file.path(out,"R_SESSIONINFO_END.txt"))
    files <- list.files(out,full.names=TRUE,recursive=TRUE)
    files <- files[!dir.exists(files) & files != status_path]
    output_hashes <- bulk_input_hashes(files)
    bulk_tsv(data.frame(path=names(output_hashes),sha256=unname(output_hashes)),file.path(out,"OUTPUT_SHA256.tsv"))
    status$status <- "ANALYSIS_COMPLETE"; status$finish_utc <- bulk_now()
    status$input_sha256 <- as.list(hashes); status$model_audits <- fits
    status$contrasts <- contrasts; status$n_contrasts <- length(contrasts)
    status$models_all_new_fit <- TRUE; status$checkpoint_model_reuse <- FALSE
    status$original60_status <- input$phase0$official60_status
    status$output_hashes_manifest <- file.path(out,"OUTPUT_SHA256.tsv")
    status$output_hashes_manifest_sha256 <- bulk_sha(status$output_hashes_manifest)
    status$stage_execution_trace <- ctx$trace
    status$global_phase0_gate <- input$gate
    status$limitations <- c("Independent DESeq2 reanalysis; source publication used glmLRT.",
      "Unknown subject pairing and unverified technical batches are not invented or silently corrected.",
      "PCA flags are report-only; no candidate-dependent or significance-dependent sample removal.",
      "Missing official60 identities preclude claiming exact reproduction of that cohort.")
    bulk_json(status_path,status)
    if(tag=="main") bulk_json(file.path(root,"results/bulk/BULK_STATUS.json"),status)
    message("ANALYSIS_COMPLETE ",status_path)
    invisible(status)
  },error=function(e) {
    status$status <- "FAILED"; status$finish_utc <- bulk_now(); status$error <- conditionMessage(e)
    bulk_json(status_path,status)
    stop(e)
  })
}

bulk_main <- function() {
  a <- bulk_parse(commandArgs(trailingOnly=TRUE))
  if(is.null(a$root)) stop("--root is required")
  root <- bulk_root(a$root)
  if(a$mode=="prepare-reference") {
    cat(bulk_prepare_reference(root),"\n")
  } else if(a$mode=="run") {
    gate <- if(is.null(a$global_gate)) file.path(root,"manifest/PHASE0_GATE.json") else a$global_gate
    run_bulk_analysis(root,gate,a$selection,a$tag)
  } else stop("Unknown mode: ",a$mode)
}

if(sys.nframe()==0L) bulk_main()
