#!/usr/bin/env Rscript
# Fixed, small QC sensitivity; does not choose a cohort from DE significance.
# Any distinct sample cohort is fitted afresh using the SAME bulk_run.R logic.
# Identical sets are explicitly not refitted; they link the new main analysis.

args0 <- commandArgs(trailingOnly=FALSE)
script_arg <- args0[startsWith(args0,"--file=")]
if(length(script_arg)!=1L) stop("Cannot locate this script")
script_path <- normalizePath(sub("^--file=","",script_arg),mustWork=TRUE)
source(file.path(dirname(script_path),"bulk_run.R"),local=FALSE)

bulk_compare_sensitivity <- function(main,sensitivity,out,label) {
  rows <- list()
  for(contrast in names(main$contrasts)) {
    a <- as.data.frame(data.table::fread(main$contrasts[[contrast]]$DESEQ2,sep="\t"))
    b <- as.data.frame(data.table::fread(sensitivity$contrasts[[contrast]]$DESEQ2,sep="\t"))
    if(anyDuplicated(a$gene_id) || anyDuplicated(b$gene_id)) stop("Cannot compare duplicated feature IDs")
    keep <- c("gene_id","gene_symbol","baseMean","log2FoldChange","stat","pvalue","padj","P_NA_reason","padj_NA_reason")
    joined <- merge(a[,keep],b[,keep],by=c("gene_id","gene_symbol"),all=TRUE,suffixes=c("_main","_sensitivity"),sort=FALSE)
    finite <- is.finite(joined$log2FoldChange_main) & is.finite(joined$log2FoldChange_sensitivity)
    joined$finite_LFC_both <- finite
    joined$LFC_sign_agrees <- ifelse(finite,sign(joined$log2FoldChange_main)==sign(joined$log2FoldChange_sensitivity),NA)
    joined$main_FDR05 <- ifelse(is.na(joined$padj_main),NA,joined$padj_main<.05)
    joined$sensitivity_FDR05 <- ifelse(is.na(joined$padj_sensitivity),NA,joined$padj_sensitivity<.05)
    bulk_tsv(joined,file.path(out,paste0(label,"_",contrast,"_ALL_GENE_COMPARISON.tsv")))
    rows[[contrast]] <- data.frame(sensitivity=label,contrast=contrast,n_features_union=nrow(joined),
      n_finite_LFC_both=sum(finite),LFC_spearman=if(sum(finite)>2) stats::cor(joined$log2FoldChange_main[finite],joined$log2FoldChange_sensitivity[finite],method="spearman") else NA_real_,
      LFC_sign_agreement_fraction=if(any(finite)) mean(joined$LFC_sign_agrees[finite]) else NA_real_,
      n_main_FDR05=sum(joined$main_FDR05,na.rm=TRUE),n_sensitivity_FDR05=sum(joined$sensitivity_FDR05,na.rm=TRUE),
      n_both_FDR05=sum(joined$main_FDR05 & joined$sensitivity_FDR05,na.rm=TRUE),
      n_main_significant_sensitivity_NA=sum(joined$main_FDR05 & is.na(joined$sensitivity_FDR05),na.rm=TRUE))
  }
  data.table::rbindlist(rows)
}

sensitivity_main <- function() {
  bulk_require(c("DESeq2","fgsea","data.table","jsonlite","digest","msigdbr"))
  a <- bulk_parse(commandArgs(trailingOnly=TRUE))
  if(is.null(a$root)) stop("--root is required")
  root <- bulk_root(a$root)
  gate <- if(is.null(a$global_gate)) file.path(root,"manifest/PHASE0_GATE.json") else a$global_gate
  main_path <- file.path(root,"results/bulk/BULK_STATUS.json")
  main <- jsonlite::read_json(main_path,simplifyVector=FALSE)
  if(!identical(main$status,"ANALYSIS_COMPLETE") || !isTRUE(main$models_all_new_fit) || !identical(main$tag,"main"))
    stop("This NEW project's main analysis must be complete before fixed sensitivity fits")
  folder <- file.path(root,"results/bulk/sensitivity")
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  status_path <- file.path(folder,"BULK_SENSITIVITY_STATUS.json")
  if(file.exists(status_path)) stop("Sensitivity status already exists; preserve it and audit, do not silently repeat")
  begin <- bulk_now(); records <- list(); comparisons <- list()
  tryCatch({
    input <- bulk_load_inputs(root,gate,"main_analysis_pass")
    protocol <- input$protocol
    labels <- names(protocol$sensitivity_min_non_rrna_counts)
    if(!setequal(labels,c("qc_low_750k","qc_high_1250k"))) stop("Unexpected sensitivity protocol labels")
    main_libraries <- sort(input$samples$library_id)
    for(label in labels) {
      raw_meta <- input$all_samples
      if(!label %in% names(raw_meta)) stop("Frozen sensitivity sample indicator missing")
      libraries <- sort(raw_meta$library_id[bulk_bool(raw_meta[[label]])])
      if(identical(libraries,main_libraries)) {
        records[[label]] <- list(status="IDENTICAL_TO_NEW_MAIN_INPUT_NO_REDUNDANT_REFIT",
          selection_column=label,n_libraries=length(libraries),main_analysis_status=main_path,
          main_analysis_status_sha256=bulk_sha(main_path),is_new_fit=FALSE,
          note="The exact same current-run raw-count cohort and fixed model imply the same analysis; no old project models read.")
      } else {
        fit <- run_bulk_analysis(root,gate,label,label)
        records[[label]] <- list(status=fit$status,selection_column=label,n_libraries=length(libraries),
          status_file=file.path(root,"results/bulk/sensitivity",label,"BULK_STATUS.json"),
          added_vs_main=setdiff(libraries,main_libraries),removed_vs_main=setdiff(main_libraries,libraries),is_new_fit=TRUE)
        comparisons[[label]] <- bulk_compare_sensitivity(main,fit,folder,label)
      }
    }
    if(length(comparisons)) bulk_tsv(data.table::rbindlist(comparisons),file.path(folder,"QC_SENSITIVITY_COMPARISON_SUMMARY.tsv"))
    bulk_json(status_path,list(status="ANALYSIS_COMPLETE",start_utc=begin,finish_utc=bulk_now(),
      protocol_sha256=bulk_sha(input$paths[["protocol"]]),main_status_sha256=bulk_sha(main_path),
      rules="Predefined 0.75M and 1.25M non-rRNA count thresholds; detected>=5000 and mapping/model/GSEA unchanged; no significance-driven selection",
      analyses=records,selection_of_primary_unchanged=TRUE,old_project_models_read=FALSE))
    message("ANALYSIS_COMPLETE ",status_path)
  },error=function(e) {
    bulk_json(status_path,list(status="FAILED",start_utc=begin,finish_utc=bulk_now(),analyses=records,error=conditionMessage(e)))
    stop(e)
  })
}

if(sys.nframe()==0L) sensitivity_main()
