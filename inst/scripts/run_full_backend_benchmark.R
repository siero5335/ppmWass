#!/usr/bin/env Rscript
# Full clean-retrieval comparison plus repeated serial/parallel runtime tests.
# Eight display configurations, seven scientific distances; sparse is an OT solver.
opts <- list(repo=NULL, input_dir=NULL, output_dir=NULL, library=NULL,
             cores="4", sizes="20,40,80", replicates="3", seed="20260908",
             cluster_boot="20000", query_boot="1000")
for(arg in commandArgs(TRUE)) {
  eq <- regexpr("=",arg,fixed=TRUE)[1]
  if(eq<1) stop("Expected --key=value")
  key <- gsub("-","_",sub("^--","",substr(arg,1,eq-1)),fixed=TRUE)
  if(!key %in% names(opts)) stop("Unknown argument: ",arg)
  opts[[key]] <- substring(arg,eq+1)
}
stopifnot(!is.null(opts$repo),!is.null(opts$input_dir),!is.null(opts$output_dir),!is.null(opts$library))
for(key in c("cores","replicates","seed","cluster_boot","query_boot")) {
  opts[[key]] <- as.integer(opts[[key]]); stopifnot(is.finite(opts[[key]]),opts[[key]]>0)
}
sizes <- as.integer(strsplit(opts$sizes,",",fixed=TRUE)[[1]])
stopifnot(all(is.finite(sizes)),all(sizes>1))
.libPaths(c(normalizePath(opts$library),.libPaths())); library(ppmWass)
source(file.path(opts$repo,"inst/scripts/lib/publication_retrieval_statistics.R"))
stopifnot(requireNamespace("msentropy",quietly=TRUE))
dir.create(opts$output_dir,recursive=TRUE,showWarnings=FALSE)
configs <- data.frame(configuration=c("ppm_wasserstein","ppm_wasserstein_sparse","composite",
            "entropy_weighted","entropy_unweighted","cosine","weighted_cosine","hellinger"),
            distance_method=c("ppm_wasserstein","ppm_wasserstein", "composite",
            "entropy_weighted","entropy_unweighted","cosine","weighted_cosine","hellinger"),
            ot_method=c("exact","exact_sparse",rep("exact",6)),stringsAsFactors=FALSE)
p0 <- eihrms_default_params(); p0$distance_method <- "ppm_wasserstein"; p0$tol_ppm <- 15
p0$class_detection_ppm <- 15; p0$use_parallel <- TRUE; p0$n_cores <- opts$cores
input_paths <- file.path(opts$input_dir,c("RECETOX_merged.msp","HREI-MSDB.msp"))
stopifnot(all(file.exists(input_paths)))
run_spec <- list(opts=opts,params=p0,configs=configs,input_md5=tools::md5sum(input_paths),
                 version=as.character(packageVersion("ppmWass")),
                 source_md5=tools::md5sum(c(file.path(opts$repo,"inst/scripts/run_full_backend_benchmark.R"),
                   list.files(file.path(opts$repo,"R"),full.names=TRUE))))
spec_file <- file.path(opts$output_dir,"run_spec.rds")
if(file.exists(spec_file)) stopifnot(identical(readRDS(spec_file),run_spec)) else saveRDS(run_spec,spec_file)
write.csv(configs,file.path(opts$output_dir,"configurations.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(opts$output_dir,"sessionInfo.txt"))
writeLines(capture.output(Sys.info()),file.path(opts$output_dir,"machine.txt"))
subset_sp <- function(sp,idx) {
  sp$df_spec <- sp$df_spec[idx,,drop=FALSE]
  for(nm in c("frag_list","loss_list","ri","loss_typ_list","mref_confidence","loss_anchor_list",
              "loss_pair_list","loss_anchor_typ_list","loss_pair_typ_list")) if(!is.null(sp[[nm]])) sp[[nm]] <- sp[[nm]][idx]
  sp
}
write_table <- function(x,path) write.csv(x,path,row.names=FALSE)
for(dataset in c("RECETOX","HREI-MSDB")) {
  dest <- file.path(opts$output_dir,gsub("-","_",dataset)); dir.create(dest,showWarnings=FALSE)
  spectra_file <- file.path(dest,"spectra.rds")
  if(file.exists(spectra_file)) sp <- readRDS(spectra_file) else {
    input <- input_paths[if(dataset=="RECETOX") 1 else 2]
    if(dataset=="HREI-MSDB") {
      # Match the publication loader: preserve DB records with duplicate names.
      lines <- readLines(input,warn=FALSE); starts <- grep("^Name:",lines)
      ends <- c(starts[-1]-1L,length(lines))
      for(i in seq_along(starts)) {
        block <- lines[starts[i]:ends[i]]; id <- grep("^DB#:",block,value=TRUE)
        id <- if(length(id)) trimws(sub("^DB#:[[:space:]]*","",id[1])) else as.character(i)
        lines[starts[i]] <- paste0(lines[starts[i]]," __HREI_DB",id)
      }
      input <- file.path(dest,"unique_names.msp"); writeLines(lines,input,useBytes=TRUE)
    }
    raw <- build_spectra_from_msp(input,p0,require_ri=FALSE,progress=FALSE)
    df <- raw$df_spec
    valid <- !is.na(df$RI) & !is.na(df$inchikey) & nchar(df$inchikey)>=14 &
      vapply(raw$frag_list,function(x) is.matrix(x) && nrow(x)>=2 && any(x[,2]!=0),logical(1))
    sp <- subset_sp(raw,which(valid)); saveRDS(sp,spectra_file)
    write_table(data.frame(dataset,raw_records=nrow(raw$msp_metadata),prepared_records=nrow(raw$df_spec),
                         retrieval_records=nrow(sp$df_spec)),file.path(dest,"qc.csv"))
  }
  keys <- sp$df_spec$inchikey; names(keys) <- sp$df_spec$id
  message(dataset,": full retrieval on ",length(keys)," spectra")
  for(k in seq_len(nrow(configs))) {
    cfg <- configs[k,]; tag <- cfg$configuration
    dfile <- file.path(dest,paste0("distance_",tag,".rds"))
    p <- p0; p$distance_method <- cfg$distance_method; p$ot_method <- cfg$ot_method
    if(!file.exists(dfile)) {
      gc(FALSE)
      elapsed <- system.time(d <- compute_distance_matrix(sp$frag_list,sp$loss_list,p,progress=FALSE))[["elapsed"]]
      stopifnot(all(is.finite(d)),identical(rownames(d),names(keys)),all(diag(d)==0))
      saveRDS(d,dfile)
      write_table(data.frame(dataset,configuration=tag,n=length(keys),cores=opts$cores,elapsed_seconds=elapsed),
                  file.path(dest,paste0("full_runtime_",tag,".csv")))
      message(dataset," ",tag," full matrix: ",round(elapsed,3)," s")
    } else d <- readRDS(dfile)
    stats_file <- file.path(dest,paste0("metrics_",tag,".csv"))
    if(!file.exists(stats_file)) {
      pq <- per_query_metrics_tie_aware(d,keys,tie_tolerance=0,random_seed=opts$seed)
      summary <- summarize_tie_aware_metrics(pq,dataset,tag,query_boot_R=opts$query_boot,
                            cluster_boot_R=opts$cluster_boot,seed=opts$seed)
      write_table(pq,file.path(dest,paste0("per_query_",tag,".csv")))
      write_table(summary,stats_file)
      pair <- ordered_pair_discrimination_metrics(d,keys)
      write_table(cbind(dataset,configuration=tag,pair$auc),file.path(dest,paste0("auc_",tag,".csv")))
      write_table(cbind(dataset,configuration=tag,pair$fdr),file.path(dest,paste0("fdr_",tag,".csv")))
      message(dataset," ",tag," retrieval statistics complete")
    }
  }
  e <- readRDS(file.path(dest,"distance_ppm_wasserstein.rds"))
  s <- readRDS(file.path(dest,"distance_ppm_wasserstein_sparse.rds"))
  stopifnot(max(abs(e-s))<1e-8)
  write_table(data.frame(dataset,cells=length(e),max_abs_error=max(abs(e-s)),
               mean_abs_error=mean(abs(e-s))),file.path(dest,"solver_agreement.csv"))
  sensitive <- list()
  for(tol in c(0,1e-12,1e-8)) for(tag in c("exact","exact_sparse")) {
    pq <- per_query_metrics_tie_aware(if(tag=="exact") e else s,keys,tie_tolerance=tol,random_seed=opts$seed)
    cols <- c("top1_fractional","rr_fractional","p_at_1_fractional","average_precision_fractional")
    sensitive[[length(sensitive)+1L]] <- data.frame(dataset,backend=tag,tolerance=tol,
      metric=cols,estimate=vapply(cols,function(nm) mean(pq[[nm]],na.rm=TRUE),numeric(1)))
  }
  write_table(do.call(rbind,sensitive),file.path(dest,"solver_tie_sensitivity.csv"))
  # All-method runtime comparison uses the same square query-library workloads.
  # Every configuration is measured separately; only within-run workers overlap.
  rtfile <- file.path(dest,"runtime_trials.csv")
  trials <- if(file.exists(rtfile)) read.csv(rtfile) else data.frame()
  for(n in sizes[sizes<=length(keys)]) {
    set.seed(opts$seed+n); idx <- sample(seq_along(keys),n)
    small <- subset_sp(sp,idx)
    write_table(data.frame(index=idx,id=names(keys)[idx]),file.path(dest,paste0("runtime_selection_",n,".csv")))
    jobs <- expand.grid(k=seq_len(nrow(configs)),cores=unique(c(1L,opts$cores)),replicate=seq_len(opts$replicates))
    set.seed(opts$seed+n); jobs <- jobs[sample(seq_len(nrow(jobs))),]
    for(j in seq_len(nrow(jobs))) {
      cfg <- configs[jobs$k[j],]; workers <- jobs$cores[j]; rep <- jobs$replicate[j]
      if(nrow(trials) && any(trials$n==n & trials$configuration==cfg$configuration & trials$cores==workers & trials$replicate==rep)) next
      p <- p0; p$distance_method <- cfg$distance_method; p$ot_method <- cfg$ot_method
      p$use_parallel <- workers>1; p$n_cores <- workers
      invisible(compute_distance_matrix_search(small$frag_list[1:2],small$loss_list[1:2],
                  small$frag_list[1:2],small$loss_list[1:2],p,progress=FALSE))
      gc(FALSE)
      elapsed <- system.time(d <- compute_distance_matrix_search(small$frag_list,small$loss_list,
                                   small$frag_list,small$loss_list,p,progress=FALSE))[["elapsed"]]
      stopifnot(all(is.finite(d)))
      matrix_file <- file.path(dest,paste0("runtime_reference_",n,"_",cfg$configuration,".rds"))
      if(file.exists(matrix_file)) stopifnot(max(abs(d-readRDS(matrix_file)))<1e-8) else saveRDS(d,matrix_file)
      trials <- rbind(trials,data.frame(dataset,n,cells=n*n,configuration=cfg$configuration,
                           cores=workers,replicate=rep,elapsed_seconds=elapsed))
      write_table(trials,rtfile)
      message(dataset," runtime n=",n," ",cfg$configuration," cores=",workers," rep=",rep," : ",round(elapsed,3)," s")
    }
    stopifnot(max(abs(readRDS(file.path(dest,paste0("runtime_reference_",n,"_ppm_wasserstein.rds")))-
                      readRDS(file.path(dest,paste0("runtime_reference_",n,"_ppm_wasserstein_sparse.rds")))))<1e-8)
  }
}
writeLines("All eight configurations completed and numerical checks passed.",file.path(opts$output_dir,"COMPLETE.txt"))
