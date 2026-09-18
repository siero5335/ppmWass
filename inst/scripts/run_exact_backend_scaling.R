#!/usr/bin/env Rscript
# Sequentially measure solver/worker configurations without competing jobs.
opts <- list(input_dir=NULL, output_dir="exact_backend_scaling", library=NULL,
             sizes="20,40,80", cores="1", replicates="3", seed="20260908")
for (arg in commandArgs(TRUE)) {
  eq <- regexpr("=",arg,fixed=TRUE)[1]
  if(eq<1) stop("Expected --key=value: ",arg)
  key <- gsub("-","_",sub("^--","",substr(arg,1,eq-1)),fixed=TRUE)
  if(!key %in% names(opts)) stop("Unknown argument: ",arg)
  opts[[key]] <- substring(arg,eq+1)
}
parse_ints <- function(x) {
  v <- suppressWarnings(as.numeric(strsplit(x,",",fixed=TRUE)[[1]]))
  if(!length(v) || any(!is.finite(v) | v<1 | v!=floor(v))) stop("Expected positive integers")
  as.integer(v)
}
sizes <- sort(unique(parse_ints(opts$sizes))); cores <- sort(unique(parse_ints(opts$cores)))
reps <- parse_ints(opts$replicates); seed <- parse_ints(opts$seed)
stopifnot(length(reps)==1,length(seed)==1,!is.null(opts$input_dir))
if(!is.null(opts$library)) .libPaths(c(normalizePath(opts$library),.libPaths()))
library(ppmWass)
dir.create(opts$output_dir,recursive=TRUE,showWarnings=FALSE)
p <- eihrms_default_params(); p$distance_method <- "ppm_wasserstein"; p$tol_ppm <- 15
timings <- list(); agreement <- list(); selection <- list(); provenance <- list()
flush_csv <- function(rows,name) write.csv(do.call(rbind,rows),file.path(opts$output_dir,name),row.names=FALSE)
saveRDS(list(options=opts,params=p,platform=Sys.info(),physical_cores=parallel::detectCores(FALSE),
             logical_cores=parallel::detectCores(TRUE)),file.path(opts$output_dir,"settings.rds"))
writeLines(capture.output(sessionInfo()),file.path(opts$output_dir,"sessionInfo.txt"))
for(dataset in c("RECETOX_merged","HREI-MSDB")) {
  input <- file.path(opts$input_dir,paste0(dataset,".msp"))
  stopifnot(file.exists(input))
  sp <- build_spectra_from_msp(input,p,progress=FALSE)
  keep <- which(vapply(sp$frag_list,nrow,integer(1))>0)
  stopifnot(length(keep)>=3*max(sizes))
  set.seed(seed); idx <- sample(keep,3*max(sizes))
  qpool <- idx[seq_len(max(sizes))]; lpool <- idx[max(sizes)+seq_len(2*max(sizes))]
  provenance[[dataset]] <- data.frame(dataset,input=basename(input),md5=unname(tools::md5sum(input)),
                                     eligible_spectra=length(keep),seed=seed)
  selection[[dataset]] <- data.frame(dataset,role=rep(c("query","library"),c(length(qpool),length(lpool))),
                      input_index=c(qpool,lpool),id=names(sp$frag_list)[c(qpool,lpool)])
  flush_csv(provenance,"inputs.csv"); flush_csv(selection,"selected_spectra.csv")
  for(nq in sizes) {
    qi <- qpool[seq_len(nq)]; li <- lpool[seq_len(2*nq)]
    qf <- sp$frag_list[qi]; ql <- sp$loss_list[qi]; lf <- sp$frag_list[li]; ll <- sp$loss_list[li]
    values <- list()
    configs <- expand.grid(backend=c("exact","exact_sparse"),cores=cores,stringsAsFactors=FALSE)
    for(k in seq_len(nrow(configs))) {
      p$ot_method <- configs$backend[k]; p$use_parallel <- configs$cores[k]>1; p$n_cores <- configs$cores[k]
      invisible(compute_distance_matrix_search(qf[1:2],ql[1:2],lf[1:2],ll[1:2],p,progress=FALSE))
    }
    for(rep in seq_len(reps)) {
      set.seed(seed+nq+rep)
      for(k in sample(seq_len(nrow(configs)))) {
        backend <- configs$backend[k]; workers <- configs$cores[k]
        p$ot_method <- backend; p$use_parallel <- workers>1; p$n_cores <- workers
        gc(FALSE)
        elapsed <- system.time(d <- compute_distance_matrix_search(qf,ql,lf,ll,p,progress=FALSE))[["elapsed"]]
        stopifnot(all(is.finite(d)))
        key <- paste(backend,workers,sep="_")
        if(!is.null(values[[key]])) stopifnot(max(abs(d-values[[key]]))<1e-8)
        values[[key]] <- d
        timings[[length(timings)+1L]] <- data.frame(dataset,n_query=nq,n_library=2*nq,cells=2*nq*nq,
                  backend,cores=workers,replicate=rep,elapsed_seconds=elapsed)
        flush_csv(timings,"timings.csv")
        message(dataset," ",nq,"x",2*nq," ",backend," cores=",workers," rep=",rep," : ",round(elapsed,3)," s")
      }
    }
    reference <- values[[paste("exact",min(cores),sep="_")]]
    for(key in names(values)) {
      d <- values[[key]]
      for(i in seq_len(nq)) {
        dr <- outer(reference[i,],reference[i,],"-"); dd <- outer(d[i,],d[i,],"-")
        decisive <- abs(dr)>1e-8
        inversions <- sum(sign(dr[decisive])!=sign(dd[decisive]))/2
        agree <- identical(which(reference[i,]<=min(reference[i,])+1e-8),which(d[i,]<=min(d[i,])+1e-8))
        error <- max(abs(reference[i,]-d[i,]))
        stopifnot(error<1e-8,agree,inversions==0)
        agreement[[length(agreement)+1L]] <- data.frame(dataset,n_query=nq,n_library=2*nq,
                  configuration=key,query=rownames(d)[i],max_abs_error=error,
                  top1_sets_agree=agree,rank_inversions=inversions)
      }
    }
    flush_csv(agreement,"agreement.csv")
    saveRDS(values,file.path(opts$output_dir,paste0(dataset,"_",nq,"_distances.rds")))
  }
}
cat("Scaling benchmark and all agreement checks passed.\n")
