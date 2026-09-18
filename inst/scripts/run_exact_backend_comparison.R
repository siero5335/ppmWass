#!/usr/bin/env Rscript
# Compare the two exact solvers on identical preprocessed search matrices.
# Example: Rscript run_exact_backend_comparison.R --input=library.msp --output-dir=results
opts <- list(input=NULL, output_dir="exact_backend_comparison", n_query=20L,
             n_library=40L, replicates=3L, seed=20260908L, library=NULL)
for (arg in commandArgs(TRUE)) {
  bits <- strsplit(sub("^--", "", arg), "=", fixed=TRUE)[[1]]
  key <- gsub("-", "_", bits[1], fixed=TRUE)
  if (length(bits)!=2L || !key %in% names(opts)) stop("Unknown argument: ",arg)
  opts[[key]] <- if(key %in% c("n_query","n_library","replicates","seed")) as.integer(bits[2]) else bits[2]
}
stopifnot(!is.null(opts$input), file.exists(opts$input),
          opts$n_query>0, opts$n_library>0, opts$replicates>0)
if (!is.null(opts$library)) .libPaths(c(normalizePath(opts$library),.libPaths()))
library(ppmWass)
p <- eihrms_default_params(); p$distance_method <- "ppm_wasserstein"
p$tol_ppm <- 15; p$use_parallel <- FALSE
sp <- build_spectra_from_msp(opts$input,p,progress=FALSE)
keep <- which(vapply(sp$frag_list,nrow,integer(1))>0)
stopifnot(length(keep)>=opts$n_query+opts$n_library)
set.seed(opts$seed)
idx <- sample(keep,opts$n_query+opts$n_library)
qi <- idx[seq_len(opts$n_query)]; li <- idx[opts$n_query+seq_len(opts$n_library)]
qf <- sp$frag_list[qi]; ql <- sp$loss_list[qi]
lf <- sp$frag_list[li]; ll <- sp$loss_list[li]
dir.create(opts$output_dir,recursive=TRUE,showWarnings=FALSE)
write.csv(data.frame(role=c(rep("query",length(qi)),rep("library",length(li))),
                     index=c(qi,li),id=names(sp$frag_list)[c(qi,li)]),
          file.path(opts$output_dir,"selected_spectra.csv"),row.names=FALSE)
saveRDS(list(options=opts,params=p,input_md5=tools::md5sum(opts$input)),file.path(opts$output_dir,"settings.rds"))
timings <- list(); distances <- list()
# Warm both backends before alternating their order in timed runs.
for(backend in c("exact","exact_sparse")) {
  p$ot_method <- backend
  invisible(compute_distance_matrix_search(qf[1],ql[1],lf[1],ll[1],p,progress=FALSE))
}
for(rep in seq_len(opts$replicates)) {
  order <- if(rep%%2) c("exact","exact_sparse") else c("exact_sparse","exact")
  for(backend in order) {
    p$ot_method <- backend; gc(FALSE)
    elapsed <- system.time(d <- compute_distance_matrix_search(qf,ql,lf,ll,p,progress=FALSE))[["elapsed"]]
    stopifnot(all(is.finite(d)))
    distances[[backend]] <- d
    timings[[length(timings)+1L]] <- data.frame(replicate=rep,backend,
       query_count=length(qi),library_count=length(li),elapsed_seconds=elapsed)
  }
  stopifnot(max(abs(distances$exact-distances$exact_sparse))<1e-8)
}
e <- distances$exact; s <- distances$exact_sparse
checks <- lapply(seq_len(nrow(e)),function(i) {
  de <- outer(e[i,],e[i,],"-"); ds <- outer(s[i,],s[i,],"-")
  decisive <- abs(de)>1e-8
  top1 <- identical(which(e[i,]<=min(e[i,])+1e-8),which(s[i,]<=min(s[i,])+1e-8))
  inversions <- sum(sign(de[decisive])!=sign(ds[decisive]))/2
  stopifnot(top1,inversions==0)
  data.frame(query=rownames(e)[i],max_abs_error=max(abs(e[i,]-s[i,])),top1_sets_agree=top1,rank_inversions=inversions)
})
timings <- do.call(rbind,timings)
write.csv(timings,file.path(opts$output_dir,"timings.csv"),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(opts$output_dir,"agreement.csv"),row.names=FALSE)
saveRDS(distances,file.path(opts$output_dir,"distances.rds"))
writeLines(capture.output(sessionInfo()),file.path(opts$output_dir,"sessionInfo.txt"))
print(aggregate(elapsed_seconds~backend,timings,median))
cat("Maximum absolute distance error:",max(abs(e-s)),"\n")
