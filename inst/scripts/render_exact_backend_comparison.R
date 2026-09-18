#!/usr/bin/env Rscript
# Render previously completed serial measurements; never starts a benchmark.
# Args: existing comparison directory (RECETOX/, HREI/), output directory.
args <- commandArgs(TRUE)
if(length(args)!=2L) stop("Usage: Rscript render_exact_backend_comparison.R INPUT_DIR OUTPUT_DIR")
input <- args[1]; output <- args[2]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
library(ggplot2)
times <- list(); pairs <- list(); checks <- list(); provenance <- list()
for(ds in c("RECETOX","HREI")) {
  base <- file.path(input,ds)
  t <- read.csv(file.path(base,"timings.csv")); t$dataset <- ds; times[[ds]] <- t
  d <- readRDS(file.path(base,"distances.rds"))
  stopifnot(identical(dimnames(d$exact),dimnames(d$exact_sparse)))
  pairs[[ds]] <- data.frame(dataset=ds,query=rep(rownames(d$exact),ncol(d$exact)),
    candidate=rep(colnames(d$exact),each=nrow(d$exact)),
    exact=as.numeric(d$exact),exact_sparse=as.numeric(d$exact_sparse))
  a <- read.csv(file.path(base,"agreement.csv")); a$dataset <- ds; checks[[ds]] <- a
  cfg <- readRDS(file.path(base,"settings.rds"))
  provenance[[ds]] <- data.frame(dataset=ds,input=basename(cfg$options$input),
    input_md5=unname(cfg$input_md5),seed=cfg$options$seed,
    n_query=cfg$options$n_query,n_library=cfg$options$n_library,
    tol_ppm=cfg$params$tol_ppm,transition_mult=cfg$params$wasserstein_transition_mult,
    w_frag=cfg$params$w_frag,w_loss=cfg$params$w_loss,
    use_parallel=cfg$params$use_parallel,replicates=cfg$options$replicates)
}
times <- do.call(rbind,times); pairs <- do.call(rbind,pairs); checks <- do.call(rbind,checks)
stopifnot(all(is.finite(pairs$exact)),all(is.finite(pairs$exact_sparse)),
          max(abs(pairs$exact-pairs$exact_sparse))<1e-8,
          all(checks$top1_sets_agree),all(checks$rank_inversions==0))
summary <- do.call(rbind,lapply(split(times,interaction(times$dataset,times$backend)),function(x)
  data.frame(dataset=x$dataset[1],backend=x$backend[1],
             median_seconds=median(x$elapsed_seconds),min_seconds=min(x$elapsed_seconds),
             max_seconds=max(x$elapsed_seconds),replicates=nrow(x))))
speedup <- vapply(split(summary,summary$dataset),function(x)
  x$median_seconds[x$backend=="exact"]/x$median_seconds[x$backend=="exact_sparse"],numeric(1))
summary$label <- sprintf("%.3f s",summary$median_seconds)
is_sparse <- summary$backend=="exact_sparse"
summary$label[is_sparse] <- sprintf("%.3f s  |  %.2fx faster",summary$median_seconds[is_sparse],speedup[summary$dataset[is_sparse]])
summary$backend <- factor(summary$backend,levels=c("exact_sparse","exact"))
summary$dataset <- factor(summary$dataset,levels=c("RECETOX","HREI"),labels=c("RECETOX","HREI-MSDB"))
times$dataset <- factor(times$dataset,levels=c("RECETOX","HREI"),labels=c("RECETOX","HREI-MSDB"))
times$backend <- factor(times$backend,levels=c("exact_sparse","exact"))
pairs$dataset <- factor(pairs$dataset,levels=c("RECETOX","HREI"),labels=c("RECETOX","HREI-MSDB"))
style <- theme_minimal(base_size=12)+theme(
  panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18),
  plot.subtitle=element_text(color="#475569",margin=margin(b=14)),
  strip.text=element_text(face="bold",size=13),
  plot.caption=element_text(hjust=0,color="#475569",size=10,margin=margin(t=14)),
  plot.margin=margin(18,22,16,18),legend.position="none")
runtime <- ggplot(summary,aes(x=median_seconds,y=backend,fill=backend))+
  geom_col(width=.48)+
  geom_errorbar(aes(xmin=min_seconds,xmax=max_seconds),orientation="y",width=.16,color="#0F172A")+
  geom_point(data=times,aes(x=elapsed_seconds,y=backend),inherit.aes=FALSE,shape=21,size=2,fill="white",color="#0F172A")+
  geom_text(aes(x=max_seconds+.07,label=label),hjust=0,size=3.6,color="#0F172A")+
  facet_wrap(~dataset,ncol=1)+
  scale_fill_manual(values=c(exact="#334E68",exact_sparse="#008577"))+
  scale_y_discrete(labels=c(exact="exact (existing)",exact_sparse="exact_sparse"))+
  scale_x_continuous(limits=c(0,max(summary$max_seconds)+1.15),expand=expansion(mult=c(0,.02)))+
  labs(title="Exact OT backend runtime",
       subtitle="20 queries x 40 candidates | fragment + derived mass differences | serial execution",
       x="Distance-matrix elapsed time (seconds)",y=NULL,
       caption="Bars: median of 3 runs; points and whiskers: observed runs and range (not confidence intervals).\nMSP loading and spectrum preparation excluded. Measurements from the completed opt-in comparison.")+style
err <- max(abs(pairs$exact-pairs$exact_sparse))
agreement <- ggplot(pairs,aes(exact,exact_sparse))+
  geom_abline(slope=1,intercept=0,color="#94A3B8",linewidth=.7,linetype="dashed")+
  geom_point(color="#008577",size=1.7,alpha=.32)+
  facet_wrap(~dataset,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=seq(0,1,.25))+scale_y_continuous(breaks=seq(0,1,.25))+
  labs(title="Distance agreement between exact backends",
       subtitle=sprintf("1,600 query-candidate pairs | maximum absolute difference: %.2e",err),
       x="Distance: exact (existing)",y="Distance: exact_sparse",
       caption="Dashed line: equality. All 40 query Top-1 tie sets agree (tolerance 1e-8).\nNo candidate-order inversions for distance differences larger than 1e-8; overlapping points are not independent replicates.")+style
for(ext in c("png","svg")) {
  ggsave(file.path(output,paste0("exact_backend_runtime.",ext)),runtime,width=10,height=6.5,dpi=200,bg="white")
  ggsave(file.path(output,paste0("exact_backend_agreement.",ext)),agreement,width=10,height=6,dpi=200,bg="white")
}
write.csv(summary,file.path(output,"runtime_summary.csv"),row.names=FALSE)
write.csv(times,file.path(output,"timings.csv"),row.names=FALSE)
write.csv(pairs,file.path(output,"distance_pairs.csv"),row.names=FALSE)
write.csv(checks,file.path(output,"agreement.csv"),row.names=FALSE)
write.csv(do.call(rbind,provenance),file.path(output,"provenance.csv"),row.names=FALSE)
cat("Rendered completed serial comparison. Maximum error:",err,"\n")
