#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
input <- args[1]; output <- args[2]
stopifnot(file.exists(file.path(input,"COMPLETE.txt")))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
library(ggplot2)
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"lib/publication_retrieval_statistics.R"))
bind_files <- function(pattern) do.call(rbind,lapply(list.files(input,pattern=pattern,full.names=TRUE,recursive=TRUE),read.csv))
metrics <- bind_files("^metrics_.*[.]csv$")
runtime <- bind_files("^runtime_trials[.]csv$")
auc <- bind_files("^auc_.*[.]csv$")
agreement <- bind_files("^solver_agreement[.]csv$")
ties <- bind_files("^solver_tie_sensitivity[.]csv$")
order <- c("ppm_wasserstein","ppm_wasserstein_sparse","composite","entropy_weighted","entropy_unweighted","cosine","weighted_cosine","hellinger")
labels <- c("ppm-Wasserstein / exact","ppm-Wasserstein / sparse exact","Composite","Entropy (weighted)","Entropy (unweighted)","Cosine","Weighted cosine","Hellinger")
paired <- list()
for(dataset in c("RECETOX","HREI-MSDB")) {
  folder <- file.path(input,gsub("-","_",dataset))
  ref <- read.csv(file.path(folder,"per_query_ppm_wasserstein.csv"))
  for(comparator in setdiff(order,"ppm_wasserstein")) {
    cmp <- read.csv(file.path(folder,paste0("per_query_",comparator,".csv")))
    for(metric in c("top1_fractional","rr_fractional","p_at_1_fractional")) {
      cluster <- if(metric=="p_at_1_fractional") "inchikey_prefix" else "inchikey"
      diff <- paired_cluster_bootstrap_difference(ref,cmp,metric,cluster,R=20000,seed=20260908)
      paired[[length(paired)+1L]] <- data.frame(dataset,reference="ppm_wasserstein",comparator,metric,
        difference=diff[["estimate"]],ci_low=diff[["ci_low"]],ci_high=diff[["ci_high"]],
        n_queries=diff[["n_queries"]],n_clusters=diff[["n_clusters"]],cluster_unit=cluster)
    }
  }
}
write.csv(do.call(rbind,paired),file.path(output,"paired_cluster_differences.csv"),row.names=FALSE)
colors <- setNames(c("#334E68","#008577","#9966AA","#2274A5","#56B4E9","#888888","#A56F40","#D55E00"),order)
style <- theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold",size=17),plot.subtitle=element_text(color="#475569"),
  strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9),legend.position="bottom")
m <- subset(metrics,metric %in% c("top1","mrr","map"))
m$method <- factor(m$method,levels=rev(order))
m$metric <- factor(m$metric,levels=c("top1","mrr","map"),labels=c("Replicate Top-1","Replicate MRR","Connectivity mAP"))
p <- ggplot(m,aes(x=estimate,y=method,color=method))+
  geom_errorbar(aes(xmin=cluster_ci_low,xmax=cluster_ci_high),orientation="y",width=.2)+geom_point(size=2.3)+
  facet_grid(dataset~metric)+scale_color_manual(values=colors,labels=setNames(labels,order))+
  scale_y_discrete(labels=setNames(labels,order))+scale_x_continuous(limits=c(0,1))+
  labs(title="Full-library retrieval: eight configurations",subtitle="Seven distances; sparse exact is an alternative solver for ppm-Wasserstein",
       x="Fractional expected score",y=NULL,color=NULL,
       caption="All eligible spectra; self matches excluded. Error bars: 95% clustered bootstrap intervals (20,000 resamples).\nPrimary tie definition: exact floating-point equality. Solver tie-tolerance sensitivity is reported separately.")+style+theme(legend.position="none")
rt <- do.call(rbind,lapply(split(runtime,interaction(runtime$dataset,runtime$n,runtime$configuration,runtime$cores)),function(x)
  data.frame(dataset=x$dataset[1],n=x$n[1],configuration=x$configuration[1],cores=x$cores[1],
             median=median(x$elapsed_seconds),min=min(x$elapsed_seconds),max=max(x$elapsed_seconds))))
rt$configuration <- factor(rt$configuration,levels=order)
rt$execution <- ifelse(rt$cores==1,"Serial",paste0(rt$cores," workers"))
q <- ggplot(rt,aes(n,median,color=configuration,group=configuration))+
  geom_line(linewidth=.8)+geom_point(size=2)+geom_errorbar(aes(ymin=pmax(min,1e-4),ymax=pmax(max,1e-4)),width=1)+
  facet_grid(dataset~execution)+scale_y_log10()+scale_x_continuous(breaks=sort(unique(rt$n)))+
  scale_color_manual(values=colors,labels=labels)+
  labs(title="Runtime across all comparison configurations",subtitle="Identical n x n query-library workloads; configurations timed sequentially",
       x="Spectra per side (n)",y="Elapsed seconds (log scale)",color=NULL,
       caption="Points: median of 3 runs; whiskers: observed min/max, not confidence intervals.\nPreprocessed spectra; timing excludes MSP loading and representation generation.")+style+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
for(ext in c("png","svg")) {
  ggsave(file.path(output,paste0("full_retrieval.",ext)),p,width=13,height=8,dpi=180,bg="white")
  ggsave(file.path(output,paste0("full_runtime.",ext)),q,width=12,height=8,dpi=180,bg="white")
}
for(nm in c("metrics","runtime","auc","agreement","ties","rt")) write.csv(get(nm),file.path(output,paste0(nm,".csv")),row.names=FALSE)
cat("Rendered full-library retrieval and all-method runtime figures.\n")
