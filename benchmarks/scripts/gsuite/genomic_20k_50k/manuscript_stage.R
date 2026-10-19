args <- commandArgs(TRUE)
stage <- args[1]; threads <- as.integer(args[2]); out <- args[3]
writeLines(as.character(Sys.getpid()), file.path(out,"r-process.pid"))
.libPaths(c("build/task-packages/library","build/task-packages/ecosystem/library","build/r-library",.libPaths()))
Sys.setenv(OMP_NUM_THREADS=threads, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
if(requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
packages <- if(stage=="database") "gact" else c("gbase","glma","gsea","gmap","gscore","gbayes")
for(p in packages) stopifnot(requireNamespace(p,quietly=TRUE))
writeLines(c(capture.output(sessionInfo()), unlist(lapply(packages,function(p) paste(p,packageVersion(p),find.package(p))))),file.path(out,"session.txt"))
if(stage!="database") {
  fixture <- "build/benchmarks/genomic-20k-50k"
  load <- function(n) readRDS(file.path(fixture,paste0(n,".rds")))$value
  simulation <- load("data_simulation"); LD <- load("sparse_ld_and_scores")
  stat <- load("gwas_and_checkstat")$statistics[1:2]
  markers <- sprintf("m%05d",1:50000)
  regions <- split(markers,rep(sprintf("region%03d",1:500),each=100))
  genes <- list()
  for(r in 1:500) {cut<-20+r%%60; nm<-sprintf("gene%04d",2*r-c(1,0));genes[[nm[1]]]<-regions[[r]][1:cut];genes[[nm[2]]]<-regions[[r]][(cut+1):100]}
  pathways <- list(enriched=names(genes)[1:160],background=names(genes)[401:560],mixed=names(genes)[seq(1,999,5)])
  if(stage=="joint") gene_stat<-gsea::gsea_input(load("vegas"))$stat[1:2]
  if(stage=="mapping") {
    pv<-setNames(sapply(names(stat),function(trait) var(simulation$simulation$Y[if(trait=="A") 1:6000 else 6001:12000,trait])),names(stat))
    prepared<-gmap::gmap_stat(stat,LD,phenotype_variance=pv)
  }
  if(stage %in% c("total","pathway")) {
    weights<-load("bayesian_regression")$weights
    recipients<-gbase::gprep(bedfiles=paste0(simulation$prefix,".bed"),ids=sprintf("id%05d",18001:20000),nthreads=threads)
    marker_sets<-lapply(pathways,function(x) unique(unlist(genes[x],use.names=FALSE)))
  }
}
controls <- switch(stage,
 vegas=list(tail_method="saddlepoint",threads=threads),
 joint=list(burnin=200L,sampling_sweeps=800L,seeds=c(11,29,47,71),threads=threads),
 mapping=list(burnin=500L,sampling_sweeps=700L,seeds=c(11,29,47,71),eigen=list(retained_positive_mass=c(.995,.99,.95))),
 list(threads=threads))
saveRDS(list(stage=stage,threads=threads,control=controls),file.path(out,"arguments.rds"))
warn<-character(); started<-Sys.time(); clock<-proc.time()[3]
value<-tryCatch(withCallingHandlers({
 switch(stage,
 database={x<-gact::gact(version="hsa.0.0.1",dbdir=file.path(out,"database"),task="download");x},
 vegas={metadata<-setNames(lapply(names(stat),function(t) data.frame(set=names(genes),sample_size=6000,mean_minor_allele_count=sapply(genes,function(ids) {p<-simulation$frequencies[match(ids,markers)];mean(12000*pmin(p,1-p))}))),names(stat));glma::vegas(stat,LD,sets=genes,metadata=metadata,control=list(tail_method="saddlepoint"),threads=threads)},
 joint={cov<-diag(2);dimnames(cov)<-list(names(gene_stat),names(gene_stat)); patterns<-rbind(null=c(0,0),A=c(1,0),B=c(0,1),both=c(1,1));colnames(patterns)<-names(gene_stat);gsea::gsea(gene_stat,pathways,"bayesc",sampling=list(gene_covariance="independent",trait_covariance=cov),prior=list(effect_covariance=cov,patterns=patterns,pattern_weights=c(null=.985,A=.005,B=.005,both=.005)),control=controls)},
 mapping=gmap::gmap(prepared,LD,regions=regions[1:20],method="bayesr",prior=list(residual_variance=1,effect_variance=.04,weights=c(.7,.1,.1,.1),variance_multipliers=c(0,.01,.1,1)),control=controls),
 total=gscore::gscore(weights,recipients,coding="standardized",scale_scores=FALSE,control=list(threads=threads)),
 pathway=gscore::gscore(weights,recipients,sets=marker_sets,include_total=FALSE,coding="standardized",scale_scores=FALSE,control=list(threads=threads)),
 stop("Unknown stage"))
 },warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
seconds<-unname(proc.time()[3]-clock)
failed<-inherits(value,"error")
if(stage=="database" && !failed) {
 if(!is.null(value$downloads)) write.csv(value$downloads,file.path(out,"downloads.csv"),row.names=FALSE)
 database_files<-list.files(value$dbdir,recursive=TRUE,full.names=TRUE); write.csv(data.frame(file=database_files,bytes=file.info(database_files)$size),file.path(out,"database-files.csv"),row.names=FALSE)
 if(is.null(value$downloads) || any(value$downloads$status != "success")) {warn<-c(warn,"Database installation incomplete: inspect downloads.csv"); failed<-TRUE}
}
if(inherits(value,"error")) warn<-c(warn,conditionMessage(value))
serialization<-system.time(saveRDS(value,file.path(out,"output.rds")))[3]
write.csv(data.frame(stage,threads,started=format(started,tz="Europe/Copenhagen"),analysis_seconds=seconds,serialization_seconds=unname(serialization),status=if(failed) "failed" else "completed",warnings=paste(unique(warn),collapse="; ")),file.path(out,"timing.csv"),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,"output-structure.txt"))
if(failed) quit(status=1)





