args<-commandArgs(TRUE);stage<-args[1];out<-args[3];root<-dirname(out);threads<-as.integer(args[2]);stopifnot(threads>=1L)
writeLines(as.character(Sys.getpid()),file.path(out,"r-process.pid"))
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC="FALSE")
if(requireNamespace("RhpcBLASctl",quietly=TRUE))RhpcBLASctl::blas_set_num_threads(threads)
library(qgg);library(gact)
writeLines(c(capture.output(sessionInfo()),paste("qgg",packageVersion("qgg"),find.package("qgg")),paste("gact",packageVersion("gact"),find.package("gact"))),file.path(out,"session.txt"))
load<-function(s){p<-file.path(root,paste0(s,"-1"),"output.rds");if(!file.exists(p)&&nzchar(Sys.getenv("QGG_BENCH_REUSE_RUN")))p<-file.path(Sys.getenv("QGG_BENCH_REUSE_RUN"),paste0(s,"-1"),"output.rds");readRDS(p)}
db_descriptor<-"build/benchmarks/manuscript-20261004-143135/database-1/output.rds"
if(stage %in% c("qgg_glist","qgg_inputs"))GA<-readRDS(db_descriptor)
if(stage=="qgg_ld")G<-load("qgg_glist")
if(stage %in% c("qgg_vegas","qgg_joint","qgg_mapping","qgg_bayesc_fit","qgg_bayesr_fit","qgg_total_c","qgg_total_r","qgg_pathway"))f<-load("qgg_inputs")
if(stage=="qgg_total_c")fit<-load("qgg_bayesc_fit")
if(stage %in% c("qgg_total_r","qgg_pathway"))fit<-load("qgg_bayesr_fit")
set.seed(20261004)
writeLines(capture.output(lapply(c("vegas","gmap","magma","gbayes","gscore","gprep"),function(n)list(name=n,defaults=formals(getExportedValue("qgg",n))))),file.path(out,"installed-defaults.txt"))
warn<-character();started<-Sys.time();clock<-proc.time()[3]
value<-tryCatch(withCallingHandlers(switch(stage,
 qgg_download={
  GA<-readRDS(db_descriptor);GA$dirs["marker"]<-out
  GA$metadata$sources[["1000G"]]<-"10.5281/zenodo.14141438"
  GA<-gact::downloadDB(GAlist=GA,what="1000G")
  stopifnot(all(file.exists(file.path(out,paste0("g1000_eur",c(".bed",".bim",".fam"))))))
  write.csv(GA$downloads,file.path(out,"downloads.csv"),row.names=FALSE)
  list(prefix=file.path(out,"g1000_eur"),doi="10.5281/zenodo.14141438",ancestry="EUR")
 },
 qgg_glist={
  prefix<-load("qgg_download")$prefix;filtered<-file.path(out,"g1000_eur_filtered")
  qgg::writeBED(bedRead=paste0(prefix,".bed"),bimRead=paste0(prefix,".bim"),famRead=paste0(prefix,".fam"),bedWrite=paste0(filtered,".bed"),bimWrite=paste0(filtered,".bim"),famWrite=paste0(filtered,".fam"),rsids=GA$rsids)
  qgg::gprep(study="1000G EUR",bedfiles=paste0(filtered,".bed"),bimfiles=paste0(filtered,".bim"),famfiles=paste0(filtered,".fam"),ncores=1)
 },
 qgg_ld={source("benchmarks/genomic_20k_50k/qgg_ld_guard.R",local=TRUE);run_website_ld(G,out)},
 qgg_inputs={
  G<-load("qgg_ld")
  # Gene methods compute LD from BED directly; preserve the failed full-LD
  # record while allowing these independent published examples to proceed.
  if(inherits(G,"error"))G<-load("qgg_glist")
  fine_stat<-gact::getMarkerStat(GAlist=GA,studyID="GWAS1")
  fine_stat<-qgg::checkStat(Glist=G,stat=fine_stat,excludeMAF=.05,excludeMAFDIFF=.05,excludeINFO=.8,excludeCGAT=TRUE,excludeINDEL=TRUE,excludeDUPS=TRUE,excludeMHC=FALSE,excludeMISS=.05,excludeHWE=1e-12)
  vegas_stat<-setNames(lapply(c("GWAS1","GWAS2"),function(id)qgg::checkStat(Glist=G,stat=gact::getMarkerStat(GAlist=GA,studyID=id))),c("GWAS1","GWAS2"))
  marker_sets<-gact::getMarkerSets(GAlist=GA,feature="Genesplus")
  kegg<-gact::getFeatureSets(GAlist=GA,feature="KEGG",minsets=50)
  scoring_sets<-gact::getMarkerSets(GAlist=GA,feature="KEGG")
  gene_stat<-gact::getVEGAS(GAlist=GA,studyID=c("GWAS1","GWAS2"))
  genes<-c("ENSG00000099204","ENSG00000118972","ENSG00000145996","ENSG00000140718")
  stopifnot(all(genes %in% names(marker_sets)),nrow(fine_stat)>0)
  list(Glist=G,fine_stat=fine_stat,vegas_stat=vegas_stat,marker_sets=marker_sets,fine_sets=marker_sets[genes],kegg=kegg,scoring_sets=scoring_sets,gene_stat=gene_stat)
 },
 qgg_vegas={
  # The website's multi-study accessor returns a list; qgg checkStat needs a
  # data frame. Check each study, then use the public multi-trait p interface.
  shared<-Reduce(intersect,lapply(f$vegas_stat,function(s)s$rsids))
  p<-sapply(f$vegas_stat,function(s)s$p[match(shared,s$rsids)]);rownames(p)<-shared
  qgg::vegas(Glist=f$Glist,sets=f$marker_sets,p=p,verbose=FALSE)
 },
 qgg_joint=qgg::magma(stat=f$gene_stat,sets=f$kegg,method="bayesC",pi=.01,nit=1000,nburn=200),
 qgg_mapping=qgg::gmap(Glist=f$Glist,stat=f$fine_stat,sets=f$fine_sets,method="bayesR",algorithm="mcmc-eigen",eigen_threshold=c(.995,.99,.95,.9),verbose=FALSE),
 qgg_bayesc_fit=qgg::gbayes(stat=f$fine_stat,Glist=f$Glist,method="bayesC",nit=1000),
 qgg_bayesr_fit=qgg::gbayes(stat=f$fine_stat,Glist=f$Glist,method="bayesR",nit=5000,nburn=1000),
 qgg_total_c=qgg::gscore(Glist=f$Glist,stat=fit$stat,ncores=1),
 qgg_total_r=qgg::gscore(Glist=f$Glist,stat=fit$stat,ncores=1),
 qgg_pathway=qgg::gscore(Glist=f$Glist,stat=fit$stat,sets=f$scoring_sets,ncores=1),
 stop("Unknown website stage")),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
seconds<-unname(proc.time()[3]-clock);failed<-inherits(value,"error");if(failed)warn<-c(warn,conditionMessage(value))
serialization<-unname(system.time(saveRDS(value,file.path(out,"output.rds")))[3])
write.csv(data.frame(stage,threads,started=format(started,tz="Europe/Copenhagen"),analysis_seconds=seconds,serialization_seconds=serialization,status=if(failed)"failed"else"completed",warnings=paste(unique(warn),collapse="; ")),file.path(out,"timing.csv"),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,"output-structure.txt"));if(failed)quit(status=1)
