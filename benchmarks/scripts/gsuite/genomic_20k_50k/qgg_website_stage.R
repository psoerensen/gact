args<-commandArgs(TRUE);stage<-args[1];threads<-as.integer(args[2]);out<-args[3];root<-dirname(out)
writeLines(as.character(Sys.getpid()),file.path(out,"r-process.pid"))
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",OMP_DYNAMIC="FALSE")
if(requireNamespace("RhpcBLASctl",quietly=TRUE))RhpcBLASctl::blas_set_num_threads(1)
library(qgg);library(gact)
writeLines(c(capture.output(sessionInfo()),paste("qgg",packageVersion("qgg"),find.package("qgg")),paste("gact",packageVersion("gact"),find.package("gact"))),file.path(out,"session.txt"))
load<-function(s)readRDS(file.path(root,paste0(s,"-1"),"output.rds"))
if(stage!="qgg_prepare")f<-load("qgg_prepare")
controls<-list(threads=1L,seed=20261004L)
if(stage=="qgg_mapping")controls<-c(controls,list(method="bayesR",algorithm="mcmc-eigen",nit=1000L,nburn=100L,nrun=1L,checkConvergence=FALSE,eigen_threshold=c(.995,.99,.95,.9)))
if(stage=="qgg_joint")controls<-c(controls,list(method="bayesC",pi=.01,nit=1000L,nburn=200L))
saveRDS(controls,file.path(out,"arguments.rds"));set.seed(20261004)
warn<-character();started<-Sys.time();clock<-proc.time()[3]
value<-tryCatch(withCallingHandlers(switch(stage,
 qgg_prepare={
  db<-normalizePath("build/benchmarks/manuscript-20261004-143135/database-1/database/hsa.0.0.1",winslash="/")
  GA<-readRDS("build/benchmarks/manuscript-20261004-143135/database-1/output.rds")
  panel<-read.delim("../qgdata/cache/1000g_phase3/integrated_call_samples_v3.20130502.ALL.panel")
  keep<-file.path(out,"eur.keep");writeLines(panel$sample[panel$super_pop=="EUR"],keep)
  prefix<-file.path(out,"eur_chr22")
  plink<-normalizePath("../qgdata/cache/tools/plink2_20260808/plink2.exe",winslash="/")
  vcf<-normalizePath("../qgdata/cache/1000g_phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",winslash="/")
  status<-system2(plink,c("--vcf",shQuote(vcf),"--keep",shQuote(keep),"--chr","22","--from-bp","20000000","--to-bp","22200000","--max-alleles","2","--snps-only","just-acgt","--maf","0.05","--make-bed","--threads","1","--out",shQuote(prefix)),stdout=file.path(out,"plink.stdout.log"),stderr=file.path(out,"plink.stderr.log"))
  stopifnot(status==0)
  G<-qgg::gprep(bedfiles=paste0(prefix,".bed"),assembly="GRCh37",ncores=1)
  stopifnot(G$n==503,!anyDuplicated(unlist(G$rsids)))
  stats<-setNames(lapply(c("GWAS1","GWAS2"),function(id){s<-gact::getMarkerStat(GA,studyID=id,rsids=unlist(G$rsids));qgg::checkStat(G,s,excludeMAF=.05,excludeMAFDIFF=.05,excludeINFO=.8,excludeCGAT=TRUE,excludeINDEL=TRUE,excludeDUPS=TRUE,excludeMHC=FALSE,excludeMISS=.05,excludeHWE=1e-12)}),c("GWAS1","GWAS2"))
  shared<-Reduce(intersect,lapply(stats,function(s)s$rsids));stats<-lapply(stats,function(s)s[match(shared,s$rsids),])
  annotation<-readRDS(file.path(db,"gsets/GRCh37_genesplus_annotation.rds"))
  genes<-annotation$EnsemblID[annotation$chr==22 & annotation$start>=20000000 & annotation$stop<=22200000]
  sets<-readRDS(file.path(db,"gsets/GRCh37.ensg2rsids.rds"));sets<-sets[intersect(genes,names(sets))]
  sets<-lapply(sets,function(ids)shared[shared %in% ids]);sets<-sets[lengths(sets)>=3];stopifnot(length(sets)>0)
  kegg<-gact::getFeatureSets(GA,feature="KEGG",minsets=50)
  gene_stat<-gact::getVEGAS(GA,studyID=c("GWAS1","GWAS2"))
  posterior<-gact::getMarkerStat(GA,studyID="BLR1")
  recipient_bed<-normalizePath("build/benchmarks/manuscript-20261004-161106/bed_conversion-1/genotypes.bed",winslash="/")
  recipient_bim<-read.table(sub("[.]bed$",".bim",recipient_bed),stringsAsFactors=FALSE)
  ref_bim<-read.table(paste0(prefix,".bim"),stringsAsFactors=FALSE)
  key<-function(b)paste(b[[1]],b[[4]],pmin(b[[5]],b[[6]]),pmax(b[[5]],b[[6]]),sep=":")
  at<-match(key(recipient_bim),key(ref_bim));hit<-!is.na(at);recipient_bim[[2]][hit]<-ref_bim[[2]][at[hit]]
  bim<-file.path(out,"recipients.bim");write.table(recipient_bim,bim,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  recipients<-qgg::gprep(bedfiles=recipient_bed,bimfiles=bim,famfiles=sub("[.]bed$",".fam",recipient_bed),assembly="GRCh37",ncores=1)
  posterior<-posterior[posterior$rsids %in% ref_bim[[2]],];stopifnot(nrow(posterior)>0)
  marker_sets<-lapply(kegg,function(gs)unique(unlist(sets[intersect(gs,names(sets))],use.names=FALSE)))
  marker_sets<-lapply(marker_sets,intersect,posterior$rsids);marker_sets<-marker_sets[lengths(marker_sets)>=2]
  stopifnot(length(marker_sets)>0)
  write.csv(data.frame(gene=names(sets),markers=lengths(sets)),file.path(out,"gene-manifest.csv"),row.names=FALSE)
  list(Glist=G,stat=stats,sets=sets,kegg=kegg,gene_stat=gene_stat,posterior=posterior,recipients=recipients,score_ids=tail(recipients$ids,2000),marker_sets=marker_sets,
   provenance=list(reference_n=G$n,reference_markers=sum(G$mchr),shared_markers=length(shared),gene_count=length(sets),gene_size=summary(lengths(sets)),scoring_markers=nrow(posterior),score_recipients=2000,pathway_scores=length(marker_sets),database=db,source="website workflows, regional chr22 restriction; cached whole-genome VEGAS and BLR1 inputs"))
 },
 qgg_vegas={p<-cbind(GWAS1=f$stat$GWAS1$p,GWAS2=f$stat$GWAS2$p);rownames(p)<-f$stat$GWAS1$rsids;qgg::vegas(Glist=f$Glist,sets=f$sets,p=p,verbose=FALSE)},
 qgg_joint=qgg::magma(stat=f$gene_stat,sets=f$kegg,method="bayesC",pi=.01,nit=1000,nburn=200),
 qgg_mapping=lapply(f$stat,function(s)qgg::gmap(Glist=f$Glist,stat=s,sets=f$sets,method="bayesR",algorithm="mcmc-eigen",eigen_threshold=c(.995,.99,.95,.9),nit=1000,nburn=100,nrun=1,checkConvergence=FALSE,verbose=FALSE,seed=20261004)),
 qgg_total=qgg::gscore(Glist=f$recipients,stat=f$posterior,ids=f$score_ids,ncores=1),
 qgg_pathway=qgg::gscore(Glist=f$recipients,stat=f$posterior,sets=f$marker_sets,ids=f$score_ids,ncores=1),
 stop("Unknown qgg stage")),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
seconds<-unname(proc.time()[3]-clock);failed<-inherits(value,"error");if(failed)warn<-c(warn,conditionMessage(value))
serialization<-unname(system.time(saveRDS(value,file.path(out,"output.rds")))[3])
write.csv(data.frame(stage,threads,started=format(started,tz="Europe/Copenhagen"),analysis_seconds=seconds,serialization_seconds=serialization,status=if(failed)"failed"else"completed",warnings=paste(unique(warn),collapse="; ")),file.path(out,"timing.csv"),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,"output-structure.txt"));if(failed)quit(status=1)
