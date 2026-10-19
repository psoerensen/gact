args<-commandArgs(TRUE);stage<-args[1];threads<-as.integer(args[2]);out<-args[3];root<-dirname(out)
writeLines(as.character(Sys.getpid()),file.path(out,"r-process.pid"))
.libPaths(c("build/task-packages/library","build/task-packages/ecosystem/library","build/r-library",.libPaths()))
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",OMP_DYNAMIC="FALSE")
if(requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)
packages<-c("gsim","gbase","glma","gsea","gmap","gscore","gbayes")
for(p in packages)stopifnot(requireNamespace(p,quietly=TRUE))
writeLines(c(capture.output(sessionInfo()),sapply(packages,function(p)paste(p,packageVersion(p),find.package(p)))),file.path(out,"session.txt"))
load<-function(s) { if(s=="preparation" && nzchar(Sys.getenv("GSIM_BENCH_FIXTURE_RDS"))) readRDS(Sys.getenv("GSIM_BENCH_FIXTURE_RDS")) else if(s=="bed_conversion" && nzchar(Sys.getenv("GSIM_BENCH_COHORT_RDS"))) readRDS(Sys.getenv("GSIM_BENCH_COHORT_RDS")) else if(s=="reference_import" && nzchar(Sys.getenv("GSIM_BENCH_REFERENCE_RDS"))) readRDS(Sys.getenv("GSIM_BENCH_REFERENCE_RDS")) else readRDS(file.path(root,paste0(s,"-1"),"output.rds")) }
if(stage %in% c("vegas_real","joint_real","mapping_real","weights_real","total_real","pathway_real")) {
 fixture<-load("preparation");LD<-fixture$LD;LD$assembly<-"GRCh37";stat<-fixture$stat;genes<-fixture$genes;pathways<-fixture$pathways
}
if(stage=="joint_real") gene_stat<-gsea::gsea_input(load("vegas_real"))$stat
if(stage=="mapping_real")prepared<-gmap::gmap_stat(stat,LD,phenotype_variance=fixture$phenotype_variance)
if(stage %in% c("total_real","pathway_real")) {
 weights<-load("weights_real")$weights
 recipients<-gbase::gprep(bedfiles=fixture$bed,ids=fixture$individuals[18001:20000],rsids=fixture$markers,nthreads=1L)
 marker_sets<-lapply(pathways,function(x)unique(unlist(genes[x],use.names=FALSE)))
}
controls<-switch(stage,joint_real=list(burnin=200L,sampling_sweeps=800L,seeds=c(11,29,47,71),threads=1L),
 mapping_real=list(burnin=500L,sampling_sweeps=700L,seeds=c(11,29,47,71),eigen=list(retained_positive_mass=c(.995,.99,.95),verification_tolerance=1e-7)),
 weights_real=list(burnin=500L,sampling_sweeps=700L,seeds=c(11,29,47,71),threads=1L,residual_policy="fixed"),list(threads=1L))
saveRDS(list(stage=stage,threads=threads,control=controls,seed=20260930L),file.path(out,"arguments.rds"))
started<-Sys.time();clock<-proc.time()[3];warn<-character()
value<-tryCatch(withCallingHandlers(switch(stage,
 reference_import={
  cache<-"../qgdata/cache/1000g_phase3"
  vcf<-file.path(cache,"ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
  stopifnot(unname(tools::md5sum(vcf))=="afdf383a32fcef939ca8a468bef23759")
  map<-read.delim(gzfile(file.path(cache,"genetic_map_GRCh37_chr22.txt.gz")),check.names=FALSE)
  panel<-read.delim(file.path(cache,"integrated_call_samples_v3.20130502.ALL.panel"))
  samples<-panel$sample[panel$super_pop=="EUR"]
  reference<-gsim::gsim_import_vcf(vcf,data.frame(chromosome="22",base_pair_position=map[["Position(bp)"]],genetic_position_cm=map[["Map(cM)"]]),
   output=file.path(out,"reference"),samples=samples,chromosome="22",region=c(20000000,22200000),unsupported="skip",overwrite=FALSE)
  reference
 },
 founders={
  reference<-load("reference_import");pedigree<-gsim::gsim_pedigree(n_generations=2,animals_per_generation=10000,founder_generations=1,new_founder_probability=1,seed=20260930)
  donor_ids<-reference$sample_ids
  if(is.null(donor_ids)) donor_ids<-read.table(reference$paths$fam,stringsAsFactors=FALSE)[[2]]
  populations<-setNames(rep("EUR",length(donor_ids)),donor_ids)
  base<-gsim::gsim_simulate_founders(reference,founder_ids=as.character(pedigree$canonical_order[pedigree$canonical_order %in% pedigree$pedigree$animal[is.na(pedigree$pedigree$sire)&is.na(pedigree$pedigree$dam)]]),populations=populations,
   ancestry_weights=c(EUR=1),mutation_age=setNames(rep(1e9,length(reference$variant_ids)),reference$variant_ids),
   N=c(EUR=length(donor_ids)),Ne=c(EUR=10000),rho=c(EUR=.02),seed=20260930,output=file.path(out,"founders"),batch_size=4096,threads=1L,overwrite=FALSE)
  list(base=base,pedigree=pedigree)
 },
 bed_conversion={x<-load("founders");gsim::gsim_simulate_pedigree(x$base,x$pedigree,seed=20260931,output=file.path(out,"genotypes"),format="bed",overwrite=FALSE)},
 preparation={
  base<-load("bed_conversion");bed<-base$paths[["bed"]]
  full<-gbase::gprep(bedfiles=bed,nthreads=1L)
  maf<-unlist(full$maf,use.names=FALSE);all_markers<-unlist(full$rsids,use.names=FALSE)
  markers<-all_markers[is.finite(maf)&maf>=.01];stopifnot(length(markers)>1274,length(full$ids)==20000)
  individuals<-full$ids
  glist<-gbase::gprep(bedfiles=bed,rsids=markers,nthreads=1L)
  set.seed(20260930);m<-length(markers);effects<-matrix(0,m,2,dimnames=list(markers,c("A","B")))
  active<-runif(m)<.02;effects[active,]<-matrix(rnorm(sum(active)*2),sum(active),2)
  phenotype<-gsim::gsim(Glist=glist,nt=2,architecture="fixed",beta=effects,h2=c(.3,.4),vg=c(.3,.4),re=0,standardize_W=FALSE,scale_effects=TRUE,seed=20260931,compute_sumstats=FALSE,chunk_size=64)
  colnames(phenotype$Y)<-c("A","B")
  ref<-gbase::gprep(bedfiles=bed,rsids=markers,ids=individuals[1:18000],nthreads=1L)
  ref$assembly<-"GRCh37"
  LD<-gbase::gprep(Glist=ref,task="sparseld",out_prefix=file.path(out,"LD"),max_distance_bp=0,max_distance_variants=273,r2=0,nthreads=1L,block_size=128,overwrite=FALSE)
  md<-LD$sparseLD$marker_metadata;statistics<-list();pv<-c()
  for(j in 1:2) {
   rows<-if(j==1)1:6000 else 6001:12000;trait<-c("A","B")[j]
   selected<-gbase::gprep(bedfiles=bed,rsids=markers,ids=individuals[rows],nthreads=1L)
   y<-setNames(phenotype$Y[,j],individuals)
   association<-glma::glma(y,selected,method="linear",threads=1L,block_size=128)
   a<-association$associations
   raw<-data.frame(rsids=markers,chr=md$chromosome,pos=md$base_pair_position,ea=md$allele1,nea=md$allele2,b=a$beta,seb=a$se,p=a$p,n=6000)
   statistics[[trait]]<-suppressMessages(gbase::checkStat(LD,raw,excludeMAF=NULL,excludeMAFDIFF=NULL,excludeINFO=NULL,excludeCGAT=FALSE,excludeINDEL=FALSE,excludeDUPS=TRUE,excludeMHC=FALSE,excludeMISS=NULL,excludeHWE=NULL))
   stopifnot(nrow(statistics[[trait]])==m);pv[trait]<-var(phenotype$Y[rows,j])
  }
  starts<-round(seq(1,m-273,length.out=1000));stopifnot(!anyDuplicated(starts))
  genes<-setNames(lapply(starts,function(start)markers[start+0:273]),sprintf("gene%04d",1:1000))
  pathways<-list(enriched=names(genes)[1:160],background=names(genes)[401:560],mixed=names(genes)[seq(1,999,5)])
  region_starts<-round(seq(1,m-99,length.out=20));regions<-setNames(lapply(region_starts,function(start)markers[start+0:99]),sprintf("region%03d",1:20))
  metadata<-setNames(lapply(names(statistics),function(trait)data.frame(set=names(genes),sample_size=6000,mean_minor_allele_count=sapply(genes,function(ids)mean(12000*maf[match(ids,all_markers)])))),names(statistics))
  list(LD=LD,stat=statistics,genes=genes,pathways=pathways,regions=regions,metadata=metadata,bed=bed,markers=markers,individuals=individuals,phenotype_variance=pv,
   provenance=list(source="real 1000G EUR phased chr22; gsim reference-conditioned founders",region="22:20000000-22200000 GRCh37",raw_markers=length(all_markers),retained_markers=m,maf_min=.01,n=20000,reference_n=18000,gene_count=1000,gene_size=274,gene_definition="size-controlled contiguous windows, not biological gene annotations",seed=20260930))
 },
 vegas_real=glma::vegas(stat,LD,sets=genes,metadata=fixture$metadata,control=list(tail_method="saddlepoint",relative_psd_tolerance=1e-7),threads=1L),
 joint_real={cov<-diag(2);dimnames(cov)<-list(names(gene_stat),names(gene_stat));patterns<-rbind(null=c(0,0),A=c(1,0),B=c(0,1),both=c(1,1));colnames(patterns)<-names(gene_stat)
  gsea::gsea(gene_stat,pathways,"bayesc",sampling=list(gene_covariance="independent",trait_covariance=cov),prior=list(effect_covariance=cov,patterns=patterns,pattern_weights=c(null=.985,A=.005,B=.005,both=.005)),control=controls)},
 mapping_real=gmap::gmap(prepared,LD,regions=fixture$regions,method="bayesr",prior=list(residual_variance=1,effect_variance=.04,weights=c(.7,.1,.1,.1),variance_multipliers=c(0,.01,.1,1)),control=controls),
 weights_real={
  # Scoring measures weight application, independently of posterior fitting.
  # Recreate fixed simulated effects; do not reuse a failed genome-wide fit.
  set.seed(20260930);m<-length(fixture$markers);active<-runif(m)<.02
  effects<-matrix(0,m,2);effects[active,]<-matrix(rnorm(sum(active)*2),sum(active),2)
  md<-LD$sparseLD$marker_metadata
  list(weights=data.frame(marker=md$marker_id,chromosome=md$chromosome,position=md$base_pair_position,
   effect_allele=md$allele1,other_allele=md$allele2,allele_frequency=md$allele_frequency,
   weight=effects[,1]/sqrt(sum(effects[,1]^2))),provenance="Fixed sparse simulation-effect weights, unit L2 norm; not posterior weights")
 },
 total_real=gscore::gscore(weights,recipients,coding="standardized",scale_scores=FALSE,control=list(threads=1L)),
 pathway_real=gscore::gscore(weights,recipients,sets=marker_sets,include_total=FALSE,coding="standardized",scale_scores=FALSE,control=list(threads=1L)),
 stop("Unknown real-LD stage")),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
seconds<-unname(proc.time()[3]-clock);failed<-inherits(value,"error")
if(failed)warn<-c(warn,conditionMessage(value))
serialization<-unname(system.time(saveRDS(value,file.path(out,"output.rds")))[3])
write.csv(data.frame(stage,threads,started=format(started,tz="Europe/Copenhagen"),analysis_seconds=seconds,serialization_seconds=serialization,status=if(failed)"failed" else "completed",warnings=paste(unique(warn),collapse="; ")),file.path(out,"timing.csv"),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,"output-structure.txt"))
if(failed)quit(status=1)

