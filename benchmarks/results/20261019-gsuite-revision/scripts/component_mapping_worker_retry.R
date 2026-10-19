args<-commandArgs(TRUE);stage<-args[1];out<-args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE);writeLines(as.character(Sys.getpid()),file.path(out,'pid.txt'))
threads<-if(stage%in%c('regional_twochain','genome_twochain'))1L else 4L
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_DYNAMIC='FALSE')
.libPaths(c('C:/Users/au223366/Documents/GitHub/gsuite/build/source-repository/library','C:/Users/au223366/Documents/GitHub/gsuite/build/r-library','C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(gbase);library(glma);library(gmap);library(gbayes);library(jsonlite)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(1)
repo<-'C:/Users/au223366/Documents/GitHub'
f<-readRDS(file.path(repo,'gsuite/build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds'))
G<-gbase::gprep(study='Frozen simulated-human tutorial',bedfiles=file.path(repo,'qgdata/simulated_human_data/human.bed'),
 bimfiles=file.path(repo,'gact/Document/glma-simulated-data/human.bim'),famfiles=file.path(repo,'qgdata/simulated_human_data/human.fam'),rsids=f$stat$rsids,nthreads=threads)
G$assembly<-'simulated-human-tutorial'
stopifnot(G$n==5000,length(unlist(G$rsids,use.names=FALSE))==37991)
regions<-lapply(f$sets,function(x)intersect(f$stat$rsids,x));stopifnot(identical(as.integer(lengths(regions)),c(850L,785L,1091L,1399L)))
ydata<-read.csv(file.path(repo,'gact/Document/glma-simulated-data/example-phenotypes.csv'),colClasses=c('character','numeric'))
y<-setNames(ydata$phenotype,ydata$id)
sets<-gbase::getChromosomeSets(G)
if(stage=='genome_twochain'){
 G<-readRDS(file.path(repo,'gsuite/build/benchmarks/manuscript-20261005-231349/cmp_suite_ld-1/output.rds'))
 gbase::validate_Glist(G)
 stopifnot(identical(unlist(G$rsids,use.names=FALSE),f$stat$rsids))
}
settings<-list(stage=stage,assembly=G$assembly,individuals=G$n,markers=37991,chromosomes=22,threads=threads)
if(stage=='regional_twochain')settings<-c(settings,list(regional_sizes=lengths(regions),nit=1000,nburn=100,nrun=2))
if(stage=='genome_twochain')settings<-c(settings,list(nit=100,nburn=0,nrun=2,max_ld_entries=1e8,LD='Reused validated sparse descriptor from manuscript-20261005-231349/cmp_suite_ld-1'))
if(grepl('^component',stage))settings<-c(settings,list(chromosome_sets=22,other_controls='Installed public defaults; non-MCMC model',public_defaults=capture.output(formals(glma::glma))))
write_json(settings,file.path(out,'settings.json'),auto_unbox=TRUE,pretty=TRUE)
warnings<-character();start<-proc.time()[3]
fit<-tryCatch(withCallingHandlers(switch(stage,
 regional_twochain=gmap::gmap(Glist=G,stat=f$stat,sets=regions,method='bayesC',nit=1000,nburn=100,nrun=2),
 genome_twochain=gbayes::gbayes(Glist=G,stat=f$stat,method='bayesC',nit=100,nburn=0,nrun=2,control=list(threads=1,max_ld_entries=1e8)),
 component_gaussian=glma::glma(y=y,Glist=G,sets=sets,method='infinitesimal_loco',threads=4),
 component_prior=glma::glma(y=y,Glist=G,sets=sets,method='prior_selected_loco',threads=4),stop('Unknown stage')),
 warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds<-unname(proc.time()[3]-start);ok<-!inherits(fit,'error')
serialize<-system.time(saveRDS(fit,file.path(out,'output.rds')))[3]
writeLines(capture.output(str(fit,max.level=2)),file.path(out,'structure.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
write_json(list(stage=stage,status=if(ok)'completed'else'failed',operation_seconds=seconds,serialization_seconds=unname(serialize),warnings=warnings,error=if(ok)NULL else conditionMessage(fit),output_bytes=file.info(file.path(out,'output.rds'))$size),file.path(out,'timing.json'),auto_unbox=TRUE,pretty=TRUE)
if(!ok)quit(status=1)
