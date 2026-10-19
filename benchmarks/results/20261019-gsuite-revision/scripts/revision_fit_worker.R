args<-commandArgs(TRUE);stage<-args[1];out<-args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines(as.character(Sys.getpid()),file.path(out,'pid.txt'))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
.libPaths(c('C:/Users/au223366/Documents/GitHub/gsuite/build/source-repository/library',
 'C:/Users/au223366/Documents/GitHub/gsuite/build/r-library',
 'C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(gsea);library(jsonlite)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(1L)
f<-readRDS('work/revision-benchmark/fixture.rds')
stat<-data.frame(gene=rownames(f$stat),chromosome=as.character(f$chromosome),score=f$stat[,1])
writeLines(c(capture.output(sessionInfo()),paste('gsea',packageVersion('gsea'),find.package('gsea'))),file.path(out,'session.txt'))
settings<-list(stage=stage,genes=nrow(f$stat),pathways=length(f$sets),annotation='Reactome from installed database',
 database='hsa.0.0.1 installed 4 October 2026',seed=10,threads=1,
 nit=if(stage=='gene_joint')1000 else 5000,nburn=if(stage=='gene_joint')200 else 1000,nrun=2,maximum_preparation_bytes=1536*1024^2)
write_json(settings,file.path(out,'settings.json'),auto_unbox=TRUE,pretty=TRUE)
ctrl<-list(threads=1,maximum_preparation_bytes=1536*1024^2)
jointctrl<-c(ctrl,list(memory_limit_bytes=1536*1024^2))
warn<-character();start<-proc.time()[3]
fit<-tryCatch(withCallingHandlers(switch(stage,
 gene_single=gsea::magma(stat=f$stat[,1,drop=FALSE],sets=f$sets,method='bayesC',pi=.01,nit=5000,nburn=1000,nrun=2,seed=10,control=ctrl),
 gene_joint=gsea::magma(stat=f$stat,sets=f$sets,method='bayesC',pi=.01,nit=1000,nburn=200,nrun=2,seed=10,control=jointctrl),
 pops_rr=gsea::pops(stat=stat,sets=f$sets,method='rr',threshold=.05),
 pops_c=gsea::pops(stat=stat,sets=f$sets,method='bayesC',pi=.001,nit=5000,nburn=1000,nrun=2,seed=10,control=ctrl),
 pops_r=gsea::pops(stat=stat,sets=f$sets,method='bayesR',pi=.001,nit=5000,nburn=1000,nrun=2,seed=10,control=ctrl),
 stop('Unknown stage')),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds<-unname(proc.time()[3]-start)
ok<-!inherits(fit,'error')
if(ok && grepl('^pops',stage))stopifnot(nrow(fit$predictions)==nrow(stat),all(is.finite(fit$predictions$score)),nrow(fit$folds)==22)
if(ok && grepl('^gene',stage))stopifnot(nrow(fit$pathways)>0)
serialize<-system.time(saveRDS(fit,file.path(out,'output.rds')))[3]
writeLines(capture.output(str(fit,max.level=2)),file.path(out,'structure.txt'))
if(ok && !is.null(fit$pathways))write.csv(fit$pathways,file.path(out,'pathways.csv'),row.names=FALSE)
if(ok && !is.null(fit$folds))write.csv(fit$folds,file.path(out,'folds.csv'),row.names=FALSE)
if(ok && !is.null(fit$predictions))write.csv(fit$predictions,file.path(out,'predictions.csv'),row.names=FALSE)
write_json(list(stage=stage,operation_seconds=seconds,serialization_seconds=unname(serialize),status=if(ok)'completed'else'failed',
 warnings=unique(warn),error=if(ok)NULL else conditionMessage(fit),output_bytes=file.info(file.path(out,'output.rds'))$size),
 file.path(out,'timing.json'),auto_unbox=TRUE,pretty=TRUE)
if(!ok)quit(status=1)
