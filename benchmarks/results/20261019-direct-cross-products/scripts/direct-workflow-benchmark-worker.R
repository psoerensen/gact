args<-commandArgs(TRUE);stage<-args[1];out<-args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE);writeLines(as.character(Sys.getpid()),file.path(out,'pid.txt'))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
.libPaths(c('C:/Users/au223366/Documents/GitHub/gsuite/build/task-packages/library','C:/Users/au223366/Documents/GitHub/gsuite/build/source-repository/library','C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(gsea);library(jsonlite)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(1L)
t0<-proc.time()[['elapsed']];f<-readRDS('work/revision-benchmark/input-preparation/output.rds');loading<-proc.time()[['elapsed']]-t0
stopifnot(nrow(f$stat)==17955,length(f$sets)==1324,length(unique(f$chromosome))==22,
 as.character(packageVersion('gsea'))=='0.1.2',as.character(packageVersion('gbayes'))=='0.1.6')
stat<-data.frame(gene=rownames(f$stat),chromosome=as.character(f$chromosome),score=f$stat[,1])
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
settings<-list(stage=stage,genes=nrow(f$stat),features=length(f$sets),threads=1,chains=2,seed=10,
 sampling_sweeps=if(stage=='magma_joint')1000 else 5000,burnin=if(stage=='magma_joint')200 else 1000,
 maximum_preparation_bytes=1536*1024^2,gene_statistics='Saved GWAS1/GWAS2 scores; no recomputation')
write_json(settings,file.path(out,'settings.json'),pretty=TRUE,auto_unbox=TRUE)
ctrl<-list(threads=1,maximum_preparation_bytes=1536*1024^2)
# Process-local tracing supplies progress after each fold. Numeric methods,
# priors and inputs are unchanged; overhead is included in workflow wall time.
if(grepl('validation',stage)) {
 bench_fold<-0L
 trace('gsea',where=asNamespace('gsea'),print=FALSE,exit=quote({
  bench_result<-returnValue();bench_fold<<-bench_fold+1L
  cat('FOLD_COMPLETED',bench_fold,'PREP',bench_result$timing$preparation_seconds,
   'FIT',bench_result$timing$blr_fit_seconds,'\n');flush.console()
  jsonlite::write_json(list(completed_folds=bench_fold),file.path(out,'progress.json'),auto_unbox=TRUE)
 }))
}
warn<-character();start<-proc.time()[['elapsed']]
fit<-tryCatch(withCallingHandlers({
 if(stage=='magma_single') magma(f$stat[,1,drop=FALSE],f$sets,method='bayesC',pi=.01,nit=5000,nburn=1000,nrun=2,seed=10,control=ctrl)
 else if(stage=='magma_joint') magma(f$stat,f$sets,method='bayesC',pi=.01,nit=1000,nburn=200,nrun=2,seed=10,control=c(ctrl,list(memory_limit_bytes=1536*1024^2)))
 else if(grepl('validation',stage)) pops(stat,sets=f$sets,method=if(grepl('_c$',stage))'bayesC'else'bayesR',pi=.001,nit=5000,nburn=1000,nrun=2,seed=10,control=ctrl)
 else {
  method<-if(grepl('_c$',stage))'bayesc'else'bayesr'
  p<-list(residual_variance=.5)
  if(method=='bayesc'){p$inclusion_probability<-.001;mass<-.001;components<-2L}
  else {p$weights<-c(.999,rep(.001/3,3));p$variance_multipliers<-c(0,.01,.1,1);mass<-sum(p$weights*p$variance_multipliers);components<-4L}
  p$effect_variance<-.5/(length(f$sets)*mass)
  p$effect_variance_prior<-list(df=4,scale=p$effect_variance)
  p$residual_variance_prior<-list(df=4,scale=.5);p$weight_prior<-rep(1,components)
  gsea(stat,f$sets,method,sampling=list(gene_covariance='independent'),prior=p,
   control=c(ctrl,list(sampling_sweeps=5000,burnin=1000,thinning=1,seeds=c(10,11),
    estimate_effect_variance=TRUE,estimate_residual_variance=TRUE,estimate_weights=TRUE)))
 }
},warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds<-proc.time()[['elapsed']]-start;ok<-!inherits(fit,'error')
if(ok){
 stopifnot(identical(fit$timing$gene_statistics$computed,FALSE),is.na(fit$timing$gene_statistics$seconds))
 if(grepl('validation',stage))stopifnot(nrow(fit$folds)==22,nrow(fit$predictions)==17955,all(is.finite(fit$predictions$score)),all(fit$timing$folds$preparation_route=='direct_set_cross_products'))
 else stopifnot(all(is.finite(fit$pathways$mean)),nrow(fit$pathways)==if(stage=='magma_joint')2648 else 1324)
 write_json(fit$timing,file.path(out,'stages.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
 if(!is.null(fit$timing$folds))write.csv(fit$timing$folds,file.path(out,'fold-timings.csv'),row.names=FALSE)
 if(!is.null(fit$predictions))write.csv(fit$predictions,file.path(out,'predictions.csv'),row.names=FALSE)
 write_json(list(prior=fit$prior,control=fit$control),file.path(out,'effective-settings.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
}
serialization<-system.time(saveRDS(fit,file.path(out,'output.rds')))[['elapsed']]
write_json(list(stage=stage,status=if(ok)'completed'else'failed',operation_seconds=seconds,input_loading_seconds=loading,
 serialization_seconds=serialization,output_bytes=file.info(file.path(out,'output.rds'))$size,warnings=unique(warn),
 error=if(ok)NULL else conditionMessage(fit)),file.path(out,'timing.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
if(!ok)quit(status=1)
cat('COMPLETED',stage,seconds,'\n')
