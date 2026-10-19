args<-commandArgs(TRUE);stage<-args[1];out<-args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines(as.character(Sys.getpid()),file.path(out,'pid.txt'))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
.libPaths(c('C:/Users/au223366/Documents/GitHub/gsuite/build/source-repository/library','C:/Users/au223366/Documents/GitHub/gsuite/build/r-library','C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(gsea);library(gbase);library(jsonlite)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(1L)
load_start<-proc.time()[3]
f<-readRDS('work/revision-benchmark/input-preparation/output.rds')
loading_seconds<-unname(proc.time()[3]-load_start)
prep_start<-proc.time()[3]
stat<-data.frame(gene=rownames(f$stat),score=f$stat[,1],stringsAsFactors=FALSE)
mapped<-gbase::mapSets(lapply(f$sets,unique),rsids=stat$gene)
mapped<-mapped[lengths(mapped)>0L]
x<-matrix(0,nrow(stat),length(mapped),dimnames=list(stat$gene,names(mapped)))
for(j in seq_along(mapped))x[mapped[[j]],j]<-1
stopifnot(nrow(x)==17955,ncol(x)==1324,all(colSums(x)>0),all(colSums(x)<nrow(x)),all(is.finite(stat$score)))
feature_preparation_seconds<-unname(proc.time()[3]-prep_start)
weights<-c(.999,rep(.001/3,3));multipliers<-c(0,.01,.1,1)
effect_variance<-.5/(ncol(x)*sum(weights*multipliers))
prior<-list(residual_variance=.5,effect_variance=effect_variance,weights=weights,variance_multipliers=multipliers,
 effect_variance_prior=list(df=4,scale=effect_variance),residual_variance_prior=list(df=4,scale=.5),weight_prior=rep(1,4))
control<-list(threads=1,burnin=1000,sampling_sweeps=5000,thinning=1,seeds=c(10,11),
 estimate_effect_variance=TRUE,estimate_residual_variance=TRUE,estimate_weights=TRUE,retain_diagnostics=TRUE,
 scale_response=TRUE,scale_features=TRUE,maximum_preparation_bytes=1536*1024^2)
settings<-list(stage=stage,api='gsea::gsea',method='bayesr',models=1,held_out_folds=0,genes=nrow(x),features=ncol(x),trait=colnames(f$stat)[1],
 annotation='Reactome from existing installed hsa.0.0.1 database',input='work/revision-benchmark/input-preparation/output.rds',
 gene_association_computation='Not run; existing cached scores reused',chains=2,prior=prior,control=control,
 fit_scope='One full-data BLR API call, including native normalization/design preparation, MCMC and result assembly; no held-out validation or prediction. External feature construction and loading separately timed.')
write_json(settings,file.path(out,'settings.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
warn<-character();start<-proc.time()[3]
fit<-tryCatch(withCallingHandlers(gsea::gsea(stat=stat,sets=x,method='bayesr',sampling=list(gene_covariance='independent'),prior=prior,control=control),
 warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds<-unname(proc.time()[3]-start);ok<-!inherits(fit,'error')
if(ok)stopifnot(nrow(fit$pathways)==1324,all(is.finite(fit$pathways$mean)),all(is.finite(fit$pathways$pip)))
serialize<-system.time(saveRDS(fit,file.path(out,'output.rds')))[3]
writeLines(capture.output(str(fit,max.level=2)),file.path(out,'structure.txt'))
if(ok)write.csv(fit$pathways,file.path(out,'coefficients.csv'),row.names=FALSE)
write_json(list(stage=stage,status=if(ok)'completed'else'failed',operation_seconds=seconds,
 input_loading_seconds=loading_seconds,feature_preparation_seconds=feature_preparation_seconds,serialization_seconds=unname(serialize),
 gene_association_computed=FALSE,model_count=1,held_out_folds=0,chains=2,genes=nrow(x),features=ncol(x),
 warnings=unique(warn),error=if(ok)NULL else conditionMessage(fit),output_bytes=file.info(file.path(out,'output.rds'))$size,
 scope=settings$fit_scope),file.path(out,'timing.json'),auto_unbox=TRUE,pretty=TRUE)
if(!ok)quit(status=1)
