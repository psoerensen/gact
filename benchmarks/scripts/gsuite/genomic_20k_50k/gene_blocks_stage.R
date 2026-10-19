args<-commandArgs(TRUE);stage<-args[1];threads<-as.integer(args[2]);out<-args[3]
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
Sys.setenv(OPENBLAS_NUM_THREADS=threads,OMP_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
library(qgg)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(threads)
writeLines(c(capture.output(sessionInfo()),paste('qgg',packageVersion('qgg'),find.package('qgg'))),file.path(out,'session.txt'))
cat('START loading reference once\n');flush.console();loading<-proc.time()
f<-readRDS('build/benchmarks/manuscript-20261004-174351/qgg_inputs-1/output.rds')
old<-readRDS('build/benchmarks/manuscript-20261004-194031/qgg_vegas-4/output.rds')$genes
loading<-proc.time()-loading
cat('START one-time marker alignment\n');flush.console();alignment<-proc.time()
targets<-c(100,300,500,1000)
chosen<-lapply(targets,function(t)head(order(abs(old$m-t),old$EnsemblID),100))
ix<-unlist(chosen);genes<-old$EnsemblID[ix];expected<-old$m[ix];stopifnot(!anyDuplicated(genes))
sets<-f$marker_sets[genes];flat<-unlist(sets,use.names=FALSE)
shared<-Reduce(intersect,lapply(f$vegas_stat,function(s)s$rsids))
present<-!is.na(match(flat,shared));group<-rep(seq_along(sets),lengths(sets))
indices<-split(match(flat[present],f$Glist$rsids[[1]]),group[present])
stopifnot(!anyNA(unlist(indices)),identical(as.integer(lengths(indices)),as.integer(expected)))
G<-f$Glist;rm(f,flat,shared,present,group,sets);alignment<-proc.time()-alignment;invisible(gc())
write.csv(data.frame(step=c('input_loading','marker_alignment'),wall_seconds=c(loading[3],alignment[3]),cpu_seconds=c(sum(loading[1:2]),sum(alignment[1:2]))),file.path(out,'setup-timing.csv'),row.names=FALSE)
measure<-function(fun){a<-proc.time();v<-fun();d<-proc.time()-a;list(value=v,wall=unname(d[3]),cpu=unname(sum(d[1:2])))}
rows<-list();spectra<-list();error<-NULL;started<-Sys.time();clock<-proc.time()
tryCatch({
for(g in seq_along(targets)){
 target<-targets[g];set.seed(20261005+target);order_gene<-sample(((g-1)*100+1):(g*100))
 cat('START numerical blocks for',target,'markers\n');flush.console()
 for(j in seq_along(order_gene)){
  i<-order_gene[j];k<-indices[[i]];stopifnot(length(unique(G$chr[[1]][k]))==1L)
  extraction<-measure(function()qgg::getG(Glist=G,chr=1,cls=k,scale=TRUE));B<-extraction$value;extraction$value<-NULL
  ld<-measure(function()cor(B));R<-ld$value;ld$value<-NULL
  # Alternate decomposition order to reduce systematic warm-cache ordering.
  if(j%%2){full<-measure(function()eigen(R,symmetric=TRUE));values<-measure(function()eigen(R,symmetric=TRUE,only.values=TRUE))}
  else{values<-measure(function()eigen(R,symmetric=TRUE,only.values=TRUE));full<-measure(function()eigen(R,symmetric=TRUE))}
  stopifnot(all(is.finite(full$value$values)),max(abs(full$value$values-values$value$values))<1e-7,abs(sum(full$value$values)-length(k))<1e-6*length(k))
  spectra[[genes[i]]]<-full$value$values
  rows[[length(rows)+1]]<-data.frame(gene=genes[i],target_markers=target,markers=length(k),individuals=G$n,
   extraction_seconds=extraction$wall,ld_seconds=ld$wall,eigen_full_seconds=full$wall,eigen_values_seconds=values$wall,
   preparation_full_seconds=extraction$wall+ld$wall+full$wall,preparation_values_seconds=extraction$wall+ld$wall+values$wall,
   extraction_cpu_seconds=extraction$cpu,ld_cpu_seconds=ld$cpu,eigen_full_cpu_seconds=full$cpu,eigen_values_cpu_seconds=values$cpu,
   genotype_matrix_bytes=as.numeric(object.size(B)),ld_matrix_bytes=as.numeric(object.size(R)),eigensystem_bytes=as.numeric(object.size(full$value)),spectrum_bytes=as.numeric(object.size(full$value$values)))
  rm(B,R,full,values,extraction,ld)
  if(j%%10==0){write.csv(do.call(rbind,rows),file.path(out,'gene-progress.csv'),row.names=FALSE);cat('Completed',j,'of 100 genes; target',target,'\n');flush.console()}
 }
}
},error=function(e)error<<-e)
duration<-proc.time()-clock
value<-if(is.null(error))list(gene_timings=do.call(rbind,rows),eigenvalues=spectra,threads_requested=threads,
 timing_policy='Setup separate; component totals; group checkpoints; automatic GC only',
 marker_policy='100 complete real genes nearest each target after shared CAD/T2D QC')else error
serialization<-unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=unname(duration[3]),serialization_seconds=serialization,status=if(is.null(error))'completed'else'failed',warnings=if(is.null(error))''else conditionMessage(error)),file.path(out,'timing.csv'),row.names=FALSE)
if(!is.null(error))quit(status=1)
