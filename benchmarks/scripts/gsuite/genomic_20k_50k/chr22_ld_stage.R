args <- commandArgs(TRUE)
stage <- args[1]; threads <- as.integer(args[2]); out <- args[3]; run <- dirname(out)
writeLines(as.character(Sys.getpid()), file.path(out, 'r-process.pid'))
.libPaths(c('build/task-packages/library','build/task-packages/ecosystem/library','build/r-library',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
if(requireNamespace('RhpcBLASctl',quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(threads)
library(qgg); library(gbase)
writeLines(c(capture.output(sessionInfo()),paste('qgg',packageVersion('qgg'),find.package('qgg')),
 paste('gbase',packageVersion('gbase'),find.package('gbase')),capture.output(gbase::gbase_native_info())),file.path(out,'session.txt'))
load_stage <- function(s) readRDS(file.path(run,paste0(s,'-1'),'output.rds'))
if(stage %in% c('prepare_qgg','prepare_gbase')) prefix <- load_stage('reference_fixture')$prefix
if(stage=='ld_qgg') G <- load_stage('prepare_qgg')
if(stage=='ld_gbase') G <- load_stage('prepare_gbase')
if(stage=='validate_ld') {Q<-load_stage('ld_qgg'); B<-load_stage('ld_gbase')}
warnings <- character(); started <- Sys.time(); clock <- proc.time()[3]
value <- tryCatch(withCallingHandlers(switch(stage,
 reference_fixture={
  original <- 'build/benchmarks/manuscript-20261004-171426/qgg_glist-1/g1000_eur_filtered'
  bim <- data.table::fread(paste0(original,'.bim'),header=FALSE,select=c(1,2))
  markers <- bim[[2]][bim[[1]]==22]
  stopifnot(length(markers)>1000, !anyDuplicated(markers))
  prefix <- file.path(out,'g1000_eur_chr22')
  qgg::writeBED(bedRead=paste0(original,'.bed'),bimRead=paste0(original,'.bim'),famRead=paste0(original,'.fam'),
    bedWrite=paste0(prefix,'.bed'),bimWrite=paste0(prefix,'.bim'),famWrite=paste0(prefix,'.fam'),rsids=markers)
  stopifnot(nrow(data.table::fread(paste0(prefix,'.bim'),header=FALSE))==length(markers))
  list(prefix=normalizePath(prefix,winslash='/'),markers=length(markers),samples=503L,chromosome=22L,
       window=1000L,r2_threshold=0,description='Shared physical chromosome 22 fixture; construction excluded from LD timings')
 },
 prepare_qgg=qgg::gprep(study='1000G EUR chr22',bedfiles=paste0(prefix,'.bed'),bimfiles=paste0(prefix,'.bim'),famfiles=paste0(prefix,'.fam'),ncores=threads),
 prepare_gbase=gbase::gprep(study='1000G EUR chr22',bedfiles=paste0(prefix,'.bed'),bimfiles=paste0(prefix,'.bim'),famfiles=paste0(prefix,'.fam'),nthreads=threads),
 ld_qgg=qgg::gprep(Glist=G,task='sparseld',msize=1000,rsids=unlist(G$rsids),ldfiles=file.path(out,'chr22.ld'),overwrite=FALSE,ncores=threads),
 ld_gbase=gbase::gprep(Glist=G,task='sparseld',out_prefix=file.path(out,'chr22'),max_distance_bp=0,max_distance_variants=1000,r2=0,block_size=1000,nthreads=threads,overwrite=FALSE),
 validate_ld={
  stopifnot(identical(Q$idsLD,B$idsLD),identical(unlist(Q$rsidsLD,use.names=FALSE),unlist(B$rsidsLD,use.names=FALSE)),Q$n==503L,B$n==503L)
  markers <- unlist(B$rsidsLD,use.names=FALSE); m<-length(markers); w<-1000L
  stopifnot(file.info(Q$ldfiles)$size==m*(2*w+1)*4)
  # Read sampled rows from the saved qgg band and compare reopened gbase CSR.
  starts<-unique(pmax(1L,pmin(m-31L,as.integer(c(1,1000,m/4,m/2,3*m/4,m-31)))))
  selected<-sort(unique(unlist(lapply(starts,function(i)i+0:31))))
  R<-gbase::ld_region(B,markers[selected])
  band<-matrix(0,length(selected),length(selected))
  con<-file(Q$ldfiles,'rb')
  for(k in seq_along(selected)) {
   i<-selected[k];seek(con,where=as.double(i-1)*(2*w+1)*4,origin='start')
   row<-readBin(con,'numeric',n=2*w+1,size=4,endian='little');stopifnot(length(row)==2*w+1)
   offset<-selected-i+w+1;valid<-offset>=1 & offset<=2*w+1
   band[k,valid]<-row[offset[valid]]
  }
  close(con)
  # gbase uses an implicit unit diagonal; qgg readers also restore that diagonal.
  stored_diagonal<-diag(band);diag(band)<-1
  max_delta<-max(abs(R-band));stopifnot(is.finite(max_delta),max_delta<2e-5)
  q_scores<-Q$ldscores[[1]]; b_scores<-B$ldscores[[1]]
  common<-match(names(q_scores),markers);stopifnot(!anyNA(common))
  score_delta<-max(abs(unname(q_scores)-unname(b_scores[common])))
  # Keep monomorphic-marker score differences visible rather than adjusting them away.
  list(markers=m,samples=B$n,window=w,r2=0,sampled_markers=length(selected),sampled_pairs=length(selected)^2,
       max_absolute_correlation_difference=max_delta,maximum_common_ld_score_difference=score_delta,
       qgg_ld_score_markers=length(q_scores),gbase_ld_score_markers=length(b_scores),
       missing_qgg_ld_score_markers=markers[!markers %in% names(q_scores)],
       sampled_qgg_nonunit_diagonal=sum(abs(stored_diagonal-1)>2e-5),
       qgg_ld_bytes=sum(file.info(Q$ldfiles)$size),gbase_ld_bytes=sum(file.info(unlist(B$sparseLD$paths))$size),
       gbase_nnz=B$sparseLD$nnz,
       note='Timing includes public sparse-LD construction, persistence and ordinary LD scores; input preparation and final RDS serialization are separate. Different fixed-band and CSR formats. Monomorphic qgg scores may be omitted; gbase diagonal is always one.')
 },stop('Unknown stage')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds <- unname(proc.time()[3]-clock);failed<-inherits(value,'error')
if(failed) warnings<-c(warnings,conditionMessage(value))
ld_bytes<-if(!failed&&stage=='ld_qgg') sum(file.info(value$ldfiles)$size) else if(!failed&&stage=='ld_gbase') sum(file.info(unlist(value$sparseLD$paths))$size) else NA_real_
serialization<-unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=seconds,serialization_seconds=serialization,ld_bytes,
 status=if(failed)'failed'else'completed',warnings=paste(unique(warnings),collapse='; ')),file.path(out,'timing.csv'),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,'output-structure.txt'))
if(failed)quit(status=1)
