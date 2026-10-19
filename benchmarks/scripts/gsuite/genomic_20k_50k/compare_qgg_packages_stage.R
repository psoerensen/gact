# Corresponding standalone workflows. Gene LD is computed on the fly; no
# persistent gene-LD cache. See qgg-package-comparison-protocol.md for limits.
args<-commandArgs(TRUE);stage<-args[1];out<-args[3];threads<-as.integer(args[2]);stopifnot(threads>=1L)
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
.libPaths(c('build/task-packages/library','build/task-packages/ecosystem/library','build/r-library',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(threads)
packages<-c('gbase','glma','gsea','gmap')
for(p in packages)stopifnot(requireNamespace(p,quietly=TRUE))
writeLines(c(capture.output(sessionInfo()),sapply(packages,function(p)paste(p,packageVersion(p),find.package(p)))),file.path(out,'session.txt'))
baseline<-Sys.getenv('QGG_COMPARE_RUN');stopifnot(nzchar(baseline))
f<-readRDS(file.path(Sys.getenv('QGG_BENCH_REUSE_RUN'),'qgg_inputs-1','output.rds'))
G<-f$Glist;set.seed(20261004)
prefix<-file.path(out,'temporary_gene')
owned<-paste0(prefix,c('.bed','.bim','.fam','.ld.row_ptr.u64.bin','.ld.col_idx.u32.0based.bin','.ld.values.f32.bin','.ld.markers.tsv','.ld.meta.txt'))
stopifnot(!any(file.exists(owned)))
record_bytes<-ceiling(G$n/4)
bed_source<-file(G$bedfiles,'rb')

# Copy the exact packed records, without decoding or simulating genotypes.
# Native public preparation then validates the temporary PLINK inputs.
prepare_gene<-function(k){
 stopifnot(length(k)>2,all(diff(k)>0),length(unique(G$chr[[1]][k]))==1)
 first<-min(k);last<-max(k)
 seek(bed_source,where=3+(first-1)*record_bytes,origin='start')
 packed<-readBin(bed_source,'raw',n=record_bytes*(last-first+1))
 stopifnot(length(packed)==record_bytes*(last-first+1))
 dim(packed)<-c(record_bytes,last-first+1)
 writeBin(c(as.raw(c(108,27,1)),as.vector(packed[,k-first+1,drop=FALSE])),paste0(prefix,'.bed'))
 bim<-data.frame(chr=G$chr[[1]][k],rsid=G$rsids[[1]][k],map=G$map[[1]][k],pos=G$pos[[1]][k],a1=G$a1[[1]][k],a2=G$a2[[1]][k])
 write.table(bim,paste0(prefix,'.bim'),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
 stopifnot(file.copy(G$famfiles,paste0(prefix,'.fam'),overwrite=FALSE))
 geneG<-gbase::gprep(bedfiles=paste0(prefix,'.bed'),nthreads=threads)
 stopifnot(identical(geneG$ids,G$ids),identical(geneG$rsids[[1]],G$rsids[[1]][k]),max(abs(geneG$af[[1]]-G$af[[1]][k]))<1e-12)
 geneG$assembly<-'GRCh37'
 gbase::gprep(Glist=geneG,task='sparseld',out_prefix=paste0(prefix,'.ld'),r2=0,max_distance_bp=0,max_distance_variants=0,allow_full_ld=TRUE,nthreads=threads)
}
cleanup_gene<-function(){
 # Fixed owned filenames only; never remove a directory or reference input.
 unlink(owned,recursive=FALSE);stopifnot(!any(file.exists(owned)))
}

# Build identifier mappings once, as part of the operation clock below.
warn<-character();started<-Sys.time();clock<-proc.time()[3]
value<-tryCatch(withCallingHandlers(switch(stage,
 suite_joint={
  z<-f$gene_stat;traits<-colnames(z);genes<-rownames(z)
  stat<-setNames(lapply(seq_len(ncol(z)),function(j)data.frame(gene=genes,score=z[,j])),traits)
  sets<-lapply(f$kegg,function(s)intersect(s,genes))
  patterns<-rbind(null=c(0,0),first=c(1,0),second=c(0,1),both=c(1,1));colnames(patterns)<-traits
  I<-diag(2);dimnames(I)<-list(traits,traits)
  effect<-.5/(length(sets)*.01)*I
  gsea::gsea(stat,sets,'bayesc',sampling=list(gene_covariance='independent',trait_covariance=.5*I),
   prior=list(effect_covariance=effect,patterns=patterns,pattern_weights=c(null=.99,first=.01/3,second=.01/3,both=.01/3),covariance_prior=list(df=4,scale=2*effect),pattern_prior=rep(1,4)),
   control=list(burnin=200L,sampling_sweeps=1000L,seeds=20261004L,threads=threads,estimate_covariance=TRUE,estimate_pattern_weights=TRUE))
 },
 suite_vegas={
  stopifnot(!any(vapply(f$marker_sets,anyDuplicated,integer(1))))
  shared<-Reduce(intersect,lapply(f$vegas_stat,function(s)s$rsids))
  sets<-lapply(f$marker_sets,function(s)intersect(s,shared));sets<-sets[lengths(sets)>2]
  sizes<-lengths(sets)
  old<-readRDS(file.path(baseline,'qgg_vegas-4','output.rds'))
  stopifnot(identical(names(sets),rownames(old$p)),identical(as.integer(sizes),as.integer(old$genes$m)))
  flat<-unlist(sets,use.names=FALSE)
  mapped<-match(flat,G$rsids[[1]]);stopifnot(!anyNA(mapped))
  index<-split(mapped,rep(seq_along(sets),sizes))
  multichr<-vapply(index,function(k)length(unique(G$chr[[1]][k])),integer(1))
  if(any(multichr!=1L))stop('Cannot match full within-gene LD across chromosomes: ',paste(names(sets)[multichr!=1L],collapse=', '))
  stats_index<-lapply(f$vegas_stat,function(s)split(match(flat,s$rsids),rep(seq_along(sets),sizes)))
  traits<-names(f$vegas_stat);result<-setNames(lapply(traits,function(x)vector('list',length(sets))),traits)
  timing<-vector('list',length(sets));diagnostics<-vector('list',length(sets))
  for(i in seq_along(sets)){
   begin<-proc.time()[3];ord<-order(index[[i]]);k<-index[[i]][ord]
   geneG<-prepare_gene(k);ld_seconds<-unname(proc.time()[3]-begin)
   st<-setNames(lapply(seq_along(traits),function(j){
    s<-f$vegas_stat[[j]][stats_index[[j]][[i]][ord],]
    stopifnot(identical(s$rsids,geneG$rsids[[1]]),identical(s$ea,geneG$a1[[1]]),identical(s$nea,geneG$a2[[1]]))
    data.frame(rsids=s$rsids,ea=s$ea,nea=s$nea,z=sqrt(qchisq(pmax(s$p,1e-10),1,lower.tail=FALSE)))
   }),traits)
   analysis_begin<-proc.time()[3]
   fit<-glma::vegas(st,geneG,setNames(list(geneG$rsids[[1]]),names(sets)[i]),control=list(maximum_markers=max(2048L,length(k)),relative_psd_tolerance=1e-7,tail_method='saddlepoint'),threads=threads)
   analysis_seconds<-unname(proc.time()[3]-analysis_begin)
   for(j in traits)result[[j]][[i]]<-fit$stat[[j]]
   diagnostics[[i]]<-fit$diagnostics
   timing[[i]]<-data.frame(gene=names(sets)[i],markers=length(k),ld_prepare_seconds=ld_seconds,analysis_seconds=analysis_seconds,resource_bytes=sum(file.info(owned)$size,na.rm=TRUE))
   cleanup_gene()
   if(i%%100==0||i==length(sets)){write.csv(do.call(rbind,timing[seq_len(i)]),file.path(out,'gene-progress.csv'),row.names=FALSE);cat('Completed',i,'of',length(sets),'genes\n');flush.console()}
  }
  list(stat=lapply(result,function(x)do.call(rbind,x)),gene_timings=do.call(rbind,timing),diagnostics=diagnostics,temporary_ld=TRUE)
 },
 suite_mapping={
  s<-f$fine_stat;yy<-median((s$b^2+(s$n-2)*s$seb^2)*s$n);n<-median(s$n);vy<-yy/(n-1)
  ww<-(yy/n)/(s$seb^2+s$b^2/s$n);wy<-s$b*ww
  sets<-lapply(f$fine_sets,function(x)intersect(x,s$rsids));sizes<-lengths(sets)
  flat<-unlist(sets,use.names=FALSE);gindices<-split(match(flat,G$rsids[[1]]),rep(seq_along(sets),sizes));sindices<-split(match(flat,s$rsids),rep(seq_along(sets),sizes))
  stopifnot(!anyNA(unlist(gindices)),!anyNA(unlist(sindices)))
  fits<-timing<-setNames(vector('list',length(sets)),names(sets));weights<-c(.99,.006,.003,.001);multipliers<-c(0,.01,.1,1)
  vb<-(vy*.5)/(5000*sum(weights*multipliers))
  for(i in seq_along(sets)){
   begin<-proc.time()[3];ord<-order(gindices[[i]]);k<-gindices[[i]][ord];r<-sindices[[i]][ord]
   geneG<-prepare_gene(k);ld_seconds<-unname(proc.time()[3]-begin)
   stopifnot(identical(s$rsids[r],geneG$rsids[[1]]),identical(s$ea[r],geneG$a1[[1]]),identical(s$nea[r],geneG$a2[[1]]))
   st<-data.frame(rsids=s$rsids[r],ea=s$ea[r],nea=s$nea[r],score=wy[r],information=ww[r])
   analysis_begin<-proc.time()[3]
   fits[[i]]<-gmap::gmap(st,geneG,setNames(list(geneG$rsids[[1]]),names(sets)[i]),method='bayesr',
    prior=list(residual_variance=vy*.5,effect_variance=vb,weights=weights,variance_multipliers=multipliers,effect_variance_prior=list(df=4,scale=2*vb),weight_prior=rep(1,4)),
    control=list(burnin=100L,sampling_sweeps=1000L,seeds=10L,estimate_effect_variance=TRUE,estimate_weights=TRUE,retain_diagnostics=TRUE,eigen=list(maximum_dimension=max(2048L,length(k)),verification_tolerance=1e-7,retained_positive_mass=c(.995,.99,.95,.9),select_on_convergence=FALSE)))
   timing[[i]]<-data.frame(gene=names(sets)[i],markers=length(k),ld_prepare_seconds=ld_seconds,analysis_seconds=unname(proc.time()[3]-analysis_begin),resource_bytes=sum(file.info(owned)$size,na.rm=TRUE))
   cleanup_gene();cat('Completed mapping',names(sets)[i],'\n');flush.console()
  }
  list(fits=fits,gene_timings=do.call(rbind,timing),working_summary=list(yy=yy,n=n,vy=vy,formula='qgg dense working ww/wy; not qgg eigen transform and not OLS conversion'),residual_policy='fixed: external LD does not support native residual learning',temporary_ld=TRUE)
 },stop('Unknown comparison stage')),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds<-unname(proc.time()[3]-clock);close(bed_source)
failed<-inherits(value,'error');if(failed)warn<-c(warn,conditionMessage(value))
serialization<-unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads=threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=seconds,serialization_seconds=serialization,status=if(failed)'failed'else'completed',warnings=paste(unique(warn),collapse='; ')),file.path(out,'timing.csv'),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,'output-structure.txt'));if(failed)quit(status=1)
