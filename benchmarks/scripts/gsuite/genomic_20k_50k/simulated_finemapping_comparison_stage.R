# The published simulated-human example, using its frozen qgg inputs.
# This compares corresponding BayesC workflows, not identical samplers.
args <- commandArgs(TRUE); stage <- args[1]; threads <- as.integer(args[2]); out <- args[3]
root <- dirname(out); stopifnot(threads == 1L)
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
.libPaths(c('build/task-packages/library','build/task-packages/ecosystem/library','build/r-library',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
if(requireNamespace('RhpcBLASctl',quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(threads)
packages <- c('qgg','gbase','gmap','gbayes')
for(p in packages) stopifnot(requireNamespace(p,quietly=TRUE))
set.seed(20261005)
writeLines(c(capture.output(sessionInfo()),sapply(packages,function(p)
  paste(p,packageVersion(p),find.package(p)))),file.path(out,'session.txt'))
frozen <- 'build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds'
input <- readRDS(frozen)
load_stage <- function(name) {
  path <- file.path(root,paste0(name,'-1'),'output.rds')
  for(run in c(Sys.getenv('SIM_COMPARE_REUSE'),Sys.getenv('QGG_BENCH_REUSE_RUN'))) {
    if(!file.exists(path) && nzchar(run)) {
      candidate <- file.path(run,paste0(name,'-1'),'output.rds')
      if(file.exists(candidate)) path <- candidate
    }
  }
  readRDS(path)
}
sets <- lapply(input$sets,function(s) input$stat$rsids[input$stat$rsids %in% s])
stopifnot(identical(as.integer(lengths(sets)),c(850L,785L,1091L,1399L)),length(input$rsids)==37991L)
write.csv(data.frame(region=names(sets),markers=lengths(sets),causal=input$causal),file.path(out,'regions.csv'),row.names=FALSE)
write.csv(data.frame(path=frozen,bytes=file.info(frozen)$size,md5=unname(tools::md5sum(frozen))),file.path(out,'input.csv'),row.names=FALSE)
if(stage=='cmp_qgg_credible') genome <- load_stage('cmp_qgg_genome')
if(stage=='cmp_suite_ld') G <- load_stage('cmp_suite_prepare')
if(stage=='cmp_suite_genome') G <- load_stage('cmp_suite_ld')
if(stage=='cmp_suite_credible_oracle') genome <- load_stage('cmp_suite_genome')
if(stage=='cmp_suite_regional') sourceG <- load_stage('cmp_suite_prepare')
familiar <- identical(Sys.getenv('SIM_COMPARE_QGG_INTERFACES'),'1')
swap_probability <- as.numeric(Sys.getenv('SIM_COMPARE_LD_SWAP_PROBABILITY','0'))
stopifnot(length(swap_probability)==1L,is.finite(swap_probability),
  swap_probability>=0,swap_probability<=1)
swap_settings <- list(probability=swap_probability,minimum_r2=.8,
  maximum_friends=50L,moves=1L)
write.csv(data.frame(familiar_interface=familiar,probability=swap_probability,
  minimum_r2=.8,maximum_friends=50L,moves=1L),
  file.path(out,'ld-swap-settings.csv'),row.names=FALSE)
warnings <- character(); started <- Sys.time(); clock <- proc.time()[3]
value <- tryCatch(withCallingHandlers(switch(stage,
 cmp_qgg_ld=qgg::gprep(Glist=input$Glist,task='sparseld',msize=1000,
   rsids=input$rsids,ldfiles=file.path(out,'human.ld'),overwrite=FALSE,ncores=threads),
 cmp_qgg_regional=qgg::gmap(Glist=input$Glist,stat=input$stat,sets=input$sets,
   method='bayesC',algorithm='mcmc',output='full',verbose=FALSE),
 cmp_qgg_genome=qgg::gbayes(Glist=input$Glist,stat=input$stat,
   method='bayesC',algorithm='mcmc',verbose=FALSE),
 cmp_qgg_credible=qgg::getCredibleSets(Glist=input$Glist,fit=genome,sets=input$sets),
 cmp_suite_prepare={
   bim <- data.table::fread(input$Glist$bimfiles,header=FALSE)
   invalid_map <- which(!is.finite(bim[[3]]) | bim[[3]] < 0)
   bim[[3]][invalid_map] <- 0
   local_bim <- file.path(out,'human-validated-map.bim')
   data.table::fwrite(bim,local_bim,sep='\t',col.names=FALSE,quote=FALSE)
   write.csv(data.frame(line=invalid_map,marker=bim[[2]][invalid_map],replacement_map=0),
     file.path(out,'genetic-map-adjustments.csv'),row.names=FALSE)
   G <- gbase::gprep(study='Frozen simulated-human website example',
     bedfiles=input$Glist$bedfiles,bimfiles=local_bim,famfiles=input$Glist$famfiles,
     rsids=input$rsids,nthreads=threads)
   G$assembly <- 'simulated-human-tutorial'
   stopifnot(identical(G$ids,input$Glist$ids),identical(unlist(G$rsids,use.names=FALSE),input$stat$rsids),
     identical(unlist(G$a1,use.names=FALSE),input$stat$ea),identical(unlist(G$a2,use.names=FALSE),input$stat$nea))
   G
 },
 cmp_suite_ld={
   G <- gbase::gprep(G,task='sparseld',out_prefix=file.path(out,'human'),r2=0,
     max_distance_bp=0,max_distance_variants=1000,block_size=1024,nthreads=threads)
   # Compare reopened CSR to qgg's frozen float32 band at sampled positions.
   m <- nrow(G$sparseLD$marker_metadata); selected <- unique(as.integer(c(1:8,1000:1007,19000:19007,(m-7):m)))
   csr <- gbase::ld_region(G,G$sparseLD$marker_metadata$marker_id[selected])
   band <- matrix(0,length(selected),length(selected)); con <- file(input$Glist$ldfiles,'rb')
   for(k in seq_along(selected)) {
     seek(con,where=as.double(selected[k]-1)*2001*4,origin='start')
     row <- readBin(con,'numeric',n=2001,size=4,endian='little'); stopifnot(length(row)==2001)
     offset <- selected-selected[k]+1001; valid <- offset>=1 & offset<=2001
     band[k,valid] <- row[offset[valid]]
   }
   close(con); diag(band) <- 1
   delta <- max(abs(csr-band)); stopifnot(is.finite(delta),delta<2e-5)
   attr(G,'comparison_ld_audit') <- list(sampled_markers=length(selected),maximum_absolute_difference=delta)
   G
 },
 cmp_suite_regional={
   if(familiar) {
     fitted <- gmap::gmap(Glist=sourceG,stat=input$stat,sets=sets,method='bayesC',
       algorithm='mcmc',nit=1000L,nburn=100L,seed=10L,output='full',
       control=list(retain_diagnostics=TRUE,ld_swap=swap_settings),nthreads=threads)
     fits <- lapply(names(sets),function(region) {
       single <- fitted
       single$estimates <- fitted$estimates[fitted$estimates$region==region,,drop=FALSE]
       single$stat <- fitted$stat[fitted$stat$region==region,,drop=FALSE]
       single$sets <- fitted$sets[region];single$details <- fitted$details[region]
       single$diagnostic_summary <- fitted$diagnostic_summary[fitted$diagnostic_summary$region==region,,drop=FALSE]
       single$selection$regions <- fitted$selection$regions[region]
       single
     });names(fits) <- names(sets)
     timing <- fitted$timings
     timing$markers <- lengths(sets)[match(timing$region,names(sets))]
     timing$ld_seconds <- timing$ld_total_seconds
     timing$model_and_sets_seconds <- timing$prepare_seconds+timing$binding_seconds+
       timing$sampling_seconds+timing$sets_seconds
     write.csv(timing,file.path(out,'regional-times.csv'),row.names=FALSE)
     write.csv(fitted$ld_swap,file.path(out,'ld-swap-counts.csv'),row.names=FALSE)
     list(fits=fits,timings=timing,full_fit=fitted,
       interpretation='qgg-style interface; full BED LD in memory; SPD-checked dense sampler; native qgg-compatible CS2 groups; one chain, no convergence claim')
   } else {
   s <- input$stat; n <- median(s$n)
   yy <- median((s$b^2+(s$n-2)*s$seb^2)*s$n); vy <- yy/(n-1)
   vb <- (vy*.5)/(5000*.01); fits <- timing <- setNames(vector('list',length(sets)),names(sets))
   for(i in seq_along(sets)) {
     begin <- proc.time()[3]
     gene <- gbase::gprep(bedfiles=sourceG$bedfiles,bimfiles=sourceG$bimfiles,
       famfiles=sourceG$famfiles,rsids=sets[[i]],nthreads=threads)
     gene$assembly <- 'simulated-human-tutorial'
     prefix <- file.path(out,paste0('region-',i))
     gene <- gbase::gprep(gene,task='sparseld',out_prefix=prefix,r2=0,
       max_distance_bp=0,max_distance_variants=0,allow_full_ld=TRUE,nthreads=threads)
     ld_seconds <- unname(proc.time()[3]-begin)
     r <- match(gene$rsids[[1]],s$rsids)
     stopifnot(!anyNA(r),identical(s$ea[r],gene$a1[[1]]),identical(s$nea[r],gene$a2[[1]]))
     st <- data.frame(rsids=s$rsids[r],ea=s$ea[r],nea=s$nea[r],score=s$wy[r],information=s$ww[r])
     analysis <- proc.time()[3]
     fits[[i]] <- gmap::gmap(st,gene,setNames(list(gene$rsids[[1]]),names(sets)[i]),method='bayesc',
       prior=list(residual_variance=vy*.5,effect_variance=vb,inclusion_probability=.01,
         effect_variance_prior=list(df=4,scale=.5*vb),residual_variance_prior=list(df=4,scale=.25*vy),weight_prior=c(1,1)),
       control=list(burnin=100L,sampling_sweeps=1000L,seeds=10L,retain_diagnostics=TRUE,
         estimate_effect_variance=TRUE,estimate_residual_variance=TRUE,estimate_weights=TRUE,
         eigen=list(retained_positive_mass=1,select_on_convergence=FALSE,verification_tolerance=1e-7),
         sets=list(coverage=.9,min_r2=.5)),
       normalization=list(trait=list(response_ss=yy,likelihood_dimension=n,in_sample_ld=TRUE)))
     timing[[i]] <- data.frame(region=names(sets)[i],markers=length(r),ld_seconds=ld_seconds,
       model_and_sets_seconds=unname(proc.time()[3]-analysis),ld_bytes=sum(file.info(unlist(gene$sparseLD$paths))$size))
     # Only exact owned files are removed; source BED and all outputs are retained.
     unlink(unlist(gene$sparseLD$paths),recursive=FALSE)
     cat('Completed region',names(sets)[i],'\n'); flush.console()
   }
   write.csv(do.call(rbind,timing),file.path(out,'regional-times.csv'),row.names=FALSE)
   list(fits=fits,timings=do.call(rbind,timing),working_normalization=list(yy=yy,n=n,vy=vy),
     interpretation='qgg dense working wy/ww; full regional LD; native eigen representation with no mass truncation; one chain, no convergence claim')
   }
 },
 cmp_suite_genome={
   if(familiar) {
     fit <- gbayes::gbayes(stat=input$stat,Glist=G,method='bayesC',algorithm='mcmc',
       nit=100L,nburn=0L,seed=20261005L,
       control=list(threads=threads,max_ld_entries=1e8,retain_diagnostics=TRUE,
         ld_swap=swap_settings))
     write.csv(as.data.frame(fit$ld_swap),file.path(out,'ld-swap-counts.csv'),row.names=FALSE)
     attr(fit,'comparison_units') <- list(interpretation='qgg-style interface over standardized summary likelihood; reference-LD residual adjustment .9; not qgg variable ww/wy and neighbourhood adjustment')
     fit
   } else {
   s <- input$stat; m <- nrow(s); n <- median(s$n)
   yy <- median((s$b^2+(n-2)*s$seb^2)*s$ww); vy <- yy/(n-1)
   vb <- .5/(m*.001)
   fit <- gbayes::gbayes(stat=s,Glist=G,method='bayesc',trait='simulation',
     prior=list(residual_variance=.5,effect_variance=vb,inclusion_probability=.001,
       effect_variance_prior=list(df=4,scale=.5*vb),residual_variance_prior=list(df=4,scale=.25),weight_prior=c(1,1)),
     control=list(burnin=0L,sampling_sweeps=100L,seeds=20261005L,threads=threads,
       estimate_effect_variance=TRUE,estimate_residual_variance=TRUE,estimate_weights=TRUE,
       residual_policy='reference_ld',residual_adjustment=.9,max_ld_entries=1e8,
       retain_diagnostics=TRUE))
   attr(fit,'comparison_units') <- list(qgg_phenotype_variance=vy,
     interpretation='standalone standardized summary likelihood and reference-LD residual policy; qgg retains supplied variable ww/wy and its LD-neighborhood variance adjustment')
   fit
   }
 },
 cmp_suite_credible_oracle={
   if(familiar) {
     sourceG <- load_stage('cmp_suite_prepare')
     result <- gmap::getCredibleSets(Glist=sourceG,fit=genome,sets=input$sets,
       nthreads=threads)
     result
   } else {
   # A common qgg postprocessor, not a standalone gsuite capability or runtime.
   adapter <- list(stat=data.frame(rsids=genome$estimates$marker,dm=genome$estimates$pip))
   result <- qgg::getCredibleSets(Glist=input$Glist,fit=adapter,sets=input$sets)
   attr(result,'implementation') <- 'qgg common postprocessor of standalone gbayes PIPs; standalone summary-posterior LD grouping API unavailable'
   result
   }
 },
 cmp_compare={
   q <- load_stage('cmp_qgg_regional'); b <- load_stage('cmp_suite_regional')
   rows <- lapply(seq_along(sets),function(i) {
     old <- q$stat; old <- old[old$rsids %in% sets[[i]],]
     new <- b$fits[[i]]$estimates; at <- match(new$marker,old$rsids); stopifnot(!anyNA(at))
     pip <- new$pip; truth <- input$causal[i]
     data.frame(region=names(sets)[i],markers=length(pip),qgg_causal_pip=old$dm[match(truth,old$rsids)],
       suite_causal_pip=pip[match(truth,new$marker)],pip_correlation=cor(pip,old$dm[at]),
       mean_absolute_pip_difference=mean(abs(pip-old$dm[at])),
       qgg_top=old$rsids[which.max(old$dm)],suite_top=new$marker[which.max(pip)],
       suite_status=b$fits[[i]]$diagnostic_summary$status[1])
   })
   result <- do.call(rbind,rows); write.csv(result,file.path(out,'regional-comparison.csv'),row.names=FALSE)
   qg <- load_stage('cmp_qgg_genome'); bg <- load_stage('cmp_suite_genome')
   at <- match(bg$estimates$marker,qg$stat$rsids); stopifnot(!anyNA(at))
   genome <- data.frame(marker=input$causal,qgg_pip=qg$stat$dm[match(input$causal,qg$stat$rsids)],
     suite_pip=bg$estimates$pip[match(input$causal,bg$estimates$marker)])
   write.csv(genome,file.path(out,'genome-causal-comparison.csv'),row.names=FALSE)
   list(regions=result,genome_causal=genome,genome_pip_correlation=cor(bg$estimates$pip,qg$stat$dm[at]),
     genome_mean_absolute_pip_difference=mean(abs(bg$estimates$pip-qg$stat$dm[at])),
     qgg_sets=load_stage('cmp_qgg_credible'),suite_pips_qgg_sets=load_stage('cmp_suite_credible_oracle'))
 },stop('Unknown stage')),warning=function(w) {warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
seconds <- unname(proc.time()[3]-clock); failed <- inherits(value,'error')
if(failed) warnings <- c(warnings,conditionMessage(value))
serialization <- unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=seconds,
  serialization_seconds=serialization,status=if(failed)'failed'else'completed',warnings=paste(unique(warnings),collapse='; ')),
  file.path(out,'timing.csv'),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,'output-structure.txt'))
if(failed) quit(status=1)
