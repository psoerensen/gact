# Reuse audited CAD/T2D inputs. No downloads, extra filters, or stored LD.
args <- commandArgs(TRUE)
out <- args[1L]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(normalizePath('build/task-packages/library'), .libPaths()))
library(gcorr)
stopifnot(packageVersion('gcorr')=='0.1.3',
  gcorr_native_info()$gcorr$runtime_version=='0.1.2')
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
writeLines(capture.output(c(sessionInfo(),gcorr_native_info())),file.path(out,'session.txt'))
phases <- data.frame(phase=character(),seconds=double())
timed <- function(label, expression) {
  cat(format(Sys.time()),label,'started\n');flush.console()
  began <- proc.time()[[3L]]
  value <- force(expression)
  phases <<- rbind(phases,data.frame(phase=label,seconds=proc.time()[[3L]]-began))
  write.csv(phases,file.path(out,'phases.csv'),row.names=FALSE)
  cat(format(Sys.time()),label,'completed\n');flush.console()
  value
}
fit_pair <- function(p,...) gcorr_fit_hess(p,shared_sample_size=0,
  phenotypic_correlation=0,...)
export_fit <- function(fit,label) {
  write.csv(fit$estimates,file.path(out,paste0(label,'-estimates.csv')),row.names=FALSE)
  write.csv(fit$diagnostics$conditioning,file.path(out,paste0(label,'-conditioning.csv')),row.names=FALSE)
  write.csv(fit$diagnostics$inference$traits,file.path(out,paste0(label,'-traits.csv')),row.names=FALSE)
  counts <- do.call(rbind,lapply(split(fit$estimates,
    paste(fit$estimates$component,fit$estimates$quantity,fit$estimates$trait1,
      fit$estimates$trait2,sep=':')),function(x) data.frame(
        component=x$component[1],quantity=x$quantity[1],trait1=x$trait1[1],trait2=x$trait2[1],
        rows=nrow(x),available_points=sum(is.finite(x$estimate)),
        available_se=sum(is.finite(x$standard_error)),
        negative_variance=sum(x$uncertainty_reason=='Estimated HESS sampling variance is negative'),
        indefinite_covariance=sum(x$uncertainty_reason=='Estimated HESS sampling covariance is indefinite'))))
  write.csv(counts,file.path(out,paste0(label,'inference-counts.csv')),row.names=FALSE)
  invisible(fit)
}
# Confirm that the new uncertainty policy leaves previously saved points intact.
old <- timed('read_saved_cap20',readRDS(
  'build/benchmarks/manuscript-20261007-070903/gcorr_genome-4/output.rds'))
strict <- timed('saved_cap20_strict',suppressWarnings(fit_pair(old$projection,uncertainty='strict')))
moment <- timed('saved_cap20_moment',suppressWarnings(fit_pair(old$projection)))
stopifnot(isTRUE(all.equal(strict$estimates$estimate,old$pair_fit$estimates$estimate,tolerance=1e-12)),
  isTRUE(all.equal(moment$estimates$estimate,strict$estimates$estimate,tolerance=0)))
export_fit(moment,'saved-cap20-moment')
rm(old,strict,moment);invisible(gc())
input <- timed('read_audited_inputs',readRDS(
  'build/benchmarks/manuscript-20261007-003431/gcorr_inputs-4/output.rds'))
stopifnot(nrow(input$z)==7221923L,identical(rownames(input$n),input$markers$rsids),
  identical(colnames(input$n),colnames(input$z)))
write.csv(data.frame(trait=colnames(input$n),markers=nrow(input$n),
  n_min=apply(input$n,2,min),n_max=apply(input$n,2,max),source='extracted stat$n'),
  file.path(out,'analysis-input-audit.csv'),row.names=FALSE)
regions <- input$regions
index <- match(unlist(regions,use.names=FALSE),input$markers$rsids)
stopifnot(!anyNA(index),!anyDuplicated(index),length(index)==7221923L)
statistics <- setNames(lapply(seq_len(ncol(input$z)),function(j) data.frame(
  rsids=input$markers$rsids[index],chr=input$markers$chr[index],pos=input$markers$pos[index],
  ea=input$markers$ea[index],nea=input$markers$nea[index],
  z=input$z[index,j],n=input$n[index,j],stringsAsFactors=FALSE)),colnames(input$z))
projection <- timed('automatic_preparation',gcorr_prepare_hess(statistics,input$Glist,
  regions,nthreads=4))
stopifnot(projection$rank_policy$maximum_components==11)
rm(input,statistics,index,regions);invisible(gc())
write.csv(do.call(rbind,Map(function(x,name) data.frame(region=name,
  markers=length(x$marker_ids),retained_components=x$retained_components,
  matrix_space=x$matrix_space,t(x$timings_seconds),check.names=FALSE),
  projection$regions,names(projection$regions))),file.path(out,'regional-phases.csv'),row.names=FALSE)
raw <- timed('automatic_moment_fit',fit_pair(projection))
export_fit(raw,'automatic-moment')
partial <- timed('automatic_partial_meta_gc_fit',fit_pair(projection,
  genomic_control=c(GWAS1=1,GWAS2=1.18),genomic_control_provenance=c(
    GWAS1='Unchanged supplied Z scale',
    GWAS2='Nikpay 2015 Figure 4: final additive meta-analysis lambda 1.18; partial undo only, cohort-specific GC remains')))
export_fit(partial,'automatic-partial-meta-gc')
timed('save_compact_projection',saveRDS(projection,file.path(out,'projection.rds'),compress=FALSE))
writeLines('Completed: original inputs and cap20 results preserved.',file.path(out,'completed.txt'))
