# Refit saved compact moments only; no genotype/LD/eigen preparation.
# Usage: Rscript --vanilla .../verify_hess_software_uncertainty.R OUTPUT_DIR
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
Sys.setenv(LC_ALL='C',LANG='C')
.libPaths(c(normalizePath('build/task-packages/library'),.libPaths()))
library(gcorr)
stopifnot(packageVersion('gcorr')=='0.1.4',
  gcorr_native_info()$gcorr$runtime_version=='0.1.3')
out <- args[1L]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
p <- readRDS('build/gact-revision-cleanup/hess-moment/projection.rds')
results <- list()
timings <- data.frame()
for(scenario in c('raw','partial_gc')) {
  start <- proc.time()[[3L]]
  fit <- if(scenario=='raw') gcorr_fit_hess(p,shared_sample_size=0,
    phenotypic_correlation=0) else gcorr_fit_hess(p,shared_sample_size=0,
    phenotypic_correlation=0,genomic_control=c(GWAS1=1,GWAS2=1.18),
    genomic_control_provenance=c(GWAS1='Unchanged supplied Z scale',
      GWAS2='Nikpay 2015 Figure 4: partial final-meta-layer 1.18 sensitivity only'))
  timings <- rbind(timings,data.frame(scenario=scenario,
    fit_seconds=proc.time()[[3L]]-start))
  previous <- read.csv(file.path('build/gact-revision-cleanup/hess-moment',
    if(scenario=='raw') 'automatic-moment-estimates.csv' else
      'automatic-partial-meta-gc-estimates.csv'))
  stopifnot(max(abs(fit$estimates$estimate-previous$estimate),na.rm=TRUE)<1e-12)
  write.csv(fit$estimates,file.path(out,paste0(scenario,'-estimates.csv')),row.names=FALSE)
  jack <- fit$diagnostics$genome_rg_jackknife
  write.csv(data.frame(region=names(p$regions),
    delete_region_rg=jack$delete_region_correlations),
    file.path(out,paste0(scenario,'-jackknife.csv')),row.names=FALSE)
  # Independent direct deletion and full-sample centering from exported rows.
  region <- fit$estimates[fit$estimates$component=='region',]
  total <- fit$estimates[fit$estimates$component=='total',]
  h1 <- region$estimate[region$quantity=='h2' & region$trait1=='GWAS1']
  h2 <- region$estimate[region$quantity=='h2' & region$trait1=='GWAS2']
  cov <- region$estimate[region$quantity=='cov']
  full <- sum(cov)/sqrt(sum(h1)*sum(h2))
  delete <- (sum(cov)-cov)/sqrt((sum(h1)-h1)*(sum(h2)-h2))
  se <- sqrt((length(delete)-1)*mean((delete-full)^2))
  stopifnot(max(abs(delete-jack$delete_region_correlations))<1e-11,
    abs(se-total$standard_error[total$quantity=='rg'])<1e-11)
  # Execute the unchanged pinned upstream variance function at native points
  # separately in the Python comparison driver to isolate its rank convention.
  results[[scenario]] <- cbind(scenario=scenario,total)
}
write.csv(do.call(rbind,results),file.path(out,'native-totals.csv'),row.names=FALSE)
write.csv(timings,file.path(out,'fit-timings.csv'),row.names=FALSE)
writeLines(capture.output(c(sessionInfo(),gcorr_native_info())),file.path(out,'session.txt'))
print(do.call(rbind,results));print(timings)
