# Summarize the focused serial pairs produced by time_ld_swap.ps1.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
runs <- read.csv(args[1],check.names=FALSE,dec=',')
destination <- args[2]
stopifnot(!dir.exists(destination),all(c('repeat','probability','run') %in% names(runs)))
dir.create(destination,recursive=TRUE)
frozen <- readRDS('build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds')
measures <- regions <- swaps <- causal <- diagnostics <- list()
reference <- list()
for(i in seq_len(nrow(runs))) {
  item <- runs[i,]
  process <- read.csv(file.path(item$run,'process-measurements.csv'),dec=',')
  for(stage in c('cmp_suite_regional','cmp_suite_genome')) {
    path <- file.path(item$run,paste0(stage,'-1'))
    timing <- read.csv(file.path(path,'timing.csv'))
    settings <- read.csv(file.path(path,'ld-swap-settings.csv'))
    stopifnot(timing$status=='completed',timing$threads==1L,
      settings$probability==item$probability,settings$minimum_r2==.8,
      settings$maximum_friends==50L,settings$moves==1L)
    value <- readRDS(file.path(path,'output.rds'))
    fit <- if(stage=='cmp_suite_regional') value$full_fit else value
    stopifnot(fit$control$ld_swap$probability==item$probability)
    key <- paste(stage,item$probability)
    if(is.null(reference[[key]])) reference[[key]] <- list(estimates=fit$estimates,swap=fit$ld_swap)
    else stopifnot(identical(reference[[key]]$estimates,fit$estimates),
      identical(reference[[key]]$swap,fit$ld_swap))
    counts <- read.csv(file.path(path,'ld-swap-counts.csv'))
    counts$stage <- stage; counts$probability <- item$probability; counts$repeat_id <- item[['repeat']]
    swaps[[length(swaps)+1L]] <- counts
    memory <- process[process$stage==stage,]
    measures[[length(measures)+1L]] <- data.frame(repeat_id=item[['repeat']],
      probability=item$probability,stage,analysis_seconds=timing$analysis_seconds,
      serialization_seconds=timing$serialization_seconds,
      process_wall_seconds=memory$processWallSeconds,
      peak_rss_gib=memory$peakWorkingSetBytes/1024^3,
      peak_private_gib=memory$sampledPeakPrivateBytes/1024^3,
      attempted=sum(counts$attempted),accepted=sum(counts$accepted))
    if(stage=='cmp_suite_regional') {
      phases <- read.csv(file.path(path,'regional-times.csv'))
      phases$probability <- item$probability; phases$repeat_id <- item[['repeat']]
      phases$regional_seconds <- phases$ld_total_seconds+phases$model_and_sets_seconds
      regions[[length(regions)+1L]] <- phases
    }
    conv <- fit$conv
    if(!'region' %in% names(conv)) conv$region <- 'genome'
    if(!'trait' %in% names(conv)) conv$trait <- rownames(conv)
    conv$stage <- stage; conv$probability <- item$probability; conv$repeat_id <- item[['repeat']]
    # Common fields only: regional convergence adds eigen-specific columns.
    diagnostics[[length(diagnostics)+1L]] <- conv[,c('stage','region','trait',
      'probability','repeat_id','nrun','draws',grep('^(rhat|ess|mcse|status)_',names(conv),value=TRUE))]
    est <- fit$estimates
    for(j in seq_along(frozen$causal)) {
      at <- est$marker==frozen$causal[j]
      if(stage=='cmp_suite_regional') at <- at & est$region==names(frozen$sets)[j]
      stopifnot(sum(at)==1L,is.finite(est$pip[at]))
      causal[[length(causal)+1L]] <- data.frame(stage,probability=item$probability,
        repeat_id=item[['repeat']],region=names(frozen$sets)[j],
        marker=frozen$causal[j],pip=est$pip[at])
    }
  }
}
measures <- do.call(rbind,measures);regions <- do.call(rbind,regions)
# Regional counters have region/trait labels; give genome counters the same labels.
swaps <- lapply(swaps,function(x) {
  if(!'region' %in% names(x)) x$region <- 'genome'
  if(!'trait' %in% names(x)) x$trait <- 'trait'
  x[,c('stage','probability','repeat_id','region','trait','attempted','accepted','acceptance_rate')]
})
swaps <- do.call(rbind,swaps);causal <- do.call(rbind,causal)
diagnostics <- do.call(rbind,diagnostics)
write.csv(runs,file.path(destination,'runs.csv'),row.names=FALSE)
for(name in c('measures','regions','swaps','causal','diagnostics'))
  write.csv(get(name),file.path(destination,paste0(name,'.csv')),row.names=FALSE)
summary_rows <- lapply(split(measures,list(measures$stage,measures$probability),drop=TRUE),function(x)
  data.frame(stage=x$stage[1],probability=x$probability[1],repeats=nrow(x),
    median_seconds=median(x$analysis_seconds),minimum_seconds=min(x$analysis_seconds),
    maximum_seconds=max(x$analysis_seconds),median_peak_rss_gib=median(x$peak_rss_gib),
    attempted=x$attempted[1],accepted=x$accepted[1]))
summary <- do.call(rbind,summary_rows)
write.csv(summary,file.path(destination,'summary.csv'),row.names=FALSE)
paired <- do.call(rbind,lapply(split(measures,list(measures$stage,measures$repeat_id),drop=TRUE),function(x) {
  off <- x[x$probability==0,];on <- x[x$probability==.05,]
  stopifnot(nrow(off)==1L,nrow(on)==1L)
  data.frame(stage=off$stage,repeat_id=off$repeat_id,off_seconds=off$analysis_seconds,
    on_seconds=on$analysis_seconds,difference_seconds=on$analysis_seconds-off$analysis_seconds,
    relative_difference_percent=100*(on$analysis_seconds/off$analysis_seconds-1))
}))
write.csv(paired,file.path(destination,'paired.csv'),row.names=FALSE)
region_summary <- do.call(rbind,lapply(split(regions,list(regions$region,regions$probability),drop=TRUE),function(x)
  data.frame(region=x$region[1],probability=x$probability[1],markers=x$markers[1],
    median_total_seconds=median(x$regional_seconds),median_ld_seconds=median(x$ld_total_seconds),
    median_prepare_seconds=median(x$prepare_seconds),
    median_binding_seconds=median(x$binding_seconds),
    median_native_fit_seconds=median(x$sampling_seconds),median_sets_seconds=median(x$sets_seconds))))
write.csv(region_summary,file.path(destination,'region-summary.csv'),row.names=FALSE)
pip_comparison <- do.call(rbind,lapply(c('cmp_suite_regional','cmp_suite_genome'),function(stage) {
  off <- reference[[paste(stage,0)]]$estimates
  on <- reference[[paste(stage,.05)]]$estimates
  stopifnot(identical(off$marker,on$marker))
  data.frame(stage,pip_correlation=cor(off$pip,on$pip),
    mean_absolute_pip_difference=mean(abs(off$pip-on$pip)),
    maximum_absolute_pip_difference=max(abs(off$pip-on$pip)))
}))
write.csv(pip_comparison,file.path(destination,'pip-comparison.csv'),row.names=FALSE)
lines <- c('# LD-swap timings on the frozen simulated-human example','',
  sprintf('%d serial paired repeats with alternating off/on order. One thread for OpenMP and BLAS,',
    length(unique(runs[['repeat']]))),
  'the same inputs, seeds, priors and iteration counts in each pair. Swaps off uses probability',
  'zero; swaps on uses probability 0.05, minimum r2 0.8, maximum 50 friends and up to one proposal',
  'per selected sweep. Regional fits have 100 burn-in and 1000 retained sweeps (seed 10);',
  'the genome fit has 100 sweeps and no burn-in (seed 20261005). Each fit has one chain.','',
  '| Workflow | Swaps | Median operation | Range | Median peak RSS | Attempts / accepted |',
  '| --- | --- | ---: | ---: | ---: | ---: |')
for(stage in c('cmp_suite_regional','cmp_suite_genome')) for(probability in c(0,.05)) {
  row <- summary[summary$stage==stage & summary$probability==probability,]
  lines <- c(lines,sprintf('| %s | %s | %.3f s | %.3f–%.3f s | %.3f GiB | %d / %d |',
    if(stage=='cmp_suite_regional') 'gmap, four dense regions' else 'gbayes, sparse genome',
    if(probability==0) 'off' else '5%',row$median_seconds,row$minimum_seconds,
    row$maximum_seconds,row$median_peak_rss_gib,row$attempted,row$accepted))
}
lines <- c(lines,'',
  'Runtime varied across repeats, including the LD phases that precede swapping.',
  'Whole-operation median differences therefore do not isolate the swap cost.',
  'paired.csv preserves each matched comparison; region-summary.csv includes native-fit phases.',
  'Do not infer a speed improvement or a precise percentage overhead from these short pairs.',
  'An untimed LD audit can be run afterwards with audit_ld_swap_example.R. Zero proposal counts',
  'must not be interpreted as evidence of mixing improvement; they can indicate weak LD.',
  '', '| Region | Markers | Off median | On median | On attempts / accepted |',
  '| --- | ---: | ---: | ---: | ---: |')
for(region in names(frozen$sets)) {
  off <- region_summary[region_summary$region==region & region_summary$probability==0,]
  on <- region_summary[region_summary$region==region & region_summary$probability==.05,]
  count <- swaps[swaps$region==region & swaps$probability==.05,][1,]
  lines <- c(lines,sprintf('| %s | %d | %.3f s | %.3f s | %d / %d |',region,off$markers,
    off$median_total_seconds,on$median_total_seconds,count$attempted,count$accepted))
}
lines <- c(lines,'',
  'Operation time excludes process startup, package/input loading and RDS serialization.',
  'Regional operation time includes reading the existing BED, computing full regional LD,',
  'preparing its dense likelihood, native fitting, posterior tables and credible groups.',
  'Genome sparse LD preparation is reused and excluded. Neither workflow recomputes genotype LD',
  'for swap proposals. The regional native-fit phase includes one-time LD-friend preparation',
  'when enabled; it is not a separate pure Gibbs-loop measurement. Phase totals omit some',
  'R wrapper work, so use whole-operation timings for workflow comparisons.',
  'Process peak RSS includes package loading and is measured from the Windows R process.',
  'Full per-repeat timings, memory, phases and swap counts are in the adjacent CSV files.','',
  sprintf('The four known causal-marker PIPs range from %.3f to %.3f across both workflows and settings.',
    min(causal$pip),max(causal$pip)),
  'Fixed seeds reproduce estimates and swap counts exactly within each setting across repeats.',
  'Enabling swaps changes RNG consumption and can change other PIPs. Full PIP comparisons are',
  'recorded separately. Native single-chain diagnostic statuses and ESS/MCSE values are retained',
  'in diagnostics.csv and the fit RDS files. These short examples do not establish a mixing',
  'improvement or convergence, especially the 100-sweep genome fit. Swaps remain opt-in.','',
  'These are fresh controls against the current installed packages and the corrected native',
  'reverse-edge check. Earlier qgg comparisons and source/evidence archives are unchanged.',
  'This report does not replace historical timings or assert numerical equivalence to qgg.','',
  'Reproduce from the gsuite checkout with:',
  '`powershell.exe -NoProfile -ExecutionPolicy Bypass -File benchmarks/genomic_20k_50k/time_ld_swap.ps1`.',
  'Then run `Rscript --vanilla benchmarks/genomic_20k_50k/summarize_ld_swap.R <runs.csv> <new-report-directory>`.',
  'The scripts reuse the preserved frozen example and prepared sparse-LD cache under build/benchmarks.',
  'Run `Rscript --vanilla benchmarks/genomic_20k_50k/audit_ld_swap_example.R <report-directory>`',
  'after the timing processes finish to reproduce the untimed LD audit.',
  'The build/evidence directory records source commits/diffs, installed-package SHA256 hashes,',
  'benchmark scripts and the run inventory; each run includes hardware, input hashes, session',
  'information, effective controls and full fitted results. No download, QC or simulation is repeated.',
  paste('Evidence directory:',normalizePath(dirname(args[1]),winslash='/')), '',
  'No commit or push was performed.')
writeLines(lines,file.path(destination,'README.md'))
print(summary,row.names=FALSE)
print(region_summary,row.names=FALSE)
print(pip_comparison,row.names=FALSE)
