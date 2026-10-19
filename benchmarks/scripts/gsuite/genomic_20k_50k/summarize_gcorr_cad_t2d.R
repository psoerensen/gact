# Summarize completed process-measured gcorr runs; no fitting or data acquisition.
# Arguments: pilot run, score run, regional run, validated input-stage directory, output directory.
args <- commandArgs(TRUE)
stopifnot(length(args)==5L)
pilot <- args[1];genome <- args[2];regional_root <- args[3];input_stage <- args[4];destination <- args[5]
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
assets <- file.path(destination,'data');dir.create(assets,showWarnings=FALSE)
read <- function(root,stage,file) read.csv(file.path(root,paste0(stage,'-4'),file),check.names=FALSE)
markdown <- function(x) {
  x[] <- lapply(x,function(v) {
    if(is.numeric(v) || (is.logical(v) && all(is.na(v))))
      ifelse(is.na(v),'unavailable',formatC(as.numeric(v),digits=5,format='g'))
    else ifelse(is.na(v),'',as.character(v))
  })
  paste(c(paste0('| ',paste(names(x),collapse=' | '),' |'),
    paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
    apply(x,1,function(r)paste0('| ',paste(r,collapse=' | '),' |'))),collapse='\n')
}
csv <- function(x,name)write.csv(x,file.path(assets,name),row.names=FALSE,na='NA')
counts <- read.csv(file.path(input_stage,'input-counts.csv'))
samples <- read.csv(file.path(input_stage,'sample-sizes.csv'))
csv(counts,'gcorr-input-counts.csv');csv(samples,'gcorr-sample-sizes.csv')
global_parts <- lapply(c('gcorr_ldsc','gcorr_gnova','gcorr_sumher'),function(stage) {
  x <- read(genome,stage,'estimates.csv');x$method <- sub('gcorr_','',stage);x
})
global_columns <- Reduce(union,lapply(global_parts,names))
global <- do.call(rbind,lapply(global_parts,function(x) {
  for(name in setdiff(global_columns,names(x)))x[[name]] <- NA
  x[,global_columns,drop=FALSE]
}))
csv(global,'gcorr-genome-estimates.csv')
compact <- global[global$quantity %in% c('h2','cov','rg') &
  (is.na(global$component) | global$component=='all'),
  intersect(c('method','quantity','trait1','trait2','estimate','standard_error'),names(global)),drop=FALSE]
# Keep raw CSV names authoritative; table selection adapts only to existing output labels.
if(!nrow(compact))compact <- global[,intersect(c('method','quantity','trait1','trait2','estimate','standard_error'),names(global)),drop=FALSE]
regional <- read(regional_root,'gcorr_genome','estimates.csv')
pilot_estimates <- read(pilot,'gcorr_chr22','estimates.csv')
sensitivity <- read(regional_root,'gcorr_genome','sensitivity-estimates.csv')
marginal <- read(regional_root,'gcorr_genome','marginal-estimates.csv')
phases <- read(regional_root,'gcorr_genome','regional-phases.csv')
csv(regional,'gcorr-regional-estimates.csv');csv(marginal,'gcorr-marginal-hess-estimates.csv')
csv(pilot_estimates,'gcorr-chr22-estimates.csv');csv(sensitivity,'gcorr-overlap-sensitivity.csv')
csv(phases,'gcorr-regional-preparation.csv')
score_process <- read.csv(file.path(genome,'process-measurements.csv'))
score_process <- score_process[score_process$stage %in% c('gcorr_ldsc','gcorr_gnova','gcorr_sumher'),]
regional_process <- read.csv(file.path(regional_root,'process-measurements.csv'))
regional_process <- regional_process[regional_process$stage=='gcorr_genome',]
pilot_process <- read.csv(file.path(pilot,'process-measurements.csv'))
pilot_process <- pilot_process[pilot_process$stage=='gcorr_chr22',]
stopifnot(nrow(pilot_process)==1L,nrow(score_process)==3L,nrow(regional_process)==1L)
process <- rbind(pilot_process,score_process,regional_process)
stopifnot(all(process$exitCode==0),
  all(tolower(as.character(process$timedOut))=='false'))
# The Windows launcher records decimal values using the local culture.
for(name in c('processWallSeconds','rProcessCpuSeconds')) {
  process[[name]] <- as.numeric(sub(',', '.', as.character(process[[name]]),fixed=TRUE))
  stopifnot(all(is.finite(process[[name]])))
}
process$path <- NULL
process$peak_working_set_GiB <- process$peakWorkingSetBytes/1024^3
process$sampled_peak_private_GiB <- process$sampledPeakPrivateBytes/1024^3
csv(process,'gcorr-process-measurements.csv')
timings <- do.call(rbind,lapply(seq_len(nrow(process)),function(i) {
  p <- process[i,];root <- if(p$stage=='gcorr_chr22')pilot else if(p$stage=='gcorr_genome')regional_root else genome
  directory <- file.path(root,paste0(p$stage,'-',p$threads))
  t <- read.csv(file.path(directory,'timing.csv'))
  prep <- if(p$stage %in% c('gcorr_chr22','gcorr_genome'))
    read.csv(file.path(directory,'preparation-phase.csv'))$seconds else NA_real_
  fit <- if(is.na(prep))read.csv(file.path(directory,'fit-phase.csv'))$seconds else
    sum(vapply(c('marginal_fit','rho_fit','sensitivity_fit'),function(phase)
      read.csv(file.path(directory,paste0(phase,'-phase.csv')))$seconds,numeric(1)))
  data.frame(analysis=p$stage,threads=p$threads,preparation_seconds=prep,
    fit_seconds=fit,analysis_seconds=t$analysis_seconds,
    process_wall_seconds=p$processWallSeconds,serialization_seconds=t$serialization_seconds,
    peak_working_set_GiB=p$peak_working_set_GiB,sampled_peak_private_GiB=p$sampled_peak_private_GiB)
}))
csv(timings,'gcorr-timings.csv')
fit_breakdown <- do.call(rbind,lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  root <- if(stage=='gcorr_chr22')pilot else regional_root
  do.call(rbind,lapply(c('marginal_fit','rho_fit','sensitivity_fit'),function(phase) {
    data.frame(analysis=stage,phase=phase,
      seconds=read(root,stage,paste0(phase,'-phase.csv'))$seconds)
  }))
}))
csv(fit_breakdown,'gcorr-fit-breakdown.csv')
breakdown <- read(regional_root,'gcorr_genome','preparation-overhead.csv')
preparation_total <- breakdown$seconds[breakdown$phase=='total']
phase_totals <- colSums(phases[,c('decode','product','eigen','projection')])
breakdown <- rbind(breakdown[breakdown$phase!='total',],
  data.frame(phase=names(phase_totals),seconds=unname(phase_totals)))
breakdown <- rbind(breakdown,data.frame(phase='other_checks_and_adapters',
  seconds=max(0,preparation_total-sum(breakdown$seconds))))
csv(breakdown,'gcorr-preparation-breakdown.csv')

png(file.path(assets,'gcorr-preparation-time.png'),width=1400,height=760,res=150)
phase_labels <- c(validation='R checks/matching',panel_setup='BED reader setup',
  decode='Genotype decoding',product='Gram matrix product',eigen='Eigen decomposition',
  projection='Both-trait projection',other_checks_and_adapters='Other checks/adapters')
par(mar=c(10,4,3,1));barplot(breakdown$seconds,names.arg=phase_labels[breakdown$phase],
  las=2,col='#247589',border=NA,ylab='Elapsed seconds',
  main='Full-genome regional preparation: one shared decomposition per region')
dev.off()
point <- regional[regional$component=='region' & regional$quantity=='cov',]
point$chr <- as.integer(sub('chr([0-9]+)_.*','\\1',point$region))
point$start <- as.numeric(sub('chr[0-9]+_([0-9]+)_.*','\\1',point$region))
point <- point[order(point$chr,point$start),]
png(file.path(assets,'gcorr-regional-covariance.png'),width=1600,height=760,res=150)
par(mar=c(4.5,5,3,1));plot(seq_len(nrow(point)),point$estimate,
  pch=16,cex=.55,col=c('#247589','#9b582c')[1+point$chr%%2],
  xaxt='n',xlab='Chromosome (regions in coordinate order)',
  ylab='Local genetic covariance (working scale)',
  main='CAD/T2D joint rho-HESS: labelled zero-overlap scenario')
abline(h=0,col='grey55',lty=2)
axis(1,at=tapply(seq_len(nrow(point)),point$chr,mean),labels=sort(unique(point$chr)))
dev.off()

totals <- regional[regional$component=='total',
  c('quantity','trait1','trait2','estimate','standard_error'),drop=FALSE]
missing_se <- sum(is.na(regional$standard_error))
negative <- sum(regional$component=='region' & regional$quantity=='h2' & regional$estimate<0,na.rm=TRUE)
raw_ratio <- totals$estimate[totals$quantity=='rg']
ratio_note <- if(length(raw_ratio)==1L && is.finite(raw_ratio) && abs(raw_ratio)>1)
  sprintf('**The zero-overlap total ratio is %.3f, outside [-1,1]. Its estimated total genetic covariance matrix is indefinite. This ratio cannot be interpreted as a biological genetic correlation or used for conditional correlation inference.** Under the user-provided no-overlap assumption, sample-size scaling, ancestry/reference compatibility and spectral-model assumptions still require review before interpreting these regional fits.',raw_ratio) else ''
regional_text <- paste(ratio_note,markdown(totals),sprintf('\n\nThere are %d negative regional heritability estimates. %d of %d reported regional/total standard errors are unavailable under the conditional plug-in checks. Raw point estimates remain in the downloads; these values do not supply calibrated region-level significance or correlation inference.',negative,missing_se,nrow(regional)),
  '\n\nOverlap sensitivity (full-genome totals):\n',
  markdown(sensitivity[,c('scenario','quantity','estimate','standard_error')]),sep='\n')
display_timings <- timings[,c('analysis','threads','preparation_seconds','fit_seconds','analysis_seconds','peak_working_set_GiB')]
display_timings$analysis <- c(gcorr_chr22='HESS/rho-HESS chr22',gcorr_genome='HESS/rho-HESS genome',
  gcorr_ldsc='LDSC',gcorr_gnova='GNOVA',gcorr_sumher='SumHer baseline')[display_timings$analysis]
display_timings$preparation_seconds <- ifelse(is.na(display_timings$preparation_seconds),
  'not separate',formatC(display_timings$preparation_seconds,digits=2,format='f'))
names(display_timings) <- c('Analysis','Threads','Preparation (s)','Fits (s)','Analysis total (s)','Peak working set (GiB)')
timing_text <- paste(markdown(display_timings),
  '\n\nThe regional fit column includes both marginal fits, the zero-overlap pair fit and two additional overlap sensitivity fits. Preparation overhead and numerical phase totals are shown below. These completed runs reuse a validated input cache; input extraction and saved-input loading are outside the analysis clocks.',
  sprintf('\n\nFor the full genome, the two marginal fits took %.2f seconds, the zero-overlap pair fit %.2f seconds, and the additional sensitivity fits %.2f seconds. These clocks include R adapters and result formatting.',
    fit_breakdown$seconds[fit_breakdown$analysis=='gcorr_genome' & fit_breakdown$phase=='marginal_fit'],
    fit_breakdown$seconds[fit_breakdown$analysis=='gcorr_genome' & fit_breakdown$phase=='rho_fit'],
    fit_breakdown$seconds[fit_breakdown$analysis=='gcorr_genome' & fit_breakdown$phase=='sensitivity_fit']),
  '\n\n',markdown(breakdown),sep='\n')
files <- list.files(assets,pattern='[.]csv$')
downloads <- paste(paste0('- [',files,'](data/',files,')'),collapse='\n')
downloads <- paste(downloads,
  '- [Runnable preparation helper](data/gcorr_cad_t2d_example.R)',
  '- [Disjoint marker-set score helper](data/gcorr_set_scores.R)',
  '- [Pinned EUR regional boundaries](data/EUR-fourier-ac125e47.bed)',
  '- [Analysis and timing protocol](data/gcorr-cad-t2d-protocol.md)',sep='\n')
display_counts <- counts
display_counts$qc <- c(reference_markers='Reference markers',regression_markers='Common QC markers',
  regional_markers='Regional markers',regions='Nonempty LD regions',unpartitioned_markers='Uncovered markers')[display_counts$qc]
names(display_counts) <- c('Input','Count')
display_counts$Count <- formatC(display_counts$Count,format='d',big.mark=',')
replacements <- list(INPUT_COUNTS=markdown(display_counts),GENOME_RESULTS=markdown(compact),
  REGIONAL_RESULTS=regional_text,TIMINGS=timing_text,
  PLOTS='![Preparation time by phase](data/gcorr-preparation-time.png)\n\n![Regional covariance, zero-overlap scenario](data/gcorr-regional-covariance.png)',
  DOWNLOADS=downloads)
saveRDS(replacements,file.path(destination,'tutorial-sections.rds'))
cat('EXPORTED',destination,'\n')
