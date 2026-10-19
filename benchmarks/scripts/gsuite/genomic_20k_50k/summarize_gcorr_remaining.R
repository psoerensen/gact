# Export the completed corrected GNOVA/SumHer and HESS workloads; no fitting.
# Args: corrected LDSC/input run, score/chr22 run, genome HESS run, tutorial dir.
args <- commandArgs(TRUE);stopifnot(length(args)==4L)
base_run <- args[1];score_run <- args[2];genome_run <- args[3];tutorial <- args[4]
evidence <- 'benchmarks/genomic_20k_50k/manuscript-results/gcorr-cad-t2d-20261007/remaining-analyses'
dir.create(evidence,recursive=TRUE,showWarnings=FALSE)
assets <- file.path(tutorial,'data');dir.create(assets,recursive=TRUE,showWarnings=FALSE)
read <- function(root,stage,name)read.csv(file.path(root,paste0(stage,'-4'),name),check.names=FALSE)
write <- function(x,name) {
  write.csv(x,file.path(evidence,name),row.names=FALSE,na='NA')
  write.csv(x,file.path(assets,name),row.names=FALSE,na='NA')
}
markdown <- function(x) {
  x[] <- lapply(x,function(v)if(is.numeric(v) || is.logical(v))
    ifelse(is.na(v),'unavailable',formatC(as.numeric(v),digits=6,format='g'))
    else ifelse(is.na(v),'',as.character(v)))
  paste(c(paste0('| ',paste(names(x),collapse=' | '),' |'),
    paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
    apply(x,1,function(r)paste0('| ',paste(r,collapse=' | '),' |'))),collapse='\n')
}
join <- function(parts) {
  fields <- Reduce(union,lapply(parts,names))
  do.call(rbind,lapply(parts,function(x) {
    for(name in setdiff(fields,names(x)))x[[name]] <- NA
    x[,fields,drop=FALSE]
  }))
}
roots <- c(gcorr_gnova=score_run,gcorr_sumher=score_run,
           gcorr_chr22=score_run,gcorr_genome=genome_run)
process <- join(lapply(names(roots),function(stage) {
  p <- read.csv(file.path(roots[[stage]],'process-measurements.csv'))
  p[p$stage==stage,,drop=FALSE]
}))
stopifnot(nrow(process)==4L,all(process$exitCode==0),
  all(tolower(as.character(process$timedOut))=='false'))
for(name in c('processWallSeconds','rProcessCpuSeconds'))
  process[[name]] <- as.numeric(sub(',','.',as.character(process[[name]]),fixed=TRUE))
process$path <- NULL
process$peak_working_set_GiB <- process$peakWorkingSetBytes/1024^3
process$sampled_peak_private_GiB <- process$sampledPeakPrivateBytes/1024^3
write(process,'gcorr-remaining-process-measurements.csv')
audits <- join(lapply(names(roots),function(stage) {
  x <- read(roots[[stage]],stage,'analysis-input-audit.csv');x$analysis <- stage;x
}))
expected_n <- c(GWAS1=49071.2732949614,GWAS2=40743.1524049809)
stopifnot(all(audits$markers==7221923),all(audits$n_min==audits$n_max),
  all(abs(audits$n_min-expected_n[audits$trait])<1e-8))
write(audits,'gcorr-remaining-input-audit.csv')
intercept <- read(score_run,'gcorr_gnova','ldsc-intercept-source.csv')
stopifnot(intercept$md5==unname(tools::md5sum(
  file.path(base_run,'gcorr_ldsc-4','estimates.csv'))))
# Export the value/checksum without a host-specific absolute source path.
intercept$source <- paste(basename(base_run),'gcorr_ldsc-4/estimates.csv',sep='/')
write(intercept,'gcorr-gnova-intercept-source.csv')
global <- join(lapply(c('gcorr_ldsc','gcorr_gnova','gcorr_sumher'),function(stage) {
  x <- read(if(stage=='gcorr_ldsc')base_run else score_run,stage,'estimates.csv')
  x$method <- sub('gcorr_','',stage);x
}))
write(global,'gcorr-current-score-estimates.csv')
selected <- global$quantity %in% c('h2','cov','rg') &
  (is.na(global$component) | global$component=='all')
global_display <- global[selected,c('method','quantity','trait1','trait2','estimate','standard_error')]
regional <- join(lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- read(roots[[stage]],stage,'estimates.csv')
  x$analysis <- stage;x
}))
write(regional,'gcorr-current-hess-estimates.csv')
marginals <- join(lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- read(roots[[stage]],stage,'marginal-estimates.csv');x$analysis <- stage;x
}))
write(marginals,'gcorr-current-marginal-hess-estimates.csv')
rank <- join(lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- read(roots[[stage]],stage,'rank-admission.csv');x$analysis <- stage;x
}))
stopifnot(all(rank$correction_denominator>0))
write(rank,'gcorr-current-rank-admission.csv')
configuration <- join(lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- read(roots[[stage]],stage,'projection-configuration.csv');x$analysis <- stage;x
}))
write(configuration,'gcorr-current-hess-configuration.csv')
phases <- read(genome_run,'gcorr_genome','regional-phases.csv')
stopifnot(nrow(phases)==1702L,sum(phases$markers)==7221923,
  sum(phases$retained_components)==unique(rank$retained_rank[rank$analysis=='gcorr_genome']))
write(phases,'gcorr-current-regional-preparation.csv')
spectral <- data.frame(summary=c('minimum','median','maximum'),
  retained_reference_eigenvalue_fraction=c(min(phases$retained_eigenvalue_fraction),
    median(phases$retained_eigenvalue_fraction),max(phases$retained_eigenvalue_fraction)))
stopifnot(all(is.finite(spectral$retained_reference_eigenvalue_fraction)))
write(spectral,'gcorr-current-spectral-retention.csv')
totals <- regional[regional$component=='total',
  c('analysis','quantity','trait1','trait2','estimate','standard_error')]
sum_checks <- do.call(rbind,lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- regional[regional$analysis==stage & regional$quantity %in% c('h2','cov'),]
  full <- x[x$component=='total',]
  do.call(rbind,lapply(seq_len(nrow(full)),function(i) {
    matching <- x$component=='region' & x$quantity==full$quantity[i] &
      x$trait1==full$trait1[i] & (if(is.na(full$trait2[i]))is.na(x$trait2) else x$trait2==full$trait2[i])
    summed <- sum(x$estimate[matching])
    data.frame(analysis=stage,quantity=full$quantity[i],trait1=full$trait1[i],
      trait2=full$trait2[i],region_sum=summed,total=full$estimate[i],
      difference=summed-full$estimate[i])
  }))
}))
stopifnot(all(is.finite(sum_checks$difference)),all(abs(sum_checks$difference)<1e-8))
write(sum_checks,'gcorr-current-hess-sum-checks.csv')
timings <- do.call(rbind,lapply(names(roots),function(stage) {
  p <- process[process$stage==stage,];t <- read(roots[[stage]],stage,'timing.csv')
  hess <- stage %in% c('gcorr_chr22','gcorr_genome')
  data.frame(analysis=stage,threads=4,
    preparation_seconds=if(hess)read(roots[[stage]],stage,'preparation-phase.csv')$seconds else NA_real_,
    fit_seconds=if(hess)sum(vapply(c('marginal_fit','rho_fit'),function(phase)
      read(roots[[stage]],stage,paste0(phase,'-phase.csv'))$seconds,numeric(1)))
      else read(roots[[stage]],stage,'fit-phase.csv')$seconds,
    analysis_seconds=t$analysis_seconds,serialization_seconds=t$serialization_seconds,
    process_wall_seconds=p$processWallSeconds,peak_working_set_GiB=p$peak_working_set_GiB,
    sampled_peak_private_GiB=p$sampled_peak_private_GiB)
}))
write(timings,'gcorr-remaining-timings.csv')
breakdown <- read(genome_run,'gcorr_genome','preparation-overhead.csv')
preparation_total <- breakdown$seconds[breakdown$phase=='total']
phase_totals <- colSums(phases[,c('decode','product','eigen','projection')])
breakdown <- rbind(breakdown[breakdown$phase!='total',],
  data.frame(phase=names(phase_totals),seconds=unname(phase_totals)))
breakdown <- rbind(breakdown,data.frame(phase='other_checks_and_adapters',
  seconds=max(0,preparation_total-sum(breakdown$seconds))))
write(breakdown,'gcorr-current-preparation-breakdown.csv')
pilot_phases <- read(score_run,'gcorr_chr22','regional-phases.csv')
write(pilot_phases,'gcorr-current-chr22-preparation.csv')
pilot_overhead <- read(score_run,'gcorr_chr22','preparation-overhead.csv')
pilot_sums <- colSums(pilot_phases[,c('decode','product','eigen','projection')])
pilot_breakdown <- rbind(pilot_overhead[pilot_overhead$phase!='total',],
  data.frame(phase=names(pilot_sums),seconds=unname(pilot_sums)))
write(pilot_breakdown,'gcorr-current-chr22-preparation-breakdown.csv')
quality <- do.call(rbind,lapply(c('gcorr_chr22','gcorr_genome'),function(stage) {
  x <- regional[regional$analysis==stage,]
  data.frame(analysis=stage,rows=nrow(x),
    negative_regional_h2=sum(x$component=='region' & x$quantity=='h2' & x$estimate<0,na.rm=TRUE),
    unavailable_points=sum(!is.finite(x$estimate)),unavailable_standard_errors=sum(!is.finite(x$standard_error)),
    rg_outside_unit_interval=sum(x$quantity=='rg' & is.finite(x$estimate) & abs(x$estimate)>1,na.rm=TRUE))
}))
write(quality,'gcorr-current-hess-availability.csv')
plot_name <- 'gcorr-current-preparation-time.png'
png(file.path(assets,plot_name),width=1400,height=760,res=150)
labels <- c(validation='R checks/matching',panel_setup='BED reader setup',decode='Genotype decoding',
  product='Gram product',eigen='Eigen decomposition',projection='Both-trait projection',
  other_checks_and_adapters='Other checks/adapters')
par(mar=c(10,4,3,1));barplot(breakdown$seconds,names.arg=labels[breakdown$phase],las=2,
  col='#247589',border=NA,ylab='Elapsed seconds',
  main='Corrected genome HESS preparation: shared decomposition, K <= 20')
invisible(dev.off());file.copy(file.path(assets,plot_name),file.path(evidence,plot_name),overwrite=TRUE)
cov_name <- 'gcorr-current-regional-covariance.png'
point <- regional[regional$analysis=='gcorr_genome' & regional$component=='region' & regional$quantity=='cov',]
point$chr <- as.integer(sub('chr([0-9]+)_.*','\\1',point$region))
point$start <- as.numeric(sub('chr[0-9]+_([0-9]+)_.*','\\1',point$region))
point <- point[order(point$chr,point$start),]
png(file.path(assets,cov_name),width=1600,height=760,res=150)
par(mar=c(4.5,5,3,1));plot(seq_len(nrow(point)),point$estimate,pch=16,cex=.55,
  col=c('#247589','#9b582c')[1+point$chr%%2],xaxt='n',
  xlab='Chromosome (regions in coordinate order)',ylab='Local genetic covariance (supplied n convention)',
  main='Corrected joint rho-HESS: no overlap, K <= 20')
abline(h=0,col='grey55',lty=2)
axis(1,at=tapply(seq_len(nrow(point)),point$chr,mean),
  labels=sort(unique(point$chr)),las=2,cex.axis=.8,gap.axis=0)
invisible(dev.off());file.copy(file.path(assets,cov_name),file.path(evidence,cov_name),overwrite=TRUE)
reasons <- unique(regional$uncertainty_reason[nzchar(regional$uncertainty_reason) &
  !is.na(regional$uncertainty_reason)])
regional_text <- paste('Completed corrected HESS/rho-HESS totals:',markdown(totals),
  'Chr22 uses K <= 50 and defines a chromosome-only modeled universe. Genome uses K <= 20 across all 1,702 regions. These are different spectral/model specifications.',
  'Actual retained ranks and correction denominators:',markdown(rank),
  'Reference LD spectral fraction retained by the genome cap (not a measured fraction of genetic signal):',markdown(spectral),
  'Estimate/uncertainty availability (raw points are preserved):',markdown(quality),
  if(length(reasons))paste('Native reasons for unavailable uncertainty:',paste(reasons,collapse='; ')) else '',
  'All regional h2/covariance contributions sum to their corresponding joint totals. Correlations are ratios and do not sum. A finite point estimate is not evidence of calibrated uncertainty.',sep='\n\n')
full_rg <- totals$estimate[totals$analysis=='gcorr_genome' & totals$quantity=='rg']
if(length(full_rg)==1L && !is.finite(full_rg))
  regional_text <- paste(regional_text,
    '**The corrected CAD genome total is negative, so genome rg is undefined. All HESS/rho-HESS standard errors are unavailable under the native conditional checks. These results do not support interpretable CAD/T2D genetic-correlation inference under this configuration.**',sep='\n\n')
if(length(full_rg)==1L && is.finite(full_rg) && abs(full_rg)>1)
  regional_text <- paste(regional_text,sprintf('**The genome total raw rg is %.3f, outside [-1,1]; it is not interpretable as a biological genetic correlation.**',full_rg),sep='\n\n')
sections <- readRDS(file.path(tutorial,'tutorial-sections.rds'))
sections$GENOME_RESULTS <- paste('Completed corrected score fits:',markdown(global_display),sep='\n\n')
sections$REGIONAL_RESULTS <- regional_text
old_timings <- read.csv(file.path(assets,'gcorr-corrected-timings.csv'))
old_timings$preparation_seconds <- NA_real_
old_process <- read.csv(file.path(assets,'gcorr-corrected-process-measurements.csv'))
old_timings$sampled_peak_private_GiB <- old_process$sampled_peak_private_GiB[
  match(old_timings$analysis,old_process$stage)]
timing_columns <- c('analysis','preparation_seconds','fit_seconds','analysis_seconds',
  'process_wall_seconds','peak_working_set_GiB','sampled_peak_private_GiB')
display_timings <- rbind(old_timings[,timing_columns],timings[,timing_columns])
display_timings$analysis <- c(gcorr_inputs='Input preparation',gcorr_ldsc='LDSC',
  gcorr_sets='Chromosome-set LDSC',gcorr_gnova='GNOVA',gcorr_sumher='SumHer baseline',
  gcorr_chr22='HESS/rho-HESS chr22 (cap 50)',gcorr_genome='HESS/rho-HESS genome (cap 20)')[display_timings$analysis]
names(display_timings) <- c('Analysis','Prepare (s)','Fit (s)','Analysis (s)',
  'Wall (s)','Peak WS (GiB)','Peak private (GiB)')
sections$CORRECTION_TIMINGS <- paste(markdown(display_timings),
  'HESS fits include two marginal fits and one no-overlap pair fit; the duplicate zero-overlap sensitivity table is an export of the same fit. Reading the input cache is outside analysis time but included in process wall time.',
  'Working set records resident R-process memory; sampled private memory records committed process allocation and can exceed physical RAM. Runs on this 16 GiB machine can therefore be affected by paging.',
  sprintf('In the chr22 pilot, validation took %.2f s and BED reader setup %.2f s; eigen decomposition across all 24 regions took %.2f s. Setup dominates this small batch.',
    pilot_breakdown$seconds[pilot_breakdown$phase=='validation'],
    pilot_breakdown$seconds[pilot_breakdown$phase=='panel_setup'],pilot_sums[['eigen']]),
  'Corrected genome preparation breakdown:',markdown(breakdown),
  paste0('![Corrected preparation time](data/',plot_name,')'),
  paste0('![Corrected regional covariance](data/',cov_name,')'),sep='\n\n')
current_files <- setdiff(list.files(evidence,pattern='[.]csv$'),
  'remaining-source-and-evidence-files.csv')
download_lines <- strsplit(sections$DOWNLOADS,'\n',fixed=TRUE)[[1]]
sections$DOWNLOADS <- paste(download_lines[!grepl(
  'remaining-source-and-evidence-files.csv',download_lines,fixed=TRUE)],collapse='\n')
# Keep historical files available, but label every old timing/result link.
for(name in c('gcorr-genome-estimates.csv','gcorr-regional-estimates.csv',
  'gcorr-marginal-hess-estimates.csv','gcorr-chr22-estimates.csv',
  'gcorr-overlap-sensitivity.csv','gcorr-regional-preparation.csv',
  'gcorr-process-measurements.csv','gcorr-timings.csv','gcorr-fit-breakdown.csv',
  'gcorr-preparation-breakdown.csv')) {
  sections$DOWNLOADS <- gsub(paste0('[',name,']'),
    paste0('[',name,' (superseded total-N run)]'),sections$DOWNLOADS,fixed=TRUE)
}
sections$DOWNLOADS <- paste(sections$DOWNLOADS,
  paste(paste0('- [',current_files,'](data/',current_files,')'),collapse='\n'),sep='\n')
sections$DOWNLOADS <- paste(unique(strsplit(sections$DOWNLOADS,'\n',fixed=TRUE)[[1]]),collapse='\n')
saveRDS(sections,file.path(tutorial,'tutorial-sections.rds'))
writeLines(c('# Corrected remaining CAD/T2D gcorr analyses, 7 October 2026','',
  'All methods use the same 7,221,923-marker website availability mask and extracted stat$n.',
  'Stored LD scores are reused by GNOVA and SumHer. HESS streams BED regions and discards LD/eigenvectors.',
  'HESS uses K <= 50 for chr22 and an explicit K <= 20 for the genome; package defaults are unchanged.',
  'The lower genome cap is set from input rank admission, not fitted trait signal; it changes the spectral model.',
  'Results remain conditional on the database information-size convention and statistical restrictions described in the tutorial.',
  'The unchanged native/R gcorr 0.1.1 baseline is preserved in the neighboring 6 October source-and-evidence.zip; this record supplies corrected workers, inputs audits and results.',
  '',markdown(global_display),'',regional_text,'','## Timing and memory','',markdown(timings),
  '',paste0('![Corrected preparation time](',plot_name,')'),
  '',paste0('![Corrected regional covariance](',cov_name,')'),
  '',paste('Input/LDSC run:',basename(base_run)),paste('Score/chr22 run:',basename(score_run)),
  paste('Genome run:',basename(genome_run))),file.path(evidence,'README.md'))
cat('REMAINING_ANALYSES_EXPORTED\n');print(global_display);print(totals);print(quality)
