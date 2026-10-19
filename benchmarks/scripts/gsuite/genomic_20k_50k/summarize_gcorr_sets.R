# Export a completed marker-set run into compact tutorial evidence.
# Arguments: run directory, output directory.
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
run <- args[1L];destination <- args[2L]
assets <- file.path(destination,'data');dir.create(assets,recursive=TRUE,showWarnings=FALSE)
stage <- file.path(run,'gcorr_sets-4')
process <- read.csv(file.path(run,'process-measurements.csv'))
process <- process[process$stage=='gcorr_sets',]
stopifnot(nrow(process)==1L,process$exitCode==0,
  tolower(as.character(process$timedOut))=='false')
for(n in c('processWallSeconds','rProcessCpuSeconds'))
  process[[n]] <- as.numeric(sub(',','.',as.character(process[[n]]),fixed=TRUE))
process$path <- NULL
process$peak_working_set_GiB <- process$peakWorkingSetBytes/1024^3
process$sampled_peak_private_GiB <- process$sampledPeakPrivateBytes/1024^3
estimates <- read.csv(file.path(stage,'estimates.csv'))
counts <- read.csv(file.path(stage,'set-counts.csv'))
wide <- data.frame(set=counts$set,regression_markers=counts$regression_markers,
  reference_markers=counts$reference_markers,stringsAsFactors=FALSE)
select <- function(component,quantity,trait) {
  x <- estimates[estimates$component==component & estimates$quantity==quantity &
    estimates$trait1==trait,,drop=FALSE]
  stopifnot(nrow(x)==1L);x
}
for(i in seq_len(nrow(wide))) {
  a <- select(wide$set[i],'h2','GWAS1');b <- select(wide$set[i],'h2','GWAS2')
  c <- select(wide$set[i],'cov','GWAS1');r <- select(wide$set[i],'rg','GWAS1')
  wide$h2_T2D[i] <- a$estimate;wide$se_h2_T2D[i] <- a$standard_error
  wide$h2_CAD[i] <- b$estimate;wide$se_h2_CAD[i] <- b$standard_error
  wide$covariance[i] <- c$estimate;wide$se_covariance[i] <- c$standard_error
  wide$rg[i] <- r$estimate;wide$se_rg[i] <- r$standard_error
  wide$rg_reason[i] <- if(is.na(r$estimate))as.character(r$undefined_reason) else
    if(is.na(r$standard_error))as.character(r$uncertainty_reason) else ''
}
wide <- wide[order(as.integer(sub('chr','',wide$set,fixed=TRUE))),]
totals <- estimates[estimates$component=='all' & estimates$quantity %in% c('h2','cov','rg'),]
for(q in c('h2','cov')) {
  for(t in unique(totals$trait1[totals$quantity==q])) {
    total <- totals$estimate[totals$quantity==q & totals$trait1==t]
    parts <- estimates$estimate[estimates$component!='all' &
      estimates$quantity==q & estimates$trait1==t]
    stopifnot(length(total)==1L,abs(total-sum(parts))<1e-10)
  }
}
write <- function(x,n)write.csv(x,file.path(assets,n),row.names=FALSE,na='NA')
write(estimates,'gcorr-chromosome-set-estimates.csv')
write(wide,'gcorr-chromosome-set-summary.csv')
write(counts,'gcorr-chromosome-set-counts.csv')
write(process,'gcorr-chromosome-set-process.csv')
t <- read.csv(file.path(stage,'timing.csv'))
prep <- read.csv(file.path(stage,'set_preparation-phase.csv'))$seconds
fit <- read.csv(file.path(stage,'fit-phase.csv'))$seconds
write(data.frame(preparation_seconds=prep,fit_seconds=fit,
  analysis_seconds=t$analysis_seconds,process_wall_seconds=process$processWallSeconds,
  peak_working_set_GiB=process$peak_working_set_GiB),
  'gcorr-chromosome-set-timings.csv')
png(file.path(assets,'gcorr-chromosome-set-results.png'),width=1600,height=1300,res=150)
par(mfrow=c(2,1),mar=c(4,5,3,1))
range_h <- range(c(wide$h2_T2D-1.96*wide$se_h2_T2D,
  wide$h2_T2D+1.96*wide$se_h2_T2D,wide$h2_CAD-1.96*wide$se_h2_CAD,
  wide$h2_CAD+1.96*wide$se_h2_CAD,0),na.rm=TRUE)
plot(1:22-.12,wide$h2_T2D,pch=16,col='#247589',xlim=c(.5,22.5),ylim=range_h,
  xaxt='n',xlab='Chromosome marker set',ylab='Heritability (supplied n convention)',
  main='Joint chromosome-set LDSC: T2D and CAD')
axis(1,at=1:22,labels=1:22)
points(1:22+.12,wide$h2_CAD,pch=17,col='#9b582c')
arrows(1:22-.12,wide$h2_T2D-1.96*wide$se_h2_T2D,
  1:22-.12,wide$h2_T2D+1.96*wide$se_h2_T2D,angle=90,code=3,length=.025,col='#247589')
arrows(1:22+.12,wide$h2_CAD-1.96*wide$se_h2_CAD,
  1:22+.12,wide$h2_CAD+1.96*wide$se_h2_CAD,angle=90,code=3,length=.025,col='#9b582c')
abline(h=0,lty=2,col='grey55');legend('topright',c('T2D','CAD'),pch=c(16,17),
  col=c('#247589','#9b582c'),bty='n')
valid <- is.finite(wide$rg)&is.finite(wide$se_rg)
limits <- range(c(wide$rg[valid]-1.96*wide$se_rg[valid],
  wide$rg[valid]+1.96*wide$se_rg[valid],wide$rg,-1,1),na.rm=TRUE)
plot(1:22,wide$rg,pch=16,col='#247589',xaxt='n',ylim=limits,
  xlab='Chromosome marker set',ylab='Genetic correlation (raw ratio)',
  main='CAD/T2D chromosome correlations: cross intercept fixed to zero')
axis(1,at=1:22,labels=1:22)
arrows(which(valid),wide$rg[valid]-1.96*wide$se_rg[valid],which(valid),
  wide$rg[valid]+1.96*wide$se_rg[valid],angle=90,code=3,length=.025,col='#247589')
abline(h=c(-1,0,1),lty=2,col='grey65')
dev.off()
markdown <- function(x) {
  x[] <- lapply(x,function(v)if(is.numeric(v))
    ifelse(is.na(v),'unavailable',formatC(v,format='g',digits=4)) else
    ifelse(is.na(v),'',as.character(v)))
  paste(c(paste0('| ',paste(names(x),collapse=' | '),' |'),
    paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
    apply(x,1,function(r)paste0('| ',paste(r,collapse=' | '),' |'))),collapse='\n')
}
shown <- wide[,c('set','h2_T2D','se_h2_T2D','h2_CAD','se_h2_CAD','rg','se_rg')]
names(shown) <- c('Set','T2D h2','T2D SE','CAD h2','CAD SE','rg','rg SE')
text <- paste(markdown(shown),'\n\nJoint model totals:\n',
  markdown(totals[,c('quantity','trait1','trait2','estimate','standard_error')]),
  sprintf('\n\nSet-score preparation took %.2f seconds; fitting took %.2f seconds. The analysis clock was %.2f seconds, process wall %.2f seconds, and peak R working set %.2f GiB. These are single observed runs; the configured budget is four threads, but this fitter exposes no thread argument.',
    prep,fit,t$analysis_seconds,process$processWallSeconds,process$peak_working_set_GiB),
  sprintf('\n\n%d of 22 set correlations have unavailable point estimates, %d have unavailable standard errors, and %d finite raw ratios lie outside [-1,1]. Raw negative h2 and invalid ratios remain visible; intervals in the plot are estimate +/- 1.96 SE, without clipping or a calibrated local-significance claim.',
    sum(!is.finite(wide$rg)),sum(!is.finite(wide$se_rg)),sum(abs(wide$rg)>1,na.rm=TRUE)),
  '\n\n![Chromosome-set heritability and correlation](data/gcorr-chromosome-set-results.png)',sep='\n')
saveRDS(list(text=text,totals=totals,wide=wide),file.path(destination,'set-tutorial-section.rds'))
cat('SET_EXPORT_COMPLETED: set sums match joint model totals\n')
