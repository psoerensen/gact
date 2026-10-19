# Fresh timing of the existing CAD/T2D tutorial, with a pre-computation baseline.
args <- commandArgs(TRUE)
stage <- args[[1]]; out <- args[[2]]; run <- dirname(out)
Sys.setenv(LC_ALL="C", LANG="C", OMP_NUM_THREADS="4", OMP_DYNAMIC="FALSE")
.libPaths(c(normalizePath("build/task-packages/library"), .libPaths()))
library(gcorr)
stopifnot(packageVersion("gcorr") == "0.1.6",
  gcorr_native_info()$gcorr$runtime_version == "0.1.4")
writeLines(as.character(Sys.getpid()), file.path(out,"r-process.pid"))
writeLines(capture.output(c(sessionInfo(),gcorr_native_info())), file.path(out,"session.txt"))
timed <- function(label, expression) {
  invisible(gc(reset=TRUE))
  # Let the monitor record the resident/private baseline before starting.
  writeLines(label, file.path(out,"ready"))
  deadline <- Sys.time()+30
  while (!file.exists(file.path(out,"go"))) {
    if(Sys.time()>deadline) stop("Memory monitor did not acknowledge baseline")
    Sys.sleep(.05)
  }
  started <- proc.time()[[3]]
  value <- force(expression)
  seconds <- proc.time()[[3]]-started
  write.csv(data.frame(stage=stage,phase=label,seconds=seconds),
    file.path(out,"timing.csv"), row.names=FALSE)
  writeLines("done",file.path(out,"done"))
  write.csv(gc(),file.path(out,"r-heap-memory.csv"))
  value
}
if (stage %in% c("hess_raw","hess_partial")) {
  projection <- readRDS(file.path(run,"hess_prepare","projection.rds"))
  fit <- timed("fit", if (stage=="hess_raw")
    gcorr_fit_hess(projection,shared_sample_size=0,phenotypic_correlation=0) else
    gcorr_fit_hess(projection,shared_sample_size=0,phenotypic_correlation=0,
      genomic_control=c(GWAS1=1,GWAS2=1.18),
      genomic_control_provenance=c(GWAS1="Unchanged supplied Z scale",
        GWAS2="Partial final CAD meta-analysis correction sensitivity: 1.18")))
  write.csv(fit$estimates,file.path(out,"estimates.csv"),row.names=FALSE)
} else {
  input <- readRDS("build/benchmarks/manuscript-20261007-003431/gcorr_inputs-4/output.rds")
  stopifnot(nrow(input$markers)==7221923L, identical(rownames(input$n),input$markers$rsids),
    identical(colnames(input$n),c("GWAS1","GWAS2")),
    all(abs(input$n[,1]-49071.2732949614)<1e-8),
    all(abs(input$n[,2]-40743.1524049809)<1e-8))
  stat <- setNames(lapply(seq_len(ncol(input$z)),function(j)
    data.frame(rsids=input$markers$rsids,chr=input$markers$chr,pos=input$markers$pos,
      ea=input$markers$ea,nea=input$markers$nea,z=input$z[,j],n=input$n[,j],
      stringsAsFactors=FALSE)),colnames(input$z))
  reference_n <- input$Glist$n
  if(stage=="hess_prepare") {
    Glist <- input$Glist; regions <- input$regions
    rm(input)
    projection <- timed("preparation",gcorr_prepare_hess(stat,Glist,regions,
      nthreads=4,maximum_components="auto"))
    stopifnot(length(projection$regions)==1702L,
      sum(vapply(projection$regions,`[[`,numeric(1),"retained_components"))==18696)
    saveRDS(projection,file.path(out,"projection.rds"))
    write.csv(data.frame(phase=names(projection$timings_seconds),
      seconds=unlist(projection$timings_seconds)),
      file.path(out,"preparation-phases.csv"),row.names=FALSE)
    write.csv(do.call(rbind,Map(function(x,name) data.frame(region=name,
      markers=x$marker_count,retained_components=x$retained_components,
      t(x$timings_seconds),check.names=FALSE),projection$regions,names(projection$regions))),
      file.path(out,"regional-phases.csv"),row.names=FALSE)
  } else {
    sets <- if(stage=="sets") {
      chr <- sort(unique(as.integer(input$markers$chr)))
      setNames(lapply(chr,function(c) input$markers$rsids[input$markers$chr==as.character(c)]),paste0("chr",chr))
    } else NULL
    threshold <- max(.001*max(input$n),80)
    input$Glist <- input$regions <- input$z <- input$n <- NULL
    architecture <- if(stage=="sumher") matrix(unname(input$ldscores),ncol=1,
      dimnames=list(names(input$ldscores),"baseline")) else NULL
    cross_intercept <- if(stage=="gnova") {
      estimates <- read.csv(file.path(run,"ldsc","estimates.csv"))
      estimates$estimate[estimates$quantity=="cross_trait_intercept" & !is.na(estimates$trait2)]
    } else NULL
    fit <- timed("fit",switch(stage,
      ldsc=gcorr_ldsc(stat,input$ldscores,input$blocks,
        marker_count=input$reference_marker_count,marginal_intercept="free",
        maximum_chi_square=threshold),
      gnova=gcorr_gnova(stat,input$blocks,input$ldscores,overlap_intercepts=cross_intercept),
      sumher=gcorr_sumher(stat,architecture,c(baseline=input$reference_marker_count),
        input$ldscores,input$blocks,reference_sample_size=reference_n,inflation="multiplicative"),
      sets={
        # Match names once, then retain integer membership throughout fitting.
        scores <- gcorr_set_scores(input$markers$rsids,input$ldscores,sets)
        counts <- lengths(sets)
        overlap <- diag(as.numeric(counts)); dimnames(overlap) <- list(names(sets),names(sets))
        gcorr_ldsc_partitioned(stat,scores,input$ldscores,counts,input$blocks,
          annotation_overlap=overlap,reference_sample_size=reference_n,
          reference_marker_count=sum(counts),marginal_intercept="free",overlap_intercept=0)
      },stop("Unknown timing stage")))
    write.csv(fit$estimates,file.path(out,"estimates.csv"),row.names=FALSE)
  }
}
writeLines("completed",file.path(out,"success"))
