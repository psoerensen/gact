args <- commandArgs(TRUE)
stage <- args[1]; threads <- as.integer(args[2]); out <- args[3]
writeLines(as.character(Sys.getpid()),file.path(out,"r-process.pid"))
.libPaths(c(normalizePath('build/task-packages/library',winslash='/'),.libPaths()))
library(gbase);library(gcorr)
stopifnot(gcorr_native_info()$gcorr$runtime_version=='0.1.1')
writeLines(capture.output(c(sessionInfo(),gcorr_native_info())),file.path(out,'session.txt'))
input_path <- Sys.getenv('GCORR_CAD_INPUT')
if(stage!='gcorr_inputs') {
  if(!nzchar(input_path)) input_path <- file.path(dirname(out),paste0('gcorr_inputs-',threads),'output.rds')
  input <- readRDS(input_path)
  # Old caches used total participants. An explicitly audited, aligned n cache
  # allows their genotype/LD metadata to be reused without copying it again.
  summary_n_path <- Sys.getenv('GCORR_SUMMARY_N_RDS')
  if(nzchar(summary_n_path))input$n <- readRDS(summary_n_path)$n
  if(is.null(input$n))stop('Input cache lacks extracted stat$n; rebuild it or supply SummaryNRds')
  stopifnot(identical(rownames(input$n),input$markers$rsids),
    identical(colnames(input$n),colnames(input$z)),all(is.finite(input$n)),all(input$n>0))
  input$sample_sizes <- setNames(vapply(seq_len(ncol(input$n)),function(j)
    if(all(input$n[,j]==input$n[1,j]))input$n[1,j] else NA_real_,numeric(1)),colnames(input$n))
  input$assumptions$scale <- 'Z/sqrt(extracted stat$n); database information-size convention'
  write.csv(data.frame(trait=colnames(input$n),markers=nrow(input$n),
    n_min=apply(input$n,2,min),n_max=apply(input$n,2,max),
    source='extracted stat$n',marker_selection=input$assumptions$qc),
    file.path(out,'analysis-input-audit.csv'),row.names=FALSE)
}
statistics <- function(input, index=seq_along(input$markers$rsids)) {
  setNames(lapply(seq_len(ncol(input$z)),function(j) {
    data.frame(rsids=input$markers$rsids[index],chr=input$markers$chr[index],
      pos=input$markers$pos[index],ea=input$markers$ea[index],nea=input$markers$nea[index],
      z=input$z[index,j],n=input$n[index,j],stringsAsFactors=FALSE)
  }),colnames(input$z))
}
timing <- function(label, expression) {
  began <- proc.time()[[3]];value <- force(expression)
  write.csv(data.frame(phase=label,seconds=unname(proc.time()[[3]]-began)),
    file.path(out,paste0(label,'-phase.csv')),row.names=FALSE)
  value
}
started <- Sys.time();clock <- proc.time()[[3]]
value <- switch(stage,
  gcorr_inputs={
    # Frozen caches from the existing gact website benchmark; no acquisition.
    old <- readRDS('build/benchmarks/manuscript-20261004-171426/qgg_glist-1/output.rds')
    fields <- c('bedfiles','bimfiles','famfiles','ids','n','nchr','rsids','chr','pos','a1','a2','af','maf','mchr')
    old <- old[fields]
    invisible(gc())
    ld <- readRDS('build/benchmarks/manuscript-20261005-142319/website_ldsc_inputs-1/output.rds')
    ld <- ld[c('z','n','ldscores','studies')]
    invisible(gc())
    original <- unlist(old$rsids,use.names=FALSE)
    index <- match(original,rownames(ld$z))
    maf <- unlist(old$maf,use.names=FALSE)
    ea <- unlist(old$a1,use.names=FALSE);nea <- unlist(old$a2,use.names=FALSE)
    chr <- as.character(unlist(old$chr,use.names=FALSE));pos <- unlist(old$pos,use.names=FALSE)
    # The website filters BED to GAlist$rsids, then uses all Glist markers.
    # This cached BED already implements that selection. Apply no extra MAF,
    # palindrome, indel, missingness, HWE or position-duplicate filter here.
    selected <- which(!is.na(index))
    marker <- data.frame(rsids=original[selected],chr=chr[selected],pos=pos[selected],
      ea=ea[selected],nea=nea[selected],stringsAsFactors=FALSE)
    write.csv(data.frame(check=c('reference_markers','common_database_score_rows',
      'matched_reference_summary_score_rows','retained_maf_le_0.05',
      'retained_palindromic','retained_duplicated_positions'),
      count=c(length(original),nrow(ld$z),nrow(marker),sum(maf[selected]<=.05),
        sum(paste0(marker$ea,marker$nea) %in% c('AT','TA','CG','GC')),
        sum(duplicated(paste(marker$chr,marker$pos))))),
      file.path(out,'website-filter-audit.csv'),row.names=FALSE)
    # Cached Z scores are oriented to the database marker table. Resolve that
    # orientation explicitly before comparing them with BIM allele-1 dosage.
    database_markers <- data.table::fread(
      'build/benchmarks/manuscript-20261004-143135/database-1/database/hsa.0.0.1/marker/markers.txt.gz',
      select=c('rsids','chr','pos','ea','nea'),data.table=FALSE)
    mapped <- match(marker$rsids,database_markers$rsids)
    stopifnot(!anyNA(mapped),all(as.character(database_markers$chr[mapped])==marker$chr),
      all(database_markers$pos[mapped]==marker$pos))
    db_ea <- database_markers$ea[mapped];db_nea <- database_markers$nea[mapped]
    direct <- db_ea==marker$ea & db_nea==marker$nea
    swapped <- db_ea==marker$nea & db_nea==marker$ea
    stopifnot(all(direct | swapped))
    z <- ld$z[index[selected],c('GWAS1','GWAS2'),drop=FALSE]
    n <- ld$n[index[selected],colnames(z),drop=FALSE]
    z <- z*ifelse(direct,1,-1)
    rm(database_markers,db_ea,db_nea)
    scores <- ld$ldscores[marker$rsids]
    stopifnot(identical(rownames(z),marker$rsids),identical(rownames(n),marker$rsids),
      all(is.finite(n)),all(n>0),all(is.finite(z)),all(is.finite(scores)),all(scores>0))
    # Cached qgg metadata is adapted to the public gbase Glist descriptor.
    # validate_Glist checks it against current BED/BIM/FAM before use.
    component <- rep.int(seq_along(old$rsids),lengths(old$rsids))[selected]
    g <- old[c('bedfiles','bimfiles','famfiles','ids','n','nchr')]
    for(field in c('rsids','chr','pos','a1','a2','af'))
      g[[field]] <- lapply(seq_along(old$rsids),function(i) {
        at <- which(component==i);source <- selected[at]-sum(lengths(old$rsids)[seq_len(i-1L)])
        x <- old[[field]][[i]][source];if(field=='chr')as.character(x)else x
      })
    for(field in c('bedfiles','bimfiles','famfiles'))
      g[[field]] <- normalizePath(g[[field]],winslash='/',mustWork=TRUE)
    g['study_ids'] <- list(NULL)
    g$physical_marker_counts <- as.double(lengths(old$rsids))
    g$source_marker_indices_0based <- lapply(seq_along(old$rsids),function(i)
      as.double(selected[component==i]-sum(lengths(old$rsids)[seq_len(i-1L)])-1L))
    g$sample_source_indices_0based <- as.double(seq_along(g$ids)-1L)
    g$mchr <- lengths(g$rsids)
    g$ldscores <- lapply(seq_along(g$rsids),function(i)scores[g$rsids[[i]]])
    class(g) <- c('gs_genotypes','Glist')
    rm(old,ld,original,index,maf,ea,nea,chr,pos,
      selected,component,mapped,direct,swapped)
    invisible(gc())
    gbase::validate_Glist(g)
    partition <- read.table('benchmarks/genomic_20k_50k/EUR-fourier-ac125e47.bed',header=TRUE)
    partition$chr <- sub('chr','',partition$chr,fixed=TRUE)
    labels <- rep(NA_character_,nrow(marker))
    for(c in unique(marker$chr)) {
      rows <- which(marker$chr==c);p <- partition[partition$chr==c,,drop=FALSE]
      at <- findInterval(marker$pos[rows],p$start)
      good <- at>0;good[good] <- marker$pos[rows][good]<p$stop[at[good]]
      labels[rows[good]] <- paste0('chr',c,'_',p$start[at[good]],'_',p$stop[at[good]])
    }
    keep <- !is.na(labels)
    # Boundary coverage and input provenance are explicit, not inferred later.
    regions <- split(marker$rsids[keep],labels[keep])
    block <- setNames(paste0('delete',pmin(200L,ceiling(seq_len(nrow(marker))*200/nrow(marker)))),marker$rsids)
    # Use the same extracted stat$n field as the existing qgg workflow.
    # Do not substitute study totals or infer a different effective-N formula.
    sample_sizes <- setNames(vapply(seq_len(ncol(n)),function(j)
      if(all(n[,j]==n[1,j]))n[1,j] else NA_real_,numeric(1)),colnames(n))
    write.csv(data.frame(trait=names(sample_sizes),sample_size=sample_sizes,
      source='extracted stat$n'),file.path(out,'sample-sizes.csv'),row.names=FALSE)
    write.csv(data.frame(qc=c('reference_markers','regression_markers','regional_markers','regions','unpartitioned_markers'),
      count=c(sum(g$physical_marker_counts),nrow(marker),sum(keep),length(regions),sum(!keep))),file.path(out,'input-counts.csv'),row.names=FALSE)
    list(markers=marker,z=z,n=n,sample_sizes=sample_sizes,ldscores=scores,blocks=block,
      Glist=g,regions=regions,reference_marker_count=nrow(marker),
      physical_reference_marker_count=sum(g$physical_marker_counts),
      assumptions=list(scale='Z/sqrt(extracted stat$n); database information-size convention, no liability conversion',
        score_source='stored EUR ordinary scores from the existing gact database; not recomputed',
        architecture='illustrative equal-marker baseline using legacy ordinary LD scores',
        genome_build='GRCh37: BIM rs429358 chr19:45411941 anchor verified',
        partition='Berisa/Pickrell EUR ac125e47bf7ff3e90be31f278a7b6a61daaba0dc',
        qc='website BED-to-GAlist marker selection; common reference/summary/LD-score rows; no extra QC exclusions'))
  },
  gcorr_ldsc={
    stat <- statistics(input)
    threshold <- max(.001*max(input$n),80)
    input$Glist <- NULL;input$regions <- NULL;input$z <- NULL;input$n <- NULL
    invisible(gc())
    fit <- timing('fit',gcorr_ldsc(stat,input$ldscores,input$blocks,
      marker_count=input$reference_marker_count,marginal_intercept='free',
      maximum_chi_square=threshold))
    write.csv(fit$estimates,file.path(out,'estimates.csv'),row.names=FALSE);fit
  },
  gcorr_gnova={
    ld_path <- Sys.getenv('GCORR_LDSC_ESTIMATES_CSV')
    if(!nzchar(ld_path))ld_path <- file.path(dirname(out),paste0('gcorr_ldsc-',threads),'estimates.csv')
    if(!file.exists(ld_path))stop('GNOVA needs a corrected LDSC estimates table from these inputs')
    ld_estimates <- read.csv(ld_path)
    stopifnot(setequal(ld_estimates$trait1[ld_estimates$quantity=='h2'],colnames(input$z)))
    intercept <- ld_estimates$estimate[ld_estimates$quantity=='cross_trait_intercept' &
      !is.na(ld_estimates$trait2)]
    if(length(intercept)!=1L || !is.finite(intercept))stop('No unique fitted cross intercept')
    write.csv(data.frame(source=normalizePath(ld_path,winslash='/'),
      md5=unname(tools::md5sum(ld_path)),cross_intercept=intercept),
      file.path(out,'ldsc-intercept-source.csv'),row.names=FALSE)
    stat <- statistics(input)
    input$Glist <- NULL;input$regions <- NULL;input$z <- NULL;input$n <- NULL
    invisible(gc())
    fit <- timing('fit',gcorr_gnova(stat,input$blocks,input$ldscores,
      overlap_intercepts=intercept))
    write.csv(fit$estimates,file.path(out,'estimates.csv'),row.names=FALSE);fit
  },
  gcorr_sumher={
    stat <- statistics(input);reference_n <- input$Glist$n
    architecture <- matrix(unname(input$ldscores),ncol=1L,
      dimnames=list(names(input$ldscores),'baseline'))
    input$Glist <- NULL;input$regions <- NULL;input$z <- NULL;input$n <- NULL
    invisible(gc())
    fit <- timing('fit',gcorr_sumher(stat,architecture,
      c(baseline=input$reference_marker_count),input$ldscores,input$blocks,
      reference_sample_size=reference_n,inflation='multiplicative'))
    write.csv(fit$estimates,file.path(out,'estimates.csv'),row.names=FALSE);fit
  },
  gcorr_sets={
    source('benchmarks/genomic_20k_50k/gcorr_set_scores.R')
    # Chromosomes enter only through this named marker-set list.
    chromosome <- sort(unique(as.integer(input$markers$chr)))
    sets <- setNames(lapply(chromosome,function(c)
      input$markers$rsids[input$markers$chr==as.character(c)]),paste0('chr',chromosome))
    # qgg's website set model normalizes by matched set membership.
    counts <- lengths(sets)
    stopifnot(sum(counts)==input$reference_marker_count)
    reference_n <- input$Glist$n
    input$Glist <- NULL;input$regions <- NULL;invisible(gc())
    prepared_sets <- timing('set_preparation',prepare_gcorr_set_scores(
      input$markers$rsids,input$ldscores,sets,counts))
    write.csv(data.frame(set=names(sets),regression_markers=lengths(sets),
      reference_markers=counts),file.path(out,'set-counts.csv'),row.names=FALSE)
    stat <- statistics(input)
    input$z <- NULL;input$n <- NULL;invisible(gc())
    fit <- timing('fit',gcorr_ldsc_partitioned(stat,
      prepared_sets$annotation_ld_scores,input$ldscores,prepared_sets$annotation_sums,
      input$blocks,annotation_overlap=prepared_sets$annotation_overlap,
      reference_sample_size=reference_n,
      reference_marker_count=prepared_sets$reference_marker_count,
      corrected=FALSE,marginal_intercept='free',overlap_intercept=0))
    write.csv(fit$estimates,file.path(out,'estimates.csv'),row.names=FALSE)
    # Keep compact estimates/metadata, without serializing the score matrix.
    fit$diagnostics$used_marker_ids <- NULL
    list(fit=fit,set_counts=data.frame(set=names(sets),regression_markers=lengths(sets),
      reference_markers=counts),assumptions=list(shared_samples=0,
      overlap_source='user understanding: CAD and T2D have no shared samples',
      cross_intercept=0,marginal_intercepts='estimated jointly across sets',
      ld_scores='legacy stored ordinary scores; masked by chromosome sets, no recomputation',
      normalization='matched marker counts by chromosome, as in the website set model',
      no_cross_chromosome_ld=TRUE))
  },
  gcorr_chr22=,gcorr_genome={
    regions <- if(stage=='gcorr_chr22')input$regions[grepl('^chr22_',names(input$regions))]else input$regions
    selected <- unlist(regions,use.names=FALSE)
    index <- match(selected,input$markers$rsids)
    maximum_components <- as.integer(Sys.getenv('GCORR_HESS_MAXIMUM_COMPONENTS','50'))
    stopifnot(length(maximum_components)==1L,is.finite(maximum_components),maximum_components>0L)
    rank_bound <- sum(pmin(lengths(regions),maximum_components,input$Glist$n-1L))
    # Admit the workload before doing any eigendecompositions. This conservative
    # bound depends on reference/region sizes and n, never on trait signal.
    if(rank_bound>=min(input$sample_sizes))
      stop('HESS rank bound exceeds supplied n; choose an explicit smaller HessMaximumComponents')
    write.csv(data.frame(trait=names(input$sample_sizes),sample_size=input$sample_sizes,
      regions=length(regions),markers=length(selected),maximum_components=maximum_components,
      minimum_eigenvalue=1,rank_upper_bound=rank_bound,
      minimum_correction_denominator=1-rank_bound/input$sample_sizes),
      file.path(out,'projection-configuration.csv'),row.names=FALSE)
    stat <- statistics(input,index)
    projection <- timing('preparation',gcorr_prepare_hess(stat,
      input$Glist,regions,nthreads=threads,maximum_components=maximum_components))
    rm(input,stat,index,selected);invisible(gc())
    retained_rank <- sum(vapply(projection$regions,`[[`,numeric(1),'retained_components'))
    write.csv(data.frame(trait=names(projection$sample_sizes),sample_size=projection$sample_sizes,
      retained_rank=retained_rank,correction_denominator=1-retained_rank/projection$sample_sizes),
      file.path(out,'rank-admission.csv'),row.names=FALSE)
    overhead <- unlist(projection$timings_seconds)
    write.csv(data.frame(phase=names(overhead),seconds=unname(overhead)),
      file.path(out,'preparation-overhead.csv'),row.names=FALSE)
    write.csv(do.call(rbind,Map(function(x,name) data.frame(region=name,
      markers=length(x$marker_ids),solver_dimension=x$solver_dimension,
      retained_components=x$retained_components,matrix_space=x$matrix_space,
      numerical_rank=x$numerical_rank,
      retained_eigenvalue_fraction=x$retained_positive_eigenvalue_fraction,
      workspace_bytes=x$workspace_bytes,t(x$timings_seconds),check.names=FALSE),
      projection$regions,names(projection$regions))),file.path(out,'regional-phases.csv'),row.names=FALSE)
    fits <- timing('marginal_fit',setNames(lapply(projection$traits,function(t)
      gcorr_fit_hess(projection,t)),projection$traits))
    write.csv(do.call(rbind,lapply(fits,`[[`,'estimates')),file.path(out,'marginal-estimates.csv'),row.names=FALSE)
    pair <- timing('rho_fit',gcorr_fit_hess(projection,
      shared_sample_size=0,phenotypic_correlation=0))
    write.csv(pair$estimates,file.path(out,'estimates.csv'),row.names=FALSE)
    # Retain the user-specified no-overlap assumption. Historical physical-N
    # sensitivity counts are not interchangeable with effective information N.
    scenarios <- data.frame(scenario='zero_overlap',
      shared_sample_size=0,phenotypic_correlation=0)
    sensitivity <- timing('sensitivity_fit',do.call(rbind,lapply(seq_len(nrow(scenarios)),function(i) {
      fit <- if(i==1L)pair else gcorr_fit_hess(projection,
        shared_sample_size=scenarios$shared_sample_size[i],
        phenotypic_correlation=scenarios$phenotypic_correlation[i])
      cbind(scenarios[rep(i,sum(fit$estimates$component=='total')),,drop=FALSE],
        fit$estimates[fit$estimates$component=='total',,drop=FALSE])
    })))
    write.csv(sensitivity,file.path(out,'sensitivity-estimates.csv'),row.names=FALSE)
    list(projection=projection,marginal_fits=fits,pair_fit=pair,sensitivity=sensitivity)
  },stop('Unknown stage'))
seconds <- unname(proc.time()[[3]]-clock)
serialization <- unname(system.time(saveRDS(value,file.path(out,'output.rds'),compress=FALSE))[[3]])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),
  analysis_seconds=seconds,serialization_seconds=serialization,status='completed'),
  file.path(out,'timing.csv'),row.names=FALSE)
writeLines(capture.output(str(value,max.level=1)),file.path(out,'output-structure.txt'))
