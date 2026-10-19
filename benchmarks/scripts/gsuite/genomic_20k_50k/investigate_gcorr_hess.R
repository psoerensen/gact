# Focused real-region numerical and model diagnostics. No acquisition, fitting
# campaign, persistent LD/eigenvectors, or alteration of production estimators.
# Usage: Rscript --vanilla investigate_gcorr_hess.R INPUT_RDS GENOME_RDS OUT_DIR [CAD_DATABASE_GZ]
# INPUT_RDS can be the original corrected gcorr input or the private reduced
# descriptor produced by this script. Only aggregate diagnostics are exported.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(3L,4L))
Sys.setenv(LC_ALL='C',LANG='C')
out <- args[3L]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE)

genome <- readRDS(args[2L])
p <- genome$projection
ns <- p$sample_sizes
rank <- vapply(p$regions,`[[`,numeric(1),'retained_components')
q <- t(vapply(p$regions,function(x)c(x$quadratic_forms[1,1],x$quadratic_forms[2,2],x$quadratic_forms[1,2]),numeric(3)))
Q <- colSums(q)[1:2L]*ns
K <- sum(rank)
totals <- data.frame(trait=names(ns),n=unname(ns),retained_rank=K,
  projected_energy=unname(Q),unit_noise_subtraction=K,
  energy_minus_noise=unname(Q-K),correction_denominator=unname(1-K/ns),
  h2=unname((Q-K)/(ns-K)),zero_h2_reinflation_factor=unname(K/Q))
native_totals <- genome$pair_fit$estimates
native_totals <- native_totals[native_totals$component=='total' & native_totals$quantity=='h2',]
stopifnot(max(abs(totals$h2-native_totals$estimate[match(totals$trait,native_totals$trait1)]))<1e-10)
write_table(totals,'genome-moment-decomposition')
cat('GENOME_DECOMPOSITION_VERIFIED\n');print(totals)
write_table(data.frame(trait=names(ns),n=unname(ns),regions=length(rank),
  current_cap=p$options$maximum_components,
  upstream_default_gwse_threshold=2,
  upstream_automatic_cap=pmin(50,ceiling(ns/(2*length(rank))))),
  'upstream-rank-policy-comparison')
# Matrix convention only, holding the measured moments/ranks fixed. This is
# not an execution of upstream HESS with its separate automatic rank policy.
rank_conventions <- do.call(rbind,lapply(seq_along(ns),function(j) {
  a <- rank/ns[j];b <- q[,j]-a
  upstream_total <- sum(b)+length(rank)*sum(a*b)/(1-sum(a))
  data.frame(trait=names(ns)[j],native_row_rank_total=totals$h2[j],
    pinned_upstream_column_rank_total=upstream_total,
    difference=upstream_total-totals$h2[j])
}))
write_table(rank_conventions,'joint-rank-convention-sensitivity')

# N scaling changes the denominator and the effect scale, but cannot change
# the sign of Q-K while N>K. This is a diagnostic, not a new N convention.
nscale <- do.call(rbind,lapply(c(1,2,4),function(scale)data.frame(
  trait=names(ns),diagnostic_n_multiplier=scale,n=unname(ns*scale),
  h2=unname((Q-K)/(ns*scale-K)))))
write_table(nscale,'sample-size-sensitivity')

# Explicit Z reinflation sensitivity: Q is multiplied by lambda, whereas
# cross moments would be multiplied by sqrt(lambda1*lambda2). These values
# are not corrections inferred or silently applied to the study statistics.
ldsc <- read.csv('benchmarks/genomic_20k_50k/manuscript-results/gcorr-cad-t2d-20261007/gcorr-corrected-ldsc-estimates.csv')
intercepts <- ldsc[ldsc$quantity=='intercept',]
inflation <- do.call(rbind,lapply(seq_along(ns),function(j) {
  fitted_intercept <- intercepts$estimate[match(names(ns)[j],intercepts$trait1)]
  factors <- c(1,1.05,1.10,1.18,1/fitted_intercept)
  data.frame(trait=names(ns)[j],scenario=c('unchanged','sensitivity_1.05',
    'sensitivity_1.10','sensitivity_1.18',
    'inverse_LDSC_intercept_sensitivity'),lambda=factors,
    h2=(factors*Q[j]-K)/(ns[j]-K),ldsc_intercept=fitted_intercept)
}))
write_table(inflation,'reinflation-sensitivity')

# Explain the existing strict conditional-uncertainty gate, without relaxing
# it or presenting raw plug-in variances as qualified standard errors.
h <- genome$marginal_fits
gates <- do.call(rbind,lapply(seq_along(ns),function(j) {
  tab <- h[[names(ns)[j]]]$estimates
  point <- tab$estimate[tab$component=='region' & tab$quantity=='h2']
  total <- totals$h2[j]
  v <- (1-total)/ns[j]
  raw_variance <- 2*rank*v^2+4*v*point
  data.frame(trait=names(ns)[j],regions=length(point),negative_local_h2=sum(point<0),
    total_in_unit_interval=total>=0 && total<=1,
    raw_plugin_variances_negative=sum(raw_variance<0),
    raw_plugin_variances_nonnegative=sum(raw_variance>=0),
    published_standard_errors_available=sum(!is.na(tab$standard_error)),
    policy='joint conditional gate requires all regional h2 nonnegative')
}))
write_table(gates,'uncertainty-gates')

# Choose by marker count, never by trait effects or the resulting estimates.
sizes <- vapply(p$regions,function(x)length(x$marker_ids),integer(1))
at <- unique(c(which.min(sizes),which.min(abs(sizes-300L)),which.min(abs(sizes-3000L))))
saved <- p$regions[at]
regions <- lapply(saved,`[[`,'marker_ids')
rm(genome,p,q,h);invisible(gc())
input <- readRDS(args[1L])
if (!is.null(input$direct_descriptor)) {
  descriptor <- input$direct_descriptor
} else {
  selected <- unlist(regions,use.names=FALSE)
  idx <- match(selected,input$markers$rsids)
  stopifnot(!anyNA(idx))
  g <- input$Glist
  all_ids <- unlist(g$rsids,use.names=FALSE)
  panel_at <- match(selected,all_ids)
  component <- rep.int(seq_along(g$rsids),lengths(g$rsids))[panel_at]
  source <- unlist(g$source_marker_indices_0based,use.names=FALSE)[panel_at]
  a1 <- unlist(g$a1,use.names=FALSE)[panel_at]
  a2 <- unlist(g$a2,use.names=FALSE)[panel_at]
  marker <- input$markers[idx,]
  direct <- marker$ea==a1 & marker$nea==a2
  swap <- marker$ea==a2 & marker$nea==a1
  stopifnot(all(direct | swap))
  descriptor <- list(markers=marker,z=input$z[idx,,drop=FALSE]*ifelse(direct,1,-1),
    n=ns,component=component,source_indices_0based=source,
    bedfiles=g$bedfiles,famfiles=g$famfiles,ids=g$ids,
    sample_indices_0based=g$sample_source_indices_0based,
    n_reference=g$n,regions=regions)
  saveRDS(list(direct_descriptor=descriptor),file.path(out,'private-direct-input.rds'),compress=FALSE)
}
rm(input);invisible(gc())
stopifnot(identical(descriptor$n,ns))

# Independent PLINK SNP-major two-bit decoder for this bounded oracle only.
# Production genotype operations continue to use installed gmat/gbits.
read_region <- function(ids) {
  rows <- match(ids,descriptor$markers$rsids)
  components <- unique(descriptor$component[rows])
  stopifnot(length(components)==1L,!anyNA(rows))
  component <- components[1L]
  fam <- read.table(descriptor$famfiles[component],stringsAsFactors=FALSE)
  stopifnot(identical(as.character(fam[[2L]]),descriptor$ids))
  n <- nrow(fam);stride <- ceiling(n/4)
  con <- file(descriptor$bedfiles[component],'rb');on.exit(close(con))
  stopifnot(identical(as.integer(readBin(con,'raw',3L)),c(108L,27L,1L)))
  dosage <- matrix(0,n,length(rows))
  for(j in seq_along(rows)) {
    seek(con,where=3+descriptor$source_indices_0based[rows[j]]*stride,origin='start')
    bytes <- as.integer(readBin(con,'raw',stride))
    stopifnot(length(bytes)==stride)
    codes <- as.vector(rbind(bitwAnd(bytes,3L),bitwAnd(bitwShiftR(bytes,2L),3L),
      bitwAnd(bitwShiftR(bytes,4L),3L),bitwAnd(bitwShiftR(bytes,6L),3L)))[seq_len(n)]
    dosage[,j] <- c(2,NA_real_,1,0)[codes+1L]
  }
  dosage <- dosage[descriptor$sample_indices_0based+1L,,drop=FALSE]
  stopifnot(nrow(dosage)==503L)
  list(dosage=dosage,z=descriptor$z[rows,,drop=FALSE])
}
standardize <- function(g) {
  g <- sweep(g,2,colMeans(g,na.rm=TRUE),'-');g[is.na(g)] <- 0
  norm <- sqrt(colSums(g*g));stopifnot(all(norm>0))
  sweep(g,2,norm,'/')
}

oracle <- caps <- split_noise <- list()
set.seed(20261007)
split_rows <- sample.int(503L)
for(name in names(regions)) {
  cat('DIRECT_REGION_BEGIN',name,'\n');flush.console()
  raw <- read_region(regions[[name]])
  w <- standardize(raw$dosage)
  began <- proc.time()[[3L]]
  decomposition <- svd(w,nu=0L,nv=min(dim(w)))
  lambda <- decomposition$d^2
  projection <- crossprod(decomposition$v,raw$z)
  eligible <- which(lambda>max(1,lambda[1L]*1e-10))
  make_q <- function(k) {
    keep <- head(eligible,k)
    scaled <- sweep(projection[keep,,drop=FALSE],1,sqrt(lambda[keep]),'/')
    scaled <- sweep(scaled,2,sqrt(ns),'/')
    crossprod(scaled)
  }
  exact <- make_q(20L)
  stored <- saved[[name]]$quadratic_forms
  error <- max(abs(exact-stored))
  stopifnot(error<1e-9, length(head(eligible,20L))==saved[[name]]$retained_components)
  # Marker-space dense eigendecomposition is affordable for the two small
  # regions and is independent of the production sample-space route.
  dense_error <- NA_real_
  if(ncol(w)<=503L) {
    eig <- eigen(crossprod(w),symmetric=TRUE)
    keep <- head(which(eig$values>max(1,eig$values[1L]*1e-10)),20L)
    dp <- crossprod(eig$vectors[,keep,drop=FALSE],raw$z)
    dp <- sweep(sweep(dp,1,sqrt(eig$values[keep]),'/'),2,sqrt(ns),'/')
    dense_error <- max(abs(crossprod(dp)-stored))
    stopifnot(dense_error<1e-9)
  }
  oracle[[name]] <- data.frame(region=name,markers=ncol(w),reference_samples=nrow(w),
    matrix_space=saved[[name]]$matrix_space,retained_components=saved[[name]]$retained_components,
    max_absolute_moment_error=error,dense_marker_space_error=dense_error,
    direct_seconds=unname(proc.time()[[3L]]-began))
  caps[[name]] <- do.call(rbind,lapply(c(5L,10L,12L,15L,20L,50L),function(cap) {
    keep <- head(eligible,cap);qq <- make_q(cap);k <- length(keep)
    data.frame(region=name,maximum_components=cap,retained_components=k,
      reference_trace_fraction=sum(lambda[keep])/sum(lambda),trait=names(ns),
      projected_energy=diag(qq)*ns,energy_per_component=diag(qq)*ns/k,
      h2_zero_outside=(diag(qq)*ns-k)/(ns-k))
  }))
  # One fixed disjoint reference split illustrates LD sampling error. It is
  # not a calibration campaign or an estimate of the GWAS noise covariance.
  train_g <- raw$dosage[split_rows[1:251],,drop=FALSE]
  holdout_g <- raw$dosage[split_rows[252:503],,drop=FALSE]
  variable <- apply(train_g,2,var,na.rm=TRUE)>0 &
    apply(holdout_g,2,var,na.rm=TRUE)>0
  stopifnot(all(!is.na(variable)),sum(variable)>1L)
  # This exclusion applies only to the illustrative reference split, whose
  # correlation is undefined for monomorphic columns. Full-study marker
  # selection and all projection oracles above remain unchanged.
  fit_w <- standardize(train_g[,variable,drop=FALSE])
  check_w <- standardize(holdout_g[,variable,drop=FALSE])
  se <- svd(fit_w,nu=0L,nv=min(dim(fit_w)))
  slambda <- se$d^2
  keep <- head(which(slambda>max(1,slambda[1L]*1e-10)),20L)
  directions <- se$v[,keep,drop=FALSE]
  # tr(U' R_holdout U Lambda_training^-1)/K: equals one only when
  # the held-out and fitted reference covariance agree in retained directions.
  ratio <- colSums((check_w %*% directions)^2)/slambda[keep]
  split_noise[[name]] <- data.frame(region=name,training_samples=251,
    holdout_samples=252,retained_components=length(keep),
    full_region_markers=ncol(raw$dosage),split_variable_markers=sum(variable),
    holdout_noise_energy_per_component=mean(ratio),
    interpretation='single reference-only split; not a GWAS noise correction')
  cat('DIRECT_REGION_VERIFIED',name,error,'\n');flush.console()
  write_table(do.call(rbind,oracle),'regional-projection-oracle')
  write_table(do.call(rbind,caps),'regional-truncation-sensitivity')
  write_table(do.call(rbind,split_noise),'reference-split-noise-diagnostic')
}
write_table(do.call(rbind,oracle),'regional-projection-oracle')
write_table(do.call(rbind,caps),'regional-truncation-sensitivity')
write_table(do.call(rbind,split_noise),'reference-split-noise-diagnostic')

if(length(args)==4L) {
  db <- data.table::fread(args[4L],select=c('rsids','b','seb','p'),data.table=FALSE)
  db <- db[match(descriptor$markers$rsids,db$rsids),]
  stopifnot(!anyNA(db$rsids))
  ratio_z <- db$b/db$seb
  p_z <- qnorm(db$p/2,lower.tail=FALSE)
  finite <- is.finite(p_z) & p_z>0 & abs(ratio_z)>.1
  diagnostic <- (ratio_z[finite]/p_z[finite])^2
  write_table(data.frame(markers=nrow(db),
    max_absolute_cached_vs_database_z_error=max(abs(abs(descriptor$z[,2L])-abs(ratio_z))),
    p_vs_effect_se_markers=sum(finite),median_effect_se_squared_over_p_z_squared=median(diagnostic),
    min_ratio=min(diagnostic),max_ratio=max(diagnostic)), 'cad-z-and-p-provenance')
}
cat('HESS_FOCUSED_INVESTIGATION_PASSED\n')
