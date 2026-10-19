# Bounded compatibility check of saved compact moments. No new projection,
# acquisition, LD/eigenvector storage, or change to the original timing study.
# Usage from gsuite: Rscript --vanilla verify_hess_robust.R GENOME_RDS OUT_DIR
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
Sys.setenv(LC_ALL="C",LANG="C")
.libPaths(c(normalizePath("build/task-packages/library",winslash="/"),.libPaths()))
library(gcorr)
stopifnot(as.character(packageVersion("gcorr"))=="0.1.2")
dir.create(args[2L],recursive=TRUE,showWarnings=FALSE)
study <- readRDS(args[1L])
projection <- study$projection
reference <- study$pair_fit$estimates
rm(study);invisible(gc())
traits <- projection$traits
stopifnot(identical(traits,c("GWAS1","GWAS2")))
warnings <- character()
capture <- function(expr) withCallingHandlers(expr,warning=function(w) {
  warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
})
unchanged <- capture(gcorr_fit_hess(projection,traits,
  shared_sample_size=0,phenotypic_correlation=0))
stopifnot(isTRUE(all.equal(unchanged$estimates,reference,tolerance=1e-12)))
write.csv(unchanged$diagnostics$conditioning,file.path(args[2L],"saved-cap20-conditioning.csv"),row.names=FALSE)
write.csv(unchanged$diagnostics$inference$traits,file.path(args[2L],"saved-cap20-inference.csv"),row.names=FALSE)
negative <- unchanged$diagnostics$inference$negative_regions
write.csv(data.frame(trait=rep(names(negative),lengths(negative)),
  region=unlist(negative,use.names=FALSE)),file.path(args[2L],"saved-cap20-negative-regions.csv"),row.names=FALSE)
cat("OLD_SCHEMA_ESTIMATES_UNCHANGED\n")
print(unchanged$diagnostics$conditioning)
print(unchanged$diagnostics$inference$traits)
policy <- getFromNamespace(".gcorr_hess_rank_policy","gcorr")("auto",
  projection$sample_sizes,lengths(lapply(projection$regions,`[[`,"marker_ids")))
write.csv(data.frame(cap=policy$maximum_components,regions=policy$region_count,
  conservative_total_rank=policy$conservative_rank_bound,
  smallest_sample_size=policy$smallest_sample_size,
  conservative_denominator=1-policy$conservative_rank_bound/policy$smallest_sample_size),
  file.path(args[2L],"automatic-genome-rank-preflight.csv"),row.names=FALSE)
stopifnot(policy$maximum_components==11,policy$conservative_rank_bound<=min(projection$sample_sizes)/2)
cat("AUTO_GENOME_PREFLIGHT_VERIFIED\n")
# Explicitly labelled sensitivity: not a reconstruction of all applied GC.
adjusted <- capture(gcorr_fit_hess(projection,"GWAS2",genomic_control=c(GWAS2=1.18),
  genomic_control_provenance=c(GWAS2="Sensitivity only: CAD meta lambda 1.18; complete cohort correction history unverified")))
k <- sum(vapply(projection$regions,`[[`,numeric(1),"retained_components"))
n <- projection$sample_sizes["GWAS2"]
q <- sum(vapply(projection$regions,function(x)x$quadratic_forms["GWAS2","GWAS2"],numeric(1)))
expected <- (1.18*q-k/n)/(1-k/n)
total <- adjusted$estimates[adjusted$estimates$component=="total",]
stopifnot(abs(total$estimate-expected)<1e-12)
write.csv(total,file.path(args[2L],"saved-cap20-cad-1.18-sensitivity.csv"),row.names=FALSE)
write.csv(data.frame(trait=names(adjusted$input_scaling$z_squared_multiplier),
  z_squared_multiplier=unname(adjusted$input_scaling$z_squared_multiplier),
  provenance=unname(adjusted$input_scaling$provenance)),
  file.path(args[2L],"scaling-provenance.csv"),row.names=FALSE)
writeLines(unique(warnings),file.path(args[2L],"expected-rank-warnings.txt"))
cat("EXPLICIT_CAD_SENSITIVITY_VERIFIED\n")
print(total)
cat("ROBUST_HESS_SAVED_MOMENT_CHECK_PASSED\n")
