#' PoPS-style chromosome-held-out gene ranking
#'
#' Fits gene features jointly by ridge or Bayesian regression and scores
#' every gene from a model that excluded its chromosome. Gene evidence can come
#' from MAGMA-style preparation, `glma::vegas()`, or another checked source, but
#' must already be represented by one finite score per gene.
#'
#' @param stat Data frame with unique `gene`, `chromosome`, and finite `score`
#'   columns. Chromosome labels define the held-out folds.
#' @param features Finite numeric gene-by-feature matrix with row names exactly
#'   matching `stat$gene` and unique feature names.
#' @param covariates Optional finite numeric matrix with the same gene row order,
#'   or data frame containing `gene` followed by numeric nuisance covariates.
#' @param sampling Explicit gene working covariance, using the same
#'   `list(gene_covariance=...)` contract as [gsea()]. Blocks must partition the
#'   genes and may not span chromosomes. Omitted sampling uses independent genes
#'   with sets or Bayesian methods.
#' @param control_features Feature names retained in every fold regardless of
#'   their marginal p-values.
#' @param control Marginal threshold, ridge penalty grid or Bayesian sampler
#'   settings and checked execution limits. Ridge defaults select at `p < 0.05`,
#'   retain at most 256 features per fold,
#'   allow a 1,024-feature union and choose among 25 penalties from `1e-2` to
#'   `1e10` by generalized cross-validation.
#' @param sets Familiar alternative to a feature matrix: named gene sets.
#' @param method "rr", "bayesC" or "bayesR"; all predictions exclude the target chromosome.
#' @param threshold Familiar alias for control$marginal_p_threshold.
#' @param pi,nit,nburn,nthin,nrun,seed Bayesian active probability or mixture weights,
#'   sampling sweeps, burn-in, thinning, chains and starting seed.
#' @param updateB,updateE,updatePi Learn effect variance, residual variance and mixture weights.
#' @param prior Advanced native Bayesian prior settings on prepared training coordinates.
#' @return A serializable `pops_fit` with ranked gene `predictions`, fold
#'   diagnostics, selected feature coefficients, transformations and provenance.
#'   Bayesian fits include separate preparation, BLR and prediction timings per
#'   fold; upstream caller-supplied gene-statistic time is unmeasured (NA).
#' @details Feature selection, nuisance adjustment and fitting are repeated
#'   after excluding the target chromosome. Feature centering/scaling uses only
#'   the declared feature matrix and no association response. Scores are
#'   feature-based priorities, not association statistics, p-values or causal
#'   probabilities. Locus construction and choosing the top gene within a locus
#'   are separate downstream operations.
#' @export
pops <- function(stat,features=NULL,covariates=NULL,sampling=NULL,
                 control_features=character(),control=list(),sets=NULL,
                 method="rr",threshold=NULL,pi=NULL,nit=5000L,nburn=1000L,
                 nthin=1L,nrun=2L,seed=10L,updateB=TRUE,updateE=TRUE,
                 updatePi=TRUE,prior=list()){
  workflow_start <- proc.time()[["elapsed"]]
  call <- match.call()
  method <- match.arg(tolower(method),c("rr","bayesc","bayesr"))
  if (!is.null(threshold)) {
    if (!is.null(control$marginal_p_threshold))
      .gsuite_stop("Specify either threshold or control$marginal_p_threshold, not both.")
    control$marginal_p_threshold <- threshold
  }
  if (!is.null(sets)) {
    if (!is.null(features)) .gsuite_stop("Specify either features or sets, not both.")
    sets <- gbase::mapSets(lapply(sets, unique), rsids = stat$gene, index=FALSE)
    sets <- sets[lengths(sets) > 0L]
    if (is.null(sampling)) sampling <- list(gene_covariance="independent")
    # Keep memberships throughout fitting and prediction for the simple exact
    # model. Screening, nuisance covariates and correlated genes retain the
    # existing general feature-design route.
    direct_sets <- method!="rr" && is.null(covariates) &&
      identical(sampling$gene_covariance,"independent") &&
      is.null(control$marginal_p_threshold)
    if(direct_sets) features <- sets else {
      features <- matrix(0, length(stat$gene), length(sets),
        dimnames = list(stat$gene, names(sets)))
      for (j in seq_along(sets)) features[match(sets[[j]],stat$gene), j] <- 1
    }
  }
  if(!is.data.frame(stat)||anyDuplicated(names(stat))||!all(c("gene","chromosome","score")%in%names(stat)))
    .gsuite_stop("stat must be a data frame with gene, chromosome and score columns.")
  .gsea_ids(stat$gene,"gene")
  chromosome<-as.character(stat$chromosome)
  if(anyNA(chromosome)||any(!nzchar(trimws(chromosome)))||length(unique(chromosome))<2L)
    .gsuite_stop("chromosome must contain at least two nonempty held-out groups.")
  if(!is.numeric(stat$score)||any(!is.finite(stat$score))) .gsuite_stop("score must be finite numeric values.")
  set_features <- is.list(features) && !is.data.frame(features)
  if(set_features && is.null(sets)) .gsuite_stop("Supply binary memberships through sets; features must be a numeric matrix.")
  if(set_features) .gmap_named(features,"features") else features<-.gsea_matrix(features,stat$gene,"features")
  feature_names <- if(set_features) names(features) else colnames(features)
  if(!is.character(control_features)||anyNA(control_features)||anyDuplicated(control_features)||any(!control_features%in%feature_names))
    .gsuite_stop("control_features must be unique feature names.")
  if(is.null(covariates))covariates<-matrix(numeric(),nrow(stat),0,dimnames=list(stat$gene,character())) else {
    if(is.data.frame(covariates)){
      if(!identical(covariates$gene,stat$gene)||anyDuplicated(names(covariates))) .gsuite_stop("covariates$gene must match stat exactly.")
      x<-covariates[setdiff(names(covariates),"gene")]
      if(!length(x)||!all(vapply(x,is.numeric,logical(1)))) .gsuite_stop("Covariates must be numeric; encode factors explicitly.")
      covariates<-as.matrix(x);rownames(covariates)<-stat$gene
    }
    covariates<-.gsea_matrix(covariates,stat$gene,"covariates")
  }
  if (method != "rr") return(.pops_bayesian(stat,features,covariates,sampling,
    control_features,control,method,pi,nit,nburn,nthin,nrun,seed,
    updateB,updateE,updatePi,prior,call,workflow_start))
  if (length(prior) || !is.null(pi)) .gsuite_stop("pi and prior require a Bayesian PoPS method.")
  defaults<-list(marginal_p_threshold=.05,ridge_penalties=10^seq(-2,10,length.out=25),
    maximum_features=4096L,maximum_selected_per_fold=256L,maximum_union_features=1024L,
    maximum_covariates=64L,maximum_block_size=2048L,
    maximum_preparation_bytes=512*1024^2,minimum_reciprocal_condition=1e-10)
  cc<-.gmap_settings(control,defaults,"control")
  .gmap_probability(cc$marginal_p_threshold,"marginal_p_threshold")
  if(!is.numeric(cc$ridge_penalties)||!length(cc$ridge_penalties)||any(!is.finite(cc$ridge_penalties))||any(cc$ridge_penalties<=0)) .gsuite_stop("ridge_penalties must be positive finite values.")
  for(nm in c("maximum_features","maximum_selected_per_fold","maximum_union_features","maximum_block_size")) .gmap_integer(cc[[nm]],nm)
  .gmap_integer(cc$maximum_covariates,"maximum_covariates",0L)
  if(cc$maximum_selected_per_fold>cc$maximum_union_features) .gsuite_stop("maximum_union_features must be at least maximum_selected_per_fold.")
  if(!is.numeric(cc$maximum_preparation_bytes)||length(cc$maximum_preparation_bytes)!=1L||!is.finite(cc$maximum_preparation_bytes)||cc$maximum_preparation_bytes<1) .gsuite_stop("maximum_preparation_bytes must be positive and finite.")
  if(!is.numeric(cc$minimum_reciprocal_condition)||length(cc$minimum_reciprocal_condition)!=1L||!is.finite(cc$minimum_reciprocal_condition)||cc$minimum_reciprocal_condition<=0||cc$minimum_reciprocal_condition>=1) .gsuite_stop("minimum_reciprocal_condition must be in (0, 1).")
  ss<-.gsea_sampling(sampling,stat$gene,"trait",FALSE,cc)
  native<-cpp_gsea_pops(transform(stat,chromosome=chromosome),features,covariates,control_features,ss,cc)
  predictions<-data.frame(gene=stat$gene,chromosome=chromosome,score=native$scores,stringsAsFactors=FALSE)
  predictions$genome_rank<-rank(-predictions$score,ties.method="min")
  predictions$chromosome_rank<-ave(-predictions$score,predictions$chromosome,FUN=function(x)rank(x,ties.method="min"))
  scaling<-data.frame(feature=colnames(features),mean=native$feature_means,scale=native$feature_scales,stringsAsFactors=FALSE)
  structure(list(call=match.call(),predictions=predictions,folds=native$folds,
    features=native$features,scaling=scaling,genes=stat$gene,chromosomes=chromosome,
    control_features=control_features,control=cc,sampling=ss,
    provenance=list(method="PoPS-style LOCO ridge",interpretation=native$interpretation,
      gene_evidence="Caller-supplied checked gene score",covariance=ss$gene_provenance)),
    class=c("pops_fit","list"))
}

#' @rdname pops
#' @export
print.pops_fit <- function(x,...){cat("PoPS-style LOCO ranking:",nrow(x$predictions),"genes |",nrow(x$folds),"chromosomes\n");print(utils::head(x$predictions[order(x$predictions$genome_rank),],10L),row.names=FALSE);invisible(x)}
#' @rdname pops
#' @export
summary.pops_fit <- function(object,...){list(gene_count=nrow(object$predictions),chromosome_count=nrow(object$folds),predictions=object$predictions,folds=object$folds,features=object$features,provenance=object$provenance)}
#' @rdname pops
#' @param x,object A `pops_fit`.
#' @param n Maximum genes to display.
#' @param ... Additional graphics arguments.
#' @export
plot.pops_fit <- function(x,n=20,...){.gmap_integer(n,"n");take<-utils::head(order(x$predictions$genome_rank),n);graphics::dotchart(rev(x$predictions$score[take]),labels=rev(x$predictions$gene[take]),xlab="PoPS-style priority score",...);invisible(x$predictions[take,,drop=FALSE])}
