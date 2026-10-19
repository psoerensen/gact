#' Pathway regression with familiar MAGMA-style arguments
#'
#' A thin interface to the existing competitive and Bayesian gsea engines. This
#' does not run the external MAGMA executable or manufacture MAGMA gene metadata.
#' @param stat Named numeric scores, a gene-by-trait matrix/data frame with gene
#'   row names, or a native gsea score data frame/list.
#' @param sets Named gene sets or a labelled feature matrix.
#' @param method "linear", "bayesC" or "bayesR".
#' @param type "joint" or "marginal" for classical pathway regression.
#' @param pi Starting active probability (BayesC), or mixture weights (BayesR).
#'   A scalar BayesR pi allocates active mass equally to its three nonzero scales.
#' @param nit,nburn,nthin,nrun,seed Sampling sweeps, burn-in, thinning, chains and seed.
#' @param sampling Explicit working covariance; defaults to independent genes
#'   and, for joint Bayesian traits, independent trait errors on original units.
#' @param prior,control Advanced native settings, including parameter learning.
#' @return A gsea_fit with pathways and the existing native diagnostics. Joint
#'   Bayesian fits also expose familiar b and PIP matrices.
#'   Eligible scalar binary-set fits use direct cross-products. Bayesian timing
#'   separates preparation and BLR; upstream gene-statistic time is unmeasured.
#' @export
magma <- function(stat, sets, method = "linear", type = c("joint", "marginal"),
                  pi = NULL, nit = 1000L, nburn = 100L, nthin = 1L,
                  nrun = 2L, seed = 10L, sampling = NULL,
                  prior = list(), control = list()) {
  interface_start <- proc.time()[["elapsed"]]
  type <- match.arg(type)
  method <- match.arg(tolower(method), c("linear", "bayesc", "bayesr"))
  native <- is.list(stat) && !is.data.frame(stat) && !is.numeric(stat)
  if (is.data.frame(stat) && all(c("gene", "score") %in% names(stat))) native <- TRUE
  if (native) {
    input <- if (is.data.frame(stat)) list(trait = stat) else stat
  } else {
    if (is.numeric(stat) && is.null(dim(stat)))
      stat <- matrix(stat, ncol = 1L, dimnames = list(names(stat), "trait"))
    if (is.data.frame(stat)) stat <- as.matrix(stat)
    if (!is.matrix(stat) || !is.numeric(stat) || !nrow(stat) || !ncol(stat))
      .gsuite_stop("stat must supply named numeric gene scores.")
    .gsea_ids(rownames(stat), "stat gene names")
    if (is.null(colnames(stat))) colnames(stat) <- paste0("trait", seq_len(ncol(stat)))
    .gsea_ids(colnames(stat), "stat trait names")
    input <- stats::setNames(lapply(seq_len(ncol(stat)), function(j)
      data.frame(gene = rownames(stat), score = stat[, j], stringsAsFactors = FALSE)), colnames(stat))
  }
  .gmap_named(input, "stat")
  genes <- input[[1L]]$gene
  if (!is.matrix(sets)) {
    sets <- gbase::mapSets(lapply(sets, unique), rsids = genes, index = FALSE)
    sets <- sets[lengths(sets) > 0L]
  }
  if (is.null(sampling)) sampling <- list(gene_covariance = "independent")
  if (method == "linear") {
    if (!is.null(pi)) .gsuite_stop("pi is used only for Bayesian methods.")
    if (length(prior)) .gsuite_stop("prior is used only for Bayesian methods.")
    if (is.null(control$joint)) control$joint <- type == "joint"
    return(gsea(input, sets, "competitive", sampling = sampling, control = control))
  }
  features <- if (is.matrix(sets)) ncol(sets) else length(sets)
  .gmap_integer(features, "number of features")
  .gmap_integer(nrun, "nrun")
  cc <- list(sampling_sweeps = nit, burnin = nburn, thinning = nthin,
    seeds = seed + seq_len(nrun) - 1L)
  flat <- names(match.call())[-1L]
  controls_for <- list(sampling_sweeps = "nit", burnin = "nburn",
    thinning = "nthin", seeds = c("seed", "nrun"))
  shared <- intersect(names(control), names(cc))
  if (any(vapply(controls_for[shared], function(x) any(x %in% flat), logical(1))))
    .gsuite_stop("Use flat MCMC arguments or their nested control settings, not both.")
  cc[shared] <- NULL
  control <- c(control, cc)
  if (length(input) == 1L) {
    if (method == "bayesc") {
      if (!is.null(pi) && !is.null(prior$inclusion_probability))
        .gsuite_stop("Specify either pi or prior$inclusion_probability, not both.")
      if (is.null(pi)) pi <- .01
      .gmap_probability(pi, "pi", TRUE)
      if (is.null(prior$inclusion_probability)) prior$inclusion_probability <- pi
      active <- prior$inclusion_probability
    } else {
      if (!is.null(pi) && !is.null(prior$weights))
        .gsuite_stop("Specify either pi or prior$weights, not both.")
      if (is.null(pi)) pi <- c(.99, .003, .003, .004)
      if (length(pi) == 1L) {
        .gmap_probability(pi, "pi", TRUE)
        pi <- c(1 - pi, rep(pi / 3, 3))
      }
      if (is.null(prior$weights)) prior$weights <- pi
      if (is.null(prior$variance_multipliers)) prior$variance_multipliers <- c(0, .01, .1, 1)
      active <- sum(prior$weights * prior$variance_multipliers)
    }
    if (is.null(prior$residual_variance)) prior$residual_variance <- .5
    if (is.null(prior$effect_variance)) prior$effect_variance <- .5 / (features * active)
  } else {
    traits <- names(input)
    if (is.null(pi)) pi <- .01
    .gmap_probability(pi, "pi", TRUE)
    if (is.null(sampling$trait_covariance)) {
      sampling$trait_covariance <- diag(vapply(input, function(x) stats::var(x$score), numeric(1)))
      dimnames(sampling$trait_covariance) <- list(traits, traits)
    }
    if (is.null(prior$effect_covariance)) {
      prior$effect_covariance <- diag(rep(.5 / (features * pi), length(traits)))
      dimnames(prior$effect_covariance) <- list(traits, traits)
    }
    if (is.null(prior$pattern_weights)) {
      patterns <- .gsea_patterns(traits)
      prior$pattern_weights <- apply(patterns, 1L, function(x) prod(ifelse(x == 1L, pi, 1 - pi)))
    }
  }
  interface_preparation <- proc.time()[["elapsed"]]-interface_start
  fit <- gsea(input, sets, method, sampling = sampling, prior = prior, control = control)
  fit$timing$interface_preparation_seconds <- interface_preparation
  fit$timing$preparation_seconds <- fit$timing$preparation_seconds+interface_preparation
  fit$b <- matrix(fit$pathways$mean, ncol = length(input),
    dimnames = list(if (is.matrix(sets)) colnames(sets) else names(sets), names(input)))
  fit$PIP <- matrix(fit$pathways$pip, ncol = length(input), dimnames = dimnames(fit$b))
  fit$timing$total_seconds <- proc.time()[["elapsed"]]-interface_start
  fit
}
