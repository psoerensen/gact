.gsuite_stop <- function(message, call. = FALSE) stop(message, call. = call.)

.gsuite_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    .gsuite_stop(sprintf("\u0060%s\u0060 must be TRUE or FALSE.", name))
  value
}

.gmap_named <- function(x, what) {
  if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))
    .gsuite_stop(paste0("\u0060", what, "\u0060 must be a nonempty uniquely named list."))
}

.gmap_settings <- function(x, defaults, what) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      anyDuplicated(names(x)) || any(!names(x) %in% names(defaults)))))
    .gsuite_stop(paste("Unknown or unnamed", what, "setting."))
  defaults[names(x)] <- x
  defaults
}

.gmap_integer <- function(x, name, minimum = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < minimum || x > .Machine$integer.max)
    .gsuite_stop(paste("Invalid integer control:", name))
  invisible(x)
}

.gmap_probability <- function(x, name, open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1 ||
      (open && (x == 0 || x == 1)))
    .gsuite_stop(paste("Invalid probability:", name))
  invisible(x)
}

.gmap_select_one <- function(value, available, what) {
  if (is.null(value) && length(available) == 1L) return(available)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% available)
    .gsuite_stop(paste("Select one", what, "from:", paste(available, collapse = ", ")))
  value
}

.ldlist_label <- function(value, name, required = FALSE) {
  if (is.null(value) && !required) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value)))
    .gsuite_stop(sprintf("\u0060%s\u0060 must be a nonempty string.", name))
  value
}

.gbayes_covariance <- function(x, traits, what) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      !identical(dim(x), rep(length(traits), 2L)) ||
      !identical(rownames(x), traits) || !identical(colnames(x), traits))
    .gsuite_stop(paste(what, "must be a finite matrix with row/column names in trait order."))
  x
}

.gbayes_vector <- function(x, labels, what, positive = TRUE,
                           normalized = FALSE, named = FALSE) {
  if (!is.numeric(x) || length(x) != length(labels) || any(!is.finite(x)) ||
      any(if (positive) x <= 0 else x < 0) ||
      (named && !identical(names(x), labels)) ||
      (!is.null(names(x)) && !identical(names(x), labels)) ||
      (normalized && abs(sum(x) - 1) > 1e-10))
    .gsuite_stop(paste("Invalid", what, "values or labels; use the declared order."))
  stats::setNames(x, labels)
}

.gmap_learning <- function(prior, control, method) {
  positive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  for (pair in list(c("effect_variance_prior", "estimate_effect_variance"),
                    c("residual_variance_prior", "estimate_residual_variance"))) {
    x <- prior[[pair[[1L]]]]
    if (control[[pair[[2L]]]]) {
      if (!is.list(x) || length(x) != 2L || anyDuplicated(names(x)) ||
          !setequal(names(x), c("df", "scale")) || !positive(x$df) ||
          !positive(x$scale))
        .gsuite_stop(paste("Learning requires", pair[[1L]],
          "= list(df = ..., scale = ...) with positive values."))
    } else if (!is.null(x))
      .gsuite_stop(paste(pair[[1L]], "is only used when its learning flag is TRUE."))
  }
  weights <- prior$weight_prior
  if (control$estimate_weights) {
    components <- if (method == "bayesc") 2L else length(prior$weights)
    if (!is.numeric(weights) || length(weights) != components ||
        any(!is.finite(weights)) || any(weights <= 0))
      .gsuite_stop("weight_prior must contain one positive concentration per component.")
  } else if (!is.null(weights))
    .gsuite_stop("weight_prior is only used with estimate_weights=TRUE.")
  invisible(NULL)
}

.gsea_patterns <- function(traits) {
  if (length(traits) < 2L || length(traits) > 12L)
    .gsuite_stop("Default complete patterns support 2 to 12 traits.")
  pattern <- as.matrix(expand.grid(rep(list(0:1), length(traits)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  storage.mode(pattern) <- "integer"
  colnames(pattern) <- traits
  rownames(pattern) <- apply(pattern, 1L, paste0, collapse = "")
  pattern
}

#' Gene-set and pathway analysis using native gsea
#'
#' Gene statistics and sets must already share an explicit background. Preparation
#' of SNP-to-gene evidence is separate. No genes are silently dropped or reordered.
#' @param stat Data frame with unique character gene IDs in gene and selected
#'   (logical, ORA), score (preranked/regression), or p_value (MAGMA). A named list
#'   supplies multiple traits. Bayesian traits are fitted jointly; classical
#'   methods are fitted separately, with multiplicity adjustment within each trait.
#' @param sets Named list of character gene vectors, or a finite numeric matrix
#'   with gene row names and feature column names (competitive/Bayesian only).
#' @param method Explicit choice: ora, preranked, competitive, magma, bayesc, bayesr.
#' @param covariates Optional numeric matrix with gene row names, or data frame
#'   with gene and numeric covariate columns. Shared across traits; no intercept.
#'   MAGMA uses its gene metadata instead; see Details.
#' @param sampling For regression, list(gene_covariance=...), where the value is
#'   "independent", a labelled matrix, or named list of labelled block matrices
#'   partitioning all genes. Joint Bayesian fits additionally require a labelled
#'   trait_covariance matrix in names(stat) order, on original response units.
#'   Optional gene_provenance and trait_provenance describe these working models.
#' @param prior Bayesian priors following gbayes: scalar residual_variance and
#'   effect_variance are required; joint fits require effect_covariance on the
#'   prepared scale. Inclusion/mixture priors and optional learning priors are
#'   documented in the pathway guide. Empty for classical methods.
#' @param control Named method-specific settings; unknown/inapplicable names fail.
#'   ORA/preranked/classical adjustment is none (default), BH or BY. Preranked uses
#'   seed=1, replicates=999, exponent=1, minimum_size=1, maximum_size=500, threads=1.
#'   Competitive joint=TRUE; MAGMA joint=FALSE. Bayesian controls include burnin=100,
#'   sampling_sweeps=1000, thinning=1, seeds=c(1,2), threads=1,
#'   retain_diagnostics=TRUE, scale_response=TRUE, scale_features=TRUE and explicit
#'   parameter-learning flags. See the pathway guide for limits and full settings.
#' @return A serializable gsea_fit with pathways, parameters, diagnostics, details,
#'   preparation and provenance; Bayesian fits also include chains and joint PIPs.
#'   Their timing list separates preparation and native BLR; supplied upstream
#'   gene-statistic computation is not timed and is recorded as NA.
#' @details Regression requires explicit response-scale gene covariance. MAGMA
#'   instead consumes quadratic-score correlations and applies its documented
#'   working-correlation power. Enabled MAGMA corrections require marker_count,
#'   effective_parameters, sample_size and mean_minor_allele_count columns.
#'   Defaults enable all four corrections; their raw/log columns are handled
#'   natively. Classical MAGMA defaults to outcome-only internal adjustment.
#'   Bayesian fits adjust both response and features and record scales. Priors,
#'   native parameter summaries and chains use prepared coordinates; pathway
#'   means/variances use original units. Joint residual covariance stays fixed.
#'   PIPs are not p-values. Chain spread is not an ESS/R-hat convergence check.
#' @export
gsea <- function(stat, sets, method, covariates=NULL, sampling=NULL,
                  prior=list(), control=list()) {
  workflow_start <- proc.time()[["elapsed"]]
  call <- match.call()
  if(missing(method) || !is.character(method) || length(method)!=1L || is.na(method))
    .gsuite_stop("Supply one explicit gsea method.")
  method <- match.arg(tolower(method),c("ora","preranked","competitive","magma","bayesc","bayesr"))
  bayes <- method %in% c("bayesc","bayesr")
  if(is.data.frame(stat)) stat <- list(trait=stat)
  .gmap_named(stat,"stat")
  traits <- names(stat); joint_traits <- bayes && length(stat)>1L
  column <- if(method=="ora") "selected" else if(method=="magma") "p_value" else "score"
  for(s in stat) {
    if(!is.data.frame(s) || anyDuplicated(names(s)) || !all(c("gene",column) %in% names(s)))
      .gsuite_stop(paste("stat requires gene and",column,"columns."))
    .gsea_ids(s$gene,"gene")
    if(column=="selected") {
      if(!is.logical(s$selected) || anyNA(s$selected)) .gsuite_stop("selected must be logical without missing values, across the complete background.")
    } else if(!is.numeric(s[[column]]) || any(!is.finite(s[[column]]))) .gsuite_stop(paste(column,"must be finite numeric values."))
  }
  genes <- stat[[1]]$gene
  if(any(!vapply(stat,function(s) identical(s$gene,genes),logical(1)))) .gsuite_stop("All traits must have the same gene IDs and order.")
  matrix_sets <- is.matrix(sets)
  if(matrix_sets) {
    if(!bayes && method!="competitive") .gsuite_stop("Feature matrices are supported by competitive and Bayesian methods only.")
    sets <- .gsea_matrix(sets,genes,"sets")
    feature_names <- colnames(sets)
  } else {
    .gmap_named(sets,"sets");feature_names <- names(sets)
    for(s in sets) {
      .gsea_ids(s,"set members")
      if(any(!s %in% genes)) .gsuite_stop("Every set member must be in the declared gene background.")
    }
  }
  regression <- !method %in% c("ora","preranked")
  cc <- .gsea_control(control,method,joint_traits)
  if(!regression && (!is.null(covariates)||!is.null(sampling)||length(prior)))
    .gsuite_stop("ORA/preranked do not use covariates, sampling or prior.")
  if(!bayes && length(prior)) .gsuite_stop("prior is for Bayesian methods only.")
  if(method=="magma" && !is.null(covariates)) .gsuite_stop("MAGMA uses declared gene metadata; use competitive for custom covariates.")
  if(is.null(covariates)) covariates <- matrix(numeric(),length(genes),0,dimnames=list(genes,character())) else {
    if(is.data.frame(covariates)) {
      if(!identical(covariates$gene,genes)||anyDuplicated(names(covariates))) .gsuite_stop("covariates$gene must match stat exactly.")
      x <- covariates[setdiff(names(covariates),"gene")]
      if(!length(x)||!all(vapply(x,is.numeric,logical(1)))) .gsuite_stop("Covariates must be numeric; encode factors explicitly.")
      covariates <- as.matrix(x);rownames(covariates) <- genes
    }
    covariates <- .gsea_matrix(covariates,genes,"covariates")
  }
  ss <- if(regression) .gsea_sampling(sampling,genes,traits,joint_traits,cc) else list()
  if(regression) {
    if(length(feature_names)>cc$maximum_features || ncol(covariates)>cc$maximum_covariates)
      .gsuite_stop("Design exceeds configured feature/covariate limits.")
    direct_sets <- bayes && !joint_traits && !matrix_sets && ss$independent && ncol(covariates)==0L
    if(!direct_sets && length(genes)*length(feature_names)*8>cc$maximum_preparation_bytes)
      .gsuite_stop("Feature payload exceeds maximum_preparation_bytes.")
  }
  if(method=="magma") {
    needed <- c(marker_count=cc$correct_size,effective_parameters=cc$correct_density,
                sample_size=cc$correct_sample_size,mean_minor_allele_count=cc$correct_mac)
    for(k in seq_along(stat)) for(nm in names(needed)) {
      if(needed[[nm]] && (is.null(stat[[k]][[nm]])||!is.numeric(stat[[k]][[nm]])||any(!is.finite(stat[[k]][[nm]]))||any(stat[[k]][[nm]]<=0)))
        .gsuite_stop(paste("Enabled MAGMA correction requires positive",nm))
      if(!needed[[nm]]) stat[[k]][[nm]] <- rep(0,length(genes))
    }
  }
  p <- if(bayes) .gsea_prior(prior,cc,method,traits,joint_traits) else list()
  r_preparation_seconds <- proc.time()[["elapsed"]]-workflow_start
  native <- cpp_gsea(stat,sets,covariates,method,ss,p,cc)
  pathways <- do.call(rbind,lapply(native$fits,`[[`,"pathways"));rownames(pathways)<-NULL
  fields <- c("parameters","diagnostics","details","preparation","chains")
  out <- lapply(fields,function(nm) lapply(native$fits,`[[`,nm));names(out)<-fields
  timing <- if(bayes) {
    stages <- lapply(native$fits,`[[`,"timing")
    stage_sum <- function(name) sum(vapply(stages,`[[`,numeric(1),name))
    list(gene_statistics=list(computed=FALSE,seconds=NA_real_,source="Caller-supplied gene scores; upstream computation excluded"),
      r_preparation_seconds=r_preparation_seconds,native=stages,
      preparation_seconds=r_preparation_seconds+stage_sum("preparation_seconds"),
      blr_fit_seconds=stage_sum("blr_fit_seconds"),
      posterior_summary_seconds=stage_sum("posterior_summary_seconds"),
      total_seconds=proc.time()[["elapsed"]]-workflow_start,
      fit_boundary="Native fit call including checks, sampling, accumulation and coefficient back-transformation")
  } else NULL
  structure(c(list(call=call,method=method,traits=traits,genes=genes,timing=timing,
    prior=p,control=cc,sets=sets,covariates=covariates,sampling=ss,
    pathways=pathways,joint=if(joint_traits) native$fits$joint$joint else NULL,
    provenance=list(native_version=native$native_version,joint_traits=native$joint_traits,
      parameter_units=if(bayes) "Prepared coordinates" else NULL,
      pathway_units=if(bayes) "Original response/feature units" else if(method=="magma") "Prepared probit-score regression" else "Method-specific statistics")),out),class="gsea_fit")
}

.gsea_ids <- function(x,what) {
  if(!is.character(x)||!length(x)||anyNA(x)||any(!nzchar(trimws(x)))||anyDuplicated(x))
    .gsuite_stop(paste(what,"must contain unique nonempty character IDs."))
}
.gsea_matrix <- function(x,genes,what) {
  if(!is.matrix(x)||!is.numeric(x)||!ncol(x)||!identical(rownames(x),genes)||any(!is.finite(x)))
    .gsuite_stop(paste(what,"must be a finite numeric matrix with row names in gene order."))
  .gsea_ids(colnames(x),paste(what,"column names"));storage.mode(x)<-"double";x
}
.gsea_control <- function(control,method,multi) {
  classical <- !method %in% c("bayesc","bayesr")
  defaults <- if(classical) list(adjustment="none") else list()
  if(method=="preranked") defaults <- c(defaults,list(seed=1,replicates=999,threads=1,
    exponent=1,direction="standard",minimum_size=1,maximum_size=500,maximum_scratch_bytes=256*1024^2))
  if(!method %in% c("ora","preranked")) defaults <- c(defaults,list(maximum_features=if(classical) 255 else 4096,
    maximum_covariates=128,maximum_block_size=2048,maximum_preparation_bytes=512*1024^2,minimum_reciprocal_condition=1e-10))
  if(method %in% c("competitive","magma")) defaults <- c(defaults,list(joint=method=="competitive",maximum_coefficients=256))
  if(method=="magma") defaults <- c(defaults,list(correct_size=TRUE,correct_density=TRUE,
    correct_sample_size=TRUE,correct_mac=TRUE,truncate=TRUE,residualize_internal=TRUE,
    lower_z_bound=-3,upper_z_sd=6,covariate_sd_bound=5,correlation_power=1.1))
  if(!classical) {
    defaults <- c(defaults,list(burnin=100,sampling_sweeps=1000,thinning=1,seeds=c(1,2),threads=1,
      retain_diagnostics=TRUE,scale_response=TRUE,scale_features=TRUE,pip_range_threshold=.1))
    defaults <- c(defaults,if(multi) list(estimate_covariance=FALSE,estimate_pattern_weights=FALSE,
      estimate_component_weights=FALSE,retain_pattern_probabilities=FALSE,memory_limit_bytes=256*1024^2) else
      list(estimate_effect_variance=FALSE,estimate_residual_variance=FALSE,estimate_weights=FALSE))
  }
  if(method=="magma") defaults[c("minimum_reciprocal_condition","maximum_coefficients")] <- NULL
  cc <- .gmap_settings(control,defaults,"control")
  if(method=="magma" && cc$maximum_block_size>2048) .gsuite_stop("MAGMA supports covariance blocks up to 2048 genes.")
  for(nm in names(defaults)[vapply(defaults,is.logical,logical(1))]) cc[[nm]] <- .gsuite_flag(cc[[nm]],nm)
  ints <- intersect(names(cc),c("burnin","sampling_sweeps","thinning","threads","replicates","minimum_size","maximum_size",
    "maximum_scratch_bytes","maximum_features","maximum_covariates","maximum_block_size","maximum_preparation_bytes","maximum_coefficients","memory_limit_bytes"))
  for(nm in ints) .gmap_integer(cc[[nm]],nm,minimum=if(nm %in% c("burnin","maximum_covariates")) 0 else 1)
  if(classical) cc$adjustment <- match.arg(cc$adjustment,c("none","BH","BY"))
  if(method=="preranked") {
    .gmap_integer(cc$seed,"seed",minimum=0)
    if(!is.numeric(cc$exponent)||length(cc$exponent)!=1||!is.finite(cc$exponent)||cc$exponent<0) .gsuite_stop("exponent must be finite and nonnegative.")
    cc$direction <- match.arg(cc$direction,c("standard","positive","negative"))
    if(cc$maximum_size<cc$minimum_size) .gsuite_stop("maximum_size must be at least minimum_size.")
  }
  if(!classical) {
    if(!is.numeric(cc$seeds)||!length(cc$seeds)||any(!is.finite(cc$seeds))||any(cc$seeds<0|cc$seeds>2^32-1|cc$seeds!=floor(cc$seeds))||anyDuplicated(cc$seeds)) .gsuite_stop("seeds must be distinct integers from 0 to 2^32-1.")
    .gmap_probability(cc$pip_range_threshold,"pip_range_threshold")
    if(multi && method=="bayesc" && cc$estimate_component_weights) .gsuite_stop("BayesC has one positive scale; omit component-weight learning.")
  }
  reals <- intersect(names(cc),c("minimum_reciprocal_condition","lower_z_bound","upper_z_sd","covariate_sd_bound","correlation_power"))
  for(nm in reals) if(!is.numeric(cc[[nm]])||length(cc[[nm]])!=1||!is.finite(cc[[nm]])) .gsuite_stop(paste("Invalid",nm))
  cc
}
.gsea_sampling <- function(sampling,genes,traits,multi,control) {
  if(is.null(sampling)) .gsuite_stop("Supply sampling$gene_covariance; use 'independent' only when intended.")
  ss <- .gmap_settings(sampling,list(gene_covariance=NULL,trait_covariance=NULL,
    gene_provenance="Caller-supplied gene working covariance",trait_provenance="Caller-supplied fixed residual trait covariance"),"sampling")
  for(nm in c("gene_provenance","trait_provenance")) .ldlist_label(ss[[nm]],nm,TRUE)
  value <- ss$gene_covariance
  if(identical(value,"independent")) ss$blocks <- list() else {
    blocks <- if(is.matrix(value)) list(all=value) else value
    .gmap_named(blocks,"gene_covariance blocks")
    seen <- character()
    ss$blocks <- lapply(blocks,function(b) {
      ids <- rownames(b);.gsea_ids(ids,"covariance gene IDs")
      if(!is.matrix(b)||!is.numeric(b)||!identical(colnames(b),ids)||!identical(dim(b),rep(length(ids),2L))||any(!is.finite(b))||any(!ids %in% genes)) .gsuite_stop("Covariance blocks must be finite square matrices with matching gene row/column names.")
      if(length(ids)>control$maximum_block_size) .gsuite_stop("Gene covariance block exceeds maximum_block_size.")
      seen <<- c(seen,ids)
      list(indices=as.integer(match(ids,genes)),covariance=b)
    })
    if(anyDuplicated(seen)||!setequal(seen,genes)) .gsuite_stop("Gene covariance blocks must partition the complete background exactly once.")
  }
  ss$independent <- identical(value,"independent")
  if(multi) ss$trait_covariance <- .gbayes_covariance(ss$trait_covariance,traits,"trait_covariance") else
    if(!is.null(ss$trait_covariance)) .gsuite_stop("trait_covariance is used only for joint Bayesian traits.")
  ss
}
.gsea_prior <- function(prior,cc,method,traits,multi) {
  positive <- function(x) is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>0
  if(!multi) {
    p <- list(residual_variance=NULL,effect_variance=NULL,effect_variance_prior=NULL,residual_variance_prior=NULL,weight_prior=NULL)
    p <- c(p,if(method=="bayesc") list(inclusion_probability=.5) else list(weights=c(.7,.1,.1,.1),variance_multipliers=c(0,.01,.1,1)))
    p <- .gmap_settings(prior,p,"prior")
    if(!positive(p$residual_variance)||!positive(p$effect_variance)) .gsuite_stop("Supply positive residual_variance and effect_variance on the prepared scale.")
    if(method=="bayesc") .gmap_probability(p$inclusion_probability,"inclusion_probability",open=TRUE) else {
      .gbayes_vector(p$weights,paste0("component_",seq_along(p$weights)),"weights",positive=FALSE,normalized=TRUE)
      if(length(p$weights)<2L||!is.numeric(p$variance_multipliers)||length(p$variance_multipliers)!=length(p$weights)||any(!is.finite(p$variance_multipliers))||p$variance_multipliers[1]!=0||any(diff(p$variance_multipliers)<=0)) .gsuite_stop("BayesR needs increasing variance_multipliers starting at zero, matching weights.")
    }
    .gmap_learning(p,cc,method);return(p)
  }
  p <- .gmap_settings(prior,list(effect_covariance=NULL,patterns=NULL,pattern_weights=NULL,
    component_scales=if(method=="bayesc") 1 else c(.01,.1,1),component_weights=if(method=="bayesc") 1 else rep(1/3,3),
    covariance_prior=NULL,pattern_prior=NULL,component_prior=NULL),"prior")
  p$effect_covariance <- .gbayes_covariance(p$effect_covariance,traits,"effect_covariance")
  if(is.null(p$patterns)) p$patterns <- .gsea_patterns(traits) else {
    x <- p$patterns
    if(!is.matrix(x)||!(is.numeric(x)||is.logical(x))||nrow(x)<2L||!identical(colnames(x),traits)||anyNA(x)||any(!x %in% c(0,1))||any(x[1,]!=0)||any(colSums(x)==0)) .gsuite_stop("patterns must be binary with trait columns, null first, and every trait represented.")
    codes <- apply(x,1,paste0,collapse="");if(anyDuplicated(codes)) .gsuite_stop("Patterns must be distinct.")
    if(is.null(rownames(x))) rownames(x)<-codes
    .gsea_ids(rownames(x),"pattern names");storage.mode(x)<-"integer";p$patterns<-x
  }
  if(is.null(p$pattern_weights)) p$pattern_weights <- rep(1/nrow(p$patterns),nrow(p$patterns))
  p$pattern_weights <- .gbayes_vector(p$pattern_weights,rownames(p$patterns),"pattern_weights",positive=FALSE,normalized=TRUE)
  s <- p$component_scales
  if(!is.numeric(s)||!length(s)||any(!is.finite(s))||any(s<=0)||any(diff(s)<=0)||(method=="bayesc"&&!identical(as.double(s),1))||(method=="bayesr"&&length(s)<2)) .gsuite_stop("BayesC requires scale 1; MT-BayesR uses an increasing positive scale ladder, without zero.")
  p$component_weights <- .gbayes_vector(p$component_weights,paste0("component_",seq_along(s)),"component_weights",positive=FALSE,normalized=TRUE)
  if(cc$estimate_covariance) {
    iw<-p$covariance_prior
    if(!is.list(iw)||!setequal(names(iw),c("df","scale"))||length(iw)!=2L||!positive(iw$df)||iw$df<=length(traits)-1) .gsuite_stop("covariance_prior requires df > traits-1 and a labelled scale matrix.")
    p$covariance_prior$scale <- .gbayes_covariance(iw$scale,traits,"covariance_prior scale")
  } else if(!is.null(p$covariance_prior)) .gsuite_stop("covariance_prior requires estimate_covariance=TRUE.")
  for(pair in list(c("pattern_prior","estimate_pattern_weights"),c("component_prior","estimate_component_weights"))) {
    labels <- if(pair[1]=="pattern_prior") rownames(p$patterns) else names(p$component_weights)
    if(cc[[pair[2]]]) p[[pair[1]]] <- .gbayes_vector(p[[pair[1]]],labels,pair[1]) else if(!is.null(p[[pair[1]]])) .gsuite_stop(paste(pair[1],"requires its learning flag."))
  }
  p
}


