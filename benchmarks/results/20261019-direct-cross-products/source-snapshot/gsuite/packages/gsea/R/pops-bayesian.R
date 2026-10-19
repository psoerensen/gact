# Chromosome splitting and prediction are R workflow responsibilities; all
# Bayesian cross-products, nuisance projection and sampling use native gsea/gbayes.
.pops_bayesian <- function(stat,features,covariates,sampling,control_features,
    control,method,pi,nit,nburn,nthin,nrun,seed,updateB,updateE,updatePi,prior,call,
    workflow_start=proc.time()[["elapsed"]]) {
  flat <- names(call)[-1L]
  assign_control <- function(name,value,arguments) {
    if (!is.null(control[[name]]) && any(arguments %in% flat))
      .gsuite_stop(paste("Specify flat arguments or control$",name,", not both.",sep=""))
    if (is.null(control[[name]])) control[[name]] <<- value
  }
  .gmap_integer(nrun,"nrun")
  assign_control("sampling_sweeps",nit,"nit")
  assign_control("burnin",nburn,"nburn")
  assign_control("thinning",nthin,"nthin")
  assign_control("seeds",seed+seq_len(nrun)-1L,c("seed","nrun"))
  assign_control("estimate_effect_variance",updateB,"updateB")
  assign_control("estimate_residual_variance",updateE,"updateE")
  assign_control("estimate_weights",updatePi,"updatePi")
  threshold <- control$marginal_p_threshold
  control$marginal_p_threshold <- NULL
  if (!is.null(threshold)) .gmap_probability(threshold,"threshold")
  if (is.null(sampling)) sampling <- list(gene_covariance="independent")
  cc <- .gsea_control(control,method,FALSE)
  ss <- .gsea_sampling(sampling,stat$gene,"trait",FALSE,cc)
  chromosome <- as.character(stat$chromosome)
  if (any(vapply(ss$blocks,function(b)
      length(unique(chromosome[b$indices])) != 1L,logical(1))))
    .gsuite_stop("PoPS covariance blocks must not span chromosomes.")
  folds <- unique(chromosome)
  set_features <- is.list(features) && !is.data.frame(features)
  feature_names <- if(set_features) names(features) else colnames(features)
  columns_of <- function(x) if(set_features) names(x) else colnames(x)
  feature_count <- function(x) if(set_features) length(x) else ncol(x)
  if(set_features && (!ss$independent || ncol(covariates) || !is.null(threshold)))
    .gsuite_stop("Direct PoPS sets require independent genes, no additional covariates and no marginal screening.")
  predictions <- numeric(nrow(stat)); fits <- selections <- vector("list",length(folds))
  fold_timing <- vector("list",length(folds))
  interface_preparation <- proc.time()[["elapsed"]]-workflow_start
  names(fits) <- names(selections) <- folds
  for (i in seq_along(folds)) {
    fold_start <- proc.time()[["elapsed"]]
    train <- chromosome != folds[i]; test <- !train
    genes <- stat$gene[train]
    if (sum(train) <= ncol(covariates)+2L) .gsuite_stop("Too few training genes in a PoPS fold.")
    if(set_features) {
      training_sets <- lapply(features,intersect,y=genes)
      means <- lengths(training_sets)/length(genes)
      varying <- lengths(training_sets)>0L & lengths(training_sets)<length(genes)
      x <- training_sets[varying]
    } else {
      means <- colMeans(features[train,,drop=FALSE])
      varying <- vapply(seq_len(ncol(features)),function(j)
        any(features[train,j] != means[j]),logical(1))
      x <- features[train,varying,drop=FALSE]
    }
    if (any(!varying[match(control_features,feature_names)]))
      .gsuite_stop("A control feature is constant in a training fold.")
    if (!feature_count(x)) .gsuite_stop("No variable features in a PoPS training fold.")
    fold_sampling <- list(gene_covariance="independent",gene_provenance=ss$gene_provenance)
    if (!ss$independent) {
      blocks <- ss$blocks[vapply(ss$blocks,function(b) train[b$indices[1L]],logical(1))]
      fold_sampling$gene_covariance <- lapply(blocks,`[[`,"covariance")
    }
    p <- setNames(rep(NA_real_,feature_count(x)),columns_of(x))
    if (!is.null(threshold)) {
      # Marginal tests fit in bounded chunks; no held-out response enters selection.
      width <- min(200L,sum(train)-ncol(covariates)-2L)
      for (start in seq.int(1L,ncol(x),by=width)) {
        columns <- seq.int(start,min(ncol(x),start+width-1L))
        screen <- gsea(stat[train,,drop=FALSE],x[,columns,drop=FALSE],"competitive",
          covariates=if(ncol(covariates)) covariates[train,,drop=FALSE] else NULL,sampling=fold_sampling,
          control=list(joint=FALSE,maximum_covariates=cc$maximum_covariates,
            maximum_block_size=cc$maximum_block_size,
            maximum_preparation_bytes=cc$maximum_preparation_bytes,
            minimum_reciprocal_condition=cc$minimum_reciprocal_condition))
        p[columns] <- screen$pathways$p_two_sided[match(colnames(x)[columns],screen$pathways$pathway)]
      }
      keep <- is.finite(p) & p < threshold | names(p) %in% control_features
      x <- x[,keep,drop=FALSE]; p <- p[keep]
      if (!ncol(x)) .gsuite_stop("No features selected in a PoPS training fold.")
    }
    fp <- prior
    if (!is.null(pi) && !is.null(fp[[if(method=="bayesc") "inclusion_probability" else "weights"]]))
      .gsuite_stop("Specify pi or its nested prior setting, not both.")
    if (method == "bayesc") {
      if (is.null(fp$inclusion_probability)) fp$inclusion_probability <- if(is.null(pi)) .001 else pi
      .gmap_probability(fp$inclusion_probability,"pi",TRUE)
      mass <- fp$inclusion_probability; components <- 2L
    } else {
      if (is.null(fp$weights)) {
        weights <- if(is.null(pi)) .001 else pi
        if (length(weights)==1L) {
          .gmap_probability(weights,"pi",TRUE)
          weights <- c(1-weights,rep(weights/3,3L))
        }
        fp$weights <- weights
      }
      if (is.null(fp$variance_multipliers)) fp$variance_multipliers <- c(0,.01,.1,1)
      if (length(fp$weights) != length(fp$variance_multipliers))
        .gsuite_stop("BayesR weights and variance multipliers must align.")
      mass <- sum(fp$weights*fp$variance_multipliers); components <- length(fp$weights)
    }
    if (!is.finite(mass) || mass <= 0) .gsuite_stop("Bayesian active prior mass must be positive.")
    if (is.null(fp$residual_variance)) fp$residual_variance <- .5
    if (is.null(fp$effect_variance)) fp$effect_variance <- .5/(feature_count(x)*mass)
    if (cc$estimate_effect_variance && is.null(fp$effect_variance_prior))
      fp$effect_variance_prior <- list(df=4,scale=fp$effect_variance)
    if (cc$estimate_residual_variance && is.null(fp$residual_variance_prior))
      fp$residual_variance_prior <- list(df=4,scale=fp$residual_variance)
    if (cc$estimate_weights && is.null(fp$weight_prior)) fp$weight_prior <- rep(1,components)
    selection_seconds <- proc.time()[["elapsed"]]-fold_start
    fit <- gsea(stat[train,,drop=FALSE],x,method,
      covariates=if(ncol(covariates)) covariates[train,,drop=FALSE] else NULL,
      sampling=fold_sampling,prior=fp,control=cc)
    prediction_start <- proc.time()[["elapsed"]]
    selected <- columns_of(x)
    beta <- fit$pathways$mean[match(selected,fit$pathways$pathway)]
    if(set_features) {
      heldout <- stat$gene[test]
      predicted <- rep(-sum(means[selected]*beta),length(heldout))
      for(j in seq_along(selected)) {
        at <- match(features[[selected[j]]],heldout,nomatch=0L)
        at <- at[at>0L];predicted[at] <- predicted[at]+beta[j]
      }
      predictions[test] <- predicted
    } else predictions[test] <- drop(sweep(features[test,selected,drop=FALSE],2L,
      means[selected],`-`) %*% beta)
    selections[[i]] <- data.frame(chromosome=folds[i],feature=selected,
      marginal_p=p,mean=beta,pip=fit$pathways$pip[match(selected,fit$pathways$pathway)],
      training_mean=means[selected],row.names=NULL)
    fold_timing[[i]] <- data.frame(chromosome=folds[i],
      selection_preparation_seconds=selection_seconds,
      model_preparation_seconds=fit$timing$preparation_seconds,
      blr_fit_seconds=fit$timing$blr_fit_seconds,
      posterior_summary_seconds=fit$timing$posterior_summary_seconds,
      prediction_seconds=proc.time()[["elapsed"]]-prediction_start,
      preparation_route=fit$timing$native[[1]]$preparation_route)
    # Retain posterior and diagnostics, rather than a feature matrix for each fold.
    fit$sets <- fit$covariates <- fit$genes <- NULL
    fit$sampling$blocks <- NULL
    fit$sampling$gene_covariance <- NULL
    fits[[i]] <- fit
  }
  prediction <- data.frame(gene=stat$gene,chromosome=chromosome,score=predictions)
  prediction$genome_rank <- rank(-prediction$score,ties.method="min")
  prediction$chromosome_rank <- ave(-prediction$score,chromosome,
    FUN=function(x) rank(x,ties.method="min"))
  fold_table <- data.frame(chromosome=folds,
    training_gene_count=vapply(folds,function(chr) sum(chromosome!=chr),integer(1)),
    selected_features=vapply(selections,nrow,integer(1)),method=method)
  fold_timing <- do.call(rbind,fold_timing)
  timing <- list(gene_statistics=list(computed=FALSE,seconds=NA_real_,
    source="Caller-supplied gene scores reused across folds; upstream computation excluded"),
    interface_preparation_seconds=interface_preparation,folds=fold_timing,
    preparation_seconds=interface_preparation+sum(fold_timing$selection_preparation_seconds+fold_timing$model_preparation_seconds),
    blr_fit_seconds=sum(fold_timing$blr_fit_seconds),
    posterior_summary_seconds=sum(fold_timing$posterior_summary_seconds),
    prediction_seconds=sum(fold_timing$prediction_seconds),
    total_seconds=proc.time()[["elapsed"]]-workflow_start,
    scope="Chromosome-held-out validation; sum of fold fits, not one full-data BLR fit")
  structure(list(call=call,predictions=prediction,folds=fold_table,
    timing=timing,
    features=do.call(rbind,selections),fits=fits,genes=stat$gene,chromosomes=chromosome,
    control_features=control_features,control=cc,sampling=list(
      independent=ss$independent,gene_provenance=ss$gene_provenance),
    provenance=list(method=paste("PoPS-style LOCO",method),
      interpretation="Chromosome-held-out feature predictions; not causal probabilities",
      gene_evidence="Caller-supplied checked gene score",covariance=ss$gene_provenance)),
    class=c("pops_fit","list"))
}
