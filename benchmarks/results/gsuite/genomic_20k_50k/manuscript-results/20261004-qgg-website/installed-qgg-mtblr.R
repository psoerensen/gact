function (yy = NULL, Xy = NULL, XX = NULL, n = NULL, b = NULL, 
    h2 = NULL, pi = 0.001, models = NULL, pimodels = NULL, vg = NULL, 
    vb = NULL, ve = NULL, ssb_prior = NULL, sse_prior = NULL, 
    updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, 
    nue = 4, nit = 1000, nburn = 500, nthin = 1, method = "bayesC", 
    verbose = NULL) 
{
    ww <- lapply(XX, diag)
    wy <- lapply(Xy, as.vector)
    m <- mean(sapply(wy, length))
    nt <- length(wy)
    XXvalues <- lapply(XX, function(x) {
        as.list(as.data.frame(x))
    })
    XXindices <- lapply(1:m, function(x) {
        (1:m) - 1
    })
    if (method == "bayesC") 
        method <- 4
    if (method == "bayesR") 
        method <- 5
    if (is.null(b)) 
        b <- lapply(1:nt, function(x) {
            rep(0, m)
        })
    if (is.matrix(b)) 
        b <- split(b, rep(1:ncol(b), each = nrow(b)))
    if (is.null(models)) {
        models <- rep(list(0:1), nt)
        models <- t(do.call(expand.grid, models))
        models <- split(models, rep(1:ncol(models), each = nrow(models)))
    }
    if (is.character(models)) {
        if (models == "restrictive") {
            models <- list(rep(0, nt), rep(1, nt))
        }
    }
    if (is.null(pimodels)) 
        pimodels <- c(1 - pi, rep(pi/(length(models) - 1), length(models) - 
            1))
    vy <- diag(yy/(n - 1), nt)
    if (is.null(h2)) 
        h2 <- 0.5
    if (is.null(vg)) 
        vg <- diag(diag(vy) * h2)
    if (is.null(ve)) 
        ve <- diag(diag(vy) * (1 - h2))
    if (method < 4 && is.null(vb)) 
        vb <- diag((diag(vy) * h2)/m)
    if (method >= 4 && is.null(vb)) 
        vb <- diag((diag(vy) * h2)/(m * pi))
    if (method < 4 && is.null(ssb_prior)) 
        ssb_prior <- diag(((nub - 2)/nub) * (diag(vg)/m))
    if (method >= 4 && is.null(ssb_prior)) 
        ssb_prior <- diag(((nub - 2)/nub) * (diag(vg)/(m * pi)))
    if (is.null(sse_prior)) 
        sse_prior <- diag(((nue - 2)/nue) * diag(ve))
    seed <- sample.int(.Machine$integer.max, 1)
    fit <- .Call("_qgg_mtblr", wy = wy, ww = ww, yy = yy, b = b, 
        XXvalues = XXvalues, XXindices = XXindices, B = vb, E = ve, 
        ssb_prior = split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))), 
        sse_prior = split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))), 
        models = models, pi = pimodels, nub = nub, nue = nue, 
        updateB = updateB, updateE = updateE, updatePi = updatePi, 
        n = n, nit = nit, nburn = nburn, nthin = nthin, seed = seed, 
        method = as.integer(method))
    names(fit) <- c("bm", "dm", "wy", "r", "b", "d", "o", "vbs", 
        "vgs", "ves", "covb", "covg", "cove", "vb", "vg", "ve", 
        "pi", "pim", "pitrait", "pimarker")
    trait_names <- names(yy)
    if (is.null(trait_names)) 
        trait_names <- paste0("T", 1:nt)
    variable_names <- names(XXvalues[[1]])
    if (is.null(variable_names)) 
        variable_names <- paste0("V", 1:length(XXvalues[[1]]))
    for (i in 1:7) {
        fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
        rownames(fit[[i]]) <- variable_names
        colnames(fit[[i]]) <- trait_names
    }
    for (i in 8:10) {
        fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
        rownames(fit[[i]]) <- paste0("Iter", 1:nrow(fit[[i]]))
        colnames(fit[[i]]) <- trait_names
    }
    for (i in 11:16) {
        fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
        colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
    }
    fit[[17]] <- fit[[17]][[1]]
    fit[[18]] <- fit[[18]][[1]]
    names(fit[[17]]) <- sapply(models, paste, collapse = "_")
    names(fit[[18]]) <- sapply(models, paste, collapse = "_")
    if (method == 4) {
        fit <- fit[1:18]
    }
    if (method == 5) {
        fit[[19]] <- as.matrix(as.data.frame(fit[[19]]))
        rownames(fit[[19]]) <- c("0", "0.01", "0.1", "1.0")
        colnames(fit[[19]]) <- trait_names
        fit[[20]] <- fit[[20]][[1]]
        names(fit[[20]]) <- c("0", "1")
    }
    if (sum(diag(fit$covb)) > 0) 
        fit$rb <- cov2cor(fit$covb)
    if (sum(diag(fit$covg)) > 0) 
        fit$rg <- cov2cor(fit$covg)
    if (sum(diag(fit$cove)) > 0) 
        fit$re <- cov2cor(fit$cove)
    return(fit)
}
