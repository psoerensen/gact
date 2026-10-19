function (y = NULL, X = NULL, W = NULL, stat = NULL, covs = NULL, 
    trait = NULL, fit = NULL, Glist = NULL, chr = NULL, rsids = NULL, 
    b = NULL, bm = NULL, seb = NULL, LD = NULL, n = NULL, formatLD = "dense", 
    vg = NULL, vb = NULL, ve = NULL, ssg_prior = NULL, ssb_prior = NULL, 
    sse_prior = NULL, lambda = NULL, scaleY = TRUE, h2 = NULL, 
    pi = 0.001, gamma = 1, updateB = TRUE, updateG = TRUE, updateE = TRUE, 
    updatePi = TRUE, adjustE = TRUE, models = NULL, nug = 4, 
    nub = 4, nue = 4, verbose = FALSE, msize = 100, mask = NULL, 
    checkConvergence = FALSE, GRMlist = NULL, ve_prior = NULL, 
    vg_prior = NULL, tol = 0.001, nit = 100, nburn = 0, nthin = 1, 
    nit_local = NULL, nit_global = NULL, method = "mixed", algorithm = "mcmc") 
{
    mask.rsids <- NULL
    if (!is.null(mask)) 
        mask.rsids <- unique(as.vector(apply(mask, 2, function(x) {
            as.vector(rownames(mask)[x])
        })))
    methods <- c("blup", "bayesN", "bayesA", "bayesL", "bayesC", 
        "bayesR")
    method <- match(method, methods) - 1
    if (!sum(method %in% c(0:5)) == 1) 
        stop("Method specified not valid")
    algorithms <- c("mcmc", "em-mcmc")
    algorithm <- match(algorithm, algorithms)
    if (is.na(algorithm)) 
        stop("algorithm argument specified not valid")
    if (method == 0) {
        updateB = FALSE
        updateE = FALSE
    }
    nt <- 1
    if (!is.null(y)) {
        if (is.list(y)) 
            nt <- length(y)
        if (is.matrix(y)) 
            nt <- ncol(y)
    }
    if (!is.null(stat)) {
        if (!is.data.frame(stat) && is.list(stat)) 
            nt <- ncol(stat$b)
        if (any(sapply(stat[, -c(1:5)], function(x) {
            any(!is.finite(x))
        }))) 
            stop("Some elements in stat not finite")
        if (any(sapply(stat[, -c(1:5)], function(x) {
            any(is.na(x))
        }))) 
            stop("Some elements in stat are NA's")
    }
    if (!is.null(GRMlist)) 
        analysis <- "mtmc-mixed"
    if (nt == 1 && !is.null(y) && !is.null(W) && formatLD == 
        "dense") 
        analysis <- "st-blr-individual-level-default"
    if (nt == 1 && !is.null(y) && !is.null(W) && formatLD == 
        "sparse") 
        analysis <- "st-blr-individual-level-sbayes"
    if (nt == 1 && !is.null(y) && formatLD == "sparse") 
        analysis <- "st-blr-individual-level-sparse-ld"
    if (nt == 1 && !is.null(y) && formatLD == "dense") 
        analysis <- "st-blr-individual-level-dense-ld"
    if (nt == 1 && is.null(y) && !is.null(stat) && !is.null(Glist)) 
        analysis <- "st-blr-sumstat-sparse-ld"
    if (nt > 1 && !is.null(y) && !is.null(W)) 
        analysis <- "mt-blr-individual-level-dense-ld"
    if (nt > 1 && !is.null(y) && formatLD == "sparse") 
        analysis <- "mt-blr-individual-level-sparse-ld"
    if (nt > 1 && !is.null(stat) && !is.null(Glist)) 
        analysis <- "mt-blr-sumstat-sparse-ld"
    message(paste("Type of analysis performed:", analysis))
    if (!is.null(GRMlist)) {
        fit <- bmm(y = y, X = X, W = W, GRMlist = GRMlist, vg = vg, 
            ve = ve, nug = nug, nue = nue, vg_prior = vg_prior, 
            ve_prior = ve_prior, updateG = updateB, updateE = updateE, 
            nit = nit, nburn = nburn, tol = tol, verbose = verbose)
    }
    if (nt == 1 && !is.null(y) && !is.null(W) && formatLD == 
        "dense") {
        fit <- bayes(y = y, X = X, W = W, b = b, scaled = TRUE, 
            h2 = h2, pi = pi, gamma = gamma, lambda = lambda, 
            vg = vg, vb = vb, ve = ve, nub = nub, nug = nug, 
            nue = nue, ssb_prior = ssb_prior, ssg_prior = ssg_prior, 
            sse_prior = sse_prior, updateB = updateB, updateG = updateG, 
            updateE = updateE, updatePi = updatePi, nit = nit, 
            nburn = nburn, nthin = nthin, method = method)
    }
    if (nt > 1 && !is.null(y) && !is.null(W)) {
        fit <- mtbayes(y = y, X = X, W = W, b = b, bm = bm, seb = seb, 
            LD = LD, n = n, vg = vg, vb = vb, ve = ve, ssb_prior = ssb_prior, 
            sse_prior = sse_prior, lambda = lambda, scaleY = scaleY, 
            h2 = h2, pi = pi, updateB = updateB, updateE = updateE, 
            updatePi = updatePi, models = models, nub = nub, 
            nue = nue, nit = nit, method = method, formatLD = formatLD, 
            algorithm = algorithm, verbose = verbose)
    }
    if (analysis == "st-blr-sumstat-sparse-ld") {
        if (is.data.frame(stat)) {
            nt <- 1
            rsidsLD <- unlist(Glist$rsidsLD)
            m <- length(rsidsLD)
            b <- wy <- ww <- matrix(0, nrow = length(rsidsLD), 
                ncol = nt)
            mask <- matrix(TRUE, nrow = length(rsidsLD), ncol = nt)
            rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- rsidsLD
            trait_names <- "bm"
            stat <- stat[rownames(stat) %in% rsidsLD, ]
            if (is.null(stat$ww)) 
                stat$ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
            if (is.null(stat$wy)) 
                stat$wy <- stat$b * stat$ww
            if (!is.null(stat$n)) 
                n <- as.integer(median(stat$n))
            ww[rownames(stat), 1] <- stat$ww
            wy[rownames(stat), 1] <- stat$wy
            mask[rownames(stat), 1] <- FALSE
            if (!is.null(mask.rsids)) 
                mask[rownames(mask) %in% mask.rsids, 1] <- TRUE
            if (any(is.na(wy))) 
                stop("Missing values in wy")
            if (any(is.na(ww))) 
                stop("Missing values in ww")
            b2 <- stat$b^2
            seb2 <- stat$seb^2
            yy <- (b2 + (n - 2) * seb2) * stat$ww
            yy <- median(yy)
        }
        if (!is.data.frame(stat) && is.list(stat)) {
            nt <- ncol(stat$b)
            trait_names <- colnames(stat$b)
            if (is.null(trait_names)) 
                trait_names <- paste0("T", 1:nt)
            rsidsLD <- unlist(Glist$rsidsLD)
            m <- length(rsidsLD)
            b <- wy <- ww <- matrix(0, nrow = length(rsidsLD), 
                ncol = nt)
            mask <- matrix(TRUE, nrow = length(rsidsLD), ncol = nt)
            rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- rsidsLD
            colnames(b) <- colnames(wy) <- colnames(ww) <- colnames(mask) <- trait_names
            rws <- rownames(stat$b) %in% rsidsLD
            if (is.null(stat$ww)) 
                stat$ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
            if (is.null(stat$wy)) 
                stat$wy <- stat$b * stat$ww
            if (!is.null(stat$n)) 
                n <- as.integer(apply(stat$n[rws, ], 2, median))
            ww[rownames(stat$ww[rws, ]), ] <- stat$ww[rws, ]
            wy[rownames(stat$wy[rws, ]), ] <- stat$wy[rws, ]
            mask[rownames(stat$wy[rws, ]), ] <- FALSE
            if (!is.null(mask.rsids)) 
                mask[rownames(mask) %in% mask.rsids, ] <- TRUE
            if (any(is.na(wy))) 
                stop("Missing values in wy")
            if (any(is.na(ww))) 
                stop("Missing values in ww")
            b2 <- (stat$b[rws, ])^2
            seb2 <- (stat$seb[rws, ])^2
            for (i in 1:nt) {
                seb2[, i] <- (n[i] - 2) * seb2[, i]
            }
            yy <- (b2 + seb2) * stat$ww[rws, ]
            yy <- apply(yy, 2, median)
        }
        if (is.null(chr)) 
            chromosomes <- 1:Glist$nchr
        if (!is.null(chr)) 
            chromosomes <- chr
        if (is.null(LD)) 
            LD <- vector(length = Glist$nchr, mode = "list")
        bm <- dm <- fit <- res <- vector(length = Glist$nchr, 
            mode = "list")
        names(bm) <- names(dm) <- names(fit) <- names(res) <- 1:Glist$nchr
        for (chr in chromosomes) {
            if (is.null(LD[[chr]])) {
                if (verbose) 
                  print(paste("Extract sparse LD matrix for chromosome:", 
                    chr))
                LD[[chr]] <- getSparseLD(Glist = Glist, chr = chr, 
                  onebased = FALSE)
            }
            rsidsLD <- names(LD[[chr]]$values)
            clsLD <- match(rsidsLD, Glist$rsids[[chr]])
            bmchr <- dmchr <- NULL
            for (trait in 1:nt) {
                if (verbose) 
                  print(paste("Fit", methods[method + 1], "on chromosome:", 
                    chr, "for trait", trait))
                fit[[chr]] <- sbayes_sparse(yy = yy[trait], wy = wy[rsidsLD, 
                  trait], ww = ww[rsidsLD, trait], b = b[rsidsLD, 
                  trait], LDvalues = LD[[chr]]$values, LDindices = LD[[chr]]$indices, 
                  mask = mask[rsidsLD, trait], method = method, 
                  nit = nit, nburn = nburn, nthin = nthin, n = n[trait], 
                  m = m, pi = pi, nue = nue, nub = nub, h2 = h2[trait], 
                  lambda = lambda, vb = vb, ve = ve, ssb_prior = ssb_prior, 
                  sse_prior = sse_prior, updateB = updateB, updateE = updateE, 
                  updatePi = updatePi, updateG = updateG, adjustE = adjustE, 
                  algorithm = algorithm)
                bmchr <- cbind(bmchr, fit[[chr]]$bm)
                dmchr <- cbind(dmchr, fit[[chr]]$dm)
            }
            colnames(bmchr) <- trait_names
            res[[chr]] <- data.frame(rsids = rsidsLD, chr = rep(chr, 
                length(rsidsLD)), pos = Glist$pos[[chr]][clsLD], 
                ea = Glist$a1[[chr]][clsLD], nea = Glist$a2[[chr]][clsLD], 
                eaf = Glist$af[[chr]][clsLD], bm = bmchr, dm = dmchr, 
                stringsAsFactors = FALSE)
            rownames(res[[chr]]) <- rsidsLD
            LD[[chr]]$values <- NULL
            LD[[chr]]$indices <- NULL
        }
        res <- do.call(rbind, res)
        rownames(res) <- res$rsids
        fit$stat <- res
        fit$stat$vm <- 2 * (1 - fit$stat$eaf) * fit$stat$eaf * 
            fit$stat$bm^2
        fit$method <- methods[method + 1]
        fit$ves <- lapply(fit[chromosomes], function(x) {
            x$ves
        })
        fit$vgs <- lapply(fit[chromosomes], function(x) {
            x$vgs
        })
        fit$vbs <- lapply(fit[chromosomes], function(x) {
            x$vbs
        })
        fit$ssbs <- lapply(fit[chromosomes], function(x) {
            x$ssbs
        })
        fit$pis <- lapply(fit[chromosomes], function(x) {
            x$pis
        })
        fit$pim <- lapply(fit[chromosomes], function(x) {
            x$pim
        })
        fit$param <- lapply(fit[chromosomes], function(x) {
            x$param
        })
        fit$bs <- lapply(fit[chromosomes], function(x) {
            matrix(x$bs, ncol = length(x$b), byrow = TRUE)
        })
        fit$mask <- mask
        zve <- zvg <- zvb <- zpi <- NULL
        if (checkConvergence) {
            zve <- sapply(fit$ves[chromosomes], function(x) {
                coda::geweke.diag(x[nburn:length(x)])$z
            })
            zvg <- sapply(fit$vgs[chromosomes], function(x) {
                coda::geweke.diag(x[nburn:length(x)])$z
            })
            zvb <- sapply(fit$vbs[chromosomes], function(x) {
                coda::geweke.diag(x[nburn:length(x)])$z
            })
            zpi <- sapply(fit$pis[chromosomes], function(x) {
                coda::geweke.diag(x[nburn:length(x)])$z
            })
        }
        ve <- sapply(fit$ves[chromosomes], function(x) {
            mean(x[nburn:length(x)])
        })
        vg <- sapply(fit$vgs[chromosomes], function(x) {
            mean(x[nburn:length(x)])
        })
        vb <- sapply(fit$vbs[chromosomes], function(x) {
            mean(x[nburn:length(x)])
        })
        pi <- sapply(fit$pim[chromosomes], function(x) {
            1 - x[1]
        })
        fit$conv <- data.frame(zve = zve, zvg = zvg, zvb = zvb, 
            zpi = zpi)
        fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi)
        fit$ve <- mean(ve)
        fit$vg <- sum(vg)
    }
    if (analysis == "mt-blr-sumstat-sparse-ld") {
        nt <- ncol(stat$b)
        trait_names <- colnames(stat$b)
        if (is.null(trait_names)) 
            trait_names <- paste0("T", 1:nt)
        rsidsLD <- unlist(Glist$rsidsLD)
        b <- wy <- ww <- matrix(0, nrow = length(rsidsLD), ncol = nt)
        rownames(b) <- rownames(wy) <- rownames(ww) <- rsidsLD
        colnames(b) <- colnames(wy) <- colnames(ww) <- trait_names
        rws <- rownames(stat$b) %in% rsidsLD
        if (!is.null(stat$n)) 
            n <- as.integer(colMeans(stat$n[rws, ]))
        if (!is.null(stat$ww)) 
            ww[rownames(stat$ww[rws, ]), ] <- stat$ww[rws, ]
        if (is.null(stat$ww)) 
            ww[rownames(stat$n[rws, ]), ] <- stat$n[rws, ]
        if (!is.null(stat$wy)) 
            wy[rownames(stat$wy[rws, ]), ] <- stat$wy[rws, ]
        if (is.null(stat$wy)) 
            wy[rownames(stat$b[rws, ]), ] <- stat$b[rws, ] * 
                stat$n[rws, ]
        if (any(is.na(wy))) 
            stop("Missing values in wy")
        if (any(is.na(ww))) 
            stop("Missing values in ww")
        b2 <- (stat$b[rws, ])^2
        seb2 <- (stat$seb[rws, ])^2
        for (i in 1:nt) {
            seb2[, i] <- (n[i] - 2) * seb2[, i]
        }
        yy <- (b2 + seb2) * stat$ww[rws, ]
        yy <- colMeans(yy)
        if (is.null(chr)) 
            chromosomes <- 1:Glist$nchr
        if (!is.null(chr)) 
            chromosomes <- chr
        if (is.null(LD)) 
            LD <- vector(length = Glist$nchr, mode = "list")
        bm <- dm <- fit <- res <- vector(length = Glist$nchr, 
            mode = "list")
        names(bm) <- names(dm) <- names(fit) <- names(res) <- 1:Glist$nchr
        for (chr in chromosomes) {
            if (is.null(LD[[chr]])) {
                if (verbose) 
                  print(paste("Extract sparse LD matrix for chromosome:", 
                    chr))
                LD[[chr]] <- getSparseLD(Glist = Glist, chr = chr, 
                  onebased = FALSE)
            }
            rsidsLD <- names(LD[[chr]]$values)
            clsLD <- match(rsidsLD, Glist$rsids[[chr]])
            bmchr <- NULL
            fit[[chr]] <- mt_sbayes_sparse(yy = yy, ww = ww[rsidsLD, 
                ], wy = wy[rsidsLD, ], b = b[rsidsLD, ], LDvalues = LD[[chr]]$values, 
                LDindices = LD[[chr]]$indices, n = n, nit = nit, 
                pi = pi, nue = nue, nub = nub, h2 = h2, vb = vb, 
                ve = ve, ssb_prior = ssb_prior, sse_prior = sse_prior, 
                updateB = updateB, updateE = updateE, updatePi = updatePi, 
                models = models, method = method, verbose = verbose)
            res[[chr]] <- data.frame(rsids = rsidsLD, chr = rep(chr, 
                length(rsidsLD)), pos = Glist$pos[[chr]][clsLD], 
                ea = Glist$a1[[chr]][clsLD], nea = Glist$a2[[chr]][clsLD], 
                eaf = Glist$af[[chr]][clsLD], bm = fit[[chr]]$bm, 
                dm = fit[[chr]]$dm, stringsAsFactors = FALSE)
            rownames(res[[chr]]) <- rsidsLD
            LD[[chr]]$values <- NULL
            LD[[chr]]$indices <- NULL
        }
        res <- do.call(rbind, res)
        rownames(res) <- res$rsids
        fit$stat <- res
        fit$method <- methods[method + 1]
        fit$stat$vm <- 2 * (1 - fit$stat$eaf) * fit$stat$eaf * 
            fit$stat$bm^2
    }
    return(fit)
}
