function (Glist = NULL, stat = NULL, sets = NULL, models = NULL, 
    rsids = NULL, ids = NULL, mask = NULL, lambda = NULL, vb = NULL, 
    vg = NULL, ve = NULL, vy = NULL, pi = NULL, gamma = NULL, 
    mc = 5000, h2 = 0.5, nub = 4, nug = 4, nue = 4, ssb_prior = NULL, 
    ssg_prior = NULL, sse_prior = NULL, vb_prior = NULL, vg_prior = NULL, 
    ve_prior = NULL, updateB = TRUE, updateG = TRUE, updateE = TRUE, 
    updatePi = TRUE, updateMH = FALSE, formatLD = "dense", checkLD = FALSE, 
    shrinkLD = FALSE, shrinkCor = FALSE, pruneLD = FALSE, checkConvergence = FALSE, 
    critVe = 3, critVg = 3, critVb = 3, critPi = 3, critB = 3, 
    critB1 = 0.5, critB2 = 3, verbose = FALSE, eigen_threshold = 0.995, 
    cs_threshold = 0.9, cs_r2 = 0.5, cs_method = "CS2", mh_p = 0.05, 
    mh_r2 = 0.95, nit = 1000, nburn = 100, nthin = 1, nrun = 1, 
    output = "summary", method = "bayesR", algorithm = "mcmc", 
    seed = 10) 
{
    methods <- c("blup", "bayesN", "bayesA", "bayesL", "bayesC", 
        "bayesR")
    method <- match(method, methods) - 1
    if (!sum(method %in% c(0:5)) == 1) 
        stop("method argument specified not valid")
    algorithms <- c("mcmc", "em-mcmc", "mcmc-eigen")
    algorithm <- match(algorithm, algorithms)
    if (is.na(algorithm)) 
        stop("algorithm argument specified not valid")
    if (!is.data.frame(stat)) 
        stop("Stat is not a data frame")
    if (any(sapply(stat[, -c(1:5)], function(x) {
        any(!is.finite(x))
    }))) 
        stop("Some elements in stat not finite")
    if (any(sapply(stat[, -c(1:5)], function(x) {
        any(is.na(x))
    }))) 
        stop("Some elements in stat NA")
    if (is.null(sets)) 
        stop("Please provide sets")
    m <- sum(stat$rsids %in% unlist(Glist$rsids))
    yy <- median((stat$b^2 + (stat$n - 2) * stat$seb^2) * stat$n)
    n <- median(stat$n)
    if (is.null(stat[["ww"]])) 
        stat$ww <- (yy/n)/(stat$seb^2 + stat$b^2/stat$n)
    if (is.null(stat[["wy"]])) 
        stat$wy <- stat$b * stat$ww
    if (method == 4 && is.null(pi)) 
        pi <- c(1 - 0.01, 0.01)
    if (method == 5 && is.null(pi)) 
        pi <- c(0.99, 0.006, 0.003, 0.001)
    if (method == 4 && is.null(gamma)) 
        gamma <- c(0, 1)
    if (method == 5 && is.null(gamma)) 
        gamma <- c(0, 0.01, 0.1, 1)
    vy <- yy/(n - 1)
    if (is.null(ve)) 
        ve <- vy * (1 - h2)
    if (is.null(vg)) 
        vg <- vy * h2
    if (method >= 4 && is.null(vb)) 
        vb <- vg/(mc * sum(pi * gamma))
    if (method >= 4 && is.null(ssb_prior)) 
        ssb_prior <- ((nub - 2)/nub) * (vg/(mc * sum(pi * gamma)))
    ssg_prior <- ((nug - 2)/nug) * vg
    sse_prior <- ((nue - 2)/nue) * ve
    sets <- mapSets(sets = sets, rsids = stat$rsids, index = FALSE)
    if (any(sapply(sets, function(x) {
        any(is.na(x))
    }))) 
        stop("NAs in sets detected - please remove these")
    chr <- as.numeric(unlist(Glist$chr))
    chrSets <- sapply(mapSets(sets = sets, Glist = Glist, index = TRUE), 
        function(x) unique(chr[x]))
    if (length(Glist$bedfiles) == 1) 
        chrSets <- setNames(rep(1, length(sets)), names(sets))
    lsets <- sapply(chrSets, length)
    sets <- sets[lsets == 1]
    if (any(lsets > 1)) 
        stop(paste("Following marker sets mapped to multiple chromosome:", 
            paste(which(lsets > 1), collapse = ",")))
    if (any(lsets == 0)) 
        stop(paste("Following marker sets not mapped to any chromosome:", 
            paste(which(lsets == 0), collapse = ",")))
    bm <- dm <- vector(mode = "list", length = length(sets))
    ves <- vgs <- vbs <- pis <- conv <- vector(mode = "list", 
        length = length(sets))
    bs <- ds <- prob <- vector(mode = "list", length = length(sets))
    pim <- vector(mode = "list", length = length(sets))
    logcpo <- rep(0, length(sets))
    fdr <- csets <- vector(mode = "list", length = length(sets))
    names(bm) <- names(dm) <- names(pim) <- names(sets)
    names(ves) <- names(vgs) <- names(pis) <- names(vbs) <- names(conv) <- names(sets)
    names(bs) <- names(ds) <- names(prob) <- names(sets)
    names(logcpo) <- names(fdr) <- names(csets) <- names(sets)
    attempts <- rep(1, length = length(sets))
    if (is.null(ids)) 
        ids <- Glist$idsLD
    if (is.null(ids)) 
        ids <- Glist$ids
    for (i in 1:length(sets)) {
        chr <- chrSets[[i]]
        rsids <- sets[[i]]
        rws <- match(rsids, stat$rsids)
        message(paste("Processing region:", i))
        pos <- getPos(Glist = Glist, chr = chr, rsids = rsids)
        message(paste("Region size in Mb:", round((max(pos) - 
            min(pos))/1e+06, 2)))
        if (!is.null(Glist$map)) 
            map <- getMap(Glist = Glist, chr = chr, rsids = rsids)
        if (!is.null(Glist$map)) 
            message(paste("Region size in cM:", round(max(map) - 
                min(map), 2)))
        W <- getG(Glist = Glist, chr = chr, rsids = rsids, ids = ids, 
            scale = TRUE)
        B <- crossprod(scale(W))/(nrow(W) - 1)
        if (shrinkLD) 
            B <- corpcor::cor.shrink(W)
        if (algorithm == 3) 
            eig <- eigen(B, symmetric = TRUE)
        for (j in 1:length(eigen_threshold)) {
            if (algorithm == 1) {
                LD <- NULL
                LDvalues <- as.list(as.data.frame(B))
                LDindices <- lapply(1:ncol(B), function(x) {
                  (1:nrow(B)) - 1
                })
                names(LDvalues) <- names(LDindices) <- rsids
            }
            if (algorithm == 3) {
                ww <- stat[rws, "n"]
                scaleb <- sqrt(1/(stat[rws, "n"] * stat[rws, 
                  "seb"] + stat[rws, "b"]^2))
                keep <- cumsum(eig$values)/sum(eig$values) < 
                  eigen_threshold[j]
                z <- t(eig$vectors[, keep]) %*% (stat[rws, "b"] * 
                  scaleb)
                wy <- z/sqrt(eig$values[keep])
                Q <- diag(sqrt(eig$values[keep])) %*% t(eig$vectors[, 
                  keep])
                colnames(Q) <- colnames(B)
                LD <- NULL
                LDvalues <- as.list(as.data.frame(Q))
                LDindices <- lapply(1:ncol(Q), function(x) {
                  (1:nrow(Q)) - 1
                })
                names(LDvalues) <- names(LDindices) <- rsids
            }
            n <- mean(stat[rws, "n"])
            m <- length(rsids)
            b <- rep(0, m)
            mask <- rep(FALSE, m)
            lambda <- rep(ve/vb, m)
            for (run in 1:nrun) {
                if (run > 1) {
                  seed[run] <- sample.int(.Machine$integer.max, 
                    1)
                }
                if (algorithm == 1) {
                  fit <- .Call("_qgg_sbayes_reg", yy = yy, wy = stat$wy[rws], 
                    ww = stat$ww[rws], LDvalues = LDvalues, LDindices = LDindices, 
                    b = b, lambda = lambda, mask = mask, pi = pi, 
                    gamma = gamma, vb = vb, vg = vg, ve = ve, 
                    ssb_prior = ssb_prior, ssg_prior = ssg_prior, 
                    sse_prior = sse_prior, nub = nub, nug = nug, 
                    nue = nue, mh_p = mh_p, mh_r2 = mh_r2, updateB = updateB, 
                    updateG = updateG, updateE = updateE, updatePi = updatePi, 
                    updateMH = updateMH, n = n, nit = nit, nburn = nburn, 
                    nthin = nthin, method = as.integer(method), 
                    algo = as.integer(algorithm), seed = as.integer(seed[run]))
                }
                if (algorithm == 3) {
                  fit <- .Call("_qgg_sbayes_reg_eigen", wy = wy, 
                    ww = ww, LDvalues = LDvalues, LDindices = LDindices, 
                    b = b, lambda = lambda, mask = mask, pi = pi, 
                    gamma = gamma, vb = vb, vg = vg, ve = ve, 
                    ssb_prior = ssb_prior, ssg_prior = ssg_prior, 
                    sse_prior = sse_prior, nub = nub, nug = nug, 
                    nue = nue, mh_p = mh_p, mh_r2 = mh_r2, updateB = updateB, 
                    updateG = updateG, updateE = updateE, updatePi = updatePi, 
                    updateMH = updateMH, n = n, nit = nit, nburn = nburn, 
                    nthin = nthin, method = as.integer(method), 
                    algo = as.integer(algorithm), seed = as.integer(seed[run]))
                }
                names(fit) <- c("bm", "dm", "coef", "vbs", "vgs", 
                  "ves", "pis", "pim", "r", "b", "param", "bs", 
                  "ds", "prob")
                if (nrun > 1) {
                  if (run == 1) {
                    rbm <- rep(0, length(fit$bm))
                    rdm <- rep(0, length(fit$dm))
                    rvbs <- rep(0, length(fit$vbs))
                    rvgs <- rep(0, length(fit$vgs))
                    rves <- rep(0, length(fit$ves))
                    rpis <- rep(0, length(fit$pis))
                    rpim <- rep(0, length(fit$pim))
                    rbs <- rep(0, length(fit$bs))
                    rds <- rep(0, length(fit$ds))
                    rprob <- rep(0, length(fit$prob))
                  }
                  rbm <- rbm + fit$bm/nrun
                  rdm <- rdm + fit$dm/nrun
                  rvbs <- rvbs + fit$vbs/nrun
                  rvgs <- rvgs + fit$vgs/nrun
                  rves <- rves + fit$ves/nrun
                  rpis <- rpis + fit$pis/nrun
                  rpim <- rpim + fit$pim/nrun
                  rbs <- rbs + fit$bs/nrun
                  rds <- rds + fit$ds/nrun
                  rprob <- rprob + fit$prob/nrun
                  if (run == nrun) {
                    fit$bm <- rbm
                    fit$dm <- rdm
                    fit$vbs <- rvbs
                    fit$vgs <- rvgs
                    fit$ves <- rves
                    fit$pis <- rpis
                    fit$pim <- rpim
                    fit$bs <- rbs
                    fit$ds <- rds
                    fit$prob <- rprob
                  }
                }
            }
            if (algorithm == 3) 
                fit$bm <- fit$bm/scaleb
            names(fit$bm) <- names(fit$dm) <- names(fit$b) <- names(LDvalues)
            fit$bs <- matrix(fit$bs, nrow = length(fit$bm))
            fit$ds <- matrix(fit$ds, nrow = length(fit$bm))
            fit$prob <- matrix(fit$prob, nrow = length(fit$bm))
            rownames(fit$bs) <- rownames(fit$ds) <- rownames(fit$prob) <- names(LDvalues)
            colnames(fit$bs) <- colnames(fit$ds) <- colnames(fit$prob) <- 1:(nit + 
                nburn)
            if (algorithm == 3) {
                for (k in 1:nrow(fit$bs)) {
                  fit$bs[k, ] <- fit$bs[k, ]/scaleb[k]
                }
            }
            critve <- critvg <- critvb <- critpi <- critb <- FALSE
            if (!updateE) 
                critve <- TRUE
            if (!updateG) 
                critvg <- TRUE
            if (!updateB) 
                critvb <- TRUE
            if (!updatePi) 
                critpi <- TRUE
            zve <- coda::geweke.diag(fit$ves[nburn:(nburn + nit)])$z
            zvg <- coda::geweke.diag(fit$vgs[nburn:(nburn + nit)])$z
            zvb <- coda::geweke.diag(fit$vbs[nburn:(nburn + nit)])$z
            zpi <- coda::geweke.diag(fit$pis[nburn:(nburn + nit)])$z
            zb <- coda::geweke.diag(apply(fit$bs[, nburn:(nburn + 
                nit)], 2, var))$z
            if (!is.na(zve)) 
                critve <- abs(zve) < critVe
            if (!is.na(zvg)) 
                critvg <- abs(zvg) < critVg
            if (!is.na(zvb)) 
                critvb <- abs(zvb) < critVb
            if (!is.na(zpi)) 
                critpi <- abs(zpi) < critPi
            if (!is.na(zb)) 
                critb <- abs(zb) < critB
            critb1 <- critb2 <- FALSE
            brws <- fit$dm > 0
            if (sum(brws) > 1 && all(c(critve, critvg, critvb, 
                critpi, critb))) {
                tstat <- fit$bs[brws, ]/stat[rws, "b"][brws]
                pdiv <- apply(tstat, 1, function(x) {
                  sum(x[nburn:length(x)] > -critB1 & x[nburn:length(x)] <= 
                    1 + critB1)
                })
                pdiv <- pdiv/length(nburn:ncol(tstat))
                pdiv <- pdiv[is.finite(pdiv)]
                if (any(pdiv < 0.95) && verbose) 
                  plot(pdiv)
                critb1 <- !any(pdiv < 0.95)
                if (!critb1) 
                  message(paste("Convergence not reached for critB1 "))
                critb <- critb1
            }
            if (checkLD && sum(brws) > 1 && !all(c(critve, critvg, 
                critvb, critpi, critb))) {
                bout <- checkb(B = B[brws, brws], b = stat[rws, 
                  "b"][brws], seb = stat[rws, "seb"][brws], critB = critB2, 
                  verbose = verbose)
                critb2 <- !any(bout$outliers)
                if (critb2) 
                  message(paste("Convergence not reached for critB2 "))
                critb <- critb2
            }
            converged <- critve & critvg & critvb & critpi & 
                critb
            if (!checkConvergence) 
                converged <- TRUE
            if (verbose) {
                layout(matrix(1:4, ncol = 2))
                pipsets <- splitWithOverlap(1:length(rsids), 
                  100, 99)
                pip <- fit$dm
                plot(pip, ylim = c(0, max(pip)), ylab = "PIP", 
                  xlab = "Position", frame.plot = FALSE)
                plot(-log10(stat[rws, "p"]), ylab = "-log10(P)", 
                  xlab = "Position", frame.plot = FALSE)
                hist(fit$ves, main = "Ve", xlab = "")
                plot(y = fit$bm, x = stat[rws, "b"], ylab = "Adjusted", 
                  xlab = "Marginal", frame.plot = FALSE)
                abline(h = 0, v = 0, lwd = 2, col = 2, lty = 2)
            }
            attempts[i] <- j
            if (!converged) {
                message(paste("Convergence not reached using eigen_threshold:", 
                  eigen_threshold[j]))
                criteria_names <- c("Variance of errors (critve)", 
                  "Genetic variance (critvg)", "Marker variance (critvb)", 
                  "Inclusion probability (critpi)", "Posterior mean (critb)")
                criteria_status <- c(critve, critvg, critvb, 
                  critpi, critb)
                message("Convergence criteria:")
                for (k in seq_along(criteria_names)) {
                  message(sprintf("  %s: %s", criteria_names[k], 
                    ifelse(criteria_status[k], "Met", "Not Met")))
                }
                if (checkConvergence & algorithm == 1) {
                  message("Serious convergence issues detected. Please check your summary statistics and LD reference data.")
                  message("If using in-sample LD and summary statistics, consider increasing the default values for the following parameters: critVe=3, critVg=3, critVb=3, critPi=3, critB=3.")
                  message("If using external LD and summary statistics, consider using algorithm=\"mcmc-eigen\" as a more robust method.")
                  message("Program terminated.")
                  return(NULL)
                }
            }
            if (converged) {
                if (verbose & algorithm == 1) 
                  message("Convergence reached")
                if (verbose & algorithm == 3) 
                  message(paste("Convergence reached using eigen_threshold:", 
                    eigen_threshold[j]))
                break
            }
        }
        if (verbose) 
            message("Computing FDR")
        cutoffs <- seq(0.01, 0.99, by = 0.01)
        cutoff_indices <- lapply(cutoffs, function(cutoff) fit$dm > 
            cutoff)
        bfdrs <- sapply(cutoff_indices, function(rws) {
            if (any(rws)) {
                fdrs <- rowMeans(1 - fit$prob[rws, , drop = FALSE], 
                  na.rm = TRUE)
                c(mean = mean(fdrs, na.rm = TRUE), quantile(fdrs, 
                  c(0.025, 0.975), na.rm = TRUE))
            }
            else {
                c(NA, NA, NA)
            }
        })
        bfdrs <- t(bfdrs)
        rownames(bfdrs) <- round(cutoffs, 2)
        if (verbose) 
            message("Saving results")
        bm[[i]] <- fit$bm
        dm[[i]] <- fit$dm
        pim[[i]] <- fit$pim
        ves[[i]] <- fit$ves
        vbs[[i]] <- fit$vbs
        vgs[[i]] <- fit$vgs
        pis[[i]] <- fit$pis
        conv[[i]] <- c(zve, zvg, zvb, zpi, zb)
        if (output == "full") {
            bs[[i]] <- fit$bs
            ds[[i]] <- fit$ds
            prob[[i]] <- fit$prob
        }
        fdr[[i]] <- bfdrs
        logcpo[i] <- fit$param[4]
        if (verbose) 
            message("Compute credible sets")
        if (sum(fit$dm) > cs_threshold) 
            csets[[i]] <- crs(prob = fit$dm, B = B, threshold = cs_threshold, 
                r2 = cs_r2, method = cs_method)
        names(bm[[i]]) <- names(dm[[i]]) <- rsids
    }
    fit <- NULL
    fit$bm <- bm
    fit$dm <- dm
    fit$pim <- pim
    fit$ves <- ves
    fit$vbs <- vbs
    fit$vgs <- vgs
    fit$pis <- pis
    if (output == "full") {
        fit$bs <- bs
        fit$ds <- ds
        fit$prob <- prob
    }
    fit$fdr <- fdr
    fit$logcpo <- logcpo
    fit$cs <- csets
    pip <- sapply(fit$dm, sum)
    minb <- sapply(fit$bm, min)
    maxb <- sapply(fit$bm, max)
    m <- sapply(fit$bm, length)
    bm <- unlist(unname(fit$bm))
    dm <- unlist(unname(fit$dm))
    marker <- data.frame(rsids = unlist(Glist$rsids), chr = unlist(Glist$chr), 
        pos = unlist(Glist$pos), ea = unlist(Glist$a1), nea = unlist(Glist$a2), 
        eaf = unlist(Glist$af), stringsAsFactors = FALSE)
    marker <- marker[marker$rsids %in% names(bm), ]
    fit$stat <- data.frame(marker, bm = bm[marker$rsids], dm = dm[marker$rsids], 
        stringsAsFactors = FALSE)
    fit$stat$vm <- 2 * (1 - fit$stat$eaf) * fit$stat$eaf * fit$stat$bm^2
    fit$method <- methods[method + 1]
    fit$mask <- mask
    ve <- sapply(fit$ves, function(x) {
        mean(x[nburn:length(x)])
    })
    vg <- sapply(fit$vgs, function(x) {
        mean(x[nburn:length(x)])
    })
    vb <- sapply(fit$vbs, function(x) {
        mean(x[nburn:length(x)])
    })
    pi <- sapply(fit$pis, function(x) {
        mean(x[nburn:length(x)])
    })
    ve_ci <- t(sapply(fit$ves, function(x) {
        quantile(x[nburn:length(x)], c(0.025, 0.975))
    }))
    vg_ci <- t(sapply(fit$vgs, function(x) {
        quantile(x[nburn:length(x)], c(0.025, 0.975))
    }))
    vb_ci <- t(sapply(fit$vbs, function(x) {
        quantile(x[nburn:length(x)], c(0.025, 0.975))
    }))
    pi_ci <- t(sapply(fit$pis, function(x) {
        quantile(x[nburn:length(x)], c(0.025, 0.975))
    }))
    fit$ci <- list(ve = cbind(mean = ve, ve_ci), vg = cbind(mean = vg, 
        vg_ci), vb = cbind(mean = vb, vb_ci), pi = cbind(mean = pi, 
        pi_ci))
    if (!is.null(Glist$map)) 
        map <- unlist(Glist$map)
    pos <- unlist(Glist$pos)
    sets <- lapply(fit$bm, names)
    setsindex <- mapSets(sets = sets, rsids = unlist(Glist$rsids))
    if (!is.null(Glist$map)) 
        cm <- sapply(setsindex, function(x) {
            max(map[x]) - min(map[x])
        })
    mb <- sapply(setsindex, function(x) {
        (max(pos[x]) - min(pos[x]))/1e+06
    })
    minmb <- sapply(setsindex, function(x) {
        min(pos[x])
    })
    maxmb <- sapply(setsindex, function(x) {
        max(pos[x])
    })
    chr <- unlist(Glist$chr)
    chr <- sapply(setsindex, function(x) {
        as.numeric(unique(chr[x]))[1]
    })
    b <- stat[fit$stat$rsids, "b"]
    conv <- t(as.data.frame(conv))
    colnames(conv) <- c("zve", "zvg", "zvb", "zpi", "zb")
    fit$conv <- data.frame(conv, ntrials = attempts, cutoff = eigen_threshold[attempts])
    if (is.null(Glist$map)) 
        fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi, 
            pip = pip, minb = minb, maxb = maxb, m = m, mb = mb, 
            chr = chr, minmb = minmb, maxmb = maxmb)
    if (!is.null(Glist$map)) 
        fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi, 
            pip = pip, minb = minb, maxb = maxb, m = m, mb = mb, 
            cm = cm, chr = chr, minmb = minmb, maxmb = maxmb)
    rownames(fit$conv) <- rownames(fit$post) <- names(sets)
    fit$ve <- mean(ve)
    fit$vg <- sum(vg)
    fit$b <- b
    fit$seed <- seed
    return(fit)
}
