function (stat = NULL, sets = NULL, method = "magma", type = "joint", 
    test = "one-sided", pi = 0.001, nit = 5000, nburn = 1000) 
{
    if (is.null(stat) || is.null(sets)) {
        stop("Both 'stat' and 'sets' must be provided.")
    }
    if (is.vector(stat)) 
        stat <- as.matrix(stat)
    if (is.null(rownames(stat))) 
        stop("Please provide names or rownames to stat object")
    y <- scale(stat, center = TRUE, scale = TRUE)
    sets <- mapSets(sets = sets, rsids = rownames(stat), index = FALSE)
    X <- designMatrix(sets = sets, rowids = rownames(y))
    m <- sapply(sets, length)
    if (method == "magma") {
        stat_summary <- computeStat(X = X, y = y[rownames(X), 
            ], scale = TRUE)
        bMAGMA <- solve(stat_summary$XX + diag(0.001, nrow(stat_summary$XX))) %*% 
            stat_summary$Xy
        bMARG <- (1/diag(stat_summary$XX)) * stat_summary$Xy
        sebMAGMA <- sqrt(diag(solve(stat_summary$XX + diag(0.001, 
            nrow(stat_summary$XX)))))
        sebMARG <- sqrt((1/diag(stat_summary$XX)))
        zMAGMA <- bMAGMA/sebMAGMA
        zMARG <- bMARG/sebMARG
        if (test == "two-sided") {
            pMAGMA <- pnorm(abs(zMAGMA), mean = 0, sd = 1, lower.tail = FALSE)
            pMARG <- pnorm(abs(zMARG), mean = 0, sd = 1, lower.tail = FALSE)
        }
        if (test == "one-sided") {
            pMAGMA <- pnorm(zMAGMA, mean = 0, sd = 1, lower.tail = FALSE)
            pMARG <- pnorm(zMARG, mean = 0, sd = 1, lower.tail = FALSE)
        }
        if (type == "marginal") 
            df <- data.frame(ID = names(sets), m = m, b = bMARG, 
                seb = sebMARG, z = zMARG, p = pMARG)
        if (type == "joint") 
            df <- data.frame(ID = names(sets), m = m, b = bMAGMA, 
                seb = sebMAGMA, z = zMAGMA, p = pMAGMA)
        o <- order(df$p, decreasing = FALSE)
        df[, 3:5] <- round(df[, 3:5], 4)
        rownames(df) <- NULL
        return(df[o, ])
    }
    if (method %in% c("blr", "bayesC", "bayesR")) {
        if (method == "blr") 
            method <- "bayesC"
        stat <- computeStat(X = X, y = y[rownames(X), ], scale = TRUE)
        if (ncol(y) == 1) {
            fit <- blr(yy = stat$yy, XX = stat$XX, Xy = stat$Xy, 
                n = stat$n, method = method, pi = pi, nit = nit, 
                nburn = nburn)
            df <- data.frame(ID = names(sets), m = m, b = fit$bm, 
                PIP = fit$dm)
            o <- order(df$PIP, decreasing = TRUE)
            df[, 3:4] <- round(df[, 3:4], 4)
            rownames(df) <- NULL
            return(df[o, ])
        }
        if (ncol(y) > 1) {
            stat$XX <- rep(list(stat$XX), ncol(y))
            fit <- mtblr(yy = stat$yy, XX = stat$XX, Xy = stat$Xy, 
                n = stat$n, method = method, pi = pi, nit = nit, 
                nburn = nburn)
            o <- order(rowSums(fit$dm), decreasing = TRUE)
            df <- list(feature = data.frame(ID = names(sets), 
                m = m)[o, ], b = round(fit$bm[o, ], 4), PIP = round(fit$dm[o, 
                ], 4))
            return(df)
        }
    }
}
