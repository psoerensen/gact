function (Glist = NULL, chr = NULL, r2 = 0, onebased = FALSE, 
    rsids = NULL, format = "sparse") 
{
    msize <- Glist$msize
    rsidsLD <- Glist$rsidsLD[[chr]]
    mchr <- length(rsidsLD)
    indices <- 1:mchr
    if (!onebased) 
        indices <- indices - 1
    rsids_indices <- c(rep(NA, msize), indices, rep(NA, msize))
    ld_indices <- vector(length = mchr, mode = "list")
    ld_values <- vector(length = mchr, mode = "list")
    names(ld_indices) <- names(ld_values) <- rsidsLD
    if (is.null(Glist$ldfiles[chr])) 
        stop("LD file not found")
    fnLD <- Glist$ldfiles[chr]
    bfLD <- file(fnLD, "rb")
    nld <- as.integer(msize * 2 + 1)
    for (i in 1:mchr) {
        ld <- readBin(bfLD, "numeric", n = nld, size = 4, endian = "little")
        ld[msize + 1] <- 1
        cls <- which((ld^2) > r2) + i - 1
        ld_indices[[i]] <- rsids_indices[cls]
        ld_values[[i]] <- ld[(ld^2) > r2]
    }
    close(bfLD)
    return(list(indices = ld_indices, values = ld_values, rsids = names(ld_indices), 
        onebased = onebased, msize = msize, fnLD = fnLD))
}
