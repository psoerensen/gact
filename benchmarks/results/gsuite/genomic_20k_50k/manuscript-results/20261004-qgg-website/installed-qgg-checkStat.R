function (Glist = NULL, stat = NULL, excludeMAF = 0.01, excludeMAFDIFF = 0.05, 
    excludeINFO = 0.8, excludeCGAT = TRUE, excludeINDEL = TRUE, 
    excludeDUPS = TRUE, excludeMHC = FALSE, excludeMISS = 0.05, 
    excludeHWE = 1e-12) 
{
    stat$ea <- toupper(stat$ea)
    stat$nea <- toupper(stat$nea)
    cpra <- unlist(Glist$cpra)
    rsids <- unlist(Glist$rsids)
    if (!is.data.frame(stat)) 
        stop("stat should be  a data frame")
    if (!is.null(stat$rsids)) 
        rownames(stat) <- stat$rsids
    fm_internal <- c("rsids", "chr", "pos", "ea", "nea", "eaf", 
        "b", "seb", "p", "n", "info")
    fm_external <- c("marker", "chr", "pos", "ea", "nea", "eaf", 
        "b", "seb", "p", "n", "info")
    fm_gbayes <- c("rsids", "chr", "pos", "ea", "nea", "eaf", 
        "bm", "dm", "vm")
    format <- "unknown"
    if (all(fm_internal[1:5] %in% colnames(stat))) 
        format <- "internal"
    if (all(fm_external[1:5] %in% colnames(stat))) 
        format <- "external"
    if (all(fm_gbayes %in% colnames(stat))) {
        stat <- stat[, 1:7]
        colnames(stat) <- c("rsids", "chr", "pos", "ea", "nea", 
            "eaf", "b")
        format <- "internal"
    }
    if (format == "unknown") {
        message("Column headings for stat object not found")
        message("Column headings for stat object should minimal be:")
        print(fm_external[1:5])
        message("or:")
        print(fm_external)
        stop("please revised your stat object accordingly")
    }
    if (format == "internal") 
        stat <- stat[, colnames(stat) %in% fm_internal]
    if (format == "external") 
        stat <- stat[, colnames(stat) %in% fm_external]
    if (format == "external") {
        cpra1 <- paste(stat[, "chr"], stat[, "pos"], stat[, "ea"], 
            stat[, "nea"], sep = "_")
        cpra2 <- paste(stat[, "chr"], stat[, "pos"], stat[, "nea"], 
            stat[, "ea"], sep = "_")
        mapped <- cpra1 %in% cpra | cpra2 %in% cpra
        message("Map markers based on cpra")
        message(paste("Number of markers in stat mapped to marker ids in Glist:", 
            sum(mapped)))
        message(paste("Number of markers in stat not mapped to marker ids in Glist:", 
            sum(!mapped)))
        stat <- stat[mapped, ]
        cpra1 <- cpra1[mapped]
        cpra2 <- cpra2[mapped]
        rws1 <- match(cpra1, cpra)
        rws2 <- match(cpra2, cpra)
        stat$marker[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
        stat$marker[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]
        isdup <- duplicated(stat$marker)
        if (any(isdup)) 
            message("Removing markers with duplicated ids")
        if (any(isdup)) 
            message(paste("Number of markers duplicated in stat:", 
                sum(isdup)))
        stat <- stat[!isdup, ]
        rownames(stat) <- stat$marker
    }
    marker <- data.frame(rsids = unlist(Glist$rsids), cpra = unlist(Glist$cpra), 
        chr = unlist(Glist$chr), pos = unlist(Glist$pos), ea = unlist(Glist$a1), 
        nea = unlist(Glist$a2), eaf = unlist(Glist$af), stringsAsFactors = FALSE)
    message("Filtering markers based on information in Glist:")
    message("")
    if (is.null(Glist$rsidsLD)) {
        rsids <- gfilter(Glist = Glist, excludeMAF = excludeMAF, 
            excludeMISS = excludeMISS, excludeCGAT = excludeCGAT, 
            excludeINDEL = excludeINDEL, excludeDUPS = excludeDUPS, 
            excludeHWE = excludeHWE, excludeMHC = excludeMHC)
        marker <- marker[marker$rsids %in% rsids, ]
        message("")
    }
    if (!is.null(Glist$rsidsLD)) {
        rsids <- unlist(Glist$rsidsLD)
        message(paste("Number of markers in sparse LD matrices:", 
            sum(marker$rsids %in% rsids)))
        message("")
        marker <- marker[marker$rsids %in% rsids, ]
    }
    rownames(marker) <- marker$rsids
    message("Filtering markers based on information in stat:")
    message("")
    if (!is.null(stat$rsids)) 
        marker_in_stat <- marker$rsids %in% stat$rsids
    if (!is.null(stat$marker)) 
        marker_in_stat <- marker$rsids %in% stat$marker
    message(paste("Number of markers in stat also found in bimfiles:", 
        sum(marker_in_stat)))
    message("")
    if (sum(marker_in_stat) == 0) 
        stop("No marker ids found in bimfiles")
    marker <- marker[marker_in_stat, ]
    stat <- stat[marker$rsids, ]
    aligned <- stat$ea == marker$ea
    message(paste("Number of effect alleles aligned with first allele in bimfiles:", 
        sum(aligned)))
    message(paste("Number of effect alleles not aligned with first allele in bimfiles:", 
        sum(!aligned)))
    message("")
    effect <- stat[, "b"]
    effect_allele <- stat[, "ea"]
    non_effect_allele <- stat[, "nea"]
    if (!is.null(stat$eaf)) 
        effect_allele_freq <- stat[, "eaf"]
    stat[!aligned, "b"] <- -effect[!aligned]
    stat[!aligned, "ea"] <- non_effect_allele[!aligned]
    stat[!aligned, "nea"] <- effect_allele[!aligned]
    if (!is.null(stat$eaf)) 
        stat[!aligned, "eaf"] <- 1 - effect_allele_freq[!aligned]
    if (is.null(stat$eaf)) {
        message("No effect allele frequency (eaf) provided - using eaf in Glist")
        stat$eaf <- marker$eaf
    }
    excludeMAFDIFF <- abs(marker$eaf - stat$eaf) > excludeMAFDIFF
    message(paste("Number of markers excluded by large difference between MAF difference:", 
        sum(excludeMAFDIFF)))
    message("")
    stat <- stat[!excludeMAFDIFF, ]
    marker <- marker[!excludeMAFDIFF, ]
    if (colnames(stat)[1] == "marker") 
        colnames(stat)[1] <- "rsids"
    if (is.null(stat$n)) 
        stat$n <- 4/((2 * stat$eaf * (1 - stat$eaf)) * stat$seb^2)
    if (!is.null(excludeMAF)) 
        stat <- stat[stat$eaf > excludeMAF, ]
    if (!is.null(excludeMAF)) 
        stat <- stat[stat$eaf < (1 - excludeMAF), ]
    if (!is.null(stat$info)) {
        lowINFO <- stat$info < excludeINFO
        message(paste("Number of markers excluded by low INFO score:", 
            sum(lowINFO)))
        message("")
        stat <- stat[!lowINFO, ]
    }
    clsstat <- c("rsids", "chr", "pos", "ea", "nea", "eaf", "b", 
        "seb", "p", "n", "info") %in% colnames(stat)
    clsstat <- c("rsids", "chr", "pos", "ea", "nea", "eaf", "b", 
        "seb", "p", "n", "info")[clsstat]
    return(stat[, clsstat])
}
