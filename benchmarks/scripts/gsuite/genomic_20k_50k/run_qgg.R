# One representative qgg timing replicate on the frozen gsuite simulation.
# Run from the gsuite repository root after benchmarks/genomic_20k_50k/run.R.
qgg_library <- file.path("build", "benchmarks", "genomic-20k-50k", "qgg-library")
if (dir.exists(qgg_library)) .libPaths(c(qgg_library, .libPaths()))
if (!requireNamespace("qgg", quietly = TRUE)) stop("Install qgg before running this comparison.")

out <- file.path("build", "benchmarks", "genomic-20k-50k")
results <- file.path("benchmarks", "genomic_20k_50k", "results")
simulation_path <- file.path(out, "data_simulation.rds")
if (!file.exists(simulation_path)) stop("Run the gsuite benchmark first to freeze its simulation.")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
resume <- "--resume" %in% commandArgs(TRUE)

seed <- 20260930L
n <- 20000L
m <- 50000L
region_size <- 100L
n_regions <- m %/% region_size
traits <- c("A", "B", "C")
ids <- sprintf("m%05d", seq_len(m))
individuals <- sprintf("id%05d", seq_len(n))
region <- rep(sprintf("region%03d", seq_len(n_regions)), each = region_size)
regions <- split(ids, region)
gwas_rows <- setNames(lapply(0:2, function(k) k * 6000L + seq_len(6000L)), traits)
reference_rows <- seq_len(18000L)
holdout_rows <- 18001:20000
simulation <- readRDS(simulation_path)$value

timings <- list()
record_timing <- function(name, seconds, started, status = "completed", warnings = character()) {
  timings[[name]] <<- data.frame(
    replicate = 1L, implementation = "qgg", stage = name, seconds = seconds,
    minutes = seconds / 60, started = format(started, tz = "Europe/Copenhagen"),
    status = status, warnings = paste(unique(warnings), collapse = "; "),
    stringsAsFactors = FALSE)
}
stage <- function(name, expression) {
  path <- file.path(out, paste0("qgg-current_", name, ".rds"))
  if (resume && file.exists(path)) {
    saved <- readRDS(path)
    record_timing(name, saved$seconds, saved$started, "completed", saved$warnings)
    message("REUSE|qgg|", name, "|", round(saved$seconds, 2), " seconds")
    return(saved$value)
  }
  started <- Sys.time()
  warning_text <- character()
  message("START|qgg|", name)
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    record_timing(name, as.numeric(difftime(Sys.time(), started, units = "secs")),
      started, "failed", c(warning_text, conditionMessage(e)))
    write.csv(do.call(rbind, timings), file.path(results, "qgg-timings.csv"), row.names = FALSE)
    stop(e)
  })
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  saveRDS(list(value = value, seconds = seconds, started = started,
    warnings = unique(warning_text)), path)
  record_timing(name, seconds, started, warnings = warning_text)
  message("DONE|qgg|", name, "|", round(seconds, 2), " seconds")
  value
}

Glist <- stage("glist_preparation", {
  qgg::gprep(bedfiles = paste0(simulation$prefix, ".bed"),
    bimfiles = paste0(simulation$prefix, ".bim"),
    famfiles = paste0(simulation$prefix, ".fam"),
    assembly = "artificial-20k-50k")
})

LD_Glist <- stage("sparse_ld_and_scores", {
  qgg::gprep(Glist = Glist, task = "sparseld",
    ldfiles = file.path(out, "qgg-current-LD.ld"), ids = individuals[reference_rows],
    msize = region_size, overwrite = TRUE, ncores = 1)
})

association <- stage("gwas_and_checkstat", {
  statistics <- setNames(vector("list", length(traits)), traits)
  for (trait in traits) {
    phenotype <- simulation$simulation$Y[gwas_rows[[trait]], trait]
    names(phenotype) <- individuals[gwas_rows[[trait]]]
    raw <- qgg::glma(phenotype, Glist = Glist, msize = 128L, scale = TRUE)
    statistics[[trait]] <- qgg::checkStat(LD_Glist, raw,
      excludeMAF = 0, excludeMAFDIFF = 1, excludeINFO = 0,
      excludeCGAT = FALSE, excludeINDEL = FALSE, excludeDUPS = TRUE,
      excludeMHC = FALSE, excludeMISS = 1, excludeHWE = 0)
    rownames(statistics[[trait]]) <- statistics[[trait]]$rsids
    stopifnot(nrow(statistics[[trait]]) == m)
  }
  statistics
})

correlation <- stage("ld_score_regression", {
  z <- do.call(cbind, lapply(association, function(x) setNames(x$b / x$seb, x$rsids)))
  colnames(z) <- traits
  sample_sizes <- setNames(rep(length(gwas_rows$A), length(traits)), traits)
  ordinary <- qgg::ldsc(Glist = LD_Glist, z = z, n = sample_sizes,
    what = "rg", intercept = TRUE)
  partitioned <- qgg::ldsc(Glist = LD_Glist, z = z, n = sample_sizes,
    sets = list(enriched = ids[seq_len(80L * region_size)]), residual = TRUE,
    what = "h2", method = "regression")
  list(ordinary = ordinary, partitioned = partitioned)
})

gene_sets <- local({
  value <- list()
  for (r in seq_len(n_regions)) {
    cut <- 20L + r %% 60L
    names_here <- sprintf("gene%04d", 2L * r - c(1L, 0L))
    value[[names_here[[1L]]]] <- regions[[r]][seq_len(cut)]
    value[[names_here[[2L]]]] <- regions[[r]][(cut + 1L):region_size]
  }
  value
})
pathway_sets <- list(enriched = names(gene_sets)[1:160],
  background = names(gene_sets)[401:560],
  mixed = names(gene_sets)[seq(1, 999, by = 5)])

gene_evidence <- stage("vegas", {
  p <- do.call(cbind, lapply(association, function(x) setNames(x$p, x$rsids)))
  colnames(p) <- traits
  qgg::vegas(Glist = LD_Glist, sets = gene_sets, p = p, verbose = FALSE)
})

gene_models <- stage("gene_models", {
  gene_stat <- gene_evidence$z
  rownames(gene_stat) <- gene_evidence$genes$EnsemblID
  set.seed(seed + 2L)
  feature_matrix <- matrix(runif(nrow(gene_stat) * 100L) < .1,
    nrow = nrow(gene_stat), dimnames = list(rownames(gene_stat),
      sprintf("feature%03d", seq_len(100L))))
  feature_matrix[, 1L] <- rownames(gene_stat) %in% pathway_sets$enriched
  feature_sets <- lapply(seq_len(ncol(feature_matrix)), function(j)
    rownames(feature_matrix)[feature_matrix[, j]])
  names(feature_sets) <- colnames(feature_matrix)
  magma <- lapply(seq_len(ncol(gene_stat)), function(j)
    qgg::magma(gene_stat[, j], pathway_sets, method = "magma"))
  names(magma) <- colnames(gene_stat)
  pops <- lapply(seq_len(ncol(gene_stat)), function(j) {
    set.seed(seed + 2L + j)
    qgg::pops(gene_stat[, j], feature_sets, method = "bayesC",
      pi = .02, nit = 1200L, nburn = 500L,
      updateG = FALSE, updateE = FALSE, updatePi = FALSE)
  })
  names(pops) <- colnames(gene_stat)
  list(magma = magma, pops = pops)
})

bayesian <- stage("bayesian_regression", {
  lapply(seq_along(c(11L, 29L, 47L, 71L)), function(i) {
    set.seed(c(11L, 29L, 47L, 71L)[i])
    qgg::gbayes(stat = association$A, Glist = LD_Glist,
      method = "bayesC", algorithm = "mcmc", nit = 1200L, nburn = 500L,
      h2 = .3, pi = .02, vb = .3 / (m * .02), ve = 1,
      updateG = FALSE, updateE = FALSE, updatePi = FALSE, msize = 128L)
  })
})

mapping <- stage("regional_fine_mapping", {
  lapply(association[1:2], function(stat) qgg::gmap(
    Glist = LD_Glist, stat = stat, sets = regions[1:20],
    method = "bayesC", algorithm = "mcmc-eigen",
    eigen_threshold = c(.995, .99, .95), nit = 1200L, nburn = 500L,
    nrun = 1L, h2 = .3, pi = .02, updateG = FALSE,
    updateE = FALSE, updatePi = FALSE, verbose = FALSE))
})

prediction <- stage("polygenic_scoring", {
  weights <- bayesian[[1L]]$stat
  total <- qgg::gscore(Glist = Glist, stat = weights,
    ids = individuals[holdout_rows], scaleMarker = TRUE, scaleGRS = FALSE,
    msize = 128L, ncores = 1L)
  sets <- qgg::gscore(Glist = Glist, stat = weights, sets = gene_sets[1:20],
    ids = individuals[holdout_rows], scaleMarker = TRUE, scaleGRS = FALSE,
    msize = 128L, ncores = 1L)
  list(total = total, sets = sets)
})

timing_table <- do.call(rbind, timings)
write.csv(timing_table, file.path(results, "qgg-timings.csv"), row.names = FALSE)
provenance <- data.frame(
  key = c("description", "simulation", "qgg_version", "qgg_commit",
    "individuals", "markers", "traits", "threads", "seed",
    "comparison_limit_ldsc", "comparison_limit_gene_models",
    "comparison_limit_fine_mapping"),
  value = c("One representative qgg execution on the frozen gsuite simulation",
    "Reused data_simulation.rds; simulation time is shared and not rerun",
    as.character(utils::packageVersion("qgg")),
    tryCatch(system2("git", c("-c", paste0("safe.directory=", normalizePath("../qgg", winslash = "/")),
      "-C", "../qgg", "rev-parse", "HEAD"), stdout = TRUE),
      error = function(e) NA_character_),
    n, m, length(traits), 1L, seed,
    "qgg partitioned LDSC is h2-only; ordinary LDSC estimates rg",
    "Same gene evidence and named sets; qgg PoPS consumes set-valued features",
    "qgg uses BayesC mcmc-eigen; gsuite timing uses multi_effect eigen mapping"),
  stringsAsFactors = FALSE)
write.csv(provenance, file.path(results, "qgg-provenance.csv"), row.names = FALSE)
message("BENCHMARK|qgg|completed|", normalizePath(results, winslash = "/"))
