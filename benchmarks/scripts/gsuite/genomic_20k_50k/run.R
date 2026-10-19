# One representative timing replicate. Run from a writable working directory.
updated_library <- file.path("build", "task-packages", "library")
task_library <- file.path("build", "task-packages", "ecosystem", "library")
root_library <- file.path("build", "r-library")
local_libraries <- c(updated_library, task_library, root_library)
.libPaths(c(local_libraries[dir.exists(local_libraries)], .libPaths()))

packages <- c("gbase", "gsim", "glma", "gcorr", "gbayes", "gscore", "gmap", "gsea")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

arguments <- commandArgs(TRUE)
thread_argument <- grep("^--threads=", arguments, value = TRUE)
thread_count <- if (length(thread_argument))
  as.integer(sub("^--threads=", "", thread_argument[[1L]])) else 1L
if (length(thread_count) != 1L || is.na(thread_count) || thread_count < 1L)
  stop("--threads must be one positive integer")
base_out <- file.path("build", "benchmarks", "genomic-20k-50k")
out <- if (thread_count == 1L) base_out else
  file.path("build", "benchmarks", sprintf("genomic-20k-50k-%dcores", thread_count))
results <- file.path("benchmarks", "genomic_20k_50k", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(results, recursive = TRUE, showWarnings = FALSE)
resume <- "--resume" %in% arguments
result_path <- function(stem) file.path(results, if (thread_count == 1L)
  paste0(stem, ".csv") else sprintf("%s-%dcores.csv", stem, thread_count))

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
rho <- seq(.2, .8, length.out = n_regions)
gwas_rows <- setNames(lapply(0:2, function(k) k * 6000L + seq_len(6000L)), traits)
reference_rows <- seq_len(18000L)
holdout_rows <- 18001:20000

timings <- list()
record_timing <- function(name, seconds, started, status = "completed", warnings = character()) {
  timings[[name]] <<- data.frame(
    replicate = 1L, stage = name, seconds = seconds,
    minutes = seconds / 60, started = format(started, tz = "Europe/Copenhagen"),
    status = status, warnings = paste(unique(warnings), collapse = "; "),
    stringsAsFactors = FALSE)
}
stage <- function(name, expression) {
  path <- file.path(out, paste0(name, ".rds"))
  if (resume && file.exists(path)) {
    saved <- readRDS(path)
    record_timing(name, saved$seconds, saved$started, "completed", saved$warnings)
    message("REUSE|", name, "|", round(saved$seconds, 2), " seconds")
    return(saved$value)
  }
  started <- Sys.time()
  warning_text <- character()
  message("START|", name)
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
  }), error = function(e) {
    record_timing(name, as.numeric(difftime(Sys.time(), started, units = "secs")),
      started, "failed", c(warning_text, conditionMessage(e)))
    write.csv(do.call(rbind, timings), result_path("timings"), row.names = FALSE)
    stop(e)
  })
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  saveRDS(list(value = value, seconds = seconds, started = started,
    warnings = unique(warning_text)), path)
  record_timing(name, seconds, started, warnings = warning_text)
  message("DONE|", name, "|", round(seconds, 2), " seconds")
  value
}

simulation <- if (thread_count > 1L &&
    file.exists(file.path(base_out, "data_simulation.rds"))) {
  frozen <- readRDS(file.path(base_out, "data_simulation.rds"))
  record_timing("data_simulation", frozen$seconds, frozen$started, status = "shared",
    warnings = frozen$warnings)
  frozen$value
} else stage("data_simulation", {
  set.seed(seed)
  prefix <- file.path(out, "genotypes")
  connection <- file(paste0(prefix, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), connection)
  frequencies <- numeric(m)
  dosage_sd <- numeric(m)
  for (j in seq_len(m)) {
    genotype <- rbinom(n, 2, .3)
    r <- (j - 1L) %/% region_size + 1L
    if ((j - 1L) %% region_size) {
      copied <- runif(n) < rho[r]
      genotype[copied] <- previous[copied]
    }
    previous <- genotype
    frequencies[j] <- mean(genotype[gwas_rows$A]) / 2
    dosage_sd[j] <- stats::sd(genotype[gwas_rows$A])
    codes <- c(3L, 2L, 0L)[genotype + 1L]
    writeBin(as.raw(colSums(matrix(codes, 4L) * c(1, 4, 16, 64))), connection)
    if (j %% 10000L == 0L) message("GENERATED|", j, " markers")
  }
  close(connection)
  write.table(data.frame(rep(seq_len(n_regions), each = region_size), ids, 0,
    as.character(rep(seq_len(region_size) * 1000L, n_regions)), "A", "G"),
    paste0(prefix, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(individuals, individuals, 0, 0, 0, -9),
    paste0(prefix, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  initial_Glist <- gbase::gprep(bedfiles = paste0(prefix, ".bed"))
  active <- runif(m) < ifelse(rep(seq_len(n_regions), each = region_size) <= 80, .08, .008)
  effects <- matrix(0, m, 3, dimnames = list(ids, traits))
  effect_correlation <- matrix(.6, 3, 3); diag(effect_correlation) <- 1
  effects[active, ] <- matrix(rnorm(sum(active) * 3, sd = .035), sum(active), 3) %*%
    chol(effect_correlation)
  effects[1:200, ] <- 0
  effects[50, ] <- .5
  effects[125, 1] <- .5
  effects[175, 2:3] <- .5
  generated <- gsim::gsim(Glist = initial_Glist, nt = 3, architecture = "fixed",
    beta = effects, h2 = c(.3, .4, .35), vg = c(.3, .4, .35), re = 0,
    standardize_W = FALSE, scale_effects = TRUE, seed = seed + 1L,
    compute_sumstats = FALSE, chunk_size = 64)
  colnames(generated$Y) <- colnames(generated$G) <- colnames(generated$B) <- traits
  list(prefix = prefix, simulation = generated, frequencies = frequencies,
    dosage_sd = dosage_sd)
})

Glist <- stage("glist_preparation", {
  gbase::gprep(bedfiles = paste0(simulation$prefix, ".bed"),
    nthreads = thread_count)
})

LD_Glist <- stage("sparse_ld_and_scores", {
  reference <- gbase::gprep(bedfiles = paste0(simulation$prefix, ".bed"),
    ids = individuals[reference_rows], nthreads = thread_count)
  reference <- gbase::gprep(Glist = reference, task = "sparseld",
    out_prefix = file.path(out, "LD"), max_distance_bp = 0,
    max_distance_variants = region_size, r2 = 0, nthreads = thread_count,
    block_size = 128, overwrite = TRUE)
  reference <- gbase::define_ld_blocks(reference, data.frame(
    block_id = names(regions), first_marker = vapply(regions, `[`, character(1), 1L),
    last_marker = vapply(regions, function(x) x[[length(x)]], character(1)),
    stringsAsFactors = FALSE))
  reference$assembly <- "artificial-20k-50k"
  scores <- gbase::ld_scores(reference)
  stopifnot(length(scores) == m, all(is.finite(scores)))
  reference
})

association <- stage("gwas_and_checkstat", {
  fits <- statistics <- setNames(vector("list", length(traits)), traits)
  metadata <- LD_Glist$sparseLD$marker_metadata
  for (trait in traits) {
    selected <- gbase::gprep(bedfiles = paste0(simulation$prefix, ".bed"),
      ids = individuals[gwas_rows[[trait]]], nthreads = thread_count)
    phenotype <- simulation$simulation$Y[, trait]
    names(phenotype) <- individuals
    fits[[trait]] <- glma::glma(phenotype, selected, method = "linear",
      threads = thread_count, block_size = 128)
    a <- fits[[trait]]$associations
    raw <- data.frame(rsids = ids, chr = metadata$chromosome,
      pos = metadata$base_pair_position, ea = metadata$allele1,
      nea = metadata$allele2, b = a$beta, seb = a$se, p = a$p,
      n = length(gwas_rows[[trait]]), stringsAsFactors = FALSE)
    statistics[[trait]] <- suppressMessages(gbase::checkStat(LD_Glist, raw,
      excludeMAF = NULL, excludeMAFDIFF = NULL, excludeINFO = NULL,
      excludeCGAT = FALSE, excludeINDEL = FALSE, excludeDUPS = TRUE,
      excludeMHC = FALSE, excludeMISS = NULL, excludeHWE = NULL))
    stopifnot(nrow(statistics[[trait]]) == m)
  }
  list(fits = fits, statistics = statistics)
})

correlation <- stage("ld_score_regression", {
  ordinary <- gcorr::ldsc(association$statistics, LD_Glist, what = "rg",
    block_count = 200L, overlap_intercept = 0,
    marginal_intercept = "fixed_one")
  partitioned <- gcorr::ldsc(association$statistics, LD_Glist,
    sets = list(enriched = ids[seq_len(80L * region_size)]), residual = TRUE,
    what = "rg", block_count = 200L, overlap_intercept = 0,
    marginal_intercept = "fixed_one")
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
  metadata <- setNames(lapply(traits, function(trait) data.frame(
    set = names(gene_sets), sample_size = length(gwas_rows[[trait]]),
    mean_minor_allele_count = vapply(gene_sets, function(markers) {
      p <- simulation$frequencies[match(markers, ids)]
      mean(2 * length(gwas_rows[[trait]]) * pmin(p, 1 - p))
    }, numeric(1)), stringsAsFactors = FALSE)), traits)
  glma::vegas(association$statistics, LD_Glist, sets = gene_sets,
    metadata = metadata, control = list(tail_method = "saddlepoint"))
})

gene_models <- stage("gene_models", {
  input <- gsea::gsea_input(gene_evidence)
  gene_stat <- input$stat
  magma <- gsea::gsea(gene_stat, pathway_sets, method = "magma",
    sampling = list(gene_covariance = "independent"))
  gene_names <- names(gene_sets)
  region_index <- rep(seq_len(n_regions), each = 2L)
  set.seed(seed + 2L)
  features <- matrix(rnorm(length(gene_names) * 100L), nrow = length(gene_names),
    dimnames = list(gene_names, sprintf("feature%03d", seq_len(100L))))
  features[, 1L] <- as.numeric(gene_names %in% pathway_sets$enriched)
  pops_stat <- data.frame(gene = gene_names,
    chromosome = as.character(ceiling(region_index / 25L)),
    score = gene_stat$A$score, stringsAsFactors = FALSE)
  pops <- gsea::pops(pops_stat, features,
    sampling = list(gene_covariance = "independent"),
    control_features = "feature001")
  list(magma = magma, pops = pops)
})

bayesian <- stage("bayesian_regression", {
  gbayes::gbayes(association$statistics$A, LD_Glist, method = "bayesc",
    trait = "A", prior = list(residual_variance = 1, effect_variance = .3 / (m * .02),
      inclusion_probability = .02),
    control = list(burnin = 500L, sampling_sweeps = 700L,
      seeds = c(11, 29, 47, 71), threads = thread_count, residual_policy = "fixed"))
})

mapping <- stage("regional_fine_mapping", {
  phenotype_variance <- setNames(vapply(traits, function(trait)
    var(simulation$simulation$Y[gwas_rows[[trait]], trait]), numeric(1)), traits)
  prepared <- gmap::gmap_stat(association$statistics[1:2], LD_Glist,
    phenotype_variance = phenotype_variance[1:2])
  gmap::gmap(prepared, LD_Glist, regions = regions[1:20], method = "multi_effect",
    prior = list(residual_variance = 1, effect_variance = .04),
    control = list(effects = 1L, max_sweeps = 200L,
      eigen = list(retained_positive_mass = c(.995, .99, .95))))
})

prediction <- stage("polygenic_scoring", {
  test <- gbase::gprep(bedfiles = paste0(simulation$prefix, ".bed"),
    ids = individuals[holdout_rows], nthreads = thread_count)
  gscore::gscore(bayesian$weights, test, sets = gene_sets[1:20],
    include_total = TRUE, coding = "standardized", scale_scores = FALSE,
    control = list(threads = thread_count))
})

timing_table <- do.call(rbind, timings)
write.csv(timing_table, result_path("timings"), row.names = FALSE)

hardware <- tryCatch({
  cpu <- utils::readRegistry(
    "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", "HLM")$ProcessorNameString
  memory_command <- paste("Add-Type -AssemblyName Microsoft.VisualBasic;",
    "[math]::Round((New-Object Microsoft.VisualBasic.Devices.ComputerInfo).TotalPhysicalMemory/1GB,1)")
  memory <- system2("powershell.exe", c("-NoProfile", "-Command",
    shQuote(memory_command)), stdout = TRUE, stderr = FALSE)
  list(cpu = cpu, memory_gb = sub(",", ".", memory, fixed = TRUE))
}, error = function(e) list(cpu = NA_character_, memory_gb = NA_character_))
workflow_commit <- tryCatch({
  revision <- suppressWarnings(system2("git", c("rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE))
  if (length(revision) == 1L && grepl("^[0-9a-f]{40}$", revision))
    revision else NA_character_
}, error = function(e) NA_character_)
provenance <- data.frame(
  key = c("description", "individuals", "markers", "traits", "gwas_samples_per_trait",
    "ld_reference_samples", "holdout_samples", "thread_limit", "thread_policy", "seed", "os", "cpu",
    "physical_memory_gb", "R", "workflow_commit", paste0("package_", packages)),
  value = c("One representative execution; not a formal benchmark", n, m, length(traits),
    length(gwas_rows$A), length(reference_rows), length(holdout_rows), thread_count,
    "Up to the thread limit where the public interface supports it", seed,
    paste(Sys.info()[c("sysname", "release", "version")], collapse = " "),
    hardware$cpu, hardware$memory_gb, R.version.string,
    workflow_commit,
    vapply(packages, function(package) as.character(utils::packageVersion(package)), character(1))),
  stringsAsFactors = FALSE)
write.csv(provenance, result_path("provenance"), row.names = FALSE)

summary_table <- transform(timing_table[c("stage", "seconds", "minutes", "status")],
  workload = c("20K individuals; 50K markers; 3 traits", "20K individuals; 50K markers",
    "18K-reference individuals; 50K markers", "3 x 6K-individual linear scans; 50K markers",
    "3 traits; ordinary and partitioned LDSC", "1,000 genes x 3 traits",
    "MAGMA-style pathways plus PoPS-style ranking", "BayesC; 50K markers; 4 chains",
    "20 regions x 2 traits", "2K held-out individuals; total plus 20 sets"))
write.csv(summary_table, result_path("summary"), row.names = FALSE)
message("BENCHMARK|completed|", normalizePath(results, winslash = "/"))
