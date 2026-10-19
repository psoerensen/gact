# Focused same-binary thread scaling for the native linear scan and BED scorer.
# Run from the gsuite repository root after the representative fixture exists.

updated_library <- file.path("build", "task-packages", "library")
task_library <- file.path("build", "task-packages", "ecosystem", "library")
root_library <- file.path("build", "r-library")
.libPaths(c(updated_library, task_library, root_library, .libPaths()))

for (package in c("gbase", "glma", "gscore")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Missing package: ", package)
}

Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OMP_DYNAMIC = "FALSE")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1L)

arguments <- commandArgs(TRUE)
block_argument <- grep("^--block-size=", arguments, value = TRUE)
block_size <- if (length(block_argument))
  as.integer(sub("^--block-size=", "", block_argument[[1L]])) else 128L
if (length(block_size) != 1L || is.na(block_size) || block_size < 1L)
  stop("--block-size must be one positive integer")
glma_only <- "--glma-only" %in% arguments

fixture <- file.path("build", "benchmarks", "genomic-20k-50k")
simulation_path <- file.path(fixture, "data_simulation.rds")
bayesian_path <- file.path(fixture, "bayesian_regression.rds")
if (!file.exists(simulation_path) || !file.exists(bayesian_path))
  stop("The cached representative simulation and Bayesian weights are required.")

simulation <- readRDS(simulation_path)$value
bayesian <- readRDS(bayesian_path)$value
n <- 20000L
m <- 50000L
individuals <- sprintf("id%05d", seq_len(n))
markers <- sprintf("m%05d", seq_len(m))
gwas_rows <- seq_len(6000L)
holdout_rows <- 18001:20000

# Preparation is deliberately outside the timed region.
gwas_glist <- gbase::gprep(
  bedfiles = paste0(simulation$prefix, ".bed"),
  ids = individuals[gwas_rows], nthreads = 1L)
score_glist <- gbase::gprep(
  bedfiles = paste0(simulation$prefix, ".bed"),
  ids = individuals[holdout_rows], nthreads = 1L)
phenotype <- simulation$simulation$Y[, "A"]
names(phenotype) <- individuals

regions <- split(markers, rep(sprintf("region%03d", seq_len(500L)), each = 100L))
gene_sets <- list()
for (region in seq_len(500L)) {
  cut <- 20L + region %% 60L
  names_here <- sprintf("gene%04d", 2L * region - c(1L, 0L))
  gene_sets[[names_here[[1L]]]] <- regions[[region]][seq_len(cut)]
  gene_sets[[names_here[[2L]]]] <- regions[[region]][(cut + 1L):100L]
}

thread_order <- c(1L, 2L, 4L, 2L, 4L, 1L, 4L, 1L, 2L)
repetition <- ave(seq_along(thread_order), thread_order, FUN = seq_along)
records <- vector("list", length(thread_order) * 2L)
reference_glma <- NULL
reference_score <- NULL
record_index <- 0L

time_call <- function(expression) {
  gc(FALSE)
  timing <- system.time(value <- force(expression), gcFirst = FALSE)
  list(value = value, elapsed = unname(timing[["elapsed"]]),
       user = unname(timing[["user.self"]]),
       system = unname(timing[["sys.self"]]))
}

for (run in seq_along(thread_order)) {
  threads <- thread_order[[run]]
  repeat_id <- repetition[[run]]
  message("START|glma|threads=", threads, "|repeat=", repeat_id)
  measured <- time_call(glma::glma(
    phenotype, gwas_glist, method = "linear",
    threads = threads, block_size = block_size))
  fit <- measured$value
  key <- fit$associations[c("beta", "se", "p")]
  if (is.null(reference_glma)) reference_glma <- key
  max_difference <- max(abs(as.matrix(key) - as.matrix(reference_glma)),
                        na.rm = TRUE)
  record_index <- record_index + 1L
  records[[record_index]] <- data.frame(
    kernel = "glma_linear", threads_requested = threads,
    threads_used = fit$diagnostics$threads_used, block_size = block_size,
    repetition = repeat_id, elapsed_seconds = measured$elapsed,
    user_seconds = measured$user, system_seconds = measured$system,
    native_seconds = fit$diagnostics$scan_seconds,
    max_abs_difference = max_difference)
  message("DONE|glma|threads=", threads, "|repeat=", repeat_id,
          "|native_seconds=", round(fit$diagnostics$scan_seconds, 4))
  rm(fit, key, measured)

  if (glma_only) next
  message("START|gscore|threads=", threads, "|repeat=", repeat_id)
  measured <- time_call(gscore::gscore(
    bayesian$weights, score_glist, sets = gene_sets[1:20],
    include_total = TRUE, coding = "standardized", scale_scores = FALSE,
    control = list(threads = threads)))
  scores <- measured$value$scores
  if (is.null(reference_score)) reference_score <- scores
  max_difference <- max(abs(scores - reference_score), na.rm = TRUE)
  record_index <- record_index + 1L
  records[[record_index]] <- data.frame(
    kernel = "gscore_bed_sets", threads_requested = threads,
    threads_used = NA_integer_, block_size = NA_integer_,
    repetition = repeat_id, elapsed_seconds = measured$elapsed,
    user_seconds = measured$user, system_seconds = measured$system,
    native_seconds = NA_real_, max_abs_difference = max_difference)
  message("DONE|gscore|threads=", threads, "|repeat=", repeat_id,
          "|elapsed_seconds=", round(measured$elapsed, 4))
  rm(scores, measured)
}

results <- do.call(rbind, records)
results$openblas_threads <- if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_get_num_procs() else NA_integer_
results$openmp_environment <- Sys.getenv("OMP_NUM_THREADS", unset = "")
results$openmp_dynamic <- Sys.getenv("OMP_DYNAMIC", unset = "")

output <- file.path(fixture, "focused-thread-scaling.csv")
write.csv(results, output, row.names = FALSE)
summary <- aggregate(
  cbind(elapsed_seconds, user_seconds, system_seconds, native_seconds) ~
    kernel + threads_requested,
  data = results, FUN = median, na.action = na.pass)
summary$speedup_elapsed <- NA_real_
for (kernel in unique(summary$kernel)) {
  rows <- which(summary$kernel == kernel)
  baseline <- summary$elapsed_seconds[rows[summary$threads_requested[rows] == 1L]]
  summary$speedup_elapsed[rows] <- baseline / summary$elapsed_seconds[rows]
}
summary_output <- file.path(fixture, "focused-thread-scaling-summary.csv")
write.csv(summary, summary_output, row.names = FALSE)
print(results, row.names = FALSE)
print(summary, row.names = FALSE)
message("RESULTS|", normalizePath(output, winslash = "/"))
message("SUMMARY|", normalizePath(summary_output, winslash = "/"))
