# Untimed agreement checks after the standalone performance measurements.
args <- commandArgs(TRUE)
runs <- jsonlite::fromJSON(args[[1L]])
load_fit <- function(label) {
  run <- runs[runs$label == label, , drop = FALSE]
  if (!nrow(run)) return(NULL)
  process <- read.csv(file.path(run$path, 'process-measurements.csv'))
  if (process$exitCode != 0 || isTRUE(process$timedOut)) return(NULL)
  readRDS(file.path(process$path, 'output.rds'))
}
checks <- list()
compare <- function(left, right, label, tolerance) {
  a <- load_fit(left); b <- load_fit(right)
  if (is.null(a) || is.null(b)) {
    checks[[label]] <<- 'Incomplete fit; agreement not checked'
    return(invisible(NULL))
  }
  stopifnot(identical(a$associations$marker, b$associations$marker),
            identical(a$associations$status, b$associations$status))
  error <- max(abs(a$associations$statistic - b$associations$statistic))
  stopifnot(is.finite(error), error < tolerance)
  if (!is.null(a$fit$loco)) {
    stopifnot(identical(a$fit$loco$chromosome, b$fit$loco$chromosome),
      isTRUE(all.equal(a$fit$loco$genetic_variance,
                       b$fit$loco$genetic_variance, tolerance = tolerance)),
      isTRUE(all.equal(a$fit$loco$residual_variance,
                       b$fit$loco$residual_variance, tolerance = tolerance)),
      abs(a$diagnostics$calibration_factor - b$diagnostics$calibration_factor) < tolerance)
  } else stopifnot(identical(a$fit$selected_prior, b$fit$selected_prior))
  checks[[label]] <<- paste('Passed; maximum absolute statistic difference:',
                            format(error, digits = 8))
}
compare('baseline-small-lmm-1', 'updated-small-lmm-1', 'before_after', 1e-5)
compare('updated-small-lmm-1', 'updated-small-lmm-4', 'ridge_threads', 1e-5)
compare('updated-frozen-prior-1', 'updated-frozen-prior-4', 'prior_threads', 1e-5)
writeLines(capture.output(checks), args[[2L]])
print(checks)
