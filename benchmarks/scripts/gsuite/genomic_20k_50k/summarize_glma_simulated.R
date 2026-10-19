# Compact evidence from fresh-process association benchmarks.
# Arguments: destination, fixture RDS, one or more run directories.
args <- commandArgs(TRUE); out <- args[1]; fixture_path <- args[2]; runs <- args[-c(1,2)]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fixture <- readRDS(fixture_path)
saveRDS(list(y = fixture$y, causal = fixture$causal, seed = fixture$seed),
        file.path(out, 'phenotype.rds'))
writeLines(unlist(fixture$Glist$rsids, use.names = FALSE),
           file.path(out, 'selected-markers.txt'))
excluded <- do.call(rbind, lapply(runs, function(run) {
  process <- read.csv(file.path(run, 'process-measurements.csv'), dec = ',')
  # Earlier harnesses did not record CPU time.
  if (!'rProcessCpuSeconds' %in% names(process)) process$rProcessCpuSeconds <- NA_real_
  process[process$exitCode != 0, ]
}))
write.csv(excluded, file.path(out, 'excluded-processes.csv'), row.names = FALSE)
measurements <- do.call(rbind, lapply(runs, function(run) {
  process <- read.csv(file.path(run, 'process-measurements.csv'), dec = ',')
  process <- process[process$exitCode == 0 &
    process$stage %in% c('glma_linear', 'glma_loco'), ]
  if (!nrow(process)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(process)), function(i) {
    row <- process[i, ]; timing <- read.csv(file.path(row$path, 'timing.csv'))
    fit <- readRDS(file.path(row$path, 'output.rds'))
    stopifnot(timing$status == 'completed', row$exitCode == 0,
              nrow(fit$associations) == sum(lengths(fixture$Glist$rsids)),
              all(fit$associations$status == 'ok'))
    if (row$stage == 'glma_loco') stopifnot(nrow(fit$fit$loco) == 22,
      all(is.finite(fit$fit$loco$genetic_variance)),
      all(is.finite(fit$fit$loco$residual_variance)),
      all(fit$fit$loco$genetic_variance > 0), all(fit$fit$loco$residual_variance > 0))
    data.frame(run = basename(run), method = row$stage, threads = row$threads,
      individuals = length(fixture$y), markers = nrow(fit$associations),
      fitted_chromosomes = if (is.null(fit$fit$loco)) NA_integer_ else nrow(fit$fit$loco),
      genetic_variance = if (is.null(fit$fit$genetic_variance)) NA_real_ else fit$fit$genetic_variance,
      residual_variance = if (is.null(fit$fit$residual_variance)) NA_real_ else fit$fit$residual_variance,
      calibration_factor = if (is.null(fit$diagnostics$calibration_factor)) NA_real_ else fit$diagnostics$calibration_factor,
      analysis_seconds = timing$analysis_seconds,
      native_fit_seconds = if (is.null(fit$diagnostics$fit_seconds)) NA_real_ else fit$diagnostics$fit_seconds,
      native_scan_seconds = fit$diagnostics$scan_seconds,
      other_call_seconds = timing$analysis_seconds - fit$diagnostics$scan_seconds -
        if (is.null(fit$diagnostics$fit_seconds)) 0 else fit$diagnostics$fit_seconds,
      serialization_seconds = timing$serialization_seconds,
      process_wall_seconds = row$processWallSeconds,
      process_cpu_seconds = if ('rProcessCpuSeconds' %in% names(row)) row$rProcessCpuSeconds else NA_real_,
      peak_resident_MiB = row$peakWorkingSetBytes / 1024^2,
      sampled_peak_private_MiB = row$sampledPeakPrivateBytes / 1024^2,
      path = row$path)
  }))
}))
write.csv(measurements, file.path(out, 'measurements.csv'), row.names = FALSE)
causal <- do.call(rbind, lapply(seq_len(nrow(measurements)), function(i) {
    measurement <- measurements[i, ]
    fit <- readRDS(file.path(measurement$path, 'output.rds'))
    rows <- fit$associations[match(fixture$causal, fit$associations$marker), ]
    data.frame(run = measurement$run, method = measurement$method,
               threads = measurement$threads, rows)
  }))
write.csv(causal, file.path(out, 'causal-associations.csv'), row.names = FALSE)
inputs <- c(fixture_path, fixture$Glist$bedfiles, fixture$Glist$bimfiles, fixture$Glist$famfiles)
write.csv(data.frame(path = inputs, bytes = file.info(inputs)$size,
  md5 = unname(tools::md5sum(inputs))), file.path(out, 'inputs.csv'), row.names = FALSE)
print(measurements[, setdiff(names(measurements), 'path')])
