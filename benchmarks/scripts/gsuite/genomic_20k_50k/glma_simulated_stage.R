# Individual-level association on the frozen gact simulated-human example.
# Use run_manuscript.ps1 -GlmaComparison; each fit is a fresh R process.
args <- commandArgs(TRUE)
stage <- args[1]; threads <- as.integer(args[2]); out <- args[3]
root <- dirname(out)
writeLines(as.character(Sys.getpid()), file.path(out, 'r-process.pid'))
.libPaths(c(Sys.getenv('GLMA_BENCH_LIBRARY'), 'build/task-packages/library', .libPaths()))
if (stage != 'glma_fixture') library(glma)
writeLines(c(capture.output(sessionInfo()), if (stage != 'glma_fixture') capture.output(glma::glma_native_info())),
           file.path(out, 'session.txt'))
frozen_path <- 'build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds'
glist_path <- 'build/benchmarks/manuscript-20261005-231349/cmp_suite_prepare-1/output.rds'
if (stage != 'glma_fixture') {
  fixture_path <- if (nzchar(Sys.getenv('GLMA_SIM_FIXTURE'))) Sys.getenv('GLMA_SIM_FIXTURE') else
    file.path(root, paste0('glma_fixture-', threads), 'output.rds')
  input <- readRDS(fixture_path)
}
gc()
started <- Sys.time(); clock <- proc.time()[3]; warnings <- character()
value <- tryCatch(withCallingHandlers(switch(stage,
  glma_fixture = {
    library(qgg)
    writeLines(c(capture.output(sessionInfo()),
      paste('qgg', packageVersion('qgg'), find.package('qgg'))), file.path(out, 'session.txt'))
    frozen <- readRDS(frozen_path)
    Glist <- readRDS(glist_path)
    stopifnot(identical(unname(unlist(Glist$rsids)), unname(frozen$rsids)))
    stopifnot(identical(unname(Glist$ids), unname(unlist(frozen$Glist$ids))),
              identical(as.numeric(Glist$sample_source_indices_0based),
                        as.numeric(seq_along(Glist$ids) - 1L)))
    set.seed(frozen$seed)
    sim <- qgg::gsim(Glist = frozen$Glist, rsids = frozen$rsids)
    stopifnot(identical(sim$causal, frozen$causal))
    replay <- qgg::glma(y = sim$y, Glist = frozen$Glist, rsids = frozen$rsids)
    audit <- all.equal(replay, frozen$stat, tolerance = 1e-10)
    writeLines(capture.output(audit), file.path(out, 'phenotype-replay-audit.txt'))
    stopifnot(isTRUE(audit))
    y <- as.numeric(sim$y)
    stopifnot(length(y) == length(Glist$ids), all(is.finite(y)))
    names(y) <- Glist$ids
    paths <- c(frozen_path, glist_path, Glist$bedfiles, Glist$bimfiles, Glist$famfiles)
    write.csv(data.frame(path = paths, bytes = file.info(paths)$size,
                        md5 = unname(tools::md5sum(paths))),
              file.path(out, 'source-files.csv'), row.names = FALSE)
    list(Glist = Glist, y = y, causal = sim$causal,
         replay_audit = audit, seed = frozen$seed)
  },
  glma_linear = glma::glma_linear(phenotype = input$y, Glist = input$Glist,
                                 covariates = input$covariates, threads = threads),
  glma_loco = glma::glma_infinitesimal_loco(phenotype = input$y,
                                          Glist = input$Glist,
                                          covariates = input$covariates, threads = threads),
  glma_prior = glma::glma(phenotype = input$y, Glist = input$Glist,
                         covariates = input$covariates,
                         method = 'prior_selected_loco', threads = threads),
  stop('Unknown glma benchmark stage')),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }),
  error = function(e) e)
seconds <- unname(proc.time()[3] - clock)
failed <- inherits(value, 'error')
if (failed) warnings <- c(warnings, conditionMessage(value))
serialization <- unname(system.time(saveRDS(value, file.path(out, 'output.rds')))[3])
write.csv(data.frame(stage, threads, started = format(started, tz = 'Europe/Copenhagen'),
  analysis_seconds = seconds, serialization_seconds = serialization,
  status = if (failed) 'failed' else 'completed',
  warnings = paste(unique(warnings), collapse = '; ')), file.path(out, 'timing.csv'), row.names = FALSE)
writeLines(capture.output(str(value, max.level = 3)), file.path(out, 'output-structure.txt'))
if (!failed && stage != 'glma_fixture') {
  write.csv(value$associations, file.path(out, 'associations.csv'), row.names = FALSE)
  writeLines(capture.output(value$diagnostics), file.path(out, 'diagnostics.txt'))
  if (!is.null(value$fit$loco)) write.csv(value$fit$loco, file.path(out, 'loco.csv'), row.names = FALSE)
  stopifnot(nrow(value$associations) == sum(lengths(input$Glist$rsids)),
            all(value$associations$status == 'ok'),
            all(is.finite(value$associations$p)))
}
if (failed) quit(status = 1)
