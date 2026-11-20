# Run against an installed gact. Optional arguments are marker counts.
# Rscript benchmarks/gwas-worker-memory.R 100000 500000
# Every study has an observed result for every backbone marker and one extra
# 80-byte source column. Each measurement starts fresh R worker processes.
library(gact)

peak_mib <- function() {
  if (.Platform$OS.type == "windows") {
    command <- sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())
    return(as.numeric(system2("powershell.exe", c("-NoProfile", "-Command", shQuote(command)),
                              stdout = TRUE)) / 1024^2)
  }
  value <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(sub("^VmHWM:[[:space:]]+([0-9]+).*", "\\1", value)) / 1024
}

run_worker <- function(i, root, index, input, headers, storage, peak_mib) {
  base <- file.path(root, paste0("study", i))
  full <- if (storage == "full") paste0(base, ".txt.gz") else NA_character_
  .Call(get("C_gact_gwas_stage", asNamespace("gact")), index, input,
        c(full, paste0(base, ".bin")), paste0(base, ".qc.json"),
        headers, c(headers[1:8], NA_character_), TRUE, NULL, FALSE)
  peak_mib()
}

measure <- function(counts) {
  root <- tempfile("gwas-worker-memory-")
  stopifnot(dir.create(root))
  on.exit(unlink(root, recursive = TRUE))
  native <- function(name) get(name, asNamespace("gact"))
  for (n in counts) {
    markers <- data.frame(rsids = paste0("rs", seq_len(n)), chr = 1L,
                          pos = seq_len(n), ea = "A", nea = "C", eaf = .2)
    catalog <- file.path(root, "markers.txt.gz")
    cpra <- file.path(root, "cpra.txt.gz")
    index <- file.path(root, "index.bin")
    data.table::fwrite(markers, catalog)
    data.table::fwrite(data.frame(cpra = paste(1L, markers$pos, "A", "C", sep = "_")),
                       cpra, col.names = FALSE)
    .Call(native("C_gact_gwas_index"), catalog, cpra, index)
    input <- file.path(root, "source.tsv.gz")
    headers <- c("chr", "pos", "ea", "nea", "b", "seb", "p", "eaf", "note")
    data.table::fwrite(data.frame(chr = 1L, pos = markers$pos, ea = "A", nea = "C",
      b = .2, seb = .1, p = .05, eaf = .2, note = strrep("x", 80)), input, sep = "\t")
    for (workers in c(1L, 2L)) for (storage in c("full", "z")) {
      cluster <- parallel::makePSOCKcluster(workers)
      peaks <- tryCatch({
        parallel::clusterCall(cluster, .libPaths, .libPaths())
        parallel::parLapply(cluster, seq_len(workers), run_worker,
                            root, index, input, headers, storage, peak_mib)
      }, finally = parallel::stopCluster(cluster))
      peaks <- unlist(peaks)
      cat(sprintf("%d\t%d\t%s\t%.1f\t%.1f\n", n, workers, storage, max(peaks), sum(peaks)))
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
counts <- if (length(args)) as.integer(args) else c(100000L, 500000L)
stopifnot(all(is.finite(counts) & counts > 0L))
cat("markers\tworkers\tstorage\tmax_worker_peak_MiB\tsum_worker_peaks_MiB\n")
measure(counts)
