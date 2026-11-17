.gact_cpp_gwas_preview <- function(file, rows, skip) {
  if (!grepl("\\.tsv\\.(gz|bgz)$", file, ignore.case = TRUE) || skip != 0L ||
      rows >= .Machine$integer.max) {
    stop("The C++ reader expects a gzip TSV with a header on the first line")
  }
  text <- .Call(C_gact_gwas_preview, file, as.integer(rows))
  if (length(strsplit(text, "\n", fixed = TRUE)[[1L]]) < 2L) stop("Preview has no data rows")
  data.table::fread(text = text, sep = "\t", nrows = rows, data.table = FALSE, check.names = FALSE)
}

.gact_cpp_gwas_mapping <- function(recipe) {
  fields <- c("chr", "pos", "ea", "nea", "b", "seb", "p", "eaf")
  mapping <- unlist(recipe$mapping)
  mapping <- mapping[!is.na(mapping)]
  header <- unlist(recipe$inspection$headers, use.names = FALSE)
  if (!all(c("marker", fields[-8L]) %in% names(mapping)) ||
      !all(names(mapping) %in% c("marker", fields)) ||
      anyDuplicated(header) || any(!nzchar(header)) ||
      any(grepl('[,"\t\r\n]', header)) ||
      any(c("n", "ncase", "ncontrol", "info", "z", "z_score") %in% header) ||
      !identical(recipe$effect_scale, "beta") || !identical(recipe$p_scale, "p") ||
      !is.null(recipe[["n"]])) {
    stop("Verified recipe is outside the C++ profile. Use the R engine for this study")
  }
  unname(as.character(mapping[fields]))
}

.gact_cpp_gwas_batch <- function(ga, descriptor, jobs, report, out_dir, workers) {
  marker <- ga$dirs[["marker"]]
  markers <- file.path(marker, "markers.txt.gz")
  cpra <- file.path(marker, "cpra.txt.gz")
  stopifnot(file.exists(markers), file.exists(cpra))
  index <- file.path(out_dir, "marker-index.bin")
  files <- c(markers, cpra)
  fingerprint <- tools::md5sum(files)
  ids <- data.table::fread(markers, select = "rsids")$rsids
  if (!identical(as.character(ga$rsids), as.character(ids))) stop("Marker catalog and GAlist disagree")
  rm(ids)
  .Call(C_gact_gwas_index, markers, cpra, index)
  recipe_dir <- file.path(dirname(ga$dirs[["gstat"]]), "recipes")
  dir.create(recipe_dir, showWarnings = FALSE)
  if (!dir.exists(recipe_dir)) stop("Cannot create database recipe directory")
  active <- which(lengths(jobs) > 0L)
  for (i in active) {
    recipe <- jobs[[i]]$recipe
    path <- file.path(recipe_dir, paste0(recipe$accession, "-", recipe$input_md5, ".json"))
    if (!file.exists(path)) {
      temp <- tempfile(tmpdir = recipe_dir)
      tryCatch({
        jsonlite::write_json(recipe, temp, auto_unbox = TRUE, pretty = TRUE)
        if (!file.rename(temp, path)) stop("Cannot save verified recipe: ", path)
      }, finally = unlink(temp))
    }
    stored <- jsonlite::read_json(path)
    fields <- c("accession", "input_md5", "mapping", .gact_gwas_metadata_fields)
    if (!identical(stored[fields], recipe[fields]) ||
        !identical(stored$inspection$headers, recipe$inspection$headers)) {
      stop("Existing recipe disagrees with reviewed recipe: ", path)
    }
    report$recipe[i] <- path
    jobs[[i]]$recipe_path <- path
  }

  stage_one <- function(i) {
    job <- jobs[[i]]
    tryCatch({
      file <- file.path(out_dir, sprintf("study_%05d", i), "staged.txt.gz")
      qc <- file.path(dirname(file), "staged.qc.json")
      recipe <- jsonlite::read_json(job$recipe_path)
      if (!identical(unname(tools::md5sum(job$file)), recipe$input_md5)) stop("Input checksum changed before staging")
      .Call(C_gact_gwas_stage, index, job$file, file, qc,
            unlist(recipe$inspection$headers, use.names = FALSE),
            .gact_cpp_gwas_mapping(recipe), identical(recipe$strand, "forward"))
      metrics <- jsonlite::read_json(qc, simplifyVector = TRUE)
      checksum <- unname(tools::md5sum(file))
      stopifnot(!is.na(checksum), metrics$n_final > 0,
                metrics$n_input == metrics$n_final + metrics$n_dropped,
                metrics$n_dropped == sum(unlist(metrics$dropped_by_reason)),
                metrics$n_flipped <= metrics$n_final)
      c(job, list(staged = file, qc_file = qc, qc = metrics, checksum = checksum))
    }, error = function(e) list(error = conditionMessage(e)))
  }
  count <- min(workers, length(active))
  staged <- if (.Platform$OS.type == "unix" && count > 1L) {
    parallel::mclapply(active, stage_one, mc.cores = count, mc.preschedule = FALSE)
  } else lapply(active, stage_one)
  passed <- vapply(staged, function(value) is.list(value) && !is.null(value$staged), logical(1))
  for (j in which(!passed)) {
    value <- staged[[j]]
    report$issues[active[j]] <- if (is.list(value) && length(value$error)) value$error else
      "C++ worker exited before producing a validated result"
  }
  if (!any(passed)) return(report)
  if (!identical(fingerprint, tools::md5sum(files))) stop("Marker catalog changed during staging")
  result <- tryCatch(.gact_publish_cpp_gwas(ga, descriptor, staged[passed]),
                     error = function(e) e)
  if (inherits(result, "error")) {
    report$issues[active[passed]] <- conditionMessage(result)
    return(report)
  }
  for (j in seq_along(result)) {
    i <- active[passed][j]
    item <- staged[passed][[j]]
    report$status[i] <- "imported"
    report$study_id[i] <- result[j]
    for (field in c("n_input", "n_final", "n_dropped", "n_flipped")) {
      report[[field]][i] <- item$qc[[field]]
    }
  }
  report
}

.gact_publish_cpp_gwas <- function(ga, descriptor, staged) {
  gstat <- ga$dirs[["gstat"]]
  inventory <- file.path(gstat, "GWAS_information.csv")
  stage <- file.path(gstat, ".gact-import.lock")
  if (!dir.create(stage, mode = "0700", showWarnings = FALSE)) stop("Another GWAS import holds the database lock")
  temp_rds <- tempfile(".gact-descriptor-", tmpdir = dirname(descriptor))
  backup <- file.path(stage, "previous-inventory.csv")
  published <- character()
  replaced <- committed <- FALSE
  on.exit({
    if (!committed) {
      if (replaced) {
        if (file.exists(backup)) {
          if (!file.rename(backup, inventory)) warning("Cannot restore the previous GWAS inventory")
        } else unlink(inventory)
      }
      unlink(published)
    }
    unlink(temp_rds)
    unlink(stage, recursive = TRUE)
  }, add = TRUE)
  current <- readRDS(descriptor)
  if (!identical(current$study, ga$study) || !identical(current$dirs, ga$dirs) ||
      !identical(current$rsids, ga$rsids)) stop("Database changed while studies were staging. Retry the batch")
  ids <- as.character(current$study$id)
  if (file.exists(inventory)) {
    if (!setequal(read.csv2(inventory, colClasses = "character")$id, ids) ||
        !file.copy(inventory, backup)) stop("Cannot verify and back up the current GWAS inventory")
  } else if (length(ids)) stop("GWAS inventory disappeared while studies were staging")
  last <- max(c(0L, suppressWarnings(as.integer(sub("^GWAS", "", ids)))), na.rm = TRUE)
  new_ids <- paste0("GWAS", last + seq_along(staged))
  source <- destination <- character()
  for (j in seq_along(staged)) {
    item <- staged[[j]]
    recipe <- item$recipe
    id <- new_ids[j]
    at <- length(current$study$id) + 1L
    current$study <- lapply(current$study, function(column) { length(column) <- at; column })
    recipeID <- sub("^GWAS", "REC", id)
    name <- paste0(id, ".txt")
    qc <- paste0(id, ".qc.json")
    recipe_path <- file.path("..", "recipes", basename(item$recipe_path))
    numbers <- lapply(c("n", "ncase", "ncontrol"), function(field) {
      if (is.null(recipe[[field]])) NA_real_ else as.numeric(recipe[[field]])
    })
    names(numbers) <- c("n", "ncase", "ncontrol")
    total <- numbers$ncase + numbers$ncontrol
    neff <- if (identical(recipe$type, "quantitative")) numbers$n else
      if (is.finite(total) && total > 0) numbers$ncase * numbers$ncontrol / total else NA_real_
    row <- c(list(id = id, file = name, created_at = Sys.time(), recipeID = recipeID,
                  recipe_file = recipe_path, qc_file = qc, checksum = item$checksum,
                  accession = recipe$accession, input_md5 = recipe$input_md5, neff = neff),
             numbers, recipe[intersect(names(recipe), c("source", "trait", "type", "category",
               "efo", "gender", "ancestry", "build", "reference", "comments",
               "study_n_total", "study_n_case", "study_n_control"))])
    for (field in names(row)) {
      column <- current$study[[field]]
      if (is.null(column)) column <- rep(row[[field]][NA_integer_], at)
      column[at] <- row[[field]]
      current$study[[field]] <- column
    }
    for (pair in list(c(item$staged, paste0(name, ".gz")), c(item$qc_file, qc))) {
      from <- file.path(stage, pair[2L])
      if (!file.copy(pair[1L], from)) stop("Cannot stage ", pair[2L], " for publication")
      source <- c(source, from)
      destination <- c(destination, file.path(gstat, pair[2L]))
    }
  }
  current$studies <- as.data.frame(current$study)
  rownames(current$studies) <- current$studies$id
  current$studyfiles <- stats::setNames(file.path(gstat, paste0(current$study$file, ".gz")), current$study$id)
  temp_csv <- file.path(stage, basename(inventory))
  write.csv2(current$studies, temp_csv, row.names = FALSE)
  saveRDS(current, temp_rds, compress = FALSE)
  suspendInterrupts({
    for (i in seq_along(source)) {
      .publish_gwas_import_file(source[i], destination[i])
      published <- c(published, destination[i])
    }
    .publish_gwas_import_file(temp_csv, inventory, replace = TRUE)
    replaced <- TRUE
    .publish_gwas_import_file(temp_rds, descriptor, replace = TRUE)
    committed <- TRUE
  })
  new_ids
}
