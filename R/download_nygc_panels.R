.gact_nygc_release_hashes <- function(record) {
  manifest <- system.file("extdata", paste0("nygc-", record, "-sha256.txt"),
                          package = "gact", mustWork = TRUE)
  checksums <- utils::read.table(manifest, colClasses = "character")
  stats::setNames(checksums[[1L]], checksums[[2L]])
}

.gact_sha256_file <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  as.character(openssl::sha256(connection))
}

.gact_download_nygc_archive <- function(record, archive, dest) {
  url <- paste0("https://zenodo.org/api/records/", record, "/files/", archive, "/content")
  curl::curl_download(url, dest, mode = "wb",
                      handle = curl::new_handle(connecttimeout = 30, timeout = 3000))
}

.gact_install_nygc_panels <- function(dest, doi) {
  if (length(dest) != 1L || !dir.exists(dest)) {
    stop("The GRCh38 GAlist must have an existing marker directory.")
  }
  record <- sub("^10\\.5281/zenodo\\.", "", doi)
  sha256 <- .gact_nygc_release_hashes(record)
  ancestries <- c("afr", "amr", "eas", "eur", "sas")
  triplets <- outer(c(".bed", ".bim", ".fam"), paste0("g1000_", ancestries),
                    function(extension, prefix) paste0(prefix, extension))
  filenames <- as.vector(triplets)
  targets <- file.path(dest, filenames)

  stage <- file.path(dest, ".nygc-install.lock")
  if (!dir.create(stage, showWarnings = FALSE)) stop("Cannot lock NYGC installation: ", stage)
  on.exit(unlink(stage, recursive = TRUE, expand = FALSE), add = TRUE)

  existing <- file.exists(targets)
  links <- Sys.readlink(targets)
  if (any(!is.na(links) & nzchar(links))) stop("NYGC destination contains a symbolic link.")
  for (i in which(existing)) {
    if (dir.exists(targets[i]) ||
        .gact_sha256_file(targets[i]) != sha256[[filenames[i]]]) {
      stop("Existing NYGC destination differs from the published release: ", targets[i])
    }
  }
  if (all(existing)) return(invisible(targets))

  for (i in which(colSums(matrix(existing, nrow = 3L)) != 3L)) {
    archive <- paste0("gact-nygc30x-grch38-", ancestries[i], ".tar.gz")
    archive_path <- file.path(stage, archive)
    message("Downloading GRCh38 NYGC panel: ", toupper(ancestries[i]))
    .gact_download_nygc_archive(record, archive, archive_path)
    if (.gact_sha256_file(archive_path) != sha256[[archive]]) {
      stop("SHA-256 mismatch for ", archive)
    }
    members <- triplets[, i]
    status <- utils::untar(archive_path, files = members, exdir = stage, tar = "internal")
    if (status != 0L) stop("Cannot unpack ", archive)
    for (member in members) {
      path <- file.path(stage, member)
      if (!file.exists(path) || .gact_sha256_file(path) != sha256[[member]]) {
        stop("SHA-256 mismatch for extracted file ", member)
      }
    }
    unlink(archive_path)
  }

  published <- character()
  committed <- FALSE
  on.exit(if (!committed) unlink(published, expand = FALSE), add = TRUE, after = FALSE)
  for (i in which(!existing)) {
    if (file.exists(targets[i]) || !file.rename(file.path(stage, filenames[i]), targets[i])) {
      stop("Cannot install NYGC file: ", targets[i])
    }
    published <- c(published, targets[i])
  }
  committed <- TRUE
  invisible(targets)
}
