#' Install the standalone gsuite analysis packages
#'
#' Installs the analysis packages explicitly, separately from downloading the
#' gact database. This function is never run when gact is loaded.
#' It checks the current gsuite source index and installs missing or older
#' packages in the destination library.
#' @param packages Packages to install; defaults to gbase and all six analysis packages.
#' @param lib Destination R library, defaulting to the first library in .libPaths().
#' @param repos Additional R repositories. The gsuite source repository is always
#'   included; the standard CRAN placeholder is replaced with cloud.r-project.org.
#' @param ... Additional arguments to utils::install.packages.
#' @param type Installation type. The default uses the published source packages.
#' @param available Optional package index to pass to utils::install.packages.
#'   Published gsuite source versions are always used to verify the result.
#' @return Invisibly, the requested package names after installation succeeds.
#' @details Source installation requires R's configured C++17 toolchain (Rtools
#'   matching R on Windows). Installing packages does not download the gact database.
#' @examples
#' \dontrun{
#' install_gsuite()
#' install_gsuite(packages = c("gbase", "gcorr"))
#' }
#' @export
install_gsuite <- function(packages = c("gbase", "gbayes", "gcorr", "glma",
    "gmap", "gscore", "gsea"), lib = .libPaths()[1L], repos = getOption("repos"),
    ..., type = "source", available = NULL) {
  supported <- c("gbase", "gbayes", "gcorr", "glma", "gmap", "gscore", "gsea")
  if (!is.character(packages) || !length(packages) || anyNA(packages) ||
      any(!packages %in% supported) || anyDuplicated(packages))
    stop("Choose unique standalone gsuite package names.", call. = FALSE)
  if (!is.character(repos) || anyNA(repos) || any(!nzchar(repos)))
    stop("repos must contain R repository URLs.", call. = FALSE)
  gsuite_repo <- "https://psoerensen.github.io/gtools/gsuite/repository"
  repos[repos == "@CRAN@"] <- "https://cloud.r-project.org"
  repos <- c(gsuite = gsuite_repo, repos)
  repos <- repos[!duplicated(unname(repos))]
  if (!length(repos[grepl("cran|r-project", repos, ignore.case = TRUE)]))
    repos <- c(repos, CRAN = "https://cloud.r-project.org")
  published <- utils::available.packages(repos = gsuite_repo, type = "source",
                                        ignore_repo_cache = TRUE)
  unavailable <- setdiff(packages, rownames(published))
  if (length(unavailable))
    stop("No compatible gsuite source version available for: ",
         paste(unavailable, collapse = ", "), ". Check the repository and your R version.",
         call. = FALSE)
  target <- stats::setNames(published[packages, "Version"], packages)
  pending <- function() {
    installed <- utils::installed.packages(lib.loc = lib, noCache = TRUE)
    version <- installed[match(packages, installed[, "Package"]), "Version"]
    missing <- is.na(version)
    missing[!missing] <- package_version(version[!missing]) < package_version(target[!missing])
    packages[missing]
  }
  needed <- pending()
  if (!length(needed)) return(invisible(packages))
  if (is.null(available) && identical(type, "source")) {
    extra <- utils::available.packages(repos = repos[-1L], type = "source")
    external <- !rownames(extra) %in% union(supported, rownames(published))
    available <- rbind(published, extra[external, , drop = FALSE])
  }
  utils::install.packages(needed, lib = lib, repos = repos, type = type,
                          available = available, ...)
  stale <- pending()
  if (length(stale))
    stop("Installation did not reach the published gsuite versions for: ",
         paste(paste0(stale, " (", target[stale], ")"), collapse = ", "),
         ". Inspect the installation output.", call. = FALSE)
  invisible(packages)
}
