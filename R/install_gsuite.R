#' Install the standalone gsuite analysis packages
#'
#' Installs the analysis packages explicitly, separately from downloading the
#' gact database. This function is never run when gact is loaded.
#' @param packages Packages to install; defaults to gbase and all six analysis packages.
#' @param lib Destination R library, defaulting to the first library in .libPaths().
#' @param repos Additional R repositories. The gsuite source repository is always
#'   included; the standard CRAN placeholder is replaced with cloud.r-project.org.
#' @param ... Additional arguments to utils::install.packages, such as type.
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
    "gmap", "gscore", "gsea"), lib = .libPaths()[1L], repos = getOption("repos"), ...) {
  supported <- c("gbase", "gbayes", "gcorr", "glma", "gmap", "gscore", "gsea")
  if (!is.character(packages) || !length(packages) || anyNA(packages) ||
      any(!packages %in% supported) || anyDuplicated(packages))
    stop("Choose unique standalone gsuite package names.", call. = FALSE)
  if (!is.character(repos) || anyNA(repos) || any(!nzchar(repos)))
    stop("repos must contain R repository URLs.", call. = FALSE)
  repos[repos == "@CRAN@"] <- "https://cloud.r-project.org"
  repos <- c(gsuite = "https://psoerensen.github.io/gtools/gsuite/repository", repos)
  repos <- repos[!duplicated(unname(repos))]
  if (!length(repos[grepl("cran|r-project", repos, ignore.case = TRUE)]))
    repos <- c(repos, CRAN = "https://cloud.r-project.org")
  current <- function(package) {
    nzchar(system.file(package = package, lib.loc = lib)) &&
      utils::packageVersion(package, lib.loc = lib) >=
        package_version(if (package == "gbase") "0.1.3" else if (package == "gcorr") "0.1.9" else if (package == "glma") "0.1.2" else "0.1.1")
  }
  needed <- packages[!vapply(packages, current, logical(1))]
  if (length(needed)) utils::install.packages(needed, lib = lib, repos = repos, ...)
  missing <- packages[!vapply(packages, current, logical(1))]
  if (length(missing)) stop("Installation did not produce: ",
    paste(missing, collapse = ", "), ". Inspect the installation output.", call. = FALSE)
  invisible(packages)
}
