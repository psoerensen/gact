\name{install_gsuite}
\alias{install_gsuite}
\title{Install the standalone gsuite analysis packages}
\usage{
install_gsuite(packages = c("gbase", "gbayes", "gcorr", "glma",
  "gmap", "gscore", "gsea"), lib = .libPaths()[1L],
  repos = getOption("repos"), ...)
}
\arguments{
\item{packages}{Unique standalone package names; defaults to gbase and all six analysis packages.}
\item{lib}{Destination R library.}
\item{repos}{Additional repositories. The gsuite source repository is included
and the CRAN placeholder is resolved to cloud.r-project.org.}
\item{...}{Arguments passed to \code{utils::install.packages()}, such as type.}
}
\description{Explicitly installs the standalone analysis packages. This helper
is never run when gact is loaded and does not download the database.}
\details{Source installation requires a configured C++17 toolchain, including
Rtools matching R on Windows. Checks that the requested package versions were
installed in the chosen library. Already current packages are reused.}
\value{Invisibly, the requested package names after installation succeeds.}
\examples{
\dontrun{
install_gsuite()
install_gsuite(packages = c("gbase", "gcorr"))
}
}
