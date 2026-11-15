\name{install_gsuite}
\alias{install_gsuite}
\title{Install the standalone gsuite analysis packages}
\usage{
install_gsuite(packages = c("gbase", "gbayes", "gcorr", "glma",
  "gmap", "gscore", "gsea"), lib = .libPaths()[1L],
  repos = getOption("repos"), ..., type = "source", available = NULL)
}
\arguments{
\item{packages}{Unique standalone package names; defaults to gbase and all six analysis packages.}
\item{lib}{Destination R library.}
\item{repos}{Additional repositories. The gsuite source repository is included
and the CRAN placeholder is resolved to cloud.r-project.org.}
\item{...}{Arguments passed to \code{utils::install.packages()}.}
\item{type}{Installation type. The default uses the published source packages.}
\item{available}{Optional package index to pass to \code{utils::install.packages()}.
Published gsuite source versions are always used to verify the result.}
}
\description{Explicitly installs the standalone analysis packages. This helper
is never run when gact is loaded and does not download the database.}
\details{Source installation requires a configured C++17 toolchain, including
Rtools matching R on Windows. The installer refreshes the gsuite source index,
installs missing or older packages, and checks the versions in the destination
library after installing. Installed versions at least as recent as the latest
compatible published versions are reused.}
\value{Invisibly, the requested package names after installation succeeds.}
\examples{
\dontrun{
install_gsuite()
install_gsuite(packages = c("gbase", "gcorr"))
}
}
