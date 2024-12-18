% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_gact.R
\name{getLDscoresDB}
\alias{getLDscoresDB}
\title{Retrieve Linkage Disequilibrium (LD) Scores for Specific Genetic Markers}
\usage{
getLDscoresDB(
  GAlist = NULL,
  chr = NULL,
  ancestry = "EUR",
  version = "HapMap3",
  rsids = NULL
)
}
\arguments{
\item{GAlist}{A list containing directory paths where marker data is stored. It is assumed that \code{GAlist$dirs["marker"]} points to the correct directory.}

\item{ancestry}{A character string specifying the ancestry population. The options are "EUR" (European), "EAS" (East Asian), or "SAS" (South Asian). Default is "EUR".}

\item{version}{A character string specifying the version of the reference panel. The options are "HapMap3", "1000G", or "Original". Default is "HapMap3".}

\item{rsids}{An optional character vector of rsIDs to filter the LD scores. If not specified, all LD scores will be returned.}
}
\value{
A named vector of LD scores, where the names are rsIDs.
}
\description{
This function retrieves LD scores for specified genetic markers based on the population ancestry and version of the reference panel. The function reads in marker files and extracts LD scores, allowing for optional filtering by specific rsIDs.
}
\details{
The function reads in precomputed LD score files for the selected ancestry and version, extracts the LD scores, and returns them as a named vector. If the \code{rsids} argument is provided, the returned vector will only include LD scores for the specified rsIDs.
}
\examples{
# Example usage:
GAlist <- list(dirs = c(marker = "/path/to/marker/files"))
ldscores <- getLDscoresDB(GAlist, ancestry = "EUR", version = "HapMap3", rsids = c("rs123", "rs456"))

}
