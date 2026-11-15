\name{addLDRegions}
\alias{addLDRegions}
\title{Register Berisa–Pickrell LD regions}
\usage{addLDRegions(GAlist, assembly = "GRCh37")}
\arguments{
\item{GAlist}{A gact database descriptor.}
\item{assembly}{Genome assembly of the marker catalog. Only GRCh37 is supported.}
}
\value{The updated descriptor with file paths and provenance in
\code{GAlist$ldRegions$GRCh37}.}
\description{Installs the bundled EUR and EAS Berisa–Pickrell intervals in the
marker directory. Starts are inclusive and stops are exclusive. EAS is labelled
ASN upstream. Use \code{downloadDB(GAlist, what = "ldetect")} for database setup.}
