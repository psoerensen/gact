\name{addLDRegions}
\alias{addLDRegions}
\title{Register Berisa–Pickrell LD regions}
\usage{addLDRegions(GAlist, assembly = "GRCh37")}
\arguments{
\item{GAlist}{A gact database descriptor.}
\item{assembly}{Genome assembly of the marker catalog. Only GRCh37 is supported.}
}
\value{The updated descriptor with provenance and interval tables in
\code{GAlist$ldRegions$GRCh37}.}
\description{Registers the bundled EUR and EAS Berisa–Pickrell intervals.
Starts are inclusive and stops are exclusive. EAS is labelled ASN upstream.}
