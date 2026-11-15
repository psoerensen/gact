\name{addLDPanels}
\alias{addLDPanels}
\title{Register custom sparse LD panels for marker statistics}
\usage{addLDPanels(GAlist, Glist)}
\arguments{
\item{GAlist}{A gact database descriptor.}
\item{Glist}{An ancestry-named list of aligned Glist objects with sparse LD,
for example \code{list(EUR = eur, EAS = eas, SAS = sas)}.}
}
\value{The updated descriptor with reference paths in \code{GAlist$ldPanels}.}
\description{Records the path to each panel's gmat sparse LD marker file and
the panel's assembly when supplied.}
\details{\code{writeMarkerStat()} uses each registered panel in preference
to the conventional filename for that ancestry. Register a new list to
replace the current references.}
