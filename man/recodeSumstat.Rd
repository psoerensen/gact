\name{recodeSumstat}
\alias{recodeSumstat}
\title{Align a summary-statistics file to the gact marker catalog}
\usage{recodeSumstat(file, GAlist, out_file = file)}
\arguments{
\item{file}{An existing gact-format CSV or TSV file, plain or gzipped.}
\item{GAlist}{A gact descriptor pointing to the marker catalog.}
\item{out_file}{Defaults to replacing the source after validation. Pass a
new path to keep the source. A path ending in \code{.gz} writes gzip.}
}
\value{Invisibly, the path, the numbers of reversed markers and repaired
Z fields, and the output MD5. For a registered study replaced in place,
the caller should update its recorded checksum.}
\description{
Matches variants by rsid, chromosome, position and allele pair. For a
reversal, exchanges the alleles, negates beta and Z, and complements the
effect-allele frequency. Repairs legacy Z fields with the opposite sign
and the same magnitude as beta divided by SE. Standard errors, p-values,
other fields and row order are preserved.
}
