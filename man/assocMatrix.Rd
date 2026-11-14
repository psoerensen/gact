\name{assocMatrix}
\alias{assocMatrix}
\alias{writeZMatrix}
\title{Export a GWAS association matrix}
\usage{
assocMatrix(GAlist, out_dir, Glist = NULL, stat = "z")

writeZMatrix(GAlist, out_dir, Glist = NULL)
}
\arguments{
\item{GAlist}{A gact database descriptor with marker and GWAS data installed.}
\item{out_dir}{A new output directory. Its parent must already exist.}
\item{Glist}{Optional list of Glist objects named by ancestry, for example
\code{list(EUR = eur, EAS = eas, SAS = sas)}, with gmat sparse LD prepared.
Panels must use the same genome build, strand and allele order as the gact
catalog for shared markers.}
\item{stat}{Statistic to export. The default, \code{"z"}, exports Z scores
and is the only supported value.}
}
\value{Invisibly, a list with the matrix and catalog paths, matrix dimensions,
and, when panels are supplied, the LD catalog and row map paths.}
\description{
Writes little-endian IEEE FP64 values in column-major order to
\code{z.f64.bin} and creates \code{markers.tsv} and \code{studies.tsv}.
Rows follow \code{GAlist$rsids}. Columns follow the study inventory.
}
\details{
\code{writeZMatrix()} is the compatibility entry point for Z-score export.

The writer reads the study files directly and checks coordinates and alleles
against the gact marker catalog. It derives Z from finite \code{b} and positive
finite \code{seb} when available, with finite stored Z as a fallback.
Missing and unusable scores are stored as NaN. The study catalog records the
ancestry and the count of observed and missing scores for each column.

All studies use the effect allele in the gact marker catalog at each matrix
row. With reference panels, the exporter checks that each sparse LD marker
shared with the Z matrix uses that same dosage allele. Reversed panel alleles
must be aligned in the genotype files and LD before export. For each panel,
\code{ld.tsv} identifies the LD marker file, binary map, and coverage counts.
Each \code{ld-<ancestry>.u32le.bin} contains one little-endian unsigned 32-bit
integer per LD row in disk order. Positive entries are one-based Z row
indices. Zero means the LD marker is absent from \code{GAlist$rsids}.
Gather the mapped Z scores in LD order and use the corresponding rows and
columns of the original LD. No sign adjustment is required for that slice.
}
