\name{writeMarkerStat}
\alias{writeMarkerStat}
\alias{assocMatrix}
\alias{writeZMatrix}
\title{Export a GWAS association matrix}
\usage{
writeMarkerStat(GAlist, studyID = NULL, out_dir = NULL, stat = "z")

assocMatrix(GAlist, out_dir, Glist = NULL, stat = "z")

writeZMatrix(GAlist, out_dir, Glist = NULL)
}
\arguments{
\item{GAlist}{A gact database descriptor with marker and GWAS data installed.}
\item{studyID}{Optional character vector of unique study IDs. The default
writes all studies in GAlist order.}
\item{out_dir}{A new output directory. Defaults to \code{gwas_matrix} next
to the marker directory, with a study ID suffix for one selected study or
\code{_selected} for multiple selected studies. Its parent must already exist.}
\item{Glist}{Optional ancestry-named list of aligned Glist objects with
gmat sparse LD, for the compatibility entry points.}
\item{stat}{Statistic to export. The default, \code{"z"}, exports Z scores
and is the only supported value.}
}
\value{Invisibly, a list with the matrix and catalog paths, matrix dimensions,
and, when available, the LD catalog, row map paths and region catalog path.}
\description{
Writes little-endian IEEE FP64 values in column-major order to
\code{z.f64.bin} and creates \code{markers.tsv} and \code{studies.tsv}.
Rows follow \code{GAlist$rsids}. Columns follow \code{studyID} when supplied
or the study inventory by default.
}
\details{
\code{assocMatrix()} and \code{writeZMatrix()} are compatibility entry points
for exporting the full study inventory. They take reference panels through
their \code{Glist} argument.

The writer reads the study files directly and checks coordinates and alleles
against the gact marker catalog. It derives Z from finite \code{b} and positive
finite \code{seb} when available, with finite stored Z as a fallback.
Missing and unusable scores are stored as NaN. The study catalog records the
ancestry and the count of observed and missing scores for each column.

For the EUR, EAS and SAS ancestries represented in the export,
\code{writeMarkerStat()} looks for sparse LD in the marker directory
with prefix \code{g1000_<ancestry>.subset.aligned}, using lowercase ancestry
codes. Studies can be exported when their ancestry's panel is absent.
An incomplete set of sparse LD files triggers an error.
Use \code{addLDPanels()} to supply a custom panel for an ancestry.

All studies use the effect allele in the gact marker catalog at each matrix
row. The exporter checks that each sparse LD marker
shared with the Z matrix uses that same dosage allele. Reversed panel alleles
must be aligned in the genotype files and LD before export. For each panel,
\code{ld.tsv} identifies the LD marker file, binary map, and coverage counts.
Each \code{ld-<ancestry>.u32le.bin} contains one little-endian unsigned 32-bit
integer per LD row in disk order. Positive entries are one-based Z row
indices. Zero means the LD marker is absent from \code{GAlist$rsids}.
Gather the mapped Z scores in LD order and use the corresponding rows and
columns of the original LD. No sign adjustment is required for that slice.

With LD region files installed by \code{downloadDB(GAlist, what = "ldetect")}, the export
writes \code{ld_regions.tsv} for panels with a matching assembly and ancestry.
Each entry names an interval table, a region index and a binary row-pair file.
Region IDs are zero-based and follow the interval table. Each index record
gives the offset and count in pairs, including for empty regions. A pair
contains two little-endian unsigned 32-bit integers, the zero-based Z and LD
rows, in LD order within each region. Region membership uses inclusive starts
and exclusive stops.
}
