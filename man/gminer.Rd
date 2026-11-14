\name{gminer}
\alias{gminer}
\title{Run gminer on a prepared association matrix}
\usage{gminer(GAlist, task = "impute", matrix_dir = NULL, out_dir = NULL)}
\arguments{
\item{GAlist}{A gact descriptor.}
\item{task}{Native task. Currently only \code{"impute"} is supported.}
\item{matrix_dir}{Directory written by \code{\link{assocMatrix}} with LD maps.
Defaults to \code{gwas_matrix} next to the gact marker directory.}
\item{out_dir}{New output directory. Defaults to \code{imputed} within
\code{matrix_dir}.}
}
\value{Invisibly, paths to the filled Z matrix, mask, predicted r-squared
matrix, per-study counts and task metadata.}
\description{Dispatches the prepared association matrix to the
\code{gminer_impute} executable on PATH. The native process reads the
observed Z scores and matching ancestry LD directly.}
\details{Binary outputs follow the original marker and study order. Missing
scores without matching LD or a prediction passing the quality threshold
remain missing.}
