% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_mat.R
\name{make_mat}
\alias{make_mat}
\title{Make a matrix from a numeric column of any BED file}
\usage{
make_mat(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
)
}
\arguments{
\item{bedfiles}{A vector of bedfile paths}

\item{regions}{A vector, data frame or GenomicRanges of genomic regions. See
details.}

\item{column}{The index of the data column needed for the matrix}

\item{mat_name}{What to name the matrix in the returned list}

\item{sparse}{Whether to return a sparse matrix}

\item{prealloc}{The number of rows to initialize the matrices with. If the
number of loci are approximately known, this can reduce runtime as fewer
resizes need to be made.}

\item{nthreads}{Set number of threads to use overriding the
\code{"iscream.threads"} option. See \code{?set_threads} for more information.}
}
\value{
A named list of
\itemize{
\item the matrix with the value of interest
\item a character vector of chromosomes and numeric vector of base positions
\item a character vector of the input sample BED file names
}
}
\description{
Queries the provided regions and produces a matrix along with genomic
positions.
Parallelized across files using threads from the \code{"iscream.threads"} option.
}
\details{
The input regions may be string vector in the form "chr:start-end"
or a GRanges object. If a data frame is provided, they must have "chr",
"start", and "end" columns.
}
\examples{
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
# examine the bedfiles
colnames <- c("chr", "start", "end", "beta", "coverage")
lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
# make matrix of beta values
make_mat(bedfiles, regions, column = 4)
}
