% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_meth_regions.R
\name{summarize_meth_regions}
\alias{summarize_meth_regions}
\title{Summarize CpGs  methylation information over genomic regions}
\usage{
summarize_meth_regions(
  bedfiles,
  regions,
  fun = "all",
  aligner = "biscuit",
  feature_col = NULL,
  mval = TRUE,
  set_region_rownames = FALSE,
  nthreads = NULL
)
}
\arguments{
\item{bedfiles}{A vector of bedfile paths}

\item{regions}{A vector, data frame or GenomicRanges of genomic regions. See
details.}

\item{fun}{Function(s) to apply over the region. See details.}

\item{aligner}{The aligner used to produce the BED files - one of "biscuit",
"bismark", "bsbolt".}

\item{feature_col}{If the input is a dataframe, the column to use as the
feature label instead of the genomic region string}

\item{mval}{Whether to calculate the M value (coverage \eqn{\times \beta})
or use the beta value when applying the function.}

\item{set_region_rownames}{Use the region strings as the returned data
frame's rownames. Can be useful if you have a named regions and want both
the regions strings rownames and the feature names.}

\item{nthreads}{Set number of threads to use overriding the
\code{"iscream.threads"} option. See \code{?set_threads} for more information.}
}
\value{
A data.frame
}
\description{
Run summarizing functions on the CpGs in bedfiles across genomic regions.
Parallelized across files using threads from the \code{"iscream.threads"} option.
}
\details{
The input regions may be string vector in the form "chr:start-end"
or a GRanges object. If a data frame is provided, they must have "chr",
"start", and "end" columns. If the string vector and GenomicRanges inputs
are named, the names will be used to describe each feature in the output
dataframe. If input dataframes have a feature column, set \code{feature_col} to
that column name to populate the output's feature column.

Supported \code{fun} arguments are given below. For each of these functions,
setting \code{mval = FALSE} will use the beta values instead of the M value:
\itemize{
\item Sum: \code{"sum"}
\item Mean: \code{"mean"}
\item Median: \code{"median"}
\item Standard deviation: \code{"stddev"}
\item Variance: \code{"variance"}
\item Minimum: \code{"min"}
\item Maximum: \code{"max"}
\item Range: \code{"range"}
\item No. of CpGs in the region: \code{"cpg_count"}
}

The summarizing computations are backed by the Armadillo library. See
\url{https://arma.sourceforge.net/docs.html#stats_fns} for futher details on the
supported functions
}
\examples{
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
# examine the bedfiles
colnames <- c("chr", "start", "end", "beta", "coverage")
lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_meth_regions(bedfiles, regions)
names(regions) <- c("A", "B", "C")
summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
summarize_meth_regions(bedfiles, regions, fun = "sum")
}
