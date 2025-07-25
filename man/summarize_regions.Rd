% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_regions.R
\name{summarize_regions}
\alias{summarize_regions}
\title{Summarize information over genomic regions from any BED file}
\usage{
summarize_regions(
  bedfiles,
  regions,
  columns,
  col_names = NULL,
  fun = "all",
  feature_col = NULL,
  set_region_rownames = FALSE,
  nthreads = NULL
)
}
\arguments{
\item{bedfiles}{A vector of BED file paths}

\item{regions}{A vector, data frame or GenomicRanges of genomic regions. See
details.}

\item{columns}{A vector of indices of the numeric columns to be summarized}

\item{col_names}{A vector of names to use for \code{columns} in the output}

\item{fun}{Function(s) to apply over the region. See details.}

\item{feature_col}{Column name of the input \code{regions} data frame containing
a name for each genomic region. Set only if the using a data frame as the
input regions format. See details.}

\item{set_region_rownames}{Use the region strings as the returned data
frame's rownames. Can be useful if you have a named regions and want both
the regions strings rownames and the feature names. See details.}

\item{nthreads}{Set number of threads to use overriding the
\code{"iscream.threads"} option. See \code{?set_threads} for more information.}
}
\value{
A data.frame
}
\description{
Run summarizing functions on BED file records across genomic regions.
Parallelized across files using threads from the \code{"iscream.threads"} option.
}
\section{Supported functions}{
\itemize{
\item Sum: \code{"sum"}
\item Mean: \code{"mean"}
\item Median: \code{"median"}
\item Standard deviation: \code{"stddev"}
\item Variance: \code{"variance"}
\item Minimum: \code{"min"}
\item Maximum: \code{"max"}
\item Range: \code{"range"}
\item No. of records in the region: \code{"count"}
}

The summarizing computations are backed by the Armadillo library. See
\url{https://arma.sourceforge.net/docs.html#stats_fns} for futher details on the
supported functions
}

\section{Using feature identifiers}{
\code{regions} may be string vector in the form "chr:start-end", a GRanges
object or a data frame with "chr", "start", and "end" columns. The \code{feature}
column of the output will contain a "chr:start-end" identifier for each
summarized region. To use other identifiers, like a gene name for a region
instead of the coordinates, set the names of the vector or GRanges to those
identifiers. These names will be used instead of the genomic region string
to describe each feature in the output dataframe. If \code{regions} is a data
frame make an additional column with the identifiers and pass that column
name to \code{feature_col}. See examples.
}

\examples{
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
# examine the bedfiles
colnames <- c("chr", "start", "end", "beta", "coverage")
lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_regions(bedfiles, regions, columns = c(4, 5), col_names = c("beta", "cov"))

# select functions
summarize_regions(
  bedfiles,
  regions,
  fun = c("mean", "stddev"),
  columns = c(4, 5),
  col_names = c("beta", "cov")
)

# add names to the regions
names(regions) <- c("A", "B", "C")
summarize_regions(
  bedfiles,
  regions,
  fun = "sum",
  columns = 5,
  col_names = "coverage"
)

# using `feature_col`
library(data.table)

# convert string vector to a data.table
regions_df <- data.table::as.data.table(regions) |>
_[, tstrsplit(regions, ":|-", fixed = FALSE, names = c("chr", "start", "end"))] |>
_[, start := as.integer(start)] |>
_[, feature := LETTERS[.I]][]
regions_df

summarize_regions(
  bedfiles,
  regions_df,
  fun = "sum",
  columns = 5,
  col_names = "coverage",
  feature_col = "feature"
)
}
