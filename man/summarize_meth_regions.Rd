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
\item{bedfiles}{A vector of BED file paths}

\item{regions}{A vector, data frame or GenomicRanges of genomic regions. See
details.}

\item{fun}{Function(s) to apply over the region. See details.}

\item{aligner}{The aligner used to produce the BED files - one of "biscuit",
"bismark", "bsbolt".}

\item{feature_col}{Column name of the input \code{regions} data frame containing
a name for each genomic region. Set only if the using a data frame as the
input regions format. See details.}

\item{mval}{Whether to calculate the M value (coverage \eqn{\times \beta})
or use the beta value when applying the function.}

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
Run summarizing functions on the CpGs in BED files across genomic regions.
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
# also see examples from ?summarize_regions

bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_meth_regions(bedfiles, regions)
names(regions) <- c("A", "B", "C")
summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
summarize_meth_regions(bedfiles, regions, fun = "sum")
}
