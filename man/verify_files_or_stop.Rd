% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verify_files.R
\name{verify_files_or_stop}
\alias{verify_files_or_stop}
\title{Verify that BED files are tabixed}
\usage{
verify_files_or_stop(bedfiles, verify_tabix = TRUE)
}
\arguments{
\item{bedfiles}{A vector of BED file paths}

\item{verify_tabix}{Whether to verify the presence of tabix files}
}
\value{
TRUE if all input BED files have an associated tabix index file.
FALSE if not
}
\description{
Verify that BED files are tabixed
}
\keyword{internal}
