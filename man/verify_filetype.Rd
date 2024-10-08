% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verify_files.R
\name{verify_filetype}
\alias{verify_filetype}
\title{Verify that the input bedfiles are of the type specified by the input aligner}
\usage{
verify_filetype(bedfiles, aligner)
}
\arguments{
\item{bedfiles}{A vector of bedfile paths}

\item{aligner}{The aligner chosen}
}
\value{
TRUE if all input bedfiles have an associated tabix index file.
FALSE if not
}
\description{
Verify that the input bedfiles are of the type specified by the input aligner
}
\keyword{internal}
