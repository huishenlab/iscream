# Verify that the input BED files are of the type specified by the input aligner

Verify that the input BED files are of the type specified by the input
aligner

## Usage

``` r
verify_filetype(bedfiles, aligner, stop_on_error = FALSE)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- aligner:

  The aligner chosen

- stop_on_error:

  Whether to warn or stop on aligner-filename mismatch

## Value

TRUE if all input BED files have an associated tabix index file. FALSE
if not
