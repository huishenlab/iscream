# Verify that BED files are tabixed

Verify that BED files are tabixed

## Usage

``` r
verify_files_or_stop(bedfiles, verify_tabix = TRUE)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- verify_tabix:

  Whether to verify the presence of tabix files

## Value

TRUE if all input BED files have an associated tabix index file. FALSE
if not
