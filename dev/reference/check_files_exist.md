# Check that files exist

Check that files exist

## Usage

``` r
check_files_exist(files_vec, error_file_prefix = "Bedfile")
```

## Arguments

- files_vec:

  A vector of file paths

- error_file_prefix:

  Error message prefix for 'Bedfile' vs 'Tabix file'

## Value

TRUE if all input BED files have an associated tabix index file. FALSE
if not
