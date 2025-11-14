# Query the chromosomes or seqnames from a vector of BED files

Query the chromosomes or seqnames from a vector of BED files

## Usage

``` r
query_chroms(bedfiles, nthreads = NULL)
```

## Arguments

- bedfiles:

  The vector of BED file paths

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

## Value

A vector of seqnames

## Examples

``` r
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
query_chroms(bedfiles)
#> [1] "chr1"
```
