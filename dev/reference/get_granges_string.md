# GRanges to region strings

Coerces GenomicRanges to `chr:start-end` strings with `as.character`. If
any regions have the same start and end, `as.character` returns
`chr:start` strings which are invalid for the htslib API. These are
corrected to `chr:start-start`.

## Usage

``` r
get_granges_string(gr)
```

## Arguments

- gr:

  A GRanges object

## Value

A character vector

## Examples

``` r
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  get_granges_string(GenomicRanges::GRanges(c("chr1:1-10", "chr2:15-20")))
}
#> [1] "chr1:1-10"  "chr2:15-20"
```
