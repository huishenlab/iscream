# DataFrame to region strings

Convert DataFrame to a vector of strings. Set feature names in a "name"
column

## Usage

``` r
get_df_string(regions_df, feature_col = NULL)
```

## Arguments

- regions_df:

  A data frame with "chr", "start" and "end" columns

- feature_col:

  The data frame column to use as the names of the output string vector

## Value

A character vector

## Examples

``` r
(df <- data.frame(chr = c("chr1", "chr2"), start = c(1, 5), end = c(4, 10)))
#>    chr start end
#> 1 chr1     1   4
#> 2 chr2     5  10
get_df_string(df)
#> [1] "chr1:1-4"  "chr2:5-10"
```
