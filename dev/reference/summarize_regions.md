# Summarize information over genomic regions from any BED file

Run summarizing functions on BED file records across genomic regions.
Parallelized across files using threads from the `"iscream.threads"`
option.

## Usage

``` r
summarize_regions(
  bedfiles,
  regions,
  columns,
  col_names = NULL,
  fun = "all",
  feature_col = NULL,
  nthreads = NULL
)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- regions:

  A vector, data frame or GenomicRanges of genomic regions. See details.

- columns:

  A vector of indices of the numeric columns to be summarized

- col_names:

  A vector of names to use for `columns` in the output

- fun:

  Function(s) to apply over the region. See details.

- feature_col:

  Column name of the input `regions` data frame or GRanges `mcols`
  containing a name for each genomic region. Set only if the using a
  data frame-like object as the input regions format, not for string
  vectors. See details.

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

## Value

A data.table

## Supported functions

- Sum: `"sum"`

- Mean: `"mean"`

- Median: `"median"`

- Mode: `"mode"`

- Anti-mode: `"antimode"`

- Sample standard deviation: `"stddev"` (`sstdev` in `bedtools map`)

- Population standard deviation: `"pstddev"` (`stdev` in `bedtools map`)

- Variance: `"variance"`

- Minimum: `"min"`

- Maximum: `"max"`

- Minimum of absolute values: `"absmin"`

- Maximum of absolute values: `"absmax"`

- Range: `"range"`

- First element: `"first"`

- Last element: `"last"`

- No. of records in the region: `"count"`

- No. of records in the region with unique data values:
  `"count_distinct"`

Most summarizing computations are backed by the Armadillo library. See
<https://arma.sourceforge.net/docs.html#stats_fns> for futher details on
the supported functions

## Using feature identifiers

`regions` may be a string vector in the form "chr:start-end", a GRanges
object or a data frame with "chr", "start", and "end" columns. If the
input data.frame or GRanges has a column with feature identifiers, like
gene names for a set of gene regions, pass that column's name to
`feature_col`. If `regions` is a vector, set its
[`names()`](https://rdrr.io/r/base/names.html) to those identifiers.
These will be used to populate a 'feature' column in the summary. See
examples.

## Examples

``` r
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
# examine the bedfiles
colnames <- c("chr", "start", "end", "beta", "coverage")
lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))
#> [[1]]
#> 
#> 
#> |chr  | start| end| beta| coverage|
#> |:----|-----:|---:|----:|--------:|
#> |chr1 |     0|   2|  1.0|        1|
#> |chr1 |     2|   4|  1.0|        1|
#> |chr1 |     4|   6|  0.0|        2|
#> |chr1 |     6|   8|  0.0|        1|
#> |chr1 |     8|  10|  0.5|        2|
#> |chr1 |    10|  12|  1.0|        2|
#> |chr1 |    12|  14|  1.0|        3|
#> 
#> [[2]]
#> 
#> 
#> |chr  | start| end| beta| coverage|
#> |:----|-----:|---:|----:|--------:|
#> |chr1 |     0|   2|    0|        2|
#> |chr1 |     4|   6|    1|        2|
#> |chr1 |     6|   8|    1|        1|
#> |chr1 |    10|  12|    0|        2|
#> |chr1 |    12|  14|    1|        1|
#> 
#> [[3]]
#> 
#> 
#> |chr  | start| end| beta| coverage|
#> |:----|-----:|---:|----:|--------:|
#> |chr1 |     2|   4|    1|        2|
#> |chr1 |     6|   8|    0|        2|
#> |chr1 |     8|  10|    1|        1|
#> 
#> [[4]]
#> 
#> 
#> |chr  | start| end| beta| coverage|
#> |:----|-----:|---:|----:|--------:|
#> |chr1 |     0|   2|  1.0|        1|
#> |chr1 |     2|   4|  1.0|        2|
#> |chr1 |     6|   8|  0.0|        1|
#> |chr1 |     8|  10|  0.5|        2|
#> |chr1 |    12|  14|  1.0|        1|
#> 

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_regions(bedfiles, regions, columns = c(4, 5), col_names = c("beta", "cov"))
#> [15:32:59.427266] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:32:59.427282] [iscream::summarize_regions] [info] using sum, mean, median, mode, antimode, stddev, pstddev, variance, min, max, absmin, absmax, range, first, last, count_distinct, count
#> [15:32:59.427286] [iscream::summarize_regions] [info] with columns 4, 5 as beta, cov
#>        chr start   end   file beta.sum cov.sum beta.mean cov.mean beta.median
#>     <char> <int> <int> <char>    <num>   <num>     <num>    <num>       <num>
#>  1:   chr1     1     6      a      2.0       4 0.6666667 1.333333        1.00
#>  2:   chr1     7    10      a      0.5       3 0.2500000 1.500000        0.25
#>  3:   chr1    11    14      a      2.0       5 1.0000000 2.500000        1.00
#>  4:   chr1     1     6      b      1.0       4 0.5000000 2.000000        0.50
#>  5:   chr1     7    10      b      1.0       1 1.0000000 1.000000        1.00
#>  6:   chr1    11    14      b      1.0       3 0.5000000 1.500000        0.50
#>  7:   chr1     1     6      c      1.0       2 1.0000000 2.000000        1.00
#>  8:   chr1     7    10      c      1.0       3 0.5000000 1.500000        0.50
#>  9:   chr1    11    14      c       NA      NA        NA       NA          NA
#> 10:   chr1     1     6      d      2.0       3 1.0000000 1.500000        1.00
#> 11:   chr1     7    10      d      0.5       3 0.2500000 1.500000        0.25
#> 12:   chr1    11    14      d      1.0       1 1.0000000 1.000000        1.00
#>     cov.median beta.mode cov.mode beta.antimode cov.antimode beta.stddev
#>          <num>     <num>    <num>         <num>        <num>       <num>
#>  1:        1.0         1        1             0            2   0.5773503
#>  2:        1.5         0        1             0            1   0.3535534
#>  3:        2.5         1        2             1            2   0.0000000
#>  4:        2.0         0        2             0            2   0.7071068
#>  5:        1.0         1        1             1            1   0.0000000
#>  6:        1.5         0        1             0            1   0.7071068
#>  7:        2.0         1        2             1            2   0.0000000
#>  8:        1.5         0        1             0            1   0.7071068
#>  9:         NA        NA       NA            NA           NA          NA
#> 10:        1.5         1        1             1            1   0.0000000
#> 11:        1.5         0        1             0            1   0.3535534
#> 12:        1.0         1        1             1            1   0.0000000
#>     cov.stddev beta.pstddev cov.pstddev beta.variance cov.variance beta.min
#>          <num>        <num>       <num>         <num>        <num>    <num>
#>  1:  0.5773503    0.4714045   0.4714045     0.3333333    0.3333333        0
#>  2:  0.7071068    0.2500000   0.5000000     0.1250000    0.5000000        0
#>  3:  0.7071068    0.0000000   0.5000000     0.0000000    0.5000000        1
#>  4:  0.0000000    0.5000000   0.0000000     0.5000000    0.0000000        0
#>  5:  0.0000000    0.0000000   0.0000000     0.0000000    0.0000000        1
#>  6:  0.7071068    0.5000000   0.5000000     0.5000000    0.5000000        0
#>  7:  0.0000000    0.0000000   0.0000000     0.0000000    0.0000000        1
#>  8:  0.7071068    0.5000000   0.5000000     0.5000000    0.5000000        0
#>  9:         NA           NA          NA            NA           NA       NA
#> 10:  0.7071068    0.0000000   0.5000000     0.0000000    0.5000000        1
#> 11:  0.7071068    0.2500000   0.5000000     0.1250000    0.5000000        0
#> 12:  0.0000000    0.0000000   0.0000000     0.0000000    0.0000000        1
#>     cov.min beta.max cov.max beta.absmin cov.absmin beta.absmax cov.absmax
#>       <num>    <num>   <num>       <num>      <num>       <num>      <num>
#>  1:       1      1.0       2           0          1         1.0          2
#>  2:       1      0.5       2           0          1         0.5          2
#>  3:       2      1.0       3           1          2         1.0          3
#>  4:       2      1.0       2           0          2         1.0          2
#>  5:       1      1.0       1           1          1         1.0          1
#>  6:       1      1.0       2           0          1         1.0          2
#>  7:       2      1.0       2           1          2         1.0          2
#>  8:       1      1.0       2           0          1         1.0          2
#>  9:      NA       NA      NA          NA         NA          NA         NA
#> 10:       1      1.0       2           1          1         1.0          2
#> 11:       1      0.5       2           0          1         0.5          2
#> 12:       1      1.0       1           1          1         1.0          1
#>     beta.range cov.range beta.first cov.first beta.last cov.last
#>          <num>     <num>      <num>     <num>     <num>    <num>
#>  1:        1.0         1          1         1       0.0        2
#>  2:        0.5         1          0         1       0.5        2
#>  3:        0.0         1          1         2       1.0        3
#>  4:        1.0         0          0         2       1.0        2
#>  5:        0.0         0          1         1       1.0        1
#>  6:        1.0         1          0         2       1.0        1
#>  7:        0.0         0          1         2       1.0        2
#>  8:        1.0         1          0         2       1.0        1
#>  9:         NA        NA         NA        NA        NA       NA
#> 10:        0.0         1          1         1       1.0        2
#> 11:        0.5         1          0         1       0.5        2
#> 12:        0.0         0          1         1       1.0        1
#>     beta.count_distinct cov.count_distinct count
#>                   <num>              <num> <num>
#>  1:                   2                  2     3
#>  2:                   2                  2     2
#>  3:                   1                  2     2
#>  4:                   2                  1     2
#>  5:                   1                  1     1
#>  6:                   2                  2     2
#>  7:                   1                  1     1
#>  8:                   2                  2     2
#>  9:                  NA                 NA    NA
#> 10:                   1                  2     2
#> 11:                   2                  2     2
#> 12:                   1                  1     1

# select functions
summarize_regions(
  bedfiles,
  regions,
  fun = c("mean", "stddev"),
  columns = c(4, 5),
  col_names = c("beta", "cov")
)
#> [15:32:59.467846] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:32:59.467861] [iscream::summarize_regions] [info] using mean, stddev
#> [15:32:59.467865] [iscream::summarize_regions] [info] with columns 4, 5 as beta, cov
#>        chr start   end   file beta.mean cov.mean beta.stddev cov.stddev
#>     <char> <int> <int> <char>     <num>    <num>       <num>      <num>
#>  1:   chr1     1     6      a 0.6666667 1.333333   0.5773503  0.5773503
#>  2:   chr1     7    10      a 0.2500000 1.500000   0.3535534  0.7071068
#>  3:   chr1    11    14      a 1.0000000 2.500000   0.0000000  0.7071068
#>  4:   chr1     1     6      b 0.5000000 2.000000   0.7071068  0.0000000
#>  5:   chr1     7    10      b 1.0000000 1.000000   0.0000000  0.0000000
#>  6:   chr1    11    14      b 0.5000000 1.500000   0.7071068  0.7071068
#>  7:   chr1     1     6      c 1.0000000 2.000000   0.0000000  0.0000000
#>  8:   chr1     7    10      c 0.5000000 1.500000   0.7071068  0.7071068
#>  9:   chr1    11    14      c        NA       NA          NA         NA
#> 10:   chr1     1     6      d 1.0000000 1.500000   0.0000000  0.7071068
#> 11:   chr1     7    10      d 0.2500000 1.500000   0.3535534  0.7071068
#> 12:   chr1    11    14      d 1.0000000 1.000000   0.0000000  0.0000000

# add names to the regions
names(regions) <- c("gene1", "gene2", "gene3")
summarize_regions(
  bedfiles,
  regions,
  fun = "sum",
  columns = 5,
  col_names = "coverage"
)
#> [15:32:59.475690] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:32:59.475704] [iscream::summarize_regions] [info] using sum
#> [15:32:59.475708] [iscream::summarize_regions] [info] with columns 5 as coverage
#>        chr start   end   file feature coverage.sum
#>     <char> <int> <int> <char>  <char>        <num>
#>  1:   chr1     1     6      a   gene1            4
#>  2:   chr1     7    10      a   gene2            3
#>  3:   chr1    11    14      a   gene3            5
#>  4:   chr1     1     6      b   gene1            4
#>  5:   chr1     7    10      b   gene2            1
#>  6:   chr1    11    14      b   gene3            3
#>  7:   chr1     1     6      c   gene1            2
#>  8:   chr1     7    10      c   gene2            3
#>  9:   chr1    11    14      c   gene3           NA
#> 10:   chr1     1     6      d   gene1            3
#> 11:   chr1     7    10      d   gene2            3
#> 12:   chr1    11    14      d   gene3            1

# using `feature_col`
library(data.table)

# convert string vector to a data.table
regions_df <- data.table::as.data.table(regions) |>
_[, tstrsplit(regions, ":|-", fixed = FALSE, names = c("chr", "start", "end"))] |>
_[, feature := names(regions)][]
regions_df
#>       chr  start    end feature
#>    <char> <char> <char>  <char>
#> 1:   chr1      1      6   gene1
#> 2:   chr1      7     10   gene2
#> 3:   chr1     11     14   gene3

summarize_regions(
  bedfiles,
  regions_df,
  fun = "sum",
  columns = 5,
  col_names = "coverage",
  feature_col = "feature"
)
#> [15:32:59.492448] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:32:59.492461] [iscream::summarize_regions] [info] using sum
#> [15:32:59.492464] [iscream::summarize_regions] [info] with columns 5 as coverage
#>        chr start   end   file feature coverage.sum
#>     <char> <int> <int> <char>  <char>        <num>
#>  1:   chr1     1     6      a   gene1            4
#>  2:   chr1     7    10      a   gene2            3
#>  3:   chr1    11    14      a   gene3            5
#>  4:   chr1     1     6      b   gene1            4
#>  5:   chr1     7    10      b   gene2            1
#>  6:   chr1    11    14      b   gene3            3
#>  7:   chr1     1     6      c   gene1            2
#>  8:   chr1     7    10      c   gene2            3
#>  9:   chr1    11    14      c   gene3           NA
#> 10:   chr1     1     6      d   gene1            3
#> 11:   chr1     7    10      d   gene2            3
#> 12:   chr1    11    14      d   gene3            1
```
