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
  set_region_rownames = FALSE,
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

  Column name of the input `regions` data frame containing a name for
  each genomic region. Set only if the using a data frame as the input
  regions format. See details.

- set_region_rownames:

  Use the region strings as the returned data frame's rownames. Can be
  useful if you have a named regions and want both the regions strings
  rownames and the feature names. See details.

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

## Value

A data.frame

## Supported functions

- Sum: `"sum"`

- Mean: `"mean"`

- Median: `"median"`

- Standard deviation: `"stddev"`

- Variance: `"variance"`

- Minimum: `"min"`

- Maximum: `"max"`

- Range: `"range"`

- No. of records in the region: `"count"`

The summarizing computations are backed by the Armadillo library. See
<https://arma.sourceforge.net/docs.html#stats_fns> for futher details on
the supported functions

## Using feature identifiers

`regions` may be string vector in the form "chr:start-end", a GRanges
object or a data frame with "chr", "start", and "end" columns. The
`feature` column of the output will contain a "chr:start-end" identifier
for each summarized region. To use other identifiers, like a gene name
for a region instead of the coordinates, set the names of the vector or
GRanges to those identifiers. These names will be used instead of the
genomic region string to describe each feature in the output dataframe.
If `regions` is a data frame make an additional column with the
identifiers and pass that column name to `feature_col`. See examples.

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
#> [20:40:00.285775] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.285790] [iscream::summarize_regions] [info] using sum, mean, median, stddev, variance, min, max, range, count
#> [20:40:00.285794] [iscream::summarize_regions] [info] with columns 4, 5 as beta, cov
#>       feature file beta.sum cov.sum beta.mean cov.mean beta.median cov.median
#> 1    chr1:1-6    a      2.0       4 0.6666667 1.333333        1.00        1.0
#> 2   chr1:7-10    a      0.5       3 0.2500000 1.500000        0.25        1.5
#> 3  chr1:11-14    a      2.0       5 1.0000000 2.500000        1.00        2.5
#> 4    chr1:1-6    b      1.0       4 0.5000000 2.000000        0.50        2.0
#> 5   chr1:7-10    b      1.0       1 1.0000000 1.000000        1.00        1.0
#> 6  chr1:11-14    b      1.0       3 0.5000000 1.500000        0.50        1.5
#> 7    chr1:1-6    c      1.0       2 1.0000000 2.000000        1.00        2.0
#> 8   chr1:7-10    c      1.0       3 0.5000000 1.500000        0.50        1.5
#> 9  chr1:11-14    c       NA      NA        NA       NA          NA         NA
#> 10   chr1:1-6    d      2.0       3 1.0000000 1.500000        1.00        1.5
#> 11  chr1:7-10    d      0.5       3 0.2500000 1.500000        0.25        1.5
#> 12 chr1:11-14    d      1.0       1 1.0000000 1.000000        1.00        1.0
#>    beta.stddev cov.stddev beta.variance cov.variance beta.min cov.min beta.max
#> 1    0.5773503  0.5773503     0.3333333    0.3333333        0       1      1.0
#> 2    0.3535534  0.7071068     0.1250000    0.5000000        0       1      0.5
#> 3    0.0000000  0.7071068     0.0000000    0.5000000        1       2      1.0
#> 4    0.7071068  0.0000000     0.5000000    0.0000000        0       2      1.0
#> 5    0.0000000  0.0000000     0.0000000    0.0000000        1       1      1.0
#> 6    0.7071068  0.7071068     0.5000000    0.5000000        0       1      1.0
#> 7    0.0000000  0.0000000     0.0000000    0.0000000        1       2      1.0
#> 8    0.7071068  0.7071068     0.5000000    0.5000000        0       1      1.0
#> 9           NA         NA            NA           NA       NA      NA       NA
#> 10   0.0000000  0.7071068     0.0000000    0.5000000        1       1      1.0
#> 11   0.3535534  0.7071068     0.1250000    0.5000000        0       1      0.5
#> 12   0.0000000  0.0000000     0.0000000    0.0000000        1       1      1.0
#>    cov.max beta.range cov.range count
#> 1        2        1.0         1     3
#> 2        2        0.5         1     2
#> 3        3        0.0         1     2
#> 4        2        1.0         0     2
#> 5        1        0.0         0     1
#> 6        2        1.0         1     2
#> 7        2        0.0         0     1
#> 8        2        1.0         1     2
#> 9       NA         NA        NA    NA
#> 10       2        0.0         1     2
#> 11       2        0.5         1     2
#> 12       1        0.0         0     1

# select functions
summarize_regions(
  bedfiles,
  regions,
  fun = c("mean", "stddev"),
  columns = c(4, 5),
  col_names = c("beta", "cov")
)
#> [20:40:00.305708] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.305721] [iscream::summarize_regions] [info] using mean, stddev
#> [20:40:00.305724] [iscream::summarize_regions] [info] with columns 4, 5 as beta, cov
#>       feature file beta.mean cov.mean beta.stddev cov.stddev
#> 1    chr1:1-6    a 0.6666667 1.333333   0.5773503  0.5773503
#> 2   chr1:7-10    a 0.2500000 1.500000   0.3535534  0.7071068
#> 3  chr1:11-14    a 1.0000000 2.500000   0.0000000  0.7071068
#> 4    chr1:1-6    b 0.5000000 2.000000   0.7071068  0.0000000
#> 5   chr1:7-10    b 1.0000000 1.000000   0.0000000  0.0000000
#> 6  chr1:11-14    b 0.5000000 1.500000   0.7071068  0.7071068
#> 7    chr1:1-6    c 1.0000000 2.000000   0.0000000  0.0000000
#> 8   chr1:7-10    c 0.5000000 1.500000   0.7071068  0.7071068
#> 9  chr1:11-14    c        NA       NA          NA         NA
#> 10   chr1:1-6    d 1.0000000 1.500000   0.0000000  0.7071068
#> 11  chr1:7-10    d 0.2500000 1.500000   0.3535534  0.7071068
#> 12 chr1:11-14    d 1.0000000 1.000000   0.0000000  0.0000000

# add names to the regions
names(regions) <- c("A", "B", "C")
summarize_regions(
  bedfiles,
  regions,
  fun = "sum",
  columns = 5,
  col_names = "coverage"
)
#> [20:40:00.312894] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.312908] [iscream::summarize_regions] [info] using sum
#> [20:40:00.312911] [iscream::summarize_regions] [info] with columns 5 as coverage
#>    feature file coverage.sum
#> 1        A    a            4
#> 2        B    a            3
#> 3        C    a            5
#> 4        A    b            4
#> 5        B    b            1
#> 6        C    b            3
#> 7        A    c            2
#> 8        B    c            3
#> 9        C    c           NA
#> 10       A    d            3
#> 11       B    d            3
#> 12       C    d            1

# using `feature_col`
library(data.table)

# convert string vector to a data.table
regions_df <- data.table::as.data.table(regions) |>
_[, tstrsplit(regions, ":|-", fixed = FALSE, names = c("chr", "start", "end"))] |>
_[, start := as.integer(start)] |>
_[, feature := LETTERS[.I]][]
regions_df
#>       chr start    end feature
#>    <char> <int> <char>  <char>
#> 1:   chr1     1      6       A
#> 2:   chr1     7     10       B
#> 3:   chr1    11     14       C

summarize_regions(
  bedfiles,
  regions_df,
  fun = "sum",
  columns = 5,
  col_names = "coverage",
  feature_col = "feature"
)
#> [20:40:00.332404] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.332417] [iscream::summarize_regions] [info] using sum
#> [20:40:00.332420] [iscream::summarize_regions] [info] with columns 5 as coverage
#>    feature file coverage.sum
#> 1        A    a            4
#> 2        B    a            3
#> 3        C    a            5
#> 4        A    b            4
#> 5        B    b            1
#> 6        C    b            3
#> 7        A    c            2
#> 8        B    c            3
#> 9        C    c           NA
#> 10       A    d            3
#> 11       B    d            3
#> 12       C    d            1
```
