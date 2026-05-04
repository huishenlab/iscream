# Summarize methylation information over genomic regions

Run summarizing functions on the CpG/CpH loci in BED files across
genomic regions. Parallelized across files using threads from the
`"iscream.threads"` option.

## Usage

``` r
summarize_meth_regions(
  bedfiles,
  regions,
  fun = "all",
  aligner = "biscuit",
  feature_col = NULL,
  mval = TRUE,
  nthreads = NULL
)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- regions:

  A vector, data frame or GenomicRanges of genomic regions. See details.

- fun:

  Function(s) to apply over the region. See details.

- aligner:

  The aligner used to produce the BED files - one of "biscuit",
  "bismark", "bsbolt".

- feature_col:

  Column name of the input `regions` data frame or GRanges `mcols`
  containing a name for each genomic region. Set only if the using a
  data frame-like object as the input regions format, not for string
  vectors. See details.

- mval:

  Whether to calculate the M value (coverage \\\times \beta\\) or use
  the beta value when applying the function.

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/reference/set_threads.md)
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
# also see examples from ?summarize_regions

bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_meth_regions(bedfiles, regions)
#> [15:34:04.325795] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:34:04.325825] [iscream::summarize_regions] [info] using sum, mean, median, mode, antimode, stddev, pstddev, variance, min, max, absmin, absmax, range, first, last, count_distinct, count
#> [15:34:04.325830] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
#>        chr start   end   file coverage.sum M.sum coverage.mean    M.mean
#>     <char> <int> <int> <char>        <num> <num>         <num>     <num>
#>  1:   chr1     1     6      a            4     2      1.333333 0.6666667
#>  2:   chr1     7    10      a            3     1      1.500000 0.5000000
#>  3:   chr1    11    14      a            5     5      2.500000 2.5000000
#>  4:   chr1     1     6      b            4     2      2.000000 1.0000000
#>  5:   chr1     7    10      b            1     1      1.000000 1.0000000
#>  6:   chr1    11    14      b            3     1      1.500000 0.5000000
#>  7:   chr1     1     6      c            2     2      2.000000 2.0000000
#>  8:   chr1     7    10      c            3     1      1.500000 0.5000000
#>  9:   chr1    11    14      c           NA    NA            NA        NA
#> 10:   chr1     1     6      d            3     3      1.500000 1.5000000
#> 11:   chr1     7    10      d            3     1      1.500000 0.5000000
#> 12:   chr1    11    14      d            1     1      1.000000 1.0000000
#>     coverage.median M.median coverage.mode M.mode coverage.antimode M.antimode
#>               <num>    <num>         <num>  <num>             <num>      <num>
#>  1:             1.0      1.0             1      1                 2          0
#>  2:             1.5      0.5             1      0                 1          0
#>  3:             2.5      2.5             2      2                 2          2
#>  4:             2.0      1.0             2      0                 2          0
#>  5:             1.0      1.0             1      1                 1          1
#>  6:             1.5      0.5             1      0                 1          0
#>  7:             2.0      2.0             2      2                 2          2
#>  8:             1.5      0.5             1      0                 1          0
#>  9:              NA       NA            NA     NA                NA         NA
#> 10:             1.5      1.5             1      1                 1          1
#> 11:             1.5      0.5             1      0                 1          0
#> 12:             1.0      1.0             1      1                 1          1
#>     coverage.stddev  M.stddev coverage.pstddev M.pstddev coverage.variance
#>               <num>     <num>            <num>     <num>             <num>
#>  1:       0.5773503 0.5773503        0.4714045 0.4714045         0.3333333
#>  2:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#>  3:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#>  4:       0.0000000 1.4142136        0.0000000 1.0000000         0.0000000
#>  5:       0.0000000 0.0000000        0.0000000 0.0000000         0.0000000
#>  6:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#>  7:       0.0000000 0.0000000        0.0000000 0.0000000         0.0000000
#>  8:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#>  9:              NA        NA               NA        NA                NA
#> 10:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#> 11:       0.7071068 0.7071068        0.5000000 0.5000000         0.5000000
#> 12:       0.0000000 0.0000000        0.0000000 0.0000000         0.0000000
#>     M.variance coverage.min M.min coverage.max M.max coverage.absmin M.absmin
#>          <num>        <num> <num>        <num> <num>           <num>    <num>
#>  1:  0.3333333            1     0            2     1               1        0
#>  2:  0.5000000            1     0            2     1               1        0
#>  3:  0.5000000            2     2            3     3               2        2
#>  4:  2.0000000            2     0            2     2               2        0
#>  5:  0.0000000            1     1            1     1               1        1
#>  6:  0.5000000            1     0            2     1               1        0
#>  7:  0.0000000            2     2            2     2               2        2
#>  8:  0.5000000            1     0            2     1               1        0
#>  9:         NA           NA    NA           NA    NA              NA       NA
#> 10:  0.5000000            1     1            2     2               1        1
#> 11:  0.5000000            1     0            2     1               1        0
#> 12:  0.0000000            1     1            1     1               1        1
#>     coverage.absmax M.absmax coverage.range M.range coverage.first M.first
#>               <num>    <num>          <num>   <num>          <num>   <num>
#>  1:               2        1              1       1              1       1
#>  2:               2        1              1       1              1       0
#>  3:               3        3              1       1              2       2
#>  4:               2        2              0       2              2       0
#>  5:               1        1              0       0              1       1
#>  6:               2        1              1       1              2       0
#>  7:               2        2              0       0              2       2
#>  8:               2        1              1       1              2       0
#>  9:              NA       NA             NA      NA             NA      NA
#> 10:               2        2              1       1              1       1
#> 11:               2        1              1       1              1       0
#> 12:               1        1              0       0              1       1
#>     coverage.last M.last coverage.count_distinct M.count_distinct count
#>             <num>  <num>                   <num>            <num> <num>
#>  1:             2      0                       2                2     3
#>  2:             2      1                       2                2     2
#>  3:             3      3                       2                2     2
#>  4:             2      2                       1                2     2
#>  5:             1      1                       1                1     1
#>  6:             1      1                       2                2     2
#>  7:             2      2                       1                1     1
#>  8:             1      1                       2                2     2
#>  9:            NA     NA                      NA               NA    NA
#> 10:             2      2                       2                2     2
#> 11:             2      1                       2                2     2
#> 12:             1      1                       1                1     1

# add names to the regions to populate the 'feature' column
names(regions) <- c("gene1", "gene2", "gene3")
summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#> [15:34:04.385389] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:34:04.385407] [iscream::summarize_regions] [info] using mean, stddev
#> [15:34:04.385411] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, beta
#>        chr start   end   file feature coverage.mean beta.mean coverage.stddev
#>     <char> <int> <int> <char>  <char>         <num>     <num>           <num>
#>  1:   chr1     1     6      a   gene1      1.333333 0.6666667       0.5773503
#>  2:   chr1     7    10      a   gene2      1.500000 0.2500000       0.7071068
#>  3:   chr1    11    14      a   gene3      2.500000 1.0000000       0.7071068
#>  4:   chr1     1     6      b   gene1      2.000000 0.5000000       0.0000000
#>  5:   chr1     7    10      b   gene2      1.000000 1.0000000       0.0000000
#>  6:   chr1    11    14      b   gene3      1.500000 0.5000000       0.7071068
#>  7:   chr1     1     6      c   gene1      2.000000 1.0000000       0.0000000
#>  8:   chr1     7    10      c   gene2      1.500000 0.5000000       0.7071068
#>  9:   chr1    11    14      c   gene3            NA        NA              NA
#> 10:   chr1     1     6      d   gene1      1.500000 1.0000000       0.7071068
#> 11:   chr1     7    10      d   gene2      1.500000 0.2500000       0.7071068
#> 12:   chr1    11    14      d   gene3      1.000000 1.0000000       0.0000000
#>     beta.stddev
#>           <num>
#>  1:   0.5773503
#>  2:   0.3535534
#>  3:   0.0000000
#>  4:   0.7071068
#>  5:   0.0000000
#>  6:   0.7071068
#>  7:   0.0000000
#>  8:   0.7071068
#>  9:          NA
#> 10:   0.0000000
#> 11:   0.3535534
#> 12:   0.0000000
summarize_meth_regions(bedfiles, regions, fun = "sum")
#> [15:34:04.396370] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [15:34:04.396385] [iscream::summarize_regions] [info] using sum
#> [15:34:04.396388] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
#>        chr start   end   file feature coverage.sum M.sum
#>     <char> <int> <int> <char>  <char>        <num> <num>
#>  1:   chr1     1     6      a   gene1            4     2
#>  2:   chr1     7    10      a   gene2            3     1
#>  3:   chr1    11    14      a   gene3            5     5
#>  4:   chr1     1     6      b   gene1            4     2
#>  5:   chr1     7    10      b   gene2            1     1
#>  6:   chr1    11    14      b   gene3            3     1
#>  7:   chr1     1     6      c   gene1            2     2
#>  8:   chr1     7    10      c   gene2            3     1
#>  9:   chr1    11    14      c   gene3           NA    NA
#> 10:   chr1     1     6      d   gene1            3     3
#> 11:   chr1     7    10      d   gene2            3     1
#> 12:   chr1    11    14      d   gene3            1     1
```
