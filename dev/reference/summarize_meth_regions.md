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
  set_region_rownames = FALSE,
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

  Column name of the input `regions` data frame containing a name for
  each genomic region. Set only if the using a data frame as the input
  regions format. See details.

- mval:

  Whether to calculate the M value (coverage \\\times \beta\\) or use
  the beta value when applying the function.

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

- Mod: `"mode"`

- Anti-mode: `"antimode"`

- Standard deviation: `"stddev"`

- Variance: `"variance"`

- Minimum: `"min"`

- Maximum: `"max"`

- Range: `"range"`

- First element: `"first"`

- Last element: `"last"`

- No. of records in the region: `"count"`

- No. of records in the region with unique data values: `"count_unique"`

Most summarizing computations are backed by the Armadillo library. See
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
# also see examples from ?summarize_regions

bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_meth_regions(bedfiles, regions)
#> [16:44:29.858147] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [16:44:29.858169] [iscream::summarize_regions] [info] using sum, mean, median, mode, antimode, stddev, variance, min, max, range, first, last, count_unique, count
#> [16:44:29.858173] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
#>       feature file coverage.sum M.sum coverage.mean    M.mean coverage.median
#> 1    chr1:1-6    a            4     2      1.333333 0.6666667             1.0
#> 2   chr1:7-10    a            3     1      1.500000 0.5000000             1.5
#> 3  chr1:11-14    a            5     5      2.500000 2.5000000             2.5
#> 4    chr1:1-6    b            4     2      2.000000 1.0000000             2.0
#> 5   chr1:7-10    b            1     1      1.000000 1.0000000             1.0
#> 6  chr1:11-14    b            3     1      1.500000 0.5000000             1.5
#> 7    chr1:1-6    c            2     2      2.000000 2.0000000             2.0
#> 8   chr1:7-10    c            3     1      1.500000 0.5000000             1.5
#> 9  chr1:11-14    c           NA    NA            NA        NA              NA
#> 10   chr1:1-6    d            3     3      1.500000 1.5000000             1.5
#> 11  chr1:7-10    d            3     1      1.500000 0.5000000             1.5
#> 12 chr1:11-14    d            1     1      1.000000 1.0000000             1.0
#>    M.median coverage.mode M.mode coverage.antimode M.antimode coverage.stddev
#> 1       1.0             1      1                 2          0       0.5773503
#> 2       0.5             2      1                 2          1       0.7071068
#> 3       2.5             3      3                 3          3       0.7071068
#> 4       1.0             2      2                 2          2       0.0000000
#> 5       1.0             1      1                 1          1       0.0000000
#> 6       0.5             1      1                 1          1       0.7071068
#> 7       2.0             2      2                 2          2       0.0000000
#> 8       0.5             1      1                 1          1       0.7071068
#> 9        NA            NA     NA                NA         NA              NA
#> 10      1.5             2      2                 2          2       0.7071068
#> 11      0.5             2      1                 2          1       0.7071068
#> 12      1.0             1      1                 1          1       0.0000000
#>     M.stddev coverage.variance M.variance coverage.min M.min coverage.max M.max
#> 1  0.5773503         0.3333333  0.3333333            1     0            2     1
#> 2  0.7071068         0.5000000  0.5000000            1     0            2     1
#> 3  0.7071068         0.5000000  0.5000000            2     2            3     3
#> 4  1.4142136         0.0000000  2.0000000            2     0            2     2
#> 5  0.0000000         0.0000000  0.0000000            1     1            1     1
#> 6  0.7071068         0.5000000  0.5000000            1     0            2     1
#> 7  0.0000000         0.0000000  0.0000000            2     2            2     2
#> 8  0.7071068         0.5000000  0.5000000            1     0            2     1
#> 9         NA                NA         NA           NA    NA           NA    NA
#> 10 0.7071068         0.5000000  0.5000000            1     1            2     2
#> 11 0.7071068         0.5000000  0.5000000            1     0            2     1
#> 12 0.0000000         0.0000000  0.0000000            1     1            1     1
#>    coverage.range M.range coverage.first M.first coverage.last M.last
#> 1               1       1              1       1             2      0
#> 2               1       1              1       0             2      1
#> 3               1       1              2       2             3      3
#> 4               0       2              2       0             2      2
#> 5               0       0              1       1             1      1
#> 6               1       1              2       0             1      1
#> 7               0       0              2       2             2      2
#> 8               1       1              2       0             1      1
#> 9              NA      NA             NA      NA            NA     NA
#> 10              1       1              1       1             2      2
#> 11              1       1              1       0             2      1
#> 12              0       0              1       1             1      1
#>    coverage.count_unique M.count_unique cpg_count
#> 1                      2              2         3
#> 2                      2              2         2
#> 3                      2              2         2
#> 4                      1              2         2
#> 5                      1              1         1
#> 6                      2              2         2
#> 7                      1              1         1
#> 8                      2              2         2
#> 9                     NA             NA        NA
#> 10                     2              2         2
#> 11                     2              2         2
#> 12                     1              1         1
names(regions) <- c("A", "B", "C")
summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#> [16:44:29.890085] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [16:44:29.890097] [iscream::summarize_regions] [info] using mean, stddev
#> [16:44:29.890101] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, beta
#>    feature file coverage.mean beta.mean coverage.stddev beta.stddev
#> 1        A    a      1.333333 0.6666667       0.5773503   0.5773503
#> 2        B    a      1.500000 0.2500000       0.7071068   0.3535534
#> 3        C    a      2.500000 1.0000000       0.7071068   0.0000000
#> 4        A    b      2.000000 0.5000000       0.0000000   0.7071068
#> 5        B    b      1.000000 1.0000000       0.0000000   0.0000000
#> 6        C    b      1.500000 0.5000000       0.7071068   0.7071068
#> 7        A    c      2.000000 1.0000000       0.0000000   0.0000000
#> 8        B    c      1.500000 0.5000000       0.7071068   0.7071068
#> 9        C    c            NA        NA              NA          NA
#> 10       A    d      1.500000 1.0000000       0.7071068   0.0000000
#> 11       B    d      1.500000 0.2500000       0.7071068   0.3535534
#> 12       C    d      1.000000 1.0000000       0.0000000   0.0000000
summarize_meth_regions(bedfiles, regions, fun = "sum")
#> [16:44:29.895245] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [16:44:29.895255] [iscream::summarize_regions] [info] using sum
#> [16:44:29.895258] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
#>    feature file coverage.sum M.sum
#> 1        A    a            4     2
#> 2        B    a            3     1
#> 3        C    a            5     5
#> 4        A    b            4     2
#> 5        B    b            1     1
#> 6        C    b            3     1
#> 7        A    c            2     2
#> 8        B    c            3     1
#> 9        C    c           NA    NA
#> 10       A    d            3     3
#> 11       B    d            3     1
#> 12       C    d            1     1
```
