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
# also see examples from ?summarize_regions

bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)

# make a vector of regions
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
summarize_meth_regions(bedfiles, regions)
#> [20:40:00.020264] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.020290] [iscream::summarize_regions] [info] using sum, mean, median, stddev, variance, min, max, range, count
#> [20:40:00.020295] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
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
#>    M.median coverage.stddev  M.stddev coverage.variance M.variance coverage.min
#> 1       1.0       0.5773503 0.5773503         0.3333333  0.3333333            1
#> 2       0.5       0.7071068 0.7071068         0.5000000  0.5000000            1
#> 3       2.5       0.7071068 0.7071068         0.5000000  0.5000000            2
#> 4       1.0       0.0000000 1.4142136         0.0000000  2.0000000            2
#> 5       1.0       0.0000000 0.0000000         0.0000000  0.0000000            1
#> 6       0.5       0.7071068 0.7071068         0.5000000  0.5000000            1
#> 7       2.0       0.0000000 0.0000000         0.0000000  0.0000000            2
#> 8       0.5       0.7071068 0.7071068         0.5000000  0.5000000            1
#> 9        NA              NA        NA                NA         NA           NA
#> 10      1.5       0.7071068 0.7071068         0.5000000  0.5000000            1
#> 11      0.5       0.7071068 0.7071068         0.5000000  0.5000000            1
#> 12      1.0       0.0000000 0.0000000         0.0000000  0.0000000            1
#>    M.min coverage.max M.max coverage.range M.range cpg_count
#> 1      0            2     1              1       1         3
#> 2      0            2     1              1       1         2
#> 3      2            3     3              1       1         2
#> 4      0            2     2              0       2         2
#> 5      1            1     1              0       0         1
#> 6      0            2     1              1       1         2
#> 7      2            2     2              0       0         1
#> 8      0            2     1              1       1         2
#> 9     NA           NA    NA             NA      NA        NA
#> 10     1            2     2              1       1         2
#> 11     0            2     1              1       1         2
#> 12     1            1     1              0       0         1
names(regions) <- c("A", "B", "C")
summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#> [20:40:00.041162] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.041174] [iscream::summarize_regions] [info] using mean, stddev
#> [20:40:00.041178] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, beta
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
#> [20:40:00.046921] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
#> [20:40:00.046934] [iscream::summarize_regions] [info] using sum
#> [20:40:00.046937] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M
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
