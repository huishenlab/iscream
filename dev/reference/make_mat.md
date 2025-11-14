# Make a matrix from a numeric column of BED files

Queries the provided regions and produces a matrix along with genomic
positions as a named list (`make_mat()`), a `RangedSummarizedExperiment`
(`make_mat_se()`), `GRanges` (`make_mat_gr()`). Parallelized across
files using threads from the `"iscream.threads"` option.

## Usage

``` r
make_mat(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
)

make_mat_se(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
)

make_mat_gr(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  prealloc = 10000,
  nthreads = NULL
)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- regions:

  A vector, data frame or GenomicRanges of genomic regions. See details.

- column:

  The index of the data column needed for the matrix

- mat_name:

  What to name the matrix in the returned object

- sparse:

  Whether to return a sparse matrix

- prealloc:

  The number of rows to initialize the matrices with. If the number of
  loci are approximately known, this can reduce runtime as fewer resizes
  need to be made.

- nthreads:

  Set the number of threads to use. Overrides the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

## Value

- `make_mat()`: A named list of

  - the matrix with the value of interest

  - a character vector of chromosomes and numeric vector of base
    positions

  - a character vector of the input sample BED file names

- `make_mat_gr()`: if `GenomicRanges` is available, a `GRanges`

- `make_mat_se()`: if `SummarizedExperiment` is available, a
  `RangedSummarizedExperiment`

## Details

The input regions may be string vector in the form "chr:start-end" or a
GRanges object. If a data frame is provided, they must have "chr",
"start", and "end" columns.

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
# make matrix of beta values
make_mat(bedfiles, regions, column = 4)
#> [20:39:51.778647] [iscream::query_all] [info] Querying 3 regions from 4 bedfiles
#> 
#> [20:39:51.779072] [iscream::query_all] [info] Creating metadata vectors
#> [20:39:51.779118] [iscream::query_all] [info] 7 loci found - 9993 extra rows allocated with 0 resizes
#> [20:39:51.779122] [iscream::query_all] [info] Creating dense matrix
#> $value
#>        a b c   d
#> [1,] 1.0 0 0 1.0
#> [2,] 1.0 0 1 1.0
#> [3,] 0.0 1 0 0.0
#> [4,] 0.0 1 0 0.0
#> [5,] 0.5 0 1 0.5
#> [6,] 1.0 0 0 0.0
#> [7,] 1.0 1 0 1.0
#> 
#> $pos
#> [1]  0  2  4  6  8 10 12
#> 
#> $chr
#> [1] "chr1" "chr1" "chr1" "chr1" "chr1" "chr1" "chr1"
#> 
#> $sampleNames
#> [1] "a" "b" "c" "d"
#> 
```
