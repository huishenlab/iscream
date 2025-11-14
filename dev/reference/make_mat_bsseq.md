# Make M/beta and coverage matrices from WGBS BED files

Queries the CpG/CpH loci from provided regions and produces M/beta and
coverage matrices with their genomic positions. Parallelized across
files using threads from the `"iscream.threads"` option. The output of
`make_mat_bsseq` may be used to create a BSseq object:
`do.call(BSseq, make_mat_bsseq(...))`.

## Usage

``` r
make_mat_bsseq(
  bedfiles,
  regions,
  aligner = "biscuit",
  mval = TRUE,
  merged = TRUE,
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- regions:

  A vector, data frame or GenomicRanges of genomic regions. See details.

- aligner:

  The aligner used to produce the BED files - one of "biscuit",
  "bismark", "bsbolt".

- mval:

  Whether to return M-values or beta-values with the coverage matrix.
  Defaults to M-value. Set `mval=FALSE` to get beta value matrix.

- merged:

  Whether the input strands have been merged/collapsed

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

A named list of

- coverage and either a beta- or M-value matrix

- a character vector of chromosomes and numeric vector of corresponding
  CpG base positions

- a character vector of the input sample names

## Details

The input regions may be string vector in the form "chr:start-end" or a
GRanges object. If a data frame is provided, they must have "chr",
"start", and "end" columns.

## Bitpacking limits

`make_mat_bsseq()` makes two matrices: M-value (or beta-value) and
coverage. For speed and memory efficiency these two values are bitpacked
during matrix creation so that only one matrix needs to be populated and
resized. This matrix is unpacked into the two required matrices only
after the matrix dimensions are known after querying all input files.
The two values are packed using the INT16 type, which has an upper limit
of 32,767, into one INT32. If the coverage values exceed 32,767, the
upper limit of a 16-bit signed integer, it will be capped at the limit.
Beta values will also be capped similarly, but any such beta values
would indicate a bug in the aligner that produced the data.

## Examples

``` r
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
# examine the BED files
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
mat <- make_mat_bsseq(bedfiles, regions)
#> [20:39:52.256743] [iscream::query_all] [info] Querying 3 regions from 4 bedfiles
#> 
#> [20:39:52.257104] [iscream::query_all] [info] Creating metadata vectors
#> [20:39:52.257135] [iscream::query_all] [info] 7 loci found - 9993 extra rows allocated with 0 resizes
#> [20:39:52.257139] [iscream::query_all] [info] Creating dense matrix
# for BSseq object run
if (requireNamespace("bsseq", quietly = TRUE)) {
  do.call(bsseq::BSseq, mat)
}
#> An object of type 'BSseq' with
#>   7 methylation loci
#>   4 samples
#> has not been smoothed
#> All assays are in-memory
```
