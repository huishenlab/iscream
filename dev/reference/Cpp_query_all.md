# Query all methylation info into M and coverage matrices

Query all methylation info into M and coverage matrices

## Usage

``` r
Cpp_query_all(
  bedfiles,
  regions,
  aligner,
  valInd,
  merged,
  sparse,
  prealloc,
  nthreads
)
```

## Arguments

- bedfiles:

  A vector of BED files

- regions:

  A vector of regions

- aligner:

  The aligner used to make the WGBS BED files, only for `make_mat_bsseq`

- valInd:

  The index of the data column needed for the matrix, for `make_mat`

- merged:

  Whether the input strands have been merged/collapsed

- prealloc:

  The number of rows to initialize the matrices with

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

## Value

A list of one or two matrices, chromosome, position, and filename
vectors
