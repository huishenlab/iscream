# Apply a function over BED file records within genomic features

This function should be called from
[`summarize_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md)
since there are few sanity checks on the C++ side.

## Usage

``` r
Cpp_summarize_regions(
  bedfiles,
  regions,
  fun_vec,
  col_indices,
  col_names,
  aligner,
  mval = FALSE,
  region_rownames = FALSE,
  nthreads = 1L
)
```

## Arguments

- bedfiles:

  A vector of BED file paths

- regions:

  A vector of genomic regions

- fun_vec:

  Vector of the armadillo-supported stats functions to apply over the
  CpGs in the ' regions: `"sum"`, `"mean"`, `"median"`, `"stddev"`,
  `"variance"` "`count`", `"min"`,`"max"`, and `"range"`.

- col_indices:

  A vector of genomic regions

- col_names:

  A vector of genomic regions

- mval:

  Calculates M values when TRUE, use beta values when FALSE

- region_rownames:

  Whether to set rownames to the regions strings. Not necessary if your
  regions vector is unnamed. If its names, then the "feature" column is
  set to the names and the rownames are set to the regions string

- nthreads:

  Number of cores to use. See details.

## Value

A summary data.frame

## Details

The optimal number of threads depends on the number of bedfiles, but is
set to half the available OpenMP cores. See
[`?get_threads`](https://huishenlab.github.io/iscream/dev/reference/get_threads.md)
for more details. It can be manaully set with
[`set_threads()`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md).
