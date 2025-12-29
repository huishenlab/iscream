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
  regions_df,
  aligner,
  mval = FALSE,
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
