# iscream 1.1.6

BUG FIXES

 - Fix a check for missing input files by correcting checking the filepath
   vector length

INTERNAL

 - Use `inherits(x, "data.frame")` instead of `"data.frame" %in% class(x)"`

# iscream 1.1.5

BREAKING CHANGES to `summarize_regions()` and `summarize_meth_regions()`:

 - The summary output now has position as chromosome, start, and end columns
   instead of just a `feature` column with the position string. This allows
   direct conversion to `GRanges` objects without having to use rownames. It is
   still possible to use chromosome names as the input query regions ("chr1"),
   the start and end in those cases will be `NA`.

 - The `feature` column is now populated only if the input regions are named, if
   a vector, or the `feature_col` argument is set to a column for
   `data.frame`/`GRanges` inputs.

 - The count column from `summarize_meth_regions()` is now named `count` instead
   of `cpg_count`.

 - These functions now return a `data.table` instead of a `data.frame` to allow
   for faster in-place modifications and consistency with `tabix()` output. Use
   `data.table::setDF()` to convert to a `data.frame` in-place.

BUG FIXES

 - Rownames are no longer set for the `summarize_regions()` or
   `summarize_meth_regions()` output using the input region strings as these
   were repeating and not unique. It was only possible because they were set
   from the C++ scope instead of R.

ENHANCEMENTS

 - `summarize_regions()` and `summarize_meth_regions()` now support the following
   new functions, inspired by [`bedtools map`](https://bedtools.readthedocs.io/en/latest/content/tools/map.html):

     * Population standard deviation (`pstddev`)
     * First element (`first`)
     * Last element (`last`)
     * Mode (`mode`)
     * Anti-mode (`antimode`)
     * Absolute min (`absmin`)
     * Absolute max (`absmax`)
     * Count of unique values (`count_distinct`)

- `get_granges_string()` can now extract names from its `mcols` using a
  `feature_col` argument as in `get_df_string()` - names were
  previously pulled only from `names()`. This means GRanges inputs to
  `summarize_regions()` can have the `feature_col` in its `mcols` rather than
  just as its names.

DOCUMENTATION

- Update Zenodo URLs in vignettes to the updated record:
  <https://zenodo.org/records/18089082>. The `genes.bed` file used in
  `vignette("performance"`) now contains the gene names that the vignette
  references.

INTERNAL

- Refactored `summarize_regions` to collect similar validation and regions
  parsing functions

- `tabix()` writes the input regions of interest to disk only once instead of
  doing it for every file


# iscream 1.1.4

BUG FIX

- Reduced the minimum R version from 4.5 to 4.4

# iscream 1.1.3

DOCUMENTATION

- Document how to make use of strand information with `tabix()`
- Improve wording of `tabix()`'s' `aligner` and `col.names` documentation

# iscream 1.1.2

- Fix typo in error message from thread count checks

# iscream 1.1.1

BUG FIXES

- Fix configure to correctly check htslib version when rhtslib is not used and
  not issue spurious warnings

# iscream 1.0.0

- First Bioconductor release

# iscream 0.99.9

DOCUMENTATION

- Vignettes suppress long startup messages

- Updated run-times and plots in tabix vs Rsamtools vignette

# iscream 0.99.8

BUG FIXES

- For `tabix()`, converting GRanges inputs of width 1 now returns them as
  `chr:start-start` instead of `chr:start` which was invalid for the htslib API
  and caused a block of records to be returned instead of just the locus.

# iscream 0.99.7

BREAKING CHANGES

- Functions that return GenomicRanges or SummarizedExperiment
  will fail if those packages are not loaded so that users don't receive objects
  they can't interact with

DOCUMENTATION

- `bsseq` creation now runs for the `make_mat_bsseq` example

# iscream 0.99.6

BREAKING CHANGES

- `make_bsseq_mat` changed to `make_mat_bsseq` for consistency with other matrix
  creating functions

BUG FIXES

- Fixed MacOS and clang issues (version checking, using the correct linkage
  tools) in configuration script

- Explicitly coerce path to string for bedfile prefix extraction in
  `Cpp_summarize_regions`

ENHANCEMENTS

- Added windows as supported platform but without multithreading support

- `pbapply` progress bar added to `tabix()` functions

INTERNAL

- Explicitly return result from `tabix.single` helper function.

- Removed unused C++ variables

- `file.path` used instead of pasting paths

DOCUMENTATION

- `bsseq` added to Suggests to run in examples and vignettes

- Changed all `\dontrun` to `\donttest`

- All vignettes set to BiocStyle

- BiocFileCache now used for all vignette file downloads

- Reading the genes file in the "performance" vignette fixed to read only three
  columns

- Added more information on bitpacking limits for `make_mat_bsseq`

- Allowed creation of BSseq objects in data structures vignette

- Package description updated to add information about methylation support

- Replaced all references to CpG with more general methylation

- Added CpH as supported in documentation

# iscream 0.99.5

INTERNAL

- tabix executable-dependent tests are skipped if tabix is not installed

# iscream 0.99.4

- tabix executable added to `SystemRequirements`

# iscream 0.99.1

- Updated LICENSE files for R CMD check

# iscream 0.99.0

- Pre-release for manuscript and Bioconductor review
