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
