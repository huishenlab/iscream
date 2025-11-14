# Package index

## Read data from BED files

### General BED files

- [`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
  [`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
  [`tabix_raw()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
  : Query records from tabixed BED files
- [`query_chroms()`](https://huishenlab.github.io/iscream/dev/reference/query_chroms.md)
  : Query the chromosomes or seqnames from a vector of BED files
- [`summarize_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md)
  : Summarize information over genomic regions from any BED file
- [`make_mat()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
  [`make_mat_se()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
  [`make_mat_gr()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
  : Make a matrix from a numeric column of BED files
- [`make_mat_bsseq()`](https://huishenlab.github.io/iscream/dev/reference/make_mat_bsseq.md)
  : Make M/beta and coverage matrices from WGBS BED files

### WGBS BED files

- [`summarize_meth_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_meth_regions.md)
  : Summarize methylation information over genomic regions
- [`make_mat_bsseq()`](https://huishenlab.github.io/iscream/dev/reference/make_mat_bsseq.md)
  : Make M/beta and coverage matrices from WGBS BED files

## Utilities

- [`get_granges_string()`](https://huishenlab.github.io/iscream/dev/reference/get_granges_string.md)
  : GRanges to region strings
- [`get_df_string()`](https://huishenlab.github.io/iscream/dev/reference/get_df_string.md)
  : DataFrame to region strings

### Multithreading

- [`get_threads()`](https://huishenlab.github.io/iscream/dev/reference/get_threads.md)
  : Get the number of available threads
- [`set_threads()`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  : Set the number of available threads

### Logging

- [`set_log_level()`](https://huishenlab.github.io/iscream/dev/reference/set_log_level.md)
  [`get_log_level()`](https://huishenlab.github.io/iscream/dev/reference/set_log_level.md)
  : Set and get logging level

### htslib

- [`htslib_version()`](https://huishenlab.github.io/iscream/dev/reference/htslib_version.md)
  : Get htslib version and available features
