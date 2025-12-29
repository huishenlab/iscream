# iscream compatible data structures

The examples here show how iscream output can be converted into other
data structures for further analysis.

``` r
library(iscream)
```

    ## iscream using 1 thread by default but parallelly::availableCores() detects 2 possibly available threads. See `?set_threads` for information on multithreading before trying to use more.

``` r
data_dir <- system.file("extdata", package = "iscream")
bedfiles <- list.files(data_dir, pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
```

## *[GenomicRanges](https://bioconductor.org/packages/3.22/GenomicRanges)*

``` r
if (!require("GenomicRanges", quietly = TRUE)) {
  stop("The 'GenomicRanges' package must be installed for this functionality")
}
```

`GRanges` objects can be used as the input regions to all of iscream’s
functions and can be returned by
[`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
and
[`make_mat_gr()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md).

### From `tabix` queries

[`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
returns a `GenomicRanges` object. The `regions` parameter can be a
string vector, a data frame or a `GRanges` object. If given a `GRanges`
object with metadata, those columns will be preserved in the output.

``` r
gr <- GRanges(regions)
values(gr) <- DataFrame(
    gene = c("gene1", "gene2", "gene3"),
    some_metadata = c("s1", "s2", "s3")
)
gr
```

    ## GRanges object with 3 ranges and 2 metadata columns:
    ##     seqnames    ranges strand |        gene some_metadata
    ##        <Rle> <IRanges>  <Rle> | <character>   <character>
    ##   A     chr1       1-6      * |       gene1            s1
    ##   B     chr1      7-10      * |       gene2            s2
    ##   C     chr1     11-14      * |       gene3            s3
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
tabix_gr(bedfiles[1], gr)
```

    ## GRanges object with 7 ranges and 4 metadata columns:
    ##       seqnames    ranges strand |        V1        V2        gene some_metadata
    ##          <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>   <character>
    ##   [1]     chr1       1-2      * |       1.0         1       gene1            s1
    ##   [2]     chr1       3-4      * |       1.0         1       gene1            s1
    ##   [3]     chr1       5-6      * |       0.0         2       gene1            s1
    ##   [4]     chr1       7-8      * |       0.0         1       gene2            s2
    ##   [5]     chr1      9-10      * |       0.5         2       gene2            s2
    ##   [6]     chr1     11-12      * |       1.0         2       gene3            s3
    ##   [7]     chr1     13-14      * |       1.0         3       gene3            s3
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

The `data.table` output of
[`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
can also be piped into `GRanges`, but does not preserve input metadata.

``` r
tabix(bedfiles[1], gr) |>
  makeGRangesFromDataFrame(
    starts.in.df.are.0based = TRUE,
    keep.extra.columns = TRUE
  )
```

    ## GRanges object with 7 ranges and 2 metadata columns:
    ##       seqnames    ranges strand |        V1        V2
    ##          <Rle> <IRanges>  <Rle> | <numeric> <integer>
    ##   [1]     chr1       1-2      * |       1.0         1
    ##   [2]     chr1       3-4      * |       1.0         1
    ##   [3]     chr1       5-6      * |       0.0         2
    ##   [4]     chr1       7-8      * |       0.0         1
    ##   [5]     chr1      9-10      * |       0.5         2
    ##   [6]     chr1     11-12      * |       1.0         2
    ##   [7]     chr1     13-14      * |       1.0         3
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

If the input BED file is not zero-based (e.g. Bismark coverage files),
set `zero_based = FALSE` in the
[`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
call to get the correct conversion from data frame to GenomicRanges.

### From `summarize_regions()`

[`summarize_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md)
returns a data frame with the input `regions` positions and summaries of
the BED-file data columns.

``` r
(summary <- summarize_regions(
  bedfiles,
  regions,
  column = 4,
  col_names = ("data_col"),
  fun = c("sum", "mean"))
)
```

    ## [20:57:08.315517] [iscream::summarize_regions] [info] Summarizing 3 regions from 4 bedfiles
    ## [20:57:08.315549] [iscream::summarize_regions] [info] using sum, mean
    ## [20:57:08.315555] [iscream::summarize_regions] [info] with columns 4 as data_col

    ##        chr start   end   file feature data_col.sum data_col.mean
    ##     <char> <int> <int> <char>  <char>        <num>         <num>
    ##  1:   chr1     1     6      a       A          2.0     0.6666667
    ##  2:   chr1     7    10      a       B          0.5     0.2500000
    ##  3:   chr1    11    14      a       C          2.0     1.0000000
    ##  4:   chr1     1     6      b       A          1.0     0.5000000
    ##  5:   chr1     7    10      b       B          1.0     1.0000000
    ##  6:   chr1    11    14      b       C          1.0     0.5000000
    ##  7:   chr1     1     6      c       A          1.0     1.0000000
    ##  8:   chr1     7    10      c       B          1.0     0.5000000
    ##  9:   chr1    11    14      c       C           NA            NA
    ## 10:   chr1     1     6      d       A          2.0     1.0000000
    ## 11:   chr1     7    10      d       B          0.5     0.2500000
    ## 12:   chr1    11    14      d       C          1.0     1.0000000

``` r
GRanges(summary)
```

    ## GRanges object with 12 ranges and 4 metadata columns:
    ##        seqnames    ranges strand |        file     feature data_col.sum
    ##           <Rle> <IRanges>  <Rle> | <character> <character>    <numeric>
    ##    [1]     chr1       1-6      * |           a           A          2.0
    ##    [2]     chr1      7-10      * |           a           B          0.5
    ##    [3]     chr1     11-14      * |           a           C          2.0
    ##    [4]     chr1       1-6      * |           b           A          1.0
    ##    [5]     chr1      7-10      * |           b           B          1.0
    ##    ...      ...       ...    ... .         ...         ...          ...
    ##    [8]     chr1      7-10      * |           c           B          1.0
    ##    [9]     chr1     11-14      * |           c           C           NA
    ##   [10]     chr1       1-6      * |           d           A          2.0
    ##   [11]     chr1      7-10      * |           d           B          0.5
    ##   [12]     chr1     11-14      * |           d           C          1.0
    ##        data_col.mean
    ##            <numeric>
    ##    [1]      0.666667
    ##    [2]      0.250000
    ##    [3]      1.000000
    ##    [4]      0.500000
    ##    [5]      1.000000
    ##    ...           ...
    ##    [8]          0.50
    ##    [9]            NA
    ##   [10]          1.00
    ##   [11]          0.25
    ##   [12]          1.00
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

### From `make_mat`

[`make_mat_gr()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
returns a `GRanges` object for dense matrices.

``` r
make_mat_gr(bedfiles, regions, column = 4, mat_name = "beta")
```

    ## [20:57:08.403125] [iscream::query_all] [info] Querying 3 regions from 4 bedfiles
    ## 
    ## [20:57:08.403451] [iscream::query_all] [info] Creating metadata vectors
    ## [20:57:08.403489] [iscream::query_all] [info] 7 loci found - 9993 extra rows allocated with 0 resizes
    ## [20:57:08.403493] [iscream::query_all] [info] Creating dense matrix

    ## GRanges object with 7 ranges and 4 metadata columns:
    ##       seqnames    ranges strand |         a         b         c         d
    ##          <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>
    ##   [1]     chr1         0      * |       1.0         0         0       1.0
    ##   [2]     chr1         2      * |       1.0         0         1       1.0
    ##   [3]     chr1         4      * |       0.0         1         0       0.0
    ##   [4]     chr1         6      * |       0.0         1         0       0.0
    ##   [5]     chr1         8      * |       0.5         0         1       0.5
    ##   [6]     chr1        10      * |       1.0         0         0       0.0
    ##   [7]     chr1        12      * |       1.0         1         0       1.0
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## *[SummarizedExperiment](https://bioconductor.org/packages/3.22/SummarizedExperiment)*

[`make_mat_se()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
returns a `RangedSummarizedExperiment` for both sparse and dense
matrices.

``` r
if (!require("SummarizedExperiment", quietly = TRUE)) {
  stop("The 'SummarizedExperiment' package must be installed for this functionality")
}
make_mat_se(bedfiles, regions, column = 4, mat_name = "beta", sparse = TRUE)
```

    ## [20:57:10.738360] [iscream::query_all] [info] Querying 3 regions from 4 bedfiles
    ## 
    ## [20:57:10.738735] [iscream::query_all] [info] Creating metadata vectors
    ## [20:57:10.738757] [iscream::query_all] [info] 7 loci found - 9993 extra rows allocated with 0 resizes
    ## [20:57:10.738763] [iscream::query_all] [info] Creating sparse matrix

    ## class: RangedSummarizedExperiment 
    ## dim: 7 4 
    ## metadata(0):
    ## assays(1): beta
    ## rownames: NULL
    ## rowData names(0):
    ## colnames(4): a b c d
    ## colData names(0):

### Making `BSseq` objects

A *[bsseq](https://bioconductor.org/packages/3.22/bsseq)* object is a
type of SummarizedExperiment, but it cannot handle sparse matrices:

``` r
if (!require("bsseq", quietly = TRUE)) {
  stop("The 'bsseq' package must be installed for this functionality")
}
mats <- make_mat_bsseq(bedfiles, regions, sparse = FALSE)
```

    ## [20:57:14.588136] [iscream::query_all] [info] Querying 3 regions from 4 bedfiles
    ## 
    ## [20:57:14.588514] [iscream::query_all] [info] Creating metadata vectors
    ## [20:57:14.588539] [iscream::query_all] [info] 7 loci found - 9993 extra rows allocated with 0 resizes
    ## [20:57:14.588543] [iscream::query_all] [info] Creating dense matrix

``` r
do.call(BSseq, mats)
```

    ## An object of type 'BSseq' with
    ##   7 methylation loci
    ##   4 samples
    ## has not been smoothed
    ## All assays are in-memory

## Session info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] bsseq_1.46.0                SummarizedExperiment_1.40.0
    ##  [3] Biobase_2.70.0              MatrixGenerics_1.22.0      
    ##  [5] matrixStats_1.5.0           GenomicRanges_1.62.1       
    ##  [7] Seqinfo_1.0.0               IRanges_2.44.0             
    ##  [9] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [11] generics_0.1.4              iscream_1.1.8              
    ## [13] BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] farver_2.1.2              R.utils_2.13.0           
    ##  [3] Biostrings_2.78.0         bitops_1.0-9             
    ##  [5] fastmap_1.2.0             RCurl_1.98-1.17          
    ##  [7] GenomicAlignments_1.46.0  stringfish_0.17.0        
    ##  [9] XML_3.99-0.20             digest_0.6.39            
    ## [11] lifecycle_1.0.4           statmod_1.5.1            
    ## [13] compiler_4.5.2            rlang_1.1.6              
    ## [15] sass_0.4.10               tools_4.5.2              
    ## [17] yaml_2.3.12               data.table_1.18.0        
    ## [19] rtracklayer_1.70.1        knitr_1.51               
    ## [21] S4Arrays_1.10.1           curl_7.0.0               
    ## [23] DelayedArray_0.36.0       RColorBrewer_1.1-3       
    ## [25] abind_1.4-8               BiocParallel_1.44.0      
    ## [27] HDF5Array_1.38.0          R.oo_1.27.1              
    ## [29] desc_1.4.3                grid_4.5.2               
    ## [31] beachmat_2.26.0           Rhdf5lib_1.32.0          
    ## [33] scales_1.4.0              gtools_3.9.5             
    ## [35] cli_3.6.5                 rmarkdown_2.30           
    ## [37] crayon_1.5.3              ragg_1.5.0               
    ## [39] RcppParallel_5.1.11-1     httr_1.4.7               
    ## [41] rjson_0.2.23              DelayedMatrixStats_1.32.0
    ## [43] pbapply_1.7-4             cachem_1.1.0             
    ## [45] rhdf5_2.54.1              parallel_4.5.2           
    ## [47] BiocManager_1.30.27       XVector_0.50.0           
    ## [49] restfulr_0.0.16           Matrix_1.7-4             
    ## [51] jsonlite_2.0.0            bookdown_0.46            
    ## [53] systemfonts_1.3.1         h5mread_1.2.1            
    ## [55] locfit_1.5-9.12           limma_3.66.0             
    ## [57] jquerylib_0.1.4           glue_1.8.0               
    ## [59] parallelly_1.46.0         pkgdown_2.2.0            
    ## [61] codetools_0.2-20          BiocIO_1.20.0            
    ## [63] htmltools_0.5.9           rhdf5filters_1.22.0      
    ## [65] BSgenome_1.78.0           R6_2.6.1                 
    ## [67] textshaping_1.0.4         sparseMatrixStats_1.22.0 
    ## [69] evaluate_1.0.5            lattice_0.22-7           
    ## [71] R.methodsS3_1.8.2         Rsamtools_2.26.0         
    ## [73] cigarillo_1.0.0           bslib_0.9.0              
    ## [75] Rcpp_1.1.0                SparseArray_1.10.8       
    ## [77] permute_0.9-8             xfun_0.55                
    ## [79] fs_1.6.6
