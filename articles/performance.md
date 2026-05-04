# Improving iscream performance

Abstract

A guide to improving iscream’s runtime and memory efficiency

## Setup

### Input BED files

Running this vignette requires downloading 2GB of single-cell whole
genome bisulfite sequencing (WGBS) BED files and tabix indices from this
Zenodo record: <https://zenodo.org/records/18089082>.

``` r

library("BiocFileCache")
#> Loading required package: dbplyr
cachedir <- BiocFileCache()
snmc_zip_path <- bfcrpath(cachedir, "https://zenodo.org/records/18089082/files/sc_beds.zip")
snmc_unzip <- file.path(tempdir(), "sc_beds")
unzip(snmc_zip_path, exdir = snmc_unzip)
genes_file <- bfcrpath(cachedir, "https://zenodo.org/records/18089082/files/genes.bed")
```

Select 100 human cell WGBS data from the snmC-seq2 ([Luo et al.
2018](#ref-luo2018a)) dataset:

``` r

snmc_dir <- file.path(snmc_unzip, "sc_biscuit")
bedfiles <- list.files(
  snmc_dir,
  pattern = "*.bed.gz$",
  full.names = TRUE
)[1:100]
```

### Regions

Here we’ll be using 5000 gene body regions as the input:

``` r

library(data.table)
regions <- fread(
  genes_file,
  col.names = c("chr", "start", "end", "gene")
)[1:5000]
head(regions)
#>       chr    start      end   gene
#>    <char>    <int>    <int> <char>
#> 1:   chr1  1471764  1497848 ATAD3B
#> 2:   chr1  3069167  3438621 PRDM16
#> 3:   chr1  2403963  2413797  PEX10
#> 4:   chr1 10472287 10630758  PEX14
#> 5:   chr1  2425979  2505532  PLCH2
#> 6:   chr1  9292893  9369532  SPSB1
```

## Multithreading

iscream uses multithreading to reduce querying runtime. The number of
threads can be set before or after loading the library. iscream will
read the option while loading and so can be set in `.Rprofile`

``` r

options("iscream.threads" = 8)
library(iscream)
#> iscream using 8 threads based on 'options(iscream.threads)' but parallelly::availableCores() detects 16 possibly available threads. See `?set_threads` for information on multithreading before trying to use more.
```

Queries are parallelized across the input BED files. While this reduces
runtime, it increases memory usage as data from multiple files are
loaded into memory at the same time. For more information on
multithreading see `?set_threads()`.

### `tabix()`

[`tabix()`](https://huishenlab.github.io/iscream/reference/tabix.md)
queries regions from BED files much like the *tabix* shell command, but
supports multiple files and can query them in parallel. Although fast,
[`tabix()`](https://huishenlab.github.io/iscream/reference/tabix.md) may
seem unresponsive on large queries as
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
doesn’t show progression, but a progress bar is in development.

iscream has two
[`tabix()`](https://huishenlab.github.io/iscream/reference/tabix.md)
implementations, one using the command line `tabix` executable and the
other using the htslib API. The command line tool is faster as it can
stream to a file which can be read from R while the htslib API stores
the strings in memory and is slower. iscream will look for the tabix
executable on startup and will only fall back to the htslib API if the
executable is not found. See
[`?tabix`](https://huishenlab.github.io/iscream/reference/tabix.md)
details for more information.

``` r

qt <- system.time(
  tbx_query <- tabix(bedfiles, regions, col.names = c("beta", "coverage"))
)
tbx_query
#>              chr     start       end  beta coverage            file
#>           <char>     <int>     <int> <num>    <int>          <char>
#>        1:   chr1    923949    923950 0.000        1 bisc_SRR6911624
#>        2:   chr1    923953    923954 0.000        1 bisc_SRR6911624
#>        3:   chr1    923959    923960 0.000        1 bisc_SRR6911624
#>        4:   chr1    923971    923972 0.000        1 bisc_SRR6911624
#>        5:   chr1    923973    923974 0.000        1 bisc_SRR6911624
#>       ---                                                          
#> 45733375:   chr4 190179369 190179370 0.000        2 bisc_SRR6911723
#> 45733376:   chr4 190179686 190179687 1.000        1 bisc_SRR6911723
#> 45733377:   chr4 190179687 190179688 0.500        2 bisc_SRR6911723
#> 45733378:   chr4 190179753 190179754 1.000        1 bisc_SRR6911723
#> 45733379:   chr4 190179754 190179755 0.333        3 bisc_SRR6911723
```

On 8 threads, retrieving 45,733,379 records across 100 files took 11.014
seconds.

### `summarize_regions`

To get a summary of the information of the gene bodies use
`summarize_regions`, providing the gene name column as the feature
column:

``` r

qt <- system.time(
  summary_query <- summarize_regions(
    bedfiles,
    regions,
    columns = 4,
    fun = c("sum", "mean", "range", "count"),
    col_names = "beta",
    feature_col = "gene"
  )
)
#> [15:00:22.484494] [iscream::summarize_regions] [info] Summarizing 5000 regions from 100 bedfiles
#> [15:00:22.484628] [iscream::summarize_regions] [info] using sum, mean, range, count
#> [15:00:22.484635] [iscream::summarize_regions] [info] with columns 4 as beta
head(summary_query)
#>       chr    start      end            file feature beta.sum beta.mean
#>    <char>    <int>    <int>          <char>  <char>    <num>     <num>
#> 1:   chr1  1471764  1497848 bisc_SRR6911624  ATAD3B   85.000 0.8585859
#> 2:   chr1  3069167  3438621 bisc_SRR6911624  PRDM16  723.500 0.5288743
#> 3:   chr1  2403963  2413797 bisc_SRR6911624   PEX10   15.000 0.2419355
#> 4:   chr1 10472287 10630758 bisc_SRR6911624   PEX14  198.000 0.7764706
#> 5:   chr1  2425979  2505532 bisc_SRR6911624   PLCH2  184.333 0.7228745
#> 6:   chr1  9292893  9369532 bisc_SRR6911624   SPSB1   64.000 0.8648649
#>    beta.range count
#>         <num> <num>
#> 1:          1    99
#> 2:          1  1368
#> 3:          1    62
#> 4:          1   255
#> 5:          1   255
#> 6:          1    74
```

On 8 threads, summarizing the 45,733,379 records across all 100 files
took 4.683 seconds. As with
[`tabix()`](https://huishenlab.github.io/iscream/reference/tabix.md),
runtime reduces and memory usage increases with increasing thread
counts.

### `make_mat`

The
[`make_mat()`](https://huishenlab.github.io/iscream/reference/make_mat.md)
function queries and stores every genomic location within the input
regions across input files. Unlike
[`summarize_regions()`](https://huishenlab.github.io/iscream/reference/summarize_regions.md)
the output matrix dimensions are unknown at runtime. Although usually
quite fast, if the number of records falling into the input regions are
large and there are few overlaps between files, this can take a long
time. Here, gene bodies are large and the final matrix can contain
millions of CpG loci. Further, with sparse data, chances are new loci
are found in every file.

Preallocating the number of rows, however, can drastically reduce
runtime. Since we got 45 million loci/CpGs from all the BED files with
the tabix query above, we can expect approximately between 5 and 10
million unique loci as single-cell WGBS data has lower coverage than
bulk. We already have the tabix query so we can get the unique CpG count
here and use it to preallocate the matrix, reducing the number of matrix
resizes. We’ll add 100,000 extra rows on top of the existing count to be
safe since every avoided resize cuts off at least a couple seconds from
the runtime. Making tabix queries can be a relatively quick way to
approximate the CpG count of a dataset. If you haven’t done a tabix
query of the full dataset, you can approximate how many CpGs to expect
based on CpG counts in one file and the coverage of your WGBS method.
Here we make a matrix of the beta-values in the 4th column:

``` r

suppressPackageStartupMessages(library("SummarizedExperiment"))
cpg.count <- tbx_query$start |> unique() |> length()
qt <- system.time(meth_mat <- make_mat_se(
  bedfiles,
  regions,
  column = 4,
  sparse = TRUE,
  prealloc = cpg.count + 1e5
))
#> [15:00:32.589487] [iscream::query_all] [info] Querying 5000 regions from 100 bedfiles
#> 
#> [15:01:12.454511] [iscream::query_all] [info] Creating metadata vectors
#> [15:01:12.911981] [iscream::query_all] [info] 7276107 loci found - 16250 extra rows allocated with 0 resizes
#> [15:01:17.720398] [iscream::query_all] [info] Creating sparse matrix
meth_mat
#> class: RangedSummarizedExperiment 
#> dim: 7276107 100 
#> metadata(0):
#> assays(1): value
#> rownames: NULL
#> rowData names(0):
#> colnames(100): bisc_SRR6911624 bisc_SRR6911625 ... bisc_SRR6911722
#>   bisc_SRR6911723
#> colData names(0):
```

Making this 7,276,107 x 100 matrix took 45.812 seconds.

## Session info

``` r

sessionInfo()
#> R version 4.5.0 (2025-04-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: AlmaLinux 9.6 (Sage Margay)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.utf8      LC_MONETARY=C.UTF-8    LC_MESSAGES=C.utf8    
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [3] GenomicRanges_1.62.0        Seqinfo_1.0.0              
#>  [5] IRanges_2.44.0              S4Vectors_0.48.0           
#>  [7] BiocGenerics_0.56.0         generics_0.1.4             
#>  [9] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [11] iscream_1.1.6               data.table_1.17.8          
#> [13] BiocFileCache_3.0.0         dbplyr_2.5.1               
#> 
#> loaded via a namespace (and not attached):
#>  [1] rappdirs_0.3.3        SparseArray_1.10.1    RSQLite_2.4.3        
#>  [4] lattice_0.22-7        magrittr_2.0.4        evaluate_1.0.5       
#>  [7] grid_4.5.0            fastmap_1.2.0         blob_1.2.4           
#> [10] Matrix_1.7-4          DBI_1.2.3             purrr_1.1.0          
#> [13] pbapply_1.7-4         httr2_1.2.1           abind_1.4-8          
#> [16] cli_3.6.5             rlang_1.1.6           XVector_0.50.0       
#> [19] parallelly_1.45.1     bit64_4.6.0-1         DelayedArray_0.36.0  
#> [22] withr_3.0.2           cachem_1.1.0          S4Arrays_1.10.0      
#> [25] tools_4.5.0           parallel_4.5.0        memoise_2.0.1        
#> [28] dplyr_1.1.4           filelock_1.0.3        curl_7.0.0           
#> [31] vctrs_0.6.5           R6_2.6.1              lifecycle_1.0.4      
#> [34] stringfish_0.17.0     bit_4.6.0             pkgconfig_2.0.3      
#> [37] RcppParallel_5.1.11-1 pillar_1.11.1         glue_1.8.0           
#> [40] Rcpp_1.1.0            xfun_0.54             tibble_3.3.0         
#> [43] tidyselect_1.2.1      knitr_1.50            compiler_4.5.0
```

## References

Luo, Chongyuan, Angeline Rivkin, Jingtian Zhou, et al. 2018. “Robust
Single-Cell DNA Methylome Profiling with snmC-seq2.” *Nat Commun* 9 (1):
3824. <https://doi.org/10.1038/s41467-018-06355-2>.
