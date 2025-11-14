# An introduction to iscream

## Introduction

iscream is a BED file querying package that can can retrieve records
within regions of interest from multiple
[tabixed](https://en.wikipedia.org/wiki/Tabix) BED files. It can make
queries like tabix, summarize data within regions and make matrices
across input files and regions. All operations can be done in parallel
and return R/Bioconductor data structures.

iscream was designed to be memory efficient and fast on large datasets.
This vignette demonstrates iscream’s functions on a small dataset. For
more information on performance and considerations on larger real
datasets see
[`vignette("performance")`](https://huishenlab.github.io/iscream/dev/articles/performance.md).

## Installation

iscream may be installed from <https://bioconductor.org> with

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("iscream")
```

Although iscream will compile using `Rhtslib`, it will perform best if
compiled against [htslib](https://www.htslib.org/download/) which has
itself been compiled with `libdeflate`. See
[`vignette("htslib")`](https://huishenlab.github.io/iscream/dev/articles/htslib.md)
for more information.

### Loading iscream

``` r
library(iscream)
```

    ## iscream using 1 thread by default but parallelly::availableCores() detects 4 possibly available threads. See `?set_threads` for information on multithreading before trying to use more.

``` r
set_threads(2)
```

    ## iscream now using 2 of 4 available threads.

On load, iscream will inform the user about the number if threads it is
set to use. This is configurable with either
[`set_threads()`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
or by setting `option("iscream.threads")`, which can be done in your
`.Rprofile` to set it on startup. Once loaded you use the following
querying functions.

## `tabix()`

tabix is a command-line tool that queries records from compressed BED
files. As input it receives the path to a tabixed BED file and one or
more genomic regions of interest. It queries and prints the records from
the tabixed BED file that fall within the regions of interest. This
output may be redirected to a file and read into R.

iscream’s
[`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
function works in a similar way but supports multiple files and multiple
input formats all in one function. It receives a vector of paths to
tabixed BED files and regions of interest. The input regions may be a
vector of `"chr:start-end"` strings, a data frame with `"chr"`,
`"start"`, and `"end"` columns or a `GRanges` object. It queries the BED
files and returns the records as a parsed
*[data.table](https://CRAN.R-project.org/package=data.table)*. To get a
`GRanges` object instead, use `tabix_gr`.

Using iscream eliminates context switching between the shell and R for
tabix queries since the queries are made in R and the result is an R
data structure.
*[Rsamtools](https://bioconductor.org/packages/3.22/Rsamtools)* is
another R package that can query data from BED files, but it only
supports one file at a time, only accepts `GRanges` as the input regions
and does not parse the tab-delimited strings. See
[`vignette("tabix")`](https://huishenlab.github.io/iscream/dev/articles/tabix.md)
for a comparison of iscream and Rsamtools.

### Setup

Using [`list.files()`](https://rdrr.io/r/base/list.files.html), you can
get the paths to all compressed bed files in a directory. These files
contain small regions from chromosome 1 and Y from the snmC-seq2
methylation data ([Luo et al. 2018](#ref-luo2018a)).

``` r
data_dir <- system.file("extdata", package = "iscream")
(bedfiles <- list.files(
  data_dir,
  pattern = "cell[1-5].bed.gz$",
  full.names = TRUE))
```

    ## [1] "/home/runner/work/_temp/Library/iscream/extdata/cell1.bed.gz"
    ## [2] "/home/runner/work/_temp/Library/iscream/extdata/cell2.bed.gz"
    ## [3] "/home/runner/work/_temp/Library/iscream/extdata/cell3.bed.gz"
    ## [4] "/home/runner/work/_temp/Library/iscream/extdata/cell4.bed.gz"

For demonstration, I’m using two regions.

``` r
regions <- c("chr1:184577-680065", "chrY:56877780-56882524")
```

**Note**: in this dataset the chromosome identifier is in a
“chr\[number\]” format. However, in some datasets, the chromosome ID is
just a number or a letter like ‘1’, or ‘Y’. To see the chromosome ID
format across your files you can use
[`query_chroms()`](https://huishenlab.github.io/iscream/dev/reference/query_chroms.md)
to get all unique chromosomes across all input files:

``` r
query_chroms(bedfiles)
```

    ## [1] "chr1" "chrY"

### Making queries

[`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
one file:

``` r
tabix(bedfiles[1], regions)
```

    ##        chr    start      end    V1    V2
    ##     <char>    <int>    <int> <num> <int>
    ##  1:   chr1   184577   184578     0     1
    ##  2:   chr1   191223   191224     0     1
    ##  3:   chr1   191264   191265     0     1
    ##  4:   chr1   191307   191308     0     1
    ##  5:   chr1   191491   191492     0     1
    ##  6:   chr1   191507   191508     0     1
    ##  7:   chr1   191526   191527     0     1
    ##  8:   chr1   679950   679951     1     1
    ##  9:   chr1   680064   680065     1     1
    ## 10:   chrY 56877780 56877781     0     1
    ## 11:   chrY 56878350 56878351     0     1
    ## 12:   chrY 56878351 56878352     0     1
    ## 13:   chrY 56878431 56878432     0     1
    ## 14:   chrY 56879483 56879484     0     2
    ## 15:   chrY 56879535 56879536     0     2
    ## 16:   chrY 56879539 56879540     1     1
    ## 17:   chrY 56879724 56879725     0     1
    ## 18:   chrY 56879834 56879835     0     1
    ## 19:   chrY 56879863 56879864     1     1
    ## 20:   chrY 56879915 56879916     1     1
    ## 21:   chrY 56879945 56879946     1     1
    ## 22:   chrY 56881130 56881131     1     1
    ## 23:   chrY 56881926 56881927     1     1
    ## 24:   chrY 56882523 56882524     0     1
    ##        chr    start      end    V1    V2

With multiple files, the output contains a column for the file name.

``` r
tabix(bedfiles, regions)
```

    ##        chr    start      end    V1    V2   file
    ##     <char>    <int>    <int> <num> <int> <char>
    ##  1:   chr1   184577   184578     0     1  cell1
    ##  2:   chr1   191223   191224     0     1  cell1
    ##  3:   chr1   191264   191265     0     1  cell1
    ##  4:   chr1   191307   191308     0     1  cell1
    ##  5:   chr1   191491   191492     0     1  cell1
    ## ---                                            
    ## 70:   chrY 56880518 56880519     0     1  cell4
    ## 71:   chrY 56880715 56880716     0     2  cell4
    ## 72:   chrY 56880914 56880915     0     2  cell4
    ## 73:   chrY 56881927 56881928     1     1  cell4
    ## 74:   chrY 56882345 56882346     1     1  cell4

### Using *[GenomicRanges](https://bioconductor.org/packages/3.22/GenomicRanges)*

Both
[`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
and
[`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
accept `GRanges` objects as input regions but
[`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
returns `GRanges` objects instead of data frames. `tabix_gr` will also
preserve any input `GRanges` metadata.

``` r
if (!require("GenomicRanges", quietly = TRUE)) {
  stop("The 'GenomicRanges' package must be installed for this functionality")
}
gr <- GRanges(regions)
tabix_gr(bedfiles, gr)
```

    ## GRangesList object of length 4:
    ## $cell1
    ## GRanges object with 24 ranges and 3 metadata columns:
    ##        seqnames    ranges strand |        V1        V2        file
    ##           <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
    ##    [1]     chr1    184578      * |         0         1       cell1
    ##    [2]     chr1    191224      * |         0         1       cell1
    ##    [3]     chr1    191265      * |         0         1       cell1
    ##    [4]     chr1    191308      * |         0         1       cell1
    ##    [5]     chr1    191492      * |         0         1       cell1
    ##    ...      ...       ...    ... .       ...       ...         ...
    ##   [20]     chrY  56879916      * |         1         1       cell1
    ##   [21]     chrY  56879946      * |         1         1       cell1
    ##   [22]     chrY  56881131      * |         1         1       cell1
    ##   [23]     chrY  56881927      * |         1         1       cell1
    ##   [24]     chrY  56882524      * |         0         1       cell1
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
    ## 
    ## ...
    ## <3 more elements>

### Setting column names

You can set the result data frame’s column names or the `GRanges`
`mcols` with `col.names`:

``` r
tabix_gr(bedfiles, regions, col.names = c("beta", "coverage"))
```

    ## GRangesList object of length 4:
    ## $cell1
    ## GRanges object with 24 ranges and 3 metadata columns:
    ##        seqnames    ranges strand |      beta  coverage        file
    ##           <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
    ##    [1]     chr1    184578      * |         0         1       cell1
    ##    [2]     chr1    191224      * |         0         1       cell1
    ##    [3]     chr1    191265      * |         0         1       cell1
    ##    [4]     chr1    191308      * |         0         1       cell1
    ##    [5]     chr1    191492      * |         0         1       cell1
    ##    ...      ...       ...    ... .       ...       ...         ...
    ##   [20]     chrY  56879916      * |         1         1       cell1
    ##   [21]     chrY  56879946      * |         1         1       cell1
    ##   [22]     chrY  56881131      * |         1         1       cell1
    ##   [23]     chrY  56881927      * |         1         1       cell1
    ##   [24]     chrY  56882524      * |         0         1       cell1
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
    ## 
    ## ...
    ## <3 more elements>

### Rsamtools-style output

[`tabix_raw()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)
will return an unparsed named list like
*[Rsamtools](https://bioconductor.org/packages/3.22/Rsamtools)*, but
with multi-file support:

``` r
tabix_raw(bedfiles, regions)
```

    ## $cell1
    ## $cell1$`chr1:184577-680065`
    ## [1] "chr1\t184577\t184578\t0.000\t1" "chr1\t191223\t191224\t0.000\t1"
    ## [3] "chr1\t191264\t191265\t0.000\t1" "chr1\t191307\t191308\t0.000\t1"
    ## [5] "chr1\t191491\t191492\t0.000\t1" "chr1\t191507\t191508\t0.000\t1"
    ## [7] "chr1\t191526\t191527\t0.000\t1" "chr1\t679950\t679951\t1.000\t1"
    ## [9] "chr1\t680064\t680065\t1.000\t1"
    ## 
    ## $cell1$`chrY:56877780-56882524`
    ##  [1] "chrY\t56877780\t56877781\t0.000\t1" "chrY\t56878350\t56878351\t0.000\t1"
    ##  [3] "chrY\t56878351\t56878352\t0.000\t1" "chrY\t56878431\t56878432\t0.000\t1"
    ##  [5] "chrY\t56879483\t56879484\t0.000\t2" "chrY\t56879535\t56879536\t0.000\t2"
    ##  [7] "chrY\t56879539\t56879540\t1.000\t1" "chrY\t56879724\t56879725\t0.000\t1"
    ##  [9] "chrY\t56879834\t56879835\t0.000\t1" "chrY\t56879863\t56879864\t1.000\t1"
    ## [11] "chrY\t56879915\t56879916\t1.000\t1" "chrY\t56879945\t56879946\t1.000\t1"
    ## [13] "chrY\t56881130\t56881131\t1.000\t1" "chrY\t56881926\t56881927\t1.000\t1"
    ## [15] "chrY\t56882523\t56882524\t0.000\t1"
    ## 
    ## 
    ## $cell2
    ## $cell2$`chr1:184577-680065`
    ## character(0)
    ## 
    ## $cell2$`chrY:56877780-56882524`
    ##  [1] "chrY\t56878351\t56878352\t0.000\t1" "chrY\t56879482\t56879483\t1.000\t1"
    ##  [3] "chrY\t56879483\t56879484\t1.000\t1" "chrY\t56879534\t56879535\t1.000\t1"
    ##  [5] "chrY\t56879538\t56879539\t1.000\t1" "chrY\t56880322\t56880323\t0.000\t1"
    ##  [7] "chrY\t56880404\t56880405\t1.000\t1" "chrY\t56881129\t56881130\t0.000\t1"
    ##  [9] "chrY\t56881285\t56881286\t0.000\t1" "chrY\t56881364\t56881365\t1.000\t1"
    ## [11] "chrY\t56882344\t56882345\t0.000\t1" "chrY\t56882523\t56882524\t0.000\t1"
    ## 
    ## 
    ## $cell3
    ## $cell3$`chr1:184577-680065`
    ## character(0)
    ## 
    ## $cell3$`chrY:56877780-56882524`
    ## [1] "chrY\t56879862\t56879863\t1.000\t1" "chrY\t56879914\t56879915\t0.500\t2"
    ## [3] "chrY\t56879944\t56879945\t1.000\t1" "chrY\t56879969\t56879970\t1.000\t1"
    ## [5] "chrY\t56879976\t56879977\t1.000\t1" "chrY\t56880322\t56880323\t0.000\t1"
    ## [7] "chrY\t56880359\t56880360\t0.000\t1" "chrY\t56881129\t56881130\t0.500\t2"
    ## 
    ## 
    ## $cell4
    ## $cell4$`chr1:184577-680065`
    ##  [1] "chr1\t191003\t191004\t1.000\t1" "chr1\t191307\t191308\t1.000\t1"
    ##  [3] "chr1\t191438\t191439\t1.000\t1" "chr1\t191491\t191492\t1.000\t1"
    ##  [5] "chr1\t191722\t191723\t1.000\t1" "chr1\t191747\t191748\t1.000\t1"
    ##  [7] "chr1\t191776\t191777\t1.000\t1" "chr1\t628458\t628459\t1.000\t1"
    ##  [9] "chr1\t631433\t631434\t1.000\t1" "chr1\t631436\t631437\t1.000\t1"
    ## 
    ## $cell4$`chrY:56877780-56882524`
    ##  [1] "chrY\t56877780\t56877781\t0.000\t1" "chrY\t56877810\t56877811\t0.000\t1"
    ##  [3] "chrY\t56879834\t56879835\t0.000\t1" "chrY\t56879863\t56879864\t1.000\t1"
    ##  [5] "chrY\t56879915\t56879916\t1.000\t1" "chrY\t56879945\t56879946\t0.000\t1"
    ##  [7] "chrY\t56880045\t56880046\t0.500\t2" "chrY\t56880081\t56880082\t0.000\t1"
    ##  [9] "chrY\t56880092\t56880093\t0.500\t2" "chrY\t56880112\t56880113\t1.000\t1"
    ## [11] "chrY\t56880137\t56880138\t1.000\t2" "chrY\t56880142\t56880143\t0.500\t2"
    ## [13] "chrY\t56880323\t56880324\t0.000\t1" "chrY\t56880360\t56880361\t1.000\t1"
    ## [15] "chrY\t56880405\t56880406\t1.000\t1" "chrY\t56880518\t56880519\t0.000\t1"
    ## [17] "chrY\t56880715\t56880716\t0.000\t2" "chrY\t56880914\t56880915\t0.000\t2"
    ## [19] "chrY\t56881927\t56881928\t1.000\t1" "chrY\t56882345\t56882346\t1.000\t1"

### WGBS BED files

Since iscream was originally developed to read Whole Genome Bisulfite
Sequencing data, it has the `aligner` argument for automatically setting
column names for the BISCUIT ([Zhou et al. 2024](#ref-zhou2024)),
Bismark ([Krueger and Andrews 2011](#ref-krueger2011)), and BSBolt
([Farrell et al. 2021](#ref-farrell2021)).

``` r
tabix_gr(bedfiles, regions, aligner = "biscuit")
```

    ## GRangesList object of length 4:
    ## $cell1
    ## GRanges object with 24 ranges and 3 metadata columns:
    ##        seqnames    ranges strand |      beta  coverage        file
    ##           <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
    ##    [1]     chr1    184578      * |         0         1       cell1
    ##    [2]     chr1    191224      * |         0         1       cell1
    ##    [3]     chr1    191265      * |         0         1       cell1
    ##    [4]     chr1    191308      * |         0         1       cell1
    ##    [5]     chr1    191492      * |         0         1       cell1
    ##    ...      ...       ...    ... .       ...       ...         ...
    ##   [20]     chrY  56879916      * |         1         1       cell1
    ##   [21]     chrY  56879946      * |         1         1       cell1
    ##   [22]     chrY  56881131      * |         1         1       cell1
    ##   [23]     chrY  56881927      * |         1         1       cell1
    ##   [24]     chrY  56882524      * |         0         1       cell1
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
    ## 
    ## ...
    ## <3 more elements>

## `summarize_regions()`

Using the same BED files and regions, iscream can also summarize the
data within regions. This requires providing the index of the data
columns you want summarized and the column names for the output. All
summarizing functions are run on each input column by default - see
[`?summarize_regions`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md)
for supported functions.

``` r
summarize_regions(
  bedfiles,
  regions,
  columns = c(4, 5),
  col_names = c("beta", "coverage")
)
```

    ## [20:40:22.932336] [iscream::summarize_regions] [info] Summarizing 2 regions from 4 bedfiles
    ## [20:40:22.932380] [iscream::summarize_regions] [info] using sum, mean, median, stddev, variance, min, max, range, count
    ## [20:40:22.932385] [iscream::summarize_regions] [info] with columns 4, 5 as beta, coverage

    ##                  feature  file beta.sum coverage.sum beta.mean coverage.mean
    ## 1     chr1:184577-680065 cell1      2.0            9 0.2222222      1.000000
    ## 2 chrY:56877780-56882524 cell1      6.0           17 0.4000000      1.133333
    ## 3     chr1:184577-680065 cell2       NA           NA        NA            NA
    ## 4 chrY:56877780-56882524 cell2      6.0           12 0.5000000      1.000000
    ## 5     chr1:184577-680065 cell3       NA           NA        NA            NA
    ## 6 chrY:56877780-56882524 cell3      5.0           10 0.6250000      1.250000
    ## 7     chr1:184577-680065 cell4     10.0           10 1.0000000      1.000000
    ## 8 chrY:56877780-56882524 cell4      9.5           26 0.4750000      1.300000
    ##   beta.median coverage.median beta.stddev coverage.stddev beta.variance
    ## 1        0.00               1   0.4409586       0.0000000     0.1944444
    ## 2        0.00               1   0.5070926       0.3518658     0.2571429
    ## 3          NA              NA          NA              NA            NA
    ## 4        0.50               1   0.5222330       0.0000000     0.2727273
    ## 5          NA              NA          NA              NA            NA
    ## 6        0.75               1   0.4432026       0.4629100     0.1964286
    ## 7        1.00               1   0.0000000       0.0000000     0.0000000
    ## 8        0.50               1   0.4722566       0.4701623     0.2230263
    ##   coverage.variance beta.min coverage.min beta.max coverage.max beta.range
    ## 1         0.0000000        0            1        1            1          1
    ## 2         0.1238095        0            1        1            2          1
    ## 3                NA       NA           NA       NA           NA         NA
    ## 4         0.0000000        0            1        1            1          1
    ## 5                NA       NA           NA       NA           NA         NA
    ## 6         0.2142857        0            1        1            2          1
    ## 7         0.0000000        1            1        1            1          0
    ## 8         0.2210526        0            1        1            2          1
    ##   coverage.range count
    ## 1              0     9
    ## 2              1    15
    ## 3             NA    NA
    ## 4              0    12
    ## 5             NA    NA
    ## 6              1     8
    ## 7              0    10
    ## 8              1    20

The `feature` column here contains the genomic region coordinates, but
can be set to something more informational if you have names for the
regions. You can also select the functions applied with `fun`:

``` r
names(regions) <- c("R1", "R2")
summarize_regions(
  bedfiles,
  regions,
  fun = c("mean", "sum"),
  columns = 5,
  col_names = "coverage"
)
```

    ## [20:40:23.018737] [iscream::summarize_regions] [info] Summarizing 2 regions from 4 bedfiles
    ## [20:40:23.018756] [iscream::summarize_regions] [info] using mean, sum
    ## [20:40:23.018760] [iscream::summarize_regions] [info] with columns 5 as coverage

    ##   feature  file coverage.mean coverage.sum
    ## 1      R1 cell1      1.000000            9
    ## 2      R2 cell1      1.133333           17
    ## 3      R1 cell2            NA           NA
    ## 4      R2 cell2      1.000000           12
    ## 5      R1 cell3            NA           NA
    ## 6      R2 cell3      1.250000           10
    ## 7      R1 cell4      1.000000           10
    ## 8      R2 cell4      1.300000           26

### WGBS BED files

For WGBS data specifically, you can use
[`summarize_meth_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_meth_regions.md),
which takes the `aligner` argument to correctly parse the columns.

``` r
summarize_meth_regions(
  bedfiles,
  regions,
  aligner = "biscuit",
  fun = c("mean", "sum")
)
```

    ## [20:40:23.087937] [iscream::summarize_regions] [info] Summarizing 2 regions from 4 bedfiles
    ## [20:40:23.087957] [iscream::summarize_regions] [info] using mean, sum
    ## [20:40:23.087961] [iscream::summarize_regions] [info] with columns 4, 5 as coverage, M

    ##   feature  file coverage.mean    M.mean coverage.sum M.sum
    ## 1      R1 cell1      1.000000 0.2222222            9     2
    ## 2      R2 cell1      1.133333 0.4000000           17     6
    ## 3      R1 cell2            NA        NA           NA    NA
    ## 4      R2 cell2      1.000000 0.5000000           12     6
    ## 5      R1 cell3            NA        NA           NA    NA
    ## 6      R2 cell3      1.250000 0.7500000           10     6
    ## 7      R1 cell4      1.000000 1.0000000           10    10
    ## 8      R2 cell4      1.300000 0.6000000           26    12

## `make_mat()`

[`make_mat()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
makes a matrix where every row is a locus found in at least one input
file and the columns are the input files. It returns a named list of the
matrix, coordinates, and file names you can use for analysis or other
data structures.
[`make_mat_se()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
returns a `RangedSummarizedExperiment`, and
[`make_mat_gr()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
returns a `GRanges` object.

``` r
if (!require("SummarizedExperiment", quietly = TRUE)) {
  stop("The 'SummarizedExperiment' package must be installed for this functionality")
}
(mat <- make_mat_se(bedfiles, regions, column = 4, mat_name = "beta"))
```

    ## [20:40:25.467901] [iscream::query_all] [info] Querying 2 regions from 4 bedfiles
    ## 
    ## [20:40:25.468396] [iscream::query_all] [info] Creating metadata vectors
    ## [20:40:25.468471] [iscream::query_all] [info] 62 loci found - 9938 extra rows allocated with 0 resizes
    ## [20:40:25.468480] [iscream::query_all] [info] Creating dense matrix

    ## class: RangedSummarizedExperiment 
    ## dim: 62 4 
    ## metadata(0):
    ## assays(1): beta
    ## rownames: NULL
    ## rowData names(0):
    ## colnames(4): cell1 cell2 cell3 cell4
    ## colData names(0):

``` r
head(assay(mat), 10)
```

    ##       cell1 cell2 cell3 cell4
    ##  [1,]     0     0   1.0     0
    ##  [2,]     0     0   0.5     0
    ##  [3,]     0     0   1.0     0
    ##  [4,]     0     0   1.0     0
    ##  [5,]     0     0   1.0     0
    ##  [6,]     0     0   0.0     0
    ##  [7,]     0     0   0.0     0
    ##  [8,]     0     0   0.5     0
    ##  [9,]     0     0   0.0     0
    ## [10,]     0     0   0.0     0

If you have sparse data, you can save memory with `sparse = TRUE`, but
only for
[`make_mat()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md)
and
[`make_mat_se()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md).

``` r
mat <- make_mat(bedfiles, regions, column = 4, mat_name = "beta", sparse = TRUE)
```

    ## [20:40:25.586620] [iscream::query_all] [info] Querying 2 regions from 4 bedfiles
    ## 
    ## [20:40:25.587081] [iscream::query_all] [info] Creating metadata vectors
    ## [20:40:25.587116] [iscream::query_all] [info] 62 loci found - 9938 extra rows allocated with 0 resizes
    ## [20:40:25.587128] [iscream::query_all] [info] Creating sparse matrix

``` r
head(mat$beta, 10)
```

    ## 10 x 4 sparse Matrix of class "dgCMatrix"
    ##       cell1 cell2 cell3 cell4
    ##  [1,]     .     .   1.0     .
    ##  [2,]     .     .   0.5     .
    ##  [3,]     .     .   1.0     .
    ##  [4,]     .     .   1.0     .
    ##  [5,]     .     .   1.0     .
    ##  [6,]     .     .   .       .
    ##  [7,]     .     .   .       .
    ##  [8,]     .     .   0.5     .
    ##  [9,]     .     .   .       .
    ## [10,]     .     .   .       .

### WGBS BED files

For a *[BSseq](https://bioconductor.org/packages/3.22/BSseq)* compatible
object, use `make_meth_mat` which returns a list of BSseq inputs. To
make the `BSseq` object run

``` r
if (require("bsseq", quietly = TRUE)) {
  meth_mat <- make_mat_bsseq(bedfiles, regions)
  do.call(BSseq, meth_mat)
}
```

    ## [20:40:29.453498] [iscream::query_all] [info] Querying 2 regions from 4 bedfiles
    ## 
    ## [20:40:29.453988] [iscream::query_all] [info] Creating metadata vectors
    ## [20:40:29.454024] [iscream::query_all] [info] 62 loci found - 9938 extra rows allocated with 0 resizes
    ## [20:40:29.454028] [iscream::query_all] [info] Creating dense matrix

    ## An object of type 'BSseq' with
    ##   62 methylation loci
    ##   4 samples
    ## has not been smoothed
    ## All assays are in-memory

## Further reading

While this vignette used a toy dataset,
[`vignette("performance")`](https://huishenlab.github.io/iscream/dev/articles/performance.md)
has tips on getting the best performance from iscream using a real and
larger dataset with 100 files.
[`vignette("TSS")`](https://huishenlab.github.io/iscream/dev/articles/TSS.md)
is an example of iscream’s application to examine methylation around
transcription start site on another real dataset. For more information
on the supported Bioconductor data structures and conversions see
[`vignette("data_structures")`](https://huishenlab.github.io/iscream/dev/articles/data_structures.md).

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
    ##  [5] matrixStats_1.5.0           GenomicRanges_1.62.0       
    ##  [7] Seqinfo_1.0.0               IRanges_2.44.0             
    ##  [9] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [11] generics_0.1.4              iscream_1.1.2              
    ## [13] BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] farver_2.1.2              R.utils_2.13.0           
    ##  [3] Biostrings_2.78.0         bitops_1.0-9             
    ##  [5] fastmap_1.2.0             RCurl_1.98-1.17          
    ##  [7] GenomicAlignments_1.46.0  stringfish_0.17.0        
    ##  [9] XML_3.99-0.20             digest_0.6.38            
    ## [11] lifecycle_1.0.4           statmod_1.5.1            
    ## [13] compiler_4.5.2            rlang_1.1.6              
    ## [15] sass_0.4.10               tools_4.5.2              
    ## [17] yaml_2.3.10               data.table_1.17.8        
    ## [19] rtracklayer_1.70.0        knitr_1.50               
    ## [21] S4Arrays_1.10.0           curl_7.0.0               
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
    ## [45] rhdf5_2.54.0              parallel_4.5.2           
    ## [47] BiocManager_1.30.26       XVector_0.50.0           
    ## [49] restfulr_0.0.16           Matrix_1.7-4             
    ## [51] jsonlite_2.0.0            bookdown_0.45            
    ## [53] systemfonts_1.3.1         h5mread_1.2.0            
    ## [55] locfit_1.5-9.12           limma_3.66.0             
    ## [57] jquerylib_0.1.4           glue_1.8.0               
    ## [59] parallelly_1.45.1         pkgdown_2.2.0            
    ## [61] codetools_0.2-20          BiocIO_1.20.0            
    ## [63] htmltools_0.5.8.1         rhdf5filters_1.22.0      
    ## [65] BSgenome_1.78.0           R6_2.6.1                 
    ## [67] textshaping_1.0.4         sparseMatrixStats_1.22.0 
    ## [69] evaluate_1.0.5            lattice_0.22-7           
    ## [71] R.methodsS3_1.8.2         Rsamtools_2.26.0         
    ## [73] cigarillo_1.0.0           bslib_0.9.0              
    ## [75] Rcpp_1.1.0                SparseArray_1.10.1       
    ## [77] permute_0.9-8             xfun_0.54                
    ## [79] fs_1.6.6

## References

Farrell, Colin, Michael Thompson, Anela Tosevska, Adewale Oyetunde, and
Matteo Pellegrini. 2021. “BiSulfite Bolt: A Bisulfite Sequencing
Analysis Platform.” *Gigascience* 10 (5): giab033.
<https://doi.org/10.1093/gigascience/giab033>.

Krueger, Felix, and Simon R. Andrews. 2011. “Bismark: A Flexible Aligner
and Methylation Caller for Bisulfite-Seq Applications.” *Bioinformatics*
27 (11): 1571–72. <https://doi.org/10.1093/bioinformatics/btr167>.

Luo, Chongyuan, Angeline Rivkin, Jingtian Zhou, Justin P. Sandoval,
Laurie Kurihara, Jacinta Lucero, Rosa Castanon, et al. 2018. “Robust
Single-Cell DNA Methylome Profiling with snmC-seq2.” *Nat Commun* 9 (1):
3824. <https://doi.org/10.1038/s41467-018-06355-2>.

Zhou, Wanding, Benjamin K Johnson, Jacob Morrison, Ian Beddows, James
Eapen, Efrat Katsman, Ayush Semwal, et al. 2024. “BISCUIT: An Efficient,
Standards-Compliant Tool Suite for Simultaneous Genetic and Epigenetic
Inference in Bulk and Single-Cell Studies.” *Nucleic Acids Research* 52
(6): gkae097. <https://doi.org/10.1093/nar/gkae097>.
