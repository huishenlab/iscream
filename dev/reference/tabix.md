# Query records from tabixed BED files

Query records from tabixed BED files

## Usage

``` r
tabix(bedfiles, regions, aligner = NULL, col.names = NULL, nthreads = NULL)

tabix_gr(
  bedfiles,
  regions,
  aligner = NULL,
  col.names = NULL,
  zero_based = TRUE,
  nthreads = NULL
)

tabix_raw(bedfiles, regions, nthreads = NULL)
```

## Arguments

- bedfiles:

  The BED files to be queried

- regions:

  A vector, data frame or GenomicRanges of genomic regions. See details.

- aligner:

  The aligner used to produce the BED files - one of "biscuit",
  "bismark", "bsbolt". Will set the result data.table's column names
  based on this argument.

- col.names:

  A vector of column names for the data columns of the result.table, not
  including "chr", "start", and "end". Set if your BED file is not from
  the supported aligners or is a general BED file.

- nthreads:

  Set number of threads to use overriding the `"iscream.threads"`
  option. See
  [`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
  for more information.

- zero_based:

  Whether the input BED file has a zero-based start column - used when
  coverting the result data frame to GenomicRanges.

## Value

- `tabix()`: A data frame

- `tabix_gr()`: A `GRanges` object for single files and `GRangesList`
  for multiple files. When making `GRanges`, the 0-based records from
  BED-files will be converted to 1-based with
  [`GenomicRanges::makeGRangesFromDataFrame()`](https://rdrr.io/pkg/GenomicRanges/man/makeGRangesFromDataFrame.html).
  Bismark's coverage files will not be converted as they are already
  1-based and the `ranges` slot will be only one position.

- `tabix_raw()`: A named list of raw strings from the regions in the
  style of
  [`Rsamtools::scanTabix`](https://rdrr.io/pkg/Rsamtools/man/scanTabix.html)

## Details

### Query method

'*iscream* has two methods to query records from BED files:

- the *tabix* shell executable: fast since its output can be redirected
  to a file (which
  [`data.table::fread()`](https://rdatatable.gitlab.io/data.table/reference/fread.html)
  can then read) instead of having to allocate memory and store it
  during the query

- *iscream's* tabix implementation, based on the *tabix* executable
  using *htslib*, but slower on large queries since it stores the
  records as they are found instead of writing to a file. However it's
  able to store each region's records independently instead of in a
  single file and is used in
  [`make_mat()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md),
  [`make_mat_bsseq()`](https://huishenlab.github.io/iscream/dev/reference/make_mat_bsseq.md),
  and
  [`summarize_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md).

When *iscream* is attached, it checks that the *tabix* executable is
available with [`Sys.which()`](https://rdrr.io/r/base/Sys.which.html)
and, if available, sets `options("tabix.method" = "shell")`. `tabix()`
then uses the *tabix* executable to make queries, except for
`tabix_raw()`. If *tabix* is not found, *iscream* uses its tabix
implementation. To use only *iscream's* tabix implementation, set
`options("tabix.method" = "htslib")`.

### Input region formats

The input regions format may be string vector in the form
"chr:start-end", a dataframe with "chr", "start" and "end" columns or a
`GRanges` object. Input regions must be 1-based. When using `"htslib"`
as the query method, if the input `GRanges` object of regions contains
any single locus regions where the start and end positions are the same,
iscream will notify that such regions were found and fixed as
`chr:start` format strings are invalid for the htslib API (see
[`?get_granges_string`](https://huishenlab.github.io/iscream/dev/reference/get_granges_string.md)).

## Examples

``` r
bedfiles <- system.file("extdata", package = "iscream") |>
  list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
tabix(bedfiles, regions, col.names = c("beta", "coverage"))
#>        chr start   end  beta coverage   file
#>     <char> <int> <int> <num>    <int> <char>
#>  1:   chr1     0     2   1.0        1      a
#>  2:   chr1     2     4   1.0        1      a
#>  3:   chr1     4     6   0.0        2      a
#>  4:   chr1     6     8   0.0        1      a
#>  5:   chr1     8    10   0.5        2      a
#>  6:   chr1    10    12   1.0        2      a
#>  7:   chr1    12    14   1.0        3      a
#>  8:   chr1     0     2   0.0        2      b
#>  9:   chr1     4     6   1.0        2      b
#> 10:   chr1     6     8   1.0        1      b
#> 11:   chr1    10    12   0.0        2      b
#> 12:   chr1    12    14   1.0        1      b
#> 13:   chr1     2     4   1.0        2      c
#> 14:   chr1     6     8   0.0        2      c
#> 15:   chr1     8    10   1.0        1      c
#> 16:   chr1     0     2   1.0        1      d
#> 17:   chr1     2     4   1.0        2      d
#> 18:   chr1     6     8   0.0        1      d
#> 19:   chr1     8    10   0.5        2      d
#> 20:   chr1    12    14   1.0        1      d
#>        chr start   end  beta coverage   file
if (require("GenomicRanges", quietly = TRUE)) {
  tabix_gr(bedfiles, regions, col.names = c("beta", "coverage"))
}
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following objects are masked from ‘package:data.table’:
#> 
#>     first, second
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> Attaching package: ‘IRanges’
#> The following object is masked from ‘package:data.table’:
#> 
#>     shift
#> GRangesList object of length 4:
#> $a
#> GRanges object with 7 ranges and 3 metadata columns:
#>       seqnames    ranges strand |      beta  coverage        file
#>          <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
#>   [1]     chr1       1-2      * |       1.0         1           a
#>   [2]     chr1       3-4      * |       1.0         1           a
#>   [3]     chr1       5-6      * |       0.0         2           a
#>   [4]     chr1       7-8      * |       0.0         1           a
#>   [5]     chr1      9-10      * |       0.5         2           a
#>   [6]     chr1     11-12      * |       1.0         2           a
#>   [7]     chr1     13-14      * |       1.0         3           a
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#> $b
#> GRanges object with 5 ranges and 3 metadata columns:
#>       seqnames    ranges strand |      beta  coverage        file
#>          <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
#>   [1]     chr1       1-2      * |         0         2           b
#>   [2]     chr1       5-6      * |         1         2           b
#>   [3]     chr1       7-8      * |         1         1           b
#>   [4]     chr1     11-12      * |         0         2           b
#>   [5]     chr1     13-14      * |         1         1           b
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#> $c
#> GRanges object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |      beta  coverage        file
#>          <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
#>   [1]     chr1       3-4      * |         1         2           c
#>   [2]     chr1       7-8      * |         0         2           c
#>   [3]     chr1      9-10      * |         1         1           c
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#> $d
#> GRanges object with 5 ranges and 3 metadata columns:
#>       seqnames    ranges strand |      beta  coverage        file
#>          <Rle> <IRanges>  <Rle> | <numeric> <integer> <character>
#>   [1]     chr1       1-2      * |       1.0         1           d
#>   [2]     chr1       3-4      * |       1.0         2           d
#>   [3]     chr1       7-8      * |       0.0         1           d
#>   [4]     chr1      9-10      * |       0.5         2           d
#>   [5]     chr1     13-14      * |       1.0         1           d
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
tabix_raw(bedfiles, regions)
#> $a
#> $a$`chr1:1-6`
#> [1] "chr1\t0\t2\t1.000\t1" "chr1\t2\t4\t1.000\t1" "chr1\t4\t6\t0.000\t2"
#> 
#> $a$`chr1:7-10`
#> [1] "chr1\t6\t8\t0.000\t1"  "chr1\t8\t10\t0.500\t2"
#> 
#> $a$`chr1:11-14`
#> [1] "chr1\t10\t12\t1.000\t2" "chr1\t12\t14\t1.000\t3"
#> 
#> 
#> $b
#> $b$`chr1:1-6`
#> [1] "chr1\t0\t2\t0.000\t2" "chr1\t4\t6\t1.000\t2"
#> 
#> $b$`chr1:7-10`
#> [1] "chr1\t6\t8\t1.000\t1"
#> 
#> $b$`chr1:11-14`
#> [1] "chr1\t10\t12\t0.000\t2" "chr1\t12\t14\t1.000\t1"
#> 
#> 
#> $c
#> $c$`chr1:1-6`
#> [1] "chr1\t2\t4\t1.000\t2"
#> 
#> $c$`chr1:7-10`
#> [1] "chr1\t6\t8\t0.000\t2"  "chr1\t8\t10\t1.000\t1"
#> 
#> $c$`chr1:11-14`
#> character(0)
#> 
#> 
#> $d
#> $d$`chr1:1-6`
#> [1] "chr1\t0\t2\t1.000\t1" "chr1\t2\t4\t1.000\t2"
#> 
#> $d$`chr1:7-10`
#> [1] "chr1\t6\t8\t0.000\t1"  "chr1\t8\t10\t0.500\t2"
#> 
#> $d$`chr1:11-14`
#> [1] "chr1\t12\t14\t1.000\t1"
#> 
#> 
```
