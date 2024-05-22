## iscream <a href="https://huishenlab.github.io/iscream/"><img src="man/figures/logo.png" align="right" height="138" style="float:right; height:138px;"/></a>

*Integrating Single-Cell Results for Exploring and Analyzing Methylation*

<!-- badges: start -->
[![R-CMD-check](https://github.com/huishenlab/iscream/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/huishenlab/iscream/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Analysis and visualization of Whole Genome Bisulfite Sequencing (WGBS) data
requires reading aligned sequencing data into formats that existing packages
like [BSseq](https://bioconductor.org/packages/devel/bioc/html/bsseq.html) and
[scMET](https://github.com/andreaskapou/scMET) can analyze. Getting the data
from on-disk formats like bedfiles to a matrix of methylation values can be
difficult because, with nearly 30 million CpGs, WGBS data can be quite large.

iscream aims to efficiently read alinged (sc)WGBS data into formats that can be
used by other packages. iscream uses [htslib](https://www.htslib.org/) to query
genomic regions to make matrices for BSSeq or aggregate the methylated reads for
scMET.

## Installation

### GitHub

You can install the development version from Github by cloning the repo and
running

```bash
git clone https://github.com/huishenlab/iscream
R CMD INSTALL iscream
```

You can also use the R [`devtools`](https://devtools.r-lib.org/) package:

```r
devtools::install_github("huishenlab/iscream")
```

or [`pak`](https://pak.r-lib.org/):

```r
pak::pkg_install("huishenlab/iscream")
```

### Dependencies

- htslib: check that you have the htslib dev package installed for your system.
  You can check this with `pkg-config` by running 

  ```bash
  pkg-config --exists htslib
  ```
### Usage

A user guide is available on the [package website](https://huishenlab.github.io/isceam/).
Bug reports may be submitted through GitHub issues.
