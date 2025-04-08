## iscream <a href="https://huishenlab.github.io/iscream/"><img src="man/figures/logo.png" align="right" height="138" style="float:right; height:138px;"/></a>

*Integrating Single-Cell Results for Exploring and Analyzing Methylation*

<!-- badges: start -->
[![R-CMD-check](https://github.com/huishenlab/iscream/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/huishenlab/iscream/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/huishenlab/iscream/graph/badge.svg)](https://app.codecov.io/gh/huishenlab/iscream)
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

## Dependencies

### *htslib* >= 1.17

*iscream* depends on the *htslib* header files. These may be installed with
your package manager:

- ubuntu/debian: `libhts-dev`  
- fedora/RHEL: `htslib-devel`  
- brew: `htslib`  
- nixpkgs: `htslib`
- conda: `bioconda::htslib`

or built manually: <https://www.htslib.org/download/>. See
[htslib.md](htslib.md) for more information.

The header files may also be found among your HPC modules - make sure the
`PKG_CONFIG_PATH` environment variable includes the `pkgconfig` location for
your installation of *htslib*. You can verify that the *htslib* development
libraries are installed with `pkg-config`:

```bash
pkg-config --cflags --libs htslib
```

#### *tabix*

Some *htslib* installations do not include the *tabix* executable (on Ubuntu you
need to install both *libhts-dev* and *tabix*). *iscream* will work without
*tabix*, but the `tabix()` function will be faster if the executable is
installed.

### GCC >= 9.4.0

GNU GCC must be installed for OpenMP support. This is usually installed by
default on Linux systems, but may need to be manually installed on MacOS to use
iscream with multiple threads[^1].

[^1]: Using OpenMP is also possible with Clang on MacOS
(<https://mac.r-project.org/openmp/>) but installing GCC with Homebrew may be
easier (<https://formulae.brew.sh/formula/gcc>).

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

### Usage

A user guide is available on the [package website](https://huishenlab.github.io/iscream/).
Bug reports may be submitted through GitHub issues.
