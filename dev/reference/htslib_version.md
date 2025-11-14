# Get htslib version and available features

Returns the version of htslib being used by iscream and whether features
such as libdeflate support are available. This information may not
always correspond to the htslib version used during iscream's
installation if a different htslib version is available for linking at
runtime.

## Usage

``` r
htslib_version()
```

## Value

None

## Examples

``` r
htslib_version()
#> 1.19
#> build=configure libcurl=yes S3=yes GCS=yes libdeflate=yes lzma=yes bzip2=yes plugins=yes plugin-path=/usr/local/lib/htslib:/usr/local/libexec/htslib:/usr/lib/x86_64-linux-gnu/htslib: htscodecs=1.6.0
```
