# Set the number of available threads

Sets the `"iscream.threads"` option to `n_threads`. To see how many
threads you have available see `?get_threads()`.

## Usage

``` r
set_threads(n_threads)
```

## Arguments

- n_threads:

  The number of threads to use

## Value

None. Sets the `iscream.threads` option to the requested number of
threads if available

## Details

iscream uses OpenMP to parallelize certain functions. You can use as
many threads as are available to you on your system to varying degrees
of performance improvements. The
[`get_threads()`](https://huishenlab.github.io/iscream/dev/reference/get_threads.md)
function uses
[`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
to report the number of available threads. Although OpenMP can detect
the number of available cores, on high performance computers (HPCs) with
resource allocating job schedulers like SLURM, OpenMP may detect all
available threads across the HPC and not limit itself to the cores that
were allocated to you by the scheduler. If your system administrator has
not set up any limits, this may result in your job taking resources from
other jobs. If there are limits, trying to use more threads that those
available will reduce iscream's performance. Job schedulers will
typically have an environment variable (e.g. `SLURM_CPUS_ON_NODE` with
SLURM) that gives you the actual number of available cores. Further, on
hyperthreaded systems, this count may be double that of the available
processors. Using hyperthreading does not guarantee any performance
improvement - it may be better to set the number of threads to half the
reported number.
[`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
takes HPC scheduler/CRAN/Bioconductor limits into account when reporting
the number of available threads but it may not reliably report
hyperthreading ('system' or 'nproc'). To set the number of threads
without having to call `set_threads()` in every session, put

    options(iscream.threads = [n_threads])

in your `.Rprofile` See
[`help('Rprofile')`](https://rdrr.io/r/base/Startup.html) for
information on startup options.

Functions currently using multithreading:

- [`tabix()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md),
  [`tabix_gr()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md),
  [`tabix_raw()`](https://huishenlab.github.io/iscream/dev/reference/tabix.md)

- [`query_chroms()`](https://huishenlab.github.io/iscream/dev/reference/query_chroms.md)

- [`make_mat()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md),
  [`make_mat_se()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md),
  [`make_mat_gr()`](https://huishenlab.github.io/iscream/dev/reference/make_mat.md),
  [`make_mat_bsseq()`](https://huishenlab.github.io/iscream/dev/reference/make_mat_bsseq.md)

- [`summarize_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_regions.md),
  [`summarize_meth_regions()`](https://huishenlab.github.io/iscream/dev/reference/summarize_meth_regions.md)

## Examples

``` r
(ncores <- parallelly::availableCores())
#> system 
#>      4 
# \donttest{
set_threads(ncores)
#> iscream now using 4 of 4 available threads.
# }
```
