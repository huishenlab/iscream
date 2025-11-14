# Check that the required threads are available

Check that the required threads are available

## Usage

``` r
check_thread_count(
  n_threads,
  avail_threads = availableCores(),
  opt_set = FALSE
)
```

## Arguments

- n_threads:

  The number of threads to check availability for

- avail_threads:

  The number of threads that are available on the system. Defaults to
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)

- opt_set:

  Whether the `iscream.threads` options is set

## Value

`n_threads` if the requested number of threads are available and stops
if not
