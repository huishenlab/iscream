# Get the number of available threads

Gets the number of threads iscream is currently set to use, whether the
`"iscream.threads"` option is set and how many threads are available for
use. To set the number of threads use
[`set_threads()`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
or set the `iscream.threads` option in your `~/.Rprofile`. See
[`?set_threads`](https://huishenlab.github.io/iscream/dev/reference/set_threads.md)
for more information.

## Usage

``` r
get_threads()
```

## Value

A named vector:

- `use_threads` = the number of threads iscream will use

- `opt_set` = whether the option was set by the user

- `avail_threads` = The number of available threads as reported by
  [`parallelly::availableCores`](https://parallelly.futureverse.org/reference/availableCores.html)

## Examples

``` r
get_threads()
#>   use_threads       opt_set avail_threads 
#>             1             1             4 
```
