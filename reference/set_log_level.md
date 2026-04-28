# Set and get logging level

Set and get logging level

## Usage

``` r
set_log_level(level = "info")

get_log_level()
```

## Arguments

- level:

  The logging verbosity level to use

  - `"info"`: the default that gives provides basic information about
    the number of files and regions used in a function

  - `"debug"`: more verbose about row allocations, how many CpGs were
    found in a region, filename parsing etc. This mode cannot be used on
    more than one thread as R cannot output messages from multiple
    threads without crashing.

  - `"off"`: no logging

## Value

- `set_log_level()`: None; sets the log level to the provided level

- `get_log_level()`: The current logging level as a string

## Examples

``` r
set_log_level("info")
get_log_level()
#> [1] "info"
```
