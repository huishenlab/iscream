# Validate provided logging level

Only "info" and "debug" are currently supported, with "debug" only
supported when using 1 thread

## Usage

``` r
validate_log_level(level = get_log_level(), n_threads)
```

## Arguments

- level:

  The logging level to validate

- n_threads:

  The number of threads that the next iscream function call will use

## Value

None; sets the log level to the provide level

## Examples

``` r
set_log_level("info")
```
