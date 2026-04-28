# spdlog Logging Lever Setter

A helper function to turn a logging level given as string into the
current logging level

## Usage

``` r
Cpp_set_log_level(name)
```

## Arguments

- name:

  A string with the logging level. Value understood are, in decreasing
  verbosity ‘trace’, ‘debug’, ‘info’, ‘warning’, ‘error’, ‘critical’,
  and ‘off’. Unrecognised names are equivalent to ‘off’.

## Value

Nothing is returned.
