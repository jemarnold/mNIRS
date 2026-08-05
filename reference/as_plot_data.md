# Validate and bind a list of mnirs data frames for plotting

Validate and bind a list of mnirs data frames for plotting

## Usage

``` r
as_plot_data(x, env = rlang::caller_env())
```

## Arguments

- x:

  A numeric vector.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.

## Value

For a single-element list, that element unchanged. Otherwise a row-bound
`data.frame` with an `interval` factor column, carrying attributes
`nirs_channels` (the union across elements), `time_channel`, and
`channel_map` — a named list mapping each channel to the interval names
whose source element declares it, so
[`plot.mnirs()`](https://jemarnold.github.io/mnirs/reference/plot.mnirs.md)
draws each channel only in its own panels.
