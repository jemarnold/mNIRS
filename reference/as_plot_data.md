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
