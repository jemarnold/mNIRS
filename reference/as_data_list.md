# Coerce `data` input to a named list of data frames

Coerce `data` input to a named list of data frames

## Usage

``` r
as_data_list(data, env = rlang::caller_env())
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
