# Apply grouping to intervals

Apply grouping to intervals

## Usage

``` r
apply_interval_groups(
  df_list,
  group_channels,
  metadata,
  group_intervals,
  zero_time = FALSE,
  verbose = TRUE,
  env = rlang::caller_env()
)
```

## Arguments

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
