# resolve start/end into time value vectors (no span applied)

resolve start/end into time value vectors (no span applied)

## Usage

``` r
resolve_interval(
  start,
  end,
  t_vec,
  event_vec = NULL,
  env = rlang::caller_env()
)
```

## Arguments

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
