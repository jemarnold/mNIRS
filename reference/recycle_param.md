# Recycle parameter to match number of events

Recycle an argument vector to a list or repeat the last list item to
match the number of events.

## Usage

``` r
recycle_param(
  param,
  n_events,
  group_intervals,
  verbose = TRUE,
  env = rlang::caller_env()
)
```
