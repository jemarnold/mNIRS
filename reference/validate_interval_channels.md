# Validate per-group channel selections for ensemble-averaging

Normalises `group_channels` to a list of channel-name vectors for
recycling across interval groups. `NULL` selects all `nirs_channels` for
every group. Unlike
[`validate_group_channels()`](https://jemarnold.github.io/mnirs/reference/validate_group_channels.md),
channels may repeat across list items: each item is one group's channel
selection.

## Usage

``` r
validate_interval_channels(
  group_channels,
  nirs_channels,
  data = NULL,
  env = rlang::caller_env()
)
```

## Arguments

- nirs_channels:

  A character vector giving the names of mNIRS columns to operate on.
  Must match column names in `data` exactly.

  - If `NULL` (default), the `nirs_channels` metadata attribute of
    `data` is used.

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
