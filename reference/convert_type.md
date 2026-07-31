# Coerce column types by role: nirs numeric, event integer, others detected

Coerce column types by role: nirs numeric, event integer, others
detected

## Usage

``` r
convert_type(
  data,
  nirs_channels = NULL,
  time_channel,
  event_channel = NULL,
  verbose = TRUE,
  env = rlang::caller_env()
)
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata.

- nirs_channels:

  A character vector giving the names of mNIRS columns to operate on.
  Must match column names in `data` exactly.

  - If `NULL` (default), the `nirs_channels` metadata attribute of
    `data` is used.

- time_channel:

  A character string naming the time or sample column. Must match a
  column name in `data` exactly.

  - If `NULL` (default), the `time_channel` metadata attribute of `data`
    is used.

- event_channel:

  A character string naming the event/lap column. Must match a column
  name in `data` exactly.

  - If `NULL` (default), the `event_channel` metadata attribute of
    `data` is used.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
