# Validate and Estimate Sample Rate

Validate and Estimate Sample Rate

## Usage

``` r
parse_sample_rate(
  data,
  file_header,
  time_channel,
  sample_rate = NULL,
  nirs_device = NULL,
  verbose = TRUE,
  env = rlang::caller_env()
)
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata.

- time_channel:

  A character string naming the time or sample column. Must match a
  column name in `data` exactly.

  - If `NULL` (default), the `time_channel` metadata attribute of `data`
    is used.

- sample_rate:

  A numeric sample rate in Hz.

  - If `NULL` (default), the `sample_rate` metadata attribute of `data`
    will be used if detected, or the sample rate will be estimated from
    `time_channel`.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
