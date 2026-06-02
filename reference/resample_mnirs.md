# Re-sample an *mnirs* data frame

Up- or down-sample an *"mnirs"* data frame to a new sample rate, filling
new samples via nearest-neighbour matching or interpolation.

## Usage

``` r
resample_mnirs(
  data,
  time_channel = NULL,
  sample_rate = NULL,
  resample_rate = sample_rate,
  method = c("none", "linear", "locf"),
  verbose = TRUE
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

- resample_rate:

  An *optional* sample rate (Hz) for the output data frame. If `NULL`
  (*default*) resamples to the existing `sample_rate`, which regularises
  any irregular samples without changing the rate.

- method:

  A character string specifying how new samples are filled. Default is
  *"none"*. Filling must be opted into explicitly (see *Details*):

  `"none"`

  :   Matches each new sample to the nearest original `time_channel`
      value without any interpolation, to within tolerance of half a
      sample-interval. New samples are returned as `NA`.

  `"locf"`

  :   (*"Last observation carried forward"*). Fills new and missing
      samples with the most recent valid non-`NA` value to the left, or
      the nearest valid value to the right for leading `NA`s. Safe for
      numeric, integer, and character columns.

  `"linear"`

  :   Fills new and missing samples via linear interpolation using
      [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html).
      Suitable for numeric columns only; non-numeric columns will fall
      back to `"locf"` behaviour.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
of class `"mnirs"`. Metadata are stored as attributes and can be
accessed with `attributes(data)`.

## Details

This function uses
[`replace_missing()`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)
(based on [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html))
to interpolate across new samples in the resampled data range.

### Sample rate and time channel

`time_channel` and `sample_rate` are retrieved automatically from `data`
of class *"mnirs"*, if not defined explicitly.

Otherwise, `sample_rate` will be estimated from the values in
`time_channel`. However, this may return unexpected values, and it is
safer to define `sample_rate` explicitly or retrieve it from *"mnirs"*
metadata.

### Default behaviour

When `resample_rate` is omitted, the output has the same `sample_rate`
as the input but with a regular, evenly-spaced `time_channel`. This is
useful for regularising data that contains missing or repeated samples
without changing the nominal rate.

### Column handling

Numeric columns are interpolated according to `method` (see
[`?replace_missing`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)).
Non-numeric columns (e.g. character event labels, integer lap numbers)
are always filled by last-observation-carried-forward, regardless of
`method`:

- For `method = "none"`, existing rows are matched to the nearest
  original values of `time_channel` without interpolation or filling,
  meaning newly created samples and any `NA`s in the original data are
  returned as `NA`.

- When down-sampling, numeric columns use time-weighted averaging.
  Non-numeric columns use the first valid value in each output bin.

## Examples

``` r
## read example data
data <- read_mnirs(
    file_path = example_mnirs("moxy_ramp"),
    nirs_channels = c(smo2 = "SmO2 Live"),
    time_channel = c(time = "hh:mm:ss"),
    verbose = TRUE
)
#> ! Estimated `sample_rate` = 2 Hz.
#> ℹ Define `sample_rate` explicitly to override.
#> Warning: ! Duplicate or irregular `time_channel` samples detected.
#> ℹ `time` = 211.59 and 1183.6.
#> ℹ Re-sample with `mnirs::resample_mnirs()`.

## note warning about irregular sampling
data
#> # A tibble: 2,202 × 2
#>     time  smo2
#>    <dbl> <dbl>
#>  1 0        54
#>  2 0.560    54
#>  3 1.11     54
#>  4 1.66     54
#>  5 2.21     54
#>  6 2.76     54
#>  7 3.31     57
#>  8 3.86     57
#>  9 4.41     57
#> 10 4.96     57
#> # ℹ 2,192 more rows

data_resampled <- resample_mnirs(
    data,               ## blank channels will be retrieved from metadata
    resample_rate = 2,  ## blank by default will resample to `sample_rate`
    method = "linear",  ## linear interpolation across resampled indices
    verbose = TRUE
)
#> ℹ Output is resampled at 2 Hz.

## note the altered `time` values resolving the above warning
data_resampled
#> # A tibble: 2,419 × 2
#>     time  smo2
#>    <dbl> <dbl>
#>  1   0    54  
#>  2   0.5  54  
#>  3   1    54  
#>  4   1.5  54  
#>  5   2    54  
#>  6   2.5  54  
#>  7   3    55.3
#>  8   3.5  57  
#>  9   4    57  
#> 10   4.5  57  
#> # ℹ 2,409 more rows
```
