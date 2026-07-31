# Apply a moving average filter

Apply a simple moving average smoothing filter to vector data.
`filter_ma()` is an alias of `filter_moving_average()`.

## Usage

``` r
filter_moving_average(
  x,
  t = seq_along(x),
  width = NULL,
  span = NULL,
  partial = FALSE,
  na.rm = FALSE,
  verbose = TRUE,
  ...
)

filter_ma(
  x,
  t = seq_along(x),
  width = NULL,
  span = NULL,
  partial = FALSE,
  na.rm = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A numeric vector of the response variable.

- t:

  An *optional* numeric vector of the predictor variable (e.g. time).
  Default is `seq_along(x)`.

- width:

  An integer defining the local window in number of samples centred on
  `idx`, between `[idx - floor(width/2), idx + floor(width/2)]`.

- span:

  A numeric value defining the local window time span around `idx` in
  units of `time_channel` or `t`, between `[t - span/2, t + span/2]`.

- partial:

  Logical; default is `FALSE`, only returns values where a full window
  of valid (non-`NA`) samples are available. If `TRUE`, ignores `NA` and
  processes available valid samples (see *Details*).

- na.rm:

  Logical; default is `FALSE`, propagates any `NA`s to the returned
  vector. If `TRUE`, ignores `NA`s and processes available valid samples
  within the local window. May return errors or warnings. (see
  *Details*).

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- ...:

  Additional arguments.

## Value

A numeric vector the same length as `x`.

## Details

### Rolling window

Applies a centred (symmetrical) moving average filter in a local window,
defined by either `width` as the number of samples around `idx` between
`[idx - floor(width/2), idx + floor(width/2)]`. Or by `span` as the
timespan in units of `time_channel` between `[t - span/2, t + span/2]`.

### Partial windows

The default `partial = FALSE` requires a complete number of samples
specified by `width` or `span` (estimated from the sample rate of `t`
when `span` is used). `NA` is returned if fewer samples are present in
the local window.

Setting `partial = TRUE` allows computation with only a single valid
sample, such as at edge conditions. But these values will be more
sensitive to noise and should be used with caution.

### Missing values

`na.rm` controls whether missing values (`NA`s) within each local window
are either propagated to the returned vector when `na.rm = FALSE` (the
default), or ignored before processing if `na.rm = TRUE`.

## Examples

``` r
x <- c(1, 3, 2, 5, 4, 6, 5, 7)
t <- c(0, 1, 2, 4, 5, 6, 7, 10)  ## irregular time with gaps

## width: centred window of 3 samples
filter_moving_average(x, width = 3)
#> [1]       NA 2.000000 3.333333 3.666667 5.000000 5.000000 6.000000       NA

## partial = TRUE fills edge values with a narrower window
filter_moving_average(x, width = 3, partial = TRUE)
#> [1] 2.000000 2.000000 3.333333 3.666667 5.000000 5.000000 6.000000 6.000000

## span: centred window of 2 time-units (accounts for irregular sampling)
filter_moving_average(x, t, span = 2)
#> [1] 2.0 2.0 2.5 4.5 5.0 5.0 5.5 7.0

## na.rm = FALSE (default): any NA in the window propagates to the result
x_na <- c(1, NA, 3, 4, 5, NA, 7, 8)
filter_moving_average(x_na, width = 3, na.rm = FALSE)
#> Warning: ! `x` contains internal NA's.
#> ℹ Set `na.rm = TRUE` to ignore NA's.
#> [1] NA NA NA  4 NA NA NA NA

## na.rm = TRUE: skip NAs and return the local mean of local valid values
filter_moving_average(x_na, width = 3, partial = TRUE, na.rm = TRUE)
#> [1] 1.0 2.0 3.5 4.0 4.5 6.0 7.5 7.5
```
