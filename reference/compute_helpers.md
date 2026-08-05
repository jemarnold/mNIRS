# Computes rolling local values

`compute_window_bounds()`: Compute the start and end indices of rolling
windows along a time variable `t`.

`window_sums()`: Windowed sums by cumulative-sum differencing.

`window_min_obs()`: Minimum number of samples spanned by a complete
window.

`compute_local_mean()`: Compute rolling means from window bounds.

`compute_local_fun()`: Compute a rolling function along `x` from a list
of rolling sample windows.

`median_no_na()`: Fast median for numeric vectors. Strips `NA`s and
replicates `median.default` arithmetic without S3 dispatch.

`compute_col_medians()`: Column medians of an `NA`-padded numeric matrix
via a single radix sort. `NA`s sort last per column; medians indexed
from per-column valid counts. Matches `median(w, na.rm = TRUE)`.

`compute_outliers()`: Computes a vector of local medians and logicals
indicating outliers of `x` within rolling windows defined by `width` or
`span`.

`compute_valid_neighbours()`: Compute a list of rolling window indices
along `x` to either side of `NA`s.

## Usage

``` r
compute_window_bounds(
  t,
  idx = seq_along(t),
  width = NULL,
  span = NULL,
  align = c("centre", "left", "right"),
  env = rlang::caller_env()
)

window_sums(v, bounds)

window_min_obs(width, span, t, min_n = 1L, env = rlang::caller_env())

compute_local_mean(x, bounds, na.rm = FALSE, min_obs = 1L)

compute_local_fun(x, window_idx, fn, ...)

median_no_na(w)

compute_col_medians(m)

compute_outliers(
  x,
  t,
  outlier_cutoff,
  width = NULL,
  span = NULL,
  env = rlang::caller_env()
)

compute_valid_neighbours(
  x,
  t = seq_along(x),
  width = NULL,
  span = NULL,
  verbose = TRUE,
  env = rlang::caller_env()
)
```

## Arguments

- t:

  An *optional* numeric vector of the predictor variable (e.g. time).
  Default is `seq_along(x)`.

- idx:

  A numeric vector of indices of `t` at which to calculate local
  windows. All indices of `t` by *default*, or can be used to only
  calculate for known indices, such as invalid values of `x`.

- width:

  An integer defining the local window in number of samples around `idx`
  in which to perform the operation, according to `align`.

- span:

  A numeric value defining the local window time span around `idx` in
  which to perform the operation, according to `align`. In units of
  `time_channel` or `t`.

- align:

  Window alignment as *"centre"/"center"* (the *default*), *"left"*, or
  *"right"*. Where *"left"* is *forward looking*, and *"right"* is
  *backward looking* from the current sample.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.

- v:

  A numeric vector to sum within windows. Callers should centre `v`
  first to contain floating-point cancellation error.

- bounds:

  A [`list()`](https://rdrr.io/r/base/list.html) of `start` and `end`
  window index vectors from `compute_window_bounds()`.

- min_n:

  A lower bound on the returned number of samples.

- x:

  A numeric vector of the response variable.

- min_obs:

  The minimum number of samples a window must span to return a value.
  Shorter (partial) windows return `NA`.

- window_idx:

  A list the same or shorter length as `x` with numeric vectors for the
  sample indices of local rolling windows.

- fn:

  A function to pass through for local rolling calculation.

- ...:

  Additional arguments.

- m:

  A numeric matrix with one column per rolling window, padded with `NA`
  where windows extend beyond the data.

- outlier_cutoff:

  A numeric value for the local outlier threshold, as the number of
  standard deviations from the local median.

  - Default `NULL` will not replace outliers.

  - Lower values are more sensitive and flag more outliers; higher
    values are more conservative.

  - `outlier_cutoff = 3` Pearson's 3 sigma edit rule.
    `outlier_cutoff = 2` approximates a Tukey-style 1.5\*IQR rule.
    `outlier_cutoff = 0` Tukey's median filter.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

## Value

`compute_window_bounds()`: A
[`list()`](https://rdrr.io/r/base/list.html) with `start` and `end`
integer vectors the same length as `idx`, giving the inclusive window
bounds at each index.

`window_sums()`: A numeric vector the same length as `bounds$start`.

`window_min_obs()`: An integer value.

`compute_local_mean()`: A numeric vector the same length as
`bounds$start`.

`compute_local_fun()`: A numeric vector the same length as `x`.

`median_no_na()`: A numeric value.

`compute_col_medians()`: A numeric vector of length `ncol(m)`.

`compute_outliers()`: A [`list()`](https://rdrr.io/r/base/list.html)
with vectors the same length as `x` for with numeric local medians and
logical identifying where `is_outlier`.

`compute_valid_neighbours()`: A list the same length as the `NA` values
in `x` with numeric vectors of sample indices of length `width` samples
or `span` units of time `t` for valid values neighbouring split to
either side of the invalid `NA`s.

## Details

The local rolling window can be specified by either `width` as the
number of samples, or `span` as the time span in units of `t`.
Specifying `width` is often faster than `span`.

`align` defaults to *"centre"* the local window around `idx` between
`[idx - floor((width-1)/2),` `idx + floor(width/2)]` when `width` is
specified. Even `width` values will bias `align` to *"left"*, with the
unequal sample forward of `idx`, effectively returning `NA` at the last
sample index. When `span` is specified, the local window is between
`[t - span/2, t + span/2]`.

`window_min_obs()` converts `span` to a sample count via the estimated
sample rate, less two samples to buffer irregular `t` at the start and
end of each window.

`compute_local_mean()` computes all window means in O(n) via
`window_sums()`. Values are centred first so the cumulative-sum
differencing error stays around `eps * sqrt(n) * sd`, far below
measurement resolution.
