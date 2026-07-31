# Apply a Butterworth digital filter

Apply a Butterworth digital filter to vector data with
[`signal::butter()`](https://rdrr.io/pkg/signal/man/butter.html) and
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)
which handles 'edges' better at the start and end of the data.

## Usage

``` r
filter_butterworth(
  x,
  order = 2L,
  W,
  type = c("low", "high", "stop", "pass"),
  edges = c("rev", "rep1", "none"),
  na.rm = FALSE,
  ...
)

filter_butter(
  x,
  order = 2L,
  W,
  type = c("low", "high", "stop", "pass"),
  edges = c("rev", "rep1", "none"),
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  A numeric vector.

- order:

  An integer defining the filter order (*default* `order = 2`).

- W:

  A one- or two-element numeric vector within `[0, 1]` defining the
  filter cutoff frequency(ies) as a fraction of the Nyquist frequency
  (see *Details*).

- type:

  A character string indicating the digital filter type (see *Details*).

  `"low"`

  :   For a *low-pass* filter (the *default*).

  `"high"`

  :   For a *high-pass* filter.

  `"stop"`

  :   For a *stop-band* (band-reject) filter.

  `"pass"`

  :   For a *pass-band* filter.

- edges:

  A character string indicating edge detection padding for `x`.

  `"rev"`

  :   Will pad `x` with the preceding 5% data in reverse sequence (*the
      default*).

  `"rep1"`

  :   Will pad `x` by repeating the last preceding value.

  `"none"`

  :   Will return the unpadded
      [`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)
      output.

- na.rm:

  Logical; default is `FALSE`, propagates any `NA`s to the returned
  vector. If `TRUE`, ignores `NA`s and processes available valid samples
  within the local window. May return errors or warnings. (see
  *Details*).

- ...:

  Additional arguments passed to the underlying method function. See
  *Details*.

## Value

A numeric vector the same length as `x`.

## Details

Applies a centred (two-pass symmetrical) Butterworth digital filter from
[`signal::butter()`](https://rdrr.io/pkg/signal/man/butter.html) and
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html).

Filter `type` defines how the desired signal frequencies are either
passed or rejected from the output signal. *Low-pass* and *high-pass*
filters allow only frequencies *lower* or *higher* than the cutoff
frequency `W` to be passed through as the output signal, respectively.
*Stop-band* defines a critical range of frequencies which are rejected
from the output signal. *Pass-band* defines a critical range of
frequencies which are passed through as the output signal.

The filter order (number of passes) is defined by `order`, typically in
the range `order = [1, 10]`. Higher filter order tends to capture more
rapid changes in amplitude, but also causes more distortion around those
change points in the signal. General advice is to use the lowest filter
order which sufficiently captures the desired rapid responses in the
data.

The critical (cutoff) frequency is defined by `W`, a numeric value for
*low-pass* and *high-pass* filters, or a two-element vector
`c(low, high)` defining the lower and upper bands for *stop-band* and
*pass-band* filters. `W` represents the desired fractional cutoff
frequency in the range `W = [0, 1]`, where `1` is the Nyquist frequency,
i.e., half the sample rate of the data in Hz.

Missing values (`NA`) in `x` will cause an error unless `na.rm = TRUE`.
Then `NA`s will be ignored and passed through to the returned vector.

## See also

[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html),
[`signal::butter()`](https://rdrr.io/pkg/signal/man/butter.html)

## Examples

``` r
set.seed(13)
sin <- sin(2 * pi * 1:150 / 50) * 20 + 40
noise <- rnorm(150, mean = 0, sd = 6)
noisy_sin <- sin + noise
without_edge_detection <- filter_butterworth(
    x = noisy_sin,
    order = 2,
    W = 0.1,
    edges = "none"
)
with_edge_detection <- filter_butterworth(
    x = noisy_sin,
    order = 2,
    W = 0.1,
    edges = "rep1"
)

ggplot2::ggplot(data.frame(), ggplot2::aes(x = seq_along(noise))) +
    theme_mnirs() +
    scale_colour_mnirs(name = NULL) +
    ggplot2::geom_line(ggplot2::aes(y = noisy_sin)) +
    ggplot2::geom_line(
        ggplot2::aes(y = without_edge_detection, colour = "without_edge_detection")
    ) +
    ggplot2::geom_line(
        ggplot2::aes(y = with_edge_detection, colour = "with_edge_detection")
    )
```
