# Plot *mnirs* objects

Create a base plot for data frames or lists of data frames with class
*"mnirs"*.

## Usage

``` r
# S3 method for class 'mnirs'
plot(x, points = FALSE, time_labels = FALSE, na.omit = FALSE, ...)
```

## Arguments

- x:

  Data frame or list of data frames of class *"mnirs"* (e.g. from
  [`extract_intervals()`](https://jemarnold.github.io/mnirs/reference/extract_intervals.md)).
  List input produces a faceted plot with one panel per element.

- points:

  Logical. Default is `FALSE`. If `TRUE` displays
  `ggplot2::geom_points()`. Otherwise displays `ggplot2::geom_lines()`.

- time_labels:

  Logical. Default is `FALSE`. If `TRUE` displays x-axis time values
  formatted as *"hh:mm:ss"* using
  [`format_hmmss()`](https://jemarnold.github.io/mnirs/reference/format_hmmss.md).
  Otherwise, x-axis values are displayed as numeric.

- na.omit:

  Logical. Default is `FALSE`. If `TRUE` omits missing (`NA`) and
  non-finite `c(Inf, -Inf, NaN)` from display.

- ...:

  Additional arguments.

## Value

A [ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Details

When `x` is a named list of *"mnirs"* data frames, elements are bound
into a single data frame and displayed as faceted panels via
[`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

Accepts some arguments in `...`, such as `nrow`, `ncol`, and `scales`
passed to
[`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
`n.breaks` overrides the default number of y-axis breaks. `breaks`
overrides the x-axis breaks directly.

## Examples

``` r
data <- read_mnirs(
    example_mnirs("train.red"),
    nirs_channels = c(smo2 = "SmO2"),
    time_channel = c(time = "Timestamp (seconds passed)"),
    verbose = FALSE
)

## plot time labels as "hh:mm:ss"
plot(data, time_labels = TRUE)


data_list <- extract_intervals(
    data,
    start = by_time(2452, 3168),
    span = c(-60, 120),
    verbose = FALSE
)

## plot a list of mnirs data frames as faceted panels
plot(data_list, time_labels = TRUE)
```
