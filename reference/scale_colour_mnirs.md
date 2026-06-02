# Scales for custom *mnirs* palette

Scales for custom *mnirs* palette

## Usage

``` r
scale_colour_mnirs(..., aesthetics = "colour")

scale_color_mnirs(..., aesthetics = "colour")

scale_fill_mnirs(..., aesthetics = "fill")
```

## Arguments

- ...:

  Arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

- aesthetics:

  A character vector with aesthetic(s) passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).
  *Default* is `"colour"`.

## Value

A [ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) scale
object.

## See also

[`theme_mnirs()`](https://jemarnold.github.io/mnirs/reference/theme_mnirs.md),
[`palette_mnirs()`](https://jemarnold.github.io/mnirs/reference/palette_mnirs.md)

## Examples

``` r
## plot example data
data <- read_mnirs(
    file_path = example_mnirs("moxy_ramp"),
    nirs_channels = c(smo2_left = "SmO2 Live",
                      smo2_right = "SmO2 Live(2)"),
    time_channel = c(time = "hh:mm:ss"),
    verbose = FALSE
)

ggplot2::ggplot(data, ggplot2::aes(x = time)) +
    theme_mnirs() +
    scale_colour_mnirs() +
    ggplot2::geom_line(ggplot2::aes(y = smo2_left, colour = "smo2_left")) +
    ggplot2::geom_line(ggplot2::aes(y = smo2_right, colour = "smo2_right"))
```
