# Rescale data range

Expand or reduce the range (min and max values) of data channels to a
new amplitude/dynamic range, e.g. rescale the range of NIRS data to
`c(0, 100)`.

## Usage

``` r
rescale_mnirs(
  data,
  nirs_channels = NULL,
  group_channels = c("ensemble", "distinct"),
  range,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata, a list of data frames, or a grouped data frame (see
  *Details*).

- nirs_channels:

  A character vector giving the names of mNIRS columns to operate on.
  Must match column names in `data` exactly.

  - If `NULL` (default), the `nirs_channels` metadata attribute of
    `data` is used.

- group_channels:

  Either a character string or a
  [`list()`](https://rdrr.io/r/base/list.html) of channel-name vectors
  specifying how to group `nirs_channels` (see *Details*).

  `"ensemble"`

  :   The *default*. Operate on all channels together, preserving the
      relative scaling between channels.

  `"distinct"`

  :   Operate on each channel independently, losing the relative scaling
      between channels.

  `list(c("A", "B"), c("C", "D"))`

  :   Operate on channels `A` & `B` in one group, and `C` & `D` in
      another group. Groups can be named (e.g.
      `list(smo2 = c("A", "B"))`). Each group must be non-empty and
      resulting group names must be unique.

- range:

  A numeric vector in the form `c(min, max)`, indicating the range of
  output values to which `nirs_channels` will be rescaled.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
of class *"mnirs"* with metadata available with
[`attributes()`](https://rdrr.io/r/base/attributes.html). For list or
grouped data frame input, returns a named list of *"mnirs"* tibbles, one
per interval.

## Details

`group_channels` controls how data channels are grouped to preserve
absolute or relative scaling.

- `group_channels = "ensemble"` (the *default*) rescales all
  `nirs_channels` to a common range, preserving relative scaling between
  channels.

- `group_channels = "distinct"` rescales each channel independently,
  losing relative scaling between channels.

- A [`list()`](https://rdrr.io/r/base/list.html) of channel-name vectors
  (e.g. `list(c("A", "B"), c("C", "D"))`) rescales channels `A` & `B`
  together and `C` & `D` together, preserving relative scaling within,
  but not between groups. `nirs_channels` omitted from the list are
  rescaled independently.

- Channel groups can be named (e.g. `list(smo2 = c("A", "B"))`) and
  names used as keys for per-group `range` argument.

- Channels (columns) in `data` not in `nirs_channels` are passed through
  without processing to the output data frame.

`nirs_channels` can be retrieved automatically from `data` of class
*"mnirs"* which has been processed with `{mnirs}`, if not defined
explicitly.

## Data input formats

*mnirs* processing functions accept `data` in multiple formats:

- A **single *"mnirs"* data frame** is processed and returned directly.

- A **list of *"mnirs"* data frames**: each interval is processed
  separately and returned as a named list.

- A **grouped *"mnirs"* data frame**, e.g. with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html):
  the data frame is split by grouping levels and each group is processed
  as a separate interval, returned as a named list.

## Per-channel arguments

Arguments apply globally to all `nirs_channels` by default. Relevant
arguments can instead be supplied uniquely per-channel as a named
[`list()`](https://rdrr.io/r/base/list.html), with names matching either
`nirs_channels` or list names in `group_channels`, e.g.:

    shift_mnirs(
        data,
        nirs_channels = c(A, B, C),
        group_channels = list(smo2 = c(A, B), hhb = C),
        to = list(100, C = 0),
        width = list(smo2 = 3),
        span = list(hhb = 5),
        position = "first"
    )

- A non-list value applies to every channel (the *default* behaviour).

- A [`list()`](https://rdrr.io/r/base/list.html) named by
  `nirs_channels` or `group_channels` applies per-channel / per-group
  values.

- A single unnamed value in the list will be applied to unlisted
  channels (e.g. `span = list(3, hhb = 5)` gives `hhb` 5 and every other
  channel 3). If no unnamed fallback value in the list, channels not
  named in the list will be returned un-processed (e.g.
  `span = list(hhb = 5)` will only process `hhb`).

- [`list()`](https://rdrr.io/r/base/list.html) names not matching
  `nirs_channels` or `group_channels` are warned about and ignored.

## Examples

``` r
## read example data
data <- read_mnirs(
    file_path = example_mnirs("moxy_ramp"),
    nirs_channels = c(smo2_left = "SmO2 Live",
                      smo2_right = "SmO2 Live(2)"),
    time_channel = c(time = "hh:mm:ss"),
    verbose = FALSE
) |>
    rescale_mnirs(        ## un-grouped nirs channels to rescale separately
        nirs_channels = c(smo2_left, smo2_right),
        group_channels = "distinct",
        range = c(0, 100)  ## rescale to a 0-100% functional exercise range
    )

data
#> # A tibble: 2,202 × 3
#>     time smo2_left smo2_right
#>    <dbl>     <dbl>      <dbl>
#>  1 0            54       78.2
#>  2 0.560        54       78.2
#>  3 1.11         54       75.6
#>  4 1.66         54       75.6
#>  5 2.21         54       75.6
#>  6 2.76         54       75.6
#>  7 3.31         57       76.9
#>  8 3.86         57       76.9
#>  9 4.41         57       76.9
#> 10 4.96         57       76.9
#> # ℹ 2,192 more rows

# \donttest{
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot(data, time_labels = TRUE) +
            ggplot2::geom_hline(yintercept = c(0, 100), linetype = "dotted")
    }

# }
```
