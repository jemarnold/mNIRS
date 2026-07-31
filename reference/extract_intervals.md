# Extract intervals from *mnirs* data

Extract intervals from *"mnirs"* time series data, specifying interval
start and end boundaries by time value, event label, lap number, or
sample index.

## Usage

``` r
extract_intervals(
  data,
  nirs_channels = NULL,
  time_channel = NULL,
  event_channel = NULL,
  sample_rate = NULL,
  group_intervals = c("distinct", "ensemble"),
  group_channels = NULL,
  start = NULL,
  end = NULL,
  span = list(c(-60, 60)),
  zero_time = FALSE,
  verbose = TRUE,
  event_groups = deprecated()
)
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata, or a [`list()`](https://rdrr.io/r/base/list.html) of such
  data frames.

- nirs_channels:

  A character vector of mNIRS channel names to operate on. Names must
  match column names in `data` exactly.

  - If `NULL` (default), channels are retrieved from *"mnirs"* metadata.

  - **\[deprecated\]** Passing a
    [`list()`](https://rdrr.io/r/base/list.html) for per-group channel
    selection is deprecated; use `group_channels`.

- time_channel:

  A character string naming the time or sample column. Must match a
  column name in `data` exactly.

  - If `NULL` (default), the `time_channel` metadata attribute of `data`
    is used.

- event_channel:

  An *optional* character string giving the name of an event/lap column.
  The column may contain character event labels or integer lap numbers.

  - Required when using
    [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
    or
    [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
    for `start` or `end`.

  - Retrieved from metadata if not defined explicitly.

- sample_rate:

  An *optional* numeric sample rate (Hz) used to bin time values for
  ensemble-averaging. If `NULL`, will be estimated from `time_channel`
  (see *Details*).

- group_intervals:

  Either a character string or a non-empty
  [`list()`](https://rdrr.io/r/base/list.html) of non-empty
  integer-valued numeric vectors specifying how to group intervals (see
  *Details*). Custom indices must be between `1` and the number of
  detected intervals.

  `"distinct"`

  :   The default. Extract each interval as an independent data frame.

  `"ensemble"`

  :   Ensemble-average each specified `nirs_channel` across all detected
      intervals, returning a single data frame.

  `list(c(1, 2), c(3, 4))`

  :   Ensemble-average each specified `nirs_channel` within each group
      and return one data frame per group.

- group_channels:

  A character vector or a [`list()`](https://rdrr.io/r/base/list.html)
  of character vectors selecting which `nirs_channels` are
  ensemble-averaged within each interval group (see *Details*).

  - If `NULL` (default), all `nirs_channels` are used for every group.

  - Only relevant when `group_intervals` contains *"ensemble"*-averaged
    intervals; with `group_intervals = "distinct"` no channel processing
    occurs.

- start:

  Specifies where intervals begin. Either raw values – numeric for time
  values, character for event labels, explicit integer (e.g. `2L`) for
  lap numbers – or created with
  [`by_time()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  or
  [`by_sample()`](https://jemarnold.github.io/mnirs/reference/by_time.md).

- end:

  Specifies where intervals end. Either raw values – numeric for time
  values, character for event labels, explicit integer (e.g. `2L`) for
  lap numbers – or created with
  [`by_time()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md),
  or
  [`by_sample()`](https://jemarnold.github.io/mnirs/reference/by_time.md).

- span:

  A one- or two-element numeric vector `c(before, after)` in units of
  `time_channel`, or a [`list()`](https://rdrr.io/r/base/list.html) of
  such vectors. (*default* `span = c(-60, 60)`. Applied additively to
  interval boundaries:

  - When both `start` and `end` are specified: `span[1]` shifts start
    times, `span[2]` shifts end times.

  - When only `start` or only `end` is specified: both `span[1]` and
    `span[2]` apply as a window around the event).

  - A single *positive* value is recycled to shift the end times (e.g.
    `span = 60` -\> `c(0, 60)`).

  - A single *negative* value is recycled to shift the start times (e.g.
    `span = -60` -\> `c(-60, 0)`).

- zero_time:

  Logical. Default is `FALSE`. If `TRUE`, re-calculates numeric
  `time_channel` values to start from zero within each interval data
  frame.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- event_groups:

  **\[deprecated\]** Renamed to `group_intervals` for naming consistency
  across the package.

## Value

A named [`list()`](https://rdrr.io/r/base/list.html) of
[tibbles](https://tibble.tidyverse.org/reference/tibble-package.html) of
class *"mnirs"*, each with metadata available via
[`attributes()`](https://rdrr.io/r/base/attributes.html). When `data` is
a list of data frames, results are flattened into a single-layer list
with interval names indicating nested number of data frame and interval
as `interval_<df>.<interval>`. Other interval names (e.g. `"ensemble"`
or custom group names) are suffixed as `<name>_<df>`.

## Details

### Interval specification

Interval `start` and `end` boundaries are specified using helper
functions, or by passing raw values directly:

- [`by_time()`](https://jemarnold.github.io/mnirs/reference/by_time.md):

  Time values in units of `time_channel`.

- [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md):

  Strings to match in `event_channel`. All matching occurrences are
  returned.

- [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md):

  Lap numbers to match in `event_channel`. Resolves to the first sample
  of each lap for `start`, and the last lap sample for `end`

- [`by_sample()`](https://jemarnold.github.io/mnirs/reference/by_time.md):

  Integer sample indices (row numbers).

Raw values supplied to `start`/`end` are auto-coerced:

- Numeric -\>
  [`by_time()`](https://jemarnold.github.io/mnirs/reference/by_time.md)

- Character -\>
  [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md),

- Explicit integer (e.g. `2L`) -\>
  [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md).

- Use
  [`by_sample()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
  explicitly for sample indices.

`start` and `end` can use different specification types (e.g., start by
label, end by time). When lengths differ, the shorter is recycled.

### Time span window

`span` additively expands the time span window around interval
boundaries.

- A two-value vector expands the `start` and `end`, respectively:
  `span = c(-60, 60)` expands the `start` earlier by `60`, and the `end`
  later by `60`. For example,
  `start = by_time(30), end = by_time(60), span = c(-5, 10)` returns an
  interval of `[25, 70]`.

- A single numeric value is recycled according to the sign: `span = -60`
  becomes `c(-60, 0)` to expand the `start` earlier. `span = 60` becomes
  `c(0, 60)` to expand the `end` later.

- If only `start` is specified alone, both span values expand the single
  boundary window: `start = by_time(30), span = c(-5, 60)` returns
  `[25, 90]`.

### Per-group channel selection with group_channels

When `group_intervals = "ensemble"` or a list of numeric grouped
intervals, `group_channels` can be specified as a list of column names
to override which channels are ensemble-averaged within each group. For
example, to exclude a channel from one interval:

    group_channels = list(
      c(A, B, C),
      c(A, C) ## channel "B" data are excluded from the second interval
    )

If all grouped intervals include all `nirs_channels`, `group_channels`
can be left as `NULL` (the default) and all channels are
ensemble-averaged within every group.

### Grouping intervals

`group_intervals` controls whether extracted intervals are returned as
distinct data frames or ensemble-averaged.

- `"distinct"`:

  The default. Extract each interval and return a list of independent
  data frames.

- `"ensemble"`:

  Ensemble-average each specified `nirs_channel` across all detected
  intervals and return a one-item list with a single data frame.

- `list(c(1, 2), c(3, 4))`:

  Ensemble-average each specified `nirs_channel` within each group and
  return a list with one data frame for each group. Any intervals
  detected but not specified in `group_intervals` are returned as
  distinct.

`group_intervals` lists can be named (e.g.
`list(low = c(1, 2), high = c(3, 4))`) and will pass those names to the
returned list of data frames.

When `group_intervals` is a list of numeric interval numbers, list items
in `group_channels` and `span` are recycled to the number of groups. If
lists are only partially specified, the final item is recycled forward
as needed. Extra items are ignored.

## Examples

``` r
## read example data
data <- read_mnirs(
    example_mnirs("train.red"),
    nirs_channels = c(
        smo2_left = "SmO2 unfiltered",
        smo2_right = "SmO2 unfiltered"
    ),
    time_channel = c(time = "Timestamp (seconds passed)"),
    zero_time = TRUE,
    verbose = FALSE
) |>
    ## avoid issues ensemble-averaging irregular samples
    resample_mnirs(method = "linear", verbose = FALSE)

## ensemble-average across multiple intervals
interval_list <- extract_intervals(
    data,                         ## channels recycled to all intervals by default
    nirs_channels = c(smo2_left, smo2_right),
    group_intervals = "ensemble", ## ensemble-average across two intervals
    start = by_time(368, 1084),   ## manually identified interval start times
    span = c(-20, 90),            ## include the last 180-sec of each interval (recycled)
    zero_time = TRUE              ## re-calculate common time to start from `0`
)

interval_list[[1L]]
#> # A tibble: 1,101 × 3
#>     time smo2_left smo2_right
#>    <dbl>     <dbl>      <dbl>
#>  1 -20        56.3       59.2
#>  2 -19.9      56.1       59.2
#>  3 -19.8      56.1       59.2
#>  4 -19.7      56.2       58.9
#>  5 -19.6      56.4       58.9
#>  6 -19.5      56.5       58.9
#>  7 -19.4      56.9       58.7
#>  8 -19.3      57.0       58.7
#>  9 -19.2      56.9       59.0
#> 10 -19.1      56.7       58.8
#> # ℹ 1,091 more rows

# \donttest{
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot(interval_list, time_labels = TRUE) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
  }

# }
```
