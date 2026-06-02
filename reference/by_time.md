# Specify interval boundaries by time, label, lap, or sample

Helper functions to define interval start or end boundaries for
[`extract_intervals()`](https://jemarnold.github.io/mnirs/reference/extract_intervals.md).

## Usage

``` r
by_time(...)

by_label(..., ignore_case = FALSE, fixed = FALSE)

by_lap(...)

by_sample(...)
```

## Arguments

- ...:

  Specify start or end boundaries.

  `by_time(...)`

  :   Numeric time values in units of `time_channel`.

  `by_label(...)`

  :   Character strings to match in `event_channel`. Matched as regular
      expressions by default; see `ignore_case` and `fixed`. All
      matching occurrences are returned.

  `by_lap(...)`

  :   Integer lap numbers to match in `event_channel`. For `start`,
      resolves to the first sample of each lap. For `end`, resolves to
      the last sample.

  `by_sample(...)`

  :   Integer sample indices (row numbers).

- ignore_case:

  For `by_label()`. If `TRUE`, match case-insensitive labels. Default
  `FALSE`.

- fixed:

  For `by_label()`. If `TRUE`, treat labels as fixed strings rather than
  regular expressions. Useful when labels contain regex metacharacters
  (`.`, `*`, `(`, etc.). Default `FALSE`.

## Value

An object of class `"mnirs_interval"` for use with the `start` and `end`
arguments of
[`extract_intervals()`](https://jemarnold.github.io/mnirs/reference/extract_intervals.md).

## Details

These helpers can be used explicitly for arguments `start`/`end`, or raw
values can be passed directly:

- Numeric -\> `by_time()`

- Character -\> `by_label()`,

- Explicit integer (e.g. `2L`) -\> `by_lap()`.

- Use `by_sample()` explicitly for sample indices.

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
    event_channel = c(lap = "Lap/Event"),
    zero_time = TRUE,
    verbose = FALSE
)

## start and end by time
extract_intervals(data, start = by_time(66), end = by_time(357))
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> $interval_1
#> # A tibble: 4,150 × 4
#>     time   lap smo2_left smo2_right
#>    <dbl> <int>     <dbl>      <dbl>
#>  1  6.04     1      67.6       70.7
#>  2  6.04     1      67.6       71.8
#>  3  6.24     1      67.6       71.4
#>  4  6.24     1      67.6       71.7
#>  5  6.30     1      67.6       71.6
#>  6  6.39     1      65.7       70.9
#>  7  6.48     1      67.4       72.0
#>  8  6.63     1      67.0       72.2
#>  9  6.69     1      67.6       70.8
#> 10  7.61     1      68.4       71.5
#> # ℹ 4,140 more rows
#> 

## start by lap
extract_intervals(data, start = by_lap(2, 4), span = 0)
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> $interval_1
#> # A tibble: 1 × 4
#>    time   lap smo2_left smo2_right
#>   <dbl> <int>     <dbl>      <dbl>
#> 1  65.9     2      68.4       71.6
#> 
#> $interval_2
#> # A tibble: 1 × 4
#>    time   lap smo2_left smo2_right
#>   <dbl> <int>     <dbl>      <dbl>
#> 1  788.     4      69.4       70.2
#> 

## introduce event_channel with "start" string
data$event <- NA_character_
data$event[1000] <- "start"
data <- create_mnirs_data(data, event_channel = "event")

## start by label, end by time
extract_intervals(data, start = by_label("start"), end = by_time(1500))
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> Warning: !  Interval 1 is partially outside data bounds.
#> ℹ Returning available data only.
#> $interval_1
#> # A tibble: 11,601 × 5
#>     time   lap smo2_left smo2_right event
#>    <dbl> <int>     <dbl>      <dbl> <chr>
#>  1  39.0     1      67.8       69.3 NA   
#>  2  39.0     1      67.8       69.1 NA   
#>  3  39.2     1      68.0       70.1 NA   
#>  4  39.2     1      68.0       69.8 NA   
#>  5  39.3     1      68.2       69.4 NA   
#>  6  39.4     1      68.2       69.7 NA   
#>  7  39.5     1      68.0       69.2 NA   
#>  8  39.6     1      68.5       69.2 NA   
#>  9  39.7     1      68.5       70.0 NA   
#> 10  39.9     1      68.5       69.8 NA   
#> # ℹ 11,591 more rows
#> 

## case-insensitive label match
extract_intervals(data, start = by_label("START", ignore_case = TRUE))
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> $interval_1
#> # A tibble: 1,211 × 5
#>     time   lap smo2_left smo2_right event
#>    <dbl> <int>     <dbl>      <dbl> <chr>
#>  1  39.0     1      67.8       69.3 NA   
#>  2  39.0     1      67.8       69.1 NA   
#>  3  39.2     1      68.0       70.1 NA   
#>  4  39.2     1      68.0       69.8 NA   
#>  5  39.3     1      68.2       69.4 NA   
#>  6  39.4     1      68.2       69.7 NA   
#>  7  39.5     1      68.0       69.2 NA   
#>  8  39.6     1      68.5       69.2 NA   
#>  9  39.7     1      68.5       70.0 NA   
#> 10  39.9     1      68.5       69.8 NA   
#> # ℹ 1,201 more rows
#> 

## literal-string label match (regex metacharacters treated as text)
data$event[1000] <- "lap.1"
data <- create_mnirs_data(data, event_channel = "event")
extract_intervals(data, start = by_label("lap.1", fixed = TRUE))
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> $interval_1
#> # A tibble: 1,211 × 5
#>     time   lap smo2_left smo2_right event
#>    <dbl> <int>     <dbl>      <dbl> <chr>
#>  1  39.0     1      67.8       69.3 NA   
#>  2  39.0     1      67.8       69.1 NA   
#>  3  39.2     1      68.0       70.1 NA   
#>  4  39.2     1      68.0       69.8 NA   
#>  5  39.3     1      68.2       69.4 NA   
#>  6  39.4     1      68.2       69.7 NA   
#>  7  39.5     1      68.0       69.2 NA   
#>  8  39.6     1      68.5       69.2 NA   
#>  9  39.7     1      68.5       70.0 NA   
#> 10  39.9     1      68.5       69.8 NA   
#> # ℹ 1,201 more rows
#> 

## multiple intervals by sample index
extract_intervals(data, start = by_sample(1000, 1500))
#> ℹ `nirs_channels` = c("smo2_left", "smo2_right") grouped together from
#>   metadata.
#> $interval_1
#> # A tibble: 1,211 × 5
#>     time   lap smo2_left smo2_right event
#>    <dbl> <int>     <dbl>      <dbl> <chr>
#>  1  39.0     1      67.8       69.3 NA   
#>  2  39.0     1      67.8       69.1 NA   
#>  3  39.2     1      68.0       70.1 NA   
#>  4  39.2     1      68.0       69.8 NA   
#>  5  39.3     1      68.2       69.4 NA   
#>  6  39.4     1      68.2       69.7 NA   
#>  7  39.5     1      68.0       69.2 NA   
#>  8  39.6     1      68.5       69.2 NA   
#>  9  39.7     1      68.5       70.0 NA   
#> 10  39.9     1      68.5       69.8 NA   
#> # ℹ 1,201 more rows
#> 
#> $interval_2
#> # A tibble: 1,211 × 5
#>     time   lap smo2_left smo2_right event
#>    <dbl> <int>     <dbl>      <dbl> <chr>
#>  1  88.5     2      58.7       65.2 NA   
#>  2  88.7     2      59.7       65.0 NA   
#>  3  88.7     2      59.7       65.6 NA   
#>  4  88.9     2      58.6       66.4 NA   
#>  5  88.9     2      58.6       66.4 NA   
#>  6  89.1     2      58.9       65.3 NA   
#>  7  89.1     2      58.9       64.8 NA   
#>  8  89.2     2      58.9       65.7 NA   
#>  9  89.3     2      58.9       65.0 NA   
#> 10  89.4     2      58.8       65.1 NA   
#> # ℹ 1,201 more rows
#> 
```
