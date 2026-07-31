# Read *mnirs* data from file

Import time-series data exported from common muscle NIRS (mNIRS) devices
and return a tibble of class `"mnirs"` with the selected signal channels
and metadata.

## Usage

``` r
read_mnirs(
  file_path,
  nirs_channels = NULL,
  time_channel = NULL,
  event_channel = NULL,
  sample_rate = NULL,
  add_timestamp = FALSE,
  zero_time = FALSE,
  keep_all = FALSE,
  verbose = TRUE
)
```

## Arguments

- file_path:

  Path of the data file to import. Supported file extensions include
  `".xlsx"`, `".xls"`, `".csv"`, and `".txt"`.

- nirs_channels:

  A character vector of one or more column names containing mNIRS
  signals to import. Names must match the file header exactly.

  - If `NULL` (default), `read_mnirs()` attempts to detect the device
    from the file contents and use a known `nirs_channel` name.

  - A *named* character vector can be used to rename columns on import,
    in the form `c(renamed = "original_name")`.

- time_channel:

  A character string giving the name of the time (or sample) column to
  import. The name must match the file header exactly.

  - If `NULL` (default), `read_mnirs()` attempts to identify a time-like
    column automatically (by known device defaults and/or time-formatted
    values).

  - A *named* character vector can be used to rename the column on
    import, in the form `c(time = "original_name")`.

- event_channel:

  An *optional* character string giving the name of an event/lap column
  to import. Names must match the file header exactly. A named character
  vector can be used to rename the column on import in the form
  `c(event = "original_name")`.

- sample_rate:

  An *optional* numeric sample rate in Hz. If left blank (`NULL`), the
  sample rate is estimated from `time_channel` (see *Details*).

- add_timestamp:

  A logical. Default is `FALSE`. If `TRUE` and if the source data
  contain an absolute date-time (POSIXct) time value, will add a
  `"timestamp"` column in addition to the specified `time_channel` as a
  numeric time column.

- zero_time:

  Logical. Default is `FALSE`. If `TRUE`, re-calculates numeric
  `time_channel` values to start from zero.

- keep_all:

  Logical. Default is `FALSE`. Will keep only the channels explicitly
  specified in `nirs_channels`, `time_channel`, and `event_channel`. If
  `TRUE` will keep all columns found in the file data table.

  - If no `nirs_channels` are specified and the file format is
    recognised, all columns in the file data table will be returned, as
    an exploratory option.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
of class `"mnirs"`. Metadata are stored as attributes and can be
accessed with `attributes(data)`.

## Details

### Header detection

`read_mnirs()` searches the file for a header row containing the
requested channel names. The header row does not need to be the first
row in the file.

- If duplicate column names exist, columns are matched in the order they
  appear and renamed with unique strings.

- Columns without a header name in the source file will be renamed to
  `col_*`, where `*` is the numeric column number in which they appear
  in the file (e.g. `col_6`). This applies to *Artinis Oxysoft* event
  label columns, which do not have a column header and must be
  identified manually.

### Renaming channels

A named character vector can be specified to rename `nirs_channels`,
`time_channel`, and `event_channel`, in the form
`c(renamed = "original_name")`. The `"original_name"` must match the
contents of the file data table header row exactly.

### Time parsing

`time_channel` will be converted to numeric for analysis.

- If `time_channel` is a date-time (POSIXct) format, it will be
  converted to numeric and re-based to start from 0, regardless of
  `zero_time`.

- Some devices export a sample index rather than time values. In those
  cases, if an export `sample_rate` is detected in the file metadata
  (e.g. *Artinis Oxysoft* exports), `read_mnirs()` will create or
  overwrite a `"time"` column in seconds derived from the sample index
  and the detected `sample_rate`.

### Sample rate

If `sample_rate` is not specified, it is estimated from differences in
`time_channel`. If `time_channel` is actually a sample index, as
described above, this may erroneously be estimated at 1 Hz.
`sample_rate` should be specified explicitly in this case.

### Data cleaning

Entirely empty rows and columns are removed. Invalid values (e.g.
`c(NaN, Inf)`) are standardized to `NA`. A warning will be displayed
when irregular sampling is detected (e.g. non-monotonic, repeated, or
unequal time values), if `verbose = TRUE`.

## Examples

``` r
read_mnirs(
    file_path = example_mnirs("moxy_ramp"), ## call an example data file
    nirs_channels = c(
        smo2_left = "SmO2 Live",            ## identify and rename channels
        smo2_right = "SmO2 Live(2)"
    ),
    time_channel = c(time = "hh:mm:ss"),    ## date-time format will be converted to numeric
    event_channel = NULL,                   ## leave blank if unused
    sample_rate = NULL,                     ## if blank, will be estimated from time_channel
    add_timestamp = FALSE,                  ## omit a date-time timestamp column
    zero_time = TRUE,                       ## recalculate time values from zero
    keep_all = FALSE,                       ## return only the specified data channels
    verbose = TRUE                          ## show warnings & messages
)
#> ! Estimated `sample_rate` = 2 Hz.
#> ℹ Define `sample_rate` explicitly to override.
#> Warning: ! Duplicate or irregular `time_channel` samples detected.
#> ℹ time = 211.59 and 1183.6.
#> ℹ Re-sample with `mnirs::resample_mnirs()`.
#> # A tibble: 2,202 × 3
#>     time smo2_left smo2_right
#>    <dbl>     <dbl>      <dbl>
#>  1 0            54         68
#>  2 0.560        54         68
#>  3 1.11         54         66
#>  4 1.66         54         66
#>  5 2.21         54         66
#>  6 2.76         54         66
#>  7 3.31         57         67
#>  8 3.86         57         67
#>  9 4.41         57         67
#> 10 4.96         57         67
#> # ℹ 2,192 more rows
```
