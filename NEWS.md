# mnirs 0.7.1

* `shift_mnirs()` now properly excludes partial windows at data edges, where fewer samples can bias calculation of *"min"* or *"max"* shift values on noise.

* `replace_mnirs()`, `shift_mnirs()`, and `filter_moving_average()` should gain a serious performance improvement when calculating rolling means, at the cost of negligible loss of precision on the order of ± ~1e-11.

* Lists of data frames exported from core functions now contain `class = "mnirs"` and should now `plot()` properly.

# mnirs 0.7.0

## Highlights

This minor version update includes mostly internal refactoring, but enough user-facing changes and a few breaking deprecations, so that it's more than just a patch. This version lays the foundation for the incoming (hopefully soon) *mnirs* `analyse_kinetics` suite of functions.

## Working with lists and grouped data frames

* Core processing functions (`resample_mnirs()`, `replace_mnirs()`, `filter_mnirs()`, `shift_mnirs()`, and `rescale_mnirs()`) now accept as the `data` input either a single data frame, a list of data frames, or a grouped data frame (requires `{dplyr}`). Single data frames are processed and returned directly, as previously. Listed and grouped data frames are each processed independently and returned as a named list.

* `extract_intervals()` now also accepts a list of data frames or a grouped data frame and returns one flattened, named list of extracted intervals.

``` r
replace_mnirs(
    data = list(df1, df2),
    outlier_cutoff = 3,
    span = 5
)
#> $interval_1
#> # A tibble:
#>     time  smo2    o2hb
#>    <dbl> <dbl>   <dbl>
#>  1   0    42.8 -0.0289
#>  2   0.1  42.8 -0.0524
#>  3   0.2  42.8 -0.0916
#>  4   0.3  42.9 -0.138 
#>  5   0.4  43.2 -0.205 
#> 
#> $interval_2
#> # A tibble:
#>     time  smo2  o2hb
#>    <dbl> <dbl> <dbl>
#>  1   9.9  51.7 -2.29
#>  2  10    51.7 -2.32
#>  3  10.1  51.8 -2.31
#>  4  10.2  52.2 -2.22
#>  5  10.3  52.4 -2.12
```

## Channel grouping and processing

* **BREAKING CHANGE**: `extract_intervals()`, `shift_mnirs()`, and `rescale_mnirs()` gain argument `group_channels`, which separates channel selection from channel grouping. `nirs_channels` selects channels; the new `group_channels` argument channels specifies *"distinct"*, *"ensemble"*, or custom channel groupp construction.

``` r
shift_mnirs(
    data,
    nirs_channels = c(smo2, o2hb),
    group_channels = "ensemble",
    to = 0,
    span = 5
)
```

* `filter_mnirs()` and `replace_mnirs()` now accept named lists for channel-specific processing arguments, in addition to a single global value as previously. `rescale_mnirs()` and `shift_mnirs()` accept processing arguments specified by channel or group. This way each `nirs_channel` can be processed with unique parameters.

``` r
filter_mnirs(
    data,
    nirs_channels = c(smo2, o2hb),
    method = list(smo2 = "moving_average", o2hb = "butterworth"),
    span = 5,     ## only used by "moving_average" channels
    order = 2,    ## only used by "butterworth" channels (and below args)
    W = 0.02,
    type = "low",
    na.rm = TRUE
)
```

## Function and argument renaming

* `filter_butterworth()` and `filter_moving_average()` are now the renamed canonical functions. `filter_butter()` and `filter_ma()` remain available as aliases.

* **BREAKING CHANGE**: `extract_intervals()` new argument `group_intervals` replaces deprecated `event_groups`. This is where interval grouping is specified for *"distinct"* or *"ensemble"*-averaged, or custom group list construction.

``` r
interval_list <- extract_intervals(
    data,
    nirs_channels = c(smo2_left, smo2_right),
    group_intervals = "ensemble",
    start = by_time(368, 1084),
    span = c(-20, 90)
)
#> $ensemble
#> # A tibble:
#>     time smo2_left smo2_right
#>    <dbl>     <dbl>      <dbl>
#>  1 -20        56.3       59.2
#>  2 -19.9      56.1       59.2
#>  3 -19.8      56.1       59.2
#>  4 -19.7      56.2       58.9
#>  5 -19.6      56.4       58.9
```


## Performance improvements & bug fixes

### `read_mnirs()`

* Improved reading of delimited files containing extra spaces around quoted values.
* Documented support for `.txt` files.
* Made automatic device and channel detection more efficient and improved related messages.
* `nirs_channel` now identified before type conversion, ensuring coercion to numeric.
* Added a warning when all values in `nirs_channels` become missing during numeric coercion.
* Improved handling of timestamp values; absolute date-times, time-only values, header timestamps, and POSIXct-type `time_channel` column.
* Refined `sample_rate` detection and warnings for irregular samples.

### `create_mnirs_data()`

* Grouped data now remain grouped when *mnirs* metadata is added or refreshed.

### `resample_mnirs()`

* `data` arg now accepts a list or grouped data frame and returns one processed data frame per data frame.
* Clarified documentation that downsampling numeric columns uses averages based on linear interpolation, rather than time-weighted averages.

### `replace_mnirs()`

* `data` arg now accepts a list or grouped data frame and returns one processed data frame per data frame.
* Arguments can be specified per-`nirs_channels` as named lists.
* Improved performance via local-median calculations when values are missing.
* Made fixed-width outlier detection faster by calculating rolling medians across many windows together.
* Added clearer checks for invalid window settings and unsorted or missing time values.

### `filter_mnirs()`

* `data` arg now accepts a list or grouped data frame and returns one processed data frame per data frame.
* Arguments can be specified per-`nirs_channels` as named lists.
* Filter settings such as `spar`, `W`, `fc`, `width`, and `span` are now explicit function arguments, making available options easier to discover.
* Simplified filter selection so each channel is sent directly to its chosen filter method.
* `filter_butterworth()` and `filter_moving_average()` are now the main worker functions; `filter_butter()` and `filter_ma()` remain available as aliases.
* Improved checks and messages for missing values, cutoff frequencies, window sizes, and unsupported settings.

### `shift_mnirs()`

* `data` arg now accepts a list or grouped data frame and returns one processed data frame per data frame.
* Arguments can be specified per-`nirs_channels` as named lists, or per-named group via `group_channels` (see below).
* **BREAKING CHANGE**: Channel grouping through a list to process channels together is now supplied to `group_channels`. Supplying a list to `nirs_channels` is deprecated.
* `nirs_channels` now defaults to `NULL`, allowing channels to be taken from *mnirs* metadata as elsewhere in the package.

### `rescale_mnirs()`

* `data` arg now accepts a list or grouped data frame and returns one processed data frame per data frame.
* Arguments can be specified per-`nirs_channels` as named lists, or per-named group via `group_channels` (see below).
* **BREAKING CHANGE**: Channel grouping through a list to process channels together is now supplied to `group_channels`. Supplying a list to `nirs_channels` is deprecated.
* `nirs_channels` now defaults to `NULL`, allowing channels to be taken from *mnirs* metadata as elsewhere in the package.

### `extract_intervals()`

* `data` arg now accepts a list or grouped data frame and returns *one named list* with interval names that identify their source data frame and interval sequence number.
* Arguments can be specified per-`nirs_channels` as named lists.
* **BREAKING CHANGE**: Interval grouping for ensemble-averaging is now supplied to `group_intervals`, renamed from `event_groups`. `event_groups` is now deprecated.
* **BREAKING CHANGE**: Channel grouping through a list to process channels together is now supplied to `group_channels`. Supplying a list to `nirs_channels` is deprecated.
* Added checks for interval groups and their selected channels, with clearer errors for invalid group numbers, unknown channels, and invalid time boundaries.
* Updated ensemble averaging so each interval group can use its own per-channel argument selection.
* Improved naming, zero-time handling, and preservation of *mnirs* metadata in extracted intervals.

### `plot.mnirs()`

* `data` arg accepts a list or grouped data frame and returns facetted plots per data frame (existing functionality since *mnirs 0.6.3*, documenting for clarity).
* Time (x-) axis now displays units `mm:ss` for data lasting less than one hour and `h:mm:ss` for longer data.
* Facets now follow the order of intervals in the input instead of being reordered alphabetically.
* Improved error messages for empty or invalid lists of data frames.

### Internal validation & shared functions

* Improved errors and warnings package-wide so they point to the function called by the user rather than an internal checking function.
* Added `as_data_list.R` containing shared support for processing `data` argument from one data frame, a list of data frames, or a grouped data frame. Existing *mnirs* metadata is kept for each interval.
* Added `channel_args.R` containing shared support for setting processing options separately for each `nirs_channel` or per-group (for `shift_mnirs()`, `replace_mnirs()`, and `extract_intervals()`) as a list named by `nirs_channels`. One (vector) value can still be supplied to all channels by default.
* Added clear warnings for unknown channel names and errors for conflicting values within a group.




# mnirs 0.6.5

## `read_mnirs()`

* `read_mnirs()` can now handle `c(".txt", ".tsv")` files via the same `data.table::fread()` pathway as `".csv"`. There are occasionally still odd file formats where columns will not be properly detetcted, usually where files are designed to be human-readible not machine-readable.

* `read_mnirs(nirs_channels = NULL)` now automatically returns **all** channels starting with "SmO2" (case insensitive), which is the most common NIRS channel name for wearable mNIRS devices. Previously, only the first detected nirs channel was returned.

``` r
## read an mNIRS file with two "smo2" channels
df <- read_mnirs(file_path = example_mnirs("moxy_ramp"))

attr(df, "nirs_channels")
#> [1] "SmO2 Live"    "SmO2 Live(2)"
```

## Core function updates

* `extract_intervals()` helper function `by_label()` matches character strings as regular expressions (regex) by default. Now more clearly documents this, and accepts two additional arguments for better handling character strings, which are both `FALSE` be default and must be explicitly opted into (see `?grep`):

    * `by_label(ignore_case = TRUE)` ignore case when evaluating strings.

    * `by_label(fixed = TRUE)` treat labels as fixed strings rather than regular expressions. Useful when labels contain regex metacharacters (`.`, `*`, `(`, etc.).

* `plot.mnirs()` internal update to pass additional arguments `n.breaks` to y-axis breaks, and `breaks` to x-axis breaks.

* `format_hmmss()` can now display fractional seconds values.

* `replace_mnirs()`, `resample_mnirs()`, and some other functions return more informative error message when `time_channel` has irregularly sampled values (time series must be monotonically increasing, non `NA`).


# mnirs 0.6.4

## Bug fix

* Fix a memory allocation bug with `resample_mnirs`. Downsampling 100k rows now allocates ~11MB rather than ~4.5 GB.


# mnirs 0.6.3

## Plotting improvements

* Generic `plot.mnirs()` ncan now be called on a list of *"mnirs"* data frames. Each data frame will be printed as a facet. This is primarily useful for printing a list of interval data frames exported from `extract_intervals()`:

``` r
read_mnirs() |> 
    extract_intervals() |> 
    plot()
## returns a plot with a facet for each interval
```

* To faciliate this, `extract_intervals()` now returns a list of data frames with `class = "mnirs"`. Otherwise, to manually plot a list of data frames, it will need to have `class(list) <- c("mnirs", class(list))` edited manually.

``` r
read_mnirs() |> 
    extract_intervals() |> 
    class()
#> [1] "mnirs" "list"
```

* Backend improvement: `print.mnirs()` generic created to avoid displaying extra `NextMethod()` calls when printing lists with `class = "mnirs"`.

## Modified lap extraction behaviour

* `extract_intervals()`: 

    * When specifying `start` and `end` with `by_lap()`, `start` will only refer to the first sample of the lap, and `end` the last sample.

    * With previous behaviour, `start = by_lap()` would include the entire specified lap(s). But this resulted in less control over displaying only parts of a lap.
    
    * Updated behaviour is now more consistent with `by_time`, `by_label`, and `by_sample()` methods, which reference `start` and `end` from a single sample. This allows extraction of e.g. only the first 60-sec of lap:

```r
read_mnirs() |> 
    extract_intervals(
        start = by_lap(1, 3),
        span = c(0, 60),
    )
## returns a list of two intervals with the first 60-sec of laps 1 and 3, respectively.
```

## Updated core functions

* `read_mnirs()`: Fix detection issue with *"PerfPro"* file formats, and another small bug fix to improve timestamp parsing.

* `create_mnirs_data()` now properly accepts tidy evaluation for `nirs_channels`, `time_channel`, and `event_channel`:

```r
create_mnirs_data(df, nirs_channels = c(o2hb, hhb))
```

* `resample_mnirs()`: 

    * Now resamples to an inclusive time range around existing data, rounded to the nearest resampled rate. Better handles edge cases where the last sample was being dropped in certain rounding conditions.

    * Fix an edge case error when `sample_rate` was mis-specified higher than the actual sample_rate of the data. Now more robustly fills non-numeric columns.


# mnirs 0.6.2

## Core updates

* `read_mnirs()`: 
    
    * Now properly parses `time_channel` as fractional unix values; i.e. timestamp (e.g. "hh:mm:ss") values are saved by Excel in all its infinite wisdom as numeric fractional Unix timestamps. Will now be properly coerced to numeric and POSIXct timestamp values can be returned.

    * Timestamps should now be returned in the user's local time zone.

## Core function argument changes

* `resample_mnirs()`: Update default `method = "none"`. Less opinionated default to force users to explicitly opt-in to specifying either "linear" or "locf" methods to fill/interpolate across new samples. Updated package documentation.

* `replace_invalid()`, `replace_outliers()`, `replace_missing()`: Remove `bypass_checks` arg intended for internal use only, to bypass redundant checks when calling from `replace_mnirs()`.

## Small edits

* `plot.mnirs()`: No longer coerces to long format data behind the scenes.

* `plot.mnirs()`: y-axis title changed from "signal" to "mNIRS".

* Fix lap marker inconsistency in `train.red_intervals.csv`. Updated relevant interval times in documentation.

* `README.md` & *"reading-mnirs-data.qmd"* vignette updates.

    * Update recommended core processing sequence: `read_mnirs()` -> `resample_mnirs()` -> `replace_mnirs()` -> ...

    * Update `extract_intervals()` vignette to `train.red_intervals.csv` end-interval reoxygenation events.

* Small documentation changes.


# mnirs 0.6.1

## Bug fixes

* Core functions updated to accept `nirs_channel` args as a list. 
    
    * An info message will be displayed when a list is not required, instead of erroring.
    
    * Additional info messages will be displayed for `shift_mnirs()`, `rescale_mnirs()`, and `extract_intervals()` when `nirs_channels` has not been specified as a list, nor retrieved from metadata. As a reminder of the grouping requirements in those functions.

* Core functions now properly update `nirs_channels` metadata when re-specified.
    
    * Previously, specifying `nirs_channels` in a function would only add any additional column name strings to the existing metadata rather than overwrite it. Meaning `nirs_channels` could only be removed from metadata by using `create_mnirs_data(nirs_channels = "...")`. The updated behaviour should mean channels need to be re-specified less often.

* `filter_mnirs()`: Fixed an error matching `method` arguments when left blank.

* *README* and *"Reading and Cleaning Data with {mnirs}"* vignette updated to fix typos and small edits.

# mnirs 0.6.0

## Prepare for CRAN submission

* *NEWS.md* is truncated to the most recent relevant release updates. The full *NEWS.md* remains in the `dev` branch.

* Some development functions have been omitted from the package build in expectation of CRAN review. They are still present in `dev` branch and can be installed with `pak::pak("jemarnold/mnirs")`.

### Updated core functions

* `filter_mnirs()`: Method-specific arguments (e.g. `order` for `method = "butterworth"`, or `width` for `method = "moving_average"`) removed from the generic function call. Continue to be passed to the appropriate method via `...`. Documentation and info/warning/abort messages updated.

* `filter_ma()` better separates effects of `partial` and `na.rm` args:
    
    * `partial = FALSE` by default returns NA at edges where insufficient number of samples are available compared to the specified `width` or `span`.
    
    * `partial = TRUE` calculates mean values at edges, as long as one valid non-`NA` sample is available.
    
    * `na.rm = FALSE` by default behaves as expected with `mean(na.rm = FALSE)`, propagating any `NA`s in the local window to the calculated mean with a warning.
        
        * **NOTE** This differs from the behaviour of `na.rm = FALSE` in `filter_butter()`, which errors if there are any internal `NA`s present. This has not been changed.
    
    * `na.rm = TRUE` ignores `NA`s and calculates local means as long as one valid sample is present.

* `plot.mnirs()`
    
    * Add new arg: `points = TRUE` will plot points in addition to lines, as a useful quick shortcut.
    
    * Update `na.omit` now omits non-valid `c(NA, NaN, Inf, -Inf)` values from plotting, not just `NA`.

* `palette_mnirs()` can now accept either a single numeric value specifying the number of colours to return, or any number of (valid) character colour names.


# mnirs 0.5.2

### Shiny app

* Update online shiny app hosted at <https://jemarnold-mnirs-app.share.connect.posit.cloud/> with basic reading and pre-processing functionality.

# mnirs 0.5.1

### Updated core functions

* `read_mnirs()` now has better automated channel detection logic for *"Artinis" / "Oxysoft"* file formats. `nirs_channels` and `time_channel` can be left blank to automatically *"sample"*, *"time"* (from `sample_rate`), and nirs channel *"2"*.

* `replace_outliers()` now accepts numeric `outlier_cutoff` values for more precise outlier detection thresholds (previously was integers only). Documentation also updated.


# mnirs 0.5.0

### Updated core functions

* `read_mnirs()` expands `event_channel` to work with integer *"lap"* numbers, or character event label as previous.
    
    * Should now work on more .csv file formats; previously read errors may have occured where the file contained header rows above the data table, resulting in improper detection of columns.
    
    * `event_channel` can now be specified as an integer `lap` column, in addition to a character column as previous.
    
    * Other {mnirs} functions may expect `event_channel` to be either character or integer-ish.

* `extract_intervals()`:
    
    * Function arguments `start` and `end` are used to specify one or both of a start and end point to the target interval.
    
    * Specify `start`/`end` values with helper functions `by_time()`, `by_label()`, `by_lap()`, and `by_sample()`.
    
    * Numeric values are automatically coerced to "time" values; Explicit integer values (e.g. `2L`) are coerced to "lap"; Character strings are coerced to event "label".

### Package resources

* A rough draft *{mnirs}* hex icon and package cheatsheet have been added.


# mnirs 0.4.2 and prior

* Full previous changelog is available on the github `dev` branch: https://github.com/jemarnold/mnirs/blob/dev/NEWS.md
