# Package index

## Read data

Read raw data from exported files.

- [`read_mnirs()`](https://jemarnold.github.io/mnirs/reference/read_mnirs.md)
  :

  Read *mnirs* data from file

- [`example_mnirs()`](https://jemarnold.github.io/mnirs/reference/example_mnirs.md)
  :

  Get path to *mnirs* example files

- [`create_mnirs_data()`](https://jemarnold.github.io/mnirs/reference/create_mnirs_data.md)
  :

  Create an *mnirs* data frame with metadata

## Process data

Clean, filter, and transform data.

- [`resample_mnirs()`](https://jemarnold.github.io/mnirs/reference/resample_mnirs.md)
  :

  Re-sample an *mnirs* data frame

- [`replace_mnirs()`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)
  [`replace_invalid()`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)
  [`replace_outliers()`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)
  [`replace_missing()`](https://jemarnold.github.io/mnirs/reference/replace_mnirs.md)
  :

  Replace outliers, invalid, and missing values in *mnirs* data

- [`filter_butterworth()`](https://jemarnold.github.io/mnirs/reference/filter_butterworth.md)
  [`filter_butter()`](https://jemarnold.github.io/mnirs/reference/filter_butterworth.md)
  : Apply a Butterworth digital filter

- [`filter_mnirs()`](https://jemarnold.github.io/mnirs/reference/filter_mnirs.md)
  : Filter a data frame

- [`filter_moving_average()`](https://jemarnold.github.io/mnirs/reference/filter_moving_average.md)
  [`filter_ma()`](https://jemarnold.github.io/mnirs/reference/filter_moving_average.md)
  : Apply a moving average filter

- [`shift_mnirs()`](https://jemarnold.github.io/mnirs/reference/shift_mnirs.md)
  : Shift data range

- [`rescale_mnirs()`](https://jemarnold.github.io/mnirs/reference/rescale_mnirs.md)
  : Rescale data range

## Detect intervals

Detect and extract intervals for further analysis

- [`extract_intervals()`](https://jemarnold.github.io/mnirs/reference/extract_intervals.md)
  :

  Extract intervals from *mnirs* data

- [`by_time()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
  [`by_label()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
  [`by_lap()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
  [`by_sample()`](https://jemarnold.github.io/mnirs/reference/by_time.md)
  : Specify interval boundaries by time, label, lap, or sample

## Print & plot

Tools for printing functions and pretty plotting with
[ggplot2](https://ggplot2.tidyverse.org)

- [`print(`*`<mnirs>`*`)`](https://jemarnold.github.io/mnirs/reference/print.mnirs.md)
  : Methods for mnirs objects

- [`plot(`*`<mnirs>`*`)`](https://jemarnold.github.io/mnirs/reference/plot.mnirs.md)
  :

  Plot *mnirs* objects

- [`theme_mnirs()`](https://jemarnold.github.io/mnirs/reference/theme_mnirs.md)
  :

  Custom *mnirs* ggplot2 theme

- [`palette_mnirs()`](https://jemarnold.github.io/mnirs/reference/palette_mnirs.md)
  :

  Custom *mnirs* colour palette

- [`scale_colour_mnirs()`](https://jemarnold.github.io/mnirs/reference/scale_colour_mnirs.md)
  [`scale_color_mnirs()`](https://jemarnold.github.io/mnirs/reference/scale_colour_mnirs.md)
  [`scale_fill_mnirs()`](https://jemarnold.github.io/mnirs/reference/scale_colour_mnirs.md)
  :

  Scales for custom *mnirs* palette

- [`breaks_timespan()`](https://jemarnold.github.io/mnirs/reference/breaks_timespan.md)
  : Breaks for time span data

- [`format_hmmss()`](https://jemarnold.github.io/mnirs/reference/format_hmmss.md)
  : Format time span data as h:mm:ss

## Files

Example data files included in the package

- [`moxy_intervals.csv`](https://jemarnold.github.io/mnirs/reference/moxy_intervals.csv.md)
  : 0.5 Hz Moxy onboard export
- [`train.red_intervals.csv`](https://jemarnold.github.io/mnirs/reference/train.red_intervals.csv.md)
  : 10 Hz Train.Red App export
- [`artinis_intervals.xlsx`](https://jemarnold.github.io/mnirs/reference/artinis_intervals.xlsx.md)
  : 10 Hz Artinis Oxysoft export recorded with Oxymon MKIII
- [`moxy_ramp.xlsx`](https://jemarnold.github.io/mnirs/reference/moxy_ramp.xlsx.md)
  : 2 Hz PerfPro export of Moxy data
- [`portamon-oxcap.xlsx`](https://jemarnold.github.io/mnirs/reference/portamon-oxcap.xlsx.md)
  : 10 Hz Artinis Oxysoft export recorded with Portamon
