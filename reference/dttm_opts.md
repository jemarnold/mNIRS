# Datetime format strings for POSIXct parsing

Time-only format must stay first:
[`parse_time_channel()`](https://jemarnold.github.io/mnirs/reference/parse_time_channel.md)
splits on `dttm_opts[1L]` vs `dttm_opts[-1L]` to detect an absolute
date-time series.

## Usage

``` r
dttm_opts
```
