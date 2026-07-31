# Validate and normalise channel grouping

Converts the `group_channels` argument to a named list of channel-name
vectors. String shortcuts expand against `nirs_channels`: `"ensemble"`
places all channels in one group (preserving relative scaling) and
`"distinct"` places each channel in its own group. Custom
[`list()`](https://rdrr.io/r/base/list.html) groupings may use bare
symbols or character names. Groups must be non-empty and their resulting
names must be unique. Channels omitted from a custom grouping are
processed independently, matching `group_intervals` behaviour in
[`extract_intervals()`](https://jemarnold.github.io/mnirs/reference/extract_intervals.md).

## Usage

``` r
validate_group_channels(
  nirs_channels,
  group_channels,
  data = NULL,
  env = rlang::caller_env()
)
```

## Arguments

- nirs_channels:

  Character vector of resolved channel names.

- group_channels:

  A quosure from
  [`rlang::enquo()`](https://rlang.r-lib.org/reference/enquo.html), a
  character string (`"ensemble"` or `"distinct"`), or a
  [`list()`](https://rdrr.io/r/base/list.html) of (optionally named)
  channel-name vectors.

- data:

  A data frame for parsing bare-symbol group members.

- env:

  Environment for symbol evaluation.

## Value

A uniquely named list of non-empty character vectors covering all
`nirs_channels`, each channel appearing in exactly one group.
