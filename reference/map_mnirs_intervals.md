# Apply an mnirs function over each interval of a multi-interval input

Shared entry point for `*_mnirs()` transformer functions accepting a
list of data frames or a grouped data frame. Normalises `data` to a
named list via
[`as_data_list()`](https://jemarnold.github.io/mnirs/reference/as_data_list.md),
then re-evaluates the captured user-facing call once per interval with
`data` swapped, so all arguments (including NSE channel expressions) are
forwarded verbatim.

## Usage

``` r
map_mnirs_intervals(data, call, eval_env, env = rlang::caller_env())
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata, a list of data frames, or a grouped data frame (see
  *Details*).

- call:

  The matched call from the user-facing function, re-evaluated with
  `data` swapped for each interval data frame.

- eval_env:

  Environment in which to re-evaluate `call`, i.e. the user-facing
  function's [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html),
  so NSE arguments resolve against the original caller.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.

## Value

A named list of processed *"mnirs"* data frames, one per interval.

## Data input formats

*mnirs* processing functions accept `data` in multiple formats:

- A **single *"mnirs"* data frame** is processed and returned directly.

- A **list of *"mnirs"* data frames**: each interval is processed
  separately and returned as a named list.

- A **grouped *"mnirs"* data frame**, e.g. with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html):
  the data frame is split by grouping levels and each group is processed
  as a separate interval, returned as a named list.
