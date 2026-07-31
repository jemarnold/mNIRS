# Read data table from raw data

Read data table from raw data

## Usage

``` r
read_data_table(
  data,
  header_row = 1L,
  nirs_channels = NULL,
  env = rlang::caller_env()
)
```

## Arguments

- data:

  A data frame of class *"mnirs"* containing time series data and
  metadata.

- nirs_channels:

  A character vector giving the names of mNIRS columns to operate on.
  Must match column names in `data` exactly.

  - If `NULL` (default), the `nirs_channels` metadata attribute of
    `data` is used.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
