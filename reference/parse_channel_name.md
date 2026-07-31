# Parse channel expressions for NSE

Converts quosures to character vectors, handling bare symbols, character
strings, lists, and tidyselect expressions.

## Usage

``` r
parse_channel_name(channel, data, env = rlang::caller_env())
```

## Arguments

- channel:

  A quosure from
  [`rlang::enquo()`](https://rlang.r-lib.org/reference/enquo.html).

- data:

  A data frame for tidyselect context.

- env:

  Environment for symbol evaluation (typically the quosure environment).

## Value

A character vector, list of character vectors, or `NULL`.
