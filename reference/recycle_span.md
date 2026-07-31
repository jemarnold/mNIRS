# recycle a single-element span to c(before, after) positive -\> c(0, x), negative -\> c(x, 0)

recycle a single-element span to c(before, after) positive -\> c(0, x),
negative -\> c(x, 0)

## Usage

``` r
recycle_span(span, env = rlang::caller_env())
```

## Arguments

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
