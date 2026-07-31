# coerce raw values to mnirs_interval objects

coerce raw values to mnirs_interval objects

## Usage

``` r
as_mnirs_interval(x, arg = "start", env = rlang::caller_env())
```

## Arguments

- x:

  A raw value or mnirs_interval object.

- arg:

  Name of the argument for error messages.

- env:

  The calling environment or a defused call, used to report errors and
  warnings as coming from the user-facing function rather than the
  validator.
