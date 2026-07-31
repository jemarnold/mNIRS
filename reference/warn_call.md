# trim caller call to bare function name for warning headers `env` accepts an environment or a call, e.g. from `sys.call(-1)`

trim caller call to bare function name for warning headers `env` accepts
an environment or a call, e.g. from `sys.call(-1)`

## Usage

``` r
warn_call(env = rlang::caller_env())
```
