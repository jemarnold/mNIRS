# Resolve per-channel arguments

Broadcasts global argument values across `nirs_channels`, applying
per-channel overrides where an argument is supplied as a named
[`list()`](https://rdrr.io/r/base/list.html) keyed by channel name. An
argument is treated as per-channel when it is a
[`list()`](https://rdrr.io/r/base/list.html) with at least one named
element, with at most one unnamed element acting as the fallback for
unlisted channels (e.g. `width = list(5, q = 7)` gives `q` 7 and every
other channel 5). Names must match `nirs_channels` or group names;
unrecognised names are warned about and ignored. Any other value
(unnamed vectors and fully unnamed lists) is applied globally to every
channel.

## Usage

``` r
resolve_channel_args(
  nirs_channels,
  group_channels = NULL,
  args,
  defaults = list(),
  choices = list(),
  verbose = TRUE,
  env = rlang::caller_env()
)
```

## Arguments

- nirs_channels:

  Character vector of resolved channel names.

- group_channels:

  An *optional* named list of channel-name vectors from
  [`validate_group_channels()`](https://jemarnold.github.io/mnirs/reference/validate_group_channels.md).
  When supplied, arguments are resolved per group: a group-name key or
  any member-channel key applies to the whole group, and conflicting
  member values within one group abort.

- args:

  Named list of per-channel-capable arguments. Each element is either a
  global value or a per-channel
  [`list()`](https://rdrr.io/r/base/list.html) map. A per-channel map
  may include a single unnamed element as the fallback for unlisted
  channels.

- defaults:

  Named list of fallback values per argument, used when a per-channel
  map omits a channel and supplies no unnamed fallback. Only needed for
  arguments whose formal default is not `NULL` (e.g.
  `method = "linear"`).

- choices:

  Named list of valid values for choice-type arguments (e.g.
  `list(method = c("linear", "median", "locf", "none"))`). Resolved
  values are matched per channel; a full default vector resolves to its
  first element, matching
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) behaviour.

- verbose:

  Logical. Default is `TRUE`. Display or silence (if `FALSE`) warnings
  and information messages helpful for troubleshooting. Ad global
  default can be set via `options(mnirs.verbose = FALSE)`.

- env:

  The calling environment, used to report errors as coming from the
  user-facing function (e.g.
  [`rescale_mnirs()`](https://jemarnold.github.io/mnirs/reference/rescale_mnirs.md)).

## Value

A named list with one element per channel (or per group when
`group_channels` is supplied); each element is a named list of that
channel's resolved argument values.
