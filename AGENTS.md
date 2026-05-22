# `{mnirs}` Agent Reference

**v0.6.2 | R | MIT** — workflow/dependency map. See `README.md` + vignettes for examples.

| | |
|---|---|
| **Author** | Jem Arnold |
| **Website** | https://jemarnold.github.io/mnirs/ |
| **Citation** | \<coming soon\> |

---

## 1. `"mnirs"` — data frame/{tibble} subclass

### Metadata

| Attribute | Type | Description |
|---|---|---|
| `nirs_channels` | character vector | NIRS signal column names |
| `time_channel` | character(1) | time column name |
| `event_channel` | character(1) | event/lap column name |
| `nirs_device` | character(1) | device name (auto-detected) |
| `sample_rate` | numeric(1) | Hz |
| `start_timestamp` | POSIXct | absolute start datetime |
| `interval_times` | numeric | set by `extract_intervals()` |
| `interval_span` | numeric(2) | span used in `extract_intervals()` |

`verbose` read from `getOption("mnirs.verbose", TRUE)` or pass explicitly.
Access: `attr(data, "nirs_channels")`.

---

## 2. Pipeline

```
read_mnirs()
  └─ resample_mnirs()                  # regularise time grid
      └─ replace_mnirs()               # clean invalid/outliers/NA
          └── filter_mnirs()           # smooth
              ├── shift_mnirs()        # optional: shift baseline
              ├── rescale_mnirs()      # optional: normalise range
              └── extract_intervals()
                └── analyse_kinetics()
                  └── plot() / print()
```

**Order rules:**
- `resample` → `replace`: window replacement needs regular `time_channel`
- `replace` → `filter`: `"smooth_spline"` / `"butterworth"` error on internal `NA` (unless `na.rm = TRUE`)
- `resample` → `filter("smooth_spline")`: errors on duplicated time values
- `resample` → `extract_intervals(event_groups = "ensemble")`: ensemble needs regular samples
- `extract_intervals` → `analyse_kinetics` for multiple intervals

---

## 3. Function Reference

### 3.1 Read

```r
read_mnirs(
    file_path,
    nirs_channels = NULL,   # char vec; named = rename: c(new = "old")
    time_channel  = NULL,   # char(1); named = rename: c(time = "Timestamp")
    event_channel = NULL,   # char(1); optional; named = rename
    sample_rate   = NULL,   # numeric(1) Hz; estimated if NULL
    add_timestamp = FALSE,  # add POSIXct "timestamp" column
    zero_time     = FALSE,  # rebase time[1] to 0
    keep_all      = FALSE,  # keep all columns, otherwise only specified channels
    verbose       = TRUE
)

example_mnirs(file = NULL)
## files: "moxy_ramp", "train.red", "artinis", "portamon", "vo2master"
## blank → returns all names; partial matching; pass to read_mnirs(file_path = ...)

create_mnirs_data(data, ...)  # low-level constructor; wraps df as "mnirs"
```

---

### 3.2 Resample

```r
resample_mnirs(
    data,
    time_channel  = NULL,
    sample_rate   = NULL,         # source Hz; estimated if NULL
    resample_rate = sample_rate,  # target Hz; default = source
    method = c("none", "locf", "linear"),
    verbose = TRUE
)
```

- `method = "none"`: new samples → `NA`
- Non-numeric cols always use `"locf"`
- `resample_rate = sample_rate`: regularise only, no rate change

---

### 3.3 Clean

```r
replace_mnirs(
    data,
    nirs_channels  = NULL,
    time_channel   = NULL,
    invalid_values = NULL,   # numeric vec; exact values to replace
    invalid_above  = NULL,   # numeric(1); replace x > threshold
    invalid_below  = NULL,   # numeric(1); replace x < threshold
    outlier_cutoff = NULL,   # numeric(1); Hampel MAD multiplier; NULL = skip
    width          = NULL,   # integer; rolling window in samples
    span           = NULL,   # numeric; rolling window in time units
    method = c("linear", "median", "locf", "none"),
    verbose = TRUE
)
```

Processing order: invalid → outliers → missing (NA).

Vector-level functions:
```r
replace_invalid(x, t, invalid_values, invalid_above, invalid_below,
                width, span, method = c("median", "none"), ...)
replace_outliers(x, t, outlier_cutoff = 3, width, span,
                 method = c("median", "none"), ...)
replace_missing(x, t, width, span,
                method = c("linear", "median", "locf"), ...)
```

`replace_invalid()` / `replace_outliers()` default `method = "median"`.
`replace_mnirs()` / `replace_missing()` default `"linear"`.

---

### 3.4 Filter

```r
filter_mnirs(
    data,
    nirs_channels = NULL,
    time_channel  = NULL,
    method = c("smooth_spline", "butterworth", "moving_average"),
    na.rm   = FALSE,       # bypass NAs
    verbose = TRUE,
    ...
)
```

**Aliases:**

| Canonical | Aliases |
|---|---|
| `"smooth_spline"` | `"spline"`, `"smooth spline"`, `"smooth-spline"` |
| `"butterworth"` | `"butter"` |
| `"moving_average"` | `"ma"`, `"moving average"`, `"moving-average"` |

**Method-specific arguments:**

**`"smooth_spline"`** (`stats::smooth.spline()`): `spar` (NULL = GCV auto); errors on NA/duplicated time.

**`"butterworth"`** (requires `{signal}`): `order` (default `2L`), `W`/`fc` (normalised or Hz), `type` (`"low"`, `"high"`, `"stop"`, `"pass"`), `edges` (`"rev"`, `"rep1"`, `"none"`); errors on NA.

**`"moving_average"`**: `width` (samples) or `span` (time units); `partial` (default `FALSE`).

Vector-level:
```r
filter_ma(x, t, width, span, partial = FALSE, na.rm = FALSE, ...)
filter_moving_average()    # alias for filter_ma()
filter_butter(x, order = 2L, W, type = "low", edges = "rev", na.rm = FALSE)
```

---

### 3.5 Transform (optional)

```r
shift_mnirs(
    data,
    nirs_channels = list(),  # list()/NULL → auto from metadata
    time_channel  = NULL,
    to       = NULL,         # numeric(1); target level (overrides `by`)
    by       = NULL,         # numeric(1); shift amount
    width    = NULL,         # integer; window in samples
    span     = NULL,         # numeric; window in time units
    position = c("min", "max", "first"),
    verbose  = TRUE
)

rescale_mnirs(
    data,
    nirs_channels = list(),
    range,                   # numeric(2): c(min, max)
    verbose = TRUE
)
```

**`nirs_channels` list semantics** (both functions):

| Syntax | Behaviour |
|---|---|
| `list("A", "B")` | channels processed independently |
| `list(c("A", "B"))` | shared reference (relative scaling preserved) |
| `list(c("A", "B"), c("C"))` | A+B share reference; C independent |

Default when retrieved from metadata: `list(c("A", "B"))`.

---

### 3.6 Extract Intervals

```r
extract_intervals(
    data,
    nirs_channels = NULL,
    time_channel  = NULL,
    event_channel = NULL,
    sample_rate   = NULL,
    start         = NULL,   # by_time(numeric), by_label(character), by_lap(integer), by_sample(numericfs)
    end           = NULL,   # same; NULL = use span around start
    span          = list(c(-60, 60)),  # c(before, after) per interval
    event_groups  = c("distinct", "ensemble"),
    zero_time     = FALSE,  # rebase time[1] to 0 per interval
    verbose       = TRUE
)
```

Returns named list of `"mnirs"` dfs.

- Raw coercion: numeric → `by_time()`; character → `by_label()`; explicit integer → `by_lap()`
- `start = by_lap(), end = NULL`: full lap range
- `start = by_lap(), end = by_lap()`: first start sample → last end sample
- `span = c(before, after)` recycled to boundaries; negative = before start, positive = after end
- `span = list()` recycled to intervals
- `event_groups = "distinct"`: one df per event
- `event_groups = "ensemble"`: single ensemble-averaged df (requires regular time grid)
- `event_groups = list(c(1,2), c(3,4))`: one ensemble df per group
- Names: `"ensemble"` or `"interval_1"`, `"interval_2"`, ...

---

### 3.7 Analyse Kinetics

```r
analyse_kinetics(
    data,           # "mnirs" df | named list of "mnirs" dfs | {dplyr} grouped df
    nirs_channels  = NULL,
    time_channel   = NULL,
    method = c("response_time", "peak_slope", "monoexponential", "sigmoidal"),
    direction = c("auto", "positive", "negative"),
    end_fit_span   = Inf,    # truncate fit after extreme; Inf = global extreme
    channel_args   = list(), # per nirs_channel argument overrides
    verbose        = TRUE,
    ...
)
```

#### `"mnirs_kinetics"` structure — list (from `analyse_kinetics()`)

| Element | Type | Description |
|---|---|---|
| `method` | character | method used |
| `model` | named list | per-interval per-channel model objects |
| `coefficients` | data frame | one row per channel per interval |
| `data` | named list | input `"mnirs"` dfs with `<channel>_fitted` cols |
| `interval_times` | data frame | `interval`, `interval_times` from metadata |
| `diagnostics` | data frame | `n_obs`, `r2`, `adj_r2`, `pseudo_r2`, `rmse`, `snr`, `cv_rmse` |
| `channel_args` | data frame | resolved args per channel per interval |
| `call` | call | matched function call |

**Aliases:**

| Canonical | Aliases |
|---|---|
| `"response_time"` | `"half time"`, `"response\|recovery time"`, `"HRT"` |
| `"peak_slope"` | `"peak slope"`, `"slope"` |
| `"monoexponential"` | `"monoexp"`, `"exponential"`, `"MRT"`, `"tau"` |
| `"sigmoidal"` | `"logistic"`, `"xmid"` |

**Input:** single df → `"interval_1"`; named list → each separate; grouped df → split by group.

**Per-method args `...`:**

`"response_time"`: `t0` (default `0`), `fraction` (default `0.5`; `0.632` ≈ MRT)

`"peak_slope"`: `width` or `span` (one required); `align` (`"centre"`, `"left"`, `"right"`); `partial` (default `FALSE`)

`"monoexponential"`: `use_time_delay` (default `TRUE`; 4-param → 3-param fallback); accepts `stats::nls()` args

**Per-channel overrides:**
```r
analyse_kinetics(data, nirs_channels = c(hhb, smo2), method = "peak_slope",
    span = 5, channel_args = list(smo2 = list(span = 10), hhb = list(direction = "negative")))
```
Names must match `nirs_channels` exactly; unrecognised names silently ignored.

**Coefficients:**

`"response_time"`: `interval`, `nirs_channels`, `time_channel`, `A`, `B`, `response_time`, `response_value`, `fitted`, `idx`

`"peak_slope"`: `interval`, `nirs_channels`, `time_channel`, `slope`, `intercept`, `fitted`, `<time_channel>`, `idx`

`"monoexponential"`: `interval`, `nirs_channels`, `time_channel`, `A`, `B`, `tau`, `k`, `TD`, `MRT`, `HRT`, `tau_fitted`, `MRT_fitted`, `HRT_fitted`

**Diagnostics:** `n_obs`, `r2`, `adj_r2`, `pseudo_r2`, `rmse`, `snr`, and `cv_rmse`

**Vector-level:**
```r
monoexponential(t, A, B, tau, TD)
## 3-parameter equation: `A + (B - A) * (1 - exp(-t / tau))`
## 4-parameter equation: `ifelse(t <= TD, A, A + (B - A) * (1 - exp(-(t - TD) / tau)))`

nls(x ~ SS_monoexp3(t, A, B, tau), data = df)
nls(x ~ SS_monoexp4(t, A, B, tau, TD), data = df)

peak_slope(x, t = seq_along(x), width = NULL, span = NULL,
    align = c("centre", "left", "right"),
    direction = c("auto", "positive", "negative"),
    partial = FALSE, na.rm = FALSE, verbose = TRUE)
## returns: slope, intercept, y, t, idx, fitted, window_idx, model

response_time(x, t = seq_along(x), t0 = 0, fraction = 0.5,
    end_fit_span = Inf, direction = c("auto", "positive", "negative"),
    verbose = TRUE)
## returns: A, B, response_time, response_value, fitted,
##          baseline_idx, response_idx, extreme_idx
```

---

### 3.8 Plot and Print

```r
plot.mnirs(x, points = FALSE, time_labels = FALSE, na.omit = FALSE)
## requires {ggplot2}; called via plot(data); returns ggplot2 object
## na.omit = TRUE connects lines across NAs

print(result)  # "mnirs_kinetics"; formatted coefficient table (max 10 rows)

theme_mnirs(base_size = 14, base_family = "sans",
            border = c("partial", "full"), ...)

palette_mnirs()              # all 12 named colours
palette_mnirs(4)             # first 4
palette_mnirs("red", "blue") # by name

scale_colour_mnirs(...)      # alias: scale_color_mnirs()
scale_fill_mnirs(...)
breaks_timespan(unit = "secs", n = 5)
format_hmmss(x)              # numeric seconds → "mm:ss" or "h:mm:ss"
```

---

## 4. Dependencies

| Package | Role | Type | Condition |
|---|---|---|---|
| `cli` | user messages | Import | |
| `data.table` | data manipulation | Import | |
| `lifecycle` | deprecation | Import | |
| `readxl` | XLS/XLSX | Import | |
| `rlang` | NSE/tidy eval | Import | |
| `tibble` | `"mnirs"` class | Import | |
| `tidyselect` | column selection | Import | |
| `signal` | Butterworth | Suggests | `"butterworth"` method |
| `ggplot2` | plotting | Suggests | `plot.mnirs()` |
| `scales` | axis formatting | Suggests | `plot.mnirs()` |
| `dplyr` | grouped df input | Suggests | grouped `analyse_kinetics()` |
| `zoo` | rolling functions | Suggests | internal testing |

---

## 5. Constraints

| Constraint | Detail |
|---|---|
| `example_mnirs()` returns path | pass to `read_mnirs(file_path = example_mnirs(...))` |
| Irregular samples warning from `read_mnirs()` | fires in pipe even after downstream `resample_mnirs()`; verify by inspecting output |
| Pipeline order | `resample → replace → filter`; wrong order = incorrect results |
| `replace_invalid` / `replace_outliers` default `"median"` | `replace_mnirs()` default `"linear"` |
| `"smooth_spline"` / `"butterworth"` fail on NA | use `replace_mnirs()` first or `na.rm = TRUE` |
| `"smooth_spline"` fails on duplicated time | use `resample_mnirs()` first |
| Ensemble needs regular samples | `event_groups = "ensemble"` warns if irregular |
| `nirs_channels = list()` default | auto-retrieves and groups all channels; see §3.5 |
| `channel_args` names must match exactly | unrecognised names silently ignored |
| `monoexponential` fallback | 4-param → 3-param → `"response_time"` with warnings |
| Single df to `analyse_kinetics()` | `interval = "interval_1"` in coefficients |

---

## 6. Key Source Files

| File | Contents |
|---|---|
| `R/read_mnirs.R` | `read_mnirs()`, `example_mnirs()`, `create_mnirs_data()` |
| `R/resample_mnirs.R` | `resample_mnirs()` |
| `R/replace_mnirs.R` | `replace_mnirs()`, `replace_invalid()`, `replace_outliers()`, `replace_missing()` |
| `R/filter_mnirs.R` | `filter_mnirs()`, `filter_ma()`, `filter_butter()` |
| `R/shift_mnirs.R` | `shift_mnirs()` |
| `R/rescale_mnirs.R` | `rescale_mnirs()` |
| `R/extract_intervals.R` | `extract_intervals()`, `by_time()`, `by_label()`, `by_lap()`, `by_sample()` |
| `R/analyse_kinetics.R` | `analyse_kinetics()`, `compute_diagnostics()` |
| `R/analyse_peak_slope.R` | `peak_slope()`, `rolling_slope()` |
| `R/analyse_monoexponential.R` | `monoexponential()`, `SS_monoexp3()`, `SS_monoexp4()` |
| `R/analyse_response_time.R` | `response_time()` |
| `R/plot.mnirs.R` | `plot.mnirs()`, `theme_mnirs()`, `palette_mnirs()`, scale functions |
| `R/validate_mnirs.R` | input validation |
| `R/data.R` | example file descriptions |
