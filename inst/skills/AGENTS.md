# `{mnirs}` Agent Reference

**v0.7.0 | R | MIT** — workflow/dependency map for muscle near-infrared
spectroscopy (mNIRS) processing & analysis package. 
Website: https://jemarnold.github.io/mnirs/.
See `README.md` + vignettes for examples.

Read `.xls(x)`/`.csv`/`.txt`/`.tsv` exports → resample → clean → filter → (transform)
→ extract intervals → analyse kinetics → plot.

---

## 1. `"mnirs"` — data frame/{tibble} subclass

Access metadata with `attributes(data)` or individually, eg: `attr(data, "nirs_channels")`.

| Attribute | Type | Description |
|---|---|---|
| `nirs_device` | character(1) | device name (auto-detected) |
| `nirs_channels` | character vector | NIRS signal column names |
| `time_channel` | character(1) | time column name |
| `event_channel` | character(1) | event/lap column name |
| `sample_rate` | numeric(1) | Hz |
| `start_timestamp` | POSIXct | absolute start datetime |
| `interval_times` | list | `list(start, end)`; set by `extract_intervals()` |
| `interval_span` | numeric(2) | span used in `extract_intervals()` |

`verbose` read from `getOption("mnirs.verbose", TRUE)` or passed explicitly.

---

## 2. Pipeline

```
read_mnirs()
└─ plot()                              # visualise "mnirs" data at each step
   └─ resample_mnirs()                 # regularise time grid
      └─ replace_mnirs()               # clean invalid/outliers/NA
         └─ filter_mnirs()             # smooth
            ├─ shift_mnirs()           # optional: shift baseline
            ├─ rescale_mnirs()         # optional: normalise range
            ├─ correct_blood_volume()  # optional: blood-volume normalise
            └─ extract_intervals()     # returns named list of "mnirs" dfs
               └─ analyse_kinetics()   # compute response rate & time course
                  └─ plot()            # plot "mnirs_kinetics": fit + markers
```

**Order:** diagram nesting = required sequence; failure modes + fixes in §5.
`extract_intervals` → `analyse_kinetics` (or other list-wise fns) for
multiple intervals.

### Data input formats

`resample_mnirs()`, `replace_mnirs()`, `filter_mnirs()`, `shift_mnirs()`,
`rescale_mnirs()`, `extract_intervals()`, and `analyse_kinetics()` accept
`data` as:

- **single `"mnirs"` df** → processed, returned directly.
- **list of `"mnirs"` dfs** (e.g. from `extract_intervals()`) →
  each processed separately; returns a list.
- **grouped df** (`dplyr::group_by()`, needs `{dplyr}`) → split per
  group; returns a list of dfs.

`extract_intervals()` on list input flattens results (§3.7).
`plot.mnirs()` takes a list of dfs → faceted panels.

### Per-channel arguments and grouping

The same functions accept per-channel args as a named `list()`
keyed by channel name:

- `span = list(smo2 = 10)` — only `smo2` = 10; others fall back to default (`NULL`).
- `span = list(hhb = 10, 5)` — one unnamed element is the fallback
  (`hhb` = 10, others = 5).
- Unrecognised names warned and ignored.

`shift_mnirs()`/`rescale_mnirs()` group channels via `group_channels`;
`extract_intervals()` groups intervals via `group_intervals` and per-group
channels via `group_channels` (§3.5, §3.7).

---

## 3. Function Reference

Unless noted, `nirs_channels`, `time_channel`, `event_channel` default to `NULL` → retrieved from metadata. 
Functions in §2's list accept list/grouped-df input and per-channel `list()` args.

### 3.1 Read

```r
read_mnirs(
    file_path,
    nirs_channels = NULL,   # char vec; (required) rename: c(new = "old")
    time_channel  = NULL,   # char(1); (required) rename: c(time = "Timestamp")
    event_channel = NULL,   # char(1); optional; rename
    sample_rate   = NULL,   # numeric(1) Hz; estimated if NULL
    add_timestamp = FALSE,  # add POSIXct "timestamp" column
    zero_time     = FALSE,  # rebase time[1] to 0
    keep_all      = FALSE,  # keep all columns, else only specified channels
    verbose       = TRUE
)
```

- Auto-detects device, header row, channel names, time column when `NULL`.
- Date-time `time_channel` → numeric, rebased to 0.
- Warns on irregular sampling (non-monotonic, repeated, unequal).

```r
example_mnirs(file = NULL)
## files: "moxy_ramp", "moxy_intervals", "train.red", "artinis",
##        "portamon", others in inst/extdata

create_mnirs_data(data, ...)  # low-level constructor; wraps df as "mnirs"
## ... = named metadata (nirs_channels, time_channel, sample_rate, ...)
```

---

### 3.2 Resample

```r
resample_mnirs(
    data, time_channel  = NULL,
    sample_rate   = NULL,         # source Hz; estimated if NULL
    resample_rate = sample_rate,  # target Hz; default = regularise to source
    method = c("none", "linear", "locf"),
    verbose = TRUE
)
```

- `"none"`: nearest-match to original; new samples → `NA`
- `"linear"`/`"locf"`: interpolate numeric via `stats::approx()`
- Non-numeric cols always either `"locf"` (up-sample) or first-in-bin (down-sample)
- `resample_rate = sample_rate`: regularise only, no rate change

---

### 3.3 Clean

```r
replace_mnirs(
    data, nirs_channels  = NULL, time_channel   = NULL,
    invalid_values = NULL,   # numeric vec; exact values to replace
    invalid_above  = NULL,   # numeric(1); replace x >= threshold
    invalid_below  = NULL,   # numeric(1); replace x <= threshold
    outlier_cutoff = NULL,   # numeric(1); Hampel MAD multiplier; NULL = skip
    width          = NULL,   # integer; rolling window in samples
    span           = NULL,   # numeric; rolling window in time units
    method = c("linear", "median", "locf", "none"),
    verbose = TRUE
)
```

Processing order: invalid → outliers → missing (NA).
`outlier_cutoff`: `3` = Pearson 3-sigma; `2` ≈ Tukey 1.5·IQR; `0` = median filter.

Vector-level:
```r
replace_invalid(x, t, invalid_values, invalid_above, invalid_below,
                width, span, method = c("median", "none"), ...)
replace_outliers(x, t, outlier_cutoff = 3, width, span,
                 method = c("median", "none"), ...)
replace_missing(x, t, width, span,
                method = c("linear", "median", "locf"), ...)
```

`replace_invalid()`/`replace_outliers()` default `"median"`;
`replace_mnirs()`/`replace_missing()` default `"linear"`.

---

### 3.4 Filter

```r
filter_mnirs(
    data, nirs_channels = NULL, time_channel = NULL,
    method = c("smooth_spline", "butterworth", "moving_average"),
    na.rm = FALSE, verbose = TRUE, ...,
    ## method-specific args:
    spar, order, W, fc, sample_rate, type, edges, width, span, partial
)
```

**Aliases:** `"smooth_spline"` = `"spline"`, `"smooth spline"`,
`"smooth-spline"`; `"butterworth"` = `"butter"`; `"moving_average"` =
`"ma"`, `"moving average"`, `"moving-average"`.

**Method args:**
- **`"smooth_spline"`** (`stats::smooth.spline()`): `spar` (NULL = GCV auto); errors on NA/duplicated time.
- **`"butterworth"`** (needs `{signal}`): `order` (default `2L`), `W`/`fc` (normalised or Hz), `sample_rate`, `type` (`"low"`/`"high"`/`"stop"`/`"pass"`), `edges` (`"rev"`/`"rep1"`/`"none"`); errors on NA.
- **`"moving_average"`**: `width` (samples) or `span` (time units); `partial` (default `FALSE`).

Vector-level:
```r
filter_moving_average(x, t, width, span, partial = FALSE, na.rm = FALSE, ...)
filter_ma(...)     # alias for filter_moving_average()
filter_butterworth(x, order = 2L, W, type = "low", edges = "rev", na.rm = FALSE)
```

---

### 3.5 Transform

```r
shift_mnirs(
    data, nirs_channels = NULL, time_channel = NULL,
    group_channels = c("ensemble", "distinct"),
    to = NULL,    # numeric(1); target level (overrides `by`)
    by = NULL,    # numeric(1); shift amount
    width = NULL, # integer; window in samples
    span = NULL,  # numeric; window in time units
    position = c("min", "max", "first"),
    verbose = TRUE
)

rescale_mnirs(
    data, nirs_channels = NULL,
    group_channels = c("ensemble", "distinct"),
    range,        # numeric(2): c(min, max)
    verbose = TRUE
)
```

`shift_mnirs()` moves values up/down, preserves absolute amplitude; 
`rescale_mnirs()` expands/contracts range. 

**`group_channels`:**

| Syntax | Behaviour |
|---|---|
| `"ensemble"` (default) | all channels share one reference (relative scaling preserved) |
| `"distinct"` | each channel independent |
| `list(c("A", "B"), c("C"))` | A+B share reference; C independent |

Groups can be named (`list(smo2 = c("A", "B"))`); names key per-group args.

---

### 3.6 Correct Blood Volume

```r
correct_blood_volume(
    data,
    oxy_channel   = NULL,   # O2Hb/oxy[haem] column name
    deoxy_channel = NULL,   # HHb/deoxy[haem] column name
    total_channel = NULL,   # THb/total[haem] column name (blood-volume proxy)
    verbose = TRUE
)
```

Normalises for blood-volume changes (Beever & Tripp et al, 2020) and 
accommodates negative values (using `shift_mnirs` internals).
Requires ≥2 of 3 channels; third derived (`total = oxy + deoxy`, etc.).
Only specified channels corrected; Names case-sensitive, must match exactly.
`total[haem]` → 0 definitionally after.

---

### 3.7 Extract Intervals

```r
extract_intervals(
    data,                   # "mnirs" df OR list of "mnirs" dfs
    nirs_channels = NULL, time_channel = NULL, event_channel = NULL,
    sample_rate = NULL,
    group_intervals = c("distinct", "ensemble"),
    group_channels = NULL,  # per-group channel selection (see below)
    start = NULL,  # by_time(numeric)/by_label(char)/by_lap(int)/by_sample(int)
    end   = NULL,  # same; NULL = window (span) around start
    span  = list(c(-60, 60)),  # boundaries c(before, after) per interval
    zero_time = FALSE,         # rebase time to 0 per interval
    verbose = TRUE,
    event_groups = deprecated()  # renamed → group_intervals (0.7.0)
)
## returns named list of "mnirs" dfs
```

**List input:** when `data` is a list of "mnirs" dfs, results are flattened 
into a single-layer list. Interval names become `interval_<df>.<interval>`; 
other names (`"ensemble"`, custom group names) suffixed `<name>_<df>`.

**`group_channels`** — selects which `nirs_channels` are ensemble-averaged
per interval group; only relevant when `group_intervals` ensemble-averages.

**Boundary helpers:**

```r
by_time(...)   # numeric time values
by_label(..., ignore_case = FALSE, fixed = FALSE)  # event labels; regex by default
by_lap(...)    # integer lap numbers (start = first sample, end = last)
by_sample(...) # integer row indices
```

- Raw coercion: numeric → `by_time()`; character → `by_label()`; explicit integer (`2L`) → `by_lap()`
- `start`/`end` may mix types; shorter recycled. `by_label()`/`by_lap()` need `event_channel`
- `span = c(before, after)`: negative shifts start earlier, positive shifts end later; single value recycled by sign (`60` → `c(0, 60)`, `-60` → `c(-60, 0)`)
- `span`/`group_channels` as `list()` recycled per interval group

**`group_intervals`:**

| Syntax | Behaviour |
|---|---|
| `"distinct"` (default) | one df per detected interval |
| `"ensemble"` | single ensemble-averaged df (needs regularised time grid) |
| `list(c(1, 2), c(3, 4))` | one ensemble df per group; named lists pass names through |

---

### 3.8 Analyse Kinetics

```r
analyse_kinetics(
    data,           # "mnirs" df | named list of "mnirs" dfs | grouped df
    nirs_channels = NULL, time_channel = NULL,
    method = c("response_time", "peak_slope", "monoexponential", "sigmoidal"),
    start_time = NULL,  # fit onset (t = 0); NULL = interval_times metadata, else t[1] else 0
    direction  = c("auto", "positive", "negative"),
    end_window = Inf,   # truncate fit after extreme; Inf = global extreme
    verbose = TRUE, ...,
    ## method-specific args:
    fraction, width, span, align, partial, na.rm, use_TD, shape
)
## analyze_kinetics(...)  # US-spelling alias
```

`direction` also constrains fitted-amplitude sign for
`"monoexponential"`/`"sigmoidal"`.

**Method aliases:**

| Canonical | Aliases |
|---|---|
| `"response_time"` | `"half response time"`, `"recovery time"`, `"half recovery time"`, `"half time"`, `"HRT"` |
| `"peak_slope"` | `"slope"` |
| `"monoexponential"` | `"monoexp"`, `"exp"`, `"exponential"`, `"MRT"`, `"tau"` |
| `"sigmoidal"` | `"sigmoid"`, `"logistic"`, `"gompertz"`, `"xmid"` |

**Per-method args (`...`):**
- **`"response_time"`**: `fraction` (default `0.5`; `0.632` ≈ MRT).
- **`"peak_slope"`**: `width` or `span` (one required); `align` (`"centre"`/`"left"`/`"right"`); `partial`, `na.rm` (default `FALSE`).
- **`"monoexponential"`**: `use_TD` (default `TRUE`; 4-param → 3-param fallback).
- **`"sigmoidal"`**: `shape` (`"symmetric"` default = `SSlogistic()`; `"gompertz"` early-inflection; `"gompertz_left"` late-inflection).

Per-channel overrides via inline named `list()` (names must match
`nirs_channels`):
```r
analyse_kinetics(data, nirs_channels = c(hhb, smo2), method = "peak_slope",
    span = list(smo2 = 10), direction = list(hhb = "negative"))
```

#### `"mnirs_kinetics"` return — formatted table & list internal components

| Element | Type | Description |
|---|---|---|
| `method` | character | method used |
| `model` | named list | per-interval per-channel model objects (`lm`/`nls`/`NULL`) |
| `coefficients` | tibble | one row per channel per interval |
| `data` | named list | input dfs with `<channel>_fitted` cols |
| `interval_times` | data frame | `interval`, `start_times` (+ `end_times` when present) |
| `diagnostics` | data frame | fit diagnostics |
| `channel_args` | data frame | resolved args per channel per interval |
| `call` | call | matched function call |

**Coefficients** (prefixed `interval`, `nirs_channels`, `time_channel`):

| Method | Columns |
|---|---|
| `"response_time"` | `A` baseline mean, `B` extreme, `response_time` (elapsed from `start_time`), `response_value` (observed at response), `fitted` (target `A + (B-A)*fraction`), `idx` (response row) |
| `"peak_slope"` | `slope` (units `x/t`), `intercept`, `fitted` (predicted value at peak), `peak_slope_time` (elapsed from `start_time`), `idx` (window-centre row) |
| `"monoexponential"` | `A` baseline, `B` asymptote, `tau` time constant, `k` rate constant (`1/tau`), `TD` time delay, `MRT` mean response time (`TD+tau`), `HRT` half-response time (`TD+tau·ln2`), `tau_fitted`/`MRT_fitted`/`HRT_fitted` (predicted value at each) |
| `"sigmoidal"` | `A` start asymptote, `B` end asymptote, `xmid` inflection time, `slope` (`dx/dt` at `xmid`), `xmid_fitted` (predicted value at `xmid`) |

**Diagnostics:** `n_obs`, `r2`, `adj_r2`, `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`, `bic`.

**Vector-level:**
```r
response_time(x, t = seq_along(x), start_time = 0, fraction = 0.5,
    direction = c("auto", "positive", "negative"), verbose = TRUE)
## → A, B, response_time, response_value, fitted,
##   baseline_idx, response_idx, extreme_idx

peak_slope(x, t = seq_along(x), width = NULL, span = NULL,
    align = c("centre", "left", "right"),
    direction = c("auto", "positive", "negative"),
    partial = FALSE, na.rm = FALSE, verbose = TRUE)
## → slope, intercept, y, t, idx, fitted, window_idx, model

monoexponential(t, A, B, tau, TD = NULL)
## 3-param: A + (B - A) * (1 - exp(-t / tau))
## 4-param: ifelse(t <= TD, A, A + (B - A) * (1 - exp(-(t - TD) / tau)))
nls(x ~ SSmonoexp(t, A, B, tau, TD), data = df)   # 3- or 4-param

logistic(t, A, B, xmid, slope, asym = NULL)  # 4-param symmetric / 5-param Richards
gompertz(t, A, B, xmid, slope)               # early-inflection
gompertz_left(t, A, B, xmid, slope)          # late-inflection
nls(x ~ SSlogistic(t, A, B, xmid, slope), data = df)      # 4- or 5-param (fragile)
nls(x ~ SSgompertz(t, A, B, xmid, slope), data = df)
nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = df)
```

---

### 3.9 Plot and Print

```r
plot.mnirs(x, points = FALSE, time_labels = FALSE, na.omit = FALSE, ...)
## needs {ggplot2}; via plot(data); returns ggplot2 object
## x = "mnirs" df or list of dfs (list → faceted panels)
## na.omit = TRUE drops NA/non-finite; ... = facet_wrap args, n.breaks, breaks

plot.mnirs_kinetics(x, fitted = TRUE, markers = TRUE, labels = TRUE, ...)
## needs {ggplot2}; via plot(result); returns ggplot2 object
## observed signal per channel, faceted by interval (builds on plot.mnirs)
## fitted  = dashed fitted curve (parametric methods; not "response_time")
## markers = dotted onset line + key coefficient point(s)
## labels  = per-panel annotation of key coefficient(s)
## ...     = label_size, others passed to plot.mnirs (points, time_labels, nrow, ncol, scales)

print(result)  # "mnirs_kinetics"; formatted coefficient table (max 10 rows)
print(data)    # "mnirs"; strips class, prints tibble

theme_mnirs(base_size = 14, base_family = "sans",
            border = c("partial", "full"),
            ink = "black", paper = "white", accent = "#0080ff", ...)

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
| `stats` | models, interpolation | Import | |
| `tibble` | `"mnirs"` class | Import | |
| `tidyselect` | column selection | Import | |
| `signal` | Butterworth | Suggests | `"butterworth"` method |
| `ggplot2` | plotting | Suggests | `plot.mnirs()`, theme/scales |
| `scales` | axis formatting | Suggests | `plot.mnirs()` |
| `dplyr` | grouped df input | Suggests | grouped-df dispatch |
| `knitr`, `quarto` | vignettes | Suggests | |
| `zoo`, `testthat` | testing | Suggests | |

---

## 5. Constraints

| Constraint | Detail |
|---|---|
| Irregular samples warning from `read_mnirs()` | fires in pipe before downstream `resample_mnirs()`; verify output |
| Pipeline order | `resample → replace → filter`; wrong order = wrong results |
| `"smooth_spline"`/`"butterworth"` fail on NA | `replace_mnirs()` first or `na.rm = TRUE` |
| `"smooth_spline"` fails on duplicated time | `resample_mnirs()` first |
| Ensemble needs regularised samples | `group_intervals = "ensemble"` warns if irregular |
| Per-channel `list()` args | unrecognised list names warned and ignored |
| `monoexponential` fallback | 4-param → 3-param; `NA` coefficients on convergence failure |
| `direction`-bounded fit | `"monoexponential"`/`"sigmoidal"` return `NA` coefficients if direction unsatisfiable |

---

## 6. Key Source Files

| File | Contents |
|---|---|
| `R/read_mnirs.R` | `read_mnirs()`, `example_mnirs()`, `create_mnirs_data()` |
| `R/resample_mnirs.R` | `resample_mnirs()` |
| `R/replace_mnirs.R` | `replace_mnirs()`, `replace_invalid/outliers/missing()` |
| `R/filter_mnirs.R` | `filter_mnirs()`, `filter_moving_average()`, `filter_ma()`, `filter_butterworth()` |
| `R/shift_mnirs.R` | `shift_mnirs()` |
| `R/rescale_mnirs.R` | `rescale_mnirs()` |
| `R/correct_blood_volume.R` | `correct_blood_volume()` |
| `R/extract_intervals.R` | `extract_intervals()` |
| `R/extract_interval_helpers.R` | `by_time/label/lap/sample()`, boundary resolution |
| `R/analyse_kinetics.R` | `analyse_kinetics()`/`analyze_kinetics()` + S3 dispatch |
| `R/analyse_kinetics_helpers.R` | channel/interval orchestration, `compute_diagnostics()` |
| `R/analyse_peak_slope.R` | `peak_slope()`, `slope()`, `rolling_slope()` |
| `R/analyse_monoexponential.R` | `monoexponential()`, `SSmonoexp()` |
| `R/analyse_sigmoidal.R` | `logistic()`, `gompertz()`, `gompertz_left()`, `SSlogistic/gompertz/gompertz_left()` |
| `R/analyse_response_time.R` | `response_time()` |
| `R/plot.mnirs.R` | `plot.mnirs()`, `plot.mnirs_kinetics()`, `kinetics_annotations()`, `as_plot_data()`, `theme_mnirs()`, `palette_mnirs()`, scale/format fns |
| `R/mnirs_methods.R` | `print.mnirs()`, `print.mnirs_kinetics()` |
| `R/channel_args.R` | `resolve_channel_args()` — per-channel/group arg broadcast |
| `R/as_data_list.R` | `map_mnirs_intervals()`, `as_data_list()` — list/grouped dispatch |
| `R/validate_mnirs.R` | input validation |
| `R/data.R` | example file descriptions |

---

## 7. Developer pointers

| Workflow | Location |
|---|---|
| Input validation | `R/validate_mnirs.R` (`validate_nirs_channels/time_channel/sample_rate/width_span/x_t/numeric/group_channels()`) |
| Per-channel arg broadcast | `resolve_channel_args()`, `R/channel_args.R` |
| List/grouped-df dispatch | `map_mnirs_intervals()`, `as_data_list()`, `R/as_data_list.R` |
| Kinetics orchestration + diagnostics | `R/analyse_kinetics_helpers.R` |
| Interval boundary resolution | `R/extract_interval_helpers.R` |
| Read/device detection | `R/read_mnirs_helpers.R` |
| User-facing messages | `cli_abort()`/`cli_warn()`/`cli_inform()` |
| Roxygen2 with markdown | pkgdown config in `_pkgdown.yml` |
