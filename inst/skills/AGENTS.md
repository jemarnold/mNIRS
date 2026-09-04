# `{mnirs}` Agent Reference

**v0.8.0 | R (>= 4.1) | MIT** — workflow/dependency map for muscle
near-infrared spectroscopy (mNIRS) processing & analysis package.
Website: https://jemarnold.github.io/mnirs/

Read `.xls(x)`/`.csv`/`.txt`/`.tsv`/`.ftn(2)` exports → resample → clean → filter → (transform)
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

`verbose` logical to display/suppress warning messages, controlled globally `getOption("mnirs.verbose", TRUE)` or passed explicitly.

---

## 2. Pipeline

```
read_mnirs()
├─ plot()                           # visualise "mnirs" data at each step
└─ resample_mnirs()                 # up/down-sample, regularise time grid
   └─ replace_mnirs()               # clean invalid/outliers/NA
      └─ filter_mnirs()             # smooth
         ├─ shift_mnirs()           # shift baseline, preserve amplitude
         ├─ rescale_mnirs()         # normalise range, modify amplitude
         ├─ correct_blood_volume()  # optional: normalise changes in THb
         └─ extract_intervals()     # detect & extract list of interval dfs
            └─ analyse_kinetics()   # compute response rate & time course
               └─ plot()            # plot "mnirs_kinetics": fit + markers
```

Failure modes + fixes in §5.

---

## 3. Function Reference

Data frame-level functions accept `data` as:

- single df in → single df out
- list of dfs → list, each processed separately
- grouped df (`dplyr::group_by()`, needs `{dplyr}`) → list, one per group

`extract_intervals()` flattens list input (§3.7). `plot.mnirs()` returns facet per-df.

`nirs_channels`, `time_channel`, `event_channel` and other args default to `NULL` → retrieved from metadata.

Most functions accept `{tidyselect}` channel names; `nirs_channels = "smo2"` or `smo2`; `time_channel = starts_with("time")`.

Functions also accept per-channel args as named `list()` keyed by channel:

- `span = list(smo2 = 10)` — only `smo2` = 10; others fall back to default (`NULL`)
- `span = list(hhb = 10, 5)` — one unnamed element is the fallback (`hhb` = 10, others = 5)
- unrecognised names warned and ignored

Channel/interval grouping: `group_channels` (§3.5), `group_intervals` (§3.7).


### 3.1 Read

```r
read_mnirs(
    file_path,
    nirs_channels = NULL,   # char vec; (required) rename: c(new = "old")
    time_channel  = NULL,   # char(1); (required) rename: c(time = "Timestamp")
    event_channel = NULL,   # char(1); optional; rename: c(lap = "Lap/Event")
    sample_rate   = NULL,   # numeric(1) Hz; estimated if NULL
    add_timestamp = FALSE,  # add POSIXct "timestamp" column
    zero_time     = FALSE,  # rebase time[1] to 0
    keep_all      = FALSE,  # keep all columns, else only specified channels
    verbose       = TRUE    # all functions have default verbose = TRUE
)
```

- Auto-detects device, header row, channel names, time column when `NULL`.
- Artinis Oxysoft: channels auto-named from file "Legend" metadata — traces cleaned & lowercased (`rx1_tx1_o2hb`). "(Sample number)" → `sample` (+ derived `time`); "(Event)" → `event`; trailing unnumbered label col → `labels` (select with `event_channel = c(event = "labels")`).
- PIONIRS `.ftn(2)`: `Time`, `TagLabel` event, `StO2*` channels; `Iteration`/`Tag` companion cols.
- Device companion cols (`labels`, `Iteration`, `Tag`, other legend metadata) returned only with `keep_all = TRUE`.
- Date-time `time_channel` converted → numeric, rebased to 0.
- Warns on irregular sampling (non-monotonic, repeated, unequal).
- Example files retrieved from `example_mnirs()`.
- Add class *"mnirs"* & metadata to df with `create_mnirs_data()`.

```r
example_mnirs(file = NULL)  # NULL = list all; partial matching
## "artinis_intervals", "moxy_intervals", "moxy_ramp",
## "portamon-oxcap", "train.red_intervals"

create_mnirs_data(data, ...)  # low-level constructor; wraps df as "mnirs"
## ... = named metadata (nirs_channels, time_channel, sample_rate, ...)
```

---

### 3.2 Resample

```r
resample_mnirs(
    data,
    time_channel  = NULL,
    sample_rate   = NULL,         # source Hz; estimated if NULL
    resample_rate = sample_rate,  # target Hz; default = regularise to source
    method = c("none", "linear", "locf")
)
```

- `"none"`: nearest-match to original; new samples → `NA`
- `"linear"`/`"locf"`: interpolate numeric via `stats::approx()`
- Non-numeric cols always either `"locf"` (up-sample) or first-in-bin (down-sample)
- Default `resample_rate = sample_rate`: regularise only, no rate change

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
    method = c("linear", "median", "locf", "none")
)
```

- Processing order: invalid → outliers → missing (NA).
- Recommended `outlier_cutoff`: `3` = Pearson 3-sigma; `2` ≈ Tukey 1.5·IQR.
- `width` XOR `span`; `width` takes precedence.

Vector-level:
```r
## default method = "median"
replace_invalid(x, t, invalid_values, invalid_above, invalid_below,
                width, span, method = c("median", "none"), ...)
replace_outliers(x, t, outlier_cutoff = 3, width, span,
                 method = c("median", "none"), ...)
## default method = "linear", as replace_mnirs()
replace_missing(x, t, width, span, method = c("linear", "median", "locf"), ...)
```

---

### 3.4 Filter

```r
filter_mnirs(
    data, nirs_channels = NULL, time_channel = NULL,
    method = c("smooth_spline", "butterworth", "moving_average"),
    na.rm = FALSE, ...,
    ## method-specific args:
    spar, order, W, fc, sample_rate, type, edges, width, span, partial
)
```

**Method aliases**: matching is case- and separator-insensitive; `<space>`, `-`, `_`.

**Method args:**
- **`"smooth_spline"`** (`stats::smooth.spline()`): `spar` (NULL = GCV auto); errors on NA/duplicated time.
- **`"butterworth"`** (needs `{signal}`): filter `order` (default `2L`), cutoff frequency(ies) `W`/`fc` (normalised or Hz), `sample_rate`, filter `type` (`"low"`/`"high"`/`"stop"`/`"pass"`), `edges` (`"rev"`/`"rep1"`/`"none"`); errors on NA.
- **`"moving_average"`**: rolling window in `width` (samples) XOR `span` (time units); `partial` (default `FALSE`).

Vector-level (can call `stats::smooth.spline()` directly):
```r
filter_moving_average(x, t, width, span, partial = FALSE, na.rm = FALSE, ...)
filter_butterworth(x, order = 2L, W, type = "low", edges = "rev", na.rm = FALSE)
filter_ma(...) & filter_butter(...) # function aliases
```

---

### 3.5 Transform

```r
## moves values up/down, preserves absolute amplitude
shift_mnirs(
    data, nirs_channels = NULL, time_channel = NULL,
    group_channels = c("ensemble", "distinct"), # or as custom `list()`
    to = NULL,    # numeric(1); target level (overrides `by`)
    by = NULL,    # numeric(1); shift amount, +/-ve
    width = NULL, # integer; window in samples
    span = NULL,  # numeric; window in time units
    position = c("min", "max", "first") # which reference value
)

## expands/contracts range
rescale_mnirs(
    data, nirs_channels = NULL,
    group_channels = c("ensemble", "distinct"),
    range         # numeric(2): c(min, max)
)
```

**`group_channels`:**

| Syntax | Behaviour |
|---|---|
| `"ensemble"` (default) | all channels share one reference (relative scaling preserved) |
| `"distinct"` | each channel independent (relative scaling lost) |
| `list(c("A", "B"), c("C"))` | A+B share reference; C independent |

Custom groups nameable (`list(smo2 = c("A", "B"))`); per-channel args can be keyed per-group instead.

---

### 3.6 Correct Blood Volume

```r
correct_blood_volume(
    data,
    oxy_channel   = NULL,  # O2Hb/oxy[haem] column name(s); vectors paired by position
    deoxy_channel = NULL,  # HHb/deoxy[haem] column name(s)
    total_channel = NULL   # THb/total[haem] column name(s); derived from pair if NULL
)
```

Multiple channel pairs: pass equal-length vectors with paired elements

Normalises blood-volume changes (Ryan et al 2012; Beever & Tripp et al, 2020). Handles negative values via `shift_mnirs`. Returns `total[haem]` → 0 definitionally.

---

### 3.7 Extract Intervals

```r
extract_intervals(
    data, nirs_channels = NULL, time_channel = NULL, event_channel = NULL,
    sample_rate = NULL,
    group_intervals = c("distinct", "ensemble"),
    group_channels = NULL,  # per-group channel selection (see below)
    start = NULL,  # by_time(numeric)/by_label(char)/by_lap(int)/by_sample(int)
    end   = NULL,  # same; NULL = window (span) around start
    span  = list(c(-60, 60)),    # boundaries c(start, end) per interval
    zero_time = FALSE,           # rebase time to 0 per interval
)
## returns named list of "mnirs" dfs
```

Returns flattened single-layer list of named intervals; `interval_<df>.<interval>`, or from `"ensemble"`, custom-grouped suffixed `<name>_<df>`.

**`group_channels`** (as above) — which `nirs_channels` are ensemble-averaged per interval group; only relevant when `group_intervals` ensemble-averages.

**Boundary helpers:**

```r
by_time(...)   # numeric time values
by_label(..., ignore_case = FALSE, fixed = FALSE)  # event labels; regex by default
by_lap(...)    # integer lap numbers (start = first sample, end = last)
by_sample(...) # integer row indices
```

- `by_label()`/`by_lap()` need `event_channel`.
- `start`/`end` may accept `list()` of multiple mixed `by_*()` types (e.g. `start = list(by_time(30), by_label("go"))`); resolved times concatenated in supplied order.
- `span = c(start, end)`: negative expands bound earlier before start, positive expands bound later after end; single value recycled by sign (`60` → `c(0, 60)`, `-60` → `c(-60, 0)`).
- `span`/`group_channels` as `list()` per interval group; recycled as needed.


**`group_intervals`:**

| Syntax | Behaviour |
|---|---|
| `"distinct"` (default) | one df per detected interval |
| `"ensemble"` | single ensemble-averaged df (needs regularised time grid) |
| `list(group1 = c(1, 2), c(3, 4))` | one ensemble df per group; names passed through |

---

### 3.8 Analyse Kinetics

```r
analyse_kinetics(
    data, nirs_channels = NULL, time_channel = NULL,
    method = c("response_time", "peak_slope", "monoexponential",
               "biexponential", "sigmoidal"),
    start_time = NULL,  # fit onset (t = 0); NULL = interval_times metadata, else t[1] else 0
    direction  = c("auto", "positive", "negative"),
    group_intervals = "ensemble",  # or list() of row numbers (see below)
    zero_time = FALSE,  # rebase time to 0 per interval/group; shifts interval_times
    end_window = Inf,   # truncate fit after first extreme; Inf = global extreme
    ...,
    ## method-specific (explicit formals, not via `...`, see below):
    fraction = 0.5, width = NULL, span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE, na.rm = FALSE,
    use_TD = TRUE, shape = c("symmetric", "gompertz", "gompertz_left"),
    tau_mult = 3, fix = NULL
    ## via `...`: control = NULL (nls methods; see below)
)
## analyze_kinetics(...) alias
```

**Method aliases**: matching is case- and separator-insensitive; `<space>`, `-`, `_`. Accepts common shorthand: `hrt`, `slope`, `mrt`, `tau`, `gompertz`, `xmid`, ...

`direction` also constrains fitted-amplitude sign for `"monoexponential"`/
`"biexponential"`/`"sigmoidal"`, resolved per channel, returns `NA` if unsatisfiable.

**Per-method args:**
- **`"response_time"`**: `fraction` (default `0.5`; `0.632` ≈ MRT; vectorised, e.g. `c(0.5, 0.632)` → one coefficient row per fraction).
- **`"peak_slope"`**: `width` XOR `span`; `align` (`"centre"`/`"left"`/`"right"`); `partial`, `na.rm` (default `FALSE`).
- **`"monoexponential"`**: `use_TD` (default `TRUE`; 4-param → 3-param fallback), `fix`.
- **`"biexponential"`**: `use_TD` (default `TRUE`; 6-param → 5-param fallback), `fix`. Sequential fit: fast monoexp on `end_window` window (`Inf` → first extreme + 20 time units) → full biexp with `A`/`tau`/`TD` held near stage-1 values, `B`/`B2`/`tau2` free. Falls back (warning; `model` column) to exp_drift → monoexp on fit failure, monotonic `texc`, `tau2 >= 2 × span`, or `|B2 - B| < 2 × rmse`. Coef columns = union of the chain (`NA` where n/a). Undocumented `model_fallback = FALSE` keeps raw fit.
- **`"exponential_drift"`**: `use_TD`, `tau_mult` (default `3`; drift onset `TD + tau_mult * tau`, always held constant), `fix`. `texc` = takeover point `max(onset, TD + tau × log(|B - A| / (|slope| × tau)))` (turning point when phases oppose). Falls back to monoexp on fit failure or `|slope| × (t_end - onset) < 2 × rmse` (`model` column; `model_fallback = FALSE` keeps raw fit).
- **`"sigmoidal"`**: `shape` (`"symmetric"` default = `SSlogistic()`; `"gompertz"` early-inflection (right); `"gompertz_left"` late-inflection), `fix`.
- **All nls methods**: `control` via `...` (`list()` or `nls.control()`, e.g. `list(maxiter = 200)`) merged over internal defaults (`maxiter = 500, warnOnly = TRUE` on `"port"` fits) at every `nls()` call incl. direction refits and fallbacks. Global only (not per-channel/interval). Unknown names abort.

Per-channel overrides via inline named `list()` (names must match `nirs_channels`):
```r
analyse_kinetics(data, nirs_channels = c(hhb, smo2), method = "peak_slope",
    span = list(smo2 = 10), direction = list(hhb = "negative"))
```

**`group_intervals`** (row grouping; no `"distinct"`):

| Syntax | Behaviour |
|---|---|
| `"ensemble"` (default) | all rows of each df analysed together |
| `list(trial1 = 1:12, trial2 = 13:24)` | each row group a separate interval; unnamed -> `interval_<n>`; multi-df input suffixed `<group>_<df>` |

Rows in no group dropped (message); overlapping rows warned. Per-interval args key by group names.
Main use: recursive analysis of `"mnirs_kinetics"` coefs across trials, e.g.
`analyse_kinetics(result, nirs_channels = slope, time_channel = peak_slope_time, method = "monoexp", group_intervals = list(trial1 = 1:12, trial2 = 13:24))`.

**`fix`** — hold params constant:

```r
fix = list(A = 0)                                  # fixed for all channels
fix = list(smo2 = list(A = 0))                     # per channel
fix = list(interval_1 = list(smo2 = list(A = 0)))  # per interval x channel
```

- Named list of finite numeric scalars; names must be model params; cannot fix all.
- Fixed params excluded from estimation and from `n_params`; reported at fixed value.
- `"monoexponential"`: `TD` fixable only when `use_TD = TRUE`, and doing so disables the 3-param fallback. Prefer `use_TD = FALSE` over `fix = list(TD = 0)`.

#### `"mnirs_kinetics"` return — formatted table & list internal components

- one row per channel per interval.

| Element | Type | Description |
|---|---|---|
| `method` | character | canonical method used |
| `model` | named list | per-interval per-channel model objects (`lm`/`nls`) |
| `coefficients` | data frame | model & derived coefficients |
| `data` | named list | input dfs with predicted `<channel>_fitted` cols |
| `interval_times` | data frame | `interval`, `start_times` (+ `end_times` when present) |
| `diagnostics` | data frame | fit diagnostics (AIC, BIC, RMSE, R^2, etc) |
| `channel_args` | data frame | resolved args |
| `warnings` | data frame | `type` = `"warning"`/`"error"` messages; empty if none. Captured regardless of `verbose` (suppresses console output only) |
| `call` | call | matched function call |

**Coefficients** (prefixed `interval`, `nirs_channels`, `start_time`):

Times are elapsed from `start_time`; `*_fitted` = predicted value at that point.

| Method | Columns |
|---|---|
| `"response_time"` | `fraction` (one row per value), `A` baseline mean, `B` extreme (peak/trough) value, `response_time`, `response_value` (observed), `fitted` (target `A + (B-A)*fraction`), `idx` (sample/row number at `response_value`) |
| `"peak_slope"` | `slope` (`x/t`), `intercept`, `fitted`, `peak_slope_time`, `idx` (sample/row number at `align` position) |
| `"monoexponential"` | `A` baseline, `B` asymptote, `tau`, `k` (`1/tau`), `TD` delay (if `use_TD`), `MRT` (`TD+tau`), `HRT` (`TD+tau·ln2`), `MRT_fitted`, `HRT_fitted` |
| `"biexponential"` | `A` start, `B` & `tau` fast component, `MRT` (`TD+tau`), `texc` (fitted excursion point; `NA` if monotonic), `B2` & `tau2` slow component, `TD` delay (if `use_TD`), `MRT_fitted`, `texc_fitted`; plus `model` and the exp_drift/monoexp columns (`NA` unless fallen back) |
| `"exponential_drift"` | monoexp columns + `texc` (excursion point where drift rate overtakes primary rate; ≥ drift onset `TD + tau_mult·tau`), `slope` (`dx/dt`), `tau_mult`, `texc_fitted`; plus `model` |
| `"sigmoidal"` | `A` & `B` start + end asymptotes, `xmid` inflection time (only literally *"middle"* for `shape = "symmetric"`), `slope` (`dx/dt` at `xmid`), `xmid_fitted` |

**Diagnostics:** `n_obs`, `n_params`, `r2`, `adj_r2`, `rmse`, `snr`, `cv_rmse`,
`aic`, `aicc`, `bic`. `n_params` = free params estimated, excluding `fix`, so
reduced-parameter fallback fits are distinguishable. `aic`/`bic` comparable only
within matching `n_obs` + `n_params`.

**Vector-level:**
```r
response_time(x, t = seq_along(x), start_time = 0, fraction = 0.5,
    direction = c("auto", "positive", "negative"))
## → A, B, response_time, response_value, fitted,
##   baseline_idx, response_idx, extreme_idx

peak_slope(x, t = seq_along(x), width = NULL, span = NULL,
    align = c("centre", "left", "right"),
    direction = c("auto", "positive", "negative"),
    partial = FALSE, na.rm = FALSE)
## → slope, intercept, y, t, idx, fitted, window_idx, model

monoexponential(t, A, B, tau, TD = NULL)
## 3-param: A + (B - A) * (1 - exp(-t / tau))
## 4-param: A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau))
nls(x ~ SSmonoexponential(t, A, B, tau, TD), data = df)   # 3- or 4-param

biexponential(t, A, B, tau, B2, tau2, TD = NULL)
## 5-param: A + (B-A)*(1 - exp(-t/tau)) + (B2-B)*(1 - exp(-t/tau2))
## 6-param: as above with ts = pmax(t - TD, 0) substituted for t
## both phases clocked from onset; excursion point
## texc = TD + log(r)/(1/tau - 1/tau2), r = -(B-A)*tau2/((B2-B)*tau), NA if r <= 1
nls(x ~ SSbiexponential(t, A, B, tau, B2, tau2, TD), data = df,
    algorithm = "port", lower = c(-Inf, -Inf, 0, -Inf, 0, 0))   # 5- or 6-param

logistic(t, A, B, xmid, slope, asym = NULL)  # 4-param symmetric / 5-param Richards
gompertz(t, A, B, xmid, slope)               # early-inflection
gompertz_left(t, A, B, xmid, slope)          # late-inflection
nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = df) # 4- or 5-param (fragile)
nls(x ~ SSgompertz(t, A, B, xmid, slope), data = df)
nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = df)
```

---

### 3.9 Plot

```r
## all plot/theme fns need {ggplot2}; scales/format need {scales}
plot.mnirs(x, points = FALSE, time_labels = FALSE, na.omit = FALSE, ...)
## x = "mnirs" df or list (list → faceted); na.omit = TRUE drops NA/non-finite;
## ... = facet_wrap args, n.breaks, breaks

plot.mnirs_kinetics(x, fitted = TRUE, markers = TRUE, labels = TRUE, ...)
## observed signal per channel, faceted by interval (builds on plot.mnirs)
## fitted = dashed curve (parametric methods; not "response_time");
## markers = dotted onset line + key coefficient point(s);
## labels = per-panel annotation; ... = label_size, passed to plot.mnirs()


theme_mnirs(base_size = 14, base_family = "sans", border = c("partial", "full"),
            ink = "black", paper = "white", accent = "#0080ff", ...)
palette_mnirs(...)      # no args = all 12; palette_mnirs(4); palette_mnirs("red", "blue")
scale_colour_mnirs(...) # alias scale_color_mnirs()
scale_fill_mnirs(...)
breaks_timespan(unit = "secs", n = 5)
format_hmmss(x)         # numeric seconds → "mm:ss" or "h:mm:ss"
```

---

## 4. Dependencies

`Depends: R (>= 4.1)`

| Package | Role | Type | Condition |
|---|---|---|---|
| `cli` | user messages | Import | |
| `data.table` | data manipulation | Import | |
| `lifecycle` | deprecation | Import | |
| `readxl` | XLS/XLSX | Import | |
| `rlang` | NSE/tidy eval | Import | |
| `stats` | models, interpolation | Import | |
| `tibble` | returned `"mnirs"` class df | Import | |
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
| Recommended pipeline order | `resample → replace → filter → ...`; different order = different results |
| `group_intervals = "ensemble"` needs regularised time grid | `resample → extract_intervals`; warns if irregular |
| `monoexponential`/`biexponential` convergence fallback | Weak fits (certain convergence errors) return fit with warnings |
| `direction`-bounded fit | return `NA` coefficients if unsatisfiable; verify `biexponential` in particular |
| `biexponential` identifiability | fast phase fixed by stage-1 monoexp on `end_window`; a window too long past the extreme gives a slow stage 1 and fallback to exp_drift/monoexp; check `warnings` + `model` column |

---

## 6. Source Map

| File | Contents |
|---|---|
| `R/read_mnirs.R` | `read_mnirs()`, `example_mnirs()`, `create_mnirs_data()` |
| `R/read_mnirs_helpers.R` | device/header/channel detection |
| `R/resample_mnirs.R` | `resample_mnirs()` |
| `R/replace_mnirs.R` | `replace_mnirs()`, `replace_invalid/outliers/missing()` |
| `R/filter_mnirs.R` | `filter_mnirs()`, `filter_moving_average()`/`filter_ma()`, `filter_butterworth()`/`filter_butter()` |
| `R/rolling_helpers.R` | internal: `compute_window_bounds()`, `compute_outliers()` |
| `R/shift_mnirs.R` / `R/rescale_mnirs.R` | `shift_mnirs()` / `rescale_mnirs()` |
| `R/correct_blood_volume.R` | `correct_blood_volume()` |
| `R/extract_intervals.R` | `extract_intervals()` |
| `R/extract_interval_helpers.R` | `by_time/label/lap/sample()`; boundary resolution |
| `R/analyse_kinetics.R` | `analyse_kinetics()`/`analyze_kinetics()` + S3 dispatch |
| `R/aanalyse_kinetics_helpers.R` | *(leading `aa` intentional — load order)* `method_aliases`, channel/interval orchestration, `detect_direction()`, `enforce_direction()`, `compute_diagnostics()`, `kinetics_warnings_df()` |
| `R/analyse_response_time.R` | `response_time()` |
| `R/analyse_peak_slope.R` | `peak_slope()`, `rolling_slope()` |
| `R/analyse_monoexponential.R` | `monoexponential()`, `SSmonoexponential()` |
| `R/analyse_biexponential.R` | `biexponential()`, `SSbiexponential()`, `biexp_init()`, `biexp_texc()` |
| `R/analyse_sigmoidal.R` | `logistic()`, `gompertz()`, `gompertz_left()`, `SS*()` |
| `R/plot.mnirs.R` | `plot.mnirs()`, `plot.mnirs_kinetics()`, `theme_mnirs()`, `palette_mnirs()`, scale/format fns |
| `R/mnirs_methods.R` | `print.mnirs()`, `print.mnirs_kinetics()` |
| `R/channel_args.R` | `resolve_channel_args()` — per-channel/group arg broadcast |
| `R/as_data_list.R` | `map_mnirs_intervals()`, `as_data_list()` — list/grouped dispatch |
| `R/validate_mnirs.R` | `validate_numeric/mnirs_data/nirs_channels/time_channel/event_channel/sample_rate/width_span/x_t/start_time/fix/findInt()` |
| `R/signif_trailing.R` | internal: `signif_trailing()`, `seq_range()` |
| `R/data.R` / `R/mnirs-package.R` | example file docs / package-level roxygen |

User-facing messages via `cli_abort()`/`cli_warn()`/`cli_inform()`.
Roxygen2 with markdown enabled; pkgdown config in `_pkgdown.yml`.
