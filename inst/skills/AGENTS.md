# mnirs agent map

Reference for orchestrator/subagents. Package: muscle NIRS file
import, time-series processing, interval extraction, kinetics analysis, plotting.
Target release: 0.7.0. Current authority: source + tests > generated
`man/`/`NAMESPACE` > README/vignettes.

## Non-negotiable workflow

- Read target `R/<stem>.R`, related `*_helpers.R`, matching
  `tests/testthat/test-<stem>.R`, then roxygen source. Trace callers before edit.
- Surgical diff. Preserve API, classes, attributes, conditions, call attribution,
  grouping, tidy evaluation, and optional-dependency behaviour unless task changes it.
- Root cause first. For bug fix: minimal failing test -> fix -> focused test -> wider
  checks. Always obtain maintainer approval before editing tests; bug/feature request
  alone does not authorise test changes.
- Never edit `NAMESPACE`, `man/*.Rd`, or `README.md` directly. Sources:
  roxygen in `R/*.R`; `README.Rmd`; `_pkgdown.yml`; Quarto `.qmd`.
- User-facing conditions use `cli_abort()`, `cli_warn()`, `cli_inform()`; pass correct
  `call`/`env`. Validation belongs in `R/validate_mnirs.R` or shared domain helper.
- Style: base pipe `|>`; `\() ...` only single-line lambdas; explicit `return()`;
  4-space indent; 80 columns (`air.toml`); lower-case `##` comments explain why.
- Prefer vectorised/base/tidyverse code. No `for` loops or nested `if` growth.
- Do not run devtools commands.
- Do not inspect git history. Current tree only. Preserve unrelated dirty changes.
- User-facing change: roxygen/docs + test + concise `NEWS.md` bullet when warranted.
  New public topic also needs `_pkgdown.yml` reference coverage.

## Agent coordination

- Orchestrator owns scope, plan, integration, final diff, and unresolved decisions.
- Delegate bounded inventories/static searches to subagents. Give exact files/question.
- One writer per file. Subagents return evidence as `path` + symbol/test name, no prose
  dump. Orchestrator verifies important claims in source before recording or editing.
- Parallelise independent inspection; serialise overlapping edits. Never let subagents
  alter tests without explicit maintainer authority.

## Core model and invariants

`"mnirs"` = tibble/data-frame subclass made only via `create_mnirs_data()`.
Canonical attributes (`mnirs_metadata`):

| attribute | contract |
|---|---|
| `nirs_device` | detected device; scalar character/`NULL` |
| `nirs_channels` | numeric signal column names |
| `time_channel` | numeric time; finite, sorted, strictly increasing after import repair |
| `event_channel` | optional character labels or integer-like laps |
| `sample_rate` | positive scalar Hz |
| `start_timestamp` | optional POSIXct origin |
| `interval_times` | extracted start/end metadata |
| `interval_span` | extraction span metadata |

Preserve column classes/order, row order, `mnirs` first in class, tibble/grouped class,
and metadata. Explicit channel args override metadata. Re-specified `nirs_channels`
replace, not append, metadata. `verbose` omitted ->
`getOption("mnirs.verbose", TRUE)` where supported.

Package-wide processing invariant: `time_channel` finite, sorted, strictly increasing.
Import may expose invalid vendor time so `resample_mnirs()` can repair it; reject missing,
repeated, or decreasing time before downstream algorithms that depend on ordering.

Accepted input varies by function. Processing/analysis entry points generally accept:

- one data frame -> one `mnirs` tibble;
- named list of data frames -> named list; extraction/plot may add class `"mnirs"`;
- `grouped_df` -> split groups through `as_data_list()`; requires suggested `dplyr`.

Do not assume every vector-level helper accepts list/grouped input.

## Happy path

```text
read_mnirs
  -> resample_mnirs
  -> replace_mnirs
  -> filter_mnirs
  -> [shift_mnirs | rescale_mnirs | correct_blood_volume]
  -> extract_intervals
  -> analyse_kinetics
  -> plot/print
```

Bracketed transforms optional; order depends scientific question. Do not state whole
pipeline compulsory. Typical constraints: regularise duplicates/irregular time before
ensemble/filtering; replace missing values before spline/Butterworth; preserve raw data.

## Cross-cutting machinery

- `R/validate_mnirs.R`: numeric/data/channel/time/event/sample-rate/window validation,
  `validate_fix()`, time helpers, caller-aware errors.
- `R/channel_args.R`: `parse_channel_name()` + `resolve_channel_args()` broadcast scalar
  or named-list args per signal/group; `validate_group_channels()` resolves
  `"ensemble"`, `"distinct"`, custom groups, omitted channels.
- `R/as_data_list.R`: normalises single/list/grouped inputs;
  `map_mnirs_intervals()` re-evaluates original call per frame.
- `R/aanalyse_kinetics_helpers.R`: intentionally current filename; method aliases,
  fit-window detection, interval/channel workers, diagnostics/result assembly. Do not
  rename until load-order intent checked.
- Tidy evaluation: channel args accept strings, symbols, external vectors, selected
  expressions. Capture with `enquo()`/`enquos()`; resolve against data via
  `parse_channel_name()`. Test all input forms when touching this layer.
- Per-channel/group args: scalar applies globally; named list overrides matching
  channel/group; one unnamed element is fallback; missing fallback can mean skip.

## Source index

| concern | source; main symbols |
|---|---|
| import/constructor | `R/read_mnirs.R`; `read_mnirs()`, `create_mnirs_data()`, `example_mnirs()` |
| device/file parsing | `R/read_mnirs_helpers.R`; detection, `read_file()`, table/time parsing |
| regularisation | `R/resample_mnirs.R`; `resample_mnirs()` |
| clean | `R/replace_mnirs.R`, `R/replace_helpers.R`; invalid/outlier/missing replacement, rolling windows |
| filter | `R/filter_mnirs.R`; S3 generic + spline/Butterworth/moving-average methods |
| transforms | `R/shift_mnirs.R`, `R/rescale_mnirs.R`, `R/correct_blood_volume.R` |
| intervals | `R/extract_intervals.R`, `R/extract_interval_helpers.R`; `by_*()`, boundary/group/ensemble logic |
| analysis dispatch | `R/analyse_kinetics.R`; generic, method S3 wrappers, US alias |
| analysis engine | `R/aanalyse_kinetics_helpers.R`; window, broadcast, diagnostics, result builder |
| response/slope | `R/analyse_response_time.R`, `R/analyse_peak_slope.R` |
| exponential | `R/analyse_monoexponential.R`, `R/analyse_biexponential.R` |
| sigmoid | `R/analyse_sigmoidal.R`; logistic/Gompertz families |
| plot/print | `R/plot.mnirs.R`, `R/mnirs_methods.R`; optional ggplot2/scales |
| package/data docs | `R/mnirs-package.R`, `R/data.R` |

Tests mirror source stems. Large integration suites:
`test-read_mnirs.R`, `test-extract_intervals.R`, `test-analyse_kinetics.R`.
Small cross-cutting suites: `test-channel_args.R`, `test-call-attribution.R`,
`test-validate_mnirs.R`. Fixtures: `tests/testthat/testdata/`; shipped examples:
`inst/extdata/`.

## Import path

`read_mnirs()` -> `read_file()` -> detect delimiter/header/device/channels -> parse
time/sample rate/start timestamp -> select/rename columns -> `create_mnirs_data()`.
Supported via current reader path: CSV/TXT/TSV and XLS/XLSX. Auto-detection is heuristic;
explicit named mappings `c(new = "source")` are authority. When debugging vendor files,
record extension, header rows, exact column names/types, locale decimal, timestamp form,
duplicates, missingness, and expected device/channel mapping. Never commit private data;
reduce to smallest synthetic or approved fixture.

## Processing semantics

- `resample_mnirs()`: target regular grid. `method="none"` nearest-match only;
  interpolation applies when regularising/up-sampling; numeric down-sampling contract is
  time-weighted averaging. Integer/non-numeric columns use bin/LOCF logic. Updates
  time/sample-rate metadata. Current numeric down-sampling deviation: MNIRS-004.
- `replace_mnirs()`: invalid -> local MAD outliers -> missing. `width` = samples;
  `span` = time units. Vector helpers expose stages.
- `filter_mnirs()`: S3 dispatch. Spline uses stats; Butterworth requires suggested
  `signal`; moving average uses rolling helpers. NA/duplicate constraints differ.
- `shift_mnirs()` preserves amplitude; `rescale_mnirs()` changes range. Both share
  custom channel grouping and per-group args. Treat group-key order/names as invariant.
- `correct_blood_volume()` accepts one data frame only; list/grouped input unsupported by
  design. Requires at least two of oxy/deoxy/total; derives third and applies correction
  only to selected channels.

## Interval and kinetics semantics

- Boundaries: `by_time()`, `by_label()` (regex unless `fixed=TRUE`), `by_lap()`,
  `by_sample()`. `start`/`end` recycle; `span` adjusts boundaries.
- `group_intervals="distinct"` returns each interval; `"ensemble"` averages; custom
  integer index lists group selected intervals and retain omitted ones separately.
- Kinetics generic assigns method class then S3-dispatches: `response_time`,
  `peak_slope`, `monoexponential`, `biexponential`, `sigmoidal`.
- Parametric methods use self-starting `nls`; `fix` holds global parameters. Failed
  convergence returns method-shaped NA coefficients plus warning when verbose.
- `start_time` defines onset/reference; `direction` selects fit window; `end_window`
  truncates after extreme. Keep elapsed vs absolute time frames explicit.
- Biexponential `direction` affects fit window only through `find_kinetics_idx()`; it does
  not constrain component/amplitude signs.
- `mnirs_kinetics` contains `method`, per-interval/channel `model`, `coefficients`,
  augmented `data` (`<channel>_fitted`), `interval_times`, `diagnostics`,
  `channel_args`, and normalised `call`.

## Debug checklist

1. Reproduce with minimal frame retaining class/attributes; print `str(attributes(x))`.
2. Identify layer: parse -> validation -> broadcast/group -> numeric kernel -> metadata
   reconstruction -> plot/print.
3. Compare direct vector helper, single `mnirs`, named list, grouped input as relevant.
4. Probe NA/NaN/Inf, duplicate/unsorted time, empty/constant signal, one row/channel,
   irregular sampling, boundaries, duplicate group keys, optional package absent.
5. For fits: inspect resolved `channel_args`, valid window indices, start values/fixed
   params, model object, fitted vector alignment, diagnostics parameter count.
6. Assert values + type/class/names/order/attributes + condition class/message/call.
7. Ask maintainer to run focused R-console test and report output; then wider suite/check.

Suggested maintainer commands in interactive R (agent does not execute):

```r
testthat::test_file("tests/testthat/test-<stem>.R")
testthat::test_local()
devtools::document()
devtools::check()
pkgdown::check_pkgdown()
```

CI authority: `.github/workflows/R-CMD-check.yaml` (macOS/Windows/Ubuntu; devel,
release, oldrel), `test-coverage.yaml`, `pkgdown.yaml`. Local static verification:
`git diff --check`, `git diff --stat`, targeted `rg`, generated-file consistency review.

## Dependencies and generated artefacts

Imports: cli, data.table, lifecycle, readxl, rlang, stats, tibble, tidyselect.
Suggested capability boundaries: dplyr grouped input; ggplot2/scales plots; signal
Butterworth; knitr/quarto docs; zoo tests. Guard suggested packages with
`requireNamespace()`/existing helper and actionable condition.

`README.md`, `NAMESPACE`, `man/*.Rd`, rendered vignette/site files are outputs.
`AGENTS.md` and `CLAUDE.md` excluded by `.Rbuildignore`. Ignore separate shipped
`inst/skills/AGENTS.md`; it is not development authority.

## Maintainer decisions

- Target 0.7.0; ignore `inst/skills/AGENTS.md`.
- Tests always require explicit edit approval.
- Numeric down-sampling uses time-weighted averaging.
- Time is finite, sorted, strictly increasing package-wide.
- `correct_blood_volume()` remains single-frame only.
- Biexponential `direction` controls fit window only via `find_kinetics_idx()`.

## Maintainer clarification needed

- For shift min/max windows, should edge windows be excluded (`partial=FALSE`)?
- For constant signals, should `rescale_mnirs()` leave values unchanged, map to lower
  bound, midpoint, or error?
