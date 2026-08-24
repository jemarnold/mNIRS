#' Analyse kinetics across mNIRS channels and intervals
#'
#' @description
#' Model oxygenation kinetics (response time course) with various parametric
#' and non-parametric methods.
#'
#' @param data A data frame, a list of data frames, or a grouped data frame of
#'   class *"mnirs"* containing time series data and metadata (see *Details*).
#' @param method A character string specifying the kinetics analysis method.
#'   Additional arguments must be specified for each method. See *Details*.
#'   \describe{
#'      \item{`"response_time"`}{Fractional (e.g. 50%, 63.2%, 90%) response
#'      time. Additional arguments: `fraction`. See [response_time()].}
#'      \item{`"peak_slope"`}{Peak rolling linear regression slope. Additional
#'      arguments: `width` or `span`, `align`, `partial`, `na.rm`. See
#'      [peak_slope()].}
#'      \item{`"monoexponential"`}{Monoexponential curve fit via
#'      [stats::nls()]. Additional arguments: `use_TD`. See
#'      [monoexponential()].}
#'      \item{`"biexponential"`}{Biexponential (fast + slow component) curve
#'      fit via [stats::nls()]. See [biexponential()].}
#'      \item{`"exponential_drift"`}{Monoexponential curve with a secondary
#'      linear drift, fit via [stats::nls()]. Additional arguments: `use_TD`,
#'      `tau_mult`. See [exponential_drift()].}
#'      \item{`"sigmoidal"`}{Logistic or Gompertz-family curve fit via
#'      [stats::nls()]. Additional arguments: `shape`. See [logistic()].}
#'   }
#' @param start_time A numeric value in units of `time_channel` specifying the
#'   time of response onset (effectively time = `0` of the response). If `NULL`
#'   (*default*), retrieves `interval_times` from *"mnirs"* metadata, or falls
#'   back to `0` or the first positive time value (see *Details*).
#' @param end_window A numeric value in units of `time_channel` specifying the
#'   window in which to look for the end of the kinetics fit; returns the
#'   window with no greater/lesser values within `end_window` after the first
#'   extreme value. `end_window = Inf` (*default*) returns the global extreme
#'   from the full sample range (see *Details*).
#' @param ... Additional arguments passed to the underlying method function.
#'   See *Details*.
#' @param fraction **response_time**: A numeric value in the range
#'   `[0, 1]` specifying the fractional response amplitude to detect.
#'   Defaults to `0.5` (50% response, i.e. half-response time).
#' @param width **peak_slope**: An integer defining the local window in
#'   number of samples around `idx` in which to perform the operation,
#'   according to `align`. Only one of either `width` or `span` must be defined.
#' @param span **peak_slope**: A numeric value defining the local window
#'   time span around `idx` in which to perform the operation, according
#'   to `align`. In units of `time_channel`. Only one of either `width` or
#'   `span` must be defined.
#' @param align **peak_slope**: Window alignment as *"centre"/"center"*
#'   (the *default*), *"left"*, or *"right"*. Where *"left"* is forward
#'   looking, and *"right"* is backward looking from the current
#'   sample by the `width` or `span`.
#' @param partial **peak_slope**: Logical; default is `FALSE`, requires
#'   local windows to have complete number of samples specified by
#'   `width` or `span`. If `TRUE`, processes available samples within the
#'   local window. See *Details*.
#' @param na.rm **peak_slope**: Logical; default is `FALSE`, propagates
#'   any `NA`s to the returned vector and may return errors or warnings.
#'   If `TRUE`, ignores `NA`s and processes available valid samples within
#'   the local window. (see *Details*).
#' @param use_TD **monoexponential, biexponential, exponential_drift**:
#'   Logical; default is
#'   `TRUE`, attempts to fit the model with a "time-delay" parameter `TD`
#'   before the response onset. i.e., a 4-parameter [SSmonoexponential()] model
#'   (A, B, tau, TD); or a 7-parameter [biexponential()] model (A, B1, tau1,
#'   B2, tau2, tau_mult, TD). If `use_TD = FALSE` or the fit fails (with a
#'   warning), attempts to fall back to a reduced model without `TD`.
#' @param shape **sigmoidal**: Character; the 4-parameter sigmoidal shape
#'   to fit. One of `"symmetric"` (*default*; calls [SSlogistic()]),
#'   `"gompertz"` (early-inflection; calls [SSgompertz()]), or
#'   `"gompertz_left"` (late-inflection; calls [SSgompertz_left()]).
#' @param fix **monoexponential, biexponential, exponential_drift,
#'   sigmoidal**: An *optional*
#'   named list of model parameters to hold constant during fitting, e.g.
#'   `fix = list(A = 0)` fixes the starting amplitude at `0`. Fixed
#'   parameters are excluded from estimation and reported at their fixed
#'   values. Specify per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. See *Details*.
#' @param tau_mult **biexponential, exponential_drift**: A numeric multiple
#'   of the fast time constant after `TD` at which the second phase begins
#'   (`texc = TD + tau_mult * tau1`, or `TD + tau_mult * tau`). For
#'   *"biexponential"* the default `NULL` estimates `tau_mult` as a free
#'   parameter bounded below by `1`; a supplied value is held constant during
#'   fitting. For *"exponential_drift"* the default `3` (~95% of the primary
#'   amplitude) is always held constant. Specify per-channel as a list keyed
#'   by channel name, e.g. `tau_mult = list(smo2 = 2)`. See *Details*.
#' @inheritParams validate_mnirs
#' @inheritParams find_kinetics_idx
#'
#' @details
#' ## Data input formats
#'
#' `analyse_kinetics()` accepts `data` in multiple formats:
#'
#' - A **single *"mnirs"* data frame** is processed as a single interval.
#' - A **list of *"mnirs"* data frames**: each interval is processed separately.
#' - A **grouped *"mnirs"* data frame**, e.g. with `dplyr::group_by()`: the
#'   data frame is split by grouping levels and each group is processed as a
#'   separate interval.
#'
#' Specified `nirs_channels` (or channels retrieved from *"mnirs"* metadata)
#' will be analysed and results returned as a formatted table.
#'
#' ## Response **start_time** and the baseline window
#'
#' `start_time` should be specified as the time point separating the
#' pre-response baseline (`time_channel <= start_time`) from the start of the
#' response fit window (`time_channel > start_time`). For intervals extracted
#' with [extract_intervals()], `start_time` will be retrieved from *"mnirs"*
#' metadata. Otherwise `start_time` defaults to `0` or the first positive
#' `time_channel` value.
#'
#' All methods are fitted on time *elapsed from* `start_time`, so the
#' returned time coefficients are relative to response onset (e.g. `t = 0`).
#'
#' For *"response_time"*, the baseline window before `start_time` defines the
#' mean starting amplitude `A` directly and anchors the start of the
#' `response_time` parameter. For *"peak_slope"*, `start_time` anchors the
#' start of the `peak_slope_time` parameter. For *"monoexponential"*,
#' *"biexponential"*, and *"sigmoidal"*, the baseline window before
#' `start_time` anchors the starting fitted amplitude `A` and the start of the
#' `TD` and `MRT`, or `xmid` parameters. (see respective *method* sections
#' below).
#'
#' The time-delay models (*"monoexponential"* and *"biexponential"* with
#' `use_TD = TRUE`) are flat at `A` before `TD`, so the pre-onset baseline is
#' included in the fit and anchors `A`. Their reduced forms (`use_TD = FALSE`,
#' or a `TD` fit that failed and fell back) have no such flat region and are
#' fitted only where `time_channel >= start_time`.
#'
#' ## Response **direction** and the fit **end_window**
#'
#' By default, `direction` is detected automatically as either *"positive"*
#' (upward) or *"negative"* (downward) response, and can be overwritten
#' manually. `end_window` is a time span in units of `time_channel` and
#' defines the end of the kinetics fitting window by locating the first peak
#' or trough value (depending on direction) that has no greater/lesser values
#' within the subsequent `end_window` time span. The curve fitting window
#' extends to the end of `end_window` beyond the detected peak/trough.
#'
#' For *"monoexponential"*, *"biexponential"*, *"exponential_drift"*, and
#' *"sigmoidal"* methods,
#' `direction` also constrains the sign of the fitted amplitude `B - A`, and
#' the sigmoidal `slope`. For the *"biexponential"* method, `direction`
#' constrains the sign of the overall fitted amplitude `B2 - A`. A fit that
#' cannot satisfy the requested direction returns `NA` coefficients with a
#' warning.
#'
#' ## Per-channel arguments
#'
#' Arguments apply globally to all `nirs_channels` by default. Arguments can
#' instead be uniquely supplied per-channel as a named `list()` with names
#' matching `nirs_channels`, e.g.:
#'
#' ```r
#' analyse_kinetics(
#'     data,
#'     nirs_channels = c(smo2_left, smo2_right),
#'     method = "peak_slope",
#'     span = list(10, smo2_right = 20),
#'     direction = list(smo2_left = "negative", smo2_right = "auto")
#' )
#' ```
#'
#' - A non-list value applies to every channel (the *default* behaviour).
#' - A `list()` named by `nirs_channels` applies per-channel values.
#' - A single unnamed value in the list will be applied to any unlisted channels
#'   (e.g. `span = list(10, smo2_right = 20)` gives `smo2_right` 20 and every
#'   other channel 10). If no unnamed fallback value in the list, channels not
#'   named in the list fall back to the argument's default.
#' - `list()` names not matching `nirs_channels` are warned about and
#'   ignored.
#'
#' `start_time`, `direction`, and `end_window` are per-channel capable, along
#' with the `method`-specific arguments. `fix` is itself a named `list()` of
#' model parameters, so a per-channel `fix` is supplied as a `list()` of
#' `list()`s keyed by channel name. A plain parameter list applies to every
#' channel:
#'
#' ```r
#' ## fix `A` at 0 for every channel
#' fix = list(A = 0)
#'
#' ## fix `A` per-channel, leaving other channels free
#' fix = list(smo2_left = list(A = 0), smo2_right = list(A = 5, B = 100))
#' ```
#'
#' ## Per-interval arguments
#'
#' For multi-interval input (a list of data frames or a grouped data frame),
#' the same arguments can also be supplied per-interval as a named `list()`
#' keyed by interval name (the list names, group keys, or `interval_<n>`),
#' e.g.:
#'
#' ```r
#' analyse_kinetics(
#'     data_list,
#'     method = "monoexponential",
#'     end_window = list(interval_1 = 30, interval_2 = 60),
#'     direction = list(
#'         interval_1 = list(smo2_left = "negative", "auto"),
#'         interval_2 = "positive"
#'     )
#' )
#' ```
#'
#' - A resolved per-interval value may itself be a per-channel `list()`,
#'   as `direction` above.
#' - A single unnamed value in the list is the fallback for unlisted
#'   intervals; otherwise intervals not named in the list fall back to the
#'   argument's default.
#' - `list()` names matching neither interval names nor `nirs_channels`
#'   are warned about and ignored.
#' - A per-interval `fix` is a `list()` of parameter lists keyed by
#'   interval name, e.g. `fix = list(interval_1 = list(A = 0))`; each
#'   element may itself be a per-channel map, e.g.
#'   `fix = list(interval_1 = list(smo2 = list(A = 0)))`. Triple nested
#'   `list()`s is janky, but it works for now.
#'
#' `method` currently only accepts a single value applied globally to all
#' intervals and `nirs_channels`. This is a current limitation (as of `0.7.1`)
#' and will be improved in future updates to allow more flexible kinetics
#' fitting.
#'
#' ## method = "response_time"
#'
#' Aliases:
#' `method = c("response time", "half recovery time", "half time", "HRT")`.
#'
#' A non-parametric approach (estimated directly from the observed data without
#' assuming a specific mathematical shape) to estimate the response time at
#' which a signal reaches a specified fraction of its total response amplitude
#' relative to the baseline. e.g. *half-response time* (`fraction = 0.5`) is
#' the time from response onset to attain 50% of the total amplitude change
#' and approximates the inflection point (`xmid` of a symmetrical sigmoid
#' function). `fraction = 0.632` approximates the time constant (`tau`;
#' \eqn{\tau}) parameter from a monoexponential function, or the inflection
#' point (`xmid`) of an asymmetrical left-Gompertz function. `fraction = 0.368`
#' approximates `xmid` of a right-Gompertz function.
#'
#' The target response value is: `fitted = A + (B - A) * fraction`
#'
#' Where `A` is the mean baseline value (`time_channel <= start_time`) and `B`
#' is the first local extreme (peak or trough) value with no greater extreme
#' values within `end_window`. `response_value` is the first observed sample
#' where the signal is equal to or greater/lesser than the target `fitted`
#' value. `response_time` is the elapsed time from `start_time` to
#' `response_value`. See [response_time()] for the full algorithm and
#' coefficients.
#'
#' ## method = "peak_slope"
#'
#' Aliases: `method = c("peak slope", "slope", "lm")`.
#'
#' A semi-parametric approach to estimate the maximum positive or negative
#' local linear slope of a signal using rolling least-squares regression. The
#' steepest local rate of change can be interpreted as the moment of greatest
#' mismatch between oxygen delivery and extraction. `peak_slope_time` is the
#' time from response onset `start_time` to this moment of greatest mismatch.
#'
#' The local window is defined by either `width` (number of samples) or `span`
#' (in units of `time_channel`). See [peak_slope()] for window mechanics,
#' partial-window behaviour, and the returned vector-level list.
#'
#' ## method = "monoexponential"
#'
#' Aliases: `method = c("monoexp", "exponential", "exp", "tau", "MRT")`.
#'
#' A parametric approach fitting a self-starting monoexponential function to
#' the response curve using [stats::nls()] with [SSmonoexponential()] for either
#' a 4-parameter (A, B, tau, TD) or 3-parameter (A, B, tau) model.
#'
#' Model equations:
#'
#' - 3-parameter: `A + (B - A) * (1 - exp(-t / tau))`
#' - 4-parameter: `A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau))`
#'
#' `TD` is the *time delay* from `start_time` to the onset of the exponential
#' response curve. `tau` is the *time constant* of the response. The
#' *rate constant* `k` is the reciprocal (`k = 1 / tau`). The
#' *mean response time* is the time sum `MRT = TD + tau`. See
#' [monoexponential()] for the model family and [SSmonoexponential()] for
#' self-start initialisation.
#'
#' Any parameter may be held constant with `fix`, e.g. `fix = list(A = 0)`.
#' This excludes them from the fit optimisation procedure, and effectively
#' reduces the function to a lower-parameter model. `TD` can only be fixed when
#' `use_TD = TRUE` and disables the 3-parameter fallback. It is recommended to
#' specify `use_TD = FALSE` rather than fix `TD = 0`.
#'
#' ## method = "biexponential"
#'
#' Aliases: `method = c("biexp", "double exponential")`.
#'
#' A parametric approach fitting a self-starting biexponential
#' excursion-recovery function to the response curve using [stats::nls()] with
#' [SSbiexponential()]. A *fast* component (`B1`, `tau1`) drives the initial
#' response toward the excursion point; a *slow* component (`B2`, `tau2`)
#' begins at the excursion point `texc = TD + tau_mult * tau1` and recovers
#' toward a stable plateau. The excursion point is reported as `texc` and
#' `texc_fitted`.
#'
#' Model equations:
#'
#' - 6-parameter: `A + (B1 - A) * (1 - exp(-t / tau1)) +
#'   (B2 - B1) * (1 - exp(-pmax(t - tau_mult * tau1, 0) / tau2))`
#' - 7-parameter, where `ts = pmax(t - TD, 0)`:
#'   `A + (B1 - A) * (1 - exp(-ts / tau1)) +
#'   (B2 - B1) * (1 - exp(-pmax(ts - tau_mult * tau1, 0) / tau2))`
#'
#' `A` is the starting value. `B1` & `tau1` are the asymptote and time
#' constant of the fast response. `B2` & `tau2` are the asymptote and time
#' constant of the slower response plateau as time approaches `Inf` (typically
#' `tau2 >> tau1`). All three of `A`, `B1`, and `B2` are values on the response
#' scale, consistent with the [monoexponential()] and sigmoidal asymptotes.
#' Set `use_TD = TRUE` (*default*) to specify the time-delay parameter `TD`.
#' See [biexponential()] for the model family and [SSbiexponential()] for
#' self-start initialisation.
#'
#' `tau_mult` sets the slow-phase onset in multiples of `tau1` after `TD`.
#' The default `NULL` estimates it freely, bounded below by `1` so the slow
#' phase begins no earlier than one fast time constant, and above so the
#' onset stays inside the record. A supplied value is held constant and
#' silently takes precedence over a `fix = list(tau_mult = )` entry. The slow
#' time constant `tau2` is capped at ten times the record span: a slow tail
#' far beyond the record identifies only its rate, not `tau2` and `B2`
#' separately.
#'
#' Any parameter may be held constant with `fix`, e.g. `fix = list(A = 0)`, as
#' above.
#'
#' ## method = "exponential_drift"
#'
#' Aliases: `method = c("exp_linear", "exp_drift", "monoexp_drift", "drift")`.
#'
#' A parametric approach fitting a two-phase curve using [stats::nls()] with
#' [SSexponential_drift()]: a pure [monoexponential()] primary response plus
#' a secondary linear drift beginning at the excursion point `texc` near the
#' primary asymptote.
#'
#' Model equation:
#'
#' `A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau)) +
#' slope * pmax(t - TD - tau_mult * tau, 0)`
#'
#' `A`, `B`, `tau`, `TD`, and the derived `k`, `MRT`, and `HRT` are as for
#' *"monoexponential"*. `slope` is the linear drift rate `dx/dt`. The drift
#' onset is not fitted: it is held at `tau_mult` multiples of `tau` after
#' `TD` (*default* `3`; ~95% of the primary amplitude) and reported as the
#' excursion point `texc = TD + tau_mult * tau` elapsed from `start_time`
#' (the same frame as `TD` and `MRT`). Set `use_TD = TRUE` (*default*) to
#' include the time-delay parameter `TD`.
#'
#' `A`, `B`, `tau`, `slope`, and `TD` may be held constant with `fix`, as
#' above.
#'
#' ## method = "sigmoidal"
#'
#' Aliases: `method = c("logistic", "gompertz", "xmid")`.
#'
#' A parametric approach fitting a self-starting 4-parameter sigmoidal function
#' to the response curve using [stats::nls()] in one of three shapes.
#'
#' Model equations (all 4-parameter):
#'
#' - `shape = "symmetric"` ([SSlogistic()]):
#'   `A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))`
#' - `shape = "gompertz"` ([SSgompertz()]):
#'   `A + (B - A) * exp(-exp(-k * (t - xmid)))` with `k = slope * e / (B - A)`.
#'   Early-acceleration; inflection height fixed at `A + (B - A) / e`; 36.8%
#'   of the amplitude.
#' - `shape = "gompertz_left"` ([SSgompertz_left()]):
#'   `A + (B - A) * (1 - exp(-exp(k * (t - xmid))))` with
#'   `k = slope * e / (B - A)`. Late-acceleration; inflection height fixed at
#'   `A + (B - A) * (1 - 1/e)`; 63.2% of the amplitude.
#'
#' `xmid` is the time from `start_time` to the *inflection point*; the steepest
#' point of the response. `slope` is the response rate `dx/dt` at the
#' inflection.
#'
#' A *"symmetric"* shape is the default for an unbiased fit when no
#' obvious asymmetry is expected. *"gompertz"* (right-inflection) growth is
#' appropriate for fast-onset, slow-tail responses. *"gompertz_left"* for
#' slow-onset, fast-tail responses. See [logistic()], [gompertz()], and
#' [gompertz_left()] for the model families and [SSlogistic()], [SSgompertz()],
#' and [SSgompertz_left()] for self-start initialisations.
#'
#' Parameters may be held constant with `fix`, e.g. `fix = list(A = 0)`, as
#' above.
#'
#' @returns A formatted table of results, with individual elements accessible
#'   as a structured list of class *"mnirs_kinetics"* containing:
#'
#'   \item{`method`}{The method used, e.g. `"response_time"`.}
#'   \item{`model`}{A named list of model objects (per interval, per
#'       `nirs_channel`). For `"peak_slope"`; each element is an
#'       [lm][stats::lm] object. For `"monoexponential"`,
#'       `"biexponential"`, and `"sigmoidal"`; an [nls][stats::nls] object.
#'       For `"response_time"`; `NULL`. Models are fitted on time
#'       *elapsed from* `start_time`, so [predict][stats::predict] expects
#'       `.t` in those units. The offset for each interval can be retrieved
#'       from `interval_times$start_times`.}
#'   \item{`coefficients`}{A [tibble][tibble::tibble-package] of coefficients
#'       with one row per `nirs_channel` per interval, containing
#'       method-specific parameters.}
#'   \item{`data`}{A list of the original input data frames augmented with a
#'       `*_fitted` column of fitted values for each processed `nirs_channel`.}
#'   \item{`interval_times`}{A data frame with one row per interval and
#'       numeric column `start_times` -- the resolved response onset used for
#'       fitting (the supplied `start_time`, else the [extract_intervals()]
#'       metadata, else 0 or the first positive time value) -- and `end_times`
#'       when any interval carries an end time from the metadata.}
#'   \item{`diagnostics`}{A data frame of model diagnostics (`n_obs`,
#'       `n_params`, `r2`, `adj_r2`, `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`,
#'       `bic`) with one row per `nirs_channel` per interval. `n_params`
#'       counts the free parameters estimated by the solver, excluding any
#'       held by `fix`, so a reduced-parameter fallback fit is distinguishable
#'       from a full one. Only model fits sharing an `n_obs` and `n_params`
#'       are fairly comparable by `aic`/`bic`.}
#'   \item{`channel_args`}{A data frame of the resolved arguments used for
#'       each `nirs_channel` with one row per `nirs_channel` per interval.}
#'   \item{`warnings`}{A data frame of warning and error messages captured
#'       during fitting, with columns `interval`, `nirs_channels` (`NA` for
#'       interval-level warnings), `type` (`"warning"` or `"error"`), and
#'       `message`; zero rows when none occurred. Conditions are captured
#'       regardless of `verbose`, which controls console output only.}
#'   \item{`call`}{The matched call.}
#'
#' @seealso [extract_intervals()], [response_time()], [peak_slope()],
#'   [monoexponential()], [biexponential()], [logistic()], [gompertz()],
#'   [gompertz_left()]
#'
#' @examples
#' result <- read_mnirs(
#'     file_path = example_mnirs("train.red"),
#'     nirs_channels = c(
#'         smo2_left = "SmO2 unfiltered",
#'         smo2_right = "SmO2 unfiltered"
#'     ),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     zero_time = TRUE,
#'     verbose = FALSE
#' ) |>
#'     resample_mnirs(method = "linear", verbose = FALSE) |>
#'     extract_intervals(
#'         group_intervals = "distinct",
#'         start = by_time(368, 1084),
#'         span = c(-20, 90),
#'         zero_time = TRUE,
#'         verbose = FALSE
#'     ) |>
#'     analyse_kinetics(
#'         nirs_channels = c(smo2_left, smo2_right),
#'         method = "peak_slope",
#'         span = 10,          ## 10-second rolling window
#'         direction = "auto", ## auto-detect slope direction
#'         verbose = FALSE
#'     )
#'
#' ## formatted table of results
#' result
#'
#' ## coefficients are accessible from the result list
#' result$coefficients
#'
#' ## along with diagnostics and other returned objects
#' result$diagnostics
#'
#' ## plot results
#' plot(result)
#'
#' @export
analyse_kinetics <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method = c(
        "response_time",
        "peak_slope",
        "monoexponential",
        "biexponential",
        "exponential_drift",
        "sigmoidal"
    ),
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    fraction = 0.5,
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE,
    na.rm = FALSE,
    use_TD = TRUE,
    shape = c("symmetric", "gompertz", "gompertz_left"),
    tau_mult = NULL,
    fix = NULL
) {
    ## normalise method aliases before matching
    key <- gsub("[ -]", "_", tolower(method))
    method <- unname(
        ifelse(key %in% names(method_aliases), method_aliases[key], method)
    )
    method <- match.arg(method)

    UseMethod(
        "analyse_kinetics",
        structure(data, class = c(method, "mnirs_kinetics"))
    )
}

#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.response_time <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    fraction = 0.5
) {
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_response_time,
        "response_time",
        mget(unlist(kinetics_dispatch[c("common", "response_time")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.peak_slope <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE,
    na.rm = FALSE
) {
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_peak_slope,
        "peak_slope",
        mget(unlist(kinetics_dispatch[c("common", "peak_slope")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.monoexponential <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    use_TD = TRUE,
    fix = NULL
) {
    ## TODO: pass additional stats::nls() args
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_monoexponential,
        "monoexponential",
        mget(unlist(kinetics_dispatch[c("common", "monoexponential")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.biexponential <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    use_TD = TRUE,
    tau_mult = NULL,
    fix = NULL
) {
    ## TODO: pass additional stats::nls() args
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_biexponential,
        "biexponential",
        mget(unlist(kinetics_dispatch[c("common", "biexponential")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.exponential_drift <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    use_TD = TRUE,
    tau_mult = 3,
    fix = NULL
) {
    ## TODO: pass additional stats::nls() args
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_exponential_drift,
        "exponential_drift",
        mget(unlist(kinetics_dispatch[c("common", "exponential_drift")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.sigmoidal <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    shape = c("symmetric", "gompertz", "gompertz_left"),
    fix = NULL
) {
    ## TODO: pass additional stats::nls() args
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    return(analyse_kinetics_intervals(
        data,
        analyse_logistic,
        "sigmoidal",
        mget(unlist(kinetics_dispatch[c("common", "sigmoidal")])),
        enquo(nirs_channels),
        enquo(time_channel),
        verbose,
        match.call(),
        sys.call(-1)
    ))
}


#' @rdname analyse_kinetics
#' @export
analyze_kinetics <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method = c(
        "response_time",
        "peak_slope",
        "monoexponential",
        "biexponential",
        "exponential_drift",
        "sigmoidal"
    ),
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...
) {
    call <- match.call()
    call[[1L]] <- quote(analyse_kinetics)
    eval(call, envir = parent.frame())
}
