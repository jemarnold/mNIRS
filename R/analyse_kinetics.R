#' Analyse kinetics across mNIRS channels and intervals
#'
#' @description
#' Fit parametric or non-parametric kinetics for each `nirs_channel` within an
#' *"mnirs"* data frame, a list of data frames, or a grouped data frame.
#' 
#' Note the `method`-specific arguments below.
#'
#' @param data A data frame of class *"mnirs"* containing time series data and
#'   metadata, a list of data frames, or a grouped data frame (see *Details*).
#' @param method A character string specifying the kinetics analysis method.
#'   Additional arguments must be specified for each method. See *Details*.
#'   \describe{
#'      \item{`"response_time"`}{Fractional (e.g. 50%, 63.2%, 90%) response
#'      time. Additional arguments: `fraction`. See [response_time()].}
#'      \item{`"peak_slope"`}{Peak local linear regression slope. Additional
#'      arguments: `width` or `span`, `align`, `partial`, `na.rm`. See
#'      [peak_slope()].}
#'      \item{`"monoexponential"`}{Monoexponential curve fit via
#'      [stats::nls()]. Additional arguments: `use_time_delay`. See
#'      [monoexponential()].}
#'      \item{`"logistic"`}{`<under development>`.}
#'   }
#' @param t0 A numeric value specifying the start of the kinetics response
#'   in units of `time_channel`. Observations where `time_channel <= t0`
#'   define the pre-response baseline window. If `NULL` (*default*), retrieves
#'   `interval_times` from *"mnirs"* metadata, or falls back to `0`.
#' @param ... Additional arguments passed to the underlying method function.
#'   See *Details*.
#' @param fraction **response_time**: A numeric value in the range
#'   `[0, 1]` specifying the fractional response amplitude to detect.
#'   Defaults to `0.5` (50% response, i.e. half-response time).
#' @param width **peak_slope**: An integer defining the local window in
#'   number of samples around `idx` in which to perform the operation,
#'   according to `align`.
#' @param span **peak_slope**: A numeric value defining the local window
#'   time span around `idx` in which to perform the operation, according
#'   to `align`. In units of `time_channel`.
#' @param align **peak_slope**: Window alignment as *"centre"/"center"*
#'   (the *default*), *"left"*, or *"right"*. Where *"left"* is *forward
#'   looking*, and *"right"* is *backward looking* from the current
#'   sample.
#' @param partial **peak_slope**: Logical; default is `FALSE`, requires
#'   local windows to have complete number of samples specified by
#'   `width` or `span`. If `TRUE`, processes available samples within the
#'   local window. See *Details*.
#' @param na.rm **peak_slope**: Logical; default is `FALSE`, propagates
#'   any `NA`s to the returned vector. If `TRUE`, ignores `NA`s and
#'   processes available valid samples within the local window. May
#'   return errors or warnings. (see *Details*).
#' @param use_time_delay **monoexponential**: Logical; default is `TRUE`
#'   to attempt to fit a 4-parameter [SSmonoexp()] model (A, B, tau, TD)
#'   with a time delay. If the 4-parameter fit fails, or if
#'   `use_time_delay = FALSE`, attempts to fit a reduced 3-parameter
#'   [SSmonoexp()] model (A, B, tau).
#' @param shape **logistic**: Character; the 4-parameter sigmoidal shape
#'   to fit. One of `"symmetric"` (*default*; calls [SSlogistic()]),
#'   `"gompertz"` (early-inflection; calls [SSgompertz()]), or
#'   `"gompertz_left"` (late-inflection; calls [SSgompertz_left()]).
#' @inheritParams validate_mnirs
#' @inheritParams find_kinetics_idx
#'
#' @details
#' ## Data input formats
#'
#' `analyse_kinetics()` accepts `data` in multiple formats:
#'
#' - A **single *"mnirs"* data frame** is processed as a single interval.
#' - A **list of *"mnirs"* data frames** -- each interval is processed
#'   separately.
#' - A **grouped *"mnirs"* data frame**, e.g. with `dplyr::group_by()` --
#'   the data frame is split by grouping levels and each group is processed
#'   as a separate interval.
#'
#' Specified `nirs_channels` (or channels retrieved from *"mnirs"* metadata)
#' will be analysed and results returned as a formatted table.
#'
#' ## Response onset (t0) and the baseline window
#'
#' `t0` should be specified as the time point separating the pre-response
#' baseline (`time_channel <= t0`) from the start of the response window
#' (`time_channel > t0`). For intervals extracted with [extract_intervals()],
#' `t0` will retrieve `interval_times` from *"mnirs"* metadata. Otherwise
#' `t0` defaults to `0`.
#'
#' For `"response_time"`, the baseline window before `t0` defines the mean
#' starting amplitude `A` directly and anchors the start of the `response_time`
#' parameter. For `"monoexponential"`, `t0` anchors the response onset for the
#' time delay (`TD`) parameter, which provides the baseline window to fit `A`
#' in a 4-parameter model (see *method = "monoexponential"* section below).
#' For `"peak_slope"`, `t0` anchors the response onset for the
#' `peak_slope_time` parameter.
#'
#' ## Response direction and the end of the fitting window
#'
#' By default, `direction` is detected automatically as either *"positive"*
#' (upward response) or *"negative"* (downward response), and can be
#' overwritten manually. The end of the fitting window is set by locating the
#' first peak (positive) or trough (negative) that has no greater/lesser value
#' within a subsequnt window defined by `end_fit_span`: a time span in units of
#' `time_channel`. The fitting window extends to the end of `end_fit_span`
#' beyond the first local extreme peak/trough.
#'
#' ## method = "response_time"
#'
#' Aliases:
#' `method = c("response time", "half recovery time", "half time", "HRT")`.
#'
#' A non-parametric approach (estimated directly from the observed data without
#' assuming a specific mathetmatical shape) to estimate the response time at
#' which a signal reaches a specified fraction of its total response amplitude
#' relative to the baseline. e.g. *half-response time* (`fraction = 0.5`) is
#' the time from response onset to attain 50% of the total amplitude change.
#' `fraction = 0.632` approximates the time constant (`tau`; \eqn{\tau})
#' parameter from a monoexponential function.
#'
#' The target response value is:
#'
#' `fitted = A + (B - A) * fraction`
#'
#' Where `A` is the mean baseline value (`time_channel <= t0`) and `B` is the
#' first local extreme (peak or trough) value after `t0`. `response_time` is
#' the elapsed time from `t0` to `response_value`; the first observed sample
#' where the signal is equal to or greater/lesser than the target `fitted`
#' value. See [response_time()] for the full algorithm and coefficients results.
#'
#' ## method = "peak_slope"
#'
#' Aliases: `method = c("peak slope", "slope")`.
#'
#' A semi-parametric approach to estimate the maximum positive or negative
#' local linear slope of a signal using rolling least-squares regression. The
#' steepest local rate of change can be interpreted as the moment of greatest
#' mismatch between oxygen delivery and extraction. `peak_slope_time` is the
#' time from response onset `t0` to this moment of greatest mismatch.
#'
#' The local window is defined by either `width` (number of samples) or `span`
#' (in units of `time_channel`). See [peak_slope()] for window mechanics,
#' partial-window behaviour, and the returned vector-level list.
#'
#' ## method = "monoexponential"
#'
#' Aliases: `method = c("monoexp", "exponential", "tau", "MRT")`.
#'
#' A parametric approach fitting a self-starting monoexponential function to
#' the response curve using [stats::nls()] with [SSmonoexp()] for either
#' a 4-parameter (A, B, tau, TD) or 3-parameter (A, B, tau) model.
#'
#' Model equations:
#'
#' - 3-parameter: `A + (B - A) * (1 - exp(-t / tau))`
#' - 4-parameter: `ifelse(t <= TD, A, A + (B - A) * (1 - exp(-(t - TD) / tau)))`
#'
#' `tau` is the *time constant* of the response. The *rate constant* `k` can be
#' derived as the reciprocal (`k = 1 / tau`). The *mean response time*
#' `MRT = TD + tau` and the *half-response time* `HRT = TD + tau * log(2)`
#' can also be derived. See [monoexponential()] for the model family and
#' [SSmonoexp()] for self-start initialisation.
#'
#' Other arguments (`...`) passed to [stats::nls()] are `<NOT YET IMPLEMENTED>`.
#'
#' ## method = "logistic"
#'
#' Aliases: `method = c("sigmoidal", "xmid")`.
#'
#' A parametric approach fitting a self-starting 4-parameter sigmoidal
#' function to the response curve using [stats::nls()] in one of three
#' shapes. All three expose the same coefficients `(A, B, xmid, slope)`
#' for cross-shape and cross-channel comparability.
#'
#' Model equations (all 4-parameter):
#'
#' - `shape = "symmetric"` ([SSlogistic()]):
#'   `A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))`
#' - `shape = "gompertz"` ([SSgompertz()]):
#'   `A + (B - A) * exp(-exp(-k * (t - xmid)))`
#'   with `k = slope * e / (B - A)`. Early-acceleration; inflection
#'   height fixed at `A + (B - A) / e`.
#' - `shape = "gompertz_left"` ([SSgompertz_left()]):
#'   `A + (B - A) * (1 - exp(-exp(k * (t - xmid))))`
#'   with `k = slope * e / (B - A)`. Late-acceleration; inflection
#'   height fixed at `A + (B - A) * (1 - 1/e)`.
#'
#' `xmid` is the time at the inflection (steepest) point of the response,
#' reported relative to `t0`. `slope` is the response rate `dx/dt` at the
#' inflection. The `symmetric` shape is the default for an unbiased fit
#' when no obvious asymmetry is expected; `gompertz` is appropriate
#' for fast-onset / slow-tail responses, and `gompertz_left` for
#' slow-onset / fast-tail responses. See [logistic()], [gompertz()],
#' and [gompertz_left()] for the model families.
#'
#' Other arguments (`...`) passed to [stats::nls()] are `<NOT YET IMPLEMENTED>`.
#'
#' @returns A formatted table of printed results, with individual elements
#'   accessible as a structured list of class *"mnirs_kinetics"* containing:
#'
#'   \item{`method`}{The method used, e.g. `"response_time"`.}
#'   \item{`model`}{A named list of model objects (per interval, per
#'       `nirs_channel`). For `"peak_slope"`, each element is an
#'       [lm][stats::lm] object; for `"monoexponential"`, an
#'       [nls][stats::nls] object; for `"response_time"`, `NULL`. `NULL`
#'       for channels where fitting failed.}
#'   \item{`coefficients`}{A [tibble][tibble::tibble-package] of coefficients
#'       with one row per `nirs_channel` per interval, containing columns
#'       `interval`, `nirs_channels`, and the method-specific parameters.}
#'   \item{`data`}{A list of the original input data frames augmented with a
#'       `*_fitted` column of fitted values for each `nirs_channel`.}
#'   \item{`interval_times`}{A data frame of interval times for each
#'       `nirs_channel` per interval, supplied from [extract_intervals()] if
#'       present in the metadata.}
#'   \item{`diagnostics`}{A data frame of model diagnostics (`n_obs`, `r2`,
#'       `adj_r2`, `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`, `bic`) with one
#'       row per `nirs_channel` per interval.}
#'   \item{`channel_args`}{A data frame of the resolved arguments used for
#'       each `nirs_channel` with one row per `nirs_channel` per interval.}
#'   \item{`call`}{The matched call.}
#'
#' @seealso [extract_intervals()], [response_time()], [peak_slope()],
#'   [monoexponential()]
#'
#' @examples
#' result <- read_mnirs(
#'     example_mnirs("train.red"),
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
#' @export
analyse_kinetics <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method = c("response_time", "peak_slope", "monoexponential", "logistic"),
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...,
    fraction = 0.5,
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE,
    na.rm = FALSE,
    use_time_delay = TRUE,
    shape = c("symmetric", "gompertz", "gompertz_left")
) {
    ## normalise method aliases before matching
    method <- gsub(
        "^HRT$|^(?:((half[ _-])?(response|recovery)[ _-]time)|half[ _-]time)$",
        "response_time",
        method,
        ignore.case = TRUE
    )
    method <- gsub(
        "^peak[ _-]slope$|^slope$",
        "peak_slope",
        method,
        ignore.case = TRUE
    )
    method <- gsub(
        "^monoexp$|^exponential$|^MRT$|^tau$",
        "monoexponential",
        method,
        ignore.case = TRUE
    )
    method <- gsub(
        "^sigmoidal$|^xmid$",
        "logistic",
        method,
        ignore.case = TRUE
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
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...,
    fraction = 0.5
) {
    ## resolve global verbose option when caller omits the argument
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        result <- analyse_response_time(
            data = data_list[[.i]],
            nirs_channels = !!enquo(nirs_channels),
            time_channel = !!enquo(time_channel),
            t0 = t0,
            fraction = fraction,
            direction = direction,
            end_fit_span = end_fit_span,
            verbose = verbose,
            interval_name = names(data_list)[[.i]],
            bypass_checks = TRUE,
            env = env
        )
    })

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(
        data_list,
        result_list,
        method = "response_time",
        match.call()
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
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...,
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE,
    na.rm = FALSE
) {
    ## resolve global verbose option when caller omits the argument
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        result <- analyse_peak_slope(
            data = data_list[[.i]],
            nirs_channels = !!enquo(nirs_channels),
            time_channel = !!enquo(time_channel),
            t0 = t0,
            width = width,
            span = span,
            align = align,
            direction = direction,
            end_fit_span = end_fit_span,
            partial = partial,
            na.rm = na.rm,
            verbose = verbose,
            interval_name = names(data_list)[[.i]],
            bypass_checks = TRUE,
            env = env
        )
    })

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(
        data_list,
        result_list,
        method = "peak_slope",
        match.call()
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
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...,
    use_time_delay = TRUE
) {
    ## TODO: pass additional stats::nls() args
    ## TODO: implement `direction`
    ## resolve global verbose option when caller omits the argument
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        analyse_monoexponential(
            data = data_list[[.i]],
            nirs_channels = !!enquo(nirs_channels),
            time_channel = !!enquo(time_channel),
            use_time_delay = use_time_delay,
            t0 = t0,
            end_fit_span = end_fit_span,
            verbose = verbose,
            interval_name = names(data_list)[[.i]],
            bypass_checks = TRUE,
            env = env
        )
    })

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(
        data_list,
        result_list,
        method = "monoexponential",
        match.call()
    ))
}


#' @rdname analyse_kinetics
#' @usage NULL
#' @export
analyse_kinetics.logistic <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...,
    shape = c("symmetric", "gompertz", "gompertz_left")
) {
    ## TODO: pass additional stats::nls() args
    ## TODO: implement `direction`
    ## resolve global verbose option when caller omits the argument
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        analyse_logistic(
            data = data_list[[.i]],
            nirs_channels = !!enquo(nirs_channels),
            time_channel = !!enquo(time_channel),
            shape = shape,
            t0 = t0,
            end_fit_span = end_fit_span,
            verbose = verbose,
            interval_name = names(data_list)[[.i]],
            bypass_checks = TRUE,
            env = env
        )
    })

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(
        data_list,
        result_list,
        method = "logistic",
        match.call()
    ))
}


#' @rdname analyse_kinetics
#' @export
analyze_kinetics <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method = c("response_time", "peak_slope", "monoexponential", "logistic"),
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...
) {
    call <- match.call()
    call[[1L]] <- quote(analyse_kinetics)
    eval(call, envir = parent.frame())
}