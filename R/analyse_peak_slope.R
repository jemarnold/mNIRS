#' Calculate linear slope
#'
#' `slope()`: Calculate the linear regression slope of a numeric vector via
#' the least-squares formula.
#'
#' @inheritParams peak_slope
#'
#' @returns `slope()` returns a numeric slope value in units of `x / t`, or
#'  `NA_real_` when insufficient valid observations are present.
#'
#' @rdname rolling_slope
#' @order 2
#' @keywords internal
slope <- function(
    x,
    t = seq_along(x),
    na.rm = FALSE,
    ...,
    env = rlang::caller_env()
) {
    ## validation =================================================
    args <- list(...)
    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t, allow_na = TRUE, env = env)
    }

    if (length(x) < max(args$min_obs, 2L) || !na.rm && anyNA(x)) {
        return(NA_real_)
    }

    ## remove invalid
    ## TODO redundant with `find_kinetics_idx`
    ## but I want NAs to get to the above `na.rm` check
    complete <- which(is.finite(x) & is.finite(t))
    x <- x[complete]
    t <- t[complete]
    n <- length(x)

    ## processing =================================================
    sum_t <- sum(t)
    sum_x <- sum(x)
    sum_tx <- sum(t * x)
    sum_t2 <- sum(t^2)
    denom <- n * sum_t2 - sum_t^2 ## should not cause integer overflow
    if (is.na(denom) || denom == 0) {
        return(NA_real_)
    }
    slope_val <- (n * sum_tx - sum_t * sum_x) / denom

    if (args$intercept %||% FALSE) {
        attr(slope_val, "intercept") <- (sum_x - slope_val * sum_t) / n
    }

    return(slope_val)
}

#' Calculate rolling linear slope
#'
#' `rolling_slope()`: Compute rolling linear regression slopes within a local
#' window along a numeric vector.
#'
#' @inheritParams peak_slope
#'
#' @details
#' See [peak_slope()] for details on window specification (`width`, `span`,
#' `align`), partial windows, and direction detection.
#'
#' Additional arguments (`...`) accepted:
#'
#' \describe{
#'   \item{`bypass_checks`}{Logical; if `TRUE`, skips input validation.
#'   Intended for internal use when checks have already been performed
#'   upstream.}
#'   \item{`min_obs`}{Integer; minimum number of valid observations required
#'   per window to return a slope. Derived from `width` or `span`, or `2L`
#'   when `partial = TRUE`.}
#'   \item{`intercept`}{Logical; if `TRUE`, `slope()` also attaches the
#'   y-intercept as `attr(slope_val, "intercept")`.}
#'   \item{`window_idx`}{Logical; if `TRUE`, the list of per-observation
#'   window indices is attached as `attr(slopes, "window_idx")`.}
#' }
#'
#' @seealso [peak_slope()]
#'
#' @returns `rolling_slope()` returns a numeric vector of rolling local slopes
#'   in units of `x / t`, the same length as `x`.
#'
#' @inheritParams validate_mnirs
#'
#' @rdname rolling_slope
#' @order 1
#' @keywords internal
rolling_slope <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    partial = FALSE,
    na.rm = FALSE,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    ## validation =================================================
    args <- list(...)
    n <- length(x)
    align <- sub("^center$", "centre", align)
    align <- match.arg(align)

    insufficient_warn <- c(
        "!" = "Insufficient valid samples detected in {.fn froll_slope}.",
        "i" = "Check length of {.arg x} and {.arg width} or {.arg span}, \\
        or specify {.arg partial} = {.val {TRUE}}."
    )

    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t, allow_na = TRUE, env = env)

        ## return NA with warning
        if (n == 0L) {
            if (verbose) {
                cli_warn(insufficient_warn, call = warn_call(env))
            }
            return(numeric(0))
        }
        if (n == 1L || all(diff(t) == 0)) {
            if (verbose) {
                cli_warn(insufficient_warn, call = warn_call(env))
            }
            return(rep(NA_real_, n))
        }

        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
        validate_width_span(width, span, verbose, env = env)
    }

    ## min_obs default to estimated width when span is specified
    min_obs <- if (partial) {
        2L
    } else {
        ## less strict span_width - 2 to allow start & end buffer
        ## with irregular t values
        max(width %||% (floor(span * estimate_sample_rate(t, env)) - 2L), 2L)
    }

    if (n < min_obs) {
        if (verbose) {
            cli_warn(insufficient_warn, call = warn_call(env))
        }
        return(rep(NA_real_, n))
    }

    ## processing =================================================
    window_idx <- compute_local_windows(
        t,
        width = width,
        span = span,
        align = align,
        env = env
    )

    if (verbose && all(lengths(window_idx) < min_obs)) {
        cli_warn(insufficient_warn, call = warn_call(env))
    }

    slopes <- vapply(window_idx, \(.idx) {
        slope(x[.idx], t[.idx], na.rm, min_obs = min_obs, bypass_checks = TRUE)
    }, numeric(1))

    if (args$window_idx %||% FALSE) {
        attr(slopes, "window_idx") <- window_idx
    }

    return(slopes)
}


#' Find peak linear slope
#'
#' Identify the maximum positive or negative local linear slope within a
#' numeric vector using rolling least-squares regression, and return a list
#' of regression parameters for the peak window. Vector-level companion to
#' [analyse_kinetics()] when `method = "peak_slope"`.
#'
#' @param ... Additional arguments.
#' @inheritParams find_kinetics_idx
#' @inheritParams compute_local_windows
#' @inheritParams replace_invalid
#' @inheritParams filter_moving_average
#'
#' @details
#' ## Rolling window
#'
#' The local rolling window is defined by either `width` (number of samples)
#' or `span` (time duration in units of `t`); one of either `width` or `span`
#' must be specified. When `align = "centre"` and
#' `width` is specified, the window spans
#' `[idx - floor((width - 1) / 2), idx + floor(width / 2)]`. Even `width`
#' values bias alignment to *"left"*, placing the unequal sample forward of
#' `idx`.
#'
#' When `span` is specified with `align = "centre"`, the window spans
#' `[t - span / 2, t + span / 2]`.
#'
#' ## Direction detection
#'
#' When `direction = "auto"`, the net slope across all of `x` determines the
#' overall trend (positive or negative), and the greatest local slope in that
#' direction is returned. If the net slope equals zero or is `NA`, the greatest
#' absolute local slope is returned.
#'
#' ## Partial windows
#'
#' The default `partial = FALSE` requires the complete number of samples
#' specified by `width` or `span`. `NA` is returned if fewer samples are
#' present in the local window.
#'
#' Setting `partial = TRUE` allows computation with as few as 2 valid samples.
#' These windows, such as at edge conditions, will be more sensitive to noise
#' and this setting should be used with caution.
#'
#' ## Missing values
#'
#' When `na.rm = FALSE` (the default), any `NA` in a window propagates `NA`
#' to the returned slope. When `na.rm = TRUE`, `NA`s are ignored and the slope
#' is computed from the remaining valid samples.
#'
#' @returns A named list containing:
#'   \item{`slope`}{The peak slope value in units of `x / t`.}
#'   \item{`intercept`}{The y-intercept of the peak local regression line.}
#'   \item{`y`}{The predicted response value at the peak slope window index.}
#'   \item{`t`}{The time value at the peak slope window index.}
#'   \item{`idx`}{The integer index position of the peak slope window.}
#'   \item{`fitted`}{A numeric vector of predicted values spanning the peak
#'   slope window.}
#'   \item{`window_idx`}{An integer vector of indices spanning the peak slope
#'   window.}
#'
#' @seealso [analyse_kinetics()], [response_time()], [monoexponential()],
#'   [rolling_slope()]
#'
#' @examples
#' x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
#' peak_slope(x, width = 5)
#'
#' x_dec <- rev(x)
#' peak_slope(x_dec, width = 5)
#'
#' @export
peak_slope <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    direction = c("auto", "positive", "negative"),
    partial = FALSE,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    args <- list(...)
    direction <- match.arg(direction)
    ## internal callers pass `env` through `...` to report conditions
    ## as coming from the user-facing function
    env <- args$env %||% environment()

    if (!(args$bypass_checks %||% FALSE)) {
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
    }

    ## calculate all rolling slopes
    slopes <- rolling_slope(
        x,
        t,
        width,
        span,
        align,
        partial,
        na.rm,
        verbose,
        window_idx = TRUE,
        bypass_checks = args$bypass_checks %||% FALSE, ## use validations
        env = env
    )

    ## pre-return NA
    na_result <- list(
        t = NA_real_,
        y = NA_real_,
        slope = NA_real_,
        intercept = NA_real_,
        idx = NA_integer_,
        fitted = NA_real_,
        window_idx = NA_integer_,
        model = NA
    )

    if (all(is.na(slopes))) {
        return(na_result)
    }

    ## detect direction from net trend, fallbackt to abs peak slope
    direction <- detect_direction(x, t, slopes, direction)

    ## manual direction calculation
    candidates <- switch(
        direction,
        positive = which(slopes > 0),
        negative = which(slopes < 0)
    )

    if (length(candidates) == 0L) {
        if (verbose) {
            cli_warn(
                c("!" = "No {direction} slopes detected."),
                call = warn_call(env)
            )
        }
        return(na_result)
    }

    ## find peak (max for positive, min for negative)
    peak_idx <- switch(
        direction,
        positive = candidates[which.max(slopes[candidates])],
        negative = candidates[which.min(slopes[candidates])]
    )

    ## get window indices at peak
    window_idx <- attr(slopes, "window_idx")[[peak_idx]]

    ## fit lm on peak window
    fit_formula <- stats::as.formula("x ~ t", env = baseenv())
    model <- stats::lm(
        fit_formula,
        data.frame(x = x[window_idx], t = t[window_idx])
    )

    slope_val <- unname(stats::coef(model)[["t"]])
    intercept_val <- unname(stats::coef(model)[["(Intercept)"]])
    t_peak <- t[peak_idx]
    y_peak <- intercept_val + slope_val * t_peak

    return(list(
        slope = slope_val,
        intercept = intercept_val,
        y = y_peak,
        t = t_peak,
        idx = peak_idx,
        fitted = unname(stats::fitted(model)),
        window_idx = window_idx,
        model = model
    ))
}


#' Analyse peak linear slope across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "peak_slope")`. Computes the maximum local
#' linear slope for each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @inheritParams validate_mnirs
#' @inheritParams peak_slope
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `slope`, `intercept`, `y`, `peak_slope_time`, `idx`.
#'   Per-channel metadata are attached as attributes:
#'   - `"model"`: a linear regression model object via [stats::lm()].
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [peak_slope()]
#'
#' @keywords internal
analyse_peak_slope <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    start_time = NULL,
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    partial = FALSE,
    na.rm = FALSE,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    ## validation ==============================================
    validate_mnirs_data(data, env = env)
    nirs_channels <- validate_nirs_channels(
        enquo(nirs_channels), data, verbose, env = env
    )
    time_channel <- validate_time_channel(enquo(time_channel), data, env = env)
    t_vec <- data[[time_channel]]
    args <- list(...)
    ## interval label; falls back to the `data` argument name when unsupplied
    interval_name <- args$interval_name %||% deparse(substitute(data))

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(
            start_time = start_time,
            width = width,
            span = span,
            align = align,
            direction = direction,
            end_window = end_window,
            partial = partial,
            na.rm = na.rm
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose,
        env = env
    )
    ## validate resolved args once, before fitting any channel
    per_channel <- validate_kinetics_args(
        per_channel,
        data,
        t_vec,
        verbose,
        env = env
    )

    ## method-specific fit: peak rolling linear slope
    peak_slope_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        slopes <- do.call(peak_slope, c(
            list(x = x_fit, t = t_fit),
            .a,
            list(verbose = verbose, bypass_checks = TRUE, env = env),
            args
        ), quote = TRUE)

        list(
            coefs = data.frame(
                slope           = slopes$slope,
                intercept       = slopes$intercept,
                fitted          = slopes$y, ## predicted response value at idx
                peak_slope_time = slopes$t - .a$start_time,
                idx             = slopes$idx
            ),
            model = slopes$model,
            fitted_data = data.frame(
                window_idx = slopes$window_idx,
                fitted     = slopes$fitted
            ),
            diag = compute_diagnostics(
                x        = x_fit[slopes$window_idx],
                t        = t_fit[slopes$window_idx],
                fitted   = slopes$fitted,
                n_params = 2L,
                verbose  = verbose,
                env      = env
            )
        )
    }

    return(analyse_kinetics_channels(
        data,
        nirs_channels,
        time_channel,
        per_channel,
        peak_slope_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
