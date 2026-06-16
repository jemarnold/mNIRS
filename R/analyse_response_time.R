#' Compute fractional kinetics response time
#'
#' Identify the time at which a numeric signal reaches a specified fraction
#' of its total response amplitude relative to a baseline period (e.g.
#' *half-response time* at 50% fractional amplitude). Vector-level companion
#' to [analyse_kinetics()] when `method = "response_time"`.
#'
#' @param t0 A numeric value specifying the start of the kinetics response
#'   in units of `t`. Observations where `t <= t0` define the baseline window.
#'   Defaults to `0`.
#' @param fraction A numeric value in the range `[0, 1]` specifying the
#'   fractional response amplitude to detect. Defaults to `0.5` (50% response,
#'   i.e. half-response time).
#' @param ... Additional arguments.
#' @inheritParams replace_invalid
#' @inheritParams find_kinetics_idx
#' @inheritParams validate_mnirs
#'
#' @details
#' ## Method
#'
#' The target response value is computed as:
#'
#' `response_fitted = A + (B - A) * fraction`
#'
#' where `A` is the mean baseline (`t <= t0`) and `B` is the extreme (peak
#' or trough) value after `t0`. The response time is the elapsed time from
#' `t0` to the first sample where `x` reaches the `response_fitted` value
#' (above for positive direction, below for negative).
#'
#' ## Direction
#'
#' When `direction = "auto"`, the net slope across `x` determines the overall
#' trend, and the corresponding extreme (maximum for positive, minimum for
#' negative) is used as `B`. If the net slope is zero or `NA`, the direction
#' of greatest absolute change is used. When `direction = "positive"` or
#' `"negative"`, the extreme is the maximum or minimum value after `t0`,
#' respectively.
#'
#' ## Baseline
#'
#' When no observations exist where `t <= t0`, the first sample `x[1]` is
#' used as the baseline and a warning is issued. `t0` cannot exceed the
#' maximum of `t`.
#'
#' @returns A named list containing:
#'   \item{`A`}{Mean baseline value (mean of `x` where `t <= t0`).}
#'   \item{`B`}{Extreme value (maximum or minimum) after `t0`.}
#'   \item{`response_time`}{Elapsed time from `t0` to the fractional response,
#'   in units of `t`.}
#'   \item{`response_value`}{The observed signal value at the response index.}
#'   \item{`fitted`}{The predicted fractional response value
#'   (`A + (B - A) * fraction`).}
#'   \item{`baseline_idx`}{Integer indices where `t <= t0`.}
#'   \item{`response_idx`}{Integer index at the `response_value`.}
#'   \item{`extreme_idx`}{Integer index at the extreme value (`B`).}
#'
#' @seealso [analyse_kinetics()], [peak_slope()], [monoexponential()]
#'
#' @examples
#' set.seed(13)
#' t <- 0:60
#' x <- monoexponential(t, A = 20, B = 60, tau = 8, TD = 10) + rnorm(length(t), 0, 1)
#'
#' ## estimated half-response time
#' HRT <- response_time(x, t, t0 = 10, fraction = 0.5)
#'
#' ## estimated mean response time (time constant; tau ~= 63.2% amplitude)
#' MRT <- response_time(x, t, t0 = 10, fraction = 0.632)
#'
#' plot(t, x, type = "l", col = "grey60", xlab = "t", ylab = "x")
#' ## baseline mean across baseline_idx
#' segments(
#'     t[min(HRT$baseline_idx)], HRT$A,
#'     t[max(HRT$baseline_idx)], HRT$A,
#'     col = "red", lwd = 2
#' )
#' ## fraction = 0.5 (red): response_value and extreme
#' points(t[HRT$response_idx], HRT$response_value, col = "red", pch = 19)
#' points(t[HRT$extreme_idx], HRT$B, col = "red", pch = 19)
#' ## fraction = 0.632 (blue): response_value
#' points(t[MRT$response_idx], MRT$response_value, col = "blue", pch = 19)
#'
#' @export
response_time <- function(
    x,
    t = seq_along(x),
    t0 = 0,
    fraction = 0.5,
    direction = c("auto", "positive", "negative"),
    verbose = TRUE,
    ...
) {
    validate_numeric(
        fraction, 1L, c(0, 1), msg2 = "between {col_blue('[0, 1]')}."
    )
    direction <- match.arg(direction)
    args <- list(...)

    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t, allow_na = TRUE)
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
        ## detect direction from net trend, fallback to abs magnitude
        direction <- detect_direction(x, t, x, direction)
    }

    baseline_idx <- which(t <= t0)

    if (!(args$bypass_checks %||% FALSE)) {
        validate_numeric(t0, 1L)
        if (length(baseline_idx) == 0L) {
            if (verbose) {
                cli_warn(c(
                    "!" = "No observations where {.arg t} <= {.arg t0} = \\
                    {.val {t0}}.",
                    "i" = "{.code x[1]} used as response baseline."
                ))
            }
            baseline_idx <- 1L
            t0 <- t[baseline_idx]
        }
        if (t0 > t[length(t)]) {
            cli_abort(c(
                "x" = "No observations in {.arg t} before {.arg t0}.",
                "i" = "{.arg t0} must be specified within the range of \\
                {.arg t}."
            ))
        }
    }

    ## process =====================================================
    ## look for extreme after t0
    x_valid <- c(rep(NA_real_, length(baseline_idx)), x[t > t0])
    extreme_idx <- if (direction == "positive") {
        which.max(x_valid)
    } else {
        which.min(x_valid)
    }

    A <- mean(x[baseline_idx], na.rm = TRUE)
    B <- x[extreme_idx]
    response_fitted <- A + (B - A) * fraction
    compare_fn <- if (direction == "positive") `>=` else `<=`
    response_idx <- which(compare_fn(x_valid, response_fitted))[1L]

    if (is.na(response_idx)) {
        if (verbose) {
            cli_warn(c(
                "!" = "No valid {.val {direction}} extremes after {.arg t0}. \\
                Returning {.val {NA}}."
            ))
        }
        response_fitted <- NA_real_
    }

    return(list(
        A = A,
        B = B,
        response_time = t[response_idx] - t0, ## real
        response_value = x[response_idx], ## real
        fitted = response_fitted,    ## predicted
        baseline_idx = baseline_idx, ## all baseline samples
        response_idx = response_idx, ## mid sample
        extreme_idx = extreme_idx    ## end sample
    ))
}


#' Analyse fractional kinetics response time across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "response_time")`. Computes the fractional
#' response time for each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams response_time
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `response_time`, `response_value`,
#'   `fitted`, `idx`. Per-channel metadata are attached as attributes:
#'   - `"model"`: `NULL` (no parametric model is fitted).
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`, containing the baseline,
#'     response, and extreme key points.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [response_time()]
#'
#' @keywords internal
analyse_response_time <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    t0 = NULL,
    fraction = 0.5,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...
) {
    ## validation ==============================================
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    validate_mnirs_data(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, verbose)
    time_channel <- validate_time_channel(enquo(time_channel), data)
    t_vec <- data[[time_channel]]
    args <- list(...)

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(
            t0 = t0,
            fraction = fraction,
            direction = direction,
            end_fit_span = end_fit_span
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose
    )
    ## validate resolved args once, before fitting any channel
    per_channel <- validate_kinetics_args(per_channel, data, t_vec, verbose)

    ## method-specific fit: fractional response time (no model fit)
    response_time_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        response <- do.call(response_time, c(
            list(x = x_fit, t = t_fit),
            .a,
            list(verbose = verbose, bypass_checks = TRUE),
            args
        ))

        coefs <- data.frame(
            A              = response$A,
            B              = response$B,
            response_time  = response$response_time,
            response_value = response$response_value,
            fitted         = response$fitted,
            idx            = response$response_idx
        )

        ## bind baseline vec with `A`, and response and extreme scalars
        fitted_data <- data.frame(
            window_idx = c(
                response$baseline_idx,
                response$response_idx,
                response$extreme_idx
            ),
            fitted = c(
                rep(response$A, length(response$baseline_idx)),
                response$fitted,
                response$B
            )
        )
        ## omit NA idx
        fitted_data <- fitted_data[!is.na(fitted_data$window_idx), ]

        list(
            coefs = coefs,
            model = NULL,
            fitted_data = fitted_data,
            diag = compute_diagnostics(
                x        = x_fit[1L:3L], ## placeholder
                t        = t_fit[1L:3L], ## placeholder
                fitted   = c(coefs$A, coefs$fitted, coefs$B),
                n_params = 0L, ## invalid for response time method
                verbose  = verbose
            )
        )
    }

    return(analyse_kinetics_channels(
        data, nirs_channels, time_channel,
        per_channel, response_time_fit, verbose, args
    ))
}
