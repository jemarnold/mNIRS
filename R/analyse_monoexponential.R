#' Monoexponential function
#'
#' Calculate a 3- or 4-parameter monoexponential curve. This model family is
#' fit by [analyse_kinetics()] when `method = "monoexponential"`, and by
#' [stats::nls()] via the self-starting wrapper [SSmonoexp()].
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting (baseline) value of the
#'   response variable.
#' @param B A numeric parameter for the ending (asymptote) value of the
#'   response variable.
#' @param tau A numeric parameter for the time constant `tau` (\eqn{\tau})
#'   of the exponential curve, in units of the predictor variable `t`.
#' @param TD A numeric parameter for the time delay before the onset of
#'   exponential response, in units of the predictor variable `t`. If `NULL`
#'   (*default*), a 3-parameter model without time delay is used.
#'
#' @details
#' ## Model equations
#'
#' 3-parameter model:
#'   `A + (B - A) * (1 - exp(-t / tau))`
#'
#' 4-parameter model:
#'   `ifelse(t <= TD, A, A + (B - A) * (1 - exp(-(t - TD) / tau)))`
#'
#' `tau` is the time constant, equal to the reciprocal of the rate constant
#' `k` (`k = 1 / tau` in reciprocal units of `time_channel`; i.e. `sec^-1s`).
#' Common derived quantities include the mean response time `MRT = TD + tau`
#' and the half-response time `HRT = TD + tau * log(2)`.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSmonoexp()], [response_time()],
#'   [peak_slope()]
#'
#' @examples
#' ## create an exponential curve with random noise
#' set.seed(13)
#' t <- 1:60
#' x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
#'     rnorm(length(t), 0, 3)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
#' model
#'
#' y <- predict(model, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y))
#'     }
#' }
#'
#' @export
monoexponential <- function(t, A, B, tau, TD = NULL) {
    if (is.null(TD)) {
        ## 3-parameter: no time delay
        y <- A + (B - A) * (1 - exp(-t / tau))
    } else {
        ## 4-parameter: with time delay
        y <- A + (B - A) * (1 - exp(-(t - TD) / tau))
        y[t < TD] <- A
    }
    return(y)
}


#' Initiate self-starting monoexponential model
#'
#' [monoexp_init()]: Returns initial values for the parameters in a `selfStart`
#' model.
#'
#' @param mCall A matched call to the function `model`.
#' @param data A data frame with time `t` and the response variable.
#' @param LHS The left-hand side expression of the model formula.
#' @param ... Additional arguments.
#'
#' @returns [monoexp_init()]: Initial starting estimates for parameters in the
#'   model called by [SSmonoexp()].
#'
#' @keywords internal
monoexp_init <- function(mCall, data, LHS, ...) {
    ## self-start parameters for nls of monoexponential fit function
    ## uses base R `SSasymp()` initialisation approach
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## check if TD parameter exists in the call
    has_TD <- "TD" %in% names(mCall)

    ## fit linear model to log-transformed differences from estimated asymptote
    ## asymptotes from 1/5 response signal
    n_asymp <- max(1, ceiling(n / 5))
    A_init <- mean(x[seq_len(n_asymp)])
    B_init <- mean(x[seq(n - n_asymp + 1, n)])

    ## estimate rate constant via linearisation (SSasymp method)
    ## sign +ve for upward (B > A); -ve for downward (B < A).
    ## linearised: log(sign * (B - x)) ~ log(B - A) - t/tau,
    #! doesn't work as well on current tests as current implementation
    # x_shifted <- sign(B_init - A_init) * (B_init - x)

    ## estimate rate constant via linearisation (SSasymp method)
    ## log(B - x) ~ log(B - A) - t/tau
    ## use shifted x to avoid log of negative/zero
    x_shifted <- B_init - x
    x_pos <- x_shifted > 0
    x_shifted[!x_pos] <- min(x_shifted[x_pos]) / 2

    if (sum(x_pos) >= 3) {
        lm_fit <- stats::lm(log(x_shifted) ~ t)
        rate <- -coef(lm_fit)[2L]
        tau_init <- if (is.finite(rate) && rate > 0) {
            1 / rate
        } else {
            diff(range(t)) / 3
        }
    } else {
        ## fallback: tau from 63.2% rise point (sign agnostic)
        target <- A_init + 0.632 * (B_init - A_init)
        tau_init <- t[which.min(abs(x - target))]
        tau_init <- max(tau_init, diff(range(t)) / 10)
    }

    tau_init <- max(tau_init, .Machine$double.eps)

    if (has_TD) {
        ## 4-parameter: estimate time delay from derivative changepoint
        dx_dt <- abs(diff(x) / diff(t))
        td_idx <- which.max(dx_dt)
        TD_init <- max(t[td_idx] - tau_init * 0.1, 0)
        return(c(A = A_init, B = B_init, tau = tau_init, TD = TD_init))
    } else {
        ## 3-parameter: no time delay
        return(c(A = A_init, B = B_init, tau = tau_init))
    }
}


#' Self-starting monoexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [monoexponential()], for use with [stats::nls()]. Supports both the
#' 3-parameter form (A, B, tau) and the 4-parameter form (A, B, tau, TD);
#' arity is inferred from the formula passed to [stats::nls()].
#'
#' @usage
#' SSmonoexp(t, A, B, tau, TD)
#'
#' @inheritParams monoexponential
#'
#' @details
#' 3-parameter model: `x ~ SSmonoexp(t, A, B, tau)`
#'
#' 4-parameter model: `x ~ SSmonoexp(t, A, B, tau, TD)`
#'
#' The 3-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()] 
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [monoexponential()], [stats::nls()], [stats::selfStart()],
#'   [stats::SSasymp()]
#'
#' @examples
#' ## create an exponential curve with random noise
#' set.seed(13)
#' t <- 1:60
#' x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
#'     rnorm(length(t), 0, 3)
#' data <- data.frame(t, x)
#'
#' ## 4-parameter fit
#' model4 <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
#' model4
#'
#' ## 3-parameter fit on the same data
#' model3 <- nls(x ~ SSmonoexp(t, A, B, tau), data = data)
#' model3
#'
#' y4 <- predict(model4, data)
#' y3 <- predict(model3, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y4, colour = "4-param")) + 
#'             ggplot2::geom_line(ggplot2::aes(y = y3, colour = "3-param"))
#'     }
#' }
#'
#' @export
SSmonoexp <- selfStart(
    model = monoexponential,
    initial = monoexp_init,
    parameters = c("A", "B", "tau", "TD")
)


#' Analyse monoexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "monoexponential")`. Fits a monoexponential
#' curve to each `nirs_channel` within a single *"mnirs"* data frame. See
#' [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_time_delay Logical; default is `TRUE` to attempt to fit a
#'   4-parameter [SSmonoexp()] model (A, B, tau, TD) with a time delay.
#'   If the 4-parameter fit fails, or if `use_time_delay = FALSE`, attempts to
#'   fit a reduced 3-parameter [SSmonoexp()] model (A, B, tau).
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `tau`, `k`, `TD`, `MRT`, `HRT`, `tau_fitted`,
#'   `MRT_fitted`, `HRT_fitted`. Per-channel metadata are attached as
#'   attributes:
#'   - `"model"`: an [nls][stats::nls] model object, or `NULL` for channels
#'     where fitting failed.
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [monoexponential()], [SSmonoexp()]
#'
#' @keywords internal
analyse_monoexponential <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    use_time_delay = TRUE, #! better arg name?
    t0 = NULL, #! better arg name?
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    verbose = TRUE,
    ...
) {
    ## validation ==================================================
    validate_mnirs_data(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, verbose)
    time_channel <- validate_time_channel(enquo(time_channel), data)
    t_vec <- data[[time_channel]]
    args <- list(...)
    ## interval label; falls back to the `data` argument name when unsupplied
    interval_name <- args$interval_name %||% deparse(substitute(data))

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(
            use_time_delay = use_time_delay,
            t0 = t0,
            direction = direction,
            end_fit_span = end_fit_span
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose
    )
    ## validate resolved args once, before fitting any channel
    per_channel <- validate_kinetics_args(per_channel, data, t_vec, verbose)

    ## NA scaffold (method columns only) for convergence failure
    na_coefs <- data.frame(
        A = NA_real_,
        B = NA_real_,
        tau = NA_real_,
        k = NA_real_,
        TD = NA_real_,
        MRT = NA_real_,
        HRT = NA_real_,
        tau_fitted = NA_real_,
        MRT_fitted = NA_real_,
        HRT_fitted = NA_real_
    )

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, n_params, e, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        msg <- c(
            "x" = "{n_params}-parameter {.fn SSmonoexp} fit failed for \\
            {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        )
        if (n_params == 4L) {
            msg <- c(msg, "i" = "Attempting 3-parameter {.fn SSmonoexp} fit.")
        }
        cli_warn(msg)
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting monoexponential via nls
    monoexponential_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        ## derive n_params from use_time_delay for internal use
        n_params <- if (.a$use_time_delay) 4L else 3L
        fit_data <- data.frame(.x = x_fit, .t = t_fit)

        ## attempt nls fit on 4-param then fall back to 3-param on failure
        model <- NULL
        if (n_params == 4L) {
            model <- tryCatch(
                nls(.x ~ SSmonoexp(.t, A, B, tau, TD), fit_data),
                error = \(e) {
                    fit_failed_warning(.nirs, n_params, e, verbose)
                    NULL
                }
            )
            if (is.null(model)) n_params <- 3L
        }

        if (n_params == 3L) {
            model <- tryCatch(
                nls(.x ~ SSmonoexp(.t, A, B, tau), fit_data),
                error = \(e) {
                    fit_failed_warning(.nirs, n_params, e, verbose)
                    NULL
                }
            )
        }

        ## TODO: implement fallback HRT method on convergence failure
        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        fitted_vals <- stats::predict(model)
        coefs <- stats::coef(model)
        TD_arg <- if (n_params == 4L) coefs[["TD"]] - .a$t0 else NULL
        TD_val <- TD_arg %||% NA_real_
        MRT_val <- sum(TD_arg, coefs[["tau"]])
        HRT_val <- sum(TD_arg, coefs[["tau"]] * log(2))

        ## predict response at tau, MRT, and HRT using the fitted model
        fitted_params <- monoexponential(
            t = c(coefs[["tau"]], MRT_val, HRT_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            TD = TD_arg
        )

        list(
            coefs = data.frame(
                A          = coefs[["A"]],
                B          = coefs[["B"]],
                tau        = coefs[["tau"]],
                k          = 1 / coefs[["tau"]], ## time_channel units^-1
                TD         = TD_val,
                MRT        = MRT_val,
                HRT        = HRT_val,
                tau_fitted = fitted_params[[1L]],
                MRT_fitted = fitted_params[[2L]],
                HRT_fitted = fitted_params[[3L]]
            ),
            model = model,
            fitted_data = data.frame(
                window_idx = valid$idx,
                fitted     = fitted_vals
            ),
            diag = compute_diagnostics(
                x_fit, t_fit, fitted_vals, n_params, verbose
            )
        )
    }

    return(analyse_kinetics_channels(
        data,
        nirs_channels,
        time_channel,
        per_channel,
        monoexponential_fit,
        verbose,
        interval_name,
        extra_args = args
    ))
}