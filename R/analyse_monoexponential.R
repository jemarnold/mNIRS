#' Monoexponential function
#'
#' Calculate a 3- or 4-parameter monoexponential curve. This model family is
#' fit by [analyse_kinetics()] when `method = "monoexponential"`, and by
#' [stats::nls()] via the self-starting wrapper [SSmonoexponential()].
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
#' @seealso [analyse_kinetics()], [SSmonoexponential()], [response_time()],
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
#' model <- nls(x ~ SSmonoexponential(t, A, B, tau, TD), data = data)
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
#' @param ... Additional arguments, including `fixed`, a named list of
#'   user-fixed parameter values from [init_fixed()] used to seed the
#'   remaining free estimates.
#'
#' @returns [monoexp_init()]: Initial starting estimates for parameters in the
#'   model called by [SSmonoexponential()].
#'
#' @keywords internal
monoexp_init <- function(mCall, data, LHS, ...) {
    ## self-start parameters for nls of monoexponential fit function
    ## uses base R `SSasymp()` initialisation approach
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## check if TD parameter exists in the call
    has_TD <- "TD" %in% names(mCall)

    ## fit linear model to log-transformed differences from estimated asymptote
    ## asymptotes from 1/5 response signal
    n_asymp <- max(1, ceiling(n / 5))
    A_init <- fixed$A %||% mean(x[seq_len(n_asymp)])
    B_init <- fixed$B %||% mean(x[seq(n - n_asymp + 1, n)])

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

    if (sum(x_pos) >= 3) {
        x_shifted[!x_pos] <- min(x_shifted[x_pos]) / 2
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

    tau_init <- fixed$tau %||% max(tau_init, .Machine$double.eps)

    if (has_TD) {
        ## 4-parameter: estimate time delay from derivative changepoint
        dx_dt <- abs(diff(x) / diff(t))
        td_idx <- which.max(dx_dt)
        TD_init <- fixed$TD %||% max(t[td_idx] - tau_init * 0.1, 0)
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
#' SSmonoexponential(t, A, B, tau, TD)
#'
#' @inheritParams monoexponential
#'
#' @details
#' 3-parameter model: `x ~ SSmonoexponential(t, A, B, tau)`
#'
#' 4-parameter model: `x ~ SSmonoexponential(t, A, B, tau, TD)`
#'
#' The 3-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSmonoexponential(t, A = 0, B, tau)` fixes the
#'   baseline at `A = 0`. Fixed parameters are excluded from estimation and
#'   are not returned by [stats::coef()].
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
#' model4 <- nls(x ~ SSmonoexponential(t, A, B, tau, TD), data = data)
#' model4
#'
#' ## 3-parameter fit on the same data
#' model3 <- nls(x ~ SSmonoexponential(t, A, B, tau), data = data)
#' model3
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(x ~ SSmonoexponential(t, A = 10, B, tau, TD), data = data)
#' coef(model_fixed)
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
SSmonoexponential <- selfStart(
    model = monoexponential,
    initial = init_fixed(monoexp_init, c("A", "B", "tau", "TD")),
    parameters = c("A", "B", "tau", "TD")
)


#' Analyse monoexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "monoexponential")`. Fits a monoexponential
#' curve to each `nirs_channel` within a single *"mnirs"* data frame. See
#' [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; default is `TRUE` to attempt to fit a
#'   4-parameter [SSmonoexponential()] model (A, B, tau, TD) with a time delay.
#'   If the 4-parameter fit fails, or if `use_TD = FALSE`, attempts to
#'   fit a reduced 3-parameter [SSmonoexponential()] model (A, B, tau).
#' @param fix An *optional* named list of model parameters to hold
#'   constant during fitting, e.g. `fix = list(A = 0)`. Fixed parameters
#'   are excluded from estimation and reported at their fixed values.
#'   Applied globally across `nirs_channels`. `TD` is fixable only when
#'   `use_TD = TRUE` for all channels; a fixed `TD` disables the
#'   3-parameter fallback.
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
#' @seealso [analyse_kinetics()], [monoexponential()], [SSmonoexponential()]
#'
#' @keywords internal
analyse_monoexponential <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    use_TD = TRUE,
    fix = NULL,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    ## validation ==================================================
    args <- list(...)
    ## interval label; falls back to the `data` argument name when unsupplied
    interval_name <- args$interval_name %||% deparse(substitute(data))

    ## shared prologue: validate data, resolve channels/time, broadcast and
    ## validate per-channel args
    setup <- setup_kinetics_worker(
        data,
        enquo(nirs_channels),
        enquo(time_channel),
        arg_list = list(
            use_TD = use_TD,
            start_time = start_time,
            direction = direction,
            end_window = end_window
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose,
        env = env
    )
    nirs_channels <- setup$nirs_channels
    time_channel <- setup$time_channel
    per_channel <- setup$per_channel

    ## global fixed parameters bypass per-channel resolution: a named
    ## list would be misread as a channel map. TD is only fixable when
    ## every channel fits the 4-parameter model
    use_TD_all <- all(vapply(per_channel, \(.a) .a$use_TD, logical(1)))
    fix <- validate_fix(
        fix,
        c("A", "B", "tau", if (use_TD_all) "TD"),
        env = env
    )

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
    fit_failed_warning <- function(.nirs, n_params, e, retry, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        msg <- c(
            "x" = "{n_params}-parameter {.fn SSmonoexponential} fit failed \\
            for {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        )
        if (retry) {
            msg <- c(
                msg,
                "i" = "Attempting 3-parameter {.fn SSmonoexponential} fit."
            )
        }
        cli_warn(msg, call = warn_call(env))
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting monoexponential via nls
    monoexponential_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        fit_data <- data.frame(.x = x_fit, .t = t_fit)
        params <- c("A", "B", "tau", if (.a$use_TD) "TD")

        ## attempt nls fit; a failed 4-param fit falls back to the
        ## 3-param model unless TD is user-fixed
        retry <- .a$use_TD && !"TD" %in% names(fix)
        model <- tryCatch(
            nls(
                build_ss_formula(quote(SSmonoexponential), params, fix),
                fit_data
            ),
            error = \(e) {
                fit_failed_warning(.nirs, length(params), e, retry, verbose)
                NULL
            }
        )
        if (is.null(model) && retry) {
            params <- setdiff(params, "TD")
            model <- tryCatch(
                nls(
                    build_ss_formula(quote(SSmonoexponential), params, fix),
                    fit_data
                ),
                error = \(e) {
                    fit_failed_warning(.nirs, length(params), e, FALSE, verbose)
                    NULL
                }
            )
        }

        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        coefs <- full_coefs(model, params, fix)

        ## enforce direction: bounded refit on D = B - A when inverted
        free_extra <- setdiff(params, c("A", "B", names(fix)))
        enforced <- enforce_direction(
            model,
            coefs,
            fit_data,
            direction = .a$direction,
            amp_fn = quote(monoexponential),
            extra = coefs[free_extra],
            ## data-scaled tau floor: tau pinned here is a degenerate
            ## step fit, not a genuine response
            extra_lower = if ("tau" %in% free_extra) {
                c(tau = diff(range(t_fit)) * 1e-6)
            },
            fn = quote(SSmonoexponential),
            fix = fix,
            .nirs = .nirs,
            interval_name = interval_name,
            verbose = verbose,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_coefs))
        }
        model <- enforced$model
        coefs <- enforced$coefs
        fitted_vals <- stats::predict(model)

        TD_arg <- if ("TD" %in% params) coefs[["TD"]] - .a$start_time else NULL
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
                x_fit,
                t_fit,
                fitted_vals,
                n_params = length(stats::coef(model)),
                verbose,
                env
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
        extra_args = c(args, list(fix = if (length(fix) > 0L) fix else NULL)),
        env = env
    ))
}
