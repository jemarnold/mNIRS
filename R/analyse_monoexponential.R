#' Monoexponential function
#'
#' Calculate a 3- or 4-parameter monoexponential curve.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting baseline of the response
#'   variable.
#' @param B A numeric parameter for the ending asymptote of the response
#'   variable.
#' @param tau A numeric parameter for the time constant `tau` (\eqn{\tau})
#'   of the exponential curve, in units of the predictor variable `t`.
#' @param TD A numeric parameter for the time delay before the onset of
#'   exponential response, in units of the predictor variable `t`. If `NULL`
#'   (*default*), a 3-parameter model without time delay is used.
#'
#' @details
#' This model family is
#' fit by [analyse_kinetics()] when `method = "monoexponential"`, and by
#' [stats::nls()] via the self-starting wrapper [SSmonoexponential()].
#'
#' ## Model equations
#'
#' 3-parameter model: `A + (B - A) * (1 - exp(-t / tau))`
#'
#' 4-parameter model: `A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau))`
#'
#' Clamping the shifted time at zero holds the curve at the baseline `A` until
#' the onset of the response at `t = TD`.
#'
#' The rate constant `k` is the reciprocal of `tau` (`k = 1 / tau`) in
#' reciprocal units of `time_channel`; i.e. `sec^-1s`).
#'
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
#' summary(model)
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
        y <- A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau))
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
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]

    ## user-fixed parameter values constrain the free estimates
    fixed <- list(...)$fixed %||% list()

    ## check if TD parameter exists in the call
    has_TD <- "TD" %in% names(mCall)

    ## profile tau (and TD for the 4-parameter model) on a coarse grid
    ## and keep the RSS-minimising start (cf. `biexp_grid_start()`).
    ## the model is linear in A and B once tau and TD are held fixed,
    ## so the asymptotes are solved by least squares at every grid
    ## point. point estimates from derivative changepoints or
    ## log-linearisation are too sensitive to noise, overshoot, and
    ## plateau data on real NIRS signals, and can strand nls with a
    ## singular gradient
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        span <- 1
    }
    tau_grid <- fixed$tau %||%
        exp(seq(log(span / 100), log(span), length.out = 25L))
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, 0.5 * span, length.out = 21L)
    }

    ## y = B + (A - B) * e, so free A and B come from simple regression
    ## of x on e, and a single free asymptote from projection on its
    ## basis; a degenerate grid point (e.g. all points pre-onset) yields
    ## non-finite estimates and is discarded via infinite rss
    solve_AB <- function(e) {
        A <- fixed$A
        B <- fixed$B
        if (is.null(A) && is.null(B)) {
            ec <- e - mean(e)
            slope <- sum(ec * x) / sum(ec^2)
            B <- mean(x) - slope * mean(e)
            A <- B + slope
        } else if (is.null(A)) {
            A <- sum(e * (x - B * (1 - e))) / sum(e^2)
        } else if (is.null(B)) {
            B <- sum((1 - e) * (x - A * e)) / sum((1 - e)^2)
        }
        rss <- if (is.finite(A + B)) sum((x - B - (A - B) * e)^2) else Inf
        return(c(A = A, B = B, rss = rss))
    }

    grid <- expand.grid(tau = tau_grid, TD = td_grid)
    fits <- vapply(seq_len(nrow(grid)), \(.i) {
        solve_AB(exp(-(if (has_TD) pmax(t - grid$TD[.i], 0) else t) /
            grid$tau[.i]))
    }, numeric(3L))
    best <- which.min(fits["rss", ])

    return(c(
        A = fits[["A", best]],
        B = fits[["B", best]],
        tau = grid$tau[best],
        TD = if (has_TD) grid$TD[best]
    ))
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
#'   name in the formula, e.g. `x ~ SSmonoexponential(t, A = 0, B, tau)` fixes
#'   the baseline at `A = 0`. Fixed parameters are excluded from estimation
#'   and are not returned by [stats::coef()].
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
#' summary(model4)
#'
#' ## 3-parameter fit on the same data
#' model3 <- nls(x ~ SSmonoexponential(t, A, B, tau), data = data)
#' summary(model3)
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(x ~ SSmonoexponential(t, A = 10, B, tau, TD), data = data)
#' summary(model_fixed)
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
#'   Applied to every channel, or per-channel as a list of lists keyed by
#'   channel name, e.g. `fix = list(smo2 = list(A = 0))`. `TD` is fixable
#'   for channels where `use_TD = TRUE`; a fixed `TD` disables the
#'   3-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `tau`, `k`, `TD`, `MRT`, `HRT`, `MRT_fitted`,
#'   `HRT_fitted`. Per-channel metadata are attached as
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
        # fmt: skip
        arg_list = mget(
            c("use_TD", "fix", "start_time", "direction", "end_window")
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        ## TD is only fixable where that channel fits the 4-parameter model
        fix_params = \(.a) c("A", "B", "tau", if (.a$use_TD) "TD"),
        verbose = verbose,
        env = env
    )
    ## NA scaffold (method columns only) for convergence failure
    # fmt: skip
    na_cols <- c(
        "A", "B", "tau", "k", "TD", "MRT", "HRT", "MRT_fitted", "HRT_fitted"
    )

    ## method-specific fit: self-starting monoexponential via nls; a
    ## failed 4-param fit falls back to the 3-param model
    monoexponential_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        fit <- fit_td_fallback(
            x_fit,
            t_fit,
            params = c("A", "B", "tau", if (.a$use_TD) "TD"),
            .a,
            fitter = \(.data, .params, on_error) {
                formula <- build_ss_formula(
                    quote(SSmonoexponential),
                    .params,
                    .a$fix
                )
                tryCatch(nls(formula, .data), error = on_error)
            },
            fn = quote(SSmonoexponential),
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(fit$model)) {
            return(build_na_results(na_cols))
        }
        params <- fit$params
        coefs <- full_coefs(fit$model, params, .a$fix)

        ## enforce direction: bounded refit on D = B - A when inverted
        free_extra <- setdiff(params, c("A", "B", names(.a$fix)))
        enforced <- enforce_direction(
            fit$model,
            coefs,
            fit$data,
            direction = .a$direction,
            amp_fn = quote(monoexponential),
            extra = coefs[free_extra],
            ## data-scaled tau floor: tau pinned here is a degenerate
            ## step fit, not a genuine response
            extra_lower = if ("tau" %in% free_extra) {
                c(tau = diff(range(t_fit)) * 1e-6)
            },
            fn = quote(SSmonoexponential),
            fix = .a$fix,
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_cols))
        }
        coefs <- enforced$coefs

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] else NULL
        MRT_val <- sum(TD_arg, coefs[["tau"]])
        HRT_val <- sum(TD_arg, coefs[["tau"]] * log(2))

        ## predict response at tau, MRT, and HRT using the fitted model;
        ## tau shifted by TD_arg so all time points share the reported frame
        fitted_params <- monoexponential(
            t = c(MRT_val, HRT_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            TD = TD_arg
        )

        build_fit_results(
            data.frame(
                A = coefs[["A"]],
                B = coefs[["B"]],
                tau = coefs[["tau"]],
                k = 1 / coefs[["tau"]], ## time_channel units^-1
                TD = TD_arg %||% NA_real_,
                MRT = MRT_val,
                HRT = HRT_val,
                MRT_fitted = fitted_params[[1L]],
                HRT_fitted = fitted_params[[2L]]
            ),
            enforced$model,
            x_fit,
            t_fit,
            valid,
            fit$keep,
            env
        )
    }

    return(analyse_kinetics_channels(
        data,
        setup$nirs_channels,
        setup$time_channel,
        setup$per_channel,
        monoexponential_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
