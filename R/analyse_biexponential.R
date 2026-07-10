#' Biexponential function
#'
#' Calculate a 5- or 6-parameter biexponential curve: a fast component that
#' drops to a rounded nadir, followed by a slow component that recovers
#' toward a stable plateau.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting value of the response
#'   variable (the `t = 0` intercept).
#' @param B1 A numeric parameter for the amplitude of the *fast* component;
#'   the depth of the initial drop to the nadir.
#' @param tau1 A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`. Dominates the initial steep fall.
#' @param B2 A numeric parameter for the amplitude of the *slow* component;
#'   the size of the recovery from the nadir toward the plateau.
#' @param tau2 A numeric parameter for the *slow* time constant (\eqn{\tau_2}),
#'   in units of the predictor variable `t`. Typically `tau2 >> tau1`.
#' @param TD A numeric parameter for the time delay before the onset of the
#'   response, in units of the predictor variable `t`. If `NULL` (*default*),
#'   a 5-parameter model without time delay is used.
#'
#' @details
#' This model family is fit by [analyse_kinetics()] when
#' `method = "biexponential"`, and by [stats::nls()] via the self-starting
#' wrapper [SSbiexponential()].
#'
#' ## Model equations
#'
#' 5-parameter model:
#'   `A - B1 * (1 - exp(-t / tau1)) + B2 * (1 - exp(-t / tau2))`
#'
#' 6-parameter model (with time delay):
#'   `ifelse(t <= TD, A, A - B1 * (1 - exp(-(t - TD) / tau1)) +
#'     B2 * (1 - exp(-(t - TD) / tau2)))`
#'
#' The fast component `B1`/`tau1` subtracts, driving the steep initial fall to
#' the nadir; the slow component `B2`/`tau2` adds, pulling the curve back up
#' once the fast term saturates. The two components carry independent
#' amplitudes and are not exchangeable.
#'
#' The response value as `t` approaches infinity is the *plateau*,
#' `A - B1 + B2`.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSbiexponential()], [monoexponential()]
#'
#' @examples
#' ## create a nadir-recovery biexponential curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B1, tau1, B2, tau2), data = data)
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
biexponential <- function(t, A, B1, tau1, B2, tau2, TD = NULL) {
    if (is.null(TD)) {
        ## 5-parameter: no time delay
        y <- A - B1 * (1 - exp(-t / tau1)) + B2 * (1 - exp(-t / tau2))
    } else {
        ## 6-parameter: with time delay
        y <- A - B1 * (1 - exp(-(t - TD) / tau1)) +
            B2 * (1 - exp(-(t - TD) / tau2))
        y[t < TD] <- A
    }
    return(y)
}


#' Initiate self-starting biexponential model
#'
#' [biexp_init()]: Returns initial values for the parameters in a `selfStart`
#' model.
#'
#' @param mCall A matched call to the function `model`.
#' @param data A data frame with time `t` and the response variable.
#' @param LHS The left-hand side expression of the model formula.
#' @param ... Additional arguments, including `fixed`, a named list of
#'   user-fixed parameter values from [init_fixed()] used to seed the
#'   remaining free estimates.
#'
#' @returns [biexp_init()]: Initial starting estimates for parameters in the
#'   model called by [SSbiexponential()].
#'
#' @keywords internal
biexp_init <- function(mCall, data, LHS, ...) {
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## check if TD parameter exists in the call
    has_TD <- "TD" %in% names(mCall)

    ## start and plateau from first and last ceiling(n/5) values
    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    plateau <- ab$B

    ## fast amplitude from the initial drop depth; slow amplitude closes the
    ## gap to the plateau (plateau = A - B1 + B2). Guard strictly positive
    B1_init <- fixed$B1 %||% max(A_init - min(x), .Machine$double.eps)
    B2_init <- fixed$B2 %||%
        max(plateau - A_init + B1_init, .Machine$double.eps)

    ## fast/slow time constants from the observed time span
    span <- diff(range(t))
    tau1_init <- fixed$tau1 %||% max(span / 10, .Machine$double.eps)
    tau2_init <- fixed$tau2 %||% max(span / 2, .Machine$double.eps)

    start <- c(
        A = A_init,
        B1 = B1_init,
        tau1 = tau1_init,
        B2 = B2_init,
        tau2 = tau2_init
    )

    if (has_TD) {
        ## estimate time delay from the steepest derivative changepoint
        dx_dt <- abs(diff(x) / diff(t))
        TD_init <- fixed$TD %||% max(t[which.max(dx_dt)] - tau1_init * 0.1, 0)
        start <- c(start, TD = TD_init)
    }

    return(start)
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. Supports both the
#' 5-parameter form (A, B1, tau1, B2, tau2) and the 6-parameter form
#' (A, B1, tau1, B2, tau2, TD); arity is inferred from the formula passed
#' to [stats::nls()].
#'
#' @usage
#' SSbiexponential(t, A, B1, tau1, B2, tau2, TD)
#'
#' @inheritParams biexponential
#'
#' @details
#' 5-parameter model: `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2)`
#'
#' 6-parameter model: `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, TD)`
#'
#' The 5-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSbiexponential(t, A = 0, B1, tau1, B2,
#'   tau2)` fixes the baseline at `A = 0`. Fixed parameters are excluded from
#'   estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [biexponential()], [stats::nls()], [stats::selfStart()],
#'   [SSmonoexponential()]
#'
#' @examples
#' ## create a nadir-recovery biexponential curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B1, tau1, B2, tau2), data = data)
#' summary(model)
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A = 70, B1, tau1, B2, tau2), data = data
#' )
#' summary(model_fixed)
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
SSbiexponential <- selfStart(
    model = biexponential,
    initial = init_fixed(
        biexp_init, c("A", "B1", "tau1", "B2", "tau2", "TD")
    ),
    parameters = c("A", "B1", "tau1", "B2", "tau2", "TD")
)


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' nadir-recovery curve to each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; default is `TRUE` to attempt to fit a
#'   6-parameter [SSbiexponential()] model (A, B1, tau1, B2, tau2, TD) with a
#'   time delay. If the 6-parameter fit fails, or if `use_TD = FALSE`, attempts
#'   to fit a reduced 5-parameter [SSbiexponential()] model
#'   (A, B1, tau1, B2, tau2).
#' @param fix An *optional* named list of model parameters (`A`, `B1`, `tau1`,
#'   `B2`, `tau2`, `TD`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied globally across
#'   `nirs_channels`. `TD` is fixable only when `use_TD = TRUE` for all
#'   channels; a fixed `TD` disables the 5-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B1`, `tau1`, `B2`, `tau2`, `TD`, `plateau`,
#'   `nadir_time`, `nadir_value`. Per-channel metadata are attached as
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
#' @seealso [analyse_kinetics()], [biexponential()], [SSbiexponential()]
#'
#' @keywords internal
analyse_biexponential <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    use_TD = FALSE,
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
    ## every channel fits the 6-parameter model
    use_TD_all <- all(vapply(per_channel, \(.a) .a$use_TD, logical(1)))
    fix <- validate_fix(
        fix,
        c("A", "B1", "tau1", "B2", "tau2", if (use_TD_all) "TD"),
        env = env
    )

    ## NA scaffold (method columns only) for convergence failure
    na_coefs <- data.frame(
        A = NA_real_,
        B1 = NA_real_,
        tau1 = NA_real_,
        B2 = NA_real_,
        tau2 = NA_real_,
        TD = NA_real_,
        plateau = NA_real_,
        nadir_time = NA_real_,
        nadir_value = NA_real_
    )

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, n_params, e, retry, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        msg <- c(
            "x" = "{n_params}-parameter {.fn SSbiexponential} fit failed \\
            for {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        )
        if (retry) {
            msg <- c(
                msg,
                "i" = "Attempting 5-parameter {.fn SSbiexponential} fit."
            )
        }
        cli_warn(msg, call = warn_call(env))
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting biexponential via nls
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        fit_data <- data.frame(.x = x_fit, .t = t_fit)
        params <- c("A", "B1", "tau1", "B2", "tau2", if (.a$use_TD) "TD")

        ## attempt nls fit; a failed 6-param fit falls back to the
        ## 5-param model unless TD is user-fixed
        retry <- .a$use_TD && !"TD" %in% names(fix)
        model <- tryCatch(
            nls(
                build_ss_formula(quote(SSbiexponential), params, fix),
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
                    build_ss_formula(quote(SSbiexponential), params, fix),
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
        fitted_vals <- stats::predict(model)

        ## time delay reported elapsed from onset (start_time)
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] - .a$start_time else NULL
        plateau_val <- coefs[["A"]] - coefs[["B1"]] + coefs[["B2"]]

        ## nadir: interior extremum where dy/dt = 0 has a closed form.
        ## s is elapsed since model onset (TD, or 0 for the 5-param model);
        ## valid only with a genuine two-component response (positive
        ## amplitudes, distinct taus, s > 0), else the curve is monotone and
        ## the extremum sits at the onset boundary
        s <- (log(coefs[["B2"]] / coefs[["tau2"]]) -
            log(coefs[["B1"]] / coefs[["tau1"]])) /
            (1 / coefs[["tau2"]] - 1 / coefs[["tau1"]])
        interior <- is.finite(s) && s > 0 &&
            coefs[["B1"]] > 0 && coefs[["B2"]] > 0

        ## nadir_time reported elapsed from start_time, mirroring
        ## MRT = TD + tau; adding TD_arg shifts the onset-relative s into the
        ## same frame, so t - TD in biexponential() recovers s exactly
        nadir_time_val <- sum(TD_arg, if (interior) s else 0)
        nadir_value_val <- biexponential(
            t = nadir_time_val,
            A = coefs[["A"]],
            B1 = coefs[["B1"]],
            tau1 = coefs[["tau1"]],
            B2 = coefs[["B2"]],
            tau2 = coefs[["tau2"]],
            TD = TD_arg
        )

        list(
            coefs = data.frame(
                A           = coefs[["A"]],
                B1          = coefs[["B1"]],
                tau1        = coefs[["tau1"]],
                B2          = coefs[["B2"]],
                tau2        = coefs[["tau2"]],
                TD          = TD_arg %||% NA_real_,
                plateau     = plateau_val,
                nadir_time  = nadir_time_val,
                nadir_value = nadir_value_val
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
        biexp_fit,
        verbose,
        interval_name,
        extra_args = c(args, list(fix = if (length(fix) > 0L) fix else NULL)),
        env = env
    ))
}
