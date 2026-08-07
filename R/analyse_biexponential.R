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
        y <- A -
            B1 * (1 - exp(-(t - TD) / tau1)) +
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

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## check if TD parameter exists in the call
    has_TD <- "TD" %in% names(mCall)

    ## time delay from the steepest derivative changepoint, in the same
    ## absolute frame as t; the fit window shift happens downstream
    TD_init <- fixed$TD %||% max(t[which.max(abs(diff(x) / diff(t)))], 0)

    ## profile the amplitudes out of a grid of time constants. deterministic,
    ## and far better seeded than an asymptote heuristic, which averages the
    ## first fifth of the window and so straddles the whole fast drop
    grid <- biexp_grid_start(x, t, TD = if (has_TD) TD_init else 0)

    if (is.null(grid)) {
        ## fall back to the asymptote heuristic when no grid point solves
        ab <- init_asymptotes(x)
        span <- max(diff(range(t)), .Machine$double.eps)
        grid <- list(
            A = ab$A,
            B1 = max(ab$A - min(x), .Machine$double.eps),
            tau1 = span / 10,
            B2 = max(ab$B - min(x), .Machine$double.eps),
            tau2 = span / 2
        )
    }

    ## user-fixed values take precedence over the grid optimum
    start <- c(
        A = fixed$A %||% grid$A,
        B1 = fixed$B1 %||% grid$B1,
        tau1 = fixed$tau1 %||% grid$tau1,
        B2 = fixed$B2 %||% grid$B2,
        tau2 = fixed$tau2 %||% grid$tau2
    )

    if (has_TD) {
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
    initial = init_fixed(biexp_init, c("A", "B1", "tau1", "B2", "tau2", "TD")),
    parameters = c("A", "B1", "tau1", "B2", "tau2", "TD")
)


#' Fit a biexponential in ratio parameterisation
#'
#' [fit_biexp_ratio()]: Fits [biexponential()] with the time constants
#' re-parameterised as `lt1 = log(tau1)` and `lr = log(tau2 / tau1)`.
#'
#' @details
#' [stats::nls()] does not converge on the natural `tau1`/`tau2` scale for
#' this model, failing even when started at the generating parameter values.
#' Estimating the log time constant and the log ratio instead decouples the
#' two components and bounds the ratio away from the degenerate `tau2 = tau1`
#' limit, where the design matrix is singular.
#'
#' Amplitudes are unbounded: NIRS responses may have negative amplitudes, so
#' identifiability rests on the time constants alone.
#'
#' A user-fixed `tau1` or `tau2` is substituted into the re-parameterised
#' formula rather than estimated: fixing `tau1` drops `lt1` and leaves `lr`
#' free about the fixed value; fixing `tau2` drops `lr` and caps `lt1` at
#' `log(tau2 / tau_ratio)`, which enforces the ratio bound from the other side.
#'
#' @param x,t Numeric vectors of the response and predictor variables.
#' @param params Character vector of parameter names in model order.
#' @param fix Named list of fixed parameter values.
#' @param tau_ratio A numeric lower bound on `tau2 / tau1`.
#' @param on_error A function called with the [stats::nls()] error condition.
#'
#' @returns An [nls][stats::nls] model in ratio parameterisation, or `NULL`.
#'
#' @keywords internal
fit_biexp_ratio <- function(x, t, params, fix, tau_ratio, on_error) {
    has_TD <- "TD" %in% params
    span <- diff(range(t))

    ## every failure route reports its condition and yields NULL
    fail <- \(.e) {
        on_error(.e)
        NULL
    }

    ## an under-determined fit is rejected before it reaches nls: the port
    ## algorithm can iterate indefinitely on such a problem rather than
    ## returning an error
    n_free <- length(setdiff(params, names(fix)))
    if (length(x) <= n_free) {
        return(fail(simpleError(sprintf(
            "%d observation%s for %d free parameters.",
            length(x),
            if (length(x) == 1L) "" else "s",
            n_free
        ))))
    }

    seed <- biexp_grid_start(x, t, tau_ratio = tau_ratio)
    if (is.null(seed)) {
        return(fail(simpleError(
            "No starting estimates could be resolved from the response."
        )))
    }

    ## TD is profiled separably: a coarse sweep locates the basin, then a
    ## refinement around it. nls polishes the value from there, so the
    ## profile only needs to be close
    TD_seed <- fix$TD %||% 0
    if (has_TD && is.null(fix$TD)) {
        profile <- \(.grid) {
            rss <- vapply(.grid, \(.td) {
                biexp_grid_start(x, t, tau_ratio, .td)$rss %||% NA_real_
            }, numeric(1))
            if (all(is.na(rss))) NULL else .grid[[which.min(rss)]]
        }
        coarse <- seq(0, span / 3, length.out = 12L)
        TD_seed <- profile(coarse) %||% 0
        step <- coarse[[2L]]
        TD_seed <- profile(
            seq(max(TD_seed - step, 0), TD_seed + step, length.out = 7L)
        ) %||%
            TD_seed
        seed <- biexp_grid_start(x, t, tau_ratio, TD_seed) %||% seed
    }

    ## a fixed tau2 caps tau1 from above, so the seed is pulled inside that
    ## cap before it is logged
    tau1_cap <- (fix$tau2 %||% Inf) / tau_ratio
    tau1_seed <- fix$tau1 %||% min(seed$tau1, tau1_cap / 1.01)
    tau2_seed <- fix$tau2 %||% seed$tau2

    ## start, lower, upper per parameter, in model order. the ratio starts
    ## strictly inside its bound; port fails on a start value sitting exactly
    ## on a boundary
    bounds <- list(
        A = c(seed$A, -Inf, Inf),
        B1 = c(seed$B1, -Inf, Inf),
        lt1 = c(log(tau1_seed), log(span / 1000), log(tau1_cap)),
        B2 = c(seed$B2, -Inf, Inf),
        lr = c(
            log(max(tau2_seed / tau1_seed, tau_ratio * 1.01)),
            log(tau_ratio),
            Inf
        ),
        TD = c(TD_seed, 0, span)
    )

    ## drop parameters the formula holds constant; a fixed time constant
    ## removes its log-scale counterpart
    free <- setdiff(
        names(bounds),
        c(
            names(fix),
            if (!is.null(fix$tau1)) "lt1",
            if (!is.null(fix$tau2)) "lr",
            if (!has_TD) "TD"
        )
    )
    ## port matches unnamed bounds positionally, so keep them aligned
    start <- vapply(bounds[free], `[[`, numeric(1), 1L)
    lower <- vapply(bounds[free], `[[`, numeric(1), 2L)
    upper <- vapply(bounds[free], `[[`, numeric(1), 3L)

    ## substitute fixed values into the formula in place of their symbols
    sub <- \(.p) fix[[.p]] %||% as.name(.p)
    tau1_expr <- fix$tau1 %||% quote(exp(lt1))
    tau2_expr <- fix$tau2 %||% substitute(T1 * exp(lr), list(T1 = tau1_expr))
    time_expr <- if (has_TD) {
        substitute(pmax(.t - TD, 0), list(TD = sub("TD")))
    } else {
        quote(.t)
    }
    rhs <- substitute(
        A - B1 * (1 - exp(-tt / t1)) + B2 * (1 - exp(-tt / t2)),
        list(
            A = sub("A"),
            B1 = sub("B1"),
            B2 = sub("B2"),
            t1 = tau1_expr,
            t2 = tau2_expr,
            tt = time_expr
        )
    )

    model <- tryCatch(
        nls(
            stats::as.formula(call("~", quote(.x), rhs)),
            data.frame(.x = x, .t = t),
            start = start,
            algorithm = "port",
            lower = lower,
            upper = upper,
            control = stats::nls.control(maxiter = 500L)
        ),
        error = fail
    )
    return(model)
}


#' Convert ratio-parameterised coefficients to the natural scale
#'
#' @param model An [nls][stats::nls] model in ratio parameterisation.
#' @param params Character vector of natural parameter names in model order.
#' @param fix Named list of fixed parameter values.
#'
#' @returns A named numeric vector on the natural `tau1`/`tau2` scale.
#'
#' @keywords internal
ratio_to_natural <- function(model, params, fix = list()) {
    coefs <- c(stats::coef(model), unlist(fix))
    ## a fixed time constant has no log-scale counterpart in the model
    coefs[["tau1"]] <- fix$tau1 %||% exp(coefs[["lt1"]])
    coefs[["tau2"]] <- fix$tau2 %||% (coefs[["tau1"]] * exp(coefs[["lr"]]))
    return(coefs[intersect(params, names(coefs))])
}


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' nadir-recovery curve to each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit a 6-parameter
#'   [SSbiexponential()] model (A, B1, tau1, B2, tau2, TD) with a time delay.
#'   If the 6-parameter fit fails, or if `use_TD = FALSE`, attempts to fit a
#'   reduced 5-parameter [SSbiexponential()] model (A, B1, tau1, B2, tau2).
#'   The user-facing default (`TRUE`) is set by [analyse_kinetics()].
#' @param fix An *optional* named list of model parameters (`A`, `B1`, `tau1`,
#'   `B2`, `tau2`, `TD`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied globally across
#'   `nirs_channels`. `TD` is fixable only when `use_TD = TRUE` for all
#'   channels; a fixed `TD` disables the 5-parameter fallback.
#' @param tau_ratio A numeric lower bound on the ratio of the slow to the fast
#'   time constant, `tau2 / tau1`; default is `2.5`. As `tau2` approaches
#'   `tau1` the two components become indistinguishable and the fit is
#'   singular, so the ratio is bounded away from that limit. The ratio is
#'   often only weakly identified and settles on this bound, in which case it
#'   sets the separation of the fast and slow components; lower values admit
#'   more similar time constants and larger, more strongly cancelling
#'   amplitudes.
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
    tau_ratio = 2.5,
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

    ## a ratio of 1 or less admits tau2 == tau1, where the two components
    ## are indistinguishable and the design is singular
    validate_numeric(
        tau_ratio, 1, c(1, Inf),
        inclusive = "right", msg1 = "one-element", env = env
    )

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
    # fmt: skip
    na_coefs <- as.data.frame(setNames(
        rep(list(NA_real_), 9L),
        c(
            "A", "B1", "tau1", "B2", "tau2", "TD", "plateau",
            "extremum_time", "extremum_value"
        )
    ))

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, n_params, e, retry, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        msg <- c(
            "x" = "{n_params}-parameter {.fn biexponential} fit failed \\
            for {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        )
        if (retry) {
            msg <- c(
                msg,
                "i" = "Attempting 5-parameter {.fn biexponential} fit."
            )
        }
        cli_warn(msg, call = warn_call(env))
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting biexponential via nls
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        params <- c("A", "B1", "tau1", "B2", "tau2", if (.a$use_TD) "TD")

        ## the 6-param model clamps to a flat `A` before `TD`, so the
        ## pre-onset baseline anchors `A`. the 5-param model has no such
        ## region and diverges at t < 0, so it is fit from `start_time` onward
        keep_rows <- function(.params) {
            if ("TD" %in% .params) rep(TRUE, length(t_fit)) else t_fit >= 0
        }

        ## attempt nls fit; a failed 6-param fit falls back to the
        ## 5-param model unless TD is user-fixed
        retry <- .a$use_TD && !"TD" %in% names(fix)
        attempt <- \(.params, .retry) {
            ## dropping TD narrows the window, so subset per attempt
            keep <- keep_rows(.params)
            fit_biexp_ratio(
                x_fit[keep],
                t_fit[keep],
                .params,
                fix,
                tau_ratio,
                \(e) {
                    fit_failed_warning(
                        .nirs,
                        length(.params),
                        e,
                        .retry,
                        verbose
                    )
                }
            )
        }
        model <- attempt(params, retry)
        if (is.null(model) && retry) {
            params <- setdiff(params, "TD")
            model <- attempt(params, FALSE)
        }

        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        coefs <- ratio_to_natural(model, params, fix)
        fitted_vals <- stats::predict(model)
        keep <- keep_rows(params)

        ## model parameters in biexponential() argument order
        pars <- as.list(coefs[c("A", "B1", "tau1", "B2", "tau2")])

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] else NULL

        ## nadir: interior extremum where dy/dt = 0 has a closed form, elapsed
        ## since model onset (TD, or 0 for the 5-param model). valid only with
        ## a genuine two-component response (positive amplitudes, distinct
        ## taus, root > 0), else the curve is monotone and the extremum sits
        ## at the onset boundary. amplitudes may be negative, so their sign is
        ## checked before taking any logarithm
        s <- with(pars, {
            root <- if (B1 > 0 && B2 > 0) {
                (log(B2 / tau2) - log(B1 / tau1)) / (1 / tau2 - 1 / tau1)
            } else {
                NA_real_
            }
            if (is.finite(root) && root > 0) root else 0
        })

        ## nadir_time reported elapsed from start_time, mirroring
        ## MRT = TD + tau; adding TD_arg shifts the onset-relative s into the
        ## same frame, so t - TD in biexponential() recovers s exactly
        nadir_time_val <- sum(TD_arg, s)
        nadir_value_val <- do.call(
            biexponential,
            c(list(t = nadir_time_val), pars, list(TD = TD_arg))
        )

        list(
            coefs = data.frame(
                pars,
                TD = TD_arg %||% NA_real_,
                plateau = pars$A - pars$B1 + pars$B2,
                nadir_time = nadir_time_val,
                nadir_value = nadir_value_val
            ),
            model = model,
            fitted_data = data.frame(
                window_idx = valid$idx[keep],
                fitted = fitted_vals
            ),
            diag = compute_diagnostics(
                x_fit[keep],
                t_fit[keep],
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
