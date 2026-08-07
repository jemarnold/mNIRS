#' Biexponential function
#'
#' Calculate a 5- or 6-parameter biexponential curve: a fast component that
#' carries the response to a rounded turning point, followed by a slow
#' component that recovers toward a stable plateau.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting value of the response
#'   variable (the `t = 0` intercept).
#' @param B1 A numeric parameter for the amplitude of the *fast* component;
#'   the size of the initial excursion to the turning point.
#' @param tau1 A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`. Dominates the initial steep fall.
#' @param B2 A numeric parameter for the amplitude of the *slow* component;
#'   the size of the recovery from the turning point toward the plateau.
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
#' 6-parameter model (with time delay), where `ts = pmax(t - TD, 0)`:
#'   `A - B1 * (1 - exp(-ts / tau1)) + B2 * (1 - exp(-ts / tau2))`
#'
#' Clamping the shifted time at zero holds the curve at the baseline `A` until
#' the onset of the response at `t = TD`.
#'
#' The fast component `B1`/`tau1` subtracts, driving the steep initial
#' response to the excursion point; the slow component `B2`/`tau2` adds, pulling
#' the curve back once the fast term saturates. Amplitudes may be negative, in
#' which case both response components are mirrored and the turning point is
#' a maximum rather than a minimum. The two components carry independent
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
#' ## create a biexponential excursion-recovery curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B1, lt1, B2, lr), data = data)
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
        ts <- pmax(t - TD, 0)
        y <- A - B1 * (1 - exp(-ts / tau1)) + B2 * (1 - exp(-ts / tau2))
    }
    return(y)
}


#' Biexponential function in log-ratio parameterisation
#'
#' [biexponential()] with the time constants re-parameterised as
#' `lt1 = log(tau1)` and `lr = log(tau2 / tau1)`, the scale on which
#' [SSbiexponential()] is estimated.
#'
#' @inheritParams biexponential
#' @param lt1 A numeric parameter for the log *fast* time constant.
#' @param lr A numeric parameter for the log ratio of the *slow* to the
#'   *fast* time constant.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @keywords internal
biexponential_ratio <- function(t, A, B1, lt1, B2, lr, TD = NULL) {
    tau1 <- exp(lt1)
    return(biexponential(t, A, B1, tau1, B2, tau1 * exp(lr), TD))
}


#' Profile the biexponential time delay by separable least squares
#'
#' A coarse sweep of candidate delays locates the RSS basin, then a
#' refinement around it. [stats::nls()] polishes the value from there, so
#' the profile only needs to be close.
#'
#' @param x,t Numeric vectors of the response and predictor variables.
#' @param tau_ratio A numeric lower bound on `tau2 / tau1`.
#' @param span Numeric range of `t`.
#'
#' @returns A numeric time delay seed, or `NULL` when no candidate solves.
#'
#' @keywords internal
biexp_profile_TD <- function(x, t, tau_ratio, span) {
    profile <- \(.grid) {
        rss <- vapply(.grid, \(.td) {
            biexp_grid_start(x, t, tau_ratio, .td)$rss %||% NA_real_
        }, numeric(1))
        if (all(is.na(rss))) NULL else .grid[[which.min(rss)]]
    }
    coarse <- seq(0, span / 3, length.out = 12L)
    TD <- profile(coarse)
    if (is.null(TD)) {
        return(NULL)
    }
    step <- coarse[[2L]]
    fine <- profile(seq(max(TD - step, 0), TD + step, length.out = 7L))
    return(fine %||% TD)
}


#' Initiate self-starting biexponential model
#'
#' [biexp_init()]: Returns initial values for the parameters in a `selfStart`
#' model, with the time constants on the log-ratio scale (`lt1`, `lr`).
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
    has_TD <- "TD" %in% names(mCall)
    span <- diff(range(t))

    ## TD profiled by the same separable RSS sweep as the package fit path
    TD_init <- fixed$TD %||%
        if (has_TD) biexp_profile_TD(x, t, 2.5, span) %||% 0 else 0

    ## profile the amplitudes out of a grid of time constants
    grid <- biexp_grid_start(x, t, TD = TD_init)
    if (is.null(grid)) {
        stop("No starting estimates could be resolved from the response.")
    }

    ## user-fixed values take precedence over the grid optimum. the grid is
    ## confined to tau2 >= 2.5 * tau1, so lr starts at or above log(2.5); the
    ## nudge keeps it strictly inside that bound for callers fitting with
    ## algorithm = "port", which errors on a start sitting on a boundary
    start <- c(
        A = fixed$A %||% grid$A,
        B1 = fixed$B1 %||% grid$B1,
        lt1 = fixed$lt1 %||% log(grid$tau1),
        B2 = fixed$B2 %||% grid$B2,
        lr = fixed$lr %||% max(log(grid$tau2 / grid$tau1), log(2.5 * 1.01))
    )

    if (has_TD) {
        start <- c(start, TD = TD_init)
    }

    return(start)
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. The time constants are
#' estimated on the log-ratio scale, `lt1 = log(tau1)` and
#' `lr = log(tau2 / tau1)`, which expresses the `tau2 / tau1` separation of
#' the two components as a simple bound on `lr`. Supports both the
#' 5-parameter form (A, B1, lt1, B2, lr) and the 6-parameter form
#' (A, B1, lt1, B2, lr, TD); arity is inferred from the formula passed
#' to [stats::nls()].
#'
#' @usage
#' SSbiexponential(t, A, B1, lt1, B2, lr, TD)
#'
#' @inheritParams biexponential_ratio
#'
#' @details
#' 5-parameter model: `x ~ SSbiexponential(t, A, B1, lt1, B2, lr)`
#'
#' 6-parameter model: `x ~ SSbiexponential(t, A, B1, lt1, B2, lr, TD)`
#'
#' The 5-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' [stats::coef()] returns the log-scale `lt1` and `lr`; recover the natural
#'   time constants as `tau1 = exp(lt1)` and `tau2 = exp(lt1 + lr)`.
#'
#' [analyse_kinetics()] fits this model with `algorithm = "port"` and box
#'   bounds.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSbiexponential(t, A = 0, B1, lt1, B2,
#'   lr)` fixes the baseline at `A = 0`. Fixed parameters are excluded from
#'   estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [biexponential()], [stats::nls()], [stats::selfStart()],
#'   [SSmonoexponential()]
#'
#' @examples
#' ## create a biexponential excursion-recovery curve with random noise
#' set.seed(13)
#' t <- 0:120
#' x <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B1, lt1, B2, lr), data = data)
#' summary(model)
#'
#' ## convert the log-scale time constants to the natural scale
#' with(as.list(coef(model)), c(tau1 = exp(lt1), tau2 = exp(lt1 + lr)))
#'
#' ## match the analyse_kinetics() fit path: port with box bounds keeping
#' ## tau1 within the data span and tau2 / tau1 >= 2.5
#' span <- diff(range(t))
#' model_port <- nls(
#'     x ~ SSbiexponential(t, A, B1, lt1, B2, lr),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, log(span / 1000), -Inf, log(2.5))
#' )
#' summary(model_port)
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A = 70, B1, lt1, B2, lr), data = data
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
    model = biexponential_ratio,
    initial = init_fixed(biexp_init, c("A", "B1", "lt1", "B2", "lr", "TD")),
    parameters = c("A", "B1", "lt1", "B2", "lr", "TD")
)


#' Fit a biexponential in ratio parameterisation
#'
#' [fit_biexp_ratio()]: Fits [biexponential()] with the time constants
#' re-parameterised as `lt1 = log(tau1)` and `lr = log(tau2 / tau1)`.
#'
#' @details
#' Convergence rests on `algorithm = "port"` with box bounds: the log-ratio
#' scale expresses the `tau2 / tau1 >= tau_ratio` constraint as a simple box
#' bound on `lr`, keeping the fit away from the degenerate `tau2 = tau1`
#' limit where the design matrix is singular. Unbounded fits converge poorly
#' on real NIRS data regardless of scale.
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

    TD_seed <- fix$TD %||% 0
    if (has_TD && is.null(fix$TD)) {
        TD_seed <- biexp_profile_TD(x, t, tau_ratio, span) %||% 0
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

    ## translate natural-scale params and fix to the SS log-ratio scale: a
    ## fixed tau1 becomes a constant lt1; a fixed tau2 becomes an lr
    ## expression in lt1, which build_ss_formula() substitutes verbatim
    params_ss <- unname(c(
        A = "A",
        B1 = "B1",
        tau1 = "lt1",
        B2 = "B2",
        tau2 = "lr",
        TD = "TD"
    )[params])
    fix_ss <- fix[setdiff(names(fix), c("tau1", "tau2"))]
    if (!is.null(fix$tau1)) {
        fix_ss$lt1 <- log(fix$tau1)
    }
    if (!is.null(fix$tau2)) {
        fix_ss$lr <- substitute(
            LT2 - L,
            list(LT2 = log(fix$tau2), L = fix_ss$lt1 %||% quote(lt1))
        )
    }

    ## free parameters in formula order; port matches unnamed bounds
    ## positionally, so start/lower/upper stay aligned with the formula
    free <- setdiff(params_ss, names(fix_ss))
    start <- vapply(bounds[free], `[[`, numeric(1), 1L)
    lower <- vapply(bounds[free], `[[`, numeric(1), 2L)
    upper <- vapply(bounds[free], `[[`, numeric(1), 3L)

    model <- tryCatch(
        nls(
            build_ss_formula(quote(SSbiexponential), params_ss, fix_ss),
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


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' excursion-recovery curve to each `nirs_channel` within a single *"mnirs"*
#' data frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit a 6-parameter
#'   [biexponential()] model (A, B1, tau1, B2, tau2, TD) with a time delay.
#'   If the 6-parameter fit fails, or if `use_TD = FALSE`, attempts to fit a
#'   reduced 5-parameter [biexponential()] model (A, B1, tau1, B2, tau2).
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
#'   `excursion_time`, `excursion_value`. Per-channel metadata are attached as
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
            "excursion_time", "excursion_value"
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

        ## back-convert log-ratio time constants to the natural scale; a
        ## fixed time constant has no log-scale counterpart in the model
        coefs <- c(stats::coef(model), unlist(fix))
        coefs[["tau1"]] <- fix$tau1 %||% exp(coefs[["lt1"]])
        coefs[["tau2"]] <- fix$tau2 %||% (coefs[["tau1"]] * exp(coefs[["lr"]]))
        fitted_vals <- stats::predict(model)
        keep <- keep_rows(params)

        ## model parameters in biexponential() argument order
        pars <- as.list(coefs[c("A", "B1", "tau1", "B2", "tau2")])

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] else NULL

        ## interior excursion where dy/dt = 0 has a closed form, elapsed since
        ## model onset (TD, or 0 for the 5-param model). the root needs only
        ## the ratio (B2/tau2)/(B1/tau1) to be positive, i.e. amplitudes of the
        ## same sign; opposite signs make the curve monotone and the turning
        ## point sits at the onset boundary. amplitudes may be negative, in
        ## which case the excursion is a maximum rather than a minimum
        s <- with(pars, {
            ratio <- (B2 / tau2) / (B1 / tau1)
            root <- if (is.finite(ratio) && ratio > 0) {
                log(ratio) / (1 / tau2 - 1 / tau1)
            } else {
                NA_real_
            }
            if (is.finite(root) && root > 0) root else 0
        })

        ## excursion_time reported elapsed from start_time, mirroring
        ## MRT = TD + tau; adding TD_arg shifts the onset-relative s into the
        ## same frame, so t - TD in biexponential() recovers s exactly
        excursion_time_val <- sum(TD_arg, s)
        excursion_value_val <- do.call(
            biexponential,
            c(list(t = excursion_time_val), pars, list(TD = TD_arg))
        )

        list(
            coefs = data.frame(
                pars,
                TD = TD_arg %||% NA_real_,
                plateau = pars$A - pars$B1 + pars$B2,
                excursion_time = excursion_time_val,
                excursion_value = excursion_value_val
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
