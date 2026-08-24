#' Biexponential function
#'
#' Calculate a two-phase curve: a fast monoexponential response toward `B1`,
#' followed by a slow component that begins at the excursion point
#' `texc = TD + tau_mult * tau1` and carries the response toward a stable
#' plateau at `B2`.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting value of the response
#'   variable (the `t = 0` intercept).
#' @param B1 A numeric parameter for the asymptote of the initial *fast*
#'   component; the value the initial response approaches at the excursion
#'   point.
#' @param tau1 A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`. Dominates the initial steep
#'   response.
#' @param B2 A numeric parameter for the asymptote of the *slow* component;
#'   the stable plateau the response recovers toward as `t` approaches
#'   infinity.
#' @param tau2 A numeric parameter for the *slow* time constant (\eqn{\tau_2}),
#'   in units of the predictor variable `t`. Typically `tau2 >> tau1`.
#' @param tau_mult A numeric multiple of `tau1` after `TD` at which the slow
#'   component begins: excursion point; `texc = TD + tau_mult * tau1`
#'   (`TD = 0` when absent).
#' @param TD A numeric parameter for the time delay before the onset of the
#'   response, in units of the predictor variable `t`. If `NULL` (*default*),
#'   a 6-parameter model without time delay is used.
#'
#' @details
#' This model family is fit by [analyse_kinetics()] when
#' `method = "biexponential"`, and by [stats::nls()] via the self-starting
#' wrapper [SSbiexponential()].
#'
#' ## Model equations
#'
#' 6-parameter model:
#'   `A + (B1 - A) * (1 - exp(-t / tau1)) +
#'   (B2 - B1) * (1 - exp(-pmax(t - tau_mult * tau1, 0) / tau2))`
#'
#' 7-parameter model (with time delay), where `ts = pmax(t - TD, 0)`:
#'   `A + (B1 - A) * (1 - exp(-ts / tau1)) +
#'   (B2 - B1) * (1 - exp(-pmax(ts - tau_mult * tau1, 0) / tau2))`
#'
#' `A`, `B1`, and `B2` are all values on the response scale. The fast
#' component is a [monoexponential()] response from `A` toward `B1`. The slow
#' component is exactly zero before the excursion point
#' `texc = TD + tau_mult * tau1`, then follows its own exponential clock with
#' amplitude `B2 - B1` toward the plateau `B2`. The curve is continuous but
#' not differentiable at the excursion point.
#'
#' The expected response is a *fast* overshoot toward a minimum or maximum
#' near `B1`, followed by a *slow* recovery back to a stable plateau at `B2`.
#' If `B1` is between `A` and `B2`, the response curve is monotonic but still
#' two-phase. If `B1 = B2`, the curve reduces to a [monoexponential()] with
#' single time constant `tau1` and asymptote `B2`.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSbiexponential()], [monoexponential()],
#'   [exponential_drift()]
#'
#' @examples
#' ## create a biexponential excursion-recovery curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(
#'     t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, tau_mult = 2
#' ) + rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0, 1),
#'     control = nls.control(warnOnly = TRUE)
#' )
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
biexponential <- function(t, A, B1, tau1, B2, tau2, tau_mult, TD = NULL) {
    ## fast monoexponential toward B1, plus a slow component whose own
    ## clock starts at the excursion point texc = TD + tau_mult * tau1
    ts <- if (is.null(TD)) t else pmax(t - TD, 0)
    return(
        (A + (B1 - A) * (1 - exp(-ts / tau1))) +
            (B2 - B1) * (1 - exp(-pmax(ts - tau_mult * tau1, 0) / tau2))
    )
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
#'   user-fixed parameter values from [init_fixed()] used to narrow the
#'   grids.
#'
#' @returns [biexp_init()]: Initial starting estimates for parameters in the
#'   model called by [SSbiexponential()].
#'
#' @keywords internal
biexp_init <- function(mCall, data, LHS, ...) {
    ## user-fixed values narrow the grids; the amplitudes are always
    ## solved free, as this is only a seed
    fixed <- list(...)$fixed %||% list()

    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    has_TD <- "TD" %in% names(mCall)
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        span <- 1
    }

    ## profile the time constants (and TD) on a coarse grid and keep the
    ## RSS-minimising start (cf. `expdrift_init()`). the model is linear
    ## in A, B1, and B2 once tau1, tau2, tau_mult, and TD are held, so
    ## those are solved by least squares at every grid point. pairs with
    ## tau2 < tau1 or a slow onset beyond the record are dropped; a grid
    ## point whose hinge has no support solves to NA and is skipped
    tau1_grid <- fixed$tau1 %||%
        exp(seq(log(span / 100), log(span / 2), length.out = 13L))
    tau2_grid <- fixed$tau2 %||%
        exp(seq(log(span / 20), log(span * 10), length.out = 9L))
    mult_grid <- fixed$tau_mult %||% c(1, 2, 3)
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, span / 3, length.out = 11L)
    }
    grid <- expand.grid(
        tau1 = tau1_grid,
        tau2 = tau2_grid,
        tau_mult = mult_grid,
        TD = td_grid
    )
    grid <- grid[
        grid$tau2 >= grid$tau1 &
            grid$TD + grid$tau_mult * grid$tau1 < span,
        ,
        drop = FALSE
    ]

    fits <- vapply(seq_len(nrow(grid)), \(.i) {
        ts <- if (has_TD) pmax(t - grid$TD[.i], 0) else t
        e1 <- exp(-ts / grid$tau1[.i])
        h <- 1 -
            exp(-pmax(ts - grid$tau_mult[.i] * grid$tau1[.i], 0) /
                grid$tau2[.i])
        X <- cbind(e1, 1 - e1 - h, h)
        cf <- qr.coef(qr(X), x)
        c(cf, sum((x - X %*% cf)^2))
    }, numeric(4L))
    best <- which.min(fits[4L, ])
    if (length(best) == 0L) {
        stop("No starting estimates could be resolved from the response.")
    }

    return(c(
        A = fits[[1L, best]],
        B1 = fits[[2L, best]],
        tau1 = grid$tau1[best],
        B2 = fits[[3L, best]],
        tau2 = grid$tau2[best],
        tau_mult = grid$tau_mult[best],
        TD = if (has_TD) grid$TD[best]
    ))
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. Supports both the
#' 6-parameter form (A, B1, tau1, B2, tau2, tau_mult) and the 7-parameter
#' form adding a time delay TD; arity is inferred from the formula passed
#' to [stats::nls()].
#'
#' @usage
#' SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult, TD)
#'
#' @inheritParams biexponential
#'
#' @details
#' 6-parameter model:
#' `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult)`
#'
#' 7-parameter model:
#' `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult, TD)`
#'
#' The hinge at `texc = TD + tau_mult * tau1` is not differentiable, so
#' `algorithm = "port"` with the time constants bounded non-negative,
#' `tau_mult >= 1`, and `control = nls.control(warnOnly = TRUE)` is
#' recommended.
#'
#' The 6-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult = 2)`
#'   holds the slow onset at `2 * tau1`. Fixed parameters are excluded from
#'   estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [biexponential()], [stats::nls()], [stats::selfStart()],
#'   [SSmonoexponential()], [SSexponential_drift()]
#'
#' @examples
#' ## create a biexponential excursion-recovery curve with random noise
#' set.seed(13)
#' t <- 0:120
#' x <- biexponential(
#'     t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, tau_mult = 2
#' ) + rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0, 1),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' ## fix the slow onset at 2 tau1
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, tau_mult = 2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model_fixed)
#'
#' @export
SSbiexponential <- selfStart(
    model = biexponential,
    initial = init_fixed(
        biexp_init,
        c("A", "B1", "tau1", "B2", "tau2", "tau_mult", "TD")
    ),
    parameters = c("A", "B1", "tau1", "B2", "tau2", "tau_mult", "TD")
)


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' excursion-recovery curve to each `nirs_channel` within a single *"mnirs"*
#' data frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit a 7-parameter
#'   [SSbiexponential()] model (A, B1, tau1, B2, tau2, tau_mult, TD) with a
#'   time delay. If the 7-parameter fit fails, or if `use_TD = FALSE`,
#'   attempts to fit a reduced 6-parameter [SSbiexponential()] model without
#'   `TD`.
#' @param tau_mult An *optional* numeric multiple of `tau1` after `TD` at
#'   which the slow component onset is held
#'   (`texc = TD + tau_mult * tau1`). If `NULL` (*default*), `tau_mult` is a
#'   free parameter estimated from the data, bounded below by `1` so the
#'   slow phase begins no earlier than one fast time constant. A supplied
#'   value is held constant during fitting and silently takes precedence
#'   over a `fix = list(tau_mult = )` entry. Applied to every channel, or
#'   per-channel as a list keyed by channel name, e.g.
#'   `tau_mult = list(smo2 = 2)`.
#' @param fix An *optional* named list of model parameters (`A`, `B1`,
#'   `tau1`, `B2`, `tau2`, `tau_mult`, `TD`) to hold constant during
#'   fitting, e.g. `fix = list(A = 0)`. Fixed parameters are excluded from
#'   estimation and reported at their fixed values. Applied to every
#'   channel, or per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. `TD` is fixable for channels where
#'   `use_TD = TRUE`; a fixed `TD` disables the 6-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B1`, `tau1`, `B2`, `tau2`, `tau_mult`, `TD`,
#'   `texc`, `texc_fitted`. Per-channel metadata are attached as
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
    use_TD = TRUE,
    tau_mult = NULL,
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
            c("use_TD", "tau_mult", "fix", "start_time",
            "direction", "end_window")
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        ## TD is only fixable where that channel fits the 7-parameter model
        fix_params = \(.a) {
            c("A", "B1", "tau1", "B2", "tau2", "tau_mult", if (.a$use_TD) "TD")
        },
        verbose = verbose,
        env = env
    )
    ## a supplied tau_mult fixes the slow onset multiple, taking precedence
    ## over any `fix` entry; NULL leaves it free with a lower bound of 1.
    ## a zero multiple would start the slow phase at the response onset,
    ## where the two components are indistinguishable
    per_channel <- lapply(setup$per_channel, \(.a) {
        if (!is.null(.a$tau_mult)) {
            validate_numeric(
                .a$tau_mult, 1, c(0, Inf),
                inclusive = "right", msg1 = "one-element positive", env = env
            )
            .a$fix$tau_mult <- .a$tau_mult
        }
        .a
    })

    ## NA scaffold (method columns only) for convergence failure
    # fmt: skip
    na_cols <- c(
        "A", "B1", "tau1", "B2", "tau2", "tau_mult", "TD",
        "texc", "texc_fitted"
    )

    ## method-specific fit: self-starting biexponential via nls; a failed
    ## 7-param fit falls back to the 6-param model
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## the hinge is non-smooth, so port often stops short of its
        ## certificate on usable coefficients, which are kept with a
        ## warning. reads `.a$fix` at call time so a per-channel fixed
        ## tau_mult is honoured
        fitter <- \(.data, .params, on_error) {
            free <- setdiff(.params, names(.a$fix))
            span <- diff(range(.data$.t))
            formula <- build_ss_formula(
                quote(SSbiexponential),
                .params,
                .a$fix
            )

            ## explicit getInitial so the tau_mult cap can be derived
            ## from the tau1 seed before the bounds are fixed
            start <- tryCatch(
                stats::getInitial(formula, .data),
                error = on_error
            )
            if (is.null(start)) {
                return(NULL)
            }
            tau1_seed <- .a$fix$tau1 %||% start[["tau1"]]

            ## taus floored as degeneracy markers; tau2 capped so a
            ## runaway slow component cannot diverge with B2; tau_mult
            ## bounded so the slow phase starts after the fast and its
            ## onset stays inside the record
            lower <- c(
                tau1 = span * 1e-6,
                tau2 = span * 1e-6,
                tau_mult = 1,
                TD = 0
            )[free]
            names(lower) <- free
            lower[is.na(lower)] <- -Inf
            upper <- c(
                tau2 = 10 * span,
                tau_mult = max(span / tau1_seed, 1)
            )[free]
            names(upper) <- free
            upper[is.na(upper)] <- Inf

            ## seed clamped inside the bounds; a fixed tau1 can nudge the
            ## grid seed past the tau_mult cap
            start <- pmin(
                pmax(start, lower[names(start)]),
                upper[names(start)]
            )

            model <- tryCatch(
                suppressWarnings(nls(
                    formula,
                    .data,
                    start = start,
                    algorithm = "port",
                    lower = lower,
                    upper = upper,
                    control = stats::nls.control(
                        maxiter = 500L,
                        warnOnly = TRUE
                    )
                )),
                error = on_error
            )
            accept_port_fit(model, on_error)
        }

        fit <- fit_td_fallback(
            x_fit,
            t_fit,
            # fmt: skip
            params = c(
                "A", "B1", "tau1", "B2", "tau2", "tau_mult",
                if (.a$use_TD) "TD"
            ),
            .a,
            fitter,
            fn = quote(SSbiexponential),
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(fit$model)) {
            return(build_na_results(na_cols))
        }
        params <- fit$params
        coefs <- full_coefs(fit$model, params, .a$fix)
        span <- diff(range(fit$data$.t))

        ## enforce direction: bounded refit on D = B2 - A when inverted.
        ## the fit bounds are carried over; tau_mult >= 1 and TD >= 0 are
        ## structural, so only D and the tau floors mark a degenerate fit
        free_extra <- setdiff(params, c("A", "B2", names(.a$fix)))
        lower_all <- c(
            tau1 = span * 1e-6,
            tau2 = span * 1e-6,
            tau_mult = 1,
            TD = 0
        )
        upper_all <- c(
            tau2 = 10 * span,
            tau_mult = max(span / coefs[["tau1"]], 1)
        )
        enforced <- enforce_direction(
            fit$model,
            coefs,
            fit$data,
            direction = .a$direction,
            amp_fn = quote(biexponential),
            extra = coefs[free_extra],
            B_name = "B2",
            extra_lower = lower_all[intersect(names(lower_all), free_extra)],
            extra_upper = upper_all[intersect(names(upper_all), free_extra)],
            floor_params = c("D", "tau1", "tau2"),
            fn = quote(SSbiexponential),
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

        ## excursion time (texc) is reported elapsed from start_time,
        ## mirroring MRT = TD + tau
        texc_val <- sum(TD_arg, coefs[["tau_mult"]] * coefs[["tau1"]])
        texc_fitted_val <- biexponential(
            t = texc_val,
            A = coefs[["A"]],
            B1 = coefs[["B1"]],
            tau1 = coefs[["tau1"]],
            B2 = coefs[["B2"]],
            tau2 = coefs[["tau2"]],
            tau_mult = coefs[["tau_mult"]],
            TD = TD_arg
        )

        build_fit_results(
            data.frame(
                A = coefs[["A"]],
                B1 = coefs[["B1"]],
                tau1 = coefs[["tau1"]],
                B2 = coefs[["B2"]],
                tau2 = coefs[["tau2"]],
                tau_mult = coefs[["tau_mult"]],
                TD = TD_arg %||% NA_real_,
                texc = texc_val,
                texc_fitted = texc_fitted_val
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
        per_channel,
        biexp_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
