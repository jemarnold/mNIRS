#' Biexponential function
#'
#' Calculate a two-phase curve: a fast monoexponential response toward `B1`
#' and a slow monoexponential response from `B1` toward a stable plateau at
#' `B2`, both clocked from the response onset and summed.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting value of the response
#'   variable (the `t = 0` intercept).
#' @param B1 A numeric parameter for the asymptote of the *fast* component;
#'   the value the fast response alone would approach.
#' @param tau1 A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`. Dominates the initial steep
#'   response.
#' @param B2 A numeric parameter for the asymptote of the *slow* component;
#'   the stable plateau the response recovers toward as `t` approaches
#'   infinity.
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
#'   `A + (B1 - A) * (1 - exp(-t / tau1)) + (B2 - B1) * (1 - exp(-t / tau2))`
#'
#' 6-parameter model (with time delay), where `ts = pmax(t - TD, 0)`:
#'   `A + (B1 - A) * (1 - exp(-ts / tau1)) +
#'   (B2 - B1) * (1 - exp(-ts / tau2))`
#'
#' `A`, `B1`, and `B2` are all values on the response scale. The fast
#' component is a [monoexponential()] response from `A` toward `B1` with
#' amplitude `B1 - A`; the slow component runs concurrently from the same
#' onset with amplitude `B2 - B1`. The curve starts at `A`, approaches `B2`
#' as `t` grows, and is smooth throughout.
#'
#' The expected response is a *fast* excursion toward a minimum or maximum
#' short of `B1`, followed by a *slow* recovery back to a stable plateau at
#' `B2`. The turning point occurs where the two phase rates cancel:
#' `texc = TD + log(r) / (1 / tau1 - 1 / tau2)` with
#' `r = -(B1 - A) * tau2 / ((B2 - B1) * tau1)`, which exists only when the
#' amplitudes oppose in sign and the fast phase dominates at the onset
#' (`r > 1`). If `B1` is between `A` and `B2`, the response is monotonic but
#' still two-phase. If `B1 = B2`, the curve reduces to a [monoexponential()]
#' with single time constant `tau1` and asymptote `B2`.
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
#' x <- biexponential(t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1, B2, tau2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0),
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
biexponential <- function(t, A, B1, tau1, B2, tau2, TD = NULL) {
    ## both phases clocked from the onset; the slow term carries the
    ## response from B1 toward B2
    ts <- if (is.null(TD)) t else pmax(t - TD, 0)
    return(
        A + (B1 - A) * (1 - exp(-ts / tau1)) + (B2 - B1) * (1 - exp(-ts / tau2))
    )
}


## time of the curve turning point, elapsed from the fit origin; NA when
## the amplitudes share a sign or the slow phase dominates from the onset
## (monotonic response), or the time constants coincide
biexp_texc <- function(A, B1, tau1, B2, tau2, TD = NULL) {
    r <- -((B1 - A) * tau2) / ((B2 - B1) * tau1)
    if (!is.finite(r) || r <= 1 || tau1 == tau2) {
        return(NA_real_)
    }
    return(sum(TD, log(r) / (1 / tau1 - 1 / tau2)))
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
    ## in A, B1, and B2 once tau1, tau2, and TD are held, so those are
    ## solved by least squares at every grid point. pairs with tau2 close
    ## to tau1 are dropped as their bases are near-collinear
    tau1_grid <- fixed$tau1 %||%
        exp(seq(log(span / 100), log(span / 2), length.out = 13L))
    tau2_grid <- fixed$tau2 %||%
        exp(seq(log(span / 20), log(span * 10), length.out = 9L))
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, span / 3, length.out = 11L)
    }
    grid <- expand.grid(tau1 = tau1_grid, tau2 = tau2_grid, TD = td_grid)
    grid <- grid[grid$tau2 >= 2 * grid$tau1, , drop = FALSE]

    fits <- vapply(seq_len(nrow(grid)), \(.i) {
        ts <- if (has_TD) pmax(t - grid$TD[.i], 0) else t
        e1 <- exp(-ts / grid$tau1[.i])
        e2 <- exp(-ts / grid$tau2[.i])
        X <- cbind(e1, e2 - e1, 1 - e2)
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
        TD = if (has_TD) grid$TD[best]
    ))
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. Supports both the
#' 5-parameter form (A, B1, tau1, B2, tau2) and the 6-parameter form adding
#' a time delay TD; arity is inferred from the formula passed to
#' [stats::nls()].
#'
#' @usage
#' SSbiexponential(t, A, B1, tau1, B2, tau2, TD)
#'
#' @inheritParams biexponential
#'
#' @details
#' 5-parameter model:
#' `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2)`
#'
#' 6-parameter model:
#' `x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, TD)`
#'
#' The two phases are weakly identified when `tau1` and `tau2` are close, so
#' `algorithm = "port"` with the time constants bounded non-negative and
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#'
#' The 5-parameter form is recommended for small samples or when no obvious
#'   time delay is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `TD` incurs no degrees-of-freedom penalty.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSbiexponential(t, A, B1, tau1 = 5, B2, tau2)`
#'   holds the fast time constant at `5`. Fixed parameters are excluded from
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
#' x <- biexponential(t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1, B2, tau2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' ## fix the fast time constant
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A, B1, tau1 = 5, B2, tau2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model_fixed)
#'
#' @export
SSbiexponential <- selfStart(
    model = biexponential,
    initial = init_fixed(
        biexp_init,
        c("A", "B1", "tau1", "B2", "tau2", "TD")
    ),
    parameters = c("A", "B1", "tau1", "B2", "tau2", "TD")
)


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' excursion-recovery curve to each `nirs_channel` within a single *"mnirs"*
#' data frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit a 6-parameter
#'   [SSbiexponential()] model (A, B1, tau1, B2, tau2, TD) with a time
#'   delay. If the 6-parameter fit fails, or if `use_TD = FALSE`, attempts
#'   to fit a reduced 5-parameter [SSbiexponential()] model without `TD`.
#' @param fix An *optional* named list of model parameters (`A`, `B1`,
#'   `tau1`, `B2`, `tau2`, `TD`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation and
#'   reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. `TD` is fixable for channels where
#'   `use_TD = TRUE`; a fixed `TD` disables the 5-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B1`, `tau1`, `B2`, `tau2`, `TD`, `texc`,
#'   `texc_fitted`. Per-channel metadata are attached as attributes:
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
        arg_list = mget(
            c("use_TD", "fix", "start_time", "direction", "end_window")
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        ## TD is only fixable where that channel fits the 6-parameter model
        fix_params = \(.a) {
            c("A", "B1", "tau1", "B2", "tau2", if (.a$use_TD) "TD")
        },
        verbose = verbose,
        env = env
    )

    time_channel <- setup$time_channel
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- c("A", "B1", "tau1", "B2", "tau2", "TD", "texc", "texc_fitted")

    ## method-specific fit: self-starting biexponential via nls; a failed
    ## 6-param fit falls back to the 5-param model
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## the phases are weakly identified when the time constants are
        ## close, so port often stops short of its certificate on usable
        ## coefficients, which are kept with a warning
        fitter <- \(.data, .params, on_error) {
            free <- setdiff(.params, names(.a$fix))
            span <- diff(range(.data[[2L]]))
            formula <- build_ss_formula(
                quote(SSbiexponential),
                .params,
                .a$fix,
                names(.data)[[1L]],
                names(.data)[[2L]]
            )

            ## taus floored as degeneracy markers; tau2 capped so a
            ## runaway slow component cannot diverge with B2
            lower <- c(tau1 = span * 1e-6, tau2 = span * 1e-6, TD = 0)[free]
            names(lower) <- free
            lower[is.na(lower)] <- -Inf
            upper <- c(tau2 = 10 * span)[free]
            names(upper) <- free
            upper[is.na(upper)] <- Inf

            port_fit <- \(start) {
                tryCatch(
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
            }
            start <- tryCatch(
                stats::getInitial(formula, .data),
                error = on_error
            )
            if (is.null(start)) {
                return(NULL)
            }
            model <- port_fit(start)

            ## the phases are exchangeable, so the optimum may land with
            ## the fast label on the slow term; refit from the swapped
            ## start (same RSS optimum) so tau1 <= tau2 is reported
            cf <- if (is.null(model)) NULL else coef(model)
            swap <- all(c("tau1", "tau2", "B1") %in% free) &&
                !is.null(cf) &&
                isTRUE(cf[["tau1"]] > cf[["tau2"]])
            if (swap) {
                A_val <- .a$fix$A %||% cf[["A"]]
                B2_val <- .a$fix$B2 %||% cf[["B2"]]
                cf[c("tau1", "tau2", "B1")] <- c(
                    cf[["tau2"]], cf[["tau1"]], A_val + B2_val - cf[["B1"]]
                )
                model <- port_fit(cf) %||% model
            }
            accept_port_fit(model, on_error)
        }

        fit <- fit_td_fallback(
            x_fit,
            t_fit,
            params = c("A", "B1", "tau1", "B2", "tau2", if (.a$use_TD) "TD"),
            .a,
            fitter,
            fn = quote(SSbiexponential),
            .nirs = .nirs,
            time_channel = time_channel,
            interval_name = interval_name,
            env = env
        )
        if (is.null(fit$model)) {
            return(build_na_results(na_cols))
        }
        params <- fit$params
        coefs <- full_coefs(fit$model, params, .a$fix)
        span <- diff(range(fit$data[[time_channel]]))

        ## enforce direction: bounded refit on the fast-phase amplitude
        ## D = B1 - A when inverted. the fit bounds are carried over;
        ## TD >= 0 is structural, so only D and the tau floors mark a
        ## degenerate fit
        free <- setdiff(params, names(.a$fix))
        lower <- c(tau1 = span * 1e-6, tau2 = span * 1e-6, TD = 0)
        upper <- c(tau2 = 10 * span)
        enforced <- enforce_direction(
            fit$model,
            coefs,
            fit$data,
            direction = .a$direction,
            amp_fn = quote(biexponential),
            lower = lower[intersect(names(lower), free)],
            upper = upper[intersect(names(upper), free)],
            floor_params = c("D", "tau1", "tau2"),
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

        ## excursion time (texc) is the fitted turning point, reported
        ## elapsed from start_time, mirroring MRT = TD + tau; NA when the
        ## fitted response is monotonic
        texc_val <- biexp_texc(
            A = coefs[["A"]],
            B1 = coefs[["B1"]],
            tau1 = coefs[["tau1"]],
            B2 = coefs[["B2"]],
            tau2 = coefs[["tau2"]],
            TD = TD_arg
        )
        texc_fitted_val <- if (is.na(texc_val)) {
            NA_real_
        } else {
            biexponential(
                t = texc_val,
                A = coefs[["A"]],
                B1 = coefs[["B1"]],
                tau1 = coefs[["tau1"]],
                B2 = coefs[["B2"]],
                tau2 = coefs[["tau2"]],
                TD = TD_arg
            )
        }

        build_fit_results(
            data.frame(
                A = coefs[["A"]],
                B1 = coefs[["B1"]],
                tau1 = coefs[["tau1"]],
                B2 = coefs[["B2"]],
                tau2 = coefs[["tau2"]],
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
        setup$per_channel,
        biexp_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
