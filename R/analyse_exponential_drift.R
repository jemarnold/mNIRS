#' Exponential-drift function
#'
#' Calculate a two-phase curve: a primary monoexponential response with a
#' secondary linear drift beginning near the asymptote.
#'
#' @param slope A numeric parameter for the linear drift rate `dx/dt`
#'   of the secondary phase, in response units per unit of the predictor
#'   variable `t`.
#' @param texc A numeric parameter for the time at which the linear
#'   drift begins, in units of the predictor variable `t` (the same time
#'   frame as `TD`, not relative to it).
#' @inheritParams monoexponential
#'
#' @details
#' This model family is fit by [analyse_kinetics()] when
#' `method = "exponential_drift"`, and by [stats::nls()] via the
#' self-starting wrapper [SSexponential_drift()].
#'
#' ## Model equations
#'
#' 5-parameter model:
#' `A + (B - A) * (1 - exp(-t / tau)) + slope * pmax(t - texc, 0)`
#'
#' 6-parameter model:
#' `A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau)) +
#' slope * pmax(t - texc, 0)`
#'
#' The primary phase is a pure [monoexponential()] response toward the
#' asymptote `B`. The secondary phase is a linear drift of rate `slope`
#' that is exactly zero before `texc`, so unlike [biexponential()] the
#' two phases do not overlap from `t = 0`; the drift only perturbs the curve
#' after the primary response is (near) complete.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSexponential_drift()],
#'   [monoexponential()], [biexponential()]
#'
#' @examples
#' ## create an exponential curve with late linear drift and random noise
#' set.seed(13)
#' t <- 1:180
#' x <- exponential_drift(
#'     t, A = 10, B = 100, tau = 12, slope = -0.1, texc = 60, TD = 15
#' ) + rnorm(length(t), 0, 3)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSexponential_drift(t, A, B, tau, slope, texc, TD),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0, 0)
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
exponential_drift <- function(t, A, B, tau, slope, texc, TD = NULL) {
    ## primary monoexponential phase + hinge-linear secondary drift
    return(monoexponential(t, A, B, tau, TD) + slope * pmax(t - texc, 0))
}


#' Initiate self-starting exponential-drift model
#'
#' [expdrift_init()]: Returns initial values for the parameters in a
#' `selfStart` model.
#'
#' @inheritParams monoexp_init
#'
#' @returns [expdrift_init()]: Initial starting estimates for parameters in
#'   the model called by [SSexponential_drift()].
#'
#' @keywords internal
expdrift_init <- function(mCall, data, LHS, ...) {
    ## user-fixed parameter values constrain the free estimates
    fixed <- list(...)$fixed %||% list()

    ## seed the primary phase from the monoexponential grid search,
    ## ignoring the drift; the joint port fit refines all parameters
    mono <- monoexp_init(
        mCall,
        data,
        LHS,
        fixed = fixed[names(fixed) %in% c("A", "B", "tau", "TD")]
    )
    TD_seed <- if ("TD" %in% names(mono)) mono[["TD"]]

    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        span <- 1
    }

    ## drift onset texc near the primary asymptote (~95% at 3 tau), held
    ## inside the record so the hinge basis has support
    texc <- fixed$texc %||%
        max(min(sum(TD_seed, 3 * mono[["tau"]]), min(t) + 0.75 * span), 0)

    ## least-squares slope of the primary-phase residuals on the hinge
    ## basis; a degenerate (all pre-texc) basis seeds a flat drift
    e <- x -
        monoexponential(t, mono[["A"]], mono[["B"]], mono[["tau"]], TD_seed)
    h <- pmax(t - texc, 0)
    slope <- fixed$slope %||% (if (sum(h^2) > 0) sum(h * e) / sum(h^2) else 0)

    return(c(
        A = mono[["A"]],
        B = mono[["B"]],
        tau = mono[["tau"]],
        slope = slope,
        texc = texc,
        TD = TD_seed
    ))
}


#' Self-starting exponential-drift model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [exponential_drift()], for use with [stats::nls()]. Supports both the
#' 5-parameter form (A, B, tau, slope, texc) and the
#' 6-parameter form adding a time delay TD; arity is inferred from the
#' formula passed to [stats::nls()].
#'
#' @usage
#' SSexponential_drift(t, A, B, tau, slope, texc, TD)
#'
#' @inheritParams exponential_drift
#'
#' @details
#' 5-parameter model:
#' `x ~ SSexponential_drift(t, A, B, tau, slope, texc)`
#'
#' 6-parameter model:
#' `x ~ SSexponential_drift(t, A, B, tau, slope, texc, TD)`
#'
#' The hinge at `texc` is not differentiable, so
#' `algorithm = "port"` with a lower bound holding `texc` (and `tau`)
#' non-negative is recommended over the default Gauss-Newton algorithm.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSexponential_drift(t, A = 0, B, tau, slope, texc)`
#'   fixes the baseline at `A = 0`. Fixed parameters are excluded from
#'   estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [exponential_drift()], [stats::nls()], [stats::selfStart()],
#'   [SSmonoexponential()]
#'
#' @examples
#' ## create an exponential curve with late linear drift and random noise
#' set.seed(13)
#' t <- 1:180
#' x <- exponential_drift(
#'     t, A = 10, B = 100, tau = 12, slope = -0.1, texc = 60, TD = 15
#' ) + rnorm(length(t), 0, 3)
#' data <- data.frame(t, x)
#'
#' ## 6-parameter fit
#' model6 <- nls(
#'     x ~ SSexponential_drift(t, A, B, tau, slope, texc, TD),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0, 0)
#' )
#' summary(model6)
#'
#' ## 5-parameter fit on the same data
#' model5 <- nls(
#'     x ~ SSexponential_drift(t, A, B, tau, slope, texc),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0)
#' )
#' summary(model5)
#'
#' @export
SSexponential_drift <- selfStart(
    model = exponential_drift,
    initial = init_fixed(
        expdrift_init,
        c("A", "B", "tau", "slope", "texc", "TD")
    ),
    parameters = c("A", "B", "tau", "slope", "texc", "TD")
)


#' Tie the drift onset texc to the primary time constant
#'
#' Builds the `fix` value substituting `texc = TD + k * tau` into the
#' model formula, referencing the free parameters as symbols and any fixed
#' `TD`/`tau` as constants. Reduces to a numeric constant when both are
#' fixed.
#'
#' @param fix Named list of fixed parameter values on the natural scale.
#' @param params Character vector of parameter names in model order.
#' @param k A numeric multiple of `tau` after `TD` at which the drift begins.
#'
#' @returns A numeric constant or a language expression suitable for
#'   [build_ss_formula()].
#'
#' @keywords internal
fix_tied_texc <- function(fix, params, k) {
    TD_part <- if ("TD" %in% params) fix$TD %||% quote(TD) else 0
    tau_part <- fix$tau %||% quote(tau)
    expr <- call("+", TD_part, call("*", k, tau_part))
    if (is.numeric(TD_part) && is.numeric(tau_part)) {
        return(eval(expr))
    }
    return(expr)
}


#' Fit an exponential-drift model with box bounds
#'
#' [fit_expdrift()]: Fits [exponential_drift()] with [stats::nls()] using
#' `algorithm = "port"`. The hinge at `texc` is non-smooth, so the drift
#' onset is box-bounded inside the record and `tau` is held positive;
#' asymptotes and the drift slope are unbounded.
#'
#' A non-converged fit is accepted with a warning when its port stop code is
#' recognised and its coefficients are finite; otherwise it is rejected and
#' `NULL` is returned.
#'
#' @param x,t Numeric vectors of the response and predictor variables.
#' @param params Character vector of parameter names in model order.
#' @param fix Named list of fixed parameter values or expressions.
#' @param on_error A function called with the [stats::nls()] error condition,
#'   or with a warning condition when a non-converged fit is accepted.
#'
#' @returns An [nls][stats::nls] model, or `NULL`.
#'
#' @keywords internal
fit_expdrift <- function(x, t, params, fix, on_error) {
    ## every failure route reports its condition and yields NULL
    fail <- \(.e) {
        on_error(.e)
        NULL
    }

    ## reject an under-determined fit before it reaches nls
    n_free <- length(setdiff(params, names(fix)))
    if (length(x) <= n_free) {
        return(fail(simpleError(sprintf(
            "%d observation%s for %d free parameters.",
            length(x),
            if (length(x) == 1L) "" else "s",
            n_free
        ))))
    }

    span <- diff(range(t))

    ## box bounds keep the hinge inside the record and the fit away from
    ## degenerate limits; asymptotes and drift slope are unbounded.
    ## lower/upper align positionally with the free parameters in model
    ## order, as expected by port
    # fmt: skip
    bounds <- list(
        A     = c(-Inf, Inf),
        B     = c(-Inf, Inf),
        tau   = c(span * 1e-6, Inf),
        slope = c(-Inf, Inf),
        texc  = c(0, max(t)),
        TD    = c(0, span)
    )
    free <- setdiff(params, names(fix))
    lower <- vapply(bounds[free], `[[`, numeric(1), 1L)
    upper <- vapply(bounds[free], `[[`, numeric(1), 2L)

    model <- tryCatch(
        suppressWarnings(nls(
            build_ss_formula(quote(SSexponential_drift), params, fix),
            data.frame(.x = x, .t = t),
            algorithm = "port",
            lower = lower,
            upper = upper,
            control = stats::nls.control(maxiter = 500L, warnOnly = TRUE)
        )),
        error = fail
    )

    ## a non-converged fit with finite coefficients is kept with a warning;
    ## the port stop code is reported in prose either way
    if (!is.null(model) && !model$convInfo$isConv) {
        port_msg <- c(
            "7" = "Singular convergence: parameters not individually identifiable.",
            "8" = "False convergence: gradient certificate failed near a non-smooth point.",
            "9" = "Function evaluation limit reached without convergence.",
            "10" = "Iteration limit reached without convergence."
        )
        code <- as.character(model$convInfo$stopCode)
        msg <- if (code %in% names(port_msg)) {
            port_msg[[code]]
        } else {
            model$convInfo$stopMessage
        }
        if (code %in% names(port_msg) && all(is.finite(stats::coef(model)))) {
            on_error(simpleWarning(msg))
        } else {
            model <- fail(simpleError(msg))
        }
    }

    return(model)
}


#' Analyse exponential-drift kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "exponential_drift")`. Fits a two-phase
#' monoexponential + linear-drift curve to each `nirs_channel` within a
#' single *"mnirs"* data frame. See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param use_TD Logical; default is `TRUE` to attempt to fit a 6-parameter
#'   [SSexponential_drift()] model (A, B, tau, slope, texc, TD)
#'   with a time delay. If the 6-parameter fit fails, or if `use_TD = FALSE`,
#'   attempts to fit a reduced 5-parameter model without `TD`.
#' @param drift_k A numeric multiple of `tau` after `TD` at which the drift
#'   onset is tied (`texc = TD + drift_k * tau`; default is `3`,
#'   ~95% of the primary amplitude) when the freely fitted texc fails or
#'   converges against its bounds. Always applied globally.
#' @param fix An *optional* named list of model parameters (`A`, `B`, `tau`,
#'   `slope`, `texc`, `TD`) to hold constant during fitting,
#'   e.g. `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. `TD` is fixable for channels where
#'   `use_TD = TRUE`; a fixed `TD` disables the 5-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `tau`, `k`, `TD`, `MRT`, `HRT`,
#'   `MRT_fitted`, `HRT_fitted`, `slope`, `texc`. Per-channel
#'   metadata are attached as attributes:
#'   - `"model"`: an [nls][stats::nls] model object, or `NULL` for channels
#'     where fitting failed.
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [exponential_drift()],
#'   [SSexponential_drift()]
#'
#' @keywords internal
analyse_exponential_drift <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    use_TD = TRUE,
    drift_k = 3,
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
            c("A", "B", "tau", "slope", "texc", if (.a$use_TD) "TD")
        },
        verbose = verbose,
        env = env
    )
    nirs_channels <- setup$nirs_channels
    time_channel <- setup$time_channel
    per_channel <- setup$per_channel

    ## a zero multiple would start the drift at the response onset, where
    ## the drift line absorbs the primary response
    validate_numeric(
        drift_k, 1, c(0, Inf),
        inclusive = "right", msg1 = "one-element positive", env = env
    )

    ## NA scaffold (method columns only) for convergence failure
    # fmt: skip
    na_coefs <- as.data.frame(setNames(
        rep(list(NA_real_), 11L),
        c("A", "B", "tau", "k", "TD", "MRT", "HRT",
        "MRT_fitted", "HRT_fitted", "slope", "texc")
    ))

    ## method-specific fit: self-starting exponential-drift via nls
    expdrift_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        params_all <- c("A", "B", "tau", "slope", "texc", if (.a$use_TD) "TD")

        ## the TD model is flat at `A` before `TD`, so the pre-texc
        ## baseline anchors `A`. the reduced model has no such region and
        ## diverges at t < 0, so it is fit from `start_time` onward
        keep_rows <- function(.params) {
            if ("TD" %in% .params) rep(TRUE, length(t_fit)) else t_fit >= 0
        }

        retry_TD <- .a$use_TD && !"TD" %in% names(.a$fix)
        texc_free <- !"texc" %in% names(.a$fix)

        attempt <- \(.params, .fix, .retry) {
            keep <- keep_rows(.params)
            fit_expdrift(x_fit[keep], t_fit[keep], .params, .fix, \(e) {
                warn_fit_failed(
                    quote(SSexponential_drift),
                    e,
                    .nirs,
                    interval_name,
                    length(.params),
                    .retry,
                    env
                )
            })
        }

        ## mono-style TD fallback within a pathway; `fix_for` rebuilds the
        ## fix list per arity so a tied texc only references free symbols
        fit_pathway <- function(fix_for) {
            params <- params_all
            model <- attempt(params, fix_for(params), retry_TD)
            if (is.null(model) && retry_TD) {
                params <- setdiff(params, "TD")
                model <- attempt(params, fix_for(params), FALSE)
            }
            list(model = model, params = params, fix = fix_for(params))
        }

        ## a free texc converged against its bounds marks a degenerate
        ## fit: at the record end there is no drift region; at zero the
        ## drift line absorbs the primary response
        texc_pinned <- function(model) {
            cf <- stats::coef(model)
            if (!"texc" %in% names(cf)) {
                return(FALSE)
            }
            short <- 0.05 * diff(range(t_fit))
            cf[["texc"]] >= max(t_fit) - short || cf[["texc"]] <= short
        }

        ## primary pathway: freely fitted drift onset at texc
        fit <- fit_pathway(\(.p) .a$fix)

        ## fallback pathway: texc tied to `TD + drift_k * tau`
        if (texc_free && (is.null(fit$model) || texc_pinned(fit$model))) {
            if (!is.null(fit$model)) {
                warn_fit_failed(
                    quote(SSexponential_drift),
                    simpleError(sprintf(
                        "`texc` converged against its bounds; refit with `texc = TD + %s * tau`.",
                        format(drift_k)
                    )),
                    .nirs,
                    interval_name,
                    length(fit$params),
                    FALSE,
                    env
                )
            }
            fit <- fit_pathway(\(.p) {
                f <- .a$fix
                f$texc <- fix_tied_texc(.a$fix, .p, drift_k)
                f
            })
        }
        model <- fit$model
        params <- fit$params
        fix_used <- fit$fix

        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        ## a tied texc is a language fix with no constant to merge
        fix_num <- fix_used[vapply(fix_used, is.numeric, logical(1))]
        coefs <- full_coefs(model, params, fix_num)

        ## enforce direction: bounded refit on D = B - A when inverted.
        ## tau/texc/TD bounds are structural, so only D and a tau pinned
        ## at its data-scaled floor mark a degenerate fit
        keep <- keep_rows(params)
        span <- diff(range(t_fit[keep]))
        fit_data <- data.frame(.x = x_fit[keep], .t = t_fit[keep])
        free_extra <- setdiff(names(stats::coef(model)), c("A", "B"))
        # fmt: skip
        extra_bounds <- list(
            tau  = c(span * 1e-6, Inf),
            texc = c(0, max(t_fit[keep])),
            TD   = c(0, span)
        )
        extra_bounds <- extra_bounds[intersect(names(extra_bounds), free_extra)]
        enforced <- enforce_direction(
            model,
            coefs,
            fit_data,
            direction = .a$direction,
            amp_fn = quote(exponential_drift),
            extra = coefs[free_extra],
            extra_lower = vapply(extra_bounds, `[[`, numeric(1), 1L),
            extra_upper = vapply(extra_bounds, `[[`, numeric(1), 2L),
            floor_params = c("D", "tau"),
            fn = quote(SSexponential_drift),
            fix = fix_used,
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_coefs))
        }
        model <- enforced$model
        coefs <- enforced$coefs
        fitted_vals <- stats::predict(model)

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] else NULL
        TD_val <- TD_arg %||% NA_real_
        MRT_val <- sum(TD_arg, coefs[["tau"]])
        HRT_val <- sum(TD_arg, coefs[["tau"]] * log(2))

        ## a tied texc is reported from the fitted coefficients
        texc_val <- if ("texc" %in% names(coefs)) {
            coefs[["texc"]]
        } else {
            sum(TD_arg, drift_k * coefs[["tau"]])
        }

        ## predict response at MRT and HRT using the full fitted model
        fitted_params <- exponential_drift(
            t = c(MRT_val, HRT_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            slope = coefs[["slope"]],
            texc = texc_val,
            TD = TD_arg
        )

        list(
            coefs = data.frame(
                A = coefs[["A"]],
                B = coefs[["B"]],
                tau = coefs[["tau"]],
                k = 1 / coefs[["tau"]], ## time_channel units^-1
                TD = TD_val,
                MRT = MRT_val,
                HRT = HRT_val,
                MRT_fitted = fitted_params[[1L]],
                HRT_fitted = fitted_params[[2L]],
                slope = coefs[["slope"]],
                texc = texc_val
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
                env = env
            )
        )
    }

    return(analyse_kinetics_channels(
        data,
        nirs_channels,
        time_channel,
        per_channel,
        expdrift_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
