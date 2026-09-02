#' Biexponential function
#'
#' Calculate a two-phase curve: a fast monoexponential response toward `B`
#' and a slow monoexponential response from `B` toward a stable plateau at
#' `B2`, both clocked from the response onset and summed.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting value of the response
#'   variable (the `t = 0` intercept).
#' @param B A numeric parameter for the asymptote of the *fast* component;
#'   the value the fast response alone would approach.
#' @param tau A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`. Dominates the initial steep
#'   response.
#' @param B2 A numeric parameter for the asymptote of the *slow* component;
#'   the stable plateau the response recovers toward as `t` approaches
#'   infinity.
#' @param tau2 A numeric parameter for the *slow* time constant (\eqn{\tau_2}),
#'   in units of the predictor variable `t`. Typically `tau2 >> tau`.
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
#'   `A + (B - A) * (1 - exp(-t / tau)) + (B2 - B) * (1 - exp(-t / tau2))`
#'
#' 6-parameter model (with time delay), where `ts = pmax(t - TD, 0)`:
#'   `A + (B - A) * (1 - exp(-ts / tau)) +
#'   (B2 - B) * (1 - exp(-ts / tau2))`
#'
#' `A`, `B`, and `B2` are all values on the response scale. The fast
#' component is a [monoexponential()] response from `A` toward `B` with
#' amplitude `B - A`; the slow component runs concurrently from the same
#' onset with amplitude `B2 - B`. The curve starts at `A`, approaches `B2`
#' as `t` grows, and is smooth throughout.
#'
#' The expected response is a *fast* excursion toward a minimum or maximum
#' short of `B`, followed by a *slow* recovery back to a stable plateau at
#' `B2`. The turning point occurs where the two phase rates cancel:
#' `texc = TD + log(r) / (1 / tau - 1 / tau2)` with
#' `ratio = -(B - A) * tau2 / ((B2 - B) * tau)`, which exists only when
#' the amplitudes oppose in sign and the fast phase dominates at the onset
#' (`ratio > 1`). If `B` is between `A` and `B2`, the response is monotonic
#' but still two-phase. If `B = B2`, the curve reduces to a
#' [monoexponential()] with single time constant `tau` and asymptote `B2`.
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
#' x <- biexponential(t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B, tau, B2, tau2),
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
biexponential <- function(t, A, B, tau, B2, tau2, TD = NULL) {
    ## both phases clocked from the onset; the slow term carries the
    ## response from B toward B2
    ts <- if (is.null(TD)) t else pmax(t - TD, 0)
    return(
        A + (B - A) * (1 - exp(-ts / tau)) + (B2 - B) * (1 - exp(-ts / tau2))
    )
}


## time of the curve turning point, elapsed from the fit origin; NA when
## the amplitudes share a sign or the slow phase dominates from the onset
## (monotonic response), or the time constants coincide
biexp_texc <- function(A, B, tau, B2, tau2, TD = NULL) {
    r <- -((B - A) * tau2) / ((B2 - B) * tau)
    if (!is.finite(r) || r <= 1 || tau == tau2) {
        return(NA_real_)
    }
    return(sum(TD, log(r) / (1 / tau - 1 / tau2)))
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
    fixed <- list(...)$fixed %||% list()
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    return(biexp_start(tx[["y"]], tx[["x"]], fixed, "TD" %in% names(mCall)))
}


## biexponential phase separation: the largest admissible tau / tau2,
## shared by the start grid and the fit bounds
tau_ratio <- 0.98


#' Grid-profiled starting estimates for the biexponential model
#'
#' Vector-level initialiser behind [biexp_init()], called directly by the
#' kinetics worker with `tau` and `TD` held at their stage-1 values to
#' seed the slow phase. Profiles the time constants (and
#' `TD`) on a coarse grid and keeps the RSS-minimising start (cf.
#' [expdrift_start()]). The model is linear in `A`, `B`, and `B2` once
#' `tau`, `tau2`, and `TD` are held, so those are solved by least squares
#' at every grid point at once: the Gram entries of the bases `e1`,
#' `e2 - e1`, `1 - e2` for every `(tau, tau2)` pair follow from the
#' column products of the two exponential matrices, and [solve_grid3()]
#' solves the pairs in one pass. User-fixed values narrow the grids; the
#' amplitudes are always solved free, as this is only a seed. Pairs with
#' `tau / tau2 > 0.98` are dropped as their bases are near-collinear,
#' unless both time constants are fixed.
#'
#' @inheritParams monoexp_start
#'
#' @returns A named numeric vector of starting estimates in model order.
#'
#' @keywords internal
biexp_start <- function(x, t, fixed = list(), has_TD = FALSE) {
    n <- length(t)
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        span <- 1
    }
    tau_grid <- fixed$tau %||%
        exp(seq(log(span / 100), log(span / 2), length.out = 13L))
    tau2_grid <- fixed$tau2 %||%
        exp(seq(log(span / 20), log(span * 10), length.out = 9L))
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, span / 3, length.out = 11L)
    }
    n1 <- length(tau_grid)
    n2 <- length(tau2_grid)
    ok <- outer(tau_grid, tau2_grid, \(.a, .b) .b >= .a / tau_ratio)
    if (!is.null(fixed$tau) && !is.null(fixed$tau2)) {
        ok[] <- TRUE
    }
    ## expand tau- and tau2-indexed vectors over the (tau, tau2) grid
    by1 <- \(v) matrix(v, n1, n2)
    by2 <- \(v) matrix(v, n1, n2, byrow = TRUE)

    ## the response is centred for conditioning; the constant is carried
    ## by the bases (they sum to one), so the asymptotes shift back
    xm <- mean(x)
    xc <- x - xm
    sx <- sum(xc)
    xx <- sum(xc^2)
    blocks <- lapply(td_grid, \(.td) {
        ts <- if (has_TD) pmax(t - .td, 0) else t
        E1 <- exp(-outer(ts, tau_grid, `/`))
        E2 <- exp(-outer(ts, tau2_grid, `/`))
        P <- crossprod(E1, E2)
        d1 <- colSums(E1^2)
        d2 <- colSums(E2^2)
        s1 <- colSums(E1)
        s2 <- colSums(E2)
        xe1 <- drop(crossprod(E1, xc))
        xe2 <- drop(crossprod(E2, xc))
        fit <- solve_grid3(
            g11 = by1(d1),
            g12 = P - by1(d1),
            g13 = by1(s1) - P,
            g22 = by1(d1) + by2(d2) - 2 * P,
            g23 = P - by1(s1) + by2(s2 - d2),
            g33 = by2(n - 2 * s2 + d2),
            b1 = by1(xe1),
            b2 = by2(xe2) - by1(xe1),
            b3 = sx - by2(xe2),
            xx = xx
        )
        fit$rss[!ok] <- Inf
        fit
    })
    k <- which.min(vapply(blocks, \(.b) min(.b$rss), numeric(1)))
    b <- blocks[[k]]
    ij <- arrayInd(which.min(b$rss), dim(b$rss))
    if (!is.finite(b$rss[ij])) {
        stop("No starting estimates could be resolved from the response.")
    }

    return(c(
        A = b$c1[ij] + xm,
        B = b$c2[ij] + xm,
        tau = tau_grid[[ij[[1L]]]],
        B2 = b$c3[ij] + xm,
        tau2 = tau2_grid[[ij[[2L]]]],
        TD = if (has_TD) td_grid[[k]]
    ))
}


#' Biexponential model with gradient
#'
#' [biexp_core()] evaluates the curve and its partial derivatives on the
#' canonical parameters. [biexp_model()] is the model function of
#' [SSbiexponential()]: [biexponential()] plus the gradient for the
#' parameters written as bare symbols in the call (see [free_params()]),
#' so [stats::nls()] skips [stats::numericDeriv()].
#'
#' @inheritParams biexponential
#'
#' @returns [biexp_core()]: a list of the curve `val` and the partial
#'   derivatives by parameter name. [biexp_model()]: a numeric vector of
#'   predicted values with a `"gradient"` attribute when any parameter is
#'   free.
#'
#' @keywords internal
biexp_core <- function(t, A, B, tau, B2, tau2, TD = NULL) {
    has_TD <- !is.null(TD)
    ts <- if (has_TD) pmax(t - TD, 0) else t
    e1 <- exp(-ts / tau)
    e2 <- exp(-ts / tau2)
    return(list(
        val = A + (B - A) * (1 - e1) + (B2 - B) * (1 - e2),
        A = e1,
        B = e2 - e1,
        tau = -(B - A) * e1 * ts / tau^2,
        B2 = 1 - e2,
        tau2 = -(B2 - B) * e2 * ts / tau2^2,
        TD = if (has_TD) {
            -(t > TD) * ((B - A) * e1 / tau + (B2 - B) * e2 / tau2)
        }
    ))
}


#' @rdname biexp_core
#' @keywords internal
biexp_model <- function(t, A, B, tau, B2, tau2, TD = NULL) {
    g <- biexp_core(t, A, B, tau, B2, tau2, TD)
    val <- g$val
    free <- free_params(
        match.call(),
        c("A", "B", "tau", "B2", "tau2", if (!is.null(TD)) "TD")
    )
    if (length(free) > 0L) {
        attr(val, "gradient") <- do.call(cbind, g[free])
    }
    return(val)
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. Supports both the
#' 5-parameter form (A, B, tau, B2, tau2) and the 6-parameter form adding
#' a time delay TD; arity is inferred from the formula passed to
#' [stats::nls()].
#'
#' @usage
#' SSbiexponential(t, A, B, tau, B2, tau2, TD)
#'
#' @inheritParams biexponential
#'
#' @details
#' 5-parameter model:
#' `x ~ SSbiexponential(t, A, B, tau, B2, tau2)`
#'
#' 6-parameter model:
#' `x ~ SSbiexponential(t, A, B, tau, B2, tau2, TD)`
#'
#' The two phases are weakly identified when `tau` and `tau2` are close, so
#' `algorithm = "port"` with the time constants bounded non-negative and
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#' [analyse_kinetics()] instead fits the phases sequentially, holding the
#' fast phase near a monoexponential estimate.
#'
#' The model function returns the analytic gradient for the free
#' parameters as a `"gradient"` attribute, so [stats::nls()] does not
#' resort to [stats::numericDeriv()] and [stats::predict()] on a fitted
#' model carries the attribute; drop it with `as.vector()`.
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
#'   `x ~ SSbiexponential(t, A, B, tau = 5, B2, tau2)`
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
#' x <- biexponential(t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSbiexponential(t, A, B, tau, B2, tau2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' ## fix the fast time constant
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A, B, tau = 5, B2, tau2),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model_fixed)
#'
#' @export
SSbiexponential <- selfStart(
    model = biexp_model,
    initial = init_fixed(
        biexp_init,
        c("A", "B", "tau", "B2", "tau2", "TD")
    ),
    parameters = c("A", "B", "tau", "B2", "tau2", "TD")
)


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' excursion-recovery curve to each `nirs_channel` within a single *"mnirs"*
#' data frame in two stages: the fast phase as a monoexponential on the
#' `end_window` window ([fit_monoexponential()]; `Inf` resolves to 30
#' time units past the first extreme), then the full
#' [SSbiexponential()] model on the whole response with `A`, `tau`, and
#' `TD` box-bounded about their stage-1 values and `B`, `B2`, `tau2`
#' free. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit the fast phase with a
#'   time delay, giving a 6-parameter [SSbiexponential()] model (A, B,
#'   tau, B2, tau2, TD). If that fit fails, or if `use_TD = FALSE`, the
#'   reduced 5-parameter model without `TD` is fit.
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `tau`, `B2`, `tau2`, `TD`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation and
#'   reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. `TD` is fixable for channels where
#'   `use_TD = TRUE`; a fixed `TD` disables the 5-parameter fallback.
#' @param tau_flex Numeric; multiplicative half-width of the stage-2
#'   `tau` bounds about the stage-1 value,
#'   `tau * [1 / (1 + tau_flex), 1 + tau_flex]`. `tau2` is floored at
#'   the `tau` ceiling divided by `0.98` and capped at ten times the span.
#' @param TD_flex Numeric; additive half-width of the stage-2 `TD` bounds
#'   in units of `time_channel`, floored at `0`.
#' @param A_flex Numeric; additive half-width of the stage-2 `A` bounds on
#'   the response scale. `NULL` (*default*) uses twice the stage-1
#'   residual standard deviation.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `TD`, `tau`, `MRT`, `texc`, `B2`, `tau2`,
#'   `MRT_fitted`, `texc_fitted`. Per-channel metadata are attached as
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
    fix = NULL,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    tau_flex = 1 / 3,
    TD_flex = 2,
    A_flex = NULL,
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
        arg_list = mget(c(
            "use_TD", "fix", "start_time", "direction", "end_window",
            "tau_flex", "TD_flex", "A_flex"
        )),
        choices = list(direction = c("auto", "positive", "negative")),
        ## TD is only fixable where that channel fits the 6-parameter model
        fix_params = \(.a) {
            c("A", "B", "tau", "B2", "tau2", if (.a$use_TD) "TD")
        },
        verbose = verbose,
        env = env
    )

    time_channel <- setup$time_channel
    ## `end_window` bounds the stage-1 fast phase only, so the global
    ## default `Inf` (the whole response) is replaced by 30 time units past
    ## the first extreme
    per_channel <- lapply(setup$per_channel, \(.a) {
        if (is.infinite(.a$end_window)) {
            .a$end_window <- 30
        }
        .a
    })
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- kinetics_coef_cols$biexponential

    ## method-specific fit in two stages. stage 1: the fast phase as a
    ## monoexponential on the `end_window` window (`x_fit`, `t_fit`).
    ## stage 2: the biexponential on the full response, with A, tau, and
    ## TD box-bounded around their stage-1 values by the `*_flex`
    ## half-widths; B, B2, and tau2 free. a failed stage returns NA, and
    ## the fallback chain resolves the row upstream
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        fix <- .a$fix %||% list()

        ## stage 1: fixed parameters shared with the fast phase carry over
        a1 <- .a
        a1$fix <- keep_fix(fix, c("A", "B", "tau", "TD"))
        # fmt: skip
        fast <- fit_monoexponential(
            .nirs, x_fit, t_fit, a1, valid, time_channel, interval_name, env
        )
        if (is.null(fast$model)) {
            return(build_na_results(na_cols))
        }
        cf1 <- fast$coefs
        has_TD <- is.finite(cf1$TD)
        params <- c("A", "B", "tau", "B2", "tau2", if (has_TD) "TD")
        free <- setdiff(params, names(fix))

        ## stage 2 window: the full response; the TD model is flat at A
        ## before TD so the pre-onset baseline is kept, as in
        ## `fit_td_fallback()`
        t_rel <- data[[time_channel]] - .a$start_time
        idx <- which(is.finite(data[[.nirs]]) & is.finite(t_rel))
        x_full <- data[[.nirs]][idx]
        t_full <- t_rel[idx]
        keep <- has_TD | t_full >= 0
        fit_data <- setNames(
            data.frame(x_full[keep], t_full[keep]),
            fit_names(.nirs, time_channel, params)
        )
        span <- diff(range(fit_data[[2L]]))
        on_error <- \(e) {
            warn_fit_failed(
                quote(SSbiexponential),
                e,
                .nirs,
                interval_name,
                length(params),
                env = env
            )
            NULL
        }
        if (nrow(fit_data) <= length(free)) {
            on_error(simpleError(sprintf(
                "%d observations for %d free parameters.",
                nrow(fit_data),
                length(free)
            )))
            return(build_na_results(na_cols))
        }

        ## fast phase from stage 1 (user-fixed values already merged in)
        prior <- c(
            list(A = cf1$A, tau = cf1$tau),
            if (has_TD) list(TD = cf1$TD)
        )
        A_flex <- .a$A_flex %||% (2 * stats::sd(stats::residuals(fast$model)))
        tau_flex <- .a$tau_flex
        lower <- c(
            A = prior$A - A_flex,
            B = -Inf,
            tau = prior$tau / (1 + tau_flex),
            B2 = -Inf,
            ## the slow phase separates above the fast-phase ceiling
            tau2 = prior$tau * (1 + tau_flex) / tau_ratio,
            TD = if (has_TD) max(0, prior$TD - .a$TD_flex)
        )
        upper <- c(
            A = prior$A + A_flex,
            B = Inf,
            tau = prior$tau * (1 + tau_flex),
            B2 = Inf,
            ## a slow tail far beyond the record identifies only its rate
            tau2 = 10 * span,
            TD = if (has_TD) prior$TD + .a$TD_flex
        )
        ## B, B2, tau2 seeded by the start grid with the fast phase held
        model <- tryCatch(
            {
                start <- biexp_start(
                    fit_data[[1L]],
                    fit_data[[2L]],
                    utils::modifyList(fix, prior[names(prior) != "A"]),
                    has_TD
                )
                start[names(prior)] <- unlist(prior)
                suppressWarnings(nls(
                    build_ss_formula(
                        quote(SSbiexponential),
                        params,
                        fix,
                        names(fit_data)[[1L]],
                        names(fit_data)[[2L]]
                    ),
                    fit_data,
                    start = pmin(pmax(start[free], lower[free]), upper[free]),
                    algorithm = "port",
                    lower = lower[free],
                    upper = upper[free],
                    control = stats::nls.control(
                        maxiter = 500L,
                        warnOnly = TRUE
                    )
                ))
            },
            error = on_error
        )
        model <- accept_port_fit(model, on_error)
        if (is.null(model)) {
            return(build_na_results(na_cols))
        }
        coefs <- full_coefs(model, params, fix)

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if (has_TD) coefs[["TD"]] else NULL
        ## fast-phase mean response time, as for the monoexponential
        MRT_val <- sum(TD_arg, coefs[["tau"]])

        ## excursion time (texc) is the fitted turning point, reported
        ## elapsed from start_time, mirroring MRT = TD + tau; NA when the
        ## fitted response is monotonic
        texc_val <- biexp_texc(
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            B2 = coefs[["B2"]],
            tau2 = coefs[["tau2"]],
            TD = TD_arg
        )
        ## predict response at MRT and texc using the full fitted model; an
        ## NA texc (monotonic fit) propagates to NA
        fitted_params <- biexponential(
            t = c(MRT_val, texc_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            B2 = coefs[["B2"]],
            tau2 = coefs[["tau2"]],
            TD = TD_arg
        )

        build_fit_results(
            data.frame(
                A = coefs[["A"]],
                B = coefs[["B"]],
                TD = TD_arg %||% NA_real_,
                tau = coefs[["tau"]],
                MRT = MRT_val,
                texc = texc_val,
                B2 = coefs[["B2"]],
                tau2 = coefs[["tau2"]],
                MRT_fitted = fitted_params[[1L]],
                texc_fitted = fitted_params[[2L]]
            ),
            model,
            x_full,
            t_full,
            utils::modifyList(valid, list(idx = idx)),
            keep,
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
