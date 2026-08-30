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
    fixed <- list(...)$fixed %||% list()
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    return(biexp_start(tx[["y"]], tx[["x"]], fixed, "TD" %in% names(mCall)))
}



## biexponential phase separation: the largest admissible tau1 / tau2,
## shared by the start grid and the fit bounds
tau_ratio <- 0.98


#' Grid-profiled starting estimates for the biexponential model
#'
#' Vector-level initialiser behind [biexp_init()], called directly by the
#' kinetics worker on the fit window. Profiles the time constants (and
#' `TD`) on a coarse grid and keeps the RSS-minimising start (cf.
#' [expdrift_start()]). The model is linear in `A`, `B1`, and `B2` once
#' `tau1`, `tau2`, and `TD` are held, so those are solved by least squares
#' at every grid point at once: the Gram entries of the bases `e1`,
#' `e2 - e1`, `1 - e2` for every `(tau1, tau2)` pair follow from the
#' column products of the two exponential matrices, and [solve_grid3()]
#' solves the pairs in one pass. User-fixed values narrow the grids; the
#' amplitudes are always solved free, as this is only a seed. Pairs with
#' `tau1 / tau2 > 0.98` are dropped as their bases are near-collinear,
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
    tau1_grid <- fixed$tau1 %||%
        exp(seq(log(span / 100), log(span / 2), length.out = 13L))
    tau2_grid <- fixed$tau2 %||%
        exp(seq(log(span / 20), log(span * 10), length.out = 9L))
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, span / 3, length.out = 11L)
    }
    n1 <- length(tau1_grid)
    n2 <- length(tau2_grid)
    ok <- outer(tau1_grid, tau2_grid, \(.a, .b) .b >= .a / tau_ratio)
    if (!is.null(fixed$tau1) && !is.null(fixed$tau2)) {
        ok[] <- TRUE
    }
    ## expand tau1- and tau2-indexed vectors over the (tau1, tau2) grid
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
        E1 <- exp(-outer(ts, tau1_grid, `/`))
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
        B1 = b$c2[ij] + xm,
        tau1 = tau1_grid[[ij[[1L]]]],
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
#' so [stats::nls()] skips [stats::numericDeriv()]. [biexp_ratio()] is
#' the internal fitting model on `(r = tau1 / tau2, tau2)`: the phase
#' separation becomes the box bound `r <= 0.98`, closing the
#' non-identifiable valley at `tau1 -> tau2` that stalls the canonical
#' fit with a runaway `B1`.
#'
#' @param r A numeric parameter for the time-constant ratio `tau1 / tau2`.
#' @inheritParams biexponential
#'
#' @returns [biexp_core()]: a list of the curve `val` and the partial
#'   derivatives by parameter name. [biexp_model()], [biexp_ratio()]: a
#'   numeric vector of predicted values with a `"gradient"` attribute when
#'   any parameter is free.
#'
#' @keywords internal
biexp_core <- function(t, A, B1, tau1, B2, tau2, TD = NULL) {
    has_TD <- !is.null(TD)
    ts <- if (has_TD) pmax(t - TD, 0) else t
    e1 <- exp(-ts / tau1)
    e2 <- exp(-ts / tau2)
    return(list(
        val = A + (B1 - A) * (1 - e1) + (B2 - B1) * (1 - e2),
        A = e1,
        B1 = e2 - e1,
        tau1 = -(B1 - A) * e1 * ts / tau1^2,
        B2 = 1 - e2,
        tau2 = -(B2 - B1) * e2 * ts / tau2^2,
        TD = if (has_TD) {
            -(t > TD) * ((B1 - A) * e1 / tau1 + (B2 - B1) * e2 / tau2)
        }
    ))
}


#' @rdname biexp_core
#' @keywords internal
biexp_model <- function(t, A, B1, tau1, B2, tau2, TD = NULL) {
    g <- biexp_core(t, A, B1, tau1, B2, tau2, TD)
    val <- g$val
    free <- free_params(
        match.call(),
        c("A", "B1", "tau1", "B2", "tau2", if (!is.null(TD)) "TD")
    )
    if (length(free) > 0L) {
        attr(val, "gradient") <- do.call(cbind, g[free])
    }
    return(val)
}


#' @rdname biexp_core
#' @keywords internal
biexp_ratio <- function(t, A, B1, r, B2, tau2, TD = NULL) {
    g <- biexp_core(t, A, B1, r * tau2, B2, tau2, TD)
    val <- g$val
    free <- free_params(
        match.call(),
        c("A", "B1", "r", "B2", "tau2", if (!is.null(TD)) "TD")
    )
    if (length(free) > 0L) {
        ## chain rule through tau1 = r * tau2
        g$r <- g$tau1 * tau2
        g$tau2 <- g$tau2 + g$tau1 * r
        attr(val, "gradient") <- do.call(cbind, g[free])
    }
    return(val)
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
#' [analyse_kinetics()] additionally fits on the ratio `tau1 / tau2` to
#' hold the phases apart.
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
    model = biexp_model,
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

    ## the phases separate at tau1 / tau2 <= `r_max`, matching the start
    ## grid; a fit pinned at the bound wants the non-identifiable
    ## tau1 = tau2 valley and is reported as such
    r_max <- tau_ratio
    pinned <- \(cf) isTRUE(cf[["r"]] >= r_max * (1 - 1e-6))
    not_separable <- "Fast and slow phases are not separable (tau1 approaches tau2)."

    ## canonical <-> ratio coefficient vectors, keeping parameter order
    to_ratio <- \(cf) {
        cf[["tau1"]] <- cf[["tau1"]] / cf[["tau2"]]
        names(cf) <- sub("^tau1$", "r", names(cf))
        cf
    }
    to_canon <- \(cf) {
        cf[["r"]] <- cf[["r"]] * cf[["tau2"]]
        names(cf) <- sub("^r$", "tau1", names(cf))
        cf
    }

    ## fit bounds. ratio space: taus floored as degeneracy markers, r
    ## capped at the separation, tau2 capped so a runaway slow component
    ## cannot diverge with B2. canonical space: the separation is the box
    ## bound tau2 >= tau1 / r_max when tau1 is user-fixed, else only the floor
    ## (the canonical refit starts at a separated optimum)
    bounds <- \(span, ratio, fix) {
        if (ratio) {
            list(
                lower = c(r = 1e-6, tau2 = span * 1e-6, TD = 0),
                upper = c(r = r_max, tau2 = 10 * span)
            )
        } else {
            list(
                lower = c(
                    tau1 = span * 1e-6,
                    tau2 = max(span * 1e-6, fix$tau1 / r_max),
                    TD = 0
                ),
                upper = c(tau2 = 10 * span)
            )
        }
    }
    ## named box over the free parameters, unbounded where unset
    box <- \(bnd, free, fill) {
        b <- setNames(bnd[free], free)
        b[is.na(b)] <- fill
        b
    }
    port_fit <- \(formula, .data, start, bnd, on_error) {
        lower <- box(bnd$lower, names(start), -Inf)
        upper <- box(bnd$upper, names(start), Inf)
        tryCatch(
            suppressWarnings(nls(
                formula,
                .data,
                ## a grid edge can overshoot its cap by rounding
                start = pmin(pmax(start, lower), upper),
                algorithm = "port",
                lower = lower,
                upper = upper,
                control = stats::nls.control(maxiter = 500L, warnOnly = TRUE)
            )),
            error = on_error
        )
    }

    ## method-specific fit: biexponential via nls port on the ratio
    ## parameterisation, or canonical when tau1 is user-fixed (r is then
    ## not a box bound); a failed 6-param fit falls back to the 5-param
    ## model. the returned model is always canonical
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        fix <- .a$fix
        ratio <- !"tau1" %in% names(fix)
        space <- \(.params) if (ratio) sub("^tau1$", "r", .params) else .params

        ## port often stops short of its certificate on usable
        ## coefficients, which are kept with a warning
        fitter <- \(.data, .params, on_error) {
            params <- space(.params)
            free <- setdiff(params, names(fix))
            start <- tryCatch(
                biexp_start(.data[[1L]], .data[[2L]], fix, "TD" %in% .params),
                error = on_error
            )
            if (is.null(start)) {
                return(NULL)
            }
            if (ratio) {
                start <- to_ratio(start)
            }
            formula <- build_ss_formula(
                if (ratio) quote(biexp_ratio) else quote(SSbiexponential),
                params,
                fix,
                names(.data)[[1L]],
                names(.data)[[2L]]
            )
            bnd <- bounds(diff(range(.data[[2L]])), ratio, fix)
            model <- accept_port_fit(
                port_fit(formula, .data, start[free], bnd, on_error),
                on_error
            )
            ## a property of the data, so the reduced model is not retried
            if (!is.null(model) && ratio && pinned(coef(model))) {
                return(on_error(fit_final_error(not_separable)))
            }
            model
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
        coefs <- full_coefs(fit$model, space(params), fix)
        span <- diff(range(fit$data[[time_channel]]))
        bnd <- bounds(span, ratio, fix)

        ## enforce direction: bounded refit on the fast-phase amplitude
        ## D = B1 - A when inverted. the fit bounds are carried over; TD
        ## >= 0 and the separation are structural, so only D and the
        ## floors mark a degenerate fit
        free <- setdiff(space(params), names(fix))
        enforced <- enforce_direction(
            fit$model,
            coefs,
            fit$data,
            direction = .a$direction,
            amp_fn = if (ratio) quote(biexp_ratio) else quote(SSbiexponential),
            fn = "SSbiexponential",
            lower = bnd$lower[intersect(names(bnd$lower), free)],
            upper = bnd$upper[intersect(names(bnd$upper), free)],
            floor_params = c("D", if (ratio) c("r", "tau2")),
            fix = fix,
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_cols))
        }
        model <- enforced$model
        coefs <- enforced$coefs

        ## canonical re-expression from the ratio optimum, so the returned
        ## model reports (tau1, tau2): an interior stationary point
        ## converges at once, a pinned one is not separable
        if (ratio) {
            degenerate <- pinned(coefs)
            model <- if (!degenerate) {
                accept_port_fit(
                    port_fit(
                        build_ss_formula(
                            quote(SSbiexponential),
                            params,
                            fix,
                            names(fit$data)[[1L]],
                            names(fit$data)[[2L]]
                        ),
                        fit$data,
                        to_canon(coefs)[setdiff(params, names(fix))],
                        bounds(span, FALSE, fix),
                        \(e) NULL
                    ),
                    \(e) NULL
                )
            }
            if (is.null(model)) {
                warn_fit_failed(
                    quote(SSbiexponential),
                    simpleError(if (degenerate) {
                        not_separable
                    } else {
                        "Canonical re-expression of the fit failed."
                    }),
                    .nirs,
                    interval_name,
                    length(params),
                    env = env
                )
                return(build_na_results(na_cols))
            }
            coefs <- full_coefs(model, params, fix)
        }

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
            model,
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
