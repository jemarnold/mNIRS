#' Biexponential function
#'
#' Calculate a 5- or 6-parameter biexponential curve: a fast component that
#' carries the response to an excursion point, followed by a slow component
#' that recovers toward a stable plateau.
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
#'   `A + (B1 - A) * (1 - exp(-ts / tau1)) + (B2 - B1) * (1 - exp(-ts / tau2))`
#'
#' `A`, `B1`, and `B2` are all values on the response scale. The fast component
#' `B1` & `tau1` drives the steep initial response away from `A` to the
#' excursion point. The slow component `B2` & `tau2` then carries the curve
#' toward the plateau `B2` as `t` approaches infinity.
#'
#' The expected response is a *fast* overshoot to a minimum or maximum at the
#' excursion point near (but not exactly equal to) `B1`, followed by a *slow*
#' recovery back to a stable plateau at `B2`. The exact excursion point is
#' reported by [analyse_kinetics()] as `texc_fitted` at `texc`.
#' If `B1` is between `A` and `B2`, the response curve will be monotonic but
#' still two-phase. If `B1 = B2`, the curve reduces to a [monoexponential()]
#' with single time constant `tau` and asymptote `B2`.
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
#' x <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40) +
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
        y <- A +
            (B1 - A) * (1 - exp(-t / tau1)) +
            (B2 - B1) * (1 - exp(-t / tau2))
    } else {
        ## 6-parameter: with time delay
        ts <- pmax(t - TD, 0)
        y <- A +
            (B1 - A) * (1 - exp(-ts / tau1)) +
            (B2 - B1) * (1 - exp(-ts / tau2))
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
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()
    has_TD <- "TD" %in% names(mCall)
    span <- diff(range(t))

    ## profile the amplitudes out of a joint grid of time constants and
    ## candidate delays; the RSS-minimising delay is the TD seed
    td_grid <- fixed$TD %||%
        if (has_TD) seq(0, span / 3, length.out = 21L) else 0
    grid <- biexp_grid_start(x, t, TD = td_grid)
    if (is.null(grid)) {
        cli_abort(c(
            "x" = "No starting estimates could be resolved from the response."
        ))
    }
    TD_init <- grid$TD

    ## user-fixed values take precedence over the grid optimum. the ratio
    ## starts strictly inside the tau2 / tau1 >= 2.5 bound of the grid
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
#'   bounds, including a lower bound on `lr` set by its `tau_ratio` argument.
#'   Resulting model parameters may differ from [SSbiexponential()].
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
#' x <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B1, lt1, B2, lr), data = data)
#' summary(model)
#'
#' y <- predict(model, data)
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


#' Translate natural-scale fixed parameters to the log-ratio scale
#'
#' A fixed `tau1` becomes a constant `lt1`; a fixed `tau2` becomes an
#' `lr` expression in `lt1` (a constant when `tau1` is also fixed),
#' enforcing the ratio bound from the other side.
#'
#' @param fix Named list of fixed parameter values on the natural scale.
#'
#' @returns A named list of fixed values on the log-ratio scale; `lr`
#'   may be a language object.
#'
#' @keywords internal
biexp_fix_ss <- function(fix) {
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
    return(fix_ss)
}


## natural-scale -> log-ratio-scale parameter names
biexp_params_ss <- c(
    A = "A",
    B1 = "B1",
    tau1 = "lt1",
    B2 = "B2",
    tau2 = "lr",
    TD = "TD"
)


#' Structural box bounds on the log-ratio scale
#'
#' Single source for the `algorithm = "port"` bounds shared by
#' [fit_biexp_ratio()] and the [enforce_direction()] refit: `lt1` spans
#' the data down to `span / 1000` and is capped by a fixed `tau2`, `lr`
#' is bounded below by `tau_ratio`, and `TD` lies within the data span.
#'
#' @param span A numeric time span of the fitted data.
#' @param tau_ratio A numeric lower bound on `tau2 / tau1`.
#' @param fix Named list of fixed parameter values on the natural scale.
#'
#' @returns A list of named numeric vectors `lower` and `upper`.
#'
#' @keywords internal
biexp_ss_bounds <- function(span, tau_ratio, fix = list()) {
    list(
        lower = c(lt1 = log(span / 1000), lr = log(tau_ratio), TD = 0),
        upper = c(
            lt1 = log((fix$tau2 %||% Inf) / tau_ratio),
            lr = Inf,
            TD = span
        )
    )
}


#' Fit a biexponential in ratio parameterisation
#'
#' [fit_biexp_ratio()]: Fits [biexponential()] with the time constants
#' re-parameterised as `lt1 = log(tau1)` and `lr = log(tau2 / tau1)`.
#'
#' @details
#' The fit uses [stats::nls()] with `algorithm = "port"` and box bounds. On
#' the log-ratio scale the `tau2 / tau1 >= tau_ratio` separation of the two
#' components is a simple lower bound on `lr`, which holds the fit away from
#' the degenerate `tau2 = tau1` limit. The asymptotes `A`, `B1`, and `B2` are
#' unbounded, so identifiability rests on the time constants alone.
#'
#' A non-converged fit is accepted with a warning when its coefficients are
#' finite and its residual sum of squares is no worse than the grid starting
#' estimates from [biexp_grid_start()]; otherwise it is rejected and `NULL`
#' is returned (see [accept_port_fit()]).
#'
#' A user-fixed `tau1` or `tau2` is substituted into the re-parameterised
#' formula rather than estimated: fixing `tau1` drops `lt1` and leaves `lr`
#' free about the fixed value; fixing `tau2` drops `lr` and caps `lt1` at
#' `log(tau2 / tau_ratio)`, enforcing the ratio bound from the other side.
#'
#' @param x,t Numeric vectors of the response and predictor variables.
#' @param params Character vector of parameter names in model order.
#' @param fix Named list of fixed parameter values.
#' @param tau_ratio A numeric lower bound on `tau2 / tau1`.
#' @param on_error A function called with the [stats::nls()] error condition,
#'   or with a warning condition when a non-converged fit is accepted;
#'   returns `NULL` (see [fit_td_fallback()]).
#'
#' @returns An [nls][stats::nls] model in ratio parameterisation, or `NULL`.
#'
#' @keywords internal
fit_biexp_ratio <- function(x, t, params, fix, tau_ratio, on_error) {
    has_TD <- "TD" %in% params
    span <- diff(range(t))

    ## joint grid over time constants and candidate delays seeds the fit;
    ## the RSS-minimising delay is the TD seed
    td_grid <- fix$TD %||%
        if (has_TD) seq(0, span / 3, length.out = 21L) else 0
    seed <- biexp_grid_start(x, t, tau_ratio, td_grid)
    if (is.null(seed)) {
        return(on_error(simpleError(
            "No starting estimates could be resolved from the response."
        )))
    }
    TD_seed <- seed$TD

    ## a fixed tau2 caps tau1 from above; the seed is pulled inside that cap
    tau1_cap <- (fix$tau2 %||% Inf) / tau_ratio
    tau1_seed <- fix$tau1 %||% min(seed$tau1, tau1_cap / 1.01)
    tau2_seed <- fix$tau2 %||% seed$tau2

    ## start, lower, upper per parameter, in model order. the ratio starts
    ## strictly inside its lower bound
    bnd <- biexp_ss_bounds(span, tau_ratio, fix)
    bounds <- list(
        A = c(seed$A, -Inf, Inf),
        B1 = c(seed$B1, -Inf, Inf),
        lt1 = c(log(tau1_seed), bnd$lower[["lt1"]], bnd$upper[["lt1"]]),
        B2 = c(seed$B2, -Inf, Inf),
        lr = c(
            log(max(tau2_seed / tau1_seed, tau_ratio * 1.01)),
            bnd$lower[["lr"]],
            bnd$upper[["lr"]]
        ),
        TD = c(TD_seed, bnd$lower[["TD"]], bnd$upper[["TD"]])
    )

    ## translate natural-scale params and fix to the SS log-ratio scale
    params_ss <- unname(biexp_params_ss[params])
    fix_ss <- biexp_fix_ss(fix)

    ## free parameters in formula order, keeping start/lower/upper aligned
    ## with the positional bounds expected by port
    free <- setdiff(params_ss, names(fix_ss))
    start <- vapply(bounds[free], `[[`, numeric(1), 1L)
    lower <- vapply(bounds[free], `[[`, numeric(1), 2L)
    upper <- vapply(bounds[free], `[[`, numeric(1), 3L)

    model <- tryCatch(
        suppressWarnings(nls(
            build_ss_formula(quote(SSbiexponential), params_ss, fix_ss),
            data.frame(.x = x, .t = t),
            start = start,
            algorithm = "port",
            lower = lower,
            upper = upper,
            control = stats::nls.control(maxiter = 500L, warnOnly = TRUE)
        )),
        error = on_error
    )

    ## a non-converged fit is kept when demonstrably good: an RSS no worse
    ## than the grid seed
    return(accept_port_fit(
        model,
        on_error,
        ok = stats::deviance(model) <= seed$rss
    ))
}


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a biexponential
#' excursion-recovery curve to each `nirs_channel` within a single *"mnirs"*
#' data frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_TD Logical; `TRUE` attempts to fit a 6-parameter
#'   [SSbiexponential()] model (A, B1, tau1, B2, tau2, TD) with a time delay.
#'   If the 6-parameter fit fails, or if `use_TD = FALSE`, attempts to fit a
#'   reduced 5-parameter [SSbiexponential()] model (A, B1, tau1, B2, tau2).
#'   The user-facing default (`TRUE`) is set by [analyse_kinetics()].
#' @param fix An *optional* named list of model parameters (`A`, `B1`, `tau1`,
#'   `B2`, `tau2`, `TD`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`. `TD` is fixable for channels where
#'   `use_TD = TRUE`; a fixed `TD` disables the 5-parameter fallback.
#' @param tau_ratio A numeric lower bound on the ratio of the slow to the fast
#'   time constant, `tau2 / tau1`; default is `2.5`. As `tau2` approaches
#'   `tau1` the two components become indistinguishable and fit convergence can
#'   fail, so the ratio is bounded away from that limit. Lower values allow
#'   closer time constants and larger, more strongly cancelling amplitudes.
#'   Model fit is often only weakly identified and settles on this bound, in
#'   which case the weak fit may be returned with a warning.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B1`, `tau1`, `B2`, `tau2`, `TD`,
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
    ## a ratio of 1 or less admits tau2 == tau1, where the two components
    ## are indistinguishable and the design is singular
    validate_numeric(
        tau_ratio, 1, c(1, Inf),
        inclusive = "right", msg1 = "one-element", env = env
    )

    ## NA scaffold (method columns only) for convergence failure
    na_cols <- c("A", "B1", "tau1", "B2", "tau2", "TD", "texc", "texc_fitted")

    ## method-specific fit: self-starting biexponential via nls; a failed
    ## 6-param fit falls back to the 5-param model
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## reads `.a$fix` at call time, so the tau2 cap below is honoured
        fitter <- \(.data, .params, on_error) {
            # fmt: skip
            fit_biexp_ratio(
                .data$.x, .data$.t, .params, .a$fix, tau_ratio, on_error
            )
        }
        fit <- fit_td_fallback(
            x_fit,
            t_fit,
            params = c("A", "B1", "tau1", "B2", "tau2", if (.a$use_TD) "TD"),
            .a,
            fitter,
            quote(SSbiexponential),
            .nirs,
            interval_name,
            env
        )
        if (is.null(fit$model)) {
            return(build_na_results(na_cols))
        }
        params <- fit$params
        span <- diff(range(fit$data$.t))

        ## tau2 far beyond the record is not identifiable -- only the rate
        ## (B2 - B1) / tau2 is, so tau2 and B2 diverge together at near-
        ## constant RSS. a runaway slow component is profiled at the
        ## horizon cap and refit through the fixed-tau2 pathway; a failed
        ## capped refit keeps the unconstrained fit
        if (is.null(.a$fix$tau2)) {
            cf <- stats::coef(fit$model)
            tau2_fit <- (.a$fix$tau1 %||% exp(cf[["lt1"]])) * exp(cf[["lr"]])
            if (tau2_fit > 10 * span) {
                .a$fix$tau2 <- 10 * span
                capped <- fit_td_fallback(
                    x_fit,
                    t_fit,
                    params,
                    .a,
                    fitter,
                    quote(SSbiexponential),
                    .nirs,
                    interval_name,
                    env,
                    retry = FALSE
                )
                if (is.null(capped$model)) {
                    .a$fix$tau2 <- NULL
                } else {
                    fit <- capped
                }
            }
        }
        model <- fit$model

        ## enforce direction: bounded refit on D = B2 - A when the overall
        ## amplitude is inverted. the refit stays on the log-ratio scale
        ## with the fit_biexp_ratio bounds carried over; those bounds are
        ## structural, so only the D sign floor marks a degenerate fit
        fix_ss <- biexp_fix_ss(.a$fix)
        fix_num <- fix_ss[vapply(fix_ss, is.numeric, logical(1))]
        coefs_ss <- full_coefs(model, unname(biexp_params_ss[params]), fix_num)
        free_extra <- setdiff(names(stats::coef(model)), c("A", "B2"))
        bnd <- biexp_ss_bounds(span, tau_ratio, .a$fix)

        ## a slow time constant stranded above the grid ceiling reseeds
        ## inside it, so the bounded refit starts from sane geometry
        extra <- coefs_ss[free_extra]
        if ("lr" %in% free_extra) {
            extra[["lr"]] <- max(
                min(extra[["lr"]], log(span * 10) - coefs_ss[["lt1"]]),
                log(tau_ratio * 1.01)
            )
        }
        enforced <- enforce_direction(
            model,
            coefs_ss,
            fit$data,
            direction = .a$direction,
            amp_fn = quote(biexponential_ratio),
            extra = extra,
            B_name = "B2",
            extra_lower = bnd$lower[intersect(names(bnd$lower), free_extra)],
            extra_upper = bnd$upper[intersect(names(bnd$upper), free_extra)],
            floor_params = "D",
            fn = quote(SSbiexponential),
            fix = fix_ss,
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_cols))
        }

        ## back-convert log-ratio time constants to the natural scale; a
        ## fixed time constant has no log-scale counterpart in the coefs
        coefs <- enforced$coefs
        coefs[["tau1"]] <- .a$fix$tau1 %||% exp(coefs[["lt1"]])
        coefs[["tau2"]] <- .a$fix$tau2 %||%
            (coefs[["tau1"]] * exp(coefs[["lr"]]))

        ## model parameters in biexponential() argument order
        pars <- as.list(coefs[c("A", "B1", "tau1", "B2", "tau2")])

        ## TD is already elapsed from start_time, matching the fit time base
        TD_arg <- if ("TD" %in% params) coefs[["TD"]] else NULL

        ## closed-form turning point where dy/dt = 0, elapsed since model
        ## onset (TD, or 0 for the 5-param model). a real root needs `B1`
        ## beyond both `A` and `B2`; a `B1` between them gives a monotone
        ## curve whose turning point sits at the onset boundary
        s <- with(pars, {
            ratio <- ((B2 - B1) / tau2) / ((A - B1) / tau1)
            root <- if (is.finite(ratio) && ratio > 0) {
                log(ratio) / (1 / tau2 - 1 / tau1)
            } else {
                NA_real_
            }
            if (is.finite(root) && root > 0) root else 0
        })

        ## excursion time (texc) is reported elapsed from start_time, mirroring
        ## MRT = TD + tau, so the onset-relative root is shifted by TD
        texc_val <- sum(TD_arg, s)
        texc_fitted_val <- do.call(
            biexponential,
            c(list(t = texc_val), pars, list(TD = TD_arg))
        )

        build_fit_results(
            data.frame(
                pars,
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
