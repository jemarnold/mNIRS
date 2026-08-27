#' Exponential-drift function
#'
#' Calculate a two-phase curve: a primary monoexponential response with a
#' secondary linear drift beginning near the asymptote.
#'
#' @param slope A numeric parameter for the linear drift rate `dx/dt`
#'   of the secondary phase, in response units per unit of the predictor
#'   variable `t`.
#' @param tau_mult A numeric multiple of `tau` after `TD` at which the linear
#'   drift begins: excursion point; `texc = TD + tau_mult * tau` (`TD = 0`
#'   when absent).
#' @inheritParams monoexponential
#'
#' @details
#' 5-parameter model:
#' `A + (B - A) * (1 - exp(-t / tau)) + slope * pmax(t - tau_mult * tau, 0)`
#'
#' 6-parameter model:
#' `A + (B - A) * (1 - exp(-pmax(t - TD, 0) / tau)) +
#' slope * pmax(t - TD - tau_mult * tau, 0)`
#'
#' The primary phase is a [monoexponential()] response toward the asymptote
#' `B`. The secondary linear drift is exactly zero before
#' `texc = TD + tau_mult * tau`.
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
#'     t, A = 10, B = 100, tau = 12, slope = -0.5, tau_mult = 3, TD = 15
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## the drift onset multiple is held constant in the formula
#' model <- nls(
#'     x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult = 3, TD),
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
exponential_drift <- function(t, A, B, tau, slope, tau_mult, TD = NULL) {
    ## primary monoexponential phase + hinge-linear secondary drift from
    ## the excursion point texc = TD + tau_mult * tau
    texc <- sum(TD, tau_mult * tau)
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
    ## user-fixed tau, TD, and tau_mult narrow the grids; the linear
    ## parameters are always solved free, as this is only a seed
    fixed <- list(...)$fixed %||% list()

    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    has_TD <- "TD" %in% names(mCall)
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        span <- 1
    }
    tau_mult <- fixed$tau_mult %||% 3

    ## profile tau (and TD) on a coarse grid and keep the RSS-minimising
    ## start (cf. `monoexp_init()`). the model is linear in A, B, and
    ## slope once tau and TD are held, so those are solved by least
    ## squares at every grid point. tau is capped so the drift onset
    ## stays inside the record; a grid point whose hinge has no support
    ## solves to NA and is skipped
    tau_grid <- fixed$tau %||%
        exp(seq(log(span / 100), log(span / tau_mult), length.out = 25L))
    td_grid <- if (!has_TD) {
        0
    } else {
        fixed$TD %||% seq(0, 0.5 * span, length.out = 21L)
    }
    grid <- expand.grid(tau = tau_grid, TD = td_grid)

    fits <- vapply(seq_len(nrow(grid)), \(.i) {
        ts <- if (has_TD) pmax(t - grid$TD[.i], 0) else t
        e <- exp(-ts / grid$tau[.i])
        X <- cbind(e, 1 - e, pmax(t - grid$TD[.i] - tau_mult * grid$tau[.i], 0))
        cf <- qr.coef(qr(X), x)
        c(cf, sum((x - X %*% cf)^2))
    }, numeric(4L))
    best <- which.min(fits[4L, ])

    return(c(
        A = fits[[1L, best]],
        B = fits[[2L, best]],
        tau = grid$tau[best],
        slope = fits[[3L, best]],
        tau_mult = tau_mult,
        TD = if (has_TD) grid$TD[best]
    ))
}


#' Self-starting exponential-drift model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [exponential_drift()], for use with [stats::nls()]. Supports both the
#' 5-parameter form (A, B, tau, slope, tau_mult) and the
#' 6-parameter form adding a time delay TD; arity is inferred from the
#' formula passed to [stats::nls()].
#'
#' @usage
#' SSexponential_drift(t, A, B, tau, slope, tau_mult, TD)
#'
#' @inheritParams exponential_drift
#'
#' @details
#' 5-parameter model:
#' `x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult)`
#'
#' 6-parameter model:
#' `x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult, TD)`
#'
#' The hinge at `texc = TD + tau_mult * tau` is not differentiable, so
#' `algorithm = "port"` with `tau` (and `TD`) bounded non-negative and
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult = 3)`
#'   holds the drift onset at `3 * tau`. Fixed parameters are excluded from
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
#'     t, A = 10, B = 100, tau = 12, slope = -0.5, tau_mult = 4, TD = 15
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## 6-parameter fit with the drift onset held at 4 tau
#' model <- nls(
#'     x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult = 4, TD),
#'     data = data,
#'     algorithm = "port",
#'     lower = c(-Inf, -Inf, 0, -Inf, 0),
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' @export
SSexponential_drift <- selfStart(
    model = exponential_drift,
    initial = init_fixed(
        expdrift_init,
        c("A", "B", "tau", "slope", "tau_mult", "TD")
    ),
    parameters = c("A", "B", "tau", "slope", "tau_mult", "TD")
)


#' Analyse exponential-drift kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "exponential_drift")`. Fits a two-phase
#' monoexponential + linear-drift curve to each `nirs_channel` within a
#' single *"mnirs"* data frame. See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param use_TD Logical; default is `TRUE` to attempt to fit a 6-parameter
#'   [SSexponential_drift()] model with a time delay. If the 6-parameter fit
#'   fails, or if `use_TD = FALSE`, attempts to fit a reduced 5-parameter
#'   model without `TD`.
#' @param tau_mult A numeric multiple of `tau` after `TD` at which the drift
#'   onset is held (`texc = TD + tau_mult * tau`; default is `3`,
#'   ~95% of the primary amplitude). Applied to every channel, or
#'   per-channel as a list keyed by channel name, e.g.
#'   `tau_mult = list(smo2 = 2)`.
#' @param fix An *optional* named list of model parameters (`A`, `B`, `tau`,
#'   `slope`, `TD`) to hold constant during fitting, e.g. `fix = list(A = 0)`.
#'   Applied to every channel, or per-channel as a list of lists keyed by
#'   channel name, e.g. `fix = list(smo2 = list(A = 0))`. `TD` is fixable
#'   for channels where `use_TD = TRUE`; a fixed `TD` disables the
#'   5-parameter fallback.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `tau`, `k`, `TD`, `MRT`, `HRT`,
#'   `MRT_fitted`, `HRT_fitted`, `slope`, `tau_mult`, `texc`, `texc_fitted`.
#'   Per-channel metadata are attached as attributes:
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
    tau_mult = 3,
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
        ## TD is only fixable where that channel fits the 6-parameter model
        fix_params = \(.a) c("A", "B", "tau", "slope", if (.a$use_TD) "TD"),
        verbose = verbose,
        env = env
    )
    ## a channel omitted from a per-channel map takes the formal default.
    ## a zero multiple would start the drift at the response onset, where
    ## the drift line absorbs the primary response
    per_channel <- lapply(setup$per_channel, \(.a) {
        tau_mult <- .a$tau_mult %||% 3
        validate_numeric(
            tau_mult, 1, c(0, Inf),
            inclusive = "right", msg1 = "one-element positive", env = env
        )
        .a$tau_mult <- tau_mult
        .a
    })

    time_channel <- setup$time_channel
    ## NA scaffold (method columns only) for convergence failure
    # fmt: skip
    na_cols <- c(
        "A", "B", "tau", "k", "TD", "MRT", "HRT", "MRT_fitted", "HRT_fitted",
        "slope", "tau_mult", "texc", "texc_fitted"
    )

    ## method-specific fit: self-starting exponential-drift via nls; a
    ## failed 6-param fit falls back to the 5-param model
    expdrift_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## the drift onset multiple is always held constant
        .a$fix <- c(.a$fix, list(tau_mult = .a$tau_mult))

        fit <- fit_td_fallback(
            x_fit,
            t_fit,
            # fmt: skip
            params = c(
                "A", "B", "tau", "slope", "tau_mult", if (.a$use_TD) "TD"
            ),
            .a,
            fitter = \(.data, .params, on_error) {
                ## tau and TD are held non-negative; the hinge is non-smooth,
                ## so port often stops short of its certificate on usable
                ## coefficients, which are kept with a warning
                free <- setdiff(.params, names(.a$fix))
                lower <- c(
                    tau = diff(range(.data[[2L]])) * 1e-6,
                    TD = 0
                )[free]
                lower[is.na(lower)] <- -Inf
                formula <- build_ss_formula(
                    quote(SSexponential_drift),
                    .params,
                    .a$fix,
                    names(.data)[[1L]],
                    names(.data)[[2L]]
                )
                model <- tryCatch(
                    suppressWarnings(nls(
                        formula,
                        .data,
                        algorithm = "port",
                        lower = lower,
                        control = stats::nls.control(warnOnly = TRUE)
                    )),
                    error = on_error
                )
                accept_port_fit(model, on_error)
            },
            fn = quote(SSexponential_drift),
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

        ## enforce direction: bounded refit on D = B - A when inverted
        enforced <- enforce_direction(
            fit$model,
            coefs,
            fit$data,
            direction = .a$direction,
            amp_fn = quote(exponential_drift),
            ## data-scaled tau floor: tau pinned here is a degenerate
            ## step fit, not a genuine response
            lower = if (!"tau" %in% names(.a$fix)) {
                c(tau = diff(range(t_fit)) * 1e-6)
            },
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
        texc_val <- sum(TD_arg, coefs[["tau_mult"]] * coefs[["tau"]])

        ## predict response at MRT, HRT, and texc using the full fitted model
        fitted_params <- exponential_drift(
            t = c(MRT_val, HRT_val, texc_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau = coefs[["tau"]],
            slope = coefs[["slope"]],
            tau_mult = coefs[["tau_mult"]],
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
                HRT_fitted = fitted_params[[2L]],
                slope = coefs[["slope"]],
                tau_mult = coefs[["tau_mult"]],
                texc = texc_val,
                texc_fitted = fitted_params[[3L]]
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
        expdrift_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
