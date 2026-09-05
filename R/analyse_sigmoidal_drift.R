#' Sigmoidal-drift function
#'
#' Calculate a three-phase curve: a primary sigmoidal response with
#' independent linear drifts at the leading and trailing asymptotes.
#'
#' @param slope_A,slope_B Numeric parameters for the linear drift rates
#'   `dx/dt` at the starting asymptote `A` and the ending asymptote `B`, in
#'   response units per unit of the predictor variable `t`.
#' @param drift_fraction A numeric fraction of the amplitude `B - A` in
#'   `(0.5, 1)` bounding the drift regions: the leading drift applies only
#'   before the sigmoid reaches `A + (1 - drift_fraction) * (B - A)` (at
#'   `texc_A`), the trailing drift only after it reaches
#'   `A + drift_fraction * (B - A)` (at `texc_B`).
#' @param shape Character; the sigmoidal shape. One of `"symmetric"`
#'   (*default*; [logistic()]), `"gompertz"` ([gompertz()]), or
#'   `"gompertz_left"` ([gompertz_left()]).
#' @inheritParams logistic
#'
#' @details
#' Model:
#' `S(t) + slope_A * pmin(t - texc_A, 0) + slope_B * pmax(t - texc_B, 0)`
#'
#' where `S(t)` is the 4-parameter sigmoid of the given `shape` with
#' asymptotes `A` and `B`, inflection `xmid`, and inflection rate `slope`.
#' Each drift is a hinge line anchored at zero at its cutoff, so the
#' leading drift is exactly zero from `texc_A` on and the trailing drift
#' exactly zero up to `texc_B`. The cutoffs never overlap
#' (`texc_A < xmid < texc_B`), so the two drifts are independent.
#'
#' The cutoffs are where the sigmoid reaches the `1 - drift_fraction` and
#' `drift_fraction` fractions of its amplitude (`p = 1 - drift_fraction`):
#'
#' - `shape = "symmetric"`: `k = 4 * slope / (B - A)`;
#'   `texc_A = xmid - log((1 - p) / p) / k`,
#'   `texc_B = xmid + log((1 - p) / p) / k`.
#' - `shape = "gompertz"`: `k = slope * e / (B - A)`;
#'   `texc_A = xmid - log(-log(p)) / k`,
#'   `texc_B = xmid - log(-log(1 - p)) / k`.
#' - `shape = "gompertz_left"`: `k = slope * e / (B - A)`;
#'   `texc_A = xmid + log(-log(1 - p)) / k`,
#'   `texc_B = xmid + log(-log(p)) / k`.
#'
#' The Gompertz forms give a longer cutoff distance on their slow side.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSsigmoidal_drift()], [logistic()],
#'   [gompertz()], [exponential_drift()]
#'
#' @examples
#' ## create a sigmoidal curve with drifting asymptotes and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_A = 0.3, slope_B = -0.4, drift_fraction = 0.95
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## the cutoff fraction is held constant in the formula
#' model <- nls(
#'     x ~ SSsigmoidal_drift(
#'         t, A, B, xmid, slope, slope_A, slope_B, drift_fraction = 0.95
#'     ),
#'     data = data,
#'     algorithm = "port",
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
sigmoidal_drift <- function(
    t,
    A,
    B,
    xmid,
    slope,
    slope_A,
    slope_B,
    drift_fraction,
    shape = c("symmetric", "gompertz", "gompertz_left")
) {
    shape <- match.arg(shape)
    cut <- sigdrift_cutoffs(A, B, xmid, slope, drift_fraction, shape)
    ## primary sigmoid + hinge-linear drifts outside the cutoffs
    S <- switch(
        shape,
        symmetric = logistic(t, A, B, xmid, slope),
        gompertz = gompertz(t, A, B, xmid, slope),
        gompertz_left = gompertz_left(t, A, B, xmid, slope)
    )
    return(
        S + slope_A * pmin(t - cut[[1L]], 0) + slope_B * pmax(t - cut[[2L]], 0)
    )
}


#' Drift cutoff times of the sigmoidal-drift model
#'
#' The times at which a sigmoid of the given `shape` reaches the
#' `1 - drift_fraction` and `drift_fraction` fractions of its amplitude, by the
#' analytic inverse of each shape (see [sigmoidal_drift()]).
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A named numeric vector `c(texc_A, texc_B)`.
#'
#' @keywords internal
sigdrift_cutoffs <- function(A, B, xmid, slope, drift_fraction, shape) {
    p <- 1 - drift_fraction
    cut <- switch(
        shape,
        symmetric = {
            q <- log((1 - p) / p) * (B - A) / (4 * slope)
            c(xmid - q, xmid + q)
        },
        gompertz = {
            k <- slope * exp(1) / (B - A)
            xmid - log(-log(c(p, 1 - p))) / k
        },
        gompertz_left = {
            k <- slope * exp(1) / (B - A)
            xmid + log(-log(c(1 - p, p))) / k
        }
    )
    return(setNames(cut, c("texc_A", "texc_B")))
}


#' Initiate self-starting sigmoidal-drift model
#'
#' [sigdrift_init()]: Returns initial values for the parameters in a
#' `selfStart` model. The `shape` written in the model call seeds the
#' matching sigmoid (`"symmetric"` when absent).
#'
#' @inheritParams logistic_init
#'
#' @returns [sigdrift_init()]: Initial starting estimates for parameters in
#'   the model called by [SSsigmoidal_drift()].
#'
#' @keywords internal
sigdrift_init <- function(mCall, data, LHS, ...) {
    fixed <- list(...)$fixed %||% list()
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    shape <- eval(mCall[["shape"]], data) %||% "symmetric"
    return(sigdrift_start(tx[["y"]], tx[["x"]], fixed, shape))
}


#' Starting estimates for the sigmoidal-drift model
#'
#' Vector-level initialiser behind [sigdrift_init()], called directly by
#' the kinetics worker on the fit window. The sigmoid is seeded as for
#' [SSgompertz()] ([init_asymptotes()], [init_inflection()]), the cutoffs
#' resolved from that seed, and each tail's residual from the seeded model
#' regressed on time from its cutoff: the intercept corrects the asymptote
#' and the slope the drift. A second pass re-seeds the sigmoid on the
#' drift-corrected response, correcting an inflection biased by the drift.
#' A tail with fewer than two points seeds a zero drift. User-fixed values
#' are held.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param t A numeric vector of the predictor variable.
#' @param fixed A named list of user-fixed parameter values.
#' @inheritParams sigmoidal_drift
#'
#' @returns A named numeric vector of starting estimates in model order.
#'
#' @keywords internal
sigdrift_start <- function(x, t, fixed = list(), shape = "symmetric") {
    ## `[[` throughout: `$` would partial-match `slope` to `slope_A`
    p <- fixed[["drift_fraction"]] %||% 0.95
    ab <- init_asymptotes(x)

    ## a tail's residual regressed on time from its cutoff: the intercept
    ## and slope correct the asymptote and drift, else no correction
    tail <- \(i, .cut, .r) {
        b <- slope(
            .r[i],
            t[i] - .cut,
            intercept = TRUE,
            bypass_checks = TRUE,
            min_obs = 2L
        )
        if (is.na(b)) c(0, 0) else c(attr(b, "intercept"), b)
    }

    refine <- \(s) {
        ## inflection from the drift-corrected response
        xd <- x -
            s$slope_A * pmin(t - s$texc_A, 0) -
            s$slope_B * pmax(t - s$texc_B, 0)
        infl <- init_inflection(xd, t, s$A, s$B)
        xmid <- fixed[["xmid"]] %||% infl$xmid
        slope <- fixed[["slope"]] %||% infl$slope
        cut <- sigdrift_cutoffs(s$A, s$B, xmid, slope, p, shape)

        ## the residual from the seeded model is, on each tail, the
        ## asymptote error plus the drift error from the cutoff
        r <- x -
            sigmoidal_drift(
                t,
                s$A,
                s$B,
                xmid,
                slope,
                s$slope_A,
                s$slope_B,
                p,
                shape
            )
        lead <- tail(t <= cut[[1L]], cut[[1L]], r)
        trail <- tail(t >= cut[[2L]], cut[[2L]], r)
        list(
            A = fixed[["A"]] %||% (s$A + lead[[1L]]),
            B = fixed[["B"]] %||% (s$B + trail[[1L]]),
            xmid = xmid,
            slope = slope,
            slope_A = fixed[["slope_A"]] %||% (s$slope_A + lead[[2L]]),
            slope_B = fixed[["slope_B"]] %||% (s$slope_B + trail[[2L]]),
            texc_A = cut[[1L]],
            texc_B = cut[[2L]]
        )
    }

    ## cutoffs at the record ends make the initial drift correction zero
    s <- refine(refine(list(
        A = fixed[["A"]] %||% ab$A,
        B = fixed[["B"]] %||% ab$B,
        slope_A = fixed[["slope_A"]] %||% 0,
        slope_B = fixed[["slope_B"]] %||% 0,
        texc_A = min(t),
        texc_B = max(t)
    )))
    return(c(
        A = s$A,
        B = s$B,
        xmid = s$xmid,
        slope = s$slope,
        slope_A = s$slope_A,
        slope_B = s$slope_B,
        drift_fraction = p
    ))
}


#' Self-starting sigmoidal-drift model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [sigmoidal_drift()], for use with [stats::nls()]: a 4-parameter sigmoid
#' (`A`, `B`, `xmid`, `slope`) with linear drifts `slope_A` and `slope_B`
#' at its asymptotes, outside the cutoff fraction `drift_fraction`.
#'
#' @usage
#' SSsigmoidal_drift(t, A, B, xmid, slope, slope_A, slope_B, drift_fraction,
#'     shape = "symmetric")
#'
#' @inheritParams sigmoidal_drift
#'
#' @details
#' `x ~ SSsigmoidal_drift(t, A, B, xmid, slope, slope_A, slope_B,
#' drift_fraction = 0.95, shape = "gompertz")`
#'
#' `drift_fraction` should be written as a constant, and `shape` is a string
#' constant (`"symmetric"` when omitted); neither is estimated. The hinges
#' at the cutoffs are not differentiable, so `algorithm = "port"` with
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSsigmoidal_drift(t, A = 0, B, xmid, slope, slope_A, slope_B,
#'   drift_fraction = 0.95)` holds the starting asymptote at `0`. Fixed
#'   parameters are excluded from estimation and are not returned by
#'   [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [sigmoidal_drift()], [stats::nls()], [stats::selfStart()],
#'   [SSlogistic()], [SSgompertz()], [SSexponential_drift()]
#'
#' @examples
#' ## create a Gompertz curve with drifting asymptotes and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_A = 0.3, slope_B = -0.4, drift_fraction = 0.95, shape = "gompertz"
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSsigmoidal_drift(
#'         t, A, B, xmid, slope, slope_A, slope_B,
#'         drift_fraction = 0.95, shape = "gompertz"
#'     ),
#'     data = data,
#'     algorithm = "port",
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' @export
SSsigmoidal_drift <- selfStart(
    model = sigmoidal_drift,
    initial = init_fixed(
        sigdrift_init,
        c("A", "B", "xmid", "slope", "slope_A", "slope_B", "drift_fraction")
    ),
    parameters = c(
        "A",
        "B",
        "xmid",
        "slope",
        "slope_A",
        "slope_B",
        "drift_fraction"
    )
)


#' Analyse sigmoidal-drift kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "sigmoidal_drift")`. Fits a sigmoidal curve
#' with linear drifts at both asymptotes to each `nirs_channel` within a
#' single *"mnirs"* data frame. See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param drift_fraction A numeric fraction of the amplitude in `(0.5, 1)`
#'   bounding the drift regions (*default* `0.95`): the leading drift
#'   applies below `A + (1 - drift_fraction) * (B - A)` and the trailing drift
#'   above `A + drift_fraction * (B - A)`. Always held constant. Applied to
#'   every channel, or per-channel as a list keyed by channel name, e.g.
#'   `drift_fraction = list(smo2 = 0.9)`.
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `xmid`, `slope`, `slope_A`, `slope_B`) to hold constant during fitting,
#'   e.g. `fix = list(A = 0)`. Applied to every channel, or per-channel as
#'   a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`.
#' @inheritParams analyse_logistic
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams analyse_monoexponential
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `xmid`, `slope`, `slope_A`, `slope_B`,
#'   `drift_fraction`, `texc_A`, `texc_B`, `xmid_fitted`. `texc_A` and `texc_B`
#'   are the drift cutoff times (see [sigmoidal_drift()]).
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
#' @seealso [analyse_kinetics()], [sigmoidal_drift()],
#'   [SSsigmoidal_drift()]
#'
#' @keywords internal
analyse_sigmoidal_drift <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    shape = c("symmetric", "gompertz", "gompertz_left"),
    drift_fraction = 0.95,
    fix = NULL,
    control = NULL,
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
        arg_list = mget(c(
            "shape", "drift_fraction", "fix", "control", "start_time",
            "direction", "end_window"
        )),
        choices = list(
            shape = c("symmetric", "gompertz", "gompertz_left"),
            direction = c("auto", "positive", "negative")
        ),
        fix_params = c("A", "B", "xmid", "slope", "slope_A", "slope_B"),
        verbose = verbose,
        env = env
    )
    per_channel <- resolve_drift_frac(setup$per_channel, env)

    time_channel <- setup$time_channel
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- kinetics_coef_cols$sigmoidal_drift
    params <- c("A", "B", "xmid", "slope", "slope_A", "slope_B", "drift_fraction")
    fn <- quote(SSsigmoidal_drift)

    ## method-specific fit: self-starting sigmoidal-drift via nls
    sigdrift_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## the cutoff fraction is always held constant; the shape rides in
        ## the formula as a string constant beside the fixed parameters
        .a$fix <- c(.a$fix, list(drift_fraction = .a$drift_fraction))
        fix_all <- c(.a$fix, list(shape = .a$shape))
        free <- setdiff(params, names(.a$fix))
        ## columns carry the channel names so the model predicts on them
        nm <- fit_names(.nirs, time_channel, params)
        fit_data <- setNames(data.frame(x_fit, t_fit), nm)
        on_error <- \(e) warn_fit_failed(fn, e, .nirs, interval_name, env = env)

        ## the hinges are non-smooth, so port often stops short of its
        ## certificate on usable coefficients, which are kept with a warning
        model <- if (nrow(fit_data) <= length(free)) {
            on_error(simpleError(sprintf(
                "%d observation%s for %d free parameters.",
                nrow(fit_data),
                if (nrow(fit_data) == 1L) "" else "s",
                length(free)
            )))
        } else {
            formula <- build_ss_formula(
                fn,
                c(params, "shape"),
                fix_all,
                nm[[1L]],
                nm[[2L]]
            )
            tryCatch(
                {
                    start <- sigdrift_start(x_fit, t_fit, .a$fix, .a$shape)
                    embed_fit_call(suppressWarnings(nls(
                        formula,
                        fit_data,
                        start = start[free],
                        algorithm = "port",
                        control = fit_control(
                            .a$control,
                            maxiter = 500L,
                            warnOnly = TRUE
                        )
                    )))
                },
                error = on_error
            )
        }
        model <- accept_port_fit(model, on_error)
        if (is.null(model)) {
            return(build_na_results(na_cols))
        }

        coefs <- full_coefs(model, params, .a$fix)

        ## enforce direction: bounded refit on D = B - A and slope sign.
        ## data-scaled slope floor: slope pinned here is a degenerate
        ## flat fit, not a genuine response
        want <- if (.a$direction == "positive") 1 else -1
        slope_eps <- diff(range(x_fit)) / diff(range(t_fit)) * 1e-6
        slope_free <- !"slope" %in% names(.a$fix)
        enforced <- enforce_direction(
            model,
            coefs,
            fit_data,
            direction = .a$direction,
            amp_fn = quote(sigmoidal_drift),
            lower = if (slope_free) {
                c(slope = if (want > 0) slope_eps else -Inf)
            },
            upper = if (slope_free) {
                c(slope = if (want > 0) Inf else -slope_eps)
            },
            fix = fix_all,
            control = .a$control,
            .nirs = .nirs,
            interval_name = interval_name,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_cols))
        }
        model <- enforced$model
        coefs <- enforced$coefs

        ## cutoff times and xmid are elapsed from start_time, matching the
        ## fit time base
        cut <- sigdrift_cutoffs(
            coefs[["A"]],
            coefs[["B"]],
            coefs[["xmid"]],
            coefs[["slope"]],
            coefs[["drift_fraction"]],
            .a$shape
        )
        xmid_fitted <- as.numeric(
            stats::predict(
                model,
                setNames(data.frame(coefs[["xmid"]]), nm[[2L]])
            )
        )

        build_fit_results(
            data.frame(
                A = coefs[["A"]],
                B = coefs[["B"]],
                xmid = coefs[["xmid"]],
                slope = coefs[["slope"]],
                slope_A = coefs[["slope_A"]],
                slope_B = coefs[["slope_B"]],
                drift_fraction = coefs[["drift_fraction"]],
                texc_A = cut[[1L]],
                texc_B = cut[[2L]],
                xmid_fitted = xmid_fitted
            ),
            model,
            x_fit,
            t_fit,
            valid,
            env = env
        )
    }

    return(analyse_kinetics_channels(
        data,
        setup$nirs_channels,
        setup$time_channel,
        per_channel,
        sigdrift_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
