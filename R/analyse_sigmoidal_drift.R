#' Sigmoidal-drift function
#'
#' Calculate a two-phase curve: a primary sigmoidal response with a
#' secondary linear drift beginning near the ending asymptote.
#'
#' @param slope_B A numeric parameter for the linear drift rate `dx/dt`
#'   of the secondary phase at the ending asymptote `B`, in response units
#'   per unit of the predictor variable `t`.
#' @param drift_fraction A numeric fraction of the amplitude `B - A` in
#'   `(0.5, 1)` at which the linear drift begins: the drift onset is where
#'   the sigmoid reaches `A + drift_fraction * (B - A)` (see
#'   [sigdrift_onset()]).
#' @param shape Character; the sigmoidal shape. One of `"symmetric"`
#'   (*default*; [logistic()]), `"gompertz"` ([gompertz()]), or
#'   `"gompertz_left"` ([gompertz_left()]).
#' @inheritParams logistic
#'
#' @details
#' Model:
#' `S(t) + slope_B * pmax(t - onset, 0)`
#'
#' where `S(t)` is the 4-parameter sigmoid of the given `shape` with
#' asymptotes `A` and `B`, inflection `xmid`, and inflection rate `slope`.
#' The drift is a hinge line anchored at zero at the onset, so it is
#' exactly zero up to the onset.
#'
#' The onset is the analytic inverse of each shape at the `drift_fraction`
#' fraction `f` of its amplitude, `onset = xmid + u / k`:
#'
#' - `shape = "symmetric"`: `k = 4 * slope / (B - A)`;
#'   `u = log(f / (1 - f))`.
#' - `shape = "gompertz"`: `k = slope * e / (B - A)`;
#'   `u = -log(-log(f))`.
#' - `shape = "gompertz_left"`: `k = slope * e / (B - A)`;
#'   `u = log(-log(1 - f))`.
#'
#' The `"gompertz"` form places its onset furthest past `xmid` (slow
#' tail) and `"gompertz_left"` nearest (fast tail).
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSsigmoidal_drift()], [logistic()],
#'   [gompertz()], [exponential_drift()]
#'
#' @examples
#' ## create a sigmoidal curve with a late linear drift and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_B = -0.4, drift_fraction = 0.95
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## the drift onset fraction is held constant in the formula
#' model <- nls(
#'     x ~ SSsigmoidal_drift(t, A, B, xmid, slope, slope_B, drift_fraction = 0.95),
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
    slope_B,
    drift_fraction,
    shape = c("symmetric", "gompertz", "gompertz_left")
) {
    shape <- match.arg(shape)
    ## primary sigmoid + hinge-linear secondary drift from the onset
    S <- switch(shape, symmetric = logistic, gompertz = gompertz, gompertz_left)
    onset <- sigdrift_onset(A, B, xmid, slope, drift_fraction, shape)
    return(S(t, A, B, xmid, slope) + slope_B * pmax(t - onset, 0))
}


#' Rate constant of a sigmoidal shape
#'
#' The rate `k` such that the sigmoid of the given `shape` is a function of
#' `u = k * (t - xmid)`: `4 * slope / (B - A)` for `"symmetric"`, else
#' `slope * e / (B - A)`. Positive for a consistent fit, where `slope` and
#' `B - A` share a sign.
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A numeric rate in units of `1 / t`.
#'
#' @keywords internal
sigdrift_rate <- function(A, B, slope, shape) {
    return(slope * (if (shape == "symmetric") 4 else exp(1)) / (B - A))
}


#' Drift onset time of the sigmoidal-drift model
#'
#' The time at which a sigmoid of the given `shape` reaches the
#' `drift_fraction` fraction of its amplitude, by the analytic inverse of
#' each shape (see [sigmoidal_drift()]). Vectorised over the numeric
#' parameters; `shape` is a single string.
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A numeric vector of onset times.
#'
#' @keywords internal
sigdrift_onset <- function(A, B, xmid, slope, drift_fraction, shape) {
    ## a fraction outside (0, 1) has no onset
    if (any(drift_fraction <= 0 | drift_fraction >= 1, na.rm = TRUE)) {
        stop("`drift_fraction` must be a fraction of the amplitude in (0, 1).")
    }
    f <- drift_fraction
    u <- switch(
        shape,
        symmetric = log(f / (1 - f)),
        gompertz = -log(-log(f)),
        gompertz_left = log(-log1p(-f))
    )
    return(xmid + u / sigdrift_rate(A, B, slope, shape))
}


#' Excursion point of the sigmoidal-drift model
#'
#' The time past the inflection at which the drift rate overtakes the
#' decaying sigmoid rate, `|S'(t)| = |slope_B|`, floored at the drift onset
#' (see [sigdrift_onset()]): the turning point of the curve when the
#' phases oppose, or where the linear trend takes over a monotonic
#' response. A drift at least as fast as the peak sigmoid rate `slope`
#' takes over from the onset. Scalar parameters only.
#'
#' @inheritParams sigmoidal_drift
#'
#' @details
#' With `r = |slope_B / slope|` and `u = k * (t - xmid)` (see
#' [sigdrift_rate()]), the sigmoid rate relative to its peak is
#' `4 * L * (1 - L)` with `L = 1 / (1 + exp(-u))` for `"symmetric"`, solved
#' as `u = 2 * atanh(sqrt(1 - r))`; `exp(1 - u - exp(-u))` for
#' `"gompertz"`; and `exp(1 + u - exp(u))` for `"gompertz_left"`. The
#' Gompertz forms have no closed inverse and are solved by
#' [stats::uniroot()] on a bracket containing the single post-inflection
#' root.
#'
#' @returns A numeric excursion time.
#'
#' @keywords internal
sigdrift_texc <- function(A, B, xmid, slope, slope_B, drift_fraction, shape) {
    onset <- sigdrift_onset(A, B, xmid, slope, drift_fraction, shape)
    r <- abs(slope_B / slope)
    if (!is.finite(r) || r >= 1) {
        return(onset)
    }
    u <- switch(
        shape,
        symmetric = 2 * atanh(sqrt(1 - r)),
        gompertz = stats::uniroot(
            \(u) u + exp(-u) - 1 + log(r),
            c(0, 1 - log(r))
        )$root,
        gompertz_left = stats::uniroot(
            \(u) exp(u) - u - 1 + log(r),
            c(0, log(2 * (1 - log(r))))
        )$root
    )
    return(max(onset, xmid + u / sigdrift_rate(A, B, slope, shape)))
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
#' [SSgompertz()] ([init_asymptotes()], [init_inflection()]), the drift
#' onset resolved from that seed, and the residual from the seeded sigmoid
#' past the onset regressed on time from the onset: the intercept corrects
#' the asymptote `B` and the slope is the drift. A second pass re-seeds the
#' sigmoid on the drift-corrected response, correcting an inflection
#' biased by the drift. Fewer than two points past the onset seed a zero
#' drift. User-fixed values are held.
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
    ## `[[` throughout: `$` would partial-match `slope` to `slope_B`
    p <- fixed[["drift_fraction"]] %||% 0.95
    ab <- init_asymptotes(x)
    A <- fixed[["A"]] %||% ab$A

    refine <- \(s) {
        ## inflection from the drift-corrected response
        xd <- x - s$slope_B * pmax(t - s$onset, 0)
        infl <- init_inflection(xd, t, A, s$B)
        xmid <- fixed[["xmid"]] %||% infl$xmid
        slope <- fixed[["slope"]] %||% infl$slope
        onset <- sigdrift_onset(A, s$B, xmid, slope, p, shape)

        ## the residual from the seeded sigmoid past the onset is the
        ## asymptote error plus the drift from the onset
        i <- t >= onset
        r <- x[i] - sigmoidal_drift(t[i], A, s$B, xmid, slope, 0, p, shape)
        b <- slope(r, t[i] - onset, intercept = TRUE, bypass_checks = TRUE)
        list(
            B = fixed[["B"]] %||% (s$B + (attr(b, "intercept") %||% 0)),
            xmid = xmid,
            slope = slope,
            slope_B = fixed[["slope_B"]] %||% (if (is.na(b)) 0 else b),
            onset = onset
        )
    }

    ## an onset at the record end makes the initial drift correction zero
    s <- refine(refine(list(
        B = fixed[["B"]] %||% ab$B,
        slope_B = fixed[["slope_B"]] %||% 0,
        onset = max(t)
    )))
    return(c(
        A = A,
        B = s$B,
        xmid = s$xmid,
        slope = s$slope,
        slope_B = s$slope_B,
        drift_fraction = p
    ))
}


#' Self-starting sigmoidal-drift model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [sigmoidal_drift()], for use with [stats::nls()]: a 4-parameter sigmoid
#' (`A`, `B`, `xmid`, `slope`) with a linear drift `slope_B` at its ending
#' asymptote from the onset fraction `drift_fraction`.
#'
#' @usage
#' SSsigmoidal_drift(t, A, B, xmid, slope, slope_B, drift_fraction,
#'     shape = "symmetric")
#'
#' @inheritParams sigmoidal_drift
#'
#' @details
#' `x ~ SSsigmoidal_drift(t, A, B, xmid, slope, slope_B,
#' drift_fraction = 0.95, shape = "gompertz")`
#'
#' `drift_fraction` should be written as a constant, and `shape` is a string
#' constant (`"symmetric"` when omitted); neither is estimated. The hinge
#' at the onset is not differentiable, so `algorithm = "port"` with
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g.
#'   `x ~ SSsigmoidal_drift(t, A = 0, B, xmid, slope, slope_B,
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
#' ## create a Gompertz curve with a late linear drift and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_B = -0.4, drift_fraction = 0.95, shape = "gompertz"
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' model <- nls(
#'     x ~ SSsigmoidal_drift(
#'         t, A, B, xmid, slope, slope_B,
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
        c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    ),
    parameters = c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
)


#' Analyse sigmoidal-drift kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "sigmoidal_drift")`. Fits a two-phase
#' sigmoidal + linear-drift curve to each `nirs_channel` within a single
#' *"mnirs"* data frame. See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param drift_fraction A numeric fraction of the amplitude in `(0.5, 1)` at
#'   which the drift onset is held (*default* `0.95`). Always held constant.
#'   Applied to every channel, or per-channel as a list keyed by channel
#'   name, e.g. `drift_fraction = list(smo2 = 0.9)`.
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `xmid`, `slope`, `slope_B`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Applied to every channel, or per-channel as a
#'   list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`.
#' @inheritParams analyse_logistic
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams analyse_monoexponential
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `xmid`, `slope`, `texc`, `slope_B`,
#'   `drift_fraction`, `xmid_fitted`, `texc_fitted`. `texc` is the
#'   excursion point where the drift rate overtakes the decaying sigmoid
#'   rate, never before the drift onset (see [sigdrift_texc()]).
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
        fix_params = c("A", "B", "xmid", "slope", "slope_B"),
        verbose = verbose,
        env = env
    )
    per_channel <- resolve_drift_frac(setup$per_channel, env)

    time_channel <- setup$time_channel
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- kinetics_coef_cols$sigmoidal_drift
    params <- c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    fn <- quote(SSsigmoidal_drift)

    ## method-specific fit: self-starting sigmoidal-drift via nls
    sigdrift_fit <- function(.nirs, x_fit, t_fit, .a, valid) {
        ## the drift onset fraction is always held constant; the shape
        ## rides in the formula as a string constant beside the fixed
        ## parameters
        .a$fix <- c(.a$fix, list(drift_fraction = .a$drift_fraction))
        fix_all <- c(.a$fix, list(shape = .a$shape))
        free <- setdiff(params, names(.a$fix))
        ## columns carry the channel names so the model predicts on them
        nm <- fit_names(.nirs, time_channel, params)
        fit_data <- setNames(data.frame(x_fit, t_fit), nm)
        on_error <- \(e) warn_fit_failed(fn, e, .nirs, interval_name, env = env)

        ## the hinge is non-smooth, so port often stops short of its
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

        ## xmid and texc are elapsed from start_time, matching the fit
        ## time base; the coefficients are in model-argument order
        cf <- c(as.list(coefs), shape = .a$shape)
        texc_val <- do.call(sigdrift_texc, cf)
        ## predict response at xmid and texc using the full fitted model
        fitted <- do.call(sigmoidal_drift, c(list(c(cf$xmid, texc_val)), cf))

        build_fit_results(
            data.frame(
                t(coefs),
                texc = texc_val,
                xmid_fitted = fitted[[1L]],
                texc_fitted = fitted[[2L]]
            )[na_cols],
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
