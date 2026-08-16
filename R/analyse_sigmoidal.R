#' Generalised logistic function
#'
#' Calculate a 4- or 5-parameter logistic (sigmoidal) curve.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting asymptote of the response
#'   variable.
#' @param B A numeric parameter for the ending asymptote of the response
#'   variable.
#' @param xmid A numeric parameter for the time value at the inflection
#'   (steepest) point of the curve, in units of the predictor variable `t`.
#' @param slope A numeric parameter for the slope `dx/dt` of the response
#'   variable at the inflection `xmid`.
#' @param asym A numeric parameter for the asymmetry index of the curve,
#'   equal to the fraction of the response where the inflection `xmid` occurs,
#'   bounded in `c(0, 1)` for `(y(xmid) - A) / (B - A)`. `asym = 0.5`
#'   is symmetric and equivalent to the 4-parameter form. If `NULL`
#'   (*default*), a symmetric 4-parameter model is used.
#'
#' @details
#' The 4-parameter symmetric form is fit by [analyse_kinetics()] when
#' `method = "sigmoidal"` and `shape = "symmetric"` (*default*) via the
#' self-starting wrapper [SSlogistic()].
#'
#' The 5-parameter Richards form is exported for advanced use directly with
#' [stats::nls()] but is not used by [analyse_kinetics()] due to convergence
#' instability. For asymmetric responses, prefer [gompertz()] / 
#' [gompertz_left()], which are more stable.
#'
#' ## Model equations
#'
#' Logistic models are re-parameterised from a Richards generalised logistic
#'   model to be interpretable.
#'
#' 4-parameter (symmetric) model:
#'   `A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))`
#'
#' 5-parameter (asymmetric) model re-parameterised so `asym` is the
#'   inflection fraction. Internally:
#'
#'   * `y = A + (B - A) / (1 + exp(-k * (t - xmid)))^(1 / v)`
#'   * `k = 2 * slope * v / ((B - A) * asym)`
#'   * `v = -log(2) / log(asym)`
#'
#' Inflection is at `t = xmid` with `dx/dt = slope` and
#'   `y(xmid) = A + (B - A) * asym` for any `asym` in `(0, 1)`. 
#' 
#'   * `asym = 0.5 -> v = 1` the model collapses to the 4-parameter form.
#'   * `asym -> 0` gives an early-acceleration curve (inflection near `A`),
#'   * `asym -> 1` gives a late-acceleration curve (inflection near `B`).
#'   * `asym = 0.368` (`1/e`) approximates a right-inflection Gompertz curve.
#'   * `asym = 0.632` (`1 - 1/e`) approximates a left-inflection Gompertz curve.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSlogistic()], [monoexponential()]
#'
#' @examples
#' ## create a logistic curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4, asym = 0.3) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(x, t)
#'
#' model <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
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
logistic <- function(t, A, B, xmid, slope, asym = NULL) {
    if (is.null(asym)) {
        ## 4-parameter symmetric
        y <- A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))
    } else {
        ## 5-parameter Richards re-parameterised: asym is the
        ## inflection-height fraction (y(xmid) - A) / (B - A).
        ## clamp to keep log(asym) finite during nls iteration
        asym <- min(max(asym, 1e-6), 1 - 1e-6)
        v <- -log(2) / log(asym)
        k <- 2 * slope * v / ((B - A) * asym)
        y <- A + (B - A) / (1 + exp(-k * (t - xmid)))^(1 / v)
    }
    return(y)
}


#' Gompertz growth function
#'
#' Calculate 4-parameter Gompertz (asymmetric sigmoidal) curves.
#'
#' @inheritParams logistic
#'
#' @details
#' The [gompertz()] curve is asymmetric, with the inflection point `xmid` closer
#' to the starting asymptote `A`; i.e. the response accelerates faster away
#' from `A`, and more slowly approaches the ending asymptote `B`.
#'
#' The modified [gompertz_left()] curve has the inflection point closer to the
#' ending asymptote `B` (slow departure from `A`, late acceleration toward `B`).
#'
#' ## Implementation
#'
#' These models are fit by [analyse_kinetics()] when `method = "sigmoidal"` and
#' `shape = "gompertz"` or `"gompertz_left"` respectively, using [stats::nls()]
#' via the self-starting wrappers [SSgompertz()] and [SSgompertz_left()].
#'
#' ## Model equations
#'
#' Both forms are re-parameterised so that `xmid` is the time at inflection
#'   and `slope` is the response rate `dx/dt` at the inflection point.
#'
#' Gompertz (right-Gompertz; early acceleration, inflection near `A`):
#'
#'   * `k = slope * e / (B - A)`
#'   * `y = A + (B - A) * exp(-exp(-k * (t - xmid)))`
#'
#' Left-Gompertz (late acceleration, inflection near `B`):
#'
#'   * `k = slope * e / (B - A)`
#'   * `y = A + (B - A) * (1 - exp(-exp(k * (t - xmid))))`
#'
#' For both forms, `y(xmid) = A + (B - A) / e` (right) or
#'   `y(xmid) = A + (B - A) * (1 - 1/e)` (left), corresponding to
#'   inflection-height fractions of `1/e` and `1 - 1/e`, respectively.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSgompertz()], [SSgompertz_left()], 
#'   [logistic()]
#'
#' @examples
#' ## create a Gompertz curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)
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
gompertz <- function(t, A, B, xmid, slope) {
    k <- slope * exp(1) / (B - A)
    y <- A + (B - A) * exp(-exp(-k * (t - xmid)))
    return(y)
}


#' @rdname gompertz
#' @export
gompertz_left <- function(t, A, B, xmid, slope) {
    k <- slope * exp(1) / (B - A)
    y <- A + (B - A) * (1 - exp(-exp(k * (t - xmid))))
    return(y)
}


#' Initiate self-starting logistic model
#'
#' [logistic_init()]: Returns initial values for the parameters in a
#' `selfStart` model.
#'
#' @param mCall A matched call to the function `model`.
#' @param data A data frame with predictor `t` and the response variable.
#' @param LHS The left-hand side expression of the model formula.
#' @param ... Additional arguments, including `fixed`, a named list of
#'   user-fixed parameter values from [init_fixed()] used to seed the
#'   remaining free estimates.
#'
#' @returns [logistic_init()]: Initial starting estimates for parameters in
#'   the model called by [SSlogistic()].
#'
#' @keywords internal
logistic_init <- function(mCall, data, LHS, ...) {
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)
    has_asym <- "asym" %in% names(mCall)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## asymptotes from first and last ceiling(n/5) values
    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    B_init <- fixed$B %||% ab$B

    ## linearisation for 4-param: log((B - y) / (y - A)) ~ t
    lo <- min(A_init, B_init)
    hi <- max(A_init, B_init)
    eps <- (hi - lo) * 1e-3
    x_clip <- pmin(pmax(x, lo + eps), hi - eps)
    xf <- log((B_init - x_clip) / (x_clip - A_init))

    xmid_init <- NA_real_
    slope_init <- NA_real_
    finite_idx <- is.finite(xf)
    if (sum(finite_idx) >= 3L) {
        b <- slope(
            xf[finite_idx],
            t[finite_idx],
            intercept = TRUE,
            bypass_checks = TRUE,
            min_obs = 2L
        )
        a <- attr(b, "intercept")
        if (is.finite(b) && b != 0) {
            xmid_init <- -a / b
            slope_init <- -b * (B_init - A_init) / 4
        }
    }

    ## fallbacks for degenerate data
    t_range <- diff(range(t))
    if (!is.finite(xmid_init) || xmid_init < min(t) || xmid_init > max(t)) {
        xmid_init <- t[which.min(abs(x - (A_init + B_init) / 2))]
    }
    if (!is.finite(slope_init) || slope_init == 0) {
        slope_init <- if (t_range > 0) {
            (B_init - A_init) / t_range
        } else {
            sign(B_init - A_init)
        }
    }

    if (!has_asym) {
        return(c(
            A = A_init,
            B = B_init,
            xmid = fixed$xmid %||% xmid_init,
            slope = fixed$slope %||% slope_init
        ))
    }

    ## 5-param: empirical inflection from smoothed derivative
    infl <- init_inflection(x, t, A_init, B_init)
    asym_emp <- (x[infl$idx] - A_init) / (B_init - A_init)
    asym_init <- min(max(asym_emp, 0.1), 0.9)

    return(c(
        A = A_init,
        B = B_init,
        xmid = fixed$xmid %||% infl$xmid,
        slope = fixed$slope %||% infl$slope,
        asym = fixed$asym %||% asym_init
    ))
}


#' Initiate self-starting Gompertz model
#'
#' [gompertz_init()]: Returns initial values for the parameters in a
#' `selfStart` model. Used by both [SSgompertz()] and [SSgompertz_left()];
#' the symmetric logistic linearisation does not apply to Gompertz forms, so
#' initialisation is derivative-based via [init_inflection()].
#'
#' @inheritParams logistic_init
#'
#' @returns [gompertz_init()]: Initial starting estimates for parameters
#'   in the model called by [SSgompertz()] or [SSgompertz_left()].
#'
#' @keywords internal
gompertz_init <- function(mCall, data, LHS, ...) {
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    B_init <- fixed$B %||% ab$B
    infl <- init_inflection(x, t, A_init, B_init)

    return(c(
        A = A_init,
        B = B_init,
        xmid = fixed$xmid %||% infl$xmid,
        slope = fixed$slope %||% infl$slope
    ))
}


#' Estimate baseline and asymptote from the first/last quintile of `x`
#'
#' Shared helper used by self-start initialisers for logistic / Gompertz
#' model families.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param n An integer length of `x`.
#'
#' @returns A list with elements `A` (starting asymptote estimate) and `B`
#'   (ending asymptote estimate).
#'
#' @keywords internal
init_asymptotes <- function(x, n = length(x)) {
    n_asymp <- max(1L, ceiling(n / 5))
    A_init <- mean(x[seq_len(n_asymp)])
    B_init <- mean(x[seq(n - n_asymp + 1L, n)])
    return(list(A = A_init, B = B_init))
}


#' Estimate inflection point from a smoothed first derivative
#'
#' Shared helper that locates the empirical inflection (peak of
#' `|dx/dt|` after smoothing) and returns the corresponding `xmid` and
#' `slope` initial values. Falls back to the half-response point and a
#' mean-rate slope when the derivative is degenerate.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param t A numeric vector of the predictor variable.
#' @param A_init Estimated starting asymptote.
#' @param B_init Estimated ending asymptote.
#'
#' @returns A list with elements `idx` (integer index into `x`), `xmid`
#'   (numeric `t` value at the inflection), and `slope` (numeric `dx/dt`
#'   at the inflection).
#'
#' @keywords internal
init_inflection <- function(x, t, A_init, B_init) {
    dx_dt <- diff(x) / diff(t)
    win <- max(3L, 2L * (length(dx_dt) %/% 20L) + 1L)
    dx_smooth <- as.numeric(stats::filter(dx_dt, rep(1 / win, win), sides = 2L))
    dx_smooth[!is.finite(dx_smooth)] <- 0

    i_infl <- which.max(abs(dx_smooth))
    xmid_init <- t[i_infl]
    slope_init <- dx_smooth[i_infl]

    ## fallback: half-response point with mean-rate slope
    t_range <- diff(range(t))
    if (!is.finite(xmid_init) || xmid_init < min(t) || xmid_init > max(t)) {
        i_infl <- which.min(abs(x - (A_init + B_init) / 2))
        xmid_init <- t[i_infl]
    }
    if (!is.finite(slope_init) || slope_init == 0) {
        slope_init <- if (t_range > 0) {
            (B_init - A_init) / t_range
        } else {
            sign(B_init - A_init)
        }
    }

    return(list(idx = i_infl, xmid = xmid_init, slope = slope_init))
}


#' Self-starting logistic model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [logistic()], for use with [stats::nls()]. Supports both the 4-parameter
#' symmetric form (A, B, xmid, slope) and the 5-parameter asymmetric form
#' (A, B, xmid, slope, asym); arity is inferred from the formula passed to
#' [stats::nls()].
#'
#' @usage
#' SSlogistic(t, A, B, xmid, slope, asym)
#'
#' @inheritParams logistic
#'
#' @details
#' 4-parameter model: `x ~ SSlogistic(t, A, B, xmid, slope)`
#'
#' 5-parameter model: `x ~ SSlogistic(t, A, B, xmid, slope, asym)`
#'
#' The 4-parameter form is used by [analyse_kinetics()] when
#'   `method = "sigmoidal"` and `shape = "symmetric"`. The 5-parameter
#'   asymmetric form is retained for advanced/experimental use only.
#'   [stats::nls()] reads the free parameters from the formula right-hand
#'   side, so omitting `asym` incurs no degrees-of-freedom penalty.
#'
#' `analyse_kinetics()` instead dispatches to [SSgompertz()] /
#'   [SSgompertz_left()] for asymmetric shapes, which are more stable.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSlogistic(t, A = 0, B, xmid, slope)`
#'   fixes the starting asymptote at `A = 0`. Fixed parameters are excluded
#'   from estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [logistic()], [stats::nls()], [stats::selfStart()],
#'   [stats::SSfpl()]
#'
#' @examples
#' ## create a logistic curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4, asym = 0.3) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## 4-parameter fit
#' model4 <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
#' summary(model4)
#'
#' ## 5-parameter fit on the same data
#' model5 <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
#' summary(model5)
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(x ~ SSlogistic(t, A = 10, B, xmid, slope), data = data)
#' summary(model_fixed)
#'
#' y4 <- predict(model4, data)
#' y5 <- predict(model5, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y5, colour = "5-param")) +
#'             ggplot2::geom_line(ggplot2::aes(y = y4, colour = "4-param"))
#'     }
#' }
#'
#' @export
SSlogistic <- selfStart(
    model = logistic,
    initial = init_fixed(logistic_init, c("A", "B", "xmid", "slope", "asym")),
    parameters = c("A", "B", "xmid", "slope", "asym")
)


#' Self-starting Gompertz model
#'
#' Creates initial coefficient estimates for `selfStart` wrappers around
#' [gompertz()] and [gompertz_left()], for use with [stats::nls()]. Both
#' wrappers use the same 4-parameter `(A, B, xmid, slope)` interface.
#'
#' @usage
#' SSgompertz(t, A, B, xmid, slope)
#'
#' SSgompertz_left(t, A, B, xmid, slope)
#'
#' @inheritParams logistic
#'
#' @details
#' Overwrites [stats::SSgompertz()].
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSgompertz(t, A = 0, B, xmid, slope)`
#'   fixes the starting asymptote at `A = 0`. Fixed parameters are excluded
#'   from estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [gompertz()], [gompertz_left()], [SSlogistic()], [stats::nls()],
#'   [stats::selfStart()], [stats::SSgompertz()]
#'
#' @examples
#' ## create a Gompertz curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)
#' summary(model)
#' 
#' ## fix the baseline A at a known value
#' model_fixed <- nls(x ~ SSgompertz(t, A = 10, B, xmid, slope), data = data)
#' summary(model_fixed)
#'
#' ## left-Gompertz
#' set.seed(16)
#' x2 <- gompertz_left(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data2 <- data.frame(t, x = x2)
#' 
#' model_left <- nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = data2)
#' summary(model_left)
#'
#' @export
SSgompertz <- selfStart(
    model = gompertz,
    initial = init_fixed(gompertz_init, c("A", "B", "xmid", "slope")),
    parameters = c("A", "B", "xmid", "slope")
)


#' @rdname SSgompertz
#' @export
SSgompertz_left <- selfStart(
    model = gompertz_left,
    initial = init_fixed(gompertz_init, c("A", "B", "xmid", "slope")),
    parameters = c("A", "B", "xmid", "slope")
)


## shape -> (self-start model fn, amplitude-reparam fn) symbols
shape_dispatch <- list(
    symmetric = list(
        model = quote(SSlogistic),
        amp = quote(logistic)
    ),
    gompertz = list(
        model = quote(SSgompertz),
        amp = quote(gompertz)
    ),
    gompertz_left = list(
        model = quote(SSgompertz_left),
        amp = quote(gompertz_left)
    )
)


#' Analyse logistic kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "sigmoidal")`. Fits a 4-parameter sigmoidal
#' curve to each `nirs_channel` within a single *"mnirs"* data frame with
#' one of three shapes: `"symmetric"`, `"gompertz"`, or `"gompertz_left"`.
#' See [analyse_kinetics()] for user-facing documentation.
#'
#' @param shape Character; the 4-parameter sigmoidal shape to fit. One of
#'   `"symmetric"` (*default*; calls [SSlogistic()]), `"gompertz"`
#'   (early-inflection; calls [SSgompertz()]), or `"gompertz_left"`
#'   (late-inflection; calls [SSgompertz_left()]).
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `xmid`, `slope`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `xmid`, `slope`, `xmid_fitted`.
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
#' @seealso [analyse_kinetics()], [logistic()], [SSlogistic()],
#'   [gompertz()], [gompertz_left()], [SSgompertz()], [SSgompertz_left()]
#'
#' @keywords internal
analyse_logistic <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    shape = c("symmetric", "gompertz", "gompertz_left"),
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
        arg_list = list(
            shape = shape,
            fix = fix,
            start_time = start_time,
            direction = direction,
            end_window = end_window
        ),
        choices = list(
            shape = c("symmetric", "gompertz", "gompertz_left"),
            direction = c("auto", "positive", "negative")
        ),
        verbose = verbose,
        env = env
    )
    nirs_channels <- setup$nirs_channels
    time_channel <- setup$time_channel

    ## validate each channel's fixed parameters before any fitting.
    ## `[` assignment keeps an unfixed channel's `fix` key present but
    ## `NULL`, so every channel serialises the same `channel_args` columns
    per_channel <- lapply(setup$per_channel, \(.a) {
        f <- validate_fix(.a$fix, c("A", "B", "xmid", "slope"), env = env)
        .a["fix"] <- list(if (length(f) > 0L) f else NULL)
        .a
    })

    ## NA scaffold (method columns only) for convergence failure
    na_coefs <- data.frame(
        A = NA_real_,
        B = NA_real_,
        xmid = NA_real_,
        slope = NA_real_,
        xmid_fitted = NA_real_
    )

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, fn, e, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        cli_warn(c(
            "x" = "{.fn {as.character(fn)}} fit failed for \\
            {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        ), call = warn_call(env))
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting sigmoidal via nls
    logistic_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        ## resolve per-channel shape and matching self-start fn
        disp <- shape_dispatch[[.a$shape]]
        ch_fn <- disp$model
        fit_data <- data.frame(.x = x_fit, .t = t_fit)
        params <- c("A", "B", "xmid", "slope")

        ## build nls formula with any fixed params as constants
        nls_formula <- build_ss_formula(ch_fn, params, .a$fix)

        model <- tryCatch(
            nls(nls_formula, fit_data),
            error = \(e) {
                fit_failed_warning(.nirs, ch_fn, e, verbose)
                NULL
            }
        )

        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        coefs <- full_coefs(model, params, .a$fix)

        ## enforce direction: bounded refit on D = B - A and slope sign
        ## data-scaled slope floor: slope pinned here is a degenerate
        ## flat fit, not a genuine response
        want <- if (.a$direction == "positive") 1 else -1
        slope_eps <- diff(range(x_fit)) / diff(range(t_fit)) * 1e-6
        free_extra <- setdiff(params, c("A", "B", names(.a$fix)))
        extra <- coefs[free_extra]
        if ("slope" %in% free_extra) {
            extra[["slope"]] <- want * max(abs(coefs[["slope"]]), slope_eps)
        }
        enforced <- enforce_direction(
            model, coefs, fit_data,
            direction = .a$direction,
            amp_fn = disp$amp,
            extra = extra,
            extra_lower = if ("slope" %in% free_extra) {
                c(slope = if (want > 0) slope_eps else -Inf)
            },
            extra_upper = if ("slope" %in% free_extra) {
                c(slope = if (want > 0) Inf else -slope_eps)
            },
            fn = ch_fn,
            fix = .a$fix,
            .nirs = .nirs,
            interval_name = interval_name,
            verbose = verbose,
            env = env
        )
        if (is.null(enforced)) {
            return(build_na_results(na_coefs))
        }
        model <- enforced$model
        coefs <- enforced$coefs
        fitted_vals <- stats::predict(model)

        ## xmid is already elapsed from start_time, matching the fit time base
        xmid_offset <- coefs[["xmid"]]

        ## predict response at the inflection point xmid
        xmid_fitted <- as.numeric(
            stats::predict(model, data.frame(.t = coefs[["xmid"]]))
        )

        list(
            coefs = data.frame(
                A           = coefs[["A"]],
                B           = coefs[["B"]],
                xmid        = xmid_offset,
                slope       = coefs[["slope"]],
                xmid_fitted = xmid_fitted
            ),
            model = model,
            fitted_data = data.frame(
                window_idx = valid$idx,
                fitted     = fitted_vals
            ),
            diag = compute_diagnostics(
                x_fit,
                t_fit,
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
        logistic_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
