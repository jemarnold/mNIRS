#' Logistic function
#'
#' Calculate a 4- or 5-parameter logistic (sigmoidal) curve. This model family
#' is fit by [analyse_kinetics()] when `method = "logistic"` or `"sigmoidal"`,
#' and by [stats::nls()] via the self-starting wrapper [SSlogistic()].
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting asymptote of the response
#'   variable.
#' @param B A numeric parameter for the ending asymptote of the response
#'   variable.
#' @param xmid A numeric parameter for the `t` value at the inflection
#'   (steepest) point of the curve, in units of the predictor variable `t`.
#' @param slope A numeric parameter for the slope `dx/dt` of the response
#'   variable at the inflection point `xmid`.
#' @param asym A numeric parameter for the asymmetry index of the curve, 
#'   equal to the fraction of the response where the inflection point `xmid` 
#'   occurs, bounded in `c(0, 1)` for `(y(xmid) - A) / (B - A)`. `asym = 0.5` 
#'   is symmetric and equivalent to the 4-parameter form. If `NULL` (*default*), a
#'   symmetric 4-parameter model is used.
#'
#' @details
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
#'   `v = -log(2) / log(asym)`
#'   `k = 2 * slope * v / ((B - A) * asym)`
#'   `y = A + (B - A) / (1 + exp(-k * (t - xmid)))^(1 / v)`
#'
#' Inflection is at `t = xmid` with `dx/dt = slope` and
#'   `y(xmid) = A + (B - A) * asym` for any `asym` in `(0, 1)`. At
#'   `asym = 0.5`, `v = 1` and the model collapses to the 4-parameter form.
#'   `asym -> 0` gives an early-acceleration curve (inflection near `A`), 
#'   `asym -> 1` gives a late-acceleration curve (inflection near `B`).
#' 
#' `asym = 0.368` (`1/e`) approximates a (right-inflection) Gompertz curve.
#'   `asym = 0.632` (`1 - 1/e`) approximates a (left-inflection) modified
#'   Gompertz curve.
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
#' model
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


#' Initiate self-starting logistic model
#'
#' [logistic_init()]: Returns initial values for the parameters in a
#' `selfStart` model.
#'
#' @param mCall A matched call to the function `model`.
#' @param data A data frame with predictor `t` and the response variable.
#' @param LHS The left-hand side expression of the model formula.
#' @param ... Additional arguments.
#'
#' @returns [logistic_init()]: Initial starting estimates for parameters in
#'   the model called by [SSlogistic()].
#'
#' @keywords internal
logistic_init <- function(mCall, data, LHS, ...) {
    ## self-start parameters for nls of logistic fit function;
    ## uses base R `SSfpl()` linearisation approach
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    has_asym <- "asym" %in% names(mCall)

    ## asymptotes from 1/5 response signal; tail order preserves direction
    n_asymp <- max(1, ceiling(n / 5))
    A_init <- mean(x[seq_len(n_asymp)])
    B_init <- mean(x[seq(n - n_asymp + 1, n)])

    ## linearise: y = A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))
    ## => log((B - y) / (y - A)) = -4 * slope / (B - A) * (t - xmid)
    ## clip x strictly inside [min(A, B), max(A, B)] for valid log argument
    lo <- min(A_init, B_init)
    hi <- max(A_init, B_init)
    eps <- (hi - lo) * 1e-3
    x_clip <- pmin(pmax(x, lo + eps), hi - eps)
    xf <- log((B_init - x_clip) / (x_clip - A_init))

    ## init NAs, then overwrite if valid
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
        mid_x <- (A_init + B_init) / 2
        xmid_init <- t[which.min(abs(x - mid_x))]
    }
    if (!is.finite(slope_init) || slope_init == 0) {
        slope_init <- if (t_range > 0) {
            (B_init - A_init) / t_range
        } else {
            sign(B_init - A_init) * 1
        }
    }

    if (has_asym) {
        ## empirical inflection-height fraction at point of steepest change;
        ## smooth `dx_dt` with a running mean (window ~ 1/10 of n) to
        ## suppress single-point noise spikes that would otherwise dominate
        ## `which.max`. Clamp into (0.1, 0.9) to keep init away from bounds.
        dx_dt <- diff(x) / diff(t)
        win <- max(3L, 2L * (length(dx_dt) %/% 20L) + 1L)
        dx_smooth <- as.numeric(stats::filter(
            dx_dt, rep(1 / win, win), sides = 2L
        ))
        dx_smooth[!is.finite(dx_smooth)] <- 0
        i_infl <- which.max(abs(dx_smooth))
        asym_emp <- (x[i_infl] - A_init) / (B_init - A_init)
        asym_init <- min(max(asym_emp, 0.1), 0.9)
        return(c(
            A = A_init,
            B = B_init,
            xmid = xmid_init,
            slope = slope_init,
            asym = asym_init
        ))
    } else {
        return(c(A = A_init, B = B_init, xmid = xmid_init, slope = slope_init))
    }
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
#' The 4-parameter form is recommended for small samples or when no obvious
#'   asymmetry is expected, as it converges more reliably. [stats::nls()]
#'   reads the free parameters from the formula right-hand side, so omitting
#'   `asym` incurs no degrees-of-freedom penalty.
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
#' ## 5-parameter fit
#' model5 <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
#' model5
#'
#' ## 4-parameter fit on the same data
#' model4 <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
#' model4
#'
#' y5 <- predict(model5, data)
#' y4 <- predict(model4, data)
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
    initial = logistic_init,
    parameters = c("A", "B", "xmid", "slope", "asym")
)


#' Analyse logistic kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "logistic")`. Fits a logistic curve to each
#' `nirs_channel` within a single *"mnirs"* data frame. See
#' [analyse_kinetics()] for user-facing documentation.
#'
#' @param use_asym Logical; default is `TRUE` to attempt to fit a 5-parameter
#'   [SSlogistic()] model (A, B, xmid, slope, asym) with an asymmetry
#'   parameter. If the 5-parameter fit fails, or if `use_asym = FALSE`,
#'   attempts to fit a reduced 4-parameter symmetric [SSlogistic()] model
#'   (A, B, xmid, slope).
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `xmid`, `slope`, `asym`, `xmid_fitted`.
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
#' @seealso [analyse_kinetics()], [logistic()], [SSlogistic()]
#'
#' @keywords internal
analyse_logistic <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    use_asym = TRUE,
    t0 = NULL,
    direction = c("auto", "positive", "negative"),
    end_fit_span = Inf,
    channel_args = list(),
    verbose = TRUE,
    ...
) {
    ## validation ==================================================
    validate_mnirs_data(data)
    args <- list(...)
    direction <- match.arg(direction)

    if (!(args$bypass_checks %||% FALSE)) {
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
    }
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, verbose)
    time_channel <- validate_time_channel(enquo(time_channel), data)
    if (!is.logical(use_asym) || length(use_asym) != 1L) {
        cli_abort(c(
            "x" = "{.arg use_asym} must be a {.cls logical} \\
            either {.val {TRUE}} or {.val {FALSE}}."
        ))
    }
    validate_numeric(
        end_fit_span, 1, c(0, Inf), msg1 = "one-element positive"
    )
    time_vec <- data[[time_channel]]
    t0 <- validate_t0(t0, data, time_vec, verbose)
    interval_names <- args$interval_names %||% substitute(data)

    default_args <- list(
        use_asym = use_asym,
        t0 = t0,
        direction = direction,
        end_fit_span = end_fit_span,
        verbose = verbose,
        args
    )

    ## NA scaffold for convergence failure
    na_coefs <- data.frame(
        nirs_channels = NA_character_,
        time_channel = time_channel,
        A = NA_real_,
        B = NA_real_,
        xmid = NA_real_,
        slope = NA_real_,
        asym = NA_real_,
        xmid_fitted = NA_real_
    )

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, n_params, e, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        msg <- c(
            "x" = "{n_params}-parameter {.fn SSlogistic} fit failed for \\
            {.field {(.nirs)}} in {.field {interval_names}}.",
            "!" = "{conditionMessage(e)}"
        )
        if (n_params == 5L) {
            msg <- c(msg, "i" = "Attempting 4-parameter {.fn SSlogistic} fit.")
        }
        cli_warn(msg)
        return(invisible(NULL))
    }

    ## process per-channel ============================================
    results <- lapply(nirs_channels, \(.nirs) {
        all_args <- utils::modifyList(
            default_args, channel_args[[.nirs]] %||% list()
        )
        ## derive n_params from use_asym for internal use
        n_params <- if (all_args$use_asym) 5L else 4L

        ## filter for valid finite idx before first extreme + end_fit_span
        valid <- find_kinetics_idx(
            data[[.nirs]], time_vec, all_args$end_fit_span, all_args$direction
        )
        all_args$direction <- valid$direction
        x_fit <- data[[.nirs]][valid$idx]
        t_fit <- time_vec[valid$idx]

        fit_data <- data.frame(.x = x_fit, .t = t_fit)

        ## attempt nls fit on 5-param then fall back to 4-param on failure
        model <- NULL
        if (n_params == 5L) {
            model <- tryCatch(
                nls(.x ~ SSlogistic(.t, A, B, xmid, slope, asym), fit_data),
                error = \(e) {
                    fit_failed_warning(.nirs, n_params, e, verbose)
                    NULL
                }
            )
            if (is.null(model)) n_params <- 4L
        }

        if (n_params == 4L) {
            model <- tryCatch(
                nls(.x ~ SSlogistic(.t, A, B, xmid, slope), fit_data),
                error = \(e) {
                    fit_failed_warning(.nirs, n_params, e, verbose)
                    NULL
                }
            )
        }

        if (is.null(model)) {
            return(build_na_results(.nirs, na_coefs, all_args, n_params))
        }

        fitted_vals <- stats::predict(model)
        coefs <- stats::coef(model)
        asym_arg <- if (n_params == 5L) coefs[["asym"]] else NULL
        asym_val <- asym_arg %||% NA_real_
        xmid_offset <- coefs[["xmid"]] - t0

        ## predict response at the inflection point xmid
        xmid_fitted <- logistic(
            t = coefs[["xmid"]],
            A = coefs[["A"]],
            B = coefs[["B"]],
            xmid = coefs[["xmid"]],
            slope = coefs[["slope"]],
            asym = asym_arg
        )

        coefs <- data.frame(
            nirs_channels = .nirs,
            time_channel  = time_channel,
            A             = coefs[["A"]],
            B             = coefs[["B"]],
            xmid          = xmid_offset,
            slope         = coefs[["slope"]],
            asym          = asym_val,
            xmid_fitted   = xmid_fitted
        )

        diag <- compute_diagnostics(
            x_fit, t_fit, fitted_vals, n_params, verbose
        )

        list(
            coefficients = coefs,
            model = model,
            fitted_data = data.frame(
                window_idx = valid$idx,
                fitted     = fitted_vals
            ),
            diagnostics = cbind(data.frame(nirs_channels = .nirs), diag),
            channel_args = build_channel_args(.nirs, all_args)
        )
    })

    return(build_channel_results(results, nirs_channels, t0, verbose))
}
