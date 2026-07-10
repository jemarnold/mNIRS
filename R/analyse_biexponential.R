#' Biexponential function
#'
#' Calculate a 5-parameter biexponential curve: a shared-asymptote response
#' split between a fast and a slow exponential component.
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting baseline of the response
#'   variable.
#' @param B A numeric parameter for the ending asymptote of the response
#'   variable.
#' @param tau1 A numeric parameter for the *fast* time constant (\eqn{\tau_1}),
#'   in units of the predictor variable `t`.
#' @param tau2 A numeric parameter for the *slow* time constant (\eqn{\tau_2}),
#'   in units of the predictor variable `t`. By convention `tau1 <= tau2`.
#' @param prop A numeric parameter in `[0, 1]` for the *amplitude fraction*
#'   carried by the fast component `tau1`; the slow component `tau2` carries
#'   the remaining `1 - prop`.
#'
#' @details
#' This model is fit by [analyse_kinetics()] when `method = "biexponential"`,
#' and by [stats::nls()] via the self-starting wrapper [SSbiexponential()].
#'
#' ## Model equation
#'
#' `A + (B - A) * (prop * (1 - exp(-t / tau1)) +
#'   (1 - prop) * (1 - exp(-t / tau2)))`
#'
#' The two components share the overall amplitude `B - A`; `prop` divides it
#' between the fast (`tau1`) and slow (`tau2`) time constants. `tau1` and
#' `tau2` are exchangeable, so a `tau1 <= tau2` ordering is used to make the
#' parameters identifiable. `prop` is clamped to `[1e-6, 1 - 1e-6]` internally
#' to keep the fit stable during iteration.
#'
#' The amplitude-weighted *mean response time* is
#' `MRT = prop * tau1 + (1 - prop) * tau2`.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSbiexponential()], [monoexponential()]
#'
#' @examples
#' ## create a biexponential curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(t, A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B, tau1, tau2, prop), data = data)
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
biexponential <- function(t, A, B, tau1, tau2, prop) {
    ## clamp amplitude fraction to keep the fit stable during iteration
    prop <- min(max(prop, 1e-6), 1 - 1e-6)
    y <- A + (B - A) * (
        prop * (1 - exp(-t / tau1)) + (1 - prop) * (1 - exp(-t / tau2))
    )
    return(y)
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
#'   user-fixed parameter values from [init_fixed()] used to seed the
#'   remaining free estimates.
#'
#' @returns [biexp_init()]: Initial starting estimates for parameters in the
#'   model called by [SSbiexponential()].
#'
#' @keywords internal
biexp_init <- function(mCall, data, LHS, ...) {
    tx <- stats::sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## asymptotes from first and last ceiling(n/5) values (shared helper)
    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    B_init <- fixed$B %||% ab$B

    ## estimate one overall time constant via linearisation (SSasymp method):
    ## log(B - x) ~ log(B - A) - t/tau, using shifted x to avoid log of
    ## negative/zero, then split into fast/slow components for identifiability
    x_shifted <- B_init - x
    x_pos <- x_shifted > 0
    tau_overall <- if (sum(x_pos) >= 3L) {
        x_shifted[!x_pos] <- min(x_shifted[x_pos]) / 2
        rate <- -coef(stats::lm(log(x_shifted) ~ t))[2L]
        if (is.finite(rate) && rate > 0) 1 / rate else diff(range(t)) / 3
    } else {
        diff(range(t)) / 3
    }

    return(c(
        A = A_init,
        B = B_init,
        tau1 = fixed$tau1 %||% max(tau_overall * 0.5, .Machine$double.eps),
        tau2 = fixed$tau2 %||% max(tau_overall * 2, .Machine$double.eps),
        prop = fixed$prop %||% 0.5
    ))
}


#' Self-starting biexponential model
#'
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [biexponential()], for use with [stats::nls()]. Fits the 5-parameter form
#' (A, B, tau1, tau2, prop).
#'
#' @usage
#' SSbiexponential(t, A, B, tau1, tau2, prop)
#'
#' @inheritParams biexponential
#'
#' @details
#' Model: `x ~ SSbiexponential(t, A, B, tau1, tau2, prop)`
#'
#' The fast/slow time constants `tau1` and `tau2` are exchangeable, so
#' [analyse_kinetics()] orders the fitted values `tau1 <= tau2` (swapping
#' `prop` accordingly) for a reproducible parameterisation.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its
#'   name in the formula, e.g. `x ~ SSbiexponential(t, A = 0, B, tau1, tau2,
#'   prop)` fixes the baseline at `A = 0`. Fixed parameters are excluded from
#'   estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [biexponential()], [stats::nls()], [stats::selfStart()],
#'   [SSmonoexponential()]
#'
#' @examples
#' ## create a biexponential curve with random noise
#' set.seed(1)
#' t <- 0:120
#' x <- biexponential(t, A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6) +
#'     rnorm(length(t), 0, 0.8)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSbiexponential(t, A, B, tau1, tau2, prop), data = data)
#' summary(model)
#'
#' ## fix the baseline A at a known value
#' model_fixed <- nls(
#'     x ~ SSbiexponential(t, A = 50, B, tau1, tau2, prop), data = data
#' )
#' summary(model_fixed)
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
SSbiexponential <- selfStart(
    model = biexponential,
    initial = init_fixed(biexp_init, c("A", "B", "tau1", "tau2", "prop")),
    parameters = c("A", "B", "tau1", "tau2", "prop")
)


#' Analyse biexponential kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "biexponential")`. Fits a 5-parameter
#' biexponential curve to each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @param fix An *optional* named list of model parameters (`A`, `B`, `tau1`,
#'   `tau2`, `prop`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied globally across
#'   `nirs_channels`.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `tau1`, `tau2`, `prop`, `MRT`, `tau1_fitted`,
#'   `tau2_fitted`, `MRT_fitted`. Per-channel metadata are attached as
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
            start_time = start_time,
            direction = direction,
            end_window = end_window
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose,
        env = env
    )
    nirs_channels <- setup$nirs_channels
    time_channel <- setup$time_channel
    per_channel <- setup$per_channel

    ## global fixed parameters bypass per-channel resolution: a named
    ## list would be misread as a channel map
    fix <- validate_fix(fix, c("A", "B", "tau1", "tau2", "prop"), env = env)

    ## NA scaffold (method columns only) for convergence failure
    na_coefs <- data.frame(
        A = NA_real_,
        B = NA_real_,
        tau1 = NA_real_,
        tau2 = NA_real_,
        prop = NA_real_,
        MRT = NA_real_,
        tau1_fitted = NA_real_,
        tau2_fitted = NA_real_,
        MRT_fitted = NA_real_
    )

    ## construct warning messages for fit failure
    fit_failed_warning <- function(.nirs, e, verbose) {
        if (!verbose) {
            return(invisible(NULL))
        }
        cli_warn(c(
            "x" = "{.fn SSbiexponential} fit failed for \\
            {.field {(.nirs)}} in {.field {interval_name}}.",
            "!" = "{conditionMessage(e)}"
        ), call = warn_call(env))
        return(invisible(NULL))
    }

    ## method-specific fit: self-starting biexponential via nls
    biexp_fit <- function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        fit_data <- data.frame(.x = x_fit, .t = t_fit)
        params <- c("A", "B", "tau1", "tau2", "prop")

        model <- tryCatch(
            nls(
                build_ss_formula(quote(SSbiexponential), params, fix),
                fit_data
            ),
            error = \(e) {
                fit_failed_warning(.nirs, e, verbose)
                NULL
            }
        )

        if (is.null(model)) {
            return(build_na_results(na_coefs))
        }

        coefs <- full_coefs(model, params, fix)

        ## enforce direction: bounded refit on D = B - A when inverted.
        ## data-scaled tau floors detect degenerate step fits; prop is
        ## clamped inside the model so needs no bound
        free_extra <- setdiff(params, c("A", "B", names(fix)))
        tau_free <- intersect(c("tau1", "tau2"), free_extra)
        enforced <- enforce_direction(
            model,
            coefs,
            fit_data,
            direction = .a$direction,
            amp_fn = quote(biexponential),
            extra = coefs[free_extra],
            extra_lower = if (length(tau_free) > 0L) {
                setNames(
                    rep(diff(range(t_fit)) * 1e-6, length(tau_free)),
                    tau_free
                )
            },
            fn = quote(SSbiexponential),
            fix = fix,
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

        ## enforce fast-first convention (tau1 <= tau2) when both are free;
        ## swapping relabels the exchangeable components, so prop flips too
        if (
            all(c("tau1", "tau2") %in% free_extra) &&
                coefs[["tau1"]] > coefs[["tau2"]]
        ) {
            coefs[c("tau1", "tau2")] <- coefs[c("tau2", "tau1")]
            coefs[["prop"]] <- 1 - coefs[["prop"]]
        }
        ## report the effective (clamped) amplitude fraction
        coefs[["prop"]] <- min(max(coefs[["prop"]], 1e-6), 1 - 1e-6)

        MRT_val <- coefs[["prop"]] * coefs[["tau1"]] +
            (1 - coefs[["prop"]]) * coefs[["tau2"]]

        ## predict response at tau1, tau2, and MRT (elapsed from onset)
        fitted_params <- biexponential(
            t = c(coefs[["tau1"]], coefs[["tau2"]], MRT_val),
            A = coefs[["A"]],
            B = coefs[["B"]],
            tau1 = coefs[["tau1"]],
            tau2 = coefs[["tau2"]],
            prop = coefs[["prop"]]
        )

        list(
            coefs = data.frame(
                A           = coefs[["A"]],
                B           = coefs[["B"]],
                tau1        = coefs[["tau1"]],
                tau2        = coefs[["tau2"]],
                prop        = coefs[["prop"]],
                MRT         = MRT_val,
                tau1_fitted = fitted_params[[1L]],
                tau2_fitted = fitted_params[[2L]],
                MRT_fitted  = fitted_params[[3L]]
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
        biexp_fit,
        verbose,
        interval_name,
        extra_args = c(args, list(fix = if (length(fix) > 0L) fix else NULL)),
        env = env
    ))
}
