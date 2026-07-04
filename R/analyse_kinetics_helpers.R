#' Detect the direction of a response signal
#'
#' Resolves whether a signal is predominantly increasing (`"positive"`) or
#' decreasing (`"negative"`) by computing the net slope of `x` over `t`.
#' Used internally to disambiguate peak (maximum) from trough (minimum)
#' detection when `direction = "auto"`.
#'
#' @param fallback A numeric vector (*defaults* to `x`) used to resolve
#'   direction when the net slope of `x` is zero or `NA`. The absolute
#'   maximum and minimum of `fallback` are compared; if `abs(max) >= abs(min)`,
#'   `"positive"` is returned.
#' @param direction A character string specifying the kinetics direction to
#'   detect when `"auto"` (*default*). When `"positive"` or `"negative"`
#'   returns unchanged.
#' @inheritParams replace_invalid
#'
#' @returns A character string: `"positive"` or `"negative"`.
#'
#' @keywords internal
detect_direction <- function(
    x,
    t = seq_along(x),
    fallback = x,
    direction = c("auto", "positive", "negative")
) {
    direction <- match.arg(direction)
    if (direction == "auto") {
        if (all(!is.finite(x))) {
            return("positive")
        }
        net_slope <- slope(x, t, na.rm = TRUE, bypass_checks = TRUE)

        direction <- if (is.na(net_slope) || net_slope == 0) {
            ## fallback to abs magnitude comparison when net slope is zero/NA
            max_pos <- abs(max(fallback, na.rm = TRUE))
            min_neg <- abs(min(fallback, na.rm = TRUE))
            if (max_pos >= min_neg) "positive" else "negative"
        } else if (net_slope > 0) {
            "positive"
        } else {
            "negative"
        }
    }

    return(direction)
}


#' Find valid model-fitting indices up to the first extreme
#'
#' Filters `x` and `t` to valid finite values, locates the first valid peak
#' (maximum) or trough (minimum) where `t >= 0`, and returns the integer
#' indices of all finite observations up to `end_window` past that extreme.
#'
#' @param end_window A numeric value in units of `time_channel` or `t`
#'   specifying the forward-looking window used to check for subsequent
#'   greater/lesser values than the candidate extreme. `end_window = Inf`
#'   (*default*) returns the global extreme from the full range of `x`.
#' @param direction A character string specifying the response direction to
#'   detect -- `"auto"` (*default*), `"positive"`, or `"negative"`. See
#'   *Details*.
#' @inheritParams replace_invalid
#'
#' @details
#' ## Direction detection
#'
#' When `direction = "auto"`, the net slope across all of `x` is computed
#' to determine the overall trend. If positive, the function searches for
#' a peak (maximum). If negative, it searches for a trough (minimum). When
#' the net slope is zero or `NA`, the direction is determined by comparing
#' `abs(max(x))` to `abs(min(x))`, with ties defaulting to `"positive"`.
#'
#' ## Negative time handling
#'
#' Only samples where `t >= 0` are used for detecting the extreme, allowing
#' pre-baseline (negative time) data to be excluded from the search.
#' However, indices where `t < 0` are included in the returned vector
#' provided they are finite.
#'
#' @returns A named list with three elements:
#'   \describe{
#'     \item{`direction`}{Character; the resolved direction used --
#'       `"positive"` (peak) or `"negative"` (trough).}
#'     \item{`extreme`}{Integer or `NULL`; the index of the
#'       first qualifying peak or trough in original `x`
#'       space, or `NULL` if no qualifying extreme was found
#'       (monotonic, horizontal, or degenerate input).}
#'     \item{`idx`}{Integer vector of all valid finite indices,
#'       truncated at `t[extreme] + end_window`.}
#'   }
#'
#' @inheritParams validate_mnirs
#' @keywords internal
find_kinetics_idx <- function(
    x,
    t = seq_along(x),
    end_window = Inf,
    direction = c("auto", "positive", "negative"),
    ...,
    env = rlang::caller_env()
) {
    ## validation =================================================
    args <- list(...)
    n <- length(x)
    direction <- match.arg(direction)

    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t, allow_na = TRUE, env = env)
        validate_numeric(
            end_window, 1, c(0, Inf), msg1 = "one-element positive",
            env = env
        )
    }

    ## all finite indices (including t < 0)
    which_valid <- which(is.finite(x) & is.finite(t))
    ## subset where t >= 0 for extreme detection
    which_positive <- which_valid[t[which_valid] >= 0]
    x_valid <- x[which_positive]
    t_valid <- t[which_positive]
    n_valid <- length(x_valid)

    invalid_return <- list(
        direction = "positive",
        extreme = NULL,
        idx = which_valid
    )

    ## early returns for degenerate inputs
    if (n_valid < 2L) {
        return(invalid_return)
    }

    ## direction detection, fallback to abs magnitude
    direction <- detect_direction(x, t, x, direction)
    extreme_fn <- if (direction == "positive") max else min
    compare_fn <- if (direction == "positive") `>=` else `<=`
    which_fn <- if (direction == "positive") which.max else which.min
    invalid_return$direction <- direction

    ## monotonic: if global extreme is the last x value
    ## horizontal: if all x values equal
    if (which_fn(x_valid) == n_valid || all(x_valid == x_valid[1L])) {
        return(invalid_return)
    }

    ## process ==================================================
    ## bin by end_window, find extreme per bin, then check forward window

    ## expand end_window to global extreme
    if (end_window == Inf) {
        end_window <- t_valid[length(t_valid)]
    }

    ## ensure end_window covers at least one adjacent sample when
    ## end_window is smaller than the minimum time step
    t_diff <- diff(t_valid)
    mod_span <- max(end_window, min(t_diff[t_diff > 0]))

    ## break t_valid into mod_span-width bins
    bin_breaks <- seq(t_valid[1L], t_valid[n_valid] + mod_span, mod_span)
    bin_idx <- findInt_mnirs(
        t_valid,
        bin_breaks,
        rightmost.closed = TRUE,
        env = env
    )

    ## extreme index per bin (first occurrence for ties)
    bin_extreme_idx <- unname(tapply(seq_len(n_valid), bin_idx, \(.idx) {
        .idx[which_fn(x_valid[.idx])]
    }, simplify = TRUE))

    ## find first bin-extreme where no forward mod_span value exceeds
    ## add tolerance for where mod_span loses floating point precision
    ## this is idx of x_valid, not `x`
    extreme_idx <- Find(\(.idx) {
        in_window <- t_valid >= t_valid[.idx] &
            t_valid <= t_valid[.idx] + mod_span + .Machine$double.eps^0.5
        compare_fn(x_valid[.idx], extreme_fn(x_valid[in_window]))
    }, bin_extreme_idx)

    if (!is.null(extreme_idx)) {
        ## map extreme back to original index space
        orig_extreme <- which_positive[extreme_idx]
        t_cutoff <- t[orig_extreme] + end_window
        truncated <- which_valid[t[which_valid] <= t_cutoff]

        return(list(
            direction = direction,
            extreme = orig_extreme,
            idx = truncated
        ))
    }

    ## fallback: no qualifying extreme found
    return(invalid_return)
}


#' Gather per-interval `mnirs_kinetics` into results structure
#'
#' Shared helper for `analyse_kinetics.*` methods. Takes a list of
#' per-interval (per-data frame) kinetics results data frames (each carrying
#' `"fitted_data"`, `"channel_args"`, and `"diagnostics"` attributes) and the
#' original `data_list`. Interval names are read from the `interval` column
#' that each result data frame already carries.
#'
#' @param data_list Named list of original interval data frames.
#' @param result_list List of per-interval result data frames with attributes.
#'
#' @returns A named list with: `coefficients`, `model`, `data`,
#'   `interval_times`, `diagnostics`, `channel_args`.
#' @keywords internal
build_kinetics_results <- function(
    data_list,
    result_list,
    method,
    call
) {
    ## interval label per result, taken from each result's `interval` column
    interval_names <- vapply(result_list, \(.r) unique(.r$interval), "")

    ## tag each channel-level attr df with its interval, then bind rows
    flatten_attr <- function(attr_name) {
        dfs <- Map(\(.res, .interval) {
            df <- attr(.res, attr_name)
            df$interval <- .interval
            return(df)
        }, result_list, interval_names)
        df <- do.call(rbind, dfs)
        df <- df[, c("interval", setdiff(names(df), "interval"))]
        rownames(df) <- NULL
        return(df)
    }

    channel_args <- flatten_attr("channel_args")
    diagnostics <- flatten_attr("diagnostics")

    ## augment `data_list` dfs with `<nirs_channels>_fitted` columns
    fitted_data_list <- Map(\(.df, .result) {
        ## extract fitted columns from "fitted_data" per `nirs_channel`
        fitted_data <- attr(.result, "fitted_data")
        fitted_cols <- lapply(fitted_data, \(.pred) {
            fitted_vec <- rep(NA_real_, nrow(.df))
            fitted_vec[.pred$window_idx] <- .pred$fitted
            fitted_vec
        })
        names(fitted_cols) <- paste0(names(fitted_data), "_fitted")
        ## agument `<nirs_channels>_fitted` columns to df
        augmented <- cbind(.df, as.data.frame(fitted_cols))
        
        ## metadata ==================================================
        metadata <- attributes(.df)
        metadata$nirs_channels <- unique(.result$nirs_channels)
        metadata$time_channel <- unique(.result$time_channel)
        create_mnirs_data(augmented, metadata)
    }, data_list, result_list)

    ## extract interval_times from each data_list attributes, if exist
    interval_times_df <- data.frame(interval = names(data_list))
    interval_times_df$interval_times <- lapply(data_list, \(.df) {
        interval_times <- attr(.df, "interval_times")
        if (is.null(interval_times)) NA_real_ else unlist(interval_times)
    })
    ## add `class = "mnirs"` for `plot.mnirs`
    class(fitted_data_list) <- c("mnirs", class(fitted_data_list))

    ## extract per-interval model lists (named by nirs_channel)
    model_list <- lapply(result_list, attr, "model")
    names(model_list) <- interval_names

    ## normalise call: function name to generic, method to canonical form
    call[[1L]] <- quote(analyse_kinetics)
    if ("method" %in% names(call)) {
        call$method <- method
    }

    ## combine scalar coefficients & relocate interval col to col[1]
    coefs <- do.call(rbind, result_list)
    coefs <- coefs[, c("interval", setdiff(names(coefs), "interval"))]
    rownames(coefs) <- NULL

    return(structure(
        list(
            method = method,
            model = model_list,
            coefficients = coefs,
            data = fitted_data_list,
            interval_times = interval_times_df,
            diagnostics = diagnostics,
            channel_args = channel_args,
            call = call
        ),
        class = "mnirs_kinetics"
    ))
}


#' Process kinetics fits across NIRS channels
#'
#' Shared per-channel skeleton for all `analyse_kinetics()` methods. Resolves
#' the fitting window for each channel via [find_kinetics_idx()], delegates
#' the method-specific fit to `fit_fn`, and assembles the standard attributed
#' result. Method workers supply a validated `per_channel` argument list and a
#' `fit_fn`; everything else is common.
#'
#' @param data A single *"mnirs"* data frame.
#' @param nirs_channels Character vector of resolved channel names.
#' @param time_channel Character; resolved time column name.
#' @param per_channel Named list (one element per channel) of resolved and
#'   validated argument lists from [resolve_channel_args()] and
#'   [validate_kinetics_args()].
#' @param fit_fn A function `(.nirs, x_fit, t_fit, .a, valid, verbose)`
#'   returning a list with `coefs` (1-row data frame of method coefficients,
#'   *without* `nirs_channels`/`time_channel`), `model`, `fitted_data`
#'   (`window_idx`/`fitted`), and `diag` (1-row data frame from
#'   [compute_diagnostics()]).
#' @param interval_name Character; the interval name recorded in the `interval`
#'   column of the returned coefficients.
#' @param extra_args Named list of additional arguments recorded in the
#'   `channel_args` result attribute.
#' @inheritParams validate_mnirs
#'
#' @returns A `data.frame` of coefficients (columns `interval`,
#'   `nirs_channels`, `time_channel`, and method parameters), one row per
#'   channel, with attributes `"model"` and `"fitted_data"` (named lists by
#'   channel) and `"diagnostics"` and `"channel_args"` (data frames, one row
#'   per channel).
#'
#' @keywords internal
analyse_kinetics_channels <- function(
    data,
    nirs_channels,
    time_channel,
    per_channel,
    fit_fn,
    verbose = TRUE,
    interval_name = NA_character_,
    extra_args = list(),
    env = rlang::caller_env()
) {
    t_vec <- data[[time_channel]]

    ## per-channel fit; collect parallel pieces keyed by channel
    fits <- setNames(
        lapply(nirs_channels, \(.nirs) {
        .a <- per_channel[[.nirs]]

        ## filter for valid finite idx before first extreme + end_window
        valid <- find_kinetics_idx(
            data[[.nirs]], t_vec, .a$end_window, .a$direction, env = env
        )
        .a$direction <- valid$direction
        x_fit <- data[[.nirs]][valid$idx]
        t_fit <- t_vec[valid$idx]

        ## method-specific fit; coefs/diag carry method columns only
        fit <- fit_fn(.nirs, x_fit, t_fit, .a, valid, verbose)

        ## serialise resolved args: NULL to NA, list() to its deparse(),
        ## dropping internal-only args, so they fit a flat data frame row
        arg_row <- lapply(c(.a, extra_args), \(.x) {
            if (is.null(.x)) NA else if (is.list(.x)) deparse(.x) else .x
        })
        arg_row[c("verbose", "bypass_checks")] <- NULL

        list(
            coefficients = cbind(
                data.frame(
                    interval      = interval_name,
                    nirs_channels = .nirs,
                    time_channel  = time_channel
                ),
                fit$coefs
            ),
            model        = fit$model,
            fitted_data  = fit$fitted_data,
            diagnostics  = cbind(data.frame(nirs_channels = .nirs), fit$diag),
            channel_args = data.frame(nirs_channels = .nirs, arg_row)
        )
    }),
        nirs_channels
    )

    ## assemble single attributed df (consumed by build_kinetics_results)
    result <- structure(
        do.call(rbind, lapply(fits, `[[`, "coefficients")),
        model        = lapply(fits, `[[`, "model"),
        fitted_data  = lapply(fits, `[[`, "fitted_data"),
        diagnostics  = do.call(rbind, lapply(fits, `[[`, "diagnostics")),
        channel_args = do.call(rbind, lapply(fits, `[[`, "channel_args"))
    )

    ## warn when time coefficients are negative (response before start_time)
    if (verbose) {
        check_cols <- intersect(
            c("TD", "tau", "response_time", "peak_slope_time"),
            names(result)
        )

        if (any(unlist(result[check_cols]) < 0, na.rm = TRUE)) {
            cli_warn(c(
                "!" = "Negative {.arg time_channel} coefficients imply the \\
                response occured before {.arg start_time}. This may \\
                indicate a poorly fitted or misparameterised model.",
                "i" = "Check {.arg time_channel} and {.arg start_time} \\
                values, or consider using a different \\
                {.fn analyse_kinetics} method."
            ), call = warn_call(env))
        }
    }

    return(result)
}


#' Validate resolved per-channel kinetics arguments
#'
#' Validates each channel's resolved argument list once, before any fitting,
#' so an invalid argument fails fast rather than after an expensive fit on an
#' earlier channel. Validation is keyed on which arguments are present.
#' Mutating validators are applied and written back:
#' [validate_start_time()] clamps `start_time`, and `align` is matched to its
#' choices. Verbose hints are emitted for the first channel only to avoid
#' repeating identical messages.
#'
#' @param per_channel Named list of resolved argument lists, one per channel.
#' @param t_vec Numeric vector of `time_channel` values.
#' @inheritParams validate_mnirs
#'
#' @returns The `per_channel` list with mutating validators applied.
#'
#' @keywords internal
validate_kinetics_args <- function(
    per_channel,
    data,
    t_vec,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    chans <- names(per_channel)
    setNames(
        lapply(chans, \(.nirs) {
        .a <- per_channel[[.nirs]]
        ## emit verbose hints once, for the first channel only
        .v <- verbose && .nirs == chans[[1L]]

        if (
            !is.null(.a$use_TD) &&
                (!is.logical(.a$use_TD) ||
                    length(.a$use_TD) != 1L)
            ) {
            cli_abort(c(
                "x" = "{.arg use_TD} must be a {.cls logical} \\
                either {.val {TRUE}} or {.val {FALSE}}."
            ), call = env)
        }
        ## always resolve: user value > interval_times metadata > 0
        .a$start_time <- validate_start_time(
            .a$start_time, data, t_vec, .v, env = env
        )
        if (!is.null(.a$end_window)) {
            validate_numeric(
                .a$end_window, 1, c(0, Inf), 
                msg1 = "one-element positive", env = env
            )
        }
        if (!is.null(.a$fraction)) {
            validate_numeric(
                .a$fraction, 1L, c(0, 1), 
                msg2 = "between {col_blue('[0, 1]')}.", env = env
            )
        }
        if (any(c("width", "span") %in% names(.a))) {
            validate_width_span(.a$width, .a$span, .v, env = env)
        }
        if (!is.null(.a$align)) {
            .a$align <- sub("^center$", "centre", .a$align[[1L]])
            .a$align <- rlang::arg_match0(
                .a$align, c("centre", "left", "right"), arg_nm = "align",
                error_call = env
            )
        }
        .a
    }),
        chans
    )
}


#' Build a standardised NA result for a failed channel
#'
#' Returns the method-specific `coefs`/`model`/`fitted_data`/`diag` list
#' expected by [analyse_kinetics_channels()] when a model fit fails,
#' populated with `NA`/`NULL` values.
#'
#' @param na_coefs A template 1-row `data.frame` of `NA` method coefficients
#'   (*without* `nirs_channels`/`time_channel`, which are added upstream).
#'
#' @returns A named list with elements `coefs`, `model`, `fitted_data`,
#'   and `diag`.
#'
#' @keywords internal
build_na_results <- function(na_coefs) {
    na_diag <- data.frame(
        n_obs = 0L,
        r2 = NA_real_,
        adj_r2 = NA_real_,
        rmse = NA_real_,
        snr = NA_real_,
        cv_rmse = NA_real_,
        aic = NA_real_,
        aicc = NA_real_,
        bic = NA_real_
    )
    return(list(
        coefs = na_coefs,
        model = NULL,
        fitted_data = data.frame(window_idx = NA_integer_, fitted = NA_real_),
        diag = na_diag
    ))
}


#' Refit an amplitude-reparameterised model with direction box bounds
#'
#' Enforces the response direction on a converged but inverted nls fit by
#' refitting with amplitude `D = B - A` bounded to the requested sign via
#' `nls(algorithm = "port")`. Sigmoid models divide by `D`, so its
#' magnitude is floored strictly above zero.
#'
#' @param amp_fn Symbol; exported model fn taking `(t, A, B, ...)`.
#' @param fit_data Data frame with columns `.x` and `.t`.
#' @param A,D0 Numeric start values; `D0` sign gives the requested
#'   direction.
#' @param extra Named numeric start values for remaining free params, in
#'   `amp_fn` argument order after `B`.
#' @param extra_lower,extra_upper Named numeric bound overrides for
#'   `extra` params. Sign-floor bounds should be data-scaled small values
#'   (not `.Machine$double.eps`) so pinned-floor degeneracy is detectable.
#'
#' @returns A named list `list(model, coefs)` with `coefs` in
#'   `(A, B, ...)` space, or `NULL` when the refit fails or any
#'   coefficient is pinned at a sign-floor bound (degenerate flat fit).
#'
#' @keywords internal
refit_direction <- function(
    amp_fn,
    fit_data,
    A,
    D0,
    extra,
    extra_lower = NULL,
    extra_upper = NULL
) {
    want <- sign(D0)
    D_eps <- diff(range(fit_data$.x)) * 1e-6

    ## build rhs: amp_fn(.t, A, A + D, <extra names>)
    rhs <- as.call(c(
        amp_fn,
        quote(.t),
        quote(A),
        call("+", quote(A), quote(D)),
        lapply(names(extra), as.name)
    ))
    nls_formula <- stats::as.formula(call("~", quote(.x), rhs))

    ## box bounds: D sign-constrained with strictly positive magnitude
    ## floor (sigmoid models divide by D); extras unbounded unless
    ## overridden
    start <- c(A = A, D = D0, extra)
    lower <- c(A = -Inf, D = if (want > 0) D_eps else -Inf,
               setNames(rep(-Inf, length(extra)), names(extra)))
    upper <- c(A = Inf, D = if (want > 0) Inf else -D_eps,
               setNames(rep(Inf, length(extra)), names(extra)))
    lower[names(extra_lower)] <- extra_lower
    upper[names(extra_upper)] <- extra_upper

    model <- tryCatch(
        nls(nls_formula, fit_data, start = start, lower = lower,
            upper = upper, algorithm = "port"),
        error = \(e) NULL
    )

    if (is.null(model)) {
        return(NULL)
    }

    ## any coefficient pinned at a sign-floor bound (e.g. D, slope, tau)
    ## indicates a degenerate flat fit: no genuine response in the
    ## requested direction
    cf <- coef(model)
    floors <- abs(ifelse(is.finite(lower), lower, upper))
    if (any(is.finite(floors) & abs(cf) <= 2 * floors)) {
        return(NULL)
    }

    coefs <- as.list(cf)
    coefs$B <- coefs$A + coefs$D
    coefs$D <- NULL
    return(list(model = model, coefs = coefs))
}


#' Warn when a fit converged against the requested direction
#'
#' Emitted when the direction-bounded refit from [refit_direction()]
#' fails or degenerates, before returning [build_na_results()].
#'
#' @param fn Symbol or character; the self-start fn named in the warning.
#' @param .nirs Character; the channel name.
#' @param direction Character; the requested direction.
#' @param interval_name Character; the interval label.
#' @inheritParams validate_mnirs
#'
#' @keywords internal
wrong_direction_warning <- function(
    fn,
    .nirs,
    direction,
    interval_name,
    verbose,
    env
) {
    if (!verbose) {
        return(invisible(NULL))
    }
    cli_warn(c(
        "x" = "{.fn {as.character(fn)}} fit for {.field {(.nirs)}} in \\
        {.field {interval_name}} could not satisfy \\
        {.code direction = {.val {direction}}}.",
        "i" = "Returning {.val {NA}} coefficients."
    ), call = warn_call(env))
    return(invisible(NULL))
}


#' Compute model diagnostics
#'
#' @param fitted A numeric vector of the predicted values.
#' @param n_params Integer; total number of estimated coefficients in the
#'   model (default `1L`). For linear models pass the number of regression
#'   coefficients (e.g. `2L` for `lm(x ~ t)`). For non-linear models
#'   (`"monoexponential"`, `"sigmoidal"`), pass the number of free parameters
#'   fit by the solver.
#' @inheritParams peak_slope
#' @inheritParams validate_mnirs
#'
#' @details
#'
#' ## r2
#'
#' Squared Pearson correlation between observed and fitted values. Equals
#'   the classic `1 - SSres / SStot` for OLS linear fits (matches
#'   `summary(lm)$r.squared`); a bounded `[0, 1]` pseudo-R² for non-linear
#'   fits such as `"monoexponential"` and `"sigmoidal"`.
#'
#' ## adj_r2
#'
#' Adjusted `R^2` penalised by `n_params`. Appropriate for OLS linear models;
#'   interpret with caution for non-linear fits.
#'
#' ## aic, aicc, bic
#'
#' Information criteria derived from a Gaussian log-likelihood with the
#'   maximum-likelihood residual variance `sigma_hat^2 = SSres / n_obs`. The
#'   effective parameter count is `k = n_params + 1` (the `+1` accounts for
#'   the estimated residual variance). Values match `stats::AIC()` and
#'   `stats::BIC()` for `lm` and `nls` fits. `aicc` is the small-sample
#'   correction and is `NA` when `n_obs - k - 1 <= 0`.
#'
#' @returns A 1-row `data.frame` with columns `n_obs`, `r2`, `adj_r2`,
#'   `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`, and `bic`.
#'
#' @keywords internal
compute_diagnostics <- function(
    x,
    t,
    fitted,
    n_params = 1L,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    n_obs <- length(fitted)

    return_na <- data.frame(
        n_obs = n_obs,
        r2 = NA_real_,
        adj_r2 = NA_real_,
        rmse = NA_real_,
        snr = NA_real_,
        cv_rmse = NA_real_,
        aic = NA_real_,
        aicc = NA_real_,
        bic = NA_real_
    )

    if (n_params < 1L || n_obs < 2L) {
        return(return_na)
    }

    if (length(x) != length(t) || length(x) != n_obs) {
        if (verbose) {
            cli_warn(c(
                "!" = "{.arg x}, {.arg t}, and {.arg fitted} must be \\
                {.cls numeric} vectors of equal lengths to return model \\
                diagnostics."
            ), call = warn_call(env))
        }
        return(return_na)
    }

    ## residuals and sums of squares (reused throughout)
    x_mean <- mean(x)
    resid <- x - fitted
    ss_res <- sum(resid^2)
    ss_tot <- sum((x - x_mean)^2)
    ss_fit <- sum((fitted - mean(fitted))^2)

    ## R²: squared Pearson correlation; equals 1 - SSres/SStot for OLS,
    ## bounded [0, 1] pseudo-R² for non-linear fits
    r2 <- if (ss_tot == 0 || ss_fit == 0) {
        NA_real_
    } else {
        stats::cor(x, fitted)^2
    }

    ## adjusted R²: penalised by n_params; valid for OLS linear models
    adj_r2 <- if (is.na(r2) || n_obs <= n_params) {
        NA_real_
    } else {
        1 - (1 - r2) * (n_obs - 1L) / (n_obs - n_params)
    }

    ## RMSE
    rmse <- sqrt(ss_res / n_obs)

    ## SNR: signal variance to residual variance, in dB
    snr <- if (ss_tot == 0 || ss_res == 0) {
        NA_real_
    } else {
        10 * log10(ss_tot / ss_res)
    }

    ## CV-RMSE: RMSE normalised by the absolute mean of observed values
    cv_rmse <- if (x_mean == 0) NA_real_ else rmse / abs(x_mean)

    ## information criteria from Gaussian log-likelihood; k includes the
    ## estimated residual variance, so matches stats::AIC()/BIC() for
    ## lm and nls fits
    k <- n_params + 1L
    if (ss_res <= 0) {
        aic <- aicc <- bic <- NA_real_
    } else {
        log_lik <- -n_obs / 2 * (log(2 * pi) + log(ss_res / n_obs) + 1)
        aic <- -2 * log_lik + 2 * k
        bic <- -2 * log_lik + log(n_obs) * k
        aicc <- if (n_obs - k - 1L <= 0L) {
            NA_real_
        } else {
            aic + 2 * k * (k + 1L) / (n_obs - k - 1L)
        }
    }

    return(data.frame(
        n_obs = n_obs,
        r2 = r2,
        adj_r2 = adj_r2,
        rmse = rmse,
        snr = snr,
        cv_rmse = cv_rmse,
        aic = aic,
        aicc = aicc,
        bic = bic
    ))
}


#' Update a model object with Fixed coefficients
#'
#' Re-fit a model with fixed coefficients provided as additional arguments.
#' Fixed coefficients are not modified when optimising for best fit.
#'
#' @param model An existing model object from `lm`, `nls`, `glm`, and others.
#' @param data An *optional* data frame to supply manually if original data
#'   frame is unavailable from a different parent environment.
#' @param ... Named model coefficients to fix.
#' @inheritParams validate_mnirs
#'
#' @details
#' If no fixed coefficients are supplied, or if a coefficient does not exist
#'   in the model, the model will be returned unchanged (with a warning).
#'
#' The function cannot update if all model coefficients are supplied as fixed,
#'   and will abort.
#'
#' @returns An updated model object with remaining free coefficients.
#'
#' @keywords internal
fix_coefs <- function(
    model,
    data = NULL,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    current_coefs <- coef(model)
    fixed_coefs <- list(...)
    fixed_names <- names(fixed_coefs)
    current_names <- names(current_coefs)

    ## validate coefs
    invalid <- setdiff(fixed_names, current_names)
    if (verbose && length(invalid) > 0) {
        cli_warn(c(
            "x" = "Unknown model coefficient{?s}: {.val {invalid}}.",
            "i" = "Returning model with known coefficients."
        ), call = warn_call(env))
    }

    ## extract data from the model environment
    if (is.null(data)) {
        data <- tryCatch(
            eval(model$call$data, envir = environment(stats::formula(model))),
            error = \(e) {
                ## fallback: try parent frames
                eval(model$call$data, envir = parent.frame(3))
            }
        )

        if (is.null(data)) {
            cli_abort(c(
                "x" = "Cannot retrieve original model data frame."
            ), call = env)
        }
    }

    ## get coef list from model and update in place from fixed coefs
    ## remove fixed coef from the start list
    start_coefs <- current_coefs[!current_names %in% fixed_names]

    if (length(start_coefs) == 0) {
        cli_abort(c(
            "x" = "Cannot update the model if all parameters are fixed. \\
            Nothing to estimate."
        ), call = env)
    }

    ## substitute fixed params into model_formula
    new_formula <- do.call(substitute, list(stats::formula(model), fixed_coefs))

    ## update the model
    return(stats::update(
        model,
        formula = new_formula,
        start = start_coefs,
        data = data
    ))
}
