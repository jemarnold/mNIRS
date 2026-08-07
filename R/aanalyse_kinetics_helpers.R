## canonical method for each accepted alias, matched case- and
## separator-insensitively (" ", "-", "_")
method_aliases <- c(
    response_time = "response_time",
    half_response_time = "response_time",
    recovery_time = "response_time",
    half_recovery_time = "response_time",
    half_time = "response_time",
    hrt = "response_time",
    peak_slope = "peak_slope",
    slope = "peak_slope",
    lm = "peak_slope",
    exp = "monoexponential",
    exponential = "monoexponential",
    mrt = "monoexponential",
    tau = "monoexponential",
    biexponential = "biexponential",
    biexp = "biexponential",
    double_exponential = "biexponential",
    logistic = "sigmoidal",
    gompertz = "sigmoidal",
    xmid = "sigmoidal"
) 

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
#' @param direction A character string specifying the response direction to
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
#' @param direction A character string specifying the response direction
#'   `"positive"`, or `"negative"`, or detect with `"auto"` (*default*). See
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
            end_window, 1, c(0, Inf), msg1 = "one-element positive", env = env
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
    bin_idx <- validate_findInt(
        t_valid, bin_breaks, rightmost.closed = TRUE, env = env
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
#' original `data_list`. Interval names are taken from `names(data_list)`.
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
    ## bind channel-level attr dfs; each already carries `interval` col[1]
    flatten_attr <- function(attr_name) {
        df <- do.call(rbind, lapply(result_list, attr, attr_name))
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
        ## augment `<nirs_channels>_fitted` columns to df
        augmented <- cbind(.df, as.data.frame(fitted_cols))
        
        ## metadata ==================================================
        metadata <- attributes(.df)
        metadata$nirs_channels <- unique(.result$nirs_channels)
        metadata$time_channel <- unique(.result$time_channel)
        create_mnirs_data(augmented, metadata)
    }, data_list, result_list)

    ## start_times: the resolved fit onset (user value > interval_times
    ## metadata > first time), uniform across channels within an interval
    it_df <- data.frame(
        interval = names(data_list),
        start_times = vapply(result_list, \(.r) {
            attr(.r, "channel_args")$start_time[[1L]]
        }, numeric(1), USE.NAMES = FALSE)
    )
    ## end_times sourced from interval_times metadata when any is present
    it_meta <- lapply(data_list, \(.df) {
        it <- attr(.df, "interval_times")
        if (is.null(it)) NA_real_ else unlist(it)
    })
    if (any(lengths(it_meta) >= 2L)) {
        it_df$end_times <- vapply(
            it_meta, `[`, numeric(1), 2L, USE.NAMES = FALSE
        )
    }
    ## add `class = "mnirs"` for `plot.mnirs`
    class(fitted_data_list) <- c("mnirs", class(fitted_data_list))

    ## extract per-interval model lists (named by nirs_channel)
    model_list <- lapply(result_list, attr, "model")
    names(model_list) <- names(data_list)

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
            interval_times = it_df,
            diagnostics = diagnostics,
            channel_args = channel_args,
            call = call
        ),
        class = "mnirs_kinetics"
    ))
}


#' Run a kinetics worker over each interval and collate results
#'
#' Shared skeleton for `analyse_kinetics.*` methods: normalises `data`
#' to a named list of interval data frames, calls the method worker
#' once per interval, and collates results via
#' [build_kinetics_results()].
#'
#' @param data A data frame, list of data frames, or grouped data frame.
#' @param worker Function; the interval-level worker, e.g.
#'   [analyse_monoexponential()].
#' @param method Character; the canonical method name.
#' @param worker_args Named list of method-specific arguments passed
#'   to `worker`.
#' @param nirs_quo,time_quo Quosures of the caller's `nirs_channels`
#'   and `time_channel` arguments, captured in the method frame.
#' @param call The matched call from the user-facing method.
#' @param env The call recorded for condition reporting.
#' @inheritParams validate_mnirs
#'
#' @returns An *"mnirs_kinetics"* object from
#'   [build_kinetics_results()].
#'
#' @keywords internal
analyse_kinetics_intervals <- function(
    data,
    worker,
    method,
    worker_args,
    nirs_quo,
    time_quo,
    verbose,
    call,
    env
) {
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        rlang::inject(worker(
            data = data_list[[.i]],
            nirs_channels = !!nirs_quo,
            time_channel = !!time_quo,
            !!!worker_args,
            verbose = verbose,
            interval_name = names(data_list)[[.i]],
            bypass_checks = TRUE,
            env = env
        ))
    })

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(data_list, result_list, method, call))
}


## canonical method -> method-specific worker arg names.
kinetics_dispatch <- list(
    common = c("start_time", "direction", "end_window"),
    response_time = c("fraction"),
    peak_slope = c("width", "span", "align", "partial", "na.rm"),
    monoexponential = c("use_TD", "fix"),
    biexponential = c("use_TD", "fix", "tau_ratio"),
    sigmoidal = c("shape", "fix")
)


#' Shared validation prologue for `analyse_<method>()` workers
#'
#' Runs the identical per-interval setup shared by every kinetics worker:
#' validates `data`, resolves `nirs_channels` and `time_channel`,
#' broadcasts global arguments across channels via [resolve_channel_args()],
#' and validates the resolved per-channel arguments via
#' [validate_kinetics_args()].
#'
#' @param data A single *"mnirs"* data frame.
#' @param nirs_quo,time_quo Quosures of the worker's `nirs_channels` and
#'   `time_channel` arguments (captured with `enquo()` in the worker frame).
#' @param arg_list Named list of the method's per-channel-capable arguments.
#' @param choices Named list of valid values for choice-type arguments,
#'   passed to [resolve_channel_args()].
#' @inheritParams validate_mnirs
#'
#' @returns A named list with `nirs_channels`, `time_channel`, `t_vec`, and
#'   `per_channel`.
#'
#' @keywords internal
setup_kinetics_worker <- function(
    data,
    nirs_quo,
    time_quo,
    arg_list,
    choices = list(),
    verbose = TRUE,
    env = rlang::caller_env()
) {
    validate_mnirs_data(data, env = env)
    nirs_channels <- validate_nirs_channels(nirs_quo, data, env)
    time_channel <- validate_time_channel(time_quo, data, env = env)
    t_vec <- data[[time_channel]]

    ## broadcast global args (applying per-channel list() overrides), then
    ## validate the resolved args once, before fitting any channel
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = arg_list,
        choices = choices,
        verbose = verbose,
        env = env
    )
    per_channel <- validate_kinetics_args(
        per_channel,
        data,
        t_vec,
        verbose,
        env = env
    )

    return(list(
        nirs_channels = nirs_channels,
        time_channel = time_channel,
        t_vec = t_vec,
        per_channel = per_channel
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
#'   [compute_diagnostics()]). `t_fit` is always time elapsed from
#'   `start_time`, so time coefficients need no further offset, and
#'   `window_idx` must index the original data frame rows.
#' @param interval_name Character; the interval name recorded in the `interval`
#'   column of the returned coefficients, `diagnostics`, and `channel_args`.
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
        nm = nirs_channels,
        lapply(nirs_channels, \(.nirs) {
            .a <- per_channel[[.nirs]]

            ## filter for valid finite idx before first extreme + end_window;
            ## data columns and `end_window` are already validated upstream
            valid <- find_kinetics_idx(
                data[[.nirs]],
                t_vec,
                .a$end_window,
                .a$direction,
                bypass_checks = TRUE,
                env = env
            )
            .a$direction <- valid$direction
            x_fit <- data[[.nirs]][valid$idx]
            ## fit on time elapsed from onset
            t_fit <- t_vec[valid$idx] - .a$start_time

            ## method-specific fit; coefs/diag carry method columns only
            fit <- fit_fn(.nirs, x_fit, t_fit, .a, valid, verbose)

            ## serialise resolved args: NULL to NA, list() to its deparse(),
            ## dropping internal-only args, so they fit a flat data frame row
            arg_row <- lapply(c(.a, extra_args), \(.x) {
                if (is.null(.x)) NA else if (is.list(.x)) deparse(.x) else .x
            })
            arg_row[c("verbose", "bypass_checks", "interval_name")] <- NULL

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
                diagnostics  = cbind(
                    data.frame(interval = interval_name, nirs_channels = .nirs),
                    fit$diag
                ),
                channel_args = data.frame(
                    interval      = interval_name,
                    nirs_channels = .nirs,
                    arg_row
                )
            )
        })
    )

    ## assemble single attributed df (consumed by build_kinetics_results)
    result <- structure(
        do.call(rbind, lapply(fits, `[[`, "coefficients")),
        model = lapply(fits, `[[`, "model"),
        fitted_data = lapply(fits, `[[`, "fitted_data"),
        diagnostics = do.call(rbind, lapply(fits, `[[`, "diagnostics")),
        channel_args = do.call(rbind, lapply(fits, `[[`, "channel_args"))
    )

    ## warn when time coefficients are negative (response before start_time)
    if (verbose) {
        check_cols <- intersect(
            c("TD", "tau", "tau1", "tau2", "response_time", "peak_slope_time"),
            names(result)
        )

        if (any(unlist(result[check_cols]) < 0, na.rm = TRUE)) {
            cli_warn(c(
                "!" = "Negative {.arg time_channel} coefficients imply the \\
                response occurred before {.arg start_time}. This may \\
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
        n_params = NA_integer_,
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


#' Resolve fixed parameters from a self-start model call
#'
#' Classifies each self-start model parameter in a matched call as free
#' (written as its own bare symbol) or fixed (written as any other
#' value, e.g. `A = 0`). Fixed expressions are evaluated for use as
#' initialisation seeds; values that cannot be resolved to a finite
#' numeric scalar return `NULL` and seeds fall back to data-driven
#' estimates.
#'
#' @param mCall A matched call to the `selfStart` model.
#' @param params Character vector of the model parameter names.
#' @param data A data frame with the model variables.
#'
#' @returns A named list of fixed parameter values, empty when no
#'   parameters are fixed.
#'
#' @keywords internal
resolve_fixed_params <- function(mCall, params, data) {
    present <- params[params %in% names(mCall)]
    fixed <- present[
        !vapply(present, \(.p) {
            identical(mCall[[.p]], as.name(.p))
        }, logical(1))
    ]
    return(lapply(setNames(nm = fixed), \(.p) {
        val <- tryCatch(eval(mCall[[.p]], data), error = \(e) NULL)
        if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
            val
        } else {
            NULL
        }
    }))
}


#' Wrap a self-start initialiser to support fixed parameters
#'
#' Decorates a `selfStart` `initial` function so parameters supplied as
#' values in the model formula (e.g. `SSmonoexponential(t, A = 0, B, tau)`) are
#' excluded from the returned start vector. [stats::nls()] reads the
#' free parameters from the names of that vector, so excluded
#' parameters are treated as constants in the formula. Fixed values are
#' forwarded to the wrapped initialiser as a `fixed` list argument to
#' seed the remaining free estimates.
#'
#' @param init A `selfStart` initial function `(mCall, data, LHS, ...)`.
#' @param params Character vector of the model parameter names.
#'
#' @returns A function suitable for the `initial` argument of
#'   [stats::selfStart()].
#'
#' @keywords internal
init_fixed <- function(init, params) {
    function(mCall, data, LHS, ...) {
        fixed <- resolve_fixed_params(mCall, params, data)
        start <- init(mCall, data, LHS, fixed = fixed, ...)
        return(start[setdiff(names(start), names(fixed))])
    }
}


#' Build a self-start model formula with optional fixed parameters
#'
#' Constructs `.x ~ fn(.t, ...)` with each free parameter as a bare
#' symbol and each fixed parameter substituted as its constant value.
#'
#' @param fn Symbol; the self-start model function.
#' @param params Character vector of parameter names in `fn` argument
#'   order.
#' @param fix Named list of fixed parameter values.
#'
#' @returns A two-sided [formula][stats::formula] on `.x` and `.t`.
#'
#' @keywords internal
build_ss_formula <- function(fn, params, fix = list()) {
    args <- lapply(setNames(nm = params), \(.p) fix[[.p]] %||% as.name(.p))
    rhs <- as.call(c(fn, quote(.t), args))
    return(stats::as.formula(call("~", quote(.x), rhs)))
}


#' Combine fitted and fixed coefficients into the full parameter vector
#'
#' @param model An [nls][stats::nls] model of the free parameters.
#' @param params Character vector of parameter names in model order.
#' @param fix Named list of fixed parameter values.
#'
#' @returns A named numeric vector ordered by `params` containing the
#'   fitted coefficients with fixed values merged in.
#'
#' @keywords internal
full_coefs <- function(model, params, fix = list()) {
    coefs <- c(stats::coef(model), unlist(fix))
    return(coefs[intersect(params, names(coefs))])
}


#' Grid search starting estimates for a biexponential fit
#'
#' Deterministic seed for [SSbiexponential()] and [fit_biexp_ratio()]. The
#' biexponential is linear in `A`, `B1` and `B2` once `tau1` and `tau2` are
#' held fixed, so the amplitudes are profiled out by least squares at every
#' point of a log-spaced grid of time constants rather than searched.
#'
#' @details
#' The grid is confined to the ridge `tau2 >= tau_ratio * tau1`. As `tau2`
#' approaches `tau1` the two exponential basis vectors coincide, the design
#' matrix becomes singular, and the amplitudes diverge with cancelling signs;
#' the ridge excludes that degenerate region.
#'
#' Grid limits scale with the span of `t`, so the seed does not depend on the
#' time units of the data. One log-spaced grid serves both components, so the
#' saturating basis `1 - exp(-t / tau)` is built once and its Gram matrix
#' supplies every candidate pair: each pair is then solved in closed form from
#' its own 3x3 sub-block, with no further pass over the data.
#'
#' No sign constraint is placed on `B1` or `B2`. NIRS response amplitudes may
#' legitimately be negative, so identifiability is enforced through the time
#' constants alone. This is not a `direction` mechanism, which steers only the
#' fit window.
#'
#' @param x A numeric vector of the response variable.
#' @param t A numeric vector of the predictor variable (time).
#' @param tau_ratio A numeric lower bound on `tau2 / tau1`.
#' @param TD A numeric time delay subtracted from `t` before fitting.
#' @param n_tau Integer resolution of the time constant grid.
#'
#' @returns A named list of starting estimates (`A`, `B1`, `tau1`, `B2`,
#'   `tau2`) with the attained `rss`, or `NULL` if no grid point is usable.
#'
#' @keywords internal
biexp_grid_start <- function(x, t, tau_ratio = 2.5, TD = 0, n_tau = 64L) {
    span <- diff(range(t))
    if (!is.finite(span) || span <= 0) {
        return(NULL)
    }

    ## response onset is flat before TD, so the delay shifts the time base
    t_fit <- pmax(t - TD, 0)
    taus <- exp(seq(log(span / 100), log(span * 10), length.out = n_tau))
    basis <- vapply(taus, \(.tau) 1 - exp(-t_fit / .tau), numeric(length(t_fit)))

    ## candidate (tau1, tau2) pairs on the ridge
    pairs <- which(outer(taus * tau_ratio, taus, "<="), arr.ind = TRUE)
    if (nrow(pairs) == 0L) {
        return(NULL)
    }
    i <- pairs[, 1L]
    j <- pairs[, 2L]

    ## normal equations for the design [1, -u_i, u_j]; signs match
    ## biexponential(), so the solved amplitudes are on the reported scale
    gram <- crossprod(basis)
    diag_gram <- diag(gram)
    su <- colSums(basis)
    xu <- as.vector(crossprod(basis, x))

    a11 <- length(t_fit)
    a12 <- -su[i]
    a13 <- su[j]
    a22 <- diag_gram[i]
    a23 <- -gram[cbind(i, j)]
    a33 <- diag_gram[j]
    v1 <- sum(x)
    v2 <- -xu[i]
    v3 <- xu[j]

    ## cofactors of the symmetric 3x3 system, evaluated across all pairs at
    ## once; a per-pair solve() would dominate the cost of the search
    c11 <- a22 * a33 - a23^2
    c12 <- a13 * a23 - a12 * a33
    c13 <- a12 * a23 - a13 * a22
    c22 <- a11 * a33 - a13^2
    c23 <- a12 * a13 - a11 * a23
    c33 <- a11 * a22 - a12^2
    det_m <- a11 * c11 + a12 * c12 + a13 * c13

    A <- (c11 * v1 + c12 * v2 + c13 * v3) / det_m
    B1 <- (c12 * v1 + c22 * v2 + c23 * v3) / det_m
    B2 <- (c13 * v1 + c23 * v2 + c33 * v3) / det_m

    ## at the least-squares solution the residual sum of squares reduces to
    ## x'x - b'X'x. a non-positive determinant marks a rank-deficient pair
    rss <- sum(x^2) - (A * v1 + B1 * v2 + B2 * v3)
    rss[!is.finite(rss) | det_m <= 0] <- NA_real_
    if (all(is.na(rss))) {
        return(NULL)
    }

    best <- which.min(rss)
    return(list(
        A = A[[best]],
        B1 = B1[[best]],
        tau1 = taus[[i[[best]]]],
        B2 = B2[[best]],
        tau2 = taus[[j[[best]]]],
        rss = max(rss[[best]], 0)
    ))
}


#' Enforce the requested direction on a converged parametric fit
#'
#' Checks the sign of the fitted amplitude `B - A` against the
#' requested `direction`. Satisfied fits are returned unchanged.
#' Inverted fits are refit with amplitude `D = B - A` sign-bounded via
#' `nls(algorithm = "port")`; sigmoid models divide by `D`, so its
#' magnitude is floored strictly above zero. A refit that fails or
#' pins any coefficient at a sign-floor bound (degenerate flat fit)
#' warns and returns `NULL`. User-fixed parameters in `fix` are held
#' constant in the refit: a fixed `A` or `B` is substituted into the
#' amplitude reparameterisation; when both asymptotes are fixed the
#' amplitude sign is predetermined, so an inverted fit cannot be
#' refit and returns `NULL`.
#'
#' @param model A converged [nls][stats::nls] model object.
#' @param coefs Named numeric coefficient vector including `A`, `B`
#'   (fixed values merged in, e.g. from [full_coefs()]).
#' @param fit_data Data frame with columns `.x` and `.t`.
#' @param direction Character; resolved `"positive"` or `"negative"`.
#' @param amp_fn Symbol; exported model fn taking `(t, A, B, ...)`.
#' @param extra Named numeric start values for remaining free params.
#' @param extra_lower,extra_upper Named numeric bound overrides for
#'   `extra` params. Sign-floor bounds should be data-scaled small
#'   values (not `.Machine$double.eps`) so pinned-floor degeneracy is
#'   detectable.
#' @param fn Symbol or character; the self-start fn named in the
#'   warning.
#' @param .nirs Character; the channel name.
#' @param interval_name Character; the interval label.
#' @param fix Named list of user-fixed parameter values.
#' @inheritParams validate_mnirs
#'
#' @returns A named list `list(model, coefs)` with `coefs` a named
#'   numeric vector in `(A, B, ...)` space including fixed values, or
#'   `NULL` when the direction cannot be satisfied (caller returns
#'   [build_na_results()]).
#'
#' @keywords internal
enforce_direction <- function(
    model,
    coefs,
    fit_data,
    direction,
    amp_fn,
    extra,
    extra_lower = NULL,
    extra_upper = NULL,
    fn,
    .nirs,
    interval_name,
    fix = list(),
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## keep the unconstrained fit when direction is already satisfied
    want <- if (direction == "positive") 1 else -1
    if (sign(coefs[["B"]] - coefs[["A"]]) == want) {
        return(list(model = model, coefs = coefs))
    }

    direction_failed <- function() {
        if (verbose) {
            cli_warn(c(
                "x" = "{.fn {as.character(fn)}} fit for \\
                {.field {(.nirs)}} in {.field {interval_name}} could \\
                not satisfy {.code direction = {.val {direction}}}.",
                "i" = "Returning {.val {NA}} coefficients."
            ), call = warn_call(env))
        }
        return(NULL)
    }

    ## both asymptotes fixed: amplitude sign is predetermined and
    ## contradicts the requested direction; no refit possible
    A_fixed <- fix[["A"]]
    B_fixed <- fix[["B"]]
    if (!is.null(A_fixed) && !is.null(B_fixed)) {
        return(direction_failed())
    }

    ## refit start: amplitude D seeded in the requested direction
    D0 <- want *
        max(abs(coefs[["B"]] - coefs[["A"]]), diff(range(fit_data$.x)) * 0.1)
    D_eps <- diff(range(fit_data$.x)) * 1e-6

    ## build rhs on amplitude D = B - A, substituting any fixed
    ## asymptote: A free `amp_fn(.t, A, A + D, ...)`, B fixed
    ## `amp_fn(.t, B - D, B, ...)`. Remaining params ride as named
    ## args: free as symbols, fixed as constants
    A_expr <- if (is.null(B_fixed)) {
        A_fixed %||% quote(A)
    } else {
        call("-", B_fixed, quote(D))
    }
    B_expr <- B_fixed %||% call("+", A_expr, quote(D))
    tail_args <- c(
        fix[setdiff(names(fix), c("A", "B"))],
        setNames(lapply(names(extra), as.name), names(extra))
    )
    rhs <- as.call(c(amp_fn, quote(.t), A_expr, B_expr, tail_args))
    nls_formula <- stats::as.formula(call("~", quote(.x), rhs))

    ## box bounds: D sign-constrained with strictly positive magnitude
    ## floor (sigmoid models divide by D); extras unbounded unless
    ## overridden. A drops out when either asymptote is fixed
    A_free <- is.null(A_fixed) && is.null(B_fixed)
    start <- c(if (A_free) c(A = coefs[["A"]]), c(D = D0), extra)
    lower <- c(
        if (A_free) c(A = -Inf),
        c(D = if (want > 0) D_eps else -Inf),
        setNames(rep(-Inf, length(extra)), names(extra))
    )
    upper <- c(
        if (A_free) c(A = Inf),
        c(D = if (want > 0) Inf else -D_eps),
        setNames(rep(Inf, length(extra)), names(extra))
    )
    lower[names(extra_lower)] <- extra_lower
    upper[names(extra_upper)] <- extra_upper

    refit <- tryCatch(
        nls(
            nls_formula,
            fit_data,
            start = start,
            lower = lower,
            upper = upper,
            algorithm = "port"
        ),
        error = \(e) NULL
    )

    ## any coefficient pinned at a sign-floor bound (e.g. D, slope,
    ## tau) indicates a degenerate flat fit: no genuine response in
    ## the requested direction
    cf <- if (is.null(refit)) NULL else coef(refit)
    floors <- abs(ifelse(is.finite(lower), lower, upper))
    if (is.null(cf) || any(is.finite(floors) & abs(cf) <= 2 * floors)) {
        return(direction_failed())
    }

    ## back-transform to full (A, B, ...) space including fixed values
    A_val <- if (is.null(B_fixed)) {
        A_fixed %||% cf[["A"]]
    } else {
        B_fixed - cf[["D"]]
    }
    out <- c(
        A = unname(A_val),
        B = unname(A_val + cf[["D"]]),
        cf[setdiff(names(cf), c("A", "D"))],
        unlist(fix[setdiff(names(fix), c("A", "B"))])
    )
    return(list(model = refit, coefs = out))
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
#' @returns A 1-row `data.frame` with columns `n_obs`, `n_params`, `r2`,
#'   `adj_r2`, `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`, and `bic`.
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
        n_params = n_params,
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
        n_params = n_params,
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
            "x" = "Unknown model coefficient{?s}: {.field {invalid}}.",
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
