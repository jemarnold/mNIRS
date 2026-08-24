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
    exponential_drift = "exponential_drift",
    exponential_linear = "exponential_drift",
    exp_drift = "exponential_drift",
    exp_linear = "exponential_drift",
    monoexp_drift = "exponential_drift",
    monoexp_linear = "exponential_drift",
    linear_drift = "exponential_drift",
    drift = "exponential_drift",
    logistic = "sigmoidal",
    gompertz = "sigmoidal",
    xmid = "sigmoidal"
)

#' Detect the direction of a response signal
#'
#' Resolves whether a signal responds upward (`"positive"`) or downward
#' (`"negative"`) by comparing the excursions of `x` above and below its
#' initial baseline, taken as the median of the earliest samples ordered
#' by `t`. The dominant excursion captures the primary response direction
#' even when a fast initial component partially recovers over most of the
#' record (e.g. biexponential drop-recovery), where a net slope would
#' misreport the trend. Used internally to disambiguate peak (maximum)
#' from trough (minimum) detection when `direction = "auto"`.
#'
#' @param fallback A numeric vector (*defaults* to `x`) used to resolve
#'   direction when the excursions above and below baseline tie (e.g.
#'   flat or symmetric data). The absolute maximum and minimum of
#'   `fallback` are compared; if `abs(max) >= abs(min)`, `"positive"` is
#'   returned.
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
    if (direction != "auto") {
        return(direction)
    }
    ## order finite samples by time; baseline = median of the earliest
    ## samples, shrinking to a single sample for short vectors
    ok <- which(is.finite(x) & is.finite(t))
    if (length(ok) == 0L) {
        return("positive")
    }
    x_ord <- x[ok][order(t[ok])]
    n_base <- max(1L, min(5L, length(x_ord) %/% 10L))
    x0 <- median(x_ord[seq_len(n_base)])

    ## dominant excursion from baseline: (max - x0) vs (x0 - min);
    ## fallback abs magnitude comparison breaks ties, itself tied positive
    net <- max(x_ord) + min(x_ord) - 2 * x0
    if (net == 0) {
        net <- abs(max(fallback, na.rm = TRUE)) -
            abs(min(fallback, na.rm = TRUE))
    }
    return(if (net >= 0) "positive" else "negative")
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
#' When `direction = "auto"`, the excursions of `x` above and below its
#' initial baseline (the median of the earliest samples) are compared via
#' [detect_direction()]. If the upward excursion dominates, the function
#' searches for a peak (maximum); if the downward excursion dominates, a
#' trough (minimum). When the excursions tie, the direction is determined
#' by comparing `abs(max(x))` to `abs(min(x))`, with ties defaulting to
#' `"positive"`.
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
    direction <- detect_direction(x, t, direction = direction)
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
#' @returns A named list with: `method`, `model`, `coefficients`, `data`,
#'   `interval_times`, `diagnostics`, `channel_args`, `warnings`, `call`.
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
    ## fallback keeps result dfs built without the attribute working
    warnings <- do.call(rbind, lapply(result_list, attr, "warnings")) %||%
        kinetics_warnings_df()
    rownames(warnings) <- NULL

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
        augmented <- cbind(.df, as.data.frame(fitted_cols, check.names = FALSE))
        
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
            warnings = warnings,
            call = call
        ),
        class = "mnirs_kinetics"
    ))
}


#' Resolve per-interval arguments
#'
#' Peels the interval layer from `worker_args` before per-interval worker
#' dispatch, mirroring the per-channel convention of
#' [resolve_channel_args()]. An argument is treated as an interval map when
#' it is a `list()` with at least one named key matching an interval name,
#' no key matching a channel name (channel maps keep their existing
#' per-channel meaning), and at most one unnamed element acting as the
#' fallback for unlisted intervals. `fix` must additionally be a list of
#' lists, so a plain parameter list (e.g. `fix = list(A = 0)`) stays global.
#'
#' Resolved values may themselves be per-channel maps, which pass untouched
#' to [resolve_channel_args()] downstream. Intervals omitted from a map with
#' no unnamed fallback resolve to `NULL`, falling through to the argument's
#' default, matching omitted-channel behaviour.
#'
#' @param worker_args Named list of method-specific arguments.
#' @param interval_names Character vector of interval names from
#'   [as_data_list()].
#' @param chan_names Character vector of resolved channel names, used only
#'   to give channel keys precedence over interval keys.
#' @inheritParams validate_mnirs
#'
#' @returns A named list with one element per interval; each element is the
#'   `worker_args` list resolved for that interval.
#'
#' @keywords internal
resolve_interval_args <- function(
    worker_args,
    interval_names,
    chan_names,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## `worker_args` comes from `mget()` of a named `unlist()` vector, so
    ## its names attribute can itself carry names, which diverts `vapply()`
    ## result naming; flatten once so name lookups are reliable
    names(worker_args) <- unname(names(worker_args))

    ## classify interval maps; channel keys win so existing per-channel
    ## lists are never reinterpreted, and a plain `fix` parameter list
    ## stays global
    is_map <- vapply(names(worker_args), \(.nm) {
        .a <- worker_args[[.nm]]
        keys <- names(.a) %||% rep("", length(.a))
        is_arg_map(.a) &&
            !any(keys %in% chan_names) &&
            any(keys %in% interval_names) &&
            (.nm != "fix" || all(vapply(.a, is.list, logical(1))))
    }, logical(1))

    ## warn once per interval-map arg: unrecognised interval names are
    ## ignored; omitted intervals with no unnamed fallback fall back to
    ## the argument's default
    if (verbose) {
        lapply(names(worker_args)[is_map], \(.nm) {
            keys <- names(worker_args[[.nm]])
            omitted <- if (all(nzchar(keys))) setdiff(interval_names, keys)
            warn_map_keys(
                .nm,
                unknown = setdiff(keys[nzchar(keys)], interval_names),
                omitted = omitted,
                what = "interval",
                match_hint = "interval names",
                env = env
            )
        })
    }

    ## peel each interval's values: named hit > unnamed fallback > NULL
    return(lapply(setNames(nm = interval_names), \(.int) {
        lapply(setNames(nm = names(worker_args)), \(.nm) {
            .a <- worker_args[[.nm]]
            if (!is_map[[.nm]]) {
                return(.a)
            }
            unnamed <- .a[!nzchar(names(.a))]
            fallback <- if (length(unnamed) > 0L) unnamed[[1L]]
            return(.a[[.int]] %||% fallback)
        })
    }))
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

    ## resolve channel names for interval-map classification only; genuine
    ## channel errors surface later in the worker with proper context
    chan_names <- tryCatch(
        validate_nirs_channels(nirs_quo, data_list[[1L]], env),
        error = \(e) character()
    )
    interval_args <- resolve_interval_args(
        worker_args,
        names(data_list),
        chan_names,
        verbose,
        env
    )

    ## iterate over each interval
    result_list <- lapply(seq_along(data_list), \(.i) {
        rlang::inject(worker(
            data = data_list[[.i]],
            nirs_channels = !!nirs_quo,
            time_channel = !!time_quo,
            !!!interval_args[[.i]],
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
    biexponential = c("use_TD", "tau_mult", "fix"),
    exponential_drift = c("use_TD", "tau_mult", "fix"),
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
#' `fix` is itself a named list, so it is classified before resolution: a
#' plain list of parameter values applies globally to every channel, while a
#' list whose elements are all lists is a per-channel map keyed by channel
#' name. The resolved `fix` is validated per channel against `fix_params`
#' when supplied.
#'
#' @param data A single *"mnirs"* data frame.
#' @param nirs_quo,time_quo Quosures of the worker's `nirs_channels` and
#'   `time_channel` arguments (captured with `enquo()` in the worker frame).
#' @param arg_list Named list of the method's per-channel-capable arguments.
#' @param choices Named list of valid values for choice-type arguments,
#'   passed to [resolve_channel_args()].
#' @param fix_params An *optional* character vector of fixable model
#'   parameter names, or a function of a channel's resolved argument list
#'   returning that vector (for models whose fixable parameters depend on
#'   another argument, e.g. `use_TD`). When supplied, each channel's
#'   resolved `fix` is validated via [validate_fix()].
#' @inheritParams validate_mnirs
#'
#' @returns A named list with `nirs_channels`, `time_channel`, and
#'   `per_channel`.
#'
#' @keywords internal
setup_kinetics_worker <- function(
    data,
    nirs_quo,
    time_quo,
    arg_list,
    choices = list(),
    fix_params = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    validate_mnirs_data(data, env = env)
    nirs_channels <- validate_nirs_channels(nirs_quo, data, env)
    time_channel <- validate_time_channel(time_quo, data, env = env)
    t_vec <- data[[time_channel]]

    ## `fix` is a named list of parameters, which `resolve_channel_args()`
    ## would misread as a channel map. Only a list of lists is per-channel;
    ## a plain parameter list is withheld and broadcast after resolution
    fix <- arg_list$fix
    fix_map <- is.list(fix) &&
        length(fix) > 0L &&
        all(vapply(fix, is.list, logical(1)))
    if (!fix_map) {
        arg_list$fix <- NULL
    }

    ## broadcast global args (applying per-channel list() overrides), then
    ## validate the resolved args once, before fitting any channel
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = arg_list,
        choices = choices,
        verbose = verbose,
        env = env
    )

    if (!fix_map && !is.null(fix)) {
        per_channel <- lapply(per_channel, \(.a) c(.a, list(fix = fix)))
    }
    per_channel <- validate_kinetics_args(
        per_channel,
        data,
        t_vec,
        verbose,
        env = env
    )

    ## validate each channel's fixed parameters against its own model.
    ## `[` assignment keeps an unfixed channel's `fix` key present but
    ## `NULL`, so every channel serialises the same `channel_args` columns
    if (!is.null(fix_params)) {
        per_channel <- lapply(per_channel, \(.a) {
            .p <- if (is.function(fix_params)) fix_params(.a) else fix_params
            f <- validate_fix(.a$fix, .p, env = env)
            .a["fix"] <- list(if (length(f) > 0L) f else NULL)
            .a
        })
    }

    return(list(
        nirs_channels = nirs_channels,
        time_channel = time_channel,
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
#' @param fit_fn A function `(.nirs, x_fit, t_fit, .a, valid)`
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
#'   channel), `"diagnostics"` and `"channel_args"` (data frames, one row
#'   per channel), and `"warnings"` (data frame of conditions captured
#'   during fitting, regardless of `verbose`; zero rows when none fire).
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

    ## collect conditions signalled during fitting; fit-path emitters signal
    ## unconditionally and this single handler governs console emission, so
    ## capture is independent of `verbose`
    warning_rows <- list()
    .nirs_active <- NA_character_
    record <- function(w) {
        warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
            interval = interval_name,
            nirs_channels = .nirs_active,
            type = if (inherits(w, "mnirs_fit_error")) "error" else "warning",
            message = gsub("\n", " - ", cli::ansi_strip(conditionMessage(w)))
        )
    }

    result <- withCallingHandlers(
        {
            ## per-channel fit; collect parallel pieces keyed by channel
            fits <- setNames(
                nm = nirs_channels,
                lapply(nirs_channels, \(.nirs) {
                .nirs_active <<- .nirs
                .a <- per_channel[[.nirs]]

                ## filter for valid finite idx before first extreme + end_window;
                ## data columns and `end_window` are already validated upstream.
                ## fit on time elapsed from onset so extreme detection and
                ## end_window truncation ignore pre-onset baseline (t < 0)
                t_rel <- t_vec - (.a$start_time %||% 0)
                valid <- find_kinetics_idx(
                    data[[.nirs]],
                    t_rel,
                    .a$end_window,
                    .a$direction,
                    bypass_checks = TRUE,
                    env = env
                )
                .a$direction <- valid$direction
                x_fit <- data[[.nirs]][valid$idx]
                t_fit <- t_rel[valid$idx]

                ## method-specific fit; coefs/diag carry method columns only
                fit <- fit_fn(.nirs, x_fit, t_fit, .a, valid)

                ## serialise resolved args: NULL to NA, list() to its deparse(),
                ## dropping internal-only args, so they fit a flat data frame row
                arg_row <- lapply(c(.a, extra_args), \(.x) {
                    if (is.null(.x)) {
                        NA
                    } else if (is.list(.x)) {
                        ## deparse() wraps beyond its default width, which would
                        ## expand the single-row data frame
                        paste(deparse(.x), collapse = "")
                    } else {
                        .x
                    }
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

            ## interval-level conditions from here on
            .nirs_active <- NA_character_

            ## assemble single attributed df (consumed by build_kinetics_results)
            result <- structure(
                do.call(rbind, lapply(fits, `[[`, "coefficients")),
                model = lapply(fits, `[[`, "model"),
                fitted_data = lapply(fits, `[[`, "fitted_data"),
                diagnostics = do.call(rbind, lapply(fits, `[[`, "diagnostics")),
                channel_args = do.call(rbind, lapply(fits, `[[`, "channel_args"))
            )

            ## warn when time coefficients are negative (response before start_time)
            # fmt: skip
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

            result
        },
        warning = \(w) {
            record(w)
            if (!verbose) rlang::cnd_muffle(w)
        },
        message = \(m) if (!verbose) rlang::cnd_muffle(m)
    )

    attr(result, "warnings") <- if (length(warning_rows) == 0L) {
        kinetics_warnings_df()
    } else {
        do.call(rbind, warning_rows)
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
#'   (*without* `nirs_channels`/`time_channel`, which are added upstream),
#'   or a character vector of their column names.
#'
#' @returns A named list with elements `coefs`, `model`, `fitted_data`,
#'   and `diag`.
#'
#' @keywords internal
build_na_results <- function(na_coefs) {
    ## a character vector names the NA columns
    if (is.character(na_coefs)) {
        na_coefs <- as.data.frame(
            setNames(rep(list(NA_real_), length(na_coefs)), na_coefs)
        )
    }
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


#' Warn on a failed or non-converged kinetics model fit
#'
#' Shared warning for the nls-based kinetics workers. Error conditions
#' report a failed fit with an optional hint that the reduced
#' (`n_params - 1`) model is attempted next; warning-class conditions
#' report a fit accepted despite non-convergence.
#'
#' @param fn Symbol or character; the model fn named in the message.
#' @param e The captured condition object.
#' @param .nirs Character; the channel name.
#' @param interval_name Character; the interval label.
#' @param n_params Integer or `NULL`; parameter count prefixed to the fn
#'   name.
#' @param retry Logical; hint that the reduced model fit is attempted next.
#' @inheritParams validate_mnirs
#'
#' @returns `invisible(NULL)`, invoked for its warning side effect.
#'
#' @keywords internal
warn_fit_failed <- function(
    fn,
    e,
    .nirs,
    interval_name,
    n_params = NULL,
    retry = FALSE,
    env = rlang::caller_env()
) {
    ## optional parameter count prefixes the model fn name
    label <- paste0(
        if (!is.null(n_params)) "{n_params}-parameter ",
        "{.fn {as.character(fn)}}"
    )
    where <- " for {.field {(.nirs)}} in {.field {interval_name}}:"
    msg <- if (inherits(e, "warning")) {
        c(
            "!" = paste0(label, " fit warning", where),
            "i" = "{conditionMessage(e)}",
            "i" = "Fit accepted. Consider a reduced \\
            {.fn analyse_kinetics} method."
        )
    } else {
        c(
            "x" = paste0(label, " fit failed", where),
            "!" = "{conditionMessage(e)}",
            if (retry) {
                c("i" = "Attempting {n_params - 1L}-parameter \\
                {.fn {as.character(fn)}} fit.")
            }
        )
    }
    ## classed so the capture handler can distinguish caught fit errors
    cli_warn(
        msg,
        class = if (!inherits(e, "warning")) "mnirs_fit_error",
        call = warn_call(env)
    )
    return(invisible(NULL))
}


#' Fit a self-start model with time-delay fallback
#'
#' Shared attempt skeleton for the nls-based kinetics workers. The TD
#' model is flat at `A` before `TD`, so the pre-onset baseline anchors
#' `A`; the reduced model has no such region and diverges at `t < 0`, so
#' it is fit from `start_time` onward. An under-determined attempt is
#' rejected before it reaches `fitter`, and a failed TD fit falls back to
#' the reduced model without `TD` unless `TD` is user-fixed. Every failure
#' is reported through [warn_fit_failed()].
#'
#' @param x_fit,t_fit Numeric vectors of the channel fit window.
#' @param params Character vector of parameter names in model order,
#'   including `TD` when the channel fits the TD model.
#' @param .a The channel's resolved argument list (`use_TD`, `fix`).
#' @param fitter A function `(.data, .params, on_error)` fitting `.params`
#'   to a data frame of `.x`/`.t` and returning an [nls][stats::nls] model
#'   or `NULL`. `on_error(e)` reports the condition `e` and returns `NULL`,
#'   so it doubles as a [tryCatch()] error handler.
#' @param retry Logical; attempt the reduced model when the TD fit fails.
#' @inheritParams warn_fit_failed
#'
#' @returns A list with `model` (or `NULL`), the `params` actually fit,
#'   the logical row filter `keep`, and the fit `data` frame.
#'
#' @keywords internal
fit_td_fallback <- function(
    x_fit,
    t_fit,
    params,
    .a,
    fitter,
    fn,
    .nirs,
    interval_name,
    env,
    retry = .a$use_TD && !"TD" %in% names(.a$fix)
) {
    attempt <- \(.params, .retry) {
        ## dropping TD narrows the window, so subset per attempt
        keep <- "TD" %in% .params | t_fit >= 0
        data <- data.frame(.x = x_fit[keep], .t = t_fit[keep])
        on_error <- \(e) {
            warn_fit_failed(
                fn,
                e,
                .nirs,
                interval_name,
                length(.params),
                .retry,
                env
            )
            NULL
        }
        n_free <- length(setdiff(.params, names(.a$fix)))
        model <- if (nrow(data) <= n_free) {
            on_error(simpleError(sprintf(
                "%d observation%s for %d free parameters.",
                nrow(data),
                if (nrow(data) == 1L) "" else "s",
                n_free
            )))
        } else {
            fitter(data, .params, on_error)
        }
        list(model = model, params = .params, keep = keep, data = data)
    }
    fit <- attempt(params, retry)
    if (is.null(fit$model) && retry) {
        fit <- attempt(setdiff(params, "TD"), FALSE)
    }
    return(fit)
}


#' Accept or reject a non-converged port fit
#'
#' [stats::nls()] with `algorithm = "port"` and `warnOnly = TRUE` returns
#' a model whose stop certificate failed. It is kept with a warning when
#' `ok` holds and its coefficients are finite; otherwise it is reported as
#' an error and dropped. The port stop code is reported in prose either
#' way.
#'
#' @param model An [nls][stats::nls] model or `NULL`.
#' @param on_error A reporting function; see [fit_td_fallback()].
#' @param ok Logical; a further acceptance condition, e.g. an RSS no worse
#'   than the starting estimates. Evaluated only for a non-converged fit.
#'
#' @returns `model` or `NULL`.
#'
#' @keywords internal
accept_port_fit <- function(model, on_error, ok = TRUE) {
    if (is.null(model) || model$convInfo$isConv) {
        return(model)
    }
    port_msg <- c(
        "7" = "Singular convergence: parameters not individually identifiable.",
        "8" = "False convergence: gradient certificate failed near a non-smooth point.",
        "9" = "Function evaluation limit reached without convergence.",
        "10" = "Iteration limit reached without convergence."
    )
    code <- as.character(model$convInfo$stopCode)
    known <- code %in% names(port_msg)
    msg <- if (known) port_msg[[code]] else model$convInfo$stopMessage
    if (known && ok && all(is.finite(stats::coef(model)))) {
        on_error(simpleWarning(msg))
        return(model)
    }
    return(on_error(simpleError(msg)))
}


#' Assemble a fitted channel result
#'
#' Counterpart of [build_na_results()] for a successful fit: the
#' `coefs`/`model`/`fitted_data`/`diag` list expected by
#' [analyse_kinetics_channels()], with fitted values and diagnostics
#' derived from `model` on the rows in `keep`.
#'
#' @param coefs A 1-row `data.frame` of method coefficients.
#' @param model A fitted model supporting [stats::predict()] and
#'   [stats::coef()].
#' @param valid The [find_kinetics_idx()] result for the channel.
#' @param keep Logical row filter of the fit window used by `model`.
#' @inheritParams fit_td_fallback
#' @inheritParams validate_mnirs
#'
#' @returns A named list with elements `coefs`, `model`, `fitted_data`,
#'   and `diag`.
#'
#' @keywords internal
build_fit_results <- function(
    coefs,
    model,
    x_fit,
    t_fit,
    valid,
    keep = TRUE,
    env = rlang::caller_env()
) {
    fitted_vals <- stats::predict(model)
    return(list(
        coefs = coefs,
        model = model,
        fitted_data = data.frame(
            window_idx = valid$idx[keep],
            fitted = fitted_vals
        ),
        diag = compute_diagnostics(
            x_fit[keep],
            t_fit[keep],
            fitted_vals,
            n_params = length(stats::coef(model)),
            env = env
        )
    ))
}


#' Zero-row kinetics warnings scaffold
#'
#' Stable column template for captured fit conditions, so binding and the
#' returned `warnings` element keep consistent columns when none fire.
#'
#' @returns A zero-row `data.frame` with columns `interval`,
#'   `nirs_channels`, `type`, and `message`.
#'
#' @keywords internal
kinetics_warnings_df <- function() {
    return(data.frame(
        interval = character(),
        nirs_channels = character(),
        type = character(),
        message = character()
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


#' Enforce the requested direction on a converged parametric fit
#'
#' Checks the direction of the fitted curve within the data span,
#' via the dominant excursion from its baseline
#' ([detect_direction()]), against the requested `direction`. The
#' asymptote sign `B - A` is not used directly: a weakly identified
#' slow component can strand the asymptote far beyond the record and
#' invert it relative to the observed response. Satisfied fits are
#' returned unchanged. Inverted fits are refit with amplitude `D = B - A`
#' sign-bounded via `nls(algorithm = "port")`; sigmoid models divide
#' by `D`, so its magnitude is floored strictly above zero. An
#' accepted refit is then re-expressed in the original
#' parameterisation, refit from its own optimum, so the returned
#' model reports consistent coefficient names (`A`, `B`, ...) rather
#' than `D`. A refit that fails, pins a sign-floored coefficient
#' (degenerate flat fit), or loses the requested sign on
#' re-expression warns and returns `NULL`. User-fixed parameters in
#' `fix` are held constant in the refit: a fixed `A` or `B` is
#' substituted into the amplitude reparameterisation; when both
#' asymptotes are fixed the amplitude sign is predetermined, so an
#' inverted fit cannot be refit and returns `NULL`.
#'
#' @param model A converged [nls][stats::nls] model object.
#' @param coefs Named numeric coefficient vector including `A` and
#'   the `B_name` asymptote (fixed values merged in, e.g. from
#'   [full_coefs()]).
#' @param fit_data Data frame with columns `.x` and `.t`.
#' @param direction Character; resolved `"positive"` or `"negative"`.
#' @param amp_fn Symbol; exported model fn taking `t`, `A`, and the
#'   `B_name` asymptote as named arguments.
#' @param B_name Character; name of the asymptote coefficient whose
#'   difference from `A` defines the directed amplitude (*default*
#'   `"B"`; `"B2"` for the biexponential).
#' @param extra Named numeric start values for remaining free params.
#' @param mirror Character; names of `extra` params that are
#'   intermediate asymptotes (e.g. the biexponential `B1`), seeded on
#'   the requested side of `A` so the refit does not start inverted.
#' @param extra_lower,extra_upper Named numeric bound overrides for
#'   `extra` params. Sign-floor bounds should be data-scaled small
#'   values (not `.Machine$double.eps`) so pinned-floor degeneracy is
#'   detectable.
#' @param floor_params Character; names of refit coefficients subject
#'   to the pinned-floor degeneracy check. `NULL` (*default*) checks
#'   every finite bound; restrict when other bounds are structural
#'   (e.g. the biexponential time-constant bounds).
#' @param fn Symbol or character; the self-start fn named in the
#'   warning.
#' @param .nirs Character; the channel name.
#' @param interval_name Character; the interval label.
#' @param fix Named list of user-fixed parameter values; values may
#'   be language objects expressed in other free parameters.
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
    B_name = "B",
    mirror = NULL,
    extra_lower = NULL,
    extra_upper = NULL,
    floor_params = NULL,
    fn,
    .nirs,
    interval_name,
    fix = list(),
    env = rlang::caller_env()
) {
    ## keep the unconstrained fit when direction is already satisfied.
    ## direction is judged on the fitted curve within the data span, so
    ## an asymptote stranded beyond the record cannot invert the check
    want <- if (direction == "positive") 1 else -1
    if (detect_direction(stats::fitted(model), fit_data$.t) == direction) {
        return(list(model = model, coefs = coefs))
    }

    direction_failed <- function() {
        cli_warn(c(
            "x" = "{.fn {as.character(fn)}} fit for \\
            {.field {(.nirs)}} in {.field {interval_name}} could \\
            not satisfy {.code direction = {.val {direction}}}.",
            "i" = "Returning {.val {NA}} coefficients."
        ), call = warn_call(env))
        return(NULL)
    }

    ## shared bounded refit; params unbounded except named overrides.
    ## the refit may start far from its constrained optimum, so it gets
    ## the same iteration budget as the primary fit
    port_fit <- function(rhs, start, lower = NULL, upper = NULL) {
        lwr <- setNames(rep(-Inf, length(start)), names(start))
        upr <- -lwr
        lwr[names(lower)] <- lower
        upr[names(upper)] <- upper
        ## non-smooth models (biexponential onset kink) can stall PORT
        ## with false convergence; a finite-coefficient stall is accepted
        ## as in the primary fit
        model <- tryCatch(
            suppressWarnings(nls(
                stats::as.formula(call("~", quote(.x), rhs)),
                fit_data,
                start = start,
                lower = lwr,
                upper = upr,
                algorithm = "port",
                control = stats::nls.control(maxiter = 500L, warnOnly = TRUE)
            )),
            error = \(e) NULL
        )
        return(accept_port_fit(model, \(e) NULL))
    }

    ## both asymptotes fixed: amplitude sign is predetermined and
    ## contradicts the requested direction; no refit possible
    A_fixed <- fix[["A"]]
    B_fixed <- fix[[B_name]]
    if (!is.null(A_fixed) && !is.null(B_fixed)) {
        return(direction_failed())
    }

    ## refit start: amplitude D seeded in the requested direction. a
    ## weakly identified fit can strand the asymptote far beyond the
    ## data, so the seed magnitude is confined to the observed scale
    x_span <- diff(range(fit_data$.x))
    D0 <- want *
        min(max(abs(coefs[[B_name]] - coefs[["A"]]), x_span * 0.1), x_span)
    D_eps <- x_span * 1e-6

    ## intermediate asymptotes seeded on the requested side of A; an
    ## inverted seed drags the refit into a degenerate tau-floor basin
    m <- intersect(mirror, names(extra))
    extra[m] <- coefs[["A"]] + want * abs(extra[m] - coefs[["A"]])

    ## build rhs on amplitude D = B - A, substituting any fixed
    ## asymptote: A free `amp_fn(.t, A, A + D, ...)`, B fixed
    ## `amp_fn(.t, B - D, B, ...)`. Remaining params ride as named
    ## args: free as symbols, fixed as constants or expressions
    A_expr <- if (is.null(B_fixed)) {
        A_fixed %||% quote(A)
    } else {
        call("-", B_fixed, quote(D))
    }
    B_expr <- B_fixed %||% call("+", A_expr, quote(D))
    tail_args <- c(
        fix[setdiff(names(fix), c("A", B_name))],
        setNames(lapply(names(extra), as.name), names(extra))
    )
    rhs <- as.call(c(
        amp_fn,
        quote(.t),
        setNames(list(A_expr, B_expr), c("A", B_name)),
        tail_args
    ))

    ## box bounds: D sign-constrained with strictly positive magnitude
    ## floor (sigmoid models divide by D); extras unbounded unless
    ## overridden. A drops out when either asymptote is fixed
    A_free <- is.null(A_fixed) && is.null(B_fixed)
    start <- c(if (A_free) c(A = coefs[["A"]]), c(D = D0), extra)
    d_bnd <- c(D = want * D_eps)
    lower <- c(if (want > 0) d_bnd, extra_lower)
    upper <- c(if (want < 0) d_bnd, extra_upper)

    refit <- port_fit(rhs, start, lower, upper)
    if (is.null(refit)) {
        return(direction_failed())
    }

    ## a coefficient pinned at a sign-floor bound (e.g. D, slope,
    ## tau) indicates a degenerate flat fit: no genuine response in
    ## the requested direction. finite lower bounds take precedence;
    ## `floor_params` restricts the check where other bounds are
    ## structural rather than sign floors
    cf <- coef(refit)
    bnd <- c(lower, upper)
    bnd <- bnd[is.finite(bnd)]
    chk <- intersect(floor_params %||% names(cf), names(bnd))
    if (any(abs(cf[chk]) <= 2 * abs(bnd[chk]))) {
        return(direction_failed())
    }

    ## back-transform to full (A, B, ...) space including fixed values
    A_val <- if (is.null(B_fixed)) {
        A_fixed %||% cf[["A"]]
    } else {
        B_fixed - cf[["D"]]
    }
    B_val <- A_val + cf[["D"]]

    ## re-express in the original parameterisation so the returned
    ## model reports consistent coefficient names rather than D. the
    ## accepted refit is an interior local minimum, so a bounded
    ## refit started there converges in place
    rhs2 <- as.call(c(
        amp_fn,
        quote(.t),
        setNames(
            list(A_fixed %||% quote(A), B_fixed %||% as.name(B_name)),
            c("A", B_name)
        ),
        tail_args
    ))
    start2 <- c(
        if (is.null(A_fixed)) c(A = unname(A_val)),
        if (is.null(B_fixed)) setNames(unname(B_val), B_name),
        cf[setdiff(names(cf), c("A", "D"))]
    )
    refit2 <- port_fit(rhs2, start2, extra_lower, extra_upper)
    if (is.null(refit2)) {
        return(direction_failed())
    }

    ## merge fixed values back into the reported coefficients; a
    ## language fix value is an expression in a free parameter and
    ## has no constant to merge. re-expression must hold the
    ## requested within-span direction, else it cannot be satisfied
    fix_num <- fix[vapply(fix, is.numeric, logical(1))]
    cf2 <- coef(refit2)
    out <- c(cf2, unlist(fix_num[setdiff(names(fix_num), names(cf2))]))
    if (detect_direction(stats::fitted(refit2), fit_data$.t) != direction) {
        return(direction_failed())
    }
    return(list(model = refit2, coefs = out))
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
#'   `summary(lm)$r.squared`); a bounded `[0, 1]` pseudo-`R^2` for non-linear
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
        cli_warn(c(
            "!" = "{.arg x}, {.arg t}, and {.arg fitted} must be \\
            {.cls numeric} vectors of equal lengths to return model \\
            diagnostics."
        ), call = warn_call(env))
        return(return_na)
    }

    ## residuals and sums of squares (reused throughout)
    x_mean <- mean(x)
    resid <- x - fitted
    ss_res <- sum(resid^2)
    ss_tot <- sum((x - x_mean)^2)
    ss_fit <- sum((fitted - mean(fitted))^2)

    ## R^2: squared Pearson correlation; equals 1 - SSres/SStot for OLS,
    ## bounded [0, 1] pseudo-R^2 for non-linear fits
    r2 <- if (ss_tot == 0 || ss_fit == 0) {
        NA_real_
    } else {
        stats::cor(x, fitted)^2
    }

    ## adjusted R^2: penalised by n_params; valid for OLS linear models
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
