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
    exp_lin = "exponential_drift",
    monoexp_drift = "exponential_drift",
    monoexp_linear = "exponential_drift",
    linear_drift = "exponential_drift",
    drift = "exponential_drift",
    logistic = "sigmoidal",
    gompertz = "sigmoidal",
    xmid = "sigmoidal"
)


## canonical method -> method-specific worker arg names.
kinetics_dispatch <- list(
    common = c("start_time", "direction", "end_window"),
    response_time = c("fraction"),
    peak_slope = c("width", "span", "align", "partial", "na.rm"),
    monoexponential = c("use_TD", "fix"),
    biexponential = c("use_TD", "fix", "tau1_flex", "TD_flex", "A_flex"),
    exponential_drift = c("use_TD", "tau_mult", "fix"),
    sigmoidal = c("shape", "fix")
)


## coefficients reported as time points elapsed from `start_time`
# fmt: skip
kinetics_time_coefs <- c(
    "response_time", "peak_slope_time", "TD", "MRT", "HRT", "texc", "xmid"
)


## canonical method -> interval worker name, resolved at call time so
## file collation order is irrelevant
kinetics_workers <- c(
    response_time = "analyse_response_time",
    peak_slope = "analyse_peak_slope",
    monoexponential = "analyse_monoexponential",
    biexponential = "analyse_biexponential",
    exponential_drift = "analyse_exponential_drift",
    sigmoidal = "analyse_logistic"
)


## nested reductions: a full method and the reduced method it collapses
## to when the extra parameters are not supported by the data (see
## `reduce_kinetics()`). `coef_map` copies reduced coefficient columns
## into the full schema (full column <- reduced column); unmapped full
## columns are NA. `fix_map` carries user-fixed parameters over to the
## reduced fit (full parameter -> reduced parameter); unmapped are
## dropped. `args` overrides reduced worker arguments so the comparator
## fits the same window as the full model. `accept` is a vectorised
## shape test on the full coefficient data frame (NULL for the F-test
## alone), explained by `reject_msg`
kinetics_reductions <- list(
    biexponential = list(
        to = "monoexponential",
        coef_map = c(
            A = "A", B1 = "B", B2 = "B", tau1 = "tau", TD = "TD", MRT = "MRT",
            MRT_fitted = "MRT_fitted"
        ),
        fix_map = c(A = "A", B2 = "B", TD = "TD"),
        ## the biexponential `end_window` bounds its fast phase only; the
        ## full model spans the whole response
        args = list(end_window = Inf),
        ## excursion-recovery only: the fitted turning point must exist
        accept = \(cf) is.finite(cf$texc),
        reject_msg = "Fitted response is monotonic ({.field texc} is \\
        {.val {NA}})."
    )
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
    ## bind channel-level attr dfs across intervals; each carries `interval`
    ## col[1]. result dfs built without `warnings` fall back to zero rows
    bind_attr <- function(attr_name) {
        df <- do.call(rbind, lapply(result_list, attr, attr_name)) %||%
            kinetics_warnings_df()
        rownames(df) <- NULL
        return(df)
    }
    channel_args <- bind_attr("channel_args")

    ## combine scalar coefficients; `interval` to col[1] and the resolved
    ## fit onset `start_time` (matched from `channel_args` by interval and
    ## channel, since coefs may hold multiple rows per channel) after
    ## `nirs_channels`
    coefs <- do.call(rbind, result_list)
    key <- \(.df) paste(.df$interval, .df$nirs_channels)
    coefs$start_time <- channel_args$start_time[match(key(coefs), key(channel_args))]
    lead_cols <- c("interval", "nirs_channels", "start_time")
    coefs <- coefs[, c(lead_cols, setdiff(names(coefs), lead_cols))]
    rownames(coefs) <- NULL

    ## start_times: resolved fit onset (user value > interval_times metadata
    ## > first time), uniform across channels within an interval;
    ## end_times sourced from interval_times metadata when any is present.
    ## recursive inputs repeat interval names across source channels, so
    ## dedupe to one row per interval
    first <- !duplicated(names(data_list))
    interval_times_df <- data.frame(
        interval = names(data_list)[first],
        start_times = coefs$start_time[!duplicated(coefs$interval)]
    )
    end_times <- vapply(data_list[first], \(.df) {
        unlist(attr(.df, "interval_times"))[2L] %||% NA_real_
    }, numeric(1), USE.NAMES = FALSE)
    if (!all(is.na(end_times))) {
        interval_times_df$end_times <- end_times
    }

    ## augment `data_list` dfs with `<nirs_channels>_fitted` columns,
    ## carrying the original metadata forward
    fitted_data_list <- Map(\(.df, .result) {
        fitted_data <- attr(.result, "fitted_data")
        fitted_cols <- lapply(fitted_data, \(.pred) {
            replace(rep(NA_real_, nrow(.df)), .pred$window_idx, .pred$fitted)
        })
        names(fitted_cols) <- paste0(names(fitted_data), "_fitted")
        metadata <- attributes(.df)
        metadata$nirs_channels <- unique(.result$nirs_channels)
        metadata$time_channel <- attr(.result, "time_channel")
        create_mnirs_data(
            cbind(.df, as.data.frame(fitted_cols, check.names = FALSE)),
            metadata
        )
    }, data_list, result_list)
    ## add `class = "mnirs"` for `plot.mnirs`
    class(fitted_data_list) <- c("mnirs", class(fitted_data_list))

    ## normalise call: function name to generic, method to canonical form
    call[[1L]] <- quote(analyse_kinetics)
    if ("method" %in% names(call)) {
        call$method <- method
    }

    ## per-interval model lists keyed by channel; intervals repeated across
    ## source channels (recursive inputs) concatenate into one list
    model_list <- setNames(lapply(result_list, attr, "model"), names(data_list))
    models <- lapply(
        split(model_list, factor(names(data_list), unique(names(data_list)))),
        \(.m) do.call(c, unname(.m))
    )

    return(structure(
        list(
            method = method,
            model = models,
            coefficients = coefs,
            data = fitted_data_list,
            interval_times = interval_times_df,
            diagnostics = bind_attr("diagnostics"),
            channel_args = channel_args,
            warnings = bind_attr("warnings"),
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


#' Split interval data frames into sample groups
#'
#' Applies the `group_intervals` argument of [analyse_kinetics()] to a
#' named list of data frames. `"ensemble"` returns `data_list` unchanged;
#' a `list()` of sample (row) indices subsets every data frame into one
#' interval per group. Samples in no group are dropped and samples in
#' several groups are warned about. Group names become interval names
#' (`interval_<n>` when unnamed), suffixed `<group>_<df>` when
#' `data_list` holds more than one data frame. Row-subset intervals no
#' longer correspond to their `interval_times`/`interval_span` metadata,
#' so those attributes are dropped.
#'
#' @param data_list Named list of data frames from [as_data_list()].
#' @param group_intervals `"ensemble"` or a `list()` of integer-valued
#'   sample index vectors; see [analyse_kinetics()].
#' @inheritParams validate_mnirs
#'
#' @returns A named list of data frames, one per group per data frame.
#'
#' @keywords internal
split_kinetics_groups <- function(
    data_list,
    group_intervals,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    if (identical(group_intervals, "ensemble")) {
        return(data_list)
    }
    if (!is.list(group_intervals)) {
        cli_abort(c(
            "x" = "{.arg group_intervals} must be {.val ensemble} or a \\
            {.cls list} of sample indices.",
            "i" = "e.g. {.code group_intervals = list(trial1 = 1:10, \\
            trial2 = 11:20)}."
        ), call = env)
    }
    ## the same groups apply to every df, so validate against the shortest
    n_samples <- min(vapply(data_list, nrow, integer(1)))
    validate_interval_groups(group_intervals, n_samples, env)

    ids <- unlist(group_intervals, use.names = FALSE)
    nms <- names(group_intervals) %||% character(length(group_intervals))
    unnamed <- !nzchar(nms)
    nms[unnamed] <- paste0("interval_", which(unnamed))

    if (verbose) {
        n_drop <- n_samples - length(unique(ids))
        if (n_drop > 0L) {
            cli_inform(c(
                "i" = "{n_drop} sample{?s} not specified by \\
                {.arg group_intervals} excluded from analysis."
            ), call = env)
        }
        dupes <- unique(ids[duplicated(ids)])
        if (length(dupes) > 0L) {
            cli_warn(c(
                "!" = "Duplicates detected of {qty(length(dupes))} \\
                sample{?s} {.field {dupes}} across {.arg group_intervals}.",
                "i" = "Re-specify {.arg group_intervals} to remove duplicates."
            ), call = warn_call(env))
        }
    }

    ## row-subset each df per group; interval boundary metadata no longer
    ## describes the subset, so only channel/device metadata carries over
    single <- length(data_list) == 1L
    out <- lapply(names(data_list), \(.nm) {
        .df <- data_list[[.nm]]
        attrs <- attributes(.df)
        metadata <- attrs[intersect(mnirs_metadata, names(attrs))]
        groups <- lapply(group_intervals, \(.g) {
            structure(
                create_mnirs_data(.df[.g, , drop = FALSE], metadata),
                interval_times = NULL,
                interval_span = NULL
            )
        })
        setNames(groups, if (single) nms else paste0(nms, "_", .nm))
    })
    return(unlist(out, recursive = FALSE))
}


#' Run a kinetics worker over each interval and collate results
#'
#' Shared skeleton for `analyse_kinetics.*` methods: normalises `data`
#' to a named list of interval data frames, splits sample groups via
#' [split_kinetics_groups()], calls the method worker once per interval,
#' and collates results via [build_kinetics_results()].
#'
#' @param data A data frame, list of data frames, or grouped data frame.
#' @param worker Function; the interval-level worker, e.g.
#'   [analyse_monoexponential()].
#' @param method Character; the canonical method name.
#' @param worker_args Named list of method-specific arguments passed
#'   to `worker`.
#' @param nirs_quo,time_quo Quosures of the caller's `nirs_channels`
#'   and `time_channel` arguments, captured in the method frame.
#' @param group_intervals `"ensemble"` or a `list()` of sample index
#'   vectors; see [split_kinetics_groups()].
#' @param zero_time Logical; if `TRUE`, rebases each interval's
#'   `time_channel` to start from zero, shifting `interval_times` metadata
#'   by the same offset.
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
    group_intervals,
    zero_time,
    verbose,
    call,
    env
) {
    ## normalise input to named list of data frames
    data_list <- as_data_list(data, env = env)

    ## recursive coef input: time-point coefs are elapsed from each interval's
    ## onset, so shift the chosen `time_channel` to absolute time
    recursive <- inherits(data, "mnirs_kinetics")
    if (recursive) {
        ## source channels of the input coefs, qualifying the fitted channel
        ## names after the worker loop below
        src_channels <- names(data_list)
        tc <- validate_time_channel(time_quo, data_list[[1L]], env = env)
        if (tc %in% kinetics_time_coefs) {
            data_list <- lapply(data_list, \(.df) {
                .df[[tc]] <- .df[[tc]] + .df$start_time
                .df
            })
            if (verbose) {
                cli_inform(c(
                    "i" = "{.arg time_channel} = {.val {tc}} converted to \\
                    absolute time from {.field start_time}."
                ), call = env)
            }
        }
    }

    ## split sample groups after the time conversion so group names key
    ## the per-interval args below
    data_list <- split_kinetics_groups(data_list, group_intervals, verbose, env)

    ## rebase each interval to its first time sample; interval_times metadata
    ## (vector, or list for ensembles) shifts by the same offset so the
    ## `start_time` fallback stays aligned
    if (zero_time) {
        tc <- validate_time_channel(time_quo, data_list[[1L]], env = env)
        data_list <- lapply(data_list, \(.df) {
            t0 <- min(.df[[tc]], na.rm = TRUE)
            it <- attr(.df, "interval_times")
            attr(.df, "interval_times") <- if (is.list(it)) {
                lapply(it, `-`, t0)
            } else if (!is.null(it)) {
                it - t0
            }
            zero_offset_data(.df, tc, t0)
        })
    }

    ## resolve channel names for interval-map classification only; genuine
    ## channel errors surface later in the worker with proper context
    chan_names <- tryCatch(
        validate_nirs_channels(nirs_quo, data_list[[1L]], env),
        error = \(e) character()
    )
    # fmt: skip
    interval_args <- resolve_interval_args(
        worker_args, names(data_list), chan_names, verbose, env
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

    ## recursive coef input: move the source channel qualifier from the
    ## interval names (`trial1_hhb`) to the fitted channel names
    ## (`hhb_slope`), so intervals group across source channels
    if (recursive) {
        src <- rep(src_channels, each = length(data_list) %/% length(src_channels))
        intervals <- mapply(\(.nm, .src) {
            sfx <- paste0("_", .src)
            if (endsWith(.nm, sfx)) substr(.nm, 1L, nchar(.nm) - nchar(sfx)) else .nm
        }, names(data_list), src, USE.NAMES = FALSE)
        relabelled <- Map(\(.res, .df, .src, .int) {
            raw <- names(attr(.res, "fitted_data"))
            pref <- paste(.src, raw, sep = "_")
            requalify <- \(.d) {
                .d$interval <- rep(.int, nrow(.d))
                .d$nirs_channels <- pref[match(.d$nirs_channels, raw)]
                .d
            }
            .res <- requalify(.res)
            attr(.res, "diagnostics") <- requalify(attr(.res, "diagnostics"))
            attr(.res, "channel_args") <- requalify(attr(.res, "channel_args"))
            attr(.res, "warnings") <- requalify(attr(.res, "warnings"))
            names(attr(.res, "fitted_data")) <- pref
            names(attr(.res, "model")) <- pref
            names(.df)[match(raw, names(.df))] <- pref
            list(res = .res, df = .df)
        }, result_list, data_list, src, intervals)
        result_list <- setNames(lapply(relabelled, `[[`, "res"), intervals)
        data_list <- setNames(lapply(relabelled, `[[`, "df"), intervals)
    }

    ## collate and return mnirs_kinetics object
    return(build_kinetics_results(data_list, result_list, method, call))
}


#' Run a kinetics method with its nested reduced-model fallback
#'
#' Wraps [analyse_kinetics_intervals()] for methods listed in
#' `kinetics_reductions`: the full model is fit first, then the reduced
#' comparator on the same intervals with arguments carried over by
#' [reduced_worker_args()], and each channel fit is resolved by
#' [reduce_kinetics()]. Methods without a reduction, or `reduce = FALSE`,
#' return the full fit with the `model` coefficient column only.
#'
#' @param method Character; the canonical full method name.
#' @param reduce Logical; attempt the reduction.
#' @inheritParams analyse_kinetics_intervals
#'
#' @returns An *"mnirs_kinetics"* object.
#'
#' @keywords internal
analyse_kinetics_reduced <- function(
    data,
    method,
    worker_args,
    nirs_quo,
    time_quo,
    group_intervals,
    zero_time,
    verbose,
    call,
    env,
    reduce = TRUE
) {
    full <- analyse_kinetics_intervals(
        data,
        get(kinetics_workers[[method]], mode = "function"),
        method,
        worker_args,
        nirs_quo,
        time_quo,
        group_intervals,
        zero_time,
        verbose,
        call,
        env
    )
    spec <- kinetics_reductions[[method]]
    if (!isTRUE(reduce) || is.null(spec)) {
        full$coefficients <- add_model_col(full$coefficients, method)
        return(full)
    }
    ## comparator conditions are captured in its `warnings`, not emitted
    reduced <- analyse_kinetics_intervals(
        data,
        get(kinetics_workers[[spec$to]], mode = "function"),
        spec$to,
        reduced_worker_args(full, worker_args, spec),
        nirs_quo,
        time_quo,
        group_intervals,
        zero_time,
        verbose = FALSE,
        call,
        env
    )
    return(reduce_kinetics(full, reduced, spec, verbose, env))
}


## `model` column after `start_time`: the method each row's coefficients
## come from
add_model_col <- function(coefs, model) {
    coefs$model <- model
    lead <- c("interval", "nirs_channels", "start_time", "model")
    return(coefs[, c(lead, setdiff(names(coefs), lead))])
}


#' Worker arguments for a reduced comparator fit
#'
#' Derives the argument list of the reduced method in `spec` from the
#' full method's arguments and result, so both fits share a window and
#' time-delay structure and the extra-sum-of-squares F-test in
#' [reduce_kinetics()] is nested. Arguments the reduced method does not
#' take are dropped. `use_TD` mirrors the model each channel actually
#' resolved to (`TD` finite), as a per-interval map of per-channel
#' values; a channel whose full fit failed keeps the user's setting.
#' `fix` is carried over through `spec$fix_map`, recursing into
#' per-interval and per-channel maps. `spec$args` overrides arguments
#' outright (e.g. the comparator window).
#'
#' @param full The full method's *"mnirs_kinetics"* result.
#' @param worker_args Named list of the full method's arguments.
#' @param spec A `kinetics_reductions` entry.
#'
#' @returns A named list of worker arguments for the reduced method.
#'
#' @keywords internal
reduced_worker_args <- function(full, worker_args, spec) {
    ## `mget()` of a named vector nests names; flatten for lookups
    names(worker_args) <- unname(names(worker_args))
    keep <- unlist(kinetics_dispatch[c("common", spec$to)], use.names = FALSE)
    args <- worker_args[intersect(names(worker_args), keep)]
    args[names(spec$args)] <- spec$args

    coefs <- full$coefficients
    intervals <- names(full$data)
    chans <- unique(coefs$nirs_channels)

    if ("use_TD" %in% keep) {
        first_coef <- setdiff(
            names(coefs),
            c("interval", "nirs_channels", "start_time")
        )[[1L]]
        td <- ifelse(
            is.na(coefs[[first_coef]]),
            full$channel_args$use_TD %||% TRUE,
            is.finite(coefs$TD)
        )
        by_int <- split(setNames(td, coefs$nirs_channels), coefs$interval)
        by_int <- by_int[intervals]
        uniform <- all(vapply(by_int, identical, logical(1), by_int[[1L]]))
        ## an interval named after a channel would be read as a channel
        ## map; fall back to a per-channel map, else the user's setting.
        ## recursive results relabel interval/channel names after fitting,
        ## so maps cannot key the reduced run (duplicated intervals):
        ## mirror only a uniform resolution
        args$use_TD <- if (anyDuplicated(intervals) > 0L) {
            if (length(unique(td)) == 1L) td[[1L]] else worker_args$use_TD
        } else if (!any(intervals %in% chans)) {
            lapply(by_int, as.list)
        } else if (uniform) {
            as.list(by_int[[1L]])
        } else {
            worker_args$use_TD
        }
    }

    if ("fix" %in% keep && !is.null(args$fix)) {
        fix <- map_fix(args$fix, spec$fix_map)
        args$fix <- if (length(fix) > 0L) fix
    }
    return(args)
}


## rename a fixed-parameter list through `fix_map`, recursing into maps
## (lists of lists); an emptied leaf stays `list()` so nested maps keep
## their keys
map_fix <- function(fix, fix_map) {
    if (length(fix) > 0L && all(vapply(fix, is.list, logical(1)))) {
        return(lapply(fix, map_fix, fix_map))
    }
    fix <- fix[intersect(names(fix), names(fix_map))]
    return(setNames(fix, fix_map[names(fix)]))
}


#' Resolve a full fit against its nested reduced comparator
#'
#' Row-wise (per interval and channel) model selection between the full
#' method's result and the reduced method's result fit by
#' [analyse_kinetics_reduced()]. The full fit is kept when it succeeded,
#' passes the shape test `spec$accept`, and the extra-sum-of-squares
#' F-test against the reduced fit rejects the reduction at `p < 0.05`.
#' The test needs both fits on the same observations with spare degrees
#' of freedom; otherwise the shape test alone decides. A failed full fit
#' falls back to the reduced fit; a row where both failed is left as is.
#'
#' Reduced rows carry the reduced coefficients, model, fitted values,
#' and diagnostics in the full method's schema via `spec$coef_map`, and
#' `model` names the method each row comes from. The reduction is
#' warned about and recorded in `warnings` along with any conditions the
#' reduced fit raised for that channel.
#'
#' @param full,reduced The *"mnirs_kinetics"* results of the full and
#'   reduced methods on the same data.
#' @param spec A `kinetics_reductions` entry.
#' @inheritParams validate_mnirs
#'
#' @returns `full`, resolved.
#'
#' @keywords internal
reduce_kinetics <- function(
    full,
    reduced,
    spec,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    alpha <- 0.05
    cf <- full$coefficients
    cr <- reduced$coefficients
    key <- \(.d) paste(.d$interval, .d$nirs_channels)
    if (!identical(key(cf), key(cr))) {
        cli_abort(c(
            "x" = "Internal: reduced fit rows do not align with the full fit."
        ), call = env)
    }
    id_cols <- c("interval", "nirs_channels", "start_time")
    coef_cols <- setdiff(names(cf), id_cols)
    full_ok <- !is.na(cf[[coef_cols[[1L]]]])
    red_ok <- !is.na(cr[[setdiff(names(cr), id_cols)[[1L]]]])
    shape_ok <- if (is.null(spec$accept)) {
        rep(TRUE, nrow(cf))
    } else {
        spec$accept(cf) %in% TRUE
    }

    ## extra-sum-of-squares F-test from the diagnostics (rss = rmse^2 * n)
    df_f <- full$diagnostics
    df_r <- reduced$diagnostics
    df1 <- df_f$n_params - df_r$n_params
    df2 <- df_f$n_obs - df_f$n_params
    comparable <- full_ok &
        red_ok &
        df_f$n_obs == df_r$n_obs &
        df1 >= 1L &
        df2 >= 1L
    comparable[is.na(comparable)] <- FALSE
    rss_f <- df_f$rmse^2 * df_f$n_obs
    rss_r <- df_r$rmse^2 * df_r$n_obs
    f_stat <- ((rss_r - rss_f) / df1) / (rss_f / df2)
    p <- rep(NA_real_, nrow(cf))
    p[comparable] <- ifelse(
        rss_f[comparable] == 0,
        0,
        stats::pf(
            f_stat[comparable],
            df1[comparable],
            df2[comparable],
            lower.tail = FALSE
        )
    )

    keep_full <- (full_ok & shape_ok & ifelse(comparable, p < alpha, TRUE)) %in%
        TRUE
    reduce <- !keep_full & (full_ok | red_ok)
    idx <- which(reduce)
    cf$model <- ifelse(reduce, spec$to, full$method)

    if (length(idx) > 0L) {
        ## coefficients, diagnostics: reduced values in the full schema
        cf[idx, coef_cols] <- NA
        cf[idx, names(spec$coef_map)] <- cr[idx, unname(spec$coef_map)]
        full$diagnostics[idx, ] <- df_r[idx, ]

        ## model, fitted values: per-interval swaps of the reduced channels;
        ## `[<-` keeps the list classes and data frame metadata
        swap <- \(.int, .x, .y, .cols) {
            ch <- cf$nirs_channels[idx][cf$interval[idx] == .int]
            ## recursive results hold one source channel per data frame, so
            ## swap only the columns each element declares
            cols <- intersect(.cols(ch), names(.x))
            if (length(cols) > 0L) {
                .x[cols] <- .y[cols]
            }
            .x
        }
        full$model[] <- Map(
            swap, names(full$model), full$model, reduced$model,
            MoreArgs = list(.cols = identity)
        )
        full$data[] <- Map(
            swap, names(full$data), full$data, reduced$data,
            MoreArgs = list(.cols = \(.ch) paste0(.ch, "_fitted"))
        )

        ## warn per reduced row; recorded alongside the reduced fit's own
        ## conditions for those channels
        fn_full <- paste0("SS", full$method)
        fn_red <- paste0("SS", spec$to)
        reasons <- vapply(idx, \(.i) {
            if (!full_ok[[.i]]) {
                "Fit failed."
            } else if (!shape_ok[[.i]]) {
                cli::format_inline(spec$reject_msg)
            } else {
                cli::format_inline(
                    "Extra-sum-of-squares F({df1[.i]}, {df2[.i]}) = \\
                    {signif(f_stat[.i], 3)}, p = {signif(p[.i], 3)} >= \\
                    {alpha}."
                )
            }
        }, character(1))
        heads <- vapply(idx, \(.i) {
            cli::format_inline(
                "{.fn {fn_full}} fit for {.field {cf$nirs_channels[.i]}} \\
                in {.field {cf$interval[.i]}} reduced to {.fn {fn_red}}."
            )
        }, character(1))
        if (verbose) {
            Map(\(.h, .r) {
                cli_warn(c("!" = .h, "i" = .r), call = warn_call(env))
            }, heads, reasons)
        }
        rw <- reduced$warnings
        rw <- rw[key(rw) %in% key(cf)[idx], , drop = FALSE]
        full$warnings <- rbind(
            full$warnings,
            rw,
            data.frame(
                interval = cf$interval[idx],
                nirs_channels = cf$nirs_channels[idx],
                type = "warning",
                message = cli::ansi_strip(paste(heads, reasons, sep = " - "))
            )
        )
        rownames(full$warnings) <- NULL
    }

    full$coefficients <- add_model_col(cf, cf$model)
    return(full)
}


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

        ## the nls formula is built on the channel names, so a channel
        ## named after a model parameter is aliased by `fit_names()`
        params <- unique(unlist(lapply(per_channel, \(.a) {
            if (is.function(fix_params)) fix_params(.a) else fix_params
        })))
        clash <- intersect(c(nirs_channels, time_channel), c(params, "D"))
        if (verbose && length(clash) > 0L) {
            alias <- paste0(".", clash)
            cli_inform(c(
                "i" = "{cli::qty(clash)}Channel name{?s} {.val {clash}} \\
                collide{?s/} with model parameter{?s}. Renamed: \\
                {.val {alias}}."
            ))
        }
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
#'   *without* `interval`/`nirs_channels`), `model`, `fitted_data`
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
#'   `nirs_channels`, and method parameters), one row per channel, with
#'   attributes `"time_channel"` (the resolved time column name), `"model"`
#'   and `"fitted_data"` (named lists by channel), `"diagnostics"` and
#'   `"channel_args"` (data frames, one row per channel), and `"warnings"`
#'   (data frame of conditions captured during fitting, regardless of
#'   `verbose`; zero rows when none fire).
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

                ## filter valid finite idx before first extreme + end_window;
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
                ## dropping internal-only args, so they fit a flat df row
                arg_row <- lapply(c(.a, extra_args), \(.x) {
                    if (is.null(.x)) {
                        NA
                    } else if (is.list(.x)) {
                        ## deparse() wraps beyond its default width, which would
                        ## expand the single-row data frame
                        paste(deparse(.x), collapse = "")
                    } else if (length(.x) > 1L) {
                        ## collapse vector args (e.g. multiple `fraction`
                        ## values) to fit the single-row data frame
                        paste(.x, collapse = ", ")
                    } else {
                        .x
                    }
                })
                arg_row[c("verbose", "bypass_checks", "interval_name")] <- NULL

                list(
                    coefficients = cbind(
                        data.frame(
                            interval = interval_name,
                            nirs_channels = .nirs
                        ),
                        fit$coefs
                    ),
                    model = fit$model,
                    fitted_data = fit$fitted_data,
                    diagnostics = cbind(
                        data.frame(
                            interval = interval_name,
                            nirs_channels = .nirs
                        ),
                        fit$diag
                    ),
                    channel_args = data.frame(
                        interval = interval_name,
                        nirs_channels = .nirs,
                        arg_row
                    )
                )
            })
            )

            ## interval-level conditions from here on
            .nirs_active <- NA_character_

            ## assemble single attributed df (consumed build_kinetics_results)
            result <- structure(
                do.call(rbind, lapply(fits, `[[`, "coefficients")),
                time_channel = time_channel,
                model = lapply(fits, `[[`, "model"),
                fitted_data = lapply(fits, `[[`, "fitted_data"),
                diagnostics = do.call(rbind, lapply(fits, `[[`, "diagnostics")),
                channel_args = do.call(
                    rbind, lapply(fits, `[[`, "channel_args")
                )
            )

            ## warn when time coefficients are negative
            ## (response before start_time)
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
                .a$fraction, Inf, c(0, 1),
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
#'   (*without* `interval`/`nirs_channels`, which are added upstream),
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
                # fmt: skip
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
#'   to a data frame with the response and time columns named per
#'   [fit_names()] and returning an [nls][stats::nls] model or `NULL`.
#'   `on_error(e)` reports
#'   the condition `e` and returns `NULL`, so it doubles as a [tryCatch()]
#'   error handler.
#' @param time_channel Character; resolved time column name.
#' @param retry Logical; attempt the reduced model when the TD fit fails.
#'   A condition of class `"mnirs_fit_final"` (see [fit_final_error()])
#'   reports a failure the reduced model cannot resolve and skips the retry.
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
    time_channel,
    interval_name,
    env,
    retry = .a$use_TD && !"TD" %in% names(.a$fix)
) {
    attempt <- \(.params, .retry) {
        ## dropping TD narrows the window, so subset per attempt. columns
        ## carry the channel names so the model predicts on them; the full
        ## `params` alias identically across both attempts
        keep <- "TD" %in% .params | t_fit >= 0
        data <- setNames(
            data.frame(x_fit[keep], t_fit[keep]),
            fit_names(.nirs, time_channel, params)
        )
        on_error <- \(e) {
            if (inherits(e, "mnirs_fit_final")) {
                retry <<- FALSE
            }
            warn_fit_failed(
                fn,
                e,
                .nirs,
                interval_name,
                length(.params),
                .retry && retry,
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


#' Fit error the reduced model cannot resolve
#'
#' An error condition for [fit_td_fallback()] `on_error` handlers that
#' describes the data rather than the attempt (e.g. inseparable
#' biexponential phases), so the reduced-model retry is skipped.
#'
#' @param message Character; the condition message.
#'
#' @returns A condition of class `"mnirs_fit_final"` and `"simpleError"`.
#'
#' @keywords internal
fit_final_error <- function(message) {
    return(structure(
        class = c("mnirs_fit_final", "simpleError", "error", "condition"),
        list(message = message, call = NULL)
    ))
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


#' Free parameters of a self-start model call
#'
#' A parameter written as a bare symbol in the model call is free (fitted
#' by [stats::nls()]); one written as a constant or expression is fixed.
#' Used by the model functions to return gradient columns for the free
#' parameters only, in call order, as [stats::nls()] indexes the
#' `"gradient"` attribute by position.
#'
#' @param mCall A matched call to the model function.
#' @param params Character vector of the model parameter names.
#'
#' @returns A character vector; the subset of `params` that are free.
#'
#' @keywords internal
free_params <- function(mCall, params) {
    return(params[vapply(params, \(.p) is.name(mCall[[.p]]), logical(1))])
}


#' Batched 3-parameter least squares over a grid
#'
#' Solves the normal equations of a 3-column linear model at every grid
#' point at once, given the Gram entries `g_ij = <c_i, c_j>` and
#' right-hand sides `b_i = <c_i, x>` as equal-shaped arrays (one element
#' per grid point). Used by the self-start initialisers to profile the
#' non-linear parameters on a grid without a per-point decomposition.
#'
#' @param g11,g12,g13,g22,g23,g33 Gram entries of the three basis columns.
#' @param b1,b2,b3 Inner products of the basis columns with the response.
#' @param xx The response sum of squares `<x, x>`.
#'
#' @returns A list with the coefficient arrays `c1`, `c2`, `c3` and the
#'   residual sum of squares `rss`, which is `Inf` where the system is
#'   singular.
#'
#' @keywords internal
solve_grid3 <- function(g11, g12, g13, g22, g23, g33, b1, b2, b3, xx) {
    ## symmetric 3x3 inverse by adjugate, vectorised over grid points
    a11 <- g22 * g33 - g23^2
    a12 <- g13 * g23 - g12 * g33
    a13 <- g12 * g23 - g13 * g22
    a22 <- g11 * g33 - g13^2
    a23 <- g12 * g13 - g11 * g23
    a33 <- g11 * g22 - g12^2
    det <- g11 * a11 + g12 * a12 + g13 * a13
    c1 <- (a11 * b1 + a12 * b2 + a13 * b3) / det
    c2 <- (a12 * b1 + a22 * b2 + a23 * b3) / det
    c3 <- (a13 * b1 + a23 * b2 + a33 * b3) / det
    ## at the least-squares optimum rss = <x, x> - <b, c>
    rss <- xx - (b1 * c1 + b2 * c2 + b3 * c3)
    rss[!is.finite(rss + c1 + c2 + c3)] <- Inf
    return(list(c1 = c1, c2 = c2, c3 = c3, rss = rss))
}


#' Build a self-start model formula with optional fixed parameters
#'
#' Constructs `x ~ fn(t, ...)` with each free parameter as a bare
#' symbol and each fixed parameter substituted as its constant value.
#'
#' @param fn Symbol; the self-start model function.
#' @param params Character vector of parameter names in `fn` argument
#'   order.
#' @param fix Named list of fixed parameter values.
#' @param x,t Character; the response and time column names (see
#'   [fit_names()]), so the returned model predicts on the original
#'   channel names.
#'
#' @returns A two-sided [formula][stats::formula] on `x` and `t`.
#'
#' @keywords internal
build_ss_formula <- function(fn, params, fix = list(), x, t) {
    args <- lapply(setNames(nm = params), \(.p) fix[[.p]] %||% as.name(.p))
    rhs <- as.call(c(fn, as.name(t), args))
    return(stats::as.formula(call("~", as.name(x), rhs)))
}


#' Alias fit column names that collide with model parameters
#'
#' [stats::nls()] formula symbols must be disjoint: a name cannot be
#' both a data column and a parameter. The fit data frame carries the
#' channel names so the model predicts on them, so a channel named
#' after a model parameter (or `D`, the amplitude used by
#' [enforce_direction()]) is prefixed with `.` in the fit and in the
#' stored model formula. Coefficients and results keep the original
#' names.
#'
#' @param x,t Character; the response and time channel names.
#' @param params Character vector of the model parameter names.
#'
#' @returns A length-2 character vector of the response and time column
#'   names to fit on.
#'
#' @keywords internal
fit_names <- function(x, t, params) {
    nm <- c(x, t)
    clash <- nm %in% c(params, "D")
    nm[clash] <- paste0(".", nm[clash])
    return(nm)
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
#' Direction is the sign of the primary amplitude `D = B - A` (`B1 - A`
#' for the biexponential), where the primary asymptote follows `A` in
#' model parameter order. A fit is satisfied when `D` and every free
#' parameter lie inside the refit box (`D` sign-constrained; `lower`/
#' `upper` for the rest, e.g. a sigmoid slope sign floor) and is
#' returned unchanged. Otherwise the model is refit on `D` via
#' `nls(algorithm = "port")` with `D` sign-bounded and its magnitude
#' floored strictly above zero (sigmoid models divide by `D`), then
#' re-expressed in the original parameterisation from that optimum so
#' the returned model reports consistent coefficient names. A refit
#' that fails, pins a sign-floored coefficient (degenerate flat fit),
#' or loses the requested sign on re-expression warns and returns
#' `NULL`. Parameters in `fix` are held constant: a fixed `A` or `B` is
#' substituted into the amplitude reparameterisation; with both fixed
#' the amplitude sign is predetermined and no refit is possible.
#'
#' @param model A converged [nls][stats::nls] model object.
#' @param coefs Named numeric coefficient vector in model parameter
#'   order with fixed values merged in (see [full_coefs()]).
#' @param fit_data Data frame with the response in the first column and
#'   time in the second; the refit formula is built on those names.
#' @param direction Character; resolved `"positive"` or `"negative"`.
#' @param amp_fn Symbol; model fn taking `t` and the model parameters as
#'   named arguments. A self-start fn returning a `"gradient"` attribute
#'   (free symbols, asymptotes first) makes both refits analytic; a plain
#'   fn falls through to [stats::numericDeriv()].
#' @param fn Character; the self-start fn named in the warning
#'   (*default* `SS<amp_fn>`, or `amp_fn` itself when already prefixed).
#' @param lower,upper Named numeric bounds for free parameters other
#'   than the asymptotes. Sign-floor bounds should be data-scaled
#'   small values (not `.Machine$double.eps`) so pinned-floor
#'   degeneracy is detectable.
#' @param floor_params Character; names of refit coefficients subject
#'   to the pinned-floor degeneracy check. `NULL` (*default*) checks
#'   every finite bound; restrict when other bounds are structural
#'   (e.g. the biexponential time-constant bounds).
#' @param fix Named list of user-fixed parameter values.
#' @param .nirs Character; the channel name.
#' @param interval_name Character; the interval label.
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
    fn = sub("^(SS)?", "SS", as.character(amp_fn)),
    lower = NULL,
    upper = NULL,
    floor_params = NULL,
    fix = list(),
    .nirs,
    interval_name,
    env = rlang::caller_env()
) {
    ## the primary asymptote follows A in model parameter order (B, or
    ## B1 for the biexponential); direction is the sign of D = B - A
    B_name <- names(coefs)[2]
    want <- if (direction == "positive") 1 else -1
    x_sym <- as.name(names(fit_data)[[1L]])
    t_sym <- as.name(names(fit_data)[[2L]])
    x_span <- diff(range(fit_data[[1L]]))
    extra <- coefs[setdiff(names(coefs), c("A", B_name, names(fix)))]

    ## refit box: D sign-constrained with a strictly positive magnitude
    ## floor (sigmoid models divide by D); other params unbounded
    ## unless overridden
    lwr <- setNames(rep(-Inf, 2 + length(extra)), c("A", "D", names(extra)))
    upr <- -lwr
    lwr[["D"]] <- if (want > 0) x_span * 1e-6 else -Inf
    upr[["D"]] <- if (want < 0) -x_span * 1e-6 else Inf
    lwr[names(lower)] <- lower
    upr[names(upper)] <- upper

    ## keep the fit when its coefficients already sit inside the box
    cf0 <- c(D = coefs[[B_name]] - coefs[["A"]], extra)
    if (all(cf0 >= lwr[names(cf0)] & cf0 <= upr[names(cf0)])) {
        return(list(model = model, coefs = coefs))
    }

    fail <- \() {
        cli_warn(c(
            "x" = "{.fn {fn}} fit for {.field {(.nirs)}} in \\
            {.field {interval_name}} could not satisfy \\
            {.code direction = {.val {direction}}}.",
            "i" = "Returning {.val {NA}} coefficients."
        ), call = warn_call(env))
        return(NULL)
    }

    ## both asymptotes fixed: amplitude sign is predetermined; no refit
    A_fix <- fix[["A"]]
    B_fix <- fix[[B_name]]
    if (!is.null(A_fix) && !is.null(B_fix)) {
        return(fail())
    }

    ## bounded refit on `start`, boxed by name. the start may be far
    ## from its constrained optimum, so it gets the same iteration
    ## budget as the primary fit; a finite-coefficient PORT stall
    ## (non-smooth onset kink) is accepted as in the primary fit
    port_fit <- \(rhs, start) {
        l <- lwr[names(start)]
        u <- upr[names(start)]
        l[is.na(l)] <- -Inf
        u[is.na(u)] <- Inf
        model <- tryCatch(
            suppressWarnings(nls(
                stats::as.formula(call("~", x_sym, rhs)),
                fit_data,
                start = start,
                lower = l,
                upper = u,
                algorithm = "port",
                control = stats::nls.control(maxiter = 500L, warnOnly = TRUE)
            )),
            error = \(e) NULL
        )
        return(accept_port_fit(model, \(e) NULL))
    }

    ## model rhs on the asymptote expressions; other params ride as
    ## named args, fixed as constants and free as symbols
    rhs <- \(A_expr, B_expr) {
        as.call(c(
            amp_fn,
            t_sym,
            setNames(list(A_expr, B_expr), c("A", B_name)),
            fix[setdiff(names(fix), c("A", B_name))],
            lapply(setNames(nm = names(extra)), as.name)
        ))
    }

    ## seeds: D in the requested direction, confined to the observed
    ## scale so a stranded asymptote cannot seed a runaway; params
    ## outside the box (e.g. an inverted slope) are mirrored into it
    D0 <- want * min(max(abs(cf0[["D"]]), x_span * 0.1), x_span)
    ex_l <- lwr[names(extra)]
    ex_u <- upr[names(extra)]
    bad <- extra < ex_l | extra > ex_u
    extra[bad] <- pmin(pmax(-extra[bad], ex_l[bad]), ex_u[bad])

    ## refit on amplitude D, substituting a fixed asymptote: A free
    ## `B = A + D`, A fixed `B = A_fix + D`, B fixed `A = B_fix - D`.
    ## the model gradient over its free symbols (asymptotes first) maps
    ## onto (A, D) by the chain rule; nls resolves `amp_D` from the
    ## formula environment, this frame
    A_free <- is.null(A_fix) && is.null(B_fix)
    amp <- eval(amp_fn)
    amp_D <- function(t, D, A = NULL, ...) {
        A_val <- if (is.null(B_fix)) A_fix %||% A else B_fix - D
        B_val <- B_fix %||% (A_val + D)
        val <- amp(t, A_val, B_val, ...)
        g <- attr(val, "gradient")
        if (is.null(g)) {
            return(val)
        }
        attr(val, "gradient") <- cbind(
            if (A_free) cbind(A = g[, 1L] + g[, 2L]),
            D = if (is.null(B_fix)) g[, 2L] else -g[, 1L],
            g[, -(1:2), drop = FALSE]
        )
        return(val)
    }
    refit <- port_fit(
        as.call(c(
            quote(amp_D),
            t_sym,
            D = quote(D),
            if (A_free) list(A = quote(A)),
            fix[setdiff(names(fix), c("A", B_name))],
            lapply(setNames(nm = names(extra)), as.name)
        )),
        c(if (A_free) c(A = coefs[["A"]]), D = D0, extra)
    )
    if (is.null(refit)) {
        return(fail())
    }

    ## a coefficient pinned at a sign floor (D, slope, tau) marks a
    ## degenerate flat fit: no genuine response in the requested
    ## direction. `floor_params` excludes structural bounds
    cf <- coef(refit)
    bnd <- c(lwr, upr)
    bnd <- bnd[is.finite(bnd)]
    chk <- intersect(floor_params %||% names(cf), names(bnd))
    if (any(abs(cf[chk]) <= 2 * abs(bnd[chk]))) {
        return(fail())
    }

    ## re-express in the original parameterisation from the interior
    ## optimum, so the returned model reports (A, B, ...) rather than D
    A_val <- if (is.null(B_fix)) A_fix %||% cf[["A"]] else B_fix - cf[["D"]]
    refit2 <- port_fit(
        rhs(A_fix %||% quote(A), B_fix %||% as.name(B_name)),
        c(
            if (is.null(A_fix)) c(A = A_val),
            if (is.null(B_fix)) setNames(A_val + cf[["D"]], B_name),
            cf[names(extra)]
        )
    )
    if (is.null(refit2)) {
        return(fail())
    }
    out <- full_coefs(refit2, names(coefs), fix)
    if (want * (out[[B_name]] - out[["A"]]) <= 0) {
        return(fail())
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
#'   `adj_r2`, `rmse`, `cv_rmse`, `snr`, `aic`, `aicc`, and `bic`.
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
        cv_rmse = NA_real_,
        snr = NA_real_,
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
        cv_rmse = cv_rmse,
        snr = snr,
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
