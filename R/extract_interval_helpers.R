#' Specify interval boundaries by time, label, lap, or sample
#'
#' Helper functions to define interval start or end boundaries for
#' [extract_intervals()].
#'
#' @param ... Specify start or end boundaries.
#'   \describe{
#'     \item{`by_time(...)`}{Numeric time values in units of `time_channel`.}
#'     \item{`by_label(...)`}{Character strings to match in `event_channel`.
#'     Matched as regular expressions by default; see `ignore_case` and
#'     `fixed`. All matching occurrences are returned.}
#'     \item{`by_lap(...)`}{Integer lap numbers to match in `event_channel`.
#'     For `start`, resolves to the first sample of each lap. For `end`,
#'     resolves to the last sample.}
#'     \item{`by_sample(...)`}{Integer sample indices (row numbers).}
#'   }
#' @param ignore_case For `by_label()`. If `TRUE`, match case-insensitive
#'   labels. Default `FALSE`.
#' @param fixed For `by_label()`. If `TRUE`, treat labels as fixed strings
#'   rather than regular expressions. Useful when labels contain regex
#'   metacharacters (`.`, `*`, `(`, etc.). Default `FALSE`.
#'
#' @details
#' These helpers can be used explicitly for arguments `start`/`end`, or raw
#' values can be passed directly:
#'   - Numeric -> [by_time()]
#'   - Character -> [by_label()],
#'   - Explicit integer (e.g. `2L`) -> [by_lap()].
#'   - Use [by_sample()] explicitly for sample indices.
#'
#' Multiple specification types can be combined for a single boundary with
#' `list()` (e.g. `list(by_time(30), by_label("go"))`). Resolved
#' boundary times are concatenated in the order supplied. Combined
#' specifications must use the `by_` helpers directly: raw values are
#' ignored with a warning.
#'
#' @returns An object of class `"mnirs_interval"` for use with the `start`
#'   and `end` arguments of [extract_intervals()].
#'
#' @examples
#' ## read example data
#' data <- read_mnirs(
#'     example_mnirs("train.red"),
#'     nirs_channels = c(
#'         smo2_left = "SmO2 unfiltered",
#'         smo2_right = "SmO2 unfiltered"
#'     ),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     event_channel = c(lap = "Lap/Event"),
#'     zero_time = TRUE,
#'     verbose = FALSE
#' )
#'
#' ## start and end by time
#' extract_intervals(data, start = by_time(66), end = by_time(357))
#'
#' ## start by lap
#' extract_intervals(data, start = by_lap(2, 4), span = 0)
#'
#' ## combine multiple specification types
#' extract_intervals(
#'     data, 
#'     start = list(by_lap(2), by_time(400)), 
#'     end = by_sample(1500)
#' )
#'
#' ## simulate event_channel with character label match 
#' data$event <- NA_character_
#' data$event[c(1000, 1001)] <- c("start", "lap.1")
#' data <- create_mnirs_data(data, event_channel = "event")
#'
#' ## case-insensitive label match
#' extract_intervals(data, start = by_label("START", ignore_case = TRUE))
#'
#' ## literal-string label match (regex metacharacters treated as text)
#' extract_intervals(data, start = by_label("lap.1", fixed = TRUE))
#'
#' @export
by_time <- function(...) {
    by_time <- c(...)
    validate_numeric(by_time)
    structure(
        list(type = "time", by_time = by_time),
        class = "mnirs_interval"
    )
}


#' @rdname by_time
#' @export
by_label <- function(..., ignore_case = FALSE, fixed = FALSE) {
    by_label <- c(...)
    if (!is.character(by_label) || length(by_label) == 0L) {
        cli_abort(c(
            "x" = "{.fn by_label} must be a valid {.cls character} vector."
        ))
    }
    ## collect any invalid logical flags, report together
    is_bad_flag <- function(x) {
        !is.logical(x) || length(x) != 1L || is.na(x)
    }
    bad <- c(
        if (is_bad_flag(ignore_case)) "ignore_case",
        if (is_bad_flag(fixed)) "fixed"
    )
    if (length(bad) > 0L) {
        cli_abort(c(
            "x" = "{.fn by_label} {cli::qty(bad)}arg{?s} {.arg {bad}} \\
            must be {.val {TRUE}} or {.val {FALSE}}."
        ))
    }

    structure(
        list(
            type = "label",
            by_label = by_label,
            ignore_case = ignore_case,
            fixed = fixed
        ),
        class = "mnirs_interval"
    )
}


#' @rdname by_time
#' @export
by_lap <- function(...) {
    by_lap <- c(...)
    validate_numeric(by_lap, range = c(1, Inf), integer = TRUE)
    structure(
        list(type = "lap", by_lap = as.integer(by_lap)),
        class = "mnirs_interval"
    )
}


#' @rdname by_time
#' @export
by_sample <- function(...) {
    by_sample <- c(...)
    validate_numeric(by_sample, range = c(1, Inf), integer = TRUE)
    structure(
        list(type = "sample", by_sample = as.integer(by_sample)),
        class = "mnirs_interval"
    )
}


#' coerce raw values to mnirs_interval objects
#' @param x A raw value or mnirs_interval object.
#' @param arg Name of the argument for error messages.
#' @inheritParams validate_mnirs
#' @keywords internal
as_mnirs_interval <- function(x, arg = "start", env = rlang::caller_env()) {
    if (is.null(x) || inherits(x, c("mnirs_interval", "mnirs_interval_list"))) {
        return(x)
    }
    ## bare lists combine multiple by_* specs
    if (is.list(x)) {
        ## c() flattens spec components into a plain named list
        if ("type" %in% names(x)) {
            cli_abort(c(
                "x" = "{.arg {arg}}: specifications combined with {.fn c}.",
                "i" = "Combine multiple {.fn by_*} specifications with \\
                {.fn list}."
            ), call = env)
        }
        ## splice nested containers flat, drop NULL elements
        specs <- unlist(lapply(x, \(.x) {
                if (inherits(.x, "mnirs_interval_list")) {
                    return(unclass(.x))
                }
                list(.x)
            }), recursive = FALSE)
        specs <- specs[lengths(specs) > 0L]
        ## combined specs must use by_* directly: warn & ignore raw values
        valid <- vapply(specs, inherits, logical(1), "mnirs_interval")
        if (any(!valid)) {
            cli_warn(c(
                "!" = "{.arg {arg}}: ignoring {sum(!valid)} element{?s}.",
                "i" = "Multiple {.arg start} and/or {.arg end} values must \\
                be specified with {.fn by_*}."
            ), call = warn_call(env))
        }
        specs <- specs[valid]
        if (length(specs) == 0L) {
            return(NULL)
        }
        if (length(specs) == 1L) {
            return(specs[[1L]])
        }
        return(structure(specs, class = "mnirs_interval_list"))
    }
    ## integer before numeric — integers are also numeric in R
    if (is.integer(x)) {
        return(by_lap(x))
    }
    if (is.numeric(x)) {
        return(by_time(x))
    }
    if (is.character(x)) {
        return(by_label(x))
    }
    cli_abort(c(
        "x" = "{.arg {arg}} must be {.cls numeric}, {.cls integer}, \\
        {.cls character}, or a {.fn by_time}, {.fn by_sample}, \\
        {.fn by_label}, {.fn by_lap} specification."
    ), call = env)
}


#' Validate per-group channel selections for ensemble-averaging
#'
#' Normalises `group_channels` to a list of channel-name vectors for
#' recycling across interval groups. `NULL` selects all `nirs_channels` for
#' every group. Unlike [validate_group_channels()], channels may repeat
#' across list items: each item is one group's channel selection.
#'
#' @inheritParams validate_mnirs
#' @keywords internal
validate_interval_channels <- function(
    group_channels,
    nirs_channels,
    data = NULL,
    env = rlang::caller_env()
) {
    ## parse tidy eval input
    if (rlang::is_quosure(group_channels)) {
        group_channels <- parse_channel_name(
            group_channels,
            data,
            rlang::quo_get_env(group_channels)
        )
    }
    if (is.null(group_channels)) {
        return(list(nirs_channels))
    }
    if (!is.list(group_channels)) {
        group_channels <- list(group_channels)
    }

    ## group members must be known channels
    members <- unlist(group_channels, use.names = FALSE)
    unknown <- setdiff(members, nirs_channels)
    if (!is.character(members) || length(unknown) > 0L) {
        cli_abort(c(
            "x" = "{.arg group_channels}: channel{?s} {.field {unknown}} \\
            not recognised.",
            "i" = "Grouped channel names must match {.arg nirs_channels} \\
            exactly."
        ), call = env)
    }
    return(group_channels)
}


#' recycle a single-element span to c(before, after)
#' positive -> c(0, x), negative -> c(x, 0)
#' @inheritParams validate_mnirs
#' @keywords internal
recycle_span <- function(span, env = rlang::caller_env()) {
    if (is.numeric(span) && length(span) == 2L) {
        return(span)
    }
    if (!is.numeric(span) || length(span) != 1L) {
        cli_abort(c(
            "x" = "{.arg span} must be a one- or two-element \\
            {.cls numeric} vector."
        ), call = env)
    }
    ## positive shifts end, negative shifts start
    return(if (span >= 0) c(0, span) else c(span, 0))
}


#' resolve a single mnirs_interval object to time values
#' @inheritParams validate_mnirs
#' @keywords internal
find_interval_time <- function(
    interval,
    t_vec,
    event_vec = NULL,
    position = c("first", "last"),
    env = rlang::caller_env()
) {
    ## multi-spec container: resolve each spec, concatenate in supplied order
    if (inherits(interval, "mnirs_interval_list")) {
        return(unlist(lapply(interval, \(.x) {
            find_interval_time(.x, t_vec, event_vec, position, env)
        }), use.names = FALSE))
    }
    switch(
        interval$type,
        time = interval$by_time,
        sample = t_vec[interval$by_sample],
        label = {
            ## OR per-label matches; grepl(fixed = TRUE) can't take a|b
            matches <- which(Reduce(`|`, lapply(interval$by_label, \(.p) {
                grepl(.p,
                    event_vec,
                    ignore.case = interval$ignore_case,
                    fixed = interval$fixed
                )
            })))
            if (length(matches) == 0L) {
                cli_abort(c(
                    "x" = "No events detected matching \\
                    {.val {interval$by_label}}.",
                    "i" = "Check {.arg event_channel} contents; see \\
                    {.arg ignore_case} and {.arg fixed} in {.fn by_label}."
                ), call = env)
            }
            t_vec[matches]
        },
        lap = {
            ## convert once outside loop
            event_int <- as.integer(event_vec)
            vapply(interval$by_lap, \(lap_val) {
                matches <- which(event_int == lap_val)
                if (length(matches) == 0L) {
                    cli_abort(c(
                        "x" = "No samples found for lap {.val {lap_val}}.",
                        "i" = "Check that {.arg event_channel} contains \\
                        lap numbers."
                    ), call = env)
                }
                idx <- if (position == "first") 1L else length(matches)
                t_vec[matches[idx]]
            }, numeric(1))
        }
    )
}


#' resolve start/end into time value vectors (no span applied)
#' @inheritParams validate_mnirs
#' @keywords internal
resolve_interval <- function(
    start,
    end,
    t_vec,
    event_vec = NULL,
    env = rlang::caller_env()
) {
    has_start <- !is.null(start)
    has_end <- !is.null(end)
    interval <- start %||% end

    ## single-boundary non-lap: reference point for span
    if (xor(has_start, has_end)) {
        start_time <- find_interval_time(
            interval,
            t_vec,
            event_vec,
            position = if (has_start) "first" else "last",
            env = env
        )
        return(list(
            start_time = start_time,
            end_time = NULL,
            has_start = has_start,
            has_end = has_end
        ))
    }

    ## both boundaries specified
    start_time <- find_interval_time(start, t_vec, event_vec, "first", env)
    end_time <- find_interval_time(end, t_vec, event_vec, "last", env)

    ## warn and truncate if lengths differ
    n_start <- length(start_time)
    n_end <- length(end_time)
    if (n_start != n_end) {
        n_intervals <- min(n_start, n_end)
        cli_warn(c(
            "!" = "Unequal lengths for {.arg start} ({.val {n_start}}) \\
            and {.arg end} ({.val {n_end}}).",
            "i" = "Returning {.val {n_intervals}} paired interval{?s}."
        ), call = warn_call(env))
        start_time <- start_time[seq_len(n_intervals)]
        end_time <- end_time[seq_len(n_intervals)]
    }

    return(list(
        start_time = start_time,
        end_time = end_time,
        has_start = TRUE,
        has_end = TRUE
    ))
}


#' Recycle parameter list to target length
#' @inheritParams validate_mnirs
#' @keywords internal
recycle_to_length <- function(
    param,
    n,
    name = c("event", "group"),
    verbose = TRUE,
    env = rlang::caller_env(),
    arg = "values"
) {
    n_param <- length(param)

    if (n_param == n) {
        return(param)
    }

    if (n_param > n) {
        if (verbose) {
            cli_inform(c(
                "!" = "{.arg {arg}} exceeds the number of \\
                {name}s by {.val {n_param - n}}.",
                "i" = "Extra values are ignored."
            ), call = env)
        }
        return(param[seq_len(n)])
    }

    ## n_param < n:  recycle last element forward
    if (verbose && n_param > 1L) {
        cli_inform(c(
            "i" = "{.arg {arg}} recycled to meet \\
            {.val {n - n_param}} unspecified {name}{qty(n - n_param)}{?s}."
        ), call = env)
    }
    return(param[c(seq_len(n_param), rep(n_param, n - n_param))])
}


## reject malformed group indices before they reach vector subscripting
validate_interval_groups <- function(
    group_intervals,
    n_intervals,
    env = rlang::caller_env()
) {
    if (length(group_intervals) == 0L) {
        cli_abort(c(
            "x" = "{.arg group_intervals} must contain at least one group."
        ), call = env)
    }

    spec <- group_intervals[[1L]][1L]
    if (
        length(group_intervals) == 1L &&
            isTRUE(spec %in% c("distinct", "ensemble"))
    ) {
        return(invisible())
    }

    group_labels <- names(group_intervals) %||%
        character(length(group_intervals))
    unnamed <- is.na(group_labels) | !nzchar(group_labels)
    group_labels[unnamed] <- paste("position", which(unnamed))

    empty <- lengths(group_intervals) == 0L
    if (any(empty)) {
        cli_abort(c(
            "x" = "{.arg group_intervals}: empty group{?s} \\
            {.field {group_labels[empty]}}.",
            "i" = "Each group must contain at least one interval index."
        ), call = env)
    }

    ## is_integerish(finite = TRUE) rejects non-numeric, NA, Inf, & fractions
    valid <- vapply(
        group_intervals,
        \(.g) rlang::is_integerish(.g, finite = TRUE) && is.null(dim(.g)),
        logical(1)
    )
    if (!all(valid)) {
        cli_abort(c(
            "x" = "{.arg group_intervals}: group{?s} \\
            {.field {group_labels[!valid]}} contains missing integer indices."
        ), call = env)
    }

    oob <- lapply(group_intervals, \(.g) .g[.g < 1L | .g > n_intervals])
    out_of_range <- lengths(oob) > 0L
    if (any(out_of_range)) {
        invalid <- unique(unlist(oob, use.names = FALSE))
        cli_abort(c(
            "x" = "{.arg group_intervals} contains out-of-range interval \\
            indices in group{qty(sum(out_of_range))}{?s} \\
            {.field {group_labels[out_of_range]}}: {.field {invalid}}.",
            "i" = "Indices must be between {.val 1} and \\
            {.val {n_intervals}}."
        ), call = env)
    }

    return(invisible())
}


#' Recycle parameter to match number of events
#'
#' Recycle an argument vector to a list or repeat the last list item to match
#' the number of events.
#'
#' @keywords internal
recycle_param <- function(
    param,
    n_events,
    group_intervals,
    verbose = TRUE,
    env = rlang::caller_env(),
    arg = "values"
) {
    ## flatten nested lists to single-depth list
    param <- if (is.list(param)) {
        lapply(param, \(.x) if (is.list(.x)) unlist(.x) else .x)
    } else {
        list(param)
    }

    ## standard recycling to events for "distinct" or "ensemble"
    if (!is.numeric(group_intervals[[1L]])) {
        return(recycle_to_length(param, n_events, "event", verbose, env, arg))
    }

    ## custom grouping: recycle per group, then map to event order;
    ## ungrouped events take the last group's value
    n_groups <- length(group_intervals)
    param <- recycle_to_length(param, n_groups, "group", verbose, env, arg)
    ids <- unlist(group_intervals, use.names = FALSE)
    keep <- ids <= n_events
    event_to_group <- rep(n_groups, n_events)
    event_to_group[ids[keep]] <- rep(
        seq_len(n_groups),
        lengths(group_intervals)
    )[keep]
    return(param[event_to_group])
}


#' apply span to resolved times and build interval_spec data frame
#' @inheritParams validate_mnirs
#' @keywords internal
apply_span <- function(
    interval_list,
    t_vec,
    span,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## extract span values (recycled by recycle_param)
    span_before <- vapply(span, `[`, numeric(1), 1L)
    span_after <- vapply(span, `[`, numeric(1), 2L)
    event_times <- interval_list$start_time

    if (interval_list$has_start && interval_list$has_end) {
        ## span[1] shifts starts, span[2] shifts ends
        start_times <- event_times + span_before
        end_times <- interval_list$end_time + span_after
        ## original boundary times for metadata
        interval_times <- Map(c, event_times, interval_list$end_time)
    } else {
        ## single-boundary: both span values on reference
        start_times <- event_times + span_before
        end_times <- event_times + span_after
        interval_times <- as.list(event_times)
    }

    ## data time bounds
    t_min <- t_vec[1L]
    t_max <- t_vec[length(t_vec)]

    ## check for entirely out of bounds intervals
    entirely_oob <- end_times < t_min | start_times > t_max
    if (any(entirely_oob)) {
        oob_ids <- which(entirely_oob)
        n_oob <- qty(length(oob_ids))
        cli_abort(c(
            "x" = "{n_oob} Interval{?s} {.field {oob_ids}} {n_oob} \\
            {?is/are} entirely outside data bounds.",
            "i" = "Intervals must be specified within existing data bounds."
        ), call = env)
    }

    ## clamp partial out of bounds to data range
    partial_oob <- start_times < t_min | end_times > t_max
    if (any(partial_oob)) {
        start_times <- pmax(t_min, start_times)
        end_times <- pmin(t_max, end_times)

        if (verbose) {
            oob_ids <- which(partial_oob)
            n_oob <- qty(length(oob_ids))
            cli_warn(c(
                "!" = "{n_oob} Interval{?s} {.field {oob_ids}} {n_oob} \\
                {?is/are} partially outside data bounds.",
                "i" = "Returning available data only."
            ), call = warn_call(env))
        }
    }

    result <- data.frame(
        span_before = span_before,
        span_after = span_after,
        start_times = start_times,
        end_times = end_times,
        stringsAsFactors = FALSE
    )
    result$interval_times <- interval_times
    return(result)
}


#' Extract interval data by time range
#' @keywords internal
extract_df_list <- function(
    data,
    t_vec,
    interval_spec,
    group_channels ## list of per-interval channel vectors
) {
    n_vec <- seq_len(nrow(interval_spec))
    df_list <- lapply(n_vec, \(.i) {
        time_range <- t_vec >= interval_spec$start_times[.i] &
            t_vec <= interval_spec$end_times[.i]
        interval_data <- data[time_range, , drop = FALSE]

        ## return interval_data with metadata
        create_mnirs_data(
            interval_data,
            nirs_channels = group_channels[[.i]], ## overwrite for interval data
            interval_times = interval_spec$interval_times[[.i]],
            interval_span = c(
                interval_spec$span_before[.i],
                interval_spec$span_after[.i]
            )
        )
    })

    names(df_list) <- sprintf("interval_%d", n_vec)
    return(df_list)
}


#' Recalculate time_channel values with zero offset at event time (t0)
#' @keywords internal
zero_offset_data <- function(data, time_channel, t0) {
    ## zero time channel to start_time `t0`
    data[[time_channel]] <- data[[time_channel]] - t0
    return(data)
}


#' Ensemble average multiple intervals
#'
#' `group_channels` is a character vector applied to every interval, or a
#' list of per-interval channel vectors; channels excluded from an interval
#' do not contribute to that channel's ensemble-mean.
#'
#' @inheritParams validate_mnirs
#' @keywords internal
ensemble_intervals <- function(
    df_list,
    group_channels,
    metadata,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    time_channel <- metadata$time_channel
    sample_rate <- metadata$sample_rate
    if (!is.list(group_channels)) {
        group_channels <- rep(list(group_channels), length(df_list))
    }
    nirs_channels <- unique(unlist(group_channels, use.names = FALSE))

    ## extract data & metadata from df_list
    interval_data <- Map(\(.df, .channels) {
        interval_times <- attr(.df, "interval_times")
        ## zero-offset to start time (first element of interval_times);
        ## ensemble-averaging time values makes no sense
        t0 <- interval_times[[1L]]
        ## NA-out channels excluded from this interval so binned means
        ## only average member channels
        excluded <- setdiff(nirs_channels, .channels)
        if (length(excluded) > 0L) {
            .df[excluded] <- NA_real_
        }
        list(
            data = zero_offset_data(.df, attr(.df, "time_channel"), t0),
            interval_times = interval_times - t0, ## zero_time recalc
            interval_span = attr(.df, "interval_span")
        )
    }, df_list, group_channels)

    ## stack interval data frames
    df_long <- do.call(rbind, lapply(interval_data, `[[`, "data"))
    ## resample times to nearest estimated sample rate for binned ensembling
    time_resampled <- round(df_long[[time_channel]] * sample_rate) / sample_rate
    unique_times <- sort(unique(time_resampled))
    bins <- match(time_resampled, unique_times)

    ## warn if any time samples have only one value, implying irregular
    ## samples and may result in alternating samples instead of ensemble-means
    if (verbose && any(tabulate(bins) < 2L)) {
        cli_warn(c(
            "!" = "Duplicate or irregular {.arg time_channel} samples \\
            detected after ensemble-averaging.",
            "i" = "Re-sample with {.fn mnirs::resample_mnirs} before \\
            ensemble-averaging.",
            "i" = "Check your resulting data for inconsistent results."
        ), call = warn_call(env))
    }

    ## channel-wise binned means as one matrix operation;
    ## sums / non-NA counts matches colMeans(na.rm = TRUE) per bin;
    ## rows with missing time values are excluded from binning
    valid <- which(!is.na(bins))
    channel_matrix <- as.matrix(df_long[valid, nirs_channels, drop = FALSE])
    channel_sums <- rowsum(channel_matrix, bins[valid], na.rm = TRUE)
    channel_counts <- rowsum(1L - is.na(channel_matrix), bins[valid])
    result <- data.frame(
        setNames(list(unique_times), time_channel),
        as.data.frame(channel_sums / channel_counts)
    )

    ## return with metadata
    return(create_mnirs_data(
        result,
        nirs_device = metadata$nirs_device,
        nirs_channels = nirs_channels,
        time_channel = time_channel,
        event_channel = metadata$event_channel,
        sample_rate = sample_rate,
        start_timestamp = metadata$start_timestamp,
        interval_times = lapply(interval_data, `[[`, "interval_times"),
        interval_span = lapply(interval_data, `[[`, "interval_span")
    ))
}

#' Zero-offset time values and add metadata
#' @keywords internal
preserve_metadata <- function(data, metadata, zero_time = FALSE) {
    if (zero_time) {
        ## zero-offset a single interval's time and interval_times attr
        interval_times <- attr(data, "interval_times")
        t0 <- interval_times[[1L]]
        data <- zero_offset_data(data, metadata$time_channel, t0)
        attr(data, "interval_times") <- interval_times - t0
    }

    return(create_mnirs_data(
        data,
        nirs_device = metadata$nirs_device,
        nirs_channels = unique(attr(data, "nirs_channels")),
        time_channel = metadata$time_channel,
        event_channel = metadata$event_channel,
        sample_rate = metadata$sample_rate,
        start_timestamp = metadata$start_timestamp,
        interval_times = attr(data, "interval_times"),
        interval_span = attr(data, "interval_span")
    ))
}


#' Normalise custom interval grouping to a complete named list
#'
#' Adds intervals missing from a custom `group_intervals` list as
#' single-interval groups, warns on duplicates, and names groups by
#' user-supplied names with `interval_<ids>` fallback.
#'
#' @inheritParams validate_mnirs
#' @keywords internal
normalise_interval_groups <- function(
    group_intervals,
    n_intervals,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    grouped_ids <- unlist(group_intervals, use.names = FALSE)

    ## add ungrouped ids as individual groups, sort by position
    ungrouped_ids <- setdiff(seq_len(n_intervals), grouped_ids)
    if (length(ungrouped_ids) > 0L) {
        group_intervals <- c(group_intervals, as.list(ungrouped_ids))
        group_intervals <- group_intervals[
            order(vapply(group_intervals, median, numeric(1), na.rm = TRUE))
        ]
        if (verbose) {
            cli_inform(c(
                "!" = "Intervals not specified by {.arg group_intervals}.",
                "i" = "Ungrouped intervals included as discrete."
            ), call = env)
        }
    }

    ## warn for duplicated intervals across groups
    dupes <- grouped_ids[duplicated(grouped_ids)]
    if (verbose && length(dupes) > 0L) {
        cli_warn(c(
            "!" = "Duplicates detected of {qty(length(dupes))} \\
            interval{?s} {.field {dupes}}.",
            "i" = "Re-specify {.arg group_intervals} to remove duplicates."
        ), call = warn_call(env))
    }

    ## user-supplied group names win; fall back to member indices
    id_names <- vapply(group_intervals, \(.g) {
        paste0("interval_", paste(.g, collapse = "_"))
    }, character(1))
    names <- names(group_intervals) %||% character(length(group_intervals))
    return(setNames(group_intervals, ifelse(nzchar(names), names, id_names)))
}


#' Apply grouping to intervals
#' @inheritParams validate_mnirs
#' @keywords internal
apply_interval_groups <- function(
    df_list,
    group_channels,
    metadata,
    group_intervals,
    zero_time = FALSE,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    n_intervals <- length(df_list)
    spec <- group_intervals[[1L]][1L]

    ## normalise grouping spec to a named list of interval-index vectors
    groups <- if (n_intervals == 1L || identical(spec, "distinct")) {
        setNames(as.list(seq_len(n_intervals)), names(df_list))
    } else if (identical(spec, "ensemble")) {
        list(ensemble = seq_len(n_intervals))
    } else {
        normalise_interval_groups(group_intervals, n_intervals, verbose, env)
    }

    ## single-interval groups pass through with metadata; multi-interval
    ## groups are ensemble-averaged across their member intervals
    return(lapply(groups, \(.g) {
        if (length(.g) == 1L) {
            return(preserve_metadata(df_list[[.g]], metadata, zero_time))
        }
        ensemble_intervals(
            df_list[.g], group_channels[.g], metadata, verbose, env
        )
    }))
}
