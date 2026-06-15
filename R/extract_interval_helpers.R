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
#' ## introduce event_channel with "start" string
#' data$event <- NA_character_
#' data$event[1000] <- "start"
#' data <- create_mnirs_data(data, event_channel = "event")
#'
#' ## start by label, end by time
#' extract_intervals(data, start = by_label("start"), end = by_time(1500))
#'
#' ## case-insensitive label match
#' extract_intervals(data, start = by_label("START", ignore_case = TRUE))
#'
#' ## literal-string label match (regex metacharacters treated as text)
#' data$event[1000] <- "lap.1"
#' data <- create_mnirs_data(data, event_channel = "event")
#' extract_intervals(data, start = by_label("lap.1", fixed = TRUE))
#'
#' ## multiple intervals by sample index
#' extract_intervals(data, start = by_sample(1000, 1500))
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
#' @keywords internal
as_mnirs_interval <- function(x, arg = "start") {
    if (is.null(x) || inherits(x, "mnirs_interval")) {
        return(x)
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
    ))
}


#' recycle a single-element span to c(before, after)
#' positive -> c(0, x), negative -> c(x, 0)
#' @keywords internal
recycle_span <- function(span) {
    if (is.numeric(span) && length(span) == 2L) {
        return(span)
    }
    if (!is.numeric(span) || length(span) != 1L) {
        cli_abort(c(
            "x" = "{.arg span} must be a one- or two-element \\
            {.cls numeric} vector."
        ))
    }
    ## positive shifts end, negative shifts start
    return(if (span >= 0) c(0, span) else c(span, 0))
}


#' resolve a single mnirs_interval object to time values
#' @keywords internal
find_interval_time <- function(
    interval,
    time_vec,
    event_vec = NULL,
    position = c("first", "last")
) {
    switch(
        interval$type,
        time = interval$by_time,
        sample = time_vec[interval$by_sample],
        label = {
            ## OR per-label matches; grepl(fixed = TRUE) can't take a|b
            matches <- which(Reduce(`|`, lapply(interval$by_label, \(.p) {
                grepl(
                    .p, event_vec, ignore.case = interval$ignore_case, fixed = interval$fixed
                )
            })))
            if (length(matches) == 0L) {
                cli_abort(c(
                    "x" = "No events detected matching \\
                    {.val {interval$by_label}}.",
                    "i" = "Check {.arg event_channel} contents; see \\
                    {.arg ignore_case} and {.arg fixed} in {.fn by_label}."
                ))
            }
            time_vec[matches]
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
                    ))
                }
                idx <- if (position == "first") 1L else length(matches)
                time_vec[matches[idx]]
            }, numeric(1))
        }
    )
}


#' resolve start/end into time value vectors (no span applied)
#' @keywords internal
resolve_interval <- function(
    start,
    end,
    time_vec,
    event_vec = NULL
) {
    has_start <- !is.null(start)
    has_end <- !is.null(end)
    interval <- start %||% end

    ## single-boundary non-lap: reference point for span
    if (xor(has_start, has_end)) {
        start_time <- find_interval_time(
            interval,
            time_vec,
            event_vec,
            position = if (has_start) "first" else "last"
        )
        return(list(
            start_time = start_time,
            end_time = NULL,
            has_start = has_start,
            has_end = has_end
        ))
    }

    ## both boundaries specified
    start_time <- find_interval_time(start, time_vec, event_vec, "first")
    end_time <- find_interval_time(end, time_vec, event_vec, "last")

    ## warn and truncate if lengths differ
    n_start <- length(start_time)
    n_end <- length(end_time)
    if (n_start != n_end) {
        n_intervals <- min(n_start, n_end)
        cli_warn(c(
            "!" = "Unequal lengths for {.arg start} ({col_blue(n_start)}) \\
            and {.arg end} ({col_blue(n_end)}).",
            "i" = "Returning {col_blue(n_intervals)} paired interval{?s}."
        ))
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
#' @keywords internal
recycle_to_length <- function(
    param,
    n,
    name = c("event", "group"),
    verbose = TRUE
) {
    n_param <- length(param)

    if (n_param == n) {
        return(param)
    }

    if (n_param > n) {
        if (verbose) {
            cli_inform(c(
                "!" = "{.arg {substitute(param)}} exceeds the number of \\
                {name}s by {.val {n_param - n}}.",
                "i" = "Extra values are ignored."
            ))
        }
        return(param[seq_len(n)])
    }

    ## n_param < n:  recycle last element forward
    if (verbose && n_param > 1L) {
        cli_inform(c(
            "i" = "{.arg {substitute(param)}} recycled to meet \\
            {.val {n - n_param}} unspecified {name}{qty(n - n_param)}{?s}."
        ))
    }
    return(param[c(seq_len(n_param), rep(n_param, n - n_param))])
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
    verbose = TRUE
) {
    ## flatten nested lists to single-depth list
    param <- if (is.list(param)) {
        lapply(param, \(.x) if (is.list(.x)) unlist(.x) else .x)
    } else {
        list(param)
    }

    ## custom grouping: recycle per group, then map to event order
    if (is.numeric(group_intervals[[1L]])) {
        n_groups <- length(group_intervals)
        groups_unlisted <- unlist(group_intervals)

        ## recycle param to number of groups
        param <- recycle_to_length(param, n_groups, "group", verbose)

        ## create mapping:  event_id -> group_id
        ## rep(1:n_groups, lengths(group_intervals)) gives group index per event in group_intervals
        group_for_event <- rep(seq_len(n_groups), lengths(group_intervals))

        ## build lookup:  position i holds group index for event i (NA if ungrouped)
        event_to_group <- integer(n_events)
        valid_events <- groups_unlisted[groups_unlisted <= n_events]
        valid_groups <- group_for_event[groups_unlisted <= n_events]
        event_to_group[valid_events] <- valid_groups

        ## fill ungrouped (zero) positions with last group index
        event_to_group[event_to_group == 0L] <- n_groups

        ## index into param by group assignment
        return(param[event_to_group])
    }

    ## standard recycling for "distinct" or "ensemble"
    return(recycle_to_length(param, n_events, "event", verbose))
}


#' apply span to resolved times and build interval_spec data frame
#' @keywords internal
apply_span <- function(
    interval_list,
    time_vec,
    span,
    verbose = TRUE
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
    t_min <- time_vec[1L]
    t_max <- time_vec[length(time_vec)]

    ## check for entirely out of bounds intervals
    entirely_oob <- end_times < t_min | start_times > t_max
    if (any(entirely_oob)) {
        oob_ids <- which(entirely_oob)
        n_oob <- qty(length(oob_ids))
        cli_abort(c(
            "x" = "{n_oob} Interval{?s} {.val {oob_ids}} {n_oob} \\
            {?is/are} entirely outside data bounds.",
            "i" = "Intervals must be specified within existing data bounds."
        ))
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
                "!" = "{n_oob} Interval{?s} {.val {oob_ids}} {n_oob} \\
                {?is/are} partially outside data bounds.",
                "i" = "Returning available data only."
            ))
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
    time_vec,
    interval_spec,
    nirs_channels ## as list
) {
    n_vec <- seq_len(nrow(interval_spec))
    df_list <- lapply(n_vec, \(.i) {
        time_range <- time_vec >= interval_spec$start_times[.i] &
            time_vec <= interval_spec$end_times[.i]
        interval_data <- data[time_range, , drop = FALSE]

        ## return interval_data with metadata
        create_mnirs_data(
            interval_data,
            nirs_channels = nirs_channels[[.i]], ## overwrite for interval data
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
#' @keywords internal
ensemble_intervals <- function(
    df_list,
    nirs_channels,
    metadata,
    verbose = TRUE
) {
    time_channel <- metadata$time_channel
    sample_rate <- metadata$sample_rate
    ## extract data & metadata from df_list
    interval_data <- lapply(df_list, \(.df) {
        interval_times <- attr(.df, "interval_times")
        time_channel <- attr(.df, "time_channel")
        ## zero-offset to start time (first element of interval_times)
        t0 <- interval_times[[1L]]
        ## return data & metadata from each interval
        list(
            ## ensemble-average time values makes no sense, so return zero-offset
            data = zero_offset_data(.df, time_channel, t0),
            interval_times = interval_times - t0, ## zero_time recalc
            interval_span = attr(.df, "interval_span")
        )
    })

    ## stack interval data frames
    df_long <- do.call(rbind, lapply(interval_data, `[[`, "data"))
    ## resample times to nearest estimated sample rate for binned ensembling
    time_resampled <- round(df_long[[time_channel]] * sample_rate) / sample_rate
    ## split row indices by unique time
    time_groups <- split(seq_len(nrow(df_long)), time_resampled)
    unique_times <- as.numeric(names(time_groups))

    ## warn if any time samples have only one value, implying irregular
    ## samples and may result in alternating samples instead of ensemble-means
    if (verbose && any(lengths(time_groups) < 2L)) {
        cli_warn(c(
            "!" = "Duplicate or irregular {.arg time_channel} samples \\
            detected after ensemble-averaging.",
            "i" = "Re-sample with {.fn mnirs::resample_mnirs} before \\
            ensemble-averaging.",
            "i" = "Check your resulting data for inconsistent results."
        ))
    }

    col_n <- length(nirs_channels)
    ## nirs_channel-wise means per unique time matrix operation
    ## nirs_channel must be vectorised and exist in the interval_data
    channel_matrix <- as.matrix(df_long[, nirs_channels, drop = FALSE])
    result_matrix <- vapply(time_groups, \(.idx) {
        colMeans(channel_matrix[.idx, , drop = FALSE], na.rm = TRUE)
    }, numeric(col_n))

    ## vapply returns channels x times (or vector if 1 channel)
    ## coerce to times x channels data frame
    result_df <- if (col_n == 1L) {
        data.frame(setNames(list(result_matrix), nirs_channels))
    } else {
        as.data.frame(t(result_matrix))
    }
    result <- data.frame(setNames(list(unique_times), time_channel), result_df)

    ## return with metadata
    return(
        create_mnirs_data(
            result,
            nirs_device = metadata$nirs_device,
            nirs_channels = unique(nirs_channels),
            time_channel = time_channel,
            event_channel = metadata$event_channel,
            sample_rate = sample_rate,
            start_timestamp = metadata$start_timestamp,
            interval_times = lapply(interval_data, `[[`, "interval_times"),
            interval_span = lapply(interval_data, `[[`, "interval_span")
        )
    )
}

#' Zero-offset time values and add metadata
#' @keywords internal
preserve_metadata <- function(data, metadata, zero_time = FALSE) {
    if (zero_time) {
        ## Zero-offset a single interval's time and interval_times attr
        interval_times <- attr(data, "interval_times")
        t0 <- interval_times[[1L]]
        data <- zero_offset_data(data, metadata$time_channel, t0)
        attr(data, "interval_times") <- interval_times - t0
        data
    }

    return(
        create_mnirs_data(
            data,
            nirs_device = metadata$nirs_device,
            nirs_channels = unique(attr(data, "nirs_channels")),
            time_channel = metadata$time_channel,
            event_channel = metadata$event_channel,
            sample_rate = metadata$sample_rate,
            start_timestamp = metadata$start_timestamp,
            interval_times = attr(data, "interval_times"),
            interval_span = attr(data, "interval_span")
        )
    )
}


#' Apply grouping to intervals
#' @keywords internal
apply_interval_groups <- function(
    df_list,
    nirs_channels,
    metadata,
    group_intervals,
    zero_time = FALSE,
    verbose = TRUE
) {
    time_channel <- metadata$time_channel
    n_intervals <- length(df_list)

    ## return distinct intervals
    if (n_intervals == 1L || group_intervals[[1L]][1L] == "distinct") {
        return(lapply(df_list, \(.x) {
            preserve_metadata(.x, metadata, zero_time)
        }))
    }

    ## return ensembled intervals
    if (group_intervals[[1L]][1L] == "ensemble") {
        return(list(
            ensemble = ensemble_intervals(
                df_list,
                unique(unlist(nirs_channels)),
                metadata,
                verbose
            )
        ))
    }

    ## custom grouping ===================================
    grouped_ids <- unlist(group_intervals)
    ungrouped_ids <- setdiff(seq_len(n_intervals), grouped_ids)

    ## add ungrouped ids as individual groups, sort by position
    if (length(ungrouped_ids) > 0) {
        group_intervals <- c(group_intervals, as.list(ungrouped_ids))
        group_intervals <- group_intervals[
            order(vapply(group_intervals, median, numeric(1), na.rm = TRUE))
        ]
        if (verbose) {
            cli_inform(c(
                "!" = "Intervals not specified by {.arg group_intervals}.",
                "i" = "Ungrouped intervals included as discrete."
            ))
        }
    }

    ## warn for duplicated intervals across groups
    dupes <- grouped_ids[duplicated(grouped_ids)]
    if (verbose && length(dupes) > 0) {
        cli_warn(c(
            "!" = "Duplicates detected of {qty(length(dupes))} \\
            interval{?s} {.val {dupes}}.",
            "i" = "Re-specify {.arg group_intervals} to remove duplicates."
        ))
    }

    ## process each group
    result <- lapply(group_intervals, \(.g) {
        if (length(.g) == 1L) {
            return(preserve_metadata(df_list[[.g]], metadata, zero_time))
        }
        ## multi-interval group: ensemble-average
        group_nirs <- unique(unlist(nirs_channels[.g]))
        ensemble_intervals(df_list[.g], group_nirs, metadata, verbose)
    })

    names(result) <- vapply(group_intervals, \(.g) {
        paste0("interval_", paste(.g, collapse = "_"))
    }, character(1))

    return(result)
}
