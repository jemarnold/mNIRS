#' Extract intervals from *{mnirs}* data
#'
#' Extract intervals from *"mnirs"* time series data, specifying interval
#' start and end boundaries by time value, event label, lap number, or sample
#' index.
#'
#' @param data A data frame of class *"mnirs"* containing time series data and
#'   metadata, or a `list()` of such data frames.
#'
#' @param nirs_channels A character vector or a `list()` of character vectors
#'   of mNIRS channel names to operate on within each interval (see *Details*).
#'   Names must match column names in `data` exactly.
#'   - Must only be specified when `group_intervals` contains *"ensemble"*-
#'     averaged intervals. If `group_intervals = "distinct"` no channel
#'     processing occurs.
#'   - If `NULL` (default), channels are retrieved from *"mnirs"* metadata.
#'
#' @param event_channel An *optional* character string giving the name of an
#'   event/lap column. The column may contain character event labels or integer
#'   lap numbers.
#'   - Required when using [by_label()] or [by_lap()] for `start` or `end`.
#'   - Retrieved from metadata if not defined explicitly.
#'
#' @param sample_rate An *optional* numeric sample rate (Hz) used to bin time
#'   values for ensemble-averaging. If `NULL`, will be estimated from
#'   `time_channel` (see *Details*).
#'
#' @param group_intervals Either a character string or a `list()` of integer
#'   vectors specifying how to group intervals (see *Details*).
#'   \describe{
#'     \item{`"distinct"`}{The default. Extract each interval as an independent
#'     data frame.}
#'     \item{`"ensemble"`}{Ensemble-average each specified `nirs_channel` across
#'     all detected intervals, returning a single data frame.}
#'     \item{`list(c(1, 2), c(3, 4))`}{Ensemble-average each specified
#'     `nirs_channel` within each group and return one data frame per group.}
#'   }
#'
#' @param start Specifies where intervals begin. Either raw values -- numeric
#'   for time values, character for event labels, explicit integer (e.g. `2L`)
#'   for lap numbers -- or created with [by_time()], [by_label()], [by_lap()],
#'   or [by_sample()].
#'
#' @param end Specifies where intervals end. Either raw values -- numeric for
#'   time values, character for event labels, explicit integer (e.g. `2L`)
#'   for lap numbers -- or created with [by_time()], [by_label()], [by_lap()],
#'   or [by_sample()].
#'
#' @param span A one- or two-element numeric vector `c(before, after)` in units
#'   of `time_channel`, or a `list()` of such vectors. (*default*
#'   `span = c(-60, 60)`. Applied additively to interval boundaries:
#'   - When both `start` and `end` are specified: `span[1]` shifts start times,
#'     `span[2]` shifts end times.
#'   - When only `start` or only `end` is specified: both `span[1]` and
#'     `span[2]` apply as a window around the event).
#'   - A single *positive* value is recycled to shift the end times (e.g.
#'     `span = 60` -> `c(0, 60)`).
#'   - A single *negative* value is recycled to shift the start times (e.g.
#'     `span = -60` -> `c(-60, 0)`).
#'
#' @param zero_time Logical. Default is `FALSE`. If `TRUE`, re-calculates
#'   numeric `time_channel` values to start from zero within each interval
#'   data frame.
#'
#' @param event_groups `r lifecycle::badge("deprecated")` Renamed to
#'   `group_intervals` for naming consistency across the package.
#'
#' @inheritParams validate_mnirs
#'
#' @details
#' ## Interval specification
#'
#' Interval `start` and `end` boundaries are specified using helper functions,
#' or by passing raw values directly:
#'
#' \describe{
#'   \item{[by_time()]}{Time values in units of `time_channel`.}
#'   \item{[by_label()]}{Strings to match in `event_channel`.
#'   All matching occurrences are returned.}
#'   \item{[by_lap()]}{Lap numbers to match in `event_channel`. Resolves to the
#'   first sample of each lap for `start`, and the last lap sample for `end`}
#'   \item{[by_sample()]}{Integer sample indices (row numbers).}
#' }
#'
#' Raw values supplied to `start`/`end` are auto-coerced:
#'   - Numeric -> [by_time()]
#'   - Character -> [by_label()],
#'   - Explicit integer (e.g. `2L`) -> [by_lap()].
#'   - Use [by_sample()] explicitly for sample indices.
#'
#' `start` and `end` can use different specification types (e.g., start by
#' label, end by time). When lengths differ, the shorter is recycled.
#'
#' ## Time span window
#'
#' `span` additively expands the time span window around interval boundaries.
#'
#' - A two-value vector expands the `start` and `end`, respectively:
#'   `span = c(-60, 60)` expands the `start` earlier by `60`, and the `end`
#'   later by `60`. For example,
#'   `start = by_time(30), end = by_time(60), span = c(-5, 10)` returns an
#'   interval of `[25, 70]`.
#' - A single numeric value is recycled according to the sign: `span = -60`
#'   becomes `c(-60, 0)` to expand the `start` earlier. `span = 60` becomes
#'   `c(0, 60)` to expand the `end` later.
#' - If only `start` is specified alone, both span values expand the single
#'   boundary window: `start = by_time(30), span = c(-5, 60)` returns
#'   `[25, 90]`.
#'
#' ## Per-interval nirs_channels for ensemble-averaging
#'
#' When `group_intervals = "ensemble"` or a list of numeric grouped intervals,
#' `nirs_channels` can be specified as a list of column names to override
#' ensemble-averaging across interval. For example, to exclude a channel
#' in one interval:
#'
#' ```r
#' nirs_channels = list(
#'   c(A, B, C),
#'   c(A, C) ## channel "B" is excluded
#' )
#' ```
#'
#' If all grouped intervals can include all `nirs_channels`, or if
#' `group_intervals = "distinct"`, a single `nirs_channels` character vector can
#' be supplied and recycled to all groups, or left as `NULL` for channels to
#' be taken from *"mnirs"* metadata.
#'
#' ## Grouping intervals
#'
#' `group_intervals` controls whether extracted intervals are returned as distinct
#' data frames or ensemble-averaged.
#'
#' \describe{
#'    \item{`"distinct"`}{The default. Extract each interval and return a
#'    list of independent data frames.}
#'    \item{`"ensemble"`}{Ensemble-average each specified `nirs_channel` across
#'    all detected intervals and return a one-item list with a single data
#'    frame.}
#'    \item{`list(c(1, 2), c(3, 4))`}{Ensemble-average each specified
#'    `nirs_channel` within each group and return a list with one data frame
#'    for each group. Any intervals detected but not specified in
#'    `group_intervals` are returned as distinct.}
#' }
#'
#' `group_intervals` lists can be named (e.g.
#' `list(low = c(1, 2), high = c(3, 4))`) and will pass those names to the
#' returned list of data frames.
#'
#' When `group_intervals` is a list of numeric interval numbers, list items in
#' `nirs_channels` and `span` are recycled to the number of groups. If lists
#' are only partially specified, the final item is recycled forward as needed.
#' Extra items are ignored.
#'
#' @returns A named `list()` of [tibbles][tibble::tibble-package] of class
#'   *"mnirs"*, each with metadata available via `attributes()`. When `data`
#'   is a list of data frames, results are flattened into a single-layer
#'   list with interval names indicating nested number of data frame and
#'   interval as `interval_<df>.<interval>`. Other interval names (e.g.
#'   `"ensemble"` or custom group names) are suffixed as `<name>_<df>`.
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
#'     zero_time = TRUE,
#'     verbose = FALSE
#' ) |>
#'     ## avoid issues ensemble-averaging irregular samples
#'     resample_mnirs(method = "linear", verbose = FALSE)
#'
#' ## ensemble-average across multiple intervals
#' interval_list <- extract_intervals(
#'     data,                         ## channels recycled to all intervals by default
#'     nirs_channels = c(smo2_left, smo2_right),
#'     group_intervals = "ensemble", ## ensemble-average across two intervals
#'     start = by_time(368, 1084),   ## manually identified interval start times
#'     span = c(-20, 90),            ## include the last 180-sec of each interval (recycled)
#'     zero_time = TRUE              ## re-calculate common time to start from `0`
#' )
#'
#' interval_list[[1L]]
#'
#' \donttest{
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot(interval_list, time_labels = TRUE) +
#'       ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
#'   }
#' }
#'
#' @export
extract_intervals <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    event_channel = NULL,
    sample_rate = NULL,
    group_intervals = c("distinct", "ensemble"),
    start = NULL,
    end = NULL,
    span = list(c(-60, 60)),
    zero_time = FALSE,
    verbose = TRUE,
    event_groups = deprecated()
) {
    ## list input → recurse per df, flatten nested interval lists
    if (!is.data.frame(data)) {
        nested <- map_mnirs_intervals(data, match.call(), parent.frame())
        ## rename intervals `interval_<df>.<interval>` before flattening;
        ## other names (e.g. "ensemble", <custom>) get `<name>_<df>` suffix
        nested <- lapply(seq_along(nested), \(.i) {
            .x <- nested[[.i]]
            names(.x) <- ifelse(
                startsWith(names(.x), "interval_"),
                sub("^interval_", paste0("interval_", .i, "."), names(.x)),
                paste0(names(.x), "_", .i)
            )
            .x
        })
        result <- unlist(nested, recursive = FALSE)
        class(result) <- c("mnirs", class(result))
        return(result)
    }

    ## validation ==============================================
    if (lifecycle::is_present(event_groups)) {
        lifecycle::deprecate_stop(
            when = "0.7.0",
            what = "extract_intervals(event_groups)",
            with = "extract_intervals(group_intervals)"
        )
    }
    validate_mnirs_data(data)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    nirs_channels <- validate_nirs_channels(
        enquo(nirs_channels), data, verbose, as_list = TRUE
    )
    time_channel <- validate_time_channel(enquo(time_channel), data)
    ## avoid floating point precision issues downstream with findInterval()
    t_vec <- round(data[[time_channel]], 6)
    ## estimate sample_rate for appropriate time binning
    sample_rate <- validate_sample_rate(
        data, time_channel, sample_rate, verbose
    )

    ## coerce raw values to mnirs_interval objects
    start <- as_mnirs_interval(start, "start")
    end <- as_mnirs_interval(end, "end")

    ## validate interval specs
    if (is.null(c(start, end))) {
        cli_abort(c(
            "x" = "No interval specification provided.",
            "i" = "Specify {.arg start} and/or {.arg end} using \\
            {.fn by_time}, {.fn by_sample}, {.fn by_label}, or {.fn by_lap}."
        ))
    }

    ## resolve event_channel if by_label or by_lap is used
    uses_event_channel <- (inherits(start, "mnirs_interval") &&
        start$type %in% c("label", "lap")) ||
        (inherits(end, "mnirs_interval") && end$type %in% c("label", "lap"))

    ## report conditions raised in handlers/lambdas from this function
    env <- environment()
    if (uses_event_channel) {
        event_channel <- tryCatch(
            validate_event_channel(event_channel, data, required = TRUE),
            error = function(e) {
                cli_abort(c(
                    "x" = "{.arg event_channel} is required when using \\
                    {.fn by_label} or {.fn by_lap}.",
                    "i" = "Specify column name containing {.cls character} \\
                    event labels or {.cls integer} lap numbers."
                ), call = env)
            }
        )
        event_vec <- data[[event_channel]]
    } else {
        event_channel <- validate_event_channel(
            event_channel, data, required = FALSE
        )
        event_vec <- NULL
    }

    if (
        verbose &&
            group_intervals[1L] != "distinct" &&
            is.null(attr(data, "nirs_channels")) &&
            !is.list(nirs_channels)
    ) {
        cli_inform(c(
            "!" = "{.fn extract_intervals} accepts {.arg nirs_channels} = \\
            {col_blue('list()')} for channel grouping when \\
            ensemble-averaging. See `?extract_intervals`."
        ))
    }

    ## expand parameters ====================================
    span <- if (is.list(span)) span else list(span)
    span <- lapply(span, recycle_span, env = env)
    group_intervals <- if (is.list(group_intervals)) {
        group_intervals
    } else {
        list(group_intervals)
    }

    ## resolve interval boundary times ==========================
    interval_list <- resolve_interval(start, end, t_vec, event_vec)
    n_events <- length(interval_list$start_time)

    ## recycle params to match number of intervals
    nirs_channels <- recycle_param(
        nirs_channels,
        n_events,
        group_intervals,
        verbose
    )
    span <- recycle_param(span, n_events, group_intervals, verbose)

    ## apply span and build interval spec ======================
    interval_spec <- apply_span(interval_list, t_vec, span, verbose)

    ## extract interval data ===================================
    df_list <- extract_df_list(data, t_vec, interval_spec, nirs_channels)

    ## apply grouping logic ====================================
    result <- apply_interval_groups(
        df_list,
        nirs_channels,
        metadata,
        group_intervals,
        zero_time,
        verbose
    )

    ## add class "mnirs" ========================================
    class(result) <- c("mnirs", class(result))

    return(result)
}
