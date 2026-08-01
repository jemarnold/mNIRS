#' Shift data range
#'
#' Move the range of data channels in a data frame up or down, while preserving
#' the absolute amplitude/dynamic range of each channel, and the relative
#' scaling across channels. e.g. shift the minimum data value to zero for all
#' positive values, or shift the mean of the first time span in a recording
#' to zero.
#'
#' @param to A numeric value in units of `nirs_channels` to which the data
#'   channels will be shifted, e.g. shift the minimum value to zero.
#' @param by A numeric value in units of `nirs_channels` by which the data
#'   channels will be shifted, e.g. shift all values up by 10 units.
#' @param position Indicates where the reference values will be shifted from.
#'   \describe{
#'      \item{`"min"`}{(The *default*) will shift the minimum value(s) `to`
#'      or `by` the specified value.}
#'      \item{`"max"`}{Will shift the maximum value(s) `to` or `by` the
#'      specified values.}
#'      \item{`"first"`}{Will shift first value(s) `to` or `by` the specified
#'      values.}
#'   }
#' @inheritParams map_mnirs_intervals
#' @inheritParams validate_mnirs
#' @inheritParams replace_mnirs
#' @inheritParams rescale_mnirs
#'
#' @inheritSection map_mnirs_intervals Data input formats
#'
#' @details
#' `group_channels` controls how data channels are grouped to preserve
#'   absolute or relative scaling (see [rescale_mnirs()]).
#'
#' - `group_channels = "ensemble"` (the *default*) shifts all
#'   `nirs_channels` to a common value, preserving relative scaling
#'   between channels.
#'
#' - `group_channels = "distinct"` shifts each channel independently,
#'   losing relative scaling between channels.
#'
#' - A `list()` of channel-name vectors (e.g. `list(c("A", "B"), c("C", "D"))`)
#'   shifts channels `A` & `B` together and `C` & `D` together, preserving
#'   relative scaling within, but not between groups. `nirs_channels`
#'   omitted from the list are rescaled independently.
#'
#' - Channel groups can be named (e.g. `list(smo2 = c("A", "B"))`) and names
#'   used as keys for per-group arguments.
#'
#' Only one of either `to` or `by` and one of either `width` or `span` should
#'   be defined for each `group_channels`. If both of either pairing are
#'   defined, `to` will be preferred over `by`, and `width` will be preferred
#'   over `span`.
#'
#' - Channels (columns) in `data` not in `nirs_channels` are passed
#'   through without processing to the output data frame.
#'
#' `nirs_channels` and `time_channel` can be retrieved automatically from
#'   `data` of class *"mnirs"* which has been processed with `{mnirs}`,
#'   if not defined explicitly.
#' 
#' When `position` is *"min"* or *"max"*, only full windows of `width` or 
#'   `span` are considered, to avoid bias from noise at edge conditions with
#'   partial samples.
#'
#' @section Per-channel arguments:
#'
#' Arguments apply globally to all `nirs_channels` by default. Relevant
#' arguments can instead be supplied uniquely per-channel as a named `list()`,
#' with names matching either `nirs_channels` or list names in
#' `group_channels`, e.g.:
#'
#' ```r
#' shift_mnirs(
#'     data,
#'     nirs_channels = c(A, B, C),
#'     group_channels = list(smo2 = c(A, B), hhb = C),
#'     to = list(100, C = 0),
#'     width = list(smo2 = 3),
#'     span = list(hhb = 5),
#'     position = "first"
#' )
#' ```
#'
#' - A non-list value applies to every channel (the *default* behaviour).
#' - A `list()` named by `nirs_channels` or `group_channels` applies
#'   per-channel / per-group values.
#' - A single unnamed value in the list will be applied to unlisted channels
#'   (e.g. `span = list(3, hhb = 5)` gives `hhb` 5 and every other channel 3).
#'   If no unnamed fallback value in the list, channels not named in the list
#'   will be returned un-processed (e.g. `span = list(hhb = 5)` will only
#'   process `hhb`).
#' - `list()` names not matching `nirs_channels` or `group_channels` are
#'   warned about and ignored.
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class *"mnirs"* with metadata
#'   available with `attributes()`. For list or grouped data frame input,
#'   returns a named list of *"mnirs"* tibbles, one per interval.
#'
#' @examples
#' ## read example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2_left = "SmO2 Live",
#'                       smo2_right = "SmO2 Live(2)"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' ) |>
#'     shift_mnirs(        ## un-grouped nirs channels to shift separately
#'         nirs_channels = c(smo2_left, smo2_right),
#'         group_channels = "distinct",
#'         to = 0,         ## NIRS values will be shifted to zero
#'         span = 120,     ## shift the *first* 120 sec of data to zero
#'         position = "first"
#'     )
#'
#' data
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         plot(data, time_labels = TRUE) +
#'             ggplot2::geom_hline(yintercept = 0, linetype = "dotted")
#'     }
#' }
#'
#' @export
shift_mnirs <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    group_channels = c("ensemble", "distinct"),
    to = NULL,
    by = NULL,
    width = NULL,
    span = NULL,
    position = c("min", "max", "first"),
    verbose = TRUE
) {
    ## list or grouped input → normalise to named list, recurse per interval
    if (inherits(data, "grouped_df") || !is.data.frame(data)) {
        return(map_mnirs_intervals(data, match.call(), parent.frame()))
    }

    ## validation =============================================
    validate_mnirs_data(data)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    nirs_parsed <- parse_channel_name(enquo(nirs_channels), data)
    if (is.list(nirs_parsed)) {
        lifecycle::deprecate_stop(
            when = "0.7.0",
            what = I(paste(
                "Passing a `list()` to `nirs_channels` for channel grouping"
            )),
            with = I("the `group_channels` argument")
        )
    }
    nirs_channels <- validate_nirs_channels(nirs_parsed, data)
    group_channels <- validate_group_channels(
        nirs_channels,
        enquo(group_channels),
        data
    )
    time_channel <- validate_time_channel(enquo(time_channel), data)
    t_vec <- data[[time_channel]]

    ## broadcast global args, applying per-channel/per-group overrides
    per_group <- resolve_channel_args(
        nirs_channels,
        args = list(
            to = to,
            by = by,
            width = width,
            span = span,
            position = position
        ),
        defaults = list(position = "min"),
        choices = list(position = c("min", "max", "first")),
        group_channels = group_channels,
        verbose = verbose
    )

    ## per-group shifts ==============================================
    env <- environment()
    shifted <- lapply(names(group_channels), \(.key) {
        .a <- per_group[[.key]]
        .cols <- group_channels[[.key]]
        ## verbose messages emitted once for the first channel/group
        .v <- verbose && .key == names(group_channels)[[1L]]

        ## do nothing condition
        if (is.null(c(.a$to, .a$by))) {
            cli_abort(c(
                "x" = "Shift value undefined",
                "i" = "One of {.arg to} or {.arg by} must be defined."
            ), call = env)
        }
        validate_numeric(.a$to, 1, msg1 = "one-element", env = env)
        validate_numeric(.a$by, 1, msg1 = "one-element", env = env)
        if (!is.null(.a$to) && !is.null(.a$by)) {
            .a$by <- NULL
            if (.v) {
                cli_inform(c(
                    "i" = "Shift {.arg to} = {.val {(.a$to)}} \\
                    overrides {.arg by}."
                ), call = env)
            }
        }

        ## shift_by does not require calculating reference positions
        if (!is.null(.a$by)) {
            return(data[.cols] + .a$by)
        }

        ## calculate shift_to reference values =======================
        validate_width_span(
            .a$width, .a$span, .v, "for `shift_mnirs()`.", env = env
        )

        if (.a$position == "first") {
            ## for span, take data <= first time value + span, assuming
            ## sorted time_channel (sum == last index)
            first_width <- .a$width %||% sum(t_vec <= t_vec[1L] + .a$span)
            ## drop = FALSE to avoid reducing single channels to vector
            shift_values <- colMeans(
                data[seq_len(first_width), .cols, drop = FALSE],
                na.rm = TRUE
            )
        } else {
            ## find local windows within width/span centred around idx
            window_idx <- compute_local_windows(
                t = t_vec, width = .a$width, span = .a$span, env = env
            )

            ## exclude partial windows to avoid bias on noise
            ## min_obs default to estimated width when span is specified
            ## less strict span_width - 2 to allow start & end buffer
            ## with irregular t values
            min_obs <- max(
                .a$width %||%
                    (floor(.a$span * estimate_sample_rate(t_vec, env)) - 2L),
                1L
            )
            window_idx <- window_idx[lengths(window_idx) >= min_obs]

            ## no complete window: reference would be Inf/-Inf
            if (length(window_idx) == 0L) {
                cli_abort(c(
                    "x" = "Insufficient valid samples detected.",
                    "i" = "{.arg width} or {.arg span} must be smaller than \\
                    the range of {.arg data}."
                ), call = env)
            }

            shift_fun <- match.fun(.a$position)
            ## compute min or max along local means per channel
            shift_values <- vapply(data[.cols], \(.x) {
                shift_fun(
                    compute_local_fun(.x, window_idx, mean, na.rm = TRUE),
                    na.rm = TRUE
                )
            }, numeric(1))
        }

        ## single reference value shared across the channel group
        group_shift <- switch(
            .a$position,
            min = min(shift_values),
            max = max(shift_values),
            first = mean(shift_values)
        )
        data[.cols] - group_shift + .a$to
    })

    ## apply shifts in group order =====================================
    data[unlist(group_channels, use.names = FALSE)] <- do.call(cbind, shifted)

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)
    metadata$time_channel <- time_channel

    return(create_mnirs_data(data, metadata))
}
