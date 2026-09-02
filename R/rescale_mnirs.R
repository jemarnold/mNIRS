#' Rescale data range
#'
#' Expand or reduce the range (min and max values) of data channels to a new
#' amplitude/dynamic range, e.g. rescale the range of NIRS data to `c(0, 100)`.
#'
#' @usage
#' rescale_mnirs(
#'   data,
#'   nirs_channels = NULL,
#'   group_channels = c("ensemble", "distinct"),
#'   range,
#'   verbose = TRUE
#' )
#'
#' @param group_channels Either a character string or a `list()` of
#'   channel-name vectors specifying how to group `nirs_channels`
#'   (see *Details*).
#'   \describe{
#'      \item{`"ensemble"`}{The *default*. Operate on all channels together,
#'      preserving the relative scaling between channels.}
#'      \item{`"distinct"`}{Operate on each channel independently, losing
#'      the relative scaling between channels.}
#'      \item{`list(c("A", "B"), c("C", "D"))`}{Operate on channels `A` & `B`
#'      in one group, and `C` & `D` in another group. Groups can be named
#'      (e.g. `list(smo2 = c("A", "B"))`). Each group must be non-empty
#'      and resulting group names must be unique.}
#'   }
#' @param range A numeric vector in the form `c(min, max)`, indicating the
#'   range of output values to which `nirs_channels` will be rescaled.
#' @inheritParams map_mnirs_intervals
#' @inheritParams validate_mnirs
#'
#' @inheritSection map_mnirs_intervals Data input formats
#'
#' @details
#' `group_channels` controls how data channels are grouped to preserve
#'   absolute or relative scaling.
#'
#' - `group_channels = "ensemble"` (the *default*) rescales all
#'   `nirs_channels` to a common range, preserving relative scaling
#'   between channels.
#'
#' - `group_channels = "distinct"` rescales each channel independently,
#'   losing relative scaling between channels.
#'
#' - A `list()` of channel-name vectors (e.g. `list(c("A", "B"), c("C", "D"))`)
#'   rescales channels `A` & `B` together and `C` & `D` together, preserving
#'   relative scaling within, but not between groups. `nirs_channels`
#'   omitted from the list are rescaled independently.
#'
#' - Channel groups can be named (e.g. `list(smo2 = c("A", "B"))`) and names
#'   used as keys for per-group `range` argument.
#'
#' - Channels (columns) in `data` not in `nirs_channels` are passed
#'   through without processing to the output data frame.
#'
#' `nirs_channels` can be retrieved automatically from `data` of class
#'   *"mnirs"* which has been processed with `{mnirs}`, if not defined
#'   explicitly.
#'
#' @inheritSection shift_mnirs Per-channel arguments
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
#'     rescale_mnirs(        ## un-grouped nirs channels to rescale separately
#'         nirs_channels = c(smo2_left, smo2_right),
#'         group_channels = "distinct",
#'         range = c(0, 100)  ## rescale to a 0-100% functional exercise range
#'     )
#'
#' data
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         plot(data, time_labels = TRUE) +
#'             ggplot2::geom_hline(yintercept = c(0, 100), linetype = "dotted")
#'     }
#' }
#'
#' @export
rescale_mnirs <- function(
    data,
    nirs_channels = NULL,
    group_channels = c("ensemble", "distinct"),
    range,
    verbose = TRUE
) {
    ## list or grouped input -> normalise to named list, recurse per interval
    if (inherits(data, "grouped_df") || !is.data.frame(data)) {
        return(map_mnirs_intervals(data, match.call(), parent.frame()))
    }

    ## validate =================================
    validate_mnirs_data(data, ncol = 1)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    nirs_parsed <- parse_channel_name(enquo(nirs_channels), data)
    if (is.list(nirs_parsed)) {
        lifecycle::deprecate_stop(
            when = "0.7.0",
            what = I(paste(
                "Passing a `list()` to `nirs_channels` for",
                "channel grouping"
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

    ## broadcast global args, applying per-channel/per-group overrides
    per_group <- resolve_channel_args(
        nirs_channels,
        args = list(range = range),
        group_channels = group_channels,
        verbose = verbose
    )

    ## rescale range per group ================================
    ## report conditions raised in the lambda from this function
    env <- environment()
    rescaled <- lapply(names(group_channels), \(.key) {
        .range <- per_group[[.key]]$range
        .cols <- group_channels[[.key]]
        validate_numeric(
            .range, 2,
            msg1 = "two-element",
            msg2 = "between {col_blue('range[1], range[2]]')}.",
            env = env
        )

        group_data <- as.matrix(data[, .cols, drop = FALSE])
        group_data_range <- range(group_data, na.rm = TRUE)
        range_diff <- diff(group_data_range)

        if (range_diff != 0) {
            group_data <- (group_data - group_data_range[1]) /
                range_diff * diff(.range) + .range[1]
        }
        group_data
    })

    ## apply re-scaling in group order
    data[unlist(group_channels, use.names = FALSE)] <- do.call(cbind, rescaled)

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)

    return(create_mnirs_data(data, metadata))
}
