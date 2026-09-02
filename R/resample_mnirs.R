#' Re-sample an *{mnirs}* data frame
#'
#' Up- or down-sample an *"mnirs"* data frame to a new sample rate, filling
#' new samples via nearest-neighbour matching or interpolation.
#'
#' @param resample_rate An *optional* sample rate (Hz) for the output data
#'   frame. If `NULL` (*default*) resamples to the existing `sample_rate`,
#'  which regularises any irregular samples without changing the rate.
#'
#' @param method A character string specifying how new samples are filled.
#'   Default is *"none"*. Filling must be opted into explicitly (see *Details*):
#'   \describe{
#'     \item{`"none"`}{Matches each new sample to the nearest original
#'     `time_channel` value without any interpolation, to within tolerance of
#'     half a sample-interval. New samples are returned as `NA`.}
#'     \item{`"locf"`}{(*"Last observation carried forward"*). Fills new and
#'     missing samples with the most recent valid non-`NA` value to the left,
#'     or the nearest valid value to the right for leading `NA`s. Safe for
#'     numeric, integer, and character columns.}
#'     \item{`"linear"`}{Fills new and missing samples via linear interpolation
#'     using [stats::approx()]. Suitable for numeric columns only; non-numeric
#'     columns will fall back to `"locf"` behaviour.}
#'   }
#'
#' @inheritParams map_mnirs_intervals
#' @inheritParams validate_mnirs
#'
#' @inheritSection map_mnirs_intervals Data input formats
#'
#' @details
#' This function uses [replace_missing()] (based on [stats::approx()]) to
#' interpolate across new samples in the resampled data range.
#'
#' ## Sample rate and time channel
#'
#' `time_channel` and `sample_rate` are retrieved automatically from `data`
#' of class *"mnirs"*, if not defined explicitly.
#'
#' Otherwise, `sample_rate` will be estimated from the values in `time_channel`.
#' However, this may return unexpected values, and it is safer to define
#' `sample_rate` explicitly or retrieve it from *"mnirs"* metadata.
#'
#' ## Default behaviour
#'
#' When `resample_rate` is omitted, the output has the same `sample_rate` as
#' the input but with a regular, evenly-spaced `time_channel`. This is useful
#' for regularising data that contains missing or repeated samples without
#' changing the nominal rate.
#'
#' ## Column handling
#'
#' Numeric columns are interpolated according to `method` (see
#' `?replace_missing`). Non-numeric columns (e.g. character event labels,
#' integer lap numbers) are always filled by last-observation-carried-forward,
#' regardless of `method`:
#'
#' - For `method = "none"`, existing rows are matched to the nearest original
#'   values of `time_channel` without interpolation or filling, meaning newly
#'   created samples and any `NA`s in the original data are returned as `NA`.
#' - When down-sampling, numeric columns use linear interpolation averaging.
#'   Non-numeric columns use the first valid value in each output bin.
#'
#' @returns A [tibble][tibble::tibble-package] of class `"mnirs"`. Metadata are
#'   stored as attributes and can be accessed with `attributes(data)`. For
#'   list or grouped data frame input, returns a named list of *"mnirs"*
#'   tibbles, one per interval.
#'
#' @examples
#' ## read example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2 = "SmO2 Live"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = TRUE
#' )
#'
#' ## note warning about irregular sampling
#' data
#'
#' data_resampled <- resample_mnirs(
#'     data,               ## blank channels will be retrieved from metadata
#'     resample_rate = 2,  ## blank by default will resample to `sample_rate`
#'     method = "linear",  ## linear interpolation across resampled indices
#'     verbose = TRUE
#' )
#'
#' ## note the altered `time` values resolving the above warning
#' data_resampled
#'
#' @export
resample_mnirs <- function(
    data,
    time_channel = NULL,
    sample_rate = NULL,
    resample_rate = sample_rate, ## placeholder indicating default condition
    method = c("none", "linear", "locf"),
    verbose = TRUE
) {
    ## list or grouped input -> normalise to named list, recurse per interval
    if (inherits(data, "grouped_df") || !is.data.frame(data)) {
        return(map_mnirs_intervals(data, match.call(), parent.frame()))
    }

    ## validation ====================================
    validate_mnirs_data(data)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    time_channel <- validate_time_channel(enquo(time_channel), data)
    sample_rate <- validate_sample_rate(
        data, time_channel, sample_rate, verbose
    )
    validate_numeric(
        resample_rate, 1, c(0, Inf), FALSE, msg1 = "one-element positive"
    )
    ## assign default resample_rate as sample_rate
    if (is.null(resample_rate)) {
        resample_rate <- sample_rate
    }
    method <- match.arg(method)

    ## calculate resampling parameters  ===================================
    resample_time <- 1 / resample_rate
    cols <- names(data)
    t_vec <- data[[time_channel]]
    n_vec <- length(t_vec)
    t_range <- range(t_vec, na.rm = TRUE) * resample_rate
    t_range <- c(floor(t_range[1L]), round(t_range[2L])) / resample_rate
    t_out <- seq(t_range[1L], t_range[2L], by = resample_time)
    n_out <- length(t_out)

    ## un-tibble cheap list view of data
    data_list <- unclass(data)

    ## identify numeric cols; integers and `time_channel` excluded
    numeric_cols <- vapply(data_list, \(.x) {
        is.numeric(.x) && !is.integer(.x)
    }, logical(1))
    numeric_cols[time_channel] <- FALSE
    numeric_names <- cols[numeric_cols]
    non_numeric_names <- cols[!numeric_cols & cols != time_channel]

    ## preallocate output column list in input order
    out <- vector("list", length(cols))
    names(out) <- cols
    out[[time_channel]] <- t_out

    ## nearest-match index for "none" method =============================
    if (method == "none") {
        ## bracket each resampled time within sorted t_vec, then pick
        ## the closer of the left/right neighbours. `left.open = TRUE`
        ## ensures that on duplicate t_vec values the *first* duplicate
        ## wins (matches the previous "first match" tie-breaking).
        left <- validate_findInt(
            t_out, t_vec, all.inside = TRUE, left.open = TRUE
        )
        right <- pmin(left + 1L, n_vec)
        ## branchless nearest-pick; left wins on tie
        idx <- left
        pick_right <- t_out - t_vec[left] > t_vec[right] - t_out
        idx[pick_right] <- right[pick_right]
        ## tol = half resample_time
        match_ok <- abs(t_out - t_vec[idx]) < resample_time * 0.5
        idx[!match_ok] <- NA
    }

    ## interpolate numeric columns (locf/linear only) ==================
    if (method != "none" && length(numeric_names) > 0) {
        approx_method <- if (method == "locf") "constant" else method
        out[numeric_names] <- lapply(data_list[numeric_names], \(.x) {
            stats::approx(
                x = t_vec,
                y = .x,
                xout = t_out,
                method = approx_method,
                rule = 2,
                f = 0,
                ties = list("ordered", mean)
            )$y
        })
    }

    ## index-based fill for non-numeric columns ===========================
    idx_names <- if (method == "none") {
        cols[cols != time_channel]
    } else {
        non_numeric_names
    }

    if (length(idx_names) > 0) {
        if (method == "none") {
            ## tolerance-matched index computed above
            out[idx_names] <- lapply(data_list[idx_names], \(.x) .x[idx])
        } else if (n_out < n_vec) {
            ## down-sample: assign each original sample to its output bin,
            ## take first non-NA value falling in each bin (or NA if all NA)
            bin <- pmin(pmax(validate_findInt(t_vec, t_out), 1L), n_out)
            out[idx_names] <- lapply(data_list[idx_names], \(.x) {
                ok <- which(!is.na(.x))
                .x[ok[match(seq_len(n_out), bin[ok])]]
            })
        } else {
            ## up-sample / regularise: assign each output sample to
            ## most recent original sample at or before current
            idx <- validate_findInt(t_out, t_vec, rightmost.closed = FALSE)
            idx[idx == 0L] <- 1L
            out[idx_names] <- lapply(data_list[idx_names], \(.x) .x[idx])
        }
    }

    ## assemble data.frame from completed column list
    result <- structure(
        out,
        class = "data.frame",
        row.names = .set_row_names(n_out)
    )

    ## Metadata =================================
    metadata$time_channel <- time_channel
    metadata$sample_rate <- resample_rate

    if (verbose) {
        cli_inform(c("i" = "Output is resampled at {.val {resample_rate}} Hz."))
    }

    return(create_mnirs_data(result, metadata))
}
