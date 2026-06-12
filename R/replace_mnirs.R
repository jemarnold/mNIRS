#' Replace outliers, invalid, and missing values in *{mnirs}* data
#'
#' Detect and replace local outliers, specified invalid values, and missing
#' `NA` values across `nirs_channels` within an *"mnirs"* data frame.
#' `replace_mnirs()` operates on a data frame, extending the vectorised 
#' functions:.
#'
#' @param invalid_values A numeric vector of invalid values to be replaced,
#'   e.g. `invalid_values = c(0, 100, 102.3)`. Default `NULL` will not replace
#'   invalid values.
#'
#' @param invalid_above,invalid_below Numeric values each specifying cutoff
#'   values, above or below which (respectively) will be replaced, *inclusive*
#'   of the specified cutoff values.
#'
#' @param outlier_cutoff A numeric value for the local outlier threshold, as
#'   the number of standard deviations from the local median.
#'   - Default `NULL` will not replace outliers.
#'   - Lower values are more sensitive and flag more outliers; higher values
#'     are more conservative.
#'   - `outlier_cutoff = 3` Pearson's 3 sigma edit rule.
#'     `outlier_cutoff = 2` approximates a Tukey-style 1.5*IQR rule.
#'     `outlier_cutoff = 0` Tukey's median filter.
#'
#' @param width An integer defining the local window in number of samples
#'   centred on `idx`, between
#'   `[idx - floor(width/2), idx + floor(width/2)]`.
#'
#' @param span A numeric value defining the local window time span around
#'   `idx` in units of `time_channel` or `t`, between
#'   `[t - span/2, t + span/2]`.
#'
#' @param method A character string indicating how to handle `NA`
#'   replacement (see *Details*):
#'   \describe{
#'     \item{`"linear"`}{Replaces `NA`s via linear interpolation (the
#'     *default*) using [stats::approx()].}
#'     \item{`"median"`}{Replaces `NA`s with the local median of valid
#'     values within a centred window defined by `width` or `span`.}
#'     \item{`"locf"`}{*"Last observation carried forward"*. Replaces
#'     `NA`s with the most recent valid value to the left for trailing
#'     `NA`s or to the right for leading `NA`s, using
#'     [stats::approx()].}
#'     \item{`"none"`}{Returns `NA`s without replacement.}
#'   }
#'
#' @inheritParams validate_mnirs
#'
#' @details
#' ## Automatic channel detection
#'
#' `nirs_channels` and `time_channel` are retrieved automatically from
#' *"mnirs"* metadata if not specified explicitly. Columns in `data` not
#' listed in `nirs_channels` are passed through unprocessed.
#'
#' ## The rolling window
#'
#' `replace_outliers()` and `replace_missing()` (when `method = "median"`)
#' operate over a local rolling window for outlier detection and median
#' interpolation. The window is specified by either `width` as the number
#' of samples, or `span` as the time span in units of `time_channel`.
#' A partial window is calculated at the edges of the data.
#'
#' @returns `replace_mnirs()` return a [tibble][tibble::tibble-package] of 
#' class `"mnirs"` with metadata available via `attributes()`.
#'
#' @examples
#' ## vectorised operations
#' x <- c(1, 999, 3, 4, 999, 6)
#' replace_invalid(x, invalid_values = 999, width = 3, method = "median")
#'
#' (x_na <- replace_outliers(x, outlier_cutoff = 3, width = 3, method = "none"))
#'
#' replace_missing(x_na, method = "linear")
#'
#' ## read example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2 = "SmO2 Live"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' )
#'
#' ## clean data
#' data_clean <- replace_mnirs(
#'     data,                  ## channels retrieved from metadata
#'     invalid_values = 0,    ## known invalid values in the data
#'     invalid_above = 90,    ## remove data spikes above 90
#'     outlier_cutoff = 3,    ## Pearson's 3 sigma edit rule
#'     width = 7,            ## window for outlier detection and interpolation
#'     method = "linear"      ## linear interpolation over NAs
#' )
#'
#' \donttest{
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     ## plot original and show where values have been replaced
#'     ## ignore warning about replacing the existing colour scale
#'     plot(data, time_labels = TRUE) +
#'       ggplot2::scale_colour_manual(
#'         name = NULL,
#'         breaks = c("smo2", "replaced"),
#'         values = palette_mnirs(2)
#'       ) +
#'       ggplot2::geom_point(
#'         data = data[data_clean$smo2 != data$smo2, ],
#'         ggplot2::aes(y = smo2, colour = "replaced"),
#'         na.rm = TRUE
#'       ) +
#'       ggplot2::geom_line(
#'         data = {
#'           data_clean[!is.na(data$smo2), "smo2"] <- NA
#'           data_clean
#'         },
#'         ggplot2::aes(y = smo2, colour = "replaced"),
#'         linewidth = 1, na.rm = TRUE
#'       )
#'   }
#' }
#'
#' @rdname replace_mnirs
#' @order 1
#' @export
replace_mnirs <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    invalid_values = NULL,
    invalid_above = NULL,
    invalid_below = NULL,
    outlier_cutoff = NULL,
    width = NULL,
    span = NULL,
    method = c("linear", "median", "locf", "none"),
    verbose = TRUE
) {
    ## validation ====================================
    method <- match.arg(method)
    check_conditions <- c(
        !is.null(c(invalid_values, invalid_above, invalid_below)),
        !is.null(outlier_cutoff),
        method != "none"
    )
    ## do nothing condition
    if (!any(check_conditions)) {
        cli_abort(c(
            "x" = "No replacement criteria specified",
            "i" = "At least one of {.arg invalid_values}, \\
            {.arg invalid_above}, {.arg invalid_below}, \\
            {.arg outlier_cutoff}, or {.arg method} must be specified."
        ))
    }
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }

    validate_mnirs_data(data)
    metadata <- attributes(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, verbose)
    time_channel <- validate_time_channel(enquo(time_channel), data)
    time_vec <- data[[time_channel]]

    if (check_conditions[2L] || method == "median") {
        validate_width_span(width, span, verbose, "for `replace_mnirs()`.")
    }

    ## remove invalid, outliers, and NA ==============================
    data[nirs_channels] <- lapply(data[nirs_channels], \(.x) {
        if (check_conditions[1L]) {
            .x <- replace_invalid(
                x = .x,
                t = time_vec,
                invalid_values = invalid_values,
                invalid_above = invalid_above,
                invalid_below = invalid_below,
                method = "none",
                bypass_checks = TRUE
            )
        }
        if (check_conditions[2L]) {
            .x <- replace_outliers(
                x = .x,
                t = time_vec,
                width = width,
                span = span,
                method = "none",
                outlier_cutoff = outlier_cutoff,
                bypass_checks = TRUE
            )
        }
        if (check_conditions[3L]) {
            .x <- replace_missing(
                x = .x,
                t = time_vec,
                width = width,
                span = span,
                method = method,
                bypass_checks = TRUE
            )
        }
        .x
    })

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)
    metadata$time_channel <- time_channel

    return(create_mnirs_data(data, metadata))
}


#' Replace invalid values
#'
#' `replace_invalid()` detects specified invalid values or range cutoffs in a 
#' numeric vector and replace them with the local median value or `NA`.
#'
#' @param x A numeric vector of the response variable.
#' @param t An *optional* numeric vector of the predictor variable (e.g. time). 
#'   Default is `seq_along(x)`.
#' @inheritParams replace_mnirs
#'
#' @details
#' ## Replace invalid values with with replace_invalid()
#' 
#' Specific `invalid_values` can be replaced, such as `c(0, 100, 102.3)`.
#' Data ranges can be replaced with cutoff values specified by `invalid_above` 
#' and `invalid_below`, where any values higher or lower than the specified 
#' cutoff values (respectively) will be replaced, *inclusive* of the cutoff 
#' values themselves.
#'
#' @returns `replace_invalid()` returns a numeric vector the same length as 
#' `x` with invalid values replaced.
#'
#' @rdname replace_mnirs
#' @order 2
#' @export
replace_invalid <- function(
    x,
    t = seq_along(x),
    invalid_values = NULL,
    invalid_above = NULL,
    invalid_below = NULL,
    width = NULL,
    span = NULL,
    method = c("median", "none"),
    verbose = TRUE,
    ...
) {
    ## validate ===============================================
    args <- list(...)
    if (is.null(c(invalid_values, invalid_above, invalid_below))) {
        cli_abort(c(
            "x" = "No replacement criteria specified",
            "i" = "At least one of {.arg invalid_values}, \\
            {.arg invalid_above}, or {.arg invalid_below} must be specified."
        ))
    }
    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t)
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
    }
    validate_numeric(invalid_values)
    validate_numeric(invalid_above, 1, msg1 = "one-element")
    validate_numeric(invalid_below, 1, msg1 = "one-element")
    method <- match.arg(method)

    ## process ========================================================
    ## fill invalid indices with NA
    invalid_idx <- c(
        which(x %in% invalid_values),
        which(x >= invalid_above),
        which(x <= invalid_below)
    )

    if (length(invalid_idx) == 0) {
        return(x)
    }

    y <- x
    y[invalid_idx] <- NA_real_

    if (method == "median") {
        if (!(args$bypass_checks %||% FALSE)) {
            validate_width_span(width, span, verbose, "for median replacement.")
        }

        window_idx <- compute_local_windows(t, invalid_idx, width, span)
        local_medians <- compute_local_fun(y, window_idx, median_nona)
        ## if method = "median"
        ## invalid_values removed to NA first,
        ## so returns local median excluding idx
        y[invalid_idx] <- local_medians
    }

    return(y)
}


#' Replace local outliers
#'
#' `replace_outliers()` detects local outliers in a numeric vector using a Hampel filter and replaces with the local median value or `NA`.
#'
#' @inheritParams replace_invalid
#'
#' @details
#' ## Outlier detection with replace_outliers()
#'
#' Rolling local medians are computed across `x` within a window defined
#' by `width` (number of samples) or `span` (time span in units of `t`).
#'
#' Outliers are detected with robust median absolute deviation (MAD),
#' adapted from `pracma::hampel()`. Deviations equal to or less than the
#' smallest absolute time series difference in `x` are excluded, to avoid
#' flagging negligible differences where local data have minimal or zero
#' variation.
#'
#' ## Replacement behaviour
#'
#' Values of `x` outside the local bounds defined by `outlier_cutoff` are
#' identified as outliers and either replaced with the local median
#' (`method = "median"`, the *default*) or set to `NA` (`method = "none"`).
#'
#' Existing `NA` values in `x` are *not* replaced. They are passed
#' through to the returned vector. See [replace_missing()].
#'
#' ## Choosing outlier_cutoff
#'
#' `outlier_cutoff` is the number of (MAD-normalised) standard deviations
#' from the local median. Higher values are more conservative; lower
#' values flag more outliers.
#'   - `outlier_cutoff = 3` -- Pearson's 3 sigma edit rule (default).
#'   - `outlier_cutoff = 2` -- approximately Tukey-style 1.5*IQR rule.
#'   - `outlier_cutoff = 0` -- Tukey's median filter (every point
#'     replaced by local median).
#'
#' @returns `replace_outliers()` returns a numeric vector the same length as 
#' `x` with outliers replaced.
#'
#' @rdname replace_mnirs
#' @order 3
#' @export
replace_outliers <- function(
    x,
    t = seq_along(x),
    outlier_cutoff = 3,
    width = NULL,
    span = NULL,
    method = c("median", "none"),
    verbose = TRUE,
    ...
) {
    ## validate ===============================================
    args <- list(...)
    if (!(args$bypass_checks %||% FALSE)) {
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
        validate_x_t(x, t)
        validate_width_span(width, span, verbose, "for `replace_outliers()`.")
    }
    validate_numeric(
        outlier_cutoff, 1, c(0, Inf), msg1 = "one-element positive"
    )
    method <- match.arg(method)

    ## process =====================================================
    outlier_stats <- compute_outliers(x, t, outlier_cutoff, width, span)
    local_medians <- outlier_stats$local_medians
    is_outlier <- outlier_stats$is_outlier

    ## fill outliers with median or NA
    y <- x
    y[is_outlier] <- if (method == "median") {
        local_medians[is_outlier]
    } else {
        NA_real_
    }
    
    return(y)
}

#' Replace missing values
#'
#' `replace_missing()` detects missing (`NA`) values in a numeric vector and
#' replaces via interpolation.
#'
#' @param ... Additional arguments.
#'
#' @inheritParams replace_invalid
#'
#' @details
#' ## Interpolation with replace_missing()
#'
#' `method = "linear"` and `method = "locf"` use [stats::approx()] with
#' `rule = 2`, so leading `NA`s are filled by *"nocb"*
#' (*"next observation carried backward"*) and trailing `NA`s by *"locf"*.
#'
#' `method = "median"` calculates the local median of valid (non-`NA`)
#' values to either side of `NA`s, within a window defined by `width`
#' (number of samples) or `span` (time span in units of `t`). Sequential
#' `NA`s are all replaced by the same median value.
#'
#' ## Edge behaviour for method = "median"
#'
#' If there are no valid values within `span` to one side of the `NA`,
#' the median of the other side is used (i.e. for leading and trailing
#' `NA`s). If there are no valid values within either side, the first
#' valid sample on either side is used (equivalent to
#' `replace_missing(x, width = 1)`).
#'
#' @returns `replace_missing()` returns a numeric vector the same length as
#' `x` with missing values replaced.
#'
#' @rdname replace_mnirs
#' @order 4
#' @export
replace_missing <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    method = c("linear", "median", "locf"),
    verbose = TRUE,
    ...
) {
    ## validate ===============================================
    args <- list(...)
    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t)
    }
    method <- match.arg(method)
    if (method == "locf") {
        method <- "constant" ## swap for approx method arg
    }

    ## process ==============================================
    if (method %in% c("linear", "constant")) {
        y <- stats::approx(
            x = t,
            y = x,
            xout = args$xout %||% t, ## = t unless explicitly specified by hidden option
            method = method, ## c("linear", "constant")
            rule = 2, ## fill leading and trailing `NA`s
            f = 0, ## locf if method = "constant"
            ties = list("ordered", mean) ## assume ordered, take mean of ties
        )$y
    } else if (method == "median") {
        if (!(args$bypass_checks %||% FALSE)) {
            if (missing(verbose)) {
                verbose <- getOption("mnirs.verbose", default = TRUE)
            }
            validate_width_span(width, span, verbose, "for median replacement.")
        }
        ## median of width or span VALID values to either side of sequential NAs
        y <- x
        na_idx <- which(is.na(x))
        window_idx <- compute_valid_neighbours(x, t, width, span, verbose)
        local_medians <- compute_local_fun(x, window_idx, median_nona)
        y[na_idx] <- local_medians
    }

    return(y)
}
