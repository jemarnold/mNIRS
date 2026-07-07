#' Filter a data frame
#'
#' @description
#' Apply digital filtering/smoothing to numeric vector data within a data frame
#' using either:
#'   1. A cubic smoothing spline.
#'   2. A Butterworth digital filter.
#'   3. A simple moving average.
#'
#' Note the `method`-specific arguments below.
#'
#' @param method A character string indicating how to filter the data.
#'   Additional arguments must be specified for each method. See *Details*.
#'   \describe{
#'      \item{`"smooth_spline"`}{Fits a cubic smoothing spline. Additional
#'      arguments: `spar`.}
#'      \item{`"butterworth"`}{Uses a centred Butterworth digital filter.
#'      Additional arguments: `order`, `W` or `fc`, `sample_rate`, `type`,
#'      `edges`. See [filter_butterworth()].}
#'      \item{`"moving_average"`}{Uses a centred moving average filter.
#'      Additional arguments: `width` or `span`, `partial`. See
#'      [filter_moving_average()].}
#'   }
#' @param na.rm Logical; default is `FALSE`, propagates any `NA`s to the
#'   returned vector. If `TRUE`, ignores `NA`s and processes available valid
#'   samples within the local window. May return errors or warnings. (see
#'   *Details*).
#' @param ... Additional arguments passed to the underlying method function.
#'   See *Details*.
#' @param spar **smooth_spline**: A numeric smoothing parameter passed to
#'   [stats::smooth.spline()]. If `NULL` (*default*), automatically
#'   determined via penalised log likelihood.
#' @param order **butterworth**: An integer defining the filter order
#'   (*default* `order = 2`).
#' @param W **butterworth**: A one- or two-element numeric vector within
#'   `[0, 1]` defining the filter cutoff frequency(ies) as a fraction of
#'   the Nyquist frequency (see *Details*). One of either `W` or `fc` must be
#'   specified.
#' @param fc **butterworth**: A one- or two-element numeric vector defining
#'   the filter absolute cutoff frequency in Hz. Used with `sample_rate` to
#'   compute `W`. One of either `W` or `fc` must be specified.
#' @param sample_rate **butterworth**: A numeric sample rate in Hz. Will
#'   be taken from metadata or estimated from `time_channel` if not
#'   defined.
#' @param type **butterworth**: A character string specifying filter type,
#'   one of: `c("low", "high", "stop", "pass")` (`"low"` is the
#'   *default*).
#' @param edges **butterworth**: A character string specifying the edge
#'   padding, one of: `c("rev", "rep1", "none")` (`"rev"` is the
#'   *default*). See [filter_butterworth()].
#' @param width **moving_average**: An integer number of samples within
#'   the local window. One of either `width` or `span` must be specified.
#' @param span **moving_average**: A numeric time duration in units of
#'   `time_channel` within the local window. One of either `width` or
#'   `span` must be specified.
#' @param partial **moving_average**: Logical; default is `FALSE`, only
#'   returns values where a full window of valid (non-`NA`) samples are
#'   available. If `TRUE`, ignores `NA` and processes available valid samples
#'   (see *Details*).
#' @inheritParams map_mnirs_intervals
#' @inheritParams validate_mnirs
#'
#' @inheritSection map_mnirs_intervals Data input formats
#'
#' @details
#' ## method = "smooth_spline"
#'
#' Aliases: `method = c("smooth spline", "spline")`
#'
#' Applies a non-parametric cubic smoothing spline from
#' [stats::smooth.spline()]. Smoothing is defined by the parameter `spar`,
#' which can be left as `NULL` and automatically determined via penalised
#' log likelihood. This usually works well for responses occurring on the
#' order of minutes or longer. `spar` can be specified typically, but not
#' necessarily, in the range `spar = [0, 1]`.
#'
#' ## method = "butterworth"
#'
#' Aliases: `method = c("butter")`
#'
#' Applies a centred (two-pass symmetrical) Butterworth digital filter
#' from [signal::butter()] and [signal::filtfilt()].
#'
#' Filter `type` defines how the desired signal frequencies are either
#' passed or rejected from the output signal. *Low-pass* and *high-pass*
#' filters allow only frequencies *lower* or *higher* than the cutoff
#' frequency, respectively to be passed through to the output signal.
#' *Stop-band* defines a critical range of frequencies which are rejected
#' from the output signal. *Pass-band* defines a critical range of
#' frequencies which are passed through as the output signal.
#'
#' The filter order (number of passes) is defined by `order`, typically
#' in the range `order = [1, 10]`. Higher filter order tends to capture
#' more rapid changes in amplitude, but also causes more distortion
#' around those change points in the signal. General advice is to use
#' the lowest filter order which sufficiently captures the desired rapid
#' responses in the data.
#'
#' The critical (cutoff) frequency can be defined by `W`, a numeric value
#' for *low-pass* and *high-pass* filters, or a two-element vector
#' `c(low, high)` defining the lower and upper bands for *stop-band*
#' and *pass-band* filters. `W` represents the desired fractional cutoff
#' frequency in the range `W = [0, 1]`, where `1` is the Nyquist
#' frequency, i.e., half the `sample_rate` of the data in Hz.
#'
#' Alternatively, the cutoff frequency can be defined by `fc` and
#' `sample_rate` together. `fc` represents the desired cutoff frequency
#' directly in Hz, and `sample_rate` is the sample rate of the recorded data
#' in Hz. Where `W = fc / (sample_rate / 2)`.
#'
#' Only one of either `W` or `fc` should be defined. If both are
#' defined, `W` will be preferred over `fc`.
#'
#' ## method = "moving_average"
#'
#' Aliases: `method = c("moving average", "ma")`
#'
#' Applies a centred (symmetrical) moving average filter in a local
#' window, defined by either `width` as the number of samples around
#' `idx` between `[idx - floor(width/2), idx + floor(width/2)]`. Or by
#' `span` as the timespan in units of `time_channel` between
#' `[t - span/2, t + span/2]`.
#'
#' ## Missing values
#'
#' Missing values (`NA`) in `nirs_channels` will cause an error for
#' `method = "smooth_spline"` or `"butterworth"`, unless `na.rm = TRUE`.
#' Then `NA`s will be ignored and passed through to the returned data.
#'
#' For `method = "moving_average"`, `na.rm` controls whether `NA`s within
#' each local window are either propagated to the returned vector when
#' `na.rm = FALSE` (the default), or ignored before processing if
#' `na.rm = TRUE`.
#'
#' @inheritSection replace_mnirs Per-channel arguments
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class *"mnirs"* with metadata
#'   available with `attributes()`. For list or grouped data frame input,
#'   returns a named list of *"mnirs"* tibbles, one per interval.
#'
#' @examples
#' ## read example data and clean for outliers
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2 = "SmO2 Live"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' ) |>
#'     replace_mnirs(
#'         invalid_values = c(0, 100),
#'         outlier_cutoff = 3,
#'         width = 7,
#'         verbose = FALSE
#'     )
#'
#' data
#'
#' data_filtered <- filter_mnirs(
#'     data,                   ## blank channels will be retrieved from metadata
#'     method = "butterworth", ## Butterworth digital filter is a common choice
#'     order = 2,              ## filter order number
#'     W = 0.02,               ## filter fractional critical frequency `[0, 1]`
#'     type = "low",           ## specify a "low-pass" filter
#'     na.rm = TRUE            ## explicitly ignore NAs
#' )
#'
#' ## note the smoothed `smo2` values
#' data_filtered
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ## plot filtered data and add the raw data back to the plot to compare
#'         plot(data_filtered, time_labels = TRUE) +
#'             ggplot2::geom_line(
#'                 data = data,
#'                 ggplot2::aes(y = smo2, colour = "smo2"), alpha = 0.4
#'             )
#'     }
#' }
#'
#' @rdname filter_mnirs
#' @export
filter_mnirs <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method = c("smooth_spline", "butterworth", "moving_average"),
    na.rm = FALSE,
    verbose = TRUE,
    ...,
    spar = NULL,
    order = 2L,
    W = NULL,
    fc = NULL,
    sample_rate = NULL,
    type = c("low", "high", "stop", "pass"),
    edges = c("rev", "rep1", "none"),
    width = NULL,
    span = NULL,
    partial = FALSE
) {
    ## list or grouped input → normalise to named list, recurse per interval
    if (inherits(data, "grouped_df") || !is.data.frame(data)) {
        return(map_mnirs_intervals(data, match.call(), parent.frame()))
    }

    ## validation ====================================
    validate_mnirs_data(data)
    ## normalise method aliases before matching
    method <- gsub(
        "^(ma|moving[ _-]average)$",
        "moving_average",
        method,
        ignore.case = TRUE
    )
    method <- gsub(
        "^(spline|smooth[ _-]spline)$",
        "smooth_spline",
        method,
        ignore.case = TRUE
    )
    method <- match.arg(method)

    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }

    UseMethod(
        "filter_mnirs",
        structure(data, class = c(method, "mnirs_filtered"))
    )
}


#' @rdname filter_mnirs
#' @usage NULL
#' @export
filter_mnirs.smooth_spline <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    ## validation ==========================================
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    metadata <- attributes(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, env)
    time_channel <- validate_time_channel(enquo(time_channel), data, env = env)

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(spar = list(...)$spar),
        verbose = verbose,
        env = env
    )

    ## processing ==========================================
    t_vec <- data[[time_channel]]

    if (anyDuplicated(t_vec)) {
        cli_abort(c(
            "x" = "{.arg time_channel} has duplicated or irregular samples.",
            "i" = "Re-sample first with {.fn mnirs::resample_mnirs}."
        ), call = env)
    }

    data[nirs_channels] <- lapply(nirs_channels, \(.x) {
        spar <- per_channel[[.x]]$spar
        validate_numeric(
            spar, 1, c(0, Inf), FALSE, msg1 = "one-element positive", env = env
        )
        x <- data[[.x]]
        ## handle NAs
        handle_na <- na.rm && anyNA(x)
        if (handle_na) {
            na_info <- preserve_na(x)
            x <- na_info$x_valid
            t_vec <- t_vec[!na_info$na_idx]
        } else if (anyNA(x)) {
            cli_abort(c(
                "x" = "{.arg nirs_channels} = {.field {.x}} contains internal \\
                {.val {NA}}'s.",
                "i" = "Set {.arg na.rm = TRUE} to ignore {.val {NA}}'s."
            ), call = env)
        }

        spline_model <- stats::smooth.spline(x = t_vec, y = x, spar = spar)

        if (is.null(spar) && verbose) {
            cli_inform(c(
                "i" = "{.arg nirs_channels} = {.field {.x}}: \\
                `smooth.spline(spar = {.val {round(spline_model$spar, 3)}})`"
            ), call = env)
        }

        if (handle_na) {
            restore_na(spline_model$y, na_info)
        } else {
            spline_model$y
        }
    })

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)

    return(create_mnirs_data(data, metadata))
}


#' @rdname filter_mnirs
#' @usage NULL
#' @export
filter_mnirs.butterworth <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    ## validation ==========================================
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    metadata <- attributes(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, env)
    time_channel <- validate_time_channel(enquo(time_channel), data, env = env)
    args <- list(...)
    sample_rate <- args$sample_rate
    sample_rate <- validate_sample_rate(
        data, time_channel, sample_rate, verbose, env = env
    )

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(
            order = args$order %||% 2L,
            W = args$W,
            fc = args$fc,
            type = args$type %||% "low",
            edges = args$edges %||% "rev"
        ),
        defaults = list(order = 2L, type = "low", edges = "rev"),
        choices = list(
            type = c("low", "high", "stop", "pass"),
            edges = c("rev", "rep1", "none")
        ),
        verbose = verbose,
        env = env
    )

    ## processing ==========================================
    data[nirs_channels] <- Map(\(.nirs, .a) {
        ## verbose validator hints emitted once for the first channel
        .v <- verbose && .nirs == nirs_channels[[1L]]

        if (is.null(c(.a$W, .a$fc))) {
            cli_abort(c(
                "x" = "Cutoff frequency undefined.",
                "i" = "One of {.arg W} or {.arg fc} must be defined for a \\
                Butterworth filter."
            ), call = env)
        }

        fc_n <- if (.a$type %in% c("low", "high")) 1 else 2
        ## order & W are validated in filter_butterworth
        validate_numeric(
            .a$fc, fc_n, c(0, Inf), inclusive = FALSE,
            msg1 = paste0(fc_n, "-element positive"), env = env
        )

        if (!is.null(.a$W) && !is.null(.a$fc)) {
            .a$fc <- NULL
            if (.v) {
                cli_inform(c(
                    "i" = "{.val Butterworth} parameter {.arg W} = \\
                    {.val {(.a$W)}} overrides {.arg fc}."
                ), call = env)
            }
        }

        if (is.null(.a$W) && !is.null(.a$fc) && !is.null(sample_rate)) {
            nq <- sample_rate * 0.5 ## nyquist frequency
            .a$W <- .a$fc / nq
            if (.a$W > 1 | .a$W <= 0) {
                cli_abort(c(
                    "x" = "{.arg fc} must be between {.val {0}} and half \\
                    the {.arg sample_rate} ({.val {signif(nq, 3)}} Hz)"
                ), call = env)
            }
        }

        filter_butterworth(
            data[[.nirs]], .a$order, .a$W, .a$type, .a$edges, na.rm, env = env
        )
    }, nirs_channels, per_channel[nirs_channels])

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)
    metadata$time_channel <- time_channel
    metadata$sample_rate <- sample_rate

    return(create_mnirs_data(data, metadata))
}


#' @rdname filter_mnirs
#' @usage NULL
#' @export
filter_mnirs.moving_average <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    method,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    ## validation ==========================================
    ## report conditions as coming from the user-facing generic
    env <- sys.call(-1)
    metadata <- attributes(data)
    nirs_channels <- validate_nirs_channels(enquo(nirs_channels), data, env)
    time_channel <- validate_time_channel(enquo(time_channel), data, env = env)
    args <- list(...)

    ## broadcast global args, applying any per-channel list() overrides
    per_channel <- resolve_channel_args(
        nirs_channels,
        args = list(
            width = args$width,
            span = args$span,
            partial = args$partial %||% FALSE
        ),
        defaults = list(partial = FALSE),
        verbose = verbose,
        env = env
    )

    ## processing ==========================================
    t_vec <- data[[time_channel]]

    data[nirs_channels] <- Map(\(.nirs, .a) {
        filter_moving_average(
            x = data[[.nirs]],
            t = t_vec,
            width = .a$width,
            span = .a$span,
            partial = .a$partial,
            na.rm = na.rm,
            verbose = verbose,
            bypass_checks = TRUE,
            env = env
        )
    }, nirs_channels, per_channel[nirs_channels])

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_channels)
    metadata$time_channel <- time_channel

    return(create_mnirs_data(data, metadata))
}


#' Apply a moving average filter
#'
#' Apply a simple moving average smoothing filter to vector data.
#' `filter_ma()` is an alias of `filter_moving_average()`.
#'
#' 
#' @param partial Logical; default is `FALSE`, only returns values where a full
#'   window of valid (non-`NA`) samples are available. If `TRUE`, ignores `NA`
#'   and processes available valid samples (see *Details*).
#' @inheritParams replace_invalid
#' @inheritParams filter_mnirs
#'
#' @details
#' ## Rolling window
#'
#' Applies a centred (symmetrical) moving average filter in a local
#' window, defined by either `width` as the number of samples around
#' `idx` between `[idx - floor(width/2), idx + floor(width/2)]`. Or by
#' `span` as the timespan in units of `time_channel` between
#' `[t - span/2, t + span/2]`.
#'
#' ## Partial windows
#'
#' The default `partial = FALSE` requires a complete number of samples
#' specified by `width` or `span` (estimated from the sample rate of `t` when
#' `span` is used). `NA` is returned if fewer samples are present in the
#' local window.
#'
#' Setting `partial = TRUE` allows computation with only a single valid sample,
#' such as at edge conditions. But these values will be more sensitive to
#' noise and should be used with caution.
#'
#' ## Missing values
#'
#' `na.rm` controls whether missing values (`NA`s) within each local window are
#' either propagated to the returned vector when `na.rm = FALSE` (the default),
#' or ignored before processing if `na.rm = TRUE`.
#'
#' @returns A numeric vector the same length as `x`.
#'
#' @examples
#' x <- c(1, 3, 2, 5, 4, 6, 5, 7)
#' t <- c(0, 1, 2, 4, 5, 6, 7, 10)  ## irregular time with gaps
#'
#' ## width: centred window of 3 samples
#' filter_moving_average(x, width = 3)
#'
#' ## partial = TRUE fills edge values with a narrower window
#' filter_moving_average(x, width = 3, partial = TRUE)
#'
#' ## span: centred window of 2 time-units (accounts for irregular sampling)
#' filter_moving_average(x, t, span = 2)
#'
#' ## na.rm = FALSE (default): any NA in the window propagates to the result
#' x_na <- c(1, NA, 3, 4, 5, NA, 7, 8)
#' filter_moving_average(x_na, width = 3, na.rm = FALSE)
#'
#' ## na.rm = TRUE: skip NAs and return the local mean of local valid values
#' filter_moving_average(x_na, width = 3, partial = TRUE, na.rm = TRUE)
#'
#' @rdname filter_moving_average
#' @export
filter_moving_average <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    partial = FALSE,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    ## validation ===========================================
    ## internal callers pass `env` through `...` to report conditions
    ## as coming from the user-facing function
    env <- list(...)$env %||% environment()
    if (!(list(...)$bypass_checks %||% FALSE)) {
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
    }
    validate_x_t(x, t, env = env)
    validate_width_span(
        width, span, verbose, "for a moving average filter.", env = env
    )

    ## handle NAs
    if (verbose && !na.rm && anyNA(x)) {
        cli_warn(c(
            "!" = "{.arg x} contains internal {.val {NA}}'s.",
            "i" = "Set {.arg na.rm = TRUE} to ignore {.val {NA}}'s."
        ), call = warn_call(env))
    }

    ## processing ==============================================
    window_idx <- compute_local_windows(
        t,
        width = width,
        span = span,
        env = env
    )

    if (!partial) {
        ## min_obs default to estimated width when span is specified
        ## less strict span_width - 2 to allow start & end buffer
        ## with irregular t values
        min_obs <- max(
            width %||% (floor(span * estimate_sample_rate(t, env)) - 2L),
            1L
        )

        ## error if fewer valid samples than min_obs
        if (sum(is.finite(x)) < min_obs) {
            cli_abort(c(
                "x" = "Insufficient valid samples detected.",
                "i" = "{.arg width} or {.arg span} must be smaller than \\
                the range of {.arg x} when {.arg partial} = {.val {FALSE}}."
            ), call = env)
        }

        which_partial <- lengths(window_idx) < min_obs
    }

    y <- vapply(window_idx, \(.idx) mean(x[.idx], na.rm = na.rm), numeric(1))

    if (!partial) {
        ## exclude incomplete windows (at edges)
        y[which_partial] <- NA_real_
    }
    ## NaN to NA
    y[!is.finite(y)] <- NA_real_
    return(y)
}


#' @rdname filter_moving_average
#' @export
filter_ma <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    partial = FALSE,
    na.rm = FALSE,
    verbose = TRUE,
    ...
) {
    ## report conditions as coming from this wrapper unless overridden
    env <- list(...)$env %||% environment()
    filter_moving_average(
        x = x,
        t = t,
        width = width,
        span = span,
        partial = partial,
        na.rm = na.rm,
        verbose = verbose,
        env = env,
        ...
    )
}


#' Apply a Butterworth digital filter
#'
#' Apply a Butterworth digital filter to vector data with [signal::butter()]
#' and [signal::filtfilt()] which handles 'edges' better at the start and end
#' of the data.
#'
#' @param x A numeric vector.
#' @param order An integer defining the filter order (*default* `order = 2`).
#' @param W A one- or two-element numeric vector within `[0, 1]` defining the
#'   filter cutoff frequency(ies) as a fraction of the Nyquist frequency
#'   (see *Details*).
#' @param type A character string indicating the digital filter type (see
#'   *Details*).
#'   \describe{
#'      \item{`"low"`}{For a *low-pass* filter (the *default*).}
#'      \item{`"high"`}{For a *high-pass* filter.}
#'      \item{`"stop"`}{For a *stop-band* (band-reject) filter.}
#'      \item{`"pass"`}{For a *pass-band* filter.}
#'   }
#' @param edges A character string indicating edge detection padding for `x`.
#'   \describe{
#'      \item{`"rev"`}{Will pad `x` with the preceding 5% data in reverse
#'      sequence (*the default*).}
#'      \item{`"rep1"`}{Will pad `x` by repeating the last preceding value.}
#'      \item{`"none"`}{Will return the unpadded [signal::filtfilt()] output.}
#'   }
#' @inheritParams filter_mnirs
#'
#' @details
#' Applies a centred (two-pass symmetrical) Butterworth digital filter from
#'   [signal::butter()] and [signal::filtfilt()].
#'
#' Filter `type` defines how the desired signal frequencies are either
#'   passed or rejected from the output signal. *Low-pass* and *high-pass*
#'   filters allow only frequencies *lower* or *higher* than the cutoff
#'   frequency `W` to be passed through as the output signal, respectively.
#'   *Stop-band* defines a critical range of frequencies which are rejected
#'   from the output signal. *Pass-band* defines a critical range of
#'   frequencies which are passed through as the output signal.
#'
#' The filter order (number of passes) is defined by `order`, typically in
#'   the range `order = [1, 10]`. Higher filter order tends to capture more
#'   rapid changes in amplitude, but also causes more distortion around
#'   those change points in the signal. General advice is to use the
#'   lowest filter order which sufficiently captures the desired rapid
#'   responses in the data.
#'
#' The critical (cutoff) frequency is defined by `W`, a numeric value for
#'   *low-pass* and *high-pass* filters, or a two-element vector
#'   `c(low, high)` defining the lower and upper bands for *stop-band* and
#'   *pass-band* filters. `W` represents the desired fractional cutoff
#'   frequency in the range `W = [0, 1]`, where `1` is the Nyquist
#'   frequency, i.e., half the sample rate of the data in Hz.
#'
#' Missing values (`NA`) in `x` will cause an error unless `na.rm = TRUE`.
#'   Then `NA`s will be ignored and passed through to the returned vector.
#'
#' @returns A numeric vector the same length as `x`.
#'
#' @seealso [signal::filtfilt()], [signal::butter()]
#'
#' @examples
#' set.seed(13)
#' sin <- sin(2 * pi * 1:150 / 50) * 20 + 40
#' noise <- rnorm(150, mean = 0, sd = 6)
#' noisy_sin <- sin + noise
#' without_edge_detection <- filter_butterworth(
#'     x = noisy_sin,
#'     order = 2,
#'     W = 0.1,
#'     edges = "none"
#' )
#' with_edge_detection <- filter_butterworth(
#'     x = noisy_sin,
#'     order = 2,
#'     W = 0.1,
#'     edges = "rep1"
#' )
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' ggplot2::ggplot(data.frame(), ggplot2::aes(x = seq_along(noise))) +
#'     theme_mnirs() +
#'     scale_colour_mnirs(name = NULL) +
#'     ggplot2::geom_line(ggplot2::aes(y = noisy_sin)) +
#'     ggplot2::geom_line(
#'         ggplot2::aes(y = without_edge_detection, colour = "without_edge_detection")
#'     ) +
#'     ggplot2::geom_line(
#'         ggplot2::aes(y = with_edge_detection, colour = "with_edge_detection")
#'     )
#'
#' @export
filter_butterworth <- function(
    x,
    order = 2L,
    W,
    type = c("low", "high", "stop", "pass"),
    edges = c("rev", "rep1", "none"),
    na.rm = FALSE,
    ...
) {
    ## validation ============================================
    ## internal callers pass `env` through `...` to report conditions
    ## as coming from the user-facing function
    env <- list(...)$env %||% environment()
    check_installed("signal", "to use Butterworth digital filter")
    validate_numeric(x, env = env)
    validate_numeric(
        order, 1, c(1, Inf), integer = TRUE, msg1 = "one-element positive",
        env = env
    )
    type <- match.arg(type)
    W_n <- if (type %in% c("low", "high")) 1 else 2
    validate_numeric(
        W, W_n, c(0, 1), inclusive = FALSE,
        msg1 = paste0(W_n, "-element positive"),
        msg2 = "between {col_blue('[0, 1]')}.",
        env = env
    )
    edges <- match.arg(edges)

    ## processing ==============================================
    ## handle NAs
    handle_na <- na.rm && anyNA(x)
    if (handle_na) {
        na_info <- preserve_na(x)
        x <- na_info$x_valid
    } else if (anyNA(x)) {
        cli_abort(c(
            "x" = "{.arg x} contains internal {.val {NA}}'s.",
            "i" = "Set {.arg na.rm = TRUE} to ignore {.val {NA}}'s."
        ), call = env)
    }

    if (edges == "none") {
        y <- signal::filtfilt(signal::butter(n = order, W, type), x = x)
    } else {
        x_n <- length(x)
        pad <- max(1L, x_n %/% 20) ## 5% padded length

        padded <- switch(
            edges,
            "rev" = c(x[pad:1], x, x[x_n:(x_n - pad + 1)]),
            "rep1" = c(rep(x[1], pad), x, rep(x[x_n], pad))
        )

        y_padded <- signal::filtfilt(
            signal::butter(n = order, W, type),
            x = padded
        )
        y <- y_padded[(pad + 1):(pad + x_n)]
    }

    ## return NAs to original positions in y
    if (handle_na) {
        return(restore_na(y, na_info))
    } else {
        return(y)
    }
}


#' Apply a Butterworth digital filter
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `filter_butter()` was renamed to [filter_butterworth()]
#'
#' @inheritParams filter_butterworth
#'
#' @returns A numeric vector the same length as `x`.
#'
#' @keywords internal
#' @export
filter_butter <- function(
    x,
    order = 2L,
    W,
    type = c("low", "high", "stop", "pass"),
    edges = c("rev", "rep1", "none"),
    na.rm = FALSE,
    ...
) {
    lifecycle::deprecate_warn(
        when = "0.7.0",
        what = "filter_butter()",
        with = "filter_butterworth()"
    )
    filter_butterworth(
        x = x,
        order = order,
        W = W,
        type = type,
        edges = edges,
        na.rm = na.rm,
        ...
    )
}