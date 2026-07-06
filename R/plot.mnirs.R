#' Plot *{mnirs}* objects
#'
#' Create a base plot for data frames or lists of data frames with class
#' *"mnirs"*.
#'
#' @param x Data frame or list of data frames of class *"mnirs"* (e.g. from
#'   [extract_intervals()]). List input produces a faceted plot with one
#'   panel per element.
#' @param points Logical. Default is `FALSE`. If `TRUE` displays
#'   `ggplot2::geom_points()`. Otherwise displays `ggplot2::geom_lines()`.
#' @param time_labels Logical. Default is `FALSE`. If `TRUE` displays x-axis
#'   time values formatted as *"h:mm:ss"* using [format_hmmss()]. Otherwise,
#'   x-axis values are displayed as numeric.
#' @param na.omit Logical. Default is `FALSE`. If `TRUE` omits missing (`NA`)
#'   and non-finite `c(Inf, -Inf, NaN)` from display.
#' @param ... Additional arguments.
#'
#' @details
#' When `x` is a named list of *"mnirs"* data frames, elements are bound into a
#' single data frame and displayed as faceted panels via
#' [ggplot2::facet_wrap()].
#'
#' Accepts some arguments in `...`, such as `nrow`, `ncol`, and `scales`
#' passed to [ggplot2::facet_wrap()]. `n.breaks` overrides the default number
#' of y-axis breaks. `breaks` overrides the x-axis breaks directly.
#'
#' @returns A [ggplot2][ggplot2::ggplot()] object.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' data <- read_mnirs(
#'     example_mnirs("train.red"),
#'     nirs_channels = c(smo2 = "SmO2"),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     verbose = FALSE
#' )
#'
#' ## plot time labels as "h:mm:ss"
#' plot(data, time_labels = TRUE)
#'
#' data_list <- extract_intervals(
#'     data,
#'     start = by_time(2452, 3168),
#'     span = c(-60, 120),
#'     verbose = FALSE
#' )
#'
#' ## plot a list of mnirs data frames as faceted panels
#' plot(data_list, time_labels = TRUE)
#'
#' @export
plot.mnirs <- function(
    x,
    points = FALSE,
    time_labels = FALSE,
    na.omit = FALSE,
    ...
) {
    check_installed(c("ggplot2"), reason = "to plot mNIRS data")

    args <- list(...)

    ## handle list of mnirs data frames
    if (is.list(x) && !is.data.frame(x)) {
        x <- as_plot_data(x)
    }

    nirs_channels <- attr(x, "nirs_channels")
    time_channel <- attr(x, "time_channel")

    ## pre-compute conditionals
    x_name <- if (time_labels) {
        ## match format_hmmss() output: "mm:ss" below one hour
        units <- if (max(abs(x[[time_channel]]), na.rm = TRUE) < 3600) {
            "(mm:ss)"
        } else {
            "(h:mm:ss)"
        }
        paste(time_channel, units)
    } else {
        ggplot2::waiver()
    }
    x_breaks <- if (!is.null(args[["breaks"]])) {
        args[["breaks"]]
    } else if (time_labels) {
        breaks_timespan()
    } else if (rlang::is_installed("scales")) {
        scales::breaks_pretty()
    } else {
        ggplot2::waiver()
    }
    x_labels <- if (time_labels) {
        format_hmmss
    } else {
        ggplot2::waiver()
    }
    y_breaks <- if (rlang::is_installed("scales")) {
        scales::breaks_pretty(n = args[["n.breaks"]] %||% 5)
    } else {
        ggplot2::waiver()
    }

    ## build base plot with axis configuration
    plot <- ggplot2::ggplot(x) +
        ggplot2::aes(x = .data[[time_channel]]) +
        theme_mnirs() +
        ggplot2::labs(x = x_name, y = "mNIRS") +
        ggplot2::scale_x_continuous(
            breaks = x_breaks,
            labels = x_labels,
            expand = ggplot2::expansion(mult = 0.01)
        ) +
        ggplot2::scale_y_continuous(
            breaks = y_breaks,
            expand = ggplot2::expansion(mult = 0.01)
        ) +
        scale_colour_mnirs(
            name = NULL,
            guide = ggplot2::guide_legend(
                override.aes = list(linewidth = 1)
            )
        )

    ## add one geom per channel
    layers <- lapply(nirs_channels, function(ch) {
        ch_data <- if (na.omit) x[is.finite(x[[ch]]), ] else x
        ch_aes <- ggplot2::aes(y = .data[[ch]], colour = ch)
        c(
            list(ggplot2::geom_line(ch_aes, data = ch_data)),
            if (points) {
                list(ggplot2::geom_point(ch_aes, data = ch_data, size = 3))
            }
        )
    })

    ## facet when plotting multiple mnirs data frames
    if ("interval" %in% names(x)) {
        facet_args <- intersect(
            names(args),
            names(formals(ggplot2::facet_wrap))
        )
        scales_arg <- args[["scales"]] %||% "free_x"
        plot <- plot +
            do.call(
            ggplot2::facet_wrap,
            c(
                list(facets = ~interval, scales = scales_arg),
                args[setdiff(facet_args, "scales")]
            )
        )
    }

    return(plot + layers)
}


#' Plot *{mnirs}* kinetics results
#'
#' Create a default plot for an *"mnirs_kinetics"* object returned from
#' [analyse_kinetics()]. Observed signals are drawn per `nirs_channel`, faceted
#' by interval, with the fitted response overlaid and the key kinetics
#' coefficient(s) annotated per panel.
#'
#' @param x An *"mnirs_kinetics"* object from [analyse_kinetics()].
#' @param fitted Logical. Default is `TRUE`; overlays a dashed fitted curve for
#'   parametric methods (`"peak_slope"`, `"monoexponential"`, `"sigmoidal"`).
#'   `"response_time"` has no fitted curve.
#' @param markers Logical. Default is `TRUE`; draws a dotted vertical line at
#'   the response onset (`start_time`) and key coefficient points.
#' @param labels Logical. Default is `TRUE`; annotates each panel with the key
#'   coefficient value(s) for the fitted method.
#' @param ... Additional arguments passed to [plot.mnirs()], such as `points`,
#'   `time_labels`, `nrow`, `ncol`, or `scales`.
#'
#' @returns A [ggplot2][ggplot2::ggplot()] object.
#'
#' @seealso [analyse_kinetics()], [plot.mnirs()]
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' result <- read_mnirs(
#'     example_mnirs("train.red"),
#'     nirs_channels = c(smo2 = "SmO2"),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     zero_time = TRUE,
#'     verbose = FALSE
#' ) |>
#'     resample_mnirs(method = "linear", verbose = FALSE) |>
#'     extract_intervals(
#'         group_intervals = "distinct",
#'         start = by_time(368, 1084),
#'         span = c(-20, 90),
#'         zero_time = TRUE,
#'         verbose = FALSE
#'     ) |>
#'     analyse_kinetics(
#'         method = "peak_slope",
#'         span = 10,
#'         verbose = FALSE
#'     )
#'
#' plot(result)
#'
#' @export
plot.mnirs_kinetics <- function(
    x,
    fitted = TRUE,
    markers = TRUE,
    labels = TRUE,
    ...
) {
    check_installed("ggplot2", reason = "to plot mNIRS data")

    ## open white marker for every kinetics key point
    key_point <- function(mapping, data, ...) {
        ggplot2::geom_point(
            mapping,
            data = data,
            size = 3,
            shape = 21,
            stroke = 1,
            fill = "white",
            show.legend = FALSE,
            ...
        )
    }

    ## observed signal + facet + theme via existing plot.mnirs
    p <- plot(x$data, ...)

    ## bound frame: interval factor (when >1) + <channel>_fitted columns
    plot_data <- as_plot_data(x$data)
    nirs <- attr(plot_data, "nirs_channels")
    time_channel <- attr(plot_data, "time_channel")
    faceted <- "interval" %in% names(plot_data)

    ## attach the resolved onset to each row for method-aware fitted overlay
    onset <- x$interval_times[c("interval", "start_times")]
    plot_data <- if (faceted) {
        merge(plot_data, onset, by = "interval", sort = FALSE)
    } else {
        transform(plot_data, start_times = onset$start_times[[1L]])
    }

    ## fitted overlay ==========================================
    ## parametric methods only: continuous dashed fitted curve in the
    ## channel colour. response_time has no curve; its points are markers.
    if (fitted && x$method != "response_time") {
        p <- p +
            lapply(nirs, \(.ch) {
            fcol <- paste0(.ch, "_fitted")
            ggplot2::geom_line(
                ggplot2::aes(y = .data[[fcol]], colour = .ch),
                data = plot_data[is.finite(plot_data[[fcol]]), , drop = FALSE],
                linetype = "dashed",
                linewidth = 1,
                show.legend = FALSE
            )
        })
    }

    ## per-interval-per-channel markers and labels ============
    ann <- kinetics_annotations(x)
    if (!faceted) {
        ann$interval <- NULL
    }

    if (markers) {
        ## dotted onset line at the resolved start_time per interval
        vdata <- unique(
            ann[intersect(c("interval", "start_times"), names(ann))]
        )
        p <- p +
            ggplot2::geom_vline(
                ggplot2::aes(xintercept = start_times),
                data = vdata,
                linetype = "dotted",
                colour = "grey50"
            )
        if (x$method == "response_time") {
            ## response and extreme (fitted values after the onset) as points
            p <- p +
                lapply(nirs, \(.ch) {
                fcol <- paste0(.ch, "_fitted")
                post <- is.finite(plot_data[[fcol]]) &
                    plot_data[[time_channel]] > plot_data$start_times
                key_point(
                    ggplot2::aes(y = .data[[fcol]], colour = .ch),
                    plot_data[post, , drop = FALSE]
                )
            })
            ## baseline as a single point at the onset (start_time, A)
            base_pts <- merge(
                x$coefficients[c("interval", "nirs_channels", "A")],
                x$interval_times[c("interval", "start_times")],
                by = "interval",
                sort = FALSE
            )
            if (!faceted) {
                base_pts$interval <- NULL
            }
            p <- p +
                key_point(
                    ggplot2::aes(
                        x = start_times,
                        y = A,
                        colour = nirs_channels
                    ),
                    base_pts,
                    inherit.aes = FALSE
                )
        } else {
            ## single key-point marker for parametric methods
            p <- p +
                key_point(
                    ggplot2::aes(x = xval, y = yval, colour = nirs_channels),
                    ann[is.finite(ann$xval), , drop = FALSE],
                    inherit.aes = FALSE
                )
        }
    }

    if (labels) {
        ## anchor to the top or bottom-right corner per row; staggered vjust
        ## keeps multi-channel labels from overlapping
        p <- p +
            ggplot2::geom_text(
                ggplot2::aes(
                    y = yval_corner,
                    label = label,
                    colour = nirs_channels,
                    vjust = vjust
                ),
                data = ann,
                x = Inf,
                hjust = 1.05,
                size = 4.5,
                show.legend = FALSE,
                inherit.aes = FALSE
            )
    }

    return(p)
}


#' Build per-panel kinetics marker and label annotations
#'
#' Maps a fitted `mnirs_kinetics` method to its key coefficient marker
#' (`xval`, `yval`) and formatted `label`, one row per `nirs_channel`
#' per interval, for [plot.mnirs_kinetics()]. Marker x-coordinates are the
#' resolved onset plus the method's time coefficient. Labels sit in the corner
#' the fitted curve vacates: `yval_corner` is the top (`Inf`) when the signal
#' falls (`A > B`) and the bottom (`-Inf`) when it rises, with `vjust`
#' staggering stacked labels by channel rank within each interval-corner.
#'
#' @param x An *"mnirs_kinetics"* object from [analyse_kinetics()].
#'
#' @returns A `data.frame` with columns `interval`, `nirs_channels`,
#'   `start_times`, `xval`, `yval`, `label`, `yval_corner`, and `vjust`.
#'
#' @keywords internal
kinetics_annotations <- function(x) {
    ## format numbers as display strings, preserving NA
    fmt <- function(v) {
        out <- rep("NA", length(v))
        ok <- !is.na(v)
        out[ok] <- signif_trailing(v[ok], 3L, "signif")
        return(out)
    }

    coefs <- merge(
        x$coefficients,
        x$interval_times[c("interval", "start_times")],
        by = "interval",
        sort = FALSE
    )

    ## per-method: time offset (x), fitted value (y), and label formatter
    spec <- switch(
        x$method,
        response_time = list(
            offset = "response_time",
            y = NA_character_,
            label = sprintf("response time = %s s", fmt(coefs$response_time))
        ),
        peak_slope = list(
            offset = "peak_slope_time",
            y = "fitted",
            label = sprintf("slope = %s", fmt(coefs$slope))
        ),
        monoexponential = list(
            offset = "MRT",
            y = "MRT_fitted",
            label = sprintf(
                "MRT = %s s\ntau = %s",
                fmt(coefs$MRT),
                fmt(coefs$tau)
            )
        ),
        sigmoidal = list(
            offset = "xmid",
            y = "xmid_fitted",
            label = sprintf(
                "xmid = %s s\nslope = %s",
                fmt(coefs$xmid),
                fmt(coefs$slope)
            )
        )
    )

    ann <- data.frame(
        interval = coefs$interval,
        nirs_channels = coefs$nirs_channels,
        start_times = coefs$start_times,
        xval = coefs$start_times + coefs[[spec$offset]],
        yval = if (is.na(spec$y)) NA_real_ else coefs[[spec$y]],
        label = spec$label,
        stringsAsFactors = FALSE
    )

    ## place labels in the corner the fitted curve vacates: top-right when
    ## the signal falls (A > B), bottom-right when it rises. peak_slope has
    ## no A/B, so fall back to the local slope sign as the trend proxy
    rises <- if (all(c("A", "B") %in% names(coefs))) {
        coefs$B > coefs$A
    } else {
        coefs$slope > 0
    }
    ann$yval_corner <- ifelse(rises, -Inf, Inf)

    ## stagger stacked labels by channel rank within each interval-corner,
    ## growing down from the top and up from the bottom
    rank <- stats::ave(
        seq_len(nrow(ann)),
        ann$interval,
        ann$yval_corner,
        FUN = seq_along
    ) -
        1
    ann$vjust <- ifelse(rises, -0.4 - rank * 1.4, 1.4 + rank * 1.4)

    return(ann)
}


#' Validate and bind a list of mnirs data frames for plotting
#' @inheritParams validate_mnirs
#' @keywords internal
as_plot_data <- function(x, env = rlang::caller_env()) {
    if (length(x) == 0L) {
        cli_abort(c(
            "x" = "{.fn plot.mnirs} must contain at least one \\
            {col_blue('\"mnirs\"')} data frame."
        ), call = env)
    }

    if (any(!vapply(x, is.data.frame, logical(1)))) {
        cli_abort(c(
            "x" = "{.fn plot.mnirs} must contain all {col_blue('\"mnirs\"')} \\
            data frames."
        ), call = env)
    }

    ## validate time_channel is consistent across elements
    time_channels <- vapply(x, \(.df) {
        attr(.df, "time_channel") %||% NA_character_
    }, character(1))

    if (anyNA(time_channels)) {
        cli_abort(c(
            "x" = "All elements of {.fn plot.mnirs} must have a \\
            {.field time_channel} attribute."
        ), call = env)
    }
    if (length(unique(time_channels)) > 1L) {
        cli_abort(c(
            "x" = "All elements of {.fn plot.mnirs} must share the same \\
            {.field time_channel}.",
            "i" = "Found: {.val {unique(time_channels)}}."
        ), call = env)
    }

    ## auto-name unnamed list elements
    if (is.null(names(x))) {
        names(x) <- paste0("interval_", seq_along(x))
    }

    ## length-1 list: unwrap to single data frame
    if (length(x) == 1L) {
        return(x[[1L]])
    }

    ## union of nirs_channels across all elements
    nirs_channels <- unique(unlist(
        lapply(x, attr, "nirs_channels"),
        use.names = FALSE
    ))

    ## pad each element with NA for any missing nirs_channels
    x <- lapply(x, \(.df) {
        .df[setdiff(nirs_channels, names(.df))] <- NA_real_
        .df
    })

    ## add interval column to each element, then row-bind
    x <- Map(\(.df, .nm) {
        .df[["interval"]] <- .nm
        .df
    }, x, names(x))
    plot_data <- do.call(rbind, unname(x))
    plot_data[["interval"]] <- factor(
        plot_data[["interval"]],
        levels = unique(plot_data[["interval"]])
    )
    attr(plot_data, "nirs_channels") <- nirs_channels
    attr(plot_data, "time_channel") <- time_channels[[1L]]

    return(plot_data)
}


#' Custom *{mnirs}* ggplot2 theme
#'
#' A `[ggplot2][ggplot2::ggplot2-package]` theme for display.
#'
#' @param base_size Base font size, given in pts.
#' @param base_family Base font family.
#' @param border Define either a *partial* or *full* border around plots.
#' @param ink Colour for text and lines. *Default* is *"black"*.
#' @param paper Background colour. *Default* is *"white"*.
#' @param accent Accent colour for highlights. *Default* is *"#0080ff"*.
#' @param ... Additional arguments to add to `[ggplot2::theme()]`.
#'
#' @details
#' - `axis.title = element_text(face = "bold")` by *default* Modify to *"plain"*.
#'
#' - `panel.grid.major` & `panel.grid.major` set to blank. Modify to
#'   `= element_line()` for visible grid lines.
#'
#' - `legend.position = "top"` by *default* Modify `"none"` to remove legend
#'   entirely.
#'
#' - `border = "partial"` uses `panel.border = element_blank()` and
#'   `axis.line = element_line()`.
#'
#' - `border = "full"` uses `panel.border = element_rect(colour = "black",`
#'   `linewidth = 1)` and `axis.line = element_line()`.
#'
#' - `base_family = "sans"` by *default*.
#'
#' @returns A [ggplot2][ggplot2::ggplot()] theme object.
#'
#' @seealso [palette_mnirs()], [scale_colour_mnirs()]
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' ## plot example data
#' read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2_left = "SmO2 Live",
#'                       smo2_right = "SmO2 Live(2)"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' ) |>
#'     plot(time_labels = TRUE)
#'
#' @export
theme_mnirs <- function(
    base_size = 14,
    base_family = "sans",
    border = c("partial", "full"),
    ink = "black",
    paper = "white",
    accent = "#0080ff",
    ...
) {
    check_installed("ggplot2", reason = "to plot mNIRS data")
    border <- match.arg(border)
    half_line <- base_size * 0.5

    if (border == "partial") {
        panel.border <- ggplot2::element_blank()
        axis.line <- ggplot2::element_line()
    } else {
        panel.border <- ggplot2::element_rect(colour = "black", linewidth = 1)
        axis.line <- ggplot2::element_blank()
    }

    ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        ink = ink,
        paper = paper,
        accent = accent
    ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = ggplot2::rel(1.2),
                lineheight = 1.1
            ),
            plot.subtitle = ggplot2::element_text(lineheight = 1.1),
            plot.caption = ggplot2::element_text(colour = "grey50"),
            panel.border = panel.border,
            axis.line = axis.line,
            axis.title = ggplot2::element_text(face = "bold"),
            strip.background = ggplot2::element_rect(fill = "grey95"),
            strip.text = ggplot2::element_text(
                margin = ggplot2::margin_auto(t = half_line * 0.5)
            ),
            plot.margin = ggplot2::margin_part(r = base_size, unit = "pt"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "top",
            legend.justification = "right",
            legend.direction = "horizontal",
            legend.margin = ggplot2::margin_auto(t = 1),
            legend.box.spacing = ggplot2::unit(half_line * 0.5, "pt")
        ) +
        ggplot2::theme(...)
}


#' Custom *{mnirs}* colour palette
#'
#' @param ... Either a single numeric specifying the number of colours to
#'   return, or character strings specifying colour names. If empty, all
#'   colours are returned.
#'
#' @returns Named (when selecting by name) or unnamed character vector of
#'   hex colours.
#'
#' @seealso [theme_mnirs()], [scale_colour_mnirs()]
#'
#' @examplesIf rlang::is_installed("scales")
#' scales::show_col(palette_mnirs())
#' scales::show_col(palette_mnirs(2))
#' scales::show_col(palette_mnirs("red", "orange"))
#'
#' @export
palette_mnirs <- function(...) {
    colours <- c(
        ## NIRS location codes
        `light blue` = "#0080ff", ## "VL"
        `dark red` = "#ba2630", ## "FCR"
        `light green` = "#5b8c52", ## "BB" "#7dbf70" alt
        `pink` = "#ff80ff", ## "VM"
        `orange` = "#ff7f00", ## "SCM"
        `dark blue` = "#00468Bff", ## "TA"
        `light red` = "#db5555", ## "ECR"
        `green` = "#42B540FF", ## "DL"
        `purple` = "#9f79ee", ## "RF"
        `brown` = "#8b4726", ## "PS"
        `blue` = "#0000ff", ## "HHb"
        `red` = "#ED0000FF" ## "O2Hb"
    )

    dots <- list(...)

    if (length(dots) == 0L) {
        return(colours)
    }

    ## numeric -> subset by count
    if (length(dots) == 1L && is.numeric(dots[[1L]])) {
        n <- dots[[1L]]
        validate_numeric(n, 1, c(1, Inf), msg1 = "one-element positive")
        if (n <= length(colours)) {
            return(unname(colours[seq_len(n)]))
        }
        ## interpolate if more colours needed, but this probably won't look good!
        return(grDevices::colorRampPalette(colours)(n))
    }

    ## character args -> subset by name
    names <- unlist(dots)

    if (!is.character(names)) {
        ## covers condition of multiple numeric vals
        cli_abort(c(
            "x" = "{.fn palette_mnirs} expects a single numeric value \\
            for the number of colours to return, or character colour names."
        ))
    }

    idx <- match(tolower(names), tolower(names(colours)))
    if (anyNA(idx)) {
        cli_abort(c(
            "x" = "{.fn palette_mnirs} unrecognised colour name{?s}: \\
            {.val {names[is.na(idx)]}}.",
            "i" = "Valid names: {.val {names(colours)}}."
        ))
    }
    return(colours[idx])
}


#' Scales for custom *{mnirs}* palette
#'
#' @param ... Arguments passed to `ggplot2::discrete_scale()`.
#' @param aesthetics A character vector with aesthetic(s) passed to
#'   `ggplot2::discrete_scale()`. *Default* is `"colour"`.
#'
#' @returns A [ggplot2][ggplot2::ggplot()] scale object.
#'
#' @seealso [theme_mnirs()], [palette_mnirs()]
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' ## plot example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2_left = "SmO2 Live",
#'                       smo2_right = "SmO2 Live(2)"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' )
#'
#' ggplot2::ggplot(data, ggplot2::aes(x = time)) +
#'     theme_mnirs() +
#'     scale_colour_mnirs() +
#'     ggplot2::geom_line(ggplot2::aes(y = smo2_left, colour = "smo2_left")) +
#'     ggplot2::geom_line(ggplot2::aes(y = smo2_right, colour = "smo2_right"))
#'
#' @rdname scale_colour_mnirs
#' @export
scale_colour_mnirs <- function(..., aesthetics = "colour") {
    check_installed("ggplot2", reason = "to plot mNIRS data")

    ggplot2::discrete_scale(
        aesthetics = aesthetics,
        palette = palette_mnirs,
        na.value = "grey10",
        ...
    )
}

#' @rdname scale_colour_mnirs
#' @export
scale_color_mnirs <- scale_colour_mnirs

#' @rdname scale_colour_mnirs
#' @export
scale_fill_mnirs <- function(..., aesthetics = "fill") {
    check_installed("ggplot2", reason = "to plot mNIRS data")

    ggplot2::discrete_scale(
        aesthetics = aesthetics,
        palette = palette_mnirs,
        na.value = "grey10",
        ...
    )
}


#' Breaks for time span data
#'
#' Pretty time span breaks for plotting in units of 5, 15, 30, 60 sec, etc.
#' Modified from [scales::breaks_timespan()].
#'
#' @param unit The time unit used to interpret numeric data input (*defaults*
#'   to *"secs"*).
#' @param n Desired number of breaks. You may get slightly more or fewer breaks
#'   than requested.
#'
#' @returns Returns a function for generating breaks.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' x <- 0:120
#' y <- sin(2 * pi * x / 15) + rnorm(length(x), 0, 0.2)
#'
#' ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
#'     theme_mnirs() +
#'     ggplot2::scale_x_continuous(breaks = breaks_timespan()) +
#'     ggplot2::geom_line()
#'
#' @export
breaks_timespan <- function(
    unit = c("secs", "mins", "hours", "days", "weeks"),
    n = 5
) {
    unit <- match.arg(unit)
    force(n)
    function(x) {
        x <- as.numeric(as.difftime(x, units = unit), units = "secs")
        range <- range(x, na.rm = TRUE)
        diff <- range[2L] - range[1L]

        ## scale of time range
        ## define nice steps for each unit
        if (diff <= 5 * 60) {
            scale <- 1 ## sec
            nice_steps <- c(1, 2, 5, 10, 15, 20, 30, 60, 120)
        } else if (diff <= 5 * 3600) {
            scale <- 60 ## min
            nice_steps <- c(1, 2, 5, 10, 15, 20, 30, 60, 120) * 60
        } else if (diff <= 5 * 86400) {
            scale <- 3600 ## hr
            nice_steps <- c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 24) * 3600
        } else {
            scale <- 86400 ## days
            nice_steps <- c(1, 7, 28) * 86400
        }

        ## scale to scale units
        range_scaled <- range / scale
        scaled_steps <- nice_steps / scale

        ## find optimal step size from nice_steps
        target_step <- diff(range_scaled) / n
        best_step <- scaled_steps[which.min(abs(scaled_steps - target_step))]

        ## generate breaks
        breaks_scaled <- seq(
            floor(range_scaled[1L] / best_step) * best_step,
            ceiling(range_scaled[2L] / best_step) * best_step,
            by = best_step
        )

        ## convert back to seconds
        round(as.numeric(as.difftime(breaks_scaled * scale, units = "secs")))
    }
}


#' Format time span data as h:mm:ss
#'
#' Convert numeric time span data to `h:mm:ss` format for pretty plotting.
#' Inspired by [ggplot2::scale_x_time()].
#'
#' @param x A numeric vector.
#'
#' @details
#' If all values are less than 3600 (1 hour), then format is returned as
#'   `mm:ss`. If any value is greater than 3600, format is returned as
#'   `h:mm:ss` with leading zeroes.
#'
#' @returns A character vector the same length as `x`.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' x <- 0:120
#' y <- sin(2 * pi * x / 15) + rnorm(length(x), 0, 0.2)
#'
#' ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
#'     theme_mnirs() +
#'     ggplot2::scale_x_continuous(
#'         breaks = breaks_timespan(),
#'         labels = format_hmmss
#'     ) +
#'     ggplot2::geom_line()
#'
#' @export
format_hmmss <- function(x) {
    # validate_numeric(x)
    x <- as.numeric(x)
    ## logical whether to handle NAs
    handle_na <- anyNA(x)

    if (handle_na) {
        na_info <- preserve_na(x)
        x <- na_info$x_valid
    }

    sign <- ifelse(x < 0, "-", "")
    hrs <- as.integer(abs(x) %/% 3600)
    mins <- as.integer((abs(x) %% 3600) %/% 60)
    secs <- abs(x) %% 60

    ## use fractional seconds format when any non-integer present
    if (any(secs %% 1 != 0, na.rm = TRUE)) {
        secs_fmt <- "%05.2f"
    } else {
        secs_fmt <- "%02d"
        secs <- as.integer(secs)
    }

    hmmss_string <- if (any(hrs > 0, na.rm = TRUE)) {
        sprintf(paste0("%s%d:%02d:", secs_fmt), sign, hrs, mins, secs)
    } else {
        sprintf(paste0("%s%02d:", secs_fmt), sign, mins, secs)
    }

    ## return y to original x length with NAs if handled
    if (handle_na) {
        return(restore_na(hmmss_string, na_info))
    } else {
        return(hmmss_string)
    }
}
