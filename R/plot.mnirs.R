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
    channel_map <- attr(x, "channel_map")

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

    ## sort facets by appearance
    if ("interval" %in% names(x) && !is.factor(x[["interval"]])) {
        x[["interval"]] <- factor(x[["interval"]], unique(x[["interval"]]))
    }

    ## build base plot with axis configuration
    plot <- ggplot2::ggplot(x) +
        ggplot2::aes(x = .data[[time_channel]]) +
        theme_mnirs() +
        ggplot2::labs(x = x_name, y = "mNIRS") +
        ggplot2::scale_x_continuous(
            breaks = x_breaks,
            labels = x_labels,
            expand = ggplot2::expansion(mult = 0.02)
        ) +
        ggplot2::scale_y_continuous(
            breaks = y_breaks,
            expand = ggplot2::expansion(mult = 0.02)
        ) +
        scale_colour_mnirs(
            name = NULL,
            guide = ggplot2::guide_legend(
                override.aes = list(linewidth = 1)
            )
        )

    ## add one geom per channel, restricted to the panels declaring it
    layers <- lapply(nirs_channels, function(channel) {
        keep <- if (is.null(channel_map)) {
            TRUE
        } else {
            x[["interval"]] %in% channel_map[[channel]]
        }
        if (na.omit) {
            keep <- keep & is.finite(x[[channel]])
        }
        ch_data <- x[keep, , drop = FALSE]
        ch_aes <- ggplot2::aes(y = .data[[channel]], colour = channel)
        c(
            list(ggplot2::geom_line(ch_aes, data = ch_data)),
            if (points) {
                list(ggplot2::geom_point(ch_aes, data = ch_data, size = 2))
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
#'   parametric methods (`"peak_slope"`, `"monoexponential"`,
#'   `"biexponential"`, `"sigmoidal"`). `"response_time"` has no fitted curve.
#' @param markers Logical. Default is `TRUE`; draws a dotted vertical line at
#'   the response onset (`start_time`) and key coefficient points.
#' @param labels Logical. Default is `TRUE`; annotates each panel with the key
#'   coefficient value(s) for the fitted method.
#' @param ... Additional arguments.
#'
#' @details
#' Accepts some arguments in `...`, such as `label_size` passed to
#' [ggplot2::geom_text()]. Also accepts args passed to [plot.mnirs()], such as
#' `points`, `time_labels`, `nrow`, `ncol`, or `scales`.
#'
#' A method with no annotation spec in [kinetics_annotations()] plots the
#' observed signal and fitted curve only, without markers or labels.
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

    ## channels with a fitted column present in the bound frame
    fit_ch <- nirs[paste0(nirs, "_fitted") %in% names(plot_data)]

    ## appearance-ordered facet levels: ggplot2 unions interval levels across
    ## layers, so any character interval column re-sorts facets alphabetically
    ## (interval_10 before _2). factor the source frames once; all overlay
    ## frames derive from them and inherit the factor through merges
    if (faceted) {
        lvls <- unique(x$interval_times$interval)
        x$interval_times$interval <- factor(x$interval_times$interval, lvls)
        x$coefficients$interval <- factor(x$coefficients$interval, lvls)
    }

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
        curved <- x$method != "peak_slope"
        p <- p +
            lapply(fit_ch, \(.ch) {
            fcol <- paste0(.ch, "_fitted")
            d <- plot_data[is.finite(plot_data[[fcol]]), , drop = FALSE]
            ## curved fits: re-predict on a dense time grid so fitted lines
            ## plot smoothly when an interval has < 100 fit-window samples
            sp <- split(d, if (faceted) d$interval else rep_len(1L, nrow(d)), drop = TRUE)
            if (curved && any(vapply(sp, nrow, 0L) < 100L)) {
                mods <- if (faceted) x$model[names(sp)] else x$model[1L]
                d <- do.call(rbind, Map(\(.d, .m) {
                    r <- range(.d[[time_channel]])
                    t <- seq(r[1L], r[2L], length.out = max(100L, nrow(.d)))
                    ## replicate the first row so interval and start_times
                    ## carry over without rebuilding the frame
                    dd <- .d[rep(1L, length(t)), , drop = FALSE]
                    dd[[time_channel]] <- t
                    ## time symbol from the model formula: a channel
                    ## colliding with a model parameter is fit aliased
                    m <- .m[[.ch]]
                    ## self-start models predict with a gradient attribute
                    dd[[fcol]] <- as.vector(stats::predict(
                        m,
                        newdata = setNames(
                            data.frame(t - .d$start_times[[1L]]),
                            as.character(stats::formula(m)[[3L]][[2L]])
                        )
                    ))
                    return(dd)
                }, sp, mods))
            }
            ggplot2::geom_line(
                ggplot2::aes(y = .data[[fcol]], colour = .ch),
                data = d,
                linetype = "dashed",
                linewidth = 1,
                show.legend = FALSE
            )
        })
    }

    ## model component overlay ================================
    ## undocumented `components = TRUE`: reconstruct the model terms from
    ## natural-scale coefficients over the fitted rows. comp1 is the
    ## primary monoexponential; comp2 is the secondary term: the
    ## biexponential slow phase (B1 to B2) clocked from the onset, or the
    ## exponential_drift linear drift from (texc, texc_fitted)
    comp_methods <- c("biexponential", "exponential_drift")
    if (isTRUE(list(...)[["components"]]) && x$method %in% comp_methods) {
        p <- p +
            lapply(fit_ch, \(.ch) {
            fcol <- paste0(.ch, "_fitted")
            d <- plot_data[is.finite(plot_data[[fcol]]), , drop = FALSE]
            if (nrow(d) == 0L) {
                return(NULL)
            }
            ## coefficient row aligned to each fitted row; overlay frame
            ## holds only time, interval, and component columns, so user
            ## channel names can never collide with coefficient names
            cf <- x$coefficients[x$coefficients$nirs_channels == .ch, ]
            co <- cf[
                if (faceted) match(d$interval, cf$interval) else rep(1L, nrow(d)),
            ]
            t_rel <- d[[time_channel]] - co$start_time
            ## TD NA marks a fit with no time delay; those fits only keep
            ## rows from the onset, so TD = 0 is equivalent
            TD <- ifelse(is.finite(co$TD), co$TD, 0)
            cd <- d[c(time_channel, if (faceted) "interval")]
            ## terms follow the model that fit each row; a coefficient
            ## absent from the schema reads as NA
            model <- co$model %||% rep(x$method, nrow(co))
            g <- \(.nm) co[[.nm]] %||% NA_real_
            biexp <- model == "biexponential"
            cd$comp1 <- ifelse(
                biexp,
                monoexponential(t_rel, g("A"), g("B1"), g("tau1"), TD),
                monoexponential(t_rel, g("A"), g("B"), g("tau"), TD)
            )
            cd$comp2 <- ifelse(
                biexp,
                monoexponential(t_rel, g("B1"), g("B2"), g("tau2"), TD),
                ifelse(
                    model == "exponential_drift" & t_rel >= g("texc"),
                    g("texc_fitted") + g("slope") * (t_rel - g("texc")),
                    NA_real_
                )
            )

            comp_line <- \(.col, .d = cd) ggplot2::geom_line(
                ggplot2::aes(y = .data[[.col]], colour = .ch),
                data = .d,
                linetype = "dotted",
                linewidth = 0.5,
                show.legend = FALSE
            )
            list(
                comp_line("comp1"),
                comp_line("comp2", cd[is.finite(cd$comp2), , drop = FALSE])
            )
        })
    }

    ## per-interval-per-channel markers and labels ============
    ## methods with no annotation spec plot the fitted curve alone
    ann <- kinetics_annotations(x)

    if (is.null(ann)) {
        return(p)
    }

    if (!faceted) {
        ann$interval <- NULL
    }

    if (markers) {
        ## dotted onset line at the resolved start_time per interval
        p <- p +
            ggplot2::geom_vline(
                ggplot2::aes(xintercept = .data$start_times),
                data = x$interval_times,
                linetype = "dotted",
                colour = "grey50"
            )

        if (x$method == "response_time") {
            ## response and extreme (fitted values after onset) as points
            p <- p +
                lapply(fit_ch, \(.ch) {
                fcol <- paste0(.ch, "_fitted")
                post <- is.finite(plot_data[[fcol]]) &
                    plot_data[[time_channel]] > plot_data$start_times
                key_point(
                    ggplot2::aes(y = .data[[fcol]], colour = .ch),
                    plot_data[post, , drop = FALSE]
                )
            })
            ## baseline as a single point at the onset (start_time, A);
            ## unique() drops duplicate rows from multiple fractions
            base_pts <- unique(x$coefficients[
                c("interval", "nirs_channels", "start_time", "A")
            ])
            if (!faceted) {
                base_pts$interval <- NULL
            }
            p <- p +
                key_point(
                    ggplot2::aes(
                        x = .data$start_time,
                        y = .data$A,
                        colour = .data$nirs_channels
                    ),
                    base_pts,
                    inherit.aes = FALSE
                )
        } else {
            ## single key-point marker for parametric methods
            p <- p +
                key_point(
                    ggplot2::aes(
                        x = .data$xval,
                        y = .data$yval,
                        colour = .data$nirs_channels
                    ),
                    ann[is.finite(ann$xval), , drop = FALSE],
                    inherit.aes = FALSE
                )
        }
    }

    if (labels) {
        ## one text row per label line anchored at the panel corner;
        ## `vjust` stacks lines inward so channels do not overlap
        p <- p +
            ggplot2::geom_text(
                ggplot2::aes(
                    x = .data$xval,
                    y = .data$yval,
                    label = .data$label,
                    colour = .data$nirs_channels,
                    vjust = .data$vjust
                ),
                data = ann[nzchar(ann$label), , drop = FALSE],
                hjust = 1.05,
                size = list(...)[["label_size"]] %||% 3.5,
                show.legend = FALSE,
                inherit.aes = FALSE
            )
    }

    return(p)
}


#' Build per-panel kinetics marker and label annotations
#'
#' Maps a fitted `mnirs_kinetics` method to its key coefficient markers
#' (`xval`, `yval`) and formatted label lines for [plot.mnirs_kinetics()].
#' Marker rows are one per `nirs_channel` per key point per interval, with
#' x-coordinates the resolved onset plus the method's time coefficient.
#' Label rows are one per label line, anchored (`xval = Inf`,
#' `yval = -Inf` or `Inf`) at the corner of the panel's right edge vacated by
#' the fitted data, judged from the coefficient direction. `vjust` stacks
#' the lines inward from the corner in channel order within each interval.
#'
#' @param x An *"mnirs_kinetics"* object from [analyse_kinetics()].
#'
#' @returns A `data.frame` with columns `interval`, `nirs_channels`,
#'   `xval`, `yval`, `label`, and `vjust`. Marker rows have an empty `label`
#'   and `NA` `vjust`; label rows have infinite `xval`/`yval`. Rows are
#'   annotated by the model that fit them (the `model` coefficient column
#'   where the method has a fallback chain, else the method).
#'   `NULL` for a method with no annotation spec, in which case
#'   [plot.mnirs_kinetics()] draws the fitted curve alone.
#'
#' @keywords internal
kinetics_annotations <- function(x) {
    coefs <- x$coefficients

    ## one label line per coefficient, `NA` when missing (e.g. `TD` for
    ## channels fitted without a time delay). `keep` drops lines that are
    ## redundant, e.g. `MRT` equals `tau` without `TD`. values show 1 decimal
    ## at most and 3 significant figures below that, whole numbers in full;
    ## "fg" without the "#" flag drops trailing zeros
    line <- \(f, v, keep = TRUE, decimals = 1L) {
        v <- signif_whole(round(v, decimals), 3L)
        ifelse(
            is.na(v) | !keep,
            NA_character_,
            sprintf(f, trimws(formatC(v, digits = 3L, format = "fg")))
        )
    }

    ## per-channel list of label lines, omitting `NA` lines
    label <- \(...) {
        apply(cbind(...), 1L, \(l) l[!is.na(l)], simplify = FALSE)
    }

    ## per-method: time offsets (x), fitted values (y), and label lines.
    ## `offset`/`y` are parallel vectors of coefficient names, one marker
    ## point per pair
    annotation_spec <- \(method, coefs) switch(
        method,
        response_time = list(
            offset = "response_time",
            y = "fitted",
            ## fraction-specific labels, e.g. "50% response = 7.9 s";
            ## outer sprintf resolves the percentage, leaving `%s` for line()
            label = label(line(
                sprintf("%g%%%% response = %%s s", coefs$fraction * 100),
                coefs$response_time
            ))
        ),
        peak_slope = list(
            offset = "peak_slope_time",
            y = "fitted",
            label = label(
                line("slope = %s /s", coefs$slope, decimals = Inf),
                line("time = %s s", coefs$peak_slope_time)
            )
        ),
        monoexponential = list(
            offset = "MRT",
            y = "MRT_fitted",
            label = label(
                line("TD = %s s", coefs$TD),
                line("tau = %s s", coefs$tau),
                line("MRT = %s s", coefs$MRT, keep = !is.na(coefs$TD))
            )
        ),
        biexponential = list(
            offset = c("MRT", "texc"),
            y = c("MRT_fitted", "texc_fitted"),
            label = label(
                line("TD = %s s", coefs$TD),
                line("tau1 = %s s", coefs$tau1),
                line("texc = %s s", coefs$texc),
                line("tau2 = %s s", coefs$tau2)
            )
        ),
        exponential_drift = list(
            offset = c("MRT", "texc"),
            y = c("MRT_fitted", "texc_fitted"),
            label = label(
                line("TD = %s s", coefs$TD),
                line("tau = %s s", coefs$tau),
                line("texc = %s s", coefs$texc),
                line("slope = %s /s", coefs$slope, decimals = Inf)
            )
        ),
        sigmoidal = list(
            offset = "xmid",
            y = "xmid_fitted",
            label = label(
                line("slope = %s /s", coefs$slope, decimals = Inf),
                line("xmid = %s s", coefs$xmid)
            )
        )
    )

    ## rows are annotated by the model that fit them: the per-row `model`
    ## where the method has a fallback chain, else the method. marker rows
    ## carry no label or stacking; label rows anchor at the panel corner.
    ## `row` keeps coefficient order across the model groups
    models <- coefs$model %||% rep(x$method, nrow(coefs))
    ann <- do.call(rbind, lapply(split(seq_len(nrow(coefs)), models), \(.i) {
        cf <- coefs[.i, , drop = FALSE]
        spec <- annotation_spec(models[[.i[[1L]]]], cf)
        if (is.null(spec)) {
            return(NULL)
        }
        markers <- Map(\(off, y) {
            data.frame(
                row = .i,
                interval = cf$interval,
                nirs_channels = cf$nirs_channels,
                xval = cf$start_time + cf[[off]],
                yval = cf[[y]],
                label = "",
                stringsAsFactors = FALSE
            )
        }, spec$offset, spec$y)
        ## channels with all-NA fits contribute no label lines
        n <- lengths(spec$label)
        labels <- data.frame(
            row = rep(.i, n),
            interval = rep(cf$interval, n),
            nirs_channels = rep(cf$nirs_channels, n),
            xval = rep(Inf, sum(n)),
            yval = rep(NA_real_, sum(n)),
            label = as.character(unlist(spec$label)),
            stringsAsFactors = FALSE
        )
        rbind(do.call(rbind, markers), labels)
    }))

    ## methods without an annotation spec degrade to a curve-only plot
    if (is.null(ann)) {
        return(NULL)
    }
    ## marker rows first, then label rows in coefficient order
    is_lab <- nzchar(ann$label)
    ann <- rbind(ann[!is_lab, ], ann[is_lab, ][order(ann$row[is_lab]), ])
    is_lab <- nzchar(ann$label)

    ## signed response direction: fitted slope sign (peak_slope, sigmoidal),
    ## otherwise plateau minus baseline; the plateau is `B2` or `B`,
    ## whichever the row's model reports
    dir <- if (is.null(coefs[["A"]])) {
        coefs[["slope"]]
    } else {
        plateau <- Reduce(
            \(.x, .y) ifelse(is.na(.x), .y, .x),
            coefs[intersect(c("B2", "B"), names(coefs))]
        )
        plateau - coefs[["A"]]
    }

    ## one corner per panel: labels anchor to the right edge, so use the half
    ## the fitted responses vacate. sign-sum majority across channels decides;
    ## ties and all-NA fits fall back to the top corner
    # fmt: skip
    rises <- stats::ave(
        sign(dir),
        coefs$interval,
        FUN = \(s) sum(s, na.rm = TRUE)
    ) > 0
    ann$yval[is_lab] <- ifelse(rises[ann$row[is_lab]], -Inf, Inf)

    ## stack lines inward from the corner, half a line-gap from the border,
    ## keeping top-to-bottom order in both corners. `vjust` is in single-line
    ## text heights (~0.7 font size), so 1.6 approximates geom_text's 1.2
    ## lineheight
    ann$vjust <- NA_real_
    if (any(is_lab)) {
        lab <- ann[is_lab, ]
        idx <- stats::ave(seq_along(lab$label), lab$interval, FUN = seq_along)
        rev_idx <- stats::ave(idx, lab$interval, FUN = rev)
        ann$vjust[is_lab] <- ifelse(
            lab$yval < 0, 0.8 - 1.6 * rev_idx, 0.2 + 1.6 * idx
        )
    }
    ann$row <- NULL
    rownames(ann) <- NULL
    return(ann)
}


#' Validate and bind a list of mnirs data frames for plotting
#' @inheritParams validate_mnirs
#'
#' @returns For a single-element list, that element unchanged. Otherwise a
#'   row-bound `data.frame` with an `interval` factor column, carrying
#'   attributes `nirs_channels` (the union across elements), `time_channel`,
#'   and `channel_map` — a named list mapping each channel to the interval
#'   names whose source element declares it, so [plot.mnirs()] draws each
#'   channel only in its own panels.
#'
#' @keywords internal
as_plot_data <- function(x, env = rlang::caller_env()) {
    if (length(x) == 0L) {
        cli_abort(c(
            "x" = "{.fn plot.mnirs} must contain at least one \\
            {.cls mnirs} data frame."
        ), call = env)
    }

    if (any(!vapply(x, is.data.frame, logical(1)))) {
        cli_abort(c(
            "x" = "{.fn plot.mnirs} must contain all {.cls mnirs} \\
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

    ## declared channels per element, before padding obscures them
    declared <- lapply(x, \(.df) attr(.df, "nirs_channels") %||% character(0))

    ## union of nirs_channels across all elements
    nirs_channels <- unique(unlist(declared, use.names = FALSE))

    ## invert to channel -> intervals declaring it, so each geom draws only
    ## in the panels whose source element declares that channel
    channel_map <- lapply(nirs_channels, \(.ch) {
        names(declared)[vapply(declared, \(.d) .ch %in% .d, logical(1))]
    })
    names(channel_map) <- nirs_channels

    ## pad each element with NA for any column it lacks, so elements with
    ## asymmetrical channels (and their derived `_fitted` columns) row-bind
    all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))
    all_cols <- union(all_cols, nirs_channels)
    x <- lapply(x, \(.df) {
        .df[setdiff(all_cols, names(.df))] <- NA_real_
        .df[all_cols]
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
    attr(plot_data, "channel_map") <- channel_map

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
#' scales::show_col(palette_mnirs("red", "blue", "green"))
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
    return(unname(colours[idx]))
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
