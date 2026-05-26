skip_if_not_installed("ggplot2")

## theme_mnirs() ============================================================
test_that("theme_mnirs returns a ggplot2 theme object", {
    theme_obj <- theme_mnirs()
    expect_s3_class(theme_obj, "theme")
    expect_s3_class(theme_obj, "gg")
})

test_that("theme_mnirs border argument works correctly", {
    partial <- theme_mnirs(border = "partial")
    full <- theme_mnirs(border = "full")
    expect_s3_class(partial$panel.border, "element_blank")
    expect_s3_class(full$panel.border, "element_rect")
})

test_that("theme_mnirs accepts custom colours", {
    custom <- theme_mnirs(ink = "red", paper = "blue", accent = "#ff0000")
    expect_s3_class(custom, "theme")
})

## palette_mnirs() ========================================
test_that("palettes returns correct colour vector", {
    all_colours <- palette_mnirs()
    expect_type(all_colours, "character")
    expect_length(all_colours, 12)
    expect_all_true(grepl("^#[0-9A-Fa-f]{6}", all_colours))
})

test_that("palettes subset by number works", {
    expect_length(palette_mnirs(3), 3)
    expect_length(palette_mnirs(1), 1)
    expect_error(palette_mnirs(2:4), "valid.*numeric")
    expect_error(palette_mnirs(0), "valid.*numeric")
})

test_that("palettes subset by name works", {
    red <- palette_mnirs("red")
    expect_named(red, "red")
    expect_equal(red[["red"]], "#ED0000FF")
    expect_error(palette_mnirs("invalid"), "unrecognised")

    multi <- palette_mnirs("red", "blue")
    expect_length(multi, 2)
    expect_named(multi, c("red", "blue"))
    expect_error(palette_mnirs("red", "invalid"), "unrecognised")

    ## mixed types error
    expect_error(palette_mnirs(TRUE), "expects")
})

test_that("palette_mnirs interpolates when n > 12", {
    many <- palette_mnirs(20)
    expect_length(many, 20)
})

## scale_colour_mnirs() ==================================================
test_that("scale_color_mnirs is an alias for scale_colour_mnirs", {
    expect_identical(scale_color_mnirs, scale_colour_mnirs)
})

test_that("scale_*_mnirs returns a ggplot2 Scale object", {
    expect_s3_class(scale_colour_mnirs(), "Scale")
    expect_s3_class(scale_colour_mnirs(), "ScaleDiscrete")

    expect_s3_class(scale_fill_mnirs(), "Scale")
    expect_s3_class(scale_fill_mnirs(), "ScaleDiscrete")
})

test_that("scale_colour_mnirs uses correct aesthetics", {
    expect_equal(scale_colour_mnirs()$aesthetics, "colour")
    expect_equal(scale_fill_mnirs()$aesthetics, "fill")
})

test_that("scale functions pass through additional arguments", {
    scale <- scale_colour_mnirs(name = "Test")
    expect_equal(scale$name, "Test")
})

test_that("scale functions use palette_mnirs", {
    # Extract the palette function and call it
    expect_equal(scale_colour_mnirs()$palette(5), palette_mnirs(5))
    expect_equal(scale_fill_mnirs()$palette(5), palette_mnirs(5))

    # Test with no argument (all colours)
    expect_equal(scale_colour_mnirs()$palette(), palette_mnirs())

    # Test with character argument
    expect_equal(
        scale_colour_mnirs()$palette("light blue", "dark red"),
        palette_mnirs("light blue", "dark red")
    )

    ## test with na.value
    expect_equal(scale_colour_mnirs()$na.value, "grey10")
    expect_equal(scale_fill_mnirs()$na.value, "grey10")
})

test_that("scale functions work in ggplot2 plots", {
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt, colour = factor(cyl))) +
        ggplot2::geom_point() +
        scale_colour_mnirs()

    expect_s3_class(p, "gg")
    expect_s3_class(p$scales$get_scales("colour"), "ScaleDiscrete")
})

## breaks_timespan() ==================================================
test_that("breaks_timespan returns a function", {
    breaks_fn <- breaks_timespan()
    expect_type(breaks_fn, "closure")
})

test_that(" breaks_timespan corresponds to nice_steps for each scale level", {
    # Test scale = 1 (diff <= 5 * 60)
    nice_steps_sec <- c(1, 2, 5, 10, 15, 20, 30, 60, 120)
    x_sec <- c(0, 150) # 2.5 min range
    breaks_sec <- breaks_timespan("secs", n = 5)(x_sec)
    steps_sec <- unique(diff(breaks_sec))
    expect_true(all(steps_sec %in% nice_steps_sec))
    expect_type(breaks_sec, "double")
    expect_true(all(breaks_sec >= 0 & breaks_sec <= 150))
    expect_true(all.equal(length(breaks_sec), 5, tolerance = 2, scale = 1))

    # Test scale = 60 (5 * 60 < diff <= 5 * 3600)
    nice_steps_min <- c(1, 2, 5, 10, 15, 20, 30, 60, 120) * 60
    x_min <- c(0, 7200) # 2 hour range
    breaks_min <- breaks_timespan("secs", n = 5)(x_min)
    steps_min <- unique(diff(breaks_min))
    expect_true(all(steps_min %in% nice_steps_min))
    expect_type(breaks_min, "double")
    expect_true(all(breaks_min >= 0 & breaks_min <= 7200))
    expect_true(all.equal(length(breaks_min), 5, tolerance = 2, scale = 1))

    # Test scale = 3600 (5 * 3600 < diff <= 5 * 86400)
    nice_steps_hr <- c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 24) * 3600
    x_hr <- c(0, 86400) # 1 day range
    breaks_hr <- breaks_timespan("secs", n = 5)(x_hr)
    steps_hr <- unique(diff(breaks_hr))
    expect_true(all(steps_hr %in% nice_steps_hr))
    expect_type(breaks_hr, "double")
    expect_true(all(breaks_hr >= 0 & breaks_hr <= 86400))
    expect_true(all.equal(length(breaks_hr), 5, tolerance = 2, scale = 1))

    # Test scale = 86400 (diff > 5 * 86400)
    nice_steps_day <- c(1, 7, 28) * 86400
    x_day <- c(0, 86400 * 28) # 28 day range
    breaks_day <- breaks_timespan("secs", n = 5)(x_day)
    steps_day <- unique(diff(breaks_day))
    expect_true(all(steps_day %in% nice_steps_day))
    expect_type(breaks_day, "double")
    expect_true(all(breaks_day >= 0 & breaks_day <= 86400 * 28))
    expect_true(all.equal(length(breaks_day), 5, tolerance = 2, scale = 1))
})

## format_hmmss() ==================================================
test_that("format_hmmss formats time correctly", {
    ## seconds
    expect_equal(format_hmmss(0), "00:00")
    expect_equal(format_hmmss(30), "00:30")
    expect_equal(format_hmmss(90), "01:30")
    expect_equal(format_hmmss(3599), "59:59")

    ## hours
    expect_equal(format_hmmss(3600), "1:00:00")
    expect_equal(format_hmmss(3661), "1:01:01")
    expect_equal(format_hmmss(7325), "2:02:05")

    ## negative
    expect_equal(format_hmmss(-30), "-00:30")
    expect_equal(format_hmmss(-3661), "-1:01:01")
})

test_that("format_hmmss handles vectors", {
    result <- format_hmmss(c(0, 30, 90, 3600))
    expect_length(result, 4)
    expect_equal(result, c("0:00:00", "0:00:30", "0:01:30", "1:00:00"))
})

test_that("format_hmmss handles NA values", {
    result <- format_hmmss(c(30, NA, 90))
    expect_length(result, 3)
    expect_true(is.na(result[2]))
})

test_that("format_hmmss handles fractional seconds", {
    ## scalar fractional
    expect_equal(format_hmmss(2.1), "00:02.10")
    expect_equal(format_hmmss(90.5), "01:30.50")
    expect_equal(format_hmmss(-2.1), "-00:02.10")

    ## fractional crossing hours -> uses h:mm:ss.ff
    expect_equal(format_hmmss(3661.25), "1:01:01.25")

    ## vector with at least one fractional uses fractional format throughout
    expect_equal(
        format_hmmss(seq(2, 7, 2.5)),
        c("00:02.00", "00:04.50", "00:07.00")
    )

    ## fractional vector with NA
    result <- format_hmmss(c(2.1, NA, 4.2))
    expect_length(result, 3)
    expect_true(is.na(result[2]))
    expect_equal(result[c(1, 3)], c("00:02.10", "00:04.20"))
})


## as_plot_data() =============================================
# Helper to create mock mNIRS object
mock_mnirs <- function() {
    df <- data.frame(
        time = 1:10,
        HHb = c(1:6, NA, NA, 9:10),
        O2Hb = c(rep(2, 3), NA, NA, rep(2, 5))
    )
    structure(
        df,
        class = c("mnirs", "data.frame"),
        nirs_channels = c("HHb", "O2Hb"),
        time_channel = "time"
    )
}

test_that("as_plot_data errors on invalid lists", {
    ## empty list
    expect_error(as_plot_data(list()), "at least one")
    ## not a df
    expect_error(as_plot_data(list(mock_mnirs(), "not_a_df")), "must contain all")
})

test_that("as_plot_data errors when element missing time_channel attribute", {
    bad <- structure(
        data.frame(time = 1:3, HHb = 1:3),
        class = c("mnirs", "data.frame"),
        nirs_channels = "HHb"
    )
    expect_error(as_plot_data(list(bad)), "time_channel attribute")
})

test_that("as_plot_data errors when elements have differing time_channel", {
    a <- mock_mnirs()
    b <- structure(
        data.frame(t = 1:10, HHb = 1:10),
        class = c("mnirs", "data.frame"),
        nirs_channels = "HHb",
        time_channel = "t"
    )
    expect_error(as_plot_data(list(a, b)), "same.*time_channel")
})

test_that("as_plot_data unwraps single-element list", {
    x <- mock_mnirs()
    result <- as_plot_data(list(x))
    expect_identical(result, x)
})

test_that("as_plot_data row-binds named list with interval factor", {
    a <- mock_mnirs()
    b <- mock_mnirs()
    result <- as_plot_data(list(pre = a, post = b))
    expect_true("interval" %in% names(result))
    expect_s3_class(result[["interval"]], "factor")
    expect_equal(levels(result[["interval"]]), c("pre", "post"))
    expect_equal(nrow(result), nrow(a) + nrow(b))
    expect_equal(attr(result, "time_channel"), "time")
    expect_equal(attr(result, "nirs_channels"), c("HHb", "O2Hb"))
})

test_that("as_plot_data auto-names unnamed list with sequential integers", {
    a <- mock_mnirs()
    b <- mock_mnirs()
    result <- as_plot_data(list(a, b))
    expect_equal(levels(result[["interval"]]), c("interval_1", "interval_2"))
})

test_that("as_plot_data unions nirs_channels across elements", {
    a <- structure(
        data.frame(time = 1:3, HHb = 1:3),
        class = c("mnirs", "data.frame"),
        nirs_channels = "HHb",
        time_channel = "time"
    )
    b <- structure(
        data.frame(time = 1:3, O2Hb = 4:6),
        class = c("mnirs", "data.frame"),
        nirs_channels = "O2Hb",
        time_channel = "time"
    )
    result <- as_plot_data(x = list(a, b))

    # plot(result)
    expect_equal(attr(result, "nirs_channels"), c("HHb", "O2Hb"))
})

## plot.mnirs() ===============================================
test_that("na.omit removes rows with any NA in nirs_channels", {
    x <- mock_mnirs()

    ## na.omit = FALSE: all 10 rows retained per channel
    p1 <- plot(x, points = TRUE)
    expect_equal(nrow(p1$layers[[1L]]$data), 10L) ## HHb line
    expect_equal(nrow(p1$layers[[3L]]$data), 10L) ## O2Hb line

    p2 <- plot(x, points = TRUE, na.omit = TRUE)
    expect_equal(nrow(p2$layers[[1L]]$data), 8L) ## HHb line
    expect_equal(nrow(p2$layers[[3L]]$data), 8L) ## O2Hb line
})

test_that("time_labels controls x-axis name and formatting", {
    x <- mock_mnirs()

    # With time_labels = FALSE (default)
    p1 <- plot(x)
    expect_true(ggplot2::is_waiver(p1$scales$get_scales("x")$name))
    expect_true(ggplot2::is_waiver(p1$scales$get_scales("x")$labels))

    # With time_labels = TRUE
    p2 <- plot(x, time_labels = TRUE)
    expect_equal(p2$labels$x, "time (h:mm:ss)")
    expect_false(ggplot2::is_waiver(p2$scales$get_scales("x")$labels))
})

test_that("plot.mnirs groups and facets", {
    x <- mock_mnirs()

    # Add a grouping column that we want to facet by
    x$group <- rep(c("A", "B"), each = 5)

    # Create plot
    p <- plot(x)

    # Check that the group column exists in plot data
    expect_true("group" %in% names(p$data))

    # Check that group values are correctly repeated
    expect_equal(p$data$group, x$group)

    # Verify faceting works without error
    expect_no_error(p + ggplot2::facet_wrap(~group))
})

test_that("plot.mnirs works with extract_intervals and faceting", {
    x <- mock_mnirs()

    # Simulate extract_intervals output with interval column
    x$interval <- factor(rep(1:2, each = 5))

    p <- plot(x)

    # Verify interval column preserved as factor
    expect_true("interval" %in% names(p$data))
    expect_s3_class(p$data$interval, "factor")

    # Test faceting works
    p_facet <- p + ggplot2::facet_wrap(~interval)
    expect_s3_class(p_facet, "ggplot")

    # Build plot to ensure no errors during rendering
    expect_no_error(ggplot2::ggplot_build(p_facet))
})

test_that("plot.mnirs uses waiver() for breaks when scales is unavailable", {
    x <- mock_mnirs()

    with_mocked_bindings(
        is_installed = function(pkg, ...) FALSE,
        .package = "rlang",
        {
            p <- plot(x, time_labels = FALSE)
            x_scale <- p$scales$get_scales("x")
            y_scale <- p$scales$get_scales("y")
            ## both break functions fall back to waiver()
            expect_true(ggplot2::is_waiver(x_scale$breaks))
            expect_true(ggplot2::is_waiver(y_scale$breaks))
        }
    )
})

test_that("plot.mnirs works on lists", {
    df_list <- read_mnirs(
        file_path = example_mnirs("train.red"),
        nirs_channels = c(
            smo2_left = "SmO2",
            smo2_right = "SmO2 unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        event_channel = c(lap = "Lap/Event"),
        verbose = FALSE
    ) |>
        resample_mnirs(method = "linear", verbose = FALSE) |>
        extract_intervals(
            start = by_lap(3, 5),
            span = c(-30, 120),
            zero_time = TRUE,
            verbose = FALSE
        )

    expect_type(df_list, "list")
    expect_s3_class(df_list, "mnirs")
    expect_length(df_list, 2)

    ## visual check
    p <- plot(df_list)
    expect_s3_class(p, "ggplot")

    ## facet wrap present for multi-element list
    facet_layers <- Filter(\(l) inherits(l, "FacetWrap"), list(p$facet))
    expect_length(facet_layers, 1)

    ## renders without error
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot.mnirs() returns ggplot2 warnings for missing values", {
    a <- structure(
        data.frame(time = 1:3, HHb = 1:3),
        class = c("mnirs", "data.frame"),
        nirs_channels = "HHb",
        time_channel = "time"
    )
    b <- structure(
        data.frame(time = 1:3, O2Hb = 4:6),
        class = c("mnirs", "data.frame"),
        nirs_channels = "O2Hb",
        time_channel = "time"
    )
    result <- as_plot_data(x = list(a, b))

    w <- tryCatch(
        print(plot(result)),
        warning = \(w) conditionMessage(w)
    )
    expect_match(w, "Removed.*containing missing")
})

test_that("n.breaks overrides y-axis break count", {
    skip_if_not_installed("scales")
    x <- mock_mnirs()

    ## default n = 5
    p_default <- plot(x)
    y_breaks_default <- p_default$scales$get_scales("y")$breaks
    expect_type(y_breaks_default, "closure")

    ## n.breaks alters break count on built scale
    p_custom <- plot(x, n.breaks = 10)
    built <- ggplot2::ggplot_build(p_custom)
    built_default <- ggplot2::ggplot_build(p_default)
    y_breaks_custom <- built$layout$panel_params[[1L]]$y$breaks
    y_breaks_def <- built_default$layout$panel_params[[1L]]$y$breaks
    expect_gt(
        length(y_breaks_custom[!is.na(y_breaks_custom)]),
        length(y_breaks_def[!is.na(y_breaks_def)])
    )
})

test_that("breaks overrides x-axis breaks directly", {
    x <- mock_mnirs()
    custom <- seq(2.1, 8.4, 2.1)

    p <- plot(x, breaks = custom)
    built <- ggplot2::ggplot_build(p)
    x_breaks <- built$layout$panel_params[[1L]]$x$breaks
    expect_equal(x_breaks[!is.na(x_breaks)], custom)
})

test_that("breaks overrides x-axis with time_labels = TRUE", {
    x <- mock_mnirs()
    custom <- seq(2.1, 8.4, 2.1)

    p <- plot(x, time_labels = TRUE, breaks = custom)
    built <- ggplot2::ggplot_build(p)
    x_breaks <- built$layout$panel_params[[1L]]$x$breaks
    expect_equal(x_breaks[!is.na(x_breaks)], custom)
})

test_that("plot.mnirs moxy.perfpro works", {
    df <- read_mnirs(
        file_path = example_mnirs("moxy_ramp"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )

    ## visual check
    plot <- plot(df, na.omit = TRUE, time_labels = TRUE, n.breaks = 8)
    expect_s3_class(plot, "ggplot")
})
