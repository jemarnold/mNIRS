skip_if_not_installed("ggplot2")

## fixtures: real mnirs_kinetics objects per method ======================
## Build via analyse_kinetics() so tests exercise the true object shape.
## Synthetic-data patterns mirror those in test-analyse_kinetics.R.

## deterministic sine data (per create_kinetics_data); columns are fixed
## smo2_left / smo2_right so channels must be drawn from those names
make_sine <- function(n = 50, sample_rate = 10) {
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    df <- data.frame(
        time = t,
        smo2_left = sin(t) * 10 + 50,
        smo2_right = cos(t) * 10 + 50
    )
    create_mnirs_data(
        df,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = "time",
        sample_rate = sample_rate,
        interval_times = 0
    )
}

## ramp-then-plateau (per create_response_time_data)
make_ramp <- function(A = 0, B = 20, channels = "smo2") {
    x <- c(rep(A, 5), seq(A, B, length.out = 10), rep(B, 5))
    t <- seq_along(x) - 5 ## t = 0 at end of baseline
    df <- setNames(data.frame(t, x), c("time", channels[1]))
    for (ch in channels[-1]) df[[ch]] <- x
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = 1,
        interval_times = 0
    )
}

## monoexponential rise/fall (per create_monoexp_data)
make_monoexp <- function(A = 50, B = 80, channels = "smo2", n = 60) {
    set.seed(13)
    t <- seq(0, n - 1, length.out = n)
    df <- setNames(
        data.frame(t, monoexponential(t, A, B, 5, 5) + rnorm(n, 0, 0.5)),
        c("time", channels[1])
    )
    for (ch in channels[-1]) {
        df[[ch]] <- monoexponential(t, A + 5, B + 5, 5, 5) + rnorm(n, 0, 0.5)
    }
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = 1,
        interval_times = 0
    )
}

## sigmoidal (per create_sigmoidal_data)
make_sigmoidal <- function(channels = "smo2", n = 60) {
    set.seed(13)
    t <- seq(0, n - 1, length.out = n)
    df <- setNames(
        data.frame(t, logistic(t, 10, 100, 30, 4) + rnorm(n, 0, 2)),
        c("time", channels[1])
    )
    for (ch in channels[-1]) {
        df[[ch]] <- logistic(t, 15, 105, 30, 4) + rnorm(n, 0, 2)
    }
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = 1,
        interval_times = 0
    )
}

## wrap single data frame or a 2-element list for faceted vs single panels
as_input <- function(d, faceted) if (faceted) list(A = d, B = d) else d

kin_peak_slope <- function(channels = "smo2_left", faceted = FALSE) {
    analyse_kinetics(
        as_input(make_sine(), faceted),
        nirs_channels = channels,
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )
}

kin_response_time <- function(channels = "smo2", faceted = FALSE) {
    analyse_kinetics(
        as_input(make_ramp(channels = channels), faceted),
        nirs_channels = channels,
        method = "response_time",
        direction = "positive",
        verbose = FALSE
    )
}

kin_monoexp <- function(A = 50, B = 80, channels = "smo2", faceted = FALSE) {
    analyse_kinetics(
        as_input(make_monoexp(A, B, channels), faceted),
        nirs_channels = channels,
        method = "monoexponential",
        use_TD = FALSE,
        verbose = FALSE
    )
}

kin_sigmoidal <- function(channels = "smo2", faceted = FALSE) {
    analyse_kinetics(
        as_input(make_sigmoidal(channels), faceted),
        nirs_channels = channels,
        method = "sigmoidal",
        shape = "symmetric",
        verbose = FALSE
    )
}

## geom class of each plot layer, e.g. "GeomLine", "GeomPoint"
layer_geoms <- function(p) {
    vapply(p$layers, \(l) class(l$geom)[[1L]], character(1))
}


## kinetics_annotations() ================================================
test_that("kinetics_annotations returns expected structure", {
    x <- kin_peak_slope(
        channels = c("smo2_left", "smo2_right"),
        faceted = TRUE
    )
    ann <- kinetics_annotations(x)

    expect_s3_class(ann, "data.frame")
    expect_named(ann, c(
        "interval", "nirs_channels", "start_times",
        "xval", "yval", "label", "yval_corner", "vjust"
    ))
    ## one row per channel per interval (2 channels x 2 intervals)
    expect_equal(nrow(ann), 4L)
})

test_that("kinetics_annotations xval is onset plus method offset", {
    ## offset column differs by method
    ps <- kin_peak_slope()
    expect_equal(
        kinetics_annotations(ps)$xval,
        ps$interval_times$start_times + ps$coefficients$peak_slope_time
    )

    me <- kin_monoexp()
    expect_equal(
        kinetics_annotations(me)$xval,
        me$interval_times$start_times + me$coefficients$MRT
    )

    sg <- kin_sigmoidal()
    expect_equal(
        kinetics_annotations(sg)$xval,
        sg$interval_times$start_times + sg$coefficients$xmid
    )

    rt <- kin_response_time()
    expect_equal(
        kinetics_annotations(rt)$xval,
        rt$interval_times$start_times + rt$coefficients$response_time
    )
})

test_that("kinetics_annotations formats method-specific labels", {
    expect_match(
        kinetics_annotations(kin_response_time())$label,
        "^response time = .+ s$"
    )
    expect_match(kinetics_annotations(kin_peak_slope())$label, "^slope = ")
    expect_match(
        kinetics_annotations(kin_monoexp())$label,
        "MRT = .+ s\ntau = "
    )
    expect_match(
        kinetics_annotations(kin_sigmoidal())$label,
        "xmid = .+ s\nslope = "
    )
})

test_that("kinetics_annotations places label in the vacated corner", {
    ## rising signal (B > A) -> bottom corner (-Inf)
    rise <- kin_monoexp(A = 50, B = 80)
    expect_true(all(kinetics_annotations(rise)$yval_corner == -Inf))

    ## falling signal (A > B) -> top corner (Inf)
    fall <- kin_monoexp(A = 80, B = 50)
    expect_true(all(kinetics_annotations(fall)$yval_corner == Inf))
})

test_that("kinetics_annotations peak_slope corner follows slope sign", {
    ## peak_slope has no A/B, so trend proxy is the local slope sign
    x <- kin_peak_slope()
    ann <- kinetics_annotations(x)
    expect_equal(ann$yval_corner, ifelse(x$coefficients$slope > 0, -Inf, Inf))
})

test_that("kinetics_annotations staggers stacked labels by channel rank", {
    ## two rising channels share the bottom corner -> distinct vjust per rank
    x <- kin_monoexp(channels = c("smo2_left", "smo2_right"))
    ann <- kinetics_annotations(x)
    expect_equal(length(unique(ann$vjust)), 2L)
    expect_equal(ann$vjust, c(-0.4, -1.8))
})

test_that("kinetics_annotations formats NA coefficients as 'NA'", {
    ## too few points for a 3-param fit -> NA coefficients
    d <- make_monoexp(n = 3)
    x <- suppressWarnings(analyse_kinetics(
        d,
        nirs_channels = "smo2",
        method = "monoexponential",
        use_TD = FALSE,
        verbose = FALSE
    ))
    ann <- kinetics_annotations(x)
    expect_match(ann$label, "MRT = NA s")
    expect_match(ann$label, "tau = NA")
})


## plot.mnirs_kinetics() =================================================
test_that("plot.mnirs_kinetics returns a ggplot and renders for each method", {
    objs <- list(
        peak_slope = kin_peak_slope(),
        response_time = kin_response_time(),
        monoexponential = kin_monoexp(),
        sigmoidal = kin_sigmoidal()
    )
    lapply(objs, function(x) {
        p <- plot(x)
        expect_s3_class(p, "ggplot")
        expect_no_error(ggplot2::ggplot_build(p))
    })
})

test_that("fitted overlay draws a dashed line per channel for parametric fits", {
    x <- kin_peak_slope(channels = c("smo2_left", "smo2_right"))
    geoms <- layer_geoms(plot(x, markers = FALSE, labels = FALSE))
    ## base signal line per channel (2) + fitted dashed line per channel (2)
    expect_equal(sum(geoms == "GeomLine"), 4L)
})

test_that("fitted = FALSE drops the parametric fitted layers", {
    x <- kin_sigmoidal()
    n_on <- length(plot(x, markers = FALSE, labels = FALSE)$layers)
    n_off <- length(
        plot(x, fitted = FALSE, markers = FALSE, labels = FALSE)$layers
    )
    expect_gt(n_on, n_off)
})

test_that("response_time has no fitted curve or points", {
    x <- kin_response_time()
    geoms <- layer_geoms(plot(x, markers = FALSE, labels = FALSE))
    ## points come from markers, not fitted; disabled here
    expect_false("GeomPoint" %in% geoms)
    ## only the base channel line; no extra dashed fitted line
    expect_equal(sum(geoms == "GeomLine"), 1L)
})

test_that("markers toggles the onset vline and parametric key-point", {
    x <- kin_peak_slope()
    on <- layer_geoms(plot(x, fitted = FALSE, labels = FALSE, markers = TRUE))
    off <- layer_geoms(plot(x, fitted = FALSE, labels = FALSE, markers = FALSE))
    expect_true("GeomVline" %in% on)
    expect_false("GeomVline" %in% off)
    ## parametric method adds a key-point marker
    expect_true("GeomPoint" %in% on)
})

test_that("markers adds the vline and key-points for response_time", {
    x <- kin_response_time()
    on <- layer_geoms(plot(x, fitted = FALSE, labels = FALSE, markers = TRUE))
    off <- layer_geoms(
        plot(x, fitted = FALSE, labels = FALSE, markers = FALSE)
    )
    expect_true("GeomVline" %in% on)
    expect_false("GeomVline" %in% off)
    ## response_time key points (baseline, response, extreme) are markers
    expect_true("GeomPoint" %in% on)
    expect_false("GeomPoint" %in% off)
})

test_that("labels toggles the geom_text layer", {
    x <- kin_monoexp()
    on <- layer_geoms(plot(x, fitted = FALSE, markers = FALSE, labels = TRUE))
    off <- layer_geoms(plot(x, fitted = FALSE, markers = FALSE, labels = FALSE))
    expect_true("GeomText" %in% on)
    expect_false("GeomText" %in% off)
})

test_that("plot.mnirs_kinetics facets multi-interval objects only", {
    multi <- plot(kin_peak_slope(faceted = TRUE))
    single <- plot(kin_peak_slope(faceted = FALSE))
    expect_length(Filter(\(l) inherits(l, "FacetWrap"), list(multi$facet)), 1L)
    expect_length(Filter(\(l) inherits(l, "FacetWrap"), list(single$facet)), 0L)
})

test_that("plot.mnirs_kinetics forwards ... to plot.mnirs", {
    p <- plot(kin_peak_slope(), time_labels = TRUE)
    expect_equal(p$labels$x, "time (mm:ss)")
})

test_that("plot.mnirs_kinetics sets label size", {
    p <- plot(kin_peak_slope(), label_size = 5)
    text_layer <- p$layers[[which(layer_geoms(p) == "GeomText")[1L]]]
    expect_equal(text_layer$aes_params$size, 5)
})
