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
    for (ch in channels[-1]) {
        df[[ch]] <- x
    }
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

## biexponential excursion-recovery (per test-analyse_biexponential.R)
make_biexp <- function(A = 70, B = 45, B2 = 60, channels = "smo2", n = 121) {
    set.seed(1)
    t <- seq(0, n - 1, length.out = n)
    df <- setNames(
        data.frame(t, biexponential(t, A, B, 5, B2, 40) + rnorm(n, 0, 0.5)),
        c("time", channels[1])
    )
    for (ch in channels[-1]) {
        df[[ch]] <- biexponential(t, A + 5, B + 5, 5, B2 + 5, 40) +
            rnorm(n, 0, 0.5)
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

## monoexponential plus linear drift from the onset -tau * log(1 - 0.98) =
## 31.3 (per create_expdrift_data); no TD so the response sits at t = 0
make_expdrift <- function(channels = "smo2", n = 120) {
    set.seed(42)
    t <- seq(0, n - 1, length.out = n)
    df <- setNames(
        data.frame(
            t,
            exponential_drift(t, 70, 40, 8, 0.2, 0.98) + rnorm(n, 0, 0.3)
        ),
        c("time", channels[1])
    )
    for (ch in channels[-1]) {
        df[[ch]] <- exponential_drift(t, 75, 45, 8, 0.2, 0.98) +
            rnorm(n, 0, 0.3)
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

kin_biexp <- function(
    A = 70,
    B = 45,
    B2 = 60,
    channels = "smo2",
    faceted = FALSE
) {
    analyse_kinetics(
        as_input(make_biexp(A, B, B2, channels), faceted),
        nirs_channels = channels,
        method = "biexponential",
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

kin_expdrift <- function(channels = "smo2", faceted = FALSE) {
    analyse_kinetics(
        as_input(make_expdrift(channels), faceted),
        nirs_channels = channels,
        method = "exponential_drift",
        use_TD = FALSE,
        verbose = TRUE
    )
}

## geom class of each plot layer, e.g. "GeomLine", "GeomPoint"
layer_geoms <- function(p) {
    vapply(p$layers, \(l) class(l$geom)[[1L]], character(1))
}

## dotted line layers drawn by `components = TRUE` (vline is not GeomLine)
comp_layers <- function(p) {
    Filter(\(l) {
        inherits(l$geom, "GeomLine") &&
            identical(l$aes_params$linetype, "dotted")
    }, p$layers)
}


## kinetics_annotations() ================================================
## marker rows have finite `xval`; label rows have non-empty `label`
ann_markers <- function(ann) ann[is.finite(ann$xval), ]
ann_labels <- function(ann) ann[nzchar(ann$label), ]

test_that("kinetics_annotations returns expected structure", {
    x <- kin_peak_slope(
        channels = c("smo2_left", "smo2_right"),
        faceted = TRUE
    )
    ann <- kinetics_annotations(x)

    expect_s3_class(ann, "data.frame")
    expect_named(ann, c(
        "interval", "nirs_channels", "xval", "yval", "label", "vjust"
    ))
    ## one marker row per channel per interval (2 channels x 2 intervals)
    markers <- ann_markers(ann)
    expect_equal(nrow(markers), 4L)
    expect_true(all(markers$label == ""))
    expect_true(all(is.na(markers$vjust)))
    ## one label row per line (slope + time) per channel per interval
    labels <- ann_labels(ann)
    expect_equal(nrow(labels), 8L)
    expect_true(all(is.infinite(labels$xval)))
    expect_true(all(is.infinite(labels$yval)))
})

test_that("kinetics_annotations xval is onset plus method offset", {
    ## offset column differs by method
    ps <- kin_peak_slope()
    expect_equal(
        ann_markers(kinetics_annotations(ps))$xval,
        ps$interval_times$start_times + ps$coefficients$peak_slope_time
    )

    me <- kin_monoexp()
    expect_equal(
        ann_markers(kinetics_annotations(me))$xval,
        me$interval_times$start_times + me$coefficients$MRT
    )

    be <- kin_biexp()
    expect_equal(
        ann_markers(kinetics_annotations(be))$xval,
        c(
            be$interval_times$start_times + be$coefficients$MRT,
            be$interval_times$start_times + be$coefficients$texc
        )
    )

    ed <- kin_expdrift()
    expect_equal(
        ann_markers(kinetics_annotations(ed))$xval,
        c(
            ed$interval_times$start_times + ed$coefficients$MRT,
            ed$interval_times$start_times + ed$coefficients$texc
        )
    )

    sg <- kin_sigmoidal()
    expect_equal(
        ann_markers(kinetics_annotations(sg))$xval,
        sg$interval_times$start_times + sg$coefficients$xmid
    )

    rt <- kin_response_time()
    expect_equal(
        ann_markers(kinetics_annotations(rt))$xval,
        rt$interval_times$start_times + rt$coefficients$response_time
    )
})

test_that("kinetics_annotations formats method-specific labels", {
    expect_match(
        ann_labels(kinetics_annotations(kin_response_time()))$label,
        "\\% response = .+ s$"
    )
    expect_match(
        ann_labels(kinetics_annotations(kin_peak_slope()))$label[[1L]],
        "^slope = "
    )
    ## `use_TD = FALSE` -> `MRT` is redundant with `tau` and is omitted
    expect_match(
        ann_labels(kinetics_annotations(kin_monoexp()))$label,
        "^tau = .+ s$"
    )
    ## one row per line, in coefficient order
    biexp <- ann_labels(kinetics_annotations(kin_biexp()))$label
    expect_equal(sub(" = .*", "", biexp), c("tau", "texc", "tau2"))
    expect_match(biexp, " s$")
    sigm <- ann_labels(kinetics_annotations(kin_sigmoidal()))$label
    expect_equal(sub(" = .*", "", sigm), c("slope", "xmid"))
    expect_match(sigm[[1L]], " /s$")
})

test_that("kinetics_annotations rounds time to 1 decimal, slope to 3 sigfig", {
    ## `kinetics_annotations()` reads only `method` and `coefficients`
    coefs <- data.frame(
        interval = "a",
        nirs_channels = "smo2",
        start_time = 0,
        A = 0,
        B = 1,
        TD = 6.58,
        tau = 2494.4,
        MRT = 2500.98,
        MRT_fitted = 1
    )
    ann <- kinetics_annotations(list(
        method = "monoexponential",
        coefficients = coefs
    ))
    expect_equal(
        ann_labels(ann)$label,
        c("TD = 6.6 s", "tau = 2494 s", "MRT = 2501 s")
    )

    coefs <- data.frame(
        interval = "a",
        nirs_channels = "smo2",
        start_time = 0,
        slope = 0.012345,
        peak_slope_time = 12.345,
        fitted = 1
    )
    ann <- kinetics_annotations(list(method = "peak_slope", coefficients = coefs))
    expect_equal(
        ann_labels(ann)$label,
        c("slope = 0.0123 /s", "time = 12.3 s")
    )
})

test_that("kinetics_annotations biexponential marks MRT and the fitted excursion", {
    x <- kin_biexp()
    markers <- ann_markers(kinetics_annotations(x))
    expect_equal(nrow(markers), 2L)
    expect_equal(
        markers$yval,
        c(x$coefficients$MRT_fitted, x$coefficients$texc_fitted)
    )
    ## both points sit below the baseline for a downward response
    expect_true(all(markers$yval < x$coefficients$A))
})

test_that("kinetics_annotations places label in the vacated corner", {
    ## rising signal (B > A) -> bottom corner (-Inf)
    rise <- kin_monoexp(A = 50, B = 80)
    expect_true(all(ann_labels(kinetics_annotations(rise))$yval == -Inf))

    ## falling signal (A > B) -> top corner (Inf)
    fall <- kin_monoexp(A = 80, B = 50)
    expect_true(all(ann_labels(kinetics_annotations(fall))$yval == Inf))
})

test_that("kinetics_annotations biexponential corner follows the plateau", {
    ## net trend is the plateau B2 against the baseline A

    ## fall-recover: plateau below baseline -> falls -> top corner (Inf)
    fall <- kin_biexp(B = 45, B2 = 55)
    expect_true(all(fall$coefficients$B2 < fall$coefficients$A))
    expect_true(all(ann_labels(kinetics_annotations(fall))$yval == Inf))

    ## rise-overshoot: plateau above baseline -> rises -> bottom corner (-Inf)
    rise <- kin_biexp(B = 95, B2 = 85)
    expect_true(all(rise$coefficients$B2 > rise$coefficients$A))
    expect_true(all(ann_labels(kinetics_annotations(rise))$yval == -Inf))
})

test_that("kinetics_annotations peak_slope corner follows the slope sign", {
    ## no asymptote to trend on, so direction comes from the fitted slope
    x <- kin_peak_slope()
    labels <- ann_labels(kinetics_annotations(x))
    expect_true(all(labels$yval == ifelse(x$coefficients$slope > 0, -Inf, Inf)))
})

test_that("kinetics_annotations mixed-direction channels share one corner", {
    ## opposite-direction channels still get a single corner per panel, with
    ## labels staggered rather than split; a tie falls back to the top corner
    set.seed(13)
    t <- 0:59
    df <- data.frame(
        time = t,
        up = monoexponential(t, 50, 80, 5, 5) + rnorm(60, 0, 0.5),
        down = monoexponential(t, 80, 50, 5, 5) + rnorm(60, 0, 0.5)
    )
    d <- create_mnirs_data(
        df,
        nirs_channels = c("up", "down"),
        time_channel = "time",
        sample_rate = 1,
        interval_times = 0
    )
    x <- analyse_kinetics(
        d,
        nirs_channels = c("up", "down"),
        method = "monoexponential",
        use_TD = FALSE,
        verbose = FALSE
    )
    labels <- ann_labels(kinetics_annotations(x))
    expect_true(all(labels$yval == Inf))
    expect_equal(labels$vjust, c(1.8, 3.4))
})

test_that("kinetics_annotations stacks label lines inward from the corner", {
    ## two rising channels share the bottom corner; single-line labels stack
    ## 1.6 text heights apart in channel order, half a gap above the border,
    ## first channel on top
    x <- kin_monoexp(channels = c("smo2_left", "smo2_right"))
    labels <- ann_labels(kinetics_annotations(x))
    expect_equal(labels$nirs_channels, c("smo2_left", "smo2_right"))
    expect_equal(labels$vjust, c(-2.4, -0.8))

    ## multi-line labels continue the same stack across channels
    x <- kin_biexp(channels = c("smo2_left", "smo2_right"))
    labels <- ann_labels(kinetics_annotations(x))
    expect_equal(labels$vjust, 0.2 + 1.6 * seq_len(nrow(labels)))
})

test_that("kinetics_annotations omits NA coefficients from labels", {
    ## too few points for a 3-param fit -> NA coefficients -> empty label
    d <- make_monoexp(n = 3)
    x <- suppressWarnings(analyse_kinetics(
        d,
        nirs_channels = "smo2",
        method = "monoexponential",
        use_TD = FALSE,
        verbose = FALSE
    ))
    ann <- kinetics_annotations(x)
    expect_equal(nrow(ann_labels(ann)), 0L)
    ## the marker row remains with an NA position
    expect_equal(ann$label, "")
    expect_true(is.na(ann$xval))
})

test_that("kinetics_annotations adds a marker row per key point", {
    ## exponential_drift: MRT and texc markers per channel
    x <- kin_expdrift(channels = c("smo2_left", "smo2_right"))
    markers <- ann_markers(kinetics_annotations(x))
    expect_equal(nrow(markers), 4L)

    texc <- markers[3:4, ]
    expect_equal(texc$nirs_channels, x$coefficients$nirs_channels)
    expect_equal(
        texc$xval,
        x$interval_times$start_times + x$coefficients$texc
    )
    expect_equal(texc$yval, x$coefficients$texc_fitted)

    ## single key-point methods have one marker per channel
    expect_equal(nrow(ann_markers(kinetics_annotations(kin_monoexp()))), 1L)
})


## plot.mnirs_kinetics() =================================================
test_that("plot.mnirs_kinetics returns a ggplot and renders for each method", {
    objs <- list(
        peak_slope = kin_peak_slope(),
        response_time = kin_response_time(),
        monoexponential = kin_monoexp(),
        biexponential = kin_biexp(),
        sigmoidal = kin_sigmoidal()
    )
    lapply(objs, \(x) {
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

test_that("fitted overlay densifies to 100 points only when < 100 samples", {
    dashed_data <- \(p) {
        Filter( \(l) {
            inherits(l$geom, "GeomLine") &&
                identical(l$aes_params$linetype, "dashed")
        }, p$layers)[[1L]]$data
    }
    ## 60 samples < 100: curve re-predicted on a 100-point grid
    x <- kin_monoexp()
    p <- plot(x, markers = FALSE, labels = FALSE)
    expect_equal(nrow(dashed_data(p)), 100L)
    ## 121 samples >= 100: observed fit-window rows kept as-is
    x <- kin_biexp()
    p <- plot(x, markers = FALSE, labels = FALSE)
    expect_equal(nrow(dashed_data(p)), sum(is.finite(x$data[[1L]]$smo2_fitted)))
})

test_that("fitted overlay and key points use a darker shade of the channel colour", {
    x <- kin_peak_slope()
    ## layers: observed line, dashed fitted line, onset vline, key point
    p <- plot(x, labels = FALSE)
    lum <- \(i) sum(grDevices::col2rgb(ggplot2::layer_data(p, i)$colour[1L]))
    expect_lt(lum(2L), lum(1L))
    expect_lt(lum(4L), lum(1L))
    ## legend keys keep the observed channel colour
    expect_equal(
        ggplot2::layer_data(p, 1L)$colour[1L],
        unname(palette_mnirs(1L))
    )
})

test_that("fitted = FALSE drops the parametric fitted layers", {
    x <- kin_sigmoidal()
    n_on <- length(plot(x, markers = FALSE, labels = FALSE)$layers)
    n_off <- length(
        plot(x, fitted = FALSE, markers = FALSE, labels = FALSE)$layers
    )
    expect_gt(n_on, n_off)
})

test_that("biexponential draws the fitted curve and an excursion key-point", {
    x <- kin_biexp()
    geoms <- layer_geoms(plot(x, labels = FALSE))
    ## base signal line + dashed fitted line
    expect_equal(sum(geoms == "GeomLine"), 2L)
    ## onset vline + single excursion marker
    expect_true("GeomVline" %in% geoms)
    expect_equal(sum(geoms == "GeomPoint"), 1L)

    ## fitted = FALSE drops the dashed overlay only
    off <- layer_geoms(plot(x, fitted = FALSE, labels = FALSE))
    expect_equal(sum(off == "GeomLine"), 1L)
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


## components overlay ====================================================
test_that("components overlays the biexponential model terms", {
    x <- kin_biexp()
    p <- plot(x, components = TRUE, markers = FALSE, labels = FALSE)
    comps <- comp_layers(p)
    expect_length(comps, 2L)

    ## fast (A -> B) and slow (B -> B2) terms sum to the fitted curve
    d <- comps[[1L]]$data
    fitted <- x$data[[1L]]$smo2_fitted
    expect_true(all.equal(
        d$comp1 + d$comp2 - x$coefficients$B,
        fitted[is.finite(fitted)],
        tolerance = 1e-6,
        scale = 1
    ))

    ## omitted by default
    expect_length(comp_layers(plot(x, markers = FALSE, labels = FALSE)), 0L)
})

test_that("components draws the exponential_drift drift line from the drift onset", {
    x <- kin_expdrift()
    p <- plot(x, components = TRUE, markers = FALSE, labels = FALSE)
    comps <- comp_layers(p)
    expect_length(comps, 2L)

    ## monoexponential term spans the full fit window
    d1 <- comps[[1L]]$data
    expect_equal(nrow(d1), sum(is.finite(x$data[[1]]$smo2_fitted)))

    # drift term restricted to t >= the drift onset, linear at the fitted slope
    d2 <- comps[[2L]]$data
    expect_lt(nrow(d2), nrow(d1))
    t_rel <- d2$time - x$interval_times$start_times[[1L]]
    cf <- x$coefficients
    onset <- expdrift_onset(cf$tau, cf$drift_frac, if (is.finite(cf$TD)) cf$TD)
    expect_true(all(t_rel >= onset))
    expect_equal(diff(d2$comp2), rep(x$coefficients$slope, nrow(d2) - 1L))
})

test_that("components tolerates channels named after coefficients", {
    ## a channel named `slope` collides with an exponential_drift
    ## coefficient; overlay frame must not join coefficient columns
    single <- kin_expdrift(channels = "slope")
    p1 <- plot(single, components = TRUE, markers = FALSE, labels = FALSE)
    expect_length(comp_layers(p1), 2L)
    expect_no_error(ggplot2::ggplot_build(p1))

    ## drift line reflects the fitted slope, not the data channel
    d2 <- comp_layers(p1)[[2L]]$data
    expect_equal(diff(d2$comp2), rep(single$coefficients$slope, nrow(d2) - 1L))

    faceted <- kin_expdrift(channels = "slope", faceted = TRUE)
    p2 <- plot(faceted, components = TRUE, markers = FALSE, labels = FALSE)
    expect_length(comp_layers(p2), 2L)
    expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("components draws per channel and per facet", {
    x <- kin_biexp(channels = c("smo2_left", "smo2_right"), faceted = TRUE)
    p <- plot(x, components = TRUE, markers = FALSE, labels = FALSE)
    comps <- comp_layers(p)
    ## two terms per channel
    expect_length(comps, 4L)
    ## overlay rows carry the interval for facet placement
    expect_true("interval" %in% names(comps[[1L]]$data))
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("components is ignored for methods without a component spec", {
    x <- kin_monoexp()
    expect_length(comp_layers(plot(x, components = TRUE)), 0L)
    expect_equal(
        length(plot(x, components = TRUE)$layers),
        length(plot(x)$layers)
    )
})

test_that("components skips channels with a failed fit", {
    ## direction mismatch fails the fit -> NA coefficients
    x <- suppressWarnings(analyse_kinetics(
        make_biexp(),
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        direction = "positive",
        verbose = FALSE
    ))
    expect_true(is.na(x$coefficients$A))
    p <- plot(x, components = TRUE, markers = FALSE, labels = FALSE)
    expect_length(comp_layers(p), 0L)
    expect_no_error(ggplot2::ggplot_build(p))
})
