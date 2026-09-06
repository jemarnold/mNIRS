## sigmoidal_drift() ================================================

shapes <- c("symmetric", "gompertz", "gompertz_left")
sigmoid_fn <- list(
    symmetric = logistic,
    gompertz = gompertz,
    gompertz_left = gompertz_left
)

test_that("sigmoidal_drift() reduces to the sigmoid without drift", {
    t <- seq(0, 100, by = 0.5)
    lapply(shapes, \(.s) {
        expect_equal(
            sigmoidal_drift(t, 10, 100, 40, 4, 0, 0.95, shape = .s),
            sigmoid_fn[[.s]](t, 10, 100, 40, 4)
        )
    })
})

test_that("sigdrift_onset() is where the sigmoid reaches the fraction", {
    p <- 0.95
    check <- \(A, B, slope) {
        lapply(shapes, \(.s) {
            onset <- sigdrift_onset(A, B, 40, slope, p, .s)
            expect_gt(onset, 40)
            expect_equal(sigmoid_fn[[.s]](onset, A, B, 40, slope), A + p * (B - A))
        })
    }
    check(10, 100, 4)
    ## falling response
    check(100, 10, -4)

    ## the Gompertz form reaches furthest on its slow tail, gompertz_left
    ## nearest on its fast tail
    g <- sigdrift_onset(10, 100, 40, 4, p, "gompertz")
    s <- sigdrift_onset(10, 100, 40, 4, p, "symmetric")
    gl <- sigdrift_onset(10, 100, 40, 4, p, "gompertz_left")
    expect_gt(g, s)
    expect_gt(s, gl)

    ## vectorised over the fraction
    expect_equal(
        sigdrift_onset(10, 100, 40, 4, c(0.9, 0.95), "symmetric"),
        c(sigdrift_onset(10, 100, 40, 4, 0.9, "symmetric"), s)
    )
    expect_error(sigdrift_onset(10, 100, 40, 4, 3, "symmetric"), "in \\(0, 1\\)")
})

test_that("sigdrift_texc() is where the drift rate overtakes the sigmoid", {
    p <- 0.95
    lapply(shapes, \(.s) {
        onset <- sigdrift_onset(10, 100, 40, 4, p, .s)
        ## a drift at least as fast as the peak rate takes over at the onset
        expect_equal(sigdrift_texc(10, 100, 40, 4, 4, p, .s), onset)
        expect_equal(sigdrift_texc(10, 100, 40, 4, -5, p, .s), onset)
        ## a fast drift relative to the rate at the onset is floored there
        expect_equal(sigdrift_texc(10, 100, 40, 4, 2, p, .s), onset)

        ## a slow drift takes over past the onset, where |S'| = |slope_B|
        texc <- sigdrift_texc(10, 100, 40, 4, -0.04, p, .s)
        expect_gt(texc, onset)
        h <- 1e-4
        dS <- diff(sigmoid_fn[[.s]](texc + c(-h, h), 10, 100, 40, 4)) / (2 * h)
        expect_equal(dS, 0.04, tolerance = 1e-4)

        ## falling response
        texc_neg <- sigdrift_texc(100, 10, 40, -4, 0.04, p, .s)
        expect_equal(texc_neg, texc)
    })
})

test_that("sigmoidal_drift() drift is hinged at the onset", {
    t <- seq(0, 100, by = 0.5)
    lapply(shapes, \(.s) {
        onset <- sigdrift_onset(10, 100, 40, 4, 0.95, .s)
        base <- sigmoidal_drift(t, 10, 100, 40, 4, 0, 0.95, shape = .s)
        drift <- sigmoidal_drift(t, 10, 100, 40, 4, -0.4, 0.95, shape = .s)

        ## zero before the onset, linear after
        expect_equal(drift[t <= onset], base[t <= onset])
        expect_equal(
            drift[t > onset] - base[t > onset],
            -0.4 * (t[t > onset] - onset)
        )
    })
})


## SSsigmoidal_drift() ==============================================

test_that("SSsigmoidal_drift() recovers parameters for every shape", {
    t <- seq(0, 119)
    lapply(shapes, \(.s) {
        set.seed(13)
        x <- sigmoidal_drift(t, 10, 100, 40, 4, -0.4, 0.95, shape = .s) +
            rnorm(length(t), 0, 1)
        data <- data.frame(t, x)
        formula <- substitute(
            x ~ SSsigmoidal_drift(
                t, A, B, xmid, slope, slope_B,
                drift_fraction = 0.95, shape = .s
            ),
            list(.s = .s)
        )
        model <- suppressWarnings(nls(
            formula,
            data = data,
            algorithm = "port",
            control = nls.control(warnOnly = TRUE)
        ))
        cf <- coef(model)
        expect_named(cf, c("A", "B", "xmid", "slope", "slope_B"))
        expect_true(all.equal(cf[["A"]], 10, tolerance = 3, scale = 1))
        expect_true(all.equal(cf[["B"]], 100, tolerance = 3, scale = 1))
        expect_true(all.equal(cf[["xmid"]], 40, tolerance = 2, scale = 1))
        expect_true(all.equal(cf[["slope"]], 4, tolerance = 1, scale = 1))
        expect_true(all.equal(cf[["slope_B"]], -0.4, tolerance = 0.1, scale = 1))
    })
})

test_that("SSsigmoidal_drift() excludes fixed parameters from estimation", {
    set.seed(13)
    t <- seq(0, 119)
    x <- sigmoidal_drift(t, 10, 100, 40, 4, -0.4, 0.95) +
        rnorm(length(t), 0, 1)
    data <- data.frame(t, x)

    model <- suppressWarnings(nls(
        x ~ SSsigmoidal_drift(
            t, A = 10, B, xmid, slope, slope_B, drift_fraction = 0.95
        ),
        data = data,
        algorithm = "port",
        control = nls.control(warnOnly = TRUE)
    ))
    expect_named(coef(model), c("B", "xmid", "slope", "slope_B"))
    expect_true(all.equal(coef(model)[["B"]], 100, tolerance = 3, scale = 1))
    ## predicts on the fixed value and the default symmetric shape
    expect_equal(
        as.numeric(predict(model, data.frame(t = 0:5))),
        sigmoidal_drift(0:5, 10, coef(model)[["B"]], coef(model)[["xmid"]],
            coef(model)[["slope"]], coef(model)[["slope_B"]], 0.95)
    )
})


## sigdrift_start() =================================================

test_that("sigdrift_start() seeds near the truth and holds fixed values", {
    set.seed(13)
    t <- seq(0, 119)
    x <- sigmoidal_drift(t, 10, 100, 40, 4, -0.4, 0.95) +
        rnorm(length(t), 0, 1)

    start <- sigdrift_start(x, t)
    expect_named(
        start,
        c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    )
    expect_true(all.equal(start[["A"]], 10, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["B"]], 100, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["xmid"]], 40, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["slope_B"]], -0.4, tolerance = 0.15, scale = 1))
    expect_equal(start[["drift_fraction"]], 0.95)

    fixed <- sigdrift_start(x, t, fixed = list(A = 12, slope_B = -0.5))
    expect_equal(fixed[["A"]], 12)
    expect_equal(fixed[["slope_B"]], -0.5)
})

test_that("sigdrift_start() seeds zero drift without tail support", {
    ## five points inside the seeded onset, so the tail has fewer than two
    t <- seq(30, 50, by = 5)
    x <- logistic(t, 10, 100, 40, 4)
    start <- sigdrift_start(x, t)
    expect_equal(start[["slope_B"]], 0)
})


## analyse_sigmoidal_drift() ========================================

## helper: rising sigmoid with a negative trailing drift; symmetric onset
## at about 57
create_sigdrift_data <- function(
    A = 10,
    B = 100,
    xmid = 40,
    slope = 4,
    slope_B = -0.4,
    drift_fraction = 0.95,
    shape = "symmetric",
    n = 120,
    sample_rate = 1,
    noise_sd = 1,
    channels = "smo2",
    seed = 42
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    df <- data.frame(time = t)
    ## successive channels are offset by 5 units
    df[channels] <- lapply(seq_along(channels) - 1L, \(.i) {
        sigmoidal_drift(
            t, A + 5 * .i, B + 5 * .i, xmid, slope, slope_B,
            drift_fraction, shape
        ) + rnorm(n, 0, noise_sd)
    })

    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate
    )
}


test_that("analyse_sigmoidal_drift() returns correct structure and recovers parameters", {
    result <- analyse_sigmoidal_drift(
        create_sigdrift_data(),
        nirs_channels = "smo2",
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "A", "B", "xmid", "slope", "texc",
        "slope_B", "drift_fraction", "xmid_fitted", "texc_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    model <- attr(result, "model")$smo2
    expect_s3_class(model, "nls")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_equal(attr(result, "channel_args")$drift_fraction, 0.95)
    expect_equal(attr(result, "channel_args")$shape, "symmetric")

    ## the onset fraction is held, never estimated
    expect_named(coef(model), c("A", "B", "xmid", "slope", "slope_B"))
    expect_true(all.equal(result$A, 10, tolerance = 3, scale = 1))
    expect_true(all.equal(result$B, 100, tolerance = 3, scale = 1))
    expect_true(all.equal(result$xmid, 40, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, 4, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope_B, -0.4, tolerance = 0.1, scale = 1))
    expect_equal(result$drift_fraction, 0.95)
    expect_true(attr(result, "diagnostics")$r2 > 0.99)
    expect_equal(attr(result, "diagnostics")$n_params, 5L)

    ## derived columns follow the fitted coefficients
    expect_equal(
        result$texc,
        sigdrift_texc(
            result$A, result$B, result$xmid, result$slope, result$slope_B,
            0.95, "symmetric"
        )
    )
    expect_gt(result$texc, result$xmid)
    expect_equal(
        c(result$xmid_fitted, result$texc_fitted),
        as.numeric(predict(model, data.frame(time = c(result$xmid, result$texc))))
    )
})

test_that("analyse_sigmoidal_drift() fits every shape", {
    lapply(shapes, \(.s) {
        result <- analyse_sigmoidal_drift(
            create_sigdrift_data(shape = .s),
            nirs_channels = "smo2",
            shape = .s,
            verbose = FALSE
        )
        expect_equal(attr(result, "channel_args")$shape, .s)
        expect_true(all.equal(result$xmid, 40, tolerance = 2, scale = 1))
        expect_true(all.equal(result$slope_B, -0.4, tolerance = 0.1, scale = 1))
        ## the stored model predicts with its shape
        model <- attr(result, "model")$smo2
        expect_equal(
            as.numeric(predict(model, data.frame(time = c(0, 40, 119)))),
            sigmoidal_drift(
                c(0, 40, 119), result$A, result$B, result$xmid, result$slope,
                result$slope_B, 0.95, shape = .s
            )
        )
    })
})

test_that("analyse_sigmoidal_drift() drift_fraction resolves per channel", {
    data <- create_sigdrift_data(channels = c("smo2", "hhb"))

    result <- analyse_sigmoidal_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_fraction = list(smo2 = 0.98, hhb = 0.9),
        verbose = FALSE
    )
    expect_equal(result$drift_fraction, c(0.98, 0.9))
    expect_equal(attr(result, "channel_args")$drift_fraction, c(0.98, 0.9))
    ## a larger fraction pushes the onset, and so texc, further from xmid
    expect_gt(result$texc[[1L]] - result$xmid[[1L]],
        result$texc[[2L]] - result$xmid[[2L]])
    expect_false("drift_fraction" %in% names(coef(attr(result, "model")$smo2)))

    ## an omitted channel takes the formal default
    result_part <- analyse_sigmoidal_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_fraction = list(smo2 = 0.98),
        verbose = FALSE
    )
    expect_equal(result_part$drift_fraction, c(0.98, 0.95))
})

test_that("analyse_sigmoidal_drift() validates drift_fraction", {
    data <- create_sigdrift_data()
    lapply(list(0.5, 1, 0, 1.2, "0.95", c(0.9, 0.95)), \(.f) {
        expect_error(
            analyse_sigmoidal_drift(
                data, nirs_channels = "smo2", drift_fraction = .f, verbose = FALSE
            ),
            "drift_fraction.*must be a valid one-element"
        )
    })
})

test_that("analyse_sigmoidal_drift() holds fixed parameters", {
    data <- create_sigdrift_data()

    result <- analyse_sigmoidal_drift(
        data,
        nirs_channels = "smo2",
        fix = list(A = 10, slope_B = -0.4),
        verbose = FALSE
    )
    expect_equal(result$A, 10)
    expect_equal(result$slope_B, -0.4)
    expect_named(coef(attr(result, "model")$smo2), c("B", "xmid", "slope"))
    expect_equal(attr(result, "diagnostics")$n_params, 3L)
    expect_true(all.equal(result$B, 100, tolerance = 3, scale = 1))

    ## the onset fraction is not a fixable parameter
    expect_error(
        analyse_sigmoidal_drift(
            data, nirs_channels = "smo2", fix = list(drift_fraction = 0.9)
        ),
        "not recognised"
    )
})

test_that("analyse_sigmoidal_drift() fails on too few observations", {
    custom_name <- create_sigdrift_data(n = 5, noise_sd = 0.1)
    expect_warning(
        result <- analyse_sigmoidal_drift(custom_name, "smo2"),
        "fit failed for.*smo2.*custom_name.*5 observations for 5 free"
    )
    expect_true(all(is.na(result[c("A", "slope_B", "texc", "xmid_fitted")])))
    expect_null(attr(result, "model")$smo2)
})


## direction ========================================================

test_that("analyse_sigmoidal_drift() direction = 'negative' matches auto on falling data", {
    data <- create_sigdrift_data(A = 100, B = 10, slope = -4)

    result_auto <- analyse_sigmoidal_drift(
        data, nirs_channels = "smo2", direction = "auto", verbose = FALSE
    )
    result_neg <- analyse_sigmoidal_drift(
        data, nirs_channels = "smo2", direction = "negative", verbose = FALSE
    )

    expect_equal(result_auto$A, result_neg$A)
    expect_equal(result_auto$slope, result_neg$slope)
    expect_true(result_auto$B < result_auto$A)
    expect_true(result_auto$slope < 0)
    expect_gt(result_auto$texc, result_auto$xmid)
})

test_that("analyse_sigmoidal_drift() direction = 'positive' rejects falling fit", {
    data <- create_sigdrift_data(A = 100, B = 10, slope = -4)

    expect_warning(
        result <- analyse_sigmoidal_drift(
            data, nirs_channels = "smo2", direction = "positive"
        ),
        "satisfy"
    )
    expect_true(all(is.na(result[c("A", "B", "slope", "slope_B", "texc")])))
})


## model fallback ===================================================

test_that("analyse_kinetics() keeps a supported drift", {
    data <- create_sigdrift_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        verbose = FALSE
    )
    forced <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        model_fallback = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(result$method, "sigmoidal_drift")
    expect_equal(
        names(cf)[1:4], c("interval", "nirs_channels", "start_time", "model")
    )
    expect_equal(cf$model, "sigmoidal_drift")
    expect_true(all.equal(cf$slope_B, -0.4, tolerance = 0.1, scale = 1))
    expect_false(any(grepl("fell back to", result$warnings$message)))
    expect_equal(cf, forced$coefficients)
})

test_that("analyse_kinetics() falls back from a negligible drift", {
    data <- create_sigdrift_data(slope_B = 0, shape = "gompertz")

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "sigmoidal_drift",
            shape = "gompertz"
        ),
        "fell back to the sigmoidal model"
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_equal(cf$model, "sigmoidal")
    expect_named(coef(model), c("A", "B", "xmid", "slope"))
    expect_equal(cf$xmid, coef(model)[["xmid"]])
    expect_true(all(is.na(
        cf[c("slope_B", "drift_fraction", "texc", "texc_fitted")]
    )))
    expect_false(is.na(cf$xmid_fitted))
    expect_equal(result$diagnostics$n_params, 4L)
    expect_equal(
        result$data[[1L]]$smo2_fitted,
        as.vector(predict(model))
    )
    ## the reduced fit keeps the shape
    expect_equal(result$channel_args$shape, "gompertz")
    expect_true(any(grepl("Drift amplitude", result$warnings$message)))

    ## the raw fit is kept on request
    forced <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        shape = "gompertz",
        model_fallback = FALSE,
        verbose = FALSE
    )
    expect_equal(forced$coefficients$model, "sigmoidal_drift")
    expect_false(is.na(forced$coefficients$slope_B))
})

test_that("sigmoidal_drift fallback resolves per channel with fix carried", {
    data <- create_sigdrift_data(channels = c("smo2", "hhb"))
    data$hhb <- create_sigdrift_data(slope_B = 0, seed = 1)$smo2

    result <- analyse_kinetics(
        data,
        nirs_channels = c(smo2, hhb),
        method = "sigmoidal_drift",
        fix = list(A = 10),
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$nirs_channels, c("smo2", "hhb"))
    expect_equal(cf$model, c("sigmoidal_drift", "sigmoidal"))
    expect_equal(cf$A, c(10, 10))
    expect_named(coef(result$model[[1L]]$hhb), c("B", "xmid", "slope"))
    fell <- result$warnings[grepl("fell back to", result$warnings$message), ]
    expect_equal(fell$nirs_channels, "hhb")
})


## analyse_kinetics() dispatch ======================================

test_that("analyse_kinetics.sigmoidal_drift dispatches via method aliases", {
    data <- create_sigdrift_data()
    aliases <- c("sigmoid_drift", "sig-drift", "Logistic Drift", "gompertz_drift")
    lapply(aliases, \(.m) {
        result <- analyse_kinetics(
            data, nirs_channels = "smo2", method = .m, verbose = FALSE
        )
        expect_equal(result$method, "sigmoidal_drift")
    })
})

test_that("analyse_kinetics.sigmoidal_drift passes shape, drift_fraction, and fix", {
    result <- analyse_kinetics(
        create_sigdrift_data(shape = "gompertz_left"),
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        shape = "gompertz_left",
        drift_fraction = 0.9,
        fix = list(B = 100),
        verbose = FALSE
    )
    cf <- result$coefficients
    expect_equal(cf$model, "sigmoidal_drift")
    expect_equal(cf$B, 100)
    expect_equal(cf$drift_fraction, 0.9)
    expect_equal(result$channel_args$shape, "gompertz_left")
    expect_equal(result$channel_args$drift_fraction, 0.9)
    expect_false("B" %in% names(coef(result$model[[1L]]$smo2)))
})

test_that("print and plot methods handle sigmoidal_drift", {
    result <- analyse_kinetics(
        create_sigdrift_data(),
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        verbose = FALSE
    )
    out <- capture.output(print(result))
    expect_match(out[[2L]], "Sigmoidal-Linear Drift")
    ## the held fraction is not displayed
    expect_false(any(grepl("drift_fraction", out)))

    p <- plot(result, components = TRUE)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    comps <- Filter(\(l) {
        inherits(l$geom, "GeomLine") &&
            identical(l$aes_params$linetype, "dotted")
    }, p$layers)
    expect_length(comps, 2L)

    cf <- result$coefficients
    onset <- sigdrift_onset(
        cf$A, cf$B, cf$xmid, cf$slope, cf$drift_fraction, "symmetric"
    )
    t_rel <- comps[[2L]]$data$time - cf$start_time
    ## drift line starts at the onset, at the fitted drift rate
    expect_true(all(t_rel >= onset))
    expect_equal(
        diff(comps[[2L]]$data$comp2),
        rep(cf$slope_B, nrow(comps[[2L]]$data) - 1L)
    )
    ## sigmoid plus the drift term recovers the fitted curve
    d <- comps[[1L]]$data
    t1 <- d$time - cf$start_time
    expect_equal(
        d$comp1 + cf$slope_B * pmax(t1 - onset, 0),
        result$data[[1L]]$smo2_fitted[is.finite(result$data[[1L]]$smo2_fitted)]
    )
})


## integration tests ================================================

test_that("analyse_sigmoidal_drift() converges on real dataset", {
    skip("Manual fit convergence check")

    intervals <- readRDS("tests/testthat/testdata/5-1_intervals_short.rds")
    deoxy <- intervals[grepl("^deoxy", names(intervals))]
    reoxy <- intervals[grepl("^reoxy", names(intervals))]

    deoxy_channels <- c("smo2_left_vl", "smo2_right_vl")
    reoxy_channels <- c("SmO2 Live")

    results <- analyse_kinetics(
        deoxy,
        # nirs_channels = c(smo2_left_vl, smo2_right_vl),
        method = "sig-lin",
        # end_window = 30
    )
    warnings()
    plot(results)

    ## end-to-end path: window detection and held drift onset.
    ## start_time = 0 anchors the fit at the interval onset
    results <- lapply(intervals, \(df) {
        analyse_exponential_drift(
            df,
            nirs_channels = nirs_channels,
            start_time = 0,
            use_TD = TRUE,
            verbose = FALSE
        )
    })

    coefs <- do.call(rbind, results)
    success <- tapply(!is.na(coefs$tau), coefs$nirs_channels, mean)
    success
    expect_true(all(success >= 1.0))

    TD_success <- tapply(!is.na(coefs$TD), coefs$nirs_channels, mean)
    TD_success
    expect_true(all(TD_success >= 0.8))

    ## converged fits keep the drift onset inside the record
    ok <- !is.na(coefs$tau)
    expect_true(all(coefs$texc[ok] >= 0))

    r2 <- unlist(lapply(results, \(x) attr(x, "diagnostics")$r2))
    mean(r2, na.rm = TRUE)
    expect_true(mean(r2 > 0.9, na.rm = TRUE) >= 0.8)
})
