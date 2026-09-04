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
            sigmoidal_drift(t, 10, 100, 40, 4, 0, 0, 0.95, shape = .s),
            sigmoid_fn[[.s]](t, 10, 100, 40, 4)
        )
    })
})

test_that("sigdrift_cutoffs() are where the sigmoid reaches the fraction", {
    p <- 0.95
    check <- \(A, B, slope) {
        lapply(shapes, \(.s) {
            cut <- sigdrift_cutoffs(A, B, 40, slope, p, .s)
            expect_named(cut, c("texc_A", "texc_B"))
            expect_lt(cut[[1L]], 40)
            expect_gt(cut[[2L]], 40)
            expect_equal(
                sigmoid_fn[[.s]](cut, A, B, 40, slope),
                A + c(1 - p, p) * (B - A),
                ignore_attr = TRUE
            )
        })
    }
    check(10, 100, 4)
    ## falling response
    check(100, 10, -4)

    ## the Gompertz forms reach further on their slow side
    g <- sigdrift_cutoffs(10, 100, 40, 4, p, "gompertz")
    expect_lt(40 - g[[1L]], g[[2L]] - 40)
    gl <- sigdrift_cutoffs(10, 100, 40, 4, p, "gompertz_left")
    expect_gt(40 - gl[[1L]], gl[[2L]] - 40)
    s <- sigdrift_cutoffs(10, 100, 40, 4, p, "symmetric")
    expect_equal(40 - s[[1L]], s[[2L]] - 40)
})

test_that("sigmoidal_drift() drifts are hinged at the cutoffs and independent", {
    t <- seq(0, 100, by = 0.5)
    lapply(shapes, \(.s) {
        cut <- sigdrift_cutoffs(10, 100, 40, 4, 0.95, .s)
        base <- sigmoidal_drift(t, 10, 100, 40, 4, 0, 0, 0.95, shape = .s)
        lead <- sigmoidal_drift(t, 10, 100, 40, 4, 0.3, 0, 0.95, shape = .s)
        trail <- sigmoidal_drift(t, 10, 100, 40, 4, 0, -0.4, 0.95, shape = .s)
        both <- sigmoidal_drift(t, 10, 100, 40, 4, 0.3, -0.4, 0.95, shape = .s)

        ## each drift is zero outside its region, linear inside
        expect_equal(lead[t >= cut[[1L]]], base[t >= cut[[1L]]])
        expect_equal(
            lead[t < cut[[1L]]] - base[t < cut[[1L]]],
            0.3 * (t[t < cut[[1L]]] - cut[[1L]])
        )
        expect_equal(trail[t <= cut[[2L]]], base[t <= cut[[2L]]])
        expect_equal(
            trail[t > cut[[2L]]] - base[t > cut[[2L]]],
            -0.4 * (t[t > cut[[2L]]] - cut[[2L]])
        )
        ## the regions never overlap, so the drifts add independently
        expect_equal(both, lead + trail - base)
    })
})


## SSsigmoidal_drift() ==============================================

test_that("SSsigmoidal_drift() recovers parameters for every shape", {
    t <- seq(0, 119)
    lapply(shapes, \(.s) {
        set.seed(13)
        x <- sigmoidal_drift(
            t, 10, 100, 40, 4, 0.3, -0.4, 0.95, shape = .s
        ) + rnorm(length(t), 0, 1)
        data <- data.frame(t, x)
        formula <- substitute(
            x ~ SSsigmoidal_drift(
                t, A, B, xmid, slope, slope_A, slope_B,
                drift_frac = 0.95, shape = .s
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
        expect_named(cf, c("A", "B", "xmid", "slope", "slope_A", "slope_B"))
        expect_true(all.equal(cf[["A"]], 10, tolerance = 3, scale = 1))
        expect_true(all.equal(cf[["B"]], 100, tolerance = 3, scale = 1))
        expect_true(all.equal(cf[["xmid"]], 40, tolerance = 2, scale = 1))
        expect_true(all.equal(cf[["slope"]], 4, tolerance = 1, scale = 1))
        expect_true(all.equal(cf[["slope_A"]], 0.3, tolerance = 0.1, scale = 1))
        expect_true(all.equal(cf[["slope_B"]], -0.4, tolerance = 0.1, scale = 1))
    })
})

test_that("SSsigmoidal_drift() excludes fixed parameters from estimation", {
    set.seed(13)
    t <- seq(0, 119)
    x <- sigmoidal_drift(t, 10, 100, 40, 4, 0.3, -0.4, 0.95) +
        rnorm(length(t), 0, 1)
    data <- data.frame(t, x)

    model <- suppressWarnings(nls(
        x ~ SSsigmoidal_drift(
            t, A = 10, B, xmid, slope, slope_A, slope_B, drift_frac = 0.95
        ),
        data = data,
        algorithm = "port",
        control = nls.control(warnOnly = TRUE)
    ))
    expect_named(coef(model), c("B", "xmid", "slope", "slope_A", "slope_B"))
    expect_true(all.equal(coef(model)[["B"]], 100, tolerance = 3, scale = 1))
    ## predicts on the fixed value and the default symmetric shape
    expect_equal(
        as.numeric(predict(model, data.frame(t = 0:5))),
        sigmoidal_drift(0:5, 10, coef(model)[["B"]], coef(model)[["xmid"]],
            coef(model)[["slope"]], coef(model)[["slope_A"]],
            coef(model)[["slope_B"]], 0.95)
    )
})


## sigdrift_start() =================================================

test_that("sigdrift_start() seeds near the truth and holds fixed values", {
    set.seed(13)
    t <- seq(0, 119)
    x <- sigmoidal_drift(t, 10, 100, 40, 4, 0.3, -0.4, 0.95) +
        rnorm(length(t), 0, 1)

    start <- sigdrift_start(x, t)
    expect_named(
        start,
        c("A", "B", "xmid", "slope", "slope_A", "slope_B", "drift_frac")
    )
    expect_true(all.equal(start[["A"]], 10, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["B"]], 100, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["xmid"]], 40, tolerance = 5, scale = 1))
    expect_true(all.equal(start[["slope_A"]], 0.3, tolerance = 0.15, scale = 1))
    expect_true(all.equal(start[["slope_B"]], -0.4, tolerance = 0.15, scale = 1))
    expect_equal(start[["drift_frac"]], 0.95)

    fixed <- sigdrift_start(x, t, fixed = list(A = 12, slope_B = -0.5))
    expect_equal(fixed[["A"]], 12)
    expect_equal(fixed[["slope_B"]], -0.5)
})

test_that("sigdrift_start() seeds zero drift without tail support", {
    ## five points inside the seeded cutoffs, so neither tail has two
    t <- seq(30, 50, by = 5)
    x <- logistic(t, 10, 100, 40, 4)
    start <- sigdrift_start(x, t)
    expect_equal(start[["slope_A"]], 0)
    expect_equal(start[["slope_B"]], 0)
})


## analyse_sigmoidal_drift() ========================================

## helper: rising sigmoid with a positive leading and negative trailing
## drift; symmetric cutoffs at about 23 and 57
create_sigdrift_data <- function(
    A = 10,
    B = 100,
    xmid = 40,
    slope = 4,
    slope_A = 0.3,
    slope_B = -0.4,
    drift_frac = 0.95,
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
            t, A + 5 * .i, B + 5 * .i, xmid, slope, slope_A, slope_B,
            drift_frac, shape
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
        "interval", "nirs_channels", "A", "B", "xmid", "slope", "slope_A",
        "slope_B", "drift_frac", "texc_A", "texc_B", "xmid_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    model <- attr(result, "model")$smo2
    expect_s3_class(model, "nls")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_equal(attr(result, "channel_args")$drift_frac, 0.95)
    expect_equal(attr(result, "channel_args")$shape, "symmetric")

    ## the cutoff fraction is held, never estimated
    expect_named(
        coef(model),
        c("A", "B", "xmid", "slope", "slope_A", "slope_B")
    )
    expect_true(all.equal(result$A, 10, tolerance = 3, scale = 1))
    expect_true(all.equal(result$B, 100, tolerance = 3, scale = 1))
    expect_true(all.equal(result$xmid, 40, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, 4, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope_A, 0.3, tolerance = 0.1, scale = 1))
    expect_true(all.equal(result$slope_B, -0.4, tolerance = 0.1, scale = 1))
    expect_equal(result$drift_frac, 0.95)
    expect_true(attr(result, "diagnostics")$r2 > 0.99)
    expect_equal(attr(result, "diagnostics")$n_params, 6L)

    ## derived columns follow the fitted coefficients
    cut <- sigdrift_cutoffs(
        result$A, result$B, result$xmid, result$slope, 0.95, "symmetric"
    )
    expect_equal(result$texc_A, cut[[1L]])
    expect_equal(result$texc_B, cut[[2L]])
    expect_lt(result$texc_A, result$xmid)
    expect_gt(result$texc_B, result$xmid)
    expect_equal(
        result$xmid_fitted,
        as.numeric(predict(model, data.frame(time = result$xmid)))
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
        expect_true(all.equal(result$slope_A, 0.3, tolerance = 0.1, scale = 1))
        expect_true(all.equal(result$slope_B, -0.4, tolerance = 0.1, scale = 1))
        ## the stored model predicts with its shape
        model <- attr(result, "model")$smo2
        expect_equal(
            as.numeric(predict(model, data.frame(time = c(0, 40, 119)))),
            sigmoidal_drift(
                c(0, 40, 119), result$A, result$B, result$xmid, result$slope,
                result$slope_A, result$slope_B, 0.95, shape = .s
            )
        )
    })
})

test_that("analyse_sigmoidal_drift() drift_frac resolves per channel", {
    data <- create_sigdrift_data(channels = c("smo2", "hhb"))

    result <- analyse_sigmoidal_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_frac = list(smo2 = 0.98, hhb = 0.9),
        verbose = FALSE
    )
    expect_equal(result$drift_frac, c(0.98, 0.9))
    expect_equal(attr(result, "channel_args")$drift_frac, c(0.98, 0.9))
    ## a larger fraction pushes the cutoffs further from xmid
    expect_lt(result$texc_A[[1L]] - result$xmid[[1L]],
        result$texc_A[[2L]] - result$xmid[[2L]])
    expect_gt(result$texc_B[[1L]] - result$xmid[[1L]],
        result$texc_B[[2L]] - result$xmid[[2L]])
    expect_false("drift_frac" %in% names(coef(attr(result, "model")$smo2)))

    ## an omitted channel takes the formal default
    result_part <- analyse_sigmoidal_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_frac = list(smo2 = 0.98),
        verbose = FALSE
    )
    expect_equal(result_part$drift_frac, c(0.98, 0.95))
})

test_that("analyse_sigmoidal_drift() validates drift_frac", {
    data <- create_sigdrift_data()
    lapply(list(0.5, 1, 0, 1.2, "0.95", c(0.9, 0.95)), \(.f) {
        expect_error(
            analyse_sigmoidal_drift(
                data, nirs_channels = "smo2", drift_frac = .f, verbose = FALSE
            ),
            "drift_frac.*must be a valid one-element"
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
    expect_named(coef(attr(result, "model")$smo2), c("B", "xmid", "slope", "slope_A"))
    expect_equal(attr(result, "diagnostics")$n_params, 4L)
    expect_true(all.equal(result$B, 100, tolerance = 3, scale = 1))

    ## the cutoff fraction is not a fixable parameter
    expect_error(
        analyse_sigmoidal_drift(
            data, nirs_channels = "smo2", fix = list(drift_frac = 0.9)
        ),
        "not recognised"
    )
})

test_that("analyse_sigmoidal_drift() fails on too few observations", {
    custom_name <- create_sigdrift_data(n = 5, noise_sd = 0.1)
    expect_warning(
        result <- analyse_sigmoidal_drift(custom_name, "smo2"),
        "fit failed for.*smo2.*custom_name.*5 observations for 6 free"
    )
    expect_true(all(is.na(result[c("A", "slope_A", "texc_B", "xmid_fitted")])))
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
    expect_lt(result_auto$texc_A, result_auto$xmid)
    expect_gt(result_auto$texc_B, result_auto$xmid)
})

test_that("analyse_sigmoidal_drift() direction = 'positive' rejects falling fit", {
    data <- create_sigdrift_data(A = 100, B = 10, slope = -4)

    expect_warning(
        result <- analyse_sigmoidal_drift(
            data, nirs_channels = "smo2", direction = "positive"
        ),
        "satisfy"
    )
    expect_true(all(is.na(result[c("A", "B", "slope", "slope_A", "texc_A")])))
})


## model fallback ===================================================

test_that("analyse_kinetics() keeps supported drifts", {
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
    expect_true(all.equal(cf$slope_A, 0.3, tolerance = 0.1, scale = 1))
    expect_false(any(grepl("fell back to", result$warnings$message)))
    expect_equal(cf, forced$coefficients)
})

test_that("analyse_kinetics() keeps the model when one drift is supported", {
    result <- analyse_kinetics(
        create_sigdrift_data(slope_A = 0),
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        verbose = FALSE
    )
    cf <- result$coefficients
    expect_equal(cf$model, "sigmoidal_drift")
    expect_true(all.equal(cf$slope_A, 0, tolerance = 0.1, scale = 1))
    expect_true(all.equal(cf$slope_B, -0.4, tolerance = 0.1, scale = 1))
    expect_false(any(grepl("fell back to", result$warnings$message)))
})

test_that("analyse_kinetics() falls back from negligible drifts", {
    data <- create_sigdrift_data(slope_A = 0, slope_B = 0, shape = "gompertz")

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
        cf[c("slope_A", "slope_B", "drift_frac", "texc_A", "texc_B")]
    )))
    expect_false(is.na(cf$xmid_fitted))
    expect_equal(result$diagnostics$n_params, 4L)
    expect_equal(
        result$data[[1L]]$smo2_fitted,
        as.vector(predict(model))
    )
    ## the reduced fit keeps the shape
    expect_equal(result$channel_args$shape, "gompertz")
    expect_true(any(grepl("Drift amplitudes", result$warnings$message)))

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
    expect_false(is.na(forced$coefficients$slope_A))
})

test_that("sigmoidal_drift fallback resolves per channel with fix carried", {
    data <- create_sigdrift_data(channels = c("smo2", "hhb"))
    data$hhb <- create_sigdrift_data(slope_A = 0, slope_B = 0, seed = 1)$smo2

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

test_that("analyse_kinetics.sigmoidal_drift passes shape, drift_frac, and fix", {
    result <- analyse_kinetics(
        create_sigdrift_data(shape = "gompertz_left"),
        nirs_channels = "smo2",
        method = "sigmoidal_drift",
        shape = "gompertz_left",
        drift_frac = 0.9,
        fix = list(B = 100),
        verbose = FALSE
    )
    cf <- result$coefficients
    expect_equal(cf$model, "sigmoidal_drift")
    expect_equal(cf$B, 100)
    expect_equal(cf$drift_frac, 0.9)
    expect_equal(result$channel_args$shape, "gompertz_left")
    expect_equal(result$channel_args$drift_frac, 0.9)
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
    expect_false(any(grepl("drift_frac", out)))

    p <- plot(result, components = TRUE)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    comps <- Filter(\(l) {
        inherits(l$geom, "GeomLine") &&
            identical(l$aes_params$linetype, "dotted")
    }, p$layers)
    expect_length(comps, 3L)

    cf <- result$coefficients
    t_rel <- comps[[2L]]$data$time - cf$start_time
    ## leading drift line ends at the cutoff, at the fitted drift rate
    expect_true(all(t_rel <= cf$texc_A))
    expect_equal(
        diff(comps[[2L]]$data$comp2),
        rep(cf$slope_A, nrow(comps[[2L]]$data) - 1L)
    )
    ## sigmoid plus both drift terms recovers the fitted curve
    d <- comps[[1L]]$data
    t1 <- d$time - cf$start_time
    expect_equal(
        d$comp1 + cf$slope_A * pmin(t1 - cf$texc_A, 0) +
            cf$slope_B * pmax(t1 - cf$texc_B, 0),
        result$data[[1L]]$smo2_fitted[is.finite(result$data[[1L]]$smo2_fitted)]
    )
})
