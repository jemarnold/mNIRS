## biexponential() ==================================================
test_that("biexponential() returns correct vector length", {
    t <- 0:120
    result <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40
    )

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("biexponential() starts at A and approaches the plateau", {
    t <- 0:500
    A <- 70
    B1 <- 40
    B2 <- 60
    result <- biexponential(
        t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40
    )

    expect_equal(result[1], A)
    expect_true(
        all.equal(result[length(result)], B2, tolerance = 0.01,
            scale = 1)
    )
})

test_that("biexponential() slow component is active from the onset", {
    t <- 0:120
    A <- 70
    B1 <- 40
    tau1 <- 5
    B2 <- 60
    tau2 <- 40
    result <- biexponential(t, A = A, B1 = B1, tau1 = tau1, B2 = B2, tau2 = tau2)

    ## the fast and slow terms run concurrently from t = 0, so the curve
    ## differs from the pure fast monoexponential for every t > 0
    fast_only <- monoexponential(t, A = A, B = B1, tau = tau1)
    expect_true(all(result[t > 0] != fast_only[t > 0]))
    expect_equal(
        result,
        fast_only + (B2 - B1) * (1 - exp(-t / tau2))
    )
})

test_that("biexponential() drops to a nadir below A and below the plateau", {
    t <- 0:120
    A <- 70
    B1 <- 40
    B2 <- 60
    result <- biexponential(
        t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40
    )

    ## interior minimum below both endpoints (nadir-recovery shape)
    expect_true(min(result) < A)
    expect_true(min(result) < B2)
    expect_true(which.min(result) > 1 && which.min(result) < length(result))
})

test_that("biexponential() TD form is flat before the delay", {
    t <- 0:120
    TD <- 15
    result <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40,
        TD = TD
    )

    expect_true(all(result[t < TD] == 70))
    expect_equal(result[t == TD], 70)
})

test_that("biexponential() reduces to the monoexponential when B1 = B2", {
    t <- 0:120
    ## the slow term vanishes and the model is the exact monoexponential
    expect_equal(
        biexponential(
            t, A = 70, B1 = 40, tau1 = 5, B2 = 40, tau2 = 50
        ),
        monoexponential(t, A = 70, B = 40, tau = 5)
    )
})


## SSbiexponential() ================================================
test_that("SSbiexponential() converges on known parameters", {
    set.seed(1)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B1", "tau1", "B2", "tau2"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], 70, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B1"]], 40, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["B2"]], 60, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["tau1"]], 5, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["tau2"]], 40, tolerance = 15, scale = 1))
})

test_that("SSbiexponential() fits the 6-parameter TD form", {
    set.seed(4)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40,
        TD = 10
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, TD),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_named(
        coef(model), c("A", "B1", "tau1", "B2", "tau2", "TD")
    )
    expect_true(all.equal(coef(model)[["TD"]], 10, tolerance = 3, scale = 1))
})

test_that("SSbiexponential() predict() returns correct length", {
    set.seed(2)
    t <- 0:120
    x <- biexponential(t, 70, 40, 5, 60, 40) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_length(predict(model, data), nrow(data))
})

test_that("SSbiexponential() fixes A at a constant", {
    set.seed(3)
    t <- 0:120
    x <- biexponential(
        t, A = 0, B1 = -25, tau1 = 5, B2 = -10, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    suppressWarnings(
        model <- nls(
            x ~ SSbiexponential(t, A = 0, B1, tau1, B2, tau2),
            data = data,
            algorithm = "port",
            lower = c(-Inf, 0, -Inf, 0),
            control = nls.control(warnOnly = TRUE)
        )
    )

    expect_named(coef(model), c("B1", "tau1", "B2", "tau2"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 0)
})

test_that("SSbiexponential() gradient matches numericDeriv for the free parameters", {
    t <- seq(-10, 120, by = 0.5)
    env <- list2env(list(
        t = t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40, TD = 3
    ))
    chk <- function(expr, pars) {
        an <- attr(eval(expr, env), "gradient")
        nd <- attr(numericDeriv(expr, pars, env), "gradient")
        expect_identical(colnames(an), pars)
        expect_equal(unname(an), unname(nd), tolerance = 1e-5)
    }
    chk(
        quote(SSbiexponential(t, A, B1, tau1, B2, tau2, TD)),
        c("A", "B1", "tau1", "B2", "tau2", "TD")
    )
    ## a constant in the formula contributes no column
    chk(
        quote(SSbiexponential(t, A, B1, tau1 = 5, B2, tau2)),
        c("A", "B1", "B2", "tau2")
    )
    ## no free parameter, no gradient; the exported fn stays plain
    expect_null(attr(SSbiexponential(t, 70, 40, 5, 60, 40), "gradient"))
    expect_null(attr(biexponential(t, 70, 40, 5, 60, 40), "gradient"))
})

test_that("biexp_start() matches a per-point least-squares grid search", {
    set.seed(8)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    start <- biexp_start(x, t, has_TD = TRUE)
    expect_named(start, c("A", "B1", "tau1", "B2", "tau2", "TD"))

    ## the linear coefficients at the chosen grid point are the lm solution
    ts <- pmax(t - start[["TD"]], 0)
    e1 <- exp(-ts / start[["tau1"]])
    e2 <- exp(-ts / start[["tau2"]])
    cf <- lm.fit(cbind(e1, e2 - e1, 1 - e2), x)$coefficients
    expect_equal(unname(start[c("A", "B1", "B2")]), unname(cf), tolerance = 1e-8)
    expect_true(start[["tau2"]] >= start[["tau1"]] / tau_ratio)
})

test_that("SSbiexponential() fixes tau1 in the formula", {
    set.seed(5)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1 = 5, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_named(coef(model), c("A", "B1", "B2", "tau2"))
    expect_true(all.equal(coef(model)[["B1"]], 40, tolerance = 4, scale = 1))
})


## analyse_biexponential() ==========================================

## helper: create excursion-recovery test data with known parameters
create_biexp_data <- function(
    A = 70,
    B1 = 40,
    tau1 = 5,
    B2 = 60,
    tau2 = 40,
    n = 120,
    sample_rate = 1,
    noise_sd = 0.5,
    channels = "smo2",
    seed = 42
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    x <- biexponential(t, A, B1, tau1, B2, tau2) +
        rnorm(n, 0, noise_sd)

    df <- setNames(data.frame(t, x), c("time", channels[1]))
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- biexponential(
                t, A + 5, B1 + 5, tau1, B2 + 5, tau2
            ) + rnorm(n, 0, noise_sd)
        }
    }

    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate
    )
}


test_that("analyse_biexponential() returns correct structure", {
    data <- create_biexp_data()
    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        verbose = TRUE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "A", "B1", "tau1",
        "MRT", "texc", "B2", "tau2", "TD", "MRT_fitted", "texc_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    expect_type(attr(result, "model"), "list")
    expect_true(inherits(attr(result, "model")$smo2, "nls"))
    expect_s3_class(attr(result, "fitted_data")$smo2, "data.frame")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_s3_class(attr(result, "diagnostics"), "data.frame")
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_s3_class(attr(result, "channel_args"), "data.frame")
})

test_that("analyse_biexponential() recovers known parameters", {
    A <- 70
    B1 <- 40
    B2 <- 60
    data <- create_biexp_data(A = A, B1 = B1, B2 = B2, noise_sd = 0.3)
    
    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 3, scale = 1))
    ## the fast asymptote absorbs some slow phase over the stage-1 window
    expect_true(all.equal(result$B1, B1, tolerance = 5, scale = 1))
    expect_true(all.equal(result$B2, B2, tolerance = 3, scale = 1))
    ## no TD: MRT is the fast time constant
    expect_equal(result$MRT, result$tau1)
    ## excursion sits inside the window, below the starting value
    expect_true(result$texc > 0)
    expect_true(result$texc_fitted < result$A)
    ## good fit
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
})

test_that("analyse_biexponential() texc is the fitted turning point", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    ## texc_fitted is the model prediction at the fitted turning point
    expect_equal(
        result$texc_fitted,
        biexponential(
            result$texc, result$A, result$B1, result$tau1,
            result$B2, result$tau2
        )
    )
    expect_true(result$texc > 0)
    expect_true(result$tau1 < result$tau2)
    ## the turning point sits near the fitted-curve minimum
    fitted <- attr(result, "fitted_data")$smo2$fitted
    expect_true(all.equal(result$texc_fitted, min(fitted), tolerance = 1,
        scale = 1))
})

test_that("analyse_biexponential() reports tau1 <= tau2", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )
    expect_lte(result$tau1, result$tau2)

    set.seed(6)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.3)
    data <- data.frame(t, x)
    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )
    expect_true(coef(model)[["tau1"]] < coef(model)[["tau2"]])
})

test_that("analyse_biexponential() reports NA texc for a monotonic fit", {
    ## B1 between A and B2: the fitted curve is monotonic, so there is no
    ## interior turning point
    set.seed(7)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 55, tau1 = 5, B2 = 40, tau2 = 40) +
        rnorm(length(t), 0, 0.3)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(is.na(result$texc))
    expect_true(is.na(result$texc_fitted))
    expect_false(is.na(result$A))
    expect_false(is.na(result$B1))
    expect_false(is.na(result$tau1))
    expect_false(is.na(result$B2))
    expect_false(is.na(result$tau2))
})

test_that("analyse_biexponential() uses start_time correctly", {
    start_time <- 12
    data <- create_biexp_data(noise_sd = 0.3)
    data$time <- data$time + start_time

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, 70, tolerance = 3, scale = 1))
    ## plateau = B2
    expect_true(all.equal(result$B2, 60, tolerance = 3, scale = 1))
})

test_that("analyse_biexponential() works with multiple channels", {
    nirs_channels <- c("smo2_left", "smo2_right")
    data <- create_biexp_data(channels = nirs_channels)

    result <- analyse_biexponential(
        data,
        nirs_channels = nirs_channels,
        verbose = FALSE
    )

    expect_equal(nrow(result), 2L)
    expect_equal(result$nirs_channels, nirs_channels)
    expect_named(attr(result, "fitted_data"), nirs_channels)
})

test_that("analyse_biexponential() returns NA for failed fit", {
    ## too few observations for the model
    custom_name <- create_biexp_data(n = 4, noise_sd = 0.1)

    expect_warning(
        result <- analyse_biexponential(
            custom_name,
            nirs_channels = "smo2",
            use_TD = FALSE
        ),
        "fit failed for.*smo2.*custom_name"
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$tau1))
    expect_true(is.na(result$MRT))
    expect_true(is.na(result$texc_fitted))
})

test_that("analyse_biexponential() suppresses fit-failure warning when verbose = FALSE", {
    custom_name <- create_biexp_data(n = 4, noise_sd = 0.1)

    expect_no_warning(
        analyse_biexponential(
            custom_name,
            nirs_channels = "smo2",
            verbose = FALSE
        )
    )
})

test_that("analyse_biexponential() end_window bounds the fast phase only", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        end_window = 40,
        use_TD = FALSE,
        verbose = FALSE
    )

    ## stage 2 spans the full response regardless of end_window
    fitted_data <- attr(result, "fitted_data")$smo2
    expect_equal(nrow(fitted_data), nrow(data))
    expect_equal(attr(result, "diagnostics")$n_obs, nrow(data))

    ## the stage-1 window drives tau1
    full <- analyse_biexponential(
        data, nirs_channels = "smo2", use_TD = FALSE, verbose = FALSE
    )
    expect_false(isTRUE(all.equal(result$tau1, full$tau1)))
})

test_that("analyse_biexponential() bounds the fast phase about the stage-1 fit", {
    data <- create_biexp_data(noise_sd = 0.3)

    ## stage 1 is the monoexponential fit on the same window
    fast <- analyse_monoexponential(
        data, nirs_channels = "smo2", end_window = 20, use_TD = FALSE,
        verbose = FALSE
    )
    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        end_window = 20,
        use_TD = FALSE,
        tau1_flex = 0.1,
        A_flex = 0.5,
        verbose = FALSE
    )

    expect_true(result$tau1 >= fast$tau / 1.1 - 1e-8)
    expect_true(result$tau1 <= fast$tau * 1.1 + 1e-8)
    expect_true(abs(result$A - fast$A) <= 0.5 + 1e-8)
    ## the slow phase separates above the fast-phase ceiling
    expect_true(result$tau2 >= fast$tau * 1.1 / tau_ratio - 1e-8)

    ## flex args pass through `...` of analyse_kinetics() and are recorded
    ca <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        end_window = 20,
        tau1_flex = list(smo2 = 0.1),
        TD_flex = 1,
        verbose = FALSE
    )$channel_args
    expect_equal(ca$tau1_flex, 0.1)
    expect_equal(ca$TD_flex, 1)
    expect_true(is.na(ca$A_flex))
    ca <- analyse_kinetics(
        data, nirs_channels = "smo2", method = "biexponential", verbose = FALSE
    )$channel_args
    expect_equal(ca$tau1_flex, 1 / 3)
    expect_equal(ca$TD_flex, 2)
})

test_that("analyse_biexponential() fits a mirrored rise-overshoot response", {
    ## inverted kinetics: the response rises to a peak before settling back.
    ## the asymptote ordering flips (B1 above A, B2 below B1) and is
    ## recovered from the data
    set.seed(31)
    t <- 0:119
    x <- biexponential(
        t, A = 30, B1 = 55, tau1 = 10, B2 = 40, tau2 = 50
    ) + rnorm(120, 0, 0.3)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$B1, 55, tolerance = 5, scale = 1))
    expect_true(all.equal(result$B2, 40, tolerance = 5, scale = 1))
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
    ## a mirrored response reports a genuine interior excursion: a maximum
    ## above the baseline
    expect_true(result$texc > 0)
    expect_true(result$texc_fitted > result$A)
})


## time delay =======================================================

test_that("analyse_biexponential() use_TD = FALSE forces the 5-param fit", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_named(coef(attr(result, "model")$smo2),
        c("A", "B1", "tau1", "B2", "tau2"))
    expect_true(is.na(result$TD))
})

test_that("analyse_biexponential() falls back to the 5-parameter fit", {
    ## the stage-1 TD fit fails on seven observations and retries without
    ## TD, so stage 2 fits the 5-parameter model
    data <- create_biexp_data(n = 7, noise_sd = 0.1)

    warns <- capture_warnings(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = TRUE
        )
    )

    expect_match(warns[[1L]], "4-parameter `SSmonoexponential")
    expect_match(warns[[1L]], "Attempting")
    ## TD is absent from the reduced model whether or not the retry converges
    expect_true(is.na(result$TD))
})

test_that("analyse_biexponential() texc is elapsed from start_time", {
    ## 6-parameter TD fit with a non-zero start_time: texc must be
    ## measured from start_time, not from the model's internal (TD) onset
    set.seed(11)
    start_time <- 20
    TD <- 8
    t <- start_time + 0:119
    x <- biexponential(
        t - start_time, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40,
        TD = TD
    ) + rnorm(120, 0, 0.15)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        use_TD = TRUE,
        verbose = FALSE
    )

    ## must fit the 6-param model, not fall back to the 5-param form
    expect_named(coef(attr(result, "model")$smo2),
        c("A", "B1", "tau1", "B2", "tau2", "TD"))

    ## MRT is the fast-phase mean response time from the fit onset
    expect_equal(result$MRT, result$TD + result$tau1)
    ## texc includes the TD offset from the fit onset
    expect_true(result$texc > TD)
})


## direction ========================================================

test_that("analyse_biexponential() direction steers the fit window", {
    ## an excursion in the requested direction fits without a refit
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        direction = "negative",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_false(is.na(result$A))
    expect_equal(attr(result, "channel_args")$direction, "negative")
})

test_that("analyse_biexponential() direction = 'negative' matches auto on falling data", {
    data <- create_biexp_data(noise_sd = 0.3)

    result_auto <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        direction = "auto",
        verbose = FALSE
    )
    result_neg <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        direction = "negative",
        verbose = FALSE
    )

    ## matching direction leaves the unconstrained fit untouched
    expect_equal(result_auto$A, result_neg$A)
    expect_equal(result_auto$B2, result_neg$B2)
    expect_equal(result_auto$tau1, result_neg$tau1)
    expect_true(result_auto$B2 < result_auto$A)
})

test_that("analyse_biexponential() direction = 'positive' rejects a negative response", {
    ## genuinely falling drop-recovery: no positive response exists
    ## within the data span, so the bounded refit degenerates
    data <- create_biexp_data(noise_sd = 0.3)

    expect_warning(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = FALSE,
            direction = "positive",
            verbose = TRUE
        ),
        "satisfy"
    )

    ## never returns a within-span negative fit against requested direction
    expect_true(is.na(result$A))
    expect_true(is.na(result$B2))
    expect_true(is.na(result$tau1))
})

test_that("analyse_biexponential() suppresses direction warning when verbose = FALSE", {
    data <- create_biexp_data(noise_sd = 0.3)

    expect_no_warning(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = FALSE,
            direction = "positive",
            verbose = FALSE
        )
    )
    expect_true(is.na(result$A))
})

test_that("analyse_biexponential() direction with both asymptotes fixed returns NA", {
    ## falling data with both asymptotes fixed contradicts the requested
    ## positive direction: no refit possible
    data <- create_biexp_data(noise_sd = 0.3)

    expect_warning(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = FALSE,
            fix = list(A = 70, B2 = 60),
            direction = "positive",
            verbose = TRUE
        ),
        "satisfy"
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$tau1))
})

test_that("analyse_kinetics() passes direction to biexponential method", {
    data <- create_biexp_data(noise_sd = 0.3)

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "biexponential",
            use_TD = FALSE,
            direction = "positive"
        ),
        "satisfy"
    )

    expect_true(is.na(result$coefficients$A))
})


## fixed parameters =================================================

test_that("analyse_biexponential() fix holds parameters constant", {
    data <- create_biexp_data(A = 0, B1 = -20, B2 = -10, noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(A = 0),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$A, 0)

    ## fixed A excluded from the fitted model coefficients
    expect_named(
        coef(attr(result, "model")$smo2),
        c("B1", "tau1", "B2", "tau2")
    )
    expect_equal(attr(result, "channel_args")$fix, "list(A = 0)")
    expect_false(is.na(attr(result, "diagnostics")$adj_r2))
})

test_that("analyse_biexponential() fix holds tau1 constant", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(tau1 = 12),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$tau1, 12)
    expect_named(
        coef(attr(result, "model")$smo2), c("A", "B1", "B2", "tau2")
    )
})

test_that("analyse_biexponential() fix holds tau2 constant on the sequential fit", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(tau2 = 40),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$tau2, 40)
    expect_named(
        coef(attr(result, "model")$smo2), c("A", "B1", "tau1", "B2")
    )
    expect_true(all.equal(result$tau1, 5, tolerance = 2, scale = 1))
})

test_that("analyse_biexponential() returned model is canonical and converged", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        verbose = FALSE
    )
    model <- attr(result, "model")$smo2

    expect_named(coef(model), c("A", "B1", "tau1", "B2", "tau2", "TD"))
    expect_true(model$convInfo$isConv)
    expect_true(result$tau2 >= result$tau1 / tau_ratio)
    ## coefficient table mirrors the model
    expect_equal(result$tau1, coef(model)[["tau1"]])
    ## the self-start model predicts with a gradient attribute
    pred <- predict(model, newdata = data.frame(time = c(0, 10)))
    expect_length(pred, 2L)
    expect_true(all(is.finite(pred)))
})

test_that("analyse_biexponential() caps a runaway tau2 at 10x the span", {
    ## fast dip onto a linear ramp: the slow limb has no finite time
    ## constant, so tau2 runs to its upper bound
    set.seed(42)
    t <- 0:59
    x <- 70 - 25 * exp(-t / 5) + 0.15 * t + rnorm(60, 0, 0.3)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- suppressWarnings(analyse_biexponential(
        data, nirs_channels = "smo2", use_TD = FALSE, verbose = FALSE
    ))

    ## pinned at the port upper bound rather than left to diverge
    expect_true(result$tau2 <= 10 * diff(range(t)) + 1e-6)
})

test_that("analyse_biexponential() fix holds TD constant", {
    set.seed(21)
    TD <- 8
    t <- 0:119
    x <- biexponential(
        t, A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40,
        TD = TD
    ) + rnorm(120, 0, 0.2)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = TRUE,
        fix = list(TD = TD),
        verbose = FALSE
    )

    expect_equal(result$TD, TD)
    ## a fixed TD is excluded from estimation and disables the 5-param retry
    expect_named(
        coef(attr(result, "model")$smo2),
        c("A", "B1", "tau1", "B2", "tau2")
    )
})

test_that("analyse_biexponential() TD is only fixable when use_TD = TRUE", {
    data <- create_biexp_data()

    expect_error(
        analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = FALSE,
            fix = list(TD = 0)
        ),
        "not recognised"
    )
})

test_that("analyse_biexponential() fix resolves per channel", {
    ## ch2 asymptotes are each + 5 by construction
    data <- create_biexp_data(
        A = 0, B1 = -20, B2 = -10, noise_sd = 0.3,
        channels = c("ch1", "ch2")
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = c("ch1", "ch2"),
        use_TD = FALSE,
        fix = list(ch1 = list(A = 0), ch2 = list(A = 5)),
        verbose = FALSE
    )

    expect_equal(result$A, c(0, 5))

    models <- attr(result, "model")
    expect_named(coef(models$ch1), c("B1", "tau1", "B2", "tau2"))
    expect_named(coef(models$ch2), c("B1", "tau1", "B2", "tau2"))

    ca <- attr(result, "channel_args")
    expect_equal(ca$fix, c("list(A = 0)", "list(A = 5)"))
})

test_that("analyse_biexponential() validates fix argument", {
    data <- create_biexp_data()

    ## unnamed list
    expect_error(
        analyse_biexponential(data, nirs_channels = "smo2", fix = list(0)),
        "uniquely named"
    )
    ## unknown parameter name
    expect_error(
        analyse_biexponential(data, nirs_channels = "smo2", fix = list(Q = 1)),
        "not recognised"
    )
    ## cannot fix every parameter
    expect_error(
        analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = TRUE,
            fix = list(
                A = 70, B1 = 40, tau1 = 5, B2 = 60, tau2 = 40, TD = 0
            )
        ),
        "Nothing to estimate"
    )
})


## dispatch =========================================================

test_that("analyse_kinetics() dispatches to the biexponential method", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs_kinetics")
    expect_equal(result$method, "biexponential")
    expect_true(all(
        c("B1", "tau1", "B2", "tau2", "texc") %in%
            names(result$coefficients)
    ))
})

test_that("analyse_kinetics() resolves the 'biexp' alias", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexp",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$method, "biexponential")
})

test_that("analyse_kinetics() passes fix to the biexponential method", {
    data <- create_biexp_data(A = 0, B1 = -20, B2 = -10, noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        fix = list(A = 0),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$coefficients$A, 0)
    expect_named(
        coef(result$model[[1L]]$smo2),
        c("B1", "tau1", "B2", "tau2")
    )
})


## monoexponential reduction ========================================

## monotonic two-phase response: B1 between A and B2, no turning point
create_monotonic_data <- function(seed = 7, TD = NULL, t = 0:119) {
    set.seed(seed)
    x <- biexponential(t, A = 70, B1 = 55, tau1 = 5, B2 = 40, tau2 = 40, TD) +
        rnorm(length(t), 0, 0.3)
    create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )
}

test_that("analyse_kinetics() reduces a monotonic biexponential fit", {
    data <- create_monotonic_data()

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "biexponential",
            use_TD = FALSE
        ),
        "reduced to"
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_equal(names(cf)[1:4], c("interval", "nirs_channels", "start_time", "model"))
    expect_equal(cf$model, "monoexponential")
    expect_named(coef(model), c("A", "B", "tau"))
    expect_equal(cf$A, coef(model)[["A"]])
    expect_equal(cf$B1, coef(model)[["B"]])
    expect_equal(cf$B2, coef(model)[["B"]])
    expect_equal(cf$tau1, coef(model)[["tau"]])
    ## reduced MRT carries over from the monoexponential (no TD: tau)
    expect_equal(cf$MRT, cf$tau1)
    expect_true(is.na(cf$tau2))
    expect_true(is.na(cf$texc))
    expect_true(is.na(cf$texc_fitted))
    ## diagnostics and fitted values follow the reduced model
    expect_equal(result$diagnostics$n_params, 3L)
    expect_equal(
        result$data[[1L]]$smo2_fitted,
        as.vector(predict(model))
    )
    ## recorded with the reason, regardless of verbose
    msgs <- result$warnings$message
    expect_true(any(grepl("reduced to", msgs)))
    expect_true(any(grepl("monotonic", msgs)))
})

test_that("analyse_kinetics() reduces a monoexponential response", {
    set.seed(3)
    t <- 0:120
    x <- monoexponential(t, A = 70, B = 40, tau = 8) + rnorm(length(t), 0, 0.3)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$coefficients$model, "monoexponential")
    expect_true(all.equal(result$coefficients$tau1, 8, tolerance = 1, scale = 1))
    expect_true(inherits(result$model[[1L]]$smo2, "nls"))
})

test_that("analyse_kinetics() keeps a supported excursion-recovery fit", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        verbose = FALSE
    )
    forced <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        force_biexponential = TRUE,
        verbose = FALSE
    )

    expect_equal(result$coefficients$model, "biexponential")
    expect_true(is.finite(result$coefficients$texc))
    expect_named(
        coef(result$model[[1L]]$smo2),
        c("A", "B1", "tau1", "B2", "tau2", "TD")
    )
    expect_false(any(grepl("reduced to", result$warnings$message)))
    ## the raw fit is unchanged by the comparison
    expect_equal(result$coefficients, forced$coefficients)
    expect_equal(result$diagnostics, forced$diagnostics)
})

test_that("force_biexponential = TRUE keeps a monotonic fit", {
    data <- create_monotonic_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        force_biexponential = TRUE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "biexponential")
    expect_true(is.na(cf$texc))
    expect_false(is.na(cf$tau2))
    expect_false(cf$B1 == cf$B2)
    expect_false(any(grepl("reduced to", result$warnings$message)))
})

test_that("reduced fit matches the time-delay structure and window", {
    data <- create_monotonic_data(seed = 11, TD = 10, t = -20:120)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        verbose = FALSE
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_equal(cf$model, "monoexponential")
    expect_equal(is.finite(cf$TD), "TD" %in% names(coef(model)))
    ## fitted rows are the reduced model's window
    expect_equal(
        result$diagnostics$n_obs,
        sum(is.finite(result$data[[1L]]$smo2_fitted))
    )
    expect_equal(result$diagnostics$n_params, length(coef(model)))
})

test_that("reduced comparator fits the full response, not end_window", {
    data <- create_monotonic_data(seed = 11, TD = 10, t = -20:120)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        end_window = 30,
        verbose = FALSE
    )

    expect_equal(result$coefficients$model, "monoexponential")
    expect_gte(result$diagnostics$n_obs, sum(data$time >= 0))
    expect_equal(result$channel_args$end_window, 30)
})

test_that("reduction carries fixed parameters over", {
    data <- create_monotonic_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        fix = list(A = 70, B2 = 40),
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "monoexponential")
    expect_equal(cf$A, 70)
    expect_equal(cf$B1, 40)
    expect_equal(cf$B2, 40)
    expect_named(coef(result$model[[1L]]$smo2), "tau")

    ## parameters without a monoexponential counterpart are dropped
    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        fix = list(tau1 = 5),
        verbose = FALSE
    )

    expect_equal(result$coefficients$model, "monoexponential")
    expect_named(coef(result$model[[1L]]$smo2), c("A", "B", "tau"))
})

test_that("reduction resolves per interval", {
    data <- list(
        excursion = create_biexp_data(noise_sd = 0.3),
        monotonic = create_monotonic_data()
    )

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$interval, c("excursion", "monotonic"))
    expect_equal(cf$model, c("biexponential", "monoexponential"))
    expect_named(
        coef(result$model$excursion$smo2), c("A", "B1", "tau1", "B2", "tau2")
    )
    expect_named(coef(result$model$monotonic$smo2), c("A", "B", "tau"))
    expect_true(is.finite(cf$texc[[1L]]))
    expect_true(is.na(cf$texc[[2L]]))
    reduced <- result$warnings[grepl("reduced to", result$warnings$message), ]
    expect_equal(reduced$interval, "monotonic")
    expect_equal(reduced$nirs_channels, "smo2")
})

test_that("reduction resolves per channel", {
    data <- create_biexp_data(noise_sd = 0.3, channels = c("smo2", "hhb"))
    data$hhb <- create_monotonic_data()$smo2

    result <- analyse_kinetics(
        data,
        nirs_channels = c(smo2, hhb),
        method = "biexponential",
        use_TD = list(smo2 = TRUE, hhb = FALSE),
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$nirs_channels, c("smo2", "hhb"))
    expect_equal(cf$model, c("biexponential", "monoexponential"))
    expect_true(is.finite(cf$TD[[1L]]))
    expect_true(is.na(cf$TD[[2L]]))
    expect_true(inherits(result$model[[1L]]$smo2, "nls"))
    expect_true(inherits(result$model[[1L]]$hhb, "nls"))
    expect_equal(
        sum(is.finite(result$data[[1L]]$hhb_fitted)),
        result$diagnostics$n_obs[[2L]]
    )
})

test_that("a row where both fits fail is left as is", {
    data <- create_mnirs_data(
        data.frame(time = 0:2, smo2 = c(70, 60, 55)),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "biexponential")
    expect_true(is.na(cf$A))
    expect_null(result$model[[1L]]$smo2)
    expect_false(any(grepl("reduced to", result$warnings$message)))
})

test_that("reduce_kinetics() applies the shape test and F-test", {
    data <- create_monotonic_data()
    args <- list(
        start_time = NULL, direction = "auto", end_window = Inf,
        use_TD = FALSE, fix = NULL
    )
    run <- \(worker, method) {
        analyse_kinetics_intervals(
            data, worker, method, args,
            rlang::quo(smo2), rlang::quo(NULL),
            "ensemble", FALSE, FALSE, quote(f()), quote(f())
        )
    }
    full <- run(analyse_biexponential, "biexponential")
    reduced <- run(analyse_monoexponential, "monoexponential")
    spec <- kinetics_reductions$biexponential

    ## the shape test rejects the monotonic fit
    out <- reduce_kinetics(full, reduced, spec, verbose = FALSE)
    expect_equal(out$coefficients$model, "monoexponential")
    expect_match(out$warnings$message, "monotonic", all = FALSE)

    ## without it the F-test decides on its own
    spec$accept <- NULL
    out <- reduce_kinetics(full, reduced, spec, verbose = FALSE)
    expect_true(out$coefficients$model %in% c("biexponential", "monoexponential"))
    if (out$coefficients$model == "monoexponential") {
        expect_match(out$warnings$message, "F\\(", all = FALSE)
    }
})

test_that("map_fix() renames through nested maps", {
    fix_map <- kinetics_reductions$biexponential$fix_map

    expect_equal(map_fix(list(A = 1, B2 = 2, tau1 = 3), fix_map), list(A = 1, B = 2))
    expect_equal(map_fix(list(tau1 = 3), fix_map), setNames(list(), character()))
    expect_equal(
        map_fix(list(smo2 = list(B2 = 2), hhb = list(tau2 = 4)), fix_map),
        list(smo2 = list(B = 2), hhb = setNames(list(), character()))
    )
})


## integration tests ================================================

test_that("analyse_biexponential() converges on real dataset", {
    skip("Manual fit convergence check")

    intervals <- readRDS(test_path("testdata/5-1_intervals_short.rds"))
    deoxy <- intervals[grepl("^deoxy", names(intervals))]

    results <- analyse_kinetics(
        deoxy,
        # nirs_channels = c(smo2_left_vl, smo2_right_vl),
        method = "biexponential",
        end_window = 30
    )
    tibble(results$warnings)
    # plot(results)
    # plot(results, components = TRUE, scales = "free")

    coefs <- results$coefficients
    ## deoxy_2 is an excursion-recovery, deoxy_3 a monotonic drop
    expect_equal(
        coefs$model,
        c("biexponential", "biexponential", "monoexponential", "monoexponential")
    )
    expect_true(all(!is.na(coefs$tau1)))

    ## reoxy
    reoxy <- intervals[grepl("^reoxy", names(intervals))]

    results <- analyse_kinetics(
        reoxy,
        nirs_channels = "SmO2 Live",
        method = "biexponential",
        verbose = FALSE
    )
    # plot(results)
    # plot(results, components = TRUE, scales = "free")

    coefs <- results$coefficients
    (success <- mean(!is.na(coefs$tau1)))
    expect_true(all(success >= 1.0))

    (texc_success <- mean(!is.na(coefs$texc)))
    expect_true(all(success >= 0.75))
})
