## biexponential() ==================================================
test_that("biexponential() returns correct vector length", {
    t <- 0:120
    result <- biexponential(
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40
    )

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("biexponential() starts at A and approaches the plateau", {
    t <- 0:500
    A <- 70
    B <- 40
    B2 <- 60
    result <- biexponential(
        t, A = A, B = B, tau = 5, B2 = B2, tau2 = 40
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
    B <- 40
    tau <- 5
    B2 <- 60
    tau2 <- 40
    result <- biexponential(t, A = A, B = B, tau = tau, B2 = B2, tau2 = tau2)

    ## the fast and slow terms run concurrently from t = 0, so the curve
    ## differs from the pure fast monoexponential for every t > 0
    fast_only <- monoexponential(t, A = A, B = B, tau = tau)
    expect_true(all(result[t > 0] != fast_only[t > 0]))
    expect_equal(
        result,
        fast_only + (B2 - B) * (1 - exp(-t / tau2))
    )
})

test_that("biexponential() drops to a nadir below A and below the plateau", {
    t <- 0:120
    A <- 70
    B <- 40
    B2 <- 60
    result <- biexponential(
        t, A = A, B = B, tau = 5, B2 = B2, tau2 = 40
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
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40,
        TD = TD
    )

    expect_true(all(result[t < TD] == 70))
    expect_equal(result[t == TD], 70)
})

test_that("biexponential() reduces to the monoexponential when B = B2", {
    t <- 0:120
    ## the slow term vanishes and the model is the exact monoexponential
    expect_equal(
        biexponential(
            t, A = 70, B = 40, tau = 5, B2 = 40, tau2 = 50
        ),
        monoexponential(t, A = 70, B = 40, tau = 5)
    )
})


## SSbiexponential() ================================================
test_that("SSbiexponential() converges on known parameters", {
    set.seed(1)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B, tau, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "tau", "B2", "tau2"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], 70, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B"]], 40, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["B2"]], 60, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["tau"]], 5, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["tau2"]], 40, tolerance = 15, scale = 1))
})

test_that("SSbiexponential() fits the 6-parameter TD form", {
    set.seed(4)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40,
        TD = 10
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B, tau, B2, tau2, TD),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_named(
        coef(model), c("A", "B", "tau", "B2", "tau2", "TD")
    )
    expect_true(all.equal(coef(model)[["TD"]], 10, tolerance = 3, scale = 1))
})

test_that("SSbiexponential() predict() returns correct length", {
    set.seed(2)
    t <- 0:120
    x <- biexponential(t, 70, 40, 5, 60, 40) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B, tau, B2, tau2),
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
        t, A = 0, B = -25, tau = 5, B2 = -10, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    suppressWarnings(
        model <- nls(
            x ~ SSbiexponential(t, A = 0, B, tau, B2, tau2),
            data = data,
            algorithm = "port",
            lower = c(-Inf, 0, -Inf, 0),
            control = nls.control(warnOnly = TRUE)
        )
    )

    expect_named(coef(model), c("B", "tau", "B2", "tau2"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 0)
})

test_that("SSbiexponential() gradient matches numericDeriv for the free parameters", {
    t <- seq(-10, 120, by = 0.5)
    env <- list2env(list(
        t = t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40, TD = 3
    ))
    chk <- function(expr, pars) {
        an <- attr(eval(expr, env), "gradient")
        nd <- attr(numericDeriv(expr, pars, env), "gradient")
        expect_identical(colnames(an), pars)
        expect_equal(unname(an), unname(nd), tolerance = 1e-5)
    }
    chk(
        quote(SSbiexponential(t, A, B, tau, B2, tau2, TD)),
        c("A", "B", "tau", "B2", "tau2", "TD")
    )
    ## a constant in the formula contributes no column
    chk(
        quote(SSbiexponential(t, A, B, tau = 5, B2, tau2)),
        c("A", "B", "B2", "tau2")
    )
    ## no free parameter, no gradient; the exported fn stays plain
    expect_null(attr(SSbiexponential(t, 70, 40, 5, 60, 40), "gradient"))
    expect_null(attr(biexponential(t, 70, 40, 5, 60, 40), "gradient"))
})

test_that("biexp_start() matches a per-point least-squares grid search", {
    set.seed(8)
    t <- 0:120
    x <- biexponential(t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    start <- biexp_start(x, t, has_TD = TRUE)
    expect_named(start, c("A", "B", "tau", "B2", "tau2", "TD"))

    ## the linear coefficients at the chosen grid point are the lm solution
    ts <- pmax(t - start[["TD"]], 0)
    e1 <- exp(-ts / start[["tau"]])
    e2 <- exp(-ts / start[["tau2"]])
    cf <- lm.fit(cbind(e1, e2 - e1, 1 - e2), x)$coefficients
    expect_equal(unname(start[c("A", "B", "B2")]), unname(cf), tolerance = 1e-8)
    expect_true(start[["tau2"]] >= start[["tau"]] / tau_ratio)
})

test_that("SSbiexponential() fixes tau in the formula", {
    set.seed(5)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B, tau = 5, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_named(coef(model), c("A", "B", "B2", "tau2"))
    expect_true(all.equal(coef(model)[["B"]], 40, tolerance = 4, scale = 1))
})


## analyse_biexponential() ==========================================

## helper: create excursion-recovery test data with known parameters
create_biexp_data <- function(
    A = 70,
    B = 40,
    tau = 5,
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
    x <- biexponential(t, A, B, tau, B2, tau2) +
        rnorm(n, 0, noise_sd)

    df <- setNames(data.frame(t, x), c("time", channels[1]))
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- biexponential(t, A + 5, B + 5, tau, B2 + 5, tau2) +
                rnorm(n, 0, noise_sd)
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
        "interval", "nirs_channels", "A", "B", "TD", "tau",
        "MRT", "texc", "B2", "tau2", "MRT_fitted", "texc_fitted"
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
    B <- 40
    B2 <- 60
    data <- create_biexp_data(A = A, B = B, B2 = B2, noise_sd = 0.3)
    
    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 3, scale = 1))
    ## the fast asymptote absorbs some slow phase over the stage-1 window
    expect_true(all.equal(result$B, B, tolerance = 5, scale = 1))
    expect_true(all.equal(result$B2, B2, tolerance = 3, scale = 1))
    ## no TD: MRT is the fast time constant
    expect_equal(result$MRT, result$tau)
    ## excursion sits inside the window, below the starting value
    expect_true(result$texc > 0)
    expect_true(result$texc_fitted < result$A)
    ## good fit
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
})

test_that("analyse_biexponential() texc is the fitted excursion point", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    ## texc_fitted is the model prediction at the fitted excursion point
    expect_equal(
        result$texc_fitted,
        biexponential(
            result$texc, result$A, result$B, result$tau,
            result$B2, result$tau2
        )
    )
    expect_true(result$texc > 0)
    expect_true(result$tau < result$tau2)
    ## the excursion point sits near the fitted-curve minimum
    fitted <- attr(result, "fitted_data")$smo2$fitted
    expect_true(all.equal(result$texc_fitted, min(fitted), tolerance = 1,
        scale = 1))
})

test_that("analyse_biexponential() reports tau <= tau2", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )
    expect_lte(result$tau, result$tau2)

    set.seed(6)
    t <- 0:120
    x <- biexponential(t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.3)
    data <- data.frame(t, x)
    model <- nls(
        x ~ SSbiexponential(t, A, B, tau, B2, tau2),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )
    expect_true(coef(model)[["tau"]] < coef(model)[["tau2"]])
})

test_that("analyse_biexponential() reports NA texc for a monotonic fit", {
    ## B between A and B2: the fitted curve is monotonic, so there is no
    ## interior excursion point
    set.seed(7)
    t <- 0:120
    x <- biexponential(t, A = 70, B = 55, tau = 5, B2 = 40, tau2 = 40) +
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
    expect_false(is.na(result$B))
    expect_false(is.na(result$tau))
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
    expect_true(is.na(result$tau))
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

    ## the stage-1 window drives tau
    full <- analyse_biexponential(
        data, nirs_channels = "smo2", use_TD = FALSE, verbose = FALSE
    )
    expect_false(isTRUE(all.equal(result$tau, full$tau)))
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
        tau_flex = 0.1,
        A_flex = 0.5,
        verbose = FALSE
    )

    expect_true(result$tau >= fast$tau / 1.1 - 1e-8)
    expect_true(result$tau <= fast$tau * 1.1 + 1e-8)
    expect_true(abs(result$A - fast$A) <= 0.5 + 1e-8)
    ## the slow phase separates above the fast-phase ceiling
    expect_true(result$tau2 >= fast$tau * 1.1 / tau_ratio - 1e-8)

    ## flex args pass through `...` of analyse_kinetics() and are recorded
    ca <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        end_window = 20,
        tau_flex = list(smo2 = 0.1),
        TD_flex = 1,
        verbose = FALSE
    )$channel_args
    expect_equal(ca$tau_flex, 0.1)
    expect_equal(ca$TD_flex, 1)
    expect_true(is.na(ca$A_flex))
    ca <- analyse_kinetics(
        data, nirs_channels = "smo2", method = "biexponential", verbose = FALSE
    )$channel_args
    expect_equal(ca$tau_flex, 1 / 3)
    expect_equal(ca$TD_flex, 2)
})

test_that("analyse_biexponential() fits a mirrored rise-overshoot response", {
    ## inverted kinetics: the response rises to a peak before settling back.
    ## the asymptote ordering flips (B above A, B2 below B) and is
    ## recovered from the data
    set.seed(31)
    t <- 0:119
    x <- biexponential(
        t, A = 30, B = 55, tau = 10, B2 = 40, tau2 = 50
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

    expect_true(all.equal(result$B, 55, tolerance = 5, scale = 1))
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
        c("A", "B", "tau", "B2", "tau2"))
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
        t - start_time, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40,
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
        c("A", "B", "tau", "B2", "tau2", "TD"))

    ## MRT is the fast-phase mean response time from the fit onset
    expect_equal(result$MRT, result$TD + result$tau)
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
    expect_equal(result_auto$tau, result_neg$tau)
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
    expect_true(is.na(result$tau))
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
    expect_true(is.na(result$tau))
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
    data <- create_biexp_data(A = 0, B = -20, B2 = -10, noise_sd = 0.3)

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
        c("B", "tau", "B2", "tau2")
    )
    expect_equal(attr(result, "channel_args")$fix, "list(A = 0)")
    expect_false(is.na(attr(result, "diagnostics")$adj_r2))
})

test_that("analyse_biexponential() fix holds tau constant", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(tau = 12),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$tau, 12)
    expect_named(
        coef(attr(result, "model")$smo2), c("A", "B", "B2", "tau2")
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
        coef(attr(result, "model")$smo2), c("A", "B", "tau", "B2")
    )
    expect_true(all.equal(result$tau, 5, tolerance = 2, scale = 1))
})

test_that("analyse_biexponential() returned model is canonical and converged", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        verbose = FALSE
    )
    model <- attr(result, "model")$smo2

    expect_named(coef(model), c("A", "B", "tau", "B2", "tau2", "TD"))
    expect_true(model$convInfo$isConv)
    expect_true(result$tau2 >= result$tau / tau_ratio)
    ## coefficient table mirrors the model
    expect_equal(result$tau, coef(model)[["tau"]])
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
        t, A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40,
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
        c("A", "B", "tau", "B2", "tau2")
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
        A = 0, B = -20, B2 = -10, noise_sd = 0.3,
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
    expect_named(coef(models$ch1), c("B", "tau", "B2", "tau2"))
    expect_named(coef(models$ch2), c("B", "tau", "B2", "tau2"))

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
                A = 70, B = 40, tau = 5, B2 = 60, tau2 = 40, TD = 0
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
        c("B", "tau", "B2", "tau2", "texc") %in%
            names(result$coefficients)
    ))
    ## the model carries its fit data in the call for update()/insight
    model <- result$model[[1L]]$smo2
    fit_data <- eval(model$call$data, envir = baseenv())
    expect_s3_class(fit_data, "data.frame")
    expect_named(fit_data, c("smo2", "time"))
    expect_equal(nrow(fit_data), length(stats::fitted(model)))
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
    data <- create_biexp_data(A = 0, B = -20, B2 = -10, noise_sd = 0.3)

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
        c("B", "tau", "B2", "tau2")
    )
})


## model fallback chain =============================================

## monotonic two-phase response: B between A and B2, no excursion point
create_monotonic_data <- function(seed = 7, TD = NULL, t = 0:119) {
    set.seed(seed)
    x <- biexponential(t, A = 70, B = 55, tau = 5, B2 = 40, tau2 = 40, TD) +
        rnorm(length(t), 0, 0.3)
    create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2",
        time_channel = "time",
        sample_rate = 1
    )
}

## fast drop then a linear recovery: a slow phase no record can resolve
create_linear_tail_data <- function(seed = 5, t = 0:119) {
    set.seed(seed)
    # fmt: skip
    x <- exponential_drift(
        t, A = 70, B = 40, tau = 5, slope = 0.1, tau_mult = 3
    ) +
        rnorm(length(t), 0, 0.3)
    create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2",
        time_channel = "time",
        sample_rate = 1
    )
}

fallback_models <- c("exponential_drift", "monoexponential")
biexp_only <- c("B2", "tau2")

test_that("analyse_kinetics() falls back from a monotonic biexponential fit", {
    data <- create_monotonic_data()

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "biexponential",
            use_TD = FALSE
        ),
        "fell back to"
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_equal(names(cf)[1:4], c("interval", "nirs_channels", "start_time", "model"))
    expect_true(cf$model %in% fallback_models)
    ## the union schema carries every model's columns, `_fitted` last
    expect_equal(names(cf)[-(1:4)], kinetics_chain_cols("biexponential"))
    ## the fallback row reports its own model's coefficients
    expect_true(all(is.na(cf[biexp_only])))
    expect_equal(cf$A, coef(model)[["A"]])
    expect_equal(cf$B, coef(model)[["B"]])
    expect_equal(cf$tau, coef(model)[["tau"]])
    expect_equal(cf$MRT, cf$tau)
    ## diagnostics and fitted values follow the fallback model
    expect_equal(result$diagnostics$n_params, length(coef(model)))
    expect_equal(
        result$data[[1L]]$smo2_fitted,
        as.vector(predict(model))
    )
    ## recorded with the reason, regardless of verbose
    msgs <- result$warnings$message
    expect_true(any(grepl("fell back to", msgs)))
    expect_true(any(grepl("monotonic", msgs)))
})

test_that("analyse_kinetics() falls back to a monoexponential response", {
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
    cf <- result$coefficients

    expect_equal(cf$model, "monoexponential")
    expect_true(all.equal(cf$tau, 8, tolerance = 1, scale = 1))
    expect_true(all(is.na(cf[c(biexp_only, "slope", "texc")])))
    expect_true(inherits(result$model[[1L]]$smo2, "nls"))
    ## both fallbacks are recorded
    expect_equal(sum(grepl("fell back to", result$warnings$message)), 2L)
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
        model_fallback = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "biexponential")
    expect_true(is.finite(cf$texc))
    expect_true(is.na(cf$slope))
    ## the fast phase shares the primary-phase columns, not its own pair
    expect_false(any(c("B1", "tau1") %in% names(cf)))
    expect_named(
        coef(result$model[[1L]]$smo2),
        c("A", "B", "tau", "B2", "tau2", "TD")
    )
    expect_false(any(grepl("fell back to", result$warnings$message)))
    ## the raw fit is unchanged by the triggers
    expect_equal(cf, forced$coefficients)
    expect_equal(result$diagnostics, forced$diagnostics)
})

test_that("model_fallback = FALSE keeps a monotonic fit", {
    data <- create_monotonic_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        model_fallback = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "biexponential")
    expect_true(is.na(cf$texc))
    expect_false(is.na(cf$tau2))
    expect_false(cf$B == cf$B2)
    expect_false(any(grepl("fell back to", result$warnings$message)))
})

test_that("a slow phase beyond the record falls back to exponential_drift", {
    data <- create_linear_tail_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$model, "exponential_drift")
    expect_true(all.equal(cf$slope, 0.1, tolerance = 0.05, scale = 1))
    expect_true(is.finite(cf$texc))
    expect_named(coef(result$model[[1L]]$smo2), c("A", "B", "tau", "slope"))
    expect_true(all(is.na(cf[biexp_only])))
    msgs <- result$warnings$message
    expect_true(any(grepl("fell back to", msgs)))
    expect_true(any(grepl("tau2 exceeds", msgs)))
})

test_that("fallback fits the full response, not end_window", {
    data <- create_monotonic_data(seed = 11, TD = 10, t = -20:120)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        end_window = 30,
        verbose = FALSE
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_true(cf$model %in% fallback_models)
    expect_equal(is.finite(cf$TD), "TD" %in% names(coef(model)))
    expect_gte(result$diagnostics$n_obs, sum(data$time >= 0))
    expect_equal(
        result$diagnostics$n_obs,
        sum(is.finite(result$data[[1L]]$smo2_fitted))
    )
    ## resolved arguments are the fallback fit's
    expect_equal(result$channel_args$end_window, Inf)
})

test_that("fallback carries fixed parameters over", {
    data <- create_monotonic_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        fix = list(A = 70, B = 55),
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_true(cf$model %in% fallback_models)
    expect_equal(cf$A, 70)
    expect_equal(cf$B, 55)
    expect_false(any(c("A", "B") %in% names(coef(result$model[[1L]]$smo2))))

    ## parameters without a counterpart are dropped
    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        fix = list(tau2 = 40),
        verbose = FALSE
    )

    expect_true(result$coefficients$model %in% fallback_models)
    expect_true(all(c("A", "B", "tau") %in% names(coef(result$model[[1L]]$smo2))))
})

test_that("fallback resolves per interval", {
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
    expect_equal(cf$model[[1L]], "biexponential")
    expect_true(cf$model[[2L]] %in% fallback_models)
    expect_named(
        coef(result$model$excursion$smo2), c("A", "B", "tau", "B2", "tau2")
    )
    expect_true(is.finite(cf$texc[[1L]]))
    expect_true(is.na(cf$tau2[[2L]]))
    fell <- result$warnings[grepl("fell back to", result$warnings$message), ]
    expect_true(all(fell$interval == "monotonic"))
    expect_true(all(fell$nirs_channels == "smo2"))
})

test_that("fallback resolves per channel", {
    data <- create_biexp_data(noise_sd = 0.3, channels = c("smo2", "hhb"))
    data$hhb <- create_monotonic_data()$smo2

    ## the per-channel map is keyed to the fallback channel only
    expect_no_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = c(smo2, hhb),
            method = "biexponential",
            use_TD = list(smo2 = TRUE, hhb = FALSE),
            fix = list(smo2 = list(A = 70), hhb = list(A = 70)),
            verbose = FALSE
        )
    )
    cf <- result$coefficients

    expect_equal(cf$nirs_channels, c("smo2", "hhb"))
    expect_equal(cf$model[[1L]], "biexponential")
    expect_true(cf$model[[2L]] %in% fallback_models)
    expect_true(is.finite(cf$TD[[1L]]))
    expect_true(is.na(cf$TD[[2L]]))
    expect_equal(cf$A, c(70, 70))
    expect_true(inherits(result$model[[1L]]$smo2, "nls"))
    expect_true(inherits(result$model[[1L]]$hhb, "nls"))
    expect_equal(result$diagnostics$nirs_channels, c("smo2", "hhb"))
    expect_equal(result$channel_args$nirs_channels, c("smo2", "hhb"))
    expect_equal(
        sum(is.finite(result$data[[1L]]$hhb_fitted)),
        result$diagnostics$n_obs[[2L]]
    )
})

test_that("a row where every fit fails reports the last method", {
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

    expect_equal(cf$model, "monoexponential")
    expect_true(is.na(cf$A))
    expect_null(result$model[[1L]]$smo2)
    expect_equal(sum(grepl("fell back to", result$warnings$message)), 2L)
})

test_that("keep_fix() filters through nested maps", {
    keep <- kinetics_fallbacks$biexponential$fix_keep

    expect_equal(keep_fix(list(A = 1, B = 2, tau2 = 3), keep), list(A = 1, B = 2))
    expect_equal(keep_fix(list(tau2 = 3), keep), setNames(list(), character()))
    expect_equal(
        keep_fix(list(smo2 = list(B = 2), hhb = list(tau2 = 4)), keep),
        list(smo2 = list(B = 2), hhb = setNames(list(), character()))
    )
})

test_that("kinetics_chain_cols() unions the chain with `_fitted` columns last", {
    expect_equal(
        kinetics_chain_cols("biexponential"),
        c(
            "A", "B", "TD", "tau", "MRT", "texc", "B2", "tau2", "k", "HRT",
            "slope", "tau_mult", "MRT_fitted", "texc_fitted", "HRT_fitted"
        )
    )
    expect_equal(
        kinetics_chain_cols("monoexponential"),
        kinetics_coef_cols$monoexponential
    )
})

test_that("bind_union() pads differing columns and keeps a fixed order", {
    a <- data.frame(x = 1, y = "a")
    b <- data.frame(y = "b", z = TRUE)

    out <- bind_union(list(a, b))
    expect_named(out, c("x", "y", "z"))
    expect_equal(out$x, c(1, NA))
    expect_equal(out$z, c(NA, TRUE))
    expect_named(bind_union(list(a, b), c("z", "y", "x", "w")), c("z", "y", "x", "w"))
    expect_equal(bind_union(list(a[0, ], b))$y, "b")
    expect_null(bind_union(list(a[0, ])))
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
    expect_equal(coefs$model[1:2], c("biexponential", "biexponential"))
    expect_true(all(coefs$model[3:4] %in% fallback_models))
    expect_true(all(!is.na(coefs$A)))

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
    (success <- mean(!is.na(coefs$tau)))
    expect_true(all(success >= 1.0))

    (texc_success <- mean(!is.na(coefs$texc)))
    expect_true(all(success >= 0.75))
})
