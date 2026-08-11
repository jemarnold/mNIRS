## biexponential() ==================================================
test_that("biexponential() returns correct vector length", {
    t <- 0:120
    result <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("biexponential() starts at A and approaches the plateau", {
    t <- 0:500
    A <- 70
    B1 <- 45
    B2 <- 60
    result <- biexponential(t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40)

    expect_equal(result[1], A)
    expect_true(
        all.equal(result[length(result)], B2, tolerance = 0.01,
            scale = 1)
    )
})

test_that("biexponential() drops to a nadir below A and below the plateau", {
    t <- 0:120
    A <- 70
    B1 <- 45
    B2 <- 60
    result <- biexponential(t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40)

    ## interior minimum below both endpoints (nadir-recovery shape)
    expect_true(min(result) < A)
    expect_true(min(result) < B2)
    expect_true(which.min(result) > 1 && which.min(result) < length(result))
})

test_that("biexponential() TD form is flat before the delay", {
    t <- 0:120
    TD <- 15
    result <- biexponential(
        t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = TD
    )

    expect_true(all(result[t < TD] == 70))
    expect_equal(result[t == TD], 70)
})

test_that("biexponential() reduces to the monoexponential when B1 = B2", {
    t <- 0:120
    ## the slow term vanishes and the model is the exact monoexponential
    expect_equal(
        biexponential(t, A = 70, B1 = 40, tau1 = 5, B2 = 40, tau2 = 50),
        monoexponential(t, A = 70, B = 40, tau = 5)
    )
})


## SSbiexponential() ================================================
test_that("SSbiexponential() converges on known parameters", {
    set.seed(1)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B1, lt1, B2, lr), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B1", "lt1", "B2", "lr"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], 70, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B1"]], 45, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["B2"]], 60, tolerance = 4, scale = 1))
    ## natural time constants recovered from the log-ratio scale
    expect_true(all.equal(exp(coefs[["lt1"]]), 5, tolerance = 2, scale = 1))
    expect_true(all.equal(
        exp(coefs[["lt1"]] + coefs[["lr"]]), 40, tolerance = 10, scale = 1
    ))
})

test_that("SSbiexponential() fits the 6-parameter TD form", {
    set.seed(4)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = 10
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, lt1, B2, lr, TD), data = data
    )

    expect_named(coef(model), c("A", "B1", "lt1", "B2", "lr", "TD"))
    expect_true(all.equal(coef(model)[["TD"]], 10, tolerance = 3, scale = 1))
})

test_that("SSbiexponential() predict() returns correct length", {
    set.seed(2)
    t <- 0:120
    x <- biexponential(t, 70, 45, 5, 60, 40) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B1, lt1, B2, lr), data = data)

    expect_length(predict(model, data), nrow(data))
})

test_that("SSbiexponential() fixes A at a constant", {
    set.seed(3)
    t <- 0:120
    x <- biexponential(t, A = 0, B1 = -25, tau1 = 5, B2 = -10, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A = 0, B1, lt1, B2, lr), data = data
    )

    expect_named(coef(model), c("B1", "lt1", "B2", "lr"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 0)
})

test_that("SSbiexponential() with port matches the package fit path", {
    set.seed(5)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)
    span <- diff(range(t))

    ## the port + bounds incantation documented in ?SSbiexponential
    model <- nls(
        x ~ SSbiexponential(t, A, B1, lt1, B2, lr),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, log(span / 1000), -Inf, log(2.5))
    )

    mnirs_data <- create_mnirs_data(
        setNames(data, c("time", "smo2")),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )
    result <- analyse_biexponential(
        mnirs_data, nirs_channels = "smo2", use_TD = FALSE, verbose = FALSE
    )

    ## both routes reach the same optimum
    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], result$A, tolerance = 0.1, scale = 1))
    expect_true(all.equal(
        exp(coefs[["lt1"]]), result$tau1, tolerance = 0.1, scale = 1
    ))
    expect_true(all.equal(
        exp(coefs[["lt1"]] + coefs[["lr"]]), result$tau2,
        tolerance = 0.5, scale = 1
    ))
})


## biexp_grid_start() ===============================================
test_that("biexp_grid_start() recovers known time constants", {
    t <- 0:200
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 10, B2 = 60, tau2 = 60)

    seed <- biexp_grid_start(x, t)

    expect_named(seed, c("A", "B1", "tau1", "B2", "tau2", "rss"))
    ## grid resolution is coarse by design; nls polishes from the seed
    expect_true(seed$tau1 > 10 / 1.5 && seed$tau1 < 10 * 1.5)
    expect_true(seed$tau2 > 60 / 1.5 && seed$tau2 < 60 * 1.5)
    expect_true(all.equal(seed$A, 70, tolerance = 5, scale = 1))
    expect_true(all.equal(seed$B1, 45, tolerance = 5, scale = 1))
    expect_true(all.equal(seed$B2, 60, tolerance = 5, scale = 1))
    ## noiseless data: the profiled residual is a small share of total SS
    expect_true(seed$rss >= 0)
    expect_true(seed$rss < sum((x - mean(x))^2) * 0.01)
})

test_that("biexp_grid_start() confines the grid to the ratio ridge", {
    t <- 0:200
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 10, B2 = 60, tau2 = 60)

    default <- biexp_grid_start(x, t)
    bounded <- biexp_grid_start(x, t, tau_ratio = 8)

    expect_true(default$tau2 >= default$tau1 * 2.5)
    expect_true(bounded$tau2 >= bounded$tau1 * 8)
})

test_that("biexp_grid_start() is invariant to the time units of t", {
    t <- 0:200
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 10, B2 = 60, tau2 = 60)

    ## grid limits scale with the span of t, so seeding in minutes rather
    ## than seconds rescales the time constants and leaves amplitudes alone
    secs <- biexp_grid_start(x, t)
    mins <- biexp_grid_start(x, t / 60)

    expect_equal(mins$tau1, secs$tau1 / 60)
    expect_equal(mins$tau2, secs$tau2 / 60)
    expect_equal(mins$A, secs$A)
    expect_equal(mins$B1, secs$B1)
})

test_that("biexp_grid_start() returns NULL for a degenerate predictor", {
    ## no span to grid over
    expect_null(biexp_grid_start(1:5, rep(0, 5)))
    ## non-finite span
    expect_null(biexp_grid_start(1:5, c(0, NA, 2, 3, 4)))
})


## analyse_biexponential() ==========================================

## helper: create excursion-recovery test data with known parameters
create_biexp_data <- function(
    A = 70,
    B1 = 45,
    tau1 = 15,
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
    x <- biexponential(t, A, B1, tau1, B2, tau2) + rnorm(n, 0, noise_sd)

    df <- setNames(data.frame(t, x), c("time", channels[1]))
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- biexponential(t, A + 5, B1 + 5, tau1, B2 + 5, tau2) +
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
    # plot(data)
    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        verbose = TRUE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "time_channel",
        "A", "B1", "tau1", "B2", "tau2", "TD",
        "excursion_time", "excursion_value"
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
    B1 <- 45
    B2 <- 60
    result <- analyse_biexponential(
        create_biexp_data(A = A, B1 = B1, B2 = B2, noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B1, B1, tolerance = 3, scale = 1))
    expect_true(all.equal(result$B2, B2, tolerance = 3, scale = 1))
    ## excursion sits inside the window, below the starting value
    expect_true(result$excursion_time > 0)
    expect_true(result$excursion_value < result$A)
    ## good fit
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
})

test_that("analyse_biexponential() excursion matches the fitted curve minimum", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    fitted <- attr(result, "fitted_data")$smo2$fitted
    ## closed-form excursion agrees with the numeric minimum of the fit
    expect_true(all.equal(result$excursion_value, min(fitted), tolerance = 0.1,
        scale = 1))
})

test_that("analyse_biexponential() uses start_time correctly", {
    start_time <- 12
    data <- create_biexp_data(noise_sd = 0.3)
    data$time <- data$time + start_time
    # plot(data)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, 70, tolerance = 2, scale = 1))
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
            nirs_channels = "smo2"
        ),
        "fit failed for.*smo2.*custom_name"
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$tau1))
    expect_true(is.na(result$excursion_value))
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

test_that("analyse_biexponential() validates tau_ratio", {
    data <- create_biexp_data()

    ## a ratio of 1 admits the singular tau2 == tau1 design
    expect_error(
        analyse_biexponential(data, nirs_channels = "smo2", tau_ratio = 1),
        "tau_ratio"
    )
    expect_error(
        analyse_biexponential(data, nirs_channels = "smo2", tau_ratio = "2"),
        "tau_ratio"
    )
    expect_error(
        analyse_biexponential(
            data, nirs_channels = "smo2", tau_ratio = c(2, 3)
        ),
        "tau_ratio"
    )
})

test_that("analyse_biexponential() end_window truncates the fit window", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        end_window = 40,
        use_TD = FALSE,
        verbose = FALSE
    )

    fitted_data <- attr(result, "fitted_data")$smo2
    expect_true(nrow(fitted_data) < nrow(data))
    expect_equal(attr(result, "diagnostics")$n_obs, nrow(fitted_data))
})

test_that("analyse_biexponential() fits a mirrored rise-overshoot response", {
    ## inverted kinetics: the response rises to a peak before settling back.
    ## the asymptote ordering flips (B1 above A, B2 below B1) and is
    ## recovered from the data
    set.seed(31)
    t <- 0:119
    x <- biexponential(t, A = 30, B1 = 55, tau1 = 10, B2 = 40, tau2 = 50) +
        rnorm(120, 0, 0.3)
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
    ## the closed-form root needs B1 beyond both A and B2, whichever the
    ## direction, so a mirrored response reports a genuine interior
    ## excursion: a maximum above the baseline
    fitted <- attr(result, "fitted_data")$smo2$fitted
    expect_true(result$excursion_time > 0)
    expect_true(result$excursion_value > result$A)
    expect_true(all.equal(result$excursion_value, max(fitted), tolerance = 0.1,
        scale = 1))
})


test_that("analyse_biexponential() excursion falls back on a monotone fit", {
    ## a B1 between A and B2 makes the curve monotone: no interior root
    ## exists and the turning point sits at the onset boundary
    set.seed(37)
    t <- 0:119
    x <- biexponential(t, A = 70, B1 = 45, tau1 = 5, B2 = 30, tau2 = 50) +
        rnorm(120, 0, 0.3)
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

    skip_if(is.na(result$tau1), "fit did not converge")
    expect_true((result$A - result$B1) * (result$B2 - result$B1) < 0)
    expect_equal(result$excursion_time, 0)
    expect_equal(result$excursion_value, result$A)
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

    ## the stored model is fit in ratio parameterisation (lt1, lr); the
    ## reported coefficients are converted back to the natural tau scale
    expect_named(coef(attr(result, "model")$smo2),
        c("A", "B1", "lt1", "B2", "lr"))
    expect_true(is.na(result$TD))
})

test_that("analyse_biexponential() falls back to the 5-parameter fit", {
    ## six observations under-determine the 6-parameter model but not the
    ## 5-parameter model, so the TD fit is rejected and the retry announced
    data <- create_biexp_data(n = 6, noise_sd = 0.1)

    warns <- capture_warnings(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = TRUE
        )
    )

    expect_match(warns[[1L]], "6-parameter")
    expect_match(warns[[1L]], "Attempting")
    ## TD is absent from the reduced model whether or not the retry converges
    expect_true(is.na(result$TD))
})

test_that("analyse_biexponential() reports natural time constants", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    ## tau1/tau2 recover exp(lt1) and exp(lt1 + lr) from the stored model
    coefs <- coef(attr(result, "model")$smo2)
    expect_equal(result$tau1, unname(exp(coefs[["lt1"]])))
    expect_equal(result$tau2, unname(exp(coefs[["lt1"]] + coefs[["lr"]])))
})

test_that("analyse_biexponential() separates the fast and slow components", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    ## the ratio bound keeps the two components identifiable; without it the
    ## time constants collapse together and the amplitudes cancel
    expect_true(result$tau2 >= result$tau1 * 2.5 - 1e-6)
    expect_true(result$tau1 > 0)
})

test_that("analyse_biexponential() tau_ratio bounds the time constants", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        tau_ratio = 5,
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(result$tau2 >= result$tau1 * 5 - 1e-6)
})

test_that("analyse_biexponential() excursion_time is elapsed from start_time", {
    ## 6-parameter TD fit with a non-zero start_time: excursion_time must be
    ## measured from start_time, not from the model's internal (TD) onset
    set.seed(11)
    start_time <- 20
    TD <- 8
    t <- start_time + 0:119
    x <- biexponential(
        t - start_time, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = TD
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
        c("A", "B1", "lt1", "B2", "lr", "TD"))

    ## excursion_time equals the elapsed time from start_time to the fitted
    ## minimum; the old onset-relative value omits the ~TD offset
    fd <- attr(result, "fitted_data")$smo2
    t_at_min <- data$time[fd$window_idx[which.min(fd$fitted)]]
    expect_true(all.equal(result$excursion_time, t_at_min - start_time,
        tolerance = 1.5, scale = 1))
    expect_true(result$excursion_time > TD)
})

test_that("analyse_biexponential() excursion_value matches the fitted minimum with TD", {
    ## origin-invariant evaluation: excursion_value must equal the fitted-curve
    ## minimum even when TD and a non-zero start_time shift the time frame
    set.seed(12)
    start_time <- 20
    TD <- 8
    t <- start_time + 0:119
    x <- biexponential(
        t - start_time, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = TD
    ) + rnorm(120, 0, 0.15)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        use_TD = TRUE
    )

    expect_named(coef(attr(result, "model")$smo2),
        c("A", "B1", "lt1", "B2", "lr", "TD"))

    ## closed-form excursion agrees with the numeric minimum of the fit
    fitted <- attr(result, "fitted_data")$smo2$fitted
    expect_true(all.equal(result$excursion_value, min(fitted), tolerance = 0.15,
        scale = 1))
    expect_true(result$excursion_value < result$A)
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

test_that("analyse_biexponential() refits a monotone decline in-direction", {
    ## monoexponential-shaped data has no rebound: the unconstrained fit can
    ## mirror the components (B1 above A) to mimic a sigmoidal decline. the
    ## direction refit bounds B1 onto the response side, yielding a valid
    ## two-phase decline rather than NA
    set.seed(41)
    t <- 0:119
    x <- monoexponential(t, A = 70, B = 40, tau = 20) + rnorm(120, 0, 0.3)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    expect_no_warning(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            use_TD = FALSE,
            verbose = TRUE
        )
    )

    expect_false(is.na(result$tau1))
    expect_true(result$B1 < result$A)
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
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
        coef(attr(result, "model")$smo2), c("B1", "lt1", "B2", "lr")
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
    ## a fixed time constant drops its log-scale counterpart from the model
    expect_named(coef(attr(result, "model")$smo2), c("A", "B1", "B2", "lr"))
    ## the ratio bound still applies about the fixed value
    expect_true(result$tau2 >= 12 * 2.5 - 1e-6)
})

test_that("analyse_biexponential() fix holds tau2 constant", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(tau2 = 50),
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_equal(result$tau2, 50)
    expect_named(coef(attr(result, "model")$smo2), c("A", "B1", "lt1", "B2"))
    ## a fixed tau2 caps tau1 at tau2 / tau_ratio from the other side
    expect_true(result$tau1 <= 50 / 2.5 + 1e-6)
})

test_that("analyse_biexponential() fix holds TD constant", {
    set.seed(21)
    TD <- 8
    t <- 0:119
    x <- biexponential(
        t, A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = TD
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
        coef(attr(result, "model")$smo2), c("A", "B1", "lt1", "B2", "lr")
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
            fix = list(A = 70, B1 = 45, tau1 = 5, B2 = 60, tau2 = 40, TD = 0)
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
        c("B1", "tau1", "B2", "tau2", "excursion_time") %in%
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
        coef(result$model[[1L]]$smo2), c("B1", "lt1", "B2", "lr")
    )
})

test_that("analyse_kinetics() passes use_TD to the biexponential method", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_named(
        coef(result$model[[1L]]$smo2), c("A", "B1", "lt1", "B2", "lr")
    )
})


## integration tests ================================================

test_that("analyse_biexponential() converges on real dataset", {
    skip("Manual fit convergence check")

    ## 5 real-world 5-1 min work-recovery intervals, 4 channels each    
    intervals <- readRDS(test_path("testdata/5-1_intervals_short.rds"))
    analyse_kinetics(
        intervals,
        method = "biexp",
        use_TD = TRUE,
        # verbose = FALSE
    # )
    # ) |> coef()
    ) |> plot(label = FALSE)

    # lapply(intervals, \(.df) {
    #     analyse_biexponential(
    #         .df,
    #         nirs_channels = smo2_left_rf,
    #         time_channel = time,
    #         use_TD = TRUE,
    #     )
    # })


    ## fit one signal at a time across all intervals; report convergence
    ## success rate per signal to flag regressions on real interval data.
    fit_5param <- function(signal) {
        vapply(intervals, \(df) {
            data <- data.frame(t = df$time, x = df[[signal]])
            model <- tryCatch(
                analyse_biexponential(
                    data,
                    nirs_channels = x,
                    time_channel = t,
                    use_TD = FALSE,
                ),
                error = \(e) NULL
            )
            as.integer(!is.null(model))
        }, integer(1L))
    }

    vl_left_success <- mean(fit_5param("smo2_left_vl"))
    vl_left_success
    vl_right_success <- mean(fit_5param("smo2_right_vl"))
    vl_right_success
    rf_left_success <- mean(fit_5param("smo2_left_rf"))
    rf_left_success
    rf_right_success <- mean(fit_5param("smo2_right_rf"))
    rf_right_success

    expect_true(vl_left_success >= 1.0)
    expect_true(vl_right_success >= 1.0)
    expect_true(rf_left_success >= 1.0)
    expect_true(rf_right_success >= 1.0)

    fit_6param <- function(signal) {
        vapply(intervals, \(df) {
            data <- data.frame(t = df$time, x = df[[signal]])
            model <- tryCatch(
                analyse_biexponential(
                    data,
                    nirs_channels = x,
                    time_channel = t,
                    use_TD = TRUE,
                ),
                error = \(e) NULL,
                warning = \(w) NULL
            )
            as.integer(!is.null(model))
        }, integer(1L))
    }

    vl_left_success <- mean(fit_6param("smo2_left_vl"))
    vl_left_success
    vl_right_success <- mean(fit_6param("smo2_right_vl"))
    vl_right_success
    rf_left_success <- mean(fit_6param("smo2_left_rf"))
    rf_left_success
    rf_right_success <- mean(fit_6param("smo2_right_rf"))
    rf_right_success

    expect_true(vl_left_success >= 0.8)
    expect_true(vl_right_success >= 0.8)
    expect_true(rf_left_success >= 0.8)
    expect_true(rf_right_success >= 0.8)
})

test_that("analyse_biexponential() converges on real dataset", {
    skip("Manual fit convergence check")

    intervals <- readRDS(test_path("testdata/5-1_intervals_short.rds"))
    nirs_channels <- c(
        "smo2_left_vl", "smo2_right_vl", "smo2_left_rf", "smo2_right_rf"
    )

    ## end-to-end path: window detection, ratio-parameterised fit and
    ## back-conversion. start_time = 0 anchors the fit at the interval onset
    results <- lapply(intervals, \(df) {
        analyse_biexponential(
            df,
            nirs_channels = nirs_channels,
            start_time = 0,
            use_TD = TRUE,
            verbose = FALSE
        )
    })

    coefs <- do.call(rbind, results)
    success <- tapply(!is.na(coefs$tau1), coefs$nirs_channels, mean)
    success
    expect_true(all(success >= 1.0))

    TD_success <- tapply(!is.na(coefs$TD), coefs$nirs_channels, mean)
    TD_success
    expect_true(all(TD_success >= 0.8))

    ## converged fits should describe the desaturation-recovery shape and
    ## keep the two components separated by the default ratio bound. the
    ## excursion only lies below the baseline for a downward excursion, so
    ## the direction check is restricted to that ordering
    ok <- !is.na(coefs$tau1)
    expect_true(all(coefs$tau2[ok] >= coefs$tau1[ok] * 2.5 - 1e-6))
    down <- ok & coefs$B1 < coefs$A & coefs$B2 > coefs$B1
    expect_true(all(coefs$excursion_value[down] <= coefs$A[down]))

    r2 <- unlist(lapply(results, \(x) attr(x, "diagnostics")$r2))
    mean(r2, na.rm = TRUE)
    expect_true(mean(r2 > 0.9, na.rm = TRUE) >= 0.8)
})
