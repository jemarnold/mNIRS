## biexponential() ==================================================
test_that("biexponential() returns correct vector length", {
    t <- 0:120
    result <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("biexponential() starts at A and approaches the plateau", {
    t <- 0:500
    A <- 70
    B1 <- 25
    B2 <- 15
    result <- biexponential(t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40)

    expect_equal(result[1], A)
    ## plateau = A - B1 + B2
    expect_true(
        all.equal(result[length(result)], A - B1 + B2, tolerance = 0.01,
            scale = 1)
    )
})

test_that("biexponential() drops to a nadir below A and below the plateau", {
    t <- 0:120
    A <- 70
    B1 <- 25
    B2 <- 15
    result <- biexponential(t, A = A, B1 = B1, tau1 = 5, B2 = B2, tau2 = 40)
    plateau <- A - B1 + B2

    ## interior minimum below both endpoints (nadir-recovery shape)
    expect_true(min(result) < A)
    expect_true(min(result) < plateau)
    expect_true(which.min(result) > 1 && which.min(result) < length(result))
})

test_that("biexponential() TD form is flat before the delay", {
    t <- 0:120
    TD <- 15
    result <- biexponential(
        t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40, TD = TD
    )

    expect_true(all(result[t < TD] == 70))
    expect_equal(result[t == TD], 70)
})


## SSbiexponential() ================================================
test_that("SSbiexponential() converges on known parameters", {
    set.seed(1)
    t <- 0:120
    x <- biexponential(t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B1, tau1, B2, tau2), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B1", "tau1", "B2", "tau2"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], 70, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B1"]], 25, tolerance = 4, scale = 1))
    expect_true(all.equal(coefs[["B2"]], 15, tolerance = 4, scale = 1))
})

test_that("SSbiexponential() fits the 6-parameter TD form", {
    set.seed(4)
    t <- 0:120
    x <- biexponential(
        t, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40, TD = 10
    ) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A, B1, tau1, B2, tau2, TD), data = data
    )

    expect_named(coef(model), c("A", "B1", "tau1", "B2", "tau2", "TD"))
    expect_true(all.equal(coef(model)[["TD"]], 10, tolerance = 3, scale = 1))
})

test_that("SSbiexponential() predict() returns correct length", {
    set.seed(2)
    t <- 0:120
    x <- biexponential(t, 70, 25, 5, 15, 40) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B1, tau1, B2, tau2), data = data)

    expect_length(predict(model, data), nrow(data))
})

test_that("SSbiexponential() fixes A at a constant", {
    set.seed(3)
    t <- 0:120
    x <- biexponential(t, A = 0, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A = 0, B1, tau1, B2, tau2), data = data
    )

    expect_named(coef(model), c("B1", "tau1", "B2", "tau2"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 0)
})


## analyse_biexponential() ==========================================

## helper: create nadir-recovery test data with known parameters
create_biexp_data <- function(
    A = 70,
    B1 = 25,
    tau1 = 15,
    B2 = 15,
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
            df[[ch]] <- biexponential(t, A + 5, B1, tau1, B2, tau2) +
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
        "plateau", "nadir_time", "nadir_value"
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
    B1 <- 25
    B2 <- 15
    result <- analyse_biexponential(
        create_biexp_data(A = A, B1 = B1, B2 = B2, noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B1, B1, tolerance = 3, scale = 1))
    expect_true(all.equal(result$B2, B2, tolerance = 3, scale = 1))
    ## plateau = A - B1 + B2
    expect_true(all.equal(result$plateau, result$A - result$B1 + result$B2,
        tolerance = 1e-6, scale = 1))
    ## nadir sits inside the window, below the starting value
    expect_true(result$nadir_time > 0)
    expect_true(result$nadir_value < result$A)
    ## good fit
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
})

test_that("analyse_biexponential() nadir matches the fitted curve minimum", {
    result <- analyse_biexponential(
        create_biexp_data(noise_sd = 0.3),
        nirs_channels = "smo2",
        use_TD = FALSE,
        verbose = FALSE
    )

    fitted <- attr(result, "fitted_data")$smo2$fitted
    ## closed-form nadir agrees with the numeric minimum of the fit
    expect_true(all.equal(result$nadir_value, min(fitted), tolerance = 0.1,
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
    expect_true(all.equal(result$plateau, 60, tolerance = 3, scale = 1))
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
    expect_true(is.na(result$plateau))
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

test_that("analyse_biexponential() nadir_time is elapsed from start_time", {
    ## 6-parameter TD fit with a non-zero start_time: nadir_time must be
    ## measured from start_time, not from the model's internal (TD) onset
    set.seed(11)
    start_time <- 20
    TD <- 8
    t <- start_time + 0:119
    x <- biexponential(
        t - start_time, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40, TD = TD
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

    ## nadir_time equals the elapsed time from start_time to the fitted
    ## minimum; the old onset-relative value omits the ~TD offset
    fd <- attr(result, "fitted_data")$smo2
    t_at_min <- data$time[fd$window_idx[which.min(fd$fitted)]]
    expect_true(all.equal(result$nadir_time, t_at_min - start_time,
        tolerance = 1.5, scale = 1))
    expect_true(result$nadir_time > TD)
})

test_that("analyse_biexponential() nadir_value matches the fitted minimum with TD", {
    ## origin-invariant evaluation: nadir_value must equal the fitted-curve
    ## minimum even when TD and a non-zero start_time shift the time frame
    set.seed(12)
    start_time <- 20
    TD <- 8
    t <- start_time + 0:119
    x <- biexponential(
        t - start_time, A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40, TD = TD
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

    expect_named(coef(attr(result, "model")$smo2),
        c("A", "B1", "tau1", "B2", "tau2", "TD"))

    ## closed-form nadir agrees with the numeric minimum of the fit
    fitted <- attr(result, "fitted_data")$smo2$fitted
    expect_true(all.equal(result$nadir_value, min(fitted), tolerance = 0.15,
        scale = 1))
    expect_true(result$nadir_value < result$A)
})


## direction ========================================================

test_that("analyse_biexponential() direction steers the fit window", {
    ## direction is used for window detection only; the non-monotone model
    ## is not rejected by a mismatched direction
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


## fixed parameters =================================================

test_that("analyse_biexponential() fix holds parameters constant", {
    data <- create_biexp_data(A = 0, B1 = 20, B2 = 10, noise_sd = 0.3)

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
        coef(attr(result, "model")$smo2), c("B1", "tau1", "B2", "tau2")
    )
    expect_equal(attr(result, "channel_args")$fix, "list(A = 0)")
    expect_false(is.na(attr(result, "diagnostics")$adj_r2))
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
            fix = list(A = 70, B1 = 25, tau1 = 5, B2 = 15, tau2 = 40, TD = 0)
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
        c("B1", "tau1", "B2", "tau2", "plateau", "nadir_time") %in%
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
    data <- create_biexp_data(A = 0, B1 = 20, B2 = 10, noise_sd = 0.3)

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
        coef(result$model[[1L]]$smo2), c("B1", "tau1", "B2", "tau2")
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
        coef(result$model[[1L]]$smo2), c("A", "B1", "tau1", "B2", "tau2")
    )
})
