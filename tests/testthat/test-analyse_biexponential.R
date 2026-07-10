## biexponential() ==================================================
test_that("biexponential() returns correct vector length", {
    t <- 0:120
    result <- biexponential(t, A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("biexponential() starts at baseline A and approaches asymptote B", {
    t <- 0:500
    A <- 50
    B <- 80
    result <- biexponential(t, A = A, B = B, tau1 = 5, tau2 = 40, prop = 0.6)

    expect_equal(result[1], A)
    expect_true(
        all.equal(result[length(result)], B, tolerance = 0.01, scale = 1)
    )
})

test_that("biexponential() handles rising (B > A) and falling (B < A)", {
    t <- 0:120
    rising <- biexponential(t, A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6)
    falling <- biexponential(t, A = 80, B = 50, tau1 = 5, tau2 = 40, prop = 0.6)

    expect_true(all(diff(rising) >= 0))
    expect_true(all(diff(falling) <= 0))
})

test_that("biexponential() clamps prop to the unit interval", {
    t <- 0:120
    ## out-of-range prop is clamped, so equivalent to a boundary value
    high <- biexponential(t, 50, 80, tau1 = 5, tau2 = 40, prop = 5)
    one <- biexponential(t, 50, 80, tau1 = 5, tau2 = 40, prop = 1 - 1e-6)

    expect_equal(high, one)
    expect_true(all(is.finite(high)))
})


## SSbiexponential() ================================================
test_that("SSbiexponential() converges on known parameters", {
    set.seed(1)
    t <- 0:120
    x <- biexponential(t, A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B, tau1, tau2, prop), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "tau1", "tau2", "prop"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], 50, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B"]], 80, tolerance = 2, scale = 1))
})

test_that("SSbiexponential() predict() returns correct length", {
    set.seed(2)
    t <- 0:120
    x <- biexponential(t, 50, 80, 5, 40, 0.6) + rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(x ~ SSbiexponential(t, A, B, tau1, tau2, prop), data = data)

    expect_length(predict(model, data), nrow(data))
})

test_that("SSbiexponential() fixes A at a constant", {
    set.seed(3)
    t <- 0:120
    x <- biexponential(t, A = 0, B = 80, tau1 = 5, tau2 = 40, prop = 0.6) +
        rnorm(length(t), 0, 0.5)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSbiexponential(t, A = 0, B, tau1, tau2, prop), data = data
    )

    expect_named(coef(model), c("B", "tau1", "tau2", "prop"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 0)
})


## analyse_biexponential() ==========================================

## helper: create biexponential test data with known parameters
create_biexp_data <- function(
    A = 50,
    B = 80,
    tau1 = 5,
    tau2 = 40,
    prop = 0.6,
    n = 120,
    sample_rate = 1,
    noise_sd = 0.5,
    channels = "smo2",
    seed = 42
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    x <- biexponential(t, A, B, tau1, tau2, prop) + rnorm(n, 0, noise_sd)

    df <- setNames(data.frame(t, x), c("time", channels[1]))
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- biexponential(t, A + 5, B + 5, tau1, tau2, prop) +
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
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "time_channel",
        "A", "B", "tau1", "tau2", "prop", "MRT",
        "tau1_fitted", "tau2_fitted", "MRT_fitted"
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
    A <- 50
    B <- 80
    result <- analyse_biexponential(
        create_biexp_data(A = A, B = B, noise_sd = 0.3),
        nirs_channels = "smo2",
        verbose = FALSE
    )

    expect_equal(result$A, A, tolerance = 2)
    expect_equal(result$B, B, tolerance = 2)
    ## fast-first convention: tau1 <= tau2
    expect_true(result$tau1 <= result$tau2)
    expect_true(result$prop >= 0 && result$prop <= 1)
    ## MRT is the amplitude-weighted mean of the two time constants
    expect_equal(
        result$MRT,
        result$prop * result$tau1 + (1 - result$prop) * result$tau2,
        tolerance = 1e-6
    )
    ## good fit
    expect_true(attr(result, "diagnostics")$r2 > 0.9)
})

test_that("analyse_biexponential() uses start_time correctly", {
    start_time <- 12
    data <- create_biexp_data(noise_sd = 0.3)
    data$time <- data$time + start_time

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        verbose = FALSE
    )

    expect_equal(result$A, 50, tolerance = 2)
    expect_equal(result$B, 80, tolerance = 2)
    expect_true(result$tau1 <= result$tau2)
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
    ## too few observations for a 5-param model
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
    expect_true(is.na(result$MRT))
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


## direction ========================================================

test_that("analyse_biexponential() direction = 'positive' rejects falling fit", {
    data <- create_biexp_data(A = 80, B = 50, noise_sd = 0.3)

    expect_warning(
        result <- analyse_biexponential(
            data,
            nirs_channels = "smo2",
            direction = "positive",
            verbose = TRUE
        ),
        "satisfy"
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$B))
    expect_true(is.na(result$tau1))
})

test_that("analyse_biexponential() per-channel direction overrides", {
    data <- create_biexp_data(
        A = 80, B = 50, noise_sd = 0.3, channels = c("ch1", "ch2")
    )

    result <- analyse_biexponential(
        data,
        nirs_channels = c("ch1", "ch2"),
        direction = list("auto", ch2 = "positive"),
        verbose = FALSE
    )

    expect_false(is.na(result$A[result$nirs_channels == "ch1"]))
    expect_true(is.na(result$A[result$nirs_channels == "ch2"]))

    ca <- attr(result, "channel_args")
    expect_equal(ca$direction, c("negative", "positive"))
})


## fixed parameters =================================================

test_that("analyse_biexponential() fix holds parameters constant", {
    data <- create_biexp_data(A = 0, B = 30, noise_sd = 0.3)

    result <- analyse_biexponential(
        data,
        nirs_channels = "smo2",
        fix = list(A = 0),
        verbose = FALSE
    )

    expect_equal(result$A, 0)
    expect_equal(result$B, 30, tolerance = 2)

    ## fixed A excluded from the fitted model coefficients
    expect_named(coef(attr(result, "model")$smo2), c("B", "tau1", "tau2", "prop"))
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
            fix = list(A = 50, B = 80, tau1 = 5, tau2 = 40, prop = 0.6)
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
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs_kinetics")
    expect_equal(result$method, "biexponential")
    expect_true(all(
        c("tau1", "tau2", "prop", "MRT") %in% names(result$coefficients)
    ))
})

test_that("analyse_kinetics() resolves the 'biexp' alias", {
    data <- create_biexp_data(noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexp",
        verbose = FALSE
    )

    expect_equal(result$method, "biexponential")
})

test_that("analyse_kinetics() passes fix to the biexponential method", {
    data <- create_biexp_data(A = 0, B = 30, noise_sd = 0.3)

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "biexponential",
        fix = list(A = 0),
        verbose = FALSE
    )

    expect_equal(result$coefficients$A, 0)
    expect_named(
        coef(result$model[[1L]]$smo2), c("B", "tau1", "tau2", "prop")
    )
})
