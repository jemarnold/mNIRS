## monoexponential() =================================================
test_that("monoexponential() returns correct vector length", {
    t <- 1:60
    result <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("monoexponential() returns baseline A before TD", {
    t <- 0:30
    A <- 10
    TD <- 15
    result <- monoexponential(t, A = A, B = 100, tau = 8, TD = TD)

    expect_true(all(result[t < TD] == A))
})

test_that("monoexponential() approaches asymptote B", {
    t <- 1:200
    B <- 100
    result <- monoexponential(t, A = 10, B = B, tau = 8, TD = 15)

    # At t >> TD + tau, should approach B
    expect_true(
        all.equal(result[length(result)], B, tolerance = 0.01, scale = 1)
    )
})

test_that("monoexponential() handles rising curves (B > A)", {
    t <- 1:60
    result <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15)

    expect_true(result[20] > result[15])
    expect_true(result[30] > result[20])
})

test_that("monoexponential() handles falling curves (B < A)", {
    t <- 1:60
    result <- monoexponential(t, A = 100, B = 10, tau = 8, TD = 15)

    expect_true(result[20] < result[15])
    expect_true(result[30] < result[20])
})

test_that("monoexponential() tau determines rate correctly", {
    t <- 1:100
    TD <- 15
    tau <- 10

    result <- monoexponential(t, A = 0, B = 100, tau = tau, TD = TD)

    # At t = TD + tau, should be ~63.2% of amplitude
    idx <- which(t == TD + tau)
    expect_true(all.equal(result[idx], 63.2, tolerance = 1, scale = 1))
})

test_that("monoexponential() handles zero and negative TD", {
    t <- 0:60
    result <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 0)
    
    expect_true(all(!is.na(result)))
    expect_true(result[1] == 10)
    expect_true(result[10] > 10)
    
    t <- -10:60
    result <- monoexponential(t, A = 10, B = 100, tau = 8, TD = -5)
    expect_true(all(!is.na(result)))
    expect_all_equal(result[1:6], 10)
    expect_true(result[7] > 10)
})


## SSmonoexp() ========================================================
test_that("SSmonoexp() converges on known parameters", {
    set.seed(13)
    t <- 1:60
    A_true <- 10
    B_true <- 100
    TD_true <- 15
    tau_true <- 8

    x <- monoexponential(t, A_true, B_true, tau_true, TD_true) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    # ggplot2::ggplot(data, ggplot2::aes(t, x)) +
    #     theme_mnirs() +
    #     ggplot2::geom_point() #+
        # ggplot2::geom_line(ggplot2::aes(y = y))

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "tau", "TD"))

    coefs <- coef(model)
    expect_true(
        all.equal(coefs[["A"]], A_true, tolerance = A_true * 0.1, scale = 1)
    )
    expect_true(
        all.equal(coefs[["B"]], B_true, tolerance = B_true * 0.1, scale = 1)
    )
    expect_true(all.equal(
            coefs[["tau"]], tau_true, tolerance = tau_true * 0.1, scale = 1
    ))
    expect_true(
        all.equal(coefs[["TD"]], TD_true, tolerance = TD_true * 0.1, scale = 1)
    )
})

test_that("SSmonoexp() handles falling exponentials", {
    set.seed(456)
    t <- 1:60
    x <- monoexponential(t, A = 100, B = 10, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    coefs <- coef(model)

    expect_true(coefs[["A"]] > coefs[["B"]])
    expect_s3_class(model, "nls")
})

test_that("SSmonoexp() predict() returns correct length", {
    set.seed(202)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    predictions <- predict(model, data)

    expect_length(predictions, nrow(data))
})

## SSmonoexp() ========================================================
test_that("SSmonoexp() converges on known parameters", {
    set.seed(13)
    t <- 1:60-1
    A_true <- 10
    B_true <- 100
    TD_true <- 0
    tau_true <- 8

    x <- monoexponential(t, A_true, B_true, tau_true, TD_true) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)
    
    model <- nls(x ~ SSmonoexp(t, A, B, tau), data = data)
    
    ## visual check
    # y <- fitted(model)
    # ggplot2::ggplot(data, ggplot2::aes(t, x)) +
    #     theme_mnirs() +
    #     ggplot2::geom_point() +
    # ggplot2::geom_line(ggplot2::aes(y = y))

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "tau"))

    coefs <- coef(model)
    expect_true(
        all.equal(coefs[["A"]], A_true, tolerance = 3, scale = 1)
    )
    expect_true(
        all.equal(coefs[["B"]], B_true, tolerance = 1, scale = 1)
    )
    expect_true(all.equal(
            coefs[["tau"]], tau_true, tolerance = 1, scale = 1
        ))
    expect_disjoint(names(coefs), "TD")
})

test_that("SSmonoexp() handles falling exponentials", {
    set.seed(456)
    t <- 1:60
    x <- monoexponential(t, A = 100, B = 10, tau = 8, TD = NULL) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau), data = data)
    coefs <- coef(model)

    expect_true(coefs[["A"]] > coefs[["B"]])
    expect_s3_class(model, "nls")
    expect_disjoint(names(coefs), "TD")
})

test_that("SSmonoexp() handles data with TD near zero", {
    set.seed(101)
    t <- 1:60
    x <- monoexponential(t, 10, 100, 8, 1) + rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    # ggplot2::ggplot(data, ggplot2::aes(t, x)) +
    #     theme_mnirs() +
    #     ggplot2::geom_point()

    ## TODO SSmonoexp4() fails for this test. Would a better initialisation succeed?
    expect_no_error(
        model <- nls(x ~ SSmonoexp(t, A, B, tau), data = data)
    )

    expect_s3_class(model, "nls")
    expect_true(all.equal(coef(model)[["tau"]], 8, tolerance = 1, scale = 1))
    expect_disjoint(names(coef(model)), "TD")
})

## OxCap modelling ===================================================
test_that("SSmonoexp() handles OxCap with few data points same as SSasymp", {
    # fmt: skip
    df <- data.frame(
        time = c(
            0, 15.1, 30.8, 45.8, 61, 76, 91.2, 106.3, 121.1, 135.9, 151.1, 
            165.8, 181, 196.1, 211.1, 225.8
        ),
        slope = c(
            2.5762, 1.3261, 0.8707, 0.733, 0.4198, 0.3614, 0.2806, 0.3069, 
            0.3728, 0.4752, 0.5505, 0.5538, 0.5384, 0.5036, 0.4911, 0.5088
        )
    )

    expect_no_error(
        model <- nls(
            slope ~ SSmonoexp(time, A, B, tau),
            data = df,
        )
    )

    model_asym <- nls(
        slope ~ SSasymp(time, Asym, R0, lrc),
        data = df,
    )
    asym_tau <- exp(-coef(model_asym)[["lrc"]])

    expect_true(all.equal(
        coef(model)[["A"]],
        coef(model_asym)[["R0"]],
        tolerance = 0.1,
        scale = 1
    ))
    expect_true(all.equal(
        coef(model)[["B"]],
        coef(model_asym)[["Asym"]],
        tolerance = 0.1,
        scale = 1
    ))
    expect_true(all.equal(
        coef(model)[["tau"]], asym_tau, tolerance = 0.1, scale = 1
    ))
})

test_that("SSmonoexp() handles OxCap with few data points better than SSasymp", {
    # fmt: skip
    df <- data.frame(
        time = c(0, 15.9, 30.9, 46, 60.9, 76, 91.5, 106.3, 121.2, 136.2, 151, 
            166.3, 181, 195.8, 211, 225.6),
        slope = c(
            2.5868, 1.1626, 0.6287, 0.3786, 0.162, 0.2219, 0.173, 0.1864, 
            0.3341, 0.2669, 0.3361, 0.4534, 0.3756, 0.4536, 0.3664, 0.4137
        )
      )

    expect_no_error(
        model <- nls(
            slope ~ SSmonoexp(time, A, B, tau),
            data = df,
        )
    )

    expect_error(
        model_asym <- nls(
            slope ~ SSasymp(time, Asym, R0, lrc),
            data = df,
        ),
        "singular gradient"
    )
})


## analyse_monoexponential() ===========================================

## helper: create monoexponential test data with known parameters
create_monoexp_data <- function(
    A = 50,
    B = 80,
    tau = 25,
    TD = 0,
    n = 60,
    sample_rate = 1,
    noise_sd = 0.5,
    channels = "smo2",
    seed = 42
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    x <- monoexponential(t, A, B, tau, TD) + rnorm(n, 0, noise_sd)

    df <- setNames(
        data.frame(t, x),
        c("time", channels[1])
    )
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- monoexponential(
                t,
                A + 5,
                B + 5,
                tau,
                TD
            ) +
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


test_that("analyse_monoexponential() returns correct structure", {
    data <- create_monoexp_data()

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = FALSE,
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "nirs_channels", "time_channel", 
        "A", "B", "tau", "k", "TD", "MRT", "HRT",
        "tau_fitted", "MRT_fitted", "HRT_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    expect_type(attr(result, "model"), "list")
    expect_true(inherits(attr(result, "model")$smo2, "nls"))
    expect_type(attr(result, "fitted_data"), "list")
    ## predicted data frame for each `nirs_channels`
    expect_s3_class(attr(result, "fitted_data")$smo2, "data.frame")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_s3_class(attr(result, "diagnostics"), "data.frame")
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_s3_class(attr(result, "channel_args"), "data.frame")
    expect_equal(nrow(attr(result, "channel_args")), 1L)
})

test_that("analyse_monoexponential() validates use_time_delay argument", {
    data <- create_monoexp_data()

    expect_error(
        analyse_monoexponential(
            data,
            nirs_channels = "smo2",
            use_time_delay = "yes"
        ),
        "use_time_delay.*logical"
    )

    expect_error(
        analyse_monoexponential(
            data,
            nirs_channels = "smo2",
            use_time_delay = c(TRUE, FALSE)
        ),
        "use_time_delay.*logical"
    )
})

test_that("analyse_monoexponential() recovers 3-param known parameters", {
    A <- 50
    B <- 80
    tau <- 25

    data <- create_monoexp_data(
        A = A, B = B, tau = tau, n = 100, noise_sd = 0.3
    )

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = FALSE,
        verbose = FALSE
    )

    expect_equal(result$A, A, tolerance = 1)
    expect_equal(result$B, B, tolerance = 1)
    expect_equal(result$tau, tau, tolerance = 1)
    expect_true(is.na(result$TD))
    expect_equal(result$k, 1 / tau, tolerance = 1)
    expect_equal(result$HRT, tau * log(2), tolerance = 1)
})

test_that("analyse_monoexponential() recovers 4-param known parameters", {
    A <- 50
    B <- 80
    tau <- 25
    TD <- 10

    data <- create_monoexp_data(
        A = A, B = B, tau = tau, TD = TD, n = 100, noise_sd = 0.3
    )

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = TRUE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$tau, tau, tolerance = 1, scale = 1))
    expect_true(all.equal(result$TD, TD, tolerance = 1, scale = 1))
    expect_true(all.equal(result$k, 1 / tau, tolerance = 1, scale = 1))
    expect_true(
        all.equal(result$HRT, tau * log(2) + TD, tolerance = 1, scale = 1)
    )
})

test_that("analyse_monoexponential() uses t0 correctly", {
    A <- 50
    B <- 80
    tau <- 25
    TD <- 10
    t0 <- 12

    data <- create_monoexp_data(
        A = A, B = B, tau = tau, TD = TD, n = 100, noise_sd = 0.3
    )
    data$time <- data$time + t0

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = TRUE,
        t0 = t0,
        verbose = FALSE
    )

    ## visual check
    # plot(data) +
    #     # ggplot2::coord_cartesian(xlim = c(0, NA)) +
    #     # ggplot2::geom_vline(xintercept = t0) +
    #     ggplot2::geom_line(
    #         ggplot2::aes(y = attributes(result)$fitted$smo2$fitted)
    #     )

    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$tau, tau, tolerance = 1, scale = 1))
    expect_true(all.equal(result$TD, TD, tolerance = 1, scale = 1))
    expect_true(all.equal(result$MRT, TD + tau, tolerance = 1, scale = 1))
    expect_true(all.equal(result$k, 1 / tau, tolerance = 1, scale = 1))
    expect_true(
        all.equal(result$HRT, tau * log(2) + TD, tolerance = 1, scale = 1)
    )

    data <- create_monoexp_data(
        A = A, B = B, tau = tau, TD = TD, n = 100, noise_sd = 0.3
    )
    analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = TRUE,
        t0 = 5,
        verbose = FALSE
    )
})

test_that("analyse_monoexponential() t0 edge cases", {
    A <- 50
    B <- 80
    tau <- 25
    TD <- 10
    t0 <- 12

    data <- create_monoexp_data(
        A = A, B = B, tau = tau, TD = TD, n = 100, noise_sd = 0.3
    )
    data$time <- data$time + t0

    ## t0 specified before time start, falls forward to t[1L]
    expect_warning(
        result <- analyse_monoexponential(
            data,
            nirs_channels = "smo2",
            use_time_delay = TRUE,
            t0 = 0,
            verbose = TRUE
        ),
        "No observations.*t0 =.*0"
    )

    expect_true(all.equal(result$TD, TD, tolerance = 1, scale = 1))
    expect_true(all.equal(result$MRT, TD + tau, tolerance = 1, scale = 1))
    expect_true(
        all.equal(result$HRT, tau * log(2) + TD, tolerance = 1, scale = 1)
    )

    ## t0 beyond time range
    expect_error(
        analyse_monoexponential(
            data,
            nirs_channels = "smo2",
            use_time_delay = TRUE,
            t0 = max(data$time) + 10,
            verbose = TRUE
        ),
        "No observations.*before.*t0"
    )
})

test_that("analyse_monoexponential() falls back from 4-param to 3-param", {
    A = 50
    B = 80
    tau = 25
    TD = 0

    ## short series with small TD makes 4-param hard to converge
    data <- create_monoexp_data(
        A = A, B = B, tau = tau, TD = TD, n = 10, noise_sd = 2, seed = 101
    )

    expect_warning(
        result <- analyse_monoexponential(
            data,
            nirs_channels = "smo2",
            use_time_delay = TRUE,
            verbose = TRUE
        ),
        "SSmonoexp.*fit failed"
    )

    ## should still return a valid result via 3-param fallback
    ## albeit with poor fit
    expect_s3_class(result, "data.frame")
    expect_true(is.na(result$TD))
    expect_false(is.na(result$tau))
})


test_that("analyse_monoexponential() returns NA for failed fit", {
    ## only 3 observations for a 3-param model
    custom_name <- create_monoexp_data(n = 3, noise_sd = 0.1)

    expect_warning(
        result <- analyse_monoexponential(
            custom_name,
            nirs_channels = "smo2",
            use_time_delay = FALSE
        ),
        "fit failed for.*smo2.*custom_name" ## call custom interval name
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$tau))
    expect_true(is.na(result$k))
})

test_that("analyse_monoexponential() suppresses fit-failure warning when verbose = FALSE", {
    ## only 3 observations for a 3-param model — guaranteed fit failure
    custom_name <- create_monoexp_data(n = 3, noise_sd = 0.1)

    expect_no_warning(
        analyse_monoexponential(
            custom_name,
            nirs_channels = "smo2",
            use_time_delay = FALSE,
            verbose = FALSE
        )
    )
})

test_that("analyse_monoexponential() works with multiple channels", {
    nirs_channels <- c("smo2_left", "smo2_right")
    data <- create_monoexp_data(channels = nirs_channels)

    result <- analyse_monoexponential(
        data,
        nirs_channels = nirs_channels,
        use_time_delay = FALSE
    )

    expect_equal(nrow(result), 2L)
    expect_equal(result$nirs_channels, nirs_channels)

    fitted_data <- attr(result, "fitted_data")
    expect_length(fitted_data, 2L)
    expect_named(fitted_data, nirs_channels)
})

test_that("analyse_monoexponential() channel_args override defaults", {
    data <- create_monoexp_data(channels = c("ch1", "ch2"), TD = 15)

    result <- analyse_monoexponential(
        data,
        nirs_channels = c("ch1", "ch2"),
        use_time_delay = FALSE,
        channel_args = list(ch2 = list(use_time_delay = TRUE))
    )

    ## ch1 has no TD (3-param), ch2 has TD (4-param)
    expect_equal(is.na(result$TD), c(TRUE, FALSE))

    ca <- attr(result, "channel_args")
    ch1_row <- ca[ca$nirs_channels == "ch1", ]
    ch2_row <- ca[ca$nirs_channels == "ch2", ]
    expect_false(ch1_row$use_time_delay)
    expect_true(ch2_row$use_time_delay)
})

test_that("analyse_monoexponential() fitted_data attribute is well-formed", {
    data <- create_monoexp_data(n = 100, tau = 25, TD = 15, noise_sd = 0.3)

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = TRUE
    )

    fitted_data <- attr(result, "fitted_data")
    expect_named(fitted_data, "smo2")
    expect_named(fitted_data$smo2, c("window_idx", "fitted"))
    expect_type(fitted_data$smo2$fitted, "double")
    ## fitted should correlate and agree well with original data
    expect_true(cor(data$smo2, fitted_data$smo2$fitted) > 0.9)
    expect_all_true(abs(data$smo2 - fitted_data$smo2$fitted) <= 3)

    ## visual check
    # library(ggplot2)
    # plot(data) +
    #     geom_line(aes(y = fitted_data$smo2$fitted))
})

test_that("analyse_monoexponential() diagnostics contain expected columns", {
    data <- create_monoexp_data()

    result <- analyse_monoexponential(
        data,
        nirs_channels = "smo2",
        use_time_delay = FALSE
    )

    diag <- attr(result, "diagnostics")
    expect_true(all(
        c("nirs_channels", "n_obs", "r2", "adj_r2", "rmse") %in% names(diag)
    ))
    expect_true(diag$r2 > 0.9)
})

## integration tests =============================================

test_that("extract model coefs", {
    set.seed(1111)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    tau <- coef(model)[["tau"]]
    expect_true(all.equal(tau, 8, tolerance = 1, scale = 1))

})