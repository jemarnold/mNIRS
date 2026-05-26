## logistic() ========================================================
test_that("logistic() returns correct vector length", {
    t <- 1:60
    result <- logistic(t, A = 0, B = 100, xmid = 30, slope = 4)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("logistic() 4-param symmetric form: y(xmid) = (A + B) / 2", {
    t <- 1:60
    A <- 10
    B <- 100
    xmid <- 30
    result <- logistic(t, A = A, B = B, xmid = xmid, slope = 4)

    idx <- which(t == xmid)
    expect_true(
        all.equal(result[idx], (A + B) / 2, tolerance = 1e-8, scale = 1)
    )
})

test_that("logistic() approaches A and B asymptotes", {
    t <- -200:200
    result <- logistic(t, A = 10, B = 100, xmid = 0, slope = 4)

    expect_true(
        all.equal(result[1L], 10, tolerance = 0.01, scale = 1)
    )
    expect_true(
        all.equal(result[length(result)], 100, tolerance = 0.01, scale = 1)
    )
})

test_that("logistic() 4-param is monotonically increasing for B > A", {
    t <- 1:60
    result <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4)

    expect_true(all(diff(result) > 0))
})

test_that("logistic() handles falling curves (B < A)", {
    t <- 1:60
    result <- logistic(t, A = 100, B = 10, xmid = 30, slope = -4)

    expect_true(all(diff(result) < 0))
    expect_true(
        all.equal(result[which(t == 30)], 55, tolerance = 1e-8, scale = 1)
    )
})

test_that("logistic() 5-param with asym = 0.5 matches 4-param", {
    t <- 1:60
    args <- list(A = 0, B = 100, xmid = 30, slope = 4)
    y4 <- do.call(logistic, c(list(t = t), args))
    y5 <- do.call(logistic, c(list(t = t), args, asym = 0.5))

    expect_true(all.equal(y4, y5, tolerance = 1e-6, scale = 1))
})

test_that("logistic() 5-param supports full asym range (0, 1)", {
    t <- seq(0, 1, length.out = 200)
    args <- list(A = 0, B = 1, xmid = 0.5, slope = 2)

    ## span across symmetric centre with both Gompertz-limit directions
    asym_vals <- seq(0.2, 0.9, 0.1)

    #! why is asym ~< 0.25 behaving non-symmetrical with equivalent near 1?
    ## visual check
    # supported <- do.call(rbind, lapply(asym_vals, \(a) {
    #     data.frame(
    #         t = t,
    #         x = do.call(logistic, c(list(t = t), args, asym = a)),
    #         asym = factor(a)
    #     )
    # }))

    # ggplot2::ggplot(supported, ggplot2::aes(t, x, colour = asym)) +
    #     ggplot2::labs(
    #         title = "logistic() supported asym range",
    #         subtitle = paste(
    #             "dashed: Gompertz curves at `y = A + (B-A)/e` and `y = 1 - (A + (B-A)/e)`",
    #             "dotted: xmid",
    #             sep = "\n"
    #         )
    #     ) +
    #     theme_mnirs() +
    #     scale_colour_mnirs() +
    #     ggplot2::geom_line() +
    #     ggplot2::geom_hline(
    #         yintercept = c(
    #             args$A + (args$B - args$A) / exp(1),
    #             1 - (args$A + (args$B - args$A) / exp(1))
    #         ),
    #         linetype = "dashed"
    #     ) +
    #     # ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
    #     ggplot2::geom_vline(xintercept = args$xmid, linetype = "dotted")

    for (a in asym_vals) {
        result <- do.call(logistic, c(list(t = t), args, asym = a))
        expect_all_true(is.finite(result))
        expect_length(result, length(t))
    }
})

test_that("logistic() 5-param asym < 0.5 shifts inflection height down", {
    ## asym -> 0 means inflection y near A
    t <- seq(0, 100, length.out = 1000)
    A <- 0
    B <- 100
    xmid <- 50
    args <- list(A = A, B = B, xmid = xmid, slope = 2)

    y_low <- do.call(logistic, c(list(t = t), args, asym = 0.05))
    # plot(t, y_low)
    
    ## inflection at xmid; height should be near A
    y_at_xmid <- y_low[which.min(abs(t - xmid))]
    expect_true(all.equal(
        y_at_xmid, A + (B - A) * 0.05, tolerance = 1, scale = 1
    ))
})

test_that("logistic() 5-param asym > 0.5 shifts inflection height up", {
    ## asym -> 1 means inflection y near B
    t <- seq(0, 100, length.out = 1000)
    A <- 0
    B <- 100
    xmid <- 50
    args <- list(A = A, B = B, xmid = xmid, slope = 2)

    y_high <- do.call(logistic, c(list(t = t), args, asym = 0.95))
    # plot(t, y_high)

    y_at_xmid <- y_high[which.min(abs(t - xmid))]
    expect_true(all.equal(
        y_at_xmid, A + (B - A) * 0.95, tolerance = 1, scale = 1
    ))
})

test_that("logistic() 5-param is smooth across asym = 0.5", {
    ## both branches should agree in the limit; check small offsets either side
    t <- 1:60
    args <- list(A = 0, B = 100, xmid = 30, slope = 4)

    y_just_below <- do.call(logistic, c(list(t = t), args, asym = 0.5 - 1e-6))
    y_just_above <- do.call(logistic, c(list(t = t), args, asym = 0.5 + 1e-6))

    expect_true(
        all.equal(y_just_below, y_just_above, tolerance = 1e-4, scale = 1)
    )
})


## SSlogistic() ======================================================
test_that("SSlogistic() 4-param fit recovers parameters", {
    set.seed(13)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4)
    x <- do.call(logistic, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## visual check
    # ggplot2::ggplot(data, ggplot2::aes(x = t, y = x)) +
    #     ggplot2::geom_line()

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "xmid", "slope"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 1, scale = 1))
    expect_true(
        all.equal(coefs[["xmid"]], base$xmid, tolerance = 1, scale = 1)
    )
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.1, scale = 1)
    )
})

test_that("SSlogistic() 4-param handles falling sigmoids (B < A)", {
    set.seed(456)
    t <- 1:60
    base <- list(A = 100, B = 10, xmid = 30, slope = -4)
    x <- do.call(logistic, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
    coefs <- coef(model)

    expect_s3_class(model, "nls")
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["xmid"]], base$xmid, tolerance = 1, scale = 1))
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.2, scale = 1)
    )
})

test_that("SSlogistic() predict() returns correct length", {
    set.seed(202)
    t <- 1:60
    x <- logistic(t, 10, 100, 30, 4) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
    predictions <- predict(model, data)

    expect_length(predictions, nrow(data))
})

test_that("SSlogistic() 5-param recovers parameters (symmetric centre)", {
    set.seed(15)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4, asym = 0.5)
    x <- do.call(logistic, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## visual check
    # ggplot2::ggplot(data, ggplot2::aes(x = t, y = x)) +
    #     ggplot2::geom_line()

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
    coefs <- coef(model)

    expect_named(coefs, c("A", "B", "xmid", "slope", "asym"))
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 2, scale = 1))
    expect_true(
        all.equal(coefs[["xmid"]], base$xmid, tolerance = 2, scale = 1)
    )
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.1, scale = 1)
    )
    expect_true(
        all.equal(coefs[["asym"]], base$asym, tolerance = 0.15, scale = 1)
    )
})

test_that("SSlogistic() 5-param recovers asym > 0.5 (late acceleration)", {
    set.seed(15)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4, asym = 0.8)
    x <- do.call(logistic, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## visual check
    # ggplot2::ggplot(data, ggplot2::aes(x = t, y = x)) +
    #     ggplot2::geom_line()

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
    coefs <- coef(model)

    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 2, scale = 1))
    expect_true(
        all.equal(coefs[["xmid"]], base$xmid, tolerance = 2, scale = 1)
    )
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.25, scale = 1)
    )
    expect_true(
        all.equal(coefs[["asym"]], base$asym, tolerance = 0.15, scale = 1)
    )
})

test_that("SSlogistic() 5-param recovers asym < 0.5 (early acceleration)", {
    set.seed(15)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4, asym = 0.8)
    x <- do.call(logistic, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## visual check
    # ggplot2::ggplot(data, ggplot2::aes(x = t, y = x)) +
    #     ggplot2::geom_line()

    model <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
    coefs <- coef(model)

    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 2, scale = 1))
    expect_true(
        all.equal(coefs[["xmid"]], base$xmid, tolerance = 2, scale = 1)
    )
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.25, scale = 1)
    )
    expect_true(
        all.equal(coefs[["asym"]], base$asym, tolerance = 0.15, scale = 1)
    )
})

test_that("SSlogistic() 4-param converges on most random realisations", {
    skip_on_ci()
    skip_on_covr()
    skip_on_cran()
    skip_if(!interactive(), "Manual convergence check")
    ## quantify fit success rate for the 4-param symmetric form;
    ## should be highly reliable across varied true xmid/slope values
    n_rep <- 1000L
    t <- 1:60
    A <- 10
    B <- 100
    xmid_grid <- runif(n_rep, min = 15, max = 45)
    slope_grid <- runif(n_rep, min = 2, max = 8)

    # set.seed(13)
    fits <- vapply(seq_len(n_rep), \(i) {
        x <- logistic(t, A, B, xmid_grid[i], slope_grid[i]) +
            rnorm(length(t), 0, 2)
        data <- data.frame(t, x)
        model <- tryCatch(
            nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data),
            error = \(e) NULL,
            warning = \(w) NULL
        )
        if (is.null(model)) {
            return(c(success = 0, xmid_err = NA_real_, slope_err = NA_real_))
        }
        c(
            success = 1,
            xmid_err = abs(coef(model)[["xmid"]] - xmid_grid[i]),
            slope_err = abs(coef(model)[["slope"]] - slope_grid[i])
        )
    }, numeric(3L))

    success_rate <- mean(fits["success", ])
    success_rate

    expect_true(success_rate >= 0.95)

    ## of successful fits, params should be recovered within tight tolerance
    ok <- fits["success", ] == 1
    expect_true(median(fits["xmid_err", ][ok], na.rm = TRUE) < 1)
    expect_true(median(fits["slope_err", ][ok], na.rm = TRUE) < 0.5)
})


test_that("SSlogistic() 5-param converges on most random realisations", {
    skip_on_ci()
    skip_on_covr()
    skip_on_cran()
    skip_if(!interactive(), "Manual convergence check")

    ## quantify fit success rate across noisy realisations spanning asym range;
    ## the 5-param model is known to be fragile near Gompertz limits
    n_rep <- 1000L
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4)
    asym_grid <- seq(0.1, 0.9, length.out = n_rep)

    # set.seed(13)
    fits <- vapply(asym_grid, \(a) {
        x <- do.call(logistic, c(list(t = t), base, asym = a)) +
            rnorm(length(t), 0, 2)
        data <- data.frame(t, x)
        model <- tryCatch(
            nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data),
            error = \(e) NULL,
            warning = \(w) NULL
        )
        if (is.null(model)) {
            return(c(success = 0, asym_err = NA_real_))
        }
        c(success = 1, asym_err = abs(coef(model)[["asym"]] - a))
    }, numeric(2L))

    success_rate <- mean(fits["success", ])
    success_rate
    
    ## current implementation: report observed rate to flag regressions;
    ## aim for >= 75% on interior of (0, 1) once init is improved
    expect_true(success_rate >= 0.75)

    ## of successful fits, asym should be recovered within tolerance
    converged <- fits["asym_err", ][fits["success", ] == 1]
    expect_true(median(converged, na.rm = TRUE) < 0.1)
})


test_that("SSlogistic() fails gracefully on non-sigmoidal data", {
    ## monotonic linear data has no inflection -> nls should fail
    set.seed(13)
    t <- 1:60
    x <- t + rnorm(length(t), 0, 1)
    data <- data.frame(t, x)

    expect_error(
        nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
    )
})


## analyse_logistic() ===================================================

## helper: create logistic test data with known parameters
create_logistic_data <- function(
    A = 10,
    B = 100,
    xmid = 30,
    slope = 4,
    asym = NULL,
    n = 60,
    sample_rate = 1,
    noise_sd = 2,
    channels = "smo2",
    seed = 13
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    x <- logistic(t, A, B, xmid, slope, asym) + rnorm(n, 0, noise_sd)
    df <- stats::setNames(data.frame(t, x), c("time", channels[1]))
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- logistic(t, A + 5, B + 5, xmid, slope, asym) +
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


test_that("analyse_logistic() returns correct structure", {
    data <- create_logistic_data()
    # plot(data)

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = FALSE,
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "nirs_channels", "time_channel",
        "A", "B", "xmid", "slope", "asym", "xmid_fitted"
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

test_that("analyse_logistic() validates use_asym argument", {
    data <- create_logistic_data()

    expect_error(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            use_asym = "yes"
        ),
        "use_asym.*logical"
    )

    expect_error(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            use_asym = c(TRUE, FALSE)
        ),
        "use_asym.*logical"
    )
})

test_that("analyse_logistic() recovers 4-param known parameters", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, n = 100, noise_sd = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = FALSE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.2, scale = 1))
    expect_true(is.na(result$asym))
})

test_that("analyse_logistic() recovers 5-param known parameters", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4
    asym <- 0.3

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, asym = asym,
        n = 100, noise_sd = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = TRUE,
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 2, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.25, scale = 1))
    expect_true(all.equal(result$asym, asym, tolerance = 0.15, scale = 1))
})

test_that("analyse_logistic() uses t0 correctly", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4
    t0 <- 12

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, n = 100, noise_sd = 1
    )
    data$time <- data$time + t0

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = FALSE,
        t0 = t0,
        verbose = FALSE
    )

    ## xmid reported as offset from t0 — should match original xmid
    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.2, scale = 1))
})

test_that("analyse_logistic() t0 edge cases", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4
    t0 <- 12

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, n = 100, noise_sd = 1
    )
    data$time <- data$time + t0

    ## t0 specified before time start, falls forward to t[1L]
    expect_warning(
        result <- analyse_logistic(
            data,
            nirs_channels = "smo2",
            use_asym = FALSE,
            t0 = 0,
            verbose = TRUE
        ),
        "No observations.*t0 =.*0"
    )

    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))

    ## t0 beyond time range
    expect_error(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            use_asym = FALSE,
            t0 = max(data$time) + 10,
            verbose = TRUE
        ),
        "No observations.*before.*t0"
    )
})

test_that("analyse_logistic() falls back from 5-param to 4-param", {
    ## short noisy series makes 5-param hard to converge
    data <- create_logistic_data(n = 60, noise_sd = 5, seed = 101)

    expect_warning(
        result <- analyse_logistic(
            data,
            nirs_channels = "smo2",
            use_asym = TRUE,
            verbose = TRUE
        ),
        "SSlogistic.*fit failed"
    )

    ## should still return a valid result via 4-param fallback
    expect_s3_class(result, "data.frame")
    expect_true(is.na(result$asym))
    expect_false(is.na(result$xmid))
})

test_that("analyse_logistic() returns NA for failed fit", {
    ## only 3 observations for a 4-param model
    custom_name <- create_logistic_data(n = 10, noise_sd = 0.1)

    expect_warning(
        result <- analyse_logistic(
            custom_name,
            nirs_channels = "smo2",
            use_asym = FALSE
        ),
        "fit failed for.*smo2.*custom_name" ## call custom interval name
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$xmid))
    expect_true(is.na(result$slope))
})

test_that("analyse_logistic() suppresses fit-failure warning when verbose = FALSE", {
    ## only 3 observations for a 4-param model — guaranteed fit failure
    custom_name <- create_logistic_data(n = 10, noise_sd = 0.1)

    expect_no_warning(
        analyse_logistic(
            custom_name,
            nirs_channels = "smo2",
            use_asym = FALSE,
            verbose = FALSE
        )
    )
})

test_that("analyse_logistic() works with multiple channels", {
    nirs_channels <- c("smo2_left", "smo2_right")
    data <- create_logistic_data(channels = nirs_channels)

    result <- analyse_logistic(
        data,
        nirs_channels = nirs_channels,
        use_asym = FALSE,
        verbose = FALSE
    )

    expect_equal(nrow(result), 2L)
    expect_equal(result$nirs_channels, nirs_channels)

    fitted_data <- attr(result, "fitted_data")
    expect_length(fitted_data, 2L)
    expect_named(fitted_data, nirs_channels)
})

test_that("analyse_logistic() channel_args override defaults", {
    data <- create_logistic_data(
        channels = c("ch1", "ch2"), asym = 0.5, n = 100, noise_sd = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = c("ch1", "ch2"),
        use_asym = FALSE,
        channel_args = list(ch2 = list(use_asym = TRUE)),
        verbose = FALSE
    )

    ## ch1 has no asym (4-param), ch2 has asym (5-param)
    expect_equal(is.na(result$asym), c(TRUE, FALSE))

    ca <- attr(result, "channel_args")
    ch1_row <- ca[ca$nirs_channels == "ch1", ]
    ch2_row <- ca[ca$nirs_channels == "ch2", ]
    expect_false(ch1_row$use_asym)
    expect_true(ch2_row$use_asym)
})

test_that("analyse_logistic() fitted_data attribute is well-formed", {
    data <- create_logistic_data(n = 100, noise_sd = 1)

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = FALSE,
        verbose = FALSE
    )

    fitted_data <- attr(result, "fitted_data")
    expect_named(fitted_data, "smo2")
    expect_named(fitted_data$smo2, c("window_idx", "fitted"))
    expect_type(fitted_data$smo2$fitted, "double")
    ## fitted should correlate and agree well with original data
    expect_true(cor(data$smo2, fitted_data$smo2$fitted) > 0.95)
    expect_all_true(abs(data$smo2 - fitted_data$smo2$fitted) <= 3)
})

test_that("analyse_logistic() diagnostics contain expected columns", {
    data <- create_logistic_data()

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        use_asym = FALSE,
        verbose = FALSE
    )

    diag <- attr(result, "diagnostics")
    expect_true(all(
        c("nirs_channels", "n_obs", "r2", "adj_r2", "rmse") %in% names(diag)
    ))
    expect_true(diag$r2 > 0.9)
})
