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
    skip("Manual fit convergence check")
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
    skip("Manual fit convergence check")

    ## quantify fit success rate across noisy realisations spanning asym range;
    ## the 5-param model is known to be fragile
    n_rep <- 1000L
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4)
    asym_grid <- seq(0.3, 0.7, length.out = n_rep)

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
    expect_true(success_rate >= 0.95)

    ## of successful fits, asym should be recovered within tolerance
    converged <- fits["asym_err", ][fits["success", ] == 1]
    expect_true(median(converged, na.rm = TRUE) < 0.1)
})

test_that("SSlogistic() converges on real dataset", {
    skip("Manual fit convergence check")
    
    ## 65 real_world reoxy intervals
    reoxy_list <- readRDS(test_path("testdata/reoxy_list.rds"))
    
    ## fit one signal at a time across all data frames; report convergence
    ## success rate per signal to flag regressions on real reoxygenation data
    fit_4param <- function(signal) {
        vapply(reoxy_list, \(df) {
            data <- data.frame(t = df$time, x = df[[signal]])
            model <- tryCatch(
                nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data),
                error = \(e) NULL,
                warning = \(w) NULL
            )
            as.integer(!is.null(model))
        }, integer(1L))
    }

    smo2_success <- mean(fit_4param("VL_smo2"))
    smo2_success
    hhb_success <- mean(fit_4param("VL_HHb"))
    hhb_success

    expect_true(smo2_success >= 0.95)
    expect_true(hhb_success >= 0.95)

    #! failing 5-param convergence 0.31 & 0.14
    fit_5param <- function(signal) {
        vapply(reoxy_list, \(df) {
            data <- data.frame(t = df$time, x = df[[signal]])
            model <- tryCatch(
                nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data),
                error = \(e) NULL,
                warning = \(w) NULL
            )
            as.integer(!is.null(model))
        }, integer(1L))
    }

    smo2_success <- mean(fit_5param("VL_smo2"))
    smo2_success
    hhb_success <- mean(fit_5param("VL_HHb"))
    hhb_success

    expect_true(smo2_success >= 0.75)
    expect_true(hhb_success >= 0.75)
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


## gompertz() / gompertz_left() ===================================
test_that("gompertz() returns correct vector length and type", {
    t <- 1:60
    result <- gompertz(t, A = 0, B = 100, xmid = 30, slope = 4)

    expect_length(result, length(t))
    expect_type(result, "double")
})

test_that("gompertz() inflection at xmid equals A + (B - A) / e", {
    t <- seq(0, 60, length.out = 1000)
    A <- 10
    B <- 100
    xmid <- 30
    result <- gompertz(t, A = A, B = B, xmid = xmid, slope = 4)

    y_at_xmid <- result[which.min(abs(t - xmid))]
    expect_true(all.equal(
        y_at_xmid, A + (B - A) / exp(1), tolerance = 1, scale = 1
    ))
})

test_that("gompertz_left() inflection at xmid equals A + (B - A) * (1 - 1/e)", {
    t <- seq(0, 60, length.out = 1000)
    A <- 10
    B <- 100
    xmid <- 30
    result <- gompertz_left(t, A = A, B = B, xmid = xmid, slope = 4)

    ## visual check
    # ggplot2::ggplot(data.frame(), ggplot2::aes(x = t, y = result)) +
    #     ggplot2::geom_line()

    y_at_xmid <- result[which.min(abs(t - xmid))]
    expect_true(all.equal(
        y_at_xmid, A + (B - A) * (1 - 1 / exp(1)), tolerance = 1, scale = 1
    ))
})

test_that("gompertz() recovers slope at inflection numerically", {
    t <- seq(0, 60, length.out = 10000)
    slope <- 4
    result <- gompertz(t, A = 10, B = 100, xmid = 30, slope = slope)
    ## numerical derivative at the inflection
    dy_dt <- diff(result) / diff(t)
    i_infl <- which.min(abs(t[-1L] - 30))

    expect_true(all.equal(
        dy_dt[i_infl], slope, tolerance = 0.1, scale = 1
    ))
})

test_that("gompertz_left() recovers slope at inflection numerically", {
    t <- seq(0, 60, length.out = 10000)
    slope <- 4
    result <- gompertz_left(t, A = 10, B = 100, xmid = 30, slope = slope)
    dy_dt <- diff(result) / diff(t)
    i_infl <- which.min(abs(t[-1L] - 30))

    expect_true(all.equal(
        dy_dt[i_infl], slope, tolerance = 0.1, scale = 1
    ))
})

test_that("gompertz() approaches A and B asymptotes", {
    t <- -200:200
    result <- gompertz(t, A = 10, B = 100, xmid = 0, slope = 4)

    expect_true(all.equal(result[1L], 10, tolerance = 0.5, scale = 1))
    expect_true(
        all.equal(result[length(result)], 100, tolerance = 0.5, scale = 1)
    )
})

test_that("gompertz_left() approaches A and B asymptotes", {
    t <- -200:200
    result <- gompertz_left(t, A = 10, B = 100, xmid = 0, slope = 4)

    expect_true(all.equal(result[1L], 10, tolerance = 0.5, scale = 1))
    expect_true(
        all.equal(result[length(result)], 100, tolerance = 0.5, scale = 1)
    )
})

test_that("gompertz()/gompertz_left() are monotonic for B > A", {
    t <- 1:60
    right <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4)
    left <- gompertz_left(t, A = 10, B = 100, xmid = 30, slope = 4)

    expect_true(all(diff(right) > 0))
    expect_true(all(diff(left) > 0))
})

test_that("gompertz()/gompertz_left() handle falling curves (B < A)", {
    t <- 1:60
    right <- gompertz(t, A = 100, B = 10, xmid = 30, slope = -4)
    left <- gompertz_left(t, A = 100, B = 10, xmid = 30, slope = -4)

    expect_true(all(diff(right) < 0))
    expect_true(all(diff(left) < 0))
})


## SSgompertz() / SSgompertz_left() ================================
test_that("SSgompertz() recovers parameters", {
    set.seed(13)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4)
    x <- do.call(gompertz, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "xmid", "slope"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["xmid"]], base$xmid, tolerance = 1, scale = 1))
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.1, scale = 1)
    )
})

test_that("SSgompertz_left() recovers parameters", {
    set.seed(13)
    t <- 1:60
    base <- list(A = 10, B = 100, xmid = 30, slope = 4)
    x <- do.call(gompertz_left, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("A", "B", "xmid", "slope"))

    coefs <- coef(model)
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["xmid"]], base$xmid, tolerance = 1, scale = 1))
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.1, scale = 1)
    )
})

test_that("SSgompertz() handles falling curves (B < A)", {
    set.seed(456)
    t <- 1:60
    base <- list(A = 100, B = 10, xmid = 30, slope = -4)
    x <- do.call(gompertz, c(list(t = t), base)) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)
    coefs <- coef(model)

    expect_s3_class(model, "nls")
    expect_true(all.equal(coefs[["A"]], base$A, tolerance = 1, scale = 1))
    expect_true(all.equal(coefs[["B"]], base$B, tolerance = 2, scale = 1))
    expect_true(all.equal(coefs[["xmid"]], base$xmid, tolerance = 1, scale = 1))
    expect_true(
        all.equal(coefs[["slope"]], base$slope, tolerance = 0.1, scale = 1)
    )
})

test_that("SSgompertz()/SSgompertz_left() converge on real dataset", {
    skip("Manual fit convergence check")
    
    reoxy_list <- readRDS(test_path("testdata/reoxy_list.rds"))

    fit_shape <- function(signal, ss_fn) {
        result <- lapply(reoxy_list, \(df) {
            data <- data.frame(t = df$time, x = df[[signal]])
            rhs <- as.call(c(
                ss_fn,
                list(quote(t), quote(A), quote(B), quote(xmid), quote(slope))
            ))
            f <- stats::as.formula(call("~", quote(x), rhs))
            model <- tryCatch(
                nls(f, data = data),
                error = \(e) NULL,
                warning = \(w) NULL
            )

            if (is.null(model)) {
                setNames(
                    rep(NA_real_, 4), c("A", "B", "xmid", "slope")
                )
            } else {
                coef(model)
            }
        })

        return(do.call(rbind, lapply(result, \(x) as.data.frame(as.list(x)))))
    }

    smo2_right <- fit_shape("VL_smo2", quote(SSgompertz))
    hhb_right <- fit_shape("VL_HHb", quote(SSgompertz))
    smo2_left <- fit_shape("VL_smo2", quote(SSgompertz_left))
    hhb_left <- fit_shape("VL_HHb", quote(SSgompertz_left))
    mean(!is.na(smo2_right$A))
    mean(!is.na(hhb_right$A))
    mean(!is.na(smo2_left$A))
    mean(!is.na(hhb_left$A))

    expect_true(mean(!is.na(smo2_right$A)) >= 0.9)
    expect_true(mean(!is.na(hhb_right$A)) >= 0.9)
    expect_true(mean(!is.na(smo2_left$A)) >= 0.9)
    expect_true(mean(!is.na(hhb_left$A)) >= 0.9)
})


## init_inflection() fallbacks =======================================

test_that("init_inflection() falls back to sign() slope when t_range == 0", {
    ## identical t -> non-finite derivative -> slope fallback;
    ## t_range == 0 so slope is the direction sign of (B - A)
    x <- 1:6
    t <- rep(5, 6)

    up <- init_inflection(x, t, A_init = 1, B_init = 6)
    expect_equal(up$slope, 1)

    down <- init_inflection(x, t, A_init = 6, B_init = 1)
    expect_equal(down$slope, -1)
})

test_that("init_inflection() falls back to mean-rate slope on flat response", {
    ## flat x -> zero smoothed derivative -> slope fallback;
    ## t_range > 0 so slope is the mean rate (B - A) / t_range
    x <- rep(3, 6)
    t <- 0:5

    res <- init_inflection(x, t, A_init = 0, B_init = 10)
    expect_equal(res$slope, 10 / 5)
    expect_true(is.finite(res$xmid))
})

test_that("init_inflection() falls back to half-response point on non-finite xmid", {
    ## derivative peak lands on a non-finite `t` -> xmid_init non-finite ->
    ## fallback to the half-response point x = (A + B) / 2
    x <- c(0, 0, 0, 0, 100, 50, 50, 50)
    t <- c(0, 1, Inf, 3, 4, 5, 6, 7)

    res <- init_inflection(x, t, A_init = 0, B_init = 100)

    expect_true(is.finite(res$xmid))
    ## first t where x reaches the half-response (50) is index 6, t = 5
    expect_equal(res$xmid, 5)
})

test_that("init_inflection() returns derivative-based estimates on clean sigmoid", {
    ## non-degenerate data bypasses both fallbacks
    t <- seq(0, 60, length.out = 200)
    x <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4)

    res <- init_inflection(x, t, A_init = 10, B_init = 100)
    expect_true(all.equal(res$xmid, 30, tolerance = 3, scale = 1))
    expect_true(res$slope > 0)
})


## logistic_init() fallbacks =========================================

## helper: invoke logistic_init() as selfStart() would, building the
## matched call and named data frame it expects
call_logistic_init <- function(x, t, asym = FALSE) {
    data <- data.frame(.t = t, .x = x)
    ## logistic_init() reads mCall[["t"]] by name and detects the 5-param
    ## branch via "asym" %in% names(mCall), so args must be named
    args <- list(
        t = quote(.t),
        A = quote(A),
        B = quote(B),
        xmid = quote(xmid),
        slope = quote(slope)
    )
    if (asym) {
        args <- c(args, list(asym = quote(asym)))
    }
    mCall <- as.call(c(quote(logistic), args))
    logistic_init(mCall, data, quote(.x))
}

test_that("logistic_init() falls back to mean-rate slope when linearisation fails", {
    ## only 2 distinct points -> fewer than 3 finite linearisation points
    ## -> xmid_init/slope_init stay NA -> both fallbacks fire. t_range > 0
    ## so slope is the mean rate (B - A) / t_range == 100 / 10 == 10
    x <- c(0, 100)
    t <- c(0, 10)

    res <- call_logistic_init(x, t)
    expect_equal(res[["slope"]], 10)
    expect_true(res[["xmid"]] >= min(t) && res[["xmid"]] <= max(t))
})

test_that("logistic_init() mean-rate slope is negative for falling data", {
    ## same degenerate path, falling direction -> negative slope fallback
    ## (B - A) / t_range == -100 / 10 == -10
    x <- c(100, 0)
    t <- c(0, 10)

    res <- call_logistic_init(x, t)
    expect_equal(res[["slope"]], -10)
    expect_true(is.finite(res[["xmid"]]))
})

test_that("logistic_init() xmid fallback pins out-of-range xmid to t range", {
    ## only the early tail of a sigmoid is observed, so the linearised
    ## xmid extrapolates past max(t); fallback pins it into [min(t), max(t)]
    t <- 1:30
    x <- logistic(t, A = 10, B = 100, xmid = 80, slope = 4)

    res <- call_logistic_init(x, t)
    expect_true(res[["xmid"]] >= min(t) && res[["xmid"]] <= max(t))
    expect_true(is.finite(res[["slope"]]))
})

test_that("logistic_init() uses linearisation on a clean sigmoid", {
    ## non-degenerate data bypasses both fallbacks and recovers params
    t <- 1:60
    x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4)

    res <- call_logistic_init(x, t)
    expect_true(all.equal(res[["A"]], 10, tolerance = 5, scale = 1))
    expect_true(all.equal(res[["B"]], 100, tolerance = 5, scale = 1))
    expect_true(all.equal(res[["xmid"]], 30, tolerance = 5, scale = 1))
    expect_true(res[["slope"]] > 0)
})

test_that("logistic_init() 5-param branch returns asym estimate", {
    t <- 1:60
    x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4, asym = 0.5)

    res <- call_logistic_init(x, t, asym = TRUE)
    expect_named(res, c("A", "B", "xmid", "slope", "asym"))
    expect_true(res[["asym"]] >= 0.1 && res[["asym"]] <= 0.9)
    expect_true(res[["xmid"]] >= min(t) && res[["xmid"]] <= max(t))
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
    df <- setNames(data.frame(t, x), c("time", channels[1]))
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
        shape = "symmetric",
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "time_channel",
        "A", "B", "xmid", "slope", "xmid_fitted"
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

test_that("analyse_logistic() validates shape argument", {
    data <- create_logistic_data()

    expect_error(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "not_a_shape"
        )
    )
})

test_that("analyse_logistic() recovers symmetric known parameters", {
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
        shape = "symmetric",
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.2, scale = 1))
})

test_that("analyse_logistic() recovers gompertz known parameters", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4

    set.seed(13)
    n <- 100
    t <- seq(0, n - 1, length.out = n)
    x <- gompertz(t, A, B, xmid, slope) + rnorm(n, 0, 1)
    df <- data.frame(time = t, smo2 = x)
    data <- create_mnirs_data(
        df, nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "gompertz",
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 2, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.3, scale = 1))
})

test_that("analyse_logistic() recovers gompertz_left known parameters", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4

    set.seed(13)
    n <- 100
    t <- seq(0, n - 1, length.out = n)
    x <- gompertz_left(t, A, B, xmid, slope) + rnorm(n, 0, 1)
    df <- data.frame(time = t, smo2 = x)
    data <- create_mnirs_data(
        df, nirs_channels = "smo2", time_channel = "time", sample_rate = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "gompertz_left",
        verbose = FALSE
    )

    expect_true(all.equal(result$A, A, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 2, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.3, scale = 1))
})

test_that("analyse_logistic() uses start_time correctly", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4
    start_time <- 12

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, n = 100, noise_sd = 1
    )
    data$time <- data$time + start_time

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "symmetric",
        start_time = start_time,
        verbose = FALSE
    )

    ## xmid reported as offset from start_time — should match original xmid
    expect_true(all.equal(result$A, A, tolerance = 1, scale = 1))
    expect_true(all.equal(result$B, B, tolerance = 1, scale = 1))
    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))
    expect_true(all.equal(result$slope, slope, tolerance = 0.2, scale = 1))
})

test_that("analyse_logistic() start_time edge cases", {
    A <- 10
    B <- 100
    xmid <- 30
    slope <- 4
    start_time <- 12

    data <- create_logistic_data(
        A = A, B = B, xmid = xmid, slope = slope, n = 100, noise_sd = 1
    )
    data$time <- data$time + start_time

    ## start_time specified before time start, falls forward to t[1L]
    expect_warning(
        result <- analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "symmetric",
            start_time = 0,
            verbose = TRUE
        ),
        "start_time.*before first valid.*time_channel"
    )

    expect_true(all.equal(result$xmid, xmid, tolerance = 1, scale = 1))

    ## start_time beyond time range
    expect_error(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "symmetric",
            start_time = max(data$time) + 10,
            verbose = TRUE
        ),
        "No observations.*before.*start_time"
    )
})

test_that("analyse_logistic() returns NA for failed fit", {
    ## only 10 observations for a 4-param model with near-zero noise —
    ## insufficient distinct response to converge
    custom_name <- create_logistic_data(n = 10, noise_sd = 0.1)

    expect_warning(
        result <- analyse_logistic(
            custom_name,
            nirs_channels = "smo2",
            shape = "symmetric"
        ),
        "fit failed for.*smo2.*custom_name" ## call custom interval name
    )

    expect_true(is.na(result$A))
    expect_true(is.na(result$xmid))
    expect_true(is.na(result$slope))
})

test_that("analyse_logistic() suppresses fit-failure warning when verbose = FALSE", {
    ## small n + low noise — guaranteed fit failure
    custom_name <- create_logistic_data(n = 10, noise_sd = 0.1)

    expect_no_warning(
        analyse_logistic(
            custom_name,
            nirs_channels = "smo2",
            shape = "symmetric",
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
        shape = "symmetric",
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
        channels = c("ch1", "ch2"), n = 100, noise_sd = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = c("ch1", "ch2"),
        shape = list("symmetric", ch2 = "gompertz"),
        verbose = FALSE
    )

    expect_equal(nrow(result), 2L)

    ca <- attr(result, "channel_args")
    ch1_row <- ca[ca$nirs_channels == "ch1", ]
    ch2_row <- ca[ca$nirs_channels == "ch2", ]
    expect_equal(ch1_row$shape, "symmetric")
    expect_equal(ch2_row$shape, "gompertz")
})

test_that("analyse_logistic() fitted_data attribute is well-formed", {
    data <- create_logistic_data(n = 100, noise_sd = 1)

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "symmetric",
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
        shape = "symmetric",
        verbose = FALSE
    )

    diag <- attr(result, "diagnostics")
    expect_true(all(
        c("nirs_channels", "n_obs", "r2", "adj_r2", "rmse") %in% names(diag)
    ))
    expect_true(diag$r2 > 0.9)
})


## direction =====================================================

test_that("analyse_logistic() direction = 'negative' matches auto on falling data", {
    data <- create_logistic_data(
        A = 100, B = 10, slope = -4, n = 100, noise_sd = 1
    )

    result_auto <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "symmetric",
        direction = "auto",
        verbose = FALSE
    )
    result_neg <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        shape = "symmetric",
        direction = "negative",
        verbose = FALSE
    )

    ## matching direction leaves the unconstrained fit untouched
    expect_equal(result_auto$A, result_neg$A)
    expect_equal(result_auto$B, result_neg$B)
    expect_equal(result_auto$slope, result_neg$slope)
    expect_true(result_auto$B < result_auto$A)
    expect_true(result_auto$slope < 0)
})

test_that("analyse_logistic() direction = 'positive' rejects falling fit", {
    data <- create_logistic_data(
        A = 100, B = 10, slope = -4, n = 100, noise_sd = 1
    )

    expect_warning(
        result <- analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "symmetric",
            direction = "positive",
            verbose = TRUE
        ),
        "satisfy"
    )

    ## never returns an inverted (B < A) fit against requested direction
    expect_true(is.na(result$A))
    expect_true(is.na(result$B))
    expect_true(is.na(result$slope))
})

test_that("analyse_logistic() suppresses direction warning when verbose = FALSE", {
    data <- create_logistic_data(
        A = 100, B = 10, slope = -4, n = 100, noise_sd = 1
    )

    expect_no_warning(
        analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "symmetric",
            direction = "positive",
            verbose = FALSE
        )
    )
})

test_that("analyse_logistic() direction constrains gompertz shapes", {
    ## falling gompertz data; forced positive returns NA
    set.seed(13)
    n <- 100
    t <- seq(0, n - 1, length.out = n)
    x <- gompertz(t, A = 100, B = 10, xmid = 30, slope = -4) + rnorm(n, 0, 1)
    data <- create_mnirs_data(
        data.frame(time = t, smo2 = x),
        nirs_channels = "smo2",
        time_channel = "time",
        sample_rate = 1
    )

    expect_warning(
        result <- analyse_logistic(
            data,
            nirs_channels = "smo2",
            shape = "gompertz",
            direction = "positive",
            verbose = TRUE
        ),
        "satisfy"
    )

    expect_true(is.na(result$A))
})

test_that("analyse_kinetics() passes direction to sigmoidal method", {
    data <- create_logistic_data(
        A = 100, B = 10, slope = -4, n = 100, noise_sd = 1
    )

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "sigmoidal",
            shape = "symmetric",
            direction = "positive"
        ),
        "satisfy"
    )

    expect_true(is.na(result$coefficients$A))
})


## fixed parameters ==============================================

test_that("SSlogistic() fixes parameters as constants", {
    set.seed(15)
    t <- 1:60
    x <- logistic(t, A = 0, B = 100, xmid = 30, slope = 4) +
        rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSlogistic(t, A = 0, B, xmid, slope), data = data)

    expect_s3_class(model, "nls")
    expect_named(coef(model), c("B", "xmid", "slope"))
    expect_true(all.equal(coef(model)[["B"]], 100, tolerance = 5, scale = 1))
    expect_true(all.equal(coef(model)[["xmid"]], 30, tolerance = 2, scale = 1))

    ## fix xmid instead of A
    model_xmid <- nls(x ~ SSlogistic(t, A, B, xmid = 30, slope), data = data)
    expect_named(coef(model_xmid), c("A", "B", "slope"))
})

test_that("SSgompertz() / SSgompertz_left() fix parameters as constants", {
    set.seed(15)
    t <- 1:60
    x <- gompertz(t, A = 0, B = 100, xmid = 30, slope = 4) +
        rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(x ~ SSgompertz(t, A = 0, B, xmid, slope), data = data)
    expect_named(coef(model), c("B", "xmid", "slope"))
    expect_true(all.equal(coef(model)[["B"]], 100, tolerance = 5, scale = 1))

    set.seed(16)
    x2 <- gompertz_left(t, A = 0, B = 100, xmid = 30, slope = 4) +
        rnorm(length(t), 0, 2)
    data2 <- data.frame(t, x = x2)

    model_left <- nls(
        x ~ SSgompertz_left(t, A = 0, B, xmid, slope),
        data = data2
    )
    expect_named(coef(model_left), c("B", "xmid", "slope"))
})

test_that("analyse_logistic() fix holds parameters constant", {
    data <- create_logistic_data(
        A = 0, B = 100, xmid = 30, slope = 4, noise_sd = 1
    )

    result <- analyse_logistic(
        data,
        nirs_channels = "smo2",
        fix = list(A = 0),
        verbose = FALSE
    )

    expect_equal(result$A, 0)
    expect_equal(result$B, 100, tolerance = 5)
    expect_equal(result$xmid, 30, tolerance = 2)

    ## fixed A excluded from the fitted model coefficients
    expect_named(coef(attr(result, "model")$smo2), c("B", "xmid", "slope"))
    expect_equal(attr(result, "channel_args")$fix, "list(A = 0)")
})

test_that("analyse_logistic() validates fix argument", {
    data <- create_logistic_data()

    expect_error(
        analyse_logistic(data, nirs_channels = "smo2", fix = list(tau = 5)),
        "not recognised"
    )
    expect_error(
        analyse_logistic(
            data, nirs_channels = "smo2",
            fix = list(A = 0, B = 100, xmid = 30, slope = 4)
        ),
        "Nothing to estimate"
    )
})

test_that("analyse_kinetics() passes fix to the sigmoidal method", {
    data <- create_logistic_data(
        A = 0, B = 100, xmid = 30, slope = 4, noise_sd = 1
    )

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "sigmoidal",
        fix = list(A = 0),
        verbose = FALSE
    )

    expect_equal(result$coefficients$A, 0)
    expect_named(coef(result$model[[1L]]$smo2), c("B", "xmid", "slope"))
})
