## exponential_drift() ==============================================
test_that("exponential_drift() is monoexponential before texc and linear after", {
    t <- 0:120
    texc <- 5 * 8
    mono <- monoexponential(t, A = 70, B = 40, tau = 8)
    result <- exponential_drift(t, 70, 40, 8, slope = 0.05, tau_mult = 5)

    ## the hinge is exactly zero before texc = tau_mult * tau
    expect_equal(result[t <= texc], mono[t <= texc])
    expect_equal(result[t > texc], mono[t > texc] + 0.05 * (t[t > texc] - texc))

    ## TD form: flat at A before TD, hinge shifts by TD
    result_TD <- exponential_drift(t, 70, 40, 8, 0.05, 5, TD = 15)
    expect_true(all(result_TD[t < 15] == 70))
    expect_equal(
        result_TD,
        monoexponential(t, 70, 40, 8, 15) + 0.05 * pmax(t - 15 - texc, 0)
    )

    ## no drift reduces to the monoexponential
    expect_equal(exponential_drift(t, 70, 40, 8, 0, 5), mono)
})


## SSexponential_drift() ============================================
test_that("SSexponential_drift() fits the 6-parameter TD form", {
    set.seed(13)
    t <- 1:180
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, tau_mult = 4, TD = 15
    ) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## the hinge is non-smooth, so port may stop short of its convergence
    ## certificate on usable coefficients
    model <- nls(
        x ~ SSexponential_drift(t, A, B, tau, slope, tau_mult = 4, TD),
        data = data,
        algorithm = "port",
        lower = c(-Inf, -Inf, 0, -Inf, 0),
        control = nls.control(warnOnly = TRUE)
    )

    expect_s3_class(model, "nls")
    coefs <- coef(model)
    expect_named(coefs, c("A", "B", "tau", "slope", "TD"))
    expect_true(all.equal(coefs[["A"]], 10, tolerance = 3, scale = 1))
    expect_true(all.equal(coefs[["B"]], 100, tolerance = 3, scale = 1))
    expect_true(all.equal(coefs[["tau"]], 12, tolerance = 3, scale = 1))
    expect_true(all.equal(coefs[["slope"]], -0.5, tolerance = 0.1, scale = 1))
    expect_true(all.equal(coefs[["TD"]], 15, tolerance = 3, scale = 1))
})

test_that("SSexponential_drift() gradient matches numericDeriv for the free parameters", {
    ## sample points off the hinge, where the one-sided derivative is exact
    t <- seq(-10, 120, by = 0.5) + 0.1
    env <- list2env(list(
        t = t, A = 70, B = 40, tau = 8, slope = -0.2, tau_mult = 3, TD = 3
    ))
    chk <- function(expr, pars) {
        an <- attr(eval(expr, env), "gradient")
        nd <- attr(numericDeriv(expr, pars, env), "gradient")
        expect_identical(colnames(an), pars)
        expect_equal(unname(an), unname(nd), tolerance = 1e-5)
    }
    chk(
        quote(SSexponential_drift(t, A, B, tau, slope, tau_mult = 3, TD)),
        c("A", "B", "tau", "slope", "TD")
    )
    chk(
        quote(SSexponential_drift(t, A, B, tau, slope, tau_mult)),
        c("A", "B", "tau", "slope", "tau_mult")
    )
    expect_null(
        attr(exponential_drift(t, 70, 40, 8, -0.2, 3), "gradient")
    )
})

test_that("expdrift_start() matches a per-point least-squares grid search", {
    set.seed(8)
    t <- 0:150
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, tau_mult = 4, TD = 15
    ) + rnorm(length(t), 0, 2)
    start <- expdrift_start(x, t, fixed = list(tau_mult = 4), has_TD = TRUE)
    expect_named(start, c("A", "B", "tau", "slope", "tau_mult", "TD"))

    ## the linear coefficients at the chosen grid point are the lm solution
    e <- exp(-pmax(t - start[["TD"]], 0) / start[["tau"]])
    h <- pmax(t - start[["TD"]] - 4 * start[["tau"]], 0)
    cf <- lm.fit(cbind(e, 1 - e, h), x)$coefficients
    expect_equal(
        unname(start[c("A", "B", "slope")]), unname(cf), tolerance = 1e-8
    )
})

test_that("SSexponential_drift() fits the 5-parameter form with a fixed A", {
    set.seed(3)
    t <- 0:150
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, tau_mult = 4
    ) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSexponential_drift(t, A = 10, B, tau, slope, tau_mult = 4),
        data = data,
        algorithm = "port",
        lower = c(-Inf, 0, -Inf),
        control = nls.control(warnOnly = TRUE)
    )

    expect_named(coef(model), c("B", "tau", "slope"))
    expect_equal(unname(predict(model, data.frame(t = 0))[1]), 10)
    expect_true(all.equal(coef(model)[["B"]], 100, tolerance = 3, scale = 1))
})


## analyse_exponential_drift() ======================================

## helper: falling primary response with a late positive drift starting
## at texc = TD + tau_mult * tau = 37
create_expdrift_data <- function(
    A = 70,
    B = 40,
    tau = 8,
    slope = 0.2,
    tau_mult = 4,
    TD = 5,
    n = 120,
    sample_rate = 1,
    noise_sd = 0.3,
    channels = "smo2",
    seed = 42
) {
    set.seed(seed)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    df <- data.frame(time = t)
    ## successive channels are offset by 5 units
    df[channels] <- lapply(seq_along(channels) - 1L, \(.i) {
        exponential_drift(t, A + 5 * .i, B + 5 * .i, tau, slope, tau_mult, TD) +
            rnorm(n, 0, noise_sd)
    })

    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate
    )
}


test_that("analyse_exponential_drift() returns correct structure and recovers parameters", {
    result <- analyse_exponential_drift(
        create_expdrift_data(),
        nirs_channels = "smo2",
        tau_mult = 4,
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "A", "B", "tau", "k", "TD", "MRT", "HRT",
        "texc", "slope", "tau_mult", "MRT_fitted", "HRT_fitted", "texc_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    expect_s3_class(attr(result, "model")$smo2, "nls")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_equal(attr(result, "channel_args")$tau_mult, 4)

    ## the onset multiple is held, never estimated
    expect_named(
        coef(attr(result, "model")$smo2),
        c("A", "B", "tau", "slope", "TD")
    )
    expect_true(all.equal(result$A, 70, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B, 40, tolerance = 3, scale = 1))
    expect_true(all.equal(result$tau, 8, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, 0.2, tolerance = 0.05, scale = 1))
    expect_true(all.equal(result$TD, 5, tolerance = 3, scale = 1))
    expect_equal(result$tau_mult, 4)
    expect_true(attr(result, "diagnostics")$r2 > 0.9)

    ## derived columns follow the fitted coefficients
    expect_equal(result$k, 1 / result$tau)
    expect_equal(result$MRT, result$TD + result$tau)
    expect_equal(result$HRT, result$TD + result$tau * log(2))
    expect_equal(result$texc, result$TD + 4 * result$tau)
    fitted_at <- \(.t) {
        exponential_drift(
            .t, result$A, result$B, result$tau, result$slope, 4, result$TD
        )
    }
    expect_equal(result$MRT_fitted, fitted_at(result$MRT))
    expect_equal(result$texc_fitted, fitted_at(result$texc))
})

test_that("analyse_exponential_drift() use_TD = FALSE fits the 5-param model from start_time", {
    ## the reduced model has no flat region, so pre-onset rows are dropped
    start_time <- 20
    data <- create_expdrift_data(TD = 0)
    data$time <- data$time + start_time

    result <- analyse_exponential_drift(
        data,
        nirs_channels = "smo2",
        start_time = start_time,
        use_TD = FALSE,
        tau_mult = 4,
        verbose = FALSE
    )

    expect_named(
        coef(attr(result, "model")$smo2),
        c("A", "B", "tau", "slope")
    )
    expect_true(is.na(result$TD))
    expect_equal(result$MRT, result$tau)
    expect_equal(result$HRT, result$tau * log(2))
    expect_equal(result$texc, 4 * result$tau)
    expect_equal(
        attr(result, "diagnostics")$n_obs, sum(data$time >= start_time)
    )
})

test_that("analyse_exponential_drift() falls back and then fails on too few observations", {
    ## five observations under-determine the 5-free-parameter TD model but
    ## not the reduced model, so the TD fit is rejected and the retry
    ## announced
    warns <- capture_warnings(
        result <- analyse_exponential_drift(
            create_expdrift_data(n = 5, noise_sd = 0.1),
            nirs_channels = "smo2"
        )
    )
    expect_match(warns[[1L]], "6-parameter")
    expect_match(warns[[1L]], "Attempting")
    expect_true(is.na(result$TD))

    ## too few observations for either model
    custom_name <- create_expdrift_data(n = 3, noise_sd = 0.1)
    expect_warning(
        result <- analyse_exponential_drift(custom_name, "smo2"),
        "fit failed for.*smo2.*custom_name.*3 observations for 5 free"
    ) |>
        expect_warning(
            "fit failed for.*smo2.*custom_name.*3 observations for 4 free"
        )
    expect_true(all(is.na(result[c("A", "tau", "slope", "texc_fitted")])))
    expect_null(attr(result, "model")$smo2)
})

test_that("analyse_exponential_drift() tau_mult resolves per channel", {
    data <- create_expdrift_data(channels = c("smo2", "hhb"))

    result <- analyse_exponential_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        tau_mult = list(smo2 = 2, hhb = 4),
        verbose = FALSE
    )
    expect_equal(result$tau_mult, c(2, 4))
    expect_equal(result$texc, result$TD + c(2, 4) * result$tau)
    expect_equal(attr(result, "channel_args")$tau_mult, c(2, 4))
    expect_false("tau_mult" %in% names(coef(attr(result, "model")$smo2)))

    ## an omitted channel takes the formal default
    result_part <- analyse_exponential_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        tau_mult = list(smo2 = 2),
        verbose = FALSE
    )
    expect_equal(result_part$tau_mult, c(2, 3))
})

test_that("analyse_exponential_drift() validates tau_mult", {
    data <- create_expdrift_data()

    ## a zero multiple starts the drift at the response onset
    for (bad in list(0, -1, "3", c(2, 3))) {
        expect_error(
            analyse_exponential_drift(data, nirs_channels = "smo2", tau_mult = bad),
            "tau_mult"
        )
    }
})

test_that("analyse_exponential_drift() fix holds parameters constant", {
    ## no drift: the curve at texc is the primary response alone
    result <- analyse_exponential_drift(
        create_expdrift_data(slope = 0),
        nirs_channels = "smo2",
        fix = list(slope = 0),
        verbose = FALSE
    )
    expect_equal(result$slope, 0)
    expect_named(coef(attr(result, "model")$smo2), c("A", "B", "tau", "TD"))
    expect_equal(
        result$texc_fitted,
        monoexponential(result$texc, result$A, result$B, result$tau, result$TD)
    )

    ## a fixed TD is excluded from estimation and disables the 5-param retry
    data <- create_expdrift_data()
    result <- analyse_exponential_drift(
        data,
        nirs_channels = "smo2",
        fix = list(TD = 5),
        verbose = FALSE
    )
    expect_equal(result$TD, 5)
    expect_named(coef(attr(result, "model")$smo2), c("A", "B", "tau", "slope"))
    expect_equal(result$MRT, 5 + result$tau)

    ## TD is only fixable when use_TD = TRUE; tau_mult is never fixable
    expect_error(
        analyse_exponential_drift(
            data, nirs_channels = "smo2", use_TD = FALSE, fix = list(TD = 0)
        ),
        "not recognised"
    )
    expect_error(
        analyse_exponential_drift(
            data, nirs_channels = "smo2", fix = list(tau_mult = 4)
        ),
        "not recognised"
    )
})

test_that("analyse_exponential_drift() enforces direction", {
    data <- create_expdrift_data()

    ## matching direction leaves the unconstrained fit untouched
    result_auto <- analyse_exponential_drift(
        data, nirs_channels = "smo2", direction = "auto", verbose = FALSE
    )
    result_neg <- analyse_exponential_drift(
        data, nirs_channels = "smo2", direction = "negative", verbose = FALSE
    )
    cols <- c("A", "B", "tau", "slope")
    expect_equal(result_auto[cols], result_neg[cols])
    expect_true(result_auto$B < result_auto$A)
    expect_equal(attr(result_neg, "channel_args")$direction, "negative")

    ## no positive primary response exists, so the bounded refit degenerates
    expect_warning(
        result_pos <- analyse_exponential_drift(
            data, nirs_channels = "smo2", direction = "positive"
        ),
        "satisfy"
    )
    expect_true(all(is.na(result_pos[c("A", "B", "slope", "texc")])))
})


## integration tests ================================================

test_that("analyse_exponential_drift() converges on real dataset", {
    skip("Manual fit convergence check")

    intervals <- readRDS(test_path("testdata/5-1_intervals_short.rds"))
    deoxy <- intervals[grepl("^deoxy", names(intervals))]

    analyse_kinetics(
        deoxy,
        method = "exp-drift",
    ) |>
        plot()

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
