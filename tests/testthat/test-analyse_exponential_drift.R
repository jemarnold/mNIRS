## exponential_drift() ==============================================
test_that("exponential_drift() is monoexponential before the onset and linear after", {
    t <- 0:120
    ## the onset is where the primary reaches the fraction of its amplitude
    onset <- expdrift_onset(8, 0.99)
    expect_equal(onset, -8 * log(0.01))
    expect_equal(monoexponential(onset, 70, 40, 8), 70 + 0.99 * (40 - 70))
    expect_equal(expdrift_onset(8, 0.95, TD = 15), 15 - 8 * log(0.05))
    mono <- monoexponential(t, A = 70, B = 40, tau = 8)
    result <- exponential_drift(t, 70, 40, 8, slope = 0.05, drift_frac = 0.99)

    ## the hinge is exactly zero before the onset
    expect_equal(result[t <= onset], mono[t <= onset])
    expect_equal(
        result[t > onset], mono[t > onset] + 0.05 * (t[t > onset] - onset)
    )

    ## TD form: flat at A before TD, hinge shifts by TD
    result_TD <- exponential_drift(t, 70, 40, 8, 0.05, 0.99, TD = 15)
    expect_true(all(result_TD[t < 15] == 70))
    expect_equal(
        result_TD,
        monoexponential(t, 70, 40, 8, 15) + 0.05 * pmax(t - 15 - onset, 0)
    )

    ## no drift reduces to the monoexponential
    expect_equal(exponential_drift(t, 70, 40, 8, 0, 0.99), mono)
})


## SSexponential_drift() ============================================
test_that("SSexponential_drift() fits the 6-parameter TD form", {
    set.seed(13)
    t <- 1:180
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, drift_frac = 0.98, TD = 15
    ) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    ## the hinge is non-smooth, so port may stop short of its convergence
    ## certificate on usable coefficients
    model <- nls(
        x ~ SSexponential_drift(t, A, B, tau, slope, drift_frac = 0.98, TD),
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
        t = t, A = 70, B = 40, tau = 8, slope = -0.2, drift_frac = 0.95, TD = 3
    ))
    chk <- function(expr, pars) {
        an <- attr(eval(expr, env), "gradient")
        nd <- attr(numericDeriv(expr, pars, env), "gradient")
        expect_identical(colnames(an), pars)
        expect_equal(unname(an), unname(nd), tolerance = 1e-5)
    }
    chk(
        quote(SSexponential_drift(t, A, B, tau, slope, drift_frac = 0.95, TD)),
        c("A", "B", "tau", "slope", "TD")
    )
    chk(
        quote(SSexponential_drift(t, A, B, tau, slope, drift_frac)),
        c("A", "B", "tau", "slope", "drift_frac")
    )
    expect_null(
        attr(exponential_drift(t, 70, 40, 8, -0.2, 0.95), "gradient")
    )
})

test_that("expdrift_start() matches a per-point least-squares grid search", {
    set.seed(8)
    t <- 0:150
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, drift_frac = 0.98, TD = 15
    ) + rnorm(length(t), 0, 2)
    start <- expdrift_start(x, t, fixed = list(drift_frac = 0.98), has_TD = TRUE)
    expect_named(start, c("A", "B", "tau", "slope", "drift_frac", "TD"))

    ## the linear coefficients at the chosen grid point are the lm solution
    e <- exp(-pmax(t - start[["TD"]], 0) / start[["tau"]])
    h <- pmax(t - expdrift_onset(start[["tau"]], 0.98, start[["TD"]]), 0)
    cf <- lm.fit(cbind(e, 1 - e, h), x)$coefficients
    expect_equal(
        unname(start[c("A", "B", "slope")]), unname(cf), tolerance = 1e-8
    )
})

test_that("SSexponential_drift() fits the 5-parameter form with a fixed A", {
    set.seed(3)
    t <- 0:150
    x <- exponential_drift(
        t, A = 10, B = 100, tau = 12, slope = -0.5, drift_frac = 0.98
    ) + rnorm(length(t), 0, 2)
    data <- data.frame(t, x)

    model <- nls(
        x ~ SSexponential_drift(t, A = 10, B, tau, slope, drift_frac = 0.98),
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
## at the onset TD - tau * log(1 - drift_frac) = 36.3
create_expdrift_data <- function(
    A = 70,
    B = 40,
    tau = 8,
    slope = 0.2,
    drift_frac = 0.98,
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
        exponential_drift(
            t, A + 5 * .i, B + 5 * .i, tau, slope, drift_frac, TD
        ) +
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
        drift_frac = 0.98,
        verbose = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_named(result, c(
        "interval", "nirs_channels", "A", "B", "TD", "tau", "k", "MRT", "HRT",
        "texc", "slope", "drift_frac", "MRT_fitted", "HRT_fitted", "texc_fitted"
    ))
    expect_equal(nrow(result), 1L)

    ## attributes
    expect_s3_class(attr(result, "model")$smo2, "nls")
    expect_named(attr(result, "fitted_data")$smo2, c("window_idx", "fitted"))
    expect_equal(nrow(attr(result, "diagnostics")), 1L)
    expect_equal(attr(result, "channel_args")$drift_frac, 0.98)

    ## the onset fraction is held, never estimated
    expect_named(
        coef(attr(result, "model")$smo2),
        c("A", "B", "tau", "slope", "TD")
    )
    expect_true(all.equal(result$A, 70, tolerance = 2, scale = 1))
    expect_true(all.equal(result$B, 40, tolerance = 3, scale = 1))
    expect_true(all.equal(result$tau, 8, tolerance = 2, scale = 1))
    expect_true(all.equal(result$slope, 0.2, tolerance = 0.05, scale = 1))
    expect_true(all.equal(result$TD, 5, tolerance = 3, scale = 1))
    expect_equal(result$drift_frac, 0.98)
    expect_true(attr(result, "diagnostics")$r2 > 0.9)

    ## derived columns follow the fitted coefficients
    expect_equal(result$k, 1 / result$tau)
    expect_equal(result$MRT, result$TD + result$tau)
    expect_equal(result$HRT, result$TD + result$tau * log(2))
    ## the drift takes over before the onset, so texc is the onset itself
    expect_equal(result$texc, expdrift_onset(result$tau, 0.98, result$TD))
    fitted_at <- \(.t) {
        exponential_drift(
            .t, result$A, result$B, result$tau, result$slope, 0.98, result$TD
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
        drift_frac = 0.98,
        verbose = FALSE
    )

    expect_named(
        coef(attr(result, "model")$smo2),
        c("A", "B", "tau", "slope")
    )
    expect_true(is.na(result$TD))
    expect_equal(result$MRT, result$tau)
    expect_equal(result$HRT, result$tau * log(2))
    expect_equal(result$texc, expdrift_onset(result$tau, 0.98))
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

test_that("analyse_exponential_drift() drift_frac resolves per channel", {
    data <- create_expdrift_data(channels = c("smo2", "hhb"))

    result <- analyse_exponential_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_frac = list(smo2 = 0.85, hhb = 0.98),
        verbose = FALSE
    )
    expect_equal(result$drift_frac, c(0.85, 0.98))
    ## texc is the turning point when past the onset (smo2), else the onset
    onset <- expdrift_onset(result$tau, c(0.85, 0.98), result$TD)
    takeover <- result$TD +
        result$tau * log((result$A - result$B) / (result$slope * result$tau))
    expect_equal(result$texc, pmax(onset, takeover))
    expect_gt(result$texc[[1L]], onset[[1L]])
    expect_true(all.equal(
        result$texc_fitted[[1L]],
        min(attr(result, "fitted_data")$smo2$fitted),
        tolerance = 1, scale = 1
    ))
    expect_equal(attr(result, "channel_args")$drift_frac, c(0.85, 0.98))
    expect_false("drift_frac" %in% names(coef(attr(result, "model")$smo2)))

    ## an omitted channel takes the formal default
    result_part <- analyse_exponential_drift(
        data,
        nirs_channels = c("smo2", "hhb"),
        drift_frac = list(smo2 = 0.85),
        verbose = FALSE
    )
    expect_equal(result_part$drift_frac, c(0.85, 0.95))
})

test_that("analyse_exponential_drift() texc is the takeover point of a monotonic drift", {
    ## drift continues in the direction of the primary response: no turning
    ## point, so texc is where the drift rate exceeds the primary rate
    result <- analyse_exponential_drift(
        create_expdrift_data(slope = -0.2),
        nirs_channels = "smo2",
        drift_frac = 0.85,
        verbose = FALSE
    )
    expect_true(result$slope < 0)
    expect_gt(result$texc, expdrift_onset(result$tau, 0.85, result$TD))
    expect_equal(
        result$texc,
        result$TD +
            result$tau * log((result$A - result$B) / (-result$slope * result$tau))
    )
})

test_that("analyse_exponential_drift() validates drift_frac", {
    data <- create_expdrift_data()

    ## the drift must start past the half-response and before the asymptote
    for (bad in list(0.5, 1, 0, -0.1, "0.95", c(0.9, 0.95))) {
        expect_error(
            analyse_exponential_drift(
                data, nirs_channels = "smo2", drift_frac = bad
            ),
            "drift_frac.*must be a valid one-element"
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
    ## no takeover without drift: texc is the onset
    expect_equal(result$texc, expdrift_onset(result$tau, 0.95, result$TD))
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

    ## TD is only fixable when use_TD = TRUE; drift_frac is never fixable
    expect_error(
        analyse_exponential_drift(
            data, nirs_channels = "smo2", use_TD = FALSE, fix = list(TD = 0)
        ),
        "not recognised"
    )
    expect_error(
        analyse_exponential_drift(
            data, nirs_channels = "smo2", fix = list(drift_frac = 0.98)
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


## model fallback ===================================================

test_that("analyse_kinetics() keeps a supported drift", {
    data <- create_expdrift_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "exponential_drift",
        drift_frac = 0.98,
        verbose = FALSE
    )
    forced <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "exponential_drift",
        drift_frac = 0.98,
        model_fallback = FALSE,
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(names(cf)[1:4], c("interval", "nirs_channels", "start_time", "model"))
    expect_equal(cf$model, "exponential_drift")
    expect_true(all.equal(cf$slope, 0.2, tolerance = 0.05, scale = 1))
    expect_false(any(grepl("fell back to", result$warnings$message)))
    expect_equal(cf, forced$coefficients)
})

test_that("analyse_kinetics() falls back from a negligible drift", {
    data <- create_expdrift_data(slope = 0)

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "exponential_drift"
        ),
        "fell back to"
    )
    cf <- result$coefficients
    model <- result$model[[1L]]$smo2

    expect_equal(cf$model, "monoexponential")
    expect_named(coef(model), c("A", "B", "tau", "TD"))
    expect_equal(cf$tau, coef(model)[["tau"]])
    expect_true(all(is.na(cf[c("slope", "drift_frac", "texc", "texc_fitted")])))
    expect_equal(result$diagnostics$n_params, 4L)
    expect_equal(
        result$data[[1L]]$smo2_fitted,
        as.vector(predict(model))
    )
    msgs <- result$warnings$message
    expect_true(any(grepl("Drift amplitude", msgs)))

    ## the raw fit is kept on request
    forced <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "exponential_drift",
        model_fallback = FALSE,
        verbose = FALSE
    )
    expect_equal(forced$coefficients$model, "exponential_drift")
    expect_false(is.na(forced$coefficients$slope))
})

test_that("exponential_drift fallback resolves per channel with fix carried", {
    data <- create_expdrift_data(channels = c("smo2", "hhb"))
    data$hhb <- create_expdrift_data(slope = 0, seed = 1)$smo2

    result <- analyse_kinetics(
        data,
        nirs_channels = c(smo2, hhb),
        method = "exponential_drift",
        drift_frac = list(smo2 = 0.98, hhb = 0.98),
        fix = list(A = 70),
        verbose = FALSE
    )
    cf <- result$coefficients

    expect_equal(cf$nirs_channels, c("smo2", "hhb"))
    expect_equal(cf$model, c("exponential_drift", "monoexponential"))
    expect_equal(cf$A, c(70, 70))
    expect_named(coef(result$model[[1L]]$hhb), c("B", "tau", "TD"))
    fell <- result$warnings[grepl("fell back to", result$warnings$message), ]
    expect_equal(fell$nirs_channels, "hhb")
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
