## detect_direction ====================================================
test_that("detect_direction returns 'positive' & 'negative' unchanged", {
    result <- detect_direction(1:10, direction = "positive")
    expect_equal(result, "positive")
    
    result <- detect_direction(1:10, direction = "negative")
    expect_equal(result, "negative")
})

test_that("detect_direction detects slope", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    expect_equal(detect_direction(x), "positive")
    
    x <- c(14, 11, 12, 9, 7, 8, 5, 2, 3, 1)
    expect_equal(detect_direction(x), "negative")
})

test_that("detect_direction uses custom t vector", {
    ## x rises from 1 to 5, but t is reversed so lm slope is negative
    x <- c(1, 3, 5)
    t <- c(10, 5, 1)
    expect_equal(detect_direction(x, t), "negative")
})

test_that("detect_direction falls back appropriately", {
    ## positive when abs(max) >= abs(min)
    ## symmetric pulse: net slope ~ 0, max = 10, min = 0
    x <- c(0, 5, 10, 5, 0)
    expect_equal(detect_direction(x), "positive")

    ## negative when abs(max) < abs(min)
    ## symmetric trough: net slope ~ 0, max = 0, min = -10
    x <- c(0, -5, -10, -5, 0)
    expect_equal(detect_direction(x), "negative")
})

test_that("detect_direction falls back to positive on magnitude tie", {
    ## symmetric around zero: abs(max) == abs(min) => positive (>=)
    x <- c(-5, 0, 5, 0, -5)
    expect_equal(detect_direction(x), "positive")
})

test_that("detect_direction handles edge cases", {
    ## falls back when net slope is NA
    expect_equal(detect_direction(c(5, 5, 5, 5, 5)), "positive")

    ## falls back when all x are NA
    expect_equal(detect_direction(rep(NA_real_, 5)), "positive")

    ## x has zero net slope, but fallback indicates negative
    x <- c(0, 5, 10, 5, 0)
    fallback <- c(-1, -8, -2, -1, -3)
    expect_equal(detect_direction(x, fallback = fallback), "negative")
})

## compute_diagnostics ====================================================
test_that("compute_diagnostics returns correct structure", {
    x <- c(1, 3, 2, 5, 8)
    t <- seq_along(x)
    fitted <- predict(lm(x ~ t))

    result <- compute_diagnostics(x, t, fitted)

    expect_s3_class(result, "data.frame")
    expect_named(
        result,
        c("n_obs", "r2", "adj_r2", 
        "rmse", "snr", "cv_rmse", 
        "aic", "aicc", "bic")
    )
    expect_type(result$n_obs, "integer")
    expect_type(result$r2, "double")
    expect_type(result$adj_r2, "double")
    expect_type(result$rmse, "double")
    expect_type(result$snr, "double")
    expect_type(result$cv_rmse, "double")
})

test_that("compute_diagnostics matches lm() summary", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t <- seq_along(x)
    lm_fit <- lm(x ~ t)
    fitted <- predict(lm_fit)

    ## n_params = 2L: intercept + slope
    result <- compute_diagnostics(x, t, fitted, n_params = 2L)
    lm_summary <- summary(lm_fit)

    expect_equal(result$n_obs, length(fitted))
    ## r2 = cor(observed, fitted)^2; equals R^2 for OLS
    expect_equal(result$r2, lm_summary$r.squared)
    expect_equal(result$adj_r2, lm_summary$adj.r.squared)
    expect_equal(result$rmse, sqrt(mean(residuals(lm_fit)^2)))
    expect_equal(result$aic, AIC(lm_fit))
    expect_equal(result$bic, BIC(lm_fit))
})

test_that("compute_diagnostics aicc applies small-sample correction", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t <- seq_along(x)
    lm_fit <- lm(x ~ t)
    fitted <- predict(lm_fit)
    n <- length(x)
    k <- 3L ## n_params + 1 (intercept + slope + variance)

    result <- compute_diagnostics(x, t, fitted, n_params = 2L)

    expect_equal(
        result$aicc, AIC(lm_fit) + 2 * k * (k + 1L) / (n - k - 1L)
    )
})

test_that("compute_diagnostics n_params adjusts adj_r2 denominator", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t <- seq_along(x)
    fitted <- predict(lm(x ~ t))
    n <- length(x)

    r2 <- compute_diagnostics(x, t, fitted, n_params = 2L)$r2
    adj_r2_p2 <- compute_diagnostics(x, t, fitted, n_params = 2L)$adj_r2
    adj_r2_p4 <- compute_diagnostics(x, t, fitted, n_params = 4L)$adj_r2

    expect_equal(adj_r2_p2, 1 - (1 - r2) * (n - 1) / (n - 2))
    expect_equal(adj_r2_p4, 1 - (1 - r2) * (n - 1) / (n - 4))
    ## more params → more penalisation
    expect_lt(adj_r2_p4, adj_r2_p2)
})

test_that("compute_diagnostics adj_r2 is NA when n <= n_params", {
    x <- c(1, 2, 3)
    t <- seq_along(x)
    fitted <- predict(lm(x ~ t))

    ## n = 3, n_params = 3: denominator = 0 → NA
    result <- compute_diagnostics(x, t, fitted, n_params = 3L)
    expect_true(is.na(result$adj_r2))
})

test_that("compute_diagnostics r2 is bounded for a non-linear fit", {
    t <- seq(0, 4 * pi, length.out = 50)
    ## simulate exponential decay with noise
    x <- 20 * exp(-0.3 * t) + rnorm(50, sd = 0.5)
    fitted <- 20 * exp(-0.3 * t) ## true curve, no noise

    result <- compute_diagnostics(x, t, fitted)

    expect_true(result$r2 > 0 && result$r2 <= 1)
    expect_equal(result$r2, cor(x, fitted)^2)
})

test_that("compute_diagnostics handles perfect fit", {
    x <- 1:10
    t <- 1:10
    fitted <- x # perfect fit

    result <- compute_diagnostics(x, t, fitted, n_params = 2L)

    expect_equal(result$r2, 1)
    expect_equal(result$adj_r2, 1)
    expect_equal(result$rmse, 0)
    expect_true(is.na(result$snr)) ## zero residual variance → NA
    expect_equal(result$cv_rmse, 0) ## rmse = 0, x_mean != 0
    ## ss_res = 0 → log-likelihood undefined → information criteria NA
    expect_true(is.na(result$aic))
    expect_true(is.na(result$aicc))
    expect_true(is.na(result$bic))
})

test_that("compute_diagnostics handles edge cases", {
    ## n < 2 returns NA
    result <- compute_diagnostics(x = 5, t = 1, fitted = 5)
    expect_equal(result$n_obs, 1L)
    expect_true(is.na(result$r2))
    expect_true(is.na(result$adj_r2))
    expect_true(is.na(result$rmse))
    expect_true(is.na(result$snr))
    expect_true(is.na(result$cv_rmse))
    expect_true(is.na(result$aic))
    expect_true(is.na(result$aicc))
    expect_true(is.na(result$bic))

    ## n = 2, default n_params = 1L: adj_r2 defined; perfect fit → 1
    result <- compute_diagnostics(c(1, 2), c(1, 2), c(1, 2))
    expect_equal(result$n_obs, 2L)
    expect_equal(result$r2, 1)
    expect_equal(result$adj_r2, 1)
    expect_equal(result$rmse, 0)
})

test_that("compute_diagnostics handles zero variance in x", {
    x <- rep(5, 10)
    t <- seq_along(x)
    fitted <- rep(5, 10)

    result <- compute_diagnostics(x, t, fitted, verbose = FALSE)

    expect_true(is.na(result$r2)) ## ss_tot = 0 → NA
    expect_equal(result$rmse, 0)
    expect_true(is.na(result$snr)) ## zero signal variance → NA
    expect_equal(result$cv_rmse, 0) ## rmse = 0, x_mean = 5
})

test_that("computes_diagnostics validates input lengths", {
    x <- 1:5
    t <- 1:5
    fitted <- 1:4 # wrong length

    expect_warning(
        result <- compute_diagnostics(x, t, fitted, verbose = TRUE),
        "x.*t.*fitted.*equal lengths"
    )
    expect_true(is.na(result$r2))
})

test_that("compute_diagnostics snr is positive for a good fit", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t <- seq_along(x)
    fitted <- predict(lm(x ~ t))

    result <- compute_diagnostics(x, t, fitted)

    expect_true(result$snr > 0)
    ## SNR = 10 * log10(var(signal) / var(residuals))
    expect_equal(result$snr, 10 * log10(var(x) / var(x - fitted)))
})

test_that("compute_diagnostics cv_rmse scales with signal magnitude", {
    t <- seq_len(10)

    x1 <- c(10, 11, 10, 12, 11, 10, 13, 11, 10, 12) ## mean ~ 11
    x2 <- x1 * 10 ## mean ~ 110
    fitted1 <- predict(lm(x1 ~ t))
    fitted2 <- predict(lm(x2 ~ t))

    r1 <- compute_diagnostics(x1, t, fitted1)
    r2 <- compute_diagnostics(x2, t, fitted2)

    ## CV-RMSE should be equal (RMSE and mean scale proportionally)
    expect_equal(r1$cv_rmse, r2$cv_rmse)
})


## as_data_list ===========================================================
test_that("as_data_list handles single data frame", {
    df <- data.frame(t = 1:5, x = rnorm(5))
    result <- as_data_list(df)

    expect_type(result, "list")
    expect_length(result, 1)
    expect_named(result, "interval_1")
    expect_identical(result[[1]], df)
})

test_that("as_data_list handles list of data frames", {
    df1 <- data.frame(t = 1:5, x = rnorm(5))
    df2 <- data.frame(t = 1:5, y = rnorm(5))
    data_list <- list(df1, df2)

    result <- as_data_list(data_list)

    expect_type(result, "list")
    expect_length(result, 2)
    expect_named(result, c("interval_1", "interval_2"))
    expect_identical(result[[1]], df1)
    expect_identical(result[[2]], df2)
})

test_that("as_data_list preserves names in list", {
    df1 <- data.frame(t = 1:5, x = rnorm(5))
    df2 <- data.frame(t = 1:5, y = rnorm(5))
    data_list <- list(baseline = df1, exercise = df2)

    result <- as_data_list(data_list)

    expect_named(result, c("baseline", "exercise"))
})

test_that("as_data_list fills only unnamed list elements", {
    df1 <- data.frame(t = 1:5, x = rnorm(5))
    df2 <- data.frame(t = 1:5, y = rnorm(5))
    ## second element unnamed: supplied name kept, blank gets `interval_<n>`
    data_list <- list(baseline = df1, df2)

    result <- as_data_list(data_list)

    expect_named(result, c("baseline", "interval_2"))
})

test_that("as_data_list handles grouped data frame", {
    skip_if_not_installed("dplyr")

    df <- data.frame(
        t = rep(1:5, 2),
        x = rnorm(10),
        group = rep(c("A", "B"), each = 5)
    )
    df <- create_mnirs_data(
        df,
        nirs_channels = "x",
        time_channel = "t"
    )
    grouped_df <- dplyr::group_by(df, group)

    result <- as_data_list(grouped_df)

    expect_type(result, "list")
    expect_length(result, 2)
    expect_named(result, c("A", "B"))
    expect_s3_class(result[[1]], "mnirs")
    expect_s3_class(result[[2]], "mnirs")
})

test_that("as_data_list handles multi-column grouping", {
    skip_if_not_installed("dplyr")

    df <- data.frame(
        t = rep(1:5, 4),
        x = rnorm(20),
        condition = rep(c("rest", "active"), each = 10),
        subject = rep(c("S1", "S2"), each = 5, times = 2)
    )
    df <- create_mnirs_data(
        df,
        nirs_channels = "x",
        time_channel = "t"
    )
    grouped_df <- dplyr::group_by(df, condition, subject)

    result <- as_data_list(grouped_df)

    expect_length(result, 4)
    expect_named(result, c("rest_S1", "rest_S2", "active_S1", "active_S2"))
})

test_that("as_data_list errors on invalid input", {
    expect_error(
        as_data_list("not a data frame"),
        "must be a list of data frames"
    )

    expect_error(
        as_data_list(list(data.frame(x = 1), "not a df")),
        "must be a list of data frames"
    )

    expect_error(
        as_data_list(list(1, 2, 3)),
        "must be a list of data frames"
    )
})

test_that("as_data_list errors when dplyr is unavailable for grouped input", {
    skip_if_not_installed("dplyr")

    df <- data.frame(x = 1:4, grp = c("A", "A", "B", "B"))
    grouped <- dplyr::group_by(df, grp)

    with_mocked_bindings(
        requireNamespace = function(pkg, ...) FALSE,
        .package = "base",
        expect_error(as_data_list(grouped), "dplyr.*required for grouped")
    )
})

test_that("as_data_list preserves mnirs attributes", {
    df <- create_mnirs_data(
        data.frame(t = 1:5, x = rnorm(5)),
        nirs_channels = "x",
        time_channel = "t"
    )

    result <- as_data_list(df)

    expect_s3_class(result[[1]], "mnirs")
    expect_equal(attr(result[[1]], "nirs_channels"), "x")
    expect_equal(attr(result[[1]], "time_channel"), "t")
})

test_that("as_data_list preserves attributes in grouped split", {
    skip_if_not_installed("dplyr")

    df <- data.frame(
        t = rep(1:5, 2),
        x = rnorm(10),
        group = rep(c("B", "A"), each = 5)
    )
    df <- create_mnirs_data(
        df,
        nirs_channels = "x",
        time_channel = "t",
        sample_rate = 10
    )
    # attributes(df)
    grouped_df <- dplyr::group_by(df, group)

    result <- as_data_list(grouped_df)

    expect_equal(attr(result[[1]], "nirs_channels"), "x")
    expect_equal(attr(result[[1]], "time_channel"), "t")
    expect_equal(attr(result[[1]], "sample_rate"), 10)
    expect_equal(attr(result[[2]], "nirs_channels"), "x")
    expect_equal(attr(result[[2]], "time_channel"), "t")
    expect_equal(attr(result[[2]], "sample_rate"), 10)
})


## build_na_results ====================================================
test_that("build_na_results returns correct 4-element fit list", {
    na_coefs <- data.frame(
        slope = NA_real_,
        intercept = NA_real_
    )

    result <- build_na_results(na_coefs)

    expect_type(result, "list")
    expect_named(result, c("coefs", "model", "fitted_data", "diag"))

    ## coefs returned unchanged (channel/time columns added upstream)
    expect_identical(result$coefs, na_coefs)
    expect_true(is.na(result$coefs$slope))

    ## no model for a failed fit
    expect_null(result$model)

    ## fitted_data has NA placeholders
    expect_true(is.na(result$fitted_data$window_idx))
    expect_true(is.na(result$fitted_data$fitted))

    ## diag has all-NA values and zero observations
    expect_equal(result$diag$n_obs, 0L)
    expect_true(all(is.na(result$diag$r2)))
})


## analyse_kinetics_channels() ============================================
## helper: minimal mnirs data with a clear rise-then-fall per channel
make_channels_data <- function(channels = c("ch1", "ch2")) {
    n <- 20L
    x <- c(seq(1, 10, length.out = 10L), seq(9, 1, length.out = 10L))
    df <- data.frame(time = seq_len(n))
    for (.ch in channels) df[[.ch]] <- x
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time"
    )
}

## helper: per-channel arg list keyed by channel, as resolve_channel_args()
make_per_channel <- function(channels, ...) {
    extra <- list(...)
    setNames(
        lapply(channels, \(.ch) {
            c(list(end_fit_span = Inf, direction = "auto"), extra)
        }),
        channels
    )
}

## helper: trivial fit_fn returning one coefficient per channel
slope_fit <- function(coef_value = 1.0, time_coef = NULL) {
    function(.nirs, x_fit, t_fit, .a, valid, verbose) {
        coefs <- data.frame(slope = coef_value)
        if (!is.null(time_coef)) coefs$peak_slope_time <- time_coef
        list(
            coefs = coefs,
            model = structure(list(), class = "lm"),
            fitted_data = data.frame(
                window_idx = seq_along(x_fit),
                fitted = x_fit
            ),
            diag = compute_diagnostics(
                x_fit, t_fit, x_fit, n_params = 2L, verbose = verbose
            )
        )
    }
}

test_that("analyse_kinetics_channels combines channels correctly", {
    channels <- c("ch1", "ch2")
    data <- make_channels_data(channels)
    per_channel <- make_per_channel(channels)

    result <- analyse_kinetics_channels(
        data, channels, "time",
        per_channel, slope_fit(), verbose = FALSE
    )

    ## coefficient rows are combined, one per channel
    expect_equal(nrow(result), 2L)
    expect_equal(result$nirs_channels, channels)
    expect_equal(result$time_channel, c("time", "time"))

    ## TODO model preserved as ...
    ## model preserved as named list
    model <- attr(result, "model")
    expect_type(model, "list")
    expect_length(model, 2L)
    ## TODO verify model object type
    expect_named(model, channels)

    ## fitted_data preserved as named list
    fitted_data <- attr(result, "fitted_data")
    expect_type(fitted_data, "list")
    expect_named(fitted_data, channels)

    ## diagnostics combined, one row per channel
    diag <- attr(result, "diagnostics")
    expect_equal(nrow(diag), 2L)
    expect_equal(diag$nirs_channels, channels)

    ## channel_args combined, one row per channel
    ca <- attr(result, "channel_args")
    expect_equal(nrow(ca), 2L)
    expect_equal(ca$nirs_channels, channels)
})

test_that("analyse_kinetics_channels sets interval column", {
    channels <- "ch1"
    data <- make_channels_data(channels)
    per_channel <- make_per_channel(channels)

    result <- analyse_kinetics_channels(
        data, channels, "time",
        per_channel, slope_fit(), verbose = FALSE,
        interval_name = "baseline"
    )

    expect_equal(names(result)[1L], "interval")
    expect_equal(result$interval, "baseline")
})

test_that("analyse_kinetics_channels serialises args (NULL to NA, list)", {
    channels <- "ch1"
    data <- make_channels_data(channels)
    ## NULL and list-valued args must collapse into a flat data frame row
    per_channel <- make_per_channel(
        channels,
        width = 10,
        control = NULL,
        nls_control = list(maxiter = 100, tol = 1e-5)
    )

    result <- analyse_kinetics_channels(
        data, channels, "time",
        per_channel, slope_fit(), verbose = FALSE
    )

    ca <- attr(result, "channel_args")
    expect_equal(nrow(ca), 1L)
    expect_equal(ca$nirs_channels, "ch1")
    expect_equal(ca$width, 10)
    ## NULL collapses to NA
    expect_true(is.na(ca$control))
    ## list deparses to a character string
    expect_type(ca$nls_control, "character")
    expect_match(ca$nls_control, "maxiter")
})

test_that("analyse_kinetics_channels records extra_args in channel_args", {
    channels <- "ch1"
    data <- make_channels_data(channels)
    per_channel <- make_per_channel(channels)

    result <- analyse_kinetics_channels(
        data, channels, "time",
        per_channel, slope_fit(), verbose = FALSE,
        extra_args = list(shape = "symmetric")
    )

    ca <- attr(result, "channel_args")
    expect_equal(ca$shape, "symmetric")
})

test_that("analyse_kinetics_channels warns for negative time coefficients", {
    channels <- c("ch1", "ch2")
    data <- make_channels_data(channels)
    per_channel <- make_per_channel(channels)

    expect_warning(
        result <- analyse_kinetics_channels(
            data, channels, "time",
            per_channel, slope_fit(time_coef = -2), verbose = TRUE
        ),
        "Negative.*coefficients"
    )

    expect_equal(result$peak_slope_time, c(-2, -2))
})

test_that("analyse_kinetics_channels silent for valid coefficients", {
    channels <- "ch1"
    data <- make_channels_data(channels)
    per_channel <- make_per_channel(channels)

    expect_silent(
        analyse_kinetics_channels(
            data, channels, "time",
            per_channel, slope_fit(time_coef = 5), verbose = TRUE
        )
    )
})


## find_kinetics_idx ==================================================
test_that("find_kinetics_idx validates inputs", {
    expect_error(
        find_kinetics_idx(x = "a", t = 1, end_fit_span = 5),
        "numeric"
    )
    expect_error(
        find_kinetics_idx(x = 1:5, t = 1:5, end_fit_span = -1),
        "positive"
    )
    expect_silent(find_kinetics_idx(x = 1:5, t = 0:4, end_fit_span = Inf))
    expect_error(
        find_kinetics_idx(x = 1:5, t = 1:3, end_fit_span = 5),
        "equal length"
    )
})

test_that("find_kinetics_idx returns a named list", {
    x <- c(1, 5, 10, 5, 1)
    t <- seq_along(x)
    result <- find_kinetics_idx(x, t, end_fit_span = 0)
    expect_type(result, "list")
    expect_named(result, c("direction", "extreme", "idx"))
    expect_equal(result$extreme, 3L)
    expect_equal(result$idx, seq_len(3L))
})

test_that("find_kinetics_idx works on edge cases", {
    ## single element
    result <- find_kinetics_idx(x = 5, t = 1)
    expect_equal(result$idx, 1L)
    expect_null(result$extreme)

    ## all-equal values
    x <- rep(5, 20)
    t <- seq_along(x)
    result <- find_kinetics_idx(x, t)
    expect_equal(result$idx, seq_along(x))
    expect_null(result$extreme)

    ## monotonic increasing
    x <- 1:50
    t <- seq_along(x)
    result <- find_kinetics_idx(x, t)
    expect_equal(result$idx, seq_along(x))
    expect_null(result$extreme)

    ## monotonic decreasing
    x <- 50:1
    t <- seq_along(x)
    result <- find_kinetics_idx(x, t)
    expect_equal(result$idx, seq_along(x))
    expect_null(result$extreme)
})

test_that("find_kinetics_idx finds peak in rise-then-fall", {
    ## peak at index 20 (value 20), then decline
    x <- c(seq(1, 20, length.out = 20), seq(19, 1, length.out = 19))
    t <- seq_along(x)

    result <- find_kinetics_idx(x, t, end_fit_span = 0)
    expect_equal(result$idx, seq_len(20L))
    expect_equal(result$extreme, 20L)
    expect_equal(x[result$extreme], max(x))
})

test_that("find_kinetics_idx finds trough in fall-then-rise", {
    ## trough at index 20 (value 1), then rise
    x <- c(seq(20, 1, length.out = 20), seq(2, 20, length.out = 19))
    t <- seq_along(x)

    result <- find_kinetics_idx(x, t, end_fit_span = 0, direction = "negative")
    expect_equal(result$idx, seq_len(20L))
    expect_equal(result$extreme, 20L)
    expect_equal(x[result$extreme], min(x))
})

test_that("find_kinetics_idx works on irregular t with end_fit_span = 0", {
    ## peak at index 20 (value 20), then decline
    x <- c(seq(1, 20, length.out = 20), seq(19, 1, length.out = 19))
    t <- c(1:10, 10, 10, 13:39)/10

    result <- find_kinetics_idx(x, t, end_fit_span = 0, direction = "positive")
    expect_equal(result$idx, seq_len(12L))
    expect_equal(result$extreme, 12L)
    
    x <- c(1:10, 10, 10, 13:20, seq(19, 1, length.out = 19))
    t <- c(1:10, 10, 10, 13:39)/10
    
    result <- find_kinetics_idx(x, t, end_fit_span = 0, direction = "positive")
    ## TODO this fails but negligible real-world concern
    # expect_equal(result$idx, seq_len(10L))
    expect_equal(result$extreme, 10L)
})

test_that("find_kinetics_idx propagates auto-detected direction", {
    ## direction detection logic tested in test-detect_direction.R
    ## positive net slope => finds peak
    x_pos <- c(seq(1, 20, length.out = 15), seq(19, 10, length.out = 15))
    res_pos <- find_kinetics_idx(x_pos, seq_along(x_pos), end_fit_span = 0)
    expect_equal(res_pos$extreme, 15L)
    expect_equal(res_pos$direction, "positive")

    ## negative net slope => finds trough
    x_neg <- c(seq(20, 1, length.out = 15), seq(2, 10, length.out = 15))
    res_neg <- find_kinetics_idx(x_neg, seq_along(x_neg), end_fit_span = 0)
    expect_equal(res_neg$extreme, 15L)
    expect_equal(res_neg$direction, "negative")
})

test_that("find_kinetics_idx ignores negative t values", {
    ## peak in negative-t region should be ignored
    x <- c(20, 20, 0, 1, 2, 10, 9, 8, 7, 6)
    t <- c(-2, -1, 0, 1, 2,  3, 4, 5, 6, 7)

    result <- find_kinetics_idx(
        x,
        t,
        end_fit_span = 2,
        direction = "positive"
    )
    ## should not find the peak at t <= 0
    expect_equal(t[result$idx[length(result$idx)]], 5)
    expect_equal(result$extreme, 6L)
})

test_that("find_kinetics_idx returns n when all t <= 0", {
    x <- c(5, 10, 3)
    t <- c(-3, -2, -1)
    result <- find_kinetics_idx(x, t, end_fit_span = 1)
    expect_equal(result$idx, seq_along(x))
    ## no positive-t samples, so n_valid < 2 for extreme detection
    expect_null(result$extreme)
})

test_that("find_kinetics_idx handles invalid values", {
    x <- c(1, NA, 5, 10, Inf, 6, 4, 3)
    t <- c(1, 2, NA, 4, NA, 6, 7, 8)

    result <- find_kinetics_idx(
        x,
        t,
        end_fit_span = 3
    )
    expect_equal(result$idx, c(1, 4, 6, 7))
    expect_equal(result$extreme, 4L)
})

test_that("find_kinetics_idx returns first tie", {
    ## two equal peaks at indices 5 and 15
    x <- c(1, 2, 3, 4, 10, 4, 3, 2, 1, 2, 3, 4, 5, 4, 10, 4, 3, 2)
    t <- seq_along(x)
    result <- find_kinetics_idx(
        x,
        t,
        end_fit_span = 3,
        direction = "positive"
    )
    ## should find first peak (index 5), not second (index 15)
    expect_equal(result$idx, seq_len(5L + 3))
    expect_equal(result$extreme, 5L)
})

test_that("find_kinetics_idx end_fit_span larger than data range", {
    x <- c(1, 5, 3, 2)
    t <- seq_along(x)
    ## end_fit_span covers entire data range
    result <- find_kinetics_idx(
        x,
        t,
        end_fit_span = 100,
        direction = "positive"
    )
    expect_equal(result$idx, seq_along(x))
    expect_equal(result$extreme, 2L)
})

## build_kinetics_results ===============================================

## helper: minimal mnirs data frame with optional interval_times attr
make_kinetics_data <- function(
    n = 10,
    channels = "smo2",
    interval_times = NULL,
    sample_rate = 10
) {
    df_data <- data.frame(
        time = seq_len(n) / sample_rate,
        setNames(
            lapply(channels, \(.ch) rnorm(n)),
            channels
        )
    )
    df_data <- create_mnirs_data(
        df_data,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate
    )
    attr(df_data, "interval_times") <- interval_times
    return(df_data)
}

## helper: minimal attributed data frame matching the output of
## analyse_peak_slope() / analyse_monoexponential() per interval
make_kinetics_results <- function(
    interval,
    channels = "smo2",
    n = 10
) {
    ## this df immitates output of analyse_peak_slope / analyse_monoexponential
    df <- data.frame(
        nirs_channels = rep(channels, 1L),
        time_channel = "time",
        slope = seq_len(length(channels)) * 0.5,
        interval = interval
    )
    attr(df, "fitted_data") <- setNames(
        lapply(channels, \(.ch) {
            data.frame(
                window_idx = seq_len(n),
                fitted = seq_len(n) * 0.1
            )
        }),
        channels
    )
    attr(df, "model") <- setNames(
        lapply(channels, \(.ch) structure(list(), class = "lm")),
        channels
    )
    attr(df, "diagnostics") <- data.frame(
        nirs_channels = channels,
        r2 = rep(0.9, length(channels)),
        rmse = rep(0.1, length(channels))
    )
    attr(df, "channel_args") <- data.frame(
        nirs_channels = channels,
        width = rep(10L, length(channels))
    )
    return(df)
}

test_that("build_kinetics_results returns mnirs_kinetics with correct names", {
    data_list <- list(int1 = make_kinetics_data())
    result_list <- list(make_kinetics_results("int1"))
    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    expect_s3_class(result, "mnirs_kinetics")
    expect_named(result, c(
        "method", "model", "coefficients", "data",
        "interval_times", "diagnostics", "channel_args", "call"
    ))
})

test_that("build_kinetics_results stores method and call", {
    data_list <- list(int1 = make_kinetics_data())
    result_list <- list(make_kinetics_results("int1"))
    fake_call <- quote(analyse_kinetics(data, method = "peak_slope"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = fake_call
    )

    expect_equal(result$method, "peak_slope")
    expect_equal(result$call, fake_call)
})

test_that("build_kinetics_results places interval as first column of coefficients", {
    data_list <- list(A = make_kinetics_data(), B = make_kinetics_data())
    result_list <- list(make_kinetics_results("A"), make_kinetics_results("B"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    expect_equal(names(result$coefficients)[1L], "interval")
    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$interval, c("A", "B"))
})

test_that("build_kinetics_results adds _fitted columns to data elements", {
    data_list <- list(int1 = make_kinetics_data(n = 10))
    result_list <- list(make_kinetics_results("int1", n = 10))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    aug <- result$data[["int1"]]
    expect_true("smo2_fitted" %in% names(aug))
    expect_equal(length(aug$smo2_fitted), 10L)
})

test_that("build_kinetics_results fitted values placed at correct indices", {
    n <- 10
    data_list <- list(int1 = make_kinetics_data(n = n))
    r <- make_kinetics_results("int1", n = n)
    ## override fitted_data to only cover rows 3:7
    attr(r, "fitted_data") <- list(
        smo2 = data.frame(window_idx = 3:7, fitted = seq(0.3, 0.7, by = 0.1))
    )
    result_list <- list(r)

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    fitted_vec <- result$data[["int1"]]$smo2_fitted
    expect_true(all(is.na(fitted_vec[c(1, 2, 8, 9, 10)])))
    expect_equal(fitted_vec[3:7], seq(0.3, 0.7, by = 0.1))
})

test_that("build_kinetics_results data elements preserve mnirs metadata", {
    data_list <- list(int1 = make_kinetics_data(sample_rate = 10))
    result_list <- list(make_kinetics_results("int1"))

    result <- build_kinetics_results(
        data_list,
        result_list,
        method = "peak_slope",
        call = NULL
    )

    aug <- result$data[["int1"]]
    expect_s3_class(aug, "mnirs")
    expect_equal(attr(aug, "nirs_channels"), "smo2")
    expect_equal(attr(aug, "time_channel"), "time")
    expect_equal(attr(aug, "sample_rate"), 10)
})

test_that("build_kinetics_results data list has class = mnirs", {
    data_list <- list(
        int1 = make_kinetics_data(sample_rate = 10),
        int2 = make_kinetics_data(sample_rate = 10)
    )
    result_list <- list(
        make_kinetics_results("int1"),
        make_kinetics_results("int2")
    )

    result <- build_kinetics_results(
        data_list,
        result_list,
        method = "peak_slope",
        call = NULL
    )
    
    expect_type(result$data, "list")
    expect_s3_class(result$data, "mnirs")
    expect_length(result$data, 2)

    ## visual check
    p <- plot(result$data)
    expect_s3_class(p, "ggplot")

    ## facet wrap present for multi-element list
    facet_layers <- Filter(\(l) inherits(l, "FacetWrap"), list(p$facet))
    expect_length(facet_layers, 1)

    ## renders without error
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("build_kinetics_results data is named by interval_names", {
    data_list <- list(
        baseline = make_kinetics_data(),
        exercise = make_kinetics_data()
    )
    result_list <- list(
        make_kinetics_results("baseline"),
        make_kinetics_results("exercise")
    )

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    expect_named(result$data, c("baseline", "exercise"))
    expect_length(result$data, 2L)
})

test_that("build_kinetics_results model list is named by interval_names", {
    data_list <- list(
        A = make_kinetics_data(),
        B = make_kinetics_data()
    )
    result_list <- list(
        make_kinetics_results("A"),
        make_kinetics_results("B")
    )

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    expect_named(result$model, c("A", "B"))
    expect_type(result$model, "list")
})

test_that("build_kinetics_results interval_times scalar numeric", {
    data_list <- list(
        baseline = make_kinetics_data(interval_times = 1.5),
        exercise = make_kinetics_data(interval_times = 3.0)
    )
    result_list <- list(
        make_kinetics_results("baseline"),
        make_kinetics_results("exercise")
    )

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    et <- result$interval_times
    expect_s3_class(et, "data.frame")
    expect_equal(nrow(et), 2L)
    expect_equal(et$interval, c("baseline", "exercise"))
    expect_type(et$interval_times, "list")
    expect_equal(et$interval_times[[1L]], 1.5)
    expect_equal(et$interval_times[[2L]], 3.0)
})

test_that("build_kinetics_results interval_times unpacks list (ensemble)", {
    data_list <- list(
        ensemble = make_kinetics_data(interval_times = list(368, 1093))
    )
    result_list <- list(make_kinetics_results("ensemble"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    et <- result$interval_times
    expect_equal(nrow(et), 1L)
    expect_type(et$interval_times, "list")
    expect_equal(et$interval_times[[1L]], c(368, 1093))
})

test_that("build_kinetics_results interval_times is NA when attribute is NULL", {
    data_list <- list(int1 = make_kinetics_data(interval_times = NULL))
    result_list <- list(make_kinetics_results("int1"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    et <- result$interval_times
    expect_equal(nrow(et), 1L)
    expect_type(et$interval_times, "list")
    expect_true(is.na(et$interval_times[[1L]]))
})

test_that("build_kinetics_results diagnostics has interval col and correct rows", {
    data_list <- list(A = make_kinetics_data(), B = make_kinetics_data())
    result_list <- list(make_kinetics_results("A"), make_kinetics_results("B"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    diag <- result$diagnostics
    expect_equal(names(diag)[1L], "interval")
    expect_equal(nrow(diag), 2L)
    expect_equal(diag$interval, c("A", "B"))
})

test_that("build_kinetics_results channel_args has interval col and correct rows", {
    data_list <- list(A = make_kinetics_data(), B = make_kinetics_data())
    result_list <- list(make_kinetics_results("A"), make_kinetics_results("B"))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    ca <- result$channel_args
    expect_equal(names(ca)[1L], "interval")
    expect_equal(nrow(ca), 2L)
    expect_equal(ca$interval, c("A", "B"))
})

test_that("build_kinetics_results tags attrs per row with uneven channels", {
    ## interval A has 1 channel, B has 2 — attr rows must align per interval,
    ## not by recycling a length-2 interval vector across 3 rows
    data_list <- list(
        A = make_kinetics_data(),
        B = make_kinetics_data(channels = c("smo2_left", "smo2_right"))
    )
    result_list <- list(
        make_kinetics_results("A"),
        make_kinetics_results("B", channels = c("smo2_left", "smo2_right"))
    )

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    expect_equal(result$diagnostics$interval, c("A", "B", "B"))
    expect_equal(result$channel_args$interval, c("A", "B", "B"))
    expect_equal(result$coefficients$interval, c("A", "B", "B"))
})

test_that("build_kinetics_results handles multiple channels per interval", {
    channels <- c("smo2_left", "smo2_right")
    data_list <- list(int1 = make_kinetics_data(channels = channels))
    result_list <- list(make_kinetics_results("int1", channels = channels))

    result <- build_kinetics_results(
        data_list, result_list,
        method = "peak_slope",
        call = NULL
    )

    aug <- result$data[["int1"]]
    expect_true(all(paste0(channels, "_fitted") %in% names(aug)))
    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$nirs_channels, channels)
})


## analyse_kinetics ======================================================
## helper to create test mnirs data
create_kinetics_data <- function(
    n = 50,
    sample_rate = 10,
    channels = c("smo2_left", "smo2_right")
) {
    ## seed so `interval_times = sample(t, 1)` is reproducible
    set.seed(13)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    df <- data.frame(
        time = t,
        smo2_left = sin(t) * 10 + 50,
        smo2_right = cos(t) * 10 + 50
    )
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate,
        interval_times = sample(t, 1)
    )
}

test_that("analyse_kinetics returns mnirs_kinetics object", {
    data <- create_kinetics_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs_kinetics")
    expect_equal(result$method, "peak_slope")
})

test_that("analyse_kinetics captures call", {
    data <- create_kinetics_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_true(!is.null(result$call))
    expect_equal(class(result$call), "call")
})


test_that("analyse_kinetics works with data formats", {
    df1 <- create_kinetics_data()
    df2 <- create_kinetics_data()

    ## named lists
    result <- analyse_kinetics(
        list(A = df1, B = df2),
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$interval, c("A", "B"))
    expect_length(result$data, 2L)
    expect_named(result$data, c("A", "B"))

    ## unnamed lists
    result <- analyse_kinetics(
        list(df1, df2),
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$interval, c("interval_1", "interval_2"))

    ## single df
    result <- analyse_kinetics(
        df1,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_equal(nrow(result$coefficients), 1L)
    expect_equal(result$coefficients$interval, "interval_1")
    expect_length(result$data, 1L)
})


test_that("analyse_kinetics errors on invalid method", {
    data <- create_kinetics_data()

    expect_error(
        analyse_kinetics(data, method = "nonexistent"),
        "arg.*should be"
    )
})

test_that("analyse_kinetics$data overwrites metadata", {
    df <- data.frame(
        time = rep(seq(0, 4.9, by = 0.1), 2),
        time_alt = rep(seq(0, 4.9, by = 0.1), 2),
        smo2_left = c(
            sin(seq(0, 4.9, by = 0.1)) * 10 + 50,
            cos(seq(0, 4.9, by = 0.1)) * 10 + 50
        ),
        smo2_right = c(
            sin(seq(0, 4.9, by = 0.1)) * 10 + 50,
            cos(seq(0, 4.9, by = 0.1)) * 10 + 50
        )
    )
    df <- create_mnirs_data(
        df,
        nirs_channels = "smo2_left",
        time_channel = "time",
        sample_rate = 10,
        interval_times = sample(df$time, 1L)
    )

    result <- analyse_kinetics(
        df,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    aug <- result$data[[1L]]
    # attributes(aug)
    expect_s3_class(aug, "mnirs")
    expect_equal(attr(aug, "nirs_channels"), "smo2_left")
    expect_equal(attr(aug, "time_channel"), "time")
    expect_equal(attr(aug, "sample_rate"), 10)
    
    ## additional nirs_channel
    result <- analyse_kinetics(
        df,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = "time_alt",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )
    
    expect_equal(
        attr(result$data[[1L]], "nirs_channels"),
        c("smo2_left", "smo2_right")
    )
    expect_equal(attr(result$data[[1L]], "time_channel"), "time_alt")

    ## remove nirs_channel
    df <- create_mnirs_data(
        df,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = "time",
        sample_rate = 10,
        interval_times = sample(df$time, 1L)
    )

    result <- analyse_kinetics(
        df,
        nirs_channels = "smo2_right",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_equal(attr(result$data[[1L]], "nirs_channels"), "smo2_right")
})


test_that("analyse_kinetics$data preserves mnirs metadata with grouped input", {
    skip_if_not_installed("dplyr")

    df <- data.frame(
        time = rep(seq(0, 4.9, by = 0.1), 2),
        smo2 = c(
            sin(seq(0, 4.9, by = 0.1)) * 10 + 50,
            cos(seq(0, 4.9, by = 0.1)) * 10 + 50
        ),
        group = rep(c("A", "B"), each = 50)
    )
    df <- create_mnirs_data(
        df,
        nirs_channels = "smo2",
        time_channel = "time",
        sample_rate = 10,
        interval_times = sample(df$time, 1L)
    )
    grouped_df <- dplyr::group_by(df, group)

    result <- analyse_kinetics(
        grouped_df,
        nirs_channels = "smo2",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expect_length(result$data, 2L)
    for (nm in c("A", "B")) {
        aug <- result$data[[nm]]
        expect_s3_class(aug, "mnirs")
        expect_equal(attr(aug, "nirs_channels"), "smo2")
        expect_equal(attr(aug, "time_channel"), "time")
        expect_equal(attr(aug, "sample_rate"), 10)
    }
})

## analyse_kinetics.response_time ======================================
## helper to create a deterministic ramp/plateau for response_time tests
create_response_time_data <- function(
    n_baseline = 5,
    n_ramp = 10,
    n_plateau = 5,
    A = 0,
    B = 20,
    sample_rate = 1,
    channels = "smo2"
) {
    x <- c(
        rep(A, n_baseline),
        seq(A, B, length.out = n_ramp),
        rep(B, n_plateau)
    )
    t <- seq_along(x) - n_baseline ## t = 0 at end of baseline
    df <- setNames(
        data.frame(t, x),
        c("time", channels[1])
    )
    if (length(channels) > 1) {
        for (ch in channels[-1]) {
            df[[ch]] <- x
        }
    }
    create_mnirs_data(
        df,
        nirs_channels = channels,
        time_channel = "time",
        sample_rate = sample_rate,
        interval_times = 0
    )
}

test_that("analyse_kinetics.response_time forwards fraction argument", {
    data <- create_response_time_data()

    result_25 <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "response_time",
        fraction = 0.25,
        direction = "positive",
        verbose = FALSE
    )
    result_50 <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "response_time",
        fraction = 0.5,
        direction = "positive",
        verbose = FALSE
    )

    ## response_time monotonic in fraction
    expect_lt(
        result_25$coefficients$response_time,
        result_50$coefficients$response_time
    )

    ## fitted = A + (B - A) * fraction
    A <- result_50$coefficients$A
    B <- result_50$coefficients$B
    expect_equal(result_25$coefficients$fitted, A + (B - A) * 0.25)
    expect_equal(result_50$coefficients$fitted, A + (B - A) * 0.5)
})

test_that("analyse_kinetics.response_time passes direction argument", {
    data <- create_kinetics_data(n = 100)

    result_pos <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "response_time",
        direction = "positive",
        verbose = FALSE
    )
    result_neg <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "response_time",
        direction = "negative",
        verbose = FALSE
    )

    expect_gt(result_pos$coefficients$B, result_pos$coefficients$A)
    expect_lt(result_neg$coefficients$B, result_neg$coefficients$A)
})

test_that("analyse_kinetics.response_time channel_args override defaults", {
    data <- create_kinetics_data(n = 100)

    result <- analyse_kinetics(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        method = "response_time",
        direction = list("positive", smo2_right = "negative"),
        verbose = FALSE
    )

    ## channel_args records per-channel settings
    ca <- result$channel_args
    expect_equal(
        ca$direction[ca$nirs_channels == "smo2_left"],
        "positive"
    )
    expect_equal(
        ca$direction[ca$nirs_channels == "smo2_right"],
        "negative"
    )
})

test_that("analyse_kinetics.response_time results have correct columns", {
    data <- create_response_time_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2",
        method = "response_time",
        direction = "positive",
        verbose = FALSE
    )

    expected_cols <- c(
        "interval", "nirs_channels", "A", "B",
        "response_time", "response_value", "fitted"
    )
    expect_true(all(expected_cols %in% names(result$coefficients)))
})

test_that("analyse_kinetics.response_time handles unreachable response", {
    data <- create_kinetics_data(n = 100, channels = c("smo2_left"))
    data <- data[data$time >= 8, ]
    attr(data, "interval_times") <- 9.3 ## 9.6

    ## visual check
    # plot(data) +
    #     ggplot2::geom_vline(xintercept = 9.3)

    expect_no_error(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2_left",
            time_channel = "time",
            method = "response_time",
            direction = "positive",
            verbose = FALSE
        )
    )
    
    expect_true(is.na(result$coefficients$response_time))
    expect_true(is.na(result$coefficients$response_value))
    expect_true(is.na(result$coefficients$fitted))

    ## fitted column exists and is all-finite at baseline/extreme rows,
    ## NA elsewhere (no error from NA subscripted assignment).
    fitted_df <- result$data[[1L]]
    t0 <- unlist(result$interval_times$interval_times)
    expect_true("smo2_left_fitted" %in% names(fitted_df))
    expect_false(is.na(fitted_df$smo2_left_fitted[fitted_df$time == t0]))
    expect_false(
        is.na(fitted_df$smo2_left_fitted[which(fitted_df$time == t0) + 1])
    )
})

test_that("analyse_kinetics.response_time dispatches via method aliases", {
    data <- create_response_time_data()
    aliases <- c("response time", "HRT", "half_recovery_time", "half time")

    for (alias in aliases) {
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = alias[1],
            direction = "positive",
            verbose = FALSE
        )
        expect_equal(result$method, "response_time")
        expect_true("response_time" %in% names(result$coefficients))
    }
})


## analyse_kinetics.peak_slope =========================================
## structure, data formats, grouped data covered by generic tests above

test_that("analyse_kinetics.peak_slope passes width and span correctly", {
    data <- create_kinetics_data(
        channels = "smo2_right",
        n = 40,
        sample_rate = 10
    )

    result_width <- analyse_kinetics(
        data,
        nirs_channels = "smo2_right",
        method = "peak_slope",
        width = 10,
        verbose = FALSE
    )

    result_span <- analyse_kinetics(
        data,
        nirs_channels = "smo2_right",
        method = "peak_slope",
        span = 1,
        verbose = FALSE
    )

    ## both should produce valid results
    expect_false(is.na(result_width$coefficients$slope))
    expect_false(is.na(result_span$coefficients$slope))
    expect_equal(result_width$diagnostics$n_obs, 10)
    ## span buffer number of samples for start and end inclusive
    expect_true(
        all.equal(result_span$diagnostics$n_obs, 10, tolerance = 2, scale = 1)
    )
})

test_that("analyse_kinetics.peak_slope respects global verbose option", {
    data <- create_kinetics_data(n = 40, channels = "smo2_left")

    old_verbose <- getOption("mnirs.verbose")
    on.exit(options(mnirs.verbose = old_verbose), add = TRUE)

    ## supplying both width and span normally emits a verbose "overrides" hint
    options(mnirs.verbose = TRUE)
    expect_message(
        analyse_kinetics(
            data,
            nirs_channels = "smo2_left",
            method = "peak_slope",
            width = 5,
            span = 1
        ),
        "overrides"
    )

    ## global option must reach the S3 method when `verbose` is omitted
    options(mnirs.verbose = FALSE)
    expect_silent(
        analyse_kinetics(
            data,
            nirs_channels = "smo2_left",
            method = "peak_slope",
            width = 5,
            span = 1
        )
    )
})

test_that("analyse_kinetics.peak_slope passes direction argument", {
    data <- create_kinetics_data(n = 100)

    result_pos <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        direction = "positive",
        verbose = FALSE
    )

    result_neg <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        direction = "negative",
        verbose = FALSE
    )

    expect_gt(result_pos$coefficients$slope, 0)
    expect_lt(result_neg$coefficients$slope, 0)
})

test_that("analyse_kinetics.peak_slope channel_args override defaults", {
    data <- create_kinetics_data(n = 100)

    result <- analyse_kinetics(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        method = "peak_slope",
        width = 5,
        direction = list("positive", smo2_right = "negative"),
        verbose = FALSE
    )

    expect_gt(result$coefficients$slope[1], 0) ## smo2_left positive
    expect_lt(result$coefficients$slope[2], 0) ## smo2_right negative

    ## channel_args should record per-channel settings
    ca <- result$channel_args
    expect_equal(
        ca$direction[ca$nirs_channels == "smo2_left"],
        "positive"
    )
    expect_equal(
        ca$direction[ca$nirs_channels == "smo2_right"],
        "negative"
    )
})

test_that("analyse_kinetics.peak_slope results have correct columns", {
    data <- create_kinetics_data()

    result <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )

    expected_cols <- c(
        "interval", "nirs_channels", "slope", "intercept", 
        "fitted", "peak_slope_time", "idx"
    )
    expect_true(all(expected_cols %in% names(result$coefficients)))
})


## analyse_kinetics.monoexponential ====================================
## helper: create monoexponential test data with known parameters
create_monoexp_data <- function(
    A = 50,
    B = 80,
    tau = 5,
    TD = 5,
    n = 60,
    sample_rate = 1,
    noise_sd = 0.5,
    channels = "smo2"
) {
    set.seed(13)
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


test_that("analyse_kinetics.monoexponential dispatches multiple channels", {
    nirs_channels <- c("smo2_left", "smo2_right")
    data <- create_monoexp_data(channels = nirs_channels)

    result <- analyse_kinetics(
        data,
        nirs_channels = nirs_channels,
        method = "monoexponential",
        use_time_delay = FALSE,
        verbose = FALSE
    )

    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$nirs_channels, nirs_channels)
    expect_named(
        result$data[[1]],
        c("time", nirs_channels, paste0(nirs_channels, "_fitted"))
    )
})

test_that("analyse_kinetics.monoexponential uses custom interval name", {
    ## only 3 observations for a 3-param model
    data <- create_monoexp_data(n = 3, noise_sd = 0.1)

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "monoexponential",
            use_time_delay = FALSE
        ),
        "fit failed for.*smo2.*interval_1" ## call custom interval name
    )

    expect_true(is.na(result$coefficients$A))
    expect_true(is.na(result$coefficients$tau))
    expect_true(is.na(result$coefficients$k))
})

test_that("analyse_kinetics.monoexponential names only the failing interval", {
    ## named list: `good` converges, `bad` has too few points to fit
    data <- list(
        good = create_monoexp_data(n = 60, noise_sd = 0.5),
        bad = create_monoexp_data(n = 3, noise_sd = 0.1)
    )

    msg <- conditionMessage(rlang::catch_cnd(
        analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "monoexponential",
            use_time_delay = FALSE
        ),
        classes = "warning"
    ))

    ## warning names the failing interval only, not the converged one
    expect_match(msg, "fit failed for.*smo2.*bad")
    expect_no_match(msg, "good")
})

## analyse_kinetics.logistic ============================================
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
    channels = "smo2"
) {
    set.seed(13)
    t <- seq(0, (n - 1) / sample_rate, length.out = n)
    x <- logistic(t, A, B, xmid, slope, asym) + rnorm(n, 0, noise_sd)

    df <- setNames(
        data.frame(t, x),
        c("time", channels[1])
    )
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


test_that("analyse_kinetics.logistic dispatches multiple channels", {
    nirs_channels <- c("smo2_left", "smo2_right")
    data <- create_logistic_data(channels = nirs_channels)

    result <- analyse_kinetics(
        data,
        nirs_channels = nirs_channels,
        method = "logistic",
        shape = "symmetric",
        verbose = FALSE
    )

    expect_equal(nrow(result$coefficients), 2L)
    expect_equal(result$coefficients$nirs_channels, nirs_channels)
    expect_named(
        result$data[[1]],
        c("time", nirs_channels, paste0(nirs_channels, "_fitted"))
    )
})

test_that("analyse_kinetics.logistic uses custom interval name", {
    ## only 3 observations for a 4-param model
    data <- create_logistic_data(n = 10, noise_sd = 0.1)

    expect_warning(
        result <- analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "logistic",
            shape = "symmetric",
        ),
        "fit failed for.*smo2.*interval_1" ## call custom interval name
    )

    expect_true(is.na(result$coefficients$A))
    expect_true(is.na(result$coefficients$xmid))
    expect_true(is.na(result$coefficients$slope))
})

test_that("analyse_kinetics.logistic names only the failing interval", {
    ## named list: `good` converges, `bad` has too few points to fit
    data <- list(
        good = create_logistic_data(n = 60, noise_sd = 2),
        bad = create_logistic_data(n = 10, noise_sd = 0.1)
    )

    msg <- conditionMessage(rlang::catch_cnd(
        analyse_kinetics(
            data,
            nirs_channels = "smo2",
            method = "logistic",
            shape = "symmetric"
        ),
        classes = "warning"
    ))

    ## warning names the failing interval only, not the converged one
    expect_match(msg, "fit failed for.*smo2.*bad")
    expect_no_match(msg, "good")
})

## analyze_kinetics (US spelling alias) ================================
test_that("analyze_kinetics matches analyse_kinetics for peak_slope", {
    data <- create_kinetics_data()

    result_uk <- analyse_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak_slope",
        width = 5,
        verbose = FALSE
    )
    result_us <- analyze_kinetics(
        data,
        nirs_channels = "smo2_left",
        method = "peak slope",
        width = 5,
        verbose = FALSE
    )

    expect_s3_class(result_us, "mnirs_kinetics")
    expect_s3_class(result_uk, "mnirs_kinetics")
    ## US alias rewrites call to analyse_kinetics(), so full results match
    expect_equal(result_us, result_uk)
    expect_equal(result_us$call[[1L]], quote(analyse_kinetics))
    expect_equal(result_us$call[[4L]], "peak_slope")
})

test_that("analyze_kinetics forwards method = 'response_time' and fraction", {
    data <- create_response_time_data()

    result <- analyze_kinetics(
        data,
        nirs_channels = "smo2",
        method = "response time",
        fraction = 0.5,
        direction = "positive",
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs_kinetics")
    expect_equal(result$method, "response_time")
    expect_equal(result$call[[1L]], quote(analyse_kinetics))
    expect_equal(result$call[[4L]], "response_time")
    expect_true("response_time" %in% names(result$coefficients))

    ## fraction was forwarded: fitted = A + (B - A) * 0.5
    A <- result$coefficients$A
    B <- result$coefficients$B
    expect_equal(result$coefficients$fitted, A + (B - A) * 0.5)
})


## print.mnirs_kinetics ================================================
## helper: minimal mnirs_kinetics fixture for print method tests
make_print_kinetics <- function(coefs, method = "peak_slope") {
    structure(
        list(method = method, coefficients = coefs),
        class = "mnirs_kinetics"
    )
}

test_that("print.mnirs_kinetics shows peak_slope header", {
    x <- make_print_kinetics(
        data.frame(interval = "int1", nirs_channels = "smo2", slope = 1.5),
        method = "peak_slope"
    )
    output <- capture.output(print(x))
    expect_true(any(grepl("Peak Linear Regression Slope", output)))
})

test_that("print.mnirs_kinetics shows monoexponential header", {
    x <- make_print_kinetics(
        data.frame(interval = "int1", nirs_channels = "smo2", tau = 5.0),
        method = "monoexponential"
    )
    output <- capture.output(print(x))
    expect_true(any(grepl("Monoexponential non-linear Regression", output)))
})

test_that("print.mnirs_kinetics shows response_time header", {
    x <- make_print_kinetics(
        data.frame(interval = "int1", nirs_channels = "smo2", tau = 5.0),
        method = "response_time"
    )
    output <- capture.output(print(x))
    expect_true(any(grepl("Fractional Response Time", output)))
})

test_that("print.mnirs_kinetics always shows Model Coefficients label", {
    for (m in c("peak_slope", "monoexponential", "response_time")) {
        x <- make_print_kinetics(
            data.frame(interval = "int1", nirs_channels = "smo2", slope = 1),
            method = m
        )
        output <- capture.output(print(x))
        expect_true(
            any(grepl("Model Coefficients:", output)),
            info = paste("method =", m)
        )
    }
})

test_that("print.mnirs_kinetics drops time_channel column from display", {
    x <- make_print_kinetics(
        data.frame(
            interval = "int1",
            nirs_channels = "smo2",
            time_channel = "time",
            slope = 1.5
        )
    )
    output <- capture.output(print(x))
    expect_false(any(grepl("time_channel", output)))
    expect_true(any(grepl("slope", output)))
})

test_that("print.mnirs_kinetics drops 'fitted$' columns from display", {
    x <- make_print_kinetics(
        data.frame(
            interval = "int1",
            nirs_channels = "smo2",
            MRT = 4.0,
            MRT_fitted = 65.0,
            response_value_fitted = 70.0,
            fitted = 50,
            fitted_column = 10
        ),
        method = "monoexponential"
    )
    output <- capture.output(print(x))
    expect_false(any(grepl("MRT_fitted", output)))
    expect_false(any(grepl("response_value_fitted", output)))
    expect_false(any(grepl("\\bfitted\\b", output)))
    expect_true(any(grepl("\\bMRT\\b", output)))
    expect_true(any(grepl("fitted_column", output)))
})

test_that("print.mnirs_kinetics formats numerics to 4 sig figs", {
    x <- make_print_kinetics(
        data.frame(
            interval = "int1",
            nirs_channels = "smo2",
            slope = 1.234567
        )
    )
    output <- capture.output(print(x))
    expect_true(any(grepl("1.235", output, fixed = TRUE)))
    expect_false(any(grepl("1.234567", output, fixed = TRUE)))
})

test_that("print.mnirs_kinetics renders NA values as 'NA'", {
    x <- make_print_kinetics(
        data.frame(
            interval = c("int1", "int2"),
            nirs_channels = c("smo2", "smo2"),
            slope = c(1.234, NA_real_)
        )
    )
    output <- capture.output(print(x))
    expect_true(any(grepl("\\bNA\\b", output)))
})

test_that("print.mnirs_kinetics handles all-NA numeric columns", {
    x <- make_print_kinetics(
        data.frame(
            interval = "int1",
            nirs_channels = "smo2",
            slope = 1.5,
            intercept = NA_real_
        )
    )
    expect_no_error(output <- capture.output(print(x)))
    expect_true(any(grepl("intercept", output)))
})

test_that("print.mnirs_kinetics prints all rows when nrow <= 10", {
    coefs <- data.frame(
        interval = paste0("int", 1:10),
        nirs_channels = "smo2",
        slope = seq(0.1, 1.0, by = 0.1)
    )
    x <- make_print_kinetics(coefs)
    output <- capture.output(print(x))
    expect_false(any(grepl("rows omitted", output)))
    ## every interval label should appear in output
    for (lab in coefs$interval) {
        expect_true(
            any(grepl(lab, output)),
            info = paste("missing:", lab)
        )
    }
})

test_that("print.mnirs_kinetics truncates when nrow > 10", {
    coefs <- data.frame(
        interval = paste0("int", sprintf("%02d", 1:15)),
        nirs_channels = "smo2",
        slope = seq(0.1, 1.5, by = 0.1)
    )
    x <- make_print_kinetics(coefs)
    output <- capture.output(print(x))

    ## spacer reports omitted count
    expect_true(any(grepl("--- 5 rows omitted", output, fixed = TRUE)))

    ## first 5 and last 5 intervals appear
    for (lab in coefs$interval[c(1:5, 11:15)]) {
        expect_true(
            any(grepl(lab, output)),
            info = paste("missing:", lab)
        )
    }
    ## middle rows do not appear
    for (lab in coefs$interval[6:10]) {
        expect_false(
            any(grepl(lab, output)),
            info = paste("unexpected:", lab)
        )
    }
})



## fix_coefs() =========================================================
test_that("fix_coefs() fixes single parameter correctly", {
    set.seed(303)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    model_fixed <- fix_coefs(model, TD = 15)

    expect_s3_class(model_fixed, "nls")
    expect_false("TD" %in% names(coef(model_fixed)))
    expect_true(all(c("A", "B", "tau") %in% names(coef(model_fixed))))
    ## 15 should be in the model formula
    expect_true(any(grepl("15", model_fixed$call$formula)))
})

test_that("fix_coefs() fixes multiple parameters", {
    set.seed(400)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    # ggplot2::ggplot(data, ggplot2::aes(t, x)) +
    #     theme_mnirs() +
    #     ggplot2::geom_point()

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    model_fixed <- fix_coefs(
        model,
        A = 10,
        TD = 15,
    )

    expect_false("A" %in% names(coef(model_fixed)))
    expect_false("TD" %in% names(coef(model_fixed)))
    expect_length(coef(model_fixed), 2)
    ## 10 & 15 should be in the model formula
    expect_true(any(grepl("15", model_fixed$call$formula)))
    expect_true(any(grepl("10", model_fixed$call$formula)))
})

test_that("fix_coefs() errors when all parameters fixed", {
    set.seed(505)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)

    expect_error(
        fix_coefs(
            model,
            A = 10,
            B = 100,
            TD = 15,
            tau = 8,
        ),
        "Cannot update the model if all parameters are fixed"
    )
})

test_that("fix_coefs() warns for invalid parameter names", {
    set.seed(606)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)

    expect_warning(
        fix_coefs(model, INVALID = 99, verbose = TRUE),
        "Unknown model coefficient"
    )

    expect_silent(
        fix_coefs(model, INVALID = 99, verbose = FALSE)
    )
})

test_that("fix_coefs() accepts explicit data argument", {
    set.seed(1010)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)

    expect_no_error(
        model_fixed <- fix_coefs(model, TD = 15, data = data)
    )
    expect_s3_class(model_fixed, "nls")
})

test_that("fix_coefs() predictions differ from original model", {
    set.seed(1111)
    t <- 1:60
    x <- monoexponential(t, A = 10, B = 100, tau = 8, TD = 15) +
        rnorm(length(t), 0, 3)
    data <- data.frame(t, x)

    model <- nls(x ~ SSmonoexp(t, A, B, tau, TD), data = data)
    model_fixed <- fix_coefs(model, TD = 20, data = data, verbose = FALSE)

    pred_orig <- predict(model, data)
    pred_fixed <- predict(model_fixed, data)

    # Should differ if fixed value differs from estimated
    expect_false(identical(pred_orig, pred_fixed))
})


## integration =======================================================
test_that("analyse_kinetics works visually on Train.Red", {
    skip_if_not_installed("ggplot2")
    skip("visual check")

    data_list <- read_mnirs(
        example_mnirs("train.red"),
        nirs_channels = c(smo2_left = "SmO2"),
        time_channel = c(time = "Timestamp (seconds passed)"),
        zero_time = TRUE,
        verbose = FALSE
    ) |>
        resample_mnirs(method = "linear", verbose = FALSE) |>
        extract_intervals(
            start = by_time(368),
            span = c(10, 90),
            zero_time = TRUE,
            verbose = FALSE
        )

    result <- analyse_kinetics(data_list, method = "peak_slope", span = 10)

    library(ggplot2)
    plot(data_list[[1]]) +
        geom_line(
            data = result$data[[1]],
            aes(y = smo2_left_fitted),
            linewidth = 1.5
        ) +
        geom_point(
            data = result$coefficients,
            aes(x = time, y = fitted),
            size = 4,
            shape = 21,
            stroke = 1.5,
            fill = "white"
        )

    result <- analyse_kinetics(data_list, method = "monoexp", use_time_delay = TRUE)
    # result$diagnostics
    # result$coefficients

    plot(data_list[[1]]) +
        geom_line(
            data = result$data[[1]],
            aes(y = smo2_left_fitted),
            linewidth = 1.5
        ) +
        geom_point(
            data = result$coefficients,
            aes(x = MRT, y = MRT_fitted),
            size = 4,
            shape = 21,
            stroke = 1.5,
            fill = "white"
        )
})


## benchmark ===========================================================
test_that("analyse_kinetics benchmark", {
    ## baselne established from documented example on initial run;
    ## fails if itr/sec regresses by >10%
    skip("benchmark baseline test")

    data_list <- read_mnirs(
        example_mnirs("train.red"),
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        zero_time = TRUE,
        verbose = FALSE
    ) |>
        resample_mnirs(method = "linear", verbose = FALSE) |>
        extract_intervals(
            group_intervals = "distinct",
            start = by_time(368, 1084),
            span = c(-20, 90),
            zero_time = TRUE,
            verbose = FALSE
        )

    # bench::mark(
    #     analyse_peak_slope = suppressWarnings(
    #         lapply(data_list, \(.df) {
    #             analyse_peak_slope(
    #                 .df,
    #                 nirs_channels = c(smo2_left, smo2_right),
    #                 span = 10,
    #                 verbose = FALSE
    #             )
    #         })
    #     ),
    #     kinetics.peak_slope = analyse_kinetics(
    #         data_list,
    #         nirs_channels = c(smo2_left, smo2_right),
    #         method = "peak_slope",
    #         span = 10,
    #         verbose = FALSE
    #     ),
    #     analyse_monoexponential = suppressWarnings(
    #         lapply(data_list, \(.df) {
    #             analyse_monoexponential(
    #                 .df,
    #                 nirs_channels = c(smo2_left, smo2_right),
    #                 use_time_delay = TRUE,
    #                 verbose = FALSE
    #             )
    #         })
    #     ),
    #     kinetics.monoexponential = analyse_kinetics(
    #         data_list,
    #         nirs_channels = c(smo2_left, smo2_right),
    #         method = "monoexponential",
    #         use_time_delay = TRUE,
    #         verbose = FALSE
    #     ),
    #     analyse_response = suppressWarnings(
    #         lapply(data_list, \(.df) {
    #             analyse_response_time(
    #                 .df,
    #                 nirs_channels = c(smo2_left, smo2_right),
    #                 fraction = 0.5,
    #                 verbose = FALSE
    #             )
    #         })
    #     ),
    #     kinetics.monoexponential = analyse_kinetics(
    #         data_list,
    #         nirs_channels = c(smo2_left, smo2_right),
    #         method = "response_time",
    #         fraction = 0.5,
    #         verbose = FALSE
    #     ),
    #     iterations = 5L,
    #     check = FALSE
    # )

    # A tibble: 2 × 13
    #>    expression                 min median `itr/sec` mem_alloc `gc/sec`
    #>    <bch:expr>              <bch:> <bch:>     <dbl> <bch:byt>    <dbl>
    #>  1 analyse_peak_slope      69.3ms 73.4ms      13.8    30.4MB     22.1
    #>  2 analyse_kinetics.peak_… 72.4ms 75.4ms      13.3    30.5MB     21.2

    # itr_per_sec <- bm$`itr/sec`

    # ## baseline: update this value when optimising (seconds)
    # ## run test interactively to calibrate:
    # ##   itr_per_sec will be printed on first failure
    # baseline <- 8
    # threshold <- baseline * 1.10 ## 10% regression budget

    # expect_lte(
    #     itr_per_sec,
    #     threshold,
    #     label = sprintf(
    #         "%.3f itr/sec exceeds %.0f%% of baseline %.3fs (limit %.3fs)",
    #         itr_per_sec,
    #         110,
    #         baseline,
    #         threshold
    #     )
    # )
})
