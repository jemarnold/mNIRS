## response_time ====================================================
test_that("response_time returns correct structure", {
    x <- c(rep(10, 5), seq(10, 60, length.out = 15), rep(60, 5))
    t <- seq_along(x)

    result <- response_time(x, t, start_time = 5, response_fraction = 0.5)

    expect_type(result, "list")
    expect_named( result, c(
        "A", "B", "response_time", "response_value",
        "fitted", "baseline_idx", "response_idx", "extreme_idx"
    ))
    expect_all_equal(lengths(result)[names(result) != "baseline_idx"], 1)
    expect_equal(lengths(result)["baseline_idx"], 5, ignore_attr = TRUE)
    expect_type(result$A, "double")
    expect_type(result$B, "double")
    expect_type(result$response_time, "double")
    expect_type(result$response_value, "double")
    expect_type(result$fitted, "double")
    expect_type(result$baseline_idx, "integer")
    expect_type(result$response_idx, "integer")
    expect_type(result$extreme_idx, "integer")
})

test_that("response_time computes correct values", {
    ## baseline idx 0-5, ramp from 6 to 16 plateau 16-20
    x <- c(rep(0, 5), seq(0, 20, length.out = 10), rep(20, 5))
    t <- seq_along(x) - 1  ## t = 0 at index 1
    result <- response_time(x, t, start_time = 0, response_fraction = 0.5)

    ## A = mean baseline (t <= 0 = first element = 0)
    expect_equal(result$A, 0)
    ## B = max value = 20
    expect_equal(result$B, 20)
    ## response time = half of 18
    expect_equal(result$response_time, 10)
    ## response_idx time+1 since time = 0 at idx = 1
    expect_equal(result$response_idx, 11)
    ## fitted = A + response_fraction * (B - A) = 10
    expect_gte(result$response_value, 10) ## 11.1 real x value
    expect_equal(result$fitted, 10)
    ## extreme_idx first plateau idx
    expect_equal(result$extreme_idx, 15)
})

test_that("response_time validates inputs", {
    expect_error(
        response_time(1:5, response_fraction = 1.5),
        "response_fraction.*valid.*numeric"
    )
    expect_error(
        response_time(1:5, response_fraction = -0.1),
        "response_fraction.*valid.*numeric"
    )
    expect_silent(response_time(1:5, t = 0:4, end_window = Inf))
    expect_error(
        response_time(x = "a"),
        "x.*valid.*numeric"
    )
    expect_error(
        response_time(1:5, t = 1:3),
        "equal length"
    )
})

test_that("response_time respects `bypass_checks`", {
    ## expect NOT warning "No observations.*t.*<=.*start_time"
    expect_warning(
        response_time(
            1:5, start_time = -5, bypass_checks = TRUE
        ),
        "No valid.*auto"
    )
    ## expect NOT error "No observations in.*t.*before.*sstart_time"
    expect_warning(
        response_time(
            1:5, start_time = 10, bypass_checks = TRUE
        ),
        "No valid.*auto"
    )
    
})

test_that("response_time respects response_fraction parameter", {
    x <- c(rep(0, 5), seq(0, 100, length.out = 20))
    t <- seq_along(x) - 1

    result_50 <- response_time(x, t, response_fraction = 0.5)
    result_25 <- response_time(x, t, response_fraction = 0.25)
    result_75 <- response_time(x, t, response_fraction = 0.75)

    ## 50% response time should fall between 25% and 75%
    expect_lt(result_25$response_time, result_50$response_time)
    expect_lt(result_50$response_time, result_75$response_time)

    ## fitted values reflect response_fraction
    expect_equal(result_50$fitted, 50)
    expect_equal(result_25$fitted, 25)
    expect_equal(result_75$fitted, 75)
})

test_that("response_time propagates auto-detected direction", {
    ## positive signal
    x_pos <- c(rep(0, 5), seq(0, 20, length.out = 10))
    result_pos <- response_time(x_pos, verbose = FALSE)
    expect_gt(result_pos$B, result_pos$A)

    ## negative signal
    x_neg <- c(rep(0, 5), seq(0, -20, length.out = 10))
    result_neg <- response_time(x_neg, verbose = FALSE)
    expect_lt(result_neg$B, result_neg$A)
})

test_that("response_time respects manual direction", {
    ## strongly positive signal, forced negative → uses which.min
    x <- c(rep(0, 5), seq(0, 20, length.out = 10), seq(20, 15, length.out = 5))
    t <- seq_along(x)

    result_pos <- response_time(
        x, t, start_time = 5, direction = "positive", verbose = FALSE
    )
    result_neg <- response_time(
        x, t, start_time = 5, direction = "negative", verbose = FALSE
    )

    ## positive: extreme = max; negative: extreme = min
    expect_equal(result_pos$B, max(x))
    expect_equal(result_neg$B, min(x))
    expect_equal(result_neg$response_time, 1) ## start_time + 1L idx
    expect_equal(result_neg$response_value, min(x))
    expect_equal(result_neg$fitted, min(x))
    expect_equal(result_neg$response_idx, 6) ## start_time + 1
})

test_that("response_time warns when no baseline observations", {
    x <- seq(0, 20, length.out = 20)
    t <- seq_along(x)

    ## start_time is below min(t), so no baseline
    expect_warning(
        result <- response_time(x, t, start_time = 0),
        "No observations where `t`"
    )
    ## falls back to x[1] as baseline
    expect_equal(result$baseline_idx, 1L)
})

test_that("response_time warns when extreme precedes start_time", {
    ## peak at index 3, but start_time = 5 → baseline spans indices 1:5 including peak
    x <- c(0, 5, 20, 10, 5, 1, 1, 1, 1, 1)
    t <- seq_along(x)

    expect_warning(
        result <- response_time(x, t, start_time = 5, direction = "positive"),
        "No valid.*extremes"
    )
    
    expect_true(is.na(result$response_time))
    expect_true(is.na(result$response_value))
    expect_true(is.na(result$response_idx))
    expect_true(is.na(result$fitted))
})

test_that("response_time errors when start_time exceeds max(t)", {
    x <- 1:10
    t <- seq_along(x)

    expect_error(
        response_time(x, t, start_time = 11),
        "No observations.*before"
    )
})


test_that("response_time response_time is relative to start_time", {
    x <- c(rep(10, 10), seq(10, 50, length.out = 20))
    t <- seq(0, 290, by = 10)  ## 0:290 seconds

    result_0 <- response_time(x, t, start_time = 0, verbose = FALSE)
    result_100 <- response_time(x, t, start_time = 100, verbose = FALSE)

    ## response_time = t[response_idx] - start_time
    expect_equal(result_0$response_time, t[result_0$response_idx] - 0)
    expect_equal(result_100$response_time, t[result_100$response_idx] - 100)
})

test_that("response_time handles NA in x", {
    ## NAs in signal don't prevent result
    x <- c(rep(0, 5), NA, seq(0, 20, length.out = 10), NA)
    t <- seq_along(x)

    result <- response_time(x, t, start_time = 5, verbose = FALSE)

    expect_false(is.na(result$A))
    expect_false(is.na(result$B))
})


## analyse_response_time =================================================
test_that("analyse_response_time returns correct structure", {
    x <- c(rep(10, 5), seq(10, 60, length.out = 15), rep(60, 5))
    q <- c(rep(20, 5), seq(20, 80, length.out = 15), rep(80, 5))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    result <- analyse_response_time(df, start_time = 5, verbose = FALSE)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
    expect_named(result, c(
        "interval", "nirs_channels", "response_fraction", "A", "B", 
        "response_time", "response_value", "fitted", "idx"
    ))

    expect_type(result$nirs_channels, "character")
    expect_equal(result$nirs_channels, c("x", "q"))
    expect_type(result$A, "double")
    expect_type(result$B, "double")
    expect_type(result$response_time, "double")
    expect_type(result$response_value, "double")
    expect_type(result$fitted, "double")
    expect_type(result$idx, "integer")

    ## attributes
    expect_type(attr(result, "model"), "list")
    expect_null(attr(result, "model")$x)  ## no model for response_time method
    expect_type(attr(result, "fitted_data"), "list")
    expect_s3_class(attr(result, "fitted_data")$x, "data.frame")
    expect_equal(
        colnames(attr(result, "fitted_data")$x),
        c("window_idx", "fitted")
    )
    expect_s3_class(attr(result, "channel_args"), "data.frame")
    expect_s3_class(attr(result, "diagnostics"), "data.frame")
})

test_that("analyse_response_time computes correct values", {
    ## known signal: baseline = 0 (t <= 0), ramp to 100 over 100s
    x <- c(rep(0, 10), seq(0, 100, length.out = 91))
    t <- seq(0, 100, by = 1)

    df <- create_mnirs_data(
        data.frame(t, x),
        nirs_channels = "x",
        time_channel = "t"
    )
    result <- analyse_response_time(df, start_time = 0, response_fraction = 0.5)

    ## A = baseline mean ~ 0
    expect_equal(result$A, 0)
    ## B = max ~ 100
    expect_equal(result$B, 100)
    ## fitted = 50
    expect_equal(result$fitted, 50)
    ## response_time at half point + 5 unit baseline
    expect_equal(result$response_time, 55)
    ## fitted = half point
    expect_equal(result$fitted, 50)
    expect_equal(result$response_value, result$fitted)
})

test_that("analyse_response_time validates data structure", {
    expect_error(analyse_response_time(1:5), "data frame")

    df <- data.frame(t = 1:10, x = 1:10)
    expect_error(analyse_response_time(df), "nirs_channels")

    df_mnirs <- create_mnirs_data(
        data.frame(t = 1:10, x = 1:10),
        nirs_channels = "x",
        time_channel = "t"
    )
    attr(df_mnirs, "time_channel") <- NULL
    expect_error(
        analyse_response_time(df_mnirs, verbose = FALSE),
        "time_channel"
    )
})

test_that("analyse_response_time processes multiple channels independently", {
    ## x: strong increase; q: strong decrease
    x <- c(rep(0, 5), seq(0, 50, length.out = 20))
    q <- c(rep(50, 5), seq(50, 0, length.out = 20))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    result <- analyse_response_time(df, start_time = 5, verbose = FALSE)

    ## x: B > A; q: B < A
    expect_gt(
        result$B[result$nirs_channels == "x"],
        result$A[result$nirs_channels == "x"]
    )
    expect_lt(
        result$B[result$nirs_channels == "q"],
        result$A[result$nirs_channels == "q"]
    )
})

test_that("analyse_response_time channel_args override defaults", {
    x <- c(rep(0, 5), seq(0, 50, length.out = 20))
    q <- c(rep(0, 5), seq(0, 50, length.out = 20))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    result <- analyse_response_time(
        df,
        start_time = 5,
        response_fraction = list(0.5, q = 0.25),
        verbose = FALSE
    )

    ## x: response_fraction = 0.5; q: response_fraction = 0.25 → q response time < x response time
    expect_lt(
        result$response_time[result$nirs_channels == "q"],
        result$response_time[result$nirs_channels == "x"]
    )
    expect_equal(attr(result, "channel_args")$response_fraction, c(0.5, 0.25))
})

test_that("vector response_fraction returns one coefficient row per response_fraction", {
    x <- c(rep(0, 5), seq(0, 100, length.out = 20))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x),
        nirs_channels = "x",
        time_channel = "t"
    )
    result <- analyse_kinetics(
        df,
        method = "response_time",
        start_time = 5,
        response_fraction = c(0.5, 0.632),
        verbose = FALSE
    )
    coefs <- result$coefficients

    expect_equal(coefs$response_fraction, c(0.5, 0.632))
    expect_equal(coefs$fitted, coefs$A + (coefs$B - coefs$A) * coefs$response_fraction)
    ## higher response_fraction reached later
    expect_lt(coefs$response_time[1], coefs$response_time[2])
    ## vector arg collapsed into the single channel_args row
    expect_equal(result$channel_args$response_fraction, "0.5, 0.632")
})

test_that("analyse_response_time fitted_data contains baseline, response, extreme", {
    x <- c(rep(0, 5), seq(0, 20, length.out = 15), rep(20, 5))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x),
        nirs_channels = "x",
        time_channel = "t"
    )
    result <- analyse_response_time(df, start_time = 5, verbose = FALSE)

    fd <- attr(result, "fitted_data")$x

    ## fitted_data has window_idx and fitted columns
    expect_named(fd, c("window_idx", "fitted"))
    ## baseline A values are constant, response and extreme differ
    n_baseline <- length(which(t <= 5))
    baseline_fitted <- fd$fitted[seq_len(n_baseline)]
    expect_all_true(baseline_fitted == result$A)
})

test_that("analyse_response_time diagnostics structure", {
    x <- c(rep(0, 5), seq(0, 20, length.out = 15), rep(20, 5))
    t <- seq_along(x)

    df <- create_mnirs_data(
        data.frame(t, x),
        nirs_channels = "x",
        time_channel = "t"
    )
    result <- analyse_response_time(df, start_time = 5, verbose = FALSE)

    diag <- attr(result, "diagnostics")

    expect_s3_class(diag, "data.frame")
    expect_equal(nrow(diag), 1L)
    expect_equal(diag$nirs_channels, "x")
    expect_all_true(
        is.na(unlist(
            diag[!names(diag) %in% c(
                "interval", "nirs_channels", "n_obs", "n_params"
            )]
        ))
    )
})


## integration ====================================================
test_that("analyse_response_time works visually on Moxy", {
    skip_if_not_installed("ggplot2")
    skip("visual check")
    
    data <- read_mnirs(
        example_mnirs("moxy_ramp"),
        nirs_channels = c(
            smo2_left = "SmO2 Live",
            smo2_right = "SmO2 Live(2)"
        ),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        extract_intervals(
            start = by_time(878),
            span = c(-15, 120),
            verbose = FALSE,
            zero_time = FALSE
        ) |> 
        _[[1L]]

    t_vec <- data$time
    response_list <- response_time(
        data$smo2_left, t_vec, 878, response_fraction = 0.5
    )
    # response_time(data$smo2_right, t_vec, response_fraction = 0.5)

    result <- analyse_response_time(
        data, 
        nirs_channels = smo2_right,
        start_time = 878,
        response_fraction = 0.5
    )

    # result
    # attributes(result)
    hrt_data <- attr(result, "fitted_data")$smo2_right
    

    library(ggplot2)
    plot(data) +
        annotate(
            "point",
            x = t_vec[
                unlist(response_list[c("response_idx", "extreme_idx")])
            ],
            y = unlist(response_list[c("response_value", "B")]),
            size = 3,
            shape = 21,
            stroke = 1
        ) +
        annotate(
            "line",
            x = t_vec[response_list$baseline_idx],
            y = response_list$A,
        ) +
        annotate(
            "point",
            x = t_vec[rev(rev(hrt_data$window_idx)[1:3])],
            y = unlist(result[c("A", "response_value", "B")]),
            size = 3,
            shape = 21,
            stroke = 1
        ) +
        annotate(
            "line",
            x = rev(t_vec[hrt_data$window_idx])[-(1:2)],
            y = result$A,
        )
})
