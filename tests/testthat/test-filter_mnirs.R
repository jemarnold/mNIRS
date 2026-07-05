## filter_moving_average() ===========================================
test_that("filter_moving_average() returns expected smoothed values", {
    x <- c(1, 2, 3, 4, 5)

    result <- filter_moving_average(x, width = 3, verbose = FALSE, na.rm = TRUE)
    expect_equal(result, c(NA, 2, 3, 4, NA))
    result <- filter_moving_average(x, width = 3, verbose = FALSE, partial = TRUE)
    expect_equal(result, c(1.5, 2, 3, 4, 4.5))

    # Width-based with single width (floor(1/2) = 0, so just x itself)
    result <- filter_moving_average(x, width = 1, verbose = FALSE)
    expect_equal(result, x)
    
    ## filter_ma alias for filter_moving_average
    result <- filter_ma(
        x,
        width = 3,
        verbose = FALSE,
        partial = FALSE
    )
    expect_equal(result, c(NA, 2, 3, 4, NA))
    result <- filter_ma(
        x,
        width = 3,
        verbose = FALSE,
        partial = TRUE
    )
    expect_equal(result, c(1.5, 2, 3, 4, 4.5))
})

test_that("filter_moving_average() handles custom time vectors", {
    x <- c(10, 20, 30, 40, 50, 60)
    t <- c(0, 1, 3, 6, 10, 11)

    # Span-based with time gaps
    result <- filter_moving_average(x, t, span = 2, verbose = FALSE)
    expect_equal(result, c(15, 15, 30, 40, 55, 55))
})

test_that("filter_moving_average() handles NA values correctly", {
    x <- c(1, NA, 3, 4, 5)

    result <- filter_moving_average(
        x,
        width = 3,
        partial = FALSE,
        na.rm = FALSE,
        verbose = FALSE
    )
    expect_equal(result, c(rep(NA, 3), 4, NA))

    result <- filter_moving_average(
        x,
        width = 3,
        partial = TRUE,
        na.rm = TRUE,
    )
    expect_equal(result, c(1, 2, 3.5, 4, 4.5))

    result <- filter_moving_average(
        x,
        width = 3,
        partial = FALSE, ## min 3 valid obs at edges
        na.rm = TRUE, ## mean excludes NA
    )
    expect_equal(result, c(NA, 2, 3.5, 4, NA))

    result <- filter_moving_average(
        x,
        width = 3,
        partial = TRUE,
        na.rm = FALSE,
        verbose = FALSE
    )
    expect_equal(result, c(rep(NA, 3), 4, 4.5))

    x <- c(1, NA, NA, NA, NA, 6, 7, 8)
    result <- filter_moving_average(x, width = 3, partial = TRUE, na.rm = TRUE)
    expect_equal(result, c(1, 1, NA, NA, 6, 6.5, 7, 7.5))
})

test_that("filter_moving_average() validates x as numeric", {
    expect_error(
        filter_moving_average(c("a", "b", "c"), width = 3, verbose = FALSE),
        "x.*numeric"
    )
})

test_that("filter_moving_average() validates t as numeric", {
    expect_error(
        filter_moving_average(
            1:5,
            t = letters[1:5],
            width = 3,
            verbose = FALSE
        ),
        "t.*numeric"
    )
})

test_that("filter_moving_average() validates x and t have same length", {
    expect_error(
        filter_moving_average(1:5, t = 1:3, width = 3, verbose = FALSE),
        "equal length"
    )
})

test_that("filter_moving_average() validates width & span constraints", {
    x <- 1:10

    ## width or span required
    expect_error(
        filter_moving_average(x, verbose = FALSE),
        "width.*or.*span.*must be defined"
    )

    # Width must be positive integer
    expect_error(
        filter_moving_average(x, width = -1, verbose = FALSE),
        "width.*integer"
    )

    ## width = 0
    expect_error(
        filter_moving_average(x, width = 0, verbose = FALSE),
        "width.*integer"
    )
    ## width > length(x)
    expect_true(all(
        filter_moving_average(x, width = 21, partial = TRUE) == mean(x)
    ))

    # Span must be positive
    expect_error(
        filter_moving_average(x, span = -1, verbose = FALSE),
        "span.*numeric"
    )

    ## span = 0
    expect_equal(filter_moving_average(x, span = 0, verbose = FALSE), x)
    ## span = 0.5 takes nearest 1 sample
    expect_equal(filter_moving_average(x, span = 0.5, verbose = FALSE), x)
    ## span > length(x)
    expect_true(all(
        filter_moving_average(x, span = 21, partial = TRUE) == mean(x)
    ))
})

test_that("filter_moving_average() warns when both width and span provided", {
    ## should produce warning
    expect_message(
        filter_moving_average(1:5, width = 3, span = 2),
        "width.*span"
    )

    old_verbose <- getOption("mnirs.verbose")
    on.exit(options(mnirs.verbose = old_verbose), add = TRUE)
    options(mnirs.verbose = FALSE)

    ## should not produce warning
    expect_silent(
        result <- filter_moving_average(1:5, width = 3, span = 2)
    )
})

test_that("filter_moving_average() handles edge cases", {
    # Single value
    expect_equal(
        filter_moving_average(5, width = 1, verbose = FALSE),
        5,
        ignore_attr = TRUE
    )

    # Two values
    expect_equal(
        filter_moving_average(c(1, 2), width = 1, verbose = FALSE),
        c(1, 2)
    )

    # All NA
    expect_error(
        filter_moving_average(c(NA, NA, NA), width = 3, verbose = FALSE),
        "x.*numeric"
    )
})

test_that("filter_moving_average with insufficient width returns NA with warning", {
    t <- 1:10
    x <- c(1, NA, 3, NA, 5, NA, 7, NA, 9, NA)

    expect_all_equal(filter_moving_average(x, t, width = 2), NA_real_) |>
        expect_warning("Set.*na.rm = TRUE")

    x <- c(NA, 1, NA)
    expect_error(
        filter_moving_average(x, width = 3, na.rm = TRUE),
        "Insufficient.*samples"
    )
})

## filter_butterworth() =========================================
test_that("filter_butterworth validates inputs correctly", {
    expect_error(
        filter_butterworth(x = "not_numeric", order = 1, W = 0.1),
        "x.*numeric"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = 0, W = 0.1),
        "n.*integer"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = -1, W = 0.1),
        "n.*integer"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = 1.5, W = 0.1),
        "n.*integer"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = 1, W = 1.5),
        "W.*numeric"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = 1, W = 0),
        "W.*numeric"
    )

    expect_error(
        filter_butterworth(x = 1:10, order = 1, W = c(0.1, 0.3), type = "low"),
        "W.*numeric"
    )
})

test_that("filter_butterworth returns correct output structure", {
    x <- rnorm(100)
    result <- filter_butterworth(x, order = 2, W = 0.1)

    expect_type(result, "double")
    expect_length(result, length(x))
    expect_false(anyNA(result))
})

test_that("filter_butterworth handles different edge options", {
    set.seed(42)
    x <- sin(2 * pi * 1:50 / 10) + rnorm(50, 0, 0.5)

    result_none <- filter_butterworth(x, order = 2, W = 0.1, edges = "none")
    result_rev <- filter_butterworth(x, order = 2, W = 0.1, edges = "rev")
    result_rep1 <- filter_butterworth(x, order = 2, W = 0.1, edges = "rep1")

    expect_length(result_none, length(x))
    expect_length(result_rev, length(x))
    expect_length(result_rep1, length(x))

    # Results should differ at edges
    expect_false(identical(result_none, result_rev))
    expect_false(identical(result_none, result_rep1))
    expect_false(identical(result_rev, result_rep1))
})

test_that("filter_butterworth handles NA values", {
    x_with_na <- c(1:5, NA, 7:10)

    # NA propagates through filtering
    expect_error(
        filter_butterworth(x_with_na, order = 1, W = 0.1, na.rm = FALSE),
        "x.*NA"
    )

    result <- filter_butterworth(x_with_na, order = 1, W = 0.1, na.rm = TRUE)
    expect_true(is.na(result[6]))

    # All NA input
    expect_error(
        filter_butterworth(rep(NA_real_, 5), order = 1, W = 0.1, na.rm = TRUE),
        "x.*valid.*numeric"
    )

    # Leading/trailing NAs
    x_edge_na <- c(NA, NA, 3:8, NA, NA)
    result_edge <- filter_butterworth(x_edge_na, order = 1, W = 0.1, na.rm = TRUE)
    expect_true(all(is.na(result_edge[c(1:2, 9:10)])))
    expect_length(result_edge, length(x_edge_na))
})

test_that("filter_butterworth handles different filter types", {
    set.seed(123)
    x <- rnorm(100)

    low <- filter_butterworth(x, order = 2, W = 0.1, type = "low")
    high <- filter_butterworth(x, order = 2, W = 0.9, type = "high")
    stop <- filter_butterworth(x, order = 2, W = c(0.3, 0.7), type = "stop")
    pass <- filter_butterworth(x, order = 2, W = c(0.3, 0.7), type = "pass")

    expect_length(low, length(x))
    expect_length(high, length(x))
    expect_length(stop, length(x))
    expect_length(pass, length(x))

    ## wrong W elements
    expect_error(
        filter_butterworth(x, order = 2, W = c(0.3, 0.7), type = "low"),
        "W.*1-element.*numeric"
    )
    expect_error(
        filter_butterworth(x, order = 2, W = 0.9, type = "pass"),
        "W.*2-element.*numeric"
    )
})

test_that("filter_butterworth handles edge cases", {
    # Very short vector
    x_short <- c(1, 2, 3)
    result <- filter_butterworth(x_short, order = 1, W = 0.5)
    expect_length(result, 3)

    # Single value
    x_single <- 5
    result_single <- filter_butterworth(x_single, order = 1, W = 0.5)
    expect_length(result_single, 1)

    # Constant signal
    x_const <- rep(10, 50)
    result_const <- filter_butterworth(x_const, order = 1, W = 0.5)
    expect_true(all(abs(result_const - 10) < 1e-10))
})

test_that("filter_butterworth smooths noisy signal", {
    set.seed(999)
    sin_wave <- sin(2 * pi * 1:100 / 20)
    x <- seq_along(sin_wave)
    noisy <- sin_wave + rnorm(100, 0, 0.3)
    filtered <- filter_butterworth(noisy, order = 2, W = 0.2)

    # Filtered signal should be smoother (lower variance)
    expect_lt(var(diff(filtered)), var(diff(noisy)))
})

test_that("filter_butterworth works visually", {
    skip("visual check for filter_butterworth")
    set.seed(999)
    sin_wave <- sin(2 * pi * 1:100 / 20)
    x <- seq_along(sin_wave)
    noisy <- sin_wave + rnorm(100, 0, 0.3)
    filtered <- filter_butterworth(noisy, order = 2, W = 0.2)

    ggplot2::ggplot(tibble()) +
        ggplot2::aes(x = x) +
        ggplot2::geom_line(ggplot2::aes(y = sin_wave, colour = "sin")) +
        ggplot2::geom_line(ggplot2::aes(y = noisy, colour = "noisy")) +
        ggplot2::geom_line(ggplot2::aes(y = filtered, colour = "filtered"))
})


## moxy_data for filter_mnirs() ==============================================
moxy_data <- read_mnirs(
    file_path = example_mnirs("moxy_ramp.xlsx"),
    nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
    time_channel = c(time = "hh:mm:ss"),
    verbose = FALSE
) |>
    resample_mnirs(method = "linear", verbose = FALSE)

## Input validation tests =======================================
test_that("filter_mnirs validates input data", {
    expect_error(
        filter_mnirs(data.frame(x = 1:10)),
        "`data`.*data.*frame"
    )
})

test_that("filter_mnirs method aliases work", {
    result_ma <- filter_mnirs(
        moxy_data, method = "moving_average", width = 5, verbose = FALSE
    )
    result_spl <- filter_mnirs(
        moxy_data, method = "smooth_spline", verbose = FALSE
    )

    expect_equal(
        filter_mnirs(
            moxy_data, method = "moving average", width = 5, verbose = FALSE
        ),
        result_ma
    )
    expect_equal(
        filter_mnirs(
            moxy_data, method = "Moving Average", width = 5, verbose = FALSE
        ),
        result_ma
    )
    expect_equal(
        filter_mnirs(
            moxy_data, method = "MA", width = 5, verbose = FALSE
        ),
        result_ma
    )
    expect_equal(
        filter_mnirs(moxy_data, verbose = FALSE),
        result_spl
    )
    expect_equal(
        filter_mnirs(moxy_data, method = "spline", verbose = FALSE),
        result_spl
    )
    expect_equal(
        filter_mnirs(moxy_data, method = "smooth spline", verbose = FALSE),
        result_spl
    )
    expect_equal(
        filter_mnirs(moxy_data, method = "Smooth Spline", verbose = FALSE),
        result_spl
    )
})

test_that("filter_mnirs validates method argument", {
    expect_error(
        filter_mnirs(moxy_data, method = "invalid", verbose = FALSE),
        "should be one of"
    )
})

## Smooth spline tests ==============================================
test_that("smooth_spline filters data correctly", {
    result <- filter_mnirs(moxy_data, method = "smooth_spline", verbose = FALSE)

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))
    expect_true(all(c("time", "smo2_left", "smo2_right") %in% names(result)))

    # Filtered data should be smoother (lower variance)
    expect_lt(var(result$smo2_left), var(moxy_data$smo2_left))
    expect_lt(var(result$smo2_right), var(moxy_data$smo2_right))
})

test_that("smooth_spline respects spar parameter", {
    result_low <- filter_mnirs(
        moxy_data,
        method = "smooth_spline",
        spar = 0.1,
        verbose = FALSE
    )
    result_high <- filter_mnirs(
        moxy_data,
        method = "smooth_spline",
        spar = 0.9,
        verbose = FALSE
    )

    # Higher spar = more smoothing
    expect_lt(var(result_high$smo2_left), var(result_low$smo2_left))
    expect_lt(var(result_high$smo2_right), var(result_low$smo2_right))
})

test_that("smooth_spline handles NAs", {
    moxy_data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        resample_mnirs(method = "none", verbose = FALSE)

    result <- filter_mnirs(
        moxy_data,
        method = "smooth_spline",
        na.rm = TRUE,
        verbose = FALSE
    )

    expect_true(anyNA(result$smo2_left))
    expect_true(anyNA(result$smo2_right))
    expect_equal(
        which(is.na(result$smo2_right)),
        which(is.na(moxy_data$smo2_right))
    )

    expect_error(
        filter_mnirs(
            moxy_data, method = "smooth_spline", na.rm = FALSE,
            verbose = FALSE
        ),
        "NA"
    )
})

test_that("smooth_spline errors with irregular samples", {
    moxy_data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )

    expect_error(
        filter_mnirs(
            moxy_data, method = "smooth_spline", verbose = FALSE
        ),
        "irregular samples"
    )
})

test_that("smooth_spline validates spar parameter", {
    expect_error(
        filter_mnirs(
            moxy_data, method = "smooth_spline", spar = -1,
            verbose = FALSE
        ),
        "spar.*numeric"
    )

    expect_error(
        filter_mnirs(
            moxy_data, method = "smooth_spline", spar = c(0.1, 0.5),
            verbose = FALSE
        ),
        "spar.*numeric"
    )
})

## Butterworth tests ==============================================
test_that("butterworth low-pass filter works", {
    result <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        W = 0.1,
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))

    # Low-pass should reduce high-frequency variance
    expect_lt(var(diff(result$smo2_left)), var(diff(moxy_data$smo2_left)))
    expect_lt(var(diff(result$smo2_right)), var(diff(moxy_data$smo2_right)))
})

test_that("butterworth accepts fc instead of W", {
    sr <- attr(moxy_data, "sample_rate")

    result_W <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        W = 0.1,
        verbose = FALSE
    )

    result_fc <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        fc = 0.1 * sr * 0.5,
        verbose = FALSE
    )

    expect_equal(result_W$smo2_left, result_fc$smo2_left)
})

test_that("butterworth prefers W when both W and fc specified", {
    expect_message(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            W = 0.1,
            fc = 5
        ),
        ".*W.*fc"
    )
})

test_that("butterworth errors without W or fc", {
    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            verbose = FALSE
        ),
        "W.*fc.*must be defined"
    )
})

test_that("butterworth errors for invalid fc", {
    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            fc = 100,
            verbose = FALSE
        ),
        "must be between.*0.*half"
    )
})

test_that("butterworth validates filter type", {
    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "invalid",
            W = 0.1,
            verbose = FALSE
        ),
        "must be one of"
    )
})

test_that("butterworth handles different filter types", {
    result_low <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        W = 0.1,
        verbose = FALSE
    )

    result_high <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "high",
        order = 2,
        W = 0.1,
        verbose = FALSE
    )

    result_stop <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "stop",
        order = 2,
        W = c(0.05, 0.15),
        verbose = FALSE
    )

    result_pass <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "pass",
        order = 2,
        W = c(0.05, 0.15),
        verbose = FALSE
    )

    expect_s3_class(result_low, "mnirs")
    expect_s3_class(result_high, "mnirs")
    expect_s3_class(result_stop, "mnirs")
    expect_s3_class(result_pass, "mnirs")
})

test_that("butterworth validates W length for stop/pass filters", {
    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "stop",
            order = 2,
            W = 0.5,
            verbose = FALSE
        ),
        "W.*numeric"
    )

    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "pass",
            order = 2,
            fc = 5,
            verbose = FALSE
        ),
        "fc.*numeric"
    )
})

test_that("butterworth handles NAs with na.rm = TRUE", {
    moxy_data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        resample_mnirs(method = "none", verbose = FALSE)

    result <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        W = 0.1,
        na.rm = TRUE,
        verbose = FALSE
    )

    expect_true(anyNA(result$smo2_left))
    expect_true(anyNA(result$smo2_right))
    expect_equal(
        which(is.na(result$smo2_right)),
        which(is.na(moxy_data$smo2_right))
    )

    expect_error(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            W = 0.1,
            na.rm = FALSE,
            verbose = FALSE
        ),
        "NA"
    )
})

## Moving average tests ==============================================
test_that("moving_average with width works", {
    result <- filter_mnirs(
        moxy_data,
        method = "moving_average",
        width = 10,
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))

    # Should smooth the data
    expect_lt(
        var(diff(result$smo2_left), na.rm = TRUE),
        var(diff(moxy_data$smo2_left), na.rm = TRUE)
    )
    expect_lt(
        var(diff(result$smo2_right), na.rm = TRUE),
        var(diff(moxy_data$smo2_right), na.rm = TRUE)
    )
})

test_that("moving_average with span works", {
    result <- filter_mnirs(
        moxy_data,
        method = "moving_average",
        span = 1,
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))
})

test_that("moving_average with larger width = more smoothing", {
    result_narrow <- filter_mnirs(
        moxy_data,
        method = "moving_average",
        width = 5,
        verbose = FALSE
    )

    result_wide <- filter_mnirs(
        moxy_data,
        method = "moving_average",
        width = 20,
        verbose = FALSE
    )

    expect_lt(
        var(result_wide$smo2_left, na.rm = TRUE),
        var(result_narrow$smo2_left, na.rm = TRUE)
    )
})

## Channel selection tests ==============================================
test_that("filter_mnirs respects nirs_channels argument", {
    result <- filter_mnirs(
        moxy_data,
        nirs_channels = "smo2_right",
        method = "smooth_spline",
        spar = 0.5,
        verbose = FALSE
    )

    # channel1 should be filtered
    expect_false(identical(result$smo2_right, moxy_data$smo2_right))
    expect_lt(var(diff(result$smo2_right)), var(diff(moxy_data$smo2_right)))

    expect_true(identical(result$smo2_left, moxy_data$smo2_left))
    expect_equal(var(diff(result$smo2_left)), var(diff(moxy_data$smo2_left)))

    expect_equal(attr(result, "nirs_channels"), "smo2_right")

    expect_false(isTRUE(all.equal(
        attr(result, "nirs_channels"),
        attr(moxy_data, "nirs_channels")
    )))
})

## Metadata preservation tests ==============================================
test_that("filter_mnirs preserves and updates metadata", {
    original_meta <- attributes(moxy_data)
    result <- filter_mnirs(moxy_data, method = "smooth_spline", verbose = FALSE)
    result_meta <- attributes(result)

    expect_equal(result_meta$nirs_channels, original_meta$nirs_channels)
    expect_equal(result_meta$time_channel, original_meta$time_channel)
})

test_that("butterworth updates sample_rate in metadata", {
    result <- filter_mnirs(
        moxy_data,
        method = "butterworth",
        type = "low",
        order = 2,
        W = 0.1,
        sample_rate = 10,
        verbose = FALSE
    )

    expect_equal(attr(result, "sample_rate"), 10)
    expect_false(attr(result, "sample_rate") == attr(moxy_data, "sample_rate"))
})

## verbose output tests ==============================================
test_that("verbose output works", {
    expect_message(
        filter_mnirs(moxy_data, method = "smooth_spline", verbose = TRUE),
        "spar ="
    ) |>
        expect_message("spar =")

    expect_silent(
        filter_mnirs(moxy_data, method = "smooth_spline", verbose = FALSE)
    )

    expect_message(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            W = 0.1,
            fc = 0.1,
            verbose = TRUE
        ),
        "W.*overrides.*fc"
    )

    expect_silent(
        filter_mnirs(
            moxy_data,
            method = "butterworth",
            type = "low",
            order = 2,
            W = 0.1,
            fc = 0.1,
            verbose = FALSE
        )
    )
})

## Integration tests ==============================================
test_that("filter_mnirs works with pipe", {
    result <- moxy_data |>
        filter_mnirs(method = "smooth_spline", verbose = FALSE)

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))
    expect_lt(var(diff(result$smo2_left)), var(diff(moxy_data$smo2_left)))
})

test_that("multiple filtering operations can be chained", {
    result <- moxy_data |>
        filter_mnirs(
            method = "butterworth",
            type = "high",
            order = 1,
            W = 0.01,
            verbose = FALSE
        ) |>
        filter_mnirs(method = "smooth_spline", spar = 0.3, verbose = FALSE)

    expect_s3_class(result, "mnirs")
    expect_equal(nrow(result), nrow(moxy_data))
})

test_that("filter_mnirs works visually on Moxy data", {
    skip("visual check for filter_mnirs() on Moxy data")

    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2 = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        resample_mnirs(method = "linear", verbose = FALSE) |>
        dplyr::mutate(
            dplyr::across(
                smo2,
                \(.x) filter_butterworth(.x, 2, 0.05),
                .names = "{.col}_filtfilt"
            )
        )

    data_filt <- filter_mnirs(
        data,
        nirs_channels = "smo2",
        time_channel = NULL,
        sample_rate = NULL,
        method = "moving_average", #c("smooth_spline", "butterworth", "moving_average"),
        # spar = 0.5,
        # type = c("low"),
        # order = 2,
        # W = 0.05,
        # fc = NULL,
        width = NULL,
        span = 30,
        # na.rm = FALSE,
        verbose = TRUE
    )

    ## visual check
    plot(data) +
        ggplot2::scale_colour_manual(
            breaks = c("smo2", "filt_vec", "filt_mnirs"),
            values = palette_mnirs()
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = smo2_filtfilt, colour = "filt_vec"),
            linewidth = 1
        ) +
        ggplot2::geom_line(
            data = data_filt,
            ggplot2::aes(y = smo2, colour = "filt_mnirs")
        )
})


## multi-interval input ================================================
test_that("filter_mnirs processes a named list of data frames", {
    make_df <- \(val) {
        create_mnirs_data(
            data.frame(time = 1:5, ch1 = rep(val, 5)),
            nirs_channels = "ch1",
            time_channel = "time"
        )
    }
    data_list <- list(a = make_df(50), b = make_df(60))

    result <- filter_mnirs(
        data_list,
        method = "moving_average",
        width = 3,
        verbose = FALSE
    )

    expect_type(result, "list")
    expect_named(result, c("a", "b"))
    expect_s3_class(result$a, "mnirs")
    ## constant input returns constant means; partial edge windows are NA
    expect_equal(result$a$ch1, c(NA, 50, 50, 50, NA))
    expect_equal(result$b$ch1, c(NA, 60, 60, 60, NA))
})

test_that("filter_mnirs processes grouped data frames", {
    skip_if_not_installed("dplyr")

    df <- create_mnirs_data(
        data.frame(
            time = rep(1:5, 2),
            ch1 = rep(c(50, 60), each = 5),
            group = rep(c("A", "B"), each = 5)
        ),
        nirs_channels = "ch1",
        time_channel = "time"
    )
    grouped_df <- dplyr::group_by(df, group)

    result <- filter_mnirs(
        grouped_df,
        method = "moving_average",
        width = 3,
        verbose = FALSE
    )

    expect_named(result, c("A", "B"))
    expect_s3_class(result$A, "mnirs")
    expect_equal(result$A$ch1, c(NA, 50, 50, 50, NA))
    expect_equal(result$B$ch1, c(NA, 60, 60, 60, NA))
})
