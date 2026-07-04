test_that("resample_mnirs upsamples correctly", {
    data <- data.frame(
        time = c(0, 1, 2),
        value = c(10, 20, 30)
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 1,
        resample_rate = 4,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(nrow(result), 9) # 0 to 2 at 0.25 intervals
    expect_equal(result$time, seq(0, 2, by = 0.25))
    expect_equal(result$value, seq(10, 30, 2.5))

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 1,
        resample_rate = 4,
        method = "none", ## return only existing samples
        verbose = FALSE
    )
    expect_equal(nrow(result), 9) # 0 to 2 at 0.25 intervals
    expect_equal(result$time, seq(0, 2, by = 0.25))
    expect_equal(result$value[1], 10)
    expect_equal(result$value[5], 20)
    expect_equal(result$value[9], 30)
    expect_true(all(is.na(result$value[c(2:4, 6:8)])))
})

test_that("resample_mnirs downsamples correctly", {
    data <- data.frame(
        time = seq(0, 2, by = 0.1),
        value = seq(10, 30, length.out = 21)
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 10,
        resample_rate = 1,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(nrow(result), 3) # 0, 1, 2
    expect_equal(result$time, c(0, 1, 2))
    expect_equal(result$value, seq(10, by = 10, len = 3))

    result <- resample_mnirs(
        data,
        "time",
        10,
        1,
        method = "none",
        verbose = FALSE
    )

    expect_equal(nrow(result), 3) # 0, 1, 2
    expect_equal(result$time, c(0, 1, 2))
    expect_equal(result$value[1], 10)
    expect_equal(result$value[2], 20)
    expect_equal(result$value[3], 30)
})

test_that("resample_mnirs handles resample_rate == sample_rate", {
    data <- data.frame(time = 1:3, value = c(10, 20, 30))
    expect_equal(
        resample_mnirs(data, "time", 1, 1, verbose = FALSE),
        data,
        ignore_attr = TRUE
    )
    expect_equal(
        resample_mnirs(data, "time", 1, verbose = FALSE),
        data,
        ignore_attr = TRUE
    )
    expect_equal(
        resample_mnirs(data, "time", 1, resample_rate = NULL, verbose = FALSE),
        data,
        ignore_attr = TRUE
    )

    ## silent sample_rate and resample_rate
    attr(data, "sample_rate")
    expect_equal(
        resample_mnirs(data, "time", verbose = FALSE),
        data,
        ignore_attr = TRUE
    )

    ## uneven sampling
    set.seed(13)
    data <- data.frame(time = 1:3 + rnorm(3, 0, 0.1), value = c(10, 20, 30))
    result <- resample_mnirs(data, "time", 1, 1, verbose = FALSE)
    expect_true(all.equal(result$value, data$value, tolerance = 1, scale = 1))

    result <- resample_mnirs(data, "time", 1, verbose = FALSE)
    expect_true(all.equal(result$value, data$value, tolerance = 1, scale = 1))

    result <- resample_mnirs(
        data,
        "time",
        1,
        1,
        method = "none",
        verbose = FALSE
    )
    expect_equal(result$value, data$value)

    result <- resample_mnirs(data, "time", 1, method = "none", verbose = FALSE)
    expect_equal(result$value, data$value)
})

test_that("resample_mnirs handles repeated samples", {
    data <- data.frame(
        time = c(1:2, 2, 3:9, 9, 10:17, 17, 18:21) / 10 + 0.1,
        value = seq(10, by = 1, len = 24)
    )

    result <- resample_mnirs(data, "time", 10, 10, "linear", verbose = FALSE)
    expect_equal(range(result$time), floor(range(data$time) * 10) / 10)
    expect_equal(result$value[2], mean(data$value[2:3]))
    expect_equal(result$value[9], mean(data$value[10:11]))
    expect_equal(result$value[17], mean(data$value[19:20]))

    result <- resample_mnirs(
        data, "time", 10, 10, "none", verbose = FALSE
    )
    expect_equal(range(result$time), floor(range(data$time) * 10) / 10)
    expect_equal(result$value[2], data$value[2])
    expect_equal(result$value[9], data$value[10])
    expect_equal(result$value[17], data$value[19])
})

test_that("resample_mnirs handles missing samples", {
    data <- data.frame(
        time = c(1, 3:8, 10:16, 18:21) / 10 + 0.1,
        value = seq(10, by = 1, len = 18)
    )

    result <- resample_mnirs(
        data,
        "time",
        10,
        10,
        method = "linear",
        verbose = FALSE
    )
    expect_equal(range(result$time), floor(range(data$time) * 10) / 10)
    expect_equal(result$value[2], mean(data$value[1:2]))
    expect_equal(result$value[9], mean(data$value[7:8]))
    expect_equal(result$value[17], mean(data$value[14:15]))

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 10,
        method = "none",
        verbose = FALSE
    )
    expect_equal(range(result$time), floor(range(data$time) * 10) / 10)
    expect_true(all(is.na(result$value[c(2, 9, 17)])))
})

test_that("resample_mnirs handles multiple numeric columns", {
    data <- data.frame(
        time = c(0, 1, 2),
        value1 = c(10, 20, 30),
        value2 = c(5, 15, 25)
    )

    result <- resample_mnirs(data, "time", 1, 2, verbose = FALSE)

    expect_equal(ncol(result), 3)
    expect_true(all(c("time", "value1", "value2") %in% names(result)))
    expect_equal(result$value1[1], 10)
    expect_equal(result$value1[3], 20)
    expect_equal(result$value1[5], 30)
    expect_equal(result$value2[1], 5)
    expect_equal(result$value2[3], 15)
    expect_equal(result$value2[5], 25)
})

test_that("resample_mnirs handles categorical columns", {
    data <- data.frame(
        time = c(0, 1, 2),
        value = c(10, 20, 30),
        category = c("A", "B", "C")
    )

    result <- resample_mnirs(data, "time", 1, 4, "linear", verbose = FALSE)
    expect_true("category" %in% names(result))
    expect_equal(result$category[1], "A") # forward fill
    expect_equal(result$category[5], "B") # at t=1
    expect_equal(result$category[9], "C") # at t=2

    result <- resample_mnirs(data, "time", 1, 0.5, "linear", verbose = FALSE)
    expect_true("category" %in% names(result))
    expect_equal(result$category[1], "A") # forward fill
    expect_equal(result$category[2], "C") # at t=2

    result <- resample_mnirs(data, "time", 1, 4, "none", verbose = FALSE)
    expect_true("category" %in% names(result))
    expect_equal(result$category[1], "A") # forward fill
    expect_equal(result$category[5], "B") # at t=1
    expect_equal(result$category[9], "C") # at t=2
    expect_true(all(is.na(result$category[c(2:4, 6:8)])))

    data <- data.frame(
        time = seq(0, 4, 0.5),
        category = NA
    )
    data$category[c(1, 2, 5, 8, 9)] <- c("A", "B", "C", "D", "E")

    result <- resample_mnirs(
        data, "time", 2, 1, "linear", verbose = FALSE
    )
    expect_equal(result$category, c("A", NA, "C", "D", "E"))
})

test_that("non-numeric columns: down-sampling all NA in interval", {
    data <- tibble(
        time = seq(0, 4, by = 0.5),
        value = 1:9,
        category = c("A", NA, NA, NA, "C", NA, NA, NA, "E")
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 2,
        resample_rate = 1,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(result$category, c("A", NA, "C", NA, "E"))

    data <- tibble(
        time = seq(0.45, 4.95, 0.5),
        event = NA
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 2,
        resample_rate = 1,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(nrow(result), 6L)
})

test_that("non-numeric columns: method = 'NA' tolerance matching", {
    data <- tibble(
        time = c(0, 0.5, 1.0, 1.5, 2.0),
        value = 1:5,
        category = c("A", "B", "C", "D", "E")
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 2,
        resample_rate = 1,
        method = "none",
        verbose = FALSE
    )

    expect_equal(
        result$category,
        c("A", "C", "E")
    )
})

test_that("non-numeric columns: character and factor types", {
    data <- tibble(
        time = seq(0, 2, by = 0.5),
        value = 1:5,
        char_col = c("A", "B", "C", "D", "E"),
        factor_col = factor(c("low", "high", "low", "high", "low"))
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 2,
        resample_rate = 1,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(result$char_col, c("A", "C", "E"))
    expect_equal(
        as.character(result$factor_col),
        c("low", "low", "low")
    )
    expect_s3_class(result$factor_col, "factor")
})

test_that("resample_mnirs forward-fills when rates exceed actual rate", {
    data <- data.frame(
        time = c(0, 1, 2), ## actual sample_rate == 1
        value = c(10, 20, 30),
        category = c("A", "B", "C")
    )

    result <- resample_mnirs(
        data,
        time_channel = "time",
        sample_rate = 5,
        resample_rate = 2,
        method = "linear",
        verbose = FALSE
    )

    expect_equal(nrow(result), length(seq(0, 2, by = 1 / 2)))
    expect_false(anyNA(result$category))
    expect_equal(result$category, c("A", "A", "B", "B", "C"))
})


test_that("resample_mnirs handles edge cases", {
    ## Intentional that this should return error, even though technically
    ## resampling from one sample to one sample should be ok
    ## It just doesn't make sense for most applications
    # Single row
    data <- data.frame(time = 1, value = 10)
    expect_error(resample_mnirs(data, "time", 1), "valid.*numeric")
})


test_that("resample_mnirs works on Moxy", {
    file_path <- example_mnirs("moxy_ramp.xlsx")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(smo2 = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )[1:15, ]

    ## works with metadata
    expect_message(
        result <- resample_mnirs(df, resample_rate = 1),
        "Output is resampled at .*1.*Hz"
    )
    expect_equal(result$time, 0:8)
    expect_s3_class(result, "mnirs")

    ## time-weighted average
    df2 <- df |>
        dplyr::mutate(
            diff = c(diff(time), diff(time)[length(diff(time))]),
            time = round(time * 1) / 1,
        ) |>
        dplyr::summarise(
            .by = time,
            smo2 = stats::weighted.mean(smo2, diff)
        )

    ## expect close enough to time-weighted average
    expect_true(all.equal(
        result, df2, tolerance = 1, scale = 1, check.attributes = FALSE
    ))

    ## should overwrite metadata
    df3 <- resample_mnirs(
        df,
        sample_rate = 2,
        resample_rate = 1.1,
        verbose = FALSE
    )
    expect_equal(attributes(df3)$sample_rate, 1.1)

    ## method = "none"
    result <- resample_mnirs(
        df,
        sample_rate = 2,
        resample_rate = 4,
        method = "none",
        verbose = FALSE
    )

    expect_true(all(result$smo2[!is.na(result$smo2)] == df$smo2))
    ## differences between original and interpolated time samples
    ## should be less than half of the resample_rate = 4 Hz = 0.125 sec
    expect_true(
        all(diff(sort(c(result$time[is.na(result$smo2)], df$time))) > 0.125)
    )
})


test_that("resample_mnirs updates metadata correctly", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        sample_rate = 2,
        verbose = FALSE
    ) |>
        resample_mnirs(
            resample_rate = 1,
            verbose = FALSE
        )
    expect_equal(attr(data, "nirs_channels"), c("smo2_left", "smo2_right"))
    expect_equal(attr(data, "time_channel"), "time")
    expect_equal(attr(data, "sample_rate"), 1)
})


test_that("resample_mnirs works visually on moxy data", {
    skip_if_not_installed("ggplot2")
    skip("visual check of ggplots")
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2 = "SmO2 Live"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )

    (p <- plot(data) +
        ggplot2::scale_colour_manual(
            breaks = c("smo2", "resample"),
            values = palette_mnirs(2)
        ) +
        ggplot2::ylim(0, 100) +
        ggplot2::xlim(0, 100))

    p +
        ggplot2::geom_line(
            data = resample_mnirs(data, method = "none"),
            ggplot2::aes(y = smo2, colour = "resample")
        )

    p +
        ggplot2::geom_line(
            data = resample_mnirs(data, method = "linear"),
            ggplot2::aes(y = smo2, colour = "resample")
        )

    p +
        ggplot2::geom_line(
            data = resample_mnirs(data, method = "locf"),
            ggplot2::aes(y = smo2, colour = "resample")
        )
})


## multi-interval input ================================================
test_that("resample_mnirs processes a named list of data frames", {
    make_df <- \(vals) {
        create_mnirs_data(
            data.frame(time = c(1, 2, 4, 5), ch1 = vals),
            nirs_channels = "ch1",
            time_channel = "time",
            sample_rate = 1
        )
    }
    data_list <- list(
        a = make_df(c(50, 51, 53, 54)),
        b = make_df(c(60, 61, 63, 64))
    )

    result <- resample_mnirs(data_list, method = "linear", verbose = FALSE)

    expect_type(result, "list")
    expect_named(result, c("a", "b"))
    expect_s3_class(result$a, "mnirs")
    ## gap at t = 3 regularised and filled by interpolation
    expect_equal(result$a$time, 1:5)
    expect_equal(result$a$ch1, c(50, 51, 52, 53, 54))
    expect_equal(result$b$ch1, c(60, 61, 62, 63, 64))
})

test_that("resample_mnirs processes grouped data frames", {
    skip_if_not_installed("dplyr")

    df <- create_mnirs_data(
        data.frame(
            time = rep(c(1, 2, 4, 5), 2),
            ch1 = c(50, 51, 53, 54, 60, 61, 63, 64),
            group = rep(c("A", "B"), each = 4)
        ),
        nirs_channels = "ch1",
        time_channel = "time",
        sample_rate = 1
    )
    grouped_df <- dplyr::group_by(df, group)

    result <- resample_mnirs(grouped_df, method = "linear", verbose = FALSE)

    expect_named(result, c("A", "B"))
    expect_s3_class(result$A, "mnirs")
    expect_equal(result$A$ch1, c(50, 51, 52, 53, 54))
    expect_equal(result$B$ch1, c(60, 61, 62, 63, 64))
})
