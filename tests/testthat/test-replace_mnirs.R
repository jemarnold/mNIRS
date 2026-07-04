## preserve_na() -------------------------------------------------------
test_that("preserve_na() correctly stores NA positions and clean values", {
    x <- c(1, NA, 3, NA, 5)
    na_info <- preserve_na(x)

    expect_type(na_info, "list")
    expect_named(na_info, c("x_valid", "x_length", "na_idx"))
    expect_equal(na_info$x_valid, x[!is.na(x)])
    expect_equal(na_info$x_length, 5)
    expect_equal(na_info$na_idx, is.na(x))
})

test_that("preserve_na() handles vector with no NA", {
    x <- 1:5
    na_info <- preserve_na(x)

    expect_equal(na_info$x_valid, x)
    expect_equal(na_info$x_length, 5)
    expect_equal(na_info$na_idx, is.na(x))
})

test_that("preserve_na() handles vector with all NA", {
    x <- rep(NA_real_, 5)
    na_info <- preserve_na(x)

    expect_length(na_info$x_valid, 0)
    expect_equal(na_info$x_length, 5)
    expect_equal(na_info$na_idx, is.na(x))
})


## restore_na() --------------------------------------------------------
test_that("restore_na() correctly restores NA positions", {
    x <- c(1, NA, 3, NA, 5)
    na_info <- preserve_na(x)
    y <- na_info$x_valid * 2 ## y_valid = c(2, 6, 10)
    result <- restore_na(y, na_info)

    expect_length(result, 5)
    expect_equal(result, x * 2)
})

test_that("restore_na() handles no NA in original vector", {
    x <- 1:5
    na_info <- preserve_na(x)
    y <- na_info$x_valid * 2
    result <- restore_na(y, na_info)

    expect_equal(result, x * 2)
})

test_that("restore_na() handles all NA in original vector", {
    x <- rep(NA_real_, 5)
    na_info <- preserve_na(x)
    y <- na_info$x_valid # no clean values to process
    result <- restore_na(y, na_info)

    expect_length(result, 5)
    expect_true(all(is.na(result)))
})

test_that("preserve_na() and restore_na() work together correctly", {
    x <- c(1, NA, 3, NA, 5, 6, NA, 8)
    na_info <- preserve_na(x)

    # Simulate processing
    y_processed <- na_info$x_valid + 10
    result <- restore_na(y_processed, na_info)

    expect_length(result, length(x))
    expect_equal(result[!is.na(x)], x[!is.na(x)] + 10)
    expect_true(is.na(result[2]))
    expect_true(is.na(result[4]))
    expect_true(is.na(result[7]))
})


## replace_outliers() --------------------------------------------------
test_that("replace_outliers() returns unchanged vector with no outliers", {
    x <- 1:20
    result <- replace_outliers(x, width = 5)
    expect_equal(result, x)

    ## span
    result <- replace_outliers(x, span = 3)
    expect_equal(result, x)
})

test_that("replace_outliers() detects and replaces outliers with median", {
    x <- c(1:4, 100, 6:10) # 100 is clear outlier
    result <- replace_outliers(x, width = 5)

    expect_type(result, "double")
    expect_length(result, length(x))
    expect_lt(result[5], 100) # outlier replaced
    expect_equal(result[5], median(x[c(3:7)])) # outlier replaced
    expect_false(anyNA(result))

    ## span
    result <- replace_outliers(x, span = 3)

    expect_type(result, "double")
    expect_length(result, length(x))
    expect_lt(result[5], 100) # outlier replaced
    expect_equal(result[5], median(x[c(3:7)])) # outlier replaced
    expect_false(anyNA(result))
})

test_that("replace_outliers() detects and replaces outliers with NA", {
    x <- c(1:10, 100, 11:20)
    result <- replace_outliers(x, method = "none", width = 5)

    expect_length(result, length(x))
    expect_true(is.na(result[11])) # outlier replaced with NA

    ## span
    result <- replace_outliers(x, method = "none", span = 3)

    expect_length(result, length(x))
    expect_true(is.na(result[11])) # outlier replaced with NA
})

test_that("replace_outliers() handles NA values in input", {
    x <- c(1:9, NA, 100, 12:15, NA, 17:20)
    result <- replace_outliers(x, width = 5)

    median(c(8, 9, NA), na.rm = TRUE)

    expect_length(result, length(x))
    expect_true(result[11] != 100) # outlier replaced
    expect_true(all(is.na(result[c(10, 16)]))) # original NA preserved

    ## span
    result <- replace_outliers(x, span = 4)

    expect_length(result, length(x))
    expect_true(result[11] != 100) # outlier replaced
    expect_true(all(is.na(result[c(10, 16)]))) # original NA preserved
})

test_that("replace_outliers() respects outlier_cutoff threshold", {
    x <- c(1:10, 15, 11:20) # mild outlier

    strict <- replace_outliers(x, width = 5, outlier_cutoff = 1)
    lenient <- replace_outliers(x, width = 5, outlier_cutoff = 5)

    expect_true(strict[11] != 15) # detected with strict threshold
    expect_equal(lenient[11], 15) # not detected with lenient threshold

    ## span
    strict <- replace_outliers(x, span = 3, outlier_cutoff = 1)
    lenient <- replace_outliers(x, span = 3, outlier_cutoff = 5)

    expect_true(strict[11] != 15) # detected with strict threshold
    expect_equal(lenient[11], 15) # not detected with lenient threshold

    ## numeric outlier_cutoff
    expect_equal(replace_outliers(x, span = 3, outlier_cutoff = 2.5)[11], 11)

    ## Tukey's median filter.
    x <- c(1:10, 15, 11:20, 20.1) ## reduce min diff to 0.1 to avoid modification
    result <- replace_outliers(x, width = 5, outlier_cutoff = 0)
    medians <- vapply(seq_along(x), \(.idx) {
        median(x[pmax(1, .idx - 2):pmin(length(x), .idx + 2)])
    }, numeric(1))
    ## added tolerance to accomodate the robust variance threshold
    expect_true(all.equal(result, medians, tolerance = 0.11, scale = 1))
})

test_that("replace_outliers() validates inputs correctly", {
    x <- 1:10

    expect_error(replace_outliers(x = "text", width = 5), "x.*?numeric") # non-numeric x
    expect_error(
        replace_outliers(x, method = "none", width = -1),
        "width.*?integer"
    ) # negative width
    expect_error(
        replace_outliers(x, method = "none", width = 5, outlier_cutoff = -1),
        "outlier_cutoff.*?positive"
    ) # negative outlier_cutoff

    ## halfes all NA
    x <- rep(NA_real_, 10)
    expect_error(replace_outliers(x, width = 5), "x.*?numeric") ## x is all NA
    ## handles all same values
    expect_equal(replace_outliers(c(1, 1, 1), width = 1), rep(1, 3))

    ## span
    x <- 1:10
    expect_error(replace_outliers("text", span = 3), "x.*?numeric") # non-numeric x
    expect_error(
        replace_outliers(x, method = "none", span = -1),
        "span.*?numeric"
    ) # negative span

    ## haldes all NA
    x <- rep(NA_real_, 10)
    expect_error(replace_outliers(x, span = 3), "x.*?numeric") ## x is all NA
})

test_that("replace_outliers() errors when x and t have different lengths", {
    expect_error(
        replace_outliers(c(1, 2, 3), t = c(1, 2), outlier_cutoff = 3, span = 1),
        "equal length"
    )
})

## replace_invalid() ==================================================
test_that("replace_invalid() returns expected structure", {
    x <- c(1, 999, 3, 4, 999, 6)
    result <- replace_invalid(x, invalid_values = 999, width = 1)
    expect_length(result, length(x))
    expect_type(result, "double")
})

test_that("replace_invalid() replaces invalid values with NA when method = 'NA'", {
    x <- c(1, 999, 3, 4, 999, 6)
    result <- replace_invalid(x, invalid_values = 999, method = "none")
    expect_equal(result, c(1, NA, 3, 4, NA, 6))
})

test_that("replace_invalid() replaces invalid values with local median", {
    x <- c(1, 2, 999, 4, 999, 6)
    result <- replace_invalid(x, invalid_values = 999, width = 3)
    expect_equal(result, c(1, 2, 3, 4, 5, 6))
})

test_that("replace_invalid() handles multiple invalid values", {
    x <- c(1, 999, 3, -1, 5, 999)

    result <- replace_invalid(
        x,
        invalid_values = c(999, -1),
        width = 3
    )

    expect_equal(result, c(1, 2, 3, 4, 5, 5))
})

test_that("replace_invalid() uses custom time vector", {
    x <- c(1, 999, 3, 4)
    t <- c(0, 1, 10, 11)
    result <- replace_invalid(x, t, invalid_values = 999, span = 2)
    expect_equal(result, c(1, 1, 3, 4))
})

test_that("replace_invalid() handles non-integer span argument", {
    x <- c(1, 999, 3, 4, 999, 6)
    t <- c(0, 1, 2, 10, 11, 12)
    result <- replace_invalid(x, t, invalid_values = 999, span = 3)
    expect_equal(result, c(1, 2, 3, 4, 5, 6))
})

test_that("replace_invalid() errors when x and t have different lengths", {
    expect_error(
        replace_invalid(c(1, 2, 3), t = c(1, 2), invalid_values = 999),
        "equal length"
    )
})

test_that("replace_invalid() errors when args are not numeric", {
    expect_error(
        replace_invalid(c("a", "b"), invalid_values = 999),
        "`x`.*numeric"
    )

    expect_error(
        replace_invalid(c(1, 2), t = c("a", "b"), invalid_values = 999),
        "`t`.*numeric"
    )

    expect_error(
        replace_invalid(c(1, 2), invalid_values = "a"),
        "invalid.*numeric"
    )
})

test_that("replace_invalid() handles edge cases", {
    ## no invalid values present
    x <- c(1, 2, 3, 4)
    expect_equal(
        replace_invalid(x, invalid_values = 999, width = 1),
        x,
        ignore_attr = TRUE
    )
    expect_equal(
        replace_invalid(x, invalid_values = 999, width = 1, method = "median"),
        x,
        ignore_attr = TRUE
    )

    ## all invalid values
    x <- c(999, 999, 999)
    expect_equal(
        replace_invalid(x, invalid_values = 999, width = 1),
        rep(NA_real_, 3)
    )
    expect_equal(
        replace_invalid(x, invalid_values = 999, width = 1, method = "median"),
        rep(NA_real_, 3)
    )
})


test_that("replace_invalid handles range filtering", {
    x <- c(0, 5, 10, 20, 25, 50.5, 100, 150, 35, 40)

    # Test invalid_below only (inclusive)
    result <- replace_invalid(x, invalid_below = 5, method = "none")
    expect_setequal(result[1:2], NA_real_)
    expect_equal(result[3:10], x[3:10])

    # Test invalid_above only (inclusive)
    result <- replace_invalid(x, invalid_above = 100, method = "none")
    expect_setequal(result[7:8], NA_real_)
    expect_equal(result[c(1:5, 9:10)], x[c(1:5, 9:10)])

    # Test both invalid_below and invalid_above
    result <- replace_invalid(
        x,
        invalid_below = 5,
        invalid_above = 100,
        method = "none"
    )
    expect_setequal(result[c(1:2, 7:8)], NA_real_)
    expect_equal(result[c(3:5, 9:10)], x[c(3:5, 9:10)])

    ## combines exact values and ranges
    result <- replace_invalid(
        x,
        invalid_values = 50.5,
        invalid_below = 5,
        invalid_above = 100,
        method = "none"
    )
    expect_setequal(result[c(1:2, 6, 7:8)], NA_real_)
    expect_equal(result[c(3:5, 9:10)], x[c(3:5, 9:10)])
})

test_that("replace_invalid errors when no invalid criteria specified", {
    x <- c(1, 2, 3)

    expect_error(
        replace_invalid(
            x,
            invalid_values = NULL,
            invalid_above = NULL,
            invalid_below = NULL
        ),
        "At least one of.*must be specified"
    )
})

test_that("replace_invalid handles overlapping conditions", {
    x <- c(-999, 0, 5, -999, 10)

    # Value meets both invalid_values and invalid_below criteria
    result <- replace_invalid(
        x,
        invalid_values = -999,
        invalid_below = 5,
        method = "none"
    )
    expect_setequal(is.na(result[1:4]), TRUE)
})


test_that("replace_invalid() works correctly", {
    x <- c(1, 2, 3, 16, 5, 6, 7)
    expect_equal(
        replace_invalid(x, invalid_values = 16, method = "none")[4],
        NA_real_
    )
    expect_equal(
        replace_invalid(x, invalid_values = 16, width = 5, method = "median")[
            4
        ],
        median(x[c(1:3, 5:7)])
    )

    ## no invalid
    x_valid <- 1:7
    result_clean <- replace_invalid(x_valid, invalid_values = 16, width = 5)
    expect_equal(result_clean, x_valid)

    ## edge cases
    expect_equal(replace_invalid(c(1), invalid_values = 16, method = "none"), 1)
    expect_equal(
        sum(replace_invalid(rep(1, 7), invalid_values = 16, method = "none")),
        7
    )
    expect_error(
        replace_invalid(rep(NA, 7), invalid_values = 16),
        "`x` must be a"
    )

    ## NA handling
    x_na <- c(1, 2, 3, NA, 5, 6, 7)
    expect_equal(
        replace_invalid(x_na, invalid_values = 35, method = "none"),
        x_na
    )
    expect_equal(
        replace_invalid(
            x_na,
            invalid_values = 35,
            width = 3,
            method = "median"
        ),
        x_na
    )
    expect_true(
        is.na(replace_invalid(x_na, invalid_values = 35, method = "none")[4])
    )

    x_na <- c(1, 2, NA, 35, 5, 6, 7)
    expect_true(
        all(is.na(replace_invalid(x_na, invalid_values = 35, method = "none")[
            3:4
        ]))
    )

    expect_equal(
        replace_invalid(
            x_na,
            invalid_values = 35,
            width = 3,
            method = "median"
        )[4],
        median(x_na[c(2:3, 5:6)], na.rm = TRUE)
    )
})


## replace_missing() =================================================
test_that("replace_missing() validates inputs", {
    # x must be numeric
    expect_error(
        replace_missing(c("a", NA, "c")),
        "x.*numeric"
    )

    # t must be numeric
    expect_error(
        replace_missing(c(1, NA, 3), t = c("a", "b", "c")),
        "t.*numeric"
    )

    # x and t must be equal length
    expect_error(
        replace_missing(c(1, NA, 3), t = c(1, 2)),
        "x.*t.*numeric.*equal length"
    )

    ## missing width/span
    expect_error(
        replace_missing(c(1, NA, 3), method = "median"),
        "width.*span.*must be defined"
    )
})

test_that("replace_missing() returns unchanged vector when no NAs", {
    x <- c(1, 2, 3, 4, 5)
    expect_equal(replace_missing(x), x)
    expect_equal(replace_missing(x, method = "locf"), x)
    expect_equal(replace_missing(x, method = "median", width = 1), x)
})

test_that("replace_missing() replaces NAs with linear method", {
    x <- c(1, 2, NA, 4, 5)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(1, 2, 3, 4, 5))

    # Multiple consecutive NAs
    x <- c(1, NA, NA, 4)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(1, 2, 3, 4))

    # Leading NA left extrapolation
    x <- c(NA, NA, 2, 3, 4)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(2, 2, 2, 3, 4))

    # Trailing NA right extrapolation = "locf"
    x <- c(1, 2, 3, NA, NA)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(1, 2, 3, 3, 3))

    # Both leading and trailing NAs
    x <- c(NA, 2, NA, 4, NA)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(2, 2, 3, 4, 4))
})

test_that("replace_missing() replaces NAs with locf method", {
    x <- c(1, 2, NA, NA, 5)
    result <- replace_missing(x, method = "locf")
    expect_equal(result, c(1, 2, 2, 2, 5))

    # Leading NA (next observation carried backward)
    x <- c(NA, NA, 3, 4)
    result <- replace_missing(x, method = "locf")
    expect_equal(result, c(3, 3, 3, 4))

    # Trailing NA (last observation carried forward)
    x <- c(1, 2, NA, NA)
    result <- replace_missing(x, method = "locf")
    expect_equal(result, c(1, 2, 2, 2))
})

test_that("replace_missing() replaces NAs with median method", {
    x <- c(1, 2, NA, 4, 5)
    result <- replace_missing(x, method = "median", width = 1)
    expect_equal(result, c(1, 2, 3, 4, 5)) # median(2, 4) = 3

    # width = 3
    x <- c(1, 2, 3, NA, 5, 6, 7)
    result <- replace_missing(x, method = "median", width = 3)
    expect_equal(result, c(1, 2, 3, 4, 5, 6, 7)) # median(2,3,5,6) = 4

    # Multiple NAs
    x <- c(1, NA, 3, NA, 5)
    result <- replace_missing(x, method = "median", width = 1)
    expect_equal(result, c(1, 2, 3, 4, 5))
})

test_that("replace_missing() replaces NAs with span method", {
    x <- c(1, 2, NA, 4, 5)
    t <- c(0, 1, 2, 3, 4)
    expect_equal(
        replace_missing(x, t = t, method = "median", span = 1),
        c(1, 2, 3, 4, 5)
    )
    expect_equal(
        replace_missing(x, t = t, method = "median", span = 0.5),
        c(1, 2, 3, 4, 5)
    )
})

test_that("replace_missing() replaces NaN", {
    x <- c(1, NaN, 3, NaN, 5)
    result <- replace_missing(x, method = "linear")
    expect_equal(result, c(1, 2, 3, 4, 5))
    expect_false(any(is.nan(result)))
})

test_that("replace_missing() works with custom time vector", {
    x <- c(10, NA, 30)
    t <- c(0, 5, 10)
    result <- replace_missing(x, t = t, method = "linear")
    expect_equal(result, c(10, 20, 30))
})

test_that("replace_missing() handles all NAs", {
    x <- c(NA, NA, NA)

    expect_error(replace_missing(x, method = "linear"), "x.*numeric")
})


## replace_mnirs() =================================================
test_that("replace_mnirs outlier removal skipped when outlier_cutoff = NULL", {
    data <- data.frame(
        time = 1:5,
        ch1 = c(50, 51, 200, 53, 54)
    )
    class(data) <- c("mnirs", "data.frame")
    attr(data, "nirs_channels") <- "ch1"
    attr(data, "time_channel") <- "time"

    ## Store original to verify no outlier processing occurred
    original_ch1 <- data$ch1

    result <- replace_mnirs(
        data,
        invalid_values = c(999),
        method = "none",
        verbose = FALSE
    )

    ## Data should be unchanged except for invalid value processing
    expect_equal(result$ch1, original_ch1)
})

test_that("replace_mnirs outlier removal processes when outlier_cutoff provided", {
    data <- data.frame(
        time = 1:10,
        ch1 = c(50, 51, 52, 200, 54, 55, 56, 57, 58, 59),
        ch2 = c(60, 61, 62, 63, 64, 65, 250, 67, 68, 69)
    )
    class(data) <- c("mnirs", "data.frame")
    attr(data, "nirs_channels") <- c("ch1", "ch2")
    attr(data, "time_channel") <- "time"

    result <- replace_mnirs(
        data,
        outlier_cutoff = 3,
        width = 5,
        method = "none",
        verbose = FALSE
    )

    ## Outlier (200) should be replaced with NA
    expect_true(is.na(result$ch1[4]))
    expect_equal(result$ch1[-4], data$ch1[-4])
    expect_true(is.na(result$ch2[7]))
    expect_equal(result$ch2[-7], data$ch2[-7])
})

test_that("replace_mnirs do nothing condition throws error appropriately", {
    data <- data.frame(
        time = 1:5,
        ch1 = 50:54
    )
    class(data) <- c("mnirs", "data.frame")
    attr(data, "nirs_channels") <- "ch1"
    attr(data, "time_channel") <- "time"

    expect_error(
        replace_mnirs(
            data,
            invalid_values = NULL,
            width = NULL,
            span = NULL,
            method = "none",
            verbose = FALSE
        ),
        "must be specified"
    )
})

test_that("replace_mnirs updates metadata correctly", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        replace_mnirs(
            invalid_values = c(0, 100),
            span = 7,
            method = "linear",
            verbose = FALSE
        )
    expect_equal(attr(data, "nirs_channels"), c("smo2_left", "smo2_right"))
    expect_equal(attr(data, "time_channel"), "time")
    expect_equal(attr(data, "sample_rate"), 2)

    ## overwrites nirs_channels
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        replace_mnirs(
            nirs_channels = "smo2_left",
            invalid_values = c(0, 100),
            span = 7,
            method = "linear",
            verbose = FALSE
        )
    expect_equal(attr(data, "nirs_channels"), "smo2_left")
})

test_that("replace_mnirs processes a named list of data frames", {
    make_df <- \(vals) {
        create_mnirs_data(
            data.frame(time = 1:5, ch1 = vals),
            nirs_channels = "ch1",
            time_channel = "time"
        )
    }
    data_list <- list(
        a = make_df(c(50, 999, 52, 53, 54)),
        b = make_df(c(60, 61, 999, 63, 64))
    )

    result <- replace_mnirs(
        data_list,
        invalid_values = 999,
        method = "linear",
        verbose = FALSE
    )

    expect_type(result, "list")
    expect_named(result, c("a", "b"))
    expect_s3_class(result$a, "mnirs")
    expect_s3_class(result$b, "mnirs")
    expect_equal(result$a$ch1, c(50, 51, 52, 53, 54))
    expect_equal(result$b$ch1, c(60, 61, 62, 63, 64))
})

test_that("replace_mnirs names unnamed list intervals", {
    df <- create_mnirs_data(
        data.frame(time = 1:5, ch1 = c(50, 999, 52, 53, 54)),
        nirs_channels = "ch1",
        time_channel = "time"
    )

    ## explicit `nirs_channels` forwarded to each interval
    result <- replace_mnirs(
        list(df, df),
        nirs_channels = "ch1",
        invalid_values = 999,
        method = "linear",
        verbose = FALSE
    )

    expect_named(result, c("interval_1", "interval_2"))
    expect_equal(result$interval_1$ch1, c(50, 51, 52, 53, 54))
})

test_that("replace_mnirs processes grouped data frames", {
    skip_if_not_installed("dplyr")

    df <- create_mnirs_data(
        data.frame(
            time = rep(1:5, 2),
            ch1 = c(50, 999, 52, 53, 54, 60, 61, 999, 63, 64),
            group = rep(c("A", "B"), each = 5)
        ),
        nirs_channels = "ch1",
        time_channel = "time",
        sample_rate = 1
    ) |> 
        dplyr::group_by(group)

    result <- replace_mnirs(
        df,
        invalid_values = 999,
        method = "linear",
        verbose = FALSE
    )

    expect_type(result, "list")
    expect_named(result, c("A", "B"))
    expect_s3_class(result$A, "mnirs")
    expect_equal(result$A$ch1, c(50, 51, 52, 53, 54))
    expect_equal(result$B$ch1, c(60, 61, 62, 63, 64))
    ## metadata copied down to each interval
    expect_equal(attr(result$A, "sample_rate"), 1)
    expect_equal(attr(result$B, "sample_rate"), 1)
})

test_that("replace_mnirs errors on list with non-data frame elements", {
    df <- create_mnirs_data(
        data.frame(time = 1:5, ch1 = 50:54),
        nirs_channels = "ch1",
        time_channel = "time"
    )

    expect_error(
        replace_mnirs(
            list(df, "not a df"),
            invalid_values = 999,
            verbose = FALSE
        ),
        "must be a list of data frames"
    )
})

test_that("replace_mnirs global verbose works", {
    expect_warning(
        read_mnirs(
            file_path = example_mnirs("moxy_ramp.xlsx"),
            nirs_channels = c(smo2 = "SmO2 Live"),
            time_channel = c(time = "hh:mm:ss")
        ) |> 
            replace_mnirs(),
        "irregular.*detected"
    ) |>
        expect_message("Estimated.*sample_rate.*2")

    old_verbose <- getOption("mnirs.verbose")
    on.exit(options(mnirs.verbose = old_verbose), add = TRUE)
    options(mnirs.verbose = FALSE)

    expect_silent(
        read_mnirs(
            file_path = example_mnirs("moxy_ramp.xlsx"),
            nirs_channels = c(smo2 = "SmO2 Live"),
            time_channel = c(time = "hh:mm:ss")
        ) |>
            replace_mnirs()
    )
})

test_that("replace_mnirs works visually on moxy data", {
    skip_if_not_installed("ggplot2")
    skip("visual check of ggplots")
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2 = "SmO2 Live"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )

    plot(data) +
        ggplot2::ylim(0, 100)

    replace_mnirs(
        data,
        invalid_values = c(0, 100),
        span = 7,
        method = "none"
    ) |>
        plot() +
        ggplot2::ylim(0, 100)

    replace_mnirs(
        data,
        invalid_values = c(0, 100),
        span = 7,
        method = "median"
    ) |>
        plot() +
        ggplot2::ylim(0, 100)

    replace_mnirs(
        data,
        invalid_values = c(0, 100),
        span = 7,
        method = "locf"
    ) |>
        plot() +
        ggplot2::ylim(0, 100)

    replace_mnirs(
        data,
        invalid_values = c(0, 100),
        span = 7,
        method = "linear"
    ) |>
        plot() +
        ggplot2::ylim(0, 100)
})
