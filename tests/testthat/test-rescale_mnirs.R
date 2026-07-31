test_that("rescale_mnirs rescales single channel correctly", {
    data <- tibble(A = c(0, 50, 100), B = c(10, 20, 30))

    result <- rescale_mnirs(data, nirs_channels = "A", range = c(0, 1))

    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, data$B) # unchanged
})

test_that("rescale_mnirs preserves relative scaling across grouped channels", {
    data <- tibble(A = c(0, 50, 100), B = c(25, 25, 50))

    result <- rescale_mnirs(data, nirs_channels = c("A", "B"), range = c(0, 1))

    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, c(0.25, 0.25, 0.5))
})

test_that("rescale_mnirs handles multiple separate groups", {
    data <- tibble(A = c(10, 100), B = c(0, 50), C = c(10, 200))

    result <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B", "C"),
        group_channels = list(c("A", "B"), "C"),
        range = c(0, 10)
    )

    expect_equal(result$A, c(1, 10))
    expect_equal(result$B, c(0, 5))
    # C scaled independently
    expect_equal(result$C, c(0, 10))
})

test_that("rescale_mnirs handles negative ranges", {
    data <- tibble(A = c(0, 50, 100))

    result <- rescale_mnirs(data, nirs_channels = "A", range = c(-1, 1))

    expect_equal(result$A, c(-1, 0, 1))
})

test_that("rescale_mnirs handles NA values", {
    data <- tibble(A = c(0, NA, 100))

    result <- rescale_mnirs(data, nirs_channels = "A", range = c(0, 1))

    expect_equal(result$A, c(0, NA, 1))
})

test_that("rescale_mnirs aborts on invalid channel selection", {
    data <- tibble(A = c(0, 50, 100))

    ## no channels detected from data or metadata
    expect_error(
        rescale_mnirs(data, nirs_channels = NULL, range = c(0, 1)),
        "nirs_channels.*not detected"
    )

    expect_error(
        rescale_mnirs(data, nirs_channels = "doesn't exist", range = c(0, 1)),
        "nirs_channels.*match exactly"
    )
})

test_that("rescale_mnirs aborts on invalid range", {
    data <- tibble(A = c(0, 50, 100))

    ## range must be a two-element numeric vector
    expect_error(
        rescale_mnirs(data, nirs_channels = "A", range = c(0, 1, 2)),
        "range.*numeric"
    )
})

test_that("rescale_mnirs passes through a constant channel", {
    data <- tibble(A = c(50, 50, 50), B = c(0, 100, 200))

    ## distinct: constant channel unchanged, varying channel rescales
    distinct <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B"),
        group_channels = "distinct",
        range = c(0, 100)
    )
    expect_equal(distinct$A, data$A)
    expect_equal(distinct$B, data$B / 2)

    ## grouped: pooled range is non-zero, so the constant channel is
    ## rescaled together with its group rather than passed through
    grouped <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B"),
        range = c(0, 100)
    )
    expect_equal(grouped$A, data$A / 2)
    expect_equal(grouped$B, data$B / 2)
})

test_that("rescale_mnirs passes through an all-constant group", {
    data <- tibble(A = c(50, 50, 50), B = c(20, 20, 20))

    ## distinct grouped constant values return unchanged
    result <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B"),
        group_channels = "distinct",
        range = c(0, 100)
    )
    expect_equal(result$A, data$A)
    expect_equal(result$B, data$B)
})

test_that("rescale_mnirs accepts a single-column data frame", {
    ## rescale requires only one column, unlike shift's two
    data <- tibble(A = c(0, 50, 100))

    result <- rescale_mnirs(data, nirs_channels = "A", range = c(0, 1))
    expect_equal(result$A, c(0, 0.5, 1))
})

test_that("rescale_mnirs errors on deprecated nirs_channels = list()", {
    data <- tibble(A = c(0, 50, 100))

    ## passing a list() for channel grouping is deprecated; use group_channels
    expect_error(
        rescale_mnirs(data, nirs_channels = list("A"), range = c(0, 1)),
        "group_channels"
    )
})

test_that("rescale_mnirs applies per-group range keyed by group name", {
    data <- tibble(
        A = c(0, 50, 100),
        B = c(0, 25, 50),
        C = c(0, 10, 20)
    )
    result <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B", "C"),
        group_channels = list(smo2 = c("A", "B"), hhb = "C"),
        range = list(smo2 = c(0, 1), hhb = c(0, 10))
    )

    ## smo2 group pooled to [0, 1]; hhb group to [0, 10]
    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, c(0, 0.25, 0.5))
    expect_equal(result$C, c(0, 5, 10))
})

test_that("rescale_mnirs applies per-channel range with unnamed fallback", {
    data <- tibble(
        A = c(0, 50, 100),
        B = c(0, 50, 100)
    )
    result <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B"),
        group_channels = "distinct",
        range = list(c(0, 100), A = c(0, 1))
    )

    ## A takes its own range; B takes the unnamed fallback
    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, c(0, 50, 100))
})

test_that("rescale_mnirs aborts on intra-group range conflict", {
    data <- tibble(A = c(0, 50, 100), B = c(0, 25, 50))

    expect_error(
        rescale_mnirs(
            data,
            nirs_channels = c("A", "B"),
            group_channels = list(g = c("A", "B")),
            range = list(A = c(0, 1), B = c(0, 2))
        ),
        "conflicting"
    )
})

test_that("rescale_mnirs updates metadata correctly", {
    data <- tibble(A = c(50, 50, 50), B = c(0, 100, 200))

    result <- rescale_mnirs(
        data,
        nirs_channels = c("A", "B"),
        group_channels = "distinct",
        range = c(0, 100),
        verbose = FALSE
    )

    expect_true(all(c("A", "B") %in% attr(result, "nirs_channels")))
})

test_that("rescale_mnirs works on Moxy", {
    file_path <- example_mnirs("moxy_ramp.xlsx")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        dplyr::mutate(
            dplyr::across(
                dplyr::matches("smo2"),
                \(.x) {
                    replace_invalid(
                        .x,
                        invalid_values = c(0, 100),
                        method = "none"
                    )
                }
            )
        )

    result <- rescale_mnirs(
        df,
        nirs_channels = c("smo2_left", "smo2_right"),
        range = c(0, 100)
    )

    # plot(df) + 
    #     ggplot2::ylim(0, 100) + 
    #     ggplot2::geom_hline(yintercept = c(0, 100))
    # plot(result) +
    #     ggplot2::ylim(0, 100) +
    #     ggplot2::geom_hline(yintercept = c(0, 100))

    ## check grouping together: min & max value should come from smo2_right
    expect_false(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))

    ## check grouped apart
    result <- rescale_mnirs(
        df,
        nirs_channels = c("smo2_left", "smo2_right"),
        group_channels = "distinct",
        range = c(0, 100)
    )

    ## check grouping together: min value should come from smo2_right
    expect_true(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))
})

test_that("rescale_mnirs works on Train.Red", {
    file_path <- example_mnirs("train.red_intervals.csv")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered",
            o2hb_left = "O2HB unfiltered",
            o2hb_right = "O2HB unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        verbose = FALSE
    )

    result <- rescale_mnirs(
        df,
        nirs_channels = c("smo2_left", "smo2_right", "o2hb_left", "o2hb_right"),
        group_channels = list(
            "smo2_left",
            "smo2_right",
            c("o2hb_left", "o2hb_right")
        ),
        range = c(0, 100),
    )

    # plot(df) +
    #     ggplot2::ylim(0, 100) +
    #     ggplot2::geom_hline(yintercept = c(0, 100))
    # plot(result) +
    #     ggplot2::ylim(0, 100) +
    #     ggplot2::geom_hline(yintercept = c(0, 100))

    ## distinct smo2 channels both reach their own min and max
    expect_true(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))

    ## o2hb channels pooled: group min and max come from opposite channels
    o2hb_min <- which.min(c(min(df$o2hb_left), min(df$o2hb_right)))
    o2hb_max <- which.max(c(max(df$o2hb_left), max(df$o2hb_right)))
    o2hb_result <- list(result$o2hb_left, result$o2hb_right)
    expect_true(any(o2hb_result[[o2hb_min]] == 0, na.rm = TRUE))
    expect_true(any(o2hb_result[[o2hb_max]] == 100, na.rm = TRUE))
})


## multi-interval input ================================================
test_that("rescale_mnirs processes a named list of data frames", {
    make_df <- \(vals) {
        create_mnirs_data(
            data.frame(time = 1:5, ch1 = vals),
            nirs_channels = "ch1",
            time_channel = "time"
        )
    }
    data_list <- list(a = make_df(50:54), b = make_df(60:64))

    result <- rescale_mnirs(
        data_list,
        group_channels = "distinct",
        range = c(0, 100),
        verbose = FALSE
    )

    expect_type(result, "list")
    expect_named(result, c("a", "b"))
    expect_s3_class(result$a, "mnirs")
    ## each interval rescaled to its own range
    expect_equal(result$a$ch1, c(0, 25, 50, 75, 100))
    expect_equal(result$b$ch1, c(0, 25, 50, 75, 100))
})

test_that("rescale_mnirs processes grouped data frames", {
    skip_if_not_installed("dplyr")

    df <- create_mnirs_data(
        data.frame(
            time = rep(1:5, 2),
            ch1 = c(50:54, 60:64),
            group = rep(c("A", "B"), each = 5)
        ),
        nirs_channels = "ch1",
        time_channel = "time"
    )
    grouped_df <- dplyr::group_by(df, group)

    result <- rescale_mnirs(
        grouped_df,
        group_channels = "distinct",
        range = c(0, 100),
        verbose = FALSE
    )

    expect_named(result, c("A", "B"))
    expect_s3_class(result$A, "mnirs")
    expect_equal(result$A$ch1, c(0, 25, 50, 75, 100))
    expect_equal(result$B$ch1, c(0, 25, 50, 75, 100))
})
