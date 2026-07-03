## helper to create test data with metadata
create_test_data <- function(
    time_max = 10,
    sample_rate = 10,
    add_metadata = TRUE
) {
    time <- seq(0, time_max, by = 1 / sample_rate)
    nrow <- length(time)
    data <- tibble(
        time = time,
        nirs1 = rnorm(nrow, 50, 5),
        nirs2 = rnorm(nrow, 60, 5),
        nirs3 = rnorm(nrow, 80, 5),
        event = c(1, rep(NA, nrow - 2), 2),
    )
    class(data) <- c("mnirs", class(data))

    if (add_metadata) {
        attr(data, "time_channel") <- "time"
        attr(data, "nirs_channels") <- c("nirs1", "nirs2")
        attr(data, "event_channel") <- "event"
        attr(data, "sample_rate") <- sample_rate
    }

    return(data)
}

test_that("shift_mnirs requires either to or by", {
    data <- tibble(
        time = 1:10,
        ch1 = rnorm(10)
    )

    expect_error(
        shift_mnirs(data, nirs_channels = "ch1", time_channel = "time"),
        "Shift value undefined"
    )
})

test_that("shift_mnirs shifts by constant correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = "distinct",
        by = 5,
        verbose = FALSE
    )

    expect_equal(result$ch1, data$ch1 + 5)
    expect_equal(result$ch2, data$ch2 + 5)
})

test_that("shift_mnirs position = 'min' preserves scaling per grouping", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    ## ensemble: both channels shifted by the shared minimum (ch1's)
    ensemble <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )
    expect_equal(ensemble$ch1, 0:9)
    expect_equal(ensemble$ch2, 10:19)
    ## relative scaling between channels preserved
    expect_equal(ensemble$ch2[1] - ensemble$ch1[1], 10)

    ## distinct: each channel shifted independently to its own minimum
    distinct <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = "distinct",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )
    expect_equal(distinct$ch1, 0:9)
    expect_equal(distinct$ch2, 0:9)
    ## relative scaling lost between channels
    expect_equal(distinct$ch2[1] - distinct$ch1[1], 0)
})

test_that("shift_mnirs handles position = 'max' correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10
    )

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        to = 100,
        width = 1,
        position = "max",
        verbose = FALSE
    )

    ## max is 10, shift to 100 means add 90
    expect_equal(result$ch1, 91:100)
})

test_that("shift_mnirs handles position = 'first' with span and width", {
    data <- tibble(
        time = seq(0, 9, by = 1),
        ch1 = c(rep(10, 3), 1:7)
    )

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        to = 0,
        position = "first",
        span = 2,
        verbose = FALSE
    )

    ## mean of first 3 values (time 0-2) is 10
    expect_equal(result$ch1, c(rep(0, 3), -9:-3))

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        to = 0,
        position = "first",
        width = 3,
        verbose = FALSE
    )
    expect_equal(result$ch1, c(rep(0, 3), -9:-3))
})

test_that("shift_mnirs preserves non-channel columns", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        other = letters[1:10]
    )

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        by = 5,
        width = 1,
        verbose = FALSE
    )

    ## selected channel shifted; unselected numeric and character pass through
    expect_equal(result$ch1, 6:15)
    expect_equal(result$ch2, 11:20)
    expect_equal(result$other, letters[1:10])
    expect_equal(result$time, 1:10)
})

test_that("shift_mnirs updates metadata correctly", {
    data <- tibble(time = 1:10, ch1 = 1:10, ch2 = 1:10)

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        by = 5,
        width = 1,
        verbose = FALSE
    )

    expect_true("ch1" %in% attr(result, "nirs_channels"))
    expect_equal(attr(result, "time_channel"), "time")
})

test_that("shift_mnirs handles multiple channel groups", {
    data <- data.frame(
        time = 1:2,
        ch1 = c(10, 20),
        ch2 = c(15, 25),
        ch3 = c(5, 35)
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2", "ch3"),
        time_channel = time,
        group_channels = list(c("ch1", "ch2"), "ch3"),
        to = 0,
        width = 1
    )

    ## grouped together: shared min comes from ch1; ch3 in its own group
    expect_true(any(result$ch1 == 0, na.rm = TRUE))
    expect_false(any(result$ch2 == 0, na.rm = TRUE))
    expect_true(any(result$ch3 == 0, na.rm = TRUE))
    ## ch1 & ch2 shifted together, preserving relative scaling
    expect_equal(result$ch1 - result$ch2, data$ch1 - data$ch2)
    ## ch3 shifted independently of ch1
    expect_false(isTRUE(all.equal(
        result$ch1 - result$ch3,
        data$ch1 - data$ch3
    )))
})

test_that("shift_mnirs handles NA values correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = c(NA, 2:10)
    )

    result <- shift_mnirs(
        data,
        nirs_channels = "ch1",
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    ## min of non-NA values is 2; NA preserved
    expect_true(is.na(result$ch1[1]))
    expect_equal(result$ch1[2:10], 0:8)
})

test_that("shift_mnirs errors on deprecated nirs_channels = list()", {
    data <- tibble(time = 1:5, ch1 = 1:5)

    ## passing a list() for channel grouping is deprecated; use group_channels
    expect_error(
        shift_mnirs(
            data,
            nirs_channels = list("ch1"),
            time_channel = "time",
            by = 1
        ),
        "group_channels"
    )
})

test_that("shift_mnirs aborts on invalid channel selection", {
    data <- data.frame(time = 1:3, value = c(10, 20, 30))

    ## no channels detected from data or metadata
    expect_error(
        shift_mnirs(data, nirs_channels = NULL, by = 5),
        "not detected"
    )

    ## named column that does not exist
    data <- create_test_data()
    expect_error(
        shift_mnirs(
            data,
            nirs_channels = nonexistent,
            time_channel = time,
            to = 0,
            width = 5,
            verbose = FALSE
        ),
        "not detected"
    )
})

test_that("shift_mnirs validates position argument", {
    data <- data.frame(
        time = 1:3,
        ch1 = c(10, 20, 30)
    )
    expect_error(
        shift_mnirs(
            data,
            "ch1",
            "time",
            to = 0,
            width = 1,
            position = "invalid"
        ),
        "must be one of"
    )
})


## per-channel & per-group argument overrides ==============================

test_that("shift_mnirs applies per-group `to` keyed by group name", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        ch3 = 21:30
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2", "ch3"),
        time_channel = "time",
        group_channels = list(smo2 = c("ch1", "ch2"), hhb = "ch3"),
        to = list(smo2 = 0, hhb = 5),
        width = 1,
        position = "min",
        verbose = FALSE
    )

    ## smo2 group shifted so its shared min (ch1's 1) lands at 0
    expect_equal(result$ch1[1], 0)
    expect_equal(result$ch2[1], 10)
    ## hhb group shifted so ch3's min (21) lands at 5
    expect_equal(result$ch3[1], 5)
})

test_that("shift_mnirs applies per-group `to` keyed by member channel", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = list(smo2 = c("ch1", "ch2")),
        to = list(ch2 = 3),
        width = 1,
        position = "min",
        verbose = FALSE
    )

    ## a member key applies to the whole group: shared min (1) shifted to 3
    expect_equal(result$ch1[1], 3)
    expect_equal(result$ch2[1], 13)
})

test_that("shift_mnirs applies per-channel `by` with unnamed fallback", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        ch3 = 21:30
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2", "ch3"),
        time_channel = "time",
        group_channels = "distinct",
        by = list(5, ch1 = 7),
        verbose = FALSE
    )

    ## ch1 takes its own value; ch2 & ch3 take the unnamed fallback
    expect_equal(result$ch1, data$ch1 + 7)
    expect_equal(result$ch2, data$ch2 + 5)
    expect_equal(result$ch3, data$ch3 + 5)
})

test_that("shift_mnirs applies per-group `width` override", {
    data <- tibble(
        time = 1:6,
        ch1 = c(10, 10, 1, 1, 1, 1),
        ch2 = c(20, 20, 20, 2, 2, 2)
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = "distinct",
        to = 0,
        position = "first",
        width = list(ch1 = 2, ch2 = 3),
        verbose = FALSE
    )

    ## each channel uses its own first-window length
    expect_equal(result$ch1[1], data$ch1[1] - mean(data$ch1[1:2]))
    expect_equal(result$ch2[1], data$ch2[1] - mean(data$ch2[1:3]))
})

test_that("shift_mnirs applies exclusive args per channel", {
    data <- tibble(
        time = seq(0, len = 6, by = 0.5),
        ch1 = c(10, 10, 1, 1, 1, 1),
        ch2 = c(20, 20, 20, 2, 2, 2)
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = "distinct",
        to = list(ch2 = 0),
        by = list(ch1 = -10),
        position = "first",
        width = list(ch1 = 2),
        span = list(ch2 = 2),
        verbose = FALSE
    )

    expect_all_equal(result$ch1[1:2], 0)
    expect_equal(result$ch2, data$ch2 - mean(data$ch2[1:5]))
})

test_that("shift_mnirs applies per-group `position` override", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2"),
        time_channel = "time",
        group_channels = list(g1 = "ch1", g2 = "ch2"),
        to = 0,
        width = 1,
        position = list(g1 = "min", g2 = "max"),
        verbose = FALSE
    )

    ## g1 shifts its minimum (1) to 0; g2 shifts its maximum (20) to 0
    expect_equal(result$ch1, 0:9)
    expect_equal(result$ch2, -9:0)
})

test_that("shift_mnirs aborts on intra-group argument conflict", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )
    expect_error(
        shift_mnirs(
            data,
            nirs_channels = c("ch1", "ch2"),
            time_channel = "time",
            group_channels = list(g = c("ch1", "ch2")),
            to = list(ch1 = 0, ch2 = 5),
            width = 1,
            position = "min",
            verbose = FALSE
        ),
        "conflicting"
    )
})

test_that("shift_mnirs informs when channel omitted from args", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    #! FIX LOTS OF EDGE CASES?
    shift_mnirs(
        data,
        nirs_channels = c(ch1, ch2),
        time_channel = time,
        group_channels = "distinct",
        to = list(ch1 = 0),
        width = list(ch1 = 1),
        position = "min"
    ) |> 
        expect_warning("`to`:.*ch2.*not specified") |>
        expect_warning("`width`:.*ch2.*not specified") |>
        expect_error("Shift value undefined")

    ## unrecognised names warned and ignored; all channels omitted from
    ## `to` leaves the shift undefined
    shift_mnirs(
        data,
        nirs_channels = c(ch1, ch2),
        time_channel = time,
        group_channels = "distinct",
        to = list(ch3 = 0),
        width = 1,
        position = "min"
    ) |>
        expect_warning("`to`:.*ch3.*not recognised") |>
        expect_warning("`to`:.*not specified") |>
        expect_error("Shift value undefined")
})

test_that("shift_mnirs `to` overrides `by` once, for the first group", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    ## both groups supplied to & by; message must fire once (first group only)
    messages <- capture_messages(
        shift_mnirs(
            data,
            nirs_channels = c("ch1", "ch2"),
            time_channel = "time",
            group_channels = list(g1 = "ch1", g2 = "ch2"),
            to = list(g1 = 0, g2 = 10),
            by = list(g1 = 100, g2 = 200),
            width = 1,
            position = "min",
            verbose = TRUE
        )
    )
    expect_length(grep("overrides", messages), 1L)
})

test_that("shift_mnirs processes omitted channel as its own group", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        ch3 = 31:40
    )
    result <- shift_mnirs(
        data,
        nirs_channels = c("ch1", "ch2", "ch3"),
        time_channel = "time",
        group_channels = list(c("ch1", "ch2")),
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    ## ch3, omitted from the custom group, shifted to its own minimum
    expect_equal(result$ch3[1], 0)
    ## ch1 & ch2 share their group minimum
    expect_equal(result$ch1[1], 0)
    expect_equal(result$ch2[1], 10)
})


## unevenly sampled data ==================================================
test_that("shift_mnirs position = 'first' handles uneven sampling", {
    ## 0.5 Hz target (2 second intervals) with deterministic jitter
    data <- data.frame(
        time = c(0.05, 1.93, 4.08, 5.92, 8.10, 9.95),
        ch1 = c(10, 5, 20, 15, 25, 30)
    )

    ## 4-second span window from the first sample
    result <- shift_mnirs(
        data,
        "ch1",
        time_channel = "time",
        to = 0,
        span = 4,
        position = "first",
        verbose = FALSE
    )

    first_window_idx <- which(data$time <= data$time[1] + 4)
    expected_mean <- mean(data$ch1[first_window_idx])
    expect_equal(result$ch1[1], data$ch1[1] - expected_mean)

    ## width = 3 averages the first three samples regardless of timing
    result <- shift_mnirs(
        data,
        "ch1",
        time_channel = "time",
        to = 0,
        width = 3,
        position = "first",
        verbose = FALSE
    )
    expect_equal(result$ch1[1], data$ch1[1] - mean(data$ch1[1:3]))
})

test_that("shift_mnirs position = 'min' shifts argmin window to `to`", {
    data <- data.frame(
        time = c(0.05, 1.93, 4.08, 5.92, 8.10, 9.95),
        ch1 = c(10, 5, 20, 15, 25, 30)
    )

    result <- shift_mnirs(
        data,
        "ch1",
        time_channel = "time",
        to = 0,
        span = 4,
        position = "min",
        verbose = FALSE
    )

    ## recreate the local-window means and find the argmin window
    t <- data$time
    window_means <- vapply(seq_along(t), \(.i) {
        ## centred span = 4 window
        idx <- which(t >= t[.i] - 2 & t <= t[.i] + 2)
        mean(data$ch1[idx])
    }, numeric(1))
    argmin <- which.min(window_means)
    argmin_idx <- which(t >= t[argmin] - 2 & t <= t[argmin] + 2)

    ## the minimum local-window mean is shifted to `to` = 0
    expect_equal(mean(result$ch1[argmin_idx]), 0)
    expect_equal(result$ch1, data$ch1 - window_means[argmin])
})


## tidy eval integration ==================================================

test_that("shift_mnirs accepts channels as strings, symbols, vectors, select", {
    data <- create_test_data()
    channels <- c("nirs1", "nirs2")
    time_col <- "time"
    expected_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))

    ## each input form selects the same channels; results must be identical
    results <- list(
        strings = shift_mnirs(
            data,
            nirs_channels = c("nirs1", "nirs2"),
            time_channel = "time",
            to = 0,
            width = 5,
            position = "first",
            verbose = FALSE
        ),
        symbols = shift_mnirs(
            data,
            nirs_channels = c(nirs1, nirs2),
            time_channel = time,
            to = 0,
            width = 5,
            position = "first",
            verbose = FALSE
        ),
        external = shift_mnirs(
            data,
            nirs_channels = channels,
            time_channel = time_col,
            to = 0,
            width = 5,
            position = "first",
            verbose = FALSE
        ),
        tidyselect = shift_mnirs(
            data,
            nirs_channels = tidyselect::starts_with("nirs"),
            time_channel = time,
            to = 0,
            width = 5,
            position = "first",
            verbose = FALSE
        )
    )

    lapply(results, \(.r) expect_s3_class(.r, "mnirs"))
    ## strings, symbols, external vector all select nirs1 & nirs2
    lapply(results[c("strings", "symbols", "external")], \(.r) {
        expect_equal(.r$nirs1[1], data$nirs1[1] - expected_mean)
        expect_equal(.r$nirs2[1], data$nirs2[1] - expected_mean)
    })
    expect_true(all(
        c("nirs1", "nirs2", "nirs3") %in% names(results$tidyselect)
    ))
    ## tidyselect picks up nirs3 too, changing the shared first mean
    select_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5], data$nirs3[1:5]))
    expect_equal(results$tidyselect$nirs1[1], data$nirs1[1] - select_mean)
})

test_that("shift_mnirs resolves group_channels from varied input forms", {
    data <- create_test_data()
    groups <- list(c("nirs1", "nirs2"))

    ## list of strings, list of symbols, and an external list object all work
    expect_s3_class(
        shift_mnirs(
            data,
            c("nirs1", "nirs2"),
            "time",
            group_channels = groups,
            to = 0,
            width = 5,
            verbose = FALSE
        ),
        "mnirs"
    )
    expect_s3_class(
        shift_mnirs(
            data,
            c(nirs1, nirs2),
            time,
            group_channels = list(c(nirs1, nirs2)),
            to = 0,
            width = 5,
            verbose = FALSE
        ),
        "mnirs"
    )
    ## mixed quoted/unquoted members
    result <- shift_mnirs(
        data,
        c("nirs1", nirs2),
        "time",
        group_channels = list("nirs1", nirs2),
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))

    ## "distinct" shifts each channel to its own first mean
    distinct <- shift_mnirs(
        data,
        c(nirs1, nirs2),
        time,
        group_channels = "distinct",
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_equal(distinct$nirs1[1], data$nirs1[1] - mean(data$nirs1[1:5]))
    expect_equal(distinct$nirs2[1], data$nirs2[1] - mean(data$nirs2[1:5]))
})

test_that("shift_mnirs() works with tidyselect in group_channels", {
    data <- data.frame(
        time = 1:10,
        smo2_left = runif(10, 50, 70),
        smo2_right = runif(10, 50, 70),
        thb = runif(10, 12, 14)
    )
    data <- create_mnirs_data(
        data,
        nirs_channels = c("smo2_left", "smo2_right", "thb"),
        time_channel = "time"
    )

    result <- shift_mnirs(
        data,
        nirs_channels = c(smo2_left, smo2_right, thb),
        time_channel = time,
        group_channels = list(tidyselect::starts_with("smo2"), thb),
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("smo2_left", "smo2_right", "thb") %in% names(result)))
})

test_that("shift_mnirs() uses metadata when channels NULL", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
})

test_that("shift_mnirs() preserves grouping with external group_channels", {
    data <- data.frame(
        time = 1:10,
        nirs1 = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
        nirs2 = c(15, 25, 35, 45, 55, 65, 75, 85, 95, 105),
        nirs3 = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95)
    )
    data <- create_mnirs_data(
        data,
        nirs_channels = c("nirs1", "nirs2", "nirs3"),
        time_channel = "time"
    )

    groups <- list(c("nirs1", "nirs2"), "nirs3")
    result <- shift_mnirs(
        data,
        nirs_channels = c("nirs1", "nirs2", "nirs3"),
        time_channel = "time",
        group_channels = groups,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )

    ## nirs1 & nirs2 grouped: shifted by the same amount
    ## nirs3 independent: shifted separately
    expect_equal(
        result$nirs1[1] - result$nirs2[1],
        data$nirs1[1] - data$nirs2[1]
    )
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1]),
        c(data$nirs1[1] - first_mean, data$nirs2[1] - first_mean)
    )
    expect_equal(result$nirs3[1], data$nirs3[1] - mean(c(data$nirs3[1:5])))
})


## integration tests ================================================
test_that("shift_mnirs works on Moxy", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
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

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = NULL,
        to = 0,
        by = NULL,
        span = 0,
        position = c("min", "max", "first"),
        verbose = FALSE
    )

    # plot(data) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))
    # plot(data_shifted) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))

    ## check grouping together: min value should come from smo2_right
    expect_false(any(data_shifted$smo2_left == 0, na.rm = TRUE))
    expect_true(any(data_shifted$smo2_right == 0, na.rm = TRUE))
    ## check both shifted together maintaining relative scaling
    expect_equal(
        data_shifted$smo2_left - data_shifted$smo2_right,
        data$smo2_left - data$smo2_right
    )
})

test_that("shift_mnirs(position = 'first') works on Moxy", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        zero_time = TRUE,
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
            ),
            dplyr::across(
                dplyr::matches("smo2"),
                \(.x) replace_missing(.x, )
            )
        )

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = NULL,
        to = 0,
        by = NULL,
        span = 120,
        position = "first",
        verbose = FALSE
    )

    first_mean <- data_shifted |>
        dplyr::filter(time <= 120) |>
        dplyr::summarise(
            mean = mean(c(smo2_left, smo2_right), na.rm = TRUE)
        ) |>
        dplyr::pull(mean)

    expect_equal(first_mean, 0)
})

test_that("shift_mnirs works on Train.Red", {
    data <- read_mnirs(
        file_path = example_mnirs("train.red_intervals.csv"),
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered",
            o2hb_left = "O2HB unfiltered",
            o2hb_right = "O2HB unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        verbose = FALSE
    )

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = c("smo2_left", "smo2_right", "o2hb_left", "o2hb_right"),
        time_channel = NULL,
        group_channels = list(
            "smo2_left",
            "smo2_right",
            c("o2hb_left", "o2hb_right")
        ),
        to = 0,
        by = NULL,
        span = 0,
        position = "min",
        verbose = FALSE
    )

    # plot(data) + ggplot2::ylim(0, 100)
    # plot(data_shifted) +
    #     ggplot2::ylim(0, 100) +
    #     ggplot2::geom_hline(yintercept = c(0))

    ## check grouping together: min value should come from each group
    expect_true(any(data_shifted$smo2_left == 0, na.rm = TRUE))
    expect_true(any(data_shifted$smo2_right == 0, na.rm = TRUE))
    expect_true(any(data_shifted$o2hb_left == 0, na.rm = TRUE))
    expect_false(any(data_shifted$o2hb_right == 0, na.rm = TRUE))
    ## check both shifted together maintaining relative scaling
    expect_equal(
        data_shifted$o2hb_left - data_shifted$o2hb_right,
        data$o2hb_left - data$o2hb_right
    )
    ## check both shifted independently
    expect_false(isTRUE(all.equal(
        data_shifted$smo2_left - data_shifted$smo2_right,
        data$smo2_left - data$smo2_right
    )))
})
