## helper to create test data with metadata
create_test_data <- function(
    time_max = 10,
    sample_rate = 10
) {
    time <- seq(0, time_max, by = 1 / sample_rate)
    nrow <- length(time)
    data <- tibble(
        time = time,
        nirs1 = rnorm(nrow, 50, 5),
        nirs2 = rnorm(nrow, 60, 5),
        event = c(1, rep(NA, nrow - 2), 2),
    )
    class(data) <- c("mnirs", class(data))
    attr(data, "time_channel") <- "time"
    attr(data, "nirs_channels") <- c("nirs1", "nirs2")
    attr(data, "event_channel") <- "event"
    attr(data, "sample_rate") <- sample_rate

    return(data)
}


## conditions attribute to the top-level user-facing function ============
test_that("read_mnirs errors attribute to read_mnirs", {
    ## raised deep in read_file()
    err <- expect_error(
        read_mnirs("nonexistent_file.csv", verbose = FALSE),
        "File not found"
    )
    expect_equal(rlang::call_name(conditionCall(err)), "read_mnirs")
})

test_that("filter_mnirs method errors attribute to the generic", {
    data <- create_test_data()

    ## method -> filter_moving_average() -> validate_width_span() chain
    err <- expect_error(
        filter_mnirs(data, method = "moving_average", verbose = FALSE),
        "Window size undefined"
    )
    expect_equal(rlang::call_name(conditionCall(err)), "filter_mnirs")

    ## resolve_channel_args() -> arg_match0(error_call = <call>)
    err <- expect_error(
        filter_mnirs(
            data, method = "butterworth", W = 0.02,
            type = "bogus", verbose = FALSE
        ),
        "type"
    )
    expect_equal(rlang::call_name(conditionCall(err)), "filter_mnirs")

    ## per-channel `method` passes, leaving the error for `width`
    err <- expect_error(
        filter_mnirs(
            data,
            method = list(nirs1 = "moving_average", nirs2 = "smooth_spline"),
            width = list(nirs1 = -1),
            spar = list(nirs2 = 0.5),
            verbose = FALSE
        ),
        "width"
    )
    expect_equal(rlang::call_name(conditionCall(err)), "filter_mnirs")
})

test_that("filter_moving_average called directly attributes to itself", {
    err <- expect_error(filter_moving_average(1:10), "Window size undefined")
    expect_equal(rlang::call_name(conditionCall(err)), "filter_moving_average")

    err <- expect_error(filter_moving_average(1:10), "Window size undefined")
    expect_equal(
        rlang::call_name(conditionCall(err)), "filter_moving_average"
    )
})

test_that("extract_intervals errors attribute to extract_intervals", {
    data <- create_test_data()

    ## raised in apply_span() for out-of-bounds intervals
    err <- expect_error(
        extract_intervals(
            data, start = by_time(1e6), span = 10, verbose = FALSE
        ),
        "outside data bounds"
    )
    expect_equal(rlang::call_name(conditionCall(err)), "extract_intervals")
})

test_that("warnings attribute to the top-level user-facing function", {
    skip_if_not_installed("signal")
    data <- create_test_data()

    ## validate_sample_rate() warns when sample_rate disagrees with estimate
    wrn <- expect_warning(
        filter_mnirs(
            data, method = "butterworth", W = 0.02, sample_rate = 999
        ),
        "inconsistent"
    )
    expect_equal(rlang::call_name(conditionCall(wrn)), "filter_mnirs")
})
