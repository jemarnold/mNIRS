## as_data_list() =====================================================
test_that("as_data_list() errors when dplyr is not installed", {
    skip_if_not_installed("dplyr")

    data <- dplyr::group_by(
        data.frame(time = 1:4, smo2 = c(60, 61, 62, 63), grp = c("a", "a", "b", "b")),
        grp
    )

    ## simulate dplyr being unavailable
    local_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "dplyr")) FALSE else TRUE
        },
        .package = "base"
    )

    expect_error(
        as_data_list(data),
        "is required for grouped data frame input"
    )
})

test_that("as_data_list() wraps a single data frame in a length-1 list", {
    data <- data.frame(time = 1:4, smo2 = c(60, 61, 62, 63))
    result <- as_data_list(data)

    expect_type(result, "list")
    expect_length(result, 1)
    expect_named(result, "interval_1")
    expect_equal(result[[1]], data)
})


## map_mnirs_intervals() ==============================================
test_that("map_mnirs_intervals() returns a list of class 'mnirs'", {
    data <- create_mnirs_data(
        data.frame(time = 1:5, ch1 = rep(50, 5)),
        nirs_channels = "ch1",
        time_channel = "time"
    )
    data_list <- list(a = data, b = data)

    ## each transformer routes list input through `map_mnirs_intervals()`;
    ## the container must carry "mnirs" so `plot()` dispatches to `plot.mnirs`
    results <- list(
        filter = filter_mnirs(
            data_list, method = "moving_average", width = 3, verbose = FALSE
        ),
        resample = resample_mnirs(
            data_list, method = "linear", verbose = FALSE
        ),
        rescale = rescale_mnirs(data_list, range = c(0, 1), verbose = FALSE),
        shift = shift_mnirs(data_list, to = 0, width = 5, verbose = FALSE),
        replace = replace_mnirs(
            data_list, invalid_below = 0, method = "linear", verbose = FALSE
        )
    )

    for (.result in results) {
        expect_type(.result, "list")
        expect_s3_class(.result, "mnirs")
        expect_named(.result, c("a", "b"))
        expect_s3_class(.result$a, "mnirs")
    }
})

test_that("plot() dispatches to plot.mnirs() on transformer list output", {
    skip_if_not_installed("ggplot2")

    data <- create_mnirs_data(
        data.frame(time = 1:5, ch1 = rep(50, 5)),
        nirs_channels = "ch1",
        time_channel = "time"
    )
    result <- filter_mnirs(
        list(a = data, b = data),
        method = "moving_average",
        width = 3,
        verbose = FALSE
    )

    ## an unclassed list falls through to `graphics::plot.default()`
    p <- plot(result)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
})
