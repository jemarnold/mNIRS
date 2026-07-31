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
