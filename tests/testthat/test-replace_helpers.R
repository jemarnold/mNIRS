## Test compute_local_windows() =========================================
test_that("compute_local_windows returns correct structure", {
    t <- 1:10
    width <- 3
    result <- compute_local_windows(t, width = width, span = NULL)

    expect_type(result, "list")
    expect_length(result, length(t))
    expect_true(all(vapply(result, is.numeric, logical(1))))

    ## with span
    t <- seq(0, 10, by = 0.5)
    span <- 1
    result <- compute_local_windows(t, width = NULL, span = span)

    expect_type(result, "list")
    expect_length(result, length(t))
})

test_that("compute_local_windows respects width boundaries", {
    t <- 1:10
    width <- 3
    result <- compute_local_windows(t, width = width)

    # First element: can only look forward
    expect_equal(result[[1]], 1:2)
    # Last element: can only look backward
    expect_equal(result[[10]], 9:10)
    # Middle element: symmetrical window
    expect_equal(result[[5]], c(4:6))

    ## even width align centreed
    width <- 2
    result <- compute_local_windows(t, width = width, align = "centre")

    # left-biased forward looking
    expect_equal(result[[1]], 1:2)
    # Last element only at idx
    expect_equal(result[[10]], 10)
    # left-biased forward looking
    expect_equal(result[[5]], c(5:6))

    ## with span
    t <- seq(0, 10, by = 0.5)
    span <- 2
    result <- compute_local_windows(t, span = span)

    # First element: can only look forward
    expect_equal(result[[1]], 1:3)
    expect_equal(t[result[[1]]], t[1:3])
    # Last element: can only look backward
    expect_equal(result[[21]], 19:21)
    expect_equal(t[result[[21]]], t[19:21])
    # Middle element: symmetrical window
    expect_equal(result[[10]], c(8:12))
    expect_equal(t[result[[10]]], t[c(8:12)])
})

test_that("compute_local_windows align = 'center'", {
    t <- 1:10
    width <- 2
    result <- compute_local_windows(t, width = width, align = "center")

    expect_equal(result[[1]], 1:2)
    expect_equal(result[[10]], 10)
    expect_equal(result[[5]], c(5:6))
})

test_that("compute_local_windows width and span >= length(x) works", {
    t <- 1:10
    result <- compute_local_windows(t, width = 20)
    expect_equal(result, rep(list(t), length(t)))

    result <- compute_local_windows(t, span = 20)
    expect_equal(result, rep(list(t), length(t)))
})

test_that("compute_local_windows works with width = 1 and span = 0", {
    t <- 1:10
    result <- compute_local_windows(t, width = 1)
    expect_equal(result, as.list(t))

    result <- compute_local_windows(t, span = 0)
    expect_equal(result, as.list(t))

    compute_local_windows(t, width = 1, span = 0)
})

## test compute_valid_neighbours() ===========================
test_that("compute_valid_neighbours() works with width", {
    x <- c(1, 2, NA, 4, 5, NA, 7)

    # width = 1
    result <- compute_valid_neighbours(x, width = 3)
    expect_length(result, 2) # Two NAs
    expect_equal(result[[1]], c(2, 4))
    expect_equal(result[[2]], c(5, 7))

    # width = 3
    x <- c(1, 2, 3, NA, 5, 6, 7)
    result <- compute_valid_neighbours(x, width = 4)
    expect_equal(result[[1]], c(2, 3, 5, 6))

    # Edge case: NA at start
    x <- c(NA, 20, 30, 40)
    result <- compute_valid_neighbours(x, width = 3)
    expect_equal(result[[1]], 2)

    # Edge case: NA at end
    x <- c(10, 20, 30, NA)
    result <- compute_valid_neighbours(x, width = 3)
    expect_equal(result[[1]], 3)
})

test_that("compute_valid_neighbours() works with span", {
    t <- 0:4
    x <- c(1, 2, NA, 4, 5)

    # span = 1 (includes t[2] = 1 and t[4] = 3)
    result <- compute_valid_neighbours(x, t = t, span = 1)
    expect_length(result, 1)
    expect_equal(result[[1]], c(2, 4))

    # span = 0.5 (no values of t within t ± span should return same as width = 1)
    result <- compute_valid_neighbours(x, t = t, span = 0.5)
    expect_equal(result[[1]], c(2, 4))
    expect_equal(result, compute_valid_neighbours(x, width = 3))

    # Multiple NAs
    x <- c(1, NA, 3, NA, 5)
    result <- compute_valid_neighbours(x, t = t, span = 1)
    expect_length(result, 2)
    expect_equal(result[[1]], c(1, 3))
    expect_equal(result[[2]], c(3, 5))

    # Edge case: NA at start
    x <- c(NA, 2, 3, 4)
    expect_equal(compute_valid_neighbours(x, span = 1)[[1]], 2)
    expect_equal(compute_valid_neighbours(x, span = 0.5)[[1]], 2)

    # Edge case: NA at end
    x <- c(1, 2, 3, NA)
    expect_equal(compute_valid_neighbours(x, span = 1)[[1]], 3)
    expect_equal(compute_valid_neighbours(x, span = 0.5)[[1]], 3)
})

test_that("compute_valid_neighbours width and span >= length(x) works", {
    x <- c(1, 2, NA, 4, 5)
    result <- compute_valid_neighbours(x, width = 8)
    expect_equal(result, list(x[!is.na(x)]))

    result <- compute_valid_neighbours(x, span = 8)
    expect_equal(result, list(x[!is.na(x)]))
})

test_that("compute_valid_neighbours works with width = 1 and span = 0", {
    x <- c(1, 2, NA, 4, 5)
    result <- compute_valid_neighbours(x, width = 1)
    expect_equal(result[[1]], c(2, 4))

    result <- compute_valid_neighbours(x, span = 0)
    expect_equal(result[[1]], c(2, 4))
})

## Test compute_local_fun() =========================================
test_that("compute_local_fun returns correct length", {
    x <- c(1, 2, 3, 4, 5)
    window_idx <- compute_local_windows(x, width = 1)
    result <- compute_local_fun(x, window_idx, fn = median)

    expect_type(result, "double")
    expect_length(result, length(x))
})

test_that("compute_local_fun calculates correct medians", {
    x <- c(10, 20, 30, 40, 50)
    # Manual windows: each looks at neighbours
    window_idx <- compute_local_windows(x, width = 3)
    result <- compute_local_fun(x, window_idx, fn = median)

    expect_equal(result[1], median(x[1:2])) # median of x[2]
    expect_equal(result[2], median(x[1:3])) # median of x[c(1,3)]
    expect_equal(result[5], median(x[4:5])) # median of x[4]
})

test_that("compute_local_fun handles NA values", {
    ## this shouldn't happen with handle_na, but just in case
    x <- c(1, NA, 3, 4, 5)
    window_idx <- compute_local_windows(x, width = 3)
    result <- compute_local_fun(x, window_idx, fn = median, na.rm = TRUE)

    expect_false(all(is.na(result[1:3])))
    expect_equal(result[1], median(x[1]))
    expect_equal(result[2], median(x[c(1, 3)]))
    expect_equal(result[3], median(x[c(3, 4)]))
})

test_that("compute_local_fun() handles empty windows", {
    x <- c(1, 2, 3)
    window_idx <- list(integer(0))
    result <- compute_local_fun(x, window_idx, fn = median, na.rm = TRUE)
    expect_true(is.na(result))
})

test_that("compute_local_fun() handles single values", {
    x <- c(10, 20, 30)
    window_idx <- list(1, 2, 3)
    result <- compute_local_fun(x, window_idx, fn = median)
    expect_equal(result, x)
})


## Test compute_outliers() ==============================================
test_that("compute_outliers returns list of logical and numeric vectors", {
    x <- c(1, 2, 3, 100, 5)
    t <- 1:5
    result <- compute_outliers(x, t, outlier_cutoff = 3, width = 3)

    expect_type(result$local_medians, "double")
    expect_type(result$is_outlier, "logical")
    expect_length(result$local_medians, length(x))
    expect_length(result$is_outlier, length(x))
    expect_true(result$is_outlier[4L]) # 100 should be flagged
    expect_false(any(result$is_outlier[-4L])) # Normal values should not be flagged
})

test_that("compute_outliers threshold sensitivity via outlier_cutoff", {
    x <- c(1, 2, 3, 10, 5)
    t <- 1:5

    # Strict threshold
    strict <- compute_outliers(x, t, outlier_cutoff = 1, width = 3)
    # Lenient threshold
    lenient <- compute_outliers(x, t, outlier_cutoff = 10, width = 3)

    expect_true(sum(strict$is_outlier) > sum(lenient$is_outlier))
    expect_equal(strict$local_medians, lenient$local_medians)
})

test_that("compute_outliers handles no outliers", {
    x <- 1:5
    t <- 1:5
    result <- compute_outliers(x, t, outlier_cutoff = 3, width = 3)

    expect_true(all(!result$is_outlier))
})

test_that("compute_outliers handles NA", {
    ## this shouldn't happen with handle_na, but just in case
    x <- c(1, 2, NA, 100, 5)
    t <- 1:5
    result <- compute_outliers(x, t, outlier_cutoff = 3, width = 3)

    expect_true(all(!result$is_outlier))
})
