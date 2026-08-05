## Test compute_window_bounds() (width/span boundaries) ================
test_that("compute_window_bounds respects width and span boundaries", {
    t <- 1:10

    ## even width: left-biased forward looking; "center" spelling accepted
    result <- compute_window_bounds(t, width = 2, align = "center")
    expect_equal(result$start, 1:10)
    expect_equal(result$end, pmin(1:10 + 1, 10))

    ## span: windows within t ± span/2
    t <- seq(0, 10, by = 0.5)
    result <- compute_window_bounds(t, span = 2)
    expect_equal(result$start[c(1, 10, 21)], c(1L, 8L, 19L))
    expect_equal(result$end[c(1, 10, 21)], c(3L, 12L, 21L))
})

test_that("compute_window_bounds handles degenerate width and span", {
    t <- 1:10

    ## window covering all of t
    expect_equal(compute_window_bounds(t, width = 20), list(start = rep(1L, 10), end = rep(10L, 10)))
    expect_equal(compute_window_bounds(t, span = 20), list(start = rep(1L, 10), end = rep(10L, 10)))

    ## single-sample windows
    expect_equal(compute_window_bounds(t, width = 1), list(start = 1:10, end = 1:10))
    expect_equal(compute_window_bounds(t, span = 0), list(start = 1:10, end = 1:10))
})

## Test compute_valid_neighbours() ======================================
test_that("compute_valid_neighbours works with width", {
    x <- c(1, 2, NA, 4, 5, NA, 7)
    result <- compute_valid_neighbours(x, width = 3)
    expect_length(result, 2) ## one window per NA
    expect_equal(result[[1]], c(2, 4))
    expect_equal(result[[2]], c(5, 7))

    ## wider window splits samples to either side
    x <- c(1, 2, 3, NA, 5, 6, 7)
    result <- compute_valid_neighbours(x, width = 4)
    expect_equal(result[[1]], c(2, 3, 5, 6))

    ## NA at start / end
    expect_equal(compute_valid_neighbours(c(NA, 2, 3, 4), width = 3)[[1]], 2)
    expect_equal(compute_valid_neighbours(c(1, 2, 3, NA), width = 3)[[1]], 3)
})

test_that("compute_valid_neighbours works with span", {
    t <- 0:4
    x <- c(1, 2, NA, 4, 5)
    result <- compute_valid_neighbours(x, t = t, span = 1)
    expect_equal(result[[1]], c(2, 4))

    ## no valid samples within span: falls back to bracketing pair
    result <- compute_valid_neighbours(x, t = t, span = 0.5)
    expect_equal(result, compute_valid_neighbours(x, width = 3))

    ## multiple NAs
    x <- c(1, NA, 3, NA, 5)
    result <- compute_valid_neighbours(x, t = t, span = 1)
    expect_equal(result[[1]], c(1, 3))
    expect_equal(result[[2]], c(3, 5))

    ## NA at start / end
    expect_equal(compute_valid_neighbours(c(NA, 2, 3, 4), span = 0.5)[[1]], 2)
    expect_equal(compute_valid_neighbours(c(1, 2, 3, NA), span = 0.5)[[1]], 3)
})

test_that("compute_valid_neighbours handles degenerate width and span", {
    x <- c(1, 2, NA, 4, 5)

    ## window covering all valid samples
    expect_equal(compute_valid_neighbours(x, width = 8), list(c(1, 2, 4, 5)))
    expect_equal(compute_valid_neighbours(x, span = 8), list(c(1, 2, 4, 5)))

    ## minimal windows return nearest neighbours
    expect_equal(compute_valid_neighbours(x, width = 1)[[1]], c(2, 4))
    expect_equal(compute_valid_neighbours(x, span = 0)[[1]], c(2, 4))
})

## Test compute_local_fun() =============================================
test_that("compute_local_fun calculates rolling medians", {
    x <- c(10, 20, 30, 40, 50)
    bounds <- compute_window_bounds(x, width = 3)
    window_idx <- Map(`:`, bounds$start, bounds$end)
    result <- compute_local_fun(x, window_idx, fn = median)

    expect_equal(result[1], median(x[1:2]))
    expect_equal(result[2], median(x[1:3]))
    expect_equal(result[5], median(x[4:5]))

    ## NA passthrough via ...
    x <- c(1, NA, 3, 4, 5)
    result <- compute_local_fun(x, window_idx, fn = median, na.rm = TRUE)
    expect_equal(result[1:3], c(1, 2, 3.5))
})

test_that("compute_local_fun handles empty and single-sample windows", {
    x <- c(10, 20, 30)
    expect_true(is.na(compute_local_fun(x, list(integer(0)), median, na.rm = TRUE)))
    expect_equal(compute_local_fun(x, list(1, 2, 3), median), x)
})

## Test compute_outliers() ==============================================
test_that("compute_outliers flags outliers with local medians", {
    x <- c(1, 2, 3, 100, 5)
    t <- 1:5
    result <- compute_outliers(x, t, outlier_cutoff = 3, width = 3)

    expect_type(result$local_medians, "double")
    expect_length(result$local_medians, length(x))
    expect_equal(result$is_outlier, c(FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("compute_outliers threshold sensitivity via outlier_cutoff", {
    x <- c(1, 2, 3, 10, 5)
    t <- 1:5
    strict <- compute_outliers(x, t, outlier_cutoff = 1, width = 3)
    lenient <- compute_outliers(x, t, outlier_cutoff = 10, width = 3)

    expect_gt(sum(strict$is_outlier), sum(lenient$is_outlier))
    expect_equal(strict$local_medians, lenient$local_medians)
})

test_that("compute_outliers handles clean data and NA", {
    t <- 1:5
    result <- compute_outliers(1:5, t, outlier_cutoff = 3, width = 3)
    expect_false(any(result$is_outlier))

    ## this shouldn't happen with handle_na, but just in case
    result <- compute_outliers(c(1, 2, NA, 100, 5), t, outlier_cutoff = 3, width = 3)
    expect_false(any(result$is_outlier))
})

## Test compute_window_bounds() =========================================
test_that("compute_window_bounds width bounds clamp at edges", {
    t <- 1:10
    result <- compute_window_bounds(t, width = 3)
    expect_equal(result$start, pmax(1:10 - 1, 1))
    expect_equal(result$end, pmin(1:10 + 1, 10))

    ## left = forward looking; right = backward looking
    result <- compute_window_bounds(t, width = 3, align = "left")
    expect_equal(result$start, 1:10)
    expect_equal(result$end, pmin(1:10 + 2, 10))

    result <- compute_window_bounds(t, width = 3, align = "right")
    expect_equal(result$start, pmax(1:10 - 2, 1))
    expect_equal(result$end, 1:10)
})

test_that("compute_window_bounds span windows stay within t ± span/2", {
    t <- seq(0, 10, by = 0.5)
    result <- compute_window_bounds(t, span = 2)
    expect_true(all(t[result$start] >= t - 1 & t[result$end] <= t + 1))
})

## Test window_sums() ===================================================
test_that("window_sums matches per-window sums", {
    set.seed(7)
    x <- rnorm(50) + 100 ## offset data exercises differencing error
    bounds <- compute_window_bounds(seq_along(x), width = 5)
    reference <- vapply(seq_along(x), \(.i) {
        sum(x[bounds$start[.i]:bounds$end[.i]])
    }, numeric(1))
    expect_equal(window_sums(x, bounds), reference)

    ## logical input sums to valid counts
    expect_equal(
        window_sums(is.finite(c(1, NA, 3)), compute_window_bounds(1:3, width = 3)),
        c(1, 2, 1)
    )
})

## Test compute_local_mean() ============================================
test_that("compute_local_mean matches naive rolling mean", {
    x <- c(10, 20, 30, 40, 50)
    bounds <- compute_window_bounds(seq_along(x), width = 3)
    reference <- vapply(seq_along(x), \(.i) {
        mean(x[bounds$start[.i]:bounds$end[.i]])
    }, numeric(1))
    expect_equal(compute_local_mean(x, bounds), reference)
})

test_that("compute_local_mean handles NA and Inf", {
    x <- c(1, 2, NA, 4, 5)
    bounds <- compute_window_bounds(seq_along(x), width = 3)

    ## NA propagates unless dropped
    result <- compute_local_mean(x, bounds, na.rm = FALSE)
    expect_equal(is.na(result), c(FALSE, TRUE, TRUE, TRUE, FALSE))
    expect_equal(compute_local_mean(x, bounds, na.rm = TRUE), c(1.5, 1.5, 3, 4.5, 4.5))

    ## all-NA window returns NA, not NaN
    result <- compute_local_mean(c(NA, NA, NA, 4, 5), bounds, na.rm = TRUE)
    expect_true(is.na(result[1]) && !is.nan(result[1]))

    ## Inf treated as missing: later windows unaffected
    result <- compute_local_mean(c(1, Inf, 3, 4, 5), bounds, na.rm = TRUE)
    expect_equal(result[4:5], c(4, 4.5))
    expect_true(all(is.finite(result)))

    ## no finite values at all: NA for every window, regardless of na.rm
    expect_identical(
        compute_local_mean(rep(NA_real_, 5), bounds, na.rm = TRUE),
        rep(NA_real_, 5)
    )
    expect_identical(
        compute_local_mean(c(Inf, -Inf, NaN, NA, Inf), bounds),
        rep(NA_real_, 5)
    )

    ## guard returns one NA per window, not per element of x
    short_bounds <- compute_window_bounds(seq_along(x), idx = 1:2, width = 3)
    expect_identical(
        compute_local_mean(rep(NA_real_, 5), short_bounds),
        rep(NA_real_, 2)
    )
})

test_that("compute_local_mean min_obs excludes partial windows", {
    x <- as.numeric(1:10)
    bounds <- compute_window_bounds(seq_along(x), width = 5)
    result <- compute_local_mean(x, bounds, min_obs = 5L)

    ## edge windows span fewer than 5 samples
    expect_true(all(is.na(result[c(1:2, 9:10)])))
    expect_equal(result[3:8], as.numeric(3:8))
})

test_that("compute_local_mean floating point error stays below resolution", {
    ## cumsum differencing loses relative precision over long vectors;
    ## absolute error must stay far below NIRS measurement resolution
    set.seed(42)
    n <- 1e5
    bounds <- compute_window_bounds(seq_len(n), width = 15)
    check_idx <- c(1:100, seq(1000, n, by = 1000), (n - 99):n)
    naive <- function(x) {
        vapply(check_idx, \(.i) {
            mean(x[bounds$start[.i]:bounds$end[.i]])
        }, numeric(1))
    }

    ## realistic NIRS magnitudes
    x <- runif(n, 60, 70) + rnorm(n, sd = 0.1)
    expect_lt(max(abs(compute_local_mean(x, bounds)[check_idx] - naive(x))), 1e-8)

    ## large offset with small signal: worst case for cancellation error
    x <- 1e4 + rnorm(n, sd = 0.01)
    expect_lt(max(abs(compute_local_mean(x, bounds)[check_idx] - naive(x))), 1e-6)
})

## Test window_min_obs() ================================================
test_that("window_min_obs returns width or converted span", {
    t <- seq(0, 99, by = 0.5) ## sample rate = 2

    expect_equal(window_min_obs(width = 5, span = NULL, t), 5)
    ## min_n floor applies
    expect_equal(window_min_obs(width = 1, span = NULL, t, min_n = 2L), 2)

    ## span converted via sample rate, less two-sample buffer
    expect_equal(window_min_obs(width = NULL, span = 5, t), 8)
    expect_equal(window_min_obs(width = NULL, span = 1, t, min_n = 2L), 2)
})

## Test median_no_na() ==================================================
test_that("median_no_na matches median with na.rm", {
    expect_identical(median_no_na(c(3, 1, 2)), median(c(3, 1, 2)))
    ## even-n mean-of-pair identical to median.default on doubles
    expect_identical(median_no_na(c(0.1, 0.2, 0.3, 0.7)), median(c(0.1, 0.2, 0.3, 0.7)))
    expect_identical(
        median_no_na(c(1, NA, 3, 2)),
        median(c(1, NA, 3, 2), na.rm = TRUE)
    )
    expect_identical(median_no_na(numeric(0)), NA_real_)
    expect_identical(median_no_na(c(NA_real_, NA_real_)), NA_real_)
})

## Test compute_col_medians() ============================================
test_that("compute_col_medians matches per-column median", {
    set.seed(1)
    m <- matrix(runif(35, 60, 70), nrow = 5)
    ## NA padding with mixed even/odd valid counts per column
    m[1, 2] <- NA
    m[1:2, 3] <- NA
    m[1:4, 4] <- NA
    m[, 5] <- NA

    result <- compute_col_medians(m)
    expect_equal(result, unname(apply(m, 2, median, na.rm = TRUE)))
    expect_true(is.na(result[5])) ## all-NA column
})
