## slope ===========================================================
test_that("slope returns correct structure", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    result <- slope(x)

    expect_type(result, "double")
    expect_equal(length(result), 1)
    expect_gt(result, 0)
})

test_that("slope handles negative values correctly", {
    x <- c(-5, -3, -1, 1, 3, 5)
    expect_equal(slope(x), 2)

    ## negative slope
    x_desc <- c(10, 8, 6, 4, 2)
    expect_equal(slope(x_desc), -2)

    ## mixed with irregular t
    x <- c(-2, 0, 2, 4)
    t <- c(1, 2.1, 3.3, 4.4)
    expect_lt(slope(x, t), 2)
})

test_that("slope calculates edge cases correctly", {
    ## identity = 1
    expect_equal(slope(1:5, 1:5), 1)
    ## single value = NA
    expect_true(is.na(slope(5)))
    ## all identical values = 0
    expect_true(slope(rep(5, 10)) == 0)
    ## all NA_real_ = NA
    expect_true(is.na(slope(x = rep(NA_real_, 5))))
    ## all t values identical = NA
    expect_true(is.na(slope(1:5, rep(1, 5))))
    ## character vectors error
    expect_error(slope(x = letters[1:4]), "valid.*numeric")
    expect_error(slope(1:4, t = letters[1:4]), "valid.*numeric")
    ## unequal length = error
    expect_error(slope(x = 1:5, t = 1:3), "equal length")
    ## all invalid = NA
    expect_error(slope(list()), "valid.*numeric")
    expect_true(is.na(slope(numeric(0))))
    expect_error(slope(NULL), "valid.*numeric")
    expect_true(is.na(slope(rep(NA_real_, 4))))
    expect_true(is.na(slope(rep(NaN, 4))))
    expect_true(is.na(slope(rep(Inf, 4))))

    ## NaN & Inf propagated by default
    x <- c(1, NaN, 3, Inf, 5)
    expect_equal(slope(x), NA_real_)
    ## NaN & Inf ignored with `na.rm = TRUE`
    expect_equal(slope(x, na.rm = TRUE), 1)
})

test_that("slope handles NA values correctly", {
    x_na <- c(1, NA, 3, 4, 5)

    # NA propagated by default, ignored with `na.rm = TRUE`
    expect_equal(slope(x_na), NA_real_)
    expect_equal(slope(x_na, na.rm = TRUE), 1)

    # NA in t
    t_na <- c(1, 2, NA, 4, 5)
    expect_equal(
        slope(1:5, t_na, na.rm = TRUE),
        coef(lm(1:5 ~ t_na))[[2L]]
    )
})

test_that("slope returns same as lm model", {
    ## irregular time
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t_irreg <- sort(runif(length(x), 0, 10))
    expect_equal(slope(x, t_irreg), unname(coef(lm(x ~ t_irreg))[2L]))

    x <- c(1, 3, NA, 5, 8, 7, 9, 12, NA, NA, NA, 17, 18)
    expect_equal(slope(x, na.rm = TRUE), coef(lm(x ~ seq_along(x)))[[2L]])

    ## boundary NAs
    x <- c(NA, 3, NA, 5, 8, 7, 9, 12, NA, NA, NA, 17, NA)
    expect_equal(slope(x, na.rm = TRUE), coef(lm(x ~ seq_along(x)))[[2L]])
})

test_that("slope works with min_obs", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    expect_gt(slope(x, partial = TRUE), 1)
    expect_gt(slope(x, min_obs = length(x)), 1)
    ## min_obs > length(x) = NA
    expect_true(is.na(slope(x, min_obs = 20)))
    ## n < min_obs = NA
    expect_true(is.na(slope(x[1:2], min_obs = 3, partial = TRUE)))
    ## min_obs < 2 ## silent argument, no warning
    # expect_warning(
    #     expect_gt(slope(x, min_obs = 1, verbose = TRUE), 1),
    #     "set to .*2"
    # )
    # expect_silent(
    #     expect_gt(slope(x, min_obs = 1, verbose = FALSE), 1)
    # )
    # expect_warning(
    #     expect_gt(slope(x, min_obs = 1, bypass_checks = FALSE), 1),
    #     "min_obs.*2"
    # )
    ## bypass doesn't warn about min_obs change
    # expect_silent(
    #     expect_gt(slope(x, min_obs = 1, bypass_checks = TRUE), 1)
    # )
})

test_that("slope handles large integer sequences without overflow", {
    n <- 1e6
    x <- rep(1, n)
    expect_equal(slope(x), 0)

    x <- as.double(seq_len(n))
    expect_equal(slope(x), 1)
})

test_that("slope respects bypass_checks", {
    x <- c(1, 3, 2, 5, 8)
    expect_error(
        slope(x, t = 1:4, bypass_checks = FALSE),
        "numeric.*equal length"
    )
    ## bypass doesn't validate length, returns bad result and warning
    expect_warning(slope(x, t = 1:4, bypass_checks = TRUE), "not a multiple")
})

test_that("slope returns intercept", {
    x <- c(1, 3, 2, 5, 8)
    expect_equal(attr(slope(x, intercept = TRUE), "intercept"), -1)
})


## rolling_slope ===========================================================
test_that("rolling_slope returns correct structure", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    result <- rolling_slope(x, width = 3)

    expect_type(result, "double")
    expect_equal(length(result), length(x))

    ## attribute for window_idx
    expect_equal(
        attr(rolling_slope(x, width = 3, window_idx = TRUE), "window_idx"),
        compute_local_windows(t = seq_along(x), width = 3)
    )
})

test_that("rolling_slope handles negative values correctly", {
    x <- c(-5, -3, -1, 1, 3, 5)
    expect_all_equal(rolling_slope(x, width = 3, partial = TRUE), 2)
    expect_all_equal(rolling_slope(x, span = 3, partial = TRUE), 2)

    ## negative slope
    x_desc <- c(10, 8, 6, 4, 2)
    expect_all_equal(rolling_slope(x_desc, width = 3, partial = TRUE), -2)
    expect_all_equal(rolling_slope(x_desc, span = 2, partial = TRUE), -2)
})

test_that("rolling_slope works with width", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    expect_message(rolling_slope(x, width = 3, span = 3), "width.*span")

    ## width = 3, centre aligned
    result <- rolling_slope(x, width = 3, align = "centre")
    expect_equal(length(result), length(x))
    expect_all_true(is.na(result[c(1, 10)])) ## boundary < width
    expect_all_true(is.na(result[-c(1, 10)]) >= 0)
    expect_equal(result[5], slope(x[4:6])) # centre window

    ## width = 3, left aligned
    result_left <- rolling_slope(x, width = 3, align = "left")
    expect_equal(length(result_left), length(x))
    expect_all_true(is.na(result_left[9:10])) ## boundary < width
    expect_all_true(is.na(result_left[-c(9:10)]) >= 0)
    expect_equal(result_left[5], slope(x[5:7]))
    
    ## width = 3, right aligned
    result_right <- rolling_slope(x, width = 3, align = "right")
    expect_equal(length(result_right), length(x))
    expect_all_true(is.na(result_right[1:2])) ## boundary < width
    expect_all_true(is.na(result_right[-c(1:2)]) >= 0)
    expect_equal(result_right[5], slope(x[3:5]))
})

test_that("rolling_slope works with span", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    ## TODO span < 7 centred effectively returns partial = TRUE because min_obs-2 and centred
    ## span = 7, centre aligned
    result <- rolling_slope(x, span = 7, align = "centre")
    expect_equal(length(result), length(x))
    expect_all_true(is.na(result[c(1, 10)])) ## boundary < width
    expect_all_true(is.na(result[-c(1, 10)]) >= 0)
    expect_equal(result[5L], slope(x[2:8])) # centre window

    ## span = 7, left aligned
    result_left <- rolling_slope(x, span = 7, align = "left")
    expect_equal(length(result_left), length(x))
    expect_all_true(is.na(result_left[7:10])) ## boundary < width
    expect_all_true(is.na(result_left[-c(9:10)]) >= 0)
    expect_equal(result_left[3L], slope(x[3:10]))

    ## span = 7, right aligned
    result_right <- rolling_slope(x, span = 7, align = "right")
    expect_equal(length(result_right), length(x))
    expect_all_true(is.na(result_right[1:4])) ## boundary < width
    expect_all_true(is.na(result_right[-c(1:2)]) >= 0)
    expect_equal(result_right[7L], slope(x[1:7]))
})

test_that("rolling_slope handles invalid width/span", {
    expect_all_true(is.na(rolling_slope(1:5, width = 1))) |>
        expect_warning("Insufficient valid")
    expect_all_true(is.na(rolling_slope(1:5, span = 0))) |>
        expect_warning("Insufficient valid")

    expect_error(rolling_slope(1:5), "width.*span.*must be defined")
    expect_error(rolling_slope(1:5, span = NA), "valid.*numeric")
    expect_error(rolling_slope(1:5, width = NA), "valid.*integer")
    expect_error(rolling_slope(1:5, span = -1), "valid.*numeric")
    expect_error(rolling_slope(1:5, width = -1), "valid.*integer")

    ## invalid align
    expect_error(
        rolling_slope(1:5, width = 3, align = "invalid"),
        "arg.*should be one of"
    )
})

test_that("rolling_slope works with partial", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    expect_all_true(rolling_slope(x, width = 3, partial = TRUE) >= 0)
    expect_all_true(rolling_slope(x, span = 3, partial = TRUE) >= 0)
    expect_all_true(is.na(
        rolling_slope(x, width = 3, partial = FALSE)[c(1, 10)]
    ))
    expect_all_true(is.na(
        rolling_slope(x, span = 7, partial = FALSE)[c(1, 10)]
    ))

    ## n < 2 = NA regardless of partial
    expect_true(is.na(rolling_slope(x[1L], width = 3, partial = TRUE))) |>
        expect_warning("Insufficient valid")
})

test_that("rolling_slope works with align", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    result_right <- rolling_slope(x, width = 3, align = "right")
    result_centre <- rolling_slope(x, width = 3, align = "centre")
    result_left <- rolling_slope(x, width = 3, align = "left")

    ## non-NA values are identical
    expect_true(all.equal(
        result_right[!is.na(result_right)],
        result_centre[!is.na(result_centre)]
    ))
    expect_true(all.equal(
        result_right[!is.na(result_right)],
        result_left[!is.na(result_left)]
    ))

    ## "center" spelling alias accepted without error
    expect_equal(
        rolling_slope(x, width = 3, align = "center"),
        result_centre
    )
})

test_that("rolling_slope handles width larger than data", {
    x <- c(1, 3, 2, 5, 8)

    ## width > n
    expect_warning(
        result <- rolling_slope(x, width = 10, partial = FALSE),
        "Insufficient valid"
    )
    expect_all_true(is.na(result))

    ## width > n but partial = TRUE
    result <- rolling_slope(x, width = 10, partial = TRUE)
    expect_equal(length(result), length(x))
    expect_equal(var(result), 0) ## all same slope for align = "centre"
})

test_that("rolling_slope calculates identity", {
    expect_all_equal(rolling_slope(1:5, 1:5, width = 3)[-c(1, 5)], 1)
    expect_all_equal(rolling_slope(1:5, 1:5, span = 3), 1)
})

test_that("rolling_slope calculates all identical values", {
    expect_all_equal(rolling_slope(rep(5, 5), width = 3)[-c(1, 5)], 0)
    expect_all_equal(rolling_slope(rep(5, 5), span = 3), 0)

    ## all t values identical = NA
    expect_all_true(is.na(rolling_slope(1:5, rep(1, 5), width = 3))) |>
        expect_warning("Insufficient valid")
    expect_all_true(is.na(rolling_slope(1:7, rep(1, 7), span = 3))) |>
        expect_warning("Insufficient valid")
})

test_that("rolling_slope calculates invalid", {
    ## single value
    expect_true(is.na(rolling_slope(x = 5, width = 3))) |>
        expect_warning("Insufficient valid")
    ## non numeric NA errors
    expect_error(rolling_slope(x = rep(NA, 5), width = 3), "valid.*numeric")
    ## NA_real_ = NA_real_
    expect_all_true(is.na(rolling_slope(x = rep(NA_real_, 5), width = 3)))
    ## empty input returns empty output for input == output length consistency
    expect_equal(length(rolling_slope(x = numeric(0), width = 3)), 0) |>
        expect_warning("Insufficient valid")
    ## invalid type = error
    expect_error(rolling_slope(list(), width = 3), "valid.*numeric")
    expect_error(rolling_slope(NULL, width = 3), "valid.*numeric")
    expect_equal(rolling_slope(c(1, 2, NULL, 3), width = 3), c(NA, 1, NA))
    ## NA_real_, NaN, Inf = NA_real_
    expect_all_true(is.na(rolling_slope(rep(NaN, 4), width = 3)))
    expect_all_true(is.na(rolling_slope(rep(Inf, 4), width = 3)))
    ## charater = error
    expect_error(rolling_slope(x = letters[1:4]), "valid.*numeric")
    expect_error(rolling_slope(1:4, t = letters[1:4]), "valid.*numeric")
    ## unequal length = error
    expect_error(rolling_slope(1:5, 1:3), "equal length")
})

test_that("rolling_slope handles NaN & Inf with na.rm", {
    x <- c(1, 2, NaN, 4, 5, Inf, 7, 8)
    ## reconcile different results between methods?
    ## differences are that centred `span` is divided by two
    ## so 2/1 = ±1 sample on each side
    ## whereas `width` is left-biased forward looking c(idx, idx+1)
    ## answers are different because questions are different
    expect_equal(
        rolling_slope(x, width = 3, na.rm = FALSE),
        c(rep(NA, 4), rep(1, 3), NA)
    )
    expect_equal(
        rolling_slope(x, span = 3, na.rm = FALSE),
        c(1, rep(NA, 3), rep(1, 4))
    )

    expect_equal(
        rolling_slope(x, width = 3, na.rm = TRUE),
        c(NA, rep(1, 6), NA)
    )
    expect_all_equal(
        rolling_slope(x, span = 3, na.rm = TRUE),
        1
    )
})

test_that("rolling_slope span returns same as lm model", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
    span <- 3

    result <- rolling_slope(x, t, span = span, partial = TRUE)
    lm_result <- vapply(t, \(.x) {
        start_t <- max(t[1], .x - span / 2)
        end_t <- min(t[length(x)], .x + span / 2)
        window_idx <- which(t >= start_t & t <= end_t)
        # ## exception from lm() NA handling for all NA
        if (sum(!is.na(x[window_idx])) > 1) {
            coef(lm(x[window_idx] ~ t[window_idx]))[[2L]]
        } else {NA}
    }, numeric(1))

    expect_equal(result, lm_result)

    ## width roll_lm with NA edges
    width <- 3
    result <- rolling_slope(x, width = width)
    lm_result <- vapply(t, \(.x) {
        start_idx <- .x - floor(width / 2)
        end_idx <- .x + floor(width / 2)
        window_idx <- start_idx:end_idx
        ## exception from lm() NA handling for all NA
        if (all(window_idx >= 1 & window_idx <= length(x))) {
            coef(lm(x[window_idx] ~ t[window_idx]))[[2L]]
        } else {
            NA
        }        
    }, numeric(1))

    expect_equal(result, lm_result)

    ## visual confirmation
    # library(ggplot2)
    # tibble(x = t, y = .env$x) |>
    #     ggplot(aes(x, y)) +
    #     scale_x_continuous(breaks = scales::breaks_pretty(n=10)) +
    #     geom_line() +
    #     geom_point(aes(y = result, colour = "rolling_slope"),
    #                shape = 21, size = 5, stroke = 1) +
    #     geom_point(aes(y = lm_result, colour = "lm"),
    #                shape = 21, size = 3, stroke = 1)
})

test_that("rolling_slope returns same as zoo::rollapply()", {
    skip_if_not_installed("zoo")
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    coef.fn <- \(x) coef(lm(x ~ seq_along(x)))[[2L]]

    ## align = "centre"
    rolling_slope_centre <- rolling_slope(
        x, width = 3, partial = TRUE, align = "centre"
    )
    rollapply_centre <- zoo::rollapply(
        x,
        FUN = coef.fn,
        width = 3,
        align = "center",
        partial = TRUE
    )
    expect_equal(rollapply_centre, rolling_slope_centre)
    rolling_slope_centre <- rolling_slope(
        x, width = 3, span = 0, partial = TRUE, align = "centre", verbose = FALSE
    )
    expect_equal(rollapply_centre, rolling_slope_centre)
    
    ## align = "left"
    rolling_slope_left <- rolling_slope(
        x, width = 3, partial = TRUE, align = "left"
    )
    rollapply_left <- zoo::rollapply(
        x,
        FUN = coef.fn,
        width = 3,
        align = "left",
        partial = TRUE
    )
    expect_equal(rollapply_left, rolling_slope_left)
    
    ## align = "right"
    rolling_slope_right <- rolling_slope(
        x, width = 3, partial = TRUE, align = "right"
    )
    rollapply_right <- zoo::rollapply(
        x,
        FUN = coef.fn,
        width = 3,
        align = "right",
        partial = TRUE
    )
    expect_equal(rollapply_right, rolling_slope_right)
})

test_that("rolling_slope handles span in units of t", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    t <- seq_along(x)
    t_5hz <- t / 5

    expect_equal(
        rolling_slope(x, t, span = 5, partial = TRUE) * 5,
        rolling_slope(x, t_5hz, span = 1, partial = TRUE)
    )

    t <- c(0.5, 1.5, 3, 3.5, 4, 5.5, 7, 8, 10, 11)
    n <- length(x)
    span <- 3
    result <- rolling_slope(x, t, span = span, partial = TRUE)

    lm_result <- vapply(seq_len(n), \(.x) {
        current_t <- t[.x]
        start_t <- max(t[1L], current_t - span / 2)
        end_t <- min(t[n], current_t + span / 2)
        window_idx <- which(t >= start_t & t <= end_t)
        
        coef(lm(x[window_idx] ~ t[window_idx]))[[2L]]
    }, numeric(1))

    expect_equal(result, lm_result)

    x <- c(1, 3, 2, 5, 8)
    ## span & min_obs > n = NA
    result <- rolling_slope(x, span = 10, verbose = FALSE)
    expect_length(result, length(x))
    expect_all_true(is.na(result))
    ## min_obs < n (partial) = identical slope values
    result <- rolling_slope(x, span = 10, partial = TRUE)
    expect_length(result, length(x))
    expect_equal(var(result), 0)
})

test_that("rolling_slope handles irregular time series with span", {
    x <- c(1, 3, 2, 5, 8)
    t <- c(0, 1, 3, 6, 10)
    span <- 4
    idx_list <- list()
    coef_result <- rep(NA, length(t))
    slope_result <- rep(NA, length(t))
    for (i in seq_along(t)) {
        t_idx <- t[i]
        idx <- t >= t_idx & t <= t_idx + span
        idx_list[[i]] <- which(idx)
        coef_result[i] <- coef(lm(x[idx] ~ t[idx]))[[2L]]
        slope_result[i] <- slope(x[idx], t[idx])
    }

    result <- rolling_slope(
        x,
        t,
        span = 4,
        align = "left",
        partial = TRUE
    )
    expect_equal(length(result), length(x))
    expect_equal(result, coef_result)
    expect_equal(result, slope_result)
})

## peak_slope ===========================================================
test_that("peak_slope returns correct structure", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    result <- peak_slope(x, width = 5)

    expect_type(result, "list")
    expect_named(
        result,
        c("slope", "intercept", "y", "t", "idx", "fitted", "window_idx", "model")
    )
    expect_type(result$slope, "double")
    expect_type(result$intercept, "double")
    expect_type(result$fitted, "double")
    expect_true(is.numeric(result$t))
    expect_type(result$idx, "integer")
    expect_type(result$window_idx, "integer")
    expect_length(result$window_idx, 5L)
})

test_that("peak_slope propagates auto-detected direction", {
    ## direction detection logic tested in test-detect_direction.R
    x_pos <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    expect_gt(peak_slope(x_pos, width = 5, verbose = FALSE)$slope, 0)

    x_neg <- rev(x_pos)
    expect_lt(peak_slope(x_neg, width = 5, verbose = FALSE)$slope, 0)
})

test_that("peak_slope respects manual direction", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    result_pos <- peak_slope(
        x,
        width = 5,
        direction = "positive",
        verbose = FALSE
    )
    expect_gt(result_pos$slope, 0)

    result_neg <- peak_slope(
        x,
        width = 2,
        direction = "negative",
        verbose = FALSE
    )
    expect_lt(result_neg$slope, 0)
})

test_that("peak_slope returns NA when no slopes in specified direction", {
    ## perfectly increasing = no negative slopes
    x <- 1:10
    expect_warning(
        result <- peak_slope(x, width = 3, direction = "negative"),
        "No negative slopes"
    )

    expect_all_true(is.na(result))

    ## perfectly decreasing = no positive slopes
    x <- 10:1
    expect_warning(
        peak_slope(x, width = 3, direction = "positive"),
        "No positive slopes"
    )
})

test_that("peak_slope calculates correct intercept and fitted", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)
    result <- peak_slope(x, width = 5, verbose = FALSE)

    ## verify fitted = intercept + slope * t
    expect_equal(result$y, result$intercept + result$slope * result$t)

    ## verify against lm()
    window_x <- x[result$window_idx]
    window_t <- result$window_idx
    lm_fit <- lm(window_x ~ window_t)
    expect_equal(result$slope, coef(lm_fit)[[2L]])
    expect_equal(result$intercept, coef(lm_fit)[[1L]])
})

test_that("peak_slope finds correct peak location", {
    ## construct data with known peak slope location
    x <- c(1, 2, 3, 10, 20, 30, 31, 32, 33)
    result <- peak_slope(x, width = 3)

    ## peak should be at idx 5 (centre of steepest section)
    expect_equal(result$idx, 5L)
    expect_equal(result$t, 5)
    expect_equal(result$window_idx, 4:6)
})

test_that("peak_slope handles NA values", {
    x <- c(1, 3, NA, 5, 8, 7, 9, 12, NA, NA, NA, 17, 18)

    result_na <- peak_slope(x, width = 7, partial = TRUE, na.rm = FALSE)
    expect_all_true(is.na(result_na))
    result_na <- peak_slope(x, width = 7, partial = TRUE, na.rm = TRUE)
    expect_all_false(is.na(result_na))
})

test_that("peak_slope handles edge cases", {
    ## single value
    expect_all_true(is.na(peak_slope(5, width = 3))) |>
        expect_warning("Insufficient valid")

    ## two values
    expect_equal(peak_slope((1:2), width = 2, partial = TRUE)$slope, 1)

    ## all identical values = slope 0, no positive/negative
    expect_all_true(is.na(peak_slope(rep(5, 10), width = 3))) |>
        expect_warning("No positive slopes")

    ## width > n
    x <- c(1, 3, 2, 5, 8)
    expect_all_true(is.na(peak_slope(x, width = 10))) |>
        expect_warning("Insufficient valid")

    ## width > n but partial = TRUE
    expect_all_false(is.na(peak_slope(x, width = 10, partial = TRUE)))
})

test_that("peak_slope works with span", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    result <- peak_slope(x, span = 5, verbose = FALSE)
    expect_all_false(is.na(result))
    expect_gt(result$slope, 0)

    ## irregular t
    t_irreg <- c(0.5, 1.5, 3, 3.5, 4, 5.5, 7, 8, 10, 11)
    result_irreg <- peak_slope(x, t_irreg, span = 3, verbose = FALSE)
    expect_all_false(is.na(result_irreg))
    expect_gt(result_irreg$slope, result$slope)
})

test_that("peak_slope align options work", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    result_centre <- peak_slope(x, width = 5, align = "centre", verbose = FALSE)
    result_left <- peak_slope(x, width = 5, align = "left", verbose = FALSE)
    result_right <- peak_slope(x, width = 5, align = "right", verbose = FALSE)

    expect_false(is.na(result_centre$slope))
    expect_false(is.na(result_left$slope))
    expect_false(is.na(result_right$slope))

    ## should return different peak idx
    expect_gt(result_centre$idx, result_left$idx)
    expect_gt(result_right$idx, result_centre$idx)
    ## should return same length window_idx
    expect_length(result_centre$window_idx, 5L)
    expect_length(result_left$window_idx, 5L)
    expect_length(result_right$window_idx, 5L)
})

test_that("peak_slope verbose messages work", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 14)

    expect_silent(peak_slope(rep(5, 5), width = 5, verbose = FALSE))
    ## auto-detection message
    expect_warning(
        peak_slope(rep(5, 5), width = 5, verbose = TRUE),
        "No.*slopes detected"
    )
})


test_that("peak_slope fitted values match window", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    results <- peak_slope(x, t, width = 5)

    fitted <- results$fitted
    window_idx <- results$window_idx

    expect_length(fitted, length(window_idx))

    # Verify fitted values match slope * t + intercept
    expect_equal(
        fitted,
        results$intercept + results$slope * t[window_idx]
    )
    ## verify fitted values match lm predictions
    expect_equal(
        fitted,
        predict(lm(x[window_idx] ~ t[window_idx])),
        ignore_attr = TRUE
    )
})


## analyse_peak_slope ====================================================
test_that("analyse_peak_slope returns correct structure", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- c(1, 3, NA, 5, 8, 7, 9, 12, NA, NA, NA, 17, 18)
    t <- as.numeric(seq_along(x))

    df <- create_mnirs_data(
        data = data.frame(time = t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "time"
    )
    results <- analyse_peak_slope(df, width = 5, verbose = FALSE)
    # attributes(results)

    expect_s3_class(results, "data.frame")
    expect_equal(nrow(results), 2)
    expect_named(results, c(
        "interval", "nirs_channels", "time_channel",
        "slope", "intercept", "fitted", "peak_slope_time", "idx"
    ))

    expect_type(results$nirs_channels, "character")
    expect_equal(results$nirs_channels, c("x", "q"))
    expect_type(results$slope, "double")
    expect_type(results$intercept, "double")
    expect_type(results$fitted, "double")
    expect_type(results$peak_slope_time, "double")
    expect_type(results$idx, "integer")

    ## metadata carried as attributes
    expect_type(attr(results, "model"), "list")
    expect_true(inherits(attr(results, "model")$x, "lm"))
    expect_type(attr(results, "fitted_data"), "list")
    expect_s3_class(attr(results, "fitted_data")$x, "data.frame")
    expect_equal(
        colnames(attr(results, "fitted_data")$x),
        c("window_idx", "fitted")
    )
    expect_s3_class(attr(results, "channel_args"), "data.frame")
    expect_s3_class(attr(results, "diagnostics"), "data.frame")
})

test_that("analyse_peak_slope handles edge cases", {
    ## single value per channel errors because validate_nirs_channels wants >=2 valid values
    df <- create_mnirs_data(
        data = data.frame(t = 1, x = 5, q = 10),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    expect_error(
        analyse_peak_slope(df, width = 3, verbose = FALSE),
        "valid.*numeric"
    )

    ## all identical values = NA (no positive/negative slopes)
    df <- create_mnirs_data(
        data = data.frame(t = 1:10, x = rep(5, 10), q = rep(10, 10)),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    results <- analyse_peak_slope(df, width = 3, verbose = FALSE)
    expect_all_true(is.na(results$slope))

    ## all NA_real_ values errors because validate_nirs_channels wants >=2 valid values
    df <- create_mnirs_data(
        data = data.frame(t = 1:5, x = rep(NA_real_, 5), q = rep(NA_real_, 5)),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    expect_error(
        analyse_peak_slope(df, width = 3, verbose = FALSE),
        "valid.*numeric"
    )

    ## all t values identical = NA
    df <- create_mnirs_data(
        data = data.frame(t = rep(1, 5), x = 1:5, q = 1:5),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    results <- analyse_peak_slope(df, width = 3, verbose = FALSE)
    expect_all_true(is.na(results$slope))

    ## NaN & Inf removed, valid slopes calculated
    x <- c(1, NaN, 3, Inf, 5, 7, 9)
    q <- c(2, 4, Inf, 6, NaN, 8, 10)
    df <- create_mnirs_data(
        data = data.frame(t = seq_along(x), x = x, q = q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )
    # plot.mnirs(df, points = TRUE, na.omit = TRUE)

    results <- analyse_peak_slope(
        df,
        width = 3,
        partial = TRUE,
        verbose = FALSE
    )
    
    expect_all_false(is.na(results$slope))
    expect_all_true(results$slope > 0)
    expect_equal(attr(results, "fitted_data")$x$window_idx, c(3, 4, 5))
    expect_equal(attr(results, "fitted_data")$x$fitted, c(5, 7, 9))
    expect_equal(attr(results, "fitted_data")$q$window_idx, c(1, 2))
    expect_equal(attr(results, "fitted_data")$q$fitted, c(2, 4))
})

test_that("analyse_peak_slope validates data structure", {
    ## non-data frame input errors
    expect_error(analyse_peak_slope(1:5, width = 3), "data frame")

    ## missing nirs_channels
    df <- data.frame(t = 1:5, x = 1:5)
    expect_error(analyse_peak_slope(df, width = 3), "nirs_channels")

    ## missing time_channel
    df <- create_mnirs_data(
        data = data.frame(t = 1:5, x = 1:5),
        nirs_channels = "x",
        time_channel = "t"
    )
    attr(df, "time_channel") <- NULL
    expect_error(
        analyse_peak_slope(df, width = 3, verbose = FALSE),
        "time_channel"
    )
})

test_that("analyse_peak_slope handles width/span validation", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    df <- create_mnirs_data(
        data = data.frame(t = seq_along(x), x = x),
        nirs_channels = "x",
        time_channel = "t"
    )

    ## missing both width and span
    expect_error(analyse_peak_slope(df, verbose = FALSE), "width.*span")

    ## width > n without partial = NA
    results <- analyse_peak_slope(df, width = 50, verbose = FALSE)
    expect_all_true(is.na(results$slope))

    ## width > n with partial = valid slope
    results <- analyse_peak_slope(
        df,
        width = 50,
        partial = TRUE,
        verbose = FALSE
    )
    expect_all_false(is.na(results$slope))
})

test_that("analyse_peak_slope processes multiple channels independently", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- rev(x) # decreasing trend
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )

    results <- analyse_peak_slope(df, width = 5, verbose = FALSE)

    # x has positive trend, q has negative trend
    expect_gt(results$slope[results$nirs_channels == "x"], 0)
    expect_lt(results$slope[results$nirs_channels == "q"], 0)

    # Different peak locations
    expect_false(
        results$idx[results$nirs_channels == "x"] ==
            results$idx[results$nirs_channels == "q"]
    )
})

test_that("analyse_peak_slope channel_args override defaults", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    data <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )

    # Different widths per channel
    results <- analyse_peak_slope(
        data,
        width = list(3, q = 7),
        verbose = FALSE
    )

    ## different widths produce different window sizes
    expect_length(attr(results, "fitted_data")$x$window_idx, 3) # x uses width = 3
    expect_length(attr(results, "fitted_data")$q$window_idx, 7) # q uses width = 7
})

test_that("analyse_peak_slope channel_args sets different directions", {
    x <- c(1, 3, 5, 7, 9, 8, 6, 4, 2) # mixed trend
    q <- c(1, 3, 5, 7, 9, 8, 6, 4, 2)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )

    results <- analyse_peak_slope(
        df,
        width = 3,
        direction = list("positive", q = "negative"),
        verbose = FALSE
    )

    expect_gt(results$slope[results$nirs_channels == "x"], 0)
    expect_lt(results$slope[results$nirs_channels == "q"], 0)
})

test_that("analyse_peak_slope handles single channel", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = "x",
        time_channel = "t"
    )

    results <- analyse_peak_slope(df, width = 5, verbose = FALSE)

    expect_equal(nrow(results), 1)
    expect_equal(results$nirs_channels, "x")
    expect_no_match(results$nirs_channels, "q")
    expect_contains(names(attr(results, "fitted_data")), "x")
    expect_no_match(names(attr(results, "fitted_data")), "q")
})

test_that("analyse_peak_slope preserves channel order", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    r <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q, r),
        nirs_channels = c("r", "x", "q"),
        time_channel = "t"
    )

    results <- analyse_peak_slope(df, width = 5, verbose = FALSE)
    expect_equal(results$nirs_channels, c("r", "x", "q"))
    expect_equal(names(attr(results, "fitted_data")), c("r", "x", "q"))
})

test_that("analyse_peak_slope fitted values match window", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x),
        nirs_channels = "x",
        time_channel = "t"
    )

    results <- analyse_peak_slope(df, width = 5, verbose = FALSE)

    fitted <- attr(results, "fitted_data")$x$fitted
    window_idx <- attr(results, "fitted_data")$x$window_idx

    expect_length(fitted, length(window_idx))

    ## verify fitted values match slope * t + intercept
    expect_equal(
        fitted,
        results$intercept + results$slope * t[window_idx]
    )
    ## verify fitted values match lm predictions
    expect_equal(
        fitted,
        predict(lm(x[window_idx] ~ t[window_idx])),
        ignore_attr = TRUE
    )
})

test_that("analyse_peak_slope channel_args stores global arguments", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )

    results <- analyse_peak_slope(
        df,
        width = 5,
        direction = "positive",
        verbose = FALSE
    )

    ## channel_args is a data.frame attribute with one row per channel
    ca <- attr(results, "channel_args")
    expect_s3_class(ca, "data.frame")
    expect_equal(nrow(ca), 2)
    expect_true("nirs_channels" %in% names(ca))

    ## each row contains the arguments used
    x_args <- ca[ca$nirs_channels == "x", ]
    expect_true(all(c(
        "nirs_channels", "width", "span", "align",
        "direction", "partial"
    ) %in% names(x_args)))
    expect_equal(x_args$width, 5)
    expect_true(is.na(x_args$span))
    expect_equal(x_args$align, "centre")
    expect_equal(x_args$direction, "positive")
    expect_false(x_args$partial)
})

test_that("analyse_peak_slope channel_args stores overrides", {
    x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    q <- rev(x)
    t <- seq_along(x)

    df <- create_mnirs_data(
        data = data.frame(t, x, q),
        nirs_channels = c("x", "q"),
        time_channel = "t"
    )

    results <- analyse_peak_slope(
        df,
        width = list(3, q = 7),
        direction = list("positive", q = "negative"),
        verbose = FALSE
    )

    ## x uses defaults
    ca <- attr(results, "channel_args")
    x_args <- ca[ca$nirs_channels == "x", ]
    expect_equal(x_args$width, 3)
    expect_equal(x_args$align, "centre")
    expect_equal(x_args$direction, "positive")
    expect_false(x_args$partial)

    ## q uses overrides
    q_args <- ca[ca$nirs_channels == "q", ]
    expect_equal(q_args$width, 7)
    expect_equal(q_args$align, "centre")
    expect_equal(q_args$direction, "negative")
    expect_false(q_args$partial)
})


## integration ============================================================
test_that("rolling_slope works visually", {
    # x <- c(1, 3, 2, 5, 8, 7, 9, 12, 11, 15, 14, 17, 18)
    x <- c(1, 3, NA, 5, 8, 7, 9, 12, NA, NA, NA, 17, 18)
    t <- seq_along(x)
    window_idx <- compute_local_windows(t, width = 3)

    ## partial = TRUE, na.rm = TRUE
    slopes <- vapply(window_idx, \(.idx) {
        df <- data.frame(t = t[.idx], x = x[.idx])
        tryCatch(coef(lm(x ~ t, df))[[2L]], error = \(e) NA_real_)
        }, numeric(1))

    expect_equal(
        rolling_slope(x, t, width = 3, partial = TRUE, na.rm = TRUE),
        slopes
    )

    ## partial = FALSE, na.rm = FALSE
    slopes <- vapply(window_idx, \(.idx) {
        df <- data.frame(t = t[.idx], x = x[.idx])
        if (length(df$t) < 3 || anyNA(df$x)) return(NA_real_)
        tryCatch(coef(lm(x ~ t, df))[[2L]], error = \(e) NA_real_)
        }, numeric(1))

    expect_equal(
        rolling_slope(x, t, width = 3, partial = FALSE, na.rm = FALSE),
        slopes
    )

    ## partial = FALSE, na.rm = TRUE
    slopes <- vapply(window_idx, \(.idx) {
        df <- data.frame(t = t[.idx], x = x[.idx])
        if (length(df$t) < 3) return(NA_real_)
        tryCatch(coef(lm(x ~ t, df))[[2L]], error = \(e) NA_real_)
        }, numeric(1))

    expect_equal(
        rolling_slope(x, t, width = 3, partial = FALSE, na.rm = TRUE),
        slopes
    )

    ## visual confirmation
    # library(ggplot2)
    # fitted_data <- lapply(window_idx, \(.idx) {
    #     df <- data.frame(t = t[.idx], x = x[.idx])
    #     fitted <- tryCatch(predict(lm(x ~ t, df), df), error = \(e) NA_real_)
    #     df$fitted <- if (!is.null(attr(fitted, "non-estim"))) NA_real_ else fitted
    #     df
    # })[2:12]

    # ggplot(data.frame(t, x), aes(t, x)) +
    #     theme_mnirs() +
    #     geom_line(linewidth = 2) +
    #     lapply(fitted_data, \(.df) {
    #         list(
    #             geom_line(
    #                 data = .df,
    #                 aes(t, fitted),
    #                 colour = "dodgerblue"
    #             ),
    #             geom_point(
    #                 data = .df,
    #                 aes(t, fitted),
    #                 colour = c("dodgerblue", "red", "dodgerblue"),
    #                 shape = c(19, 18, 19),
    #                 size = c(2, 5, 2)
    #             )
    #         )
    #     })
})

test_that("rolling_slope works visually on example_mnirs() data", {
    skip("visual check for rolling_slope")

    data <- read_mnirs(
        example_mnirs("moxy_ramp"),
        nirs_channels = c(
            smo2_right = "SmO2 Live(2)"
        ),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        filter_mnirs(method = "butterworth", order = 4, W = 0.05)

    data$slope <- rolling_slope(
        x = data$smo2_right,
        t = data$time,
        span = 10,
        partial = FALSE,
    )

    library(ggplot2)
    plot(data) +
        ggplot2::geom_hline(
            yintercept = 50,
            linetype = "dotted",
            colour = "red"
        ) +
        ggplot2::geom_line(ggplot2::aes(y = 50 + data$slope * 10), colour = "red")

    data <- read_mnirs(
        example_mnirs("train.red"),
        nirs_channels = c(
            smo2_right = "SmO2 unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        verbose = FALSE
    ) |>
        filter_mnirs(method = "butterworth", order = 4, W = 0.05)

    data$slope <- rolling_slope(
        x = data$smo2_right,
        t = data$time,
        span = 10,
        partial = FALSE
    )

    plot(data) +
        ggplot2::geom_hline(
            yintercept = 65,
            linetype = "dotted",
            colour = "red"
        ) +
        ggplot2::geom_line(ggplot2::aes(y = 65 + data$slope * 10), colour = "red")
})


## benchmark ===========================================================
test_that("analyse_frollslope benchmark", {
    ## baselne established from documented example on initial run;
    ## fails if itr/sec regresses by >10%
    skip("benchmark analyse_peak_slope baseline test")

    df <- read_mnirs(
        example_mnirs("train.red"),
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        zero_time = TRUE,
        verbose = FALSE
    ) |>
        resample_mnirs(method = "linear", verbose = FALSE) |>
        extract_intervals(
            group_intervals = "distinct",
            start = by_time(368, 1084),
            span = c(-20, 90),
            zero_time = TRUE,
            verbose = FALSE
        )
    # df <- df[[1L]]

    # bench::mark(
    #     slope = slope(df$smo2_left, df$time, na.rm = TRUE),
    #     check = TRUE,
    #     iterations = 500
    # )

    # bench::mark(
    #     rolling_slope = rolling_slope(
    #         df[[1L]]$smo2_left,
    #         df[[1L]]$time,
    #         span = 10,
    #         align = "right",
    #         partial = FALSE,
    #     ),
    #     check = TRUE,
    #     iterations = 50
    # )

    # bench::mark(
    #     analyse_peak_slope = analyse_peak_slope(
    #         df[[1L]],
    #         nirs_channels = c(smo2_left),
    #         span = 10,
    #         partial = FALSE,
    #         verbose = FALSE
    #     ),
    #     iterations = 20,
    #     check = FALSE
    # )

    # bench::mark(
    #     peak_slope = analyse_kinetics(
    #         df[[1L]],
    #         nirs_channels = c(smo2_left),
    #         method = "peak_slope",
    #         span = 10,
    #         partial = FALSE,
    #         verbose = FALSE
    #     ),
    #     monoexponential = analyse_kinetics(
    #         df[[1L]],
    #         nirs_channels = c(smo2_left),
    #         method = "monoexp",
    #         use_time_delay = TRUE,
    #         partial = FALSE,
    #         verbose = FALSE
    #     ),
    #     check = FALSE,
    #     iterations = 5
    # )
})
