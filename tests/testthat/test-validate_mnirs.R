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


## validate_numeric() ==================================
test_that("validate_numeric accepts valid numeric values", {
    expect_silent(validate_numeric(5))
    expect_silent(validate_numeric(c(1, 2, 3)))
    expect_silent(validate_numeric(1.5))
})

test_that("validate_numeric checks element count", {
    x <- c(1, 2)
    expect_error(validate_numeric(x, elements = 1), "numeric")
    expect_silent(validate_numeric(5, elements = 1))
    expect_silent(validate_numeric(c(1, 2, 3), elements = 3))
    ## fewer elements than specified fails
    expect_error(validate_numeric(1, elements = 2), "numeric")
})

test_that("validate_numeric checks range with inclusive bounds", {
    expect_silent(validate_numeric(5, range = c(0, 10)))
    expect_silent(validate_numeric(0, range = c(0, 10)))
    expect_silent(validate_numeric(10, range = c(0, 10)))
    expect_error(validate_numeric(-1, range = c(0, 10)), "numeric")
    expect_error(validate_numeric(11, range = c(0, 10)), "numeric")
})

test_that("validate_numeric checks range with exclusive bounds", {
    expect_silent(validate_numeric(5, range = c(0, 10), inclusive = FALSE))
    expect_error(
        validate_numeric(0, range = c(0, 10), inclusive = FALSE),
        "numeric"
    )
    expect_error(
        validate_numeric(10, range = c(0, 10), inclusive = FALSE),
        "numeric"
    )
})

test_that("validate_numeric rejects non-numeric input", {
    expect_error(validate_numeric("text"), "numeric")
    expect_error(validate_numeric(TRUE), "numeric")
    expect_error(validate_numeric(list(1, 2))) ## "not implemented for type `list`"
})

test_that("validate_numeric handles NA, NaN, Inf", {
    expect_silent(validate_numeric(c(1, NA)))
    expect_error(validate_numeric(NA), "numeric")
    expect_error(validate_numeric(NaN), "numeric")
    expect_silent(validate_numeric(Inf))
    expect_error(validate_numeric(Inf, range = c(0, 1)), "numeric")
    expect_silent(validate_numeric(c(1, NA, 3), elements = 2))
})

test_that("validate_numeric handles NULL", {
    expect_silent(validate_numeric(NULL))
})

test_that("validate_numeric reports errors from the caller's `env`", {
    ## a named caller env attributes the abort to the user-facing function
    caller <- function() {
        validate_numeric("text", env = rlang::current_env())
    }
    err <- expect_error(caller(), "numeric")
    expect_equal(rlang::call_name(conditionCall(err)), "caller")
})

## validate_mnirs_data() ========================================
test_that("validate_mnirs_data() accepts valid data frames", {
    data <- create_test_data()
    expect_silent(validate_mnirs_data(data))
})

test_that("validate_mnirs_data() rejects non-data frames", {
    expect_error(
        validate_mnirs_data(list(a = 1, b = 2)),
        "must be a data frame"
    )
    expect_error(validate_mnirs_data(c(1, 2, 3)), "must be a data frame")
})

test_that("validate_mnirs_data() rejects data frames with < 2 columns", {
    expect_error(validate_mnirs_data(data.frame(x = 1:10)), "at least")
})


## parse_channel_name() ============================================
test_that("parse_channel_name() returns NULL for NULL input", {
    data <- data.frame(a = 1, b = 2)
    result <- parse_channel_name(rlang::quo(NULL), data)
    expect_null(result)
})

test_that("parse_channel_name() handles already-evaluated character", {
    data <- data.frame(a = 1, b = 2)
    channel <- "a"
    result <- parse_channel_name(channel = rlang::enquo(channel), data)
    expect_equal(result, "a")
})

test_that("parse_channel_name() handles already-evaluated character vector", {
    data <- data.frame(a = 1, b = 2, c = 3)
    ch <- c("a", "b")
    result <- parse_channel_name(rlang::enquo(ch), data)
    expect_equal(result, c("a", "b"))
})

test_that("parse_channel_name() handles already-evaluated list", {
    data <- data.frame(a = 1, b = 2, c = 3)
    ch <- list(c("a", "b"), "c")
    result <- parse_channel_name(rlang::enquo(ch), data)
    expect_equal(result, list(c("a", "b"), "c"))
})

test_that("parse_channel_name() handles bare symbol matching column", {
    data <- data.frame(smo2 = 1, thb = 2)
    env <- rlang::current_env()
    expr <- rlang::quo(smo2)
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, "smo2")
})

test_that("parse_channel_name() handles bare symbol as external object", {
    data <- data.frame(a = 1, b = 2, c = 3)
    channels <- c("a", "b")
    result <- parse_channel_name(rlang::enquo(channels), data)
    expect_equal(result, c("a", "b"))
})

test_that("parse_channel_name() handles external list object", {
    data <- data.frame(a = 1, b = 2, c = 3)
    channels <- list(c("a", "b"), "c")
    result <- parse_channel_name(rlang::enquo(channels), data)
    expect_equal(result, list(c("a", "b"), "c"))
})

test_that("parse_channel_name() handles c() with bare symbols", {
    data <- data.frame(smo2 = 1, thb = 2, hhr = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(c(smo2, thb))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, c("smo2", "thb"))
})

test_that("parse_channel_name() handles list() with bare symbols", {
    data <- data.frame(smo2 = 1, thb = 2, hhr = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(list(c(smo2, thb), hhr))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, list(c("smo2", "thb"), "hhr"))
})

test_that("parse_channel_name() handles list() with mixed input", {
    data <- data.frame(a = 1, b = 2, c = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(list(c("a", "b"), c))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, list(c("a", "b"), "c"))
})

test_that("parse_channel_name() handles tidyselect starts_with()", {
    data <- data.frame(smo2_left = 1, smo2_right = 2, thb = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(tidyselect::starts_with("smo2"))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, c("smo2_left", "smo2_right"))
})

test_that("parse_channel_name() handles tidyselect ends_with()", {
    data <- data.frame(smo2_left = 1, thb_left = 2, smo2_right = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(tidyselect::ends_with("left"))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, c("smo2_left", "thb_left"))
})

test_that("parse_channel_name() handles tidyselect matches()", {
    data <- data.frame(left_smo2 = 1, smo2_right = 2, thb = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(tidyselect::matches("smo2"))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, c("left_smo2", "smo2_right"))
})

test_that("parse_channel_name() handles tidyselect in list()", {
    data <- data.frame(smo2_a = 1, smo2_b = 2, thb = 3)
    env <- rlang::current_env()
    expr <- rlang::quo(list(tidyselect::starts_with("smo2"), thb))
    result <- parse_channel_name(expr, data, env = env)
    expect_equal(result, list(c("smo2_a", "smo2_b"), "thb"))
})

test_that("parse_channel_name() returns non-existent colname symbol", {
    data <- data.frame(a = 1, b = 2)
    env <- rlang::new_environment()
    channel <- rlang::new_quosure(rlang::sym("nonexistent"), env = env)
    result <- parse_channel_name(channel, data, env = env)
    ## should pass through nonexistent colname to validate_nirs_channels
    expect_equal(result, "nonexistent")
})

test_that("parse_channel_name returns NULL on logical channel", {
    data <- data.frame(a = 1, b = 2)
    env <- rlang::new_environment()
    expect_null(parse_channel_name(rlang::quo(TRUE), data, env))
    list_null <- parse_channel_name(rlang::quo(list(TRUE)), data, env)
    expect_true(is.list(list_null))
    expect_null(unlist(list_null))
})



## validate_nirs_channels() ========================================
test_that("validate_nirs_channels() uses metadata when NULL", {
    data <- create_test_data()
    result <- validate_nirs_channels(NULL, data, verbose = FALSE)
    expect_equal(result, c("nirs1", "nirs2"))
})

test_that("validate_nirs_channels() uses explicit channels when provided", {
    data <- create_test_data()
    result <- validate_nirs_channels("nirs1", data)
    expect_equal(result, "nirs1")
})

test_that("validate_nirs_channels() works with nirs_channels = list()", {
    data <- create_test_data()
    nirs_vec <- c("nirs1", "nirs2")
    result <- validate_nirs_channels(nirs_vec, data)
    expect_equal(result, nirs_vec)

    result <- validate_nirs_channels(enquo(nirs_vec), data)
    expect_equal(result, nirs_vec)

    attr(data, "nirs_channels") <- nirs_vec
    expect_message(
        result <- validate_nirs_channels(
            NULL, data, verbose = TRUE, as_list = TRUE
        ),
        "`nirs_channels`.*grouped"
    )
    expect_equal(result, nirs_vec)

    nirs_list <- list(c("nirs1", "nirs2"), "nirs3")
    result <- validate_nirs_channels(nirs_list, data, as_list = TRUE)
    expect_equal(result, nirs_list)

    result <- validate_nirs_channels(enquo(nirs_list), data, as_list = TRUE)
    expect_equal(result, nirs_list)
})

test_that("validate_nirs_channels() errors when not in metadata or provided", {
    data <- create_test_data(add_metadata = FALSE)
    expect_error(validate_nirs_channels(NULL, data), "not detected in metadata")
})

test_that("validate_nirs_channels() errors when columns don't exist", {
    data <- create_test_data()
    expect_error(validate_nirs_channels("nonexistent", data), "match exactly")
})

test_that("validate_nirs_channels() errors when col symbols don't exist", {
    data <- data.frame(a = 1, b = 2)
    env <- rlang::new_environment()
    channel <- rlang::new_quosure(rlang::sym("nonexistent"), env = env)
    expect_error(validate_nirs_channels(channel, data), "match exactly")
})

test_that("validate_nirs_channels() errors for non-numeric channels", {
    data <- create_test_data()
    data$nirs1 <- as.character(data$nirs1)
    expect_error(
        validate_nirs_channels(c("nirs1", "nirs2"), data),
        "must contain valid.*numeric"
    )
})

test_that("validate_nirs_channels() errors when < 2 valid values", {
    data <- create_test_data()
    data$nirs1 <- c(1, rep(NA, nrow(data) - 1))
    expect_error(
        validate_nirs_channels(c("nirs1", "nirs2"), data),
        "must contain valid.*numeric"
    )
})

test_that("validate_nirs_channels() silently flattens list on default path", {
    data <- create_test_data()
    nirs_list <- list(c("nirs1", "nirs2"), "nirs3")
    expect_no_message(
        result <- validate_nirs_channels(
            nirs_list, data, verbose = TRUE, as_list = FALSE
        )
    )
    expect_equal(result, c("nirs1", "nirs2", "nirs3"))
})


## validate_time_channel() ========================================
test_that("validate_time_channel() uses metadata when NULL", {
    data <- create_test_data()
    result <- validate_time_channel(NULL, data)
    expect_equal(result, "time")
})

test_that("validate_time_channel() uses explicit channel when provided", {
    data <- create_test_data()
    data$time_new <- data$time
    result <- validate_time_channel("time_new", data)
    expect_equal(result, "time_new")
})

test_that("validate_time_channel() errors when not in metadata or provided", {
    data <- create_test_data(add_metadata = FALSE)
    expect_error(validate_time_channel(NULL, data), "not detected in metadata")
})

test_that("validate_time_channel() errors when column doesn't exist", {
    data <- create_test_data()
    expect_error(validate_time_channel("nonexistent", data), "match exactly")
})

test_that("validate_time_channel() errors for non-numeric channel", {
    data <- create_test_data()
    data$time <- as.character(data$time)
    expect_error(validate_time_channel("time", data), "must contain valid.*numeric")
})

test_that("validate_time_channel() errors when < 2 valid values", {
    data <- create_test_data()
    data$time <- c(1, rep(NA, nrow(data) - 1))
    expect_error(validate_time_channel("time", data), "must contain valid.*numeric")
})


## validate_event_channel() ========================================
test_that("validate_event_channel() uses metadata when NULL", {
    data <- create_test_data()
    result <- validate_event_channel(NULL, data)
    expect_equal(result, "event")
})

test_that("validate_event_channel() uses explicit channel when provided", {
    data <- create_test_data()
    data$lap <- data$event
    result <- validate_event_channel("lap", data)
    expect_equal(result, "lap")
})

test_that("validate_event_channel() errors when not in metadata or provided", {
    data <- create_test_data(add_metadata = FALSE)
    expect_error(validate_event_channel(NULL, data), "not detected in metadata")

    ## no metadata, required = FALSE
    expect_equal(validate_event_channel(NULL, data, required = FALSE), NULL)
})

test_that("validate_event_channel() errors when column doesn't exist", {
    data <- create_test_data()
    expect_error(validate_event_channel("nonexistent", data), "match exactly")
})

test_that("validate_event_channel() errors when all NA", {
    data <- create_test_data()
    data_one <- data[1:(nrow(data) - 1), ] ## one valid should work
    expect_equal(validate_event_channel("event", data_one), "event")
    data_na <- data[2:(nrow(data) - 1), ]
    expect_error(
        validate_event_channel("event", data_na),
        "must contain valid"
    )
})

test_that("validate_event_channel() works for quosures", {
    data <- create_test_data()
    result <- validate_event_channel(rlang::quo("event"), data)
    expect_equal(result, "event")
})

test_that("validate_event_channel() accepts integer column", {
    data <- create_test_data()
    ## create_test_data() uses c(1, rep(NA, nrow - 2), 2) — integerish
    data$lap <- c(1L, rep(NA_integer_, nrow(data) - 2L), 2L)
    result <- validate_event_channel("lap", data)
    expect_equal(result, "lap")
})

test_that("validate_event_channel() error conditions", {
    data <- create_test_data()

    ## errors for non-integerish numeric column
    data$bad_col <- c(1.5, rep(NA, nrow(data) - 2), 2.5)
    expect_error(validate_event_channel("bad_col", data), "integer")

    ## errors when integer column is all NA
    data$lap <- rep(NA_integer_, nrow(data))
    expect_error(validate_event_channel("lap", data), "must contain valid")
})



## within() ===============================================================
test_that("within() handles basic inclusive range (default)", {
    expect_equal(within(5, c(1, 10)), TRUE)
    expect_equal(within(1, c(1, 10)), TRUE)
    expect_equal(within(10, c(1, 10)), TRUE)
    expect_equal(within(0, c(1, 10)), FALSE)
    expect_equal(within(11, c(1, 10)), FALSE)
})

test_that("within() handles vectorised inputs", {
    expect_equal(
        within(c(0, 1, 5, 10, 11), c(1, 10)),
        c(FALSE, TRUE, TRUE, TRUE, FALSE)
    )
    expect_equal(
        within(1:10, c(3, 7)),
        c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
    )
})

test_that("within() handles non-numeric `vec` values", {
    expect_error(
        within(1:5, c("A", "B")),
        "valid.*numeric"
    )
})

test_that("within() handles inclusive/exclusive", {
    ## both exclusive
    expect_false(within(1, c(1, 10), inclusive = FALSE))
    expect_false(within(10, c(1, 10), inclusive = FALSE))
    expect_true(within(5, c(1, 10), inclusive = FALSE))
    expect_equal(
        within(c(1, 5, 10), c(1, 10), inclusive = FALSE),
        c(FALSE, TRUE, FALSE)
    )

    ## left inclusive
    expect_true(within(1, c(1, 10), inclusive = "left"))
    expect_false(within(10, c(1, 10), inclusive = "left"))
    expect_equal(
        within(c(1, 5, 10), c(1, 10), inclusive = "left"),
        c(TRUE, TRUE, FALSE)
    )

    ## right inclusive
    expect_false(within(1, c(1, 10), inclusive = "right"))
    expect_true(within(10, c(1, 10), inclusive = "right"))
    expect_equal(
        within(c(1, 5, 10), c(1, 10), inclusive = "right"),
        c(FALSE, TRUE, TRUE)
    )
})

test_that("within() detects positive non-zero values", {
    expect_true(within(0, c(0, Inf)))
    expect_true(within(Inf, c(0, Inf)))
    expect_false(within(0, c(0, Inf), inclusive = FALSE))
    expect_true(within(1, c(0, Inf), inclusive = FALSE))
    expect_false(within(-1, c(0, Inf)))
    expect_false(within(-1, c(0, Inf), inclusive = FALSE))
    expect_equal(
        within(c(-1, 0, 0.001, 1), c(0, Inf), inclusive = FALSE),
        c(FALSE, FALSE, TRUE, TRUE)
    )
})

test_that("within() handles NA values", {
    expect_error(within(NA, c(1, 10)), "valid.*numeric")
    expect_equal(within(NA_real_, c(1, 10)), NA)
    expect_equal(within(c(1, NA, 5), c(1, 10)), c(TRUE, NA, TRUE))
    expect_false(within(5, c(NA, 10)))
    expect_false(within(5, c(1, NA)))
})

test_that("within() handles infinite values", {
    expect_true(within(Inf, c(1, Inf)), TRUE)
    expect_false(within(Inf, c(1, Inf), inclusive = FALSE), FALSE)
    expect_true(within(-Inf, c(-Inf, 10)), TRUE)
    expect_false(within(-Inf, c(-Inf, 10), inclusive = FALSE), FALSE)
})

test_that("within() handles negative ranges", {
    expect_true(within(-5, c(-10, -1)))
    expect_true(within(-10, c(-10, -1)))
    expect_false(within(-11, c(-10, -1)))
})

test_that("within() handles degenerate ranges where left equals right", {
    expect_true(within(5, c(5, 5)), TRUE)
    expect_false(within(5, c(5, 5), inclusive = FALSE))
    expect_false(within(5, c(5, 5), inclusive = "left"))
    expect_false(within(5, c(5, 5), inclusive = "right"))
})

test_that("within is equivalent to dplyr::between()", {
    expect_equal(
        within(c(0, 5, 10, 15), c(1, 10)),
        dplyr::between(c(0, 5, 10, 15), 1, 10)
    )

    expect_equal(within(NA_real_, c(1, 10)), dplyr::between(NA_real_, 1, 10))
    expect_equal(
        within(c(1, NA, 5), c(1, 10)),
        dplyr::between(c(1, NA, 5), 1, 10)
    )
})


## estimate_sample_rate() ======================================
test_that("estimate_sample_rate works correctly", {
    # Regular 100 Hz data
    time_100hz <- seq(0, 1, by = 0.01)
    expect_equal(estimate_sample_rate(time_100hz), 100)

    # 50 Hz data
    time_50hz <- seq(0, 1, by = 0.02)
    expect_equal(estimate_sample_rate(time_50hz), 50)

    ## irregular sampling
    set.seed(13)
    x <- seq(0, 1, by = 0.02)
    x <- x + rnorm(length(x), 0, 0.001)
    expect_true(estimate_sample_rate(x) %in% c(49.5, 50, 50.5))
    expect_true(
        all.equal(estimate_sample_rate(x), 50, tolerance = 1, scale = 1)
    )
    expect_true(
        all.equal(estimate_sample_rate(x), 50, tolerance = 1, scale = 1)
    )

    # With NAs in diffs (not in x directly, as diff removes one element)
    time_with_gaps <- c(0, 0.01, NA, 0.03, 0.04)
    expect_type(estimate_sample_rate(time_with_gaps), "double")

    # Edge case: single diff value
    expect_equal(estimate_sample_rate(c(0, 0.01)), 100)
    expect_error(estimate_sample_rate(NA), "numeric")

    ## edge case sample rate undetectable returns NULL
    expect_error(
        estimate_sample_rate(c(1, 1, 1, 1, 1)),
        "Unable to estimate"
    )
    expect_error(
        estimate_sample_rate(Inf),
        "Unable to estimate"
    )
})


## validate_sample_rate() ========================================
test_that("validate_sample_rate() uses metadata when NULL", {
    data <- create_test_data(sample_rate = 15)
    result <- validate_sample_rate(data, "time", NULL)
    expect_equal(result, 15)
})

test_that("validate_sample_rate() uses explicit value with warning when provided", {
    data <- create_test_data(sample_rate = 10)
    expect_equal(validate_sample_rate(data, "time", 20), 20) |>
        expect_warning("appears to be inconsistent")
    ## uses explicit without warning
    expect_equal(validate_sample_rate(data, "time", 20, verbose = FALSE), 20) |>
        expect_silent()
})

test_that("validate_sample_rate() does not warn for integer time_channel", {
    data <- create_test_data(time_max = 100, sample_rate = 1)
    expect_silent(validate_sample_rate(data, "time", 5, FALSE))
})

test_that("validate_sample_rate() estimates from time_channel when NULL", {
    rate = 9.5
    data <- create_test_data(
        time_max = 10,
        sample_rate = rate,
        add_metadata = FALSE
    )
    result <- validate_sample_rate(data, "time", NULL, verbose = FALSE)
    expect_equal(result, round(rate))

    rate = 0.8
    data <- create_test_data(
        time_max = 10,
        sample_rate = rate,
        add_metadata = FALSE
    )
    result <- validate_sample_rate(data, "time", NULL, verbose = FALSE)
    expect_equal(result, round(rate))

    rate = 11
    data <- create_test_data(
        time_max = 10,
        sample_rate = rate,
        add_metadata = FALSE
    )
    result <- validate_sample_rate(data, "time", NULL, verbose = FALSE)
    expect_equal(result, 10)

    rate = 44
    data <- create_test_data(
        time_max = 10,
        sample_rate = rate,
        add_metadata = FALSE
    )
    result <- validate_sample_rate(data, "time", NULL, verbose = FALSE)
    expect_equal(result, 50)

    rate = 98
    data <- create_test_data(
        time_max = 10,
        sample_rate = rate,
        add_metadata = FALSE
    )
    result <- validate_sample_rate(data, "time", NULL, verbose = FALSE)
    expect_equal(result, 100)
})

test_that("validate_sample_rate() shows message when estimating", {
    data <- create_test_data(add_metadata = FALSE)
    expect_message(
        validate_sample_rate(data, "time", NULL, verbose = TRUE),
        "Estimated"
    )
})

test_that("validate_sample_rate() errors for non-single non-numeric non-positive", {
    data <- create_test_data()
    expect_error(
        validate_sample_rate(data, "time", "10", FALSE),
        "must be .*numeric"
    )
    expect_error(
        validate_sample_rate(data, "time", NA, FALSE),
        "must be .*numeric"
    )
    expect_error(
        validate_sample_rate(data, "time", 0, FALSE),
        "must be .*positive"
    )
    expect_error(
        validate_sample_rate(data, "time", -5, FALSE),
        "must be .*positive"
    )
    expect_error(
        validate_sample_rate(data, "time", c(10, 20), FALSE),
        "must be .*one-element"
    )
})

## validate_width_span ==============================
test_that("validate_width_span() validates inputs", {
    expect_error(
        validate_width_span(width = NULL, span = NULL),
        "width.*span.*must be defined"
    )

    expect_message(
        validate_width_span(width = 2, span = 1),
        "width.*overrides.*span"
    )

    expect_error(
        validate_width_span(width = -1),
        "width.*valid.*integer"
    )

    expect_error(
        validate_width_span(width = 1.5),
        "width.*valid.*integer"
    )

    expect_error(
        validate_width_span(span = -1),
        "span.*valid.*numeric"
    )
})

test_that("validate_width_span() reports message", {
    expect_error(
        validate_width_span(
            width = NULL, span = NULL, msg = "for this function."
        ),
        "width.*span.*must be defined.*for this function."
    )
})

## validate_x_t =================================
test_that("validate_x_t() validates inputs", {
    expect_error(validate_x_t(x = 1:10, t = 1:5), "numeric.*equal length")
    expect_error(validate_x_t(x = NULL, t = 1:10), "valid.*numeric")
    expect_error(validate_x_t(x = 1:10, t = NULL), "valid.*numeric")
    expect_error(validate_x_t(x = NA_real_, t = 1), "valid.*numeric")
    expect_silent(validate_x_t(x = NA_real_, t = 1, allow_na = TRUE))
})


## validate_start_time() ==============================
test_that("validate_start_time() uses explicit start_time when provided", {
    data <- create_test_data()
    attr(data, "interval_times") <- 3
    expect_equal(validate_start_time(5, data, data[["time"]]), 5)
})

test_that("validate_start_time() falls back to interval_times attribute", {
    data <- create_test_data()
    attr(data, "interval_times") <- 3
    expect_equal(validate_start_time(NULL, data, data[["time"]]), 3)
})

test_that("validate_start_time() falls back to zero when no start_time or attribute", {
    data <- create_test_data()
    expect_equal(validate_start_time(NULL, data, data[["time"]]), 0)
})

test_that("validate_start_time() rejects non-numeric or multi-element start_time", {
    data <- create_test_data()
    expect_error(validate_start_time("5", data, data[["time"]]), "numeric")
    expect_error(validate_start_time(c(1, 2), data, data[["time"]]), "numeric")
    expect_error(validate_start_time(NA, data, data[["time"]]), "numeric")
})

test_that("validate_start_time() warns and resets when no observations <= start_time", {
    data <- create_test_data()
    expect_warning(
        result <- validate_start_time(-1, data, data[["time"]]),
        "No observations.*time_channel"
    )
    expect_equal(result, data[["time"]][1L])

    ## verbose = FALSE suppresses warning, still resets
    expect_silent(
        result <- validate_start_time(-1, data, data[["time"]], verbose = FALSE)
    )
    expect_equal(result, data[["time"]][1L])
})

test_that("validate_start_time() errors when start_time exceeds time range", {
    data <- create_test_data(time_max = 10)
    expect_error(
        validate_start_time(100, data, data[["time"]]),
        "No observations.*before"
    )
})

## findInt_mnirs() ==================================
test_that("findInt_mnirs provides informative error", {
    expect_error(findInt_mnirs(1, c(3, 1, 2)), "time_channel.*sorted")
    expect_error(findInt_mnirs(1, c(1, NA, 2)), "time_channel.*sorted")
    ## valid input passes through to findInterval()
    expect_equal(findInt_mnirs(2.5, 1:5), findInterval(2.5, 1:5))
})