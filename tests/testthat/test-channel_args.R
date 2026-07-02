## resolve_channel_args ====================================================

test_that("resolve_channel_args broadcasts global values to all channels", {
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = 5, method = "linear"),
        verbose = FALSE
    )

    expect_named(out, c("smo2", "o2hb"))
    expect_equal(out$smo2$width, 5)
    expect_equal(out$o2hb$width, 5)
    expect_equal(out$smo2$method, "linear")
    expect_equal(out$o2hb$method, "linear")
})

test_that("resolve_channel_args applies per-channel overrides", {
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = list(smo2 = 10), span = list(o2hb = 30)),
        verbose = FALSE
    )

    expect_equal(out$smo2$width, 10)
    expect_equal(out$o2hb$span, 30)
    ## omitted channel falls back to the formal default (NULL)
    expect_null(out$smo2$span)
    expect_null(out$o2hb$width)
})

test_that("resolve_channel_args unnamed value sets the global fallback", {
    out <- resolve_channel_args(
        c("smo2", "o2hb", "hhb"),
        args = list(width = list(5, o2hb = 7)),
        verbose = FALSE
    )

    expect_equal(out$smo2$width, 5)
    expect_equal(out$o2hb$width, 7)
    expect_equal(out$hhb$width, 5)
})

test_that("resolve_channel_args lists without channel names stay global", {
    ## purely unnamed multi-element list is a global value
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = list(5, 7)),
        verbose = FALSE
    )
    expect_equal(out$smo2$width, list(5, 7))
    expect_equal(out$o2hb$width, list(5, 7))

    ## single unnamed element with no channel key is global, not a fallback
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = list(10)),
        verbose = FALSE
    )
    expect_equal(out$smo2$width, list(10))
    expect_equal(out$o2hb$width, list(10))
})

test_that("resolve_channel_args informs on omitted channels", {
    ## unnamed channel omits channel and produces a warning
    resolve_channel_args(
        nirs_channels = c("smo2", "o2hb"),
        args = list(width = list(smo2 = 5)),
        verbose = TRUE
    ) |>
        expect_warning("width.*o2hb.*not specified")

    ## an unnamed fallback covers omitted channels: no warning
    expect_no_warning(
        resolve_channel_args(
            c("smo2", "o2hb"),
            args = list(width = list(10, smo2 = 5)),
            verbose = TRUE
        )
    )
})

test_that("resolve_channel_args defaults apply when channel omitted", {
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(method = list(smo2 = "median")),
        defaults = list(method = "linear"),
        choices = list(method = c("linear", "median", "locf", "none")),
        verbose = FALSE
    )

    expect_equal(out$smo2$method, "median")
    expect_equal(out$o2hb$method, "linear")
})

test_that("resolve_channel_args choices match like match.arg", {
    ## a full default vector resolves to its first element
    out <- resolve_channel_args(
        "smo2",
        args = list(method = c("linear", "median", "locf", "none")),
        choices = list(method = c("linear", "median", "locf", "none")),
        verbose = FALSE
    )
    expect_equal(out$smo2$method, "linear")

    ## invalid choice values abort
    expect_error(
        resolve_channel_args(
            "smo2",
            args = list(method = "wrong"),
            choices = list(method = c("linear", "median")),
            verbose = FALSE
        ),
        "method.*must be one of.*linear"
    )
})

test_that("resolve_channel_args warns on non-channel list names", {
    ## named list with no recognised channel names is a per-channel map:
    ## unknown names warned and ignored, channels fall back to defaults
    resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = list(typo = 7)),
        verbose = TRUE
    ) |>
        expect_warning("typo.*not recognised") |>
        expect_warning("width.*not specified")

    out <- suppressWarnings(
        resolve_channel_args(
            c("smo2", "o2hb"),
            args = list(width = list(typo = 7)),
            verbose = TRUE
        )
    )
    expect_null(out$smo2$width)
    expect_null(out$o2hb$width)

    ## unnamed vectors are always global
    out <- resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(invalid_values = c(0, 100)),
        verbose = FALSE
    )
    expect_equal(out$smo2$invalid_values, c(0, 100))
    expect_equal(out$o2hb$invalid_values, c(0, 100))
})

test_that("resolve_channel_args warns on partly matching list names", {
    resolve_channel_args(
        c("smo2", "o2hb"),
        args = list(width = list(smo2 = 5, typo = 7)),
        verbose = TRUE
    ) |>
        expect_warning("typo.*not recognised") |>
        expect_warning("width.*o2hb.*not specified")

    ## unrecognised names are ignored; valid keys resolve per channel
    out <- suppressWarnings(
        resolve_channel_args(
            c("smo2", "o2hb"),
            args = list(width = list(smo2 = 5, typo = 7)),
            verbose = TRUE
        )
    )
    expect_equal(out$smo2$width, 5)
    expect_null(out$o2hb$width)
})

test_that("resolve_channel_args resolves per group with group_channels", {
    groups <- list(oxy = c("smo2_l", "smo2_r"), hhb = "hhb")

    ## group-name key applies to the whole group
    out <- resolve_channel_args(
        c("smo2_l", "smo2_r", "hhb"),
        args = list(to = list(oxy = 0, hhb = 5)),
        group_channels = groups,
        verbose = FALSE
    )
    expect_named(out, c("oxy", "hhb"))
    expect_equal(out$oxy$to, 0)
    expect_equal(out$hhb$to, 5)

    ## member-channel key also applies to the whole group
    out <- resolve_channel_args(
        c("smo2_l", "smo2_r", "hhb"),
        args = list(to = list(smo2_r = 1)),
        group_channels = groups,
        verbose = FALSE
    )
    expect_equal(out$oxy$to, 1)
    expect_null(out$hhb$to)
})

test_that("resolve_channel_args aborts on intra-group conflicts", {
    expect_error(
        resolve_channel_args(
            c("smo2_l", "smo2_r"),
            args = list(to = list(smo2_l = 0, smo2_r = 5)),
            group_channels = list(oxy = c("smo2_l", "smo2_r")),
            verbose = FALSE
        ),
        "conflicting.*group_channels"
    )

    ## equal member values do not conflict
    out <- resolve_channel_args(
        c("smo2_l", "smo2_r"),
        args = list(to = list(smo2_l = 0, smo2_r = 0)),
        group_channels = list(oxy = c("smo2_l", "smo2_r")),
        verbose = FALSE
    )
    expect_equal(out$oxy$to, 0)
})

test_that("resolve_channel_args reports errors from the caller's `env`", {
    ## a named caller env attributes the abort to the user-facing function
    caller <- function() {
        resolve_channel_args(
            c("smo2_l", "smo2_r"),
            args = list(to = list(smo2_l = 0, smo2_r = 5)),
            group_channels = list(smo2 = c("smo2_l", "smo2_r")),
            verbose = FALSE
        )
    }
    err <- expect_error(caller(), "conflicting")
    expect_equal(rlang::call_name(conditionCall(err)), "caller")

    ## invalid choices are also attributed to the caller
    caller_choice <- function() {
        resolve_channel_args(
            "smo2",
            args = list(method = "wrong"),
            choices = list(method = c("linear", "median")),
            verbose = FALSE
        )
    }
    err <- expect_error(caller_choice(), "linear")
    expect_equal(rlang::call_name(conditionCall(err)), "caller_choice")
})


## validate_group_channels =================================================

test_that("validate_group_channels expands string shortcuts", {
    channels <- c("smo2", "o2hb", "hhb")

    ensemble <- validate_group_channels(channels, "ensemble")
    expect_length(ensemble, 1L)
    expect_equal(ensemble[[1L]], channels)

    distinct <- validate_group_channels(channels, "distinct")
    expect_length(distinct, 3L)
    expect_named(distinct, channels)

    ## default argument vector resolves to "ensemble"
    default <- validate_group_channels(channels, c("ensemble", "distinct"))
    expect_length(default, 1L)
})

test_that("validate_group_channels normalises custom lists", {
    channels <- c("smo2_l", "smo2_r", "hhb", "thb")

    out <- validate_group_channels(
        channels,
        list(smo2 = c("smo2_l", "smo2_r"), "hhb")
    )

    ## unnamed groups keyed by first member; omitted channels appended
    ## as their own distinct groups
    expect_named(out, c("smo2", "hhb", "thb"))
    expect_equal(out$smo2, c("smo2_l", "smo2_r"))
    expect_equal(out$hhb, "hhb")
    expect_equal(out$thb, "thb")
})

test_that("validate_group_channels aborts on invalid groups", {
    channels <- c("smo2", "o2hb")

    ## unknown members
    expect_error(
        validate_group_channels(channels, list(c("smo2", "typo"))),
        "typo.*not recognised"
    )

    ## duplicated members across groups
    expect_error(
        validate_group_channels(channels, list("smo2", c("smo2", "o2hb"))),
        "more than one"
    )

    ## invalid shortcut string
    expect_error(
        validate_group_channels(channels, "wrong"),
        "one of.*ensemble.*distinct.*wrong"
    )
})

test_that("validate_group_channels reports errors from the caller's `env`", {
    channels <- c("smo2", "o2hb")

    ## the abort is attributed to the user-facing function via `env`
    caller <- function() {
        validate_group_channels(channels, list(c("smo2", "typo")))
    }
    err <- expect_error(caller(), "typo.*not recognised")
    expect_equal(rlang::call_name(conditionCall(err)), "caller")
})
