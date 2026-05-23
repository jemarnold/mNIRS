## example data ==========================================
make_df <- function() {
    data.frame(
        time  = 1:4,
        oxy   = c(2, 5, 6, 7),
        deoxy = c(1, 2, 4, 5),
        total = c(3, 7, 10, 12)
    )
}


test_that("returns mnirs tibble with same dimensions and names", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_s3_class(result, "tbl_df")
    expect_equal(dim(result), dim(df))
    expect_named(result, names(df))
})


test_that("errors when fewer than two channels specified", {
    df <- make_df()
    expect_error(
        correct_blood_volume(df, oxy_channel = "oxy"),
        "At least two of"
    )
    expect_error(
        correct_blood_volume(df),
        "At least two of"
    )
})


test_that("errors when channel name not found in data", {
    df <- make_df()
    expect_error(
        correct_blood_volume(
            df,
            oxy_channel   = "oxy",
            deoxy_channel = "bad_col"
        ),
        "not found.*case-sensitive"
    )
})


test_that("verbose returns message", {
    df <- make_df()
    expect_message(
        result <- correct_blood_volume(
            df,
            oxy_channel   = "oxy",
            deoxy_channel = "deoxy",
            total_channel = "total",
            verbose       = TRUE
        ),
        "channels have been corrected"
    )
})


test_that("first sample is zero (cumsum starts at 0)", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$oxy[1L],   0)
    expect_equal(result$deoxy[1L], 0)
    expect_equal(result$total[1L], 0)
})


test_that("corrected total is zero at every sample when supplied", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$total, rep(0, nrow(df)))
})


test_that("corrected oxy + corrected deoxy is zero at every sample", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$oxy + result$deoxy, rep(0, nrow(df)))
})


test_that("corrected channels match hand-calculated values", {
    ## input: oxy = (2, 5, 6, 7), deoxy = (1, 2, 4, 5), total = (3, 7, 10, 12)
    ## min = 1; shift = max(0, -1) + eps = eps  (negligible)
    ## beta[-1] = oxy[-1] / total[-1] = (5/7, 6/10, 7/12)
    ## diff(oxy)   = (3, 1, 1)
    ## diff(deoxy) = (1, 2, 1)
    ## diff(total) = (4, 3, 2)
    ##
    ## corrected oxy increments = diff(oxy) - beta[-1] * diff(total):
    ##   3 - (5/7)*4   =  1/7
    ##   1 - (6/10)*3  = -4/5
    ##   1 - (7/12)*2  = -1/6
    ## corrected oxy = cumsum(c(0, 1/7, -4/5, -1/6))
    ##              = (0, 1/7, -23/35, -173/210)
    ##
    ## corrected deoxy increments = diff(deoxy) - (1 - beta[-1]) * diff(total):
    ##   1 - (2/7)*4   = -1/7
    ##   2 - (4/10)*3  =  4/5
    ##   1 - (5/12)*2  =  1/6
    ## corrected deoxy = cumsum(c(0, -1/7, 4/5, 1/6))
    ##                = (0, -1/7, 23/35, 173/210)
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$oxy,   c(0, 1/7, -23/35, -173/210))
    expect_equal(result$deoxy, c(0, -1/7,  23/35,  173/210))
})


test_that("constant ratio across channels gives all-zero correction", {
    ## oxy:deoxy:total ratio constant; each diff exactly cancelled
    ## oxy = (0, 2, 4, 6), deoxy = (0, 1, 2, 3), total = (0, 3, 6, 9)
    ## beta[-1] = (2/3, 2/3, 2/3)
    ## corrected oxy   incr = 2 - (2/3)*3 = 0
    ## corrected deoxy incr = 1 - (1/3)*3 = 0
    df <- data.frame(
        oxy   = c(0, 2, 4, 6),
        deoxy = c(0, 1, 2, 3),
        total = c(0, 3, 6, 9)
    )
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$oxy,   rep(0, nrow(df)))
    expect_equal(result$deoxy, rep(0, nrow(df)))
    expect_equal(result$total, rep(0, nrow(df)))
})


test_that("constant signals produce zero-valued corrections", {
    ## diff of any constant series is zero, so corrections collapse to zero
    df <- data.frame(
        oxy   = rep(10, 5L),
        deoxy = rep(5,  5L),
        total = rep(15, 5L)
    )
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(result$oxy,   rep(0, nrow(df)))
    expect_equal(result$deoxy, rep(0, nrow(df)))
    expect_equal(result$total, rep(0, nrow(df)))
})


test_that("ensemble shift handles negative values", {
    ## oxy = (-1, 0, 1), deoxy = (-2, -1, 0), total = (-3, -1, 1)
    ## min = -3; shift = 3 + eps
    ## shifted oxy = (2, 3, 4), deoxy = (1, 2, 3), total = (0, 2, 4) (+ eps)
    ## beta[-1] = oxy[-1] / total[-1] = (3/2, 1)
    ## diff(oxy) = (1, 1), diff(deoxy) = (1, 1), diff(total) = (2, 2)
    ## corrected oxy   incr = 1 - (3/2)*2 = -2;  1 - 1*2 = -1
    ## corrected deoxy incr = 1 - (-1/2)*2 = 2;  1 - 0*2 =  1
    ## corrected oxy   = (0, -2, -3)
    ## corrected deoxy = (0,  2,  3)
    df <- data.frame(
        time  = 1:3,
        oxy   = c(-1,  0, 1),
        deoxy = c(-2, -1, 0),
        total = c(-3, -1, 1)
    )
    # df <- create_mnirs_data(
    #     df,
    #     nirs_channels = c(oxy, deoxy, total),
    #     time_channel = time
    # )
    # plot(df)
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    # plot(result)

    expect_equal(result$oxy,   c(0, -2, -3))
    expect_equal(result$deoxy, c(0,  2,  3))
    expect_equal(result$total, c(0,  0,  0))
})


test_that("derives total from oxy + deoxy when total not specified", {
    ## same as hand-calculated case but total derived internally
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        verbose       = FALSE
    )
    ## unspecified total column passes through unchanged
    expect_equal(result$total, df$total)
    ## corrected oxy and deoxy match hand-calculated values
    expect_equal(result$oxy,   c(0,  1/7, -23/35, -173/210))
    expect_equal(result$deoxy, c(0, -1/7,  23/35,  173/210))
})


test_that("derives oxy from total - deoxy when oxy not specified", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    ## unspecified oxy column passes through unchanged
    expect_equal(result$oxy, df$oxy)
    ## corrected deoxy and total match hand-calculated values
    expect_equal(result$deoxy, c(0, -1/7, 23/35, 173/210))
    expect_equal(result$total, rep(0, nrow(df)))
})


test_that("derives deoxy from total - oxy when deoxy not specified", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        total_channel = "total",
        verbose       = FALSE
    )
    ## unspecified deoxy column passes through unchanged
    expect_equal(result$deoxy, df$deoxy)
    ## corrected oxy and total match hand-calculated values
    expect_equal(result$oxy,   c(0, 1/7, -23/35, -173/210))
    expect_equal(result$total, rep(0, nrow(df)))
})


test_that("non-NIRS columns pass through unchanged", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        deoxy_channel = "deoxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_identical(result$time, df$time)
})


test_that("metadata nirs_channels reflects specified channels only", {
    df <- make_df()
    result <- correct_blood_volume(
        df,
        oxy_channel   = "oxy",
        total_channel = "total",
        verbose       = FALSE
    )
    expect_equal(attr(result, "nirs_channels"), c("oxy", "total"))
})

## integration ==============================================
test_that("correct_blood_volume preserves metadata", {
    data <- read_mnirs(
        file_path = example_mnirs("artinis"),
        nirs_channels = c(o2hb = 2, hhb = 3),
        time_channel = c(sample = 1),
        verbose = FALSE
    ) |>
        ## derive thb
        dplyr::mutate(thb = o2hb + hhb) |>
        ## update nirs_channels with thb
        create_mnirs_data(nirs_channels = c(thb, o2hb, hhb))

    expect_message(
        result <- correct_blood_volume(
            data,
            oxy_channel = "o2hb",
            deoxy_channel = "hhb",
            total_channel = "thb",
        ),
        "o2hb.*hhb.*thb.*corrected"
    )

    expect_equal(attr(result, "nirs_device"), attr(data, "nirs_device"))
    expect_setequal(attr(result, "nirs_channels"), attr(data, "nirs_channels"))
    expect_equal(attr(result, "time_channel"), attr(data, "time_channel"))
    expect_equal(attr(result, "sample_rate"), attr(data, "sample_rate"))

    ## overwrite nirs_channels
    result <- correct_blood_volume(
        data,
        oxy_channel = "o2hb",
        deoxy_channel = "hhb",
        total_channel = NULL,
        verbose = FALSE
    )

    expect_equal(
        setdiff(attr(data, "nirs_channels"), attr(result, "nirs_channels")),
        "thb"
    )
})
