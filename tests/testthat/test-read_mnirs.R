## example_mnirs() ====================================================
test_that("example_mnirs() returns all files when file = NULL", {
    files <- example_mnirs()
    expect_type(files, "character")
    expect_true(length(files) > 0)
})

test_that("example_mnirs() returns valid path for exact match", {
    # Assumes at least one file exists in extdata
    all_files <- example_mnirs()
    skip_if(length(all_files) == 0, "No example files available")

    path <- example_mnirs(all_files[1])
    expect_type(path, "character")
    expect_true(file.exists(path))
    expect_match(path, all_files[1], fixed = TRUE)
})

test_that("example_mnirs() returns valid path for partial match", {
    skip_if_not(
        any(grepl("moxy_ramp", example_mnirs(), fixed = TRUE)),
        "moxy_ramp.xlsx not available"
    )

    path <- example_mnirs("moxy_ramp")
    expect_true(file.exists(path))
    expect_match(path, "moxy_ramp", fixed = TRUE)
})

test_that("example_mnirs() errors on multiple partial matches", {
    skip_if_not(
        sum(grepl("moxy", example_mnirs(), fixed = TRUE)) > 1,
        "Multiple moxy files not available"
    )

    expect_error(example_mnirs("moxy"), "Multiple files match")
})

test_that("example_mnirs() errors on non-existent file", {
    expect_error(
        example_mnirs("nonexistent_file_xyz"),
        "'arg' should be one of"
    )
})

test_that("example_mnirs() does not show files with `~`", {
    # Create test directory structure
    test_dir <- file.path(tempdir(), "epl_test")
    dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)

    # Create test files
    file.create(file.path(test_dir, "data.csv"))
    file.create(file.path(test_dir, "~temp.csv"))
    file.create(file.path(test_dir, "results.xlsx"))

    result <- list.files(test_dir)
    expect_length(result, 3)
    expect_true(any(grepl("^~", result)))

    result <- list.files(test_dir, pattern = "^[^~]")

    expect_length(result, 2)
    expect_false(any(grepl("^~", result)))
    expect_true(all(c("data.csv", "results.xlsx") %in% result))
})


## read_file() ==================================================

test_that("read_file() reads moxy (perfpro) xlsx files correctly", {
    file_path <- example_mnirs("moxy_ramp")
    skip_if(!grepl("\\.xls(x)?$", file_path, ignore.case = TRUE))

    result <- read_file(file_path)

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 5)
    expect_gt(nrow(result), 0)
    expect_all_true(unlist(lapply(result, is.character)))
})

test_that("read_file() reads moxy .csv files correctly", {
    file_path <- example_mnirs("moxy_intervals")
    skip_if(!grepl("\\.csv$", file_path, ignore.case = TRUE))

    result <- read_file(file_path)

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 7)
    expect_gt(nrow(result), 0)
    expect_all_true(unlist(lapply(result, is.character)))
})

test_that("read_file() reads explicit comma train.red csv correctly", {
    file_path <- example_mnirs("train.red")
    skip_if(!grepl("\\.csv$", file_path, ignore.case = TRUE))

    result <- read_file(file_path)

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 12)
    expect_gt(nrow(result), 0)
    expect_all_true(unlist(lapply(result, is.character)))
})

test_that("read_file() reads train.red files correctly", {
    file_path <- test_path("testdata/train.red-mre.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    result <- read_file(file_path)

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 23)
    expect_gt(nrow(result), 0)
    expect_all_true(unlist(lapply(result, is.character)))
})

test_that("read_file() reads vo2master files correctly", {
    file_path <- test_path("testdata/vo2master.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    result <- read_file(file_path)

    expect_s3_class(result, "data.frame")
    expect_gt(ncol(result), 12)
    expect_gt(nrow(result), 0)
    expect_all_true(unlist(lapply(result, is.character)))
})

test_that("read_file() errors", {
    expect_error(
        read_file("nonexistent_file.xlsx"),
        "File not found"
    )

    temp_file <- tempfile(fileext = ".docx")
    writeLines("test", temp_file)
    on.exit(unlink(temp_file))

    expect_error(
        read_file(temp_file),
        "Unsupported file type"
    )
})

test_that("read_file() errors when Excel file cannot be opened", {
    file_path <- example_mnirs("moxy_ramp")
    skip_if(!grepl("\\.xls(x)?$", file_path, ignore.case = TRUE))

    local_mocked_bindings(
        read_excel = function(...) stop("cannot be opened"),
        .package = "readxl"
    )

    expect_error(
        read_file(file_path),
        "File cannot be opened"
    )
})

test_that("read_file() re-throws other Excel read errors", {
    file_path <- example_mnirs("moxy_ramp")
    skip_if(!grepl("\\.xls(x)?$", file_path, ignore.case = TRUE))

    local_mocked_bindings(
        read_excel = function(...) stop("some other error"),
        .package = "readxl"
    )

    expect_error(
        read_file(file_path),
        "some other error"
    )
})


## detect_mnirs_device() ===============================================
test_that("detect_mnirs_device works on example files", {
    ## xlsx files
    expect_equal(
        read_file(example_mnirs("moxy_ramp")) |>
            detect_mnirs_device(),
        list(
            nirs_device = "Moxy",
            header_row = 6
        )
    )

    expect_equal(
        read_file(example_mnirs("artinis_intervals")) |>
            detect_mnirs_device(),
        list(
            nirs_device = "Artinis",
            header_row = 38
        )
    )

    ## csv files
    expect_equal(
        read_file(example_mnirs("train.red")) |>
            detect_mnirs_device(),
        list(
            nirs_device = "Train.Red",
            header_row = 41
        )
    )
})

test_that("detect_mnirs_device works on internal example files", {
    skip_on_ci()
    skip_on_covr()
    file_path <- test_path("testdata/train.red-mre.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    expect_equal(
        read_file(file_path) |>
            detect_mnirs_device(),
        list(
            nirs_device = "Train.Red",
            header_row = 521
        )
    )

    file_path <- test_path("testdata/vo2master.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    expect_equal(
        read_file(file_path) |>
            detect_mnirs_device(),
        list(
            nirs_device = "VO2master",
            header_row = 1
        )
    )
})

test_that("detect_mnirs_device() returns NULL when no match", {
    data <- data.frame(
        V1 = c("Unknown", "device", "data"),
        V2 = c("header", "col1", "val1")
    )

    expect_equal(
        detect_mnirs_device(data),
        list(nirs_device = NULL, header_row = 1)
    )
})

test_that("detect_mnirs_device() requires 'oxysoft' for Artinis match", {
    ## bare numeric rows match Artinis regex, but no 'oxysoft' above
    data <- data.frame(
        V1 = c("some", "1", "2"),
        V2 = c("header", "2", "3"),
        V3 = c("info", "3", "4")
    )

    expect_equal(
        detect_mnirs_device(data),
        list(nirs_device = NULL, header_row = 1)
    )

    ## 'oxysoft' present above numeric row — case-insensitive
    data <- data.frame(
        V1 = c("Exported from OxySoft", "more", "1", "2"),
        V2 = c("v1", "info", "2", "3"),
        V3 = c("", "", "3", "4")
    )

    expect_equal(
        detect_mnirs_device(data),
        list(nirs_device = "Artinis", header_row = 3)
    )
})


## resolve_channels() ==================================================
## helpers: device list as from detect_mnirs_device(); user list as in read_mnirs()
test_device <- function(nirs_device = NULL, header_row = 1L) {
    list(nirs_device = nirs_device, header_row = header_row)
}
test_user <- function(nirs = NULL, time = NULL, event = NULL) {
    lapply(list(time = time, event = event, nirs = nirs), name_channels)
}

test_that("resolve_channels() returns user channels when provided", {
    raw <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw,
        test_device("Moxy"),
        test_user(nirs = c(smo2 = "SmO2 Live"), time = c(time = "hh:mm:ss")),
        verbose = FALSE
    )

    expect_named(result, c("time", "event", "extra", "nirs"))
    expect_equal(result$nirs, c(smo2 = "SmO2 Live"))
    expect_equal(result$time, c(time = "hh:mm:ss"))
    expect_null(result$event)
    expect_null(result$extra)

    ## NULL device
    result <- resolve_channels(
        raw,
        test_device(NULL),
        test_user(nirs = "O2Hb", time = "Time"),
        verbose = FALSE
    )

    expect_equal(result$nirs, c(O2Hb = "O2Hb"))
    expect_equal(result$time, c(Time = "Time"))
})

test_that("resolve_channels() applies device default time with user nirs", {
    raw <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw,
        test_device("Moxy"),
        test_user(nirs = c(smo2 = "SmO2 Live")),
        verbose = FALSE
    )

    expect_equal(result$time, device_patterns$Moxy$time_channel)

    ## user time_channel overrides device default
    result <- resolve_channels(
        raw,
        test_device("Moxy"),
        test_user(time = c(time = "custom_time")),
        verbose = FALSE
    )

    expect_equal(result$nirs, "SmO2 Live")
    expect_equal(result$time, c(time = "custom_time"))
})

test_that("resolve_channels() detects known channels for device", {
    raw <- data.frame(
        V1 = c("meta", "hh:mm:ss", "00:00:01"),
        V2 = c("meta", "SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw,
        test_device("Moxy", 2L),
        test_user(),
        verbose = FALSE
    )

    expect_equal(result$nirs, "SmO2 Live")
    expect_equal(result$time, device_patterns$Moxy$time_channel)
})

test_that("resolve_channels() detects multiple SmO2 channels", {
    raw <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("SmO2 (1)", "55"),
        V3 = c("SmO2 (2)", "60"),
        V4 = c("HR", "120"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw, test_device("PerfPro"), test_user(), verbose = FALSE
    )

    expect_equal(result$nirs, c("SmO2 (1)", "SmO2 (2)"))
})

test_that("resolve_channels() matches SmO2 case-insensitively", {
    raw <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("smo2 raw", "55"),
        V3 = c("SMO2_LIVE", "60"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw, test_device("Moxy"), test_user(), verbose = FALSE
    )

    expect_equal(result$nirs, c("smo2 raw", "SMO2_LIVE"))
})

test_that("resolve_channels() drops redundant unfiltered/Averaged SmO2", {
    raw <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("SmO2", "55"),
        V3 = c("SmO2 unfiltered", "60"),
        V4 = c("SmO2 Averaged", "60"),
        V5 = c("Lap/Event", "1"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw, test_device("Train.Red"), test_user(), verbose = FALSE
    )

    expect_equal(result$nirs, "SmO2")
    expect_equal(result$event, "Lap/Event")

    ## device event default dropped when absent from header
    expect_null(resolve_channels(
        raw[-5L], test_device("Train.Red"), test_user(), verbose = FALSE
    )$event)
})

test_that("resolve_channels() parses Artinis legend block", {
    ## Artinis names channels in the legend block, not the header row
    raw <- data.frame(
        V1 = c("OxySoft export of", "Legend", "Column", "1", "2", "3", "4", NA, "1"),
        V2 = c(NA, NA, "Trace (Measurement)", "(Sample number)",
               "VL O2Hb", "VL HHb", "(Event)", NA, "2"),
        V3 = c(rep(NA, 8), "3"),
        V4 = c(rep(NA, 8), "4"),
        V5 = NA_character_,
        stringsAsFactors = FALSE
    )
    device <- test_device("Artinis", 9L)

    result <- resolve_channels(raw, device, test_user(), verbose = FALSE)

    expect_equal(result$nirs, c(vl_o2hb = "2", vl_hhb = "3"))
    expect_equal(result$time, c(sample = "1"))
    expect_equal(result$event, c(event = "4"))
    ## unnumbered trailing column mapped to "labels"
    expect_equal(result$extra, c(labels = "col_5"))

    ## `labels` accompanies an auto-detected event only
    result <- resolve_channels(
        raw, device, test_user(event = c(ev = "4")), verbose = FALSE
    )
    expect_equal(result$event, c(ev = "4"))
    expect_null(result$extra)

    ## unless keep_all
    result <- resolve_channels(
        raw, device, test_user(event = c(ev = "4")), keep_all = TRUE,
        verbose = FALSE
    )
    expect_equal(result$extra, c(labels = "col_5"))

    ## legend entries for user-claimed columns are dropped
    result <- resolve_channels(
        raw, device, test_user(time = c(t = "4")), verbose = FALSE
    )
    expect_equal(result$time, c(t = "4"))
    expect_null(result$event)

    ## `event_channel = "labels"` aliases the unnumbered label column
    result <- resolve_channels(
        raw, device, test_user(event = c(event = "labels")), verbose = FALSE
    )
    expect_equal(result$event, c(event = "col_5"))
    expect_null(result$extra)
})

test_that("resolve_channels() Artinis falls back without legend", {
    raw <- data.frame(
        V1 = c("1", "2", "3"),
        V2 = c("4", "5", "6"),
        stringsAsFactors = FALSE
    )
    device <- test_device("Artinis", 2L)

    ## no legend: time defaults to sample = "1", nirs cannot be resolved
    expect_error(
        resolve_channels(raw, device, test_user(), verbose = FALSE),
        "cannot be determined"
    )

    result <- resolve_channels(
        raw, device, test_user(nirs = c(O2Hb = 2)), verbose = FALSE
    )

    expect_equal(result$nirs, c(O2Hb = "2"))
    expect_equal(result$time, c(sample = "1"))
})

test_that("resolve_channels() detects SmO2 for unknown device", {
    raw <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("SmO2", "55"),
        stringsAsFactors = FALSE
    )

    result <- resolve_channels(
        raw, test_device(NULL), test_user(), verbose = FALSE
    )

    expect_equal(result$nirs, "SmO2")
    expect_null(result$time)
})

test_that("resolve_channels() detects known channels for PerfPro", {
    file_path <- test_path("testdata/perfpro-mre.xlsx")
    skip_if_not(file.exists(file_path), "testdata not available")

    raw <- read_file(file_path)
    device <- detect_mnirs_device(raw)

    expect_match(device$nirs_device, "PerfPro")
    expect_equal(device$header_row, 3)

    result <- resolve_channels(raw, device, test_user(), verbose = FALSE)

    expect_equal(result$nirs, c("SmO2 (1614)", "SmO2 (1615)"))
    expect_equal(result$time, "Time")
})

test_that("resolve_channels() errors when no SmO2 columns found", {
    raw <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("HR", "120"),
        stringsAsFactors = FALSE
    )

    ## NULL device + no SmO2 columns
    expect_error(
        resolve_channels(raw, test_device(NULL), test_user(), verbose = FALSE),
        "cannot be determined"
    )

    ## known device + no SmO2 columns
    expect_error(
        resolve_channels(raw, test_device("Moxy"), test_user(), verbose = FALSE),
        "cannot be determined"
    )
})

test_that("resolve_channels() verbose messages for detection", {
    raw <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    expect_message(
        resolve_channels(raw, test_device("Moxy"), test_user(), verbose = TRUE),
        "Moxy.*detected"
    )

    ## unknown device labelled "Unknown"
    expect_message(
        resolve_channels(raw, test_device(NULL), test_user(), verbose = TRUE),
        "Unknown.*detected"
    )

    ## no message when user provides channels
    expect_no_message(
        resolve_channels(
            raw, test_device("Moxy"), test_user(nirs = "SmO2 Live"),
            verbose = TRUE
        )
    )
})


## parse_oxysoft_legend() ==============================================
test_that("clean_channel_names() collapses non-alphanumerics and lowercases", {
    expect_equal(clean_channel_names("VL O2Hb"), "vl_o2hb")
    expect_equal(clean_channel_names("Rx1 - Tx1 tHb"), "rx1_tx1_thb")
    expect_equal(clean_channel_names("(TSI)"), "tsi")
    expect_equal(clean_channel_names("  a  b  "), "a_b")
    expect_equal(clean_channel_names("***"), "")
})

test_that("parse_oxysoft_legend() parses full legend with labels column", {
    data <- data.frame(
        V1 = c("Legend", "Column", "1", "2", "3", "4", "5", NA, "1"),
        V2 = c(NA, "Trace (Measurement)", "(Sample number)", "Rx1 - Tx1 tHb",
               "Rx1 - Tx1 O2Hb", "(TSI)", "(Event)", NA, "2"),
        V3 = c(rep(NA, 8), "3"),
        V4 = c(rep(NA, 8), "4"),
        V5 = c(rep(NA, 8), "5"),
        V6 = NA_character_,
        stringsAsFactors = FALSE
    )

    result <- parse_oxysoft_legend(data, header_row = 9L)

    expect_named(result, c("time", "event", "extra", "nirs"))
    expect_equal(result$time, c(sample = "1"))
    expect_equal(
        result$nirs,
        c(rx1_tx1_thb = "2", rx1_tx1_o2hb = "3")
    )
    expect_equal(result$event, c(event = "5"))
    ## unknown parenthesised trace kept, unnumbered trailing col -> labels
    expect_equal(result$extra, c(tsi = "4", labels = "col_6"))
})

test_that("parse_oxysoft_legend() omits labels without trailing column", {
    data <- data.frame(
        V1 = c("Legend", "Column", "1", "2", "3", NA, "1"),
        V2 = c(NA, "Trace (Measurement)", "(Sample number)", "O2Hb",
               "(Event)", NA, "2"),
        V3 = c(rep(NA, 6), "3"),
        stringsAsFactors = FALSE
    )

    result <- parse_oxysoft_legend(data, header_row = 7L)

    expect_equal(result$nirs, c(o2hb = "2"))
    expect_equal(result$event, c(event = "3"))
    expect_null(result$extra)
})

test_that("parse_oxysoft_legend() returns NULL when legend missing or malformed", {
    ## no legend row
    data <- data.frame(
        V1 = c("meta", "1"),
        V2 = c("meta", "2"),
        stringsAsFactors = FALSE
    )
    expect_null(parse_oxysoft_legend(data, header_row = 2L))

    ## empty trace name
    data <- data.frame(
        V1 = c("Legend", "Column", "1", "2", "1"),
        V2 = c(NA, "Trace (Measurement)", "(Sample number)", NA, "2"),
        stringsAsFactors = FALSE
    )
    expect_null(parse_oxysoft_legend(data, header_row = 5L))

    ## legend ids not present in the header row
    data <- data.frame(
        V1 = c("Legend", "Column", "1", "2", "A"),
        V2 = c(NA, "Trace (Measurement)", "(Sample number)", "O2Hb", "B"),
        stringsAsFactors = FALSE
    )
    expect_null(parse_oxysoft_legend(data, header_row = 5L))
})


## find_header_row() ===================================================
test_that("find_header_row() finds the row containing all channels", {
    raw <- data.frame(
        V1 = c("meta1", "meta2", "O2Hb", "10", "20"),
        V2 = c("meta1", "meta2", "HHb", "5", "15"),
        V3 = c("meta1", "meta2", "Time", "0.1", "0.2"),
        stringsAsFactors = FALSE
    )

    expect_equal(find_header_row(raw, c("O2Hb", "HHb")), 3L)
    ## `start` row tried first
    expect_equal(find_header_row(raw, c("O2Hb", "HHb"), start = 3L), 3L)
    ## numeric channel ids coerced to character
    expect_equal(find_header_row(raw, c(oxy = "O2Hb"), start = 1L), 3L)
})

test_that("find_header_row() matches renamed duplicate headers", {
    raw <- data.frame(
        V1 = c("meta", "SmO2", "10"),
        V2 = c("meta", "SmO2", "5"),
        V3 = c("meta", "Time", "0.1"),
        stringsAsFactors = FALSE
    )

    expect_equal(find_header_row(raw, c(left = "SmO2", right = "SmO2_1")), 2L)
})

test_that("find_header_row() errors when channels not found", {
    raw <- data.frame(
        V1 = c("header", "WrongChannel", "10"),
        V2 = c("header", "Time", "0.1"),
        stringsAsFactors = FALSE
    )

    expect_error(find_header_row(raw, "O2Hb"), "Channel names not detected")
})

test_that("find_header_row() is case sensitive", {
    raw <- data.frame(
        V1 = c("meta", "o2hb", "10"),
        V2 = c("meta", "time", "0.1"),
        stringsAsFactors = FALSE
    )

    expect_error(find_header_row(raw, "O2Hb"), "case sensitive")
})


## detect_time_channel() ==============================================
test_that("detect_time_channel finds time column by name", {
    df <- data.frame(time = 1:5, value = 6:10)
    expect_equal(detect_time_channel(df, verbose = FALSE), "time")

    df <- data.frame(Time = 1:5, value = 6:10)
    expect_equal(detect_time_channel(df, verbose = FALSE), "Time")

    df <- tibble("hh:mm:ss" = 1:5, value = 6:10)
    expect_equal(detect_time_channel(df, verbose = FALSE), "hh:mm:ss")

    df <- data.frame("hms" = 1:5, value = 6:10)
    expect_equal(detect_time_channel(df, verbose = FALSE), "hms")

    df <- data.frame(duration = 1:5, value = 6:10)
    expect_equal(detect_time_channel(df, verbose = FALSE), "duration")
})

test_that("detect_time_channel finds character time format", {
    df <- data.frame(
        value = 1:5,
        string_col = c(
            "12:30:45",
            "12:30:46",
            "12:30:47",
            "12:30:48",
            "12:30:49"
        )
    )
    expect_equal(detect_time_channel(df, verbose = FALSE), "string_col")

    # Test with H:MM format
    df <- data.frame(
        value = 1:5,
        string_col = c("1:30", "1:31", "1:32", "1:33", "1:34")
    )
    expect_equal(detect_time_channel(df, verbose = FALSE), "string_col")
})

test_that("detect_time_channel handles NA and empty values in character column", {
    df <- data.frame(
        value = 1:5,
        time_str = c(NA, "", "12:30:46", "12:30:47", "12:30:48")
    )
    expect_equal(detect_time_channel(df, verbose = FALSE), "time_str")
})

test_that("detect_time_channel errors when no time column found", {
    df <- data.frame(x = 1:5, y = 6:10)
    expect_error(
        detect_time_channel(df, verbose = FALSE),
        "time_channel.*not detected"
    )
})

test_that("detect_time_channel verbose messages work", {
    df <- data.frame(time = 1:5, value = 6:10)
    expect_message(
        detect_time_channel(df, verbose = TRUE),
        "Detected.*time_channel"
    )
})


## rename_duplicates() ================================================
test_that("rename_duplicates() handles duplicate strings", {
    x <- c("O2Hb", "HHb", "O2Hb", "Time")
    result <- rename_duplicates(x)

    expect_equal(result, c("O2Hb", "HHb", "O2Hb_1", "Time"))
    expect_true(all(!duplicated(result)))
})

test_that("rename_duplicates() handles empty strings", {
    x <- c("O2Hb", "", "HHb", "")
    result <- rename_duplicates(x)

    expect_equal(result, c("O2Hb", "col_2", "HHb", "col_4"))
})

test_that("rename_duplicates() handles NA values", {
    x <- c("O2Hb", NA, "HHb", NA)
    result <- rename_duplicates(x)

    expect_equal(result, c("O2Hb", "col_2", "HHb", "col_4"))
})

test_that("rename_duplicates() handles NULL input", {
    expect_null(rename_duplicates(NULL))
})

test_that("rename_duplicates() handles all unique values", {
    x <- c("O2Hb", "HHb", "Time")
    result <- rename_duplicates(x)

    expect_equal(result, x)
})

test_that("rename_duplicates() handles multiple duplicates", {
    x <- c("O2Hb", "O2Hb", "O2Hb")
    result <- rename_duplicates(x)

    expect_equal(result, c("O2Hb", "O2Hb_1", "O2Hb_2"))
})

## name_channels() ====================================================
test_that("name_channels() names unnamed vector", {
    x <- c("O2Hb", "HHb", "Time")
    result <- name_channels(x)

    expect_equal(names(result), c("O2Hb", "HHb", "Time"))
    expect_equal(as.character(result), x)
})

test_that("name_channels() preserves existing names", {
    x <- c(oxy = "O2Hb", deoxy = "HHb", time = "Time")
    result <- name_channels(x)

    expect_equal(names(result), c("oxy", "deoxy", "time"))
    expect_equal(as.character(result), c("O2Hb", "HHb", "Time"))
})

test_that("name_channels() fills in missing names", {
    x <- c(oxy = "O2Hb", "HHb", time = "Time")
    result <- name_channels(x)

    expect_equal(names(result), c("oxy", "HHb", "time"))
})

test_that("name_channels() handles NA names", {
    x <- c("O2Hb", "HHb")
    names(x) <- c(NA, "deoxy")
    result <- name_channels(x)

    expect_equal(names(result), c("O2Hb", "deoxy"))
})

test_that("name_channels() handles all empty names", {
    x <- c("O2Hb", "HHb", "Time")
    names(x) <- c("", "", "")
    result <- name_channels(x)

    expect_equal(names(result), c("O2Hb", "HHb", "Time"))
})

test_that("name_channels() coerces numeric ids to character and passes NULL", {
    expect_equal(name_channels(c(HHb = 2, 3)), c(HHb = "2", "3" = "3"))
    expect_null(name_channels(NULL))
})


## match_channels() ===================================================
test_that("match_channels() keeps list order and names channels", {
    result <- match_channels(
        list(time = c(time = "Time"), nirs = c(oxy = "O2Hb", deoxy = "HHb")),
        data_names = c("O2Hb", "HHb", "Time"),
        verbose = FALSE
    )

    expect_named(result, c("time", "nirs"))
    expect_equal(result$time, c(time = "Time"))
    expect_equal(result$nirs, c(oxy = "O2Hb", deoxy = "HHb"))
})

test_that("match_channels() names unnamed channels and drops NULL", {
    result <- match_channels(
        list(time = "Time", event = NULL, nirs = "O2Hb"),
        data_names = c("O2Hb", "Time"),
        verbose = FALSE
    )

    expect_named(result, c("time", "nirs"))
    expect_equal(result$time, c(Time = "Time"))
    expect_equal(result$nirs, c(O2Hb = "O2Hb"))
})

test_that("match_channels() handles un-renamed duplicate channels", {
    expect_warning(
        result <- match_channels(
            list(time = "Time", nirs = c("O2Hb", "O2Hb")),
            data_names = c("O2Hb", "O2Hb_1", "Time"),
            verbose = TRUE
        ),
        "Duplicate channel names"
    )

    expect_equal(result$nirs, c(O2Hb = "O2Hb", O2Hb_1 = "O2Hb_1"))
})

test_that("match_channels() handles renamed duplicate data columns", {
    expect_no_warning(
        result <- match_channels(
            list(time = "Time", nirs = c(oxy1 = "O2Hb", oxy2 = "O2Hb")),
            data_names = c("O2Hb", "O2Hb_1", "Time"),
            verbose = TRUE
        )
    )

    expect_equal(result$nirs, c(oxy1 = "O2Hb", oxy2 = "O2Hb_1"))
})

test_that("match_channels() makes names unique across roles", {
    expect_warning(
        result <- match_channels(
            list(time = c(smo2 = "Time"), nirs = c(smo2 = "SmO2")),
            data_names = c("SmO2", "Time"),
            verbose = TRUE
        ),
        "smo2 = smo2_1"
    )

    expect_equal(result$time, c(smo2 = "Time"))
    expect_equal(result$nirs, c(smo2_1 = "SmO2"))
})

test_that("match_channels() errors when channel not found", {
    expect_error(
        match_channels(
            list(time = "Time", nirs = "HHb"),
            data_names = c("O2Hb", "Time")
        ),
        "Channel names not detected"
    )
})

test_that("match_channels() suppresses warnings with verbose", {
    expect_silent(
        match_channels(
            list(time = "Time", nirs = c(o2hb = "O2Hb", o2hb = "O2Hb")),
            data_names = c("O2Hb", "O2Hb_1", "Time"),
            verbose = FALSE
        )
    )
})

test_that("read_mnirs() selects, orders, and renames columns", {
    file_path <- tempfile(fileext = ".csv")
    on.exit(unlink(file_path))
    writeLines(
        c(
            "O2Hb,HHb,Time,custom,Extra",
            "10,5,0.1,x,a",
            "20,15,0.2,y,b"
        ),
        file_path
    )

    ## default keep_all = FALSE: time first, then nirs
    df <- read_mnirs(
        file_path,
        nirs_channels = c(oxy = "O2Hb", deoxy = "HHb"),
        time_channel = c(time = "Time"),
        verbose = FALSE
    )

    expect_equal(names(df), c("time", "oxy", "deoxy"))
    expect_equal(attr(df, "nirs_channels"), c("oxy", "deoxy"))
    expect_equal(attr(df, "time_channel"), "time")

    ## keep_all: custom names take priority over clashing data names
    df <- read_mnirs(
        file_path,
        nirs_channels = c(custom = "O2Hb"),
        time_channel = "Time",
        keep_all = TRUE,
        verbose = FALSE
    )

    expect_equal(names(df), c("Time", "custom", "HHb", "custom_1", "Extra"))
    expect_equal(df$custom, c(10, 20))
    expect_equal(df$custom_1, c("x", "y"))
})


## convert_type() ================================================
test_that("convert_type() applies unopinionated typing to data columns", {
    data <- tibble(
        time = c(1, 2, 3),
        lap = c(1.0, NA_real_, 2.0),  ## event -> integer laps
        B = c(10, 20, 30),            ## whole numbers -> double
        x = c(10.5, 11.0, 11.5),      ## fractional -> double
    )

    result <- convert_type(data, list(time = "time", event = "lap"))

    expect_type(result$time, "double")   ## time left unchanged
    expect_type(result$lap, "integer")
    expect_equal(result$lap, c(1L, NA_integer_, 2L))
    expect_type(result$B, "double")     ## unopinionated: whole -> duble
    expect_type(result$x, "double")
})


test_that("convert_type() standardises empty and 'NA' strings to NA", {
    data <- data.frame(
        time = c("1", "2", "3"),
        A = c("a", "", "b"),
        B = c("x", "NA", "y"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, list(time = "time"))

    expect_equal(result$A, c("a", NA_character_, "b"))
    expect_equal(result$B, c("x", NA_character_, "y"))
})

test_that("convert_type() standardises Inf/NaN to NA in numeric cols", {
    data <- data.frame(
        time = c("1", "2", "3", "4"),
        A = c("1.5", "Inf", "-Inf", "NaN"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, list(time = "time"))

    expect_type(result$A, "double")
    expect_equal(result$A, c(1.5, NA_real_, NA_real_, NA_real_))
})

test_that("convert_type() standardises non-finite integers to NA", {
    ## event and whole-number cols resolve to integer via type.convert
    data <- tibble(
        time = c("1", "2", "3"),
        lap = c("1", NA, "2"),
        B = c("10", "20", "30"),
    )

    result <- convert_type(data, list(time = "time", event = "lap"))

    expect_type(result$lap, "integer")
    expect_equal(result$lap, c(1L, NA_integer_, 2L))
    expect_type(result$B, "double")
})

test_that("convert_type() preserves valid numeric values", {
    data <- data.frame(
        time = c("1", "2"),
        A = c("0", "-0"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, list(time = "time"))
    expect_equal(result$A, c(0, 0))
})

test_that("convert_type() forces nirs_channels to numeric", {
    data <- tibble(
        time = c("1", "2", "3"),
        smo2 = c("55", "< 0.5", "60"),  ## qualifier in a signal col
        other = c("55", "< 0.5", "60"),  ## same values, not a signal
        label = c("Start", "Lap", "Stop")
    )

    result <- convert_type(data, list(time = "time", nirs = "smo2"))

    ## specified `nirs_channels` coerced to numeric
    expect_type(result$smo2, "double")
    expect_equal(result$smo2, c(55, NA_real_, 60))
    ## partly-numeric col type.converted without opinion
    expect_type(result$other, "character")
    expect_equal(result$other, data$other)
    ## non-signal text col unaffected
    expect_type(result$label, "character")
    expect_equal(result$label, c("Start", "Lap", "Stop"))
})

test_that("convert_type() warns per nirs channel coerced to all NA", {
    data <- tibble(
        time = c("1", "2", "3"),
        smo2 = c("a", "b", "c"),   ## non-numeric -> all NA
        hhb = c("x", "y", "z"),    ## non-numeric -> all NA
        o2hb = c("55", "z", "60"), ## partly numeric -> kept
    )

    expect_warning(
        result <- convert_type(
            data,
            list(time = "time", nirs = c("smo2", "hhb", "o2hb"))
        ),
        "smo2"
    ) |> 
        expect_warning("hhb")
    expect_equal(result$smo2, rep(NA_real_, 3))
    expect_equal(result$o2hb, c(55, NA_real_, 60))
})

test_that("convert_type() all-NA warning respects verbose = FALSE", {
    data <- tibble(
        time = c("1", "2"),
        smo2 = c("a", "b"),
    )

    expect_no_warning(
        convert_type(
            data, list(time = "time", nirs = "smo2"), verbose = FALSE
        )
    )
})


## remove_empty_rows_cols() ===========================================
test_that("remove_empty_rows_cols() removes empty rows & cols", {
    data <- data.frame(
        A = c("1", "", "3"),
        B = c("x", "", "z"),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)
    expect_equal(nrow(result), 2)
    expect_equal(result$A, c("1", "3"))

    data <- data.frame(
        A = c("1", "2"),
        B = c("", ""),
        C = c("x", "y"),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)

    expect_equal(ncol(result), 2)
    expect_equal(names(result), c("A", "C"))

    data <- data.frame(
        A = c("1", "", "3"),
        B = c("", "", ""),
        C = c("x", "", "z"),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)
    expect_equal(nrow(result), 2)
    expect_equal(result$A, c("1", "3"))
    expect_equal(ncol(result), 2)
    expect_equal(names(result), c("A", "C"))
})

test_that("remove_empty_rows_cols() handles NA values", {
    data <- data.frame(
        A = c("1", NA, "3"),
        B = c(NA, NA, NA),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)

    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), 1)
})

test_that("remove_empty_rows_cols() preserves non-empty data", {
    data <- data.frame(
        A = c("1", "2", "3"),
        B = c("x", "y", "z"),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)

    expect_equal(result, data)
})

test_that("remove_empty_rows_cols() handles all empty data", {
    data <- data.frame(
        A = c("", ""),
        B = c("", ""),
        stringsAsFactors = FALSE
    )

    result <- remove_empty_rows_cols(data)

    expect_equal(nrow(result), 0)
    expect_equal(ncol(result), 0)
})

## extract_start_timestamp() ==========================================
test_that("extract_start_timestamp() returns NULL when no timestamps present", {
    data <- data.frame(
        V1 = c("device", "sensor"),
        V2 = c("model_x", "ch1"),
        stringsAsFactors = FALSE
    )
    expect_null(extract_start_timestamp(data))
})

test_that("extract_start_timestamp() detects ISO 8601 timestamp", {
    data <- data.frame(
        V1 = c("Start", "2025-03-01T08:30:00"),
        V2 = c("Device", "Moxy"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_s3_class(result, "POSIXct")
    expect_equal(
        result,
        as.POSIXct("2025-03-01T08:30:00", format = "%Y-%m-%dT%H:%M:%OS")
    )
})

test_that("extract_start_timestamp() detects yyyy-mm-dd HH:MM:SS format", {
    data <- data.frame(
        V1 = c("Start", "2025-06-15 14:22:10"),
        V2 = c("Device", "Moxy"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_s3_class(result, "POSIXct")
    expect_equal(result, as.POSIXct("2025-06-15 14:22:10"))
})

test_that("extract_start_timestamp() returns earliest timestamp when multiple present", {
    data <- data.frame(
        V1 = c("2025-01-01 10:00:00", "2025-01-01 11:00:00"),
        V2 = c("2025-01-01 09:00:00", "metadata"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_equal(result, as.POSIXct("2025-01-01 09:00:00"))
})

test_that("extract_start_timestamp() compares mixed formats chronologically", {
    data <- data.frame(
        V1 = c("2025-01-01 00:00:00", "31/12/2024 23:00:00")
    )

    expect_equal(
        extract_start_timestamp(data),
        as.POSIXct("31/12/2024 23:00:00", format = "%d/%m/%Y %H:%M:%OS")
    )
})

test_that("extract_start_timestamp() ignores NA and empty strings", {
    data <- data.frame(
        V1 = c(NA, ""),
        V2 = c("2025-05-10T07:00:00", "invalid:time"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_s3_class(result, "POSIXct")
    expect_equal(
        result,
        as.POSIXct("2025-05-10T07:00:00", format = "%Y-%m-%dT%H:%M:%OS")
    )
})

test_that("extract_start_timestamp() works with real example file header", {
    file_header <- read_file(example_mnirs("moxy_ramp.xlsx"))[1:20, ]
    result <- extract_start_timestamp(file_header)
    expect_s3_class(result, "POSIXct")
})


## parse_time_channel() ================================================
test_that("parse_time_channel() returns a list of time, timestamp, start_timestamp", {
    result <- parse_time_channel(c(0, 1, 2))

    expect_type(result, "list")
    expect_named(result, c("time", "timestamp", "start_timestamp"))
    expect_type(result$time, "double")
})

test_that("parse_time_channel() coerces numeric-string time to numeric", {
    result <- parse_time_channel(c("0", "1", "2"))

    expect_type(result$time, "double")
    expect_equal(result$time, c(0, 1, 2))
})

test_that("parse_time_channel() preserves numeric time (zero_time = FALSE)", {
    result <- parse_time_channel(c(10.5, 20.5, 30.5), zero_time = FALSE)

    expect_equal(result$time, c(10.5, 20.5, 30.5))
    expect_null(result$timestamp)
    expect_null(result$start_timestamp)
})

test_that("parse_time_channel() recalculates numeric time from zero", {
    result <- parse_time_channel(c(10, 20, 30), zero_time = TRUE)

    expect_equal(result$time, c(0, 10, 20))
})

test_that("parse_time_channel() parses ISO 8601 character timestamps to numeric", {
    result <- parse_time_channel(c("2025-01-01T10:00:00", "2025-01-01T10:00:01"))

    expect_type(result$time, "double")
    expect_equal(result$time, c(0, 1))
})

test_that("parse_time_channel() parses date-time formats to numeric", {
    formats <- list(
        c("2025-01-01 10:00:00", "2025-01-01 10:00:01"),
        c("2025/01/01 10:00:00", "2025/01/01 10:00:01"),
        c("01-01-2025 10:00:00", "01-01-2025 10:00:01"),
        c("01/01/2025 10:00:00", "01/01/2025 10:00:01")
    )

    for (fmt in formats) {
        expect_equal(parse_time_channel(fmt)$time, c(0, 1))
    }
})

test_that("parse_time_channel() parses time-only H:MM:SS character format", {
    result <- parse_time_channel(c("10:00:00", "10:00:01"))

    expect_equal(result$time, c(0, 1))
})

test_that("parse_time_channel() handles milliseconds in timestamps", {
    result <- parse_time_channel(
        c("2025-01-01T10:00:00.123", "2025-01-01T10:00:01.456")
    )

    expect_equal(result$time[1], 0)
    expect_true(result$time[2] > 1 & result$time[2] < 2)
})

test_that("parse_time_channel() converts POSIXct to numeric seconds from zero", {
    t0 <- as.POSIXct("2025-01-01 10:00:00")
    result <- parse_time_channel(t0 + c(0, 1))

    expect_type(result$time, "double")
    expect_equal(result$time, c(0, 1))
    expect_s3_class(result$timestamp, "POSIXct")
    expect_equal(as.numeric(result$timestamp), as.numeric(t0 + c(0, 1)))

    ## POSIXct is always relative — always starts from 0
    result <- parse_time_channel(t0 + c(0, 1), zero_time = FALSE)
    expect_equal(result$time, c(0, 1))
})

test_that("parse_time_channel() returns start_timestamp from POSIXct time", {
    t0 <- as.POSIXct("2025-03-15 08:00:00")
    result <- parse_time_channel(t0 + c(0, 1, 2))

    ## start_timestamp is extracted from the POSIXct series when not in header
    expect_equal(result$start_timestamp, t0, ignore_attr = TRUE)
    expect_equal(result$timestamp[1L], result$start_timestamp, ignore_attr = TRUE)
})

test_that("parse_time_channel() start_timestamp is first sample when non-monotonic", {
    t0 <- as.POSIXct("2025-01-01 10:00:00")
    result <- parse_time_channel(t0 + c(5, 0, 10))

    ## start_timestamp + time must reconstruct the original timestamps
    expect_equal(result$start_timestamp, t0 + 5, ignore_attr = TRUE)
    expect_equal(as.numeric(result$timestamp), as.numeric(t0 + c(5, 0, 10)))
})

test_that("parse_time_channel() does not force header for dated time series", {
    t0 <- as.POSIXct("2025-03-15 08:00:00")

    result <- parse_time_channel(
        t0 + 0:2,
        start_timestamp = stop("header timestamp was forced")
    )

    expect_equal(result$start_timestamp, t0, ignore_attr = TRUE)
})

test_that("parse_time_channel() dated time series takes priority over header", {
    t0 <- as.POSIXct("2025-03-15 08:00:00")
    header_start <- as.POSIXct("2024-01-01 00:00:00")

    result <- parse_time_channel(
        format(t0 + 0:2, "%Y-%m-%d %H:%M:%S"),
        start_timestamp = header_start
    )

    expect_equal(result$start_timestamp, t0)
    expect_equal(result$timestamp, t0 + 0:2, ignore_attr = TRUE)
})

test_that("parse_time_channel() header anchors time-only series", {
    header_start <- as.POSIXct("2025-03-15 10:00:00")

    result <- parse_time_channel(
        c("10:00:00", "10:00:01"),
        start_timestamp = header_start
    )

    expect_equal(result$start_timestamp, header_start)
    expect_equal(result$timestamp, header_start + 0:1)
})

test_that("parse_time_channel() header anchors fractional-day series", {
    header_start <- as.POSIXct("2025-03-15 12:00:00")

    result <- parse_time_channel(
        c(0.5, 0.5 + 1 / 86400),
        start_timestamp = header_start
    )

    expect_equal(result$start_timestamp, header_start)
    expect_equal(result$timestamp, header_start + 0:1, tolerance = 1e-6)
})

test_that("parse_time_channel() header reconstructs absolute timestamps", {
    start_ts <- as.POSIXct("2025-06-01 09:00:00")

    result <- parse_time_channel(c(0, 1, 2), start_timestamp = start_ts)

    expect_s3_class(result$timestamp, "POSIXct")
    expect_equal(as.numeric(result$timestamp), as.numeric(start_ts + c(0, 1, 2)))
    expect_equal(result$start_timestamp, start_ts)
})

test_that("parse_time_channel() returns NULL timestamps when none available", {
    result <- parse_time_channel(c(10, 20, 30), start_timestamp = NULL)

    expect_null(result$timestamp)
    expect_null(result$start_timestamp)
})

test_that("read_mnirs() add_timestamp inserts POSIXct column after time_channel", {
    file_path <- tempfile(fileext = ".csv")
    on.exit(unlink(file_path))
    writeLines(
        c(
            "SmO2,recorded_at",
            "55,2025-03-15 08:00:00",
            "56,2025-03-15 08:00:01"
        ),
        file_path
    )

    df <- read_mnirs(
        file_path,
        nirs_channels = "SmO2",
        time_channel = c(time = "recorded_at"),
        add_timestamp = TRUE,
        verbose = FALSE
    )

    expect_equal(names(df), c("time", "timestamp", "SmO2"))
    expect_s3_class(df$timestamp, "POSIXct")
    expect_equal(df$timestamp[1L], attr(df, "start_timestamp"), ignore_attr = TRUE)

    ## no timestamp available — column not added silently
    writeLines(c("SmO2,time", "55,10", "56,11"), file_path)
    df <- read_mnirs(
        file_path,
        nirs_channels = "SmO2",
        time_channel = "time",
        add_timestamp = TRUE,
        verbose = FALSE
    )

    expect_false("timestamp" %in% names(df))
    expect_null(attr(df, "start_timestamp"))
})

test_that("parse_time_channel() works on fraction-of-day", {
    ## Moxy.csv saved as excel will coerce date-time
    ## to numeric fraction-of-day
    file_path <- test_path("testdata/moxy-occlusion.xlsx")
    skip_if_not(file.exists(file_path), "testdata not available")

    data <- suppressMessages(readxl::read_excel(file_path)[-(1:2), 1:2])
    x <- as.numeric(data[[2L]])

    result <- parse_time_channel(x, start_timestamp = NULL)

    expect_equal(class(result$time), "numeric")
    expect_equal(result$time[1L], 0)
    expect_equal(median(diff(result$time)), 2)
    expect_equal(class(result$timestamp), c("POSIXct", "POSIXt"))
    expect_equal(class(result$start_timestamp), c("POSIXct", "POSIXt"))
    ## should return today's date, local time zone, precise timestamp
    expect_equal(as.Date(result$start_timestamp), Sys.Date())
    expect_equal(format(result$start_timestamp, "%Z"), format(Sys.time(), "%Z"))
    expect_equal(format(result$start_timestamp, "%H:%M:%OS"), "13:52:59")
})

test_that("read_mnirs() returns local time zone start_timestamp", {
    perfpro <- test_path("testdata/perfpro-mre.xlsx")
    moxy_occl <- test_path("testdata/moxy-occlusion.xlsx")
    skip_if_not(file.exists(perfpro), "testdata not available")
    skip("run in local time PDT/PST zone")

    file_list <- c(
        perfpro, ## today's date, 0:00:00
        moxy_occl, ## today's date, 13:52:59
        example_mnirs("moxy_intervals"), ## today's date, 13:17:13
        example_mnirs("moxy_ramp") ## today's date 0:29:00.41
    )

    timestamp_list <- lapply(file_list, \(.file) {
        df <- read_mnirs(
            .file, add_timestamp = TRUE, zero_time = TRUE, verbose = FALSE
        )
        start_timestamp <- attr(df, "start_timestamp")
        expect_true(format(start_timestamp, "%Z") %in% c("PDT", "PST"))
        start_timestamp
    })
    timestamp_list <- do.call(c, timestamp_list)
    expect_equal(as.Date(timestamp_list), rep(Sys.Date(), 4))

    expect_equal(
        format(timestamp_list, "%H:%M:%OS2"),
        c("00:00:00.00", "13:52:59.00", "13:17:13.00", "00:29:00.41")
    )
})

## oxysoft_sample_rate() ==============================================
test_that("oxysoft_sample_rate() reads export sample rate from header", {
    raw <- read_file(example_mnirs("artinis_intervals"))
    header <- raw[seq_len(detect_mnirs_device(raw)$header_row), ]

    expect_equal(oxysoft_sample_rate(header), 10)
})

test_that("read_mnirs() derives Artinis time from sample index", {
    df <- read_mnirs(
        example_mnirs("artinis_intervals"),
        nirs_channels = c(HHb = 2, O2Hb = 3),
        time_channel = c(sample = 1),
        verbose = FALSE
    )

    expect_equal(names(df)[1:2], c("time", "sample"))
    expect_equal(attr(df, "time_channel"), "time")
    expect_equal(attr(df, "sample_rate"), 10)
    expect_equal(df$time, df$sample / 10)

    ## existing "time" name is kept unique
    df <- read_mnirs(
        example_mnirs("artinis_intervals"),
        nirs_channels = c(HHb = 2),
        time_channel = c(time = 1),
        verbose = FALSE
    )

    expect_equal(names(df)[1:2], c("time_1", "time"))
    expect_equal(attr(df, "time_channel"), "time_1")
})

## detect_irregular_samples() =========================================
test_that("detect_irregular_samples returns invisibly with no irregularities", {
    x <- seq(0, 100, by = 1)
    expect_invisible(detect_irregular_samples(x, "time", verbose = TRUE))
    expect_invisible(detect_irregular_samples(x, "time", verbose = FALSE))
})

test_that("detect_irregular_samples detects duplicated samples", {
    x <- c(0, 1, 2, 2, 3, 4)
    expect_warning(
        detect_irregular_samples(x, "time"),
        "irregular.*detected"
    )
    expect_warning(
        detect_irregular_samples(x, "time"),
        "time.*=.*2"
    )
})

test_that("detect_irregular_samples detects unordered samples", {
    x <- c(0, 1, 3, 2, 4, 5)
    expect_warning(
        detect_irregular_samples(x, "time"),
        "irregular.*detected"
    )
    expect_warning(
        detect_irregular_samples(x, "time"),
        "time.*=.*2"
    )
})

test_that("detect_irregular_samples detects large gaps (>= 3600)", {
    x <- c(0, 1, 2, 3602, 3603)
    expect_warning(
        detect_irregular_samples(x, "time"),
        "irregular.*detected"
    )
    expect_warning(
        detect_irregular_samples(x, "time"),
        "time.*=.*2"
    )
})

test_that("detect_irregular_samples shows first 3 when > 5 irregularities", {
    x <- c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
    expect_warning(
        detect_irregular_samples(x, "time"),
        "and 3 more"
    )
})

test_that("detect_irregular_samples shows all when <= 5 irregularities", {
    x <- c(0, 1, 1, 2, 2, 3)
    w <- expect_warning(detect_irregular_samples(x, "time"))
    expect_false(grepl("more", w$message))
})

test_that("detect_irregular_samples respects verbose = FALSE", {
    x <- c(0, 1, 1, 2, 2)
    expect_invisible(detect_irregular_samples(x, "time", verbose = FALSE))
    expect_no_warning(detect_irregular_samples(x, "time", verbose = FALSE))
})

test_that("detect_irregular_samples handles multiple irregularity types", {
    x <- c(0, 1, 1, 3, 2, 3606)
    expect_warning(
        detect_irregular_samples(x, "time"),
        "irregular.*detected"
    )
})

test_that("detect_irregular_samples uses correct time_channel name", {
    x <- c(0, 1, 1, 2)
    expect_warning(
        detect_irregular_samples(x, "my_time_col"),
        "my_time_col"
    )
})

## read_mnirs() =======================================================
test_that("read_mnirs preserves dated time series timestamps", {
    file_path <- tempfile(fileext = ".csv")
    on.exit(unlink(file_path))
    writeLines(
        c(
            "recorded_at,SmO2",
            "2025-03-15 08:00:00,55",
            "2025-03-15 08:00:01,56"
        ),
        file_path
    )
    expected <- as.POSIXct(c(
        "2025-03-15 08:00:00",
        "2025-03-15 08:00:01"
    ))

    result <- read_mnirs(
        file_path,
        nirs_channels = "SmO2",
        time_channel = "recorded_at",
        add_timestamp = TRUE,
        verbose = FALSE
    )

    expect_s3_class(attr(result, "start_timestamp"), "POSIXct")
    expect_equal(
        attr(result, "start_timestamp"), min(expected), ignore_attr = TRUE
    )
    expect_equal(result$timestamp, expected)
    expect_equal(result$recorded_at, c(0, 1))
})

## read_mnirs() auto-detection =========================================
test_that("read_mnirs auto-detects Moxy channels when nirs_channels = NULL", {
    file_path <- example_mnirs("moxy_ramp")

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = NULL,
            time_channel = NULL,
            verbose = TRUE
        ),
        "Moxy.*detected"
    ) |>
        expect_warning("irregular") |>
        expect_message("Estimated.*sample_rate.*2")

    expect_s3_class(df, "mnirs")
    expect_equal(attr(df, "nirs_device"), "Moxy")
    expect_equal(attr(df, "nirs_channels"), c("SmO2 Live", "SmO2 Live(2)"))
    expect_equal(attr(df, "time_channel"), device_patterns$Moxy$time_channel)
})

test_that("read_mnirs auto-detects Train.Red channels when nirs_channels = NULL", {
    file_path <- example_mnirs("train.red")

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = NULL,
            time_channel = NULL,
            verbose = TRUE
        ),
        "Train.Red.*detected"
    ) |>
        expect_warning("Duplicate channel names") |>
        expect_warning("irregular") |>
        expect_message("Estimated.*sample_rate.*10")

    expect_s3_class(df, "mnirs")
    expect_equal(attr(df, "nirs_device"), "Train.Red")
    expect_equal(attr(df, "nirs_channels"), c("SmO2", "SmO2_1"))
    expect_equal(
        attr(df, "time_channel"),
        device_patterns$Train.Red$time_channel
    )
})

test_that("read_mnirs accepts renamed duplicate channels as originals", {
    df <- read_mnirs(
        example_mnirs("train.red"),
        nirs_channels = c(smo2_left = "SmO2", smo2_right = "SmO2_1"),
        verbose = FALSE
    )
    expect_equal(attr(df, "nirs_channels"), c("smo2_left", "smo2_right"))
})

test_that("read_mnirs auto-detects Artinis channels when nirs_channels = NULL", {
    file_path <- example_mnirs("artinis_intervals")

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = NULL,
            time_channel = NULL,
            verbose = TRUE
        ),
        "Artinis.*detected"
    ) |>
        expect_message("Oxysoft.*sample_rate.*10")

    expect_s3_class(df, "mnirs")
    expect_equal(attr(df, "nirs_device"), "Artinis")
    expect_equal(attr(df, "nirs_channels"), c("vl_o2hb", "vl_hhb"))
    expect_equal(attr(df, "time_channel"), "time")
    expect_equal(attr(df, "event_channel"), "event")
    expect_true(all(c("sample", "time", "vl_o2hb", "vl_hhb", "event") %in% names(df)))
})

test_that("read_mnirs keep_all = FALSE returns only specified columns by default", {
    file_path <- example_mnirs("moxy_ramp")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(smo2 = "SmO2 Live"),
        time_channel = c(time = "hh:mm:ss"),
        # keep_all = FALSE,
        verbose = FALSE
    )

    ## default keep_all = FALSE returns only smo2 + time
    expect_equal(ncol(df), 2)
    expect_true("smo2" %in% names(df))
    expect_true("time" %in% names(df))
})

test_that("read_mnirs keep_all = TRUE returns all columns", {
    file_path <- example_mnirs("moxy_ramp")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(smo2 = "SmO2 Live"),
        time_channel = c(time = "hh:mm:ss"),
        keep_all = TRUE,
        verbose = FALSE
    )

    ## default keep_all = TRUE returns more columns than just smo2 + time
    expect_gt(ncol(df), 2)
    expect_true("smo2" %in% names(df))
    expect_true("time" %in% names(df))
})

test_that("read_mnirs returns all columns when auto-detecting nirs_channels", {
    file_path <- example_mnirs("moxy_ramp")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = NULL,
        time_channel = NULL,
        keep_all = FALSE, ## explicitly FALSE should be overridden
        verbose = FALSE
    )

    ## default keep_all = TRUE returns more columns than just smo2 + time
    expect_gt(ncol(df), 2)
    expect_true("SmO2 Live" %in% names(df))
    expect_true("hh:mm:ss" %in% names(df))
})


## moxy ===============================================================
test_that("read_mnirs moxy .xlsx works with timestamp", {
    expect_warning(
        df <- read_mnirs(
            file_path = example_mnirs("moxy_ramp.xlsx"),
            nirs_channels = c(
                smo2_left = "SmO2 Live",
                smo2_right = "SmO2 Live(2)"
            ),
            time_channel = c(time = "hh:mm:ss"),
            add_timestamp = TRUE,
            verbose = TRUE
        ),
        "irregular"
    ) |>
        expect_message("Estimated.*sample_rate.*2")

    expect_s3_class(df, "mnirs")
    expect_s3_class(df, "data.frame")
    expect_true(all(
        c("time", "timestamp", "smo2_left", "smo2_right") %in% names(df)
    ))
    expect_equal(df$time[1], 0)
    expect_equal(class(df$time), "numeric")
    expect_true(any(class(df$timestamp) %in% "POSIXct"))

    ## check that time diffs should be 0 < Δ < 1 with proper POSIXct import
    expect_gt(sum(diff(df$time[1:100]) < 1 & diff(df$time[1:100]) > 0), 0)
    expect_lt(sum(diff(df$time[1:100]) %in% c(0, 1)), 99)

    expect_true(all(
        c("nirs_device", "nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))

    expect_equal(attr(df, "nirs_device"), "Moxy")
    expect_equal(attr(df, "sample_rate"), 2)
})

test_that("read_mnirs moxy .csv works converting time to numeric", {
    expect_message(
        df <- read_mnirs(
            file_path = example_mnirs("moxy_intervals"),
            nirs_channels = c(smo2_left = "SmO2 Live", thb = "THb"),
            time_channel = c(time = "hh:mm:ss"),
            add_timestamp = FALSE,
            verbose = TRUE
        ),
        "Estimated.*sample_rate.*0.5"
    )

    expect_equal(class(df$time), "numeric")
    expect_false(class(df$time) %in% "POSIXct")
    expect_false(c("timestamp") %in% names(df))

    ## check that time diffs should be Δ >= 2 with proper POSIXct import
    expect_gt(sum(diff(df$time[1:100]) >= 2), 0)

    expect_true(all(
        c("nirs_device", "nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))

    expect_equal(attr(df, "nirs_device"), "Moxy")
    expect_equal(attr(df, "sample_rate"), 0.5)
})


test_that("read_mnirs moxy invalid channel names", {
    file_path <- example_mnirs("moxy_ramp.xlsx")

    ## invalid channel names
    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(""),
            time_channel = c(time = "hh:mm:ss"),
            verbose = FALSE
        ),
        "not detected"
    )

    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(smo2_left = "smo2_doesnt_exist"),
            time_channel = c(time = "hh:mm:ss"),
            verbose = FALSE
        ),
        "not detected"
    )

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 Live",
                smo2_right = "SmO2 Live(2)"
            ),
            time_channel = NULL,
            verbose = TRUE
        ),
        "Estimated.*sample_rate.*2"
    ) |>
        expect_warning("irregular")

    ## device default time_channel used
    expect_equal(attr(df, "time_channel"), "hh:mm:ss")

    ## duplicate input names are renamed
    expect_warning(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(smo2 = "SmO2 Live", smo2 = "SmO2 Live(2)"),
            time_channel = c(time = "hh:mm:ss"),
            verbose = TRUE
        ),
        "Duplicate"
    ) |>
        expect_warning("irregular") |>
        expect_message("Estimated.*sample_rate.*2")

    expect_true(all(c("smo2", "smo2_1") %in% names(df)))
})


## train.red ========================================================
test_that("read_mnirs train.red works", {
    file_path <- example_mnirs("train.red_intervals.csv")

    expect_no_message(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            verbose = FALSE
        )
    )

    expect_warning(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            verbose = TRUE
        ),
        "irregular"
    ) |>
        expect_message("Estimated.*sample_rate.*10")

    expect_s3_class(df, "mnirs")
    expect_s3_class(df, "data.frame")
    expect_true(all(
        c("time", "smo2_left", "smo2_right") %in% names(df)
    ))
    expect_equal(class(df$time), "numeric")
    expect_gte(df$time[1], 0)

    ## check that time diffs should be 0 < Δ < 1 with proper POSIXct import
    expect_gt(sum(diff(df$time[1:100]) < 1 & diff(df$time[1:100]) > 0), 0)
    expect_lt(sum(diff(df$time[1:100]) %in% c(0, 1)), 99)

    expect_true(all(
        c("nirs_device", "nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))

    expect_equal(attr(df, "nirs_device"), "Train.Red")
    expect_equal(attr(df, "sample_rate"), 10)
})

test_that("read_mnirs external train.red mre works", {
    file_path <- test_path("testdata/train.red-mre.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    expect_no_message(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            verbose = FALSE
        )
    )

    expect_warning(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            verbose = TRUE
        ),
        "irregular"
    ) |>
        expect_message("Estimated.*sample_rate.*10")

    expect_s3_class(df, "mnirs")
    expect_s3_class(df, "data.frame")
    expect_true(all(
        c("time", "smo2_left", "smo2_right") %in% names(df)
    ))
    expect_equal(class(df$time), "numeric")
    expect_gte(df$time[1], 0)

    ## check that time diffs should be 0 < Δ < 1 with proper POSIXct import
    expect_gt(sum(diff(df$time[1:100]) < 1 & diff(df$time[1:100]) > 0), 0)
    expect_lt(sum(diff(df$time[1:100]) %in% c(0, 1)), 99)

    expect_true(all(
        c("nirs_device", "nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))

    expect_equal(attr(df, "nirs_device"), "Train.Red")
    expect_equal(attr(df, "sample_rate"), 10)
    expect_s3_class(attr(df, "start_timestamp"), "POSIXct")
})

test_that("read_mnirs coerces integerish event_channel to integer", {
    ## train.red "Lap/Event" column contains lap numbers (doubles from CSV)
    data <- read_mnirs(
        file_path = example_mnirs("train.red_intervals.csv"),
        nirs_channels = c(smo2_left = "SmO2 unfiltered"),
        time_channel = c(time = "Timestamp (seconds passed)"),
        event_channel = c(lap = "Lap/Event"),
        verbose = FALSE
    )
    # rlang::is_integerish(df$lap)

    expect_type(data$lap, "integer")
})

test_that("read_mnirs auto-detects Train.Red event_channel", {
    df <- read_mnirs(example_mnirs("train.red_intervals.csv"), verbose = FALSE)

    expect_equal(attr(df, "event_channel"), "Lap/Event")
    expect_type(df[["Lap/Event"]], "integer")
})

test_that("read_mnirs train.red works with zero_time", {
    file_path <- example_mnirs("train.red_intervals.csv")

    expect_equal(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            zero_time = TRUE,
            verbose = FALSE
        )$time[1],
        0
    )
})

test_that("read_mnirs train.red invalid channel names", {
    file_path <- example_mnirs("train.red_intervals.csv")

    ## invalid channel names
    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(" "),
            time_channel = c(time = "Timestamp (seconds passed)"),
        ),
        "not detected"
    )

    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(smo2_left = "smo2_doesnt_exist"),
            time_channel = c(time = "Timestamp (seconds passed)"),
        ),
        "not detected"
    )

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_left = "SmO2 unfiltered",
                smo2_right = "SmO2 unfiltered"
            ),
            time_channel = NULL,
            verbose = TRUE
        ),
        "Estimated.*sample_rate.*10"
    ) |>
        expect_warning("irregular")

    expect_equal(attr(df, "time_channel"), "Timestamp (seconds passed)")

    ## duplicate input names are renamed
    expect_warning(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2 = "SmO2 unfiltered",
                smo2 = "SmO2 unfiltered"
            ),
            time_channel = c(time = "Timestamp (seconds passed)"),
            verbose = TRUE
        ),
        "Duplicate"
    ) |>
        expect_warning("irregular") |>
        expect_message("Estimated.*sample_rate.*10")

    expect_true(all(c("smo2", "smo2_1") %in% names(df)))
})


## oxysoft =======================================================
test_that("read_mnirs oxysoft works", {
    file_path <- example_mnirs("artinis_intervals")

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(HHb = 2, O2Hb = 3),
            time_channel = c(sample = 1),
            verbose = TRUE
        ),
        "Oxysoft.*sample_rate.*10"
    )

    expect_s3_class(df, "mnirs")
    expect_s3_class(df, "data.frame")
    ## "event" auto-detected from the legend
    expect_true(all(
        c("time", "HHb", "O2Hb", "event") %in% names(df)
    ))
    expect_equal(attr(df, "event_channel"), "event")
    expect_equal(class(df$time), "numeric")
    expect_gte(df$time[1], 0)
    expect_equal(df$sample[1:10] / 10, df$time[1:10])

    expect_true(all.equal(diff(df$time[1:100]), rep(0.1, 99)))

    expect_true(all(
        c("nirs_device", "nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))

    expect_equal(attr(df, "nirs_device"), "Artinis")
    expect_equal(attr(df, "sample_rate"), 10)
    expect_equal(attr(df, "time_channel"), "time")
})

test_that("read_mnirs Oxysoft Portamon works", {
    file_path <- example_mnirs("portamon-oxcap")

    expect_equal(
        read_file(file_path) |>
            detect_mnirs_device(),
        list(
            nirs_device = "Artinis",
            header_row = 42
        )
    )

    df <- read_mnirs(
        file_path = example_mnirs("portamon-oxcap.xlsx"),
        nirs_channels = c(thb = 2, hhb = 3, o2hb = 4),
        time_channel = NULL,
        event_channel = c(event = "col_6"),
        verbose = FALSE
    )

    expect_true(all(
        c("sample", "time", "event", "thb", "hhb", "o2hb") %in% names(df)
    ))

    ## auto detect channels from the legend
    df2 <- read_mnirs(
        file_path,
        nirs_channels = NULL,
        time_channel = NULL,
        verbose = FALSE
    )

    expect_true(all(
        c("sample", "time", "rx1_tx1_thb", "rx1_tx1_hhb", "rx1_tx1_o2hb",
          "event", "labels") %in% names(df2)
    ))
    expect_equal(
        attr(df2, "nirs_channels"),
        c("rx1_tx1_thb", "rx1_tx1_hhb", "rx1_tx1_o2hb")
    )
    expect_equal(attr(df2, "event_channel"), "event")
    expect_true("Occlusion" %in% df2$labels)

    for (d in list(df, df2)) {
        expect_s3_class(d, "mnirs")
        expect_s3_class(d, "data.frame")

        expect_equal(class(d$time), "numeric")
        expect_equal(d$time[1], 0)
        expect_equal(d$sample[1:10] / 10, d$time[1:10])

        expect_true(all.equal(diff(d$time[1:100]), rep(0.1, 99)))

        expect_true(all(
            c(
                "nirs_device",
                "nirs_channels",
                "time_channel",
                "sample_rate"
            ) %in% names(attributes(d))
        ))

        expect_equal(attr(d, "nirs_device"), "Artinis")
        expect_equal(attr(d, "sample_rate"), 10)
        expect_equal(attr(d, "time_channel"), "time")
    }
})

test_that("read_mnirs Oxysoft event_channel = 'labels' aliases label column", {
    df <- read_mnirs(
        example_mnirs("portamon-oxcap"),
        nirs_channels = c(thb = 2, hhb = 3, o2hb = 4),
        event_channel = c(event = "labels"),
        verbose = FALSE
    )

    expect_equal(attr(df, "event_channel"), "event")
    expect_true("Occlusion" %in% df$event)
    expect_false("labels" %in% names(df))
})

test_that("read_mnirs Oxysoft edge case channel names", {
    file_path <- example_mnirs("artinis_intervals")

    old_verbose <- getOption("mnirs.verbose")
    on.exit(options(mnirs.verbose = old_verbose), add = TRUE)
    options(mnirs.verbose = FALSE)

    ## invalid channel names
    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(""),
            time_channel = c(sample = 1),
        ),
        "not detected"
    )

    expect_error(
        read_mnirs(
            file_path = file_path,
            nirs_channels = c(smo2_left = "smo2_doesnt_exist"),
            time_channel = c(sample = 1),
        ),
        "not detected"
    )

    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(HHb = 2, O2Hb = 3),
            time_channel = NULL,
            verbose = TRUE
        ),
        "Oxysoft.*sample_rate.*10"
    )

    ## detected as "sample" then updated to "time" automatically
    expect_equal(attr(df, "time_channel"), "time")

    ## duplicate input names are renamed
    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(HHb = 2, HHb = 3),
            time_channel = c(sample = 1),
            verbose = TRUE
        ),
        "Oxysoft.*sample_rate.*10"
    ) |>
        expect_warning("Duplicate")

    expect_true(all(c("HHb", "HHb_1") %in% names(df)))
})

## VO2master app ========================================================
test_that("read_mnirs VO2master with ',' decimals returns numeric", {
    file_path <- test_path("testdata/vo2master.csv")
    skip_if_not(file.exists(file_path), "testdata not available")

    nirs_channels <- c(
        smo2_1 = "SmO2[%]",
        smo2_2 = "SmO2 -  2[%]",
        smo2_3 = "SmO2 -  3[%]"
    )
    time_channel <- c(time = "Time[s]")

    df_raw <- tibble::as_tibble(
        data.table::fread(
            file_path,
            header = TRUE,
            colClasses = "character"
        )
    )

    expect_all_true(vapply(df_raw, is.character, logical(1L)))

    ## should convert decimal "," to numeric
    df <- convert_type(df_raw, list(time = time_channel))
    expect_all_true(vapply(df[, -c(1:2)], is.numeric, logical(1L)))

    ## integrated test
    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(
                smo2_1 = "SmO2[%]",
                smo2_2 = "SmO2 -  2[%]",
                smo2_3 = "SmO2 -  3[%]"
            ),
            time_channel = c(time = "Time[s]"),
            verbose = TRUE
        ),
        "Estimated.*sample_rate.*1"
    )

    expect_type(df$time, "double")
    expect_equal(sum(diff(df$time[1:100]) == 1), 99)
    ## smo2 should be numeric from "27,90"
    expect_type(df$smo2_1, "double")
    expect_type(df$smo2_2, "double")
    expect_type(df$smo2_3, "double")

    expect_true(all(
        c("nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))
    expect_equal(attr(df, "sample_rate"), 1)
})

## PerfPro ========================================================
test_that("read_mnirs PerfPro", {
    file_path <- test_path("testdata/perfpro-mre.xlsx")
    skip_if_not(file.exists(file_path), "testdata not available")

    raw <- read_file(file_path)
    device <- detect_mnirs_device(raw)

    expect_equal(device$nirs_device, "PerfPro")
    expect_equal(device$header_row, 3)

    expect_message(
        channels <- resolve_channels(raw, device, test_user(), verbose = TRUE),
        "PerfPro"
    )

    expect_equal(channels$time, "Time")
    expect_equal(channels$nirs, c("SmO2 (1614)", "SmO2 (1615)"))

    ## integrated test
    expect_message(
        df <- read_mnirs(
            file_path = file_path,
            nirs_channels = c(smo2 = "SmO2 (1614)"),
            time_channel = c(time = "Time"),
            add_timestamp = TRUE
        ),
        "Estimated.*sample_rate.*2"
    )

    expect_type(df$time, "double")
    expect_true(inherits(df$timestamp, "POSIXct"))
    expect_equal(df$time[1L], 0)
    expect_true(all(diff(df$time[1:25]) < 1))
    ## smo2 should be numeric from "27,90"
    expect_type(df$smo2, "double")

    expect_true(all(
        c("nirs_channels", "time_channel", "sample_rate") %in%
            names(attributes(df))
    ))
    expect_equal(attr(df, "sample_rate"), 2)
})

## create_mnirs_data() ==================================================

test_that("create_mnirs_data edge cases", {
    ## error when data isn't a dataframe
    vec <- c(1, 2)
    expect_error(create_mnirs_data(vec), "must be a data frame")

    ## attach metadata unlisted
    df <- tibble(A = 1:2, B = letters[1:2])
    df_meta <- create_mnirs_data(
        df,
        nirs_channels = c("A", "B"),
        sample_rate = 1
    )
    expect_null(attr(df, "sample_rate"))
    expect_equal(attr(df_meta, "sample_rate"), 1)
    expect_equal(attr(df_meta, "nirs_channels"), c("A", "B"))
})

test_that("create_mnirs_data accepts NSE for *_channels", {
    df <- tibble(A = 1:2, B = 3:4, C = 5:6, lap = c("a", "b"))

    ## bare symbols
    df_sym <- create_mnirs_data(
        df,
        nirs_channels = B,
        time_channel = A,
        event_channel = lap
    )
    expect_equal(attr(df_sym, "nirs_channels"), "B")
    expect_equal(attr(df_sym, "time_channel"), "A")
    expect_equal(attr(df_sym, "event_channel"), "lap")

    ## c() expression
    df_c <- create_mnirs_data(df, nirs_channels = c(B, C))
    expect_equal(attr(df_c, "nirs_channels"), c("B", "C"))

    ## tidyselect helper
    df_sel <- create_mnirs_data(
        df, nirs_channels = tidyselect::starts_with("B")
    )
    expect_equal(attr(df_sel, "nirs_channels"), "B")

    ## external character vector
    chans <- c("B", "C")
    df_ext <- create_mnirs_data(df, nirs_channels = chans)
    expect_equal(attr(df_ext, "nirs_channels"), c("B", "C"))

    ## list form still works (single list argument)
    meta <- list(nirs_channels = c("B", "C"), sample_rate = 2)
    df_list <- create_mnirs_data(df, meta)
    expect_equal(attr(df_list, "nirs_channels"), c("B", "C"))
    expect_equal(attr(df_list, "sample_rate"), 2)
})

test_that("create_mnirs_data preserves grouping", {
    skip_if_not_installed("dplyr")
    
    df <- tibble(g = c("a", "a", "b"), A = 1:3, B = 4:6)
    grouped <- dplyr::group_by(df, g)

    df_grp <- create_mnirs_data(grouped, nirs_channels = "B", sample_rate = 1)

    ## grouping survives and mnirs class is intact
    expect_s3_class(df_grp, "mnirs")
    expect_true(dplyr::is_grouped_df(df_grp))
    expect_equal(dplyr::group_vars(df_grp), "g")
    expect_equal(attr(df_grp, "nirs_channels"), "B")

    ## ungrouped input stays ungrouped
    df_ungrp <- create_mnirs_data(df, nirs_channels = "B")
    expect_false(dplyr::is_grouped_df(df_ungrp))
})
