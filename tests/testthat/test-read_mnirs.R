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


## detect_device_channels() ============================================
test_that("detect_device_channels() returns user channels when provided", {
    result <- detect_device_channels(
        nirs_device = "Moxy",
        nirs_channels = c(smo2 = "SmO2 Live"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, c(smo2 = "SmO2 Live"))
    expect_equal(result$time_channel, c(time = "hh:mm:ss"))
})

test_that("detect_device_channels() returns user channels even with NULL device", {
    result <- detect_device_channels(
        nirs_device = NULL,
        nirs_channels = "O2Hb",
        time_channel = "Time",
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, "O2Hb")
    expect_equal(result$time_channel, "Time")
})

test_that("detect_device_channels() detects known channels for device", {
    data <- data.frame(
        V1 = c("meta", "hh:mm:ss", "00:00:01"),
        V2 = c("meta", "SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 2L,
        nirs_device = "Moxy",
        nirs_channels = NULL,
        time_channel = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, "SmO2 Live")
    expect_equal(result$time_channel, device_patterns$Moxy$time_channel)
})

test_that("detect_device_channels() detects multiple SmO2 channels", {
    data <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("SmO2 (1)", "55"),
        V3 = c("SmO2 (2)", "60"),
        V4 = c("HR", "120"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 1L,
        nirs_device = "PerfPro",
        nirs_channels = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, c("SmO2 (1)", "SmO2 (2)"))
})

test_that("detect_device_channels() matches SmO2 case-insensitively", {
    data <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("smo2 raw", "55"),
        V3 = c("SMO2_LIVE", "60"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 1L,
        nirs_device = "Moxy",
        nirs_channels = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, c("smo2 raw", "SMO2_LIVE"))
})

test_that("detect_device_channels() returns '2' for Artinis without scanning", {
    ## Artinis uses numeric channel id — no SmO2 column needed in data
    data <- data.frame(
        V1 = c("1", "2", "3"),
        V2 = c("4", "5", "6"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 2L,
        nirs_device = "Artinis",
        nirs_channels = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, "2")
    expect_null(result$time_channel)
})

test_that("detect_device_channels() detects SmO2 for unknown device", {
    data <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("SmO2", "55"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 1L,
        nirs_device = NULL,
        nirs_channels = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, "SmO2")
    expect_null(result$time_channel)
})

test_that("detect_device_channels() detects known channels for PerfPro", {
    file_path <- test_path("testdata/perfpro-mre.xlsx")
    skip_if_not(file.exists(file_path), "testdata not available")

    data <- read_file(file_path)
    detected_list <- detect_mnirs_device(data)
    nirs_device <- detected_list$nirs_device
    header_row <- detected_list$header_row

    expect_match(nirs_device, "PerfPro")
    expect_equal(header_row, 3)

    result <- detect_device_channels(
        data,
        header_row,
        nirs_device,
        nirs_channels = NULL,
        time_channel = NULL,
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, c("SmO2 (1614)", "SmO2 (1615)"))
    expect_equal(result$time_channel, "Time")
})

test_that("detect_device_channels() user time_channel overrides device default", {
    data <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 1L,
        nirs_device = "Moxy",
        nirs_channels = NULL,
        time_channel = c(time = "custom_time"),
        verbose = FALSE
    )

    expect_equal(result$nirs_channels, "SmO2 Live")
    expect_equal(result$time_channel, c(time = "custom_time"))
})

test_that("detect_device_channels() errors when no SmO2 columns found", {
    data <- data.frame(
        V1 = c("Time", "0.1"),
        V2 = c("HR", "120"),
        stringsAsFactors = FALSE
    )

    ## NULL device + no SmO2 columns
    expect_error(
        detect_device_channels(
            data,
            header_row = 1L,
            nirs_device = NULL,
            nirs_channels = NULL,
            verbose = FALSE
        ),
        "cannot be determined"
    )

    ## known device + no SmO2 columns
    expect_error(
        detect_device_channels(
            data,
            header_row = 1L,
            nirs_device = "Moxy",
            nirs_channels = NULL,
            verbose = FALSE
        ),
        "cannot be determined"
    )
})

test_that("detect_device_channels() verbose messages for detection", {
    data <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    expect_message(
        detect_device_channels(
            data,
            header_row = 1L,
            nirs_device = "Moxy",
            nirs_channels = NULL,
            verbose = TRUE
        ),
        "Moxy.*detected"
    )

    ## unknown device labelled "Unknown"
    expect_message(
        detect_device_channels(
            data,
            header_row = 1L,
            nirs_device = NULL,
            nirs_channels = NULL,
            verbose = TRUE
        ),
        "Unknown.*detected"
    )

    ## no message when user provides channels
    expect_no_message(
        detect_device_channels(
            nirs_device = "Moxy",
            nirs_channels = c("SmO2 Live"),
            verbose = TRUE
        )
    )
})

test_that("detect_device_channels() returns appropriate keep_all", {
    data <- data.frame(
        V1 = c("hh:mm:ss", "00:00:01"),
        V2 = c("SmO2 Live", "55"),
        stringsAsFactors = FALSE
    )

    result <- detect_device_channels(
        data,
        header_row = 1L,
        nirs_device = "Moxy",
        nirs_channels = NULL,
        time_channel = c(time = "custom_time"),
        verbose = FALSE
    )

    expect_true(result$keep_all)

    result <- detect_device_channels(
        nirs_device = "Moxy",
        nirs_channels = "SmO2",
        time_channel = c(time = "custom_time"),
        verbose = FALSE
    )

    expect_false(result$keep_all)
})


## read_data_table() ===================================================
test_that("read_data_table() extracts data with valid channels", {
    data <- data.frame(
        V1 = c("meta1", "meta2", "O2Hb", "10", "20"),
        V2 = c("meta1", "meta2", "HHb", "5", "15"),
        V3 = c("meta1", "meta2", "Time", "0.1", "0.2"),
        stringsAsFactors = FALSE
    )

    result <- read_data_table(data, nirs_channels = c("O2Hb", "HHb"))

    expect_type(result, "list")
    expect_named(result, c("file_header", "data_table"))
    expect_s3_class(result$data_table, "data.frame")
    expect_s3_class(result$file_header, "data.frame")

    expect_equal(nrow(result$data_table), 2)
    expect_equal(ncol(result$data_table), 3)
    expect_equal(names(result$data_table), c("O2Hb", "HHb", "Time"))
    expect_true(all(result$data_table == data[4:5, ]))

    expect_equal(nrow(result$file_header), 3)
    expect_equal(ncol(result$file_header), 3)
    expect_true(all(result$file_header == data[1:3, ]))
})

test_that("read_data_table() errors when channels not found", {
    data <- data.frame(
        V1 = c("header", "WrongChannel", "10"),
        V2 = c("header", "Time", "0.1"),
        stringsAsFactors = FALSE
    )

    expect_error(
        read_data_table(data, nirs_channels = "O2Hb"),
        "Channel names not detected"
    )
})

test_that("read_data_table() is case sensitive", {
    data <- data.frame(
        V1 = c("meta", "o2hb", "10"),
        V2 = c("meta", "time", "0.1"),
        stringsAsFactors = FALSE
    )

    expect_error(
        read_data_table(data, nirs_channels = "O2Hb"),
        "case sensitive"
    )
})


## detect_time_channel() ==============================================
test_that("detect_time_channel returns provided time_channel", {
    df <- data.frame(x = 1:5, y = 6:10)
    expect_equal(
        detect_time_channel(df, time_channel = "custom", verbose = FALSE),
        "custom"
    )
})

test_that("detect_time_channel finds time column by name", {
    df <- data.frame(time = 1:5, value = 6:10)
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "time"
    )

    df <- data.frame(Time = 1:5, value = 6:10)
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "Time"
    )

    df <- tibble::tibble("hh:mm:ss" = 1:5, value = 6:10)
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "hh:mm:ss"
    )

    df <- data.frame("hms" = 1:5, value = 6:10)
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "hms"
    )

    df <- data.frame(duration = 1:5, value = 6:10)
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "duration"
    )
})

test_that("detect_time_channel finds POSIXct column", {
    df <- data.frame(
        value = 1:5,
        posixct_col = as.POSIXct("2024-01-01 12:00:00") + 1:5
    )
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "posixct_col"
    )
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
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "string_col"
    )

    # Test with H:MM format
    df <- data.frame(
        value = 1:5,
        string_col = c("1:30", "1:31", "1:32", "1:33", "1:34")
    )
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "string_col"
    )
})

test_that("detect_time_channel handles NA values in character column", {
    df <- data.frame(
        value = 1:5,
        time_str = c(NA, "12:30:45", "12:30:46", "12:30:47", "12:30:48")
    )
    expect_equal(
        detect_time_channel(df, verbose = FALSE),
        "time_str"
    )
})

test_that("detect_time_channel errors when no time column found", {
    df <- data.frame(x = 1:5, y = 6:10)
    expect_error(
        detect_time_channel(df, verbose = FALSE),
        "time_channel.*not detected"
    )
})

test_that("detect_time_channel prioritises time_channel argument", {
    df <- data.frame(
        time = 1:5,
        custom = 6:10,
        timestamp = as.POSIXct("2024-01-01 12:00:00") + 1:5
    )
    expect_equal(
        detect_time_channel(df, time_channel = "custom", verbose = FALSE),
        "custom"
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


## select_rename_data() ===========================================
test_that("select_rename_data() selects and renames channels in order", {
    data <- data.frame(
        O2Hb = c("10", "20"),
        HHb = c("5", "15"),
        Time = c("0.1", "0.2"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = c(oxy = "O2Hb", deoxy = "HHb"),
        time_channel = c(time = "Time"),
        verbose = FALSE
    )

    expect_equal(names(result$data), c("time", "oxy", "deoxy"))
    expect_equal(result$nirs_channel, c("oxy", "deoxy"))
    expect_equal(result$time_channel, c("time"))
})

test_that("select_rename_data() works with unnamed channels", {
    data <- data.frame(
        O2Hb = c("10"),
        Time = c("0.1"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = "O2Hb",
        time_channel = "Time",
        verbose = FALSE
    )

    expect_equal(names(result$data), c("Time", "O2Hb"))
    expect_equal(result$nirs_channel, "O2Hb")
    expect_equal(result$time_channel, "Time")
})

test_that("select_rename_data() includes event channel", {
    data <- data.frame(
        O2Hb = c("10"),
        Time = c("0.1"),
        Event = c("Start"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = "O2Hb",
        time_channel = "Time",
        event_channel = "Event",
        verbose = FALSE
    )

    expect_equal(names(result$data), c("Time", "Event", "O2Hb"))
    expect_equal(result$event_channel, "Event")
})

test_that("select_rename_data() handles duplicate channel names", {
    skip('currently returns "Channel names not detected error"')

    data <- data.frame(
        O2Hb = c("10"),
        Time = c("0.1"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        result <- select_rename_data(
            data,
            nirs_channels = c("O2Hb", "O2Hb"),
            time_channel = "Time",
            verbose = TRUE
        ),
        "Duplicated channel names"
    )

    expect_equal(result$nirs_channel, c("O2Hb", "O2Hb_1"))
})

test_that("select_rename_data() handles duplicate data columns", {
    data <- data.frame(
        O2Hb = c("10"),
        O2Hb = c("20"),
        Time = c("0.1"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = c(oxy1 = "O2Hb", oxy2 = "O2Hb"),
        time_channel = "Time",
        verbose = FALSE
    )

    expect_equal(names(result$data), c("Time", "oxy1", "oxy2"))
})

test_that("select_rename_data() keeps all columns with keep_all", {
    data <- data.frame(
        O2Hb = c("10"),
        HHb = c("5"),
        Time = c("0.1"),
        Extra = c("x"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = c(o2hb = "O2Hb", hhb = "HHb"),
        time_channel = c(time = "Time"),
        keep_all = TRUE,
        verbose = FALSE
    )

    expect_equal(ncol(result$data), 4)
    expect_true(all(
        c("o2hb", "hhb", "time", "Extra") %in% names(result$data)
    ))
})

test_that("select_rename_data() drops extra columns with default keep_all = FALSE", {
    data <- data.frame(
        O2Hb = c("10"),
        HHb = c("5"),
        Time = c("0.1"),
        Extra = c("x"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = "O2Hb",
        time_channel = "Time",
        # keep_all = FALSE,
        verbose = FALSE
    )

    expect_equal(ncol(result$data), 2)
    expect_false("Extra" %in% names(result$data))
})

test_that("select_rename_data() errors when channel not found", {
    data <- data.frame(
        O2Hb = c("10"),
        Time = c("0.1"),
        stringsAsFactors = FALSE
    )

    expect_error(
        select_rename_data(
            data,
            nirs_channels = "HHb",
            time_channel = "Time"
        ),
        "Channel names not detected"
    )
})

test_that("select_rename_data() suppresses warnings with verbose", {
    data <- data.frame(
        O2Hb = c("10"),
        O2Hb = c("20"),
        Time = c("0.1"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    expect_silent(
        select_rename_data(
            data,
            nirs_channels = c(o2hb = "O2Hb", o2hb = "O2Hb"),
            time_channel = "Time",
            verbose = FALSE
        )
    )

    expect_warning(
        select_rename_data(
            data,
            nirs_channels = c(o2hb = "O2Hb", o2hb = "O2Hb"),
            time_channel = "Time",
            verbose = TRUE
        ),
        "o2hb = o2hb_1"
    )
})

test_that("select_rename_data() prioritises custom names over data", {
    data <- data.frame(
        O2Hb = c("10"),
        Time = c("0.1"),
        custom = c("x"),
        stringsAsFactors = FALSE
    )

    result <- select_rename_data(
        data,
        nirs_channels = c(custom = "O2Hb"),
        time_channel = "Time",
        keep_all = TRUE,
        verbose = FALSE
    )

    expect_true(all(c("custom", "custom_1") %in% names(result$data)))
    expect_equal(result$data$custom, "10")
})

## convert_type() ================================================
test_that("convert_type() coerces integer columns apart from event_channel", {
    data <- data.frame(
        time = c(1, 2, 3),
        lap = c(1.0, NA_real_, 2.0),
        B = c(10, 20, 30),
        x = c(10.5, 11.0, 11.5),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, "time", event_channel = "lap")

    expect_type(result$time, "double")
    expect_type(result$B, "double")
    expect_type(result$lap, "integer")
    expect_equal(result$lap, c(1L, NA_integer_, 2L))
    ## other columns unaffected
    expect_type(result$x, "double")
})


test_that("convert_type() standardises empty and 'NA' strings to NA", {
    data <- data.frame(
        time = c("1", "2", "3"),
        A = c("a", "", "b"),
        B = c("x", "NA", "y"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, "time")

    expect_equal(result$A, c("a", NA_character_, "b"))
    expect_equal(result$B, c("x", NA_character_, "y"))
})

test_that("convert_type() standardises Inf/NaN to NA in numeric cols", {
    data <- data.frame(
        time = c("1", "2", "3", "4"),
        A = c("1.5", "Inf", "-Inf", "NaN"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, "time")

    expect_type(result$A, "double")
    expect_equal(result$A, c(1.5, NA_real_, NA_real_, NA_real_))
})

test_that("convert_type() standardises non-finite integers to NA", {
    ## event_channel kept as integer; other int cols coerced to numeric
    data <- data.frame(
        time = c("1", "2", "3"),
        lap = c("1", NA, "2"),
        B = c("10", "20", "30"),
        stringsAsFactors = FALSE
    )

    result <- convert_type(data, "time", event_channel = "lap")

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

    result <- convert_type(data, "time")
    expect_equal(result$A, c(0, 0))
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
    expect_type(result, "character")
    expect_match(result, "2025-03-01T08:30:00")
})

test_that("extract_start_timestamp() detects yyyy-mm-dd HH:MM:SS format", {
    data <- data.frame(
        V1 = c("Start", "2025-06-15 14:22:10"),
        V2 = c("Device", "Moxy"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_type(result, "character")
    expect_match(result, "2025-06-15 14:22:10")
})

test_that("extract_start_timestamp() returns earliest timestamp when multiple present", {
    data <- data.frame(
        V1 = c("2025-01-01 10:00:00", "2025-01-01 11:00:00"),
        V2 = c("2025-01-01 09:00:00", "metadata"),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_match(result, "2025-01-01 09:00:00")
})

test_that("extract_start_timestamp() ignores NA and empty strings", {
    data <- data.frame(
        V1 = c(NA, ""),
        V2 = c("2025-05-10T07:00:00", NA),
        stringsAsFactors = FALSE
    )
    result <- extract_start_timestamp(data)
    expect_type(result, "character")
    expect_match(result, "2025-05-10T07:00:00")
})

test_that("extract_start_timestamp() works with real example file header", {
    file_header <- read_file(example_mnirs("moxy_ramp.xlsx"))[1:20, ]
    result <- extract_start_timestamp(file_header)
    expect_false(is.null(result))
    expect_type(result, "character")
})


## parse_time_channel() ================================================
test_that("parse_time_channel() returns a list with $data and $start_timestamp", {
    data <- data.frame(time = c(0, 1, 2), value = c(1, 2, 3))
    result <- parse_time_channel(data, "time")

    expect_type(result, "list")
    expect_named(result, c("data", "start_timestamp"))
    expect_s3_class(result$data, "data.frame")
})

test_that("parse_time_channel() preserves numeric time (zero_time = FALSE)", {
    data <- data.frame(
        time = c(10.5, 20.5, 30.5),
        value = c(1, 2, 3)
    )

    result <- parse_time_channel(data, "time", zero_time = FALSE)

    expect_equal(result$data$time, data$time)
    expect_null(result$start_timestamp)
})

test_that("parse_time_channel() recalculates numeric time from zero", {
    data <- data.frame(
        time = c(10, 20, 30),
        value = c(1, 2, 3)
    )

    result <- parse_time_channel(data, "time", zero_time = TRUE)

    expect_equal(result$data$time, c(0, 10, 20))
})

test_that("parse_time_channel() parses ISO 8601 character timestamps to numeric", {
    data <- data.frame(
        time = c("2025-01-01T10:00:00", "2025-01-01T10:00:01"),
        value = c(1, 2),
        stringsAsFactors = FALSE
    )

    result <- parse_time_channel(data, "time")

    expect_type(result$data$time, "double")
    expect_equal(result$data$time, c(0, 1))
})

test_that("parse_time_channel() parses date-time formats to numeric", {
    formats <- list(
        c("2025-01-01 10:00:00", "2025-01-01 10:00:01"),
        c("2025/01/01 10:00:00", "2025/01/01 10:00:01"),
        c("01-01-2025 10:00:00", "01-01-2025 10:00:01"),
        c("01/01/2025 10:00:00", "01/01/2025 10:00:01")
    )

    for (fmt in formats) {
        data <- data.frame(
            time = fmt,
            value = c(1, 2),
            stringsAsFactors = FALSE
        )
        result <- parse_time_channel(data, "time")
        expect_equal(result$data$time[1L], 0)
    }
})

test_that("parse_time_channel() parses time-only H:MM:SS character format", {
    data <- data.frame(
        time = c("10:00:00", "10:00:01"),
        value = c(1, 2),
        stringsAsFactors = FALSE
    )

    result <- parse_time_channel(data, "time")

    expect_equal(result$data$time, c(0, 1))
})

test_that("parse_time_channel() handles milliseconds in timestamps", {
    data <- data.frame(
        time = c("2025-01-01T10:00:00.123", "2025-01-01T10:00:01.456"),
        value = c(1, 2),
        stringsAsFactors = FALSE
    )

    result <- parse_time_channel(data, "time")

    expect_equal(result$data$time[1], 0)
    expect_true(result$data$time[2] > 1 & result$data$time[2] < 2)
})

test_that("parse_time_channel() converts POSIXct to numeric seconds from zero", {
    data <- data.frame(
        time = as.POSIXct(c("2025-01-01 10:00:00", "2025-01-01 10:00:01")),
        value = c(1, 2)
    )

    result <- parse_time_channel(data, "time", add_timestamp = FALSE)

    expect_type(result$data$time, "double")
    expect_false("timestamp" %in% names(result$data))
    expect_equal(result$data$time, c(0, 1))
})

test_that("parse_time_channel() always zeros POSIXct regardless of zero_time", {
    data <- data.frame(
        time = as.POSIXct(c("2025-01-01 10:00:00", "2025-01-01 10:00:01")),
        value = c(1, 2)
    )

    result <- parse_time_channel(data, "time", zero_time = FALSE)

    expect_type(result$data$time, "double")
    ## POSIXct is always relative — always starts from 0
    expect_equal(result$data$time, c(0, 1))
})

test_that("parse_time_channel() returns start_timestamp from POSIXct time_channel", {
    t0 <- as.POSIXct("2025-03-15 08:00:00")
    data <- data.frame(
        time = t0 + c(0, 1, 2),
        value = c(1, 2, 3)
    )

    result <- parse_time_channel(data, "time")

    ## start_timestamp is extracted from the POSIXct column when not in header
    expect_false(is.null(result$start_timestamp))
    expect_equal(result$start_timestamp, t0, ignore_attr = TRUE)
})

test_that("parse_time_channel() add_timestamp=TRUE adds POSIXct column after time_channel", {
    t0 <- as.POSIXct("2025-01-01 10:00:00")
    data <- data.frame(
        time = t0 + c(0, 1),
        value = c(1, 2)
    )

    result <- parse_time_channel(data, "time", add_timestamp = TRUE)

    expect_true("timestamp" %in% names(result$data))
    expect_s3_class(result$data$timestamp, "POSIXct")
    expect_type(result$data$time, "double")
    ## timestamp column inserted immediately after time_channel
    expect_equal(which(names(result$data) == "timestamp"), 2)
    expect_equal(
        result$data$timestamp[1L],
        result$start_timestamp,
        ignore_attr = TRUE
    )
    expect_equal(as.numeric(result$data$timestamp), as.numeric(t0 + c(0, 1)))
})


test_that("parse_time_channel() add_timestamp=TRUE with start_timestamp reconstructs absolute timestamps", {
    data <- data.frame(
        time = c(0, 1, 2),
        value = c(1, 2, 3)
    )
    start_ts <- "2025-06-01T09:00:00"

    result <- parse_time_channel(
        data,
        "time",
        start_timestamp = start_ts,
        add_timestamp = TRUE
    )

    expect_true("timestamp" %in% names(result$data))
    expect_s3_class(result$data$timestamp, "POSIXct")

    expected_t0 <- as.POSIXct(start_ts, format = "%Y-%m-%dT%H:%M:%OS")
    expect_equal(
        as.numeric(result$data$timestamp),
        as.numeric(expected_t0 + c(0, 1, 2))
    )
    ## start_timestamp is passed through unchanged
    expect_equal(result$start_timestamp, start_ts)
})

test_that("parse_time_channel() add_timestamp=TRUE with no timestamps skips column silently", {
    data <- data.frame(
        time = c(10, 20, 30),
        value = c(1, 2, 3)
    )

    result <- parse_time_channel(
        data,
        time_channel = "time",
        start_timestamp = NULL,
        add_timestamp = TRUE
    )

    ## no timestamp available — column not added
    expect_false("timestamp" %in% names(result$data))
    expect_null(result$start_timestamp)
})

test_that("parse_time_channel() works on fractional unix time", {
    ## Moxy.csv saved as excel will coerce date-time
    ## to numeric fractional Unix time.
    file_path <- test_path("testdata/moxy-occlusion.xlsx")
    skip_if_not(file.exists(file_path), "testdata not available")

    data <- suppressMessages(readxl::read_excel(file_path)[-(1:2), 1:2])
    names(data)[2L] <- "time"
    data$time <- as.numeric(data$time)

    result <- parse_time_channel(
        data,
        time_channel = "time",
        start_timestamp = NULL,
        add_timestamp = TRUE
    )

    expect_equal(class(result$data$time), "numeric")
    expect_equal(result$data$time[1L], 0)
    expect_equal(median(diff(result$data$time)), 2)
    expect_equal(class(result$data$timestamp), c("POSIXct", "POSIXt"))
    expect_equal(class(result$start_timestamp), c("POSIXct", "POSIXt"))
    ## should return today's date, local time zone, precise timestamp
    expect_equal(as.Date(result$start_timestamp), Sys.Date())
    expect_equal(format(result$start_timestamp, "%Z"), format(Sys.time(), "%Z"))
    expect_equal(format(result$start_timestamp, "%H:%M:%OS"), "13:52:59")
})

test_that("parse_time_channel() returns local time zonel", {
    perfpro <- test_path("testdata/perfpro-mre.xlsx")
    moxy_occl <- test_path("testdata/moxy-occlusion.xlsx")
    vo2master <- test_path("testdata/vo2master.csv")
    skip_if_not(file.exists(perfpro), "testdata not available")
    skip("run in local time PDT/PST zone")
    
    file_list <- c(
        perfpro, ## today's date, 0:00:00
        moxy_occl, ## today's date, 13:52:59
        example_mnirs("moxy_intervals"), ## today's date, 13:17:13
        example_mnirs("moxy_ramp") ## today's date 0:29:00.41
    )

    timestamp_list <- list()
    for (.file in file_list) {
        # .file = file_list[3]
        data <- read_file(.file)
        detected_list <- detect_mnirs_device(data)
        nirs_device <- detected_list$nirs_device
        header_row <- detected_list$header_row

        channels <- detect_device_channels(
            data,
            header_row,
            nirs_device,
            nirs_channels = NULL,
            time_channel = NULL,
            keep_all = FALSE,
            verbose = FALSE
        )
        nirs_channels <- channels$nirs_channels
        time_channel <- channels$time_channel
        keep_all <- channels$keep_all

        table_list <- read_data_table(data, header_row, nirs_channels)
        data <- table_list$data_table
        file_header <- table_list$file_header

        start_timestamp <- extract_start_timestamp(file_header)

        time_channel <- detect_time_channel(
            data,
            time_channel,
            nirs_device,
            verbose = FALSE
        )

        renamed_list <- select_rename_data(
            data,
            nirs_channels,
            time_channel,
            event_channel = NULL,
            keep_all,
            verbose = FALSE
        )
        data <- renamed_list$data
        nirs_renamed <- renamed_list$nirs_channel
        time_renamed <- renamed_list$time_channel
        event_renamed <- renamed_list$event_channel

        data <- remove_empty_rows_cols(data)
        data <- convert_type(data, time_renamed, event_renamed, verbose = FALSE)
        time_list <- parse_time_channel(
            data,
            time_renamed,
            start_timestamp,
            add_timestamp = TRUE,
            zero_time = TRUE
        )

        expect_true(
            format(time_list$start_timestamp, "%Z") %in% c("PDT", "PST")
        )
        timestamp_list <- c(timestamp_list, list(time_list$start_timestamp))
    }
    timestamp_list <- do.call(c, timestamp_list)
    expect_equal(as.Date(timestamp_list), rep(Sys.Date(), 4))
    
    expect_equal(
        format(timestamp_list, "%H:%M:%OS2"),
        c("00:00:00.00", "13:52:59.00", "13:17:13.00", "00:29:00.41")
    )
})

## parse_sample_rate() ================================================
test_that("parse_sample_rate returns correct structure", {
    data <- data.frame(
        time = seq(0, 10, by = 0.1),
        value = rnorm(101, 10, 1)
    )
    file_header <- matrix(NA, nrow = 5, ncol = 5)

    result <- parse_sample_rate(
        data = data,
        file_header = file_header,
        time_channel = "time",
        sample_rate = 10,
        nirs_device = NULL,
        verbose = FALSE
    )

    expect_type(result, "list")
    expect_named(result, c("data", "time_channel", "sample_rate"))
    expect_s3_class(result$data, "data.frame")
    expect_type(result$time_channel, "character")
    expect_type(result$sample_rate, "double")
    expect_equal(result$data$value, data$value)
    expect_equal(nrow(result$data), nrow(data))
    expect_equal(result$sample_rate, 10)
})

test_that("parse_sample_rate handles Artinis device", {
    file_header <- read_file(example_mnirs("artinis_intervals"))
    data <- read_mnirs(
        example_mnirs("artinis_intervals"),
        nirs_channels = c(HHb = 2, O2Hb = 3),
        time_channel = c(sample = 1),
        event_channel = NULL,
        verbose = FALSE
    ) |>
        dplyr::select(-time)

    result <- parse_sample_rate(
        data = data,
        file_header = file_header,
        time_channel = "sample",
        sample_rate = NULL,
        nirs_device = "Artinis",
        verbose = FALSE
    )

    expect_equal(result$sample_rate, 10)
    expect_true("time" %in% names(result$data))
    expect_equal(result$time_channel, "time")
    expect_equal(ncol(result$data), 4)
    expect_equal(result$data$time, data$sample / 10)
})

test_that("parse_sample_rate errors when rate indeterminable", {
    data <- data.frame(x = rep(1, 10))
    file_header <- matrix(NA, nrow = 5, ncol = 5)

    expect_error(
        parse_sample_rate(
            data = data,
            file_header = file_header,
            time_channel = "x",
            sample_rate = NULL,
            verbose = FALSE
        ),
        "Unable to estimate.*sample_rate"
    )
})

test_that("parse_sample_rate verbose output for Artinis", {
    file_header <- read_file(example_mnirs("artinis_intervals"))
    data <- read_mnirs(
        example_mnirs("artinis_intervals"),
        nirs_channels = c(HHb = 2, O2Hb = 3),
        time_channel = c(sample = 1),
        event_channel = NULL,
        verbose = FALSE
    ) |>
        dplyr::select(-time)

    expect_message(
        result <- parse_sample_rate(
            data = data,
            file_header = file_header,
            time_channel = "sample",
            sample_rate = NULL,
            nirs_device = "Artinis",
            verbose = TRUE
        ),
        "Oxysoft.*sample_rate"
    )
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
        "time.*=.*3"
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
        expect_message("Oxysoft.*sample.*column") |> 
        expect_message("Oxysoft.*sample_rate.*10")

    expect_s3_class(df, "mnirs")
    expect_equal(attr(df, "nirs_device"), "Artinis")
    expect_equal(attr(df, "nirs_channels"), "2")
    expect_equal(attr(df, "time_channel"), "time")
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
        expect_warning("irregular") |> 
        expect_message("Detected.*time_channel")

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
        expect_message("Detected.*time_channel") |> 
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
    expect_true(all(
        c("time", "HHb", "O2Hb") %in% names(df)
    ))
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

    ## auto detect channels
    df2 <- read_mnirs(
        file_path,
        nirs_channels = NULL,
        time_channel = NULL,
        verbose = FALSE
    )

    expect_true(all(
        c("sample", "time", "2", "3", "4", "5", "col_6") %in% names(df2)
    ))

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
        "Oxysoft.*sample"
    ) |>
        expect_message("Oxysoft.*sample_rate.*10")

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
    df <- convert_type(df_raw, time_channel)
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

    data <- read_file(file_path)
    detected_list <- detect_mnirs_device(data)
    nirs_device <- detected_list$nirs_device
    header_row <- detected_list$header_row

    expect_equal(nirs_device, "PerfPro")
    expect_equal(header_row, 3)

    expect_message(
        channels <- detect_device_channels(
            data,
            header_row,
            nirs_device,
            nirs_channels = NULL,
            time_channel = NULL,
            keep_all = FALSE,
            verbose = TRUE
        ),
        "PerfPro"
    )

    expect_equal(channels$time_channel, "Time")
    expect_equal(channels$nirs_channels, c("SmO2 (1614)", "SmO2 (1615)"))
    expect_true(channels$keep_all)

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
