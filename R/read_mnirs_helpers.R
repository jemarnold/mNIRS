#' Read raw data frame from file path
#' @inheritParams validate_mnirs
#' @keywords internal
read_file <- function(file_path, env = rlang::caller_env()) {
    ## validation: check file exists
    if (!file.exists(file_path)) {
        cli_abort(c(
            "x" = "File not found. Check that file exists.",
            "i" = "{.arg file_path} = {.path {file_path}}"
        ), call = env)
    }

    ## import data_raw from either excel or csv
    if (grepl("\\.(csv|tsv|txt)$", file_path, ignore.case = TRUE)) {
        ## sample lines for separator and column count detection
        ## strip whitespace inside quoted fields and around line edges
        lines <- trimws(
            gsub('\\s*"\\s*', '"', readLines(file_path, warn = FALSE))
        )
        nrows <- length(lines)

        ## detect separator: comma vs tab from first 10 lines
        head_lines <- lines[seq_len(min(10L, nrows))]
        sep <- if (any(grepl("\t", head_lines))) "\t" else ","

        ## find the max number of separators from the end of the data file
        tail_lines <- lines[seq(to = nrows, by = 1L, len = min(50, nrows))]
        n_seps <- max(lengths(gregexpr(sep, tail_lines, fixed = TRUE)))

        ## pad the first line so fread infers the correct column count
        lines <- c(strrep(sep, n_seps), lines)

        ## read with explicit sep and column count to handle
        ## irregular header rows with fewer columns than data
        data_raw <- data.table::fread(
            text = lines,
            header = FALSE,
            fill = Inf,
            sep = sep,
            colClasses = "character",
        )[-1, ]
    } else if (grepl("\\.xls(x)?$", file_path, ignore.case = TRUE)) {
        ## report error when file is open and cannot be accessed by readxl
        data_raw <- tryCatch(
            readxl::read_excel(
                path = file_path,
                col_names = FALSE,
                col_types = "text",
                .name_repair = "minimal"
            ),
            error = \(e) {
                if (grepl("cannot be opened", e$message)) {
                    cli_abort(c(
                        "x" = "File cannot be opened.",
                        "i" = "Check file is not in use by another \\
                        application.",
                        "i" = "{e$message}"
                    ), call = env)
                } else {
                    cli_abort(e$message, call = env)
                }
            }
        )
    } else {
        ## validation: check file types
        cli_abort(c(
            "x" = "Unsupported file type.",
            "i" = "{.arg file_path} = {.path {file_path}}",
            "i" = "Currently only {.var .csv}, {.var .txt}, and {.var .xls(x)} supported."
        ), call = env)
    }

    return(data.frame(data_raw))
}


#' Known channel names and detection patterns for supported mNIRS devices
#' @keywords internal
device_patterns <- list(
    Artinis = list(
        ## keep `time_channel = NULL` to force `detect_time_channel` message
        ## also requires matching "oxysoft" string
        time_channel = NULL,
        pattern = "^(\\d+ )+(\\d+|NA) ?$",
        fixed = FALSE
    ),
    Train.Red = list(
        time_channel = "Timestamp (seconds passed)",
        pattern = c("Timestamp (seconds passed)", "SmO2"),
        fixed = TRUE
    ),
    Moxy = list(
        time_channel = "hh:mm:ss",
        pattern = c("hh:mm:ss", "SmO2 Live"),
        fixed = TRUE
    ),
    VO2master = list(
        time_channel = "Time[s]",
        pattern = c("Time[s]", "SmO2[%]"),
        fixed = TRUE
    ),
    ## matches `nirs_channels` with "SmO2" below
    PerfPro = list(
        time_channel = "Time",
        pattern = c("Time.*SmO2"),
        fixed = FALSE
    )
)


#' Detect mnirs device from file metadata
#' @keywords internal
detect_mnirs_device <- function(data) {
    ## collapse each row to a single string for pattern matching
    data_strings <- do.call(paste, c(data, list(sep = " ")))

    ## first row index where all of a device's patterns match; NA if none
    first_rows <- vapply(device_patterns, \(.d) {
        hits <- Reduce(`&`, lapply(
            .d$pattern, grepl, x = data_strings, fixed = .d$fixed
        ))
        which(hits)[1L]
    }, integer(1L))

    if (all(is.na(first_rows))) {
        return(list(nirs_device = NULL, header_row = 1L))
    }

    ## earliest matching row; ties resolved by device order via which.min
    matched_row <- min(first_rows, na.rm = TRUE)
    device_name <- names(first_rows)[which.min(first_rows)]

    ## require "oxysoft" match for Artinis pattern
    if (identical(device_name, "Artinis")) {
        above_strings <- data_strings[seq_len(matched_row)]
        if (!any(grepl("oxysoft", above_strings, ignore.case = TRUE))) {
            return(list(nirs_device = NULL, header_row = 1L))
        }
    }

    return(list(nirs_device = device_name, header_row = matched_row))
}


#' Detect known channels for a device
#' @inheritParams validate_mnirs
#' @keywords internal
detect_device_channels <- function(
    data,
    header_row = 1L,
    nirs_device = NULL,
    nirs_channels = NULL,
    time_channel = NULL,
    keep_all = FALSE,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## user-specified channels always take priority
    if (!is.null(nirs_channels)) {
        return(list(
            ## if `time_channel = NULL` defined at `detect_time_channel`
            time_channel = time_channel,
            nirs_channels = nirs_channels,
            keep_all = keep_all
        ))
    }

    ## scan header row for cols starting with "SmO2" ignore case, ignore NA
    header <- Filter(Negate(is_empty), as.character(data[header_row, ]))

    nirs_channels <- if (identical(nirs_device, "Artinis")) {
        ## Artinis Oxysoft exports use numeric channel ids, not "SmO2"
        header[startsWith(header, "2")]
    } else {
        header[startsWith(toupper(header), "SMO2")]
    }

    ## drop redundant unfiltered Train.Red, and Averaged Moxy channels
    nirs_channels <- nirs_channels[
        !grepl("unfiltered|Averaged", nirs_channels, ignore.case = TRUE)
    ]

    if (length(nirs_channels) == 0L) {
        cli_abort(c(
            "x" = "{.arg nirs_channels} cannot be determined automatically.",
            "i" = "Define {.arg nirs_channels} explicitly."
        ), call = env)
    }

    ## check for NULL
    nirs_device <- nirs_device %||% "Unknown"

    ch_list <- list(
        ## priority: user `time_channel` -> device default -> NULL
        time_channel = time_channel %||%
            device_patterns[[nirs_device]]$time_channel,
        nirs_channels = nirs_channels,
        keep_all = TRUE ## return all cols to view potential nirs_channels
    )

    if (verbose) {
        cli_inform(c(
            "!" = "{.val {nirs_device}} file format detected. \\
            {.arg nirs_channels} set to {.field {ch_list$nirs_channels}}.",
            "i" = "Override by specifying {.arg nirs_channels} explicitly."
        ), call = env)
    }

    return(ch_list)
}


#' Read data table from raw data
#' @inheritParams validate_mnirs
#' @keywords internal
read_data_table <- function(
    data,
    header_row = 1L,
    nirs_channels = NULL,
    env = rlang::caller_env()
) {
    nrows <- nrow(data)
    ## find the first row where ALL nirs_channels match
    ## start with `header_row` passed from `detect_mnirs_device()`
    header_row <- Find(\(.i) {
        all(nirs_channels %in% data[.i, ])
    }, c(header_row, seq_len(nrows)))

    ## validation: all channels must be detected to extract the data frame
    ## return error if channels string is detected at multiple rows
    if (length(header_row) == 0) {
        cli_abort(c(
            "x" = "Channel names not detected.",
            "i" = "Column names are case sensitive and must match exactly."
        ), call = env)
    }

    ## extract the data_table, and name by header row
    table_rows <- (header_row + 1L):nrows
    data_table <- setNames(data[table_rows, ], data[header_row, ])
    file_header <- data[seq_len(header_row), ]

    return(list(
        file_header = file_header,
        data_table = data_table
    ))
}


#' Datetime format strings for POSIXct parsing
#' @keywords internal
dttm_opts <- c(
    "%H:%M:%OS",
    "%Y-%m-%dT%H:%M:%OS",
    "%Y-%m-%dT%H:%M:%OS%z",
    "%Y-%m-%d %H:%M:%OS",
    "%Y/%m/%d %H:%M:%OS",
    "%d-%m-%Y %H:%M:%OS",
    "%d/%m/%Y %H:%M:%OS"
)


#' Extract earliest POSIXct value from file header metadata
#' @keywords internal
extract_start_timestamp <- function(file_header) {
    header_values <- unlist(file_header, use.names = FALSE)
    header_values <- header_values[!is_empty(header_values)]
    ## all dttm_opts contain %H:%M, so candidates must contain ":"
    header_values <- header_values[grepl(":", header_values, fixed = TRUE)]

    ## search for POSIXct values, return the earliest time value
    ## vulnerable to invalid timestamps
    parsed <- which(
        !is.na(vapply(header_values, \(.x) {
        as.POSIXct(.x, tryFormats = dttm_opts, optional = TRUE)
    }, numeric(1L)))
    )

    if (length(parsed) == 0L) {
        return(NULL)
    }

    ## return the earliest character string timestamp, assuming start time
    return(min(header_values[parsed]))
}


#' Detect time_channel from header row
#' @inheritParams validate_mnirs
#' @keywords internal
detect_time_channel <- function(
    data,
    time_channel = NULL,
    nirs_device = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    if (!is.null(time_channel)) {
        return(time_channel)
    }

    ## return default sample column for Artinis Oxysoft
    ## when `nirs_channels` defined but `time_channel = NULL`
    if (!is.null(nirs_device) && nirs_device == "Artinis") {
        if (verbose) {
            cli_inform(c(
                "!" = "Oxysoft {.val sample} column detected."
            ), call = env)
        }
        return(c(sample = "1"))
    }

    col_names <- names(data)

    ## match column names to possible time column names
    time_regex <- "time|duration|hms|h+:m+:s+"
    time_idx <- grep(time_regex, col_names, ignore.case = TRUE)[1L]

    ## find name of POSIXct column
    if (is.na(time_idx)) {
        time_idx <- Position(\(.col) inherits(.col, "POSIXct"), data)
    }

    ## find name of character column with time format strings
    if (is.na(time_idx)) {
        time_idx <- Position(\(.col) {
            is.character(.col) && {
                val <- .col[which(!is.na(.col))[1L]]
                !is.na(val) && grepl("^\\d{1,2}:\\d{2}(:\\d{2})?", val)
            }
        }, data)
    }

    if (!is.na(time_idx)) {
        if (verbose) {
            cli_inform(c(
                "!" = "Detected {.arg time_channel} = \\
                {.field {col_names[time_idx]}}."
            ), call = env)
        }
        return(col_names[time_idx])
    }

    cli_abort(c(
        "x" = "{.arg time_channel} not detected.",
        "i" = "Define {.arg time_channel} explicitly."
    ), call = env)
}


#' Detect empty or NA strings
#' @keywords internal
is_empty <- function(x) {
    is.na(x) | x == ""
}


#' Rename duplicate strings in a vector with `make.unique()`
#' @keywords internal
rename_duplicates <- function(x) {
    if (is.null(x)) {
        return(x)
    }
    ## rename blank strings
    empty <- which(is_empty(x))
    x[empty] <- paste0("col_", empty)

    return(make.unique(x, sep = "_"))
}


#' Force names on character strings
#' @keywords internal
name_channels <- function(x) {
    ## if channels not named, names from object
    names <- names(x) %||% character(length(x))
    empty_names <- is_empty(names)
    names[empty_names] <- as.character(x)[empty_names]
    return(setNames(x, names))
}


#' Select data table columns and rename from channels, handling duplicates
#' @inheritParams validate_mnirs
#' @keywords internal
select_rename_data <- function(
    data,
    nirs_channels,
    time_channel,
    event_channel = NULL,
    keep_all = FALSE,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## ensure all channel inputs are named (name = original_col_name)
    ch_list <- list(
        time_channel = time_channel,
        event_channel = event_channel,
        nirs_channels = nirs_channels
    ) |>
        lapply(\(.x) if (is.null(.x)) .x else name_channels(.x))

    ## original column names (values) mapped to user names (names)
    ## rename_duplicates makes user-facing names unique
    orig_names <- lapply(ch_list, \(.x) {
        if (is.null(.x)) NULL else rename_duplicates(as.character(.x))
    })
    user_names <- lapply(ch_list, \(.x) {
        if (is.null(.x)) NULL else rename_duplicates(names(.x))
    })

    ## flat vectors for column matching
    orig_vec <- unlist(orig_names, use.names = FALSE)
    user_vec <- unlist(user_names, use.names = FALSE)

    ## de-duplicate data column names
    data_names <- rename_duplicates(names(data))

    ## check channels exist in data
    missing <- setdiff(orig_vec, data_names)
    if (length(missing) > 0L) {
        cli_abort(c(
            "x" = "Channel names not detected.",
            "i" = "Column names are case sensitive and must match exactly."
        ), call = env)
    }

    ## keep all columns or only specified channels
    selected_cols <- if (keep_all) {
        c(orig_vec, setdiff(data_names, orig_vec))
    } else {
        orig_vec
    }

    ## select and rename: channel columns get user names,
    ## remaining columns keep de-duplicated data names
    result <- setNames(data, data_names)[, selected_cols, drop = FALSE]
    channel_idx <- match(orig_vec, names(result))

    ## resolve clashes: user names take priority over data names
    all_names <- rename_duplicates(c(user_vec, names(result)))
    names(result) <- all_names[!all_names %in% user_vec]
    names(result)[channel_idx] <- user_vec

    ## warn if any channels were renamed from their input names
    renamed <- user_vec != unlist(lapply(ch_list, names), use.names = FALSE)
    if (verbose && any(renamed)) {
        old <- unlist(lapply(ch_list, names), use.names = FALSE)[renamed]
        new <- user_vec[renamed]
        cli_warn(c(
            "!" = "Duplicate channel names detected.",
            "i" = "Renamed: {.field {paste(old, new, sep = ' = ')}}",
            "i" = "Unique channel names can be defined explicitly."
        ), call = warn_call(env))
    }

    return(list(
        data = result,
        nirs_channel = user_names$nirs_channels,
        time_channel = user_names$time_channel,
        event_channel = user_names$event_channel
    ))
}


#' Remove Empty Rows and Columns
#' @keywords internal
remove_empty_rows_cols <- function(data) {
    data <- data[rowSums(!is_empty(data)) > 0, , drop = FALSE]
    return(data[, colSums(!is_empty(data)) > 0, drop = FALSE])
}


#' Coerce column types by role: nirs numeric, event integer, others detected
#' @inheritParams validate_mnirs
#' @keywords internal
convert_type <- function(
    data,
    nirs_channels = NULL,
    time_channel,
    event_channel = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## convert decimal "," to "."
    is_char <- vapply(data, is.character, logical(1L))
    char_cols <- setdiff(names(data)[is_char], time_channel)
    data[char_cols] <- lapply(
        data[char_cols], gsub, pattern = ",", replacement = ".", fixed = TRUE
    )

    ## coerce by role, then standardise NA once:
    ## - time_channel left unchanged for parse_time_channel
    ## - nirs_channels  -> numeric
    ## - event_channel  -> integer or character
    ## - other columns  -> integer to numeric, otherwise unopinionated
    data[] <- Map(\(.x, .nm) {
        if (.nm %in% event_channel) {
            .x <- utils::type.convert(
                .x, na.strings = c("NA", ""), as.is = TRUE
            )
        } else if (.nm %in% nirs_channels) {
            .x <- suppressWarnings(as.numeric(.x))
        } else if (!.nm %in% time_channel) {
            .x <- utils::type.convert(
                .x, na.strings = c("NA", ""), as.is = TRUE
            )
            if (is.integer(.x)) {
                .x <- as.numeric(.x)
            }
        }
        ## standardise empty/NaN/Inf to NA by resulting type
        if (is.character(.x)) {
            .x[.x %in% c("", "NA")] <- NA_character_
        } else {
            .x[!is.finite(.x)] <- NA
        }
        .x
    }, data, names(data))

    ## warn per channel when nirs values are all coerced to NA
    if (verbose) {
        nirs_cols <- intersect(nirs_channels, names(data))
        all_na <- vapply(data[nirs_cols], \(.x) all(is.na(.x)), logical(1L))
        lapply(nirs_cols[all_na], \(.nm) {
            cli_warn(c(
                "!" = "Channel {.field {(.nm)}} values coerced to {.val {NA}}.",
                "i" = "Check the source data contents should be numeric values."
            ), call = warn_call(env))
        })
    }

    return(data)
}


#' Parse time_channel character or dttm to numeric
#' @keywords internal
parse_time_channel <- function(
    data,
    time_channel,
    start_timestamp = NULL,
    add_timestamp = FALSE,
    zero_time = FALSE
) {
    t_vec <- data[[time_channel]]

    ## character time -> numeric (sample index / seconds) where numeric-like,
    ## otherwise parse date-time formats to POSIXct
    if (is.character(t_vec)) {
        num <- suppressWarnings(as.numeric(t_vec))
        t_vec <- if (any(!is.na(num))) {
            num
        } else {
            as.POSIXct(t_vec, tryFormats = dttm_opts, optional = TRUE)
        }
    }

    ## fractional unix time to POSIXct coerced to local time zone
    if (is.numeric(t_vec) && all(t_vec <= 1, na.rm = TRUE)) {
        t_vec <- as.POSIXct(
            as.character(as.POSIXct(Sys.Date(), "UTC")),
            tz = Sys.timezone()
        ) +
            t_vec * 86400
    }

    ## recalculate numeric time to start from zero
    if (zero_time && is.numeric(t_vec)) {
        t_vec <- t_vec - t_vec[1L]
    }

    ## preserve POSIXct timestamp and convert to numeric seconds
    timestamp_vec <- NULL
    if (inherits(t_vec, "POSIXct")) {
        timestamp_vec <- t_vec
        t_vec <- as.numeric(difftime(t_vec, t_vec[1L], units = "secs"))
    }

    data[[time_channel]] <- t_vec

    ## add_timestamp preserves or adds POSIXct/dttm column
    if (add_timestamp) {
        ## add "timestamp" col after `time_channel` position
        col_names <- names(data)
        time_idx <- match(time_channel, col_names)
        data_names <- append(col_names, "timestamp", time_idx)
        data$timestamp <- NA_real_
        data <- data[data_names]

        ## if neither header start_timestamp or timestamp_vec exist
        ## then return NULL and don't append column
        if (!is.null(start_timestamp)) {
            start_time <- as.POSIXct(
                start_timestamp,
                tryFormats = dttm_opts,
                optional = TRUE
            )
            data$timestamp <- start_time + t_vec
        } else if (!is.null(timestamp_vec)) {
            data$timestamp <- timestamp_vec
        } else {
            ## column removed if no timestamp detected
            data$timestamp <- NULL
        }
    }

    ## extract earliest POSIXct value as start_timestamp metadata
    ## if not already passed from file header
    if (is.null(start_timestamp) && !is.null(timestamp_vec)) {
        start_timestamp <- min(timestamp_vec, na.rm = TRUE)
    }

    return(list(data = data, start_timestamp = start_timestamp))
}


#' Validate and Estimate Sample Rate
#' @inheritParams validate_mnirs
#' @keywords internal
parse_sample_rate <- function(
    data,
    file_header,
    time_channel,
    sample_rate = NULL,
    nirs_device = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## if Oxysoft, sample_rate will be detected = 1
    ## extract and overwrite with exported sample_rate
    ## create new "time" column at col_idx behind `time_channel`
    if (!is.null(nirs_device) && nirs_device == "Artinis") {
        pos <- which(file_header == "Export sample rate", arr.ind = TRUE)
        sample_rate <- as.numeric(file_header[pos[1L], pos[2L] + 1L])

        col_names <- names(data)
        time_idx <- match(time_channel, col_names)
        data_names <- append(col_names, "time", time_idx)
        data_names <- rename_duplicates(data_names)
        t_vec <- data[[time_channel]] / sample_rate
        time_channel <- setdiff(data_names, col_names)
        data[[time_channel]] <- t_vec
        data <- data[data_names]

        if (verbose) {
            cli_inform(c(
                "!" = "Oxysoft {.arg sample_rate} = {.val {sample_rate}} Hz.",
                "i" = "{.arg time_channel} = {.field {time_channel}} added to \\
                the data frame, in {.cls seconds}."
            ), call = env)
        }
    }

    ## validate priority user input sample_rate
    ## metadata check will be skipped
    ## will estimate from time_channel (time_channel)
    ## will error on unable to estimate sample_rate
    sample_rate <- validate_sample_rate(
        data, time_channel, sample_rate, verbose, env = env
    )

    return(list(
        data = data,
        time_channel = time_channel,
        sample_rate = sample_rate
    ))
}


#' Report warnings for unbalanced time_channel samples
#' @inheritParams validate_mnirs
#' @keywords internal
detect_irregular_samples <- function(
    x,
    time_channel,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    if (!verbose) {
        return(invisible())
    }

    ## flag duplicated, unordered, or big-gap samples
    diffs <- diff(x)
    irregular <- duplicated(x)
    irregular[-1L] <- irregular[-1L] | diffs < 0 | diffs >= 3600

    ## skip if no irregular samples
    if (!any(irregular)) {
        return(invisible())
    }

    irregular_vec <- unique(round(x[irregular], 6))
    n_total <- length(irregular_vec)

    info_msg <- if (n_total > 5L) {
        ## more than 5: print the first three plus remaining count
        "{.field {time_channel}} = {.val {irregular_vec[seq_len(3L)]}} \\
        and {n_total - 3L} more."
    } else {
        ## 5 or fewer: print all
        "{.field {time_channel}} = {.val {irregular_vec}}."
    }

    cli_warn(c(
        "!" = "Duplicate or irregular {.arg time_channel} samples detected.",
        "i" = info_msg,
        "i" = "Re-sample with {.fn mnirs::resample_mnirs}."
    ), call = warn_call(env))

    return(invisible())
}
