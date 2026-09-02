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

    ## import data_raw from pionirs .ftn(2), csv, or excel
    if (grepl("\\.ftn2?$", file_path, ignore.case = TRUE)) {
        ## pionirs tab-separated export: header on row 1, regular columns
        data_raw <- data.table::fread(
            file_path,
            header = FALSE,
            colClasses = "character"
        )
    } else if (grepl("\\.(csv|tsv|txt)$", file_path, ignore.case = TRUE)) {
        ## sample lines for separator and column count detection
        lines <- readLines(file_path, warn = FALSE)
        nrows <- length(lines)

        ## detect separator: comma vs tab from first 10 lines
        head_lines <- lines[seq_len(min(10L, nrows))]
        sep <- if (any(grepl("\t", head_lines, fixed = TRUE))) "\t" else ","

        ## find the max number of separators from the end of the data file
        tail_lines <- lines[seq(to = nrows, by = 1L, len = min(50, nrows))]
        n_seps <- max(lengths(gregexpr(sep, tail_lines, fixed = TRUE)))

        ## pad the first line so fread infers the correct column count
        ## from the data table rather than narrower metadata rows
        data_raw <- data.table::fread(
            text = c(strrep(sep, n_seps), lines),
            header = FALSE,
            fill = Inf,
            sep = sep,
            colClasses = "character",
            strip.white = TRUE
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
            "i" = "Currently only {.var .csv}, {.var .txt}, {.var .xls(x)}, \\
            and {.var .ftn(2)} supported."
        ), call = env)
    }

    return(as.data.frame(data_raw))
}


#' Known channel names and detection patterns for supported mNIRS devices
#'
#' Per device: `pattern` strings which must all match one header row
#' (`fixed` regex flag); default `time_channel` and `event_channel` names;
#' `extra_channels` companion columns (e.g. sample index, numeric tag)
#' returned with `keep_all = TRUE`.
#' @keywords internal
device_patterns <- list(
    Artinis = list(
        ## channels resolved from the legend block by `parse_oxysoft_legend()`
        ## also requires matching "oxysoft" string
        time_channel = NULL,
        pattern = "^(\\d+ )+(\\d+|NA) ?$",
        fixed = FALSE
    ),
    Train.Red = list(
        time_channel = "Timestamp (seconds passed)",
        event_channel = "Lap/Event",
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
    PerfPro = list(
        time_channel = "Time",
        pattern = c("Time.*SmO2"),
        fixed = FALSE
    ),
    PIONIRS = list(
        time_channel = "Time",
        event_channel = "TagLabel",
        extra_channels = c("Iteration", "Tag"),
        pattern = c("Iteration Time", "StO2", "TagLabel"),
        fixed = TRUE
    )
)


#' Detect mnirs device from file metadata
#' @keywords internal
detect_mnirs_device <- function(data, chunk = 200L) {
    ## collapse rows to single strings for pattern matching
    row_strings <- \(.rows) {
        do.call(paste, c(data[.rows, , drop = FALSE], list(sep = " ")))
    }

    ## header rows sit near the top: scan in chunks and stop at the first
    ## match rather than pattern-matching every data row
    n <- nrow(data)
    #! test replace `for` loop
    for (from in seq.int(1L, max(n, 1L), by = chunk)) {
        strings <- row_strings(seq.int(from, min(from + chunk - 1L, n)))

        ## first row index where all of a device's patterns match; NA if none
        first_rows <- vapply(device_patterns, \(.d) {
            hits <- Reduce(`&`, lapply(
                .d$pattern, grepl, x = strings, fixed = .d$fixed
            ))
            which(hits)[1L]
        }, integer(1L))

        if (all(is.na(first_rows))) {
            next
        }

        ## earliest matching row; ties resolved by device order via which.min
        matched_row <- min(first_rows, na.rm = TRUE) + from - 1L
        device_name <- names(first_rows)[which.min(first_rows)]

        ## require "oxysoft" match at or above the row for Artinis pattern
        if (
            identical(device_name, "Artinis") &&
                !any(grepl(
                    "oxysoft",
                    row_strings(seq_len(matched_row)),
                    ignore.case = TRUE
                ))
        ) {
            break
        }

        return(list(nirs_device = device_name, header_row = matched_row))
    }

    return(list(nirs_device = NULL, header_row = 1L))
}


#' Parse channel names from the Oxysoft "Legend" metadata block
#'
#' Legend rows above the numeric header row map column ids to trace names.
#' Returns a channel list of named mappings `c(new_name = "original_col")`,
#' or `NULL` when the legend is missing or malformed.
#' @param raw A raw character data frame from `read_file()`.
#' @param header_row Integer row index of the numeric data table header.
#' @keywords internal
parse_oxysoft_legend <- function(raw, header_row) {
    ## legend entries occupy the first two columns above the header row
    col1 <- trimws(as.character(raw[[1L]][seq_len(header_row - 1L)]))
    legend_row <- which(col1 == "Legend")[1L]
    if (is.na(legend_row) || legend_row + 1L > header_row - 1L) {
        return(NULL)
    }

    ## entry rows: first consecutive run of integer ids after the legend row,
    ## skipping the "Column / Trace (Measurement)" title row
    ids <- suppressWarnings(as.integer(col1))
    candidates <- seq(legend_row + 1L, header_row - 1L)
    runs <- rle(!is.na(ids[candidates]))
    run <- which(runs$values)[1L]
    if (is.na(run)) {
        return(NULL)
    }
    entry_rows <- candidates[
        sum(runs$lengths[seq_len(run - 1L)]) + seq_len(runs$lengths[run])
    ]

    cols <- as.character(ids[entry_rows])
    traces <- trimws(as.character(raw[[2L]][entry_rows]))

    ## legend must be complete and its ids present in the header row
    header_cells <- trimws(as.character(raw[header_row, ]))
    if (any(is_empty(traces)) || !all(cols %in% header_cells)) {
        return(NULL)
    }

    ## build mapping, NULL when no ids match
    named_map <- \(.ids, .names) {
        if (length(.ids) == 0L) {
            return(NULL)
        }
        return(setNames(.ids, rep_len(.names, length(.ids))))
    }

    ## classify traces: parenthesised entries are metadata columns
    paren <- grepl("^\\(.*\\)$", traces)
    time_id <- cols[traces == "(Sample number)"]
    event_id <- cols[traces == "(Event)"]
    extra <- paren & !cols %in% c(time_id, event_id)

    ## unnumbered trailing columns hold event label text; blank header cells
    ## become `col_<idx>` downstream via `rename_duplicates()`
    label_idx <- which(
        is_empty(header_cells) & seq_along(header_cells) > max(ids[entry_rows])
    )

    return(list(
        time = named_map(time_id, "sample"),
        extra = named_map(cols[extra], clean_channel_names(traces[extra])),
        event = named_map(event_id, "event"),
        labels = named_map(
            paste0("col_", label_idx, recycle0 = TRUE),
            "labels"
        ),
        nirs = named_map(cols[!paren], clean_channel_names(traces[!paren]))
    ))
}


#' Extract the export sample rate from Oxysoft header metadata
#' @param header A character data frame of the file rows above the data table.
#' @keywords internal
oxysoft_sample_rate <- function(header) {
    pos <- which(header == "Export sample rate", arr.ind = TRUE)
    return(as.numeric(header[pos[1L], pos[2L] + 1L]))
}


#' Resolve channels from user input, device defaults, or the Oxysoft legend
#'
#' User-specified channels take priority. Otherwise `nirs` channels are read
#' from the Oxysoft legend (Artinis) or header cells starting with `SmO2`;
#' `time` falls back to the device default, and is detected later by
#' `detect_time_channel()` when still `NULL`; `event` falls back to the
#' device default when present in the header row. Device companion columns
#' (`extra`, `labels`) are returned only with `keep_all = TRUE`.
#' @param raw A raw character data frame from `read_file()`.
#' @param device Output of `detect_mnirs_device()`.
#' @param user A list of user-specified `time`, `event`, and `nirs` channels,
#'   each a named character vector `c(new = "original")` or `NULL`.
#' @inheritParams validate_mnirs
#' @returns A list of `time`, `extra`, `event`, `labels`, and `nirs` channel
#'   mappings, each a named `c(new = "original")` vector or `NULL`. List
#'   order sets output column order.
#' @keywords internal
resolve_channels <- function(
    raw,
    device,
    user,
    keep_all = FALSE,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    nirs_device <- device$nirs_device %||% "Unknown"

    auto <- if (identical(nirs_device, "Artinis")) {
        legend <- parse_oxysoft_legend(raw, device$header_row) %||%
            list(time = c(sample = "1"))
        ## user `event_channel = "labels"` aliases the unnumbered label column
        if (
            identical(unname(user$event), "labels") &&
                length(legend$labels) == 1L
        ) {
            user$event[] <- legend$labels
        }
        legend
    } else {
        ## header cells starting with "SmO2" or "StO2", ignoring case; drop
        ## redundant Train.Red unfiltered and Moxy Averaged channels
        header <- Filter(
            Negate(is_empty),
            as.character(raw[device$header_row, ])
        )
        dev <- device_patterns[[nirs_device]]
        list(
            time = dev$time_channel,
            extra = intersect(dev$extra_channels, header),
            event = intersect(dev$event_channel, header),
            nirs = header[
                grepl("^S[mt]O2", header, ignore.case = TRUE) &
                    !grepl("unfiltered|Averaged", header, ignore.case = TRUE)
            ]
        )
    }

    ## drop auto entries for columns the user has already claimed;
    ## empty roles become NULL
    auto <- lapply(auto, \(.x) {
        .x <- .x[!.x %in% unlist(user)]
        if (length(.x) > 0L) name_channels(.x)
    })

    ## user-specified channels take priority over detected
    channels <- list(
        time = user$time %||% auto$time,
        extra = if (keep_all) auto$extra,
        event = user$event %||% auto$event,
        labels = if (keep_all) auto$labels,
        nirs = user$nirs %||% auto$nirs
    )

    if (is.null(channels$nirs)) {
        cli_abort(c(
            "x" = "{.arg nirs_channels} cannot be determined automatically.",
            "i" = "Define {.arg nirs_channels} explicitly."
        ), call = env)
    }

    if (verbose && is.null(user$nirs)) {
        cli_inform(c(
            "!" = "{.val {nirs_device}} file format detected. \\
            {.arg nirs_channels} set to {.field {names(channels$nirs)}}.",
            "i" = "Override by specifying {.arg nirs_channels} explicitly."
        ), call = env)
    }

    return(channels)
}


#' Find the header row containing all `nirs_channels`
#' @param raw A raw character data frame from `read_file()`.
#' @param nirs_channels Character vector of original column names.
#' @param start Integer row index to try first, from `detect_mnirs_device()`.
#' @inheritParams validate_mnirs
#' @keywords internal
find_header_row <- function(
    raw,
    nirs_channels,
    start = 1L,
    env = rlang::caller_env()
) {
    ## match against uniquified names so renamed duplicates (e.g. `SmO2_1`)
    ## and blank columns (`col_5`) resolve as they appear in the data table
    header_row <- Find(\(.i) {
        all(nirs_channels %in% rename_duplicates(as.character(raw[.i, ])))
    }, c(start, seq_len(nrow(raw))))

    if (is.null(header_row)) {
        cli_abort(c(
            "x" = "Channel names not detected.",
            "i" = "Column names are case sensitive and must match exactly."
        ), call = env)
    }

    return(header_row)
}


#' Datetime format strings for POSIXct parsing
#'
#' Time-only format must stay first: `parse_dttm()` treats `dttm_opts[1L]`
#' as a relative time of day and the rest as absolute date-times.
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


#' Detect the first `dttm_opts` format matching a character vector
#'
#' Tested on the first non-empty value only, in UTC (local time zone
#' parsing is slow on Windows).
#' @returns A format string, or `NULL` when none match.
#' @keywords internal
detect_dttm_format <- function(x) {
    x1 <- x[!is_empty(x)][1L]
    if (is.na(x1)) {
        return(NULL)
    }
    return(Find(\(.f) !is.na(strptime(x1, .f, tz = "UTC")), dttm_opts))
}


#' Convert `H:MM(:SS.fff)` strings to seconds of day
#' @keywords internal
hms_to_seconds <- function(x) {
    ## empty strings split to nothing; keep them NA rather than 0
    x[is_empty(x)] <- NA_character_
    parts <- strsplit(x, ":", fixed = TRUE)
    return(suppressWarnings(vapply(parts, \(.p) {
        sum(as.numeric(.p) * c(3600, 60, 1)[seq_along(.p)])
    }, numeric(1L))))
}


#' Parse character date-times with one `dttm_opts` format to local POSIXct
#'
#' Time-only strings are anchored to today's local midnight, matching the
#' Excel fraction-of-day convention.
#' @keywords internal
parse_dttm <- function(x, fmt) {
    if (identical(fmt, dttm_opts[1L])) {
        return(as.POSIXct(format(Sys.Date())) + hms_to_seconds(x))
    }
    return(as.POSIXct(x, format = fmt))
}


#' Extract earliest POSIXct value from file header metadata
#' @keywords internal
extract_start_timestamp <- function(file_header) {
    header_values <- unlist(file_header, use.names = FALSE)
    header_values <- header_values[!is_empty(header_values)]
    ## all dttm_opts contain %H:%M, so candidates must contain ":"
    header_values <- header_values[grepl(":", header_values, fixed = TRUE)]
    if (length(header_values) == 0L) {
        return(NULL)
    }

    ## header formats may differ: screen candidates per format in UTC
    ## (cheap), then parse the hits in local time
    parsed <- unlist(lapply(dttm_opts, \(.f) {
        hits <- header_values[!is.na(strptime(header_values, .f, tz = "UTC"))]
        if (length(hits) > 0L) as.numeric(parse_dttm(hits, .f))
    }))
    parsed <- parsed[!is.na(parsed)]

    if (length(parsed) == 0L) {
        return(NULL)
    }

    return(as.POSIXct(min(parsed)))
}


#' Detect time_channel from column names or time-formatted values
#' @inheritParams validate_mnirs
#' @keywords internal
detect_time_channel <- function(
    data,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    col_names <- names(data)

    ## match column names to possible time column names
    time_regex <- "time|duration|hms|h+:m+:s+"
    time_idx <- grep(time_regex, col_names, ignore.case = TRUE)[1L]

    ## otherwise first column whose first value is a time format string
    if (is.na(time_idx)) {
        time_idx <- Position(\(.col) {
            val <- .col[which(!is_empty(.col))[1L]]
            !is.na(val) && grepl("^\\d{1,2}:\\d{2}(:\\d{2})?", val)
        }, data)
    }

    if (is.na(time_idx)) {
        cli_abort(c(
            "x" = "{.arg time_channel} not detected.",
            "i" = "Define {.arg time_channel} explicitly."
        ), call = env)
    }

    if (verbose) {
        cli_inform(c(
            "!" = "Detected {.arg time_channel} = \\
            {.field {col_names[time_idx]}}."
        ), call = env)
    }

    return(col_names[time_idx])
}


#' Detect empty or NA strings
#' @keywords internal
is_empty <- function(x) {
    is.na(x) | x == ""
}


#' Clean legend trace names to syntactic column names
#' @keywords internal
clean_channel_names <- function(x) {
    ## collapse non-alphanumeric runs to "_", trim edges, lowercase
    x <- gsub("[^[:alnum:]]+", "_", x)
    return(tolower(gsub("^_+|_+$", "", x)))
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
#'
#' Returns a named character vector `c(new = "original")`; unnamed elements
#' are named by their value. `NULL` passes through.
#' @keywords internal
name_channels <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    names <- names(x) %||% character(length(x))
    x <- as.character(x)
    empty_names <- is_empty(names)
    names[empty_names] <- x[empty_names]
    return(setNames(x, names))
}


#' Select, rename, and order channel columns
#'
#' Original names are made unique to match `rename_duplicates(names(data))`;
#' duplicated new names are made unique with a warning. Channel names take
#' priority over clashing names of other data columns. Columns are ordered
#' by role, followed by all remaining columns when `keep_all = TRUE`.
#' @param data The named character data table.
#' @param channels A list of named `c(new = "original")` channel mappings by
#'   role from `resolve_channels()`; `NULL` roles are dropped.
#' @inheritParams validate_mnirs
#' @returns A list of the selected `data` and `channels` as new names by
#'   role.
#' @keywords internal
select_channels <- function(
    data,
    channels,
    keep_all = FALSE,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    channels <- Filter(Negate(is.null), channels)
    roles <- rep(names(channels), lengths(channels))
    orig <- rename_duplicates(unlist(channels, use.names = FALSE))
    new_in <- unlist(lapply(channels, names), use.names = FALSE)
    new <- rename_duplicates(new_in)
    col_idx <- match(orig, names(data))

    if (anyNA(col_idx)) {
        cli_abort(c(
            "x" = "Channel names not detected.",
            "i" = "Column names are case sensitive and must match exactly."
        ), call = env)
    }

    renamed <- new != new_in
    if (verbose && any(renamed)) {
        cli_warn(c(
            "!" = "Duplicate channel names detected.",
            "i" = "Renamed: {.field {paste(new_in[renamed], new[renamed], sep = ' = ')}}",
            "i" = "Unique channel names can be defined explicitly."
        ), call = warn_call(env))
    }

    ## rename remaining columns first so channel names win any clash
    names(data) <- rename_duplicates(c(new, names(data)))[-seq_along(new)]
    names(data)[col_idx] <- new
    data <- data[c(col_idx, if (keep_all) setdiff(seq_along(data), col_idx))]

    return(list(
        data = data,
        channels = split(new, factor(roles, levels = names(channels)))
    ))
}


#' Remove Empty Rows and Columns
#' @keywords internal
remove_empty_rows_cols <- function(data) {
    filled <- !is_empty(as.matrix(data))
    return(data[rowSums(filled) > 0, colSums(filled) > 0, drop = FALSE])
}


#' Coerce column types by role: nirs numeric, event integer, others detected
#' @param channels A list of `time`, `event`, and `nirs` column names.
#' @inheritParams validate_mnirs
#' @keywords internal
convert_type <- function(
    data,
    channels,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## convert decimal "," to "."
    is_char <- vapply(data, is.character, logical(1L))
    char_cols <- setdiff(names(data)[is_char], channels$time)
    data[char_cols] <- lapply(
        data[char_cols], gsub, pattern = ",", replacement = ".", fixed = TRUE
    )

    ## coerce by role, then standardise NA once:
    ## - time           left unchanged for parse_time_channel
    ## - nirs           -> numeric
    ## - event          -> integer or character
    ## - other columns  -> integer to numeric, otherwise unopinionated
    data[] <- Map(\(.x, .nm) {
        if (.nm %in% channels$event) {
            ## "-" is the pionirs placeholder for no event
            .x <- utils::type.convert(
                .x, na.strings = c("NA", "", "-"), as.is = TRUE
            )
        } else if (.nm %in% channels$nirs) {
            .x <- suppressWarnings(as.numeric(.x))
        } else if (!.nm %in% channels$time) {
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
        nirs_cols <- intersect(channels$nirs, names(data))
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


#' Parse time_channel character or dttm to numeric seconds
#'
#' @param x The time column vector: numeric, character, or POSIXct.
#' @param start_timestamp Optional POSIXct from the file header, evaluated
#'   lazily only when `x` is not an absolute date-time series.
#' @param zero_time Logical; re-base numeric time to start from zero.
#' @returns A list of `time` (numeric seconds), `timestamp` (POSIXct
#'   vector or `NULL`), and `start_timestamp` (POSIXct or `NULL`).
#' @keywords internal
parse_time_channel <- function(x, start_timestamp = NULL, zero_time = FALSE) {
    dated <- inherits(x, "POSIXct")

    ## character time -> numeric (sample index / seconds) where numeric-like,
    ## otherwise parse the detected date-time format to POSIXct; a time-only
    ## series is relative and must be anchored by a header timestamp
    if (is.character(x)) {
        num <- suppressWarnings(as.numeric(x))
        fmt <- if (all(is.na(num))) detect_dttm_format(x)
        if (is.null(fmt)) {
            x <- num
        } else {
            x <- parse_dttm(x, fmt)
            dated <- !identical(fmt, dttm_opts[1L])
        }
    }

    ## excel fraction-of-day time to POSIXct from today's local midnight
    if (is.numeric(x) && all(x >= 0 & x <= 1, na.rm = TRUE)) {
        x <- as.POSIXct(format(Sys.Date())) + x * 86400
    }

    ## recalculate numeric time to start from zero
    if (zero_time && is.numeric(x)) {
        x <- x - x[1L]
    }

    ## preserve POSIXct timestamp and convert to numeric seconds
    timestamp <- NULL
    if (inherits(x, "POSIXct")) {
        timestamp <- x
        x <- as.numeric(difftime(x, x[1L], units = "secs"))
    }

    ## a dated series anchors itself; `!dated` short-circuits so the lazy
    ## header argument is never forced. Otherwise a header timestamp anchors
    ## the relative series, else fall back to the parsed series.
    ## first sample, not earliest, so `start_timestamp + time` matches the
    ## zeroed time when samples are out of order
    if (!dated && !is.null(start_timestamp)) {
        timestamp <- start_timestamp + x
    } else {
        start_timestamp <- timestamp[1L]
    }

    return(list(
        time = x,
        timestamp = timestamp,
        start_timestamp = start_timestamp
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
