#' Read *{mnirs}* data from file
#'
#' Import time-series data exported from common muscle NIRS (mNIRS) devices and
#' return a tibble of class `"mnirs"` with the selected signal channels and
#' metadata.
#'
#' @param file_path Path of the data file to import. Supported file extensions
#'   include `".xls(x)"`, `".csv"`, `".txt"`, and *PIONIRS* `".ftn(2)"`.
#'
#' @param nirs_channels A character vector of one or more column names
#'   containing mNIRS signals to import. Names must match the file header
#'   exactly.
#'
#'   - If `NULL` (default), `read_mnirs()` attempts to detect the device from
#'     the file contents and use a known `nirs_channel` name.
#'   - A *named* character vector can be used to rename columns on import, in
#'     the form `c(renamed = "original_name")`.
#'
#' @param time_channel A character string giving the name of the time
#'   (or sample) column to import. The name must match the file header exactly.
#'
#'   - If `NULL` (default), `read_mnirs()` attempts to identify a time-like
#'     column automatically (by known device defaults and/or time-formatted
#'     values).
#'   - A *named* character vector can be used to rename the column on import,
#'     in the form `c(time = "original_name")`.
#'
#' @param event_channel An *optional* character string giving the name of an
#'   event/lap column to import. Names must match the file header exactly.
#'   A named character vector can be used to rename the column on import in
#'   the form `c(event = "original_name")`.
#'
#' @param sample_rate An *optional* numeric sample rate in Hz. If left blank
#'   (`NULL`), the sample rate is estimated from `time_channel` (see *Details*).
#'
#' @param add_timestamp A logical. Default is `FALSE`. If `TRUE` and if the
#'   source data contain an absolute date-time (POSIXct) time value, will add
#'   a `"timestamp"` column in addition to the specified `time_channel` as a
#'   numeric time column.
#'
#' @param zero_time Logical. Default is `FALSE`. If `TRUE`, re-calculates
#'   numeric `time_channel` values to start from zero.
#'
#' @param keep_all Logical. Default is `FALSE`. Will keep only the channels
#'   explicitly specified in `nirs_channels`, `time_channel`, and
#'   `event_channel`. If `TRUE` will keep all columns found in the file
#'   data table.
#'
#'   - If no `nirs_channels` are specified and the file format is recognised,
#'     all columns in the file data table will be returned, as an exploratory
#'     option.
#'
#' @param verbose Logical. Default is `TRUE`. Display or silence (if `FALSE`)
#'   warnings and information messages helpful for troubleshooting. Ad
#'   global default can be set via `options(mnirs.verbose = FALSE)`.
#'
#' @details
#' ## Header detection
#' `read_mnirs()` searches the file for a header row containing the requested
#' channel names. The header row does not need to be the first row in the file.
#'
#' - If duplicate column names exist (e.g. for *Train.Red* files), they are
#'   made unique with a numbered suffix (`_*`), and can be renamed accordingly:
#'   `nirs_channels = c(smo2_left = "smo2", smo2_right = "smo2_1")`.
#' - Unnamed columns in the source file will be renamed to `col_*`, where `*`
#'   is the ordered column number in the file (e.g. `col_6`; *Artinis Oxysoft*
#'   files have an exception to this renaming convention. See below).
#'
#' ## Artinis Oxysoft exports
#' *Artinis Oxysoft* files have numbered data columns and channel names in a
#' "Legend" metadata block. `read_mnirs()` can detect and rename columns
#' automatically:
#'
#' - `nirs_channels` names become clean lower-case column names with
#'   underscores (e.g. `"Rx1 - Tx1 O2Hb"` becomes `rx1_tx1_o2hb`).
#' - `"(Sample number)"` column is renamed `sample`, and a `time` column
#'   in seconds is derived from the export sample rate (see *Time parsing*).
#' - `"(Event)"` column is renamed `event` and set as `event_channel`.
#' - *Oxysoft* exports a trailing un-numbered column containing optional
#'   event label text. This is renamed `labels` if present, or dropped when
#'   empty. It can be selected as the event column explicitly with
#'   `event_channel = c(event = "labels")`.
#'
#' Explicit `nirs_channels`, `time_channel`, and `event_channel` renaming
#' (as below) overrides the automatic detection names.
#'
#' ## Renaming channels
#' All `channels` can be renamed by specifying a named character vector in
#' the form `c(renamed = "original_name")`. The `"original_name"` must match
#' the file header row exactly.
#'
#' ## Time parsing
#' If `time_channel` is `NULL`, it is resolved from the known device default,
#' or by detecting a time-like column name (e.g. `"time"`, `"hh:mm:ss"`), or
#' by detecting a column with time-formatted (POSIXct-like) values.
#'
#' - If `time_channel` is a date-time (POSIXct) format, it will be converted
#'   to numeric and re-based to start from 0, regardless of `zero_time`.
#' - Some devices export a sample index rather than time values. In those
#'   cases, if an export `sample_rate` is detected in the file metadata (e.g.
#'   *Artinis Oxysoft* exports), `read_mnirs()` will create a *"time"* column
#'   in seconds derived from the sample index and the detected `sample_rate`.
#'
#' ## Sample rate
#' If `sample_rate` is not specified, it is estimated from differences in
#' `time_channel`. If `time_channel` is actually a sample index, as described
#' above, this may erroneously be estimated at 1 Hz. `sample_rate` should be
#' specified explicitly in this case.
#'
#' ## Data cleaning
#' Entirely empty rows and columns are removed. Invalid values (e.g.
#' `c(NaN, Inf)`) are standardized to `NA`. A warning will be displayed when
#' irregular sampling is detected (e.g. non-monotonic, repeated, or unequal
#' time values), if `verbose = TRUE`. In this case, it is recommended to use
#' `resample_mnirs()` to standardise the time grid to the desired
#' `sample_rate`.
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class `"mnirs"`. Metadata are stored
#'   as attributes and can be accessed with `attributes(data)`.
#'
#' @examples
#' read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"), ## call an example data file
#'     nirs_channels = c(
#'         smo2_left = "SmO2 Live",            ## identify and rename channels
#'         smo2_right = "SmO2 Live(2)"
#'     ),
#'     time_channel = c(time = "hh:mm:ss"),    ## date-time format will be converted to numeric
#'     event_channel = NULL,                   ## leave blank if unused
#'     sample_rate = NULL,                     ## if blank, will be estimated from time_channel
#'     add_timestamp = FALSE,                  ## omit a date-time timestamp column
#'     zero_time = TRUE,                       ## recalculate time values from zero
#'     keep_all = FALSE,                       ## return only the specified data channels
#'     verbose = TRUE                          ## show warnings & messages
#' )
#'
#' @export
read_mnirs <- function(
    file_path,
    nirs_channels = NULL,
    time_channel = NULL,
    event_channel = NULL,
    sample_rate = NULL,
    add_timestamp = FALSE,
    zero_time = FALSE,
    keep_all = FALSE,
    verbose = TRUE
) {
    ## global options overrides implicit but not explicit `verbose`
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }

    ## import raw character table & detect device
    raw <- read_file(file_path)
    device <- detect_mnirs_device(raw)
    nirs_device <- device$nirs_device

    ## user channels as `c(new = "original")`
    user <- lapply(
        list(time = time_channel, event = event_channel, nirs = nirs_channels),
        name_channels
    )

    ## `keep_all` coerced to `TRUE` when `nirs_channels` are auto-detected,
    ## otherwise respect user supplied `keep_all`
    keep_all <- keep_all || is.null(nirs_channels)

    ## resolve channels from user input, device defaults, or Oxysoft legend
    channels <- resolve_channels(raw, device, user, keep_all, verbose)
    header_row <- find_header_row(raw, channels$nirs, device$header_row)

    ## name data table by header row
    ## blank/duplicate names made unique
    data <- setNames(
        raw[-seq_len(header_row), ],
        rename_duplicates(as.character(raw[header_row, ]))
    )
    channels$time <- channels$time %||% detect_time_channel(data, verbose)

    ## select and rename channels
    ## names take priority over clashing data names
    channels <- match_channels(channels, names(data), verbose)
    original <- unlist(channels, use.names = FALSE)
    new <- unlist(lapply(channels, names), use.names = FALSE)
    col_idx <- match(original, names(data))
    names(data) <- rename_duplicates(c(new, names(data)))[-seq_along(new)]
    names(data)[col_idx] <- new
    data <- data[c(col_idx, if (keep_all) setdiff(seq_along(data), col_idx))]
    channels <- lapply(channels, names)

    ## pionirs exports a sample index and numeric tag: keep "Iteration"
    ## beside the time column and "Tag" before the event column when
    ## retained as extra columns rather than named channels
    if (identical(nirs_device, "PIONIRS")) {
        extra <- setdiff(names(data), new)
        relocate <- \(.df, .col, ...) {
            tibble::add_column(
                .df[names(.df) != .col], "{.col}" := .df[[.col]], ...
            )
        }
        if ("Iteration" %in% extra) {
            data <- relocate(data, "Iteration", .after = channels$time)
        }
        if ("Tag" %in% extra && !is.null(channels$event)) {
            data <- relocate(data, "Tag", .before = channels$event[1L])
        }
    }

    ## remove empty rows and columns
    ## drop metadata for an empty event column
    data <- remove_empty_rows_cols(data)
    event <- intersect(channels$event, names(data))
    channels["event"] <- list(if (length(event) > 0L) event)

    ## convert decimal "," to "." and column types by role
    data <- convert_type(data, channels, verbose)

    ## time to numeric, header timestamp parsed lazily only when the
    ## time series lacks an absolute date-time
    time_list <- parse_time_channel(
        data[[channels$time]],
        extract_start_timestamp(raw[seq_len(header_row), ]),
        zero_time
    )
    data[[channels$time]] <- time_list$time
    if (add_timestamp && !is.null(time_list$timestamp)) {
        data <- tibble::add_column(
            data,
            timestamp = time_list$timestamp,
            .after = channels$time
        )
    }

    ## Oxysoft exports a sample index: derive a "time" column in seconds
    ## from the export sample rate, placed in front of the sample column
    if (identical(nirs_device, "Artinis")) {
        sample_rate <- oxysoft_sample_rate(raw[seq_len(header_row), ])
        time_new <- make.unique(c(names(data), "time"), sep = "_")[
            ncol(data) + 1L
        ]
        data <- tibble::add_column(
            data,
            "{time_new}" := data[[channels$time]] / sample_rate,
            .before = channels$time
        )
        channels$time <- time_new

        if (verbose) {
            cli_inform(c(
                "!" = "Oxysoft {.arg sample_rate} = {.val {sample_rate}} Hz.",
                "i" = "{.arg time_channel} = {.field {time_new}} added to \\
                the data frame, in {.cls seconds}."
            ))
        }
    }

    ## validate or estimate sample rate; warn for irregular samples
    sample_rate <- validate_sample_rate(
        data, channels$time, sample_rate, verbose
    )
    detect_irregular_samples(data[[channels$time]], channels$time, verbose)

    ## assign metadata to attributes(data)
    metadata <- list(
        nirs_device = nirs_device,
        nirs_channels = channels$nirs,
        time_channel = channels$time,
        event_channel = channels$event,
        sample_rate = sample_rate,
        start_timestamp = time_list$start_timestamp,
        verbose = verbose
    )

    return(create_mnirs_data(data, metadata))
}

#' Metadata names  of class `"mnirs"`, retrieved with `attr()`
#' @keywords internal
mnirs_metadata <- c(
    "nirs_device",
    "nirs_channels",
    "time_channel",
    "event_channel",
    "sample_rate",
    "start_timestamp",
    "interval_times",
    "interval_span"
)


#' Create an *{mnirs}* data frame with metadata
#'
#' Manually add class `"mnirs"` and metadata to an existing data frame.
#'
#' @param data A data frame with existing metadata (accessed with
#'   `attributes(data)`).
#'
#' @param ... Additional arguments with metadata to add to the data frame.
#'   Can be either seperate named arguments or a list of named values.
#'   - nirs_device
#'   - nirs_channels
#'   - time_channel
#'   - event_channel
#'   - sample_rate
#'   - start_timestamp
#'   - interval_times
#'   - interval_span
#'
#' @details
#' Typically will only be called internally, but can be used to inject
#'   *{mnirs}* metadata into any data frame.
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class `"mnirs"`. Metadata are stored
#'   as attributes and can be accessed with `attributes(data)`.
#'
#' @examples
#' data <- data.frame(
#'     A = 1:3,
#'     B = seq(10, 30, 10),
#'     C = seq(11, 33, 11)
#' )
#'
#' attributes(data)
#'
#' ## inject metadata
#' nirs_data <- create_mnirs_data(
#'     data,
#'     nirs_channels = c("B", "C"),
#'     time_channel = "A",
#'     sample_rate = 1
#' )
#'
#' attributes(nirs_data)
#'
#' @export
create_mnirs_data <- function(data, ...) {
    validate_mnirs_data(data, 1L)

    ## tidy eval ========================================================
    ## capture quosures so bare symbols / tidyselect resolve against `data`
    dots <- rlang::enquos(...)
    args <- Map(\(.q, .nm) {
        if (.nm %in% c("nirs_channels", "time_channel", "event_channel")) {
            parse_channel_name(.q, data)
        } else {
            rlang::eval_tidy(.q)
        }
    }, dots, names(dots) %||% rep("", length(dots)))

    ## overwrite existing attributes and add from incoming metadata
    ## incoming metadata from `...` can be either listed or un-listed
    incoming_metadata <- if (length(args) == 1L && is.list(args[[1L]])) {
        args[[1L]]
    } else {
        args
    }

    #! check missing `utils` dependency
    metadata <- utils::modifyList(attributes(data), incoming_metadata)

    ## preserve grouping: `new_tibble()` resets class, so re-add `grouped_df`
    grp <- if (inherits(data, "grouped_df")) "grouped_df"

    nirs_data <- tibble::new_tibble(
        data,
        class = c("mnirs", grp),
        nirs_device = metadata$nirs_device,
        nirs_channels = metadata$nirs_channels,
        time_channel = metadata$time_channel,
        event_channel = metadata$event_channel,
        sample_rate = metadata$sample_rate,
        start_timestamp = metadata$start_timestamp,
        interval_times = metadata$interval_times,
        interval_span = metadata$interval_span,
    )

    tibble::validate_tibble(nirs_data)

    return(nirs_data)
}


#' Get path to *{mnirs}* example files
#'
#' @param file Name of file as character string. If `NULL`, returns a vector
#' of all available file names.
#'
#' @returns
#' A file path character string for selected example files stored in this
#'   package.
#'
#' @examples
#' ## lists all files
#' example_mnirs()
#'
#' ## partial matching will error if matches multiple
#' try(example_mnirs("moxy"))
#'
#' example_mnirs("moxy_ramp")
#'
#' @export
example_mnirs <- function(file = NULL) {
    dir_files <- list.files(
        system.file("extdata", package = "mnirs"),
        pattern = "^[^~]" ## exclude open files
    )

    if (is.null(file)) {
        return(dir_files)
    }

    matches <- grep(file, dir_files, fixed = TRUE, value = TRUE)
    if (length(matches) > 1L) {
        cli_abort(c(
            "x" = "Multiple files match {.val {file}}:",
            "i" = "Matching files: {.val {matches}}"
        ))
    }

    file <- match.arg(file, choices = dir_files)
    return(system.file("extdata", file, package = "mnirs", mustWork = TRUE))
}
