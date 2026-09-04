#' Read *{mnirs}* data from file
#'
#' Import time-series data exported from common muscle NIRS (mNIRS) devices and
#' return a [tibble][tibble::tibble-package] (data frame) of class `"mnirs"`
#' with the specified signal channels and metadata.
#'
#' @param file_path Path of the data file to import. Supported file extensions
#'   include `".xls(x)"`, `".csv"`, `".txt"`, and `".ftn(2)"`.
#'
#' @param nirs_channels A character vector of one or more column names
#'   containing mNIRS signals to import. Names must match the file contents
#'   exactly.
#'
#'   - If `NULL` (*default*), `read_mnirs()` attempts to automatically detect
#'     known `nirs_channel` names from the file contents.
#'   - A named character vector is used to rename columns, in the form
#'     `c(renamed = "original_name")`.
#'
#' @param time_channel A single character vector for the time (or sample)
#'   column name. Must match the file contents exactly.
#'
#'   - If `NULL` (*default*), `read_mnirs()` attempts to automatically detect
#'     a time-like column from known device defaults, and/or time-formatted
#'     values.
#'   - A named character vector is used to rename the column, e.g.
#'     `c(time = "original_name")`.
#'
#' @param event_channel An *optional* single character vector for the event
#'   or lap column name. Must match the file contents exactly. A named
#'   character vector is used to rename the column, e.g.
#'   `c(event = "original_name")`.
#'
#' @param sample_rate An *optional* numeric sample rate in Hz. If `NULL`
#'   (*default*), the sample rate is estimated from `time_channel` (see
#'   *Details*).
#'
#' @param add_timestamp Logical. Default is `FALSE`. If `TRUE` and the
#'   source data contain date-time (POSIXct) values, will add a `"timestamp"`
#'   column in addition to the specified `time_channel` as a numeric time
#'   column.
#'
#' @param zero_time Logical. Default is `FALSE`. If `TRUE`, re-calculates
#'   numeric `time_channel` values to start from zero.
#'
#' @param keep_all Logical. `FALSE` (*default*) will only keep the channels
#'   explicitly specified in `channels`. If `TRUE`, will keep all columns
#'   found in the file data table.
#'
#'   - If no `channels` are specified and the NIRS device file format is
#'     recognised, then all columns in the file data table will be returned
#'     to allow exploration of the file.
#'
#' @param verbose Logical. `TRUE` (*default*) will display, and `FALSE` will
#'   silence warnings and information messages helpful for troubleshooting.
#'   Global default can be set via `options(mnirs.verbose = FALSE)`.
#'
#' @details
#' ## Header detection
#' `read_mnirs()` searches the file for a header row containing the requested
#' channel names. The header row does not need to be the first row in the file.
#'
#' - If duplicate column names exist, they are made unique with a numbered
#'   suffix (e.g. `*_1`), and can be renamed accordingly:
#'   `nirs_channels = c(smo2_left = "smo2", smo2_right = "smo2_1")`.
#' - Unnamed columns containing data in the source file will be renamed to
#'   `col_n`, where `n` is the ordered column number in the file (e.g.
#'   `col_6`). *Artinis Oxysoft* files are an exception to this renaming
#'   convention. See *Artinis Oxysoft exports* below).
#'
#' ## Renaming channels
#' All `channels` can be renamed with a named character vector in the form
#' `c(renamed = "original_name")`. The `"original_name"` must match
#' the file contents header row exactly.
#'
#' ## Artinis Oxysoft exports
#' *Artinis Oxysoft* files have numbered data columns, with a "Legend"
#' metadata block with channel names. `read_mnirs()` can detect and rename
#' these channels automatically:
#'
#' - `nirs_channels` names become clean lower-case column names with
#'   underscores (e.g. `"Rx1 - Tx1 O2Hb"` becomes `rx1_tx1_o2hb`). Channels
#'   should still be explicitly renamed pointing to their number, e.g.
#'   `nirs_channels = c(o2hb = 2, hhb = 3)`.
#' - `"(Sample number)"` column is renamed `sample`, and a `time` column
#'   in seconds is automatically derived from the export sample rate.
#' - `"(Event)"` column is renamed `event` and set as `event_channel`.
#' - *Oxysoft* exports a trailing un-numbered column containing optional
#'   event label text. This is renamed `labels` and returned with
#'   `keep_all = TRUE`, or dropped when empty. It can be selected as the
#'   event column explicitly with `event_channel = c(event = "labels")`.
#'
#' Explicit `nirs_channels`, `time_channel`, and `event_channel` renaming
#' (as above) overrides automatically detected names.
#'
#' ## PIONIRS exports
#' *PIONIRS* `.ftn(2)` files are detected with `"Time"` as `time_channel`,
#' `"TagLabel"` as `event_channel`, and `StO2` channels as `nirs_channels`.
#' The `"Iteration"` sample index and numeric `"Tag"` companion columns are
#' returned beside `time_channel` and `event_channel` with `keep_all = TRUE`.
#'
#' ## Time parsing
#' If `time_channel` is left as `NULL`, it can be resolved from a known
#' NIRS device default, or by detecting a time-like column name (e.g.
#' `"time"`, `"hh:mm:ss"`), or by detecting a column with date-time formatted
#' (POSIXct-like) values.
#'
#' If `time_channel` is a date-time (POSIXct) format, it will be converted
#' to numeric and re-based to start from `0`, regardless of `zero_time`.
#'
#' ## Sample rate
#' If `sample_rate` is not specified, it is estimated from differences in
#' `time_channel`. When irregular time sampling is detected, the estimated
#' median `sample_rate` will be approximated as common known recording rate
#' (e.g. an estimated rate of `11` may be rounded to `10 Hz`).
#'
#' If `time_channel` is specified as a sample index (e.g.
#' *Artinis Oxysoft "sample"* or *PIONIRS "Iterations"*), `sample_rate` will
#' be mis-estimated as `1 Hz`. `sample_rate` should be specified explicitly in
#' this case.
#'
#' ## Data cleaning
#' Entirely empty rows and columns are removed. Invalid values (e.g.
#' `c(NaN, Inf, "-")`) are standardized to `NA`. A warning is displayed
#' (respecting `verbose`) when irregular sampling is detected (e.g.
#' non-monotonic, repeated, or unequal `time_channel` values). In this case,
#' it is recommended to use `resample_mnirs()` to standardise the time grid to
#' the desired `sample_rate`.
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
    ## as named `c(new = "original")` by role; role order = column order
    channels <- resolve_channels(raw, device, user, keep_all, verbose)
    header_row <- find_header_row(raw, channels$nirs, device$header_row)
    header <- raw[seq_len(header_row), ]

    ## name data table by header row
    ## blank/duplicate names made unique
    data <- setNames(
        raw[-seq_len(header_row), ],
        rename_duplicates(as.character(raw[header_row, ]))
    )
    channels$time <- channels$time %||%
        name_channels(detect_time_channel(data, verbose))

    ## select, rename, and order channel columns
    selected <- select_channels(data, channels, keep_all, verbose)
    data <- selected$data
    channels <- selected$channels

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
        extract_start_timestamp(header),
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
        sample_rate <- oxysoft_sample_rate(header)
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
