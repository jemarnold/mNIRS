#' Validate `{mnirs}` parameters
#'
#' Resolve and validate *{mnirs}* metadata and perform basic data quality
#' checks.
#'
#' @param data A data frame of class *"mnirs"* containing time series data and
#'   metadata.
#'
#' @param nirs_channels A character vector giving the names of mNIRS columns to
#'   operate on. Must match column names in `data` exactly.
#'   - If `NULL` (default), the `nirs_channels` metadata attribute of `data` is
#'     used.
#'
#' @param time_channel A character string naming the time or sample column. 
#'   Must match a column name in `data` exactly.
#'   - If `NULL` (default), the `time_channel` metadata attribute of `data` is
#'     used.
#'
#' @param event_channel A character string naming the event/lap column. Must 
#'   match a column name in `data` exactly.
#'   - If `NULL` (default), the `event_channel` metadata attribute of `data` is
#'     used.
#'
#' @param as_list Logical. Default is `FALSE`. If `nirs_channels` is specified
#'   as a list, it will be coerced to a flat character vector and an
#'   information message is displayed (when `verbose = TRUE`). If `TRUE`,  
#'   `nirs_channels` is returned as-is, i.e. as a list for callers which 
#'   require it. 
#'
#' @param required Logical. Default is `TRUE`. `event_channel` must be
#'   present or detected in metadata. If `FALSE`, `event_channel` may be `NULL`.
#'
#' @param x A numeric vector.
#'
#' @param sample_rate A numeric sample rate in Hz.
#'   - If `NULL` (default), the `sample_rate` metadata attribute of `data` will
#'     be used if detected, or the sample rate will be estimated from
#'     `time_channel`.
#'
#' @param elements An integer. Default is `Inf`. The number of numeric elements
#'   expected in `x`.
#'
#' @param range A two-element numeric vector giving the valid range for `x`.
#'
#' @param inclusive A character vector specifying which boundaries of `range`
#'   are included. Any of `"left"`, `"right"` (default is both). Use `FALSE` to
#'   exclude both endpoints.
#'
#' @param integer Logical. Default is `FALSE`. If `TRUE`, validate `x` as
#'   integer-like values using [rlang::is_integerish()]. Otherwise tested as a
#'   numeric value.
#' 
#' @param allow_na Logical. Default is `FALSE`. If `TRUE`, allows pass through
#'   of `NA` to the returned numeric/integer vector.
#'
#' @param msg1,msg2 A character string appended to the [cli::cli_abort()]
#'   message when numeric validation fails.
#'
#' @param env The calling environment or a defused call, used to report
#'   errors and warnings as coming from the user-facing function rather
#'   than the validator.
#'
#' @inheritParams read_mnirs
#'
#' @details
#' `validate_mnirs()` is an internal documentation topic for a set of
#' validators used throughout the package. These validators:
#'
#' - Prefer explicit user-supplied arguments.
#' - Fall back to *"mnirs"* metadata attributes when available.
#' - Fail fast with informative [cli::cli_abort()] messages when values are
#'   missing or invalid.
#'
#' @returns Returns the validated object (e.g. a resolved `time_channel`
#'   string), or invisibly returns `NULL` for successful validations. On
#'   failure, an error is thrown via [cli::cli_abort()].
#'
#' @name validate_mnirs
#' @keywords internal
NULL

#' validate_numeric abort message construction
#' @keywords internal
abort_validation <- function(
    name,
    integer = FALSE,
    msg1 = "",
    msg2 = "",
    env = rlang::caller_env()
) {
    type <- if (integer) "integer" else "numeric"

    cli_abort(c(
        "x" = paste0(
            "{.arg {name}} must be a valid ", msg1, " {.cls {type}} ", msg2
        )
    ), call = env)
}


#' @rdname validate_mnirs
validate_numeric <- function(
    x,
    elements = Inf,
    range = NULL,
    inclusive = c("left", "right"),
    integer = FALSE,
    allow_na = FALSE,
    msg1 = "",
    msg2 = "",
    env = rlang::caller_env()
) {
    ## pass through NULL
    if (is.null(x)) {
        return(invisible(NULL))
    }

    name <- rlang::as_label(substitute(x))

    ## cheap early type check
    if (!is.numeric(x)) {
        abort_validation(name, integer, msg1, msg2, env)
    }

    ## valid elements length — skip NA scan when allow_na = TRUE
    if (!allow_na) {
        valid <- !is.na(x)
        n_valid <- sum(valid)
        if (n_valid == 0L) abort_validation(name, integer, msg1, msg2, env)
    } else {
        n_valid <- length(x)
    }

    ## elements check
    if (is.finite(elements) && n_valid != elements) {
        abort_validation(name, integer, msg1, msg2, env)
    }

    ## subset once for range/integer checks
    needs_subset <- !is.null(range) || integer
    if (needs_subset && !allow_na) {
        x_valid <- if (n_valid < length(x)) x[valid] else x
    }

    ## range check
    if (!is.null(range) && !all(within(x_valid, range, inclusive))) {
        abort_validation(name, integer, msg1, msg2, env)
    }
    ## expensive integer check
    if (integer && !rlang::is_integerish(x_valid)) {
        abort_validation(name, integer, msg1, msg2, env)
    }

    return(invisible())
}


#' @rdname validate_mnirs
validate_mnirs_data <- function(
    data, 
    ncol = 2L, 
    env = rlang::caller_env()
) {
    ## validate is a data frame with at least two columns
    if (!is.data.frame(data) || length(data) < ncol) {
        cli_abort(c(
            "x" = "{.arg data} must be a data frame with at least \\
            {.val {ncol}} column{?s}."
        ), call = env)
    }

    return(invisible())
}

#' Parse channel expressions for NSE
#'
#' Converts quosures to character vectors, handling bare symbols, character
#' strings, lists, and tidyselect expressions.
#'
#' @param channel A quosure from `rlang::enquo()`.
#' @param data A data frame for tidyselect context.
#' @param env Environment for symbol evaluation (typically the quosure
#'   environment).
#'
#' @returns A character vector, list of character vectors, or `NULL`.
#'
#' @keywords internal
parse_channel_name <- function(
    channel,
    data,
    env = rlang::caller_env()
) {
    if (rlang::quo_is_null(channel)) {
        return(NULL)
    }

    channel_raw <- rlang::quo_get_expr(channel)

    ## already-evaluated list or character
    if (is.list(channel_raw) || is.character(channel_raw)) {
        return(channel_raw)
    }

    ## bare symbol: check if column name, otherwise evaluate
    if (rlang::quo_is_symbol(channel)) {
        sym_name <- rlang::as_name(channel)
        if (sym_name %in% names(data)) {
            return(sym_name)
        }

        ## external object: evaluate and return directly
        result <- tryCatch(
            rlang::eval_tidy(channel, env = env),
            error = \(e) NULL
        )
        if (is.list(result) || is.character(result)) {
            return(result)
        }
        return(sym_name)
    }

    ## list() call: recurse on each element
    if (rlang::is_call(channel_raw, "list")) {
        result <- lapply(rlang::call_args(channel_raw), \(.arg) {
            parse_channel_name(rlang::new_quosure(.arg, env = env), data, env)
        })
        return(unname(result))
    }

    ## evaluate: tidyselect first, then fallback to direct evaluation
    ## handles c(), tidyselect helpers, symbols, and external objects
    tryCatch(
        unname(names(tidyselect::eval_select(channel, data))),
        error = \(e) {
            result <- rlang::eval_tidy(channel, env = env)
            if (is.list(result) || is.character(result)) result else NULL
        }
    )
}


#' @rdname validate_mnirs
validate_nirs_channels <- function(
    nirs_channels,
    data,
    verbose = FALSE, ## hint only emitted on the as_list grouping path
    as_list = FALSE,
    env = rlang::caller_env()
) {
    ## parse tidy eval input
    if (rlang::is_quosure(nirs_channels)) {
        nirs_channels <- parse_channel_name(
            nirs_channels, data, rlang::quo_get_env(nirs_channels)
        )
    }
    nirs_unlisted <- unlist(nirs_channels)

    ## if not defined, check metadata
    if (is.null(nirs_unlisted) || length(nirs_unlisted) == 0) {
        nirs_channels <- attr(data, "nirs_channels") ## should be vector
        nirs_unlisted <- nirs_channels
        if (verbose && as_list && !is.null(nirs_unlisted)) {
            cli_inform(c(
                "i" = "{.arg nirs_channels} = \\
                {col_blue({deparse(nirs_unlisted)})} \\
                grouped together from metadata."
            ), call = env)
        }
    }

    ## if still not defined, return error
    if (is.null(nirs_unlisted)) {
        cli_abort(c(
            "x" = "{.arg nirs_channels} not detected in metadata.",
            "i" = "Check your data attributes or define \\
            {.arg nirs_channels} explicitly."
        ), call = env)
    }

    ## validate exists in data
    if (!is.character(nirs_unlisted) || !all(nirs_unlisted %in% names(data))) {
        cli_abort(c(
            "x" = "{.arg nirs_channels} not detected in {.arg data}.",
            "i" = "Channel names are case-sensitive and must match exactly."
        ), call = env)
    }

    ## validate is numeric and has >=2 valid values
    invalid_channels <- vapply(data[nirs_unlisted], \(.x) {
        !is.numeric(.x) || sum(is.finite(.x)) < 2
    }, logical(1))

    if (sum(invalid_channels) > 0) {
        cli_abort(c(
            "x" = "{.arg nirs_channels} must contain valid {.cls numeric} data."
        ), call = env)
    }

    ## preserve list grouping for callers that need it (as_list = TRUE)
    if (as_list) {
        return(nirs_channels)
    }

    ## default: return a flat character vector of channel names
    return(nirs_unlisted)
}


#' @rdname validate_mnirs
validate_time_channel <- function(
    time_channel,
    data,
    env = rlang::caller_env()
) {
    ## parse tidy eval input
    if (rlang::is_quosure(time_channel)) {
        time_channel <- parse_channel_name(
            time_channel, data, rlang::quo_get_env(time_channel)
        )
    }

    ## if not defined, check metadata
    if (is.null(time_channel)) {
        time_channel <- attr(data, "time_channel")
    }

    ## if still not defined, return error
    if (is.null(time_channel)) {
        cli_abort(c(
            "x" = "{.arg time_channel} not detected in metadata.",
            "i" = "Check your data attributes or define \\
            {.arg time_channel} explicitly."
        ), call = env)
    }

    ## validate exists in data
    if (!is.character(time_channel) || !time_channel %in% names(data)) {
        cli_abort(c(
            "x" = "{.arg time_channel} not detected in {.arg data}.",
            "i" = "Channel names are case-sensitive and must match exactly."
        ), call = env)
    }

    ## validate is numeric and has >=2 valid values
    if (
        !is.numeric(data[[time_channel]]) ||
            sum(is.finite(data[[time_channel]])) < 2
    ) {
        cli_abort(c(
            "x" = "{.arg time_channel} must contain valid {.cls numeric} data."
        ), call = env)
    }

    return(time_channel)
}


#' @rdname validate_mnirs
validate_event_channel <- function(
    event_channel,
    data,
    required = TRUE,
    env = rlang::caller_env()
) {
    ## parse tidy eval input
    if (rlang::is_quosure(event_channel)) {
        event_channel <- parse_channel_name(
            event_channel, data, rlang::quo_get_env(event_channel)
        )
    }
    ## if not defined, check metadata
    if (is.null(event_channel)) {
        event_channel <- attr(data, "event_channel")
    }

    ## if still not defined, return error
    if (is.null(event_channel) && required) {
        cli_abort(c(
            "x" = "{.arg event_channel} not detected in metadata.",
            "i" = "Check your data attributes or define {.arg event_channel} \\
            explicitly."
        ), call = env)
    } else if (is.null(event_channel) && !required) {
        ## return event_channel = NULL if not required
        return(event_channel)
    }

    ## validate exists in data
    if (!is.character(event_channel) || !event_channel %in% names(data)) {
        cli_abort(c(
            "x" = "{.arg event_channel} not detected in {.arg data}.",
            "i" = "Channel names are case-sensitive and must match exactly."
        ), call = env)
    }

    ## validate column type: must be character or integerish
    col <- data[[event_channel]]
    if (!is.character(col) && !rlang::is_integerish(col)) {
        cli_abort(c(
            "x" = "{.arg event_channel} must contain valid {.cls character} \\
            event labels or {.cls integer} lap numbers."
        ), call = env)
    }

    ## check for empty column — character columns also check for empty strings
    valid_values <- if (is.character(col)) {
        !is.na(col) & nzchar(col)
    } else {
        !is.na(col)
    }
    if (sum(valid_values) == 0) {
        cli_abort(c(
            "x" = "{.arg event_channel} must contain valid {.cls character} \\
            event labels or {.cls integer} lap numbers."
        ), call = env)
    }

    return(event_channel)
}


#' @rdname validate_mnirs
estimate_sample_rate <- function(x, env = rlang::caller_env()) {
    ## estimate samples per second
    sample_rate_raw <- 1 / median(diff(x), na.rm = TRUE)

    if (!is.finite(sample_rate_raw) || sample_rate_raw == 0) {
        cli_abort(c(
            "x" = "Unable to estimate {.arg sample_rate}.",
            "i" = "Check that {.arg time_channel} values are consistent.",
            "i" = "Set {.arg sample_rate} = {.cls numeric}."
        ), call = env)
    }

    pretty_vals <- c(
        0.25, 0.5, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 50, 60, 75, 100
    )
    return(pretty_vals[which.min(abs(pretty_vals - sample_rate_raw))])
}


#' @rdname validate_mnirs
validate_sample_rate <- function(
    data,
    time_channel,
    sample_rate,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## if not defined, check metadata
    sample_rate <- sample_rate %||% attr(data, "sample_rate")

    ## skip estimation when sample_rate ~ 1 (integer time_channel, the
    ## common cheap case); estimate when NULL or verbose check needed
    near_one <- \(x) isTRUE(all.equal(1, x, tol = 1e-3, scale = 1))
    if (is.null(sample_rate) || (verbose && !near_one(sample_rate))) {
        ## time_channel must be validated before this
        t_vec <- as.numeric(data[[time_channel]])
        ## will error if unable to estimate sample_rate
        sample_rate_est <- estimate_sample_rate(t_vec, env)
    }

    ## if still not defined, use estimated sample_rate
    if (is.null(sample_rate)) {
        sample_rate <- sample_rate_est
        if (verbose) {
            cli_inform(c(
                "!" = "Estimated {.arg sample_rate} = {.val {sample_rate}} Hz.",
                "i" = "Define {.arg sample_rate} explicitly to override."
            ), call = env)
        }
    }

    ## validate has one numeric value
    validate_numeric(
        sample_rate, 1, c(0, Inf), FALSE, 
        msg1 = "one-element positive", env = env
    )

    ## warn when user-provided sample_rate disagrees with the estimate
    if (
        verbose && !near_one(sample_rate) &&
        !isTRUE(all.equal(sample_rate_est, sample_rate, tol = 0.5, scale = 1))
    ) {
        cli_warn(c(
            "!" = "`sample_rate = {.val {sample_rate}}` appears to be \\
            inconsistent with {.arg time_channel}. Estimated \\
            `sample_rate = {.val {sample_rate_est}}`.",
            "i" = "Check that your sample rate and {.arg time_channel} \\
            values are consistent."
        ), call = env)
    }

    return(sample_rate)
}

#' @rdname validate_mnirs
validate_width_span <- function(
    width = NULL,
    span = NULL,
    verbose = TRUE,
    msg = "",
    env = rlang::caller_env()
) {
    if (is.null(c(width, span))) {
        cli_abort(c(
            "x" = "Window size undefined",
            "i" = paste(
                "One of {.arg width} or {.arg span} must be defined", msg
            )
        ), call = env)
    }
    validate_numeric(
        width, 1, c(1, Inf), integer = TRUE, 
        msg1 = "one-element positive", env = env
    )
    validate_numeric(
        span, 1, c(0, Inf), msg1 = "one-element positive", env = env
    )
    if (verbose && !is.null(width) && !is.null(span)) {
        cli_inform(c(
            "i" = "{.arg width} = {.val {width}} overrides {.arg span}."
        ), call = env)
    }
}


#' @rdname validate_mnirs
validate_x_t <- function(
    x,
    t,
    allow_na = FALSE,
    env = rlang::caller_env()
) {
    ## exclude NULL by defaulting to allow_na character
    x <- x %||% character()
    t <- t %||% character()
    validate_numeric(x, allow_na = allow_na, env = env)
    validate_numeric(t, allow_na = allow_na, env = env)
    if (length(x) != length(t)) {
        cli_abort(c(
            "x" = "{.arg x} and {.arg t} must be {.cls numeric} vectors \\
            of equal length."
        ), call = env)
    }
}


#' Validate t0
#' @keywords internal
validate_t0 <- function(
    t0,
    data,
    t_vec,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## fall back to metadata or zero
    ## TODO fall back to first non-zero time value rather than 0?
    t0 <- t0 %||% attr(data, "interval_times") %||% 0
    validate_numeric(t0, 1L, env = env)

    if (length(which(t_vec <= t0)) == 0L) {
        if (verbose) {
            cli_warn(c(
                "!" = "No observations where {.arg time_channel} <= \\
                `t0 = {.val {t0}}`.",
                "i" = "All samples included in response."
            ), call = env)
        }
        t0 <- t_vec[1L]
    }
    if (t0 > t_vec[length(t_vec)]) {
        cli_abort(c(
            "x" = "No observations in {.arg time_channel} before {.arg t0}.",
            "i" = "{.arg t0} must be specified within the range of \\
            {.arg time_channel}."
        ), call = env)
    }

    return(t0)
}


#' wrap findInterval: informative 'time_channel' error message
#' @keywords internal
findInt_mnirs <- function(..., env = rlang::caller_env()) {
    withCallingHandlers(findInterval(...), error = \(e) {
        cli_abort(c(
            "x" = "Duplicate or irregular {.arg time_channel} samples \\
            detected.",
            "i" = "{.arg time_channel} must be sorted without repeating, \\
            decreasing, or {.val {NA}}s."
        ), call = env)
    })
}