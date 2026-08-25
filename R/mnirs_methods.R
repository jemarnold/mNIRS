#' Methods for mnirs objects
#'
#' Generic methods for objects of class `"mnirs"`.
#'
#' @param x Object of class `"mnirs"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @returns
#' \item{`print`}{Returns `x` without class attributes.}
#'
#' @examples
#' x <- read_mnirs(
#'     example_mnirs("train.red"),
#'     nirs_channels = c(smo2 = "SmO2"),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     verbose = FALSE
#' ) |>
#'     resample_mnirs(method = "linear", verbose = FALSE) |>
#'     extract_intervals(
#'         start = by_time(2452, 3168),
#'         span = c(-60, 120),
#'         verbose = FALSE
#'     )
#'
#' print(x)
#'
#' @export
print.mnirs <- function(x, ...) {
    class(x) <- setdiff(class(x), "mnirs")
    print(x, ...)
    return(invisible(x))
}


#' Methods for mnirs_kinetics objects
#'
#' Generic methods for objects returned from [analyse_kinetics()].
#'
#' @param x Object of class `"mnirs_kinetics"` returned from
#'   [analyse_kinetics()].
#' @param ... Additional arguments.
#'
#' @returns
#' \item{`print`}{Returns a model summary}
#'
#' @examples
#' result <- read_mnirs(
#'     example_mnirs("train.red"),
#'     nirs_channels = c(smo2 = "SmO2"),
#'     time_channel = c(time = "Timestamp (seconds passed)"),
#'     zero_time = TRUE,
#'     verbose = FALSE
#' ) |>
#'     resample_mnirs(method = "linear", verbose = FALSE) |>
#'     extract_intervals(
#'         group_intervals = "distinct", ## return each interval distinctly
#'         start = by_time(368, 1084),
#'         span = c(-20, 90),
#'         zero_time = TRUE,
#'         verbose = FALSE
#'     ) |>
#'     analyse_kinetics(
#'         method = "peak_slope",
#'         span = 10,
#'         verbose = FALSE
#'     )
#'
#' print(result)
#'
#' @export
print.mnirs_kinetics <- function(x, ...) {
    coefs <- x$coefficients
    ## remove columns unnecessary for display
    drop_start_time <- isTRUE(all(coefs$start_time == 0))
    drop_cols <- c(
        if (drop_start_time) "start_time",
        "tau_mult",
        grep("fitted$", names(coefs), value = TRUE)
    )
    coefs <- coefs[, !names(coefs) %in% drop_cols, drop = FALSE]
    nrows <- nrow(coefs)

    ## prep numeric values for display; ignore all na cols
    numeric_cols <- vapply(coefs, \(.x) {
        is.numeric(.x) && !all(is.na(.x))
    }, logical(1))
    ## prep only non-NA values
    coefs[, numeric_cols] <- lapply(coefs[, numeric_cols, drop = FALSE], \(.x) {
        .x[!is.na(.x)] <- signif_trailing(.x[!is.na(.x)], 4L, "signif")
        .x[is.na(.x)] <- "NA"
        .x
    })
    ## display titles
    method_labels <- c(
        response_time = "Fractional Response Time",
        peak_slope = "Peak Linear Response Rate",
        monoexponential = "Monoexponential One-Phase Kinetics",
        exponential_drift = "Exponential-Linear Two-Phase Kinetics",
        biexponential = "Biexponential Two-Phase Kinetics",
        sigmoidal = "Sigmoidal Inflection Kinetics"
    )

    cat("\n")
    cat(method_labels[[x$method]])
    cat("\n")
    cat("    Model Coefficients:")
    cat("\n")
    if (nrows <= 10) {
        print(coefs)
    } else {
        start_ids <- 1:5
        end_idx <- (nrows - 4):nrows
        ## get first and last rows
        display <- rbind(coefs[start_ids, ], coefs[end_idx, ])

        ## format data frame with row numbers
        output <- utils::capture.output(print(display, row.names = FALSE))

        ## insert row numbers and spacer
        header <- output[1L]
        lines <- output[-1L]

        ## add row numbers to data lines
        top_lines <- paste0(sprintf("%6d:", start_ids), " ", lines[1:5])
        bottom_lines <- paste0(sprintf("%6d:", end_idx), " ", lines[6:10])

        ## create spacer with column count
        spacer <- sprintf("    --- %d rows omitted", nrows - 10)

        ## print output
        cat(sprintf("        %s", header), "\n", sep = "")
        cat(top_lines, sep = "\n")
        cat(spacer, "\n")
        cat(bottom_lines, sep = "\n")

        ## TODO add instructions how to retrieve returned object?
    }
    cat("\n\n")

    return(invisible(x))
}
