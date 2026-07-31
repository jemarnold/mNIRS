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
    invisible(x)
}

