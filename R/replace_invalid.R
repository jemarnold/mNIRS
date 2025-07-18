#' Replace Invalid Values
#'
#' Detect specific values such as `c(0, 100)` in vector data and replaces
#' with `NA` or the local median value.
#'
#' @param x A numeric vector
#' @param values A numeric vector of values to be replaced, e.g.
#'  `values = c(0, 100)`.
#' @param width A numeric scalar for the window length of `(2 · width + 1)` samples.
#' @param return Indicates whether outliers should be replaced with `NA`
#'  (*default*) or the local `"median"` value.
#'
#' @details
#' Useful to overwrite known invalid/nonsense values, such as `0`, `100`, or `102.3`.
#'
#' *TODO: allow for overwriting all values greater or less than known values.*
#'
#' @seealso [pracma::hampel()]
#'
#' @examples
#' set.seed(13)
#' (x <- sample.int(10, 20, replace = TRUE))
#' (y <- replace_invalid(x, values = c(1, 10), width = 5))
#'
#' @return A numeric vector of filtered data.
#'
#' @export
replace_invalid <- function(
        x,
        values,
        width,
        return = c("NA", "median")
) {

    return <- match.arg(return)

    ## validation: `x` must be a numeric vector
    if (!is.numeric(x)) {
        cli::cli_abort("{.arg x} must be a {.cls numeric} vector.")
    }

    ## validation: `values` must be a numeric vector
    if (!is.numeric(values)) {
        cli::cli_abort(paste(
            "{.arg values} must be a {.cls numeric} vector."))
    }

    ## if `return = "median"` then return local median values
    if (return == "median") {

        ## validation: `width` must be a numeric scalar
        if (!is.numeric(width) | length(width) > 1) {
            cli::cli_abort(paste("{.arg width} must be a {.cls numeric} scalar."))
        }

        ## validation: `width` must be shorter than half length(x)
        if (width >= ceiling(length(x)/2)) {
            cli::cli_abort(paste(
                "{.arg width} must be less than half the length of {.arg x}."))
        }

        y <- x
        n <- length(x)
        for (i in 1:n) {
            # Calculate the window bounds, ensuring they stay within vector limits
            start_idx <- max(1, i - width)
            end_idx <- min(n, i + width)
            x0 <- median(x[start_idx:end_idx])
            if (x[i] %in% values) {
                y[i] <- if (return == "median") {x0} else {NA_real_}
            }
        }

    } else {
        ## if `return = "NA"` then simply overwrite to `NA`
        y <- x
        n <- length(x)
        for (i in 1:n) {
            if (x[i] %in% values) {
                y[i] <- NA_real_
            }
        }

    }

    return(y)
}
