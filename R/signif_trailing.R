#' Format numbers for display as character strings
#'
#' `signif_trailing()` converts numeric values to character strings to
#' preserve trailing zeroes
#'
#' @param x A numeric vector.
#' @param digits An integer specifying the number of *decimal places* or
#'   *significant figures* to preserve. Negative `digits` values will round to
#'   the nearest whole value of `10^(digits)`.
#' @param format Indicates how to treat `digits`. Either the desired number of
#'   decimal places (`format = "digits"`, the *default*) or significant figures
#'   after the decimal place (`format = "signif"`).
#' @param trim Logical; if `TRUE` (the *default*), caps `digits` at the
#'   number of decimal places or significant figures observed in `x`. If
#'   `FALSE`, uses the exact `digits` value.
#'
#' @details
#' `signif_trailing()`
#'
#' - Negative `digits` round to the respective integer place, e.g.
#'   `signif_trailing(123, digits = -1)` returns `"120"`.
#'
#' @returns `signif_trailing()` returns a character vector of formatted numbers
#'   the same length as `x`.
#'
#' @seealso [formatC()], [round()], [signif()]
#'
#' @rdname signif_trailing
#' @order 1
#' @keywords internal
signif_trailing <- function(
    x,
    digits = 2L,
    format = c("digits", "signif"),
    trim = TRUE
) {
    format <- match.arg(format)

    if (format == "digits") {
        validate_numeric(x)
        validate_numeric(digits, 1, c(-Inf, Inf), FALSE, TRUE)
        if (trim) digits <- min(digits, count_decimals(x))
        formatC_x <- round(x, digits)
        formatC_format <- "f"
    } else {
        ## if whole digits >= sig figs, return rounded whole number
        ## x & digits validated by `signif_whole`
        if (trim) digits <- min(digits, count_sigfigs(x))
        formatC_x <- signif_whole(x, digits)
        formatC_format <- "fg"
    }

    result <- formatC(
        x = formatC_x,
        digits = max(0, digits),
        format = formatC_format,
        flag = "#"
    )

    ## remove trailing `.` or "NA"
    return(gsub("\\.$|NA| ", "", result))
}


#' Count maximum decimal places across a numeric vector
#'
#' Returns the largest number of decimal places present in any finite,
#' non-NA element of `x`. Used internally by `signif_trailing()` for
#' `format = "digits"`.
#'
#' @param x A numeric vector.
#' @returns A single non-negative integer.
#'
#' @keywords internal
count_decimals <- function(x) {
    ## handle non-finite elements where nchar shouldn't apply
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(0L)
    }

    ## format with enough precision to capture true decimal places
    txt <- format(x, scientific = FALSE, trim = TRUE)
    decimals <- vapply(strsplit(txt, ".", fixed = TRUE), \(.s) {
        if (length(.s) < 2L) 0L else nchar(.s[[2L]])
    }, integer(1))

    return(max(decimals))
}


#' Count maximum significant figures across a numeric vector
#'
#' Returns the largest number of significant figures present in any finite,
#' non-NA element of `x`. Used internally by `signif_trailing()` for
#' `format = "signif"`.
#'
#' @param x A numeric vector.
#' @returns A single positive integer (minimum 1).
#'
#' @keywords internal
count_sigfigs <- function(x) {
    ## handle non-finite elements where nchar shouldn't apply
    x <- x[is.finite(x) & x != 0]
    if (length(x) == 0L) {
        return(1L)
    }

    ## format with enough precision, strip sign and leading/trailing zeros
    txt <- format(abs(x), scientific = FALSE, trim = TRUE)
    sigfigs <- vapply(txt, \(.s) {
        ## remove decimal point
        .s <- gsub(".", "", .s, fixed = TRUE)
        ## remove leading zeros
        .s <- sub("^0+", "", .s)
        ## remove trailing zeros only for values without a decimal
        ## (integers represented exactly)
        nchar(.s)
    }, integer(1), USE.NAMES = FALSE)

    return(max(sigfigs))
}


#' Round numbers to significant figures or whole value
#'
#' `signif_whole()` rounds numeric values to a specified number of significant
#' figures, or the nearest whole value if the number of digits of `x` are
#' greater than `digits`.
#'
#' @details
#' Decimal rounding is based on the "banker's rounding" default behaviour of
#'   `signif()` and `round()`, where `signif(123.45, 4)` or `round(123.45, 1)`
#'   each return `123.4`.
#'
#' `signif_whole()`
#'
#' - Negative `digits` round to the nearest whole value as if `digits = 0`,
#'   e.g. `signif_whole(123, digits = -5)` still returns `123`.
#'
#' @returns `signif_whole()` returns a numeric vector the same length as `x`.
#'
#' @rdname signif_trailing
#' @order 2
#' @keywords internal
signif_whole <- function(x, digits = 5L) {
    validate_numeric(x)
    validate_numeric(digits, 1, c(-Inf, Inf), FALSE, TRUE)

    ## if whole digits >= sig figs, return rounded whole number
    whole_digits <- floor(log10(abs(x))) + 1
    whole_digits[!is.finite(whole_digits)] <- 1 ## handle 0, NA, Inf
    should_round <- whole_digits >= digits

    ## vectorised ifelse evaluation
    return(ifelse(should_round, round(x), signif(x, digits)))
}


#' Format p-values for display
#'
#' `signif_pvalue()` displays p-values as either formatted numeric strings
#' or significance symbols.
#'
#' @param format Indicates how to treat `digits`. Either the desired 
#'   significance criteria over which to display the absolute p value 
#'   (`format = "digits"`, the *default*), or the smallest significance
#'   criteria to print as less than (`format = "signif"`).
#' @param display Specifies output type, either *"value"* (the *default*) for
#'   formatted numbers or *"symbol"* for significance symbols.
#' @param symbol Character string specifying the significance symbol.
#'   *Default* is "*".
#' @param symbol_repeat Logical indicating whether to repeat symbols for
#'   different significance levels. Default is *FALSE*.
#' @param alpha A numeric value specifying significance threshold.
#'   *Default* is `0.05`.
#'
#' @details
#' `signif_pvalue()`
#'
#' - When `format = "digits"` and e.g. `digits = 3`, `x` is rounded to 3
#'   decimal places, or shown as *"p < 0.001"* below a 3-decimal place 
#'   significance threshold.
#' - `digits = 1` with `format = "digits"` displays *"p < `alpha`"*,
#'   e.g. *"p < 0.05"*.
#' - When `format = "signif"`, `digits` sets the lowest threshold
#'   (e.g. `digits = 3` gives thresholds `alpha`, `0.01`, `0.001`). Values
#'   below `alpha` show the nearest threshold above them, e.g. `p = 0.04`
#'   gives *"p < 0.05"*; `p = 0.009` gives *"p < 0.01"*.
#' - When `display = "symbol"`, if `symbol_repeat = TRUE`: Uses repeated
#'   symbols based on thresholds
#'   `(0.001 = "***", 0.01 = "**", alpha = "*", ns = "")`.
#' - If `symbol_repeat = FALSE`: Shows one symbol `"*"` for p < alpha,
#'   otherwise empty string.
#'
#' @returns `signif_pvalue()` returns a character vector of formatted p-values
#'   or significance symbols the same length as `x`.
#'
#' @rdname signif_trailing
#' @order 3
#' @keywords internal
signif_pvalue <- function(
    x,
    digits = 3L,
    format = c("digits", "threshold"),
    display = c("value", "symbol"),
    symbol = "*",
    symbol_repeat = FALSE,
    alpha = 0.05
) {
    validate_numeric(x)
    validate_numeric(digits, 1, c(0, Inf), FALSE, TRUE)
    format <- match.arg(format)
    display <- match.arg(display)
    validate_numeric(
        alpha, 1, c(0, 1), FALSE, 
        msg1 = "one-element", msg2 = "between {col_blue('[0, 1]')}"
    )

    x[is.infinite(x)] <- NA_real_

    if (display == "symbol" && symbol_repeat) {
        return(strrep(symbol, 3L - findInt_mnirs(x, c(0.001, 0.01, alpha))))
    }
    if (display == "symbol") {
        return(ifelse(x >= alpha, "", symbol))
    }

    if (format == "digits") {
        threshold <- 10^-digits
        if (digits == 1) {
            digits <- -floor(log10(alpha))
            threshold <- alpha
        }
        return(ifelse(
            x < threshold,
            sprintf("p < %.*f", digits, threshold),
            paste0("p = ", signif_trailing(x, digits, format, FALSE))
        ))
    }

    ## format = "threshold": find nearest conventional threshold above x;
    ## digits defines the lowest threshold (e.g. digits = 3 -> 0.001)
    alpha_exp <- -floor(log10(alpha))
    thresholds <- unique(c(alpha, 10^-(alpha_exp:digits)))
    nearest <- vapply(x, \(.p) {
        above <- thresholds[thresholds > .p]
        if (length(above) == 0L) NA_real_ else min(above)
    }, numeric(1))
    
    # compute formatted_threshold only for non-NA nearest values
    formatted_threshold <- rep(NA_character_, length(nearest))
    valid <- !is.na(nearest)
    if (any(valid)) {
        fd <- as.integer(-floor(log10(nearest[valid])))
        formatted_threshold[valid] <- sprintf("%.*f", fd, nearest[valid])
    }

    return(ifelse(
        is.na(x) | x >= alpha,
        paste0("p = ", signif_trailing(x, digits, "digits", FALSE)),
        paste0("p < ", formatted_threshold)
    ))
}


#' Generate numeric sequence from range of a vector
#'
#' Creates a numeric sequence spanning the range of input vector with either a
#' specified step size or a desired output length.
#'
#' @param x A numeric vector.
#' @param by A numeric step size for the output sequence. *Default* is `1`. 
#'   Sign determines order of returned vector (negative `by` returns a
#'   descending sequence). `direction` takes precedence over `by` sign.
#' @param length.out A positive integer giving the desired length of the
#'   sequence. *Default* is `NULL`. If supplied, takes precedence over `by`.
#' @param direction Order of returned vector. Either `"up"` for ascending or
#'   `"down"` for descending. If supplied, takes precedence over the `by` sign.
#'
#' @details
#' The output vector will likely be a different length than the input `x`.
#'
#' @returns A numeric vector spanning the range of the input `x`.
#'
#' @seealso [seq()], [range()]
#'
#' @keywords internal
seq_range <- function(
    x,
    by = 1,
    length.out = NULL,
    direction = c("up", "down")
) {
    if (!is.null(length.out) && !missing(by)) {
        cli_inform(c("i" = "{.arg length.out} overrides {.arg by}."))
    }
    ## explicit `direction` wins; otherwise infer from sign of `by`
    direction <- if (missing(direction)) {
        if (by < 0) "down" else "up"
    } else {
        match.arg(direction)
    }
    x_range <- range(x, na.rm = TRUE)
    if (direction == "down") x_range <- rev(x_range)
    if (!is.null(length.out)) {
        return(seq(x_range[1L], x_range[2L], length.out = length.out))
    }
    by <- if (direction == "down") -abs(by) else abs(by)
    return(seq(x_range[1L], x_range[2L], by = by))
}


#' Wrap vector elements
#'
#' Rotates vector elements by moving the first `n` elements to the end.
#'
#' @param x A vector.
#' @param n An integer specifying number of elements to move from start to end.
#'   *Default* is `0` (no wrapping). Negative values move elements from end
#'   to start.
#'
#' @details
#' The function:
#' - Returns `x` unchanged if `n = 0`.
#' - Moves first `n` elements to the end of the vector.
#' - For negative `n`, effectively moves elements from end to start.
#' - If `n` is larger than `length(x)`, positions wrap around.
#'
#' @returns A vector with all the same elements as `x`.
#'
#' @keywords internal
wrap <- function(x, n = 0L) {
    validate_numeric(n, 1, integer = TRUE)

    if (length(x) == 0L) {
        return(x)
    }

    n <- n %% length(x)
    if (n == 0L) {
        return(x)
    }

    return(c(x[-(1:n)], x[1:n]))
}


#' Detect if numeric values fall within range of a vector
#'
#' Vectorised check for `x %in% vec`, inclusive or exclusive of left and right
#' boundary values, specified independently.
#'
#' @param x A numeric vector.
#' @param vec A numeric vector from which `left` and `right` boundary values
#'   for `x` will be taken.
#' @param inclusive A character vector to specify which of `left` and/or
#'   `right` boundary values should be included in the range, or both (the
#'   default), or excluded if `FALSE`.
#'
#' @details
#' `inclusive = FALSE` can be used to test for positive non-zero values:
#'   `within(x, c(0, Inf), inclusive = FALSE)`.
#'
#' @returns A logical vector the same length as `x`.
#'
#' @seealso [dplyr::between()]
#'
#' @keywords internal
within <- function(x, vec, inclusive = c("left", "right")) {
    if (!is.numeric(x)) {
        abort_validation(substitute(x))
    }
    if (!is.numeric(vec)) {
        abort_validation(substitute(vec))
    }
    inclusive <- match.arg(
        as.character(inclusive), ## force FALSE to character
        choices = c("left", "right", "FALSE"),
        several.ok = TRUE
    )

    ## extract bounds from vec
    left <- min(vec, na.rm = TRUE)
    right <- max(vec, na.rm = TRUE)

    if ("FALSE" %in% inclusive) {
        return(x > left & x < right)
    }

    left_op <- if ("left" %in% inclusive) `>=` else `>`
    right_op <- if ("right" %in% inclusive) `<=` else `<`

    return(left_op(x, left) & right_op(x, right))
}
