#' Computes rolling local values
#'
#' `compute_window_bounds()`: Compute the start and end indices of rolling
#' windows along a time variable `t`.
#'
#' @param idx A numeric vector of indices of `t` at which to calculate local
#'   windows. All indices of `t` by *default*, or can be used to only calculate
#'   for known indices, such as invalid values of `x`.
#' @param width An integer defining the local window in number of samples
#'   around `idx` in which to perform the operation, according to `align`.
#' @param span A numeric value defining the local window time span around `idx`
#'   in which to perform the operation, according to `align`. In units of
#'   `time_channel` or `t`.
#' @param align Window alignment as *"centre"/"center"* (the *default*),
#'   *"left"*, or *"right"*. Where *"left"* is *forward looking*, and *"right"*
#'   is *backward looking* from the current sample.
#' @inheritParams replace_invalid
#'
#' @returns
#' `compute_window_bounds()`: A `list()` with `start` and `end` integer vectors
#'   the same length as `idx`, giving the inclusive window bounds at each index.
#'
#' @details
#' The local rolling window can be specified by either `width` as the number of
#'   samples, or `span` as the time span in units of `t`. Specifying `width`
#'   is often faster than `span`.
#'
#' `align` defaults to *"centre"* the local window around `idx` between
#'   `[idx - floor((width-1)/2),` `idx + floor(width/2)]` when `width` is
#'   specified. Even `width` values will bias `align` to *"left"*, with the
#'   unequal sample forward of `idx`, effectively returning `NA` at the last
#'   sample index. When `span` is specified, the local window is between
#'   `[t - span/2, t + span/2]`.
#'
#' @rdname compute_helpers
#' @inheritParams validate_mnirs
#' @keywords internal
compute_window_bounds <- function(
    t,
    idx = seq_along(t),
    width = NULL,
    span = NULL,
    align = c("centre", "left", "right"),
    env = rlang::caller_env()
) {
    align <- sub("^center$", "centre", align)
    align <- match.arg(align)
    n <- length(t)

    if (!is.null(width)) {
        ## right = backward looking; left = forward looking
        ## centre = left-biased
        offsets <- switch(
            align,
            centre = c(-floor((width - 1L) / 2L), floor(width / 2L)),
            left = c(0L, width - 1L),
            right = c(-(width - 1L), 0L)
        )
        start_idx <- pmax.int(1L, idx + offsets[1L])
        end_idx <- pmin.int(n, idx + offsets[2L])
    } else {
        offsets <- switch(
            align,
            centre = c(-0.5, 0.5),
            left = c(0, 1),
            right = c(-1, 0)
        ) *
            span
        start_idx <- validate_findInt(
            t[idx] + offsets[1L],
            t,
            left.open = TRUE,
            env = env
        ) +
            1L
        end_idx <- validate_findInt(t[idx] + offsets[2L], t, env = env)
    }

    ## inclusive of x[i] for detect outliers
    return(list(start = start_idx, end = end_idx))
}


#' @description
#' `window_sums()`: Windowed sums by cumulative-sum differencing.
#'
#' @param v A numeric vector to sum within windows. Callers should centre
#'   `v` first to contain floating-point cancellation error.
#' @param bounds A `list()` of `start` and `end` window index vectors from
#'   `compute_window_bounds()`.
#'
#' @returns
#' `window_sums()`: A numeric vector the same length as `bounds$start`.
#'
#' @rdname compute_helpers
#' @keywords internal
window_sums <- function(v, bounds) {
    cs <- cumsum(c(0, v))
    return(cs[bounds$end + 1L] - cs[bounds$start])
}


#' @description
#' `window_min_obs()`: Minimum number of samples spanned by a complete window.
#'
#' @param min_n A lower bound on the returned number of samples.
#'
#' @returns
#' `window_min_obs()`: An integer value.
#'
#' @details
#' `window_min_obs()` converts `span` to a sample count via the estimated
#'   sample rate, less two samples to buffer irregular `t` at the start and
#'   end of each window.
#'
#' @rdname compute_helpers
#' @keywords internal
window_min_obs <- function(
    width,
    span,
    t,
    min_n = 1L,
    env = rlang::caller_env()
) {
    ## `span` only converted when `width` is undefined
    return(
        max(width %||% (floor(span * estimate_sample_rate(t, env)) - 2L), min_n)
    )
}


#' @description
#' `compute_local_mean()`: Compute rolling means from window bounds.
#'
#' @param min_obs The minimum number of samples a window must span to return a
#'   value. Shorter (partial) windows return `NA`.
#' @inheritParams replace_invalid
#'
#' @returns
#' `compute_local_mean()`: A numeric vector the same length as `bounds$start`.
#'
#' @details
#' `compute_local_mean()` computes all window means in O(n) via
#'   [window_sums()]. Values are centred first so the cumulative-sum
#'   differencing error stays around `eps * sqrt(n) * sd`, far below
#'   measurement resolution.
#'
#' @rdname compute_helpers
#' @keywords internal
compute_local_mean <- function(x, bounds, na.rm = FALSE, min_obs = 1L) {
    ## non-finite treated as missing: Inf would otherwise return NaN
    valid <- is.finite(x)
    if (!any(valid)) {
        return(rep(NA_real_, length(bounds$start)))
    }

    ## centre so cumsum differencing error stays ~eps * sqrt(n) * sd
    offset <- mean(x[valid])
    x0 <- x - offset
    x0[!valid] <- 0

    n_window <- bounds$end - bounds$start + 1L
    n_valid <- window_sums(valid, bounds)
    y <- window_sums(x0, bounds) / n_valid + offset

    ## empty (all-NA) windows and partial windows below min_obs return NA
    y[n_valid == 0L | n_window < min_obs] <- NA_real_
    ## NA propagates unless na.rm
    if (!na.rm) {
        y[n_valid < n_window] <- NA_real_
    }
    return(y)
}


#' @description
#' `compute_local_fun()`: Compute a rolling function along `x` from a list of
#' rolling sample windows.
#'
#' @param window_idx A list the same or shorter length as `x` with numeric
#'   vectors for the sample indices of local rolling windows.
#' @param fn A function to pass through for local rolling calculation.
#' @param ... Additional arguments.
#'
#' @returns
#' `compute_local_fun()`: A numeric vector the same length as `x`.
#'
#' @rdname compute_helpers
#' @keywords internal
compute_local_fun <- function(x, window_idx, fn, ...) {
    vapply(seq_along(window_idx), \(.i) {
        fn(x[window_idx[[.i]]], ...)
    }, numeric(1))
}


#' @description
#' `median_no_na()`: Fast median for numeric vectors. Strips `NA`s and
#' replicates `median.default` arithmetic without S3 dispatch.
#'
#' @returns
#' `median_no_na()`: A numeric value.
#'
#' @rdname compute_helpers
#' @keywords internal
median_no_na <- function(w) {
    if (anyNA(w)) {
        w <- w[!is.na(w)]
    }
    n <- length(w)
    if (n == 0L) {
        return(NA_real_)
    }
    half <- (n + 1L) %/% 2L
    if (n %% 2L == 1L) {
        sort.int(w, partial = half)[half]
    } else {
        mean(sort.int(w, partial = half + 0L:1L)[half + 0L:1L])
    }
}


#' @description
#' `compute_col_medians()`: Column medians of an `NA`-padded numeric matrix
#' via a single radix sort. `NA`s sort last per column; medians indexed
#' from per-column valid counts. Matches `median(w, na.rm = TRUE)`.
#'
#' @param m A numeric matrix with one column per rolling window, padded
#'   with `NA` where windows extend beyond the data.
#'
#' @returns
#' `compute_col_medians()`: A numeric vector of length `ncol(m)`.
#'
#' @rdname compute_helpers
#' @keywords internal
compute_col_medians <- function(m) {
    w <- nrow(m)
    nv <- w - colSums(is.na(m))
    ## sort all windows at once: by column, then value, NAs last
    o <- order(col(m), m, na.last = TRUE, method = "radix")
    ms <- matrix(m[o], nrow = w)
    half <- pmax((nv + 1L) %/% 2L, 1L) ## guard nv == 0
    cols <- seq_len(ncol(m))
    lo <- ms[cbind(half, cols)]
    hi <- ms[cbind(half + (nv %% 2L == 0L), cols)] ## even nv: mean of pair
    med <- (lo + hi) / 2
    med[nv == 0L] <- NA_real_
    return(med)
}


#' @description
#' `compute_outliers()`: Computes a vector of local medians and logicals
#' indicating outliers of `x` within rolling windows defined by `width`
#' or `span`.
#'
#' @returns
#' `compute_outliers()`: A `list()` with vectors the same length as `x` for
#' with numeric local medians and logical identifying where `is_outlier`.
#'
#' @rdname compute_helpers
#' @keywords internal
compute_outliers <- function(
    x,
    t,
    outlier_cutoff,
    width = NULL,
    span = NULL,
    env = rlang::caller_env()
) {
    n <- length(x)
    L <- 1.4826 ## 1 / qnorm(0.75): MAD at the 75% percentile of |Z|
    # MAD = median(|x - median(x)|) within each window
    ## median of absolute local residuals from the local median
    if (!is.null(width)) {
        ## width: fixed-size windows, vectorised over an NA-padded matrix.
        ## padding out-of-range cells with NA is equivalent to the clamped
        ## partial edge windows because medians ignore NA.
        ## same offsets as compute_window_bounds(align = "centre")
        offsets <- (-floor((width - 1L) / 2L)):(floor(width / 2L))
        idx <- outer(offsets, seq_len(n), `+`)
        oob <- idx < 1L | idx > n
        idx[oob] <- 1L
        m <- matrix(x[idx], nrow = width)
        m[oob] <- NA_real_

        local_meds <- compute_col_medians(m)
        local_mad <- compute_col_medians(abs(m - rep(local_meds, each = width)))
    } else {
        ## span: variable-size windows, per-window loop over bounds
        bounds <- compute_window_bounds(t, span = span, env = env)
        local_stats <- vapply(seq_len(n), \(.i) {
            w <- x[bounds$start[.i]:bounds$end[.i]]
            local_median <- median_no_na(w)
            local_mad <- median_no_na(abs(w - local_median))

            c(local_median, local_mad)
        }, numeric(2))

        local_meds <- local_stats[1L, ]
        local_mad <- local_stats[2L, ]
    }

    ## robust variance threshold based on minimum sample difference
    abs_diffs <- abs(diff(x[!is.na(x)]))
    smallest_var <- suppressWarnings(min(abs_diffs[abs_diffs > 1e-5]))

    ## logical outlier positions
    abs_dev <- abs(x - local_meds)
    is_outlier <- abs_dev > smallest_var &
        abs_dev > (L * outlier_cutoff * local_mad)
    ## NAs from is_outlier check should return FALSE
    is_outlier[is.na(is_outlier)] <- FALSE

    ## return list of vectors w/ local logicals and medians
    return(list(
        local_medians = local_meds,
        is_outlier = is_outlier
    ))
}


#' @description
#' `compute_valid_neighbours()`: Compute a list of rolling window indices along
#' `x` to either side of `NA`s.
#'
#' @returns
#' `compute_valid_neighbours()`: A list the same length as the `NA` values in
#'   `x` with numeric vectors of sample indices of length `width` samples or
#'   `span` units of time `t` for valid values neighbouring split to either
#'   side of the invalid `NA`s.
#'
#' @rdname compute_helpers
#' @keywords internal
compute_valid_neighbours <- function(
    x,
    t = seq_along(x),
    width = NULL,
    span = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    na_idx <- which(is.na(x))
    valid_idx <- which(!is.na(x))
    n_valid <- length(valid_idx)
    n_na <- length(na_idx)

    if (!is.null(width)) {
        ## find position to the left of each NA in valid_idx sequence
        pos <- findInterval(na_idx, valid_idx)
        half_width <- floor(width / 2L)

        window_idx <- vector("list", n_na)
        for (i in seq_len(n_na)) {
            ## extract width samples before and after
            left <- max(1L, pos[i] - half_width + 1L):pos[i]
            right <- min(n_valid, pos[i] + 1L):min(n_valid, pos[i] + half_width)
            window_idx[[i]] <- valid_idx[sort(unique(c(left, right)))]
        }
        return(window_idx)
    }

    ## pre-compute for span approach
    t_valid <- t[valid_idx]
    t_na <- t[na_idx]
    half_span <- span * 0.5

    ## build per-NA valid neighbours with binary search on sorted `t_valid`
    ## falls back to naerest bracketing pair when no valid samples within `span`
    window_idx <- lapply(seq_len(n_na), \(.i) {
        lo <- validate_findInt(
            t_na[.i] - half_span, t_valid, left.open = TRUE, env = env
        ) + 1L
        hi <- validate_findInt(t_na[.i] + half_span, t_valid, env = env)
        if (lo > hi) {
            pos <- findInterval(na_idx[.i], valid_idx)
            return(unique(valid_idx[c(pos, min(n_valid, pos + 1L))]))
        }
        valid_idx[lo:hi]
    })

    ## window of valid values exclusive around `x`
    return(window_idx)
}


#' Preserve and restore NA information within a vector
#'
#' `preserve_na()` stores `NA` vector positions and extracts valid non-`NA`
#' values for later restoration with `restore_na()`.
#'
#' @param x A vector containing missing `NA` values.
#'
#' @returns
#' `preserve_na()` returns a list `na_info` with components:
#'   - `na_info$x_valid`: A vector with `NA` values removed.
#'   - `na_info$x_length`: A numeric value of the original input vector length.
#'   - `na_info$na_idx`: A logical vector preserving `NA` positions.
#'
#' `restore_na()` returns a vector `y` the same length as the original
#'   input vector `x` with `NA` values restored to their original positions.
#'
#' @keywords internal
preserve_na <- function(x) {
    na_info <- list(
        x_valid = x[!is.na(x)],
        x_length = length(x),
        na_idx = is.na(x)
    )
    return(na_info)
}


#' Preserve and restore NA information within a vector
#'
#' `restore_na()` restores `NA` values to their original vector positions
#' after processing valid non-`NA` values returned from `preserve_na()`.
#'
#' @param y A vector of valid non-`NA` values returned from `preserve_na()`.
#' @param na_info A list returned from `preserve_na()`.
#'
#' @rdname preserve_na
#' @keywords internal
restore_na <- function(y, na_info) {
    if (all(!na_info$na_idx)) {
        return(y)
    }
    ## fill original length of NAs
    result <- rep(NA, na_info$x_length)
    if (all(na_info$na_idx)) {
        return(result)
    }
    ## replace non-NA with processed output values
    result[!na_info$na_idx] <- y
    return(result)
}
