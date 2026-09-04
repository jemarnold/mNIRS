#' Correct for blood volume changes
#'
#' Normalises mNIRS channels for the effects of blood volume changes, following
#' the sample-wise iterative method of *Beever & Tripp et al, 2020*.
#'
#' @param oxy_channel A character vector naming the `oxy[haem]` (oxygenated
#'   haemoglobin and myoglobin; *O2Hb*) column(s) in `data`. Must match exactly.
#' @param deoxy_channel A character vector naming the `deoxy[haem]`
#'   (deoxygenated haemoglobin and myoglobin; *HHb*) column(s) in `data`.
#'   Must match exactly.
#' @param total_channel A character vector naming the `total[haem]` (total
#'   haemoglobin and myoglobin; *THb*; proxy for blood volume) column(s) in
#'   `data`. Must match exactly.
#' @inheritParams map_mnirs_intervals
#' @inheritParams validate_mnirs
#'
#' @inheritSection map_mnirs_intervals Data input formats
#'
#' @details
#' ## Specify NIRS component channels
#'
#' At least two of `oxy_channel`, `deoxy_channel`, and `total_channel` must
#' be specified to calculate the blood volume correction factor. Best practice
#' is to specify all existing channels in `data`. Missing channels are derived
#' from the specified pair before the correction is applied.
#'
#' - `total` = `oxy + deoxy`
#' - `oxy`   = `total - deoxy`
#' - `deoxy` = `total - oxy`
#'
#' Multiple channel pairs can be corrected in one call by passing equal-length
#' vectors, with each element number forming a pair (e.g.
#' `oxy_channel = c(o2hb_1, o2hb_2), deoxy_channel = c(hhb_1, hhb_2)`).
#'
#' *NOTE*: the returned data frame will *ONLY* include corrected values for
#' the specified channels. Non-specified channels will remain uncorrected and
#' will therefore no longer be comparable to corrected channels. Best practice
#' is to specify all existing channels in `data`.
#'
#' ## Compute blood volume correction
#'
#' If any NIRS channels have negative values, all specified channels will be
#' ensemble-shifted by a common offset so that all channels contain only
#' positive values. Relative scaling across channels is preserved. This is
#' modified from the method in *Beever & Tripp et al, 2020* to properly
#' calculate `total[haem]` and the blood volume correction factor `beta` when
#' there are negative NIRS values.
#'
#' The correction factor `beta` is effectively the single-channel fractional
#' (%) oxygen saturation used to normalise `oxy[haem]` and `deoxy[haem]`
#' relative to an adjusted invariant `total[haem]`. This is computed as the
#' cumulative sum of adjusted incremental differences:
#'
#' \deqn{\Delta\text{O2Hb}_c = \Delta\text{O2Hb} - \beta \cdot \Delta\text{THb}}
#' \deqn{\Delta\text{HHb}_c = \Delta\text{HHb} - (1 - \beta) \cdot \Delta\text{THb}}
#'
#' After correction, `total[haem]` is zero (blood volume changes are
#' normalised).
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class *"mnirs"* with blood
#'   volume-corrected channels written back to the specified columns, and with
#'   metadata available with `attributes()`. For list or grouped data frame
#'   input, returns a named list of *"mnirs"* tibbles, one per interval.
#'
#' @references
#' Beever AT, Tripp TR, Zhang J, MacInnis MJ (2020) Nirs-Derived Skeletal
#' Muscle Oxidative Capacity Is Correlated with Aerobic Fitness and
#' Independent of Sex. J Appl Physiol (1985).
#' \doi{10.1152/japplphysiol.00017.2020}
#'
#' Ryan TE, Erickson ML, Brizendine JT, et al. (2012) Noninvasive Evaluation of
#' Skeletal Muscle Mitochondrial Capacity with near-Infrared Spectroscopy:
#' Correcting for Blood Volume Changes. J Appl Physiol (1985).
#' \doi{10.1152/japplphysiol.00319.2012}
#'
#' @examples
#' data <- read_mnirs(
#'     file_path = example_mnirs("artinis"),
#'     nirs_channels = c(o2hb = 2, hhb = 3),
#'     time_channel = c(sample = 1),
#'     verbose = FALSE,
#' )
#'
#' plot(data)
#'
#' result <- correct_blood_volume(
#'     data,
#'     oxy_channel = "o2hb",
#'     deoxy_channel = "hhb", ## thb will be derived from o2hb + hhb
#' )
#'
#' plot(result)
#'
#' @export
correct_blood_volume <- function(
    data,
    oxy_channel = NULL,
    deoxy_channel = NULL,
    total_channel = NULL,
    verbose = TRUE
) {
    ## list or grouped input -> normalise to named list, recurse per interval
    if (inherits(data, "grouped_df") || !is.data.frame(data)) {
        return(map_mnirs_intervals(data, match.call(), parent.frame()))
    }

    ## validation =====================================================
    validate_mnirs_data(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }

    ## collect and parse NSE channel names
    channels <- list(
        oxy = parse_channel_name(enquo(oxy_channel), data),
        deoxy = parse_channel_name(enquo(deoxy_channel), data),
        total = parse_channel_name(enquo(total_channel), data)
    )

    ## check if columns exist in data
    specified <- !vapply(channels, is.null, logical(1))

    ## require at least two channels to derive the third
    if (sum(specified) < 2L) {
        cli_abort(c(
            "x" = "At least two of {.arg oxy_channel}, {.arg deoxy_channel}, \\
            {.arg total_channel} are required."
        ))
    }

    ## channels are paired by position, so specified args must align
    if (length(unique(lengths(channels[specified]))) > 1L) {
        cli_abort(c(
            "x" = "{.arg oxy_channel}, {.arg deoxy_channel}, and \\
            {.arg total_channel} must have the same length.",
            "i" = "Channels are paired by position."
        ))
    }

    ## confirm specified column names exist in data (case-sensitive)
    missing_cols <- setdiff(unlist(channels), names(data))
    if (length(missing_cols) > 0L) {
        cli_abort(c(
            "x" = "Column{?s} {.field {missing_cols}} not found in \\
            {.arg data}.",
            "i" = "Channel names are case-sensitive and must match exactly."
        ))
    }

    ## processing ====================================================
    ## pair channels by position; ensemble shift applied within each pair
    pairs <- lapply(seq_along(channels[[which(specified)[1L]]]), \(.i) {
        vapply(channels[specified], `[[`, "", .i)
    })

    corrected <- lapply(pairs, \(.p) {
        ## pull specified channels; unspecified channels as NULL
        oxy <- if (specified[["oxy"]]) data[[.p[["oxy"]]]]
        deoxy <- if (specified[["deoxy"]]) data[[.p[["deoxy"]]]]
        total <- if (specified[["total"]]) data[[.p[["total"]]]]

        ## prefer user-specified channels, derive if NULL
        total <- total %||% (oxy + deoxy)
        oxy <- oxy %||% (total - deoxy)
        deoxy <- deoxy %||% (total - oxy)

        ## ensemble shift channels to min value == 0
        ## preserves relative scaling with all positive
        shift <- max(0, -min(oxy, deoxy, total, na.rm = TRUE)) +
            .Machine$double.eps
        oxy <- oxy + shift
        deoxy <- deoxy + shift
        total <- total + shift

        ## sample-wise blood volume correction factor
        ## (Beever & Tripp et al, 2020)
        ## dropping the first sample to align with incremental diffs
        beta <- (oxy / total)[-1L]
        diff_total <- diff(total)

        ## corrected signals are the cumulative sum of adjusted incremental
        ## diffs; first sample starts at zero
        ## total reduces to zero by construction once blood volume is
        ## normalised
        list(
            oxy = cumsum(c(0, diff(oxy) - beta * diff_total)),
            deoxy = cumsum(c(0, diff(deoxy) - (1 - beta) * diff_total)),
            total = double(length(total))
        )[specified]
    })

    ## write corrected values back to user-specified columns only
    ## metadata: over-write `nirs_channels` with user-specified col names
    nirs_channels <- unlist(pairs, use.names = FALSE)
    data[nirs_channels] <- unlist(corrected, recursive = FALSE)

    if (verbose) {
        cli_inform(c(
            "i" = "{.field {nirs_channels}} channels have been corrected \\
            for changes in blood volume."
        ))
    }

    return(create_mnirs_data(data, nirs_channels = nirs_channels))
}
