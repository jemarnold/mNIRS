#' Resolve per-channel arguments
#'
#' Broadcasts global argument values across `nirs_channels`, applying
#' per-channel overrides where an argument is supplied as a named `list()`
#' keyed by channel name. An argument is treated as per-channel only when
#' it is a `list()` whose named elements all match `nirs_channels` or group
#' names, with at most one unnamed element acting as the fallback for
#' unlisted channels (e.g. `width = list(5, q = 7)` gives `q` 7 and every
#' other channel 5). Any other value (including unnamed vectors and named
#' lists such as `control = list(maxiter = 100)`) is applied globally to
#' every channel.
#'
#' @param nirs_channels Character vector of resolved channel names.
#' @param args Named list of per-channel-capable arguments. Each element is
#'   either a global value or a per-channel `list()` map. A per-channel map
#'   may include a single unnamed element as the fallback for unlisted
#'   channels.
#' @param defaults Named list of fallback values per argument, used when a
#'   per-channel map omits a channel and supplies no unnamed fallback. Only
#'   needed for arguments whose formal default is not `NULL` (e.g.
#'   `method = "linear"`).
#' @param choices Named list of valid values for choice-type arguments
#'   (e.g. `list(method = c("linear", "median", "locf", "none"))`).
#'   Resolved values are matched per channel; a full default vector
#'   resolves to its first element, matching [match.arg()] behaviour.
#' @param group_channels An *optional* named list of channel-name vectors
#'   from [validate_group_channels()]. When supplied, arguments are
#'   resolved per group: a group-name key or any member-channel key
#'   applies to the whole group, and conflicting member values within one
#'   group abort.
#' @param call The calling environment, used to report errors as coming
#'   from the user-facing function (e.g. [rescale_mnirs()]).
#' @inheritParams validate_mnirs
#'
#' @returns A named list with one element per channel (or per group when
#'   `group_channels` is supplied); each element is a named list of that
#'   channel's resolved argument values.
#'
#' @keywords internal
resolve_channel_args <- function(
    nirs_channels,
    args,
    defaults = list(),
    choices = list(),
    group_channels = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## units of resolution: single channels, or groups when supplied
    units <- group_channels %||% setNames(as.list(nirs_channels), nirs_channels)
    valid_ch <- c(nirs_channels, names(units))

    ## a list is a per-channel map only when at least one element names a
    ## known channel or group, every named element is a known name, and at
    ## most one element is unnamed (the fallback for unlisted channels)
    is_per_channel <- function(.a) {
        if (!is.list(.a) || is.data.frame(.a)) {
            return(FALSE)
        }
        named <- nzchar(names(.a))
        any(named) && sum(!named) <= 1L && all(names(.a)[named] %in% valid_ch)
    }

    ## warn when a list arg looks like a mistyped per-channel map
    if (verbose) {
        invisible(lapply(names(args), \(.nm) {
            .a <- args[[.nm]]
            names <- names(.a)
            if (
                is.list(.a) && !is.data.frame(.a) && length(names) > 0L &&
                    !is_per_channel(.a) && any(names %in% valid_ch)
            ) {
                unknown <- setdiff(names[nzchar(names)], valid_ch)
                cli_warn(c(
                    "!" = "{.arg {(.nm)}} treated as a global value: \\
                    name{?s} {.val {unknown}} not found in \\
                    {.arg nirs_channels}.",
                    "i" = "Per-channel lists must be named by \\
                    {.arg nirs_channels}, with only one unnamed default value."
                ), call = env)
            }
        }))
    }

    ## inform when a per-channel list omits channels and gives no unnamed
    ## fallback; those channels fall back to the argument's formal default
    if (verbose) {
        invisible(lapply(names(args), \(.nm) {
            .a <- args[[.nm]]
            if (!is_per_channel(.a) || any(!nzchar(names(.a)))) {
                return(invisible())
            }
            omitted <- names(units)[vapply(names(units), \(.key) {
                length(intersect(names(.a), c(.key, units[[.key]]))) == 0L
            }, logical(1))]
            if (length(omitted) > 0L) {
                cli_inform(c(
                    "i" = "{.arg {(.nm)}}: channel{?s} {.val {omitted}} \\
                    not listed and will not be processed."
                ))
            }
        }))
    }

    ## resolve one argument for one channel/group
    resolve_one <- function(.a, .nm, .key, .members) {
        if (is_per_channel(.a)) {
            ## the lone unnamed element is the fallback for unlisted channels
            unnamed <- .a[!nzchar(names(.a))]
            default_val <- if (length(unnamed)) unnamed[[1L]] else NULL
            ## group-name key preferred, then member-channel keys
            hits <- .a[intersect(names(.a), c(.key, .members))]
            if (length(hits) > 1L && length(unique(hits)) > 1L) {
                cli_abort(c(
                    "x" = "{.arg {(.nm)}} has conflicting values within \\
                    {.arg group_channels} = {.val {(.key)}}.",
                    "i" = "Grouped channels must share one value per \\
                    argument."
                ), call = env)
            }
            .a <- if (length(hits) > 0L) {
                hits[[1L]] %||% default_val %||% defaults[[.nm]]
            } else {
                default_val %||% defaults[[.nm]]
            }
        }
        ## match choice-type args; a full default vector resolves to its
        ## first element, matching `match.arg()` behaviour
        if (!is.null(choices[[.nm]]) && !is.null(.a)) {
            if (identical(.a, choices[[.nm]])) {
                .a <- .a[[1L]]
            }
            .a <- rlang::arg_match0(
                .a, choices[[.nm]], arg_nm = .nm, error_call = env
            )
        }
        return(.a)
    }

    out <- lapply(names(units), \(.key) {
        setNames(
            lapply(names(args), \(.nm) {
                resolve_one(args[[.nm]], .nm, .key, units[[.key]])
            }),
            names(args)
        )
    })
    return(setNames(out, names(units)))
}


#' Validate and normalise channel grouping
#'
#' Converts the `group_channels` argument to a named list of channel-name
#' vectors. String shortcuts expand against `nirs_channels`: `"ensemble"`
#' places all channels in one group (preserving relative scaling) and
#' `"distinct"` places each channel in its own group. Custom `list()`
#' groupings may use bare symbols or character names; channels omitted
#' from a custom grouping are processed independently, matching
#' `group_intervals` behaviour in [extract_intervals()].
#'
#' @param nirs_channels Character vector of resolved channel names.
#' @param group_channels A quosure from `rlang::enquo()`, a character
#'   string (`"ensemble"` or `"distinct"`), or a `list()` of channel-name
#'   vectors. Lists may be named; unnamed groups are keyed by their first
#'   member.
#' @param data A data frame for parsing bare-symbol group members.
#' @param env Environment for symbol evaluation.
#' @inheritParams validate_mnirs
#'
#' @returns A named list of character vectors covering all
#'   `nirs_channels`, each channel appearing in exactly one group.
#'
#' @keywords internal
validate_group_channels <- function(
    nirs_channels,
    group_channels,
    data = NULL,
    verbose = TRUE,
    env = rlang::caller_env()
) {
    ## parse tidy eval input, preserving group names from list() calls
    if (rlang::is_quosure(group_channels)) {
        expr <- rlang::quo_get_expr(group_channels)
        if (rlang::is_call(expr, "list")) {
            quo_env <- rlang::quo_get_env(group_channels)
            group_channels <- lapply(rlang::call_args(expr), \(.arg) {
                parse_channel_name(
                    rlang::new_quosure(.arg, env = quo_env), data, env
                )
            })
        } else {
            group_channels <- parse_channel_name(group_channels, data, env)
        }
    }

    ## string shortcuts expand against nirs_channels
    if (is.character(group_channels)) {
        shortcut <- rlang::arg_match0(
            arg = group_channels[[1L]],
            values = c("ensemble", "distinct"),
            arg_nm = "group_channels"
        )
        group_channels <- if (shortcut == "ensemble") {
            list(nirs_channels)
        } else {
            as.list(nirs_channels)
        }
    }

    if (!is.list(group_channels)) {
        cli_abort(c(
            "x" = "{.arg group_channels} must be {.val ensemble}, \\
            {.val distinct}, or a {.cls list} of channel names."
        ), call = env)
    }

    ## group members must be known channels
    members <- unlist(group_channels, use.names = FALSE)
    unknown <- setdiff(members, nirs_channels)
    if (!is.character(members) || length(unknown) > 0L) {
        cli_abort(c(
            "x" = "{.arg group_channels} contains unknown \\
            channel{?s}: {.val {unknown}}.",
            "i" = "Group members must match {.arg nirs_channels} exactly."
        ), call = env)
    }

    ## each channel may belong to one group only
    if (anyDuplicated(members) > 0L) {
        dupes <- unique(members[duplicated(members)])
        cli_abort(c(
            "x" = "Channel{?s} {.val {dupes}} assigned to more than one \\
            group in {.arg group_channels}.",
            "i" = "Each channel may belong to one group only."
        ), call = env)
    }

    ## channels omitted from custom groups are processed independently,
    ## matching `group_intervals` behaviour for unspecified intervals
    group_channels <- c(
        group_channels, as.list(setdiff(nirs_channels, members))
    )

    ## name unnamed groups by their first member
    names <- names(group_channels) %||% rep("", length(group_channels))
    first_members <- vapply(group_channels, `[[`, "", 1L)
    names[!nzchar(names)] <- first_members[!nzchar(names)]

    return(setNames(group_channels, names))
}
