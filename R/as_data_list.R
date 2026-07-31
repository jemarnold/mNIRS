#' Coerce `data` input to a named list of data frames
#' @inheritParams validate_mnirs
#' @keywords internal
as_data_list <- function(data, env = rlang::caller_env()) {
    ## grouped data frame → split by groups
    if (inherits(data, "grouped_df")) {
        if (!requireNamespace("dplyr", quietly = TRUE)) {
            cli_abort(c(
                "x" = "{.pkg dplyr} is required for grouped data frame input.",
                "i" = "Install with {.code install.packages(\"dplyr\")}."
            ), call = env)
        }

        ## refactor grouping variables to order of appearance
        group_vars <- dplyr::group_vars(data)
        df_grp <- data |>
            dplyr::ungroup() |>
            dplyr::mutate(
                dplyr::across(dplyr::all_of(group_vars), \(.x) {
                    factor(.x, levels = unique(.x))
                })
            ) |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)))
        keys <- do.call(paste, c(dplyr::group_keys(df_grp), list(sep = "_")))
        data_list <- dplyr::group_split(df_grp, .keep = TRUE)

        ## copy mnirs metadata down to each df in the list
        data_list <- lapply(data_list, \(.df) {
            create_mnirs_data(.df, attributes(data)[mnirs_metadata])
        })
        names(data_list) <- keys

        return(data_list)
    }

    ## single data frame → length-1 list
    if (is.data.frame(data)) {
        return(setNames(list(data), "interval_1"))
    }

    ## list of data frames — validate
    if (!is.list(data) || !all(vapply(data, is.data.frame, logical(1)))) {
        cli_abort(
            "{.arg data} must be a list of data frames, or a single grouped \\
            or ungrouped data frame.",
            call = env
        )
    }

    ## name any unnamed intervals `interval_<n>`, keeping supplied names
    nms <- names(data) %||% rep("", length(data))
    empty <- !nzchar(nms)
    nms[empty] <- paste0("interval_", seq_along(data))[empty]
    names(data) <- nms

    return(data)
}


#' Apply an mnirs function over each interval of a multi-interval input
#'
#' Shared entry point for `*_mnirs()` transformer functions accepting a list
#' of data frames or a grouped data frame. Normalises `data` to a named list
#' via [as_data_list()], then re-evaluates the captured user-facing call once
#' per interval with `data` swapped, so all arguments (including NSE channel
#' expressions) are forwarded verbatim.
#'
#' @param data A data frame of class *"mnirs"* containing time series data
#'   and metadata, a list of data frames, or a grouped data frame (see
#'   *Details*).
#' @param call The matched call from the user-facing function, re-evaluated
#'   with `data` swapped for each interval data frame.
#' @param eval_env Environment in which to re-evaluate `call`, i.e. the
#'   user-facing function's `parent.frame()`, so NSE arguments resolve
#'   against the original caller.
#' @inheritParams validate_mnirs
#'
#' @section Data input formats:
#'
#' *{mnirs}* processing functions accept `data` in multiple formats:
#'
#' - A **single *"mnirs"* data frame** is processed and returned directly.
#' - A **list of *"mnirs"* data frames**: each interval is processed
#'   separately and returned as a named list.
#' - A **grouped *"mnirs"* data frame**, e.g. with `dplyr::group_by()`: the
#'   data frame is split by grouping levels and each group is processed as
#'   a separate interval, returned as a named list.
#'
#' @returns A named list of processed *"mnirs"* data frames, one per
#'   interval.
#'
#' @keywords internal
map_mnirs_intervals <- function(
    data,
    call,
    eval_env,
    env = rlang::caller_env()
) {
    data_list <- as_data_list(data, env = env)
    return(lapply(data_list, \(.df) {
        call$data <- .df
        eval(call, envir = eval_env)
    }))
}
