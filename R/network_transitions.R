#' Classify edges between two network states
#'
#' Merges two endpoint tables on their key columns and labels every edge as
#' gained, lost, sign-stable or sign-flipped. This is the shared reduction behind
#' the state-transition, donor-deletion and cell-budget summaries.
#'
#' @param previous,next_state Endpoint tables with the key columns and a signed
#'   direction column.
#' @param keys Key columns identifying an edge.
#' @param direction_column Column holding the signed direction in both tables.
#' @param labels Named character vector with the class names for
#'   \code{gained}, \code{lost}, \code{stable} and \code{flip}.
#' @return A \code{data.table} with the key columns,
#'   \code{previous_direction}, \code{next_direction} and \code{topology_class};
#'   edges present in only one state have \code{NA} on the other side.
#' @export
classify_edge_transitions <- function(previous, next_state, keys,
                                      direction_column = "direction",
                                      labels = c(
                                        gained = "gained", lost = "lost",
                                        stable = "shared_sign_stable",
                                        flip = "shared_sign_flip"
                                      )) {
  if (!all(c("gained", "lost", "stable", "flip") %in% names(labels))) {
    stop("labels must name gained, lost, stable and flip.", call. = FALSE)
  }
  prepare <- function(table, label) {
    if (!is.data.frame(table)) {
      stop(label, " must be a data frame or data.table.", call. = FALSE)
    }
    table <- data.table::as.data.table(data.table::copy(table))
    missing <- setdiff(c(keys, direction_column), names(table))
    if (length(missing)) {
      stop(label, " lacks column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    table <- table[, c(keys, direction_column), with = FALSE]
    if (anyNA(table[, keys, with = FALSE]) || anyDuplicated(table[, keys, with = FALSE])) {
      stop(label, " has missing or duplicated edge keys.", call. = FALSE)
    }
    table
  }
  previous <- prepare(previous, "previous")
  next_state <- prepare(next_state, "next_state")
  data.table::setnames(previous, direction_column, "previous_direction")
  data.table::setnames(next_state, direction_column, "next_direction")
  merged <- merge(previous, next_state, by = keys, all = TRUE, sort = FALSE)
  data.table::set(
    merged, j = "topology_class",
    value = data.table::fifelse(
      is.na(merged$previous_direction), labels[["gained"]],
      data.table::fifelse(
        is.na(merged$next_direction), labels[["lost"]],
        data.table::fifelse(
          merged$previous_direction == merged$next_direction,
          labels[["stable"]], labels[["flip"]]
        )
      )
    )
  )
  merged[]
}
