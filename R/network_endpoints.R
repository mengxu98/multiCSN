#' Endpoint tables of a layered multi-omic network
#'
#' The five endpoint tables written by \code{\link{fit_layered_network}} (and by
#' the manuscript pipeline) share one contract: a file name, a key column set and
#' a signed weight. This function returns that contract so writers, readers and
#' validators cannot drift apart.
#'
#' @return A data frame with columns \code{endpoint} and \code{file} plus a list
#'   column \code{keys} holding the key columns of each endpoint.
#' @export
network_endpoints <- function() {
  spec <- data.frame(
    endpoint = c("TF-gene", "TF-region", "region-gene", "factorized-chain", "mediated-TF-gene"),
    file = c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv", "mediated_tf_gene.tsv"),
    stringsAsFactors = FALSE
  )
  spec$keys <- list(
    c("regulator", "target"),
    c("regulator", "region"),
    c("region", "target"),
    c("regulator", "region", "target"),
    c("regulator", "target")
  )
  spec
}

#' Read layered network endpoints from a result directory
#'
#' @param root Directory holding the endpoint files of one fitted network.
#' @param endpoints Endpoints to read; defaults to all five.
#' @param direction Replace the signed weight column with a \code{direction}
#'   column holding its sign (finite, non-zero), which is what the transition and
#'   stability helpers consume.
#' @param weight_column Name of the signed weight column.
#' @param validate Check that key columns are present, complete and unique.
#' @return A named list of \code{data.table}s, one per endpoint.
#' @export
read_network_endpoints <- function(root, endpoints = NULL, direction = FALSE,
                                   weight_column = "weight", validate = TRUE) {
  spec <- network_endpoints()
  if (is.null(endpoints)) {
    endpoints <- spec$endpoint
  }
  unknown <- setdiff(endpoints, spec$endpoint)
  if (length(unknown)) {
    stop("Unknown endpoint(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (!dir.exists(root)) {
    stop("Result directory does not exist: ", root, call. = FALSE)
  }
  output <- vector("list", length(endpoints))
  names(output) <- endpoints
  for (index in seq_along(endpoints)) {
    endpoint <- endpoints[[index]]
    spec_row <- spec[spec$endpoint == endpoint, , drop = FALSE]
    keys <- spec_row$keys[[1L]]
    path <- file.path(root, spec_row$file)
    if (!file.exists(path)) {
      stop("Missing ", endpoint, " endpoint: ", path, call. = FALSE)
    }
    table <- data.table::fread(path)
    columns <- if (isTRUE(direction)) c(keys, weight_column) else keys
    missing <- setdiff(columns, names(table))
    if (length(missing)) {
      stop(
        endpoint, " endpoint lacks column(s): ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    if (isTRUE(validate) &&
      (anyNA(table[, keys, with = FALSE]) ||
        anyDuplicated(table[, keys, with = FALSE]))) {
      stop(endpoint, " endpoint has missing or duplicated keys.", call. = FALSE)
    }
    if (isTRUE(direction)) {
      values <- suppressWarnings(as.numeric(table[[weight_column]]))
      if (any(!is.finite(values) | values == 0)) {
        stop(endpoint, " endpoint weights must be finite and non-zero.", call. = FALSE)
      }
      table <- table[, keys, with = FALSE]
      data.table::set(table, j = "direction", value = sign(values))
    }
    output[[index]] <- table
  }
  output
}
