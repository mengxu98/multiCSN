#' Endpoint tables of a layered multi-omic network
#'
#' The five endpoint tables written by \code{\link{fit_layered_network}} (and by
#' the manuscript pipeline) share one contract: a file name, a key column set and
#' a signed weight. This function returns that contract so writers, readers and
#' validators cannot drift apart.
#'
#' @return A data frame with columns \code{endpoint} and \code{file} plus a list
#'   column \code{keys} holding the key columns of each endpoint.
#' @export network_endpoints write_network_endpoints
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

write_network_endpoints <- function(fit, root,
                                    format = c("tsv", "csv", "tsv.gz", "csv.gz"),
                                    overwrite = FALSE) {
  format <- match.arg(format)
  spec <- network_endpoints()
  tables <- c(
    fit[sub("[.]tsv$", "", spec$file)],
    list(
      tf_gene_accounting = fit$accounting$tf_gene,
      tf_region_accounting = fit$accounting$tf_region,
      region_gene_accounting = fit$accounting$region_gene,
      endpoint_summary = fit$summary
    )
  )
  if (any(vapply(tables, is.null, logical(1)))) {
    stop("Incomplete layered fit.", call. = FALSE)
  }
  for (i in seq_len(nrow(spec))) {
    check_endpoint_keys(tables[[i]], spec$keys[[i]], spec$endpoint[[i]])
  }
  tables$tf_region_accounting <- as.data.frame(tables$tf_region_accounting)
  names(tables$tf_region_accounting)[names(tables$tf_region_accounting) == "region"] <- "target"
  paths <- file.path(root, paste0(names(tables), ".", format))
  if (!isTRUE(overwrite) && any(file.exists(paths))) {
    stop("Output files already exist.", call. = FALSE)
  }
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(tables)) {
    write_endpoint_table(tables[[i]], paths[[i]], format)
  }
  invisible(stats::setNames(paths, names(tables)))
}

#' Read layered network endpoints from a result directory
#'
#' @param root Directory holding the endpoint files of one fitted network.
#' @param endpoints Endpoints to read; defaults to all five.
#' @param direction Replace the signed weight column with a \code{direction}
#'   column holding its sign (finite, non-zero).
#' @param format File format: tsv, csv, tsv.gz or csv.gz.
#' @param weight_column Name of the signed weight column.
#' @param validate Check that key columns are present, complete and unique.
#' @return A named list of \code{data.table}s, one per endpoint.
#' @export
read_network_endpoints <- function(root, endpoints = NULL, direction = FALSE,
                                   weight_column = "weight", validate = TRUE,
                                   format = c("tsv", "csv", "tsv.gz", "csv.gz")) {
  format <- match.arg(format)
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
    path <- file.path(root, sub("[.]tsv$", paste0(".", format), spec_row$file))
    sep <- if (startsWith(format, "csv")) "," else "\t"
    table <- read_endpoint_table(path, sep, keys)
    columns <- if (isTRUE(direction)) c(keys, weight_column) else keys
    check_endpoint_keys(table, keys, endpoint, columns, validate)
    if (isTRUE(direction)) {
      values <- suppressWarnings(as.numeric(table[[weight_column]]))
      if (any(!is.finite(values) | values == 0)) {
        stop(endpoint, " endpoint weights must be finite and non-zero.", call. = FALSE)
      }
      table <- table[, keys, with = FALSE]
      data.table::setDT(table)
      data.table::set(table, j = "direction", value = sign(values))
    }
    output[[index]] <- table
  }
  output
}

check_endpoint_keys <- function(table, keys, endpoint, columns = keys, validate = TRUE) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    stop(endpoint, " endpoint lacks column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (isTRUE(validate)) {
    key_table <- as.data.frame(table)[, keys, drop = FALSE]
    if (anyNA(key_table) || anyDuplicated(key_table)) {
      stop(endpoint, " endpoint has missing or duplicated keys.", call. = FALSE)
    }
  }
  invisible(NULL)
}

write_endpoint_table <- function(table, path, format) {
  sep <- if (startsWith(format, "csv")) "," else "\t"
  data.table::fwrite(table, path, sep = sep, quote = "auto", na = "NA",
                    compress = if (endsWith(format, ".gz")) "gzip" else "none")
  invisible(path)
}

read_endpoint_table <- function(path, sep = "\t", keys = character()) {
  data.table::as.data.table(utils::read.delim(
    path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = stats::setNames(rep("character", length(keys)), keys)
  ))
}
