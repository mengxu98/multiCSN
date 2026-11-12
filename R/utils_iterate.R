map_list <- function(.x, .f, ...) {
  .f <- as_function(.f)
  out <- lapply(.x, .f, ...)
  names(out) <- names(.x)
  out
}

map_typed <- function(.x, .f, type, ...) {
  .f <- as_function(.f)
  out <- vapply(.x, function(value) {
    result <- .f(value, ...)
    if (length(result) != 1L) {
      stop("Result must be length 1, not ", length(result), ".", call. = FALSE)
    }
    coerce_scalar(result, type)
  }, type, USE.NAMES = FALSE)
  names(out) <- names(.x)
  out
}

coerce_scalar <- function(result, type) {
  if (is.factor(result) && !is.character(type)) {
    result <- unclass(result)
    attributes(result) <- NULL
  }
  missing_logical <- is.logical(result) && is.na(result)
  if (is.character(type)) {
    if (is.character(result) || missing_logical) {
      return(as.character(result))
    }
    stop("Can't coerce a ", typeof(result), " result to a character.", call. = FALSE)
  }
  if (is.logical(type)) {
    if (is.logical(result)) {
      return(result)
    }
    if (is.numeric(result) && (is.na(result) || unclass(result) %in% c(0, 1))) {
      return(as.logical(unclass(result)))
    }
    stop("Can't coerce a ", typeof(result), " result to a logical.", call. = FALSE)
  }
  if (is.numeric(result) || is.logical(result)) {
    return(as.numeric(unclass(result)))
  }
  stop("Can't coerce a ", typeof(result), " result to a double.", call. = FALSE)
}

map_chr_list <- function(.x, .f, ...) map_typed(.x, .f, character(1), ...)

map_lgl_list <- function(.x, .f, ...) map_typed(.x, .f, logical(1), ...)

map_dbl_list <- function(.x, .f, ...) map_typed(.x, .f, numeric(1), ...)

imap_list <- function(.x, .f, ...) {
  .f <- as_function(.f)
  index <- if (is.null(names(.x))) seq_along(.x) else names(.x)
  out <- Map(.f, .x, index, MoreArgs = list(...))
  names(out) <- names(.x)
  out
}

map2_list <- function(.x, .y, .f, ...) {
  .f <- as_function(.f)
  out <- Map(.f, .x, .y, MoreArgs = list(...))
  names(out) <- names(.x)
  out
}

bind_rows_id <- function(x, .id = NULL) {
  if (is.null(.id)) {
    return(dplyr::bind_rows(x))
  }
  dplyr::bind_rows(x, .id = .id)
}

map_dfr_list <- function(.x, .f, ..., .id = NULL) {
  bind_rows_id(map_list(.x, .f, ...), .id = .id)
}

imap_dfr_list <- function(.x, .f, ..., .id = NULL) {
  bind_rows_id(imap_list(.x, .f, ...), .id = .id)
}

pmap_dfr_list <- function(.l, .f, ..., .id = NULL) {
  .f <- as_function(.f)
  rows <- do.call(Map, c(list(.f), as.list(.l), list(MoreArgs = list(...))))
  names(rows) <- names(.l[[1L]])
  bind_rows_id(rows, .id = .id)
}

reduce_list <- function(.x, .f) {
  if (!length(.x)) {
    stop("Must supply `.init` when `.x` is empty.", call. = FALSE)
  }
  Reduce(as_function(.f), .x)
}

list_c <- function(x) {
  unlist(unname(x), recursive = FALSE, use.names = TRUE)
}

as_function <- function(.f) {
  if (inherits(.f, "formula")) {
    body <- .f[[length(.f)]]
    environment <- environment(.f)
    return(eval(call("function", as.pairlist(alist(.x = , ... = )), body), environment))
  }
  .f
}
