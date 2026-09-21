# Internal helpers for cost-safe sparse row blocks.
#
# `x[rows, , drop = FALSE]` on a *C*ompressed sparse matrix costs O(nnz) per
# call, because the requested rows are scattered across every column. The
# layered fit performs one such subset per target -- tens of thousands of times
# per layer -- so this single pattern dominates the runtime at real scale.
# Reading the same entries from the row-compressed slots costs O(nnz of the
# requested rows) and returns an identical block.

.row_compressed_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!methods::is(x, "CsparseMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "CsparseMatrix")
  }
  methods::as(x, "RsparseMatrix")
}

.row_compressed_empty <- function(n_rows, n_columns) {
  methods::new("dgCMatrix",
    i = integer(), p = integer(n_columns + 1L), x = numeric(),
    Dim = c(as.integer(n_rows), as.integer(n_columns))
  )
}

#' Extract columns of selected rows from a row-compressed sparse matrix
#'
#' Returns the identical block as `t(x[rows, , drop = FALSE])`.
#' @noRd
.row_compressed_block <- function(row_compressed, rows) {
  rows <- as.integer(rows)
  if (!length(rows)) {
    return(.row_compressed_empty(ncol(row_compressed), 0L))
  }
  p <- row_compressed@p
  starts <- p[rows] + 1L
  lengths <- p[rows + 1L] - p[rows]
  positions <- rep.int(starts, lengths) + sequence(lengths) - 1L
  methods::new("dgCMatrix",
    i = as.integer(row_compressed@j[positions]),
    p = as.integer(c(0L, cumsum(lengths))),
    x = as.numeric(row_compressed@x[positions]),
    Dim = c(as.integer(ncol(row_compressed)), length(rows))
  )
}

#' Column indices of the stored entries of one row
#' @noRd
.row_compressed_nonzero <- function(row_compressed, row) {
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  keep <- row_compressed@x[positions] != 0
  as.integer(row_compressed@j[positions][keep]) + 1L
}

#' Column indices of the positive entries of one row
#' @noRd
.row_compressed_positive <- function(row_compressed, row) {
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  keep <- row_compressed@x[positions] > 0
  as.integer(row_compressed@j[positions][keep]) + 1L
}

#' Dense representation of one row (`as.numeric(x[row, ])`)
#' @noRd
.row_compressed_dense <- function(row_compressed, row, length_out) {
  values <- numeric(length_out)
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  values[row_compressed@j[positions] + 1L] <- row_compressed@x[positions]
  values
}
