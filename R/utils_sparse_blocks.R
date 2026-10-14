row_compressed_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!methods::is(x, "CsparseMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "CsparseMatrix")
  }
  methods::as(x, "RsparseMatrix")
}

row_compressed_empty <- function(n_rows, n_columns) {
  methods::new("dgCMatrix",
    i = integer(), p = integer(n_columns + 1L), x = numeric(),
    Dim = c(as.integer(n_rows), as.integer(n_columns))
  )
}

row_compressed_block <- function(row_compressed, rows) {
  rows <- as.integer(rows)
  if (anyNA(rows) || any(rows < 1L | rows > nrow(row_compressed))) {
    stop("Row indices must be within the sparse matrix.", call. = FALSE)
  }
  if (!length(rows)) {
    return(row_compressed_empty(ncol(row_compressed), 0L))
  }
  if (methods::is(row_compressed, "dgRMatrix")) {
    return(gather_sparse_block(
      row_compressed@p, row_compressed@j, row_compressed@x,
      rows, ncol(row_compressed)
    ))
  }
  p <- row_compressed@p
  starts <- p[rows] + 1L
  lengths <- p[rows + 1L] - p[rows]
  if (sum(as.double(lengths)) > .Machine$integer.max) {
    stop("Sparse block exceeds the compressed matrix index limit.", call. = FALSE)
  }
  positions <- rep.int(starts, lengths) + sequence(lengths) - 1L
  block <- methods::new("dgCMatrix")
  block@i <- as.integer(row_compressed@j[positions])
  block@p <- as.integer(c(0L, cumsum(lengths)))
  block@x <- as.numeric(row_compressed@x[positions])
  block@Dim <- c(as.integer(ncol(row_compressed)), length(rows))
  block
}

column_compressed_block <- function(column_compressed, columns) {
  columns <- as.integer(columns)
  if (anyNA(columns) || any(columns < 1L | columns > ncol(column_compressed))) {
    stop("Column indices must be within the sparse matrix.", call. = FALSE)
  }
  if (methods::is(column_compressed, "dgCMatrix")) {
    return(gather_sparse_block(
      column_compressed@p, column_compressed@i, column_compressed@x,
      columns, nrow(column_compressed)
    ))
  }
  p <- column_compressed@p
  starts <- p[columns] + 1L
  lengths <- p[columns + 1L] - p[columns]
  if (sum(as.double(lengths)) > .Machine$integer.max) {
    stop("Sparse block exceeds the compressed matrix index limit.", call. = FALSE)
  }
  positions <- rep.int(starts, lengths) + sequence(lengths) - 1L
  block <- methods::new("dgCMatrix")
  block@i <- as.integer(column_compressed@i[positions])
  block@p <- as.integer(c(0L, cumsum(lengths)))
  block@x <- as.numeric(column_compressed@x[positions])
  block@Dim <- as.integer(c(nrow(column_compressed), length(columns)))
  block
}

row_compressed_nonzero <- function(row_compressed, row) {
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  keep <- row_compressed@x[positions] != 0
  as.integer(row_compressed@j[positions][keep]) + 1L
}

row_compressed_positive <- function(row_compressed, row) {
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  keep <- row_compressed@x[positions] > 0
  as.integer(row_compressed@j[positions][keep]) + 1L
}

row_compressed_dense <- function(row_compressed, row, length_out) {
  values <- numeric(length_out)
  p <- row_compressed@p
  positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
  values[row_compressed@j[positions] + 1L] <- row_compressed@x[positions]
  values
}
