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
  if (!length(rows)) {
    return(row_compressed_empty(ncol(row_compressed), 0L))
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

column_compressed_block <- function(column_compressed, columns) {
  columns <- as.integer(columns)
  p <- column_compressed@p
  starts <- p[columns] + 1L
  lengths <- p[columns + 1L] - p[columns]
  positions <- rep.int(starts, lengths) + sequence(lengths) - 1L
  methods::new("dgCMatrix",
    i = as.integer(column_compressed@i[positions]),
    p = as.integer(c(0L, cumsum(lengths))),
    x = as.numeric(column_compressed@x[positions]),
    Dim = as.integer(c(nrow(column_compressed), length(columns)))
  )
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
