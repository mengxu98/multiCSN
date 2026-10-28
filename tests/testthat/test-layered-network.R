test_that("row-compressed blocks reproduce sparse row subsetting", {
  set.seed(11)
  x <- Matrix::rsparsematrix(40, 12, density = 0.25)
  dimnames(x) <- list(paste0("r", seq_len(40)), paste0("c", seq_len(12)))
  x[3, ] <- 0
  x[5, 1] <- 0
  row_compressed <- multiCSN:::row_compressed_matrix(x)
  rows <- c(2L, 3L, 5L, 17L)
  reference <- Matrix::t(x[rows, , drop = FALSE])
  block <- multiCSN:::row_compressed_block(row_compressed, rows)
  expect_identical(block@i, reference@i)
  expect_identical(block@p, reference@p)
  expect_identical(block@x, reference@x)
  expect_identical(block@Dim, reference@Dim)
  expect_identical(
    multiCSN:::row_compressed_nonzero(row_compressed, 3L), integer()
  )
  expect_identical(
    multiCSN:::row_compressed_nonzero(row_compressed, 5L),
    as.integer(which(x[5, ] != 0))
  )
  expect_equal(
    multiCSN:::row_compressed_dense(row_compressed, 7L, ncol(x)),
    as.numeric(x[7, ])
  )
  expect_equal(dim(multiCSN:::row_compressed_block(row_compressed, integer())), c(12L, 0L))
})

test_that("column-compressed blocks reproduce sparse column subsetting", {
  set.seed(13)
  x <- Matrix::rsparsematrix(25, 18, density = 0.2)
  dimnames(x) <- list(paste0("r", seq_len(25)), paste0("c", seq_len(18)))
  x[, 3] <- 0
  for (columns in list(integer(), 3L, 1:5, c(14L, 2L, 14L, 18L))) {
    reference <- x[, columns, drop = FALSE]
    block <- multiCSN:::column_compressed_block(x, columns)
    expect_identical(block@i, reference@i)
    expect_identical(block@p, reference@p)
    expect_identical(block@x, reference@x)
    expect_identical(block@Dim, reference@Dim)
    design <- matrix(rnorm(25 * 4), nrow = 25)
    expect_identical(
      unname(as.matrix(Matrix::crossprod(design, block))),
      unname(as.matrix(Matrix::crossprod(design, reference)))
    )
  }
})

test_that("sparse block assembly preserves slots and rejects invalid indices", {
  set.seed(1014)
  x <- methods::as(Matrix::rsparsematrix(70, 80, density = 0.1), "dgCMatrix")
  x[, 3] <- 0
  x[2, 2] <- NA_real_
  x[4, 5] <- Inf
  rows <- multiCSN:::row_compressed_matrix(x)
  for (selected in list(integer(), 3L, c(65L, 3L, 2L, 4L, 65L, 1L))) {
    for (by_row in c(FALSE, TRUE)) {
      source <- if (by_row) rows else x
      indices <- if (by_row) source@j else source@i
      first <- source@p[selected] + 1L
      sizes <- source@p[selected + 1L] - source@p[selected]
      positions <- rep.int(first, sizes) + sequence(sizes) - 1L
      expected <- methods::new("dgCMatrix",
        i = as.integer(indices[positions]),
        p = as.integer(c(0L, cumsum(sizes))),
        x = as.numeric(source@x[positions]),
        Dim = as.integer(c(if (by_row) ncol(source) else nrow(source), length(selected)))
      )
      actual <- if (by_row) {
        multiCSN:::row_compressed_block(source, selected)
      } else {
        multiCSN:::column_compressed_block(source, selected)
      }
      expect_identical(actual, expected)
      expect_true(methods::validObject(actual))
    }
  }
  for (bad in list(NA_integer_, 0L, -1L, 81L)) {
    expect_error(multiCSN:::row_compressed_block(rows, bad), "Row indices")
    expect_error(multiCSN:::column_compressed_block(x, bad), "Column indices")
  }
})

test_that("layered candidate construction matches the reference selection", {
  set.seed(23)
  peaks <- paste0("p", seq_len(12))
  genes <- c("g1", "g2", "g3")
  cells <- paste0("c", seq_len(30))
  peak_tf_gate <- Matrix::Matrix(
    matrix(as.numeric(runif(length(peaks) * 3) > 0.6),
      nrow = length(peaks),
      dimnames = list(peaks, c("tf1", "tf2", "tf3"))
    ),
    sparse = TRUE
  )
  peaks2gene <- Matrix::Matrix(
    matrix(rep(c(1, 0, 1, 0), length.out = length(peaks) * length(genes)),
      nrow = length(genes), dimnames = list(genes, peaks)
    ),
    sparse = TRUE
  )
  gene_by_cell <- Matrix::Matrix(
    matrix(rpois(length(genes) * length(cells), 1),
      nrow = length(genes),
      dimnames = list(genes, cells)
    ),
    sparse = TRUE
  )
  candidate <- multiCSN:::layered_candidate_functions(
    gene_by_cell = gene_by_cell, peaks2gene = peaks2gene,
    peak_tf_gate = peak_tf_gate, features = genes, regulators = colnames(peak_tf_gate),
    min_detected = 0L, endpoint = "tf_gene"
  )
  for (index in seq_along(genes)) {
    target <- genes[[index]]
    regions <- colnames(peaks2gene)[as.logical(peaks2gene[target, ])]
    regions <- intersect(regions, rownames(peak_tf_gate))
    expected <- if (length(regions)) {
      which(Matrix::colSums(peak_tf_gate[regions, , drop = FALSE] != 0) > 0)
    } else {
      integer()
    }
    expected <- setdiff(expected, match(target, colnames(peak_tf_gate)))
    expect_identical(candidate(index), as.integer(expected))
  }
  region_candidate <- multiCSN:::layered_candidate_functions(
    gene_by_cell = gene_by_cell, peaks2gene = peaks2gene,
    peak_tf_gate = peak_tf_gate, features = genes, regulators = colnames(peak_tf_gate),
    min_detected = 0L, endpoint = "tf_region"
  )
  for (row in seq_len(nrow(peak_tf_gate))) {
    expect_identical(
      region_candidate(row),
      as.integer(which(as.numeric(peak_tf_gate[row, ]) > 0))
    )
  }
})

test_that("region-gene layer matches the naive per-target reference", {
  set.seed(31)
  peaks <- paste0("p", seq_len(9))
  genes <- paste0("g", seq_len(4))
  cells <- paste0("c", seq_len(40))
  set.seed(7)
  peak_by_cell <- Matrix::Matrix(
    matrix(rpois(length(peaks) * length(cells), 0.8),
      nrow = length(peaks),
      dimnames = list(peaks, cells)
    ),
    sparse = TRUE
  )
  gene_by_cell <- Matrix::Matrix(
    matrix(rpois(length(genes) * length(cells), 1.5),
      nrow = length(genes),
      dimnames = list(genes, cells)
    ),
    sparse = TRUE
  )
  peaks2gene <- Matrix::Matrix(
    matrix(rep(c(1, 0, 1), length.out = length(peaks) * length(genes)),
      nrow = length(genes), dimnames = list(genes, peaks)
    ),
    sparse = TRUE
  )
  settings <- list(
    response_chunk = 1024L, cores = 1L, max_support_size = NULL,
    min_improvement = 1e-10, verbose = FALSE
  )
  fitted <- multiCSN:::fit_layered_region_gene(
    peak_by_cell, gene_by_cell, genes, peaks2gene, settings
  )
  naive <- lapply(genes, function(target) {
    candidates <- colnames(peaks2gene)[as.logical(peaks2gene[target, ])]
    candidates <- intersect(candidates, rownames(peak_by_cell))
    empty <- data.frame(
      region = character(), target = character(), standardized_beta = numeric(),
      deletion_delta_bic = numeric(), stringsAsFactors = FALSE
    )
    if (!length(candidates)) {
      return(empty)
    }
    n_obs <- ncol(peak_by_cell)
    x <- Matrix::t(peak_by_cell[candidates, , drop = FALSE])
    means <- as.numeric(Matrix::colSums(x)) / n_obs
    centered_ss <- pmax(0, as.numeric(Matrix::colSums(x^2)) - n_obs * means^2)
    scales <- sqrt(centered_ss / (n_obs - 1L))
    variable <- is.finite(scales) & scales > 0
    candidates <- candidates[variable]
    if (!length(candidates)) {
      return(empty)
    }
    x <- x[, variable, drop = FALSE]
    means <- means[variable]
    scales <- scales[variable]
    gram <- (as.matrix(Matrix::crossprod(x)) - n_obs * tcrossprod(means)) /
      tcrossprod(scales)
    gram <- (gram + t(gram)) / 2
    diag(gram) <- n_obs - 1
    y <- as.numeric(gene_by_cell[target, , drop = TRUE])
    y_mean <- mean(y)
    y_scale <- sqrt(sum((y - y_mean)^2) / (n_obs - 1L))
    if (!is.finite(y_scale) || y_scale <= 0) {
      return(empty)
    }
    xty <- matrix(
      (as.numeric(Matrix::crossprod(x, y)) - n_obs * means * y_mean) /
        (scales * y_scale),
      ncol = 1L
    )
    fit <- inferCSN::fit_greedy_l0_batch(
      gram = gram, xty = xty, response_ss = n_obs - 1,
      candidates = list(seq_along(candidates)), n_obs = n_obs,
      max_support_size = NULL, min_improvement = 1e-10
    )
    if (!length(fit$predictor_index)) {
      return(empty)
    }
    data.frame(
      region = candidates[fit$predictor_index], target = target,
      standardized_beta = as.numeric(fit$standardized_beta),
      deletion_delta_bic = as.numeric(fit$deletion_delta_bic),
      stringsAsFactors = FALSE
    )
  })
  reference <- do.call(rbind, naive)
  rownames(reference) <- NULL
  expect_equal(fitted$edges, reference)
  expect_identical(
    fitted$accounting$n_selected,
    as.integer(c(
      nrow(naive[[1]]), nrow(naive[[2]]), nrow(naive[[3]]), nrow(naive[[4]])
    ))
  )
  expect_identical(fitted$accounting$target, genes)
})

test_that("checkpoints round-trip without changing the layered fit", {
  set.seed(41)
  cells <- paste0("c", seq_len(20))
  responses <- paste0("t", seq_len(6))
  design <- matrix(rnorm(length(cells) * 3),
    nrow = length(cells),
    dimnames = list(cells, paste0("tf", seq_len(3)))
  )
  design <- sweep(design, 2L, colMeans(design), "-")
  design <- sweep(design, 2L, sqrt(colSums(design^2) / (nrow(design) - 1L)), "/")
  response <- Matrix::Matrix(
    matrix(rpois(length(responses) * length(cells), 2),
      nrow = length(responses),
      dimnames = list(responses, cells)
    ),
    sparse = TRUE
  )
  settings <- list(
    response_chunk = 2L, cores = 1L, exclude_response_alias = FALSE,
    max_support_size = NULL, min_improvement = 1e-10
  )
  candidates <- function(index) seq_len(ncol(design))
  cache <- file.path(tempdir(), "layered-checkpoints")
  unlink(cache, recursive = TRUE)
  first <- multiCSN:::fit_layered_shared_design(
    design, Matrix::t(response), responses, candidates, "TF-gene", settings, cache
  )
  expect_true(length(list.files(cache, recursive = TRUE)) > 0L)
  second <- multiCSN:::fit_layered_shared_design(
    design, Matrix::t(response), responses, candidates, "TF-gene", settings, cache
  )
  expect_identical(first$edges, second$edges)
  expect_identical(first$accounting, second$accounting)
})

test_that("checkpoint identity rejects changed inputs and legacy parts", {
  root <- tempfile("layered-identity-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  inputs <- list(
    cells = "c1", features = "g1", regulators = "tf1",
    tf_region_names = "p1", gene_by_cell = matrix(1, 1, 1),
    all_peak_by_cell = matrix(1, 1, 1),
    peaks2gene = matrix(1, 1, 1), peak_tf_gate = matrix(1, 1, 1)
  )
  settings <- list(min_detected = 0L, cores = 1L, checkpoint_dir = root)
  multiCSN:::layered_prepare_checkpoint(root, inputs, settings)
  expect_silent(multiCSN:::layered_prepare_checkpoint(root, inputs, settings))
  parallel_settings <- settings
  parallel_settings$cores <- 2L
  expect_silent(multiCSN:::layered_prepare_checkpoint(root, inputs, parallel_settings))
  changed <- inputs
  changed$gene_by_cell[1, 1] <- 2
  expect_error(
    multiCSN:::layered_prepare_checkpoint(root, changed, settings),
    "different input"
  )
  changed_settings <- settings
  changed_settings$min_detected <- 20L
  expect_error(
    multiCSN:::layered_prepare_checkpoint(root, inputs, changed_settings),
    "different input"
  )
  unlink(file.path(root, "layered_checkpoint_identity.rds"))
  dir.create(file.path(root, "tf_gene_parts"))
  saveRDS(list(edges = data.frame()), file.path(root, "tf_gene_parts", "chunk_000001_000001.rds"))
  expect_error(
    multiCSN:::layered_prepare_checkpoint(root, inputs, settings),
    "no identity"
  )
})

test_that("parallel worker errors retain the original message", {
  expect_error(
    multiCSN:::layered_lapply(1:2, function(i) stop("target failed"), cores = 2L),
    "target failed"
  )
})

test_that("PSOCK workers return ordered results", {
  previous <- getOption("multicsn.parallel_backend")
  on.exit(options(multicsn.parallel_backend = previous), add = TRUE)
  options(multicsn.parallel_backend = "psock")
  observed <- multiCSN:::layered_lapply(1:3, function(i) i * i, cores = 2L)
  expect_identical(unname(observed), list(1L, 4L, 9L))
})

test_that("fit_layered_network validates its inputs", {
  x <- Matrix::Matrix(matrix(1:12, 3, 4, dimnames = list(paste0("g", 1:3), paste0("c", 1:4))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(x)
  expect_error(
    fit_layered_network(object, verbose = FALSE),
    "chromatin accessibility assay"
  )
  expect_error(
    fit_layered_network(object, min_detected = -1L, verbose = FALSE),
    "non-negative"
  )
  expect_error(
    fit_layered_network(object, response_chunk = 0L, verbose = FALSE),
    "positive integer"
  )
  expect_error(
    fit_layered_network(object, target_chunk = NA_integer_, verbose = FALSE),
    "target_chunk.*positive integer"
  )
})

test_that("native sparse gathering preserves explicit values and logical fallback", {
  for (logical_values in c(FALSE, TRUE)) {
    x <- Matrix::sparseMatrix(
      i = c(1L, 3L, 2L, 4L), j = c(1L, 1L, 3L, 5L),
      x = if (logical_values) c(TRUE, FALSE, NA, TRUE) else c(0, -0, NaN, -Inf),
      dims = c(6L, 7L)
    )
    rows <- multiCSN:::row_compressed_matrix(x)
    for (selected in list(integer(), 2L, c(5L, 1L, 3L, 5L, 2L))) {
      for (by_row in c(FALSE, TRUE)) {
        source <- if (by_row) rows else x
        indices <- if (by_row) source@j else source@i
        lengths <- source@p[selected + 1L] - source@p[selected]
        positions <- rep.int(source@p[selected] + 1L, lengths) + sequence(lengths) - 1L
        expected <- methods::new("dgCMatrix",
          i = as.integer(indices[positions]),
          p = as.integer(c(0L, cumsum(lengths))), x = as.numeric(source@x[positions]),
          Dim = as.integer(c(if (by_row) ncol(source) else nrow(source), length(selected)))
        )
        actual <- if (by_row) multiCSN:::row_compressed_block(source, selected) else multiCSN:::column_compressed_block(source, selected)
        expect_identical(actual, expected)
        expect_identical(serialize(actual, NULL), serialize(expected, NULL))
      }
    }
  }
})

test_that("layered binding preserves all internal table schemas and empty parts", {
  schemas <- c("empty_tf_gene_edges", "empty_region_gene_edges", "empty_shared_accounting", "empty_region_gene_accounting")
  for (schema in schemas) {
    empty <- get(schema, envir = asNamespace("multiCSN"))
    part <- empty()[rep(NA_integer_, 3L), , drop = FALSE]
    for (name in names(part)) {
      part[[name]] <- switch(typeof(empty()[[name]]),
        character = c("x", "y", NA_character_),
        integer = c(2L, 0L, NA_integer_),
        double = c(1.5, Inf, NA_real_)
      )
    }
    parts <- list(first = part, NULL, empty(), last = part[3:1, , drop = FALSE])
    expected <- do.call(rbind, parts[c(1L, 4L)])
    rownames(expected) <- NULL
    expect_identical(multiCSN:::layered_rbind(parts, empty), expected)
    expect_identical(multiCSN:::layered_rbind(list(NULL, empty()), empty), empty())
  }
})

test_that("sparse gathering does not touch the random state and checks overflow", {
  set.seed(10143)
  x <- methods::as(Matrix::rsparsematrix(8, 9, density = 0.2), "dgCMatrix")
  seed <- .Random.seed
  multiCSN:::column_compressed_block(x, c(4L, 2L))
  expect_identical(.Random.seed, seed)
  expect_error(
    multiCSN:::gather_sparse_block(c(0L, .Machine$integer.max), integer(), numeric(), c(1L, 1L), 1L),
    "compressed matrix index limit"
  )
})

test_that("native candidate unions preserve sparse value and ordering semantics", {
  row_compressed_positive <-
    function(row_compressed, row) {
      p <- row_compressed@p
      positions <- p[[row]] + seq_len(p[[row + 1L]] - p[[row]])
      keep <- row_compressed@x[positions] > 0
      as.integer(row_compressed@j[positions][keep]) + 1L
    }
  set.seed(10414)
  for (replicate in seq_len(30)) {
    gate <- Matrix::rsparsematrix(17, 9, density = 0.35)
    domain <- Matrix::rsparsematrix(11, 19, density = 0.3)
    gate@x[seq_len(4)] <- c(0, NA_real_, NaN, Inf)
    domain@x[seq_len(4)] <- c(0, NA_real_, NaN, -Inf)
    g <- multiCSN:::row_compressed_matrix(gate)
    d <- multiCSN:::row_compressed_matrix(domain)
    mapping <- sample(c(seq_len(17), NA_integer_), 19, replace = TRUE)
    for (row in seq_len(nrow(domain))) {
      links <- mapping[multiCSN:::row_compressed_nonzero(d, row)]
      links <- links[!is.na(links)]
      expected <- sort(unique(unlist(lapply(links, function(i) {
        row_compressed_positive(g, i)
      }), use.names = FALSE)))
      actual <- multiCSN:::sparse_candidate_union(
        d@p, d@j, d@x, g@p, g@j, g@x, mapping, row, ncol(g)
      )
      expect_identical(actual, as.integer(expected))
    }
    for (row in seq_len(nrow(gate))) {
      expect_identical(
        multiCSN:::positive_sparse_row(g@p, g@j, g@x, row),
        row_compressed_positive(g, row)
      )
    }
  }
})

test_that("checkpoint identity detects package code changes without version changes", {
  package_root <- tempfile("package-identity-")
  dir.create(package_root)
  on.exit(unlink(package_root, recursive = TRUE), add = TRUE)
  dir.create(file.path(package_root, "R"))
  dir.create(file.path(package_root, "libs"))
  writeLines("Version: 0.0.1", file.path(package_root, "DESCRIPTION"))
  writeLines("useDynLib(multiCSN)", file.path(package_root, "NAMESPACE"))
  writeBin(as.raw(1:5), file.path(package_root, "R", "multiCSN.rdb"))
  binary <- file.path(package_root, "libs", "multiCSN.so")
  writeBin(as.raw(1:9), binary)
  original <- multiCSN:::layered_package_identity
  testthat::local_mocked_bindings(
    layered_package_identity = function(root) original(package_root),
    .package = "multiCSN"
  )
  root <- tempfile("checkpoint-identity-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  inputs <- list(cells = "c", features = "g", regulators = "r")
  settings <- list(cores = 1L)
  multiCSN:::layered_prepare_checkpoint(root, inputs, settings)
  expect_silent(multiCSN:::layered_prepare_checkpoint(root, inputs, settings))
  writeBin(as.raw(9:1), binary)
  expect_error(multiCSN:::layered_prepare_checkpoint(root, inputs, settings), "different input")
  writeBin(as.raw(1:9), binary)
  expect_silent(multiCSN:::layered_prepare_checkpoint(root, inputs, settings))
  added <- file.path(package_root, "libs", "extra.dll")
  writeBin(as.raw(1:3), added)
  expect_error(multiCSN:::layered_prepare_checkpoint(root, inputs, settings), "different input")
  unlink(added)
  unlink(binary)
  expect_error(multiCSN:::layered_prepare_checkpoint(root, inputs, settings), "different input")
  writeBin(as.raw(1:9), binary)
  writeBin(as.raw(5:1), file.path(package_root, "R", "multiCSN.rdb"))
  expect_error(multiCSN:::layered_prepare_checkpoint(root, inputs, settings), "different input")
})

test_that("shared Gram reuse preserves complete shared-design output", {
  set.seed(91014)
  design <- scale(matrix(rnorm(300), 30, 10))
  colnames(design) <- paste0("r", seq_len(10))
  response <- Matrix::Matrix(matrix(rnorm(240), 30, 8), sparse = TRUE)
  settings <- list(
    response_chunk = 3L, cores = 1L, exclude_response_alias = TRUE,
    max_support_size = NULL, min_improvement = 1e-10
  )
  candidates <- function(i) seq_len(ncol(design))
  reference <- multiCSN:::fit_layered_shared_design(
    design, response, paste0("g", seq_len(8)), candidates, "TF-gene", settings
  )
  reused <- multiCSN:::fit_layered_shared_design(
    design, response, paste0("g", seq_len(8)), candidates, "TF-gene", settings,
    gram = crossprod(design)
  )
  expect_identical(reused, reference)
})
