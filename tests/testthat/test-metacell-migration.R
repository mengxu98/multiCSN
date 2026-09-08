test_that("weighted metacells match independent cell-group weighted means", {
  skip_if_not_installed("BiocNeighbors")
  set.seed(100)
  x <- matrix(rexp(40 * 6), 40, 6,
    dimnames = list(paste0("c", 1:40), paste0("g", 1:6))
  )
  groups <- NULL
  aggregate_original <- meta_cell_ge
  local_mocked_bindings(meta_cell_ge = function(ge, groups, ...) {
    assign("groups", groups, envir = parent.env(environment()))
    aggregate_original(ge, groups, ...)
  })
  weights <- seq_len(nrow(x))
  observed <- meta_cells(x,
    gamma = 4, genes_use = colnames(x), pc_num = 2,
    weights = weights
  )
  expect_length(groups, nrow(x))
  expected <- t(vapply(sort(unique(groups)), function(g) {
    cells <- which(groups == g)
    colSums(x[cells, , drop = FALSE] * weights[cells]) / sum(weights[cells])
  }, numeric(ncol(x))))
  expect_equal(unname(as.matrix(observed)), unname(expected), tolerance = 1e-12)
  unweighted <- meta_cells(x, gamma = 4, genes_use = colnames(x), pc_num = 2)
  unit_weighted <- meta_cells(x,
    gamma = 4, genes_use = colnames(x), pc_num = 2,
    weights = rep(1, nrow(x))
  )
  expect_equal(unit_weighted, unweighted, tolerance = 0)
})

test_that("subsampling strategies keep cell-by-gene shape and neighborhood means", {
  skip_if_not_installed("BiocNeighbors")
  set.seed(71)
  x <- matrix(rexp(40 * 6), 40, 6,
    dimnames = list(paste0("c", 1:40), paste0("g", 1:6))
  )
  set.seed(19)
  seeds <- sample(seq_len(nrow(x)), 20)
  expect_equal(subsampling(x, method = "sample", ratio = .5, seed = 19, verbose = FALSE), x[seeds, ])
  d <- as.matrix(dist(x))
  diag(d) <- Inf
  expected <- t(vapply(seeds, function(i) {
    neighbors <- order(d[i, ])[1:3]
    colMeans(x[c(i, neighbors), , drop = FALSE])
  }, numeric(ncol(x))))
  observed <- subsampling(x, method = "pseudobulk", ratio = .5, seed = 19, k = 3, verbose = FALSE)
  expect_equal(unname(observed), unname(expected), tolerance = 1e-12)
  m <- subsampling(x,
    method = "meta_cells", ratio = .25, genes_use = colnames(x), pc_num = 2,
    verbose = FALSE
  )
  expect_equal(dim(m), c(10L, 6L))
  expect_identical(colnames(m), colnames(x))
  sparse <- subsampling(Matrix::Matrix(x, sparse = TRUE), "sample", .5,
    seed = 19, verbose = FALSE
  )
  expect_s4_class(sparse, "sparseMatrix")
  expect_equal(as.matrix(sparse), x[seeds, ])
})

test_that("aggregate_assay still summarizes predefined metadata groups", {
  counts <- matrix(seq_len(24), 4, 6,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:6))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$group <- c("A", "A", "B", "B", "B", "A")
  result <- aggregate_assay(object, group_name = "group", assay = "RNA", layer = "counts")
  expected <- rbind(A = rowMeans(counts[, c(1, 2, 6)]), B = rowMeans(counts[, 3:5]))
  expect_equal(as.matrix(result@assays$RNA@misc$summary$group), expected)
  expect_equal(
    SeuratObject::LayerData(result, assay = "RNA", layer = "counts"),
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  )
})
