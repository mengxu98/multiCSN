test_that("dynamic filtering uses inferCSN and preserves reusable object output", {
  set.seed(42)
  time <- seq(0, 1, length.out = 60)
  x <- matrix(rpois(8 * 60, 4), 8, dimnames = list(paste0("g", 1:8), paste0("c", 1:60)))
  x[1, ] <- rpois(60, exp(time * 3))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(x, sparse = TRUE))
  object$time <- time
  SeuratObject::LayerData(object, layer = "data") <- log1p(Matrix::Matrix(x, sparse = TRUE))
  config <- list(fit_method = "pretsa", padjust_threshold = 1, n_candidates = 2)
  out <- get_trend_genes(object, "time", config, rownames(object), 1, verbose = FALSE)
  reference <- inferCSN::select_trend_features(t(log1p(x)), time,
    padjust_threshold = 1, n_candidates = 2
  )
  expect_identical(out$genes, reference$features)
  expect_identical(get_trend_genes(out$seurat, "time", config, rownames(object), 1, verbose = FALSE), out$genes)
  expect_equal(colnames(out$seurat@tools$DynamicFeatures_time$raw_matrix), c("pseudotime", rownames(object)))
  expect_identical(get_trend_genes(object, "time", NULL, rownames(object), 1), rownames(object))
})
