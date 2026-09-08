test_that("sparse formula fitting preserves predictor names and no-intercept columns", {
  set.seed(19)
  d <- data.frame(x = rnorm(60), z = rnorm(60))
  d$y <- 2 * d$x + rnorm(60, sd = 0.02)
  expect_identical(fit_srm2(y ~ x, d, verbose = FALSE)$coefficients$variable, "x")
  expect_identical(fit_srm2(y ~ 0 + x, d, verbose = FALSE)$coefficients$variable, "x")
  expect_error(fit_srm2(y ~ x, d, n_folds = 3), "unused argument")
})

test_that("Seurat inference forwards support limits to both branches", {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1:12, 3, 4,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:4))
  ), sparse = TRUE))
  static_limit <- dynamic_limit <- NULL
  local_mocked_bindings(
    .csn_celltypes = function(...) "A",
    get_attribute = function(...) "g1",
    fit_models = function(object, ...) {
      static_limit <<- list(...)$max_support_size
      object
    },
    .process_csn = function(object, ...) object,
    export_csn = function(...) NULL,
    state_dynamic = function(object, ...) {
      dynamic_limit <<- list(...)$max_support_size
      object
    },
    .package = "multiCSN"
  )
  inferCSN(object, max_support_size = 1L, verbose = FALSE)
  expect_identical(static_limit, 1L)
  inferCSN(object, pseudotime = "time", max_support_size = 2L, verbose = FALSE)
  expect_identical(dynamic_limit, 2L)
})


test_that("support limits are only accepted by the sparse backend", {
  set.seed(4)
  d <- data.frame(y = rnorm(30), x = rnorm(30))
  expect_error(fit_model(y ~ x, d, method = "glm", max_support_size = 1), "only supported")
  expect_true(is.list(fit_model(y ~ x, d, method = "glm")))
})

test_that("greedy_l0 dispatch matches the sparse formula wrapper", {
  set.seed(19)
  d <- data.frame(x = rnorm(60), z = rnorm(60))
  d$y <- 2 * d$x + rnorm(60, sd = 0.02)
  expect_equal(
    fit_model(y ~ x + z, d, method = "greedy_l0", verbose = FALSE),
    fit_srm2(y ~ x + z, d, verbose = FALSE)
  )
  expect_error(fit_model(y ~ x, d, method = "srm"), "arg")
})
