test_that("KNN agrees with direct distances including ties and small samples", {
  skip_if_not_installed("BiocNeighbors")
  fixtures <- list(
    matrix(c(0, 0, 1, 0, -1, 0, 0, 2, 0, 0), ncol = 2, byrow = TRUE),
    matrix(c(0, 1), ncol = 1),
    matrix(0, nrow = 1, ncol = 2)
  )
  for (x in fixtures) {
    for (k in unique(c(1L, max(1L, nrow(x) - 1L)))) {
      observed <- suppressWarnings(multiCSN:::build_nn2(x, k = k))
      distances <- unname(as.matrix(dist(x)))
      diag(distances) <- Inf
      for (i in seq_len(nrow(x))) {
        valid <- which(!is.na(observed$idx[i, ]))
        expected_n <- min(k, nrow(x) - 1L)
        expect_length(valid, expected_n)
        if (!expected_n) next
        idx <- observed$idx[i, valid]
        expect_false(i %in% idx)
        expect_length(unique(idx), expected_n)
        expect_equal(observed$dist[i, valid], distances[i, idx], tolerance = 1e-12)
        expect_equal(sort(observed$dist[i, valid]), sort(distances[i, ])[seq_len(expected_n)], tolerance = 1e-12)
      }
    }
  }
})

test_that("KNN has exact neighbor identities when distances are not tied", {
  skip_if_not_installed("BiocNeighbors")
  x <- cbind(c(0, 0.3, 1.7, 4.2, 9.1), c(0, 1.1, 0.2, 2.4, 3))
  d <- as.matrix(dist(x))
  diag(d) <- Inf
  expected <- t(apply(d, 1, order))[, 1:2, drop = FALSE]
  one <- multiCSN:::build_nn2(x, k = 2, mode = "out")
  expect_equal(one$idx, unname(expected))
  expect_equal(igraph::ecount(one$graph_knn), 2L * nrow(x))
  expect_false(any(igraph::which_loop(one$graph_knn)))
  two <- thisutils::run_biocneighbors_knn(x, k = 2, exclude_self = TRUE, n_threads = 2)
  expect_identical(one$idx, two$idx)
  expect_identical(one$dist, two$dist)
})

test_that("metacell aggregation agrees with explicit weighted group means", {
  x <- matrix(seq_len(24), nrow = 4)
  groups <- c(2L, 1L, 2L, 3L, 1L, 3L)
  weights <- c(1, 2, 3, 4, 2, 1)
  for (weighted in c(FALSE, TRUE)) {
    w <- if (weighted) weights else rep(1, ncol(x))
    expected <- sapply(sort(unique(groups)), function(g) {
      keep <- which(groups == g)
      as.vector(x[, keep, drop = FALSE] %*% w[keep] / sum(w[keep]))
    })
    observed <- multiCSN:::meta_cell_ge(x, groups, weights = if (weighted) weights else NULL)
    expect_equal(as.matrix(observed), expected)
  }
})
