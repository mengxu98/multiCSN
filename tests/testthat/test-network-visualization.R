test_that("heatmap conversion preserves filtered values and numeric gene order", {
  weights <- matrix(seq_len(9), 3, dimnames = list(
    c("g10", "g2", "g1"), c("g10", "g2", "g1")
  ))
  weights["g2", "g1"] <- -4
  converted <- network_to_heatmap_matrix(
    weights,
    regulators = c("g10", "g2"), targets = c("g10", "g1"),
    switch_matrix = FALSE
  )
  expect_identical(rownames(converted), c("g2", "g10"))
  expect_identical(colnames(converted), c("g1", "g10"))
  expect_equal(converted, weights[c("g2", "g10"), c("g1", "g10")])
})

test_that("network ranking excludes zero and missing weights before node filtering", {
  edges <- data.frame(
    regulator = c("g10", "g2", "g2", "g10"),
    target = c("g1", "g1", "zero", "missing"),
    weight = c(-2, 3, 0, NA_real_)
  )
  ranked <- compute_betweenness_degree(
    list(state = edges),
    regulators = "g10", targets = c("g1", "zero", "missing")
  )$state
  expect_setequal(ranked$gene, c("g10", "g1"))
  expect_identical(ranked$gene[ranked$is_regulator], "g10")
})
