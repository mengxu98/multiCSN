test_that("module coverage counts unique reachable genes within each module", {
  x <- plot_reachability_community_coverage(
    list(A = c("g1", "g1", "g2"), B = "g3", Empty = character()),
    list(M1 = c("g1", "g2"), M2 = c("g2", "g3"))
  )
  expect_equal(x$A, c(1, .5))
  expect_equal(x$B, c(0, .5))
  expect_equal(x$Empty, c(0, 0))
  expect_identical(rownames(x), c("M1", "M2"))
})

test_that("module paths use weighted directed routes and mark unreachable nodes", {
  edges <- data.frame(
    from = c("A", "B", "A", "D"),
    to = c("B", "C", "C", "E"), edge_length = c(1, 1, 5, 1)
  )
  x <- suppressWarnings(find_paths_to(edges, "C"))
  expect_equal(x["A", "path_length"], 2)
  expect_identical(x["A", "path"], "A--B--C")
  expect_equal(x["C", "path_length"], 0)
  expect_true(is.infinite(x["D", "path_length"]))
  expect_true(is.na(x["D", "path"]))
})
