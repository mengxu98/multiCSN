test_that("roundoff-scale evidence ties are grouped without chaining", {
  expect_equal(network_ordinal_weight(c(1, -1, 1), c(354, 354 + 1e-11, 10)),
    c(.75, -.75, .25),
    tolerance = 0
  )

  e <- c(1, 1 - 1.5e-12, 1 - 3e-12)
  expect_equal(network_ordinal_weight(rep(1, 3), e), c(.75, .75, .25), tolerance = 0)
  expect_equal(network_ordinal_weight(c(1, 1), c(1, 1 + 1e-8)), c(.25, .75))
})

test_that("numerically tied strongest paths cannot invent a confident sign", {
  chains <- data.frame(
    regulator = "TF", target = "G", region = c("a", "b"),
    chain_delta_bic = c(5, 5 + 1e-13), weight = c(.5, -.5)
  )
  expect_equal(nrow(strongest_mediated_projection(chains)), 0L)
  chains$weight <- .5
  projection <- strongest_mediated_projection(chains)
  expect_equal(projection$region, "a")
  expect_equal(projection$strongest_path_ties, 2L)

  expect_equal(projection$chain_delta_bic, max(chains$chain_delta_bic), tolerance = 0)
})
