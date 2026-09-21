# Reference implementation mirroring the documented projection semantics:
# per regulator-target pair keep the strongest evidence, treat values within
# 1e-12 * (1 + |max|) as ties, drop pairs whose ties disagree in sign and use the
# lexicographically first region as representative.
reference_strongest_projection <- function(chains) {
  if (!nrow(chains)) {
    return(data.frame(
      regulator = character(), target = character(), region = character(),
      chain_delta_bic = numeric(), strongest_path_ties = integer(),
      strongest_path_sign_consistent = logical(), weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  groups <- split(
    seq_len(nrow(chains)),
    interaction(chains$regulator, chains$target, drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(index) {
    group <- chains[index, , drop = FALSE]
    best <- max(group$chain_delta_bic)
    tied <- group[abs(group$chain_delta_bic - best) <= 1e-12 * (1 + abs(best)), , drop = FALSE]
    directions <- unique(sign(tied$weight))
    chosen <- tied[order(tied$region), , drop = FALSE][1L, , drop = FALSE]
    data.frame(
      regulator = chosen$regulator, target = chosen$target, region = chosen$region,
      chain_delta_bic = best,
      chain_direction = if (length(directions) == 1L) directions else NA_real_,
      strongest_path_ties = nrow(tied),
      strongest_path_sign_consistent = length(directions) == 1L,
      stringsAsFactors = FALSE
    )
  })
  projection <- do.call(rbind, rows)
  projection <- projection[projection$strongest_path_sign_consistent, , drop = FALSE]
  projection$weight <- multiCSN:::.signed_ordinal_from_evidence(
    projection$chain_direction, projection$chain_delta_bic
  )
  projection$chain_direction <- NULL
  projection <- projection[order(
    -abs(projection$weight), projection$regulator, projection$target, projection$region
  ), c(
    "regulator", "target", "region", "chain_delta_bic", "strongest_path_ties",
    "strongest_path_sign_consistent", "weight"
  ), drop = FALSE]
  rownames(projection) <- NULL
  projection
}

make_chains <- function(regulator, region, target, evidence, weight) {
  data.frame(
    regulator = regulator, region = region, target = target,
    chain_delta_bic = evidence, weight = weight, stringsAsFactors = FALSE
  )
}

test_that("projection matches the reference implementation on random chains", {
  set.seed(11)
  pairs <- expand.grid(
    regulator = paste0("tf", 1:6), target = paste0("g", 1:25),
    stringsAsFactors = FALSE
  )
  rows <- do.call(rbind, lapply(seq_len(300L), function(i) {
    pair <- pairs[sample.int(nrow(pairs), 1L), , drop = FALSE]
    make_chains(
      pair$regulator, sprintf("chr1:%d-%d", i * 100L, i * 100L + 50L), pair$target,
      stats::runif(1, 1, 10), sample(c(-1, 1), 1L) * stats::runif(1, 0.1, 1)
    )
  }))
  chains <- rows[!duplicated(rows[, c("regulator", "region", "target")]), , drop = FALSE]
  expect_gt(sum(duplicated(chains[, c("regulator", "target")])), 0L)
  expect_identical(
    strongest_mediated_projection(chains),
    reference_strongest_projection(chains)
  )
})

test_that("projection handles ties, sign conflicts and tolerance boundaries", {
  chains <- make_chains(
    regulator = c("tf1", "tf1", "tf1", "tf2", "tf2", "tf2", "tf3", "tf3"),
    region = c(
      "chr1:300-350", "chr1:100-150", "chr1:200-250",
      "chr1:100-150", "chr1:200-250", "chr1:300-350",
      "chr1:100-150", "chr1:200-250"
    ),
    target = c("g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3"),
    evidence = c(5, 5, 5, 9, 9 - 1e-13, 4, 3, 3 + 1e-12 * (1 + 3)),
    weight = c(0.5, 0.5, 0.5, 0.4, -0.4, 0.2, 0.1, 0.1)
  )
  projection <- strongest_mediated_projection(chains)
  expect_identical(projection, reference_strongest_projection(chains))
  tf1 <- projection[projection$regulator == "tf1", , drop = FALSE]
  expect_identical(tf1$region, "chr1:100-150")
  expect_identical(tf1$strongest_path_ties, 3L)
  expect_false("tf2" %in% projection$regulator)
  tf3 <- projection[projection$regulator == "tf3", , drop = FALSE]
  expect_identical(tf3$strongest_path_ties, 2L)
  expect_equal(tf3$chain_delta_bic, 3 + 1e-12 * (1 + 3))
})

test_that("projection validates its input and handles empty tables", {
  empty <- make_chains(character(), character(), character(), numeric(), numeric())
  projection <- strongest_mediated_projection(empty)
  expect_identical(nrow(projection), 0L)
  expect_identical(
    names(projection),
    c(
      "regulator", "target", "region", "chain_delta_bic", "strongest_path_ties",
      "strongest_path_sign_consistent", "weight"
    )
  )
  duplicated <- make_chains(
    c("tf1", "tf1"), c("chr1:1-2", "chr1:1-2"), c("g1", "g1"), c(1, 2), c(1, 1)
  )
  expect_error(strongest_mediated_projection(duplicated), "Invalid factorized chain table")
  missing_column <- make_chains("tf1", "chr1:1-2", "g1", 1, 1)[, -5]
  expect_error(strongest_mediated_projection(missing_column), "Invalid factorized chain table")
})
