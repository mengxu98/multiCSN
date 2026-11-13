brute_average_precision <- function(labels, scores) {
  thresholds <- sort(unique(scores), decreasing = TRUE)
  recall_previous <- 0
  total <- 0
  for (threshold in thresholds) {
    selected <- scores >= threshold
    recall <- sum(labels & selected) / sum(labels)
    total <- total + (recall - recall_previous) * sum(labels & selected) / sum(selected)
    recall_previous <- recall
  }
  total
}

brute_expected_precision <- function(labels, scores, k, draws = 4000L) {
  set.seed(1)
  mean(replicate(draws, {
    ordering <- order(-scores, stats::runif(length(scores)))
    mean(labels[ordering][seq_len(k)])
  }))
}

test_that("tie-aware ranking metrics match reference definitions", {
  set.seed(5)
  for (iteration in 1:20) {
    n <- sample(10:60, 1)
    scores <- round(stats::runif(n), 1)
    labels <- stats::runif(n) < 0.3
    labels[1:2] <- c(TRUE, FALSE)
    expect_equal(tie_aware_average_precision(labels, scores), brute_average_precision(labels, scores))
    pairs <- outer(scores[labels], scores[!labels], ">") + 0.5 * outer(scores[labels], scores[!labels], "==")
    expect_equal(tie_aware_auroc(labels, scores), mean(pairs))
    k <- sample(seq_len(n), 1)
    expect_equal(
      expected_precision_at_k(labels, scores, k),
      brute_expected_precision(labels, scores, k),
      tolerance = 0.05
    )
    combined <- ranking_metrics(scores, labels, top_k = k)
    expect_equal(combined$auroc, tie_aware_auroc(labels, scores))
    expect_equal(combined$average_precision, tie_aware_average_precision(labels, scores))
    expect_equal(combined$topk$precision, expected_precision_at_k(labels, scores, k))
  }
})

test_that("ranking metrics are invariant to the order of tied scores", {
  scores <- c(3, 2, 2, 2, 1, 0)
  labels <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  swapped <- c(1, 3, 2, 4, 5, 6)
  expect_identical(
    ranking_metrics(scores, labels, top_k = 2:4)[c("auroc", "average_precision", "auprc", "topk")],
    ranking_metrics(scores[swapped], labels[swapped], top_k = 2:4)[c("auroc", "average_precision", "auprc", "topk")]
  )
  expect_true(is.na(tie_aware_auroc(c(TRUE, TRUE), c(1, 2))))
  expect_true(is.na(ranking_metrics(1:3, c(FALSE, FALSE, FALSE))$auroc))
  expect_identical(tie_aware_auroc(c(TRUE, NA, FALSE), c(1, 5, NaN)), NA_real_)
})

test_that("signed ordinal weights rank magnitudes with shared ties", {
  expect_equal(signed_ordinal_weight(c(-2, 1, 0, 2, NA)), c(-0.75, 0.25, 0, 0.75, 0))
  expect_identical(signed_ordinal_weight(c(0, NA)), c(0, 0))
})

test_that("edge adapters accept common layouts", {
  expect_identical(
    as_edge_table(data.frame(Gene1 = "a", Gene2 = "b", EdgeWeight = "0.5")),
    data.frame(regulator = "a", target = "b", weight = 0.5)
  )
  expect_identical(
    as_edge_table(data.frame(Gene2 = "b", Gene1 = "a"), weight = FALSE)$regulator,
    "a"
  )
  expect_identical(as_edge_table(data.frame(x = "a", y = "b"))$weight, 1)
  m <- matrix(c(0, 0.2, -0.7, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_identical(
    weight_matrix_to_edges(m),
    data.frame(regulator = c("B", "A"), target = c("A", "B"), weight = c(-0.7, 0.2))
  )
})

test_that("evaluate_network scores the complete gene-pair universe", {
  truth <- data.frame(regulator = c("A", "A", "B"), target = c("B", "C", "C"))
  perfect <- data.frame(regulator = c("A", "A", "B"), target = c("B", "C", "C"), weight = c(3, 2, 1))
  scores <- evaluate_network(perfect, truth)
  expect_identical(scores$AUROC, 1)
  expect_identical(scores$AP, 1)
  expect_identical(scores$EarlyPrecision, 1)
  expect_identical(scores$GoldEdges, 3L)
  expect_equal(scores$EPR, 2)
  reversed <- data.frame(regulator = c("C", "C", "B"), target = c("A", "B", "A"), weight = 1)
  expect_identical(evaluate_network(reversed, truth)$AUROC, 0)
  empty <- evaluate_network(perfect, data.frame(regulator = "A", target = "A"))
  expect_true(is.na(empty$AUROC))
})

test_that("network evaluation rejects missing and empty truth endpoints", {
  prediction <- data.frame(regulator = "A", target = "B", weight = 1)
  for (bad in c(NA_character_, "")) {
    expect_error(evaluate_network(prediction, data.frame(regulator = bad, target = "B")),
                 "Ground-truth edges require")
    expect_error(evaluate_network(prediction, data.frame(regulator = "A", target = bad)),
                 "Ground-truth edges require")
  }
})

test_that("ranking_metrics rejects invalid scores and labels before coercion", {
  for (scores in list(
    numeric(), c(1, NA_real_), c(1, Inf), c("1", "2"),
    c(TRUE, FALSE), matrix(1:2), c(1 + 1i, 2)
  )) {
    labels <- rep(c(FALSE, TRUE), length.out = length(scores))
    expect_error(ranking_metrics(scores, labels), "scores.*finite numeric vector")
  }
  for (labels in list(
    c(0, 2), c(0, 0.5), c(0, NA_real_), c(0, Inf),
    c("0", "1"), matrix(0:1), c(0, 1 + 1i)
  )) {
    expect_error(ranking_metrics(1:2, labels), "labels.*logical.*0 or 1")
  }
  expect_error(ranking_metrics(1:2, TRUE), "same length")
})

test_that("ranking_metrics rejects invalid budgets even when a class is absent", {
  for (labels in list(c(FALSE, TRUE), c(TRUE, TRUE))) {
    for (top_k in list(
      integer(), 0, -1, 1.5, NA_integer_, NaN, Inf,
      "1", TRUE, matrix(1L), c(1L, 0L), .Machine$integer.max + 1
    )) {
      expect_error(ranking_metrics(1:2, labels, top_k), "top_k.*positive integers")
    }
  }
  expect_identical(ranking_metrics(1:2, 0:1, c(1, 3, 1))$topk$top_k, c(1L, 2L, 1L))
  expect_identical(ranking_metrics(1:2, c(TRUE, TRUE), 3)$topk$top_k, 3L)
})
