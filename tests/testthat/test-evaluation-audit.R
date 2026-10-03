test_that("asymmetric confusion counts retain their row and column meaning", {
  pred_data <- list(predictor_binary = c(1, 0, 0), true_label = c(1, 1, 0))
  expect_equal(binary_counts_from_pred_data(pred_data), c(TN = 1, FP = 0, FN = 1, TP = 1))
  truth <- data.frame(regulator = c("a", "b"), target = c("b", "c"))
  pred <- data.frame(regulator = c("a", "a"), target = c("b", "c"), weight = c(3, 2))
  expect_equal(calculate_precision(pred, truth)$metrics$Value, 1)
  expect_equal(calculate_recall(pred, truth)$metrics$Value, 0.5)
})

test_that("EPR uses a fixed candidate pool and actual tie-expanded predictions", {
  truth <- data.frame(regulator = c("a", "b", "c"), target = c("b", "c", "a"))
  short <- data.frame(regulator = c("x", "a"), target = c("y", "b"), weight = c(100, -5))
  expect_equal(calculate_epr(short, truth)$metrics$Value, 2)
  tied <- data.frame(
    regulator = c("a", "b", "a", "c"),
    target = c("b", "c", "c", "b"), weight = c(5, 4, 4, 4)
  )
  expect_equal(calculate_epr(tied, truth)$metrics$Value, 1)
  expect_equal(calculate_epr(tied[c(4, 3, 2, 1), ], truth)$metrics$Value, 1)
  short$weight <- 0
  expect_equal(calculate_epr(short, truth)$metrics$Value, 0)
  expect_equal(calculate_epr(short[0, ], truth)$metrics$Value, 0)
})

test_that("edge endpoint punctuation does not merge distinct relations", {
  truth <- data.frame(regulator = c("a-b", "a"), target = c("c", "b-c"))
  pred <- transform(truth, weight = c(2, 1))
  expect_equal(calculate_si(pred, truth)$metrics$Value, 2)
  expect_equal(calculate_ji(pred[1, ], truth)$metrics$Value, 0.5)
  pred <- data.frame(regulator = "a|||b", target = "c", weight = 2)
  truth <- data.frame(regulator = c("a|||b", "a"), target = c("c", "b|||c"))
  gold <- prepare_calculate_metrics(pred, truth)
  expect_equal(sum(gold$label), 2)
  expect_equal(sum(gold$weight > 0), 1)
})

test_that("feedback motif counts use the actual cycle list", {
  skip_if_not_installed("igraph")
  truth <- data.frame(regulator = c("a", "b", "c"), target = c("b", "c", "a"))
  pred <- transform(truth[1:2, ], weight = 1)
  metric <- calculate_motif_ratios(pred, truth)$metrics
  expect_equal(metric$Value[metric$Metric == "FBLRatio"], 0)
})

test_that("Youden ties select the highest threshold deterministically", {
  scores <- 6:1
  labels <- c(1, 0, 1, 0, 1, 0)
  expect_equal(select_best_binary_threshold(scores, labels)$threshold, 6)
})

test_that("signed EPR distinguishes missing truth from empty predictions", {
  truth <- data.frame(
    regulator = c("a", "b", "c"), target = c("b", "c", "a"),
    type = c("+", "-", "-")
  )
  pred <- data.frame(regulator = character(), target = character(), weight = numeric())
  expect_equal(calculate_signed_epr(pred, truth)$metrics$Value, c(0, 0))
  truth$type <- "+"
  expect_true(is.na(calculate_signed_epr(pred, truth)$metrics$Value[2]))
})
