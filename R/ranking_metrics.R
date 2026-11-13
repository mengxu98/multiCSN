#' Tie-aware ranking metrics
#'
#' Ranking metrics for a score vector against binary labels. Scores tied at a
#' threshold are treated as one block, so every metric is invariant to the
#' order of tied predictions.
#'
#' \code{ranking_metrics()} returns the threshold average precision
#' (\code{average_precision}), the trapezoidal area under the precision-recall
#' curve (\code{auprc}), the AUROC with average ranks for ties and expected
#' precision, recall and F1 for each \code{top_k}, where a top-k boundary that
#' falls inside a tie block takes the block's expected positives.
#'
#' @param scores Numeric prediction scores; larger is more confident.
#' @param labels Logical or 0/1 truth labels.
#' @param top_k Edge budgets for the expected top-k precision.
#' @param k Edge budget for \code{expected_precision_at_k()}.
#' @return \code{ranking_metrics()} returns a list with \code{auroc},
#'   \code{average_precision}, \code{auprc}, a \code{topk} data frame and the
#'   descending score \code{order}. The other functions return one number;
#'   missing labels and non-finite scores are dropped and \code{NA} is returned
#'   when a class is absent.
#' @examples
#' scores <- c(0.9, 0.8, 0.8, 0.3, 0.1)
#' labels <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
#' ranking_metrics(scores, labels, top_k = 2)$topk
#' tie_aware_auroc(labels = labels, scores = scores)
#' tie_aware_average_precision(labels = labels, scores = scores)
#' expected_precision_at_k(labels = labels, scores = scores, k = 2)
#' @name ranking_metrics
NULL

#' @rdname ranking_metrics
#' @export
ranking_metrics <- function(scores, labels,
                            top_k = c(25L, 50L, 100L, 200L, 500L, 1000L)) {
  scores <- suppressWarnings(as.numeric(scores))
  labels <- as.logical(labels)
  if (length(scores) != length(labels)) {
    stop("scores and labels must have the same length.", call. = FALSE)
  }
  if (length(scores) == 0L || any(!is.finite(scores)) || anyNA(labels)) {
    stop("Scores must be finite and labels must be non-missing.", call. = FALSE)
  }
  n_positive <- sum(labels)
  n_negative <- length(labels) - n_positive
  if (n_positive == 0L || n_negative == 0L) {
    return(list(
      auroc = NA_real_,
      average_precision = NA_real_,
      auprc = NA_real_,
      topk = data.frame(
        top_k = as.integer(top_k),
        cutoff_score = NA_real_,
        precision = NA_real_,
        recall = NA_real_,
        f1 = NA_real_
      ),
      order = order(-scores, seq_along(scores))
    ))
  }
  ord <- order(-scores, seq_along(scores))
  ordered_scores <- scores[ord]
  ordered_labels <- labels[ord]
  block_end <- c(which(diff(ordered_scores) != 0), length(ordered_scores))
  true_selected <- cumsum(ordered_labels)[block_end]
  recall <- true_selected / n_positive
  precision <- true_selected / block_end
  average_precision <- sum(c(recall[1L], diff(recall)) * precision)
  pr_recall <- c(0, recall)
  pr_precision <- c(1, precision)
  auprc <- sum(
    diff(pr_recall) *
      (utils::head(pr_precision, -1L) + utils::tail(pr_precision, -1L)) / 2
  )
  ranks <- rank(scores, ties.method = "average")
  auroc <- (sum(ranks[labels]) - n_positive * (n_positive + 1) / 2) /
    (as.double(n_positive) * n_negative)
  topk <- lapply(as.integer(top_k), function(k) {
    k <- min(k, length(scores))
    cutoff_score <- ordered_scores[[k]]
    above <- ordered_scores > cutoff_score
    tied <- ordered_scores == cutoff_score
    expected_tp <- sum(ordered_labels[above]) +
      (k - sum(above)) * mean(ordered_labels[tied])
    prec <- expected_tp / k
    rec <- expected_tp / n_positive
    data.frame(
      top_k = k,
      cutoff_score = cutoff_score,
      precision = prec,
      recall = rec,
      f1 = if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
    )
  })
  list(
    auroc = auroc,
    average_precision = average_precision,
    auprc = auprc,
    topk = do.call(rbind, topk),
    order = ord
  )
}

ranking_inputs <- function(labels, scores) {
  labels <- as.logical(labels)
  scores <- suppressWarnings(as.numeric(scores))
  if (length(labels) != length(scores)) {
    stop("labels and scores must have the same length.", call. = FALSE)
  }
  keep <- !is.na(labels) & is.finite(scores)
  list(labels = labels[keep], scores = scores[keep])
}

#' @rdname ranking_metrics
#' @export
tie_aware_average_precision <- function(labels, scores) {
  input <- ranking_inputs(labels, scores)
  labels <- input$labels
  scores <- input$scores
  positives <- sum(labels)
  if (!length(labels) || positives == 0L) {
    return(NA_real_)
  }
  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  scores <- scores[ordering]
  block_end <- cumsum(rle(scores)$lengths)
  block_positives <- diff(c(0, cumsum(labels)[block_end]))
  precision <- cumsum(block_positives) / block_end
  sum((block_positives / positives) * precision)
}

#' @rdname ranking_metrics
#' @export
tie_aware_auroc <- function(labels, scores) {
  input <- ranking_inputs(labels, scores)
  labels <- input$labels
  positives <- sum(labels)
  negatives <- length(labels) - positives
  if (positives == 0L || negatives == 0L) {
    return(NA_real_)
  }
  ranks <- rank(input$scores, ties.method = "average")
  (sum(ranks[labels]) - positives * (positives + 1) / 2) /
    (as.double(positives) * negatives)
}

#' @rdname ranking_metrics
#' @export
expected_precision_at_k <- function(labels, scores, k) {
  input <- ranking_inputs(labels, scores)
  labels <- input$labels
  scores <- input$scores
  k <- min(as.integer(k), length(labels))
  if (!length(labels) || is.na(k) || k <= 0L) {
    return(NA_real_)
  }
  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  scores <- scores[ordering]
  block_end <- cumsum(rle(scores)$lengths)
  block_size <- diff(c(0L, block_end))
  block_positives <- diff(c(0, cumsum(labels)[block_end]))
  remaining <- k
  expected_true <- 0
  for (index in seq_along(block_size)) {
    if (remaining <= 0L) break
    selected <- min(remaining, block_size[[index]])
    expected_true <- expected_true +
      selected * block_positives[[index]] / block_size[[index]]
    remaining <- remaining - selected
  }
  expected_true / k
}

#' Signed ordinal edge weights
#'
#' Converts signed coefficients into signed ranks in (0, 1]: magnitudes are
#' ranked from strongest to weakest, magnitudes within \code{tolerance} share
#' one rank, and the sign of the coefficient is kept. Zero and non-finite
#' coefficients receive weight 0. Use \code{\link{network_ordinal_weight}} when
#' the rank should follow deletion evidence rather than the coefficient size.
#'
#' @param values Signed coefficients.
#' @param tolerance Magnitudes closer than this are tied.
#' @return Numeric vector of signed ordinal weights.
#' @examples
#' signed_ordinal_weight(c(-2, 1, 0, 2, 0.5))
#' @export
signed_ordinal_weight <- function(values, tolerance = 1e-12) {
  values <- suppressWarnings(as.numeric(values))
  out <- numeric(length(values))
  keep <- is.finite(values) & values != 0
  if (!any(keep)) {
    return(out)
  }
  kept_index <- which(keep)
  magnitude <- abs(values[kept_index])
  ord <- order(-magnitude, kept_index)
  sorted_magnitude <- magnitude[ord]
  group <- integer(length(ord))
  group[[1]] <- 1L
  anchor <- sorted_magnitude[[1]]
  if (length(ord) > 1L) {
    for (i in 2:length(ord)) {
      if (abs(sorted_magnitude[[i]] - anchor) >= tolerance) {
        group[[i]] <- group[[i - 1L]] + 1L
        anchor <- sorted_magnitude[[i]]
      } else {
        group[[i]] <- group[[i - 1L]]
      }
    }
  }
  n_groups <- max(group)
  out[kept_index[ord]] <- sign(values[kept_index[ord]]) *
    (1 - (group - 0.5) / n_groups)
  out
}

#' Standardize an edge table
#'
#' Accepts the common layouts of predicted networks and ground truths
#' (\code{regulator}/\code{target}, \code{Gene1}/\code{Gene2}, or the first two
#' columns) and returns \code{regulator}, \code{target} and, for predictions,
#' a numeric \code{weight}. The weight column is the first of \code{weight},
#' \code{score}, \code{cor_value}, \code{importance}, \code{corr} and
#' \code{value}, then the third column, otherwise 1.
#'
#' @param table A data frame of edges.
#' @param weight Return a \code{weight} column.
#' @return A data frame.
#' @examples
#' as_edge_table(data.frame(Gene1 = "A", Gene2 = "B", EdgeWeight = 0.4))
#' @export
as_edge_table <- function(table, weight = TRUE) {
  if (is.null(table) || ncol(table) < 2L) {
    stop("An edge table needs at least two columns.", call. = FALSE)
  }
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  table_names <- tolower(colnames(table))
  if (all(c("regulator", "target") %in% table_names)) {
    columns <- c(match("regulator", table_names), match("target", table_names))
  } else if (all(c("gene1", "gene2") %in% table_names)) {
    columns <- c(match("gene1", table_names), match("gene2", table_names))
  } else if (!isTRUE(weight) && identical(table_names[1:2], c("gene2", "gene1"))) {
    columns <- c(2L, 1L)
  } else {
    columns <- c(1L, 2L)
  }
  if (!isTRUE(weight)) {
    out <- table[, columns, drop = FALSE]
    colnames(out) <- c("regulator", "target")
    return(out)
  }
  candidates <- c("weight", "score", "cor_value", "importance", "corr", "value")
  weight_column <- match(candidates, table_names)
  weight_column <- weight_column[!is.na(weight_column)][1L]
  if (is.na(weight_column)) {
    if (ncol(table) >= 3L) {
      weight_column <- 3L
    } else {
      table$weight <- 1
      weight_column <- ncol(table)
    }
  }
  out <- table[, c(columns, weight_column), drop = FALSE]
  colnames(out) <- c("regulator", "target", "weight")
  out$regulator <- as.character(out$regulator)
  out$target <- as.character(out$target)
  out$weight <- suppressWarnings(as.numeric(out$weight))
  out[stats::complete.cases(out[, c("regulator", "target")]), , drop = FALSE]
}

#' Convert a square weight matrix to an edge table
#'
#' @param weight_matrix Target-by-regulator weight matrix (rows are targets,
#'   columns are regulators).
#' @param remove_self Drop self-edges.
#' @return A data frame with \code{regulator}, \code{target} and \code{weight},
#'   sorted by decreasing absolute weight.
#' @examples
#' m <- matrix(c(0, 0.2, -0.7, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
#' weight_matrix_to_edges(m)
#' @export
weight_matrix_to_edges <- function(weight_matrix, remove_self = TRUE) {
  mat <- as.matrix(weight_matrix)
  if (is.null(rownames(mat))) {
    rownames(mat) <- seq_len(nrow(mat))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- seq_len(ncol(mat))
  }
  index <- which(!is.na(mat), arr.ind = TRUE)
  if (nrow(index) == 0) {
    return(data.frame(
      regulator = character(0), target = character(0), weight = numeric(0)
    ))
  }
  network_table <- data.frame(
    regulator = colnames(mat)[index[, 2]],
    target = rownames(mat)[index[, 1]],
    weight = as.numeric(mat[index]),
    stringsAsFactors = FALSE
  )
  if (remove_self) {
    network_table <- network_table[
      network_table$regulator != network_table$target, ,
      drop = FALSE
    ]
  }
  network_table <- network_table[
    order(abs(network_table$weight), decreasing = TRUE), ,
    drop = FALSE
  ]
  rownames(network_table) <- NULL
  network_table
}

#' BEELINE-style network evaluation
#'
#' Scores a predicted network against a ground truth on the complete directed
#' gene-pair universe of the ground-truth genes (self-pairs excluded).
#' Predictions are ranked by absolute weight and pairs without a prediction
#' share score 0. Ranking metrics are tie-aware. The binary summaries
#' (accuracy, precision, recall, F1, Jaccard index and signal index) use the
#' Youden-optimal score threshold, and early precision keeps every non-zero
#' prediction tied at the ground-truth edge count, dividing by the number
#' actually selected.
#'
#' @param network_table Predicted edges, see \code{\link{as_edge_table}}.
#' @param ground_truth Ground-truth edges with non-missing, non-empty endpoints.
#' @return A named list with \code{AUROC}, \code{AUPRC}, \code{AP},
#'   \code{ACC}, \code{Precision}, \code{Recall}, \code{F1}, \code{JI},
#'   \code{SI}, \code{EarlyPrecision}, \code{EPR}, \code{TopK},
#'   \code{EarlySelectedEdges}, \code{EarlyTruePositives},
#'   \code{PositiveRate} and \code{GoldEdges}.
#' @examples
#' truth <- data.frame(regulator = c("A", "A", "B"), target = c("B", "C", "C"))
#' prediction <- data.frame(
#'   regulator = c("A", "B", "C", "A"),
#'   target = c("B", "C", "A", "C"),
#'   weight = c(0.9, 0.5, 0.4, 0.1)
#' )
#' str(evaluate_network(prediction, truth))
#' @export
evaluate_network <- function(network_table, ground_truth) {
  gt <- as_edge_table(ground_truth, weight = FALSE)
  if (anyNA(gt$regulator) || anyNA(gt$target) ||
      any(!nzchar(gt$regulator) | !nzchar(gt$target))) {
    stop("Ground-truth edges require non-missing, non-empty endpoints.", call. = FALSE)
  }
  gt <- gt[gt$regulator != gt$target, , drop = FALSE]
  gt <- unique(gt)
  if (!nrow(gt)) {
    return(list(
      AUROC = NA_real_, AUPRC = NA_real_, AP = NA_real_, ACC = NA_real_,
      Precision = NA_real_, Recall = NA_real_, F1 = NA_real_, JI = NA_real_,
      SI = 0, EarlyPrecision = NA_real_, EPR = NA_real_, TopK = 0L,
      EarlySelectedEdges = 0L, EarlyTruePositives = 0L,
      PositiveRate = NA_real_, GoldEdges = 0L
    ))
  }
  pred <- as_edge_table(network_table)
  gt_genes <- sort(unique(c(gt$regulator, gt$target)))
  pred <- pred[
    pred$regulator %in% gt_genes & pred$target %in% gt_genes &
      pred$regulator != pred$target, ,
    drop = FALSE
  ]
  true_edge_ids <- unique(metric_edge_ids(gt$regulator, gt$target))
  pred <- pred[is.finite(pred$weight), , drop = FALSE]
  pred$weight <- abs(pred$weight)
  pred$edge_id <- metric_edge_ids(pred$regulator, pred$target)
  pred <- pred[order(-pred$weight), , drop = FALSE]
  pred <- pred[!duplicated(pred$edge_id), , drop = FALSE]
  pred$label <- as.integer(pred$edge_id %in% true_edge_ids)

  n_genes <- length(gt_genes)
  n_true <- length(true_edge_ids)
  n_possible <- as.double(n_genes) * as.double(max(n_genes - 1L, 0L))
  n_negative <- n_possible - n_true
  explicit_pos <- sum(pred$label)
  missing_pos <- n_true - explicit_pos
  missing_neg <- n_negative - (nrow(pred) - explicit_pos)
  if (missing_pos < 0 || missing_neg < 0) {
    stop("Network evaluation encountered an invalid edge accounting state.", call. = FALSE)
  }
  weights <- unique(pred$weight)
  block <- match(pred$weight, weights)
  grouped <- data.frame(
    weight = weights,
    pos = tabulate(block[pred$label == 1L], length(weights)),
    neg = tabulate(block[pred$label == 0L], length(weights))
  )
  if (missing_pos > 0 || missing_neg > 0) {
    zero <- which(grouped$weight == 0)
    if (length(zero)) {
      grouped$pos[zero[[1]]] <- grouped$pos[zero[[1]]] + missing_pos
      grouped$neg[zero[[1]]] <- grouped$neg[zero[[1]]] + missing_neg
    } else {
      grouped <- rbind(grouped, data.frame(weight = 0, pos = missing_pos, neg = missing_neg))
    }
  }
  grouped <- grouped[order(-grouped$weight), , drop = FALSE]
  tp_cum <- cumsum(grouped$pos)
  fp_cum <- cumsum(grouped$neg)

  positive_rate <- safe_metric_divide(n_true, n_possible, default = NA_real_)
  top_k <- as.integer(n_true)
  selected <- select_early_precision_edges(pred$weight, top_k)
  top_n <- length(selected)
  top_hits <- sum(pred$label[selected])
  early_precision <- safe_metric_divide(top_hits, top_n, default = 0)
  epr <- safe_metric_divide(early_precision, positive_rate, default = NA_real_)

  if (n_true == 0 || n_negative <= 0) {
    auroc <- NA_real_
  } else {
    ascending <- grouped[order(grouped$weight), , drop = FALSE]
    neg_before <- c(0, utils::head(cumsum(ascending$neg), -1))
    auroc <- safe_metric_divide(
      sum(ascending$pos * (neg_before + 0.5 * ascending$neg)),
      n_true * n_negative,
      default = NA_real_
    )
  }
  precision_desc <- ifelse((tp_cum + fp_cum) > 0, tp_cum / (tp_cum + fp_cum), 0)
  recall_desc <- tp_cum / n_true
  ap <- sum(grouped$pos / n_true * precision_desc)
  reverse <- rev(seq_along(precision_desc))
  pr_precision <- c(precision_desc[reverse], 1)
  pr_recall <- c(recall_desc[reverse], 0)
  auprc <- abs(sum(
    diff(pr_recall) *
      (utils::head(pr_precision, -1) + utils::tail(pr_precision, -1)) / 2
  ))

  if (n_negative <= 0) {
    tp_best <- 0
    fp_best <- 0
  } else {
    best <- which.max(as.double(tp_cum) * n_negative - as.double(fp_cum) * n_true)
    tp_best <- tp_cum[[best]]
    fp_best <- fp_cum[[best]]
  }
  tn_best <- n_negative - fp_best
  fn_best <- n_true - tp_best
  precision_best <- safe_metric_divide(tp_best, tp_best + fp_best, default = 0)
  recall_best <- safe_metric_divide(tp_best, tp_best + fn_best, default = 0)
  f1_best <- if ((precision_best + recall_best) == 0) {
    0
  } else {
    2 * precision_best * recall_best / (precision_best + recall_best)
  }
  acc_best <- safe_metric_divide(tp_best + tn_best, n_true + n_negative, default = 0)
  union_n <- (tp_best + fp_best) + n_true - tp_best
  ji_best <- if (union_n == 0) 0 else tp_best / union_n

  list(
    AUROC = as.numeric(auroc),
    AUPRC = as.numeric(auprc),
    AP = as.numeric(ap),
    ACC = as.numeric(acc_best),
    Precision = as.numeric(precision_best),
    Recall = as.numeric(recall_best),
    F1 = as.numeric(f1_best),
    JI = as.numeric(ji_best),
    SI = as.integer(tp_best),
    EarlyPrecision = early_precision,
    EPR = as.numeric(epr),
    TopK = top_k,
    EarlySelectedEdges = top_n,
    EarlyTruePositives = as.integer(top_hits),
    PositiveRate = positive_rate,
    GoldEdges = as.integer(n_true)
  )
}

select_early_precision_edges <- function(weights, top_k) {
  eligible <- which(is.finite(weights) & weights > 0)
  if (!length(eligible) || top_k <= 0) {
    return(integer(0))
  }
  ranked <- eligible[order(-weights[eligible])]
  cutoff <- weights[ranked[[min(as.integer(top_k), length(ranked))]]]
  eligible[weights[eligible] >= cutoff]
}
