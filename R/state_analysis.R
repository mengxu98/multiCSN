normalize_regulator_flag <- function(x) {
  if (is.null(x)) {
    return(logical(0))
  }
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  x_chr <- tolower(trimws(as.character(x)))
  !is.na(x_chr) & x_chr %in% c("true", "t", "1", "yes", "y")
}

choose_rank_value_column <- function(df, preferred = NULL) {
  preferred <- preferred[!is.na(preferred) & nzchar(preferred)]
  for (nm in preferred) {
    if (nm %in% colnames(df) && is.numeric(df[[nm]])) {
      return(nm)
    }
  }
  numeric_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
  numeric_cols <- setdiff(
    numeric_cols,
    c(
      "order_key",
      "n_targets_from",
      "n_targets_to",
      "n_regulators",
      "n_cells"
    )
  )
  numeric_cols[[1]] %ss% NULL
}

aggregate_target_weights <- function(df, weight_column = "weight") {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      target = character(0),
      x = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  target <- as.character(df$target)
  weight <- abs(as.numeric(df[[weight_column]]))
  keep <- !is.na(target) & nzchar(target) & !is.na(weight)
  if (!any(keep)) {
    return(data.frame(
      target = character(0),
      x = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- stats::aggregate(
    weight[keep],
    by = list(target = target[keep]),
    FUN = sum
  )
  out$target <- as.character(out$target)
  out
}

calculate_set_jaccard <- function(set1, set2) {
  set1 <- unique(as.character(set1))
  set2 <- unique(as.character(set2))
  set1 <- set1[!is.na(set1) & nzchar(set1)]
  set2 <- set2[!is.na(set2) & nzchar(set2)]
  union_n <- length(union(set1, set2))
  if (union_n == 0) {
    return(NA_real_)
  }
  length(intersect(set1, set2)) / union_n
}

compare_state_similarity <- function(
  object,
  network = DefaultNetwork(object),
  type = c("networks", "regulators", "targets"),
  state_ids = NULL,
  weight_column = "weight"
) {
  type <- match.arg(type)
  state_tbl <- summarize_networks(object, network = network)
  state_tbl <- state_tbl[state_tbl$keep & state_tbl$inferable, , drop = FALSE]
  state_tbl <- state_tbl[
    state_tbl$state_id %in% resolve_network_states(object, network = network), ,
    drop = FALSE
  ]
  if (!is.null(state_ids)) {
    state_tbl <- state_tbl[state_tbl$state_id %in% state_ids, , drop = FALSE]
  }
  state_tbl <- state_tbl[order(state_tbl$order_key), , drop = FALSE]
  if (nrow(state_tbl) < 2) {
    stop(
      "compare_state_similarity: need at least two kept inferable states.",
      call. = FALSE
    )
  }

  state_ids <- as.character(state_tbl$state_id)
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = state_ids,
    weight_cutoff = NULL,
    inferable_only = FALSE
  )

  feature_sets <- lapply(state_ids, function(sid) {
    net <- nets[[sid]]
    if (is.null(net) || nrow(net) == 0) {
      return(character(0))
    }
    if (identical(type, "networks")) {
      return(unique(paste(
        as.character(net$regulator),
        as.character(net$target),
        sep = "\t"
      )))
    }
    if (identical(type, "regulators")) {
      return(unique(as.character(net$regulator)))
    }
    unique(as.character(net$target))
  })
  names(feature_sets) <- state_ids

  sim_mat <- matrix(
    NA_real_,
    nrow = length(state_ids),
    ncol = length(state_ids)
  )
  rownames(sim_mat) <- state_ids
  colnames(sim_mat) <- state_ids
  for (i in seq_along(state_ids)) {
    for (j in seq_along(state_ids)) {
      sim_mat[i, j] <- calculate_set_jaccard(
        feature_sets[[i]],
        feature_sets[[j]]
      )
    }
  }

  counts <- data.frame(
    state_id = state_ids,
    feature_count = vapply(feature_sets, length, integer(1)),
    mean_similarity = vapply(
      seq_along(state_ids),
      function(i) {
        vals <- sim_mat[i, setdiff(seq_along(state_ids), i), drop = TRUE]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0) {
          return(NA_real_)
        }
        mean(vals)
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    similarity = sim_mat,
    counts = counts,
    state_table = state_tbl,
    feature_sets = feature_sets,
    type = type
  )
}

infer_feature_network_kind <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL
) {
  state_tbl <- tryCatch(
    summarize_networks(
      object,
      network = network,
      celltypes = celltypes
    ),
    error = function(e) NULL
  )
  if (!is.null(state_tbl) && nrow(state_tbl) > 0) {
    state_cols <- c("state_type", "keep", "inferable", "support_score")
    if (all(state_cols %in% colnames(state_tbl))) {
      stable_n <- sum(
        !is.na(state_tbl$state_type) & state_tbl$state_type == "stable"
      )
      if (stable_n > 0) {
        return("dynamic")
      }
    }
  }
  "static"
}

resolve_network_states <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  inferable_only = TRUE
) {
  networks <- GetNetwork(object, network = network)
  available_states <- names(networks) %ss% character(0)

  if (!is.null(celltypes)) {
    return(intersect(as.character(celltypes), available_states))
  }

  states <- tryCatch(
    summarize_networks(object, network = network),
    error = function(e) NULL
  )
  if (is.null(states) || nrow(states) == 0) {
    return(available_states)
  }

  if (isTRUE(inferable_only)) {
    states <- states[states$keep & states$inferable, , drop = FALSE]
  }

  intersect(as.character(states$state_id), available_states)
}

export_state_networks <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  weight_cutoff = NULL,
  inferable_only = FALSE
) {
  state_ids <- resolve_network_states(
    object,
    network = network,
    celltypes = celltypes,
    inferable_only = inferable_only
  )
  if (length(state_ids) == 0) {
    return(list())
  }

  res <- export_csn(
    object,
    active_network = network,
    celltypes = state_ids,
    weight_cutoff = weight_cutoff
  )

  if (is.null(res)) {
    return(list())
  }
  if (is.data.frame(res)) {
    res <- list(res)
    names(res) <- state_ids[[1]]
    return(res)
  }
  if (is.null(names(res))) {
    names(res) <- state_ids[seq_along(res)]
  }
  res
}

make_metric_long_table <- function(tbl, metrics) {
  pieces <- lapply(metrics, function(metric) {
    if (!metric %in% colnames(tbl)) {
      stop(
        sprintf("Metric '%s' not found in summary table.", metric),
        call. = FALSE
      )
    }
    data.frame(
      state_id = tbl$state_id,
      state_type = tbl$state_type,
      metric = metric,
      value = tbl[[metric]],
      order_key = tbl$order_key,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

prepare_feature_plot_table <- function(df, feature_col, value_col) {
  split_df <- split(df, df$state_id)
  out <- lapply(split_df, function(x) {
    x <- x[order(x[[value_col]], decreasing = TRUE), , drop = FALSE]
    x$feature_label <- paste(x$state_id, x[[feature_col]], sep = "___")
    x$feature_label <- factor(x$feature_label, levels = rev(x$feature_label))
    x
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

get_ordered_state_levels <- function(state_ids) {
  state_ids <- unique(as.character(state_ids))
  state_ids[order(vapply(state_ids, state_order_key, numeric(1)))]
}

get_shared_state_regulators <- function(
  regulator_table,
  state_ids = NULL,
  min_state_fraction = 1
) {
  if (is.null(regulator_table) || nrow(regulator_table) == 0) {
    return(character(0))
  }

  state_ids <- state_ids %ss% get_ordered_state_levels(regulator_table$state_id)
  if (length(state_ids) == 0) {
    return(character(0))
  }

  present <- stats::aggregate(
    as.integer(regulator_table$state_id %in% state_ids),
    by = list(gene = regulator_table$gene, state_id = regulator_table$state_id),
    FUN = max
  )
  wide <- stats::xtabs(x ~ gene + state_id, data = present)
  state_presence <- rowSums(wide > 0)
  keep_n <- ceiling(length(state_ids) * min_state_fraction)
  names(state_presence)[state_presence >= keep_n]
}

assign_cells_to_states <- function(
  object,
  pseudotime_column = NULL,
  keep_only = TRUE,
  inferable_only = FALSE
) {
  states <- summarize_networks(
    object,
    pseudotime_column = pseudotime_column
  )
  if (isTRUE(keep_only)) {
    states <- states[states$keep, , drop = FALSE]
  }
  if (isTRUE(inferable_only)) {
    states <- states[states$inferable, , drop = FALSE]
  }
  if (nrow(states) == 0) {
    return(data.frame())
  }

  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )

  pieces <- lapply(states$state_id, function(sid) {
    cells <- density_result$cells[[sid]]
    if (is.null(cells) || length(cells) == 0) {
      return(NULL)
    }
    row <- states[states$state_id == sid, , drop = FALSE]
    data.frame(
      cell = cells,
      state_id = sid,
      state_type = row$state_type[[1]],
      order_key = row$order_key[[1]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

fetch_tf_expression_data <- function(
  object,
  tf,
  pseudotime_column = NULL,
  assay = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "fetch_tf_expression_data: object must be a Seurat object.",
      call. = FALSE
    )
  }
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  assay <- assay %ss% Seurat::DefaultAssay(object)

  expr_mat <- Seurat::GetAssayData(
    object,
    assay = assay,
    layer = "data"
  )
  if (!tf %in% rownames(expr_mat)) {
    stop(
      sprintf("TF '%s' not found in assay '%s'.", tf, assay),
      call. = FALSE
    )
  }

  cells <- colnames(object)
  meta <- object@meta.data[cells, , drop = FALSE]
  require_meta_column(meta, pseudotime_column, what = "pseudotime_column")
  expr <- as.numeric(expr_mat[tf, cells, drop = TRUE])

  out <- data.frame(
    cell = cells,
    pseudotime = meta[[pseudotime_column]],
    expression = expr,
    stringsAsFactors = FALSE
  )

  state_map <- assign_cells_to_states(
    object,
    pseudotime_column = pseudotime_column,
    keep_only = TRUE,
    inferable_only = FALSE
  )
  if (nrow(state_map) > 0) {
    idx <- match(out$cell, state_map$cell)
    out$state_id <- state_map$state_id[idx]
    out$state_type <- state_map$state_type[idx]
    out$order_key <- state_map$order_key[idx]
  } else {
    out$state_id <- NA_character_
    out$state_type <- NA_character_
    out$order_key <- NA_real_
  }

  out <- out[is.finite(out$pseudotime), , drop = FALSE]
  out
}

extract_tf_target_table <- function(
  object,
  tf,
  network = DefaultNetwork(object),
  celltypes = NULL,
  weight_column = "weight"
) {
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = celltypes,
    inferable_only = FALSE
  )
  if (length(nets) == 0) {
    return(data.frame())
  }

  out <- purrr::imap_dfr(nets, function(df, sid) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df <- df[df$regulator == tf, , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    if (!weight_column %in% colnames(df)) {
      stop(
        sprintf(
          "extract_tf_target_table: weight_column '%s' not found.",
          weight_column
        ),
        call. = FALSE
      )
    }

    agg_weight <- tapply(df[[weight_column]], df$target, sum, na.rm = TRUE)
    agg_abs <- tapply(abs(df[[weight_column]]), df$target, sum, na.rm = TRUE)
    agg_n <- tapply(df[[weight_column]], df$target, length)

    data.frame(
      state_id = sid,
      target = names(agg_weight),
      weight = as.numeric(agg_weight),
      abs_weight = as.numeric(agg_abs[names(agg_weight)]),
      n_edges = as.numeric(agg_n[names(agg_weight)]),
      stringsAsFactors = FALSE
    )
  })

  if (nrow(out) == 0) {
    return(out)
  }
  out$order_key <- vapply(
    as.character(out$state_id),
    state_order_key,
    numeric(1)
  )
  out <- out[order(out$order_key, -out$abs_weight), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Internal helper: integrated regulator prioritization
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param assay RNA assay used for TF expression summaries.
#' @param peak_assay Chromatin assay used for motif-accessibility summaries.
#' @param method Regulator ranking method.
#' @param require_all_states Require TFs to be present in all inferable states.
#' @param min_state_fraction Minimum fraction of states where TF must be present when \code{require_all_states = FALSE}.
#' @param min_targets Minimum number of targets in at least one state.
#'
#' @return A scored data.frame of key TF candidates.
#' @noRd
.prioritize_regulators_integrated <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  method = c("page_rank", "degree_distribution"),
  require_all_states = TRUE,
  min_state_fraction = 1,
  min_targets = 10
) {
  if (!methods::is(object, "Seurat")) {
    stop("prioritize_features: object must be a Seurat object.", call. = FALSE)
  }
  method <- match.arg(method)
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  assay <- assay %ss%
    (tryCatch(Params(object)$rna_assay, error = function(e) NULL) %ss%
      Seurat::DefaultAssay(object))
  peak_assay <- resolve_peak_assay_name(object, peak_assay = peak_assay)

  reg_tbl <- rank_features(
    object,
    type = "regulators",
    network = network,
    method = method,
    top_n = NULL
  )
  if (nrow(reg_tbl) == 0) {
    return(data.frame())
  }
  score_col <- choose_rank_value_column(
    reg_tbl,
    preferred = c(
      "page_rank",
      "rank_value",
      "degree",
      "out_degree",
      "in_degree"
    )
  )
  if (is.null(score_col)) {
    stop(
      "prioritize_features: unable to identify regulator score column.",
      call. = FALSE
    )
  }

  state_ids <- get_ordered_state_levels(reg_tbl$state_id)
  if (isTRUE(require_all_states)) {
    shared_tfs <- get_shared_state_regulators(
      reg_tbl,
      state_ids = state_ids,
      min_state_fraction = 1
    )
  } else {
    shared_tfs <- get_shared_state_regulators(
      reg_tbl,
      state_ids = state_ids,
      min_state_fraction = min_state_fraction
    )
  }
  reg_tbl <- reg_tbl[reg_tbl$gene %in% shared_tfs, , drop = FALSE]
  if (nrow(reg_tbl) == 0) {
    return(data.frame())
  }

  tf_target_counts <- get_state_regulator_target_counts(
    object,
    network = network,
    state_ids = state_ids
  )
  tf_target_summary <- stats::aggregate(
    tf_target_counts$n_targets,
    by = list(gene = tf_target_counts$regulator),
    FUN = function(x) max(x, na.rm = TRUE)
  )
  colnames(tf_target_summary)[2] <- "max_targets"
  tf_keep <- tf_target_summary$gene[
    tf_target_summary$max_targets >= min_targets
  ]
  reg_tbl <- reg_tbl[reg_tbl$gene %in% tf_keep, , drop = FALSE]
  if (nrow(reg_tbl) == 0) {
    return(data.frame())
  }

  genes <- unique(as.character(reg_tbl$gene))
  centrality_mat <- matrix(0, nrow = length(genes), ncol = length(state_ids))
  rownames(centrality_mat) <- genes
  colnames(centrality_mat) <- state_ids
  for (sid in state_ids) {
    sub <- reg_tbl[reg_tbl$state_id == sid, , drop = FALSE]
    if (nrow(sub) == 0) {
      next
    }
    agg <- stats::aggregate(
      sub[[score_col]],
      by = list(gene = sub$gene),
      FUN = function(x) max(x, na.rm = TRUE)
    )
    centrality_mat[agg$gene, sid] <- agg$x
  }

  expression_mat <- vapply(
    state_ids,
    function(sid) {
      vals <- compute_state_gene_expression(
        object,
        state_id = sid,
        genes = genes,
        pseudotime_column = pseudotime_column,
        assay = assay
      )
      out <- rep(0, length(genes))
      names(out) <- genes
      out[names(vals)] <- vals
      out
    },
    numeric(length(genes))
  )
  rownames(expression_mat) <- genes

  motif_mat <- compute_state_tf_motif_accessibility(
    object,
    state_ids = state_ids,
    tfs = genes,
    pseudotime_column = pseudotime_column,
    peak_assay = peak_assay
  )
  motif_mat <- motif_mat[genes, state_ids, drop = FALSE]

  rewiring <- compare_state_rewiring(object, network = network)
  rew_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df[df$regulator %in% genes, , drop = FALSE]
  })
  rew_summary <- if (nrow(rew_tbl) > 0) {
    stats::aggregate(
      rew_tbl$rewiring_score,
      by = list(gene = rew_tbl$regulator),
      FUN = function(x) mean(x, na.rm = TRUE)
    )
  } else {
    data.frame(gene = genes, x = 0, stringsAsFactors = FALSE)
  }
  colnames(rew_summary)[2] <- "rewiring_mean"

  tf_tbl <- data.frame(
    gene = genes,
    centrality_peak = apply(centrality_mat, 1, max, na.rm = TRUE),
    centrality_dynamic = apply(centrality_mat, 1, function(x) {
      diff(range(x, na.rm = TRUE))
    }),
    motif_accessibility_peak = apply(motif_mat, 1, max, na.rm = TRUE),
    expression_peak = apply(expression_mat, 1, max, na.rm = TRUE),
    state_specificity = vapply(
      genes,
      function(g) specificity_index(centrality_mat[g, ]),
      numeric(1)
    ),
    anchor_state = state_ids[max.col(
      centrality_mat[genes, , drop = FALSE],
      ties.method = "first"
    )],
    stringsAsFactors = FALSE
  )
  tf_tbl$rewiring_mean <- rew_summary$rewiring_mean[match(
    tf_tbl$gene,
    rew_summary$gene
  )]
  tf_tbl$rewiring_mean[is.na(tf_tbl$rewiring_mean)] <- 0
  tf_tbl$max_targets <- tf_target_summary$max_targets[match(
    tf_tbl$gene,
    tf_target_summary$gene
  )]
  tf_tbl$max_targets[is.na(tf_tbl$max_targets)] <- 0

  tf_tbl$centrality_peak_z <- safe_zscore(tf_tbl$centrality_peak)
  tf_tbl$centrality_dynamic_z <- safe_zscore(tf_tbl$centrality_dynamic)
  tf_tbl$rewiring_mean_z <- safe_zscore(tf_tbl$rewiring_mean)
  tf_tbl$motif_accessibility_peak_z <- safe_zscore(
    tf_tbl$motif_accessibility_peak
  )
  tf_tbl$expression_peak_z <- safe_zscore(tf_tbl$expression_peak)
  tf_tbl$state_specificity_z <- safe_zscore(tf_tbl$state_specificity)

  tf_tbl$key_tf_score <-
    0.30 *
    tf_tbl$centrality_peak_z +
    0.15 * tf_tbl$centrality_dynamic_z +
    0.20 * tf_tbl$rewiring_mean_z +
    0.15 * tf_tbl$motif_accessibility_peak_z +
    0.10 * tf_tbl$expression_peak_z +
    0.10 * tf_tbl$state_specificity_z

  tf_tbl$anchor_order_key <- vapply(
    tf_tbl$anchor_state,
    state_order_key,
    numeric(1)
  )
  tf_tbl$type <- "regulators"
  tf_tbl$strategy <- "integrated"
  tf_tbl$network_kind <- "dynamic"
  tf_tbl$feature <- tf_tbl$gene
  tf_tbl$feature_score <- tf_tbl$key_tf_score
  tf_tbl$evidence_used <- if (is.null(peak_assay) || nrow(motif_mat) == 0) {
    "network,rewiring,specificity,expression"
  } else {
    "network,rewiring,specificity,expression,motif"
  }
  tf_tbl <- tf_tbl[
    order(tf_tbl$key_tf_score, decreasing = TRUE), ,
    drop = FALSE
  ]
  rownames(tf_tbl) <- NULL
  tf_tbl
}

#' @title Select top prioritized TFs
#' @inheritParams prioritize_features
#' @param pseudotime_column Optional pseudotime column.
#' @param assay RNA assay used for TF expression summaries.
#' @param peak_assay Chromatin assay used for motif-accessibility summaries.
#' @param method Regulator ranking method.
#' @param require_all_states Require TFs to be present in all inferable states.
#' @param min_state_fraction Minimum fraction of states where TF must be present
#'   when \code{require_all_states = FALSE}.
#' @param min_targets Minimum number of targets in at least one state.
#' @param top_n Number of TFs returned.
#' @return Character vector of TFs.
#' @export
select_key_tfs <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  method = c("page_rank", "degree_distribution"),
  require_all_states = TRUE,
  min_state_fraction = 1,
  min_targets = 10,
  top_n = 6
) {
  tbl <- prioritize_features(
    object = object,
    type = "regulators",
    strategy = "integrated",
    network = network,
    pseudotime_column = pseudotime_column,
    assay = assay,
    peak_assay = peak_assay,
    method = method,
    require_all_states = require_all_states,
    min_state_fraction = min_state_fraction,
    min_targets = min_targets
  )
  tbl$gene[seq_len(min(top_n, nrow(tbl)))]
}

#' @title Prioritize network features using network or integrated evidence
#'
#' @param object A \code{Seurat} object.
#' @param type Feature type, \code{"regulators"} or \code{"targets"}.
#' @param strategy Prioritization strategy, \code{"network"} or \code{"integrated"}.
#' @param network Dynamic or static network name. Defaults to \code{DefaultNetwork(object)}.
#' @param ... Additional arguments forwarded to the underlying ranking or integrated-prioritization routine.
#'
#' @return A data.frame of prioritized features.
#' @export
prioritize_features <- function(
  object,
  type = c("regulators", "targets"),
  strategy = c("network", "integrated"),
  network = DefaultNetwork(object),
  ...
) {
  type <- match.arg(type)
  strategy <- match.arg(strategy)

  if (identical(strategy, "network")) {
    return(rank_features(
      object = object,
      type = type,
      network = network,
      ...
    ))
  }

  network_kind <- infer_feature_network_kind(object, network = network)
  if (!identical(network_kind, "dynamic")) {
    stop(
      "prioritize_features: integrated strategy is currently only supported for dynamic networks.",
      call. = FALSE
    )
  }

  if (identical(type, "regulators")) {
    .prioritize_regulators_integrated(
      object = object,
      network = network,
      ...
    )
  } else {
    .prioritize_targets_integrated(
      object = object,
      network = network,
      ...
    )
  }
}

#' Internal helper: integrated target prioritization
#'
#' @param object A \code{Seurat} object.
#' @param tfs Optional TF vector. If \code{NULL}, derived from \code{tf_scores}.
#' @param tf_scores Optional output of \code{prioritize_features(type = "regulators", strategy = "integrated")}.
#' @param network Dynamic network name.
#' @param pseudotime_column Optional pseudotime column.
#' @param assay RNA assay used for target-expression support.
#' @param peak_assay Chromatin assay used for peak support.
#' @param weight_column Edge-weight column.
#' @param dynamic_only Restrict candidate targets to dynamic genes when available.
#' @param dynamic_padjust_threshold Adjusted P-value threshold for dynamic genes.
#' @param peak_to_gene_method Peak-to-gene linking method.
#' @param upstream,downstream,only_tss Peak-to-gene window parameters.
#' @param top_targets_per_tf Optional truncation after scoring.
#'
#' @return A data.frame of scored TF-target pairs.
#' @noRd
.prioritize_targets_integrated <- function(
  object,
  tfs = NULL,
  tf_scores = NULL,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  weight_column = "weight",
  dynamic_only = TRUE,
  dynamic_padjust_threshold = 0.05,
  peak_to_gene_method = "Signac",
  upstream = 100000,
  downstream = 0,
  only_tss = FALSE,
  top_targets_per_tf = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop("prioritize_features: object must be a Seurat object.", call. = FALSE)
  }
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  assay <- assay %ss%
    (tryCatch(Params(object)$rna_assay, error = function(e) NULL) %ss%
      Seurat::DefaultAssay(object))
  peak_assay <- resolve_peak_assay_name(object, peak_assay = peak_assay)

  if (is.null(tf_scores) || nrow(tf_scores) == 0) {
    tf_scores <- .prioritize_regulators_integrated(
      object = object,
      network = network,
      pseudotime_column = pseudotime_column,
      assay = assay,
      peak_assay = peak_assay
    )
  }
  if (!is.null(tfs)) {
    tf_scores <- tf_scores[tf_scores$gene %in% tfs, , drop = FALSE]
  }
  if (nrow(tf_scores) == 0) {
    return(data.frame())
  }

  state_ids <- get_ordered_state_levels(tf_scores$anchor_state)
  nets_all <- export_state_networks(
    object,
    network = network,
    inferable_only = TRUE
  )
  edge_tbl <- purrr::imap_dfr(nets_all, function(df, sid) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    if (!weight_column %in% colnames(df)) {
      stop(
        "prioritize_features: weight_column not found in exported network.",
        call. = FALSE
      )
    }
    if (!all(c("regulator", "target") %in% colnames(df))) {
      stop(
        "prioritize_features: exported network must contain regulator and target columns.",
        call. = FALSE
      )
    }
    edge_input <- df[, c("regulator", "target", weight_column), drop = FALSE]
    colnames(edge_input)[3] <- "weight"
    edge_input$state_id <- sid
    edge_input$regulator <- as.character(edge_input$regulator)
    edge_input$target <- as.character(edge_input$target)
    edge_input$weight <- suppressWarnings(as.numeric(edge_input$weight))
    edge_input <- edge_input[
      !is.na(edge_input$regulator) &
        nzchar(edge_input$regulator) &
        !is.na(edge_input$target) &
        nzchar(edge_input$target) &
        !is.na(edge_input$weight),
      c("state_id", "regulator", "target", "weight"),
      drop = FALSE
    ]
    if (nrow(edge_input) == 0) {
      return(NULL)
    }
    edge_input$abs_weight <- abs(edge_input$weight)
    agg <- stats::aggregate(
      abs_weight ~ state_id + regulator + target,
      data = edge_input,
      FUN = sum,
      na.rm = TRUE
    )
    agg
  })
  if (nrow(edge_tbl) == 0) {
    return(data.frame())
  }
  edge_total <- stats::aggregate(
    edge_tbl$abs_weight,
    by = list(regulator = edge_tbl$regulator, target = edge_tbl$target),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(edge_total)[3] <- "total_abs_weight"

  dynamic_genes <- get_dynamic_feature_genes(
    object,
    pseudotime_column = pseudotime_column,
    padjust_threshold = dynamic_padjust_threshold
  )

  anchor_targets <- lapply(seq_len(nrow(tf_scores)), function(i) {
    tf <- tf_scores$gene[[i]]
    sid <- tf_scores$anchor_state[[i]]
    sub <- edge_tbl[
      edge_tbl$regulator == tf & edge_tbl$state_id == sid, ,
      drop = FALSE
    ]
    if (isTRUE(dynamic_only) && length(dynamic_genes) > 0) {
      sub <- sub[sub$target %in% dynamic_genes, , drop = FALSE]
    }
    unique(sub$target)
  })
  names(anchor_targets) <- tf_scores$gene
  candidate_targets <- unique(unlist(anchor_targets, use.names = FALSE))
  if (length(candidate_targets) == 0) {
    return(data.frame())
  }

  gene_peak_map <- build_gene_peak_map(
    object,
    genes = candidate_targets,
    peak_assay = peak_assay,
    method = peak_to_gene_method,
    upstream = upstream,
    downstream = downstream,
    only_tss = only_tss
  )
  peak_tf <- build_peak_tf_incidence(object, tfs = tf_scores$gene)
  if (!is.null(gene_peak_map) && !is.null(peak_tf)) {
    common_peaks <- intersect(rownames(gene_peak_map), rownames(peak_tf))
    gene_peak_map <- gene_peak_map[common_peaks, , drop = FALSE]
    peak_tf <- peak_tf[common_peaks, tf_scores$gene, drop = FALSE]
  }

  target_rank_tbl <- rank_features(
    object,
    type = "targets",
    network = network,
    top_n = NULL,
    weight_column = weight_column
  )

  state_peak_acc_cache <- list()
  if (!is.null(gene_peak_map) && nrow(gene_peak_map) > 0) {
    peak_names <- rownames(gene_peak_map)
    for (sid in unique(tf_scores$anchor_state)) {
      acc <- compute_state_peak_accessibility(
        object,
        state_id = sid,
        peaks = peak_names,
        pseudotime_column = pseudotime_column,
        peak_assay = peak_assay
      )
      acc <- acc[peak_names]
      acc[is.na(acc)] <- 0
      state_peak_acc_cache[[sid]] <- acc
    }
  }

  out <- lapply(seq_len(nrow(tf_scores)), function(i) {
    tf <- tf_scores$gene[[i]]
    sid <- tf_scores$anchor_state[[i]]
    sub <- edge_tbl[
      edge_tbl$regulator == tf & edge_tbl$state_id == sid, ,
      drop = FALSE
    ]
    if (isTRUE(dynamic_only) && length(dynamic_genes) > 0) {
      sub <- sub[sub$target %in% dynamic_genes, , drop = FALSE]
    }
    if (nrow(sub) == 0) {
      return(NULL)
    }
    sub <- merge(sub, edge_total, by = c("regulator", "target"), all.x = TRUE)
    colnames(sub)[colnames(sub) == "regulator"] <- "tf"
    sub$target_specificity <- sub$abs_weight /
      pmax(sub$total_abs_weight, .Machine$double.eps)

    expr_vals <- compute_state_gene_expression(
      object,
      state_id = sid,
      genes = sub$target,
      pseudotime_column = pseudotime_column,
      assay = assay
    )
    sub$target_expression <- expr_vals[sub$target]
    sub$target_expression[is.na(sub$target_expression)] <- 0

    tgt_rank_sub <- target_rank_tbl[
      target_rank_tbl$state_id == sid, ,
      drop = FALSE
    ]
    sub$target_weighted_indegree <- tgt_rank_sub$weighted_indegree[match(
      sub$target,
      tgt_rank_sub$target
    )]
    sub$target_weighted_indegree[is.na(sub$target_weighted_indegree)] <- 0

    sub$peak_support <- 0
    sub$n_support_peaks <- 0
    if (
      !is.null(gene_peak_map) &&
        nrow(gene_peak_map) > 0 &&
        !is.null(peak_tf) &&
        tf %in% colnames(peak_tf)
    ) {
      genes_use <- intersect(sub$target, colnames(gene_peak_map))
      if (length(genes_use) > 0) {
        motif_idx <- which(as.vector(peak_tf[, tf]) > 0)
        if (length(motif_idx) > 0) {
          gp_sub <- gene_peak_map[motif_idx, genes_use, drop = FALSE]
          peak_acc <- state_peak_acc_cache[[sid]][rownames(gp_sub)]
          support_sum <- as.numeric(Matrix::crossprod(peak_acc, gp_sub))
          support_n <- as.numeric(Matrix::colSums(gp_sub != 0))
          support_mean <- support_sum / pmax(support_n, 1)
          names(support_mean) <- colnames(gp_sub)
          names(support_n) <- colnames(gp_sub)
          sub$peak_support <- support_mean[sub$target]
          sub$peak_support[is.na(sub$peak_support)] <- 0
          sub$n_support_peaks <- support_n[sub$target]
          sub$n_support_peaks[is.na(sub$n_support_peaks)] <- 0
        }
      }
    }
    sub$dynamic_supported <- if (length(dynamic_genes) > 0) {
      sub$target %in% dynamic_genes
    } else {
      TRUE
    }
    sub$peak_supported <- sub$n_support_peaks > 0
    sub$anchor_state <- sid
    sub
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }

  split_out <- split(out, out$tf)
  split_out <- lapply(split_out, function(df) {
    df$abs_weight_z <- safe_zscore(df$abs_weight)
    df$target_specificity_z <- safe_zscore(df$target_specificity)
    df$target_weighted_indegree_z <- safe_zscore(df$target_weighted_indegree)
    df$target_expression_z <- safe_zscore(df$target_expression)
    df$peak_support_z <- safe_zscore(df$peak_support)
    df$target_score <-
      0.35 *
      df$abs_weight_z +
      0.20 * df$target_specificity_z +
      0.15 * df$target_weighted_indegree_z +
      0.10 * df$target_expression_z +
      0.20 * df$peak_support_z
    df <- df[order(df$target_score, decreasing = TRUE), , drop = FALSE]
    if (!is.null(top_targets_per_tf) && nrow(df) > top_targets_per_tf) {
      df <- df[seq_len(top_targets_per_tf), , drop = FALSE]
    }
    df
  })
  out <- do.call(rbind, split_out)
  out$type <- "targets"
  out$strategy <- "integrated"
  out$network_kind <- "dynamic"
  out$feature <- out$target
  out$feature_score <- out$target_score
  out$evidence_used <- if (!all(out$n_support_peaks == 0)) {
    "network,specificity,expression,peak"
  } else {
    "network,specificity,expression"
  }
  rownames(out) <- NULL
  out
}

#' @title Select top prioritized targets for chosen TFs
#' @inheritParams prioritize_features
#' @param tfs Optional TF vector. If \code{NULL}, derived from \code{tf_scores}.
#' @param tf_scores Optional output of
#'   \code{prioritize_features(type = "regulators", strategy = "integrated")}.
#' @param pseudotime_column Optional pseudotime column.
#' @param assay RNA assay used for target-expression support.
#' @param peak_assay Chromatin assay used for peak support.
#' @param top_n Number of targets kept per TF.
#' @return A data.frame of top targets per TF.
#' @export
select_key_targets <- function(
  object,
  tfs = NULL,
  tf_scores = NULL,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  top_n = 10
) {
  prioritize_features(
    object = object,
    type = "targets",
    strategy = "integrated",
    tfs = tfs,
    tf_scores = tf_scores,
    network = network,
    pseudotime_column = pseudotime_column,
    assay = assay,
    peak_assay = peak_assay,
    top_targets_per_tf = top_n
  )
}

#' @title Plot integrated multiome key-TF scores
#' @param tf_scores Output of \code{prioritize_features(type = "regulators", strategy = "integrated")}.
#' @param top_n Number of TFs shown.
#' @return A patchwork/ggplot object.
#' @export
plot_key_tf_multiome <- function(tf_scores, top_n = 12) {
  if (is.null(tf_scores) || nrow(tf_scores) == 0) {
    stop("plot_key_tf_multiome: tf_scores is empty.", call. = FALSE)
  }
  plot_df <- tf_scores[seq_len(min(top_n, nrow(tf_scores))), , drop = FALSE]
  plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene))

  p_rank <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = gene, y = key_tf_score, fill = key_tf_score)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "Integrated key-TF ranking",
      x = NULL,
      y = "KeyTFScore",
      fill = "Score"
    )

  component_cols <- c(
    "centrality_peak_z",
    "centrality_dynamic_z",
    "rewiring_mean_z",
    "motif_accessibility_peak_z",
    "expression_peak_z",
    "state_specificity_z"
  )
  component_labels <- c(
    centrality_peak_z = "Centrality peak",
    centrality_dynamic_z = "Centrality dynamic",
    rewiring_mean_z = "Rewiring",
    motif_accessibility_peak_z = "Motif accessibility",
    expression_peak_z = "TF expression",
    state_specificity_z = "State specificity"
  )
  heat_df <- do.call(
    rbind,
    lapply(component_cols, function(comp) {
      data.frame(
        gene = plot_df$gene,
        component = component_labels[[comp]],
        value = plot_df[[comp]],
        stringsAsFactors = FALSE
      )
    })
  )
  heat_df$component <- factor(
    heat_df$component,
    levels = component_labels[component_cols]
  )
  p_heat <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = component, y = gene, fill = value)
  ) +
    ggplot2::geom_tile(color = "grey92") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "Score components", x = NULL, y = NULL, fill = "z")

  p <- p_rank | p_heat
  print(p)
  invisible(p)
}

#' @title Plot integrated key-target ranking with ATAC support
#' @param target_scores Output of \code{prioritize_features(type = "targets", strategy = "integrated")}.
#' @param top_n Number of targets shown per TF.
#' @return A \code{ggplot} object.
#' @export
plot_key_target_multiome <- function(target_scores, top_n = 10) {
  if (is.null(target_scores) || nrow(target_scores) == 0) {
    stop("plot_key_target_multiome: target_scores is empty.", call. = FALSE)
  }
  split_df <- split(target_scores, target_scores$tf)
  plot_df <- lapply(split_df, function(df) {
    df <- df[order(df$target_score, decreasing = TRUE), , drop = FALSE]
    df <- df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
    df$target_label <- paste(df$tf, df$target, sep = "___")
    df$target_label <- factor(df$target_label, levels = rev(df$target_label))
    df
  })
  plot_df <- do.call(rbind, plot_df)
  rownames(plot_df) <- NULL

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = target_label,
      y = target_score,
      color = peak_support,
      size = pmax(abs_weight, 1e-8),
      shape = peak_supported
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = target_label, y = 0, yend = target_score),
      color = "grey75"
    ) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::facet_grid(tf ~ ., scales = "free_y", space = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_color_gradient(low = "grey85", high = "#B2182B") +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Integrated key-target ranking",
      x = NULL,
      y = "TargetScore",
      color = "Peak support",
      size = "Edge weight",
      shape = "ATAC-supported"
    )
  print(p)
  invisible(p)
}

build_state_density_plot <- function(
  object,
  pseudotime_column = NULL,
  group_column = NULL,
  palette = NULL,
  palette_name = "Chinese"
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "build_state_density_plot: object must be a Seurat object.",
      call. = FALSE
    )
  }

  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )
  state_windows <- density_result$windows
  state_windows <- state_windows[state_windows$keep, , drop = FALSE]
  prepared <- prepare_density_input(
    object = object,
    pseudotime_column = pseudotime_column,
    group_column = group_column %ss% density_result$params$group_column
  )
  plot_data <- prepared$meta[, c(".group", pseudotime_column), drop = FALSE]
  colnames(plot_data) <- c("cluster", "pseudotime")
  plot_data <- plot_data[
    !is.na(plot_data$cluster) & plot_data$cluster != "", ,
    drop = FALSE
  ]
  if (nrow(plot_data) == 0) {
    stop(
      "build_state_density_plot: no cells with valid group labels.",
      call. = FALSE
    )
  }

  group_order <- density_result$group_order
  plot_data$cluster <- factor(plot_data$cluster, levels = group_order)
  if (is.null(palette)) {
    palette <- .multicsn_palette_colors(group_order, palette = palette_name)
  } else if (is.null(names(palette))) {
    palette <- stats::setNames(
      rep(palette, length.out = length(group_order)),
      group_order
    )
  }
  palette <- palette[group_order]

  max_dens <- max(vapply(
    split(plot_data$pseudotime, plot_data$cluster),
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        return(0)
      }
      max(stats::density(x)$y)
    },
    numeric(1)
  ))

  stable_windows <- state_windows[
    state_windows$state_type == "stable", ,
    drop = FALSE
  ]
  transition_windows <- state_windows[
    state_windows$state_type == "transition", ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pseudotime))
  if (nrow(stable_windows) > 0) {
    stable_windows$fill_group <- vapply(
      stable_windows$parent_groups,
      function(x) strsplit(x, "\\|", fixed = FALSE)[[1]][1],
      character(1)
    )
    p <- p +
      ggplot2::geom_rect(
        data = stable_windows,
        ggplot2::aes(
          xmin = left,
          xmax = right,
          ymin = 0,
          ymax = Inf,
          fill = fill_group
        ),
        alpha = 0.15,
        inherit.aes = FALSE
      )
  }
  if (nrow(transition_windows) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = transition_windows,
        ggplot2::aes(xmin = left, xmax = right, ymin = 0, ymax = Inf),
        fill = "grey60",
        alpha = 0.18,
        inherit.aes = FALSE
      )
  }

  line_data <- state_windows[
    order(vapply(
      as.character(state_windows$state_id),
      state_order_key,
      numeric(1)
    )), ,
    drop = FALSE
  ]
  n_states <- nrow(line_data)
  n_spaces <- n_states + 2L
  line_data$y <- max_dens * (n_spaces - seq_len(n_states)) / n_spaces
  line_data$label_x <- (line_data$left + line_data$right) / 2
  line_data$label_y <- line_data$y + max_dens * 0.02

  p <- p +
    ggplot2::geom_density(
      ggplot2::aes(color = cluster, fill = cluster),
      alpha = 0.55
    ) +
    ggplot2::geom_vline(
      data = density_result$boundaries,
      ggplot2::aes(xintercept = boundary),
      linetype = "dashed",
      color = "grey40",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = line_data,
      ggplot2::aes(x = left, xend = right, y = y, yend = y),
      color = "grey20",
      linewidth = 0.6,
      arrow = ggplot2::arrow(angle = 15, type = "closed"),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = line_data,
      ggplot2::aes(x = label_x, y = label_y, label = state_id),
      inherit.aes = FALSE,
      size = 3.2,
      vjust = 0
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x = pseudotime_column,
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      panel.grid = ggplot2::element_blank()
    )

  p
}

get_state_regulator_target_counts <- function(
  object,
  network = DefaultNetwork(object),
  state_ids = NULL
) {
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = state_ids,
    inferable_only = TRUE
  )
  if (length(nets) == 0) {
    return(data.frame())
  }

  out <- purrr::imap_dfr(nets, function(df, sid) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    tgt_counts <- tapply(df$target, df$regulator, function(x) length(unique(x)))
    data.frame(
      state_id = sid,
      regulator = names(tgt_counts),
      n_targets = as.numeric(tgt_counts),
      stringsAsFactors = FALSE
    )
  })
  if (nrow(out) == 0) {
    return(out)
  }
  out$order_key <- vapply(
    as.character(out$state_id),
    state_order_key,
    numeric(1)
  )
  out <- out[order(out$order_key, -out$n_targets), , drop = FALSE]
  rownames(out) <- NULL
  out
}

compute_state_gene_expression <- function(
  object,
  state_id,
  genes,
  pseudotime_column = NULL,
  assay = NULL
) {
  if (length(genes) == 0) {
    return(setNames(numeric(0), character(0)))
  }

  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  assay <- assay %ss% Seurat::DefaultAssay(object)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )
  state_cells <- density_result$cells[[state_id]]
  state_cells <- intersect(state_cells, colnames(object))
  if (length(state_cells) == 0) {
    return(stats::setNames(rep(0, length(genes)), genes))
  }

  expr_mat <- Seurat::GetAssayData(
    object,
    assay = assay,
    layer = "data"
  )
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) {
    return(numeric(0))
  }

  vals <- Matrix::rowMeans(expr_mat[genes, state_cells, drop = FALSE])
  vals <- as.numeric(vals)
  names(vals) <- genes
  vals
}

safe_zscore <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1 || stats::sd(x, na.rm = TRUE) == 0 || all(!is.finite(x))) {
    out <- rep(0, length(x))
  } else {
    out <- as.numeric(scale(x))
    out[!is.finite(out)] <- 0
  }
  out
}

specificity_index <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  if (length(x) <= 1) {
    return(0)
  }
  x <- x - min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    return(0)
  }
  x <- x / xmax
  sum(1 - x) / (length(x) - 1)
}

resolve_peak_assay_name <- function(object, peak_assay = NULL) {
  params <- tryCatch(Params(object), error = function(e) list())
  peak_assay <- peak_assay %ss% params$peak_assay
  assay_names <- names(object@assays)
  if (!is.null(peak_assay) && peak_assay %in% assay_names) {
    return(peak_assay)
  }
  preferred <- assay_names[tolower(assay_names) %in% c("atac", "peaks", "peak")]
  preferred[[1]] %ss% NULL
}

get_assay_matrix_safe <- function(
  object,
  assay,
  layer_primary = "data",
  layer_fallback = "counts"
) {
  mat <- tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer_primary),
    error = function(e) NULL
  )
  if (!is.null(mat)) {
    return(mat)
  }
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer_fallback),
    error = function(e) {
      Seurat::GetAssayData(object, assay = assay, slot = layer_primary)
    }
  )
}

compute_state_peak_accessibility <- function(
  object,
  state_id,
  peaks,
  pseudotime_column = NULL,
  peak_assay = NULL
) {
  if (length(peaks) == 0) {
    return(setNames(numeric(0), character(0)))
  }

  peak_assay <- resolve_peak_assay_name(object, peak_assay = peak_assay)
  if (is.null(peak_assay)) {
    return(stats::setNames(rep(NA_real_, length(peaks)), peaks))
  }
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )
  state_cells <- density_result$cells[[state_id]]
  state_cells <- intersect(state_cells, colnames(object))
  if (length(state_cells) == 0) {
    return(stats::setNames(rep(0, length(peaks)), peaks))
  }

  peak_mat <- get_assay_matrix_safe(
    object,
    assay = peak_assay,
    layer_primary = "data",
    layer_fallback = "counts"
  )
  peaks <- intersect(peaks, rownames(peak_mat))
  if (length(peaks) == 0) {
    return(numeric(0))
  }

  vals <- Matrix::rowMeans(peak_mat[peaks, state_cells, drop = FALSE])
  vals <- as.numeric(vals)
  names(vals) <- peaks
  vals
}

get_dynamic_feature_genes <- function(
  object,
  pseudotime_column = NULL,
  padjust_threshold = 0.05
) {
  pseudotime_column <- pseudotime_column %ss%
    tryCatch(
      resolve_pseudotime_column(object, pseudotime_column),
      error = function(e) NULL
    )
  if (is.null(pseudotime_column)) {
    return(character(0))
  }
  tool_name <- paste0("DynamicFeatures_", pseudotime_column)
  tool <- tryCatch(object@tools[[tool_name]], error = function(e) NULL)
  df <- tool$DynamicFeatures %ss% NULL
  if (is.null(df) || nrow(df) == 0) {
    return(character(0))
  }
  padj <- if ("padjust" %in% colnames(df)) df$padjust else df$pvalue
  rownames(df)[is.finite(padj) & padj < padjust_threshold]
}

build_peak_tf_incidence <- function(object, tfs = NULL) {
  regions <- NetworkRegions(object)
  motif2tf <- NetworkTFs(object)
  if (is.null(regions@motifs) || is.null(motif2tf)) {
    return(NULL)
  }

  peaks2motif <- NULL
  if (is.list(regions@motifs)) {
    motif_parts <- lapply(regions@motifs, function(x) {
      if (is.null(x) || is.null(x$motifs)) {
        return(NULL)
      }
      x$motifs@data
    })
    motif_parts <- motif_parts[!vapply(motif_parts, is.null, logical(1))]
    if (length(motif_parts) == 0) {
      return(NULL)
    }
    all_peaks <- unique(unlist(
      lapply(motif_parts, rownames),
      use.names = FALSE
    ))
    all_motifs <- unique(unlist(
      lapply(motif_parts, colnames),
      use.names = FALSE
    ))
    peaks2motif <- Matrix::Matrix(
      0,
      nrow = length(all_peaks),
      ncol = length(all_motifs),
      sparse = TRUE
    )
    rownames(peaks2motif) <- all_peaks
    colnames(peaks2motif) <- all_motifs
    for (mat in motif_parts) {
      peaks2motif[rownames(mat), colnames(mat)] <- pmax(
        peaks2motif[rownames(mat), colnames(mat)],
        mat
      )
    }
  } else {
    peaks2motif <- regions@motifs@data
  }
  if (
    is.null(peaks2motif) || nrow(peaks2motif) == 0 || ncol(peaks2motif) == 0
  ) {
    return(NULL)
  }

  tf_use <- intersect(colnames(motif2tf), tfs %ss% colnames(motif2tf))
  motif_use <- intersect(colnames(peaks2motif), rownames(motif2tf))
  if (length(tf_use) == 0 || length(motif_use) == 0) {
    return(NULL)
  }
  peak_tf <- peaks2motif[, motif_use, drop = FALSE] %*%
    motif2tf[motif_use, tf_use, drop = FALSE]
  peak_tf@x[] <- 1
  peak_tf
}

compute_state_tf_motif_accessibility <- function(
  object,
  state_ids,
  tfs,
  pseudotime_column = NULL,
  peak_assay = NULL
) {
  peak_tf <- build_peak_tf_incidence(object, tfs = tfs)
  if (is.null(peak_tf)) {
    out <- Matrix::Matrix(
      0,
      nrow = length(tfs),
      ncol = length(state_ids),
      sparse = FALSE
    )
    rownames(out) <- tfs
    colnames(out) <- state_ids
    return(out)
  }

  peak_assay <- resolve_peak_assay_name(object, peak_assay = peak_assay)
  if (is.null(peak_assay)) {
    out <- Matrix::Matrix(
      0,
      nrow = length(tfs),
      ncol = length(state_ids),
      sparse = FALSE
    )
    rownames(out) <- tfs
    colnames(out) <- state_ids
    return(out)
  }

  peak_names <- rownames(peak_tf)
  state_peak_acc <- lapply(state_ids, function(sid) {
    acc <- compute_state_peak_accessibility(
      object,
      state_id = sid,
      peaks = peak_names,
      pseudotime_column = pseudotime_column,
      peak_assay = peak_assay
    )
    acc[peak_names] %ss% rep(0, length(peak_names))
  })
  names(state_peak_acc) <- state_ids

  support_mat <- vapply(
    state_ids,
    function(sid) {
      acc <- state_peak_acc[[sid]]
      names(acc) <- peak_names
      numer <- as.numeric(Matrix::crossprod(acc[peak_names], peak_tf))
      denom <- as.numeric(Matrix::colSums(peak_tf))
      out <- numer / pmax(denom, 1)
      names(out) <- colnames(peak_tf)
      out[tfs]
    },
    numeric(length(tfs))
  )
  rownames(support_mat) <- tfs
  support_mat
}

build_gene_peak_map <- function(
  object,
  genes,
  peak_assay = NULL,
  method = "Signac",
  upstream = 100000,
  downstream = 0,
  only_tss = FALSE
) {
  peak_assay <- resolve_peak_assay_name(object, peak_assay = peak_assay)
  if (is.null(peak_assay) || length(genes) == 0) {
    return(NULL)
  }
  gene_annot <- Signac::Annotation(.csn_get_assay(object, peak_assay))
  if (is.null(gene_annot)) {
    return(NULL)
  }
  gene_annot <- gene_annot[gene_annot$gene_name %in% genes]
  gene_annot <- gene_annot[!duplicated(gene_annot$gene_name)]
  if (length(gene_annot) == 0) {
    return(NULL)
  }
  peak_ranges <- Signac::StringToGRanges(rownames(.csn_get_assay(
    object,
    peak_assay
  )))
  peak_gene <- find_peaks_near_genes(
    peaks = peak_ranges,
    genes = gene_annot,
    method = method,
    upstream = upstream,
    downstream = downstream,
    only_tss = only_tss
  )
  peak_gene[, intersect(colnames(peak_gene), genes), drop = FALSE]
}

build_state_subnetwork_data <- function(
  object,
  state_id,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  top_tfs = 10,
  top_targets_per_tf = 12,
  max_targets = 60
) {
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = state_id,
    inferable_only = FALSE
  )
  df <- nets[[state_id]]
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  tf_counts <- tapply(df$target, df$regulator, function(x) length(unique(x)))
  tf_counts <- sort(tf_counts, decreasing = TRUE)
  tf_keep <- names(tf_counts)[seq_len(min(top_tfs, length(tf_counts)))]
  if (length(tf_keep) == 0) {
    return(NULL)
  }

  edge_list <- lapply(tf_keep, function(tf) {
    sub <- df[df$regulator == tf, , drop = FALSE]
    sub <- sub[order(abs(sub$weight), decreasing = TRUE), , drop = FALSE]
    sub <- sub[
      seq_len(min(max(top_targets_per_tf * 2, top_targets_per_tf), nrow(sub))), ,
      drop = FALSE
    ]
    sub
  })
  edge_df <- do.call(rbind, edge_list)
  if (is.null(edge_df) || nrow(edge_df) == 0) {
    return(NULL)
  }

  target_strength <- tapply(
    abs(edge_df$weight),
    edge_df$target,
    sum,
    na.rm = TRUE
  )
  target_shared_n <- tapply(edge_df$regulator, edge_df$target, function(x) {
    length(unique(x))
  })
  target_score <- target_strength * (1 + 0.75 * pmax(target_shared_n - 1, 0))
  target_summary <- data.frame(
    target = names(target_strength),
    strength = as.numeric(target_strength),
    shared_n = as.numeric(target_shared_n[names(target_strength)]),
    score = as.numeric(target_score[names(target_strength)]),
    stringsAsFactors = FALSE
  )
  target_summary <- target_summary[
    order(-target_summary$shared_n, -target_summary$score), ,
    drop = FALSE
  ]
  target_keep <- target_summary$target[seq_len(min(
    max_targets,
    nrow(target_summary)
  ))]
  edge_df <- edge_df[edge_df$target %in% target_keep, , drop = FALSE]
  if (nrow(edge_df) == 0) {
    return(NULL)
  }

  edge_df <- do.call(
    rbind,
    lapply(split(edge_df, edge_df$regulator), function(sub) {
      sub <- sub[order(abs(sub$weight), decreasing = TRUE), , drop = FALSE]
      sub[seq_len(min(top_targets_per_tf, nrow(sub))), , drop = FALSE]
    })
  )
  rownames(edge_df) <- NULL

  all_genes <- unique(c(edge_df$regulator, edge_df$target))
  expr <- compute_state_gene_expression(
    object,
    state_id = state_id,
    genes = all_genes,
    pseudotime_column = pseudotime_column,
    assay = assay
  )
  expr_z <- if (length(expr) > 1 && stats::sd(expr) > 0) {
    as.numeric(scale(expr))
  } else {
    rep(0, length(expr))
  }
  names(expr_z) <- names(expr)

  node_df <- data.frame(
    name = all_genes,
    node_type = ifelse(all_genes %in% tf_keep, "TF", "Target"),
    expression = expr[all_genes],
    expression_z = expr_z[all_genes],
    stringsAsFactors = FALSE
  )
  reg_sizes <- sort(tf_counts[tf_keep], decreasing = TRUE)
  node_df$node_size <- ifelse(
    node_df$node_type == "TF",
    as.numeric(reg_sizes[node_df$name]),
    as.numeric(target_strength[node_df$name] %ss% 1)
  )
  node_df$shared_n <- ifelse(
    node_df$node_type == "TF",
    NA_real_,
    as.numeric(target_shared_n[node_df$name] %ss% 1)
  )
  node_df$label <- ifelse(
    node_df$node_type == "TF" |
      (!is.na(node_df$shared_n) & node_df$shared_n >= 2) |
      node_df$node_size >=
        stats::quantile(node_df$node_size, probs = 0.93, na.rm = TRUE),
    node_df$name,
    NA_character_
  )

  list(edges = edge_df, nodes = node_df, tf_counts = tf_counts)
}

build_state_subnetwork_layout <- function(edge_df, node_df) {
  tf_nodes <- node_df$name[node_df$node_type == "TF"]
  target_nodes <- node_df$name[node_df$node_type == "Target"]

  if (length(tf_nodes) == 0) {
    node_df$x <- 0
    node_df$y <- 0
    return(node_df)
  }

  tf_angles <- seq(0, 2 * pi, length.out = length(tf_nodes) + 1)[
    -(length(tf_nodes) + 1)
  ]
  tf_angles <- tf_angles + pi / 2
  names(tf_angles) <- tf_nodes

  target_membership <- lapply(target_nodes, function(target) {
    sub <- edge_df[edge_df$target == target, , drop = FALSE]
    regs <- unique(sub$regulator)
    regs <- regs[regs %in% tf_nodes]
    reg_w <- abs(sub$weight[match(regs, sub$regulator)])
    reg_w[!is.finite(reg_w)] <- 1
    reg_w[reg_w <= 0] <- .Machine$double.eps
    mean_angle <- atan2(
      sum(sin(tf_angles[regs]) * reg_w),
      sum(cos(tf_angles[regs]) * reg_w)
    )
    data.frame(
      name = target,
      regulator_key = paste(sort(regs), collapse = "|"),
      primary_tf = if (length(regs) > 0) {
        regs[[which.max(reg_w)]]
      } else {
        NA_character_
      },
      shared_n = length(regs),
      mean_angle = mean_angle,
      stringsAsFactors = FALSE
    )
  })
  target_info <- do.call(rbind, target_membership)

  tf_target_n <- tapply(edge_df$target, edge_df$regulator, function(x) {
    length(unique(x))
  })
  tf_target_n <- as.numeric(tf_target_n[tf_nodes] %ss% rep(1, length(tf_nodes)))
  tf_target_n[!is.finite(tf_target_n)] <- 1
  tf_load <- tf_target_n / max(tf_target_n, na.rm = TRUE)
  tf_scale <- 0.72 + 0.22 * sqrt(tf_load)
  tf_df <- data.frame(
    name = tf_nodes,
    x = cos(tf_angles) * 1.75 * tf_scale,
    y = sin(tf_angles) * 1.05 * tf_scale,
    stringsAsFactors = FALSE
  )

  circular_delta <- function(a, b) {
    atan2(sin(a - b), cos(a - b))
  }

  place_arc_group <- function(
    nodes,
    center_angle,
    radius_x,
    radius_y,
    span = 0.22,
    rings = 3
  ) {
    if (length(nodes) == 0) {
      return(NULL)
    }
    rings <- max(1, rings)
    group_size <- ceiling(length(nodes) / rings)
    res <- vector("list", length(nodes))
    idx <- 1L
    for (ring_idx in seq_len(rings)) {
      start_idx <- idx
      end_idx <- min(length(nodes), start_idx + group_size - 1L)
      if (start_idx > length(nodes)) {
        break
      }
      ring_nodes <- nodes[start_idx:end_idx]
      local_span <- span * (1 + 0.1 * (ring_idx - 1))
      offsets <- if (length(ring_nodes) == 1L) {
        0
      } else {
        seq(-local_span, local_span, length.out = length(ring_nodes))
      }
      ring_x <- radius_x + 0.28 * (ring_idx - 1)
      ring_y <- radius_y + 0.16 * (ring_idx - 1)
      for (j in seq_along(ring_nodes)) {
        angle <- center_angle + offsets[[j]]
        res[[idx]] <- data.frame(
          name = ring_nodes[[j]],
          x = cos(angle) * ring_x,
          y = sin(angle) * ring_y,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
    do.call(rbind, res)
  }

  target_layout <- list()

  core_targets <- target_info[target_info$shared_n >= 3, , drop = FALSE]
  if (nrow(core_targets) > 0) {
    core_targets <- core_targets[
      order(core_targets$mean_angle, core_targets$name), ,
      drop = FALSE
    ]
    core_rings <- max(1L, min(4L, ceiling(nrow(core_targets) / 8)))
    core_group <- split(
      core_targets$name,
      rep(seq_len(core_rings), length.out = nrow(core_targets))
    )
    core_idx <- 1L
    for (ring_idx in seq_along(core_group)) {
      ring_nodes <- core_group[[ring_idx]]
      angles <- if (length(ring_nodes) == 1L) {
        0
      } else {
        seq(0, 2 * pi, length.out = length(ring_nodes) + 1)[
          -(length(ring_nodes) + 1)
        ] +
          0.25 * ring_idx
      }
      radius_x <- 0.58 + 0.34 * (ring_idx - 1)
      radius_y <- 0.38 + 0.18 * (ring_idx - 1)
      for (j in seq_along(ring_nodes)) {
        target_layout[[core_idx]] <- data.frame(
          name = ring_nodes[[j]],
          x = cos(angles[[j]]) * radius_x,
          y = sin(angles[[j]]) * radius_y,
          stringsAsFactors = FALSE
        )
        core_idx <- core_idx + 1L
      }
    }
  }

  bridge_targets <- target_info[target_info$shared_n == 2, , drop = FALSE]
  if (nrow(bridge_targets) > 0) {
    bridge_groups <- split(bridge_targets, bridge_targets$regulator_key)
    for (grp in bridge_groups) {
      regs <- strsplit(grp$regulator_key[[1]], "\\|", fixed = FALSE)[[1]]
      regs <- regs[regs %in% tf_nodes]
      if (length(regs) != 2) {
        next
      }
      center_angle <- atan2(
        sum(sin(tf_angles[regs])),
        sum(cos(tf_angles[regs]))
      )
      grp <- grp[
        order(abs(circular_delta(grp$mean_angle, center_angle)), grp$name), ,
        drop = FALSE
      ]
      placed <- place_arc_group(
        nodes = grp$name,
        center_angle = center_angle,
        radius_x = 1.45,
        radius_y = 0.95,
        span = 0.16,
        rings = max(1L, min(3L, ceiling(nrow(grp) / 5)))
      )
      if (!is.null(placed)) {
        target_layout <- c(target_layout, split(placed, seq_len(nrow(placed))))
      }
    }
  }

  unique_targets <- target_info[target_info$shared_n <= 1, , drop = FALSE]
  if (nrow(unique_targets) > 0) {
    unique_groups <- split(unique_targets, unique_targets$primary_tf)
    tf_order <- names(sort(tf_angles))
    tf_angles_ordered <- tf_angles[tf_order]
    prev_angles <- c(
      tail(tf_angles_ordered, 1) - 2 * pi,
      head(tf_angles_ordered, -1)
    )
    next_angles <- c(
      tail(tf_angles_ordered, -1),
      head(tf_angles_ordered, 1) + 2 * pi
    )
    tf_sector <- pmin(
      0.5,
      0.38 *
        pmin(
          abs(tf_angles_ordered - prev_angles),
          abs(next_angles - tf_angles_ordered)
        )
    )
    names(tf_sector) <- tf_order

    for (tf in names(unique_groups)) {
      grp <- unique_groups[[tf]]
      if (nrow(grp) == 0 || is.na(tf) || !tf %in% names(tf_angles)) {
        next
      }
      grp <- grp[order(grp$name), , drop = FALSE]
      placed <- place_arc_group(
        nodes = grp$name,
        center_angle = tf_angles[[tf]],
        radius_x = 1.68,
        radius_y = 1.00,
        span = 0.9 * (tf_sector[[tf]] %ss% 0.24),
        rings = max(1L, min(3L, ceiling(nrow(grp) / 7)))
      )
      if (!is.null(placed)) {
        target_layout <- c(target_layout, split(placed, seq_len(nrow(placed))))
      }
    }
  }

  placed_names <- vapply(target_layout, function(x) x$name[[1]], character(1))
  missing_targets <- setdiff(target_nodes, placed_names)
  if (length(missing_targets) > 0) {
    fallback_angles <- seq(0, 2 * pi, length.out = length(missing_targets) + 1)[
      -(length(missing_targets) + 1)
    ]
    fallback_df <- data.frame(
      name = missing_targets,
      x = cos(fallback_angles) * 1.2,
      y = sin(fallback_angles) * 0.75,
      stringsAsFactors = FALSE
    )
    target_layout <- c(
      target_layout,
      split(fallback_df, seq_len(nrow(fallback_df)))
    )
  }

  target_df <- if (length(target_layout) > 0) {
    do.call(rbind, target_layout)
  } else {
    data.frame(
      name = character(0),
      x = numeric(0),
      y = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  layout_df <- rbind(tf_df, target_df)
  out <- merge(node_df, layout_df, by = "name", all.x = TRUE, sort = FALSE)
  out[match(node_df$name, out$name), , drop = FALSE]
}

build_state_subnetwork_plot <- function(
  object,
  state_id,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  top_tfs = 10,
  top_targets_per_tf = 12,
  max_targets = 60,
  layout = "manual"
) {
  subdat <- build_state_subnetwork_data(
    object,
    state_id = state_id,
    network = network,
    pseudotime_column = pseudotime_column,
    assay = assay,
    top_tfs = top_tfs,
    top_targets_per_tf = top_targets_per_tf,
    max_targets = max_targets
  )
  if (is.null(subdat)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = paste("No network for", state_id),
          size = 5
        )
    )
  }

  layout_nodes <- build_state_subnetwork_layout(
    edge_df = subdat$edges,
    node_df = subdat$nodes
  )
  layout_nodes$plot_size <- ifelse(
    layout_nodes$node_type == "TF",
    layout_nodes$node_size * 0.46,
    layout_nodes$node_size * 1.22
  )

  edge_abs <- abs(subdat$edges$weight)
  edge_q <- stats::quantile(
    edge_abs,
    probs = c(0.1, 0.95),
    na.rm = TRUE,
    names = FALSE
  )
  edge_span <- max(edge_q[[2]] - edge_q[[1]], .Machine$double.eps)
  edge_scaled <- pmin(pmax((edge_abs - edge_q[[1]]) / edge_span, 0), 1)
  edge_scaled[!is.finite(edge_scaled)] <- 0.5
  edge_width_plot <- 0.16 + 0.58 * sqrt(edge_scaled)
  edge_alpha_plot <- 0.16 + 0.30 * edge_scaled

  graph <- igraph::graph_from_data_frame(
    d = transform(
      subdat$edges[, c("regulator", "target", "weight"), drop = FALSE],
      weight = pmax(abs(weight), .Machine$double.eps),
      edge_width_plot = edge_width_plot,
      edge_alpha_plot = edge_alpha_plot
    ),
    directed = TRUE,
    vertices = layout_nodes
  ) |>
    tidygraph::as_tbl_graph()

  if (identical(layout, "manual")) {
    p <- ggraph::ggraph(graph, layout = "manual", x = x, y = y)
  } else {
    p <- ggraph::ggraph(graph, layout = layout)
  }

  p <- p +
    ggraph::geom_edge_arc(
      ggplot2::aes(width = edge_width_plot, alpha = edge_alpha_plot),
      strength = 0.55,
      fold = TRUE,
      n = 80,
      color = "grey60",
      lineend = "round",
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      data = function(x) x[x$node_type == "Target", , drop = FALSE],
      ggplot2::aes(fill = expression_z, size = plot_size),
      shape = 21,
      color = "grey58",
      stroke = 0.22,
      alpha = 0.98
    ) +
    ggraph::geom_node_point(
      data = function(x) x[x$node_type == "TF", , drop = FALSE],
      ggplot2::aes(fill = expression_z, size = plot_size),
      shape = 23,
      color = "grey35",
      stroke = 0.28,
      alpha = 1
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      repel = TRUE,
      size = 3.25,
      family = "",
      color = "grey15"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#6BAED6",
      mid = "white",
      high = "#FB6A4A"
    ) +
    ggplot2::scale_size_continuous(range = c(3.0, 10.2)) +
    ggraph::scale_edge_width_identity() +
    ggraph::scale_edge_alpha_identity() +
    ggplot2::theme_void() +
    ggplot2::labs(title = state_id, fill = "Expr z", size = "Node size") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(4, 4, 4, 4),
      legend.position = "none"
    )

  p
}

build_state_tf_barplot <- function(
  object,
  state_id,
  network = DefaultNetwork(object),
  top_tfs = 8
) {
  tf_counts_tbl <- get_state_regulator_target_counts(
    object,
    network = network,
    state_ids = state_id
  )
  tf_counts_tbl <- tf_counts_tbl[
    tf_counts_tbl$state_id == state_id, ,
    drop = FALSE
  ]
  if (nrow(tf_counts_tbl) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = paste("No TF counts for", state_id),
          size = 4
        )
    )
  }

  tf_counts_tbl <- tf_counts_tbl[
    seq_len(min(top_tfs, nrow(tf_counts_tbl))), ,
    drop = FALSE
  ]
  tf_counts_tbl$regulator <- factor(
    tf_counts_tbl$regulator,
    levels = rev(tf_counts_tbl$regulator)
  )

  ggplot2::ggplot(
    tf_counts_tbl,
    ggplot2::aes(x = regulator, y = n_targets)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = regulator, y = 0, yend = n_targets),
      color = "#D9A05B"
    ) +
    ggplot2::geom_point(color = "#D98C2B", size = 2.2) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Number of targets") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8)
    )
}

#' @title Summarize per-state network architecture
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param celltypes Optional subset of states.
#' @param weight_cutoff Optional absolute-weight cutoff for exported edges.
#'
#' @return A data.frame summarizing cells, edges, regulators, and targets per state.
#' @export
summarize_state_networks <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  weight_cutoff = NULL
) {
  summarize_networks(
    object = object,
    network = network,
    pseudotime_column = pseudotime_column,
    celltypes = celltypes,
    weight_cutoff = weight_cutoff
  )
}

#' @title Compare rewiring between static cell type networks
#'
#' @param object A \code{Seurat} object.
#' @param network Static network name. Defaults to \code{DefaultNetwork(object)}.
#' @param celltypes Optional subset of cell types.
#' @param weight_column Weight column used when exporting networks.
#'
#' @return A list with pair summaries and per-TF rewiring tables.
#' @export
compare_celltype_rewiring <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  weight_column = "weight"
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "compare_celltype_rewiring: object must be a Seurat object.",
      call. = FALSE
    )
  }

  celltypes <- resolve_network_states(
    object,
    network = network,
    celltypes = celltypes,
    inferable_only = FALSE
  )
  if (length(celltypes) < 2) {
    stop(
      "compare_celltype_rewiring: need at least two cell types.",
      call. = FALSE
    )
  }

  nets <- export_state_networks(
    object,
    network = network,
    celltypes = celltypes,
    inferable_only = FALSE
  )

  tf_tables <- list()
  pair_summary <- list()
  idx <- 1L
  for (i in seq_len(length(celltypes) - 1)) {
    for (j in seq.int(i + 1, length(celltypes))) {
      celltype_from <- celltypes[[i]]
      celltype_to <- celltypes[[j]]
      a <- nets[[celltype_from]]
      b <- nets[[celltype_to]]
      if (is.null(a) || is.null(b) || nrow(a) == 0 || nrow(b) == 0) {
        next
      }
      if (!weight_column %in% colnames(a) || !weight_column %in% colnames(b)) {
        stop(
          "compare_celltype_rewiring: weight_column not found in exported networks.",
          call. = FALSE
        )
      }

      edges_a <- paste(a$regulator, a$target, sep = "\t")
      edges_b <- paste(b$regulator, b$target, sep = "\t")
      union_edges <- union(edges_a, edges_b)
      edge_jaccard <- if (length(union_edges) == 0) {
        1
      } else {
        length(intersect(edges_a, edges_b)) / length(union_edges)
      }

      tfs <- sort(unique(c(a$regulator, b$regulator)))
      shared_targets_global <- intersect(unique(a$target), unique(b$target))
      tf_df <- purrr::map_dfr(tfs, function(tf) {
        a_tf <- a[a$regulator == tf, c("target", weight_column), drop = FALSE]
        b_tf <- b[b$regulator == tf, c("target", weight_column), drop = FALSE]
        a_tf$target <- as.character(a_tf$target)
        b_tf$target <- as.character(b_tf$target)
        ta <- unique(a_tf$target)
        tb <- unique(b_tf$target)
        union_n <- length(union(ta, tb))
        ji <- if (union_n == 0) {
          NA_real_
        } else {
          length(intersect(ta, tb)) / union_n
        }
        wa <- aggregate_target_weights(a_tf, weight_column = weight_column)
        wb <- aggregate_target_weights(b_tf, weight_column = weight_column)
        compare_targets <- union(wa$target, wb$target)
        if (length(shared_targets_global) > 0) {
          compare_targets_shared <- intersect(
            compare_targets,
            shared_targets_global
          )
          if (length(compare_targets_shared) > 0) {
            compare_targets <- compare_targets_shared
          }
        }
        va <- stats::setNames(numeric(length(compare_targets)), compare_targets)
        vb <- stats::setNames(numeric(length(compare_targets)), compare_targets)
        if (nrow(wa) > 0 && length(compare_targets) > 0) {
          idx_a <- match(wa$target, compare_targets, nomatch = 0L)
          keep_a <- idx_a > 0
          va[idx_a[keep_a]] <- wa$x[keep_a]
        }
        if (nrow(wb) > 0 && length(compare_targets) > 0) {
          idx_b <- match(wb$target, compare_targets, nomatch = 0L)
          keep_b <- idx_b > 0
          vb[idx_b[keep_b]] <- wb$x[keep_b]
        }
        weighted_ji <- if (length(compare_targets) == 0) {
          NA_real_
        } else {
          denom <- sum(pmax(va, vb), na.rm = TRUE)
          if (denom <= 0) NA_real_ else sum(pmin(va, vb), na.rm = TRUE) / denom
        }
        data.frame(
          celltype_from = celltype_from,
          celltype_to = celltype_to,
          regulator = tf,
          n_targets_from = length(ta),
          n_targets_to = length(tb),
          regulon_jaccard = ji,
          weighted_regulon_jaccard = weighted_ji,
          rewiring_score = if (is.na(weighted_ji)) {
            NA_real_
          } else {
            1 - weighted_ji
          },
          binary_rewiring_score = if (is.na(ji)) NA_real_ else 1 - ji,
          stringsAsFactors = FALSE
        )
      })

      key <- paste(celltype_from, celltype_to, sep = "__")
      tf_tables[[key]] <- tf_df[
        order(tf_df$rewiring_score, decreasing = TRUE), ,
        drop = FALSE
      ]
      pair_summary[[idx]] <- data.frame(
        celltype_from = celltype_from,
        celltype_to = celltype_to,
        n_edges_from = nrow(a),
        n_edges_to = nrow(b),
        edge_jaccard = edge_jaccard,
        edge_rewiring_score = 1 - edge_jaccard,
        mean_tf_regulon_jaccard = mean(tf_df$regulon_jaccard, na.rm = TRUE),
        mean_tf_weighted_regulon_jaccard = mean(
          tf_df$weighted_regulon_jaccard,
          na.rm = TRUE
        ),
        mean_tf_rewiring_score = mean(tf_df$rewiring_score, na.rm = TRUE),
        mean_tf_binary_rewiring_score = mean(
          tf_df$binary_rewiring_score,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  list(
    celltype_pairs = do.call(rbind, pair_summary) %ss% data.frame(),
    tf_rewiring = tf_tables
  )
}

#' @title Plot cell type rewiring heatmap
#'
#' @param object A \code{Seurat} object.
#' @param network Static network name. Defaults to \code{DefaultNetwork(object)}.
#' @param celltypes Optional subset of cell types.
#' @param top_tfs Number of TFs to display.
#' @param metric Heatmap value. Supported values include
#'   \code{"rewiring_score"}, \code{"weighted_regulon_jaccard"},
#'   \code{"binary_rewiring_score"}, and \code{"regulon_jaccard"}.
#' @param weight_column Weight column used when exporting networks.
#'
#' @return A \code{ggplot} object.
#' @export
plot_celltype_specificity_heatmap <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  top_tfs = 20,
  metric = c("rewiring_score", "regulon_jaccard"),
  weight_column = "weight"
) {
  metric <- match.arg(metric)
  rewiring <- compare_celltype_rewiring(
    object,
    network = network,
    celltypes = celltypes,
    weight_column = weight_column
  )

  tf_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df$pair_id <- nm
    df$pair_label <- paste(df$celltype_from, df$celltype_to, sep = " vs ")
    df
  })
  if (nrow(tf_tbl) == 0) {
    stop(
      "plot_celltype_specificity_heatmap: no rewiring table available.",
      call. = FALSE
    )
  }

  tf_mean <- stats::aggregate(
    tf_tbl[[metric]],
    by = list(regulator = tf_tbl$regulator),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  tf_mean <- tf_mean[order(tf_mean$x, decreasing = TRUE), , drop = FALSE]
  top_regulators <- tf_mean$regulator[seq_len(min(top_tfs, nrow(tf_mean)))]

  plot_df <- tf_tbl[tf_tbl$regulator %in% top_regulators, , drop = FALSE]
  pair_levels <- unique(paste(
    rewiring$celltype_pairs$celltype_from,
    rewiring$celltype_pairs$celltype_to,
    sep = " vs "
  ))
  plot_df$pair_label <- factor(plot_df$pair_label, levels = pair_levels)
  plot_df$regulator <- factor(plot_df$regulator, levels = rev(top_regulators))

  low_color <- if (metric == "rewiring_score") "white" else "#08306B"
  high_color <- if (metric == "rewiring_score") "#B22222" else "white"

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pair_label, y = regulator, fill = !!rlang::sym(metric))
  ) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::scale_fill_gradient(
      low = low_color,
      high = high_color,
      na.value = "grey95"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Cell type pair", y = "TF", fill = metric)

  print(p)
  invisible(p)
}

#' @title Summarize transition-specific TF rewiring and target gain/loss
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param state_pairs Optional two-column data.frame with columns \code{state_from}
#'   and \code{state_to}. If \code{NULL}, adjacent ordered states are compared.
#' @param weight_column Weight column used when exporting networks.
#' @param top_tfs Number of rewired TFs kept per transition.
#' @param top_targets Number of gain/loss targets kept per TF and transition.
#'
#' @return A list with transition summary, TF summary, and target summary tables.
#' @export
summarize_transition_analysis <- function(
  object,
  network = DefaultNetwork(object),
  state_pairs = NULL,
  weight_column = "weight",
  top_tfs = 10,
  top_targets = 10
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "summarize_transition_analysis: object must be a Seurat object.",
      call. = FALSE
    )
  }

  rewiring <- compare_state_rewiring(
    object,
    network = network,
    state_pairs = state_pairs,
    weight_column = weight_column
  )
  tf_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df$pair_id <- nm
    df$transition_id <- paste(df$state_from, df$state_to, sep = " -> ")
    if (!is.null(top_tfs) && nrow(df) > top_tfs) {
      df <- df[seq_len(top_tfs), , drop = FALSE]
    }
    df
  })

  nets <- export_state_networks(
    object,
    network = network,
    inferable_only = FALSE
  )
  target_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    state_from <- unique(df$state_from)[1]
    state_to <- unique(df$state_to)[1]
    net_from <- nets[[state_from]]
    net_to <- nets[[state_to]]
    if (
      is.null(net_from) ||
        is.null(net_to) ||
        nrow(net_from) == 0 ||
        nrow(net_to) == 0
    ) {
      return(NULL)
    }
    if (
      !weight_column %in% colnames(net_from) ||
        !weight_column %in% colnames(net_to)
    ) {
      stop(
        "summarize_transition_analysis: weight_column not found in exported networks.",
        call. = FALSE
      )
    }
    tf_keep <- df$regulator[seq_len(min(top_tfs, nrow(df)))]

    gain_loss_tbl <- purrr::map_dfr(tf_keep, function(tf) {
      from_tbl <- net_from[
        net_from$regulator == tf,
        c("target", weight_column),
        drop = FALSE
      ]
      to_tbl <- net_to[
        net_to$regulator == tf,
        c("target", weight_column),
        drop = FALSE
      ]
      colnames(from_tbl)[2] <- "weight_from"
      colnames(to_tbl)[2] <- "weight_to"
      merged <- merge(from_tbl, to_tbl, by = "target", all = TRUE)
      merged$weight_from[is.na(merged$weight_from)] <- 0
      merged$weight_to[is.na(merged$weight_to)] <- 0
      merged$abs_from <- abs(as.numeric(merged$weight_from))
      merged$abs_to <- abs(as.numeric(merged$weight_to))
      merged$delta_abs_weight <- merged$abs_to - merged$abs_from
      merged$transition_direction <- ifelse(
        merged$abs_to > merged$abs_from,
        "gain",
        ifelse(merged$abs_to < merged$abs_from, "loss", "stable")
      )
      merged <- merged[merged$transition_direction != "stable", , drop = FALSE]
      if (nrow(merged) == 0) {
        return(NULL)
      }
      merged <- merged[
        order(abs(merged$delta_abs_weight), decreasing = TRUE), ,
        drop = FALSE
      ]
      if (!is.null(top_targets) && nrow(merged) > top_targets) {
        merged <- merged[seq_len(top_targets), , drop = FALSE]
      }
      merged$regulator <- tf
      merged$state_from <- state_from
      merged$state_to <- state_to
      merged$transition_id <- paste(state_from, state_to, sep = " -> ")
      merged
    })
    gain_loss_tbl
  })

  list(
    transition_summary = rewiring$state_pairs,
    tf_summary = tf_tbl,
    target_summary = target_tbl
  )
}

.summarize_states_base <- function(
  object,
  pseudotime_column = NULL,
  network = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop("summarize_states: object must be a Seurat object.", call. = FALSE)
  }
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )
  windows <- state_windows_from_result(
    object,
    density_result,
    pseudotime_column
  )
  metrics <- density_result$metrics
  windows$n_cells <- metrics$n_cells[match(windows$state_id, metrics$state_id)]
  windows$support_score <- metrics$support_score[match(
    windows$state_id,
    metrics$state_id
  )]
  windows$order_key <- vapply(
    as.character(windows$state_id),
    state_order_key,
    numeric(1)
  )
  windows <- windows[order(windows$order_key), , drop = FALSE]
  rownames(windows) <- NULL
  windows
}

#' @title Summarize per-state network architecture
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param celltypes Optional subset of states.
#' @param weight_cutoff Optional absolute-weight cutoff for exported edges.
#'
#' @return A data.frame summarizing cells, edges, regulators, and targets per state.
#' @export
summarize_networks <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  weight_cutoff = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop("summarize_networks: object must be a Seurat object.", call. = FALSE)
  }

  states <- .summarize_states_base(
    object,
    pseudotime_column = pseudotime_column,
    network = network
  )
  if (!is.null(celltypes)) {
    states <- states[states$state_id %in% celltypes, , drop = FALSE]
  }

  states$n_edges <- NA_integer_
  states$n_regulators <- NA_integer_
  states$n_targets <- NA_integer_
  states$has_network <- FALSE

  nets <- export_state_networks(
    object,
    network = network,
    celltypes = states$state_id,
    weight_cutoff = weight_cutoff,
    inferable_only = FALSE
  )

  for (sid in names(nets)) {
    df <- nets[[sid]]
    idx <- match(sid, states$state_id)
    if (is.na(idx) || is.null(df) || nrow(df) == 0) {
      next
    }
    states$n_edges[[idx]] <- nrow(df)
    states$n_regulators[[idx]] <- length(unique(df$regulator))
    states$n_targets[[idx]] <- length(unique(df$target))
    states$has_network[[idx]] <- TRUE
  }

  states
}

#' @title Rank network features across ordered states
#'
#' @param object A \code{Seurat} object.
#' @param type Feature type, \code{"regulators"} or \code{"targets"}.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param celltypes Optional subset of states.
#' @param method Gene-rank method, \code{"page_rank"} or \code{"degree_distribution"}.
#' @param top_n Number of rows to keep per state. Use \code{NULL} for all.
#' @param weight_column Weight column used when \code{type = "targets"}.
#'
#' @return A data.frame combining ranked features for the selected states.
#' @export
rank_features <- function(
  object,
  type = c("regulators", "targets"),
  network = DefaultNetwork(object),
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  top_n = 20,
  weight_column = "weight"
) {
  if (!methods::is(object, "Seurat")) {
    stop("rank_features: object must be a Seurat object.", call. = FALSE)
  }

  type <- match.arg(type)
  method <- match.arg(method)
  network_kind <- infer_feature_network_kind(
    object,
    network = network,
    celltypes = celltypes
  )

  if (identical(network_kind, "dynamic") && identical(type, "regulators")) {
    celltypes <- resolve_network_states(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = TRUE
    )
    if (length(celltypes) == 0) {
      return(data.frame())
    }

    object <- calculate_gene_rank(
      object,
      network = network,
      celltypes = celltypes,
      method = method
    )
    ranks <- GeneRanks(object, network = network, celltypes = celltypes)
    if (is.null(ranks) || length(ranks) == 0) {
      return(data.frame())
    }

    out <- purrr::imap_dfr(ranks, function(rank_obj, sid) {
      if (is.null(rank_obj) || is.null(rank_obj$df)) {
        return(NULL)
      }
      df <- rank_obj$df
      reg_col <- if ("regulator" %in% colnames(df)) {
        "regulator"
      } else if ("is_regulator" %in% colnames(df)) {
        "is_regulator"
      } else {
        NULL
      }
      if (!is.null(reg_col)) {
        keep_reg <- normalize_regulator_flag(df[[reg_col]])
        df <- df[keep_reg, , drop = FALSE]
      }
      if (nrow(df) == 0) {
        return(NULL)
      }
      if (!is.null(top_n) && nrow(df) > top_n) {
        df <- df[seq_len(top_n), , drop = FALSE]
      }
      df$state_id <- sid
      df$rank_method <- rank_obj$method %ss% method
      df$feature <- df$gene
      df$type <- type
      df$network_kind <- network_kind
      df$feature_score <- df[[choose_rank_value_column(
        df,
        preferred = c(
          "page_rank",
          "rank_value",
          "degree",
          "out_degree",
          "in_degree"
        )
      )]]
      df
    })
  } else if (identical(network_kind, "dynamic")) {
    nets <- export_state_networks(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = TRUE
    )
    if (length(nets) == 0) {
      return(data.frame())
    }

    out <- purrr::imap_dfr(nets, function(df, sid) {
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      if (!weight_column %in% colnames(df)) {
        stop(
          "rank_features: weight_column not found in exported network.",
          call. = FALSE
        )
      }

      abs_weight <- abs(df[[weight_column]])
      weighted_indegree <- tapply(abs_weight, df$target, sum, na.rm = TRUE)
      mean_abs_weight <- tapply(abs_weight, df$target, mean, na.rm = TRUE)
      max_abs_weight <- tapply(abs_weight, df$target, max, na.rm = TRUE)
      n_regulators <- tapply(abs_weight, df$target, length)

      agg <- data.frame(
        target = names(weighted_indegree),
        weighted_indegree = as.numeric(weighted_indegree),
        mean_abs_weight = as.numeric(mean_abs_weight[names(weighted_indegree)]),
        max_abs_weight = as.numeric(max_abs_weight[names(weighted_indegree)]),
        n_regulators = as.numeric(n_regulators[names(weighted_indegree)]),
        stringsAsFactors = FALSE
      )
      agg <- agg[
        order(agg$weighted_indegree, agg$n_regulators, decreasing = TRUE), ,
        drop = FALSE
      ]
      if (!is.null(top_n) && nrow(agg) > top_n) {
        agg <- agg[seq_len(top_n), , drop = FALSE]
      }
      agg$state_id <- sid
      agg$feature <- agg$target
      agg$type <- type
      agg$network_kind <- network_kind
      agg$feature_score <- agg$weighted_indegree
      agg
    })
  } else if (identical(type, "regulators")) {
    celltypes <- resolve_network_states(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = FALSE
    )
    if (length(celltypes) == 0) {
      return(data.frame())
    }

    object <- calculate_gene_rank(
      object,
      network = network,
      celltypes = celltypes,
      method = method
    )
    ranks <- GeneRanks(object, network = network, celltypes = celltypes)
    nets <- export_state_networks(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = FALSE
    )
    if (is.null(ranks) || length(ranks) == 0 || length(nets) == 0) {
      return(data.frame())
    }

    edge_tbl <- purrr::imap_dfr(nets, function(df, celltype_id) {
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      if (!weight_column %in% colnames(df)) {
        stop(
          "rank_features: weight_column not found in exported network.",
          call. = FALSE
        )
      }
      data.frame(
        celltype_id = celltype_id,
        regulator = as.character(df$regulator),
        target = as.character(df$target),
        abs_weight = abs(as.numeric(df[[weight_column]])),
        stringsAsFactors = FALSE
      )
    })
    if (nrow(edge_tbl) == 0) {
      return(data.frame())
    }

    out <- purrr::imap_dfr(ranks, function(rank_obj, celltype_id) {
      if (is.null(rank_obj) || is.null(rank_obj$df)) {
        return(NULL)
      }
      df <- rank_obj$df
      reg_col <- if ("regulator" %in% colnames(df)) {
        "regulator"
      } else if ("is_regulator" %in% colnames(df)) {
        "is_regulator"
      } else {
        NULL
      }
      if (!is.null(reg_col)) {
        keep_reg <- normalize_regulator_flag(df[[reg_col]])
        df <- df[keep_reg, , drop = FALSE]
      }
      if (nrow(df) == 0) {
        return(NULL)
      }

      score_col <- choose_rank_value_column(
        df,
        preferred = c(
          "page_rank",
          "rank_value",
          "degree",
          "out_degree",
          "in_degree"
        )
      )
      if (is.null(score_col)) {
        return(NULL)
      }

      sub_edges <- edge_tbl[edge_tbl$celltype_id == celltype_id, , drop = FALSE]
      target_n <- tapply(sub_edges$abs_weight, sub_edges$regulator, length)
      weight_sum <- tapply(
        sub_edges$abs_weight,
        sub_edges$regulator,
        sum,
        na.rm = TRUE
      )

      out_df <- data.frame(
        gene = as.character(df$gene),
        centrality_score = as.numeric(df[[score_col]]),
        n_targets = as.numeric(target_n[as.character(df$gene)]),
        outgoing_weight = as.numeric(weight_sum[as.character(df$gene)]),
        celltype_id = celltype_id,
        rank_method = rank_obj$method %ss% method,
        stringsAsFactors = FALSE
      )
      out_df$n_targets[is.na(out_df$n_targets)] <- 0
      out_df$outgoing_weight[is.na(out_df$outgoing_weight)] <- 0
      out_df
    })
    if (nrow(out) == 0) {
      return(out)
    }

    centrality_mat <- xtabs(centrality_score ~ gene + celltype_id, data = out)
    target_mat <- xtabs(n_targets ~ gene + celltype_id, data = out)
    weight_mat <- xtabs(outgoing_weight ~ gene + celltype_id, data = out)
    genes <- rownames(centrality_mat)
    spec_tbl <- data.frame(
      gene = genes,
      centrality_specificity = vapply(
        genes,
        function(g) specificity_index(centrality_mat[g, ]),
        numeric(1)
      ),
      target_specificity = vapply(
        genes,
        function(g) specificity_index(target_mat[g, ]),
        numeric(1)
      ),
      weight_specificity = vapply(
        genes,
        function(g) specificity_index(weight_mat[g, ]),
        numeric(1)
      ),
      anchor_celltype = colnames(centrality_mat)[max.col(
        centrality_mat[genes, , drop = FALSE],
        ties.method = "first"
      )],
      stringsAsFactors = FALSE
    )
    out <- merge(out, spec_tbl, by = "gene", all.x = TRUE)
    out$edge_specificity <- out$outgoing_weight /
      pmax(
        rowSums(weight_mat)[match(out$gene, rownames(weight_mat))],
        .Machine$double.eps
      )
    out <- split(out, out$celltype_id)
    out <- lapply(out, function(df) {
      df$centrality_score_z <- safe_zscore(df$centrality_score)
      df$n_targets_z <- safe_zscore(df$n_targets)
      df$outgoing_weight_z <- safe_zscore(df$outgoing_weight)
      df$edge_specificity_z <- safe_zscore(df$edge_specificity)
      df$centrality_specificity_z <- safe_zscore(df$centrality_specificity)
      df$target_specificity_z <- safe_zscore(df$target_specificity)
      df$weight_specificity_z <- safe_zscore(df$weight_specificity)
      df$tf_specificity_score <-
        0.30 *
        df$centrality_score_z +
        0.15 * df$n_targets_z +
        0.15 * df$outgoing_weight_z +
        0.15 * df$edge_specificity_z +
        0.15 * df$centrality_specificity_z +
        0.10 * df$weight_specificity_z
      df <- df[df$anchor_celltype == df$celltype_id, , drop = FALSE]
      df <- df[
        order(df$tf_specificity_score, decreasing = TRUE), ,
        drop = FALSE
      ]
      if (!is.null(top_n) && nrow(df) > top_n) {
        df <- df[seq_len(top_n), , drop = FALSE]
      }
      df$type <- type
      df$network_kind <- network_kind
      df$feature <- df$gene
      df$feature_score <- df$tf_specificity_score
      df
    })
    out <- do.call(rbind, out)
  } else {
    celltypes <- resolve_network_states(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = FALSE
    )
    nets <- export_state_networks(
      object,
      network = network,
      celltypes = celltypes,
      inferable_only = FALSE
    )
    if (length(nets) == 0) {
      return(data.frame())
    }

    tf_tbl <- rank_features(
      object,
      type = "regulators",
      network = network,
      celltypes = celltypes,
      method = method,
      top_n = NULL,
      weight_column = weight_column
    )
    tf_key <- if (nrow(tf_tbl) > 0) {
      paste(tf_tbl$celltype_id, tf_tbl$gene, sep = "\t")
    } else {
      character(0)
    }

    out <- purrr::imap_dfr(nets, function(df, celltype_id) {
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      if (!weight_column %in% colnames(df)) {
        stop(
          "rank_features: weight_column not found in exported network.",
          call. = FALSE
        )
      }
      df$abs_weight <- abs(as.numeric(df[[weight_column]]))
      df$is_specific_tf <- paste(celltype_id, df$regulator, sep = "\t") %in%
        tf_key

      weighted_indegree <- tapply(df$abs_weight, df$target, sum, na.rm = TRUE)
      mean_abs_weight <- tapply(df$abs_weight, df$target, mean, na.rm = TRUE)
      max_abs_weight <- tapply(df$abs_weight, df$target, max, na.rm = TRUE)
      n_regulators <- tapply(df$abs_weight, df$target, length)
      specific_tf_weight <- tapply(
        df$abs_weight * df$is_specific_tf,
        df$target,
        sum,
        na.rm = TRUE
      )
      specific_tf_n <- tapply(df$is_specific_tf, df$target, sum, na.rm = TRUE)

      agg <- data.frame(
        target = names(weighted_indegree),
        weighted_indegree = as.numeric(weighted_indegree),
        mean_abs_weight = as.numeric(mean_abs_weight[names(weighted_indegree)]),
        max_abs_weight = as.numeric(max_abs_weight[names(weighted_indegree)]),
        n_regulators = as.numeric(n_regulators[names(weighted_indegree)]),
        specific_tf_weight = as.numeric(specific_tf_weight[names(
          weighted_indegree
        )]),
        n_specific_tfs = as.numeric(specific_tf_n[names(weighted_indegree)]),
        celltype_id = celltype_id,
        stringsAsFactors = FALSE
      )
      agg$specific_tf_weight[is.na(agg$specific_tf_weight)] <- 0
      agg$n_specific_tfs[is.na(agg$n_specific_tfs)] <- 0
      agg
    })
    if (nrow(out) == 0) {
      return(out)
    }

    indegree_mat <- xtabs(weighted_indegree ~ target + celltype_id, data = out)
    specific_tf_mat <- xtabs(
      specific_tf_weight ~ target + celltype_id,
      data = out
    )
    targets <- rownames(indegree_mat)
    spec_tbl <- data.frame(
      target = targets,
      indegree_specificity = vapply(
        targets,
        function(g) specificity_index(indegree_mat[g, ]),
        numeric(1)
      ),
      specific_tf_support_specificity = vapply(
        targets,
        function(g) specificity_index(specific_tf_mat[g, ]),
        numeric(1)
      ),
      anchor_celltype = colnames(indegree_mat)[max.col(
        indegree_mat[targets, , drop = FALSE],
        ties.method = "first"
      )],
      stringsAsFactors = FALSE
    )
    out <- merge(out, spec_tbl, by = "target", all.x = TRUE)
    out$target_specificity <- out$weighted_indegree /
      pmax(
        rowSums(indegree_mat)[match(out$target, rownames(indegree_mat))],
        .Machine$double.eps
      )
    out <- split(out, out$celltype_id)
    out <- lapply(out, function(df) {
      df$weighted_indegree_z <- safe_zscore(df$weighted_indegree)
      df$n_regulators_z <- safe_zscore(df$n_regulators)
      df$specific_tf_weight_z <- safe_zscore(df$specific_tf_weight)
      df$n_specific_tfs_z <- safe_zscore(df$n_specific_tfs)
      df$target_specificity_z <- safe_zscore(df$target_specificity)
      df$indegree_specificity_z <- safe_zscore(df$indegree_specificity)
      df$specific_tf_support_specificity_z <- safe_zscore(
        df$specific_tf_support_specificity
      )
      df$target_specificity_score <-
        0.30 *
        df$weighted_indegree_z +
        0.10 * df$n_regulators_z +
        0.20 * df$specific_tf_weight_z +
        0.10 * df$n_specific_tfs_z +
        0.15 * df$target_specificity_z +
        0.15 * df$indegree_specificity_z
      df <- df[df$anchor_celltype == df$celltype_id, , drop = FALSE]
      df <- df[
        order(df$target_specificity_score, decreasing = TRUE), ,
        drop = FALSE
      ]
      if (!is.null(top_n) && nrow(df) > top_n) {
        df <- df[seq_len(top_n), , drop = FALSE]
      }
      df$type <- type
      df$network_kind <- network_kind
      df$feature <- df$target
      df$feature_score <- df$target_specificity_score
      df
    })
    out <- do.call(rbind, out)
  }

  if (nrow(out) == 0) {
    return(out)
  }
  out$order_key <- vapply(
    as.character(out$state_id),
    state_order_key,
    numeric(1)
  )
  out <- out[order(out$order_key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @title Compare adjacent state rewiring
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param state_pairs Optional two-column data.frame with columns \code{state_from}
#'   and \code{state_to}. If \code{NULL}, adjacent ordered states are compared.
#' @param weight_column Weight column used when exporting networks.
#'
#' @return A list with pair summaries and per-TF rewiring tables.
#' @export
compare_state_rewiring <- function(
  object,
  network = DefaultNetwork(object),
  state_pairs = NULL,
  weight_column = "weight"
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "compare_state_rewiring: object must be a Seurat object.",
      call. = FALSE
    )
  }

  states <- summarize_networks(object, network = network)
  states <- states[states$keep & states$inferable, , drop = FALSE]
  states <- states[
    states$state_id %in% resolve_network_states(object, network = network), ,
    drop = FALSE
  ]
  if (nrow(states) < 2) {
    stop(
      "compare_state_rewiring: need at least two kept inferable states.",
      call. = FALSE
    )
  }

  if (is.null(state_pairs)) {
    state_pairs <- data.frame(
      state_from = states$state_id[-nrow(states)],
      state_to = states$state_id[-1],
      stringsAsFactors = FALSE
    )
  }

  req_states <- unique(as.character(unlist(state_pairs, use.names = FALSE)))
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = req_states,
    weight_cutoff = NULL,
    inferable_only = FALSE
  )
  if (any(!req_states %in% names(nets))) {
    stop(
      "compare_state_rewiring: some requested states do not have exported networks.",
      call. = FALSE
    )
  }

  tf_tables <- list()
  pair_summary <- purrr::pmap_dfr(state_pairs, function(state_from, state_to) {
    a <- nets[[state_from]]
    b <- nets[[state_to]]
    if (is.null(a) || is.null(b) || nrow(a) == 0 || nrow(b) == 0) {
      return(NULL)
    }
    if (!weight_column %in% colnames(a) || !weight_column %in% colnames(b)) {
      stop(
        "compare_state_rewiring: weight_column not found in exported networks.",
        call. = FALSE
      )
    }

    edges_a <- paste(a$regulator, a$target, sep = "\t")
    edges_b <- paste(b$regulator, b$target, sep = "\t")
    union_edges <- union(edges_a, edges_b)
    edge_jaccard <- if (length(union_edges) == 0) {
      1
    } else {
      length(intersect(edges_a, edges_b)) / length(union_edges)
    }

    tfs <- sort(unique(c(a$regulator, b$regulator)))
    shared_targets_global <- intersect(unique(a$target), unique(b$target))
    tf_df <- purrr::map_dfr(tfs, function(tf) {
      a_tf <- a[a$regulator == tf, c("target", weight_column), drop = FALSE]
      b_tf <- b[b$regulator == tf, c("target", weight_column), drop = FALSE]
      a_tf$target <- as.character(a_tf$target)
      b_tf$target <- as.character(b_tf$target)
      ta <- unique(a_tf$target)
      tb <- unique(b_tf$target)
      union_n <- length(union(ta, tb))
      ji <- if (union_n == 0) NA_real_ else length(intersect(ta, tb)) / union_n

      wa <- aggregate_target_weights(a_tf, weight_column = weight_column)
      wb <- aggregate_target_weights(b_tf, weight_column = weight_column)
      compare_targets <- union(wa$target, wb$target)
      if (length(shared_targets_global) > 0) {
        compare_targets_shared <- intersect(
          compare_targets,
          shared_targets_global
        )
        if (length(compare_targets_shared) > 0) {
          compare_targets <- compare_targets_shared
        }
      }
      va <- stats::setNames(numeric(length(compare_targets)), compare_targets)
      vb <- stats::setNames(numeric(length(compare_targets)), compare_targets)
      if (nrow(wa) > 0 && length(compare_targets) > 0) {
        idx_a <- match(wa$target, compare_targets, nomatch = 0L)
        keep_a <- idx_a > 0
        va[idx_a[keep_a]] <- wa$x[keep_a]
      }
      if (nrow(wb) > 0 && length(compare_targets) > 0) {
        idx_b <- match(wb$target, compare_targets, nomatch = 0L)
        keep_b <- idx_b > 0
        vb[idx_b[keep_b]] <- wb$x[keep_b]
      }
      weighted_ji <- if (length(compare_targets) == 0) {
        NA_real_
      } else {
        denom <- sum(pmax(va, vb), na.rm = TRUE)
        if (denom <= 0) NA_real_ else sum(pmin(va, vb), na.rm = TRUE) / denom
      }
      data.frame(
        state_from = state_from,
        state_to = state_to,
        regulator = tf,
        n_targets_from = length(ta),
        n_targets_to = length(tb),
        regulon_jaccard = ji,
        weighted_regulon_jaccard = weighted_ji,
        rewiring_score = if (is.na(weighted_ji)) NA_real_ else 1 - weighted_ji,
        binary_rewiring_score = if (is.na(ji)) NA_real_ else 1 - ji,
        stringsAsFactors = FALSE
      )
    })

    tf_tables[[paste(state_from, state_to, sep = "__")]] <<-
      tf_df[order(tf_df$rewiring_score, decreasing = TRUE), , drop = FALSE]

    data.frame(
      state_from = state_from,
      state_to = state_to,
      n_edges_from = nrow(a),
      n_edges_to = nrow(b),
      edge_jaccard = edge_jaccard,
      edge_rewiring_score = 1 - edge_jaccard,
      mean_tf_regulon_jaccard = mean(tf_df$regulon_jaccard, na.rm = TRUE),
      mean_tf_weighted_regulon_jaccard = mean(
        tf_df$weighted_regulon_jaccard,
        na.rm = TRUE
      ),
      mean_tf_rewiring_score = mean(tf_df$rewiring_score, na.rm = TRUE),
      mean_tf_binary_rewiring_score = mean(
        tf_df$binary_rewiring_score,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })

  list(state_pairs = pair_summary, tf_rewiring = tf_tables)
}

#' @title Plot state-level network summary metrics
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param celltypes Optional subset of states.
#' @param metrics Metrics to display.
#' @param weight_cutoff Optional absolute-weight cutoff for exported edges.
#'
#' @return A \code{ggplot} object.
#' @export
plot_network_summary <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  metrics = c("n_cells", "n_edges", "n_regulators", "n_targets"),
  weight_cutoff = NULL
) {
  tbl <- summarize_networks(
    object,
    network = network,
    pseudotime_column = pseudotime_column,
    celltypes = celltypes,
    weight_cutoff = weight_cutoff
  )
  if (nrow(tbl) == 0) {
    stop("plot_network_summary: no states available.", call. = FALSE)
  }

  plot_df <- make_metric_long_table(tbl, metrics)
  plot_df <- plot_df[!is.na(plot_df$value), , drop = FALSE]
  state_levels <- unique(as.character(tbl$state_id))
  plot_df$state_id <- factor(
    as.character(plot_df$state_id),
    levels = state_levels
  )
  plot_df$metric <- factor(plot_df$metric, levels = metrics)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = state_id, y = value, fill = state_type)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::labs(x = "State", y = NULL, fill = "State type")

  print(p)
  invisible(p)
}

#' @title Plot top network features per state
#'
#' @param object A \code{Seurat} object.
#' @param type Feature type, \code{"regulators"} or \code{"targets"}.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param celltypes Optional subset of states.
#' @param method Gene-rank method used when \code{type = "regulators"}.
#' @param top_n Number of features to display per state.
#' @param value_column Optional score column to plot.
#' @param weight_column Weight column used when \code{type = "targets"}.
#' @param data Optional precomputed ranked-feature table used instead of \code{object}.
#' @param view View mode, \code{"ranking"} for per-state rankings or \code{"dynamics"}
#'   for per-state feature-score heatmaps (\code{type = "regulators"} only).
#' @param scale_value Score scaling for \code{view = "dynamics"}: \code{"zscore"}
#'   to standardize each feature across states, or \code{"raw"} for raw scores.
#'
#' @return A \code{ggplot} object.
#' @export
plot_features <- function(
  object = NULL,
  data = NULL,
  type = c("regulators", "targets"),
  network = NULL,
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  top_n = 10,
  value_column = NULL,
  weight_column = "weight",
  view = c("ranking", "dynamics"),
  scale_value = c("zscore", "raw")
) {
  type <- match.arg(type)
  method <- match.arg(method)
  view <- match.arg(view)
  scale_value <- match.arg(scale_value)
  if (is.null(data)) {
    if (is.null(object)) {
      stop("plot_features: provide either object or data.", call. = FALSE)
    }
    network <- network %ss% DefaultNetwork(object)
    tbl <- rank_features(
      object,
      type = type,
      network = network,
      celltypes = celltypes,
      method = method,
      top_n = if (view == "dynamics") NULL else top_n,
      weight_column = weight_column
    )
  } else {
    tbl <- data
  }
  if (nrow(tbl) == 0) {
    stop(sprintf("plot_features: no ranked %s available.", type), call. = FALSE)
  }

  network_kind <- unique(as.character(tbl$network_kind))
  network_kind <- network_kind[!is.na(network_kind) & nzchar(network_kind)]
  network_kind <- network_kind[[1]] %ss% "dynamic"

  if (identical(view, "dynamics")) {
    if (!identical(type, "regulators")) {
      stop(
        "plot_features: view='dynamics' only supports type='regulators'.",
        call. = FALSE
      )
    }
    if (!identical(network_kind, "dynamic")) {
      stop(
        "plot_features: view='dynamics' only supports dynamic networks.",
        call. = FALSE
      )
    }
    value_column <- value_column %ss%
      choose_rank_value_column(
        tbl,
        preferred = c(
          "feature_score",
          "page_rank",
          "rank_value",
          "degree",
          "out_degree",
          "in_degree"
        )
      )
    if (is.null(value_column)) {
      stop(
        "plot_features: unable to identify a score column for dynamics view.",
        call. = FALSE
      )
    }
    plot_df <- tbl[,
      c("state_id", "gene", value_column, "order_key"),
      drop = FALSE
    ]
    colnames(plot_df)[colnames(plot_df) == value_column] <- "score"
    if (scale_value == "zscore") {
      plot_df$score_plot <- ave(plot_df$score, plot_df$gene, FUN = function(x) {
        if (length(x) < 2 || stats::sd(x, na.rm = TRUE) == 0) {
          return(rep(0, length(x)))
        }
        as.numeric(scale(x))
      })
    } else {
      plot_df$score_plot <- plot_df$score
    }
    tf_rank <- stats::aggregate(
      abs(plot_df$score_plot),
      by = list(gene = plot_df$gene),
      FUN = function(x) max(x, na.rm = TRUE)
    )
    tf_rank <- tf_rank[order(tf_rank$x, decreasing = TRUE), , drop = FALSE]
    top_genes <- tf_rank$gene[seq_len(min(top_n, nrow(tf_rank)))]
    plot_df <- plot_df[plot_df$gene %in% top_genes, , drop = FALSE]
    state_levels <- unique(as.character(plot_df$state_id)[order(
      plot_df$order_key
    )])
    plot_df$state_id <- factor(plot_df$state_id, levels = unique(state_levels))
    plot_df$gene <- factor(plot_df$gene, levels = rev(top_genes))

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = state_id, y = gene, fill = score_plot)
    ) +
      ggplot2::geom_tile(color = "grey90") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = "State",
        y = "TF",
        fill = if (scale_value == "zscore") "z-score" else value_column
      )
    print(p)
    return(invisible(p))
  }

  preferred <- if (identical(type, "regulators")) {
    c(
      "feature_score",
      "tf_specificity_score",
      "page_rank",
      "rank_value",
      "degree",
      "out_degree",
      "in_degree"
    )
  } else {
    c(
      "feature_score",
      "target_specificity_score",
      "weighted_indegree",
      "mean_abs_weight",
      "max_abs_weight"
    )
  }
  feature_col <- if (identical(type, "regulators")) "gene" else "target"
  point_color <- if (identical(type, "regulators")) "#B22222" else "#1F78B4"

  value_column <- value_column %ss%
    choose_rank_value_column(tbl, preferred = preferred)
  if (is.null(value_column)) {
    stop(
      sprintf("plot_features: unable to identify a score column for %s.", type),
      call. = FALSE
    )
  }

  tbl <- prepare_feature_plot_table(
    tbl,
    feature_col = feature_col,
    value_col = value_column
  )
  facet_col <- if (identical(network_kind, "dynamic")) {
    "state_id"
  } else {
    "celltype_id"
  }
  tbl[[facet_col]] <- factor(
    as.character(tbl[[facet_col]]),
    levels = unique(as.character(tbl[[facet_col]]))
  )

  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = feature_label, y = !!rlang::sym(value_column))
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = feature_label, yend = !!rlang::sym(value_column)),
      y = 0,
      color = "grey70"
    ) +
    ggplot2::geom_point(color = point_color, size = 2) +
    ggplot2::facet_grid(
      stats::as.formula(paste(facet_col, "~ .")),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::labs(x = NULL, y = value_column)

  print(p)
  invisible(p)
}

#' @title Plot adjacent-state TF rewiring heatmap
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param type Plot type. \code{"rewiring"} shows TF-by-transition rewiring;
#'   \code{"networks"}, \code{"regulators"}, and \code{"targets"} show
#'   HARNexus-style pairwise state similarity based on edges, regulator sets,
#'   or target-gene sets.
#' @param plot_type Plot form. Default is \code{"heatmap"}; \code{"venn"} and
#'   \code{"upset"} are supported for similarity plots.
#' @param weight.by Weight column used when exporting networks.
#' @param group.by Reserved for future grouping support in similarity plots.
#' @param palette Palette name passed to
#'   \code{thisplot::palette_colors()} via \code{.multicsn_palette_colors()}.
#' @param palcolor Optional named color vector overriding state/group colors.
#'
#' @return A drawable plot object. When \pkg{ComplexHeatmap} is available,
#'   the function returns a grid grob captured from the drawn heatmap;
#'   otherwise it returns a \code{ggplot} object.
#' @export
plot_rewiring <- function(
  object,
  network = DefaultNetwork(object),
  type = c("rewiring", "networks", "regulators", "targets"),
  plot_type = c("heatmap", "venn", "upset"),
  weight.by = "weight",
  group.by = NULL,
  palette = "Chinese",
  palcolor = NULL
) {
  type <- match.arg(type)
  plot_type <- match.arg(plot_type)
  metric <- "rewiring_score"
  top_tfs <- 20L

  if (!identical(type, "rewiring")) {
    sim <- compare_state_similarity(
      object,
      network = network,
      type = type,
      weight_column = weight.by
    )
    sim_mat <- sim$similarity
    state_tbl <- sim$state_table
    feature_sets <- sim$feature_sets
    state_ids <- colnames(sim_mat)

    tf_counts <- state_tbl$n_regulators[match(state_ids, state_tbl$state_id)]
    gene_counts <- state_tbl$n_targets[match(state_ids, state_tbl$state_id)]
    cell_counts <- state_tbl$n_cells[match(state_ids, state_tbl$state_id)]

    state_colors <- character(0)
    if (!is.null(palcolor) && !is.null(names(palcolor))) {
      state_colors <- palcolor[state_ids]
    }
    if (length(state_colors) == 0 || any(is.na(state_colors))) {
      state_colors <- .multicsn_palette_colors(
        state_ids,
        palette = palette
      )
    }
    state_colors <- state_colors[state_ids]

    if (!identical(plot_type, "heatmap")) {
      if (!requireNamespace("thisplot", quietly = TRUE)) {
        stop(
          "plot_rewiring: plot_type='venn'/'upset' requires the thisplot package.",
          call. = FALSE
        )
      }
      all_features <- sort(unique(unlist(feature_sets, use.names = FALSE)))
      set_df <- as.data.frame(
        stats::setNames(
          lapply(state_ids, function(sid) all_features %in% feature_sets[[sid]]),
          state_ids
        ),
        stringsAsFactors = FALSE
      )
      rownames(set_df) <- all_features
      p <- thisplot::StatPlot(
        meta.data = set_df,
        stat.by = state_ids,
        plot_type = plot_type,
        palette = palette,
        palcolor = palcolor,
        title = switch(type,
          networks = "Overlap of network edge sets",
          regulators = "Overlap of regulator sets",
          targets = "Overlap of target-gene sets",
          "Set overlap"
        )
      )
      print(p)
      return(invisible(p))
    }

    if (
      requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("circlize", quietly = TRUE)
    ) {
      cell_counts[!is.finite(cell_counts)] <- 0
      names(cell_counts) <- state_ids
      tf_counts[!is.finite(tf_counts)] <- 0
      gene_counts[!is.finite(gene_counts)] <- 0
      state_fill_gp <- grid::gpar(fill = unname(state_colors), col = "black")
      positive_axis_at <- function(x) {
        xmax <- max(x, na.rm = TRUE)
        if (!is.finite(xmax) || xmax <= 0) {
          return(numeric(0))
        }
        brks <- pretty(c(0, xmax), n = 3)
        brks <- brks[brks > 0 & brks <= xmax]
        if (length(brks) == 0) {
          brks <- xmax
        }
        brks <- sort(unique(brks))
        brks
      }
      cells_axis <- positive_axis_at(cell_counts)
      tf_axis <- positive_axis_at(tf_counts)
      gene_axis <- positive_axis_at(gene_counts)

      top_anno <- ComplexHeatmap::columnAnnotation(
        "Cells count" = ComplexHeatmap::anno_barplot(
          cell_counts,
          baseline = 0,
          bar_width = 0.5,
          gp = state_fill_gp,
          border = TRUE,
          height = grid::unit(1.8, "cm"),
          axis_param = list(
            at = cells_axis,
            labels = format(cells_axis, scientific = FALSE, big.mark = ",")
          )
        ),
        " " = ComplexHeatmap::anno_simple(
          seq_along(state_ids),
          col = structure(
            unname(state_colors),
            names = as.character(seq_along(state_ids))
          ),
          height = grid::unit(0.15, "cm")
        ),
        height = grid::unit(2, "cm"),
        gap = grid::unit(0.2, "cm"),
        annotation_name_gp = grid::gpar(fontsize = 7)
      )

      right_anno <- ComplexHeatmap::rowAnnotation(
        "TFs count" = ComplexHeatmap::anno_barplot(
          tf_counts,
          baseline = 0,
          gp = state_fill_gp,
          border = TRUE,
          width = grid::unit(1.6, "cm"),
          axis_param = list(
            at = tf_axis,
            labels = format(tf_axis, scientific = FALSE, big.mark = ",")
          )
        ),
        "Target genes count" = ComplexHeatmap::anno_barplot(
          gene_counts,
          baseline = 0,
          gp = state_fill_gp,
          border = TRUE,
          width = grid::unit(1.6, "cm"),
          axis_param = list(
            at = gene_axis,
            labels = format(gene_axis, scientific = FALSE, big.mark = ",")
          )
        ),
        " " = ComplexHeatmap::anno_simple(
          seq_along(state_ids),
          col = structure(
            unname(state_colors),
            names = as.character(seq_along(state_ids))
          ),
          width = grid::unit(0.15, "cm")
        ),
        annotation_name_gp = grid::gpar(fontsize = 7)
      )

      col_fun <- circlize::colorRamp2(
        c(0, 0.25, 0.5, 0.75, 1),
        grDevices::colorRampPalette(c("white", "#2177B8"))(5)
      )
      heatmap_size <- grid::unit(max(3.2, length(state_ids) * 0.8), "cm")
      ht <- ComplexHeatmap::Heatmap(
        sim_mat,
        name = "Jaccard\nsimilarity",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        na_col = "white",
        show_row_names = TRUE,
        show_column_names = TRUE,
        rect_gp = grid::gpar(col = "gray50", lwd = 0.6),
        width = heatmap_size,
        height = heatmap_size,
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- sim_mat[i, j]
          if (is.na(val)) {
            grid::grid.rect(
              x = x,
              y = y,
              width = width,
              height = height,
              gp = grid::gpar(fill = "white", col = NA, lwd = 0)
            )
          } else {
            grid::grid.rect(
              x = x,
              y = y,
              width = width,
              height = height,
              gp = grid::gpar(fill = NA, col = "gray80", lwd = 0.5)
            )
            text_col <- if (val > 0.5) "white" else "gray40"
            grid::grid.text(
              sprintf("%.2f", val),
              x,
              y,
              gp = grid::gpar(fontsize = 8, col = text_col)
            )
          }
        },
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 45,
        top_annotation = top_anno,
        right_annotation = right_anno,
        heatmap_legend_param = list(
          title = "Jaccard\nsimilarity",
          title_gp = grid::gpar(fontsize = 10),
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
          legend_height = grid::unit(4, "cm")
        ),
        row_title = NULL
      )
      legend_list <- list()
      legend_list[[length(legend_list) + 1L]] <- ComplexHeatmap::Legend(
        title = "States",
        title_gp = grid::gpar(fontsize = 10),
        at = names(state_colors),
        labels = names(state_colors),
        legend_gp = grid::gpar(fill = unname(state_colors)),
        ncol = if (length(state_colors) > 9) 2 else 1
      )
      grob <- grid::grid.grabExpr(
        ComplexHeatmap::draw(
          ht,
          annotation_legend_list = legend_list,
          newpage = FALSE
        ),
        wrap = TRUE
      )
      grid::grid.newpage()
      grid::grid.draw(grob)
      attr(grob, "heatmap_object") <- ht
      return(invisible(grob))
    }

    sim_df <- as.data.frame(as.table(sim_mat), stringsAsFactors = FALSE)
    colnames(sim_df) <- c("state_from", "state_to", "similarity")
    p <- ggplot2::ggplot(
      sim_df,
      ggplot2::aes(x = state_to, y = state_from, fill = similarity)
    ) +
      ggplot2::geom_tile(color = "gray80", linewidth = 0.35) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf("%.2f", similarity),
          color = similarity > 0.5
        ),
        size = 3
      ) +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::colorRampPalette(c("white", "#2177B8"))(5),
        values = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0, 1)
      ) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "white", "FALSE" = "gray40"),
        guide = "none"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Jaccard\nsimilarity")
    print(p)
    return(invisible(p))
  }

  if (!identical(plot_type, "heatmap")) {
    stop(
      "plot_rewiring: plot_type='venn'/'upset' is only supported for type='networks', 'regulators', or 'targets'.",
      call. = FALSE
    )
  }

  rewiring <- compare_state_rewiring(
    object,
    network = network,
    weight_column = weight.by
  )

  tf_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df$pair_id <- nm
    df$pair_label <- paste(df$state_from, df$state_to, sep = " -> ")
    df
  })
  if (nrow(tf_tbl) == 0) {
    stop("plot_rewiring: no rewiring table available.", call. = FALSE)
  }

  tf_mean <- stats::aggregate(
    tf_tbl[[metric]],
    by = list(regulator = tf_tbl$regulator),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  tf_mean <- tf_mean[order(tf_mean$x, decreasing = TRUE), , drop = FALSE]
  top_regulators <- tf_mean$regulator[seq_len(min(top_tfs, nrow(tf_mean)))]

  plot_df <- tf_tbl[tf_tbl$regulator %in% top_regulators, , drop = FALSE]
  pair_levels <- unique(paste(
    rewiring$state_pairs$state_from,
    rewiring$state_pairs$state_to,
    sep = " -> "
  ))

  plot_df <- merge(
    expand.grid(
      regulator = top_regulators,
      pair_label = pair_levels,
      stringsAsFactors = FALSE
    ),
    plot_df,
    by = c("regulator", "pair_label"),
    all.x = TRUE,
    sort = FALSE
  )

  regulator_order <- stats::aggregate(
    plot_df[[metric]],
    by = list(regulator = plot_df$regulator),
    FUN = function(x) max(x, na.rm = TRUE)
  )
  regulator_order$x[!is.finite(regulator_order$x)] <- -Inf
  regulator_order <- regulator_order[
    order(regulator_order$x, decreasing = TRUE), ,
    drop = FALSE
  ]

  plot_df$pair_label <- factor(plot_df$pair_label, levels = pair_levels)
  plot_df$regulator <- factor(
    plot_df$regulator,
    levels = rev(regulator_order$regulator)
  )

  fill_title <- switch(metric,
    rewiring_score = "Rewiring score",
    weighted_regulon_jaccard = "Weighted Jaccard",
    binary_rewiring_score = "Binary rewiring",
    regulon_jaccard = "Binary Jaccard",
    metric
  )

  if (
    requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)
  ) {
    plot_df$pair_label_chr <- as.character(plot_df$pair_label)
    plot_df$regulator_chr <- as.character(plot_df$regulator)
    plot_df$metric_value <- plot_df[[metric]]
    plot_df$metric_missing <- is.na(plot_df$metric_value)
    heatmap_mat <- stats::xtabs(
      stats::as.formula("metric_value ~ regulator_chr + pair_label_chr"),
      data = plot_df
    )
    heatmap_mat <- heatmap_mat[
      rev(regulator_order$regulator),
      pair_levels,
      drop = FALSE
    ]
    missing_mask <- stats::xtabs(
      metric_missing ~ regulator_chr + pair_label_chr,
      data = plot_df
    )
    heatmap_mat[missing_mask > 0] <- NA_real_

    row_summary <- stats::aggregate(
      plot_df[[metric]],
      by = list(regulator = plot_df$regulator),
      FUN = function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) {
          return(0)
        }
        mean(x)
      }
    )
    row_summary$x[!is.finite(row_summary$x)] <- 0
    row_summary <- row_summary[
      match(rev(regulator_order$regulator), row_summary$regulator), ,
      drop = FALSE
    ]
    row_bar <- matrix(row_summary$x, ncol = 1)
    rownames(row_bar) <- row_summary$regulator

    pair_summary <- rewiring$state_pairs
    pair_summary$pair_label <- paste(
      pair_summary$state_from,
      pair_summary$state_to,
      sep = " -> "
    )
    pair_bar_col <- if (
      metric %in% c("rewiring_score", "binary_rewiring_score")
    ) {
      "mean_tf_rewiring_score"
    } else {
      "mean_tf_weighted_regulon_jaccard"
    }
    if (metric == "regulon_jaccard") {
      pair_bar_col <- "mean_tf_regulon_jaccard"
    }
    if (!pair_bar_col %in% colnames(pair_summary)) {
      pair_bar_col <- "edge_rewiring_score"
    }
    pair_summary <- pair_summary[
      match(pair_levels, pair_summary$pair_label), ,
      drop = FALSE
    ]
    pair_bar <- pair_summary[[pair_bar_col]]
    pair_bar[!is.finite(pair_bar)] <- 0
    names(pair_bar) <- pair_levels

    top_anno <- ComplexHeatmap::columnAnnotation(
      "Mean score" = ComplexHeatmap::anno_barplot(
        pair_bar,
        baseline = 0,
        bar_width = 0.5,
        gp = grid::gpar(fill = "#E59586", col = NA),
        border = FALSE,
        height = grid::unit(2, "cm"),
        axis_param = list(
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
        )
      ),
      annotation_name_gp = grid::gpar(fontsize = 7)
    )

    right_anno <- ComplexHeatmap::rowAnnotation(
      "Mean score" = ComplexHeatmap::anno_barplot(
        row_bar,
        baseline = 0,
        gp = grid::gpar(fill = "#E59586", col = NA),
        border = FALSE,
        width = grid::unit(2, "cm"),
        axis_param = list(
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
        )
      ),
      annotation_name_gp = grid::gpar(fontsize = 7)
    )

    col_fun <- circlize::colorRamp2(
      c(0, 0.25, 0.5, 0.75, 1),
      grDevices::colorRampPalette(c("white", "#E59586"))(5)
    )
    ht <- ComplexHeatmap::Heatmap(
      heatmap_mat,
      name = fill_title,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      col = col_fun,
      na_col = "white",
      rect_gp = grid::gpar(col = "gray50", lwd = 0.6),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- heatmap_mat[i, j]
        if (is.na(val)) {
          grid::grid.rect(
            x = x,
            y = y,
            width = width,
            height = height,
            gp = grid::gpar(fill = "white", col = NA, lwd = 0)
          )
        } else {
          grid::grid.rect(
            x = x,
            y = y,
            width = width,
            height = height,
            gp = grid::gpar(fill = NA, col = "gray80", lwd = 0.5)
          )
        }
      },
      top_annotation = top_anno,
      right_annotation = right_anno,
      heatmap_legend_param = list(
        title = fill_title,
        title_gp = grid::gpar(fontsize = 10),
        at = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
        legend_height = grid::unit(4, "cm")
      )
    )
    grob <- grid::grid.grabExpr(
      ComplexHeatmap::draw(ht, newpage = FALSE),
      wrap = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(grob)
    attr(grob, "heatmap_object") <- ht
    return(invisible(grob))
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pair_label, y = regulator, fill = !!rlang::sym(metric))
  ) +
    ggplot2::geom_tile(
      color = "#E6E6E6",
      linewidth = 0.35,
      width = 0.95,
      height = 0.95
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#E59586",
      na.value = "white",
      limits = c(0, 1)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "grey35",
        linewidth = 0.6
      ),
      legend.key.height = grid::unit(1.8, "cm"),
      legend.key.width = grid::unit(0.55, "cm")
    ) +
    ggplot2::labs(x = "Adjacent state pair", y = "TF", fill = fill_title)

  print(p)
  invisible(p)
}

#' @title Plot single-TF case study across ordered states
#'
#' @param object A \code{Seurat} object.
#' @param tf A transcription factor name.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param assay Assay used to fetch expression.
#' @param method Gene-rank method.
#' @param weight_column Weight column used when exporting networks.
#' @param top_targets Number of top targets to display in the target heatmap.
#' @param ranked_features Optional precomputed regulator ranking table; when
#'   \code{NULL}, \code{rank_features(type = "regulators")} is called.
#'
#' @return A patchwork/ggplot object.
#' @export
plot_tf_case_study <- function(
  object,
  tf,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  method = c("page_rank", "degree_distribution"),
  weight_column = "weight",
  top_targets = 15,
  ranked_features = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop("plot_tf_case_study: object must be a Seurat object.", call. = FALSE)
  }

  method <- match.arg(method)
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)

  reg_tbl <- ranked_features %ss%
    rank_features(
      object,
      type = "regulators",
      network = network,
      method = method,
      top_n = NULL
    )
  reg_tbl <- reg_tbl[reg_tbl$gene == tf, , drop = FALSE]
  if (nrow(reg_tbl) == 0) {
    stop(
      sprintf(
        "plot_tf_case_study: TF '%s' was not found in regulator rankings.",
        tf
      ),
      call. = FALSE
    )
  }

  score_col <- choose_rank_value_column(
    reg_tbl,
    preferred = c(
      "page_rank",
      "rank_value",
      "degree",
      "out_degree",
      "in_degree"
    )
  )
  if (is.null(score_col)) {
    stop(
      "plot_tf_case_study: unable to identify regulator score column.",
      call. = FALSE
    )
  }

  state_levels <- get_ordered_state_levels(reg_tbl$state_id)
  reg_tbl$state_id <- factor(reg_tbl$state_id, levels = state_levels)

  p_centrality <- ggplot2::ggplot(
    reg_tbl,
    ggplot2::aes(x = state_id, y = !!rlang::sym(score_col), group = 1)
  ) +
    ggplot2::geom_line(color = "#B22222", linewidth = 0.7) +
    ggplot2::geom_point(color = "#B22222", size = 2.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste0(tf, ": regulator activity"),
      x = "State",
      y = score_col
    )

  expr_df <- fetch_tf_expression_data(
    object,
    tf = tf,
    pseudotime_column = pseudotime_column,
    assay = assay
  )
  if (nrow(expr_df) == 0) {
    stop("plot_tf_case_study: no expression data available.", call. = FALSE)
  }
  if (all(is.na(expr_df$state_id))) {
    expr_df$state_id <- "unassigned"
  }
  expr_df$state_id <- factor(
    as.character(expr_df$state_id),
    levels = unique(c(state_levels, "unassigned"))
  )

  p_expr <- ggplot2::ggplot(
    expr_df,
    ggplot2::aes(x = pseudotime, y = expression, color = state_id)
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggplot2::geom_smooth(
      se = FALSE,
      linewidth = 0.8,
      color = "black",
      method = "loess"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0(tf, ": expression along pseudotime"),
      x = pseudotime_column,
      y = "Expression"
    ) +
    ggplot2::theme(legend.position = "bottom")

  rew <- compare_state_rewiring(
    object,
    network = network,
    weight_column = weight_column
  )
  rew_df <- purrr::imap_dfr(rew$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df <- df[df$regulator == tf, , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$pair_label <- paste(df$state_from, df$state_to, sep = " -> ")
    df
  })
  if (nrow(rew_df) == 0) {
    rew_df <- data.frame(
      pair_label = "NA",
      rewiring_score = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  p_rewire <- ggplot2::ggplot(
    rew_df,
    ggplot2::aes(x = pair_label, y = rewiring_score)
  ) +
    ggplot2::geom_col(fill = "#1F78B4", width = 0.7, na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste0(tf, ": adjacent-state rewiring"),
      x = "State pair",
      y = "Rewiring score"
    )

  target_tbl <- extract_tf_target_table(
    object,
    tf = tf,
    network = network,
    weight_column = weight_column
  )

  if (nrow(target_tbl) > 0) {
    target_rank <- stats::aggregate(
      target_tbl$abs_weight,
      by = list(target = target_tbl$target),
      FUN = function(x) sum(x, na.rm = TRUE)
    )
    target_rank <- target_rank[
      order(target_rank$x, decreasing = TRUE), ,
      drop = FALSE
    ]
    top_target_ids <- target_rank$target[seq_len(min(
      top_targets,
      nrow(target_rank)
    ))]
    target_plot_df <- target_tbl[
      target_tbl$target %in% top_target_ids, ,
      drop = FALSE
    ]
    target_plot_df$target <- factor(
      target_plot_df$target,
      levels = rev(top_target_ids)
    )
    target_plot_df$state_id <- factor(
      as.character(target_plot_df$state_id),
      levels = state_levels
    )

    p_targets <- ggplot2::ggplot(
      target_plot_df,
      ggplot2::aes(x = state_id, y = target, fill = abs_weight)
    ) +
      ggplot2::geom_tile(color = "grey90") +
      ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ) +
      ggplot2::labs(
        title = paste0(tf, ": target strength by state"),
        x = "State",
        y = "Target",
        fill = "Abs weight"
      )
  } else {
    p_targets <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0("No targets found for ", tf),
        size = 5
      )
  }

  p <- patchwork::wrap_plots(
    p_centrality,
    p_expr,
    p_rewire,
    p_targets,
    ncol = 2
  ) +
    patchwork::plot_annotation(title = paste0("TF case study: ", tf))

  print(p)
  invisible(p)
}

#' @title Plot paper-style ordered-state overview
#'
#' @param object A \code{Seurat} object.
#' @param network Dynamic network name. Defaults to \code{DefaultNetwork(object)}.
#' @param pseudotime_column Optional pseudotime column.
#' @param assay Assay used to summarize state-level expression in subnetworks.
#' @param state_ids Optional subset of states. Defaults to kept inferable states.
#' @param top_tfs_network Number of TFs shown in each state subnetwork.
#' @param top_targets_per_tf Number of top targets retained per TF in each subnetwork.
#' @param max_targets Maximum union of targets shown per state subnetwork.
#' @param top_tfs_bar Number of TFs shown in the bottom target-count summary.
#' @param palette Optional named palette for the density panel.
#'
#' @return A patchwork/ggplot object.
#' @export
plot_state_overview_figure <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  state_ids = NULL,
  top_tfs_network = 10,
  top_targets_per_tf = 12,
  max_targets = 60,
  top_tfs_bar = 8,
  palette = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "plot_state_overview_figure: object must be a Seurat object.",
      call. = FALSE
    )
  }

  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  states <- summarize_networks(
    object,
    pseudotime_column = pseudotime_column,
    network = network
  )
  states <- states[states$keep & states$inferable, , drop = FALSE]
  states <- states[
    states$state_id %in% resolve_network_states(object, network = network), ,
    drop = FALSE
  ]
  if (!is.null(state_ids)) {
    states <- states[states$state_id %in% state_ids, , drop = FALSE]
  }
  if (nrow(states) == 0) {
    stop(
      "plot_state_overview_figure: no inferable states available.",
      call. = FALSE
    )
  }
  states <- states[order(states$order_key), , drop = FALSE]

  p_density <- build_state_density_plot(
    object,
    pseudotime_column = pseudotime_column,
    palette = palette
  )

  subnetwork_plots <- lapply(states$state_id, function(sid) {
    build_state_subnetwork_plot(
      object,
      state_id = sid,
      network = network,
      pseudotime_column = pseudotime_column,
      assay = assay,
      top_tfs = top_tfs_network,
      top_targets_per_tf = top_targets_per_tf,
      max_targets = max_targets
    )
  })
  bar_plots <- lapply(states$state_id, function(sid) {
    build_state_tf_barplot(
      object,
      state_id = sid,
      network = network,
      top_tfs = top_tfs_bar
    )
  })
  names(subnetwork_plots) <- states$state_id
  names(bar_plots) <- states$state_id

  p_mid <- patchwork::wrap_plots(subnetwork_plots, nrow = 1)
  p_bottom <- patchwork::wrap_plots(bar_plots, nrow = 1)

  p <- p_density /
    p_mid /
    p_bottom +
    patchwork::plot_layout(heights = c(1.1, 2.3, 1.2)) +
    patchwork::plot_annotation(
      title = "Ordered-state GRN overview",
      subtitle = paste(states$state_id, collapse = "  ->  ")
    )

  print(p)
  invisible(p)
}

compute_pairwise_state_edge_overlap <- function(
  object,
  network = DefaultNetwork(object),
  state_ids = NULL
) {
  nets <- export_state_networks(
    object,
    network = network,
    celltypes = state_ids,
    inferable_only = TRUE
  )
  if (length(nets) < 2) {
    stop(
      "compute_pairwise_state_edge_overlap: need at least two inferable state networks.",
      call. = FALSE
    )
  }

  state_ids <- names(nets)
  state_ids <- state_ids[order(vapply(state_ids, state_order_key, numeric(1)))]
  edge_sets <- lapply(nets[state_ids], function(df) {
    if (is.null(df) || nrow(df) == 0) {
      character(0)
    } else {
      unique(paste(df$regulator, df$target, sep = "\t"))
    }
  })

  out <- lapply(state_ids, function(s1) {
    lapply(state_ids, function(s2) {
      e1 <- edge_sets[[s1]]
      e2 <- edge_sets[[s2]]
      union_edges <- union(e1, e2)
      inter_edges <- intersect(e1, e2)
      jaccard <- if (length(union_edges) == 0) {
        1
      } else {
        length(inter_edges) / length(union_edges)
      }
      data.frame(
        state_from = s1,
        state_to = s2,
        n_shared_edges = length(inter_edges),
        n_union_edges = length(union_edges),
        edge_jaccard = jaccard,
        stringsAsFactors = FALSE
      )
    })
  })
  out <- do.call(rbind, unlist(out, recursive = FALSE))
  out$state_from <- factor(out$state_from, levels = state_ids)
  out$state_to <- factor(out$state_to, levels = rev(state_ids))
  out
}

build_state_edge_overlap_plot <- function(
  object,
  network = DefaultNetwork(object),
  state_ids = NULL
) {
  overlap_df <- compute_pairwise_state_edge_overlap(
    object,
    network = network,
    state_ids = state_ids
  )
  ggplot2::ggplot(
    overlap_df,
    ggplot2::aes(x = state_from, y = state_to, fill = edge_jaccard)
  ) +
    ggplot2::geom_tile(color = "grey88") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", edge_jaccard)),
      size = 3.2,
      color = "grey10"
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#2166AC",
      limits = c(0, 1)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Pairwise edge overlap",
      x = NULL,
      y = NULL,
      fill = "Jaccard"
    )
}

build_adjacent_rewiring_summary_plot <- function(rewiring_result) {
  pair_df <- rewiring_result$state_pairs
  if (is.null(pair_df) || nrow(pair_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No adjacent-state rewiring summary",
          size = 5
        )
    )
  }

  pair_df$pair_label <- factor(
    paste(pair_df$state_from, pair_df$state_to, sep = " -> "),
    levels = paste(pair_df$state_from, pair_df$state_to, sep = " -> ")
  )
  long_df <- rbind(
    data.frame(
      pair_label = pair_df$pair_label,
      metric = "Edge rewiring",
      value = pair_df$edge_rewiring_score,
      stringsAsFactors = FALSE
    ),
    data.frame(
      pair_label = pair_df$pair_label,
      metric = "Mean TF rewiring",
      value = pair_df$mean_tf_rewiring_score,
      stringsAsFactors = FALSE
    )
  )
  long_df$metric <- factor(
    long_df$metric,
    levels = c("Edge rewiring", "Mean TF rewiring")
  )

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = pair_label, y = value, group = metric, color = metric)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(
      values = c("Edge rewiring" = "#B2182B", "Mean TF rewiring" = "#2166AC")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "top"
    ) +
    ggplot2::labs(
      title = "Adjacent-state rewiring summary",
      x = NULL,
      y = "Rewiring score",
      color = NULL
    )
}

build_rewiring_heatmap_plot <- function(
  object,
  network = DefaultNetwork(object),
  state_pairs = NULL,
  top_tfs = 20,
  metric = c("rewiring_score", "regulon_jaccard"),
  weight_column = "weight"
) {
  metric <- match.arg(metric)
  rewiring <- compare_state_rewiring(
    object,
    network = network,
    state_pairs = state_pairs,
    weight_column = weight_column
  )
  tf_tbl <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df$pair_id <- nm
    df$pair_label <- paste(df$state_from, df$state_to, sep = " -> ")
    df
  })
  if (nrow(tf_tbl) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No TF rewiring table available",
          size = 5
        )
    )
  }

  tf_mean <- stats::aggregate(
    tf_tbl[[metric]],
    by = list(regulator = tf_tbl$regulator),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  tf_mean <- tf_mean[order(tf_mean$x, decreasing = TRUE), , drop = FALSE]
  top_regulators <- tf_mean$regulator[seq_len(min(top_tfs, nrow(tf_mean)))]
  plot_df <- tf_tbl[tf_tbl$regulator %in% top_regulators, , drop = FALSE]
  pair_levels <- unique(paste(
    rewiring$state_pairs$state_from,
    rewiring$state_pairs$state_to,
    sep = " -> "
  ))
  plot_df$pair_label <- factor(plot_df$pair_label, levels = pair_levels)
  plot_df$regulator <- factor(plot_df$regulator, levels = rev(top_regulators))

  low_color <- if (metric == "rewiring_score") "white" else "#08306B"
  high_color <- if (metric == "rewiring_score") "#B22222" else "white"

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pair_label, y = regulator, fill = !!rlang::sym(metric))
  ) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::scale_fill_gradient(
      low = low_color,
      high = high_color,
      na.value = "grey95"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "TF rewiring across adjacent states",
      x = NULL,
      y = "TF",
      fill = metric
    )
}

build_tf_dynamics_heatmap_plot <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  top_tfs = 25,
  value_column = NULL,
  scale_value = c("zscore", "raw")
) {
  method <- match.arg(method)
  scale_value <- match.arg(scale_value)
  tbl <- rank_features(
    object,
    type = "regulators",
    network = network,
    celltypes = celltypes,
    method = method,
    top_n = NULL
  )
  if (nrow(tbl) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No regulator ranking available",
          size = 5
        )
    )
  }

  value_column <- value_column %ss%
    choose_rank_value_column(
      tbl,
      preferred = c(
        "page_rank",
        "rank_value",
        "degree",
        "out_degree",
        "in_degree"
      )
    )
  if (is.null(value_column)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No regulator score column",
          size = 5
        )
    )
  }

  plot_df <- tbl[,
    c("state_id", "gene", value_column, "order_key"),
    drop = FALSE
  ]
  colnames(plot_df)[colnames(plot_df) == value_column] <- "score"
  if (scale_value == "zscore") {
    plot_df$score_plot <- ave(plot_df$score, plot_df$gene, FUN = function(x) {
      if (length(x) < 2 || stats::sd(x, na.rm = TRUE) == 0) {
        return(rep(0, length(x)))
      }
      as.numeric(scale(x))
    })
  } else {
    plot_df$score_plot <- plot_df$score
  }

  tf_rank <- stats::aggregate(
    abs(plot_df$score_plot),
    by = list(gene = plot_df$gene),
    FUN = function(x) max(x, na.rm = TRUE)
  )
  tf_rank <- tf_rank[order(tf_rank$x, decreasing = TRUE), , drop = FALSE]
  top_genes <- tf_rank$gene[seq_len(min(top_tfs, nrow(tf_rank)))]
  plot_df <- plot_df[plot_df$gene %in% top_genes, , drop = FALSE]
  state_levels <- unique(as.character(plot_df$state_id)[order(
    plot_df$order_key
  )])
  plot_df$state_id <- factor(plot_df$state_id, levels = unique(state_levels))
  plot_df$gene <- factor(plot_df$gene, levels = rev(top_genes))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = state_id, y = gene, fill = score_plot)
  ) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "TF activity dynamics",
      x = "State",
      y = "TF",
      fill = if (scale_value == "zscore") "z-score" else value_column
    )
}

#' @title Plot manuscript Figure 2 rewiring overview
#' @inheritParams plot_rewiring
#' @param state_ids Optional subset of states.
#' @param top_tfs Number of rewired TFs kept per transition in the heatmap.
#' @param metric Rewiring metric used for the heatmap, \code{"rewiring_score"}
#'   or \code{"regulon_jaccard"}.
#' @param weight_column Weight column used when exporting networks.
#' @return A patchwork/ggplot object.
#' @export
plot_rewiring_overview_figure <- function(
  object,
  network = DefaultNetwork(object),
  state_ids = NULL,
  top_tfs = 20,
  metric = c("rewiring_score", "regulon_jaccard"),
  weight_column = "weight"
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "plot_rewiring_overview_figure: object must be a Seurat object.",
      call. = FALSE
    )
  }
  metric <- match.arg(metric)
  state_pairs_df <- NULL
  if (!is.null(state_ids)) {
    state_tbl <- summarize_networks(object, network = network)
    state_tbl <- state_tbl[state_tbl$keep & state_tbl$inferable, , drop = FALSE]
    state_tbl <- state_tbl[
      state_tbl$state_id %in% resolve_network_states(object, network = network), ,
      drop = FALSE
    ]
    state_tbl <- state_tbl[state_tbl$state_id %in% state_ids, , drop = FALSE]
    state_tbl <- state_tbl[order(state_tbl$order_key), , drop = FALSE]
    if (nrow(state_tbl) < 2) {
      stop(
        "plot_rewiring_overview_figure: need at least two states in state_ids.",
        call. = FALSE
      )
    }
    state_pairs_df <- data.frame(
      state_from = state_tbl$state_id[-nrow(state_tbl)],
      state_to = state_tbl$state_id[-1],
      stringsAsFactors = FALSE
    )
    state_ids <- state_tbl$state_id
  }
  rewiring <- compare_state_rewiring(
    object,
    network = network,
    state_pairs = state_pairs_df,
    weight_column = weight_column
  )
  p_overlap <- build_state_edge_overlap_plot(
    object,
    network = network,
    state_ids = state_ids
  )
  p_summary <- build_adjacent_rewiring_summary_plot(rewiring)
  p_heatmap <- build_rewiring_heatmap_plot(
    object,
    network = network,
    state_pairs = state_pairs_df,
    top_tfs = top_tfs,
    metric = metric,
    weight_column = weight_column
  )

  p <- (p_overlap | p_summary) /
    p_heatmap +
    patchwork::plot_layout(heights = c(1.0, 1.35)) +
    patchwork::plot_annotation(title = "Ordered-state GRN rewiring overview")
  print(p)
  invisible(p)
}

#' @title Plot manuscript Figure 3 key-TF summary
#' @inheritParams plot_tf_case_study
#' @param tfs Optional TF vector. If \code{NULL}, select automatically.
#' @param top_tfs_heatmap Number of TFs shown in the global TF-dynamics heatmap.
#' @param n_case_tfs Number of TFs selected automatically when \code{tfs=NULL}.
#' @param feature_scores Optional integrated key-TF score table; when \code{NULL},
#'   recomputed with \code{prioritize_features(strategy = "integrated")}.
#' @return A patchwork/ggplot object.
#' @export
plot_key_tf_summary_figure <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  tfs = NULL,
  method = c("page_rank", "degree_distribution"),
  weight_column = "weight",
  top_tfs_heatmap = 25,
  n_case_tfs = 6,
  feature_scores = NULL
) {
  if (!methods::is(object, "Seurat")) {
    stop(
      "plot_key_tf_summary_figure: object must be a Seurat object.",
      call. = FALSE
    )
  }
  method <- match.arg(method)
  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  assay <- assay %ss% Seurat::DefaultAssay(object)

  reg_tbl <- rank_features(
    object,
    type = "regulators",
    network = network,
    method = method,
    top_n = NULL
  )
  feature_scores <- feature_scores %ss%
    prioritize_features(
      object,
      type = "regulators",
      strategy = "integrated",
      network = network,
      pseudotime_column = pseudotime_column,
      assay = assay,
      method = method
    )
  rewiring <- compare_state_rewiring(
    object,
    network = network,
    weight_column = weight_column
  )
  if (is.null(tfs) || length(tfs) == 0) {
    tfs <- feature_scores$gene[seq_len(min(n_case_tfs, nrow(feature_scores)))]
  }
  tfs <- unique(as.character(tfs))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  if (length(tfs) == 0) {
    stop(
      "plot_key_tf_summary_figure: no TFs available for summary plotting.",
      call. = FALSE
    )
  }

  state_tbl <- summarize_networks(
    object,
    pseudotime_column = pseudotime_column,
    network = network
  )
  state_tbl <- state_tbl[state_tbl$keep & state_tbl$inferable, , drop = FALSE]
  state_tbl <- state_tbl[
    state_tbl$state_id %in% resolve_network_states(object, network = network), ,
    drop = FALSE
  ]
  state_tbl <- state_tbl[order(state_tbl$order_key), , drop = FALSE]
  state_levels <- state_tbl$state_id

  p_tf_heatmap <- build_tf_dynamics_heatmap_plot(
    object,
    network = network,
    method = method,
    top_tfs = top_tfs_heatmap
  )

  score_col <- choose_rank_value_column(
    reg_tbl,
    preferred = c(
      "page_rank",
      "rank_value",
      "degree",
      "out_degree",
      "in_degree"
    )
  )
  reg_sel <- reg_tbl[reg_tbl$gene %in% tfs, , drop = FALSE]
  reg_sel$state_id <- factor(
    as.character(reg_sel$state_id),
    levels = state_levels
  )
  reg_sel$gene <- factor(reg_sel$gene, levels = tfs)
  p_centrality <- ggplot2::ggplot(
    reg_sel,
    ggplot2::aes(
      x = state_id,
      y = !!rlang::sym(score_col),
      group = gene,
      color = gene
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = "Selected TF centrality trajectories",
      x = "State",
      y = score_col,
      color = "TF"
    )

  expr_mat <- vapply(
    state_levels,
    function(sid) {
      vals <- compute_state_gene_expression(
        object,
        state_id = sid,
        genes = tfs,
        pseudotime_column = pseudotime_column,
        assay = assay
      )
      out <- vals[tfs]
      out[is.na(out)] <- 0
      out
    },
    numeric(length(tfs))
  )
  rownames(expr_mat) <- tfs
  colnames(expr_mat) <- state_levels
  expr_plot <- as.data.frame(as.table(expr_mat), stringsAsFactors = FALSE)
  colnames(expr_plot) <- c("gene", "state_id", "expression")
  expr_plot$gene <- factor(expr_plot$gene, levels = rev(tfs))
  expr_plot$state_id <- factor(expr_plot$state_id, levels = state_levels)
  expr_plot$expression_z <- ave(
    expr_plot$expression,
    expr_plot$gene,
    FUN = function(x) {
      if (length(x) < 2 || stats::sd(x, na.rm = TRUE) == 0) {
        return(rep(0, length(x)))
      }
      as.numeric(scale(x))
    }
  )
  p_expression <- ggplot2::ggplot(
    expr_plot,
    ggplot2::aes(x = state_id, y = gene, fill = expression_z)
  ) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Selected TF expression by state",
      x = "State",
      y = "TF",
      fill = "Expr z"
    )

  rew_sel <- purrr::imap_dfr(rewiring$tf_rewiring, function(df, nm) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df <- df[df$regulator %in% tfs, , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$pair_label <- paste(df$state_from, df$state_to, sep = " -> ")
    df
  })
  if (nrow(rew_sel) > 0) {
    pair_levels <- unique(paste(
      rewiring$state_pairs$state_from,
      rewiring$state_pairs$state_to,
      sep = " -> "
    ))
    rew_sel$pair_label <- factor(rew_sel$pair_label, levels = pair_levels)
    rew_sel$regulator <- factor(rew_sel$regulator, levels = rev(tfs))
    p_rewiring_sel <- ggplot2::ggplot(
      rew_sel,
      ggplot2::aes(x = pair_label, y = regulator, fill = rewiring_score)
    ) +
      ggplot2::geom_tile(color = "grey90") +
      ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = "Selected TF rewiring",
        x = "Adjacent state pair",
        y = "TF",
        fill = "Rewiring"
      )
  } else {
    p_rewiring_sel <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "No selected-TF rewiring table",
        size = 5
      )
  }

  p <- (p_tf_heatmap | p_centrality) /
    (p_expression | p_rewiring_sel) +
    patchwork::plot_layout(widths = c(1.15, 1), heights = c(1.15, 1)) +
    patchwork::plot_annotation(
      title = "Key transcription-factor dynamics and rewiring",
      subtitle = paste("Selected TFs:", paste(tfs, collapse = ", "))
    )
  print(p)
  invisible(p)
}

#' @title Plot state-level network summary metrics
#'
#' @inheritParams plot_network_summary
#'
#' @return A \code{ggplot} object.
#' @export
plot_state_network_summary <- function(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  metrics = c("n_cells", "n_edges", "n_regulators", "n_targets"),
  weight_cutoff = NULL
) {
  plot_network_summary(
    object = object,
    network = network,
    pseudotime_column = pseudotime_column,
    celltypes = celltypes,
    metrics = metrics,
    weight_cutoff = weight_cutoff
  )
}

#' @title Plot adjacent-state TF rewiring heatmap
#'
#' @inheritParams plot_rewiring
#'
#' @return A \code{ggplot} object.
#' @export
plot_state_rewiring <- function(
  object,
  network = DefaultNetwork(object),
  type = c("rewiring", "networks", "regulators", "targets"),
  plot_type = c("heatmap", "venn", "upset"),
  weight.by = "weight",
  group.by = NULL,
  palette = "Chinese",
  palcolor = NULL
) {
  type <- match.arg(type)
  plot_type <- match.arg(plot_type)
  plot_rewiring(
    object = object,
    network = network,
    type = type,
    plot_type = plot_type,
    weight.by = weight.by,
    group.by = group.by,
    palette = palette,
    palcolor = palcolor
  )
}
