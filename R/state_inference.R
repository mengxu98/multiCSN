get_dynamic_genes <- function(
  seurat_obj,
  pseudotime_column,
  dynamic_features,
  genes_all,
  cores,
  assay = NULL,
  verbose = TRUE
) {
  if (is.null(dynamic_features) || !is.list(dynamic_features)) {
    return(genes_all)
  }
  n_candidates <- if (!is.null(dynamic_features$n_candidates)) {
    dynamic_features$n_candidates
  } else {
    NULL
  }
  padjust_threshold <- if (!is.null(dynamic_features$padjust_threshold)) {
    dynamic_features$padjust_threshold
  } else {
    0.05
  }
  fit_method <- if (!is.null(dynamic_features$fit_method)) {
    dynamic_features$fit_method
  } else {
    "pretsa"
  }

  tool_name <- paste0("DynamicFeatures_", pseudotime_column)
  tool <- seurat_obj@tools[[tool_name]]

  if (!is.null(tool) && !is.null(tool$DynamicFeatures)) {
    cached_features <- NULL
    if (!is.null(tool$raw_matrix) && ncol(tool$raw_matrix) > 1) {
      cached_features <- setdiff(colnames(tool$raw_matrix), "pseudotime")
    } else if (!is.null(tool$DynamicFeatures) && nrow(tool$DynamicFeatures) > 0) {
      cached_features <- rownames(tool$DynamicFeatures)
    }
    cached_features <- unique(as.character(cached_features))
    current_features <- unique(as.character(genes_all))
    same_candidates <- !is.null(cached_features) &&
      length(cached_features) == length(current_features) &&
      identical(sort(cached_features), sort(current_features))

    if (isTRUE(same_candidates)) {
      df <- tool$DynamicFeatures
      padj <- if ("padjust" %in% names(df)) df$padjust else df$pvalue
      dyn <- rownames(df)[padj < padjust_threshold]
      dyn <- intersect(dyn, genes_all)
      if (!is.null(n_candidates) && length(dyn) > n_candidates) {
        ord <- order(padj[dyn], decreasing = FALSE)
        dyn <- dyn[ord][seq_len(n_candidates)]
      }
      if (length(dyn) > 0) {
        thisutils::log_message(
          "Using {.val {length(dyn)}} pre-computed dynamic genes out of {.val {length(genes_all)}} candidates (padjust < {padjust_threshold}, {if (!is.null(n_candidates)) paste0('top ', n_candidates) else 'all'}, method={fit_method})",
          verbose = verbose
        )
        return(dyn)
      }
    } else {
      thisutils::log_message(
        "Dynamic gene cache ignored because candidate features changed ({.val {length(cached_features %ss% character(0))}} cached vs {.val {length(current_features)}} current). Recomputing dynamic features.",
        verbose = verbose
      )
    }
  }

  seurat_obj <- scop::RunDynamicFeatures(
    seurat_obj,
    lineages = pseudotime_column,
    features = genes_all,
    n_candidates = n_candidates,
    fit_method = fit_method,
    layer = "data",
    assay = assay,
    cores = cores,
    verbose = verbose
  )
  tool <- seurat_obj@tools[[tool_name]]
  df <- tool$DynamicFeatures
  padj <- if ("padjust" %in% names(df)) df$padjust else df$pvalue
  dyn <- rownames(df)[padj < padjust_threshold]
  dyn <- intersect(dyn, genes_all)
  if (!is.null(n_candidates) && length(dyn) > n_candidates) {
    ord <- order(padj[dyn], decreasing = FALSE)
    dyn <- dyn[ord][seq_len(n_candidates)]
  }
  if (length(dyn) == 0) {
    stop(
      sprintf(
        "No dynamic genes passed padjust_threshold = %s within %d candidates. Consider using a higher padjust_threshold.",
        padjust_threshold,
        length(genes_all)
      ),
      call. = FALSE
    )
  }
  thisutils::log_message(
    "{.val {length(dyn)}} dynamic genes out of {.val {length(genes_all)}} candidates (padjust < {padjust_threshold}, {if (!is.null(n_candidates)) paste0('top ', n_candidates) else 'all'}, method={fit_method})",
    verbose = verbose
  )
  list(genes = dyn, seurat = seurat_obj)
}

require_meta_column <- function(meta, col, what = "column") {
  if (is.null(col) || !nzchar(col)) {
    stop(sprintf("%s must be a non-empty string.", what), call. = FALSE)
  }
  if (!col %in% colnames(meta)) {
    stop(sprintf("%s '%s' not found in Seurat metadata.", what, col), call. = FALSE)
  }
  invisible(TRUE)
}

state_id_stable <- function(i) {
  paste0("state_", i)
}

state_id_transition <- function(i, j) {
  paste0("state_", i, "_", j)
}

sort_state_ids <- function(ids) {
  nums <- regmatches(ids, gregexpr("\\d+", ids))
  keys <- t(vapply(
    nums,
    function(n) {
      n <- as.integer(n)
      if (length(n) == 1) c(n, 0L) else c(n[1], n[2])
    },
    integer(2)
  ))
  ids[order(keys[, 1], keys[, 2])]
}

build_state_attributes <- function(
  cells,
  genes,
  peaks = NULL,
  window = NULL,
  state_type = NULL,
  parent_groups = NULL,
  inferable = NULL
) {
  if (is.null(peaks)) {
    peaks <- character(0)
  }
  list(
    genes = data.frame(gene = as.character(genes)),
    peaks = data.frame(peak = as.character(peaks)),
    cells = as.character(cells),
    window = as.numeric(window %ss% c(NA_real_, NA_real_)),
    state_type = state_type %ss% NA_character_,
    parent_groups = as.character(parent_groups %ss% NA_character_),
    inferable = inferable %ss% NA
  )
}

infer_state_skeleton <- function(meta, pseudotime_col, group_col = NULL) {
  require_meta_column(meta, pseudotime_col, what = "pseudotime_column")
  pt <- meta[[pseudotime_col]]
  keep <- !is.na(pt) & is.finite(pt)
  meta <- meta[keep, , drop = FALSE]
  pt <- meta[[pseudotime_col]]

  if (!is.null(group_col) && nzchar(group_col) && group_col %in% colnames(meta)) {
    grp <- meta[[group_col]]
  } else if (".group" %in% colnames(meta)) {
    grp <- meta[[".group"]]
  } else if ("celltype" %in% colnames(meta)) {
    grp <- meta[["celltype"]]
  } else if ("seurat_clusters" %in% colnames(meta)) {
    grp <- meta[["seurat_clusters"]]
  } else {
    grp <- rep("group_1", nrow(meta))
  }

  grp <- as.character(grp)
  df <- data.frame(group = grp, pseudotime = pt, stringsAsFactors = FALSE)
  df <- df[!is.na(df$pseudotime), , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No cells with non-NA pseudotime found.", call. = FALSE)
  }

  med <- stats::aggregate(pseudotime ~ group, df, stats::median)
  med <- med[order(med$pseudotime), , drop = FALSE]
  group_order <- as.character(med$group)

  list(
    meta = meta,
    group_col = group_col %ss% groups_col(meta),
    groups = group_order
  )
}

groups_col <- function(meta) {
  for (c in c("celltype", "seurat_clusters", "cluster", "group")) {
    if (c %in% colnames(meta)) {
      return(c)
    }
  }
  NULL
}

#' @title Infer dynamic state networks from density-partitioned pseudotime states
#'
#' @param object A \code{Seurat} object carrying CSNObject state.
#' @param pseudotime_column Metadata column containing pseudotime values.
#' @param method Network inference method.
#' @param penalty Penalty passed to \code{inferCSN()} for \code{Network}.
#' @param r_squared_threshold R-squared threshold.
#' @param regulators Optional regulator set.
#' @param targets Optional target set.
#' @param cores Number of cores.
#' @param verbose Logical.
#' @param dynamic_features Optional dynamic feature filtering config.
#' @param group_column Optional grouping column reused by \code{density_points()}.
#' @param posterior_threshold Maximum pairwise posterior dominance allowed inside a
#'   transition window; lower values make transition windows narrower.
#' @param n_grid Number of grid points used for kernel density evaluation.
#' @param n_boot Number of bootstrap resamples used to assess boundary stability.
#' @param bandwidth Bandwidth passed to \code{stats::density()}.
#' @param overwrite_density_points Logical; when \code{TRUE}, recompute and overwrite
#'   a cached density partition.
#' @param min_cells Ignored. Retained for backwards compatibility; transition
#'   states are determined automatically.
#' @param ... Additional parameters forwarded to \code{inferCSN()} on
#'   \code{Network} objects.
#'
#' @return A \code{Seurat} object with a dynamic network stored in CSNObject state.
#' @export
state_dynamic <- function(
  object,
  pseudotime_column,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  penalty = "L0",
  r_squared_threshold = 0,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  verbose = TRUE,
  dynamic_features = NULL,
  group_column = NULL,
  posterior_threshold = 0.8,
  n_grid = 512,
  n_boot = 50,
  bandwidth = "nrd0",
  overwrite_density_points = FALSE,
  min_cells = NULL,
  ...
) {
  method <- match.arg(method)

  seurat_obj <- .get_multicsn_seurat(object, require_state = TRUE)
  meta <- seurat_obj@meta.data
  require_meta_column(meta, pseudotime_column, what = "pseudotime_column")

  if (!is.null(min_cells) && isTRUE(verbose)) {
    thisutils::log_message(
      "Argument `min_cells` is ignored in density-based dynamic state inference; transition states are determined automatically.",
      verbose = verbose,
      message_type = "warning"
    )
  }

  density_result <- get_density_points_result(
    object,
    pseudotime_column = pseudotime_column,
    group_column = group_column,
    recompute = overwrite_density_points,
    posterior_threshold = posterior_threshold,
    n_grid = n_grid,
    n_boot = n_boot,
    bandwidth = bandwidth,
    verbose = verbose
  )


  rna_assay <- Params(object)$rna_assay %ss% "RNA"
  expr <- Seurat::GetAssayData(seurat_obj, assay = rna_assay, layer = "data")
  expr <- Matrix::as.matrix(expr)

  pt <- meta[[pseudotime_column]]
  keep_cells <- colnames(expr) %in% rownames(meta)[!is.na(pt) & is.finite(pt)]
  expr <- expr[, keep_cells, drop = FALSE]
  meta <- meta[colnames(expr), , drop = FALSE]

  if (verbose) {
    thisutils::log_message(
      "Dynamic state skeleton groups (ordered by median {pseudotime_column}): {.val {density_result$group_order}}",
      verbose = verbose
    )
  }

  windows_all <- density_result$windows
  metrics_all <- density_result$metrics
  if (verbose) {
    n_stable <- sum(windows_all$state_type == "stable")
    n_trans <- sum(windows_all$state_type == "transition" & windows_all$keep)
    thisutils::log_message(
      "Dynamic states: {.val {n_stable}} stable, {.val {n_trans}} transition states",
      verbose = verbose
    )
    for (sid in sort_state_ids(windows_all$state_id[windows_all$keep])) {
      row <- windows_all[match(sid, windows_all$state_id), , drop = FALSE]
      thisutils::log_message(
        "  state {sid} ({row$state_type[[1]]}): cells={length(density_result$cells[[sid]])}, {pseudotime_column}=[{format(row$left[[1]], digits = 3)}, {format(row$right[[1]], digits = 3)}], inferable={isTRUE(row$inferable[[1]])}",
        verbose = verbose
      )
    }
  }

  genes_all <- rownames(expr)
  targets_base <- if (is.null(targets)) genes_all else intersect(as.character(targets), genes_all)
  dyn_res <- get_dynamic_genes(
    seurat_obj = seurat_obj,
    pseudotime_column = pseudotime_column,
    dynamic_features = dynamic_features,
    genes_all = targets_base,
    cores = cores,
    assay = rna_assay,
    verbose = verbose
  )
  if (is.list(dyn_res)) {
    targets_use <- dyn_res$genes
    object <- dyn_res$seurat
    seurat_obj <- object
    expr <- Matrix::as.matrix(Seurat::GetAssayData(seurat_obj, assay = rna_assay, layer = "data"))
    expr <- expr[, colnames(expr) %in% rownames(meta), drop = FALSE]
  } else {
    targets_use <- dyn_res
  }
  targets_use <- intersect(targets_use, targets_base)
  if (length(targets_use) == 0) {
    targets_use <- targets_base
  }

  regulators_use <- regulators
  if (is.null(regulators_use)) {
    regulators_use <- .multicsn_get_tfs(object)
  }
  if (is.null(regulators_use)) {
    regulators_use <- character(0)
  }
  regulators_use <- intersect(as.character(regulators_use), genes_all)

  windows_keep <- windows_all[windows_all$keep & windows_all$inferable, , drop = FALSE]
  if (nrow(windows_keep) == 0) {
    stop("state_dynamic: no inferable states were produced by density_points().", call. = FALSE)
  }

  networks <- list()
  attrs <- .multicsn_get_attributes(object)

  for (sid in sort_state_ids(windows_keep$state_id)) {
    cells <- intersect(density_result$cells[[sid]], colnames(expr))
    if (length(cells) < 2) {
      next
    }
    row <- windows_keep[match(sid, windows_keep$state_id), , drop = FALSE]
    genes_use <- union(targets_use, regulators_use)
    mat_cells_genes <- t(expr[genes_use, cells, drop = FALSE])

    thisutils::log_message(
      "Inferring dynamic state network: {.val {sid}} ({row$state_type[[1]]}, {.val {length(cells)}} cells)",
      verbose = verbose
    )

    net <- initiate_object(
      mat_cells_genes,
      regulators = regulators_use,
      targets = targets_use
    )
    net <- inferCSN(
      net,
      penalty = penalty,
      r_squared_threshold = r_squared_threshold,
      cores = cores,
      verbose = verbose,
      method = method,
      ...
    )

    pt_s <- meta[cells, pseudotime_column]
    metric_row <- metrics_all[match(sid, metrics_all$state_id), , drop = FALSE]
    parent_groups <- strsplit(row$parent_groups[[1]], "\\|")[[1]]

    net@params$state_id <- sid
    net@params$state_type <- row$state_type[[1]]
    net@params$state_window <- c(row$left[[1]], row$right[[1]])
    net@params$pseudotime_range <- c(min(pt_s, na.rm = TRUE), max(pt_s, na.rm = TRUE))
    net@params$parent_groups <- parent_groups
    net@params$atac_supported <- FALSE
    net@params$support_score <- metric_row$support_score[[1]] %ss% NA_real_
    net@params$network_type <- "dynamic"
    net@params$pseudotime_column <- pseudotime_column
    net@params$state_strategy <- "density_intersection"
    net@params$inferable <- TRUE

    networks[[sid]] <- net
    attrs[[sid]] <- build_state_attributes(
      cells = cells,
      genes = targets_use,
      window = c(row$left[[1]], row$right[[1]]),
      state_type = row$state_type[[1]],
      parent_groups = parent_groups,
      inferable = TRUE
    )
  }

  dyn_name <- paste0(method, "_network")
  dyn_networks <- .multicsn_get_networks(object, network = dyn_name) %ss% list()
  for (sid in names(networks)) {
    dyn_networks[[sid]] <- networks[[sid]]
  }

  object <- .multicsn_set_networks(object, dyn_name, dyn_networks)
  object <- .multicsn_set_attributes(object, attrs)
  object
}
