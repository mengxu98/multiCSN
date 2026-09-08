#' @title Partition pseudotime states by density intersections
#'
#' @param object A \code{Seurat} object or a metadata \code{data.frame}.
#' @param pseudotime_column Column name containing pseudotime values.
#' @param group_column Optional metadata column containing ordered groups.
#'   When \code{NULL}, the function reuses \code{celltype},
#'   \code{seurat_clusters}, \code{cluster}, or \code{group} when available.
#' @param meta_data Deprecated alias of \code{object} for backwards compatibility.
#' @param store Logical. When \code{object} is a \code{Seurat} object, store the
#'   partition result in \code{object@@tools}.
#' @param overwrite Logical. Overwrite an existing stored partition.
#' @param posterior_threshold Maximum pairwise posterior dominance allowed inside
#'   a transition window. Lower values make transition windows narrower.
#' @param n_grid Number of grid points used for kernel density evaluation.
#' @param n_boot Number of bootstrap resamples used to assess boundary stability.
#' @param bandwidth Bandwidth passed to \code{stats::density()}.
#' @param plot Logical. When \code{TRUE}, call \code{plot_states()} after
#'   computing the partition.
#' @param verbose Logical.
#'
#' @return If \code{object} is a \code{Seurat} object, returns the updated
#'   object. Otherwise returns the partition result list.
#' @export
#'
#' @examples
#' test_data <- rbind(
#'   data.frame(cluster = "cluster1", pseudotime = rnorm(500, mean = 1, sd = 1)),
#'   data.frame(cluster = "cluster2", pseudotime = rnorm(500, mean = 2, sd = 1)),
#'   data.frame(cluster = "cluster3", pseudotime = rnorm(500, mean = 3, sd = 1))
#' )
#' density_points(test_data, pseudotime_column = "pseudotime", group_column = "cluster")
density_points <- function(
  object = NULL,
  pseudotime_column = "pseudotime",
  group_column = NULL,
  meta_data = NULL,
  store = TRUE,
  overwrite = FALSE,
  posterior_threshold = 0.8,
  n_grid = 512,
  n_boot = 50,
  bandwidth = "nrd0",
  plot = FALSE,
  verbose = TRUE
) {
  if (is.null(object)) {
    object <- meta_data
  }
  if (is.null(object)) {
    stop(
      "density_points: please supply a Seurat object or metadata data.frame.",
      call. = FALSE
    )
  }

  if (methods::is(object, "Seurat")) {
    tool_name <- density_points_tool_name(pseudotime_column)
    if (!overwrite && !is.null(object@tools[[tool_name]])) {
      if (isTRUE(plot)) {
        plot_states(object, pseudotime_column = pseudotime_column)
      }
      return(object)
    }
    prepared <- prepare_density_input(
      object = object,
      pseudotime_column = pseudotime_column,
      group_column = group_column
    )
    result <- infer_state_partition(
      meta = prepared$meta,
      pseudotime_column = pseudotime_column,
      group_column = prepared$group_column,
      posterior_threshold = posterior_threshold,
      n_grid = n_grid,
      n_boot = n_boot,
      bandwidth = bandwidth,
      verbose = verbose
    )
    if (isTRUE(store)) {
      object@tools[[tool_name]] <- result
    }
    if (isTRUE(plot)) {
      plot_states(object, pseudotime_column = pseudotime_column)
    }
    return(object)
  }

  if (!is.data.frame(object)) {
    stop(
      "density_points: object must be a Seurat object or metadata data.frame.",
      call. = FALSE
    )
  }

  prepared <- prepare_density_input(
    meta_data = object,
    pseudotime_column = pseudotime_column,
    group_column = group_column
  )
  infer_state_partition(
    meta = prepared$meta,
    pseudotime_column = pseudotime_column,
    group_column = prepared$group_column,
    posterior_threshold = posterior_threshold,
    n_grid = n_grid,
    n_boot = n_boot,
    bandwidth = bandwidth,
    verbose = verbose
  )
}

#' @title Plot pseudotime states
#'
#' @param object A \code{Seurat} object.
#' @param pseudotime_column Optional pseudotime column. If missing, the function
#'   reuses the value stored in the active dynamic network or in
#'   \code{DensityPoints_<pseudotime_column>}.
#' @param group_column Optional grouping column for density curves.
#' @param palette Optional colors.
#' @param palette_name Palette name used by \code{.multicsn_palette_colors()}.
#'
#' @return A data.frame of kept state windows.
#' @export
plot_states <- function(
  object,
  pseudotime_column = NULL,
  group_column = NULL,
  palette = NULL,
  palette_name = "Chinese"
) {
  if (!methods::is(object, "Seurat")) {
    stop("plot_states: object must be a Seurat object.", call. = FALSE)
  }

  pseudotime_column <- resolve_pseudotime_column(object, pseudotime_column)
  density_result <- get_density_points_result(
    object,
    pseudotime_column,
    recompute = FALSE
  )
  state_windows <- state_windows_from_result(
    object,
    density_result,
    pseudotime_column
  )
  state_windows <- state_windows[state_windows$keep, , drop = FALSE]
  if (nrow(state_windows) == 0) {
    stop("plot_states: no kept states available.", call. = FALSE)
  }

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
    stop("plot_states: no cells with valid group labels.", call. = FALSE)
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

  line_data <- state_windows
  line_data <- line_data[
    order(vapply(
      as.character(line_data$state_id),
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
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      data = density_result$boundaries,
      ggplot2::aes(xintercept = boundary),
      linetype = "dashed",
      color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = line_data,
      ggplot2::aes(x = left, xend = right, y = y, yend = y),
      color = "black",
      linewidth = 0.5,
      arrow = ggplot2::arrow(angle = 15, type = "closed"),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = line_data,
      ggplot2::aes(x = label_x, y = label_y, label = state_id),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(
      x = pseudotime_column,
      y = "Density",
      color = "cluster",
      fill = "cluster"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  print(p)
  state_windows
}

#' @title Legacy dynamic windowing wrapper
#'
#' @param meta_data Metadata data.frame.
#' @param pseudotime_column Pseudotime column.
#' @param group_column Group column.
#' @return Named list of state cells.
#' @export
dynamic_windowing <- function(
  meta_data,
  pseudotime_column = "pseudotime",
  group_column = "cluster"
) {
  res <- density_points(
    object = meta_data,
    pseudotime_column = pseudotime_column,
    group_column = group_column,
    store = FALSE,
    verbose = FALSE
  )
  res$cells[res$windows$state_id[res$windows$keep]]
}

density_points_tool_name <- function(pseudotime_column) {
  paste0("DensityPoints_", pseudotime_column)
}

resolve_pseudotime_column <- function(object, pseudotime_column = NULL) {
  if (!methods::is(object, "Seurat")) {
    return(pseudotime_column)
  }
  if (!is.null(pseudotime_column) && nzchar(pseudotime_column)) {
    return(pseudotime_column)
  }
  active_network <- tryCatch(
    .multicsn_resolve_network(
      object,
      network = NULL,
      preferred = "dynamic",
      verbose = FALSE,
      caller = "plot_states"
    ),
    error = function(e) NULL
  )
  if (!is.null(active_network) && nzchar(active_network)) {
    nets_active <- GetNetwork(object, network = active_network)
    first_net <- if (!is.null(nets_active) && length(nets_active) > 0) {
      nets_active[[1]]
    } else {
      NULL
    }
    if (!is.null(first_net) && methods::is(first_net, "Network")) {
      pt <- first_net@params$pseudotime_column %ss% NULL
      if (!is.null(pt) && nzchar(pt)) {
        return(pt)
      }
    }
  }
  tool_names <- names(object@tools)
  density_tools <- grep("^DensityPoints_", tool_names, value = TRUE)
  if (length(density_tools) > 0) {
    return(sub("^DensityPoints_", "", density_tools[[1]]))
  }
  stop("Cannot determine pseudotime_column.", call. = FALSE)
}

get_density_points_result <- function(
  object,
  pseudotime_column,
  group_column = NULL,
  recompute = FALSE,
  ...
) {
  tool_name <- density_points_tool_name(pseudotime_column)
  cached <- object@tools[[tool_name]]
  if (!recompute && !is.null(cached)) {
    return(cached)
  }
  object2 <- density_points(
    object = object,
    pseudotime_column = pseudotime_column,
    group_column = group_column,
    store = TRUE,
    overwrite = recompute,
    ...
  )
  object2@tools[[tool_name]]
}

prepare_density_input <- function(
  object = NULL,
  meta_data = NULL,
  pseudotime_column,
  group_column = NULL
) {
  if (!is.null(object) && methods::is(object, "Seurat")) {
    meta <- object@meta.data
    if (is.null(group_column) || !nzchar(group_column)) {
      group_column <- groups_col(meta)
      if (is.null(group_column)) {
        idents <- tryCatch(
          as.character(Seurat::Idents(object)),
          error = function(e) NULL
        )
        if (!is.null(idents)) {
          meta$.group <- idents
          group_column <- ".group"
        }
      }
    }
  } else {
    meta <- meta_data %ss% object
  }
  require_meta_column(meta, pseudotime_column, what = "pseudotime_column")
  if (is.null(group_column) || !nzchar(group_column)) {
    group_column <- groups_col(meta)
  }
  if (
    is.null(group_column) ||
      !nzchar(group_column) ||
      !group_column %in% colnames(meta)
  ) {
    stop("density_points: unable to determine group_column.", call. = FALSE)
  }

  meta <- meta[
    !is.na(meta[[pseudotime_column]]) & is.finite(meta[[pseudotime_column]]), ,
    drop = FALSE
  ]
  meta <- meta[
    !is.na(meta[[group_column]]) & as.character(meta[[group_column]]) != "", ,
    drop = FALSE
  ]
  if (nrow(meta) == 0) {
    stop(
      "density_points: no cells with finite pseudotime and valid groups.",
      call. = FALSE
    )
  }
  meta$.group <- as.character(meta[[group_column]])
  list(meta = meta, group_column = group_column)
}

infer_state_partition <- function(
  meta,
  pseudotime_column,
  group_column,
  posterior_threshold = 0.8,
  n_grid = 512,
  n_boot = 50,
  bandwidth = "nrd0",
  verbose = TRUE
) {
  sk <- infer_state_skeleton(
    meta,
    pseudotime_col = pseudotime_column,
    group_col = group_column
  )
  meta <- sk$meta
  groups <- sk$groups
  pt <- meta[[pseudotime_column]]
  names(pt) <- rownames(meta)
  global_min <- min(pt, na.rm = TRUE)
  global_max <- max(pt, na.rm = TRUE)
  if (length(groups) < 2) {
    stop(
      "density_points: at least two ordered groups are required.",
      call. = FALSE
    )
  }

  group_stats <- purrr::map_dfr(groups, function(g) {
    vals <- pt[meta$.group == g]
    dens <- density_on_grid(
      vals,
      from = global_min,
      to = global_max,
      n_grid = n_grid,
      bandwidth = bandwidth
    )
    data.frame(
      group = g,
      n_cells = length(vals),
      pseudotime_median = stats::median(vals),
      mode = dens$mode,
      t_min = min(vals),
      t_max = max(vals),
      stringsAsFactors = FALSE
    )
  })

  if (verbose) {
    thisutils::log_message(
      "Dynamic state skeleton groups (ordered by median {pseudotime_column}): {.val {groups}}",
      verbose = verbose
    )
  }

  pair_list <- vector("list", length(groups) - 1L)
  for (i in seq_len(length(groups) - 1L)) {
    left_group <- groups[[i]]
    right_group <- groups[[i + 1L]]
    left_values <- pt[meta$.group == left_group]
    right_values <- pt[meta$.group == right_group]
    pair_list[[i]] <- compute_pair_metrics(
      left_group = left_group,
      right_group = right_group,
      left_values = left_values,
      right_values = right_values,
      global_min = global_min,
      global_max = global_max,
      posterior_threshold = posterior_threshold,
      n_grid = n_grid,
      n_boot = n_boot,
      bandwidth = bandwidth
    )
  }

  boundaries <- do.call(rbind, lapply(pair_list, function(x) x$boundary_row))
  if (nrow(boundaries) > 1) {
    boundaries$boundary <- cummax(boundaries$boundary)
  }
  adjusted_intervals <- adjust_transition_intervals(
    boundaries,
    global_min,
    global_max
  )
  boundaries$left_transition <- adjusted_intervals$left
  boundaries$right_transition <- adjusted_intervals$right
  boundaries$keep <- adjusted_intervals$keep

  window_res <- build_windows_and_cells(
    meta = meta,
    groups = groups,
    pseudotime_column = pseudotime_column,
    boundaries = boundaries,
    global_min = global_min,
    global_max = global_max
  )

  metrics <- window_res$metrics
  stable_counts <- metrics$n_cells[match(
    state_id_stable(seq_along(groups)),
    metrics$state_id
  )]
  names(stable_counts) <- state_id_stable(seq_along(groups))

  for (i in seq_len(nrow(boundaries))) {
    tid <- state_id_transition(i, i + 1L)
    midx <- match(tid, metrics$state_id)
    if (is.na(midx)) {
      next
    }
    n_left <- metrics$n_left[[midx]]
    n_right <- metrics$n_right[[midx]]
    n_eff <- if ((n_left + n_right) > 0) {
      4 * n_left * n_right / (n_left + n_right)
    } else {
      0
    }
    ref <- min(
      stable_counts[[state_id_stable(i)]],
      stable_counts[[state_id_stable(i + 1L)]]
    )
    inferable <- boundaries$keep[[i]] &&
      is.finite(ref) &&
      ref > 0 &&
      (n_eff >= max(3, 0.08 * ref))
    boundaries$inferable[[i]] <- inferable
    metrics$inferable[[midx]] <- inferable
    metrics$support_score[[midx]] <- boundaries$support_score[[i]]
  }
  stable_idx <- metrics$state_type == "stable"
  metrics$inferable[stable_idx] <- metrics$n_cells[stable_idx] >= 10
  metrics$support_score[stable_idx] <- NA_real_

  windows <- window_res$windows
  windows$keep <- TRUE
  windows$inferable <- metrics$inferable[match(
    windows$state_id,
    metrics$state_id
  )]
  transition_idx <- windows$state_type == "transition"
  windows$keep[transition_idx] <- boundaries$keep[match(
    windows$state_id[transition_idx],
    state_id_transition(
      seq_len(nrow(boundaries)),
      seq_len(nrow(boundaries)) + 1L
    )
  )]
  windows$inferable[transition_idx] <- boundaries$inferable[match(
    windows$state_id[transition_idx],
    state_id_transition(
      seq_len(nrow(boundaries)),
      seq_len(nrow(boundaries)) + 1L
    )
  )]

  cells <- window_res$cells
  metrics$keep <- windows$keep[match(metrics$state_id, windows$state_id)]
  metrics$n_cells <- vapply(cells[metrics$state_id], length, integer(1))

  if (verbose) {
    n_stable <- sum(windows$state_type == "stable")
    n_trans_keep <- sum(windows$state_type == "transition" & windows$keep)
    thisutils::log_message(
      "Dynamic states: {.val {n_stable}} stable, {.val {n_trans_keep}} transition states",
      verbose = verbose
    )
    for (sid in sort_state_ids(windows$state_id[windows$keep])) {
      row <- windows[match(sid, windows$state_id), , drop = FALSE]
      thisutils::log_message(
        "  {sid} ({row$state_type[[1]]}): cells={length(cells[[sid]])}, {pseudotime_column}=[{format(row$left[[1]], digits = 3)}, {format(row$right[[1]], digits = 3)}]",
        verbose = verbose
      )
    }
  }

  list(
    params = list(
      pseudotime_column = pseudotime_column,
      group_column = group_column,
      strategy = "density_intersection",
      posterior_threshold = posterior_threshold,
      n_grid = n_grid,
      n_boot = n_boot,
      bandwidth = bandwidth
    ),
    group_order = groups,
    group_stats = group_stats,
    boundaries = boundaries,
    windows = windows,
    cells = cells,
    metrics = metrics
  )
}

compute_pair_metrics <- function(
  left_group,
  right_group,
  left_values,
  right_values,
  global_min,
  global_max,
  posterior_threshold,
  n_grid,
  n_boot,
  bandwidth
) {
  left_density <- density_on_grid(
    left_values,
    global_min,
    global_max,
    n_grid,
    bandwidth
  )
  right_density <- density_on_grid(
    right_values,
    global_min,
    global_max,
    n_grid,
    bandwidth
  )
  grid <- left_density$x
  mode_left <- left_density$mode
  mode_right <- right_density$mode
  boundary <- find_density_boundary(
    grid = grid,
    left_density = left_density$y,
    right_density = right_density$y,
    mode_left = mode_left,
    mode_right = mode_right
  )
  post_left <- left_density$y /
    pmax(left_density$y + right_density$y, .Machine$double.eps)
  interval <- find_transition_interval(
    grid = grid,
    posterior_left = post_left,
    boundary = boundary,
    threshold = posterior_threshold
  )

  overlap_score <- clip01(trapz_integral(
    grid,
    pmin(left_density$y, right_density$y)
  ))
  entropy_p <- pmin(pmax(post_left, .Machine$double.eps), 1 - .Machine$double.eps)
  entropy_vals <- -(entropy_p * log(entropy_p) + (1 - entropy_p) * log(1 - entropy_p)) / log(2)
  pair_support <- left_density$y + right_density$y
  between_modes <- grid >= min(mode_left, mode_right) &
    grid <= max(mode_left, mode_right)
  denom <- trapz_integral(grid[between_modes], pair_support[between_modes])
  in_interval <- grid >= interval$left & grid <= interval$right
  entropy_score <- if (denom > 0 && any(in_interval)) {
    trapz_integral(
      grid[in_interval],
      entropy_vals[in_interval] * pair_support[in_interval]
    ) /
      denom
  } else {
    0
  }
  n_left_interval <- sum(
    left_values >= interval$left & left_values <= interval$right
  )
  n_right_interval <- sum(
    right_values >= interval$left & right_values <= interval$right
  )
  balance_score <- if ((n_left_interval + n_right_interval) > 0) {
    4 *
      n_left_interval *
      n_right_interval /
      (n_left_interval + n_right_interval)^2
  } else {
    0
  }
  boot <- bootstrap_boundary_stability(
    left_values = left_values,
    right_values = right_values,
    global_min = global_min,
    global_max = global_max,
    n_grid = n_grid,
    bandwidth = bandwidth,
    n_boot = n_boot
  )
  width_score <- clip01(
    (interval$right - interval$left) /
      max(abs(mode_right - mode_left), .Machine$double.eps)
  )
  support_score <- clip01(
    0.30 *
      overlap_score +
      0.25 * clip01(entropy_score) +
      0.20 * balance_score +
      0.15 * boot$stability +
      0.10 * width_score
  )
  keep <- (interval$right > interval$left) &&
    (support_score >= 0.45) &&
    (n_left_interval > 0) &&
    (n_right_interval > 0)

  list(
    boundary_row = data.frame(
      boundary_id = paste0(left_group, "__", right_group),
      left_group = left_group,
      right_group = right_group,
      mode_left = mode_left,
      mode_right = mode_right,
      boundary = boundary,
      left_transition = interval$left,
      right_transition = interval$right,
      overlap_score = overlap_score,
      entropy_score = clip01(entropy_score),
      balance_score = balance_score,
      stability_score = boot$stability,
      width_score = width_score,
      support_score = support_score,
      keep = keep,
      inferable = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

build_windows_and_cells <- function(
  meta,
  groups,
  pseudotime_column,
  boundaries,
  global_min,
  global_max
) {
  pt <- meta[[pseudotime_column]]
  names(pt) <- rownames(meta)
  windows <- list()
  cells <- list()
  metrics <- list()
  current_left <- global_min

  for (i in seq_along(groups)) {
    stable_id <- state_id_stable(i)
    stable_right <- if (i <= nrow(boundaries) && isTRUE(boundaries$keep[[i]])) {
      boundaries$left_transition[[i]]
    } else if (i <= nrow(boundaries)) {
      boundaries$boundary[[i]]
    } else {
      global_max
    }
    stable_right <- max(stable_right, current_left)
    group_i <- groups[[i]]
    if (i < length(groups)) {
      stable_cells <- rownames(meta)[
        meta$.group == group_i & pt >= current_left & pt < stable_right
      ]
    } else {
      stable_cells <- rownames(meta)[
        meta$.group == group_i & pt >= current_left & pt <= stable_right
      ]
    }
    windows[[stable_id]] <- data.frame(
      state_id = stable_id,
      state_type = "stable",
      left = current_left,
      right = stable_right,
      parent_groups = group_i,
      stringsAsFactors = FALSE
    )
    cells[[stable_id]] <- stable_cells
    metrics[[stable_id]] <- data.frame(
      state_id = stable_id,
      state_type = "stable",
      n_cells = length(stable_cells),
      n_left = length(stable_cells),
      n_right = 0,
      support_score = NA_real_,
      keep = TRUE,
      inferable = FALSE,
      stringsAsFactors = FALSE
    )

    if (i <= nrow(boundaries)) {
      transition_id <- state_id_transition(i, i + 1L)
      trans_left <- if (isTRUE(boundaries$keep[[i]])) {
        boundaries$left_transition[[i]]
      } else {
        boundaries$boundary[[i]]
      }
      trans_right <- if (isTRUE(boundaries$keep[[i]])) {
        boundaries$right_transition[[i]]
      } else {
        boundaries$boundary[[i]]
      }
      trans_left <- max(trans_left, stable_right)
      trans_right <- max(trans_right, trans_left)
      left_group <- groups[[i]]
      right_group <- groups[[i + 1L]]
      trans_cells <- rownames(meta)[
        meta$.group %in%
          c(left_group, right_group) &
          pt >= trans_left &
          pt <= trans_right
      ]
      n_left <- sum(meta[trans_cells, ".group", drop = TRUE] == left_group)
      n_right <- sum(meta[trans_cells, ".group", drop = TRUE] == right_group)
      windows[[transition_id]] <- data.frame(
        state_id = transition_id,
        state_type = "transition",
        left = trans_left,
        right = trans_right,
        parent_groups = paste(left_group, right_group, sep = "|"),
        stringsAsFactors = FALSE
      )
      cells[[transition_id]] <- trans_cells
      metrics[[transition_id]] <- data.frame(
        state_id = transition_id,
        state_type = "transition",
        n_cells = length(trans_cells),
        n_left = n_left,
        n_right = n_right,
        support_score = NA_real_,
        keep = FALSE,
        inferable = FALSE,
        stringsAsFactors = FALSE
      )
      current_left <- trans_right
    }
  }

  list(
    windows = do.call(rbind, windows),
    cells = cells,
    metrics = do.call(rbind, metrics)
  )
}

adjust_transition_intervals <- function(boundaries, global_min, global_max) {
  left <- boundaries$left_transition
  right <- boundaries$right_transition
  keep <- boundaries$keep
  if (length(left) == 0) {
    return(list(left = left, right = right, keep = keep))
  }
  left <- pmax(left, global_min)
  right <- pmin(right, global_max)
  for (i in seq_len(length(left) - 1L)) {
    split <- (boundaries$boundary[[i]] + boundaries$boundary[[i + 1L]]) / 2
    if (isTRUE(keep[[i]]) && right[[i]] > split) {
      right[[i]] <- split
    }
    if (isTRUE(keep[[i + 1L]]) && left[[i + 1L]] < split) {
      left[[i + 1L]] <- split
    }
  }
  keep <- keep & (right > left)
  list(left = left, right = right, keep = keep)
}

state_windows_from_result <- function(
  object,
  density_result,
  pseudotime_column,
  prefer_network = FALSE
) {
  if (!isTRUE(prefer_network)) {
    return(density_result$windows)
  }
  nets <- tryCatch(
    {
      active_network <- .multicsn_resolve_network(
        object,
        network = NULL,
        preferred = "dynamic",
        verbose = FALSE,
        caller = "plot_states"
      )
      GetNetwork(object, network = active_network)
    },
    error = function(e) NULL
  )
  if (is.list(nets) && length(nets) > 0) {
    state_ids <- names(nets)
    windows <- purrr::map_dfr(state_ids, function(sid) {
      net <- nets[[sid]]
      if (is.null(net) || !methods::is(net, "Network")) {
        return(NULL)
      }
      window <- net@params$state_window %ss% net@params$pseudotime_range
      if (is.null(window) || length(window) != 2) {
        row <- density_result$windows[
          density_result$windows$state_id == sid, ,
          drop = FALSE
        ]
        if (nrow(row) == 0) {
          return(NULL)
        }
        window <- c(row$left[[1]], row$right[[1]])
      }
      row <- density_result$windows[
        density_result$windows$state_id == sid, ,
        drop = FALSE
      ]
      if (nrow(row) == 0) {
        row <- data.frame(
          state_id = sid,
          state_type = net@params$state_type %ss% NA_character_,
          left = as.numeric(window[1]),
          right = as.numeric(window[2]),
          parent_groups = paste(
            net@params$parent_groups %ss% NA_character_,
            collapse = "|"
          ),
          keep = TRUE,
          inferable = TRUE,
          stringsAsFactors = FALSE
        )
      } else {
        row$left <- as.numeric(window[1])
        row$right <- as.numeric(window[2])
        row$keep <- TRUE
      }
      row
    })
    if (nrow(windows) > 0) {
      return(windows)
    }
  }
  density_result$windows
}

density_on_grid <- function(
  values,
  from,
  to,
  n_grid = 512,
  bandwidth = "nrd0"
) {
  values <- values[is.finite(values)]
  if (length(values) < 2 || stats::sd(values) == 0) {
    center <- if (length(values) == 0) mean(c(from, to)) else values[[1]]
    x <- seq(from, to, length.out = n_grid)
    y <- stats::dnorm(x, mean = center, sd = max((to - from) / 100, 1e-3))
  } else {
    den <- stats::density(
      values,
      from = from,
      to = to,
      n = n_grid,
      bw = bandwidth
    )
    x <- den$x
    y <- den$y
  }
  y <- y / trapz_integral(x, y)
  list(x = x, y = y, mode = x[which.max(y)][1])
}

find_density_boundary <- function(
  grid,
  left_density,
  right_density,
  mode_left,
  mode_right
) {
  lo <- min(mode_left, mode_right)
  hi <- max(mode_left, mode_right)
  idx <- which(grid >= lo & grid <= hi)
  x <- grid[idx]
  d <- left_density[idx] - right_density[idx]
  if (length(x) < 2) {
    return(mean(c(mode_left, mode_right)))
  }
  roots <- c()
  total <- c()
  for (k in seq_len(length(x) - 1L)) {
    if (!is.finite(d[[k]]) || !is.finite(d[[k + 1L]])) {
      next
    }
    if (d[[k]] == 0) {
      roots <- c(roots, x[[k]])
      total <- c(total, left_density[idx][[k]] + right_density[idx][[k]])
    } else if (d[[k]] * d[[k + 1L]] < 0) {
      root <- x[[k]] - d[[k]] * (x[[k + 1L]] - x[[k]]) / (d[[k + 1L]] - d[[k]])
      roots <- c(roots, root)
      total <- c(
        total,
        approx(
          x = x[c(k, k + 1L)],
          y = left_density[idx][c(k, k + 1L)] +
            right_density[idx][c(k, k + 1L)],
          xout = root,
          rule = 2
        )$y
      )
    }
  }
  if (length(roots) == 0) {
    return(mean(c(mode_left, mode_right)))
  }
  roots[[which.min(total)]]
}

find_transition_interval <- function(
  grid,
  posterior_left,
  boundary,
  threshold = 0.8
) {
  dominance <- pmax(posterior_left, 1 - posterior_left)
  eligible <- dominance <= threshold
  if (!any(eligible)) {
    return(list(left = boundary, right = boundary))
  }
  boundary_idx <- which.min(abs(grid - boundary))
  if (!eligible[[boundary_idx]]) {
    candidate_idx <- which(eligible)
    if (length(candidate_idx) == 0) {
      return(list(left = boundary, right = boundary))
    }
    boundary_idx <- candidate_idx[[which.min(abs(
      candidate_idx - boundary_idx
    ))]]
  }
  left_idx <- boundary_idx
  right_idx <- boundary_idx
  while (left_idx > 1 && eligible[[left_idx - 1L]]) {
    left_idx <- left_idx - 1L
  }
  while (right_idx < length(grid) && eligible[[right_idx + 1L]]) {
    right_idx <- right_idx + 1L
  }
  list(left = grid[[left_idx]], right = grid[[right_idx]])
}

bootstrap_boundary_stability <- function(
  left_values,
  right_values,
  global_min,
  global_max,
  n_grid,
  bandwidth,
  n_boot
) {
  if (n_boot <= 1 || length(left_values) < 3 || length(right_values) < 3) {
    return(list(stability = 0))
  }
  boundaries <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    left_boot <- sample(left_values, length(left_values), replace = TRUE)
    right_boot <- sample(right_values, length(right_values), replace = TRUE)
    left_density <- density_on_grid(
      left_boot,
      global_min,
      global_max,
      n_grid,
      bandwidth
    )
    right_density <- density_on_grid(
      right_boot,
      global_min,
      global_max,
      n_grid,
      bandwidth
    )
    boundaries[[b]] <- find_density_boundary(
      grid = left_density$x,
      left_density = left_density$y,
      right_density = right_density$y,
      mode_left = left_density$mode,
      mode_right = right_density$mode
    )
  }
  boundaries <- boundaries[is.finite(boundaries)]
  if (length(boundaries) < 2) {
    return(list(stability = 0))
  }
  spread <- stats::IQR(boundaries)
  ref <- abs(stats::median(right_values) - stats::median(left_values))
  list(stability = clip01(1 - spread / max(ref, .Machine$double.eps)))
}

trapz_integral <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) {
    return(0)
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}


clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

state_order_key <- function(sid) {
  nums <- as.numeric(regmatches(sid, gregexpr("[0-9]+", sid))[[1]])
  if (length(nums) == 0) {
    return(0)
  }
  if (length(nums) == 1) {
    return(nums[[1]])
  }
  nums[[1]] + 0.5
}
