#' @include setClass.R
#' @include setGenerics.R
#' @title Calculate gene rank
#'
#' @param object Network object
#' @param ... Additional arguments
#'
#' @return Data frame with gene ranks
#'
#' @rdname calculate_gene_rank
#' @export
setGeneric(
  name = "calculate_gene_rank",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("calculate_gene_rank")
  }
)

#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @return The Network object with gene ranks stored in \code{params$gene_rank}.
#' @rdname calculate_gene_rank
setMethod(
  "calculate_gene_rank",
  signature(object = "Network"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE,
           method = c("page_rank", "degree_distribution"),
           ...) {
    method <- match.arg(method)
    if (method == "page_rank") {
      ranks <- calculate_page_rank(object, regulators, targets, directed, ...)
    } else {
      ranks <- calculate_degree_distribution(object, regulators, targets, directed, ...)
    }
    object@params$gene_rank <- list(df = ranks, method = method)
    return(object)
  }
)

#' @param network Name of the network to use.
#' @param celltypes Character vector specifying which states/cell types to include; NULL for all.
#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @return The Seurat object with gene ranks stored in each Network's \code{params$gene_rank}.
#' @rdname calculate_gene_rank
#' @export
setMethod(
  "calculate_gene_rank",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           regulators = NULL,
           targets = NULL,
           directed = FALSE,
           method = c("page_rank", "degree_distribution"),
           ...) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "active",
      verbose = TRUE,
      caller = "calculate_gene_rank"
    )
    method <- match.arg(method)
    nets <- GetNetwork(object, network = network, celltypes = celltypes)
    if (is.null(nets) || length(nets) == 0) {
      stop("No networks found for the specified network or cell types.")
    }
    cell_names <- names(nets)
    if (is.null(cell_names)) cell_names <- seq_along(nets)
    for (i in seq_along(nets)) {
      net <- nets[[i]]
      if (is.null(net)) next
      net <- calculate_gene_rank(net, regulators = regulators, targets = targets, directed = directed, method = method, ...)
      object <- .multicsn_set_network_entry(object, network, cell_names[i], net)
    }
    return(object)
  }
)

#' @rdname calculate_gene_rank
#' @export
setMethod(
  "calculate_gene_rank",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)

#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @rdname calculate_gene_rank
#' @export
setMethod(
  "calculate_gene_rank",
  signature(object = "data.frame"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE,
           method = c("page_rank", "degree_distribution"),
           ...) {
    method <- match.arg(method)
    if (method == "page_rank") {
      return(
        calculate_page_rank(
          object,
          regulators,
          targets,
          directed,
          ...
        )
      )
    } else {
      return(
        calculate_degree_distribution(
          object,
          regulators,
          targets,
          directed,
          ...
        )
      )
    }
  }
)

#' @title Calculate PageRank
#'
#' @param object Network object
#' @param ... Additional arguments
#' @return Data frame with PageRank scores
#' @export
#' @rdname calculate_page_rank
setGeneric(
  name = "calculate_page_rank",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("calculate_page_rank")
  }
)

#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @rdname calculate_page_rank
#' @export
setMethod(
  "calculate_page_rank",
  signature(object = "Network"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- as.data.frame(object@network)
    .calculate_page_rank(network_table, directed)
  }
)

#' @param regulators Character vector, regulators to include
#' @param targets Character vector, targets to include
#' @param directed Logical, whether the network is directed
#' @rdname calculate_page_rank
#' @export
setMethod(
  "calculate_page_rank",
  signature(object = "data.frame"),
  function(object,
           regulators = NULL,
           targets = NULL,
           directed = FALSE) {
    network_table <- inferCSN::network_format(
      object,
      regulators,
      targets,
      abs_weight = FALSE
    )
    .calculate_page_rank(network_table, directed)
  }
)

.calculate_page_rank <- function(network_table, directed) {
  network <- igraph::graph_from_data_frame(
    network_table,
    directed = directed
  )
  page_rank_res <- data.frame(
    igraph::page_rank(network, directed = directed)$vector
  )
  colnames(page_rank_res) <- c("page_rank")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "page_rank")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$page_rank,
    decreasing = TRUE
  ), ]
  page_rank_res$regulator <- ifelse(
    page_rank_res$gene %in% unique(network_table$regulator),
    "TRUE", "FALSE"
  )
  rownames(page_rank_res) <- NULL

  return(page_rank_res)
}

#' @title calculate_degree_distribution
#'
#' @param object The input object.
#' @param ... Parameters for other methods.
#'
#' @rdname calculate_degree_distribution
#'
#' @export
#' @export
setGeneric(
  name = "calculate_degree_distribution",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("calculate_degree_distribution")
  }
)

#' @param regulators regulators
#' @param targets targets
#' @param directed directed
#' @rdname calculate_degree_distribution
#' @export
setMethod(
  "calculate_degree_distribution",
  signature(object = "Network"),
  function(object, regulators = NULL, targets = NULL, directed = TRUE, ...) {
    network_table <- as.data.frame(object@network)

    .degree_distribution(network_table, directed)
  }
)

#' @param regulators regulators
#' @param targets targets
#' @param directed directed
#' @rdname calculate_degree_distribution
#' @export
setMethod(
  "calculate_degree_distribution",
  signature(object = "data.frame"),
  function(object, regulators = NULL, targets = NULL, directed = TRUE, ...) {
    network_table <- inferCSN::network_format(
      object,
      regulators,
      targets,
      abs_weight = FALSE
    )
    .degree_distribution(network_table, directed)
  }
)

.calculate_power_fit <- function(deg) {
  degree_freq <- table(deg)
  k <- as.numeric(names(degree_freq))
  pk <- as.numeric(degree_freq) / sum(degree_freq)

  if (length(unique(k)) > 1) {
    fit <- stats::lm(log(pk) ~ log(k))
    return(summary(fit)$r.squared)
  }
  return(0)
}

.degree_distribution <- function(
  network_table, directed
) {
  network <- igraph::graph_from_data_frame(
    network_table,
    directed = directed
  )

  total_degrees <- igraph::degree(network, mode = "total")

  result <- data.frame(
    gene = names(total_degrees),
    degree = as.numeric(total_degrees)
  )

  degree_freq <- table(total_degrees)
  result$P_k <- degree_freq[match(result$degree, names(degree_freq))] / sum(degree_freq)

  if (directed) {
    in_degrees <- igraph::degree(network, mode = "in")
    out_degrees <- igraph::degree(network, mode = "out")

    in_power_law_score <- .calculate_power_fit(in_degrees)
    out_power_law_score <- .calculate_power_fit(out_degrees)

    result$in_degree <- in_degrees[result$gene]
    result$out_degree <- out_degrees[result$gene]
    result$in_power_law_fit <- in_power_law_score
    result$out_power_law_fit <- out_power_law_score
    result$rank_value <- (result$in_degree * in_power_law_score +
      result$out_degree * out_power_law_score) / 2
  } else {
    power_law_score <- .calculate_power_fit(total_degrees)
    result$power_law_fit <- power_law_score
    result$rank_value <- result$degree * power_law_score
  }

  result <- result[order(result$rank_value, decreasing = TRUE), ]
  result$regulator <- result$gene %in% unique(network_table$regulator)
  rownames(result) <- NULL

  return(result)
}

#' @title Plot gene ranks and network properties
#'
#' @param object Network object
#' @param ... Other params
#' @return Combined ggplot object
#' @export
#' @rdname plot_gene_rank
setGeneric(
  name = "plot_gene_rank",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_gene_rank")
  }
)

#' @param method Character, ranking method: "page_rank" or "degree_distribution"
#' @param weight_cutoff Numeric, threshold for edge weight filtering
#' @param compare_random Logical, whether to compare with randomized network
#' @param top_n Integer, number of top genes in centrality plot
#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "Network"),
  function(object,
           method = c("page_rank", "degree_distribution"),
           weight_cutoff = 0.1,
           compare_random = TRUE,
           top_n = 30,
           ...) {
    method <- match.arg(method)
    network_table <- tryCatch(export_csn(object), error = function(e) NULL)
    if (is.null(network_table) || nrow(network_table) == 0) {
      stop("No network edges found. Run inferCSN and ensure network has edges.")
    }
    network_table <- as.data.frame(network_table)
    if (!"weight" %in% colnames(network_table) || !is.numeric(network_table$weight)) {
      stop("Network table must have numeric 'weight' column.")
    }
    network_table <- network_table[abs(network_table$weight) >= weight_cutoff, , drop = FALSE]

    g_orig <- igraph::graph_from_data_frame(network_table, directed = TRUE)
    degrees_orig <- igraph::degree(g_orig, mode = "total")

    degree_freq <- table(degrees_orig)
    df_orig <- data.frame(
      k = as.numeric(names(degree_freq)),
      P_k = as.numeric(degree_freq) / sum(degree_freq)
    )

    dist_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 10),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        aspect.ratio = 1
      )

    p1 <- ggplot(df_orig, aes(x = log(k), y = log(P_k))) +
      geom_point(size = 1) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        color = "grey50",
        linewidth = 0.3,
        formula = y ~ x
      ) +
      dist_theme +
      labs(
        x = "log k",
        y = "log P(k)",
        title = "Degree distribution",
        tag = "a"
      )

    if (nrow(df_orig) > 1) {
      model <- stats::lm(log(P_k) ~ log(k), data = df_orig)
      r2 <- summary(model)$r.squared
      p1 <- p1 + annotate("text",
        x = max(log(df_orig$k)) - 0.1,
        y = max(log(df_orig$P_k)) - 0.1,
        label = sprintf("R\u00b2 = %.2f", r2),
        size = 3,
        color = "steelblue",
        hjust = 1
      )
    }

    stored <- object@params$gene_rank
    if (is.null(stored) || stored$method != method) {
      object <- calculate_gene_rank(object, method = method)
      ranks <- object@params$gene_rank$df
    } else {
      ranks <- stored$df
    }
    rank_col <- if (method == "page_rank") "page_rank" else "rank_value"

    centrality_df <- ranks %>%
      dplyr::arrange(dplyr::desc(.data[[rank_col]])) %>%
      head(top_n) %>%
      dplyr::mutate(centrality = .data[[rank_col]] / max(.data[[rank_col]]))

    cent_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
      )

    p_cent <- ggplot(
      centrality_df,
      aes(x = centrality, y = stats::reorder(gene, centrality))
    ) +
      geom_point(color = "steelblue", size = 2) +
      cent_theme +
      labs(
        x = if (method == "page_rank") "PageRank centrality" else "Degree centrality",
        y = NULL,
        tag = "c"
      )

    if (compare_random) {
      n_edges <- igraph::ecount(g_orig)
      all_nodes <- unique(c(network_table$regulator, network_table$target))

      random_edges <- data.frame(
        regulator = sample(all_nodes, n_edges, replace = TRUE),
        target = sample(all_nodes, n_edges, replace = TRUE),
        weight = network_table$weight
      )

      g_random <- igraph::graph_from_data_frame(random_edges, directed = TRUE)
      degrees_random <- igraph::degree(g_random, mode = "total")

      degree_freq_random <- table(degrees_random)
      df_random <- data.frame(
        k = as.numeric(names(degree_freq_random)),
        P_k = as.numeric(degree_freq_random) / sum(degree_freq_random)
      )

      p2 <- ggplot(df_random, aes(x = log(k), y = log(P_k))) +
        geom_point(size = 1) +
        geom_smooth(
          method = "lm",
          se = FALSE,
          color = "grey50",
          linewidth = 0.3,
          formula = y ~ x
        ) +
        dist_theme +
        labs(
          x = "log k", y = "log P(k)",
          title = "Degree distribution\nof randomized network",
          tag = "b"
        )

      if (nrow(df_random) > 1) {
        model_random <- stats::lm(log(P_k) ~ log(k), data = df_random)
        r2_random <- summary(model_random)$r.squared
        p2 <- p2 + annotate("text",
          x = max(log(df_random$k)) - 0.1,
          y = max(log(df_random$P_k)) - 0.1,
          label = sprintf("R\u00b2 = %.2f", r2_random),
          size = 3,
          color = "steelblue",
          hjust = 1
        )
      }
    }

    if (compare_random) {
      dist_plots <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1))
      return(dist_plots | p_cent + patchwork::plot_layout(widths = c(1, 1.5)))
    }
    return(p1 | p_cent + patchwork::plot_layout(widths = c(1, 1.5)))

    return(
      list(
        distribution = if (compare_random) p1 / p2 else p1,
        centrality = p_cent
      )
    )
  }
)

#' @param top_n Integer, number of top genes to show in centrality plot.
#' @param network Name of the network to use; defaults to the active network.
#' @param celltypes Character vector specifying which states/cell types to include; NULL for all.
#' @param combine Logical. When \code{TRUE}, combine per-state plots with
#'   \code{patchwork::wrap_plots}; otherwise return a list of plots.
#' @param nrow,ncol Integer or \code{NULL}. Rows/columns for \code{patchwork::wrap_plots} when \code{combine=TRUE}.
#' @param byrow Logical. Fill plots by row in \code{wrap_plots}; default \code{TRUE}.
#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           method = c("page_rank", "degree_distribution"),
           weight_cutoff = 0.1,
           compare_random = TRUE,
           combine = TRUE,
           top_n = 30,
           nrow = NULL,
           ncol = NULL,
           byrow = TRUE,
           ...) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "active",
      verbose = TRUE,
      caller = "plot_gene_rank"
    )
    method <- match.arg(method)
    ranks_list <- GeneRanks(object, network = network, celltypes = celltypes)
    need_compute <- is.null(ranks_list) || any(sapply(ranks_list, function(r) is.null(r) || is.null(r$df)))
    if (need_compute) {
      object <- calculate_gene_rank(object, network = network, celltypes = celltypes, method = method, ...)
      ranks_list <- GeneRanks(object, network = network, celltypes = celltypes)
    }
    cell_names <- names(ranks_list)
    if (is.null(cell_names)) cell_names <- seq_along(ranks_list)
    plot_list <- list()
    for (i in seq_along(ranks_list)) {
      ranks <- ranks_list[[i]]
      if (is.null(ranks) || is.null(ranks$df) || nrow(ranks$df) == 0) next
      p <- .plot_rank_centrality(ranks$df, ranks$method, top_n, title = cell_names[i])
      plot_list[[cell_names[i]]] <- p
    }
    if (length(plot_list) == 0) stop("No gene ranks to plot.")
    if (combine && length(plot_list) > 1) {
      wrap_args <- list(plot_list, nrow = nrow, ncol = ncol, byrow = byrow)
      return(do.call(patchwork::wrap_plots, wrap_args))
    }
    if (combine && length(plot_list) == 1) {
      return(plot_list[[1]])
    }
    plot_list
  }
)

#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)

.plot_rank_centrality <- function(ranks_df, method, top_n = 30, title = NULL) {
  rank_col <- if (method == "page_rank") "page_rank" else "rank_value"
  if (!rank_col %in% colnames(ranks_df)) rank_col <- colnames(ranks_df)[2]
  centrality_df <- ranks_df %>%
    dplyr::arrange(dplyr::desc(.data[[rank_col]])) %>%
    head(top_n) %>%
    dplyr::mutate(centrality = .data[[rank_col]] / max(.data[[rank_col]], na.rm = TRUE))
  gene_col <- if ("gene" %in% colnames(centrality_df)) "gene" else colnames(centrality_df)[1]
  p <- ggplot(centrality_df, aes(x = .data[["centrality"]], y = stats::reorder(.data[[gene_col]], .data[["centrality"]]))) +
    geom_point(color = "steelblue", size = 2) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10)
    ) +
    labs(x = if (method == "page_rank") "PageRank" else "Degree centrality", y = NULL, title = title)
  p
}

#' @rdname plot_gene_rank
#' @export
setMethod(
  "plot_gene_rank",
  signature(object = "data.frame"),
  function(object,
           method = c("page_rank", "degree_distribution"),
           weight_cutoff = 0.1,
           compare_random = TRUE,
           top_n = 30) {
    if ("gene" %in% colnames(object) && ("page_rank" %in% colnames(object) || "rank_value" %in% colnames(object))) {
      method_used <- if ("page_rank" %in% colnames(object)) "page_rank" else "degree_distribution"
      .plot_rank_centrality(object, method_used, top_n)
    } else if ("regulator" %in% colnames(object) && "target" %in% colnames(object) && "weight" %in% colnames(object)) {
      network_obj <- new("Network", network = object)
      plot_gene_rank(network_obj, method, weight_cutoff, compare_random)
    } else {
      stop("data.frame must have either (gene, page_rank/rank_value) for ranks or (regulator, target, weight) for network.")
    }
  }
)
