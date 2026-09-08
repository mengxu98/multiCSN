#' @title quick plot of dynamic networks
#'
#' @param network the result of running epochGRN
#' @param regulators regulators
#' @param only_TFs plot only regulator network
#' @param network_order which epochs or transitions to plot
#' @param weight_threshold weight_threshold
#'
#' @return plot
#'
#' @export
plot_dynamic_network <- function(
  network,
  regulators,
  only_TFs = TRUE,
  network_order = NULL,
  weight_threshold = NULL
) {
  plot_list <- list()

  if (!is.null(network_order)) {
    network <- network[network_order]
  }

  for (i in 1:length(network)) {
    network_table <- network[[i]]

    colnames(network_table) <- c("regulator", "target", "weight")
    if (only_TFs) {
      network_table <- network_table[network_table$regulator %in% regulators, ]
    }
    if (!is.null(weight_threshold)) {
      network_table <- network_table[abs(network_table$weight) > weight_threshold, ]
    }
    network_table$interaction <- "Activation"
    network_table$interaction[network_table$weight < 0] <- "Repression"

    net <- igraph::graph_from_data_frame(
      network_table[, c("regulator", "target", "interaction")],
      directed = FALSE
    )
    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    network_plot_data <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    network_plot_data$is_regulator <- as.character(network_plot_data$name %in% regulators)

    cols <- c("Activation" = "blue", "Repression" = "red")
    plot_list[[i]] <- ggplot() +
      ggnetwork::geom_edges(
        data = network_plot_data,
        aes(x = x, y = y, xend = xend, yend = yend, color = interaction),
        size = 0.75,
        curvature = 0.1,
        alpha = 0.6
      ) +
      ggnetwork::geom_nodes(
        data = network_plot_data,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "darkgray",
        size = 6,
        alpha = 0.5
      ) +
      ggnetwork::geom_nodes(
        data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "#8C4985",
        size = 6,
        alpha = 0.8
      ) +
      scale_color_manual(values = cols) +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "FALSE", ],
        aes(x = x, y = y, label = name),
        size = 2,
        color = "#5A8BAD"
      ) +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
        aes(x = x, y = y, label = name),
        size = 3.5,
        color = "black"
      ) +
      theme_blank() +
      ggtitle(names(network)[i])

    common_legend <- get_legend(plot_list[[i]])
    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }

  plot_list$legend <- common_legend
  do.call(gridExtra::grid.arrange, plot_list)
}


#' Plot the dynamic differential network but colored by communities and optionally faded by igraph::betweenness
#'
#' @param network the dynamic network
#' @param regulators regulators
#' @param top_edges top_edges
#' @param only_TFs whether or not to only plot regulators and exclude non-regulators
#' @param network_order the network_order in which to plot epochs, or which epochs to plot
#' @param communities community assignments or the result of running find_commumities.
#' The names in this object should match the names of the epoch networks in network.
#' If NULL, it will be automatically run.
#' @param compute_betweenness whether or not to fade nodes by igraph::betweenness
#'
#' @return plot
#'
#' @export
plot_detail_network <- function(
  network,
  regulators,
  top_edges = NULL,
  only_TFs = TRUE,
  network_order = NULL,
  communities = NULL,
  compute_betweenness = TRUE
) {
  plot_list <- list()

  if (!is.null(communities)) {
    if (names(communities) != names(network)) {
      names(communities) <- names(network)
    }
  }

  if (!is.null(network_order)) {
    network <- network[network_order]
    if (!is.null(communities)) {
      communities <- communities[network_order]
    }
  }

  if (any(sapply(network, function(x) {
    ("interaction" %in% names(x))
  }) == FALSE)) {
    network <- add_interaction_type(network)
  }

  for (i in seq_len(length(network))) {
    message(
      "\rPlotting for ", i, "/", length(network), " networks.",
      appendLF = FALSE
    )

    network_table <- network[[i]]

    if (!is.null(top_edges)) {
      network_table <- network_table[1:top_edges, ]
    }

    if (only_TFs) {
      network_table <- network_table[network_table$target %in% regulators, ]
    }
    if (nrow(network_table) == 0) {
      next
    }

    network_table <- network_table[, c("regulator", "target", "weight", "interaction")]
    network_table$weight <- abs(network_table$weight)

    net <- igraph::graph_from_data_frame(
      network_table,
      directed = FALSE
    )

    if (compute_betweenness) {
      betweenness <- igraph::betweenness(
        net,
        directed = FALSE,
        normalized = TRUE
      )
      betweenness <- as.data.frame(betweenness)
      colnames(betweenness) <- "betweenness"
      betweenness$gene <- rownames(betweenness)
      if (!is.null(communities)) {
        communities <- communities
      } else {
        communities <- as.data.frame(
          as.table(
            igraph::membership(igraph::cluster_louvain(net))
          )
        )
      }
      colnames(communities) <- c("gene", "communities")
      vtx_features <- merge(betweenness, communities, by = "gene", all = TRUE)
    } else {
      if (!is.null(communities)) {
        communities <- communities
      } else {
        communities <- as.data.frame(
          as.table(
            igraph::membership(igraph::cluster_louvain(net))
          )
        )
      }
      colnames(communities) <- c("gene", "communities")
      vtx_features <- communities
    }

    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]

    network_plot_data <- ggnetwork::ggnetwork(
      net,
      layout = layout_ordered,
      cell.jitter = 0
    )

    network_plot_data$is_regulator <- as.character(
      network_plot_data$name %in% regulators
    )
    if (compute_betweenness) {
      network_plot_data$betweenness <- vtx_features$betweenness[match(network_plot_data$name, vtx_features$gene)]
    }
    network_plot_data$communities <- as.factor(
      vtx_features$communities[match(network_plot_data$name, vtx_features$gene)]
    )

    cols <- c("Activation" = "blue", "Repression" = "red")
    num_cols2 <- length(unique(network_plot_data$communities))
    if (num_cols2 <= 8) {
      cols2 <- RColorBrewer::brewer.pal(num_cols2, "Set1")
    } else {
      cols2 <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Set1")
      )(num_cols2)
    }
    names(cols2) <- unique(network_plot_data$communities)
    cols <- c(cols, cols2)

    p <- ggplot() +
      ggnetwork::geom_edges(
        data = network_plot_data,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend,
          linewidth = weight,
          color = interaction
        ),
        curvature = 0.1,
        alpha = 0.6
      ) +
      scale_color_manual(values = cols) +
      ggnetwork::geom_nodes(
        data = network_plot_data[network_plot_data$is_regulator == "FALSE", ],
        aes(x = x, y = y),
        color = "darkgray",
        alpha = 0.5,
        size = 6
      )
    if (compute_betweenness) {
      p <- p +
        ggnetwork::geom_nodes(
          data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
          aes(x = x, y = y, color = communities, alpha = betweenness),
          size = 6
        )
    } else {
      p <- p +
        ggnetwork::geom_nodes(
          data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
          aes(x = x, y = y, color = communities),
          alpha = .5,
          size = 6
        )
    }
    p <- p +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
        aes(x = x, y = y, label = name),
        size = 2.5,
        color = "#8C4985"
      ) +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "FALSE", ],
        aes(x = x, y = y, label = name),
        size = 2.5,
        color = "#5A8BAD"
      ) +
      theme_blank() +
      ggtitle(names(network)[i]) +
      theme(legend.position = "none")

    plot_list[[i]] <- p
  }
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  do.call(gridExtra::grid.arrange, plot_list)
}

#' quick plot of top regulators in dynamic networks
#'
#' @param network_list the result of running epochGRN
#' @param gene_ranks the result of running compute_pagerank
#' @param regulators regulators
#' @param targets targets
#' @param regulators_num number of top regulators to plot
#' @param targets_num number of top targets to plot for each regulator
#' @param network_order which epochs or transitions to plot
#' @param method Choose method to rank regulators for plot.
#'  Defaule set to `weight`, mean choose regulators rely weights in network infer by \link[inferCSN]{inferCSN}.
#'  Also can choose `page_rank`, then will ranking regulators rely the result of \link[igraph]{page_rank}.
#'
#' @return plot
#'
#' @export
plot_top_features <- function(
  network_list,
  gene_ranks = NULL,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  targets_num = NULL,
  network_order = NULL,
  method = "weight"
) {
  plot_list <- list()

  if (!is.null(network_order)) {
    network_list <- network_list[network_order]
  }

  common_legend <- NULL
  for (i in seq_len(length(network_list))) {
    epoch <- names(network_list)[i]
    network_table <- network_list[[i]]
    thisutils::log_message(
      paste0("Plotting for ", i, "/", length(network_list), " network: ", epoch)
    )

    if (is.null(regulators)) {
      if (method == "weight") {
        top_regulators <- network_table$regulator
      } else if (method == "page_rank") {
        if (!is.null(gene_ranks)) {
          rank <- gene_ranks[[epoch]]
        } else {
          rank <- calculate_gene_rank(network_table)
        }
        top_regulators <- rownames(rank[rank$is_regulator == TRUE, ])
      }
    } else {
      top_regulators <- regulators
    }

    top_regulators <- top_regulators[!duplicated(top_regulators)]

    if (!is.null(regulators_num)) {
      top_regulators <- top_regulators[1:min(regulators_num, length(top_regulators))]
    }

    network_table_sub <- network_table[network_table$regulator %in% top_regulators, ]
    top_regulators <- unique(network_table_sub$regulator)
    if (length(top_regulators) == 0) {
      message("No network found.")
      next
    }

    network_table_plot <- purrr::map_dfr(
      top_regulators, function(x) {
        if (is.null(targets)) {
          targets <- as.character(
            network_table[network_table$regulator == x, "target"]
          )
          if (method == "weight") {
            targets <- network_table[network_table$target %in% targets, "target"]
          } else if (method == "page_rank") {
            rank_targets <- rank[targets, ]
            rank_targets <- rank_targets[order(
              rank_targets$page_rank,
              decreasing = TRUE
            ), ]
            targets <- rownames(rank_targets)
          }
        }
        targets <- targets[!duplicated(targets)]
        if (!is.null(targets_num)) {
          targets <- targets[1:min(targets_num, length(targets))]
        }

        add <- network_table[network_table$regulator == x, ]
        add[add$target %in% targets, ]
      }
    )

    network_table_plot <- add_interaction_type(network_table_plot)
    network_table_plot$weight <- abs(network_table_plot$weight)
    net <- igraph::graph_from_data_frame(
      network_table_plot,
      directed = FALSE
    )

    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    network_plot_data <- ggnetwork::ggnetwork(net, layout = layout_ordered, cell.jitter = 0)

    network_plot_data$is_regulator <- as.character(network_plot_data$name %in% top_regulators)

    cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

    plot_list[[i]] <- ggplot() +
      ggnetwork::geom_edges(
        data = network_plot_data,
        aes(
          x = x, y = y,
          xend = xend, yend = yend,
          color = interaction
        ),
        linewidth = 0.75,
        curvature = 0.1,
        alpha = 0.6
      ) +
      ggnetwork::geom_nodes(
        data = network_plot_data[network_plot_data$is_regulator == "FALSE", ],
        aes(x = x, y = y),
        color = "darkgray",
        size = 6,
        alpha = 0.8
      ) +
      ggnetwork::geom_nodes(
        data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
        aes(x = x, y = y),
        color = "#8C4985",
        size = 6,
        alpha = 0.8
      ) +
      scale_color_manual(values = cols) +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "FALSE", ],
        aes(x = x, y = y, label = name),
        size = 2,
        color = "#5A8BAD"
      ) +
      geom_nodelabel_repel(
        data = network_plot_data[network_plot_data$is_regulator == "TRUE", ],
        aes(x = x, y = y, label = name),
        size = 2.5,
        color = "#8C4985"
      ) +
      theme_blank() +
      ggtitle(names(network_list)[i])

    if (is.null(common_legend) | length(unique(network_table_plot$interaction)) == 2) {
      common_legend <- get_legend(plot_list[[i]])
    }

    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }
  plot_list$legend <- common_legend

  patchwork::wrap_plots(plot_list)
}

#' quick plot of top regulators given targets in dynamic networks based on reconstruction weight, colored by expression and interaction type
#'
#' @param network the result of running epochGRN
#' @param targets targets
#' @param epochs_list result of running assign_epochs
#' @param weight_column column name containing reconstruction weights to use
#' @param gene_ranks gene_ranks
#' @param regulators_num number of top regulators to plot
#' @param network_order which epochs or transitions to plot
#' @param fixed_layout whether or not to fix node positions across epoch networks
#' @param layout_alg layout algorithm if fixed_layout. Defaults to FR. Ignored if fixed_layout==FALSE.
#' @param declutter if TRUE, will only label nodes with active interactions in given network
#'
#' @return plot
#'
#' @export
plot_targets_with_top_regulators_detail <- function(
  network,
  targets,
  epochs_list,
  weight_column = "weight",
  gene_ranks = NULL,
  regulators_num = 5,
  network_order = NULL,
  fixed_layout = TRUE,
  layout_alg = "fr",
  declutter = TRUE
) {
  plot_list <- list()

  if (!is.null(network_order)) {
    network <- network[network_order]
  }


  ktgraph <- list()
  for (i in 1:length(network)) {
    mean_expression <- epochs_list[[i]]$mean_expression
    epoch <- names(network)[i]
    network_table <- network[[epoch]]


    tgs <- as.character(network_table[network_table$target %in% targets, "target"])
    if (length(tgs) == 0) {
      next
    }

    network_table$interaction <- "Activation"
    network_table$interaction[network_table$weight < 0] <- "Repression"

    edges_to_keep <- data.frame(regulator = character(), target = character())
    if (weight_column == "page_rank") {
      for (target in tgs) {
        if (is.null(gene_ranks)) {
          stop("Need to supply gene_ranks.")
        }
        rank <- gene_ranks[[epoch]]
        regs_of_targets <- as.character(network_table[network_table$target == target, "regulator"])
        rank_regs <- rank[regs_of_targets, ]
        rank_regs <- rank_regs[order(rank_regs$page_rank, decreasing = TRUE), ]
        top_regs <- rownames(rank_regs)[1:regulators_num]

        edges <- network_table[network_table$target == target, ]
        edges <- edges[edges$regulator %in% top_regs, c("regulator", "target", "interaction")]

        edges_to_keep <- rbind(edges_to_keep, data.frame(regulator = as.character(edges$regulator), target = as.character(edges$target), interaction = as.character(edges$interaction)))


        if (nrow(edges_to_keep) == 1) {
          edges_to_keep <- rbind(edges_to_keep, data.frame(regulator = NA, target = NA))
        }
      }
    } else {
      for (target in tgs) {
        edges <- network_table[network_table$target == target, ]
        edges <- edges[order(edges[, weight_column], decreasing = TRUE), ]
        edges <- edges[1:regulators_num, c("regulator", "target", "interaction")]

        edges_to_keep <- rbind(edges_to_keep, data.frame(regulator = as.character(edges$regulator), target = as.character(edges$target), interaction = as.character(edges$interaction)))
      }
    }

    ktgraph[[epoch]] <- edges_to_keep
  }


  if (fixed_layout) {
    agg <- dplyr::bind_rows(ktgraph)[, c("regulator", "target", "interaction")]
    agg <- agg[!duplicated(agg), ]
    agg <- igraph::graph_from_data_frame(agg, directed = FALSE)
    agg <- igraph::delete_vertices(agg, v = igraph::V(agg)$name[is.na(igraph::V(agg)$name)])


    if (layout_alg == "sugiyama") {
      layout <- igraph::layout_with_sugiyama(agg)
      layout <- layout$layout
    } else if (layout_alg == "mds") {
      layout <- igraph::layout_with_mds(agg)
    } else if (layout_alg == "lefttoright") {
      layout <- igraph::layout_with_mds(agg)
      agg_vtcs <- data.frame(name = igraph::V(agg)$name)

      temp_splits <- list()
      for (i in 1:length(ktgraph)) {
        temp_splits[[i]] <- union(ktgraph[[i]]$regulator, ktgraph[[i]]$target)
      }

      rownames(agg_vtcs) <- agg_vtcs$name
      agg_vtcs$avg_epoch <- NA
      for (row in rownames(agg_vtcs)) {
        agg_vtcs[row, "avg_epoch"] <- mean(which(sapply(temp_splits, FUN = function(x) row %in% x)))
      }

      layout[, 1] <- jitter(agg_vtcs$avg_epoch, factor = 5)
    } else {
      layout <- igraph::layout_with_fr(agg)
    }

    rownames(layout) <- igraph::V(agg)$name
  }


  for (i in 1:length(network)) {
    epoch <- names(network)[i]
    edges_to_keep <- ktgraph[[epoch]]


    if (fixed_layout) {
      net <- igraph::graph_from_data_frame(edges_to_keep[, c("regulator", "target", "interaction")], directed = FALSE)

      addvtcs <- igraph::V(agg)$name[!(igraph::V(agg)$name %in% igraph::V(net)$name)]
      net <- igraph::add_vertices(net, length(addvtcs), attr = list(name = addvtcs))
    } else {
      net <- igraph::graph_from_data_frame(edges_to_keep[, c("regulator", "target", "interaction")], directed = FALSE)
    }


    net <- igraph::delete_vertices(net, v = igraph::V(net)$name[is.na(igraph::V(net)$name)])
    net <- igraph::delete_vertices(net, v = igraph::V(net)$name[igraph::V(net)$name == "NA"])


    expression_from <- mean_expression[mean_expression$epoch == strsplit(epoch, split = "..", fixed = TRUE)[[1]][1], ]
    expression_to <- mean_expression[mean_expression$epoch == strsplit(epoch, split = "..", fixed = TRUE)[[1]][2], ]
    if (strsplit(epoch, split = "..", fixed = TRUE)[[1]][1] == strsplit(epoch, split = "..", fixed = TRUE)[[1]][2]) {
      igraph::V(net)$expression <- expression_from$mean_expression[match(igraph::V(net)$name, expression_from$gene)]
    } else {
      igraph::V(net)$expression <- ifelse(igraph::V(net)$name %in% edges_to_keep$regulator, expression_from$mean_expression[match(igraph::V(net)$name, expression_from$gene)], expression_to$mean_expression[match(igraph::V(net)$name, expression_to$gene)])
    }


    if (fixed_layout) {
      layout_ordered <- layout[igraph::V(net)$name, ]
      network_plot_data <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    } else {
      layout <- igraph::layout_with_fr(net)
      rownames(layout) <- igraph::V(net)$name
      layout_ordered <- layout[igraph::V(net)$name, ]
      network_plot_data <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    }


    network_plot_data$type <- "regulator"
    network_plot_data$type[network_plot_data$name %in% targets] <- "known target"
    network_plot_data$type <- factor(network_plot_data$type, levels = c("regulator", "known target"))

    network_plot_data <- network_plot_data[!(is.na(network_plot_data$name)), ]

    cols <- c("Activation" = "blue", "Repression" = "red")

    plot_list[[i]] <- ggplot() +
      geom_edges(data = network_plot_data, aes(x = x, y = y, xend = xend, yend = yend, color = interaction), size = 0.75, curvature = 0.1, alpha = .6) +
      ggnetwork::geom_nodes(data = network_plot_data, aes(x = x, y = y, xend = xend, yend = yend, shape = type, size = expression, alpha = expression), color = "black") +
      scale_color_manual(values = cols) +
      theme_blank() +
      ggtitle(names(network)[i])

    if (declutter) {
      keep <- union(edges_to_keep$regulator, edges_to_keep$target)
      plot_list[[i]] <- plot_list[[i]] + geom_nodelabel_repel(data = network_plot_data[network_plot_data$name %in% keep, ], aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD")
    } else {
      plot_list[[i]] <- plot_list[[i]] + geom_nodelabel_repel(data = network_plot_data, aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD")
    }

    common_legend <- get_legend(plot_list[[i]])

    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }


  plot_list[sapply(plot_list, is.null)] <- NULL

  plot_list$legend <- common_legend
  do.call(gridExtra::grid.arrange, plot_list)
}


#' Updated plot of top regulators given targets in dynamic networks based on a weight column.
#' Top regulators computed for each epoch, but maintained in plot across epochs if present in epoch subnetwork.
#'
#'
#' @param network the result of running epochGRN
#' @param targets targets
#' @param epochs result of running assign_epochs
#' @param weight_column column name containing reconstruction weights to use
#' @param gene_ranks gene_ranks
#' @param regulators_num number of top regulators to plot
#' @param network_order which epochs or transitions to plot
#' @param fixed_layout whether or not to fix node positions across epoch networks
#' @param declutter if TRUE, will only label nodes with active interactions in given network
#' @param show_expression if TRUE, size and shade of node indicates mean expression in a given epoch.
#' @param node_size node_size
#' @param label_size label_size
#' @param title_size title_size
#' @param legend_size legend_size
#'
#' @return plot
#'
#' @export
plot_targets_and_regulators <- function(
  network,
  targets,
  epochs = NULL,
  weight_column = "weight",
  gene_ranks = NULL,
  regulators_num = 5,
  network_order = NULL,
  fixed_layout = TRUE,
  declutter = TRUE,
  show_expression = TRUE,
  node_size = 5,
  label_size = 2,
  title_size = 10,
  legend_size = 8
) {
  plot_list <- list()


  if (!is.null(network_order)) {
    network <- network[network_order]
  }


  ktgraph <- list()
  for (i in 1:length(network)) {
    if (!is.null(epochs) & ("mean_expression" %in% names(epochs[[i]]))) {
      message("expression")
      mean_expression <- epochs[[i]]$mean_expression
    } else {
      show_expression <- FALSE
    }
    epoch <- names(network)[i]
    network_table <- network[[epoch]]


    tgs <- as.character(network_table[network_table$target %in% targets, "target"])
    if (length(tgs) == 0) {
      next
    }

    network_table$interaction <- "Activation"
    network_table$interaction[network_table$weight < 0] <- "Repression"

    edges_to_keep <- data.frame(regulator = character(), target = character())
    if (!is.null(gene_ranks) & (weight_column %in% colnames(gene_ranks[[1]]))) {
      for (target in tgs) {
        rank <- gene_ranks[[epoch]]
        regs_of_targets <- as.character(network_table[network_table$target == target, "regulator"])
        rank_regs <- rank[regs_of_targets, ]
        rank_regs <- rank_regs[order(rank_regs[, weight_column], decreasing = TRUE), ]
        top_regs <- rownames(rank_regs)[1:regulators_num]

        edges <- network_table[network_table$target == target, ]
        edges <- edges[edges$regulator %in% top_regs, c("regulator", "target", "interaction")]

        edges_to_keep <- rbind(
          edges_to_keep,
          data.frame(
            regulator = as.character(edges$regulator),
            target = as.character(edges$target),
            interaction = as.character(edges$interaction)
          )
        )


        if (nrow(edges_to_keep) == 1) {
          edges_to_keep <- rbind(edges_to_keep, data.frame(regulator = NA, target = NA))
        }
      }
    } else if (weight_column %in% colnames(network_table)) {
      for (target in tgs) {
        edges <- network_table[network_table$target == target, ]
        edges <- edges[order(edges[, weight_column], decreasing = TRUE), ]
        edges <- edges[1:regulators_num, c("regulator", "target", "interaction")]

        edges_to_keep <- rbind(edges_to_keep, data.frame(regulator = as.character(edges$regulator), target = as.character(edges$target), interaction = as.character(edges$interaction)))
      }
    } else {
      if (is.null(gene_ranks)) {
        stop("Need to supply gene_ranks.")
      } else {
        stop("invalid weight_column")
      }
    }

    ktgraph[[epoch]] <- edges_to_keep
  }


  agg <- dplyr::bind_rows(ktgraph)[, c("regulator", "target", "interaction")]
  agg <- agg[!duplicated(agg), ]
  agg_net <- igraph::graph_from_data_frame(agg, directed = FALSE)
  agg_net <- igraph::delete_vertices(agg_net, v = igraph::V(agg_net)$name[is.na(igraph::V(agg_net)$name)])

  if (fixed_layout) {
    layout <- igraph::layout_with_fr(agg_net)
    rownames(layout) <- igraph::V(agg_net)$name
  }


  for (i in 1:length(network)) {
    epoch <- names(network)[i]
    edges_to_keep <- ktgraph[[epoch]]


    to_add <- agg[!(paste(agg$regulator, agg$target) %in% paste(edges_to_keep$regulator, edges_to_keep$target)), ]
    to_add <- to_add[(paste(to_add$regulator, to_add$target)) %in% (paste(network[[epoch]]$regulator, network[[epoch]]$target)), ]
    edges_to_keep <- rbind(edges_to_keep, to_add)


    if (fixed_layout) {
      net <- igraph::graph_from_data_frame(edges_to_keep[, c("regulator", "target", "interaction")], directed = FALSE)

      addvtcs <- igraph::V(agg_net)$name[!(igraph::V(agg_net)$name %in% igraph::V(net)$name)]
      net <- igraph::add_vertices(net, length(addvtcs), attr = list(name = addvtcs))
    } else {
      net <- igraph::graph_from_data_frame(edges_to_keep[, c("regulator", "target", "interaction")], directed = FALSE)
    }


    net <- igraph::delete_vertices(net, v = igraph::V(net)$name[is.na(igraph::V(net)$name)])
    net <- igraph::delete_vertices(net, v = igraph::V(net)$name[igraph::V(net)$name == "NA"])


    if (show_expression) {
      expression_from <- mean_expression[mean_expression$epoch == strsplit(epoch, split = "..", fixed = TRUE)[[1]][1], ]
      expression_to <- mean_expression[mean_expression$epoch == strsplit(epoch, split = "..", fixed = TRUE)[[1]][2], ]
      if (strsplit(epoch, split = "..", fixed = TRUE)[[1]][1] == strsplit(epoch, split = "..", fixed = TRUE)[[1]][2]) {
        igraph::V(net)$expression <- expression_from$mean_expression[match(igraph::V(net)$name, expression_from$gene)]
      } else {
        igraph::V(net)$expression <- ifelse(igraph::V(net)$name %in% edges_to_keep$regulator, expression_from$mean_expression[match(igraph::V(net)$name, expression_from$gene)], expression_to$mean_expression[match(igraph::V(net)$name, expression_to$gene)])
      }
    }


    if (fixed_layout) {
      layout_ordered <- layout[igraph::V(net)$name, ]
      network_plot_data <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    } else {
      layout <- igraph::layout_with_fr(net)
      rownames(layout) <- igraph::V(net)$name
      layout_ordered <- layout[igraph::V(net)$name, ]
      network_plot_data <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    }


    network_plot_data$type <- "regulator"
    network_plot_data$type[network_plot_data$name %in% targets] <- "known target"
    network_plot_data$type <- factor(network_plot_data$type, levels = c("regulator", "known target"))

    network_plot_data <- network_plot_data[!(is.na(network_plot_data$name)), ]

    cols <- c("Activation" = "blue", "Repression" = "red")

    plot_list[[i]] <- ggplot() +
      geom_edges(
        data = network_plot_data,
        aes(x = x, y = y, xend = xend, yend = yend, color = interaction),
        size = 0.75,
        curvature = 0.1,
        alpha = 0.6
      )

    if (show_expression) {
      plot_list[[i]] <- plot_list[[i]] +
        ggnetwork::geom_nodes(
          data = network_plot_data,
          aes(x = x, y = y, xend = xend, yend = yend, shape = type, size = expression, alpha = expression),
          color = "black"
        )
    } else {
      plot_list[[i]] <- plot_list[[i]] +
        ggnetwork::geom_nodes(
          data = network_plot_data,
          aes(x = x, y = y, xend = xend, yend = yend, shape = type),
          color = "darkgray",
          size = node_size,
          alpha = 0.8
        ) +
        ggnetwork::geom_nodes(
          data = network_plot_data[network_plot_data$type == "regulator", ],
          aes(x = x, y = y, xend = xend, yend = yend, shape = type),
          color = "#8C4985",
          size = node_size,
          alpha = 0.8
        )
    }

    plot_list[[i]] <- plot_list[[i]] +
      scale_color_manual(values = cols) +
      theme_blank() +
      ggtitle(names(network)[i])

    if (declutter) {
      keep <- union(edges_to_keep$regulator, edges_to_keep$target)
      plot_list[[i]] <- plot_list[[i]] +
        geom_nodelabel_repel(
          data = network_plot_data[network_plot_data$name %in% keep, ],
          aes(x = x, y = y, label = name),
          size = label_size,
          color = "#5A8BAD"
        )
    } else {
      plot_list[[i]] <- plot_list[[i]] +
        geom_nodelabel_repel(
          data = network_plot_data,
          aes(x = x, y = y, label = name),
          size = label_size,
          color = "#5A8BAD"
        )
    }

    plot_list[[i]] <- plot_list[[i]] + theme(legend.key.size = unit(legend_size / 8, "lines"), legend.key.width = unit(legend_size / 8, "lines"), legend.title = element_text(size = legend_size), legend.text = element_text(size = legend_size), plot.title = element_text(size = title_size))

    common_legend <- get_legend(plot_list[[i]])

    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }


  plot_list[sapply(plot_list, is.null)] <- NULL

  plot_list$legend <- common_legend
  do.call(gridExtra::grid.arrange, plot_list)
}


get_legend <- function(plot_list) {
  tmp <- ggplot_gtable(ggplot_build(plot_list))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' plots results of findDynGenes
#'
#' @inheritParams hm_dyn_epoch
#' @param geneAnn geneAnn
#' @param dynRes dynRes
#'
#' @return heatmaps
#'
#' @export
hm_dyn_clust <- function(
  matrix,
  dynRes,
  geneAnn,
  row_cols,
  limits = c(0, 10),
  toScale = FALSE,
  fontsize_row = 4
) {
  meta_data <- dynRes$cells
  t1 <- meta_data$pseudotime
  names(t1) <- as.vector(meta_data$cell_name)
  grps <- as.vector(meta_data$group)
  names(grps) <- as.vector(meta_data$cell_name)

  ord1 <- sort(t1)

  matrix <- matrix[, names(ord1)]
  grps <- grps[names(ord1)]

  genes <- rownames(geneAnn)


  missingGenes <- setdiff(genes, rownames(matrix))
  if (length(missingGenes) > 0) {
    cat("Missing genes: ", paste0(missingGenes, collapse = ","), "\n")
    genes <- intersect(genes, rownames(matrix))
  }

  value <- matrix[genes, ]
  if (toScale) {
    if (class(value)[1] != "matrix") {
      value <- t(scale(Matrix::t(value)))
    } else {
      value <- t(scale(t(value)))
    }
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]
  groupNames <- unique(grps)
  cells <- names(grps)

  xcol <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames
  names(row_cols) <- unique(as.vector(geneAnn$epoch))
  anno_colors <- list(group = xcol, epoch = row_cols)
  xx <- data.frame(group = as.factor(grps))
  rownames(xx) <- cells

  geneX <- as.data.frame(geneAnn[, "epoch"])
  rownames(geneX) <- rownames(geneAnn)
  colnames(geneX) <- "epoch"

  val_col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(25)

  pheatmap::pheatmap(
    value,
    cluster_cols = FALSE, cluster_rows = FALSE, color = val_col,
    show_colnames = FALSE, annotation_row = geneX,
    annotation_col = xx,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    annotation_colors = anno_colors,
    fontsize_row = fontsize_row
  )
}

#' plots results of findDynGenes
#'
#' @param matrix expression matrix
#' @param dynRes result of running findDynGenes
#' @param cluster cluster
#' @param topX topX
#' @param cRow cRow
#' @param cCol cCol
#' @param limits limits
#' @param toScale toScale
#' @param fontsize_row fontsize_row
#' @param geneAnn geneAnn
#' @param anno_colors anno_colors
#' @param show_rownames show_rownames
#' @param filename filename
#' @param width width
#' @param height height
#'
#' @return heatmap list
#'
#' @export
hm_dyn <- function(
  matrix,
  dynRes,
  cluster = "cluster",
  topX = 25,
  cRow = FALSE,
  cCol = FALSE,
  limits = c(0, 10),
  toScale = FALSE,
  fontsize_row = 4,
  geneAnn = FALSE,
  anno_colors = NULL,
  show_rownames = TRUE,
  filename = NA,
  width = NA,
  height = NA
) {
  colors <- NULL
  if (!is.null(anno_colors)) {
    colors <- anno_colors
  }

  meta_data <- dynRes$meta_data
  t1 <- meta_data$pseudotime
  names(t1) <- as.vector(meta_data$cells)
  grps <- as.vector(meta_data[, cluster])
  names(grps) <- as.vector(meta_data$cells)

  ord1 <- sort(t1)

  matrix <- matrix[, names(ord1)]
  grps <- grps[names(ord1)]

  genes <- dynRes$genes[1:topX]


  missingGenes <- setdiff(genes, rownames(matrix))
  if (length(missingGenes) > 0) {
    cat("Missing genes: ", paste0(missingGenes, collapse = ","), "\n")
    genes <- intersect(genes, rownames(matrix))
  }


  peakTime <- apply(matrix[genes, ], 1, which.max)
  genes_ordered <- names(sort(peakTime))

  value <- matrix[genes_ordered, ]
  if (toScale) {
    if (class(value)[1] != "matrix") {
      value <- t(scale(Matrix::t(value)))
    } else {
      value <- t(scale(t(value)))
    }
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]
  groupNames <- unique(grps)
  cells <- names(grps)

  xcol <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 11, name = "Paired"))
  )(length(groupNames))
  names(xcol) <- groupNames
  anno_colors <- list(
    group = xcol,
    pseudotime = viridis::magma(length(ord1) / 2, direction = -1)
  )
  xx <- data.frame(group = as.factor(grps), pseudotime = ord1)

  if (!is.null(meta_data$epoch)) {
    epochs <- as.vector(meta_data$epoch)
    names(epochs) <- as.vector(meta_data$cells)
    epochs <- epochs[names(ord1)]
    xx <- cbind(xx, data.frame(epoch = epochs))
  }

  rownames(xx) <- cells
  val_col <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
  )(25)
  if (is.data.frame(geneAnn)) {
    cat("right spot\n")
    if (is.null(colors)) {
      pheatmap::pheatmap(
        value,
        cluster_rows = cRow,
        cluster_cols = cCol,
        color = val_col,
        show_colnames = FALSE,
        annotation_row = geneAnn,
        show_rownames = show_rownames,
        annotation_col = xx,
        annotation_names_col = FALSE,
        annotation_colors = anno_colors,
        fontsize_row = fontsize_row,
        filename = filename,
        width = width,
        height = height
      )
    } else {
      pheatmap::pheatmap(
        value,
        cluster_rows = cRow,
        cluster_cols = cCol,
        color = val_col,
        show_colnames = FALSE,
        annotation_row = geneAnn,
        show_rownames = show_rownames,
        annotation_col = xx,
        annotation_names_col = FALSE,
        annotation_colors = colors,
        fontsize_row = fontsize_row,
        filename = filename,
        width = width,
        height = height
      )
    }
  } else {
    if (is.null(colors)) {
      pheatmap::pheatmap(
        value,
        cluster_rows = cRow,
        cluster_cols = cCol,
        color = val_col,
        show_colnames = FALSE,
        annotation_names_row = FALSE,
        show_rownames = show_rownames,
        annotation_col = xx,
        annotation_names_col = FALSE,
        annotation_colors = anno_colors,
        fontsize_row = fontsize_row,
        border_color = NA,
        filename = filename,
        width = width,
        height = height
      )
    } else {
      pheatmap::pheatmap(
        value,
        cluster_rows = cRow,
        cluster_cols = cCol,
        color = val_col,
        show_colnames = FALSE,
        annotation_names_row = FALSE,
        show_rownames = show_rownames,
        annotation_col = xx,
        annotation_names_col = FALSE,
        annotation_colors = colors,
        fontsize_row = fontsize_row,
        border_color = NA,
        filename = filename,
        width = width,
        height = height
      )
    }
  }
}


#' heatmap
#'
#' @param matrix expression matrix
#' @param epochRes result of running findDynGenes
#' @param row_cols row_cols
#' @param limits limits
#' @param toScale toScale
#' @param fontsize_row fontsize_row
#' '
#' @return heatmaps
#'
#' @export
hm_dyn_epoch <- function(
  matrix,
  epochRes,
  row_cols,
  limits = c(0, 10),
  toScale = FALSE,
  fontsize_row = 4
) {
  meta_data <- epochRes$cells
  t1 <- meta_data$pseudotime
  names(t1) <- as.vector(meta_data$cell_name)
  grps <- as.vector(meta_data$epoch)
  names(grps) <- as.vector(meta_data$cell_name)
  ord1 <- sort(t1)
  matrix <- matrix[, names(ord1)]
  grps <- grps[names(ord1)]


  geneTab <- epochRes$genes
  genes <- rownames(geneTab)


  missingGenes <- setdiff(genes, rownames(matrix))
  if (length(missingGenes) > 0) {
    cat("Missing genes: ", paste0(missingGenes, collapse = ","), "\n")
    genes <- intersect(genes, rownames(matrix))
  }


  genes_ordered <- rownames(geneTab)[order(geneTab$peakTime)]
  value <- matrix[genes_ordered, ]

  if (toScale) {
    if (class(value)[1] != "matrix") {
      value <- t(scale(Matrix::t(value)))
    } else {
      value <- t(scale(t(value)))
    }
  }
  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]
  groupNames <- unique(grps)
  cells <- names(grps)

  xcol <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Paired")))(length(groupNames))
  names(xcol) <- groupNames

  names(row_cols) <- unique(as.vector(geneTab$epoch))
  anno_colors <- list(group = xcol, epoch = row_cols)
  xx <- data.frame(group = as.factor(grps))
  rownames(xx) <- cells

  geneX <- as.data.frame(geneTab[, "epoch"])
  rownames(geneX) <- rownames(geneTab)
  colnames(geneX) <- "epoch"

  cat("OK\n")

  val_col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(25)

  neps <- length(unique(meta_data$epoch))

  c_gaps <- rep(0, neps - 1)
  epTally <- table(meta_data$epoch)
  c_gaps[1] <- epTally[1]
  for (i in seq(2, length(epTally) - 1)) {
    c_gaps[i] <- c_gaps[i - 1] + epTally[i]
    cat(c_gaps[i], "\n")
  }

  g_gaps <- length(unique(meta_data$epoch))
  epTally <- table(geneTab$epoch)
  g_gaps[1] <- epTally[1]
  for (i in seq(2, length(epTally) - 1)) {
    g_gaps[i] <- g_gaps[i - 1] + epTally[i]
    cat(g_gaps[i], "\n")
  }

  pheatmap::pheatmap(
    value,
    cluster_cols = FALSE, cluster_rows = FALSE, color = val_col,
    show_colnames = FALSE, annotation_row = geneX,
    annotation_col = xx,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    annotation_colors = anno_colors,
    fontsize_row = fontsize_row, gaps_row = g_gaps
  )
}
