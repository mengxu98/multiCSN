#' @include setClass.R
#' @include setGenerics.R

#' @title Plot dynamic networks
#'
#' @param network_table The weight data table of network.
#' @param regulators A character vector of regulators to include.
#' @param targets A character vector of targets to include.
#' @param legend_position The position of legend.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' data(example_matrix, package = "inferCSN")
#' network_table <- inferCSN::inferCSN(example_matrix)
#' plot_static_networks(
#'   network_table,
#'   regulators = "g1"
#' )
#' plot_static_networks(
#'   network_table,
#'   targets = "g1"
#' )
#' plot_static_networks(
#'   network_table,
#'   regulators = "g2",
#'   targets = "g3"
#' )
plot_static_networks <- function(
  network_table,
  regulators = NULL,
  targets = NULL,
  legend_position = "right"
) {
  network_table <- inferCSN::network_format(
    network_table,
    regulators = regulators,
    targets = targets
  )

  net <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "weight", "Interaction")],
    directed = FALSE
  )

  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name, ]
  regulator_network <- ggnetwork(
    net,
    layout = layout_ordered,
    cell.jitter = 0
  )

  regulator_network$is_regulator <- as.character(
    regulator_network$name %in% regulators
  )
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  g <- ggplot() +
    geom_edges(
      data = regulator_network,
      aes(
        x = x, y = y,
        xend = xend, yend = yend,
        size = weight,
        color = Interaction
      ),
      size = 0.75,
      curvature = 0.1,
      alpha = .6
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y),
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y),
      color = "#8C4985",
      size = 6,
      alpha = .8
    ) +
    scale_color_manual(values = cols) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y, label = name),
      size = 2,
      color = "#5A8BAD"
    ) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y, label = name),
      size = 3.5,
      color = "black"
    ) +
    theme_blank() +
    theme(legend.position = legend_position)

  return(g)
}

#' @title Plot contrast networks
#'
#' @md
#' @param network_table The weight data table of network.
#' @param degree_value Degree value to filter nodes.
#' Default is `0`.
#' @param weight_value Weight value to filter edges.
#' Default is `0`.
#' @param legend_position The position of legend.
#' Default is `"bottom"`.
#'
#' @return
#' A ggplot2 object.
#' @export
#'
#' @examples
#' data(example_matrix, package = "inferCSN")
#' network_table <- inferCSN::inferCSN(example_matrix)
#' plot_contrast_networks(network_table[1:50, ])
plot_contrast_networks <- function(
  network_table,
  degree_value = 0,
  weight_value = 0,
  legend_position = "bottom"
) {
  graph <- inferCSN::network_format(network_table) |>
    tidygraph::as_tbl_graph() |>
    dplyr::mutate(
      degree = tidygraph::centrality_degree(mode = "out")
    ) |>
    dplyr::filter(degree > degree_value) |>
    tidygraph::activate(edges)

  g <- ggraph(graph, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      aes(
        colour = Interaction,
        filter = weight > weight_value,
        edge_width = weight
      ),
      arrow = arrow(length = unit(3, "mm")),
      start_cap = square(3, "mm"),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_width(range = c(0, 1)) +
    facet_edges(~Interaction) +
    geom_node_point(aes(size = degree), colour = "#A1B7CE") +
    geom_node_text(aes(label = name), repel = TRUE) +
    coord_fixed() +
    theme_graph(
      base_family = "serif",
      foreground = "steelblue",
      fg_text_colour = "white"
    ) +
    theme(legend.position = legend_position)

  return(g)
}

#' @title Plot dynamic networks
#'
#' @param object Network table (data.frame with regulator, target, weight, celltype) or CSNObject.
#' @param ... Arguments passed to methods.
#' @rdname plot_dynamic_networks
#' @export
setGeneric(
  name = "plot_dynamic_networks",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_dynamic_networks")
  }
)

#' @param object CSNObject with dynamic network.
#' @param network Name of the network; default uses \code{DefaultNetwork(object)}.
#' @param celltypes Character vector of states/cell types; \code{NULL} for all.
#' @param weight_cutoff Numeric threshold for edge weight filtering.
#' @rdname plot_dynamic_networks
#' @export
setMethod(
  "plot_dynamic_networks",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           weight_cutoff = NULL,
           celltypes_order = NULL,
           ntop = 10,
           filter_by_ntop = TRUE,
           ...) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "dynamic",
      verbose = TRUE,
      caller = "plot_dynamic_networks"
    )
    nets <- GetNetwork(object, network = network)
    if (is.null(nets) || !is.list(nets)) {
      stop("No dynamic network found. Ensure object has a network with multiple states.")
    }
    if (is.null(celltypes)) celltypes <- names(nets)
    celltypes <- intersect(celltypes, names(nets))
    if (length(celltypes) == 0) stop("No matching cell types/states in network.")
    dyn_list <- lapply(celltypes, function(sid) {
      net <- nets[[sid]]
      df <- tryCatch(export_csn(net, weight_cutoff = weight_cutoff), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      df <- df[, c("regulator", "target", "weight"), drop = FALSE]
      df$celltype <- sid
      df
    })
    dyn_list <- dyn_list[!sapply(dyn_list, is.null)]
    if (length(dyn_list) == 0) stop("No network edges to plot.")
    dyn_table <- dplyr::bind_rows(dyn_list)
    if (is.null(celltypes_order)) celltypes_order <- celltypes
    plot_dynamic_networks(dyn_table, celltypes_order = celltypes_order, ntop = ntop, filter_by_ntop = filter_by_ntop, ...)
  }
)

#' @rdname plot_dynamic_networks
#' @export
setMethod(
  "plot_dynamic_networks",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)

#' @param object Data.frame with columns regulator, target, weight, celltype.
#' @param celltypes_order The order of cell types.
#' @param ntop The number of top regulators to show (labels and edge filtering).
#' When \code{filter_by_ntop=TRUE}, only edges from top \code{ntop} regulators are plotted to avoid huge networks.
#' Default is `10`.
#' @param filter_by_ntop Logical. If \code{TRUE} (default), filter edges to top \code{ntop} regulators per celltype.
#' @param title The title of figure.
#' Default is `NULL`.
#' @param theme_type The theme of figure.
#' Could be `"theme_void"`, `"theme_blank"`, or `"theme_facet"`.
#' Default is `"theme_void"`.
#' @param plot_type The type of figure.
#' Could be `"ggplot"`, `"animate"`, or `"ggplotly"`.
#' Default is `"ggplot"`.
#' @param layout The layout of figure.
#' Could be `"fruchtermanreingold"` or `"kamadakawai"`.
#' Default is `"fruchtermanreingold"`.
#' @param nrow The number of rows of figure (for \code{facet_wrap} when \code{combine=TRUE}).
#' Default is `2`.
#' @param ncol Integer or \code{NULL}. Number of columns for \code{facet_wrap}.
#' @param byrow Logical. If \code{TRUE}, fill facets by row (\code{as.table=TRUE} in \code{facet_wrap}).
#' @param combine Logical. If \code{TRUE} (default), return a single faceted plot; if \code{FALSE}, return a named list of plots per celltype.
#' @param figure_save Whether to save the figure file.
#' Default is `FALSE`.
#' @param figure_name The name of figure file.
#' Default is `NULL`.
#' @param figure_width The width of figure.
#' Default is `6`.
#' @param figure_height The height of figure.
#' Default is `6`.
#' @param seed The seed random use to plot network.
#' Default is `1`.
#'
#' @return
#' A dynamic figure object.
#' @export
#'
#' @examples
#' data(example_matrix, package = "inferCSN")
#' network <- inferCSN::inferCSN(example_matrix, verbose = FALSE)
#' network$celltype <- rep(
#'   c("cluster5", "cluster3", "cluster2", "cluster1", "cluster6"),
#'   length.out = nrow(network)
#' )
#'
#' celltypes_order <- c(
#'   "cluster5", "cluster3",
#'   "cluster2", "cluster1",
#'   "cluster6"
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order[1:3]
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order,
#'   plot_type = "ggplotly"
#' )
#' @rdname plot_dynamic_networks
#' @export
setMethod(
  "plot_dynamic_networks",
  signature(object = "data.frame"),
  function(object,
           celltypes_order = NULL,
           ntop = 10,
           filter_by_ntop = TRUE,
           title = NULL,
           theme_type = "theme_void",
           plot_type = "ggplot",
           layout = "fruchtermanreingold",
           nrow = 2,
           ncol = NULL,
           byrow = TRUE,
           combine = TRUE,
           figure_save = FALSE,
           figure_name = NULL,
           figure_width = 6,
           figure_height = 6,
           seed = 1) {
    network_table <- as.data.frame(object)
    if (ncol(network_table) < 4) stop("object must have columns: regulator, target, weight, celltype")
    names(network_table)[1:4] <- c("regulator", "target", "weight", "celltype")
    if (is.null(celltypes_order)) celltypes_order <- unique(network_table$celltype)
    network_table$regulator <- as.character(network_table$regulator)
    network_table$target <- as.character(network_table$target)
    celltypes_list <- unique(intersect(celltypes_order, network_table$celltype))


    if (filter_by_ntop) {
      network_table <- purrr::map_dfr(
        celltypes_list,
        .f = function(x) {
          sub <- network_table[which(network_table$celltype == x), ]
          if (nrow(sub) == 0) {
            return(sub)
          }
          top_regs <- sub |>
            dplyr::group_by(regulator) |>
            dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
            dplyr::arrange(dplyr::desc(n)) |>
            utils::head(ntop) |>
            dplyr::pull(regulator)
          sub[sub$regulator %in% top_regs, , drop = FALSE]
        }
      )
    } else {
      network_table <- purrr::map_dfr(
        celltypes_list,
        .f = function(x) network_table[which(network_table$celltype == x), ]
      )
    }

    nodes <- unique(
      c(network_table$regulator, network_table$target)
    )
    dnodes <- data.frame(id = 1:length(nodes), label = nodes)
    edges <- dplyr::left_join(
      network_table,
      dnodes,
      by = c("regulator" = "label")
    ) |>
      dplyr::rename(from = id) |>
      dplyr::left_join(
        dnodes,
        by = c("target" = "label")
      ) |>
      dplyr::rename(to = id) |>
      dplyr::select(from, to, weight, celltype)
    edges$Interaction <- ifelse(
      edges$weight > 0, "Activation", "Repression"
    )
    edges$weight <- abs(edges$weight)


    edges <- edges[edges$from != edges$to, , drop = FALSE]
    if (nrow(edges) == 0) {
      stop("No edges remaining after removing self-loops. Check network_table.")
    }


    dedges <- edges |>
      dplyr::group_by(from, to, celltype) |>
      dplyr::summarise(
        weight = mean(weight, na.rm = TRUE),
        Interaction = dplyr::first(Interaction),
        .groups = "drop"
      )
    dnodes$label <- gsub("\\.", "-", dnodes$label)
    network_data <- network::network(
      dedges,
      vertex.attr = dnodes,
      matrix.type = "edgelist",
      ignore.eval = FALSE,
      directed = TRUE,
      multiple = TRUE
    )

    set.seed(seed)
    layout <- match.arg(
      layout,
      c("fruchtermanreingold", "kamadakawai")
    )
    ggnetwork_data <- ggnetwork(
      network_data,
      arrow.size = 0.1,
      arrow.gap = 0.015,
      by = "celltype",
      weights = "weight",
      layout = layout
    )

    nodes_data <- purrr::map_dfr(
      celltypes_list,
      .f = function(x) {
        nodes_data_celltype <- network_table[which(network_table$celltype == x), ] |>
          dplyr::group_by(
            regulator
          ) |>
          dplyr::summarise(
            targets_num = dplyr::n()
          ) |>
          dplyr::arrange(
            dplyr::desc(targets_num)
          ) |>
          as.data.frame()
        nodes_data_celltype$label_genes <- as.character(
          nodes_data_celltype$regulator
        )
        if (nrow(nodes_data_celltype) > ntop) {
          cf <- nodes_data_celltype$targets_num[ntop]
          nodes_data_celltype$label_genes[which(nodes_data_celltype$targets_num < cf)] <- ""
        } else if (nrow(nodes_data_celltype) == 0) {
          return()
        }
        nodes_data_celltype$celltype <- x

        return(nodes_data_celltype)
      }
    )

    names(nodes_data)[1] <- "label"
    ggnetwork_data <- merge(
      ggnetwork_data,
      nodes_data,
      by = c("label", "celltype"),
      all.x = T
    )
    ggnetwork_data$targets_num[which(is.na(ggnetwork_data$targets_num))] <- 0
    ggnetwork_data$label_genes[which(is.na(ggnetwork_data$label_genes))] <- ""
    ggnetwork_data$celltype <- factor(
      ggnetwork_data$celltype,
      levels = celltypes_order
    )
    cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

    plot_type <- match.arg(
      plot_type,
      c("ggplot", "animate", "ggplotly")
    )
    p <- ggplot(ggnetwork_data, aes(x, y, xend = xend, yend = yend))
    if (plot_type == "ggplotly") {
      p <- p + geom_edges(
        aes(color = Interaction),
        linewidth = 0.7,
        arrow = arrow(length = unit(3, "pt"), type = "closed")
      )
    } else {
      p <- p + geom_edges(
        aes(color = Interaction, alpha = weight),
        linewidth = 0.7,
        arrow = arrow(length = unit(3, "pt"), type = "closed")
      )
    }
    p <- p +
      geom_nodes(
        aes(size = targets_num),
        color = "darkgray",
        alpha = 0.9
      ) +
      geom_nodetext(
        aes(label = label_genes),
        color = "black"
      ) +
      theme(aspect.ratio = 2, legend.position = "bottom") +
      scale_color_manual(values = cols)

    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }

    theme_type <- match.arg(
      theme_type,
      c("theme_void", "theme_facet", "theme_blank")
    )
    p <- switch(
      EXPR = theme_type,
      "theme_void" = p + theme_void(),
      "theme_blank" = p + theme_blank(),
      "theme_facet" = p + theme_facet()
    )
    if (plot_type == "ggplot") {
      if (combine) {
        p <- p +
          facet_wrap(~celltype, nrow = nrow, ncol = ncol, as.table = byrow)
        if (figure_save) {
          if (is.null(figure_name)) {
            figure_name <- "networks.pdf"
          }
          ggsave(
            figure_name,
            p,
            width = figure_width,
            height = figure_height
          )
        }
      } else {
        plot_list <- lapply(celltypes_list, function(ct) {
          p_sub <- p %+% ggnetwork_data[ggnetwork_data$celltype == ct, ]
          p_sub + ggtitle(ct)
        })
        names(plot_list) <- celltypes_list
        return(plot_list)
      }
    }

    if (plot_type == "animate") {
      p <- p + gganimate::transition_states(states = celltype)
      p <- gganimate::animate(
        p,
        render = gganimate::gifski_renderer()
      )
      if (figure_save) {
        if (is.null(figure_name)) {
          figure_name <- "networks.gif"
        }
        gganimate::anim_save(figure_name, animation = p)
      }
    }

    if (plot_type == "ggplotly") {
      if (combine) {
        p <- p +
          facet_wrap(~celltype, nrow = nrow, ncol = ncol, as.table = byrow)
        p <- plotly::ggplotly(p)
      } else {
        plot_list <- lapply(celltypes_list, function(ct) {
          p_sub <- p %+% ggnetwork_data[ggnetwork_data$celltype == ct, ]
          plotly::ggplotly(p_sub + ggtitle(ct))
        })
        names(plot_list) <- celltypes_list
        return(plot_list)
      }
    }

    return(p)
  }
)
