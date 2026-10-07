#' @include setClass.R
#' @include setGenerics.R

#' @title Plot static networks
#'
#' @param network_table The weight data table of network.
#' @param regulators A character vector of regulators to include.
#' @param targets A character vector of targets to include.
#' @param legend_position The position of legend. Default is \code{"right"}.
#' @param theme_use Theme function name or choice: \code{"theme_this"},
#' \code{"theme_blank"}, or \code{"theme_void"}.
#' Default is \code{"theme_this"}.
#' @param cols Optional named color vector for Interaction (\code{Activation} and \code{Repression}).
#' Default is \code{c("Activation" = "#1F78B4", "Repression" = "#E31A1C")}.
#' @param palette Palette name used by \code{thisplot::palette_colors()}. Default is \code{"Chinese"}.
#' @param palcolor Optional named color vector overriding default node colors.
#' @param aspect.ratio Aspect ratio of panel. Default is \code{1}.
#' @param label.size Label text size. Default is \code{3.5}.
#' @param label.fg Label foreground text color. Default is \code{"black"}.
#' @param label.bg Label outline/background color. Default is \code{"white"}.
#' @param label.bg.r Label outline radius. Default is \code{0.12}.
#' @param tf_nodes Optional character vector of focal TF nodes. When supplied,
#' these nodes are drawn as diamonds and other nodes as circles.
#' @param node_state Optional node-state mapping. Supply either a named character
#' vector or a data.frame with columns \code{name} and \code{state}.
#' @param state_colors Optional named vector of colors for \code{node_state}.
#' @param node_size Optional named numeric vector of node sizes.
#' @param label_nodes Optional character vector of nodes to label. Default labels all nodes.
#' @param seed Optional integer seed for the force-directed layout.
#' @param edge_alpha Edge alpha. Default is \code{0.7}.
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
  legend_position = "right",
  theme_use = c("theme_this", "theme_blank", "theme_void"),
  cols = NULL,
  palette = "Chinese",
  palcolor = NULL,
  aspect.ratio = 1,
  label.size = 3.5,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.12,
  tf_nodes = NULL,
  node_state = NULL,
  state_colors = NULL,
  node_size = NULL,
  label_nodes = NULL,
  seed = NULL,
  edge_alpha = 0.7
) {
  theme_use <- if (is.character(theme_use)) theme_use[1] else "theme_this"
  network_table <- inferCSN::network_format(
    network_table,
    regulators = regulators,
    targets = targets
  )

  use_tf_shapes <- !is.null(tf_nodes)
  if (is.null(tf_nodes)) {
    tf_nodes <- regulators
  }
  tf_nodes <- intersect(as.character(tf_nodes), network_table$regulator)
  tf_nodes <- unique(c(tf_nodes, setdiff(regulators, tf_nodes)))

  state_map <- NULL
  if (!is.null(node_state)) {
    if (is.data.frame(node_state)) {
      if (!all(c("name", "state") %in% names(node_state))) {
        stop(
          "`node_state` data.frame must contain `name` and `state` columns.",
          call. = FALSE
        )
      }
      state_map <- stats::setNames(
        as.character(node_state$state),
        as.character(node_state$name)
      )
    } else if (is.character(node_state) && !is.null(names(node_state))) {
      state_map <- node_state
    } else {
      stop(
        "`node_state` must be a named character vector or a data.frame with `name` and `state` columns.",
        call. = FALSE
      )
    }
  }

  if (!is.null(node_size) && is.null(names(node_size))) {
    stop("`node_size` must be a named numeric vector.", call. = FALSE)
  }

  net <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "weight", "Interaction")],
    directed = FALSE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }
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
  node_network <- regulator_network[
    !duplicated(regulator_network[, c("x", "y", "name")]), ,
    drop = FALSE
  ]
  label_names <- if (is.null(label_nodes)) {
    node_network$name
  } else {
    intersect(as.character(label_nodes), node_network$name)
  }
  label_network <- node_network[node_network$name %in% label_names, , drop = FALSE]

  if (is.null(cols)) {
    cols <- c("Activation" = "#1F78B4", "Repression" = "#E31A1C")
  }

  node_colors <- thisplot::palette_colors(
    c("Regulator", "Target"),
    palette = palette, palcolor = palcolor
  )
  reg_col <- node_colors[["Regulator"]]
  tgt_col <- node_colors[["Target"]]

  hybrid <- use_tf_shapes || !is.null(state_map)

  add_static_labels <- function(plot, data, size, fontface) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot + ggrepel::geom_text_repel(
        data = data,
        aes(x = x, y = y, label = name),
        size = size,
        fontface = fontface,
        color = label.fg,
        bg.color = label.bg,
        bg.r = label.bg.r,
        segment.colour = NA,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    } else {
      plot + geom_nodetext(
        data = data,
        aes(x = x, y = y, label = name),
        size = size,
        fontface = fontface,
        color = label.fg
      )
    }
  }

  g <- ggplot() +
    geom_edges(
      data = regulator_network,
      aes(
        x = x, y = y,
        xend = xend, yend = yend,
        linewidth = weight,
        color = Interaction
      ),
      linewidth = 0.65,
      curvature = 0.08,
      alpha = edge_alpha
    )

  if (hybrid) {
    node_rows <- node_network
    node_rows$is_tf <- node_rows$name %in% tf_nodes
    node_rows$node_type <- ifelse(node_rows$is_tf, "Focal TF", "Target")
    node_rows$node_size <- ifelse(node_rows$is_tf, 6.4, 3.0)
    if (!is.null(node_size)) {
      mapped_size <- unname(node_size[node_rows$name])
      use_size <- !is.na(mapped_size)
      node_rows$node_size[use_size] <- mapped_size[use_size]
    }
    if (!is.null(state_map)) {
      node_rows$peak_state <- unname(state_map[node_rows$name])
      state_levels <- if (!is.null(state_colors) && !is.null(names(state_colors))) {
        names(state_colors)
      } else {
        unique(stats::na.omit(node_rows$peak_state))
      }
      node_rows$peak_state <- factor(node_rows$peak_state, levels = state_levels)
      if (is.null(state_colors)) {
        state_colors <- thisplot::palette_colors(
          state_levels,
          palette = palette,
          palcolor = palcolor
        )
        if (is.null(names(state_colors))) {
          names(state_colors) <- state_levels
        }
      } else if (is.null(names(state_colors))) {
        names(state_colors) <- state_levels
      }
      group_col <- "Peak state"
      node_rows[[group_col]] <- as.character(node_rows$peak_state)
      group_colors <- state_colors
    } else {
      group_col <- "Node type"
      node_rows[[group_col]] <- node_rows$node_type
      group_colors <- c("Focal TF" = reg_col, "Target" = tgt_col)
    }
    if (!use_tf_shapes) {
      node_rows$node_type <- "Node"
      group_colors <- c("Node" = reg_col)
    }
    node_rows$label_size <- ifelse(
      node_rows$is_tf, label.size, max(1, label.size - 0.7)
    )
    node_rows$label_face <- ifelse(node_rows$is_tf, "bold", "plain")
    node_rows <- node_rows[order(node_rows$is_tf), , drop = FALSE]
    return(thisplot::NetworkPlot(
      edge = as.data.frame(network_table[, c("regulator", "target", "weight", "Interaction")]),
      node = as.data.frame(node_rows[, c(
        "name", group_col, "node_type", "node_size", "label_size", "label_face"
      )]),
      from = "regulator", to = "target", weight = "weight",
      edge_group = "Interaction",
      edge_palcolor = cols,
      node_group = group_col,
      node_palcolor = group_colors,
      node_shape = "node_type",
      node_shape_values = c("Focal TF" = 23, "Target" = 21, "Node" = 21),
      node_size = "node_size",
      layout = "fr",
      seed = seed,
      edge_curvature = 0.08,
      edge_alpha = edge_alpha,
      label = TRUE,
      label_nodes = label_names,
      label.size = label.size,
      label_size = "label_size",
      label_face = "label_face",
      label.fg = label.fg,
      label.bg = label.bg,
      label.bg.r = label.bg.r,
      aspect.ratio = aspect.ratio,
      legend.position = legend_position,
      theme_use = theme_use
    ))
  } else {
    g <- g +
      geom_nodes(
        data = regulator_network[regulator_network$is_regulator == "FALSE", ],
        aes(x = x, y = y),
        color = "grey50",
        fill = tgt_col,
        shape = 21,
        stroke = 0.35,
        size = 3.5,
        alpha = 0.8
      ) +
      geom_nodes(
        data = regulator_network[regulator_network$is_regulator == "TRUE", ],
        aes(x = x, y = y),
        color = "black",
        fill = reg_col,
        shape = 21,
        stroke = 0.5,
        size = 6,
        alpha = 0.95
      ) +
      scale_color_manual(values = cols)
    g <- add_static_labels(g, label_network, label.size, "bold")
  }

  if (identical(theme_use, "theme_this")) {
    g <- g + thisplot::theme_this(aspect.ratio = aspect.ratio) +
      theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = legend_position
      )
  } else if (identical(theme_use, "theme_blank")) {
    g <- g + thisplot::theme_blank(add_coord = FALSE) +
      theme(aspect.ratio = aspect.ratio, legend.position = legend_position)
  } else {
    g <- g + theme_void() +
      theme(aspect.ratio = aspect.ratio, legend.position = legend_position)
  }
  g <- g + theme(panel.border = element_blank())

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
#' @param cols Named vector of colors for edge interactions. Default uses \code{c("Activation" = "#1F78B4", "Repression" = "#E31A1C")}.
#' @param regulator_color Color of regulator nodes. Default is \code{"#8C4985"}.
#' @param target_color Color of target nodes. Default is \code{"#B0C4DE"}.
#' @param label.size Node label font size. Default is \code{3.5}.
#' @param base_family Font family for graph theme. Default is \code{""}.
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
  cols = NULL,
  regulator_color = "#8C4985",
  target_color = "#B0C4DE",
  label.size = 3.5,
  base_family = "",
  legend_position = "bottom"
) {
  if (is.null(cols)) {
    cols <- c("Activation" = "#1F78B4", "Repression" = "#E31A1C")
  }

  formatted <- inferCSN::network_format(network_table)
  formatted$abs_weight <- abs(as.numeric(formatted$weight))
  graph <- formatted |>
    tidygraph::as_tbl_graph() |>
    dplyr::mutate(
      degree = tidygraph::centrality_degree(mode = "out"),
      is_regulator = name %in% formatted$regulator
    )
  graph <- graph |> dplyr::filter(degree > degree_value)
  graph <- graph |>
    tidygraph::activate(edges) |>
    dplyr::filter(weight > weight_value)

  g <- ggraph(graph, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      aes(
        edge_colour = Interaction,
        edge_width = abs_weight,
        edge_alpha = abs_weight
      ),
      arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_width(range = c(0.4, 1.3)) +
    scale_edge_alpha(range = c(0.4, 0.9)) +
    scale_edge_colour_manual(values = cols) +
    geom_node_point(aes(size = degree, color = is_regulator)) +
    scale_color_manual(
      values = c(`TRUE` = regulator_color, `FALSE` = target_color),
      labels = c(`TRUE` = "Regulator", `FALSE` = "Target"),
      name = "Node type"
    ) +
    geom_node_text(aes(label = name), repel = TRUE, size = label.size, fontface = "bold") +
    facet_edges(~Interaction) +
    coord_fixed() +
    guides(
      edge_colour = guide_legend(title = "Interaction", order = 1),
      color = guide_legend(title = "Node type", order = 2),
      size = guide_legend(title = "Degree", order = 3),
      edge_width = guide_legend(title = "Edge weight", order = 4),
      edge_alpha = "none"
    ) +
    theme_void(base_family = base_family) +
    theme(
      plot.margin = margin(10, 20, 10, 20),
      strip.text = element_text(size = 12, face = "bold", margin = margin(b = 8)),
      legend.position = legend_position,
      legend.box = "vertical",
      legend.margin = margin(t = 8)
    )

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
    network <- multicsn_resolve_network(
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
    stop_csnobject_runtime()
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
#' @param theme_use Theme function name or choice: \code{"theme_this"}, \code{"theme_blank"}, or \code{"theme_void"}.
#' Default is \code{"theme_this"}.
#' @param theme_type Legacy theme choice: \code{"theme_void"},
#' \code{"theme_blank"}, or \code{"theme_facet"}.
#' @param plot_type The type of figure.
#' Could be `"ggplot"`, `"animate"`, or `"ggplotly"`.
#' Default is `"ggplot"`.
#' @param layout The layout of figure.
#' Could be `"anchor"` (temporally aligned coordinates across states),
#' `"fruchtermanreingold"` or `"kamadakawai"`.
#' Default is `"anchor"`.
#' @param align_nodes Logical. If \code{TRUE} (default), anchor node positions across all states
#' based on the union graph to ensure temporal coherence and avoid jumping nodes.
#' @param cols Optional named color vector for Interaction (\code{Activation} and \code{Repression}).
#' Default is \code{c("Activation" = "#1F78B4", "Repression" = "#E31A1C")}.
#' @param palette Palette name used by \code{thisplot::palette_colors()}. Default is \code{"Chinese"}.
#' @param palcolor Optional named color vector overriding default node colors.
#' @param regulator_color Fill color for regulator nodes. Default is \code{NULL} (resolved via \code{palette}).
#' @param target_color Fill color for target nodes. Default is \code{NULL} (resolved via \code{palette}).
#' @param aspect.ratio Numeric aspect ratio for panels. Default is \code{1}.
#' @param aspect_ratio Alias for \code{aspect.ratio}.
#' @param label.size Label text size. Default is \code{3.5}.
#' @param label.fg Label foreground text color. Default is \code{"black"}.
#' @param label.bg Label outline/background color. Default is \code{"white"}.
#' @param label.bg.r Label outline radius. Default is \code{0.12}.
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
           theme_use = "theme_this",
           theme_type = NULL,
           plot_type = "ggplot",
           layout = "anchor",
           align_nodes = TRUE,
           cols = NULL,
           palette = "Chinese",
           palcolor = NULL,
           regulator_color = NULL,
           target_color = NULL,
           aspect.ratio = 1,
           aspect_ratio = NULL,
           label.size = 3.5,
           label.fg = "black",
           label.bg = "white",
           label.bg.r = 0.12,
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

    if (!is.null(aspect_ratio)) {
      aspect.ratio <- aspect_ratio
    }
    if (!is.null(theme_type)) {
      theme_use <- theme_type
    }

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
    layout <- if (is.character(layout)) layout[1] else "anchor"
    if (!layout %in% c("anchor", "fruchtermanreingold", "kamadakawai")) {
      layout <- "anchor"
    }

    if (identical(layout, "anchor") || isTRUE(align_nodes)) {
      union_g <- igraph::graph_from_data_frame(
        dedges[, c("from", "to", "weight"), drop = FALSE],
        directed = FALSE,
        vertices = data.frame(name = dnodes$id)
      )
      if (identical(layout, "kamadakawai")) {
        layout_mat <- igraph::layout_with_kk(union_g)
      } else {
        layout_mat <- igraph::layout_with_fr(union_g)
      }
      layout_mat[, 1] <- scales::rescale(layout_mat[, 1], to = c(0, 1))
      layout_mat[, 2] <- scales::rescale(layout_mat[, 2], to = c(0, 1))
      rownames(layout_mat) <- as.character(dnodes$id)
      layout_param <- layout_mat
    } else {
      layout_param <- layout
    }

    ggnetwork_data <- ggnetwork(
      network_data,
      arrow.size = 0.1,
      arrow.gap = 0.015,
      by = "celltype",
      weights = "weight",
      layout = layout_param
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
      all.x = TRUE
    )
    ggnetwork_data$targets_num[which(is.na(ggnetwork_data$targets_num))] <- 0
    ggnetwork_data$label_genes[which(is.na(ggnetwork_data$label_genes))] <- ""
    ggnetwork_data$celltype <- factor(
      ggnetwork_data$celltype,
      levels = celltypes_order
    )
    ggnetwork_data$is_regulator <- as.character(
      ggnetwork_data$label %in% network_table$regulator
    )

    if (is.null(cols)) {
      cols <- c("Activation" = "#1F78B4", "Repression" = "#E31A1C")
    }

    node_colors <- thisplot::palette_colors(
      c("Regulator", "Target"),
      palette = palette, palcolor = palcolor
    )
    reg_col <- if (!is.null(regulator_color)) {
      regulator_color
    } else {
      node_colors[["Regulator"]]
    }
    tgt_col <- if (!is.null(target_color)) {
      target_color
    } else {
      node_colors[["Target"]]
    }

    plot_type <- match.arg(
      plot_type,
      c("ggplot", "animate", "ggplotly")
    )
    p <- ggplot(ggnetwork_data, aes(x, y, xend = xend, yend = yend))
    if (plot_type == "ggplotly") {
      p <- p + geom_edges(
        aes(color = Interaction),
        linewidth = 0.65,
        arrow = arrow(length = unit(3, "pt"), type = "closed")
      )
    } else {
      p <- p + geom_edges(
        aes(color = Interaction, alpha = weight),
        linewidth = 0.65,
        arrow = arrow(length = unit(3.5, "pt"), type = "closed")
      )
    }

    node_df_fun <- function(is_reg) {
      function(df) {
        sub_df <- df[df$is_regulator == as.character(is_reg), , drop = FALSE]
        sub_df[!duplicated(sub_df[, c("x", "y", "label")]), , drop = FALSE]
      }
    }
    label_df_fun <- function(df) {
      sub_df <- df[nzchar(as.character(df$label_genes)), , drop = FALSE]
      if (nrow(sub_df) == 0) {
        return(sub_df)
      }
      sub_df[!duplicated(sub_df[, c("x", "y", "label_genes")]), , drop = FALSE]
    }

    p <- p +
      geom_nodes(
        data = node_df_fun(FALSE),
        aes(size = pmax(targets_num, 1)),
        color = "grey55",
        fill = tgt_col,
        shape = 21,
        stroke = 0.35,
        alpha = 0.75
      ) +
      geom_nodes(
        data = node_df_fun(TRUE),
        aes(size = pmax(targets_num, 3)),
        color = "grey15",
        fill = reg_col,
        shape = 21,
        stroke = 0.5,
        alpha = 0.95
      )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df_fun,
          aes(label = label_genes),
          color = label.fg,
          bg.color = label.bg,
          bg.r = label.bg.r,
          size = label.size,
          fontface = "bold",
          max.overlaps = 50,
          box.padding = 0.35,
          point.padding = 0.25
        )
    } else {
      p <- p +
        geom_nodetext(
          data = label_df_fun,
          aes(label = label_genes),
          color = label.fg,
          size = label.size,
          fontface = "bold"
        )
    }

    p <- p +
      scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
      scale_color_manual(values = cols)

    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }

    theme_use <- if (is.character(theme_use)) theme_use[1] else "theme_this"
    if (identical(theme_use, "theme_facet")) {
      p <- p + ggnetwork::theme_facet()
    } else if (!is.null(theme_type) && identical(theme_use, "theme_blank")) {
      p <- p + ggnetwork::theme_blank()
    } else if (identical(theme_use, "theme_this")) {
      p <- p + thisplot::theme_this(aspect.ratio = aspect.ratio) +
        theme(
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom"
        )
    } else if (identical(theme_use, "theme_blank")) {
      p <- p + thisplot::theme_blank(add_coord = FALSE) +
        theme(aspect.ratio = aspect.ratio, legend.position = "bottom")
    } else {
      p <- p + theme_void() +
        theme(aspect.ratio = aspect.ratio, legend.position = "bottom")
    }
    p <- p + theme(panel.border = element_blank())

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
