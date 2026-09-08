#' Function to score targets of effectors
#'
#' Adds columns for mean binding score, maximum binding score, and percent frequency target is a hit based on threshold
#'
#' @param aList list of dataframes containing binding data for each effector
#' @param threshold binding score threshold to compute hit frequency
#'
#' @return updated list of dataframes
#'
#' @export
score_targets <- function(aList, threshold = 50) {
  aList <- lapply(aList, function(x) {
    x$STRING <- NULL
    x
  })


  aList <- lapply(aList, function(x) {
    if (!is.null(x$Target_genes)) {
      rownames(x) <- x$Target_genes
      x$Target_genes <- NULL
    }
    x
  })


  aList <- lapply(aList, function(x) {
    x$mean_score <- x[, grepl("Average", colnames(x))]
    x
  })

  aList <- lapply(aList, function(x) {
    x <- x[, !grepl("Average", colnames(x))]
  })


  aList <- lapply(aList, function(x) {
    x$max_score <- apply(x, 1, max)
    x
  })


  aList <- lapply(aList, function(x) {
    columns <- colnames(x)[!(colnames(x) %in% c("Target_genes", "STRING", "mean_score", "max_score"))]
    n_pos <- apply(x[, columns], 1, function(y) {
      sum(y > threshold)
    })
    x$percent_freq <- n_pos / length(columns)
    x
  })

  aList
}


#' Finds binding targets given list of dataframes containing binding info for effectors
#'
#' @param aList the result of running score_targets
#' @param column column name used in either ranking or thresholding the possible targets
#' @param by_rank if TRUE will return top n_targets; if FALSE will use thresholds
#' @param n_targets the number of targets to return for each effector if by_rank=TRUE
#' @param threshold threshold to use if by_rank=FALSE
#'
#' @return list of binding targets for list of effectors
#'
#' @export
find_targets <- function(
  aList,
  column = "max_score",
  by_rank = FALSE,
  n_targets = 2000,
  threshold = NULL
) {
  if (by_rank) {
    aList <- lapply(aList, function(x) {
      x <- x[order(x[, column], decreasing = TRUE), ]
      n_targets <- min(n_targets, nrow(x))
      x[1:n_targets, ]
    })
  } else if (!is.null(threshold)) {
    aList <- lapply(aList, function(x) {
      x[(x[, column] > threshold), ]
    })
  } else {
    stop("Insufficient parameters.")
  }

  res <- lapply(aList, function(x) {
    rownames(x)
  })

  res
}


#' Computes mean expression of groups of genes
#'
#' @param matrix genes-by-cells expression matrix
#' @param module_list list containing grouped genes, each element is a module of genes
#'
#' @return module expression
#'
#' @export
mean_module_expression <- function(
  matrix,
  module_list
) {
  res <- data.frame(matrix(ncol = ncol(matrix), nrow = length(module_list)))
  colnames(res) <- colnames(matrix)
  rownames(res) <- names(module_list)
  for (m in names(module_list)) {
    genes <- intersect(module_list[[m]], rownames(matrix))
    exp_sub <- matrix[genes, ]
    if (length(genes) == 1) {
      res[m, ] <- matrix[genes, ]
    } else {
      res[m, ] <- colMeans(exp_sub)
    }
  }
  res
}

#' Computes mean expression of groups of genes in a dynamic network
#'
#' @param matrix genes-by-cells expression matrix
#' @param community_list list of dataframes with community assginemnts. Each dataframe is the result of running subnets.
#'
#' @return subnetwork expression
#'
#' @export
mean_subnetwork_expression <- function(
  matrix,
  community_list
) {
  subnets <- c()
  for (e in names(community_list)) {
    subnets <- c(subnets, paste(e, unique(as.character(community_list[[e]]$communities)), sep = "_"))
  }


  res <- data.frame(matrix(ncol = ncol(matrix), nrow = length(subnets)))
  colnames(res) <- colnames(matrix)
  rownames(res) <- subnets
  for (n in subnets) {
    e <- strsplit(n, "_")[[1]][1]
    c <- strsplit(n, "_")[[1]][2]

    genes <- community_list[[e]]$gene[as.character(community_list[[e]]$communities) == as.character(c)]
    exp_sub <- matrix[as.character(genes), ]

    if (length(genes) == 1) {
      res[n, ] <- matrix[as.character(genes), ]
    } else {
      res[n, ] <- colMeans(exp_sub)
    }
  }

  res
}

#' Function to assign nodes to communities via Louvain clustering
#'
#' @param df dataframe containing a static network
#' @param tfs TFs
#' @param tfonly if TRUE will limit network to TFs only
#'
#' @return dataframe containing community assignments for each gene
#'
#' @export
subnets <- function(
  df,
  tfs,
  tfonly = TRUE
) {
  if (tfonly) {
    df <- df[df$target %in% tfs, ]
  }

  df <- df[, c("regulator", "target", "weight", "interaction")]

  net <- igraph::graph_from_data_frame(df, directed = FALSE)

  c <- as.data.frame(
    as.table(igraph::membership(igraph::cluster_louvain(net)))
  )
  colnames(c) <- c("gene", "communities")

  c
}

#' @title rough_hierarchy
#'
#' @description
#'  returns rough roots in the network, rough roots selected as those connected to most number of nodes
#'
#' @param network_table The weight data table of network.
#' @param abs_weight Whether to use absolute weight values.
#' @param directed Whether the network is directed.
#'
#' @return list
#' @export
#'
#' @examples
#' data("example_matrix", package = "inferCSN")
#' network_table <- inferCSN::inferCSN(example_matrix)
#' rough_hierarchy(network_table)
rough_hierarchy <- function(
  network_table,
  abs_weight = TRUE,
  directed = TRUE
) {
  if (!igraph::is.igraph(network_table)) {
    if (abs_weight) {
      network_table$weight <- abs(network_table$weight)
    }
    network_table <- igraph::graph_from_data_frame(
      network_table,
      directed = directed
    )
  }

  distances <- data.frame(
    igraph::distances(network_table, mode = "out")
  )
  distances$num_connected <- rowSums(distances != Inf)

  roots <- rownames(distances)[distances$num_connected == max(distances$num_connected)]

  return(
    list(
      roots = roots,
      num_paths = max(distances$num_connected)
    )
  )
}


#' Function to return shortest path from 1 regulator to 1 target in a static network
#'
#' @param network_table The weight data table of network.
#' @param regulator The starting gene.
#' @param target The end gene.
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#' @param compare_to_average if TRUE will compute normalized against average path length
#'
#' @return shortest path, distance, normalized distance, and action
#'
#' @export
#'
#' @examples
#' data("example_matrix", package = "inferCSN")
#' network_table <- inferCSN::inferCSN(example_matrix)
#' static_shortest_path(
#'   network_table,
#'   regulator = "g1",
#'   target = "g2"
#' )
static_shortest_path <- function(
  network_table,
  regulator,
  target,
  weight_column = "weight",
  compare_to_average = FALSE
) {
  network_table$normalized_score <- thisutils::normalization(
    network_table[, weight_column],
    method = "maximum"
  )
  network_table$edge_length <- 1 - network_table$normalized_score

  ig <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "edge_length", "weight")],
    directed = TRUE
  )

  distance <- igraph::distances(
    ig,
    v = regulator,
    to = target,
    mode = "out",
    weights = igraph::E(ig)$edge_length
  )[1, 1]

  if (!is.finite(distance)) {
    stop("No path")
  }


  path <- suppressWarnings(
    igraph::shortest_paths(
      ig,
      from = regulator,
      to = target,
      mode = "out",
      output = "both",
      weights = igraph::E(ig)$edge_length
    )
  )

  vpath <- path$vpath[[1]]$name

  if (compare_to_average) {
    avg_path_length <- igraph::mean_distance(ig, directed = TRUE)
  }


  if (compare_to_average) {
    list(
      path = vpath,
      distance = distance,
      distance_over_average = distance / avg_path_length,
      action = ifelse((sum(path$epath[[1]]$corr < 0) %% 2) == 0, 1, -1)
    )
  } else {
    list(
      path = vpath,
      distance = distance,
      action = ifelse((sum(path$epath[[1]]$corr < 0) %% 2) == 0, 1, -1)
    )
  }
}


#' Function to return shortest path from 1 regulator to 1 target in a dynamic network
#'
#' @param network_table a dyanmic network
#' @param regulator the starting regulator
#' @param target the end regulator
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#' @param compare_to_average if TRUE will compute normalized against average path length
#'
#' @return shortest path, distance, normalized distance, and action
#'
#' @export
dynamic_shortest_path <- function(
  network_table,
  regulator,
  target,
  weight_column = "weight",
  compare_to_average = FALSE
) {
  network_table <- do.call("rbind", network_table)
  network_table <- network_table[!duplicated(network_table[, c("regulator", "target")]), ]

  return(
    static_shortest_path(
      network_table,
      regulator,
      target,
      weight_column,
      compare_to_average
    )
  )
}


#' Function to return shortest path from multiple TFs to multiple targets in a dynamic network
#'
#' @param ... Arguments for other methods.
#' @param regulators the starting TFs
#' @param targets the end TFs
#' @param weight_column column name in network_table with edge weights that will be converted to distances
#'
#' @return dataframe with shortest path, distance, normalized distance, and action
#'
#' @export
setGeneric(
  name = "dynamic_shortest_path_multiple",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("dynamic_shortest_path_multiple")
  }
)

#' @param object A dynamic network list (each element a data.frame with regulator/target/weight and optional corr).
#' @rdname dynamic_shortest_path_multiple
#' @export
setMethod(
  "dynamic_shortest_path_multiple",
  signature(object = "list"),
  function(object,
           regulators,
           targets,
           weight_column = "weight") {
    network_table <- object
    res <- data.frame(
      regulator = character(),
      target = character(),
      path = character(),
      distance = numeric(),
      action = numeric()
    )

    for (regulator in regulators) {
      for (target in targets) {
        tryCatch(
          {
            path <- dynamic_shortest_path(
              network_table,
              regulator = regulator,
              target = target,
              weight_column = weight_column,
              compare_to_average = FALSE
            )
            res <- rbind(
              res,
              data.frame(
                regulator = regulator,
                target = target,
                path = paste(path[["path"]], collapse = "--"),
                distance = path[["distance"]],
                action = path[["action"]]
              )
            )
          },
          error = function(e) {}
        )
      }
    }

    network_table <- do.call("rbind", network_table)
    network_table$normalized_score <- thisutils::normalization(
      network_table[, weight_column],
      method = "maximum"
    )
    network_table$edge_length <- 1 - network_table$normalized_score
    ig <- igraph::graph_from_data_frame(
      network_table[, c("regulator", "target", "edge_length", "corr")],
      directed = TRUE
    )

    avg_path_length <- igraph::mean_distance(ig, directed = TRUE)
    res$distance_over_average <- res$distance / avg_path_length

    return(res)
  }
)

#' @param object A Network object (single network).
#' @param regulators Starting genes/TFs.
#' @param targets End genes.
#' @param weight_column Column name for weights.
#' @rdname dynamic_shortest_path_multiple
#' @export
setMethod(
  "dynamic_shortest_path_multiple",
  signature(object = "Network"),
  function(object,
           regulators,
           targets,
           weight_column = "weight") {
    network_table <- tryCatch(export_csn(object), error = function(e) NULL)
    if (is.null(network_table) || nrow(network_table) == 0) {
      stop("No network edges found.")
    }
    network_table <- as.data.frame(network_table)
    if (!"corr" %in% colnames(network_table)) {
      network_table$corr <- network_table[[weight_column]]
    }
    net_list <- list(static = network_table[, c("regulator", "target", weight_column, "corr"), drop = FALSE])
    res <- dynamic_shortest_path_multiple(net_list, regulators = regulators, targets = targets, weight_column = weight_column)
    object@params$shortest_paths <- list(
      method = "dynamic_shortest_path_multiple",
      regulators = regulators,
      targets = targets,
      weight_column = weight_column,
      df = res
    )
    object
  }
)

#' @param object A CSNObject with a dynamic network (states in object@networks[[network]]).
#' @param network Name of the dynamic network.
#' @param celltypes States/cell types to include; NULL for all.
#' @param weight_cutoff Optional cutoff passed to export_csn.
#' @param corr_column Column name for correlation; if absent, fall back to weight.
#' @rdname dynamic_shortest_path_multiple
#' @export
setMethod(
  "dynamic_shortest_path_multiple",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           regulators,
           targets,
           weight_column = "weight",
           weight_cutoff = NULL,
           corr_column = "mean_corr") {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "dynamic",
      verbose = TRUE,
      caller = "dynamic_shortest_path_multiple"
    )
    nets <- GetNetwork(object, network = network)
    if (is.null(nets) || !is.list(nets) || length(nets) == 0) {
      stop("No networks found for the specified network.")
    }
    if (is.null(celltypes)) celltypes <- names(nets)
    celltypes <- intersect(celltypes, names(nets))
    if (length(celltypes) == 0) stop("No matching states/cell types.")

    net_list <- lapply(celltypes, function(sid) {
      net <- nets[[sid]]
      df <- tryCatch(export_csn(net, weight_cutoff = weight_cutoff), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      df <- as.data.frame(df)
      if (!"corr" %in% colnames(df)) {
        df$corr <- if (corr_column %in% colnames(df)) df[[corr_column]] else df[[weight_column]]
      }
      df[, c("regulator", "target", weight_column, "corr"), drop = FALSE]
    })
    net_list <- net_list[!sapply(net_list, is.null)]
    if (length(net_list) == 0) stop("No edges to compute shortest paths.")
    names(net_list) <- celltypes[!sapply(net_list, is.null)]

    res <- dynamic_shortest_path_multiple(net_list, regulators = regulators, targets = targets, weight_column = weight_column)
    state <- .read_multicsn_state(object, init = TRUE)
    shortest_paths <- state$shortest_paths %ss% list()
    shortest_paths[[network]] <- list(
      method = "dynamic_shortest_path_multiple",
      celltypes = names(net_list),
      regulators = regulators,
      targets = targets,
      weight_column = weight_column,
      weight_cutoff = weight_cutoff,
      df = res
    )
    state$shortest_paths <- shortest_paths
    .write_multicsn_state(object, state)
  }
)

#' @rdname dynamic_shortest_path_multiple
#' @export
setMethod(
  "dynamic_shortest_path_multiple",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)

#' Adds an extra column to the result of dynamic_shortest_path_multiple that predicts overall action based on correlation between "from" and "to"
#'
#' @param spDF the result of running dynamic_shortest_path_multiple
#' @param matrix corresponding genes-by-cells expression matrix
#'
#' @return shortest paths dataframe with added action_by_corr column
#'
#' @export
cor_and_add_action <- function(
  spDF,
  matrix
) {
  genes <- union(spDF$from, spDF$to)
  matrix <- matrix[genes, ]

  corrs <- apply(spDF[, c("from", "to")], 1, function(x) {
    stats::cor(matrix[x[1], ], matrix[x[2], ])
  })
  spDF$action_by_corr <- corrs
  spDF$action_by_corr <- ifelse(spDF$action_by_corr < 0, -1, 1)

  spDF
}


#' Useful plotting function to plot heatmap of module expression across time with pre-split matrix
#'
#' @param expList list of expression matrices
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param toScale whether or not to scale expression
#' @param limits limits on expression
#' @param smooth whether or not to smooth expression across pseudotime for cleaner plotting
#' @param order_by name in expList that is used to order rows in the heatmap
#' @param thresh_on threshold expression is considered on, used in ordering the rows
#' @param fontSize heatmap font size
#' @param anno_colors annotation colors
#'
#' @return pheatmap
#'
#' @export
heatmap_by_treatment_group <- function(
  expList,
  meta_data,
  pseudotime_column = "pseudotime",
  toScale = T,
  limits = c(0, 5),
  smooth = TRUE,
  order_by = "WAG",
  thresh_on = 0.02,
  fontSize = 8,
  anno_colors = NA
) {
  if (!is.list(expList)) {
    expList <- list(A = expList)
  }
  expList <- lapply(expList, function(x) {
    st <- meta_data[colnames(x), ]
    st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
    x[, rownames(st)]
  })
  if (smooth) {
    expList <- lapply(expList, function(x) {
      st <- data.frame(
        cell_name = rownames(meta_data[colnames(x), ]),
        pseudotime = meta_data[colnames(x), pseudotime_column]
      )
      rownames(st) <- st$cell_name
      expression_ksmooth(x, st, bandwith = 0.05)
    })
  }

  expList <- lapply(expList, function(x) {
    x[rowSums(is.na(x)) == 0, ]
  })


  peak_time <- apply(expList[[order_by]], 1, function(x) {
    ifelse(any(x > thresh_on), mean(which(x > thresh_on)[1:10]), length(x) + 1)
  })
  peak_time[is.na(peak_time)] <- ncol(expList[[order_by]]) + 1


  genes_ordered <- names(sort(peak_time))

  matrix <- do.call(cbind, expList)
  meta_data <- meta_data[colnames(matrix), ]

  meta_data$cell_name <- rownames(meta_data)
  if ("epoch" %in% colnames(meta_data)) {
    col_ann <- meta_data[, c("treatment", pseudotime_column, "epoch")]
  } else {
    col_ann <- meta_data[, c("treatment", pseudotime_column)]
  }

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


  gaps <- sapply(expList, ncol)
  gaps <- cumsum(gaps)
  gaps <- gaps[-length(gaps)]

  pheatmap::pheatmap(
    value,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = F,
    annotation_colors = anno_colors,
    annotation_col = col_ann,
    border_color = NA,
    gaps_col = gaps,
    annotation_names_row = F,
    fontsize = fontSize
  )
}

#' Useful plotting function to plot heatmap with pre-split matrix
#'
#' @param expList list of expression matrices
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param toScale whether or not to scale expression
#' @param limits limits on expression
#' @param smooth whether or not to smooth expression across pseudotime for cleaner plotting
#' @param fontSize heatmap font size
#' @param anno_colors annotation colors
#'
#' @return pheatmap
#'
#' @export
plot_heatmap_by_treatment <- function(
  expList,
  meta_data,
  pseudotime_column = "latent_time",
  toScale = T,
  limits = c(0, 5),
  smooth = TRUE,
  anno_colors = NA,
  fontSize = 8
) {
  if (!is.list(expList)) {
    expList <- list(A = expList)
  }

  expList <- lapply(expList, function(x) {
    st <- meta_data[colnames(x), ]
    st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
    x[, rownames(st)]
  })
  if (smooth) {
    expList <- lapply(expList, function(x) {
      st <- data.frame(cell_name = rownames(meta_data[colnames(x), ]), pseudotime = meta_data[colnames(x), pseudotime_column])
      rownames(st) <- st$cell_name
      expression_ksmooth(x, st, bandwith = 0.03)
    })
  }


  rows <- data.frame(subnet = rownames(expList[[1]]), network = unlist(lapply(strsplit(rownames(expList[[1]]), "_"), "[[", 1)))
  genes_ordered <- c()
  for (n in unique(rows$network)) {
    peak_time <- apply(expList[[1]][startsWith(rownames(expList[[1]]), n), ], 1, function(x) {
      mean(order(x, decreasing = TRUE)[1:30])
    })
    genes_ordered <- c(genes_ordered, names(sort(peak_time)))
  }

  matrix <- do.call(cbind, expList)
  meta_data <- meta_data[colnames(matrix), ]

  meta_data$cell_name <- rownames(meta_data)
  col_ann <- meta_data[, c("treatment", pseudotime_column)]

  row_ann <- data.frame(subnet = rownames(matrix), network = unlist(lapply(strsplit(rownames(matrix), "_"), "[[", 1)))
  rownames(row_ann) <- row_ann$subnet
  row_ann$subnet <- NULL
  row_ann$network <- factor(row_ann$network, levels = unique(row_ann$network))

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


  gaps <- sapply(expList, ncol)
  gaps <- cumsum(gaps)
  gaps <- gaps[-length(gaps)]

  rowgaps <- table(row_ann$network)
  rowgaps <- cumsum(rowgaps)
  rowgaps <- rowgaps[-length(rowgaps)]

  pheatmap::pheatmap(
    value,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = F,
    annotation_colors = anno_colors,
    annotation_col = col_ann,
    annotation_row = row_ann,
    border_color = NA,
    gaps_col = gaps,
    gaps_row = rowgaps,
    annotation_names_row = F,
    fontsize = fontSize
  )
}

#' Function that orders genes based on peak expression
#'
#' @param matrix expression matrix
#' @param meta_data sample table
#' @param pseudotime_column column in sample table containing pseudotime annotation
#' @param smooth whether or not to smooth expression across pseudotime
#'
#' @return ordered genes
#'
#' @export
order_genes <- function(
  matrix,
  meta_data,
  pseudotime_column,
  smooth = TRUE
) {
  st <- meta_data[colnames(matrix), ]
  st <- st[order(st[, pseudotime_column], decreasing = FALSE), ]
  matrix <- matrix[, rownames(st)]

  if (smooth) {
    st <- data.frame(
      cell_name = rownames(st),
      pseudotime = st[, pseudotime_column]
    )
    rownames(st) <- st$cell_name
    matrix <- expression_ksmooth(matrix, st, bandwith = 0.05)
  }

  peak_time <- apply(matrix, 1, which.max)
  genes_ordered <- names(sort(peak_time))

  return(genes_ordered)
}

find_paths_to <- function(network_table, module) {
  modnet_ig <- igraph::graph_from_data_frame(network_table, directed = TRUE)
  mods <- igraph::V(modnet_ig)$name
  mods <- mods[mods != module]

  res <- igraph::distances(modnet_ig, to = module, weights = igraph::E(modnet_ig)$edge_length, mode = "out")
  colnames(res) <- c("path_length")

  paths <- c()
  for (mod in rownames(res)) {
    path <- igraph::shortest_paths(modnet_ig, from = mod, to = module, mode = "out", weights = igraph::E(modnet_ig)$edge_length)$vpath[[1]]$name
    paths <- c(paths, paste(path, collapse = "--"))
  }

  res <- cbind(res, data.frame(path = paths))

  res$path[is.infinite(res$path_length)] <- NA

  res
}

plot_reachability_community_coverage <- function(
  reachableList,
  communities
) {
  reachableList <- lapply(reachableList, function(x) {
    x[!duplicated(x)]
  })

  res <- data.frame(matrix(ncol = length(reachableList), nrow = length(communities)))
  colnames(res) <- names(reachableList)
  rownames(res) <- names(communities)

  for (m in names(communities)) {
    pct_covered <- sapply(reachableList, function(x) {
      sum(x %in% communities[[m]]) / length(communities[[m]])
    })
    res[m, ] <- pct_covered
  }

  res
}
