#' @include setClass.R
#' @include setGenerics.R
NULL


.edge_uniqueness_list <- function(grnDFs, tfs, weight_column) {
  if (length(grnDFs) < 2) {
    stop("edge_uniqueness requires at least two networks in the list.")
  }
  genes <- unlist(lapply(grnDFs, function(x) union(x$target, x$regulator)), use.names = FALSE)
  genes <- unique(genes)
  tfs <- intersect(genes, tfs)
  if (length(tfs) == 0) {
    stop("No regulators in tfs appear in the networks.")
  }

  adj_list <- sapply(names(grnDFs), function(x) NULL)

  for (df in names(adj_list)) {
    graph <- igraph::graph_from_data_frame(grnDFs[[df]][, c("regulator", "target", weight_column)], directed = TRUE)
    addvtcs <- setdiff(genes, igraph::V(graph)$name)
    graph <- igraph::add_vertices(graph, length(addvtcs), name = addvtcs)
    adj <- igraph::as_adjacency_matrix(graph, attr = weight_column)

    adj <- adj[tfs, ]
    adj <- t(as.matrix(adj))

    adj_list[[df]] <- adj
  }

  full_df <- as.data.frame(as.table(as.matrix(adj_list[[1]])))
  full_df <- full_df[, 1:2]
  colnames(full_df) <- c("target", "regulator")

  for (df in names(adj_list)) {
    add <- as.data.frame(as.table(as.matrix(adj_list[[df]])))
    colnames(add) <- c("target", "regulator", df)

    full_df <- merge(full_df, add, by = c("target", "regulator"))
  }


  full_df$min <- apply(full_df[, !names(full_df) %in% c("target", "regulator")], 1, FUN = min)
  full_df$max <- apply(full_df[, !names(full_df) %in% c("target", "regulator")], 1, FUN = max)

  full_df$diff <- full_df$max - full_df$min


  full_df <- full_df[order(full_df$diff, decreasing = TRUE), ]

  full_df
}

#' @title Edge uniqueness across multiple GRNs
#' @name edge_uniqueness
#' @param object Named list of GRN data.frames (regulator, target, weight) or CSNObject.
#' @param ... Passed to methods.
#' @rdname edge_uniqueness
#' @export
setGeneric(
  name = "edge_uniqueness",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("edge_uniqueness")
  }
)

#' @param tfs Transcription factors to restrict rows of adjacency.
#' @param weight_column Column name for edge weights.
#' @rdname edge_uniqueness
#' @export
setMethod(
  "edge_uniqueness",
  signature(object = "list"),
  function(object, tfs, weight_column) {
    .edge_uniqueness_list(object, tfs, weight_column)
  }
)

#' @param network Name of dynamic network in \code{object@networks}.
#' @param celltypes States to include; \code{NULL} for all.
#' @param weight_cutoff Passed to \code{export_csn}.
#' @param tfs If \code{NULL}, uses \code{object@metadata$tfs}.
#' @rdname edge_uniqueness
#' @export
setMethod(
  "edge_uniqueness",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           tfs = NULL,
           weight_column = "weight",
           weight_cutoff = NULL) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "dynamic",
      verbose = TRUE,
      caller = "edge_uniqueness"
    )
    nets <- GetNetwork(object, network = network)
    if (is.null(nets) || !is.list(nets) || length(nets) < 2) {
      stop("edge_uniqueness requires a dynamic network with at least two states.")
    }
    if (is.null(celltypes)) celltypes <- names(nets)
    celltypes <- intersect(celltypes, names(nets))
    if (length(celltypes) < 2) {
      stop("Need at least two states with networks for edge_uniqueness.")
    }

    grnDFs <- lapply(celltypes, function(sid) {
      df <- tryCatch(export_csn(nets[[sid]], weight_cutoff = weight_cutoff), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      df <- as.data.frame(df)
      if (!weight_column %in% colnames(df)) {
        stop("weight_column '", weight_column, "' not found in exported network.")
      }
      df[, c("regulator", "target", weight_column), drop = FALSE]
    })
    ok <- !sapply(grnDFs, is.null)
    grnDFs <- grnDFs[ok]
    names(grnDFs) <- celltypes[ok]
    if (length(grnDFs) < 2) stop("Fewer than two non-empty networks after export.")

    if (is.null(tfs)) {
      tfs <- .multicsn_get_tfs(object)
    }
    if (is.null(tfs) || length(tfs) == 0) {
      stop("Supply tfs= or initialize TFs with `find_motifs()`.")
    }

    res <- .edge_uniqueness_list(grnDFs, tfs, weight_column)
    state <- .read_multicsn_state(object, init = TRUE)
    edge_uniqueness <- state$edge_uniqueness %ss% list()
    edge_uniqueness[[network]] <- list(
      method = "edge_uniqueness",
      celltypes = names(grnDFs),
      tfs = tfs,
      weight_column = weight_column,
      weight_cutoff = weight_cutoff,
      df = res
    )
    state$edge_uniqueness <- edge_uniqueness
    .write_multicsn_state(object, state)
  }
)

#' @rdname edge_uniqueness
#' @export
setMethod(
  "edge_uniqueness",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)


#' Compute a dynamic difference network
#'
#'
#'
#' @param edgeDF the result of running edge_uniqueness
#' @param epochs list of epoch gene assignments
#' @param condition condition of interest, should be one of the treatment (aka network name) or column names in edgeDF
#' @param type "on" or "off" specifies either finding edges that are active in the condition network but off others or inactive in condition network but active in others
#' @param diff_thresh edge difference threshold to determine if edge is uniquely on or off
#' @param condition_thresh edge weight theshold in condition network to keep or filter out in difference network
#'
#' @return list of dataframes representing the dynamic difference network
#'
#' @export
#'
dynamic_difference_network <- function(
  edgeDF,
  epochs,
  condition,
  type,
  diff_thresh = 3,
  condition_thresh = 6
) {
  edgeDF <- edgeDF[edgeDF$diff != 0, ]


  epochs <- lapply(epochs, function(x) {
    x$mean_expression <- NULL
    x
  })

  temp <- vector(mode = "list", length = length(names(epochs[[1]])))
  names(temp) <- names(epochs[[1]])

  for (e in names(epochs)) {
    temp <- Map(c, temp, epochs[[e]])
  }

  epochs <- sapply(temp, unique)

  dynamic_edges <- list()
  for (epoch in names(epochs)) {
    dynamic_edges[[epoch]] <- edgeDF[edgeDF$regulator %in% epochs[[epoch]], ]
  }

  edgeDF <- dynamic_edges


  conditions <- names(edgeDF[[1]])[!(names(edgeDF[[1]]) %in% c("target", "regulator", "min", "max", "diff"))]
  if (!(condition %in% conditions)) {
    stop("condition not in network.")
  }
  diffnet <- lapply(edgeDF, function(x) x[, c("target", "regulator", condition, "diff")])
  diffnet <- lapply(diffnet, function(x) x[x$diff > diff_thresh, ])

  if (type == "on") {
    diffnet <- Map(function(x, y) x[x$regulator %in% y, ], diffnet, epochs)
    diffnet <- lapply(diffnet, function(x) x[x[, condition] >= condition_thresh, ])
  } else if (type == "off") {
    diffnet <- lapply(diffnet, function(x) x[x[, condition] < condition_thresh, ])
  } else {
    stop("type should be 'on' or 'off'.")
  }

  diffnet
}


#' Perform community detection on a dynamic network
#'
#' @param diffnet diffnet
#' @param method community detection method, currently only "louvain"
#' @param use_weights whether or not to use edge weights (for weighted graphs)
#' @param weight_column if using weights, name of the column containing edge weights
#'
#' @return community assignments of nodes in the dynamic network
#'
#' @export
diffnet_community_detection <- function(
  diffnet,
  method = "louvain",
  use_weights = FALSE,
  weight_column = NULL
) {
  graphs <- lapply(diffnet, function(x) {
    g <- igraph::graph_from_data_frame(x, directed = FALSE)
    g
  })

  if (use_weights) {
    weights <- igraph::edge_attr(x, weight_column)
  } else {
    weights <- NA
  }

  if (method == "louvain") {
    communities <- lapply(graphs, function(x) {
      c <- igraph::cluster_louvain(x, weights = weights)
      c
    })
  }

  communities
}


#' Adds interaction type to dynamic differential network
#'
#' @param diffnet diffnet
#' @param type "on" or "off" depending on the type of differential network. If "on" will assign type based on grnDF_on. Otherwise interaction assigned from grnDF_offlist.
#' @param grnDF_on the static network in which diffnet edges are active
#' @param grnDF_offlist list of static networks in which diffnet edges are inactive
#'
#' @return community assignments of nodes in the dynamic network
#'
#' @export
add_type <- function(
  diffnet,
  type,
  grnDF_on,
  grnDF_offlist
) {
  diffnet <- diffnet[sapply(diffnet, function(x) dim(x)[1]) > 0]
  fun <- function(x, y) {
    merge(x, y, by = c("target", "regulator"), all.x = TRUE)
  }
  if (type == "on") {
    added <- lapply(diffnet, fun, grnDF_on[, c("target", "regulator", "corr")])
    added <- lapply(added, transform, interaction = ifelse(corr < 0, "repression", "activation"))
  } else if (type == "off") {
    added <- diffnet
    for (i in 1:length(grnDF_offlist)) {
      added <- lapply(added, fun, grnDF_offlist[[i]][, c("target", "regulator", "corr")])
    }
    added <- lapply(added, function(x) {
      cols <- colnames(x)[grepl("corr", colnames(x))]
      x$interaction <- NA
      x[(rowSums((x[, cols] >= 0) | (is.na(x[, cols])), na.rm = TRUE)) == length(cols), "interaction"] <- "activation"
      x[(rowSums((x[, cols] < 0) | (is.na(x[, cols])), na.rm = TRUE)) == length(cols), "interaction"] <- "repression"
      x
    })
  } else {
    stop("invalid type.")
  }

  added
}


#' Computes frobenius distance in a pairwise manner between two sets of networks
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param weight_column	column name containing edge weights
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of frobenius distances
#'
#' @export
compute_frobenius_distance <- function(
  netlist1,
  netlist2,
  weight_column = "weight",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
) {
  if (compare_within_netlist1) {
    df1 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()
    for (i in 1:nrow(df1)) {
      net1 <- netlist1[[df1$X1[i]]]
      net2 <- netlist1[[df1$X2[i]]]

      net1 <- reshape2::acast(
        net1,
        regulator ~ target,
        value.var = weight_column
      )
      net1[is.na(net1)] <- 0
      net1 <- t(net1)

      net2 <- reshape2::acast(
        net2,
        regulator ~ target,
        value.var = weight_column
      )
      net2[is.na(net2)] <- 0
      net2 <- t(net2)

      rows <- union(rownames(net1), rownames(net2))
      cols <- union(colnames(net1), colnames(net2))


      missing <- setdiff(rows, rownames(net1))
      addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
      rownames(addnet1) <- missing
      colnames(addnet1) <- colnames(net1)
      net1 <- rbind(net1, addnet1)

      missing <- setdiff(rows, rownames(net2))
      addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
      rownames(addnet2) <- missing
      colnames(addnet2) <- colnames(net2)
      net2 <- rbind(net2, addnet2)

      missing <- setdiff(cols, colnames(net1))
      addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
      rownames(addnet1) <- rownames(net1)
      colnames(addnet1) <- missing
      net1 <- cbind(net1, addnet1)

      missing <- setdiff(cols, colnames(net2))
      addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
      rownames(addnet2) <- rownames(net2)
      colnames(addnet2) <- missing
      net2 <- cbind(net2, addnet2)


      net1 <- net1[rows, cols]
      net2 <- net2[rows, cols]

      diff <- net1 - net2
      score <- c(score, norm(diff, "F"))
    }

    df1$score <- score
  }


  if (compare_within_netlist2) {
    df2 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()

    for (i in 1:nrow(df2)) {
      net1 <- netlist2[[df2$X1[i]]]
      net2 <- netlist2[[df2$X2[i]]]

      net1 <- reshape2::acast(
        net1,
        regulator ~ target,
        value.var = weight_column
      )
      net1[is.na(net1)] <- 0
      net1 <- t(net1)

      net2 <- reshape2::acast(
        net2,
        regulator ~ target,
        value.var = weight_column
      )
      net2[is.na(net2)] <- 0
      net2 <- t(net2)

      rows <- union(rownames(net1), rownames(net2))
      cols <- union(colnames(net1), colnames(net2))


      missing <- setdiff(rows, rownames(net1))
      addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
      rownames(addnet1) <- missing
      colnames(addnet1) <- colnames(net1)
      net1 <- rbind(net1, addnet1)

      missing <- setdiff(rows, rownames(net2))
      addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
      rownames(addnet2) <- missing
      colnames(addnet2) <- colnames(net2)
      net2 <- rbind(net2, addnet2)

      missing <- setdiff(cols, colnames(net1))
      addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
      rownames(addnet1) <- rownames(net1)
      colnames(addnet1) <- missing
      net1 <- cbind(net1, addnet1)

      missing <- setdiff(cols, colnames(net2))
      addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
      rownames(addnet2) <- rownames(net2)
      colnames(addnet2) <- missing
      net2 <- cbind(net2, addnet2)


      net1 <- net1[rows, cols]
      net2 <- net2[rows, cols]

      diff <- net1 - net2
      score <- c(score, norm(diff, "F"))
    }

    df2$score <- score
  }


  df3 <- data.frame(
    expand.grid(1:length(netlist1), 1:length(netlist2))
  )
  colnames(df3) <- c("X1", "X2")
  score <- c()

  for (i in 1:nrow(df3)) {
    net1 <- netlist1[[df3$X1[i]]]
    net2 <- netlist2[[df3$X2[i]]]

    net1 <- reshape2::acast(
      net1,
      regulator ~ target,
      value.var = weight_column
    )
    net1[is.na(net1)] <- 0
    net1 <- t(net1)

    net2 <- reshape2::acast(
      net2,
      regulator ~ target,
      value.var = weight_column
    )
    net2[is.na(net2)] <- 0
    net2 <- t(net2)

    rows <- union(rownames(net1), rownames(net2))
    cols <- union(colnames(net1), colnames(net2))


    missing <- setdiff(rows, rownames(net1))
    addnet1 <- matrix(0, nrow = length(missing), ncol = ncol(net1))
    rownames(addnet1) <- missing
    colnames(addnet1) <- colnames(net1)
    net1 <- rbind(net1, addnet1)

    missing <- setdiff(rows, rownames(net2))
    addnet2 <- matrix(0, nrow = length(missing), ncol = ncol(net2))
    rownames(addnet2) <- missing
    colnames(addnet2) <- colnames(net2)
    net2 <- rbind(net2, addnet2)

    missing <- setdiff(cols, colnames(net1))
    addnet1 <- matrix(0, nrow = nrow(net1), ncol = length(missing))
    rownames(addnet1) <- rownames(net1)
    colnames(addnet1) <- missing
    net1 <- cbind(net1, addnet1)

    missing <- setdiff(cols, colnames(net2))
    addnet2 <- matrix(0, nrow = nrow(net2), ncol = length(missing))
    rownames(addnet2) <- rownames(net2)
    colnames(addnet2) <- missing
    net2 <- cbind(net2, addnet2)


    net1 <- net1[rows, cols]
    net2 <- net2[rows, cols]

    diff <- net1 - net2
    score <- c(score, norm(diff, "F"))
  }

  df3$score <- score


  if (compare_within_netlist1) {
    df1$group <- "netlist1"
  }
  if (compare_within_netlist2) {
    df2$group <- "netlist2"
  }
  df3$group <- "cross_comparison"

  if (compare_within_netlist1 & compare_within_netlist2) {
    df <- rbind(df1, df2)
    df <- rbind(df, df3)
  } else if (compare_within_netlist2) {
    df <- rbind(df2, df3)
  } else if (compare_within_netlist1) {
    df <- rbind(df1, df3)
  } else {
    df <- df3
  }

  df
}


#' Computes Jaccard similarity between top regulators in two sets of networks
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param n_regs the number of top regulators to compare from each network
#' @param method method to find top regulators. Currently only supports "pagerank"
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of Jaccard similarities of top regulators
#'
#' @export
#'
compute_JI_topregs <- function(
  netlist1,
  netlist2,
  n_regs = 15,
  method = "pagerank",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
) {
  if (compare_within_netlist1) {
    df1 <- data.frame(t(utils::combn(1:length(netlist1), 2)))
    score <- c()

    for (i in 1:nrow(df1)) {
      net1 <- netlist1[[df1$X1[i]]]
      net2 <- netlist1[[df1$X2[i]]]

      net1_ranks <- compute_pagerank(list(X = net1))
      net2_ranks <- compute_pagerank(list(X = net2))

      net1_topregs <- net1_ranks$X$gene[1:n_regs]
      net2_topregs <- net2_ranks$X$gene[1:n_regs]

      ji <- length(
        intersect(net1_topregs, net2_topregs)
      ) / length(union(net1_topregs, net2_topregs))

      score <- c(score, ji)
    }

    df1$jaccard <- score
  }


  if (compare_within_netlist2) {
    df2 <- data.frame(t(utils::combn(1:length(netlist2), 2)))
    score <- c()

    for (i in 1:nrow(df2)) {
      net1 <- netlist2[[df2$X1[i]]]
      net2 <- netlist2[[df2$X2[i]]]

      net1_ranks <- compute_pagerank(list(X = net1))
      net2_ranks <- compute_pagerank(list(X = net2))

      net1_topregs <- net1_ranks$X$gene[1:n_regs]
      net2_topregs <- net2_ranks$X$gene[1:n_regs]

      ji <- length(
        intersect(net1_topregs, net2_topregs)
      ) / length(union(net1_topregs, net2_topregs))

      score <- c(score, ji)
    }

    df2$jaccard <- score
  }


  df3 <- data.frame(expand.grid(1:length(netlist1), 1:length(netlist2)))
  colnames(df3) <- c("X1", "X2")
  score <- c()

  for (i in 1:nrow(df3)) {
    net1 <- netlist1[[df3$X1[i]]]
    net2 <- netlist2[[df3$X2[i]]]

    net1_ranks <- compute_pagerank(list(X = net1))
    net2_ranks <- compute_pagerank(list(X = net2))

    net1_topregs <- net1_ranks$X$gene[1:n_regs]
    net2_topregs <- net2_ranks$X$gene[1:n_regs]

    ji <- length(
      intersect(net1_topregs, net2_topregs)
    ) / length(union(net1_topregs, net2_topregs))

    score <- c(score, ji)
  }

  df3$jaccard <- score


  if (compare_within_netlist1) {
    df1$group <- "netlist1"
  }
  if (compare_within_netlist2) {
    df2$group <- "netlist2"
  }
  df3$group <- "cross_comparison"

  if (compare_within_netlist1 & compare_within_netlist2) {
    df <- rbind(df1, df2)
    df <- rbind(df, df3)
  } else if (compare_within_netlist2) {
    df <- rbind(df2, df3)
  } else if (compare_within_netlist1) {
    df <- rbind(df1, df3)
  } else {
    df <- df3
  }

  df
}


#' Computes Jaccard similarity between top regulators in two sets of networks across a range of top X regulators
#'
#' @param netlist1 list of grnDFs
#' @param netlist2 list of grnDFs
#' @param n_regs a vector indicating which values of top regulators to scan across
#' @param func func
#' @param method method to find top regulators. Currently only supports "pagerank"
#' @param weight_column column name in grnDFs containing edge weights
#' @param compare_within_netlist1 whether or not to do pairwise comparisons between networks in netlist1
#' @param compare_within_netlist2 whether or not to do pairwise comparisons between networks in netlist2
#'
#' @return dataframe of Jaccard similarities of top regulators
#'
#' @export
JI_across_topregs <- function(
  netlist1,
  netlist2,
  n_regs = 3:15,
  func = "mean",
  method = "pagerank",
  weight_column = "zscore",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
) {
  if (func == "mean") {
    res <- data.frame(
      group = character(),
      jaccard = numeric(),
      n_regs = numeric()
    )
    for (i in n_regs) {
      ji <- compute_JI_topregs(
        netlist1,
        netlist2,
        n_regs = i,
        method = method,
        compare_within_netlist1 = compare_within_netlist1,
        compare_within_netlist2 = compare_within_netlist2
      )
      ji <- ji[, c("group", "jaccard")]
      ji <- aggregate(. ~ group, ji, mean)
      ji$n_regs <- i

      res <- rbind(res, ji)
    }
  }
  res
}


compute_betweenness <- function(
  grnDF,
  tfs = NULL,
  normalized = TRUE
) {
  if (is.null(tfs)) {
    tfs <- unique(grnDF$regulator)
  }

  grnDF <- grnDF[, c("regulator", "target", "weight")]
  g <- igraph::graph_from_data_frame(grnDF, directed = FALSE)

  betweenness <- igraph::betweenness(
    g,
    tfs,
    directed = FALSE,
    normalized = normalized
  )

  betweenness
}


#' Computes betweenness and degree of each regulator for each network in a list of networks
#'
#' @param netlist netlist
#' @param weight_column column name in grnDFs containing edge weights
#' @param normalized whether or not to normalize degree and betweenness
#'
#' @return dataframe listing network, betweenness of each regulator, degree of each regulator
#'
#' @export
biglist_compute_betweenness_degree <- function(
  netlist,
  weight_column = "zscore",
  normalized = TRUE
) {
  res <- data.frame(
    betweenness = numeric(),
    degree = numeric(),
    network = character(),
    regulator = character()
  )

  if (is.null(names(netlist))) {
    names(netlist) <- seq(1:length(netlist))
  }

  for (name in names(netlist)) {
    net <- netlist[[name]]
    g <- igraph::graph_from_data_frame(
      net[, c("regulator", "target")],
      directed = FALSE
    )
    betweenness <- compute_betweenness(
      net,
      normalized = normalized
    )
    degree <- degree(
      g,
      unique(net$regulator),
      mode = "all",
      normalized = normalized
    )

    df <- cbind(betweenness, degree = degree[names(betweenness)])
    df <- as.data.frame(df)
    df$network <- name
    df$regulator <- rownames(df)

    res <- rbind(res, df)
  }

  res
}


#' Plot the dynamic differential network
#'
#' @param grn the dynamic network
#' @param tfs TFs
#' @param only_TFs whether or not to only plot TFs and exclude non-regulators
#' @param order the order in which to plot epochs, or which epochs to plot
#'
#' @return plot
#'
#' @export
#'
plot_dyn_diffnet <- function(
  grn,
  tfs,
  only_TFs = TRUE,
  order = NULL
) {
  g <- list()

  if (!is.null(order)) {
    grn <- grn[order]
  }

  for (i in 1:length(grn)) {
    df <- grn[[i]]

    if (only_TFs) {
      df <- df[df$target %in% tfs, ]
    }
    if (nrow(df) == 0) {
      next
    }

    net <- igraph::graph_from_data_frame(
      df[, c("regulator", "target", "interaction")],
      directed = FALSE
    )
    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    tfnet <- ggnetwork::ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    tfnet$is_regulator <- as.character(tfnet$name %in% tfs)

    cols <- c("activation" = "blue", "repression" = "red")
    g[[i]] <- ggplot() +
      geom_edges(
        data = tfnet,
        aes(x = x, y = y, xend = xend, yend = yend, color = interaction),
        size = 0.75,
        curvature = 0.1,
        alpha = .6
      ) +
      scale_color_manual(values = cols) +
      geom_nodes(
        data = tfnet,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "darkgray",
        size = 6,
        alpha = .5
      ) +
      geom_nodes(
        data = tfnet[tfnet$is_regulator == "TRUE", ],
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "#8C4985",
        size = 6,
        alpha = .8
      ) +
      geom_nodelabel_repel(
        data = tfnet,
        aes(x = x, y = y, label = name),
        size = 2.5,
        color = "#5A8BAD"
      ) +
      theme_blank() +
      ggtitle(names(grn)[i])

    g[[i]] <- g[[i]] + theme(legend.position = "none")
  }
  g <- g[!sapply(g, is.null)]
  do.call(gridExtra::grid.arrange, g)
}


#' Plot the dynamic differential network but colored by communities and optionally faded by betweenness
#'
#' @param grn the dynamic network
#' @param tfs TFs
#' @param only_TFs whether or not to only plot TFs and exclude non-regulators
#' @param order the order in which to plot epochs, or which epochs to plot
#' @param compute_betweenness whether or not to fade nodes by betweenness
#'
#' @return plot
#'
#' @export
plot_diffnet_detail <- function(
  grn,
  tfs,
  only_TFs = TRUE,
  order = NULL,
  compute_betweenness = TRUE
) {
  g <- list()

  if (!is.null(order)) {
    grn <- grn[order]
  }

  if (any(sapply(grn, function(x) {
    ("interaction" %in% names(x))
  }) == FALSE)) {
    grn <- add_interaction_type(grn)
  }

  for (i in 1:length(grn)) {
    df <- grn[[i]]

    if (only_TFs) {
      df <- df[df$target %in% tfs, ]
    }
    if (nrow(df) == 0) {
      next
    }

    df <- df[, c("regulator", "target", "weight", "interaction")]

    net <- igraph::graph_from_data_frame(df, directed = FALSE)

    if (compute_betweenness) {
      b <- igraph::betweenness(net, directed = FALSE, normalized = TRUE)
      b <- as.data.frame(b)
      colnames(b) <- "betweenness"
      b$gene <- rownames(b)
      b <- b[, c("gene", "betweenness")]
      c <- as.data.frame(as.table(igraph::membership(igraph::cluster_louvain(net))))
      colnames(c) <- c("gene", "communities")
      vtx_features <- merge(b, c, by = "gene", all = TRUE)
    } else {
      c <- as.data.frame(as.table(igraph::membership(igraph::cluster_louvain(net))))
      colnames(c) <- c("gene", "communities")
      vtx_features <- c
    }

    layout <- igraph::layout_with_fr(net)
    rownames(layout) <- igraph::V(net)$name
    layout_ordered <- layout[igraph::V(net)$name, ]
    tfnet <- ggnetwork(net, layout = layout_ordered, cell.jitter = 0)
    tfnet$is_regulator <- as.character(tfnet$name %in% tfs)
    if (compute_betweenness) {
      tfnet$betweenness <- vtx_features$betweenness[match(tfnet$name, vtx_features$gene)]
    }
    tfnet$communities <- as.factor(vtx_features$communities[match(tfnet$name, vtx_features$gene)])

    cols <- c("activation" = "blue", "repression" = "red")
    num_cols2 <- length(unique(tfnet$communities))
    if (num_cols2 <= 8) {
      cols2 <- RColorBrewer::brewer.pal(num_cols2, "Set1")
    } else {
      cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_cols2)
    }
    names(cols2) <- unique(tfnet$communities)
    cols <- c(cols, cols2)

    if (compute_betweenness) {
      g[[i]] <- ggplot() +
        geom_edges(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend, color = interaction), size = 0.75, curvature = 0.1, alpha = .6) +
        scale_color_manual(values = cols) +
        geom_nodes(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgray", size = 6, alpha = .5) +
        geom_nodes(data = tfnet[tfnet$is_regulator == "TRUE", ], aes(x = x, y = y, xend = xend, yend = yend, color = communities, alpha = betweenness + .5), size = 6) +
        geom_nodelabel_repel(data = tfnet, aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD") +
        theme_blank() +
        ggtitle(names(grn)[i])
    } else {
      g[[i]] <- ggplot() +
        geom_edges(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend, color = interaction), size = 0.75, curvature = 0.1, alpha = .6) +
        scale_color_manual(values = cols) +
        geom_nodes(data = tfnet, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgray", size = 6, alpha = .5) +
        geom_nodes(data = tfnet[tfnet$is_regulator == "TRUE", ], aes(x = x, y = y, xend = xend, yend = yend, color = communities), alpha = .5, size = 6) +
        geom_nodelabel_repel(data = tfnet, aes(x = x, y = y, label = name), size = 2.5, color = "#5A8BAD") +
        theme_blank() +
        ggtitle(names(grn)[i])
    }

    g[[i]] <- g[[i]] + theme(legend.position = "none")
  }
  g <- g[!sapply(g, is.null)]
  do.call(gridExtra::grid.arrange, g)
}
