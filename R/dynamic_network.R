#' Define states
#'
#' @param dynamic_object A trajectory result with a named vector of gene P values in \code{genes} and a cell table in \code{cells}, or a named list of such results. The cell table contains \code{pseudotime} and cell identifiers; group-based partitioning also requires \code{group}.
#' @param matrix genes-by-cells expression matrix, or a list of expression matrices per path. If list, names should match names of dynamic_object.
#' @param method method to define states. Either "pseudotime", "cell_order", "group", "con_similarity", "kmeans", "hierarchical"
#' @param num_states number of states to define. Ignored when explicit pseudotime cuts or group assignments define the states.
#' @param pseudotime_cuts vector of pseudotime cutoffs. If NULL, cuts are set to max(pseudotime)/num_states.
#' @param group_assignments a list of vectors where names(assignment) are state names, and vectors contain groups belonging to corresponding state
#' @param p_value p_value
#' @param winSize winSize
#'
#' @return updated list of dynamic_object with state column included in dynamic_object$cells
#' @export
define_states <- function(
  dynamic_object,
  matrix,
  method = "pseudotime",
  num_states = 2,
  pseudotime_cuts = NULL,
  group_assignments = NULL,
  p_value = 0.05,
  winSize = 2
) {
  dynamic_object <- normalize_state_data(dynamic_object)
  if (!is.list(dynamic_object[[1]])) {
    dynamic_object <- list(dynamic_object)
    matrix <- list(matrix)
  }

  new_dynRes <- list()
  for (path in 1:length(dynamic_object)) {
    if (!is.null(names(dynamic_object))) {
      path <- names(dynamic_object)[path]
    }
    path_dyn <- dynamic_object[[path]]
    path_dyn$cells$state <- NA


    if (method == "pseudotime") {
      if (is.null(pseudotime_cuts)) {
        pseudotime_cuts <- seq(
          min(path_dyn$cells$pseudotime),
          max(path_dyn$cells$pseudotime),
          (max(path_dyn$cells$pseudotime) - min(path_dyn$cells$pseudotime)) / num_states
        )
        pseudotime_cuts <- pseudotime_cuts[-length(pseudotime_cuts)]
        pseudotime_cuts <- pseudotime_cuts[-1]
      }

      path_dyn <- split_states_by_pseudotime(path_dyn, pseudotime_cuts)
    }


    if (method == "cell_order") {
      t1 <- path_dyn$cells$pseudotime
      names(t1) <- as.vector(path_dyn$cells$cell_name)

      chunk_size <- floor(length(t1) / num_states)

      for (i in 1:num_states) {
        if (i == num_states) {
          cells_in_state <- names(t1)[(1 + ((i - 1) * chunk_size)):length(t1)]
        } else {
          cells_in_state <- names(t1)[(1 + ((i - 1) * chunk_size)):(i * chunk_size)]
        }
        path_dyn$cells$state[path_dyn$cells$cell_name %in% cells_in_state] <- paste0("state_", i)
      }
    }


    if (method == "group") {
      if (is.null(group_assignments)) {
        stop("Must provide group_assignments for group method.")
      }

      path_dyn <- split_states_by_group(path_dyn, group_assignments)
    }


    if (method == "con_similarity") {
      cuts <- find_cuts_by_similarity(matrix[[path]], path_dyn, winSize = winSize, p_value = p_value)
      path_dyn <- split_states_by_pseudotime(path_dyn, cuts)
    }


    if (method == "kmeans") {
      cuts <- find_cuts_by_clustering(
        matrix[[path]],
        path_dyn,
        num_states = num_states,
        method = "kmeans",
        p_value = p_value
      )
      path_dyn <- split_states_by_pseudotime(path_dyn, cuts)
    }


    if (method == "hierarchical") {
      cuts <- find_cuts_by_clustering(matrix[[path]], path_dyn, num_states = num_states, method = "hierarchical", p_value = p_value)
      path_dyn <- split_states_by_pseudotime(path_dyn, cuts)
    }


    new_dynRes[[path]] <- path_dyn
  }

  if (length(new_dynRes) == 1) {
    new_dynRes <- new_dynRes[[1]]
  }
  new_dynRes
}

#' Splits data into states
#'
#' Splits data into states by assigning cells to states
#'
#' @param dynamic_object result of running findDynGenes or compileDynGenes
#' @param cuts vector of pseudotime cutoffs
#' @param state_names names of resulting states, must have length of length(cuts)+1
#'
#' @return updated dynamic_object with state column included in dynamic_object$cells
#' @export
split_states_by_pseudotime <- function(dynamic_object, cuts, state_names = NULL) {
  dynamic_object <- normalize_state_data(dynamic_object)
  state_names <- normalize_state_labels(state_names)
  sampTab <- dynamic_object$cells

  if (max(cuts) > max(sampTab$pseudotime)) {
    stop("Cuts must be within pseudotime.")
  }

  if (!is.null(state_names) & (length(state_names) != length(cuts) + 1)) {
    stop("Length of state_names must be equal to 1+length(cuts).")
  }

  if (is.null(state_names)) {
    state_names <- paste0("state_", seq_len(length(cuts) + 1))
  }

  cuts <- c(-0.1, cuts, max(sampTab$pseudotime))
  sampTab$state <- NA
  for (i in 2:length(cuts)) {
    sampTab$state[(cuts[i - 1] < sampTab$pseudotime) & (sampTab$pseudotime <= cuts[i])] <- state_names[i - 1]
  }

  dynamic_object$cells <- sampTab
  dynamic_object
}

#' Splits data into states manually
#'
#' Splits data into states given group assignment
#'
#' @param dynamic_object result of running findDynGenes or compileDynGenes
#' @param assignment a list of vectors where names(assignment) are state names, and vectors contain groups belonging to corresponding state
#'
#' @return updated dynamic_object with state column included in dynamic_object$cells
#' @export
split_states_by_group <- function(dynamic_object, assignment) {
  dynamic_object <- normalize_state_data(dynamic_object)
  assignment <- normalize_state_data(assignment)
  sampTab <- dynamic_object$cells
  sampTab$state <- NA

  for (e in names(assignment)) {
    sampTab$state[sampTab$group %in% assignment[[e]]] <- e
  }

  dynamic_object$cells <- sampTab
  dynamic_object
}

#' @title find_cuts_by_similarity
#'
#' @description
#'  Returns cuts to define states via sliding window comparison
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object result of running findDynGenes or define_states
#' @param winSize number of cells to each side to compare each cell to
#' @param limit_to vector of genes on which to base state cuts, for example, limiting to TFs
#' @param p_value pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into states
#' @export
find_cuts_by_similarity <- function(
  matrix,
  dynamic_object,
  winSize = 2,
  limit_to = NULL,
  p_value = 0.05
) {
  matrix <- matrix[names(dynamic_object$genes[dynamic_object$genes < p_value]), ]

  if (!is.null(limit_to)) {
    limit_to <- intersect(limit_to, rownames(matrix))
    matrix <- matrix[limit_to, ]
  }


  xcorr <- stats::cor(matrix[, rownames(dynamic_object$cells)])
  pShift <- rep(0, nrow(xcorr) - 1)
  end <- nrow(xcorr)

  for (i in 2:end - 1) {
    past <- seq(i - 1:max(1, i - winSize))
    pastandfuture <- seq(i + 1:min(end, i + winSize))
    closer_to <- which.max(c(mean(xcorr[i, past]), mean(xcorr[i, pastandfuture])))
    pShift[i] <- c(0, 1)[closer_to]
  }

  consecutive_diffs <- diff(pShift)
  cuts_index <- which(consecutive_diffs == 1)


  cat("Cut points: ", dynamic_object$cells[cuts_index, ]$pseudotime, "\n")


  dynamic_object$cells[cuts_index, ]$pseudotime
}


#' Returns cuts to define states
#'
#' Returns cuts to define states via clustering
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object result of running findDynGenes or define_states
#' @param num_states the number of states
#' @param limit_to vector of genes on which to base state cuts, for example, limiting to TFs
#' @param method what clustering method to use, either 'kmeans' or 'hierarchical'
#' @param p_value pval threshold if gene is dynamically expressed
#'
#' @return vector of  pseudotimes at which to cut data into states
#' @export
find_cuts_by_clustering <- function(
  matrix,
  dynamic_object,
  num_states,
  limit_to = NULL,
  method = "kmeans",
  p_value = 0.05
) {
  matrix <- matrix[names(dynamic_object$genes[dynamic_object$genes < p_value]), ]


  if (!is.null(limit_to)) {
    limit_to <- intersect(limit_to, rownames(matrix))
    matrix <- matrix[limit_to, ]
  }

  if (method == "kmeans") {
    clustering <- stats::kmeans(t(matrix), num_states, iter.max = 100)$cluster

    cuts <- c()
    for (cluster in unique(clustering)) {
      cells <- names(clustering)[clustering == cluster]
      cuts <- c(cuts, max(dynamic_object$cells$pseudotime[rownames(dynamic_object$cells) %in% cells]))
    }
    cuts <- sort(cuts)
    cuts <- cuts[1:length(cuts) - 1]
  } else if (method == "hierarchical") {
    clustering <- stats::hclust(stats::dist(t(matrix)), method = "centroid")
    clusterCut <- stats::cutree(clustering, 3)
    cuts <- c()
    for (cluster in unique(clusterCut)) {
      cells <- names(clusterCut)[clusterCut == cluster]
      cuts <- c(cuts, max(dynamic_object$cells$pseudotime[rownames(dynamic_object$cells) %in% cells]))
    }
    cuts <- sort(cuts)
    cuts <- cuts[1:length(cuts) - 1]
  } else {
    stop("method must be either 'kmeans' or 'hierarchical'.")
  }

  cat("Cut points: ", cuts, "\n")
  cuts
}

#' Assigns genes to states
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_states
#' @param method method of assigning state genes, either "active_expression" (looks for active expression in state) or "DE" (looks for differentially expressed genes per state)
#' @param p_value pval threshold if gene is dynamically expressed
#' @param pThresh_DE pval if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into state with max expression.
#'
#' @return states a list detailing genes active in each state
#' @export
#'
assign_genes_to_states <- function(
  matrix,
  dynamic_object,
  method = "active_expression",
  p_value = 0.05,
  pThresh_DE = 0.05,
  active_thresh = 0.33,
  toScale = FALSE,
  forceGenes = TRUE
) {
  dynamic_object <- normalize_state_data(dynamic_object)
  if (active_thresh < 0 | active_thresh > 1) {
    stop("active_thresh must be between 0 and 1.")
  }


  exp <- Matrix::as.matrix(
    matrix[intersect(rownames(matrix), names(dynamic_object$genes[dynamic_object$genes < p_value])), ]
  )

  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }


  state_names <- unique(dynamic_object$cells$state)
  states <- vector("list", length(state_names))
  names(states) <- state_names

  navg <- ceiling(ncol(exp) * 0.05)


  thresholds <- data.frame(
    gene = rownames(exp),
    thresh = rep(0, length(rownames(exp)))
  )
  rownames(thresholds) <- thresholds$gene
  for (gene in rownames(exp)) {
    profile <- exp[gene, ][order(exp[gene, ], decreasing = FALSE)]
    bottom <- mean(profile[1:navg])
    top <- mean(profile[(length(profile) - navg):length(profile)])

    thresh <- ((top - bottom) * active_thresh) + bottom
    thresholds[gene, "thresh"] <- thresh
  }

  mean_expression <- data.frame(
    gene = character(),
    state = numeric(),
    mean_expression = numeric()
  )
  if (method == "active_expression") {
    for (state in names(states)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$state == state, "cells"]
      chunk <- exp[, chunk_cells]

      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      states[[state]] <- rownames(chunk_df[chunk_df$active, ])

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          state = rep(state, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  } else {
    for (state in names(states)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$state == state, "cells"]
      chunk <- exp[, chunk_cells]

      background <- exp[, !(colnames(exp) %in% chunk_cells)]

      diffres <- data.frame(gene = character(), mean_diff = double(), pval = double())
      for (gene in rownames(exp)) {
        t <- stats::t.test(chunk[gene, ], background[gene, ])
        ans <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, ans)
      }

      diffres$padj <- stats::p.adjust(diffres$pval, method = "BH")
      diffres <- diffres[diffres$mean_diff > 0, ]

      states[[state]] <- diffres$gene[diffres$padj < pThresh_DE]


      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      states[[state]] <- intersect(states[[state]], rownames(chunk_df[chunk_df$active, ]))

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          state = rep(state, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  }


  if (forceGenes) {
    assignedGenes <- unique(unlist(states))
    orphanGenes <- setdiff(rownames(exp), assignedGenes)
    message("There are ", length(orphanGenes), " orphan genes\n")
    for (oGene in orphanGenes) {
      xdat <- mean_expression[mean_expression$gene == oGene, ]
      state_id <- xdat[which.max(xdat$mean_expression), ]$state
      states[[state_id]] <- append(states[[state_id]], oGene)
    }
  }

  states$mean_expression <- mean_expression

  states
}

#' Divides grnDF into states, filters interactions between genes not in same or consecutive states
#'
#' @param grnDF result of GRN reconstruction
#' @param states result of running assign_genes_to_states
#' @param state_network dataframe outlining higher level state connectivity (i.e. state transition network).
#' If NULL, will assume states is ordered linear trajectory
#'
#' @return list of GRNs across states and transitions
#' @export
split_network_by_states <- function(
  grnDF,
  states,
  state_network = NULL
) {
  states <- normalize_state_data(states)
  if (!is.null(state_network)) {
    state_network$from <- normalize_state_labels(state_network$from)
    state_network$to <- normalize_state_labels(state_network$to)
  }
  states$mean_expression <- NULL
  all_dyngenes <- unique(unlist(states, use.names = FALSE))


  if (is.null(state_network)) {
    state_network <- data.frame(from = character(), to = character())
    for (i in 1:(length(names(states)) - 1)) {
      df <- data.frame(from = c(names(states)[i]), to = c(names(states)[i + 1]))
      state_network <- rbind(state_network, df)
    }
  }


  state_network <- rbind(
    state_network,
    data.frame(from = names(states), to = names(states))
  )


  state_network$name <- paste(state_network[, 1], state_network[, 2], sep = "..")
  GRN <- vector("list", nrow(state_network))
  names(GRN) <- state_network$name

  print(state_network)

  for (t in 1:nrow(state_network)) {
    from <- as.character(state_network[t, 1])
    to <- as.character(state_network[t, 2])


    temp <- grnDF[grnDF$regulator %in% states[[from]], ]


    if (from != to) {
      remove_tgs_in_both <- intersect(states[[to]], states[[from]])
      remove_tgs_in_neither <- intersect(setdiff(all_dyngenes, states[[from]]), setdiff(all_dyngenes, states[[to]]))

      temp <- temp[!(temp$TG %in% remove_tgs_in_both), ]
      temp <- temp[!(temp$TG %in% remove_tgs_in_neither), ]
    }


    GRN[[state_network[t, "name"]]] <- temp
  }

  GRN
}


#' Assigns genes to states just based on which mean is maximal
#'
#' @param matrix expression matrix
#' @param dynamic_object result of running findDynGenes
#' @param num_states num_states
#' @param pThresh pThresh
#' @param toScale toScale
#' @param key_word key_word
#'
#' @return data.frame of dynamically expressed genes, cluster, peakTime, ordered by peaktime
#' @export
assign_genes_to_states_simple <- function(
  matrix,
  dynamic_object,
  num_states = 3,
  pThresh = 0.01,
  toScale = FALSE,
  key_word = "state_"
) {
  dynamic_object <- normalize_state_data(dynamic_object)
  exp <- matrix[names(dynamic_object$genes[dynamic_object$genes < pThresh]), ]

  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }

  navg <- ceiling(ncol(exp) * 0.05)


  thresholds <- data.frame(gene = rownames(exp), thresh = rep(0, length(rownames(exp))))
  rownames(thresholds) <- thresholds$gene
  for (gene in rownames(exp)) {
    profile <- exp[gene, ][order(exp[gene, ], decreasing = FALSE)]
    bottom <- mean(profile[1:navg])
    top <- mean(profile[(length(profile) - navg):length(profile)])

    thresh <- ((top - bottom) * 0.33) + bottom
    thresholds[gene, "thresh"] <- thresh
  }


  t1 <- dynamic_object$cells$pseudotime
  names(t1) <- as.vector(dynamic_object$cells$cell_name)
  sort(t1, decreasing = FALSE)
  exp <- exp[, names(t1)]

  mean_expression <- data.frame(gene = character(), state = character(), mean_expression = numeric())


  state_names <- paste0(rep(key_word, num_states), seq(1:num_states))
  states <- vector("list", length(state_names))
  names(states) <- state_names


  ptmax <- max(dynamic_object$cells$pseudotime)
  ptmin <- min(dynamic_object$cells$pseudotime)
  chunk_size <- (ptmax - ptmin) / num_states

  cell_states <- rep("", length(names(t1)))
  names(cell_states) <- names(t1)

  for (i in 1:length(states)) {
    lower_bound <- ptmin + ((i - 1) * chunk_size)
    upper_bound <- ptmin + (i * chunk_size)
    chunk_cells <- rownames(dynamic_object$cells[dynamic_object$cells$pseudotime >= lower_bound & dynamic_object$cells$pseudotime <= upper_bound, ])
    chunk <- exp[, chunk_cells]

    chunk_df <- data.frame(means = rowMeans(chunk))
    chunk_df <- cbind(chunk_df, thresholds)
    chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

    states[[i]] <- rownames(chunk_df[chunk_df$active, ])
    genesPeakTimes <- apply(chunk, 1, which.max)
    gpt <- as.vector(dynamic_object$cells[chunk_cells, ][genesPeakTimes, ]$pseudotime)

    mean_expression <- rbind(
      mean_expression,
      data.frame(
        gene = rownames(chunk),
        state = rep(state_names[i], length(rownames(chunk))), mean_expression = rowMeans(chunk),
        peakTime = gpt
      )
    )
    cell_states[chunk_cells] <- state_names[i]
  }


  genes <- unique(as.vector(mean_expression$gene))
  cat("n genes: ", length(genes), "\n")
  gene_states <- rep("", length(genes))
  gene_state_peak_time <- rep(0, length(genes))
  state_mean <- rep(0, length(genes))

  names(gene_states) <- genes
  names(gene_state_peak_time) <- genes
  names(state_mean) <- genes
  for (gene in genes) {
    x <- mean_expression[mean_expression$gene == gene, ]
    xi <- which.max(x$mean_expression)
    gene_states[[gene]] <- as.vector(x[xi, ]$state)
    gene_state_peak_time[[gene]] <- as.vector(x[xi, ]$peakTime)
    state_mean[[gene]] <- max(x$mean_expression)
  }

  geneDF <- data.frame(
    gene = genes,
    state = gene_states,
    peakTime = gene_state_peak_time,
    state_mean = state_mean,
    pval = dynamic_object$genes[genes]
  )
  cells2 <- dynamic_object$cells[names(t1), ]
  cells2$state <- cell_states

  list(genes = geneDF, cells = cells2)
}
