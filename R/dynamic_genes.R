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
define_states_new <- function(
  dynamic_object,
  matrix,
  method = "pseudotime",
  num_states = 2,
  pseudotime_cuts = NULL,
  group_assignments = NULL,
  p_value = 0.05,
  winSize = 2
) {
  define_states(
    dynamic_object = dynamic_object, matrix = matrix, method = method,
    num_states = num_states, pseudotime_cuts = pseudotime_cuts,
    group_assignments = group_assignments, p_value = p_value, winSize = winSize
  )
}

#' Assigns genes to states
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_states
#' @param method method of assigning state genes, either "active_expression" (looks for active expression in state) or "DE" (looks for differentially expressed genes per state)
#' @param p_value p value threshold if gene is dynamically expressed
#' @param pThresh_DE p value if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into state with max expression.
#'
#' @return states a list detailing genes active in each state
#' @export
assign_genes_to_states_new1 <- function(
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
  dynamic_object <- dynamic_object[[2]]

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


  state_names <- unique(dynamic_object$cells$group)
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
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == state, "cells"]
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
        res <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, res)
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
assign_network <- function(
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
  message("Starting.")

  dynamic_object <- dynamic_object[[2]]

  dynamic_genes <- dynamic_object$genes
  dynamic_genes <- dynamic_genes[dynamic_genes$dynamic_genes < p_value, ]$genes
  exp <- Matrix::as.matrix(
    matrix[intersect(rownames(matrix), dynamic_genes), ]
  )

  if (toScale) {
    if (class(exp)[1] != "matrix") {
      exp <- t(scale(Matrix::t(exp)))
    } else {
      exp <- t(scale(t(exp)))
    }
  }


  state_names <- unique(dynamic_object$cells$group)
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
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == state, "cells"]

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
        res <- data.frame(gene = gene, mean_diff = (t$coefficient[1] - t$coefficient[2]), pval = t$p.value)
        diffres <- rbind(diffres, res)
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
