#' Define epochs
#'
#' @param dynamic_object results of running findDynGenes, or a list of results of running findDynGenes per path. If list, names should match names of matrix.
#' @param matrix genes-by-cells expression matrix, or a list of expression matrices per path. If list, names should match names of dynamic_object.
#' @param method method to define epochs. Either "pseudotime", "cell_order", "group", "con_similarity", "kmeans", "hierarchical"
#' @param num_epochs number of epochs to define. Ignored if epoch_transitions, pseudotime_cuts, or group_assignments are provided.
#' @param pseudotime_cuts vector of pseudotime cutoffs. If NULL, cuts are set to max(pseudotime)/num_epochs.
#' @param group_assignments a list of vectors where names(assignment) are epoch names, and vectors contain groups belonging to corresponding epoch
#' @param p_value p_value
#' @param winSize winSize
#'
#' @return updated list of dynamic_object with epoch column included in dynamic_object$cells
#' @export
define_epochs_new <- function(
  dynamic_object,
  matrix,
  method = "pseudotime",
  num_epochs = 2,
  pseudotime_cuts = NULL,
  group_assignments = NULL,
  p_value = 0.05,
  winSize = 2
) {
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
    path_dyn$cells$epoch <- NA


    if (method == "pseudotime") {
      if (is.null(pseudotime_cuts)) {
        pseudotime_cuts <- seq(
          min(path_dyn$cells$pseudotime),
          max(path_dyn$cells$pseudotime),
          (max(path_dyn$cells$pseudotime) - min(path_dyn$cells$pseudotime)) / num_epochs
        )
        pseudotime_cuts <- pseudotime_cuts[-length(pseudotime_cuts)]
        pseudotime_cuts <- pseudotime_cuts[-1]
      }

      path_dyn <- split_epochs_by_pseudotime(path_dyn, pseudotime_cuts)
    }


    if (method == "cell_order") {
      t1 <- path_dyn$cells$pseudotime
      names(t1) <- as.vector(path_dyn$cells$cell_name)

      chunk_size <- floor(length(t1) / num_epochs)

      for (i in 1:num_epochs) {
        if (i == num_epochs) {
          cells_in_epoch <- names(t1)[(1 + ((i - 1) * chunk_size)):length(t1)]
        } else {
          cells_in_epoch <- names(t1)[(1 + ((i - 1) * chunk_size)):(i * chunk_size)]
        }
        path_dyn$cells$epoch[path_dyn$cells$cell_name %in% cells_in_epoch] <- paste0("epoch", i)
      }
    }


    if (method == "group") {
      if (is.null(group_assignments)) {
        stop("Must provide group_assignments for group method.")
      }

      path_dyn <- split_epochs_by_group(path_dyn, group_assignments)
    }


    if (method == "con_similarity") {
      cuts <- find_cuts_by_similarity(matrix[[path]], path_dyn, winSize = winSize, p_value = p_value)
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }


    if (method == "kmeans") {
      cuts <- find_cuts_by_clustering(
        matrix[[path]],
        path_dyn,
        num_epochs = num_epochs,
        method = "kmeans",
        p_value = p_value
      )
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }


    if (method == "hierarchical") {
      cuts <- find_cuts_by_clustering(matrix[[path]], path_dyn, num_epochs = num_epochs, method = "hierarchical", p_value = p_value)
      path_dyn <- split_epochs_by_pseudotime(path_dyn, cuts)
    }


    new_dynRes[[path]] <- path_dyn
  }

  if (length(new_dynRes) == 1) {
    new_dynRes <- new_dynRes[[1]]
  }
  new_dynRes
}

#' Assigns genes to epochs
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_epochs
#' @param method method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)
#' @param p_value p value threshold if gene is dynamically expressed
#' @param pThresh_DE p value if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.
#'
#' @return epochs a list detailing genes active in each epoch
#' @export
assign_epochs_new1 <- function(
  matrix,
  dynamic_object,
  method = "active_expression",
  p_value = 0.05,
  pThresh_DE = 0.05,
  active_thresh = 0.33,
  toScale = FALSE,
  forceGenes = TRUE
) {
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


  epoch_names <- unique(dynamic_object$cells$group)
  epochs <- vector("list", length(epoch_names))
  names(epochs) <- epoch_names

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
    epoch = numeric(),
    mean_expression = numeric()
  )
  if (method == "active_expression") {
    for (epoch in names(epochs)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == epoch, "cells"]
      chunk <- exp[, chunk_cells]

      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- rownames(chunk_df[chunk_df$active, ])

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  } else {
    for (epoch in names(epochs)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
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

      epochs[[epoch]] <- diffres$gene[diffres$padj < pThresh_DE]


      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- intersect(epochs[[epoch]], rownames(chunk_df[chunk_df$active, ]))

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  }


  if (forceGenes) {
    assignedGenes <- unique(unlist(epochs))
    orphanGenes <- setdiff(rownames(exp), assignedGenes)
    message("There are ", length(orphanGenes), " orphan genes\n")
    for (oGene in orphanGenes) {
      xdat <- mean_expression[mean_expression$gene == oGene, ]
      ep <- xdat[which.max(xdat$mean_expression), ]$epoch
      epochs[[ep]] <- append(epochs[[ep]], oGene)
    }
  }

  epochs$mean_expression <- mean_expression

  epochs
}

#' Assigns genes to epochs
#'
#' @param matrix genes-by-cells expression matrix
#' @param dynamic_object individual path result of running define_epochs
#' @param method method of assigning epoch genes, either "active_expression" (looks for active expression in epoch) or "DE" (looks for differentially expressed genes per epoch)
#' @param p_value pval threshold if gene is dynamically expressed
#' @param pThresh_DE pval if gene is differentially expressed. Ignored if method is active_expression.
#' @param active_thresh value between 0 and 1. Percent threshold to define activity
#' @param toScale whether or not to scale the data
#' @param forceGenes whether or not to rescue orphan dyanmic genes, forcing assignment into epoch with max expression.
#'
#' @return epochs a list detailing genes active in each epoch
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


  epoch_names <- unique(dynamic_object$cells$group)
  epochs <- vector("list", length(epoch_names))
  names(epochs) <- epoch_names

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
    epoch = numeric(),
    mean_expression = numeric()
  )
  if (method == "active_expression") {
    for (epoch in names(epochs)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$group == epoch, "cells"]

      chunk <- exp[, chunk_cells]

      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- rownames(chunk_df[chunk_df$active, ])

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  } else {
    for (epoch in names(epochs)) {
      chunk_cells <- dynamic_object$cells[dynamic_object$cells$epoch == epoch, "cells"]
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

      epochs[[epoch]] <- diffres$gene[diffres$padj < pThresh_DE]


      chunk_df <- data.frame(means = rowMeans(chunk))
      chunk_df <- cbind(chunk_df, thresholds)
      chunk_df$active <- (chunk_df$means >= chunk_df$thresh)

      epochs[[epoch]] <- intersect(epochs[[epoch]], rownames(chunk_df[chunk_df$active, ]))

      mean_expression <- rbind(
        mean_expression,
        data.frame(
          gene = rownames(chunk),
          epoch = rep(epoch, length(rownames(chunk))),
          mean_expression = rowMeans(chunk)
        )
      )
    }
  }


  if (forceGenes) {
    assignedGenes <- unique(unlist(epochs))
    orphanGenes <- setdiff(rownames(exp), assignedGenes)
    message("There are ", length(orphanGenes), " orphan genes\n")
    for (oGene in orphanGenes) {
      xdat <- mean_expression[mean_expression$gene == oGene, ]
      ep <- xdat[which.max(xdat$mean_expression), ]$epoch
      epochs[[ep]] <- append(epochs[[ep]], oGene)
    }
  }

  epochs$mean_expression <- mean_expression

  epochs
}
