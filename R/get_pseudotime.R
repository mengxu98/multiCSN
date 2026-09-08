#' Get pseudotime information
#' @param object The input data (matrix, Seurat, etc.)
#' @param ... Arguments for other methods
#' @return Pseudotime information corresponding to cells
#' @export
get_pseudotime <- function(object, ...) {
  UseMethod("get_pseudotime")
}

#' @param meta_data Input meta data.
#' @param embeddings Embeddings information
#' @param cluster_column The column used for slingshot.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#' @param seed Random seed for reproducibility. Default is 1.
#' @param start_cluster The start cluster.
#' @param end_cluster The end cluster.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @rdname get_pseudotime
#' @export
#' @method get_pseudotime default
get_pseudotime.default <- function(
  object,
  meta_data = NULL,
  embeddings = NULL,
  cluster_column = "cluster",
  cores = 1,
  seed = 1,
  start_cluster = NULL,
  end_cluster = NULL,
  verbose = TRUE,
  ...
) {
  if (is.null(meta_data)) {
    cells <- rownames(object)
    if (is.null(cells)) {
      cells <- paste0("cell_", rep(1:nrow(object)))
      rownames(object) <- cells
    }
    meta_data <- data.frame(
      cells = cells,
      cluster = "cluster"
    )
    rownames(meta_data) <- cells
  } else {
    if (!cluster_column %in% colnames(meta_data)) {
      meta_data$cluster <- "cluster"
    } else {
      meta_data$cluster <- meta_data[, cluster_column]
    }
  }

  if (is.null(embeddings)) {
    set.seed(seed)
    embeddings <- suppressMessages(
      uwot::umap(object, n_threads = cores)
    )
  }

  thisutils::log_message("Running `slingshot`", verbose = verbose)
  pseudotime_res <- slingshot::slingshot(
    embeddings,
    clusterLabels = meta_data$cluster,
    start.clus = start_cluster,
    end.clus = end_cluster
  ) |> slingshot::slingPseudotime()
  pseudotime_res <- apply(
    pseudotime_res, 2, function(x) {
      thisutils::normalization(x)
    }
  ) |> as.data.frame()
  colnames(pseudotime_res) <- paste0("pseudotime_slingshot", 1:ncol(pseudotime_res))
  pseudotime_res <- cbind.data.frame(cluster = meta_data$cluster, pseudotime_res)

  return(pseudotime_res)
}

#' @param assay The assay used for slingshot.
#' @param layer The layer used for slingshot.
#' @param reduction The reduction used for slingshot.
#' @rdname get_pseudotime
#' @export
#' @method get_pseudotime Seurat
get_pseudotime.Seurat <- function(
  object,
  assay = "RNA",
  layer = "data",
  cluster_column = "cluster",
  reduction = "umap",
  start_cluster = NULL,
  end_cluster = NULL,
  ...
) {
  embeddings <- Seurat::Embeddings(object, reduction = reduction)

  result <- get_pseudotime(
    Seurat::GetAssay(object, layer = layer),
    meta_data = object@meta.data,
    embeddings = embeddings,
    cluster_column = cluster_column,
    start_cluster = start_cluster,
    end_cluster = end_cluster,
    ...
  )
  result <- result[colnames(object), ]
  object <- Seurat::AddMetaData(object, result)

  return(object)
}
