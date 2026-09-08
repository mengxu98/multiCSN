.csn_celltypes <- function(object) {
  .multicsn_get_celltypes(object)
}

.csn_get_assay <- function(object, assay = NULL) {
  if (methods::is(object, "CSNObject")) {
    return(Seurat::GetAssay(object@data, assay = assay))
  }
  if (methods::is(object, "Seurat")) {
    return(Seurat::GetAssay(object, assay = assay))
  }
  stop(".csn_get_assay: object must be a Seurat or CSNObject object.", call. = FALSE)
}

.csn_layer_data <- function(object, ...) {
  if (methods::is(object, "CSNObject")) {
    return(SeuratObject::LayerData(object@data, ...))
  }
  if (methods::is(object, "Seurat")) {
    return(SeuratObject::LayerData(object, ...))
  }
  stop(".csn_layer_data: object must be a Seurat or CSNObject object.", call. = FALSE)
}

#' Aggregate Seurat assay over groups
#' @param object The csn object.
#' @param group_name A character vector indicating the metadata column to aggregate over.
#' @param fun The summary function to be applied to each group.
#' @param assay The assay to summarize.
#' @param layer The layer to summarize.
#' @return A Seurat object.
#' @export
aggregate_assay <- function(
  object,
  group_name,
  fun = "mean",
  assay = "RNA",
  layer = "data"
) {
  ass_mat <- Matrix::t(
    Seurat::GetAssayData(
      object,
      assay = assay,
      layer = layer
    )
  )
  groups <- as.character(object@meta.data[[group_name]])
  agg_mat <- thisutils::aggregate_matrix(
    ass_mat,
    groups = groups,
    fun = fun
  )
  if (is.null(object@assays[[assay]]@misc$summary)) {
    object@assays[[assay]]@misc$summary <- list()
  }
  object@assays[[assay]]@misc$summary[[group_name]] <- agg_mat
  return(object)
}
