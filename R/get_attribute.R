#' @include setGenerics.R
NULL

.get_attribute_impl <- function(
  object,
  celltypes = NULL,
  active_network = NULL,
  attribute = c(
    "genes",
    "tfs",
    "peaks",
    "regulators",
    "targets",
    "celltypes",
    "cells",
    "modules",
    "coefficients"
  ),
  ...
) {
  if (is.null(active_network)) {
    active_network <- DefaultNetwork(object)
  }

  nets_active <- GetNetwork(object, network = active_network)
  is_dynamic_network <- FALSE
  if (is.list(nets_active) && length(nets_active) > 0) {
    first_net <- nets_active[[1]]
    if (!is.null(first_net) &&
      methods::is(first_net, "Network") &&
      !is.null(first_net@params$network_type) &&
      identical(first_net@params$network_type, "dynamic")) {
      is_dynamic_network <- TRUE
    }
  }

  if (is.null(celltypes)) {
    if (is_dynamic_network && !is.null(nets_active) && length(nets_active) > 0) {
      celltypes <- names(nets_active)
    } else {
      celltypes <- .csn_celltypes(object)
    }
  }

  attribute <- match.arg(attribute)
  attrs_state <- .multicsn_get_attributes(object)
  attributes <- switch(
    EXPR = attribute,
    "genes" = lapply(
      celltypes,
      function(c) attrs_state[[c]]$genes$gene
    ),
    "peaks" = lapply(
      celltypes,
      function(c) attrs_state[[c]]$peaks$peak
    ),
    "regulators" = {
      networks <- GetNetwork(
        object,
        network = active_network,
        celltypes = celltypes
      )
      lapply(
        networks,
        function(net) {
          coeffs <- methods::slot(net, "coefficients")
          modules_obj <- methods::slot(net, "modules")
          modules_meta <- methods::slot(modules_obj, "meta")
          if (nrow(modules_meta) == 0) {
            unique(coeffs$tf)
          } else {
            unique(modules_meta$tf)
          }
        }
      )
    },
    "targets" = {
      networks <- GetNetwork(
        object,
        network = active_network,
        celltypes = celltypes
      )
      lapply(
        networks,
        function(net) {
          coeffs <- methods::slot(net, "coefficients")
          modules_obj <- methods::slot(net, "modules")
          modules_meta <- methods::slot(modules_obj, "meta")
          if (nrow(modules_meta) == 0) {
            unique(coeffs$target)
          } else {
            unique(modules_meta$target)
          }
        }
      )
    },
    "cells" = lapply(
      celltypes,
      function(c) attrs_state[[c]]$cells
    ),
    "modules" = lapply(
      celltypes,
      function(c) GetNetwork(object, network = active_network, celltypes = c)[[1]]@modules
    ),
    "coefficients" = lapply(
      celltypes,
      function(c) GetNetwork(object, network = active_network, celltypes = c)[[1]]@coefficients
    )
  )
  if (attribute == "tfs") {
    return(.multicsn_get_tfs(object))
  }
  if (attribute == "celltypes") {
    if (is_dynamic_network && !is.null(nets_active) && length(nets_active) > 0) {
      return(names(nets_active))
    }
    return(.csn_celltypes(object))
  }

  attributes <- purrr::set_names(attributes, celltypes)
  if (length(attributes) == 1) {
    return(attributes[[1]])
  }

  attributes
}

#' @param celltypes A character vector specifying the celltypes to get attributes for.
#' @param active_network A character string specifying the active network.
#' @param attribute Attribute to get: genes, tfs, peaks, regulators, targets, celltypes, cells, modules, coefficients.
#' @return A character vector or list of regulatory genes
#' @rdname get_attribute
#' @export
setMethod(
  "get_attribute",
  "Seurat",
  .get_attribute_impl
)

#' @rdname get_attribute
#' @export
setMethod(
  "get_attribute",
  "CSNObject",
  function(object, ...) {
    .get_attribute_impl(object, ...)
  }
)
