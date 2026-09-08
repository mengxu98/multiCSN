#' @include setClass.R
NULL

#' @title Get any attribute from a CSNObject object
#'
#' @param object A CSNObject object
#' @param ... Additional arguments
#'
#' @return A character vector of regulatory genes
#' @export
#'
#' @rdname get_attribute
setGeneric(
  name = "get_attribute",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("get_attribute")
  }
)

#' @title Get active network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname DefaultNetwork
setGeneric(
  name = "DefaultNetwork",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("DefaultNetwork")
  }
)

#' @title Get summary of seurat assay
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname GetAssaySummary
setGeneric(
  name = "GetAssaySummary",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("GetAssaySummary")
  }
)

#' @title Get network
#'
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname GetNetwork
setGeneric(
  name = "GetNetwork",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("GetNetwork")
  }
)

#' @title Get network graph
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname NetworkGraph
setGeneric(
  name = "NetworkGraph",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("NetworkGraph")
  }
)

#' @title Get gene ranks from Network/CSNObject
#' @param object The input data.
#' @param ... Arguments for other methods
#' @rdname GeneRanks
#' @export
setGeneric(
  name = "GeneRanks",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("GeneRanks")
  }
)

#' @title Get stored shortest-path results from Network/CSNObject
#' @param object The input data.
#' @param ... Arguments for other methods
#' @rdname ShortestPaths
#' @export
setGeneric(
  name = "ShortestPaths",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("ShortestPaths")
  }
)

#' @title Get stored edge_uniqueness results from CSNObject
#' @param object The input data.
#' @param ... Arguments for other methods
#' @rdname EdgeUniqueness
#' @export
setGeneric(
  name = "EdgeUniqueness",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("EdgeUniqueness")
  }
)

#' @title Get stored extract_genes_transition results from CSNObject
#' @param object The input data.
#' @param ... Arguments for other methods
#' @rdname ExtractGenesTransition
#' @export
setGeneric(
  name = "ExtractGenesTransition",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("ExtractGenesTransition")
  }
)

#' @title Get network modules
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname NetworkModules
setGeneric(
  name = "NetworkModules",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("NetworkModules")
  }
)

#' @title Get network parameters
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname NetworkParams
setGeneric(
  name = "NetworkParams",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("NetworkParams")
  }
)

#' @title Get network regions
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname NetworkRegions
setGeneric(
  name = "NetworkRegions",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("NetworkRegions")
  }
)

#' @title Get network TFs
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname NetworkTFs
setGeneric(
  name = "NetworkTFs",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("NetworkTFs")
  }
)

#' @title Get parameters
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname Params
setGeneric(
  name = "Params",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("Params")
  }
)

#' @title Embed perturbation
#'
#' @param object The input data.
#'
#' @export
#'
#' @rdname embedPerturbation
setGeneric(
  name = "embedPerturbation",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("embedPerturbation")
  }
)

#' @title Export CSN
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname export_csn
setGeneric(
  name = "export_csn",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("export_csn")
  }
)

#' @title Get network graph
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname get_network_graph
setGeneric(
  name = "get_network_graph",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("get_network_graph")
  }
)

#' @title Get TF network
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname get_tf_network
setGeneric(
  name = "get_tf_network",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("get_tf_network")
  }
)

#' @title Initiate object
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname initiate_object
setGeneric(
  name = "initiate_object",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("initiate_object")
  }
)

#' @title Metrics
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname metrics
setGeneric(
  name = "metrics",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("metrics")
  }
)

#' @title Plot perturbation
#'
#' @param object The input data.
#'
#' @export
#'
#' @rdname plotPerturbation
setGeneric(
  name = "plotPerturbation",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plotPerturbation")
  }
)

#' @title Plot perturbation trajectory
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname plotPerturbationTrajectory
setGeneric(
  name = "plotPerturbationTrajectory",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plotPerturbationTrajectory")
  }
)

#' @title Plot goodness of fit
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname plot_gof
setGeneric(
  name = "plot_gof",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_gof")
  }
)

#' @title Plot module metrics
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname plot_module_metrics
setGeneric(
  name = "plot_module_metrics",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_module_metrics")
  }
)

#' @title Plot network graph
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname plot_network_graph
setGeneric(
  name = "plot_network_graph",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_network_graph")
  }
)

#' @title Plot TF network
#'
#' @param object The input data.
#' @param ... Arguments for other methods
#'
#' @export
#'
#' @rdname plot_tf_network
setGeneric(
  name = "plot_tf_network",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("plot_tf_network")
  }
)

#' @title Predict perturbation
#'
#' @param object The input data.
#'
#' @export
#'
#' @rdname predictPerturbation
setGeneric(
  name = "predictPerturbation",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("predictPerturbation")
  }
)
#' @title Initiate the \code{CSNObject} object
#' @param object The input data, a seurat object.
#' @param ... Arguments for other methods
#' @rdname initiate_object
#' @export initiate_object
setGeneric(
  "initiate_object",
  signature = "object",
  function(object, ...) {
    standardGeneric("initiate_object")
  }
)

#' @title Get network regions
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname NetworkRegions
#' @export NetworkRegions
setGeneric(
  "NetworkRegions",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkRegions")
  }
)

#' @title Get network modules
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname NetworkModules
#' @export NetworkModules
setGeneric(
  "NetworkModules",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkModules")
  }
)

#' @title Get network parameters
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname NetworkParams
#' @export NetworkParams
setGeneric(
  "NetworkParams",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkParams")
  }
)

#' @title Get network graph
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname NetworkGraph
#' @export NetworkGraph
setGeneric(
  "NetworkGraph",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkGraph")
  }
)

#' @title Get network
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname GetNetwork
#' @export GetNetwork
setGeneric(
  "GetNetwork",
  signature = "object",
  function(object, ...) {
    standardGeneric("GetNetwork")
  }
)

#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname Params
#' @export Params
setGeneric(
  "Params",
  signature = "object",
  function(object, ...) {
    standardGeneric("Params")
  }
)

#' @title Get network TFs
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname NetworkTFs
#' @export NetworkTFs
setGeneric(
  "NetworkTFs",
  signature = "object",
  function(object, ...) {
    standardGeneric("NetworkTFs")
  }
)

#' @title Get any attribute from a CSNObject
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname get_attribute
#' @export get_attribute
setGeneric(
  "get_attribute",
  signature = "object",
  function(object, ...) {
    standardGeneric("get_attribute")
  }
)

#' @title Get active network name
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname DefaultNetwork
#' @export DefaultNetwork
setGeneric(
  "DefaultNetwork",
  signature = "object",
  function(object, ...) {
    standardGeneric("DefaultNetwork")
  }
)

#' @param value Name of the active network to set.
#' @rdname DefaultNetwork
#' @export DefaultNetwork<-
setGeneric(
  "DefaultNetwork<-",
  signature = "object",
  function(object, ..., value) {
    standardGeneric("DefaultNetwork<-")
  }
)

#' @title Get metrics
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname metrics
#' @export metrics
setGeneric(
  "metrics",
  signature = "object",
  function(object, ...) {
    standardGeneric("metrics")
  }
)

#' @title Get summary of seurat assay
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname GetAssaySummary
#' @export GetAssaySummary
setGeneric(
  "GetAssaySummary",
  signature = "object",
  function(object, ...) {
    standardGeneric("GetAssaySummary")
  }
)

#' @title Export network from CSN object
#' @description Export network data from a CSN object with optional filtering
#' @param object A CSNObject object
#' @param ... Additional arguments
#' @return A list of network data frames by cell type
#' @export
setGeneric(
  "export_csn",
  function(object, ...) {
    standardGeneric("export_csn")
  }
)
