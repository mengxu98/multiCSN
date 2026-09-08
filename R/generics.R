#' @title Plot goodness-of-fit
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname plot_gof
#' @export
setGeneric("plot_gof", function(object, ...) standardGeneric("plot_gof"))

#' @title Plot module metrics
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname plot_module_metrics
#' @export
setGeneric("plot_module_metrics", function(object, ...) standardGeneric("plot_module_metrics"))

#' @title Get network graph
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname get_network_graph
#' @export
setGeneric("get_network_graph", function(object, ...) standardGeneric("get_network_graph"))

#' @title Plot network graph
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname plot_network_graph
#' @export
setGeneric("plot_network_graph", function(object, ...) standardGeneric("plot_network_graph"))

#' @title Get TF network
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname get_tf_network
#' @export
setGeneric("get_tf_network", function(object, ...) standardGeneric("get_tf_network"))

#' @title Plot TF network
#' @param object CSNObject or Network
#' @param ... Additional arguments
#' @rdname plot_tf_network
#' @export
setGeneric("plot_tf_network", function(object, ...) standardGeneric("plot_tf_network"))
