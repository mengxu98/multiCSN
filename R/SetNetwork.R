#' @include setClass.R
#' @title Set Network object
#' @description Set or replace a network in a CSNObject.
#' @param object CSNObject
#' @param network Network object to set
#' @param name Character, name of the network (optional)
#' @return Updated CSNObject
#' @export
SetNetwork <- function(object, network, name = NULL) {
  if (methods::is(object, "Seurat")) {
    name <- name %ss% DefaultNetwork(object)
    if (is.null(name) || !nzchar(name)) {
      stop("No active network specified")
    }
    return(.multicsn_set_networks(object, name, network))
  }

  if (is.null(name)) {
    name <- object@active_network
  }

  if (is.null(name)) {
    stop("No active network specified")
  }

  object@networks[[name]] <- network
  object@active_network <- name
  return(object)
}
