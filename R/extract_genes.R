#' @include setClass.R
#' @include setGenerics.R
NULL

extract_genes <- function(
  network_list,
  regulators = NULL,
  targets = NULL,
  top_regulators_num = 10,
  top_targets_num = 10,
  common_regulators = FALSE
) {
  if (common_regulators) {
    regulators_list <- purrr::map(
      network_list, function(x) {
        regulators <- x$regulator
        regulators[!duplicated(regulators)]
      }
    )
    regulators_common <- purrr::reduce(
      regulators_list,
      .f = intersect
    )
  }

  if (is.list(network_list)) {
    genes_list <- purrr::map(
      network_list, function(x) {
        regulators_list <- x$regulator
        regulators_list <- regulators_list[!duplicated(regulators_list)]

        if (common_regulators) {
          if (is.null(regulators)) {
            regulators <- regulators_common[1:min(top_regulators_num, length(regulators_common))]
          }
        } else {
          if (is.null(regulators)) {
            regulators <- regulators_list[1:min(top_regulators_num, length(regulators_list))]
          }
        }

        targets_list <- purrr::map_dfr(
          regulators, function(g) {
            x_g <- inferCSN::network_format(
              x,
              regulators = g,
              targets = targets,
              abs_weight = FALSE
            )
            x_g[1:top_targets_num, ]
          }
        )
        targets <- unique(targets_list$target)

        return(list(
          regulators = regulators,
          targets = targets
        ))
      }
    )
  }

  return(genes_list)
}


setdiff_genes <- function(list) {
  genes_list <- list()
  for (i in seq_len(length(list))) {
    genes1 <- list[[i]]
    nums_list <- seq_len(length(list))[-i]
    list_new <- list[nums_list]
    genes2 <- purrr::list_c(list_new)
    genes_list[[i]] <- unique(setdiff(genes1, genes2))
  }

  return(genes_list)
}


.extract_genes_transition_list <- function(
  network_list,
  regulators = NULL,
  targets = NULL,
  transition_networks = NULL,
  top_regulators_num = 5,
  top_targets_num = 5,
  common_regulators = TRUE,
  common_targets = TRUE
) {
  if (is.null(transition_networks)) {
    network_names <- sort(names(network_list))
    transition_networks <- purrr::map(
      seq_len(((length(network_names) - 1) / 2)),
      .f = function(x) {
        x <- 2 * x - 1
        network_names[x:(x + 2)]
      }
    )
  }

  transition_networks_list <- purrr::map(
    transition_networks,
    .f = function(x) {
      network_list_sub <- network_list[x]
      transition_networks_list_sub <- extract_genes(
        network_list_sub,
        common_regulators = common_regulators,
        top_regulators_num = top_regulators_num,
        top_targets_num = top_targets_num
      )
      if (common_regulators && common_targets) {
        targets_list <- purrr::map(
          transition_networks_list_sub, function(x) {
            targets <- x$target
            targets[!duplicated(targets)]
          }
        )
        targets_common <- purrr::reduce(
          targets_list,
          .f = intersect
        )

        targets_setdiff <- setdiff_genes(targets_list)
        names(targets_setdiff) <- names(targets_list)

        transition_networks_list_sub <- purrr::map2(
          transition_networks_list_sub,
          targets_setdiff,
          .f = function(x, y) {
            list(regulators = x$regulators, targets = y)
          }
        )
        transition_networks_list_sub[["intersect"]] <- list(
          regulators = transition_networks_list_sub[[1]][1],
          targets = targets_common
        )
      }

      return(transition_networks_list_sub)
    }
  )
  names(transition_networks_list) <- paste0(
    "transition_network",
    seq_len(length(transition_networks))
  )

  transition_networks_list
}

#' @title Extract transition-specific regulators and targets
#' @name extract_genes_transition
#' @param object Named list of network tables or CSNObject.
#' @param ... Passed to methods.
#' @rdname extract_genes_transition
#' @export
setGeneric(
  name = "extract_genes_transition",
  signature = c("object"),
  def = function(object, ...) {
    standardGeneric("extract_genes_transition")
  }
)

#' @param regulators Optional regulator subset; kept for backward compatibility
#'   and currently not used in the ranking.
#' @param targets Optional target subset; kept for backward compatibility
#'   and currently not used in the ranking.
#' @param transition_networks Optional list of state triplets defining transitions;
#'   \code{NULL} derives triplets from consecutive states.
#' @param top_regulators_num Number of top regulators kept per state.
#' @param top_targets_num Number of top targets kept per regulator.
#' @param common_regulators Logical; when \code{TRUE}, restrict regulators to those
#'   shared across the states of a transition.
#' @param common_targets Logical; when \code{TRUE}, split targets into
#'   transition-specific sets and a shared \code{intersect} set.
#' @rdname extract_genes_transition
#' @export
setMethod(
  "extract_genes_transition",
  signature(object = "list"),
  function(object,
           regulators = NULL,
           targets = NULL,
           transition_networks = NULL,
           top_regulators_num = 5,
           top_targets_num = 5,
           common_regulators = TRUE,
           common_targets = TRUE) {
    .extract_genes_transition_list(
      object,
      regulators,
      targets,
      transition_networks,
      top_regulators_num,
      top_targets_num,
      common_regulators,
      common_targets
    )
  }
)

#' @param network Name of dynamic network in \code{object@networks}.
#' @param celltypes States to include; \code{NULL} for all (sorted by name for transition triplets).
#' @param weight_cutoff Passed to \code{export_csn}.
#' @rdname extract_genes_transition
#' @export
setMethod(
  "extract_genes_transition",
  signature(object = "Seurat"),
  function(object,
           network = NULL,
           celltypes = NULL,
           weight_cutoff = NULL,
           regulators = NULL,
           targets = NULL,
           transition_networks = NULL,
           top_regulators_num = 5,
           top_targets_num = 5,
           common_regulators = TRUE,
           common_targets = TRUE) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "dynamic",
      verbose = TRUE,
      caller = "extract_genes_transition"
    )
    nets <- GetNetwork(object, network = network)
    if (is.null(nets) || !is.list(nets) || length(nets) == 0) {
      stop("No networks found for the specified network.")
    }
    if (is.null(celltypes)) celltypes <- names(nets)
    celltypes <- intersect(celltypes, names(nets))
    network_list <- lapply(celltypes, function(sid) {
      df <- tryCatch(export_csn(nets[[sid]], weight_cutoff = weight_cutoff), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      as.data.frame(df)
    })
    ok <- !sapply(network_list, is.null)
    network_list <- network_list[ok]
    names(network_list) <- celltypes[ok]
    res <- .extract_genes_transition_list(
      network_list,
      regulators,
      targets,
      transition_networks,
      top_regulators_num,
      top_targets_num,
      common_regulators,
      common_targets
    )

    state <- .read_multicsn_state(object, init = TRUE)
    extract_res <- state$extract_genes_transition %ss% list()
    extract_res[[network]] <- list(
      method = "extract_genes_transition",
      celltypes = names(network_list),
      transition_networks = transition_networks,
      top_regulators_num = top_regulators_num,
      top_targets_num = top_targets_num,
      common_regulators = common_regulators,
      common_targets = common_targets,
      result = res
    )
    state$extract_genes_transition <- extract_res
    .write_multicsn_state(object, state)
  }
)

#' @rdname extract_genes_transition
#' @export
setMethod(
  "extract_genes_transition",
  signature(object = "CSNObject"),
  function(object, ...) {
    .stop_csnobject_runtime()
  }
)
