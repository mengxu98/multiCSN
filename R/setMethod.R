#' @include setClass.R
#' @include setGenerics.R
NULL

#' @title Get GRN inference parameters
#' @rdname Params
#' @export
setMethod(
  f = "Params",
  signature = "Seurat",
  definition = function(object, ...) {
    .multicsn_get(object, "params", default = list())
  }
)

setMethod(
  f = "Params",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@params)
  }
)

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#' @rdname GetNetwork
#' @export
setMethod(
  f = "GetNetwork",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    network <- .multicsn_resolve_network(
      object,
      network = network,
      celltypes = celltypes,
      preferred = "active",
      verbose = FALSE,
      caller = "GetNetwork"
    )
    nets <- .multicsn_get_networks(object, network = network)
    if (is.null(nets)) {
      return(NULL)
    }
    if (!is.null(celltypes)) {
      out <- lapply(
        celltypes,
        function(x) nets[[x]]
      )
      names(out) <- celltypes
      return(out)
    }
    nets
  }
)

setMethod(
  f = "GetNetwork",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    if (!is.null(celltypes)) {
      return(
        lapply(
          celltypes,
          function(x) {
            object@networks[[network]][[x]]
          }
        )
      )
    }
    return(object@networks[[network]])
  }
)

#' @rdname NetworkTFs
#' @export
setMethod(
  f = "NetworkTFs",
  signature = "Seurat",
  definition = function(object, ...) {
    NetworkRegions(object)@motifs2tfs
  }
)

setMethod(
  f = "NetworkTFs",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@regions@motifs2tfs)
  }
)

#' @rdname NetworkRegions
#' @export
setMethod(
  f = "NetworkRegions",
  signature = "Seurat",
  definition = function(object, ...) {
    .multicsn_get_regions(object)
  }
)

setMethod(
  f = "NetworkRegions",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@regions)
  }
)

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#' @rdname NetworkModules
#' @export
setMethod(
  f = "NetworkModules",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(networks)) {
      return(NULL)
    }
    lapply(networks, function(net) NetworkModules(net))
  }
)

setMethod(
  f = "NetworkModules",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkModules(net)
        )
      )
    }
    return(
      lapply(
        networks,
        function(net) NetworkModules(net)
      )
    )
  }
)

#' @rdname NetworkModules
#' @export
setMethod(
  f = "NetworkModules",
  signature = "Network",
  definition = function(object, ...) {
    return(object@modules)
  }
)

#' @rdname NetworkParams
#' @export
setMethod(
  f = "NetworkParams",
  signature = "Network",
  definition = function(object, ...) {
    return(object@params)
  }
)

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#' @rdname NetworkParams
#' @export
setMethod(
  f = "NetworkParams",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(networks)) {
      return(NULL)
    }
    lapply(networks, function(net) NetworkParams(net))
  }
)

setMethod(
  f = "NetworkParams",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkParams(net)
        )
      )
    }
    return(
      lapply(
        networks, function(net) NetworkParams(net)
      )
    )
  }
)

#' @param network network
#' @param graph graph
#' @param celltypes cell types to analyze, NULL for all cell types
#' @rdname NetworkGraph
#' @export
setMethod(
  f = "NetworkGraph",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        graph = "module_graph",
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(networks)) {
      return(NULL)
    }
    lapply(networks, function(net) NetworkGraph(net, graph = graph))
  }
)

setMethod(
  f = "NetworkGraph",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        graph = "module_graph",
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(celltypes)) {
      return(
        lapply(
          networks,
          function(net) NetworkGraph(net, graph = graph)
        )
      )
    }
    return(
      lapply(
        networks,
        function(net) NetworkGraph(net, graph = graph)
      )
    )
  }
)

#' @rdname NetworkGraph
#' @export
setMethod(
  f = "NetworkGraph",
  signature = "Network",
  definition = function(object,
                        graph = "module_graph",
                        ...) {
    if (!graph %in% names(object@graphs)) {
      stop(
        paste0(
          "The requested graph '", graph, "' does not exist. ",
          "Try (re-)running `get_network_graph()`."
        )
      )
    }
    return(object@graphs[[graph]])
  }
)

#' @param network Name of the network.
#' @param celltypes Cell types to get ranks for; NULL for all.
#' @rdname GeneRanks
#' @export
setMethod(
  f = "GeneRanks",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(object, network = network, celltypes = celltypes)
    if (is.null(networks)) {
      return(NULL)
    }
    out <- lapply(networks, function(net) GeneRanks(net, ...))
    if (!is.null(names(networks))) names(out) <- names(networks)
    out
  }
)

setMethod(
  f = "GeneRanks",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    networks <- GetNetwork(object, network = network, celltypes = celltypes)
    if (is.null(networks)) {
      return(NULL)
    }
    out <- lapply(networks, function(net) GeneRanks(net, ...))
    if (!is.null(names(networks))) names(out) <- names(networks)
    out
  }
)

#' @rdname GeneRanks
#' @export
setMethod(
  f = "GeneRanks",
  signature = "Network",
  definition = function(object, ...) {
    object@params$gene_rank
  }
)

#' @param network Name of the network.
#' @param celltypes Cell types/states to get paths for; NULL for all.
#' @rdname ShortestPaths
#' @export
setMethod(
  f = "ShortestPaths",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    state <- .read_multicsn_state(object)
    res <- state$shortest_paths %ss% NULL
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    out <- res[[network]]
    if (!is.null(celltypes) && is.list(out) && !is.null(names(out))) {
      out <- out[intersect(celltypes, names(out))]
    }
    out
  }
)

setMethod(
  f = "ShortestPaths",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    res <- tryCatch(object@metadata$shortest_paths, error = function(e) NULL)
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    out <- res[[network]]
    if (!is.null(celltypes) && is.list(out) && !is.null(names(out))) {
      out <- out[intersect(celltypes, names(out))]
    }
    out
  }
)

#' @rdname ShortestPaths
#' @export
setMethod(
  f = "ShortestPaths",
  signature = "Network",
  definition = function(object, ...) {
    object@params$shortest_paths
  }
)

#' @param network Name of the network.
#' @rdname EdgeUniqueness
#' @export
setMethod(
  f = "EdgeUniqueness",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        ...) {
    state <- .read_multicsn_state(object)
    res <- state$edge_uniqueness %ss% NULL
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    res[[network]]
  }
)

setMethod(
  f = "EdgeUniqueness",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        ...) {
    res <- tryCatch(object@metadata$edge_uniqueness, error = function(e) NULL)
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    res[[network]]
  }
)

#' @param network Name of the stored network whose transition results are returned;
#'   defaults to \code{DefaultNetwork(object)}.
#' @rdname ExtractGenesTransition
#' @export
setMethod(
  f = "ExtractGenesTransition",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        ...) {
    state <- .read_multicsn_state(object)
    res <- state$extract_genes_transition %ss% NULL
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    res[[network]]
  }
)

setMethod(
  f = "ExtractGenesTransition",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        ...) {
    res <- tryCatch(object@metadata$extract_genes_transition, error = function(e) NULL)
    if (is.null(res) || is.null(res[[network]])) {
      return(NULL)
    }
    res[[network]]
  }
)

#' @rdname DefaultNetwork
#' @export
setMethod(
  f = "DefaultNetwork",
  signature = "Seurat",
  definition = function(object, ...) {
    .multicsn_get_active_network(object)
  }
)

setMethod(
  f = "DefaultNetwork",
  signature = "CSNObject",
  definition = function(object, ...) {
    return(object@active_network)
  }
)

#' @param value Name of the network to activate.
#' @rdname DefaultNetwork
#' @export
setReplaceMethod(
  f = "DefaultNetwork",
  signature = "Seurat",
  definition = function(object, value) {
    .multicsn_set_active_network(object, value)
  }
)

setReplaceMethod(
  f = "DefaultNetwork",
  signature = "CSNObject",
  definition = function(object, value) {
    value <- value[1] %ss% ""
    if (is.na(value) || !nzchar(value)) {
      stop("DefaultNetwork<-: value must be a non-empty character string")
    }
    if (!value %in% names(object@networks)) {
      stop(
        sprintf(
          "DefaultNetwork<-: network '%s' not found in object@networks",
          value
        )
      )
    }
    object@active_network <- value
    return(object)
  }
)

#' @param group_name group_name
#' @param assay assay
#' @param verbose verbose
#' @rdname GetAssaySummary
#' @export
setMethod(
  f = "GetAssaySummary",
  signature = "Seurat",
  definition = function(object,
                        group_name,
                        assay = NULL,
                        verbose = TRUE,
                        ...) {
    if (is.null(assay)) {
      assay <- object@active.assay
    }
    smry <- Seurat::Misc(
      object[[assay]]
    )$summary[[group_name]]
    if (is.null(smry)) {
      thisutils::log_message(
        "summary of '", group_name, "' does not yet exist",
        verbose = verbose,
        message_type = "warning"
      )
      thisutils::log_message(
        "Summarizing information for '", group_name, "'",
        verbose = verbose
      )
      object <- aggregate_assay(
        object,
        group_name = group_name,
        assay = assay
      )
      smry <- GetAssaySummary(
        object,
        group_name = group_name,
        assay = assay,
        verbose = verbose
      )
    }
    return(smry)
  }
)

#' @rdname GetAssaySummary
#' @export
setMethod(
  f = "GetAssaySummary",
  signature = "CSNObject",
  definition = function(object,
                        group_name,
                        assay = NULL,
                        verbose = TRUE,
                        ...) {
    return(
      GetAssaySummary(
        object@data,
        group_name = group_name,
        assay = assay,
        verbose = verbose
      )
    )
  }
)

#' @title Get fitted coefficients
#' @param object Network object
#' @param ... Other parameters
#' @method coef Network
#' @return Return the fitted coefficients
#' @export
coef.Network <- function(object, ...) {
  return(object@coefficients)
}

#' @title Get fitted coefficients
#' @param object CSNObject object
#' @param network Name of the network to use.
#' @param celltypes Celltypes to plot.
#' @param ... Other parameters
#' @method coef CSNObject
#' @return Return the fitted coefficients
#' @export
coef.CSNObject <- function(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  ...
) {
  networks <- GetNetwork(
    object,
    network = network,
    celltypes = celltypes
  )
  if (is.null(celltypes)) {
    return(
      lapply(
        networks, function(net) {
          net@coefficients
        }
      )
    )
  }
  return(networks@coefficients)
}

#' @param network network
#' @param celltypes cell types to analyze, NULL for all cell types
#' @rdname metrics
#' @export
setMethod(
  f = "metrics",
  signature = "Seurat",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    celltypes_all <- get_attribute(
      object,
      attribute = "celltypes"
    )
    celltypes <- intersect(
      celltypes %ss% celltypes_all,
      celltypes_all
    )
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )
    if (is.null(networks)) {
      return(NULL)
    }

    res <- lapply(
      networks,
      function(net) {
        net@metrics
      }
    )
    names(res) <- celltypes

    res
  }
)

setMethod(
  f = "metrics",
  signature = "CSNObject",
  definition = function(object,
                        network = DefaultNetwork(object),
                        celltypes = NULL,
                        ...) {
    celltypes_all <- get_attribute(
      object,
      attribute = "celltypes"
    )
    celltypes <- intersect(
      celltypes %ss% celltypes_all,
      celltypes_all
    )
    networks <- GetNetwork(
      object,
      network = network,
      celltypes = celltypes
    )

    res <- lapply(
      networks,
      function(net) {
        net@metrics
      }
    )
    names(res) <- celltypes

    return(res)
  }
)

#' @rdname metrics
#' @export
setMethod(
  f = "metrics",
  signature = "Network",
  definition = function(object, celltypes = NULL, ...) {
    return(object@metrics)
  }
)

#' @title Print Network objects
#' @param x x
#' @param ... other parameters
#' @rdname print
#' @export
#' @method print Network
print.Network <- function(x, ...) {
  coeffs <- methods::slot(x, "coefficients")
  modules_obj <- methods::slot(x, "modules")
  modules_meta <- methods::slot(modules_obj, "meta")

  if (nrow(modules_meta) == 0) {
    n_genes <- length(unique(coeffs$target))
    n_tfs <- length(unique(coeffs$tf))
  } else {
    n_genes <- length(unique(modules_meta$target))
    n_tfs <- length(unique(modules_meta$tf))
  }
  cat(paste0(
    "A Network object\n", "with ", n_tfs, " TFs and ",
    n_genes, " targets"
  ))
}

setMethod(
  "show",
  signature = "Network",
  function(object) {
    print(object)
  }
)

#' @title Print Modules objects
#' @rdname print
#' @export
#' @method print Modules
print.Modules <- function(x, ...) {
  n_mods <- length(x@features$genes_pos)
  cat(paste0(
    "An Modules object with ", n_mods, " TF modules"
  ))
}

setMethod(
  "show",
  signature = "Modules",
  function(object) {
    print(object)
  }
)

#' @title Print Regions objects
#' @rdname print
#' @export
#' @method print Regions
print.Regions <- function(x, ...) {
  n_regs <- length(x@ranges)
  n_peaks <- length(unique(x@peaks))
  cat(paste0(
    "An Regions object\n", "with ", n_regs, " candidate genomic regions ",
    "in ", n_peaks, " peaks"
  ))
}

setMethod(
  "show",
  signature = "Regions",
  function(object) {
    print(object)
  }
)

#' @export
DefaultAssay.CSNObject <- function(object, ...) {
  Seurat::DefaultAssay(object@data, ...)
}

#' @export
"DefaultAssay<-.CSNObject" <- function(object, ..., value) {
  Seurat::DefaultAssay(object@data, ...) <- value
  object
}

setMethod(
  f = "$",
  signature = "CSNObject",
  definition = function(x, name) {
    name <- as.character(name)
    if (identical(name, "celltypes")) {
      return(.csn_celltypes(x))
    }
    if (identical(name, "metadata")) {
      return(x@metadata)
    }
    md <- x@metadata
    if (!is.null(md) && !is.null(names(md)) && name %in% names(md)) {
      return(md[[name]])
    }
    if (!is.null(x@data) &&
      methods::is(x@data, "Seurat") &&
      !is.null(x@data@meta.data) &&
      name %in% colnames(x@data@meta.data)) {
      return(x@data@meta.data[[name]])
    }
    stop(
      sprintf("Unknown metadata field '%s' for CSNObject", name),
      call. = FALSE
    )
  }
)

#' @title Print CSNObject objects
#' @param x x
#' @param ... other parameters
#' @rdname print
#' @export
#' @method print CSNObject
print.CSNObject <- function(x, ...) {
  seurat <- x@data
  assays <- tryCatch(names(seurat@assays), error = function(e) character(0))
  active_assay <- tryCatch(seurat@active.assay, error = function(e) NA_character_)
  n_cells <- tryCatch(ncol(seurat), error = function(e) NA_integer_)
  n_features <- tryCatch(nrow(seurat[[active_assay]]), error = function(e) NA_integer_)
  celltypes <- .csn_celltypes(x)

  active_network <- DefaultNetwork(x)
  active_network <- active_network[1] %ss% ""
  if (is.na(active_network)) {
    active_network <- ""
  }
  network_names <- names(x@networks)
  has_motifs <- !is.null(x@regions@motifs)
  n_tfs <- tryCatch(ncol(x@regions@motifs2tfs), error = function(e) 0L)
  if (is.null(n_tfs) || length(n_tfs) == 0) {
    n_tfs <- 0L
  }

  cli::cli_h2("CSNObject")
  cli::cli_h3("Seurat")
  cli::cli_ul(c(
    "cells: {n_cells}",
    "features (active assay): {n_features} ({active_assay})",
    if (length(assays) > 0) {
      "assays: {assays}"
    } else {
      "assays: <none>"
    },
    "celltypes: {length(celltypes)}"
  ))

  cli::cli_h3("Regions")
  cli::cli_ul(c(
    "motifs: {ifelse(has_motifs, 'yes', 'no')}",
    "TF map columns: {n_tfs}"
  ))

  cli::cli_h3("Networks")
  cli::cli_ul(c(
    "active network: {if (nzchar(active_network)) active_network else '<none>'}",
    if (length(network_names) > 0) {
      "available networks: {network_names}"
    } else {
      "available networks: <none>"
    }
  ))

  if (nzchar(active_network) && active_network %in% names(x@networks)) {
    nets <- x@networks[[active_network]]
    if (is.list(nets) && length(nets) > 0) {
      ct_use <- intersect(names(nets), celltypes %ss% names(nets))
      if (length(ct_use) == 0) {
        ct_use <- names(nets)
      }
      cli::cli_h3("Network summary ({active_network})")
      ct_width <- max(nchar(ct_use))
      n_edges_vec <- integer(length(ct_use))
      out_reg_vec <- integer(length(ct_use))
      out_tar_vec <- integer(length(ct_use))
      is_empty <- logical(length(ct_use))
      for (i in seq_along(ct_use)) {
        ct <- ct_use[i]
        net <- nets[[ct]]
        if (is.null(net) || !methods::is(net, "Network")) {
          is_empty[i] <- TRUE
          next
        }
        is_empty[i] <- FALSE
        network_tbl <- tryCatch(export_csn(net), error = function(e) NULL)
        if (!is.null(network_tbl) && nrow(network_tbl) > 0) {
          n_edges_vec[i] <- nrow(network_tbl)
          out_reg_vec[i] <- length(unique(network_tbl$regulator))
          out_tar_vec[i] <- length(unique(network_tbl$target))
        } else {
          n_edges_vec[i] <- 0L
          out_reg_vec[i] <- 0L
          out_tar_vec[i] <- 0L
        }
      }
      w_edges <- max(1L, nchar(as.character(max(n_edges_vec))))
      w_out_reg <- max(1L, nchar(as.character(max(out_reg_vec))))
      w_out_tar <- max(1L, nchar(as.character(max(out_tar_vec))))
      items <- character(length(ct_use))
      for (i in seq_along(ct_use)) {
        ct <- ct_use[i]
        if (is_empty[i]) {
          items[i] <- sprintf("%-*s: %s", ct_width, ct, "<empty>")
        } else {
          items[i] <- sprintf(
            "%-*s: %*d edges (%*d regulators -> %*d targets)",
            ct_width,
            ct,
            w_edges, n_edges_vec[i],
            w_out_reg,
            out_reg_vec[i],
            w_out_tar,
            out_tar_vec[i]
          )
        }
      }
      if (length(items) > 0) {
        cli::cli_ul(items)
      }
    }
  }
}

setMethod(
  "show",
  signature = "CSNObject",
  function(object) {
    print(object)
  }
)

.process_Network <- function(
  object,
  r_squared_threshold = 0
) {
  metrics <- methods::slot(object, "metrics")
  metrics <- metrics[metrics$r_squared >= r_squared_threshold, ]
  targets <- unique(metrics$target)
  coefficients <- methods::slot(object, "coefficients")
  coefficients <- coefficients[coefficients$target %in% targets, ]

  if (is.null(coefficients) || nrow(coefficients) == 0) {
    methods::slot(object, "network") <- data.frame(
      regulator = character(),
      target = character(),
      weight = numeric(),
      mean_weight = numeric(),
      mean_corr = numeric(),
      n_regions = integer()
    )
    return(object)
  }

  coefficients_renamed <- dplyr::rename(
    coefficients,
    regulator = tf
  )
  aggregated_edges <- dplyr::group_by(
    coefficients_renamed,
    regulator,
    target
  )
  aggregated_edges <- dplyr::summarise(
    aggregated_edges,
    sum_weight = sum(coefficient),
    mean_weight = mean(coefficient),
    mean_corr = mean(corr),
    n_regions = dplyr::n(),
    .groups = "drop"
  )

  aggregated_edges <- dplyr::group_by(aggregated_edges, target)
  aggregated_edges <- dplyr::mutate(
    aggregated_edges,
    weight = thisutils::normalization(sum_weight, method = "unit_vector")
  )
  aggregated_edges <- dplyr::ungroup(aggregated_edges)
  aggregated_edges <- as.data.frame(aggregated_edges)

  aggregated_edges <- dplyr::select(
    aggregated_edges,
    regulator,
    target,
    weight,
    sum_weight,
    mean_weight,
    mean_corr,
    n_regions
  )

  methods::slot(object, "network") <- aggregated_edges
  return(object)
}

.process_csn <- function(
  object,
  r_squared_threshold = 0
) {
  active_network <- DefaultNetwork(object)
  networks <- GetNetwork(object, network = active_network)

  if (is.null(networks) || length(networks) == 0) {
    return(object)
  }

  for (celltype in names(networks)) {
    network <- networks[[celltype]]
    if (!is.null(network) && methods::is(network, "Network")) {
      network <- .process_Network(
        network,
        r_squared_threshold = r_squared_threshold
      )
      if (methods::is(object, "Seurat")) {
        object <- .multicsn_set_network_entry(object, active_network, celltype, network)
      } else {
        object@networks[[active_network]][[celltype]] <- network
      }
    }
  }

  return(object)
}

#' @param active_network Character string specifying which network to export
#' @param celltypes Character vector of cell types to export
#' @param weight_cutoff Numeric threshold for filtering edges by absolute weight
#'
#' @examples
#' \dontrun{
#' data(pbmcmultiome_sub, package = "scop")
#' data(motifs)
#' data(motif2tf)
#'
#' object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
#' object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
#' object <- initiate_object(
#'   object,
#'   group.by = "CellType",
#'   rna_assay = "RNA",
#'   peak_assay = "peaks",
#'   verbose = FALSE
#' )
#' genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
#' object <- find_motifs(
#'   object,
#'   pfm = motifs,
#'   motif_tfs = motif2tf,
#'   genome = genome,
#'   backend = "motifmatchr",
#'   verbose = FALSE
#' )
#' object <- inferCSN(object, cores = 2, verbose = FALSE)
#'
#' networks <- export_csn(object)
#' head(networks[[1]])
#' }
#' @rdname export_csn
#' @export
setMethod(
  "export_csn",
  signature = "Seurat",
  function(object,
           active_network = NULL,
           celltypes = NULL,
           weight_cutoff = NULL,
           ...) {
    active_network <- active_network %ss% DefaultNetwork(object)
    net_all <- GetNetwork(object, network = active_network)
    if (is.null(net_all)) {
      return(NULL)
    }

    celltypes_all <- names(net_all)
    celltypes <- intersect(celltypes %ss% celltypes_all, celltypes_all)

    res <- lapply(
      net_all[celltypes],
      function(x) {
        export_csn(x, weight_cutoff = weight_cutoff)
      }
    ) |>
      purrr::set_names(celltypes)

    if (length(celltypes) == 1) {
      return(res[[1]])
    }

    res
  }
)

setMethod(
  "export_csn",
  signature = "CSNObject",
  function(object,
           active_network = NULL,
           celltypes = NULL,
           weight_cutoff = NULL,
           ...) {
    active_network <- active_network %ss% DefaultNetwork(object)
    net_all <- object@networks[[active_network]]

    celltypes_all <- names(net_all)
    celltypes <- intersect(celltypes %ss% celltypes_all, celltypes_all)

    res <- lapply(
      net_all[celltypes],
      function(x) {
        export_csn(x, weight_cutoff = weight_cutoff)
      }
    ) |>
      purrr::set_names(celltypes)

    if (length(celltypes) == 1) {
      return(res[[1]])
    }

    return(res)
  }
)

#' @rdname export_csn
#' @export
setMethod(
  "export_csn",
  signature = "Network",
  function(object,
           weight_cutoff = NULL,
           ...) {
    network <- methods::slot(object, "network")

    if (is.null(network)) {
      return(NULL)
    }

    network <- inferCSN::network_format(network, abs_weight = FALSE)

    if (!is.null(weight_cutoff)) {
      network <- network[abs(network$weight) >= weight_cutoff, ]
    }

    return(network)
  }
)
