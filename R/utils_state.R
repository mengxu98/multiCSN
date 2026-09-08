#' @include setClass.R
NULL

.MULTICSN_STATE_SLOT <- "CSNObject"
.MULTICSN_SCHEMA_VERSION <- "2.0.0"

.empty_regions <- function() {
  methods::new(
    "Regions",
    motifs = NULL,
    motifs2tfs = NULL,
    ranges = GenomicRanges::GRanges(),
    peaks = numeric(0)
  )
}

.multicsn_state_template <- function() {
  list(
    schema_version = .MULTICSN_SCHEMA_VERSION,
    params = list(),
    active_network = character(0),
    networks = list(),
    regions = .empty_regions(),
    attributes = list(),
    tfs = NULL,
    summary = NULL,
    shortest_paths = NULL,
    edge_uniqueness = NULL,
    extract_genes_transition = NULL
  )
}

.normalize_multicsn_state <- function(state) {
  template <- .multicsn_state_template()
  if (is.null(state)) {
    return(template)
  }

  if (!is.list(state)) {
    stop("CSNObject state must be a list.", call. = FALSE)
  }

  for (nm in names(template)) {
    if (is.null(state[[nm]])) {
      state[[nm]] <- template[[nm]]
    }
  }

  state$schema_version <- state$schema_version %ss% .MULTICSN_SCHEMA_VERSION
  state$params <- state$params %ss% list()
  state$active_network <- state$active_network %ss% character(0)
  state$networks <- state$networks %ss% list()
  state$attributes <- state$attributes %ss% list()
  state$tfs <- state$tfs %ss% NULL
  state$summary <- state$summary %ss% NULL

  if (is.null(state$regions) || !methods::is(state$regions, "Regions")) {
    state$regions <- .empty_regions()
  }

  state
}

.get_multicsn_seurat <- function(object, require_state = FALSE) {
  if (methods::is(object, "Seurat")) {
    if (require_state) {
      state <- tryCatch(
        Seurat::Misc(object, .MULTICSN_STATE_SLOT),
        error = function(e) NULL
      )
      if (is.null(state)) {
        stop(
          "This Seurat object has not been initialized for CSNObject state. ",
          "Please run `initiate_object()` first.",
          call. = FALSE
        )
      }
    }
    return(object)
  }
  stop("Expected a Seurat object.", call. = FALSE)
}

.read_multicsn_state <- function(object, init = FALSE) {
  if (methods::is(object, "Seurat")) {
    state <- tryCatch(
      Seurat::Misc(object, .MULTICSN_STATE_SLOT),
      error = function(e) NULL
    )
    if (is.null(state) && init) {
      state <- .multicsn_state_template()
    }
    return(.normalize_multicsn_state(state))
  }
  stop("Expected a Seurat object.", call. = FALSE)
}

.write_multicsn_state <- function(object, state) {
  if (!methods::is(object, "Seurat")) {
    stop("CSNObject state can only be written to Seurat objects.", call. = FALSE)
  }

  normalized_state <- .normalize_multicsn_state(state)
  withCallingHandlers(
    {
      Seurat::Misc(object, .MULTICSN_STATE_SLOT) <- normalized_state
    },
    warning = function(w) {
      if (grepl("Overwriting miscellanous data for CSNObject", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  object
}

.multicsn_get <- function(object, field, default = NULL, init = FALSE) {
  state <- .read_multicsn_state(object, init = init)
  state[[field]] %ss% default
}

.multicsn_set <- function(object, field, value) {
  state <- .read_multicsn_state(object, init = TRUE)
  state[[field]] <- value
  .write_multicsn_state(object, state)
}

.multicsn_get_active_network <- function(object) {
  active <- .multicsn_get(object, "active_network", default = character(0))
  active[1] %ss% ""
}

.multicsn_set_active_network <- function(object, value) {
  value <- value[1] %ss% ""
  if (is.na(value) || !nzchar(value)) {
    stop("DefaultNetwork<-: value must be a non-empty character string", call. = FALSE)
  }
  nets <- .multicsn_get(object, "networks", default = list())
  if (!value %in% names(nets)) {
    stop(
      sprintf("DefaultNetwork<-: network '%s' not found in CSNObject state", value),
      call. = FALSE
    )
  }
  .multicsn_set(object, "active_network", value)
}

.multicsn_get_networks <- function(object, network = NULL) {
  nets <- .multicsn_get(object, "networks", default = list())
  network <- network %ss% .multicsn_get_active_network(object)
  if (is.null(network) || !nzchar(network)) {
    return(NULL)
  }
  nets[[network]]
}

.multicsn_set_networks <- function(object, network, value) {
  network <- network[1] %ss% ""
  nets <- .multicsn_get(object, "networks", default = list(), init = TRUE)
  nets[[network]] <- value
  object <- .multicsn_set(object, "networks", nets)
  if (nzchar(network)) {
    object <- .multicsn_set(object, "active_network", network)
  }
  object
}

.multicsn_set_network_entry <- function(object, network, celltype, value) {
  nets <- .multicsn_get_networks(object, network = network)
  nets <- nets %ss% list()
  nets[[celltype]] <- value
  .multicsn_set_networks(object, network, nets)
}

.multicsn_get_attributes <- function(object) {
  .multicsn_get(object, "attributes", default = list())
}

.multicsn_set_attributes <- function(object, value) {
  .multicsn_set(object, "attributes", value)
}

.multicsn_get_regions <- function(object) {
  .multicsn_get(object, "regions", default = .empty_regions())
}

.multicsn_set_regions <- function(object, value) {
  .multicsn_set(object, "regions", value)
}

.multicsn_get_tfs <- function(object) {
  .multicsn_get(object, "tfs", default = NULL)
}


.multicsn_get_celltypes <- function(object, active_network = NULL) {
  nets_active <- .multicsn_get_networks(object, network = active_network)
  is_dynamic_network <- FALSE

  if (is.list(nets_active) && length(nets_active) > 0) {
    first_net <- nets_active[[1]]
    if (!is.null(first_net) &&
      methods::is(first_net, "Network") &&
      identical(tryCatch(first_net@params$network_type, error = function(e) NULL), "dynamic")) {
      is_dynamic_network <- TRUE
    }
  }

  if (is_dynamic_network && length(nets_active) > 0) {
    return(names(nets_active))
  }

  attrs <- .multicsn_get_attributes(object)
  attr_names <- names(attrs) %ss% character(0)
  if (length(attr_names) > 0) {
    return(attr_names)
  }

  seurat <- .get_multicsn_seurat(object, require_state = FALSE)
  idents <- tryCatch(Seurat::Idents(seurat), error = function(e) NULL)
  if (!is.null(idents)) {
    lvls <- levels(idents)
    if (!is.null(lvls) && length(lvls) > 0) {
      return(lvls)
    }
    uniq <- unique(as.character(idents))
    uniq <- uniq[!is.na(uniq) & nzchar(uniq)]
    if (length(uniq) > 0) {
      return(uniq)
    }
  }

  character(0)
}

.multicsn_palette_colors <- function(values, palette = "Chinese", fallback = NULL) {
  vals <- as.character(values)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  vals_unique <- unique(vals)
  if (length(vals_unique) == 0) {
    return(character(0))
  }

  if (requireNamespace("thisplot", quietly = TRUE)) {
    cols <- tryCatch(
      thisplot::palette_colors(vals_unique, palette = palette),
      error = function(e) NULL
    )
    if (!is.null(cols) && length(cols) > 0) {
      if (is.null(names(cols)) || any(!nzchar(names(cols)))) {
        names(cols) <- vals_unique[seq_len(min(length(cols), length(vals_unique)))]
      }
      return(cols)
    }
  }

  if (is.null(fallback)) {
    fallback <- c(
      "#0066ff", "#0099ff", "#66cc33", "#66ff33", "#ccff00", "#ffff99",
      "#ffff33", "#ffcc33", "#ff9933", "#ff6633", "#cc3300", "#ff3333"
    )
  }
  cols <- rep(fallback, length.out = length(vals_unique))
  stats::setNames(cols, vals_unique)
}

.multicsn_network_type <- function(nets) {
  if (!is.list(nets) || length(nets) == 0) {
    return("unknown")
  }
  first_net <- nets[[1]]
  if (!is.null(first_net) &&
    methods::is(first_net, "Network") &&
    identical(tryCatch(first_net@params$network_type, error = function(e) NULL), "dynamic")) {
    return("dynamic")
  }
  "static"
}

.multicsn_resolve_network <- function(
  object,
  network = NULL,
  celltypes = NULL,
  preferred = c("active", "static", "dynamic", "any"),
  verbose = TRUE,
  caller = "function"
) {
  preferred <- match.arg(preferred)
  nets_all <- .multicsn_get(object, "networks", default = list())
  net_names <- names(nets_all) %ss% character(0)

  if (!is.null(network) && nzchar(network[1] %ss% "")) {
    return(network[1])
  }

  if (length(net_names) == 0) {
    stop("No networks found in Seurat::Misc(object, 'CSNObject').", call. = FALSE)
  }

  if (!is.null(celltypes) && length(celltypes) > 0) {
    candidates <- net_names[vapply(
      net_names,
      function(nm) all(celltypes %in% (names(nets_all[[nm]]) %ss% character(0))),
      logical(1)
    )]
    if (length(candidates) == 1) {
      thisutils::log_message(
        sprintf("%s: using network '%s' inferred from celltypes/states.", caller, candidates[[1]]),
        verbose = verbose
      )
      return(candidates[[1]])
    }
  }

  typed <- split(
    net_names,
    vapply(net_names, function(nm) .multicsn_network_type(nets_all[[nm]]), character(1))
  )
  active <- .multicsn_get_active_network(object)

  choose_first <- function(x) x[[1]] %ss% NULL

  chosen <- switch(preferred,
    active = if (nzchar(active)) active else choose_first(net_names),
    static = choose_first(typed$static) %ss% if (nzchar(active)) active else choose_first(net_names),
    dynamic = choose_first(typed$dynamic) %ss% if (nzchar(active)) active else choose_first(net_names),
    any = if (nzchar(active)) active else choose_first(net_names)
  )

  if (is.null(chosen) || !nzchar(chosen)) {
    stop("Unable to resolve a network automatically.", call. = FALSE)
  }

  thisutils::log_message(
    sprintf("%s: using network '%s'.", caller, chosen),
    verbose = verbose
  )
  chosen
}

.stop_csnobject_runtime <- function() {
  stop(
    "CSNObject runtime support has been removed. ",
    "Please use a Seurat object with state stored in Seurat::Misc(object, 'CSNObject').",
    call. = FALSE
  )
}

#' Show a compact summary of CSNObject state stored in a Seurat object
#'
#' @param object A \code{Seurat} object initialized by \code{multiCSN}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{object}.
#' @export
CSNObject <- function(object, ...) {
  object <- .get_multicsn_seurat(object, require_state = TRUE)
  state <- .read_multicsn_state(object)
  assays <- tryCatch(names(object@assays), error = function(e) character(0))
  active_assay <- tryCatch(object@active.assay, error = function(e) NA_character_)
  n_cells <- tryCatch(ncol(object), error = function(e) NA_integer_)
  n_features <- tryCatch(nrow(object[[active_assay]]), error = function(e) NA_integer_)
  celltypes <- .multicsn_get_celltypes(object)

  active_network <- .multicsn_get_active_network(object)
  network_names <- names(state$networks) %ss% character(0)
  has_motifs <- !is.null(state$regions@motifs)
  n_tfs <- tryCatch(ncol(state$regions@motifs2tfs), error = function(e) 0L)
  if (is.null(n_tfs) || length(n_tfs) == 0) {
    n_tfs <- 0L
  }

  cli::cli_h2("CSNObject")
  cli::cli_h3("Seurat")
  cli::cli_ul(c(
    "cells: {n_cells}",
    "features (active assay): {n_features} ({active_assay})",
    if (length(assays) > 0) "assays: {assays}" else "assays: <none>",
    "celltypes/states: {length(celltypes)}"
  ))

  cli::cli_h3("Regions")
  cli::cli_ul(c(
    "motifs: {ifelse(has_motifs, 'yes', 'no')}",
    "TF map columns: {n_tfs}"
  ))

  cli::cli_h3("Networks")
  cli::cli_ul(c(
    "active network: {if (nzchar(active_network)) active_network else '<none>'}",
    if (length(network_names) > 0) "available networks: {network_names}" else "available networks: <none>"
  ))

  if (nzchar(active_network) && active_network %in% network_names) {
    nets <- state$networks[[active_network]]
    if (is.list(nets) && length(nets) > 0) {
      ct_use <- names(nets)
      items <- vapply(
        ct_use,
        function(ct) {
          net <- nets[[ct]]
          if (is.null(net) || !methods::is(net, "Network")) {
            return(sprintf("%s: <empty>", ct))
          }
          network_tbl <- tryCatch(export_csn(net), error = function(e) NULL)
          if (is.null(network_tbl) || nrow(network_tbl) == 0) {
            return(sprintf("%s: 0 edges", ct))
          }
          sprintf(
            "%s: %d edges (%d regulators -> %d targets)",
            ct,
            nrow(network_tbl),
            length(unique(network_tbl$regulator)),
            length(unique(network_tbl$target))
          )
        },
        character(1)
      )
      cli::cli_h3("Network summary ({active_network})")
      cli::cli_ul(items)
    }
  }

  invisible(object)
}
