#' @include utils_sparse_blocks.R setClass.R
NULL

# ---- empty tables ----------------------------------------------------------

.empty_tf_gene_edges <- function() {
  data.frame(
    regulator = character(), target = character(),
    standardized_beta = numeric(), deletion_delta_bic = numeric(),
    stringsAsFactors = FALSE
  )
}

.empty_tf_region_edges <- function() {
  data.frame(
    regulator = character(), region = character(),
    standardized_beta = numeric(), deletion_delta_bic = numeric(),
    stringsAsFactors = FALSE
  )
}

.empty_region_gene_edges <- function() {
  data.frame(
    region = character(), target = character(),
    standardized_beta = numeric(), deletion_delta_bic = numeric(),
    stringsAsFactors = FALSE
  )
}

.empty_shared_accounting <- function() {
  data.frame(
    endpoint = character(), target = character(), n_candidates = integer(),
    n_selected = integer(), rss = numeric(), bic = numeric(), status = character(),
    stringsAsFactors = FALSE
  )
}

.empty_region_gene_accounting <- function() {
  data.frame(
    endpoint = character(), target = character(), n_candidates = integer(),
    n_selected = integer(), status = character(), seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

.layered_rbind <- function(parts, empty) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  parts <- parts[vapply(parts, nrow, integer(1)) > 0L]
  if (!length(parts)) {
    return(empty())
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# ---- numerics --------------------------------------------------------------

.layered_standardize_design <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  means <- colMeans(x)
  centered <- sweep(x, 2L, means, "-")
  scales <- sqrt(colSums(centered^2) / (nrow(centered) - 1L))
  variable <- is.finite(scales) & scales > 0
  list(
    matrix = sweep(centered[, variable, drop = FALSE], 2L, scales[variable], "/"),
    variable = variable, means = means, scales = scales
  )
}

.layered_response_statistics <- function(response, n_obs) {
  means <- as.numeric(Matrix::colSums(response)) / n_obs
  sumsq <- as.numeric(Matrix::colSums(response^2))
  centered_ss <- pmax(0, sumsq - n_obs * means^2)
  scales <- sqrt(centered_ss / (n_obs - 1L))
  list(means = means, scales = scales, variable = is.finite(scales) & scales > 0)
}

.layered_affine_alias <- function(standardized_crossproduct, n_obs, tolerance = 1e-10) {
  correlation <- abs(as.numeric(standardized_crossproduct)) / (n_obs - 1)
  is.finite(correlation) & correlation >= 1 - tolerance
}

.layered_lapply <- function(x, fun, cores = 1L) {
  cores <- as.integer(cores)
  if (cores > 1L && .Platform$OS.type != "windows") {
    return(parallel::mclapply(x, fun, mc.cores = cores, mc.preschedule = FALSE))
  }
  lapply(x, fun)
}

.layered_checkpoint <- function(root, name, ...) {
  if (is.null(root)) {
    return(NULL)
  }
  file.path(root, name, ...)
}

# ---- candidate builders ----------------------------------------------------

.layered_candidate_functions <- function(gene_by_cell, peaks2gene, peak_tf_gate,
                                         features, regulators, min_detected,
                                         endpoint = c("tf_gene", "tf_region")) {
  endpoint <- match.arg(endpoint)
  peaks2gene_rows <- .row_compressed_matrix(peaks2gene)
  gate_rows <- .row_compressed_matrix(peak_tf_gate)
  if (identical(endpoint, "tf_region")) {
    return(function(response_index) {
      .row_compressed_positive(gate_rows, response_index)
    })
  }
  gate_row_of_region <- match(colnames(peaks2gene), rownames(peak_tf_gate))
  gene_positions <- match(features, rownames(peaks2gene))
  self_of <- match(features, regulators)
  target_detected <- Matrix::rowSums(gene_by_cell[features, , drop = FALSE] != 0)
  tf_detected <- if (min_detected > 0L) {
    Matrix::rowSums(gene_by_cell[regulators, , drop = FALSE] != 0)
  } else {
    NULL
  }
  function(target_index) {
    if (min_detected > 0L && target_detected[[target_index]] < min_detected) {
      return(integer())
    }
    gene_row <- gene_positions[[target_index]]
    candidates <- integer()
    if (!is.na(gene_row)) {
      rows <- gate_row_of_region[.row_compressed_nonzero(peaks2gene_rows, gene_row)]
      rows <- rows[!is.na(rows)]
      if (length(rows)) {
        candidate_sets <- lapply(rows, function(row) .row_compressed_positive(gate_rows, row))
        candidates <- sort(unique(unlist(candidate_sets, use.names = FALSE)))
      }
    }
    self <- self_of[[target_index]]
    if (!is.na(self)) {
      candidates <- setdiff(candidates, self)
    }
    if (min_detected > 0L) {
      candidates <- candidates[tf_detected[candidates] >= min_detected]
    }
    as.integer(candidates)
  }
}

# ---- layer fitters ---------------------------------------------------------

.fit_layered_shared_design <- function(design, response, response_names,
                                       candidate_function, endpoint, settings,
                                       checkpoint_root = NULL) {
  n_obs <- nrow(design)
  gram <- crossprod(design)
  response_stats <- .layered_response_statistics(response, n_obs)
  chunk_size <- as.integer(settings$response_chunk)
  chunk_starts <- seq.int(1L, length(response_names), by = chunk_size)
  part_root <- gsub("-", "_", tolower(endpoint))
  chunk_results <- .layered_lapply(seq_along(chunk_starts), function(chunk_index) {
    first <- chunk_starts[[chunk_index]]
    last <- min(length(response_names), first + chunk_size - 1L)
    indices <- first:last
    cache_path <- .layered_checkpoint(
      checkpoint_root, paste0(part_root, "_parts"),
      sprintf("chunk_%06d_%06d.rds", first, last)
    )
    if (!is.null(cache_path) && file.exists(cache_path)) {
      return(readRDS(cache_path))
    }
    raw_xty <- as.matrix(Matrix::crossprod(
      design, response[, indices, drop = FALSE]
    ))
    xty <- raw_xty
    variable <- response_stats$variable[indices]
    if (any(variable)) {
      xty[, variable] <- sweep(
        raw_xty[, variable, drop = FALSE], 2L,
        response_stats$scales[indices][variable], "/"
      )
    }
    if (any(!variable)) xty[, !variable] <- 0
    candidates <- lapply(indices, candidate_function)
    if (isTRUE(settings$exclude_response_alias)) {
      candidates <- lapply(seq_along(candidates), function(j) {
        cand <- candidates[[j]]
        cand[!.layered_affine_alias(xty[cand, j], n_obs)]
      })
    }
    fit <- inferCSN::fit_greedy_l0_batch(
      gram = gram, xty = xty,
      response_ss = ifelse(variable, n_obs - 1, 0),
      candidates = candidates, n_obs = n_obs,
      max_support_size = settings$max_support_size,
      min_improvement = settings$min_improvement
    )
    edge_table <- if (length(fit$predictor_index)) {
      data.frame(
        feature = colnames(design)[fit$predictor_index],
        target = response_names[indices][fit$target_index],
        standardized_beta = as.numeric(fit$standardized_beta),
        deletion_delta_bic = as.numeric(fit$deletion_delta_bic),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        feature = character(), target = character(),
        standardized_beta = numeric(), deletion_delta_bic = numeric(),
        stringsAsFactors = FALSE
      )
    }
    accounting <- data.frame(
      endpoint = endpoint, target = response_names[indices],
      n_candidates = as.integer(fit$candidate_size),
      n_selected = as.integer(fit$support_size),
      rss = as.numeric(fit$rss), bic = as.numeric(fit$bic),
      status = as.character(fit$status), stringsAsFactors = FALSE
    )
    result <- list(edges = edge_table, accounting = accounting)
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, cache_path)
    }
    result
  }, cores = settings$cores)
  failures <- vapply(
    chunk_results, function(x) inherits(x, "try-error") || is.null(x$edges), logical(1)
  )
  if (any(failures)) {
    stop("Layered fit failed for ", sum(failures), " ", endpoint, " chunk(s).", call. = FALSE)
  }
  edges <- .layered_rbind(
    lapply(chunk_results, `[[`, "edges"),
    function() .empty_tf_gene_edges()
  )
  accounting <- .layered_rbind(
    lapply(chunk_results, `[[`, "accounting"),
    .empty_shared_accounting
  )
  list(edges = edges, accounting = accounting)
}

.fit_layered_region_gene <- function(peak_by_cell, gene_by_cell, features, peaks2gene,
                                     settings, checkpoint_root = NULL) {
  n_obs <- ncol(gene_by_cell)
  peak_rows <- .row_compressed_matrix(peak_by_cell)
  domain_rows <- .row_compressed_matrix(peaks2gene)
  gene_rows <- .row_compressed_matrix(gene_by_cell)
  domain_to_peak_row <- match(colnames(peaks2gene), rownames(peak_by_cell))
  genes_to_domain_row <- match(features, rownames(peaks2gene))
  genes_to_gene_row <- match(features, rownames(gene_by_cell))
  results <- .layered_lapply(seq_along(features), function(target_index) {
    target <- features[[target_index]]
    cache_path <- .layered_checkpoint(
      checkpoint_root, "region_gene_parts", sprintf("target_%06d.rds", target_index)
    )
    if (!is.null(cache_path) && file.exists(cache_path)) {
      return(readRDS(cache_path))
    }
    target_started <- proc.time()[["elapsed"]]
    domain_row <- genes_to_domain_row[[target_index]]
    candidate_rows <- if (is.na(domain_row)) {
      integer()
    } else {
      rows <- domain_to_peak_row[.row_compressed_nonzero(domain_rows, domain_row)]
      rows[!is.na(rows)]
    }
    candidates <- rownames(peak_by_cell)[candidate_rows]
    status <- "no_domain_peak"
    selected <- 0L
    part <- data.frame(
      region = character(), target = character(),
      standardized_beta = numeric(), deletion_delta_bic = numeric(),
      stringsAsFactors = FALSE
    )
    if (length(candidates)) {
      candidate_block <- .row_compressed_block(peak_rows, candidate_rows)
      means <- as.numeric(Matrix::colSums(candidate_block)) / n_obs
      centered_ss <- pmax(
        0, as.numeric(Matrix::colSums(candidate_block^2)) - n_obs * means^2
      )
      scales <- sqrt(centered_ss / (n_obs - 1L))
      variable <- is.finite(scales) & scales > 0
      candidates <- candidates[variable]
      candidate_rows <- candidate_rows[variable]
      if (length(candidates)) {
        candidate_block <- candidate_block[, variable, drop = FALSE]
        means <- means[variable]
        scales <- scales[variable]
        raw_gram <- as.matrix(Matrix::crossprod(candidate_block))
        gram <- (raw_gram - n_obs * tcrossprod(means)) / tcrossprod(scales)
        gram <- (gram + t(gram)) / 2
        diag(gram) <- n_obs - 1
        y <- .row_compressed_dense(gene_rows, genes_to_gene_row[[target_index]], n_obs)
        y_mean <- mean(y)
        y_scale <- sqrt(sum((y - y_mean)^2) / (n_obs - 1L))
        if (is.finite(y_scale) && y_scale > 0) {
          raw_xty <- as.numeric(Matrix::crossprod(candidate_block, y))
          xty <- matrix(
            (raw_xty - n_obs * means * y_mean) / (scales * y_scale), ncol = 1L
          )
          fit <- inferCSN::fit_greedy_l0_batch(
            gram = gram, xty = xty, response_ss = n_obs - 1,
            candidates = list(seq_along(candidates)), n_obs = n_obs,
            max_support_size = settings$max_support_size,
            min_improvement = settings$min_improvement
          )
          status <- as.character(fit$status[[1]])
          selected <- length(fit$predictor_index)
          if (selected) {
            part <- data.frame(
              region = candidates[fit$predictor_index], target = target,
              standardized_beta = as.numeric(fit$standardized_beta),
              deletion_delta_bic = as.numeric(fit$deletion_delta_bic),
              stringsAsFactors = FALSE
            )
          }
        } else {
          status <- "constant_response"
        }
      } else {
        status <- "no_variable_domain_peak"
      }
    }
    result <- list(
      edges = part,
      accounting = data.frame(
        endpoint = "region-gene", target = target, n_candidates = length(candidates),
        n_selected = selected, status = status,
        seconds = proc.time()[["elapsed"]] - target_started,
        stringsAsFactors = FALSE
      )
    )
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, cache_path)
    }
    result
  }, cores = settings$cores)
  failures <- vapply(
    results, function(x) inherits(x, "try-error") || is.null(x$accounting), logical(1)
  )
  if (any(failures)) {
    stop("Layered region-gene fit failed for ", sum(failures), " target(s).", call. = FALSE)
  }
  list(
    edges = .layered_rbind(lapply(results, `[[`, "edges"), .empty_region_gene_edges),
    accounting = .layered_rbind(
      lapply(results, `[[`, "accounting"), .empty_region_gene_accounting
    )
  )
}

# ---- input preparation -----------------------------------------------------

.layered_inputs <- function(object, cell_group, cells, regulators, targets,
                            min_peak_cells, upstream, downstream, tf_region_scope,
                            renormalize, sort_regulators, verbose) {
  params <- Params(object)
  rna_assay <- params$rna_assay
  peak_assay <- params$peak_assay
  if (is.null(peak_assay)) {
    stop("A layered multi-omic fit requires a chromatin accessibility assay.", call. = FALSE)
  }
  catalogue <- if (is.null(cell_group)) NULL else {
    get_attribute(object, celltypes = cell_group, attribute = "cells")
  }
  selected <- if (is.null(cells)) catalogue else as.character(cells)
  if (is.null(selected) || !length(selected)) {
    stop("No cells selected for the layered fit.", call. = FALSE)
  }
  if (!is.null(catalogue)) {
    outside <- setdiff(selected, catalogue)
    if (length(outside)) {
      stop(
        "Selected cells are not part of the '", cell_group, "' catalogue: ",
        paste(utils::head(outside, 3L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (isTRUE(renormalize)) {
    thisutils::log_message("Re-normalizing RNA and chromatin for the selected cells", verbose = verbose)
    object <- subset(object, cells = selected)
    object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
    object <- Signac::RunTFIDF(object, assay = peak_assay, verbose = FALSE)
  }
  gene_layer <- get_layer_data(object, assay = rna_assay, layer = "data")
  peak_layer <- get_layer_data(object, assay = peak_assay, layer = "data")
  peaks_requested <- get_attribute(object, celltypes = cell_group, attribute = "peaks")
  if (is.null(peaks_requested)) {
    peaks_requested <- rownames(peak_layer)
  }
  targets_requested <- targets
  if (is.null(targets_requested)) {
    targets_requested <- get_attribute(object, celltypes = cell_group, attribute = "genes")
  }
  if (is.null(targets_requested)) {
    targets_requested <- rownames(gene_layer)
  }
  regulators_requested <- regulators
  if (is.null(regulators_requested)) {
    regulators_requested <- get_attribute(object, attribute = "tfs")
  }
  gene_by_cell <- gene_layer[, selected, drop = FALSE]
  peak_by_cell <- peak_layer[peaks_requested, selected, drop = FALSE]
  gene_annotation <- get_peak_annotation(object, peak_assay)
  features <- sort(unique(intersect(
    intersect(unique(gene_annotation$gene_name), targets_requested),
    rownames(gene_by_cell)
  )))
  regulators <- intersect(regulators_requested, rownames(gene_by_cell))
  if (isTRUE(sort_regulators)) {
    regulators <- sort(regulators)
  }
  if (!length(features) || !length(regulators)) {
    stop("The layered fit requires both target genes and transcription factors.", call. = FALSE)
  }
  regions <- NetworkRegions(object)
  peaks2motif <- regions@motifs@data[peaks_requested, , drop = FALSE]
  motif2tf <- regions@motifs2tfs[, regulators, drop = FALSE]
  peak_ranges <- parse_peak_ranges(peaks_requested)
  names(peak_ranges) <- peaks_requested
  thisutils::log_message(
    "Selecting candidate regulatory regions near ", length(features), " genes",
    verbose = verbose
  )
  peaks_near_gene <- find_peaks_near_genes(
    peaks = peak_ranges, method = "Signac",
    genes = gene_annotation[gene_annotation$gene_name %in% features, , drop = FALSE],
    upstream = upstream, downstream = downstream, only_tss = FALSE
  )
  peaks2gene <- collapse_peak_gene_domains(peaks_near_gene)
  observed <- Matrix::rowSums(peak_by_cell != 0) >= min_peak_cells
  domain_peak <- Matrix::colSums(peaks2gene != 0) > 0 &
    observed[colnames(peaks2gene)]
  peaks2gene <- peaks2gene[, domain_peak, drop = FALSE]
  all_peak_by_cell <- peak_by_cell
  peak_by_cell <- peak_by_cell[colnames(peaks2gene), , drop = FALSE]
  motif_peak <- Matrix::rowSums(
    peaks2motif[colnames(peaks2gene), , drop = FALSE] != 0
  ) > 0
  tf_region_names <- if (identical(tf_region_scope, "domain_peaks")) {
    colnames(peaks2gene)[motif_peak]
  } else {
    peaks_requested[observed & Matrix::rowSums(peaks2motif != 0) > 0]
  }
  motif_use <- intersect(colnames(peaks2motif), rownames(motif2tf))
  peak_tf_gate <- peaks2motif[tf_region_names, motif_use, drop = FALSE] %*%
    motif2tf[motif_use, regulators, drop = FALSE]
  peak_tf_gate[is.na(peak_tf_gate)] <- 0
  peak_tf_gate[peak_tf_gate < 0] <- 0
  peak_tf_gate <- methods::as(peak_tf_gate, "generalMatrix")
  if (!identical(tf_region_scope, "domain_peaks")) {
    keep <- Matrix::rowSums(peak_tf_gate != 0) > 0
    tf_region_names <- tf_region_names[keep]
    peak_tf_gate <- peak_tf_gate[tf_region_names, , drop = FALSE]
  }
  thisutils::log_message(
    "Layered fit input: ", length(selected), " cells, ", length(features),
    " target genes, ", length(regulators), " transcription factors, ",
    length(tf_region_names), " candidate regions",
    verbose = verbose
  )
  list(
    object = object, cells = selected, features = features, regulators = regulators,
    gene_by_cell = gene_by_cell, peak_by_cell = peak_by_cell,
    all_peak_by_cell = all_peak_by_cell, peaks2gene = peaks2gene,
    peak_tf_gate = peak_tf_gate, tf_region_names = tf_region_names
  )
}

# ---- driver ----------------------------------------------------------------

.fit_layered_core <- function(inputs, settings) {
  gene_by_cell <- inputs$gene_by_cell
  features <- inputs$features
  regulators <- inputs$regulators
  tf_design_info <- .layered_standardize_design(
    Matrix::t(gene_by_cell[regulators, , drop = FALSE])
  )
  tf_design <- tf_design_info$matrix
  regulators <- colnames(tf_design)
  peak_tf_gate <- inputs$peak_tf_gate[, regulators, drop = FALSE]
  tf_gene_candidate <- .layered_candidate_functions(
    gene_by_cell = gene_by_cell, peaks2gene = inputs$peaks2gene,
    peak_tf_gate = peak_tf_gate, features = features, regulators = regulators,
    min_detected = settings$min_detected, endpoint = "tf_gene"
  )
  tf_region_candidate <- .layered_candidate_functions(
    gene_by_cell = gene_by_cell, peaks2gene = inputs$peaks2gene,
    peak_tf_gate = peak_tf_gate, features = features, regulators = regulators,
    min_detected = settings$min_detected, endpoint = "tf_region"
  )
  thisutils::log_message("Fitting the TF-gene layer", verbose = settings$verbose)
  tf_gene <- .fit_layered_shared_design(
    tf_design, Matrix::t(gene_by_cell[features, , drop = FALSE]), features,
    tf_gene_candidate, "TF-gene", settings, settings$checkpoint_dir
  )
  names(tf_gene$edges)[names(tf_gene$edges) == "feature"] <- "regulator"
  thisutils::log_message("Fitting the TF-region layer", verbose = settings$verbose)
  tf_region <- .fit_layered_shared_design(
    tf_design,
    Matrix::t(inputs$all_peak_by_cell[inputs$tf_region_names, , drop = FALSE]),
    inputs$tf_region_names, tf_region_candidate, "TF-region", settings,
    settings$checkpoint_dir
  )
  names(tf_region$edges)[names(tf_region$edges) == "feature"] <- "regulator"
  names(tf_region$edges)[names(tf_region$edges) == "target"] <- "region"
  names(tf_region$accounting)[names(tf_region$accounting) == "target"] <- "region"
  thisutils::log_message("Fitting the region-gene layer", verbose = settings$verbose)
  region_gene <- .fit_layered_region_gene(
    inputs$peak_by_cell, gene_by_cell, features, inputs$peaks2gene, settings,
    settings$checkpoint_dir
  )
  tf_gene$edges$weight <- network_ordinal_weight(
    sign(tf_gene$edges$standardized_beta), tf_gene$edges$deletion_delta_bic
  )
  tf_region$edges$weight <- network_ordinal_weight(
    sign(tf_region$edges$standardized_beta), tf_region$edges$deletion_delta_bic
  )
  region_gene$edges$weight <- network_ordinal_weight(
    sign(region_gene$edges$standardized_beta), region_gene$edges$deletion_delta_bic
  )
  thisutils::log_message("Assembling factorized chains", verbose = settings$verbose)
  chains <- factorized_chain_support(
    tf_region$edges, region_gene$edges, tf_gene$edges
  )
  mediated <- strongest_mediated_projection(chains)
  list(
    tf_gene = tf_gene$edges, tf_region = tf_region$edges,
    region_gene = region_gene$edges, triplets = chains, mediated_tf_gene = mediated,
    accounting = list(
      tf_gene = tf_gene$accounting, tf_region = tf_region$accounting,
      region_gene = region_gene$accounting
    ),
    summary = data.frame(
      endpoint = c("TF-gene", "TF-region", "region-gene", "factorized-chain", "mediated-TF-gene"),
      selected = c(
        nrow(tf_gene$edges), nrow(tf_region$edges), nrow(region_gene$edges),
        nrow(chains), nrow(mediated)
      ),
      stringsAsFactors = FALSE
    )
  )
}

.layered_store <- function(object, fit, network_name, cell_group, settings) {
  entry <- function(coefficients, metrics, regulators, targets, endpoint) {
    parameters <- settings
    parameters$endpoint <- endpoint
    methods::new(
      "Network",
      params = parameters,
      regulators = as.character(regulators),
      targets = as.character(targets),
      metrics = as.data.frame(metrics),
      coefficients = as.data.frame(coefficients),
      network = as.data.frame(coefficients)
    )
  }
  layers <- list(
    tf_gene = list(
      edges = fit$tf_gene, accounting = fit$accounting$tf_gene,
      regulators = unique(fit$tf_gene$regulator), targets = unique(fit$tf_gene$target)
    ),
    tf_region = list(
      edges = fit$tf_region, accounting = fit$accounting$tf_region,
      regulators = unique(fit$tf_region$regulator), targets = unique(fit$tf_region$region)
    ),
    region_gene = list(
      edges = fit$region_gene, accounting = fit$accounting$region_gene,
      regulators = unique(fit$region_gene$region), targets = unique(fit$region_gene$target)
    )
  )
  for (layer in names(layers)) {
    value <- layers[[layer]]
    layer_name <- paste0(network_name, "_", layer)
    network <- entry(
      value$edges, value$accounting, value$regulators, value$targets, layer_name
    )
    object <- .multicsn_set_network_entry(
      object, network = layer_name, celltype = cell_group, value = network
    )
  }
  object
}

#' Fit a layered multi-omic gene regulatory network
#'
#' Fits the layered multiCSN model on single-cell multiome data: a
#' transcription-factor-to-gene layer, a transcription-factor-to-region layer and
#' a region-to-gene layer, each solved as independent per-target Greedy-l0
#' problems, followed by the factorized chain projection that links the three
#' layers. By default only the selected cells are used and the expression layers
#' stored in the object are taken as-is.
#'
#' @param object A Seurat object with chromatin accessibility and gene expression
#'   assays plus the multiCSN state (see \code{\link{initiate_object}}).
#' @param cell_group Catalogue key (cell type, state or timepoint) used to select
#'   cells, peaks and target genes.
#' @param cells Optional character vector of exact cells to fit; when supplied the
#'   fit uses exactly these cells and validates them against \code{cell_group}.
#' @param regulators,targets Optional character vectors restricting the
#'   transcription factors and target genes.
#' @param min_detected Minimum number of cells in which a target gene or a
#'   transcription factor must be detected; \code{0} disables the filter.
#' @param min_peak_cells Minimum number of cells in which a candidate regulatory
#'   region must be detected.
#' @param upstream,downstream Window in base pairs used to assign candidate
#'   regions to genes.
#' @param tf_region_scope Candidate regions of the TF-region layer:
#'   \code{"motif_peaks"} fits every motif-carrying accessible region, while
#'   \code{"domain_peaks"} restricts the layer to regions inside gene domains.
#' @param exclude_response_alias Drop candidates whose standardized cross-product
#'   with the response is an affine alias (|r| >= 1 - 1e-10).
#' @param renormalize Re-normalize RNA counts and re-run TF-IDF on the selected
#'   cells before fitting.
#' @param response_chunk Number of targets fitted per chunk.
#' @param max_support_size,max_improvement Greedy-l0 stopping rules forwarded to
#'   \code{inferCSN::fit_greedy_l0_batch()}.
#' @param sort_regulators Sort the transcription factors before fitting.
#' @param cores Number of forked workers used per chunk (ignored on Windows).
#' @param checkpoint_dir Optional directory storing per-chunk checkpoints so an
#'   interrupted fit can resume; \code{NULL} keeps everything in memory.
#' @param store Store the three layers as \code{Network} entries in the returned
#'   object.
#' @param network_name Name of the stored network entries.
#' @param verbose Print progress messages.
#' @return A list with \code{tf_gene}, \code{tf_region}, \code{region_gene},
#'   \code{triplets}, \code{mediated_tf_gene}, \code{accounting}, \code{summary},
#'   \code{settings}, \code{cells} and, when \code{store = TRUE}, \code{object}.
#' @export
fit_layered_network <- function(object,
                                cell_group = NULL,
                                cells = NULL,
                                regulators = NULL,
                                targets = NULL,
                                min_detected = 0L,
                                min_peak_cells = 3L,
                                upstream = 1e5,
                                downstream = 0L,
                                tf_region_scope = c("motif_peaks", "domain_peaks"),
                                exclude_response_alias = FALSE,
                                renormalize = FALSE,
                                response_chunk = 1024L,
                                max_support_size = NULL,
                                min_improvement = 1e-10,
                                sort_regulators = FALSE,
                                cores = 1L,
                                checkpoint_dir = NULL,
                                store = FALSE,
                                network_name = "multiCSN",
                                verbose = TRUE) {
  tf_region_scope <- match.arg(tf_region_scope)
  min_detected <- as.integer(min_detected)
  min_peak_cells <- as.integer(min_peak_cells)
  response_chunk <- as.integer(response_chunk)
  cores <- as.integer(cores)
  if (is.na(min_detected) || min_detected < 0L) {
    stop("`min_detected` must be a non-negative integer.", call. = FALSE)
  }
  if (is.na(min_peak_cells) || min_peak_cells < 1L) {
    stop("`min_peak_cells` must be a positive integer.", call. = FALSE)
  }
  if (is.na(response_chunk) || response_chunk < 1L) {
    stop("`response_chunk` must be a positive integer.", call. = FALSE)
  }
  if (is.na(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  settings <- list(
    cell_group = cell_group, min_detected = min_detected,
    min_peak_cells = min_peak_cells, upstream = upstream, downstream = downstream,
    tf_region_scope = tf_region_scope,
    exclude_response_alias = isTRUE(exclude_response_alias),
    renormalize = isTRUE(renormalize), response_chunk = response_chunk,
    max_support_size = max_support_size, min_improvement = min_improvement,
    sort_regulators = isTRUE(sort_regulators), cores = cores,
    checkpoint_dir = checkpoint_dir, verbose = isTRUE(verbose)
  )
  inputs <- .layered_inputs(
    object = object, cell_group = cell_group, cells = cells,
    regulators = regulators, targets = targets, min_peak_cells = min_peak_cells,
    upstream = upstream, downstream = downstream, tf_region_scope = tf_region_scope,
    renormalize = renormalize, sort_regulators = sort_regulators, verbose = verbose
  )
  fit <- .fit_layered_core(inputs, settings)
  result <- c(
    fit,
    list(settings = settings, cells = inputs$cells, targets = inputs$features,
         regulators = inputs$regulators)
  )
  if (isTRUE(store)) {
    result$object <- .layered_store(
      inputs$object, fit, network_name, cell_group, settings
    )
  }
  result
}
