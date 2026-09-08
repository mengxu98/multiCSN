#' @include setClass.R
#' @include setGenerics.R
#' @include utils_object.R
NULL

#' Internal helper to normalize regulators/targets specification
#' into a per-celltype named list.
#'
#' @param x Regulators/targets specification: a character vector, or a named list
#'   of character vectors with names matching \code{celltypes}.
#' @param celltypes Cell types the specification is expanded to.
#' @param arg_name Argument name used in error messages.
parse_regulators_targets <- function(x, celltypes, arg_name = "targets") {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.atomic(x) && !is.list(x)) {
    return(
      purrr::set_names(
        rep(list(as.character(x)), length(celltypes)),
        celltypes
      )
    )
  }

  if (is.list(x)) {
    if (length(x) == 1L) {
      return(
        purrr::set_names(
          rep(list(as.character(x[[1]])), length(celltypes)),
          celltypes
        )
      )
    }
    if (is.null(names(x)) ||
      length(x) != length(celltypes) ||
      !identical(sort(names(x)), sort(as.character(celltypes)))) {
      stop(
        sprintf(
          "When %s is a list with multiple elements, length and names must match celltypes.",
          arg_name
        )
      )
    }
    x[] <- lapply(x, as.character)
    return(x[as.character(celltypes)])
  }

  stop(sprintf("%s must be a vector or list.", arg_name))
}

#' @param object A Seurat object containing gene expression and/or chromatin accessibility data.
#' @param regulators Regulators to consider for CSN inference.
#' @param targets Targets to consider for CSN inference.
#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "matrix",
  definition = function(object,
                        regulators = character(0),
                        targets = character(0),
                        ...) {
    if (length(regulators) == 0) {
      regulators <- colnames(object)
    } else {
      regulators <- intersect(
        regulators, colnames(object)
      )
    }
    if (length(targets) == 0) {
      targets <- colnames(object)
    } else {
      targets <- intersect(
        targets, colnames(object)
      )
    }

    object <- methods::new(
      Class = "Network",
      data = object,
      regulators = regulators,
      targets = targets
    )

    return(object)
  }
)

#' @param object A Seurat object containing gene expression and/or chromatin accessibility data.
#' @param regulators Regulators to consider for CSN inference. Can be a character
#'   vector applied to all cell types, or a named list (names must match cell
#'   types) of character vectors per cell type.
#' @param targets Targets to consider for CSN inference. Same format rules as
#'   \code{regulators}.
#' @param celltype Selected celltype(s) to analyze. If NULL, all cells will be used.
#' @param group.by Column name in Seurat meta.data for cell type information.
#' @param filter_mode Filter mode for identifying cell-type specific features.
#' @param filter_by When filter_mode is "variable", specify "aggregate" or "celltype".
#' @param rna_assay Name of the gene expression assay.
#' @param rna_min_pct,rna_logfc_threshold,rna_test_method,n_variable_genes RNA params.
#' @param peak_assay Name of the chromatin accessibility assay. If NULL, only RNA.
#' @param peak_min_pct,peak_logfc_threshold,peak_test_method,n_variable_peaks Peak params.
#' @param regions Candidate regions (GRanges or data frame). If NULL, all peaks.
#' @param regions_extend Base pairs to extend regions.
#' @param exclude_exons Whether to consider exons for binding site inference.
#' @param only_pos Only return positive markers.
#' @param verbose Print progress messages.
#' @param p_value Significance threshold for filtering features.
#' @param ... Additional arguments passed to marker detection functions.
#' @return CSNObject object.
#'
#' @examples
#' data(pbmcmultiome_sub, package = "scop")
#'
#' object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
#' object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
#'
#' object <- initiate_object(
#'   object,
#'   group.by = "CellType",
#'   rna_assay = "RNA",
#'   peak_assay = "peaks",
#'   verbose = FALSE
#' )
#' object
#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "Seurat",
  definition = function(object,
                        celltype = NULL,
                        regulators = NULL,
                        targets = NULL,
                        group.by = NULL,
                        filter_mode = c(
                          "celltype_specific",
                          "variable",
                          "celltype",
                          "unfiltered"
                        ),
                        filter_by = c("celltype", "aggregate"),
                        rna_assay = "RNA",
                        rna_min_pct = 0.1,
                        rna_logfc_threshold = 0.25,
                        rna_test_method = "wilcox",
                        n_variable_genes = 2000,
                        peak_assay = NULL,
                        peak_min_pct = 0.05,
                        peak_logfc_threshold = 0.1,
                        peak_test_method = "LR",
                        n_variable_peaks = "q5",
                        regions = NULL,
                        regions_extend = 0,
                        exclude_exons = TRUE,
                        only_pos = TRUE,
                        verbose = TRUE,
                        p_value = 0.05,
                        ...) {
    if (!is.null(group.by) && group.by != "aggregate") {
      if (!group.by %in% colnames(object@meta.data)) {
        stop(sprintf("Column '%s' not found in object metadata", group.by))
      }
      Seurat::Idents(object) <- group.by
    } else {
      Seurat::Idents(object) <- "aggregate"
    }

    if (!is.null(celltype)) {
      if (!all(celltype %in% Seurat::Idents(object))) {
        stop("Some specified clusters not found in data")
      }
      object <- subset(object, idents = celltype)
      if (length(celltype) == 1) {
        Seurat::Idents(object) <- celltype
      }
    }

    filter_mode <- match.arg(filter_mode)
    filter_by <- match.arg(filter_by)
    celltypes <- unique(Seurat::Idents(object))

    Seurat::DefaultAssay(object) <- rna_assay

    use_custom_targets <- !is.null(targets)
    regulators_list <- NULL
    targets_list <- NULL

    all_genes <- rownames(
      Seurat::GetAssayData(
        object,
        assay = rna_assay,
        layer = "data"
      )
    )

    if (!is.null(regulators)) {
      regulators_list <- parse_regulators_targets(
        regulators,
        celltypes,
        "regulators"
      )
      regulators_list <- lapply(
        regulators_list,
        function(r) intersect(r, all_genes)
      )
    }

    if (use_custom_targets) {
      thisutils::log_message(
        "Using user-specified targets; skipping gene marker detection.",
        verbose = verbose
      )
      targets_list <- parse_regulators_targets(
        targets,
        celltypes,
        "targets"
      )
      targets_list <- lapply(
        targets_list,
        function(t) intersect(t, all_genes)
      )

      genes_markers <- purrr::map_dfr(
        celltypes,
        function(ct) {
          feats <- targets_list[[ct]]
          if (is.null(feats) || length(feats) == 0) {
            return(.create_empty_feature_df(is_peak = FALSE))
          }
          .create_feature_df(
            features = feats,
            celltype = ct,
            infinite_logfc = TRUE,
            is_peak = FALSE
          )
        }
      )
    } else {
      genes_markers <- process_features(
        object = object,
        celltypes = celltypes,
        mode = filter_mode,
        by = filter_by,
        is_peak = FALSE,
        group.by = group.by,
        params = list(
          assay = rna_assay,
          min_pct = rna_min_pct,
          logfc_threshold = rna_logfc_threshold,
          test_method = rna_test_method,
          n_features = n_variable_genes,
          only_pos = only_pos,
          p_value = p_value
        ),
        verbose = verbose,
        ...
      )
    }

    if (is.null(peak_assay)) {
      thisutils::log_message(
        "No peak assay found. Specify 'peak_assay' in 'initiate_object'.",
        verbose = verbose,
        message_type = "warning"
      )
    }
    peaks_markers <- NULL
    if (!is.null(peak_assay)) {
      peak_filter_mode <- filter_mode
      if (identical(filter_mode, "celltype_specific")) {
        peak_filter_mode <- "variable"
        thisutils::log_message(
          "Using variable peak filtering when `filter_mode = 'celltype_specific'`.",
          verbose = verbose
        )
      }
      Seurat::DefaultAssay(object) <- peak_assay
      object$nCount_peaks <- Matrix::colSums(
        Seurat::GetAssayData(object, assay = peak_assay, layer = "data") > 0
      )

      peaks_markers <- process_features(
        object = object,
        celltypes = celltypes,
        mode = peak_filter_mode,
        by = filter_by,
        is_peak = TRUE,
        params = list(
          assay = peak_assay,
          min_pct = peak_min_pct,
          logfc_threshold = peak_logfc_threshold,
          test_method = peak_test_method,
          n_features = n_variable_peaks,
          only_pos = only_pos,
          latent.vars = "nCount_peaks",
          p_value = p_value
        ),
        verbose = verbose,
        ...
      )
    }

    regions_obj <- process_regions(
      object = object,
      regions = regions,
      regions_extend = regions_extend,
      peak_assay = peak_assay,
      exclude_exons = exclude_exons,
      verbose = verbose,
      peaks_markers = peaks_markers,
      celltypes = celltypes,
      ...
    )

    attributes <- process_attributes(
      object = object,
      celltypes = celltypes,
      genes_markers = genes_markers,
      peaks_markers = peaks_markers,
      regions_obj = regions_obj,
      verbose = verbose,
      p_value = p_value
    )

    params <- list(
      peak_assay = peak_assay,
      rna_assay = rna_assay,
      exclude_exons = exclude_exons,
      filter_mode = filter_mode,
      filter_by = filter_by
    )
    if (!is.null(regulators_list)) {
      params$regulators <- regulators_list
    }

    summary_df <- create_summary(attributes, celltypes, peak_assay)
    thisutils::log_message("Summary of cell-type specific features:", verbose = verbose)
    thisutils::log_message(summary_df, verbose = verbose)

    object <- .write_multicsn_state(
      object,
      list(
        schema_version = .MULTICSN_SCHEMA_VERSION,
        params = params,
        active_network = character(0),
        networks = list(),
        regions = regions_obj,
        attributes = attributes,
        tfs = NULL,
        summary = create_summary(
          attributes,
          celltypes,
          peak_assay
        )
      )
    )

    return(object)
  }
)

#' @rdname initiate_object
#' @export
setMethod(
  f = "initiate_object",
  signature = "CSNObject",
  definition = function(object, ...) {
    .stop_csnobject_runtime()
  }
)

process_features <- function(
  object,
  celltypes,
  mode,
  by,
  is_peak,
  group.by = NULL,
  params,
  verbose,
  ...
) {
  switch(
    EXPR = mode,
    "celltype_specific" = process_celltype_specific_genes(
      object = object,
      celltypes = celltypes,
      group.by = group.by,
      params = params,
      verbose = verbose
    ),
    "celltype" = process_celltype_features(
      object,
      celltypes,
      is_peak,
      params,
      verbose,
      ...
    ),
    "variable" = process_variable_features(
      object,
      celltypes,
      by,
      is_peak,
      params,
      verbose,
      ...
    ),
    "unfiltered" = create_feature_matrix(
      rownames(Seurat::GetAssayData(object, assay = params$assay, layer = "data")),
      celltypes,
      is_peak = is_peak
    )
  )
}

process_celltype_features <- function(
  object,
  celltypes,
  is_peak,
  params,
  verbose,
  ...
) {
  if (length(celltypes) <= 1) {
    return(
      process_variable_features(
        object, celltypes, "aggregate", is_peak, params, verbose, ...
      )
    )
  }

  feature_type <- if (is_peak) "peaks" else "genes"
  thisutils::log_message(
    "Finding celltype-specific ", feature_type,
    verbose = verbose
  )

  purrr::map_dfr(
    celltypes, function(x) {
      markers <- suppressWarnings(
        Seurat::FindMarkers(
          object = object,
          ident.1 = x,
          assay = params$assay,
          min.pct = params$min_pct,
          logfc.threshold = params$logfc_threshold,
          test.use = params$test_method,
          only.pos = params$only_pos,
          verbose = FALSE,
          ...
        )
      )

      if (!is.null(markers) && nrow(markers) > 0) {
        markers$celltype <- x
        markers[[if (is_peak) "peak" else "gene"]] <- rownames(markers)
      }

      return(markers)
    }
  )
}

process_variable_features <- function(
  object,
  celltypes,
  by,
  is_peak,
  params,
  verbose,
  ...
) {
  if (by == "aggregate") {
    object <- find_variable_features(
      object,
      is_peak,
      params
    )
    var_features <- Seurat::VariableFeatures(object)
    res <- create_feature_matrix(
      var_features,
      celltypes,
      is_peak = is_peak
    )
    return(res)
  }

  res <- purrr::map_dfr(
    celltypes, function(x) {
      cell_names <- colnames(object)[Seurat::Idents(object) == x]
      cell_subset <- object[, cell_names]
      cell_subset <- find_variable_features(
        cell_subset,
        is_peak,
        params
      )
      var_features <- Seurat::VariableFeatures(cell_subset)
      create_feature_matrix(var_features, x, is_peak = is_peak)
    }
  )
  return(res)
}

process_celltype_specific_genes <- function(
  object,
  celltypes,
  group.by,
  params,
  verbose
) {
  if (is.null(group.by) || identical(group.by, "aggregate")) {
    stop(
      "`filter_mode = 'celltype_specific'` requires a valid `group.by` ",
      "column in Seurat metadata."
    )
  }
  gene_sets <- get_celltype_specific_genes(
    seurat = object,
    group.by = group.by,
    assay = params$assay,
    verbose = verbose
  )

  if (length(gene_sets) == 0) {
    return(.create_empty_feature_df(is_peak = FALSE))
  }

  res <- purrr::map_dfr(
    celltypes,
    function(x) {
      features <- gene_sets[[as.character(x)]]
      if (is.null(features) || length(features) == 0) {
        return(.create_empty_feature_df(is_peak = FALSE))
      }
      .create_feature_df(features, x, is_peak = FALSE)
    }
  )

  return(res)
}


find_variable_features <- function(
  object,
  is_peak,
  params
) {
  if (is_peak) {
    Signac::FindTopFeatures(
      object,
      min.cutoff = params$n_features,
      verbose = FALSE
    )
  } else {
    Seurat::FindVariableFeatures(
      object,
      nfeatures = params$n_features,
      verbose = FALSE
    )
  }
}

.create_feature_df <- function(
  features,
  celltype,
  infinite_logfc = TRUE,
  is_peak = FALSE
) {
  df <- data.frame(
    avg_log2FC = rep(if (infinite_logfc) Inf else 1, length(features)),
    p_val_adj = rep(0, length(features)),
    celltype = celltype,
    stringsAsFactors = FALSE
  )

  if (is_peak) {
    df$peak <- features
  } else {
    df$gene <- features
  }

  return(df)
}

.create_empty_feature_df <- function(is_peak = FALSE) {
  df <- data.frame(
    avg_log2FC = numeric(0),
    p_val_adj = numeric(0),
    celltype = character(0),
    stringsAsFactors = FALSE
  )

  if (is_peak) {
    df$peak <- character(0)
  } else {
    df$gene <- character(0)
  }

  return(df)
}

create_feature_matrix <- function(
  features,
  celltypes,
  is_peak = FALSE
) {
  if (length(celltypes) == 1) {
    celltypes <- list(celltypes)
  }
  purrr::map_dfr(
    celltypes,
    ~ .create_feature_df(features, .x, is_peak = is_peak)
  )
}

process_attributes <- function(
  object,
  celltypes,
  genes_markers,
  peaks_markers,
  regions_obj = NULL,
  verbose = TRUE,
  p_value = 0.05
) {
  res <- purrr::map(
    celltypes, function(x) {
      thisutils::log_message(
        "Processing results for {.val {x[1]}}",
        verbose = verbose
      )

      sig_genes <- filter_significant_features(
        genes_markers, x, p_value
      )
      cells <- colnames(object)[Seurat::Idents(object) == x]

      if (!is.null(peaks_markers)) {
        sig_peaks <- filter_significant_features(peaks_markers, x, p_value)

        if (!is.null(regions_obj) && length(regions_obj@ranges) > 0) {
          if (!is.null(sig_peaks) && nrow(sig_peaks) > 0) {
            peak_names <- sig_peaks$peak
            peak_ranges <- Signac::StringToGRanges(peak_names)
            overlaps <- IRanges::findOverlaps(peak_ranges, regions_obj@ranges)
            overlapping_indices <- unique(S4Vectors::queryHits(overlaps))

            if (length(overlapping_indices) > 0) {
              sig_peaks <- sig_peaks[overlapping_indices, ]
              n_peaks_final <- nrow(sig_peaks)
            } else {
              sig_peaks <- NULL
              n_peaks_final <- 0
            }
          } else {
            n_peaks_final <- 0
          }
        } else {
          n_peaks_final <- if (is.null(sig_peaks)) 0 else nrow(sig_peaks)
        }
      } else {
        sig_peaks <- NULL
        n_peaks_final <- 0
      }

      list(
        celltype = x,
        cells = cells,
        n_cells = length(cells),
        genes = sig_genes,
        n_genes = if (is.null(sig_genes)) 0 else nrow(sig_genes),
        peaks = sig_peaks,
        n_peaks = n_peaks_final
      )
    }
  ) |>
    purrr::set_names(celltypes)

  return(res)
}

filter_significant_features <- function(
  markers,
  celltype,
  p_value = 0.05
) {
  sig_features <- markers[
    markers$celltype == celltype &
      markers$p_val_adj < p_value,
  ]
  if (!is.null(sig_features) && nrow(sig_features) > 0) {
    sig_features[order(sig_features$avg_log2FC, decreasing = TRUE), ]
  } else {
    sig_features
  }

  return(sig_features)
}

process_regions <- function(
  object,
  regions = NULL,
  regions_extend = 0,
  peak_assay = NULL,
  exclude_exons = TRUE,
  verbose = TRUE,
  peaks_markers = NULL,
  celltypes = NULL,
  ...
) {
  if (is.null(peak_assay)) {
    res <- methods::new(
      Class = "Regions",
      ranges = GenomicRanges::GRanges(
        seqnames = character(0),
        ranges = IRanges::IRanges(
          start = integer(0),
          end = integer(0)
        )
      ),
      peaks = numeric(0),
      motifs = NULL
    )

    return(res)
  }

  thisutils::log_message(
    "Processing candidate regions",
    verbose = verbose
  )

  gene_annot <- Signac::Annotation(object[[peak_assay]])
  if (is.null(gene_annot)) {
    stop("Please provide a gene annotation for the ChromatinAssay.")
  }

  peak_ranges <- Signac::StringToGRanges(
    rownames(Seurat::GetAssay(object, assay = peak_assay))
  )

  if (!is.null(regions)) {
    if (is.data.frame(regions)) {
      regions <- GenomicRanges::GRanges(
        seqnames = regions[, 1],
        ranges = IRanges::IRanges(
          start = regions[, 2],
          end = regions[, 3]
        )
      )
    }
    if (!inherits(regions, "GRanges")) {
      thisutils::log_message(
        "Regions must be a GRanges object or a data frame with three columns: chrom, start, end.",
        message_type = "error"
      )
    }

    if (regions_extend > 0) {
      regions_extended <- GenomicRanges::resize(
        regions,
        width = GenomicRanges::width(regions) + 2 * regions_extend,
        fix = "center"
      )
      thisutils::log_message(
        "Extended regions by ", regions_extend, " bp on each side",
        verbose = verbose
      )
    } else {
      regions_extended <- regions
    }

    cand_ranges <- regions_extended

    if (!is.null(peaks_markers) && !is.null(celltypes)) {
      significant_peaks <- unique(peaks_markers$peak[!is.na(peaks_markers$peak)])
      if (length(significant_peaks) > 0) {
        peak_indices <- which(names(peak_ranges) %in% significant_peaks)
        if (length(peak_indices) > 0) {
          peak_ranges_filtered <- peak_ranges[peak_indices]
          thisutils::log_message(
            "Filtered to ", length(peak_ranges_filtered), " cell-type-specific peaks",
            verbose = verbose
          )

          cand_olaps <- IRanges::findOverlaps(regions_extended, peak_ranges_filtered)
          if (length(cand_olaps) > 0) {
            cand_ranges <- IRanges::pintersect(
              peak_ranges_filtered[S4Vectors::subjectHits(cand_olaps)],
              regions_extended[S4Vectors::queryHits(cand_olaps)]
            )
            thisutils::log_message(
              "Found ", length(cand_ranges), " candidate regions overlapping with specified regions and significant peaks",
              if (regions_extend > 0) paste0(" (extended by ", regions_extend, " bp)") else "",
              verbose = verbose
            )
          } else {
            thisutils::log_message(
              "No overlaps found between specified regions and significant peaks",
              verbose = verbose,
              message_type = "warning"
            )
            cand_ranges <- GenomicRanges::GRanges()
          }
        } else {
          thisutils::log_message(
            "No significant peaks found, using provided regions directly",
            verbose = verbose
          )
        }
      } else {
        thisutils::log_message(
          "No significant peaks found, using provided regions directly",
          verbose = verbose
        )
      }
    } else {
      thisutils::log_message(
        "Using provided regions directly (", length(cand_ranges), " regions)",
        if (regions_extend > 0) paste0(" (extended by ", regions_extend, " bp)") else "",
        verbose = verbose
      )
    }
  } else {
    if (!is.null(peaks_markers) && !is.null(celltypes)) {
      significant_peaks <- unique(peaks_markers$peak[!is.na(peaks_markers$peak)])
      if (length(significant_peaks) > 0) {
        peak_indices <- which(names(peak_ranges) %in% significant_peaks)
        if (length(peak_indices) > 0) {
          peak_ranges <- peak_ranges[peak_indices]
          thisutils::log_message(
            "Filtered to ", length(peak_ranges), " cell-type-specific peaks",
            verbose = verbose
          )
        }
      }
    }

    cand_ranges <- peak_ranges
  }

  if (exclude_exons) {
    exon_ranges <- gene_annot[gene_annot$type == "exon", ]
    names(exon_ranges@ranges) <- NULL
    exon_ranges <- IRanges::intersect(
      exon_ranges,
      exon_ranges
    )
    exon_ranges <- GenomicRanges::GRanges(
      seqnames = exon_ranges@seqnames,
      ranges = exon_ranges@ranges
    )
    cand_ranges <- GenomicRanges::subtract(
      cand_ranges,
      exon_ranges,
      ignore.strand = TRUE
    ) |> unlist()
  }

  peak_overlaps <- IRanges::findOverlaps(
    cand_ranges,
    peak_ranges
  )
  peak_matches <- S4Vectors::subjectHits(peak_overlaps)

  res <- methods::new(
    Class = "Regions",
    ranges = cand_ranges,
    peaks = peak_matches,
    motifs = NULL
  )

  return(res)
}


create_summary <- function(
  attributes,
  celltypes,
  peak_assay
) {
  res <- do.call(
    rbind,
    lapply(
      celltypes,
      function(x) {
        attr <- attributes[[x]]
        if (!is.null(peak_assay)) {
          data.frame(
            celltype = x,
            cells = attr$n_cells,
            genes = attr$n_genes,
            peaks = attr$n_peaks
          )
        } else {
          data.frame(
            celltype = x,
            cells = attr$n_cells,
            genes = attr$n_genes
          )
        }
      }
    )
  )

  return(res)
}
