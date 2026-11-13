#' @include setClass.R
#' @include inferCSN.R

#' @param seed Stored on \code{Network} objects for provenance.
#' @rdname inferCSN
#' @method inferCSN Network
#' @export
#'
#' @examples
#' \dontrun{
#' data("example_matrix", package = "inferCSN")
#' object <- initiate_object(example_matrix)
#' object <- inferCSN(object)
#' }
setMethod(
  f = inferCSN,
  signature = signature(object = "Network"),
  definition = function(object,
                        regulators = NULL,
                        targets = NULL,
                        max_support_size = NULL,
                        cores = 1,
                        verbose = TRUE,
                        seed = 1,
                        r_squared_threshold = 0,
                        method = "greedy_l0",
                        gene_cor_threshold = 0,
                        ...) {
    method <- match.arg(method, "greedy_l0")
    matrix <- thisutils::as_matrix(object@data)
    if (is.null(regulators)) {
      regulators <- object@regulators
      if (length(regulators) == 0) {
        regulators <- colnames(matrix)
      }
    } else {
      regulators <- intersect(
        regulators,
        colnames(matrix)
      )
    }
    if (is.null(targets)) {
      targets <- object@targets
      if (length(targets) == 0) {
        targets <- colnames(matrix)
      }
    } else {
      targets <- intersect(
        targets,
        colnames(matrix)
      )
    }

    model_fits <- fit_models(
      object = matrix,
      regulators = regulators,
      targets = targets,
      gene_cor_threshold = gene_cor_threshold,
      method = method,
      scale = FALSE,
      cores = cores,
      verbose = verbose,
      max_support_size = max_support_size,
      ...
    )
    model_fits <- model_fits[!map_lgl_list(model_fits, is.null)]
    if (length(model_fits) == 0) {
      thisutils::log_message(
        "fitting model failed for all genes.",
        verbose = verbose,
        message_type = "warning"
      )
    }

    coefficients <- map_dfr_list(
      model_fits,
      function(x) x$coefficients,
      .id = "target"
    )
    coefficients <- format_coefs(
      coefficients,
      variable = NULL,
      adjust_method = "fdr"
    )
    corrs <- map_dfr_list(
      model_fits,
      function(x) x$corr,
      .id = "target"
    )
    if (nrow(coefficients) > 0) {
      coefficients <- suppressMessages(
        dplyr::left_join(coefficients, corrs)
      )
    }
    object@metrics <- map_dfr_list(
      model_fits,
      function(x) x$metrics,
      .id = "target"
    )

    object@coefficients <- coefficients
    object@regulators <- regulators
    object@targets <- targets
    object@params <- list(
      method = method,
      seed = seed,
      gene_cor_threshold = gene_cor_threshold,
      r_squared_threshold = r_squared_threshold,
      cores = cores,
      verbose = verbose
    )
    object <- process_Network(
      object,
      r_squared_threshold = r_squared_threshold
    )
    object@network <- export_csn(object)

    return(object)
  }
)

#' @param celltypes Character vector of cell types to infer networks for.
#' @param r_squared_threshold Threshold of \eqn{R^2} used when formatting \code{Network} / \code{Seurat} results.
#' @param network_name Name of the stored network.
#' @param peak_to_gene_method Method to link peaks to genes. One of \code{Signac} or \code{GREAT}.
#' @param upstream Distance upstream of the gene to consider as a potential regulatory region.
#' @param downstream Distance downstream of the gene to consider as a potential regulatory region.
#' @param extend Distance added to the basal regulatory region when \code{peak_to_gene_method = "GREAT"}.
#' @param only_tss Measure distance from the TSS (\code{TRUE}) or from the gene body (\code{FALSE}).
#' @param peak_to_gene_domains \code{GenomicRanges} object with regulatory regions for each gene.
#' @param gene_cor_threshold Threshold for TF-target gene correlation.
#' @param peak_cor_threshold Threshold for binding peak-target gene correlation.
#' @param aggregate_rna_col Metadata column used to aggregate RNA.
#' @param aggregate_peaks_col Metadata column used to aggregate peaks.
#' @param method For matrix input, \code{"greedy_l0"} (default) or the optional
#'   L0Learn penalties \code{"L0"}, \code{"L0L1"}, and \code{"L0L2"}.
#'   The \code{Network} and \code{Seurat} methods support only \code{"greedy_l0"}.
#' @param interaction_term Interaction between TF and binding site: \code{"+"}, \code{":"}, or \code{"*"}.
#' @param scale Whether to z-transform expression and accessibility matrices.
#' @param adjust_method Method for adjusting p-values.
#' @param pseudotime_column Metadata column for pseudotime. If set, runs state-dynamic inference.
#' @param dynamic_features Optional dynamic-gene filter when \code{pseudotime_column} is set.
#'   A list with \code{fit_method} (\code{"pretsa"} by default or \code{"gam"}),
#'   \code{padjust_threshold} (default 0.05), and optional \code{n_candidates}
#'   (default NULL, no cap). Compatible cached statistics are reused; otherwise
#'   \code{inferCSN::select_trend_features()} computes the candidates.
#'
#' @rdname inferCSN
#' @method inferCSN Seurat
#' @export
setMethod(
  f = inferCSN,
  signature = signature(object = "Seurat"),
  definition = function(object,
                        pseudotime = NULL,
                        regulators = NULL,
                        targets = NULL,
                        max_support_size = NULL,
                        cores = 1,
                        verbose = TRUE,
                        r_squared_threshold = 0,
                        celltypes = NULL,
                        network_name = paste0(method, "_network"),
                        peak_to_gene_method = c("Signac", "GREAT"),
                        upstream = 100000,
                        downstream = 0,
                        extend = 1000000,
                        only_tss = FALSE,
                        peak_to_gene_domains = NULL,
                        gene_cor_threshold = 0.1,
                        peak_cor_threshold = 0.,
                        aggregate_rna_col = NULL,
                        aggregate_peaks_col = NULL,
                        method = "greedy_l0",
                        interaction_term = ":",
                        adjust_method = "fdr",
                        scale = FALSE,
                        pseudotime_column = NULL,
                        dynamic_features = NULL,
                        ...) {
    method <- match.arg(method, "greedy_l0")
    if (is.null(pseudotime_column) &&
      is.character(pseudotime) &&
      length(pseudotime) == 1L) {
      pseudotime_column <- pseudotime
    }
    if (!is.null(pseudotime_column)) {
      return(
        state_dynamic(
          object = object,
          pseudotime_column = pseudotime_column,
          method = method,
          r_squared_threshold = r_squared_threshold,
          regulators = regulators,
          targets = targets,
          cores = cores,
          verbose = verbose,
          dynamic_features = dynamic_features,
          max_support_size = max_support_size,
          ...
        )
      )
    }


    peak_to_gene_method <- match.arg(peak_to_gene_method)

    if (is.null(celltypes)) {
      celltypes <- csn_celltypes(object)
    }

    for (celltype in celltypes) {
      thisutils::log_message(
        "Inferring network for: {.val {celltype}}",
        verbose = verbose
      )
      celltype_genes <- get_attribute(
        object,
        celltypes = celltype,
        attribute = "genes"
      )

      celltype_genes <- targets %ss% celltype_genes

      if (length(celltype_genes) == 0) {
        thisutils::log_message(
          "No targets found for {.val {celltype}}, skipping",
          verbose = verbose,
          message_type = "warning"
        )
        next
      }

      object <- fit_models(
        object = object,
        regulators = regulators,
        targets = celltype_genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        gene_cor_threshold = gene_cor_threshold,
        peak_cor_threshold = peak_cor_threshold,
        aggregate_rna_col = aggregate_rna_col,
        aggregate_peaks_col = aggregate_peaks_col,
        method = method,
        interaction_term = interaction_term,
        adjust_method = adjust_method,
        scale = scale,
        verbose = verbose,
        cores = cores,
        celltype = celltype,
        max_support_size = max_support_size,
        ...
      )
    }

    thisutils::log_message(
      "Network inference summary:",
      verbose = verbose
    )

    object <- process_csn(
      object,
      r_squared_threshold = r_squared_threshold
    )
    for (celltype in celltypes) {
      network <- export_csn(
        object,
        celltypes = celltype
      )
      if (!is.null(network)) {
        thisutils::log_message(
          sprintf(
            "   %s: %d edges (%d regulators -> %d targets)",
            celltype,
            nrow(network),
            length(unique(network$regulator)),
            length(unique(network$target))
          ),
          verbose = verbose,
          cli_model = FALSE
        )
      } else {
        thisutils::log_message(
          celltype, ": no edges inferred",
          verbose = verbose,
          message_type = "warning"
        )
      }
    }

    return(object)
  }
)

#' @rdname inferCSN
#' @method inferCSN CSNObject
#' @export
setMethod(
  f = inferCSN,
  signature = signature(object = "CSNObject"),
  definition = function(object, ...) {
    stop_csnobject_runtime()
  }
)
