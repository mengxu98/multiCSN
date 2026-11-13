#' @title inferring cell-type specific gene regulatory network
#'
#' @description Infers cell-specific gene regulatory networks. Matrix methods come from
#' \pkg{inferCSN}; \code{Network} and \code{Seurat} methods are implemented here.
#'
#' @param object Expression matrix, \code{Network}, \code{Seurat}, or \code{CSNObject}.
#'   Matrix input is handled by \code{\link[inferCSN:inferCSN]{inferCSN::inferCSN()}}.
#' @param pseudotime Optional pseudotime vector or branch matrix for matrix input.
#'   For \code{Seurat} objects, a single column name is treated as \code{pseudotime_column}.
#' @param regulators,targets Optional gene subsets.
#' @param max_support_size Optional support-size limit for matrix inference.
#' @param lag_fraction Fractional lag used when \code{lag_steps} is \code{NULL}.
#' @param lag_steps Optional integer lag for matrix inference.
#' @param cores Number of workers.
#' @param verbose Whether to report progress.
#' @param ... Additional method arguments.
#'
#' @docType methods
#' @rdname inferCSN
#' @return A regulator-target edge table for matrix input, a \code{Network} object, or a \code{Seurat} object.
#' @details The Seurat method with chromatin peaks fits TF-by-region interaction
#' predictors through \code{fit_models}; projecting those coefficients does not
#' fit independent TF-gene, TF-region and region-gene models. For factorized
#' analyses, fit the three endpoint responses independently and use
#' \code{factorized_chain_support} to assemble their selected supports.
#' @export
#'
#' @examples
#' data("example_matrix", package = "inferCSN")
#' network_table_1 <- inferCSN(
#'   example_matrix
#' )
#'
#' network_table_2 <- inferCSN(
#'   example_matrix,
#'   cores = 2
#' )
#'
#' head(network_table_1)
#'
#' identical(
#'   network_table_1,
#'   network_table_2
#' )
#'
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g0")
#' )
#'
#' \dontrun{
#' data("example_ground_truth", package = "inferCSN")
#' calculate_metrics(
#'   network_table_1,
#'   example_ground_truth,
#'   return_plot = TRUE
#' )
#' }
#'
#' # multiome input: motif scanning comes first
#' \dontrun{
#' # Supply a Seurat object containing RNA and chromatin accessibility assays.
#' data(motifs)
#'
#' object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
#' object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
#' object <- initiate_object(
#'   object,
#'   group.by = "CellType",
#'   rna_assay = "RNA",
#'   peak_assay = "peaks",
#'   verbose = FALSE
#' )
#'
#' genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
#' object <- find_motifs(
#'   object,
#'   pfm = motifs,
#'   genome = genome,
#'   backend = "motifmatchr",
#'   verbose = FALSE
#' )
#'
#' object <- inferCSN(object, cores = 2, verbose = FALSE)
#' networks <- export_csn(object)
#' head(networks[[1]])
#' }
inferCSN <- methods::getGeneric(
  "inferCSN",
  where = asNamespace("inferCSN")
)
