#' @title inferring cell-type specific gene regulatory network
#'
#' @inheritParams inferCSN::inferCSN
#' @param ... Parameters for other methods.
#'
#' @docType methods
#' @rdname inferCSN
#' @return A data table of regulator-target regulatory relationships or a CSNObject.
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
#'
#' genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
#' object <- find_motifs(
#'   object,
#'   pfm = motifs,
#'   motif_tfs = motif2tf,
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
