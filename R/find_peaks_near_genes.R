#' Find peaks or regions near gene body or TSS
#'
#' @param peaks A \code{GRanges} object with peak regions.
#' @param genes A \code{GRanges} object with gene coordinates.
#' @param method Character specifying the method to link peak overlapping motif regions to nearby genes.
#' One of 'Signac' or 'GREAT'.
#' @param upstream Integer defining the distance upstream of the gene/TSS to consider.
#' @param downstream Integer defining the distance downstream of the gene/TSS to consider.
#' @param extend Integer defining the distance from the upstream and downstream of the
#' basal regulatory region. Used only when method is 'GREAT'.
#' @param sep Vector of separators to use for genomic string.
#' @param only_tss Logical value. Measure distance from the TSS (\code{TRUE})
#' or from the entire gene body (\code{FALSE}).
#' @param verbose Logical value. Display messages
#' @return A sparse binary Matrix with gene/peak matches.
#' @export
find_peaks_near_genes <- function(
  peaks,
  genes,
  sep = c("-", "-"),
  method = c("Signac", "GREAT"),
  upstream = 100000,
  downstream = 0,
  extend = 1000000,
  only_tss = FALSE,
  verbose = TRUE
) {
  method <- match.arg(method)

  if (method == "Signac") {
    if (only_tss) {
      genes <- IRanges::resize(x = genes, width = 1, fix = "start")
    }
    genes_extended <- suppressWarnings(
      expr = Signac::Extend(
        genes,
        upstream = upstream,
        downstream = downstream
      )
    )
    overlaps <- IRanges::findOverlaps(
      query = peaks,
      subject = genes_extended,
      type = "any",
      select = "all"
    )
    hit_matrix <- Matrix::sparseMatrix(
      i = S4Vectors::queryHits(overlaps),
      j = S4Vectors::subjectHits(overlaps),
      x = 1,
      dims = c(
        length(peaks),
        length(genes_extended)
      )
    )
    rownames(hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name
  } else if (method == "GREAT") {
    utils::data(EnsDb.Hsapiens.v93.annot.UCSC.hg38, envir = environment())
    gene_annot_use <- EnsDb.Hsapiens.v93.annot.UCSC.hg38[
      which(EnsDb.Hsapiens.v93.annot.UCSC.hg38$gene_name %in% genes$gene_name),
    ]
    gene_annot_tss <- dplyr::select(
      tibble::as_tibble(gene_annot_use), "seqnames",
      "start" = "tss", "end" = "tss", "strand"
    )
    tss <- GenomicRanges::GRanges(gene_annot_tss)
    basal_reg <- suppressWarnings(
      expr = Signac::Extend(
        tss,
        upstream = upstream,
        downstream = downstream
      )
    )
    basal_overlaps <- suppressWarnings(
      IRanges::findOverlaps(
        query = peaks,
        subject = basal_reg,
        type = "any",
        select = "all",
        minoverlap = 2
      )
    )
    peak_all <- Signac::GRangesToString(grange = peaks, sep = sep)
    basal_peak_mapped_idx <- S4Vectors::queryHits(basal_overlaps)
    peak_unmapped_idx <- setdiff(seq_along(peak_all), basal_peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)
    gene_bound <- GenomicRanges::GRanges(gene_annot_use)
    body_overlaps <- IRanges::findOverlaps(
      query = peak_unmapped_region,
      subject = gene_bound,
      type = "any",
      select = "all",
      minoverlap = 2
    )
    body_peak_mapped_idx <- peak_unmapped_idx[S4Vectors::queryHits(body_overlaps)]
    peak_mapped_idx <- c(basal_peak_mapped_idx, body_peak_mapped_idx)
    peak_unmapped_idx <- setdiff(seq_along(peak_all), peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)
    extend_reg <- suppressWarnings(
      expr = Signac::Extend(
        basal_reg,
        upstream = extend,
        downstream = extend
      )
    )
    extended_overlaps <- suppressWarnings(
      IRanges::findOverlaps(
        query = peak_unmapped_region,
        subject = extend_reg,
        type = "any",
        select = "all",
        minoverlap = 2
      )
    )
    extended_peak_mapped_idx <- peak_unmapped_idx[S4Vectors::queryHits(extended_overlaps)]
    hit_matrix <- Matrix::sparseMatrix(
      i = c(
        basal_peak_mapped_idx,
        body_peak_mapped_idx,
        extended_peak_mapped_idx
      ),
      j = c(
        S4Vectors::subjectHits(basal_overlaps),
        S4Vectors::subjectHits(body_overlaps),
        S4Vectors::subjectHits(extended_overlaps)
      ),
      x = 1,
      dims = c(length(peaks), length(basal_reg))
    )
    rownames(hit_matrix) <- peak_all
    colnames(hit_matrix) <- c(basal_reg$gene_name)
  }

  return(hit_matrix)
}
