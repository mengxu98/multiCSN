#' Build a consensus peak set across samples
#'
#' Merges the peak sets of several samples into one non-overlapping universe.
#' Overlapping and adjacent intervals are collapsed with
#' \code{GenomicRanges::reduce()}, optionally restricted to a set of sequence
#' levels and annotated with a \code{Seqinfo} (for example the UCSC hg38
#' chromosome lengths), which is what downstream count projection needs.
#'
#' @param peaks A \code{GRanges}, a character vector of \code{"chr:start-end"}
#'   regions, or a list of either (one element per sample).
#' @param seqlevels Optional character vector of sequence levels to keep.
#' @param seqinfo Optional \code{Seqinfo} used to annotate the kept sequence
#'   levels.
#' @return A named \code{GRanges} whose names are \code{"chr:start-end"}.
#' @export
consensus_peak_set <- function(peaks, seqlevels = NULL, seqinfo = NULL) {
  if (is.list(peaks) && !methods::is(peaks, "GRanges")) {
    ranges <- lapply(peaks, as_peak_ranges)
  } else {
    ranges <- list(as_peak_ranges(peaks))
  }
  if (!length(ranges)) stop("At least one peak set is required.", call. = FALSE)
  combined <- do.call(c, unname(ranges))
  if (!length(combined)) stop("The peak sets are empty.", call. = FALSE)
  consensus <- GenomicRanges::reduce(combined, ignore.strand = TRUE)
  if (!is.null(seqlevels)) {
    keep <- intersect(as.character(seqlevels), GenomeInfoDb::seqlevels(consensus))
    if (!length(keep)) stop("No sequence level of the peak sets is retained.", call. = FALSE)
    consensus <- GenomeInfoDb::keepSeqlevels(consensus, keep, pruning.mode = "coarse")
  }
  if (!is.null(seqinfo)) {
    levels_used <- GenomeInfoDb::seqlevels(consensus)
    missing <- setdiff(levels_used, GenomeInfoDb::seqlevels(seqinfo))
    if (length(missing)) {
      stop("Seqinfo lacks sequence levels: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    GenomeInfoDb::seqinfo(consensus) <- seqinfo[levels_used]
  }
  names(consensus) <- Signac::GRangesToString(grange = consensus, sep = c(":", "-"))
  consensus
}

#' Map sample peaks into a consensus peak set
#'
#' @param peaks A \code{GRanges} or \code{"chr:start-end"} character vector.
#' @param consensus Consensus ranges from \code{\link{consensus_peak_set}}.
#' @param type Overlap type forwarded to \code{GenomicRanges::findOverlaps()}.
#'   The default requires each peak to fall inside one consensus region.
#' @param require_one_to_one Fail when a peak maps to zero or several consensus
#'   regions; that is the condition under which count projection is well defined.
#' @return Integer vector of consensus indices, one per input peak.
#' @export
map_peaks_to_consensus <- function(peaks, consensus, type = "within",
                                   require_one_to_one = TRUE) {
  ranges <- as_peak_ranges(peaks)
  if (!methods::is(consensus, "GRanges")) {
    consensus <- as_peak_ranges(consensus)
  }
  hits <- GenomicRanges::findOverlaps(ranges, consensus, type = type, ignore.strand = TRUE)
  index <- rep(NA_integer_, length(ranges))
  query <- S4Vectors::queryHits(hits)
  subject <- S4Vectors::subjectHits(hits)
  if (anyDuplicated(query)) {
    if (require_one_to_one) {
      stop(
        "Peaks map to several consensus regions; consensus aggregation would ",
        "double-count entries.",
        call. = FALSE
      )
    }
    keep <- !duplicated(query)
    query <- query[keep]
    subject <- subject[keep]
  }
  index[query] <- subject
  if (require_one_to_one && (anyNA(index) || !length(index))) {
    stop("Peaks do not map one-to-one into the consensus regions.", call. = FALSE)
  }
  index
}

#' Project a peak count matrix onto consensus regions
#'
#' @param counts Peak-by-cell count matrix.
#' @param peaks Peak identifiers of \code{counts} (defaults to its row names).
#' @param consensus Consensus ranges from \code{\link{consensus_peak_set}}.
#' @param type,require_one_to_one Passed to \code{\link{map_peaks_to_consensus}}.
#' @return Count matrix with one row per consensus region.
#' @export
project_counts_to_consensus <- function(counts, peaks = rownames(counts), consensus,
                                        type = "within", require_one_to_one = TRUE) {
  if (!length(peaks) || nrow(counts) != length(peaks)) {
    stop("Peak identifiers and count rows disagree.", call. = FALSE)
  }
  index <- map_peaks_to_consensus(
    peaks, consensus, type = type, require_one_to_one = require_one_to_one
  )
  projection <- Matrix::sparseMatrix(
    i = index, j = seq_along(index), x = 1,
    dims = c(length(consensus), length(index)),
    dimnames = list(names(consensus), as.character(peaks))
  )
  projected <- projection %*% counts
  rownames(projected) <- names(consensus)
  colnames(projected) <- colnames(counts)
  projected
}

as_peak_ranges <- function(x) {
  if (methods::is(x, "GRanges")) {
    return(x)
  }
  if (!is.character(x) || !length(x)) {
    stop("Peaks must be GRanges or 'chr:start-end' strings.", call. = FALSE)
  }
  ranges <- Signac::StringToGRanges(regions = x, sep = c(":", "-"))
  ranges
}
