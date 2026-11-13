#' @include setClass.R
#' @include setGenerics.R
NULL

#' @title Scan for motifs in candidate regions
#' @param object The input data, a csn object.
#' @param ... Arguments for other methods
#' @rdname find_motifs
#' @export find_motifs
setGeneric(
  "find_motifs",
  signature = "object",
  function(object, ...) {
    standardGeneric("find_motifs")
  }
)

#' @param pfm A PFMatrixList object with position weight matrices. Defaults to
#'   the bundled \code{\link{motifs}}.
#' @param genome A BSgenome object.
#' @param motif_tfs A data frame matching motifs with TFs. First column motif,
#'   second TF. Defaults to the bundled \code{\link{motif2tf}}.
#' @param backend Motif scanning backend. `"signac"` uses `Signac::AddMotifs()`;
#'   `"motifmatchr"` uses `motifmatchr::matchMotifs()` directly and wraps the
#'   result into a Signac `Motif` object.
#' @param verbose Display messages.
#' @return A CSNObject object with updated motif info.
#'
#' @examples
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
#' }
#' @rdname find_motifs
#' @export
setMethod(
  f = "find_motifs",
  signature = "Seurat",
  definition = function(object,
                        pfm = NULL,
                        genome,
                        motif_tfs = NULL,
                        backend = c("signac", "motifmatchr"),
                        verbose = TRUE,
                        ...) {
    backend <- match.arg(backend)
    if (is.null(pfm)) {
      pfm <- bundled_dataset("motifs")
    }
    params <- Params(object)
    if (is.null(params$peak_assay)) {
      thisutils::log_message(
        "Skipping motif finding as no peak assay was found.",
        verbose = verbose,
        message_type = "warning"
      )
      return(object)
    }

    thisutils::log_message(
      "Adding TF information",
      verbose = verbose
    )
    check_r("motifmatchr", strict = TRUE, verbose = FALSE)
    if (!is.null(motif_tfs)) {
      motif2tf <- motif_tfs
    } else {
      motif2tf <- bundled_dataset("motif2tf")
    }

    motif2tf <- motif_tf_matrix(motif2tf)

    assay_genes <- rownames(csn_get_assay(object, params$rna_assay))
    if (!is.null(params$regulators)) {
      regulators_union <- unique(
        unlist(params$regulators, use.names = FALSE)
      )
      tfs_use <- intersect(
        intersect(regulators_union, assay_genes),
        colnames(motif2tf)
      )
    } else {
      tfs_use <- intersect(
        assay_genes,
        colnames(motif2tf)
      )
    }

    if (length(tfs_use) == 0) {
      stop(
        "None of the provided TFs were found in the dataset. ",
        "Consider providing a custom motif-to-TF map as 'motif_tfs'"
      )
    }
    object <- multicsn_set(object, "tfs", tfs_use)
    regions <- multicsn_get_regions(object)
    regions@motifs2tfs <- motif2tf[, tfs_use]
    object <- multicsn_set_regions(object, regions)

    celltypes <- get_attribute(
      object,
      attribute = "celltypes"
    )

    celltype_peaks <- map_list(
      celltypes,
      ~ get_attribute(
        object,
        celltypes = .x,
        attribute = "peaks"
      )
    ) |>
      stats::setNames(celltypes)

    peaks_use <- unique(unlist(celltype_peaks, use.names = FALSE))
    if (length(peaks_use) == 0) {
      thisutils::log_message(
        "No peaks found for motif scanning.",
        verbose = verbose,
        message_type = "warning"
      )
      regions <- multicsn_get_regions(object)
      regions@motifs <- NULL
      object <- multicsn_set_regions(object, regions)
      return(object)
    }

    thisutils::log_message(
      sprintf(
        "Processing motifs once for %d unique peaks across %d celltypes",
        length(peaks_use),
        length(celltypes)
      ),
      verbose = verbose
    )
    result <- process_peak_motifs(
      peaks = peaks_use,
      label = "all selected peaks",
      genome = genome,
      pfm = pfm,
      backend = backend,
      verbose = verbose
    )
    regions <- multicsn_get_regions(object)
    regions@motifs <- if (is.null(result)) NULL else result$motifs
    object <- multicsn_set_regions(object, regions)
    log_message(
      "Motifs processed",
      message_type = "success",
      verbose = verbose
    )

    return(object)
  }
)

#' @rdname find_motifs
#' @export
setMethod(
  f = "find_motifs",
  signature = "CSNObject",
  definition = function(object, ...) {
    stop_csnobject_runtime()
  }
)

process_peak_motifs <- function(
  peaks,
  label,
  genome,
  pfm,
  backend = c("signac", "motifmatchr"),
  verbose
) {
  backend <- match.arg(backend)
  thisutils::log_message(
    "Processing motifs for ", label,
    verbose = verbose
  )

  if (length(peaks) == 0) {
    thisutils::log_message(
      "no significant peaks found for ", label,
      verbose = verbose,
      message_type = "warning"
    )
    return(NULL)
  }

  validate_genome_seqinfo(genome)

  peak_ranges <- string_to_ranges(peaks)
  motif_pos <- suppressWarnings(
    if (identical(backend, "signac")) {
      Signac::AddMotifs(
        object = peak_ranges,
        genome = genome,
        pfm = pfm,
        verbose = FALSE
      )
    } else {
      create_motif_object_motifmatchr(
        features = peak_ranges,
        genome = genome,
        pfm = pfm
      )
    }
  )

  list(
    motifs = motif_pos,
    n_motifs = ncol(motif_pos),
    peaks = rownames(motif_pos),
    n_peaks = nrow(motif_pos),
    genome = class(genome)[1]
  )
}

create_motif_object_motifmatchr <- function(features, genome, pfm) {
  motif_ix <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = features,
    genome = genome,
    out = "scores",
    p.cutoff = 5e-5,
    bg = "subject",
    w = 7
  )

  motif_matrix <- motifmatchr::motifMatches(motif_ix)
  motif_matrix <- as(motif_matrix, "CsparseMatrix")
  rownames(motif_matrix) <- ranges_to_string(features)

  if (is.null(names(pfm))) {
    warning(
      "No 'names' attribute found in PFMatrixList. Extracting names from individual entries.",
      immediate. = TRUE
    )
    colnames(motif_matrix) <- vapply(
      X = pfm,
      FUN = methods::slot,
      FUN.VALUE = "character",
      "name"
    )
  } else {
    colnames(motif_matrix) <- names(pfm)
  }

  Signac::CreateMotifObject(
    data = motif_matrix,
    pwm = pfm,
    positions = NULL
  )
}

validate_genome_seqinfo <- function(genome) {
  GenomeInfoDb::seqinfo(genome)
  invisible(TRUE)
}

motif_tf_matrix <- function(motif2tf) {
  pairs <- unique(data.frame(
    motif = motif2tf[[1]], tf = motif2tf[[2]], stringsAsFactors = FALSE
  ))
  motifs <- unique(pairs$motif)
  tfs <- unique(pairs$tf)
  incidence <- matrix(0, length(motifs), length(tfs), dimnames = list(motifs, tfs))
  incidence[cbind(match(pairs$motif, motifs), match(pairs$tf, tfs))] <- 1
  Matrix::Matrix(incidence, sparse = TRUE)
}

bundled_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "multiCSN", envir = environment)
  environment[[name]]
}

#' Prepare a multiome object for layered network inference
#'
#' Builds the candidate catalogue used by \code{\link{fit_layered_network}}:
#' it initializes the multiCSN state without marker filtering (every RNA gene
#' is a candidate target), aligns the gene annotation with the peak chromosome
#' names and scans the selected peaks for motifs.
#'
#' @param object A Seurat object with normalized RNA and TF-IDF chromatin assays.
#' @param group.by Metadata column with the cell groups to fit, or
#'   \code{"aggregate"} for one group with all cells.
#' @param genome A BSgenome object matching the peak coordinates.
#' @param rna_assay,peak_assay Assay names.
#' @param targets Candidate target genes; defaults to every RNA gene.
#' @param pfm,motif_tfs Motif matrices and motif-to-TF map; default to the
#'   bundled \code{\link{motifs}} and \code{\link{motif2tf}}.
#' @param backend Motif scanning backend, see \code{\link{find_motifs}}.
#' @param exclude_exons Remove exonic peaks from the candidates.
#' @param align_annotation Add peak chromosomes missing from the gene
#'   annotation to its sequence levels.
#' @param verbose Print progress messages.
#' @return The object with candidate genes, peaks, transcription factors and
#'   motif matches.
#' @examples
#' \dontrun{
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' object <- prepare_candidates(object, group.by = "cell_type", genome = genome)
#' fit <- fit_layered_network(object, cell_group = "HSC", preset = "state")
#' }
#' @export
prepare_candidates <- function(object,
                               group.by = "aggregate",
                               genome,
                               rna_assay = "RNA",
                               peak_assay = "ATAC",
                               targets = NULL,
                               pfm = NULL,
                               motif_tfs = NULL,
                               backend = c("motifmatchr", "signac"),
                               exclude_exons = TRUE,
                               align_annotation = TRUE,
                               verbose = TRUE) {
  backend <- match.arg(backend)
  if (!methods::is(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  missing_assays <- setdiff(c(rna_assay, peak_assay), SeuratObject::Assays(object))
  if (length(missing_assays)) {
    stop("Missing assay(s): ", paste(missing_assays, collapse = ", "), call. = FALSE)
  }
  check_layered_assay_layers(object, c(rna_assay, peak_assay))
  if (is.null(targets)) {
    targets <- rownames(object[[rna_assay]])
  }
  if (isTRUE(align_annotation)) {
    object <- align_annotation_seqlevels(object, peak_assay)
  }
  object <- initiate_object(
    object,
    group.by = group.by, targets = targets,
    filter_mode = "unfiltered", filter_by = "aggregate",
    rna_assay = rna_assay, peak_assay = peak_assay,
    regions = NULL, regions_extend = 0, exclude_exons = exclude_exons,
    verbose = verbose
  )
  find_motifs(
    object,
    pfm = pfm, motif_tfs = motif_tfs, genome = genome,
    backend = backend, verbose = verbose
  )
}

align_annotation_seqlevels <- function(object, peak_assay) {
  assay <- object[[peak_assay]]
  if (!methods::is(assay, "ChromatinAssay")) {
    return(object)
  }
  annotation <- Signac::Annotation(assay)
  if (is.null(annotation)) {
    stop("The chromatin assay needs a gene annotation.", call. = FALSE)
  }
  chromosomes <- GenomeInfoDb::seqlevels(parse_peak_ranges(rownames(assay)))
  missing <- setdiff(chromosomes, GenomeInfoDb::seqlevels(annotation))
  if (length(missing)) {
    genome <- unique(stats::na.omit(GenomeInfoDb::genome(annotation)))
    GenomeInfoDb::seqlevels(annotation) <- c(GenomeInfoDb::seqlevels(annotation), missing)
    if (length(genome) == 1L) {
      GenomeInfoDb::genome(annotation) <- genome
    }
    Signac::Annotation(object[[peak_assay]]) <- annotation
  }
  object
}
