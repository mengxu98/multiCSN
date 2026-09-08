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

#' @param pfm A PFMatrixList object with position weight matrices.
#' @param genome A BSgenome object.
#' @param motif_tfs A data frame matching motifs with TFs. First column motif, second TF.
#' @param backend Motif scanning backend. `"signac"` uses `Signac::AddMotifs()`;
#'   `"motifmatchr"` uses `motifmatchr::matchMotifs()` directly and wraps the
#'   result into a Signac `Motif` object.
#' @param verbose Display messages.
#' @return A CSNObject object with updated motif info.
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
#' }
#' @rdname find_motifs
#' @export
setMethod(
  f = "find_motifs",
  signature = "Seurat",
  definition = function(object,
                        pfm,
                        genome,
                        motif_tfs = NULL,
                        backend = c("signac", "motifmatchr"),
                        verbose = TRUE,
                        ...) {
    backend <- match.arg(backend)
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
    check_r(
      "motifmatchr",
      verbose = FALSE
    )
    if (!is.null(motif_tfs)) {
      motif2tf <- motif_tfs
    } else {
      utils::data(motif2tf, envir = environment())
    }

    motif2tf <- motif2tf |>
      dplyr::select("motif" = 1, "tf" = 2) |>
      dplyr::distinct() |>
      dplyr::mutate(val = 1) |>
      tidyr::pivot_wider(
        names_from = "tf",
        values_from = val,
        values_fill = 0
      ) |>
      tibble::column_to_rownames("motif") |>
      as.matrix() |>
      Matrix::Matrix(sparse = TRUE)

    assay_genes <- rownames(.csn_get_assay(object, params$rna_assay))
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
    object <- .multicsn_set(object, "tfs", tfs_use)
    regions <- .multicsn_get_regions(object)
    regions@motifs2tfs <- motif2tf[, tfs_use]
    object <- .multicsn_set_regions(object, regions)

    celltypes <- get_attribute(
      object,
      attribute = "celltypes"
    )

    celltype_peaks <- purrr::map(
      celltypes,
      ~ get_attribute(
        object,
        celltypes = .x,
        attribute = "peaks"
      )
    ) |>
      purrr::set_names(celltypes)

    peaks_use <- unique(unlist(celltype_peaks, use.names = FALSE))
    if (length(peaks_use) == 0) {
      thisutils::log_message(
        "No peaks found for motif scanning.",
        verbose = verbose,
        message_type = "warning"
      )
      regions <- .multicsn_get_regions(object)
      regions@motifs <- NULL
      object <- .multicsn_set_regions(object, regions)
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
    result <- .process_peak_motifs(
      peaks = peaks_use,
      label = "all selected peaks",
      genome = genome,
      pfm = pfm,
      backend = backend,
      verbose = verbose
    )
    regions <- .multicsn_get_regions(object)
    regions@motifs <- if (is.null(result)) NULL else result$motifs
    object <- .multicsn_set_regions(object, regions)
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
    .stop_csnobject_runtime()
  }
)

.process_peak_motifs <- function(
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

  .validate_genome_seqinfo(genome)

  peak_ranges <- Signac::StringToGRanges(peaks)
  motif_pos <- suppressWarnings(
    if (identical(backend, "signac")) {
      Signac::AddMotifs(
        object = peak_ranges,
        genome = genome,
        pfm = pfm,
        verbose = FALSE
      )
    } else {
      .create_motif_object_motifmatchr(
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

.create_motif_object_motifmatchr <- function(features, genome, pfm) {
  motif_ix <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = features,
    genome = genome,
    out = "scores"
  )

  motif_matrix <- motifmatchr::motifMatches(motif_ix)
  motif_matrix <- as(motif_matrix, "CsparseMatrix")
  rownames(motif_matrix) <- Signac::GRangesToString(features, sep = c("-", "-"))

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

.validate_genome_seqinfo <- function(genome) {
  seqinfo_ok <- tryCatch(
    {
      GenomeInfoDb::seqinfo(genome)
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  if (seqinfo_ok) {
    return(invisible(TRUE))
  }

  has_seqinfo_pkg <- "package:Seqinfo" %in% search()
  fix_hint <- if (has_seqinfo_pkg) {
    " Detach it using detach('package:Seqinfo', unload = TRUE) and retry."
  } else {
    ""
  }

  stop(
    paste0(
      "GenomeInfoDb::seqinfo(genome) failed for the provided BSgenome object. ",
      "This usually indicates a conflicting 'seqinfo' generic in the R session",
      " (often from the 'Seqinfo' package).",
      fix_hint
    ),
    call. = FALSE
  )
}
