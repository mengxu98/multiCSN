#' Read assay layer data
#' @param object A Seurat object or a legacy CSNObject.
#' @param ... Arguments passed to SeuratObject::LayerData.
#' @return An assay layer matrix.
#' @export
get_layer_data <-
  function(object, ...) {
    if (methods::is(object, "CSNObject")) {
      return(SeuratObject::LayerData(object@data, ...))
    }
    if (methods::is(object, "Seurat")) {
      return(SeuratObject::LayerData(object, ...))
    }
    stop("object must be a Seurat or CSNObject object.", call. = FALSE)
  }

#' Read peak assay annotations
#' @param object A Seurat object or a legacy CSNObject.
#' @param peak_assay Name of the peak assay.
#' @return Genomic annotations, or NULL if unavailable.
#' @export
get_peak_annotation <-
  function(object, peak_assay) {
    assay_obj <- tryCatch(.csn_get_assay(object, peak_assay), error = function(e) NULL)
    annotation <- NULL
    if (!is.null(assay_obj) && methods::is(assay_obj, "ChromatinAssay")) {
      annotation <- tryCatch(Signac::Annotation(assay_obj), error = function(e) NULL)
    }
    if (!is.null(annotation)) {
      return(annotation)
    }
    if (methods::is(object, "CSNObject")) {
      object <- object@data
    }
    object@misc$multiCSN_annotations[[peak_assay]] %ss% object@misc$multiCSN_input$annotation
  }

#' Parse peak identifiers
#' @param peaks Character peak identifiers in chromosome:start-end, chromosome-start-end or chromosome_start_end form.
#' @param sep Retained for compatibility with peak identifier callers.
#' @return Named genomic ranges.
#' @export
parse_peak_ranges <-
  function(peaks, sep = c("-", "-")) {
    peaks <- as.character(peaks)
    if (length(peaks) == 0) {
      return(GenomicRanges::GRanges())
    }
    patterns <- c("^(.+):([0-9]+)-([0-9]+)$", "^(.+)-([0-9]+)-([0-9]+)$", "^(.+)_([0-9]+)_([0-9]+)$")
    parse_one <- function(peak) {
      for (pattern in patterns) {
        match <- regmatches(peak, regexec(pattern, peak, perl = TRUE))[[1]]
        if (length(match) == 4L) {
          return(match[2:4])
        }
      }
      character(0)
    }
    parsed <- lapply(peaks, parse_one)
    invalid <- which(lengths(parsed) != 3L)
    if (length(invalid)) {
      examples <- paste(utils::head(peaks[invalid], 3L), collapse = ", ")
      stop("Could not parse peak identifier(s) as chromosome/start/end: ", examples, call. = FALSE)
    }
    parsed <- do.call(rbind, parsed)
    starts <- suppressWarnings(as.integer(parsed[, 2]))
    ends <- suppressWarnings(as.integer(parsed[, 3]))
    invalid_coordinates <- !is.finite(starts) | !is.finite(ends) | starts < 1L | ends < starts
    if (any(invalid_coordinates)) {
      examples <- paste(utils::head(peaks[invalid_coordinates], 3L), collapse = ", ")
      stop("Peak coordinates must be finite positive integers with end >= start: ", examples, call. = FALSE)
    }
    ranges <- GenomicRanges::GRanges(seqnames = parsed[, 1], ranges = IRanges::IRanges(
      start = starts,
      end = ends
    ))
    names(ranges) <- peaks
    ranges
  }

#' Map peak sequence names to a genome
#' @param peak_ranges Genomic ranges to map.
#' @param genome A genome object with sequence levels and lengths.
#' @return A list of mapped ranges, support flags and original/mapped sequence names.
#' @export
map_peak_ranges_to_genome <-
  function(peak_ranges, genome) {
    genome_levels <- GenomeInfoDb::seqlevels(genome)
    original_seqnames <- as.character(GenomicRanges::seqnames(peak_ranges))
    resolve_seqname <- function(seqname) {
      if (seqname %in% genome_levels) {
        return(seqname)
      }
      direct_candidates <- unique(c(paste0("chr", seqname), sub("^chr", "", seqname)))
      direct_hit <- direct_candidates[direct_candidates %in% genome_levels]
      if (length(direct_hit) == 1L) {
        return(direct_hit)
      }
      accession <- sub("\\.([0-9]+)$", "v\\1", seqname)
      accession_pattern <- paste0("(^|_)", accession, "(_|$)")
      accession_hit <- genome_levels[grepl(accession_pattern, genome_levels)]
      if (length(accession_hit) == 1L) {
        return(accession_hit)
      }
      seqname
    }
    mapped_seqnames <- vapply(original_seqnames, resolve_seqname, FUN.VALUE = character(1), USE.NAMES = FALSE)
    mapped_ranges <- GenomicRanges::GRanges(
      seqnames = mapped_seqnames, ranges = GenomicRanges::ranges(peak_ranges),
      strand = GenomicRanges::strand(peak_ranges)
    )
    names(mapped_ranges) <- names(peak_ranges)
    S4Vectors::mcols(mapped_ranges) <- S4Vectors::mcols(peak_ranges)
    supported <- mapped_seqnames %in% genome_levels
    genome_lengths <- GenomeInfoDb::seqlengths(genome)
    sequence_lengths <- unname(genome_lengths[mapped_seqnames])
    supported <- supported & as.integer(IRanges::start(mapped_ranges)) >= 1L & is.finite(sequence_lengths) &
      as.integer(IRanges::end(mapped_ranges)) <= sequence_lengths
    list(ranges = mapped_ranges, supported = supported, original_seqnames = original_seqnames, mapped_seqnames = mapped_seqnames)
  }

#' Collapse repeated gene annotations into peak-gene domains
#' @param peak_annotation_matrix Sparse peak-by-annotation matrix, with annotation gene names as column names.
#' @return A gene-by-peak sparse matrix with duplicate gene annotations summed.
#' @export
collapse_peak_gene_domains <-
  function(peak_annotation_matrix) {
    if (!methods::is(peak_annotation_matrix, "sparseMatrix")) {
      stop("peak_annotation_matrix must be a sparse Matrix.", call. = FALSE)
    }
    peak_names <- rownames(peak_annotation_matrix)
    annotation_gene <- colnames(peak_annotation_matrix)
    if (is.null(peak_names) || is.null(annotation_gene) || anyNA(peak_names) || anyNA(annotation_gene) ||
      any(!nzchar(peak_names)) || any(!nzchar(annotation_gene))) {
      stop("Peak and annotation gene names must be complete.", call. = FALSE)
    }
    gene_names <- unique(annotation_gene)
    nonzero <- Matrix::summary(peak_annotation_matrix)
    if (!nrow(nonzero)) {
      return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(
        length(gene_names),
        length(peak_names)
      ), dimnames = list(gene_names, peak_names)))
    }
    methods::as(Matrix::sparseMatrix(
      i = match(annotation_gene[nonzero$j], gene_names), j = nonzero$i,
      x = nonzero$x, dims = c(length(gene_names), length(peak_names)), dimnames = list(
        gene_names,
        peak_names
      )
    ), "generalMatrix")
  }
