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
    assay_obj <- tryCatch(csn_get_assay(object, peak_assay), error = function(e) NULL)
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
#' @return Named genomic ranges.
#' @export
parse_peak_ranges <-
  function(peaks, sep = NULL) {
    if (!is.null(sep) && (!is.character(sep) || length(sep) != 2L ||
      anyNA(sep) || any(!nzchar(sep)))) {
      stop("`sep` must contain two non-empty character separators.", call. = FALSE)
    }
    peaks <- as.character(peaks)
    if (length(peaks) == 0) {
      return(GenomicRanges::GRanges())
    }
    fast <- if (is.null(sep)) split_peak_identifiers(peaks) else split_peak_identifiers_with_sep(peaks, sep)
    chromosomes <- fast$chromosome
    starts <- fast$start
    ends <- fast$end
    pending <- which(is.na(chromosomes) | is.na(starts) | is.na(ends))
    matched_pattern <- if (is.null(sep)) rep(FALSE, length(peaks)) else fast$matched
    for (index in if (is.null(sep)) pending else integer()) {
      for (pattern in .peak_identifier_patterns) {
        matched <- regmatches(peaks[index], regexec(pattern, peaks[index], perl = TRUE))[[1L]]
        if (length(matched) == 4L) {
          matched_pattern[index] <- TRUE
          chromosomes[index] <- matched[[2L]]
          starts[index] <- suppressWarnings(as.integer(matched[[3L]]))
          ends[index] <- suppressWarnings(as.integer(matched[[4L]]))
          break
        }
      }
    }
    pending <- which(!matched_pattern & (is.na(chromosomes) | is.na(starts) | is.na(ends)))
    invalid <- pending
    if (length(invalid)) {
      examples <- paste(utils::head(peaks[invalid], 3L), collapse = ", ")
      stop("Could not parse peak identifier(s) as chromosome/start/end: ", examples, call. = FALSE)
    }
    invalid_coordinates <- !is.finite(starts) | !is.finite(ends) | starts < 1L | ends < starts
    if (any(invalid_coordinates)) {
      examples <- paste(utils::head(peaks[invalid_coordinates], 3L), collapse = ", ")
      stop("Peak coordinates must be finite positive integers with end >= start: ", examples, call. = FALSE)
    }
    ranges <- GenomicRanges::GRanges(seqnames = chromosomes, ranges = IRanges::IRanges(
      start = starts,
      end = ends
    ))
    names(ranges) <- peaks
    ranges
  }

.peak_identifier_patterns <- c(
  "^(.+):([0-9]+)-([0-9]+)$",
  "^(.+)-([0-9]+)-([0-9]+)$",
  "^(.+)_([0-9]+)_([0-9]+)$"
)

split_peak_identifiers <- function(peaks) {
  n <- length(peaks)
  chromosome <- rep(NA_character_, n)
  start <- rep(NA_integer_, n)
  end <- rep(NA_integer_, n)
  take <- function(index, pattern, coordinate_delimiter) {
    matched <- regexpr(pattern, peaks[index], perl = TRUE)
    valid <- which(matched > 1L)
    if (!length(valid)) {
      return(invisible(NULL))
    }
    chosen <- index[valid]
    boundary <- as.integer(matched[valid])
    suffix <- substring(peaks[chosen], boundary + 1L)
    coordinates <- data.table::tstrsplit(suffix, coordinate_delimiter, fixed = TRUE)
    chromosome[chosen] <<- substr(peaks[chosen], 1L, boundary - 1L)
    start[chosen] <<- suppressWarnings(as.integer(coordinates[[1L]]))
    end[chosen] <<- suppressWarnings(as.integer(coordinates[[2L]]))
    invisible(NULL)
  }
  take(seq_len(n), ":[0-9]+-[0-9]+$", "-")
  remaining <- which(is.na(chromosome))
  if (length(remaining)) take(remaining, "-[0-9]+-[0-9]+$", "-")
  remaining <- which(is.na(chromosome))
  if (length(remaining)) take(remaining, "_[0-9]+_[0-9]+$", "_")
  list(chromosome = chromosome, start = start, end = end)
}

split_peak_identifiers_with_sep <- function(peaks, sep) {
  last_separator <- function(values, separator) {
    vapply(gregexpr(separator, values, fixed = TRUE), function(matches) {
      positions <- matches[!is.na(matches) & matches > 0L]
      if (length(positions)) positions[[length(positions)]] else NA_integer_
    }, integer(1))
  }
  end_boundary <- last_separator(peaks, sep[[2L]])
  prefix <- substr(peaks, 1L, end_boundary - 1L)
  start_boundary <- last_separator(prefix, sep[[1L]])
  chromosomes <- substr(prefix, 1L, start_boundary - 1L)
  start_text <- substring(prefix, start_boundary + nchar(sep[[1L]]))
  end_text <- substring(peaks, end_boundary + nchar(sep[[2L]]))
  matched <- !is.na(chromosomes) & nzchar(chromosomes) &
    !is.na(start_text) & grepl("^[0-9]+$", start_text) &
    !is.na(end_text) & grepl("^[0-9]+$", end_text)
  chromosomes[!matched] <- NA_character_
  list(
    chromosome = chromosomes,
    start = suppressWarnings(as.integer(start_text)),
    end = suppressWarnings(as.integer(end_text)),
    matched = matched
  )
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

string_to_ranges <- function(regions, sep = NULL) {
  ranges <- parse_peak_ranges(regions, sep = sep)
  names(ranges) <- NULL
  GenomeInfoDb::sortSeqlevels(ranges)
}

ranges_to_string <- function(grange, sep = c("-", "-")) {
  paste0(
    as.character(GenomicRanges::seqnames(grange)), sep[[1]],
    GenomicRanges::start(grange), sep[[2]], GenomicRanges::end(grange)
  )
}
