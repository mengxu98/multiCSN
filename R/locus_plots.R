#' Normalize an ATAC assay to a Signac ChromatinAssay
#'
#' @param object A \code{Seurat} object.
#' @param peak_assay Name of the assay to normalize. Defaults to
#'   \code{Params(object)$peak_assay} or \code{"ATAC"} when present.
#' @param new_assay Optional assay name to write to. Defaults to overwriting
#'   \code{peak_assay}.
#' @param fragments Optional fragment object/path passed to
#'   \code{Signac::CreateChromatinAssay()}.
#' @param annotation Optional gene annotation \code{GRanges}. Defaults to the
#'   current assay annotation when available.
#' @param genome Optional genome metadata forwarded to
#'   \code{Signac::CreateChromatinAssay()}.
#' @param sep Separator used in peak names.
#'
#' @return A \code{Seurat} object with a \code{ChromatinAssay}.
#' @export
normalize_chromatin_assay <- function(
  object,
  peak_assay = NULL,
  new_assay = NULL,
  fragments = NULL,
  annotation = NULL,
  genome = NULL,
  sep = c("-", "-")
) {
  if (!methods::is(object, "Seurat")) {
    stop("normalize_chromatin_assay: object must be a Seurat object.", call. = FALSE)
  }

  peak_assay <- peak_assay %ss% (tryCatch(Params(object)$peak_assay, error = function(e) NULL) %ss% if ("ATAC" %in% names(object@assays)) "ATAC" else NULL)
  if (is.null(peak_assay) || !peak_assay %in% names(object@assays)) {
    stop("normalize_chromatin_assay: unable to identify peak assay.", call. = FALSE)
  }
  new_assay <- new_assay %ss% peak_assay

  assay_obj <- object[[peak_assay]]
  annotation <- annotation %ss% tryCatch(Signac::Annotation(assay_obj), error = function(e) NULL)

  if (inherits(assay_obj, "ChromatinAssay")) {
    if (!is.null(annotation)) {
      Signac::Annotation(assay_obj) <- annotation
    }
    if (!is.null(fragments)) {
      Signac::Fragments(assay_obj) <- fragments
    }
    object[[new_assay]] <- assay_obj
    return(object)
  }

  counts <- tryCatch(
    Seurat::GetAssayData(object, assay = peak_assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts) || nrow(counts) == 0) {
    stop("normalize_chromatin_assay: counts layer is required to create a ChromatinAssay.", call. = FALSE)
  }

  peak_ranges <- tryCatch(
    Signac::StringToGRanges(rownames(counts), sep = sep),
    error = function(e) NULL
  )
  if (is.null(peak_ranges) || length(peak_ranges) != nrow(counts)) {
    stop("normalize_chromatin_assay: failed to parse peak coordinates from rownames.", call. = FALSE)
  }

  chrom_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    ranges = peak_ranges,
    sep = sep,
    genome = genome,
    fragments = fragments,
    annotation = annotation
  )

  data_layer <- tryCatch(
    Seurat::GetAssayData(object, assay = peak_assay, layer = "data"),
    error = function(e) NULL
  )
  if (!is.null(data_layer) && all(dim(data_layer) == dim(counts))) {
    chrom_assay <- SeuratObject::SetAssayData(
      object = chrom_assay,
      layer = "data",
      new.data = data_layer
    )
  }

  object[[new_assay]] <- chrom_assay
  object
}

.resolve_peak_assay_for_locus <- function(object, peak_assay = NULL) {
  assay_names <- names(object@assays)
  chromatin_assays <- assay_names[vapply(object@assays, function(x) inherits(x, "ChromatinAssay"), logical(1))]
  peak_assay <- peak_assay %ss% if ("ATAC" %in% chromatin_assays) {
    "ATAC"
  } else {
    NULL
  }
  peak_assay <- peak_assay %ss% (tryCatch(Params(object)$peak_assay, error = function(e) NULL))
  if (!is.null(peak_assay) && peak_assay %in% assay_names && !inherits(object[[peak_assay]], "ChromatinAssay") && length(chromatin_assays) > 0) {
    peak_assay <- if ("ATAC" %in% chromatin_assays) "ATAC" else chromatin_assays[[1]]
  }
  peak_assay <- peak_assay %ss% if (length(chromatin_assays) > 0) chromatin_assays[[1]] else NULL
  peak_assay <- peak_assay %ss% if ("ATAC" %in% assay_names) "ATAC" else NULL
  if (is.null(peak_assay) || !peak_assay %in% names(object@assays)) {
    stop("Unable to find a peak assay for locus plotting.", call. = FALSE)
  }
  peak_assay
}

.build_gene_peak_map_local <- function(
  object,
  genes,
  peak_assay = NULL,
  upstream = 100000,
  downstream = 0,
  only_tss = FALSE
) {
  peak_assay <- .resolve_peak_assay_for_locus(object, peak_assay)
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(NULL)
  }

  gene_annot <- Signac::Annotation(object[[peak_assay]])
  if (is.null(gene_annot)) {
    stop("Peak assay is missing genomic annotations.", call. = FALSE)
  }
  gene_annot <- gene_annot[gene_annot$gene_name %in% genes]
  gene_annot <- gene_annot[!duplicated(gene_annot$gene_name)]
  if (length(gene_annot) == 0) {
    return(NULL)
  }

  peak_ranges <- Signac::StringToGRanges(rownames(object[[peak_assay]]))
  peak_gene <- find_peaks_near_genes(
    peaks = peak_ranges,
    genes = gene_annot,
    method = "Signac",
    upstream = upstream,
    downstream = downstream,
    only_tss = only_tss
  )
  peak_gene[, intersect(colnames(peak_gene), genes), drop = FALSE]
}

#' Plot TF-supported peak overlap across ordered states
#'
#' @param object A \code{Seurat} object.
#' @param key_target_tbl Output of \code{\link{prioritize_features}} with \code{type = "targets"} and \code{strategy = "integrated"}.
#' @param tfs TFs to visualize.
#' @param pseudotime_column Pseudotime column used for state assignment.
#' @param peak_assay Peak assay name.
#' @param top_targets_per_tf Number of target genes retained per TF.
#'
#' @return A \code{ggplot} object.
#' @export
plot_tf_peak_overlap_scenic_style <- function(
  object,
  key_target_tbl,
  tfs,
  pseudotime_column,
  peak_assay = NULL,
  top_targets_per_tf = 6
) {
  peak_assay <- .resolve_peak_assay_for_locus(object, peak_assay)
  tfs <- unique(as.character(tfs))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  if (length(tfs) == 0 || is.null(key_target_tbl) || nrow(key_target_tbl) == 0) {
    stop("plot_tf_peak_overlap_scenic_style: tfs and key_target_tbl are required.", call. = FALSE)
  }

  target_tbl <- key_target_tbl[key_target_tbl$tf %in% tfs, , drop = FALSE]
  if (nrow(target_tbl) == 0) {
    stop("plot_tf_peak_overlap_scenic_style: no TF-target rows remained after filtering.", call. = FALSE)
  }
  target_tbl <- do.call(rbind, lapply(split(target_tbl, target_tbl$tf), function(df) {
    df <- df[order(df$target_score, decreasing = TRUE), , drop = FALSE]
    df[seq_len(min(top_targets_per_tf, nrow(df))), , drop = FALSE]
  }))
  rownames(target_tbl) <- NULL

  gene_peak_map <- .build_gene_peak_map_local(
    object,
    genes = unique(target_tbl$target),
    peak_assay = peak_assay
  )
  peak_tf <- build_peak_tf_incidence(object, tfs = tfs)
  if (is.null(gene_peak_map) || is.null(peak_tf)) {
    stop("plot_tf_peak_overlap_scenic_style: unable to construct peak support matrices.", call. = FALSE)
  }

  common_peaks <- intersect(rownames(gene_peak_map), rownames(peak_tf))
  gene_peak_map <- gene_peak_map[common_peaks, , drop = FALSE]
  peak_tf <- peak_tf[common_peaks, tfs, drop = FALSE]

  tf_peak_sets <- lapply(tfs, function(tf) {
    genes_use <- unique(target_tbl$target[target_tbl$tf == tf])
    genes_use <- intersect(genes_use, colnames(gene_peak_map))
    if (length(genes_use) == 0 || !tf %in% colnames(peak_tf)) {
      return(character(0))
    }
    gp_sub <- gene_peak_map[, genes_use, drop = FALSE]
    motif_peaks <- rownames(peak_tf)[as.vector(peak_tf[, tf]) > 0]
    linked_peaks <- rownames(gp_sub)[Matrix::rowSums(gp_sub != 0) > 0]
    intersect(motif_peaks, linked_peaks)
  })
  names(tf_peak_sets) <- tfs
  tf_peak_sets <- tf_peak_sets[lengths(tf_peak_sets) > 0]
  if (length(tf_peak_sets) == 0) {
    stop("plot_tf_peak_overlap_scenic_style: no TF-supported peak sets available.", call. = FALSE)
  }

  state_levels <- get_ordered_state_levels(unique(as.character(key_target_tbl$anchor_state)))
  group_sets <- list()
  tf_names <- names(tf_peak_sets)
  for (tf in tf_names) {
    other_peaks <- unique(unlist(tf_peak_sets[setdiff(tf_names, tf)], use.names = FALSE))
    group_sets[[paste0(tf, " only")]] <- setdiff(tf_peak_sets[[tf]], other_peaks)
  }
  if (length(tf_names) >= 2) {
    pair_idx <- utils::combn(tf_names, 2, simplify = FALSE)
    for (pair in pair_idx) {
      pair_peaks <- Reduce(intersect, tf_peak_sets[pair])
      if (length(pair_peaks) > 0) {
        group_sets[[paste(pair, collapse = " + ")]] <- pair_peaks
      }
    }
  }
  group_sets <- group_sets[lengths(group_sets) > 0]
  if (length(group_sets) == 0) {
    stop("plot_tf_peak_overlap_scenic_style: no overlap groups available.", call. = FALSE)
  }

  heat_df <- do.call(rbind, lapply(names(group_sets), function(group_name) {
    peaks_use <- group_sets[[group_name]]
    do.call(rbind, lapply(state_levels, function(sid) {
      acc <- compute_state_peak_accessibility(
        object,
        state_id = sid,
        peaks = peaks_use,
        pseudotime_column = pseudotime_column,
        peak_assay = peak_assay
      )
      data.frame(
        group = group_name,
        state_id = sid,
        mean_accessibility = mean(acc, na.rm = TRUE),
        n_peaks = length(peaks_use),
        stringsAsFactors = FALSE
      )
    }))
  }))

  heat_df$group <- factor(heat_df$group, levels = rev(names(group_sets)))
  heat_df$state_id <- factor(heat_df$state_id, levels = state_levels)

  ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = state_id, y = group, fill = mean_accessibility)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_peaks),
      size = 3,
      color = "black"
    ) +
    ggplot2::scale_fill_gradient(low = "white", high = "#08519C") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "TF peak-set overlap panel",
      subtitle = "Fill: mean state accessibility; text: number of supporting peaks",
      x = "State",
      y = "Peak-set group",
      fill = "ATAC"
    )
}

.ensure_state_group_column <- function(
  object,
  pseudotime_column,
  column_name = NULL,
  include_unassigned = FALSE
) {
  column_name <- column_name %ss% paste0("codex_state_", pseudotime_column)
  if (column_name %in% colnames(object@meta.data)) {
    return(list(object = object, column = column_name))
  }

  state_df <- assign_cells_to_states(
    object,
    pseudotime_column = pseudotime_column,
    keep_only = TRUE,
    inferable_only = FALSE
  )
  object[[column_name]] <- "unassigned"
  if (!is.null(state_df) && nrow(state_df) > 0) {
    idx <- match(state_df$cell, colnames(object))
    idx <- idx[!is.na(idx)]
    object@meta.data[idx, column_name] <- as.character(
      state_df$state_id[match(colnames(object)[idx], state_df$cell)]
    )
    state_levels <- get_ordered_state_levels(unique(as.character(state_df$state_id)))
    object[[column_name]][, 1] <- factor(
      object[[column_name]][, 1],
      levels = c(state_levels, "unassigned")
    )
  }

  if (isFALSE(include_unassigned)) {
    keep <- as.character(object[[column_name]][, 1]) != "unassigned"
    keep[is.na(keep)] <- FALSE
    object <- object[, keep, drop = FALSE]
  }

  list(object = object, column = column_name)
}

.build_signac_links_for_target <- function(
  object,
  target_gene,
  key_target_tbl,
  tfs = NULL,
  pseudotime_column,
  peak_assay = NULL
) {
  peak_assay <- .resolve_peak_assay_for_locus(object, peak_assay)
  target_gene <- as.character(target_gene[[1]])
  target_tbl <- key_target_tbl[key_target_tbl$target == target_gene, , drop = FALSE]
  if (!is.null(tfs)) {
    target_tbl <- target_tbl[target_tbl$tf %in% tfs, , drop = FALSE]
  }
  if (nrow(target_tbl) == 0) {
    stop(sprintf("No target records found for %s.", target_gene), call. = FALSE)
  }

  tfs_use <- unique(target_tbl$tf)
  state_levels <- get_ordered_state_levels(unique(as.character(target_tbl$anchor_state)))
  gene_peak_map <- .build_gene_peak_map_local(
    object,
    genes = target_gene,
    peak_assay = peak_assay
  )
  peak_tf <- build_peak_tf_incidence(object, tfs = tfs_use)
  if (is.null(gene_peak_map) || is.null(peak_tf)) {
    stop("Unable to construct peak-to-gene or peak-to-TF mappings.", call. = FALSE)
  }

  common_peaks <- intersect(rownames(gene_peak_map), rownames(peak_tf))
  gene_peak_map <- gene_peak_map[common_peaks, target_gene, drop = FALSE]
  peak_tf <- peak_tf[common_peaks, tfs_use, drop = FALSE]
  peak_keep <- rownames(gene_peak_map)[Matrix::rowSums(gene_peak_map != 0) > 0]
  if (length(peak_keep) == 0) {
    stop(sprintf("No linked peaks found for %s.", target_gene), call. = FALSE)
  }

  gene_annot <- Signac::Annotation(object[[peak_assay]])
  gene_gr <- gene_annot[gene_annot$gene_name == target_gene]
  gene_gr <- gene_gr[1]
  gene_mid <- floor((GenomicRanges::start(gene_gr) + GenomicRanges::end(gene_gr)) / 2)

  link_df <- do.call(rbind, lapply(seq_len(nrow(target_tbl)), function(i) {
    tf <- target_tbl$tf[[i]]
    tf_score <- as.numeric(target_tbl$target_score[[i]] %ss% 0)
    tf_peak_support <- as.numeric(target_tbl$peak_support[[i]] %ss% 0)
    peaks_tf <- rownames(peak_tf)[as.vector(peak_tf[, tf, drop = FALSE]) > 0]
    peaks_tf <- intersect(peaks_tf, peak_keep)
    if (length(peaks_tf) == 0) {
      return(NULL)
    }

    do.call(rbind, lapply(peaks_tf, function(pk) {
      acc_by_state <- vapply(state_levels, function(sid) {
        acc <- compute_state_peak_accessibility(
          object,
          state_id = sid,
          peaks = pk,
          pseudotime_column = pseudotime_column,
          peak_assay = peak_assay
        )
        as.numeric(acc[[1]] %ss% 0)
      }, numeric(1))
      peak_ranges <- Signac::StringToGRanges(pk)
      peak_mid <- floor((GenomicRanges::start(peak_ranges) + GenomicRanges::end(peak_ranges)) / 2)
      data.frame(
        peak = pk,
        tf = tf,
        chr = as.character(GenomicRanges::seqnames(peak_ranges)),
        peak_start = GenomicRanges::start(peak_ranges),
        peak_end = GenomicRanges::end(peak_ranges),
        peak_mid = peak_mid,
        gene = target_gene,
        gene_mid = gene_mid,
        score = tf_score + 0.25 * tf_peak_support + max(acc_by_state, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
  if (is.null(link_df) || nrow(link_df) == 0) {
    stop(sprintf("No TF-supported links found for %s.", target_gene), call. = FALSE)
  }

  link_df_raw <- link_df
  link_df <- stats::aggregate(
    score ~ peak + chr + peak_start + peak_end + peak_mid + gene + gene_mid,
    data = link_df_raw,
    FUN = max
  )
  tf_map <- stats::aggregate(
    tf ~ peak + chr + peak_start + peak_end + peak_mid + gene + gene_mid,
    data = link_df_raw,
    FUN = function(x) paste(sort(unique(as.character(x))), collapse = ", ")
  )
  tf_n <- stats::aggregate(
    tf ~ peak + chr + peak_start + peak_end + peak_mid + gene + gene_mid,
    data = link_df_raw,
    FUN = function(x) length(unique(as.character(x)))
  )
  link_df <- merge(link_df, tf_map, by = c("peak", "chr", "peak_start", "peak_end", "peak_mid", "gene", "gene_mid"), all.x = TRUE, sort = FALSE)
  link_df <- merge(link_df, tf_n, by = c("peak", "chr", "peak_start", "peak_end", "peak_mid", "gene", "gene_mid"), all.x = TRUE, sort = FALSE)
  colnames(link_df)[colnames(link_df) == "tf.x"] <- "tf_label"
  colnames(link_df)[colnames(link_df) == "tf.y"] <- "tf_n"
  link_df <- link_df[order(-link_df$score, -link_df$tf_n, abs(link_df$peak_mid - link_df$gene_mid)), , drop = FALSE]

  links_gr <- GenomicRanges::GRanges(
    seqnames = link_df$chr,
    ranges = IRanges::IRanges(
      start = pmin(link_df$peak_mid, link_df$gene_mid),
      end = pmax(link_df$peak_mid, link_df$gene_mid)
    )
  )
  links_gr$score <- link_df$score
  links_gr$gene <- link_df$gene
  links_gr$peak <- link_df$peak
  links_gr$peak_start <- link_df$peak_start
  links_gr$peak_end <- link_df$peak_end
  links_gr$target_gene <- target_gene
  links_gr$tf_label <- link_df$tf_label
  links_gr$tf_n <- link_df$tf_n

  list(links = links_gr, gene = gene_gr, link_table = link_df)
}

.build_peak_pseudocoverage_plot <- function(
  object,
  region,
  group_column,
  peak_assay,
  title = "Pseudo-coverage",
  n_bins = 200,
  normalize_by_group = TRUE
) {
  peak_ranges <- Signac::StringToGRanges(rownames(object[[peak_assay]]))
  region_gr <- Signac::StringToGRanges(region)
  hits <- GenomicRanges::findOverlaps(peak_ranges, region_gr)
  if (length(hits) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No peaks in region", size = 5)
    )
  }

  peak_idx <- unique(S4Vectors::queryHits(hits))
  peaks_use <- rownames(object[[peak_assay]])[peak_idx]
  peak_sub <- peak_ranges[peak_idx]
  counts <- Seurat::GetAssayData(object, assay = peak_assay, layer = "counts")[peaks_use, , drop = FALSE]
  groups <- as.character(object@meta.data[[group_column]])
  groups[is.na(groups) | !nzchar(groups)] <- "unassigned"
  group_levels <- unique(groups)
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)
  bin_edges <- seq(region_start, region_end, length.out = n_bins + 1)
  peak_mid <- floor((GenomicRanges::start(peak_sub) + GenomicRanges::end(peak_sub)) / 2)
  bin_id <- cut(
    peak_mid,
    breaks = bin_edges,
    include.lowest = TRUE,
    labels = FALSE
  )
  bin_centers <- floor((bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2)

  cov_df <- do.call(rbind, lapply(group_levels, function(grp) {
    cell_idx <- which(groups == grp)
    if (length(cell_idx) == 0) {
      return(NULL)
    }
    vals <- Matrix::rowMeans(counts[, cell_idx, drop = FALSE])
    bin_vals <- tapply(vals, bin_id, sum, na.rm = TRUE)
    y <- numeric(length(bin_centers))
    if (length(bin_vals) > 0) {
      y[as.integer(names(bin_vals))] <- as.numeric(bin_vals)
    }
    if (isTRUE(normalize_by_group) && max(y, na.rm = TRUE) > 0) {
      y <- y / max(y, na.rm = TRUE)
    }
    data.frame(
      group = grp,
      x = bin_centers,
      value = y,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(cov_df) || nrow(cov_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No grouped counts available", size = 5)
    )
  }

  cov_df$group <- factor(cov_df$group, levels = rev(group_levels))
  ggplot2::ggplot(cov_df, ggplot2::aes(x = x, y = value)) +
    ggplot2::geom_area(fill = "#9ECAE1", alpha = 0.9, linewidth = 0) +
    ggplot2::geom_line(color = "#4A4A4A", linewidth = 0.25) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Genomic position",
      y = if (isTRUE(normalize_by_group)) "Scaled peak counts" else "Mean peak counts"
    )
}

.build_gene_model_plot <- function(
  gene_gr,
  region,
  gene_label = NULL
) {
  region_gr <- Signac::StringToGRanges(region)
  gene_label <- gene_label %ss% as.character(gene_gr$gene_name[[1]] %ss% "gene")
  gene_df <- data.frame(
    start = GenomicRanges::start(gene_gr),
    end = GenomicRanges::end(gene_gr),
    gene = gene_label,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(gene_df) +
    ggplot2::geom_segment(
      ggplot2::aes(x = start, xend = end, y = 0, yend = 0),
      linewidth = 1.2,
      color = "#238B45"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = (start + end) / 2, y = 0),
      size = 2,
      color = "#238B45"
    ) +
    ggplot2::annotate(
      "text",
      x = mean(c(gene_df$start, gene_df$end)),
      y = -0.15,
      label = gene_label,
      size = 4,
      color = "#238B45"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(GenomicRanges::start(region_gr), GenomicRanges::end(region_gr)),
      ylim = c(-0.25, 0.25),
      clip = "off"
    ) +
    ggplot2::theme_void()
}

.build_link_arc_plot <- function(
  link_table,
  region,
  title = "Links",
  top_n = 8,
  show_tf_labels = TRUE
) {
  region_gr <- Signac::StringToGRanges(region)
  if (is.null(link_table) || nrow(link_table) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No links available", size = 5)
    )
  }

  arc_df <- link_table
  arc_df <- arc_df[arc_df$peak_mid >= GenomicRanges::start(region_gr) & arc_df$peak_mid <= GenomicRanges::end(region_gr), , drop = FALSE]
  if (nrow(arc_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No links inside region", size = 5)
    )
  }
  arc_df <- arc_df[order(arc_df$peak_mid), , drop = FALSE]
  if (!is.null(top_n) && nrow(arc_df) > top_n) {
    arc_df <- arc_df[order(-arc_df$score, -arc_df$tf_n, abs(arc_df$peak_mid - arc_df$gene_mid)), , drop = FALSE]
    arc_df <- arc_df[seq_len(top_n), , drop = FALSE]
    arc_df <- arc_df[order(arc_df$peak_mid), , drop = FALSE]
  }
  arc_df$curve_y <- seq(0.12, 0.35, length.out = nrow(arc_df))
  arc_df$label_y <- 0.505

  label_df <- NULL
  if (isTRUE(show_tf_labels)) {
    label_source <- arc_df[order(-arc_df$score, -arc_df$tf_n), , drop = FALSE]
    label_df <- label_source[!duplicated(label_source$tf_label), , drop = FALSE]
    label_df <- label_df[nzchar(label_df$tf_label), , drop = FALSE]
    label_df$label_y <- seq(0.505, 0.535, length.out = nrow(label_df))
  }

  p <- ggplot2::ggplot(arc_df) +
    ggplot2::geom_curve(
      ggplot2::aes(
        x = peak_mid,
        y = curve_y,
        xend = gene_mid,
        yend = curve_y,
        color = score
      ),
      curvature = 0.22,
      linewidth = 0.85,
      alpha = 0.95
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = peak_start, xmax = peak_end, ymin = 0.40, ymax = 0.49),
      fill = "#9ECAE1",
      alpha = 0.85,
      color = NA
    ) +
    ggplot2::coord_cartesian(
      xlim = c(GenomicRanges::start(region_gr), GenomicRanges::end(region_gr)),
      ylim = c(0, 0.52),
      clip = "off"
    ) +
    ggplot2::scale_color_gradient(low = "#9ECAE1", high = "#08519C") +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, color = "Score")
  if (isTRUE(show_tf_labels) && !is.null(label_df) && nrow(label_df) > 0) {
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = peak_mid, y = label_y, label = tf_label),
      size = 2.5,
      angle = 90,
      vjust = -0.1,
      hjust = 0,
      color = "grey20"
    )
  }
  p
}

#' Plot a target-gene locus using Signac when available, with robust fallbacks
#'
#' @param object A \code{Seurat} object.
#' @param target_gene Target gene symbol.
#' @param key_target_tbl Output of \code{\link{prioritize_features}} with \code{type = "targets"} and \code{strategy = "integrated"}.
#' @param tfs Optional TF subset used to build links.
#' @param pseudotime_column Pseudotime column used for state assignment.
#' @param assay RNA assay used for expression tracks when \code{CoveragePlot()}
#'   is available.
#' @param peak_assay Peak assay name.
#' @param extend_upstream,extend_downstream Genomic extension around the gene.
#' @param auto_extend Whether to expand the plotting window to include all
#'   inferred linked peaks for the target gene.
#' @param include_unassigned Whether to keep unassigned cells in grouped tracks.
#' @param top_links Maximum number of links shown in fallback link plots.
#' @param show_tf_labels Whether to display TF labels above linked peaks in the
#'   fallback link panel.
#' @param normalize_tracks Whether to scale each grouped pseudo-coverage track
#'   to a maximum of 1 in fallback mode.
#'
#' @return A patchwork/ggplot object.
#' @export
plot_target_locus_signac <- function(
  object,
  target_gene,
  key_target_tbl,
  tfs = NULL,
  pseudotime_column,
  assay = "RNA",
  peak_assay = NULL,
  extend_upstream = 5000,
  extend_downstream = 5000,
  auto_extend = TRUE,
  include_unassigned = FALSE,
  top_links = 8,
  show_tf_labels = TRUE,
  normalize_tracks = TRUE
) {
  peak_assay <- .resolve_peak_assay_for_locus(object, peak_assay)
  state_info <- .ensure_state_group_column(
    object = object,
    pseudotime_column = pseudotime_column,
    include_unassigned = include_unassigned
  )
  object <- state_info$object

  link_res <- .build_signac_links_for_target(
    object = object,
    target_gene = target_gene,
    key_target_tbl = key_target_tbl,
    tfs = tfs,
    pseudotime_column = pseudotime_column,
    peak_assay = peak_assay
  )

  region_chr <- as.character(GenomicRanges::seqnames(link_res$gene))
  if (isTRUE(auto_extend) && nrow(link_res$link_table) > 0) {
    region_start <- min(
      GenomicRanges::start(link_res$gene),
      link_res$link_table$peak_start,
      na.rm = TRUE
    ) - extend_upstream
    region_end <- max(
      GenomicRanges::end(link_res$gene),
      link_res$link_table$peak_end,
      na.rm = TRUE
    ) + extend_downstream
  } else {
    region_start <- GenomicRanges::start(link_res$gene) - extend_upstream
    region_end <- GenomicRanges::end(link_res$gene) + extend_downstream
  }
  gene_region <- paste0(
    region_chr,
    "-",
    max(1, region_start),
    "-",
    region_end
  )

  fragments_ok <- FALSE
  fragments_obj <- tryCatch(Signac::Fragments(object[[peak_assay]]), error = function(e) NULL)
  if (!is.null(fragments_obj) && length(fragments_obj) > 0) {
    fragments_ok <- TRUE
  }
  chromatin_assay_ok <- inherits(object[[peak_assay]], "ChromatinAssay")

  if (isTRUE(chromatin_assay_ok)) {
    old_links <- tryCatch(Signac::Links(object[[peak_assay]]), error = function(e) NULL)
    Signac::Links(object[[peak_assay]]) <- link_res$links
    on.exit(
      {
        if (is.null(old_links)) {
          Signac::Links(object[[peak_assay]]) <- GenomicRanges::GRanges()
        } else {
          Signac::Links(object[[peak_assay]]) <- old_links
        }
      },
      add = TRUE
    )
  }

  if (isTRUE(fragments_ok) && isTRUE(chromatin_assay_ok)) {
    return(
      Signac::CoveragePlot(
        object = object,
        region = gene_region,
        assay = peak_assay,
        expression.assay = assay,
        annotation = TRUE,
        peaks = TRUE,
        links = TRUE,
        group.by = state_info$column,
        extend.upstream = 0,
        extend.downstream = 0
      ) +
        patchwork::plot_annotation(
          title = paste0("Signac locus plot: ", target_gene),
          subtitle = "CoveragePlot with annotation, peaks, and custom peak-gene links"
        )
    )
  }

  p_cov <- .build_peak_pseudocoverage_plot(
    object = object,
    region = gene_region,
    group_column = state_info$column,
    peak_assay = peak_assay,
    title = paste0("Signac locus plot: ", target_gene),
    normalize_by_group = normalize_tracks
  )
  if (isTRUE(chromatin_assay_ok)) {
    p_links <- tryCatch(
      Signac::LinkPlot(
        object = object,
        region = gene_region,
        assay = peak_assay,
        min.cutoff = 0
      ),
      error = function(e) NULL
    )
    p_annot <- tryCatch(
      Signac::AnnotationPlot(
        object = object[[peak_assay]],
        region = gene_region
      ),
      error = function(e) NULL
    )
  } else {
    p_links <- NULL
    p_annot <- NULL
  }

  if (is.null(p_links) || is.null(p_annot)) {
    p_links <- .build_link_arc_plot(
      link_table = link_res$link_table,
      region = gene_region,
      title = "Peak-gene links",
      top_n = top_links,
      show_tf_labels = show_tf_labels
    )
    p_annot <- .build_gene_model_plot(
      gene_gr = link_res$gene,
      region = gene_region,
      gene_label = target_gene
    )
  }

  p_cov / p_links / p_annot + patchwork::plot_layout(heights = c(4.2, 1.9, 1.0))
}
