#' Pseudobulk a perturbation screen into a paired design
#'
#' Sums single-cell counts of one or more modalities per
#' replicate x condition sample, which is the input of the paired effect model
#' used to validate an inferred network against a perturbation screen. The
#' control label keeps its own sample in every replicate, so the resulting
#' design is paired and complete by construction.
#'
#' @param cells A matrix (features x cells) or a named list of matrices, one per
#'   modality (for example RNA and ATAC). All modalities must share the cells.
#' @param sample_info Data frame with \code{cell_id}, \code{replicate},
#'   \code{timepoint} and \code{target} columns.
#' @param targets Perturbation targets to aggregate; defaults to every target of
#'   \code{sample_info} except the control label, in order of appearance.
#' @param control_label Label identifying the unperturbed control cells.
#' @param roles Named character vector mapping \code{control} and
#'   \code{perturbed} to the role labels stored in the sample metadata.
#' @param require_complete Fail when a replicate lacks a control or target
#'   sample.
#' @param verbose Print progress messages.
#' @return A list with \code{sample_metadata} and \code{pseudobulk} (one matrix
#'   per modality, columns ordered exactly like \code{sample_metadata$sample_id}).
#' @export
pseudobulk_perturbation <- function(cells, sample_info, targets = NULL,
                                    control_label = "NT",
                                    roles = c(
                                      control = "network_training_control",
                                      perturbed = "heldout_perturbation"
                                    ),
                                    require_complete = TRUE, verbose = TRUE) {
  if (!is.list(cells) || is.matrix(cells) || methods::is(cells, "Matrix")) {
    cells <- list(cells)
    names(cells) <- "counts"
  }
  if (is.null(names(cells)) || any(!nzchar(names(cells))) || anyDuplicated(names(cells))) {
    stop("Modalities must be a named, unique list.", call. = FALSE)
  }
  required <- c("cell_id", "replicate", "timepoint", "target")
  if (!is.data.frame(sample_info) || !all(required %in% names(sample_info))) {
    stop("sample_info must have columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  sample_info$cell_id <- as.character(sample_info$cell_id)
  if (anyDuplicated(sample_info$cell_id)) {
    stop("sample_info contains duplicated cell identifiers.", call. = FALSE)
  }
  if (is.null(targets)) {
    targets <- setdiff(unique(as.character(sample_info$target)), control_label)
  }
  group_levels <- c(control_label, as.character(targets))
  if (anyDuplicated(group_levels) || length(group_levels) < 3L) {
    stop("At least one control and two distinct targets are required.", call. = FALSE)
  }
  for (modality in names(cells)) {
    missing <- setdiff(sample_info$cell_id, colnames(cells[[modality]]))
    if (length(missing)) {
      stop(
        "Modality '", modality, "' lacks ", length(missing),
        " annotated cell(s), e.g. ", missing[[1L]],
        call. = FALSE
      )
    }
  }
  replicates <- unique(as.character(sample_info$replicate))
  pseudobulk <- lapply(names(cells), function(x) NULL)
  names(pseudobulk) <- names(cells)
  sample_rows <- vector("list", length(replicates))
  for (index in seq_along(replicates)) {
    replicate <- replicates[[index]]
    block <- sample_info[as.character(sample_info$replicate) == replicate, , drop = FALSE]
    timepoints <- unique(as.character(block$timepoint))
    if (length(timepoints) != 1L) {
      stop("Replicate '", replicate, "' maps to several timepoints.", call. = FALSE)
    }
    group_index <- match(as.character(block$target), group_levels)
    if (require_complete && anyNA(group_index)) {
      stop("Replicate '", replicate, "' carries an undeclared target label.", call. = FALSE)
    }
    if (require_complete && !all(seq_along(group_levels) %in% group_index)) {
      stop(
        "Replicate '", replicate, "' lacks a sample for: ",
        paste(group_levels[!seq_along(group_levels) %in% group_index], collapse = ", "),
        call. = FALSE
      )
    }
    design <- Matrix::sparseMatrix(
      i = seq_len(nrow(block)), j = group_index, x = 1,
      dims = c(nrow(block), length(group_levels)),
      dimnames = list(block$cell_id, group_levels)
    )
    for (modality in names(cells)) {
      aggregate <- cells[[modality]][, block$cell_id, drop = FALSE] %*% design
      colnames(aggregate) <- paste(replicate, group_levels, sep = "__")
      pseudobulk[[modality]] <- if (is.null(pseudobulk[[modality]])) {
        aggregate
      } else {
        cbind(pseudobulk[[modality]], aggregate)
      }
    }
    sample_rows[[index]] <- data.frame(
      sample_id = paste(replicate, group_levels, sep = "__"),
      replicate = replicate,
      timepoint = timepoints,
      target = group_levels,
      role = ifelse(group_levels == control_label, roles[["control"]], roles[["perturbed"]]),
      cells = as.integer(table(factor(block$target, levels = group_levels))),
      stringsAsFactors = FALSE
    )
    if (isTRUE(verbose)) {
      thisutils::log_message("Aggregated replicate ", replicate, verbose = TRUE)
    }
  }
  sample_metadata <- do.call(rbind, sample_rows)
  rownames(sample_metadata) <- sample_metadata$sample_id
  for (modality in names(pseudobulk)) {
    if (!identical(colnames(pseudobulk[[modality]]), sample_metadata$sample_id)) {
      stop("Pseudobulk columns and sample metadata disagree for ", modality, call. = FALSE)
    }
  }
  list(sample_metadata = sample_metadata, pseudobulk = pseudobulk)
}

#' Fit paired perturbation effects on pseudobulk counts
#'
#' Fits one quasi-likelihood negative-binomial model per timepoint and target
#' with a paired design (\code{~ replicate + condition}) and reports per-feature
#' log fold changes and FDR values for the perturbed-minus-control contrast.
#'
#' @param pseudobulk Named list of pseudobulk matrices (see
#'   \code{\link{pseudobulk_perturbation}}).
#' @param sample_metadata Sample metadata with \code{sample_id}, \code{replicate},
#'   \code{timepoint}, \code{target}, \code{cells}.
#' @param targets,timepoints Subsets to fit; defaults to every target and
#'   timepoint present.
#' @param control_label Label of the control samples.
#' @param effect_definition Label stored in the effect tables.
#' @param min_replicates Minimum number of paired replicates per contrast.
#' @param output_dir Optional directory receiving one \code{<timepoint>__<target>.tsv.gz}
#'   table per contrast; when \code{NULL} the tables are returned in memory.
#' @param output_subdirs Optional named vector mapping modality names to the
#'   subdirectory created below \code{output_dir} (defaults to the modality name).
#' @param verbose Print progress messages.
#' @return A list with \code{effects} (per modality, one table or path per
#'   contrast) and \code{summary} (one row per modality x timepoint x target).
#' @export
fit_perturbation_effects <- function(pseudobulk, sample_metadata, targets = NULL,
                                     timepoints = NULL,
                                     control_label = "NT",
                                     effect_definition = "heldout_perturbed_minus_NT_control",
                                     min_replicates = 3L, output_dir = NULL,
                                     output_subdirs = NULL,
                                     verbose = TRUE) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for fit_perturbation_effects().", call. = FALSE)
  }
  required <- c("sample_id", "replicate", "timepoint", "target")
  if (!is.data.frame(sample_metadata) || !all(required %in% names(sample_metadata))) {
    stop("sample_metadata must have columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (is.null(names(pseudobulk)) || any(!nzchar(names(pseudobulk)))) {
    stop("pseudobulk must be a named list of matrices.", call. = FALSE)
  }
  if (!is.null(output_subdirs)) {
    if (is.null(names(output_subdirs)) || any(!nzchar(names(output_subdirs)))) {
      stop("output_subdirs must be a named vector of modality directories.", call. = FALSE)
    }
  }
  sample_metadata$sample_id <- as.character(sample_metadata$sample_id)
  if (is.null(targets)) {
    targets <- setdiff(unique(as.character(sample_metadata$target)), control_label)
  }
  if (is.null(timepoints)) {
    timepoints <- unique(as.character(sample_metadata$timepoint))
  }
  tasks <- expand.grid(
    timepoint = as.character(timepoints), target = as.character(targets),
    stringsAsFactors = FALSE
  )
  contrasts <- split(tasks, seq_len(nrow(tasks)))
  summary_rows <- list()
  effects <- lapply(names(pseudobulk), function(x) list())
  names(effects) <- names(pseudobulk)
  for (modality in names(pseudobulk)) {
    counts <- pseudobulk[[modality]]
    for (task_index in seq_along(contrasts)) {
      task <- contrasts[[task_index]]
      keep <- as.character(sample_metadata$timepoint) == task$timepoint &
        as.character(sample_metadata$target) %in% c(control_label, task$target)
      info <- sample_metadata[keep, , drop = FALSE]
      info <- info[order(info$replicate, info$target), , drop = FALSE]
      paired <- table(info$replicate)
      if (nrow(info) != 2L * length(unique(info$replicate)) ||
        length(unique(info$replicate)) < min_replicates || any(paired != 2L)) {
        stop(
          "Incomplete paired pseudobulk design for ", task$target, " at ", task$timepoint,
          call. = FALSE
        )
      }
      info$condition <- factor(
        ifelse(info$target == control_label, "control", "perturbed"),
        levels = c("control", "perturbed")
      )
      info$replicate <- factor(info$replicate)
      design <- stats::model.matrix(~ replicate + condition, data = info)
      y <- edgeR::DGEList(counts = counts[, info$sample_id, drop = FALSE])
      tested <- edgeR::filterByExpr(y, design = design)
      if (!any(tested)) {
        stop("No tested features for ", task$target, " at ", task$timepoint, call. = FALSE)
      }
      y <- y[tested, , keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y, method = "TMM")
      y <- edgeR::estimateDisp(y, design, robust = TRUE)
      fit <- edgeR::glmQLFit(y, design, robust = TRUE)
      test <- edgeR::glmQLFTest(fit, coef = ncol(design))
      result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
      result$feature <- rownames(result)
      result <- result[, c("feature", setdiff(colnames(result), "feature")), drop = FALSE]
      result$target <- task$target
      result$timepoint <- task$timepoint
      result$effect_definition <- effect_definition
      slug <- paste(gsub(" ", "_", task$timepoint, fixed = TRUE), task$target, sep = "__")
      output_path <- NULL
      if (!is.null(output_dir)) {
        subdir <- if (is.null(output_subdirs)) modality else {
          value <- output_subdirs[[modality]]
          if (is.null(value) || !nzchar(value)) modality else value
        }
        modality_dir <- file.path(output_dir, subdir)
        dir.create(modality_dir, recursive = TRUE, showWarnings = FALSE)
        output_path <- file.path(modality_dir, paste0(slug, ".tsv.gz"))
        data.table::fwrite(result, output_path, sep = "\t", quote = FALSE, na = "NA")
      }
      effects[[modality]][[slug]] <- if (is.null(output_path)) result else output_path
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        modality = modality,
        timepoint = task$timepoint,
        target = task$target,
        paired_replicates = length(unique(info$replicate)),
        control_cells = sum(info$cells[info$target == control_label]),
        perturbation_cells = sum(info$cells[info$target == task$target]),
        release_features = nrow(counts),
        tested_features = nrow(result),
        significant_fdr_0_05 = sum(result$FDR <= 0.05, na.rm = TRUE),
        significant_up_fdr_0_05 = sum(result$FDR <= 0.05 & result$logFC > 0, na.rm = TRUE),
        significant_down_fdr_0_05 = sum(result$FDR <= 0.05 & result$logFC < 0, na.rm = TRUE),
        output_file = if (is.null(output_path)) NA_character_ else output_path,
        stringsAsFactors = FALSE
      )
      if (isTRUE(verbose)) {
        thisutils::log_message(
          "Fitted ", modality, " contrast for ", task$target, " at ", task$timepoint,
          verbose = TRUE
        )
      }
    }
  }
  list(effects = effects, summary = do.call(rbind, summary_rows))
}
