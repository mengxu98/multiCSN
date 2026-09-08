#' Rank signed network edges by deletion evidence
#' @param coefficient Signed coefficients.
#' @param deletion_delta_bic Nonnegative deletion delta-BIC values, one per coefficient.
#' @details Descending evidence is grouped against each group's maximum using
#'   absolute difference <= 1e-12 * (1 + abs(maximum)). This numerical rule does
#'   not change fitted support, coefficients or raw deletion evidence.
#' @return Signed ordinal weights; invalid or zero coefficients receive zero.
#' @export
network_ordinal_weight <-
  function(coefficient, deletion_delta_bic) {
    coefficient <- suppressWarnings(as.numeric(coefficient))
    deletion_delta_bic <- suppressWarnings(as.numeric(deletion_delta_bic))
    if (length(coefficient) != length(deletion_delta_bic)) {
      stop("Coefficient and deletion-evidence lengths differ.", call. = FALSE)
    }
    out <- numeric(length(coefficient))
    keep <- is.finite(coefficient) & coefficient != 0 & is.finite(deletion_delta_bic) & deletion_delta_bic >=
      0
    if (!any(keep)) {
      return(out)
    }
    kept_index <- which(keep)
    magnitude <- deletion_delta_bic[kept_index]
    ord <- order(magnitude, decreasing = TRUE, kept_index)
    sorted_magnitude <- magnitude[ord]
    group <- integer(length(ord))
    group[[1]] <- 1L
    group_anchor <- sorted_magnitude[[1]]
    if (length(ord) > 1L) {
      for (i in 2:length(ord)) {
        if (abs(sorted_magnitude[[i]] - group_anchor) > 1e-12 * (1 + abs(group_anchor))) {
          group[[i]] <- group[[i - 1L]] + 1L
          group_anchor <- sorted_magnitude[[i]]
        } else {
          group[[i]] <- group[[i - 1L]]
        }
      }
    }
    n_groups <- max(group)
    mapped <- 1 - (group - 0.5) / n_groups
    out[kept_index[ord]] <- sign(coefficient[kept_index[ord]]) * mapped
    out
  }

.empty_factorized_chain_table <-
  function() {
    data.frame(
      regulator = character(), target = character(), region = character(), tf_region_beta = numeric(),
      tf_region_delta_bic = numeric(), region_gene_beta = numeric(), region_gene_delta_bic = numeric(),
      chain_delta_bic = numeric(), chain_direction = numeric(), weight = numeric(), backbone_selected = logical(),
      backbone_sign_concordant = logical(), stringsAsFactors = FALSE
    )
  }

.signed_ordinal_from_evidence <-
  function(direction, evidence) {
    direction <- sign(as.numeric(direction))
    evidence <- as.numeric(evidence)
    if (length(direction) != length(evidence) || any(!is.finite(direction) | direction == 0) || any(!is.finite(evidence) |
      evidence < 0)) {
      stop("Invalid direction or evidence supplied to the ordinal transform.", call. = FALSE)
    }
    if (!length(direction)) {
      return(numeric())
    }
    network_ordinal_weight(direction, evidence)
  }

.validate_factorized_endpoint <-
  function(x, keys, beta_column = "standardized_beta", evidence_column = "deletion_delta_bic", label = "endpoint") {
    required <- c(keys, beta_column, evidence_column)
    if (!is.data.frame(x) || !all(required %in% names(x))) {
      stop(label, " is missing required columns: ", paste(setdiff(required, names(x)), collapse = ", "),
        call. = FALSE
      )
    }
    if (anyDuplicated(as.data.frame(x)[, keys, drop = FALSE])) {
      stop(label, " contains duplicated endpoint keys.", call. = FALSE)
    }
    beta <- as.numeric(x[[beta_column]])
    evidence <- as.numeric(x[[evidence_column]])
    if (any(!is.finite(beta) | beta == 0) || any(!is.finite(evidence) | evidence < 0)) {
      stop(label, " contains invalid coefficients or deletion evidence.", call. = FALSE)
    }
    invisible(TRUE)
  }

#' Assemble factorized TF-region-target support
#' @param tf_region TF-region edges with regulator, region, standardized_beta and deletion_delta_bic columns.
#' @param region_gene Region-gene edges with region, target, standardized_beta and deletion_delta_bic columns.
#' @param tf_gene Optional TF-gene backbone with regulator, target and the same evidence columns.
#' @return A chain table ranked by minimum endpoint deletion evidence, with backbone concordance flags.
#' @export
factorized_chain_support <-
  function(tf_region, region_gene, tf_gene = NULL) {
    .validate_factorized_endpoint(tf_region, c("regulator", "region"), label = "TF-region endpoint")
    .validate_factorized_endpoint(region_gene, c("region", "target"), label = "region-gene endpoint")
    tf_region <- data.frame(
      regulator = as.character(tf_region$regulator), region = as.character(tf_region$region),
      tf_region_beta = as.numeric(tf_region$standardized_beta), tf_region_delta_bic = as.numeric(tf_region$deletion_delta_bic),
      stringsAsFactors = FALSE
    )
    region_gene <- data.frame(
      region = as.character(region_gene$region), target = as.character(region_gene$target),
      region_gene_beta = as.numeric(region_gene$standardized_beta), region_gene_delta_bic = as.numeric(region_gene$deletion_delta_bic),
      stringsAsFactors = FALSE
    )
    chains <- merge(tf_region, region_gene, by = "region", sort = FALSE)
    if (!nrow(chains)) {
      return(.empty_factorized_chain_table())
    }
    chains$chain_delta_bic <- pmin(chains$tf_region_delta_bic, chains$region_gene_delta_bic)
    chains$chain_direction <- sign(chains$tf_region_beta * chains$region_gene_beta)
    chains$weight <- .signed_ordinal_from_evidence(chains$chain_direction, chains$chain_delta_bic)
    chains$backbone_selected <- FALSE
    chains$backbone_sign_concordant <- FALSE
    if (!is.null(tf_gene)) {
      .validate_factorized_endpoint(tf_gene, c("regulator", "target"), label = "TF-gene endpoint")
      backbone <- data.frame(
        regulator = as.character(tf_gene$regulator), target = as.character(tf_gene$target),
        backbone_direction = sign(as.numeric(tf_gene$standardized_beta)), stringsAsFactors = FALSE
      )
      chains <- merge(chains, backbone, by = c("regulator", "target"), all.x = TRUE, sort = FALSE)
      chains$backbone_selected <- is.finite(chains$backbone_direction)
      chains$backbone_sign_concordant <- chains$backbone_selected & chains$backbone_direction == chains$chain_direction
      chains$backbone_direction <- NULL
    }
    chains <- chains[order(-abs(chains$weight), chains$regulator, chains$region, chains$target), names(.empty_factorized_chain_table()),
      drop = FALSE
    ]
    rownames(chains) <- NULL
    if (anyDuplicated(chains[, c("regulator", "region", "target")])) {
      stop("The assembled factorized network contains duplicated keys.", call. = FALSE)
    }
    chains
  }

#' Project strongest sign-consistent mediated paths
#' @param chains Output of factorized_chain_support.
#' @details Paths within 1e-12 * (1 + abs(maximum evidence)) of the maximum
#'   are numerical ties. Opposite-sign ties are excluded. The representative
#'   region is lexicographically first; pair evidence remains the maximum.
#' @return One strongest path per regulator-target pair; inconsistent direction ties are excluded.
#' @export
strongest_mediated_projection <-
  function(chains) {
    required <- c("regulator", "region", "target", "chain_delta_bic", "weight")
    keys <- c("regulator", "region", "target")
    if (!is.data.frame(chains) || !all(required %in% names(chains)) || anyDuplicated(as.data.frame(chains)[,
      keys,
      drop = FALSE
    ])) {
      stop("Invalid factorized chain table.", call. = FALSE)
    }
    if (!nrow(chains)) {
      return(data.frame(
        regulator = character(), target = character(), region = character(), chain_delta_bic = numeric(),
        strongest_path_ties = integer(), strongest_path_sign_consistent = logical(), weight = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    groups <- split(seq_len(nrow(chains)), interaction(chains$regulator, chains$target,
      drop = TRUE,
      lex.order = TRUE
    ))
    rows <- lapply(groups, function(index) {
      group <- chains[index, , drop = FALSE]
      best <- max(group$chain_delta_bic)
      tied <- group[abs(group$chain_delta_bic - best) <= 1e-12 * (1 + abs(best)), , drop = FALSE]
      directions <- unique(sign(tied$weight))
      chosen <- tied[order(tied$region), , drop = FALSE][1L, , drop = FALSE]
      data.frame(
        regulator = chosen$regulator, target = chosen$target, region = chosen$region, chain_delta_bic = best,
        chain_direction = if (length(directions) == 1L) {
          directions
        } else {
          NA_real_
        }, strongest_path_ties = nrow(tied), strongest_path_sign_consistent = length(directions) ==
          1L, stringsAsFactors = FALSE
      )
    })
    projection <- do.call(rbind, rows)
    projection <- projection[projection$strongest_path_sign_consistent, , drop = FALSE]
    projection$weight <- .signed_ordinal_from_evidence(projection$chain_direction, projection$chain_delta_bic)
    projection$chain_direction <- NULL
    projection <- projection[order(
      -abs(projection$weight), projection$regulator, projection$target,
      projection$region
    ), c(
      "regulator", "target", "region", "chain_delta_bic", "strongest_path_ties",
      "strongest_path_sign_consistent", "weight"
    ), drop = FALSE]
    rownames(projection) <- NULL
    projection
  }
