# -*- coding: utf-8 -*-

#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib multiCSN, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#'
#' @description
#' Inferring cell-specific gene regulatory networks from single-cell multiome data.
#'
#' @author Meng Xu (Maintainer), \email{mengxu98@qq.com}
#'
#' @source \url{https://mengxu98.github.io/multiCSN/}
#'
#' @md
#' @docType package
#' @name multiCSN-package
"_PACKAGE"

utils::globalVariables(
  c(
    "abs_weight",
    "boundary",
    "component",
    "curve_y",
    "edge_jaccard",
    "end",
    "expression_z",
    "feature_label",
    "fill_group",
    "gene_mid",
    "key_tf_score",
    "label",
    "label_x",
    "label_y",
    "left",
    "n_targets",
    "pair_label",
    "peak_end",
    "peak_mid",
    "peak_start",
    "peak_support",
    "peak_supported",
    "plot_size",
    "rewiring_score",
    "right",
    "score",
    "score_plot",
    "similarity",
    "start",
    "state_from",
    "state_id",
    "state_to",
    "state_type",
    "target_label",
    "target_score",
    "tf_label",
    "value"
  )
)

#' @title The logo of multiCSN
#'
#' @description
#' The multiCSN logo, using ASCII or Unicode characters
#' Use [cli::ansi_strip] to get rid of the colors.
#'
#' @md
#' @param unicode Unicode symbols on UTF-8 platforms.
#' Default is [cli::is_utf8_output].
#'
#' @return
#' A character vector with class `multicsn_logo`.
#'
#' @references
#' \url{https://github.com/tidyverse/tidyverse/blob/main/R/logo.R}
#'
#' @export
#' @examples
#' multicsn_logo()
multicsn_logo <- function(unicode = cli::is_utf8_output()) {
  logo <- c(
    "          0          1        2             3     4
                           ____  _  ___________ _   __
          ____ ___  __  __/ / /_(_)/ ____/ ___// | / /
         / __ `__ ./ / / / / __/ // /    .__ ./  |/ /
        / / / / / / /_/ / / /_/ // /___ ___/ / /|  /
       /_/ /_/ /_/.__,_/_/.__/_/ .____//____/_/ |_/
      5               6      7        8          9"
  )

  hexa <- c("*", ".", "o", "*", ".", "o", "*", ".", "o", "*")
  if (unicode) {
    hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]
  }

  cols <- c(
    "red", "yellow", "green", "magenta", "cyan",
    "yellow", "green", "white", "magenta", "cyan"
  )

  col_hexa <- mapply(
    function(x, y) cli::make_ansi_style(y)(x),
    hexa, cols,
    SIMPLIFY = FALSE
  )

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(
    cli::col_blue(logo),
    class = "multicsn_logo"
  )
}

#' @title Print logo
#'
#' @param x Input information.
#' @param ... Other parameters.
#'
#' @return Print the ASCII logo
#'
#' @method print multicsn_logo
#'
#' @export
#'
print.multicsn_logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  verbose <- thisutils::get_verbose()
  if (isTRUE(verbose)) {
    version <- utils::packageVersion(pkgname)
    date <- utils::packageDate(pkgname)
    url <- utils::packageDescription(
      pkgname,
      fields = "URL"
    )

    msg <- paste0(
      cli::col_grey(strrep("-", 60)),
      "\n",
      cli::col_blue("Version: ", version, " (", date, " update)"),
      "\n",
      cli::col_blue("Website: ", cli::style_italic(url)),
      "\n\n",
      cli::col_grey("This message can be suppressed by:"),
      "\n",
      cli::col_grey("  suppressPackageStartupMessages(library(multiCSN))"),
      "\n",
      cli::col_grey("  or options(log_message.verbose = FALSE)"),
      "\n",
      cli::col_grey(strrep("-", 60))
    )

    packageStartupMessage(multicsn_logo())
    packageStartupMessage(msg)
  }
}
