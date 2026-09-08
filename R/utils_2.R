#' @title Format network coefficients
#'
#' @param coefficients A data frame with coefficients
#' @param variable variable
#' @param adjust_method adjust_method
#'
#' @return A data frame.
#'
#' @export
format_coefs <- function(
  coefficients,
  variable = ":",
  adjust_method = "fdr"
) {
  if (dim(coefficients)[1] == 0) {
    return(coefficients)
  }

  if ("pval" %in% colnames(coefficients)) {
    coefficients$padj <- stats::p.adjust(
      coefficients$pval,
      method = adjust_method
    )
  }

  if (is.null(variable)) {
    coefficients <- coefficients %>%
      dplyr::mutate(
        tf = variable,
        region = "1"
      ) %>%
      dplyr::select(tf, target, region, variable, coefficient)

    return(coefficients)
  }

  term_pattern <- paste0("(.+)", variable, "(.+)")
  region_pattern <- "[\\d\\w]+_\\d+_\\d+"
  coefs_use <- coefficients %>%
    dplyr::filter(!variable %in% c("(Intercept)", "Intercept")) %>%
    dplyr::mutate(
      tf_ = stringr::str_replace(variable, term_pattern, "\\1"),
      region_ = stringr::str_replace(variable, term_pattern, "\\2")
    ) %>%
    dplyr::mutate(
      tf = ifelse(stringr::str_detect(tf_, region_pattern), region_, tf_),
      region = ifelse(!stringr::str_detect(tf_, region_pattern), region_, tf_)
    ) %>%
    dplyr::select(-region_, -tf_) %>%
    dplyr::mutate(
      region = stringr::str_replace_all(region, "_", "-"),
      tf = stringr::str_replace_all(tf, "_", "-"),
      target = stringr::str_replace_all(target, "_", "-")
    ) %>%
    dplyr::select(tf, target, region, variable, tidyselect::everything())

  return(coefs_use)
}
