#' @title Fit model
#'
#' @description Fits a sparse regression model through [inferCSN::fit_greedy_l0()].
#'
#' @param formula A model formula.
#' @param data A data frame containing the response and predictors.
#' @param method The sole supported method, `"greedy_l0"`.
#' @param max_support_size Maximum selected support size.
#' @param ... Additional arguments passed to [fit_srm2()].
#'
#' @return A list containing the model, fit metrics, and selected coefficients.
#' @export
#' @examples
#' data("example_matrix", package = "inferCSN")
#' df <- as.data.frame(example_matrix)
#' fit_model(g1 ~ ., data = df)
fit_model <- function(formula, data, method = "greedy_l0",
                      max_support_size = NULL, ...) {
  match.arg(method, "greedy_l0")
  fit_srm2(formula, data, max_support_size = max_support_size, ...)
}

#' @title Fit a sparse regression model
#'
#' @description
#' Fits a sparse regression model using [inferCSN::fit_greedy_l0()].
#'
#' @param formula An object of class `formula` describing the model.
#' @param data A data frame containing the variables in the model.
#' @param verbose Whether to show warning messages.
#' @param max_support_size,min_improvement See [inferCSN::fit_greedy_l0()].
#'
#' @return A list containing model, fit metrics, and nonzero coefficients.
#' @export
#' @examples
#' data("example_matrix", package = "inferCSN")
#' df <- as.data.frame(example_matrix)
#' fit_srm2(g1 ~ ., data = df)
fit_srm2 <- function(
  formula,
  data,
  verbose = TRUE,
  max_support_size = NULL,
  min_improvement = 1e-10
) {
  model_frame <- stats::model.frame(formula, data = data)
  model_mat <- stats::model.matrix(formula, data = model_frame)
  model_mat <- model_mat[, colnames(model_mat) != "(Intercept)", drop = FALSE]
  if (!ncol(model_mat)) {
    stop("At least one predictor is required.", call. = FALSE)
  }
  response <- stats::model.response(model_frame)

  result <- inferCSN::fit_greedy_l0(
    x = model_mat,
    y = response,
    max_support_size = max_support_size,
    min_improvement = min_improvement,
    verbose = verbose
  )

  nonzero <- result$coefficients$coefficient != 0
  result$metrics <- tibble::tibble(
    r_squared = result$metrics$r_squared
  )
  result$coefficients <- tibble::tibble(
    variable = result$coefficients$variable[nonzero],
    coefficient = result$coefficients$coefficient[nonzero]
  )

  result
}
