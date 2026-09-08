#' @title Plot coefficients
#'
#' @md
#' @param data Input data with `regulator` and `weight` columns.
#' @param style Plotting style: `"binary"` or `"continuous"`. Both use
#' [thisplot::StatPlot] signed value bars; `"binary"` uses the two endpoint
#' colors and `"continuous"` includes a neutral midpoint color.
#' @param positive_color Color for positive weights.
#' Default is `"#3d67a2"`.
#' @param negative_color Color for negative weights.
#' Default is `"#c82926"`.
#' @param neutral_color Color for weights near zero (used in `"continuous"`).
#' Default is `"#cccccc"`.
#' @param bar_width Width of the bars.
#' Default is `0.7`.
#' @param text_size Size of the text for weight values.
#' Default is `3`.
#' @param show_values Whether to show weight values on bars.
#' Default is `TRUE`.
#' @param ... Additional arguments passed to [thisplot::StatPlot].
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_matrix, package = "inferCSN")
#' network_table <- inferCSN::inferCSN(example_matrix, targets = "g1")
#' plot_coefficient(network_table)
#' plot_coefficient(network_table, style = "binary")
#' }
plot_coefficient <- function(
  data,
  style = c("continuous", "binary"),
  positive_color = "#3d67a2",
  negative_color = "#c82926",
  neutral_color = "#cccccc",
  bar_width = 0.7,
  text_size = 3,
  show_values = TRUE,
  ...
) {
  if (!requireNamespace("thisplot", quietly = TRUE)) {
    thisutils::log_message(
      "{.fn plot_coefficient} requires {.pkg thisplot}",
      message_type = "error"
    )
  }
  style <- match.arg(style)
  data <- as.data.frame(data)
  if (!all(c("regulator", "weight") %in% colnames(data))) {
    thisutils::log_message(
      "{.arg data} must contain {.field regulator} and {.field weight} columns",
      message_type = "error"
    )
  }
  palcolor <- if (identical(style, "binary")) {
    c(negative_color, positive_color)
  } else {
    c(negative_color, neutral_color, positive_color)
  }
  thisplot::StatPlot(
    data,
    stat.by = "regulator",
    value.by = "weight",
    stat_type = "value",
    plot_type = "bar",
    flip = TRUE,
    palette = "RdBu",
    palcolor = palcolor,
    bar_width = bar_width,
    label = show_values,
    label.size = text_size,
    xlab = "Regulator",
    ylab = "Weight",
    value_legend_title = "Weight",
    grid_major = FALSE,
    ...
  )
}

#' @title Plot coefficients for multiple targets
#'
#' @param data Input data.
#' @param targets Targets to plot.
#' Default is `NULL`.
#' @param nrow Number of rows for the plot.
#' Default is `NULL`.
#' @param ... Other arguments passed to [plot_coefficient].
#'
#' @return A patchwork of ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_matrix, package = "inferCSN")
#' network_table <- inferCSN::inferCSN(
#'   example_matrix,
#'   targets = c("g1", "g2", "g3")
#' )
#' plot_coefficients(network_table, show_values = FALSE)
#' plot_coefficients(network_table, targets = "g1")
#' }
plot_coefficients <- function(
  data,
  targets = NULL,
  nrow = NULL,
  ...
) {
  if (is.null(targets)) {
    targets <- unique(data$target)
  }
  p_list <- lapply(targets, function(target) {
    plot_coefficient(data[data$target == target, ], ...)
  })
  if (!is.null(nrow)) {
    patchwork::wrap_plots(p_list, nrow = nrow)
  } else {
    patchwork::wrap_plots(p_list)
  }
}
