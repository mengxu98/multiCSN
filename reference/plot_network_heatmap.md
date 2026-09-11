# Plot a network heatmap

Plot a network heatmap

## Usage

``` r
plot_network_heatmap(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE,
  show_names = FALSE,
  show_names_position = c("outer", "all"),
  rect_color = NA,
  heatmap_size_lock = TRUE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  ncol = NULL,
  nrow = NULL,
  performance_metrics = "none",
  performance_ground_truth = NULL,
  truth_cell_border = TRUE,
  truth_cell_border_legend = TRUE,
  truth_match_border_color = "#2E7D32",
  truth_mismatch_border_color = "#F2C94C",
  truth_cell_border_lwd = 1.2,
  border_color = "gray",
  heatmap_palette = NULL,
  heatmap_palcolor = NULL,
  heatmap_color = c("#1966ad", "white", "#bb141a"),
  row_anno_palette = "Set1",
  row_anno_palcolor = NULL,
  col_anno_palette = "Set2",
  col_anno_palcolor = NULL,
  row_anno = NULL,
  column_anno = NULL,
  anno_width = 1,
  anno_height = 1,
  row_anno_type = c("boxplot", "barplot", "histogram", "density", "lines", "points",
    "horizon"),
  column_anno_type = c("boxplot", "barplot", "histogram", "density", "lines", "points"),
  legend_name = NULL,
  row_title = "Regulators"
)
```

## Arguments

- network_table:

  Network edge table.

- regulators, targets:

  Nodes to include.

- switch_matrix:

  Convert edge tables to matrices.

- show_names, heatmap_size_lock:

  Logical display controls.

- show_names_position:

  Name placement for multiple heatmaps.

- heatmap_size, heatmap_height, heatmap_width:

  Heatmap dimensions.

- heatmap_title:

  Heatmap title.

- ncol, nrow:

  Layout dimensions for multiple heatmaps.

- performance_metrics:

  Metrics appended to titles.

- performance_ground_truth:

  Ground-truth network for metrics.

- truth_cell_border, truth_cell_border_legend:

  Ground-truth border controls.

- truth_match_border_color, truth_mismatch_border_color:

  Border colors.

- truth_cell_border_lwd:

  Border width.

- border_color, rect_color:

  Border and cell colors.

- heatmap_palcolor, heatmap_palette:

  Heatmap colors.

- heatmap_color:

  Alias for \`heatmap_palcolor\`. Used when \`heatmap_palcolor\` is
  \`NULL\`.

- row_anno_palette, row_anno_palcolor:

  Row annotation colors.

- col_anno_palette, col_anno_palcolor:

  Column annotation colors.

- row_anno, column_anno:

  Enable annotations. Default is on for a single network and off for a
  list of networks.

- anno_width, anno_height:

  Annotation dimensions.

- row_anno_type, column_anno_type:

  Annotation plot types.

- legend_name, row_title:

  Legend and row titles.

## Value

A heatmap object.

## Examples

``` r
data(example_matrix, package = "inferCSN")
data("example_ground_truth", package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix, verbose = FALSE)

# both heatmaps must show the same genes, so align them on the shared set
regulators <- intersect(
  example_ground_truth$regulator,
  network_table$regulator
)
targets <- intersect(
  example_ground_truth$target,
  network_table$target
)

p1 <- plot_network_heatmap(
  example_ground_truth[, 1:3],
  regulators = regulators,
  targets = targets,
  heatmap_title = "Ground truth",
  legend_name = "Ground truth"
)
p2 <- plot_network_heatmap(
  network_table,
  regulators = regulators,
  targets = targets,
  heatmap_title = "inferCSN",
  legend_name = "inferCSN"
)
ComplexHeatmap::draw(p1 + p2)
```
