# Plot adjacent-state TF rewiring heatmap

Plot adjacent-state TF rewiring heatmap

## Usage

``` r
plot_rewiring(
  object,
  network = DefaultNetwork(object),
  type = c("rewiring", "networks", "regulators", "targets"),
  plot_type = c("heatmap", "venn", "upset"),
  weight.by = "weight",
  group.by = NULL,
  palette = "Chinese",
  palcolor = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- type:

  Plot type. `"rewiring"` shows TF-by-transition rewiring; `"networks"`,
  `"regulators"`, and `"targets"` show HARNexus-style pairwise state
  similarity based on edges, regulator sets, or target-gene sets.

- plot_type:

  Plot form. Default is `"heatmap"`; `"venn"` and `"upset"` are
  supported for similarity plots.

- weight.by:

  Weight column used when exporting networks.

- group.by:

  Reserved for future grouping support in similarity plots.

- palette:

  Palette name passed to
  [`thisplot::palette_colors()`](https://mengxu98.github.io/thisplot/reference/palette_colors.html)
  via `.multicsn_palette_colors()`.

- palcolor:

  Optional named color vector overriding state/group colors.

## Value

A drawable plot object. When ComplexHeatmap is available, the function
returns a grid grob captured from the drawn heatmap; otherwise it
returns a `ggplot` object.
