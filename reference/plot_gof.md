# Plot goodness of fit

Plot goodness of fit

Plot goodness-of-fit

Plot goodness-of-fit metrics.

## Usage

``` r
plot_gof(object, ...)

plot_gof(object, ...)

# S4 method for class 'Seurat'
plot_gof(object, network = NULL, celltypes = NULL, point_size = 0.5, ...)
```

## Arguments

- object:

  An object.

- ...:

  Additional arguments.

- network:

  Name of the network to use.

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

- point_size:

  Float indicating the point size.

## Value

A ggplot2 object.
