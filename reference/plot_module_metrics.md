# Plot module metrics

Plot module metrics

Plot module metrics

Plot module metrics number of genes, number of peaks and number of TFs
per gene.

## Usage

``` r
plot_module_metrics(object, ...)

plot_module_metrics(object, ...)

# S4 method for class 'Seurat'
plot_module_metrics(object, network = NULL, celltypes = NULL, ...)
```

## Arguments

- object:

  CSNObject or Network

- ...:

  Additional arguments

- network:

  Name of the network to use.

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

## Value

A ggplot2 object.
