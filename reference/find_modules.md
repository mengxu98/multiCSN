# Find TF modules in regulatory network

Find TF modules in regulatory network

## Usage

``` r
find_modules(object, ...)

# S4 method for class 'Network'
find_modules(
  object,
  rsq_thresh = 0.1,
  nvar_thresh = 10,
  min_genes_per_module = 5,
  verbose = TRUE,
  ...
)

# S4 method for class 'Seurat'
find_modules(
  object,
  network = NULL,
  rsq_thresh = 0.1,
  nvar_thresh = 10,
  min_genes_per_module = 5,
  verbose = TRUE,
  ...
)

# S4 method for class 'CSNObject'
find_modules(object, ...)
```

## Arguments

- object:

  An object.

- ...:

  Arguments for other methods

- rsq_thresh:

  Float indicating the \\R^2\\ threshold on the adjusted p-value.

- nvar_thresh:

  Integer indicating the minimum number of variables in the model.

- min_genes_per_module:

  Integer indicating the minimum number of genes in a module.

- verbose:

  Print messages.

- network:

  Name of the network to use.

## Value

A Network object.

A CSNObject object
