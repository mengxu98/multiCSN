# Find TF modules in regulatory network

Find TF modules in regulatory network

## Usage

``` r
find_modules(object, ...)

# S4 method for class 'Network'
find_modules(
  object,
  p_thresh = 0.05,
  rsq_thresh = 0.1,
  nvar_thresh = 10,
  min_genes_per_module = 5,
  xgb_method = c("tf", "target"),
  xgb_top = 50,
  verbose = TRUE,
  ...
)

# S4 method for class 'Seurat'
find_modules(
  object,
  network = NULL,
  p_thresh = 0.05,
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

- p_thresh:

  Float indicating the significance threshold on the adjusted p-value.

- rsq_thresh:

  Float indicating the \\R^2\\ threshold on the adjusted p-value.

- nvar_thresh:

  Integer indicating the minimum number of variables in the model.

- min_genes_per_module:

  Integer indicating the minimum number of genes in a module.

- xgb_method:

  Method to get modules from xgb models `tf` - Choose top targets for
  each TF. `target` - Choose top TFs for each target gene.

- xgb_top:

  Interger indicating how many top targets/TFs to return.

- verbose:

  Print messages.

- network:

  Name of the network to use.

## Value

A Network object.

A CSNObject object
