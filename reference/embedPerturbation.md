# Embed perturbation

Embed perturbation

Embed perturbation results

## Usage

``` r
embedPerturbation(
  object,
  reduction_method = "umap",
  dims = 2,
  scale = TRUE,
  cores = 1,
  seed = 1
)

embedPerturbation(
  object,
  reduction_method = "umap",
  dims = 2,
  scale = TRUE,
  cores = 1,
  seed = 1
)

# S4 method for class 'Seurat'
embedPerturbation(
  object,
  reduction_method = "umap",
  dims = 2,
  scale = TRUE,
  cores = 1,
  seed = 1
)

# S4 method for class 'Network'
embedPerturbation(
  object,
  reduction_method = "umap",
  dims = 2,
  scale = TRUE,
  cores = 1,
  seed = 1
)
```

## Arguments

- object:

  A Network or CSNObject object

- reduction_method:

  Character, dimensionality reduction method ("umap", "tsne", "pca")

- dims:

  Number of dimensions to compute

- scale:

  Whether to scale the data before reduction

- cores:

  Number of cores for parallel processing

- seed:

  Random seed for reproducibility

## Value

Updated object with embedded perturbation results
