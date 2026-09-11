# Get pseudotime information

Get pseudotime information

## Usage

``` r
get_pseudotime(object, ...)

# Default S3 method
get_pseudotime(
  object,
  meta_data = NULL,
  embeddings = NULL,
  cluster_column = "cluster",
  cores = 1,
  seed = 1,
  start_cluster = NULL,
  end_cluster = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'Seurat'
get_pseudotime(
  object,
  assay = "RNA",
  layer = "data",
  cluster_column = "cluster",
  reduction = "umap",
  start_cluster = NULL,
  end_cluster = NULL,
  ...
)
```

## Arguments

- object:

  The input data (matrix, Seurat, etc.)

- ...:

  Arguments for other methods

- meta_data:

  Input meta data.

- embeddings:

  Embeddings information

- cluster_column:

  The column used for slingshot.

- cores:

  Number of cores to use for parallel processing. Default is 1.

- seed:

  Random seed for reproducibility. Default is 1.

- start_cluster:

  The start cluster.

- end_cluster:

  The end cluster.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

- assay:

  The assay used for slingshot.

- layer:

  The layer used for slingshot.

- reduction:

  The reduction used for slingshot.

## Value

Pseudotime information corresponding to cells
