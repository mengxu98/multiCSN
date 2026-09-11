# Calculate gene rank

Calculate gene rank

## Usage

``` r
calculate_gene_rank(object, ...)

# S4 method for class 'Network'
calculate_gene_rank(
  object,
  regulators = NULL,
  targets = NULL,
  directed = FALSE,
  method = c("page_rank", "degree_distribution"),
  ...
)

# S4 method for class 'Seurat'
calculate_gene_rank(
  object,
  network = NULL,
  celltypes = NULL,
  regulators = NULL,
  targets = NULL,
  directed = FALSE,
  method = c("page_rank", "degree_distribution"),
  ...
)

# S4 method for class 'CSNObject'
calculate_gene_rank(object, ...)

# S4 method for class 'data.frame'
calculate_gene_rank(
  object,
  regulators = NULL,
  targets = NULL,
  directed = FALSE,
  method = c("page_rank", "degree_distribution"),
  ...
)
```

## Arguments

- object:

  Network object

- ...:

  Additional arguments

- regulators:

  Character vector, regulators to include

- targets:

  Character vector, targets to include

- directed:

  Logical, whether the network is directed

- method:

  Character, ranking method: "page_rank" or "degree_distribution"

- network:

  Name of the network to use.

- celltypes:

  Character vector specifying which states/cell types to include; NULL
  for all.

## Value

Data frame with gene ranks

The Network object with gene ranks stored in `params$gene_rank`.

The Seurat object with gene ranks stored in each Network's
`params$gene_rank`.
