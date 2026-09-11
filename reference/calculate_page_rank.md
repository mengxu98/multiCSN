# Calculate PageRank

Calculate PageRank

## Usage

``` r
calculate_page_rank(object, ...)

# S4 method for class 'Network'
calculate_page_rank(
  object,
  regulators = NULL,
  targets = NULL,
  directed = FALSE
)

# S4 method for class 'data.frame'
calculate_page_rank(
  object,
  regulators = NULL,
  targets = NULL,
  directed = FALSE
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

## Value

Data frame with PageRank scores
