# Function to return shortest path from 1 regulator to 1 target in a dynamic network

Function to return shortest path from 1 regulator to 1 target in a
dynamic network

## Usage

``` r
dynamic_shortest_path(
  network_table,
  regulator,
  target,
  weight_column = "weight",
  compare_to_average = FALSE
)
```

## Arguments

- network_table:

  a dyanmic network

- regulator:

  the starting regulator

- target:

  the end regulator

- weight_column:

  column name in network_table with edge weights that will be converted
  to distances

- compare_to_average:

  if TRUE will compute normalized against average path length

## Value

shortest path, distance, normalized distance, and action
