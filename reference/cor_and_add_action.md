# Adds an extra column to the result of dynamic_shortest_path_multiple that predicts overall action based on correlation between "from" and "to"

Adds an extra column to the result of dynamic_shortest_path_multiple
that predicts overall action based on correlation between "from" and
"to"

## Usage

``` r
cor_and_add_action(spDF, matrix)
```

## Arguments

- spDF:

  the result of running dynamic_shortest_path_multiple

- matrix:

  corresponding genes-by-cells expression matrix

## Value

shortest paths dataframe with added action_by_corr column
