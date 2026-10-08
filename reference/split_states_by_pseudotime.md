# Splits data into states

Splits data into states by assigning cells to states

## Usage

``` r
split_states_by_pseudotime(dynamic_object, cuts, state_names = NULL)
```

## Arguments

- dynamic_object:

  result of running findDynGenes or compileDynGenes

- cuts:

  vector of pseudotime cutoffs

- state_names:

  names of resulting states, must have length of length(cuts)+1

## Value

updated dynamic_object with state column included in
dynamic_object\$cells
