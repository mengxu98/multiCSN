# Splits data into states manually

Splits data into states given group assignment

## Usage

``` r
split_states_by_group(dynamic_object, assignment)
```

## Arguments

- dynamic_object:

  result of running findDynGenes or compileDynGenes

- assignment:

  a list of vectors where names(assignment) are state names, and vectors
  contain groups belonging to corresponding state

## Value

updated dynamic_object with state column included in
dynamic_object\$cells
