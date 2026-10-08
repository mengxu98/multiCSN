# Divides grnDF into states, filters interactions between genes not in same or consecutive states

Divides grnDF into states, filters interactions between genes not in
same or consecutive states

## Usage

``` r
split_network_by_states(grnDF, states, state_network = NULL)
```

## Arguments

- grnDF:

  result of GRN reconstruction

- states:

  result of running assign_genes_to_states

- state_network:

  dataframe outlining higher level state connectivity (i.e. state
  transition network). If NULL, will assume states is ordered linear
  trajectory

## Value

list of GRNs across states and transitions
