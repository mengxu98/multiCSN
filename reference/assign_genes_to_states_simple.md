# Assigns genes to states just based on which mean is maximal

Assigns genes to states just based on which mean is maximal

## Usage

``` r
assign_genes_to_states_simple(
  matrix,
  dynamic_object,
  num_states = 3,
  pThresh = 0.01,
  toScale = FALSE,
  key_word = "state_"
)
```

## Arguments

- matrix:

  expression matrix

- dynamic_object:

  result of running findDynGenes

- num_states:

  num_states

- pThresh:

  pThresh

- toScale:

  toScale

- key_word:

  key_word

## Value

data.frame of dynamically expressed genes, cluster, peakTime, ordered by
peaktime
