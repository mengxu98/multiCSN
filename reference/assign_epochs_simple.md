# Assigns genes to epochs just based on which mean is maximal

Assigns genes to epochs just based on which mean is maximal

## Usage

``` r
assign_epochs_simple(
  matrix,
  dynamic_object,
  num_epochs = 3,
  pThresh = 0.01,
  toScale = FALSE,
  key_word = "epoch"
)
```

## Arguments

- matrix:

  expression matrix

- dynamic_object:

  result of running findDynGenes

- num_epochs:

  num_epochs

- pThresh:

  pThresh

- toScale:

  toScale

- key_word:

  key_word

## Value

data.frame of dynamically expressed genes, cluster, peakTime, ordered by
peaktime
