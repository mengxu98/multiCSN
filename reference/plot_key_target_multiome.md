# Plot integrated key-target ranking with ATAC support

Plot integrated key-target ranking with ATAC support

## Usage

``` r
plot_key_target_multiome(target_scores, top_n = 10)
```

## Arguments

- target_scores:

  Output of
  `prioritize_features(type = "targets", strategy = "integrated")`.

- top_n:

  Number of targets shown per TF.

## Value

A `ggplot` object.
