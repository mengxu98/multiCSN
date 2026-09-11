# Plot integrated multiome key-TF scores

Plot integrated multiome key-TF scores

## Usage

``` r
plot_key_tf_multiome(tf_scores, top_n = 12)
```

## Arguments

- tf_scores:

  Output of
  `prioritize_features(type = "regulators", strategy = "integrated")`.

- top_n:

  Number of TFs shown.

## Value

A patchwork/ggplot object.
