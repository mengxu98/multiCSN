# Prioritize network features using network or integrated evidence

Prioritize network features using network or integrated evidence

## Usage

``` r
prioritize_features(
  object,
  type = c("regulators", "targets"),
  strategy = c("network", "integrated"),
  network = DefaultNetwork(object),
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Feature type, `"regulators"` or `"targets"`.

- strategy:

  Prioritization strategy, `"network"` or `"integrated"`.

- network:

  Dynamic or static network name. Defaults to `DefaultNetwork(object)`.

- ...:

  Additional arguments forwarded to the underlying ranking or
  integrated-prioritization routine.

## Value

A data.frame of prioritized features.
