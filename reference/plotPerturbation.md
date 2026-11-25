# Plot perturbation

Plot perturbation

Visualize perturbation results

## Usage

``` r
plotPerturbation(
  object,
  point_size = 1,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5
)

plotPerturbation(
  object,
  point_size = 1,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5
)

# S4 method for class 'Network'
plotPerturbation(
  object,
  point_size = 1,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5
)

# S4 method for class 'Seurat'
plotPerturbation(
  object,
  point_size = 1,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5
)

# S4 method for class 'CSNObject'
plotPerturbation(
  object,
  point_size = 1,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5
)
```

## Arguments

- object:

  A Network or CSNObject object

- point_size:

  Numeric, size of points

- alpha:

  Point transparency

- show_trajectory:

  Logical, whether to show trajectory arrows

- arrow_size:

  Size of trajectory arrows

## Value

ggplot object
