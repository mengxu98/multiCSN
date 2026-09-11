# Plot perturbation

Plot perturbation

Visualize perturbation results

## Usage

``` r
plotPerturbation(
  object,
  label_points = FALSE,
  point_size = 1,
  max_overlaps = 10,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5,
  grid_n = 25
)

plotPerturbation(
  object,
  label_points = FALSE,
  point_size = 1,
  max_overlaps = 10,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5,
  grid_n = 25
)

# S4 method for class 'Network'
plotPerturbation(
  object,
  label_points = FALSE,
  point_size = 1,
  max_overlaps = 10,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5,
  grid_n = 25
)

# S4 method for class 'Seurat'
plotPerturbation(
  object,
  label_points = FALSE,
  point_size = 1,
  max_overlaps = 10,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5,
  grid_n = 25
)

# S4 method for class 'CSNObject'
plotPerturbation(
  object,
  label_points = FALSE,
  point_size = 1,
  max_overlaps = 10,
  alpha = 0.6,
  show_trajectory = TRUE,
  arrow_size = 0.5,
  grid_n = 25
)
```

## Arguments

- object:

  A Network or CSNObject object

- label_points:

  Logical, whether to add point labels

- point_size:

  Numeric, size of points

- max_overlaps:

  Maximum number of overlapping labels to allow

- alpha:

  Point transparency

- show_trajectory:

  Logical, whether to show trajectory arrows

- arrow_size:

  Size of trajectory arrows

- grid_n:

  Number of grid points for trajectory arrows

## Value

ggplot object
