# Predict perturbation

Predict perturbation

## Usage

``` r
predictPerturbation(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 5,
  scale_factor = 0.1,
  verbose = TRUE
)

predictPerturbation(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 5,
  scale_factor = 0.1,
  verbose = TRUE
)

# S4 method for class 'Network'
predictPerturbation(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 3,
  scale_factor = 0.1,
  verbose = TRUE
)

# S4 method for class 'Seurat'
predictPerturbation(
  object,
  perturb_tfs,
  use_weight = FALSE,
  n_iter = 3,
  scale_factor = 0.1,
  verbose = TRUE
)

# S4 method for class 'CSNObject'
predictPerturbation(object, perturb_tfs, use_weight = FALSE)
```

## Arguments

- object:

  A Network or CSNObject object

- perturb_tfs:

  Named vector of TF perturbation states

- use_weight:

  Logical, whether to use network weights instead of coefficients

- n_iter:

  Number of iterations for perturbation simulation

- scale_factor:

  Scale factor for perturbation effect (0-1)

- verbose:

  Whether to show progress messages

## Value

Updated object with perturbation predictions
