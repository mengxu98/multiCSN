# Cross-validation for regularized generalized linear models

Performs cross-validation to select optimal regularization parameters
for GLMnet models.

## Usage

``` r
fit_cvglmnet(formula, data, family = gaussian, alpha = 0.5, ...)
```

## Arguments

- formula:

  An object of class *`formula`* with a symbolic description of the
  model to be fitted.

- data:

  A *`data.frame`* containing the variables in the model.

- family:

  A description of the error distribution and link function. See
  *[`stats::family`](https://rdrr.io/r/stats/family.html)* for details.

- alpha:

  The elasticnet mixing parameter, between 0 and 1. See
  *[`glmnet::glmnet`](https://rdrr.io/pkg/glmnet/man/glmnet.html)* for
  details.

- ...:

  Additional parameters passed to
  *[`glmnet::cv.glmnet`](https://rdrr.io/pkg/glmnet/man/cv.glmnet.html)*.

## Value

A list containing two data frames:

- *`metrics`* - Cross-validated performance metrics

- *`coefficients`* - Coefficients selected by cross-validation
