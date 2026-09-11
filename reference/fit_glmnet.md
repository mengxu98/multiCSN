# Fit regularized generalized linear model

Fits a regularized generalized linear model using elastic net penalties
via glmnet.

## Usage

``` r
fit_glmnet(formula, data, family = gaussian, alpha = 0.5, ...)
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

  The elasticnet mixing parameter, between 0 and 1:

  - *`alpha = 1`*: Lasso penalty

  - *`alpha = 0`*: Ridge penalty

  - *`0 < alpha < 1`*: Elastic net penalty

- ...:

  Additional parameters passed to
  *[`glmnet::glmnet`](https://rdrr.io/pkg/glmnet/man/glmnet.html)*.

## Value

A list containing two data frames:

- *`metrics`* - Goodness of fit measures including lambda and R-squared

- *`coefficients`* - Regularized coefficient estimates
