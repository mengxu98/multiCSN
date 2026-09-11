# Fit generalized linear model

Fits a generalized linear model using standard maximum likelihood
estimation.

## Usage

``` r
fit_glm(formula, data, family = gaussian, ...)
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

- ...:

  Additional parameters passed to
  *[`stats::glm`](https://rdrr.io/r/stats/glm.html)*.

## Value

A list containing two data frames:

- *`metrics`* - Goodness of fit measures including R-squared

- *`coefficients`* - Fitted coefficients with standard errors and
  p-values
