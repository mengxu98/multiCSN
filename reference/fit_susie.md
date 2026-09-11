# Fit a SuSiE regression model

Fits a Sum of Single Effects (SuSiE) regression model using the susieR
package.

## Usage

``` r
fit_susie(formula, data, ...)
```

## Arguments

- formula:

  An object of class *`formula`* with a symbolic description of the
  model to be fitted.

- data:

  A *`data.frame`* containing the variables in the model.

- ...:

  Additional parameters passed to
  *[`susieR::susie`](https://rdrr.io/pkg/susieR/man/susie.html)*.

## Value

A list containing two data frames:

- *`metrics`* - Goodness of fit measures

- *`coefficients`* - Fitted coefficients
