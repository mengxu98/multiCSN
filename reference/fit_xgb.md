# Fit a gradient boosting regression model with XGBoost

Fits a gradient boosted tree model using the XGBoost implementation.

## Usage

``` r
fit_xgb(
  formula,
  data,
  params = list(max_depth = 3, eta = 0.01, objective = "reg:squarederror"),
  nrounds = 1000,
  nthread = -1,
  ...
)
```

## Arguments

- formula:

  An object of class *`formula`* with a symbolic description of the
  model to be fitted.

- data:

  A *`data.frame`* containing the variables in the model.

- params:

  A list of XGBoost parameters including:

  - *`max_depth`* - Maximum tree depth

  - *`eta`* - Learning rate

  - *`objective`* - Loss function to optimize

- nrounds:

  Maximum number of boosting iterations.

- nthread:

  Number of parallel threads used (-1 for all cores).

- ...:

  Additional parameters passed to
  *[`xgboost::xgb.train`](https://rdrr.io/pkg/xgboost/man/xgb.train.html)*.

## Value

A list containing two data frames:

- *`metrics`* - Model performance metrics

- *`coefficients`* - Feature importance measures
