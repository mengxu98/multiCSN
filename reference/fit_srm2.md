# Fit a sparse regression model

Fits a sparse regression model using
[`inferCSN::fit_greedy_l0()`](https://mengxu98.github.io/inferCSN/reference/fit_greedy_l0.html).

## Usage

``` r
fit_srm2(
  formula,
  data,
  verbose = TRUE,
  max_support_size = NULL,
  min_improvement = 1e-10
)
```

## Arguments

- formula:

  An object of class *`formula`* with a symbolic description of the
  model to be fitted.

- data:

  A *`data.frame`* containing the variables in the model.

- verbose:

  If `TRUE`, show warning messages.

- max_support_size, min_improvement:

  See
  [inferCSN::fit_greedy_l0](https://mengxu98.github.io/inferCSN/reference/fit_greedy_l0.html).

## Value

A list containing two data frames:

- *`metrics`* - Goodness of fit measures

- *`coefficients`* - Fitted coefficients with sparse structure

## Examples

``` r
data("example_matrix", package = "inferCSN")
df <- as.data.frame(example_matrix)
fit_srm2(g1 ~ ., data = df)
#> $model
#> $model$algorithm
#> [1] "greedy_l0"
#> 
#> $model$support
#> [1] "g5" "g6"
#> 
#> $model$bic
#> [1] -111.8727
#> 
#> $model$rss
#> [1] 29.79499
#> 
#> $model$n_obs
#> [1] 100
#> 
#> 
#> $metrics
#> # A tibble: 1 × 1
#>   r_squared
#>       <dbl>
#> 1     0.699
#> 
#> $coefficients
#> # A tibble: 2 × 2
#>   variable coefficient
#>   <chr>          <dbl>
#> 1 g5            -0.492
#> 2 g6             0.601
#> 
```
