# Fit model

Fits a sparse regression model through \[inferCSN::fit_greedy_l0()\].

## Usage

``` r
fit_model(formula, data, method = "greedy_l0", max_support_size = NULL, ...)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame containing the response and predictors.

- method:

  The sole supported method, \`"greedy_l0"\`.

- max_support_size:

  Maximum selected support size.

- ...:

  Additional arguments passed to \[fit_srm2()\].

## Value

A list containing the model, fit metrics, and selected coefficients.

## Examples

``` r
data("example_matrix", package = "inferCSN")
df <- as.data.frame(example_matrix)
fit_model(g1 ~ ., data = df)
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
