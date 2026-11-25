# Signed ordinal edge weights

Converts signed coefficients into signed ranks in (0, 1\]: magnitudes
are ranked from strongest to weakest, magnitudes within `tolerance`
share one rank, and the sign of the coefficient is kept. Zero and
non-finite coefficients receive weight 0. Use
[`inferCSN::deletion_rank_weight()`](https://mengxu98.github.io/inferCSN/reference/deletion_rank_weight.html)
when the rank should follow deletion evidence rather than the
coefficient size.

## Usage

``` r
signed_ordinal_weight(values, tolerance = 1e-12)
```

## Arguments

- values:

  Signed coefficients.

- tolerance:

  Magnitudes closer than this are tied.

## Value

Numeric vector of signed ordinal weights.

## Examples

``` r
signed_ordinal_weight(c(-2, 1, 0, 2, 0.5))
#> [1] -0.8333333  0.5000000  0.0000000  0.8333333  0.1666667
```
