# Rank signed network edges by deletion evidence

Rank signed network edges by deletion evidence

## Usage

``` r
network_ordinal_weight(coefficient, deletion_delta_bic)
```

## Arguments

- coefficient:

  Signed coefficients.

- deletion_delta_bic:

  Nonnegative deletion delta-BIC values, one per coefficient.

## Value

Signed ordinal weights; invalid or zero coefficients receive zero.

## Details

Descending evidence is grouped against each group's maximum using
absolute difference \<= 1e-12 \* (1 + abs(maximum)). This numerical rule
does not change fitted support, coefficients or raw deletion evidence.
