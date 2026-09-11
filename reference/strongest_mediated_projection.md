# Project strongest sign-consistent mediated paths

Project strongest sign-consistent mediated paths

## Usage

``` r
strongest_mediated_projection(chains)
```

## Arguments

- chains:

  Output of factorized_chain_support.

## Value

One strongest path per regulator-target pair; inconsistent direction
ties are excluded.

## Details

Paths within 1e-12 \* (1 + abs(maximum evidence)) of the maximum are
numerical ties. Opposite-sign ties are excluded. The representative
region is lexicographically first; pair evidence remains the maximum.
