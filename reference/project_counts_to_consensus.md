# Project a peak count matrix onto consensus regions

Project a peak count matrix onto consensus regions

## Usage

``` r
project_counts_to_consensus(
  counts,
  peaks = rownames(counts),
  consensus,
  type = "within",
  require_one_to_one = TRUE
)
```

## Arguments

- counts:

  Peak-by-cell count matrix.

- peaks:

  Peak identifiers of `counts` (defaults to its row names).

- consensus:

  Consensus ranges from
  [`consensus_peak_set`](https://mengxu98.github.io/multiCSN/reference/consensus_peak_set.md).

- type, require_one_to_one:

  Passed to
  [`map_peaks_to_consensus`](https://mengxu98.github.io/multiCSN/reference/map_peaks_to_consensus.md).

## Value

Count matrix with one row per consensus region.
