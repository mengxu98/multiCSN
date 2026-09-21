# Map sample peaks into a consensus peak set

Map sample peaks into a consensus peak set

## Usage

``` r
map_peaks_to_consensus(
  peaks,
  consensus,
  type = "within",
  require_one_to_one = TRUE
)
```

## Arguments

- peaks:

  A `GRanges` or `"chr:start-end"` character vector.

- consensus:

  Consensus ranges from
  [`consensus_peak_set`](https://mengxu98.github.io/multiCSN/reference/consensus_peak_set.md).

- type:

  Overlap type forwarded to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).
  The default requires each peak to fall inside one consensus region.

- require_one_to_one:

  Fail when a peak maps to zero or several consensus regions; that is
  the condition under which count projection is well defined.

## Value

Integer vector of consensus indices, one per input peak.
