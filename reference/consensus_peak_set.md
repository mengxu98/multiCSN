# Build a consensus peak set across samples

Merges the peak sets of several samples into one non-overlapping
universe. Overlapping and adjacent intervals are collapsed with
[`GenomicRanges::reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html),
optionally restricted to a set of sequence levels and annotated with a
`Seqinfo` (for example the UCSC hg38 chromosome lengths), which is what
downstream count projection needs.

## Usage

``` r
consensus_peak_set(peaks, seqlevels = NULL, seqinfo = NULL)
```

## Arguments

- peaks:

  A `GRanges`, a character vector of `"chr:start-end"` regions, or a
  list of either (one element per sample).

- seqlevels:

  Optional character vector of sequence levels to keep.

- seqinfo:

  Optional `Seqinfo` used to annotate the kept sequence levels.

## Value

A named `GRanges` whose names are `"chr:start-end"`.
