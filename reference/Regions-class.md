# The Regions class

The Regions object stores the genomic regions that are considered by the
model. It stores their genomic positions, how they map to the peaks in
the Seurat object and motif matches.

## Slots

- `motifs`:

  A `Motifs` object with matches of TF motifs.

- `motifs2tfs`:

  tfs.

- `ranges`:

  A `GenomicRanges` object.

- `peaks`:

  A numeric vector with peak indices for each region.
