# The CSNObject class

The CSNObject object is an extended `Seurat` object for the storage and
analysis of celltype-specific gene regulatory networks.

## Slots

- `data`:

  Seurat object.

- `metadata`:

  A list with metadata about the object.

- `regions`:

  A
  [`Regions`](https://mengxu98.github.io/multiCSN/reference/Regions-class.md)
  object containing information about the genomic regions included in
  the network.

- `networks`:

  A
  [`Network`](https://mengxu98.github.io/multiCSN/reference/Network-class.md)
  object containing the inferred regulatory network and information
  about the model fit.

- `params`:

  A list storing parameters for network inference.

- `active_network`:

  A string indicating the active network.
