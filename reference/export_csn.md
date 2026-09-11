# Export CSN

Export network data from a CSN object with optional filtering

## Usage

``` r
export_csn(object, ...)

export_csn(object, ...)

# S4 method for class 'Seurat'
export_csn(
  object,
  active_network = NULL,
  celltypes = NULL,
  weight_cutoff = NULL,
  ...
)

# S4 method for class 'Network'
export_csn(object, weight_cutoff = NULL, ...)
```

## Arguments

- object:

  A CSNObject object

- ...:

  Additional arguments

- active_network:

  Character string specifying which network to export

- celltypes:

  Character vector of cell types to export

- weight_cutoff:

  Numeric threshold for filtering edges by absolute weight

## Value

A list of network data frames by cell type

## Examples

``` r
if (FALSE) { # \dontrun{
data(pbmcmultiome_sub, package = "scop")
data(motifs)
data(motif2tf)

object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
object <- initiate_object(
  object,
  group.by = "CellType",
  rna_assay = "RNA",
  peak_assay = "peaks",
  verbose = FALSE
)
genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
object <- find_motifs(
  object,
  pfm = motifs,
  motif_tfs = motif2tf,
  genome = genome,
  backend = "motifmatchr",
  verbose = FALSE
)
object <- inferCSN(object, cores = 2, verbose = FALSE)

networks <- export_csn(object)
head(networks[[1]])
} # }
```
