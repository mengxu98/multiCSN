# Scan for motifs in candidate regions

Scan for motifs in candidate regions

## Usage

``` r
find_motifs(object, ...)

# S4 method for class 'Seurat'
find_motifs(
  object,
  pfm,
  genome,
  motif_tfs = NULL,
  backend = c("signac", "motifmatchr"),
  verbose = TRUE,
  ...
)

# S4 method for class 'CSNObject'
find_motifs(object, ...)
```

## Arguments

- object:

  The input data, a csn object.

- ...:

  Arguments for other methods

- pfm:

  A PFMatrixList object with position weight matrices.

- genome:

  A BSgenome object.

- motif_tfs:

  A data frame matching motifs with TFs. First column motif, second TF.

- backend:

  Motif scanning backend. \`"signac"\` uses \`Signac::AddMotifs()\`;
  \`"motifmatchr"\` uses \`motifmatchr::matchMotifs()\` directly and
  wraps the result into a Signac \`Motif\` object.

- verbose:

  Display messages.

## Value

A CSNObject object with updated motif info.

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
} # }
```
