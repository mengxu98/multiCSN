# Prepare a multiome object for layered network inference

Builds the candidate catalogue used by
[`fit_layered_network`](https://mengxu98.github.io/multiCSN/reference/fit_layered_network.md):
it initializes the multiCSN state without marker filtering (every RNA
gene is a candidate target), aligns the gene annotation with the peak
chromosome names and scans the selected peaks for motifs.

## Usage

``` r
prepare_candidates(
  object,
  group.by = "aggregate",
  genome,
  rna_assay = "RNA",
  peak_assay = "ATAC",
  targets = NULL,
  pfm = NULL,
  motif_tfs = NULL,
  backend = c("motifmatchr", "signac"),
  exclude_exons = TRUE,
  align_annotation = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with normalized RNA and TF-IDF chromatin assays.

- group.by:

  Metadata column with the cell groups to fit, or `"aggregate"` for one
  group with all cells.

- genome:

  A BSgenome object matching the peak coordinates.

- rna_assay, peak_assay:

  Assay names.

- targets:

  Candidate target genes; defaults to every RNA gene.

- pfm, motif_tfs:

  Motif matrices and motif-to-TF map; default to the bundled
  [`motifs`](https://mengxu98.github.io/multiCSN/reference/motifs.md)
  and
  [`motif2tf`](https://mengxu98.github.io/multiCSN/reference/motif2tf.md).

- backend:

  Motif scanning backend, see
  [`find_motifs`](https://mengxu98.github.io/multiCSN/reference/find_motifs.md).

- exclude_exons:

  Remove exonic peaks from the candidates.

- align_annotation:

  Add peak chromosomes missing from the gene annotation to its sequence
  levels.

- verbose:

  Print progress messages.

## Value

The object with candidate genes, peaks, transcription factors and motif
matches.

## Details

Both assays require a single normalized layer named `data`. Join split
layers with
[`SeuratObject::JoinLayers()`](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
before preparing candidates. RNA values must use the normalization
intended for the fit; chromatin values must be TF-IDF normalized. This
function does not normalize either assay.

Every call rebuilds the multiCSN candidate catalogue and clears networks
from the previous catalogue on the returned object. It sets active
identities from `group.by`, selects the RNA default assay and, when
requested, extends annotation sequence levels. Other assays and metadata
are retained. Keep the returned object and reuse it for fits of the same
catalogue.

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
object <- prepare_candidates(object, group.by = "cell_type", genome = genome)
fit <- fit_layered_network(object, cell_group = "HSC", preset = "state")
} # }
```
