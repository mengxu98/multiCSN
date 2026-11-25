# Normalize an ATAC assay to a Signac ChromatinAssay

Normalize an ATAC assay to a Signac ChromatinAssay

## Usage

``` r
normalize_chromatin_assay(
  object,
  peak_assay = NULL,
  new_assay = NULL,
  fragments = NULL,
  annotation = NULL,
  genome = NULL,
  sep = c("-", "-")
)
```

## Arguments

- object:

  A `Seurat` object.

- peak_assay:

  Name of the assay to normalize. Defaults to
  `Params(object)$peak_assay` or `"ATAC"` when present.

- new_assay:

  Optional assay name to write to. Defaults to overwriting `peak_assay`.

- fragments:

  Optional fragment object/path passed to
  [`Signac::CreateChromatinAssay()`](https://stuartlab.org/signac/reference/CreateChromatinAssay.html).

- annotation:

  Optional gene annotation `GRanges`. Defaults to the current assay
  annotation when available.

- genome:

  Optional genome metadata forwarded to
  [`Signac::CreateChromatinAssay()`](https://stuartlab.org/signac/reference/CreateChromatinAssay.html).

- sep:

  Character vector of two literal separators used in peak names.
  Explicit values control coordinate parsing. If omitted, the supported
  identifier forms are detected automatically.

## Value

A `Seurat` object with a `ChromatinAssay`.
