# Collapse repeated gene annotations into peak-gene domains

Collapse repeated gene annotations into peak-gene domains

## Usage

``` r
collapse_peak_gene_domains(peak_annotation_matrix)
```

## Arguments

- peak_annotation_matrix:

  Sparse peak-by-annotation matrix, with annotation gene names as column
  names.

## Value

A gene-by-peak sparse matrix with duplicate gene annotations summed.
