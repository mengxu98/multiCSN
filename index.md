# **multiCSN**

## **Introduction**

[multiCSN](https://mengxu98.github.io/multiCSN/) is an R package for
**infer**ring **C**ell-**S**pecific gene regulatory **N**etwork from
single-cell omics data.

![multiCSN workflow
diagram](https://raw.githubusercontent.com/mengxu98/figures/main/multiCSN/multiCSN.svg#gh-light-mode-only)

![multiCSN workflow
diagram](https://raw.githubusercontent.com/mengxu98/figures/main/multiCSN/multiCSN-dark.svg#gh-dark-mode-only)

## **Installation**

Install the development version from
[GitHub](https://github.com/mengxu98/multiCSN) use
[pak](https://github.com/r-lib/pak):

``` r

if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/multiCSN")
```

## **Usage**

### **Layered multi-omic networks**

For a Seurat object with normalized RNA and TF-IDF chromatin assays,
build the candidate catalogue once and fit one network per cell group.
The bundled motifs are used by default.

``` r

library(multiCSN)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

object <- prepare_candidates(object, group.by = "cell_type", genome = genome)

fit <- fit_layered_network(
  object,
  cell_group = "Erythroblast",
  preset = "state",
  output_dir = "networks/erythroblast"
)
```

`preset = "state"` fits cell states that share one normalization;
`preset = "timepoint"` re-normalizes each group and requires target
genes and TFs to be detected in at least 20 cells. Arguments you pass
explicitly override the preset. With `output_dir` the network tables are
written only after the fit completes, and an existing directory is never
overwritten.

Each assay needs one normalized `data` layer. Join split RNA layers
before candidate preparation with
`object <- SeuratObject::JoinLayers(object, assay = "RNA", layers = "data", new = "data")`.
Per-group normalization also needs joined `counts` layers. Candidate
preparation rebuilds the catalogue, so prepare once and then fit each
group from that object.

For an in-memory result, omit `output_dir`. The default `store = FALSE`
keeps the input object unchanged and returns the complete layer and path
tables, accounting, effective settings and selected cell and feature
IDs. It does not include a Seurat object. Set `store = TRUE` to also
return an object with stored networks.

Read, check and compare the saved networks:

``` r

validate_network_endpoints("networks/erythroblast")
networks <- read_network_endpoints("networks/erythroblast")

early <- read_network_endpoints("networks/proerythroblast", "TF-gene", direction = TRUE, simplify = TRUE)
late <- read_network_endpoints("networks/erythroblast", "TF-gene", direction = TRUE, simplify = TRUE)
classify_edge_transitions(early, late, from = "Proerythroblast", to = "Erythroblast")
```

### **Network evaluation**

[`evaluate_network()`](https://mengxu98.github.io/multiCSN/reference/evaluate_network.md)
scores a predicted edge table against a ground truth with BEELINE-style
metrics, and
[`ranking_metrics()`](https://mengxu98.github.io/multiCSN/reference/ranking_metrics.md),
[`tie_aware_auroc()`](https://mengxu98.github.io/multiCSN/reference/ranking_metrics.md),
[`tie_aware_average_precision()`](https://mengxu98.github.io/multiCSN/reference/ranking_metrics.md)
and
[`expected_precision_at_k()`](https://mengxu98.github.io/multiCSN/reference/ranking_metrics.md)
give tie-aware ranking metrics for any score vector.

``` r

evaluate_network(prediction, ground_truth)$AUPRC
```

### **Single-omic networks**

``` r

library(multiCSN)
data("example_matrix", package = "inferCSN")

network <- inferCSN(example_matrix, verbose = FALSE)

pseudotime <- setNames(seq(0, 1, length.out = nrow(example_matrix)), rownames(example_matrix))
ordered_network <- inferCSN(example_matrix, pseudotime = pseudotime, verbose = FALSE)
```

The RNA matrix has cells in rows and genes in columns. This public
method uses the independent `inferCSN` engine. Use the layered workflow
above for paired RNA and chromatin assays. See
[`?fit_layered_network`](https://mengxu98.github.io/multiCSN/reference/fit_layered_network.md)
for the meanings of layer coefficients, deletion evidence, ordinal
weights and mediated paths.

See the [function
reference](https://mengxu98.github.io/multiCSN/reference/index.html) for
the full interface.
