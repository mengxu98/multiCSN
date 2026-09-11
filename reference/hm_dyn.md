# plots results of findDynGenes

plots results of findDynGenes

## Usage

``` r
hm_dyn(
  matrix,
  dynRes,
  cluster = "cluster",
  topX = 25,
  cRow = FALSE,
  cCol = FALSE,
  limits = c(0, 10),
  toScale = FALSE,
  fontsize_row = 4,
  geneAnn = FALSE,
  anno_colors = NULL,
  show_rownames = TRUE,
  filename = NA,
  width = NA,
  height = NA
)
```

## Arguments

- matrix:

  expression matrix

- dynRes:

  result of running findDynGenes

- cluster:

  cluster

- topX:

  topX

- cRow:

  cRow

- cCol:

  cCol

- limits:

  limits

- toScale:

  toScale

- fontsize_row:

  fontsize_row

- geneAnn:

  geneAnn

- anno_colors:

  anno_colors

- show_rownames:

  show_rownames

- filename:

  filename

- width:

  width

- height:

  height

## Value

heatmap list
