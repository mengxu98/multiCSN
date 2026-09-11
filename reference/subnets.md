# Function to assign nodes to communities via Louvain clustering

Function to assign nodes to communities via Louvain clustering

## Usage

``` r
subnets(df, tfs, tfonly = TRUE)
```

## Arguments

- df:

  dataframe containing a static network

- tfs:

  TFs

- tfonly:

  if TRUE will limit network to TFs only

## Value

dataframe containing community assignments for each gene
