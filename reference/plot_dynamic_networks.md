# Plot dynamic networks

Plot dynamic networks

## Usage

``` r
plot_dynamic_networks(object, ...)

# S4 method for class 'Seurat'
plot_dynamic_networks(
  object,
  network = NULL,
  celltypes = NULL,
  weight_cutoff = NULL,
  celltypes_order = NULL,
  ntop = 10,
  filter_by_ntop = TRUE,
  ...
)

# S4 method for class 'CSNObject'
plot_dynamic_networks(object, ...)

# S4 method for class 'data.frame'
plot_dynamic_networks(
  object,
  celltypes_order = NULL,
  ntop = 10,
  filter_by_ntop = TRUE,
  title = NULL,
  theme_type = "theme_void",
  plot_type = "ggplot",
  layout = "fruchtermanreingold",
  nrow = 2,
  ncol = NULL,
  byrow = TRUE,
  combine = TRUE,
  figure_save = FALSE,
  figure_name = NULL,
  figure_width = 6,
  figure_height = 6,
  seed = 1
)
```

## Arguments

- object:

  Data.frame with columns regulator, target, weight, celltype.

- ...:

  Arguments passed to methods.

- network:

  Name of the network; default uses `DefaultNetwork(object)`.

- celltypes:

  Character vector of states/cell types; `NULL` for all.

- weight_cutoff:

  Numeric threshold for edge weight filtering.

- celltypes_order:

  The order of cell types.

- ntop:

  The number of top regulators to show (labels and edge filtering). When
  `filter_by_ntop=TRUE`, only edges from top `ntop` regulators are
  plotted to avoid huge networks. Default is \`10\`.

- filter_by_ntop:

  Logical. If `TRUE` (default), filter edges to top `ntop` regulators
  per celltype.

- title:

  The title of figure. Default is \`NULL\`.

- theme_type:

  The theme of figure. Could be \`"theme_void"\`, \`"theme_blank"\`, or
  \`"theme_facet"\`. Default is \`"theme_void"\`.

- plot_type:

  The type of figure. Could be \`"ggplot"\`, \`"animate"\`, or
  \`"ggplotly"\`. Default is \`"ggplot"\`.

- layout:

  The layout of figure. Could be \`"fruchtermanreingold"\` or
  \`"kamadakawai"\`. Default is \`"fruchtermanreingold"\`.

- nrow:

  The number of rows of figure (for `facet_wrap` when `combine=TRUE`).
  Default is \`2\`.

- ncol:

  Integer or `NULL`. Number of columns for `facet_wrap`.

- byrow:

  Logical. If `TRUE`, fill facets by row (`as.table=TRUE` in
  `facet_wrap`).

- combine:

  Logical. If `TRUE` (default), return a single faceted plot; if
  `FALSE`, return a named list of plots per celltype.

- figure_save:

  Whether to save the figure file. Default is \`FALSE\`.

- figure_name:

  The name of figure file. Default is \`NULL\`.

- figure_width:

  The width of figure. Default is \`6\`.

- figure_height:

  The height of figure. Default is \`6\`.

- seed:

  The seed random use to plot network. Default is \`1\`.

## Value

A dynamic figure object.

## Examples

``` r
data(example_matrix, package = "inferCSN")
network <- inferCSN::inferCSN(example_matrix, verbose = FALSE)
network$celltype <- rep(
  c("cluster5", "cluster3", "cluster2", "cluster1", "cluster6"),
  length.out = nrow(network)
)

celltypes_order <- c(
  "cluster5", "cluster3",
  "cluster2", "cluster1",
  "cluster6"
)

plot_dynamic_networks(
  network,
  celltypes_order = celltypes_order
)


plot_dynamic_networks(
  network,
  celltypes_order = celltypes_order[1:3]
)


plot_dynamic_networks(
  network,
  celltypes_order = celltypes_order,
  plot_type = "ggplotly"
)

{"x":{"data":[{"x":[0.33219424671294684,0.81772459646872142],"y":[0,0.12642705410946112],"text":"x: 0.3321942<br />y: 0.0000000<br />xend: 0.817724596<br />yend: 0.126427054<br />Interaction: Activation","type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.14769496527857964,0.64596452221429412,null,0.83224055579712186,0.34671020604134734],"y":[0.88251048555948375,0.99665065693185118,null,0.130206858969004,0.003779804859542879],"text":["x: 0.1476950<br />y: 0.8825105<br />xend: 0.645964522<br />yend: 0.996650657<br />Interaction: Activation","x: 0.1476950<br />y: 0.8825105<br />xend: 0.645964522<br />yend: 0.996650657<br />Interaction: Activation",null,"x: 0.8322406<br />y: 0.1302069<br />xend: 0.346710206<br />yend: 0.003779805<br />Interaction: Activation","x: 0.8322406<br />y: 0.1302069<br />xend: 0.346710206<br />yend: 0.003779805<br />Interaction: Activation"],"type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.66058580468238981,0.16231624774667536,null,0.14769496527857964,0.0043429937724203949],"y":[1,0.88585982862763257,null,0.88251048555948375,0.40860261732119535],"text":["x: 0.6605858<br />y: 1.0000000<br />xend: 0.162316248<br />yend: 0.885859829<br />Interaction: Activation","x: 0.6605858<br />y: 1.0000000<br />xend: 0.162316248<br />yend: 0.885859829<br />Interaction: Activation",null,"x: 0.1476950<br />y: 0.8825105<br />xend: 0.004342994<br />yend: 0.408602617<br />Interaction: Activation","x: 0.1476950<br />y: 0.8825105<br />xend: 0.004342994<br />yend: 0.408602617<br />Interaction: Activation"],"type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0,0.14335197150615925],"y":[0.39424509646651024,0.86815296470479864],"text":"x: 0.0000000<br />y: 0.3942451<br />xend: 0.143351972<br />yend: 0.868152965<br />Interaction: Activation","type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[1,0.67041124530335205,null,0.83224055579712186,0.99503507535131974],"y":[0.60847171988143267,0.98866594880001224,null,0.130206858969004,0.59431723213962706],"text":["x: 1.0000000<br />y: 0.6084717<br />xend: 0.670411245<br />yend: 0.988665949<br />Interaction: Repression","x: 1.0000000<br />y: 0.6084717<br />xend: 0.670411245<br />yend: 0.988665949<br />Interaction: Repression",null,"x: 0.8322406<br />y: 0.1302069<br />xend: 0.995035075<br />yend: 0.594317232<br />Interaction: Repression","x: 0.8322406<br />y: 0.1302069<br />xend: 0.995035075<br />yend: 0.594317232<br />Interaction: Repression"],"type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(255,0,102,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1,0.83720548044580212,null,0.66058580468238981,0.99017455937903776,null,0,0.32252883387926751],"y":[0.60847171988143267,0.14436134671080963,null,1,0.61980577108142043,null,0.39424509646651024,0.011470823630173488],"text":["x: 1.0000000<br />y: 0.6084717<br />xend: 0.837205480<br />yend: 0.144361347<br />Interaction: Repression","x: 1.0000000<br />y: 0.6084717<br />xend: 0.837205480<br />yend: 0.144361347<br />Interaction: Repression",null,"x: 0.6605858<br />y: 1.0000000<br />xend: 0.990174559<br />yend: 0.619805771<br />Interaction: Repression","x: 0.6605858<br />y: 1.0000000<br />xend: 0.990174559<br />yend: 0.619805771<br />Interaction: Repression",null,"x: 0.0000000<br />y: 0.3942451<br />xend: 0.322528834<br />yend: 0.011470824<br />Interaction: Repression","x: 0.0000000<br />y: 0.3942451<br />xend: 0.322528834<br />yend: 0.011470824<br />Interaction: Repression"],"type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(255,0,102,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,0.0096654128336793321],"y":[0,0.38277427283633675],"text":"x: 0.3321942<br />y: 0.0000000<br />xend: 0.009665413<br />yend: 0.382774273<br />Interaction: Repression","type":"scatter","mode":"lines","line":{"width":2.6456692913385824,"color":"rgba(255,0,102,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["x: 0.3321942<br />y: 0.0000000<br />targets_num: 1","x: 1.0000000<br />y: 0.6084717<br />targets_num: 1","x: 0.6605858<br />y: 1.0000000<br />targets_num: 0","x: 0.1476950<br />y: 0.8825105<br />targets_num: 0","x: 0.0000000<br />y: 0.3942451<br />targets_num: 0","x: 0.8322406<br />y: 0.1302069<br />targets_num: 1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,169,169,1)","opacity":0.90000000000000002,"size":[22.677165354330711,22.677165354330711,3.7795275590551185,3.7795275590551185,3.7795275590551185,22.677165354330711],"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(169,169,169,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["x: 0.3321942<br />y: 0.0000000<br />targets_num: 0","x: 1.0000000<br />y: 0.6084717<br />targets_num: 1","x: 0.6605858<br />y: 1.0000000<br />targets_num: 1","x: 0.1476950<br />y: 0.8825105<br />targets_num: 0","x: 0.0000000<br />y: 0.3942451<br />targets_num: 1","x: 0.8322406<br />y: 0.1302069<br />targets_num: 0"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,169,169,1)","opacity":0.90000000000000002,"size":[3.7795275590551185,22.677165354330711,22.677165354330711,3.7795275590551185,22.677165354330711,3.7795275590551185],"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(169,169,169,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["x: 0.3321942<br />y: 0.0000000<br />targets_num: 0","x: 1.0000000<br />y: 0.6084717<br />targets_num: 0","x: 0.6605858<br />y: 1.0000000<br />targets_num: 0","x: 0.1476950<br />y: 0.8825105<br />targets_num: 1","x: 0.0000000<br />y: 0.3942451<br />targets_num: 0","x: 0.8322406<br />y: 0.1302069<br />targets_num: 1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,169,169,1)","opacity":0.90000000000000002,"size":[3.7795275590551185,3.7795275590551185,3.7795275590551185,22.677165354330711,3.7795275590551185,22.677165354330711],"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(169,169,169,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["x: 0.3321942<br />y: 0.0000000<br />targets_num: 0","x: 1.0000000<br />y: 0.6084717<br />targets_num: 0","x: 0.6605858<br />y: 1.0000000<br />targets_num: 1","x: 0.1476950<br />y: 0.8825105<br />targets_num: 1","x: 0.0000000<br />y: 0.3942451<br />targets_num: 0","x: 0.8322406<br />y: 0.1302069<br />targets_num: 0"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,169,169,1)","opacity":0.90000000000000002,"size":[3.7795275590551185,3.7795275590551185,22.677165354330711,22.677165354330711,3.7795275590551185,3.7795275590551185],"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(169,169,169,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["x: 0.3321942<br />y: 0.0000000<br />targets_num: 1","x: 1.0000000<br />y: 0.6084717<br />targets_num: 0","x: 0.6605858<br />y: 1.0000000<br />targets_num: 0","x: 0.1476950<br />y: 0.8825105<br />targets_num: 0","x: 0.0000000<br />y: 0.3942451<br />targets_num: 1","x: 0.8322406<br />y: 0.1302069<br />targets_num: 0"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,169,169,1)","opacity":0.90000000000000002,"size":[22.677165354330711,3.7795275590551185,3.7795275590551185,3.7795275590551185,22.677165354330711,3.7795275590551185],"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(169,169,169,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["g1","g2","","","","g6"],"hovertext":["x: 0.3321942<br />y: 0.0000000<br />label_genes: g1","x: 1.0000000<br />y: 0.6084717<br />label_genes: g2","x: 0.6605858<br />y: 1.0000000<br />label_genes: ","x: 0.1476950<br />y: 0.8825105<br />label_genes: ","x: 0.0000000<br />y: 0.3942451<br />label_genes: ","x: 0.8322406<br />y: 0.1302069<br />label_genes: g6"],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["","g2","g3","","g5",""],"hovertext":["x: 0.3321942<br />y: 0.0000000<br />label_genes: ","x: 1.0000000<br />y: 0.6084717<br />label_genes: g2","x: 0.6605858<br />y: 1.0000000<br />label_genes: g3","x: 0.1476950<br />y: 0.8825105<br />label_genes: ","x: 0.0000000<br />y: 0.3942451<br />label_genes: g5","x: 0.8322406<br />y: 0.1302069<br />label_genes: "],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["","","","g4","","g6"],"hovertext":["x: 0.3321942<br />y: 0.0000000<br />label_genes: ","x: 1.0000000<br />y: 0.6084717<br />label_genes: ","x: 0.6605858<br />y: 1.0000000<br />label_genes: ","x: 0.1476950<br />y: 0.8825105<br />label_genes: g4","x: 0.0000000<br />y: 0.3942451<br />label_genes: ","x: 0.8322406<br />y: 0.1302069<br />label_genes: g6"],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["","","g3","g4","",""],"hovertext":["x: 0.3321942<br />y: 0.0000000<br />label_genes: ","x: 1.0000000<br />y: 0.6084717<br />label_genes: ","x: 0.6605858<br />y: 1.0000000<br />label_genes: g3","x: 0.1476950<br />y: 0.8825105<br />label_genes: g4","x: 0.0000000<br />y: 0.3942451<br />label_genes: ","x: 0.8322406<br />y: 0.1302069<br />label_genes: "],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.33219424671294684,1,0.66058580468238981,0.14769496527857964,0,0.83224055579712186],"y":[0,0.60847171988143267,1,0.88251048555948375,0.39424509646651024,0.130206858969004],"text":["g1","","","","g5",""],"hovertext":["x: 0.3321942<br />y: 0.0000000<br />label_genes: g1","x: 1.0000000<br />y: 0.6084717<br />label_genes: ","x: 0.6605858<br />y: 1.0000000<br />label_genes: ","x: 0.1476950<br />y: 0.8825105<br />label_genes: ","x: 0.0000000<br />y: 0.3942451<br />label_genes: g5","x: 0.8322406<br />y: 0.1302069<br />label_genes: "],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":27.68949771689498,"r":0,"b":0,"l":0},"paper_bgcolor":"rgba(0,0,0,0)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,0.32191780821917804],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":0,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y2","title":"","hoverformat":".2f"},"yaxis":{"domain":[0.52968036529680362,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":0,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.32191780821917804,"y0":0,"y1":11.68949771689498,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.34474885844748859,"x1":0.65525114155251141,"y0":0,"y1":11.68949771689498,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.67808219178082185,"x1":1,"y0":0,"y1":11.68949771689498,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.32191780821917804,"y0":0,"y1":11.68949771689498,"yanchor":0.47031963470319632,"ysizemode":"pixel"},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.34474885844748859,"x1":0.65525114155251141,"y0":0,"y1":11.68949771689498,"yanchor":0.47031963470319632,"ysizemode":"pixel"}],"annotations":[{"text":"cluster5","x":0.16095890410958902,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster3","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster2","x":0.83904109589041087,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster1","x":0.16095890410958902,"y":0.47031963470319632,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster6","x":0.5,"y":0.47031963470319632,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"xaxis2":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":0,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.34474885844748859,0.65525114155251141],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y2","title":"","hoverformat":".2f"},"xaxis3":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":0,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.67808219178082185,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis2":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":0,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0,0.47031963470319632],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"title":{"text":"Interaction<br />targets_num","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1eca6c51c979":{"x":{},"y":{},"xend":{},"yend":{},"colour":{},"type":"scatter"},"1eca7ccddc18":{"x":{},"y":{},"xend":{},"yend":{},"size":{}},"1eca682adcc":{"x":{},"y":{},"xend":{},"yend":{},"label":{}}},"cur_data":"1eca6c51c979","visdat":{"1eca6c51c979":["function (y) ","x"],"1eca7ccddc18":["function (y) ","x"],"1eca682adcc":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
