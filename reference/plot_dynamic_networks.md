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
  theme_use = "theme_this",
  theme_type = NULL,
  plot_type = "ggplot",
  layout = "anchor",
  align_nodes = TRUE,
  cols = NULL,
  palette = "Chinese",
  palcolor = NULL,
  regulator_color = NULL,
  target_color = NULL,
  aspect.ratio = 1,
  aspect_ratio = NULL,
  label.size = 3.5,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.12,
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

- theme_use:

  Theme function name or choice: `"theme_this"`, `"theme_blank"`, or
  `"theme_void"`. Default is `"theme_this"`.

- theme_type:

  Legacy theme choice: `"theme_void"`, `"theme_blank"`, or
  `"theme_facet"`.

- plot_type:

  The type of figure. Could be \`"ggplot"\`, \`"animate"\`, or
  \`"ggplotly"\`. Default is \`"ggplot"\`.

- layout:

  The layout of figure. Could be \`"anchor"\` (temporally aligned
  coordinates across states), \`"fruchtermanreingold"\` or
  \`"kamadakawai"\`. Default is \`"anchor"\`.

- align_nodes:

  Logical. If `TRUE` (default), anchor node positions across all states
  based on the union graph to ensure temporal coherence and avoid
  jumping nodes.

- cols:

  Optional named color vector for Interaction (`Activation` and
  `Repression`). Default is
  `c("Activation" = "#1F78B4", "Repression" = "#E31A1C")`.

- palette:

  Palette name used by
  [`thisplot::palette_colors()`](https://mengxu98.github.io/thisplot/reference/palette_colors.html).
  Default is `"Chinese"`.

- palcolor:

  Optional named color vector overriding default node colors.

- regulator_color:

  Fill color for regulator nodes. Default is `NULL` (resolved via
  `palette`).

- target_color:

  Fill color for target nodes. Default is `NULL` (resolved via
  `palette`).

- aspect.ratio:

  Numeric aspect ratio for panels. Default is `1`.

- aspect_ratio:

  Alias for `aspect.ratio`.

- label.size:

  Label text size. Default is `3.5`.

- label.fg:

  Label foreground text color. Default is `"black"`.

- label.bg:

  Label outline/background color. Default is `"white"`.

- label.bg.r:

  Label outline radius. Default is `0.12`.

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
#> Warning: geom_GeomTextRepel() has yet to be implemented in plotly.
#>   If you'd like to see this geom implemented,
#>   Please open an issue with your example code at
#>   https://github.com/ropensci/plotly/issues
#> Warning: geom_GeomTextRepel() has yet to be implemented in plotly.
#>   If you'd like to see this geom implemented,
#>   Please open an issue with your example code at
#>   https://github.com/ropensci/plotly/issues
#> Warning: geom_GeomTextRepel() has yet to be implemented in plotly.
#>   If you'd like to see this geom implemented,
#>   Please open an issue with your example code at
#>   https://github.com/ropensci/plotly/issues
#> Warning: geom_GeomTextRepel() has yet to be implemented in plotly.
#>   If you'd like to see this geom implemented,
#>   Please open an issue with your example code at
#>   https://github.com/ropensci/plotly/issues
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width
#> Warning: Aspect ratios aren't yet implemented, but you can manually set a suitable height/width

{"x":{"data":[{"x":[0,0.26468930604579755],"y":[0.44947179223390687,0.02693097255202298],"text":"x: 0.0000000<br />y: 0.44947179<br />xend: 0.264689306<br />yend: 0.0269309726<br />Interaction: Activation","type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(31,120,180,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.70341031487267713,0.99218658158723638,null,0.27265229416232589,0.0079629881165283467],"y":[1,0.52676508514466902,null,0.014219134751923799,0.43675995443380772],"text":["x: 0.7034103<br />y: 1.00000000<br />xend: 0.992186582<br />yend: 0.5267650851<br />Interaction: Activation","x: 0.7034103<br />y: 1.00000000<br />xend: 0.992186582<br />yend: 0.5267650851<br />Interaction: Activation",null,"x: 0.2726523<br />y: 0.01421913<br />xend: 0.007962988<br />yend: 0.4367599544<br />Interaction: Activation","x: 0.2726523<br />y: 0.01421913<br />xend: 0.007962988<br />yend: 0.4367599544<br />Interaction: Activation"],"type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(31,120,180,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1,0.71122373328544075,null,0.70341031487267713,0.27665230405322505],"y":[0.51396076975074467,0.98719568460607565,null,1,0.97195146153806722],"text":["x: 1.0000000<br />y: 0.51396077<br />xend: 0.711223733<br />yend: 0.9871956846<br />Interaction: Activation","x: 1.0000000<br />y: 0.51396077<br />xend: 0.711223733<br />yend: 0.9871956846<br />Interaction: Activation",null,"x: 0.7034103<br />y: 1.00000000<br />xend: 0.276652304<br />yend: 0.9719514615<br />Interaction: Activation","x: 0.7034103<br />y: 1.00000000<br />xend: 0.276652304<br />yend: 0.9719514615<br />Interaction: Activation"],"type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(31,120,180,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.26168459747923606,0.68844260829868809],"y":[0.97096771364609136,0.99901625210802414],"text":"x: 0.2616846<br />y: 0.97096771<br />xend: 0.688442608<br />yend: 0.9990162521<br />Interaction: Activation","type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(31,120,180,1)","dash":"solid"},"hoveron":"points","name":"Activation","legendgroup":"Activation","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"x":[0.66721983195133527,0.99184747850305466,null,0.27265229416232589,0.65222956260785647],"y":[0,0.50136964616141233,null,0.014219134751923799,0.00054020830237822698],"text":["x: 0.6672198<br />y: 0.00000000<br />xend: 0.991847479<br />yend: 0.5013696462<br />Interaction: Repression","x: 0.6672198<br />y: 0.00000000<br />xend: 0.991847479<br />yend: 0.5013696462<br />Interaction: Repression",null,"x: 0.2726523<br />y: 0.01421913<br />xend: 0.652229563<br />yend: 0.0005402083<br />Interaction: Repression","x: 0.2726523<br />y: 0.01421913<br />xend: 0.652229563<br />yend: 0.0005402083<br />Interaction: Repression"],"type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(227,26,28,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.66721983195133527,0.28764256350580469,null,1,0.67537235344828073,null,0.26168459747923606,0.0067274604729195553],"y":[0,0.013678926449545572,null,0.51396076975074467,0.012591123589332343,null,0.97096771364609136,0.46287855454198501],"text":["x: 0.6672198<br />y: 0.00000000<br />xend: 0.287642564<br />yend: 0.0136789264<br />Interaction: Repression","x: 0.6672198<br />y: 0.00000000<br />xend: 0.287642564<br />yend: 0.0136789264<br />Interaction: Repression",null,"x: 1.0000000<br />y: 0.51396077<br />xend: 0.675372353<br />yend: 0.0125911236<br />Interaction: Repression","x: 1.0000000<br />y: 0.51396077<br />xend: 0.675372353<br />yend: 0.0125911236<br />Interaction: Repression",null,"x: 0.2616846<br />y: 0.97096771<br />xend: 0.006727460<br />yend: 0.4628785545<br />Interaction: Repression","x: 0.2616846<br />y: 0.97096771<br />xend: 0.006727460<br />yend: 0.4628785545<br />Interaction: Repression"],"type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(227,26,28,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,0.2549571370063165],"y":[0.44947179223390687,0.95756095133801322],"text":"x: 0.0000000<br />y: 0.44947179<br />xend: 0.254957137<br />yend: 0.9575609513<br />Interaction: Repression","type":"scatter","mode":"lines","line":{"width":2.4566929133858273,"color":"rgba(227,26,28,1)","dash":"solid"},"hoveron":"points","name":"Repression","legendgroup":"Repression","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"visible":false,"showlegend":false,"xaxis":null,"yaxis":null,"hoverinfo":"text","frame":null},{"x":[0,0.66721983195133527,1,0.70341031487267713,0.26168459747923606,0.27265229416232589],"y":[0.44947179223390687,0,0.51396076975074467,1,0.97096771364609136,0.014219134751923799],"text":["x: 0.0000000<br />y: 0.44947179<br />pmax(targets_num, 3): 3","x: 0.6672198<br />y: 0.00000000<br />pmax(targets_num, 3): 3","x: 1.0000000<br />y: 0.51396077<br />pmax(targets_num, 3): 3","x: 0.7034103<br />y: 1.00000000<br />pmax(targets_num, 3): 3","x: 0.2616846<br />y: 0.97096771<br />pmax(targets_num, 3): 3","x: 0.2726523<br />y: 0.01421913<br />pmax(targets_num, 3): 3"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(176,213,223,1)","opacity":0.94999999999999996,"size":20.138937164395049,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(38,38,38,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null},{"showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"showlegend":false,"xaxis":"x","yaxis":"y2","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":39.910336239103373,"r":7.3059360730593621,"b":10.958904109589042,"l":10.958904109589042},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243},"xaxis":{"domain":[0,0.32191780821917804],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y2","title":"","hoverformat":".2f"},"yaxis":{"domain":[0.53735990037359904,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"transparent","width":0.66417600664176002,"linetype":"none"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.32191780821917804,"y0":0,"y1":24.574512245745126,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"transparent","width":0.66417600664176002,"linetype":"none"},"yref":"paper","xref":"paper","layer":"below","x0":0.34474885844748859,"x1":0.65525114155251141,"y0":0,"y1":24.574512245745126,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"transparent","width":0.66417600664176002,"linetype":"none"},"yref":"paper","xref":"paper","layer":"below","x0":0.67808219178082185,"x1":1,"y0":0,"y1":24.574512245745126,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"transparent","width":0.66417600664176002,"linetype":"none"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.32191780821917804,"y0":0,"y1":24.574512245745126,"yanchor":0.46264009962640101,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"transparent","width":0.66417600664176002,"linetype":"none"},"yref":"paper","xref":"paper","layer":"below","x0":0.34474885844748859,"x1":0.65525114155251141,"y0":0,"y1":24.574512245745126,"yanchor":0.46264009962640101,"ysizemode":"pixel"}],"annotations":[{"text":"cluster5","x":0.16095890410958902,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":16.604400166044005},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster3","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":16.604400166044005},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster2","x":0.83904109589041087,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":16.604400166044005},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster1","x":0.16095890410958902,"y":0.46264009962640101,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":16.604400166044005},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"cluster6","x":0.5,"y":0.46264009962640101,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":16.604400166044005},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"xaxis2":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"domain":[0.34474885844748859,0.65525114155251141],"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y2","title":"","hoverformat":".2f"},"xaxis3":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"domain":[0.67808219178082185,1],"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis2":{"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"domain":[0,0.46264009962640101],"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"orientation":"h","x":0.5,"y":-0.14999999999999999,"xanchor":"center","yanchor":"top","title":{"text":"Interaction","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1ff323fc32ea":{"x":{},"y":{},"xend":{},"yend":{},"colour":{},"type":"scatter"},"1ff31ef1c9f4":{"x":{},"y":{},"xend":{},"yend":{},"size":{}},"1ff3355f4f1c":{"x":{},"y":{},"xend":{},"yend":{},"size":{}},"1ff369a631f":{"x":{},"y":{},"xend":{},"yend":{},"label":{}}},"cur_data":"1ff323fc32ea","visdat":{"1ff323fc32ea":["function (y) ","x"],"1ff31ef1c9f4":["function (y) ","x"],"1ff3355f4f1c":["function (y) ","x"],"1ff369a631f":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
