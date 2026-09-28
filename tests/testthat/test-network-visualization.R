test_that("static network labels contain one row per node", {
  edges <- data.frame(
    regulator = c("A", "A", "A"), target = c("B", "C", "D"),
    weight = c(1, -2, 3)
  )
  plot <- plot_static_networks(edges)
  labels <- plot$layers[[length(plot$layers)]]$data
  expect_equal(nrow(labels), length(unique(c(edges$regulator, edges$target))))
  expect_equal(length(unique(labels$name)), nrow(labels))
  expect_s3_class(plot$theme$panel.border, "element_blank")
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

test_that("contrast network preserves directed degree and magnitude filtering", {
  edges <- data.frame(
    regulator = c("A", "A", "B", "C"), target = c("B", "C", "C", "A"),
    weight = c(1, -2, 3, 1)
  )
  plot <- plot_contrast_networks(edges)
  expect_setequal(plot$data$name, c("A", "B", "C"))
  built <- ggplot2::ggplot_build(plot)
  expect_setequal(unique(built$data[[1]]$edge_colour), c("#1F78B4", "#E31A1C"))
  thresholded <- ggplot2::ggplot_build(plot_contrast_networks(edges, weight_value = 2.5))
  expect_true(all(thresholded$data[[1]]$edge_colour == "#1F78B4"))
})

test_that("dynamic network honors palette and legacy theme choices", {
  edges <- data.frame(
    regulator = c("A", "A", "B"), target = c("B", "C", "C"),
    weight = c(1, -2, 3), celltype = c("s1", "s1", "s2")
  )
  plot <- plot_dynamic_networks(
    edges, palcolor = c(Regulator = "red", Target = "blue")
  )
  fills <- vapply(plot$layers[2:3], function(layer) layer$aes_params$fill, character(1))
  expect_identical(unname(fills), c("blue", "red"))
  expect_s3_class(plot$theme$panel.border, "element_blank")
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  legacy <- plot_dynamic_networks(edges, theme_type = "theme_facet")
  expect_s3_class(legacy$theme$panel.border, "element_blank")
  expect_s3_class(ggplot2::ggplot_build(legacy), "ggplot_built")
})

test_that("static network supports state fills and focal TF diamonds", {
  edges <- data.frame(
    regulator = c("A", "A", "B"),
    target = c("B", "C", "C"),
    weight = c(1, -2, 3)
  )
  plot <- plot_static_networks(
    edges,
    tf_nodes = "A",
    node_state = c(A = "HSC", B = "ProE", C = "ProE"),
    state_colors = c(HSC = "#3569A8", ProE = "#E4A33C"),
    seed = 42
  )
  built <- ggplot2::ggplot_build(plot)
  shapes <- unique(unlist(lapply(built$data, function(x) {
    if ("shape" %in% names(x)) x$shape else NULL
  })))
  fills <- unique(unlist(lapply(built$data, function(x) {
    if ("fill" %in% names(x)) x$fill else NULL
  })))
  expect_true(all(c(21, 23) %in% shapes))
  expect_true(all(c("#3569A8", "#E4A33C") %in% fills))
})

test_that("static network can label only selected nodes", {
  edges <- data.frame(
    regulator = c("A", "A", "B"),
    target = c("B", "C", "C"),
    weight = c(1, -2, 3)
  )
  plot <- plot_static_networks(
    edges,
    tf_nodes = "A",
    node_state = c(A = "HSC", B = "ProE", C = "ProE"),
    state_colors = c(HSC = "#3569A8", ProE = "#E4A33C"),
    label_nodes = "A",
    seed = 42
  )
  built <- ggplot2::ggplot_build(plot)
  label_values <- unique(unlist(lapply(built$data, function(layer) {
    if ("label" %in% names(layer)) {
      as.character(layer$label)
    } else {
      NULL
    }
  })))
  expect_true("A" %in% label_values)
  expect_false(any(c("B", "C") %in% label_values))
})

test_that("static network still draws unlabelled nodes", {
  edges <- data.frame(
    regulator = c("A", "A", "B"),
    target = c("B", "C", "C"),
    weight = c(1, -2, 3)
  )
  plot <- plot_static_networks(
    edges,
    tf_nodes = "A",
    node_state = c(A = "HSC", B = "ProE", C = "ProE"),
    state_colors = c(HSC = "#3569A8", ProE = "#E4A33C"),
    label_nodes = "A",
    seed = 42
  )
  built <- ggplot2::ggplot_build(plot)
  node_layer <- which(vapply(
    plot$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)
  ))[[1L]]
  expect_equal(nrow(built$data[[node_layer]]), length(unique(c(edges$regulator, edges$target))))
})
