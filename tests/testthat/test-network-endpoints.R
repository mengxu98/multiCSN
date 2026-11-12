empty <- function(x) suppressWarnings(unlink(x, recursive = TRUE))

write_layer <- function(root, file, table) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(table, file.path(root, file), sep = "\t", quote = FALSE, na = "NA")
}

test_that("read_network_endpoints reads and validates a result directory", {
  root <- file.path(tempdir(), "endpoint-reader")
  empty(root)
  write_layer(root, "tf_gene.tsv", data.frame(
    regulator = c("tf1", "tf2"), target = c("g1", "g2"),
    standardized_beta = c(0.5, -0.25), deletion_delta_bic = c(3, 4), weight = c(0.75, -0.5)
  ))
  write_layer(root, "tf_region.tsv", data.frame(
    regulator = "tf1", region = "chr1:1-2", standardized_beta = 0.5,
    deletion_delta_bic = 3, weight = 0.75
  ))
  write_layer(root, "region_gene.tsv", data.frame(
    region = "chr1:1-2", target = "g1", standardized_beta = 0.5,
    deletion_delta_bic = 3, weight = 0.75
  ))
  write_layer(root, "triplets.tsv", data.frame(
    regulator = "tf1", region = "chr1:1-2", target = "g1",
    chain_delta_bic = 3, weight = 0.75
  ))
  write_layer(root, "mediated_tf_gene.tsv", data.frame(
    regulator = "tf1", target = "g1", region = "chr1:1-2",
    chain_delta_bic = 3, weight = 0.75
  ))
  tables <- read_network_endpoints(root)
  expect_identical(names(tables), network_endpoints()$endpoint)
  expect_identical(sort(names(tables$`TF-gene`)), sort(c(
    "regulator", "target", "standardized_beta", "deletion_delta_bic", "weight"
  )))
  directed <- read_network_endpoints(root, endpoints = "TF-gene", direction = TRUE)
  expect_identical(names(directed$`TF-gene`), c("regulator", "target", "direction"))
  expect_identical(directed$`TF-gene`$direction, c(1, -1))
  expect_error(
    read_network_endpoints(root,
      endpoints = "TF-gene", direction = TRUE,
      weight_column = "missing_column"
    ),
    "lacks column"
  )
  expect_error(read_network_endpoints(root, endpoints = "nope"), "Unknown endpoint")
  expect_error(read_network_endpoints(file.path(root, "absent")), "does not exist")

  write_layer(root, "region_gene.tsv", data.frame(
    region = c("chr1:1-2", "chr1:1-2"), target = c("g1", "g1"),
    standardized_beta = c(0.5, 0.5), deletion_delta_bic = c(3, 3), weight = c(0.5, 0.5)
  ))
  expect_error(
    read_network_endpoints(root, endpoints = "region-gene"),
    "missing or duplicated keys"
  )
  write_layer(root, "region_gene.tsv", data.frame(
    region = "chr1:1-2", target = "g1", standardized_beta = 0.5,
    deletion_delta_bic = 3, weight = 0
  ))
  expect_error(
    read_network_endpoints(root, endpoints = "region-gene", direction = TRUE),
    "finite and non-zero"
  )
  empty(root)
})

test_that("classify_edge_transitions labels shared, gained, lost and flipped edges", {
  previous <- data.frame(
    regulator = c("tf1", "tf2", "tf3"), target = c("g1", "g1", "g1"),
    direction = c(1, -1, 1)
  )
  next_state <- data.frame(
    regulator = c("tf1", "tf2", "tf4"), target = c("g1", "g1", "g1"),
    direction = c(1, 1, -1)
  )
  transitions <- classify_edge_transitions(previous, next_state, c("regulator", "target"))
  expect_identical(
    names(transitions),
    c("regulator", "target", "previous_direction", "next_direction", "topology_class")
  )
  expect_identical(nrow(transitions), 4L)
  class_of <- function(regulator) {
    transitions$topology_class[transitions$regulator == regulator]
  }
  expect_identical(class_of("tf1"), "shared_sign_stable")
  expect_identical(class_of("tf2"), "shared_sign_flip")
  expect_identical(class_of("tf3"), "lost")
  expect_identical(class_of("tf4"), "gained")
  expect_true(is.na(transitions$next_direction[transitions$regulator == "tf3"]))
  expect_true(is.na(transitions$previous_direction[transitions$regulator == "tf4"]))
  relabelled <- classify_edge_transitions(
    previous, next_state, c("regulator", "target"),
    labels = c(gained = "new", lost = "gone", stable = "kept", flip = "reversed")
  )
  expect_identical(sort(unique(relabelled$topology_class)), c("gone", "kept", "new", "reversed"))
})

test_that("classify_edge_transitions validates its inputs", {
  previous <- data.frame(regulator = "tf1", target = "g1", direction = 1)
  next_state <- data.frame(regulator = "tf1", target = "g1", direction = 1)
  expect_error(
    classify_edge_transitions(previous, next_state, c("regulator", "missing")),
    "column"
  )
  expect_error(
    classify_edge_transitions(previous, next_state, c("regulator", "target"),
      direction_column = "absent"
    ),
    "column"
  )
  expect_error(
    classify_edge_transitions(
      rbind(previous, previous), next_state,
      c("regulator", "target")
    ),
    "missing or duplicated edge keys"
  )
  expect_error(
    classify_edge_transitions(previous, next_state, c("regulator", "target"),
      labels = c(a = "a", b = "b", c = "c", d = "d")
    ),
    "labels must name"
  )
})

fake_layered_fit <- function() {
  edge <- function(...) data.frame(..., stringsAsFactors = FALSE)
  fit <- list(
    tf_gene = edge(regulator = c("tf1", "tf2"), target = c("g1", "g2"), weight = c(0.75, -0.25)),
    tf_region = edge(regulator = "tf1", region = "chr1-1-2", weight = 0.5),
    region_gene = edge(region = "chr1-1-2", target = "g1", weight = 0.5),
    triplets = edge(regulator = "tf1", region = "chr1-1-2", target = "g1", weight = 0.5),
    mediated_tf_gene = edge(regulator = "tf1", target = "g1", region = "chr1-1-2", weight = 0.5),
    accounting = list(
      tf_gene = edge(endpoint = "TF-gene", target = c("g1", "g2"), n_selected = 1L),
      tf_region = edge(endpoint = "TF-region", region = "chr1-1-2", n_selected = 1L),
      region_gene = edge(endpoint = "region-gene", target = "g1", n_selected = 1L)
    )
  )
  fit$summary <- edge(
    endpoint = network_endpoints()$endpoint,
    selected = c(2L, 1L, 1L, 1L, 1L)
  )
  fit
}

test_that("written endpoints round-trip in every format and validate", {
  fit <- fake_layered_fit()
  for (format in c("tsv", "csv", "tsv.gz", "csv.gz")) {
    root <- tempfile(paste0("endpoints-", format, "-"))
    write_network_endpoints(fit, root, format = format)
    expect_true(validate_network_endpoints(root))
    tables <- read_network_endpoints(root)
    expect_equal(as.data.frame(tables$`TF-gene`), fit$tf_gene)
    single <- read_network_endpoints(root, endpoints = "TF-gene", direction = TRUE, simplify = TRUE)
    expect_identical(single$direction, c(1, -1))
    unlink(root, recursive = TRUE)
  }
})

test_that("validate_network_endpoints reports incomplete or inconsistent output", {
  root <- tempfile("endpoints-invalid-")
  write_network_endpoints(fake_layered_fit(), root)
  expect_false(validate_network_endpoints(file.path(root, "absent"), error = FALSE))
  summary <- data.table::fread(file.path(root, "endpoint_summary.tsv"))
  summary$selected[[1]] <- 5L
  data.table::fwrite(summary, file.path(root, "endpoint_summary.tsv"), sep = "\t")
  expect_error(validate_network_endpoints(root), "count mismatch")
  file.create(file.path(root, "tf_region_accounting.tsv"))
  expect_error(validate_network_endpoints(root), "Incomplete network output")
  unlink(root, recursive = TRUE)
})

test_that("classify_edge_transitions labels the transition when states are named", {
  previous <- data.frame(regulator = c("a", "b"), target = "x", direction = c(1, -1))
  next_state <- data.frame(regulator = c("b", "c"), target = "x", direction = c(1, 1))
  merged <- classify_edge_transitions(previous, next_state, from = "early", to = "late")
  expect_identical(unique(merged$transition), "early -> late")
  expect_identical(unique(merged$previous_state), "early")
  expect_identical(sort(merged$topology_class), sort(c("lost", "shared_sign_flip", "gained")))
  expect_error(classify_edge_transitions(previous, next_state, from = "early"), "single state names")
})
