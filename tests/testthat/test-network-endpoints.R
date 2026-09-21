empty <- function(x) suppressWarnings(unlink(x, recursive = TRUE))

write_layer <- function(root, file, table) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(table, file.path(root, file), sep = "\t", quote = FALSE, na = "NA")
}

test_that("endpoint spec matches the layered fit output contract", {
  spec <- network_endpoints()
  expect_identical(nrow(spec), 5L)
  expect_identical(
    spec$file,
    c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv", "mediated_tf_gene.tsv")
  )
  expect_identical(spec$keys[[1]], c("regulator", "target"))
  expect_identical(spec$keys[[3]], c("region", "target"))
  expect_identical(spec$keys[[4]], c("regulator", "region", "target"))
})

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
    read_network_endpoints(root, endpoints = "TF-gene", direction = TRUE,
                           weight_column = "missing_column"),
    "lacks column"
  )
  expect_error(read_network_endpoints(root, endpoints = "nope"), "Unknown endpoint")
  expect_error(read_network_endpoints(file.path(root, "absent")), "does not exist")

  write_layer(root, "region_gene.tsv", data.frame(
    region = c("chr1:1-2", "chr1:1-2"), target = c("g1", "g1"),
    standardized_beta = c(0.5, 0.5), deletion_delta_bic = c(3, 3), weight = c(0.5, 0.5)
  ))
  expect_error(read_network_endpoints(root, endpoints = "region-gene"),
               "missing or duplicated keys")
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
    "lacks column"
  )
  expect_error(
    classify_edge_transitions(previous, next_state, c("regulator", "target"),
                              direction_column = "absent"),
    "lacks column"
  )
  expect_error(
    classify_edge_transitions(rbind(previous, previous), next_state,
                              c("regulator", "target")),
    "missing or duplicated edge keys"
  )
  expect_error(
    classify_edge_transitions(previous, next_state, c("regulator", "target"),
                              labels = c(a = "a", b = "b", c = "c", d = "d")),
    "labels must name"
  )
})
