state_partition_fixture <- function() {
  cells <- paste0("c", 1:6)
  list(
    matrix = rbind(g1 = c(5, 4, 3, 1, 1, 1), g2 = c(1, 1, 1, 3, 4, 5), g3 = rep(1, 6)),
    object = list(
      genes = c(g1 = 0.001, g2 = 0.001, g3 = 0.001),
      cells = data.frame(cell_name = cells, cells = cells,
        pseudotime = seq(0, 1, length.out = 6), group = rep(c("A", "B", "C"), each = 2),
        row.names = cells)
    )
  )
}

test_that("state partitions preserve cut boundaries and custom group labels", {
  f <- state_partition_fixture()
  colnames(f$matrix) <- rownames(f$object$cells)
  observed <- define_states(f$object, f$matrix, pseudotime_cuts = 0.5)
  expect_identical(observed$cells$state, rep(c("state_1", "state_2"), each = 3))
  expect_false("epoch" %in% names(observed$cells))
  expect_identical(observed$genes, f$object$genes)
  expect_identical(define_states_new(f$object, f$matrix, pseudotime_cuts = 0.5), observed)
  ordered <- define_states(f$object, f$matrix, method = "cell_order", num_states = 3)
  expect_identical(ordered$cells$state, rep(paste0("state_", 1:3), each = 2))
  grouped <- define_states(f$object, f$matrix, method = "group",
    group_assignments = list(Early = "A", Late = c("B", "C")))
  expect_identical(grouped$cells$state, c(rep("Early", 2), rep("Late", 4)))
  boundary <- split_states_by_pseudotime(f$object, 0.4, state_names = c("Early", "Late"))
  expect_identical(boundary$cells$state, c(rep("Early", 3), rep("Late", 3)))
  expect_error(split_states_by_pseudotime(f$object, 0.5, "one"), "state_names")
})

test_that("gene assignments and transition networks use the same state identifiers", {
  f <- state_partition_fixture()
  colnames(f$matrix) <- rownames(f$object$cells)
  object <- define_states(f$object, f$matrix, pseudotime_cuts = 0.5)
  assigned <- assign_genes_to_states(f$matrix, object, forceGenes = FALSE)
  expect_identical(assigned$state_1, c("g1", "g3"))
  expect_identical(assigned$state_2, c("g2", "g3"))
  expect_equal(assigned$mean_expression$mean_expression, c(4, 1, 1, 1, 4, 1))
  expect_identical(assigned$mean_expression$state, rep(c("state_1", "state_2"), each = 3))
  edges <- data.frame(regulator = c("g1", "g1", "g3", "g2"),
    TG = c("g2", "g3", "g2", "g1"), weight = c(0.2, 0.4, 0.6, 0.8))
  invisible(capture.output(networks <- split_network_by_states(edges, assigned)))
  expect_setequal(names(networks), c("state_1..state_2", "state_1..state_1", "state_2..state_2"))
  expect_identical(networks[["state_1..state_2"]], edges[c(1, 3), ])
  object$cells$epoch <- sub("state_", "epoch", object$cells$state)
  object$cells$state <- NULL
  expect_identical(assign_genes_to_states(f$matrix, object, forceGenes = FALSE), assigned)
  simple <- suppressMessages(assign_genes_to_states_simple(f$matrix, object, num_states = 2))
  expect_identical(simple$genes$state, c("state_1", "state_2", "state_1"))
  expect_identical(simple$cells$state, rep(c("state_1", "state_2"), each = 3))
  expect_false(anyDuplicated(names(simple$cells)) > 0L)
})

test_that("saved assignments migrate without changing custom labels or input data", {
  legacy <- list(
    genes = data.frame(gene = c("g1", "g2"), epoch = c("1", "2"), mean_expression = c(4, 5)),
    cells = data.frame(epoch = factor(c("epoch1", "epoch2"))),
    networks = list(epoch1..epoch2 = data.frame(weight = c(0.2, 0.5)))
  )
  migrated <- multiCSN:::normalize_state_data(legacy)
  expect_identical(as.character(migrated$cells$state), c("state_1", "state_2"))
  expect_identical(migrated$genes$state, c("state_1", "state_2"))
  expect_identical(names(migrated$networks), "state_1..state_2")
  expect_identical(migrated$networks[[1]]$weight, legacy$networks[[1]]$weight)
  expect_identical(multiCSN:::normalize_state_data(migrated), migrated)
  expect_true("epoch" %in% names(legacy$cells))
  custom <- data.frame(state = c("HSC", "MK/E", "Ery", NA))
  expect_identical(multiCSN:::normalize_state_data(custom), custom)
  expect_error(multiCSN:::normalize_state_data(data.frame(epoch = "epoch1", state = "state_2")),
    "assignments disagree")
  expect_error(multiCSN:::normalize_state_data(list(epoch1 = "g1", state_1 = "g2")),
    "identifiers collide")
})

test_that("heatmap annotations consume saved state assignments", {
  skip_if_not_installed("viridis")
  f <- state_partition_fixture()
  colnames(f$matrix) <- rownames(f$object$cells)
  meta <- f$object$cells
  meta$epoch <- rep(c("epoch1", "epoch2"), each = 3)
  plot <- hm_dyn(f$matrix, list(meta_data = meta, genes = rownames(f$matrix)),
    cluster = "group", topX = 3)
  expect_s4_class(plot, "Heatmap")
  expect_true("state" %in% names(plot@top_annotation@anno_list))
  expect_false("epoch" %in% names(plot@top_annotation@anno_list))
  expect_identical(plot@matrix, f$matrix[c("g1", "g3", "g2"), ])
})
