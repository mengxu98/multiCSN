test_that("candidate preparation rejects split RNA data before motif scanning", {
  object <- multiome_fixture()
  object[["RNA"]] <- split(object[["RNA"]], f = object$celltype)
  expect_error(
    prepare_candidates(object, group.by = "celltype", genome = NULL, verbose = FALSE),
    "Assay 'RNA' requires one joined 'data' layer.*JoinLayers"
  )
})

test_that("layered fitting rejects split RNA data before using candidates", {
  object <- initiate_object(
    multiome_fixture(), group.by = "celltype", filter_mode = "unfiltered",
    filter_by = "aggregate", rna_assay = "RNA", peak_assay = "ATAC",
    verbose = FALSE
  )
  object[["RNA"]] <- split(object[["RNA"]], f = object$celltype)
  expect_error(
    fit_layered_network(object, cell_group = "A", preset = "state", verbose = FALSE),
    "Assay 'RNA' requires one joined 'data' layer.*JoinLayers"
  )
})

test_that("per-group normalization requires joined counts", {
  object <- initiate_object(
    multiome_fixture(), group.by = "celltype", filter_mode = "unfiltered",
    filter_by = "aggregate", rna_assay = "RNA", peak_assay = "ATAC",
    verbose = FALSE
  )
  object[["RNA"]] <- split(object[["RNA"]], f = object$celltype, layers = "counts")
  expect_error(
    fit_layered_network(object, cell_group = "A", preset = "timepoint", verbose = FALSE),
    "Assay 'RNA' requires one joined 'counts' layer.*JoinLayers"
  )
})
