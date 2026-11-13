test_that("prepare_candidates and fit_layered_network run end to end", {
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
  object <- prepare_candidates(
    multiome_fixture(), group.by = "celltype", genome = genome, verbose = FALSE
  )
  expect_setequal(get_attribute(object, attribute = "celltypes"), c("A", "B"))
  expect_true(length(get_attribute(object, attribute = "tfs")) > 0L)

  root <- tempfile("layered-")
  output <- file.path(root, "A")
  fit <- fit_layered_network(
    object, cell_group = "A", preset = "state", output_dir = output, verbose = FALSE
  )
  expect_true(validate_network_endpoints(output))
  expect_identical(list.files(root, all.files = TRUE, no.. = TRUE), "A")
  expect_true(all(file.exists(fit$paths)))
  expect_identical(fit$settings$tf_region_scope, "domain_peaks")
  expect_true(fit$settings$sort_regulators)
  expect_error(
    fit_layered_network(object, cell_group = "A", output_dir = output, verbose = FALSE),
    "will not be overwritten"
  )
  in_memory <- fit_layered_network(
    object, cell_group = "A",
    min_detected = 0L, tf_region_scope = "domain_peaks", sort_regulators = TRUE,
    verbose = FALSE
  )
  expect_identical(in_memory$tf_gene, fit$tf_gene)
  expect_identical(in_memory$triplets, fit$triplets)
  positional <- fit_layered_network(
    object, "A", in_memory$cells, NULL, NULL, 0L, 3L, 1e5, 0L, "domain_peaks",
    FALSE, FALSE, 1024L, 256L, NULL, 1e-10, TRUE, 1L, NULL, FALSE, "multiCSN", FALSE
  )
  expect_identical(positional$cells, in_memory$cells)
  expect_identical(positional$tf_gene, in_memory$tf_gene)
  expect_false(positional$settings$verbose)

  timepoint <- fit_layered_network(
    object, cell_group = "B", preset = "timepoint", min_detected = 5L, verbose = FALSE
  )
  expect_identical(timepoint$settings$min_detected, 5L)
  expect_identical(timepoint$settings$tf_region_scope, "motif_peaks")
  expect_true(timepoint$settings$renormalize)
  empty <- fit_layered_network(
    object, cell_group = "A", preset = "state", min_peak_cells = 1000L,
    verbose = FALSE
  )
  expect_true(all(empty$summary$selected == 0L))
  expect_identical(nrow(empty$accounting$tf_region), 0L)
  empty_motif <- fit_layered_network(
    object, cell_group = "A", tf_region_scope = "motif_peaks",
    min_peak_cells = 1000L, verbose = FALSE
  )
  expect_true(all(empty_motif$summary$selected == 0L))
  empty_path <- file.path(root, "empty")
  write_network_endpoints(empty, empty_path)
  expect_true(validate_network_endpoints(empty_path))
  constant_tfs <- object
  data <- SeuratObject::LayerData(constant_tfs, assay = "RNA", layer = "data")
  data[get_attribute(constant_tfs, attribute = "tfs"), ] <- 1
  constant_tfs <- SeuratObject::SetAssayData(constant_tfs, assay = "RNA", layer = "data", new.data = data)
  no_tf_fit <- fit_layered_network(
    constant_tfs, cell_group = "A", preset = "state", verbose = FALSE
  )
  expect_identical(nrow(no_tf_fit$tf_gene), 0L)
  expect_identical(nrow(no_tf_fit$tf_region), 0L)
  expect_identical(nrow(no_tf_fit$triplets), 0L)
  expect_true(nrow(no_tf_fit$region_gene) > 0L)
  expect_true(all(no_tf_fit$accounting$tf_gene$status == "no_variable_candidate"))
  expect_true(all(no_tf_fit$accounting$tf_gene$n_selected == 0L))
  no_tf_path <- file.path(root, "no-variable-tf")
  write_network_endpoints(no_tf_fit, no_tf_path)
  expect_true(validate_network_endpoints(no_tf_path))
  unlink(root, recursive = TRUE)
})

test_that("presets are validated", {
  expect_error(multiCSN:::layered_preset("weekly"), "should be one of")
  expect_identical(multiCSN:::layered_preset("timepoint")$min_detected, 20L)
})
