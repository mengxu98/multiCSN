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

  timepoint <- fit_layered_network(
    object, cell_group = "B", preset = "timepoint", min_detected = 5L, verbose = FALSE
  )
  expect_identical(timepoint$settings$min_detected, 5L)
  expect_identical(timepoint$settings$tf_region_scope, "motif_peaks")
  expect_true(timepoint$settings$renormalize)
  unlink(root, recursive = TRUE)
})

test_that("presets are validated", {
  expect_error(multiCSN:::layered_preset("weekly"), "should be one of")
  expect_identical(multiCSN:::layered_preset("timepoint")$min_detected, 20L)
})
