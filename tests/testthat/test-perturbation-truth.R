make_screen <- function(replicates = paste0("rep", 1:3), targets = c("A", "B"),
                        cells_per_group = 40L, features = 60L, seed = 5) {
  set.seed(seed)
  groups <- expand.grid(
    replicate = replicates, target = c("NT", targets), stringsAsFactors = FALSE
  )
  cell_ids <- unlist(lapply(seq_len(nrow(groups)), function(i) {
    paste(groups$replicate[[i]], groups$target[[i]], seq_len(cells_per_group), sep = "_")
  }))
  sample_info <- data.frame(
    cell_id = cell_ids,
    replicate = rep(groups$replicate, each = cells_per_group),
    target = rep(groups$target, each = cells_per_group),
    timepoint = "day 7",
    stringsAsFactors = FALSE
  )
  counts <- matrix(
    stats::rnbinom(length(cell_ids) * features, mu = 5, size = 2),
    nrow = features, dimnames = list(paste0("f", seq_len(features)), cell_ids)
  )
  signal <- sample_info$target == "A"
  counts[1, signal] <- counts[1, signal] + 60
  list(counts = Matrix::Matrix(counts, sparse = TRUE), sample_info = sample_info)
}

test_that("pseudobulk aggregation sums cells into paired samples", {
  screen <- make_screen()
  aggregated <- pseudobulk_perturbation(
    list(RNA = screen$counts), screen$sample_info, targets = c("A", "B"), verbose = FALSE
  )
  metadata <- aggregated$sample_metadata
  expect_identical(nrow(metadata), 9L)
  expect_identical(
    metadata$sample_id,
    c(
      "rep1__NT", "rep1__A", "rep1__B", "rep2__NT", "rep2__A", "rep2__B",
      "rep3__NT", "rep3__A", "rep3__B"
    )
  )
  expect_identical(unique(metadata$cells), 40L)
  expect_identical(
    unique(metadata$role),
    c("network_training_control", "heldout_perturbation")
  )
  expect_identical(colnames(aggregated$pseudobulk$RNA), metadata$sample_id)
  # column sums equal totals of the contributing cells
  raw <- screen$counts[, screen$sample_info$cell_id, drop = FALSE]
  expect_equal(
    as.numeric(Matrix::colSums(aggregated$pseudobulk$RNA[, "rep1__A", drop = FALSE])),
    sum(raw[, screen$sample_info$replicate == "rep1" & screen$sample_info$target == "A"])
  )
  expect_equal(
    as.numeric(aggregated$pseudobulk$RNA[1, "rep2__A"]),
    sum(raw[1, screen$sample_info$replicate == "rep2" & screen$sample_info$target == "A"])
  )
})

test_that("pseudobulk aggregation validates its inputs", {
  screen <- make_screen()
  expect_error(
    pseudobulk_perturbation(screen$counts, screen$sample_info[, c("cell_id", "replicate")],
                            verbose = FALSE),
    "sample_info must have columns"
  )
  missing_cells <- screen$sample_info
  expect_error(
    pseudobulk_perturbation(
      list(RNA = screen$counts[, -1]), missing_cells, targets = c("A", "B"), verbose = FALSE
    ),
    "lacks 1 annotated cell"
  )
  duplicated_cells <- rbind(screen$sample_info, screen$sample_info[1, , drop = FALSE])
  expect_error(
    pseudobulk_perturbation(
      list(RNA = screen$counts), duplicated_cells, targets = c("A", "B"), verbose = FALSE
    ),
    "duplicated cell identifiers"
  )
  incomplete <- screen$sample_info[!(screen$sample_info$replicate == "rep3" &
    screen$sample_info$target == "B"), , drop = FALSE]
  expect_error(
    pseudobulk_perturbation(
      list(RNA = screen$counts[, incomplete$cell_id]), incomplete,
      targets = c("A", "B"), verbose = FALSE
    ),
    "lacks a sample for: B"
  )
  # one replicate carrying two timepoints is rejected
  mixed_timepoint <- screen$sample_info
  rep2 <- mixed_timepoint$replicate == "rep2"
  renamed <- which(rep2)[seq_len(sum(rep2) %/% 2L)]
  original_ids <- mixed_timepoint$cell_id[renamed]
  mixed_timepoint$cell_id[renamed] <- paste0(original_ids, "_d9")
  mixed_timepoint$timepoint[renamed] <- "day 9"
  extra <- screen$counts[, original_ids, drop = FALSE]
  colnames(extra) <- mixed_timepoint$cell_id[renamed]
  mixed_counts <- cbind(screen$counts, extra)
  expect_error(
    pseudobulk_perturbation(
      list(RNA = mixed_counts), mixed_timepoint, targets = c("A", "B"), verbose = FALSE
    ),
    "maps to several timepoints"
  )
})

test_that("paired effect fitting recovers the spiked perturbation effect", {
  skip_if_not_installed("edgeR")
  screen <- make_screen()
  aggregated <- pseudobulk_perturbation(
    list(RNA = screen$counts), screen$sample_info, targets = c("A", "B"), verbose = FALSE
  )
  output <- file.path(tempdir(), "perturbation-effects")
  unlink(output, recursive = TRUE)
  fitted <- fit_perturbation_effects(
    aggregated$pseudobulk, aggregated$sample_metadata, targets = c("A", "B"),
    output_dir = output, verbose = FALSE
  )
  expect_identical(nrow(fitted$summary), 2L)
  expect_identical(
    sort(fitted$summary$target), c("A", "B")
  )
  expect_identical(unique(fitted$summary$paired_replicates), 3L)
  expect_identical(unique(fitted$summary$control_cells), 120L)
  expect_identical(unique(fitted$summary$perturbation_cells), 120L)
  files <- list.files(output, pattern = "[.]tsv[.]gz$", recursive = TRUE)
  expect_identical(length(files), 2L)
  effect <- data.table::fread(file.path(output, "RNA", "day_7__A.tsv.gz"))
  expect_identical(names(effect)[1], "feature")
  expect_identical(unique(effect$target), "A")
  expect_identical(unique(effect$effect_definition), "heldout_perturbed_minus_NT_control")
  expect_true(effect$logFC[effect$feature == "f1"] > 1)
  expect_true(effect$FDR[effect$feature == "f1"] < 0.05)
  expect_gt(fitted$summary$significant_up_fdr_0_05[fitted$summary$target == "A"], 0L)
  # the same call without output_dir keeps the tables in memory
  in_memory <- fit_perturbation_effects(
    aggregated$pseudobulk, aggregated$sample_metadata, targets = "A", verbose = FALSE
  )
  expect_true(is.data.frame(in_memory$effects$RNA[["day_7__A"]]))
  expect_identical(in_memory$summary$output_file, NA_character_)
})
