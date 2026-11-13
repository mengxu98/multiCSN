test_that("checkpoint hashes use the existing SHA-256 byte contract", {
  object <- list(value = c(NA_real_, Inf, pi), text = "基因", matrix = Matrix::Diagonal(3))
  bytes <- serialize(object, NULL, xdr = FALSE, version = 3L)
  expected <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  expect_identical(multiCSN:::layered_sha256(object), expected)
  if ("sha256sum" %in% getNamespaceExports("tools")) {
    expect_identical(expected, unname(tools::sha256sum(bytes = bytes)))
  }
  root <- tempfile("package-hashes-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  writeBin(as.raw(0:255), file.path(root, "DESCRIPTION"))
  file.create(file.path(root, "NAMESPACE"))
  expected_files <- vapply(
    c("DESCRIPTION", "NAMESPACE"),
    function(path) digest::digest(file = file.path(root, path), algo = "sha256"),
    character(1)
  )
  expect_identical(multiCSN:::layered_package_identity(root), expected_files)
})

test_that("explicit peak separators preserve literal coordinate boundaries", {
  examples <- list(
    list(peaks = c("chr1:100:200", "chr2:300:400"), sep = c(":", ":")),
    list(peaks = c("chr1:100-200", "chr2:300-400"), sep = c(":", "-")),
    list(peaks = c("chr1.100|200", "chr2.300|400"), sep = c(".", "|")),
    list(peaks = c("chr1::100++200", "chr2::300++400"), sep = c("::", "++")),
    list(peaks = c("chrUn-1-100-200", "chr2-300-400"), sep = c("-", "-"))
  )
  for (example in examples) {
    ranges <- parse_peak_ranges(example$peaks, sep = example$sep)
    expect_identical(names(ranges), example$peaks)
    expect_identical(as.integer(IRanges::start(ranges)), c(100L, 300L))
    expect_identical(as.integer(IRanges::end(ranges)), c(200L, 400L))
  }
  contig <- parse_peak_ranges("chr1:part:100:200", sep = c(":", ":"))
  expect_identical(as.character(GenomeInfoDb::seqnames(contig)), "chr1:part")
  expect_error(parse_peak_ranges("chr1:100:200", sep = c(":", "")), "sep")
  expect_error(parse_peak_ranges("chr1:100:200", sep = ":"), "sep")
  expect_error(parse_peak_ranges("chr1:100:200", sep = c(":", NA)), "sep")
  expect_error(parse_peak_ranges("chr1:100:200", sep = c("-", "-")), "Could not parse")
  expect_error(parse_peak_ranges("chr1:200:100", sep = c(":", ":")), "end >= start")
})

test_that("normalizing a plain peak assay honors sep and preserves counts", {
  counts <- Matrix::Matrix(
    matrix(c(1, 3, 2, 4, 5, 2, 3, 1), nrow = 2,
      dimnames = list(c("chr1:100:200", "chr2:300:400"), paste0("c", 1:4))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts, assay = "ATAC")
  object <- SeuratObject::SetAssayData(object, assay = "ATAC", layer = "data", new.data = counts / 2)
  normalized <- suppressWarnings(normalize_chromatin_assay(object, peak_assay = "ATAC", sep = c(":", ":")))
  ranges <- Signac::granges(normalized[["ATAC"]])
  expect_identical(as.character(GenomeInfoDb::seqnames(ranges)), c("chr1", "chr2"))
  expect_identical(as.integer(IRanges::start(ranges)), c(100L, 300L))
  expect_identical(as.integer(IRanges::end(ranges)), c(200L, 400L))
  actual <- SeuratObject::LayerData(normalized, assay = "ATAC", layer = "counts")
  expect_identical(colnames(actual), colnames(counts))
  expect_equal(unname(as.matrix(actual)), unname(as.matrix(counts)))
  actual_data <- SeuratObject::LayerData(normalized, assay = "ATAC", layer = "data")
  expect_equal(unname(as.matrix(actual_data)), unname(as.matrix(counts / 2)))
  automatic_counts <- counts
  rownames(automatic_counts) <- c("chr1:100-200", "chr2:300-400")
  automatic_object <- SeuratObject::CreateSeuratObject(automatic_counts, assay = "ATAC")
  automatic <- suppressWarnings(normalize_chromatin_assay(automatic_object, peak_assay = "ATAC"))
  expect_equal(SeuratObject::LayerData(automatic, assay = "ATAC", layer = "counts"), actual)
})
