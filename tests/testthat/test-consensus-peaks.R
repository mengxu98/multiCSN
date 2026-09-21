test_that("consensus peaks merge overlaps and adjacent intervals", {
  peaks <- list(
    c("chr1:100-200", "chr1:500-600", "chr2:100-150"),
    c("chr1:150-250", "chr1:601-700", "chr2:100-150")
  )
  consensus <- consensus_peak_set(peaks)
  expect_s4_class(consensus, "GRanges")
  expect_identical(names(consensus), c("chr1:100-250", "chr1:500-700", "chr2:100-150"))
  expect_identical(
    as.character(GenomeInfoDb::seqnames(consensus)), c("chr1", "chr1", "chr2")
  )
  limited <- consensus_peak_set(peaks, seqlevels = "chr2")
  expect_identical(names(limited), "chr2:100-150")
  seqinfo <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), seqlengths = c(1000, 1000))
  annotated <- consensus_peak_set(peaks, seqlevels = c("chr1", "chr2"), seqinfo = seqinfo)
  expect_identical(as.integer(GenomeInfoDb::seqlengths(annotated)), c(1000L, 1000L))
  expect_identical(length(annotated), 3L)
  expect_identical(
    as.character(GenomeInfoDb::seqlevels(annotated)), c("chr1", "chr2")
  )
  expect_error(
    consensus_peak_set(list("chr9:1-10"), seqlevels = "chr1"),
    "No sequence level"
  )
})

test_that("peak mapping and count projection are one-to-one", {
  consensus <- consensus_peak_set(list(c("chr1:100-200", "chr1:300-400")))
  index <- map_peaks_to_consensus(c("chr1:120-180", "chr1:301-399"), consensus)
  expect_identical(index, c(1L, 2L))
  expect_error(
    map_peaks_to_consensus("chr1:90-210", consensus),
    "one-to-one"
  )
  counts <- Matrix::Matrix(
    matrix(1:6, nrow = 3, dimnames = list(c("p1", "p2", "p3"), c("c1", "c2"))),
    sparse = TRUE
  )
  wide <- consensus_peak_set(list(c("chr1:100-200", "chr1:300-400")))
  expect_error(
    project_counts_to_consensus(counts, c("chr1:120-180", "chr1:90-210", "chr1:301-399"), wide),
    "one-to-one"
  )
  map <- list(c("chr1:100-150", "chr1:300-350"), c("chr1:140-200", "chr1:340-400"))
  wide <- consensus_peak_set(map)
  expect_identical(names(wide), c("chr1:100-200", "chr1:300-400"))
  counts <- Matrix::Matrix(
    matrix(c(1:4, 5:8), nrow = 4,
           dimnames = list(c("chr1:100-150", "chr1:300-350", "chr1:140-200", "chr1:340-400"),
                           c("c1", "c2"))),
    sparse = TRUE
  )
  projected <- project_counts_to_consensus(counts, rownames(counts), wide)
  expect_identical(rownames(projected), names(wide))
  expect_equal(
    as.matrix(projected),
    matrix(c(1 + 3, 2 + 4, 5 + 7, 6 + 8), nrow = 2,
           dimnames = list(names(wide), c("c1", "c2")))
  )
})
