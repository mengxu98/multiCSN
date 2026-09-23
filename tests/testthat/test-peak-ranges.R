test_that("parse_peak_ranges keeps every documented peak identifier form", {
  peaks <- c(
    "1:2-4", "chr1-5-8", "chr2_1_20",
    "chr1_KI270706v1_random:125391-126274",
    "GL000191-1-1000"
  )
  ranges <- parse_peak_ranges(peaks)
  expect_s4_class(ranges, "GRanges")
  expect_identical(names(ranges), peaks)
  expect_identical(
    as.character(GenomeInfoDb::seqnames(ranges)),
    c("1", "chr1", "chr2", "chr1_KI270706v1_random", "GL000191")
  )
  expect_identical(as.integer(IRanges::start(ranges)), c(2L, 5L, 1L, 125391L, 1L))
  expect_identical(as.integer(IRanges::end(ranges)), c(4L, 8L, 20L, 126274L, 1000L))
})

test_that("parse_peak_ranges treats the last delimiter as the coordinate boundary", {
  ranges <- parse_peak_ranges(c("chr1:2:3-4", "chr1-2-3-4"))
  expect_identical(as.character(GenomeInfoDb::seqnames(ranges)), c("chr1:2", "chr1-2"))
  expect_identical(as.integer(IRanges::start(ranges)), c(3L, 3L))
  expect_identical(as.integer(IRanges::end(ranges)), c(4L, 4L))
})

test_that("parse_peak_ranges validates coordinates and reports bad identifiers", {
  expect_error(parse_peak_ranges("chr1:5-2"), "end >= start")
  expect_error(parse_peak_ranges("chr1:0-2"), "positive integers")
  expect_error(parse_peak_ranges("chr1:5_8"), "Could not parse")
  expect_error(parse_peak_ranges("chr1:abc-def"), "Could not parse")
  expect_error(
    parse_peak_ranges(c("chr1:1-2", "chr1:5-2")),
    "chr1:5-2"
  )
  empty <- parse_peak_ranges(character(0))
  expect_identical(length(empty), 0L)
})

test_that("the split fast path and the pattern fallback agree", {
  # coordinates beyond the integer range keep the coordinate error, not the
  # "could not parse" error
  expect_error(parse_peak_ranges("chr1:21474836470-21474836480"), "finite positive")
  # forms that only look like an identifier must reach the same error
  expect_error(parse_peak_ranges(c("chr1:1-2", "chr1:5_8")), "Could not parse")
  expect_error(parse_peak_ranges(c("chr1-5_2", "scaffold-7:0:0")), "Could not parse")
  # several delimiters in one sequence name resolve to the last delimiter
  ranges <- parse_peak_ranges(c(
    "chr1_KI270706v1_random:125391-126274", "chrUn-1_2-3-4", "GL000191-1-1000"
  ))
  expect_identical(
    as.character(GenomeInfoDb::seqnames(ranges)),
    c("chr1_KI270706v1_random", "chrUn-1_2", "GL000191")
  )
  expect_identical(as.integer(IRanges::start(ranges)), c(125391L, 3L, 1L))
  expect_identical(as.integer(IRanges::end(ranges)), c(126274L, 4L, 1000L))
  # a large mixed batch stays consistent with the per-item results
  batch <- c("chr1:1-2", "chr2_3_4", "chr3-5-6", "chr4:7:8-9")
  expect_identical(
    as.character(GenomeInfoDb::seqnames(parse_peak_ranges(batch))),
    c("chr1", "chr2", "chr3", "chr4:7")
  )
})

test_that("parse_peak_ranges is vectorised across many peaks", {
  set.seed(3)
  starts <- sample(1:1e6, 5000L, TRUE)
  peaks <- paste0(
    "chr", sample(1:22, 5000L, TRUE), ":",
    starts, "-", starts + sample(1:999, 5000L, TRUE)
  )
  ranges <- parse_peak_ranges(peaks)
  expect_identical(length(ranges), 5000L)
  expect_identical(as.character(GenomeInfoDb::seqnames(ranges)), sub(":.*$", "", peaks))
  expect_true(all(as.integer(IRanges::end(ranges)) >= as.integer(IRanges::start(ranges))))
})

test_that("one unusually delimited peak does not disable the fast path", {
  ordinary <- paste0("chr1:", seq_len(100L), "-", seq_len(100L) + 1L)
  mixed <- c(ordinary, "chr1:extra:101-102")
  parsed <- multiCSN:::.split_peak_identifiers(mixed)
  expect_false(anyNA(parsed$chromosome))
  expect_identical(parsed$chromosome, c(rep("chr1", 100L), "chr1:extra"))
  expect_identical(as.integer(IRanges::start(parse_peak_ranges(mixed))),
                   c(seq_len(100L), 101L))
})

test_that("mixed delimiter batches match the original peak grammar", {
  set.seed(7)
  chromosomes <- sample(c("chr1", "chrUn-1", "chrA_2", "chrB:3"), 400L, TRUE)
  starts <- sample.int(100000L, 400L)
  ends <- starts + sample.int(1000L, 400L, TRUE)
  separators <- sample(c(":", "-", "_"), 400L, TRUE)
  peaks <- ifelse(
    separators == ":", paste0(chromosomes, ":", starts, "-", ends),
    paste0(chromosomes, separators, starts, separators, ends)
  )
  reference <- t(vapply(peaks, function(peak) {
    for (pattern in multiCSN:::.peak_identifier_patterns) {
      matched <- regmatches(peak, regexec(pattern, peak, perl = TRUE))[[1L]]
      if (length(matched) == 4L) return(matched[2:4])
    }
    stop("The fixture contains an invalid peak identifier")
  }, character(3)))
  parsed <- parse_peak_ranges(peaks)
  expect_identical(as.character(GenomeInfoDb::seqnames(parsed)), unname(reference[, 1L]))
  expect_identical(as.integer(IRanges::start(parsed)), as.integer(reference[, 2L]))
  expect_identical(as.integer(IRanges::end(parsed)), as.integer(reference[, 3L]))
})
