test_that("deletion evidence determines ordinal weights and chain projections", {
  expect_equal(
    network_ordinal_weight(c(1, -2, 0, NA, 1), c(5, 5, 100, 1, 0)),
    c(0.75, -0.75, 0, 0, 0.25)
  )
  tr <- data.frame(
    regulator = c("A", "A", "B"), region = c("r1", "r2", "r1"),
    standardized_beta = c(1, -1, 1), deletion_delta_bic = c(6, 6, 4)
  )
  rg <- data.frame(
    region = c("r1", "r2"), target = "T",
    standardized_beta = 1, deletion_delta_bic = 5
  )
  tg <- data.frame(
    regulator = c("A", "B"), target = "T",
    standardized_beta = 1, deletion_delta_bic = 3
  )
  chains <- factorized_chain_support(tr, rg, tg)
  expect_equal(chains$chain_delta_bic, c(5, 5, 4))
  expect_equal(chains$weight, c(0.75, -0.75, 0.25))
  expect_identical(chains$backbone_sign_concordant, c(TRUE, FALSE, TRUE))
  projected <- strongest_mediated_projection(chains)
  expect_identical(projected$regulator, "B")
  expect_equal(projected$weight, 0.5)
  expect_error(factorized_chain_support(rbind(tr, tr[1, ]), rg), "duplicated")
  expect_equal(nrow(factorized_chain_support(tr[FALSE, ], rg)), 0L)
})

test_that("peak adapters preserve coordinates and collapse repeated annotations", {
  peaks <- parse_peak_ranges(c("1:2-4", "chr1-5-8", "chr2_1_20"))
  expect_equal(as.integer(IRanges::start(peaks)), c(2L, 5L, 1L))
  genome <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), seqlengths = c(10, 10))
  mapped <- map_peak_ranges_to_genome(peaks, genome)
  expect_identical(mapped$mapped_seqnames, c("chr1", "chr1", "chr2"))
  expect_identical(mapped$supported, c(TRUE, TRUE, FALSE))
  expect_identical(names(mapped$ranges), names(peaks))
  expect_error(parse_peak_ranges("chr1:5-2"), "coordinates")
  x <- Matrix::sparseMatrix(
    i = c(1, 1, 2), j = 1:3, x = c(1, 2, 3),
    dims = c(2, 3), dimnames = list(c("p1", "p2"), c("g1", "g1", "g2"))
  )
  collapsed <- collapse_peak_gene_domains(x)
  expect_equal(as.matrix(collapsed), matrix(c(3, 0, 0, 3), 2, 2,
    dimnames = list(c("g1", "g2"), c("p1", "p2"))
  ))
})

test_that("Seurat layer and annotation adapters preserve their inputs", {
  x <- Matrix::Matrix(matrix(1:12, 3, 4, dimnames = list(paste0("g", 1:3), paste0("c", 1:4))), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(x)
  expect_equal(
    get_layer_data(object, assay = "RNA", layer = "counts"),
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  )
  annotation <- parse_peak_ranges("chr1:2-4")
  object@misc$multiCSN_annotations <- list(ATAC = annotation)
  expect_equal(get_peak_annotation(object, "ATAC"), annotation)
})
