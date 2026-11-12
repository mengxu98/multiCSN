multiome_fixture <- function(seed = 7) {
  set.seed(seed)
  motif2tf <- multiCSN:::bundled_dataset("motif2tf")
  tfs <- sample(unique(motif2tf$tf), 25)
  targets <- paste0("GENE", 1:40)
  genes <- c(tfs, targets)
  cells <- paste0("cell", seq_len(240))
  chromosome <- sample(paste0("chr", 1:3), length(genes), TRUE)
  start <- sample(seq(2e6, 4e7, by = 2e5), length(genes))
  annotation <- GenomicRanges::GRanges(
    chromosome, IRanges::IRanges(start, width = 5000),
    strand = sample(c("+", "-"), length(genes), TRUE),
    gene_name = genes, gene_id = paste0("ENSG", seq_along(genes)),
    tx_id = paste0("ENST", seq_along(genes)),
    gene_biotype = "protein_coding", type = "exon"
  )
  GenomeInfoDb::genome(annotation) <- "hg38"
  peaks <- unique(unlist(lapply(seq_along(genes), function(i) {
    offset <- start[[i]] + sample(-60000:3000, 4)
    paste0(chromosome[[i]], "-", offset, "-", offset + 400)
  })))
  rna <- matrix(stats::rpois(length(genes) * length(cells), 2),
    length(genes),
    dimnames = list(genes, cells)
  )
  for (j in 1:15) {
    rna[targets[[j]], ] <- stats::rpois(length(cells), 1 + rna[tfs[[(j %% 5) + 1]], ])
  }
  atac <- matrix(stats::rbinom(length(peaks) * length(cells), 3, 0.3),
    length(peaks),
    dimnames = list(peaks, cells)
  )
  object <- SeuratObject::CreateSeuratObject(
    methods::as(rna, "dgCMatrix"),
    meta.data = data.frame(celltype = rep(c("A", "B"), each = 120), row.names = cells)
  )
  object[["ATAC"]] <- suppressWarnings(Signac::CreateChromatinAssay(
    methods::as(atac, "dgCMatrix"),
    ranges = parse_peak_ranges(peaks), annotation = annotation, genome = "hg38"
  ))
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  Signac::RunTFIDF(object, assay = "ATAC", verbose = FALSE)
}
