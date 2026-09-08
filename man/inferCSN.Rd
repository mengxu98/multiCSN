% Please edit documentation in R/inferCSN.R, R/inferCSN-methods.R
\docType{methods}
\name{inferCSN}
\alias{inferCSN}
\alias{inferCSN,Network-method}
\alias{inferCSN,Seurat-method}
\alias{inferCSN,CSNObject-method}
\title{inferring cell-type specific gene regulatory network}
\usage{
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

\S4method{inferCSN}{Network}(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  r_squared_threshold = 0,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  gene_cor_threshold = 0,
  ...
)

\S4method{inferCSN}{Seurat}(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  r_squared_threshold = 0,
  celltypes = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  alpha = 0.5,
  family = "gaussian",
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  pseudotime_column = NULL,
  dynamic_features = NULL,
  ...
)

\S4method{inferCSN}{CSNObject}(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{Expression matrix, \code{Network}, \code{Seurat}, or \code{CSNObject}.
Matrix input is handled by \code{\link[inferCSN:inferCSN]{inferCSN::inferCSN()}}.}

\item{pseudotime}{Optional pseudotime vector or branch matrix for matrix input.
For \code{Seurat} objects, a single column name is treated as \code{pseudotime_column}.}

\item{regulators, targets}{Optional gene subsets.}

\item{max_support_size}{Optional support-size limit for matrix inference.}

\item{lag_fraction}{Fractional lag used when \code{lag_steps} is \code{NULL}.}

\item{lag_steps}{Optional integer lag for matrix inference.}

\item{cores}{Number of workers.}

\item{verbose}{Whether to report progress.}

\item{...}{Additional method arguments.}

\item{penalty, cross_validation, seed, n_folds}{Stored on \code{Network} objects for provenance. Fitting uses \code{method}.}

\item{r_squared_threshold}{Threshold of \eqn{R^2} used when formatting \code{Network} / \code{Seurat} results.}

\item{method}{Fitting backend for \code{Network} / \code{Seurat} methods:
\code{"greedy_l0"} (greedy L0 via \code{\link[inferCSN:fit_greedy_l0]{inferCSN::fit_greedy_l0()}}),
\code{"glm"}, \code{"glmnet"}, \code{"cv.glmnet"}, \code{"xgb"}, or \code{"susie"}.}

\item{gene_cor_threshold}{Threshold for TF-target gene correlation.}

\item{celltypes}{Character vector of cell types to infer networks for.}

\item{network_name}{Name of the stored network.}

\item{peak_to_gene_method}{Method to link peaks to genes. One of \code{Signac} or \code{GREAT}.}

\item{upstream}{Distance upstream of the gene to consider as a potential regulatory region.}

\item{downstream}{Distance downstream of the gene to consider as a potential regulatory region.}

\item{extend}{Distance added to the basal regulatory region when \code{peak_to_gene_method = "GREAT"}.}

\item{only_tss}{Measure distance from the TSS (\code{TRUE}) or from the gene body (\code{FALSE}).}

\item{peak_to_gene_domains}{\code{GenomicRanges} object with regulatory regions for each gene.}

\item{peak_cor_threshold}{Threshold for binding peak-target gene correlation.}

\item{aggregate_rna_col}{Metadata column used to aggregate RNA.}

\item{aggregate_peaks_col}{Metadata column used to aggregate peaks.}

\item{alpha}{Elastic-net mixing parameter. See \code{\link[glmnet]{glmnet}}.}

\item{family}{Error distribution and link function. See \code{\link[stats]{family}}.}

\item{interaction_term}{Interaction between TF and binding site: \code{"+"}, \code{":"}, or \code{"*"}.}

\item{adjust_method}{Method for adjusting p-values.}

\item{scale}{Whether to z-transform expression and accessibility matrices.}

\item{pseudotime_column}{Metadata column for pseudotime. If set, runs state-dynamic inference.}

\item{dynamic_features}{Optional dynamic-gene filter when \code{pseudotime_column} is set.}
}
\value{
A regulator-target edge table for matrix input, a \code{Network} object, or a \code{Seurat} object.
}
\details{
The Seurat method with chromatin peaks fits TF-by-region interaction predictors
through \code{fit_models}; projecting those coefficients does not fit independent
TF-gene, TF-region and region-gene models. For factorized analyses, fit the three
endpoint responses independently and use \code{factorized_chain_support} to
assemble their selected supports.
}
\description{
Infers cell-specific gene regulatory networks. Matrix methods come from
\pkg{inferCSN}; \code{Network} and \code{Seurat} methods are implemented here.
}
\examples{
data("example_matrix", package = "inferCSN")
network_table_1 <- inferCSN(
  example_matrix
)

network_table_2 <- inferCSN(
  example_matrix,
  cores = 2
)

head(network_table_1)

identical(
  network_table_1,
  network_table_2
)

inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g0")
)

\dontrun{
data("example_ground_truth", package = "inferCSN")
calculate_metrics(
  network_table_1,
  example_ground_truth,
  return_plot = TRUE
)
}
\dontrun{
data("example_matrix", package = "inferCSN")
object <- initiate_object(example_matrix)
object <- inferCSN(object)
}
\dontrun{
# multiome input: motif scanning comes first
data(pbmcmultiome_sub, package = "scop")
data(motifs)
data(motif2tf)

object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
object <- initiate_object(
  object,
  group.by = "CellType",
  rna_assay = "RNA",
  peak_assay = "peaks",
  verbose = FALSE
)

genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
object <- find_motifs(
  object,
  pfm = motifs,
  motif_tfs = motif2tf,
  genome = genome,
  backend = "motifmatchr",
  verbose = FALSE
)

object <- inferCSN(object, cores = 2, verbose = FALSE)
networks <- export_csn(object)
head(networks[[1]])
}
}
