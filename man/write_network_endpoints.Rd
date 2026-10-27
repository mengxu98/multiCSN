\name{write_network_endpoints}
\alias{write_network_endpoints}
\title{Write layered network result tables}
\usage{
write_network_endpoints(
  fit,
  root,
  format = c("tsv", "csv", "tsv.gz", "csv.gz"),
  overwrite = FALSE
)
}
\arguments{
\item{fit}{The result returned by \code{fit_layered_network()}.}
\item{root}{Output directory.}
\item{format}{File format: tsv, csv, tsv.gz or csv.gz.}
\item{overwrite}{Replace existing output files. Defaults to \code{FALSE}.}
}
\description{
Writes the five network endpoint tables, three per-target accounting tables
and the endpoint summary. TF-region accounting uses \code{target} as the
response column name. No Seurat object or checkpoint files are written.
}
\value{
An invisible named character vector of output paths.
}
\details{
Use the same \code{format} with \code{read_network_endpoints()} to read the
five network tables. Accounting and summary tables can be read as ordinary
delimited tables. Empty networks retain their column headers.
}
