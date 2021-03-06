% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embed_cells.R
\name{embed_cells}
\alias{embed_cells}
\title{Embed Cells}
\usage{
embed_cells(
  input,
  lem,
  embedding_key = NULL,
  method = "tSNE",
  comp_select = NULL,
  tSNE_perp = 30,
  iterations = 500,
  write_colData = T,
  colData_labs = NULL,
  seed = 100,
  print_progress = T
)
}
\arguments{
\item{input}{the input sce}

\item{lem}{a LEM within reducedDims() must be provided}

\item{embedding_key}{The name of the embedding within metadata(reducedDim(input)) slot of SCE}

\item{method}{the method for embedding. tSNE or OTHER?}

\item{comp_select}{a vector of indices to use for embedding. NULL uses all, seq(1:10) would use first 10 sampleFactors withing the reducedDim(lem)}

\item{tSNE_perp}{number of cells expressed above threshold for a given gene, 10-100 recommended.}

\item{iterations}{The number of iterations for tSNE to perform.}

\item{write_colData}{whether or not to write the embedding to colData}

\item{colData_labs}{a vector of length 2 for the names of the X and Y colnames within colData()}

\item{seed}{For tSNE, the seed. Can set to NULL if desired.}

\item{print_progress}{will print progress if TRUE}
}
\description{
This function will embed cells using sampleFactors() within reducedDim()
}
\details{

}
\examples{

}
