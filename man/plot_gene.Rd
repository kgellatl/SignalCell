% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gene.R
\name{plot_gene}
\alias{plot_gene}
\title{Plot on a gene or gene}
\usage{
plot_gene(
  input,
  gene,
  facet_by = NULL,
  assay = NULL,
  coords = NULL,
  title = NULL,
  colors = c("gray", "blue", "red", "yellow"),
  size = 1.5,
  ncol = 2,
  text_sizes = c(20, 10, 5, 10, 5, 5),
  xlab = "x",
  ylab = "y"
)
}
\arguments{
\item{input}{The input data}

\item{gene}{the gene or genes to plot}

\item{facet_by}{If ONE gene is being plotted it can be faceted. You cannot facet multiple gene}

\item{assay}{The assay to operate on. Will default to get_def_assay(input)}

\item{coords}{A vector of character lem, x, and u. See set_xy function.}

\item{title}{The title}

\item{colors}{Colors for the cells}

\item{size}{The size of the points}

\item{ncol}{controls the number of columns if faceting}

\item{text_sizes}{the text sizes on the plot}

\item{xlab}{character x label}

\item{ylab}{character y label}

\item{theme}{the ggplot theme}
}
\description{
This will plot gene information onto a 2d  plot
}
\details{
Utilize information stored in colData to control the plot display.
}
\examples{
plot_tsne_gene(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)


}
