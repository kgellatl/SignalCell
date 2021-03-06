% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_violin.R
\name{plot_violin}
\alias{plot_violin}
\title{Plot Violin}
\usage{
plot_violin(
  input,
  gene,
  color_by,
  facet_by = NULL,
  assay = NULL,
  title = NULL,
  jitter = T,
  plot_mean = T,
  number_labels = T,
  sig = 2,
  size = 1,
  ncol = 2,
  text_sizes = c(20, 10, 5, 10, 5, 5)
)
}
\arguments{
\item{input}{The input data}

\item{gene}{the gene or genes to plot}

\item{color_by}{colData() column to facet plot by}

\item{facet_by}{If ONE gene is being plotted it can be faceted. You cannot facet multiple genes}

\item{assay}{The assay to operate on. Will default to get_def_assay(input)}

\item{title}{The title}

\item{jitter}{Whether or not to add jitter points to violin}

\item{plot_mean}{whether or not to plot the mean value along with a secondary axis}

\item{number_labels}{Add counts and fraction expression to plot}

\item{sig}{number of decimals for number_labels}

\item{size}{The size of the points}

\item{ncol}{controls the number of columns if faceting}

\item{text_sizes}{the text sizes on the plot}
}
\description{
This will plot gene expression via violin plot
}
\details{
Utilize information stored in colData to control the plot display.
}
\examples{
plot_violin(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)


}
