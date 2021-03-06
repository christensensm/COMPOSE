% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.volcano.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Volcano plot for sgRNA fold-change results}
\usage{
plot_volcano(
  results,
  print = T,
  save = F,
  sig.cutoff = 0.05,
  contrast = NULL,
  top = NULL
)
}
\arguments{
\item{results}{a calc.DESeq2.L2FC object}

\item{print}{logical - do you want to print the  plot to plots}

\item{save}{logical - do you want to save the  plot to pdf}

\item{sig.cutoff}{adjusted p-value cutoff for sgRNA fold-changes}

\item{top}{numerical - number of top guides to label (by p-value)}
}
\description{
Visualization of fold-changes and p-values for sgRNAs
}
\examples{

...
}
\seealso{
\code{\link{ggplot}, \link{ggrepel}} which this function uses to plot
}
