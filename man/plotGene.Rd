% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGene.R
\name{plotGene}
\alias{plotGene}
\title{Plot splicing ratio boxplot by condition for each exon in a given gene}
\usage{
plotGene(
  data,
  gene,
  condition = NULL,
  title = NULL,
  ylab = "exon relative abundance",
  xlab = "Exons"
)
}
\arguments{
\item{data}{expression dataset. Must be of class `RangedSummarizedExperiment`.}

\item{condition}{a `factor` representing the condition of each individual.}

\item{title}{plot title. Optional.}

\item{ylab}{y-axis label. Optional.}

\item{xlab}{x-axix label. Optional.}
}
\description{
Plot splicing ratio boxplot by condition for each exon in a given gene
}
\examples{
data(adipose)
plotGene(adipose.chr16, gene="ENSG00000197165.10", condition="inv16p11.2")

}
