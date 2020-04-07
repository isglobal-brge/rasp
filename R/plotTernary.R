#' Plot proportional exon expression
#'
#' @param x expression dataset. Can be of class `DEXSeqDataSet`, `SummarizedExperiment` or
#'  `RangedSummarizedExperiment`. A `matrix` or `data.frame` is also accepted.
#' @param gene name of the gene to plot.
#' @param group a `factor` representing the condition of each individual.
#'   Required if `data` is a `matrix` or `data.fame`, ignored otherwise.
#' @param geneCol `character` or `numeric` vector specifying the columns of `x` with that contain the gene names.
#'   Required if `x` is a `matrix` or `data.fame`, ignored otherwise.
#' @param transCol `character` or `numeric` vector specifying the columns of `x` with that contain the transcript names.
#'   Required if `x` is a `matrix` or `data.fame`, ignored otherwise.
#' @param transcripts `numeric` vector of length 3 with the indices of the transcripts to ploth.
#' @examples
#' data(adipose)
#' plotTernary(adipose.chr16, "ENSG00000003249.13", transcripts = c(2, 4, 7))
#' 
#' data(YRI)
#' plotTernary(YRI, "ENSG00000160741")
# @import SummarizedExperiment
# @import DirichletReg
#' @export
plotTernary <- function(x, gene, group, geneCol = 1, transCol = 2, transcripts, ...) {
  if (missing(transcripts)) transcripts <- 1:3
  else if (!inherits(transcripts, "numeric")) stop("'transcripts' must be a numeric index")
  else if (length(transcripts) != 3) stop("Please specify 3 transcripts")

  if (inherits(x, "RangedSummarizedExperiment")) {
    xx <- SummarizedExperiment::assay(x)[grep(gene, rownames(x)), ]
  } else {
    xx <- x[x[, 1] == gene, ]
    rownames(xx) <- xx[, transCol] 
    xx <- xx[, -c(geneCol, transCol)]
  }
  
  if (nrow(xx) < 3) stop("A 'gene' must have at least 3 transcripts to be plotted")

  ratio.dr <- DirichletReg::DR_data(t(xx))
  if (missing(group))
    DirichletReg::plot.DirichletRegData(ratio.dr, dims = transcripts, ...)
  else {
    if (length(group) != ncol(xx)) stop("'group' length does not match the number of samples in 'x'")
    if (nlevels(group) > 5) warning("Too many levels in 'group'. Only the first 5 will be plotted")
    
    mycol <- c("red", "blue", "darkgreen", "yellow", "orange")[1:nlevels(group)]
    ok <- as.character(factor(group, labels = mycol))
    DirichletReg::plot.DirichletRegData(ratio.dr, col = ok, dims = transcripts, ...)
    legend("topleft", levels(group), col = mycol, pch = 16)
  }       
}