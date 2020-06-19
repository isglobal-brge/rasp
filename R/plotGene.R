#' Plot splicing ratio boxplot by condition for each exon in a given gene
#'
#' @param data expression dataset. Must be of class `RangedSummarizedExperiment`.
#' @gene gene to be plotted.
#' @param condition a `factor` representing the condition of each individual.
#' @param title plot title. Optional.
#' @param ylab y-axis label. Optional.
#' @param xlab x-axix label. Optional.
#' @examples
#' data(adipose)
#' plotGene(adipose.chr16, gene="ENSG00000197165.10", condition="inv16p11.2")
#' 
#' @import ggplot2 
#' @import tibble tibble
#' @export
#' 

plotGene <- function(data, gene, condition = NULL,  title = NULL, 
                     ylab = "exon relative abundance", 
                     xlab="Exons"){

    sel <- which(rownames(data) %in% gene)
    if(length(sel)==0)
        stop("This gene is not in your RSE object")
    exonCounts <- SummarizedExperiment::assay(data[sel, ])
    exonNames <- 1:nrow(exonCounts)
    rownames(exonCounts) <- exonNames

    pheno <- SummarizedExperiment::colData(data)
    ii <- grep(condition, colnames(pheno))
    if (length(ii)!=1)
        stop("'condition' is not a variable in the metadata")
    condition <- pheno[, ii]
    if(!inherits(condition, "factor"))
        condition <- as.factor(condition)


    non.all.zero <- !apply(exonCounts, 2, function(x) all(x == 0))
    exonCounts <- exonCounts[, non.all.zero]
    condition <- condition[non.all.zero]
    xx <- apply(exonCounts, 2, function(y) y / sum(y))
    plotDat <- tibble(ratio = c(xx),
                          exon = rep(exonNames, ncol(xx)),
                          condition = rep(condition, each = nrow(xx)))
    gg <- ggplot(plotDat, aes(x = condition, y = ratio, fill=condition)) + 
        geom_boxplot() + 
        facet_grid(. ~ exon) + 
        labs(y=ylab, x=xlab) +  theme(axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank()) +
        ggtitle(title)
    gg
}
