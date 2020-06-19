#' Plot splicing ratios of exon expression stored in a RangeSummarizedExperiment object
#'
#' @param data expression dataset. Must be of class `RangedSummarizedExperiment`.
#' @param gene name of the gene to plot.
#' @param condition a `factor` representing the condition of each individual.
#' @param jitterWidth amount of graphical jitter to use in the plot.
#' @param inds `character` or `numeric` vector. Plot only the specified individuals (samples).
#' @param exons `character` or `numeric` vector. Plot only the specified exons.
#' @examples
#' data(adipose)
#' plotSplicingRatios("ENSG00000003249.13", adipose.chr16, "inv16p11.2")
#' 
#' @import ggplot2 
#' @import tibble tibble
#' @export
#' 

plotSplicingRatios <- function(data, gene, condition, jitterWidth = 0, inds, exons) {
    sel <- which(rownames(data) %in% gene)
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
   
    if (!missing(exons)) {
        if (any(!exons%in%exonNames))
            stop("This gene has not these exons")
        else
            exonCounts <- exonCounts[exons, ]
    }
    if (!missing(inds)) {
        if (is.character(inds)) {
            exonCounts <- exonCounts[, colnames(exonCounts) %in% inds]
            condition <- condition[colnames(exonCounts) %in% inds]
        }
        else {
            exonCounts <- exonCounts[, inds]
            condition <- condition[inds]
        }
    }
    
    non.all.zero <- !apply(exonCounts, 2, function(x) all(x == 0))
    exonCounts <- exonCounts[, non.all.zero]
    condition <- condition[non.all.zero]
    xx <- apply(exonCounts, 2, function(y) y / sum(y))
    plotDat <- tibble(ratio = c(xx),
                          exon = rep(rownames(xx), ncol(xx)),
                          ind = rep(colnames(xx), each = nrow(xx)),
                          condition = rep(condition, each = nrow(xx)))

    gg <- ggplot(data = plotDat, aes(x = ind, y = ratio, col = exon)) +
        geom_point(size = 1, shape = 16, position = position_jitter(height = 0, 
                                                                      width = jitterWidth)) +
        facet_grid(~ condition) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "gray93"),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank()) +
             scale_color_discrete("Exons") +
        ggtitle(gene) + ylab("Splicing Ratio") + xlab("Condition")
    
    gg
}
