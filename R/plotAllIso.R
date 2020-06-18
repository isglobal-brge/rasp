#' Plot proportional exon expression
#'
#' @param gene name of the gene to plot.
#' @param data expression dataset. Can be of class `DEXSeqDataSet`, `SummarizedExperiment` or
#'  `RangedSummarizedExperiment`. A `matrix` or `data.frame` is also accepted.
#' @param group a `factor` representing the condition of each individual.
#'   Required if `data` is a `matrix` or `data.fame`, ignored otherwise.
#' @param gene_col name of the column that contains the gene IDs.
#' @param jitterWidth amount of graphical jitter to use in the plot.
#' @param inds `character` or `numeric` vector. Plot only the specified individuals (samples).
#' @param isos `character` or `numeric` vector. Plot only the specified exons.
#' @examples
#' data(adipose)
#' plotAllIso("ENSG00000003249.13", adipose.chr16, "inv16p11.2")
#' 
#' data(YRI)
#' plotAllIso("ENSG00000215915", YRI, attr(YRI, "gender"))
#' 
#' @import ggplot2
#' @export
#'  

plotAllIso <- function(gene, data, group, gene_col = 'gene_id', jitterWidth = 0, inds, isos) {
    if (is.data.frame(data))
        x <- data[data[[gene_col]] == gene, ]
    else if (inherits(data, "RangedSummarizedExperiment") |
             inherits(data, "SummarizedExperiment")) {
        x <- SummarizedExperiment::assay(data)[rownames(data) == gene, ]
        group <- SummarizedExperiment::colData(data)[, group]
    }
    
    if (!missing(isos)) {
        if (is.character(isos))
            x <- x[x$transcript_id %in% isos, ]
        else
            x <- x[isos, ]
    }
    
    if(!missing(inds)) {
        if(is.character(inds)) {
            x <- x[, colnames(x) %in% inds]
            group <- group[colnames(x) %in% inds]
        } else {
            x <- x[,c(1,2,inds + 2)]
            group <- group[inds]
        }
    }
    
    # TODO: generalize code below.
    if (is.data.frame(data)) xx <- apply(x[, -(1:2)], 2, function(y) y / sum(y)) # Fix for YRI datset (first two colums are not counts).
    else xx <- apply(x, 2, function(y) y / sum(y))
    plotDat <- data.frame(ratio = c(xx),
                          iso = if ('transcript_id' %in% names(x)) rep(x$transcript_id, ncol(xx))
                                else rep(1:nrow(xx), ncol(xx)),
                          ind = rep(colnames(xx), each = nrow(xx)),
                          group = rep(group, each = nrow(xx)))
    gg <- ggplot(data = plotDat, aes(x = ind, y = ratio, col = as.factor(iso))) +
         geom_point(size = 2.5, shape = 16, position = position_jitter(height = 0, 
                                                                       width = jitterWidth)) +
        facet_grid(~group, scales = "free_x", space = "free_x") +
        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill = "gray93"),
                              axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        geom_errorbar(aes(ymin = 0, ymax = 1), width = 0, linetype = 2, alpha = 0.1) + 
        scale_x_discrete("") + scale_y_continuous("Splicing ratios") + 
        scale_color_discrete("Isoforms") +  ggtitle(gene)
    gg
}