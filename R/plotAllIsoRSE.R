#' Plot proportional exon expression stored in a RangeSummarizedExperiment object
#'
#' @param gene name of the gene to plot.
#' @param data expression dataset. Must be of class `RangedSummarizedExperiment`.
#' @param group a `factor` representing the condition of each individual.
#' @param jitterWidth amount of graphical jitter to use in the plot.
#' @param inds `character` or `numeric` vector. Plot only the specified individuals (samples).
#' @param isos `character` or `numeric` vector. Plot only the specified exons.
#' @examples
#' data(adipose)
#' plotAllIsoRSE("ENSG00000003249.13", adipose.chr16, "inv16p11.2")
#' @import ggplot2
#' @export
plotAllIsoRSE <- function(gene, data, group, jitterWidth = 0, inds, isos) {
    sel <- which(rownames(data) %in% gene)
    x <- SummarizedExperiment::assay(data[sel, ])
    rownames(x) <- paste("exon", 1:nrow(x))
    group <- SummarizedExperiment::colData(data)[, group]
   
     if (!missing(isos)) {
        if (is.character(isos))
            x <- x[x$transcript_id %in% isos, ]
        else
            x <- x[isos, ]
    }
    if (!missing(inds)) {
        if (is.character(inds)) {
            x <- x[, colnames(x) %in% inds]
            group <- group[colnames(x) %in% inds]
        }
        else {
            x <- x[, inds]
            group <- group[inds]
        }
    }
    
    non.all.zero <- !apply(x, 2, function(x) all(x == 0))
    x <- x[, non.all.zero]
    group <- group[non.all.zero]
    xx <- apply(x, 2, function(y) y / sum(y))
    plotDat <- data.frame(ratio = c(xx),
                          iso = rep(rownames(xx), ncol(xx)),
                          ind = rep(colnames(xx), each = nrow(xx)),
                          group = rep(group, each = nrow(xx)))
    gg <- ggplot(data = plotDat, aes(x = ind, y = ratio, col = iso))
    gg <- gg + geom_point(size = 2.5, shape = 16,
                          position = position_jitter(height = 0, width = jitterWidth))
    gg <- gg + facet_grid(~ group, scales = "free_x", space = "free_x")
    gg <- gg + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "gray93"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5))
    gg <- gg + geom_errorbar(aes(ymin = 0, ymax = 1),
                             width = 0, linetype = 2, alpha = 0.1)
    gg <- gg + scale_x_discrete("")
    gg <- gg + scale_y_continuous("Splicing ratios")
    gg <- gg + scale_color_discrete("Isoforms")
    gg <- gg + ggtitle(gene)
    gg
}
