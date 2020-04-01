plotAllIso <- function(gene, data, group, jitterWidth = 0,
                       inds, isos){
    require(ggplot2)
    if (is.data.frame(data))
        x <- data[data$gene_id == gene,]
    else if (inherits(data, "RangedSummarizedExperiment") |
             inherits(data, "SummarizedExperiment")){
        x <- lapply(split(data, gene), SummarizedExperiment::assay)[[1]]
        group <- colData(data)[,group]
    }
    if(!missing(isos)){
        if(is.character(isos))
            x <- x[x$transcript_id%in%isos,]
        else
            x <- x[isos,]
    }
    if(!missing(inds)){
        if(is.character(inds)){
            x <- x[,colnames(x)%in%inds]
            group <- group[colnames(x)%in%inds]
        }
        else{
            x <- x[,c(1,2,inds+2)]
            group <- group[inds]
        }
    }
    xx <- apply(x[,3:ncol(x)], 2, function(y)
        y/sum(y))
    plotDat <- data.frame(ratio = c(xx),
                          iso = rep(x$transcript_id, ncol(xx)),
                          ind = rep(colnames(xx), each = nrow(xx)),
                          group = rep(group, each = nrow(xx)))
    gg <- ggplot(data = plotDat, aes(x = ind, y = ratio, col = iso))
    gg <- gg + geom_point(size = 2.5, shape = 16,
                          position = position_jitter(height = 0,
                              width = jitterWidth))
    gg <- gg + facet_grid(~group, scales = "free_x", space = "free_x")
    gg <- gg + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(
                         fill = "gray93"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5))
    gg <- gg + geom_errorbar(aes(ymin = 0, ymax = 1),
                             width = 0, linetype = 2,
                             alpha = 0.1)
    gg <- gg + scale_x_discrete("")
    gg <- gg + scale_y_continuous("Splicing ratios")
    gg <- gg + scale_color_discrete("Isoforms")
    gg <- gg + ggtitle(gene)
    gg
}
