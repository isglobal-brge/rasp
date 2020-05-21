#' Plot proportional exon expression
#'
#' @param gene name of the gene to plot.
#' @param data expression dataset. Can be of class `DEXSeqDataSet`, `SummarizedExperiment` or
#'  `RangedSummarizedExperiment`. A `matrix` or `data.frame` is also accepted.
#' @param condition a `factor` representing the condition of each individual.
#'   Required if `data` is a `matrix` or `data.fame`, ignored otherwise.
#' @param fileName base name of the output PDF file.
#' @param factor `character` or `numeric` vector. Plot only the specified individuals (samples).
#' @param title `character` or `numeric` vector. Plot only the specified exons.
#' @param ylab1
#' @param ylab2
#' @param col
#' @param exons
#' @examples
#' data(adipose)
#' genePlot("ENSG00000003249.13", adipose.chr16, "inv16p11.2")
#' 
#' data(YRI)
#' genePlot("ENSG00000215915", YRI, attr(YRI, "gender"))
#' @export
genePlot <- function(gene, data, condition = NULL, fileName = NULL, factor = 0.5, title = FALSE,
                     ylab1 = "gene expression (counts)", ylab2 = "exon relative abundance",
                     col = blues9[3:9], exons = NULL) {
    if (is.data.frame(data))
        data <- data[data$gene_id == gene, ]
    else if (inherits(data, "RangedSummarizedExperiment") |
             inherits(data, "SummarizedExperiment")) {
        data <- SummarizedExperiment::assay(data)[rownames(data) == gene, ]
    }
    # if (is.null(condition)) condition <- pData(data)$condition
    if (!inherits(condition, "factor")) condition <- as.factor(condition)
    if (nlevels(condition) != 2) stop("Only two levels are accepted for 'condition'")

    if (is.null(fileName)) fileName <- paste0(sub(".pdf$", "", gene), ".pdf")
    # ids <- geneIDs(data)
    # geneix <- grep(pattern = gene, ids)
    # if (!length(geneix)) stop("gene not found")
    # exonCounts <- counts(data)[geneix,]
    exonCounts <- counts(data)[gene, ]
    exonNames <- sapply(strsplit(names(ids[geneix]), ":"), '[', 2)

    if(!is.null(exons)){
        if(!any(exons%in%exonNames))
            stop("Could not find selected 'exons'")
        geneix <- geneix[exonNames%in%exons]
        exonNames <- exons
    }
    nExons <- length(geneix)
    meanCond1 <- apply(exonCounts[, condition %in% levels(condition)[1]], 1, mean)
    meanCond2 <- apply(exonCounts[, condition %in% levels(condition)[2]], 1, mean)
    geneSum <- apply(exonCounts, 2, sum)
    exonProps <- apply(exonCounts, 2, function(x) x / sum(x))
    meanDat <- data.frame(mean = c(meanCond1, meanCond2),
                         condition = c(rep(levels(condition)[1], length(meanCond1)),
                             rep(levels(condition)[2], length(meanCond2))))
    propDat <- as.data.frame(t(exonProps))
    names(propDat) <- paste0("Exon",1:nExons)

    require(grDevices)
    unsaturate <- function(col, factor) {
        rgb <- col2rgb(col)
        r <- rgb["red", ]
        g <- rgb["green", ]
        b <- rgb["blue", ]
        conv <- as.list(as.data.frame(t(rgb2hsv(r, g, b))))
        conv[[2]] <- pmin(1, conv[[2]] * factor)
        do.call(hsv, conv)
    }
    dcols <- col
    lcols <- unsaturate(dcols, factor)
    nCond1 <- sum(condition %in% levels(condition)[1])
    nCond2 <- sum(condition %in% levels(condition)[2])
    
    pdf(file = fileName, width = ifelse(nExons <= 2, 10, 15) + nExons/2)
    oldMar <- par()$mar
    par(mar = c(8, 0, 4, 0))
    par(las = 2)
    ws <- ifelse(nExons <= 2, c(2/5, 3/5), c(2/(nExons+1), (nExons-1)/(nExons+1)))
    layout(matrix(c(1,2), nrow = 1), widths = ws)
    bp <- boxplot(mean ~ condition, data = meanDat, col = c("lightgrey", "white"),
                  frame = FALSE, axes = FALSE, par(pch = 20))
    axis(4)
    text(par("usr")[2],(par("usr")[3]+par("usr")[4])/2 , labels = ylab1,  srt = 90, xpd = NA, adj = c(0.5, 7), cex = 1.2)
    lines(1:2, rep(par("usr")[3] - (par("usr")[4]-par("usr")[3])/100, 2),  xpd = TRUE)
    text(1.5,  par("usr")[3], labels = gene, srt = 90, xpd = TRUE, adj = c(1.5,.5))
    par(mar = c(8, 5, 4, 6))
    for (i in 1:nExons) {
        if (i == 1) {
            boxplot(propDat[,i]~condition, col = c(lcols[i], dcols[i]),
                    xlim = c(0, 3*nExons),
                    frame = FALSE, at = c(1+(i-1)*3, 2+(i-1)*3),
                    axes = FALSE, ylim = c(0,1), par(pch = 20))
            lines(1:2, c(-.05,-.05), xpd = TRUE)
            text(1.5,  par("usr")[3], labels = exonNames[1], srt = 90, xpd = TRUE, adj = c(1.5,.5))
        } else {
            boxplot(propDat[,i]~condition, col = c(lcols[i], dcols[i]),
                    add = TRUE, 
                    frame = FALSE, at = c(1+(i-1)*3, 2+(i-1)*3),
                    axes = FALSE, par(pch = 20))
            lines(c(1+(i-1)*3, 2+(i-1)*3), c(-.05,-.05), xpd = TRUE)
            text((3+(i-1)*6)/2, par("usr")[3], labels = exonNames[i], srt = 90, xpd = TRUE, adj = c(1.5,.5))
        }
    }
    
    if (title) mtext(paste0(gene, ": ", levels(condition)[1], " vs. ", levels(condition)[2]), side = 3, padj = 1, col = "darkblue", font = 2, las = 1, cex = 1.5, line = 3, at = 2)
    text(par("usr")[2],(par("usr")[3]+par("usr")[4])/2 , labels = ylab2,  srt = 90, xpd = NA, adj = c(0.5, 7), cex = 1.2)
    axis(4)
    dev.off()
    on.exit(par(mar = oldMar, las = 0))
}
