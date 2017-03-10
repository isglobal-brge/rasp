genePlot <- function(data, gene, condition = NULL, fileName = NULL, factor = 0.5, title = FALSE,
                     ylab1 = "gene expression (counts)", ylab2 = "exon relative abundance", col = blues9[3:9], exons = NULL){
    if(!inherits(data, "ExonCountSet"))
        stop("'data' must be of class 'ExonCountSet'")
    if(is.null(condition))
        condition <- pData(data)$condition
    if(nlevels(condition) != 2)
        stop("At the moment, only two levels are accepted for 'condition'")
    if(!inherits(condition, "factor"))
        condition <- as.factor(condition)
    if(is.null(fileName))
        fileName <- paste0(gene,".pdf")
    conditionLevs <- levels(condition)
    ids <- geneIDs(data)
    geneix <- grep(pattern = gene, ids)
    if(length(geneix) == 0)
        stop("gene not found")
    exonCounts <- counts(data)[geneix,]
    exonNames <- unlist(lapply(strsplit(names(ids[geneix]), ":"),
                               function(x)
                               x[2]
                               ))
    if(!is.null(exons)){
        if(!any(exons%in%exonNames))
            stop("Could not find selected 'exons'")
        geneix <- geneix[exonNames%in%exons]
        exonNames <- exons
    }
    nExons <- length(geneix)
    meanCond1 <- apply(exonCounts[,condition%in%conditionLevs[1]], 1, mean)
    meanCond2 <- apply(exonCounts[,condition%in%conditionLevs[2]], 1, mean)
    geneSum <- apply(exonCounts, 2, sum)
    exonProps <- apply(exonCounts, 2, function(x)
                       x/sum(x))
    meanDat <- data.frame(mean = c(meanCond1, meanCond2),
                         condition = c(rep(conditionLevs[1], length(meanCond1)),
                             rep(conditionLevs[2], length(meanCond2))))
    propDat <- as.data.frame(t(exonProps))
    names(propDat) <- paste0("Exon",1:nExons)
                                        #    cols <- c("blue", "green", "grey", "cyan", "salmon",
                                        #              "goldenrod", "seagreen")
                                        #    lcols <- paste0("light",cols)
                                        #    dcols <- paste0("dark",cols)
    require(grDevices)
    unsaturate <- function(col, factor) {
        rgb <- col2rgb(col)
        r <- rgb["red",]
        g <- rgb["green",]
        b <- rgb["blue",]
        conv <- as.list(as.data.frame(t(rgb2hsv(r, g, b))))
        conv[[2]] <- pmin(1, conv[[2]] * factor)
        do.call(hsv, conv)
    }
    dcols <- col
    lcols <- unsaturate(dcols, factor)
    nCond1 <- sum(condition%in%conditionLevs[1])
    nCond2 <- sum(condition%in%conditionLevs[2])
    pdf(file = fileName, width = ifelse(nExons <= 2, 10, 15) + nExons/2)
    oldMar <- par()$mar
    par(mar = c(8, 0, 4, 0))
    par(las = 2)
    ws <- ifelse(nExons <= 2, c(2/5, 3/5), c(2/(nExons+1), (nExons-1)/(nExons+1)))
#    ws <- c(/nExons, (nExons-1)/nExons)
    layout(matrix(c(1,2), nrow = 1), widths = ws)
    bp <- boxplot(mean ~ condition, data = meanDat, col = c("lightgrey", "white"),
                  frame = FALSE, axes = FALSE, par(pch = 20))
#    axis(1, tick = FALSE, labels = conditionLevs, at = c(1,2))
    axis(4)
    text(par("usr")[2],(par("usr")[3]+par("usr")[4])/2 , labels = ylab1,  srt = 90, xpd = NA, adj = c(0.5, 7), cex = 1.2)
    lines(1:2, rep(par("usr")[3] - (par("usr")[4]-par("usr")[3])/100, 2),  xpd = TRUE)
    text(1.5,  par("usr")[3], labels = gene, srt = 90, xpd = TRUE, adj = c(1.5,.5))
#    text(1:2, par("usr")[3], labels = conditionLevs, srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex = .9)
    par(mar = c(8, 5, 4, 6))
    for(i in 1:nExons){
        if(i == 1){
            boxplot(propDat[,i]~condition, col = c(lcols[i], dcols[i]),
                    xlim = c(0, 3*nExons),
                    frame = FALSE, at = c(1+(i-1)*3, 2+(i-1)*3),
                    axes = FALSE, ylim = c(0,1), par(pch = 20))
            lines(1:2, c(-.05,-.05), xpd = TRUE)
#            text(1:2, par("usr")[3], labels = paste(exonNames[1], conditionLevs, sep = "|"), srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex = .9)
            text(1.5,  par("usr")[3], labels = exonNames[1], srt = 90, xpd = TRUE, adj = c(1.5,.5))
        }
        else{
            boxplot(propDat[,i]~condition, col = c(lcols[i], dcols[i]),
                    add = TRUE, 
                    frame = FALSE, at = c(1+(i-1)*3, 2+(i-1)*3),
                    axes = FALSE, par(pch = 20))
            lines(c(1+(i-1)*3, 2+(i-1)*3), c(-.05,-.05), xpd = TRUE)
#            lines(c(i+1.5*(i-1), i+1.5*(i-1)+1), c(-.02,-.02), xpd = TRUE)
            text((3+(i-1)*6)/2, par("usr")[3], labels = exonNames[i], srt = 90, xpd = TRUE, adj = c(1.5,.5))
#            text(i+1.5*(i-1)+1/2, par("usr")[3], labels = exonNames[i], srt = 90, xpd = TRUE, adj = c(1.5,.5))
#            text(c(1+(i-1)*3, 2+(i-1)*3), par("usr")[3], labels = paste(exonNames[i], conditionLevs, sep = "|"), srt = 90, adj = c(1,.5), xpd = TRUE, cex = .9)
        }
    }
    if(title)
        mtext(paste0(gene, ": ", conditionLevs[1], " vs. ", conditionLevs[2]), side = 3, padj = 1, col = "darkblue", font = 2, las = 1, cex = 1.5, line = 3, at = 2)
    text(par("usr")[2],(par("usr")[3]+par("usr")[4])/2 , labels = ylab2,  srt = 90, xpd = NA, adj = c(0.5, 7), cex = 1.2)
    axis(4)
    dev.off()
    on.exit(par(mar = oldMar, las = 0))
}
