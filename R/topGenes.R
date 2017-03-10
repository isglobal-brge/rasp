topGenes <- function (x, n = 6L, sort.by = "pvalue", pval.adjust.cutoff = 1, print = TRUE, ...) 
{
    if (!inherits(x, "rasp")) 
        stop("x must be an object of class 'rasp'")
    sort.by <- match.arg(sort.by, c("pvalue"))
    nn <- colnames(x) 
    if (print) 
        cat(paste("Comparison of groups:", nn[1], "vs", nn[2], "\n"))
    orderedIdx <- NA
    sortstr <- ""
    orderedIdx <- order(x[, "pvalue"])
    sortstr <- "Corrected p-value"
    x <- x[orderedIdx, ]
    x <- x[x[,"padjust"] <= pval.adjust.cutoff & !is.na(x[,"p.homogeneity"]),]
    if (print) {
        topstr <- "Showing"
        if (nrow(x) > n) 
            topstr <- paste("Showing top", n)
        cat(paste(topstr, "genes ranked by", sortstr, "\n"))
        cat(paste("Maximum adjusted P-value of", pval.adjust.cutoff, 
            "\n"))
        print.data.frame(as.data.frame(x)[1:min(nrow(x), n), 
            ], ...)
    }
    invisible(as.data.frame(x)[1:min(nrow(x), n), ])
}

