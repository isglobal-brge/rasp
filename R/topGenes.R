#' Get the top ranking genes.
#' 
#' By default, `topGenes` will return *and* print the 6 genes with best p-value.
#'
#' @param x object returned by `rasp`.
#' @param n number of top genes to return.
#' @param sort.by column to be used as a sorting key.
#' @param pval.adjust.cutoff return only genes with a `padjust` value lower than this.
#' @param print whether to print the output, in addition to only returning.
#' @param ... additional arguments passed to `print.data.frame`.
#' @return a `data.frame`
#'
#' @export
topGenes <- function (x, n = 6L, sort.by = "pvalue", pval.adjust.cutoff = 1, print = TRUE, ...) {
    if (!inherits(x, "rasp")) stop("x must be an object of class 'rasp'")
    if (pval.adjust.cutoff < 0 || pval.adjust.cutoff > 1) stop("'pval.adjust.cutoff' must be between 0 and 1")
    # sort.by <- match.arg(sort.by, c("pvalue"))
    orderedIdx <- NA
    sortstr <- ""
    orderedIdx <- order(x[, sort.by])
    sortstr <- "Corrected p-value"
    x <- x[orderedIdx, ]
    x <- x[x[, "padjust"] <= pval.adjust.cutoff ,]
    if (print) {
        topstr <- "Showing"
        if (nrow(x) > n) 
            topstr <- paste("Showing top", n)
        cat(paste(topstr, "genes ranked by", sortstr, "\n"))
        cat(paste("Maximum adjusted P-value of", pval.adjust.cutoff, "\n"))
        print.data.frame(as.data.frame(x)[1:min(nrow(x), n), ], ...)
    }
    invisible(as.data.frame(x)[1:min(nrow(x), n), ])
}
