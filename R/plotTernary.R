plotTernary <- function(x, gene, group, geneCol=1, transCol=2, transcripts, ...) {
 
  if (missing(transcripts))
    dims <- c(1,2,3)
  else
    dims <- transcripts

  if (inherits(x, "ExonCountSet")) {
   x <- counts(x)
   xx <- x[grep(gene, rownames(x)),]
   rownames(xx) <- unlist(lapply(strsplit(rownames(xx), "\\:"), function(x) x[2]))
  }
  else {
    xx <- x[x[,1]==gene, ]
    rownames(xx) <- xx[,transCol] 
    xx <- xx[, -c(geneCol, transCol)]
  }
  total <- colSums(xx)
  ratio <- sweep(xx, 2, FUN="/", total)
  ratio.dr <- DR_data(t(ratio))
  if (missing(group))
    plot(ratio.dr, dims=dims, ...)
  else {
    mycol <- c("red", "blue", "darkgreen", "yellow", "orange")[1:nlevels(group)]
    ok <- as.character(factor(group, labels=mycol))
    plot(ratio.dr, col=ok, dims=dims, ...)
    legend("topleft", levels(group), col=mycol, pch=16)
  }       
}
