rasp <- function(x, y, expressionCols, geneidCol, test = c("asymptotic", "permutation"), filter = 0.1, filterFreq = 0.05, transformation = FALSE, Nperms = 1e5-1, type = "median", pairwise = FALSE, cores = parallel::detectCores(), maxTrans = 100, testGroup = FALSE, verbose = FALSE) {
    # TODO: add support for formula notation.
    if (0) next
    
    # TODO: Remove support for ExonCountSet, since is deprecated.
    if (class(x) == "ExonCountSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor") 

        xx <- split(x, droplevels(as.factor(geneIDs(x))))
    
    } else if (class(x) == "DEXSeqDataSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor")
        
        xx <- split(x, droplevels(as.factor(geneIDs(x))))
        
    } else if (class(x) == "RangedSummarizedExperiment") {
      # if(missing(y)) {
      #   warning("'y' was not specified, the design of 'x' will be used instead")
      #   y <- design(x)
      # } else if (!is.factor(y)) stop("'y' should be a factor")
      if (!is.factor(y)) stop("'y' should be a factor")
      
      xx <- lapply(split(x, names(SummarizedExperiment::rowRanges(x))), SummarizedExperiment::assay)
      
    } else {
        if (!is.factor(y)) stop("'y' should be a factor") 
        xx <- split(x[, expressionCols], droplevels(as.factor(x[, geneidCol])))
    }
    
    # Compute over all genes.
    # pb <- txtProgressBar(max = length(xx), style = 3)
    
    # ans <- list()
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    ans <- foreach::`%dopar%`(foreach::foreach (nm = 1:length(xx)), {
      # setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      if (all(xx[[nm]] == 0)) NA # mlm crashes on empty (0 count) and single-exon genes.
      else if (nrow(xx[[nm]]) == 1) NA # mlm crashes on single exon genes.
      else mlm::mlm(t(xx[[nm]]) ~ y)$aov.tab["y", "Pr(>F)"]
    })
    
    parallel::stopCluster(cl)

    # setTxtProgressBar(pb, length(xx))
    
    # Post-process results.
    out <- do.call(rbind, ans) # TODO: substitute by unlist.
    # nlev <- nlevels(y)
    # pvals <- out[, -c(1:(nlev + 1)), drop = FALSE]
    pvals <- out[, 1, drop = FALSE] # TODO: remove this, just call pvals the previous 'out' variable.
    if(ncol(pvals) > 1) pvals[, 2:ncol(pvals)] <- t(apply(pvals[,2:ncol(pvals)], 1, p.adjust, method = "holm"))
    pvals.adj <- apply(pvals, 2, p.adjust, "BH")
    # colnames(pvals.adj) <- gsub("pvalue", "padjust", colnames(pvals.adj))
    out <- as.data.frame(cbind(out, pvals.adj))
    names(out) <- c("pvalue", "padjust")
    rownames(out) <- names(xx)
    class(out) <- c("data.frame", "rasp")
    
    close(pb)
    
    out
}