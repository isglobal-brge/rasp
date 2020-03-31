# expressionCols: character or integer vector specifying the columns with expr data (samples).
# geneidCol: name or number of the column encoding the gene IDs.

rasp <- function(x, y, expressionCols, geneidCol, cores = parallel::detectCores()) {
    if (class(formula) == "formula") {
      next
    }
    
    # TODO: Maybe remove support for ExonCountSet, since is deprecated.
    if (class(x) == "ExonCountSet") {
        warning("'ExonCountSet' is deprecated, please use a 'DEXSeqDataSet' object.")
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor") 

        x <- split(x, droplevels(as.factor(geneIDs(x))))
    
    } else if (class(x) == "DEXSeqDataSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor")
        
        x <- split(x, droplevels(as.factor(geneIDs(x))))
        
    } else if (class(x) == "SummarizedExperiment" ||
               class(x) == "RangedSummarizedExperiment") {
      # if(missing(y)) {
      #   warning("'y' was not specified, the design of 'x' will be used instead")
      #   y <- design(x)
      # } else if (!is.factor(y)) stop("'y' should be a factor")
      if (!is.factor(y)) stop("'y' should be a factor")
      
      x <- lapply(split(x, names(SummarizedExperiment::rowRanges(x))), SummarizedExperiment::assay)
      
    } else {
        if (!is.factor(y)) stop("'y' should be a factor")
        x <- split(x[, expressionCols], droplevels(as.factor(x[, geneidCol])))
    }
    
    # Compute over all genes.
    # pb <- txtProgressBar(max = length(x), style = 3)
    
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    ans <- foreach::`%dopar%`(foreach::foreach (nm = names(x)), {
      # setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      if (all(x[[nm]] == 0)) NA # mlm crashes on empty (0 count) and single-exon genes.
      else if (nrow(x[[nm]]) == 1) NA # mlm crashes on single exon genes.
      else mlm::mlm(t(x[[nm]]) ~ y)$aov.tab["y", "Pr(>F)"]
    })
    
    parallel::stopCluster(cl)

    # setTxtProgressBar(pb, length(x))
    
    # Post-process results.
    # out <- do.call(rbind, ans) # TODO: substitute by unlist.
    pvals <- as.matrix(unlist(ans))
    # nlev <- nlevels(y)
    # pvals <- out[, -c(1:(nlev + 1)), drop = FALSE]
    # pvals <- out[, 1, drop = FALSE] # TODO: remove this, just call pvals the previous 'out' variable.
    # if(ncol(pvals) > 1) pvals[, 2:ncol(pvals)] <- t(apply(pvals[,2:ncol(pvals)], 1, p.adjust, method = "holm"))
    pvals.adj <- apply(pvals, 2, p.adjust, "BH")
    # colnames(pvals.adj) <- gsub("pvalue", "padjust", colnames(pvals.adj))
    # out <- as.data.frame(cbind(out, pvals.adj))
    out <- as.data.frame(cbind(pvals, pvals.adj))
    names(out) <- c("pvalue", "padjust")
    rownames(out) <- names(x)
    class(out) <- c("data.frame", "rasp")
    
    # close(pb)
    
    out
}