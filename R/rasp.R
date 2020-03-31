# expressionCols: character or integer vector specifying the columns with expr data (samples).
# geneidCol: name or number of the column encoding the gene IDs.

rasp <- function(x, y, expressionCols, geneidCol, cores = parallel::detectCores()) {
    # Prepare data list.
    if (class(x) == "DEXSeqDataSet") {
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
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    ans <- foreach::`%dopar%`(foreach::foreach (nm = names(x)), {
      if (all(x[[nm]] == 0)) NA # mlm crashes on empty (0 count) and single-exon genes.
      else if (nrow(x[[nm]]) == 1) NA # mlm crashes on single exon genes.
      else mlm::mlm(t(x[[nm]]) ~ y)$aov.tab["y", "Pr(>F)"]
    })

    parallel::stopCluster(cl)

    # Post-process results.
    pvals <- as.matrix(unlist(ans))
    pvals.adj <- apply(pvals, 2, p.adjust, "BH")
    out <- cbind(pvals, pvals.adj)
    rownames(out) <- names(x)
    colnames(out) <- c("pvalue", "padjust")
    class(out) <- c("rasp", "matrix")

    out
}