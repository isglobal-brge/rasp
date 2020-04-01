# expressionCols: character or integer vector specifying the columns with expr data (samples).
# geneidCol: name or number of the column encoding the gene IDs.

rasp <- function(formula, x, group, expressionCols, geneidCol, filterInd = 0.1 ,
                 filterExon = 0.05, 
                 transform = "none",
                 cores = parallel::detectCores() -1, ...) {
    
    
    # Prepare data list.
    if (class(x) == "DEXSeqDataSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor")

        x <- split(x, droplevels(as.factor(geneIDs(x))))

    } else if (class(x) == "SummarizedExperiment" ||
               class(x) == "RangedSummarizedExperiment") {
      
      vars <- all.vars(formula)
      
      if(any(!vars%in%colnames(colData(x))))
        stop("variables in the formula should be in the metadata (e.g. 'colData()')")
      
      data <- colData(x)[, vars, drop=FALSE]
      colnames(data)[1] <- "group"
      
      if (!is.factor(data$group)) stop(paste(vars[1], "should be a factor"))

      x <- lapply(split(x, BiocGenerics::rownames(x)), 
                  SummarizedExperiment::assay)

    } else {
      if (missing(group))
        stop("Providing count data as a data.frame requires a grouping variable to be passed through 'group' argument")
      data <- data.frame(group=group)
      if (!is.factor(data$group)) stop("'y' should be a factor")
      x <- split(x[, expressionCols], droplevels(as.factor(x[, geneidCol])))
    }

    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    ans <- foreach::`%dopar%`(foreach::foreach (nm = names(x)), {
      if (all(x[[nm]] == 0)) NA # mlm crashes on empty (0 count) and single-exon genes.
      else if (nrow(x[[nm]]) == 1) NA # mlm crashes on single exon genes.
      else testRasp(t(x[[nm]]), data=data, 
                    filterInd = filterInd , filterExon = filterExon, 
                    transform = transform, ...)
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