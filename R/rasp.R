#' rasp
#' 
#' Differential alternative splicing test.
#' 
#' Test for differential alternative splicing across two or more conditions.
#' 
#' @param formula a `formula` with a symbolic description of the model to be fitted.
#' @param x expression dataset. Can be of class `DEXSeqDataSet`, `SummarizedExperiment` or
#'  `RangedSummarizedExperiment`. A `matrix` or `data.frame` is also accepted.
#' @param group a `factor` representing the condition of each individual.
#'   Required if `x` is a `matrix` or `data.fame`, ignored otherwise.
#' @param expressionCols `character` or `numeric` vector specifying the columns of `x` with that contain
#'   the expression data (i.e. the different samples).
#'   Required if `x` is a `matrix` or `data.fame`, ignored otherwise.
#' @param geneidCol name or index of the column of `x` encoding the gene IDs.
#'   Required if `x` is a `matrix` or `data.fame`, ignored otherwise.
#' @param filterInd filter individuals (samples) with mean relative frequency < `filterInd`.
#' @param filterExon filter exons with mean relative frequency < `filterFreq`.
#' @param transform argument for `mlm`. Can be any of `none`, `sqrt` or `log`.
#' @param cores `integer`. Number of parallel processes to use.
#' @param ... additional arguments passed to `mlm`.
#' @return A `rasp` class `matrix` containing the p-values (also adjusted for multiple testing).
#'   `NA` will be returned for single exon genes and genes not passing filtering criteria.
#' @examples
#' data(adipose)
#' rasp(~ inv16p11.2, adipose.chr16)
#' 
#' data("YRI")
#' rasp(x = YRI, group = attr(YRI, "gender"), expressionCols = 3:71, geneidCol = "gene_id")
#' @export
rasp <- function(formula, x, group, expressionCols, geneidCol,
                 filterInd = 0.1,
                 filterExon = 0.05, 
                 transform = "none",
                 cores = 1 , ...) {
    
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
      
      if (any(!vars %in% colnames(SummarizedExperiment::colData(x))))
        stop("variables in the formula should be in the metadata (e.g. 'colData()')")
      
      data <- SummarizedExperiment::colData(x)[, vars, drop=FALSE]
      colnames(data)[1] <- "group"
      
      if (!is.factor(data$group)) stop(paste(vars[1], "should be a factor"))

      x <- lapply(split(x, rownames(x)), SummarizedExperiment::assay)

    } else {
      if (missing(group))
        stop("Providing count data as a data.frame requires a grouping variable to be passed through 'group' argument")
      data <- data.frame(group=group)
      if (!is.factor(data$group)) stop("'group' should be a factor")
      if (missing(expressionCols) | missing(geneidCol))
        stop("arguments 'geneidCol' and 'expressionCols' must be provided")
      x <- split(x[, expressionCols], droplevels(as.factor(x[, geneidCol])))
    }

    if (length(x) > 1){
      cl <- parallel::makeCluster(cores)
      
      
      doSNOW::registerDoSNOW(cl)
      pb <- txtProgressBar(max = length(x), style = 3)
      opts <- list(progress = function(i) setTxtProgressBar(pb, i))
    
       # ans <- foreach::`%dopar%`(foreach::foreach (i = 1:length(x),
       #                                           .export = "testRasp",
       #                                           .options.snow = opts), {
       # nm <- names(x)[i]
       # if (all(x[[nm]] == 0)) NA # mlm crashes on empty (0 count) genes.
       # else if (nrow(x[[nm]]) == 1) NA # mlm crashes on single exon genes.
       # else testRasp(t(x[[nm]]), data = data, 
       #               filterInd = filterInd,
       #               filterExon = filterExon, 
       #               transform = transform, ...)
       # })

      ans <- foreach (i = 1:length(x),.export = "testRasp",.options.snow = opts) %dopar% {
        nm <- names(x)[i]
        if (all(x[[nm]] == 0)) NA # mlm crashes on empty (0 count) genes.
        else if (nrow(x[[nm]]) == 1) NA # mlm crashes on single exon genes.
        else testRasp(t(x[[nm]]), data = data, 
                      filterInd = filterInd,
                      filterExon = filterExon, 
                      transform = transform, ...)
        }
      
     close(pb)
     parallel::stopCluster(cl)
    } else{
      ans <- testRasp(t(x[[1]]), data=data, 
               filterInd = filterInd,
               filterExon = filterExon, 
               transform = transform, ...)
    }

    # Post-process results.
    pvals <- as.matrix(unlist(ans))
    if (nrow(pvals) > 1){
      nexons <- sapply(x, nrow)
      pvals.adj <- apply(pvals, 2, p.adjust, "BH")
      out <- cbind(nexons, pvals, pvals.adj)
      rownames(out) <- names(x)
      colnames(out) <- c("n.exons", "pvalue", "padjust")
    }
    else{
      nexons <- nrow(x[[1]])
      out <- cbind(nexons, pvals)
      rownames(out) <- names(x)
      colnames(out) <- c("n.exons", "pvalue")
    }
      
    class(out) <- c("rasp", "matrix")
    out
}
