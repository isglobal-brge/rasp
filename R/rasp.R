rasp <- function(x, y, expressionCols, geneidCol, test = c("asymptotic", "permutation"), filter = 0.1, filterFreq = 0.05, transformation = FALSE, Nperms = 1e5-1, type = "median", pairwise = FALSE, mc.cores = 1, maxTrans = 100, testGroup = FALSE, verbose = FALSE) {
    filt <- TRUE # TODO: what does this undefined variable do?
    
    # TODO: Remove support for ExonCountSet, since is deprecated.
    if (class(x) == "ExonCountSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor") 

        x2 <- x
        xx <- split(x, droplevels(as.factor(geneIDs(x2)[filt])))
    
    } else if (class(x) == "DEXSeqDataSet") {
        if(missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        } else if (!is.factor(y)) stop("'y' should be a factor")
        
        x2 <- x
        xx <- split(x, droplevels(as.factor(geneIDs(x2)[filt])))
        
    } else if (class(x) == "RangedSummarizedExperiment") {
      # if(missing(y)) {
      #   warning("'y' was not specified, the design of 'x' will be used instead")
      #   y <- design(x)
      # } else if (!is.factor(y)) stop("'y' should be a factor")
      if (!is.factor(y)) stop("'y' should be a factor")
      
      x2 <- x
      xx <- lapply(split(x, names(rowRanges(x))[filt]), assay)
      
    } else {
        if (!is.factor(y)) stop("'y' should be a factor") 
        xx <- split(x[, expressionCols], droplevels(as.factor(x[, geneidCol])))
    }
    
    test <- match.arg(test)
    if (test == "asymptotic") {
        # TODO: are these ever used?
        n <- ncol(x)
        g <- nlevels(y)
    }

    # assign("aux", 0, .GlobalEnv)
    pb <- txtProgressBar(max = length(xx), style = 3)
    masterDesc <- try(get("masterDescriptor", envir = getNamespace("parallel")), TRUE)
    
    # TODO: take these functions outside.
    f <- function(x, group, test, Nperms, nc, coreID, type){
        masterDesc <- get("masterDescriptor", envir = getNamespace("parallel"))
        if (masterDesc() == coreID) {
            # auxx <- get("aux", envir = .GlobalEnv)
            assign("aux", auxx + 1, .GlobalEnv)
            setTxtProgressBar(pb, nc * aux)
        }
        testRasp(x = x, group = group, test = test, type = type, Nperms = Nperms,
                 transformation=transformation, multipleTesting=TRUE, maxTrans=maxTrans,
                 pairwise = pairwise, testGroup = testGroup, filter = filter)
    }
    f2 <- function(nm, group, test, type, Nperms, transformation, verb){
        if (verb) print(nm)
        # aux <<- aux + 1
        # setTxtProgressBar(pb, aux)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
        testRasp(x = xx[[nm]], group = group, test = test, type = type, Nperms = Nperms,
                 transformation=transformation, multipleTesting=TRUE, maxTrans=maxTrans,
                 pairwise = pairwise, testGroup = testGroup, filter = filter)
    }
    
    if (mc.cores > 1) {
        if (class(masterDesc) == "try-error")
            stop("It appears you are trying to use multiple cores from Windows, this is not possible")
        mclapp <- try(get("mclapply", envir = getNamespace("parallel")), TRUE)
        detectCor <- try(get("detectCores", envir = getNamespace("parallel")),TRUE)
        nAvailableCores <- detectCor()
        coreID <- mclapp(as.list(1:mc.cores), function(x) masterDesc(),
                         mc.cores = mc.cores)[[1]]
        ans <- mclapp(xx, f, group=y, test=test, Nperms=Nperms,
                      nc = mc.cores, coreID = coreID, type = type, mc.cores = mc.cores)
    } else {
        ans <- lapply(names(xx), f2, group=y, test=test, type = type, Nperms=Nperms,
                      transformation=transformation, verb = verbose)
    }
    setTxtProgressBar(pb, length(xx))
    
    out <- do.call(rbind, ans) 
    nlev <- nlevels(y)
    pvals <- out[, -c(1:(nlev + 1)), drop = FALSE]
    if(ncol(pvals) > 1) pvals[, 2:ncol(pvals)] <- t(apply(pvals[,2:ncol(pvals)], 1, p.adjust, method = "holm"))
    pvals.adj <- apply(pvals, 2, p.adjust, "BH")
    colnames(pvals.adj) <- gsub("pvalue", "padjust", colnames(pvals.adj))
    out <- as.data.frame(cbind(out, pvals.adj))
    class(out) <- c("data.frame", "rasp")
    
    close(pb)
    
    out
}
