\name{rasp}
\alias{rasp}
\title{
  Alternative splicing test
}
\description{
This function implements two different methods of testing for differential
splicing across two or more conditions. 
}
\usage{
rasp(x, y, expressionCols, geneidCol, test = c("asymptotic", "permutation"), filter = 0.1, transformation = FALSE, Nperms = 1e+05 - 1, type = "median", pairwise = FALSE, mc.cores = 1, nb.mont = 1e+06 - 1, maxTrans = 100)
}
\arguments{
  \item{x}{
    Either an object of class 'ExonCountSet' or a 'data.frame' with both
    gene and transcript ids and their respective expression. For the
    latter the 'expressionCols' and 'geneidCol' columns are mandatory.
}
  \item{y}{
    Object of class 'factor' representing the condition of each
    individual. 
}
  \item{expressionCols}{
    Integer vector representing the numbers of the columns which contain
    the expression. Only necessary when 'x' is of class 'data.frame'.
}
\item{geneidCol}{
  Integer representing the number of the column which contains the
  gene ids. Only necessary when 'x' is of class 'data.frame'.
}
  \item{test}{
    The testing method to be used, "asymptotic" or
    "permutation". Default is 'asymptotic'.
}
  \item{filter}{
    Numeric value between 0 and 1. All transcripts-exons with a
    proportion of non-zero values greater will be filtered out. Default
    is '0.1' which would filter all transcripts with more than a 90\% of zeros.
}
  \item{transformation}{
    Logical indicating whether the values on the boundary should be
    transformed or not. Default is 'FALSE'.
}
  \item{Nperms}{
    Number of permutations to perform when test is equal to
    "permutation". Default is 1e05 - 1.
}
  \item{type}{
    The type to be used in the beta dispersion model. Default is "median".
}
  \item{pairwise}{
    Logical indicating whether pairwise comparisons of group means
    should be performed or not in the permutation test. Default is 'FALSE'.
}
  \item{mc.cores}{
    Number of cores to be used during the computations (NOT AVAILABLE
  FOR WINDOWS USERS).
}
  \item{nb.mont}{
    Number of Montecarlo iterations to be performed when test is equal
    to "asymptotic". Default is 1e06 - 1.
}
  \item{maxTrans}{
    Maximum number of exons- transcripts per gene to be analyzed.
    }
}
\details{
}
\value{
}
\references{

}
\author{

}
\note{

}



\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, expressionCols, geneidCol, test = c("asymptotic", 
    "permutation"), filter = 0.1, transformation = FALSE, Nperms = 1e+05 - 
    1, type = "median", pairwise = FALSE, mc.cores = 1, nb.mont = 1e+06 - 
    1, maxTrans = 100) 
{
    if (class(x) == "ExonCountSet") {
        if (missing(y)) {
            warning("'y' was not specified, the design of 'x' will be used instead")
            y <- design(x)
        }
        else {
            if (!is.factor(y)) 
                stop("'y' should be a factor")
        }
        x2 <- x
        m <- apply(counts(x) > 0, 1, mean)
        filt <- m > filter
        x <- as.data.frame(counts(x)[filt, ])
        xx <- split(x, droplevels(as.factor(geneIDs(x2)[filt])))
    }
    else {
        if (!is.factor(y)) 
            stop("'y' should be a factor")
        m <- apply(x[, expressionCols] > 0, 1, mean)
        x <- x[m > filter, ]
        xx <- split(x[, expressionCols], droplevels(as.factor(x[, 
            geneidCol])))
    }
    test <- match.arg(test)
    if (test == "asymptotic") {
        n <- ncol(x)
        g <- nlevels(y)
    }
    assign("aux", 0, .GlobalEnv)
    N <- length(xx)
    pb <- txtProgressBar(min = 0, max = N, initial = 0, style = 3)
    masterDesc <- try(get("masterDescriptor", envir = getNamespace("parallel")), 
        TRUE)
    f <- function(x, group, test, Nperms, nc, coreID, type) {
        masterDesc <- get("masterDescriptor", envir = getNamespace("parallel"))
        if (masterDesc() == coreID) {
            auxx <- get("aux", envir = .GlobalEnv)
            assign("aux", auxx + 1, .GlobalEnv)
            setTxtProgressBar(pb, nc * aux)
        }
        testRasp(x = x, group = group, test = test, type = type, 
            Nperms = Nperms, transformation = transformation, 
            nb.mont = nb.mont, , multipleTesting = TRUE, randomChisqN = randomChisqN, 
            randomChisqD = randomChisqD, maxTrans = maxTrans)
    }
    f2 <- function(x, group, test, type, Nperms, transformation) {
        aux <<- aux + 1
        setTxtProgressBar(pb, aux)
        testRasp(x = x, group = group, test = test, type = type, 
            Nperms = Nperms, transformation = transformation, 
            nb.mont = nb.mont, multipleTesting = TRUE, randomChisqN = randomChisqN, 
            randomChisqD = randomChisqD, maxTrans = maxTrans)
    }
    if (mc.cores > 1) {
        if (class(masterDesc) == "try-error") 
            stop("It appears you are trying to use multiple cores from Windows, this is not possible")
        mclapp <- try(get("mclapply", envir = getNamespace("parallel")), 
            TRUE)
        detectCor <- try(get("detectCores", envir = getNamespace("parallel")), 
            TRUE)
        nAvailableCores <- detectCor()
        coreID <- mclapp(as.list(1:mc.cores), function(x) masterDesc(), 
            mc.cores = mc.cores)[[1]]
        randomChisqN <- do.call(rbind, mclapp(1:maxTrans, function(i) rchisq(nb.mont, 
            df = g - 1), mc.cores = mc.cores))
        randomChisqD <- do.call(rbind, mclapp(1:maxTrans, function(i) rchisq(nb.mont, 
            df = n - g), mc.cores = mc.cores))
        ans <- mclapp(xx, f, group = y, test = test, Nperms = Nperms, 
            Nc = mc.cores, coreID = coreID, type = type, mc.cores = mc.cores)
    }
    else {
        randomChisqN <- do.call(rbind, lapply(1:maxTrans, function(i) rchisq(nb.mont, 
            df = g - 1)))
        randomChisqD <- do.call(rbind, lapply(1:maxTrans, function(i) rchisq(nb.mont, 
            df = n - g)))
        ans <- lapply(xx, f2, group = y, test = test, type = type, 
            Nperms = Nperms, transformation = transformation)
    }
    setTxtProgressBar(pb, N)
    out <- do.call(rbind, ans)
    nlev <- nlevels(y)
    pvals <- out[, -c(1:(nlev + 1)), drop = FALSE]
    pvals.adj <- apply(pvals, 2, function(x) p.adjust(x, "BH"))
    if (nlev == 2) 
        colnames(pvals.adj) <- "padjust"
    else colnames(pvals.adj) <- "padjust"
    out <- cbind(out, pvals.adj)
    close(pb)
    out
  }
}
