testRasp <- function(x, data, filterInd, filterExon, transform, ...) {
    # Check arguments.
    if (filterInd < 0 | filterInd > 1) stop("'filterInd' must be between 0 and 1")
    if (filterExon < 0 | filterExon > 1) stop("'filterExon' must be between 0 and 1")
  
    if (nrow(x) != length(data$group))
        stop("number of individuals doesn't match")
    if (any(is.na(data$group)))
        stop("'group' variable contains missing values")
  
    # Filter individuals (samples).
    if (filterInd) {
        m <- rowMeans(x > 0)
        nfilt <- sum(m < filterInd, na.rm = TRUE)
        if (nfilt) {
            ind.keep <- m >= filterInd
            x <- x[ind.keep, , drop = FALSE]
            warning(nfilt, " individuals filtered out given 'filterInd' argument")
            # TODO: if all individuals excluded, return NA.
        }
        else ind.keep <- TRUE
    }
  
    # Filter exons.
    if (filterExon) {
        ## Filter exons with mean relative frequency < 'filterFreq'
        auxa1 <- rowSums(x)
        auxR <- sweep(x, 1, auxa1, FUN="/")
        mrf <- colMeans(auxR)
        nffilt <- sum(mrf < filterExon, na.rm = TRUE)
        if(nffilt > 0){
            x <- x[, mrf >= filterExon]
            warning(nffilt, " exons filtered out given 'filterExon' argument")
            # TODO: if all exons excluded, return NA.
        }
    }
    
    data <- data[ind.keep, , drop=FALSE]
    
    # Compute p-values.
    mod <- try(mlm::mlm(x ~ ., data=data, transform=transform, ...), TRUE)
    if (!inherits(mod, "try-error")) mod$aov.tab["group", "Pr(>F)"]
    else NA
}