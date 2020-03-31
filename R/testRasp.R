testRasp <- function(x, data, filterInd, filterExon, transform, ...) {
    
    if(nrow(x)!=length(data$group))
        stop("number of individuals doesn't match")
    if(any(is.na(data$group)))
        stop("'group' variable contains missing values")
    if(filterInd > 0){
        m <- apply(x > 0, 1, mean)
        nfilt <- sum(m < filterInd, na.rm = TRUE)
        if(nfilt > 0){
            ind.keep <- m >= filterInd
            x <- x[ind.keep,]
            warning(nfilt, " individuals filtered out given 'filterInd' argument")
        }
        else{
            ind.keep <- rep(TRUE, nrow(x))
        }
    }
    if(filterExon > 0){
        ## Filter exons with mean relative frequency < 'filterFreq'
        auxa1 <- rowSums(x)
        auxR <- sweep(x, 1, auxa1, FUN="/")
        mrf <- colMeans(auxR)
        nffilt <- sum(mrf < filterExon, na.rm = TRUE)
        if(nffilt > 0){
            x <- x[, mrf >= filter]
            warning(nffilt, " exons filtered out given 'filterExon' argument")
        }
    }
    
    data <- data[ind.keep, , drop=FALSE]
    
    mod <- try(mlm::mlm(x ~ ., data=data, transform=transform, ...), TRUE)
    if (!inherits(mod, "try-error"))
        out <- mod$aov.tab["group", "Pr(>F)"]
    else
        out <- NA
    return(out)
}
