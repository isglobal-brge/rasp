testRasp <- function(x, data, filterInd, filterExon, transform, ...) {

    # Check arguments.
    if (filterInd < 0 | filterInd > 1) stop("'filterInd' must be between 0 and 1")
    if (filterExon < 0 | filterExon > 1) stop("'filterExon' must be between 0 and 1")
    
    if (nrow(x) != length(data$group))
      stop("number of individuals doesn't match")
    if (any(is.na(data$group)))
      warning("WARNING : 'group' variable contains missing values")
    
    # Filter individuals (samples).
    if (filterInd > 0){
      m <- apply(x > 0, 1, mean)
      nfilt <- sum(m < filterInd, na.rm = TRUE)
      if (nfilt > 0){
        if (nfilt == nrow(x)) return(NA)
        else {
          ind.keep <- m >= filterInd
          x <- x[ind.keep, , drop=FALSE]
          warning(nfilt, " individuals filtered out given 'filterInd' argument")
          if (nrow(x) == 1) return(NA)
        }
      }
      else ind.keep <- TRUE
    }
    
    # Filter exons.
    if (filterExon) {
      # Filter exons with mean relative frequency < 'filterFreq'
      auxa1 <- rowSums(x)
      auxR <- sweep(x, 1, auxa1, FUN="/")
      mrf <- colMeans(auxR)
      nffilt <- sum(mrf < filterExon, na.rm = TRUE)
      if(nffilt > 0){
        if(nffilt == ncol(x)) return(NA)
        else{
          x <- x[, mrf >= filterExon, drop=FALSE]
          warning(nffilt, " exons filtered out given 'filterExon' argument")
          if (ncol(x) == 1) return(NA)
        }
      }
    }
    
    
    data <- data[ind.keep, , drop=FALSE]
    
    ss <- rowSums(x)
    x <- sweep(x, 1, FUN="/", ss) 
    
    # Compute p-values.
    mod <- try(mlm::mlm(x ~ ., data=data, transform=transform, ...), TRUE)
    if (!inherits(mod, "try-error")) 
      ans <- mod$aov.tab["group", "Pr(>F)"]
    else return(NA) 
    
    
  
  return(ans)

}
