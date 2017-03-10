SS <- function(d, group, lev){
    dd <- as.matrix(d)
    ddg <- dd[group == lev, group == lev]
    ddg2 <- ddg^2
    sum(ddg2)/(2*nrow(ddg))
}

fSt <- function(ss, group, lev){
    ss <- unlist(ss)
    (nlevels(group)-1)*ss[lev]/sum(ss[levels(group)!=lev])
}

permTestSGi <- function(d, group, lev, NpermsSG){
    levs <- levels(group)
    unlist(lapply(1:NpermsSG, function(i){
        gp <- sample(group)
        ssp <- lapply(levs, function(lev)
            SS(d, gp, lev))
        names(ssp) <- levs
        fSt(ssp, gp, lev)
    }))
}

permTestSG <- function(d, group, NpermsSG){
    levs <- levels(group)
    ss <- lapply(levs, function(lev)
        SS(d, group, lev))
    names(ss) <- levs
    unlist(lapply(levs, function(lev){
        F0 <- fSt(ss, group, lev)
        Fp <- permTestSGi(d, group, lev, NpermsSG)
        pval <- sum(Fp < F0)/NpermsSG
        return(pval)
    }))
}
 
testRasp <- function(x, group, test="asymptotic", type="median", transformation=FALSE,
                     Nperms=1e5-1, 
                     maxTrans=100, pairwise = FALSE,
                     NpermsSG = 1e3,
                     testGroup = FALSE,
                     filter = 0.1,
                     filterFreq = 0.05,
                     ...) {
    if(ncol(x)!=length(group))
        stop("number of individuals doesn't match")
    if(any(is.na(group)))
        stop("'group' variable contains missing values")
    if(filter > 0){
        m <- apply(x > 0, 1, mean)
        nfilt <- sum(m < filter, na.rm = TRUE)
        if(nfilt > 0){
            x <- x[m >= filter,]
            warning(nfilt, " exons filtered out due to absolute frequency filter")
        }
    }
    if(filterFreq > 0){
        ## Filter exons with mean relative frequency < 'filterFreq'
        auxa1<-colSums(x)
        auxR <- sweep(x, 2, auxa1, FUN="/")
        mrf <- apply(auxR, 1, mean)
        nffilt <- sum(mrf < filterFreq, na.rm = TRUE)
        if(nffilt > 0){
            x <- x[mrf >= filter,]
            warning(nffilt, " exons filtered out due to relative frequency filter")
        }
    }
    all.zero <- which(apply(x, 2, function(x) all(x==0)))
    if (length(all.zero) >= 1) {
        x <- x[,-all.zero]
        group <- group[-all.zero]
        if(any(table(group) <= 1)){
            warning("At least two non-null individuals are needed for each of the group levels")
            return(rep(NA, nlevels(group)+2))
        }
    }
    group <- as.factor(group)
    n <- ncol(x)
    nS <- nrow(x)
    method <- charmatch(test, c("asymptotic", "permutation"), nomatch=NA)
    if (is.na(method))
        stop(" 'test' argument should be 'asymptotic' or 'permutation'")

    nlev <- nlevels(group)
    
    if(nS == 1 || nS>maxTrans) {
        if(testGroup)
            out <- rep(NA, 2*(nlev + 1))
        else
            out <- rep(NA, nlev + 2)
    }
    else {
        a1<-colSums(x)
        R <- sweep(x, 2, a1, FUN="/")
        if (transformation) {
            R <- suppressWarnings(t(DR_data(t(R))))
        }
        ## Compute Hellinger distance
        d <- DistH(R)
        d.dist <- as.dist(d)

        if (method == 1)
            {
                ## When most of the individuals are located in the edges of the simplex                ## and n=2, betadisper sometimes crashes
        ## Compare dispersions
                mod <- try(betadisper(d.dist, group, type=type), silent=T)
                if(!inherits(mod,"try-error")) {
                    temp <- try(anova(mod), silent=TRUE)
                    if(!inherits(temp,"try-error")) 
                        p.homog <- temp[1,5]
                    else
                        p.homog <- NA   
                }
                else {
                    p.homog <- NA         
                }

                ## Compare means
                ##     Implement our own SS computation for the sake of efficiency:
                ssTotal <- SS(d.dist, " ", " ")
                ssResi <- lapply(levels(group), SS, d = d.dist, group = group)
                ssRes <- sum(unlist(ssResi))
                SSg <- ssTotal - ssRes
                Chisq <- SSg/(nlev-1)
                e <- try(eigenG(d), TRUE)
                eigenStats  <- c(length (e), sum(e>0), sum(e<0))
                if (eigenStats [3] > 0)  # Maybe not necessary. Depending on Hellinger curvature 
                    e <- abs(e)
                if (length(e) < nS)
                    e <- c(e, rep(0, nS - length(e))) 
                pvalue <- imhof(Chisq*(nlev-1), lambda=e, h=rep(nlev-1,nS))$Qq
                if(pvalue < 0)
                    pvalue <- liu(Chisq*(nlev-1), lambda=e, h=rep(nlev-1,nS))
            }
        if (method == 2) {
                                        # Compare dispersion
            mod <- try(betadisper(d.dist, group, type=type), silent=T)
            if(!inherits(mod,"try-error")) {
                perm   <- permutest(mod, pairwise = pairwise, control = permControl(nperm = 999))
                if (perm[1]$tab[1,6] < 0.01){
                    perm   <- permutest(mod, pairwise = pairwise, control = permControl(nperm = Nperms))
                }   
                p.homog <- perm[1]$tab[1,6]
            }
            else{
                p.homog <- NA
            }   
            
                                        # Compare means 
            mod.ad <- try(adonis(d.dist ~ group, permutations=999), silent=T)
            if(!inherits(mod.ad,"try-error")) {
                if (mod.ad$aov.tab[1,6] < 0.01 ){
                    mod.ad <- try(adonis(d.dist ~ group, permutations=Nperms), silent=T)
                } 
                pvalue <- mod.ad$aov.tab[1,6]
            }   
            else{
                pvalue <- NA
            }   
        }
        
        ## Calculate d-bar for each population
        if(!inherits(mod,"try-error")) {   
            d.bar <- tapply(mod$distances, mod$group, mean)
        }
        else {
            d.bar <- rep(NA, nlev) 
        }   
        names(p.homog) <- "p.homogeneity"
        names(pvalue) <- "pvalue"
        
        if(nlev > 2 & testGroup){ ## Calculate group-specific p-values
            resPT <- permTestSG(d, group, NpermsSG)
            out <- c(d.bar, p.homog, pvalue, resPT)
            names(out)[(length(out)-nlevels(group)+1):length(out)] <-
                paste0("pvalue_", levels(group))
        }
        else
            out <- c(d.bar, p.homog, pvalue)
    }
    
    out
}
