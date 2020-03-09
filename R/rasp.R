rasp <- function(x, formula, expressionCols, geneidCol, filterFreq=0.05, transform = "none", type="II",
                 contrasts = NULL, mc.cores=1){

  filt <- TRUE # to be supplied  ....  related to filterFreq 
  
  if(inherits(x, "DEXSeqDataSet") | inherits(x, "RangedSummarizedExperiment")){
    
    if (missing(formula) & inherits(x, "DEXSeqDataSet")){
        warning("'formula' was not specified, the design of 'x' will be used instead")
        y <- model.matrix(design(x), data=colData(x))
    }
      
    else if(!missing(formula) & inherits(x, "DEXSeqDataSet"))  {
       y <- model.matrix(formula, data=colData(x))
    }
    else{
     stop ("formula must be provided in the argument 'formula'")        
    }
    x2 <- x
    xx <- split(x, droplevels(as.factor(geneIDs(x2)[filt])))
}

    
masterDesc <- try(get("masterDescriptor",envir = getNamespace("parallel")), TRUE)
if (mc.cores>1){
    if (class(masterDesc) == "try-error")
        stop("It appears you are trying to use multiple cores from Windows, this is not possible")
    mclapp <- try(get("mclapply", envir = getNamespace("parallel")), TRUE)
    detectCor <- try(get("detectCores", envir = getNamespace("parallel")),TRUE)
    nAvailableCores <- detectCor()
    coreID <- mclapp(as.list(1:mc.cores), function(x) masterDesc(),
                     mc.cores = mc.cores)[[1]]
    xx<-xx[mclapp(xx,nrow,mc.cores=mc.cores)>1]
    xx<-mclapp(xx,t,mc.cores=mc.cores)
    n<-length(xx)
    1
    cl <- parallel::makeForkCluster(mc.cores)
    doParallel::registerDoParallel(cl)
    p.values<-foreach(i = 1:n, .combine="rbind" ) %dopar%{
        ans <- mlm(xx[[i]]~y,type=type)
        ans[[2]][1,6]
    }
    parallel::stopCluster(cl)
}
else {
    xx<-xx[lapply(xx,nrow)>1]
    xx<-lapply(xx,t)
    n<-length(xx)
    p.values<-numeric(n)
    for (i in 1:n) {
        ans <- mlm(xx[[i]] ~ y, type=type, transform=transform, contrast=contrast, ...)
        p.values[i]<-ans[[2]][1,6]
    }
}
colnames(p.values)<-"p.value"
out<-cbind("gene_id"=names(xx), p.values)
class(out)<-c("rasp", "matrix")
out
}

