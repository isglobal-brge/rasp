
DistH <- function(data) {
    k <- ncol(data)
    interdist <- matrix(data=0,nrow=k,ncol=k)
    for(i in 2:k) {
      for(j in 1:(i-1)) {
        interdist[i,j] <- hellingerDist(data[,i],data[,j])
      }
    }
    return(interdist)
}


hellingerDist <-
function (x1, x2) {
   a <- (sqrt(x1) - sqrt(x2))
#   b <- sqrt(sum(a*a)) 
   b <- sqrt(a%*%a)	
   return(b)
}


