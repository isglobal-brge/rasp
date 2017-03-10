eigenG <- function (x,tol=10^-12) {

   A <- (- 0.5) * x^2
   A <-  A + t(A)
   
   n <- ncol(A)

   I <- diag (1, nrow=n) 
   J <- matrix(1/n, ncol=n, nrow=n)
   Aux <- I - J
   
   G <- Aux %*% A %*% Aux

   e <- eigen(G, symmetric=T, only.values=T)$values
   index <- abs(e) > tol

   ans <- e[index]/n
   ans

}


