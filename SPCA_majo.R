#####################################################################
##### Layered PCA by majorization
#####  2018/4/5
#####
#####     X = FA' = F(A1+A2+...+Am)'
####      X : n*p data matrix.
####      F : n*r score matrix
####      Am: p*q matrix of loading layer (sparsest)
####      m : number of layers   r : number of components
#####################################################################
source("https://raw.githubusercontent.com/nyamashita/MultipleStarts/master/MultipleStart_nosnow.R")

LPCA_majo <- function(X, #data matrix; objects*variable
                      r, #number of components
                      m, #number of layers
                      itemax=1000, #max of iteration
                      eps=1e-7){

  N <- nrow(X)
  P <- ncol(X)
  
  #objective function to be minimized
  LS <- function(F,A){
    sum((X - F%*%t(A))^2)
  }
  
  
  #initial values: PCA of random matrix Xmat
  Xmat <- matrix(rnorm(N*P),N,P)
  res <- svd(Xmat)
  F <- sqrt(N)*res$u[,1:r]
  A <- res$v[,1:r]%*%diag(res$d[1:r])/sqrt(N)
  Alayers <- array(0,c(P,r,m))
  for(l in 1:m){
    Alayers[,,l] <- A/m
  }
  
  history <- c()
  
  #iteration start
  for(ite in 1:itemax){
    
    #F-step: svd of sqrt(N)*XA
    res <- svd(sqrt(N)*X%*%A)
    F <- sqrt(N)*res$u[,1:r]%*%t(res$v[,1:r])
    
    #A-step: majorization
    for(l in 1:m){
      X_tilde <- X - (F%*%t(A) - F%*%t(Alayers[,,l]))
      Q_l <- t(X_tilde)%*%F/N
      Anew <- matrix(0,P,r)
      for(row in 1:P){ #for rach row of layer
        Anew[row,which.max(Q_l[row,]^2)] <- Q_l[row,which.max(Q_l[row,]^2)]
      }
      Alayers[,,l] <- Anew
      A <- apply(Alayers,c(1,2),sum)
    }
    
    history[ite] <- LS(F,A)
    
    if(ite>1){
      if((history[ite-1] - history[ite]) < eps){break}
     }
  }
  
  list(F=F,
       A=A,
       Alayers=Alayers,
       LOSS_MIN = min(history),
       HISTORY = history)

}

#example
#X <- matrix(rnorm(1000), 100, 10)
#res <- MULTIPLE_STARTS("LPCA_majo(X,2,3)",100)