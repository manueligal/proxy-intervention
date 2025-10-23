#Likelihood function for causal parametrization
Lred <- function(theta,tab_source,tab_target){
  k_Y <- dim(tab_source)[1]
  k_X <- dim(tab_source)[2]
  k_W <- k_U <- dim(tab_source)[3]
  k_E <- dim(tab_source)[4]
  
  #Dimension of the parameter that is optimized
  par_dim2 <- k_U*(k_E+k_W+1+k_X+k_W*k_Y*k_X)
  
  #P(U|E)
  PUE <- apply(matrix(theta[1:(k_E*k_U)],k_U,k_E),2,softmax)
  
  #P(W|U)
  PWU <- apply(matrix(theta[(k_E*k_U+1):(k_U*(k_E+k_W))],k_W,k_U),2,softmax)
  
  #Q(U)
  QU <- softmax(theta[(k_U*(k_E+k_W)+1):(k_U*(k_E+k_W+1))])
  
  #P(X|U)
  PXU <- apply(matrix(theta[(k_U*(k_E+k_W+1)+1):(k_U*(k_E+k_W+1+k_X))],k_X,k_U),2,softmax)
  
  #P(Y|U,W,X)
  PYUWX <- apply(array(theta[(k_U*(k_E+k_W+1+k_X)+1):par_dim2],c(k_Y,k_U,k_W,k_X)),2:4,softmax)

  #Q(W)
  QW <- as.numeric(PWU%*%QU)
  
  #Part of the log-likelihood corresponding to the target domain
  L <- -sum(tab_target*log(QW))
  
  PYXWE <- array(NA,c(k_Y,k_X,k_W,k_E))
  for(E in 1:k_E){
    #All the values in this loop correspond to a specific domain
    #P(U)
    PU <- PUE[,E]
    
    for(x in 1:k_X){
      #P(X=x)
      Px <- as.numeric(PXU[x,]%*%PU)
      
      #P(U|X=x)
      PUx <- PXU[x,]*PU/Px
      
      #P(Y,W,X=x)
      for(s in 1:k_W){
        PYXWE[,x,s,E] <- PYUWX[,,s,x]%*%diag(PUx)%*%t(PWU)[,s]*Px
      }
    }
  }
  
  #Part of the log-likelihood corresponding to the source domains
  L <- L-sum(tab_source*log(PYXWE))
  
  return(ifelse(is.finite(L),L,1e10))
}

#Function h to define the reduced estimator
h <- function(eta,k_W,k_E){
  #Construction of the matrices Q(W), P(W|E,x) and P(y|E,x) from eta
  QW <- c(eta[1:(k_W-1)]/eta[k_W],1-sum(eta[1:(k_W-1)]/eta[k_W]))

  PWEx_num <- matrix(eta[(k_W+1):(k_W+(k_W-1)*k_E)],k_W-1,k_E)
  PWEx_den <- eta[(k_W+k_W*k_E+1):(k_W+(k_W+1)*k_E)]
  PWEx0 <- PWEx_num%*%diag(noNaN(1/PWEx_den))
  PWEx <- rbind(PWEx0,1-colSums(PWEx0))

  PyEx_num <- eta[(k_W+(k_W-1)*k_E+1):(k_W+k_W*k_E)]
  PyEx <- PyEx_num*noNaN(1/PWEx_den)
  
  #Estimation of the causal effect
  if(kappa(PWEx)<1e14){
    estim_effect <- as.numeric(PyEx%*%pseudosolve(PWEx)%*%QW)
  }else{
    estim_effect <- 1
  }
  
  return(estim_effect)
}
