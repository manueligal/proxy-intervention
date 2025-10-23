source('Lfunction.R')

#Estimate the causal effect using the reduced parameterization
ReducedEffectEstimation <- function(Sample,y,doX){
  k_W <- max(Sample[,'W'])
  k_Y <- max(Sample[,'Y'])
  
  #Selection of the observations where X=doX
  SampledoX <- Sample[Sample[,'X']==doX,]
  
  #Estimation of the matrices
  QW <- as.vector(table(factor(Sample[which(Sample[,'E']==(k_E+1)),'W'],levels=1:k_W))/sum(Sample[,'E']==(k_E+1)))
  
  PWEx_n0 <- table(factor(SampledoX[,'W'],levels=1:k_W),SampledoX[,'E'])
  PWEx_n <- PWEx_n0[,-which(colnames(PWEx_n0)==as.character(k_E+1))]
  PWEx <- sweep(PWEx_n,2,colSums(PWEx_n),'/')
  
  PYEx_n0 <- table(factor(SampledoX[,'Y'],levels=1:k_Y),SampledoX[,'E'])
  PYEx_n <- PYEx_n0[,-which(colnames(PYEx_n0)==as.character(k_E+1))]
  PyEx <- sweep(PYEx_n,2,colSums(PYEx_n),'/')[y,]
  
  #Estimation of the causal effect
  if(kappa(PWEx)<1e14){
    estim_effect <- min(max(as.numeric(PyEx%*%pseudosolve(PWEx)%*%QW),0),1)
  }else{
    estim_effect <- NA
  }
  
  return(list(estim_effect=estim_effect,PyEx=PyEx,PWEx=PWEx,QW=QW))
}

#Estimate the causal effect using the causal parameterization
CausalEffectEstimation <- function(Sample,y,doX,Nseeds){
  #Calculate the dimensions
  k_Y <- max(Sample[,'Y'])
  k_U <- k_W <- max(Sample[,'W'])
  k_E <- max(Sample[,'E'])-1
  k_X <- max(Sample[,'X'])
  
  #Dimension of the parameter that is optimized
  par_dim2 <- k_U*(k_E+k_W+1+k_X+k_W*k_Y*k_X)
  
  #Table of the number of observations n(y,x,w,e)
  tab_source <- table(Sample[,'Y'],Sample[,'X'],Sample[,'W'],Sample[,'E'])[,,,1:k_E]
  tab_target  <- table(Sample[Sample[,'E']==(k_E+1),'W'])
  
  theta_unr_opt <- NA
  L_opt <- 1e10
  for(seed in 1:Nseeds){
    opt <- optim(runif(par_dim2),Lred,method='L-BFGS-B',tab_source=tab_source,tab_target=tab_target,control=list(maxit=5e4))
    if(opt$value<L_opt){
      L_opt <- opt$value
      theta_unr_opt <- opt$par
    }
  }
  
  #Transformation to keep the components in [0,1]
  pars  <- theta_unr_opt
  
  #P(U|E) optimal
  PUEo <- apply(matrix(pars[1:(k_E*k_U)],k_U,k_E),2,softmax)
  
  #P(W|U) optimal
  PWUo <- apply(matrix(pars[(k_E*k_U+1):(k_U*(k_E+k_W))],k_W,k_U),2,softmax)

  #Q(U) optimal
  QUo <- softmax(pars[(k_U*(k_E+k_W)+1):(k_U*(k_E+k_W+1))])
  
  #P(X|U) optimal
  PXUo <- apply(matrix(pars[(k_U*(k_E+k_W+1)+1):(k_U*(k_E+k_W+1+k_X))],k_X,k_U),2,softmax)
  
  #P(y|U,W,doX) optimal
  PYUWXo <- apply(array(pars[(k_U*(k_E+k_W+1+k_X)+1):par_dim2],c(k_Y,k_U,k_W,k_X)),2:4,softmax)
  PyUWxo <- PYUWXo[y,,,doX]
  
  #Estimation of the causal effect
  estim_effect <- as.numeric(diag(PyUWxo%*%PWUo)%*%QUo)
  
  return(list(estim_effect=estim_effect,pars=pars,L=L_opt))
}
