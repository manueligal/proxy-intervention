library(ggplot2)
library(reshape2)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)

#Parameters for the number of samples and sample size
M <- 10
N <- 5
n <- 2e4

#Number of categories
k_W <- 2
k_E <- 3
k_Y <- 2
k_X <- 2

#Study the causal effect q(Y=y|do(X=doX)). The value in () is in binary notation
doX <- (1)+1
y <- (1)+1

set.seed(126)
results <- data.frame(matrix(NA,M*N,7))
colnames(results) <- c('Oracle','Causal','Reduced','NoAdj','NoAdj*','WAdj','WAdj*')

for(m in 1:M){
  #Select only those cases with certain level of confounding
  diff <- 0
  while(diff<0.1){
    matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)
    effect <- as.numeric(matrices$effect)
      
    QyUx <- diag(matrices$PYUWX[y,,,doX]%*%matrices$PWU)
    Qx <- as.numeric(matrices$PXU[doX,]%*%matrices$QU)
    obs_effect <- QyUx%*%(matrices$PXU[doX,]*matrices$QU)/Qx
      
    diff <- abs(effect-obs_effect)
  }
  
  for(j in 1:N){
    #Generation of the sample
    Sample <- SampleGeneration(matrices,n)
    SampledoX <- Sample[Sample[,'X']==doX,]
    
    #Estimate for the reduced parameterization
    results[(m-1)*N+j,'Reduced'] <- ReducedEffectEstimation(Sample,y,doX)$estim_effect-effect
    
    #Estimate for the causal parameterization
    results[(m-1)*N+j,'Causal'] <- CausalEffectEstimation(Sample,y,doX,3)$estim_effect-effect
    
    #Estimate from the intervention distribution (oracle)
    Sample_oracle <- Sample
    Sample_oracle[,'X'] <- rep(doX,n)
    Sample_oracle[,'Y'] <- mapply(substitution2,Sample_oracle[,'U'],Sample_oracle[,'W'],Sample_oracle[,'X'],MoreArgs=list(P=matrices$PYUWX))
    results[(m-1)*N+j,'Oracle'] <- mean(Sample_oracle[Sample_oracle[,'E']==(k_E+1),'Y']==y)-effect
    
    #Estimate using the pooled conditional distribution
    results[(m-1)*N+j,'NoAdj'] <- mean(SampledoX[SampledoX[,'E']!=(k_E+1),'Y']==y)-effect
    
    #Estimate using the conditional distribution in the target domain
    results[(m-1)*N+j,'NoAdj*'] <- mean(SampledoX[SampledoX[,'E']==(k_E+1),'Y']==y)-effect
    
    #Estimate using the (not valid) adjustment set W
    #With the pooled data from source domains
    PYWx_n_p <- table(SampledoX[SampledoX[,'E']!=(k_E+1),'Y'],SampledoX[SampledoX[,'E']!=(k_E+1),'W'])
    PyWx_p <- sweep(PYWx_n_p,2,colSums(PYWx_n_p),'/')[y,]
    QW_p <- as.vector(table(Sample[which(Sample[,'E']!=(k_E+1)),'W'])/sum(Sample[,'E']!=(k_E+1)))
    results[(m-1)*N+j,'WAdj'] <- PyWx_p%*%QW_p-effect
    
    #With the target domain data
    PYWx_n <- table(SampledoX[SampledoX[,'E']==(k_E+1),'Y'],SampledoX[SampledoX[,'E']==(k_E+1),'W'])
    PyWx <- sweep(PYWx_n,2,colSums(PYWx_n),'/')[y,]
    QW <- as.vector(table(Sample[which(Sample[,'E']==(k_E+1)),'W'])/sum(Sample[,'E']==(k_E+1)))
    results[(m-1)*N+j,'WAdj*'] <- PyWx%*%QW-effect
  }
}

results['sample'] <- 1:(N*M)
results_melt <- melt(results,id.vars='sample',variable.name='method',value.name='error')

ggplot(results_melt,aes(x=method,y=abs(error))) +
  geom_boxplot() +
  theme_bw() +
  labs(x='Method',y='Absolute error') +
  ylim(0,0.35)
