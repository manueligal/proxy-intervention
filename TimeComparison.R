source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)
library(ggplot2)

#Generation parameters
sample_sizes <- c(1e2,1e3,1e4,1e5)  
M <- 5 
N <- 2   

#Number of categories
k_E <- 2  
k_W <- 2  
k_X <- 2  
k_Y <- 2  

#Study the causal effect q(Y=y|do(X=doX)). The value in () is in binary notation
doX <- (1)+1
y <- (1)+1

sample_list <- list()
set.seed(127)

for(t in 1:length(sample_sizes)){
  n <- sample_sizes[t]
  sample_list0 <- list()
  for(m in 1:M){
    matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)
    
    for(j in 1:N){
      Sample <- SampleGeneration(matrices,n)
      SampledoX <- Sample[Sample[,'X']==doX,]
      
      Sample_oracle <- Sample
      Sample_oracle[,'X'] <- rep(doX,n)
      Sample_oracle[,'Y'] <- mapply(substitution2,Sample_oracle[,'U'],Sample_oracle[,'W'],Sample_oracle[,'X'],MoreArgs=list(P=matrices$PYUWX))
      
      sample_list0[[(m-1)*N+j]] <- list(Sample=Sample,SampledoX=SampledoX,Sample_oracle=Sample_oracle)
    }
  }
  sample_list[[t]] <- sample_list0
}

timer_method <- function(method,size_num){
  for(sample_num in 1:M*N){
    Sample <- sample_list[[size_num]][[sample_num]]$Sample
    SampledoX <- sample_list[[size_num]][[sample_num]]$SampledoX
    Sample_oracle <- sample_list[[size_num]][[sample_num]]$Sample_oracle
    
    if(method=='oracle'){
      estim <- mean(Sample_oracle[Sample_oracle[,'E']==(k_E+1),'Y']==y)
    }else if(method=='causal'){
      estim <- CausalEffectEstimation(Sample,y,doX,3)$estim_effect
    }else if(method=='reduced'){
      estim <- ReducedEffectEstimation(Sample,y,doX)$estim_effect
    }else if(method=='noadj'){
      estim <- mean(SampledoX[SampledoX[,'E']!=(k_E+1),'Y']==y)
    }else if(method=='noadj*'){
      estim <- mean(SampledoX[SampledoX[,'E']==(k_E+1),'Y']==y)
    }else if(method=='wadj'){
      PYWx_n_p <- table(SampledoX[SampledoX[,'E']!=(k_E+1),'Y'],SampledoX[SampledoX[,'E']!=(k_E+1),'W'])
      PyWx_p <- sweep(PYWx_n_p,2,colSums(PYWx_n_p),'/')[y,]
      QW_p <- as.vector(table(Sample[which(Sample[,'E']!=(k_E+1)),'W'])/sum(Sample[,'E']!=(k_E+1)))
      estim <- PyWx_p%*%QW_p
    }else if(method=='wadj*'){
      PYWx_n <- table(SampledoX[SampledoX[,'E']==(k_E+1),'Y'],SampledoX[SampledoX[,'E']==(k_E+1),'W'])
      PyWx <- sweep(PYWx_n,2,colSums(PYWx_n),'/')[y,]
      QW <- as.vector(table(Sample[which(Sample[,'E']==(k_E+1)),'W'])/sum(Sample[,'E']==(k_E+1)))
      estim <- PyWx%*%QW
    }
  }
}

methods <- c('oracle','causal','reduced','noadj','noadj*','wadj','wadj*')
times <- data.frame(method=rep(methods,length(sample_sizes)),time=rep(0,length(sample_sizes)*length(methods)),n=rep(sample_sizes,each=length(methods)))

for(t in 1:length(sample_sizes)){
  for(method_num in 1:length(methods)){
    method <- methods[method_num]
    times[(t-1)*length(methods)+method_num,'time'] <- as.numeric(system.time(replicate(5,timer_method(method,t)))[1])
  }
}

ggplot(times,aes(x=n,y=time)) +
  geom_point(aes(color=method)) +
  geom_line(aes(color=method)) +
  theme_bw() +
  scale_color_manual('',breaks=methods,,values=c('red','orange','black','gray','blue','purple','green')) +
  labs(x='Sample size',y='Time (s)',color='Method') +
  scale_x_continuous(breaks=sample_sizes,trans='log10') +
  scale_y_continuous(breaks=c(0,1,5,20,80,240),trans=scales::pseudo_log_trans(base=10, sigma=0.25))
