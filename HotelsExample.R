library(ggplot2)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)
source('utils/FunctionsCI.R',chdir=TRUE)

#Load the data
train <- read.csv('../data/train.csv', header=TRUE)

#Transform the price into a categorical variable
train$price_usd_cat <- as.numeric(cut(train$price_usd,breaks=c(0,75,125,175,225,1e8)))

#Obtain the matrix of the sample of (E,W,X,Y)
sub <- train[,c('prop_id','price_usd_cat','position','click_bool','random_bool')]
sub <- sub[!is.na(sub$price_usd_cat),]
sub$click_bool <- sub$click_bool+1

#Divide into the observational and interventional datasets
obs <- sub[sub$random_bool==0,1:4]
exp <- sub[sub$random_bool==1,1:4]

#Keep only hotels with a minimum sample size
#Observational dataset
hotels_obs <- table(obs$prop_id)
hotels_selection_obs <- as.double(rownames(hotels_obs[hotels_obs>2000]))

#Randomized dataset
hotels_exp <- table(exp$prop_id)
hotels_selection_exp <- as.double(rownames(hotels_exp[hotels_exp>1500]))

#Select the source domains
hotels_kept <- setdiff(hotels_selection_obs,hotels_selection_exp)
obs_selection <- obs[obs$prop_id%in%hotels_kept,]
obs_selection$prop_id <- factor(obs_selection$prop_id)
levels(obs_selection$prop_id) <- 1:length(hotels_kept)
obs_selection$prop_id <- as.numeric(obs_selection$prop_id)

#Selection of parameters
k_E <- length(hotels_kept)
k_Y <- 2
k_X <- 2
k_W <- length(unique(obs$price_usd))
doX <- (1)+1
y <- (1)+1


#Estimate the causal effect for some holdout target domains
results <- data.frame(matrix(NA,length(hotels_selection_exp),12))
colnames(results) <- c('estim','ci1','ci2','exp','ci3','ci4','cond','ci5','ci6','proxy')

for(i_id in 1:length(hotels_selection_exp)){
  id <- hotels_selection_exp[i_id]
  
  #Observational dataset corresponding to the target domain
  target <- obs[obs$prop_id==id,]
  target$prop_id <- k_E+1
  Sample <- rbind(obs_selection,target)
  Sample$position <- (Sample$position==1)*1+1
  colnames(Sample) <- c('E','W','X','Y')

  #Estimation of the causal effect
  estim_effect <- ReducedEffectEstimation(Sample,y,doX)$estim_effect
  results[i_id,'estim'] <- estim_effect
  
  #Calculation of the causal effect from the interventional data
  target_exp <- exp[exp$prop_id==id,]
  target_exp$position <- (target_exp$position==1)*1+1
  colnames(target_exp) <- c('E','W','X','Y')
  estim_exp <- mean(target_exp[target_exp$X==doX,'Y']==y)
  n_exp <- length(target_exp[target_exp$X==doX,'Y'])
  results[i_id,'exp'] <- estim_exp

  #Calculation of the causal effect using the conditional distribution
  target$position <- (target$position==1)*1+1
  colnames(target) <- c('E','W','X','Y')
  estim_cond <- mean(target[target$X==doX,'Y']==y)
  n_cond <- length(target[target$X==doX,'Y'])
  results[i_id,'cond'] <- estim_cond
  
  #Calculation of the causal effect using proxy adjustment
  SampledoX <- Sample[Sample[,'X']==doX,]
  PYWx_n <- table(SampledoX[SampledoX[,'E']==(k_E+1),'Y'],factor(SampledoX[SampledoX[,'E']==(k_E+1),'W'],levels=1:k_W))
  PyWx <- sweep(PYWx_n,2,colSums(PYWx_n),'/')[y,]
  QW <- as.vector(table(Sample[which(Sample[,'E']==(k_E+1)),'W'])/sum(Sample[,'E']==(k_E+1)))
  results[i_id,'proxy'] <- noNaN(PyWx)%*%QW

  sigma <- ReducedSigma(Sample,y,doX)
  results[i_id,c('ci1','ci2')] <- estim_effect+c(-1,1)*qnorm(0.975)*sigma
  
  results[i_id,c('ci3','ci4')] <- estim_exp+c(-1,1)*qnorm(0.975)*sqrt(estim_exp*(1-estim_exp)/n_exp)
  
  results[i_id,c('ci5','ci6')] <- estim_cond+c(-1,1)*qnorm(0.975)*sqrt(estim_cond*(1-estim_cond)/n_cond)
}

results2 <- cbind(id=1:nrow(results),results)
results_melt0 <- cbind(results2[,c(1,2:4)],method=rep('Reduced',length(hotels_selection_exp)))
results_melt1 <- cbind(results2[,c(1,5:7)],method=rep('Oracle',length(hotels_selection_exp)))
results_melt2 <- cbind(results2[,c(1,8:10)],method=rep('NoAdj*',length(hotels_selection_exp)))
results_melt3 <- cbind(results2[,c(1,11:13)],method=rep('WAdj*',length(hotels_selection_exp)))
colnames(results_melt1)[2:4] <- c('estim','ci1','ci2')
colnames(results_melt2)[2:4] <- c('estim','ci1','ci2')
colnames(results_melt3)[2:4] <- c('estim','ci1','ci2')
results_melt <- rbind(results_melt0,results_melt1,results_melt2,results_melt3)
results_melt$target <- rep(1:length(hotels_selection_exp),4)
results_melt$method <- factor(results_melt$method,levels=c('Reduced','Oracle','NoAdj*','WAdj*'))
results_melt[c('ci1','ci2')] <- apply(results_melt[c('ci1','ci2')],1:2,clamp)

#Pooled source domains data estimation
estim_pooled <- mean(Sample[(Sample$E!=(k_E+1))&(Sample$X==doX),'Y']==y)
n_pooled <- length(Sample[(Sample$E!=(k_E+1))&(Sample$X==doX),'Y'])
CI_pooled <- estim_pooled+c(-1,1)*qnorm(0.975)*sqrt(estim_pooled*(1-estim_pooled)/n_pooled)

PYWx_n_p <- table(SampledoX[SampledoX[,'E']!=(k_E+1),'Y'],SampledoX[SampledoX[,'E']!=(k_E+1),'W'])
PyWx_p <- sweep(PYWx_n_p,2,colSums(PYWx_n_p),'/')[y,]
QW_p <- as.vector(table(Sample[which(Sample[,'E']!=(k_E+1)),'W'])/sum(Sample[,'E']!=(k_E+1)))
WAdj_pooled <- as.numeric(PyWx_p%*%QW_p)

#Estimate both effects for different atomic interventions
K <- 20
results_k <- data.frame(matrix(NA,K,9))
colnames(results_k) <- c('estim1','ci1','ci2','exp','ci3','ci4','cond','ci5','ci6')

#Dataset corresponding to the target domain
id <- hotels_selection_exp[1]
target <- obs[obs$prop_id==id,]
target$prop_id <- k_E+1
target_exp <- exp[exp$prop_id==id,]

for(k in 1:K){
  Sample <- rbind(obs_selection,target)
  Sample$position <- (Sample$position==k)*1+1
  colnames(Sample) <- c('E','W','X','Y')
  
  #Estimation using the reduced parameterization
  estim_effect1 <- ReducedEffectEstimation(Sample,y=2,doX=2)$estim_effect
  results_k[k,'estim1'] <- estim_effect1

  sigma1 <- ReducedSigma(Sample,y,doX)
  results_k[k,c('ci1','ci2')] <- sapply(estim_effect1+c(-1,1)*qnorm(0.975)*sigma1,clamp)
  
  #Calculation of the causal effect from the interventional data
  Sample_exp <- target_exp
  Sample_exp$position <- (Sample_exp$position==k)*1+1
  colnames(Sample_exp) <- c('E','W','X','Y')
  estim_exp <- mean(Sample_exp[Sample_exp$X==doX,'Y']==y)
  n_exp <- length(Sample_exp[Sample_exp$X==doX,'Y'])
  results_k[k,'exp'] <- estim_exp
  
  #Calculation of the causal effect using the conditional distribution
  target_c <- target
  target_c$position <- (target$position==k)*1+1
  colnames(target_c) <- c('E','W','X','Y')
  estim_cond <- mean(target_c[target_c$X==doX,'Y']==y)
  n_cond <- length(target_c[target_c$X==doX,'Y'])
  results_k[k,'cond'] <- estim_cond
  
  if(n_exp>0){
    results_k[k,c('ci3','ci4')] <- as.numeric(binom.test(estim_exp*n_exp,n_exp)$conf.int)
    results_k[k,c('ci5','ci6')] <- as.numeric(binom.test(estim_cond*n_cond,n_cond)$conf.int)
  }
}

results_k_melt1 <- cbind(results_k[,c(4:6,1:3)],method=rep('Reduced',K))
results_k_melt2 <- cbind(results_k[,c(4:6,7:9)],method=rep('NoAdj*',K))
colnames(results_k_melt2)[4:6] <- c('estim1','ci1','ci2')
results_k_melt <- rbind(results_k_melt1,results_k_melt2)
results_k_melt$target <- rep(1:K,2)
results_k_melt$method <- factor(results_k_melt$method,levels=c('Reduced','NoAdj*'))

#We calculate the proportion of clicks across all positions
average_clicks <- mean(target_exp$click_bool==y)

ggplot(data=results_melt,aes(x=target)) +
  geom_hline(aes(yintercept=estim_pooled,color='NoAdj')) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=CI_pooled[1],ymax=CI_pooled[2]),fill='gray',alpha=0.015) +
  geom_hline(aes(yintercept=WAdj_pooled,color='WAdj')) +
  geom_point(aes(y=estim,color=method),position=position_dodge(0.6)) +
  geom_errorbar(aes(ymin=(ci1),ymax=(ci2),color=method),width=1,position=position_dodge(0.6)) +
  theme_bw() +
  theme(strip.background=element_blank(),strip.text.x=element_blank()) +
  labs(x='Target Domain ID',y=expression(widehat(q)[n]~'(Y=1|do(X=1))')) +
  scale_color_manual('',breaks=c('Reduced','NoAdj*','NoAdj','Oracle','WAdj*','WAdj'),values=c('Reduced'='black','Oracle'='red','NoAdj*'='blue','NoAdj'='darkgray','WAdj*'='green','WAdj'='purple')) +
  coord_cartesian(ylim=c(-0.025,0.425))

ggplot(data=results_k_melt,aes(x=target)) +
  geom_hline(aes(yintercept=average_clicks,linetype='Avg click proportion'),color='red',alpha=0.5) +
  geom_point(aes(y=estim1,color=method),position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(ci1),ymax=(ci2),color=method),position=position_dodge(0.5),width=1) +
  geom_point(aes(y=exp,color='Oracle')) +
  geom_errorbar(aes(ymin=(ci3),ymax=(ci4),color='Oracle'),width=0.5) +
  theme_bw() +
  labs(x='Position x',y=expression(widehat(q)[n]~'(Y=1|do(X=x))')) +
  scale_color_manual('',breaks=c('Reduced','NoAdj*','Oracle'),values=c('Reduced'='black','Oracle'='red','NoAdj*'='blue')) +
  scale_linetype_manual('',values=c('Avg click proportion'='dashed')) +
  coord_cartesian(ylim=c(-0.025,0.5))

#Calculation of the MAE for the different methods
oracle <- results$exp
print(paste('Reduced:',mean(abs(results$estim-oracle))))
print(paste('NoAdj:',mean(abs(estim_pooled-oracle))))
print(paste('NoAdj*:',mean(abs(results$cond-oracle))))
print(paste('WAdj:',mean(abs(WAdj_pooled-oracle))))
print(paste('WAdj*:',mean(abs(results$proxy-oracle))))