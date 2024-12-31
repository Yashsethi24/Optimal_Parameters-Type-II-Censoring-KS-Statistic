library('sde')

#Function for calculating KS statistic value
Type<-function(n,itr,mu,sd)
{
  B<-array(0,dim = c(itr,2))                             
  
  for(i in 1:itr)
  { 
    x0<-rlnorm(n)
    x<-exp(log(x0)*sd+mu)
    d<-log(x)
    estm<-sum(d)/n
    estv<-(sum((d-estm)^2))/(n-1)
    estsd<-sqrt(estv)
    
  
    B[i,1]<-sqrt(n)*(ks.test(x,"plnorm",estm,estsd)$statistic)
    
    x0<-rlnorm(n)
    d0<-log(x0)
    estm0<-sum(d0)/n
    estv0<-(sum((d0-estm0)^2))/(n-1)
    estsd0<-sqrt(estv0)
    
    
    B[i,2]<-sqrt(n)*(ks.test(x0,"plnorm",estm0,estsd0)$statistic)   # null distribution 
    
  
    
  }
  
  B
  
}


###########################################################


#Main program begins:

mu<-7                                     
sd<-2                                     
itr<-10000                               # Number of iterations
fnitr<-100                                 # final  No. of iterations
qt<-c(0.95)
N<-500  # Sample size


ksdiff<-array(0,dim = c(fnitr,2))          
qtdiff<-array(0,dim = c(fnitr,length(qt))) 
qtval<-array(0,dim = c(fnitr,2*length(qt)))


####################Quantile analysis##########################

for(i in 1 : fnitr) {
  
  
  kstestval<-Type(N,itr,mu,sd)                   #computed values of KS statistic 
  qtval[i,1]<-quantile(kstestval[,1],qt)
  qtval[i,2]<-quantile(kstestval[,2],qt)
  ksqt<-wilcox.test(kstestval[,1],kstestval[,2],alternative ="two.sided")
  ksdiff[i,]<-c(ksqt$statistic, ksqt$p.value)
  
  cat(i, ";  ")
}

############## Plotting#################


hist(ksdiff[,2], probability = T)          # pvalue

#hist(qtdiff[,1], probability = T)          # quantile difference



