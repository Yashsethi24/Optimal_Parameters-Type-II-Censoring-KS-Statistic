# Function to calculate Kolmogorov statistic test value
library('sde')
Type<-function(N,itr,co,T0)       #T0 stopping time for Type1 censoring
{ 
  Bfin<-array(0, dim=c(itr,3))
  B<-c()                          #For complete data
  Bt1<-c()                        #For Type I
  Bt2<-c()                        #For Type II
  for(i in 1:itr)
  {
    d<-c()                        #Data  
    ecm<-c()                      #empirical lower for complete data
    ecp<-c()                      #empirical upper for complete data
    cd<-c()                       #Cumulative 
    dm<-c()                       #D-
    dp<-c()                       #D+
    D<-c()                        #D
    
    ecmt1<-c()                      #empirical lower for Type I censored data
    ecpt1<-c()                      #empirical upper for Type I censored data
    cdt1<-c()                       #Cumulative for Type I censored
    dmt1<-c()                       #D- for Type I censored
    dpt1<-c()                       #D+ for Type I censored
    Dt1<-c()                        #D  for Type I censored
    
    ecmt2<-c()                    #empirical lower for Type II censored data
    ecpt2<-c()                    #empirical upper for Type II censored data
    cdt2<-c()                     #Cumulative for Type II censored
    dmt2<-c()                     #D- for Type II censored
    dpt2<-c()                     #D+ for Type II censored
    Dt2<-c()                      #D  for Type II censored
    n<-N                          #No. of obs
    censor_obj<-co                #censoring objects
    censor_time<-T0               #censoring time
    d<-sort(runif(n,0,1))         #Data
    ecm<-((1:n)-1)/n   
    ecp<-(1:n)/n
    cd=d
    dp=abs(ecp-cd)
    dm=abs(ecm-cd)
    D<-max(max(dm),max(dp))
    
    g<-which(d<=censor_time)       #Index of Non-Censored data of Type I
    non_censoredt1<-sort(d[g])    #Non-Censored data of Type I
    lt1<-length(non_censoredt1)   #Length of Type I data
    ecmt1<-ecm[1:lt1]              
    ecpt1<-ecp[1:lt1]
    cdt1<-non_censoredt1
    dpt1<-abs(cdt1-ecpt1)          
    dmt1<-abs(cdt1-ecmt1)
    d3t1<-abs(censor_time-lt1/n) 
    Dt1<-max(dpt1,dmt1,d3t1)
    
    non_censoredt2<-d[1:censor_obj] #Non-Censored data
    lt2<-length(non_censoredt2)
    ecmt2<-ecm[1:lt2]  
    ecpt2<-ecp[1:lt2]             
    cdt2<-non_censoredt2         
    dpt2<-abs(cdt2-ecpt2)          
    dmt2<-abs(cdt2-ecmt2)
    Dt2<-max(dpt2,dmt2)
    B[i]<-sqrt(n)*D
    Bt1[i]<-sqrt(n)*Dt1
    Bt2[i]<-sqrt(n)*Dt2
  }
  
  Bfin[,1]<- B
  Bfin[,2]<- Bt1
  Bfin[,3]<- Bt2 
 
  Bfin 
}



#Function to calculate Kolmogrov distribution value 
kolcomp2<-function(N, brsize, itr, co, T0){
  bbmax<-array(0, dim = c(itr,4))
  
  for (i in 1:itr){
    br<-(BBridge(x=0, y=0, t0=0, T=1, N=brsize))
    brabs<- abs(br)
    rab<-rbeta(1,shape1 =  co,shape2 =  (N-co+1)) 
    u<-seq(0,1,length=(brsize+1) )
    bbmax[i,1]<-max(brabs)
    bbmax[i,2]<-max(brabs[1:floor(T0*brsize)])
    bbmax[i,3]<-max(brabs[1:floor(rab*brsize)]) # max(brabs[1:floor((co/N)*brsize)])
    bbmax[i,4]<-max(abs(sqrt((co-1)/N)*br+u*((co-1)/N-rab)*sqrt(N)), abs(co/N-rab))
    
   
  }
  
  bbmax
}  


###########################################################


#Main program begins:
N<-200                                #Sample size
itr<-10000                            # No. of iterations
brsize<-N*100                         # Grid size of Brownian bridge
fnitr<-1000                           # final  No. of iterations
qt<-c(0.95)                           #c( seq(0.1,0.9,by=0.1), seq(0.91,0.99,by=0.01)) # probability 
co<-80                                #Censored objects
T0<-0.4                               #Censoring Time


ksdiff0<-array(0,dim = c(fnitr,2))  # 2 sample  KS test values   for complete data
ksdiff1<-array(0,dim = c(fnitr,2))  # 2 sample  KS test values   for Type I
ksdiff2a<-array(0,dim = c(fnitr,2)) # 2 sample  KS test values   using (co/N)
ksdiff2b<-array(0,dim = c(fnitr,2)) # 2 sample  KS test values   using beta 
ksdiff2c<-array(0,dim = c(fnitr,2)) # 2 sample  KS test values   using function beta 

qtdiff0<-array(0,dim = c(fnitr,length(qt))) # quantile difference for complete data
qtdiff1<-array(0,dim = c(fnitr,length(qt))) # quantile difference for Type I
qtdiff2a<-array(0,dim = c(fnitr,length(qt))) # quantile difference using (co/N)
qtdiff2b<-array(0,dim = c(fnitr,length(qt))) # quantile difference using beta
qtdiff2c<-array(0,dim = c(fnitr,length(qt))) # quantile difference using  function of beta
qtval<-array(0,dim = c(fnitr,8*length(qt)))

####################Quantile analysis##########################

for(i in 1 : fnitr) {
  brval<-kolcomp2(N, brsize, itr, co, T0 )       #Kolmogorv distribution values  4 values 
  kstestval<-Type(N,itr,co,T0)                   #computed values of KS statistic 3 values 
  
  qtval[i,1]<-quantile(brval[,1],qt)
  qtval[i,2]<-quantile(brval[,2],qt)
  qtval[i,3]<-quantile(brval[,3],qt)
  qtval[i,4]<-quantile(brval[,4],qt)
  qtval[i,5]<-quantile(kstestval[,1],qt)
  qtval[i,6]<-quantile(kstestval[,3],qt)
  qtval[i,7]<-qtval[i,6]
  qtval[i,8]<-qtval[i,6]
  
  qtdiff0[i, ]<-qtval[i,5]-qtval[i,1]
  qtdiff2a[i, ]<-qtval[i,6]-qtval[i,2]
  qtdiff2b[i, ]<-qtval[i,7]-qtval[i,3]
  qtdiff2c[i, ]<-qtval[i,8]-qtval[i,4]
  
  ksqt0<-(ks.test(kstestval[,1],brval[,1],alternative = "two.sided"))
  ksdiff0[i,]<-c(ksqt0$statistic, ksqt0$p.value)
  
  ksqt2a<-(ks.test(kstestval[,2],brval[,2],alternative = "two.sided"))
  ksdiff1[i,]<-c(ksqt2a$statistic, ksqt2a$p.value)
  
  ksqt2b<-(ks.test(kstestval[,3],brval[,3],alternative = "two.sided"))
  ksdiff2b[i,]<-c(ksqt2b$statistic, ksqt2b$p.value)
  
  ksqt2c<-(ks.test(kstestval[,3],brval[,4],alternative = "two.sided"))
  ksdiff2c[i,]<-c(ksqt2c$statistic, ksqt2c$p.value)
  
 cat(i, date() ,"\n ")
}

write.csv (qtval, "qtval.csv"   )
cat("=============Complete ===========","\n")
hist(ksdiff0[,1],  probability = T) # D value for Complete Data
hist(ksdiff0[,2], probability = T) # pvalue  for Complete Data
hist(qtdiff0[,1], probability = T) #  quantile difference for Complete Data
print(var(qtdiff0[,1]))
print(median(qtdiff0[,1]))

cat("========= Type IIa===============","\n")
hist(ksdiff2a[,1],  probability = T) # D value for Type II Data
hist(ksdiff2a[,2], probability = T) # pvalue  for Type II Data
hist(qtdiff2a[,1], probability = T) #  quantile difference for Type II Data
print(var(qtdiff2a[,1]))
print(median(qtdiff2a[,1]))

cat("========= Type IIb===============","\n")
hist(ksdiff2b[,1],  probability = T) # D value for Type II Data
hist(ksdiff2b[,2], probability = T) # pvalue for Type II Data
hist(qtdiff2b[,1], probability = T) #  quantile difference for Type II Data
print(var(qtdiff2b[,1]))
print(median(qtdiff2b[,1]))


cat("========= Type I===============","\n")
hist(ksdiff2c[,1],  probability = T) # D value for Type II Data
hist(ksdiff2c[,2], probability = T) # pvalue   for Type II Data
hist(qtdiff2c[,1], probability = T) #  quantile difference for Type II Data
print(var(qtdiff2c[,1]))
print(median(qtdiff2c[,1]))

plot(ecdf(qtdiff0[,1]))
lines(ecdf(qtdiff1[,1]),col=2)
lines(ecdf(qtdiff2a[,1]),col=3)
lines(ecdf(qtdiff2b[,1]),col=4)
