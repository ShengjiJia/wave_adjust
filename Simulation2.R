#######get data, please modify the file path if necessary
CGHdata <- read.csv("C:/Users/PC/Desktop/我的文档/Research/Projects/change points(wave adjust)/CGHdataset.csv")

library(grpreg)
library(KernSmooth)
library(np)
library(DNAcopy)
library(imputeTS)

#################################define functions
estimateSigma<-function (Y, h = 10) {  
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

localMax<-function (y, span = 5) {  
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

MultiScan<-function(y, h=20){
  m=dim(y)[1]
  n=dim(y)[2]
  s=matrix(0,nr=m,nc=n)
  sigma=1:m
  for (i in 1:m){
    sigma[i]=estimateSigma(y[i,], h = max(3, 2*floor(log(n))))
    s[i,]=(localDiagnostic(y[i,],h=h))^2*h/(2*sigma[i])
  }
  S=apply(s,2,sum)
  index=localMax(S,span=2*h)
  return(data.frame(S[index],index))
}

threshold<-function(y, alpha=0.05, h=20) {
  m=dim(y)[1]
  n=dim(y)[2]
  empirical=NULL
  for (k in 1:100) {
    Y=matrix(rnorm(n*m),nr=m)
    s=matrix(0,nr=m,nc=n)
    for (i in 1:m){ s[i,]=(localDiagnostic(Y[i,],h=h))^2*h/2 }
    S=apply(s,2,sum)
    index=localMax(S,span=2*h)
    empirical=c(empirical, S[index])
  }
  return(quantile(empirical, probs=1-alpha))
}

TruePositive<-function(estimate, true, error=3){         
  #calculate the number of true positives
  Num=0
  if(length(estimate)>0){
    for(l in 1:length(true)){
      if(min(abs(estimate-true[l]))<=error)
        Num=Num+1
    }
  }
  return(Num)
}

FalsePositive<-function(estimate, true, error=3){         
  #calculate the number of flase positives
  Num=0
  if(length(estimate)>0){
    for(l in 1:length(estimate)){
      if(min(abs(estimate[l]-true))>error)
        Num=Num+1
    }
  }
  return(Num)
}

screening<-function(x, y, h=10){       #local linear smoothing (rightlimit-leftlimit)
  n=length(x)
  xx=1:(n+2*h)
  yy=c(rep(y[1],h),y,rep(y[n],h))       #create data outside the boundaries
  right=rep(0,n)     #rightlimit for xx[1:n]
  left=rep(0,n)       #leftlimit for xx[(1+2*h):(n+2*h)]
  for (i in 1:n){
    model=locpoly(xx[i:(i+2*h)],yy[i:(i+2*h)],kernel="epanech",bandwidth=h,gridsize=1+2*h)
    right[i]=model$y[1]
    left[i]=model$y[1+2*h]
  }
  L=c(rep(0,h),right[(2*h+2):n]-left[1:(n-2*h-1)],rep(0,h+1))
  return(L)
}


##############################real wave estimation
y=CGHdata[2:2301,19]               
x=1:length(y)
y[which(abs(y)>2)]=NA     #delete outliers
y=na_ma(y,k=5,weighting="linear")     #imputation
h=10
x1=matrix(1,nrow=length(y),ncol=length(y))
for(i in 1:(length(y)-1)) {x1[i,(i+1):length(y)]=0} 
X1=x1
for(i in 2:length(y)) {X1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=length(y))$y }      
y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=length(y))$y
gSCAD=grpreg(X=X1[,2:length(y)], y=y0, penalty="grSCAD",family="gaussian",lambda=0.005) 
re=y-x1%*%gSCAD$beta
rwave=10*locpoly(x,re,kernel="epanech",bandwidth=20,gridsize=2000)$y
#################################Simulation2
set.seed(18)
n=2000
d=5
J=20
x=1:n
tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)
h=10
W=function(t,x,h){diag(1/h*3/4*(1-((x-t)/h)^2)*((1-((x-t)/h)^2)>0))}
D=function(t,x,h){matrix(c(rep(1,times=length(x)),(x-t)/h),nrow=length(x),ncol=2)}
s=function(t,x,h){t(c(1,0))%*%solve(t(D(t,x,h))%*%W(t,x,h)%*%D(t,x,h))%*%t(D(t,x,h))%*%W(t,x,h)}
S=NULL                   #smoothing matrix 
for(i in 1:length(x)){S=rbind(S,s(x[i],x,h))}
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) x1[i,(i+1):n]=0       
XX=kronecker(cbind(x1[,1],(diag(1,n)-S)%*%x1[ ,-1]), diag(1, d))                   
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,d))  
num=matrix(0,nr=5,nc=100)       #number of estimated change-points
num1=matrix(0,nr=5,nc=100) 
TP=matrix(0,nr=5,nc=100)        #number of true positives
TP1=matrix(0,nr=5,nc=100)  
FP=matrix(0,nr=5,nc=100)        #number of flase positives
FP1=matrix(0,nr=5,nc=100)

for(j in 1:100){
  ##############################data generation for Scenario I
  beta=matrix(sample(c(-1,-0.5,0,0.5,1), size=J*d, replace = TRUE), nrow=d)
  theta=runif(1,min=0,max=2*pi)
  signal=matrix(rep(0, n*d), nrow=d)
  y=matrix(rep(0, n*d), nrow=d)
  wave=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    wave[k, ]=sample(c(-0.5,-0.3,0.3,0.5),size=1)*sin(pi*x/50+theta)
    signal[k, ]=0
    for(i in 1:J) {signal[k, ]<-signal[k, ]+beta[k,i]*(x>tau[i])}
    y[k, ]=signal[k, ]+wave[k, ]+rnorm(n,mean=0,sd=0.2)  
  }
  ##############################method1 CBS
  CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
  estimate1=sort(unique(CBS$output[,4]))
  estimate1=estimate1[-length(estimate1)]
  if(length(which(diff(estimate1)<5))>0) estimate1=estimate1[-which(diff(estimate1)<5)]
  num[1,j]=length(estimate1)
  TP[1,j]=TruePositive(estimate1,tau)
  FP[1,j]=FalsePositive(estimate1,tau)
  ##############################method2 SaRa
  sara=MultiScan(y,h=20)
  thres=threshold(y,alpha=0.05,h=20)
  estimate2=sara$index[which(sara$S.index.>thres)]
  num[2,j]=length(estimate2)
  TP[2,j]=TruePositive(estimate2,tau)
  FP[2,j]=FalsePositive(estimate2,tau)
  ##############################method3 gSCAD+adjusted  
  y0=as.vector(t((diag(1,n)-S)%*%t(y)))            
  gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=40) 
  opt=which.min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #minimize gBIC
  estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[d*(1:(length(which(gSCAD$beta[,opt]!=0))/d))]/d-1)[-1])
  if(length(which(diff(estimate3)<5))>0)  estimate3=estimate3[-which(diff(estimate3)<5)]
  num[3,j]=length(estimate3)
  TP[3,j]=TruePositive(estimate3,tau)
  FP[3,j]=FalsePositive(estimate3,tau)
  #############################wave estimation
  WAVE=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    jumpsize=screening(x,y[k,])[estimate3]
    fit=0
    for(i in 1:length(estimate3)) {fit<-fit+jumpsize[i]*(x>estimate3[i])}
    WAVE[k,]=S%*%(y[k,]-fit)
  }
  ##############################method4 CBS-adjusted
  yy=y-WAVE
  CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
  estimate4=sort(unique(CBS$output[,4]))
  estimate4=estimate4[-length(estimate4)]
  if(length(which(diff(estimate4)<5))>0) estimate4=estimate4[-which(diff(estimate4)<5)]
  num[4,j]=length(estimate4)
  TP[4,j]=TruePositive(estimate4,tau)
  FP[4,j]=FalsePositive(estimate4,tau)
  ##############################method5 SaRa-adjusted
  sara=MultiScan(yy,h=20)
  thres=threshold(yy,alpha=0.05,h=20)
  estimate5=sara$index[which(sara$S.index.>thres)]
  num[5,j]=length(estimate5)
  TP[5,j]=TruePositive(estimate5,tau)
  FP[5,j]=FalsePositive(estimate5,tau)
}

set.seed(23)
for(j in 1:100){
  ##############################data generation for Scenario II 
  beta=matrix(sample(c(-1,-0.5,0,0.5,1), size=J*d, replace = TRUE), nrow=d)
  signal=matrix(rep(0, n*d), nrow=d)
  y=matrix(rep(0, n*d), nrow=d)
  wave=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    wave[k, ]=sample(c(-1,-1.5,1,1.5),size=1)*rwave         #realwave
    signal[k, ]=0
    for(i in 1:J) {signal[k, ]<-signal[k, ]+beta[k,i]*(x>tau[i])}
    y[k, ]=signal[k, ]+wave[k, ]+rnorm(n,mean=0,sd=0.2)  
  }
  ##############################method1 CBS
  CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
  estimate1=sort(unique(CBS$output[,4]))
  estimate1=estimate1[-length(estimate1)]
  if(length(which(diff(estimate1)<5))>0) estimate1=estimate1[-which(diff(estimate1)<5)]
  num1[1,j]=length(estimate1)
  TP1[1,j]=TruePositive(estimate1,tau)
  FP1[1,j]=FalsePositive(estimate1,tau)
  ##############################method2 SaRa
  sara=MultiScan(y,h=20)
  thres=threshold(y,alpha=0.05,h=20)
  estimate2=sara$index[which(sara$S.index.>thres)]
  num1[2,j]=length(estimate2)
  TP1[2,j]=TruePositive(estimate2,tau)
  FP1[2,j]=FalsePositive(estimate2,tau)
  ##############################method3 gSCAD+adjusted  
  y0=as.vector(t((diag(1,n)-S)%*%t(y)))            
  gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=40) 
  opt=which.min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #minimize gBIC
  estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[d*(1:(length(which(gSCAD$beta[,opt]!=0))/d))]/d-1)[-1])
  if(length(which(diff(estimate3)<5))>0)  estimate3=estimate3[-which(diff(estimate3)<5)]
  num1[3,j]=length(estimate3)
  TP1[3,j]=TruePositive(estimate3,tau)
  FP1[3,j]=FalsePositive(estimate3,tau)
  #############################wave estimation
  WAVE=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    jumpsize=screening(x,y[k,])[estimate3]
    fit=0
    for(i in 1:length(estimate3)) {fit<-fit+jumpsize[i]*(x>estimate3[i])}
    WAVE[k,]=S%*%(y[k,]-fit)
  }
  ##############################method4 CBS-adjusted
  yy=y-WAVE
  CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
  estimate4=sort(unique(CBS$output[,4]))
  estimate4=estimate4[-length(estimate4)]
  if(length(which(diff(estimate4)<5))>0) estimate4=estimate4[-which(diff(estimate4)<5)]
  num1[4,j]=length(estimate4)
  TP1[4,j]=TruePositive(estimate4,tau)
  FP1[4,j]=FalsePositive(estimate4,tau)
  ##############################method5 SaRa-adjusted
  sara=MultiScan(yy,h=20)
  thres=threshold(yy,alpha=0.05,h=20)
  estimate5=sara$index[which(sara$S.index.>thres)]
  num1[5,j]=length(estimate5)
  TP1[5,j]=TruePositive(estimate5,tau)
  FP1[5,j]=FalsePositive(estimate5,tau)
}

##############################show table
rbind(apply(num,1,mean),
      apply(num,1,sd),
      apply(TP,1,mean),
      apply(TP,1,sd),
      apply(FP,1,mean),
      apply(FP,1,sd),
      apply(num1,1,mean),
      apply(num1,1,sd),
      apply(TP1,1,mean),
      apply(TP1,1,sd),
      apply(FP1,1,mean),
      apply(FP1,1,sd))
