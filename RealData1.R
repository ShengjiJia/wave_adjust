#######get data, please modify the file path if necessary
CGHdata <- read.csv("C:/Users/PC/Desktop/我的文档/Research/Projects/change points(wave adjust)/CGHdataset.csv")

library(grpreg)
library(KernSmooth)
library(DNAcopy)
library(imputeTS)
#################################define some functions
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

#################################real data 1
#index=c(11,19,45,56)
index=c(6,19,20,21)
data=CGHdata[2:2301,(1+3*index)]               
n=nrow(data)
d=ncol(data)
x=1:n
y=t(as.matrix(data))
y[which(abs(y)>2)]=NA     #delete outliers
for (i in 1:d) {y[i,]=na_ma(y[i,],k=5,weighting="linear")}        #imputation
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) {x1[i,(i+1):n]=0} 
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,d))                    

##############################method1 CBS
set.seed(18)
CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
estimate1=sort(unique(CBS$output[,4]))
estimate1=estimate1[-length(estimate1)]
if(length(which(diff(estimate1)<5))>0) estimate1=estimate1[-which(diff(estimate1)<5)]
##############################method2 SaRa
sara=MultiScan(y,h=20)
thres=threshold(y,alpha=0.05,h=20)
estimate2=sara$index[which(sara$S.index.>thres)]
##############################method3 gSCAD+adjusted  
hh=c(10,15,20)
gBIC=1:3
for(j in 1:3){
  h=hh[j]
  xx1=x1
  for(i in 2:n) {xx1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
  XX=kronecker(xx1, diag(1, d)) 
  y0=matrix(0,nrow=d,ncol=n)
  for (i in 1:d) {y0[i,]=y[i,]-locpoly(x,y[i,],kernel="epanech",bandwidth=h,gridsize=n)$y}
  y0=as.vector(y0)            
  gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=40) 
  gBIC[j]=min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #minimize gBIC for fixed h
}
h=hh[which.min(gBIC)]        #optimal h
for(i in 2:n) {x1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
XX=kronecker(x1, diag(1, d)) 
y0=matrix(0,nrow=d,ncol=n)
for (i in 1:d) {y0[i,]=y[i,]-locpoly(x,y[i,],kernel="epanech",bandwidth=h,gridsize=n)$y}
y0=as.vector(y0)            
gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=40) 
opt=which.min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #optimal lambda
estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[d*(1:(length(which(gSCAD$beta[,opt]!=0))/d))]/d-1)[-1])
if(length(which(diff(estimate3)<5))>0)  estimate3=estimate3[-which(diff(estimate3)<5)]
#############################wave estimation
WAVE=matrix(rep(0, n*d), nrow=d)
fit=matrix(rep(0, n*d), nrow=d)
for(k in 1:d){
  jumpsize=screening(x,y[k,])[estimate3]
  for(i in 1:length(estimate3)) {fit[k,]<-fit[k,]+jumpsize[i]*(x>estimate3[i])}
  WAVE[k,]=locpoly(x,y[k,]-fit[k,],kernel="epanech",bandwidth=h,gridsize=n)$y
}
##############################method4 CBS-adjusted
yy=y-WAVE
CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
estimate4=sort(unique(CBS$output[,4]))
estimate4=estimate4[-length(estimate4)]
if(length(which(diff(estimate4)<5))>0) estimate4=estimate4[-which(diff(estimate4)<5)]
##############################method5 SaRa-adjusted
sara=MultiScan(yy,h=20)
thres=threshold(yy,alpha=0.05,h=20)
estimate5=sara$index[which(sara$S.index.>thres)]

##############################show figures
par(mfrow=c(4,1))
plot(x, y[1,], ylim=c(-2,2), xlab="locations", ylab="Log 2 ratio", main="X1333-4", pch=20, col=8)
lines(x ,fit[1,]+WAVE[1,], col=2, lwd=2) 
abline(v=estimate3,lty=2)
plot(x, y[2,], xlab="locations", ylab="Log 2 ratio", main="X1533-1", pch=20, col=8)
lines(x ,fit[2,]+WAVE[2,], col=2, lwd=2) 
abline(v=estimate3,lty=2)
plot(x, y[3,], xlab="locations", ylab="Log 2 ratio", main="X1533-10", pch=20, col=8)
lines(x ,fit[3,]+WAVE[3,], col=2, lwd=2) 
abline(v=estimate3,lty=2)
plot(x, y[4,], xlab="locations", ylab="Log 2 ratio", main="X1533-13", pch=20, col=8)
lines(x ,fit[4,]+WAVE[4,], col=2, lwd=2)
abline(v=estimate3,lty=2)

par(mfrow=c(3,2))
yy=y[4,]
plot(x, yy, xlim=c(300,800), ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="original CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate1)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate1,n)[i]+1):(c(0,estimate1,n)[i+1])]), c(0,estimate1,n)[i+1]-c(0,estimate1,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(300,800),  ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="original SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate2)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate2,n)[i]+1):(c(0,estimate2,n)[i+1])]), c(0,estimate2,n)[i+1]-c(0,estimate2,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(300,800),  ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="proposed", pch=20, col=8)
lines(x ,fit[4,]+WAVE[4,], col=2, lwd=2)
plot(x, yy, xlim=c(300,800),  ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="adjusted CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate4)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate4,n)[i]+1):(c(0,estimate4,n)[i+1])]), c(0,estimate4,n)[i+1]-c(0,estimate4,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(300,800),  ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="adjusted SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate5)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate5,n)[i]+1):(c(0,estimate5,n)[i+1])]), c(0,estimate5,n)[i+1]-c(0,estimate5,n)[i]))
}
lines(fitted,col=2,lwd=2)
re=yy-fit[4,]-WAVE[4,]
qqnorm(re)
