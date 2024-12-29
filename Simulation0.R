#######get data, please modify the file path if necessary
CGHdata <- read.csv("C:/Users/PC/Desktop/我的文档/Research/Projects/change points(wave adjust)/CGHdataset.csv")

library(ncvreg)
library(grpreg)
library(KernSmooth)
library(imputeTS)

######estimate wave from real data
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
wave=10*locpoly(x,re,kernel="epanech",bandwidth=20,gridsize=500)$y
######simulation in introduction
set.seed(23)
n=500
x1=matrix(1,nrow=n,ncol=n)     #matrix X
for(i in 1:(n-1))
  x1[i,(i+1):n]=0
x=1:n
tau=c(140,185,315,355)
beta=c(1,-1,-1,1)
sigma=0.2
signal=0
J=4
for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
estimate1=NULL
estimate2=NULL
for(i in 1:100){
  e=rnorm(n,mean=0,sd=sigma) 
  ######no wave + lasso
  y=signal+e  
  cvfit<-cv.ncvreg(x1[ ,2:n], y, penalty="lasso", max.iter=100000, eps=1e-3)
  estimate=which(abs(coef(cvfit))>sigma)
  estimate1=c(estimate1, as.vector(estimate))
  ######wave + lasso
  y=signal+wave+e 
  cvfit<-cv.ncvreg(x1[ ,2:n], y, penalty="lasso", max.iter=100000, eps=1e-3)
  estimate=which(abs(coef(cvfit))>sigma)
  estimate2=c(estimate2, as.vector(estimate))
}
######show Figure1
par(mfrow=c(2,2))
plot(x,signal, type="l", lwd=2, ylim=c(-1.5,1.5), xlab="locations", ylab="signal", main="true signal")
abline(v=c(140,185,315,355), lty=2)
plot(x,wave, type="l", lwd=2, ylim=c(-1.5,1.5), xlab="locations", ylab="wave", main="wave patterns")
hist(estimate1[-c(which(estimate1<5),which(estimate1>(n-5)))],breaks=50,xlim=c(0,n),ylim=c(0,100),xlab="locations",main="Lasso estimates without waves")
hist(estimate2[-c(which(estimate2<5),which(estimate2>(n-5)))],breaks=50,xlim=c(0,n),ylim=c(0,100),xlab="locations",main="Lasso estimates with waves")
