library(grpreg)
library(KernSmooth)
library(np)

#################################define some functions
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

#################################Simulation1
set.seed(23)
n=1000       #sample size
d=5            #number of sequences
J=10           #number of change points
x=1:n         #location
tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
h=10
W=function(t,x,h){diag(1/h*3/4*(1-((x-t)/h)^2)*((1-((x-t)/h)^2)>0))}
D=function(t,x,h){matrix(c(rep(1,times=length(x)),(x-t)/h),nrow=length(x),ncol=2)}
s=function(t,x,h){t(c(1,0))%*%solve(t(D(t,x,h))%*%W(t,x,h)%*%D(t,x,h))%*%t(D(t,x,h))%*%W(t,x,h)}
S=NULL                    #smoothing matrix S
for(i in 1:length(x)){S=rbind(S,s(x[i],x,h))}
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) x1[i,(i+1):n]=0
X=kronecker(x1,diag(1, d))       
XX=kronecker(cbind(x1[,1],(diag(1,n)-S)%*%x1[ ,-1]), diag(1, d))                   
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,d))  
TP1=matrix(NA,nr=100,nc=20)                #true positive
TP2=matrix(NA,nr=100,nc=20)
TP3=matrix(NA,nr=100,nc=20)
TP4=matrix(NA,nr=100,nc=20)
FP1=matrix(NA,nr=100,nc=20)                #false positive
FP2=matrix(NA,nr=100,nc=20)
FP3=matrix(NA,nr=100,nc=20)
FP4=matrix(NA,nr=100,nc=20)
for(j in 1:100){
  ##############################data generation
  beta=matrix(sample(c(-1,-0.5,0,0.5,1), size=J*d, replace = TRUE), nrow=d)
  theta=runif(1,min=0,max=2*pi)
  signal=matrix(rep(0, n*d), nrow=d)
  y=matrix(rep(0, n*d), nrow=d)    #with wave
  y0=matrix(rep(0, n*d), nrow=d)    #without wave
  wave=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    wave[k, ]=sample(c(-0.5,-0.3,0.3,0.5),size=1)*sin(pi*x/50+theta)
    signal[k, ]=0
    for(i in 1:J) {signal[k, ]<-signal[k, ]+beta[k,i]*(x>tau[i])}
    y0[k, ]=signal[k, ]+rnorm(n,mean=0,sd=0.2)  
    y[k, ]=wave[k, ]+y0[k, ]
  }
  ##############################gLasso+no wave
  yy=as.vector(y0)                 
  gLasso=grpreg(X=X[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grLasso",family="gaussian", gmax=20)   
  for(i in 1:dim(gLasso$beta)[2]){
    estimate=as.vector((which(gLasso$beta[,i]!=0)[d*(1:(length(which(gLasso$beta[,i]!=0))/d))]/d-1)[-1])
    if(length(estimate)<=20)  TP1[j,length(estimate)]=TruePositive(estimate, tau)
    if(length(estimate)<=20)  FP1[j,length(estimate)]=FalsePositive(estimate, tau)
  }
  ############################gSCAD+no wave
  gSCAD=grpreg(X=X[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grSCAD",family="gaussian", gmax=20) 
  for(i in 1:dim(gSCAD$beta)[2]){
    estimate=as.vector((which(gSCAD$beta[,i]!=0)[d*(1:(length(which(gSCAD$beta[,i]!=0))/d))]/d-1)[-1])
    if(length(estimate)<=20)  TP2[j,length(estimate)]=TruePositive(estimate, tau)
    if(length(estimate)<=20)  FP2[j,length(estimate)]=FalsePositive(estimate, tau)
  }
  ##############################gSCAD+wave   
  yy=as.vector(y)          
  gSCAD=grpreg(X=X[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grSCAD",family="gaussian", gmax=20) 
  for(i in 1:dim(gSCAD$beta)[2]){
    estimate=as.vector((which(gSCAD$beta[,i]!=0)[d*(1:(length(which(gSCAD$beta[,i]!=0))/d))]/d-1)[-1])
    if(length(estimate)<=20)  TP3[j,length(estimate)]=TruePositive(estimate, tau)
    if(length(estimate)<=20)  FP3[j,length(estimate)]=FalsePositive(estimate, tau)
  }
  ##############################gSCAD+adjusted
  yy=as.vector(t((diag(1,n)-S)%*%t(y)))
  gSCAD=grpreg(X=XX[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grSCAD",family="gaussian", gmax=20) 
  for(i in 1:dim(gSCAD$beta)[2]){
    estimate=as.vector((which(gSCAD$beta[,i]!=0)[d*(1:(length(which(gSCAD$beta[,i]!=0))/d))]/d-1)[-1])
    if(length(estimate)<=20)  TP4[j,length(estimate)]=TruePositive(estimate, tau)
    if(length(estimate)<=20)  FP4[j,length(estimate)]=FalsePositive(estimate, tau)
  }
}

##############################show plots
par(mfrow=c(2,2))
plot(x=1:10,y=apply(TP1,2,mean,na.rm=T)[1:10],xlab="Number of detected change-points",ylab="True positives",ylim=c(-1,10),main="Scenario I, True positives", type="b")
lines(x=1:10,y=apply(TP2,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(FP1,2,mean,na.rm=T)[1:10], xlab="Number of detected change-points",ylab="False positives",ylim=c(-1,10),main="Scenario I, False positives", type="b")
lines(x=1:10,y=apply(FP2,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("topright", legend=c("SCAD","Lasso"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(TP3,2,mean,na.rm=T)[1:10], xlab="Number of detected change-points",ylab="True positives",ylim=c(-1,10),main="Scenario II, True positives", type="b")
lines(x=1:10,y=apply(TP4,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("bottomright", legend=c("adjustment","no adjustment"),pch=c(4,1),bty="n")
plot(x=1:10,y=apply(FP3,2,mean,na.rm=T)[1:10], xlab="Number of detected change-points",ylab="False positives",ylim=c(-1,10),main="Scenario II, False positives", type="b")
lines(x=1:10,y=apply(FP4,2,mean,na.rm=T)[1:10],type="b",pch=4)
legend("topright", legend=c("adjustment","no adjustment"),pch=c(4,1),bty="n")

estimate=as.vector((which(gSCAD$beta[,J]!=0)[d*(1:(length(which(gSCAD$beta[,J]!=0))/d))]/d-1)[-1])
jumpsize=screening(x,y[1,])[estimate]
fit=0
for(i in 1:length(estimate)) {fit<-fit+jumpsize[i]*(x>estimate[i])}
Wave1=S%*%(y[1,]-fit)
plot(wave[1,]~x,type="l",ylim=c(-1,1),xlab="location",ylab="f1(x)",main="true wave")
plot(Wave1~x,type="l",ylim=c(-1,1),xlab="location",ylab="f1(x)",main="estimated wave")
jumpsize=screening(x,y[2,])[estimate]
fit=0
for(i in 1:length(estimate)) {fit<-fit+jumpsize[i]*(x>estimate[i])}
Wave2=S%*%(y[2,]-fit)
plot(wave[2,]~x,type="l",ylim=c(-1,1),xlab="location",ylab="f2(x)",main="true wave")
plot(Wave2~x,type="l",ylim=c(-1,1),xlab="location",ylab="f2(x)",main="estimated wave")