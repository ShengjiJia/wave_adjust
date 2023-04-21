library(grpreg)
library(KernSmooth)
library(np)
library(DNAcopy)
library(imputeTS)

#######get data
SNPdata <- read.delim("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(wave adjust)/SNPdata.txt")

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

#################################chromosome 21
y=rbind(SNPdata$X99HI0697A.Log.R.Ratio[which(SNPdata$Chr==21)],
        SNPdata$X99HI0698C.Log.R.Ratio[which(SNPdata$Chr==21)],
        SNPdata$X99HI0700A.Log.R.Ratio[which(SNPdata$Chr==21)])              
n=ncol(y)
d=nrow(y)
x=1:n
y[which(abs(y)>1.5)]=NA     #delete outliers
for (i in 1:d) {y[i,]=na_ma(y[i,],k=5,weighting="linear")}        #imputation
h=10
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) {x1[i,(i+1):n]=0} 
for(i in 2:n) {x1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
XX=kronecker(x1, diag(1, d))                   
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,d))  

##############################method1 CBS
set.seed(18)
CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
estimate1=sort(unique(CBS$output[,4]))
estimate1=estimate1[-length(estimate1)]
if(length(which(diff(estimate1)<5))>0)  estimate1=estimate1[-which(diff(estimate1)<5)]
##############################method2 SaRa
sara=MultiScan(y,h=10)
thres=threshold(y,alpha=0.1,h=10)
estimate2=sara$index[which(sara$S.index.>thres)]
##############################method3 gSCAD+adjusted  
y0=matrix(0,nrow=d,ncol=n)
for (i in 1:d) {y0[i,]=y[i,]-locpoly(x,y[i,],kernel="epanech",bandwidth=h,gridsize=n)$y}
y0=as.vector(y0)            
gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", lambda=0.003*0.9^(1:10-1)) 
opt=which.min(log(gSCAD$loss/n)+gSCAD$df*0.3*log(n)*log(log(n))/n)   #minimize gBIC
estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[d*(1:(length(which(gSCAD$beta[,opt]!=0))/d))]/d-1)[-1])
if(length(which(diff(estimate3)<5))>0)  estimate3=estimate3[-which(diff(estimate3)<5)]
#############################wave estimation
wave=matrix(rep(0, n*d), nrow=d)
fit=matrix(rep(0, n*d), nrow=d)
for(k in 1:d){
  jumpsize=screening(x,y[k,])[estimate3]
  for(i in 1:length(estimate3)) {fit[k,]<-fit[k,]+jumpsize[i]*(x>estimate3[i])}
  wave[k,]=locpoly(x,y[k,]-fit[k,],kernel="epanech",bandwidth=h,gridsize=n)$y
}
##############################method4 CBS-adjusted
yy=y-wave
CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
estimate4=sort(unique(CBS$output[,4]))
estimate4=estimate4[-length(estimate4)]
if(length(which(diff(estimate4)<5))>0)  estimate4=estimate4[-which(diff(estimate4)<5)]
##############################method5 SaRa-adjusted
sara=MultiScan(yy,h=10)
thres=threshold(yy,alpha=0.1,h=10)
estimate5=sara$index[which(sara$S.index.>thres)]

##############################show plots
par(mfrow=c(2,3))
yy=y[2,]
plot(x, yy, xlim=c(1000,1500), ylim=c(-1,0.5), xlab="locations", ylab="Log R ratio", main="original CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate1)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate1,n)[i]+1):(c(0,estimate1,n)[i+1])]), c(0,estimate1,n)[i+1]-c(0,estimate1,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(1000,1500),  ylim=c(-1,0.5), xlab="locations", ylab="Log R ratio", main="original SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate2)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate2,n)[i]+1):(c(0,estimate2,n)[i+1])]), c(0,estimate2,n)[i+1]-c(0,estimate2,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(1000,1500),  ylim=c(-1,0.5), xlab="locations", ylab="Log R ratio", main="proposed", pch=20, col=8)
lines(x ,fit[2,]+wave[2,], col=2, lwd=2)
plot(x, yy, xlim=c(1000,1500),  ylim=c(-1,0.5), xlab="locations", ylab="Log R ratio", main="adjusted CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate4)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate4,n)[i]+1):(c(0,estimate4,n)[i+1])]), c(0,estimate4,n)[i+1]-c(0,estimate4,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, yy, xlim=c(1000,1500),  ylim=c(-1,0.5), xlab="locations", ylab="Log R ratio", main="adjusted SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate5)+1)){
  fitted=c(fitted,rep(mean(yy[(c(0,estimate5,n)[i]+1):(c(0,estimate5,n)[i+1])]), c(0,estimate5,n)[i+1]-c(0,estimate5,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x=(120:250)/10,y=0.5*sin(2*pi*(120:250)/250),xlab="locations",ylab="signals",ylim=c(-1,2),type="l",lty=2,lwd=2,col="blue")
lines(x=(120:250)/10,y=-0.2*(((120:250)/10)<14)+0.3*(((120:250)/10)>=14)*(((120:250)/10)<22)+0.1*(((120:250)/10)>=22),type="l",lty=2,lwd=2,col="red")
lines(x=(120:250)/10,y=0.5*sin(2*pi*(120:250)/250)-0.2*(((120:250)/10)<14)+0.3*(((120:250)/10)>=14)*(((120:250)/10)<22)+0.1*(((120:250)/10)>=22),lwd=2,col=1)
legend("topright", legend=c("waves","CNVs","observations"),lty=c(2,2,1),lwd=c(2,2,2),col=c("blue","red","black"),bty="n")

#################################chromosome 22
y=rbind(SNPdata$X99HI0697A.Log.R.Ratio[which(SNPdata$Chr==22)],
        SNPdata$X99HI0698C.Log.R.Ratio[which(SNPdata$Chr==22)],
        SNPdata$X99HI0700A.Log.R.Ratio[which(SNPdata$Chr==22)])              
n=ncol(y)
d=nrow(y)
x=1:n
y[which(abs(y)>1.5)]=NA     #delete outliers
for (i in 1:d) {y[i,]=na_ma(y[i,],k=5,weighting="linear")}        #imputation
h=10
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) {x1[i,(i+1):n]=0} 
for(i in 2:n) {x1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
XX=kronecker(x1, diag(1, d))                   
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,d))  

##############################method1 CBS
set.seed(18)
CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
Estimate1=sort(unique(CBS$output[,4]))
Estimate1=Estimate1[-length(Estimate1)]
if(length(which(diff(Estimate1)<5))>0)  Estimate1=Estimate1[-which(diff(Estimate1)<5)]
##############################method2 SaRa
sara=MultiScan(y,h=10)
thres=threshold(y,alpha=0.1,h=10)
Estimate2=sara$index[which(sara$S.index.>thres)]
##############################method3 gSCAD+adjusted  
y0=matrix(0,nrow=d,ncol=n)
for (i in 1:d) {y0[i,]=y[i,]-locpoly(x,y[i,],kernel="epanech",bandwidth=h,gridsize=n)$y}
y0=as.vector(y0)            
gSCAD=grpreg(X=XX[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", lambda=0.003*0.95^(1:10-1)) 
opt=which.min(log(gSCAD$loss/n)+gSCAD$df*0.3*log(n)*log(log(n))/n)   #minimize gBIC
Estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[d*(1:(length(which(gSCAD$beta[,opt]!=0))/d))]/d-1)[-1])
if(length(which(diff(Estimate3)<5))>0)  Estimate3=Estimate3[-which(diff(Estimate3)<5)]
#############################wave estimation
WAVE=matrix(rep(0, n*d), nrow=d)
FIT=matrix(rep(0, n*d), nrow=d)
for(k in 1:d){
  jumpsize=screening(x,y[k,])[Estimate3]
  for(i in 1:length(Estimate3)) {FIT[k,]<-FIT[k,]+jumpsize[i]*(x>Estimate3[i])}
  WAVE[k,]=locpoly(x,y[k,]-FIT[k,],kernel="epanech",bandwidth=h,gridsize=n)$y
}
##############################method4 CBS-adjusted
yy=y-WAVE
CBS=DNAcopy::segment(CNA(t(yy), chrom=rep(1,n), maploc=1:n))
Estimate4=sort(unique(CBS$output[,4]))
Estimate4=Estimate4[-length(Estimate4)]
if(length(which(diff(Estimate4)<5))>0)  Estimate4=Estimate4[-which(diff(Estimate4)<5)]
##############################method5 SaRa-adjusted
sara=MultiScan(yy,h=10)
thres=threshold(yy,alpha=0.1,h=10)
Estimate5=sara$index[which(sara$S.index.>thres)]

###############################show table
rbind(c(length(estimate1),length(estimate2),length(estimate3),length(estimate4),length(estimate5)),
      c(length(Estimate1),length(Estimate2),length(Estimate3),length(Estimate4),length(Estimate5)))
