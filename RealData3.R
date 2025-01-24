library(seqbbs)   #is available at https://github.com/metalhelix/seqbbs
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

#################################real data 3
test_filename <- system.file("extdata", "paper.txt", package="seqbbs")
ratios <- read.table(test_filename, header = FALSE)
y=ratios$V1               
n=length(y)     
x=1:n
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) {x1[i,(i+1):n]=0}                   

##############################method1 CBS
set.seed(18)
CBS=DNAcopy::segment(CNA(as.matrix(y), chrom=rep(1,n), maploc=1:n))
estimate1=sort(unique(CBS$output[,4]))
estimate1=estimate1[-length(estimate1)]
if(length(which(diff(estimate1)<5))>0) estimate1=estimate1[-which(diff(estimate1)<5)]
##############################method2 SaRa
D=abs(localDiagnostic(y, h=10))
index=localMax(D,span=10)
lambda=2*sqrt(2/10)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
candidate=index[which(D[index]>lambda)]
candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
bic=rep(0,length(candidate1))
for(j in 1:length(candidate1)){
  model=lm(y~x1[,candidate1[1:j]])
  bic[j]=n*log(summary(model)$sigma)+j*log(n)
}
j1=which.min(bic)
estimate2=sort(candidate1[1:j1])
##############################method3 gSCAD+adjusted  
group=NULL
for (i in 1:n) group<-c(group,rep(i-1,1))  
hh=c(10,15,20)
gBIC=1:3
for(j in 1:3){
  h=hh[j]
  xx1=x1
  for(i in 2:n) {xx1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
  y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y         
  gSCAD=grpreg(X=xx1[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=50) 
  gBIC[j]=min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #minimize gBIC for fixed h
}
h=hh[which.min(gBIC)]        #optimal h
xx1=x1
for(i in 2:n) {xx1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y           
gSCAD=grpreg(X=xx1[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grSCAD",family="gaussian", gmax=100) 
opt=which.min(log(gSCAD$deviance/n)+gSCAD$df*0.2*log(n)*log(log(n))/n)   #optimal lambda
estimate3=as.vector((which(gSCAD$beta[,opt]!=0)[-1]))-1
if(length(which(diff(estimate3)<5))>0)  estimate3=estimate3[-(which(diff(estimate3)<5)+1)]
#############################wave estimation
fit=rep(0, n)
jumpsize=screening(x,y)[estimate3]
for(i in 1:length(estimate3)) {fit<-fit+jumpsize[i]*(x>estimate3[i])}
WAVE=locpoly(x,y-fit,kernel="epanech",bandwidth=h,gridsize=n)$y
##############################method4 CBS-adjusted
yy=y-WAVE
CBS=DNAcopy::segment(CNA(as.matrix(yy), chrom=rep(1,n), maploc=1:n))
estimate4=sort(unique(CBS$output[,4]))
estimate4=estimate4[-length(estimate4)]
if(length(which(diff(estimate4)<5))>0) estimate4=estimate4[-which(diff(estimate4)<5)]
##############################method5 SaRa-adjusted
D=abs(localDiagnostic(yy, h=10))
index=localMax(D,span=10)
lambda=2*sqrt(2/10)*mad(diff(yy))/sqrt(2)        #use robust estimator of sigma
candidate=index[which(D[index]>lambda)]
candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
bic=rep(0,length(candidate1))
for(j in 1:length(candidate1)){
  model=lm(yy~x1[,candidate1[1:j]])
  bic[j]=n*log(summary(model)$sigma)+j*log(n)
}
j1=which.min(bic)
estimate5=sort(candidate1[1:j1])

##############################show figures
par(mfrow=c(2,3))
plot(x, y, xlim=c(300,800), ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="original CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate1)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate1,n)[i]+1):(c(0,estimate1,n)[i+1])]), c(0,estimate1,n)[i+1]-c(0,estimate1,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="original SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate2)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate2,n)[i]+1):(c(0,estimate2,n)[i+1])]), c(0,estimate2,n)[i+1]-c(0,estimate2,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="proposed", pch=20, col=8)
lines(x ,fit+WAVE, col=2, lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="adjusted CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate4)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate4,n)[i]+1):(c(0,estimate4,n)[i+1])]), c(0,estimate4,n)[i+1]-c(0,estimate4,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="adjusted SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate5)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate5,n)[i]+1):(c(0,estimate5,n)[i+1])]), c(0,estimate5,n)[i+1]-c(0,estimate5,n)[i]))
}
lines(fitted,col=2,lwd=2)
