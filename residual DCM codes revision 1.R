###################################################################
#                                                                 #
#  The codes for the manuscript "Surrogate Residuals for DCMs"    #
#                                                                 #
#                        Chao Cheng                               #
#                                                                 #
#                      cqplus@126.com                             #
#                                                                 #
###################################################################
# Contents:                                                       #
# Section 1: Functions                                            #
# Section 2: Numerical Examples                                   #
###################################################################

#### initialize

library("truncnorm");library("nnet");library("tmvtnorm");library("MNP")
library("miscTools");library("foreign");
library("mlogit")

##################################################################
#
# Section 1: FUNCTIONS
#
##################################################################

##############################################################################
#Surrogate residuals for discrete choice models
##############################################################################
#Arguments:
# object:  An object of class glm, mlogit, or mnp
# means : the mean structures (i.e the deterministic parts).
# ...   : Additional optional arguments. (Currently ignored.)
##############################################################################
#Output : a vectir or matrix representing the surrogate residuals
##############################################################################
#Note:
# This function is only for logit regression (based on glm function), probit 
# regresssion (based on glm function), multinominal logistic regression (based
# on mlogit package) and multinominal probit regression (based on MNP package)
# now.
##############################################################################
surrogate_resid=function(object,means,...) {UseMethod("surrogate_resid")}
surrogate_resid.default=function(object,means=NULL,...) {
  if (class(object)[1]=="glm") {
    if (is.null(means)) {
      stop("please input the data representing the mean structure (determinstic part)!")
    }
    if (object$family$link=="logit") {
      n=length(means)
      min.v=rep(-Inf,n);max.v=rep(Inf,n)
      min.v[which(object$y>0.5)]=0
      max.v[which(object$y<0.5)]=0
      out=rlogit(n=n,location=means,min=min.v,max=max.v)-means
      class(out)="surrogate_resid"
      out
    } else if (object$family$link=="probit"){
      n=length(means)
      min.v=rep(-Inf,n);max.v=rep(Inf,n)
      min.v[which(object$y>0.5)]=0
      max.v[which(object$y<0.5)]=0
      out=rprobit(n=n,location=means,min=min.v,max=max.v)-means
      class(out)="surrogate_resid"
      out
    }
  } else if (class(object)=="mlogit") {
    fitted =fitted(object,outcome=F)
    M=dim(fitted)[2];
    mean.mat = log(fitted) - log(fitted)[,1] %*% t(rep(1,M))
    N=dim(mean.mat)[1]
    outcome=object$probabilities+object$residuals
    z=rep(0,N)
    for (i in (1:N)) {
      z[i] = which(outcome[i,]>0.999999)
    }
    error=matrix(0,nrow=dim(mean.mat)[1],ncol=M)
    for (i in (1:dim(mean.mat)[1])) {
      error1 = rmgumbel(n=1,coef=mean.mat[i,],max.loc=z[i],burn.in=100)
      error1[is.infinite(error1)]=0
      error[i,] = error1 - 0.5772
    }
    colnames(error) = colnames(object$probabilities)
    class(error)="surrogate_resid"
    error
  } else if (class(object)=="mnp") {
    if (is.null(means)) {
      stop("please input the data representing the mean structure (determinstic part)!")
    }
    cov=covmnp(object)
    N=length(object$y)
    z=object$y
    error=matrix(0,nrow=N,ncol=length(object$alt)-1)
    for (i in (1:N)) {
      error[i,] = rmulti_trun(1,Mean=means[i,],Var=cov,
                              max.location=z[i])
    }
    error=error-means
    colnames(error)=object$alt[-which(object$alt==object$base)]
    class(error)="surrogate_resid"
    error
  } else {
    stop("The model is not currently available in this package.")
  }
}


##############################################################################
# Other function supporting surrogate_resid() are shown below
##############################################################################
rlogit=function(n=1,location=0,min=-Inf,max=Inf) {
  #n=1;location=0;min=-Inf;max=Inf
  u=runif(n,0,1)
  stats::qlogis(u*(stats::plogis(q=max,location=location)-stats::plogis(q=min,location=location))+stats::plogis(q=min,location=location),
                location=location)
}

rprobit=function(n=1,location=0,min=-Inf,max=Inf) {
  u=runif(n,0,1)
  stats::qnorm(u*(stats::pnorm(q=max,mean=location)-stats::pnorm(q=min,mean=location))+stats::pnorm(q=min,mean=location),
               mean=location)
}
rmulti_trun=function(n,Mean=rep(0,3),Var=diag(rep(1,3)),max.location=1) {
  #n=1000;Mean=rep(0,5);Var=rep(1,5);max.location=3
  n.dim=length(Mean)
  D = matrix(0,nrow=n.dim,ncol=n.dim)
  D[,max.location] = rep(1,n.dim)
  for (i in (1:n.dim)) {
    if (max.location != i) {
      D[i,i] = -1
    }
  }
  upper = rep(Inf, n.dim)
  lower = rep(-Inf, n.dim)  
  for (i in (1:n.dim)) {
    if (max.location != i) {
      lower[i]=0
    } else {
      lower[i]=0
    }
  }
  if (max.location==0) {
    D = diag(rep(1,n.dim))
    upper = rep(0, n.dim)
    lower = rep(-Inf, n.dim)
  }
  z = tmvtnorm::rtmvnorm(n=n,lower=lower,upper=upper,mean=Mean,sigma=Var,D=D,algorithm="gibbsR",burn.in.samples=100)
  return(z)
}
coemnp=function(x) {
  coef=coef.mnp(x)
  coef=apply(coef,2,mean)
  coef
}
covmnp = function(x) {
  Cov=cov.mnp(x)
  Cov=apply(Cov,c(1,2),mean)
  Cov
}

## The distribution function of standard gumbel distribution
pgumbel=function(x) { exp(-exp(-x)) }
qgumbel = function(x) {-log(-log(x))}
fgumbel=function(x){exp(-x-exp(-x))}
## Generate Standard Gumbel Distribution
rgumbel = function(n,a=-Inf,b=Inf) {
  if (a==-Inf && b==Inf) {
    x = runif(n,0,1)
  } else if(a==-Inf && is.infinite(b)==F) {
    x = runif(n,0,pgumbel(b))
  } else if (is.infinite(a)==F && b==Inf) {
    x = runif(n,pgumbel(a),1)
  } else {
    x = runif(n,pgumbel(a),pgumbel(b))
  }
  -log(-log(x))
}

rmgumbel = function(n,coef,max.loc=1,burn.in=50) {
  m=length(coef)
  n.order = c(max.loc,c(1:m)[-max.loc])
  n.coef = coef[n.order]
  N = n + burn.in + 1
  mcmc=matrix(0,nrow=N,ncol=m)
  mcmc[1,1] =  rgumbel(1)
  for (i in (2:m)) {
    mcmc[1,i] = rgumbel(1,b = n.coef[1] - n.coef[i] + mcmc[1,1])
  }
  for (i in (2:N)) {
    Max=max(-n.coef[1] + mcmc[i-1,2:m] + n.coef[2:m])
    mcmc[i,1] = rgumbel(1,a=Max)
    for (j in (2:m)) {
      mcmc[i,j] = rgumbel(1,b = n.coef[1] - n.coef[j] + mcmc[i,1])
    }
  }
  mcmc=mcmc[-(1:(burn.in+1)),,drop=F]
  if(max.loc==1){
    mcmc
  } else if (max.loc == m){
    reorder=c(2:m,1)
    mcmc[,reorder]
  } else {
    reorder=c(2:max.loc,1,(max.loc+1):m)
    mcmc[,reorder]
  }
}

## Pearson and Deviance Residual for Multinomial logit regression
res.mlogit=function(x,type="pearson") {
  if (type=="pearson") {
    pro=x$probabilities
    Var=pro*(1-pro)
    #Var=pro
    resp=matrix(x$model$y,ncol=dim(pro)[2],byrow=T)
    res=(resp-pro)/sqrt(Var)
    return(res)
  } else {
    pro=x$probabilities
    resp=matrix(x$model$y,ncol=dim(pro)[2],byrow=T)
    Sign=2*(resp-0.5)
    dev=resp*log(pro) + (1-resp)*log(1-pro)
    res=Sign*sqrt(-2*dev)
    return(res)
  }
}

boot.surres=function(mymodel,K=10) {
  res.mat=matrix(NA,ncol=K,nrow=dim(mymodel$residuals)[1])
  for (i in (1:K)) {
    res.mat[,i]=surrogate_resid(mymodel)[,1]
  }
  res.mat
}

#######################################################################################
#                                                                                     #
# Section 2: Numerical Examples                                                       #
#                                                                                     #
#######################################################################################

#############################################################################
# example 1
#############################################################################

## generate data
set.seed(1)
x=runif(2000,-5,5)
prob=exp(1+2*x)/(1+exp(1+2*x))
y=as.numeric(runif(1:2000)<prob)
eg1=data.frame(y=y,x=x)
m1=glm(y~x,family=binomial(link = "logit"),data=eg1)
## pearson, deviance and surrogate residuals
means = cbind(1,eg1[,2]) %*% coef(m1)
r.p=residuals(m1,type="pearson")
r.d=residuals(m1,type="deviance")
r.s=surrogate_resid(m1,means)
## plots
par(mfrow=c(3,2))
{
  plot(eg1$x,r.p,xlab="X",main="(a)",ylab="Pearson's Residual",cex.lab=1.5)
  l1=loess(r.p~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  qqnorm(r.p, xlab="Normal Distribution",ylab="Pearson's residual",main="(b)",cex.lab=1.5);
  abline(0,1,col=2,lwd=2)
  
  plot(eg1$x,r.d,xlab="X",main="(c)",ylab="Deviance Residual",cex.lab=1.5)
  l1=loess(r.d~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  qqnorm(r.d, xlab="Normal Distribution",main="(d)",ylab="Deviance Residual",cex.lab=1.5);
  abline(0,1,col=2,lwd=2)
  
  plot(eg1$x,r.s,xlab="X",main="(e)",ylab="Surrogate Residual",cex.lab=1.5)
  l1=loess(r.p~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  plot(qlogis((1:2000)/(2000+1)),sort(r.s),xlab="Logistic Distribution",
       ylab="Surrogate Residual",main="(f)",cex.lab=1.5)
  abline(0,1,col=2,lwd=2)
}

#############################################################################
# example 2
#############################################################################

## generate data
set.seed(76453)
n=2000
x=runif(n,-5,5)
x1= 1+7*x-2*x^2  + rgumbel(n)
x2= 8+2*x - 4*x^2 + rgumbel(n)
x3= 0+0*x+0*x^2 + rgumbel(n)
y=rep(3,n)
for ( i in 1:n) {
  if (x1[i]>x2[i] && x1[i]>x3[i]) {
    y[i]=1
  } else if (x2[i]>x1[i] && x2[i]>x3[i]) {
    y[i]=2
  }
}
z = y
y = as.factor(paste("y",y,sep=""))
eg2=data.frame(y=y,x=x)
mlogit.data <- mlogit.data(eg2, shape = "wide", choice = "y")
summary(mlogit(y ~ 0 | x + I(x^2),data=mlogit.data))
## the orginal model
m2=mlogit(y ~ 0 | x ,data=mlogit.data,reflevel = "y3")
## the improve model
m22=mlogit(y~0|x+I(x^2),data=mlogit.data,reflevel = "y3")
## surrogate residual
error=surrogate_resid(m22)
## plots
ylab.list=c(expression(r[1]),expression(r[2]),expression(r[3]))
ylab.list2=c(expression('Quantiles of '*r[1]),
             expression('Quantiles of '*r[2]),
             expression('Quantiles of '*r[3]))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
layout(matrix(1:6,nr=3,byrow=T))
par(mar=c(5, 5, 4, 2) + 0.1)
yname=c("y1","y2","y3")
for ( i in (1:3)) {
  plot(x,error[,yname[i]],ylab=ylab.list[i],xlab="X",main=main.list[2*(i-1)+1],ylim=c(-3,8),cex.lab=1.5);
  l1=loess(error[,yname[i]]~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  abline(v=-1,col=2,lty=2);abline(v=1,col=2,lty=2);abline(v=3.5,col=2,lty=2)
  
  plot(qgumbel((1:n)/(n+1)),sort(error[,yname[i]]+0.5772),xlab="Theoretical Quantiles",
       ylab=ylab.list2[i],main=main.list[2*i],cex.lab=1.5)
  abline(0,1,col="red")
}



# pearson and deviance residual
ylab.list=c(expression(r[1]^{(P)}),expression(r[2]^{(P)}),expression(r[3]^{(P)}))
ylab.list2=c(expression(r[1]^{(D)}),expression(r[2]^{(D)}),expression(r[3]^{(D)}))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
rp=res.mlogit(m2,type="pearson")
rd=res.mlogit(m2,type="deviance")
layout(matrix(1:6,nr=3,byrow=T))
par(mar=c(5, 5, 4, 2) + 0.1)
yname=c("y1","y2","y3")
for ( i in (1:3)) {
  plot(x,rp[,i],ylab=ylab.list[i],xlab="X",main=main.list[2*(i-1)+1],ylim=c(-3,8),cex.lab=1.5);
  l1=loess(rp[,i]~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  abline(v=-1,col=2,lty=2);abline(v=1,col=2,lty=2);abline(v=3.5,col=2,lty=2)
  
  plot(x,rd[,i],ylab=ylab.list2[i],xlab="X",main=main.list[2*i],ylim=c(-3,8),cex.lab=1.5);
  l1=loess(rd[,i]~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  abline(v=-1,col=2,lty=2);abline(v=1,col=2,lty=2);abline(v=3.5,col=2,lty=2)
}


#############################################################################
# example 3
#############################################################################

## generate data
set.seed(14567)
n=2000
x.y1 = runif(n,-5,5)
x.y2 = runif(n,-5,5)
x.y3 = runif(n,-5,5)
x1= 0 + 10*x.y1 +  rgumbel(n)
x2= 0 + 5*x.y2  + rgumbel(n)
x3= 0 - 10*x.y3 + rgumbel(n)
y=rep(3,n)
for ( i in 1:n) {
  if (x1[i]>x2[i] && x1[i]>x3[i]) {
    y[i]=1
  } else if (x2[i]>x1[i] && x2[i]>x3[i]) {
    y[i]=2
  }
}
z = y
y = as.factor(paste("y",y,sep=""))
eg3=data.frame(y,x.y1,x.y2,x.y3)
mlogit.data <- mlogit.data(eg3, shape = "wide", choice = "y",varying = 2:4)
## DCM model
m3 = mlogit(y ~ x,data=mlogit.data,reflevel = "y3")
## surrogate residual
error=surrogate_resid(m3)
means=matrix(0,ncol=3,nrow=2000)
means[,1]=coef(m3)[1] + coef(m3)[3]*x.y1
means[,2]=coef(m3)[2] + coef(m3)[3]*x.y2
means[,3]=0 + coef(m3)[3]*x.y3

layout(matrix(1:3,nr=3,byrow=T))
{
  x.list=data.frame(x.y1,x.y2,x.y3)
  plot(x.list[,1],coef(m3)[1] + coef(m3)[3]*x.y1+error[,"y1"],ylab=expression(s[1])
       ,xlab=expression(x[11]),main="(a)",cex.lab=1.5)
  uu=coef(m3)[1] + coef(m3)[3]*x.y1+error[,"y1"]
  l1=loess(uu~x.list[,1]);
  od <- order(x.list[,1]);lines(x.list[,1][od],l1$fitted[od],lwd=3,col="red")
  
  plot(x.list[,2],error[,"y2"]+coef(m3)[2] + coef(m3)[3]*x.y2,ylab=expression(s[2])
       ,xlab=expression(x[21]),main="(b)",cex.lab=1.5)
  uu=error[,"y2"]+coef(m3)[2] + coef(m3)[3]*x.y2
  l1=loess(uu~x.list[,2]);
  od <- order(x.list[,2]);lines(x.list[,2][od],l1$fitted[od],lwd=3,col="red")
  
  plot(x.list[,3],error[,"y3"]+0 + coef(m3)[3]*x.y3,ylab=expression(s[3])
       ,xlab=expression(x[31]),main="(c)",cex.lab=1.5)
  uu=error[,"y3"]+0 + coef(m3)[3]*x.y3
  l1=loess(uu~x.list[,3]);
  od <- order(x.list[,3]);lines(x.list[,3][od],l1$fitted[od],lwd=3,col="red")
}

#### pearson and deviance residual
ylab.list=c(expression(r[1]^{(P)}),expression(r[2]^{(P)}),expression(r[3]^{(P)}))
ylab.list2=c(expression(r[1]^{(D)}),expression(r[2]^{(D)}),expression(r[3]^{(D)}))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
rp=res.mlogit(m3,type="pearson")
rd=res.mlogit(m3,type="deviance")

layout(matrix(1:6,nr=3,byrow=T))
{
  x.list=data.frame(x.y1,x.y2,x.y3)
  plot(x.list[,1],coef(m3)[1] + coef(m3)[3]*x.y1+rp[,1],ylab=expression(s[1]^{(P)})
       ,xlab=expression(x[11]),main="(a)",cex.lab=1.5)
  uu=coef(m3)[1] + coef(m3)[3]*x.y1+rp[,1]
  l1=loess(uu~x.list[,1]);
  od <- order(x.list[,1]);lines(x.list[,1][od],l1$fitted[od],lwd=3,col="red")
  
  plot(x.list[,1],coef(m3)[1] + coef(m3)[3]*x.y1+rd[,1],ylab=expression(s[1]^{(D)})
       ,xlab=expression(x[11]),main="(b)",cex.lab=1.5)
  uu=coef(m3)[1] + coef(m3)[3]*x.y1+rd[,1]
  l1=loess(uu~x.list[,1]);
  od <- order(x.list[,1]);lines(x.list[,1][od],l1$fitted[od],lwd=3,col="red")
  #
  plot(x.list[,2],rp[,2]+coef(m3)[2] + coef(m3)[3]*x.y2,ylab=expression(s[2]^{(P)})
       ,xlab=expression(x[21]),main="(c)",cex.lab=1.5)
  uu=rp[,2]+coef(m3)[2] + coef(m3)[3]*x.y2
  l1=loess(uu~x.list[,2]);
  od <- order(x.list[,2]);lines(x.list[,2][od],l1$fitted[od],lwd=3,col="red")
  
  plot(x.list[,2],rd[,2]+coef(m3)[2] + coef(m3)[3]*x.y2,ylab=expression(s[2]^{(D)})
       ,xlab=expression(x[21]),main="(d)",cex.lab=1.5)
  uu=rd[,2]+coef(m3)[2] + coef(m3)[3]*x.y2
  l1=loess(uu~x.list[,2]);
  od <- order(x.list[,2]);lines(x.list[,2][od],l1$fitted[od],lwd=3,col="red")
  #
  plot(x.list[,3],rp[,3]+0 + coef(m3)[3]*x.y3,ylab=expression(s[3]^{(P)})
       ,xlab=expression(x[31]),main="(e)",cex.lab=1.5)
  uu=rp[,3]+0 + coef(m3)[3]*x.y3
  l1=loess(uu~x.list[,3]);
  od <- order(x.list[,3]);lines(x.list[,3][od],l1$fitted[od],lwd=3,col="red")
  
  plot(x.list[,3],rd[,3]+0 + coef(m3)[3]*x.y3,ylab=expression(s[3]^{(D)})
       ,xlab=expression(x[31]),main="(f)",cex.lab=1.5)
  uu=rd[,3]+0 + coef(m3)[3]*x.y3
  l1=loess(uu~x.list[,3]);
  od <- order(x.list[,3]);lines(x.list[,3][od],l1$fitted[od],lwd=3,col="red")
}

#############################################################################
# example 4
#############################################################################

## generate data
set.seed(76453)
n=2000
x=runif(n,-5,5)
z0=runif(n)>0.5
x1= -10-5*x+5*z0+3*x*z0  + rgumbel(n)
x2= -10+5*x-5*z0+6*x*z0 + rgumbel(n)
x3= 0 + rgumbel(n)
y=rep(3,n)
for ( i in 1:n) {
  if (x1[i]>x2[i] && x1[i]>x3[i]) {
    y[i]=1
  } else if (x2[i]>x1[i] && x2[i]>x3[i]) {
    y[i]=2
  }
}
z = y
y = as.factor(paste("y",y,sep=""))
eg4=data.frame(y,x=x,z0=z0)
mlogit.data <- mlogit.data(eg4, shape = "wide", choice = "y")
summary(mlogit(y ~ 0 | x+ z0,data=mlogit.data))
## DCM models
m41=mlogit(y ~ 0 | x ,data=mlogit.data[which(mlogit.data$z0==0),],reflevel = "y3")
m42=mlogit(y ~ 0 | x ,data=mlogit.data[which(mlogit.data$z0==1),],reflevel="y3")
## surrogate residuals
res1=surrogate_resid(m41)
res2=surrogate_resid(m42)
## plots
ylab.list=c(expression(s[1]),expression(s[2]),expression(s[3]))
xlab.list1=c(expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'))
xlab.list2=c(expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
ylim.list=c(-60,25,-90,50,-5,10)
layout(matrix(1:6,nr=3,byrow=T))
coef1=m41$coefficients;coef2=m42$coefficients
for ( i in (1:3)) {
  if (i==1) {
    uu = coef1[1] + coef1[3]*x[which(z0==0)]+res1[,2]
  } else if (i==2) {
    uu = coef1[2] + coef1[4]*x[which(z0==0)]+res1[,3]
  } else {
    uu = res1[,1]
  }
  plot(x[which(z0==0)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+1],
       xlab=xlab.list1[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  l1=loess(uu~x[which(z0==0)]);
  od <- order(x[which(z0==0)]);lines(x[which(z0==0)][od],l1$fitted[od],lwd=3,col="red")
  
  
  if (i==1) {
    uu = coef2[1] + coef2[3]*x[which(z0==1)]+res2[,2]
  } else if (i==2) {
    uu = coef2[2] + coef2[4]*x[which(z0==1)]+res2[,3]
  } else {
    uu = res2[,1]
  }
  plot(x[which(z0==1)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+2],
       xlab=xlab.list2[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  l1=loess(uu~x[which(z0==1)]);
  od <- order(x[which(z0==1)]);lines(x[which(z0==1)][od],l1$fitted[od],lwd=3,col="red")
}

##### pearson residual
res1=res.mlogit(m41,type="pearson")
res2=res.mlogit(m42,type="pearson")
par(mar=c(5, 5, 4, 2) + 0.1)
ylab.list=c(expression(s[1]^{(P)}),expression(s[2]^{(P)}),expression(s[3]^{(P)}))
xlab.list1=c(expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'))
xlab.list2=c(expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
ylim.list=c(-100,100,-100,100,-50,50)
layout(matrix(1:6,nr=3,byrow=T))
coef1=m41$coefficients;coef2=m42$coefficients
for ( i in (1:3)) {
  if (i==1) {
    uu = coef1[1] + coef1[3]*x[which(z0==0)]+res1[,1]
  } else if (i==2) {
    uu = coef1[2] + coef1[4]*x[which(z0==0)]+res1[,2]
  } else {
    uu = res1[,3]
  }
  plot(x[which(z0==0)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+1],
       xlab=xlab.list1[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  #l1=loess(uu~x[which(z0==0)]);
  #od <- order(x[which(z0==0)]);lines(x[which(z0==0)][od],l1$fitted[od],lwd=3,col="red")
  
  
  if (i==1) {
    uu = coef2[1] + coef2[3]*x[which(z0==1)]+res2[,1]
  } else if (i==2) {
    uu = coef2[2] + coef2[4]*x[which(z0==1)]+res2[,2]
  } else {
    uu = res2[,3]
  }
  plot(x[which(z0==1)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+2],
       xlab=xlab.list2[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  #l1=loess(uu~x[which(z0==1)]);
  #od <- order(x[which(z0==1)]);lines(x[which(z0==1)][od],l1$fitted[od],lwd=3,col="red")
}


#### deviance residual

res1=res.mlogit(m41,type="deviance")
res2=res.mlogit(m42,type="deviance")
par(mar=c(5, 5, 4, 2) + 0.1)
ylab.list=c(expression(s[1]^{(D)}),expression(s[2]^{(D)}),expression(s[3]^{(D)}))
xlab.list1=c(expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'),
             expression(x[1]*'(Control)'))
xlab.list2=c(expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'),
             expression(x[1]*'(Treatment)'))
main.list=c("(a)","(b)","(c)","(d)","(e)","(f)")
ylim.list=c(-60,25,-90,50,-5,10)
layout(matrix(1:6,nr=3,byrow=T))
coef1=m41$coefficients;coef2=m42$coefficients
for ( i in (1:3)) {
  if (i==1) {
    uu = coef1[1] + coef1[3]*x[which(z0==0)]+res1[,1]
  } else if (i==2) {
    uu = coef1[2] + coef1[4]*x[which(z0==0)]+res1[,2]
  } else {
    uu = res1[,3]
  }
  plot(x[which(z0==0)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+1],
       xlab=xlab.list1[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  #l1=loess(uu~x[which(z0==0)]);
  #od <- order(x[which(z0==0)]);lines(x[which(z0==0)][od],l1$fitted[od],lwd=3,col="red")
  
  
  if (i==1) {
    uu = coef2[1] + coef2[3]*x[which(z0==1)]+res2[,1]
  } else if (i==2) {
    uu = coef2[2] + coef2[4]*x[which(z0==1)]+res2[,2]
  } else {
    uu = res2[,3]
  }
  plot(x[which(z0==1)],uu,
       ylab=ylab.list[i],main=main.list[2*(i-1)+2],
       xlab=xlab.list2[i],ylim=ylim.list[(2*i-1):(2*i)],cex.lab=1.5);
  #l1=loess(uu~x[which(z0==1)]);
  #od <- order(x[which(z0==1)]);lines(x[which(z0==1)][od],l1$fitted[od],lwd=3,col="red")
}

##########################################################################

#############################################################################
# Bootstrapping residual
#############################################################################

## generate data
set.seed(76453)
n=50
x=runif(n,-5,5)
x1= 1+7*x-2*x^2  + rgumbel(n)
x2= 8+2*x - 4*x^2 + rgumbel(n)
x3= 0+0*x+0*x^2 + rgumbel(n)
y=rep(3,n)
for ( i in 1:n) {
  if (x1[i]>x2[i] && x1[i]>x3[i]) {
    y[i]=1
  } else if (x2[i]>x1[i] && x2[i]>x3[i]) {
    y[i]=2
  }
}
z = y
y = as.factor(paste("y",y,sep=""))
eg2=data.frame(y=y,x=x)
mlogit.data <- mlogit.data(eg2, shape = "wide", choice = "y")
summary(mlogit(y ~ 0 | x + I(x^2),data=mlogit.data))
## the orginal model
m2=mlogit(y ~ 0 | x + I(x^2),data=mlogit.data,reflevel = "y3")

## surrogate residual
set.seed(105)
main.list1=c("(a) Replication 1", "(c) Replication 2", "(e) Replication 3")
main.list2=c("(b) Replication 1", "(d) Replication 2", "(f) Replication 3")
layout(matrix(1:6,nr=3,byrow=T))
for (i in (1:3)) {
  error=surrogate_resid(m2)[,1]
  plot(x,error,ylab=expression(r[1]),xlab="X",main=main.list1[i],ylim=c(-3,8),cex.lab=1.5);
  l1=loess(error~x);
  od <- order(x);lines(x[od],l1$fitted[od],lwd=3,col="red")
  
  br=boot.surres(m2,K=10)
  mat.res=matrix(0,ncol=2,nrow=100*10)
  mat.res[,1]=rep(x,times=10)
  mat.res[,2]=as.vector(br)
  plot(mat.res[,1],mat.res[,2],ylim=c(-3,8),xlab="X",main=main.list2[i],ylab=expression(r[1]^{(B)}),cex.lab=1.5)
  l1=loess(mat.res[,2]~mat.res[,1]);
  od <- order(mat.res[,1]);lines(mat.res[,1][od],l1$fitted[od],lwd=3,col="red")
}








