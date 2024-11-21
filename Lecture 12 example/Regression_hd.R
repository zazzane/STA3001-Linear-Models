set.seed(0)
n = 50
p = 30
x = matrix(rnorm(n*p),nrow=n)

bstar = c(runif(10,0.5,1),runif(20,0,0.3))
mu = as.numeric(x%*%bstar)

par(mar=c(4.5,4.5,0.5,0.5))
hist(bstar,breaks=30,col="gray",main="",
     xlab="True coefficients")

set.seed(1)
R = 100

fit = matrix(0,R,n)
err = numeric(R)

for (i in 1:R) {
  cat(c(i,", "))
  y = mu + rnorm(n)
  ynew = mu + rnorm(n)
  
  a = lm(y~x+0)
  bls = coef(a)
  fit[i,] = x%*%bls
  err[i] = mean((ynew-fit[i,])^2)
}

prederr = mean(err)

bias = sum((colMeans(fit)-mu)^2)/n
var = sum(apply(fit,2,var))/n

bias
var
p/n
1 + bias + var
prederr




set.seed(0)
n = 50
p = 30
x = matrix(rnorm(n*p),nrow=n)

bstar = c(runif(10,0.5,1),runif(20,0,0.3))
mu = as.numeric(x%*%bstar)

par(mar=c(4.5,4.5,0.5,0.5))
hist(bstar,breaks=30,col="gray",main="",
     xlab="True coefficients")

library(MASS)

set.seed(1)
R = 100
nlam = 60
lam = seq(0,25,length=nlam)

fit.ls = matrix(0,R,n)
fit.rid  = array(0,dim=c(R,nlam,n))
err.ls = numeric(R)
err.rid = matrix(0,R,nlam)

for (i in 1:R) {
  cat(c(i,", "))
  y = mu + rnorm(n)
  ynew = mu + rnorm(n)
  
  a = lm(y~x+0)
  bls = coef(a)
  fit.ls[i,] = x%*%bls
  err.ls[i] = mean((ynew-fit.ls[i,])^2)
  
  aa = lm.ridge(y~x+0,lambda=lam)
  brid = coef(aa)
  fit.rid[i,,] = brid%*%t(x)
  err.rid[i,] = rowMeans(scale(fit.rid[i,,],center=ynew,scale=F)^2)
}

aveerr.ls = mean(err.ls)
aveerr.rid = colMeans(err.rid)

bias.ls = sum((colMeans(fit.ls)-mu)^2)/n
var.ls = sum(apply(fit.ls,2,var))/n

bias.rid = rowSums(scale(apply(fit.rid,2:3,mean),center=mu,scale=F)^2)/n
var.rid = rowSums(apply(fit.rid,2:3,var))/n

mse.ls = bias.ls + var.ls
mse.rid = bias.rid + var.rid
prederr.ls = mse.ls + 1
prederr.rid = mse.rid + 1

bias.ls
var.ls
p/n

prederr.ls
aveerr.ls

cbind(prederr.rid,aveerr.rid)

par(mar=c(4.5,4.5,0.5,0.5))
plot(lam,prederr.rid,type="l",
     xlab="Amount of shrinkage",ylab="Prediction error")
abline(h=prederr.ls,lty=2)
text(c(1,24),c(1.48,1.48),c("Low","High"))
legend("topleft",lty=c(2,1),
       legend=c("Linear regression","Ridge regression"))

par(mar=c(4.5,4.5,0.5,0.5))
plot(lam,mse.rid,type="l",ylim=c(0,max(mse.rid)),
     xlab=expression(paste(lambda)),ylab="")
lines(lam,bias.rid,col="red")
lines(lam,var.rid,col="blue")
abline(h=mse.ls,lty=2)
legend("bottomright",lty=c(2,1,1,1),
       legend=c("Linear MSE","Ridge MSE","Ridge Bias^2","Ridge Var"),
       col=c("black","black","red","blue"))


