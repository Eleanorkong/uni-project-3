setwd("~/Desktop/MT3508 project")
data<-read.table("MT3508 project data.txt",header=T)
#check relationship between day and weights
plot(data$day,data$weights,pch=19,cex=0.8,xlab="Days",ylab="Weights")

#Model 1: exp(y)=beta0+beta1*x1(first)+beta2*x2(second)+beta3*x3(third)
#AIC: 333.2507
data1=data[c(1:284),c(1:4)]
data1=na.omit(data1)

negllik1=function(trtheta,data1){
  n=dim(data1)[1]
  mu=trtheta[1]+trtheta[2]*data1$first+trtheta[3]*data1$second+trtheta[4]*data1$third
  sigma2=exp(trtheta[5])
  ll=-n/2*log(2*pi)-n/2*log(sigma2)-sum((data1$weights-mu)^2)/(2*sigma2) 
  return(-ll)
}

start=c(rep(0,4),log(var(data1$weights)))

negllik1(start,data1)

mle1.opt=optim(start,negllik1,data1=data1,hessian=TRUE,method="BFGS")

mu=c(mle1.opt$par[1:4])
mu
# [1] 136.83279729  -0.03712647   0.15555299   0.17664981

sigma2=exp(mle1.opt$par[5])
sigma2
# [1] 0.5320047

#95% CI using Hessian 
H=mle1.opt$hessian
vcv=solve(H)
beta0.ci=mle1.opt$par[1]+c(-1.96,1.96)*sqrt(vcv[1,1])
beta0.ci
# [1] 136.2225 137.4431
beta1.ci=mle1.opt$par[2]+c(-1.96,1.96)*sqrt(vcv[2,2])
beta1.ci
# [1] -0.04635698 -0.02789597
beta2.ci=mle1.opt$par[3]+c(-1.96,1.96)*sqrt(vcv[3,3])
beta2.ci
# [1] 0.1508258 0.1602802
beta3.ci=mle1.opt$par[4]+c(-1.96,1.96)*sqrt(vcv[4,4])
beta3.ci
# [1] 0.1710829 0.1822167
logsigma2.ci=mle1.opt$par[5]+c(-1.96,1.96)*sqrt(vcv[5,5])
exp(logsigma2.ci)
# [1] 0.4245808 0.6666079

#AIC 
AIC.1=2*mle1.opt$value+2*length(mle1.opt$pars)
AIC.1

#Expected values
n=dim(data1) [1]
X=matrix(c(rep(1,n),data1$first,data1$second,data1$third),ncol=4)
head(X)
E.y=X%*%mle1.opt$par[1:4]

#Residuals
library(TSA)
resids=data1$weights-E.y

#Run test on residuals
runs(resids,0)$pvalue 
plot(density(resids))

#Goodness-of-fits (KS on residuals)
ks.test(resids,"pnorm",0,sqrt(var(resids)),alternative="two.sided")

cvm.test(resids,"pnorm")
ad.test(resids,"pnorm")


#residuals against fitted
pos=(resids)>0
col=rep("blue",length(resids))
col[pos]="red"
plot(E.y,resids,col=col,pch=19)
abline(0,0,lty=2)

#QQ-plot
qqnorm(resids)
qqline(resids,lty=2)

#glm for model 1
fit1=glm(weights~first+second+third,family=gaussian,data=data1)
summary(fit1)
anova(fit1,test="F")

1-pchisq(deviance(fit1),df.residual(fit1))
plot(residuals(fit1,type="response")~predict(fit1,type="response"),xlab=expression(hat(mu)),ylab="Raw residuals")
plot(residuals(fit1,type="deviance")~predict(fit1,type="response"),xlab=expression(hat(mu)),ylab="Deviance residuals")
plot(residuals(fit1,type="pearson")~predict(fit1,type="response"),xlab=expression(hat(mu)),ylab="Pearson residuals")

#lm for model 1
fit1lm=lm(weights~first+second+third,family=gaussian,data=data1)
summary(fit1lm)
ols_plot_resid_qq(fit1lm)
ols_test_normality(fit1lm)
ols_test_correlation(fit1lm)
ols_plot_resid_fit(fit1lm)
ols_plot_resid_hist(fit1lm)
plot(fit1lm)


#Model 2: exp(y)=beta0+beta1*x1(first)+beta23*x23(second and third)
#AIC: 352.5245
data2<-read.table("MT3508 project data 2.txt",header=T)
data2=data2[c(1:284),c(1:3)]
data2=na.omit(data2)

negllik2=function(trtheta,data2){
  n=dim(data2)[1]
  mu=trtheta[1]+trtheta[2]*data2$first+trtheta[3]*data2$second
  sigma2=exp(trtheta[4])
  ll=-n/2*log(2*pi)-n/2*log(sigma2)-sum((data2$weights-mu)^2)/(2*sigma2) 
  return(-ll)
}

start=c(rep(0,3),log(var(data2$weights)))

negllik2(start,data2)

mle2.opt=optim(start,negllik2,data2=data2,hessian=TRUE,method="BFGS")

mu=c(mle2.opt$par[1:3])
mu
#[1] 137.27823975  -0.04672614   0.16500378

sigma2=exp(mle2.opt$par[4])
sigma2
# [1] 0.6045461

#95% CI using Hessian 
H=mle2.opt$hessian
vcv=solve(H)
beta0.ci=mle2.opt$par[1]+c(-1.96,1.96)*sqrt(vcv[1,1])
beta0.ci
# [1] 136.6611 137.8954
beta1.ci=mle2.opt$par[2]+c(-1.96,1.96)*sqrt(vcv[2,2])
beta1.ci
# [1] -0.05551253 -0.03793974
beta2.ci=mle2.opt$par[3]+c(-1.96,1.96)*sqrt(vcv[3,3])
beta2.ci
# [1] 0.1624701 0.1675375
logsigma2.ci=mle2.opt$par[4]+c(-1.96,1.96)*sqrt(vcv[4,4])
exp(logsigma2.ci)
# [1] 0.4824644 0.7575189

#AIC 
AIC.2=2*mle2.opt$value+2*length(mle2.opt$pars)
AIC.2

#Expected values
n=dim(data2) [1]
X=matrix(c(rep(1,n),data2$first,data2$second),ncol=3)
head(X)
E.y=X%*%mle2.opt$par[1:3]

#Residuals
library(TSA)
resids=data2$weights-E.y

#Run test
runs(resids,0)$pvalue
plot(density(resids))

#Goodness-of-fits (KS on residuals)
ks.test(resids,"pnorm",alternative="two.sided")

#residuals against fitted
pos=(resids)>0
col=rep("blue",length(resids))
col[pos]="red"
plot(E.y,resids,col=col,pch=19)
abline(0,0,lty=2)

#QQ-plot
qqnorm(resids)
qqline(resids,lty=2)

#glm for model 2
fit2=glm(weights~first+second,family=gaussian,data=data2)
summary(fit2)
anova(fit2,test="F")

1-pchisq(deviance(fit2),df.residual(fit2))
plot(residuals(fit2,type="response")~predict(fit2,type="response"),xlab=expression(hat(mu)),ylab="Raw residuals")
plot(residuals(fit2,type="deviance")~predict(fit2,type="response"),xlab=expression(hat(mu)),ylab="Deviance residuals")
plot(residuals(fit2,type="pearson")~predict(fit2,type="response"),xlab=expression(hat(mu)),ylab="Pearson residuals")

#lm for model 2
fit2lm=lm(weights~first+second,family=gaussian,data=data2)
summary(fit2lm)
ols_plot_resid_qq(fit2lm)
ols_test_normality(fit2lm)
ols_test_correlation(fit2lm)
ols_plot_resid_fit(fit2lm)
ols_plot_resid_hist(fit2lm)
plot(fit2lm)


#Model3: based on model1 but add hour
#AIC: 327.9256
data1hr=data[c(1:284),c(1:5)]
data1hr=na.omit(data1hr)

negllik1hr=function(trtheta,data1hr){
  n=dim(data1hr)[1]
  mu=trtheta[1]+trtheta[2]*data1hr$first+trtheta[3]*data1hr$second+trtheta[4]*data1hr$third+trtheta[5]*data1hr$hours
  sigma2=exp(trtheta[6])
  ll=-n/2*log(2*pi)-n/2*log(sigma2)-sum((data1hr$weights-mu)^2)/(2*sigma2) 
  return(-ll)
}

start=c(rep(0,5),log(var(data1hr$weights)))

negllik1hr(start,data1hr)

mle1hr.opt=optim(start,negllik1hr,data1hr=data1hr,hessian=TRUE,method="BFGS")

mu=c(mle1hr.opt$par[1:5])
mu

sigma2=exp(mle1hr.opt$par[6])
sigma2

#95% CI using Hessian 
H=mle1hr.opt$hessian
vcv=solve(H)
beta0.ci=mle1hr.opt$par[1]+c(-1.96,1.96)*sqrt(vcv[1,1])
beta0.ci

beta1.ci=mle1hr.opt$par[2]+c(-1.96,1.96)*sqrt(vcv[2,2])
beta1.ci

beta2.ci=mle1hr.opt$par[3]+c(-1.96,1.96)*sqrt(vcv[3,3])
beta2.ci

beta3.ci=mle1hr.opt$par[4]+c(-1.96,1.96)*sqrt(vcv[4,4])
beta3.ci

beta4.ci=mle1hr.opt$par[5]+c(-1.96,1.96)*sqrt(vcv[5,5])
beta4.ci

logsigma2.ci=mle1hr.opt$par[6]+c(-1.96,1.96)*sqrt(vcv[6,6])
exp(logsigma2.ci)

#AIC
AIC.3=2*mle1hr.opt$value+2*length(mle1hr.opt$pars)
AIC.3

#Expected values
n=dim(data1hr) [1]
X=matrix(c(rep(1,n),data1hr$first,data1hr$second,data1hr$third,data1hr$hours),ncol=5)
head(X)
E.y=X%*%mle1hr.opt$par[1:5]

#Residuals
library(TSA)
resids=data1hr$weights-E.y

#Run test
runs(resids,0)$pvalue
plot(density(resids))

#Goodness-of-fits (KS on residuals)
ks.test(resids,"pnorm",alternative="two.sided")

#residuals against fitted
pos=(resids)>0
col=rep("blue",length(resids))
col[pos]="red"
plot(E.y,resids,col=col,pch=19)
abline(0,0,lty=2)

#QQ plot
qqnorm(resids)
qqline(resids)

#glm for model 3
fit3=glm(weights~first+second+third+hours,family=gaussian,data=data1hr)
summary(fit3)

1-pchisq(deviance(fit3),df.residual(fit3))
plot(residuals(fit3,type="response")~predict(fit3,type="response"),xlab=expression(hat(mu)),ylab="Raw residuals")
plot(residuals(fit3,type="deviance")~predict(fit3,type="response"),xlab=expression(hat(mu)),ylab="Deviance residuals")
plot(residuals(fit3,type="pearson")~predict(fit3,type="response"),xlab=expression(hat(mu)),ylab="Pearson residuals")
anova(fit3,test="F")
confint(fit3)
anova(fit1,fit3,test="F")

#lm for model 3
fit3lm=lm(weights~first+second+third+hours,data=data1hr)
summary(fit3lm)
ols_plot_resid_qq(fit3lm)
ols_test_normality(fit3lm)
ols_test_correlation(fit3lm)
ols_plot_resid_fit(fit3lm)
ols_plot_resid_hist(fit3lm)
plot(fit3lm)

#LRT (model 1 and model 3)
llik1=-mle1.opt$value
llik2=-mle1hr.opt$value
LRT=2*(llik2-llik1)
LRT
pchisq(LRT,df=1,lower.tail=FALSE)


#Model4: based on model2 but add hour
#AIC: 349.178
data2hr<-read.table("MT3508 project data 2.txt",header=T)
data2hr=data2hr[c(1:284),c(1:4)]
data2hr=na.omit(data2hr)

negllik2hr=function(trtheta,data2hr){
  n=dim(data2hr)[1]
  mu=trtheta[1]+trtheta[2]*data2hr$first+trtheta[3]*data2hr$second+trtheta[4]*data2hr$hours
  sigma2=exp(trtheta[5])
  ll=-n/2*log(2*pi)-n/2*log(sigma2)-sum((data2hr$weights-mu)^2)/(2*sigma2) 
  return(-ll)
}

start=c(rep(0,4),log(var(data2hr$weights)))

negllik2hr(start,data2hr)

mle2hr.opt=optim(start,negllik2hr,data2hr=data2hr,hessian=TRUE,method="BFGS")
mle2hr.opt$counts

mu=c(mle2hr.opt$par[1:4])
mu

sigma2=exp(mle2hr.opt$par[5])
sigma2


#95% CI using Hessian 
H=mle2hr.opt$hessian
vcv=solve(H)
beta0.ci=mle2hr.opt$par[1]+c(-1.96,1.96)*sqrt(vcv[1,1])
beta0.ci

beta1.ci=mle2hr.opt$par[2]+c(-1.96,1.96)*sqrt(vcv[2,2])
beta1.ci

beta2.ci=mle2hr.opt$par[3]+c(-1.96,1.96)*sqrt(vcv[3,3])
beta2.ci

beta3.ci=mle2hr.opt$par[4]+c(-1.96,1.96)*sqrt(vcv[4,4])
beta3.ci

logsigma2.ci=mle2hr.opt$par[5]+c(-1.96,1.96)*sqrt(vcv[5,5])
exp(logsigma2.ci)

#AIC
AIC.4=2*mle2hr.opt$value+2*length(mle2hr.opt$pars)
AIC.4

#Expected values 
n=dim(data2hr) [1]
X=matrix(c(rep(1,n),data2hr$first,data2hr$second,data2hr$hours),ncol=4)
head(X)
E.y=X%*%mle2hr.opt$par[1:4]

#Residuals
library(TSA)
resids=data2hr$weights-E.y

#Run test
runs(resids,0)$pvalue
plot(density(resids))

#Goodness-of-fits (KS on residuals)
ks.test(resids,"pnorm",alternative="two.sided")

#residuals vs fitted 
pos=(resids)>0
col=rep("blue",length(resids))
col[pos]="red"
plot(E.y,resids,col=col,pch=19)
abline(0,0,lty=2)

#QQ plot
qqnorm(resids)
qqline(resids)

#glm for model 4
fit4=glm(weights~first+second+hours,family=gaussian,data=data2hr)
summary(fit4)
anova(fit4,test="F")

1-pchisq(deviance(fit4),df.residual(fit4))
plot(residuals(fit4,type="response")~predict(fit4,type="response"),xlab=expression(hat(mu)),ylab="Raw residuals")
plot(residuals(fit4,type="deviance")~predict(fit4,type="response"),xlab=expression(hat(mu)),ylab="Deviance residuals")
plot(residuals(fit4,type="pearson")~predict(fit4,type="response"),xlab=expression(hat(mu)),ylab="Pearson residuals")

#lm for model 4
fit4lm=lm(weights~first+second+hours,data=data1hr)
summary(fit4lm)
ols_plot_resid_qq(fit4lm)
ols_test_normality(fit4lm)
ols_test_correlation(fit4lm)
ols_plot_resid_fit(fit4lm)
ols_plot_resid_hist(fit4lm)
plot(fit4lm) 

