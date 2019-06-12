
###Calculate the raw trend pre-standardisation
####convert data to log space using delta log normal distribution
deltalognorm.mean=function(x){
  n=length(x)
  m=length(x[x>0])
  pnz=m/n
  xnz=subset(x, x>0)
  mnlogNZ=mean(log(xnz))
  sdlogNZ=sd(log(xnz))
  #selogNZ=sdlogNZ/sqrt(m)
  theta=log(pnz)+mnlogNZ+sdlogNZ^2/2
  cc=(1-pnz)^(n-1)
  dd=1+(n-1)*pnz
  vartheta=((dd-cc)*(1-cc*dd)-m*(1-cc)^2)/(m*(1-cc*dd)^2)+sdlogNZ^2/m+sdlogNZ^4/(2*(m-1))
  lowCL=exp(theta-1.96*sqrt(vartheta))
  uppCL=exp(theta+1.96*sqrt(vartheta))
  return(c(exp(theta), lowCL, uppCL))
}

#Calculate Mean raw cpue across fishery
cpueagg.sub=data.frame(year=allyears, cpue=NA, lowCL=NA, uppCL=NA)
for (iyear in allyears){
  subdat.sub=subset(idata, year==iyear)
  cpueagg.sub[cpueagg.sub$year==iyear, -1]=deltalognorm.mean(subdat.sub$cpue)
}
cpueagg.sub

##export out raw cpue trend
write.csv(cpueagg.sub, paste(DirOut,"Redfishes raw data cpue trend", ".csv", sep=""), row.names = F)

##=====================================================================================================
##===========STEP 3 - Standardisation Model============================================================
##=====================================================================================================

##========================================================
## HURDLE MODEL ANALYSIS
# Binomial GLM (Cpue = 0 and Cpue > 0)
# Linear model (Cpue >0)
# Combine results of both the Binomial and Linear models
# Run bootstrap to get correct confidence intervals 
##========================================================


###Run model by calendar year
###choose model factors
modelRHS="fyear+fmonth+as.factor(block30)+as.factor(gear)"
#modelRHS="fyear+fmonth+as.factor(block30)+as.factor(gear)+as.factor(TargetBCG)"
##distinguish between zero and positive cpue for modelling
idata$ind=0
idata$ind[idata$cpue>0]=1

#####binomial GLM##### 
model.formula.glm1=paste("glm1=glm(ind~",modelRHS, ", data=idata, family=binomial(link='logit'),na.action=na.exclude)", sep="")
eval(parse(text=model.formula.glm1))
summary(glm1)
anova(glm1, test="Chi")
drop1(glm1, .~., test="Chi")
lsm1=summary(lsmeans(glm1, ~fyear, type="response"))
lsm1$fyear=as.numeric(as.character(lsm1$fyear))

idata$pred.prop=predict(glm1,type="response")

#####Linear Model (Cpue >0) #####
idatax=subset(idata, ind==1)
### log normal distribution GLM with cpue>0
model.formula.glm2=paste("glm2=lm(log(cpue)~",modelRHS, ", data=idatax)", sep="")
eval(parse(text=model.formula.glm2))
summary(glm2)
coef(glm2)
anova(glm2)#type 1 error
drop1(glm2, test="F")#type 3 error
exp(coef(glm2))
idata$pred.nonzero=exp(predict(glm2,newdata=idata))
idata$pred.hurdle=idata$pred.prop * idata$pred.nonzero
head(idata)
sigma.glm2=summary(glm2)$sigma
lsm2=summary(lsmeans(glm2, ~fyear, type="link"))
lsm2$response=exp(lsm2$lsmean)*exp(sigma.glm2^2/2)
lsm2$asymp.LCL=exp(lsm2$lower.CL)*exp(sigma.glm2^2/2)
lsm2$asymp.UCL=exp(lsm2$upper.CL)*exp(sigma.glm2^2/2)

#####Combine Binomial GLM with LM Cpue>0
lsm2$fyear=as.numeric(as.character(lsm2$fyear))
lsmB=cbind(lsm1[, c("fyear","prob", "asymp.LCL", "asymp.UCL")], 
           lsm2[, c("response", "asymp.LCL", "asymp.UCL")])
lsmB$index=lsmB$response*lsmB$prob
names(lsmB)[c(3,4, 6,7)]=c("plow", "pupp", "rlow", "rupp")
lsmB$low=lsmB$rlow*lsmB$plow
lsmB$upp=lsmB$rupp*lsmB$pupp
lsmB

###boot strap mean to obtain correct confidence intervals
library(boot)
hurdle_fn=function(data, i) {
  dat_boot=data[i, ]
  model.formula.m1=paste("m1=glm(ind~",modelRHS, ", data=dat_boot, family=binomial(link='logit'),na.action=na.exclude)", sep="")
  eval(parse(text=model.formula.m1))
  model.formula.m2=paste("m2=lm(log(cpue)~",modelRHS, ", data=subset(dat_boot, ind==1))", sep="")
  eval(parse(text=model.formula.m2))
  bin_coefs=summary(lsmeans(m1, ~fyear, type="response"))$prob
  sigma.m2=summary(m2)$sigma
  lsm_lognormal=summary(lsmeans(m2, ~fyear, type="link"))
  lognormal_coefs=exp(lsm_lognormal$lsmean)*exp(sigma.m2^2/2)
  return((bin_coefs*lognormal_coefs))
}

boot.out=boot(idata, hurdle_fn, R=1000, stype="i")

hurdle.conf.int=data.frame(fyear=allyears, index=boot.out$t0)
hurdle.conf.int$lowCL=apply(boot.out$t, 2, function(x)(quantile(x, 0.025)))
hurdle.conf.int$uppCL=apply(boot.out$t, 2, function(x)(quantile(x, 0.975)))
hurdle.conf.int


###OUTPUT RESULTS#######
file.out=T
filename.out=paste(DirOut,"Redfishes model outputs ",format(Sys.Date(), "%Y_%m_%d"), ".txt", sep="")

if (file.out) sink(file=filename.out) 
if (file.out) print("Binomial Model")
if (file.out) summary(glm1)
if (file.out) anova(glm1, test="Chi")
if (file.out) drop1(glm1, .~., test="Chi")
if (file.out) print("Linear Model")
if (file.out) summary(glm2)
if (file.out) anova(glm2)#type 1 error
if (file.out) drop1(glm2, test="F")#type 3 error
if (file.out) sink(NULL)

##export out standardised trend
write.csv(hurdle.conf.int, paste(DirOut,"Redfishes bootstrapped mean ",format(Sys.Date(), "%Y_%m_%d"), ".csv", sep=""),row.names=F)
