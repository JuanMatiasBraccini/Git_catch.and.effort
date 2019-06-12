#issue: response variable must by cpue

MOD=Pred.BaseCase
Grid.dat=Grid.data[[i]]
SPEC=SPECIES.vec[[i]]
  


impute.extract.yr.block=function(MOD,Grid.dat,SPEC)
{
  #1. Get models
  GLMbi=MOD$Bi
  GLMlog=MOD$LogN
  
  #Get coefficients
  Bi.coef=dummy.coef(GLMbi)
  Log.coef=dummy.coef(GLMlog)
  
  Bi.name=names(Bi.coef)
  Log.name=names(Log.coef)
  
  #2. Get NA coefficients requiring imputation
  Set.to.inter=c(paste(names(Log.coef$FINYEAR[1]),":",names(Log.coef$BLOCKX),sep=""),
                 paste(names(Log.coef$FINYEAR),":",names(Log.coef$BLOCKX[1]),sep=""))
  dummy=Log.coef$"FINYEAR:BLOCKX"
  Set.to.NA=subset(dummy,dummy==0 & !names(dummy)%in%Set.to.inter)
  a=ifelse(names(dummy)%in%names(Set.to.NA), NA,dummy)
  names(a)=names(dummy)
  Log.coef$"FINYEAR:BLOCKX"=a
  rm(a,dummy)
  
  
  #3. Extract year-block index for non-NA coefficients
  
    #3.1 binominal part (Maunder and Punt, 2004; Eq. 9)
  Bi.coef.no.yr.terms=Bi.coef[-match(c("(Intercept)","FINYEAR"),names(Bi.coef))]
  These=match(names(Bi.coef.no.yr.terms),names(Grid.dat$bin))
  names(These)=sapply(Grid.dat$bin[These],class)
  Grid.dat$bin$BLOCKX=names(rev(sort(table(MOD$BiData[,match("BLOCKX",names(MOD$BiData))]))))[1]
  
  biTemp <- rep(NA,length(Bi.coef.no.yr.terms)) # specify value of other terms
  for (o in 1:length(Bi.coef.no.yr.terms))
  {
    dmy=Bi.coef.no.yr.terms[[o]]
    if(names(These[o])=="numeric") biTemp[o]<- dmy*Grid.dat$bin[1,These[o]]
    if(names(These[o])=="factor")  biTemp[o]<-dmy[match(Grid.dat$bin[1,These[o]],names(dmy))]
  }
  BiSumCoeff <- Bi.coef$"(Intercept)" + Bi.coef$FINYEAR + sum(biTemp)  
  BiIndex <- exp(BiSumCoeff)/(1+exp(BiSumCoeff))  # inverse of logit link 
  
  
  #3.2 lognormal part  
  Log.coef.no.yr_blk.terms=Log.coef[-match(c("(Intercept)","FINYEAR:BLOCKX"),names(Log.coef))]
  These=match(names(Log.coef.no.yr_blk.terms),names(Grid.dat$pos))
  names(These)=sapply(Grid.dat$pos[These],class)
     
  biasCorr <- GLMlog$deviance/GLMlog$df.residual/2    # bias correction  
  

  LogTemp <- rep(NA,length(Log.coef.no.yr_blk.terms)) # specify value of other terms
  for (o in 1:length(Log.coef.no.yr_blk.terms))
  {
    dmy=Log.coef.no.yr_blk.terms[[o]]
    if(names(These[o])=="numeric") LogTemp[o]<- dmy*Grid.dat$pos[1,These[o]]
    if(names(These[o])=="factor")  LogTemp[o]<-dmy[match(Grid.dat$pos[1,These[o]],names(dmy))]
  }
  LogSumCoeff <- Log.coef$"(Intercept)" + Log.coef$"FINYEAR:BLOCKX" + sum(LogTemp) + biasCorr 
  LogIndex <- exp(LogSumCoeff)  # inverse of log 
  
  
  
  #4. Impute missing coefficients
  
}



