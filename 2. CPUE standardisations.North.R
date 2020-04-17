#--------- CPUE STANDARDISATIONS OF Northern Shark Fisheries ---------#
# notes: Heupel & McAuley 2007 raised these issues about NSF cpue:
#           1. Demersal longline replaced pelagic gillnetting in the 1990s
#           2. 2002-03 a large ex-pelagic longliner entered fishery

#   Hence, for standardisations
                # use Method==LL, 
                # consider Vessel as a model term,

# Data selection: Removed years, months, vessels, with less than Min.yrs observations
#                 Kept species with at least Min.rec.per.yr for at least Min.yrs


# Response variable: cpue (kg/hook hours)
rm(list=ls(all=TRUE))

library(tidyverse)
library(cede)
library(mgcv)
library(mvtnorm)
library(doParallel)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)   


##############--- 1. DATA SECTION ---###################

fn.in=function(NM)read.csv(paste('C:/Matias/Analyses/Data_outs/',NM,sep=""))
Effort=fn.in(NM='Effort.monthly.NSF.csv')
Data=fn.in(NM='Data.monthly.NSF.csv')
Species.names=read.csv('C:/Matias/Data/Species_names_shark.only.csv')

##############--- 2. PARAMETERS SECTION ---###################
Min.rec.per.yr=10  #minimum records per year
Min.yrs=5          #minimum number of years with Min.rec.per.yr

Categorical=c("FINYEAR" ,"MONTH" ,"VESSEL")


##############--- 3. PROCEDURE SECTION ---###################

# Select unique effort records
Effort=Effort%>%
  filter(!is.na(hook.hours))%>%
  distinct(Same.return,.keep_all = T)

# Combine catch and effort for NSF
NSF.code=c('C051','C127','CL02','')  #fishery codes for NSF. Note, the '' code is form Rory's Table81.d data, which has no code
Data=Data%>%
    filter(METHOD%in%c('LL') & Reporter=="good" & FisheryCode%in%NSF.code
           & SPECIES<50000 & !SPECIES%in%c(22999,31000))%>%
  group_by(Same.return,FINYEAR,MONTH,VESSEL,SPECIES,BLOCKX,LAT,LONG,zone)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  left_join(Effort%>%dplyr::select(Same.return,hook.days,hook.hours),by="Same.return")%>%
  mutate(cpue=LIVEWT.c/hook.hours,
         yr=as.numeric(substr(FINYEAR,1,4)))%>%
  filter(!is.na(hook.hours))

# Define species and years to analyse
TAB=with(Data,table(SPECIES,FINYEAR))
TAB[TAB<Min.rec.per.yr]=0
TAB[TAB>=Min.rec.per.yr]=1
Keep.sp=names(which(rowSums(TAB)>=Min.yrs))
Keep.yrs=names(which(colSums(TAB)>=1))
Data=Data%>%filter(SPECIES%in%as.numeric(Keep.sp) & FINYEAR%in%Keep.yrs)


# One row per record
Data.w=Data%>%ungroup() %>%
      dplyr::select(-LIVEWT.c) %>%
      spread(SPECIES,cpue,fill=0)%>%
      data.frame
colnames(Data.w)[grepl("X",colnames(Data.w))]=gsub("X","",colnames(Data.w)[grepl("X",colnames(Data.w))])


#Nominal delta lognormal
Store.nominal.cpue=vector('list',length(Keep.sp))
names(Store.nominal.cpue)=Keep.sp
fn.nominal.delta.log=function(d)
{
  #remove years, months, vessels, with less than Min.yrs observations
  Kip.yr=table(d$FINYEAR)
  Kip.yr=names(Kip.yr[Kip.yr>=Min.yrs])
  Kip.mn=table(d$MONTH) 
  Kip.mn=as.integer(names(Kip.mn[Kip.mn>=Min.yrs]))
  Kip.vs=table(d$VESSEL)
  Kip.vs=names(Kip.vs[Kip.vs>=Min.yrs])
  
  d <- d%>%
    filter(FINYEAR%in%Kip.yr)%>%
    filter(MONTH%in%Kip.mn)%>%
    filter(VESSEL%in%Kip.vs)
  
  out = d %>%
    group_by(FINYEAR) %>%
    summarise(n = length(cpue),
              m = length(cpue[cpue>0]),
              mean.lognz = mean(log(cpue[cpue>0])),
              sd.lognz = sd(log(cpue[cpue>0]))) %>%
    mutate(p.nz = m/n,
           theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
           c = (1-p.nz)^(n-1),
           d = 1+(n-1)*p.nz,
           vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
             sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
           mean = exp(theta),
           lowCL = exp(theta - 1.96*sqrt(vartheta)),
           uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
    as.data.frame%>%
    dplyr::select(FINYEAR,mean,lowCL,uppCL)
  return(out)
}
for(s in 1:length(Keep.sp))
{
  this=Keep.sp[s]
  Store.nominal.cpue[[s]]=fn.nominal.delta.log(
            d=Data.w[,-match(Keep.sp[-match(this,Keep.sp)],
                    colnames(Data.w))]%>%rename(cpue=!!this))
}


# Standardised log normal
fn.delta=function(d,Formula.bi.gam,Formula.gam)   
{
  Bi <- d %>%mutate(catch.pos=as.numeric(cpue>0))
  
  #remove months, year-vessel with low observations
  Kip=Bi%>%filter(catch.pos>0)
  Kip.mn=table(Kip$MONTH) 
  Kip.mn=as.integer(names(Kip.mn[Kip.mn>=Min.yrs]))
  Kip.yr=table(Kip$FINYEAR)
  Kip.yr=names(Kip.yr[Kip.yr>=Min.yrs])
  Kip.vs=table(Kip$VESSEL)
  Kip.vs=names(Kip.vs[Kip.vs>=Min.yrs])
  
  Bi <- Bi%>%
    filter(FINYEAR%in%Kip.yr)%>%
    filter(MONTH%in%Kip.mn)%>%
    filter(VESSEL%in%Kip.vs)
  
  Kip=Bi%>%filter(catch.pos>0)      #some years have only 1 vessel, so remove those
  Kip.yr.vsl=table(Kip$FINYEAR,Kip$VESSEL)
  Kip.yr.vsl[Kip.yr.vsl>0]=1
  Kip.yr=names(which(rowSums(Kip.yr.vsl)>1))

  Bi <- Bi%>%filter(FINYEAR%in%Kip.yr)
  
  #get pos data
  d <- Bi%>%
    filter(catch.pos>0)%>%
    mutate(ln.cpue=log(cpue))
  
  id.fctr=which(labels(terms(Formula.gam))%in%Categorical)
  Bi=makecategorical(labels(terms(Formula.gam))[id.fctr],Bi)
  d=makecategorical(labels(terms(Formula.gam))[id.fctr],d)
  
  res.gam_bi <-gam(Formula.bi.gam,data=Bi, family="binomial",method="REML")
  res.gam <-gam(Formula.gam,data=d,method="REML")
  
  return(list(res.gam=res.gam,res.gam_bi=res.gam_bi,DATA=d,DATA_bi=Bi))
  
}
Stand.out=vector('list',length(Keep.sp))
names(Stand.out)=Keep.sp
Pred=Stand.out
for(s in 1:length(Keep.sp))
{
  this=Keep.sp[s]
  Form.gam=formula(ln.cpue~FINYEAR+MONTH+VESSEL+ s(LAT,LONG))
  if(this%in%c("18003","18007","18013")) Form.gam=formula(ln.cpue~FINYEAR+MONTH+VESSEL)
  Stand.out[[s]]=fn.delta(
                    d=Data.w[,-match(Keep.sp[-match(this,Keep.sp)],colnames(Data.w))]%>%rename(cpue=!!this),
                    Formula.bi.gam=formula(catch.pos~FINYEAR+MONTH+VESSEL+log(hook.hours)+ s(LAT,LONG)),
                    Formula.gam=Form.gam)
}


# Extract year and calculate uncertainty 
fn.MC.delta.cpue=function(BiMOD,MOD,BiData,PosData,niter,pred.term,ALL.terms)
{
  #don't use binomial part if no contrast (applicable to monthly records where catch is mostly positive)
  if(class(BiMOD)[1]=="gam") Bi.cont=quantile(summary(BiMOD)$se,probs=.6)
  if(class(BiMOD)[1]=="glm") Bi.cont=quantile(summary(BiMOD)$coefficients[,2],probs=.6)  
  
  #get terms
  ALL.terms=ALL.terms[which(ALL.terms%in%colnames(BiData))]
  Covar.bi=as.matrix(vcov(BiMOD))
  Covar.pos=as.matrix(vcov(MOD))
  dummy.BiMOD=BiMOD
  dummy.MOD=MOD
  knstnt.terms=ALL.terms[-match(pred.term,ALL.terms)]
  id.fctr=knstnt.terms[which(knstnt.terms%in%Categorical)]
  id.cont=knstnt.terms[which(!knstnt.terms%in%Categorical)]
  newdata.pos=matrix(nrow=1,ncol=length(knstnt.terms))
  colnames(newdata.pos)=c(id.fctr,id.cont) 
  newdata.pos=as.data.frame(newdata.pos)
  newdata.bi=newdata.pos
  for(ii in 1:ncol(newdata.pos))
  {
    if(colnames(newdata.pos)[ii]%in%id.fctr)
    {
      id=match(colnames(newdata.bi)[ii],names(BiData))
      if(!is.na(id))
      {
        dummy=sort(table(BiData[,id]))
        newdata.bi[,ii]= factor(names(dummy[length(dummy)]),levels(BiData[,id]))
      }
      id=match(colnames(newdata.pos)[ii],names(PosData))
      if(!is.na(id))
      {
        dummy=sort(table(PosData[,id]))
        newdata.pos[,ii]= factor(names(dummy[length(dummy)]),levels(PosData[,id]))
      }
    }
    if(colnames(newdata.pos)[ii]%in%id.cont)
    {
      id=match(colnames(newdata.bi)[ii],names(BiData))
      newdata.bi[,ii]= mean(BiData[,id])
      
      id=match(colnames(newdata.pos)[ii],names(PosData))
      newdata.pos[,ii]= mean(PosData[,id])
    }
  }
  nms.coef=names(unlist(dummy.coef(MOD)))
  pred.var=sapply(strsplit(nms.coef[grepl(pred.term, nms.coef)], paste(pred.term,".",sep="")), "[", 2)
  
  pred.dat.pos=data.frame(factor(pred.var,levels=pred.var))
  colnames(pred.dat.pos)=pred.term
  newdata.pos=cbind(pred.dat.pos,newdata.pos)
  
  newdata.bi=cbind(pred.dat.pos,newdata.bi)

  set.seed(999)
  
  Bi.pars.rand=rmvnorm(niter,mean=coef(BiMOD),sigma=Covar.bi)
  Pos.pars.rand=rmvnorm(niter,mean=coef(MOD),sigma=Covar.pos)
  
  MC.preds=matrix(nrow=niter,ncol=nrow(newdata.bi))
  
  
  for(n in 1:niter)
  {
    #Binomial part
    if(Bi.cont<5) dummy.BiMOD$coefficients=Bi.pars.rand[n,] else
      dummy.BiMOD=BiMOD
    newdata.bi$Pred.bi=predict(dummy.BiMOD,newdata=newdata.bi, type="response")
    
    #Positive part
    dummy.MOD$coefficients=Pos.pars.rand[n,]
    a=predict(dummy.MOD,newdata=newdata.pos, type="response",se.fit=T)
    newdata.pos$Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    
    dummy=left_join(newdata.bi,newdata.pos,by=pred.term)%>%
      #mutate(Index=Pred)%>%
      mutate(Index=Pred.bi*Pred)%>%
      dplyr::select(Index)%>%
      unlist
    
    MC.preds[n,]=dummy
  }
  
  #Get summary stats
  MEAN=colMeans(MC.preds,na.rm=T)
  SD=apply(MC.preds,2,sd,na.rm=T)
  LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T))
  UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T))
  
  Stats=cbind(pred.dat.pos,data.frame(MEAN=MEAN,SD=SD,LOW=LOW,UP=UP))
  Stats=Stats[order(Stats[,1]),]
  Stats=Stats%>%
        rename(Mean=MEAN, LOW.CI=LOW, UP.CI=UP)%>%
        mutate(CV=SD/Mean)
  rownames(Stats)=NULL
  return(Stats)
}
Niter=500   #MC interations
system.time({
  for(s in 1:length(Keep.sp))
  {
    Pred[[s]]=fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$res.gam_bi,
                               MOD=Stand.out[[s]]$res.gam,
                               BiData=Stand.out[[s]]$DATA_bi,
                               PosData=Stand.out[[s]]$DATA,
                               niter=Niter,
                               pred.term='FINYEAR',
                               ALL.terms=c(Categorical,"LAT","LONG","hook.hours"))
  }
})


##############--- 4. EXPORT DATA ---###################
setwd("C:/Matias/Analyses/Data_outs")

for (s in 1:length(Keep.sp))
{
  NM=Species.names$Name[match(Keep.sp[s],Species.names$SPECIES)]
  write.csv(Pred[[s]],paste(NM,".annual.abundance.NSF.csv",sep=""),row.names=F) 
}


##############--- 5. PLOT INDEX ---###################
setwd("C:/Matias/Analyses/Catch and effort/Outputs/NSF")
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))


Plot.cpue=function(cpue.stand,cpue.nom,Nrmlzd)
{
  tcl=.15
  cpue.stand=cpue.stand%>%
    mutate(Yrs=as.numeric(substr(FINYEAR,1,4)))
  cpue.nom=cpue.nom%>%
    mutate(Yrs=as.numeric(substr(FINYEAR,1,4)))%>%
    filter(Yrs%in%cpue.stand$Yrs)
  Yrs=cpue.stand$Yrs
  if(Nrmlzd=="YES")
  {
    cpue.stand=cpue.stand%>%
      mutate(Mn=mean(Mean),
             LOW.CI=LOW.CI/Mn,
              UP.CI=UP.CI/Mn,
              Mean=Mean/Mn)
    
    cpue.nom=cpue.nom%>%
      mutate(Mn=mean(mean),
             lowCL=lowCL/Mn,
             uppCL=uppCL/Mn,
             mean=mean/Mn)
  }
  ymax=max(c(cpue.stand$UP.CI,cpue.nom$uppCL),na.rm=T)

  
  plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
  with(cpue.stand,
       {
         points(Yrs, Mean, pch=19, lty=2, col="black",cex=1.5)
         arrows(x0=Yrs, y0=LOW.CI, 
                x1=Yrs, y1=UP.CI, 
                code=3, angle=90, length=0.05, col="black",lwd=1.5)
       })
  with(cpue.nom,
       {
         points(Yrs+tcl, mean, pch=19, lty=2, col="grey60",cex=1.5)
         arrows(x0=Yrs+tcl, y0=lowCL, 
                x1=Yrs+tcl, y1=uppCL, 
                code=3, angle=90, length=0.05, col="grey60",lwd=1.5)
       }) 
}

tiff(file="Figure 1. Standardised and nominal.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
smart.par(n.plots=length(Keep.sp),MAR=c(1,2.5,2.5,1),OMA=c(2.5,2,.1,.1),MGP=c(.1, 0.5, 0))

for(s in 1:length(Keep.sp))
{
  Plot.cpue(cpue.stand=Pred[[s]],cpue.nom=Store.nominal.cpue[[s]],Nrmlzd="YES")
  NM=Species.names$Name[match(Keep.sp[s],Species.names$SPECIES)]
  mtext(NM,side=3,line=.35,font=1,las=0,cex=1.35)
}
legend('topright',c("Standardised","Nominal"),bty='n',col=c("black","grey60"),pch=19,cex=1.5)
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()
