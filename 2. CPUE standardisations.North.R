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
NSF.code=c('C051','C127','CL02')  #fishery codes for NSF
Data=Data%>%
    filter(METHOD%in%c('LL') & Reporter=="good" & FisheryCode%in%NSF.code
           & SPECIES<50000 & !SPECIES%in%c(22999,31000))%>%
  group_by(Same.return,FINYEAR,MONTH,VESSEL,SPECIES,BLOCKX,LAT,LONG,zone)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  left_join(Effort%>%select(Same.return,hook.days,hook.hours),by="Same.return")%>%
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
      select(-LIVEWT.c) %>%
      spread(SPECIES,cpue,fill=0)%>%
      data.frame
colnames(Data.w)[grepl("X",colnames(Data.w))]=gsub("X","",colnames(Data.w)[grepl("X",colnames(Data.w))])


#Nominal delta lognormal
Store.nominal.cpue=vector('list',length(Keep.sp))
names(Store.nominal.cpue)=Keep.sp
fn.nominal.delta.log=function(d)
{
  Bi <- d %>%mutate(catch.pos=as.numeric(cpue>0))
  
  #remove years, months, vessels, with less than Min.yrs observations
  Kip=Bi%>%filter(catch.pos>0)
  Kip.yr=table(Kip$FINYEAR)
  Kip.yr=names(Kip.yr[Kip.yr>=Min.yrs])
  Kip.mn=table(Kip$MONTH) 
  Kip.mn=as.integer(names(Kip.mn[Kip.mn>=Min.yrs]))
  Kip.vs=table(Kip$VESSEL)
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
  return(d)
}
for(s in 1:length(Keep.sp))
{
  this=Keep.sp[s]
  Store.nominal.cpue[[s]]=fn.nominal.delta.log(d=Data.w[,-match(Keep.sp[-match(this,Keep.sp)],
                                                                colnames(Data.w))]%>%rename(cpue=!!this))
}


# Standardised log normal
fn.delta=function(d,Formula.bi.gam,Formula.gam)   
{
  Bi <- d %>%mutate(catch.pos=as.numeric(cpue>0))
  
  #remove years, months, vessels, with low observations
  Kip=Bi%>%filter(catch.pos>0)
  Kip.yr=table(Kip$FINYEAR)
  Kip.yr=names(Kip.yr[Kip.yr>=Min.yrs])
  Kip.mn=table(Kip$MONTH) 
  Kip.mn=as.integer(names(Kip.mn[Kip.mn>=Min.yrs]))
  Kip.vs=table(Kip$VESSEL)
  Kip.vs=names(Kip.vs[Kip.vs>=Min.yrs])
  
  Bi <- Bi%>%
      filter(FINYEAR%in%Kip.yr)%>%
      filter(MONTH%in%Kip.mn)%>%
      filter(VESSEL%in%Kip.vs)
  
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
Store.cpue=vector('list',length(Keep.sp))
names(Store.cpue)=Keep.sp
for(s in 1:length(Keep.sp))
{
  this=Keep.sp[s]
  Form.gam=formula(ln.cpue~FINYEAR+MONTH+VESSEL+ s(LAT,LONG))
  if(this%in%c("18003","18013")) Form.gam=formula(ln.cpue~FINYEAR+MONTH+VESSEL)
  Store.cpue[[s]]=fn.delta(
                    d=Data.w[,-match(Keep.sp[-match(this,Keep.sp)],colnames(Data.w))]%>%rename(cpue=!!this),
                    Formula.bi.gam=formula(catch.pos~FINYEAR+MONTH+VESSEL+ s(LAT,LONG)),
                    Formula.gam=Form.gam)
}


# Extract year and calculate uncertainty     #ACA



##############--- 4. EXPORT DATA ---###################



##############--- 5. PLOT INDEX ---###################

#Store.nominal.cpue