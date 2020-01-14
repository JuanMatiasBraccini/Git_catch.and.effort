# MISSING: Spatial closures and marine parks shape files and blockx within closures
#           calculate overlap


# Header ---------------------------------------------------------
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)
library(tidyverse)
library(mgcv)

# 2. Data section ---------------------------------------------------------
setwd('C:/Matias/Analyses/Data_outs')

  #TDGDLF catch and effort used for cpue standardisation
Data.monthly=read.csv("Data.monthly.GN.csv")
Data.daily=read.csv("Data.daily.GN.csv")
Effort.monthly=read.csv("Effort.monthly.csv")
Effort.daily=read.csv("Effort.daily.csv")

  #NSF catch and effort
Data.monthly.NSF=read.csv("Data.monthly.NSF.csv")
Data.daily.NSF=read.csv("Data.daily.NSF.csv")
Effort.monthly.NSF=read.csv("Effort.monthly.NSF.csv")
Effort.daily.NSF=read.csv("Effort.daily.NSF.csv")


  #Shape files. Spatial closures and marine parks       #MISSING

# 3. Parameters section ---------------------------------------------------------
Min.obs=100

# 4. Procedure section ---------------------------------------------------------

Fin.yr.mon=paste(1975:2005,substr(1976:2006,3,4),sep='-')

# 4.1 Manipulate data 
  #Catch
Data.monthly=Data.monthly%>%
                filter(SPECIES<=24900 & !SPECIES==22999 &    #sharks only
                         METHOD=="GN" & FINYEAR%in%Fin.yr.mon)%>%  
                mutate(lat.corner=LAT,
                       long.corner=LONG,
                       SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  #Keep species with a minimum of Min.obs
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return,VESSEL,Bioregion,zone,BLOCKX,LONG,LAT,
                       lat.corner,long.corner,SNAME,SPECIES,Sch.or.DogS,
                       FINYEAR,YEAR.c,MONTH,LIVEWT.c)%>%
                data.frame
            
Data.daily=Data.daily%>%
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="GN")%>%
                mutate(lat.corner=-(abs(as.numeric(substr(block10,1,2))+
                                          10*(as.numeric(substr(block10,3,3)))/60)),
                       long.corner=100+as.numeric(substr(block10,4,5))+
                                        10*(as.numeric(substr(block10,6,6)))/60,
                       SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return.SNo,Same.return,VESSEL,Bioregion,zone,block10,
                       BLOCKX,lat.corner,long.corner,SHOTS,
                       SNAME,SPECIES,FINYEAR,YEAR.c,MONTH,LIVEWT.c)%>%
                data.frame
                
Data.monthly.NSF=Data.monthly.NSF%>%
                filter(SPECIES<=24900 & !SPECIES==22999 
                       & METHOD=="LL" & FINYEAR%in%Fin.yr.mon)%>%
                mutate(lat.corner=LAT,
                       long.corner=LONG,
                       SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
        dplyr::select(Same.return,VESSEL,Bioregion,zone,BLOCKX,LONG,LAT,
                      lat.corner,long.corner,SNAME,SPECIES,FINYEAR,
                      YEAR.c,MONTH,LIVEWT.c)%>%
               data.frame
  
Data.daily.NSF=Data.daily.NSF%>%
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="LL")%>%
                mutate(lat.corner=-(abs(as.numeric(substr(block10,1,2))+
                                          10*(as.numeric(substr(block10,3,3)))/60)),
                       long.corner=100+as.numeric(substr(block10,4,5))+
                         10*(as.numeric(substr(block10,6,6)))/60,
                       SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return.SNo,Same.return,VESSEL,Bioregion,zone,block10,
                       BLOCKX,lat.corner,long.corner,SHOTS,
                       SNAME,SPECIES,FINYEAR,YEAR.c,MONTH,LIVEWT.c)%>%
                data.frame

  #Effort
Effort.monthly=Effort.monthly%>%
                      filter(Same.return%in%unique(Data.monthly$Same.return))%>%
               dplyr::select(Same.return,Km.Gillnet.Hours.c)%>%
                      group_by(Same.return)%>%
                      summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c))%>%
                      rename(Effort=Km.Gillnet.Hours.c)

Effort.daily=Effort.daily%>%
                      filter(Same.return.SNo%in%unique(Data.daily$Same.return.SNo))%>%
               dplyr::select(Same.return.SNo,Km.Gillnet.Hours.c)%>%
                      group_by(Same.return.SNo)%>%
                      summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c))%>%
                      rename(Effort=Km.Gillnet.Hours.c)

Effort.monthly.NSF=Effort.monthly.NSF%>%
                      filter(Same.return%in%unique(Data.monthly.NSF$Same.return))%>%
               dplyr::select(Same.return,hook.hours)%>%
                      group_by(Same.return)%>%
                      summarise(hook.hours=max(hook.hours))%>%
                      rename(Effort=hook.hours)
 
Effort.daily.NSF=Effort.daily.NSF%>%
                      filter(Same.return.SNo%in%unique(Data.daily.NSF$Same.return.SNo))%>%
               dplyr::select(Same.return.SNo,hook.hours)%>%
                      group_by(Same.return.SNo)%>%
                      summarise(hook.hours=max(hook.hours))%>%
                      rename(Effort=hook.hours)


# 4.2 Analyse data  
#note: this function puts data in wide form (one observation per record,
#       adds effort and performs GAM)
fn.ann=function(ktch,efrt,Joint)
{
  #add effort   
  d=left_join(ktch,efrt,by=Joint)%>%
    filter(!is.na(Effort))
  
  #wide form
  other.stuff=d[!duplicated(d[,match(Joint,names(d))]),]%>%
                dplyr::select(-c(SNAME,LIVEWT.c))
   other.stuff=other.stuff[,-match("SPECIES",names(other.stuff))] 
 
    
  ktch.wide=d%>%
        group_by(.dots=c(Joint,'SPECIES'))%>%
        summarise(LIVEWT.c=sum(LIVEWT.c))%>%
        dplyr::select(Joint,SPECIES,LIVEWT.c)%>%
               spread(SPECIES,LIVEWT.c,fill=0)
  Sp=colnames(ktch.wide)[-1]
  
  ktch.wide=ktch.wide%>%
              left_join(other.stuff,by=Joint)
  
  #run GAM on presence/absence with effort as offset
  Store=vector('list',length(Sp))
  names(Store)=Sp
  for(s in 1:length(Sp))
  {
    Drp=Sp[-s]
    dd=ktch.wide[,-match(Drp,names(ktch.wide))]%>%
            rename(LIVEWT.c=!!Sp[s])%>%
            mutate(LIVEWT.c=ifelse(LIVEWT.c>0,1,0),
                   log.effort=log(Effort))
    
    #binominal gam
    mod <-gam(LIVEWT.c~s(MONTH,bs='cc', k = 12)+s(long.corner,lat.corner)+offset(log.effort),
              data=dd,family=binomial,method="REML")
    Store[[s]]=list(mod=mod,dat=dd)
  }
  return(Store)
}
system.time({
  Res.monthly=fn.ann(ktch=Data.monthly,
                     efrt=Effort.monthly,
                     Joint='Same.return')
  Res.daily=fn.ann(ktch=Data.daily,
                   efrt=Effort.daily,
                   Joint='Same.return.SNo')
  Res.monthly.NSF=fn.ann(ktch=Data.monthly.NSF,
                         efrt=Effort.monthly.NSF,
                         Joint='Same.return')
  Res.daily.NSF=fn.ann(ktch=Data.daily.NSF,
                       efrt=Effort.daily.NSF,
                       Joint='Same.return.SNo')
})    #takes 525 secs


# 4.3 Predict spatial occurrence
fn.pred=function(dd)
{
  Store=vector('list',length(dd))
  names(Store)=names(dd)
  for(i in 1:length(dd))
  {
    mod=dd[[i]]$mod
    dat=dd[[i]]$dat
    Mn=sort(table(dat$MONTH))
    new.gam=data.frame(log.effort=mean(dat$log.effort),
                       MONTH=as.numeric(names(Mn[length(Mn)])),
                       long.corner=dat$long.corner,
                       lat.corner=dat$lat.corner)%>%
      mutate(long.lat=paste(long.corner,lat.corner))%>%
      distinct(long.lat,.keep_all = T)
    new.gam$Prob=predict(mod,newdata=new.gam,type='response')
    Store[[i]]=new.gam
  }
  return(Store)
  
}
system.time({
  Pred.monthly=fn.pred(dd=Res.monthly)
  Pred.daily=fn.pred(dd=Res.daily)
  Pred.monthly.NSF=fn.pred(dd=Res.monthly.NSF)
  Pred.daily.NSF=fn.pred(dd=Res.daily.NSF)
})


# 4.4 Calculate spatial overlap


# 5. Outputs section ---------------------------------------------------------

# 5.1 Export spatial overlap table


# 5.2 Map spatial overlap table

