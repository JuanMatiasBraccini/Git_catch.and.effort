# MISSING: Spatial closures and marine parks shape files and blockx within closures
#           calculate overlap. 
#         Note that Figure 3 of Brikmanis paper has the MAPs
#           used in their study, make sure these are accounted for


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
######dummy
Closed.Blks=c(20000:25000,31000:33000)
BLKS=sort(c(unique(Data.monthly$BLOCKX),
            unique(Data.monthly.NSF$BLOCKX),
            unique(Data.daily$BLOCKX),
            unique(Data.daily.NSF$BLOCKX)))
Closures=rbind(Data.daily%>%dplyr::select(BLOCKX,block10),
              Data.daily.NSF%>%dplyr::select(BLOCKX,block10))%>%
                    distinct(block10,.keep_all = T)
Misn.blk=data.frame(BLOCKX=BLKS[which(!BLKS%in%Closures$BLOCKX)])%>%
  mutate(block10=BLOCKX*10)
Closures=rbind(Closures,Misn.blk)%>%
              mutate(Closed=ifelse(BLOCKX%in%Closed.Blks,"YES","NO"))

###


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
                       BLOCKX=dat$BLOCKX,
                       long.corner=dat$long.corner,
                       lat.corner=dat$lat.corner)%>%
      mutate(long.lat=paste(long.corner,lat.corner))
    if("block10" %in% colnames(dat)) new.gam$block10=dat$block10
    new.gam=new.gam%>%distinct(long.lat,.keep_all = T)
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
})  #takes 8 secs


# 4.4 Convert Prob of occurrence into categories
fn.cat=function(dd)
{
  for( i in 1:length(dd))
  {
    dd[[i]]=dd[[i]]%>%mutate(Suitability=
                               ifelse(Prob< .3,'Low',
                               ifelse(Prob>= .3 & Prob<.6,'Suitable',
                               ifelse(Prob>=.6,'High',NA))),
                             Suitability=factor(Suitability,
                                                levels=c('Low',
                                                         'Suitable',
                                                         'High')))
  }
  return(dd)
}
system.time({
  Pred.monthly=fn.cat(dd=Pred.monthly)
  Pred.daily=fn.cat(dd=Pred.daily)
  Pred.monthly.NSF=fn.cat(dd=Pred.monthly.NSF)
  Pred.daily.NSF=fn.cat(dd=Pred.daily.NSF)
})  #takes 1 secs


# 4.5 Calculate spatial overlap   
#note:  calculate the number of grid cells within
#       each IUCN zone for each habitat category (low, suitable, high)
fn.overlap=function(dd)
{
  Overlap=vector('list',length(dd))
  names(Overlap)=names(dd)
  for( i in 1:length(dd))
  {
    Clsd=Closures%>%filter(BLOCKX%in%unique(dd[[i]]$BLOCKX))
    if("block10" %in% colnames(dd[[i]])) 
    {
      dd[[i]]=dd[[i]]%>%left_join(Clsd,by=c("block10"))
    }else
    {
      dd[[i]]=dd[[i]]%>%left_join(Clsd,by=c("BLOCKX"))
    }
    
    Overlap[[i]]=with(dd[[i]],table(Closed,Suitability,useNA = 'ifany'))
  }
  return(Overlap)
}

Overlap.monthly=fn.overlap(dd=Pred.monthly)
Overlap.daily=fn.overlap(dd=Pred.daily)
Overlap.monthly.NSF=fn.overlap(dd=Pred.monthly.NSF)
Overlap.daily.NSF=fn.overlap(dd=Pred.daily.NSF)

#Deje aca: how to output this???

# 5. Outputs section ---------------------------------------------------------
setwd('C:\\Matias\\Analyses\\Catch and effort\\Outputs\\Spatial protection')

# 5.1 Summary of available data
Tab.sp.name= Data.monthly%>%
           dplyr::select(SNAME,SPECIES)%>%
                  bind_rows(Data.daily%>%
           dplyr::select(SNAME,SPECIES))%>%
                  bind_rows(Data.monthly.NSF%>%
           dplyr::select(SNAME,SPECIES))%>%
                  bind_rows(Data.daily.NSF%>%
           dplyr::select(SNAME,SPECIES))%>%
                  distinct(SPECIES,.keep_all = T)

fn.tab1=function(dd)
{
  Tbl=as.data.frame(table(dd$SPECIES))
  colnames(Tbl)=c("SPECIES","N")
  Tbl$SPECIES=as.character(Tbl$SPECIES)
  return(Tbl)
}
N.monthly=fn.tab1(dd=Data.monthly)
N.daily=fn.tab1(dd=Data.daily)
N.monthly.NSF=fn.tab1(dd=Data.monthly.NSF)
N.daily.NSF=fn.tab1(dd=Data.daily.NSF) 
N.tot=full_join(N.monthly,N.daily,by="SPECIES")%>%
  full_join(N.monthly.NSF,by="SPECIES")%>%
  full_join(N.daily.NSF,by="SPECIES")%>%
  replace(is.na(.), 0)%>%
  mutate(Number.of.occurrences=N.x+N.y+N.x.x+N.y.y,
         SPECIES=as.numeric(SPECIES))%>%
  dplyr::select(SPECIES,Number.of.occurrences)%>%
  full_join(Tab.sp.name,by="SPECIES")%>%
  select(SNAME,Number.of.occurrences)

write.csv(N.tot,'Table1_Number of records per species.csv',row.names = F)

Total.records=length(unique(Data.monthly$Same.return))+
              length(unique(Data.daily$Same.return.SNo))+
              length(unique(Data.monthly.NSF$Same.return))+
              length(unique(Data.daily.NSF$Same.return.SNo))

write.csv(Total.records,'Table1_Total records.csv',row.names = F)

# 5.2 Export spatial overlap table


# 5.3 Map spatial overlap table

