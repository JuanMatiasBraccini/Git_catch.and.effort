# Header ---------------------------------------------------------
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)
library(tidyverse)

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


  #Spatial closures and marine parks

# 3. Parameters section ---------------------------------------------------------
Min.obs=100

# 4. Procedure section ---------------------------------------------------------

Fin.yr.mon=paste(1975:2005,substr(1976:2006,3,4),sep='-')

# 4.1 Manipulate data 
  #Catch
Data.monthly=Data.monthly%>%
                filter(SPECIES<=24900 & !SPECIES==22999 &    #sharks only
                         METHOD=="GN" & FINYEAR%in%Fin.yr.mon)%>%  
                mutate(SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  #Keep species with a minimum of Min.obs
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return,VESSEL,Bioregion,zone,BLOCKX,LONG,LAT,
                       SNAME,SPECIES,Sch.or.DogS,FINYEAR,YEAR.c,MONTH,LIVEWT.c)
            

Data.daily=Data.daily%>%
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="GN")%>%
                mutate(SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return.SNo,Same.return,VESSEL,Bioregion,zone,block10,
                       BLOCKX,LONG,LongMin,LAT,LatMin,SHOTS,
                       SNAME,SPECIES,FINYEAR,YEAR.c,MONTH,LIVEWT.c)
                
Data.monthly.NSF=Data.monthly.NSF%>%
                filter(SPECIES<=24900 & !SPECIES==22999 
                       & METHOD=="LL" & FINYEAR%in%Fin.yr.mon)%>%
                mutate(SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
        dplyr::select(Same.return,VESSEL,Bioregion,zone,BLOCKX,LONG,LAT,
                      SNAME,SPECIES,FINYEAR,YEAR.c,MONTH,LIVEWT.c)
  
Data.daily.NSF=Data.daily.NSF%>%
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="LL")%>%
                mutate(SNAME=tolower(SNAME),
                       SNAME=ifelse(SPECIES==18003,'shark, bronze whaler',
                             ifelse(SPECIES==18023,'shark, spinner',SNAME)))%>%
                group_by(SPECIES)%>%  
                mutate(count = n())%>%
                filter(count>=Min.obs)%>%
         dplyr::select(Same.return.SNo,Same.return,VESSEL,Bioregion,zone,block10,
                       BLOCKX,LONG,LongMin,LAT,LatMin,SHOTS,
                       SNAME,SPECIES,FINYEAR,YEAR.c,MONTH,LIVEWT.c)

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
  for(s in 1:length(Sp))
  {
    Drp=Sp[-s]
    dd=ktch.wide[,-match(Drp,names(ktch.wide))]%>%
            rename(LIVEWT.c=!!Sp[s])%>%
            mutate(LIVEWT.c=ifelse(LIVEWT.c>0,1,0),
                   log.effort=log(Effort))
    
    #ACA create binominal gam with effort as offset....
    mod <-gam(Formula.gam,data=dd,method="REML")
  }
}
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




