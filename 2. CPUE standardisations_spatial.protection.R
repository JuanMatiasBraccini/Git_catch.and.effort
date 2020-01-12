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

# Manipulate data to get into right format
Data.monthly=Data.monthly%>%
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="GN")%>%   #sharks only
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
                filter(SPECIES<=24900 & !SPECIES==22999 & METHOD=="LL")%>%
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
  
#Put data in spread form (one observation per record) and add effort (unique and sum)
#ACA


