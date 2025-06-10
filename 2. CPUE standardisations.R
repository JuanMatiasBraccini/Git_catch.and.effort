#--------- CPUE STANDARDISATIONS OF DIFFERENT DATASETS ---------#

#NOTE:  This script standardises the catch and effort data for the main TDGDLF shark species,

#       Irrelevant. Model Term not selected:
#             If want to update SOI and Mean Freo Sealevel each year, run "Get.SOI.Freo.R"  in C:\Matias\Data\Oceanography
#             If want to update Temperature, run "SST.r"   in C:\Matias\Data\Oceanography


# HEADER -----------------------------------------------------------------------
rm(list=ls(all=TRUE))

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


#library(glmmADMB)
library(lattice)
library(bbmle) #for AIC comparisons
library(tweedie)
library(pscl) #zero inflated GLMs
library(ggplot2)
library(lme4) #mixed models
#library(nlme) 
#detach("package:lme4")  # detach lme4 (not compatible with nlme)
library(pvclust) #cluster anlyses
library(cluster) 
library(fpc)
library(coefplot)
#install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source")
#library(coefplot2)
require(statmod) # Provides tweedie family functions
library(VGAM) #zero truncated models
library(Hmisc)#for error bars
library(plotrix)
library(qpcR)  #AKAIKE stuff
library(zoo)    		#needed for filling in NAs
library(lsmeans)
library(mvtnorm)      #for multivariate normal pdf
library(lunar)     #moon phases
library(MASS)
library(stringr)
library(dplyr)
library(tidyr)
library(corrplot)
library(cluster)
library(factoextra) #for plotting
#library(cede)     #Malcolm Haddon's
library(gridExtra)
#library(glmulti)  #model selection
library(fitdistrplus)  #select distribution
library(coefplot)  #visualise coefficients
library(emmeans)  #for model predictions
library(doParallel)
library(tibble)
library(cluster)
library(factoextra) #for plotting
library(mgcv)
library(data.table)
library(PBSmapping)
library(ggpubr)
library(broom)  #nice display of model fit
library(mgcViz)
library(grid)
library(tictoc)
library(ggcorrplot)
library(correlation)
library(ggrepel)
library(rlist)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240,dplyr.summarise.inform = FALSE)   


setwd(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other"))
source("Compare.error.structure.R")
source("Deviance.explained.R")
source("Sorting.objects.R")
#source("MS.Office.outputs.R")
setwd(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics"))
source("fn.fig.R")
source("Nominal_cpue_functions.R")
source(handl_OneDrive("R\\Sort ls by size.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))

source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))  #my themes

source(handl_OneDrive('Analyses/Catch and effort/Git_catch.and.effort/CPUE_Auxiliary functions.R'))

# DATA SECTION -----------------------------------------------------------------------
setwd(handl_OneDrive('Analyses/Data_outs'))
Data.daily.original=fread("Data.daily.original.csv",data.table=FALSE)
Data.monthly.GN=fread("Data.monthly.GN.csv",data.table=FALSE) 
Data.daily.GN=fread("Data.daily.GN.csv",data.table=FALSE)
Effort.daily=fread("Effort.daily.csv",data.table=FALSE)
Effort.monthly=fread("Effort.monthly.csv",data.table=FALSE)
Mesh.monthly=fread("Mesh.monthly.csv",data.table=FALSE)
Mesh.size=fread("Mesh.size.csv",data.table=FALSE)
Mangmnt.TDGDLF=fread("Mangmnt.TDGDLF.csv",data.table=FALSE)


lst <- strsplit(Data.daily.GN$Same.return.SNo, "\\s+")
Data.daily.GN$SNo <- sapply(lst ,'[', 1)
Data.daily.GN$DSNo <- sapply(lst, '[', 2)
Data.daily.GN$TSNo <- sapply(lst, '[', 3)
rm(lst)

#Block10 locations
BlOCK_10=fread(handl_OneDrive("Data/Mapping/Blocks_10NM.csv"),data.table=FALSE)
names(BlOCK_10)=c("block10","LAT","LONG")
Metro_BlOCK_10=subset(BlOCK_10, LAT>(-33) & LAT<=(-31) & LONG<116)

#Southern Oscillation Index (index of La Niña, El Niño events)
#note: negative values of SOI below −7 indicate El Niño episodes
#      positive values of SOI above +7 are typical of La Niña episode
SOI=fread(handl_OneDrive("Data/Oceanography/SOI.csv"),data.table=FALSE)

#Mean Freo sea level (index of Leeuwin Current)
Freo=fread(handl_OneDrive("Data/Oceanography/Freo_mean_sea_level.csv"),data.table=FALSE)  

#SST
SST=fread(handl_OneDrive("Data/Oceanography/SST.csv"),data.table=FALSE) 

#Fishable areas       
Depth.range="species_specific"
#Depth.range=200
if(Depth.range==200)
{
  Whis.fishArea=fread(handl_OneDrive("Data/Catch and Effort/FishableArea/BLOCKX_whiskery_200.csv"),data.table=FALSE)
  Gum.fishArea=fread(handl_OneDrive("Data/Catch and Effort/FishableArea/BLOCKX_gummy_200.csv"),data.table=FALSE)
  Dusky.fishArea=fread(handl_OneDrive("Data/Catch and Effort/FishableArea/BLOCKX_dusky_200.csv"),data.table=FALSE)
  Sand.fishArea=fread(handl_OneDrive("Data/Catch and Effort/FishableArea/BLOCKX_sandbar_200.csv"),data.table=FALSE)  
}
if(Depth.range=="species_specific")
{
  Grab.Area=handl_OneDrive("Data/Catch and Effort/FishableArea/")
  Whis.fishArea=fread(paste(Grab.Area,"BLOCKX_whiskery_30_70.csv",sep=""),data.table=FALSE)
  Gum.fishArea=fread(paste(Grab.Area,"BLOCKX_gummy_les_70.csv",sep=""),data.table=FALSE)
  Dusky.fishArea=fread(paste(Grab.Area,"BLOCKX_dusky_less_60.csv",sep=""),data.table=FALSE)
  Sand.fishArea=fread(paste(Grab.Area,"BLOCKX_sandbar_30_120.csv",sep=""),data.table=FALSE)  
  
  Whis.fishArea_b10=fread(paste(Grab.Area,"block10.whiskery.csv",sep=""),data.table=FALSE)
  Gum.fishArea_b10=fread(paste(Grab.Area,"block10.gummy.csv",sep=""),data.table=FALSE)
  Dusky.fishArea_b10=fread(paste(Grab.Area,"block10.dusky.csv",sep=""),data.table=FALSE)
  Sand.fishArea_b10=fread(paste(Grab.Area,"block10.sandbar.csv",sep=""),data.table=FALSE)  
}


#Vessel characteristics from DoF's annual survey and questionnaires
Vessel.charac=read.csv(handl_OneDrive("Data/Fishing power/Vessel.charac.for_TDGDLF.cpue.stand.csv"))
#extracted from source("...Git_catch_and_effort/3.Vessel_characteristics.R")

#vessel monthly CFL
Monthly_CLF_vessel=read.csv(handl_OneDrive("Data/Catch and Effort/Monthly_CLF_vessel.csv"))


# PARAMETERS SECTION -----------------------------------------------------------------------

#Control if monthly cpue is exported for stock assessments
xport_monthly=FALSE # SQL server overwrite 'Reporter' column so now all records id as 'good reporters' for monthly records

#1.1 Control what parts of script are activated

#1.1.1 Data controls

#Control if doing separated analysis of monthly and daily records
Separate.monthly.daily="YES"

#Control how to aggregate daily records
#Use.Date="YES"    #Rory's approach (aggregating by DATE). this is used for getting the monthly agg.
Use.Date="NO"     # aggregating Daily records by Same.return.SNo (i.e. SNo, DsNo and TSNo). 
#                     This is more appropriate as some fishers are fishing and reporting more than 1 shot 
#                     per day so if using date, then only the max of the effort of
#                     that day is used, underestimating effort

#Control if removing closures' data or not
Remove.closure="NO"

#control if showing example of folly index
Show.folly.eg="NO"
Show.variability.cpue.eg="NO"

#control how to remove blocks with low catch
Remove.blk.by="blk_only"   #remove only blocks
#Remove.blk.by="yr_blk"    #remove by yr_block combos

#control if removing or reallocating border blocks
BOUND.BLK="REALLOCATE"
#BOUND.BLK="REMOVE"


#1.1.2 Procedure controls 

#Define is manually setting core areas based on Rory's
change.core.manually=FALSE

#Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.jpeg="NO"
Do.tiff="YES"

#Control type of model run
#Model.run="First"    # for first time doing standardisation. This allows selection of 
#     best model and sensitivity tests
Model.run="Standard"  # for running standardisations in subsequent years

#Control if calculating efficiency creep from vessel survey data
get.efficiency.creep=FALSE #cannot match ves vars with catch return (need vessel-boath, only vessel available)

#Define is testing model structure (i.e. selection of model terms)
Def.mod.Str="NO"
do.monthly.def.str=TRUE

#Control if excluding firt year of monthly data for gummy and sandbar
Exclude.yr.gummy=TRUE  #very small sample size
Exclude.yr.sandbar=TRUE

#Control what approach used for standardisation
Use.Tweedie=TRUE
Use.Delta=FALSE
Use.Qualif.level=FALSE

#Control if doing spatial standardisation
do.spatial.cpiui=TRUE

#Control if doing Stephens and McCall
do_Stephens_McCall="YES"

#Control if doing cluster analysis of daily data
do_cluster="NO"

#Control if doing PCA analysis of daily data                    
do_pca="NO"   

#Control if doing sensitivity tests
if(Model.run=="First") do.sensitivity="YES"
if(Model.run=="Standard") do.sensitivity="NO"

#Control if comparing glm and gam spatial predictions
compare.glm.gam_spatial="NO"


#1.1.3 Reporting controls

#Control if doing influence plots (Bentley et al 2012)
if(Model.run=="First")  do.influence="YES"
if(Model.run=="Standard")  do.influence="NO"


#control color of plots
#do.colors="YES"
do.colors="NO"  #grey scale
if(do.colors=="YES") what.color="cols"
if(do.colors=="NO") what.color="black"

#Control if doing exploratory analyses
if(Model.run=="First") do.Exploratory="YES" else
  do.Exploratory="NO" 

#Control if producing cpue paper figure
plot.cpue.paper.figures="NO"

#Use area weight in annual index
use.blok.area='NO'


#1.2 Criteria for keeping species for analysis
N.keep=5      #in years
Min.kg=100   #in kg
Min.annual.prop.zero=0.2  #minimum annual proportion of zero records to be selected for analysis
core.per=90

#Maximum possible effort
Net.max=8500    #Rory pers comm
Max.km.gn.h.monthly=24*20*Net.max/1000
Max.km.gn.h.daily=24*Net.max/1000

#Species definitions
Shark.species=c(5001:24900)  #note: rays are reported as a mixed so standardisation for Eagle rays (39001) is not possible
Indicator.sp=c(17001,17003,18003,18007)
Greynurse.protection='1999-00'
TARGETS.name=c("SHARK, WHISKERY","SHARK, GUMMY","SHARK, DUSKY WHALER","SHARK, THICKSKIN (SANDBAR)")
TARGETS=list(17003,17001,18003,18007)
names(TARGETS)=TARGETS.name
N.species=length(TARGETS)
SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
SPvec=c(17003,17001,18003,18007)
names(SPvec)=SPECIES.vec

reset.hammerheads.to.smooth.hh=FALSE  #assume that all hammerheads are smooth HH east of 120 E

#Zn1-Zn2 boundary blocks
Boundary.Blks=c(3416,3516,3616)  

#Criteria for selecting vessels and blocks for standardisation
MIN.ktch=100   #Minimum annual catch of target species (in kgs)
MIN.records.yr=3 #Minimum number of records of target species per vessel per year for at least the Threshold.n.yrs
MIN.obs.BLK=3   # Block. minimum number of years with observations for each block
MIN.obs.BLK.sens=0  
MIN.obs.BL10K=3 # Block. minimum number of observation per block10
Threshold.n.yrs=3   #Vessel. minimum number of years reporting the species
Threshold.n.yrs.monthly=3 
Threshold.n.yrs.daily=3 
Threshold.n.yrs.sens=0  #sensitivity
Threshold.n.vessls.per.yr=3  #Vessel. keep years with a least 3 different vessels. Not used


#qualification levels   #not used
QL_expl_ktch_prop=.9   #proportion of explained annual catch for selected record
PRP.MLCLM=0.1          #proportion of catch of target


#Qualification levels minimum number of years with positive records
Min.Vess.yr=5 #monthly
Min.Vess.yr.d=5 #daily

#Cpue correction for assumed increase in fishing power prior to 1995
# Rory McAuley: 2% annual upto 1994-95, then constant...consistent with 2.4% reported by 
#               Palomares & Pauly 2019
Fish.Pow=.02


#Minimum weights
#McAuley unpublished (TL to FL)
a.w=1.0044;b.w=13.171
a.g=1.0837;b.g=4.6424
a.d=1.1849;b.d=2.9835
a.s=1.2185;b.s=0.2133

bwt.g=0.0000004623   #DOF unpublished
bwt.w=0.0000027500   #McAuley and Simpfendorfer, 2003
bwt.d=0.0000034694   #McAuley and Simpfendorfer, 2003
bwt.s=0.0000021684   #DOF unpublished
awt.g=3.47701
awt.w=3.08059
awt.d=3.10038
awt.s=3.20688


#Min observed FL in catch (observers data)
Min.FL.w=88
Min.FL.g=43
Min.FL.d=55
Min.FL.s=47

#Size at birth (FL,cm)
Size.birth.w=25
Size.birth.g=33.5
Size.birth.d=75.3
Size.birth.s=42.5


aa=function(FL,a,b) TL=FL*a+b
ww=function(TL,a,b) TW=b*TL^a

#min weight (kg) in catch
TW.w=ww(aa(Min.FL.w,a.w,b.w),awt.w,bwt.w)
TW.g=ww(aa(Min.FL.g,a.g,b.g),awt.g,bwt.g)
TW.d=ww(aa(Min.FL.d,a.d,b.d),awt.d,bwt.d)
TW.s=ww(aa(Min.FL.s,a.s,b.s),awt.s,bwt.s)

#min weight (kg) in population
TW.w.birth=ww(aa(Size.birth.w,a.w,b.w),awt.w,bwt.w)
TW.g.birth=ww(aa(Size.birth.g,a.g,b.g),awt.g,bwt.g)
TW.d.birth=ww(aa(Size.birth.d,a.d,b.d),awt.d,bwt.d)
TW.s.birth=ww(aa(Size.birth.s,a.s,b.s),awt.s,bwt.s)

Min.weight=c(.8,.8,.8,.8)
names(Min.weight)=SPECIES.vec

#max weights
max.w.whis=ww(160,awt.w,bwt.w)*1.2  #max TL from Last and Stevens plust 20% to allow for max observations
max.w.gum=ww(175,awt.g,bwt.g)*1.2
max.w.dus=ww(365,awt.d,bwt.d)*1.2
max.w.san=ww(240,awt.s,bwt.s)*1.2


Max.weight=c(max.w.whis,max.w.gum,max.w.dus,max.w.san)
names(Max.weight)=SPECIES.vec


#Year catchability changed for whiskery sharks (Simpfendorfer 2000).    NOT USED....
#Q_change="1982-83"    
Q_change="1975-76"    #use all monthly records in one series (this bypasses the use of Q_change)

#5% percent increase in catch and effort prior 1990 
#note: Rory McAuley
Inc.per=1.05  

#McAuly core areas
Dusky.range=c(-28,120)
Sandbar.range=c(-26,118)
Whiskery.range=c(-28,129)
Gummy.range=c(116,129) 



# Remove duplicates and incomplete last year & aggregate weights-----------------------------------------------------------------------
Tab.incomplit=table(Data.daily.GN$FINYEAR,Data.daily.GN$MONTH)
Tab.incomplit[Tab.incomplit>0]=1
Incomplete.year=names(which(rowSums(Tab.incomplit)<12))
if(length(Incomplete.year)>0)
{
  Data.daily.GN=subset(Data.daily.GN,!FINYEAR==Incomplete.year)
  Effort.daily=subset(Effort.daily,!finyear==Incomplete.year)
  Data.monthly.GN=subset(Data.monthly.GN,!FINYEAR==Incomplete.year)
  Effort.monthly=subset(Effort.monthly,!FINYEAR==Incomplete.year)
}


Data.daily.GN=Data.daily.GN%>%
  mutate(Dupli=paste(Same.return.SNo,SPECIES,CONDITN,LIVEWT.c))%>% 
  distinct(Dupli,.keep_all = T)%>%
  group_by(Same.return.SNo,SPECIES)%>%
  mutate(LIVEWT=sum(LIVEWT,na.rm=T),
         LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
  mutate(Dupli=paste(Same.return.SNo,SPECIES))%>%
  distinct(Dupli,.keep_all = T)%>%
  dplyr::select(-c(LIVEWT.orgnl,LIVEWT.reap,Dupli))%>%
  ungroup()
Data.monthly.GN=Data.monthly.GN%>%
  mutate(Dupli=paste(Same.return,SPECIES,CONDITN,LIVEWT.c))%>%
  distinct(Dupli,.keep_all = T)%>%
  group_by(Same.return,SPECIES)%>%
  mutate(LIVEWT=sum(LIVEWT,na.rm=T),
         LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
  mutate(Dupli=paste(Same.return,SPECIES))%>%
  distinct(Dupli,.keep_all = T)%>%
  dplyr::select(-c(LIVEWT.orgnl,LIVEWT.reap,Dupli))%>%
  ungroup()

# BASIC MANIPULATIONS -----------------------------------------------------------------------

#Remove Nil fish caught from daily
a1=Data.daily.GN%>%filter(SPECIES==9998)
Data.daily.GN=Data.daily.GN%>%filter(!SPECIES==9998)
Effort.daily=Effort.daily%>%filter(!Same.return.SNo%in%unique(a1$Same.return.SNo))

#Reset BLOCKX to 4 digits as some have 5 digits
Data.monthly.GN$BLOCKX=as.integer(substr(Data.monthly.GN$BLOCKX,1,4))
Effort.monthly$BLOCKX=as.integer(substr(Effort.monthly$BLOCKX,1,4))
Data.daily.GN$BLOCKX=as.integer(substr(Data.daily.GN$BLOCKX,1,4))
Effort.daily$blockx=as.integer(substr(Effort.daily$blockx,1,4))


#...Effort 
# -- Daily                
Effort.daily$LAT=-abs(Effort.daily$LAT)
Block.lat.long=Effort.daily%>%distinct(blockx,.keep_all = TRUE)%>%dplyr::select(blockx,LAT,LONG)

#km gn days
Eff.daily.c.daily=aggregate(Km.Gillnet.Days.c~Same.return.SNo+vessel+finyear+month+blockx+block10,
                            data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)    
Eff.daily.daily=aggregate(Km.Gillnet.Days.inv~Same.return.SNo+vessel+finyear+month+blockx+block10,
                          data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)

#km gn hours 
Eff.daily.hour.c.daily=aggregate(Km.Gillnet.Hours.c~Same.return.SNo+vessel+finyear+month+blockx+block10,
                                 data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)    
Eff.daily.hour.daily=aggregate(Km.Gillnet.Hours.inv~Same.return.SNo+vessel+finyear+month+blockx+block10,
                               data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)  

#merge into single file
Eff.daily.c.daily=merge(Eff.daily.c.daily,Eff.daily.hour.c.daily,by=c("blockx","block10","finyear","month","vessel","Same.return.SNo"),all=T)
Eff.daily.daily=merge(Eff.daily.daily,Eff.daily.hour.daily,by=c("blockx","block10","finyear","month","vessel","Same.return.SNo"),all=T)
Eff.daily.c.daily=merge(Eff.daily.c.daily,Eff.daily.daily,by=c("blockx","block10","finyear","month","vessel","Same.return.SNo"),all.x=T)
Eff.daily.c.daily=merge(Eff.daily.c.daily,Block.lat.long,by=c("blockx"),all.x=T)  


# -- Monthly 
Monthly=subset(Effort.monthly,NETLEN.c>100 & METHOD=="GN")
Block.lat.long=Effort.monthly%>%distinct(BLOCKX,.keep_all = TRUE)%>%dplyr::select(BLOCKX,LAT,LONG)

#km gn days     
Eff.monthly.c=aggregate(Km.Gillnet.Days.c~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)
Eff.monthly=aggregate(Km.Gillnet.Days.inv~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)

#km gn hours                     
Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)
Eff.monthly.hour=aggregate(Km.Gillnet.Hours.inv~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)

#merge into single file
THESE=c("BLOCKX","FINYEAR","MONTH","VESSEL","YEAR.c")
Eff.monthly.c=merge(Eff.monthly.c,Eff.monthly.hour.c,by=THESE,all=T)
Eff.monthly=merge(Eff.monthly,Eff.monthly.hour,by=THESE,all=T)
Eff.monthly.c=merge(Eff.monthly.c,Eff.monthly,by=THESE,all.x=T)
Eff.monthly.c=merge(Eff.monthly.c,Block.lat.long,by=c("BLOCKX"),all.x=T)


#put back effort reporter variable

#daily
a=Effort.daily%>%
  dplyr::select(blockx,block10,Same.return.SNo,vessel,Eff.Reporter)%>%
  mutate(dupli=paste(blockx,block10,Same.return.SNo,vessel))%>%
  distinct(dupli,.keep_all = TRUE)%>%
  dplyr::select(-dupli)
Eff.daily.c.daily=left_join(Eff.daily.c.daily,a,by=c("blockx","block10","Same.return.SNo","vessel"))

#monthly      
a=Effort.monthly%>%
  dplyr::select(BLOCKX,FINYEAR,MONTH,VESSEL,Eff.Reporter)%>%
  mutate(dupli=paste(BLOCKX,FINYEAR,MONTH,VESSEL))%>%
  distinct(dupli,.keep_all = TRUE)%>%
  dplyr::select(-dupli)
Eff.monthly.c=left_join(Eff.monthly.c,a,by=c("BLOCKX","FINYEAR","MONTH","VESSEL"))

rm(a)


#put back effort variables
#monthly
Eff.monthly.c$dummy=with(Eff.monthly.c,paste(BLOCKX,FINYEAR,MONTH,VESSEL,
                                             Km.Gillnet.Hours.c,Km.Gillnet.Days.c))
a=Effort.monthly[,match(c("BLOCKX","FINYEAR","MONTH","VESSEL",
                          "SHOTS.c","NETLEN.c","BDAYS.c","HOURS.c",
                          "Km.Gillnet.Hours.c","Km.Gillnet.Days.c"),names(Effort.monthly))]
a$dummy=paste(a$BLOCKX,a$FINYEAR,a$MONTH,a$VESSEL,a$Km.Gillnet.Hours.c,a$Km.Gillnet.Days.c)
a=a[!duplicated(a$dummy),]
a=subset(a,dummy%in%Eff.monthly.c$dummy,select=c(dummy,SHOTS.c,NETLEN.c,BDAYS.c,HOURS.c))
Eff.monthly.c=merge(Eff.monthly.c,a,by="dummy",all.x=T)
Eff.monthly.c=Eff.monthly.c[,-match("dummy",names(Eff.monthly.c))]

#daily
a=Effort.daily[,match(c("blockx","block10","Same.return.SNo","vessel",
                        "netlen.c","hours.c","bdays.c","shots.c","nlines.c"),names(Effort.daily))]
a$dummy=paste(a$blockx,a$block10,a$Same.return.SNo,a$vessel)
a=a[!duplicated(a$dummy),]
Eff.daily.c.daily$dummy=with(Eff.daily.c.daily,paste(blockx,block10,Same.return.SNo,vessel))
a=subset(a,dummy%in%Eff.daily.c.daily$dummy,select=c(dummy,netlen.c,hours.c,bdays.c,shots.c,nlines.c))
Eff.daily.c.daily=merge(Eff.daily.c.daily,a,by="dummy",all.x=T)
Eff.daily.c.daily=Eff.daily.c.daily[,-match("dummy",names(Eff.daily.c.daily))]
rm(a)


#Merge monthly and daily-aggregated
Eff.monthly.c=Eff.monthly.c[,-match("YEAR.c",names(Eff.monthly.c))]

#monthly
Eff=Eff.monthly.c
Eff$Eff.Reporter=with(Eff,ifelse(is.na(Eff.Reporter),"good",Eff.Reporter))

#daily
Eff.daily=Eff.daily.c.daily
names(Eff.daily)[match(c("blockx","finyear","month","vessel"),names(Eff.daily))]=c("BLOCKX",
                                                                                   "FINYEAR","MONTH","VESSEL")
Eff.daily$Eff.Reporter=with(Eff.daily,ifelse(is.na(Eff.Reporter),"good",Eff.Reporter))
rm(Eff.daily.c.daily)

#Add mesh size    
#Monthly
Mesh.monthly$mesh=Mesh.monthly$MeshSizeHigh
Mesh.monthly$mesh=with(Mesh.monthly,ifelse(mesh==0,NA,mesh))
a=subset(Mesh.monthly,select=c(VESSEL,FINYEAR,MONTH,BLOCKX,METHOD,mesh))
a$dummy=with(a,paste(VESSEL,FINYEAR,MONTH,BLOCKX))
Eff$dummy=with(Eff,paste(VESSEL,FINYEAR,MONTH,BLOCKX))
a=subset(a,dummy%in%Eff$dummy)
a=a[!duplicated(a$dummy),]
a=a[,-match(c("dummy",'METHOD'),names(a))]
Eff=Eff[,-match("dummy",names(Eff))]
Eff=merge(Eff,a,by=c("VESSEL","FINYEAR","MONTH","BLOCKX"),all.x=T)

#Daily
Mesh.size$mesh=Mesh.size$mshigh
a=subset(Mesh.size,select=c(SNo,DSNo,TSNo,mesh))
a$Same.return.SNo=with(a,paste(SNo,DSNo,TSNo))
a=a[!duplicated(a$Same.return.SNo),]
a=subset(a,Same.return.SNo%in%unique(Eff.daily$Same.return.SNo),select=c(Same.return.SNo,mesh))
Eff.daily=merge(Eff.daily,a,by="Same.return.SNo",all.x=T)  
rm(Mesh.size)


#Add depth to daily
a=subset(Data.daily.original,select=c(SNo,DSNo,TSNo,depthMin,depthMax))
a$Same.return.SNo=with(a,paste(SNo,DSNo,TSNo))
a=a[!duplicated(a$Same.return.SNo),]
a$Mean.depth=(a$depthMin+a$depthMax)/2
a=subset(a,Same.return.SNo%in%unique(Eff.daily$Same.return.SNo),select=c(Same.return.SNo,Mean.depth))
Eff.daily=merge(Eff.daily,a,by="Same.return.SNo",all.x=T)
rm(Data.daily.original)


# Define 'hammerhead' records as 'smooth hammerhead' if east of 120 E (Bartes & Braccini 2021) 
if(reset.hammerheads.to.smooth.hh)
{
  Data.daily.GN=Data.daily.GN%>%
    mutate(SPECIES=ifelse(LONG>=120 & SPECIES==19000,19004,SPECIES),
           SNAME=ifelse(SPECIES==19004,"SHARK, SMOOTH HAMMERHEAD",SNAME),
           RSCommonName=ifelse(SPECIES==19004,"Smooth Hammerhead Shark",RSCommonName))
  
  Data.monthly.GN=Data.monthly.GN%>%
    mutate(SPECIES=ifelse(LONG>120 & SPECIES==19000,19004,SPECIES),
           SNAME=ifelse(SPECIES==19004,"SHARK, SMOOTH HAMMERHEAD",SNAME),
           RSCommonName=ifelse(SPECIES==19004,"Smooth Hammerhead Shark",RSCommonName))
}


# Define skates and rays as eagle rays following interviews with fishers for PA project
Data.daily.GN=Data.daily.GN%>%
  mutate(SPECIES=ifelse(LAT<=(-30.5) & SPECIES%in%c(31000,990001),39001,SPECIES),
         SNAME=ifelse(SPECIES==39001,"Eagle ray",SNAME),
         RSCommonName=ifelse(SPECIES==39001,"Eagle ray",RSCommonName))
Data.monthly.GN=Data.monthly.GN%>%
  mutate(SPECIES=ifelse(LAT<=(-30.5) & SPECIES%in%c(31000,990001),39001,SPECIES),
         SNAME=ifelse(SPECIES==39001,"Eagle ray",SNAME),
         RSCommonName=ifelse(SPECIES==39001,"Eagle ray",RSCommonName))



# SELECT SPECIES-------------------------------------------------------------------------
# note: at least 'N.keep' years of catches of at least 'Min.kg') and put into species list
#       Hammerhead species were reported by fishers as 'hammerheads' (code 19000) until 2023 and reapportioned
#       in Fishcube to one of 3 species using observed proportions so don't use reapportioned for species cpue 

#1. First Criteria: annual catch of at least 'Min.kg' for at least 'N.keep' years
A=Data.monthly.GN%>%
  filter(SPECIES%in%Shark.species) %>%
  group_by(SPECIES,FINYEAR) %>%
  summarise(Weight=round(sum(LIVEWT.c))) %>%
  spread(FINYEAR, Weight)%>%
  data.frame()

A=A %>% mutate_at(.vars = names(A)[-match("SPECIES",names(A))], function(x)(ifelse(x>=Min.kg, 1, 0)))
B=as.data.frame(A)

Anm=A$SPECIES
A=rowSums(A[,-1],na.rm=T)
names(A)=Anm
SpiSis=unique(Data.monthly.GN%>%
                mutate(RSCommonName=ifelse(SPECIES==8001,'Greynurse Shark',RSCommonName))%>%
                filter(SPECIES%in%as.numeric(names(A[A>=N.keep])))%>%
                dplyr::select(SPECIES,RSCommonName) %>%
                filter(!RSCommonName==""))
SpiSis=SpiSis[!duplicated(SpiSis$SPECIES),] 
nms=SpiSis$RSCommonName
SpiSis=SpiSis$SPECIES
names(SpiSis)=nms

SpiSis=SpiSis[!SpiSis%in%c(22999)]   

if(do.Exploratory=="YES")
{
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Species selection/Selection.pdf'), width=12)
  
  #Presence/absence
  fn.choose.sp(d=Data.monthly.GN%>%filter(SPECIES%in%Shark.species & !FINYEAR%in%unique(Data.daily.GN$FINYEAR)),
               crit='Same.return',Titl='Monthly')
  
  fn.choose.sp(d=Data.daily.GN%>%filter(SPECIES%in%Shark.species),
               crit='Same.return.SNo',Titl='Daily')
  
  #Pass selection criteria
  B%>%
    gather('Year','Pass',-SPECIES)%>%
    mutate(Year=as.numeric(substr(Year,2,5)),
           Pass=ifelse(!Pass==1,NA,Pass),
           Col=ifelse(SPECIES%in%SpiSis,"full","empty"))%>%
    left_join(Data.monthly.GN%>%distinct(SPECIES,SNAME),by='SPECIES')%>%
    ggplot(aes(Year,Pass),colour=Col)+
    geom_point(aes(bg=factor(Col)),shape=21)+
    facet_wrap(~SNAME)+
    ylab("Pass data selection criteria")+theme_PA(strx.siz=8)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    scale_fill_manual(breaks = c("full", "empty"),
                      values=c("green", "transparent"))+
    ggtitle(paste("Criteria:",N.keep,"years with at least",Min.kg,"kg per year"))
  
  #Delta log nominal cpue
  cpue=Data.monthly.GN%>%
            filter(SPECIES%in%SpiSis & !FINYEAR%in%unique(Data.daily.GN$FINYEAR))%>%
            dplyr::select(Same.return,FINYEAR,SNAME,LIVEWT.c)%>%
            spread(SNAME,LIVEWT.c,fill=0)%>%
            left_join(Effort.monthly%>%distinct(Same.return,Km.Gillnet.Hours.c),by='Same.return')%>%
            gather(SNAME,LIVEWT.c,-c(Same.return,FINYEAR,Km.Gillnet.Hours.c))%>%
            mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
            rename(season=FINYEAR)
  cpue.sp=unique(cpue$SNAME)
  DLnMean_all=vector('list',length(cpue.sp))
  for(i in 1:length(cpue.sp))DLnMean_all[[i]]=fn.delta.log(d=cpue%>%
                                                         filter(SNAME==cpue.sp[i])%>%   
                                                         group_by(season,Same.return)%>%
                                                         summarise(cpue=mean(cpue))%>%
                                                         ungroup())%>%
                                                        mutate(SNAME=cpue.sp[i])
  p=do.call(rbind,DLnMean_all)%>%
                mutate(yr=as.numeric(substr(year,1,4)))%>%
                ggplot(aes(yr,mean))+
    geom_point()+
    geom_line(linetype='dotted')+
    geom_errorbar(aes(ymin=lowCL,ymax=uppCL),alpha=0.5)+
    theme_PA(strx.siz=8)+theme(legend.position = 'top')+ylab('Delta lognormal (mean +/= 95%CI)')+
    facet_wrap(~SNAME,scales='free_y')
  print(p)

  cpue=Data.daily.GN%>%
    filter(SPECIES%in%SpiSis)%>%
    dplyr::select(Same.return.SNo,FINYEAR,SNAME,LIVEWT.c)%>%
    spread(SNAME,LIVEWT.c,fill=0)%>%
    left_join(Effort.daily%>%distinct(Same.return.SNo,Km.Gillnet.Hours.c),by='Same.return.SNo')%>%
    gather(SNAME,LIVEWT.c,-c(Same.return.SNo,FINYEAR,Km.Gillnet.Hours.c))%>%
    mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
    rename(season=FINYEAR)
  cpue.sp=unique(cpue$SNAME)
  DLnMean_all=vector('list',length(cpue.sp))
  for(i in 1:length(cpue.sp))DLnMean_all[[i]]=fn.delta.log(d=cpue%>%
                                                             filter(SNAME==cpue.sp[i])%>%   
                                                             group_by(season,Same.return.SNo)%>%
                                                             summarise(cpue=mean(cpue))%>%
                                                             ungroup())%>%
                                                      mutate(SNAME=cpue.sp[i])
  p=do.call(rbind,DLnMean_all)%>%
    mutate(yr=as.numeric(substr(year,1,4)))%>%
    ggplot(aes(yr,mean))+
    geom_point()+
    geom_line(linetype='dotted')+
    geom_errorbar(aes(ymin=lowCL,ymax=uppCL),alpha=0.5)+
    theme_PA(strx.siz=8)+theme(legend.position = 'top')+ylab('Delta lognormal (mean +/= 95%CI)')+
    facet_wrap(~SNAME,scales='free_y')
  print(p)
  
  dev.off()
}


#2. Second criteria. Minimum proportion of zero catch records
#Get proportion of 0 catch by year
Prop.ktch=vector('list',length(SpiSis))
names(Prop.ktch)=SpiSis
Prop.ktch.daily=Prop.ktch

d=Data.monthly.GN%>%
  dplyr::select(Same.return,FINYEAR,BLOCKX,SPECIES,LIVEWT.c)%>%
  group_by(Same.return,FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  spread(SPECIES,LIVEWT.c,fill=0)

d.daily=Data.daily.GN%>%
  dplyr::select(Same.return.SNo,FINYEAR,BLOCKX,SPECIES,LIVEWT.c)%>%
  group_by(Same.return.SNo,FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  spread(SPECIES,LIVEWT.c,fill=0)

for(s in 1:length(SpiSis))
{
  if(SpiSis[s]%in%colnames(d))
  {
    Bi=d[,c("Same.return","FINYEAR","BLOCKX",SpiSis[s])]%>%
      rename(Catch.Target=!!as.symbol(SpiSis[s]))
    Ag.blk=Bi%>%
      group_by(BLOCKX)%>%
      summarise(tot=sum(Catch.Target))%>%
      arrange(-tot)%>%
      mutate(cumsum=cumsum(tot),
             percentile=cumsum/sum(tot))%>%
      filter(percentile<=core.per/100)%>%
      pull(BLOCKX)
    Bi=Bi%>%filter(BLOCKX%in%Ag.blk)%>%
      mutate(catch.pos=as.numeric(Catch.Target>0))
    TAB=table(Bi$catch.pos,Bi$FINYEAR)
    Prop.ktch[[s]]=round(TAB[2,]/colSums(TAB),2)
  }
  
  if(SpiSis[s]%in%colnames(d.daily))
  {
    Bi=d.daily[,c("Same.return.SNo","FINYEAR","BLOCKX",SpiSis[s])]%>%
      rename(Catch.Target=!!as.symbol(SpiSis[s]))
    Ag.blk=Bi%>%
      group_by(BLOCKX)%>%
      summarise(tot=sum(Catch.Target))%>%
      arrange(-tot)%>%
      mutate(cumsum=cumsum(tot),
             percentile=cumsum/sum(tot))%>%
      filter(percentile<=core.per/100)%>%
      pull(BLOCKX)
    Bi=Bi%>%filter(BLOCKX%in%Ag.blk)%>%
      mutate(catch.pos=as.numeric(Catch.Target>0))
    TAB=table(Bi$catch.pos,Bi$FINYEAR)
    Prop.ktch.daily[[s]]=round(TAB[2,]/colSums(TAB),2)
  }
}
rm(d,d.daily)

Prop.ktch=do.call(rbind,Prop.ktch)
Prop.ktch.daily=do.call(rbind,Prop.ktch.daily)

write.csv(Prop.ktch,handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/Prop.records.with.catch.by.year.csv"))
write.csv(Prop.ktch.daily,handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/Prop.records.with.catch.daily.by.year.csv"))

#keep species with at least 10 of years with Min.annual.prop.zero
Prop.ktch[Prop.ktch>=Min.annual.prop.zero]=1
Prop.ktch[Prop.ktch<Min.annual.prop.zero]=0

Prop.ktch.daily[Prop.ktch.daily>=Min.annual.prop.zero]=1
Prop.ktch.daily[Prop.ktch.daily<Min.annual.prop.zero]=0

First.year.catch=apply(Prop.ktch,1, function(x) head(x[x!=0],1))
if(Min.annual.prop.zero<0.2) names(First.year.catch$'18007')='1986-87'  #cannot estimate coef for 1985
First.year.catch.daily=apply(Prop.ktch.daily,1, function(x) head(x[x!=0],1))

Annual.year.sp=names(which(rowSums(Prop.ktch)>=10))  
Annual.year.sp.daily=names(which(rowSums(Prop.ktch.daily)>=10))

SpiSis=subset(SpiSis,SpiSis%in%unique(c(Annual.year.sp,Annual.year.sp.daily)))

SP.list=as.list(SpiSis)

#remove these species; not enough positive record data to estimate glm coefficients    
SP.list=SP.list[-match(c("School Shark"),names(SP.list))]

#remove these as different species mixed up and not fishing core area for sawsharks and hammerheads
SP.list=SP.list[-match(c("Wobbegong","Common Sawshark","Hammerhead Sharks"),names(SP.list))]  

SpiSis=SpiSis[match(names(SP.list),names(SpiSis))]

First.year.catch=First.year.catch[match(unlist(SP.list),names(First.year.catch))]
First.year.catch.daily=First.year.catch.daily[match(unlist(SP.list),names(First.year.catch.daily))]

First.year.catch=First.year.catch[!is.na(names(First.year.catch))]
First.year.catch.daily=First.year.catch.daily[!is.na(names(First.year.catch.daily))]

nnn=1:length(SP.list)


# DEFINE TARGET SPECIES INDEX ----------------------------------------------
Tar.sp=match(TARGETS,SP.list)
Non.Tar.sp=nnn[which(!nnn%in%Tar.sp)]

# EFFECTIVE AREA (90% of catch) AND RASTER -----------------------------------------------------------------------
fn.scale=function(x,scaler) ((x/max(x,na.rm=T))^0.5)*scaler

Core=SP.list
pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Species core areas/cores.pdf'))
for(s in nnn)
{
  Kr=core.per
  d=Data.monthly.GN%>%filter(SPECIES%in%SP.list[[s]])%>%
    mutate(LAT=round(LAT),LONG=round(LONG))
  Nm=unique(d$SPECIES)
  Nm=ifelse(length(Nm)>3,'others',Nm)
  d=aggregate(LIVEWT.c~LAT+LONG,d,sum)
  d=d[order(-d$LIVEWT.c),]
  d$CumSum=cumsum(d$LIVEWT.c)
  d$CumSum=100*d$CumSum/max(d$CumSum)
  plot(d$LONG,d$LAT,cex=fn.scale(d$LIVEWT.c,4),pch=19,col="steelblue",
       ylab="Lat",xlab="Long",xlim=c(112,129),ylim=c(-36,-26),main=paste(Nm,names(SP.list)[s]))
  d=subset(d,CumSum<=Kr)
  Rnglat=range(d$LAT)
  Rnglon=range(d$LONG)
  polygon(c(Rnglon[1],Rnglon[2],Rnglon[2],Rnglon[1]),
          c(Rnglat[1],Rnglat[1],Rnglat[2],Rnglat[2]),border=2)
  Core[[s]]=list(Lat=Rnglat,Long=Rnglon)
}
dev.off()
if(Model.run=="First")
{
  Tab.d=Data.monthly.GN%>%filter(SPECIES%in%unlist(SP.list))%>%
    mutate(LAT=round(LAT),LONG=round(LONG))%>%
    group_by(Same.return,LAT,LONG,SPECIES)%>%
    summarise(Catch=sum(LIVEWT.c))%>%
    spread(SPECIES,Catch)
  
  theme_set(theme_pubr())
  nnn.i=2*nnn
  plot_list=vector('list',length(nnn)*2)
  for(s in nnn)
  {
    d=Tab.d[,match(c("LAT","LONG",SP.list[[s]]),names(Tab.d))]
    colnames(d)[3]="Catch"
    
    p1=ggplot(d, aes(LONG, LAT)) +
      geom_raster(aes(fill = log(Catch+1e-5)), interpolate = F)+
      #geom_contour(aes(z = log(Catch+1e-5)),linetype=1,col='black')+
      scale_fill_gradient2(low="white", high="dodgerblue4", guide="colorbar")+
      labs(title = paste(names(SP.list)[s],"   Density using zero and non zero catch"),
           fill = "log catch")
    
    dd=d%>%group_by(LAT,LONG)%>%summarise(Catch=sum(Catch,na.rm=T)/1000)%>%mutate(Catch=ifelse(Catch==0,NA,Catch))
    p2=ggplot(dd, aes(x=LONG, y=LAT,size =Catch)) + geom_point(alpha=0.7) +labs(title ="Positive catch (tonnes)") 
    
    plot_list[[nnn.i[s]-1]]=p1
    plot_list[[nnn.i[s]]]=p2
  }
  multi.page <-ggarrange(plotlist=plot_list, nrow = 2, ncol = 1)
  ggexport(multi.page, filename = handl_OneDrive("Analyses/Catch and effort/Outputs/species core areas/raster.pdf"))
}

# Fishery spatial expansion -----------------------------------------------------------------------
if(Model.run=="First")
{
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Spatial expansion_Monthly.pdf'))
  for(s in 1:length(SP.list))
  {
    print(paste('Spatial expansion Monthly----',names(SP.list)[s]))
    d.ktch=Data.monthly.GN%>%filter(SPECIES%in%SP.list[[s]])%>%
      mutate(LAT=round(LAT),LONG=round(LONG))%>%
      group_by(LAT,LONG,Same.return,YEAR.c)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
    
    d.eff=Effort.monthly%>%filter(Same.return%in%d.ktch$Same.return)%>%
      group_by(LAT,LONG,Same.return,YEAR.c)%>%
      summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c,na.rm=T),
                Km.Gillnet.Days.c=max(Km.Gillnet.Days.c,na.rm=T))
    
    d=inner_join(d.ktch,d.eff,by=c("LAT","LONG","Same.return","YEAR.c"))%>%
      mutate(cpue_days=LIVEWT.c/Km.Gillnet.Days.c,
             cpue_hours=LIVEWT.c/Km.Gillnet.Hours.c)
    
    p=d%>%
      filter(YEAR.c<=2005)%>%
      ggplot(aes(LONG,LAT,size=cpue_hours))+
      geom_point(color=2)+
      facet_wrap(~YEAR.c)+
      ggtitle(names(SP.list)[s])+
      theme(legend.position = 'top')
    
    print(p)
    
  }
  dev.off()
}

# Spatial re distribution due to climate change -----------------------------------------------------------------------
if(Model.run=="First")   
{
  #summer SST anomaly for the West Coast Bioregion
  #source: #eyeballed from Arani's 2023 WA’s Ocean Climate Report
  SST.anomaly.WCB=data.frame(year=2006:2023,
                             Anomaly=c(-.65,-.1,.7,.15,-.05,1.8,1.4,1,.1,
                                       .4,-.35,-.9,-.5,-.7,.2,.55,.8,.25))   
  dis.sp=table(Data.daily.GN$RSCommonName)
  dis.sp=subset(dis.sp,dis.sp>1e3)
  dis.sp=names(dis.sp)
  dis.sp=subset(dis.sp,!dis.sp%in%c("Sharks")) 
  Taxa=list(elasmos='elasmobranchs',
            scalies=c('demersal','est, coastal','pelagic'))
  Min.x=min(Data.daily.GN$LongFC)
  Max.x=max(Data.daily.GN$LongFC)
  Min.y=min(Data.daily.GN$LatFC)
  Max.y=max(Data.daily.GN$LatFC)
  for(q in 1:length(Taxa))
  {
    out.folder=names(Taxa)[q]
    estas=unique(Data.daily.GN%>%filter(type%in%Taxa[[q]])%>%pull(RSCommonName))
    dis.sp1=dis.sp[which(dis.sp%in%estas)]
    
    for(s in 1:length(dis.sp1))
    {
      print(paste('Spatial distribution -Climate Change - Daily----',dis.sp1[s]))
      
      d.ktch=Data.daily.GN%>%
        mutate(RSCommonName=ifelse(RSCommonName=="Smooth Hammerhead Shark" & YEAR.c<=2023,
                                   "Hammerhead Sharks",  #HH not id to species level prior 2023
                                   RSCommonName))%>%
        filter(RSCommonName%in%dis.sp1[s])%>%
        group_by(LatFC,LongFC,Same.return.SNo,YEAR.c,MONTH)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
      
      if(nrow(d.ktch)>0)
      {
        NM=capitalize(tolower(dis.sp1[s]))
        NM=ifelse(NM=="Bronze whaler","Copper shark",
                  ifelse(NM=="Dusky Whaler","Dusky shark",
                         ifelse(NM=="Wobbegong","Wobbegongs",NM)))
        
        d.eff=Effort.daily%>%filter(Same.return.SNo%in%d.ktch$Same.return.SNo)%>%
          group_by(Same.return.SNo)%>%
          summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c,na.rm=T),
                    Km.Gillnet.Days.c=max(Km.Gillnet.Days.c,na.rm=T))
        
        d=inner_join(d.ktch,d.eff,by=c("Same.return.SNo"))%>%
          mutate(cpue_days=LIVEWT.c/Km.Gillnet.Days.c,
                 cpue_hours=LIVEWT.c/Km.Gillnet.Hours.c,
                 #Yr_Mn=paste(MONTH,YEAR.c,sep='_'),
                 Yr_Mn=YEAR.c)%>%
          group_by(Yr_Mn,LatFC,LongFC)%>%
          summarise(cpue=mean(cpue_hours))
        
        p=d%>%
          left_join(SST.anomaly.WCB,by=c('Yr_Mn'='year'))%>%
          rename(SST.anomaly=Anomaly)%>%
          ggplot(aes(LongFC,LatFC,size=cpue,color=SST.anomaly))+
          geom_point(alpha=0.5)+
          facet_wrap(~Yr_Mn)+
          scale_colour_gradient2(low = "navyblue",mid = "forestgreen", high = "red")+
          theme_PA()+ylab('Latitude (°S)')+xlab('Longitude (°E)')+
          theme(legend.position = 'top')+
          xlim(Min.x,Max.x)+
          scale_y_continuous(limits = c(Min.y,Max.y),
                             breaks = pretty(Data.daily.GN$LatFC),
                             labels = abs(pretty(Data.daily.GN$LatFC)))
        print(p)
        ggsave(handl_OneDrive(paste('Analyses/Catch and effort/Outputs/Spatial_Climate change/',
                                    out.folder,'/',NM,'.tiff',sep='')),width=10,height= 10,compression="lzw")  
        
      }
      
      
    } 
  }
  
  
}

# Check whiskery shark change in catchability -----------------------------------------------------------------------
if(Model.run=="First")
{
  whis.q.period1=c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81","1981-82")
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Paper/Whiskery_targeting.pdf'))
  d=Data.monthly.GN%>%filter(!FINYEAR%in%Data.daily.GN$FINYEAR)%>%
    mutate(LAT=round(LAT),LONG=round(LONG))
  prop_period1=d%>%filter(FINYEAR%in%whis.q.period1)%>%
    group_by(BLOCKX, SPECIES)%>%
    summarise (n = sum(LIVEWT.c)) %>%
    mutate(freq_per1 = n / sum(n)) %>%
    data.frame
  
  prop_period2=d%>%filter(!FINYEAR%in%whis.q.period1)%>%
    group_by(BLOCKX, SPECIES)%>%
    summarise (n = sum(LIVEWT.c)) %>%
    mutate(freq_per2 = n / sum(n)) %>%
    data.frame
  
  prop=full_join(prop_period1,prop_period2,by=c("BLOCKX","SPECIES"))%>%
    filter(SPECIES==17003)%>%
    dplyr::select(BLOCKX,SPECIES,freq_per1,freq_per2)
  barplot(as.matrix(prop[,3:4]),beside =T,main="Prop of whiskery per block") 
  
  
  prop=d%>%filter(SPECIES%in%c(17003,18003))%>%
    group_by(SPECIES,FINYEAR)%>%
    summarise (n = sum(LIVEWT.c))
  
  prop.blk.yr=prop%>%group_by(FINYEAR)%>%
    summarise (n.blk.yr = sum(n)) 
  prop=prop%>%left_join(prop.blk.yr,by=c('FINYEAR'))%>%
    mutate(prop = n / n.blk.yr) %>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
    data.frame%>%arrange(FINYEAR)
  with(prop%>%filter(SPECIES==17003),plot(year,prop,type='b',
                                          ylim=c(0,.8),cex=1.25,
                                          ylab="Proportion of whiskery out of whiskery or dusky"))
  with(prop%>%filter(SPECIES==18003),points(year+.3,prop,type='b',pch=21
                                            ,cex=1.25,bg=3))
  legend("bottomright",c("whiskery","dusky"),pch=21,cex=1.5,pt.bg=c("white","green"),bty='n')
  dev.off()
}

# adjust core areas following Rory McAuley -----------------------------------------------------------------------
if(change.core.manually)
{
  Core$"Dusky Whaler"$Lat[2]=Dusky.range[1]
  Core$"Dusky Whaler"$Long[2]=Dusky.range[2]
  
  Core$"Sandbar Shark"$Long[2]=Sandbar.range[2]
  
  Core$"Whiskery Shark"$Lat[2]=Whiskery.range[1]
  Core$"Whiskery Shark"$Long[2]=Whiskery.range[2]
  
  Core$"Gummy Shark"$Long=Gummy.range
}

# FURTHER DATA MANIPULATIONS -----------------------------------------------------------------------

#put date back in Daily data set
get.dates=subset(Effort.daily,Same.return.SNo%in%unique(Data.daily.GN$Same.return.SNo),select=c(Same.return.SNo,date))
get.dates=get.dates[!duplicated(get.dates$Same.return.SNo),]
Data.daily.GN=Data.daily.GN%>%
  left_join(get.dates,by="Same.return.SNo")%>%
  mutate(date=as.Date(date))


#create some useful vars
Post.yrs=max(unique(sort(Data.monthly.GN$YEAR.c)))
Post.yrs=paste(2006:(Post.yrs-1),substr(2007:Post.yrs,3,4),sep="-")
Daily.l.years=sort(unique(Data.daily.GN$FINYEAR))

#Monthly
# Set all records of reapportioned returns to bad
Data.monthly.GN$Reporter.old=Data.monthly.GN$Reporter  #reset to Reporter, this was lost by FISHCUBE
Baddies=subset(Data.monthly.GN,Reporter.old=="bad")   
Baddies=unique(Baddies$Same.return)
Data.monthly.GN$Reporter.old=with(Data.monthly.GN,ifelse(Same.return%in%Baddies,'bad',Reporter.old))
Data.monthly.GN$Reporter=Data.monthly.GN$Reporter.old  

#remove small net length which correspond to non-shark gillnet
Data.monthly.GN=subset(Data.monthly.GN,NETLEN.c >100)

#remove caess data post 2005/06 in monthly records
Data.monthly.GN=subset(Data.monthly.GN,!(FINYEAR%in%Post.yrs & TYPE.DATA=="monthly"))

#remove data daily records in data monthly
Data.monthly.GN=subset(Data.monthly.GN,!FINYEAR%in%Daily.l.years)

#deal with boundary blocks in zone 1-2
if(BOUND.BLK=="REMOVE")Data.monthly.GN=subset(Data.monthly.GN,!BLOCKX%in%Boundary.Blks)
if(BOUND.BLK=="REALLOCATE")Data.monthly.GN$zone=with(Data.monthly.GN,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))

Data.monthly.GN$Same.return=with(Data.monthly.GN,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))

#Create monthly Km.Gillnet.Hours_shot.c
Eff$Km.Gillnet.Hours_shot.c=with(Eff,Km.Gillnet.Hours.c*SHOTS.c)

# Set all records of reapportioned returns to bad
if('Reporter.old'%in%names(Data.daily.GN))
{
  Baddies=subset(Data.daily.GN,Reporter.old=="bad")   
  if(nrow(Baddies)>0)
  {
    Baddies=unique(Baddies$Same.return)
    Data.daily.GN$Reporter.old=with(Data.daily.GN,ifelse(Same.return%in%Baddies,'bad',Reporter.old))
  }
  Data.daily.GN$Reporter=Data.daily.GN$Reporter.old 
}

#remove small net length which correspond to non-shark gillnet
Data.daily.GN=subset(Data.daily.GN,netlen.c >100)

#deal with boundary blocks in zone 1-2
if(BOUND.BLK=="REMOVE")Data.daily.GN=subset(Data.daily.GN,!BLOCKX%in%Boundary.Blks)
if(BOUND.BLK=="REALLOCATE")Data.daily.GN$zone=with(Data.daily.GN,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))

Data.daily.GN$Same.return=with(Data.daily.GN,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))

#Create daily Km.Gillnet.Hours_shot.c
Eff.daily$Km.Gillnet.Hours_shot.c=with(Eff.daily,Km.Gillnet.Hours.c*shots.c)


# ADD ENVIRONMENTAL VARIABLES -----------------------------------------------------------------------

# Add T and T residuals 
Lat.rng=seq(min(Data.monthly.GN$LAT),max(Data.monthly.GN$LAT),1)
Long.rng=seq(min(Data.monthly.GN$LONG),max(Data.monthly.GN$LONG),1)
SST=SST%>%filter(Lat%in%Lat.rng & Long%in%Long.rng)%>%
  rename(Temperature=value)

#Monthly
Data.monthly.GN=Data.monthly.GN %>%
  left_join(SST,by=c("YEAR.c"="year","MONTH"="month","LONG"="Long","LAT"="Lat"))%>% 
  arrange(YEAR.c,MONTH,LAT,LONG)%>%
  mutate(Temperature=ifelse(is.na(Temperature),na.approx(Temperature),Temperature))%>%
  group_by(MONTH,BLOCKX)%>%
  mutate(Temp.res=Temperature-mean(Temperature,na.rm=T))%>%   #temperature anomaly
  data.frame
#Daily
Data.daily.GN=Data.daily.GN %>%
  mutate(LAT.round=-floor(abs(LAT)),LONG.round=floor(LONG)) %>%
  left_join(SST,by=c("YEAR.c"="year","MONTH"="month","LONG.round"="Long","LAT.round"="Lat"))%>% 
  arrange(YEAR.c,MONTH,LAT.round,LONG.round)%>%
  mutate(Temperature=ifelse(is.na(Temperature),na.approx(Temperature),Temperature))%>%
  group_by(MONTH,BLOCKX)%>%
  mutate(Temp.res=Temperature-mean(Temperature,na.rm=T))%>%
  dplyr::select(-c(LONG.round,LAT.round))%>%   
  data.frame

# Add SOI, Freo and Moon (the latter to daily only) 
#SOI is already an anomaly
Freo=Freo%>%
  mutate(MeanSeaLevel=MeanSeaLevel-mean(MeanSeaLevel,na.rm=T))%>% #Fit Freo as an anomaly
  rename(Freo=MeanSeaLevel)%>%
  mutate(Freo_lag6=lag(Freo,6),
         Freo_lag12=lag(Freo,12))

      #Monthly
Data.monthly.GN=Data.monthly.GN%>%left_join(SOI,by=c("YEAR.c"="Year","MONTH"="Month"))%>%
  left_join(Freo,by=c("YEAR.c"="Year","MONTH"="Month")) 

      #Daily
Data.daily.GN=Data.daily.GN%>%left_join(SOI,by=c("YEAR.c"="Year","MONTH"="Month"))%>%
  left_join(Freo,by=c("YEAR.c"="Year","MONTH"="Month")) %>%
  mutate(Lunar=lunar.illumination(date),
         Lunar.phase=lunar.phase(date,name=T))

    #show example of temperature, temp.res and lunar data
if(Model.run=="First") 
{
  Data.daily.GN%>%
    filter((LAT==(-29) & LONG==114) | (LAT==(-30) & LONG==114) |
             (LAT==(-33) & LONG==114) | (LAT==(-34) & LONG==114))%>%
    distinct(date,YEAR.c,LAT,LONG,Temperature)%>%
    mutate(YEAR.c=as.character(YEAR.c),
           LAT=abs(LAT))%>%
    arrange(date)%>%
    ggplot(aes(date,Temperature,color=YEAR.c))+
    geom_line()+
    facet_grid(LAT~LONG)+theme_PA()+xlab('')+
    theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/example_Temperature.tiff'),
         width = 8,height = 8, dpi = 300, compression = "lzw")
  
  Data.daily.GN%>%
    filter((LAT==(-29) & LONG==114) | (LAT==(-30) & LONG==114) |
             (LAT==(-33) & LONG==114) | (LAT==(-34) & LONG==114))%>%
    distinct(date,YEAR.c,LAT,LONG,Temp.res)%>%
    mutate(YEAR.c=as.character(YEAR.c),
           LAT=abs(LAT))%>%
    arrange(date)%>%
    ggplot(aes(date,Temp.res,color=YEAR.c))+
    geom_line()+geom_hline(yintercept=0, linetype="dashed",color = "black")+
    facet_grid(LAT~LONG)+theme_PA()+xlab('')+
    theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/example_Temp.res.tiff'),
         width = 8,height = 8, dpi = 300, compression = "lzw")
  
  
  Data.daily.GN%>%
    filter(FINYEAR=='2010-11')%>%
    distinct(date,MONTH,Lunar,Lunar.phase)%>%
    arrange(date)%>%
    ggplot(aes(date,Lunar,color=Lunar.phase))+
    geom_point()+
    facet_wrap(~MONTH,scales='free_x')+theme_PA()+xlab('')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/example_Lunar phase.tiff'),
         width = 10,height = 6, dpi = 300, compression = "lzw")
  
}


#Replace 0 depth with mean of block10
Eff.daily=Eff.daily %>% 
  group_by(block10) %>%
  mutate(Mean.depth= replace(Mean.depth, Mean.depth<5, mean(Mean.depth, na.rm=TRUE)),
         Mean.depth=ifelse(Mean.depth==0,NA,Mean.depth))

#Set records with km.gn.hours > Max possible to 'bad' reporter  
Eff=Eff%>%
  mutate(Eff.Reporter=ifelse(Km.Gillnet.Hours.c>Max.km.gn.h.monthly,'bad',Eff.Reporter))   
Eff.daily=Eff.daily%>%
  mutate(Eff.Reporter=ifelse(Km.Gillnet.Hours.c>Max.km.gn.h.daily,'bad',Eff.Reporter),
         Eff.Reporter=ifelse(Mean.depth>120,'bad',Eff.Reporter))  


# CREATE SKIPPER EXPERIENCE VARIABLE ----------------------------------------------
#note: for daily, the combo 'VESSEL & MastersName' is used but for monthly, only 'VESSEL' is available
Data.daily.GN=Data.daily.GN%>%
  mutate(MastersName=case_when(MastersName=="1st trip - j. smythe - 2nd trip - n. triantaryllou"~"triantafyllou, neoclis",
                               MastersName=="andrew francis joy & m. tonkin"~'tonkin, michael',
                               MastersName%in%c('john rowbottom & chris black','c. black & mark robinson')~'black, christopher barry',
                               MastersName=='greg whetstone / mason thomas'~'whetstone,greg',
                               MastersName%in%c('brian / jason scimone','jason / brian scimone')~'scimone, brian',
                               MastersName%in%c('leighton matthews','l. matthews')~'matthews, leighton',
                               MastersName=='peter hughes & william ronald reay'~'reay, william ronald',
                               TRUE~MastersName))
if(Model.run=="First") 
{
  Table.experience_daily=Data.daily.GN%>%
    distinct(VESSEL,MastersName,Same.return.SNo,YEAR.c)%>%
    group_by(VESSEL,MastersName,YEAR.c)%>%
    tally()%>%
    rename(n.shots=n)%>%
    mutate(dummy=1)%>%
    group_by(VESSEL,MastersName)%>%
    mutate(Years.of.exp=cumsum(dummy))
  
  Table.experience_monthly=Data.monthly.GN%>%
    distinct(VESSEL,Same.return,YEAR.c)%>%
    group_by(VESSEL,YEAR.c)%>%
    tally()%>%
    rename(n.shots=n)%>%
    mutate(dummy=1)%>%
    group_by(VESSEL)%>%
    mutate(Years.of.exp=cumsum(dummy))
  
  #plot it
  Table.experience_daily%>%
    mutate(VESSEL_MastersName=paste(VESSEL, MastersName,sep='-'))%>%
    ggplot(aes(YEAR.c,VESSEL_MastersName))+
    geom_point(aes(size=n.shots),color='darkred')+theme(axis.text.y = element_text(size = 6))
  ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Efficiency creep/Experience_daily.tiff'),
         width = 6,height = 10, dpi = 300, compression = "lzw")
  
  
  A=Table.experience_monthly%>%
    group_by(VESSEL)%>%
    summarise(n=sum(n.shots))%>%
    arrange(-n)%>%
    mutate(Cumsum=cumsum(n),
           Per=Cumsum/sum(n))
  
  Table.experience_monthly=Table.experience_monthly%>%
    filter(VESSEL%in%unique(c(A%>%filter(Per<=0.8)%>%pull(VESSEL),Table.experience_daily$VESSEL)))
  
  Table.experience_monthly%>%
    ggplot(aes(YEAR.c,VESSEL))+
    geom_point(aes(size=n.shots),color='darkred')+theme(axis.text.y = element_text(size = 6))
  ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Efficiency creep/Experience_monthly.tiff'),
         width = 12,height = 10, dpi = 300, compression = "lzw")
  
}


# Extract first year of reporting by Vessel as proxy to experience
Vessel.start.rep=Data.monthly.GN%>%
  distinct(FINYEAR,VESSEL)%>%
  mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
  group_by(VESSEL)%>%
  mutate(Start.reporting=min(finyear))%>%
  ungroup()%>%
  distinct(VESSEL,Start.reporting,.keep_all = T)%>%
  dplyr::select(VESSEL,Start.reporting)

Vessel.start.rep.daily=Data.daily.GN%>%
  distinct(FINYEAR,VESSEL)%>%
  mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
  group_by(VESSEL)%>%
  mutate(Start.reporting=min(finyear))%>%
  ungroup()%>%
  distinct(VESSEL,Start.reporting,.keep_all = T)%>%
  dplyr::select(VESSEL,Start.reporting)%>%
  filter(!VESSEL%in%Vessel.start.rep$VESSEL)

Vessel.start.rep=rbind(Vessel.start.rep,Vessel.start.rep.daily)

# REMOVE DAILY RECORDS FROM MONTHLY EFFORT-----------------------------------------------------------------------
Effort.monthly=Effort.monthly%>%
  filter(FINYEAR%in%unique(Data.monthly.GN$FINYEAR))
Eff=Eff%>%
  filter(FINYEAR%in%unique(Data.monthly.GN$FINYEAR))

# EFFICIENCY_ADD SKIPPER'S YEARS OF EXPERIENCE ----------------------------------------------
#note: Vessel rego is used as proxy to skippers experience
Data.monthly.GN=Data.monthly.GN%>%
                  left_join(Vessel.start.rep,by='VESSEL')%>%
                  mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
                         Yrs.of.experience=finyear-Start.reporting,
                         Yrs.of.experience=ifelse(Yrs.of.experience==0,1,Yrs.of.experience))%>%
                  dplyr::select(-c(Start.reporting,finyear))

Data.daily.GN=Data.daily.GN%>%
                  left_join(Vessel.start.rep,by='VESSEL')%>%
                  mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
                         Yrs.of.experience=finyear-Start.reporting,
                         Yrs.of.experience=ifelse(Yrs.of.experience==0,1,Yrs.of.experience))%>%
                  dplyr::select(-c(Start.reporting,finyear))

# CREATE SPECIES DATA SETS FOR STANDARDISATIONS ----------------------------------------------

#Monthly  
#note: select species within CORE area and add effort by Same return
#      #remove Bronze Whaler due to few records and uncertain species ID prior to Daily logbooks
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)   # takes 1 sec per species
clusterEvalQ(cl, .libPaths(.libPaths()[1]))  #location of used libraries "C:/Users/myb/AppData/Local/R/win-library/4.4"
system.time({Species.list=foreach(s=nnn,.packages=c('dplyr','doParallel')) %dopar%
  {
    if(!SP.list[[s]]==18001) return(fn.cpue.data(Dat=Data.monthly.GN %>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                                                  LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                                                 EffrrT=Eff%>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                                        LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                                                 sp=SP.list[[s]]))
  }
})    #takes 3 secs per species
names(Species.list)=names(SP.list) 

#Daily 
#note: select species within CORE area and add effort by date or ID (==Same.return.SNo). Note that for catch aggregating by date
#       or by ID makes no difference but it's needed for merging with effort
system.time({Species.list.daily=foreach(s=nnn,.packages=c('dplyr','doParallel')) %dopar%
  {
    return(fn.cpue.data.daily(Dat=Data.daily.GN %>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                             LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                              EffrrT=Eff.daily%>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                           LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                              sp=SP.list[[s]]))
  }
})    #takes 13 secs per species
names(Species.list.daily)=names(SP.list) 

#Unbalanced raw data
get.raw=TRUE
if(get.raw)
{
  #Monthly
  system.time({Species.list.raw=foreach(s= 1:length(Tar.sp),.packages=c('dplyr','doParallel')) %dopar%
    {
      a=fn.cpue.data(Dat=Data.monthly.GN,
                     EffrrT=Eff,
                     sp=SP.list[[Tar.sp[s]]])%>%
        filter(SPECIES==SP.list[[Tar.sp[s]]])%>%
        dplyr::select(FINYEAR,SPECIES,LIVEWT.c,Km.Gillnet.Hours.c)
      return(a)
    }
  })
  names(Species.list.raw)=names(SP.list)[Tar.sp] 
  
  #Daily 
  system.time({Species.list.daily.raw=foreach(s= 1:length(Tar.sp),.packages=c('dplyr','doParallel')) %dopar%
    {
      a=fn.cpue.data.daily(Dat=Data.daily.GN,
                           EffrrT=Eff.daily,
                           sp=SP.list[[Tar.sp[s]]])%>%
        filter(SPECIES==SP.list[[Tar.sp[s]]])%>%
        dplyr::select(FINYEAR,SPECIES,LIVEWT.c,Km.Gillnet.Hours.c)
      return(a)
    }
  })
  names(Species.list.daily.raw)=names(SP.list)[Tar.sp] 
  
}

#Remove irrelevant variables (not used after prelim analysis)
for(s in nnn)
{
  if(!is.null(Species.list[[s]]))
  {
    Species.list[[s]] = Species.list[[s]] %>%
      dplyr::select(-c(LIVEWT,Boundary.blk,Km.Gillnet.Hours_shot.c,
                       TYPE.DATA,Sch.or.DogS,Freo_lag6,Freo_lag12,mesh,
                       NETLEN.c, BDAYS.c))
    
  }
  if(!is.null(Species.list.daily[[s]]))
  {
    Species.list.daily[[s]] = Species.list.daily[[s]] %>%
      dplyr::select(-c(LIVEWT,Km.Gillnet.Days.inv,
                       Km.Gillnet.Hours.inv,Km.Gillnet.Hours_shot.c,netlen.c,
                       bdays.c,TYPE.DATA,LIVEWT,nfish,Freo_lag6,Freo_lag12))
  }
}


# VESSELS & BLOCKS----------------------------------------------
Table.species.by.vessel=Data.daily.GN%>%
  filter(SPECIES%in%unlist(SP.list))%>%
  distinct(VESSEL,FINYEAR,RSCommonName)%>%
  mutate(N=1)%>%
  group_by(FINYEAR,RSCommonName)%>%
  summarise(N=sum(N))%>%
  spread(RSCommonName,N)

# Extract overall number of vessels reporting catch of species per species range
N.VES=matrix(rep(NA,length(nnn)),ncol=length(nnn))
colnames(N.VES)=names(SP.list)
for(i in nnn)
{
  a=b=NULL
  if(!is.null(Species.list[[i]]))a=unique(subset(Species.list[[i]],SPECIES%in%SP.list[[i]])$VESSEL)
  if(!is.null(Species.list.daily[[i]]))b=unique(subset(Species.list.daily[[i]],SPECIES%in%SP.list[[i]])$VESSEL)
  N.VES[,i]=length(unique(c(a,b)))
}
setwd(handl_OneDrive('Analyses/Catch and effort'))
hndl=paste(getwd(),"/Outputs/Paper/",sep="")
write.csv(N.VES,paste(hndl,"All.Vessels.by.species.csv",sep=""),row.names=T)


# Extract number of blocks where sharks have been caught within core (effective) area
Tol.blks=SP.list
for(i in nnn)
{
  a=b=NULL
  if(!is.null(Species.list[[i]]))a=unique(subset(Species.list[[i]],SPECIES%in%SP.list[[i]] & LIVEWT.c>0)$BLOCKX)
  if(!is.null(Species.list.daily[[i]]))b=unique(subset(Species.list.daily[[i]],SPECIES%in%SP.list[[i]] & LIVEWT.c>0)$BLOCKX)
  Tol.blks[[i]]=sort(unique(c(a,b)))
}
Blks.by.species=do.call(c,Tol.blks)
write.csv(Blks.by.species,paste(hndl,"All.Blks.by.species.csv",sep=""),row.names=T)

Table.species.by.BLOCKX=Data.daily.GN%>%
  filter(SPECIES%in%unlist(SP.list))%>%
  distinct(BLOCKX,FINYEAR,RSCommonName)%>%
  mutate(N=1)%>%
  group_by(FINYEAR,RSCommonName)%>%
  summarise(N=sum(N))%>%
  spread(RSCommonName,N)

# Block weights   
if(use.blok.area=="YES")
{
  AREA.W=vector('list',length=N.species)
  names(AREA.W)=SPECIES.vec
  AREA.W_b10=AREA.W
  
  
  #GIS approach (complied by Dale Smith)
  #BLOCKX
  names(Whis.fishArea)=names(Gum.fishArea)=names(Dusky.fishArea)=
    names(Sand.fishArea)=c("BLOCKX","Fish.Area","MaxDepth")
  Whis.fishArea$Fish.Area=with(Whis.fishArea,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Gum.fishArea$Fish.Area=with(Gum.fishArea,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Dusky.fishArea$Fish.Area=with(Dusky.fishArea,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Sand.fishArea$Fish.Area=with(Sand.fishArea,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  
  AREA.W$"Whiskery shark"=Whis.fishArea
  AREA.W$"Gummy shark"=Gum.fishArea  
  AREA.W$"Dusky shark"=Dusky.fishArea
  AREA.W$"Sandbar shark"=Sand.fishArea
  
  if(is.na(match(3521,Whis.fishArea$BLOCKX)))AREA.W$"Whiskery shark"=rbind(AREA.W$"Whiskery shark",data.frame(BLOCKX=3521,Fish.Area=0.01,MaxDepth=Whis.fishArea$MaxDepth[1]))
  if(is.na(match(3319,Dusky.fishArea$BLOCKX)))AREA.W$"Dusky shark"=rbind(AREA.W$"Dusky shark",data.frame(BLOCKX=3319,Fish.Area=0.01,MaxDepth=Dusky.fishArea$MaxDepth[1]))
  
  
  #block10
  names(Whis.fishArea_b10)=names(Gum.fishArea_b10)=names(Dusky.fishArea_b10)=
    names(Sand.fishArea_b10)=c("BLOCKX","block10","Fish.Area","MaxDepth")
  Whis.fishArea_b10$Fish.Area=with(Whis.fishArea_b10,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Gum.fishArea_b10$Fish.Area=with(Gum.fishArea_b10,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Dusky.fishArea_b10$Fish.Area=with(Dusky.fishArea_b10,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  Sand.fishArea_b10$Fish.Area=with(Sand.fishArea_b10,ifelse(Fish.Area<0.01,0.01,Fish.Area))
  
  AREA.W_b10$"Whiskery shark"=Whis.fishArea_b10
  AREA.W_b10$"Gummy shark"=Gum.fishArea_b10  
  AREA.W_b10$"Dusky shark"=Dusky.fishArea_b10
  AREA.W_b10$"Sandbar shark"=Sand.fishArea_b10
  
  #equal weights
  AREA.W.equal=AREA.W
  AREA.W_b10.equal=AREA.W_b10
  for(q in 1:length(AREA.W.equal)) AREA.W.equal[[q]]$Fish.Area=1
  for(q in 1:length(AREA.W_b10.equal)) AREA.W_b10.equal[[q]]$Fish.Area=1
  
}


# Create useful vars
FINYEAR.monthly=as.character(unique(Data.monthly.GN$FINYEAR))
FINYEAR.monthly=sort(FINYEAR.monthly)
N.yrs=length(FINYEAR.monthly)

FINYEAR.daily=as.character(unique(Data.daily.GN$FINYEAR))
FINYEAR.daily=sort(FINYEAR.daily)
N.yrs.daily=length(FINYEAR.daily) 

FINYEAR.ALL=unique(c(FINYEAR.monthly,FINYEAR.daily))
FINYEAR.ALL=sort(FINYEAR.ALL)
N.yrs.ALL=length(FINYEAR.ALL)


# Determine indicative vessels and blocks        
#steps: 1. select vessels that meet criteria (fishing for at least Threshold.n.yrs or Threshold.n.yrs.daily
#         and catching at least MIN.ktch)
#       2. for those vessels, select blocks with at least MIN.obs.BLK years of observations
BLKS.used=vector('list',length(SP.list)) 
names(BLKS.used)=names(SP.list)
BLKS.not.used=VES.used=VES.not.used=
  BLKS.used.daily=BLKS.not.used.daily=BLKS_10.used.daily=BLKS_10.not.used.daily=
  VES.used.daily=VES.not.used.daily=BLKS.used
if(Remove.blk.by=="blk_only")  
{
  for(i in nnn)  
  {
    print(paste('Calculating blocks and vessels used for ---',names(Species.list)[i]))
    NM=names(SP.list)[i]
    
    #monthly
    if(!is.null(Species.list[[i]]))
    {
      finy=sort(unique(Species.list[[i]]$FINYEAR))
      Ves.sel.BC=Threshold.n.yrs.monthly
      BLK.sel.BC=MIN.obs.BLK
      Min.ktch=MIN.ktch*1
      if(length(SP.list[[i]])<3)
      {
        if(SP.list[[i]][1]%in%Indicator.sp)
        {
          Ves.sel.BC=Threshold.n.yrs.monthly
          BLK.sel.BC=MIN.obs.BLK
          Min.ktch=MIN.ktch 
          if(18007%in%SP.list[[i]]) finy=finy[-match(c("1985-86","1986-87","1987-88"),finy)]
        }
      }
      dummy=fn.see.all.yrs.ves.blks(a=subset(Species.list[[i]],FINYEAR%in%finy),SP=SP.list[[i]],
                                    NM=NM,what=".monthly",Ves.sel.BC=Ves.sel.BC,Ves.sel.sens=Threshold.n.yrs.sens,
                                    BLK.sel.BC=BLK.sel.BC,BLK.sel.sens=MIN.obs.BLK.sens,
                                    Min.ktch=Min.ktch,MIN.rec=MIN.records.yr)
      
      if(!is.null(dummy))
      {
        BLKS.used[[i]]=dummy$Blks.BC
        BLKS.not.used[[i]]=dummy$Drop.blks
        VES.used[[i]]=dummy$Ves.BC
        VES.not.used[[i]]=dummy$Drop.ves
        write.csv(dummy$Tab.sel.ves,paste0(handl_OneDrive("Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Selected and dropped vessels and blocks/"),
                                           paste0(NM,'.monthly_annual records for selected vessels.csv',sep="")),row.names = FALSE)
      }
      
    }
    #Daily
    if(!is.null(Species.list.daily[[i]]))
    {
      Ves.sel.BC=Threshold.n.yrs.daily
      BLK.sel.BC=MIN.obs.BLK
      Min.ktch=MIN.ktch*1 
      finy=unique(Species.list.daily[[i]]$FINYEAR)
      if(length(SP.list[[i]])<3)
      {
        if(SP.list[[i]][1]%in%Indicator.sp)
        {
          Ves.sel.BC=Threshold.n.yrs.daily
          BLK.sel.BC=MIN.obs.BLK
          Min.ktch=MIN.ktch 
        }
      }
      dummy=fn.see.all.yrs.ves.blks(a=subset(Species.list.daily[[i]],FINYEAR%in%finy),SP=SP.list[[i]],
                                    NM=NM,what=".daily",Ves.sel.BC=Ves.sel.BC,Ves.sel.sens=Threshold.n.yrs.sens,
                                    BLK.sel.BC=BLK.sel.BC,BLK.sel.sens=MIN.obs.BLK.sens,
                                    Min.ktch=Min.ktch,MIN.rec=MIN.records.yr)
      if(!is.null(dummy))
      {
        BLKS.used.daily[[i]]=dummy$Blks.BC
        BLKS.not.used.daily[[i]]=dummy$Drop.blks
        BLKS_10.used.daily[[i]]=dummy$Blks.BC_10
        BLKS_10.not.used.daily[[i]]=dummy$Drop.blks_10
        VES.used.daily[[i]]=dummy$Ves.BC
        VES.not.used.daily[[i]]=dummy$Drop.ves
        write.csv(dummy$Tab.sel.ves,paste0(handl_OneDrive("Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Selected and dropped vessels and blocks/"),
                                           paste0(NM,'.daily_annual records for selected vessels.csv',sep="")),row.names = FALSE)
      }
    }
  }
  
  #compared used and not used vessels and blocks
  fn.fig(handl_OneDrive("Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Proportion blocks & vessels used & not used"),1400,2400)
  par(mfrow=c(length(nnn),2),mar=c(1,1,.1,.3),oma=c(2,2,1,.1),las=1,mgp=c(1,.5,0),xpd=T)
  for(i in nnn)
  {
    NMS.arG=c("","")
    #Monthly
    if(i==nnn[length(nnn)]) NMS.arG=c("Blocks","Vessels")
    if(is.null(BLKS.used[[i]])) plot.new()
    if(!is.null(BLKS.used[[i]])) fn.comp.used_not.used(blk.used=BLKS.used[[i]],blk.not.used=BLKS.not.used[[i]],
                                                       vsl.used=VES.used[[i]],vsl.not.used=VES.not.used[[i]],
                                                       NMS.arg=NMS.arG)
    legend('bottomleft',names(Species.list)[i],bty='n',cex=.9)
    if(i==length(nnn)) mtext("Monthly",1,1.5)
    #Daily
    if(is.null(BLKS.used.daily[[i]])) plot.new()
    if(!is.null(BLKS.used.daily[[i]])) fn.comp.used_not.used(blk.used=BLKS.used.daily[[i]],blk.not.used=BLKS.not.used.daily[[i]],
                                                             vsl.used=VES.used.daily[[i]],vsl.not.used=VES.not.used.daily[[i]],
                                                             NMS.arg=NMS.arG)
    if(i==length(nnn)) mtext("Daily",1,1.5)
    
    if(i==1)legend('center',c('used','not used'),fill=c("chartreuse3","brown1"),horiz = T,cex=1.25,bg="white")
  }
  dev.off() 
}

# EFFICIENCY_EXPAND VESSEL CHARACTERISTICS FOR SELECTED VESSELS----------------------------------------------
#select only vessels used for cpue stand
if(get.efficiency.creep)
{
  Vessel.charac=Vessel.charac%>%
    filter(BOATREGO%in%c(unique(unlist(VES.used)),unique(unlist(VES.used.daily))))
  ves.vars=names(Vessel.charac)
  ves.vars.drop=c("BOATNAME","BOATREGO","LICYEAR","PFL","SKIPPER")
  ves.vars=subset(ves.vars,!ves.vars%in%ves.vars.drop)
  
  #show history of vessel characteristics
  if(Model.run=="First")  
  {
    p.list=vector('list',length(ves.vars))
    pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Efficiency creep/Vessel_char_hist.pdf'))
    for(p in 1:length(p.list))
    {
      d=Vessel.charac%>%
        distinct(across(all_of(c('BOATNAME','BOATREGO',ves.vars[p]))))%>%
        mutate(BOATREGO.BOATNAME=paste(BOATREGO,BOATNAME,sep='_'))
      
      if(is.integer(d[,ves.vars[p]])|is.numeric(d[,ves.vars[p]]))
      {
        min.y=min(d[,ves.vars[p]],na.rm=T)-1
        p.list[[p]]=d%>%
          ggplot(aes_string(x="BOATREGO.BOATNAME",y=ves.vars[p]))+
          geom_segment( aes_string(x="BOATREGO.BOATNAME", xend="BOATREGO.BOATNAME", y=min.y, yend=ves.vars[p])) +
          geom_point( size=5, color="red", fill="orange", shape=21, stroke=2)+
          coord_flip()+ggtitle(ves.vars[p])+ylim(min.y,NA)
      }else
      {
        p.list[[p]]=d%>%
          ggplot(aes_string("BOATREGO.BOATNAME",fill=ves.vars[p]))+
          geom_bar()+
          coord_flip()+ggtitle(ves.vars[p])
      }
      print(p.list[[p]]+theme_PA())
    }
    dev.off()
  }
  
  #Expand Vessel.charac
  amend.date.built=TRUE
  Vessel.charac=Vessel.charac%>%
    mutate(DATEBUILT=case_when(BOATNAME=='CINDERELLA' & is.na(DATEBUILT)~1993,
                               BOATNAME=='COSAN II' & is.na(DATEBUILT)~2004,
                               BOATREGO=='E 061'& is.na(DATEBUILT)~2014, 
                               TRUE~DATEBUILT),
           BOATNAME=case_when(BOATNAME=='STEVE MAYREE D'~'SOUTHWESTERN',
                              BOATNAME=='PLANJAK'~'PLANJAK II',
                              BOATNAME=='FALCON 2'~'FALCON II',
                              BOATNAME=='ST. GERARD'~'ST GERARD M',
                              BOATNAME=='SOUTH WESTERN'~'SOUTHWESTERN',
                              BOATNAME=='SVET-NIKOLA'~'SVETI NIKOLA',
                              TRUE~BOATNAME))%>%
    filter(!BOATREGO%in%c('E 030','E 007')) #'E 030': 3 vessels with this rego, cannot determine period of time of each
  
  #create dummy data frame with observations for all years (expert judgement fill of missing years)
  dammy=Vessel.charac%>%
    distinct(BOATREGO,BOATNAME,LICYEAR)%>%
    group_by(BOATREGO,BOATNAME)%>%
    mutate(first.year=min(LICYEAR),
           last.year=max(LICYEAR))%>%
    ungroup()%>%
    distinct(BOATREGO,BOATNAME,.keep_all = TRUE)%>%
    mutate(first.year=case_when(BOATREGO=='G 297' & BOATNAME=='MISS-DEB-A-DEL II'~1998,
                                BOATREGO=='G 297' & BOATNAME=='ST GERARD M'~2000,
                                BOATREGO=='B 091' & BOATNAME=='PLANJAK II'~1997,
                                BOATREGO=='B 142' & BOATNAME=='SOUTHWESTERN'~2000,
                                BOATREGO=='E 061' & BOATNAME=='DESTINY'~2002,
                                BOATREGO=='E 056' & BOATNAME=='FATAL ATTRACTION'~1998,
                                TRUE~first.year),
           last.year=case_when( BOATREGO=='F 768' & BOATNAME=='CINDERELLA'~2003,
                                BOATREGO=='E 035' & BOATNAME=='BLUEBIRD II'~2002,
                                TRUE~last.year))
  if(Model.run=="First")
  {
    dammy%>% 
      gather(Yr,value,-c(BOATREGO,BOATNAME,LICYEAR))%>%
      arrange(BOATREGO,BOATNAME)%>%
      mutate(BOATREGO.BOATNAME=paste(BOATREGO,BOATNAME,sep='_'))%>%
      ggplot(aes(BOATREGO.BOATNAME,value))+
      geom_point(size=3)+geom_line()+
      coord_flip()+theme_PA()
    ggsave(handl_OneDrive('Analyses/Catch and effort/Outputs/Efficiency creep/Vessel_char_time line.tiff'),
           width = 10,height = 8,compression = "lzw")
  }
  
  if(amend.date.built)
  {
    first.year.data=Data.monthly.GN%>%
      filter(VESSEL%in%unique(Vessel.charac$BOATREGO))%>%
      distinct(VESSEL,FINYEAR)%>%
      mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
      group_by(VESSEL)%>%
      summarise(first.year.catch=min(finyear))
    dammy=left_join(dammy,first.year.data,by=c('BOATREGO'='VESSEL'))
  }
  
  dis.FINYEAR.daily=Data.daily.GN%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(VESSEL,BoatName)%>%
    summarise(last.year.ktch=max(finyear))%>%
    ungroup()%>%
    mutate(BoatName=toupper(BoatName))
  dammy=dammy%>% 
    ungroup()%>%
    left_join(dis.FINYEAR.daily,by=c('BOATNAME'='BoatName','BOATREGO'='VESSEL'))
  
  #reset first year to start at first year of catch data
  dammy=dammy%>%
    mutate(first.year=case_when(BOATREGO=='B 038' & BOATNAME=='SVETI NIKOLA'~1980,
                                BOATREGO=='B 067' & BOATNAME=='CHRISTINE NERELLE'~1980,
                                BOATREGO=='B 091' & BOATNAME=='FILLY JO'~1979,
                                BOATREGO=='B 098' & BOATNAME=='CINDERELLA'~1993,
                                BOATREGO=='B 142' & BOATNAME=='WARNBRO LADY'~1980,
                                BOATREGO=='E 034' & BOATNAME=='MARIAN'~1980,
                                BOATREGO=='E 035' & BOATNAME=='BLUEBIRD II'~1984,
                                BOATREGO=='E 045' & BOATNAME=='DOREEN'~1983,
                                BOATREGO=='E 056' & BOATNAME=='DESTINY'~1980,
                                BOATREGO=='E 061' & BOATNAME=='STORMRAKER'~1988,
                                BOATREGO=='E 067' & BOATNAME=='QUADRANT'~1998,
                                BOATREGO=='F 417' & BOATNAME=='SANTA BARBARA II'~1979,
                                BOATREGO=='F 517' & BOATNAME=='SAN MARGO'~1987,
                                BOATREGO=='F 768' & BOATNAME=='CINDERELLA'~1993,
                                BOATREGO=='G 297' & BOATNAME=='ELIZABETH MARIA II'~1985,
                                TRUE~first.year))
  #expand years
  dammy=dammy%>%
    mutate(last.year=ifelse(!is.na(last.year.ktch) & !BOATNAME%in%c('COSAN II'),last.year.ktch,last.year))%>%
    dplyr::select(BOATREGO,BOATNAME,first.year,last.year)%>%
    rowwise() %>% 
    do(data.frame(BOATREGO= .$BOATREGO,BOATNAME=.$BOATNAME, LICYEAR = .$first.year:.$last.year)) %>% 
    arrange(BOATREGO,BOATNAME,LICYEAR)
  
  #add expanded years to vessel char
  Vessel.charac.exp=full_join(Vessel.charac,
                              dammy,by=c('BOATNAME','BOATREGO','LICYEAR'))%>%
    mutate(FINYEAR=paste(LICYEAR,substr(LICYEAR+1,3,4),sep='-'),
           finyear=as.integer(substr(FINYEAR,1,4)))%>%
    arrange(BOATNAME,BOATREGO,finyear)%>%
    dplyr::select(-c(PFL,SKIPPER))
  
  
  #fill in missing years  
  Vessel.charac.exp=Vessel.charac.exp%>%
    group_by(BOATNAME,BOATREGO)%>%
    fill(ENGDERAT,ENGNUM,FRZCAP,DATEBUILT,FLYBRIDGE,
         HULLCONS,HULLNUMB,HULLTYPE,WHEELHOU,
         GPS,PLOT,COECHO,RADAR,SONAR,BWECHO, .direction = "downup")%>%
    fill(ENGPOWR,ENGSPD,BRNTNK,ICEBOX,LHTABV,LHTBLW,GROSSTON,MAXBEAM,
         MAXDRAU,REG_LENGTH, .direction = "downup")%>%
    ungroup()
  
  Vessel.charac.exp=Vessel.charac.exp%>%
    mutate(GPS=ifelse(finyear<GPS & !is.na(GPS),'N',ifelse(finyear>=GPS & !is.na(GPS),'Y',GPS)),
           PLOT=ifelse(finyear<PLOT & !is.na(PLOT),'N',ifelse(finyear>=PLOT & !is.na(PLOT),'Y',PLOT)),
           COECHO=ifelse(finyear<COECHO & !is.na(COECHO),'N',ifelse(finyear>=COECHO ,'Y',COECHO)),
           RADAR=ifelse(finyear<RADAR & !is.na(RADAR),'N',ifelse(finyear>=RADAR & !is.na(RADAR),'Y',RADAR)),
           SONAR=ifelse(finyear<SONAR & !is.na(SONAR),'N',ifelse(finyear>=SONAR & !is.na(SONAR),'Y',SONAR)),
           BWECHO=ifelse(finyear<BWECHO & !is.na(BWECHO),'N',ifelse(finyear>=BWECHO & !is.na(BWECHO),'Y',BWECHO)))%>%
    data.frame()
  
  check.filled.properly=FALSE  
  if(check.filled.properly)
  {
    x=Vessel.charac%>%
      mutate(x=paste(BOATREGO,BOATNAME))
    x1=x%>%distinct(BOATREGO,BOATNAME)%>%arrange(BOATNAME)
    x1=paste(x1$BOATREGO,x1$BOATNAME)
    a=Vessel.charac.exp%>%
      mutate(x=paste(BOATREGO,BOATNAME))
    for(i in 1:length(x1))
    {
      b=a%>%filter(x==x1[i])%>%data.frame
      b.ori=x%>%filter(x==x1[i])%>%data.frame
      p.list=vector('list',length(ves.vars))
      for(v in 1:length(ves.vars))
      {
        p.list[[v]]=b%>%
          ggplot(aes_string('LICYEAR',ves.vars[v]))+
          geom_point(color='orange')+
          geom_point(data=b.ori,aes_string('LICYEAR',ves.vars[v],color=ves.vars[v]),size=3)+
          theme_PA()+theme(legend.position = 'none')
      }
      ggarrange(plotlist = p.list,ncol=5,nrow=5)
      ggsave(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/check filled out correctly/',x1[i],'.tiff')),
             width = 12,height = 10,compression = "lzw")
    }
  }
  
  #show each vessel char tru time
  if(Model.run=="First")  
  {
    pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Efficiency creep/Vessel_char_thru time.pdf'))
    for(v in 1:length(ves.vars))
    {
      p=Vessel.charac.exp%>%
        mutate(Boat.name_rego=paste(BOATREGO,BOATNAME,sep='_'))%>%
        ggplot(aes_string('LICYEAR','Boat.name_rego', color=ves.vars[v]))+
        geom_point( size=5)+
        ggtitle(ves.vars[v])+ylab('')+
        theme_PA()+geom_vline(xintercept=2006,color='orange',linewidth=1.5,alpha=.6) 
      print(p)
    }
    dev.off()
  }
  
}

# ILLUSTRATE FOLLY EFFECT (MEAN vs SUM) AND CATCH RATE VARIABILITY---------------------------
if(Show.folly.eg=="YES")
{
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/MeanVsSum"))
  
  #5.1.1 Dummy example
  n=10
  Barcos=LETTERS[1:n]
  sipiui=10
  Consistent.cpue=data.frame(Vessel=Barcos,Effort=seq(10,100,length.out=n))
  Consistent.cpue$Catch=Consistent.cpue$Effort*sipiui
  Consistent.cpue$cpue=Consistent.cpue$Catch/Consistent.cpue$Effort
  Consistent.cpue$mean.cpue=round(mean(Consistent.cpue$cpue),2)
  Consistent.cpue$sd.mean.cpue=round(sd(Consistent.cpue$cpue),2)
  Consistent.cpue$sum.cpue=round(sum(Consistent.cpue$Catch)/sum(Consistent.cpue$Effort),2)
  
  Sum.under1=Consistent.cpue
  Sum.under1$Catch[c(4,6,8)]=0
  Sum.under1$cpue=Sum.under1$Catch/Sum.under1$Effort
  Sum.under1$mean.cpue=round(mean(Sum.under1$cpue),2)
  Sum.under1$sd.mean.cpue=round(sd(Sum.under1$cpue),2)
  Sum.under1$sum.cpue=round(sum(Sum.under1$Catch)/sum(Sum.under1$Effort),2)
  
  
  Sum.under2=Consistent.cpue
  Sum.under2$Catch[c(1,2,3)]=Sum.under2$Catch[c(1,2,3)]*10
  Sum.under2$cpue=Sum.under2$Catch/Sum.under2$Effort
  Sum.under2$mean.cpue=round(mean(Sum.under2$cpue),2)
  Sum.under2$sd.mean.cpue=round(sd(Sum.under2$cpue),2)
  Sum.under2$sum.cpue=round(sum(Sum.under2$Catch)/sum(Sum.under2$Effort),2)
  
  write.csv(Consistent.cpue,"Consistent.cpue.csv",row.names=F)
  write.csv(Sum.under1,"Sum.under1.csv",row.names=F)
  write.csv(Sum.under2,"Sum.under2.csv",row.names=F)
}
if(Show.variability.cpue.eg=="YES")
{
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs"))
  
  fn.fig("E.g.variability.cpue.presentation",2400,2000)
  par(mfcol=c(2,1),mai=c(.45,.6,.2,.1),oma=c(.8,.8,.1,.1),mgp=c(2.5,.65,0))
  
  
  #whiskery
  a=subset(Data.monthly.GN.whiskery,SPECIES==17003)
  agg1=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~VESSEL+FINYEAR,a,mean)
  All.yrS=sort(unique(agg1$FINYEAR))
  Ves1=unique(agg1$VESSEL)
  #YMAX=max(agg1[,3])
  YMAX=150
  plot(1:length(All.yrS),ylim=c(0,YMAX),xaxt='n',col="transparent",ylab="",xlab="",
       main="Annual cpue by vessel")
  axis(1,1:length(All.yrS),F,tck=-0.025)
  axis(1,seq(1,length(All.yrS),10),F,tck=-0.05)
  CL=rainbow(length(Ves1))
  for(i in 1:length(Ves1))
  {
    x=subset(agg1,VESSEL==Ves1[i])
    Mis.yr=All.yrS[which(!All.yrS%in%x$FINYEAR)]
    Mis.d=x[1:length(Mis.yr),]
    Mis.d$VESSEL=Ves1[i]
    Mis.d$FINYEAR=Mis.yr
    Mis.d[,3]=NA
    x=rbind(x,Mis.d)
    x=x[order(x$FINYEAR),]
    lines(x[,3],col=CL[i])
  }
  legend("topright","Whiskery shark",cex=2,bty='n')
  
  #gummy
  a=subset(Data.monthly.GN.gummy,SPECIES==17001)
  agg1=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~VESSEL+FINYEAR,a,mean)
  All.yrS=sort(unique(agg1$FINYEAR))
  Ves1=unique(agg1$VESSEL)
  #YMAX=max(agg1[,3])
  YMAX=150
  plot(1:length(All.yrS),ylim=c(0,YMAX),xaxt='n',col="transparent",ylab="",xlab="",
       main="")
  axis(1,1:length(All.yrS),F,tck=-0.025)
  CL=rainbow(length(Ves1))
  for(i in 1:length(Ves1))
  {
    x=subset(agg1,VESSEL==Ves1[i])
    Mis.yr=All.yrS[which(!All.yrS%in%x$FINYEAR)]
    Mis.d=x[1:length(Mis.yr),]
    Mis.d$VESSEL=Ves1[i]
    Mis.d$FINYEAR=Mis.yr
    Mis.d[,3]=NA
    x=rbind(x,Mis.d)
    x=x[order(x$FINYEAR),]
    lines(x[,3],col=CL[i])
  }
  legend("topright","Gummy shark",cex=2,bty='n')
  
  axis(1,seq(1,length(All.yrS),10),All.yrS[seq(1,length(All.yrS),10)],tck=-0.05,cex.axis=1.25)
  mtext("CPUE (kg /km.gn.day)",2,-0.75,las=3,outer=T,cex=1.75)
  mtext("Financial year",1.1,-.25,outer=T,cex=1.75)
  dev.off()
}

# CONSTRUCT WIDE DATABASE FOR STANDARDISATIONS ----------------------------------------------
#steps: 
#   1. select "Good" records (the variable "Reporter" includes good/bad catch and effort reporters)
#   2. Construct a single row for each record (i.e. 'year-month-vessel-block-gear' for monthly
#      returns and 'year-Session-vessel-block10-gear' for daily logbooks), with catch of target
#      and other species as separate columns, giving a 0 catch if no catch

DATA.list.LIVEWT.c=vector('list',length(SP.list)) 
names(DATA.list.LIVEWT.c)=names(SP.list)
DATA.list.LIVEWT.c.daily=DATA.list.LIVEWT.c
Prop.Catch=rep(NA,length(SP.list))
names(Prop.Catch)=names(SP.list)
Prop.Catch.daily=Prop.Catch

#monthly  
These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Hours.c","Km.Gillnet.Days.c",
                "zone","MONTH","BLOCKX","SHOTS.c","HOURS.c")
for(i in nnn)
{
  if(!(is.null(Species.list[[i]])))
  {
    dummy=Effort.data.fun(DATA=subset(Species.list[[i]],Reporter=="good"),
                          target=SP.list[[i]],
                          ktch="LIVEWT.c")  
    DATA.list.LIVEWT.c[[i]]=dummy$dat
    Prop.Catch[i]=dummy$prop.with.catch 
  }
  print(paste0('constructing wide Monthly dataframe for -----',names(DATA.list.LIVEWT.c)[i]))
}

#daily 
These.efforts.daily=c("FINYEAR","date","TSNo","Km.Gillnet.Hours.c","Km.Gillnet.Days.c",
                      "zone","MONTH","BLOCKX","block10","VESSEL","mesh",
                      "Same.return.SNo","nlines.c","shots.c",
                      "hours.c","BoatName","MastersName")
for(i in nnn)
{
  if(!is.null(Species.list.daily[[i]])) 
  {
    dummy=Effort.data.fun.daily(DATA=subset(Species.list.daily[[i]],Reporter=="good"),
                                target=SP.list[[i]],
                                ktch="LIVEWT.c",
                                Aggregtn="SNo")
    DATA.list.LIVEWT.c.daily[[i]]=dummy$dat
    Prop.Catch.daily[i]=dummy$prop.with.catch   
  }
  print(paste0('constructing wide Daily dataframe for -----',names(DATA.list.LIVEWT.c.daily)[i]))
}

#add vessel characteristics to data 
#note: by not having 'boatname' in Monthly data, I have to drop this var when merging,
#     this is an issue as many vessel regos were used in parallel in different boats, eg 'G 297'
#     this is was attempted to be solved in '#reset first year to start at first year of catch data' 
#     but still not enough, we need to have 'boatname' included in Monthly.returns and merge by BOATNAME-BOATREGO-FINYEAR
if(get.efficiency.creep)
{
  for(i in nnn)
  {
    print(paste('Adding vessel characteristicts to ---',names(DATA.list.LIVEWT.c)[i]))
    if(!is.null(DATA.list.LIVEWT.c[[i]]))
    {
      DATA.list.LIVEWT.c[[i]]=DATA.list.LIVEWT.c[[i]]%>%
        left_join(Vessel.charac.exp%>%
                    distinct(BOATREGO,FINYEAR,.keep_all = T)%>%
                    dplyr::select(-c(finyear,LICYEAR,BOATNAME)),
                  by=c('VESSEL'='BOATREGO','FINYEAR'))  
    }  
    if(!is.null(DATA.list.LIVEWT.c.daily[[i]]))
    {
      DATA.list.LIVEWT.c.daily[[i]]=DATA.list.LIVEWT.c.daily[[i]]%>%
        left_join(Vessel.charac.exp%>%
                    distinct(BOATREGO,FINYEAR,.keep_all = T)%>%
                    dplyr::select(-c(finyear,LICYEAR,BOATNAME)),
                  by=c('VESSEL'='BOATREGO','FINYEAR'))
    }
  }
}

# EXPORT PROPORTIONS WITH CATCH ----------------------------------------------
write.csv(Prop.Catch,paste(hndl,"Prop.records.with.catch.monthly.csv",sep=""),row.names=T)
write.csv(Prop.Catch.daily,paste(hndl,"Prop.records.with.catch.daily.csv",sep=""),row.names=T)

# GET NUMBER OF SPECIES BY DAILY RECORD ----------------------------------------------
if(Model.run=="First")  
{
  for(i in nnn)
  {
    TRGT=names(Species.list.daily)[i]
    fn.species.per.session(d=Species.list.daily[[i]],target=TRGT)
    ggsave(handl_OneDrive(paste0("Analyses/Catch and effort/Outputs/Number of species caught by Daily sessions/",TRGT,".tiff")),
           width = 6,height = 8,compression = "lzw")
  }
}

# EFFICIENCY_CATCH RATE BY FISHER TO ID FISHING EFFICIENCY CREEP ----------------------------------------------
#note: Monthly is aggregated so low incidence of 0 catch records. Not very informative
if(Model.run=="First")  
{
  #PROPORTION ZERO CATCH THRU TIME
  indis=match(c("Gummy Shark","Whiskery Shark","Dusky Whaler","Sandbar Shark"),names(DATA.list.LIVEWT.c))
  for(i in indis)
  {
    print(paste('PROPORTION ZERO CATCH THRU TIME for ----',names(DATA.list.LIVEWT.c)[i]))
    
    explained.ktch=0.90
    if(names(Species.list)[i]%in%c("Dusky Whaler","Whiskery Shark")) explained.ktch=0.8
    if(!is.null(DATA.list.LIVEWT.c[[i]]))
    {
      fn.prop.0.catch.by.fisher(d=DATA.list.LIVEWT.c[[i]]%>%
                                  mutate(Ves.var=VESSEL,
                                         time.var=YEAR.c+MONTH/12.95,
                                         cpue=Catch.Target/Km.Gillnet.Hours.c),
                                explained.ktch=explained.ktch,
                                NM=names(DATA.list.LIVEWT.c)[i],
                                series='Monthly')
    }
    
    explained.ktch=0.95
    if(!is.null(DATA.list.LIVEWT.c.daily[[i]]))
    {
      fn.prop.0.catch.by.fisher(d=DATA.list.LIVEWT.c.daily[[i]]%>%
                                  mutate(Ves.var=MastersName,
                                         time.var=YEAR.c+MONTH/12.95,
                                         cpue=Catch.Target/Km.Gillnet.Hours.c),
                                explained.ktch=explained.ktch,
                                NM=names(DATA.list.LIVEWT.c.daily)[i],
                                series='Daily')
    }
  }
  
  #CATCH BY FISHER  
  
}
# CALCULATE BLOCK CORNERS FOR GAM ----------------------------------------------
for(s in nnn)
{
  if(!is.null(DATA.list.LIVEWT.c[[s]])) DATA.list.LIVEWT.c[[s]] = DATA.list.LIVEWT.c[[s]] %>%  
      mutate(LAT10.corner=LAT, LONG10.corner=LONG)
  
  if(!is.null(DATA.list.LIVEWT.c.daily[[s]])) DATA.list.LIVEWT.c.daily[[s]] = DATA.list.LIVEWT.c.daily[[s]] %>% 
      mutate(LAT10.corner=-(abs(as.numeric(substr(block10,1,2))+10*(as.numeric(substr(block10,3,3)))/60)),
             LONG10.corner=100+as.numeric(substr(block10,4,5))+10*(as.numeric(substr(block10,6,6)))/60)
}

# REMOVE YEARS WITHOUT SANDBAR, TIGER or SPINNER SHARK CODE FROM DATA SET----------------------------------------------
#note: Dropping these years because vessels don't meet selection criteria and no positive catch
fn.sel.yrs.used=function(DD,ThrShld.n.vess)
{
  a=with(DD,table(FINYEAR,VESSEL))
  a[a>0]=1
  a=rowSums(a)
  a[a<ThrShld.n.vess]=NA
  return(names(a[which(!is.na(a))]))
}
DD=subset(DATA.list.LIVEWT.c$"Sandbar Shark",BLOCKX%in%as.numeric(BLKS.used$"Sandbar Shark"))      
DD=subset(DD,VESSEL%in%VES.used$"Sandbar Shark")
San.Yrs=fn.sel.yrs.used(DD,ThrShld.n.vess=5)
rm(DD)
DATA.list.LIVEWT.c$"Sandbar Shark"=subset(DATA.list.LIVEWT.c$"Sandbar Shark",FINYEAR%in%San.Yrs)

#no reported catch before these years for spinner and tiger
Spinner.yrs=paste(1991:2005,substr(1992:2006,3,4),sep='-')
DATA.list.LIVEWT.c$"Spinner Shark"=subset(DATA.list.LIVEWT.c$"Spinner Shark",FINYEAR%in%Spinner.yrs)

Tiger.yrs=paste(1995:2005,substr(1996:2006,3,4),sep='-')
DATA.list.LIVEWT.c$"Tiger Shark"=subset(DATA.list.LIVEWT.c$"Tiger Shark",FINYEAR%in%Tiger.yrs)


# FIX SOME SPECIES NAMES ----------------------------------------------
Nms.sp=capitalize(tolower(names(SP.list)))  
Nms.sp[match(c("Dusky whaler"),Nms.sp)]=c("Dusky shark")


# IDENTIFY FISHING ON DIFFERENT HABITATS (~TARGETING BEHAVIOUR, only applicable to Daily logbooks) ----------------------------------------------
#note:  more code in 2.CPUE standardisations_delta.R)
#      this uses all species accounting for 95% of catch
HndL.Species_targeting=handl_OneDrive("Analyses/Catch and effort/Outputs/Species targeting/")

#Density plots
if(Model.run=="First")
{
  theme_set(theme_pubr())
  nnn.i=4*nnn
  plot_list=vector('list',length(nnn)*4)
  for(s in nnn)
  {
    if(!is.null(DATA.list.LIVEWT.c.daily[[s]]))
    {
      d=DATA.list.LIVEWT.c.daily[[s]]
      with(d,plot(Catch.Target,Catch.Gummy))
      
      p1=ggplot(d, aes(x=Catch.Target, y=Catch.Gummy) ) +
        geom_bin2d(bins = 50) +
        theme_bw()+labs(title =names(DATA.list.LIVEWT.c.daily)[s])
      
      p2=ggplot(d, aes(x=Catch.Target, y=Catch.Whiskery) ) +
        geom_bin2d(bins = 50) +
        theme_bw()
      
      p3=ggplot(d, aes(x=Catch.Target, y=Catch.Dusky) ) +
        geom_bin2d(bins = 50) +
        theme_bw()
      
      
      p4=ggplot(d, aes(x=Catch.Target, y=Catch.Sandbar) ) +
        geom_bin2d(bins = 50) +
        theme_bw()
      
      plot_list[[nnn.i[s]-3]]=p1
      plot_list[[nnn.i[s]-2]]=p2
      plot_list[[nnn.i[s]-1]]=p3
      plot_list[[nnn.i[s]]]=p4
      
    }
  }
  multi.page <-ggarrange(plotlist=plot_list, nrow = 2, ncol = 2)
  ggexport(multi.page, filename = paste(HndL.Species_targeting,"Density.pdf",sep=''))
}

#Run Stephens & McCall 
if(do_Stephens_McCall=="YES")
{
  dir_plots=paste0(HndL.Species_targeting,'Stephens_McCall')
  Explained.catch.percent=0.95 #0.99 original
  Minyears.with.catch=5 #10 original
  Min.avrg.annual.catch=1 #1 original
  use.all.species.minlocs=TRUE #FALSE original
  check.dodgy.predictions=FALSE
  if(Model.run=="First") check.dodgy.predictions=TRUE
  
  Data.daily.GN.DD=Data.daily.GN%>% 
                    mutate(SNAME=case_when(SPECIES%in%c(19001:19004)~'Hammerhead Sharks',
                                             TRUE~SNAME),
                           RSCommonName=case_when(SPECIES%in%c(19001:19004)~'Hammerhead Sharks',
                                                    TRUE~RSCommonName),
                           RSSpeciesCode=case_when(SPECIES%in%c(19001:19004)~37019000,
                                                   TRUE~RSSpeciesCode),
                           RSSpeciesId=case_when(SPECIES%in%c(19001:19004)~90,
                                                 TRUE~RSSpeciesId),
                           SPECIES=case_when(SPECIES%in%c(19001:19004)~19000,
                                             TRUE~SPECIES))%>%
                          group_by(SNAME,SPECIES,FINYEAR,Same.return.SNo,zone)%>%
                          summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
                          ungroup()%>%
                          filter(!SNAME=="Nil Fish Caught")%>%
                          mutate(SNAME=case_when(SPECIES==22999~'Other sharks',
                                                 SPECIES==20000~'Dogfishes',
                                                 TRUE~SNAME))
  
  #1. Catch of target vs catch of others and calculation of minlocs
  minlocs.vec=rep(NA,length(SP.list)) #minimum number of co-occurrences for species to be considered
  if(!use.all.species.minlocs)
  {
    for( i in 1:length(SP.list))
    {
      output_dir=paste0(dir_plots,'/',names(SpiSis[i]))
      if(!dir.exists(output_dir)) dir.create(output_dir)
      print(paste('Catch of target vs catch of others for ----',names(SpiSis[i])))
      d=left_join(Data.daily.GN.DD%>%
                    filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[i]]$Same.return.SNo))%>%
                    dplyr::select(Same.return.SNo,zone,LIVEWT.c,SPECIES,SNAME),
                  Effort.daily%>%
                    filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[i]]$Same.return.SNo))%>%
                    mutate(FishingSeason=as.numeric(substr(finyear,1,4)))%>%
                    distinct(Same.return.SNo,Km.Gillnet.Hours.c,FishingSeason),
                  by='Same.return.SNo')%>%
        filter(!SPECIES==9998)
      Tab=d%>%
        group_by(SPECIES)%>%
        tally()%>%
        arrange(-n)%>%
        mutate(Cumsum=cumsum(n),
               Percent=Cumsum/sum(n))%>%
        mutate(id = row_number())
      Tab%>%
        ggplot(aes(id,Percent))+
        geom_point(alpha=0.35,size=3.5,color='orange')+
        xlab('Number of species')+ylab('Proportion of catch (by number)')+
        geom_text(data=Tab%>%filter(abs(Percent - 0.9) == min(abs(Percent - 0.9))),
                  aes(id,Percent,label=paste0(round(100*Percent),'%, ',id,
                                              ' species, minimum of ',n,' co-occurrences')), hjust = 0,color='black')+
        geom_text(data=Tab%>%filter(abs(Percent - 0.95) == min(abs(Percent - 0.95))),
                  aes(id,Percent,label=paste0(round(100*Percent),'%, ',id,
                                              ' species, minimum of ',n,' co-occurrences')), hjust = 0,color='black')+
        geom_text(data=Tab%>%filter(abs(Percent - 0.99) == min(abs(Percent - 0.99))),
                  aes(id,Percent,label=paste0(round(100*Percent),'%, ',id,
                                              ' species, minimum of ',n,' co-occurrences')), hjust = 0,color='black')+
        theme_PA()
      ggsave(paste0(output_dir,"/Cumulative catch.tiff"),width = 6, height = 6,compression="lzw")
      
      Tab=Tab%>%
        filter(Percent<=Explained.catch.percent)
      minlocs.vec[i]=min(Tab$n)
      
      d1=d%>%
        filter(SPECIES%in%Tab$SPECIES)%>%
        mutate(SNAME=case_when(SNAME=='Triggerfishes & Leatherjackets'~'Leatherjackets',
                               SNAME=='Temperate Basses & Rockcods'~'Rockcods',
                               TRUE~SNAME),
               CPUE=LIVEWT.c/Km.Gillnet.Hours.c,
               SP.ID=ifelse(SPECIES==SP.list[[i]],'Target',SNAME))%>%
        dplyr::select(CPUE,SP.ID,Same.return.SNo)%>%
        mutate(SP.ID=str_replace_all(SP.ID," ","_"))%>%
        spread(SP.ID,CPUE,fill=0)%>%
        dplyr::select(-Same.return.SNo)
      
      dis.sp=names(d1)
      dis.sp=dis.sp[-match("Target",dis.sp)]
      plot.list=vector('list',length(dis.sp))
      for(dd in 1:length(dis.sp))
      {
        plot.list[[dd]]=d1[,c('Target',dis.sp[dd])]%>%
          ggplot(aes_string(x='Target', y=dis.sp[dd]) ) +
          geom_bin2d(bins = 70) +scale_fill_continuous(type = "viridis") +
          theme_PA()
      }
      n.plts=n2mfrow(length(plot.list))
      multi.page <-ggarrange(plotlist=plot.list, nrow=n.plts[1],ncol = n.plts[2])
      ggsave(multi.page, filename = paste0(output_dir,"/Density.tiff"),width = 15, height = 15,compression="lzw")
      
    }
  }
  if(use.all.species.minlocs) minlocs.vec=rep(1,length(minlocs.vec))  
  
  #2. Plot catch by year for each species and extract relevant species
  Kept.species.Steph.Mac=vector('list',length(SP.list))
  names(Kept.species.Steph.Mac)=names(SP.list)
  for( i in 1:length(SP.list))
  {
    target.name=names(SP.list)[i]
    output_dir=paste(dir_plots,target.name,sep='/')
    if(!dir.exists(output_dir)) dir.create(output_dir)
    print(paste('Catch by year and non-collocated species for ----',target.name))
    Kept.species.Steph.Mac[[i]]=Species.catch.ranking(d=Data.daily.GN.DD%>%
                                                          filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[i]]$Same.return.SNo))%>%
                                                          filter(!SPECIES==9998),
                                                      TITL="Catch by year",
                                                      min.avrg.catch=Min.avrg.annual.catch, 
                                                      minyears=Minyears.with.catch,
                                                      minlocs=minlocs.vec[i],
                                                      target=SP.list[[i]],
                                                      drop.noncollocated=FALSE)
    print(Kept.species.Steph.Mac[[i]]$p)
    ggsave(paste0(output_dir,"/Top caught species.tiff"),width = 6,height = 12,compression = "lzw")
    
    p=fn.prop.by.shot(NMS=Kept.species.Steph.Mac[[i]]$dat%>%distinct(SNAME,SPECIES),
                      dd=Kept.species.Steph.Mac[[i]]$d_multi,
                      tar=SP.list[[i]])
    print(p)
    ggsave(paste0(output_dir,"/Top caught species_proportions.tiff"),width = 5,height = 12,compression = "lzw")
    
  }

  #3. Run Stephens and McCall
  tested.modls=c(1,2) #mod 2 has zone interactions
  Selected.model=data.frame(SNAME=names(SP.list))%>%
    mutate(Selected.model=ifelse(SNAME%in%c('Gummy Shark','Whiskery Shark',
                                            'Dusky Whaler','Sandbar Shark'),'Model_2','Model_1'))
  
  Stephens.McCall=vector('list',length(tested.modls))
  names(Stephens.McCall)=paste("Model",tested.modls,sep='_')
  system.time({for(m in 1:length(tested.modls))
  {
    MODL=tested.modls[m]
    Dummy=SP.list
    for(i in 1:length(SP.list))
    {
      target=SP.list[[i]]
      do.it=Selected.model[i,'Selected.model']
      if(paste0('Model_',MODL)==do.it) 
      {
        #Get inputs 
        d=Data.daily.GN.DD%>%
          filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[i]]$Same.return.SNo))%>%
          dplyr::select(Same.return.SNo,zone,LIVEWT.c,SPECIES,SNAME)%>%
          filter(!SPECIES==9998)
        EffrT=Effort.daily%>%
          filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[i]]$Same.return.SNo))%>%
          mutate(FishingSeason=as.numeric(substr(finyear,1,4)))%>%
          distinct(Same.return.SNo,Km.Gillnet.Hours.c,FishingSeason)
        
        
        #Set all non-co-located species to the same
        noncollocated=Kept.species.Steph.Mac[[i]]$noncollocated
        d_multi=Kept.species.Steph.Mac[[i]]$d_multi
        #bindata=Kept.species.Steph.Mac[[i]]$bindata
        indspecies_ids=Kept.species.Steph.Mac[[i]]$indspecies_ids 
        indspecies_names=Kept.species.Steph.Mac[[i]]$indspecies_names

        # add 'other' to group non-selected species and avoid 0 catch records of selected species
        id.other.sp=which(!names(d_multi)%in%c("Same.return.SNo","zone",target,indspecies_ids))
        d_multi$'99'  = rowSums(d_multi[, id.other.sp])
        indspecies_ids=c(indspecies_ids,99)
        indspecies_names=c(indspecies_names,"other fish")
        id.sp=which(!names(d_multi)%in%c("Same.return.SNo","zone"))
        bindata=d_multi[,id.sp]/d_multi[,id.sp]
        bindata[is.na(bindata)] = 0
        
        #remove grouped species
        d_multi=d_multi[,-which(!names(d_multi)%in%c("Same.return.SNo","zone",target,indspecies_ids,99))]
        bindata=bindata[,-which(!names(bindata)%in%c("Same.return.SNo","zone",target,indspecies_ids,99))]
       
        # if (length(noncollocated)>=1 & !is.na(noncollocated))
        # {
        #   dcols = which(names(d_multi)%in%noncollocated)
        #   if (length(dcols)>1)
        #   {
        #     d_multi$'99'  = rowSums(d_multi[, dcols])
        #   } else
        #   {
        #     d_multi$'99'=d_multi[, dcols]
        #   }
        #   bindata$'99' = d_multi$'99' / d_multi$'99'
        #   bindata$'99'[is.na(bindata$'99')] = 0
        #   d_multi=d_multi[,-dcols]
        #   bindata=bindata[,-which(names(bindata)%in%noncollocated)]
        # }
        
        names(bindata) = paste0("ind", names(bindata))
        target.name=names(SP.list)[i]
        
        print(paste('Apply Stephen & McCall method to------',target.name,'---model',MODL))
        
        output_dir=paste(dir_plots,target.name,sep='/')
 
        #combine effort, catch and bidata
        dat.comb=d_multi
        iidd=match(names(bindata),paste0("ind", names(dat.comb)))
        names(dat.comb)[iidd]=paste0('a_',names(dat.comb)[iidd])
        dat.comb=dat.comb%>%data.frame
        dat.comb$Catch=dat.comb[,match(paste0('a_',target),names(dat.comb))]
        dat.comb$Proportion=dat.comb$Catch/(rowSums(dat.comb[,iidd]))
        dat.comb=cbind(dat.comb,bindata)%>%
          left_join(EffrT%>%rename(Effort=Km.Gillnet.Hours.c),by='Same.return.SNo')%>%
          mutate(CPUE=Catch/Effort)
        
        #fit binomial model
        options(scipen=0)
        smfit=FitStephensMacCallModel(smdat=dat.comb, species_ids=target, indspecies_ids,indspecies_names,use_model=MODL)
        
        #reset predictions set to high but no catch other than 99
        id.dodgy=which(rowSums(smfit$data[,paste0('ind',c(target,subset(indspecies_ids,!indspecies_ids==99)))])==0 &
                         smfit$data[,paste0('ind',99)]>0)
        smfit$data[id.dodgy,'Pred']=0
        smfit$data[id.dodgy,'TargetSM']=0
        
        #plot model
        plot_width = 20;  plot_height = 20
        Plots=PlotStephensMacCallModel(smfit)
        names(Plots)=c('coefs.All_id','coefs.All_name','coefs.Sig_id','coefs.Sig_name')
        for(p in 1:length(Plots))
        {
          print(Plots[[p]])
          plot_name = paste0("filtermodel_", smfit$use_model, "_",names(Plots)[p],".png")
          ggsave(paste(output_dir, plot_name, sep="/"),width = plot_width, height = plot_height, scale=1, units = "cm")
        }
        
        #plot critical values
        Plots=PlotStephensMacCallCriticalValues(smfit, species_ids=target, species_names=target.name) 
        plot_width = 17; plot_height = 10
        names(Plots)=c('critical_values','predicted_probs','selected_proportions','selected_cpue','selected_cpue_without outliers')
        for(p in 1:length(Plots))
        {
          print(Plots[[p]])
          plot_name = paste0("filtermodel_", smfit$use_model, "_",names(Plots)[p],".png")
          ggsave(paste(output_dir, plot_name, sep="/"),width = plot_width, height = plot_height, scale=1, units = "cm")
        }
        
        #plot all cpue vs target cpue using delta lognormal
        Plots=PlotStephensMacCallTarget.vs.All_cpue(smfit, species_names=target.name)
        print(Plots)
        plot_name = paste0("filtermodel_", smfit$use_model, "_Delta_lognormal_cpue_All.vs.Taget.png")
        ggsave(paste(output_dir, plot_name, sep="/"),width = plot_width, height = plot_height, scale=1, units = "cm")
        
        #plot cpue vs pred prob
        p1=smfit$data%>%
          mutate(TargetSM=as.character(TargetSM))%>%
          ggplot(aes(CPUE,Pred,color = TargetSM))+
          geom_point()+
          coord_cartesian(xlim = c(0,quantile(smfit$data$CPUE,probs=0.95)))+
          theme_PA()
        print(p1)
        plot_name = paste0("CPUE vs pred prob_model",smfit$use_model, ".png")
        ggsave(paste(output_dir, plot_name, sep="/"),width = plot_width, height = plot_height, scale=1, units = "cm")
        
        p2=smfit$data%>%
          ggplot(aes(CPUE,Pred))+
          geom_bin2d(bins = 200)+
          coord_cartesian(xlim = c(0,quantile(smfit$data$CPUE,probs=0.95)))+
          geom_hline(yintercept=unique(smfit$data$critval), color = "red")+
          theme_PA()
        print(p2)
        plot_name = paste0("CPUE vs pred prob_density_model",smfit$use_model, ".png")
        ggsave(paste(output_dir, plot_name, sep="/"),width = plot_width, height = plot_height, scale=1, units = "cm")
        
        chek.dodgy=check.dodgy.predictions
        if(!target%in%c(18007,18003,17001,17003)) chek.dodgy=FALSE
        if(chek.dodgy)
        {
          pps=fn.untangle(dat=smfit$data,axs=6,tar.sp=target)
          print(pps$p1)
          plot_name = paste0("CPUE vs pred prob_High cpue but low(dodgy) and high(correct) pred. prob",
                             smfit$use_model, ".png")
          ggsave(paste(output_dir, plot_name, sep="/"),width = 8, height = 12, scale=1, units = "cm")
          
          if(!is.null(pps$p2))
          {
            print(pps$p2)
            plot_name = paste0("CPUE vs pred prob_Low cpue but high(dodgy) and low(correct) pred. prob",
                               smfit$use_model, ".png")
            ggsave(paste(output_dir, plot_name, sep="/"),width = 8, height = 12, scale=1, units = "cm")
          }
          rm(pps)
        }
        
        #Store predictions
        Dummy[[i]]=list(d_multi=d_multi,
                        dat.comb=dat.comb,
                        indspecies_ids=indspecies_ids,
                        indspecies_names=indspecies_names,
                        smfit.data=smfit$data%>%
                            dplyr::select(Same.return.SNo,Pred,TargetSM)%>%
                            rename(Step.MCal_target_prob=Pred,
                                   Step.MCal_target_group=TargetSM))
                          
        rm(smfit,Plots,indspecies_ids,indspecies_names,d_multi,bindata,noncollocated,dat.comb,p1,p2)
        
      }
     }
    Stephens.McCall[[m]]=Dummy
  }})  #takes 100 secs
  
  #4. Add to cpue stand data set
  for( i in 1:length(SP.list))
  {
    ID.mod=match(Selected.model[i,'Selected.model'],names(Stephens.McCall))
    DATA.list.LIVEWT.c.daily[[i]]=left_join(DATA.list.LIVEWT.c.daily[[i]],
                                            Stephens.McCall[[ID.mod]][[i]]$smfit.data,by=c("Same.return.SNo")) 
  }
  
  if(exists('Data.daily.GN.DD')) rm(Data.daily.GN.DD)
}

#Run cluster
if(do_cluster=="YES")
{
  #note: CLARA analysis as per Campbell et al 2017 on nfish as this has data at Sesssion level
  #       The CLARA (Clustering Large Applications) algorithm is an extension to the 
  #       PAM (Partitioning Around Medoids) clustering method for large data sets. It intends to 
  #       reduce the computation time in the case of large data set.
  Store.cluster=vector('list',length(SP.list)) 
  names(Store.cluster)=names(SP.list)
  
  #Using binomial data based on Stephen & MacCall
  N.clus=rep(2,length(SP.list))  #from initial optimum number
  if(Model.run=="First")out.clara=TRUE else out.clara=FALSE
  tic()
  for(i in 1:length(SP.list))
  {
    print(paste('Cluster analysis for ---------',names(SP.list)[i]))
    
    #run cluster
    ID.mod=match(Selected.model[i,'Selected.model'],names(Stephens.McCall))
    selected.species=c(SP.list[[i]],Stephens.McCall[[ID.mod]][[i]]$indspecies_ids)
    names(selected.species)=c(names(SP.list)[i],Stephens.McCall[[ID.mod]][[i]]$indspecies_names)
    a=Stephens.McCall[[ID.mod]][[i]]$d_multi%>%data.frame
    rownames(a)=a$Same.return.SNo
    colnames(a)=str_remove(colnames(a),"X")
    
    Store.cluster[[i]]=fn.cluster.Stephen.MacCall(a=a,
                                                  selected.species=selected.species,
                                                  n.clus=N.clus[i],
                                                  target=names(SP.list)[i],
                                                  check.clust.num=FALSE,
                                                  out.clara=out.clara,
                                                  apply.scale=FALSE,
                                                  do.proportion=TRUE)
    
    #add to data frame
    DATA.list.LIVEWT.c.daily[[i]]=left_join(DATA.list.LIVEWT.c.daily[[i]],
                                            Store.cluster[[i]],by=c("Same.return.SNo"))
  }
  toc()
  
  #Using catch rates of all species
  use.this.cluster=FALSE
  if(use.this.cluster)
  {
    #run cluster analysis
    N.clus=c(2,2,2,2)  #from initial optimum number
    for(i in 1:length(Tar.sp))
    {
      Store.cluster[[Tar.sp[i]]]=fn.cluster(data=Species.list.daily,
                                            TarSp=Tar.sp[i],
                                            n.clus=N.clus[i],
                                            target=names(Species.list.daily)[Tar.sp[i]],
                                            check.clustrbl="NO")
    }
    
    #add cluster to original data for use in standardisations
    for(s in Tar.sp)
    {
      DATA.list.LIVEWT.c.daily[[s]]=left_join(DATA.list.LIVEWT.c.daily[[s]]%>%
                                                mutate(Same.return.SNo.block10=paste(Same.return.SNo,block10)),
                                              Store.cluster[[s]],by=c("Same.return.SNo.block10"))%>%
        dplyr::select(-Same.return.SNo.block10)
    }
  }
  
  #Using aggregated catch by main species 
  Clus.vars=c("Catch.Target","Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar",
              "Catch.Groper","Catch.Snapper","Catch.Blue_mor")
  Tar.clus.vars=c("Catch.Whiskery","Catch.Gummy","Catch.Dusky","Catch.Sandbar")
  
  use.this.cluster=FALSE
  if(use.this.cluster)
  {
    n.clus=c(2,2,2,2)  #from initial optimum number
    scalem="YES"
    for(i in 1:length(Tar.sp))
    {
      Store.cluster[[Tar.sp[i]]]=fn.cluster(data=DATA.list.LIVEWT.c.daily,TarSp=Tar.sp[i],
                                            target=Tar.clus.vars[i],
                                            varS=Clus.vars,scaling=scalem,
                                            check.clustrbl="NO",n.clus=n.clus[i])
    }
    #add cluster to original data for use in standardisations
    for(s in Tar.sp)
    {
      DATA.list.LIVEWT.c.daily[[s]]=left_join(DATA.list.LIVEWT.c.daily[[s]],Store.cluster[[s]],by=c("Same.return.SNo"))
    }
  }
  
  rm(Store.cluster)
  
  #check cpue by cluster group
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster.boxplot",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,3,3,1),oma=c(2,1,.1,.5),las=1,mgp=c(2,.6,0))
  for(i in Tar.sp) boxplot((Catch.Target/Km.Gillnet.Hours.c)~cluster_clara,DATA.list.LIVEWT.c.daily[[i]],
                           main=names(DATA.list.LIVEWT.c.daily)[i],ylab="cpue",ylim=c(0,30))
  dev.off()
  
  for(s in 1:length(Tar.sp))
  {
    tiff(paste(HndL.Species_targeting,"Cluster/CLARA_cluster.boxplot_by_group_",Nms.sp[Tar.sp[s]],
               ".tiff",sep=''),width = 2400,
         height = 1800,units = "px", res = 300, compression = "lzw")    
    fn.compare.targeting(DAT=DATA.list.LIVEWT.c.daily[[Tar.sp[s]]],
                         Drop.var=Tar.clus.vars[s],
                         Title=Nms.sp[Tar.sp[s]])
    dev.off()
  }
}

#Run PCA
if(do_pca=="YES")
{
  #Winker et al 2014
  PercentExpl=90
  for(i in 1:length(Tar.sp))
  {
    target=paste("Proportion.",names(SP.list)[Tar.sp[i]],sep="")
    
    a=Species.list.daily[[Tar.sp[i]]]%>%
      filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[Tar.sp[i]]]$Same.return.SNo))%>%
      filter(!SPECIES%in%31000)%>%
      mutate(SPECIES=ifelse(SPECIES==19004,19000,SPECIES))
    TOP.ktch=a%>%group_by(SPECIES)%>%
      summarise(Catch=sum(LIVEWT.c))%>%
      arrange(-Catch)%>%
      mutate(Cum.ktch=cumsum(Catch),
             Quantiles=Cum.ktch/sum(Catch))
    top.sp=TOP.ktch$SPECIES[1:which.min(abs(TOP.ktch$Quantiles-95/100))]
    #a=a%>%mutate(SPECIES=ifelse(!SPECIES%in%top.sp,999999,SPECIES))%>%
    a=a%>%filter(SPECIES%in%top.sp)%>%
      group_by(Same.return.SNo,SPECIES)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T),
                Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c))%>%
      mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
      dplyr::select(SPECIES,cpue,Same.return.SNo)%>% 
      spread(SPECIES,cpue)%>%
      column_to_rownames(var = "Same.return.SNo")
    
    #proportion
    a[is.na(a)]=0
    a=a/rowSums(a,na.rm = T)
    
    #square root
    #a=a*a              
    
    #check correlation
    if(Model.run=="First")
    {
      M<-cor(a)
      fn.fig(paste(HndL.Species_targeting,"PCA/correlation_",target,sep=""),2400,2400)
      corrplot(M, method="circle")
      dev.off()
    }
    
    #run pca
    res.pca <- prcomp(a,center = T, scale. = F)    #it's already scaled
    
    #Visualize all eigenvalues, from most to least contribution
    if(Model.run=="First")
    {
      fn.fig(paste(HndL.Species_targeting,"PCA/eigenvalues_",target,sep=""),2400,2400)
      fviz_eig(res.pca)
      dev.off()
    }
    
    
    #Graph of each species proportion.
    #note: positive correlated variables point to the same side of the plot.
    #      negative correlated variables point to opposite sides of the grap
    if(Model.run=="First")
    {
      fn.fig(paste(HndL.Species_targeting,"PCA/PCA_species_",target,sep=""),2400,2400)
      fviz_pca_var(res.pca,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE )    # Avoid text overlapping
      dev.off()
    }
    
    
    #Acces the results
    Eig.val=get_eigenvalue(res.pca)
    #res.var <- get_pca_var(res.pca)
    #res.var$coord          # Coordinates
    #res.var$contrib        # Contributions to the PCs
    #res.var$cos2           # Quality of representation 
    PCA.variance=Eig.val
    keep.axis=which(PCA.variance$cumulative.variance.percent<PercentExpl)
    if(length(keep.axis)==1)keep.axis=1:2
    if(length(keep.axis)>3)keep.axis=keep.axis[1:3]
    PCA.scores=get_pca_ind(res.pca)$coord          # PCA Coordinates
    PCA.scores=data.frame(Same.return.SNo=row.names(PCA.scores),PCA.scores[,keep.axis])     #keep only relevant axes     
    row.names(PCA.scores)=NULL                                          
    PCA.axis=as.data.frame(res.pca$rotation)       # PCA axis loading. How species contribute to each axis
    
    #Plot variance explained
    if(Model.run=="First")
    {
      fn.fig(paste(HndL.Species_targeting,"PCA/Variance_",target,sep=""),2400,2400)
      plot(PCA.variance$cumulative.variance.percent,ylab="Cumulative variance explained",
           ylim=c(0,100),xlab="Axis")
      text(3,PCA.variance$cumulative.variance.percent[3],
           paste(round(PCA.variance$cumulative.variance.percent[3]),"%"),pos=3,col=2,srt=45)
      dev.off()
    }
    
    #Add pca axis to data
    DATA.list.LIVEWT.c.daily[[Tar.sp[i]]]=left_join(DATA.list.LIVEWT.c.daily[[Tar.sp[i]]],PCA.scores,by=c("Same.return.SNo"))
    
  }
  
  #check cpue by PCA dimension
  if(Model.run=="First")
  {
    fn.fig(paste(HndL.Species_targeting,"PCA/Cpue_dim.1",sep=""),2400,2400)
    par(mfcol=c(2,2),mar=c(1,3,3,1),oma=c(2,1,.1,.5),las=1,mgp=c(2,.6,0))
    for(s in Tar.sp)
    {
      with(DATA.list.LIVEWT.c.daily[[s]],plot((Catch.Target/Km.Gillnet.Hours.c)~Dim.1, 
                                              main=names(DATA.list.LIVEWT.c.daily)[s],ylab="cpue"))
      
    }
    dev.off()
  }
}

#Compare Stephens & McCall with cluster   
if(do_cluster=="YES" & do_Stephens_McCall=="YES")
{
  for(i in 1:length(SP.list))
  {
    p=DATA.list.LIVEWT.c.daily[[i]]%>%
      mutate(CPUE=Catch.Target/Km.Gillnet.Hours.c)%>%
      dplyr::select(CPUE,Step.MCal_target_prob,Step.MCal_target_group,cluster_clara)%>%
      rename(Step.MCal.prob=Step.MCal_target_prob,
             Step.MCal=Step.MCal_target_group,
             cluster=cluster_clara)

    p1=p%>%
      mutate(Step.MCal.prob=factor(round(Step.MCal.prob,1)),
             Step.MCal=paste('group',Step.MCal),
             cluster=paste('group',cluster))%>%
      gather(Method,Value,-CPUE)%>%
      ggplot(aes(Value,CPUE,fill=Method))+
      geom_boxplot()+
      ylim(0,quantile(p$CPUE,.99))+
      facet_wrap(~Method,scales='free_x')+
      theme_PA()+theme(legend.position = 'none')+xlab('')
    
    p.corr.dat=p%>%
      dplyr::select(Step.MCal,cluster)%>%
      filter(!is.na(Step.MCal) | !is.na(cluster))
    p2=ggcorrplot(p.corr.dat %>% 
                    correlation(),
                  lab = T,
                  show.diag = F)
    
    ggarrange(plotlist = list(p1,p2),ncol=1,nrow=2)

    plot_name = paste0("Cluster Vs Steph MacCal_", names(SP.list)[i],".tiff")
    ggsave(paste(HndL.Species_targeting, plot_name, sep="/"),width = 8,height = 6,compression = "lzw")
    

    
  }

}

#Remove Catch.x variables
for(i in 1:length(SP.list))
{
  drop=names(DATA.list.LIVEWT.c.daily[[i]])[grep('Catch.',names(DATA.list.LIVEWT.c.daily[[i]]))]
  drop=subset(drop,!drop%in%c("Catch.Target","Catch.Total"))
  if(length(drop)>0) DATA.list.LIVEWT.c.daily[[i]]=DATA.list.LIVEWT.c.daily[[i]]%>%dplyr::select(-all_of(drop))
}

# TABLE OF SENSITIVITY SCENARIOS ----------------------------------------------
Tab.Sensi=data.frame(Scenario=c("Base case","Nominal","All vessels & blocks","No efficiency"),
                     Standardisation=c("Yes","No",rep("Yes",2)),
                     Records_used=c("Reliable only","All",rep("Reliable only",2)),
                     Vessels_used=c("3 years","All","All","3 years"),
                     Blocks_used=c("3 years","All","All","3 years"),
                     Efficiency_increase=c(rep("Yes",3),"No"))

setwd(paste(getwd(),"/Outputs/Paper",sep=""))
# fn.word.table(WD=getwd(),TBL=Tab.Sensi,Doc.nm="Sensitivity tests",caption=NA,paragph=NA,
#               HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
#               Zebra='NO',Zebra.col='grey60',Grid.col='black',
#               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")



# COMPUTE RAW UNBALANCED INDEX ----------------------------------------------
Raw.index=vector('list',length(Species.list.raw)) 
names(Raw.index)=names(Species.list.raw)
Raw.index.daily=Raw.index
for(s in 1:length(Species.list.raw))
{
  #Monthly
  
  Raw.index[[s]]=Species.list.raw[[s]] %>%
    mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
    group_by(FINYEAR) %>%
    summarise(mean = mean(cpue),
              n = length(cpue),
              sd = sd(cpue)) %>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame%>%
    dplyr::select(-c(n,sd))
  
  #Daily
  Raw.index.daily[[s]]=Species.list.daily.raw[[s]] %>%
    mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
    group_by(FINYEAR) %>%
    summarise(mean = mean(cpue),
              n = length(cpue),
              sd = sd(cpue)) %>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame%>%
    dplyr::select(-c(n,sd))
}




# COMPUTE FOLLY AND NOMINAL INDICES FOR EXPORTING ----------------------------------------------

#--Effective
#note: As done by Rory, the Effective (ie Foly) index has all records (good  & bad reporters) from all blocks 
#       and vessels within effective area without 0 catches

DATA.list.LIVEWT.c_all_reporters=vector('list',length(SP.list)) 
names(DATA.list.LIVEWT.c_all_reporters)=names(SP.list)
DATA.list.LIVEWT.c.daily_all_reporters=DATA.list.LIVEWT.c_all_reporters

#monthly  
for ( i in Tar.sp)
{
  #create data sets 
  dummy=Effort.data.fun(Species.list[[i]],SP.list[[i]],"LIVEWT.c")  
  DATA.list.LIVEWT.c_all_reporters[[i]]=dummy$dat  
}

#daily 
for ( i in Tar.sp)
{
  #create data sets 
  dummy=Effort.data.fun.daily(Species.list.daily[[i]],SP.list[[i]],"LIVEWT.c",Aggregtn="SNo")
  DATA.list.LIVEWT.c.daily_all_reporters[[i]]=dummy$dat
}


#calculate Effective cpue (as per McAuley, ie use all records within effective area)
#note: effective is used as the conventionally used nominal cpue
Effective=vector('list',length(SP.list)) 
names(Effective)=names(SP.list)
Effective_daily=Effective
for(s in Tar.sp)   
{
  #Monthly
  Effective[[s]]=fn.out.effective(a=subset(Species.list[[s]],SPECIES==SP.list[[s]]))
  
  #Daily
  Effective_daily[[s]]=fn.out.effective(a=subset(Species.list.daily[[s]],SPECIES==SP.list[[s]]))
}


#-- Nominal   (not used)
#note: Ratio = mean(catch)/mean(effort)
#      Mean = mean(cpue)
#     LnMean= exp(mean(log(cpue))+bias corr)
#     DLnMean = exp(log(prop pos)+exp(mean(log(cpue))+bias corr)
Store_nom_cpues_monthly=vector('list',length(SP.list)) 
names(Store_nom_cpues_monthly)=names(SP.list)
Store_nom_cpues_daily=Store_nom_cpues_monthly
Hnd.ains=handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/Ainsline_different_cpues/")
suppressWarnings(for(s in nnn)
{
  #Monthly
  if(!is.null(DATA.list.LIVEWT.c[[s]]) & !is.null(BLKS.used[[s]]))
  {
    Store_nom_cpues_monthly[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c[[s]],Ktch.targt='catch.target',
                                            Effrt=c('km.gillnet.days.c','km.gillnet.hours.c'),
                                            explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
                                            cpue.units = c("kg/km gillnet day","kg/km gillnet hour"),spname=names(SP.list)[s],
                                            BLks=as.numeric(BLKS.used[[s]]),VesL=VES.used[[s]],Type="_monthly_")
  }
  
  #Daily weight
  if(!is.null(DATA.list.LIVEWT.c.daily[[s]])& !is.null(BLKS.used.daily[[s]]))
  {
    Store_nom_cpues_daily[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily[[s]],Ktch.targt='catch.target',
                                          Effrt=c('km.gillnet.days.c','km.gillnet.hours.c'),
                                          explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
                                          cpue.units = c("kg/km gillnet day","kg/km gillnet hour"),spname=names(SP.list)[s],
                                          BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Type="_daily_")
  }
})

# EVALUATE BALANCE OF DATA FOR QL ----------------------------------------------
#note: remove vessels with less than Min.Vess.yr (monthly) and Min.Vess.yr.d (daily)
#      for those vessels, keep blocks with more than Min.Vess.yr / Min.Vess.yr.d
if(!exists('BLKS.used.indi'))
{
  BLKS.used.indi=BLKS.used  #keep copy of indicators approach
  VES.used.indi=VES.used
  BLKS.used.daily.indi=BLKS.used.daily
  VES.used.daily.indi=VES.used.daily
}
hndl.kept=handl_OneDrive("Analyses/Catch and effort/Outputs/Kept_blocks_vessels/")
HndL=paste(hndl.kept,'QL_balanced_design/',sep="")
for(s in Tar.sp)
{
  #monthly
  a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=names(SP.list)[s],
                      what="monthly",MN.YR=Min.Vess.yr,pLot=T)
  BLKS.used[[s]]=a$this.blks
  VES.used[[s]]=a$this.ves
  
  #daily
  a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=names(SP.list)[s],
                      what="daily",MN.YR=Min.Vess.yr.d,pLot=T)
  BLKS.used.daily[[s]]=a$this.blks
  VES.used.daily[[s]]=a$this.ves
}

#export blocks and vessels used
for(s in nnn)
{
  write.csv(BLKS.used[[s]],paste(hndl.kept,"blocks_used_",names(SP.list)[s],"_monthly.csv",sep=""))
  write.csv(VES.used[[s]],paste(hndl.kept,"vessels_used_",names(SP.list)[s],"_monthly.csv",sep=""))
  write.csv(BLKS.used.daily[[s]],paste(hndl.kept,"blocks_used_",names(SP.list)[s],"_daily.csv",sep=""))
  write.csv(VES.used.daily[[s]],paste(hndl.kept,"vessels_used_",names(SP.list)[s],"_daily.csv",sep=""))
}

#show blocks kept
CEX=.85
SRt=45
fn.fig(paste(hndl.kept,"block_used_map_indicators",sep=""),2400, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in Tar.sp)
{
  #Monthly
  fn.show.blk(dat=BLKS.used[[s]],CEX=CEX,SRt=SRt,dat.all=BLKS.used.indi[[s]])
  if(s==Tar.sp[1])mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  fn.show.blk(dat=BLKS.used.daily[[s]],CEX=CEX,SRt=SRt,dat.all=BLKS.used.daily.indi[[s]])
  if(s==Tar.sp[1])mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",names(SP.list)[s],bty='n',cex=1.5)
}
mtext("Longitude",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Latitude",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()

fn.fig(paste(hndl.kept,"block_used_map_non_indicators",sep=""),2400, 2400)    
par(mfrow=n2mfrow(length(Non.Tar.sp)*2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in Non.Tar.sp) 
{
  #Monthly
  fn.show.blk(dat=BLKS.used[[s]],CEX=CEX,SRt=SRt,dat.all=BLKS.used.indi[[s]])
  if(s==min(Non.Tar.sp))mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  fn.show.blk(dat=BLKS.used.daily[[s]],CEX=CEX,SRt=SRt,dat.all=BLKS.used.daily.indi[[s]])
  if(s==min(Non.Tar.sp))mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",names(SP.list)[s],bty='n',cex=1.5)
}
mtext("Longitude",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Latitude",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Paper"))

# SHOW EFFECT OF USING km gn d VS km g h for gummy ----------------------------------------------
Show.gummy.hour=FALSE     
if(Show.gummy.hour)
{
  fn.fig("Appendix 2_Gummy_km.gn.d_VS_km.gn.h",1600,2400)
  par(mfrow=c(3,1),mar=c(1,3,2,1),oma=c(2,.5,.5,1.75),mgp=c(1.5,.5,0),cex.lab=1.5,las=1)
  Get.Mns(d=DATA.list.LIVEWT.c$"Gummy Shark",grp="FINYEAR",
          Vars=c("HOURS.c","SHOTS.c",'Km.Gillnet.Days.c','Km.Gillnet.Hours.c','Catch.Target'),
          LGND=c("Hours fished per day","Number of shots",
                 'km gillnet days','km gillnet hours','Catch (kg)'),
          add.arrow=c(rep(TRUE,5),FALSE))
  dev.off()
}

# EXPORT DATA FOR AINSLIE ----------------------------------------------
Exprt_Ainslie="NO"
if(Exprt_Ainslie=="YES")
{
  for(s in Tar.sp)
  {
    write.csv(DATA.list.LIVEWT.c[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c[s]),"_monthly.csv",sep=""),row.names=F)
    write.csv(DATA.list.LIVEWT.c.daily[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c.daily[s]),"_daily.csv",sep=""),row.names=F)
  }
}

# OUTPUT DATA TABLES ----------------------------------------------
#Output table with number of records available in effective area and numbers used in standardisation
do.this=FALSE  #need to update package Reporter not available anymore
if(do.this)
{
  TABle=vector('list',length(SP.list))
  names(TABle)=names(SP.list)
  for(s in Tar.sp)
  {
    #total records
    Tot.m=table(DATA.list.LIVEWT.c_all_reporters[[s]]$FINYEAR)
    Tot.d=table(DATA.list.LIVEWT.c.daily_all_reporters[[s]]$FINYEAR)
    
    Good.m=table(DATA.list.LIVEWT.c[[s]]$FINYEAR)
    Good.d=table(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR)
    
    #Used in Standardisation after selecting blocks and vessels
    DD=DATA.list.LIVEWT.c[[s]]
    DD=subset(DD,BLOCKX%in%as.numeric(BLKS.used[[s]]))      
    DD=subset(DD,VESSEL%in%VES.used[[s]])
    
    DD_daily=DATA.list.LIVEWT.c.daily[[s]]
    DD_daily=subset(DD_daily,BLOCKX%in%as.numeric(BLKS.used.daily[[s]]))      
    DD_daily=subset(DD_daily,VESSEL%in%VES.used.daily[[s]])
    
    Stand.m=table(DD$FINYEAR)
    Stand.d=table(DD_daily$FINYEAR)
    
    a=table(DD$FINYEAR,DD$VESSEL)
    a[a>0]=1
    a=rowSums(a)
    a[a<Threshold.n.vessls.per.yr]=NA
    Stand.m[which(is.na(a))]=NA
    
    a=table(DD_daily$FINYEAR,DD_daily$VESSEL)
    a[a>0]=1
    a=rowSums(a)
    a[a<Threshold.n.vessls.per.yr]=NA
    Stand.d[which(is.na(a))]=NA
    
    Tot.m=Tot.m[match(names(Good.m),names(Tot.m))]
    TaBs=data.frame(Year=c(names(Tot.m),names(Tot.d)),
                    Record=c(rep("Monthly",length(Tot.m)),rep("Daily",length(Tot.d))),
                    A_Total.num.rec.eff.area=c(Tot.m,Tot.d),
                    B_Good.num.rec.eff.area=round(c(Good.m/Tot.m,Good.d/Tot.d),2),
                    C_Used.fo.stand=round(c(Stand.m/Tot.m,Stand.d/Tot.d),2))
    
    names(TaBs)[3:ncol(TaBs)]=paste(names(TaBs)[3:ncol(TaBs)],names(DATA.list.LIVEWT.c)[s],sep="_")
    
    rm(DD,DD_daily)
    TABle[[s]]=TaBs
    
  }
  TABle=TABle[lengths(TABle) != 0]
  Table.nsamp <- Reduce(function(x, y) merge(x, y, all=T,by=c("Year", "Record")), TABle, accumulate=F)
  Table.nsamp=Table.nsamp[,c("Year", "Record",sort(names(Table.nsamp[3:ncol(Table.nsamp)])))]
  Export.tbl(WD=getwd(),Tbl=Table.nsamp,Doc.nm="Sample_sizes",caption=NA,paragph=NA,
             HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
             Zebra='NO',Zebra.col='grey60',Grid.col='black',
             Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
             HDR.names=c('Year', 'Record','Total number of records within eff. area',
                         'Proportion reliable','Proportion used in standardisation'),
             HDR.span=c(1,1,N.species,N.species,N.species),
             HDR.2nd=c("","",rep(c("Dusky","Gummy","Sandbar","Whiskery"),3)))
}


# EXPLORATORY ANALYSES -----------------------------------------------------------------------
#note: This applies cede() exploration, and checks for outliers in catch and effort
#         max monthly ktch, effort (~ 40 tonnes, ~ 5800 km gn h (@ 30 days X 24 h X 8000 m), respectively) 
#         max trip (daily kg) ktch, effort (~ 15 tonnes, ~ 1900 km gn h (@ 10 days X 24 h X 8000 m), respectively) 
#      Also does GAM exploratory stuff

#Temperature
explore.Oceanographic=FALSE
if(explore.Oceanographic)
{
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/Oceanographic variables/Temperature.pdf'))
  SST%>%
    group_by(Long,Lat,year)%>%
    summarise(Temperature=mean(Temperature))%>%
    ggplot(aes(x = Long, y = Lat)) +
    labs(x = "Long", y = "Lat", fill = "Value") +
    geom_raster(aes(fill = Temperature))+
    facet_wrap(vars(year))

  SST%>%
    filter(Lat<(-29) & Lat>=(-31))%>%
    group_by(year,month)%>%
    summarise(mean=mean(Temperature,na.rm=T),
              sd=sd(Temperature,na.rm=T),
              n=length(Temperature))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n))%>%
    ggplot(aes(x = year, y = mean)) +
    geom_point()+
    geom_smooth(method = "lm")+
    geom_errorbar(aes(ymin=lowCL, ymax=uppCL), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Year", y = "Temperature (+/- 95% CI)") +
    ggtitle("Mean temperature 29 to 31  S")+
    facet_wrap(vars(month))
  
  SST%>%
    filter(Lat<(-31))%>%
    group_by(year,month)%>%
    summarise(mean=mean(Temperature,na.rm=T),
              sd=sd(Temperature,na.rm=T),
              n=length(Temperature))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n))%>%
    ggplot(aes(x = year, y = mean)) +
    geom_point()+
    geom_smooth(method = "lm")+
    geom_errorbar(aes(ymin=lowCL, ymax=uppCL), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Year", y = "Temperature (+/- 95% CI)") +
    ggtitle("Mean temperature South of 31  S")+
    facet_wrap(vars(month))
  
  SST%>%
    filter(Lat<(-29) & month%in%c(6:9))%>%
    group_by(year)%>%
    summarise(mean=mean(Temperature,na.rm=T),
              sd=sd(Temperature,na.rm=T),
              n=length(Temperature))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n))%>%
    ggplot(aes(x = year, y = mean)) +
    geom_point()+
    geom_smooth(method = "lm")+
    geom_errorbar(aes(ymin=lowCL, ymax=uppCL), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Year", y = "Temperature (+/- 95% CI)") +
    ggtitle("June to September.  Mean temperature south of 29 S")
  
  SST%>%
    filter(Lat<(-29) & month%in%c(12,1:4))%>%
    group_by(year)%>%
    summarise(mean=mean(Temperature,na.rm=T),
              sd=sd(Temperature,na.rm=T),
              n=length(Temperature))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n))%>%
    ggplot(aes(x = year, y = mean)) +
    geom_point()+
    geom_smooth(method = "lm")+
    geom_errorbar(aes(ymin=lowCL, ymax=uppCL), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Year", y = "Temperature (+/- 95% CI)") +
    ggtitle("December to April.  Mean temperature south of 29 S")
  
  SST%>%
    filter(Lat<(-29) & Lat>=(-31))%>%
    mutate(Temperature=Temperature-mean(Temperature,na.rm=T))%>%
    group_by(year)%>%
    summarise(Anomaly=mean(Temperature,na.rm=T))%>%
    mutate(Sign=as.factor(ifelse(Anomaly<0,'neg','pos')))%>%
    ggplot(aes(x = year, y = Anomaly, fill = Sign)) +
    geom_bar(stat = "identity", show.legend = FALSE)  +
    ggtitle("Temperature anomalies.  29 to 31  S")+
    scale_fill_manual("legend", values = c("neg" = "blue", "pos" = "red"))
  
  SST%>%
    filter(Lat<(-31))%>%
    mutate(Temperature=Temperature-mean(Temperature,na.rm=T))%>%
    group_by(year)%>%
    summarise(Anomaly=mean(Temperature,na.rm=T))%>%
    mutate(Sign=as.factor(ifelse(Anomaly<0,'neg','pos')))%>%
    ggplot(aes(x = year, y = Anomaly, fill = Sign)) +
    geom_bar(stat = "identity", show.legend = FALSE)  +
    ggtitle("Temperature anomalies.  South of 31  S")+
    scale_fill_manual("legend", values = c("neg" = "blue", "pos" = "red"))
  dev.off()
}
#SOI
if(explore.Oceanographic)
{
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/Oceanographic variables/SOI.pdf'))
  SOI%>%
    mutate(year=Year+Month/13,
           Sign=as.factor(ifelse(SOI<0,'neg','pos')))%>%
    ggplot(aes(x = year, y = SOI, fill = Sign)) +
    geom_bar(stat = "identity", show.legend = FALSE)  +
    ggtitle("SOI")+
    scale_fill_manual("legend", values = c("neg" = "blue", "pos" = "red"))+
    geom_hline(yintercept=-7,linetype="dashed", color = "blue")+
    annotate("text", x = min(SOI$Year)+.75, y = -8, label = "El Niño", size = 5)+
    geom_hline(yintercept=7,linetype="dashed", color = "red")+
    annotate("text", x =min(SOI$Year)+.75 , y = 8, label = "La Niña", size = 5)
  dev.off()
}
#Freo
if(explore.Oceanographic)
{
  pdf(handl_OneDrive('Analyses/Catch and effort/Outputs/Exploratory/Oceanographic variables/Freo.pdf'))
  Freo%>%
    mutate(year=Year+Month/13,
           Sign=as.factor(ifelse(Freo<0,'neg','pos')))%>%
    ggplot(aes(x = year, y = Freo, fill = Sign)) +
    geom_bar(stat = "identity", show.legend = FALSE)  +
    ggtitle("Freo")+
    scale_fill_manual("legend", values = c("neg" = "blue", "pos" = "red"))
  dev.off()
}
if(do.Exploratory=="YES")
{
  hndl.expl=handl_OneDrive("Analyses/Catch and effort/Outputs/Exploratory/")
  
  # Standard exploratory analyses
  system.time({
    for(s in Tar.sp)
    {
      #monthly
      pdf(paste(hndl.expl,names(SP.list)[s],"_monthly.pdf",sep="")) 
      if(!is.null(Store_nom_cpues_monthly[[s]]))
      {
        dummy=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
        dummy=subset(dummy,blockx%in%BLKS.used[[s]])
        fn.expl.cede(d=dummy,PREDS=Predictors_monthly,kg=TRUE,Do.ggplts=F)
      }
      if(!is.null(DATA.list.LIVEWT.c[[s]])) fn.box.plt.year(d=subset(DATA.list.LIVEWT.c[[s]],Catch.Target>0))
      if(!is.null(Store_nom_cpues_monthly[[s]])) check.cpue(Store_nom_cpues_monthly[[s]]$QL_dat,paste("monthly",names(DATA.list.LIVEWT.c)[s]),4)
      if(names(SP.list)[s]=="Whiskery Shark") fn.pred.effect(DATA=DATA.list.LIVEWT.c[[s]],PREDS=Predictors_monthly)
      dev.off()
      
      #daily
      pdf(paste(hndl.expl,names(SP.list)[s],"_daily.pdf",sep=""))
      if(!is.null(Store_nom_cpues_daily[[s]]))
      {
        dummy=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
        dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
        fn.expl.cede(d=dummy,PREDS=Predictors_daily,kg=TRUE,Do.ggplts=T)
      }
      if(!is.null(DATA.list.LIVEWT.c.daily[[s]])) fn.box.plt.year(d=subset(DATA.list.LIVEWT.c.daily[[s]],Catch.Target>0))
      if(!is.null(Store_nom_cpues_daily[[s]])) check.cpue(Store_nom_cpues_daily[[s]]$QL_dat,paste("daily",names(DATA.list.LIVEWT.c)[s]),4) 
      if(names(SP.list)[s]=="Whiskery Shark") fn.pred.effect(DATA=DATA.list.LIVEWT.c.daily[[s]],PREDS=subset(Predictors_daily,!Predictors_daily%in%c('mesh','nlines.c')))
      dev.off()
    }
  })    #takes 2 mins
  
  # GAM exploratory analyses 
  #   notes on GAM:  
  #       Test significance for smooth terms: if you cannot draw a horizontal 
  #           line through the 95% confidence interval the smooth term is significant.
  #       Good GAM tutorial https://www.youtube.com/watch?v=q4_t8jXcQgc
  #                         https://noamross.github.io/gams-in-r-course/
  
  # gam.check() for checking model fit:
  #   . Q-Q plot, which compares the model residuals to a normal distribution. 
  #   A well-fit model's residuals will be close to a straight line. 
  # 
  #   . On bottom left is a histogram of residuals. We would expect this to have a symmetrical bell shape. 
  # 
  #   . On top-right is a plot of residual values. These should be evenly distributed around zero. 
  # 
  #   . Finally, on the bottom-right is plot of response against fitted values. 
  #   A perfect model would form a straight line. 
  #   We don't expect a perfect model, but we do expect the pattern to cluster around the 1-to-1 line.
  
  
  # concurvity() for checking colinearity among terms:
  #   . The first mode, full = TRUE, reports overall concurvity for each smooth. 
  #   Specifically, it shows how much each smooth is predetermined by all the other smooths.
  #   Since concurvity is complex, the function reports three different ways of measuring concurvity. 
  #   Each is better in some situations. 
  #   What is important is that you should always look at the worst case, 
  #   and if the value is high (say, over 0.8), inspect your model more carefully.
  # 
  #   . If any of these values from the full = TRUE mode is high, we will want to also use 
  #   the second mode, setting full = FALSE. With full = FALSE
  #   These show the degree to which each variable is predetermined by each other variable, 
  #   rather than all the other variables. This can be used to pinpoint which variables 
  #   have a close relationship.
  
  system.time({
    for(s in Tar.sp)
    {
      NM=names(DATA.list.LIVEWT.c.daily)[s]
      
      #Daily
      print(paste("------",NM,"-daily"))
      Covars=c('YEAR.c','MONTH','LAT10.corner','LONG10.corner',
               'Temperature','Temp.res','Freo','SOI','Lunar','nlines.c',
               'mesh','Mean.depth','Dim.1','Dim.2','Dim.3')
      if(do_pca=="YES") Covars=c(Covars,Targeting.vars)
      
      fn.explr.gam.rel(d=DATA.list.LIVEWT.c.daily[[s]],
                       Fktrs=c('FINYEAR','month.as.factor','VESSEL','shots.c'),
                       Covars=Covars,
                       OUT=paste(names(DATA.list.LIVEWT.c.daily)[s],'_daily',sep=''))
      
      #Monthly
      print(paste("------",NM,"-monthly"))
      fn.explr.gam.rel(d=DATA.list.LIVEWT.c[[s]],
                       Fktrs=c('FINYEAR','month.as.factor','VESSEL','shots.c'),
                       Covars=c('YEAR.c','MONTH','LAT10.corner','LONG10.corner',
                                'Temperature','Temp.res','Freo','SOI'),
                       OUT=paste(names(DATA.list.LIVEWT.c)[s],'_monthly',sep=''))
    }
  })    #takes 9 mins
  
}

# EFFICIENCY_ESTIMATE EFFORT CREEP----------------------------------------------
if(get.efficiency.creep)
{
  #Display ves.vars thru time
  if(Model.run=="First")
  {
    for(i in Tar.sp)
    {
      NM=names(DATA.list.LIVEWT.c)[i]
      print(paste("Display vessel carachteristics thru time for ------",NM))
      fun.check.ves.char.on.cpue(d=DATA.list.LIVEWT.c[[i]],NM=NM)
    }
  }
  
  #other vars are highly correlated or very small sample size
  ves.vars.selected=c("ENGDERAT","ENGNUM","ENGPOWR","FLYBRIDGE","GROSSTON",
                      "HULLCONS","HULLNUMB","HULLTYPE","MAXBEAM","WHEELHOU",
                      "GPS","COECHO","PLOT","RADAR")
  
  #Determine representativeness of vessels filling in vessel survey
  VSL.surveyed=vector('list',length(nnn))
  names(VSL.surveyed)=names(Species.list)
  dumi=dumi1=VSL.surveyed
  for(i in nnn)
  {
    if(i%in%Tar.sp)
    {
      NM=names(DATA.list.LIVEWT.c)[i]
      print(paste("How representative surveyed vessels are for ------",NM))
      x=fun.check.ves.char.representative(d=DATA.list.LIVEWT.c[[i]],NM)
      dumi[[i]]=x$p
      dumi1[[i]]=x$p1
      VSL.surveyed[[i]]=x$vesl.list
    }
  }
  
  dumi<-dumi[!sapply(dumi,is.null)]
  ggarrange(plotlist = dumi,ncol=2,nrow=2)
  ggsave(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/Vessels with vessel char_annual prop. total effort.tiff')),
         width=8,height= 6,compression="lzw")
  dumi1<-dumi1[!sapply(dumi1,is.null)]
  ggarrange(plotlist = dumi1,ncol=2,nrow=2,common.legend = TRUE)
  ggsave(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/Vessels with vessel char_cum total effort.tiff')),
         width=8,height= 6,compression="lzw")

  Boat.names_surveyed=Data.daily.GN%>%
    filter(VESSEL%in%unique(unlist(VSL.surveyed)))%>%
    distinct(VESSEL, BoatName,MastersName)%>%
    arrange(MastersName,BoatName)
  
  #Determine technology adoption trend  
  #fun.technology.adoption(d=Vessel.charac.exp%>%rename(VESSEL=BOATREGO),dis.var=c("GPS","PLOT","COECHO","RADAR"))
  
  
  
  #Model efficiency creep

}


#Explore vessel and vessel char effect on selected records for standardisation
if(Model.run=="First")
{
  for(s in Tar.sp)
  {
    NM=names(DATA.list.LIVEWT.c)[s]
    
    #Monthly
    iid=match(SpiSis[match(NM,names(SpiSis))],names(First.year.catch))
    Firs.yr.ktch=names(First.year.catch[[iid]])
    theseyears=sort(unique(DATA.list.LIVEWT.c[[s]]$FINYEAR))
    theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
    
    d=DATA.list.LIVEWT.c[[s]]%>%
      filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]] & FINYEAR%in%theseyears )%>%
      mutate(CPUE=Catch.Target/Km.Gillnet.Hours.c)
    
    #remove first years with very few positive catch observation
    if(NM=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
    if(NM=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
    x=fn.show.cpue.vsl(d)
    x$p
    ggsave(paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by vessel_monthly.tiff'))),
           width = 10,height = 8, dpi = 300, compression = "lzw")
    write.csv(x$d,paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by vessel_monthly.csv'))),
              row.names = F)
    
    ggarrange(plotlist=list(x$p1,x$p2,x$p3,x$p4,x$p5,x$p6))
    ggsave(paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by char_monthly.tiff'))),
           width = 10,height = 8, dpi = 300, compression = "lzw")
    
    
    #Daily
    iid=match(SpiSis[match(names(DATA.list.LIVEWT.c.daily)[s],names(SpiSis))],names(First.year.catch.daily))
    Firs.yr.ktch=names(First.year.catch.daily[[iid]])
    theseyears=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR))
    theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
    
    #remove first year of transition from Monthly to Daily returns due to effort misreporting
    drop.this=match("2006-07",theseyears)
    if(!is.na(drop.this))theseyears=theseyears[-drop.this]
    
    d=DATA.list.LIVEWT.c.daily[[s]]%>%
      filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]] & FINYEAR%in%theseyears)%>%
      mutate(CPUE=Catch.Target/Km.Gillnet.Hours.c)
    x=fn.show.cpue.vsl(d)
    x$p
    ggsave(paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by vessel_daily.tiff'))),
           width = 10,height = 8, dpi = 300, compression = "lzw")
    write.csv(x$d,paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by vessel_daily.csv'))),
              row.names = F)
    
    ggarrange(plotlist=list(x$p1,x$p2,x$p3,x$p4,x$p5,x$p6))
    ggsave(paste0(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Exploratory/check vessel and chars/',NM,'_cpue by char_daily.tiff'))),
           width = 10,height = 8, dpi = 300, compression = "lzw")
    
    
  }
}


# CONSTRUCT STANDARDISED ABUNDANCE INDEX----------------------------------------------
#note: standard run takes ~ 10 secs per species
source(handl_OneDrive('Analyses/Catch and effort/Git_catch.and.effort/CPUE Construct standardised abundance index.R'))

# INFLUENCE PLOTS ---------------------------------------------------------
# Bentley et al 2012; not useful for delta-MC method. check https://github.com/trophia/influ/blob/master/influ.R
if(do.influence=="YES")
{
  if(Use.Delta)
  {
    HnDll=handl_OneDrive("Analyses/Catch and effort/Outputs/Influence.plot/")
    Store.Influence=vector('list',length(SP.list)) 
    names(Store.Influence)=names(SP.list)
    Store.Influence.daily=Store.Influence
    for(s in Tar.sp)
    {
      #Monthly
      Terms=all.vars(Best.Model[[s]])[-1]
      Terms=Terms[-match("finyear",Terms)]
      Term.Type=Terms
      Term.Type=ifelse(Term.Type%in%Categorical,"CAT","Cont")
      names(Term.Type)=Terms
      Term.Type=subset(Term.Type,Term.Type=="CAT")
      Terms=names(Term.Type)
      pdf(paste(HnDll,Nms.sp[s],".monthly.CDI.pdf",sep=""))
      Store.Influence[[s]]=Influence.fn(MOD=Stand.out[[s]]$Pos$res,DAT=Stand.out[[s]]$Pos$DATA,
                                        Term.type=Term.Type,termS=Terms,add.Influence="YES",SCALER=4)
      dev.off()
      
      
      #Daily
      Terms=all.vars(Best.Model.daily.gam_delta[[s]]$pos)[-1]
      Terms=Terms[-match("finyear",Terms)]
      Term.Type=Terms
      Term.Type=ifelse(Term.Type%in%Categorical,"CAT","Cont")
      names(Term.Type)=Terms
      Term.Type=subset(Term.Type,Term.Type=="CAT")
      Terms=names(Term.Type)
      Term.Type=Term.Type[!is.na(Term.Type)]
      pdf(paste(HnDll,Nms.sp[s],".daily.CDI.pdf",sep=""))
      Store.Influence.daily[[s]]=Influence.fn(MOD=Stand.out.daily[[s]]$Pos$res.gam,
                                              DAT=Stand.out.daily[[s]]$Pos$DATA,Term.type=Term.Type,termS=Terms,
                                              add.Influence="YES",SCALER=4)
      dev.off()
    } 
    
    Store.Influence=Store.Influence[Tar.sp] 
    Store.Influence.daily=Store.Influence.daily[Tar.sp] 
    Over.inf.monthly=do.call(rbind,lapply(Store.Influence, '[[', match('Over.all.influence',names(Store.Influence[[1]]))))
    myList=vector('list',length(Store.Influence.daily))
    for(l in 1:length(myList))myList[[l]]=data.frame(t(Store.Influence.daily[[l]]$Over.all.influence))
    Over.inf.daily=as.data.frame(rbindlist(myList, fill = TRUE, use.names=TRUE))
    rownames(Over.inf.daily)=names(Store.Influence.daily)
    Over.inf.monthly=round(100*Over.inf.monthly,2)
    Over.inf.daily=round(100*Over.inf.daily,2)
    
    Over.inf.percent=as.data.frame(rbindlist(list(data.frame(Species=rownames(Over.inf.monthly),dat="Monthly",Over.inf.monthly),
                                                  data.frame(Species=rownames(Over.inf.daily),dat="Daily",Over.inf.daily)),
                                             fill = TRUE, use.names=TRUE))
    fn.word.table(WD=getwd(),TBL=Over.inf.percent,Doc.nm="Influence_table_percent",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
    
    #Compare influence of all terms
    LWD=3
    LTY.col=c("black","grey45","grey20","grey85","grey55","grey70")
    Whr=c("topleft","topleft","topleft","topright")
    Whr.d1=c("bottomleft","bottomleft","bottomleft","topright")
    Whr.d2=c("topright","topright","topright","bottom")
    
    
    fn.fig("Figure 5.All.terms.Influence",2400, 2400)
    par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
    for ( s in 1:length(Tar.sp))
    {
      #Monthly
      dumy=unlist(lapply(Store.Influence, "[[", match("Annual.Dev",names(Store.Influence[[1]]))))
      YLIM=c(min(dumy),max(dumy))
      Compare.term.infl.fun(A=Store.Influence[[s]],Whr[s],WHERE2=NA,YLIM=YLIM,spliT="NO")
      if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      #Daily
      dumy=unlist(lapply(Store.Influence.daily, "[[", match("Annual.Dev",names(Store.Influence.daily[[1]]))))
      YLIM=c(min(dumy),max(dumy))
      Compare.term.infl.fun(A=Store.Influence.daily[[s]],WHERE=Whr.d1[s],WHERE2=Whr.d2[s],YLIM=YLIM,spliT="YES")
      if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      
      mtext(Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Influence",side=2,line=1.25,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
    #Show each species monthly CDI
    for ( s in 1:length(Tar.sp))
    {
      fn.fig(paste("Figure 6.",Nms.sp[Tar.sp[s]],"_CDI",sep=''),2400, 1200)
      par(mfcol=c(2,3),mar=c(1,1,.5,1.25),oma=c(2,4,1,.2),las=1,mgp=c(1.9,.5,0))
      Fig.CDI.paper.fn(store=Store.Influence[[s]]$store,
                       SCALER=2.5,termS=c("vessel","blockx","month"))
      mtext("Financial year                   Coefficient    ",
            side=2,line=2.8,cex=1,las=3,outer=T)
      dev.off()
    }
    
  }
}

# COMPARE ALL DIFFERENT WAYS OF CALCULATING CPUES -----------------------------------------------------------------------
if(Model.run=="First")
{
  dis=c("Whiskery Shark","Gummy Shark","Dusky Whaler","Sandbar Shark")
  for(d in 1:length(dis))
  {
    NM=dis[d]
    fn.plot.all.indices(sp=NM)
    ggsave(handl_OneDrive(paste('Analyses/Catch and effort/Outputs/Compare all index types/',
                                NM,'.jpeg',sep='')),width=10,height= 10)  
    
  }
  do.dis=FALSE
  if(do.dis) source(handl_OneDrive('Analyses/Catch and effort/Git_catch.and.effort/compare cpue series.R'))
  
}  

# CONSTRUCT SPATIAL STANDARDISED INDEX ---------------------------------------------------------
#note: standard run takes ~ 7 secs per species
if(do.spatial.cpiui) source(handl_OneDrive('Analyses/Catch and effort/Git_catch.and.effort/CPUE Construct standardised abundance index_spatial.R'))
on.exit(stopCluster(cl))
#ACA

# EXPORT INDICES -----------------------------------------------------------------------
setwd(handl_OneDrive("Analyses/Data_outs"))
Sel.vars=c("finyear","response","CV","lower.CL","upper.CL","SE")
nams.Sel.vars=c("Finyear","Mean","CV","LOW.CI","UP.CI","SE")

Nms.sp=ifelse(Nms.sp=="Bronze whaler","Copper shark",Nms.sp)
Nms.sp=ifelse(Nms.sp=="Hammerhead sharks","Hammerheads",Nms.sp)

#1. Target species 

#1.1. Absolute scale 

  #zones combined with NO efficiency creep
for (s in Tar.sp)
{
  a=subset(Pred[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.monthly_no.creep.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.daily_no.creep.csv",sep=""),row.names=F) 
  
  rm(a)
}
  #zones combined with efficiency creep
for (s in Tar.sp)
{
  a=subset(Pred.creep[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.monthly.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily.creep[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.daily.csv",sep=""),row.names=F) 
  rm(a)
}
  #Spatial
if(do.spatial.cpiui)
{
  #by zones with NO efficiency creep
  for (s in 1:length(Tar.sp))
  {
    Zn=names(Zone_preds.monthly[[s]])
    for(z in 1:length(Zn))
    {
      #Standardised
      a=subset(Zone_preds.monthly[[s]][[z]]$Preds,select=Sel.vars)   
      names(a)=nams.Sel.vars
      ii=Tar.sp[s]
      if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.monthly.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
      
      a=subset(Zone_preds.daily[[s]][[z]]$Preds,select=Sel.vars)    
      names(a)=nams.Sel.vars
      write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.daily.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
      
      rm(a)
    }
  }
  
  #by zones with efficiency creep
  for (s in 1:length(Tar.sp))
  {
    Zn=names(Zone_preds.monthly[[s]])
    for(z in 1:length(Zn))
    {
      #Standardised
      a=subset(Zone_preds.monthly[[s]][[z]]$Preds.creep,select=Sel.vars)   
      names(a)=nams.Sel.vars
      ii=Tar.sp[s]
      if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.monthly",Zn[z],".csv",sep=""),row.names=F) 
      
      a=subset(Zone_preds.daily[[s]][[z]]$Preds.creep,select=Sel.vars)   
      names(a)=nams.Sel.vars
      write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.daily",Zn[z],".csv",sep=""),row.names=F) 
      
      rm(a)
    }
  }
}

#1.2.Relative scale   

  #zones combined with efficiency creep
for (s in Tar.sp)
{
  a=subset(Pred.normlzd[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.monthly_relative.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily.normlzd[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(getwd(),'/',Nms.sp[s],'/',Nms.sp[s],".annual.abundance.basecase.daily_relative.csv",sep=""),row.names=F) 
  rm(a)
}
  #spatial
if(do.spatial.cpiui)
{
  #4.22.12.6 by zones with efficiency creep
  for (s in 1:length(Tar.sp))
  {
    Zn=names(Zone_preds.monthly[[s]])
    for(z in 1:length(Zn))
    {
      a=subset(Zone_preds.monthly[[s]][[z]]$Preds.nrmlzd,select=Sel.vars)   
      names(a)=nams.Sel.vars
      ii=Tar.sp[s]
      if(xport_monthly) write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.monthly.",Zn[z],"_relative.csv",sep=""),row.names=F) 
      
      a=subset(Zone_preds.daily[[s]][[z]]$Preds.nrmlzd,select=Sel.vars)   
      names(a)=nams.Sel.vars
      write.csv(a,paste(getwd(),'/',Nms.sp[ii],'/',Nms.sp[ii],".annual.abundance.basecase.daily.",Zn[z],"_relative.csv",sep=""),row.names=F) 
      
      rm(a)
    }
  }
}


#2. Other species zones combined with efficiency creep  
Sel.vars.other=c("finyear","response","CV","lower.CL","upper.CL","SE")

  #2.1. Absolute scale with creep
for (s in nnn[-sort(Tar.sp)])
{
  nmx=Nms.sp[s]
  if(nmx=="Wobbegong") nmx="Wobbegongs"
  if(nmx=="Common sawshark") nmx="Sawsharks"
  if(!is.null(Pred.creep[[s]]))
  {
    a=subset(Pred.creep[[s]],select=Sel.vars.other)   
    names(a)=nams.Sel.vars
    if(xport_monthly) write.csv(a,paste(getwd(),'/',nmx,'/',nmx,".annual.abundance.basecase.monthly.csv",sep=""),row.names=F) 
  }
  
  if(!is.null(Pred.daily.creep[[s]]))
  {
    a=subset(Pred.daily.creep[[s]],select=Sel.vars.other)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(getwd(),'/',nmx,'/',nmx,".annual.abundance.basecase.daily.csv",sep=""),row.names=F) 
  }
  rm(a,nmx)
}
  #2.2. Relative scale
for (s in nnn[-sort(Tar.sp)])
{
  nmx=Nms.sp[s]
  if(nmx=="Wobbegong") nmx="Wobbegongs"
  if(nmx=="Common sawshark") nmx="Sawsharks"
  if(!is.null(Pred.normlzd[[s]]))
  {
    a=subset(Pred.normlzd[[s]],select=Sel.vars.other)   
    names(a)=nams.Sel.vars
    if(xport_monthly) write.csv(a,paste(getwd(),'/',nmx,'/',nmx,".annual.abundance.basecase.monthly_relative.csv",sep=""),row.names=F) 
  }
  
  if(!is.null(Pred.daily.normlzd[[s]]))
  {
    a=subset(Pred.daily.normlzd[[s]],select=Sel.vars.other)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(getwd(),'/',nmx,'/',nmx,".annual.abundance.basecase.daily_relative.csv",sep=""),row.names=F) 
  }
  rm(a,nmx)
}


# EXPLORE SMOOTH HAMMERHEAD ---------------------------------------------------------
if(reset.hammerheads.to.smooth.hh)
{
  Explore.smooth.HH=FALSE
  if(Explore.smooth.HH)
  {
    hndl.ZmuHH=handl_OneDrive('Analyses/Catch and effort/Outputs/Paper/smooth HH_explore/')
    # Explore under reporting smooth HH following CITES listing by some vessel 
    a=DATA.list.LIVEWT.c.daily$`Smooth Hammerhead Shark`
    X=a%>%
      group_by(FINYEAR,VESSEL)%>%
      summarise(Catch.Target=sum(Catch.Target))%>%
      data.frame%>%
      spread(FINYEAR,Catch.Target,fill=0)%>%
      mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE) )%>%
      filter(Total>500)%>%
      dplyr::select(-Total)
    X1=X%>%
      gather(FINYEAR,Catch,-VESSEL)%>%
      mutate(yr=as.numeric(substr(FINYEAR,1,4)))
    
    Vs.lvl= X1%>%
      group_by(VESSEL)%>%
      summarise(Sum=sum(Catch))%>%
      arrange(-Sum)
    
    X1=X1%>%mutate(VESSEL=factor(VESSEL,levels=Vs.lvl$VESSEL))
    X1%>%
      ggplot(aes(x=yr, y=VESSEL)) +
      geom_tile(aes(fill =  Catch), colour = "grey50") +
      scale_fill_viridis(direction = -1,guide = guide_colorbar()) + 
      theme_classic(base_size = 28, base_line_size = 0.5) +
      theme(axis.ticks.length=unit(-0.25, "cm"),
            axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
            panel.grid.minor = element_blank(),panel.spacing.x=unit(-1, "cm"),panel.spacing.y=unit(-1, "cm"),
            strip.background = element_blank(), 
            legend.spacing.x = unit(1, "pt"), 
            legend.key.width  = unit(5, "cm"), 
            legend.position = "bottom" ,
            legend.title = element_blank(), 
            legend.key.height = unit(1, 'cm'),
            plot.margin=unit(c(0,0,0,0), "cm"))
    ggsave(paste(hndl.ZmuHH,"Smooth HH daily catch by vessel.tiff"),width = 12,height = 12, dpi = 300, compression = "lzw")
    
    
    s=match("Smooth Hammerhead Shark",names(DATA.list.LIVEWT.c.daily)) 
    fn.bubl=function(Tab1)
    {
      N=nrow(Tab1)
      plot(1:nrow(Tab1),1:nrow(Tab1),col="transparent",xlim=c(1,ncol(Tab1)),yaxt='n',xaxt='n',xlab='',ylab='')
      for(n in 1:ncol(Tab1)) points(rep(n,N),1:N,cex=fn.scale (Tab1[,n],3),pch=19,col="steelblue")
      axis(1,1:ncol(Tab1),colnames(Tab1),las=3)
      axis(2,1:nrow(Tab1),rownames(Tab1),las=1)
    }
    d=DATA.list.LIVEWT.c.daily[[s]]%>%
      filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
    colnames(d)=tolower(colnames(d))
    
    
    
    
    fn.fig(paste(hndl.ZmuHH,"annual records distribution_vessel"),2000, 2400)   
    fn.bubl(Tab1=table(d$vessel,d$finyear))
    dev.off()
    
    fn.fig(paste(hndl.ZmuHH,"annual records distribution_month"),2000, 2400)
    fn.bubl(Tab1=table(d$month,d$finyear))
    dev.off()
    
    fn.fig(paste(hndl.ZmuHH,"annual records distribution_block10"),2000, 2400)
    fn.bubl(Tab1=table(d$block10,d$finyear))
    dev.off()
    
    Terms=Predictors_daily
    Continuous=Covariates.daily
    
    Terms=tolower(Terms)
    Continuous=tolower(Continuous)
    Factors=Terms[!Terms%in%Continuous]
    Terms=all.vars(Best.Model.daily[[s]])[-1]
    d <- d%>%
      dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
      mutate(cpue=catch.target/km.gillnet.hours.c)
    d <- makecategorical(Factors[Factors%in%Terms],d)
    
    #1. Fit model
    mod<-bam(Best.Model.daily[[s]],data=d,family='tw',method="fREML",discrete=TRUE) 
    
    #2.Make predictions
    lsmeas=pred.fun(mod=mod,biascor="NO",PRED="finyear")
    
    
    #nominal
    out.mean = d %>%
      group_by(finyear) %>%
      summarise(mean = mean(cpue),
                n = length(cpue),
                sd = sd(cpue)) %>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      data.frame
    
    out.delta.log = d %>%
      mutate(cpue.target=cpue)%>%
      group_by(finyear) %>%
      summarise(n = length(cpue.target),
                m = length(cpue.target[cpue.target>0]),
                mean.lognz = mean(log(cpue.target[cpue.target>0])),
                sd.lognz = sd(log(cpue.target[cpue.target>0]))) %>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
      as.data.frame
    
    
    NRMLIZD="YES"
    Y=as.numeric(substr(lsmeas$finyear,1,4))  
    fn.fig(paste(hndl.ZmuHH,"index"),2000, 2400)
    with(lsmeas,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(response)
           if(NRMLIZD=="NO")
           {
             YMAX=max(response+1.96*SE/Mn)
             YMIN=min(response-1.96*SE/Mn)
           }
           
           if(NRMLIZD=="YES")
           {
             YMAX=2
             YMIN=0
           }
           
           plot(Y-.1,response/Mn,col="orange",pch=19,cex=1.5,lwd=1.5,type='o',main=Nms.sp[s],xlab='financial year',
                ylab='relative cpue',ylim=c(YMIN,YMAX))
           segments(Y-.1,(response+1.96*SE)/Mn,
                    Y-.1,(response-1.96*SE)/Mn,col="orange")
         })
    with(out.mean,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y,mean/Mn,col="grey60",pch=19,type='o')
           segments(Y,(lowCL)/Mn,
                    Y,(uppCL)/Mn,col="grey60")
         })
    with(out.delta.log,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y+.15,mean/Mn,col="grey80",pch=19,type='o')
           segments(Y+.15,(lowCL)/Mn,
                    Y+.15,(uppCL)/Mn,col="grey80")
         })
    legend('topright',c("lsmeans","nominal.mean","nominal.delta.log"),
           lty=1,col=c('orange',"grey50","grey80"),bty='n')
    points(2015,1,pch="*",col=2,cex=3)
    text(2015,1,"CITES listing",col=2,srt=45,pos=4)
    dev.off()
    
    library(RColorBrewer)
    
    #Explore cpue by vessel
    d=DATA.list.LIVEWT.c.daily[[s]]%>%
      filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
    colnames(d)=tolower(colnames(d))
    Uni.ves=unique(d$vessel)
    
    ves.cpue=vector('list',length(Uni.ves))
    names(ves.cpue)=Uni.ves
    for(v in 1:length(Uni.ves))
    {
      ves.cpue[[v]] = d %>%
        filter(vessel==Uni.ves[v])%>%
        mutate(cpue.target=catch.target/km.gillnet.hours.c)%>%
        group_by(finyear) %>%
        summarise(n = length(cpue.target),
                  m = length(cpue.target[cpue.target>0]),
                  mean.lognz = mean(log(cpue.target[cpue.target>0])),
                  sd.lognz = sd(log(cpue.target[cpue.target>0]))) %>%
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
        dplyr::select(finyear,mean,lowCL,uppCL)
      
    }
    
    n <- length(Uni.ves)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    CL = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    NRMLIZD="YES"
    Y=sort(as.numeric(substr(unique(d$finyear),1,4)))  
    
    fn.fig(paste(hndl.ZmuHH,"cpue by vessel"),2000, 2400)   
    plot(Y,Y,ylim=c(0,3),ylab="CPUE",xlab="Financial year",col="transparent")
    for(v in 1:length(Uni.ves))with(ves.cpue[[v]],
                                    {
                                      yy=sort(as.numeric(substr(ves.cpue[[v]]$finyear,1,4)))
                                      if(NRMLIZD=="NO")Mn=1
                                      if(NRMLIZD=="YES")Mn=mean(mean)
                                      points(yy+.15,mean/Mn,col=CL[v],pch=19,type='o',cex=1.25)
                                      #segments(yy+.15,(lowCL)/Mn,
                                      #          yy+.15,(uppCL)/Mn,col=CL[v])
                                    })
    legend("topleft",Uni.ves,fill=CL,cex=.95,bty='n',ncol=4)
    dev.off()
    
    fn.fig(paste(hndl.ZmuHH,"prop records reporting smoothHH by vessel"),2000, 2400)   
    d2=d%>%mutate(catch.target=ifelse(catch.target>0,1,0))%>%select(vessel,finyear,catch.target)
    Yrs=as.numeric(substr(sort(unique(d2$finyear)),1,4))
    d2=tbl_df(d2)%>%
      group_by(vessel, finyear) %>%
      summarise (n = n(),
                 count=sum(catch.target)) %>%
      mutate(prop = count / n)%>%
      select(-c( n,count))%>%
      spread(finyear,prop,fill=0)%>%
      data.frame
    Tab1=as.matrix(d2[,-1])
    colnames(Tab1)=Yrs
    rownames(Tab1)=d2$vessel
    
    N=nrow(Tab1)
    plot(1:nrow(Tab1),1:nrow(Tab1),col="transparent",xlim=c(1,ncol(Tab1)),yaxt='n',xaxt='n',xlab='',ylab='')
    for(n in 1:ncol(Tab1)) points(rep(n,N),1:N,cex=Tab1[,n]*4,pch=19,col="steelblue")
    axis(1,1:ncol(Tab1),colnames(Tab1),las=3)
    axis(2,1:nrow(Tab1),rownames(Tab1),las=1)
    dev.off()
    
    
    #Refit model without dodgy vessels
    Dodgy=c("E 035")
    Uni.ves=Uni.ves[!Uni.ves%in%Dodgy]
    
    d1 <- d%>%
      filter(vessel%in%Uni.ves)%>%
      dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
      mutate(cpue=catch.target/km.gillnet.hours.c)
    d1 <- makecategorical(Factors[Factors%in%Terms],d1)
    
    #1. Fit model
    
    mod<-bam(cpue ~ finyear + s(vessel, bs = "re") +s(month, k = 12, 
                                                      bs = "cc") + s(long10.corner, lat10.corner) + s(mean.depth) ,data=d1,family='tw',method="fREML",discrete=TRUE) 
    
    #2.Make predictions
    lsmeas=pred.fun(mod=mod,biascor="NO",PRED="finyear")
    
    
    #nominal
    out.mean = d1 %>%
      group_by(finyear) %>%
      summarise(mean = mean(cpue),
                n = length(cpue),
                sd = sd(cpue)) %>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      data.frame
    
    out.delta.log = d1 %>%
      mutate(cpue.target=cpue)%>%
      group_by(finyear) %>%
      summarise(n = length(cpue.target),
                m = length(cpue.target[cpue.target>0]),
                mean.lognz = mean(log(cpue.target[cpue.target>0])),
                sd.lognz = sd(log(cpue.target[cpue.target>0]))) %>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
      as.data.frame
    
    
    NRMLIZD="YES"
    Y=as.numeric(substr(lsmeas$finyear,1,4))  
    fn.fig(paste(hndl.ZmuHH,"index_without_E035"),2000, 2400)
    with(lsmeas,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(response)
           if(NRMLIZD=="NO")
           {
             YMAX=max(response+1.96*SE/Mn)
             YMIN=min(response-1.96*SE/Mn)
           }
           
           if(NRMLIZD=="YES")
           {
             YMAX=2
             YMIN=0
           }
           
           plot(Y-.1,response/Mn,col="orange",pch=19,cex=1.5,lwd=1.5,type='o',main=Nms.sp[s],xlab='financial year',
                ylab='relative cpue',ylim=c(YMIN,YMAX))
           segments(Y-.1,(response+1.96*SE)/Mn,
                    Y-.1,(response-1.96*SE)/Mn,col="orange")
         })
    with(out.mean,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y,mean/Mn,col="grey60",pch=19,type='o')
           segments(Y,(lowCL)/Mn,
                    Y,(uppCL)/Mn,col="grey60")
         })
    with(out.delta.log,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y+.15,mean/Mn,col="grey80",pch=19,type='o')
           segments(Y+.15,(lowCL)/Mn,
                    Y+.15,(uppCL)/Mn,col="grey80")
         })
    legend('topright',c("lsmeans","nominal.mean","nominal.delta.log"),
           lty=1,col=c('orange',"grey50","grey80"),bty='n')
    points(2015,1,pch="*",col=2,cex=3)
    text(2015,1,"CITES listing",col=2,srt=45,pos=4)
    dev.off()
    
    
  }
}


# EXPLORE DUSKY AND SANDBAR UNCERTAINTY IN DAILY DATA--------------------------------------------
Explore.why.dusky.sandbar.uncertain=FALSE
if(Explore.why.dusky.sandbar.uncertain)
{
  dusk.san=names(DATA.list.LIVEWT.c.daily)
  dusk.san=match(c("Dusky Whaler","Sandbar Shark"),dusk.san)
  for(s in dusk.san)
  {
    d=DATA.list.LIVEWT.c.daily[[s]]%>%
      filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]] & FINYEAR%in%theseyears)
    Terms=Predictors_daily
    Continuous=Covariates.daily
    colnames(d)=tolower(colnames(d))
    Terms=tolower(Terms)
    Continuous=tolower(Continuous)
    Factors=Terms[!Terms%in%Continuous]
    Terms=all.vars(Best.Model.daily[[s]])[-1]
    d <- d%>%
      dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
      mutate(cpue=catch.target/km.gillnet.hours.c)
    d <- makecategorical(Factors[Factors%in%Terms],d)
    
    
    #Best.Model.daily[[s]]
    What="finyear+ vessel + month + depth + lat long"
    mod<-bam(cpue ~ finyear+ s(vessel, bs = "re")+ 
               s(month, k = 12,bs = "cc")+ s(mean.depth)+ s(long10.corner, lat10.corner),data=d,family='tw',method="fREML",discrete=TRUE)
    
    #2.Make predictions
    lsmeas=pred.fun(mod=mod,biascor="NO",PRED="finyear")
    
    
    #nominal
    out.mean = d %>%
      group_by(finyear) %>%
      summarise(mean = mean(cpue),
                n = length(cpue),
                sd = sd(cpue)) %>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      data.frame
    
    out.delta.log = d %>%
      mutate(cpue.target=cpue)%>%
      group_by(finyear) %>%
      summarise(n = length(cpue.target),
                m = length(cpue.target[cpue.target>0]),
                mean.lognz = mean(log(cpue.target[cpue.target>0])),
                sd.lognz = sd(log(cpue.target[cpue.target>0]))) %>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
      as.data.frame
    
    
    NRMLIZD="YES"
    Y=as.numeric(substr(lsmeas$finyear,1,4))  
    
    hndl.crap=handl_OneDrive('Analyses/Catch and effort/Outputs/Paper/why dusky and sandbar daily so uncertain/')
    fn.fig(paste(hndl.crap,Nms.sp[s],"_",What,sep=""),2000, 2400)
    with(lsmeas,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(response)
           if(NRMLIZD=="NO")
           {
             YMAX=max(response+1.96*SE/Mn)
             YMIN=min(response-1.96*SE/Mn)
           }
           
           if(NRMLIZD=="YES")
           {
             YMAX=2
             YMIN=0
           }
           
           plot(Y-.1,response/Mn,col="orange",pch=19,cex=1.5,lwd=1.5,type='o',main=Nms.sp[s],xlab='financial year',
                ylab='relative cpue',ylim=c(YMIN,YMAX))
           segments(Y-.1,(response+1.96*SE)/Mn,
                    Y-.1,(response-1.96*SE)/Mn,col="orange")
         })
    with(out.mean,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y,mean/Mn,col="grey60",pch=19,type='o')
           segments(Y,(lowCL)/Mn,
                    Y,(uppCL)/Mn,col="grey60")
         })
    with(out.delta.log,
         {
           if(NRMLIZD=="NO")Mn=1
           if(NRMLIZD=="YES")Mn=mean(mean)
           points(Y+.15,mean/Mn,col="grey80",pch=19,type='o')
           segments(Y+.15,(lowCL)/Mn,
                    Y+.15,(uppCL)/Mn,col="grey80")
         })
    legend('topright',c("lsmeans","nominal.mean","nominal.delta.log"),
           lty=1,col=c('orange',"grey50","grey80"),bty='n')
    dev.off()
  }
  
  
}




# Create CPUE stand paper figures -----------------------------------------------------------------------
if (plot.cpue.paper.figures=="YES")
{
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Paper"))
  source(handl_OneDrive('Analyses/Catch and effort/Git_catch.and.effort/CPUE stand paper figures.r'))
}