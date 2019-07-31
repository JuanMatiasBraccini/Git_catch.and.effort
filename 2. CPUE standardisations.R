#--------- CPUE STANDARDISATIONS OF DIFFERENT DATASETS ---------#

#NOTE:  This script standardises the catch and effort data for the 4 commercial shark species,

# missing: Vessel Characteristics from annual survey and Questionnaires (see 'Vessel characteristics from DoF's annual survey')

#       To update SOI and Mean Freo Sealevel each year, run "Get.SOI.Freo.R" 
#       To update Temperature, run "SST.r"


#Index:  #----1. DATA SECTION-----#  
#           1.1 Import data
#           1.2 Control what parts of script are activated
#	            1.2.1 Data controls
#	            1.2.2 Procedure controls
#	            1.2.3 Reporting controls

#----2. PARAMETERS SECTION-----#

#----3. FUNCTIONS SECTION-----#

#----4. PROCEDURE SECTION-----#
#           4.1 Deal with zone1-zone2 Boundary blocks to a zone
#           4.2 Extract number of vessels per species range
#           4.3 Data fixes
#           4.4 Remove NA effort
#           4.5 Create useful vars
#           4.6 Proportion of dusky and copper shark
#           4.7 Define indicative vessels and blocks 
#           4.8 Put data into a list
#           4.9 Compare nominal all records VS 'good reporters' only
#           4.10 Construct wide database for analysis
#           4.12 Drop first years of sandbar data because vessels don't meet selection 
#           4.13 Corroborate effective area 
#           4.14 Identify targeting behaviour
#           4.15 Table of sensitivity scenarios
#           4.16 Compute foly and nominal index for exporting
#           4.17 Evaluate balance of data subset based on QL
#           4.18 Show gummy monthly cpue effect of using km gn d or km g h
#           4.19 Export data to ainslie
#           4.20 Output data tables
#           4.21 Check outliers in catch and effort for removing nonsense values 
#           4.22 Construct index of abundance 
#               4.22.1 Explore data used for standardisation
#               4.22.2 Show applied effort creep
#               4.22.3 Select model structure
#               4.22.4 Run standardisation
#               4.22.5 Export deviance explained
#               4.22.6 Run sensitivity tests
#               4.22.7 Fit diagnostics
#               4.22.8 Plot base case and nominal 
#               4.22.9 Influence plots 
#               4.22.10 Show month and block effects
#               4.22.11 Construct spatial standardised catch rates
#               4.22.12 Export catch rates

#----5. REPORT SECTION-----#


rm(list=ls(all=TRUE))

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
library(cede)     #Malcolm Haddon's
library(gridExtra)
library(glmulti)  #model selection
library(fitdistrplus)  #select distribution
library(coefplot)  #visualise coefficients
library(emmeans)  #for model predictions
library(doParallel)
library(tibble)
library(cluster)
library(factoextra) #for plotting
library(mgcv)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)   



setwd("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other")
source("Compare.error.structure.R")
source("Deviance.explained.R")
source("Sorting.objects.R")
source("MS.Office.outputs.R")
setwd("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics")
source("fn.fig.R")
source("Nominal_cpue_functions.R")
source("C:\\Matias\\R\\Sort ls by size.R")



##############--- 1. DATA SECTION ---###################

setwd('C:/Matias/Analyses/Catch and effort/Data_outs')
Data.daily.original=read.csv("Data.daily.original.csv")
Data.monthly.GN=read.csv("Data.monthly.GN.csv")
Data.daily.GN=read.csv("Data.daily.GN.csv")
Effort.daily=read.csv("Effort.daily.csv")
Effort.monthly=read.csv("Effort.monthly.csv")
Mesh.monthly=read.csv("Mesh.monthly.csv")
Mesh.size=read.csv("Mesh.size.csv")

lst <- strsplit(Data.daily.GN$Same.return.SNo, "\\s+")
Data.daily.GN$SNo <- sapply(lst ,'[', 1)
Data.daily.GN$DSNo <- sapply(lst, '[', 2)
Data.daily.GN$TSNo <- sapply(lst, '[', 3)
rm(lst)

#Block10 locations
BlOCK_10=read.csv("C:/Matias/Data/Mapping/Blocks_10NM.csv")
names(BlOCK_10)=c("block10","LAT","LONG")
Metro_BlOCK_10=subset(BlOCK_10, LAT>(-33) & LAT<=(-31) & LONG<116)

#Southern Oscillation Index
SOI=read.csv("C:/Matias/Data/Oceanography/SOI.csv")

#Mean Freo sea level
Freo=read.csv("C:/Matias/Data/Oceanography/Freo_mean_sea_level.csv")  

#SST
SST=read.csv("C:/Matias/Data/Oceanography/SST.csv") 

#Fishable areas       
Depth.range="species_specific"
#Depth.range=200
if(Depth.range==200)
{
  Whis.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_whiskery_200.csv")
  Gum.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_gummy_200.csv")
  Dusky.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_dusky_200.csv")
  Sand.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_sandbar_200.csv")  
}
if(Depth.range=="species_specific")
{
  Grab.Area="C:/Matias/Data/Catch and Effort/FishableArea/"
  Whis.fishArea=read.csv(paste(Grab.Area,"BLOCKX_whiskery_30_70.csv",sep=""))
  Gum.fishArea=read.csv(paste(Grab.Area,"BLOCKX_gummy_les_70.csv",sep=""))
  Dusky.fishArea=read.csv(paste(Grab.Area,"BLOCKX_dusky_less_60.csv",sep=""))
  Sand.fishArea=read.csv(paste(Grab.Area,"BLOCKX_sandbar_30_120.csv",sep=""))  
  
  Whis.fishArea_b10=read.csv(paste(Grab.Area,"block10.whiskery.csv",sep=""))
  Gum.fishArea_b10=read.csv(paste(Grab.Area,"block10.gummy.csv",sep=""))
  Dusky.fishArea_b10=read.csv(paste(Grab.Area,"block10.dusky.csv",sep=""))
  Sand.fishArea_b10=read.csv(paste(Grab.Area,"block10.sandbar.csv",sep=""))  
  
  
}


#Vessel characteristics from DoF's annual survey                MISSING
#note: must update ""VesselGearSurveyData.csv" in C:\Matias\Data\Fishing power
#source("C:/Matias/Analyses/Catch and effort/Git_catch_and_effort/5.Vessel_characteristics.R")


#Vessel characteristics from fisher questionnaire    MISSING, incomplete survey
#read("C:\\Matias\\Data\\Fishing power\\Questionnaire responses\\Questionaire.xlsx")




#1.2 Control what parts of script are activated

#1.2.1 Data controls

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


#1.2.2 Procedure controls 

#Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.jpeg="YES"
Do.tiff="NO"



#Control type of model run
#Model.run="First"    # for first time of doing standardisation. This allows selection of 
#     best model and sensitivity tests
Model.run="Standard"  # for running standardisations in subsequent years

#Define is testing model structure (i.e. selection of model terms)
Def.mod.Str="NO"

#Control if continue exploring fit
improve="NO"

#Control if checking correlation between catch covariates
Chck.Corr="NO"

#Control if doing cluster analysis of daily data
do_cluster="NO"

#Control if doing PCA analysis of daily data
do.PCA="NO"

#Control if comparing lognormal to gamma
Compare.best.lognormal.gamma="NO"

#control if checking interactions
check.interactions="NO"

#control if combining vessel categories
Combine.ves="NO"

#Control if comparing cpue and catch as reponse vars, effort as offset, vessel effect, etc
COMPARE.RAW.etc="NO"
if(Model.run=="First") COMPARE.RAW.etc="YES"


#Control if doing sensitivity tests
if(Model.run=="First") do.sensitivity="YES"
if(Model.run=="Standard") do.sensitivity="NO"

#Control if aggregating daily
do.aggregated.daily="NO"


#Control if comparing best model from catch vs from cpue as response vars
Compare.best.catch.cpue="NO"

#Control if doing k-fold cross validation
Do.K_n.fold.test="NO"
#Do.K_n.fold.test="YES"

#Control if fitting to catch or cpue
Fit.to.what="catch"
#Fit.to.what="cpue"

#Control if comparing glm and gam spatial predictions
compare.glm.gam="NO"

#Control if extracting glm deviance table
if(Model.run=="Standard") Extract.Deviance.table="NO" 
if(Model.run=="First") Extract.Deviance.table="YES" 

#Control if showing binomial fit
do.binomial.fit="NO" 



#Control if adding interactions to model
With.interact="YES"
#control criteria for selecting indicative vessels
#Criteria.indi='all'  #criteria for selecting sensitivity of indicative vessel
Criteria.indi='subset'  #use different criteria for indicative vessels

#control criteria for indicative vessels
second.criteria="top.percent"

WHICH.VESSEL="top.vess.across.yr"   #select top vessels across years
#WHICH.VESSEL="top.vess.by.yr"      #select top vessels per year

#Control if sorting factor levels 
#note: this specifies the glm reference level of block and vessel
#Sort.levels="NO"
#Sort.levels="Habitat.area"
#Sort.levels="Most.common"
Sort.levels="Highest.catch" 

#Control new data factors
Sel.lev='most.common'



#Control if doing AMM actions
do.actions.AMM.2017="NO"

#1.2.3 Reporting controls

#Control if doing influence plots (Bentley et al 2012)
if(do.sensitivity=="NO") do.influence="NO"
if(do.sensitivity=="YES") do.influence='YES'

#control color of plots
#do.colors="YES"
do.colors="NO"  #grey scale

if(do.colors=="YES") what.color="cols"
if(do.colors=="NO") what.color="black"

#Control if doing exploratory analyses
do.Exploratory="NO"

#Control of showing relative cpue
Show.relative.index="YES"

#Control what index to export
#Export.relative="NO"
Export.relative="YES"


##############--- 2. PARAMETERS SECTION ---###################

Stand.eff=1000          #express all standardised catches as: catch per 1000 km gn d
Lg.Efrt=log(Stand.eff)

Red.th=1  #percentage reduction in deviance to accept term

#Criteria for keeping species for analysis
N.keep=5      #in years
Min.kg=100   #in kg

#Species definitions
Shark.species=5001:24900
Indicator.sp=c(17001,17003,18003,18001,18007)
Greynurse.protection='1999-00'
TARGETS.name=c("SHARK, WHISKERY","SHARK, GUMMY","SHARK, BRONZE WHALER","SHARK, THICKSKIN (SANDBAR)")
TARGETS=list(17003,17001,c(18003,18001),18007)
names(TARGETS)=TARGETS.name
N.species=length(TARGETS)
SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
SPvec=c(17003,17001,18003,18007)
names(SPvec)=SPECIES.vec

#Zn1-Zn2 boundary blocks
Boundary.Blks=c(3416,3516,3616)  

#Minimum annual catch of target species (in kgs)
MIN.ktch=100 

#Indicative blocks
  # minimum number of years with observations for each block
MIN.obs.BLK=5  
MIN.obs.BLK.sens=0  
MIN.obs.BL10K=5 # minimum number of observation per block10


#Indicative vessels 
  #minimum number of years reporting the species
Threshold.n.yrs=5
Threshold.n.yrs.monthly=5 
Threshold.n.yrs.daily=5 
Threshold.n.yrs.sens=0  #sensitivity
Threshold.n.vessls.per.yr=3  #keep years with a least 3 different vessels


#qualification levels
QL_expl_ktch_prop=.9   #proportion of explained annual catch for selected record
PRP.MLCLM=0.1          #proportion of catch of target


#Qualification levels minimum number of years with positive records
Min.Vess.yr=5 #monthly
Min.Vess.yr.d=5 #daily

#Cpue correction for assumed increase in fishing power prior to 1995
# Rory McAuley: 2% annual (i.e. 2%, 4%, 6%, etc) to 1995, then flat
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
Min.FL.w=88;Min.FL.g=43; Min.FL.d=55; Min.FL.s=47

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


##############--- 3. FUNCTIONS SECTION ---###################

fn.scale=function(x,scaler) ((x/max(x,na.rm=T))^0.5)*scaler

#functions for creating species data lists
fn.cpue.data=function(Dat,EffrrT,sp)
{
  TAB=with(subset(Dat,SPECIES%in%sp),unique(YEAR.c))
  if(length(TAB)>=N.keep)
  {
    #aggregate records by Same return (drop issues with Condition and Bioregion...)
    Dat=Dat%>%group_by(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX,Boundary.blk,SPECIES,SNAME,YEAR.c,
                       LAT,LONG,Same.return,TYPE.DATA,zone,Reporter,Sch.or.DogS,
                       Temperature,Temp.res,SOI,Freo,Freo_lag6,Freo_lag12)%>%
      summarise(LIVEWT = sum(LIVEWT),
                LIVEWT.c = sum(LIVEWT.c))%>%
      data.frame()
    
    #add effort
    Ids=match(c("LAT","LONG"),names(EffrrT))
    Dat=Dat%>%left_join(EffrrT[,-Ids],by=c("BLOCKX","FINYEAR","MONTH","VESSEL"))
    
    #consider effor reporter
    Dat$Reporter=with(Dat,ifelse(Eff.Reporter=="bad","bad",Reporter))  
    
    #remove records with NA effort
    Dat=subset(Dat,!(is.na(Km.Gillnet.Days.c) | is.na(Km.Gillnet.Hours.c)))
    
    #remove school shark or dogfish shots if not the target species
    idd=which(sp%in%c(17008,20000))
    if(length(idd)==0) Dat=subset(Dat,!(Sch.or.DogS=="Yes"))
    
    #keep records from first year with data as some species (e.g. Sandbars) 
    # didn't have a code for reporting early on
    TAB=with(subset(Dat,SPECIES%in%sp),table(YEAR.c))
    Dat=subset(Dat,YEAR.c>=as.numeric(names(TAB[1])))
    
    #for greynurse, drop years post protection
    idd=which(sp%in%8001)
    if(length(idd)==1 & length(sp)<3)
    {
      Greyn.yrs=sort(unique(Dat$FINYEAR))
      Greyn.yrs=Greyn.yrs[1:match(Greynurse.protection,Greyn.yrs)]
      Dat=subset(Dat,FINYEAR%in%Greyn.yrs)
    }
  }else
  {
    Dat=NULL
  }
  return(Dat)
}
fn.cpue.data.daily=function(Dat,EffrrT,sp)
{
  TAB=with(subset(Dat,SPECIES%in%sp),unique(YEAR.c))
  if(length(TAB)>=N.keep)
  {
    Dat=Dat%>%group_by(date,Same.return.SNo,TSNo,Same.return,FINYEAR,MONTH,
                       VESSEL,METHOD,BLOCKX,block10,SPECIES,SNAME,YEAR.c,LAT,LONG,
                       TYPE.DATA,zone,Reporter,Sch.or.DogS,ZnID,
                       Temperature,Temp.res,SOI,Freo,Freo_lag6,Freo_lag12,Lunar)%>%
      summarise(LIVEWT = sum(LIVEWT),
                LIVEWT.c = sum(LIVEWT.c),
                nfish = sum(nfish))%>%
      data.frame()
    
    #add effort
    Ids=match(c("BLOCKX","FINYEAR","MONTH","LAT","LONG","VESSEL","block10"),names(EffrrT))
    Dat=Dat%>%left_join(EffrrT[,-Ids],by=c("Same.return.SNo"))
    
    #consider effor reporter
    Dat$Reporter=with(Dat,ifelse(Eff.Reporter=="bad","bad",Reporter))  
    
    #remove records with NA effort
    Dat=subset(Dat,!(is.na(Km.Gillnet.Days.c) | is.na(Km.Gillnet.Hours.c)))
    
    #remove school shark or dogfish shots if not the target species
    idd=which(sp%in%c(17008,20000))
    if(length(idd)==0) Dat=subset(Dat,!(Sch.or.DogS=="Yes"))
  }else
  {
    Dat=NULL
  }
  return(Dat)
}

smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

fn.see.all.yrs.ves.blks=function(a,SP,NM,what,Ves.sel.BC,Ves.sel.sens,BLK.sel.BC,BLK.sel.sens,Min.ktch)
{
  All.ves=unique(as.character(a$VESSEL))
  All.blk=unique(as.character(a$BLOCKX))
  a=subset(a,Reporter=="good")
  dddd=subset(a,SPECIES%in%SP)
  if(nrow(dddd)>10)
  {
    CATCH.sp=dddd %>% group_by(YEAR.c,VESSEL)%>%
                      summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
                      spread(VESSEL, LIVEWT.c) %>%
                      arrange(YEAR.c) %>%
                      data.frame()
    Vess=names(CATCH.sp)[2:ncol(CATCH.sp)]
    Vess=chartr(".", " ", Vess)
    Yrs=CATCH.sp$YEAR.c
    Z=as.matrix(CATCH.sp[,-1])
    
    #Step 1. Select vessels with > X years of records of a minimum catch
    ZZ=Z
    ZZ[ZZ<Min.ktch]=NA
    ZZ[ZZ>=Min.ktch]=1
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    
    pdf(paste("C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Vessel_pos_records_by_yr/",paste(NM,what,sep=""),".pdf",sep="")) 
    
        #Ves.sel.BC
    par(mar=c(3,3.5,.8,.8))
    WHICh=which(Yrs.with.ktch>Ves.sel.BC)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.BC=Vess[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.BC=Ves.BC[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Ves.BC,las=1,cex.axis=.9)
    legend("top",paste("vessels with >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
    Drop.ves=All.ves[which(!All.ves%in%Ves.BC)]
    
        #Ves.sel.sens
    WHICh=which(Yrs.with.ktch>Ves.sel.sens)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.Sens=Vess[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.Sens=Ves.Sens[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Ves.Sens,las=1,cex.axis=.6)
    legend("top",paste("vessels with >=",Ves.sel.sens, "years of records"),bty='n')
    
    #plot CPUEs
    a$CPUE.km.gn.day=a$LIVEWT.c/a$Km.Gillnet.Days.c
    a$CPUE.km.gn.h=a$LIVEWT.c/a$Km.Gillnet.Hours.c
    a.mean.cpue.km.day_all=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,SPECIES%in%SP),mean)
    a.mean.cpue.km.h_all=aggregate(CPUE.km.gn.h~YEAR.c,subset(a, SPECIES%in%SP),mean)
    a.mean.cpue.km.day_Sens=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES%in%SP),mean)
    a.mean.cpue.km.h_Sens=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES%in%SP),mean)
    a.mean.cpue.km.day_BC=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES%in%SP),mean)
    a.mean.cpue.km.h_BC=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES%in%SP),mean)
    
    #kmgday
    par(mar=c(3,3,.5,4),mgp=c(2,.7,0))
    Yrs=a.mean.cpue.km.day_Sens$YEAR.c
    plot(a.mean.cpue.km.day_all$YEAR.c,a.mean.cpue.km.day_all$CPUE.km.gn.day,ylab="",xlab="")
    points(a.mean.cpue.km.day_Sens$YEAR.c,a.mean.cpue.km.day_Sens$CPUE.km.gn.day,pch=19,col=2)
    if(nrow(a.mean.cpue.km.day_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.day_BC$CPUE.km.gn.day,pch=19,col=3)
    
    #kmgday  
    par(new = T)
    plot(a.mean.cpue.km.h_all$YEAR.c,a.mean.cpue.km.h_all$CPUE.km.gn.h,ylab=NA, axes=F,xlab=NA,pch=0,cex=2)
    points(a.mean.cpue.km.h_Sens$YEAR.c,a.mean.cpue.km.h_Sens$CPUE.km.gn.h,pch=15,col=2,cex=2)
    if(nrow(a.mean.cpue.km.h_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.h_BC$CPUE.km.gn.h,pch=15,col=3,cex=2)
    axis(side = 4)
    mtext("Financial year",1,line=2)
    mtext("Nominal CPUE (Kg/km.gn.day)",2,line=2)
    mtext("Nominal CPUE (Kg/km.gn.hour)",4,line=2)
    legend("top",c("Kg/km.gn.day","Kg/km.gn.hour"),bty='n',pch=c(0,19))
    legend("bottomleft",c("all",paste(Ves.sel.sens,"y"),paste(Ves.sel.BC,"y")),bty='n',pch=c(0,19,19),col=c(1,2,3))
    
    
    # step 2. For selected vessels, plot number of blocks by year
      #Ves.BC
        #all blocks
    AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.BC) %>%
              group_by(YEAR.c,BLOCKX)%>%
              summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
              spread(BLOCKX, LIVEWT.c) %>%
              arrange(YEAR.c) %>%
              data.frame()
    BLOCs=substr(names(AA)[2:ncol(AA)],2,6)
    Yrs=AA$YEAR.c
    Z=as.matrix(AA[,-1])
    ZZ=Z
    ZZ[ZZ>0]=1
    ZZZ=ZZ
    ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
    ZZZ=ZZZ[,ID.sort]
    if(!is.matrix(ZZZ)) ZZZ=t(as.matrix(ZZZ))
    par(mar=c(3,3.5,.8,.8))
    image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
    legend("top",paste("all blocks for vessels >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
    
        #blocks with > bc records
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    WHICh=which(Yrs.with.ktch>BLK.sel.BC)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.BC=BLOCs[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.BC=Blks.BC[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Blks.BC,las=1,cex.axis=.9)
    legend("top",paste("blocks with >=",BLK.sel.BC, "years of records for vessels >=",Ves.sel.BC,"years of records and >",Min.ktch,"kg per year"),
           cex=0.75,bty='n')
    Drop.blks=All.blk[which(!All.blk%in%Blks.BC)]  
    
        #Ves.Sens
          #all blocks
    AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.Sens) %>%
      group_by(YEAR.c,BLOCKX)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
      spread(BLOCKX, LIVEWT.c) %>%
      arrange(YEAR.c) %>%
      data.frame()
    BLOCs=substr(names(AA)[2:ncol(AA)],2,6)
    Yrs=AA$YEAR.c
    Z=as.matrix(AA[,-1])
    ZZ=Z
    ZZ[ZZ>0]=1
    ZZZ=ZZ
    ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
    ZZZ=ZZZ[,ID.sort]
    if(!is.matrix(ZZZ)) ZZZ=t(as.matrix(ZZZ))
    par(mar=c(3,3.5,.8,.8))
    image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
    legend("top",paste("all blocks for vessels >=",Ves.sel.sens, "years of records"),bty='n')
    
        #blocks with > Sens records
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    WHICh=which(Yrs.with.ktch>BLK.sel.sens)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.Sens=BLOCs[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.Sens=Blks.Sens[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Blks.Sens,las=1,cex.axis=.9)
    legend("top",paste("blocks with >=",BLK.sel.sens, "years of records for vessels >=",Ves.sel.sens,"years of records of",SP),bty='n')
    dev.off()
    
    Drop.blks_10=Blks.BC_10=NULL
    if(what==".daily")
    {
      AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.BC) %>%
        group_by(YEAR.c,block10)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
        spread(block10, LIVEWT.c) %>%
        arrange(YEAR.c) %>%
        data.frame()
      AA=as.matrix(AA[,-1])
      AA[AA>0]=1
      Yrs.with.ktch=colSums(AA,na.rm=T)
      Blks.BC_10=substr(names(which(Yrs.with.ktch>BLK.sel.BC)),2,50)
      Drop.blks_10=unique(a$block10)[which(!unique(a$block10)%in%as.numeric(Blks.BC_10))]
    }
    return(list(Ves.BC=Ves.BC, Ves.Sens=Ves.Sens, Blks.BC=Blks.BC,Blks.BC_10=Blks.BC_10, Blks.Sens=Blks.Sens,
                Drop.ves=Drop.ves, Drop.blks=Drop.blks,Drop.blks_10=Drop.blks_10))
    
  }
}

#function for converting continuous var to factor
cfac=function(x,breaks=NULL)  
{
  if(is.null(breaks)) breaks=unique(quantile(x,probs = seq(0, 1, 0.1)))
  x=cut(x,breaks,include.lowest=T,right=F)
  levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
           c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
  return(x)
}

clog=function(x) log(x+0.05)   #function for applying log

#function for correlations
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)   
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#predictors effect
fn.pred.effect <- function(DATA,PREDS) 
{
  colnames(DATA)=tolower(colnames(DATA))
  PREDS=PREDS[which(PREDS%in%colnames(DATA))]
  DATA=DATA %>% mutate_each_(funs(factor(.)),PREDS[which(PREDS%in%Categorical)])%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  
  par(mfcol=c(3,2),mai=c(.6,.65,.3,.1),oma=c(.2,.2,.1,.1))
  hist(log(DATA$cpue),main="hist log(cpue)",xlab="log(cpue)",ylab="Count")
  
  boxcox(cpue+0.00001 ~ log(year.c), data = DATA,lambda = seq(0, 1, length = 10))
  legend("topright","Box Cox (should be small)",bty='n')
  
  # Cook distance to see outliers or overdisperse data (if distance >1)
  M1.1=glm(log(cpue+0.00001)~finyear,family=gaussian,data=DATA)
  plot(M1.1,which=4)
  legend("topright","outliers or overdispersion if distance >1",bty='n',cex=.85, text.col=2)
  
  plot(table(DATA$catch.target),type='h',xlab="Catch",ylab="Count",main="Catch zero inflation and right tail")
  
  #Outliers response var
  boxplot(DATA$cpue~DATA$finyear,main="Outliers in response var?",ylab="cpue (kg/km.gn.day)")
  
  #boxplot of response var and predictors
  smart.par(length(PREDS),c(2,2,2,.1),c(.1,.3,.1,.1),c(1.1,.35,0))
  for(d in 1:length(PREDS))
  {
    a=DATA[,match(c("cpue",PREDS[d]),names(DATA))]
    if(!is.factor(a[,2])) a[,2]=cfac(a[,2])
    plot(clog(a[,1])~a[,2],main=PREDS[d],ylab="",xlab="")
  }
  mtext("log(cpue)",side=2,line=-1,las=3,outer=T)
  
  #Covariate correlations.
  Covars=DATA%>%mutate(month=as.numeric(as.character(month)))%>%
    select(month,PREDS[which(!PREDS%in%Categorical)])
  pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor)
}

#functions for reshaping data
  #monthly data
Effort.data.fun=function(DATA,target,ktch)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  if(nrow(DATA)>0)
  {
    #remove nonsense lat
    DATA=subset(DATA,LAT>=(-36))
    
    #calculate effort (only need max effort)
    Match.these.eff=match(These.efforts,names(DATA))
    Effort.data1=DATA[,Match.these.eff]
    Effort.data=Effort.data1%>%
      group_by(zone,FINYEAR,Same.return,MONTH,BLOCKX,SHOTS.c)%>%
      summarise(Km.Gillnet.Days.c = max(Km.Gillnet.Days.c),
                Km.Gillnet.Hours.c = max(Km.Gillnet.Hours.c))%>%
      data.frame()
    
    #target species catch 
    #note: catch targeted at other species: pointless as these are multiple trips combined in one month
    ID=match(c(ktch),colnames(DATA))
    DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
    DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
    
    #reshape catch data
    TABLE=DATA%>%group_by(MONTH,FINYEAR,BLOCKX,VESSEL,Same.return,LAT,LONG,YEAR.c)%>%
      summarise(Catch.Target = sum(Catch.Target,na.rm=T),
                Catch.Total = sum(Catch.Total,na.rm=T))
    Enviro=DATA%>%group_by(MONTH,FINYEAR,BLOCKX)%>%
      summarise(Temperature=mean(Temperature),
                Temp.res=mean(Temp.res),
                Freo=mean(Freo))
    TABLE=TABLE%>%left_join(Enviro,by=c("FINYEAR","MONTH","BLOCKX"))    %>%
      arrange(FINYEAR,MONTH,BLOCKX) %>%
      data.frame()
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    
    #merge catch and effort
    dat=TABLE%>%left_join(Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"))
    
    #Add mesh size                                
    # d=subset(DATA,select=c(Same.return,mesh))
    # d=d[!duplicated(d$Same.return),]
    # dat=dat%>%left_join(d,by=c("Same.return"))
    
  }else
  {
    dat=DATA
    prop.with.catch=0
  }
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}

  #daily data
Effort.data.fun.daily=function(DATA,target,ktch,Aggregtn)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  if(nrow(DATA)>0)
  {
    #remove nonsense lat
    DATA=subset(DATA,LAT>=(-36))
    
    #calculate effort
    Match.these.eff=match(These.efforts.daily,names(DATA))
    Effort.data1=DATA[,Match.these.eff]
    
    #aggregate at shot level
    if(Use.Date=="NO")
    {
      #max effort by Sno, DsNo & TSNo
      Effort.data=Effort.data1%>%      
        group_by(zone,FINYEAR,Same.return.SNo,MONTH,BLOCKX,block10,shots.c)%>%
        summarise(Km.Gillnet.Days.c = max(Km.Gillnet.Days.c),
                  Km.Gillnet.Hours.c = max(Km.Gillnet.Hours.c))%>%
        data.frame()
      
      #aggregate at TSNo if required
      if(Aggregtn=="TSNo")
      {
        Effort.data$TSNo=word(Effort.data$Same.return.SNo,3)
        Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+
                                FINYEAR+TSNo+MONTH+BLOCKX,Effort.data,sum)
      }
    }
    
    #target species catch 
    ID=match(c(ktch),colnames(DATA))
    DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
    
    #catch targeted at other species
    DATA$Catch.Gummy=with(DATA,ifelse(SPECIES==17001,DATA[,ID],0))
    DATA$Catch.Whiskery=with(DATA,ifelse(SPECIES==17003,DATA[,ID],0))
    DATA$Catch.Dusky=with(DATA,ifelse(SPECIES%in%c(18003,18001),DATA[,ID],0))
    DATA$Catch.Sandbar=with(DATA,ifelse(SPECIES==18007,DATA[,ID],0))
    DATA$Catch.Groper=with(DATA,ifelse(SPECIES%in%c(384002),DATA[,ID],0))
    DATA$Catch.Snapper=with(DATA,ifelse(SPECIES%in%c(353001),DATA[,ID],0))
    DATA$Catch.Blue_mor=with(DATA,ifelse(SPECIES%in%c(377004),DATA[,ID],0))
    DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
    
    #reshape catch data
    if(Use.Date=="NO")
    {
      if(Aggregtn=="SNo") 
      {
        TABLE=DATA%>%group_by(MONTH,FINYEAR,BLOCKX,block10,VESSEL,Same.return.SNo,date,LAT,LONG,YEAR.c,Lunar)%>%
          summarise(Catch.Target = sum(Catch.Target,na.rm=T),
                    Catch.Gummy=sum(Catch.Gummy,na.rm=T),
                    Catch.Whiskery=sum(Catch.Whiskery,na.rm=T),
                    Catch.Dusky=sum(Catch.Dusky,na.rm=T),
                    Catch.Sandbar=sum(Catch.Sandbar,na.rm=T),
                    Catch.Groper=sum(Catch.Groper,na.rm=T),
                    Catch.Snapper=sum(Catch.Snapper,na.rm=T),
                    Catch.Blue_mor=sum(Catch.Blue_mor,na.rm=T),
                    Catch.Total=sum(Catch.Total,na.rm=T))
        Enviro=DATA%>%group_by(MONTH,FINYEAR,BLOCKX)%>%
          summarise(Temperature=mean(Temperature),
                    Temp.res=mean(Temp.res),
                    Freo=mean(Freo))
        TABLE=TABLE%>%left_join(Enviro,by=c("FINYEAR","MONTH","BLOCKX"))    %>%
          arrange(Same.return.SNo,FINYEAR,MONTH,BLOCKX) %>%
          data.frame()
      }
      #aggregating by trip
      if(Aggregtn=="TSNo")   
      {
        TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Dusky,Catch.Sandbar,
                              Catch.Groper,Catch.Snapper,Catch.Blue_mor,Catch.Dhufish,
                              Catch.Other.shrk,Catch.Other.scalefish,
                              Catch.non_indicators,Catch.Total)~MONTH+FINYEAR+BLOCKX+VESSEL+
                          TSNo+YEAR.c,data=DATA,sum,na.rm=T)
        xx=subset(DATA,select=c(BLOCKX,LAT,LONG))
        xx=xx[!duplicated(xx$BLOCKX),]
        xx$LAT=round(xx$LAT)
        xx$LONG=round(xx$LONG)
        TABLE=merge(TABLE,xx,by="BLOCKX",all.x=T)
        TABLE=TABLE[order(TABLE$TSNo,TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
      }
    }
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    
    #merge catch and effort
    if(Aggregtn=="SNo") dat=TABLE%>%left_join(Effort.data,by=c("Same.return.SNo","FINYEAR","MONTH","BLOCKX","block10"))
    if(Aggregtn=="TSNo") dat=TABLE%>%left_join(Effort.data,by=c("TSNo","FINYEAR","MONTH","BLOCKX","block10"))
    
    
    #Add mesh size, shots, depth and nlines for each session
    if(Aggregtn=="SNo")
    {
      d=subset(DATA,select=c(Same.return.SNo,VESSEL,mesh,nlines.c,Mean.depth))
      d=d[!duplicated(paste(d$Same.return.SNo,d$VESSEL)),]
      dat=dat%>%left_join(d,by=c("Same.return.SNo","VESSEL"))
    }
    
  }else
  {
    dat=DATA
    prop.with.catch=0
  }
  
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}

fn.box.plt.year=function(d)
{
  d$cpue=d$Catch.Target/d$Km.Gillnet.Days.c
  par(mfcol=c(3,1),mar=c(4,4,1,1),mgp=c(2,.7,0))
  boxplot(Catch.Target~FINYEAR,d,ylab="Catch (kg)")
  boxplot(Km.Gillnet.Days.c~FINYEAR,d,ylab="Km gn d")
  boxplot(cpue~FINYEAR,d,ylab="cpue")
  
  FINYRS=sort(unique(d$FINYEAR))
  for(xx in 1:length(FINYRS))    boxplot(cpue~BLOCKX,subset(d,FINYEAR==FINYRS[xx]),ylab="cpue",main=FINYRS[xx])
  
  Yrs=as.numeric(substr(sort(unique(d$FINYEAR)),1,4))
  
  Ag.blk=aggregate(cpue~FINYEAR+BLOCKX,d,mean)
  Ag.blk= reshape(Ag.blk, v.names = "cpue", idvar="BLOCKX",timevar ="FINYEAR", direction = "wide")
  blks=Ag.blk$BLOCKX
  Ag.blk=Ag.blk[,-match('BLOCKX',names(Ag.blk))]
  CL=rainbow(ncol(Ag.blk))
  plot(Yrs,Ag.blk[1,],col=CL[1],pch=19,ylim=c(0,max(Ag.blk,na.rm=T)))
  for(qq in 2:ncol(Ag.blk)) points(Yrs,Ag.blk[qq,],pch=19,col=CL[qq])
  legend("topright",paste(blks),bty='n',pch=19,col=CL)
  
  Ag.vsl=aggregate(cpue~FINYEAR+VESSEL,d,mean)
  Ag.vsl= reshape(Ag.vsl, v.names = "cpue", idvar="VESSEL",timevar ="FINYEAR", direction = "wide")
  vsls=Ag.vsl$VESSEL
  Ag.vsl=Ag.vsl[,-match('VESSEL',names(Ag.vsl))]
  CL=rainbow(ncol(Ag.vsl))
  plot(Yrs,Ag.vsl[1,],col=CL[1],pch=19,ylim=c(0,max(Ag.vsl,na.rm=T)))
  for(qq in 2:ncol(Ag.vsl)) points(Yrs,Ag.vsl[qq,],pch=19,col=CL[qq])
  legend("topright",paste(vsls),bty='n',pch=19,col=CL)
}


fn.check.eff.area=function(ALL,WHAT,Effrt,QL_prop_ktch,spname,TYPE)
{
  id=match(WHAT,names(ALL$whiskery))
  ALL.shots=unique(ALL$whiskery[,id])
  gum.shots=unique(ALL$gummy[,id])
  dus.shots=unique(ALL$dusky[,id])
  san.shots=unique(ALL$sandbar[,id])
  gum.shots=gum.shots[which(!gum.shots%in%ALL.shots)]
  dus.shots=dus.shots[which(!dus.shots%in%ALL.shots)]
  san.shots=san.shots[which(!san.shots%in%ALL.shots)]
  dat=ALL$whiskery
  if(length(gum.shots)>0) dat=rbind(dat,ALL$gummy[ALL$gummy[,id]%in%gum.shots,])
  if(length(dus.shots)>0) dat=rbind(dat,ALL$dusky[ALL$dusky[,id]%in%dus.shots,])
  if(length(san.shots)>0) dat=rbind(dat,ALL$sandbar[ALL$sandbar[,id]%in%san.shots,])

  Ktch.targt=paste("catch",spname,sep=".")
  names(dat) =  casefold(names(dat))
  
  ## some variable names
  dat$catch = dat[,match(Ktch.targt,names(dat))]
  dat$effort = dat[,match(Effrt,names(dat))]
  dat$year = dat$year.c
  dat$fymonth = factor(dat$month, levels=c(7:12, 1:6))
  dat$season = as.numeric(substring(dat$finyear, 1, 4))
  dat$smonth = factor(dat$month, levels=c(7:12, 1:6))
  
  
  
  ## calculate year-specific QL and add column to data set
  dat$prop = dat$catch / dat$catch.total
  
  qldat = CalcQL(dat, prop.catch=QL_prop_ktch)
  dat = merge(dat, qldat, all=TRUE)
  dat$target = ifelse(dat$prop > dat$ql, 1, 0)
  dat=dat[order(dat$target),]
  
  dat$col=ifelse(dat$target==0,"grey80","grey20")
  uni.rs=sort(unique(dat$finyear))
  
  fn.fig(paste("Outputs/Effective_Area/",spname,"_Qualif_level_",TYPE,sep=""),2400,2400)
  smart.par(n.plots=length(uni.rs),MAR=c(1.5,1.5,1,.5),OMA=c(3,3,.5,.5),MGP=c(1,.5,0))
  for(u in 1:length(uni.rs))
  {
    with(subset(dat,finyear==uni.rs[u]),plot(long,lat,col=col,pch=19,ylab="",xlab="",
                                             main=uni.rs[u],cex.main=1.25,ylim=c(-36,-26),xlim=c(113,129)))
  }
  legend("topright",c("targeted","not targeted"),bty='n',pch=19,col=c("grey20","grey80"),cex=1.5)
  mtext("Longitude",1,line=.7,outer=T,cex=1.5)
  mtext("Latitude",2,line=.7,outer=T,las=3,cex=1.5)
  dev.off()
  
  
  #cpue
  dat$cpue=dat$catch/dat$effort
  scaler=5
  dat$CX=mapply(fn.scale,dat$cpue,scaler)
  dat=dat[order(dat$cpue),]
  fn.fig(paste("Outputs/Effective_Area/",spname,"_Pos_cpue_",TYPE,sep=""),2400,2400)
  smart.par(n.plots=length(uni.rs),MAR=c(1.5,1.5,1,.5),OMA=c(3,3,.5,.5),MGP=c(1,.5,0))
  for(u in 1:length(uni.rs))
  {
    with(subset(dat,finyear==uni.rs[u]),plot(long,lat,cex=CX,pch=19,ylab="",xlab="",
                                             col=rgb(.1,.1,.1,.2),main=uni.rs[u],cex.main=1.25,ylim=c(-36,-26),xlim=c(113,129)))
  }
  Qnt=round(quantile(dat$cpue[dat$cpue>0],probs=c(.25,.75,1)),1)
  legend("topright",paste(Qnt),bty='n',pch=19,col=rgb(.1,.1,.1,.2),
         pt.cex=mapply(fn.scale,Qnt,scaler),cex=1.5,title="kg/km gn d")
  mtext("Longitude",1,line=.7,outer=T,cex=1.5)
  mtext("Latitude",2,line=.7,outer=T,las=3,cex=1.5)
  dev.off()
}

export.foly=function(DATA,DATA1)     #function for creating foly indices 
{
  DATA=subset(DATA,Catch.Target>0)
  DATA1=subset(DATA1,Catch.Target>0)
  
  Ktch=aggregate(Catch.Target~FINYEAR,DATA,sum)
  Ktch_d=aggregate(Catch.Target~FINYEAR,DATA1,sum)
  
  Efrt=aggregate(Km.Gillnet.Days.c~FINYEAR,DATA,sum)
  Efrt_d=aggregate(Km.Gillnet.Days.c~FINYEAR,DATA1,sum)
  
  foly=merge(Ktch,Efrt,by="FINYEAR")
  foly_d=merge(Ktch_d,Efrt_d,by="FINYEAR")
  foly$cpue=foly$Catch.Target/foly$Km.Gillnet.Days.c
  foly_d$cpue=foly_d$Catch.Target/foly_d$Km.Gillnet.Days.c
  
  return(rbind(subset(foly,select=c(FINYEAR,cpue)),subset(foly_d,select=c(FINYEAR,cpue))))
}

#calculate 4 different nominal cpues and choose data based on qualification level (90% of years catch)
fn.ainslie=function(dat,Ktch.targt,Effrt,explr,QL_prop_ktch,Prop.Malcolm,cpue.units,spname,BLks,VesL,Type)
{
  names(dat) =  casefold(names(dat))
  
  ## some variable names
  dat$catch = dat[,match(Ktch.targt,names(dat))]
  dat$year = dat$year.c
  dat$fymonth = factor(dat$month, levels=c(7:12, 1:6))
  dat$season = as.numeric(substring(dat$finyear, 1, 4))
  dat$smonth = factor(dat$month, levels=c(7:12, 1:6))
  
  CPUE.All=vector('list',length(Effrt))
  names(CPUE.All)=Effrt
  CPUE.blk_vess=CPUE.non_zero=CPUE.QL_target=CPUE.Malcolm=CPUE.All
  
  for(ef in 1:length(Effrt))
  {
    pdf(paste(Hnd.ains,spname,Type,paste(Effrt[ef]),".pdf",sep=""))
    
    dat$effort = dat[,match(Effrt[ef],names(dat))]
    
    if(explr=="YES")
    {
      table(dat$year, dat$month)
      table(dat$season, dat$smonth)
      
      all.years = sort(unique(dat$year))
      par(mfrow=c(2,2))
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.gummy, list(year, month), sum)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.whiskery, list(year, month), sum)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.dusky, list(year, month), sum)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.sandbar, list(year, month), sum)), xlab="", ylab="", las=1, main="Sandbar")
      
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.gummy/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.whiskery/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.dusky/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.sandbar/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Sandbar")
      
      
      all.seasons = sort(unique(dat$season))
      par(mfrow=c(2,2))
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.gummy, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.whiskery, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.dusky, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.sandbar, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Sandbar")
      
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.gummy/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.whiskery/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.dusky/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.sandbar/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Sandbar")
      
      
      par(mfrow=c(2,1))
      boxplot(catch.target/km.gillnet.days.c ~ year, dat)
      boxplot(catch.target/km.gillnet.days.c ~ season, dat)
      
      boxplot(catch.target ~ month, dat, main="")
      boxplot(catch.target/km.gillnet.days.c ~ month, dat)
      
    }
    
    
    ## calculate year-specific QL and add column to data set
    
    dat$prop = dat$catch.target / dat$catch.total
    
    if(ef==1)
    {
      par(mfrow=c(1,1))
      boxplot(prop ~ season, dat)
      mtext("Proportion of target species catch out of total catch",3,-1,col="red")
    } 
    qldat = CalcQL(dat, prop.catch=QL_prop_ktch)
    dat = merge(dat, qldat, all=TRUE)
    dat$target = ifelse(dat$prop > dat$ql, 1, 0)
    if(ef==1)
    {
      par(mfrow=c(2,1))
      boxplot(catch.target/km.gillnet.days.c ~ season, dat)
      mtext("Qualification levels_all",3,-1,col="red")
      boxplot(catch.target/km.gillnet.days.c ~ season, subset(dat, target==1))
      mtext("Qualification levels_target only",3,-1,col="red")
    }
    
    
    ## malcolm's targeting - proportion below which no variation exists in cpue
    
    dat$cpue.target = dat$catch / dat$effort
    all.seasons = unique(dat$season)
    if(ef==1)
    {
      smart.par(n.plots=length(all.seasons+1),MAR=c(2,2,.1,.1),OMA=c(2,2,.5,.5),MGP=c(1,.5,0))
      with(dat, plot(prop, cpue.target,  pch=16, col=rgb(1,0,0,0.1), ylab='',xlab=''))
      legend('top',"All Seasons",bty='n')
      ## hmmm this is hard to see - probably very low <0.1
      for (i in all.seasons)
      {
        with(subset(dat, season==i), plot(prop, cpue.target, pch=16, ylab='',xlab='', col=rgb(0,0,1,0.1)))
        legend('top',paste(i),bty='n')
      }
      ## still say <0.1
      mtext(paste("Cpue (",cpue.units[ef],")",sep=""),2,outer=T,las=3)
      mtext("Proportion",1,outer=T)
    }

    ## plot raw mean cpues and CIs using 4 different data sets
    
    #compare different subsets of the data
    par(mfrow=c(3,2),mar=c(3,3,1,1),oma=c(1,1,1,1),mgp=c(2,.5,0))
    
    CPUE.All[[ef]] = CalcMeanCPUE(cpuedata = dat, catch.column="catch", effort.column="effort",
                                  plot.title = paste(spname, "_All rec"), cpue.units = cpue.units[ef], 
                                  draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
    
    CPUE.blk_vess[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, blockx%in%BLks & vessel%in%VesL), catch.column="catch", 
                                       effort.column="effort",plot.title = paste(spname, "_indicative_blk_ves"),
                                       cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
    
    CPUE.non_zero[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, catch>0), catch.column="catch", 
                                       effort.column="effort",plot.title = paste(spname, "_Nonzero"),
                                       cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
    
    CPUE.QL_target[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, target==1), catch.column="catch",
                                        effort.column="effort",plot.title = paste(spname, "_QL Target"),
                                        cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
    
    if(nrow(subset(dat, prop>Prop.Malcolm))>100)CPUE.Malcolm[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, prop>Prop.Malcolm), catch.column="catch",
                                      effort.column="effort",plot.title = paste(spname, "_Malcolm_Prop>",Prop.Malcolm),
                                      cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=TRUE,PaR="NO",showLNMean="YES")
    dev.off()
  }
  
  #show proportion of records selected by year
  #note: show that prob of catching doesn't change thru time
  TAB=subset(dat,select=c(finyear,target)) %>%
    group_by(finyear,target) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n))%>%
    as.data.frame
  TAB$yr=substr(TAB$finyear,1,4)
  
  return(list(CPUE.All=CPUE.All,CPUE.blk_vess=CPUE.blk_vess,CPUE.non_zero=CPUE.non_zero,
              CPUE.QL_target=CPUE.QL_target,CPUE.Malcolm=CPUE.Malcolm,
              QL_dat=subset(dat, target==1),Prop.selected=TAB))
}

fn.sel.yrs.used=function(DD,ThrShld.n.vess)
{
  a=with(DD,table(FINYEAR,VESSEL))
  a[a>0]=1
  a=rowSums(a)
  a[a<ThrShld.n.vess]=NA
  return(names(a[which(!is.na(a))]))
}
fn.sel.yrs.used.glm=function(DD)
{
  a=with(DD,table(finyear,vessel))
  a[a>0]=1
  a=rowSums(a)
  a[a<Threshold.n.vessls.per.yr]=NA
  return(names(a[which(!is.na(a))]))
}

fn.check.balanced=function(d,SP,what,MN.YR,pLot)
{
  fn.plt=function(dd)
  {
    a=dd
    a[a>0]=1
    Orderd=rev(sort(rowSums(a)))
    dd=as.matrix(dd[match(names(Orderd),row.names(dd)),])
    Nx=c(1,ncol(dd))
    Ny=c(1,nrow(dd)+2)
    Mx=max(dd)
    par(mar=c(2,3,.8,.8),mgp=c(2,.5,0))
    plot(1,1,xlim=Nx,ylim=Ny,col="transparent",ann=F,xaxt='n',yaxt='n')
    Selected=names(Orderd[Orderd>=MN.YR])
    show.pol=match(Selected,row.names(dd))
    Nx.p=c(Nx[1]-1,Nx[2]+1)
    if(length(Selected)>0)polygon(x=c(Nx.p,rev(Nx.p)),y=c(rep(show.pol[1]-1,2),rep(show.pol[length(show.pol)],2)),col=rgb(.1,.1,.1,.25),border='transparent')
    for(i in 1:nrow(dd)) points(Nx[1]:Nx[2],rep(i,Nx[2]),pch=21,cex=fn.scale(dd[i,],2.5),bg=rgb(.1,.1,.1,.4))
    axis(1,1:ncol(dd),colnames(dd))
    axis(2,1:nrow(dd),rownames(dd),las=1,cex.axis=.5)
    mtext("number of records per year",3,1,cex=1.25)
    Lab=round(quantile(dd,probs=c(.95,.995,1)))
    legend('topright',paste(Lab),pch=21,pt.bg="grey70",
           pt.cex=fn.scale(Lab,2.5),horiz=T,title="# of records")
    
    return(Selected)
  }
  
  if(pLot)pdf(paste(HndL,paste(SP,"_",what,sep=""),".pdf",sep="")) 
  
  #First, select vessels 
  Ves.Yr=with(d,table(vessel,finyear))
  this.ves=fn.plt(Ves.Yr)
  mtext("Vessel (all)",2,line=1.65,cex=1.25)
  
  #Second, select blocks for selected vessel
  BLK.Yr=with(subset(d,vessel%in%this.ves),table(blockx,finyear))
  this.blks=fn.plt(BLK.Yr)
  mtext("Block (for selected Vessels)",2,line=1.65,cex=1.25)
  
  
  #Third, keep only selected blocks and vessels
  d=subset(d,blockx%in%this.blks)
  d=subset(d,vessel%in%this.ves)
  
  BLK.Yr=with(d,table(blockx,finyear))
  this.blks=fn.plt(BLK.Yr)
  mtext("Block (for selected block and vessel)",2,line=1.65,cex=1.25)
  
  Ves.Yr=with(d,table(vessel,finyear))
  this.ves=fn.plt(Ves.Yr)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  
  MN.Yr=with(d,table(month,finyear))
  kk=fn.plt(MN.Yr)
  mtext("Month (for selected block and vessel)",2,line=1.65,cex=1.25)

  
  MN.Ves=with(d,table(vessel,month))
  kk=fn.plt(MN.Ves)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Month",1,line=1,cex=1.1)
  
  MN.blk=with(d,table(blockx,month))
  kk=fn.plt(MN.blk)
  mtext("Block (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Month",1,line=1,cex=1.1)
  
  ves.blk=with(d,table(vessel,blockx))
  kk=fn.plt(ves.blk)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Block",1,line=1,cex=1.1)
  
  if(pLot)dev.off()
  
  
  return(list(this.blks=this.blks,this.ves=this.ves))
}

Get.Mns=function(d,grp,Vars,LGND,add.arrow)
{
  d=d[,match(c(grp,Vars),names(d))]
  d$cpue.d=d$Catch.Target/d$Km.Gillnet.Days.c
  d$cpue.h=d$Catch.Target/d$Km.Gillnet.Hours.c
  d$cpue.h_shot=d$Catch.Target/(d$Km.Gillnet.Hours_shot.c)
  for(v in 1:length(Vars))
  {
    if(Vars[v]=="HOURS.c")
    {
      ddd=subset(d,SHOTS.c==1)
      B= ddd[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
        mutate(lowCL = mean - 1.96*sd/sqrt(n),
               uppCL = mean + 1.96*sd/sqrt(n)) %>%
        as.data.frame
      B$yr=as.numeric(substr(B$FINYEAR,1,4))

      ddd=subset(d,SHOTS.c==2)
      B2= ddd[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
        mutate(lowCL = mean - 1.96*sd/sqrt(n),
               uppCL = mean + 1.96*sd/sqrt(n)) %>%
        as.data.frame
      B2$yr=as.numeric(substr(B2$FINYEAR,1,4)) 
      
      plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col=1,ylim=c(min(c(B$lowCL,B2$lowCL)),max(c(B$uppCL,B2$uppCL))))
      arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col=1)
      #if(add.arrow[v])
      #{
     #   Is=(nrow(B)-5):nrow(B)
    #    arrows(B$yr[1],B$mean[1],mean(B$yr[Is]),mean(B$mean[Is]),col=1,lwd=2)
      #  legend("bottomright",paste(round(mean(B$mean[Is])/B$mean[1],1),"fold",sep="-"),bty='n',cex=1)
      #  with(B[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.1,.8,.15)))
   #   }
      
      points(B2$yr,B2$mean,pch=19,col='grey65')
      B2$lowCL=ifelse(B2$lowCL==B2$mean,B2$lowCL*.999,B2$lowCL)
      B2$uppCL=ifelse(B2$uppCL==B2$mean,B2$uppCL*1.0001,B2$uppCL)
      arrows(x0=B2$yr, y0=B2$lowCL,x1=B2$yr, y1=B2$uppCL,code=3, angle=90, length=0.05, col='grey65')
    #  if(add.arrow[v])
    #  {
    #    Is=(nrow(B2)-5):nrow(B2)
     #   arrows(B2$yr[1],B2$mean[1],mean(B2$yr[Is]),mean(B2$mean[Is]),col='forestgreen',lwd=2)
     #   legend("topright",paste(round(mean(B2$mean[Is])/B2$mean[1],1),"fold",sep="-"),bty='n',cex=1)
    #    with(B2[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.8,.1,.15)))
    #  }
      legend("topleft",c("1 shot","2 shots"),text.col=c("black","grey65"),bty='n',cex=1.5)
      
    }
    if(Vars[v]=="SHOTS.c")
    {
      ddd=subset(d,SHOTS.c%in%c(1,2))
      
      B=ddd[,match(c(grp,Vars[v]),names(d))] %>%
        group_by(FINYEAR,SHOTS.c) %>%
        summarise (n = n()) %>%
        mutate(freq = n / sum(n))%>%
        as.data.frame
      B$yr=as.numeric(substr(B$FINYEAR,1,4))
      B1=reshape(subset(B,select=c(SHOTS.c,freq,yr)),
              v.names = "freq", idvar = "SHOTS.c",
              timevar = "yr", direction = "wide")
      barplot(as.matrix(B1[,-1]),names.arg=unique(B$yr),legend.text=c("1 shot","2 shots"))
      box()
      
    }
    if(!Vars[v]%in%c("HOURS.c","SHOTS.c")){
      B= d[,match(c(grp,Vars[v]),names(d))] %>%
        na.omit()%>%
        group_by(FINYEAR) %>%
        summarise_all(funs(mean=mean,sd=sd,n=length)) %>%
        mutate(lowCL = mean - 1.96*sd/sqrt(n),
               uppCL = mean + 1.96*sd/sqrt(n)) %>%
        as.data.frame
      B$yr=as.numeric(substr(B$FINYEAR,1,4))
      plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col=1,ylim=c(min(B$lowCL),max(B$uppCL)))
      arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col=1)
  #    if(add.arrow[v])
  #    {
   #     Is=(nrow(B)-5):nrow(B)
    #    arrows(B$yr[1],B$mean[1],mean(B$yr[Is]),mean(B$mean[Is]),col='black',lwd=2)
   #     legend("bottomright",paste(round(mean(B$mean[Is])/B$mean[1],1),"fold",sep="-"),bty='n',cex=1.5)
    #    with(B[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.1,.8,.15)))
     # }
    }

  }
  
  #cpues
  B= d[,match(c(grp,"cpue.d"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr,B$mean,xlab="",ylab="CPUE (kg gn d)",pch=19,col='black',ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='black')
  axis(2,col="black",col.axis = "black")
  
  par(new=T)
  B= d[,match(c(grp,"cpue.h"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr,B$mean,xlab="",ylab="",pch=19,col='grey50',axes=F,ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='grey50')
  
  
  B= d[,match(c(grp,"cpue.h_shot"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  arrows(x0=B$yr+.5, y0=B$lowCL,x1=B$yr+.5, y1=B$uppCL,code=3, angle=90, length=0.05, col='black')
  points(B$yr+.5,B$mean,xlab="",ylab="",pch=21,col='black',bg="white")
  axis(side = 4,col="grey50",col.axis = "grey50")
  mtext("CPUE (kg gn h)",4,line=1.5,col="grey50",cex=1)
  legend("top",c("kg/km gillnet days","kg/km gillnet hours","kg/km gillnet hours_shot"),
         bty='n',col=c("black","grey50","black"),cex=1.35,pt.bg=c("black","grey50","white"),pch=21)
}

#show blocks and vessels kept
fn.show.blk=function(dat,CEX,SRt) 
{
  dat=sort(dat)
  LAT.kept=sapply(dat, function(x) -as.numeric(substr(x, 1, 2)))
  LONG.kept=sapply(dat, function(x) 100+as.numeric(substr(x, 3, 4)))
  
  Y=-36:-26; X=seq(113,129,length.out=length(Y))
  plot(X,Y,ylab='',xlab="",col="transparent",cex.lab=1.5,cex.axis=1.25)
  for(e in 1:length(LAT.kept))
  {
    dd.y=c(LAT.kept[e]-1,LAT.kept[e]-1,LAT.kept[e],LAT.kept[e])
    dd.x=c(LONG.kept[e],LONG.kept[e]+1,LONG.kept[e]+1,LONG.kept[e])
    polygon(dd.x,dd.y,col=rgb(0, 0, 1,0.25), border=rgb(0, 0, 1,0.5))
    text(LONG.kept[e]+0.5,LAT.kept[e]-0.5,dat[e],cex=CEX,col=1,srt=SRt,font=2)
  }
}

#nice table
Export.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                    body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                    HDR.names,HDR.span,HDR.2nd)
{
  mydoc = docx(Doc.nm)  #create r object
  mydoc = addSection( mydoc, landscape = T )   #landscape table
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                     header.cell.props = cellProperties(background.color=HdR.bg), 
                     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                        font.weight="bold",font.family =Fnt.hdr), 
                     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
  #Add header
  MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
  
  #Add second header
  MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
  
  
  # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
                                 inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                 outer.vertical = borderNone(),
                                 outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  mydoc = addFlexTable( mydoc, MyFTable)   
  mydoc = addSection( mydoc, landscape = F ) 
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}


check.cpue=function(DATA,NAME,cl)   #function for checking cpue outliers
{
  Ktc.q=quantile(DATA$catch.target,probs=seq(0,1,.1))
  Eff.q=quantile(DATA$km.gillnet.hours.c,probs=seq(0,1,.1))
  CPUE.q=quantile(DATA$cpue.target,probs=seq(0,1,.1))
  par(mfcol=c(3,1),mai=c(.2,.5,.3,.1),mgp=c(2.5,.5,0))
  boxplot(catch.target~finyear,DATA,ylab="KG",main=NAME,col="grey80")
  abline(h=Ktc.q[6],col=cl,lwd=2)
  text(1,Ktc.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=Ktc.q[10],col=cl,lwd=2)
  text(1,Ktc.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(Ktc.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
  
  boxplot(km.gillnet.hours.c~finyear,DATA,ylab="km gn hr",col="grey80")
  abline(h=Eff.q[6],col=cl,lwd=2)
  text(1,Eff.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=Eff.q[10],col=cl,lwd=2)
  text(1,Eff.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(Eff.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
  
  boxplot(cpue.target~finyear,DATA,ylab="KG / km gn hr",col="grey80")
  abline(h=CPUE.q[6],col=cl,lwd=2)
  text(1,CPUE.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=CPUE.q[10],col=cl,lwd=2)
  text(1,CPUE.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(CPUE.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
}           

fn.plt=function(a,y,TITL)
{
  plot(1:nrow(a),ylim=c(0,max(a,na.rm=T)),col="transparent",ann=F,axes=F)
  CL=rainbow(ncol(a))
  for(pp in 1:ncol(a)) lines(1:nrow(a),a[,pp],col=CL[pp],lwd=4)
  axis(1,1:nrow(a),rownames(a))
  nn=seq(0,max(a,na.rm=T),length.out=5)
  axis(2,nn,round(nn))
  mtext(y,2,3,las=3,cex=1.5)
  legend("topright",colnames(a),text.col=CL,bty='n',title=TITL)
}

fn.expl.cede=function(d,PREDS,kg,Do.ggplts)    #function for exploratory analysis
{
  PREDS=PREDS[which(PREDS%in%colnames(d))]
  
  Yrs=length(unique(d$year.c))
  div=1
  if(kg) div=1000  #in tonnes
  
  fn.plt(tapsum(d,"catch.target","finyear","month",div=div),"Catch","month")
  fn.plt(tapsum(d,"catch.target","finyear","zone",div=div),"Catch","zone")
  fn.plt(tapsum(d,"km.gillnet.hours.c","finyear","zone",div=1.0),"Effort","km.gillnet.hours")
  
  #explore turn over of vessel per year
  cbv <- tapsum(d,"catch.target","vessel","year.c",div=div) # often more vessels than years
  total <- rowSums(cbv,na.rm=TRUE)
  cbv1 <- cbv[order(total),] 
  to <- turnover(cbv1)    
  yearBubble(cbv1,ylabel="sqrt(catch-per-vessel)",diam=0.125,txt=c(2,3,4,5),hline=TRUE)
  
  plot.new()
  grid.table(to)
  
  #depth bin selection
  if(!is.na(match("mean.depth",PREDS)))
  {
    par(mfcol=c(2,2),mar=c(2,2,2,.1))
    barplot(table(trunc(d$mean.depth/2) * 2),main="2 m bin")
    barplot(table(trunc(d$mean.depth/5) * 5),main="5 m bin")
    barplot(table(trunc(d$mean.depth/10) * 10),main="10 m bin")
    barplot(table(trunc(d$mean.depth/25) * 25),main="25 m bin")
    mtext("Depth categories",3,-2,outer=T,col=2)
    
    cc <- histyear(d,Lbound=0,Rbound=max(d$mean.depth)*1.1,inc=10,pickvar="mean.depth",
                   years="year.c",varlabel="Depth (m)",plots=n2mfrow(Yrs),vline=120)
    
    d$DepCat=trunc(d$mean.depth/10) * 10
    
  }
  
  if(!is.na(match("nlines.c",PREDS)))
  {
    par(mfcol=c(1,1),mar=c(1,1,2,1))
    barplot(table(d$nlines.c),main="n lines")
  }
  if(!is.na(match("mesh",PREDS)))
  {
    barplot(table(d$mesh),main="mesh")
  }
  
  #effort distribution
  outH <- histyear(d,Lbound=0,Rbound=max(d$km.gillnet.hours.c)*1.1,inc=10,pickvar="km.gillnet.hours.c",
                   years="year.c",varlabel="km.gillnet.hours",plots=n2mfrow(Yrs),vline=NA)
  
  #catch vs effort
  par(mfrow=c(1,1),mai=c(0.45,0.45,0.05,0.05),cex=0.85, 
      mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)  
  plot(d$km.gillnet.hours.c,d$catch.target,type="p",pch=16,col=rgb(1,0,0,1/5),
       ylim=c(0,max(d$catch.target)),xlab="km.gillnet.hours",ylab="Catch")
  abline(h=0.0,col="grey")
  
  #Exploration of Spatial distribution of data
  leftlong <- 113;  rightlong <- 129
  uplat <- -26;  downlat <- -36
  plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
  dd=subset(d,catch.target>0,select=c(lat,long,catch.target))
  names(dd)=c("Lat","Long","catch.target")
  addpoints(dd,intitle="Location of Positive catches")
  
  plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
  plotpolys(dd,leftlong,rightlong,uplat,downlat,gridon=1,leg="left",
            intitle="1 degree squares",mincount=2,namecatch="catch.target",textout = F)
  
  
  #cpue distribution by year and by month
  d$cpue=d$catch.target/d$km.gillnet.hours.c
  d$LnCE=log(d$cpue)
  cc=histyear(d,Lbound=min(d$LnCE)*1.1,Rbound=max(d$LnCE)*1.1,inc=0.2,pickvar="LnCE",
              years="year.c",varlabel="log(CPUE)",plots=n2mfrow(Yrs))
  cc=histyear(d,Lbound=min(d$LnCE)*1.1,Rbound=max(d$LnCE)*1.1,inc=0.2,pickvar="LnCE",
              years="month",varlabel="log(CPUE)",plots=n2mfrow(12))
  
  #cpue boxplots
  Cat=PREDS[which(PREDS%in%Categorical)]
  dd <- makecategorical(Cat,d)
  
  smart.par(length(PREDS),c(2.5,2.5,.1,.1),c(1,1,1,1),c(1.5,.5,0))
  for(pp in 1:length(PREDS))
  {
    x=dd[,match(c("cpue",PREDS[pp],"finyear"),names(dd))]
    if(!(is.factor(x[,2])|is.integer(x[,2]))) x[,2]=cut(x[,2],breaks=quantile(x[,2]))
    boxplot(x$cpue~x[,2],ylab="cpue",xlab=PREDS[pp],notch=F,varwidth=T)
    #"varwidth=T": box widths proportional to the square roots of the sample sizes
  }  
  
  if(Do.ggplts)
  {
    for(pp in 1:length(PREDS))
    {
      x=dd[,match(c("cpue",PREDS[pp],"finyear"),names(dd))]
      if(!(is.factor(x[,2])|is.integer(x[,2])))
      {
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue,color = finyear))+geom_point()
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue))+geom_density_2d()+
          xlab(PREDS[pp]) +geom_point(alpha=0.2,color="brown",size=1.5)
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue))+geom_point()+geom_smooth(aes(color = finyear)) +facet_wrap( ~ finyear)
      }
    } 
    ggplot(dd, aes(x = vessel, y = year.c)) +  geom_jitter()
    ggplot(dd, aes(x = year.c, y = cpue)) + geom_jitter(alpha = 0.6) + facet_wrap( ~ vessel) + coord_flip()
    ggplot(dd, aes(x = cpue, fill = vessel)) +geom_histogram(bins = 25)
    ggplot(dd, aes(x = month, y = cpue, color = cluster_clara)) +   geom_boxplot()
  }
  
}

fn.sel.discrete.dist=function(d)
{
  P <- fitdist(d, "pois",discrete=T)
  nb <- fitdist(d, "nbinom",discrete=T)
  a=gofstat(list(P, nb), fitnames = c("Poisson", "neg.bin"))
  return(names(sort(a$aic))[1])
}
fn.show.mod.sel=function(MODS,outs)
{
  plot(MODS, type="s")
  plot.new()
  tmp <- weightable(MODS)
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = .65)),
    colhead = list(fg_params=list(cex = .8)),
    rowhead = list(fg_params=list(cex = .8)))
  myt <- gridExtra::tableGrob(tmp[1:outs,], theme = mytheme)
  grid.draw(myt)
}
fn.modl.sel=function(RESPNS)   #function for model structure selection
{
  PREDS[id.cov]=paste("LN",PREDS[id.cov],sep="")
  PREDS=PREDS[-match(always,PREDS)]
  if(RESPNS=="LNcpue")Formula=formula(paste("LNcpue",paste(PREDS,collapse="+"),sep="~"))
  if(RESPNS=="catch")Formula=formula(paste(Response,paste(paste(PREDS,collapse="+"),
                              paste("offset(","LN",efrt,")",sep=""),sep="+"),sep="~"))
  
  res <- glmulti(Formula,data=d,level=ifelse(Inter=="MainTerm",1,ifelse(Inter=="2way",2,"3way")),
                 method="h",fitfunction=fitFun,
                 always=paste('+',paste(always,collapse="+"),sep=""),
                 crit="aicc",confsetsize=2^length(PREDS),plotty=F,report=T)
  
  return(list(res=res,BEST=res@formulas[[1]]))
}

viz.coef=function(MOD,WHAT) coefplot(MOD,coefficient=WHAT)

fn.stand=function(d,Response,RESPNS,PREDS,efrt,Formula,Formula.gam)   #function for standardisation
{
  id.fctr=which(PREDS%in%Categorical)
  d=makecategorical(PREDS[id.fctr],d) 
  d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
  id.cov=which(PREDS%in%Covariates)
  d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
  res=NULL
  if(!is.null(Formula)) res <- glm(Formula,data=d)
  res.gam=NULL
  if(!is.null(Formula.gam)) res.gam <-gam(Formula.gam,data=d,method="REML")
  return(list(res=res,res.gam=res.gam,DATA=d))
}

fn.delta=function(d,Response,PREDS,efrt,Formula,Formula.gam)   #function for standardisation
{
  ALLvars=all.vars(Formula)[-1]
  Formula.bi=as.formula(paste('catch.pos',"~",paste(paste(ALLvars,collapse="+"),"offset(LNeffort)",sep="+")))
  id.fctr=which(PREDS%in%Categorical)
  
  Bi <- d %>%mutate(catch.pos=as.numeric(catch.target>0))
  TAB=table(Bi$catch.pos,Bi$finyear)
  TAB[TAB>0]=1
  drop.yrs=names(which(TAB[2,]==0))
  if(length(drop.yrs)>0)
  {
    Bi=Bi%>%filter(!finyear%in%drop.yrs)
    d=d%>%filter(!finyear%in%drop.yrs)
  }
  Bi=makecategorical(PREDS[id.fctr],Bi)
  Bi$LNeffort=log(Bi[,match(efrt,names(Bi))])
  
  d=d%>%filter(catch.target>0)
  d=makecategorical(PREDS[id.fctr],d)
  d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
  
  res=res_bi=NULL
  if(!is.null(Formula))
  {
    res_bi <- glm(Formula.bi, data=Bi, family="binomial", maxit=100)
    res <- glm(Formula,data=d)
    while(any(is.na(coef(res)))|any(is.na(coef(res_bi))))
    {
      na.coef=names(which(is.na(coef(res_bi))))
      if(isTRUE(any(grepl('blockx', na.coef))))
      {
        drp=substr(na.coef[grepl('blockx', na.coef)],nchar('blockx')+1,15)
        Bi=subset(Bi,!blockx%in%drp)%>%
          mutate(blockx=droplevels(blockx))
      }
      if(isTRUE(any(grepl('vessel', na.coef))))
      {
        drp=substr(na.coef[grepl('vessel', na.coef)],nchar('vessel')+1,15)
        Bi=subset(Bi,!vessel%in%drp)%>%
          mutate(vessel=droplevels(vessel))
      }
      if(isTRUE(any(grepl('month', na.coef))))
      {
        drp=substr(na.coef[grepl('vessel', na.coef)],nchar('month')+1,15)
        Bi=subset(Bi,!month%in%drp)%>%
          mutate(month=droplevels(month))
      }
      res_bi <- glm(Formula.bi, data=Bi, family="binomial", maxit=100)
      
      
      na.coef=names(which(is.na(coef(res))))
      if(isTRUE(any(grepl('blockx', na.coef))))
      {
        drp=substr(na.coef[grepl('blockx', na.coef)],nchar('blockx')+1,15)
        d=subset(d,!blockx%in%drp)%>%
          mutate(blockx=droplevels(blockx))
      }
      if(isTRUE(any(grepl('vessel', na.coef))))
      {
        drp=substr(na.coef[grepl('vessel', na.coef)],nchar('vessel')+1,15)
        d=subset(d,!vessel%in%drp)%>%
          mutate(vessel=droplevels(vessel))
      }
      if(isTRUE(any(grepl('month', na.coef))))
      {
        drp=substr(na.coef[grepl('vessel', na.coef)],nchar('month')+1,15)
        d=subset(d,!month%in%drp)%>%
          mutate(month=droplevels(month))
      }
      res <- glm(Formula,data=d)
    }
  }

  res.gam=res.gam_bi=NULL
  if(!is.null(Formula.gam))
  {
    ALLvars.gam=all.vars(Formula.gam)[-1]
    ALLvars.gam=ALLvars.gam[-match(c('long10.corner','lat10.corner'),ALLvars.gam)]
    Formula.bi.gam=as.formula(paste('catch.pos',"~",paste(paste(paste(ALLvars.gam,collapse="+"),
                                                                "s(long10.corner,lat10.corner)",sep='+'),"offset(LNeffort)",sep="+")))
    res.gam_bi <-gam(Formula.bi.gam,data=Bi, family="binomial",method="REML")
    res.gam <-gam(Formula.gam,data=d,method="REML")
  }
  
  return(list(res=res,res_bi=res_bi,res.gam=res.gam,res.gam_bi=res.gam_bi,DATA=d,DATA_bi=Bi))
}
fn.MC.delta.cpue=function(BiMOD,MOD,BiData,PosData,niter,pred.term,ALL.terms)
{
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
  newdata.bi=cbind(newdata.bi,LNeffort=mean(BiData$LNeffort))
  
  set.seed(999)
  
  Bi.pars.rand=rmvnorm(niter,mean=coef(BiMOD),sigma=Covar.bi)
  Pos.pars.rand=rmvnorm(niter,mean=coef(MOD),sigma=Covar.pos)
  
  MC.preds=matrix(nrow=niter,ncol=nrow(newdata.bi))
  for(n in 1:niter)
  {
    #Binomial part
    dummy.BiMOD$coefficients=Bi.pars.rand[n,]
    newdata.bi$Pred.bi=predict(dummy.BiMOD,newdata=newdata.bi, type="response")
    
    
    #Positive part
    dummy.MOD$coefficients=Pos.pars.rand[n,]
    a=predict(dummy.MOD,newdata=newdata.pos, type="response",se.fit=T)
    newdata.pos$Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    
    dummy=left_join(newdata.bi,newdata.pos,by=pred.term)%>%
      mutate(Index=Pred.bi*Pred)%>%
      select(Index)%>%
      unlist
    
    MC.preds[n,]=dummy
  }
  
  #Get summary stats
  MEAN=colMeans(MC.preds,na.rm=T)
  SD=apply(MC.preds,2,sd,na.rm=T)
  CV=SD/MEAN
  LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T))
  UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T))
  
  Stats=cbind(pred.dat.pos,data.frame(MEAN=MEAN,SD=SD,CV=CV,LOW=LOW,UP=UP))
  Stats=Stats[order(Stats[,1]),]
  Stats=Stats%>%rename(response=MEAN, lower.CL=LOW, upper.CL=UP)
  rownames(Stats)=NULL
  return(Stats)
}

Anova.and.Dev.exp=function(MOD,SP,type,gam.extra)   #function for extracting term significance and deviance explained
{
  #Anovas
  if(class(MOD)[1]=='glm')
  {
    Anova.tab=anova(MOD, test = "Chisq")
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/MOD$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    Dev.exp=sum(Term.dev.exp)
    ANOVA=as.data.frame.matrix(Anova.tab)
    Term=data.frame(Percent.dev.exp=Term.dev.exp)
    Table=ANOVA[-1,match(c("Df","Pr(>Chi)"),names(ANOVA))]
    Term=Term[match(rownames(Term), rownames(Table)),]
    Table=cbind(Table,Term)
    names(Table)[match("Term",names(Table))]="Percent.dev.exp"
    Table$term=rownames(ANOVA)[2:nrow(ANOVA)]
    Table=Table%>%select('term','Df','Pr(>Chi)','Percent.dev.exp')
    All=Table[1,]
    All[,1:ncol(All)]=NA
    All$term="model"
    All$Percent.dev.exp=round(Dev.exp,3)
    Table=rbind(Table,All)
    vars <- c(df = "Df", 'p-value' ="Pr(>Chi)")
    Table= Table %>% mutate_at(c("Pr(>Chi)","Percent.dev.exp"), round, 3) %>%
      rename(!!vars)
  }
  
  if(class(MOD)[1]=="gam")
  {
    Anova.tab=anova(MOD)
    ANOVA=as.data.frame(Anova.tab$pTerms.table[,-2])
    s.mat=as.data.frame(Anova.tab$s.table)
    s.mat=s.mat[,-(2:3)]
    names(s.mat)=names(ANOVA)
    ANOVA=rbind(ANOVA,s.mat)
    ANOVA$term=rownames(ANOVA)
    gamo=gam.extra
    for(l in 2:length(gam.extra)) gamo[l]=gam.extra[l]-gam.extra[l-1]
    gamo=100*gamo
    gam.dev.exp=data.frame(Percent.dev.exp=gamo,term=names(gam.extra))
    ANOVA=ANOVA%>%left_join(gam.dev.exp,by="term")
    
    ANOVA = ANOVA %>% select(term, df, 'p-value', Percent.dev.exp)
    model=ANOVA[1,]
    model[,]=NA
    model$Percent.dev.exp=sum(gamo)
    model$term='model'
    ANOVA=rbind(ANOVA,model)
    Table= ANOVA %>% mutate_at(c("p-value","Percent.dev.exp"), round, 3)  %>%
      mutate_at(c("df"), round, 0)
    
  }
  Table$"p-value"=ifelse(Table$"p-value"<0.001,"<0.001",Table$"p-value")
  dummy=Table[1:2,]
  dummy[,]=NA
  dummy$term=c(type,SP)
  Table=rbind(dummy,Table)
  Table[is.na(Table)] <- ""
  rownames(Table)=NULL
  return(Table)
}

BiasCor.fn=function(Median,sigma) biasCorr <- exp(Median+(sigma^2)/2) #function for bias corrected mean in normal space

#function for model predictions
#note: this issues marginal means with accounts for unbalanced data
#see https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
pred.fun=function(MOD,biascor,PRED,Pred.type)             
{
  lsm=summary(emmeans(MOD, PRED, type=Pred.type))
  id.low=which(names(lsm)%in%c("lower.CL","asymp.LCL"))
  id.up=which(names(lsm)%in%c("upper.CL","asymp.UCL"))
  if(biascor=="YES")
  {
    sigma.glm=sqrt(summary(MOD)$dispersion) # residuals standard error
    lsm$response=exp(lsm$emmean)*exp(sigma.glm^2/2)
    lsm$lower.CL=exp(lsm[,id.low])*exp(sigma.glm^2/2)
    lsm$upper.CL=exp(lsm[,id.up])*exp(sigma.glm^2/2)
  }
  if(biascor=="NO")
  {
    if(is.na(match("response",names(lsm))))lsm$response=lsm$emmean
    lsm$lower.CL=lsm[,id.low]
    lsm$upper.CL=lsm[,id.up]
  }
  return(lsm)
}

pred.fun.spatial=function(DAT,MOD,PRED,FORM,Spatial.grid)
{
  TermS=all.vars(FORM)
  TermS=TermS[-match(c("LNcpue",PRED),TermS)]
  id.cat=TermS[which(TermS%in%Categorical)]
  NewDat=as.data.frame(matrix(nrow=1,ncol=length(id.cat)))
  names(NewDat)=id.cat
  for(ii in 1:length(id.cat))
  {
    dummy=sort(table(DAT[,match(id.cat[ii],names(DAT))]))
    NewDat[,ii]=factor(names(dummy[length(dummy)]),levels(DAT[,match(id.cat[ii],names(DAT))]))
  }
  id.cont=TermS[which(!TermS%in%Categorical)]
  if(length(id.cont)>0)
  {
    NewDat.cont=as.data.frame(matrix(nrow=1,ncol=length(id.cont)))
    names(NewDat.cont)=id.cont
    for(ii in 1:length(id.cont)) NewDat.cont[,ii]=mean(DAT[,match(id.cont[ii],names(DAT))])
    NewDat=cbind(NewDat,NewDat.cont)
  }
  NewDat=cbind(Spatial.grid,NewDat)
  PRD=predict(MOD,newdata=NewDat,type='link',se.fit=T)
  NewDat$cpue=exp(PRD$fit+(PRD$se.fit^2)/2)
  return(NewDat)
}

#function for model predictions by zone
#note: ref_grid selects the subset of blocks
pred.fun.zone=function(MOD,Subset,biascor,PRED,Pred.type)             
{
  #REFGRD=ref_grid(MOD)
  #head(REFGRD@ grid)  #see grid
  #plot(REFGRD, by = "month")
  #emmip(REFGRD, finyear ~ month)
  REFGRD=ref_grid(MOD,at=list(blockx=Subset))
  lsm=summary(emmeans(REFGRD, PRED, type=Pred.type))
  if(biascor=="YES")
  {
    sigma.glm=sqrt(summary(MOD)$dispersion) # residuals standard error
    lsm$response=exp(lsm$emmean)*exp(sigma.glm^2/2)
    lsm$lower.CL=exp(lsm$asymp.LCL)*exp(sigma.glm^2/2)
    lsm$upper.CL=exp(lsm$asymp.UCL)*exp(sigma.glm^2/2)
  }
  if(biascor=="NO")
  {
    lsm$response=lsm$emmean
    lsm$lower.CL=lsm$asymp.LCL
    lsm$upper.CL=lsm$asymp.UCL
  }
  return(lsm)
}


Plot.cpue=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar,add.lines)    #plot cpues
{
  if(inherits(cpuedata, "list")) 
  {
    if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
    if(length(cpuedata)<=3)tc=seq(-.5*0.15,.5*0.15,length.out=length(cpuedata))
    ymax = max(unlist(lapply(cpuedata, `[`, "upper.CL")),na.rm=T)
    Yrs=as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4))
    plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
    if(COL=='color')CL=c("black","forestgreen", "red","bisque3","blue2","dodgerblue")
    if(COL=='grey') CL=gray.colors(length(cpuedata),start=0.2,end=0.65)
    for(l in 1:length(cpuedata))
    {
      aaa=cpuedata[[l]]
      aaa$finyear=as.character(aaa$finyear)
      msn=Yrs[which(!Yrs%in%as.numeric(substr(aaa$finyear,1,4)))]
      if(length(msn)>0)
      {
        ad=aaa[length(msn),]
        ad[,]=NA
        ad$finyear=msn
        aaa=rbind(aaa,ad)
        aaa=aaa[order(aaa$finyear),]
      }
      
      with(aaa,
           {
             if(add.lines=="NO") points(Yrs+tc[l], response, pch=16, lty=2, col=CL[l],cex=CxS)
             if(add.lines=="YES") points(Yrs+tc[l], response, "o", pch=16, lty=2, col=CL[l],cex=CxS)
             arrows(x0=Yrs+tc[l], y0=lower.CL, 
                    x1=Yrs+tc[l], y1=upper.CL, 
                    code=3, angle=90, length=0.05, col=CL[l])
           })
      if(ADD.LGND=="YES") legend(whereLGND,names(cpuedata),bty='n',pch=16,col=CL,cex=1.45)
    }
  }
  
  if(inherits(cpuedata, "data.frame"))
  {
    ymax = max(cpuedata$upper.CL,na.rm=T)
    Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
    plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
    CL="black"
    if(add.lines=="NO")points(Yrs, cpuedata$response, pch=16, lty=2, col=CL,cex=CxS)
    if(add.lines=="YES")points(Yrs, cpuedata$response, "o", pch=16, lty=2, col=CL,cex=CxS)
    arrows(x0=Yrs, y0=cpuedata$lower.CL, 
           x1=Yrs, y1=cpuedata$upper.CL, 
           code=3, angle=90, length=0.05, col=CL)
  }
}
Plot.cpue.other=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar,All.yrs)    #plot cpues other species
{
  if(is.null(cpuedata[[1]])) plot.new()
  if(!is.null(cpuedata[[1]]))
  {
    if(inherits(cpuedata, "list")) 
    {
      if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
      if(length(cpuedata)<=3)tc=seq(-.5*0.15,.5*0.15,length.out=length(cpuedata))
      ymax = max(unlist(lapply(cpuedata, `[`, "upper.CL")),na.rm=T)
      Yrs=as.numeric(substr(All.yrs,1,4))
      plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent")
      if(COL=='color')CL=c("black","forestgreen", "red","bisque3","blue2","dodgerblue")
      if(COL=='grey') CL=gray.colors(length(cpuedata),start=0.2,end=0.65)
      for(l in 1:length(cpuedata))
      {
        aaa=cpuedata[[l]]%>%mutate(finyear=as.character(finyear))
        aaa=data.frame(finyear=All.yrs)%>%left_join(aaa,by="finyear")
        with(aaa,
             {
               points(Yrs+tc[l], response, "o", pch=16, lty=2, col=CL[l],cex=CxS)
               arrows(x0=Yrs+tc[l], y0=lower.CL, 
                      x1=Yrs+tc[l], y1=upper.CL, 
                      code=3, angle=90, length=0.05, col=CL[l])
             })
        if(ADD.LGND=="YES") legend(whereLGND,names(cpuedata),bty='n',pch=16,col=CL,cex=1.45)
      }
    }
    if(inherits(cpuedata, "data.frame"))
    {
      ymax = max(cpuedata$upper.CL,na.rm=T)
      Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
      plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent")
      CL="black"
      points(Yrs, cpuedata$response, "o", pch=16, lty=2, col=CL,cex=CxS)
      arrows(x0=Yrs, y0=cpuedata$lower.CL, 
             x1=Yrs, y1=cpuedata$upper.CL, 
             code=3, angle=90, length=0.05, col=CL)
    }
  }
}

#functions for influence plots (Bentley et al 2012)
bubble.plot=function(x,y,z,scaler,Xlab,Ylab)  
{
  xo=outer(x,rep(1,length=length(y)))
  yo=t(outer(y,rep(1,length=length(x))))
  zo=z
  for(zz in 1:nrow(zo))zo[zz,]=((zo[zz,]/max(zo[zz,]))^0.5)*scaler
  matplot(xo,yo,type="n",xlab=Xlab,ylab=Ylab,xaxt='n',yaxt='n')
  for(s in 1:length(x))
  {
    points(xo[s,],yo[s,],cex=zo[,s],pch=16,col="grey80")
    points(xo[s,],yo[s,],cex=zo[,s],pch=1,col="black")
  }
}
Influence.fn=function(MOD,DAT,Term.type,termS,add.Influence,SCALER)  
{
  termS=subset(termS,termS%in%names(Term.type))
  #extract main term coefficients for each species
  nt=length(termS)
  Store1=Store2=MatcH=COEF.list=COEF.SE.list=vector('list',nt)
  ID=c(1,grep("[:]", names(coef(MOD))))
  Cofs=coef(MOD)[-ID]
  if(class(MOD)[1]=="glm")Cofs.SE=summary(MOD)$coefficients[-ID, 2]
  if(class(MOD)[1]=="gam")Cofs.SE=summary(MOD)$se[-ID]
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      Store1[[p]]=as.character(levels(DAT[,match(termS[p],names(DAT))]))[-1]
      Store2[[p]]=paste(termS[p],Store1[[p]],sep="")
    }
  }
  for(p in 1:nt)MatcH[[p]]=if (Term.type[p]=="CAT") match(Store2[[p]],names(Cofs))
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") COEF.list[[p]]=Cofs[MatcH[[p]]]
    if (Term.type[p]=="Cont") COEF.list[[p]]=Cofs[match(termS[p],names(Cofs))]
  }
  for(p in 1:nt) 
  {
    if (Term.type[p]=="CAT")
    {
      COEF.list[[p]]=data.frame(Dummy=Store1[[p]],coef=COEF.list[[p]])
      COEF.list[[p]]$Dummy=as.character(COEF.list[[p]]$Dummy)
    }
  }
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      A=as.character(levels(DAT[,match(termS[p],names(DAT))]))[1]
      COEF.list[[p]]=rbind(COEF.list[[p]],data.frame(Dummy=A,coef=0))
    }
  }
  for(p in 1:nt) if (Term.type[p]=="CAT")colnames(COEF.list[[p]])=c(termS[p],paste("Coef.",termS[p],sep=""))
  
  #attach coefficients to data
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") DAT=merge(DAT,COEF.list[[p]],by=termS[p],all.x=T)
    if (Term.type[p]=="Cont")
    {
      DAT=cbind(DAT,COEF.list[[p]]*DAT[,match(names(COEF.list[[p]]),names(DAT))]) #coef X value
      colnames(DAT)[ncol(DAT)]=paste("Coef.",termS[p],sep="")      
    }
  }
  
  #ny
  ny=table(DAT$finyear)
  
  #calculate mean of coefficient
  Coef.vec=match(paste("Coef.",termS,sep=""),names(DAT))
  Mean.coef=Annual.Dev=vector('list',nt)
  names(Annual.Dev)=termS
  Over.all.influence=rep(NA,nt)
  names(Over.all.influence)=termS
  for(p in 1:nt) Mean.coef[[p]]=mean(DAT[,Coef.vec[p]],na.rm=T)
  
  #calculate overall and annual deviation from mean (i.e. influence)
  #- Categorical variables
  for(p in 1:nt)
  {
    dev=rep(NA,length(ny))
    for(t in 1:length(ny))
    {
      a=subset(DAT,finyear==names(ny[t]))
      dev[t]=(sum(a[,Coef.vec[p]]-Mean.coef[[p]]))/ny[t]
    }  
    
    #Store Annual deviance
    #note: exp because it's multiplicative
    Annual.Dev[[p]]=exp(dev)  
    
    #Store Overall influence of variable
    Over.all.influence[p]=exp(sum(abs(dev))/length(ny))-1    
  }
  
  #plot CDI (categorical vars only)
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      COEF=Cofs[MatcH[[p]]]
      COEF.SE=Cofs.SE[MatcH[[p]]]
      COef.nm=names(COEF)
      COEF=sort(COEF)
      COef.nm.sorted=names(COEF)
      COEF.SE=COEF.SE[match(COef.nm.sorted,names(COEF.SE))]
      x=1:length(COEF)
      if(length(x)>1)
      {
        nf <- layout(matrix(c(1,1,1,2,2,2), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
        if(add.Influence=="YES")nf <- layout(matrix(c(1,1,0,2,2,3), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
        par(mar=c(0,0,0,0),oma=c(4,6,1,1),las=1,mgp=c(1,.9,0))
        #layout.show(nf)
        
        # Coefficients
        minSE=COEF-COEF.SE
        maxSE=COEF+COEF.SE
        plot(x,COEF,xlab="",xaxt="n",ylim=c(min(minSE),max(maxSE)),cex.axis=1.25,pch=19,cex=2)
        arrows(x, minSE, x, maxSE, code=3, angle=90, length=0.1)
        axis(1,1:length(COEF),F,tcl=0.5)
        axis(1,seq(1,length(COEF),2),F,cex.axis=1.15,tcl=1)      
        mtext("Coefficient",side=2,line=4,cex=1.5,las=3)
        
        # Bubble plot of records
        DAT[,match(termS[p],names(DAT))]=as.factor(DAT[,match(termS[p],names(DAT))])
        TAb=table(DAT$finyear,DAT[,match(termS[p],names(DAT))])
        TAb=TAb[,match(Store1[[p]],colnames(TAb))]
        Prop.Rec=TAb/rowSums(TAb)
        if(colnames(Prop.Rec)[1]=="1")Prop.Rec=Prop.Rec[,-1]
        Nombres=gsub("[^[:digit:]]", "", names(COEF))
        
        if(termS[p]=="vessel") Nombres=1:length(Nombres)   #change vessel name for dummy
        rownames(Prop.Rec)=names(ny[1:length(ny)])
        colnames(Prop.Rec)=COef.nm
        Prop.Rec=Prop.Rec[,match(COef.nm.sorted,colnames(Prop.Rec))]
        bubble.plot(x,1:length(ny),Prop.Rec,scaler=SCALER,termS[p],"Financial year")
        axis(1,1:length(COEF),F,tck=-0.015)
        axis(1,seq(1,length(COEF),1),Nombres[seq(1,length(COEF),1)],cex.axis=1.15,tck=-0.025)
        axis(2,1:length(ny),F,tck=-0.015)
        axis(2,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1,tck=-0.025)   
        mtext("Financial year",side=2,line=4.35,cex=1.5,las=3)
        mtext(termS[p],side=1,line=2.5,cex=1.5)
        
        if(add.Influence=="YES")
        {
          #Influence plot
          plot(Annual.Dev[[p]],1:length(ny),type="o",pch=19,xlab="",ylab="",cex=2,cex.axis=1.25,yaxt='n')
          abline(v=1,lty=3,col=1)
          mtext("Influence",side=1,line=2.5,cex=1.5)      
          axis(2,1:length(ny),F,tcl=0.5)
          axis(2,seq(1,length(ny),2),F,tcl=1)
        }
        
      }
    }
  }
  return(list(Annual.Dev=Annual.Dev,ny=ny,Over.all.influence=Over.all.influence))
}
Compare.term.infl.fun=function(A,WHERE,YLIM)
{
  Annual.Dev=A$Annual.Dev
  ny=A$ny
  if(is.null(YLIM))YLIM=c(min(unlist(lapply(Annual.Dev,min))),max(unlist(lapply(Annual.Dev,max))))
  NamE=names(A$Over.all.influence)
  
  nt=length(Annual.Dev)  
  LTY=c(1,4,3,1,3,2)
  plot(1:length(ny),Annual.Dev[[1]],col=LTY.col[1],type="l",xlab="",ylab="",lwd=LWD,
       cex.axis=1.35,xaxt='n',ylim=YLIM)
  abline(h=1,lty=3,col=1)
  axis(1,1:length(ny),F,tck=-0.02)
  axis(1,seq(1,length(ny),2),F,tck=-0.04)
  axis(1,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1.35,tck=-0.04)
  for(p in 2:nt)lines(1:length(ny),Annual.Dev[[p]],lwd=LWD,lty=LTY[p],col=LTY.col[p])
  LEG=paste(NamE," (",round(100*A$Over.all.influence,1),"%)",sep="")
  legend(WHERE,LEG,bty='n',lty=LTY,col=LTY.col,lwd=LWD,cex=1.25,pt.cex=1.5)
}
Fig.CDI.paper.fn=function(MOD,DAT,Term.type,termS,SCALER,YLABs,CxAx)  
{
  #extract main term coefficients for each species
  nt=length(termS)
  Store1=Store2=MatcH=COEF.list=COEF.SE.list=vector('list',nt)
  ID=c(1,grep("[:]", names(coef(MOD))))
  Cofs=coef(MOD)[-ID]
  Cofs.SE=summary(MOD)$coefficients[-ID, 2]
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      Store1[[p]]=as.character(levels(DAT[,match(termS[p],names(DAT))]))[-1]
      Store2[[p]]=paste(termS[p],Store1[[p]],sep="")
    }
  }
  for(p in 1:nt)MatcH[[p]]=if (Term.type[p]=="CAT") match(Store2[[p]],names(Cofs))
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") COEF.list[[p]]=Cofs[MatcH[[p]]]
    if (Term.type[p]=="Cont") COEF.list[[p]]=Cofs[match(termS[p],names(Cofs))]
  }
  for(p in 1:nt) 
  {
    if (Term.type[p]=="CAT")
    {
      COEF.list[[p]]=data.frame(Dummy=Store1[[p]],coef=COEF.list[[p]])
      COEF.list[[p]]$Dummy=as.character(COEF.list[[p]]$Dummy)
    }
  }
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      A=as.character(levels(DAT[,match(termS[p],names(DAT))]))[1]
      COEF.list[[p]]=rbind(COEF.list[[p]],data.frame(Dummy=A,coef=0))
    }
  }
  for(p in 1:nt) if (Term.type[p]=="CAT")colnames(COEF.list[[p]])=c(termS[p],paste("Coef.",termS[p],sep=""))
  
  #attach coefficients to data
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") DAT=merge(DAT,COEF.list[[p]],by=termS[p],all.x=T)
    if (Term.type[p]=="Cont")
    {
      DAT=cbind(DAT,COEF.list[[p]]*DAT[,match(names(COEF.list[[p]]),names(DAT))]) #coef X value
      colnames(DAT)[ncol(DAT)]=paste("Coef.",termS[p],sep="")      
    }
  }
  
  #ny
  ny=table(DAT$finyear)
  
  #calculate mean of coefficient
  Coef.vec=match(paste("Coef.",termS,sep=""),names(DAT))
  Mean.coef=Annual.Dev=vector('list',nt)
  names(Annual.Dev)=termS
  Over.all.influence=rep(NA,nt)
  names(Over.all.influence)=termS
  for(p in 1:nt) Mean.coef[[p]]=mean(DAT[,Coef.vec[p]],na.rm=T)
  
  #calculate overall and annual deviation from mean (i.e. influence)
  #- Categorical variables
  for(p in 1:nt)
  {
    dev=rep(NA,length(ny))
    for(t in 1:length(ny))
    {
      a=subset(DAT,finyear==names(ny[t]))
      dev[t]=(sum(a[,Coef.vec[p]]-Mean.coef[[p]]))/ny[t]
    }  
    
    #Store Annual deviance
    #note: exp because it's multiplicative
    Annual.Dev[[p]]=exp(dev)  
    
    #Store Overall influence of variable
    Over.all.influence[p]=exp(sum(abs(dev))/length(ny))-1    
  }
  
  #plot CDI (categorical vars only)
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      COEF=Cofs[MatcH[[p]]]
      COEF.SE=Cofs.SE[MatcH[[p]]]
      COef.nm=names(COEF)
      COEF=sort(COEF)
      COef.nm.sorted=names(COEF)
      COEF.SE=COEF.SE[match(COef.nm.sorted,names(COEF.SE))]
      x=1:length(COEF)
      
      # Coefficients
      minSE=COEF-COEF.SE
      maxSE=COEF+COEF.SE
      plot(x,COEF,xlab="",ylab="",xaxt="n",ylim=c(min(minSE),max(maxSE)),cex.axis=CxAx,
           pch=19,cex=.9)
      arrows(x, minSE, x, maxSE, code=3, angle=90, length=0.1,lwd=.75)
      axis(1,1:length(COEF),F,tcl=0.25)
      axis(1,seq(1,length(COEF),2),F,tcl=.5)  
      
      # Bubble plot of records
      DAT[,match(termS[p],names(DAT))]=as.factor(DAT[,match(termS[p],names(DAT))])
      TAb=table(DAT$finyear,DAT[,match(termS[p],names(DAT))])
      TAb=TAb[,match(Store1[[p]],colnames(TAb))]
      Prop.Rec=TAb/rowSums(TAb)
      if(colnames(Prop.Rec)[1]=="1")Prop.Rec=Prop.Rec[,-1]
      Nombres=gsub("[^[:digit:]]", "", names(COEF))
      
      if(termS[p]=="vessel") Nombres=1:length(Nombres)   #change vessel name for dummy
      rownames(Prop.Rec)=names(ny[1:length(ny)])
      colnames(Prop.Rec)=COef.nm
      Prop.Rec=Prop.Rec[,match(COef.nm.sorted,colnames(Prop.Rec))]
      bubble.plot(x,1:length(ny),Prop.Rec,scaler=SCALER,"","")
      axis(1,1:length(COEF),F,tck=-0.015)
      axis(1,seq(1,length(COEF),1),Nombres[seq(1,length(COEF),1)],cex.axis=CxAx,tck=-0.025)
      axis(2,1:length(ny),F,tck=-0.015)
      axis(2,seq(1,length(ny),2),substr(names(ny)[seq(1,length(ny),2)],1,4),cex.axis=CxAx,tck=-0.025)   
      mtext(YLABs[p],side=1,line=1.75,cex=1)
    }
  }
  
}

Plot.cpue.spatial=function(cpuedata,var)
{
  if(var[1]=='blockx')cpuedata=cpuedata%>%mutate( Lat=-round(as.numeric(substr(get(var),1,2)),2),
                                                  Long=round(100+as.numeric(substr(get(var),3,4)),2))else
                                                    cpuedata=cpuedata%>%mutate( Lat=round(get(var[2]),2),
                                                                                Long=round(get(var[1]),2))
                                                  
                                                  YLIM=floor(range(Full.lat))    
                                                  XLIM=floor(range(Full.long)) 
                                                  
                                                  misn.lat=sort(Full.lat[which(!Full.lat%in%unique(cpuedata$Lat))])
                                                  misn.lon=sort(Full.long[which(!Full.long%in%unique(cpuedata$Long))])
                                                  if(length(misn.lat)>0 | length(misn.lon)>0)
                                                  {
                                                    if(var[1]=='blockx')
                                                    {
                                                      combo=expand.grid(Lon=Full.long,Lat=Full.lat)%>%
                                                        mutate(blockx=paste(abs(Lat),Lon-100,sep=''))%>%
                                                        select(blockx)
                                                      cpuedata=cpuedata%>%mutate(blockx=as.character(blockx))
                                                      cpuedata=combo%>%left_join(cpuedata,by=var)%>%
                                                        mutate( Lat=-round(as.numeric(substr(get(var),1,2)),2),
                                                                Long=round(100+as.numeric(substr(get(var),3,4)),2))
                                                      
                                                    }else
                                                    {
                                                      combo=expand.grid(Long=Full.long,Lat=Full.lat)
                                                      cpuedata=combo%>%left_join(cpuedata,by=c('Long','Lat'))
                                                    }
                                                  }
                                                  cpuedata=cpuedata%>%select(c(cpue,Lat,Long)) 
                                                  cpuedata.spread=cpuedata%>%spread(Lat,cpue)
                                                  Lon=as.numeric(cpuedata.spread$Long)
                                                  cpuedata.spread=as.matrix(cpuedata.spread[,-1]) 
                                                  LaT=as.numeric(colnames(cpuedata.spread))
                                                  brk<- quantile( c(cpuedata.spread),probs=seq(0,1,.1),na.rm=T)
                                                  YLIM[1]=YLIM[1]-0.5
                                                  YLIM[2]=YLIM[2]+0.5
                                                  XLIM[1]=XLIM[1]-0.5
                                                  XLIM[2]=XLIM[2]+0.5
                                                  image.plot(Lon,LaT,cpuedata.spread, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
                                                             lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
}

Pos.Diag.fn=function(MODEL,SPECIES,M)   #function for positive catch diagnostics
{
  RES=MODEL$residuals   #residuals
  Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
  PRED=predict(MODEL)
  
  qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="",xlab="")
  qqline(RES, col = 'grey40',lwd=1.5,lty=2)
  mtext(SPECIES,3,outer=F,line=0.25,cex=1.3)
  if(s==1) mtext("Residuals",2,outer=F,line=2,las=3,cex=M)
  if(s==2) mtext("                        Theoretical quantiles",1,outer=F,line=1.5,cex=M)
  
  hist(Std.RES,xlim=c(-5,5),ylab="",xlab="",main="",col="grey",breaks=50)
  box()
  if(s==1) mtext("Frequency",2,outer=F,line=2.5,las=3,cex=M)
  if(s==2) mtext("                      Standardised residuals",1,outer=F,line=1.5,cex=M)
  
  plot(PRED,Std.RES,ylab="",xlab="",ylim=c(-5,5))
  abline(0,0,lwd=1.5,lty=2,col='grey40')
  if(s==1) mtext("Standardised residuals",2,outer=F,line=2,las=3,cex=M)
  if(s==2) mtext("                         Fitted values",1,outer=F,line=1.5,cex=M)
  
}

fn.table.terms=function(d,PREDS)
{
  id.fctr=which(PREDS%in%Categorical)
  d=makecategorical(PREDS[id.fctr],d)
  
  store=vector('list',length(PREDS))
  names(store)=PREDS
  for(p in 1:length(store))
  {
    xx=d[,match(PREDS[p],names(d))]
    if(is.factor(xx))
    {
      levels=length(levels(xx))
      type="Categorical"
    }
    if(!is.factor(xx))
    {
      Unik=1
      type="Continuous"
      levels=""
    }
    a=data.frame(Term=PREDS[p],Type=type,Levels=levels)
    store[[p]]=a
  }
  return(do.call(rbind,store))
}


##############--- 4. PROCEDURE SECTION ---###################

#Correct weird Blocks and Reset BLOCKX to 4 digits as some have 5 digits
Data.monthly.GN$BLOCKX=with(Data.monthly.GN,ifelse(BLOCKX>50000 & !(is.na(LAT) & is.na(LONG)),
                                        paste(abs(floor(LAT)),floor(LONG)-100,0,sep=""),BLOCKX))
Data.daily.GN$BLOCKX=with(Data.daily.GN,ifelse(BLOCKX>50000 & !(is.na(LAT) & is.na(LONG)),
                                        paste(abs(floor(LAT)),floor(LONG)-100,0,sep=""),BLOCKX))
Effort.monthly$BLOCKX=with(Effort.monthly,ifelse(BLOCKX>50000 & !(is.na(LAT) & is.na(LONG)),
                                        paste(abs(floor(LAT)),floor(LONG)-100,0,sep=""),BLOCKX))
Effort.daily$blockx=with(Effort.daily,ifelse(blockx>50000 & !(is.na(LAT) & is.na(LONG)),
                                        paste(abs(floor(LAT)),floor(LONG)-100,0,sep=""),blockx))

Data.monthly.GN$BLOCKX=as.integer(substr(Data.monthly.GN$BLOCKX,1,4))
Effort.monthly$BLOCKX=as.integer(substr(Effort.monthly$BLOCKX,1,4))
Data.daily.GN$BLOCKX=as.integer(substr(Data.daily.GN$BLOCKX,1,4))
Effort.daily$blockx=as.integer(substr(Effort.daily$blockx,1,4))


#...Effort 
  # -- Daily                
Effort.daily$LAT=-abs(Effort.daily$LAT)
Block.lat.long=Effort.daily[!duplicated(Effort.daily$blockx),match(c("blockx","LAT","LONG"),names(Effort.daily))]

#km gn days  
Eff.daily.c.daily=aggregate(Km.Gillnet.Days.c~Same.return.SNo+vessel+finyear+month+blockx+block10,
                            data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)    
Eff.daily.daily=aggregate(Km.Gillnet.Days.inv~Same.return.SNo+vessel+finyear+month+blockx+block10,
                          data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)

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
Block.lat.long=Effort.monthly[!duplicated(Effort.monthly$BLOCKX),match(c("BLOCKX","LAT","LONG"),names(Effort.monthly))]

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
a=Effort.daily[,match(c("blockx","block10","Same.return.SNo","vessel","Eff.Reporter"),names(Effort.daily))]
a=a[!duplicated(paste(a$blockx,a$block10,a$Same.return.SNo,a$vessel)),]
Eff.daily.c.daily=merge(Eff.daily.c.daily,a,by=c("blockx","block10","Same.return.SNo","vessel"),all.x=T)

#monthly
a=Effort.monthly[,match(c("BLOCKX","FINYEAR","MONTH","VESSEL","Eff.Reporter"),names(Effort.monthly))]
a=a[!duplicated(paste(a$BLOCKX,a$FINYEAR,a$MONTH,a$VESSEL)),]
Eff.monthly.c=merge(Eff.monthly.c,a,by=c("BLOCKX","FINYEAR","MONTH","VESSEL"),all.x=T)

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
#note:only use meshes of 6.5 inch (165 mm) and 7 inch (178 mm) in cpue standardisation
#Monthly
Mesh.monthly$mesh=Mesh.monthly$MSHIGH
Mesh.monthly$mesh=with(Mesh.monthly,ifelse(mesh==0,NA,mesh))
#a=subset(Mesh.monthly,mesh%in%c("165","178"))
#A=with(a,table(VESSEL,mesh))  #vessels fishing with a unique mesh so extrapolate to all records for same vessel
#a=a[!duplicated(a$VESSEL),match(c("VESSEL","mesh"),names(a))]
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

#Put data into species list (keep species with at least N.keep years of catches>Min.kg)
A=Data.monthly.GN%>%filter(SPECIES%in%Shark.species) %>%
  group_by(SPECIES,FINYEAR) %>%
  summarise(Weight=round(sum(LIVEWT.c))) %>%
  spread(FINYEAR, Weight)%>%
  data.frame()
A=A %>% mutate_at(.vars = names(A)[-match("SPECIES",names(A))], function(x)(ifelse(x>=Min.kg, 1, 0)))
Anm=A$SPECIES
A=rowSums(A[,-1],na.rm=T)
names(A)=Anm
SpiSis=unique(Data.monthly.GN%>%mutate(RSCommonName=ifelse(SPECIES==8001,'Greynurse Shark',RSCommonName))%>%
                         filter(SPECIES%in%as.numeric(names(A[A>=N.keep])))%>%
                         select(SPECIES,RSCommonName) %>%
                         filter(!RSCommonName==""))
SpiSis=SpiSis[!duplicated(SpiSis$SPECIES),]
nms=SpiSis$RSCommonName
SpiSis=SpiSis$SPECIES
names(SpiSis)=nms
SpiSis=SpiSis[-match(c(22999),SpiSis)]
SP.list=as.list(SpiSis)

#remove thresher and dogfish because there is not enough positive record data to estimate glm coefficients
SP.list=SP.list[-match(c("Thresher Shark","Gulper sharks, Sleeper Sharks & Dogfishes"),names(SP.list))]

#remove blacktips because only a few years of daily available
SP.list=SP.list[-match(c("Blacktip Shark"),names(SP.list))]

#combine dusky and bronzy (Rory request due to species id) and sawsharks general with common saw
SP.list=SP.list[-match(c("Common Sawshark","SawShark",'Dusky Whaler','Bronze Whaler'),names(SP.list))]  

SP.list$'Dusky Whaler Bronze Whaler'=c(18003,18001)
SP.list$'Sawsharks'=c(23900,23002)
SP.list$All.Non.indicators=Shark.species[-match(Indicator.sp,Shark.species)]    
nnn=1:length(SP.list)


#get core area (90% of catch)
core.per=90
Core=SP.list
pdf('C:/Matias/Analyses/Catch and effort/species core areas/cores.pdf')
for(s in nnn)
{
  d=subset(Data.monthly.GN,SPECIES%in%SP.list[[s]])
  Nm=unique(d$SPECIES)
  Nm=ifelse(length(Nm)>3,'others',Nm)
  d=aggregate(LIVEWT.c~LAT+LONG,d,sum)
  d=d[order(-d$LIVEWT.c),]
  d$CumSum=cumsum(d$LIVEWT.c)
  d$CumSum=100*d$CumSum/max(d$CumSum)
  plot(d$LONG,d$LAT,cex=fn.scale(d$LIVEWT.c,4),pch=19,col="steelblue",
       ylab="Lat",xlab="Long",main=paste(Nm,names(SP.list)[s]))
  d=subset(d,CumSum<=core.per)
  Rnglat=range(d$LAT)
  Rnglon=range(d$LONG)
  polygon(c(Rnglon[1],Rnglon[2],Rnglon[2],Rnglon[1]),
          c(Rnglat[1],Rnglat[1],Rnglat[2],Rnglat[2]),border=2)
  Core[[s]]=list(Lat=Rnglat,Long=Rnglon)
}
dev.off()

#adjust core areas following McAuley and Simpfendorfer
Dusky.range=c(-28,120)
Core$"Dusky Whaler Bronze Whaler"$Lat[2]=Dusky.range[1]
Core$"Dusky Whaler Bronze Whaler"$Long[2]=Dusky.range[2]
Sandbar.range=c(-26,118)
Core$"Sandbar Shark"$Long[2]=Sandbar.range[2]
Whiskery.range=c(-28,129)
Core$"Whiskery Shark"$Lat[2]=Whiskery.range[1]
Core$"Whiskery Shark"$Long[2]=Whiskery.range[2]
Gummy.range=c(116,129)
Core$"Gummy Shark"$Long=Gummy.range


#put date back in Daily data set
get.dates=subset(Effort.daily,Same.return.SNo%in%unique(Data.daily.GN$Same.return.SNo),select=c(Same.return.SNo,date))
get.dates=get.dates[!duplicated(get.dates$Same.return.SNo),]
Data.daily.GN=merge(Data.daily.GN,get.dates,"Same.return.SNo",all.x=T)
Data.daily.GN$date=as.Date(Data.daily.GN$date)

#create some useful vars
Post.yrs=max(unique(sort(Data.monthly.GN$YEAR.c)))
Post.yrs=paste(2006:(Post.yrs-1),substr(2007:Post.yrs,3,4),sep="-")
Daily.l.years=sort(unique(Data.daily.GN$FINYEAR))

  #Monthly
# Set all records of reapportioned returns to bad
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


#Create Km.Gillnet.Hours_shot.c
Eff$Km.Gillnet.Hours_shot.c=with(Eff,Km.Gillnet.Hours.c*SHOTS.c)

  #Daily
Data.daily.GN$Reporter.old="good"

# Set all records of reapportioned returns to bad
Baddies=subset(Data.daily.GN,Reporter.old=="bad")   
if(nrow(Baddies)>0)
{
  Baddies=unique(Baddies$Same.return)
  Data.daily.GN$Reporter.old=with(Data.daily.GN,ifelse(Same.return%in%Baddies,'bad',Reporter.old))
}
Data.daily.GN$Reporter=Data.daily.GN$Reporter.old 

#remove small net length which correspond to non-shark gillnet
Data.daily.GN=subset(Data.daily.GN,netlen.c >100)

#deal with boundary blocks in zone 1-2
if(BOUND.BLK=="REMOVE")Data.daily.GN=subset(Data.daily.GN,!BLOCKX%in%Boundary.Blks)
if(BOUND.BLK=="REALLOCATE")Data.daily.GN$zone=with(Data.daily.GN,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))

Data.daily.GN$Same.return=with(Data.daily.GN,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))

#Create Km.Gillnet.Hours_shot.c
Eff.daily$Km.Gillnet.Hours_shot.c=with(Eff.daily,Km.Gillnet.Hours.c*shots.c)


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
                  group_by(MONTH,YEAR.c)%>%
                  mutate(Temp.res=Temperature/mean(Temperature,na.rm=T))

#Daily
Data.daily.GN=Data.daily.GN %>%
                  mutate(LAT.round=-floor(abs(LAT)),LONG.round=floor(LONG)) %>%
                  left_join(SST,by=c("YEAR.c"="year","MONTH"="month","LONG.round"="Long","LAT.round"="Lat"))%>% 
                  arrange(YEAR.c,MONTH,LAT.round,LONG.round)%>%
                  mutate(Temperature=ifelse(is.na(Temperature),na.approx(Temperature),Temperature))%>%
                  group_by(MONTH,YEAR.c)%>%
                  mutate(Temp.res=Temperature/mean(Temperature,na.rm=T))%>%
                  select(-c(LONG.round,LAT.round))
  


# Add SOI, Freo and Moon (the later to daily only) 
Freo=Freo%>%rename(Freo=MeanSeaLevel)%>%
            mutate(Freo_lag6=lag(Freo,6),
                   Freo_lag12=lag(Freo,12))

  #Monthly
Data.monthly.GN=Data.monthly.GN%>%left_join(SOI,by=c("YEAR.c"="Year","MONTH"="Month"))%>%
                                   left_join(Freo,by=c("YEAR.c"="Year","MONTH"="Month")) 
  


  #Daily
Data.daily.GN=Data.daily.GN%>%left_join(SOI,by=c("YEAR.c"="Year","MONTH"="Month"))%>%
                              left_join(Freo,by=c("YEAR.c"="Year","MONTH"="Month")) %>%
                              mutate(Lunar=lunar.phase(date,name=T))




#Create species data sets  

  #Monthly
#note: select species range, and add effort by Same return
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
#getDoParWorkers()
system.time({Species.list=foreach(s=nnn,.packages=c('dplyr','doParallel')) %dopar%
  {
    return(fn.cpue.data(Dat=Data.monthly.GN %>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                         LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                        EffrrT=Eff%>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                               LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                        sp=SP.list[[s]]))
  }
})
names(Species.list)=names(SP.list) 
stopCluster(cl)

  #Daily 
#note: select species range and add effort by date or ID (==Same.return.SNo). Note that for catch aggregating by date
#       or by ID makes no difference but it's needed for merging with effort
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
system.time({Species.list.daily=foreach(s=nnn,.packages=c('dplyr','doParallel')) %dopar%
  {
    return(fn.cpue.data.daily(Dat=Data.daily.GN %>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                             LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                              EffrrT=Eff.daily%>% filter(LAT>=Core[[s]]$Lat[1] & LAT<=Core[[s]]$Lat[2] &
                                                           LONG>=Core[[s]]$Long[1]& LONG<=Core[[s]]$Long[2]),
                              sp=SP.list[[s]]))
  }
})
names(Species.list.daily)=names(SP.list) 
stopCluster(cl)

#Remove variables not used after prelim analysis
for(s in nnn)
{
  if(!is.null(Species.list[[s]])) Species.list[[s]] = Species.list[[s]] %>%  select(-c(LIVEWT,Boundary.blk,Km.Gillnet.Hours_shot.c,
                              TYPE.DATA,Sch.or.DogS,SOI,Freo_lag6,Freo_lag12,mesh,
                              NETLEN.c, BDAYS.c,HOURS.c))
  if(!is.null(Species.list.daily[[s]])) Species.list.daily[[s]] = Species.list.daily[[s]] %>%  select(-c(LIVEWT,Km.Gillnet.Days.inv,
                                Km.Gillnet.Hours.inv,Km.Gillnet.Hours_shot.c,netlen.c,hours.c,
                                bdays.c,TYPE.DATA,LIVEWT,nfish,SOI,Freo_lag6,Freo_lag12))
}

#Keep vessel characteristics from vessel survey for vessels that have fished      
if(exists("TDGDLF.survey"))  TDGDLF.survey=subset(TDGDLF.survey, BOATREGO%in%
                                    unique(c(Data.monthly.GN$VESSEL,Data.daily.GN$VESSEL)))   


#4.2 Extract number of vessels reporting catch of species per species range
N.VES=matrix(rep(NA,length(nnn)),ncol=length(nnn))
colnames(N.VES)=names(SP.list)
for(i in nnn)
{
  a=b=NULL
  if(!is.null(Species.list[[i]]))a=unique(subset(Species.list[[i]],SPECIES%in%SP.list[[i]])$VESSEL)
  if(!is.null(Species.list.daily[[i]]))b=unique(subset(Species.list.daily[[i]],SPECIES%in%SP.list[[i]])$VESSEL)
  N.VES[,i]=length(unique(c(a,b)))
}
setwd('C:/Matias/Analyses/Catch and effort')
hndl=paste(getwd(),"/Outputs/Paper/",sep="")
write.csv(N.VES,paste(hndl,"All.Vessels.by.species.csv",sep=""),row.names=T)


#4.2.1 Extract number of blocks where shark has been caught within effective area
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


#4.2.2 Block weights     
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


#4.5 Create useful vars
FINYEAR.monthly=as.character(unique(Data.monthly.GN$FINYEAR))
FINYEAR.monthly=sort(FINYEAR.monthly)
N.yrs=length(FINYEAR.monthly)

FINYEAR.daily=as.character(unique(Data.daily.GN$FINYEAR))
FINYEAR.daily=sort(FINYEAR.daily)
N.yrs.daily=length(FINYEAR.daily) 

FINYEAR.ALL=unique(c(FINYEAR.monthly,FINYEAR.daily))
FINYEAR.ALL=sort(FINYEAR.ALL)
N.yrs.ALL=length(FINYEAR.ALL)


#Export data for spatial distribution paper
write.csv(Data.daily.GN,'C:/Matias/Analyses/Catch and effort/Data_outs/Data.daily.GN_for_spatial_analysis.csv',row.names=F)

#4.6 Proportion of dusky and copper shark
fn.fig("proportion of dusky and copper shark_TDGLDF",2000,2400)
par(mfcol=c(2,1),las=1,mai=c(.8,.85,.1,.1),mgp=c(2.5,.8,0))
  #monthly
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Species.list[[match("Dusky Whaler Bronze Whaler",names(Species.list))]],SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Species.list[[match("Dusky Whaler Bronze Whaler",names(Species.list))]],SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="",xlab="",pch=19,col=2,cex=1.75,cex.lab=1.5,ylim=c(0,.25))
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)
legend("topright","Monthly returns",bty='n',cex=1.5)

  #daily
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Species.list.daily[[match("Dusky Whaler Bronze Whaler",names(Species.list))]],SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Species.list.daily[[match("Dusky Whaler Bronze Whaler",names(Species.list))]],SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="",xlab="Financial year",pch=19,col=2,cex=1.75,cex.lab=1.5,ylim=c(0,.25))
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)
legend("topright","Daily logbooks",bty='n',cex=1.5)

mtext("Bronze whaler shark catch / Dusky shark catch",2,-1.5,las=3,outer=T,cex=1.75)
dev.off()


#4.7 Determine indicative vessels and blocks        
#steps: 1. select vessels that meet criteria (fishing for at least Threshold.n.yrs/Threshold.n.yrs.daily
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
    #monthly
    if(!is.null(Species.list[[i]]))
    {
      finy=sort(unique(Species.list[[i]]$FINYEAR))
      Ves.sel.BC=0
      BLK.sel.BC=0
      Min.ktch=1
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
                                    NM=names(SP.list)[i],what=".monthly",
                                    Ves.sel.BC=Ves.sel.BC,Ves.sel.sens=Threshold.n.yrs.sens,
                                    BLK.sel.BC=BLK.sel.BC,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=Min.ktch)
      if(!is.null(dummy))
      {
        BLKS.used[[i]]=dummy$Blks.BC
        BLKS.not.used[[i]]=dummy$Drop.blks
        VES.used[[i]]=dummy$Ves.BC
        VES.not.used[[i]]=dummy$Drop.ves
      }
      
    }
    #Daily
    if(!is.null(Species.list.daily[[i]]))
    {
      Ves.sel.BC=0
      BLK.sel.BC=0
      Min.ktch=1
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
                                    NM=names(SP.list)[i],what=".daily",
                                    Ves.sel.BC=Ves.sel.BC,Ves.sel.sens=Threshold.n.yrs.sens,
                                    BLK.sel.BC=BLK.sel.BC,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=Min.ktch)
      if(!is.null(dummy))
      {
        BLKS.used.daily[[i]]=dummy$Blks.BC
        BLKS.not.used.daily[[i]]=dummy$Drop.blks
        BLKS_10.used.daily[[i]]=dummy$Blks.BC_10
        BLKS_10.not.used.daily[[i]]=dummy$Drop.blks_10
        VES.used.daily[[i]]=dummy$Ves.BC
        VES.not.used.daily[[i]]=dummy$Drop.ves
      }
    }
  }
}

# Illustrate folly effect: mean Vs sum
if(Show.folly.eg=="YES")
{
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper/MeanVsSum")
  
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
  
  #5.1.2 dropping record based on cpue threshols
  
  
  #need to combine monthly and daily to re run...
  # tiff(file="Folly.Sum.Mean.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  # par(mfrow=c(4,3),las=1,mai=c(.4,.5,.2,.01),oma=c(2,.1,.5,3),mgp=c(1,.7,0))
  # fn.plot.folly(Species.list[[1]],17003,CPUE.threshold=100,CPUE.threshold1=10)
  # mtext("Whiskery",4,las=3,cex=1.5,line=1)
  # 
  # fn.plot.folly(Species.list[[2]],17001,CPUE.threshold=100,CPUE.threshold1=10)
  # mtext("Gummy",4,las=3,cex=1.5,line=1)
  # 
  # fn.plot.folly(Species.list[[3]],18003,CPUE.threshold=100,CPUE.threshold1=10)
  # mtext("Dusky",4,las=3,cex=1.5,line=1)
  # 
  # fn.plot.folly(Species.list[[4]],18007,CPUE.threshold=25,CPUE.threshold1=5)
  # mtext("CPUE (kg/km.gn.day)",2,outer=T,las=3,cex=1.5,line=-1.5)
  # mtext("Financial year",1,outer=T,cex=1.5)
  # mtext("Sandbar",4,las=3,cex=1.5,line=1)
  # dev.off()
  
}

#Show example of variability in catch rates and need for standardisation
if(Show.variability.cpue.eg=="YES")
{
  setwd("C:/Matias/Analyses/Catch and effort/Outputs")
  
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


#4.10 Construct wide database for analysis 
#steps: 
#   1. select "Good" records (the variable "Reporter" includes good/bad catch and effort reporters)
#   2. Construct a single row for each record (i.e. 'year-month-vessel-block-gear' for monthly
#      returns and 'year-Session-vessel-block10-gear' for daily logbooks), with catch of target
#      and other species as separate columns, giving a 0 catch for column "target" if no catch

DATA.list.LIVEWT.c=vector('list',length(SP.list)) 
names(DATA.list.LIVEWT.c)=names(SP.list)
DATA.list.LIVEWT.c.daily=DATA.list.LIVEWT.c
Prop.Catch=rep(NA,length(SP.list))
names(Prop.Catch)=names(SP.list)
Prop.Catch.daily=Prop.Catch

  #monthly  
These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Hours.c","Km.Gillnet.Days.c",
                "zone","MONTH","BLOCKX","SHOTS.c")
for(i in nnn)
{
  if(!is.null(Species.list[[i]]))
  {
    dummy=Effort.data.fun(subset(Species.list[[i]],Reporter=="good"),SP.list[[i]],"LIVEWT.c")  
    DATA.list.LIVEWT.c[[i]]=dummy$dat   
    Prop.Catch[i]=dummy$prop.with.catch 
  }
}


  #daily 
These.efforts.daily=c("FINYEAR","date","TSNo","Km.Gillnet.Hours.c","Km.Gillnet.Days.c",
                      "zone","MONTH","BLOCKX","block10","VESSEL","mesh",
                      "Same.return.SNo","nlines.c","shots.c")
for(i in nnn)
{
  if(!is.null(Species.list.daily[[i]]))
  {
     # #set Reporter to Bad if average weight is outside species range
    # #note: some records have dodgy nfish so class as 'bad reporters'
    # Species.list.daily[[i]]$Avrg.w=Species.list.daily[[i]]$LIVEWT.c/Species.list.daily[[i]]$nfish
    # Species.list.daily[[i]]$Reporter=with(Species.list.daily[[i]],
    #                                       ifelse((Avrg.w>Max.weight[i]| Avrg.w<Min.weight[i]) & SPECIES==SPvec[i],"bad",Reporter))                    
    dummy=Effort.data.fun.daily(subset(Species.list.daily[[i]],Reporter=="good"),SP.list[[i]],
                                ktch="LIVEWT.c",Aggregtn="SNo")
    DATA.list.LIVEWT.c.daily[[i]]=dummy$dat
    Prop.Catch.daily[i]=dummy$prop.with.catch   
  }
}

#Export proportions with Catch
write.csv(Prop.Catch,paste(hndl,"Prop.records.with.catch.monthly.csv",sep=""),row.names=T)
write.csv(Prop.Catch.daily,paste(hndl,"Prop.records.with.catch.daily.csv",sep=""),row.names=T)


#Calculate block corners for gam
for(s in nnn)
{
  if(!is.null(DATA.list.LIVEWT.c[[s]])) DATA.list.LIVEWT.c[[s]] = DATA.list.LIVEWT.c[[s]] %>%  
                                                  mutate(LAT10.corner=LAT, LONG10.corner=LONG)
      
  if(!is.null(DATA.list.LIVEWT.c.daily[[s]])) DATA.list.LIVEWT.c.daily[[s]] = DATA.list.LIVEWT.c.daily[[s]] %>% 
                mutate(LAT10.corner=-(abs(as.numeric(substr(block10,1,2))+10*(as.numeric(substr(block10,3,3)))/60)),
                       LONG10.corner=100+as.numeric(substr(block10,4,5))+10*(as.numeric(substr(block10,6,6)))/60)
}

#Define target species index
Tar.sp=match(TARGETS,SP.list)


#4.12 Drop first years of sandbar data because vessels don't meet selection criteria and no positive catch
DD=subset(DATA.list.LIVEWT.c$"Sandbar Shark",BLOCKX%in%as.numeric(BLKS.used$"Sandbar Shark"))      
DD=subset(DD,VESSEL%in%VES.used$"Sandbar Shark")
San.Yrs=fn.sel.yrs.used(DD,ThrShld.n.vess=5)
rm(DD)
DATA.list.LIVEWT.c$"Sandbar Shark"=subset(DATA.list.LIVEWT.c$"Sandbar Shark",FINYEAR%in%San.Yrs)


#4.14  Identify targeting behaviour   (more in 2.CPUE standardisations_delta.R)
if(do_cluster=="YES")
{
  #clustering analysis
  fn.cluster=function(data,TarSp,target,varS,scaling,check.clustrbl,n.clus)
  {
    a=data[[TarSp]]%>% column_to_rownames(var = "Same.return.SNo")%>%
      select(varS[-match(target,varS)])
    if(scaling=="YES")a=scale(a)
    
    #step 1. Define if data are clusterable
    if(check.clustrbl=="YES")
    {
      #random sample to reduce computation time
      ran.samp=sample(1:nrow(a),15000,replace=F)
      
      res <- get_clust_tendency(a[ran.samp,], n = nrow(a[ran.samp,])-1, graph = FALSE)
      if(1-res$hopkins_stat>0.75) clusterable="YES"else  clusterable="NO"
      print(clusterable)
    }
    
    #step 2. Determine optimum number of clusters
    if(check.clustrbl=="YES")
    {
      fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_optimal_numbers_",target,sep=""),2400,2400)
      b=fviz_nbclust(a, clara, method = "silhouette",print.summary=T)
      b+theme_classic()
      dev.off()
      num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
    }
    
    #step 3. fit clara
    if(!exists("num.clus")) num.clus=n.clus
    clara.res <- clara(a, num.clus, samples = 50, pamLike = TRUE)
    
    #step 4. visualize CLARA clusters in data scattergram
    fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster_",target,sep=""),2400,2400)
    fviz_cluster(clara.res, 
                 palette = rainbow(num.clus), # color palette
                 ellipse.type = "t", # Concentration ellipse
                 geom = "point", pointsize = 1,
                 ggtheme = theme_classic())
    dev.off()
    
    #step 5. add cluster to input data
    dd.clara <- cbind(as.data.frame(a), cluster_clara = clara.res$cluster)
    dd.clara=dd.clara%>%rownames_to_column(var = "Same.return.SNo")%>%
      select(cluster_clara,Same.return.SNo)%>%remove_rownames()
    return(dd.clara)
  }
  #note: CLARA analysis as per Campbell et al 2017 on nfish as this has data at Sesssion level
  #       The CLARA (Clustering Large Applications) algorithm is an extension to the 
  #       PAM (Partitioning Around Medoids) clustering method for large data sets. It intended to 
  #       reduce the computation time in the case of large data set.
  Clus.vars=c("Catch.Target","Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar",
              "Catch.Groper","Catch.Snapper","Catch.Blue_mor")
  Tar.clus.vars=c("Catch.Whiskery","Catch.Gummy","Catch.Dusky","Catch.Sandbar")
  n.clus=c(2,2,2,2)  #from initial optimum number
  HndL.Species_targeting="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Species_targeting/"
  scalem="YES"
  
  Store.cluster=vector('list',length(SP.list)) 
  names(Store.cluster)=names(SP.list)
  for(i in 1:length(Tar.sp))
  {
    Store.cluster[[Tar.sp[i]]]=fn.cluster(data=DATA.list.LIVEWT.c.daily,TarSp=Tar.sp[i],target=Tar.clus.vars[i],
                                          varS=Clus.vars,scaling=scalem,check.clustrbl="NO",n.clus=n.clus[i])
  }
  
  #add cluster to original data for use in standardisations
  for(i in Tar.sp)
  {
    DATA.list.LIVEWT.c.daily[[i]]=left_join(DATA.list.LIVEWT.c.daily[[i]],Store.cluster[[i]],by=c("Same.return.SNo"))
  }
  rm(Store.cluster)
  
  #check cpue by cluster group
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster.boxplot",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,3,3,1),oma=c(2,1,.1,.5),las=1,mgp=c(2,.6,0))
  for(i in Tar.sp) boxplot((Catch.Target/Km.Gillnet.Hours.c)~cluster_clara,DATA.list.LIVEWT.c.daily[[i]],
                           main=names(DATA.list.LIVEWT.c.daily)[i],ylab="cpue")
  dev.off()
  
  
  
  fn.compare.targeting=function(DAT,Drop.var,Title)
  {
    d=DAT%>% select(c(cluster_clara,tolower(Clus.vars[-match(Drop.var,Clus.vars)])))
    d[-1]=d[-1]/rowSums(d[-1])
    d=d%>%gather('ID','value',-cluster_clara)           
    ggplot(d) +labs(title = Title,x = "",y="Proportion")+
      geom_boxplot(aes(x=ID, y=value, fill=cluster_clara))
  }
  for(s in 1:length(Tar.sp))
  {
    pdf(paste(HndL.Species_targeting,"Cluster/Targeting_",Nms.sp[Tar.sp[s]],".pdf",sep="")) 
    fn.compare.targeting(DAT=DATA.list.LIVEWT.c.daily[[Tar.sp[s]]],Drop.var=Tar.clus.vars[s],Nms.sp[Tar.sp[s]])
    dev.off()
  }
}


#4.15 Table of sensitivity scenarios       
Tab.Sensi=data.frame(Scenario=c("Base case","2 years","No efficiency"),
                     Vessels_used=c("5 years","2 years","5 years"),
                     Blocks_used=c("5 years","2 years","5 years"),
                     Efficiency_increase=c(rep("Yes",2),"No"))

setwd(paste(getwd(),"/Outputs/Paper",sep=""))
fn.word.table(WD=getwd(),TBL=Tab.Sensi,Doc.nm="Sensitivity tests",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")




#4.16 Compute foly and nominal index for exporting  

  #--Foly
#note: As done by Rory, the Foly (ie Effective) index has all records (good  & bad reporters) from all blocks 
#       and vessels within effective area without 0 catches

DATA.list.LIVEWT.c_all_reporters=vector('list',length(SP.list)) 
names(DATA.list.LIVEWT.c_all_reporters)=names(SP.list)
DATA.list.LIVEWT.c.daily_all_reporters=List.foly.nom=DATA.list.LIVEWT.c_all_reporters

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


#calculate foly
for ( i in Tar.sp)
{
  List.foly.nom[[i]]=export.foly(DATA.list.LIVEWT.c_all_reporters[[i]],DATA.list.LIVEWT.c.daily_all_reporters[[i]])
}


#-- Nominal 
#note: Ratio = mean(catch)/mean(effort)
#      Mean = mean(cpue)
#     LnMean= exp(mean(log(cpue))+bias corr)
#     DLnMean = exp(log(prop pos)+exp(mean(log(cpue))+bias corr)

Store_nom_cpues_monthly=vector('list',length(SP.list)) 
names(Store_nom_cpues_monthly)=names(SP.list)
Store_nom_cpues_daily=Store_nom_cpues_monthly
Hnd.ains="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Ainsline_different_cpues/"
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

#4.17 Evaluate balance of data subset based on QL  
#note: remove vessels with less than Min.Vess.yr (monthly) and Min.Vess.yr.d (daily)
#      for those vessels, keep blocks with more than Min.Vess.yr / Min.Vess.yr.d
if(!exists('BLKS.used.indi'))
{
  BLKS.used.indi=BLKS.used  #keep copy of indicators approach
  VES.used.indi=VES.used
  BLKS.used.daily.indi=BLKS.used.daily
  VES.used.daily.indi=VES.used.daily
}
hndl.kept="C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/"
HndL=paste(hndl.kept,'QL_balanced_design/',sep="")
for(s in Tar.sp)
{
  #monthly
  a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=names(SP.list)[s],
                      what="monthly",MN.YR=Min.Vess.yr,pLot=T)
  BLKS.used[[s]]=a$this.blks
  VES.used[[s]]=a$this.ves
  write.csv(BLKS.used[[s]],paste(hndl.kept,"blocks_used_",names(SP.list)[s],"_monthly.csv",sep=""))
  
  #daily
  a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=names(SP.list)[s],
                      what="daily",MN.YR=Min.Vess.yr.d,pLot=T)
  BLKS.used.daily[[s]]=a$this.blks
  VES.used.daily[[s]]=a$this.ves
  write.csv(BLKS.used.daily[[s]],paste(hndl.kept,"blocks_used_",names(SP.list)[s],"_daily.csv",sep=""))
}

#show blocks kept
CEX=.85
SRt=45
fn.fig(paste(hndl.kept,"block_used_map",sep=""),2400, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in Tar.sp)
{
  #Monthly
  fn.show.blk(dat=BLKS.used[[s]],CEX=CEX,SRt=SRt)
  if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  fn.show.blk(dat=BLKS.used.daily[[s]],CEX=CEX,SRt=SRt)
  if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",names(SP.list)[s],bty='n',cex=1.5)
}
mtext("Longitude",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Latitude",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")


#Show proportion of selected records by year
fn.fig("Appendix4_proportion_records_selected",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in Tar.sp)
{
  #Monthly
  with(subset(Store_nom_cpues_monthly[[s]]$Prop.selected,target==1),plot(yr,freq,ylim=c(0,1),
                          las=1,pch=19,cex=1.25,ylab="",xlab=""))
  if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)

  #Daily
  with(subset(Store_nom_cpues_daily[[s]]$Prop.selected,target==1),plot(yr,freq,ylim=c(0,1),
              las=1,pch=19,cex=1.25,ylab="",xlab=""))
  if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",names(SP.list)[s],bty='n',cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Proportion of records selected",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


#4.18 Show gummy monthly cpue effect of using km gn d or km g h
if(Model.run=="First")
{
  fn.fig("Appendix5_Gummy_km.gn.d_VS_km.gn.h",2400,2400)
  par(mfrow=c(4,2),mar=c(1,3,2,1),oma=c(2,.5,.5,1.75),mgp=c(1.5,.5,0),cex.lab=1.5)
  Get.Mns(d=DATA.list.LIVEWT.c$"Gummy Shark",grp="FINYEAR",
          Vars=c("HOURS.c","SHOTS.c",'Km.Gillnet.Days.c','Km.Gillnet.Hours.c','Km.Gillnet.Hours_shot.c','Catch.Target'),
          LGND=c("Hours fished per day","Number of shots",
                 'km gillnet days','km gillnet hours','km gillnet hours shot','Catch (kg)'),
          add.arrow=c(rep(TRUE,5),FALSE))
  dev.off()
}


#4.19 Export data to ainslie
Exprt_Ainslie="NO"
if(Exprt_Ainslie=="YES")
{
  for(s in Tar.sp)
  {
    write.csv(DATA.list.LIVEWT.c[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c[s]),"_monthly.csv",sep=""),row.names=F)
    write.csv(DATA.list.LIVEWT.c.daily[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c.daily[s]),"_daily.csv",sep=""),row.names=F)
  }
}


#4.20 Output data tables  
#Output table with number of records available in effective area and numbers used in standardisation
TABle=vector('list',length(SP.list))
names(TABle)=names(SP.list)
for(s in Tar.sp)
{
  
  #total records
  Tot.m=table(DATA.list.LIVEWT.c_all_reporters[[s]]$FINYEAR)
  Tot.d=table(DATA.list.LIVEWT.c.daily_all_reporters[[s]]$FINYEAR)
  
  Good.m=table(DATA.list.LIVEWT.c[[s]]$FINYEAR)
  Good.d=table(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR)
  
  
  #Used in Standardisation after selecting blocks, vessels and Qualification levels
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
                  B_Good.num.rec.eff.area=c(Good.m,Good.d),
                  C_Used.fo.stand=c(Stand.m,Stand.d))
  
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
                       'Number of reliable records','Number of records used in stand.'),
           HDR.span=c(1,1,N.species,N.species,N.species),
           HDR.2nd=c("","",rep(c("Dusky","Gummy","Sandbar","Whiskery"),3)))


#4.22 Construct index of abundance     
ZONES=c("West","Zone1","Zone2")
Eff.vars=c("km.gillnet.hours.c")
Covariates=c("temp.res","freo","freo_lag6","freo_lag12")
Predictors_monthly=c("finyear","vessel","month","blockx",Covariates) 
Predictors_daily=c("finyear","vessel","month","block10","shots.c","lunar","mean.depth",
                   "nlines.c","mesh",Covariates)
Response="catch.target"    #note that cpue is calculated inside stand function

Categorical=c("finyear","vessel","month","blockx","block10","shots.c",
              "lunar","mean.depth","nlines.c","mesh")


#   4.22.1 Explore data used for standardisation
#note: this applies cede() exploration, and checks for outliers in catch and effort
#      max monthly ktch, effort (~ 40 tonnes, ~ 5800 km gn h (@ 30 days X 24 h X 8000 m), respectively) 
#      max trip (daily kg) ktch, effort (~ 15 tonnes, ~ 1900 km gn h (@ 10 days X 24 h X 8000 m), respectively) 

  #check potential effect of predictors
if(do.Exploratory=="YES")
{
  hndl.expl="C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/"
  for(s in nnn)
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
    if(names(SP.list)[s]=="Whiskery Shark") fn.pred.effect(DATA=DATA.list.LIVEWT.c.daily[[s]],PREDS=Predictors_daily)
    dev.off()
  }
}

  #check data properties and degrees of freedom   
if(Model.run=="First")
{
  Prop.deg.free.m=Prop.deg.free.d=Prop.deg.free.d_n=Store_nom_cpues_monthly
  for(s in Tar.sp)
  {
    dummy=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used[[s]])
    Prop.m=properties(dummy)
    d.f=Prop.m[match(Predictors_monthly,rownames(Prop.m)),]
    d.f$dummy=as.numeric(with(d.f,ifelse(Class%in%c('numeric'),1,Unique)))
    d.f.m=data.frame(Deg.F=sum(d.f$dummy),Obser=nrow(dummy))
    Prop.deg.free.m[[s]]=list(dat=d.f,deg.f=d.f.m)
    
    dummy=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    Prop.d=properties(dummy)
    d.f=Prop.d[match(Predictors_daily,rownames(Prop.d)),]
    d.f$dummy=as.numeric(with(d.f,ifelse(Class%in%c('numeric'),1,Unique)))
    d.f.d=data.frame(Deg.F=sum(d.f$dummy),Obser=nrow(dummy))
    Prop.deg.free.d[[s]]=list(dat=d.f,deg.f=d.f.d)
  }
}

  #show records dropped by data selection process (starting from good records)
if(Model.run=="First")
{
  fn.fig("show_how_records_drop",2400, 2400) 
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,1),las=1,mgp=c(1.9,.7,0))
  for(s in Tar.sp)
  {
    #Monthly
    first=table(DATA.list.LIVEWT.c_all_reporters[[s]]$FINYEAR)
    dummy=Store_nom_cpues_monthly[[s]]$QL_dat
    second=table(dummy$finyear)
    dummy=subset(dummy,vessel%in%VES.used[[s]])
    third=table(dummy$finyear)
    dummy=subset(dummy,blockx%in%BLKS.used[[s]])
    fourth=table(dummy$finyear)
    LGTXT=NULL
    if(SP.list[[s]][1]==18007)LGTXT=c("whole data set","QL","QL_sel. vess.","QL_sel. vess. & block")
    barplot(rbind(first[match(names(second),names(first))],second,third,fourth),legend.text=LGTXT,
            args.legend=list(x="topright",cex=.9,bty='n',xjust=0))
    box()
    
    #Daily
    
    first=table(DATA.list.LIVEWT.c.daily_all_reporters[[s]]$FINYEAR)
    dummy=Store_nom_cpues_daily[[s]]$QL_dat
    second=table(dummy$finyear)
    dummy=subset(dummy,vessel%in%VES.used.daily[[s]])
    third=table(dummy$finyear)
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    fourth=table(dummy$finyear)
    barplot(rbind(first,second,third,fourth))
    box()
    mtext(names(SP.list)[s],4,0.5,cex=1,las=3)
  }
  mtext("Frequency",2,-0.75,las=3,outer=T,cex=1.5)
  dev.off()
}


#   4.22.2 Show applied effort creep
Fish.pow.inc=FINYEAR.ALL
names(Fish.pow.inc)=FINYEAR.ALL
Fish.pow.inc=as.numeric(substr(Fish.pow.inc,1,4))
Fish.pow.inc=c(seq(0,0.4,by=Fish.Pow),
              rep(0.4,length(1996:Fish.pow.inc[length(Fish.pow.inc)])))
if(Model.run=="First")
{
  fn.fig("Effort_creep/Effort_creep_applied",2400,2400)
  par(las=1,mgp=c(2.5,.9,0))
  plot(1:length(Fish.pow.inc),Fish.pow.inc*100,xaxt='n',ylab="Effort creep (% increase)",
       xlab="Finacial year",type='o',pch=19,cex.lab=1.6,cex.axis=1.35,cex=1.5)
  axis(1,1:length(Fish.pow.inc),F,tck=-0.01)
  axis(1,seq(1,length(Fish.pow.inc),10),FINYEAR.ALL[seq(1,length(Fish.pow.inc),10)],tck=-0.04,cex.axis=1.35)
  segments(6,Fish.pow.inc[6]*100,6,Fish.pow.inc[5]*100,col=2,lwd=2)
  segments(5,Fish.pow.inc[5]*100,6,Fish.pow.inc[5]*100,col=2,lwd=2)
  text(9,Fish.pow.inc[5]*100,"2 %",col=2,cex=2)
  dev.off()
}
Eff.creep=data.frame(finyear=FINYEAR.ALL,effort.creep=Fish.pow.inc)


#   4.22.3 Select model structure  

    #remove predictors identified as highly correlated
Predictors_monthly=Predictors_monthly[-match(c("temp.res",'freo','freo_lag6','freo_lag12'),Predictors_monthly)]
Predictors_daily=Predictors_daily[-match(c("temp.res",'nlines.c','freo','freo_lag6','freo_lag12'),Predictors_daily)]
cNSTNT=c('finyear','vessel','month','blockx')
cNSTNT.daily=c('finyear','vessel','month','block10')

      #4.22.3.1 extract best model
Best.Model=vector('list',length(SP.list)) 
names(Best.Model)=names(SP.list)
Best.Model.daily=Store.Best.Model=Store.Best.Model.daily=Best.Model.daily.gam=Best.Model
if(Def.mod.Str=="YES")     #takes 16 minutes
{
  fitFun= function(formula, data,always="", ...) 
  {
    glm(as.formula(paste(deparse(formula), always)), data=data, ...)
  }
  efrt="km.gillnet.hours.c"
  Inter="MainTerm"
  system.time({
  for(s in Tar.sp)
  {
      #monthly
    PREDS=Predictors_monthly
    d=Store_nom_cpues_monthly[[s]]$QL_dat%>%filter(vessel%in%VES.used[[s]] & blockx%in%BLKS.used[[s]])
    d=makecategorical(PREDS[which(PREDS%in%Categorical)],d)
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    always=cNSTNT
    Store.Best.Model[[s]]=fn.modl.sel(RESPNS="LNcpue")
    best.fm=Store.Best.Model[[s]]$BEST    
    best.fm=all.vars(best.fm)
    best.fm=as.formula(paste(best.fm[1],'~',paste(c(always,best.fm[-1]),collapse='+')))
    Best.Model[[s]]=best.fm
    rm(d,id.cov,PREDS)
    
      #daily
    d=Store_nom_cpues_daily[[s]]$QL_dat%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%BLKS.used.daily[[s]])
    PREDS=Predictors_daily
    d=d %>% mutate(mean.depth=10*round(mean.depth/10),
                   nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                   mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
     d=makecategorical(PREDS[which(PREDS%in%Categorical)],d)
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    always=cNSTNT.daily
    Store.Best.Model.daily[[s]]=fn.modl.sel(RESPNS="LNcpue")
    best.fm=Store.Best.Model.daily[[s]]$BEST    
    best.fm=all.vars(best.fm)
    best.fm=as.formula(paste(best.fm[1],'~',paste(c(always,best.fm[-1]),collapse='+')))
    Best.Model.daily[[s]]=best.fm
    rm(d,id.cov,PREDS)
}
  
  #4.22.3.3. show selection outcomes
  hndl.modl.sel="C:/Matias/Analyses/Catch and effort/Outputs/Model Selection/"
  for(s in Tar.sp)
  {
    pdf(paste(hndl.modl.sel,names(SP.list)[s],"_monthly.pdf",sep=""))
    fn.show.mod.sel(MODS=Store.Best.Model[[s]]$res,outs=length(Store.Best.Model[[s]]$res@objects))
    dev.off()
    
    pdf(paste(hndl.modl.sel,names(SP.list)[s],"_daily.pdf",sep=""))
    fn.show.mod.sel(MODS=Store.Best.Model.daily[[s]]$res,outs=length(Store.Best.Model.daily[[s]]$res@objects))
    dev.off()
  }
  
  })
}   
if(Def.mod.Str=="NO")
{
  for(s in nnn[-sort(Tar.sp)])
  {
    Best.Model[[s]]=formula(LNcpue ~ finyear + vessel + month + blockx)
    Best.Model.daily[[s]]=Best.Model[[s]]
    Best.Model.daily.gam[[s]]=formula(LNcpue ~ finyear + vessel + month + s(long10.corner,lat10.corner))
  }
  
  #Monthly
  Best.Model$`Gummy Shark`=formula(LNcpue ~ finyear + vessel + month + blockx )
  Best.Model$`Whiskery Shark`=  Best.Model$`Dusky Whaler Bronze Whaler`=
  Best.Model$`Sandbar Shark`=Best.Model$`Gummy Shark`
    
  #Daily
  Best.Model.daily$`Gummy Shark`=formula(LNcpue~finyear+vessel+month+block10+shots.c+lunar+mean.depth+mesh)
  Best.Model.daily$`Whiskery Shark`=formula(LNcpue~finyear+vessel+month+block10+shots.c+lunar+mesh)
  Best.Model.daily$`Dusky Whaler Bronze Whaler`=formula(LNcpue~finyear+vessel+month+block10+shots.c+lunar+mean.depth)
  Best.Model.daily$`Sandbar Shark`=formula(LNcpue ~finyear+vessel+month+block10+shots.c+lunar+mean.depth)
  
  
  Best.Model.daily.gam$`Gummy Shark`=formula(LNcpue~finyear+vessel+month+s(long10.corner,lat10.corner)+
                                                    shots.c+lunar+mean.depth+mesh)
  Best.Model.daily.gam$`Whiskery Shark`=formula(LNcpue~finyear+vessel+month+s(long10.corner,lat10.corner)+
                                                        shots.c+lunar+mesh)
  Best.Model.daily.gam$`Dusky Whaler Bronze Whaler`=formula(LNcpue~finyear+vessel+month+s(long10.corner,lat10.corner)+
                                                                   shots.c+lunar+mean.depth)
  Best.Model.daily.gam$`Sandbar Shark`=formula(LNcpue ~finyear+vessel+month+s(long10.corner,lat10.corner)+
                                                        shots.c+lunar+mean.depth)
}

#Export table of term levels
if(Model.run=="First")
{
  terms.table=vector('list',length(Tar.sp))
  for(s in Tar.sp)   
  {
    #monthly
    DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used[[s]])
    Mon=fn.table.terms(d=DAT,PREDS=Predictors_monthly)
    Mon=cbind(Data="Monthly",Mon)
    
    #daily
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    DAT=DAT %>% mutate(mean.depth=10*round(mean.depth/10),
                   nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                   mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    Day=fn.table.terms(d=DAT,PREDS=Predictors_daily)
    Day=cbind(Data="Daily",Day)
    DAT=rbind(Mon,Day)
    terms.table[[s]]=cbind(Species=names(SP.list)[s],DAT)
  }
  Tab.Terms=do.call(rbind,terms.table)
  row.names(Tab.Terms)=NULL
  fn.word.table(WD=getwd(),TBL=Tab.Terms,Doc.nm="Terms_table",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}


#   4.22.4 Run standardisation for Target species (based on qualification levels)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
  #monthly
system.time({Stand.out=foreach(s=Tar.sp,.packages=c('dplyr','cede')) %dopar%
  {
    DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])  #selet blocks and vessels
    DAT=subset(DAT,blockx%in%BLKS.used[[s]])
    DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select years with a min number of vessels
    return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                    efrt="km.gillnet.hours.c",Formula=Best.Model[[s]],Formula.gam=NULL))
    rm(DAT)
  }
})

  #daily
system.time({Stand.out.daily=foreach(s=Tar.sp,.packages=c('dplyr','cede','mgcv')) %dopar%
  {
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
    DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                      nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                      mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",
                    PREDS=Predictors_daily,efrt="km.gillnet.hours.c",
                    Formula=NULL,Formula.gam=Best.Model.daily.gam[[s]]))
    rm(DAT)
  }
})

names(Stand.out)=names(Stand.out.daily)=names(SP.list)[Tar.sp]
stopCluster(cl)

Nms.sp=names(SP.list)  
Nms.sp[match("Dusky Whaler Bronze Whaler",Nms.sp)]=c("Dusky Shark")
Nms.sp[match("All.Non.indicators",Nms.sp)]=c("All non-ind.")
Nms.sp[match("Hammerhead Sharks",Nms.sp)]=c("Hammerhead")
Nms.sp[match("Wobbegong",Nms.sp)]=c("Wobbegongs")
Nms.sp[match("Shortfin Mako",Nms.sp)]=c("Mako")


#   4.22.5 Run sensitivity tests     
    #free up some memory
rm(Species.list.daily,Species.list,Data.daily.GN,
   Data.monthly.GN,Effort.monthly,Effort.daily,DATA.list.LIVEWT.c.daily_all_reporters,
   DATA.list.LIVEWT.c_all_reporters)
if(Model.run=="First")      #takes 7 mins
{
  sens=Tab.Sensi
  sens$Efrt.used="km.gillnet.hours.c"
  Sens.pred=vector('list',length(SP.list))
  names(Sens.pred)=names(SP.list)
  Sens.pred.daily=Sens.pred.normlzd=Sens.pred.daily.normlzd=Sens.pred
  
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  system.time({for(s in Tar.sp)   #takes 55 secs
  {
    #1. Fit models
    sens_monthly=foreach(o=1:nrow(sens),.packages=c('dplyr','cede')) %dopar%   
      {
        MiN.YR=as.numeric(substr(sens$Vessels_used[o],1,1))
        EFrT=sens$Efrt.used[o]
        a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=names(SP.list)[s],
                            what="monthly",MN.YR=MiN.YR,pLot=F)
        if(names(SP.list)[s]=="Sandbar shark") a$this.blks=subset(a$this.blks,!a$this.blks=="3517") #cannot estimate this coef for 0==1
        DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%a$this.ves)
        DAT=subset(DAT,blockx%in%a$this.blks)
        DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
        return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,efrt=EFrT,
                        Formula=Best.Model[[s]],Formula.gam=NULL))
        rm(DAT)
      }
    sens_daily=foreach(o=1:nrow(sens),.packages=c('dplyr','cede')) %dopar%
      {
        MiN.YR=as.numeric(substr(sens$Vessels_used[o],1,1))
        EFrT=sens$Efrt.used[o]
        a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=names(SP.list)[s],
                            what="daily",MN.YR=MiN.YR,pLot=F)
        
        DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%a$this.ves)
        DAT=subset(DAT,blockx%in%a$this.blks)
        DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
        DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                          nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                          mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
        return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                        efrt=EFrT,Formula=Best.Model[[s]],Formula.gam=NULL))
        rm(DAT)
      }
    names(sens_monthly)=names(sens_daily)=sens$Scenario
    
    #2. Predict years based on emmeans (formerly lsmeans) considering log bias corr if required
    dummy=vector('list',length=nrow(sens))
    names(dummy)=sens$Scenario
    dummy.daily=dummy
    for(o in 1:nrow(sens))
    {
      d=sens_monthly[[o]]$DATA   #note: need data as global for ref_grid
      dummy[[o]]=pred.fun(MOD=sens_monthly[[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
      
      d=sens_daily[[o]]$DATA
      dummy.daily[[o]]=pred.fun(MOD=sens_daily[[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
      rm(d)
    }
    
    #3. Apply efficiency creep where required     
    for(o in 1:nrow(sens))
    {
      if(sens$Efficiency_increase[o]=="Yes")
      {
        #monthly
        add.crp=Eff.creep$effort.creep[match(dummy[[o]]$finyear,Eff.creep$finyear)]
        dummy[[o]]$response=dummy[[o]]$response*(1-add.crp)
        dummy[[o]]$lower.CL=dummy[[o]]$lower.CL*(1-add.crp)
        dummy[[o]]$upper.CL=dummy[[o]]$upper.CL*(1-add.crp)
        
        #daily
        add.crp=Eff.creep$effort.creep[match(dummy.daily[[o]]$finyear,Eff.creep$finyear)]
        dummy.daily[[o]]$response=dummy.daily[[o]]$response*(1-add.crp)
        dummy.daily[[o]]$lower.CL=dummy.daily[[o]]$lower.CL*(1-add.crp)
        dummy.daily[[o]]$upper.CL=dummy.daily[[o]]$upper.CL*(1-add.crp)
      }
    }
    
    #4. Normalise
    dummy.normlzd=dummy
    dummy.daily.normlzd=dummy.daily
    for(o in 1:nrow(sens))
    {
      #monthly
      Mn=mean(dummy[[o]]$response)
      dummy.normlzd[[o]]$response=dummy[[o]]$response/Mn
      dummy.normlzd[[o]]$lower.CL=dummy[[o]]$lower.CL/Mn
      dummy.normlzd[[o]]$upper.CL=dummy[[o]]$upper.CL/Mn
      
      #daily
      Mn=mean(dummy.daily[[o]]$response)
      dummy.daily.normlzd[[o]]$response=dummy.daily[[o]]$response/Mn
      dummy.daily.normlzd[[o]]$lower.CL=dummy.daily[[o]]$lower.CL/Mn
      dummy.daily.normlzd[[o]]$upper.CL=dummy.daily[[o]]$upper.CL/Mn
    }
    
    Sens.pred[[s]]=dummy
    Sens.pred.daily[[s]]=dummy.daily
    Sens.pred.normlzd[[s]]=dummy.normlzd
    Sens.pred.daily.normlzd[[s]]=dummy.daily.normlzd
  }})
  stopCluster(cl)
  
  #Plot   
  fn.fig("Appendix_Sensitivity",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
  for(s in Tar.sp)
  {
    LgND="NO"
    if(s==5)LgND="YES"
    Plot.cpue(cpuedata=Sens.pred[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
    if(s==5) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    LgND="NO"
    Plot.cpue(cpuedata=Sens.pred.daily[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
    if(s==5) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",Nms.sp[s],bty='n',cex=1.75)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  fn.fig("Appendix_Sensitivity_nomalised",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
  for(s in Tar.sp)
  {
    LgND="NO"
    if(s==5)LgND="YES"
    Plot.cpue(cpuedata=Sens.pred.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
    if(s==5) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    LgND="NO"
    Plot.cpue(cpuedata=Sens.pred.daily.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
    if(s==5) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",Nms.sp[s],bty='n',cex=1.75)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Relative CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
}
 
# 4.22.6 Run standardisation for Other species (based on delta method due to excess zeros)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
  #monthly          takes 2 sec
system.time({Stand.out.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','cede')) %dopar%
  {
    if(!is.null(DATA.list.LIVEWT.c[[s]]) & !is.null(BLKS.used[[s]]))
    {
      DAT=DATA.list.LIVEWT.c[[s]]
      colnames(DAT)=tolower(colnames(DAT)) 
      DAT=DAT%>%filter(vessel%in%VES.used[[s]] & blockx%in%as.numeric(BLKS.used[[s]]))  #selet blocks and vessels
      return(fn.delta(d=DAT,Response="catch.target",PREDS=Predictors_monthly,
                      efrt="km.gillnet.hours.c",Formula=Best.Model[[s]],Formula.gam=NULL))
      rm(DAT)
    }
  }
})
  #daily            takes 3.5 sec
system.time({Stand.out.daily.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','cede','mgcv')) %dopar%
  {
    if(!is.null(DATA.list.LIVEWT.c.daily[[s]])& !is.null(BLKS.used.daily[[s]]))
    {
      DAT=DATA.list.LIVEWT.c.daily[[s]]
      colnames(DAT)=tolower(colnames(DAT)) 
      DAT=DAT%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%as.numeric(BLKS.used.daily[[s]]))  #selet blocks and vessels
      return(fn.delta(d=DAT,Response="catch.target",PREDS=Predictors_monthly,
                      efrt="km.gillnet.hours.c",
                      Formula=NULL,Formula.gam=Best.Model.daily.gam[[s]]))
      rm(DAT)
    }
  }
})
names(Stand.out.other)=names(Stand.out.daily.other)=names(SP.list)[nnn[-sort(Tar.sp)]]
stopCluster(cl)


#Combine target and other species lists
Stand.out=c(Stand.out,Stand.out.other)
Stand.out=Stand.out[names(SP.list)]
Stand.out.daily=c(Stand.out.daily,Stand.out.daily.other)
Stand.out.daily=Stand.out.daily[names(SP.list)]

rm(Stand.out.other,Stand.out.daily.other)

#   4.22.7 Export deviance explained
if(Model.run=="First")  #takes 6 mins
{
  #calculate gam deviance explained by each term
  gam.list=Stand.out.daily
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  system.time({dummy1=foreach(s=nnn,.packages=c('mgcv')) %dopar%
    {
      if(!is.null(gam.list[[s]]))
      {
        if(s %in% Tar.sp)
        {
          ALLvars.gam=all.vars(Best.Model.daily.gam[[s]])[-1]
          ALLvars.gam=c(1,"s(long10.corner,lat10.corner)",ALLvars.gam[-match(c('long10.corner','lat10.corner'),ALLvars.gam)])
          dev.exp=rep(NA,length(ALLvars.gam))
          names(dev.exp)=ALLvars.gam
          for(g in 1:length(ALLvars.gam))
          {
            if(g==1) added.bit=paste(ALLvars.gam[g],collapse="+") else
              added.bit=paste(ALLvars.gam[1:g],collapse="+")
            Formula.gam=as.formula(paste('LNcpue',"~",added.bit))
            res.gam <-gam(Formula.gam,data=gam.list[[s]]$DATA,method="REML")
            dev.exp[g]=summary(res.gam)$dev.expl
          }
          dev.exp=dev.exp[-1]  
          return(list(dev.exp=dev.exp))
        }
        if(!s %in% Tar.sp)
        {
          ALLvars.gam=all.vars(Best.Model.daily.gam[[s]])[-1]
          ALLvars.gam=c(1,"s(long10.corner,lat10.corner)",ALLvars.gam[-match(c('long10.corner','lat10.corner'),ALLvars.gam)])
          dev.exp=rep(NA,length(ALLvars.gam))
          names(dev.exp)=ALLvars.gam
          dev.exp.BI=dev.exp
          for(g in 1:length(ALLvars.gam))
          {
            if(g==1) added.bit=paste(ALLvars.gam[g],collapse="+") else
              added.bit=paste(ALLvars.gam[1:g],collapse="+")
            
            Formula.bi.gam=as.formula(paste('catch.pos',"~",paste(added.bit,"offset(LNeffort)",sep="+")))
            res.gam_bi <-gam(Formula.bi.gam,data=gam.list[[s]]$DATA_bi, family="binomial",method="REML")
            
            Formula.gam=as.formula(paste('LNcpue',"~",added.bit))
            res.gam <-gam(Formula.gam,data=gam.list[[s]]$DATA,method="REML")
            
            dev.exp.BI[g]=summary(res.gam_bi)$dev.expl
            dev.exp[g]=summary(res.gam)$dev.expl
          }
          dev.exp.BI=dev.exp.BI[-1]
          dev.exp=dev.exp[-1] 
          return(list(Bi=dev.exp.BI,Pos=dev.exp))
        }
      }
    }
  })
  stopCluster(cl)
  names(dummy1)=names(gam.list)
  gam.list=dummy1
  
  Dev.exp=vector('list',length(SP.list)) 
  names(Dev.exp)=names(SP.list)
  Dev.exp.daily=Dev.exp
  system.time({for(s in Tar.sp)
  {
    if(!is.null(Stand.out[[s]]$res))Dev.exp[[s]]=Anova.and.Dev.exp(MOD=Stand.out[[s]]$res,
                                                SP=Nms.sp[s],type="Monthly",gam.extra=NULL)
    if(!is.null(Stand.out.daily[[s]]$res.gam))Dev.exp.daily[[s]]=Anova.and.Dev.exp(MOD=Stand.out.daily[[s]]$res.gam,
                                                SP=Nms.sp[s],type="Daily",gam.extra=gam.list[[s]]$dev.exp)
  }})   
  Tab.Dev.Exp=rbind(do.call(rbind,Dev.exp),do.call(rbind,Dev.exp.daily))
  rownames(Tab.Dev.Exp)=NULL
  fn.word.table(WD=getwd(),TBL=Tab.Dev.Exp,Doc.nm="ANOVA_table",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  Dev.exp.bi=Dev.exp
  Dev.exp.daily.bi=Dev.exp.daily
  system.time({for(s in nnn[-sort(Tar.sp)])
  {
    if(!is.null(Stand.out[[s]]$res))
    {
      Dev.exp[[s]]=Anova.and.Dev.exp(MOD=Stand.out[[s]]$res,
                        SP=Nms.sp[s],type="Monthly.pos",gam.extra=NULL)
      Dev.exp.bi[[s]]=Anova.and.Dev.exp(MOD=Stand.out[[s]]$res_bi,
                                     SP=Nms.sp[s],type="Monthly.bi",gam.extra=NULL)
    }

    if(!is.null(Stand.out.daily[[s]]$res.gam))
    {
      Dev.exp.daily[[s]]=Anova.and.Dev.exp(MOD=Stand.out.daily[[s]]$res.gam,
                                    SP=Nms.sp[s],type="Daily.pos",gam.extra=gam.list[[s]]$Pos)
      Dev.exp.daily.bi[[s]]=Anova.and.Dev.exp(MOD=Stand.out.daily[[s]]$res.gam,
                                           SP=Nms.sp[s],type="Daily.bi",gam.extra=gam.list[[s]]$Bi)
    }

  }})   
  Tab.Dev.Exp=rbind(do.call(rbind,Dev.exp),do.call(rbind,Dev.exp.bi),
                    do.call(rbind,Dev.exp.daily),do.call(rbind,Dev.exp.daily.bi))
  rownames(Tab.Dev.Exp)=NULL
  fn.word.table(WD=getwd(),TBL=Tab.Dev.Exp,Doc.nm="ANOVA_table_other.species",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}

#   4.22.8 Fit diagnostics
if(Model.run=="First") 
{
  fn.fig("Appendix 6",2000, 2400)
  par(mfcol=c(2*3,4),las=1,mar=c(2,2,1.75,1),oma=c(1,2,.1,2),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
  for(s in Tar.sp) 
  {
    #Monthly
    Pos.Diag.fn(MODEL=Stand.out[[s]]$res,SPECIES=Nms.sp[s],M=.9)
    
    #Daily
    Pos.Diag.fn(MODEL=Stand.out.daily[[s]]$res.gam,SPECIES="",M=.9)
  }
  mtext(c("Daily logbooks                                          Monthly returns     "),4,
        outer=T,las=3,line=0,cex=1.3)
  dev.off()
}

#   4.22.8 Plot base case and nominal 

#Extract comparable nominal cpue
Selected.data='CPUE.QL_target'
Selected.effort='km.gillnet.hours.c'

id.dat=match(Selected.data,names(Store_nom_cpues_monthly$`Whiskery Shark`))
id.eff=match(Selected.effort,names(Store_nom_cpues_monthly$`Whiskery Shark`[[1]]))

Nominl=vector('list',length(SP.list)) 
names(Nominl)=names(SP.list)
Nominl.daily=Nominl
for(s in nnn) 
{
  if(!is.null(Store_nom_cpues_monthly[[s]])) Nominl[[s]]=subset(Store_nom_cpues_monthly[[s]][[id.dat]][[id.eff]],method=="LnMean")
  if(!is.null(Store_nom_cpues_daily[[s]])) Nominl.daily[[s]]=subset(Store_nom_cpues_daily[[s]][[id.dat]][[id.eff]],method=="LnMean")
  if(is.na(match("finyear",names(Nominl[[s]]))))
  {
    Nominl[[s]]$finyear=Nominl[[s]]$season
    Nominl.daily[[s]]$finyear=Nominl.daily[[s]]$season
  }
}

#Compare glm and gam spatial predictions
if(compare.glm.gam=="YES")
{
  AOV.tab=function(Anova.tab,GLM)
  {
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/GLM$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    Dev.exp=sum(Term.dev.exp)
    
    ANOVA=as.data.frame.matrix(Anova.tab)
    Term=data.frame(Percent.dev.exp=Term.dev.exp)
    Table=ANOVA[-1,match(c("Deviance","Pr(>Chi)"),names(ANOVA))]
    Term=Term[match(rownames(Term), rownames(Table)),]
    Table=cbind(Table,Term)
    names(Table)[match("Term",names(Table))]="Percent.dev.exp"
    Table= Table %>% mutate_at(c("Deviance","Percent.dev.exp"), round, 3) 
    All=Table[1,]
    rownames(Table)=rownames(ANOVA)[2:nrow(ANOVA)]
    rownames(All)="Model"
    All[,1:ncol(All)]=NA
    All$Percent.dev.exp=round(Dev.exp,3)
    return(rbind(Table,All))
  }
  fn.compare.glm.gam=function(d)
  {
    Formula=formula(LNcpue ~ finyear + vessel + blockx )
    Formula_10=formula(LNcpue ~ finyear + vessel + block10 )
    Formula.gam=formula(LNcpue ~ finyear + vessel + s(long10.corner, lat10.corner) )
    
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    
    res <- glm(Formula,data=d)
    res_10 <- glm(Formula_10,data=d)
    res.gam <-gam(Formula.gam,data=d,method="REML")
    
    Tab.aov=AOV.tab(anova(res, test = "Chisq"),res)
    row.names(Tab.aov)=paste("glm",row.names(Tab.aov),sep='...')
    Tab.aov_10=AOV.tab(anova(res_10, test = "Chisq"),res_10)
    row.names(Tab.aov_10)=paste("glm_10",row.names(Tab.aov_10),sep='...')
    
    aa=summary(res.gam)
    Tab.aov.gam=data.frame(Deviance=NA,'Pr(>Chi)'=NA,Percent.dev.exp=round(100*aa$dev.expl,2))
    row.names(Tab.aov.gam)=paste("gam","Model",sep='...')
    colnames( Tab.aov.gam)=colnames(Tab.aov)
    
    grid.table(rbind(Tab.aov,Tab.aov_10,Tab.aov.gam))
    
    
    #finyear
    BLKr=sort(table(d$blockx))
    BLKr=names(BLKr[length(BLKr)])
    BLKr10=sort(table(d$block10))
    BLKr10=names(BLKr10[length(BLKr10)])
    long10.corner=as.numeric(d[d$block10==BLKr10,match('long10.corner',names(d))][1])
    lat10.corner=as.numeric(d[d$block10==BLKr10,match('lat10.corner',names(d))][1])
    
    VSl= sort(table(d$vessel))
    VSl=names(VSl[length(VSl)])
    
    new.block=data.frame(finyear=factor(levels(d$finyear)),
                         vessel=factor(VSl,levels(d$vessel)),
                         blockx=factor(BLKr,levels(d$blockx)))
    new.block$cpue=predict(res,newdata=new.block,type='response')
    
    new.block10=data.frame(finyear=factor(levels(d$finyear)),
                           vessel=factor(VSl,levels(d$vessel)),
                           block10=factor(BLKr10,levels(d$block10)))
    new.block10$cpue=predict(res_10,newdata=new.block10,type='response')
    
    new.gam=data.frame(finyear=factor(levels(d$finyear)),
                       vessel=factor(VSl,levels(d$vessel)),
                       long10.corner=long10.corner,
                       lat10.corner=lat10.corner)
    new.gam$cpue=predict(res.gam,newdata=new.gam,type='response')
    
    new.block$cpue=new.block$cpue/mean(new.block$cpue)
    new.block10$cpue=new.block10$cpue/mean(new.block10$cpue)
    new.gam$cpue=new.gam$cpue/mean(new.gam$cpue)
    
    yrs=new.gam$finyear
    plot(1:length(yrs),new.block$cpue,xaxt='n',xlab="financial year",ylab="normalised cpue",
         cex.lab=1.5,pch=19,cex=1.5,type='o',
         ylim=c(min(c(new.block$cpue,new.block10$cpue,new.gam$cpue)),max(c(new.block$cpue,new.block10$cpue,new.gam$cpue))))
    axis(1,1:length(yrs),yrs)
    points((1:length(yrs))+.15,new.block10$cpue,pch=19,col='cyan3',cex=1.5)
    lines((1:length(yrs))+.15,new.block10$cpue,pch=19,col='cyan3')
    points((1:length(yrs))-.15,new.gam$cpue,pch=19,col='grey60',cex=1.5)
    lines((1:length(yrs))-.15,new.gam$cpue,pch=19,col='grey60')
    legend("topright",c("block","block10","gam"),bty='n',pch=19,col=c("black","cyan3","grey60"),cex=1.5)
    #space
    par(mfcol=c(3,1),mar=c(2,2,2,2),oma=c(1,1,.1,.1),mgp=c(1.5,.6,0))
    
    FINYr=sort(table(d$finyear))
    FINYr=names(FINYr[length(FINYr)])
    VSl= sort(table(d$vessel))
    VSl=names(VSl[length(VSl)])
    
    new.block=data.frame(finyear=factor(FINYr,levels(d$finyear)),
                         vessel=factor(VSl,levels(d$vessel)),
                         blockx=factor(levels(d$blockx)))
    new.block$cpue=predict(res,newdata=new.block,type='response')
    new.block=new.block%>%mutate( Lat=-as.numeric(substr(blockx,1,2)),
                                  Long=100+as.numeric(substr(blockx,3,4)))%>%
      select(-c(blockx,finyear,vessel))
    YLIM=c(min(new.block$Lat),max(new.block$Lat))
    XLIM=c(min(new.block$Long),max(new.block$Long))
    seq.lat=seq(YLIM[1],YLIM[2])
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.block$Lat))]
    seq.lon=seq(XLIM[1],XLIM[2])
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.block$Long))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(Lat=seq.lat, Long=seq.lon)
      new.block=combo%>%left_join(new.block,by=c("Lat","Long"))
    }
    new.block=new.block%>%spread(Lat,cpue)
    Lon=as.numeric(new.block$Long)
    new.block=as.matrix(new.block[,-1]) 
    LaT=as.numeric(colnames(new.block))
    brk<- quantile( c(new.block),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.block, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Blockx",bty='n',cex=1.5)
    
    
    
    new.block10=data.frame(finyear=factor(FINYr,levels(d$finyear)),
                           vessel=factor(VSl,levels(d$vessel)),
                           block10=factor(levels(d$block10)))
    new.block10$cpue=predict(res_10,newdata=new.block10,type='response')
    dummy=subset(d,select=c(block10,lat10.corner,long10.corner)) %>%
      mutate(block10=as.character(block10)) %>%
      distinct(block10,.keep_all =T)
    new.block10=new.block10%>%mutate(block10=as.character(block10))%>%
      left_join(dummy,by="block10") %>%
      select(-c(block10,finyear,vessel)) %>%
      mutate(lat10.corner=round(lat10.corner,2),
             long10.corner=round(long10.corner,2))
    YLIM=c(min(new.block10$lat10.corner),max(new.block10$lat10.corner))
    XLIM=c(min(new.block10$long10.corner),max(new.block10$long10.corner))
    #seq(112,129,length.out = 7+6*16)
    seq.lat=c(-26.83, -26.67, -26.50, -26.33, -26.17, -26.00,
              -27.83, -27.67, -27.50, -27.33, -27.17, -27.00,
              -28.83, -28.67, -28.50, -28.33, -28.17, -28.00,
              -29.83, -29.67, -29.50, -29.33, -29.17, -29.00,
              -30.83, -30.67, -30.50, -30.33, -30.17, -30.00,
              -31.83, -31.67, -31.50, -31.33, -31.17, -31.00,
              -32.83, -32.67, -32.50, -32.33, -32.17, -32.00,
              -33.83, -33.67, -33.50, -33.33, -33.17, -33.00,
              -34.83, -34.67, -34.50, -34.33, -34.17, -34.00,
              -35.83, -35.67, -35.50, -35.33, -35.17, -35.00)
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.block10$lat10.corner))]
    seq.lon=c(113.83, 113.67, 113.50, 113.33, 113.17, 113.00,
              114.83, 114.67, 114.50, 114.33, 114.17, 114.00,
              115.83, 115.67, 115.50, 115.33, 115.17, 115.00,
              116.83, 116.67, 116.50, 116.33, 116.17, 116.00,
              117.83, 117.67, 117.50, 117.33, 117.17, 117.00,
              118.83, 118.67, 118.50, 118.33, 118.17, 118.00,
              119.83, 119.67, 119.50, 119.33, 119.17, 119.00,
              120.83, 120.67, 120.50, 120.33, 120.17, 120.00,
              121.83, 121.67, 121.50, 121.33, 121.17, 121.00,
              122.83, 122.67, 122.50, 122.33, 122.17, 122.00,
              123.83, 123.67, 123.50, 123.33, 123.17, 123.00,
              124.83, 124.67, 124.50, 124.33, 124.17, 124.00,
              125.83, 125.67, 125.50, 125.33, 125.17, 125.00,
              126.83, 126.67, 126.50, 126.33, 126.17, 126.00,
              127.83, 127.67, 127.50, 127.33, 127.17, 127.00,
              128.83, 128.67, 128.50, 128.33, 128.17, 128.00,
              129.83, 129.67, 129.50, 129.33, 129.17, 129.00)
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.block10$long10.corner))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(lat10.corner=seq.lat, long10.corner=seq.lon)
      new.block10=combo%>%left_join(new.block10,by=c("lat10.corner","long10.corner"))
    }
    new.block10=new.block10%>%spread(lat10.corner,cpue)
    Lon=as.numeric(new.block10$long10.corner)
    new.block10=as.matrix(new.block10[,-1]) 
    LaT=as.numeric(colnames(new.block10))
    brk<- quantile( c(new.block10),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.block10, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Block10",bty='n',cex=1.5)
    
    
    
    ##option using lsmeans
    # dummy=subset(d,select=c(block10,lat10.corner,long10.corner)) %>%
    #   mutate(block10=as.character(block10)) %>%
    #   distinct(block10,.keep_all =T)%>%
    #   arrange(block10)
    # Long.seq=dummy$long10.corner
    # Lat.seq=dummy$lat10.corner
    # a=summary(ref_grid(res.gam, at = list(finyear=new.finyr,lat10.corner =Lat.seq ,long10.corner = Long.seq)))
    
    new.gam=d%>%select(block10,lat10.corner,long10.corner)%>%
      distinct(block10,.keep_all=T)%>%
      mutate(finyear=factor(FINYr,levels(d$finyear)),
             vessel=factor(VSl,levels(d$vessel)))
    new.gam$cpue=predict(res.gam,newdata=new.gam,type='response')
    new.gam=new.gam%>%select(-c(block10,finyear,vessel))%>%
      mutate(lat10.corner=round(lat10.corner,2),
             long10.corner=round(long10.corner,2))
    
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.gam$lat10.corner))]
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.gam$long10.corner))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(lat10.corner=seq.lat, long10.corner=seq.lon)
      new.gam=combo%>%left_join(new.gam,by=c("lat10.corner","long10.corner"))
    }
    new.gam=new.gam%>%spread(lat10.corner,cpue)
    Lon=as.numeric(new.gam$long10.corner)
    new.gam=as.matrix(new.gam[,-1]) 
    LaT=as.numeric(colnames(new.gam))
    brk<- quantile( c(new.gam),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.gam, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Gam",bty='n',cex=1.5)
    
    mtext("Longitude",1,outer=T)
    mtext("Latitude",2,outer=T)
    
    
    #explore gam
    xlpr.gam="NO"
    if(xlpr.gam=="YES")
    {
      par(mfrow = c(2,2))
      gam.check(res.gam)
      #small p-values indicate that residuals are not randomly distributed. 
      # This often means there are not enough basis functions
      
      
      par(mfrow = c(1,2))
      plot(res.gam, residuals = TRUE, pch = 1)
      plot(res.gam, residuals = TRUE, pch = 1, scheme = 2)
      #In this plot the axes represent values of our predictor variables, x1 and x2. 
      #The interior is a topographic map of predicted values. 
      #The contour lines represent points of equal predicted values, and they are labeled. 
      #The dotted lines show uncertainty in prediction; they represent how contour
      # lines would move if predictions were one standard error higher or lower.
      
      vis.gam(x = res.gam,                # GAM object
              view = c("long10.corner", "lat10.corner"),   # variables
              plot.type = "contour", too.far = 0.05)    # kind of plot 
      
    }
    
    
  }
  system.time({for(s in Tar.sp)
  {
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
    DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                      nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                      mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    
    Response="catch.target"
    RESPNS="LNcpue"
    PREDS=c("finyear","vessel","block10","blockx")
    efrt="km.gillnet.hours.c"
    
    pdf(paste(getwd(),"/Compare glm and gam/Compare.glm.vs.gam_",Nms.sp[s],".pdf",sep=""))
    fn.compare.glm.gam(d=DAT)
    dev.off()
    rm(DAT)
  }
  })
}


#Predict years (considering log bias corr if required)

  #Other species
Niter=1000   #MC interations
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
      #monthly          takes 0.3 sec per iteration
system.time({Pred.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mvtnorm')) %dopar%
  {
    if(!is.null(Stand.out[[s]]))
    {
      return(fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$res_bi,
                              MOD=Stand.out[[s]]$res,
                              BiData=Stand.out[[s]]$DATA_bi,
                              PosData=Stand.out[[s]]$DATA,
                              niter=Niter,
                              pred.term='finyear',
                              ALL.terms=Predictors_monthly))
    } 
  }
})
      #Daily              takes 0.2 sec per iteration
system.time({Pred.daily.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mvtnorm')) %do%
  {
    if(!is.null(Stand.out.daily[[s]]))
    {
      return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$res.gam_bi,
                              MOD=Stand.out.daily[[s]]$res.gam,
                              BiData=Stand.out.daily[[s]]$DATA_bi,
                              PosData=Stand.out.daily[[s]]$DATA,
                              niter=Niter,
                              pred.term='finyear',
                              ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
    } 
  }
})
names(Pred.other)=names(Pred.daily.other)=names(SP.list)[nnn[-sort(Tar.sp)]]

  #Target species
    #monthly          takes 80 sec
system.time({Pred.tar=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
  {
    d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
    return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="finyear",Pred.type="link"))
    rm(d)
  }
})
    #Daily              takes 70 sec
system.time({Pred.daily.tar=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
  {
    d=Stand.out.daily[[s]]$DATA
    return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="finyear",Pred.type="link"))
    rm(d)
  }
})
names(Pred.tar)=names(Pred.daily.tar)=names(SP.list)[Tar.sp]
stopCluster(cl) 

Pred=c(Pred.tar,Pred.other)
Pred=Pred[names(SP.list)]
Pred.daily=c(Pred.daily.tar,Pred.daily.other)
Pred.daily=Pred.daily[names(SP.list)]

  #Apply efficiency creep      
Pred.creep=Pred
Pred.daily.creep=Pred.daily
Nominl.creep=Nominl
Nominl.daily.creep=Nominl.daily
for(s in nnn)
{
  #monthly
  if(!is.null(Pred.creep[[s]]))
  {
    yrs=Nominl.creep[[s]]$finyear
    yrs=paste(yrs,substr(yrs+1,3,4),sep="-")
    Pred.creep[[s]]=subset(Pred.creep[[s]],finyear%in%yrs)
    
    add.crp=Eff.creep$effort.creep[match(Pred.creep[[s]]$finyear,Eff.creep$finyear)]
    Pred.creep[[s]]$response=Pred.creep[[s]]$response*(1-add.crp)
    Pred.creep[[s]]$lower.CL=Pred.creep[[s]]$lower.CL*(1-add.crp)
    Pred.creep[[s]]$upper.CL=Pred.creep[[s]]$upper.CL*(1-add.crp)
    
    yrs=as.character(substr(Pred.creep[[s]]$finyear,1,4))
    Nominl.creep[[s]]=subset(Nominl.creep[[s]],finyear%in%yrs)
    Nominl.creep[[s]]$response=Nominl.creep[[s]]$mean*(1-add.crp)
    Nominl.creep[[s]]$lower.CL=Nominl.creep[[s]]$lowCL*(1-add.crp)
    Nominl.creep[[s]]$upper.CL=Nominl.creep[[s]]$uppCL*(1-add.crp)
  }
  
  #daily
  if(!is.null(Pred.daily.creep[[s]]))
  {
    yrs=Nominl.daily.creep[[s]]$finyear
    yrs=paste(yrs,substr(yrs+1,3,4),sep="-")
    Pred.daily.creep[[s]]=subset(Pred.daily.creep[[s]],finyear%in%yrs)
    add.crp=Eff.creep$effort.creep[match(Pred.daily.creep[[s]]$finyear,Eff.creep$finyear)]
    Pred.daily.creep[[s]]$response=Pred.daily.creep[[s]]$response*(1-add.crp)
    Pred.daily.creep[[s]]$lower.CL=Pred.daily.creep[[s]]$lower.CL*(1-add.crp)
    Pred.daily.creep[[s]]$upper.CL=Pred.daily.creep[[s]]$upper.CL*(1-add.crp)
    
    yrs=as.character(substr(Pred.daily.creep[[s]]$finyear,1,4))
    Nominl.daily.creep[[s]]=subset(Nominl.daily.creep[[s]],finyear%in%yrs)
    Nominl.daily.creep[[s]]$response=Nominl.daily.creep[[s]]$mean*(1-add.crp)
    Nominl.daily.creep[[s]]$lower.CL=Nominl.daily.creep[[s]]$lowCL*(1-add.crp)
    Nominl.daily.creep[[s]]$upper.CL=Nominl.daily.creep[[s]]$uppCL*(1-add.crp)
  }
}

      #Plot Target
fn.fig("Figure 4.Annual_Index",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
for(s in Tar.sp)
{
  #Monthly
  Mon.dat=list(Standardised=Pred.creep[[s]],Nominal=Nominl.creep[[s]])
  LgND="NO"
  if(s==Tar.sp[1])LgND="YES"
  Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
  if(s==Tar.sp[1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  Daily.dat=list(Standardised=Pred.daily.creep[[s]],Nominal=Nominl.daily.creep[[s]])
  LgND="NO"
  Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
  if(s==Tar.sp[1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext(Nms.sp[s],4,line=1,las=3,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("CPUE (kg/ km gillnet hour)",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
dev.off()


      #Plot Target normalised    
Pred.normlzd=Pred.creep
Pred.daily.normlzd=Pred.daily.creep
Nominl.normlzd=Nominl.creep
Nominl.daily.normlzd=Nominl.daily.creep
for(s in nnn)
{
  #monthly
  if(!is.null(Pred.normlzd[[s]]))
  {
    Mn=mean(Pred.normlzd[[s]]$response)
    Pred.normlzd[[s]]$response=Pred.normlzd[[s]]$response/Mn
    Pred.normlzd[[s]]$lower.CL=Pred.normlzd[[s]]$lower.CL/Mn
    Pred.normlzd[[s]]$upper.CL=Pred.normlzd[[s]]$upper.CL/Mn
    
    Mn=mean(Nominl.normlzd[[s]]$response)
    Nominl.normlzd[[s]]$response=Nominl.normlzd[[s]]$response/Mn
    Nominl.normlzd[[s]]$lower.CL=Nominl.normlzd[[s]]$lower.CL/Mn
    Nominl.normlzd[[s]]$upper.CL=Nominl.normlzd[[s]]$upper.CL/Mn
  }
  
  #daily
  if(!is.null(Pred.daily.normlzd[[s]]))
  {
    Mn=mean(Pred.daily.normlzd[[s]]$response)
    Pred.daily.normlzd[[s]]$response=Pred.daily.normlzd[[s]]$response/Mn
    Pred.daily.normlzd[[s]]$lower.CL=Pred.daily.normlzd[[s]]$lower.CL/Mn
    Pred.daily.normlzd[[s]]$upper.CL=Pred.daily.normlzd[[s]]$upper.CL/Mn
    
    Mn=mean(Nominl.daily.normlzd[[s]]$response)
    Nominl.daily.normlzd[[s]]$response=Nominl.daily.normlzd[[s]]$response/Mn
    Nominl.daily.normlzd[[s]]$lower.CL=Nominl.daily.normlzd[[s]]$lower.CL/Mn
    Nominl.daily.normlzd[[s]]$upper.CL=Nominl.daily.normlzd[[s]]$upper.CL/Mn
  }
}
fn.fig("Figure 4.Annual_Index_normalised",2000, 2400)    
#no point in showing nominal as it is not comparable in relative terms
par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
for(s in Tar.sp)
{
  #Monthly
  Mon.dat=list(Standardised=Pred.normlzd[[s]])
  LgND="NO"
  if(s==Tar.sp[1])LgND="YES"
  Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
  if(s==Tar.sp[1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  Daily.dat=list(Standardised=Pred.daily.normlzd[[s]])
  LgND="NO"
  Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
  if(s==Tar.sp[1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext(Nms.sp[s],4,line=1,las=3,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
dev.off()


#Plot Other species normalised
fn.fig("Figure 4.Annual_Index_normalised_other species",1200, 2400)    
par(mfrow=c(length(nnn[-sort(Tar.sp)]),2),mar=c(1,1,.75,.95),oma=c(2.5,3,1,.25),las=1,mgp=c(1.9,.5,0),cex.axis=.8)
for(s in nnn[-sort(Tar.sp)])
{
  #Monthly
  Mon.dat=list(Standardised=Pred.normlzd[[s]])
  LgND="NO"
  suppressWarnings(Plot.cpue.other(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1,Yvar="finyear",All.yrs=FINYEAR.monthly))
  if(s==nnn[-sort(Tar.sp)][1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.25)
  
  #Daily
  Daily.dat=list(Standardised=Pred.daily.normlzd[[s]])
  LgND="NO"
  suppressWarnings(Plot.cpue.other(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1,Yvar="finyear",All.yrs=FINYEAR.daily))
  if(s==nnn[-sort(Tar.sp)][1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.25)
  mtext(gsub(paste("Shark", collapse="|"), "", Nms.sp[s]),4,0.1,cex=.7,las=3)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.25,outer=T)
mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.25,outer=T)
dev.off()

#boxplots of factors
if(Model.run=="First")
{
  for(s in Tar.sp)
  {
    fn.fig(paste("boxplot.daily.preds_",Nms.sp[s],sep=""),2000, 2400)  
    par(mfcol=c(3,2))  
    with(Stand.out.daily[[s]]$DATA,
         {
           Ylim=c(0,quantile(catch.target/km.gillnet.hours.c,probs=0.99))
           boxplot((catch.target/km.gillnet.hours.c)~cluster_clara,col="grey70",ylim=Ylim)
           boxplot((catch.target/km.gillnet.hours.c)~mean.depth,col="grey70",ylim=Ylim)
           boxplot((catch.target/km.gillnet.hours.c)~lunar,col="grey70",ylim=Ylim)
           boxplot((catch.target/km.gillnet.hours.c)~month,col="grey70",ylim=Ylim)
           boxplot((catch.target/km.gillnet.hours.c)~mesh,col="grey70",ylim=Ylim)
           boxplot((catch.target/km.gillnet.hours.c)~shots.c,col="grey70",ylim=Ylim)
         })
    dev.off()
  }
}


#   4.22.9 Show effects for other terms 
if(Model.run=="First")
{
  #Predict vessel based on emmeans (formerly lsmeans) considering log bias corr if required
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  #monthly          takes 80 sec
  system.time({Pred.vess=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
      return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="vessel",Pred.type="link"))
      rm(d)
    }
  })
  #Daily              takes 70 sec
  system.time({Pred.vess.daily=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out.daily[[s]]$DATA
      return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="vessel",Pred.type="link"))
      rm(d)
    }
  })
  names(Pred.vess)=names(Pred.vess.daily)=names(SP.list)[Tar.sp]
  stopCluster(cl)
  fn.fig("Figure 2.Vessel effect",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,2.5,.1,.2),las=1,mgp=c(1.9,.7,0))
  for(s in 1:length(Pred.vess))
  {
    #Monthly
    Mon.dat=Pred.vess[[s]]
    LgND="NO"
    Mon.dat$vessel=1:nrow(Mon.dat)
    if(s==1)LgND="YES"
    Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="vessel",add.lines="NO")
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    #Daily
    Daily.dat=Pred.vess.daily[[s]]
    LgND="NO"
    Daily.dat$vessel=1:nrow(Daily.dat)
    Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="vessel",add.lines="NO")
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    mtext( Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
  }
  mtext("Vessel",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CPUE (kg/ km gillnet hour)",side=2,line=0.5,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  
  #Predict month based on emmeans (formerly lsmeans) considering log bias corr if required
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
    #monthly          takes 80 sec
  system.time({Pred.month=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
      return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="month",Pred.type="link"))
      rm(d)
    }
  })
    #Daily              takes 70 sec
  system.time({Pred.month.daily=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out.daily[[s]]$DATA
      return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="month",Pred.type="link"))
      rm(d)
    }
  })
  names(Pred.month)=names(Pred.month.daily)=names(SP.list)[Tar.sp]
  stopCluster(cl)
  
  fn.fig("Figure.Monthly effect",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,2.5,.1,.2),las=1,mgp=c(1.9,.7,0))
  for(s in 1:length(Pred.vess))
  {
    #Monthly
    Mon.dat=Pred.month[[s]]
    LgND="NO"
    if(s==1)LgND="YES"
    Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month",add.lines="YES")
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    #Daily
    Daily.dat=Pred.month.daily[[s]]
    LgND="NO"
    Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month",add.lines="YES")
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    mtext( Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
  }
  mtext("Month",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CPUE (kg/ km gillnet hour)",side=2,line=0,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  
  #Predict spatial cpue (monthly blocks and daily lat/long)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
    #monthly          takes 4 sec
  system.time({Pred.spatial.monthly=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %dopar%
    {
      
      return(pred.fun.spatial(DAT=Stand.out[[s]]$DATA,MOD=Stand.out[[s]]$res,
                              PRED='blockx',FORM=Best.Model[[s]],
                              Spatial.grid=data.frame(blockx=factor(levels(Stand.out[[s]]$DATA$blockx))))
      )
      
    }
  })
    #Daily              takes 0.1 sec
  system.time({Pred.spatial.daily=foreach(s=Tar.sp,.packages=c('dplyr')) %do%
    {
      return(pred.fun.spatial(DAT=Stand.out.daily[[s]]$DATA,MOD=Stand.out.daily[[s]]$res.gam,
                              PRED=c('long10.corner','lat10.corner'),FORM=Best.Model.daily.gam[[s]],
                              Spatial.grid=Stand.out.daily[[s]]$DATA%>%select(block10,lat10.corner,long10.corner)%>%
                                distinct(block10,.keep_all=T)))
      
    }
  })
  names(Pred.spatial.monthly)=names(Pred.spatial.daily)=names(SP.list)[Tar.sp]
  stopCluster(cl)
  fn.fig("Figure 3.Spatial effect",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,2,1.5,3),oma=c(2.5,2.5,.1,2.5),las=1,mgp=c(1.9,.5,0))
  for(s in 1:length(Pred.vess))
  {
    #Monthly
    Full.long=seq(113,129)
    Full.lat=seq(-35,-26) 
    Plot.cpue.spatial(cpuedata=Pred.spatial.monthly[[s]],var='blockx')
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.35)
    #if(s%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    #axis(2,seq(-36,-26,2),rev(seq(26,36,2)),tck=-0.025,cex.axis=1.35)
    
    
    #Daily
    Full.long=apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(113:129)),rep(113:129,each=6)),1,sum)
    Full.lat=-apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(35:26)),rep(35:26,each=6)),1,sum)
    Plot.cpue.spatial(cpuedata=Pred.spatial.daily[[s]],var=c('long10.corner','lat10.corner'))
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.35)
    #if(s%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    mtext( Nms.sp[Tar.sp[s]],4,line=6,las=3,cex=1.5)
  }
  mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.2,font=1,las=0,cex=1.35,outer=T)
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=0.65,font=1,las=0,cex=1.35,outer=T)
  dev.off()
}


#   4.22.10 Influence plots     
if(do.influence=="YES")
{
  HnDll="C:/Matias/Analyses/Catch and effort/Outputs/Influence.plot/"
  Store.Influence=vector('list',length(SP.list)) 
  names(Store.Influence)=names(SP.list)
  Store.Influence.daily=Store.Influence
  for(s in Tar.sp)
  {
    #Monthly
    Terms=all.vars(Best.Model[[s]])[-1]
    Terms=Terms[-match("finyear",Terms)]
    Term.Type=Terms
    Term.Type=ifelse(Term.Type%in%Categorical,"CAT",NA)
    names(Term.Type)=Terms
    pdf(paste(HnDll,Nms.sp[s],".monthly.CDI.pdf",sep=""))
    Store.Influence[[s]]=Influence.fn(MOD=Stand.out[[s]]$res,DAT=Stand.out[[s]]$DATA,
                  Term.type=Term.Type,termS=Terms,add.Influence="YES",SCALER=4)
    dev.off()
    
    
    #Daily
    Terms=all.vars(Best.Model.daily.gam[[s]])[-1]
    Terms=Terms[-match("finyear",Terms)]
    Term.Type=Terms
    Term.Type=ifelse(Term.Type%in%Categorical,"CAT",NA)
    names(Term.Type)=Terms
    Term.Type=Term.Type[!is.na(Term.Type)]

    pdf(paste(HnDll,Nms.sp[s],".daily.CDI.pdf",sep=""))
    Store.Influence.daily[[s]]=Influence.fn(MOD=Stand.out.daily[[s]]$res.gam,
            DAT=Stand.out.daily[[s]]$DATA,Term.type=Term.Type,termS=Terms,
            add.Influence="YES",SCALER=4)
    dev.off()
  } 
  Store.Influence=Store.Influence[!sapply(Store.Influence, is.null)] 
  Store.Influence.daily=Store.Influence.daily[!sapply(Store.Influence.daily, is.null)] 
  
  Over.inf.monthly=do.call(rbind,lapply(Store.Influence, '[[', match('Over.all.influence',names(Store.Influence[[1]]))))
  #ACA. keep updating 1:N.species with Tar.sp where appropriate, and SPECIES.vec for  Nms.sp 
  Over.inf.daily=do.call(rbind,lapply(Store.Influence.daily, '[[', match('Over.all.influence',names(Store.Influence.daily[[1]]))))
  Over.inf.monthly=round(100*Over.inf.monthly,2)
  Over.inf.daily=round(100*Over.inf.daily,2)
  Over.inf.percent=rbind(data.frame(Species=rownames(Over.inf.monthly),dat="Monthly",Over.inf.monthly),
                         data.frame(Species=rownames(Over.inf.daily),dat="Daily",Over.inf.daily))
  
  fn.word.table(WD=getwd(),TBL=Over.inf.percent,Doc.nm="Influence_table_percent",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  #Compare influence of all terms
  LWD=3
  LTY.col=c("black","grey70","grey25","grey45","grey55")
  Whr=c("bottomright","bottomright","bottomright","topright")
  Whr.d=c("bottomright","bottomright","bottomright","bottom")
  
  
  fn.fig("Figure 5.All.terms.Influence",2400, 2400)
  par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
  for ( s in 1:N.species)
  {
    #Monthly
    dumy=unlist(lapply(Store.Influence, "[[", match("Annual.Dev",names(Store.Influence[[1]]))))
    YLIM=c(min(dumy),max(dumy))
    Compare.term.infl.fun(A=Store.Influence[[s]],Whr[s],YLIM=YLIM)
    if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    #Daily
    dumy=unlist(lapply(Store.Influence.daily, "[[", match("Annual.Dev",names(Store.Influence.daily[[1]]))))
    YLIM=c(min(dumy),max(dumy))
    Compare.term.infl.fun(A=Store.Influence.daily[[s]],Whr.d[s],YLIM=YLIM)
    if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    
    mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Influence",side=2,line=1.25,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  #Show dusky and sandbar CDI
  SCLR=1.75
  CDI.sp=c("Dusky shark","Sandbar shark")
  fn.fig("Figure 6.Dusky_Sandbar_CDI",2400, 1200)
  par(mfcol=c(2,6),mar=c(1,2,.5,.5),oma=c(2,1.5,1,.2),las=1,mgp=c(1.9,.5,0))
  for(s in 1:length(CDI.sp))
  {
    #Monthly
    if(CDI.sp[s]=="Dusky shark") Terms="vessel"
    if(CDI.sp[s]=="Sandbar shark") Terms=c("vessel","blockx")
    Term.Type=ifelse(Terms%in%c("vessel" , "month" , "blockx"),"CAT",NA)
    YLABs=Terms
    YLABs=ifelse(YLABs=="blockx","Block",ifelse(YLABs=="vessel","Vessel",NA))
    
    ii=match(CDI.sp[s],names(Stand.out))
    Fig.CDI.paper.fn(MOD=Stand.out[[ii]]$res,DAT=Stand.out[[ii]]$DATA,
                     Term.type=Term.Type,termS=Terms,SCALER=SCLR,YLABs,CxAx=.75) 
    
    if(CDI.sp[s]=="Sandbar shark")
    {
      Terms=c("vessel","blockx","month")
      Term.Type=ifelse(Terms%in%c("vessel" , "month" , "blockx"),"CAT",NA)
      YLABs=Terms
      YLABs=ifelse(YLABs=="blockx","Block",ifelse(YLABs=="vessel","Vessel",
                                                  ifelse(YLABs=="month","Month",NA)))
      ii=match(CDI.sp[s],names(Stand.out.daily))
      Fig.CDI.paper.fn(MOD=Stand.out.daily[[ii]]$res,DAT=Stand.out.daily[[ii]]$DATA,
                       Term.type=Term.Type,termS=Terms,SCALER=SCLR,YLABs,CxAx=.75) 
    }
    mtext(c("Financial year                   Coefficient "),side=2,line=.25,cex=1,las=3,outer=T)
    mtext(c("Dusky (monthly)          Sandbar (monthly)                                          Sandbar (daily)                                                                    "),
          3,line=-.25,cex=.75,outer=T)
  }
  dev.off()
  
}



#   4.22.11 Construct spatial standardised catch rates    #ACA: use blocks for monthly and lat long for daily
#1. fit glms to each specific zone
#2. then predict, etc
Pred.zone=Pred.zone.creep=Pred.zone.nrmlzd=Pred
Pred.daily.zone=Pred.daily.zone.creep=Pred.daily.zone.nrmlzd=Pred.daily
system.time({for(s in 1:N.species)
{
  #Monthly
  DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used[[s]])
  ZonEs=unique(DAT$zone)
  if(SPECIES.vec[s]=="Sandbar shark") ZonEs=subset(ZonEs,!ZonEs=="Zone2")   #no enough observations in zone2
  pred.temp=vector('list',length(ZonEs))
  names(pred.temp)=ZonEs
  pred.temp.crip=pred.temp.nrm=pred.temp

  for(z in 1:length(ZonEs))
  {
    #1. Fit model to zone data
    model=fn.stand(d=subset(DAT,zone==ZonEs[z]),Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                   efrt="km.gillnet.hours.c",Formula=Best.Model[[s]])

    #2. Predict for selected blocks
    d=model$DATA
    Pred.1=pred.fun(MOD=model$res,biascor="YES",PRED="finyear",Pred.type="link")

    #3. Apply creep
    Pred.1.c=Pred.1
    add.crp=Eff.creep$effort.creep[match(Pred.1.c$finyear,Eff.creep$finyear)]
    Pred.1.c$response=Pred.1.c$response*(1-add.crp)
    Pred.1.c$lower.CL=Pred.1.c$lower.CL*(1-add.crp)
    Pred.1.c$upper.CL=Pred.1.c$upper.CL*(1-add.crp)

    #4.Normalize
    Pred.1.c.norm=Pred.1.c
    Mn=mean(Pred.1.c.norm$response)
    Pred.1.c.norm$response=Pred.1.c.norm$response/Mn
    Pred.1.c.norm$lower.CL=Pred.1.c.norm$lower.CL/Mn
    Pred.1.c.norm$upper.CL=Pred.1.c.norm$upper.CL/Mn

    #5.Store
    pred.temp[[z]]=Pred.1
    pred.temp.crip[[z]]=Pred.1.c
    pred.temp.nrm[[z]]=Pred.1.c.norm
    rm(d,Pred.1,Pred.1.c,Pred.1.c.norm)
  }
  Pred.zone[[s]]=pred.temp
  Pred.zone.creep[[s]]=pred.temp.crip
  Pred.zone.nrmlzd[[s]]=pred.temp.nrm


  #Daily
  DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
  ZonEs=unique(DAT$zone)
  if(SPECIES.vec[s]=="Sandbar shark") ZonEs=subset(ZonEs,!ZonEs=="Zone2")   #no enough observations in zone2
  pred.temp=vector('list',length(ZonEs))
  names(pred.temp)=ZonEs
  pred.temp.crip=pred.temp.nrm=pred.temp

  for(z in 1:length(ZonEs))
  {
    #1. Fit model to zone data
    model=fn.stand(d=subset(DAT,zone==ZonEs[z]),Response="catch.target",RESPNS="LNcpue",
                   PREDS=Predictors_daily,efrt="km.gillnet.hours.c",Formula=Best.Model.daily[[s]])

    #2. Predict for selected blocks
    d=model$DATA
    Pred.1=pred.fun(MOD=model$res,biascor="YES",PRED="finyear",Pred.type="link")

    #3. Apply creep
    Pred.1.c=Pred.1
    add.crp=Eff.creep$effort.creep[match(Pred.1.c$finyear,Eff.creep$finyear)]
    Pred.1.c$response=Pred.1.c$response*(1-add.crp)
    Pred.1.c$lower.CL=Pred.1.c$lower.CL*(1-add.crp)
    Pred.1.c$upper.CL=Pred.1.c$upper.CL*(1-add.crp)

    #4.Normalize
    Pred.1.c.norm=Pred.1.c
    Mn=mean(Pred.1.c.norm$response)
    Pred.1.c.norm$response=Pred.1.c.norm$response/Mn
    Pred.1.c.norm$lower.CL=Pred.1.c.norm$lower.CL/Mn
    Pred.1.c.norm$upper.CL=Pred.1.c.norm$upper.CL/Mn

    #5.Store
    pred.temp[[z]]=Pred.1
    pred.temp.crip[[z]]=Pred.1.c
    pred.temp.nrm[[z]]=Pred.1.c.norm
    rm(d,Pred.1,Pred.1.c,Pred.1.c.norm)
  }
  Pred.daily.zone[[s]]=pred.temp
  Pred.daily.zone.creep[[s]]=pred.temp.crip
  Pred.daily.zone.nrmlzd[[s]]=pred.temp.nrm

}})



#   4.22.12 Export catch rates        MISSING: add other species relative..Pred.normlzd Pred.daily.normlzd
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Index")
Sel.vars=c("finyear","response","SE","lower.CL","upper.CL")
nams.Sel.vars=c("Finyear","Mean","CV","LOW.CI","UP.CI")
 
#Absolute scale 
    #4.22.12.1 zones combined with NO efficiency creep
for (s in 1:N.species)
{
  #Standardised
  a=subset(Pred[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly_no.creep.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily_no.creep.csv",sep=""),row.names=F) 

  rm(a)
}

    #4.22.12.2 zones combined with efficiency creep
List.foly.nom.creep=List.foly.nom    #first add effort creep to foly 
for(s in 1:N.species)
{
  if(names(List.foly.nom.creep)[s]=="san") List.foly.nom.creep[[s]]$Foly=subset(List.foly.nom.creep[[s]]$Foly,FINYEAR%in%San.Yrs)
  add.crp=Eff.creep$effort.creep[match(List.foly.nom.creep[[s]]$Foly$FINYEAR,Eff.creep$finyear)]
  List.foly.nom.creep[[s]]$Foly$cpue=List.foly.nom.creep[[s]]$Foly$cpue*(1-add.crp)
}
for (s in 1:N.species)
{
  #Standardised
  a=subset(Pred.creep[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily.creep[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.csv",sep=""),row.names=F) 
  

  #Unstandardised
  write.csv(List.foly.nom.creep[[s]]$Foly,paste(SPECIES.vec[s],".annual.folly.csv",sep=""),row.names=F)  
  
  a=subset(Nominl.creep[[s]],select=c("finyear","response","lower.CL","upper.CL"))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.nominal.monthly.csv",sep=""),row.names=F)    
  
  a=subset(Nominl.daily.creep[[s]],select=c("finyear","response","lower.CL","upper.CL"))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.nominal.daily.csv",sep=""),row.names=F) 

  rm(a)
}

    #4.22.12.3 by zones without efficiency creep
for (s in 1:N.species)
{
  Zn=names(Pred.zone[[s]])
  for(z in 1:length(Zn))
  {
      #Standardised
    a=subset(Pred.zone[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
  
    a=subset(Pred.daily.zone[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
  
    rm(a)
  }
}

    #4.22.12.4 by zones with efficiency creep
for (s in 1:N.species)
{
  Zn=names(Pred.zone[[s]])
  for(z in 1:length(Zn))
  {
    #Standardised
    a=subset(Pred.zone.creep[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.",Zn[z],".csv",sep=""),row.names=F) 
    
    a=subset(Pred.daily.zone.creep[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.",Zn[z],".csv",sep=""),row.names=F) 
    
    rm(a)
  }
}


#Relative scale   
    #4.22.12.5 zones combined with efficiency creep
for (s in 1:N.species)
{
  #Standardised
  a=subset(Pred.normlzd[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly_relative.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily.normlzd[[s]],select=Sel.vars)   
  names(a)=nams.Sel.vars
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily_relative.csv",sep=""),row.names=F) 
  rm(a)
}

    #4.22.12.6 by zones with efficiency creep
for (s in 1:N.species)
{
  Zn=names(Pred.zone[[s]])
  for(z in 1:length(Zn))
  {
    #Standardised
    a=subset(Pred.zone.nrmlzd[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.",Zn[z],"_relative.csv",sep=""),row.names=F) 
    
    a=subset(Pred.daily.zone.nrmlzd[[s]][[z]],select=Sel.vars)   
    names(a)=nams.Sel.vars
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.",Zn[z],"_relative.csv",sep=""),row.names=F) 
    
    rm(a)
  }
}


##############--- 5. REPORT SECTION FROM 1.Manipulate data.R---###################
plot.cpue.paper.figures="NO"
if (plot.cpue.paper.figures=="YES")
{
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
  
  #Some functions
  
  PLOT.fn=function(DATA1,DATA2,max1,max2,LEG,CL2)   #function plot vessels and blocks
  {
    plot(1:NN.monthly,DATA1,ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1,ylim=c(0,max1)
         ,cex.axis=.75,lwd=1.1)
    axis(1,at=1:NN.monthly,labels=F,tck=-0.015)
    
    par(new=T)
    plot(1:NN.monthly,DATA2,col=CL2,type='o',pch=19,axes=F,ann='F',cex=1,lwd=1.1,ylim=c(0,max2))
    #  axis(4,at=pretty(DATA2),labels=pretty(DATA2),las=2,cex.axis=.75)
    axis(4,at=seq(0,max2,40),labels=seq(0,max2,40),las=2,cex.axis=.75,col.axis=CL2)
    legend("topleft",LEG,bty='n',cex=.75)
  }
  
  fn.eff.probs=function(gear,eff.var,eff.var1,letra,CLs)
  {
    fn.wrong=function()
    {
      ID1=match(corr,names(dummy))
      dummy$VAR=ifelse(!(dummy[,ID]==dummy[,ID1]),"Wrong","OK")
      dummy$VAR=with(dummy,ifelse(is.na(VAR),"Wrong",VAR))
      Agg=table(dummy[,ID2],dummy$VAR)
      Agg=as.matrix(round(100*Agg/rowSums(Agg),4))
      if(ncol(Agg)==2)Tab=cbind(Agg[,1],Agg[,2])
      if(ncol(Agg)==1)Tab=cbind(Agg[,1],rep(0,length(Agg[,1])))
      colnames(Tab)=c("ok","wrong")
      return(Tab)
    }
    
    #monthly
    dummy=subset(Effort.monthly,!FINYEAR%in%Daily.l.years & METHOD==gear & !(BLOCKX%in%Estuaries))
    corr=paste(eff.var,".c",sep="")
    dummy=dummy[,match(c("FINYEAR","METHOD","Same.return",eff.var,corr),names(dummy))]
    dummy=dummy[!duplicated(dummy$Same.return),]
    THIS="FINYEAR"
    eff.var.1=eff.var
    ID=match(eff.var,names(dummy))
    ID2=match(THIS,names(dummy))
    
    
    Tab.m=fn.wrong()
    
    #daily
    dummy=subset(Effort.daily,method==gear & !(blockx%in%Estuaries))
    corr=paste(eff.var1,".c",sep="")
    dummy=dummy[,match(c("finyear","method","date","ID",eff.var1,corr),names(dummy))]
    dummy$ID=with(dummy,paste(date,ID))
    dummy=dummy[!duplicated(dummy$ID),]
    ID=match(eff.var1,names(dummy))
    ID2=match("finyear",names(dummy))
    
    Tab.d=fn.wrong()
    
    TAB=rbind(Tab.m,Tab.d)
    rownames(TAB)=NULL
    barplot(rbind(TAB[,1],TAB[,2]),beside=F,col=CLs,cex.axis=1.25,cex.lab=1.25,
            ylim=c(0,110),yaxs="i",xaxs="i",legend.text=paste(letra),          
            args.legend=list(x = "topleft",cex=1,bty='n',
                             fill="transparent",border='transparent'),ylab="")
    box()
    AXIS1()
    AXIS2()
    return(data.frame(OK=mean(TAB[,1]),wrong=mean(TAB[,2])))
  }
  
  fn.explore=function(DATA,NAMES)
  {
    #remove initial years when sandbar was not reported
    if(NAMES=="Sandbar shark")
    {
      DATA=subset(DATA,!(FINYEAR%in%c("1975-76","1976-77","1977-78",
                                      "1978-79","1979-80","1980-81","1981-82","1982-83","1983-84","1984-85")))
    }
    
    #Remove vessel levels that don't occur in data
    DATA$VESSEL=DATA$VESSEL[, drop=TRUE]
    
    DATA$CPUE=with(DATA,Catch.Target/Km.Gillnet.Days.c)
    DATA$SP.Target=ifelse(DATA$CPUE>0,1,0)
    
    #Tables
    TABLE1=sort(with(DATA,table(VESSEL)))
    TABLE1=data.frame(VESSEL=names(TABLE1),Count=as.numeric(TABLE1))
    
    TABLE6=aggregate(Catch.Target~VESSEL,data=DATA,mean)
    names(TABLE6)[2]="Mean Catch"
    TABLE6.1=aggregate(Catch.Target~VESSEL,data=DATA,sd)
    names(TABLE6.1)[2]="SD Catch"
    TABLE6=merge(TABLE6,TABLE1,by="VESSEL")
    TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
    TABLE6=TABLE6[order(TABLE6$Count),]
    
    
    #cumulative catch
    TABLE12=aggregate(Catch.Target~VESSEL,data=DATA,sum)
    TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
    TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
    TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
    
    TABLE13=aggregate(Catch.Target~BLOCKX,data=DATA,sum)
    TABLE13=TABLE13[order(-TABLE13$Catch.Target),]
    TABLE13$CumCatch=cumsum(TABLE13$Catch.Target)
    TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$Catch.Target),2)
    
    #cumulative records
    TABLE16=rev(sort(table(DATA$VESSEL)))
    TABLE16=data.frame(Records=as.numeric(TABLE16),VESSEL=names(TABLE16))
    TABLE16$CumRecords=cumsum(TABLE16$Records)
    TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
    
    TABLE17=rev(sort(table(DATA$BLOCKX)))
    TABLE17=data.frame(Records=as.numeric(TABLE17),BLOCKX=names(TABLE17))
    TABLE17$CumRecords=cumsum(TABLE17$Records)
    TABLE17$PerCumRecords=round(TABLE17$CumRecords*100/sum(TABLE17$Records),2)
    
    
    
    return(list(Ves.Cum.Ca=TABLE12$PerCumCatch,Block.Cum.Ca=TABLE13$PerCumCatch,
                TABLE17=TABLE17,TABLE16=TABLE16,Mean_catch_Vessel=TABLE6))
  }
  
  fun.rec.per.ves=function(DATA)
  {
    Count.ves.records=table(DATA$Count)
    Count.ves.records=c(subset(Count.ves.records,as.numeric(names(Count.ves.records))<Min.rec.ves),sum(
      subset(Count.ves.records,as.numeric(names(Count.ves.records))<Min.rec.ves)))
    names(Count.ves.records)[length(Count.ves.records)]=paste(Min.rec.ves-1,"+",sep="")
    b=barplot(Count.ves.records,xaxt='n',yaxt='n')
    axis(1,at=b,labels=F,tck=-0.016)
    axis(2,at=seq(0,300,100),labels=seq(0,300,100),tck=-0.016,cex=1.25,las=1)
    axis(1,at=b[seq(2,(length(b)-2),by=2)],labels=names(Count.ves.records)[seq(2,(length(b)-2),by=2)],
         tck=-0.032,cex=1.25)
    mtext("Number of records per vessel",side=1,line=1.2,font=1,las=0,cex=1.125,outer=F)
    mtext("Count",side=2,line=1.75,font=1,las=0,cex=1.25,outer=F)
    axis(1,at=b[length(b)],labels=paste(">",19,sep=""),tck=-0.032,cex.axis=1.15)
    
  }
  
  fn.Figure5=function(DATA,NAMES,SP)
  {
    #remove initial years when sandbar was not reported
    if(NAMES=="Sandbar shark")
    {
      DATA=subset(DATA,!(FINYEAR%in%c("1975-76","1976-77","1977-78",
                                      "1978-79","1979-80","1980-81","1981-82","1982-83","1983-84","1984-85")))
    }
    
    #Remove vessel levels that don't occur in data
    DATA$VESSEL=DATA$VESSEL[, drop=TRUE]
    
    DATA=subset(DATA,SPECIES==SP)
    DATA$CPUE=with(DATA,LIVEWT.c/Km.Gillnet.Hours.c)
    DATA$SP.Target=ifelse(DATA$CPUE>0,1,0)
    DATA$Catch.Target=DATA$LIVEWT.c
    
    #Tables
    TABLE1=sort(with(DATA,table(VESSEL)))
    
    TABLE6=aggregate(Catch.Target~VESSEL,data=DATA,mean)
    names(TABLE6)[2]="Mean Catch"
    TABLE6.1=aggregate(Catch.Target~VESSEL,data=DATA,sd)
    names(TABLE6.1)[2]="SD Catch"
    TABLE6=merge(TABLE6,data.frame(VESSEL=names(TABLE1),Count=as.numeric(TABLE1)),by="VESSEL")
    TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
    TABLE6=TABLE6[order(TABLE6$Count),]
    
    
    #cumulative catch
    TABLE12=aggregate(Catch.Target~VESSEL,data=DATA,sum)
    TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
    TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
    TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
    
    TABLE13=aggregate(Catch.Target~BLOCKX,data=DATA,sum)
    TABLE13=TABLE13[order(-TABLE13$Catch.Target),]
    TABLE13$CumCatch=cumsum(TABLE13$Catch.Target)
    TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$Catch.Target),2)
    
    #cumulative records
    TABLE16=rev(sort(table(DATA$VESSEL)))
    TABLE16=data.frame(Records=as.numeric(TABLE16),VESSEL=names(TABLE16))
    TABLE16$CumRecords=cumsum(TABLE16$Records)
    TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
    
    TABLE17=rev(sort(table(DATA$BLOCKX)))
    TABLE17=data.frame(Records=as.numeric(TABLE17),BLOCKX=names(TABLE17))
    TABLE17$CumRecords=cumsum(TABLE17$Records)
    TABLE17$PerCumRecords=round(TABLE17$CumRecords*100/sum(TABLE17$Records),2)
    
    
    
    return(list(Ves.Cum.Ca=TABLE12$PerCumCatch,Block.Cum.Ca=TABLE13$PerCumCatch,
                TABLE17=TABLE17,TABLE16=TABLE16,Mean_catch_Vessel=TABLE6))
  }
  
  #Extract monthly records by Year-Month-Vessel-Block
  
  These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Days.inv",
                  "Km.Gillnet.Days.c","zone","MONTH","BLOCKX")
  
  Effort.data.fun=function(DATA,target,ktch)
  {
    #remove record if no effort data
    ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
    ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
    if(length(ID)>0)DATA=DATA[-ID,]
    
    # remove nonsense lat
    DATA=subset(DATA,LAT>=(-36))
    
    #calculate effort
    Match.these.eff=match(These.efforts,names(DATA))
    Effort.data=DATA[,Match.these.eff]
    Effort.data=aggregate(cbind(Km.Gillnet.Days.inv,Km.Gillnet.Days.c)~zone+
                            FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,max)
    Effort.data=aggregate(cbind(Km.Gillnet.Days.inv,Km.Gillnet.Days.c)~zone+
                            FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,sum)
    
    #target species catch 
    ID=match(c(ktch),colnames(DATA))
    DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
    
    
    #catch targeted at other species
    DATA$Catch.Gummy=with(DATA,ifelse(SPECIES==17001,DATA[,ID],0))
    DATA$Catch.Whiskery=with(DATA,ifelse(SPECIES==17003,DATA[,ID],0))
    DATA$Catch.Dusky=with(DATA,ifelse(SPECIES%in%c(18003,18001),DATA[,ID],0))
    DATA$Catch.Sandbar=with(DATA,ifelse(SPECIES==18007,DATA[,ID],0))
    DATA$Catch.Scalefish=with(DATA,ifelse(SPECIES%in%188000:599001,DATA[,ID],0))
    DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
    
    #reshape catch data
    TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Scalefish,
                          Catch.Dusky,Catch.Sandbar,Catch.Total)~MONTH+
                      FINYEAR+BLOCKX+VESSEL+Same.return+LAT+LONG+
                      YEAR.c,data=DATA,sum,na.rm=T)
    TABLE=TABLE[order(TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    
    #merge catch and effort
    dat=merge(TABLE,Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"),all.x=T)
    
    #create "other shark catch" variable
    dat$Catch.other.shk=NA
    if(target[1]==17003)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Sandbar
    if(target[1]==17001)dat$Catch.other.shk=dat$Catch.Whiskery
    if(target[1]%in%c(18003,18001))dat$Catch.other.shk=dat$Catch.Whiskery+dat$Catch.Sandbar
    if(target[1]==18007)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Whiskery
    
    #recalculate 60 by 60 blocks
    dat$BLOCKX.orignl=dat$BLOCKX
    dat$BLOCKX=as.numeric(substr(dat$BLOCKX,1,4))
    
    
    return(list(dat=dat,prop.with.catch=prop.with.catch))
  }
  
  DATA.list.LIVEWT.c=vector('list',length=N.species)
  names(DATA.list.LIVEWT.c)=names(Species.list)
  
  
  #Create data sets for plotting cpue paper figures
  for ( i in 1:N.species)DATA.list.LIVEWT.c[[i]]=Effort.data.fun(Species.list[[i]],TARGETS[[i]],"LIVEWT.c")$dat
  
  #Create figures 1 to 5
  if (plot.cpue.paper.figures=="YES")
  {
    tiff(file="Figure 1. Map.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mar=c(2,2,2,2),oma=c(1,1,1,1))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
    polygon(x=c(116.5,116.5,112,112),y=c(-26.5,-33,-33,-26.5),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
    polygon(x=c(116.5,116.5,112,112),y=c(-33,-37,-37,-33),lwd=1.5,col=rgb(.3,.3,.3,alpha=.5))
    polygon(x=c(129,129,116.5,116.5),y=c(-30,-37,-37,-30),lwd=1.5,col=rgb(.7,.7,.7,alpha=.2))
    
    axis(side = 1, at =seq(LONGG[1],LONGG[length(LONGG)],length.out = 7+6*(length(LONGG)-2)), labels = F, tcl = 34,lty=3,col="grey60")
    axis(side = 4, at = seq(LATT[1],LATT[length(LATT)],length.out = 7+6*(length(LATT)-2)), labels = F,tcl =34,lty=3,col="grey30")
    axis(side = 1, at =LONGG, labels = F, tcl = 34,lty=1,col="grey30")
    axis(side = 4, at = LATT, labels = F,tcl =34,lty=1,col="grey30")
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=c(-37,-25),xlim=South.WA.long, zlim=c(-1,-300),
                                 nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
    
    par(new=T,mar=c(2,2,2,2),oma=c(1,1,1,1))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
    axis(side = 1, at =seq(112,129,2), labels = seq(112,129,2), tcl = .35,las=1,cex.axis=1.25,padj=-1.25)
    axis(side = 2, at = seq(-36,-25,2), labels = -seq(-36,-25,2),tcl = .35,las=2,cex.axis=1.25,hadj=.3)
    text(116.73,Perth[2],("Perth"),col="black", cex=1.1)
    points(115.86,-31.95,pch=19,cex=1.5)
    text(116.73,-33.55,("Bunbury"),col="black", cex=1.1)
    points(115.6,-33.55,pch=19,cex=1.5)
    text(117.7,-34.75,("Albany"),col="black", cex=1.1)
    points(117.8,-35,pch=19,cex=1.5)
    text(122,-33.62,("Esperance"),col="black", cex=1.1)
    points(121.9,-33.86,pch=19,cex=1.5)
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1.2,las=3,cex=1.75)
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.75,cex=1.75)
    
    text(113.5,-30.5,("WCDGDLF"),col="black", cex=1.4)
    text(114,-34.75,("JASDGDLF"),col="black", cex=1.4) 
    text(114,-35.55,("(Zone 1)"),col="black", cex=1.4)
    text(122,-34.75,("JASDGDLF"),col="black", cex=1.4)
    text(122,-35.55,("(Zone 2)"),col="black", cex=1.4)
    
    par(fig=c(.5,.92,.5,.92), new = T,mgp=c(.1,.4,0))
    plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
            col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
    text(134,-22.5,("Australia"),col="black", cex=2)
    dev.off()
    
    
    #Appendix 1.                
    
    #Stack-up catch plot
    these.ones=c(22999,17001,17003,18007,18003)
    names(these.ones)=c("SHARK, OTHER","SHARK, GUMMY","SHARK, WHISKERY",
                        "SHARK, THICKSKIN (SANDBAR)","SHARK, BRONZE WHALER")
    #Monthly
    FINYEAR.monthly=sort(unique(Data.monthly.GN$FINYEAR))
    Mn.yr=subset(FINYEAR.monthly,!FINYEAR.monthly %in% FINYEAR.daily)
    NN.monthly=length(Mn.yr)
    STORE=matrix(nrow=NN.monthly,ncol=length(these.ones))
    colnames(STORE)=sort(these.ones)
    for(i in 1:NN.monthly)
    {
      Dat=subset(Data.monthly.GN, FINYEAR==FINYEAR.monthly[i])    
      test=aggregate(LIVEWT~FINYEAR+Spec.old,data=Dat,sum,na.rm=T)
      test$SP=ifelse(test$Spec.old%in%these.ones,test$Spec.old,"OTHER")
      test=aggregate(LIVEWT~FINYEAR+ SP,data=test,sum,na.rm=T)
      prop=100*test[,3]/sum(test[,3])
      names(prop)=test[,2]
      if(!(names(prop)[4]=="18007"))prop=c(prop[1:3],0,prop[4])
      if(names(prop)[4]=="18007")prop=prop[1:5]
      STORE[i,]=prop
    }
    
    #Daily
    NN.daily=length(FINYEAR.daily)
    STORE.daily=matrix(nrow=NN.daily,ncol=length(these.ones))
    colnames(STORE.daily)=sort(these.ones)
    for(i in 1:NN.daily)
    {
      Dat=subset(Data.daily.GN, FINYEAR==FINYEAR.daily[i])    
      Dat$Spec.old=Dat$SPECIES
      test=aggregate(LIVEWT~FINYEAR+Spec.old,data=Dat,sum,na.rm=T)
      test$SP=ifelse(test$Spec.old%in%these.ones,test$Spec.old,"OTHER")
      test=aggregate(LIVEWT~FINYEAR+ SP,data=test,sum,na.rm=T)
      prop=100*test[,3]/sum(test[,3])
      names(prop)=test[,2]
      if(!(names(prop)[4]=="18007"))prop=c(prop[1:3],0,prop[4])
      if(names(prop)[4]=="18007")prop=prop[1:5]
      STORE.daily[i,]=prop
    }
    STORE=rbind(STORE,STORE.daily)
    
    
    #Stack-up barplot of effort problems
    COL.BAR=c("white","grey35","grey55","grey75","black")
    NN=NN.monthly+NN.daily
    AXIS1=function()axis(1,at=b,labels=F,tck=-0.016)
    AXIS2=function()axis(1,at=b[seq(1,NN,by=5)],labels=F,tck=-0.03)
    AXIS3=function()axis(1,at=b[seq(1,NN,by=5)],labels=FINYEAR.monthly[seq(1,NN,by=5)],tck=-0.035,cex.axis=1.25)
    
    tiff(file="Appendix 1. Data problems_All.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    layout(matrix(c(1,1,2:5), 3, 2, byrow = TRUE))
    par(mar=c(2,.75,.25,2),oma=c(2,4,.25,.01),las=1,mgp=c(1,.75,0))
    
    #Catch (species id problem)
    b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
              ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
              args.legend=list(x = "topleft",horiz=T,cex=1.7,pt.cex=1.7,bty='n',inset=c(0, -0.025)),ylab="")
    box()
    AXIS1()
    AXIS3()
    
    #BDAYS                                        
    print(fn.eff.probs("GN","BDAYS","bdays","",c("grey85","black"))  )
    mtext("Number of days fished per month",3,line=-1.375,cex=1.25)
    
    #HOURS
    print(fn.eff.probs("GN","HOURS","hours","",c("grey85","black")))
    mtext("Number of hours fished per day",3,line=-1.375,cex=1.25)
    
    #SHOTS
    print(fn.eff.probs("GN","SHOTS","shots","",c("grey85","black")))
    AXIS3()
    mtext("Number of shots per day",3,line=-1.375,cex=1.25)
    
    #NETLEN     
    print(fn.eff.probs("GN","NETLEN","netlen","",c("grey85","black")))
    AXIS3()
    mtext("Net length per shot",3,line=-1.375,cex=1.25)
    
    mtext("Financial year",side=1,line=0.8,font=1,las=0,cex=1.75,outer=T)
    mtext("Percentage",side=2,line=2,font=1,las=0,cex=1.75,outer=T)
    dev.off()
    
    
    tiff(file="Appendix 1. Data problems.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(2,1),mar=c(1,3,.1,.1),oma=c(2.5,.5,.1,.5),las=1,mgp=c(1.8,.6,0))
    
    #Catch (species id problem)
    b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
              ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
              args.legend=list(x = "topleft",horiz=T,cex=1.17,pt.cex=1.25,bty='n',inset=c(0, -0.025)),ylab="")
    box()
    AXIS1()
    AXIS2()
    
    # #NETLEN     
    fn.eff.probs("GN","NETLEN","netlen","",c("grey85","black"))
    AXIS3()
    
    mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.75,outer=T)
    mtext("Percentage",side=2,line=-0.75,font=1,las=0,cex=1.75,outer=T)
    dev.off()
    
    
    tiff(file="Appendix 1. Data problems_catch_only.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(1,1),mar=c(1,3,.15,.1),oma=c(3,.5,.15,.65),las=1,mgp=c(1.8,.9,0))
    
    #Catch (species id problem)
    b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
              ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
              args.legend=list(x = "topleft",horiz=T,cex=1.17,pt.cex=1.25,bty='n',inset=c(0, -0.015)),ylab="")
    box()
    AXIS1()
    AXIS2()
    AXIS3()
    mtext("Financial year",side=1,line=1.5,font=1,las=0,cex=1.75,outer=T)
    mtext("Percentage",side=2,line=-0.75,font=1,las=0,cex=1.75,outer=T)
    dev.off()
    
    
    tiff(file="Appendix 1. Data problems_effort_only.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(2,2),mar=c(1,3.2,.1,.1),oma=c(2.5,.5,.1,.65),las=1,mgp=c(1.8,.6,0))
    
    #BDAYS                                        
    print(fn.eff.probs("GN","BDAYS","bdays","(a)",c("grey85","black"))  )
    
    #HOURS
    print(fn.eff.probs("GN","HOURS","hours","(b)",c("grey85","black")))
    AXIS3()
    #SHOTS
    print(fn.eff.probs("GN","SHOTS","shots","(c)",c("grey85","black")))
    
    #NETLEN     
    print(fn.eff.probs("GN","NETLEN","netlen","(d)",c("grey85","black")))
    AXIS3()
    
    mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.75,outer=T)
    mtext("Percentage",side=2,line=-1,font=1,las=0,cex=1.75,outer=T)
    dev.off()
    
    
    
    #Appendix 2. Flowchart diagram 
    source("C:/Matias/Analyses/SOURCE_SCRIPTS/flow_chart.R")
    #text
    LABELS <- list("raw data",
                   c("Exclude records from estuaries, non-gillnet gear","and school shark and dogfish targeting"),
                   c("Valid catch composition"),
                   c("Correct catch but do not use record"," in catch and effort standardisation"),
                   "Valid effort",
                   c("Correct effort but do not use record","in catch and effort standardisation"),
                   c("Financial year > 1989-90"),
                   c("Adjust incomplete","catch and effort"),
                   c("Record within effective area"),
                   c("Do not use record in","catch and effort standardisation"),
                   c("Standardise catch and effort","Construct abundance index"),
                   c("Adjust for increase in","fishing efficiency"))
    
    #type of shape
    #Note:
    #oval: start and terminal points
    #square or round: process
    #diammond: decision (yes/no)
    
    SHAPES=c("oval","round","diamond","round","diamond","round","diamond",
             "round","diamond","round","round","oval")
    
    #Shape coordinates
    MaInX=0.725
    X.COOR=c(rep(MaInX,3),MaInX*0.35,MaInX,MaInX*0.35,MaInX,MaInX*0.35,
             MaInX,MaInX*0.35,MaInX,MaInX)
    N.labl=length(LABELS)
    Y.COOR=rep(NA,N.labl)
    Y.COOR[1]=0.975
    delta=c(rep(0.115,3),rep(0.07,8),0.115)
    for(q in 2:N.labl) Y.COOR[q]=Y.COOR[q-1]-delta[q]
    
    #Shape size
    X.size=c(.08,.19,.13,.165,.09,.165,.125,.125,.125,.125,.13,.125)
    
    #arrows
    ArROW=c(rep("Straight",2),rep(c("Side.left","Side.right"),4),"Straight","Side.back")
    
    tiff(file="Appendix 2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mar=c(0.1,0.1,0.1,0.1))
    fn.flow.chart(lab=LABELS,SHPE=SHAPES,X.coor=X.COOR,Y.coor=Y.COOR,SX=X.size,ARRw=ArROW,CEX=.9,n=8,n1=9)
    dev.off()  
    
    
    #Appendix 3. Effort dynamics (expansion and contraction)
    NN.monthly=NN
    Lat.seq=c(-26,-28,-30,-32,-34,-36)
    
    tiff(file="Appendix 3. Effort dynamics.gillnets.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=(1+length(DATA.lista)),MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(.1, 0.5, 0))
    #number of blocks and vessels per year
    PLOT.fn(BLKS.YEAR,Effort.expan$N.ves.yr,60,180,"",CL2="grey65")
    axis(1,at=seq(1,NN.monthly,5),labels=FINYEAR.monthly[seq(1,NN.monthly,5)],tck=-0.03,cex.axis=.9)
    mtext("Financial year",side=1,line=1.5,font=1,las=0,cex=1,outer=F)
    mtext("Number of blocks fished",side=2,line=1.35,font=1,las=0,cex=.85,outer=F)
    mtext("Number of licence holders fishing",side=4,line=1.3,las=3,cex=.75,outer=F,col="grey65")
    
    #effort by block per calendar year groups
    for (i in 1:length(DATA.lista))
    {
      DATA=DATA.lista[[i]][-which(duplicated(DATA.lista[[i]]$Same.return)),]
      
      fn.eff.plot(DATA,tcl.1=0,tcl.2=0,EffortBreakSS)
      mtext(Yr.range[i],side=3,line=-1.25,cex=.95)
      axis(side = 1, at =Long.seq, labels = F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
      if(i%in%c(10,7,8)) axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      if(i%in%c(3,6,9)) axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
      if(i==8) color.legend(126,-26,129,-30.5,round(EffortBreakSS,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
    }  
    mtext("Latitude (?S)",side=2,line=0.4,las=3,cex=1.25,outer=T)
    mtext("Longitude (?E)",side=1,line=0.75,cex=1.25,outer=T)  
    dev.off()
    
    
    
    #Appendix 4. Cumulative catch and vessels
    SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
    id=70
    Min.rec.ves=20  #minimum accepted number of records to keep a vessel
    
    
    Data.Summary=vector('list',length=N.species)
    
    for ( i in 1:N.species)Data.Summary[[i]]=fn.explore(DATA.list.LIVEWT.c[[i]],SPECIES.vec[i])
    Data.Fig5=vector('list',length=N.species)
    spe=c(17003,17001,18003,18007)
    for ( i in 1:N.species)Data.Fig5[[i]]=fn.Figure5(Species.list[[i]],SPECIES.vec[i],spe[i])
    
    tiff(file="Appendix 4.Cummulative.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    layout(matrix(c(1,2,1,2,3,4,3,4), 4, 2, byrow = TRUE))
    par(mar=c(3,.75,.1,3),oma=c(1,4,.1,.01),las=1,mgp=c(1,.9,0))
    line.col=c("black","grey20","grey55","grey80")
    line.type=c(1,2,1,1)
    plot(Data.Fig5[[1]]$Block.Cum.Ca[1:id],ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
    for (i in 2:N.species)lines(Data.Fig5[[i]]$Block.Cum.Ca[1:id],col=line.col[i],lwd=3,lty=line.type[i])
    axis(1,1:length(Data.Fig5[[1]]$Block.Cum.Ca[1:id]),labels=F,tck=-0.015)
    axis(1,seq(5,length(Data.Fig5[[1]]$Block.Cum.Ca[1:id]),5),
         labels=F,tck=-0.03)
    mtext("Cumulative catch (%)",side=2,line=2.75,font=1,las=0,cex=1.5,outer=F)
    
    plot(Data.Fig5[[1]]$Ves.Cum.Ca,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
    for (i in 2:N.species)lines(Data.Fig5[[i]]$Ves.Cum.Ca,col=line.col[i],lwd=3,lty=line.type[i])
    axis(1,seq(0,length(Data.Fig5[[1]]$Ves.Cum.Ca),10),labels=F,tck=-0.015)
    axis(1,seq(100,length(Data.Fig5[[1]]$Ves.Cum.Ca),100),labels=F,tck=-0.03)
    legend("bottomright",SPECIES.vec,bty='n',cex=1.5,col=line.col,lwd=3,lty=line.type)
    
    plot(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]*.99,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
    for (i in 2:N.species)lines(Data.Fig5[[i]]$TABLE17$PerCumRecords[1:id],col=line.col[i],lwd=3,lty=line.type[i])
    axis(1,1:length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),labels=F,tck=-0.015)
    axis(1,seq(5,length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),5),
         labels=seq(5,length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),5),tck=-0.03,cex.axis=1.65)
    mtext("Number of blocks",side=1,line=2.5,font=1,las=0,cex=1.5,outer=F)
    mtext("Cumulative records (%)",side=2,line=2.75,font=1,las=0,cex=1.5,outer=F)
    
    
    plot(Data.Fig5[[1]]$TABLE16$PerCumRecords*.99,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
    for (i in 2:N.species)lines(Data.Fig5[[i]]$TABLE16$PerCumRecords,col=line.col[i],lwd=3,lty=line.type[i])
    axis(1,seq(0,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),10),labels=F,tck=-0.015)
    axis(1,seq(100,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),100),
         labels=seq(100,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),100),tck=-0.03,cex.axis=1.65)
    mtext("Number of vessels",side=1,line=2.5,font=1,las=0,cex=1.5,outer=F)
    
    
    #Plot of records per vessel
    par(fig=c(0.60,1.00,.1,0.35), new = T,mgp=c(.25,.2,0),las=1)
    fun.rec.per.ves(Data.Fig5[[1]]$Mean_catch_Vessel)
    dev.off()
    
    
    #Plot effective area
    Dusky=c(X1=South.WA.long[1],X2=Dusky.range[2],Y1=South.WA.lat[1],Y2=Dusky.range[1])
    Sandbar=c(X1=South.WA.long[1],X2=Sandbar.range[2],Y1=South.WA.lat[1],Y2=Sandbar.range[1])
    Whiskery=c(X1=South.WA.long[1],X2=Whiskery.range[2],Y1=South.WA.lat[1],Y2=Whiskery.range[1])
    Gummy=c(X1=Gummy.range[1],X2=Gummy.range[2],Y1=South.WA.lat[1],Y2=-31.6)
    LISta=list(Whiskery=Whiskery,Gummy=Gummy,Dusky=Dusky,Sandbar=Sandbar)
    fn.show=function(X1,X2,Y1,Y2)
    {
      plotmap(a,b,PLATE,"grey85",South.WA.long,South.WA.lat)
      polygon(c(X1,X2,X2,X1),c(Y2,Y2,Y1,Y1),col='grey35',border="transparent")
      par(new=T)
      plotmap(a,b,PLATE,"grey85",South.WA.long,South.WA.lat)
    }
    jpeg(file="Effective_area.jpeg",width = 2400, height = 2400,units = "px", res = 300)
    par(mfcol=c(2,2),mai=c(.5,.5,.2,.1),oma=c(.1,.1,.65,.1))
    for(x in 1:length(LISta))
    {
      fn.show(X1=LISta[[x]][[1]],X2=LISta[[x]][[2]],Y1=LISta[[x]][[3]],Y2=LISta[[x]][[4]]) 
      mtext(names(LISta)[x],3,cex=2)
      At=c(LISta[[x]][[1]],LISta[[x]][[2]])
      axis(side = 1, at =At, labels = At, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      At=c(LISta[[x]][[3]],LISta[[x]][[4]])
      axis(side = 2, at =At , labels = -At,tcl = .35,las=2,cex.axis=1,hadj=.65)
    }
    dev.off()  
    
  }
}

