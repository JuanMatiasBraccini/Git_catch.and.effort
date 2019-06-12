#--------- CPUE STANDARDISATIONS OF DIFFERENT DATASETS ---------#

#NOTE:  This script standardises the catch and effort data for the 4 commercial shark species,

#       To update Mean Freo Sealevel each year, extract csv data from "http://uhslc.soest.hawaii.edu/data/download/fd"
#       and run the script "Get.Freo.R" to get the mean monthly value from the daily records


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
#           4.11 Add SOI, Freo and moon phase to SNo 
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
#               4.22.11 Show daily_nfish
#               4.22.12 Construct spatial standardised catch rates
#               4.22.13 Export catch rates



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
library(corrplot)

library('cluster')
library(factoextra) #for plotting

library(cede)
library(gridExtra)

library(glmulti)  #model selection

library(fitdistrplus)  #select distribution

library(coefplot)  #visualise coefficients

library(emmeans)  #for model predictions


options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)   


setwd("C:/Matias/Analyses/SOURCE_SCRIPTS")
source("Delta_lognormal.R")
source("Delta_gamma.R")
#source("Bootstrap_Delta_Methods.R")
source("Compare.error.structure.R")
source("Deviance.explained.R")
source("Sorting.objects.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/MS.Office.outputs.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/fn.fig.R")
source("C:\\Matias\\R\\Sort ls by size.R")
source("C:\\Matias\\Analyses\\SOURCE_SCRIPTS\\Population dynamics\\Nominal_cpue_functions.R")


#----1. DATA SECTION-----#

#1.1 Import data

setwd("C:/Matias/Analyses/Catch and effort/Data_outs")

#monthly
Data.monthly.GN.whiskery=read.csv(file ="Data.monthly.GN.whiskery.csv",stringsAsFactors=FALSE)
Data.monthly.GN.gummy=read.csv(file ="Data.monthly.GN.gummy.csv",stringsAsFactors=FALSE)
Data.monthly.GN.dusky=read.csv(file ="Data.monthly.GN.dusky.csv",stringsAsFactors=FALSE)
Data.monthly.GN.sandbar=read.csv(file ="Data.monthly.GN.sandbar.csv",stringsAsFactors=FALSE)

#daily
Data.daily.GN.whiskery=read.csv(file ="Data.daily.GN.whiskery.csv",stringsAsFactors=FALSE)
Data.daily.GN.gummy=read.csv(file ="Data.daily.GN.gummy.csv",stringsAsFactors=FALSE)
Data.daily.GN.dusky=read.csv(file ="Data.daily.GN.dusky.csv",stringsAsFactors=FALSE)
Data.daily.GN.sandbar=read.csv(file ="Data.daily.GN.sandbar.csv",stringsAsFactors=FALSE)


#Block10 locations
BlOCK_10=read.csv("C:/Matias/Data/Mapping/Blocks_10NM.csv")
names(BlOCK_10)=c("block10","LAT","LONG")
Metro_BlOCK_10=subset(BlOCK_10, LAT>(-33) & LAT<=(-31) & LONG<116)

#Southern Oscillation Index
SOI=read.csv("C:/Matias/Data/SOI.1975_2013.csv")
names(SOI)[1]="Month.soi"

#Mean Freo sea level
Freo=read.csv("C:/Matias/Data/Freo_mean_sea_level.csv")  
names(Freo)[c(1,3)]=c('Year',"Freo")

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


#Vessel characteristics
#note: must update ""VesselGearSurveyData.csv" in C:\Matias\Data\Fishing power
#source("C:/Matias/Analyses/Catch and effort/Git_catch_and_effort/4.Vessel_characteristics.R")



#1.2 Control what parts of script are activated

#1.2.1 Data controls

#Control if doing separated analysis of monthly and daily records
Separate.monthly.daily="YES"

#Control how to aggregate daily records
#Use.Date="YES"    #Rory's approach (aggregating by DATE)
Use.Date="NO"     # aggregating by Same.return.SNo (i.e. SNo, DsNo and TSNo). 
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
do.cluster="YES"

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


#Control of showing relative cpue
Show.relative.index="YES"

#Control what index to export
#Export.relative="NO"
Export.relative="YES"



#----2. PARAMETERS SECTION-----#

Stand.eff=1000          #express all standardised catches as: catch per 1000 km gn d
Lg.Efrt=log(Stand.eff)

Red.th=1  #percentage reduction in deviance to accept term


#Indicative blocks
MIN.obs.BLK=5  # minimum number of years with observations for each block
MIN.obs.BLK.sens=0  
MIN.obs.BL10K=5 # minimum number of observation per block10

MIN.ktch=100 #minimum annual catch of target species (in kgs)

#Indicative vessles
Threshold.n.yrs=10
Threshold.n.yrs.monthly=10 #minimum number of years reporting the species
Threshold.n.yrs.daily=5 
Threshold.n.yrs.sens=0  #sensitivity
Threshold.n.vessls.per.yr=5  #keep years with a least 5 vessels

Per.int.bc=.1  #interval for grouping vessels (not used)
Per.int.sens=.05   #sensitivity test (not used)
Prop.ktch.exp=.9   #keep blocks and vessels that explain 90% of catch (not used)
Min.rec.ves=10  #minimum accepted number of records to keep a vessel (not used)
Min.rec.ves.sens=0   #sensitivity (not used)

#Qualification levels minimum number of years with positive records
Min.Vess.yr=5 #monthly
Min.Vess.yr.d=5 #daily

Fish.Pow=.02  #2% annual (i.e. 2%, 4%, 6%, etc) increase in fishing power prior 1995  #RORY

SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
SPvec=c(17003,17001,18003,18007)
names(SPvec)=SPECIES.vec

POP.GRW.RATE=c(0.13,0.37,0.02,0.02)    #population growth rate from demography (used for imputing recent yrs)
names(POP.GRW.RATE)=SPECIES.vec

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



#----3. FUNCTIONS SECTION-----#

smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

fn.n.vessels=function(dat,dat1)  return(length(unique(c(dat$VESSEL,dat1$VESSEL))))

fn.n.blk.caught=function(Dat,Dat1,SPEC)  #total number of bloks where shark caught
{
  Dat=subset(Dat,SPECIES%in%SPEC & LIVEWT.c>0)
  Dat1=subset(Dat1,SPECIES%in%SPEC & LIVEWT.c>0)
  return(unique(c(Dat$BLOCKX,Dat1$BLOCKX)))
}

fn.see.all.yrs.ves.blks=function(a,SP,what,Ves.sel.BC,Ves.sel.sens,BLK.sel.BC,BLK.sel.sens,Min.ktch)
{
  
  a$BLOCKX=substr(a$BLOCKX,1,4)
  a=subset(a,select=c(MONTH,YEAR.c,BLOCKX,VESSEL,Same.return,SPECIES,LIVEWT.c,Km.Gillnet.Days.c,Km.Gillnet.Hours.c,Reporter))
  All.ves=unique(as.character(a$VESSEL))
  All.blk=unique(as.character(a$BLOCKX))
  
  a=subset(a,Reporter=="good")
  CATCH.sp=aggregate(LIVEWT.c~YEAR.c+VESSEL,subset(a,SPECIES==SP),sum)
  CATCH.sp=reshape(CATCH.sp,v.names="LIVEWT.c",idvar="YEAR.c",timevar="VESSEL",direction="wide")
  CATCH.sp=CATCH.sp[order(CATCH.sp$YEAR.c),]
  
  Vess=substr(names(CATCH.sp)[2:ncol(CATCH.sp)],10,30)
  Yrs=CATCH.sp$YEAR.c
  Z=as.matrix(CATCH.sp[,-1])
  
  #Vess with > X years of records
  ZZ=Z
  ZZ[ZZ<Min.ktch]=NA
  ZZ[ZZ>=Min.ktch]=1
  
  Yrs.with.ktch=colSums(ZZ,na.rm=T)
  
  pdf(paste("C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Vessel_pos_records_by_yr/",paste(SP,what,sep=""),".pdf",sep="")) 
  
  #Ves.sel.BC
  par(mar=c(3,3.5,.8,.8))
  WHICh=which(Yrs.with.ktch>Ves.sel.BC)
  Z.this=ZZ[,WHICh]
  Ves.BC=Vess[WHICh]
  ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
  Z.this=Z.this[,ID.sort]
  Ves.BC=Ves.BC[ID.sort]
  image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:ncol(Z.this),Ves.BC,las=1,cex.axis=.9)
  legend("top",paste("vessels with >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
  
  Drop.ves=All.ves[which(!All.ves%in%Ves.BC)]
  
  #Ves.sel.sens
  WHICh=which(Yrs.with.ktch>Ves.sel.sens)
  Z.this=ZZ[,WHICh]
  Ves.Sens=Vess[WHICh]
  ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
  Z.this=Z.this[,ID.sort]
  Ves.Sens=Ves.Sens[ID.sort]
  image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:ncol(Z.this),Ves.Sens,las=1,cex.axis=.6)
  legend("top",paste("vessels with >=",Ves.sel.sens, "years of records"),bty='n')
  
  #plot CPUEs
  a$CPUE.km.gn.day=a$LIVEWT.c/a$Km.Gillnet.Days.c
  a$CPUE.km.gn.h=a$LIVEWT.c/a$Km.Gillnet.Hours.c
  
  a.mean.cpue.km.day_all=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,SPECIES==SP),mean)
  a.mean.cpue.km.h_all=aggregate(CPUE.km.gn.h~YEAR.c,subset(a, SPECIES==SP),mean)
  
  a.mean.cpue.km.day_Sens=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
  a.mean.cpue.km.h_Sens=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
  
  a.mean.cpue.km.day_BC=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES==SP),mean)
  a.mean.cpue.km.h_BC=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES==SP),mean)
  
  #kmgday
  par(mar=c(3,3,.5,4),mgp=c(2,.7,0))
  Yrs=a.mean.cpue.km.day_Sens$YEAR.c
  plot(Yrs,a.mean.cpue.km.day_all$CPUE.km.gn.day,ylab="",xlab="")
  points(Yrs,a.mean.cpue.km.day_Sens$CPUE.km.gn.day,pch=19,col=2)
  if(nrow(a.mean.cpue.km.day_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.day_BC$CPUE.km.gn.day,pch=19,col=3)
  
  #kmgday  
  par(new = T)
  plot(Yrs,a.mean.cpue.km.h_all$CPUE.km.gn.h,ylab=NA, axes=F,xlab=NA,pch=0,cex=2)
  points(Yrs,a.mean.cpue.km.h_Sens$CPUE.km.gn.h,pch=15,col=2,cex=2)
  if(nrow(a.mean.cpue.km.h_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.h_BC$CPUE.km.gn.h,pch=15,col=3,cex=2)
  axis(side = 4)
  mtext("Financial year",1,line=2)
  mtext("Nominal CPUE (Kg/km.gn.day)",2,line=2)
  mtext("Nominal CPUE (Kg/km.gn.hour)",4,line=2)
  legend("top",c("Kg/km.gn.day","Kg/km.gn.hour"),bty='n',pch=c(0,19))
  legend("bottomleft",c("all",paste(Ves.sel.sens,"y"),paste(Ves.sel.BC,"y")),bty='n',pch=c(0,19,19),col=c(1,2,3))
  
  
  #For each vessel group, plot number of blocks by year
  #Ves.BC
  #all blocks
  AA=aggregate(LIVEWT.c~YEAR.c+BLOCKX,subset(a,SPECIES==SP & VESSEL%in%Ves.BC),sum)
  AA=reshape(AA,v.names="LIVEWT.c",idvar="YEAR.c",timevar="BLOCKX",direction="wide")
  AA=AA[order(AA$YEAR.c),]
  BLOCs=substr(names(AA)[2:ncol(AA)],10,30)
  Yrs=AA$YEAR.c
  Z=as.matrix(AA[,-1])
  ZZ=Z
  ZZ[ZZ>0]=1
  ZZZ=ZZ
  ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
  ZZZ=ZZZ[,ID.sort]
  par(mar=c(3,3.5,.8,.8))
  image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
  legend("top",paste("all blocks for vessels >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
  
  #blocks with > bc records
  Yrs.with.ktch=colSums(ZZ,na.rm=T)
  WHICh=which(Yrs.with.ktch>BLK.sel.BC)
  Z.this=ZZ[,WHICh]
  Blks.BC=BLOCs[WHICh]
  ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
  Z.this=Z.this[,ID.sort]
  Blks.BC=Blks.BC[ID.sort]
  image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:ncol(Z.this),Blks.BC,las=1,cex.axis=.9)
  legend("top",paste("blocks with >=",BLK.sel.BC, "years of records for vessels >=",Ves.sel.BC,"years of records and >",Min.ktch,"kg per year"),
         cex=0.75,bty='n')
  
  Drop.blks=All.blk[which(!All.blk%in%Blks.BC)]  
  
  #Ves.Sens
  #all blocks
  AA=aggregate(LIVEWT.c~YEAR.c+BLOCKX,subset(a,SPECIES==SP & VESSEL%in%Ves.Sens),sum)
  AA=reshape(AA,v.names="LIVEWT.c",idvar="YEAR.c",timevar="BLOCKX",direction="wide")
  AA=AA[order(AA$YEAR.c),]
  BLOCs=substr(names(AA)[2:ncol(AA)],10,30)
  Yrs=AA$YEAR.c
  Z=as.matrix(AA[,-1])
  ZZ=Z
  ZZ[ZZ>0]=1
  
  ZZZ=ZZ
  ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
  ZZZ=ZZZ[,ID.sort]
  
  par(mar=c(3,3.5,.8,.8))
  image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
  legend("top",paste("all blocks for vessels >=",Ves.sel.sens, "years of records"),bty='n')
  #blocks with > Sens records
  Yrs.with.ktch=colSums(ZZ,na.rm=T)
  WHICh=which(Yrs.with.ktch>BLK.sel.sens)
  Z.this=ZZ[,WHICh]
  Blks.Sens=BLOCs[WHICh]
  
  ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
  Z.this=Z.this[,ID.sort]
  Blks.Sens=Blks.Sens[ID.sort]
  
  
  image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
  axis(1,1:length(Yrs),Yrs)
  axis(2,1:ncol(Z.this),Blks.Sens,las=1,cex.axis=.9)
  legend("top",paste("blocks with >=",BLK.sel.sens, "years of records for vessels >=",Ves.sel.sens,"years of records of",SP),bty='n')
  dev.off()
  
  return(list(Ves.BC=Ves.BC, Ves.Sens=Ves.Sens, Blks.BC=Blks.BC, Blks.Sens=Blks.Sens, Drop.ves=Drop.ves, Drop.blks=Drop.blks))
}

#functions for reshaping data
  #monthly data
Effort.data.fun=function(DATA,target,ktch)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  #remove nonsense lat
  DATA=subset(DATA,LAT>=(-36))
  
  #calculate effort
  Match.these.eff=match(These.efforts,names(DATA))
  Effort.data1=DATA[,Match.these.eff]
  Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c,Km.Gillnet.Hours_shot.c)~zone+
                          FINYEAR+Same.return+MONTH+BLOCKX+SHOTS.c+BDAYS.c+HOURS.c+NETLEN.c,Effort.data1,max)
  # Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+      #redundant
  #                         FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,sum)
  
  Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data1,max)
  #Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data.inv,sum)

  Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return",all.x=T)

  #target species catch 
  ID=match(c(ktch),colnames(DATA))
  DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
  DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
  #note: catch targeted at other species: pointless as these are multiple trips combined in one month
  
  #reshape catch data
  TABLE=aggregate(cbind(Catch.Target,Catch.Total)~MONTH+FINYEAR+BLOCKX+VESSEL+Same.return+LAT+LONG+YEAR.c,data=DATA,sum,na.rm=T)
  TABLE=TABLE[order(TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  
  #proportion of records with target catch
  prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
  
  #merge catch and effort
  dat=merge(TABLE,Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"),all.x=T)
  
  
  
  #Add mesh size                                
  d=subset(DATA,select=c(Same.return,mesh))
  d=d[!duplicated(d$Same.return),]
  dat=merge(dat,d,by="Same.return",all.x=T)
  
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}

  #daily data
Effort.data.fun.daily=function(DATA,target,ktch,Aggregtn)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  #remove nonsense lat
  DATA=subset(DATA,LAT>=(-36))
  
  #calculate effort
  Match.these.eff=match(These.efforts.daily,names(DATA))
  Effort.data1=DATA[,Match.these.eff]
  
  #aggregate at shot level
  if(Use.Date=="NO")
  {
    #max effort by Sno, DsNo & TSNo
    Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c,Km.Gillnet.Hours_shot.c)~zone+
                            FINYEAR+Same.return.SNo+MONTH+BLOCKX+block10+shots.c+hours.c+netlen.c,Effort.data1,max)
    Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return.SNo,Effort.data1,max)
 
    
    #aggregate at TSNo if required
    if(Aggregtn=="TSNo")
    {
      Effort.data$TSNo=word(Effort.data$Same.return.SNo,3)
      Effort.data.inv$TSNo=word(Effort.data.inv$Same.return.SNo,3)
      
      Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c,Km.Gillnet.Hours_shot.c)~zone+
                              FINYEAR+TSNo+MONTH+BLOCKX,Effort.data,sum)
      Effort.data.inv=aggregate(Km.Gillnet.Days.inv~TSNo,Effort.data.inv,sum)
    }
    
    #merge as appropriate
    if(Aggregtn=="SNo")
    {
      Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return.SNo",all.x=T)
    }
    if(Aggregtn=="TSNo")
    {
      Effort.data=merge(Effort.data,Effort.data.inv,by="TSNo",all.x=T)
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
  DATA$Catch.Dhufish=with(DATA,ifelse(SPECIES%in%c(320000),DATA[,ID],0))
  
  Other.shk=c(5001:24900,25000:31000)
  Other.shk=subset(Other.shk,!Other.shk%in%c(17001,17003,18003,18001,18007))
  DATA$Catch.Other.shrk=with(DATA,ifelse(SPECIES%in%Other.shk,DATA[,ID],0))
  
  Other.scalie=188000:599001
  Other.scalie=subset(Other.scalie,!Other.scalie%in%c(384002,353001,377004,320000))
  DATA$Catch.Other.scalefish=with(DATA,ifelse(SPECIES%in%Other.scalie,DATA[,ID],0))
  
  DATA$Catch.non_indicators=with(DATA,ifelse(!SPECIES%in%c(17001,17003,18003,18001,18007),DATA[,ID],0))
  
  DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
  
  #reshape catch data
  if(Use.Date=="NO")
  {
    if(Aggregtn=="SNo") TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Dusky,Catch.Sandbar,
                                              Catch.Groper,Catch.Snapper,Catch.Blue_mor,Catch.Dhufish,
                                              Catch.Other.shrk,Catch.Other.scalefish,
                                              Catch.non_indicators,Catch.Total)~MONTH+FINYEAR+BLOCKX+block10+VESSEL+
                                          Same.return.SNo+date+LAT+LONG+YEAR.c,data=DATA,sum,na.rm=T)   
    
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
    }
  }
  
  if(Aggregtn=="SNo")TABLE=TABLE[order(TABLE$Same.return.SNo,TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  if(Aggregtn=="TSNo")TABLE=TABLE[order(TABLE$TSNo,TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  
  #proportion of records with target catch
  prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
  
  #merge catch and effort
  if(Aggregtn=="SNo")dat=merge(TABLE,Effort.data,by=c("Same.return.SNo","FINYEAR","MONTH","BLOCKX","block10"),all.x=T)
  if(Aggregtn=="TSNo")dat=merge(TABLE,Effort.data,by=c("TSNo","FINYEAR","MONTH","BLOCKX"),all.x=T)
  
  
  #Add mesh size, shots, depth and nlines for each session
  if(Aggregtn=="SNo")
  {
    d=subset(DATA,select=c(Same.return.SNo,VESSEL,mesh,nlines.c,Mean.depth))
    d=d[!duplicated(paste(d$Same.return.SNo,d$VESSEL)),]
    dat=merge(dat,d,by=c("Same.return.SNo","VESSEL"),all.x=T)
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

fn.scale=function(x,max,scaler) ((x/max)^0.5)*scaler

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
  
  
  ##==========================================================================
  ## calculate year-specific QL and add column to data set
  ##==========================================================================
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
  dat$CX=mapply(fn.scale,dat$cpue,max(dat$cpue),scaler)
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
         pt.cex=mapply(fn.scale,Qnt,max(Qnt),scaler),cex=1.5,title="kg/km gn d")
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

#calculate 4 different nominal cpues
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
    
    ##==========================================================================
    ## calculate year-specific QL and add column to data set
    ##==========================================================================
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
    
    ##==========================================================================
    ## malcolm's targeting - proportion below which no variation exists in cpue
    ##==========================================================================
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
    
    
    ##==========================================================================
    ## plot raw mean cpues and CIs using 4 different data sets
    ##==========================================================================
    
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
    
    CPUE.Malcolm[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, prop>Prop.Malcolm), catch.column="catch",
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
    polygon(x=c(Nx.p,rev(Nx.p)),y=c(rep(show.pol[1]-1,2),rep(show.pol[length(show.pol)],2)),col=rgb(.1,.1,.1,.25),border='transparent')
    for(i in 1:nrow(dd)) points(Nx[1]:Nx[2],rep(i,Nx[2]),pch=21,cex=fn.scale(dd[i,],Mx,2.5),bg=rgb(.1,.1,.1,.4))
    axis(1,1:ncol(dd),colnames(dd))
    axis(2,1:nrow(dd),rownames(dd),las=1,cex.axis=.5)
    mtext("number of records per year",3,1,cex=1.25)
    Lab=round(quantile(dd,probs=c(.95,.995,1)))
    legend('topright',paste(Lab),pch=21,pt.bg="grey70",
           pt.cex=fn.scale(Lab,Mx,2.5),horiz=T,title="# of records")
    
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
      B= d[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
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

fn.expl.cede=function(d,PREDS,kg,Do.ggplts)    #function for exploratory analysis
{
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
  if(!is.na(match("DepCat",PREDS)))
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
    #ggplot(dd, aes(x = month, y = cpue, color = cluster_clara)) +   geom_boxplot()
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
fn.modl.sel=function(Inter,RESPNS)   #function for model structure selection
{
  PREDS[id.cov]=paste("LN",PREDS[id.cov],sep="")
  
  Pos.mod=2^length(PREDS)
  if(RESPNS=="LNcpue")Formula=formula(paste("LNcpue",paste(PREDS,collapse="+"),sep="~"))
  
  if(RESPNS=="catch")Formula=formula(paste(Response,paste(paste(PREDS,collapse="+"),
                                                          paste("offset(","LN",efrt,")",sep=""),sep="+"),sep="~"))
  
  Level=ifelse(Inter=="MainTerm",1,ifelse(Inter=="2way",2,"3way"))
  MeThod="h"
  
  res <- glmulti(Formula,data=d, level=Level,method=MeThod, fitfunction=fitFun,
                 crit="aicc", confsetsize=Pos.mod,plotty = F, report = T)
  return(list(res=res,BEST=res@formulas[[1]]))
}

viz.coef=function(MOD,WHAT) coefplot(MOD,coefficient=WHAT)

fn.stand=function(d,Response,RESPNS,PREDS,efrt,Formula)   #function for standardisation
{
  id.fctr=which(PREDS%in%Categorical)
  d=makecategorical(PREDS[id.fctr],d) 
  d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
  id.cov=which(PREDS%in%Covariates)
  d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
  PREDS[id.cov]=paste("LN",PREDS[id.cov],sep="")
  if(RESPNS=="LNcpue") res <- glm(Formula,data=d)
  if(RESPNS=="catch")  res <- glm.nb(Formula,data=d)
  return(list(res=res,DATA=d))
}

Anova.and.Dev.exp=function(GLM,SP,type)   #function for extracting term significance and deviance explained
{
  #Anovas
  Anova.tab=anova(GLM, test = "Chisq")
  
  #Deviance explained
  #By each term
  n=2:length(Anova.tab$Deviance)
  Term.dev.exp=100*(Anova.tab$Deviance[n]/GLM$null.deviance)
  names(Term.dev.exp)=rownames(Anova.tab)[n]
  
  #By full model
  Dev.exp=sum(Term.dev.exp)
  
  #Combine as table
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
  Table=rbind(Table,All)
  
  Table$"Pr(>Chi)"=ifelse(Table$"Pr(>Chi)"<0.001,"<0.001",Table$"Pr(>Chi)")
  dummy=Table[1:2,]
  dummy[,]=NA
  rownames(dummy)=c(type,SP)
  Table=rbind(dummy,Table)
  Table[is.na(Table)] <- ""
  Table=cbind(data.frame(Term=rownames(Table)),Table)
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
  if(biascor=="YES")
  {
    sigma.glm=sqrt(summary(MOD)$dispersion) # residuals standard error
    lsm$response=exp(lsm$emmean)*exp(sigma.glm^2/2)
    lsm$lower.CL=exp(lsm$asymp.LCL)*exp(sigma.glm^2/2)
    lsm$upper.CL=exp(lsm$asymp.UCL)*exp(sigma.glm^2/2)
  }
  if(biascor=="NO")
  {
    if(is.na(match("response",names(lsm))))lsm$response=lsm$emmean
    lsm$lower.CL=lsm$asymp.LCL
    lsm$upper.CL=lsm$asymp.UCL
  }
  return(lsm)
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


Plot.cpue=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar)    #plot cpues
{
  if(inherits(cpuedata, "list")) 
  {
    if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
    if(length(cpuedata)<3)tc=seq(-.5*0.15,.5*0.15,length.out=length(cpuedata))
    ymax = max(unlist(lapply(cpuedata, `[`, "upper.CL")))
    Yrs=as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4))
    plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
    if(COL=='color')CL=c("black","forestgreen","dodgerblue","red","bisque3")
    if(COL=='grey') CL=gray.colors(length(cpuedata),start=0.2,end=0.65)
    for(l in 1:length(cpuedata))
    {
      with(cpuedata[[l]],
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
    ymax = max(cpuedata$upper.CL)
    Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
    plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
    CL="black"
    points(Yrs, cpuedata$response, "o", pch=16, lty=2, col=CL,cex=CxS)
    arrows(x0=Yrs, y0=cpuedata$lower.CL, 
           x1=Yrs, y1=cpuedata$upper.CL, 
           code=3, angle=90, length=0.05, col=CL)
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

Plot.cpue.spatial=function(cpuedata,scaler,colPalet,CxTxt)
{
  cpuedata$CV=round(100*((cpuedata$upper.CL-cpuedata$response)/1.96)/cpuedata$response)
  rbPal <- colorRampPalette(colPalet)
  Nbreaks=10
  cpuedata$Col <- rbPal(Nbreaks)[as.numeric(cut(cpuedata$CV,breaks = Nbreaks))]
  
  with(cpuedata,plot(LONG+0.5,LAT-0.5,ylim=c(-36,-26),xlim=c(113,129),
                     cex=((response/max(response))^0.5)*scaler,pch=19,col=Col,
                     yaxt='n',xaxt='n',ylab="",xlab=""))
  axis(1,113:129,F)
  axis(2,-36:-26,F)
  LEgnd=round(quantile(cpuedata$response,probs=c(.01,.5,.99)),2)
  LEgnd.cex=((Quant/max(cpuedata$response))^0.5)*scaler
  legend('top',paste(LEgnd),bty='n',title="CPUE (kg/km gillnet h)",pch=19,pt.cex=LEgnd.cex,cex=1.25)
  
  LEgn.err=round(quantile(cpuedata$CV,probs=c(.01,.5,.85,.99)))
  legend("topright",paste(LEgn.err,"%"),bty='n',title="CV",
         col =rbPal(Nbreaks)[as.numeric(cut(LEgn.err,breaks = Nbreaks))],pch=19,pt.cex=2)
  with(cpuedata,text(LONG+0.5,LAT-0.5,blockx,cex=CxTxt,srt=45,adj=c(-0.25,0.5)))
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
  
  plot(PRED,Std.RES,ylab="",xlab="")
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

#----4. PROCEDURE SECTION-----#


#Keep vessel characteristics from vessel survey for vessels that have fished      #missing: fill in missing vessel info and add to cpue data. Also add wheather prediction
if(exists("TDGDLF.survey"))TDGDLF.survey=subset(TDGDLF.survey, BOATREGO%in%unique(
  c(Data.monthly.GN.whiskery$VESSEL,Data.monthly.GN.gummy$VESSEL,
    Data.monthly.GN.dusky$VESSEL,Data.monthly.GN.sandbar$VESSEL,
    Data.daily.GN.whiskery$VESSEL,Data.daily.GN.gummy$VESSEL,
    Data.daily.GN.dusky$VESSEL,Data.daily.GN.sandbar$VESSEL)
))    


#4.1 Deal with zone1-zone2 Boundary blocks to a zone 
Boundary.Blks=c(34160,35160,36160)
if(BOUND.BLK=="REMOVE")
{
  Data.monthly.GN.whiskery=subset(Data.monthly.GN.whiskery,!BLOCKX%in%Boundary.Blks)
  Data.monthly.GN.gummy=subset(Data.monthly.GN.gummy,!BLOCKX%in%Boundary.Blks)
  Data.monthly.GN.dusky=subset(Data.monthly.GN.dusky,!BLOCKX%in%Boundary.Blks)
  Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,!BLOCKX%in%Boundary.Blks) 
}

if(BOUND.BLK=="REALLOCATE")
{
  #a=subset(Data.monthly.GN.whiskery,BLOCKX%in%Boundary.Blks)
  #b=subset(Data.monthly.GN.whiskery,!BLOCKX%in%Boundary.Blks)
  # RAND=sample(c("Zone2","Zone1"),nrow(a),replace=T)
  # a$zone=with(a,ifelse(BLOCKX%in%Boundary.Blks,RAND,zone))
  #a$zone="Zone1"
  #Data.monthly.GN.whiskery=rbind(b,a)
  Data.monthly.GN.whiskery$zone=with(Data.monthly.GN.whiskery,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))
  
  #a=subset(Data.monthly.GN.gummy,BLOCKX%in%Boundary.Blks)
  #b=subset(Data.monthly.GN.gummy,!BLOCKX%in%Boundary.Blks)
  # RAND=sample(c("Zone2","Zone1"),nrow(a),replace=T)
  # a$zone=with(a,ifelse(BLOCKX%in%Boundary.Blks,RAND,zone))
  #if(nrow(a)>0)a$zone="Zone1"
  #Data.monthly.GN.gummy=rbind(b,a)
  Data.monthly.GN.gummy$zone=with(Data.monthly.GN.gummy,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))
  
  
  #a=subset(Data.monthly.GN.dusky,BLOCKX%in%Boundary.Blks)
  #b=subset(Data.monthly.GN.dusky,!BLOCKX%in%Boundary.Blks)
  # RAND=sample(c("Zone2","Zone1"),nrow(a),replace=T)
  # a$zone=with(a,ifelse(BLOCKX%in%Boundary.Blks,RAND,zone))
  #a$zone="Zone1"
  #Data.monthly.GN.dusky=rbind(b,a) 
  Data.monthly.GN.dusky$zone=with(Data.monthly.GN.dusky,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))
  
  #a=subset(Data.monthly.GN.sandbar,BLOCKX%in%Boundary.Blks)
  #b=subset(Data.monthly.GN.sandbar,!BLOCKX%in%Boundary.Blks)
  # RAND=sample(c("Zone2","Zone1"),nrow(a),replace=T)
  # a$zone=with(a,ifelse(BLOCKX%in%Boundary.Blks,RAND,zone))
  #a$zone="Zone1"
  #Data.monthly.GN.sandbar=rbind(b,a) 
  Data.monthly.GN.sandbar$zone=with(Data.monthly.GN.sandbar,ifelse(BLOCKX%in%Boundary.Blks,"Zone1",zone))
  
}


#4.2 Extract number of vessels per species range
TARGETS.name=c("SHARK, WHISKERY","SHARK, GUMMY","SHARK, BRONZE WHALER","SHARK, THICKSKIN (SANDBAR)")
TARGETS=list(17003,17001,c(18003,18001),18007)
names(TARGETS)=TARGETS.name
N.species=length(TARGETS)

N.VES=matrix(rep(0,N.species),ncol=N.species)
colnames(N.VES)=c("WH","GM","Dus","San")

N.VES[,1]=fn.n.vessels(subset(Data.monthly.GN.whiskery,SPECIES==17003),
                       subset(Data.daily.GN.whiskery,SPECIES==17003))

N.VES[,2]=fn.n.vessels(subset(Data.monthly.GN.gummy,SPECIES==17001),
                       subset(Data.daily.GN.gummy,SPECIES==17001))

N.VES[,3]=fn.n.vessels(subset(Data.monthly.GN.dusky,SPECIES%in%c(18003,18001)),
                       subset(Data.daily.GN.dusky,SPECIES%in%c(18003,18001)))

N.VES[,4]=fn.n.vessels(subset(Data.monthly.GN.sandbar,SPECIES==18007 & YEAR.c>=1985),
                       subset(Data.daily.GN.sandbar,SPECIES==18007 & YEAR.c>=1985))
setwd('C:/Matias/Analyses/Catch and effort')
hndl=paste(getwd(),"/Outputs/Paper/",sep="")
write.csv(N.VES,paste(hndl,"All.Vessels.by.species.csv",sep=""),row.names=F)


#4.2.1 Extract number of blocks where shark has been caught within effective area
Tol.blks.whi=length(fn.n.blk.caught(Data.monthly.GN.whiskery,Data.daily.GN.whiskery,17003))
Tol.blks.gum=length(fn.n.blk.caught(Data.monthly.GN.gummy,Data.daily.GN.gummy,17001))
Tol.blks.dus=length(fn.n.blk.caught(Data.monthly.GN.dusky,Data.daily.GN.dusky,c(18003,18001)))
Tol.blks.san=length(fn.n.blk.caught(subset(Data.monthly.GN.sandbar,YEAR.c>=1986),Data.daily.GN.sandbar,18007))

Blks.by.species=cbind(Tol.blks.whi,Tol.blks.gum,Tol.blks.dus,Tol.blks.san)
colnames(Blks.by.species)=c("WH","GM","Dus","San")
write.csv(Blks.by.species,paste(hndl,"All.Blks.by.species.csv",sep=""),row.names=F)


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


#4.3 Data fixes
#Remove initial years for sandbar as they were not reported in fishery stats
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,
                               !(FINYEAR%in%paste(1975:1984,substr(1976:1985,3,4),sep="-")))


#4.4 Remove NA effort
#monthly
Data.monthly.GN.whiskery=subset(Data.monthly.GN.whiskery,!is.na(Km.Gillnet.Days.c))
Data.monthly.GN.gummy=subset(Data.monthly.GN.gummy,!is.na(Km.Gillnet.Days.c))
Data.monthly.GN.dusky=subset(Data.monthly.GN.dusky,!is.na(Km.Gillnet.Days.c))
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,!is.na(Km.Gillnet.Days.c))

#daily
Data.daily.GN.whiskery=subset(Data.daily.GN.whiskery,!is.na(Km.Gillnet.Days.c))
Data.daily.GN.gummy=subset(Data.daily.GN.gummy,!is.na(Km.Gillnet.Days.c))
Data.daily.GN.dusky=subset(Data.daily.GN.dusky,!is.na(Km.Gillnet.Days.c))
Data.daily.GN.sandbar=subset(Data.daily.GN.sandbar,!is.na(Km.Gillnet.Days.c))


#4.5 Create useful vars
FINYEAR.monthly=as.character(unique(Data.monthly.GN.whiskery$FINYEAR))
FINYEAR.monthly=sort(FINYEAR.monthly)
N.yrs=length(FINYEAR.monthly)

FINYEAR.daily=as.character(unique(Data.daily.GN.whiskery$FINYEAR))
FINYEAR.daily=sort(FINYEAR.daily)
N.yrs.daily=length(FINYEAR.daily) 

FINYEAR.ALL=c(FINYEAR.monthly,FINYEAR.daily)
FINYEAR.ALL=sort(FINYEAR.ALL)
N.yrs.ALL=length(FINYEAR.ALL)

new.q=FINYEAR.monthly[match(Q_change,FINYEAR.monthly):length(FINYEAR.monthly)]  #whiskery q years


#4.6 Proportion of dusky and copper shark
fn.fig("proportion of dusky and copper shark_TDGLDF",2000,2400)
par(mfcol=c(2,1),las=1,mai=c(.8,.85,.1,.1),mgp=c(2.5,.8,0))
  #monthly
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Data.monthly.GN.dusky,SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Data.monthly.GN.dusky,SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="",xlab="",pch=19,col=2,cex=1.75,cex.lab=1.5,ylim=c(0,.25))
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)
legend("topright","Monthly returns",bty='n',cex=1.5)

  #daily
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Data.daily.GN.dusky,SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Data.daily.GN.dusky,SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="",xlab="Financial year",pch=19,col=2,cex=1.75,cex.lab=1.5,ylim=c(0,.25))
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)
legend("topright","Daily logbooks",bty='n',cex=1.5)

mtext("Bronze whaler shark catch / Dusky shark catch",2,-1.5,las=3,outer=T,cex=1.75)
dev.off()


#4.7 Define indicative vessels and blocks 
#steps: 1. select vessels that meet criteria (fishing for at least Threshold.n.yrs/Threshold.n.yrs.daily
#         and catching at least MIN.ktch)
#       2. for those vessels, select blocks with at least MIN.obs.BLK years of observations
BLKS.used=vector('list',length=N.species)
names(BLKS.used)=SPECIES.vec 
VES.used=BLKS.used.daily=VES.used.daily=BLKS.used

if(Remove.blk.by=="blk_only")  
{
  #Monthly
  Sort.blks.gum=fn.see.all.yrs.ves.blks(a=Data.monthly.GN.gummy,SP=17001,what=".monthly",
                                        Ves.sel.BC=Threshold.n.yrs.monthly,Ves.sel.sens=Threshold.n.yrs.sens,
                                        BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.whi_q.change=fn.see.all.yrs.ves.blks(a=Data.monthly.GN.whiskery,SP=17003,what=".monthly",
                                                 Ves.sel.BC=Threshold.n.yrs.monthly,Ves.sel.sens=Threshold.n.yrs.sens,
                                                 BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.dus=fn.see.all.yrs.ves.blks(a=Data.monthly.GN.dusky,SP=18003,what=".monthly",
                                        Ves.sel.BC=Threshold.n.yrs.monthly,Ves.sel.sens=Threshold.n.yrs.sens,
                                        BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.san=fn.see.all.yrs.ves.blks(a=subset(Data.monthly.GN.sandbar,!FINYEAR%in%c("1985-86","1986-87","1987-88")),SP=18007,what=".monthly",
                                        Ves.sel.BC=Threshold.n.yrs.monthly,Ves.sel.sens=Threshold.n.yrs.sens,
                                        BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  
  
  #identified selected and dropped vessels and blocks     
  BLKS.used$`Gummy shark`=blok.gum=Sort.blks.gum$Blks.BC
  No.blok.gum=Sort.blks.gum$Drop.blks
  VES.used$`Gummy shark`=vsl.gum=Sort.blks.gum$Ves.BC
  No.vsl.gum=Sort.blks.gum$Drop.ves
  
  BLKS.used$`Whiskery shark`=blok.whi_q.change=Sort.blks.whi_q.change$Blks.BC
  No.blok.whi_q.change=Sort.blks.whi_q.change$Drop.blks
  VES.used$`Whiskery shark`=vsl.whi=Sort.blks.whi_q.change$Ves.BC
  No.vsl.whi=Sort.blks.whi_q.change$Drop.ves
  
  BLKS.used$`Dusky shark`=blok.dus=Sort.blks.dus$Blks.BC
  No.blok.dus=Sort.blks.dus$Drop.blks
  VES.used$`Dusky shark`=vsl.dus=Sort.blks.dus$Ves.BC
  No.vsl.dus=Sort.blks.dus$Drop.ves
  
  BLKS.used$`Sandbar shark`=blok.san=Sort.blks.san$Blks.BC
  No.blok.san=Sort.blks.san$Drop.blks
  VES.used$`Sandbar shark`=vsl.san=Sort.blks.san$Ves.BC
  No.vsl.san=Sort.blks.san$Drop.ves
  
  
  #Daily
  Sort.blks.gum.daily=fn.see.all.yrs.ves.blks(a=Data.daily.GN.gummy,SP=17001,what=".daily",
                                              Ves.sel.BC=Threshold.n.yrs.daily,Ves.sel.sens=Threshold.n.yrs.sens,
                                              BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.whi.daily=fn.see.all.yrs.ves.blks(a=Data.daily.GN.whiskery,SP=17003,what=".daily",
                                              Ves.sel.BC=Threshold.n.yrs.daily,Ves.sel.sens=Threshold.n.yrs.sens,
                                              BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.dus.daily=fn.see.all.yrs.ves.blks(a=Data.daily.GN.dusky,SP=18003,what=".daily",
                                              Ves.sel.BC=Threshold.n.yrs.daily,Ves.sel.sens=Threshold.n.yrs.sens,
                                              BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  Sort.blks.san.daily=fn.see.all.yrs.ves.blks(a=Data.daily.GN.sandbar,SP=18007,what=".daily",
                                              Ves.sel.BC=Threshold.n.yrs.daily,Ves.sel.sens=Threshold.n.yrs.sens,
                                              BLK.sel.BC=MIN.obs.BLK,BLK.sel.sens=MIN.obs.BLK.sens,Min.ktch=MIN.ktch)
  
  #identified selected and dropped vessels and blocks 
  BLKS.used.daily$`Gummy shark`=blok.gum.daily=Sort.blks.gum.daily$Blks.BC
  No.blok.gum.daily=Sort.blks.gum.daily$Drop.blks
  VES.used.daily$`Gummy shark`=vsl.gum.daily=Sort.blks.gum.daily$Ves.BC
  No.vsl.gum.daily=Sort.blks.gum.daily$Drop.ves
  
  BLKS.used.daily$`Whiskery shark`=blok.whi.daily=Sort.blks.whi.daily$Blks.BC
  No.blok.whi.daily=Sort.blks.whi.daily$Drop.blks
  VES.used.daily$`Whiskery shark`=vsl.whi.daily=Sort.blks.whi.daily$Ves.BC
  No.vsl.whi.daily=Sort.blks.whi.daily$Drop.ves
  
  BLKS.used.daily$`Dusky shark`=blok.dus.daily=Sort.blks.dus.daily$Blks.BC
  No.blok.dus.daily=Sort.blks.dus.daily$Drop.blks
  VES.used.daily$`Dusky shark`=vsl.dus.daily=Sort.blks.dus.daily$Ves.BC
  No.vsl.dus.daily=Sort.blks.dus.daily$Drop.ves
  
  BLKS.used.daily$`Sandbar shark`=blok.san.daily=Sort.blks.san.daily$Blks.BC
  No.blok.san.daily=Sort.blks.san.daily$Drop.blks
  VES.used.daily$`Sandbar shark`=vsl.san.daily=Sort.blks.san.daily$Ves.BC
  No.vsl.san.daily=Sort.blks.san.daily$Drop.ves
}


#4.8 Put data into a list
#note: this still has all records, all vessels and all blocks
Species.list=list(whiskery=Data.monthly.GN.whiskery,gummy=Data.monthly.GN.gummy,
                  dusky=Data.monthly.GN.dusky,sandbar=Data.monthly.GN.sandbar)

Species.list.daily=list(whiskery=Data.daily.GN.whiskery,gummy=Data.daily.GN.gummy,
                        dusky=Data.daily.GN.dusky,sandbar=Data.daily.GN.sandbar)


#reset blockxs to first 4 digits as some diferent 5 digit blockxs have same lat and long
for(s in 1:N.species)
{
  Species.list[[s]]$BLOCKX=as.integer(substr(Species.list[[s]]$BLOCKX,1,4))
  Species.list[[s]]$Same.return=with(Species.list[[s]],paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
  
  Species.list.daily[[s]]$BLOCKX=as.integer(substr(Species.list.daily[[s]]$BLOCKX,1,4))
  Species.list.daily[[s]]$Same.return=with(Species.list.daily[[s]],paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
  
}

rm(Data.monthly.GN.whiskery,Data.monthly.GN.gummy,Data.monthly.GN.dusky,Data.monthly.GN.sandbar,
   Data.daily.GN.whiskery,Data.daily.GN.gummy,Data.daily.GN.dusky,Data.daily.GN.sandbar)


#Create Km.Gillnet.Hours_shot.c
for(s in 1:N.species)
{
  Species.list[[s]]$Km.Gillnet.Hours_shot.c=with(Species.list[[s]],Km.Gillnet.Hours.c*SHOTS.c)
  Species.list.daily[[s]]$Km.Gillnet.Hours_shot.c=with(Species.list.daily[[s]],Km.Gillnet.Hours.c*shots.c)
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


#4.9 Compare nominal all records VS 'good reporters' only
for(i in 1:N.species)
{
  fn.fig(paste(hndl,"All records vs Good records/",names(Species.list)[i],sep=""),2400,2400)
  par(mfcol=c(2,1),mai=c(.5,1.1,.3,.1),mgp=c(2.5,.5,0))
  #Monthly
  test=subset(Species.list[[i]],SPECIES%in%TARGETS[[i]])
  D=subset(test,Reporter=="good")
  test=subset(test,FINYEAR%in%unique(D$FINYEAR))
  agg.good=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~FINYEAR,D,mean)
  agg.all=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~FINYEAR,test,mean)
  S=c(agg.all[,2],agg.good[,2])
  plot(agg.all[,2],pch=19,col=2,type='l',lwd=2,ylim=c(min(S),max(c(S))),ylab="CPUE")
  lines(agg.good[,2],lwd=2,col=3)
  legend("top",c("all records","good records"),bty='n',col=2:3,lty=1,lwd=2)
  
  #Daily
  test=subset(Species.list.daily[[i]],SPECIES%in%TARGETS[[i]])
  D=subset(test,Reporter=="good")
  agg.good=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~FINYEAR,D,mean)
  agg.all=aggregate((LIVEWT.c/Km.Gillnet.Days.c)~FINYEAR,test,mean)
  S=c(agg.all[,2],agg.good[,2])
  plot(agg.all[,2],pch=19,col=2,type='l',lwd=2,ylim=c(min(S),max(c(S))),ylab="CPUE")
  lines(agg.good[,2],lwd=2,col=3)
  
  dev.off()
}


#4.10 Construct wide database for analysis 
#steps: 
#   1. select "Good" records (the variable "Reporter" includes good/bad catch and effort reporters)
#   2. Construct a single record for each record (i.e. year-month-vessel-block-gear for monthly
#      returns and year-Session-vessel-block10-gear for daily logbooks), with catch of target
#      and other species as separate columns, giving a 0 catch for column "target" if no catch
#      where a Session is TSNo or SNo.DSNo.TSNo (see below)

#note: this still keeps all blocks and vessels from good reporters

#variables used in standardisation
#monthly
These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Days.inv",
                "Km.Gillnet.Hours.c","Km.Gillnet.Days.c","Km.Gillnet.Hours_shot.c",
                "zone","MONTH","BLOCKX","mesh","SHOTS.c","BDAYS.c","HOURS.c","NETLEN.c")

#daily
These.efforts.daily=c("FINYEAR","date","TSNo","Km.Gillnet.Days.inv",
                      "Km.Gillnet.Hours.c","Km.Gillnet.Days.c","Km.Gillnet.Hours_shot.c",
                      "zone","MONTH","BLOCKX","block10","VESSEL","mesh",
                      "Same.return.SNo","nlines.c","shots.c","netlen.c","hours.c")

  #monthly  
DATA.list.LIVEWT.c=vector('list',length=N.species)
Prop.Catch=vector(length=N.species)
names(DATA.list.LIVEWT.c)=names(Prop.Catch)=names(Species.list)
for ( i in 1:N.species)
{
  #create data sets 
  dummy=Effort.data.fun(subset(Species.list[[i]],Reporter=="good"),TARGETS[[i]],"LIVEWT.c")  
  DATA.list.LIVEWT.c[[i]]=dummy$dat   
  
  #proportion with catch
  Prop.Catch[i]=dummy$prop.with.catch   
}


  #daily 
#note: DATA.list.LIVEWT.c.daily has kg as response variable and it's aggregated by trip 
#         as this is how data collected so block10 is dropped as can be multiple block10s per trip
#     DATA.list.LIVEWT.c.daily_nfish has nfish per Session as as response variable as 
#         this is how data collected
DATA.list.LIVEWT.c.daily=vector('list',length=N.species)
Prop.Catch.daily=vector(length=N.species)
names(DATA.list.LIVEWT.c.daily)=names(Prop.Catch.daily)=names(Species.list.daily)
DATA.list.LIVEWT.c.daily_nfish=DATA.list.LIVEWT.c.daily
for ( i in 1:N.species)
{
  #set Reporter to Bad if average weight is outside species range
  #note: some records have dodgy nfish so class as 'bad reporters'
  Species.list.daily[[i]]$Avrg.w=Species.list.daily[[i]]$LIVEWT.c/Species.list.daily[[i]]$nfish
  Species.list.daily[[i]]$Reporter=with(Species.list.daily[[i]],
                                        ifelse((Avrg.w>Max.weight[i]| Avrg.w<Min.weight[i]) & SPECIES==SPvec[i],"bad",Reporter))                    
  
  
  
  #create data sets 
  #Catch in kg aggregated by TSNo (as catch is weighed at end of trip)
  dummy=Effort.data.fun.daily(subset(Species.list.daily[[i]],Reporter=="good"),TARGETS[[i]],
                              ktch="LIVEWT.c",Aggregtn="TSNo")
  DATA.list.LIVEWT.c.daily[[i]]=dummy$dat
  
  #proportion with catch
  Prop.Catch.daily[i]=dummy$prop.with.catch   
  
  
  #Catch in Numbers at SNo as number of fish is reported for each session  
  dummy=Effort.data.fun.daily(subset(Species.list.daily[[i]],Reporter=="good" & !is.na(nfish)),TARGETS[[i]],
                              ktch="nfish",Aggregtn="SNo")
  DATA.list.LIVEWT.c.daily_nfish[[i]]=dummy$dat
  
  rm(dummy)
  
}

write.csv(Prop.Catch,paste(hndl,"Prop.records.with.catch.monthly.csv",sep=""),row.names=T)
write.csv(Prop.Catch.daily,paste(hndl,"Prop.records.with.catch.daily.csv",sep=""),row.names=T)

#Remove NA effort
for ( i in 1:N.species)
{
  DATA.list.LIVEWT.c[[i]]=subset(DATA.list.LIVEWT.c[[i]],!is.na(Km.Gillnet.Days.c))
  DATA.list.LIVEWT.c.daily[[i]]=subset(DATA.list.LIVEWT.c.daily[[i]],!is.na(Km.Gillnet.Days.c))
}

#4.11 Add SOI, Freo and moon phase to SNo    
Freo$Freo.Lag6=c(rep(NA,6),Freo$Freo[1:(length(Freo$Freo)-6)])
Freo$Freo.Lag12=c(rep(NA,12),Freo$Freo[1:(length(Freo$Freo)-12)])
Freo=subset(Freo,Year>1973)
for(i in 1:N.species)
{
  #monthly
  a=SOI
  a$dummy=with(a,paste(Year.soi, Month.soi))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c[[i]]=merge(DATA.list.LIVEWT.c[[i]],a, by.x=c("YEAR.c","MONTH"),
                                by.y=c("Year.soi","Month.soi"),all.x=T)
  
  a=Freo
  a$dummy=with(a,paste(Year, Month))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c[[i]]=merge(DATA.list.LIVEWT.c[[i]],a, by.x=c("YEAR.c","MONTH"),
                                by.y=c("Year","Month"),all.x=T)
  
  #daily
  a=SOI
  a$dummy=with(a,paste(Year.soi, Month.soi))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily[[i]]=merge(DATA.list.LIVEWT.c.daily[[i]],a, by.x=c("YEAR.c","MONTH"),
                                      by.y=c("Year.soi","Month.soi"),all.x=T)
  
  a=Freo
  a$dummy=with(a,paste(Year, Month))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily[[i]]=merge(DATA.list.LIVEWT.c.daily[[i]],a, by.x=c("YEAR.c","MONTH"),
                                      by.y=c("Year","Month"),all.x=T)
  
  
  
  #daily nfish
  a=SOI
  a$dummy=with(a,paste(Year.soi, Month.soi))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily_nfish[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a, by.x=c("YEAR.c","MONTH"),
                                            by.y=c("Year.soi","Month.soi"),all.x=T)
  
  a=Freo
  a$dummy=with(a,paste(Year, Month))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily_nfish[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a, by.x=c("YEAR.c","MONTH"),
                                            by.y=c("Year","Month"),all.x=T)
  rm(a)
  
  #Add moon phase to daily records
  DATA.list.LIVEWT.c.daily_nfish[[i]]$Moon=lunar.phase(as.Date(DATA.list.LIVEWT.c.daily_nfish[[i]]$date),name=4)
}

#Explore cpue by block
HnDl="C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/Boxplots_year_blocks/"
for(i in 1:N.species)
{
  pdf(paste(HnDl,names(DATA.list.LIVEWT.c)[i],"_daily",".pdf",sep=""))
  fn.box.plt.year(d=subset(DATA.list.LIVEWT.c.daily[[i]],Catch.Target>0))
  dev.off()
  
  pdf(paste(HnDl,names(DATA.list.LIVEWT.c)[i],"_monthly",".pdf",sep=""))
  fn.box.plt.year(d=subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0))
  dev.off()
}


#4.12 Drop first years of sandbar data because vessels don't meet selection 
#   criteria and no positive catch
DD=DATA.list.LIVEWT.c$sandbar
DD=subset(DD,BLOCKX%in%as.numeric(BLKS.used$`Sandbar shark`))      
DD=subset(DD,VESSEL%in%VES.used$`Sandbar shark`)
a=with(DD,table(FINYEAR,VESSEL))
a[a>0]=1
a=rowSums(a)
a[a<Threshold.n.vessls.per.yr]=NA
San.Yrs=names(a[which(!is.na(a))])
rm(DD)
DATA.list.LIVEWT.c$sandbar=subset(DATA.list.LIVEWT.c$sandbar,FINYEAR%in%San.Yrs)


#4.13 Corroborate effective area  
if(Model.run=="First")for(s in 1:N.species)
{
  #Daily
  fn.check.eff.area(ALL=DATA.list.LIVEWT.c.daily,WHAT="TSNo",
                    Effrt='km.gillnet.days.c',
                    QL_prop_ktch=.9,
                    spname=casefold(unlist(strsplit(SPECIES.vec[s], " "))[1]),
                    TYPE='Daily')
}


#4.14  Identify targeting behaviour   
props="YES"  #response variable as proportion
SQRT="NO"   #transform
CENTRE="YES"
SCALE="YES"

use.what="catch"  #use proportional catch to be independent of abundance
#use.what="cpue"

Ktch.var=c("Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar")

agg_PCA="by_Same.return.SNo"
HndL.Species_targeting="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Species_targeting/"


  # #Hoyle et al 2015
  # Targt=c(17003,17001,18003,18007)
  # for(i in 1:N.species)  fn.targeting(Species.list.daily[[i]],SP=Targt[i])  
  

  
  #CLARA analysis as per Campbell et al 2017 on nfish as this has data at Sesssion level
  #CLARA clustering
  #The CLARA (Clustering Large Applications) algorithm is an extension to the 
  # PAM (Partitioning Around Medoids) clustering method for large data sets. It intended to 
  # reduce the computation time in the case of large data set.
  
  #As almost all partitioning algorithm, it requires the user to specify the appropriate
  # number of clusters to be produced. This can be estimated using the function fviz_nbclust.
  
  #The R function clara() [cluster package] can be used to compute CLARA algorithm. 
  # The simplified format is clara(x, k, pamLike = TRUE), 
  #  where "x" is the data and k is the number of clusters to be generated.
  
  #After, computing CLARA, the R function fviz_cluster() can be used to visualize the results. 
  # The format is fviz_cluster(clara.res), where clara.res is the CLARA results.
  
  #References
  #Kaufman, Leonard, and Peter Rousseeuw. 1990. Finding Groups in Data: An Introduction to Cluster Analysis.
  
  
  #note:Cluster is done on combined data set. Then select Same.return.SNo accordingly
  if(agg_PCA=="by_TSNo")
  {
    ALL.shots=unique(DATA.list.LIVEWT.c.daily$whiskery$Same.return.SNo)
    gum.shots=unique(DATA.list.LIVEWT.c.daily$gummy$Same.return.SNo)
    dus.shots=unique(DATA.list.LIVEWT.c.daily$dusky$Same.return.SNo)
    san.shots=unique(DATA.list.LIVEWT.c.daily$sandbar$Same.return.SNo)
    gum.shots=gum.shots[which(!gum.shots%in%ALL.shots)]
    dus.shots=dus.shots[which(!dus.shots%in%ALL.shots)]
    san.shots=san.shots[which(!san.shots%in%ALL.shots)]
    ALL=DATA.list.LIVEWT.c.daily$whiskery
    if(length(gum.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily$gummy,Same.return.SNo%in%gum.shots))
    if(length(dus.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily$dusky,Same.return.SNo%in%dus.shots))
    if(length(san.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily$sandbar,Same.return.SNo%in%san.shots))
    
    d=ALL[,match(c(Ktch.var),names(ALL))]
    d.cpue=as.data.frame(d)
    d.cpue=d.cpue[,match(Ktch.var,names(d.cpue))]
  }else
  {
    ALL.shots=unique(DATA.list.LIVEWT.c.daily_nfish$whiskery$Same.return.SNo)
    gum.shots=unique(DATA.list.LIVEWT.c.daily_nfish$gummy$Same.return.SNo)
    dus.shots=unique(DATA.list.LIVEWT.c.daily_nfish$dusky$Same.return.SNo)
    san.shots=unique(DATA.list.LIVEWT.c.daily_nfish$sandbar$Same.return.SNo)
    gum.shots=gum.shots[which(!gum.shots%in%ALL.shots)]
    dus.shots=dus.shots[which(!dus.shots%in%ALL.shots)]
    san.shots=san.shots[which(!san.shots%in%ALL.shots)]
    ALL=DATA.list.LIVEWT.c.daily_nfish$whiskery
    if(length(gum.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily_nfish$gummy,Same.return.SNo%in%gum.shots))
    if(length(dus.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily_nfish$dusky,Same.return.SNo%in%dus.shots))
    if(length(san.shots)>0) ALL=rbind(ALL,subset(DATA.list.LIVEWT.c.daily_nfish$sandbar,Same.return.SNo%in%san.shots))
    
    if(use.what=="catch")d=ALL[,match(c(Ktch.var,"Same.return.SNo"),names(ALL))]
    if(use.what=="cpue")
    {
      d=ALL[,match(c(Ktch.var,"Km.Gillnet.Days.c","Same.return.SNo"),names(ALL))]
      d[,match(Ktch.var,names(d))]=d[,match(Ktch.var,names(d))]/d$Km.Gillnet.Days.c
      d=d[,-match("Km.Gillnet.Days.c",names(d))]
    }
    
    d.cpue=as.data.frame(d[,-match("Same.return.SNo",names(d))])
    row.names(d.cpue)=d$Same.return.SNo
  }
  
  id=which(is.na(d.cpue))
  if(length(id)>0)d.cpue[id]=0
  id=which(rowSums(d.cpue)==0)
  if(length(id)>0)d.cpue=d.cpue[-id,]
  
  if(props=="YES") d.cpue=d.cpue/rowSums(d.cpue)  #convert to proportion  
  if(SQRT=="YES") d.cpue=d.cpue^0.5              #square root
  
  df=d.cpue
  if(CENTRE=="YES" & SCALE=="NO") df <- scale(df,center = TRUE, scale = FALSE)
  if(CENTRE=="YES" & SCALE=="YES") df <- scale(df,center = TRUE, scale = TRUE)
  rm(ALL.shots)
  
  #First. Determine if data is clusterable or not
  # Compute Hopkins statistic 
  #source: http://www.sthda.com/english/articles/29-cluster-validation-essentials/95-assessing-clustering-tendency-essentials/
  #fviz_pca_ind(prcomp(df), title = "PCA",geom = "point", ggtheme = theme_classic())   #any obvious clustering??
  
  test.clusterable="NO"
  if(test.clusterable=="YES")
  {
    fn.Clusterable=function(Hopkins)
    {
      if(1-Hopkins>0.75) a="clusterable"else
        a="Not clusterable"
      return(a)
    }
    ii=sample(1:nrow(df),2000,replace=F)  #random subsample to reduce computation time
    system.time({res <- get_clust_tendency(df[ii,], n = nrow(df[ii,])-1, graph = FALSE)})
    Clusterable=fn.Clusterable(res$hopkins_stat)
    write.csv(Clusterable,paste(HndL.Species_targeting,"Cluster/Clusterable.csv",sep=""))
  }
  
  #Second. Determine optimal number of clusters
  memory.limit(50000)
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_optimal_numbers",sep=""),2400,2400)
  a=fviz_nbclust(df, clara, method = "silhouette",print.summary=T)
  a+theme_classic()
  dev.off()
  
  #compute CLARA with PAM algorithm 
  num.clus=as.numeric(as.character(a$data$clusters[match(max(a$data$y),a$data$y)]))
  clara.res <- clara(df, num.clus, samples = 50, pamLike = TRUE)
  rm(a)
  # Print components of clara.res
  #print(clara.res)
  #  medoids: Objects that represent clusters
  #  clustering: a vector containing the cluster number of each object
  #  sample: labels or case numbers of the observations in the best sample, that is, 
  #         the sample used by the clara algorithm for the final partition
  
  
  
  #visualize CLARA clusters in data scattergram
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster",sep=""),2400,2400)
  fviz_cluster(clara.res, 
               palette = rainbow(num.clus), # color palette
               ellipse.type = "t", # Concentration ellipse
               geom = "point", pointsize = 1,
               ggtheme = theme_classic())
  dev.off()
  
  #add cluster to input data
  dd.clara <- cbind(as.data.frame(d.cpue), cluster_clara = clara.res$cluster)
  dd.clara$Same.return.SNo=row.names(dd.clara)
  row.names(dd.clara)=NULL
  
  #check cpue by group
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster.boxplot",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,3,3,1),oma=c(2,1,.1,.5),las=1,mgp=c(2,.6,0))
  boxplot(Catch.Whiskery~cluster_clara,dd.clara,ylab="Whiskery prop. catch")
  boxplot(Catch.Gummy~cluster_clara,dd.clara,ylab="Gummy prop. catch")
  boxplot(Catch.Dusky~cluster_clara,dd.clara,ylab="Dusky prop. catch")
  boxplot(Catch.Sandbar~cluster_clara,dd.clara,ylab="Sandbar prop. catch")
  dev.off()
  
  
  #add cluster to original data for use in standardisations
  for(i in 1:N.species)
  {
    a=subset(dd.clara,Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily_nfish[[i]]$Same.return.SNo),
             select=c(Same.return.SNo,cluster_clara))
    DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a,by="Same.return.SNo",all.x=T)
    DATA.list.LIVEWT.c.daily_nfish[[i]]$cluster_clara=as.character(DATA.list.LIVEWT.c.daily_nfish[[i]]$cluster_clara)
    DATA.list.LIVEWT.c.daily_nfish[[i]]$cluster_clara=with(DATA.list.LIVEWT.c.daily_nfish[[i]],
                                            ifelse(is.na(cluster_clara),"Other",cluster_clara))
  }
  rm(clara.res,dd.clara)
  
  do.normal.cluster="NO"
  if(do.normal.cluster=="YES")
  {
    
    # Dissimilarity matrix
    d <- dist(df, method = "euclidean")
    
    # Hierarchical clustering using Ward's method
    res.hc <- hclust(d, method = "ward.D2" )
    rm(d)  #huge object
    
    #plot groups on dendogram
    fn.fig(paste(HndL.Species_targeting,"Cluster/Standard_CLUSTER_dendogram",sep=""),2400,2400)
    #fviz_dend(res.hc, k =num.clus, k_colors = "jco",
    #          as.ggplot = TRUE, show_labels = FALSE)
    plot(res.hc, cex = 0.1, hang = -1,xlab="")
    rect.hclust(res.hc, k = num.clus, border = rainbow(num.clus))
    dev.off()
    
    #add cluster group to data
    grp <- cutree(res.hc, k = num.clus)  # Cut tree into k groups
    dd=as.data.frame(d.cpue)
    dd$Group=grp
    
    #visualize results in scatterplot
    fn.fig(paste(HndL.Species_targeting,"Cluster/Standard_CLUSTER_scatterplot",sep=""),2400,2400)
    fviz_cluster(list(data = df, cluster = grp))
    dev.off()
    
    #See weight of indicator species by group
    fn.fig(paste(HndL.Species_targeting,"Cluster/Standard_CLUSTER_boxplot",sep=""),2400,2400)
    par(mfcol=c(2,2),mar=c(1,3,3,1),oma=c(2,1,.1,.5),las=1,mgp=c(2,.6,0))
    boxplot(Catch.Whiskery~Group,dd,ylab="Whiskery")
    boxplot(Catch.Gummy~Group,dd,ylab="Gummy")
    boxplot(Catch.Dusky~Group,dd,ylab="Dusky")
    boxplot(Catch.Sandbar~Group,dd,ylab="Sandbar")
    dev.off()
    
    #add cluster to original data
    dd$Same.return.SNo=row.names(dd)
    row.names(dd)=NULL
    for(i in 1:N.species)
    {
      a=subset(dd,Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily_nfish[[i]]$Same.return.SNo),
               select=c(Same.return.SNo,Group))
      DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a,by="Same.return.SNo",all.x=T)
    }
    rm(res.hc,dd,grp)
  }
  
  rm(ALL,df,d.cpue,d)

lsos()

#4.15 Table of sensitivity scenarios       
Tab.Sensi=data.frame(Scenario=c("Base case","3 years","7 years","Shots","No efficiency"),
                     Vessels_used=c("5 years","3 years","7 years","5 years","5 years"),
                     Blocks_used=c("5 years","3 years","7 years"," 5 years","5 years"),
                     Effort=c(rep("days X net X hours",3),"days X net X hours X shots","days X net X hours"),
                     Efficiency_increase=c(rep("Yes",4),"No"))

setwd(paste(getwd(),"/Outputs/Paper",sep=""))
fn.word.table(WD=getwd(),TBL=Tab.Sensi,Doc.nm="Sensitivity tests",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")



#Add mesh, depth and nlines categories to daily_nfish  
for(s in 1:N.species)
{
  DATA.list.LIVEWT.c.daily_nfish[[s]]$DepCat=with(DATA.list.LIVEWT.c.daily_nfish[[s]],trunc(Mean.depth/10) * 10)
  DATA.list.LIVEWT.c.daily_nfish[[s]]$LineCat=with(DATA.list.LIVEWT.c.daily_nfish[[s]],
                         ifelse(nlines.c>2,"3_or_more",as.character(nlines.c)))
  DATA.list.LIVEWT.c.daily_nfish[[s]]$MeshCat=with(DATA.list.LIVEWT.c.daily_nfish[[s]],
                         ifelse(!mesh%in%c(165,178),"other",as.character(mesh)))
}


#4.16 Compute foly and nominal index for exporting  

  #--Foly
#note: As done by Rory, Foly, ie Effective, has all records (good  & bad reporters) from all blocks 
#       and vessels within effective area without 0 catches

#monthly  
DATA.list.LIVEWT.c_all_reporters=DATA.list.LIVEWT.c
for ( i in 1:N.species)
{
  #create data sets 
  dummy=Effort.data.fun(Species.list[[i]],TARGETS[[i]],"LIVEWT.c")  
  DATA.list.LIVEWT.c_all_reporters[[i]]=dummy$dat  
}

#daily 
DATA.list.LIVEWT.c.daily_all_reporters=DATA.list.LIVEWT.c.daily_all_reporters_nfish=DATA.list.LIVEWT.c.daily
for ( i in 1:N.species)
{
  #create data sets 
  dummy=Effort.data.fun.daily(Species.list.daily[[i]],TARGETS[[i]],"LIVEWT.c",Aggregtn="TSNo")
  DATA.list.LIVEWT.c.daily_all_reporters[[i]]=dummy$dat
  
  dummy=Effort.data.fun.daily(Species.list.daily[[i]],TARGETS[[i]],ktch="nfish",Aggregtn="SNo")
  DATA.list.LIVEWT.c.daily_all_reporters_nfish[[i]]=dummy$dat
  rm(dummy)
}

Foly.whis=export.foly(DATA.list.LIVEWT.c_all_reporters$whiskery,DATA.list.LIVEWT.c.daily_all_reporters$whiskery)
Foly.gum=export.foly(DATA.list.LIVEWT.c_all_reporters$gummy,DATA.list.LIVEWT.c.daily_all_reporters$gummy)
Foly.dus=export.foly(DATA.list.LIVEWT.c_all_reporters$dusky,DATA.list.LIVEWT.c.daily_all_reporters$dusky)
Foly.san=export.foly(DATA.list.LIVEWT.c_all_reporters$sandbar,DATA.list.LIVEWT.c.daily_all_reporters$sandbar)
List.foly.nom=list(whis=list(Foly=Foly.whis),gum=list(Foly=Foly.gum),dus=list(Foly=Foly.dus),san=list(Foly=Foly.san))


#-- Nominal 
#note: Ratio = mean(catch)/mean(effort)
#      Mean = mean(cpue)
#     LnMean= exp(mean(log(cpue))+bias corr)
#     DLnMean = exp(log(prop pos)+exp(mean(log(cpue))+bias corr)

QL_expl_ktch_prop=.9   #proportion of explained annual catch for selected record
PRP.MLCLM=0.1          #proportion of catch of target

Store_nom_cpues_monthly=Store_nom_cpues_daily=Store_nom_cpues_daily_nfish=DATA.list.LIVEWT.c.daily
Hnd.ains="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Ainsline_different_cpues/"

for(s in 1:N.species)
{
  #Monthly
  Store_nom_cpues_monthly[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c[[s]],Ktch.targt='catch.target',
                   Effrt=c('km.gillnet.days.c','km.gillnet.hours.c'),
                  explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
                 cpue.units = c("kg/km gillnet day","kg/km gillnet hour"),spname=SPECIES.vec[s],
                 BLks=as.numeric(BLKS.used[[s]]),VesL=VES.used[[s]],Type="_monthly_")
  
  #Daily weight
  Store_nom_cpues_daily[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily[[s]],Ktch.targt='catch.target',
                Effrt=c('km.gillnet.days.c','km.gillnet.hours.c'),
               explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
               cpue.units = c("kg/km gillnet day","kg/km gillnet hour"),spname=SPECIES.vec[s],
              BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Type="_daily_")
  
  #Daily numbers
  Store_nom_cpues_daily_nfish[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily_nfish[[s]],Ktch.targt='catch.target',
             Effrt=c('km.gillnet.days.c','km.gillnet.hours.c'),
            explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
             cpue.units = c("Individuals/km gillnet day","Individuals/km gillnet hour"),spname=SPECIES.vec[s],
          BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Type="_daily_n_")
}


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
for(s in 1:N.species)
{
  #monthly
  a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=SPECIES.vec[s],
                      what="monthly",MN.YR=Min.Vess.yr,pLot=T)
  BLKS.used[[s]]=a$this.blks
  VES.used[[s]]=a$this.ves
  write.csv(BLKS.used[[s]],paste(hndl.kept,"blocks_used_",SPECIES.vec[s],"_monthly.csv",sep=""))
  
  #daily
  a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=SPECIES.vec[s],
                      what="daily",MN.YR=Min.Vess.yr.d,pLot=T)
  BLKS.used.daily[[s]]=a$this.blks
  VES.used.daily[[s]]=a$this.ves
  write.csv(BLKS.used.daily[[s]],paste(hndl.kept,"blocks_used_",SPECIES.vec[s],"_daily.csv",sep=""))
}

#show blocks kept
CEX=.85
SRt=45
  
fn.fig(paste(hndl.kept,"block_used_map",sep=""),2400, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in 1:N.species)
{
  #Monthly
  fn.show.blk(dat=BLKS.used[[s]],CEX=CEX,SRt=SRt)
  if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  fn.show.blk(dat=BLKS.used.daily[[s]],CEX=CEX,SRt=SRt)
  if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",SPECIES.vec[s],bty='n',cex=1.5)
}
mtext("Longitude",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Latitude",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")

#Show proportion of selected records by year
fn.fig("Figure QL_proportion_records_selected",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in 1:N.species)
{
  #Monthly
  with(subset(Store_nom_cpues_monthly[[s]]$Prop.selected,target==1),plot(yr,freq,ylim=c(0,1),
                          las=1,pch=19,cex=1.25,ylab="",xlab=""))
  if(s==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)

  #Daily
  with(subset(Store_nom_cpues_daily[[s]]$Prop.selected,target==1),plot(yr,freq,ylim=c(0,1),
              las=1,pch=19,cex=1.25,ylab="",xlab=""))
  if(s==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",SPECIES.vec[s],bty='n',cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Proportion of records selected",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


#4.18 Show gummy monthly cpue effect of using km gn d or km g h
if(Model.run=="First")
{
  fn.fig("Gummy_km.gn.d_VS_km.gn.h",2400,2400)
  par(mfrow=c(4,2),mar=c(1,3,2,1),oma=c(2,.5,.5,1.75),mgp=c(1.5,.5,0),cex.lab=1.5)
  Get.Mns(d=DATA.list.LIVEWT.c$gummy,grp="FINYEAR",
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
  for(s in 1:N.species)
  {
    write.csv(DATA.list.LIVEWT.c[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c[s]),"_monthly.csv",sep=""),row.names=F)
    write.csv(DATA.list.LIVEWT.c.daily[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c.daily[s]),"_daily.csv",sep=""),row.names=F)
    write.csv(DATA.list.LIVEWT.c.daily_nfish[[s]],paste("For_Ainslie/",names(DATA.list.LIVEWT.c.daily_nfish[s]),"_daily_nfish.csv",sep=""),row.names=F)
  }
}



#4.20 Output data tables  
#Output table with number of records available in effective area and numbers used in standardisation
TABle=vector('list',N.species)
names(TABle)=SPECIES.vec
for(s in 1:N.species)
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

Table.nsamp <- Reduce(function(x, y) merge(x, y, all=T,by=c("Year", "Record")), TABle, accumulate=F)
Table.nsamp=Table.nsamp[,c("Year", "Record",sort(names(Table.nsamp[3:ncol(Table.nsamp)])))]
Export.tbl(WD=getwd(),Tbl=Table.nsamp,Doc.nm="Sample_sizes",caption=NA,paragph=NA,
           HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
           Zebra='NO',Zebra.col='grey60',Grid.col='black',
           Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
           HDR.names=c('Year', 'Record','Total number of records within eff. area',
                       'Number of records for good reporters','Number of records used in stand.'),
           HDR.span=c(1,1,N.species,N.species,N.species),
           HDR.2nd=c("","",rep(c("Dusky","Gummy","Sandbar","Whiskery"),3)))

rm(DATA.list.LIVEWT.c_all_reporters,DATA.list.LIVEWT.c.daily_all_reporters)


#4.21 Check outliers in catch and effort for removing nonsense values   
#note: max monthly ktch, effort (~ 40 tonnes, ~ 5800 km gn h (@ 30 days X 24 h X 8000 m), respectively) 
#      max trip (daily kg) ktch, effort (~ 15 tonnes, ~ 1900 km gn h (@ 10 days X 24 h X 8000 m), respectively) 
pdf("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/Check.outliers.pdf") 
for(i in 1:N.species)
{
  check.cpue(Store_nom_cpues_monthly[[i]]$QL_dat,paste("monthly",names(DATA.list.LIVEWT.c)[i]),4)
  check.cpue(Store_nom_cpues_daily[[i]]$QL_dat,paste("daily",names(DATA.list.LIVEWT.c)[i]),4)  
}
dev.off()


#4.22 Construct index of abundance     
ZONES=c("West","Zone1","Zone2")
Eff.vars=c("km.gillnet.hours.c","km.gillnet.hours_shot.c")
Covariates=c("freo")
Predictors_monthly=c("finyear","vessel","month","blockx",Covariates) 
Predictors_daily=Predictors_monthly
Predictors_daily_nfish=c(Predictors_daily,"cluster_clara","meshcat","linecat","depcat","moon")
Response="catch.target"    #note that cpue is calculated inside stand function

Categorical=c("finyear","vessel","month","blockx",
              "cluster_clara","moon","depcat","linecat","meshcat")

  

#   4.22.1 Explore data used for standardisation

  #check data properties and degrees of freedom   
if(Model.run=="First")
{
  Prop.deg.free.m=Prop.deg.free.d=Prop.deg.free.d_n=Store_nom_cpues_monthly
  for(s in 1:N.species)
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
    
    dummy=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    Prop.d_n=properties(dummy)
    d.f=Prop.d_n[match(Predictors_daily_nfish,rownames(Prop.d_n)),]
    d.f$dummy=as.numeric(with(d.f,ifelse(Class%in%c('numeric'),1,Unique)))
    d.f.d.n=data.frame(Deg.F=sum(d.f$dummy),Obser=nrow(dummy))
    Prop.deg.free.d_n[[s]]=list(dat=d.f,deg.f=d.f.d.n)
  }
}

  #apply cede() functions    
if(Model.run=="First")
{
  hndl.cede="C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/cede/"
  for(s in 1:N.species)
  {
    #monthly
    pdf(paste(hndl.cede,SPECIES.vec[s],"_monthly.pdf",sep="")) 
    dummy=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used[[s]])
    fn.expl.cede(d=dummy,PREDS=Predictors_monthly,kg=TRUE,Do.ggplts=FALSE)
    dev.off()
    
    #daily
    pdf(paste(hndl.cede,SPECIES.vec[s],"_daily.pdf",sep=""))
    dummy=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    fn.expl.cede(d=dummy,PREDS=Predictors_daily,kg=TRUE,Do.ggplts=FALSE)
    dev.off()
    
    #daily nfish
    pdf(paste(hndl.cede,SPECIES.vec[s],"_daily_nfish.pdf",sep=""))
    dummy=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    fn.expl.cede(d=dummy,PREDS=Predictors_daily_nfish,kg=FALSE,Do.ggplts=FALSE)
    dev.off()
    rm(dummy)
  }
}

  #show records dropped by data selection process (starting from good records)
if(Model.run=="First")
{
  fn.fig("show_how_records_drop",2400, 2400) 
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,1),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
    first=table(DATA.list.LIVEWT.c[[s]]$FINYEAR)
    dummy=Store_nom_cpues_monthly[[s]]$QL_dat
    second=table(dummy$finyear)
    dummy=subset(dummy,vessel%in%VES.used[[s]])
    third=table(dummy$finyear)
    dummy=subset(dummy,blockx%in%BLKS.used[[s]])
    fourth=table(dummy$finyear)
    LGTXT=NULL
    if(s==4)LGTXT=c("whole data set","QL","QL_sel. vess.","QL_sel. vess. & block")
    barplot(rbind(first,second,third,fourth),legend.text=LGTXT,
            args.legend=list(x="topright",cex=.9,bty='n',xjust=0))
    box()
    
    first=table(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR)
    dummy=Store_nom_cpues_daily[[s]]$QL_dat
    second=table(dummy$finyear)
    dummy=subset(dummy,vessel%in%VES.used.daily[[s]])
    third=table(dummy$finyear)
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    fourth=table(dummy$finyear)
    barplot(rbind(first,second,third,fourth))
    box()
    mtext(SPECIES.vec[s],4,0.5,cex=1.5,las=3)
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

      #4.22.3.1. select dist for count
Disc.dist=SPECIES.vec
names(Disc.dist)=SPECIES.vec
for(s in 1:N.species)
{
  dummy=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
  Disc.dist[s]=fn.sel.discrete.dist(dummy$catch.target)
}


      #4.22.3.2. extract best model
Best.Model=vector('list',N.species)
names(Best.Model)=SPECIES.vec
Best.Model.daily=Best.Model.daily_nfish=
Store.Best.Model=Store.Best.Model.daily=Store.Best.Model.daily_nfish=Best.Model

if(Def.mod.Str=="YES")     #takes 40 minutes
{
  efrt="km.gillnet.hours.c"
  system.time({
  for(s in 1:N.species)
  {
    fitFun=function(formula, data,...) glm(formula,data=data,...) 
    
    #monthly
    d=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    d=subset(d,blockx%in%BLKS.used[[s]])
    PREDS=Predictors_monthly
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    Store.Best.Model[[s]]=fn.modl.sel(Inter="MainTerm",RESPNS="LNcpue")
    Best.Model[[s]]=Store.Best.Model[[s]]$BEST
    rm(d)
    
    #daily
    d=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    d=subset(d,blockx%in%BLKS.used.daily[[s]])
    PREDS=Predictors_daily
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    Store.Best.Model.daily[[s]]=fn.modl.sel(Inter="MainTerm",RESPNS="LNcpue")
    Best.Model.daily[[s]]=Store.Best.Model.daily[[s]]$BEST 
    rm(d,fitFun)
    
    
    #daily nfish
    if(Disc.dist[s]=="neg.bin") fitFun=function(formula, data) glm.nb(formula,data=data)
    if(Disc.dist[s]=="Poisson") fitFun=function(formula, data,family) glm(formula,data=data,family="poisson")
    d=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    d=subset(d,blockx%in%BLKS.used.daily[[s]])
    PREDS=Predictors_daily_nfish
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    Store.Best.Model.daily_nfish[[s]]=fn.modl.sel(Inter="MainTerm",RESPNS="catch")
    Best.Model.daily_nfish[[s]]=Store.Best.Model.daily_nfish[[s]]$BEST  
    
    rm(d,fitFun)
    
  }
  
  #4.22.3.3. show selection outcomes
  hndl.modl.sel="C:/Matias/Analyses/Catch and effort/Outputs/Model Selection/"
  for(s in 1:N.species)
  {
    pdf(paste(hndl.modl.sel,SPECIES.vec[s],"_monthly.pdf",sep=""))
    fn.show.mod.sel(MODS=Store.Best.Model[[s]]$res,outs=20)
    dev.off()
    
    pdf(paste(hndl.modl.sel,SPECIES.vec[s],"_daily.pdf",sep=""))
    fn.show.mod.sel(MODS=Store.Best.Model.daily[[s]]$res,outs=20)
    dev.off()
    
    pdf(paste(hndl.modl.sel,SPECIES.vec[s],"_daily_nfish.pdf",sep=""))
    fn.show.mod.sel(MODS=Store.Best.Model.daily_nfish[[s]]$res,outs=20)
    dev.off()
  }
  
  })
}   

if(Def.mod.Str=="NO")
{
  for(s in 1:N.species)
  {
    Best.Model[[s]]=formula("LNcpue ~ finyear + vessel + month + blockx")
    Best.Model.daily[[s]]=formula("LNcpue ~ finyear + vessel + month + blockx")
    if(names(Best.Model.daily_nfish)[s]=="Whiskery shark")
    {
      Best.Model.daily_nfish[[s]]=formula("catch.target~finyear+vessel+month+blockx+
                                    cluster_clara+offset(LNkm.gillnet.hours.c)")
    }
    if(names(Best.Model.daily_nfish)[s]=="Gummy shark")
    {
      Best.Model.daily_nfish[[s]]=formula("catch.target~finyear+vessel+month+blockx+
                                    cluster_clara+offset(LNkm.gillnet.hours.c)")
    }
    if(names(Best.Model.daily_nfish)[s]=="Dusky shark")
    {
      Best.Model.daily_nfish[[s]]=formula("catch.target~finyear+vessel+month+blockx+
                                    cluster_clara+depcat+offset(LNkm.gillnet.hours.c)")
    }
    if(names(Best.Model.daily_nfish)[s]=="Sandbar shark")
    {
      Best.Model.daily_nfish[[s]]=formula("catch.target~finyear+vessel+month+blockx+
                                    cluster_clara+depcat+offset(LNkm.gillnet.hours.c)")
    }
  }
}

#visualise coefficients daily_nfish
if(Def.mod.Str=="YES")
{
  for(s in 1:N.species)
  {
    d=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    d=subset(d,blockx%in%BLKS.used.daily[[s]])
    PREDS=Predictors_daily_nfish
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    PREDS[id.cov]=paste("LN",PREDS[id.cov],sep="")
    Formula=formula(paste(Response,paste(paste(PREDS,collapse="+"),
                   paste("offset(","LN",efrt,")",sep=""),sep="+"),sep="~"))
    
    MOD=glm.nb(Formula,data=d)
    AOV=anova(MOD,test='Chisq')
    
    pdf(paste(hndl.modl.sel,"viz.coeffs_",SPECIES.vec[s],".pdf",sep=""))
    viz.coef(MOD,
             WHAT=c(paste("depcat",levels(d$depcat),sep=""),
                    paste("cluster_clara",levels(d$cluster_clara),sep=""),
                    paste("linecat",levels(d$linecat),sep=""),
                    paste("moon",levels(d$moon),sep="")))
    
    par(mfcol=c(2,2),mar=c(2,2,.25,.15))
    boxplot((d$catch.target/d$km.gillnet.hours.c)~d$depcat,cex.axis=.8)
    boxplot((d$catch.target/d$km.gillnet.hours.c)~d$cluster_clara,cex.axis=.8)
    boxplot((d$catch.target/d$km.gillnet.hours.c)~d$linecat,cex.axis=.8)
    boxplot((d$catch.target/d$km.gillnet.hours.c)~d$moon,cex.axis=.8)
    
    
    Expl.Dev=data.frame(Term=rownames(AOV),Deviance=AOV$Deviance,
                        Res.Dev=AOV$'Resid. Dev')
    Expl.Dev$Expl.Dev=paste(round((Expl.Dev$Deviance*100/Expl.Dev$Res.Dev[1]),2),"%")
    plot.new()
    grid.table(Expl.Dev)
    
    dev.off()
  }
}


#Export table of term levels
terms.table=vector('list',length(N.species))
for(s in 1:N.species)   
{
  #monthly
  DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used[[s]])
  Mon=fn.table.terms(d=DAT,PREDS=Predictors_monthly)
  Mon=cbind(Data="Monthly",Mon)
  
  #daily
  DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
  Day=fn.table.terms(d=DAT,PREDS=Predictors_daily)
  Day=cbind(Data="Daily",Day)
  
  #daily_nfish
  DAT=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
  Day_nfish=fn.table.terms(d=DAT,PREDS=Predictors_daily_nfish)
  Day_nfish=cbind(Data="Daily_nfish",Day_nfish)
  DAT=rbind(Mon,Day,Day_nfish)
 
   terms.table[[s]]=cbind(Species=SPECIES.vec[s],DAT)
}
Tab.Terms=do.call(rbind,terms.table)
row.names(Tab.Terms)=NULL
fn.word.table(WD=getwd(),TBL=Tab.Terms,Doc.nm="Terms_table",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")




#   4.22.4 Run standardisation
Stand.out=vector('list',length=N.species)
names(Stand.out)=SPECIES.vec 
Stand.out.daily=Stand.out.daily_nfish=Stand.out

system.time({for(s in 1:N.species)   #takes 20 secs
{
  #monthly
  DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used[[s]])
  Stand.out[[s]]=fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
           efrt="km.gillnet.hours.c",Formula=Best.Model[[s]])
  rm(DAT)
  
  #daily
  DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
  Stand.out.daily[[s]]=fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",
            PREDS=Predictors_daily,efrt="km.gillnet.hours.c",Formula=Best.Model.daily[[s]])
  rm(DAT)
  
  #daily_nfish
  DAT=subset(Store_nom_cpues_daily_nfish[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
  DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
  Stand.out.daily_nfish[[s]]=fn.stand(d=DAT,Response="catch.target",RESPNS="catch",
              PREDS=Predictors_daily_nfish,efrt="km.gillnet.hours.c",Formula=Best.Model.daily_nfish[[s]])
  rm(DAT)
}})


#   4.22.5 Export deviance explained
if(Model.run=="First")
{
  Dev.exp=vector('list',length=N.species)
  names(Dev.exp)=SPECIES.vec
  Dev.exp.daily=Dev.exp.daily_nfish=Dev.exp
  
  system.time({for(s in 1:N.species)
  {
    Dev.exp[[s]]=Anova.and.Dev.exp(GLM=Stand.out[[s]]$res,SP=SPECIES.vec[s],type="Monthly")
    Dev.exp.daily[[s]]=Anova.and.Dev.exp(GLM=Stand.out.daily[[s]]$res,SP=SPECIES.vec[s],type="Daily")
    Dev.exp.daily_nfish[[s]]=Anova.and.Dev.exp(GLM=Stand.out.daily_nfish[[s]]$res,SP=SPECIES.vec[s],type="Daily_nfish")
  }})   #takes 7 seconds
  
  Tab.Dev.Exp=rbind(do.call(rbind,Dev.exp),do.call(rbind,Dev.exp.daily),do.call(rbind,Dev.exp.daily_nfish))
  rownames(Tab.Dev.Exp)=NULL
  
  fn.word.table(WD=getwd(),TBL=Tab.Dev.Exp,Doc.nm="ANOVA_table",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
}



#   4.22.6 Run sensitivity tests     
if(Model.run=="First")      #takes 12 mins
{
  sens=Tab.Sensi
  sens$Efrt.used=with(sens,ifelse(Effort=="days X net X hours","km.gillnet.hours.c",
                      ifelse(Effort=="days X net X hours X shots","km.gillnet.hours_shot.c",NA)))
  Stand.out_sens=Stand.out.daily_sens=Stand.out
  
  system.time({for(s in 1:N.species)   #takes 25 secs
  {
    sens_monthly=vector('list',length=nrow(sens))
    names(sens_monthly)=sens$Scenario
    sens_daily=sens_monthly
    for(o in 1:nrow(sens))
    {
      MiN.YR=as.numeric(substr(sens$Vessels_used[o],1,1))
      EFrT=sens$Efrt.used[o]
        
      #monthly
      a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=SPECIES.vec[s],
                          what="monthly",MN.YR=MiN.YR,pLot=F)
      if(SPECIES.vec[s]=="Sandbar shark") a$this.blks=subset(a$this.blks,!a$this.blks=="3517") #cannot estimate this coef for 0==1
      DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%a$this.ves)
      DAT=subset(DAT,blockx%in%a$this.blks)
      sens_monthly[[o]]=fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",
                            PREDS=Predictors_monthly,efrt=EFrT,Formula=Best.Model[[s]])
      rm(DAT)
      
      #daily
      a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=SPECIES.vec[s],
                          what="daily",MN.YR=MiN.YR,pLot=F)

      DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%a$this.ves)
      DAT=subset(DAT,blockx%in%a$this.blks)
      sens_daily[[o]]=fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",
                       PREDS=Predictors_daily,efrt=EFrT,Formula=Best.Model.daily[[s]])
      rm(DAT)
    }
    Stand.out_sens[[s]]=sens_monthly
    Stand.out.daily_sens[[s]]=sens_daily
  }})
  
  
  #Predict years based on emmeans (formerly lsmeans) considering log bias corr if required
  Sens.pred=vector('list',length=N.species)
  names(Sens.pred)=SPECIES.vec 
  Sens.pred.daily=Sens.pred
  system.time({
    for(s in 1:N.species)   
    {
      dummy=vector('list',length=nrow(sens))
      names(dummy)=sens$Scenario
      dummy.daily=dummy
      for(o in 1:nrow(sens))
      {
        d=Stand.out_sens[[s]][[o]]$DATA   #note: need data as global for ref_grid
        dummy[[o]]=pred.fun(MOD=Stand.out_sens[[s]][[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
        
        d=Stand.out.daily_sens[[s]][[o]]$DATA
        dummy.daily[[o]]=pred.fun(MOD=Stand.out.daily_sens[[s]][[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
        rm(d)
      }
      Sens.pred[[s]]=dummy
      Sens.pred.daily[[s]]=dummy.daily
    }
  })   
  
  #Apply efficiency creep where required     
  Sens.pred.creep=Sens.pred
  Sens.pred.daily.creep=Sens.pred.daily
  for(s in 1:N.species)
  {
    for(o in 1:nrow(sens))
    {
      if(sens$Efficiency_increase[o]=="Yes")
      {
        #monthly
        add.crp=Eff.creep$effort.creep[match(Sens.pred.creep[[s]][[o]]$finyear,Eff.creep$finyear)]
        Sens.pred.creep[[s]][[o]]$response=Sens.pred.creep[[s]][[o]]$response*(1-add.crp)
        Sens.pred.creep[[s]][[o]]$lower.CL=Sens.pred.creep[[s]][[o]]$lower.CL*(1-add.crp)
        Sens.pred.creep[[s]][[o]]$upper.CL=Sens.pred.creep[[s]][[o]]$upper.CL*(1-add.crp)

        #daily
        add.crp=Eff.creep$effort.creep[match(Sens.pred.daily.creep[[s]][[o]]$finyear,Eff.creep$finyear)]
        Sens.pred.daily.creep[[s]][[o]]$response=Sens.pred.daily.creep[[s]][[o]]$response*(1-add.crp)
        Sens.pred.daily.creep[[s]][[o]]$lower.CL=Sens.pred.daily.creep[[s]][[o]]$lower.CL*(1-add.crp)
        Sens.pred.daily.creep[[s]][[o]]$upper.CL=Sens.pred.daily.creep[[s]][[o]]$upper.CL*(1-add.crp)
      }
    }
  }
  
  #Normalise    
  Sens.pred.normlzd=Sens.pred.creep
  Sens.pred.daily.normlzd=Sens.pred.daily.creep
  for(s in 1:N.species)
  {
    for(o in 1:nrow(sens))
    {
        #monthly
        Mn=mean(Sens.pred.normlzd[[s]][[o]]$response)
        Sens.pred.normlzd[[s]][[o]]$response=Sens.pred.normlzd[[s]][[o]]$response/Mn
        Sens.pred.normlzd[[s]][[o]]$lower.CL=Sens.pred.normlzd[[s]][[o]]$lower.CL/Mn
        Sens.pred.normlzd[[s]][[o]]$upper.CL=Sens.pred.normlzd[[s]][[o]]$upper.CL/Mn
        
        #daily
        Mn=mean(Sens.pred.daily.normlzd[[s]][[o]]$response)
        Sens.pred.daily.normlzd[[s]][[o]]$response=Sens.pred.daily.normlzd[[s]][[o]]$response/Mn
        Sens.pred.daily.normlzd[[s]][[o]]$lower.CL=Sens.pred.daily.normlzd[[s]][[o]]$lower.CL/Mn
        Sens.pred.daily.normlzd[[s]][[o]]$upper.CL=Sens.pred.daily.normlzd[[s]][[o]]$upper.CL/Mn
    }
  }
  
  #Plot   
  fn.fig("Appendix_Sensitivity",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
    LgND="NO"
    if(s==1)LgND="YES"
    Plot.cpue(cpuedata=Sens.pred.creep[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear")
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    LgND="NO"
    Plot.cpue(cpuedata=Sens.pred.daily.creep[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear")
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",SPECIES.vec[s],bty='n',cex=1.75)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  fn.fig("Appendix_Sensitivity_nomalised",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
    LgND="NO"
    if(s==1)LgND="YES"
    Plot.cpue(cpuedata=Sens.pred.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear")
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    LgND="NO"
    Plot.cpue(cpuedata=Sens.pred.daily.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
              COL="color",CxS=1.15,Yvar="finyear")
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",SPECIES.vec[s],bty='n',cex=1.75)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Relative CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
}


#   4.22.7 Fit diagnostics
if(Model.run=="First") 
{
  fn.fig("Appendix 6",2000, 2400)
  par(mfcol=c(2*3,4),las=1,mar=c(2,2,1.75,1),oma=c(1,2,.1,2),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
  for(s in 1:N.species) 
  {
    #Monthly
    Pos.Diag.fn(MODEL=Stand.out[[s]]$res,SPECIES=SPECIES.vec[s],M=.9)
    
    #Daily
    Pos.Diag.fn(MODEL=Stand.out.daily[[s]]$res,SPECIES="",M=.9)
  }
  mtext(c("Daily logbooks                                          Monthly returns     "),4,
        outer=T,las=3,line=0,cex=1.3)
  dev.off()
}



#   4.22.8 Plot base case and nominal 

#Extract comparable nominal cpue
Selected.data='CPUE.QL_target'
Selected.effort='km.gillnet.hours.c'

id.dat=match(Selected.data,names(Store_nom_cpues_monthly[[1]]))
id.eff=match(Selected.effort,names(Store_nom_cpues_monthly[[1]][[1]]))

Nominl=vector('list',length=N.species)
names(Nominl)=SPECIES.vec 
Nominl.daily=Nominl
for(s in 1:N.species) 
{
  Nominl[[s]]=subset(Store_nom_cpues_monthly[[s]][[id.dat]][[id.eff]],method=="LnMean")
  Nominl.daily[[s]]=subset(Store_nom_cpues_daily[[s]][[id.dat]][[id.eff]],method=="LnMean")
}

#Predict years based on emmeans (formerly lsmeans) considering log bias corr if required
Pred=vector('list',length=N.species)
names(Pred)=SPECIES.vec 
Pred.daily=Pred
system.time({
  for(s in 1:N.species)   
  {
    d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
    Pred[[s]]=pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="finyear",Pred.type="link")
      
    d=Stand.out.daily[[s]]$DATA
    Pred.daily[[s]]=pred.fun(MOD=Stand.out.daily[[s]]$res,biascor="YES",PRED="finyear",Pred.type="link")
    rm(d)
  }
})     #takes 3 minutes


#Apply efficiency creep      
Pred.creep=Pred
Pred.daily.creep=Pred.daily
Nominl.creep=Nominl
Nominl.daily.creep=Nominl.daily
for(s in 1:N.species)
{
  #monthly
  add.crp=Eff.creep$effort.creep[match(Pred.creep[[s]]$finyear,Eff.creep$finyear)]
  Pred.creep[[s]]$response=Pred.creep[[s]]$response*(1-add.crp)
  Pred.creep[[s]]$lower.CL=Pred.creep[[s]]$lower.CL*(1-add.crp)
  Pred.creep[[s]]$upper.CL=Pred.creep[[s]]$upper.CL*(1-add.crp)
  
  Nominl.creep[[s]]$response=Nominl.creep[[s]]$mean*(1-add.crp)
  Nominl.creep[[s]]$lower.CL=Nominl.creep[[s]]$lowCL*(1-add.crp)
  Nominl.creep[[s]]$upper.CL=Nominl.creep[[s]]$uppCL*(1-add.crp)
  
  
  #daily
  add.crp=Eff.creep$effort.creep[match(Pred.daily.creep[[s]]$finyear,Eff.creep$finyear)]
  Pred.daily.creep[[s]]$response=Pred.daily.creep[[s]]$response*(1-add.crp)
  Pred.daily.creep[[s]]$lower.CL=Pred.daily.creep[[s]]$lower.CL*(1-add.crp)
  Pred.daily.creep[[s]]$upper.CL=Pred.daily.creep[[s]]$upper.CL*(1-add.crp)
  
  Nominl.daily.creep[[s]]$response=Nominl.daily.creep[[s]]$mean*(1-add.crp)
  Nominl.daily.creep[[s]]$lower.CL=Nominl.daily.creep[[s]]$lowCL*(1-add.crp)
  Nominl.daily.creep[[s]]$upper.CL=Nominl.daily.creep[[s]]$uppCL*(1-add.crp)
}

fn.fig("Figure 4.Annual_Index",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
for(s in 1:N.species)
{
  #Monthly
  Mon.dat=list(Standardised=Pred.creep[[s]],Nominal=Nominl.creep[[s]])
  LgND="NO"
  if(s==1)LgND="YES"
  Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear")
  if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  Daily.dat=list(Standardised=Pred.daily.creep[[s]],Nominal=Nominl.daily.creep[[s]])
  LgND="NO"
  Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear")
  if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("CPUE (kg/ km gillnet hour)",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
dev.off()


#Normalise    
Pred.normlzd=Pred.creep
Pred.daily.normlzd=Pred.daily.creep
Nominl.normlzd=Nominl.creep
Nominl.daily.normlzd=Nominl.daily.creep
for(s in 1:N.species)
{
  #monthly
  Mn=mean(Pred.normlzd[[s]]$response)
  Pred.normlzd[[s]]$response=Pred.normlzd[[s]]$response/Mn
  Pred.normlzd[[s]]$lower.CL=Pred.normlzd[[s]]$lower.CL/Mn
  Pred.normlzd[[s]]$upper.CL=Pred.normlzd[[s]]$upper.CL/Mn
  
  Mn=mean(Nominl.normlzd[[s]]$response)
  Nominl.normlzd[[s]]$response=Nominl.normlzd[[s]]$response/Mn
  Nominl.normlzd[[s]]$lower.CL=Nominl.normlzd[[s]]$lower.CL/Mn
  Nominl.normlzd[[s]]$upper.CL=Nominl.normlzd[[s]]$upper.CL/Mn
  
  
  #daily
  Mn=mean(Pred.daily.normlzd[[s]]$response)
  Pred.daily.normlzd[[s]]$response=Pred.daily.normlzd[[s]]$response/Mn
  Pred.daily.normlzd[[s]]$lower.CL=Pred.daily.normlzd[[s]]$lower.CL/Mn
  Pred.daily.normlzd[[s]]$upper.CL=Pred.daily.normlzd[[s]]$upper.CL/Mn
  
  Mn=mean(Nominl.daily.normlzd[[s]]$response)
  Nominl.daily.normlzd[[s]]$response=Nominl.daily.normlzd[[s]]$response/Mn
  Nominl.daily.normlzd[[s]]$lower.CL=Nominl.daily.normlzd[[s]]$lower.CL/Mn
  Nominl.daily.normlzd[[s]]$upper.CL=Nominl.daily.normlzd[[s]]$upper.CL/Mn
  
  
}

fn.fig("Figure 4.Annual_Index_normalised",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
for(s in 1:N.species)
{
  #Monthly
  Mon.dat=list(Standardised=Pred.normlzd[[s]],Nominal=Nominl.normlzd[[s]])
  LgND="NO"
  if(s==1)LgND="YES"
  Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear")
  if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  Daily.dat=list(Standardised=Pred.daily.normlzd[[s]],Nominal=Nominl.daily.normlzd[[s]])
  LgND="NO"
  Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear")
  if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
dev.off()


#   4.22.9 Influence plots     
if(do.influence=="YES")
{
  HnDll="C:/Matias/Analyses/Catch and effort/Outputs/Influence.plot/"
  Store.Influence=Store.Influence.daily=Pred.creep
  for(s in 1:N.species)
  {
    #Monthly
    Terms=all.vars(Best.Model[[s]])[-1]
    Terms=Terms[-match("finyear",Terms)]
    Term.Type=ifelse(Terms%in%c("vessel" , "month" , "blockx"),"CAT",NA)
    
    pdf(paste(HnDll,SPECIES.vec[s],".monthly.CDI.pdf",sep=""))
    Store.Influence[[s]]=Influence.fn(MOD=Stand.out[[s]]$res,DAT=Stand.out[[s]]$DATA,
                  Term.type=Term.Type,termS=Terms,add.Influence="YES",SCALER=4)
    dev.off()
    
    
    #Daily
    Terms=all.vars(Best.Model.daily[[s]])[-1]
    Terms=Terms[-match("finyear",Terms)]
    Term.Type=ifelse(Terms%in%c("vessel" , "month" , "blockx"),"CAT",NA)
    
    pdf(paste(HnDll,SPECIES.vec[s],".daily.CDI.pdf",sep=""))
    Store.Influence.daily[[s]]=Influence.fn(MOD=Stand.out.daily[[s]]$res,
            DAT=Stand.out.daily[[s]]$DATA,Term.type=Term.Type,termS=Terms,
            add.Influence="YES",SCALER=4)
    dev.off()
  } 
  
  Over.inf.monthly=do.call(rbind,lapply(Store.Influence, '[[', match('Over.all.influence',names(Store.Influence[[1]]))))
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
    YLABs=termS
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
    mtext(c("Financial year       Coefficient"),side=2,line=.25,cex=1,las=3,outer=T)
    mtext(c("Dusky (monthly)          Sandbar (monthly)                                          Sandbar (daily)                                                                    "),
          3,line=-.25,cex=.75,outer=T)
  }
  dev.off()
  
}


#   4.22.10 Show month and block effects 
if(Model.run=="First")
{
  #Predict month based on emmeans (formerly lsmeans) considering log bias corr if required
  Pred.month=vector('list',length=N.species)
  names(Pred.month)=SPECIES.vec 
  Pred.month.daily=Pred.month
  system.time({
    for(s in 1:N.species)   
    {
      d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
      Pred.month[[s]]=pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="month",Pred.type="link")
      
      d=Stand.out.daily[[s]]$DATA
      Pred.month.daily[[s]]=pred.fun(MOD=Stand.out.daily[[s]]$res,biascor="YES",PRED="month",Pred.type="link")
      rm(d)
    }
})     #takes 1.5 minutes
  
  
  fn.fig("Figure 2.Monthly effect",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,2.5,.1,.2),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
    #Monthly
    Mon.dat=Pred.month[[s]]
    LgND="NO"
    if(s==1)LgND="YES"
    Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month")
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    
    #Daily
    Daily.dat=Pred.month.daily[[s]]
    LgND="NO"
    Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month")
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
  }
  mtext("Month",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CPUE (kg/ km gillnet hour)",side=2,line=0,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  
  #Predict block based on emmeans (formerly lsmeans) considering log bias corr if required
  Pred.blockx=vector('list',length=N.species)
  names(Pred.blockx)=SPECIES.vec 
  Pred.blockx.daily=Pred.blockx
  system.time({
    for(s in 1:N.species)   
    {
      d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
      Pred.blockx[[s]]=pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="blockx",Pred.type="link")
      
      d=Stand.out.daily[[s]]$DATA
      Pred.blockx.daily[[s]]=pred.fun(MOD=Stand.out.daily[[s]]$res,biascor="YES",PRED="blockx",Pred.type="link")
      rm(d)
    }
})     #takes 1.5 minutes
  
  SCLr=2.75
  CxTxt=0.0001
  #CxTxt=0.7
  Paleta=c("grey85","grey35")
  #Paleta=c("lightskyblue1","dodgerblue4")
  fn.fig("Figure 3.Block effect",2000, 2400)    
  par(mfrow=c(4,2),mar=c(1,1.5,1.5,1),oma=c(2.5,2.5,.1,1.75),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
    #Monthly
    Mon.dat=Pred.blockx[[s]]
    Mon.dat$LAT=-as.numeric(substr(Mon.dat$blockx,1,2))
    Mon.dat$LONG=100+as.numeric(substr(Mon.dat$blockx,3,4))
    Plot.cpue.spatial(cpuedata=Mon.dat,scaler=SCLr,colPalet=Paleta,CxTxt=CxTxt)
    if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    if(s%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    axis(2,seq(-36,-26,2),rev(seq(26,36,2)),tck=-0.025,cex.axis=1.35)
    
    
    #Daily
    Daily.dat=Pred.blockx.daily[[s]]
    Daily.dat$LAT=-as.numeric(substr(Daily.dat$blockx,1,2))
    Daily.dat$LONG=100+as.numeric(substr(Daily.dat$blockx,3,4))
    Plot.cpue.spatial(cpuedata=Daily.dat,scaler=SCLr,colPalet=Paleta,CxTxt=CxTxt)
    if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
    if(s%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
  }
  mtext("Longitude (ºE)",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Latitude (ºS)",side=2,line=0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
}


#ACA
#   4.22.11 Show daily_nfish
if(Model.run=="First")
{
  Pred.daily_nfish=Pred
  system.time({
    for(s in 1:N.species)   
    {
      d=Stand.out.daily_nfish[[s]]$DATA
      Pred.daily_nfish[[s]]=pred.fun(MOD=Stand.out.daily_nfish[[s]]$res,biascor="NO",PRED="finyear",Pred.type="response")
      rm(d)
    }
  })     #takes 3 minutes
  
  fn.fig("Figure 4.Annual_index_daily_nfish",2000, 2400)    
  par(mfrow=c(2,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
  for(s in 1:N.species)
  {
     LgND="NO"
    if(s==1)LgND="YES"
     Daily.dat=list(Standardised=Pred.daily_nfish[[s]])
    Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear")
    mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
}


#   4.22.12 Construct spatial standardised catch rates
#1. first glms to each specific zone
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



#   4.22.13 Export catch rates
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Index")

#4.22.13.1 zones combined without efficiency creep
for (s in 1:N.species)
{
  #Standardised
  a=subset(Pred[[s]],select=c(finyear,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly_no.creep.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily[[s]],select=c(finyear,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily_no.creep.csv",sep=""),row.names=F) 

  rm(a)
}

#4.22.13.2 zones combined with efficiency creep
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
  a=subset(Pred.creep[[s]],select=c(finyear,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.csv",sep=""),row.names=F) 
  
  a=subset(Pred.daily.creep[[s]],select=c(finyear,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.csv",sep=""),row.names=F) 
  

  #Unstandardised
  write.csv(List.foly.nom.creep[[s]]$Foly,paste(SPECIES.vec[s],".annual.folly.csv",sep=""),row.names=F)  
  
  a=subset(Nominl.creep[[s]],select=c(season,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.nominal.monthly.csv",sep=""),row.names=F)    
  
  a=subset(Nominl.daily.creep[[s]],select=c(season,response,lower.CL,upper.CL))   
  names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
  write.csv(a,paste(SPECIES.vec[s],".annual.nominal.daily.csv",sep=""),row.names=F) 

  rm(a)
}

#4.22.13.3 by zones without efficiency creep
for (s in 1:N.species)
{
  Zn=names(Pred.zone[[s]])
  for(z in 1:length(Zn))
  {
      #Standardised
    a=subset(Pred.zone[[s]][[z]],select=c(finyear,response,lower.CL,upper.CL))   
    names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
  
    a=subset(Pred.daily.zone[[s]][[z]],select=c(finyear,response,lower.CL,upper.CL))   
    names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.",Zn[z],"_no.creep.csv",sep=""),row.names=F) 
  
    rm(a)
  }
}

#4.22.13.4 by zones with efficiency creep
for (s in 1:N.species)
{
  Zn=names(Pred.zone[[s]])
  for(z in 1:length(Zn))
  {
    #Standardised
    a=subset(Pred.zone.creep[[s]][[z]],select=c(finyear,response,lower.CL,upper.CL))   
    names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.monthly.",Zn[z],".csv",sep=""),row.names=F) 
    
    a=subset(Pred.daily.zone.creep[[s]][[z]],select=c(finyear,response,lower.CL,upper.CL))   
    names(a)=c("Finyear","Mean","LOW.CI","UP.CI")
    write.csv(a,paste(SPECIES.vec[s],".annual.abundance.basecase.daily.",Zn[z],".csv",sep=""),row.names=F) 
    
    rm(a)
  }
}