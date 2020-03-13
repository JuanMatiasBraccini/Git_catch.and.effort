#--------- CPUE STANDARDISATIONS OF DIFFERENT DATASETS ---------#

#NOTE:  This script standardises the catch and effort data for the 4 commercial shark species,

#       To update Mean Freo Sealevel each year, extract csv data from "http://uhslc.soest.hawaii.edu/data/download/fd"
#       and run the script "Get.Freo.R" to get the mean monthly value from the daily records

# See Haddon's workshop, Haddon/.../data_exp_and_cpue_stand

#Index:  #----1. DATA SECTION-----#  
#           1.1 Import data
#           1.2 Control what parts of script are activated
#	            1.2.1 Data controls
#	            1.2.2 Procedure controls
#	            1.2.3 Reporting controls

#----2. PARAMETERS SECTION-----#

#----3. FUNCTIONS SECTION-----#

#----4. PROCEDURE SECTION-----#
#           4.1. Add SOI and Freo sea level data
#           4.2 Deal with zone1-zone2 Boundary blocks to a zone
#           4.3 Extract number of vessels per species range
#           4.4 Data fixes
#           4.5 Create useful vars
#           4.6 Remove NA effort
#           4.7 Explore proportion of dusky and copper shark
#           #4.9 Define indicative vessels and blocks 
#           4.10 Put data into a list
#           4.11 Compare nominal all records VS 'good reporters' only
#           4.12 Construct wide database for analysis
#           4.13 Compare wide vs long (raw) data sets
#           4.14 Corroborate effective area 
#           4.15 Identify targeting behaviour
#           4.16 Table of sensitivity scenarios
#           4.17 Compute foly index for exporting
#           4.18 Compute nominal for exporting
#           4.19 Input data tables
#           4.20 Plot cpue by year and month
#           4.21 Check outliers in catch and effort for removing nonsense values
#           4.22 Construct index of abundance
#             4.22.1 Explore data sets
#             4.22.2 Combine vessels in categories
#             4.22.3  Define model structure with Yr-Blk interactions
#             4.22.4 Write out best model
#             4.22.5 Get blocks from each time series
#             4.22.6 Show effort creep
#             4.22.7 Standardise catches
#             4.22.8 Get CI thru MC simulations
#             4.22.9 Base Case model fit diagnotics
#             4.22.10 Show Term Coefficients 
#             4.22.11 Records used
#             4.22.12 Plot data gaps 
#             4.22.13 Calculate relative standardised catch for MC simulations
#             4.22.14 Calculate relative standardised catch for sensitivity tests
#             4.22.15 Create index without effort creep
#             4.22.16 Plot indices
#             4.22.17 sensitivity tests
#             4.22.18 Influence Plots 
#             4.22.19 Compare creep and non-creep
#             4.22.20 Export relative index 
#             4.22.21  AMM actions


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


options(stringsAsFactors = FALSE)   #avoid characters being read in as factors


setwd("C:/Matias/Analyses/SOURCE_SCRIPTS")
source("Delta_lognormal.R")
source("Delta_gamma.R")
#source("Bootstrap_Delta_Methods.R")
source("Compare.error.structure.R")
source("Deviance.explained.R")
source("Sorting.objects.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/MS.Office.outputs.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/fn.fig.R")




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


#SOFAR effective cpue, plot and export
#hndLL="C:/Matias/Data/Catch and Effort/Effective cpue/"
# Gummy.Effective.SOFAR=read.csv(paste(hndLL,"Gummy.Effective.SOFAR.csv",sep=""))
# Dusky.Effective.SOFAR=read.csv(paste(hndLL,"Dusky.Effective.SOFAR.csv",sep=""))
# Whiskery.Effective.SOFAR=read.csv(paste(hndLL,"Whiskery.Effective.SOFAR.csv",sep=""))
# Sandbar.Effective.SOFAR=read.csv(paste(hndLL,"Sandbar.Effective.SOFAR.csv",sep=""))


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

#Control type of analysis done
ANALYSIS="Simple"     #no imputation
#ANALYSIS="Complex"   #imputation, bootstrapping, etc


#Control type of model run
#Model.run="First"    # for first time of doing standardisation. This allows selection of 
#     best model and sensitivity tests
Model.run="Standard"  # for running standardisations in subsequent years


#Control if including 'catch of other' and environmental covariates in model
Add.cov="NO"
#Add.cov="YES"


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

#Control whether to test model structure or not
#note: turn on only if wanting to define model structure and error
if(Model.run=="First")  Def.mod.Str.="YES"  
if(Model.run=="Standard") Def.mod.Str.="NO"

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

#Control if fitting models to dataset by zone or combined data set
#fit.to.by.zone="NO"     #fit one model to all data and predict blocks within zone
fit.to.by.zone="YES"   #fit separate model to different zones


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


if(ANALYSIS=="Complex")
{  
  #Control type of business rule for imputing missing last year of data
  Second.rule="combined" #use r if negative slope, use slope if positive
  #Second.rule="pop.growth.rate"       #use pop growth rate to project from last year (use average last three years)
  #Second.rule="average.last.three"  #use average of last 3 years (used given low number of observations in last year)
  #Second.rule="last.ob.carried.forward" 
  #Second.rule="linear.extrapolation"  #linear extrapolation
  
  impute.aver=3 #number of years used to impute mean for recent years
  impute.YRs.back=10  #number of years back to impute linear extrapolation
  
  
  #Control if extracting SE directly from glm 
  Extract.SE="NO"
  #Control if binomial part is set to same new data as positive part for predictions
  Bi.same.Pos="YES"  
  
  #Control if using same new data for all vessels as for base case
  All.ves.dat.same.BC.dat="NO"  
  #All.ves.dat.same.BC.dat="YES"
  
  #control if using same new data for predicting effort creep
  Creep.dat.same.BC.dat="YES"
  #Creep.dat.same.BC.dat="NO"
  
  #Control if fitting a yr-block interaction for sandbars
  Interaction.sandbar="NO"  
  
  #control if further exploring sandbar catch rates
  Do.sandbar.explore="NO"
  
  #Control new data factors
  Sel.lev='weighted.aver'
  #Sel.lev='most.common'
  
  #Control if standardised cpue are expressed in relative terms or not
  #RELATIVE="YES"
  RELATIVE="NO"
  
  #Control standardised cpue relative to what (note that relative to 1st year reduces uncertainty, creates a funnel)
  if(RELATIVE=="YES")
  {
    relative.to="Mean.1"
    #relative.to="First.yr"
  }
  
  #Control if plotting effort creep together
  compare.creep.sep="NO"
  
  
}


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

fn.plot.folly=function(a,spec,CPUE.threshold,CPUE.threshold1) #function for plotting folly example
{
  a=subset(a,SPECIES==spec)
  
  a=a[!duplicated(a$Same.return),]
  a$cpue=a$LIVEWT.c/a$Km.Gillnet.Days.c
  YRs=sort(unique(a$FINYEAR))
  YR=1:length(YRs)
  fn.plot=function(dat,MAIN)
  {
    a1=aggregate(Km.Gillnet.Days.c~FINYEAR,dat,sum)
    a2=aggregate(LIVEWT.c~FINYEAR,dat,sum)   
    a3=aggregate(cpue~FINYEAR,dat,mean)
    a4=aggregate(cpue~FINYEAR,dat,sd)
    a4[,2]=2*(a4[,2]/(table(dat$FINYEAR))^0.5) #CI
    LOW1=a3[,2]-a4[,2]
    UP1=a3[,2]+a4[,2]   
    YMAX=max(c(UP1,a2[,2]/a1[,2]))
    Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec <- c(LOW1, tail(UP1, 1), rev(UP1), LOW1[1])
    plot(YR,a3[,2],type='l',xlab="",ylab="",ylim=c(0,YMAX), main=MAIN,xaxt='n',cex.main=1.15,cex.axis=1.15)
    polygon(Year.Vec, Biom.Vec, col = "grey90", border = "grey90")
    lines(YR,a3[,2],type='l',lwd=2)
    lines(YR,a2[,2]/a1[,2],col=3,lwd=2)
    axis(1,YR,F,tck=-0.02)
    axis(1,YR[seq(YR[1],YR[length(YR)],5)],YRs[seq(YR[1],YR[length(YR)],5)],tck=-0.06,cex.axis=1.25)
  }
  fn.plot(a,"All records")
  
  #Subset by cpue threshold
  d=subset(a,cpue<CPUE.threshold)
  fn.plot(d,paste("Records with cpue<",CPUE.threshold,sep=""))
  
  #Subset by cpue threshold1
  d=subset(a,cpue<CPUE.threshold1)
  fn.plot(d,paste("Records with cpue<",CPUE.threshold1,sep=""))  
  if(spec==17003) legend("topright",c("mean approach","sum approach"),bty='n',lty=1,col=c(1,3),lwd=2,cex=1.45)
}

fn.see.month=function(d)
{
  d=subset(d,Catch.Target>0)
  d$cpue=d$Catch.Target/d$Km.Gillnet.Days.c
  
  Mon.ktch=aggregate(Catch.Target~MONTH,d,mean)
  Mon.ktch.sd=aggregate(Catch.Target~MONTH,d,sd)
  Mon.ktch.N=aggregate(Catch.Target~MONTH,d,length)
  Mon.ktch.SE=Mon.ktch.sd$Catch.Target/(Mon.ktch.N$Catch.Target)^0.5
  
  with(Mon.ktch,plot(MONTH,Catch.Target,ylim=c(0,max(Catch.Target+Mon.ktch.SE)),
                     ylab="",xlab="",pch=19,col=1,cex.lab=1.25))
  with(Mon.ktch,segments(MONTH,Catch.Target+Mon.ktch.SE,MONTH, Catch.Target-Mon.ktch.SE))
  
  par(new=T)
  Mon.cpue=aggregate(cpue~MONTH,d,mean)
  Mon.cpue.sd=aggregate(cpue~MONTH,d,sd)
  Mon.cpue.N=aggregate(cpue~MONTH,d,length)
  Mon.cpue.SE=Mon.cpue.sd$cpue/(Mon.cpue.N$cpue)^0.5
  with(Mon.cpue,plot(MONTH,cpue,ylim=c(0,max(cpue+Mon.cpue.SE)),pch=19,col=2,
                     ylab="",xlab="",axes=F,ann='F'))
  with(Mon.cpue,segments(MONTH,cpue+Mon.cpue.SE,MONTH, cpue-Mon.cpue.SE,col=2))
  axis(4,col=2,col.ticks=2,col.axis=2)
  legend("bottomleft",SPECIES.vec[i],bty='n',cex=1.5)
  abline(v=7,col=3,lwd=2)
}

fn.n.vessels=function(dat,dat1)
{
  return(length(unique(c(dat$VESSEL,dat1$VESSEL))))
}

#functions for removing blocks with low catch
fn.n.blk.caught=function(Dat,Dat1,SPEC)  #total number of bloks where shark caught
{
  Dat=subset(Dat,SPECIES%in%SPEC & LIVEWT.c>0)
  Dat1=subset(Dat1,SPECIES%in%SPEC & LIVEWT.c>0)
  return(unique(c(Dat$BLOCKX,Dat1$BLOCKX)))
}
fn.drop.blk=function(Dat,Tops) 
{
  ALL=unique(Dat$BLOCKX)
  ALL=ALL[-match(Tops,ALL)]
}
fn.check.blok=function(Dat,SPEC,YR,Threshold)
{
  Dat=subset(Dat,FINYEAR%in%YR & SPECIES==SPEC)
  if(nrow(Dat)>0)
  {
    Dat$cpue=with(Dat,LIVEWT.c/Km.Gillnet.Hours.c)
    CPUe=aggregate(cpue~BLOCKX,Dat,sum,na.rm=T)
    CPUe=CPUe[order(-CPUe$cpue),]
    CPUe$Cum=cumsum(CPUe$cpue)/sum(CPUe$cpue)
    id=which(CPUe$Cum<=Threshold)
    TopBlok=CPUe$BLOCKX[id]
    
    return(list(cpue=CPUe,TopBlok=TopBlok))
  }
}

fn.cum.ktch=function(dat,SP,Threshold,MinObs)
{
  Dat=subset(dat,SPECIES==SP & LIVEWT.c>0 & Reporter=="good",select=c(BLOCKX,LIVEWT.c))
  TAB=table(Dat$BLOCKX)
  TAB=TAB[TAB>MinObs]
  Dat=subset(Dat,BLOCKX%in%names(TAB))
  KTCH=aggregate(LIVEWT.c~BLOCKX,Dat,sum,na.rm=T)
  KTCH=KTCH[order(-KTCH$LIVEWT.c),]
  KTCH$Cum=round(cumsum(KTCH$LIVEWT.c)/sum(KTCH$LIVEWT.c),2)
  id=which(KTCH$Cum<=Threshold)
  TopBlok=KTCH$BLOCKX[id]
  all.blks=unique(dat$BLOCKX)
  Dropped.blks=subset(all.blks,!all.blks%in%TopBlok)
  
  return(list(TopBloks=TopBlok,Dropped.blks=Dropped.blks))
}
fn.cum.ktch.daily=function(dat,SP,Threshold,MinObs)
{
  Dat=subset(dat,SPECIES==SP & LIVEWT.c>0 & Reporter=="good",select=c(block10,LIVEWT.c))
  TAB=table(Dat$block10)
  TAB=TAB[TAB>MinObs]
  Dat=subset(Dat,block10%in%names(TAB))
  KTCH=aggregate(LIVEWT.c~block10,Dat,sum,na.rm=T)
  KTCH=KTCH[order(-KTCH$LIVEWT.c),]
  KTCH$Cum=round(cumsum(KTCH$LIVEWT.c)/sum(KTCH$LIVEWT.c),2)
  id=which(KTCH$Cum<=Threshold)
  TopBlok=KTCH$block10[id]
  all.blks=unique(dat$block10)
  Dropped.blks=subset(all.blks,!all.blks%in%TopBlok)
  return(list(TopBloks=TopBlok,Dropped.blks=Dropped.blks))
}
fn.show.blk=function(dat,dat1,show.dropped)  #Show removed blocks
{
  dat=sort(dat)
  dat1=sort(dat1)
  dat1=subset(dat1,!dat1%in%c(33151,28132,28142,29142))
  
  LAT.dropped=sapply(dat, function(x) -as.numeric(substr(x, 1, 2)))
  LONG.dropped=sapply(dat, function(x) 100+as.numeric(substr(x, 3, 4)))
  
  LAT.kept=sapply(dat1, function(x) -as.numeric(substr(x, 1, 2)))
  LONG.kept=sapply(dat1, function(x) 100+as.numeric(substr(x, 3, 4)))
  
  plot(X,Y,ylab='LATITUDE',xlab="LONGITUDE",col="transparent",cex.lab=1.5,cex.axis=1.25)
  #main=paste("N blocks=",length(dat1))
  
  #kept
  for(e in 1:length(LAT.kept))
  {
    dd.y=c(LAT.kept[e]-1,LAT.kept[e]-1,LAT.kept[e],LAT.kept[e])
    dd.x=c(LONG.kept[e],LONG.kept[e]+1,LONG.kept[e]+1,LONG.kept[e])
    polygon(dd.x,dd.y,col=rgb(0, 0, 1,0.25), border=1)
    text(LONG.kept[e]+0.5,LAT.kept[e]-0.5,dat1[e],cex=1.2,col=1,srt=45)
  }
  
  if(show.dropped=="YES")
  {
    #dropped
    for(e in 1:length(LAT.dropped))
    {
      dd.y=c(LAT.dropped[e]-1,LAT.dropped[e]-1,LAT.dropped[e],LAT.dropped[e])
      dd.x=c(LONG.dropped[e],LONG.dropped[e]+1,LONG.dropped[e]+1,LONG.dropped[e])
      polygon(dd.x,dd.y,col=rgb(1, 0, 0,0.25), border=NA)
      text(LONG.dropped[e]+0.5,LAT.dropped[e]-0.5,dat[e],cex=0.55,col=1)
    }
  }
  
}

#function for spatial expansion of effort
Effect.area=function(a)
{
  yr=sort(unique(a$FINYEAR))
  n=length(yr)
  
  agg=aggregate(Km.Gillnet.Days.c~FINYEAR+BLOCKX+LAT+LONG,a,sum)
  
  for(i in 1:n)
  {
    b=subset(agg,FINYEAR==yr[i])
    z=b$Km.Gillnet.Days.c/max(agg$Km.Gillnet.Days.c)
    plot(b$LONG,b$LAT,cex=z*2,main=yr[i],ylab="",xlab="",pch=19,ylim=c(-36,-26),
         xlim=c(113,129),cex.axis=0.8,cex.main=.85,col="steelblue4")
    legend("top",paste(round(c(max(agg$Km.Gillnet.Days.c),max(agg$Km.Gillnet.Days.c)/2)),"km.gn.d"),
           pch=19,pt.cex=c(1*2,1/2),bty='n',col="steelblue4",cex=.85)
    
  }
}


#function for cluster analysis of daily data
fn.targeting=function(dat,SP)
{
  dat=subset(dat,!SPECIES==999999,select=c(date,VESSEL,block10,LAT,LONG,SPECIES,LIVEWT.c))
  IDVar=c("date","VESSEL","block10","LAT","LONG")
  dat1=reshape(dat,v.names="LIVEWT.c",idvar=IDVar,timevar="SPECIES",direction="wide")
  names(dat1)[(length(IDVar)+1):ncol(dat1)]=substr(names(dat1)[(length(IDVar)+1):ncol(dat1)],10,30)
  Total=rowSums(dat1[(length(IDVar)+1):ncol(dat1)],na.rm=T)
  
  #1. Proportion of shot's catch and grouping of non-target sp
  
  Non.Targt=dat1[,-match(c(IDVar,as.character(Targt)),names(dat1))]
  Non.Targt=rowSums(Non.Targt,na.rm=T)
  dat1=dat1[,match(c(IDVar,as.character(Targt)),names(dat1))]
  dat1$Other=Non.Targt
  MATRIX=dat1[,(length(IDVar)+1):ncol(dat1)]
  MATRIX$Total=Total
  MATRIX=as.matrix(MATRIX[,match(c(Targt,"Other"),names(MATRIX))]/MATRIX$Total)
  MATRIX[is.na(MATRIX)]=0
  
  #2. calculate distance matrix using Euclidean distance 
  #note: this calculates the dissimilarity between samples, the closer the samples, the more similar
  #rownames(MATRIX)=with(dat1,paste(date,".",VESSEL,".",block10,sep=""))
  IID=match(as.character(SP),colnames(MATRIX))
  rownames(MATRIX)=round(MATRIX[,IID],3)   #name matrix with proportion of target species
  
  #standardize variables (scale matrix by centering the columns)
  MATRIX_s=scale(MATRIX)
  
  mydata=MATRIX_s
  
  
  #3. Cluster the distance matrix
  
  # Determine number of clusters
  #   wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  #   for (i in 2:length(wss)) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
  #   plot(1:length(wss), wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
  #   
  #Ward's method
  Dist.Mat=dist(mydata,method='euclidean')
  hc <- hclust(Dist.Mat, method="ward.D") 
  plot(hc)
  groups <- cutree(hc, k=5) # cut tree into 5 clusters
  rect.hclust(hc, k=5, border="red") # draw dendogram with red borders around the 5 clusters 
  
  #   #4.  Ward Hierarchical Clustering with Bootstrapped p values
  #   fit <- pvclust(mydata, method.hclust="ward.D", method.dist="euclidean")
  #   plot(fit) # dendogram with p values
  #   # add rectangles around groups highly supported by the data
  #   pvrect(fit, alpha=.95)
  
  #   #5.  K-Means Cluster Analysis
  #   fit <- kmeans(mydata, 5) # 5 cluster solution  
  #   aggregate(mydata,by=list(fit$cluster),FUN=mean)  # get cluster means
  #   # append cluster assignment
  #   mydata <- data.frame(mydata, fit$cluster)
  
  #6. Plotting clusters
  
  #6.1 K-Means Clustering with 5 clusters
  # fit <- kmeans(mydata, 5)
  
  # Cluster Plot against 1st 2 principal components
  
  # vary parameters for most readable graph
  #  clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,labels=2, lines=0)
  
  # Centroid Plot against 1st 2 discriminant functions
  #  plotcluster(mydata, fit$cluster)
  
  #   #4. non-metric MDS plot
  #   fit <- cmdscale(Dist.Mat,eig=TRUE, k=2) # k is the number of dim
  #   x <- fit$points[,1]
  #   y <- fit$points[,2]
  #   plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  #        main="Metric  MDS",  type="n")
  #   text(x, y, labels = row.names(MATRIX), cex=.7)
}

#check variability in vessels cpue and month cpue by year
fn.vessl.cpue.yr=function(dat)  
{
  Yr=sort(unique(dat$FINYEAR)) 
  dat$cpue=dat$Catch.Target/dat$Km.Gillnet.Days.c
  tab1=aggregate(cpue~FINYEAR+VESSEL,dat,mean)
  tab11=reshape(tab1,v.names = "cpue", idvar = "FINYEAR",timevar = "VESSEL", direction = "wide")
  plot(tab11[,2],ylim=c(0,max(tab11[,2:ncol(tab11)],na.rm=T)),type='l')
  for(s in 3:ncol(tab11))lines(tab11[,s],col=runif(1,1,12))
  
  tab1=aggregate(cpue~FINYEAR+MONTH,dat,mean)
  tab11=reshape(tab1,v.names = "cpue", idvar = "FINYEAR",timevar = "MONTH", direction = "wide")
  plot(tab11[,2],ylim=c(0,max(tab11[,2:ncol(tab11)],na.rm=T)),type='l')
  for(s in 3:ncol(tab11))lines(tab11[,s],col=runif(1,1,12))
}

#Aggregate predictions by bins
fn.bin.pred=function(dat,Int)   
{
  dat$BINS=cut(dat$Pred, seq(min(dat$Pred),max(dat$Pred),Int), include.lowest=TRUE)
  Agg.pred=aggregate(Pred~BINS,dat,mean)
  Agg.Obs=aggregate(log.cpue~BINS,dat,mean)
  
  #plot
  par(mfcol=c(2,1))
  LIMS=c(min(c(dat$log.cpue,dat$Pred)),max(c(dat$log.cpue,dat$Pred)))
  plot(dat$log.cpue,dat$Pred,ylim=LIMS,xlim=LIMS)
  lines(dat$log.cpue,dat$log.cpue,col=3)
  
  plot(Agg.Obs$log.cpue,Agg.Obs$log.cpue,type='l',lty=2,lwd=2,col="grey",ylab="",xlab="",cex.axis=1.5)
  points(Agg.pred$Pred,Agg.Obs$log.cpue,pch=19,cex=1.75)
  return(dat)
}

fn.Tabl1=function(dat)
{
  YR=sort(unique(dat$YEAR.c))
  YR=paste(YR[1],"-",YR[length(YR)],sep="")    
  BLc=length(unique(dat$BLOCKX))
  MN=sort(unique(dat$MONTH))
  MN=paste(MN[1],"-",MN[length(MN)],sep="")  
  Vsl=length(unique(dat$VESSEL))
  return(data.frame(Financial_years=YR,Spatial_blk=BLc,Month=MN,Vessel=Vsl))         
}
fn.Tabl1.daily=function(dat)
{
  YR=sort(unique(dat$YEAR.c))
  YR=paste(YR[1],"-",YR[length(YR)],sep="")    
  BLc=length(unique(dat$block10))
  MN=sort(unique(dat$MONTH))
  MN=paste(MN[1],"-",MN[length(MN)],sep="")  
  Vsl=length(unique(dat$VESSEL))
  return(data.frame(Financial_years=YR,Spatial_blk=BLc,Month=MN,Vessel=Vsl))         
}


fn.AICc=function(DAT)    #function for calculation AICc
{
  Sample.size=dim(model.matrix(DAT))[1]
  Num.pars=dim(model.matrix(DAT))[2]
  LoGLike=logLik(DAT)[1]
  return(AIC.c=-2*LoGLike+(2*Num.pars*(Sample.size/(Sample.size-Num.pars-1))))
} 
fn.Dev.Exp=function(DAT)Dsquared(DAT,adjust=T)  #function for calculating adjusted deviance explained
fit.model.to.inspect=function(DATA,VES.No.K,Formula.cpue.pos,Vars,MIN.Wght,Int) #function for fitting model for checking predictions
{
  #1. drop vessels that never reported catch of target species   
  DATA=subset(DATA,!(VESSEL%in%VES.No.K))
  #DATA=subset(DATA,!is.na(new.level.vessel))
  DATA=subset(DATA,Catch.Target>=MIN.Wght)  #remove records with nonsense catch value
  
  #2. put data in proper shape
  DataFile=DATA[,match(Vars,names(DATA))]
  
  #drop levels not occurring in data
  for(f in 1:ncol(DataFile))if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
  
  #convert to factor and sort levels
  DataFile$FINYEAR=as.factor(DataFile$FINYEAR)
  DataFile$BLOCKX=as.factor(DataFile$BLOCKX)
  DataFile$MONTH=as.factor(DataFile$MONTH)
  DataFile$VESSEL=as.factor(DataFile$VESSEL)
  #DataFile$new.level.vessel=as.factor(DataFile$new.level.vessel)
  
  #log effort
  DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c)
  
  
  # log catch other
  DataFile$log.Catch.Gummy=log(DataFile$Catch.Gummy+0.000000001)
  DataFile$log.Catch.Whiskery=log(DataFile$Catch.Whiskery+0.000000001)
  DataFile$log.Catch.Dusky=log(DataFile$Catch.Dusky+0.000000001)
  DataFile$log.Catch.Sandbar=log(DataFile$Catch.Sandbar+0.000000001)
  DataFile$log.Catch.Scalefish=log(DataFile$Catch.Scalefish+0.000000001)
  DataFile$log.Freo=log(DataFile$Freo)
  DataFile$log.Freo.Lag6=log(DataFile$Freo.Lag6)
  DataFile$log.Freo.Lag12=log(DataFile$Freo.Lag12)
  
  # Create positive (non-zero) dataset
  Id.Ktch=match("Catch.Target",names(DataFile))
  dat <- DataFile[DataFile[,Id.Ktch]>0,]  
  dat$log.Catch=log(dat$Catch.Target)  
  dat$cpue=dat$Catch.Target/dat$Km.Gillnet.Days.c
  dat$log.cpue=log(dat$cpue)
  
  #Lognormal GLM
  MODEL <- glm(Formula.cpue.pos, data=dat, family=gaussian, maxit=500)
  
  par(mfcol=c(3,1))
  Pos.Diag.fn(MODEL,SPECIES.vec[i])
  print(paste("AIC",fn.AICc(MODEL)))
  print(paste("Deviance explained",fn.Dev.Exp(MODEL)))
  dat$Pred=predict(MODEL,type="response")
  dat=fn.bin.pred(dat,Int)
  return(dat)    
}


cfac=function(x,breaks=NULL)  #function for converting continuous var to factor
{
  if(is.null(breaks)) breaks=unique(quantile(x,probs = seq(0, 1, 0.1)))
  x=cut(x,breaks,include.lowest=T,right=F)
  levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                                                 c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
  return(x)
}
clog=function(x) log(x+0.05)   #function for applying log

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
  Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+
        FINYEAR+Same.return+MONTH+BLOCKX+SHOTS.c+BDAYS.c+HOURS.c+NETLEN.c,Effort.data1,max)
  # Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+      #redundant
  #                         FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,sum)
  
  Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data1,max)
  #Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data.inv,sum)
  
  Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return,Effort.data1,max)
  #Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return,Effort.data.no.creep,sum)
  
  Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return",all.x=T)
  Effort.data=merge(Effort.data,Effort.data.no.creep,by="Same.return",all.x=T)
  
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
    Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+
                            FINYEAR+Same.return.SNo+MONTH+BLOCKX+block10+shots.c+hours.c+netlen.c,Effort.data1,max)
    Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return.SNo,Effort.data1,max)
    Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return.SNo,Effort.data1,max)
    
    
    #aggregate at TSNo if required
    if(Aggregtn=="TSNo")
    {
      Effort.data$TSNo=word(Effort.data$Same.return.SNo,3)
      Effort.data.inv$TSNo=word(Effort.data.inv$Same.return.SNo,3)
      Effort.data.no.creep$TSNo=word(Effort.data.no.creep$Same.return.SNo,3)
      
      Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+
                              FINYEAR+TSNo+MONTH+BLOCKX,Effort.data,sum)
      Effort.data.inv=aggregate(Km.Gillnet.Days.inv~TSNo,Effort.data.inv,sum)
      Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~TSNo,Effort.data.no.creep,sum)
    }
    
    #merge as appropriate
    if(Aggregtn=="SNo")
    {
      Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return.SNo",all.x=T)
      Effort.data=merge(Effort.data,Effort.data.no.creep,by="Same.return.SNo",all.x=T) 
    }
    if(Aggregtn=="TSNo")
    {
      Effort.data=merge(Effort.data,Effort.data.inv,by="TSNo",all.x=T)
      Effort.data=merge(Effort.data,Effort.data.no.creep,by="TSNo",all.x=T) 
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


qual.levl.fun=function(DAT,Prop.ktch.exp)
{
  dummy=aggregate(Catch.Target~FINYEAR,DAT,sum)
  Yrs=subset(dummy,Catch.Target>0,select=FINYEAR)
  DAT=subset(DAT,FINYEAR%in%Yrs$FINYEAR)
  ZNs=sort(unique(DAT$zone))
  Yrs=sort(unique(DAT$FINYEAR))
  STR.zn=vector('list',length(ZNs))
  STR.yrs=vector('list',length(Yrs))
  for(zn in 1:length(ZNs))
  {
    Dat=subset(DAT,zone==ZNs[zn])
    for(yrs in 1:length(Yrs))
    {
      dat=subset(Dat,FINYEAR==Yrs[yrs])
      dat$prop.target=dat$Catch.Target/dat$Catch.Total
      dat=dat[order(-dat$prop.target),]
      dat$CumSum=cumsum(dat$Catch.Target)
      dat$PropTarget=dat$CumSum/sum(dat$Catch.Target)
      dat$Qualif=dat$prop.target[which(dat$PropTarget>=Prop.ktch.exp)[1]]
      dat=subset(dat,prop.target>=unique(dat$Qualif))   #remove records outside qualification level
      STR.yrs[[yrs]]=dat
    }
    STR.zn[[zn]]=do.call(rbind,STR.yrs)
  }
  return(do.call(rbind,STR.zn))
}
plot.qual.levl.fun=function(DAT,Prop.ktch.exp)
{
  ZNs=sort(unique(DAT$zone))
  Yrs=sort(unique(DAT$FINYEAR))
  nyrs=length(Yrs)
  CL=2:4
  plot(1:nyrs,ylim=c(0,max(DAT$Qualif)),col='transparent',xlab='',xaxt='n',ylab="")
  legend('bottomleft',paste(names(QUALIF.LEVL)[i],"shark"),bty='n',cex=1.5)
  if(i==1) legend('top',ZNs,bty='n',cex=1.25,lty=1,lwd=2,col=CL)
  axis(1,1:nyrs,F,tck=-0.02)
  axis(1,seq(1,nyrs,5),Yrs[seq(1,nyrs,5)],tck=-0.04)
  
  for(q in 1:length(ZNs))
  {
    aa=subset(DAT,zone==ZNs[q])
    aa=aa[!duplicated(aa$FINYEAR),]
    lines(1:length(aa$FINYEAR),aa$Qualif,lwd=2,col=CL[q])
  }
}

fn.indicative=function(DAT,Min.yrs,n.records,second.criteria,KTCH.threshold)  #function for keeping indicative vessels
{
  #1st Criteria, vessels with at least n.records
  DAT.1=subset(DAT,Catch.Target>0)
  Vessel.n=with(DAT.1,table(VESSEL))
  Vessel.n=Vessel.n[Vessel.n>=n.records]
  
  #2nd Criteria, reporting catch of the target specie for >= Min.yrs
  DAT.1=subset(DAT.1,VESSEL%in%names(Vessel.n))
  Vessel.n.yrs=with(DAT.1,table(VESSEL,FINYEAR))
  Vessel.n.yrs=ifelse(Vessel.n.yrs>=1,1,0)
  Vessel.n.yrs=sort(rowSums(Vessel.n.yrs))
  Indic.ves=names(Vessel.n.yrs[which(Vessel.n.yrs>=Min.yrs)])
  
  
  #3rd Criteria
  DATA.indic=subset(DAT,VESSEL%in%Indic.ves)
  TABLE12=aggregate(Catch.Target~VESSEL+FINYEAR,data=DATA.indic,sum)
  a=reshape(TABLE12,v.names="Catch.Target",idvar="FINYEAR",timevar="VESSEL",direction='wide')
  FinYrs=a$FINYEAR
  a[a>0]=1
  a$FINYEAR=FinYrs    
  a$ROWSUM=rowSums(a[,2:ncol(a)],na.rm=T)
  a=subset(a,ROWSUM>1)  #drop years when only 1 or no vessel reporting catch
  dropped.yrs=FinYrs[which(!FinYrs%in%a$FINYEAR)]
  if(length(dropped.yrs)==0) dropped.yrs=NA
  FinYrs=a$FINYEAR  
  DATA.indic=subset(DATA.indic,FINYEAR%in%a$FINYEAR)
  
  #annual catch within top KTCH.threshold
  if(second.criteria=="top.percent")
  {
    Uni.yr=sort(unique(DATA.indic$FINYEAR))
    vess.store=vector("list",length(Uni.yr))
    names(vess.store)=Uni.yr
    DATS=vess.store
    for (q in 1: length(Uni.yr))
    {
      TABLE12=aggregate(Catch.Target~VESSEL,data=subset(DATA.indic,FINYEAR==Uni.yr[q]),sum)
      TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
      TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
      TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
      TABLE12$Delta=KTCH.threshold-TABLE12$PerCumCatch
      id=which(TABLE12$Delta<0)[1]
      vess.store[[q]]=TABLE12$VESSEL[1:id]
      DATS[[q]]=subset(DATA.indic,VESSEL%in%vess.store[[q]] & FINYEAR ==Uni.yr[q])
    }
    DATA.indic=do.call(rbind,DATS)
  }
  
  #median target species totalcatch >= fishery median catch
  if(second.criteria=="median.fishery")
  {
    Threshold.median.catch=aggregate(Catch.Target~VESSEL,data=DAT.1,sum)
    Threshold.median.catch=median(Threshold.median.catch$Catch.Target)
    Vessel.median.catch=aggregate(Catch.Target~VESSEL,data=DATA.indic,sum)
    Indic.ves=subset(Vessel.median.catch,Catch.Target>=Threshold.median.catch)$VESSEL
    DATA.indic=subset(DATA.indic,VESSEL%in%Indic.ves)
  }
  
  
  #How much catch explained
  prop.ktc.explained=sum(DATA.indic$Catch.Target)/sum(DAT$Catch.Target)
  
  #vessels kepth
  Indic.ves=sort(unique(DATA.indic$VESSEL))
  Kept.ves=round(length(Indic.ves)/length(unique(DAT$VESSEL)),2)
  dropped.ves=1-Kept.ves
  
  return(list(Indic.ves=Indic.ves,vess.store=vess.store,
              Prop.Ves.kept=Kept.ves,Prop.Ves.dropped=dropped.ves,
              Prop.ktch.exp.=prop.ktc.explained,dropped.yrs=dropped.yrs))
}

Yrs.by.blk=function(DATS)   #function for tables
{ 
  DATS.no=subset(DATS,Catch.Target==0)
  DATS=subset(DATS,Catch.Target>0)
  Table1=table(DATS$FINYEAR,DATS$BLOCKX)
  Table2=table(DATS$FINYEAR,DATS$VESSEL)
  Table3=table(DATS$FINYEAR,DATS$MONTH)
  Table4=table(DATS$BLOCKX)
  Table5=table(DATS$VESSEL)
  
  Unic.ves.no.K=unique(DATS.no$VESSEL)
  id=which(!Unic.ves.no.K%in%unique(DATS$VESSEL))
  Unic.ves.no.K=Unic.ves.no.K[id]
  
  Unic.blk.no.K=unique(DATS.no$BLOCKX)
  id=which(!Unic.blk.no.K%in%unique(DATS$BLOCKX))
  Unic.blk.no.K=Unic.blk.no.K[id]
  
  
  return(list(Table1=Table1,Table2=Table2,Table3=Table3,Table4=Table4,Table5=Table5,
              Table6=Unic.blk.no.K,Table7=Unic.ves.no.K))
}

plt.ktch.yr_mn=function(DATA)  #function for plotting cpue by year and month
{
  this=aggregate((Catch.Target/Km.Gillnet.Days.c)~FINYEAR+MONTH,DATA,mean)
  Reshaped=(reshape(this,v.names = "(Catch.Target/Km.Gillnet.Days.c)", idvar = "FINYEAR",
                    timevar = "MONTH", direction = "wide"))
  FinYr=sort(as.character(Reshaped[,1]))
  N.yrs=length(FinYr)
  Blk.lab=sort(unique(DATA$MONTH))
  Blk=1:length(Blk.lab)
  
  numInt=10
  couleurs=rev(heat.colors(numInt))
  CPUEBreaks=quantile(unlist(Reshaped[,2:ncol(Reshaped)]),probs=seq(0,1,1/numInt),na.rm=T)
  numberLab=2
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  Reshaped=Reshaped[order(Reshaped[,1]),]
  Reshaped=Reshaped[,-1]
  Reshaped=as.matrix(Reshaped)
  
  par(mfcol=c(1,1),mai=c(0.8,1.0,.1,1.5),oma=c(0.01,.1,.1,.01),las=1)
  image(1:length(FinYr),Blk,z=Reshaped,axes = FALSE,ylab="Month",xlab="Financial year",
        col =couleurs,breaks=CPUEBreaks,cex.lab=1.75)
  axis(1,at=1:N.yrs,labels=F,tck=-0.01)
  axis(1,at=seq(1,N.yrs,by=5),labels=FinYr[seq(1,N.yrs,by=5)],tck=-0.02,cex.axis=1.05)
  
  axis(2,at=Blk,labels=F,tck=-0.01)
  axis(2,at=seq(1,length(Blk),by=1),labels=Blk.lab[seq(1,length(Blk),by=1)],tck=-0.02,cex.axis=1.05)
  color.legend(N.yrs+4,length(Blk.lab),N.yrs+6,length(Blk.lab)*.80,round(CPUEBreaks,1),rect.col=couleurs,
               gradient="y",col=colLeg,cex=0.95)
  box()
}

fn.plot.vess.selection=function(DAT,DAT1) #function for plotting number of vessels & catch for each data set
{
  ax.fn=function(YR)
  {
    Nn=length(YR)
    axis(1,1:Nn,F,tck=-0.02)
    axis(1,seq(1,Nn,5),YR[seq(1,Nn,5)],tck=-0.04,cex.axis=1.25)
  }
  
  LWD=2
  pt.fn=function(what,what1,what2)
  {
    NN=1:length(what)
    plot(NN,what,type='l',ylim=c(0,max(what)),ylab="",xlab="",xaxt="n",lwd=LWD,cex.axis=1.25)
    
    NN1=which(!names(what)%in%names(what1))
    if(length(NN1)>0)
    {
      crap=what[NN1]
      crap[crap>0]=NA
      what1=c(what1,crap)
      what1=what1[sort(names(what1))]
    }
    lines(NN,what1,lty=2,lwd=LWD)
    
    NN1=which(!names(what)%in%names(what2))
    if(length(NN1)>0)
    {
      crap=what[NN1]
      crap[crap>0]=NA
      what2=c(what2,crap)
      what2=what2[sort(names(what2))]
    }
    lines(NN,what2,col="grey70",lwd=LWD)
    ax.fn(YRs)
  }
  
  #number of vessels
  fn.ves=function(d)
  {
    N=with(d,table(VESSEL,FINYEAR))
    N[N>0]=1
    return(colSums(N)) 
  }
  N.ves.all=fn.ves(DAT)
  N.ves.all.pos.ktch=fn.ves(subset(DAT,Catch.Target>0))
  N.ves.ind=fn.ves(DAT1)
  YRs=names(N.ves.all)
  pt.fn(N.ves.all,N.ves.all.pos.ktch,N.ves.ind)
  if(i==2) mtext("Number of vessels                      ",2,cex=1.5,line=2,las=3)
  if(i==1) legend("topright",c("All vessels","Vessels with catch","Indicative vessels"),bty='n',
                  lwd=LWD,lty=c(1,2,1),col=c("black","black","grey70"),cex=1.3)
  
  #Total catch of species
  fn.Tot.Ktch=function(d)
  {
    N=aggregate(Catch.Target~FINYEAR,d,sum)
    a=N$Catch.Target/1000
    names(a)=N$FINYEAR
    return(a) 
  }
  Tot.ktch.all=fn.Tot.Ktch(DAT)
  Tot.ktch.all.pos.ktch=fn.Tot.Ktch(subset(DAT,Catch.Target>0))
  Tot.ktch.ind=fn.Tot.Ktch(DAT1)
  pt.fn(Tot.ktch.all,Tot.ktch.all.pos.ktch,Tot.ktch.ind)
  if(i==2) mtext("Total catch of species (tonnes)                 ",2,cex=1.5,line=2.5,las=3)
  return(Prop.catch.exp=Tot.ktch.ind/Tot.ktch.all)
}

check.cpue=function(DATA,NAME,BYDAY)   #function for checking cpue outliers
{
  DATA$Catch.Target=DATA$Catch.Target
  DATA$Km.Gillnet.Days.c=DATA$Km.Gillnet.Days.c
  
  if(BYDAY=="YES")
  {
    DATA$Catch.Target=DATA$Catch.Target/30
    DATA$Km.Gillnet.Days.c=DATA$Km.Gillnet.Days.c/30
  }
  
  DATA$CPUE=DATA$Catch.Target/DATA$Km.Gillnet.Days.c
  
  Ktc.q=quantile(DATA$Catch.Target,probs=seq(0,1,.1))
  Eff.q=quantile(DATA$Km.Gillnet.Days.c,probs=seq(0,1,.1))
  CPUE.q=quantile(DATA$CPUE,probs=seq(0,1,.1))
  par(mfcol=c(2,1),mai=c(.5,1.1,.3,.1),mgp=c(2.5,.5,0))
  boxplot((Catch.Target)~FINYEAR,DATA,ylab="KG",main=NAME)
  abline(h=Ktc.q[6],col=2)
  legend("topright",paste("median=",round(Ktc.q[6],3),"       "),text.col=2,cex=1.15,bty='n')
  
  boxplot(CPUE~FINYEAR,DATA,ylab="KG / km gn hr")
  abline(h=CPUE.q[6],col=2)
  legend("topright",paste("median=",round(CPUE.q[6],3),"       "),text.col=2,cex=1.15,bty='n')
  
  return(list(Ktc.q=Ktc.q,Eff.q=Eff.q,CPUE.q=CPUE.q))
}

plot0=function(DATA,MAX,Ymax,Xmax,SPECIES)  #function for plotting 0 catches
{
  DATA.no0=subset(DATA,Catch.Target>0)
  HIST=hist(DATA.no0$Catch.Target/1000,breaks=seq(0,MAX,by=1),plot=F)
  Prop=HIST$counts/nrow(DATA)
  Catch.0=round((nrow(DATA)-nrow(DATA.no0))/nrow(DATA),2)
  
  Prop.plus=sum(Prop[(Xmax+1):length(Prop)])
  Prop=Prop[2:(Xmax)]
  #  barplot(c(Prop,Prop.plus),ylim=c(0,Ymax),las=1,names.arg=c(DATA$breaks[2:(Xmax)],paste(">",Xmax-1,sep="")),xlim=c(0,Xmax+2))
  a=barplot(c(Prop,Prop.plus),ylim=c(0,Ymax),las=1,xlim=c(0,Xmax+2))
  axis(1,at=a[1:(length(Prop)+1)],labels=F,tck=-0.032,cex.axis=.5)
  if(SPECIES%in%c("Gummy shark","Sandbar shark"))
  {
    axis(1,at=a[1:length(Prop)],labels=1:length(Prop),tck=-0.032,cex.axis=.9)
    axis(1,at=a[length(Prop)+1],labels=paste(">",Xmax-1,sep=""),tck=-0.032,cex.axis=.875)
  }
  box()
  legend("bottomright",paste("Proportion 0 catch=",round(Catch.0,2)),bty='n',cex=.9)
  
}

Zero.catch.fun2=function(DATA) #function for plotting 0 catches
{
  DATA$Target=with(DATA,ifelse(Catch.Target>0,1,0))
  DATA$NoTarget=with(DATA,ifelse(Catch.Target==0,1,0))
  test=aggregate(cbind(Target,NoTarget)~FINYEAR,data=DATA,sum,na.rm=T)
  test$sum=test$Target+test$NoTarget
  proportion=test$Target/test$sum
  return(data.frame(FINYEAR=test$FINYEAR,proportion=proportion))
}

plt.miss.yr_blk=function(DATA)   #function for plotting missing year-blocks
{
  this=aggregate((Catch.Target/Km.Gillnet.Days.c)~FINYEAR+BLOCKX,DATA,mean)
  Reshaped=(reshape(this,v.names = "(Catch.Target/Km.Gillnet.Days.c)", idvar = "FINYEAR",
                    timevar = "BLOCKX", direction = "wide"))
  FinYr=sort(as.character(Reshaped[,1]))
  N.yrs=length(FinYr)
  Blk.lab=sort(unique(DATA$BLOCKX))
  Blk=1:length(Blk.lab)
  
  numInt=20
  couleurs=rev(gray(seq(0.2,0.9,length=numInt)))  
  CPUEBreaks=quantile(unlist(Reshaped[,2:ncol(Reshaped)]),probs=seq(0,1,1/numInt),na.rm=T)
  numberLab=2
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  Reshaped=Reshaped[order(Reshaped[,1]),]
  Reshaped=Reshaped[,-1]
  Reshaped=as.matrix(Reshaped)
  
  par(mfcol=c(1,1),mai=c(0.8,1.0,.1,1),oma=c(0.01,.1,.1,.01),las=1)
  image(1:length(FinYr),Blk,z=Reshaped,axes = FALSE,ylab="Block",xlab="Financial year",
        col =couleurs,breaks=CPUEBreaks,cex.lab=1.75)
  axis(1,at=1:N.yrs,labels=F,tck=-0.01)
  axis(1,at=seq(1,N.yrs,by=5),labels=FinYr[seq(1,N.yrs,by=5)],tck=-0.02,cex.axis=1.05)
  
  axis(2,at=Blk,labels=F,tck=-0.01)
  axis(2,at=seq(1,length(Blk),by=5),labels=Blk.lab[seq(1,length(Blk),by=5)],tck=-0.02,cex.axis=1.05)
  color.legend(N.yrs+2,length(Blk.lab),N.yrs+4,length(Blk.lab)*.60,round(CPUEBreaks,2),rect.col=couleurs,
               gradient="y",col=colLeg,cex=1,align="rb")
  box()
  Missing.yr.blks=sum(is.na(Reshaped))
  return(Missing.yr.blks)
}

Check.deg.free=function(DATA)  #function for checking observations per Month-Year-Block
{
  Dat.zn.W=subset(DATA,zone=="West")
  Dat.zn.Zn1=subset(DATA,zone=="Zone1")
  Dat.zn.Zn2=subset(DATA,zone=="Zone2")
  
  Yr.Mn.W=with(Dat.zn.W,table(FINYEAR,MONTH))
  Yr.Mn.Zn1=with(Dat.zn.Zn1,table(FINYEAR,MONTH))
  Yr.Mn.Zn2=with(Dat.zn.Zn2,table(FINYEAR,MONTH))
  Yr.Blk.W=with(Dat.zn.W,table(FINYEAR,BLOCKX))
  Yr.Blk.Zn1=with(Dat.zn.Zn1,table(FINYEAR,BLOCKX))
  Yr.Blk.Zn2=with(Dat.zn.Zn2,table(FINYEAR,BLOCKX))
  Yr.Blk.Mn.W=with(Dat.zn.W,table(FINYEAR,MONTH,BLOCKX))
  Yr.Blk.Mn.Zn1=with(Dat.zn.Zn1,table(FINYEAR,MONTH,BLOCKX))
  Yr.Blk.Mn.Zn2=with(Dat.zn.Zn2,table(FINYEAR,MONTH,BLOCKX))
  return(list(Yr.Mn.W=Yr.Mn.W,Yr.Mn.Zn1=Yr.Mn.Zn1,Yr.Mn.Zn2=Yr.Mn.Zn2,Yr.Blk.W=Yr.Blk.W,Yr.Blk.Zn1=Yr.Blk.Zn1,
              Yr.Blk.Zn2=Yr.Blk.Zn2,Yr.Blk.Mn.W=Yr.Blk.Mn.W,Yr.Blk.Mn.Zn1=Yr.Blk.Mn.Zn1,Yr.Blk.Mn.Zn2=Yr.Blk.Mn.Zn2))
}

fn.explore=function(DATA,NAMES)    #function for exploring data
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
  id=match(c(Predictors,"CPUE"),names(DATA))
  id=subset(id,!is.na(id))
  DATA$SP.Target=ifelse(DATA$CPUE>0,1,0)
  
  #Tables
  TABLE1=sort(with(DATA,table(VESSEL)))
  TABLE2=sort(with(DATA,table(BLOCKX)))
  TABLE3=with(DATA,table(FINYEAR,MONTH))
  TABLE3.3=with(DATA,table(BLOCKX,FINYEAR))
  TABLE4=with(DATA,table(MONTH,BLOCKX,FINYEAR))
  TABLE5=with(DATA,table(VESSEL,FINYEAR))
  
  TABLE6=aggregate(Catch.Target~VESSEL,data=DATA,mean)
  names(TABLE6)[2]="Mean Catch"
  TABLE6.1=aggregate(Catch.Target~VESSEL,data=DATA,sd)
  names(TABLE6.1)[2]="SD Catch"
  TABLE6=merge(TABLE6,data.frame(VESSEL=names(TABLE1),Count=as.numeric(TABLE1)),by="VESSEL")
  TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
  TABLE6=TABLE6[order(TABLE6$Count),]
  
  TABLE7=aggregate(Catch.Target~BLOCKX,data=DATA,mean)
  names(TABLE7)[2]="Mean Catch"
  TABLE7.1=aggregate(Catch.Target~BLOCKX,data=DATA,sd)
  names(TABLE7.1)[2]="SD Catch"
  TABLE7=merge(TABLE7,data.frame(BLOCKX=names(TABLE2),Count=as.numeric(TABLE2)),by="BLOCKX")
  TABLE7=merge(TABLE7,TABLE7.1,by="BLOCKX")
  TABLE7=TABLE7[order(TABLE7$Count),]
  
  N=with(DATA,table(FINYEAR))
  N=c(N)
  
  TABLE8=aggregate(Catch.Target~FINYEAR,data=DATA,sum)
  #TABLE8=aggregate(Catch.Target~FINYEAR,data=DATA,mean)
  TABLE9=aggregate(Catch.Target~FINYEAR,data=DATA,sd)
  TABLE10=aggregate(CPUE~FINYEAR,data=DATA,mean)
  TABLE11=aggregate(CPUE~FINYEAR,data=DATA,sd)
  
  TABLE.eff=aggregate(Km.Gillnet.Days.c~FINYEAR,data=DATA,sum)
  #  TABLE.eff=aggregate(Km.Gillnet.Days.c~FINYEAR,data=DATA,mean)
  TABLE.eff1=aggregate(Km.Gillnet.Days.c~FINYEAR,data=DATA,sd)
  
  LEVELS=with(DATA,c(length(levels(as.factor(BLOCKX))),length(levels(VESSEL))))
  names(LEVELS)=c("Block", "Vessel")
  
  #cumulative catch
  TABLE12=aggregate(Catch.Target~VESSEL,data=DATA,sum)
  TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
  TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
  TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
  Vessels.with.catch=nrow(subset(TABLE12,Catch.Target>0))
  #TopVessels=subset(TABLE12,PerCumCatch<=CatchThreshold, select=c(VESSEL))
  
  TABLE13=aggregate(Catch.Target~BLOCKX,data=DATA,sum)
  TABLE13=TABLE13[order(-TABLE13$Catch.Target),]
  TABLE13$CumCatch=cumsum(TABLE13$Catch.Target)
  TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$Catch.Target),2)
  Blocks.with.catch=nrow(subset(TABLE13,Catch.Target>0))
  #TopBlocks.catch=subset(TABLE13,PerCumCatch<=CatchThreshold, select=c(BLOCKX))
  TABLE14=aggregate(SP.Target~BLOCKX,data=DATA,sum)
  TABLE14=merge(TABLE14,TABLE13,by='BLOCKX')
  names(TABLE14)[2]="Number.obs"
  #TABLE14$Top=ifelse(TABLE14$BLOCKX%in%TopBlocks[,1],"grey55","black")
  
  TABLE15=aggregate(SP.Target~VESSEL,data=DATA,sum)
  TABLE15=merge(TABLE15,TABLE12,by='VESSEL')
  names(TABLE15)[2]="Number.obs"
  #TABLE15$Top=ifelse(TABLE15$VESSEL%in%TopVessels[,1],"grey45","black")
  
  
  #cumulative records
  TABLE16=rev(sort(table(DATA$VESSEL)))
  TABLE16=data.frame(Records=as.numeric(TABLE16),VESSEL=names(TABLE16))
  TABLE16$CumRecords=cumsum(TABLE16$Records)
  TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
  TopVessels.records=subset(TABLE16,PerCumRecords<=RecordThreshold, select=c(VESSEL))
  TABLE16.b=subset(TABLE16,VESSEL%in%TopVessels.records$VESSEL,select=c("Records","VESSEL"))
  agg=subset(TABLE16,!(VESSEL%in%TopVessels.records$VESSEL))
  Per.grouped.rec.vesl=round(100*(sum(agg$Records)/(sum(TABLE16.b$Records)+sum(agg$Records))),1)
  TABLE16.b=rbind(TABLE16.b,data.frame(Records=sum(agg$Records),VESSEL=as.factor("PlusGroup")))
  
  Blocks.zones=unique(data.frame(BLOCKX=DATA$BLOCKX,zone=DATA$zone))
  TABLE17=rev(sort(table(DATA$BLOCKX)))
  TABLE17=data.frame(Records=as.numeric(TABLE17),BLOCKX=names(TABLE17))
  TABLE17$CumRecords=cumsum(TABLE17$Records)
  TABLE17$PerCumRecords=round(TABLE17$CumRecords*100/sum(TABLE17$Records),2)
  TABLE17.1=TABLE17
  TABLE17=merge(TABLE17,Blocks.zones,by="BLOCKX")
  TopBlocks.records=subset(TABLE17,PerCumRecords<=RecordThreshold, select=c(BLOCKX))
  TABLE17.b=subset(TABLE17,BLOCKX%in%TopBlocks.records$BLOCKX,select=c("Records","BLOCKX","zone"))
  Blocks.rec.lev=nrow(TABLE17.b)
  Blocks.rec.lev=c(nrow(TABLE17.b[TABLE17.b$zone=="West",]),nrow(TABLE17.b[TABLE17.b$zone=="Zone1",]),
                   nrow(TABLE17.b[TABLE17.b$zone=="Zone2",]))
  agg=subset(TABLE17,!(BLOCKX%in%TopBlocks.records$BLOCKX))
  agg.grouped=aggregate(Records~zone,data=agg,sum)
  agg.grouped$BLOCKX=as.factor(paste("PlusGroup.",agg.grouped$zone,sep="")) 
  
  agg$BLOCKX.new=as.factor(paste("PlusGroup.",agg$zone,sep=""))
  
  Blocks.rec.lev.grouped=c(nrow(agg[agg$zone=="West",]),nrow(agg[agg$zone=="Zone1",]),nrow(agg[agg$zone=="Zone2",]))
  names(Blocks.rec.lev.grouped)=names(Blocks.rec.lev)=ZONAS
  Per.grouped.rec.blk=round(100*(agg.grouped$Records/(aggregate(Records~zone,TABLE17.b,sum)[2]+agg.grouped$Records)),1)
  #   rownames(Per.grouped.rec.blk)=agg.grouped$zone
  #   TABLE17.b=rbind(TABLE17.b,agg.grouped[,match(c("Records","BLOCKX","zone"),names(agg.grouped))])
  
  
  # N.vessels.kept=round(100*nrow(subset(DATA,VESSEL%in%TopVessels[,1]))/nrow(DATA),0)
  #  N.blocks.kept=round(100*nrow(subset(DATA,BLOCKX%in%TopBlocks[,1]))/nrow(DATA),0)
  
  # N.vessel.and.blocks.kept=round(100*nrow(subset(DATA,BLOCKX%in%TopBlocks[,1]&VESSEL%in%TopVessels[,1]))/nrow(DATA),0)
  
  FINYEAR.monthly=as.character(unique(DATA$FINYEAR))
  a=function()axis(1,at=seq(1,nrow(TABLE8),by=1),labels=F,tck=-0.01)
  b=function()axis(1,at=seq(1,nrow(TABLE8),5),labels=FINYEAR.monthly[seq(1,nrow(TABLE8),5)],tck=-0.04)
  
  plot(1:nrow(TABLE8),TABLE8[,2]/1000,ylab="Total catch (tonnes)",xlab="Year",cex.lab=1.25,type="o",xaxt='n',pch=19)
  #segments(1:nrow(TABLE8),TABLE8[,2]-TABLE9[,2]/(N^0.5),1:nrow(TABLE8),TABLE8[,2]+TABLE9[,2]/(N^0.5))
  a();b()
  plot(1:nrow(TABLE.eff),TABLE.eff[,2]/1000,ylab="Total Effort (1000 km net days)",xlab="Year",cex.lab=1.25,type="o",xaxt='n',pch=19)
  #segments(1:nrow(TABLE.eff),TABLE.eff[,2]-TABLE.eff1[,2]/(N^0.5),1:nrow(TABLE.eff),TABLE.eff[,2]+TABLE.eff1[,2]/(N^0.5))
  a();b()
  #   plot(1:nrow(TABLE10),TABLE10[,2],ylab="Mean cpue (±SE)",xlab="Year",cex.lab=1.25,type="o",xaxt='n',pch=19)
  #   segments(1:nrow(TABLE10),TABLE10[,2]-TABLE11[,2]/(N^0.5),1:nrow(TABLE10),TABLE10[,2]+TABLE11[,2]/(N^0.5))
  #   a();b()
  pos.cpue=DATA$CPUE[DATA$CPUE>0]
  hist(log(pos.cpue),main=paste(nrow(DATA)-length(pos.cpue),"records with 0 catch (",
                                round(100*(nrow(DATA)-length(pos.cpue))/nrow(DATA),0),"%)"))
  box()
  
  
  #indicative vessels
  #1st, in the fishery for >= Threshold.n.yrs
  Vessel.n.yrs=with(DATA,table(VESSEL,FINYEAR))
  Vessel.n.yrs=ifelse(Vessel.n.yrs>=1,1,0)
  Vessel.n.yrs=rowSums(Vessel.n.yrs)
  Indic.ves=names(Vessel.n.yrs[which(Vessel.n.yrs>=Threshold.n.yrs)])
  
  DATA.indic=subset(DATA,VESSEL%in%Indic.ves)
  DATA.less.5yrs=subset(DATA,!(VESSEL%in%Indic.ves))
  
  #2nd, median total catch >= fishery median catch
  Threshold.median.catch=median(DATA$Catch.Target)
  Vessel.median.catch=aggregate(Catch.Target~VESSEL,data=DATA.indic,median)
  Indic.ves=as.character(Vessel.median.catch[which(Vessel.median.catch$Catch.Target>=Threshold.median.catch),1])
  Kept.indic=round(100*nrow(subset(DATA,VESSEL%in%Indic.ves))/nrow(DATA))
  
  #   hist(Vessel.n.yrs,breaks=35,main="",xlab="Number of years in fishery",cex.lab=1.25)
  #   box()
  #   
  #   hist(Vessel.median.catch$Catch.Target,breaks=50,main="",xlab="Median catch (kg)",cex.lab=1.25)
  #   box()
  #   legend("topright",paste("Percent records kept (indicative ves)=",Kept.indic),bty='n',cex=1.25)
  
  
  #Blocks by year-month
  #   DATA$Yr.Mn=with(DATA,paste(YEAR.c,MONTH,sep="."))
  #   Block.Yr.Mn=with(DATA,table(Yr.Mn,BLOCKX))
  #   Block.Yr.Mn=ifelse(Block.Yr.Mn>=1,1,0)
  #   Block.Yr.Mn=colSums(Block.Yr.Mn)
  
  return(list(Factor.levels=LEVELS,Summary=summary(DATA[,id]),Vessels=TABLE1,Blocks=TABLE2,
              Yr_Mn=TABLE3,Yr_Block=TABLE3.3,Mn_Blokc_Yr=TABLE4,Vessel_Yr=TABLE5,Mean_catch_Vessel=TABLE6,
              Mean_catch_Block=TABLE7,Indicative.ves=Indic.ves,
              Ves.Cum.Ca=TABLE12$PerCumCatch,Block.Cum.Ca=TABLE13$PerCumCatch,
              TABLE14=TABLE14,TABLE12=TABLE12,
              TABLE15=TABLE15,TABLE16=TABLE16,TopVessels.records=TopVessels.records,
              TABLE16.b=TABLE16.b,Per.grouped.rec.vesl=Per.grouped.rec.vesl,
              Vessels.with.catch=Vessels.with.catch,Blocks.with.catch=Blocks.with.catch,
              TABLE17=TABLE17.1,TopBlocks.records=TopBlocks.records,TABLE17.b=TABLE17.b,
              Per.grouped.rec.blk=Per.grouped.rec.blk,Blocks.rec.lev=Blocks.rec.lev,
              Blocks.rec.lev.grouped=Blocks.rec.lev.grouped,agg=agg))
}

fn.catch.cpue=function(DATA)   #function for plotting catch VS cpue
{
  DATA$CPUE=with(DATA,Catch.Target/Km.Gillnet.Days.c)
  KTCH=aggregate(Catch.Target~FINYEAR,DATA,sum,na.rm=T)
  CPUE=aggregate(CPUE~FINYEAR,DATA,mean,na.rm=T)
  max1=max(KTCH[,2]/1000,na.rm=T)*1.05
  max2=max(CPUE[,2],na.rm=T)*1.05
  FInYEAR=as.character(unique(DATA$FINYEAR))
  N=length(FInYEAR)
  plot(1:N,(KTCH[,2]/1000),ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
       ,cex.axis=1.2,lwd=1.75)
  axis(1,at=1:N,labels=F,tck=-0.015)
  axis(1,at=seq(1,N,5),labels=FInYEAR[seq(1,N,5)],tck=-0.0225,cex.axis=1.1)
  par(new=T)
  plot(1:N,CPUE[,2],col="grey70",type='l',axes=F,ann='F',lwd=2,ylim=c(0,max2))
  axis(4,at=pretty(CPUE[,2]),labels=pretty(CPUE[,2]),las=2,cex.axis=1.2)
  mtext("Year",1,1.25,outer=T,cex=1.35)
  mtext("Total catch (tonnes)",2,-1.15,outer=T,las=3,cex=1.35)
  mtext("CPUE (kg/km gn.d)",4,-1.5,outer=T,las=3,cex=1.35)
  legend("bottom",c("catch","cpue"),col=c("black","grey70"),lty=1,bty='n',lwd=2,cex=1.5)
  mtext("Increase in catch and decrease in CPUE = increase in Q or efficiency",3)
}

Wide.vs.long=function(D.wide,D.long)  #function for comparing long VS wide data
{
  if(TARGETS[[i]][1]%in%18007)
  {
    S=paste(1988:2005,substr(1989:2006,3,4),sep="-")
    D.wide=subset(D.wide,FINYEAR%in%S)
    D.long=subset(D.long,FINYEAR%in%S)
  }
  
  par(mfcol=c(3,1),mai=c(.5,1,.1,.1),mgp=c(3,0.6,0),las=1)
  long=aggregate(LIVEWT.c~FINYEAR,D.long,sum)
  wide=aggregate(Catch.Target~FINYEAR,D.wide,sum)  
  plot(long[,2],pch=19,ylab="catch")
  lines(wide[,2],col=2)
  
  long.e=aggregate(Km.Gillnet.Days.c~zone+
                     FINYEAR+Same.return+MONTH+BLOCKX,D.long,max,na.rm=T)
  long.e=aggregate(Km.Gillnet.Days.c~FINYEAR,D.long,sum,na.rm=T)
  wide.e=aggregate(Km.Gillnet.Days.c~FINYEAR,D.wide,sum,na.rm=T)  
  plot(long.e[,2],pch=19,ylab="effort")
  lines(wide.e[,2],col=2)
  
  cpue.long=long[,2]/long.e[,2]
  cpue.wide=wide[,2]/wide.e[,2]
  plot(cpue.long,pch=19,ylab="cpue")
  lines(cpue.wide,col=2)
  legend("topleft",c("wide"),bty='n',lty=1,col=2,cex=2)
  
  
}
Wide.vs.long.daily=function(D.wide,D.long,DATA)
{
  
  par(mfcol=c(3,1),mai=c(.5,1,.1,.1),mgp=c(3,0.6,0),las=1)
  
  long=aggregate(LIVEWT.c~FINYEAR,D.long,sum)
  wide=aggregate(Catch.Target~FINYEAR,D.wide,sum)  
  plot(long[,2],pch=19,ylab="catch")
  lines(wide[,2],col=2)
  
  long.e=aggregate(cbind(Km.Gillnet.Days.c,SOI,Freo,Freo.Lag6,Freo.Lag12)~zone+
                     FINYEAR+date+TSNo+MONTH+BLOCKX+block10+VESSEL,DATA,max)
  long.e=aggregate(Km.Gillnet.Days.c~zone+FINYEAR+TSNo+
                     MONTH+BLOCKX+block10+VESSEL,long.e,sum)
  long.e=aggregate(Km.Gillnet.Days.c~FINYEAR,long.e,sum,na.rm=T)
  wide.e=aggregate(Km.Gillnet.Days.c~FINYEAR,D.wide,sum,na.rm=T)  
  plot(long.e[,2],pch=19,ylab="effort")
  lines(wide.e[,2],col=2)
  
  cpue.long=long[,2]/long.e[,2]
  cpue.wide=wide[,2]/wide.e[,2]
  plot(cpue.long,pch=19,ylab="cpue")
  lines(cpue.wide,col=2)
  legend("topleft",c("wide"),bty='n',lty=1,col=2,cex=2) 
}

fun.ves.blk.zone=function(DAT)  #function for Number of vessels and blocks per zone and species
{
  V=table(DAT$VESSEL,DAT$zone)
  V=ifelse(V>0,1,V)
  n.ves=colSums(V)
  
  B=table(DAT$BLOCKX,DAT$zone)
  B=ifelse(B>0,1,B)
  n.blok=colSums(B)
  
  return(list(n.ves=n.ves,n.blok=n.blok))
}

fn.new.ves=function(DATA,Ves.No.ktch,Per.int)   #function for assigning vessels to categories
{
  #drop vessels that never reported catch of target species
  DATA=subset(DATA,!(VESSEL%in%Ves.No.ktch))
  DATA=DATA[,match(VARIABLES,names(DATA))]
  
  # Create positive (non-zero) dataset
  PosData <- subset(DATA,Catch.Target>0)
  PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c
  PosData$log.Catch=log(PosData$Catch.Target+0.00001)
  PosData$log.Effort=log(PosData$Km.Gillnet.Days.c+0.00001)
  
  #fit simple mixed model for ranking all vessels each year
  YRs=sort(unique(PosData$FINYEAR))
  STORE=vector('list',length(YRs))
  for (y in 1:length(YRs))
  {
    dat=subset(PosData,FINYEAR==YRs[y])
    dat$VESSEL=as.factor(dat$VESSEL)
    dat$BLOCKX=as.factor(dat$BLOCKX)
    simple.mixed=lme(log.Catch~BLOCKX+offset(log.Effort), data = dat, random = ~ 1 |VESSEL)
    
    #Predict catches for each vessel
    Eff.unit=mean(dat$log.Effort)
    a=names(sort(table(dat$BLOCKX)))
    a=a[length(a)]
    Pred.dat=expand.grid(log.Effort=Eff.unit,BLOCKX=factor(a,levels=levels(dat$BLOCKX)),
                         VESSEL=levels(dat$VESSEL))
    Catch.pred=predict(simple.mixed,newdata=Pred.dat,type='response')
    Catch.pred=data.frame(Pred.ktch=Catch.pred,VESSEL=names(Catch.pred))
    Catch.pred=Catch.pred[order(-Catch.pred$Pred.ktch),]
    Catch.pred$Ranking.Ktch.pred=1:nrow(Catch.pred)
    Catch.pred$FINYEAR=YRs[y]
    
    
    # extract random coefficients
    Rand.eff=ranef(simple.mixed)
    vessel.coef=data.frame(VESSEL=rownames(Rand.eff),coef=Rand.eff)
    
    #get percentiles and sort
    percentiles=unique(quantile(vessel.coef$X.Intercept.,probs=seq(0,1,Per.int)))
    names(percentiles)=0:(length(percentiles)-1)
    a=aggregate(cpue~VESSEL,dat,mean)
    b=aggregate(cpue~VESSEL,dat,sd)
    a=merge(a,b,by="VESSEL")
    names(a)[2:3]=c("Mean","SD")
    a=merge(a,vessel.coef,by="VESSEL",all.x=T)
    
    #Create new vessel levels
    a$new.level.vessel=cut(a$X.Intercept.,percentiles,include.lowest = T, 
                           right = T,labels=names(percentiles)[2:length(percentiles)])
    
    
    #Re order factor levels
    a$new.level.vessel <- factor(a$new.level.vessel,levels = names(percentiles))
    
    #Add predicted catch
    a=merge(a,Catch.pred,by="VESSEL")
    a=a[,-match(c("Mean","SD"),names(a))]
    
    STORE[[y]]=a
  }
  return(vessel.coef=do.call(rbind,STORE))
}

fn.plot.nw.ves.cat=function(dat,NM)   #function for plotting new vessel categories
{
  ktch=aggregate(Pred.ktch~new.level.vessel,dat,mean)
  ktch.sd=aggregate(Pred.ktch~new.level.vessel,dat,sd)
  n=1:length(ktch$new.level.vessel)
  plot(n,ktch$Pred.ktch,ylab="",xlab="",xaxt='n',pch=19,cex=1.25,ylim=c(0,max(ktch$Pred.ktch+ktch.sd$Pred.ktch)))
  segments(n,ktch$Pred.ktch-ktch.sd$Pred.ktch,n,ktch$Pred.ktch+ktch.sd$Pred.ktch,lwd=2)
  axis(1,n)
  legend("bottomright",NM,bty='n',cex=1.5)
}

fn.add.new.ves=function(DAT,Ves.no.kt,DAT1) #function for adding new vessel categories to data
{
  DAT=subset(DAT,!(VESSEL%in%Ves.no.kt))
  DAT1=DAT1[,match(c("VESSEL","new.level.vessel","FINYEAR"),names(DAT1))]
  DAT=merge(DAT,DAT1,by=c("VESSEL","FINYEAR"),all.x=T)
  return(DAT)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)   #function for correlations
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}



#function for calculating folly and nominal
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
export.nominal=function(DATA,DATA1)     #function for creating nominal indices 
{
  DATA=subset(DATA,Catch.Target>0)
  DATA1=subset(DATA1,Catch.Target>0)
  
  DATA$cpue=DATA$Catch.Target/DATA$Km.Gillnet.Days.c
  DATA1$cpue=DATA1$Catch.Target/DATA1$Km.Gillnet.Days.c
  Nominal=aggregate(cpue~FINYEAR,DATA,mean)
  Nominal_d=aggregate(cpue~FINYEAR,DATA1,mean)
  return(rbind(Nominal,Nominal_d))
}


#functions for defining model structure
fn.define.models=function(res.var,Interaction)  #function for defining model structures
{
  Species.mod=Species.model
  for(s in 1:N.species)
  {
    #binomial
    PREDS=Predictors[-match("Km.Gillnet.Days.c",Predictors)]    
    PREDS=PREDS[-drop.pred[s]]
    log.this=c("Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar","Catch.Scalefish",
               "Freo","Freo.Lag6","Freo.Lag12")
    id=match(log.this,PREDS)
    id=id[!is.na(id)]
    PREDS[id]=paste("log.",PREDS[id],sep="")
    bin.form=vector('list',length=N.models)
    bin.form[[1]]=paste(PREDS[match(c("FINYEAR","BLOCKX","VESSEL"),PREDS)],collapse="+")
    PREDS=PREDS[-match(c("FINYEAR","BLOCKX","VESSEL"),PREDS)]
    for(x in 2:length(bin.form)) bin.form[[x]]=paste(bin.form[[x-1]],"+",PREDS[x-1],sep="")
    for(x in 1:length(bin.form)) 
    {
      bin.form[[x]]=as.formula(paste(paste("Catch.Target", "~", collapse=NULL),
                                     bin.form[[x]],"+offset(log.Effort)")) 
    }
    
    
    #poscatch
    PREDS=Predictors[-match("Km.Gillnet.Days.c",Predictors)]
    PREDS=PREDS[-drop.pred[s]]
    id=match(log.this,PREDS)
    id=id[!is.na(id)]
    PREDS[id]=paste("log.",PREDS[id],sep="")    
    pos.log.form=vector('list',length=N.models)
    pos.log.form[[1]]=paste("FINYEAR",Interaction,"BLOCKX","+","VESSEL",sep="")
    PREDS=PREDS[-match(c("FINYEAR","BLOCKX","VESSEL"),PREDS)]
    for(x in 2:length(pos.log.form)) pos.log.form[[x]]=paste(pos.log.form[[x-1]],"+",PREDS[x-1],sep="")
    pos.form=pos.log.form
    
    if(res.var=="cpue")
    {
      for(x in 1:length(pos.log.form)) 
      {
        pos.log.form[[x]]=as.formula(paste(paste("log.cpue", "~", collapse=NULL), pos.log.form[[x]])) 
        pos.form[[x]]=as.formula(paste(paste("cpue", "~", collapse=NULL), pos.form[[x]])) 
      }
    }
    
    if(res.var=="catch")
    {
      for(x in 1:length(bin.form)) 
      {
        pos.log.form[[x]]=as.formula(paste(paste("log.Catch", "~", collapse=NULL),
                                           pos.log.form[[x]],"+offset(log.Effort)")) 
        pos.form[[x]]=as.formula(paste(paste("Catch.Target", "~", collapse=NULL),pos.form[[x]],"+offset(log.Effort)"))
      }
    }
    Species.mod[[s]]=list(bin.form=bin.form,pos.log.form=pos.log.form,pos.form=pos.form)
    
  } 
  return(Species.mod)
}
Test.models=function(DATA,No.ves,SPEC,MIN.Wght)   #function for testing models
{
  #1. Create storing objects
  form.bin=vector("list",length=length(bin.form))
  form.pos.log=vector("list",length=length(pos.log.form))
  form.pos=vector("list",length=length(pos.form))
  
  names(form.bin)=bin.form
  names(form.pos.log)=pos.log.form
  names(form.pos)=pos.form
  
  #2. Drop vessels that never reported catch of target species   
  DATA=subset(DATA,!(VESSEL%in%No.ves))
  DATA=subset(DATA,Catch.Target>=MIN.Wght| Catch.Target==0)
  
  #3. Put data in proper shape
  DataFile=DATA[,match(VARIABLES,names(DATA))]
  
  # convert to factor and sort levels
  DataFile$FINYEAR=as.factor(DataFile$FINYEAR)
  DataFile$BLOCKX=as.factor(DataFile$BLOCKX)
  DataFile$MONTH=as.factor(DataFile$MONTH)
  
  DataFile$FINYEAR=factor(DataFile$FINYEAR,levels=(names(table(DataFile$FINYEAR))))
  DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=(names(table(DataFile$BLOCKX))))
  DataFile$MONTH=factor(DataFile$MONTH,levels=(names(table(DataFile$MONTH))))
  
  # log effort
  DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c)
  
  # log catch other
  DataFile$log.Catch.Gummy=log(DataFile$Catch.Gummy+0.000000001)
  DataFile$log.Catch.Whiskery=log(DataFile$Catch.Whiskery+0.000000001)
  DataFile$log.Catch.Dusky=log(DataFile$Catch.Dusky+0.000000001)
  DataFile$log.Catch.Sandbar=log(DataFile$Catch.Sandbar+0.000000001)
  DataFile$log.Catch.Scalefish=log(DataFile$Catch.Scalefish+0.000000001)
  DataFile$log.Freo=log(DataFile$Freo)
  DataFile$log.Freo.Lag6=log(DataFile$Freo.Lag6)
  DataFile$log.Freo.Lag12=log(DataFile$Freo.Lag12)
  
  
  # create binary dataset (0/1 for CPUE)
  BiData <- DataFile
  BiData[,1] <- as.numeric(DataFile[,1]>0)
  
  # create positive (non-zero) dataset
  PosData <- DataFile[DataFile[,1]>0,]  
  PosData$log.Catch=log(PosData$Catch.Target)
  PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c
  PosData$log.cpue=log(PosData$cpue)
  
  
  #4. Loop over all possible model structures    
  #note: too big GLM so have to export
  
  # binomial
  NAME.bi=rep(NA,length(bin.form))
  for(i in 1:length(bin.form))
  {
    Bin.mod=NULL
    if(!is.null(bin.form[[i]])) Bin.mod=glm(bin.form[[i]], data=BiData, family="binomial", maxit=500)  
    NAME.bi[i]=paste(SPEC,".Bi",i,".RData",sep="")
    save(Bin.mod,file=NAME.bi[i])
    rm(Bin.mod)
  }
  
  #Lognormal
  NAME.pos=rep(NA,length(pos.log.form))
  for(i in 1:length(pos.log.form))
  {
    LogN.model=glm(pos.log.form[[i]], data=PosData, family=gaussian, maxit=500)
    NAME.pos[i]=paste(SPEC,".LogN",i,".RData",sep="")
    save(LogN.model,file=NAME.pos[i])
    rm(LogN.model)
  }
  
  #Gamma
  NAME.pos.gamma=rep(NA,length(pos.form))
  for(i in 1:length(pos.form))
  {
    Gamma.model=glm(pos.form[[i]], data=PosData, family=Gamma(link=log), maxit=500)
    NAME.pos.gamma[i]=paste(SPEC,".Gamma",i,".RData",sep="")
    save(Gamma.model,file=NAME.pos.gamma[i])
    rm(Gamma.model)
  }
  
  return(list(NAME.bi=NAME.bi,NAME.pos=NAME.pos,NAME.pos.gamma=NAME.pos.gamma))
}
load.and.select=function(dir,N.spec)   
{
  setwd(dir)
  for(j in N.spec) All.models[[j]]=Names.models(SPECIES.vec[j])
  
  BESTS=All.models
  for(j in N.spec)
  {
    NMs=All.models[[j]]
    Bi=NMs$NAME.bi
    Pos=NMs$NAME.pos
    Gam=NMs$NAME.pos.Gamma
    
    #Binomial models
    AIC.bin=Dev.bin=NULL
    for(i in 1:length(Bi))
    {
      load(Bi[[i]])
      if(!is.null(Bin.mod))
      {
        AIC.bin[i]=fn.AICc(Bin.mod)
        Dev.bin[i]=fn.Dev.Exp(Bin.mod)
        rm(Bin.mod)
      }
      
    }
    
    #Positive catch models
    #Lognormal
    for(i in 1:length(Pos))
    {
      load(Pos[[i]])
      AIC.Pos[i]=fn.AICc(LogN.model)
      Dev.pos[i]=fn.Dev.Exp(LogN.model)
      rm(LogN.model)
    }
    
    #Gamma
    for(i in 1:length(Gam))
    {
      load(Gam[[i]])
      AIC.Pos.Gamma[i]=fn.AICc(Gamma.model)
      Dev.pos.Gamma[i]=fn.Dev.Exp(Gamma.model)
      rm(Gamma.model)
    }
    
    
    
    #Store things  
    #bin
    Aic.Analysis.Bin=NULL
    if(!is.null(AIC.bin))
    {
      Aic.Analysis.Bin=fn.AIC.BIC.ratio(AIC.bin)
      Aic.Analysis.Bin$Best.Mod=bin.form[Aic.Analysis.Bin$Best.Mod]
    }
    #Pos
    #Lognormal
    Aic.Analysis.Pos=fn.AIC.BIC.ratio(AIC.Pos)
    Aic.Analysis.Pos$Best.Mod=pos.log.form[Aic.Analysis.Pos$Best.Mod]
    
    #Gamma
    Aic.Analysis.Pos.Gamma=fn.AIC.BIC.ratio(AIC.Pos.Gamma)
    Aic.Analysis.Pos.Gamma$Best.Mod=pos.form[Aic.Analysis.Pos.Gamma$Best.Mod]
    
    #Put together
    Dev.explained=list(bi=Dev.bin,pos=Dev.pos,pos.Gamma=Dev.pos.Gamma)
    BESTS[[j]]=list(AIC.bin=AIC.bin,AIC.Pos=AIC.Pos,AIC.Pos.Gamma=AIC.Pos.Gamma,
                    Dev.explained=Dev.explained,
                    Aic.Analysis.Bin=Aic.Analysis.Bin,Aic.Analysis.Pos=Aic.Analysis.Pos,
                    Aic.Analysis.Pos.Gamma=Aic.Analysis.Pos.Gamma)
    
  }
  return(BESTS)
}
Names.models=function(SPEC)  
{
  #1. Create storing objects
  form.bin=vector("list",length=length(bin.form))
  form.pos.log=vector("list",length=length(pos.log.form))
  form.pos=vector("list",length=length(pos.form))
  names(form.bin)=bin.form
  names(form.pos.log)=pos.log.form
  names(form.pos)=pos.form
  
  # binomial
  NAME.bi=rep(NA,length(bin.form))
  for(i in 1:length(bin.form))  NAME.bi[i]=paste(SPEC,".Bi",i,".RData",sep="")
  #Lognormal
  NAME.pos=rep(NA,length(pos.log.form))
  for(i in 1:length(pos.log.form))NAME.pos[i]=paste(SPEC,".LogN",i,".RData",sep="")  
  #Gamma
  NAME.pos.Gamma=rep(NA,length(pos.form))
  for(i in 1:length(pos.form))NAME.pos.Gamma[i]=paste(SPEC,".Gamma",i,".RData",sep="")
  
  
  return(list(NAME.bi=NAME.bi,NAME.pos=NAME.pos,NAME.pos.Gamma=NAME.pos.Gamma))
} 

fn.see.bl10.closure=function(D)
{
  D=subset(D, LAT>(-33) & LAT<=(-31) & LONG<116)
  blck10=with(D,table(FINYEAR,block10))
  blck=with(D,table(FINYEAR,BLOCKX))
  return(list(blk10=blck10,BLCKX=blck))
}

#Comparing cpue, catch, random effects, etc functions
Compare.all=function(RESHAPED,RESHAPED.ALL.VES)   #function for testing models
{
  
  RESHAPED=subset(RESHAPED,Catch.Target>0)
  RESHAPED$log.Effort=log(RESHAPED$Km.Gillnet.Days.c)
  RESHAPED$log.Catch=log(RESHAPED$Catch.Target)
  RESHAPED$cpue=RESHAPED$Catch.Target/RESHAPED$Km.Gillnet.Days.c
  RESHAPED$log.cpue=log(RESHAPED$cpue)
  RESHAPED$FINYEAR=as.factor(RESHAPED$FINYEAR)
  RESHAPED$BLOCKX=as.factor(RESHAPED$BLOCKX)
  RESHAPED$MONTH=as.factor(RESHAPED$MONTH)
  RESHAPED$VESSEL=as.factor(RESHAPED$VESSEL)
  
  RESHAPED.ALL.VES=subset(RESHAPED.ALL.VES,Catch.Target>0)
  RESHAPED.ALL.VES$log.Effort=log(RESHAPED.ALL.VES$Km.Gillnet.Days.c)
  RESHAPED.ALL.VES$log.Catch=log(RESHAPED.ALL.VES$Catch.Target)
  RESHAPED.ALL.VES$cpue=RESHAPED.ALL.VES$Catch.Target/RESHAPED.ALL.VES$Km.Gillnet.Days.c
  RESHAPED.ALL.VES$log.cpue=log(RESHAPED.ALL.VES$cpue)
  RESHAPED.ALL.VES$FINYEAR=as.factor(RESHAPED.ALL.VES$FINYEAR)
  RESHAPED.ALL.VES$BLOCKX=as.factor(RESHAPED.ALL.VES$BLOCKX)
  RESHAPED.ALL.VES$MONTH=as.factor(RESHAPED.ALL.VES$MONTH)
  RESHAPED.ALL.VES$VESSEL=as.factor(RESHAPED.ALL.VES$VESSEL)
  
  #Nominal
  Nominal_RESHAPED=aggregate(log.cpue~FINYEAR,RESHAPED,mean)
  Nominal_RESHAPED.ALL.VES=aggregate(log.cpue~FINYEAR,RESHAPED.ALL.VES,mean)
  
  #Lognormal GLM
  
  #RESHAPED
  MODEL_RESHAPED.logcpue <- glm(log.cpue ~ FINYEAR , data=RESHAPED, family=gaussian, maxit=500)
  MODEL_RESHAPED.offset <- glm(log.Catch ~ FINYEAR + offset(log.Effort), data=RESHAPED, family=gaussian, maxit=500)
  MODEL_RESHAPED.main <- glm(log.Catch ~ FINYEAR + log.Effort, data=RESHAPED, family=gaussian, maxit=500)
  
  MODEL_RESHAPED.logcpue_vessel <- glm(log.cpue ~ FINYEAR+VESSEL, data=RESHAPED, family=gaussian, maxit=500)
  MODEL_RESHAPED.offset_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + offset(log.Effort), data=RESHAPED, family=gaussian, maxit=500)
  MODEL_RESHAPED.main_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + log.Effort, data=RESHAPED, family=gaussian, maxit=500)
  
  #RESHAPED.ALL.VES
  MODEL_RESHAPED.ALL.VES.logcpue <- glm(log.cpue ~ FINYEAR , data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  MODEL_RESHAPED.ALL.VES.offset <- glm(log.Catch ~ FINYEAR + offset(log.Effort), data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  MODEL_RESHAPED.ALL.VES.main <- glm(log.Catch ~ FINYEAR + log.Effort, data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  
  MODEL_RESHAPED.ALL.VES.logcpue_vessel <- glm(log.cpue ~ FINYEAR+VESSEL, data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  MODEL_RESHAPED.ALL.VES.offset_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + offset(log.Effort), data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  MODEL_RESHAPED.ALL.VES.main_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + log.Effort, data=RESHAPED.ALL.VES, family=gaussian, maxit=500)
  
  
  #Plot annual cpue
  
  #RESHAPED
  Preds_RESHAPED.logcpue=lsmeans(MODEL_RESHAPED.logcpue,specs=c("FINYEAR"))
  Preds_RESHAPED.offset=lsmeans(MODEL_RESHAPED.offset,specs=c("FINYEAR"))
  Preds_RESHAPED.main=lsmeans(MODEL_RESHAPED.main,specs=c("FINYEAR"))
  
  Preds_RESHAPED.logcpue_vessel=lsmeans(MODEL_RESHAPED.logcpue_vessel,specs=c("FINYEAR"))
  Preds_RESHAPED.offset_vessel=lsmeans(MODEL_RESHAPED.offset_vessel,specs=c("FINYEAR"))
  Preds_RESHAPED.main_vessel=lsmeans(MODEL_RESHAPED.main_vessel,specs=c("FINYEAR"))
  
  #RESHAPED.ALL.VES
  Preds_RESHAPED.ALL.VES.logcpue=lsmeans(MODEL_RESHAPED.ALL.VES.logcpue,specs=c("FINYEAR"))
  Preds_RESHAPED.ALL.VES.offset=lsmeans(MODEL_RESHAPED.ALL.VES.offset,specs=c("FINYEAR"))
  Preds_RESHAPED.ALL.VES.main=lsmeans(MODEL_RESHAPED.ALL.VES.main,specs=c("FINYEAR"))
  
  Preds_RESHAPED.ALL.VES.logcpue_vessel=lsmeans(MODEL_RESHAPED.ALL.VES.logcpue_vessel,specs=c("FINYEAR"))
  Preds_RESHAPED.ALL.VES.offset_vessel=lsmeans(MODEL_RESHAPED.ALL.VES.offset_vessel,specs=c("FINYEAR"))
  Preds_RESHAPED.ALL.VES.main_vessel=lsmeans(MODEL_RESHAPED.ALL.VES.main_vessel,specs=c("FINYEAR"))
  
  
  #RESHAPED
  logcpue=summary(Preds_RESHAPED.logcpue)
  offset=summary(Preds_RESHAPED.offset)
  main=summary(Preds_RESHAPED.main)  
  logcpue_RESHAPED=as.data.frame(logcpue[c("lsmean")])
  offset_RESHAPED=as.data.frame(offset[c("lsmean")])
  main_RESHAPED=as.data.frame(main[c("lsmean")])
  
  logcpue=summary(Preds_RESHAPED.logcpue_vessel)
  offset=summary(Preds_RESHAPED.offset_vessel)
  main=summary(Preds_RESHAPED.main_vessel)  
  logcpue_RESHAPED_vessel=as.data.frame(logcpue[c("lsmean")])
  offset_RESHAPED_vessel=as.data.frame(offset[c("lsmean")])
  main_RESHAPED_vessel=as.data.frame(main[c("lsmean")])
  
  
  #RESHAPED.ALL.VES
  logcpue=summary(Preds_RESHAPED.ALL.VES.logcpue)
  offset=summary(Preds_RESHAPED.ALL.VES.offset)
  main=summary(Preds_RESHAPED.ALL.VES.main)  
  logcpue_RESHAPED.ALL.VES=as.data.frame(logcpue[c("lsmean")])
  offset_RESHAPED.ALL.VES=as.data.frame(offset[c("lsmean")])
  main_RESHAPED.ALL.VES=as.data.frame(main[c("lsmean")])
  
  logcpue=summary(Preds_RESHAPED.ALL.VES.logcpue_vessel)
  offset=summary(Preds_RESHAPED.ALL.VES.offset_vessel)
  main=summary(Preds_RESHAPED.ALL.VES.main_vessel)  
  logcpue_RESHAPED.ALL.VES_vessel=as.data.frame(logcpue[c("lsmean")])
  offset_RESHAPED.ALL.VES_vessel=as.data.frame(offset[c("lsmean")])
  main_RESHAPED.ALL.VES_vessel=as.data.frame(main[c("lsmean")])
  
  
  #_RESHAPED
  fn2(MODEL_RESHAPED.logcpue,MODEL_RESHAPED.offset,MODEL_RESHAPED.main,
      MODEL_RESHAPED.logcpue_vessel,MODEL_RESHAPED.offset_vessel,MODEL_RESHAPED.main_vessel,
      logcpue_RESHAPED,offset_RESHAPED,main_RESHAPED,logcpue_RESHAPED_vessel,offset_RESHAPED_vessel,main_RESHAPED_vessel)
  mtext("Reshaped data, indicative vessels",side=3,line=-2)
  
  #_RESHAPED.ALL.VES
  fn2(MODEL_RESHAPED.ALL.VES.logcpue,MODEL_RESHAPED.ALL.VES.offset,MODEL_RESHAPED.ALL.VES.main,
      MODEL_RESHAPED.ALL.VES.logcpue_vessel,MODEL_RESHAPED.ALL.VES.offset_vessel,MODEL_RESHAPED.ALL.VES.main_vessel,
      logcpue_RESHAPED.ALL.VES,offset_RESHAPED.ALL.VES,main_RESHAPED.ALL.VES,
      logcpue_RESHAPED.ALL.VES_vessel,offset_RESHAPED.ALL.VES_vessel,main_RESHAPED.ALL.VES_vessel)
  mtext("Reshaped data, all vessels",side=3,line=-2)
}

fn2=function(MOD.logcpue,MOD.offset,MOD.main,MOD.logcpue_vessel,MOD.offset_vessel,MOD.main_vessel,
             logcpue,offset,main,logcpue_vessel,offset_vessel,main_vessel)
{
  dev.exp=paste(" (",round(100*c(Dsquared(MOD.logcpue),Dsquared(MOD.offset),Dsquared(MOD.main),
                                 Dsquared(MOD.logcpue_vessel),Dsquared(MOD.offset_vessel),
                                 Dsquared(MOD.main_vessel)),0),"%)",sep="")  
  
  logcpue$lsmean=exp(logcpue$lsmean)
  offset$lsmean=exp(offset$lsmean)
  main$lsmean=exp(main$lsmean)
  logcpue_vessel$lsmean=exp(logcpue_vessel$lsmean)
  offset_vessel$lsmean=exp(offset_vessel$lsmean)
  main_vessel$lsmean=exp(main_vessel$lsmean)
  
  A=c(logcpue$lsmean/mean(logcpue$lsmean),offset$lsmean/mean(offset$lsmean),
      main$lsmean/mean(main$lsmean),logcpue_vessel$lsmean/mean(logcpue_vessel$lsmean),
      offset_vessel$lsmean/mean(offset_vessel$lsmean),main_vessel$lsmean/mean(main_vessel$lsmean))
  Ylim=c(min(A),max(A))
  with(logcpue,plot(lsmean/mean(lsmean),ylim=Ylim,type='l',col=1,ylab="relative catch rate",cex.lab=1.35))
  with(offset,lines(lsmean/mean(lsmean),col=2))
  with(main,lines(lsmean/mean(lsmean),col=3))
  
  with(logcpue_vessel,lines(lsmean/mean(lsmean),col=1,lty=3))
  with(offset_vessel,lines(lsmean/mean(lsmean),col=2,lty=3))
  with(main_vessel,lines(lsmean/mean(lsmean),col=3,lty=3))
  
  legend("bottomleft",
         paste(c("cpue","offset","main term","cpue_vess","offset_vess","main term_vess"),dev.exp),
         lty=c(rep(1,3),rep(3,3)),col=rep(1:3,2),bty='n',cex=1.5)
  
}

fn3=function(offset,offset.rand)
{
  offset$lsmean=exp(offset$lsmean)
  offset.rand$lsmean=exp(offset.rand$lsmean)  
  
  A=c(offset$lsmean/mean(offset$lsmean),offset.rand$lsmean/mean(offset.rand$lsmean))
  Ylim=c(min(A),max(A))
  with(offset,plot(lsmean/mean(lsmean),ylim=Ylim,type='l',col=1,ylab="relative catch rate",cex.lab=1.35))
  with(offset.rand,lines(lsmean/mean(lsmean),col=2))    
  legend("bottomleft",c("main effect","random effect"),lty=1,col=1:2,bty='n',cex=1.5)  
}
Compare.vess.random=function(RESHAPED)   #function for testing models
{   
  RESHAPED=subset(RESHAPED,Catch.Target>0)
  RESHAPED$log.Effort=log(RESHAPED$Km.Gillnet.Days.c)
  RESHAPED$log.Catch=log(RESHAPED$Catch.Target)
  RESHAPED$FINYEAR=as.factor(RESHAPED$FINYEAR)
  RESHAPED$VESSEL=as.factor(RESHAPED$VESSEL)
  
  #model
  MODEL_RESHAPED <- glm(log.Catch ~ FINYEAR +VESSEL+ offset(log.Effort), data=RESHAPED, family=gaussian, maxit=500)  
  MODEL_RESHAPED.random=lmer(log.Catch ~ FINYEAR +offset(log.Effort) +(1 |VESSEL), data=RESHAPED)
  
  #Plot annual cpue
  Preds_RESHAPED=lsmeans(MODEL_RESHAPED,specs=c("FINYEAR"))
  Preds_RESHAPED.random=lsmeans(MODEL_RESHAPED.random,specs=c("FINYEAR"))
  offset=summary(Preds_RESHAPED)
  offset_RESHAPED=as.data.frame(offset[c("lsmean")])
  offset=summary(Preds_RESHAPED.random)
  offset_RESHAPED.random=as.data.frame(offset[c("lsmean")])
  fn3(offset_RESHAPED,offset_RESHAPED.random)
  mtext("Reshaped data, indicative vessels",side=3,line=-2)  
}


#Information theory functions
Inc.Dev.Exp=function(A,B)round((100*A/B)-100) #increase in deviance explained from one model to next

fn.AIC.BIC.ratio=function(DAT)  #Function for calculating AIC and BIC weights and ratios, deviance explained and store best model
{
  MIN=min(DAT)
  Delta=DAT-MIN
  Like.model.give.dat=exp(-Delta/2)
  Weight=Like.model.give.dat/sum(Like.model.give.dat)
  id=which(Weight==max(Weight))
  names(Weight)=NULL
  Evidence.ratio=outer(Weight[id],Weight, "/")
  return(list(Best.Mod=id,Delta=Delta,Like=Like.model.give.dat,Weight=Weight,Evidence.ratio_how.much.better=Evidence.ratio))
}

fn.improved.dev=function(DAT,Red.th,NAMES)  #function for calculating improved in deviance
{
  Dev.red=DAT[1:(length(DAT)-1)]-DAT[2:length(DAT)]
  RED=DAT[1:(length(DAT)-1)]*Red.th/100
  Keep.term=ifelse(Dev.red>RED,"KEEP","DROP")
  names(Keep.term)=NAMES[-1]
  return(Keep.term=Keep.term)
}



#function for applying best model and constructing index
fn.stand.cpue=function(DataFile,Formula.cpue,Formula.cpue.pos,MESH,Effort.calc)
{
  #A. Standardise catch
  if(MESH=="YES")
  {
    DataFile=subset(DataFile,!is.na(mesh))
    ALLvars=all.vars(Formula.cpue)
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],"mesh",paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
    
    ALLvars=all.vars(Formula.cpue.pos)
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],"mesh",paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue.pos=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
  }
  
  #drop levels not occurring in data
  for(f in 1:ncol(DataFile))if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
  
  #log effort
  if(Effort.calc=="Km.Gillnet.Days")DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c)
  if(Effort.calc=="Km.Gillnet.Hours")DataFile$log.Effort=log(DataFile$Km.Gillnet.Hours.c)

  
  #Find.Blk.10=match("block10",names(DataFile))
  Find.Blk.10=NA
  
  #Create binary dataset (0/1 for CPUE)
  BiData <- DataFile
  Id.Ktch=match("Catch.Target",names(DataFile))
  BiData[,Id.Ktch] <- as.numeric(DataFile[,Id.Ktch]>0)
  
  #Keep vessels with a minimum number of observations
  a=with(BiData,table(FINYEAR,as.character(VESSEL)))
  a=colSums(a)
  a=subset(a,a>Min.rec.ves)
  BiData=subset(BiData,VESSEL%in%names(a))
  
  
  #convert to factor
  BiData$FINYEAR=factor(BiData$FINYEAR,levels=sort(unique(BiData$FINYEAR)))
  BiData$MONTH=factor(BiData$MONTH,levels=1:12)  
  BiData$BLOCKX=as.factor(BiData$BLOCKX)
  BiData$VESSEL=as.factor(BiData$VESSEL)
  if(!is.na(Find.Blk.10)) BiData$block10=as.factor(BiData$block10)
  if(MESH=="YES")BiData$mesh=as.factor(BiData$mesh)
  
  #Create positive (non-zero) dataset
  PosData<-DataFile[DataFile[,Id.Ktch]>0,]  
  PosData$log.Catch=log(PosData$Catch.Target)    
  
  #Keep vessels with a minimum number of observations
  a=with(PosData,table(FINYEAR,as.character(VESSEL)))
  a=colSums(a)
  a=subset(a,a>Min.rec.ves)
  PosData=subset(PosData,VESSEL%in%names(a))
  
  
  #convert to factor
  PosData$FINYEAR=factor(PosData$FINYEAR,levels=sort(unique(PosData$FINYEAR)))
  PosData$MONTH=factor(PosData$MONTH,levels=1:12)  
  PosData$BLOCKX=as.factor(PosData$BLOCKX)  
  PosData$VESSEL=as.factor(PosData$VESSEL)  
  if(!is.na(Find.Blk.10)) PosData$block10=as.factor(PosData$block10)
  if(MESH=="YES")PosData$mesh=as.factor(PosData$mesh)
  
  if(length(levels(PosData$BLOCKX))==1)
  {
    ALLvars=all.vars(Formula.cpue)
    ALLvars=subset(ALLvars,!ALLvars=="BLOCKX")
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
    
    ALLvars=all.vars(Formula.cpue.pos)
    ALLvars=subset(ALLvars,!ALLvars=="BLOCKX")
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue.pos=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
  }
  
  if(length(unique(PosData$FINYEAR))>1)
  {
    #binomial GLM
    GLMbi <- glm(Formula.cpue, data=BiData, family="binomial",maxit=500)
    
    #Lognormal GLM
    GLMlog <- glm(Formula.cpue.pos, data=PosData, family=gaussian, maxit=500)  
    
    
    #B. Return outputs
    LisT=list(GLMlog=GLMlog,GLMbi=GLMbi,PosData=PosData,BiData=BiData)
    
    return(LisT)
    
  }
  
}
see.BC.fit.to.data=function(GLM,RespVar)
{
  BiData=GLM$BiData
  PosData=GLM$PosData 
  Log.GLMbi=GLM$GLMbi
  GLMlog=GLM$GLMlog 
  
  BiData=subset(BiData,VESSEL%in%names(dummy.coef(Log.GLMbi)$VESSEL))
  PosData=subset(PosData,VESSEL%in%names(dummy.coef(GLMlog)$VESSEL))
  
  par(mfcol=c(3,1),mai=c(.6,.75,.3,.1),oma=c(1,1,1,1),las=1,mgp=c(2,.5,0))
  
  if(!is.null(Log.GLMbi))
  {
    BiData$Pred=round(predict(Log.GLMbi,newdata=BiData,type="response"))
    a=aggregate(Catch.Target~FINYEAR,BiData,mean)
    b=aggregate(Pred~FINYEAR,BiData,mean)
    plot(a[,2],pch=19,ylab="Prob of catch",xlab="",xaxt='n',cex.lab=1.5,ylim=c(0,1.1))
    lines(b[,2],col=2,lwd=2)
    axis(1,seq(1,nrow(a),1),F,tck=-0.02)
    axis(1,seq(1,nrow(a),5),F,tck=-0.05)
  }
  if(is.null(Log.GLMbi)) 
  {
    plot.blank("Prob of catch") 
  }
  mtext("Presence/Absence data",3,cex=1.5)
  #if(show.not.bias.corr=="NO")legend("bottomright",c("observed","predicted"),lty=c(0,1),pch=c(19,""),col=1:2,bty='n',cex=1.5)
  #if(show.not.bias.corr=="YES")legend("bottomright",c("obs","pred.bias.corr","pred.no.bias.corr"),lty=c(0,1,1),pch=c(19,"",""),col=1:3,bty='n',cex=1.5)
  
  #Pos.pred=predict(GLMlog,type="response",se.fit =T)
  
  
  PosData$Pred=predict(GLMlog,newdata=PosData,type="response")
  #   PosData$Pred=BiasCor.fn(Pos.pred$fit,Pos.pred$se.fit)
  #   PosData$Pred.no.bias.corr=exp(Pos.pred$fit)
  #   PosData$Pred.SE.up=exp(Pos.pred$fit+Pos.pred$se.fit)
  #   PosData$Pred.SE.low=exp(Pos.pred$fit-Pos.pred$se.fit)
  
  
  if(RespVar=="catch")
  {
    a=aggregate(log.Catch~FINYEAR,PosData,sum)
    b=aggregate(Pred~FINYEAR,PosData,sum)
    #     d=aggregate(Pred.no.bias.corr/1000~FINYEAR,PosData,sum)
    #     b.SE.up=aggregate(Pred.SE.up/1000~FINYEAR,PosData,sum)
    #     b.SE.low=aggregate(Pred.SE.low/1000~FINYEAR,PosData,sum)
    #     
    LABS="Log total catch"
  }
  
  if(RespVar=="cpue")
  {
    a=aggregate(cpue~FINYEAR,PosData,mean)
    b=aggregate(Pred~FINYEAR,PosData,mean)
    #     d=aggregate(Pred.no.bias.corr~FINYEAR,PosData,mean)
    #     b.SE.up=aggregate(Pred.SE.up~FINYEAR,PosData,mean)
    #     b.SE.low=aggregate(Pred.SE.low~FINYEAR,PosData,mean)
    
    LABS="CPUE (kg/km.gn.day)"
  }
  RANG=c(min(c(a[,2],b[,2])),max(c(a[,2],b[,2])))   
  # RANG=c(0,max(c(a[,2],b.SE.up[,2]))) 
  plot(a[,2],pch=19,ylab=LABS,xlab="",xaxt='n',ylim=RANG,cex.lab=1.5)
  lines(b[,2],col=2,lwd=2)
  #   if(show.not.bias.corr=="YES")lines(d[,2],col=3,lwd=2)
  #   lines(b.SE.up[,2],col=2,lty=2,lwd=2.5)
  #   lines(b.SE.low[,2],col=2,lty=2,lwd=2.5)
  mtext("Positive catch data",3,cex=1.5)
  mtext("Financial year",1,cex=1.5,line=2)
  axis(1,seq(1,nrow(a),1),F,tck=-0.02)
  axis(1,seq(1,nrow(a),5),a[,1][seq(1,nrow(a),5)],tck=-0.05)  
  
  plot(PosData$log.Catch,PosData$Pred)
  lines(PosData$log.Catch,PosData$log.Catch,col=2)
}
BiasCor.fn=function(Median,SE) biasCorr <- exp(Median+(SE^2)/2) #function for bias corrected mean in normal space
Pos.Diag.fn=function(MODEL,SPECIES)   #function for positive catch diagnostics
{
  RES=MODEL$residuals   #residuals
  Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
  PRED=predict(MODEL)
  M=1.5
  
  qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="",xlab="")
  qqline(RES, col = 'grey40',lwd=1.5,lty=2)
  mtext(SPECIES,3,outer=F,line=0.25,cex=1.3)
  if(i==1) mtext("Residuals",2,outer=F,line=2,las=3,cex=M)
  if(i==2) mtext("                        Theoretical quantiles",1,outer=F,line=2,cex=M)
  
  hist(Std.RES,xlim=c(-5,5),ylab="",xlab="",main="",col="grey",breaks=50)
  box()
  if(i==1) mtext("Frequency",2,outer=F,line=2.5,las=3,cex=M)
  if(i==2) mtext("                      Standardised residuals",1,outer=F,line=2,cex=M)
  
  plot(PRED,Std.RES,ylab="",xlab="")
  abline(0,0,lwd=1.5,lty=2,col='grey40')
  if(i==1) mtext("Standardised residuals",2,outer=F,line=2,las=3,cex=M)
  if(i==2) mtext("                         Fitted values",1,outer=F,line=2,cex=M)
  
}
fn.pred.month=function(MOD,SP)
{
  LogN=MOD$GLMlog
  This=c("MONTH2","MONTH3","MONTH4","MONTH5","MONTH6","MONTH7","MONTH8","MONTH9",
         "MONTH10","MONTH11","MONTH12") 
  COFS=coef(LogN)[match(This,names(coef(LogN)))]
  plot(2:12,COFS,type='o',cex=2,pch=19,ylab="GLM Coefficient",xlab="Month") 
  legend("top",SP,bty='n',cex=2)
}
fn.blk.zone=function(BLK,dat)   #get blocks per zone
{
  da=dat[,match(c(BLK,"zone"),names(dat))]
  da=da[!duplicated(da[,1]),]
  rownames(da)=NULL
  return(da)
}

#function for extracting term significance and deviance explained
Anova.and.Dev.exp=function(GLMbi,GLMlog)
{
  #Anovas
  Anova.tab.bi=anova(GLMbi, test = "Chisq")
  Anova.tab.pos=anova(GLMlog, test = "Chisq")
  
  #Deviance explained
  #By each term
  n=2:length(Anova.tab.bi$Deviance)
  Term.dev.exp.bi=100*(Anova.tab.bi$Deviance[n]/GLMbi$null.deviance)
  names(Term.dev.exp.bi)=rownames(Anova.tab.bi)[n]
  
  n=2:length(Anova.tab.pos$Deviance)
  Term.dev.exp.pos=100*(Anova.tab.pos$Deviance[n]/GLMlog$null.deviance)
  names(Term.dev.exp.pos)=rownames(Anova.tab.pos)[n]
  
  #By full model
  Dev.exp.bi=sum(Term.dev.exp.bi)
  Dev.exp.pos=sum(Term.dev.exp.pos)
  
  #Dev.exp.bi=fn.Dev.Exp(GLMbi)
  #Dev.exp.pos=fn.Dev.Exp(GLMlog)
  
  LisT=list(Anova.bi=Anova.tab.bi,Anova.pos=Anova.tab.pos,
            Dev.exp.bi=Dev.exp.bi,Dev.exp.pos=Dev.exp.pos,
            Term.Dev.exp.bi=Term.dev.exp.bi,Term.Dev.exp.pos=Term.dev.exp.pos)
  
  return(LisT)
  
}
Comb.anova.and.dev=function(d,NM)
{
  #Bi part
  ANOVA.bi=as.data.frame.matrix(d$Anova.bi)
  Term.bi=data.frame(Percent.dev.exp=d$Term.Dev.exp.bi)
  Table.bi=ANOVA.bi[-1,match(c("Deviance","Pr(>Chi)"),names(ANOVA.bi))]
  Term.bi=Term.bi[match(rownames(Term.bi), rownames(Table.bi)),]
  Table.bi=cbind(Table.bi,Term.bi)
  names(Table.bi)[match("Term.bi",names(Table.bi))]="Percent.dev.exp"
  All.bi=Table.bi[1,]
  rownames(All.bi)="Model"
  All.bi[,1:ncol(All.bi)]=NA
  All.bi$Percent.dev.exp=d$Dev.exp.bi
  Table.bi=rbind(Table.bi,All.bi)
  colnames(Table.bi)=paste("bi",colnames(Table.bi),sep="_")
  
  #Pos part
  ANOVA.pos=as.data.frame.matrix(d$Anova.pos)
  Term.pos=data.frame(Percent.dev.exp=d$Term.Dev.exp.pos)
  Table.pos=ANOVA.pos[-1,match(c("Deviance","Pr(>Chi)"),names(ANOVA.pos))]
  Term.pos=Term.pos[match(rownames(Term.pos), rownames(Table.pos)),]
  Table.pos=cbind(Table.pos,Term.pos)
  names(Table.pos)[match("Term.pos",names(Table.pos))]="Percent.dev.exp"
  All.pos=Table.pos[1,]
  rownames(All.pos)="Model"
  All.pos[,1:ncol(All.pos)]=NA
  All.pos$Percent.dev.exp=d$Dev.exp.pos
  Table.pos=rbind(Table.pos,All.pos)
  colnames(Table.pos)=paste("pos",colnames(Table.pos),sep="_")
  
  Table=cbind(data.frame(Species=c(NM,rep("",nrow(Table.bi)-1)),Terms=row.names(Table.bi)),Table.bi,Table.pos)
  Table[,3:ncol(Table)]=round(Table[,3:ncol(Table)],3)
  Table$"bi_Pr(>Chi)"=ifelse(Table$"bi_Pr(>Chi)"<0.001,"<0.001",Table$"bi_Pr(>Chi)")
  Table$"pos_Pr(>Chi)"=ifelse(Table$"pos_Pr(>Chi)"<0.001,"<0.001",Table$"pos_Pr(>Chi)")
  Table[is.na(Table)] <- ""
  return(Table)  
}
fn.exp.anova=function(DAT,NM,WHAT)
{
  Bi.tab=as.data.frame(anova(DAT$GLMbi, test = "Chisq"))
  Bi.tab=cbind(TERM=rownames(Bi.tab),Bi.tab)
  fn.word.table(WD=getwd(),TBL=Bi.tab,Doc.nm=paste(NM,"_Bi_",WHAT,sep=""),caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  Log.tab=anova(DAT$GLMlog, test = "Chisq")
  Log.tab=as.data.frame(anova(Stand.BaseCase[[i]]$GLMlog, test = "Chisq"))
  Log.tab=cbind(TERM=rownames(Log.tab),Log.tab)
  fn.word.table(WD=getwd(),TBL=Log.tab,Doc.nm=paste(NM,"_Pos_",WHAT,sep=""),caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}
fn.get.Significance=function(MOD,SPEC) 
{
  LogN=MOD$LogN
  Bi=MOD$Bi  
  Anova.bi=anova(Bi, test = "Chisq")
  Anova.LogN=anova(LogN, test = "Chisq")
  return(list(Bi=Anova.bi,LogN=Anova.LogN))
}

#function for sensitivity test: predicting area-weighted annual index
fn.stand.cpue.sens=function(MOD,FORMULA,Apply.area.w,Ar)
{
  GLMbi=MOD$GLMbi
  GLMlog=MOD$GLMlog
  BiData=MOD$BiData
  PosData=MOD$PosData
  Formula.cpue=FORMULA$Bi
  Formula.cpue.pos=FORMULA$Log
  
  #Get model terms
  bi.terms=all.vars(Formula.cpue)
  pos.terms=all.vars(Formula.cpue.pos)
  
  Find.Blk.10=match("block10",pos.terms)
  
  #Predict year-block
  if(is.na(Find.Blk.10))
  {
    bi.terms=bi.terms[-match(c("Catch.Target", "FINYEAR","BLOCKX"),bi.terms)]
    pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","BLOCKX"),pos.terms)]  
    Bi.dat.pred=fn.match.terms(bi.terms,BiData,GLMbi)
    Pos.dat.pred=fn.match.terms(pos.terms,PosData,GLMlog)
  }
  if(!is.na(Find.Blk.10))
  {
    bi.terms= bi.terms[-match(c("Catch.Target", "FINYEAR","block10"),bi.terms)]
    pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","block10"),pos.terms)]   
    Bi.dat.pred=fn.match.terms.daily(bi.terms,BiData,GLMbi)
    Pos.dat.pred=fn.match.terms.daily(pos.terms,PosData,GLMlog)
  }
  Bi.dat.new=Bi.dat.pred$Value
  Pos.dat.new=Pos.dat.pred$Value
  Bi.dat.new$prob=predict(GLMbi,Bi.dat.new,type="response")
  Pos.dat.new$pos=predict(GLMlog,Pos.dat.new,type="response")  
  if(is.na(Find.Blk.10))
  {
    Index=merge(Bi.dat.new[,match(c("FINYEAR","BLOCKX","prob"),names(Bi.dat.new))],
                Pos.dat.new[,match(c("FINYEAR","BLOCKX","pos"),names(Pos.dat.new))],by=c("FINYEAR","BLOCKX"))    
  }
  if(!is.na(Find.Blk.10))
  {
    Index=merge(Bi.dat.new[,match(c("FINYEAR","block10","prob"),names(Bi.dat.new))],
                Pos.dat.new[,match(c("FINYEAR","block10","pos"),names(Pos.dat.new))],by=c("FINYEAR","block10"))    
  }
  Index$Pos=exp(Index$pos)
  
  #Combine prob of catch and pos catch
  Index$Index=Index$Pos*Index$prob
  
  #Weight by area
  if(Apply.area.w=="No") Ar$Fish.Area=1
  if(is.na(Find.Blk.10)) Ar=subset(Ar,BLOCKX%in%unique(Index$BLOCKX),select=c(BLOCKX,Fish.Area))
  if(!is.na(Find.Blk.10)) Ar=subset(Ar,block10%in%unique(Index$block10),select=c(block10,Fish.Area))
  Ar$Fish.Area=Ar$Fish.Area/sum(Ar$Fish.Area)  #out of total available area
  if(is.na(Find.Blk.10)) Index=merge(Index,Ar,by="BLOCKX",all.x=T) 
  if(!is.na(Find.Blk.10)) Index=merge(Index,Ar,by="block10",all.x=T)
  Index$Index.w=Index$Index*Index$Fish.Area
  
  #Aggregate by year
  I_y.w=aggregate(Index.w~FINYEAR,Index,sum,na.rm=T)
  
  I_y.w$cpue=(I_y.w$Index.w)/exp(Lg.Efrt)   #express as catch rate  
  I_y.w$FINYEAR=as.character(I_y.w$FINYEAR)
  I_y.w=I_y.w[order(I_y.w$FINYEAR),]
  
  return(I_y.w)  
}


#functions for calculating CI thru MC simulations
get.terms=function(FORMULA)
{
  a=gsub(" ", "", unlist(strsplit(as.character(FORMULA)[3],"[+]")), fixed = TRUE)
  a=a[-(1:2)]
  IDD=match("offset(log.Effort)",a)
  a[IDD]="log.Effort"
  return(a)
}

fn.cof=function(mod,End,what)
{
  a=substr(names(coef(mod)),1,End)
  a=which(a==what)
  return(coef(mod)[a])
}
fn.match.terms=function(TERMS,DAT,MOD)
{
  YRs=unique(as.character(levels(DAT$FINYEAR)))
  BLKs=unique(as.character(levels(DAT$BLOCKX)))
  ID=match(TERMS,names(DAT))
  Clas=mapply(class,DAT[,ID])
  Value=as.data.frame(matrix(ncol=length(Clas),nrow=1))   
  names(Value)=names(Clas)
  for(w in 1:length(Clas))
  {
    if(Clas[w]=='factor') 
    {
      xx=match(names(Clas[w]),names(DAT))     
      
      if(Sel.lev=='most.common') 
      {
        Tabla=sort(table(DAT[,xx]))
        THIS=names(Tabla[length(Tabla)])        
      }
      
      if(Sel.lev=='average')
      {
        NN=names(Clas[w])
        COF=fn.cof(MOD,nchar(NN),NN)        
        names(COF)=substr(names(COF),(nchar(NN)+1),50)
        MEAN.level=mean(COF)
        THIS=names(sort(abs(COF-MEAN.level)))[1]        
      }
      Value[,w]=factor(THIS,levels=levels(DAT[,xx]))
    }
    
    if(Clas[w]=='numeric') Value[,w]=mean(DAT[,match(names(Clas[w]),names(DAT))],na.rm=T)
    if(TERMS[w]=='log.Effort') Value[,w]=Lg.Efrt
  }
  DD=cbind(expand.grid(FINYEAR=YRs,BLOCKX=BLKs),Value)
  return(list(Value=DD,Clas=Clas))
}
fn.match.terms.daily=function(TERMS,DAT,MOD)
{
  YRs=unique(as.character(levels(DAT$FINYEAR)))
  BLK_10=unique(as.character(levels(DAT$block10)))
  ID=match(TERMS,names(DAT))
  Clas=mapply(class,DAT[,ID])
  Value=as.data.frame(matrix(ncol=length(Clas),nrow=1))   
  names(Value)=names(Clas)
  for(w in 1:length(Clas))
  {   
    if(Clas[w]=='factor') 
    {
      xx=match(names(Clas[w]),names(DAT))     
      
      if(Sel.lev=='most.common') 
      {
        Tabla=sort(table(DAT[,xx]))
        THIS=names(Tabla[length(Tabla)])        
      }
      
      if(Sel.lev=='average')
      {
        NN=names(Clas[w])
        COF=fn.cof(MOD,nchar(NN),NN)        
        names(COF)=substr(names(COF),(nchar(NN)+1),50)
        MEAN.level=mean(COF)
        THIS=names(sort(abs(COF-MEAN.level)))[1]        
      }
      Value[,w]=factor(THIS,levels=levels(DAT[,xx]))
    }
    
    if(Clas[w]=='numeric') Value[,w]=mean(DAT[,match(names(Clas[w]),names(DAT))],na.rm=T)
    if(TERMS[w]=='log.Effort') Value[,w]=Lg.Efrt
  }
  DD=cbind(expand.grid(FINYEAR=YRs,block10=BLK_10),Value)
  return(list(Value=DD,Clas=Clas))
}

#Function for Monte Carlo index
fn.MC.cpue=function(MOD,Formula.cpue,Formula.cpue.pos,niter,No.MC.Bi,Area.w,BLKs,MESH)
{
  if(MESH=="YES")
  {
    ALLvars=all.vars(Formula.cpue)
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],"mesh",paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
    
    ALLvars=all.vars(Formula.cpue.pos)
    TERMS=c(ALLvars[2:(length(ALLvars)-1)],"mesh",paste("offset(",ALLvars[length(ALLvars)],")",sep=""))
    Formula.cpue.pos=as.formula(paste(ALLvars[1],"~",paste(TERMS,collapse="+")))
  }
  
  #A. Get models and data
  GLMbi=MOD$GLMbi
  GLMlog=MOD$GLMlog
  BiData=MOD$BiData
  PosData=MOD$PosData
  
  
  #B. Build confidence bounds
  
  #B.1. Draw random sample of coefficients from multivariate normal distribution
  #Binomial part
  Covar.bi=as.matrix(vcov(GLMbi))
  set.seed(999);Bi.pars.rand=rmvnorm(niter,mean=coef(GLMbi),sigma=Covar.bi)    
  
  #Positive part
  Covar.pos=as.matrix(vcov(GLMlog))
  set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(GLMlog),sigma=Covar.pos)    
  
  
  #B.2 Predict new set of parameters
  #get model terms
  bi.terms=all.vars(Formula.cpue)
  pos.terms=all.vars(Formula.cpue.pos)
  
  #Find.Blk.10=match("block10",bi.terms)
  Find.Blk.10=NA
  
  if(is.na(Find.Blk.10))
  {
    if(BLKs=="byzone")
    {
      BiData=subset(BiData,BLOCKX%in%BLks)
      PosData=subset(PosData,BLOCKX%in%BLks)
    }
    
    bi.terms=bi.terms[-match(c("Catch.Target", "FINYEAR","BLOCKX"),bi.terms)]
    pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","BLOCKX"),pos.terms)]  
    Bi.dat.pred=fn.match.terms(bi.terms,BiData,GLMbi)
    Pos.dat.pred=fn.match.terms(pos.terms,PosData,GLMlog)
    
    if(BLKs=="byzone")
    {
      Bi.dat.pred$Value=subset(Bi.dat.pred$Value,BLOCKX%in%BLks)
      Pos.dat.pred$Value=subset(Pos.dat.pred$Value,BLOCKX%in%BLks)
    }
    
  }
  
  if(!is.na(Find.Blk.10))
  {
    if(BLKs=="byzone")
    {
      BiData=subset(BiData,block10%in%BLks)
      PosData=subset(PosData,block10%in%BLks)
    }
    
    bi.terms= bi.terms[-match(c("Catch.Target", "FINYEAR","block10"),bi.terms)]
    pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","block10"),pos.terms)]   
    Bi.dat.pred=fn.match.terms.daily(bi.terms,BiData,GLMbi)
    Pos.dat.pred=fn.match.terms.daily(pos.terms,PosData,GLMlog)
    
    if(BLKs=="byzone")
    {
      Bi.dat.pred$Value=subset(Bi.dat.pred$Value,block10%in%BLks)
      Pos.dat.pred$Value=subset(Pos.dat.pred$Value,block10%in%BLks)
    }
    
  }
  
  
  #predict MC data
  dummy.GLMbi=GLMbi
  dummy.GLMlog=GLMlog
  YrSS=sort(unique(as.character(Bi.dat.pred$Value$FINYEAR)))
  MC.preds=matrix(nrow=niter,ncol=length(YrSS))
  colnames(MC.preds)=YrSS
  
  if(No.MC.Bi=="YES") Bi.dat.pred$Value$Pred.bi=predict(GLMbi,newdata=Bi.dat.pred$Value, type="response")
  
  
  for(n in 1:niter)
  {
    #Get predictions
    #Bidata
    if(!No.MC.Bi=="YES")
    {
      dummy.GLMbi$coefficients=Bi.pars.rand[n,]
      Bi.dat.pred$Value$Pred.bi=predict(dummy.GLMbi,newdata=Bi.dat.pred$Value, type="response")      
    }
    
    #Posdata
    dummy.GLMlog$coefficients=Pos.pars.rand[n,]
    a=predict(dummy.GLMlog,newdata=Pos.dat.pred$Value, type="response",se.fit=T)
    Pos.dat.pred$Value$Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    
    
    #Combine prob of catch and positive catch
    if(is.na(Find.Blk.10)) Index=merge(subset(Bi.dat.pred$Value,select=c(FINYEAR,BLOCKX,Pred.bi)),
                                       subset(Pos.dat.pred$Value,select=c(FINYEAR,BLOCKX,Pred)), by=c("FINYEAR","BLOCKX"),all=T)
    
    if(!is.na(Find.Blk.10)) Index=merge(subset(Bi.dat.pred$Value,select=c(FINYEAR,block10,Pred.bi)),
                                        subset(Pos.dat.pred$Value,select=c(FINYEAR,block10,Pred)), by=c("FINYEAR","block10"),all=T)
    
    Index$Index=Index$Pred.bi*Index$Pred
    
    
    #Calculate annual index, weighted by area
    if(is.na(Find.Blk.10)) Area.w=subset(Area.w,BLOCKX%in%unique(Index$BLOCKX),select=c(BLOCKX,Fish.Area))
    if(!is.na(Find.Blk.10)) Area.w=subset(Area.w,block10%in%unique(Index$block10),select=c(block10,Fish.Area))
    
    Area.w$Fish.Area=Area.w$Fish.Area/sum(Area.w$Fish.Area)  #out of total available area
    
    if(is.na(Find.Blk.10)) Index=merge(Index,Area.w,by="BLOCKX",all.x=T) 
    if(!is.na(Find.Blk.10)) Index=merge(Index,Area.w,by="block10",all.x=T)
    
    Index$Index.w=Index$Index*Index$Fish.Area
    I_y.w=aggregate(Index.w~FINYEAR,Index,sum, na.action = na.pass)
    
    
    #Store simulations
    MC.preds[n,]=I_y.w$Index.w
  }
  rm(a)
  
  
  #B.5 Get summary stats
  MEAN=colMeans(MC.preds,na.rm=T)
  SD=apply(MC.preds,2,sd,na.rm=T)
  CV=SD/MEAN
  LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T))
  UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T))
  
  Stats=data.frame(FINYEAR=colnames(MC.preds),MEAN=MEAN,SD=SD,CV=CV,LOW=LOW,UP=UP)
  Stats=Stats[order(Stats$FINYEAR),]
  rownames(Stats)=NULL
  return(Stats)
}


#function for showing term effect
fn.show.coef=function(MOD,SP,This,ShowLevel,WHRE,CEX,Pt.cex)
{
  a=as.data.frame(coef(summary(MOD)))
  Check=substr(rownames(a),1,nchar(This))
  That=which(Check%in%This) 
  a=a[That,] 
  NMS=as.character(substr(rownames(a),nchar(This)+1,100))
  if(ShowLevel=="NO") NMS=as.character(2:(length(NMS)+1))
  
  COFS=a$Estimate
  SE=a$'Std. Error'
  names(COFS)=names(SE)=NMS
  #   if(ShowLevel=="YES")
  #   {
  #     IDss=match(sort(as.numeric(NMS)),NMS)
  #     COFS=COFS[IDss]
  #     SE=SE[IDss]
  #     NMS=NMS[IDss]
  #   }
  
  YLIM=c(min(COFS-SE),max(COFS+SE))
  plot(1:length(NMS),COFS,cex=Pt.cex,pch=19,ylab="",xlab="",ylim=YLIM,xaxt='n',cex.axis=1.15)
  if(ShowLevel=="YES")axis(1,1:length(NMS),labels=NMS,cex.axis=CEX)
  if(ShowLevel=="NO")axis(1,1:length(NMS),labels=F)
  segments(1:length(COFS),COFS,1:length(COFS),COFS+SE,lwd=2)
  segments(1:length(COFS),COFS,1:length(COFS),COFS-SE,lwd=2)
  legend(WHRE,SP,bty='n',cex=2)
}
fn.show.coef.spatial=function(MOD,SP,This,dat)
{
  a=as.data.frame(coef(summary(MOD)))
  Check=substr(rownames(a),1,nchar(This))
  That=which(Check%in%This) 
  a=a[That,] 
  NMS=as.character(substr(rownames(a),nchar(This)+1,100))
  IDi=match(This,names(dat))
  dat=dat[!duplicated(dat[IDi]),]
  names(dat)[IDi]="BLK"
  dat$LAT=-as.numeric(substr(dat$BLK,1,2))
  dat$LONG=100+as.numeric(substr(dat$BLK,3,4))
  All.blks=unique(dat[,IDi])
  Ref.blk=All.blks[which(!as.character(All.blks)%in%NMS)]
  COFS=c(0,a$Estimate)
  MAT=data.frame(BLK=c(Ref.blk,NMS),Coef=COFS)
  MAT=merge(MAT,dat,by="BLK",all.x=T)
  all.lat=seq(min(MAT$LAT),max(MAT$LAT))
  Misin.lat=all.lat[which(!all.lat%in%MAT$LAT)]
  all.lon=seq(min(MAT$LONG),max(MAT$LONG))
  Misin.lon=all.lon[which(!all.lon%in%MAT$LONG)]
  if(length(Misin.lon)>0 & length(Misin.lat)>0)
  {
    Add=data.frame(BLK=paste(-Misin.lat,Misin.lon-100,sep=""),Coef=NA,LONG=Misin.lon, LAT=Misin.lat)
    MAT=rbind(MAT,Add)
  }
  
  if(length(Misin.lon)==0 & length(Misin.lat)>0)
  {
    Add=data.frame(BLK=paste(-Misin.lat,114-100,sep=""),Coef=NA,LONG=114, LAT=Misin.lat)
    MAT=rbind(MAT,Add)
  }
  if(length(Misin.lon)>0 & length(Misin.lat)==0)
  {
    Add=data.frame(BLK=paste(34,Misin.lon-100,sep=""),Coef=NA,LONG=Misin.lon, LAT=-34)
    MAT=rbind(MAT,Add)
  }
  
  
  dat=reshape(MAT[match(c("Coef","LONG","LAT"),names(MAT))],
              idvar="LONG",timevar="LAT",v.names="Coef", direction="wide")
  names(dat)[2:ncol(dat)]=substr(names(dat)[2:ncol(dat)],start=6,stop=15)
  dat=dat[order(dat$LONG),]
  
  Breaks=sort(quantile(unlist(as.matrix(dat[,2:ncol(dat)])),probs=seq(0,1,1/numInt),na.rm=T))
  
  Lon=dat$LONG
  d=dat[,2:ncol(dat)]
  Sort.Lat=as.character(sort(as.numeric(colnames(d))))
  d=d[,match(Sort.Lat,colnames(d))]
  d=as.matrix(d)
  Lat=as.numeric(colnames(d))
  
  Lat=Lat-.5
  Lon=Lon+.5
  
  image(Lon,Lat,d,col =couleurs,breaks=Breaks,xaxt='n',
        yaxt='n',ylab="",xlab="",main="",cex.main=1.25,ylim=c(-36,-26),xlim=c(113,129))
  box()
  axis(1,113:129,F,tck=-0.015)
  axis(2,-36:-26,F,tck=-0.015)
  if(!SP=="")legend("topright",SP,bty='n',cex=1.5)
  
}
fn.show.coef.spatial.daily=function(MOD,SP,This,dat,where)
{
  a=as.data.frame(coef(summary(MOD)))
  Check=substr(rownames(a),1,nchar(This))
  That=which(Check%in%This) 
  a=a[That,] 
  NMS=as.numeric(substr(rownames(a),nchar(This)+1,100))
  dat=dat[!duplicated(dat[,1]),1]
  COFS=c(0,a$Estimate)
  names(COFS)=c(dat[which(!dat%in%NMS)],NMS)
  MAT=data.frame(BLK=names(COFS),Coef=COFS)
  MAT$BLK=as.numeric(as.character(MAT$BLK))
  MAT=merge(MAT,BlOCK_10,by.x="BLK",by.y="block10",all.y=T)
  MAT=MAT[order(MAT$LAT),]
  dat=reshape(MAT[match(c("Coef","LONG","LAT"),names(MAT))],
              idvar="LONG",timevar="LAT",v.names="Coef", direction="wide")
  names(dat)[2:ncol(dat)]=substr(names(dat)[2:ncol(dat)],start=6,stop=15)
  dat=dat[order(dat$LONG),]
  Breaks=sort(quantile(unlist(as.matrix(dat[,2:ncol(dat)])),probs=seq(0,1,1/numInt),na.rm=T))
  Lon=dat$LONG
  d=dat[,2:ncol(dat)]
  d=as.matrix(d)
  Lat=as.numeric(colnames(d))     
  image(Lon,Lat,d,col =couleurs,breaks=Breaks,xaxt='n',
        yaxt='n',ylab="",xlab="",main="",cex.main=1.25,ylim=c(-36,-26),xlim=c(113,129))
  box()
  axis(1,113:129,F,tck=-0.015)
  axis(2,-36:-26,F,tck=-0.015)
  legend(where,SP,bty='n',cex=1.75)
}

#function for expressing in relative terms     
fn.relative=function(D)
{
  Mn=mean(D$MEAN,na.rm=T)
  D$LOW=D$LOW/Mn
  D$UP=D$UP/Mn
  D$SD=D$SD/Mn
  D$MEAN=D$MEAN/Mn
  return(D)
}
fn.relative.fol.nom=function(D,what)
{
  if(!is.null(D)) D$cpue=D$cpue/mean(D$cpue)
  return(D)
}

#function for removing effort creep post construction of index
fn.remove.eff.creep=function(d,VRS)
{
  Power=Fish.pow.inc[match(d$FINYEAR,names(Fish.pow.inc))]
  d[match(VRS,names(d))]=d[match(VRS,names(d))]/Power
  return(d)
}

#function for converting standardised catch to catch rate
fn.rate=function(d,VRS)
{
  d[match(VRS,names(d))]=d[match(VRS,names(d))]/Stand.eff
  return(d)
}

#functions for plotting index
PLOT.Index=function(Index,ADD.his.nom,FOLY,NOM,Index_q.change)
{
  if(!is.null(Index))
  {
    Index$FINYEAR=as.character(Index$FINYEAR)
    if(!is.null(Index_q.change)) Index_q.change$FINYEAR=as.character(Index_q.change$FINYEAR)
    
    #confidence bounds
    YR=1:length(Index$FINYEAR)
    Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec.mon <- c(Index$LOW, tail(Index$UP, 1),rev(Index$UP), Index$LOW[1])      
    
    
    if(!is.null(Index_q.change)) FInYEAR=c(Index$FINYEAR,Index_q.change$FINYEAR)else
      FInYEAR=Index$FINYEAR
    
    N=length(FInYEAR)
    
    MxY=max(Index$UP,na.rm=T)
    
    if(ADD.his.nom=="YES")
    {
      x=which(FOLY$FINYEAR%in%FInYEAR)
      x.nom=which(NOM$FINYEAR%in%FInYEAR)
      MxY=max(c(Index$UP,FOLY$cpue[x],NOM$cpue[x]),na.rm=T)
      if(!is.null(Index_q.change)) MxY=max(c(MxY,Index_q.change$UP),na.rm=T)
    }
    
    AX=pretty(seq(0,ceiling(MxY),by=1))
    plot(1:length(Index$MEAN),Index$MEAN,ylab="",xlab="",xaxt="n",type="l",lwd=2,cex.axis=1.2,
         col=1,ylim=c(0,MxY),yaxt='n',xlim=c(1,N))   
    polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5), border = 'transparent',lwd=1.25)
    lines(1:length(Index$MEAN),Index$MEAN,lwd=2)
    if(ADD.his.nom=="YES")
    {
      points(1:length(FOLY$FINYEAR[x]),FOLY$cpue[x],cex=1.5,pch=21,bg="grey50")
      points(1:length(NOM$FINYEAR[x.nom]),NOM$cpue[x.nom],cex=1.5,pch=21,bg="white")
    }
    
    if(!is.null(Index_q.change))
    {
      YR=match(Index_q.change$FINYEAR,FInYEAR)
      Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
      Biom.Vec.mon <- c(Index_q.change$LOW, tail(Index_q.change$UP, 1),rev(Index_q.change$UP), Index_q.change$LOW[1])      
      polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5), border = 'transparent',lwd=1.25)
      lines(YR,Index_q.change$MEAN,lwd=2)
    }
    
    axis(2,at=AX,labels=T,las=2,cex.axis=1.2)
    axis(1,at=1:N,labels=F,tck=-0.025)
    if(N>=10) AX=seq(1,N,by=5)
    if(N<10) AX=seq(1,N,by=2)
    axis(1,at=AX,labels=FInYEAR[AX],tck=-0.05,cex.axis=1.1)
    
    if(!is.null(Index_q.change))
    {
      X1=match(Index$FINYEAR,FInYEAR)[length(Index$FINYEAR)]
      X2=match(Index_q.change$FINYEAR,FInYEAR)[1]
      XN=c(match(Index$FINYEAR,FInYEAR)[length(Index$FINYEAR)],match(Index_q.change$FINYEAR,FInYEAR)[1])
      polygon(x=c(XN,rev(XN)),y=c(-0.5,-0.5,MxY*10,MxY*10),col="white",border="white")
      abline(v=X1,lwd=2)
      abline(v=X2,lwd=2)
    }
  }
  
  if(is.null(Index)) plot(1,col='transparent',ann=F,axes=F)
}
PLOT.Index.zone=function(LISTA,LISTA2,ALL.yrs)
{
  ALL.yrs=data.frame(FINYEAR=ALL.yrs)
  MX=rep(NA,length(LISTA))
  for(l in 1:length(LISTA))
  {
    if(!is.null(LISTA[[l]]))
    {
      #add missing years
      LISTA[[l]]=merge(ALL.yrs,LISTA[[l]],by="FINYEAR",all=T)   
      MX[l]=max(LISTA[[l]]$UP,na.rm=T)
      if(!is.null(LISTA2[[l]])) MX[l]=max(c(LISTA[[l]]$UP,LISTA2[[l]]$UP),na.rm=T)
    }
  }
  
  MxY=max(MX,na.rm=T)
  FINYEAR=as.character(ALL.yrs$FINYEAR)
  N=length(FINYEAR)
  inc=c(0,0.15,0.3)
  plot(1:N,rep(MxY,N),ylab="",xlab="",xaxt="n",cex.axis=1.2,col='transparent',ylim=c(0,MxY),xlim=c(1,N+0.25))            
  for(l in 1:length(LISTA))
  {
    if(!is.null(LISTA[[l]]))
    {
      Index=LISTA[[l]]
      YR=1:length(Index$FINYEAR)
      YR=YR+inc[l]
      segments(YR,Index$MEAN,YR,Index$LOW,col=1,lwd=1.5)
      segments(YR,Index$MEAN,YR,Index$UP,col=1,lwd=1.5)
      points(YR,Index$MEAN,bg=CL[match(names(LISTA)[l],names(CL))],pch=21,cex=1.5,col=1)        
    }      
  }
  
  if(!is.null(LISTA2))
  {
    for(l in 1:length(LISTA2))
    {
      if(!is.null(LISTA2[[l]]))
      {
        Index=LISTA2[[l]]
        YR=match(Index$FINYEAR,ALL.yrs$FINYEAR)
        YR=YR+inc[l]
        segments(YR,Index$MEAN,YR,Index$LOW,col=1,lwd=1.5)
        segments(YR,Index$MEAN,YR,Index$UP,col=1,lwd=1.5)
        points(YR,Index$MEAN,bg=CL[match(names(LISTA2)[l],names(CL))],pch=21,cex=1.5,col=1)        
      }      
    }
  }
  
  axis(1,at=1:N,labels=F,tck=-0.025)
  if(N>=10) AX=seq(1,N,by=5)
  if(N<10) AX=seq(1,N,by=2)
  axis(1,at=AX,labels=FINYEAR[AX],tck=-0.05,cex.axis=1.1)    
  
  if(!is.null(LISTA2))
  {
    for(l in 1:length(LISTA2))
    {
      if(!is.null(LISTA2[[l]]))
      {
        Index=LISTA[[l]]
        Index=subset(Index,!is.na(MEAN))
        Index_q.change=LISTA2[[l]]
        
        X1=match(Index$FINYEAR,FINYEAR)[length(Index$FINYEAR)]
        X2=match(Index_q.change$FINYEAR,FINYEAR)[1]
        XN=c(X1,X2)
        polygon(x=c(XN,rev(XN)),y=c(-0.5,-0.5,MxY*10,MxY*10),col="white",border="white")
        abline(v=X1,lwd=2)
        abline(v=X2,lwd=2)
      }      
    }
  }
  
}
PLOT.Index.same.plot=function(Ind.mon,Ind.day)
{
  if(!is.null(Ind.mon))
  {
    Ind.mon$FINYEAR=as.character(Ind.mon$FINYEAR)
    Ind.day$FINYEAR=as.character(Ind.day$FINYEAR)
    
    Dummy.mon=Ind.mon
    Dummy.mon[,-match("FINYEAR",names(Dummy.mon))]=NA
    Dummy.day=Ind.day
    Dummy.day[,-match("FINYEAR",names(Dummy.day))]=NA
    
    #confidence bounds
    YR=1:length(Ind.mon$FINYEAR)
    Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec.mon <- c(Ind.mon$asymp.LCL, tail(Ind.mon$asymp.UCL, 1), rev(Ind.mon$asymp.UCL), Ind.mon$asymp.LCL[1])
    
    YR=YR[length(YR)]+(1:length(Ind.day$FINYEAR))
    Year.Vec.day <- c(YR, tail(YR, 1), rev(YR), YR[1])   
    Biom.Vec.day <- c(Ind.day$asymp.LCL, tail(Ind.day$asymp.UCL, 1), rev(Ind.day$asymp.UCL), Ind.day$asymp.LCL[1])
    
    Ind.mon=rbind(Ind.mon,Dummy.day)
    Ind.day=rbind(Dummy.mon,Ind.day)
    
    FInYEAR=Ind.mon$FINYEAR
    N=length(FInYEAR)
    
    #Monthly
    MxY=max(Ind.mon$asymp.UCL,na.rm=T)
    plot(1:N,Ind.mon$response,ylab="",xlab="",xaxt="n",type="l",cex=1.25,lwd=2,cex.axis=1.2,
         col=1,ylim=c(0,MxY),yaxt='n')
    polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.7), border = 'transparent',lwd=1.25)
    AX=pretty(seq(0,ceiling(MxY),by=1))
    axis(2,at=AX,labels=T,las=2,cex.axis=1.2)
    axis(1,at=1:N,labels=F,tck=-0.015)
    axis(1,at=seq(1,N,5),labels=FInYEAR[seq(1,N,5)],tck=-0.0225,cex.axis=1.1)
    polygon(x=c(-10,match("2005-06",FInYEAR),match("2005-06",FInYEAR),-10),
            y=c(-1000,-1000,MxY*1.5,MxY*1.5),border="transparent",col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.3))
    
    #Daily
    MxY=max(Ind.day$asymp.UCL,na.rm=T)
    par(new=T)
    plot(1:N,Ind.day$response,ylab="",xlab="",xaxt="n",type="l",cex=1.25,lwd=2,cex.axis=1.2,col=1,ylim=c(0,MxY),
         yaxt='n',xaxt='n')
    polygon(Year.Vec.day, Biom.Vec.day, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.7), border = 'transparent',lwd=1.25)
    AX=pretty(seq(0,ceiling(MxY),by=1))
    axis(4,at=AX,labels=T,las=2,cex.axis=1.2)
    
    polygon(x=c(match("2005-06",FInYEAR),N*1.5,N*1.5,match("2005-06",FInYEAR)),
            y=c(-1000,-1000,MxY*1.5,MxY*1.5),border="transparent",col=rgb(red=0.1, green=0.1, blue=0.1, alpha=0.1))    
    
  }
  
  if(is.null(Ind.mon)) plot(1,col='transparent',ann=F,axes=F)
}

PLOT.km.gn._vs_km.gn.h=function(Index,HOURS)
{
  if(!is.null(Index))
  {
    Index$FINYEAR=as.character(Index$FINYEAR)
    
    #confidence bounds
    YR=1:length(Index$FINYEAR)
    Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec.mon <- c(Index$LOW, tail(Index$UP, 1),rev(Index$UP), Index$LOW[1])
    Biom.Vec.mon.hour <- c(HOURS$LOW, tail(HOURS$UP, 1),rev(HOURS$UP), HOURS$LOW[1])
    FInYEAR=Index$FINYEAR
    N=length(FInYEAR)
    MxY=max(c(Index$UP,HOURS$UP),na.rm=T)
    AX=pretty(seq(0,ceiling(MxY),by=1))
    plot(1:length(Index$MEAN),Index$MEAN,ylab="",xlab="",xaxt="n",type="l",lwd=2,cex.axis=1.2,
         col=1,ylim=c(0,MxY),yaxt='n',xlim=c(1,N))   
    polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.2), border = 'transparent',lwd=1.25)
    lines(1:length(Index$MEAN),Index$MEAN,lwd=2)
    
    polygon(Year.Vec.mon, Biom.Vec.mon.hour, col = rgb(red=0.1, green=0.5, blue=0.6, alpha=0.2), border = 'transparent',lwd=1.25)
    lines(1:length(HOURS$MEAN),HOURS$MEAN,lwd=2,col=rgb(red=0.1, green=0.5, blue=0.6))
    axis(2,at=AX,labels=T,las=2,cex.axis=1.2)
    axis(1,at=1:N,labels=F,tck=-0.025)
    if(N>=10) AX=seq(1,N,by=5)
    if(N<10) AX=seq(1,N,by=2)
    axis(1,at=AX,labels=FInYEAR[AX],tck=-0.05,cex.axis=1.1)
  }
  
  if(is.null(Index)) plot(1,col='transparent',ann=F,axes=F)
}

PLOT.CV.zone=function(LISTA,LISTA2,ALL.yrs)
{
  ALL.yrs=data.frame(FINYEAR=ALL.yrs)
  MX=rep(NA,length(LISTA))
  for(l in 1:length(LISTA))
  {
    if(!is.null(LISTA[[l]]))
    {
      #add missing years
      LISTA[[l]]=merge(ALL.yrs,LISTA[[l]],by.x="FINYEAR",by.y="FINYEAR",all=T)   
      MX[l]=max(LISTA[[l]]$CV,na.rm=T)
      if(!is.null(LISTA2[[l]])) MX[l]=max(c(LISTA[[l]]$CV,LISTA2[[l]]$CV),na.rm=T)
    }
  }
  
  MxY=max(MX,na.rm=T)
  FINYEAR=as.character(ALL.yrs$FINYEAR)
  N=length(FINYEAR)
  inc=c(0,0.15,0.3)
  plot(1:N,rep(MxY,N),ylab="",xlab="",xaxt="n",cex.axis=1.2,col='transparent',ylim=c(0,MxY),xlim=c(1,N+0.25))            
  for(l in 1:length(LISTA))
  {
    if(!is.null(LISTA[[l]]))
    {
      Index=LISTA[[l]]
      YR=1:length(Index$FINYEAR)
      YR=YR+inc[l]
      points(YR,Index$CV,bg=CL[match(names(LISTA)[l],names(CL))],pch=21,cex=1.5,col=1)        
    }      
  }
  
  if(!is.null(LISTA2))
  {
    for(l in 1:length(LISTA2))
    {
      if(!is.null(LISTA2[[l]]))
      {
        Index=LISTA2[[l]]
        YR=match(Index$FINYEAR,ALL.yrs$FINYEAR)
        YR=YR+inc[l]
        points(YR,Index$CV,bg=CL[match(names(LISTA2)[l],names(CL))],pch=21,cex=1.5,col=1)        
      }      
    }
  }
  
  axis(1,at=1:N,labels=F,tck=-0.025)
  if(N>=10) AX=seq(1,N,by=5)
  if(N<10) AX=seq(1,N,by=2)
  axis(1,at=AX,labels=FINYEAR[AX],tck=-0.05,cex.axis=1.1)    
  
  if(!is.null(LISTA2))
  {
    for(l in 1:length(LISTA2))
    {
      if(!is.null(LISTA2[[l]]))
      {
        Index=LISTA[[l]]
        Index=subset(Index,!is.na(MEAN))
        Index_q.change=LISTA2[[l]]
        
        X1=match(Index$FINYEAR,FINYEAR)[length(Index$FINYEAR)]
        X2=match(Index_q.change$FINYEAR,FINYEAR)[1]
        XN=c(X1,X2)
        polygon(x=c(XN,rev(XN)),y=c(-0.5,-0.5,MxY*10,MxY*10),col="white",border="white")
        abline(v=X1,lwd=2)
        abline(v=X2,lwd=2)
      }      
    }
  }
}

fn.show.sens=function(FOLY,STAND,STAND.d,COMBO)
{
  Mn.Yrs=as.character(STAND$"Stand_Base case"$FINYEAR)
  Nn.mon=1:length(Mn.Yrs)
  Day.Yrs=as.character(STAND.d$"Stand_Base case"$FINYEAR)
  Nn.day=1:length(Day.Yrs)
  FOLY=subset(FOLY,FINYEAR%in%c(Mn.Yrs,Day.Yrs))
  
  #Monthly
  #non standardised cpue
  R=length(STAND)
  store=1:R
  for(r in 1:R) store[r]=max(STAND[[r]]$cpue)
  Ylim=max(c(FOLY$cpue[Nn.mon],store))      
  plot(Nn.mon,FOLY$cpue[Nn.mon],ylab="",xlab="",type='l',ylim=c(0,Ylim),lty=LIN[1],
       lwd=2,xaxt='n',col=CL[1],cex.axis=1.15)
  
  #standardised cpue
  x=match(names(STAND),names(CL))
  for(r in 1:R) lines(match(STAND[[r]]$FINYEAR,Mn.Yrs),STAND[[r]]$cpue,lty=LIN[x[r]],col=CL[x[r]],lwd=2.5)
  axis(1,Nn.mon,F,tck=-0.02)
  axis(1,Nn.mon[seq(1,length(Nn.mon),5)],Mn.Yrs[seq(1,length(Mn.Yrs),5)],tck=-0.04,cex.axis=1.1)
  if(COMBO=="NO")legend("bottomleft","Monthly returns",bty='n',cex=1.25)
  if(!is.na(WHERE))
  {
    legend(WHERE,nm.SCENARIOS,col=CL,lty=LIN,bty='n',cex=1.25,lwd=2.5)
    if(COMBO=="YES") mtext("Monthly returns",3,0,cex=1.5)
  }
  
  
  #Daily
  ii=length(Nn.mon)+Nn.day
  
  #non standardised cpue
  R=length(STAND.d)
  store=1:R
  for(r in 1:R) store[r]=max(STAND.d[[r]]$cpue)
  Ylim=max(c(FOLY$cpue[ii],store))      
  plot(Nn.day,FOLY$cpue[ii],ylab="",xlab="",type='l',ylim=c(0,Ylim),lty=LIN[1],
       lwd=2,xaxt='n',col=CL[1],cex.axis=1.15)
  
  #standardised cpue
  x=match(names(STAND.d),names(CL))
  for(r in 1:R) lines(match(STAND.d[[r]]$FINYEAR,Day.Yrs),STAND.d[[r]]$cpue,lty=LIN[x[r]],col=CL[x[r]],lwd=2.5)
  axis(1,Nn.day,F,tck=-0.02)
  axis(1,Nn.day[seq(1,length(Day.Yrs),2)],Day.Yrs[seq(1,length(Day.Yrs),2)],tck=-0.04,cex.axis=1.1)
  if(COMBO=="NO")legend("bottomleft","Daily logbooks",bty='n',cex=1.25)
  if(!is.na(WHERE))if(COMBO=="YES") mtext("Daily logbooks",3,0,cex=1.5)
}

Compare.Eff.creep.No.creep=function(Index,Index_no.creep,Index_q.change,Index_q.change_no.creep,ADD.CI)
{
  Index$FINYEAR=as.character(Index$FINYEAR)
  if(!is.null(Index_q.change)) Index_q.change$FINYEAR=as.character(Index_q.change$FINYEAR)
  
  #confidence bounds
  YR=1:length(Index$FINYEAR)
  Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec.mon <- c(Index$LOW, tail(Index$UP, 1),rev(Index$UP), Index$LOW[1])      
  Biom.Vec.mon_creep<- c(Index_no.creep$LOW, tail(Index_no.creep$UP, 1),rev(Index_no.creep$UP), Index_no.creep$LOW[1])  
  
  if(!is.null(Index_q.change)) FInYEAR=c(Index$FINYEAR,Index_q.change$FINYEAR)else
    FInYEAR=Index$FINYEAR
  
  N=length(FInYEAR)
  
  MxY=max(Index_no.creep$MEAN,na.rm=T)
  if(!is.null(Index_q.change)) MxY=max(c(MxY,Index_q.change_no.creep$MEAN),na.rm=T)
  
  if(ADD.CI=="YES") 
  {
    MxY=max(Index_no.creep$UP,na.rm=T)
    if(!is.null(Index_q.change)) MxY=max(c(MxY,Index_q.change_no.creep$UP),na.rm=T)
  }
  
  
  
  AX=pretty(seq(0,ceiling(MxY),by=1))
  plot(1:length(Index$MEAN),Index$MEAN,ylab="",xlab="",xaxt="n",type="l",lwd=2,cex.axis=1.2,
       col=1,ylim=c(0,MxY),yaxt='n',xlim=c(1,N))   
  
  if(ADD.CI=="YES")polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5), border = 'transparent',lwd=1.25)
  
  lines(1:length(Index_no.creep$MEAN),Index_no.creep$MEAN,lwd=2,lty=2)
  if(ADD.CI=="YES")polygon(Year.Vec.mon, Biom.Vec.mon_creep, col = rgb(red=0.25, green=0.25, blue=0.25, alpha=0.3), border = 'transparent',lwd=1.25)
  
  if(!is.null(Index_q.change))
  {
    YR=match(Index_q.change$FINYEAR,FInYEAR)
    Year.Vec.mon <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec.mon <- c(Index_q.change$LOW, tail(Index_q.change$UP, 1),rev(Index_q.change$UP), Index_q.change$LOW[1])      
    Biom.Vec.mon_creep <- c(Index_q.change_no.creep$LOW, tail(Index_q.change_no.creep$UP, 1),rev(Index_q.change_no.creep$UP), Index_q.change_no.creep$LOW[1])      
    
    lines(YR,Index_q.change$MEAN,lwd=2)
    if(ADD.CI=="YES")polygon(Year.Vec.mon, Biom.Vec.mon, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5), border = 'transparent',lwd=1.25)
    
    lines(YR,Index_q.change_no.creep$MEAN,lwd=2,lty=2)
    if(ADD.CI=="YES")polygon(Year.Vec.mon, Biom.Vec.mon_creep, col = rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5), border = 'transparent',lwd=1.25)
    
  }
  
  
  axis(2,at=AX,labels=T,las=2,cex.axis=1.2)
  axis(1,at=1:N,labels=F,tck=-0.025)
  if(N>=10) AX=seq(1,N,by=5)
  if(N<10) AX=seq(1,N,by=2)
  axis(1,at=AX,labels=FInYEAR[AX],tck=-0.05,cex.axis=1.1)
  
  if(!is.null(Index_q.change))
  {
    X1=match(Index$FINYEAR,FInYEAR)[length(Index$FINYEAR)]
    X2=match(Index_q.change$FINYEAR,FInYEAR)[1]
    XN=c(match(Index$FINYEAR,FInYEAR)[length(Index$FINYEAR)],match(Index_q.change$FINYEAR,FInYEAR)[1])
    polygon(x=c(XN,rev(XN)),y=c(-0.5,-0.5,MxY*10,MxY*10),col="white",border="white")
    abline(v=X1,lwd=2)
    abline(v=X2,lwd=2)
  }
}


#functions for influence plots
bubble.plot=function(x,y,z,scale,Xlab,Ylab)  
{
  xo=outer(x,rep(1,length=length(y)))
  yo=t(outer(y,rep(1,length=length(x))))
  zo=z*scale
  matplot(xo,yo,type="n",xlab=Xlab,ylab=Ylab,xaxt='n',yaxt='n')
  for(s in 1:length(x))
  {
    points(xo[s,],yo[s,],cex=zo[,s],pch=16,col="grey80")
    points(xo[s,],yo[s,],cex=zo[,s],pch=1,col="black")
  }
}
Influence.fn=function(MOD,DAT,Term.type,termS,LABEL,XLIM,SCALE,add.Influence)  #Bentley et al 2012
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
  ny=table(DAT$FINYEAR)
  
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
      a=subset(DAT,FINYEAR==names(ny[t]))
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
      
      fn.fig(paste(HnDll,SPECIES.vec[i],".CDI.",termS[p],sep=""),2400, 1600)
      nf <- layout(matrix(c(1,1,1,2,2,2), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
      if(add.Influence=="YES")nf <- layout(matrix(c(1,1,0,2,2,3), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
      par(mar=c(0,0,0,0),oma=c(4,6,1,1),las=1,mgp=c(1,.9,0))
      layout.show(nf)
      
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
      Prop.Rec=NULL
      for(t in 1:length(ny))
      {
        a=subset(DAT,FINYEAR==names(ny[t]))
        Prop.rec=table(a[,match(termS[p],names(a))])
        Prop.rec=Prop.rec[-1]/sum(Prop.rec[-1])
        Prop.Rec=rbind(Prop.Rec,Prop.rec)
      }  
      if(colnames(Prop.Rec)[1]=="1")Prop.Rec=Prop.Rec[,-1]
      Nombres=gsub("[^[:digit:]]", "", names(COEF))
      if(termS[p]=="VESSEL") Nombres=1:length(Nombres)   #change vessel name for dummy
      rownames(Prop.Rec)=names(ny[1:length(ny)])
      colnames(Prop.Rec)=COef.nm
      Prop.Rec=Prop.Rec[,match(COef.nm.sorted,colnames(Prop.Rec))]
      bubble.plot(x,1:length(ny),Prop.Rec,SCALE[p],termS[p],"Financial year")
      axis(1,1:length(COEF),F,tck=-0.015)
      axis(1,seq(1,length(COEF),1),Nombres[seq(1,length(COEF),1)],cex.axis=1.15,tck=-0.025)
      mtext(LABEL[p],side=1,line=2.5,cex=1.5)
      axis(2,1:length(ny),F,tck=-0.015)
      axis(2,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1,tck=-0.025)   
      mtext("Financial year",side=2,line=4,cex=1.5,las=3)
      
      if(add.Influence=="YES")
      {
        #Influence plot
        plot(Annual.Dev[[p]],1:length(ny),type="o",pch=19,xlab="",ylab="",cex=2,cex.axis=1.25,yaxt='n',xlim=XLIM)
        abline(v=1,lty=3,col=1)
        mtext("Influence",side=1,line=2.5,cex=1.5)      
        axis(2,1:length(ny),F,tcl=0.5)
        axis(2,seq(1,length(ny),2),F,tcl=1)
        return(list(Annual.Dev=Annual.Dev,ny=ny,Over.all.influence=Over.all.influence))
      }
      
      
      dev.off()
    }
  }
  
  
}

Fig3.plot.fun=function(A,YLIM,WHERE)
{
  Annual.Dev=A$Annual.Dev
  ny=A$ny
  
  NamE=names(A$Over.all.influence)
  NamE=ifelse(NamE=="BLOCKX","Block",
              ifelse(NamE=="VESSEL","Vessel",
                     ifelse(NamE=="MONTH","Month",
                            ifelse(NamE=="log.Catch.Whiskery","Whiskery_c",
                                   ifelse(NamE=="log.Catch.Whiskery","Whiskery_c",
                                          ifelse(NamE=="log.Catch.Gummy","Gummy_c",
                                                 ifelse(NamE=="log.Catch.Dusky","Dusky_c",
                                                        ifelse(NamE=="log.Catch.Sandbar","Sandbar_c",
                                                               ifelse(NamE=="log.Catch.Scalefish","Scalefish_c",NamE)))))))))
  
  nt=length(Annual.Dev)  
  LTY=c(1,4,3,1,3,2)
  plot(1:length(ny),Annual.Dev[[1]],col=LTY.col[1],type="l",xlab="",ylab="",lwd=LWD,
       cex.axis=1.35,xaxt='n',ylim=YLIM)
  abline(h=1,lty=3,col=1)
  axis(1,1:length(ny),F,tck=-0.02)
  axis(1,seq(1,length(ny),2),F,tck=-0.04)
  axis(1,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1.35,tck=-0.04)
  for(p in 2:nt)lines(1:length(ny),Annual.Dev[[p]],lwd=LWD,lty=LTY[p],col=LTY.col[p])
  LEG=paste(NamE," (",100*round(A$Over.all.influence,2),"%)",sep="")
  legend(WHERE,LEG,bty='n',lty=LTY,col=LTY.col,lwd=LWD,cex=1.25,pt.cex=1.5)
}


#function to explore vessels with catch of study species and blocks fished to aim for balanced design
# fn.see.all.yrs.ves.blks.agg=function(a,b,SP,Ves.sel.BC,Ves.sel.sens,BLK.sel.BC,BLK.sel.sens)
# {
#   
#   a=subset(a,select=c(MONTH,YEAR.c,BLOCKX,VESSEL,Same.return,SPECIES,LIVEWT.c,Km.Gillnet.Days.c,Km.Gillnet.Hours.c))
#   b=subset(b,select=c(MONTH,YEAR.c,BLOCKX,VESSEL,Same.return,SPECIES,LIVEWT.c,Km.Gillnet.Days.c,Km.Gillnet.Hours.c))
#   
#   All=rbind(a,b)
#   CATCH.sp=aggregate(LIVEWT.c~YEAR.c+VESSEL,subset(All,SPECIES==SP),sum)
#   CATCH.sp=reshape(CATCH.sp,v.names="LIVEWT.c",idvar="YEAR.c",timevar="VESSEL",direction="wide")
#   CATCH.sp=CATCH.sp[order(CATCH.sp$YEAR.c),]
#   
#   Vess=substr(names(CATCH.sp)[2:ncol(CATCH.sp)],10,30)
#   Yrs=CATCH.sp$YEAR.c
#   
#   #All vess total catch
#   Z=as.matrix(CATCH.sp[,-1])
#   #image(x=1:length(Yrs),y=1:length(Vess),Z,xaxt='n',yaxt='n',ann=F)
#   #axis(1,1:length(Yrs),Yrs)
#   #axis(2,1:length(Vess),F)
#   #axis(2,seq(1,length(Vess),5),Vess[seq(1,length(Vess),5)],las=1,cex.axis=.6)
#   
#   
#   #Vess with > X years of records
#   ZZ=Z
#   ZZ[ZZ>0]=1
#   Yrs.with.ktch=colSums(ZZ,na.rm=T)
#   
#   pdf(paste("C:/Matias/Analyses/Catch and effort/Outputs/Vessel_pos_records_by_yr/",SP,".pdf",sep="")) 
#   
#   #Ves.sel.BC
#   par(mar=c(3,3.5,.8,.8))
#   WHICh=which(Yrs.with.ktch>Ves.sel.BC)
#   Z.this=ZZ[,WHICh]
#   Ves.BC=Vess[WHICh]
#   image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:ncol(Z.this),Ves.BC,las=1,cex.axis=.9)
#   legend("top",paste("vessels with >=",Ves.sel.BC, "years of records of",SP),bty='n')
#   
#   #Ves.sel.sens
#   WHICh=which(Yrs.with.ktch>Ves.sel.sens)
#   Z.this=ZZ[,WHICh]
#   Ves.Sens=Vess[WHICh]
#   image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:ncol(Z.this),Ves.Sens,las=1,cex.axis=.6)
#   legend("top",paste("vessels with >=",Ves.sel.sens, "years of records of",SP),bty='n')
#   
#   #plot CPUEs
#   #monthly
#   a$CPUE.km.gn.day=a$LIVEWT.c/a$Km.Gillnet.Days.c
#   a$CPUE.km.gn.h=a$LIVEWT.c/a$Km.Gillnet.Hours.c
#   
#   a.mean.cpue.km.day_all=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,SPECIES==SP),mean)
#   a.mean.cpue.km.h_all=aggregate(CPUE.km.gn.h~YEAR.c,subset(a, SPECIES==SP),mean)
#   
#   a.mean.cpue.km.day_Sens=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
#   a.mean.cpue.km.h_Sens=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
#   
#   a.mean.cpue.km.day_BC=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES==SP),mean)
#   a.mean.cpue.km.h_BC=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES==SP),mean)
#   
#   #kmgday
#   par(mar=c(3,3,.5,4),mgp=c(2,.7,0))
#   Yrs=a.mean.cpue.km.day_Sens$YEAR.c
#   plot(Yrs,a.mean.cpue.km.day_all$CPUE.km.gn.day,ylab="",xlab="")
#   points(Yrs,a.mean.cpue.km.day_Sens$CPUE.km.gn.day,pch=19,col=2)
#   if(nrow(a.mean.cpue.km.day_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.day_BC$CPUE.km.gn.day,pch=19,col=3)
#   
#   #kmgday  
#   par(new = T)
#   plot(Yrs,a.mean.cpue.km.h_all$CPUE.km.gn.h,ylab=NA, axes=F,xlab=NA,pch=0,cex=2)
#   points(Yrs,a.mean.cpue.km.h_Sens$CPUE.km.gn.h,pch=15,col=2,cex=2)
#   if(nrow(a.mean.cpue.km.h_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.h_BC$CPUE.km.gn.h,pch=15,col=3,cex=2)
#   axis(side = 4)
#   mtext("Financial year",1,line=2)
#   mtext("Nominal CPUE (Kg/km.gn.day)",2,line=2)
#   mtext("Nominal CPUE (Kg/km.gn.hour)",4,line=2)
#   legend("top",c("Kg/km.gn.day","Kg/km.gn.hour"),bty='n',pch=c(0,19))
#   legend("bottomleft",c("all",paste(Ves.sel.sens,"y"),paste(Ves.sel.BC,"y")),bty='n',pch=c(0,19,19),col=c(1,2,3))
#   
#   
#   
#   #daily
#   b$CPUE.km.gn.day=b$LIVEWT.c/b$Km.Gillnet.Days.c
#   b$CPUE.km.gn.h=b$LIVEWT.c/b$Km.Gillnet.Hours.c
#   
#   b.mean.cpue.km.day_all=aggregate(CPUE.km.gn.day~YEAR.c,subset(b,SPECIES==SP),mean)
#   b.mean.cpue.km.h_all=aggregate(CPUE.km.gn.h~YEAR.c,subset(b, SPECIES==SP),mean)
#   
#   b.mean.cpue.km.day_Sens=aggregate(CPUE.km.gn.day~YEAR.c,subset(b,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
#   b.mean.cpue.km.h_Sens=aggregate(CPUE.km.gn.h~YEAR.c,subset(b,VESSEL%in%Ves.Sens& SPECIES==SP),mean)
#   
#   b.mean.cpue.km.day_BC=aggregate(CPUE.km.gn.day~YEAR.c,subset(b,VESSEL%in%Ves.BC& SPECIES==SP),mean)
#   b.mean.cpue.km.h_BC=aggregate(CPUE.km.gn.h~YEAR.c,subset(b,VESSEL%in%Ves.BC& SPECIES==SP),mean)
#   
#   #kmgday
#   Yrs=b.mean.cpue.km.day_Sens$YEAR.c
#   plot(Yrs,b.mean.cpue.km.day_all$CPUE.km.gn.day,ylab="",xlab="")
#   points(Yrs,b.mean.cpue.km.day_Sens$CPUE.km.gn.day,pch=19,col=2)
#   points(Yrs,b.mean.cpue.km.day_BC$CPUE.km.gn.day,pch=19,col=3)
#   
#   #kmgday  
#   par(new = T)
#   plot(Yrs,b.mean.cpue.km.h_all$CPUE.km.gn.h,ylab=NA, axes=F,xlab=NA,pch=0,cex=2)
#   points(Yrs,b.mean.cpue.km.h_Sens$CPUE.km.gn.h,pch=15,cex=2,col=2)
#   points(Yrs,b.mean.cpue.km.h_BC$CPUE.km.gn.h,pch=15,col=3,cex=2)
#   axis(side = 4)
#   mtext("Financial year",1,line=2)
#   mtext("Nominal CPUE (Kg/km.gn.day)",2,line=2)
#   mtext("Nominal CPUE (Kg/km.gn.hour)",4,line=2)
#   legend("top",c("Kg/km.gn.day","Kg/km.gn.hour"),bty='n',pch=c(0,19))
#   legend("bottomleft",c("all",paste(Ves.sel.sens,"y"),paste(Ves.sel.BC,"y")),bty='n',pch=c(0,19,19),col=c(1,2,3))
#   
#   
#   
#   #For each vessel group, plot number of blocks by year
#   #Ves.BC
#   #alll blocks
#   AA=aggregate(LIVEWT.c~YEAR.c+BLOCKX,subset(All,SPECIES==SP & VESSEL%in%Ves.BC),sum)
#   AA=reshape(AA,v.names="LIVEWT.c",idvar="YEAR.c",timevar="BLOCKX",direction="wide")
#   AA=AA[order(AA$YEAR.c),]
#   BLOCs=substr(names(AA)[2:ncol(AA)],10,30)
#   Yrs=AA$YEAR.c
#   Z=as.matrix(AA[,-1])
#   ZZ=Z
#   ZZ[ZZ>0]=1
#   par(mar=c(3,3.5,.8,.8))
#   image(x=1:length(Yrs),y=1:ncol(ZZ),ZZ,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:length(BLOCs),BLOCs,las=1,cex.axis=.9)
#   legend("top",paste("all blocks for vessels >=",Ves.sel.BC, "years of records of",SP),bty='n')
#   
#   #blocks with > bc records
#   Yrs.with.ktch=colSums(ZZ,na.rm=T)
#   WHICh=which(Yrs.with.ktch>BLK.sel.BC)
#   Z.this=ZZ[,WHICh]
#   Blks.BC=BLOCs[WHICh]
#   image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:ncol(Z.this),Blks.BC,las=1,cex.axis=.9)
#   legend("top",paste("blocks with >=",BLK.sel.BC, "of records for vessels >=",Ves.sel.BC,"years of records of",SP),bty='n')
#   
#   
#   #Ves.Sens
#   #alll blocks
#   AA=aggregate(LIVEWT.c~YEAR.c+BLOCKX,subset(All,SPECIES==SP & VESSEL%in%Ves.Sens),sum)
#   AA=reshape(AA,v.names="LIVEWT.c",idvar="YEAR.c",timevar="BLOCKX",direction="wide")
#   AA=AA[order(AA$YEAR.c),]
#   BLOCs=substr(names(AA)[2:ncol(AA)],10,30)
#   Yrs=AA$YEAR.c
#   Z=as.matrix(AA[,-1])
#   ZZ=Z
#   ZZ[ZZ>0]=1
#   par(mar=c(3,3.5,.8,.8))
#   image(x=1:length(Yrs),y=1:ncol(ZZ),ZZ,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:length(BLOCs),BLOCs,las=1,cex.axis=.9)
#   legend("top",paste("all blocks for vessels >=",Ves.sel.sens, "years of records of",SP),bty='n')
#   
#   #blocks with > Sens records
#   Yrs.with.ktch=colSums(ZZ,na.rm=T)
#   WHICh=which(Yrs.with.ktch>BLK.sel.sens)
#   Z.this=ZZ[,WHICh]
#   Blks.Sens=BLOCs[WHICh]
#   image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
#   axis(1,1:length(Yrs),Yrs)
#   axis(2,1:ncol(Z.this),Blks.Sens,las=1,cex.axis=.9)
#   legend("top",paste("blocks with >=",BLK.sel.sens, "years of records for vessels >=",Ves.sel.sens,"years of records of",SP),bty='n')
#   
#   
#   
#   dev.off()
#   
#   return(list(Ves.BC=Ves.BC, Ves.Sens=Ves.Sens, Blks.BC=Blks.BC, Blks.Sens=Blks.Sens))
# }

fn.see.all.yrs.ves.blks=function(a,SP,what,Ves.sel.BC,Ves.sel.sens,BLK.sel.BC,BLK.sel.sens,Min.ktch)
{
  
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


#export simple cpue (positive catch only)
simple=function(d,SP)
{
  print("Lat and long range")
  print(range(d$LAT))
  print(range(d$LONG))
  
  print('table of method')
  print(table(d$METHOD,useNA = 'ifany'))
  
  print('table of mesh')
  print(table(d$mesh,useNA = 'ifany'))
  
  print("effort == effort no creep")
  print(nrow(d)==sum(d$Km.Gillnet.Days.c==d$Km.Gillnet.Days.c.no.creep))
  
  dd=subset(d,SPECIES%in%SP)
  dd$cpue=dd$LIVEWT.c/dd$Km.Gillnet.Days.c
  dd$cpue.hours=dd$LIVEWT.c/dd$Km.Gillnet.Hours.c
  
  CPUE.nominal=aggregate(cpue~FINYEAR,dd,mean)
  CPUE.nominal.hours=aggregate(cpue.hours~FINYEAR,dd,mean)
  CATCH=aggregate(LIVEWT.c~FINYEAR,dd,sum)
  EFRT=aggregate(Km.Gillnet.Hours.c~FINYEAR,dd,sum)
  CPUE.folly=0.38*CATCH[,2]/EFRT[,2]
  CPUE.nominal$cpue=0.38*CPUE.nominal$cpue
  CPUE.nominal.hours$cpue.hours=0.38*CPUE.nominal.hours$cpue.hours
  y=CPUE.nominal$cpue/mean(CPUE.nominal$cpue)
  y1=CPUE.nominal.hours$cpue.hours/mean(CPUE.nominal.hours$cpue.hours)
  y2=CPUE.folly/mean(CPUE.folly)
  plot(y,type='l',xaxt='n',main=unique(dd$SNAME)[1],lwd=2,ylim=c(min(c(y,y1,y2)),max(c(y,y1,y2))))
  axis(1,1:nrow(CPUE.nominal),CPUE.nominal$FINYEAR)
  lines(y1,col=2,lwd=2)
  lines(y2,col=3,lwd=2)
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

#calculate 4 different nominal cpues
fn.ainslie=function(dat,Ktch.targt,Effrt,explr,QL_prop_ktch,Prop.Malcolm,cpue.units,spname,BLks,VesL,Show.QL)
{
  source("C:\\Matias\\Analyses\\SOURCE_SCRIPTS\\Population dynamics\\Nominal_cpue_functions.R")
  
  names(dat) =  casefold(names(dat))
  
  ## some variable names
  dat$catch = dat[,match(Ktch.targt,names(dat))]
  dat$effort = dat[,match(Effrt,names(dat))]
  dat$year = dat$year.c
  dat$fymonth = factor(dat$month, levels=c(7:12, 1:6))
  dat$season = as.numeric(substring(dat$finyear, 1, 4))
  dat$smonth = factor(dat$month, levels=c(7:12, 1:6))
  
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
  
  if(Show.QL)
  {
    par(mfrow=c(1,1))
    boxplot(prop ~ season, dat)
    mtext("Proportion of target species catch out of total catch",3,-1,col="red")
  } 
  qldat = CalcQL(dat, prop.catch=QL_prop_ktch)
  dat = merge(dat, qldat, all=TRUE)
  dat$target = ifelse(dat$prop > dat$ql, 1, 0)
  if(Show.QL)
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
  
  if(Show.QL)
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
    mtext(paste("Cpue (",cpue.units,")",sep=""),2,outer=T,las=3)
    mtext("Proportion",1,outer=T)
  }
  

  ##==========================================================================
  ## plot raw mean cpues and CIs using 4 different data sets
  ##==========================================================================
  
  #compare different subsets of the data
  par(mfrow=c(3,2),mar=c(3,3,1,1),oma=c(1,1,1,1),mgp=c(2,.5,0))
  
  CPUE.All = CalcMeanCPUE(cpuedata = dat, catch.column="catch", effort.column="effort",
                    plot.title = paste(spname, "_All rec"), cpue.units = cpue.units, 
                    draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
  
  CPUE.blk_vess = CalcMeanCPUE(cpuedata = subset(dat, blockx%in%BLks & vessel%in%VesL), catch.column="catch", 
                    effort.column="effort",plot.title = paste(spname, "_indicative_blk_ves"),
                    cpue.units = cpue.units,draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")

  CPUE.non_zero = CalcMeanCPUE(cpuedata = subset(dat, catch>0), catch.column="catch", 
                    effort.column="effort",plot.title = paste(spname, "_Nonzero"),
                    cpue.units = cpue.units,draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
  
  CPUE.QL_target = CalcMeanCPUE(cpuedata = subset(dat, target==1), catch.column="catch",
                    effort.column="effort",plot.title = paste(spname, "_QL Target"),
                    cpue.units = cpue.units,draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
  
  CPUE.Malcolm = CalcMeanCPUE(cpuedata = subset(dat, prop>Prop.Malcolm), catch.column="catch",
                  effort.column="effort",plot.title = paste(spname, "_Malcolm_Prop>",Prop.Malcolm),
                  cpue.units = cpue.units,draw.plot=TRUE, show.legend=TRUE,PaR="NO",showLNMean="YES")
  
  return(list(CPUE.All=CPUE.All,CPUE.blk_vess=CPUE.blk_vess,CPUE.non_zero=CPUE.non_zero,
              CPUE.QL_target=CPUE.QL_target,CPUE.Malcolm=CPUE.Malcolm,QL_dat=subset(dat, target==1)))
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
  
  
  source("C:\\Matias\\Analyses\\SOURCE_SCRIPTS\\Population dynamics\\Nominal_cpue_functions.R")
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

Get.Mns=function(d,grp,Vars,LGND,add.arrow)
{
  d=d[,match(c(grp,Vars),names(d))]
  d$cpue.d=d$Catch.Target/d$Km.Gillnet.Days.c
  d$cpue.h=d$Catch.Target/d$Km.Gillnet.Hours.c
  for(v in 1:length(Vars))
  {
    B= d[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      as.data.frame
    B$yr=as.numeric(substr(B$FINYEAR,1,4))
    plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col='blue',ylim=c(min(B$lowCL),max(B$uppCL)))
    arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='blue')
    if(add.arrow[v])
    {
      Is=(nrow(B)-5):nrow(B)
      arrows(B$yr[1],B$mean[1],mean(B$yr[Is]),mean(B$mean[Is]),col='black',lwd=2)
      legend("bottomright",paste(round(mean(B$mean[Is])/B$mean[1],1),"fold",sep="-"),bty='n',cex=1.5)
      with(B[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.8,.1,.1,.3)))
    }
  }
  
  
  B= d[,match(c(grp,"cpue.d"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr,B$mean,xlab="",ylab="kg/km gillnet days",pch=19,col='blue',ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='blue')
  
  par(new=T)
  B= d[,match(c(grp,"cpue.h"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr,B$mean,xlab="",ylab="",pch=19,col='red',axes=F,ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='red')
  axis(side = 4,col="red",col.axis = "red")
  mtext("kg/km gillnet hours",4,line=1.5,col="red",cex=1)
  
}




#----4. PROCEDURE SECTION-----#

#export simple cpue (positive catch only)
if(Model.run=="First") 
{
  fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/Simple.cpue",2400,2400)
  par(mfcol=c(2,2),mar=c(1,1,3,1),oma=c(2,4,.1,.5),las=1,mgp=c(1.9,.6,0))
  simple(d=Data.daily.GN.whiskery,SP=17003)
  simple(d=Data.daily.GN.gummy,SP=17001)
  simple(d=Data.daily.GN.dusky,SP=c(18003,18001))
  simple(d=Data.daily.GN.sandbar,SP=18007)
  mtext("Relative cpue",side=2,las=3,line=1,outer=T,cex=1.5)
  legend("bottomleft",c('nominal cpue','nominal cpue hours','effective cpue'),
         col=1:3,lty=1,bty='n',lwd=2,cex=1)
  dev.off()
}



#Keep vessel characteristics from vessel survey for vessels that have fished      #missing: fill in missing vessel info and add to cpue data. Also add wheather prediction
TDGDLF.survey=subset(TDGDLF.survey, BOATREGO%in%unique(
  c(Data.monthly.GN.whiskery$VESSEL,Data.monthly.GN.gummy$VESSEL,
    Data.monthly.GN.dusky$VESSEL,Data.monthly.GN.sandbar$VESSEL,
    Data.daily.GN.whiskery$VESSEL,Data.daily.GN.gummy$VESSEL,
    Data.daily.GN.dusky$VESSEL,Data.daily.GN.sandbar$VESSEL)
  ))    


#4.2 Deal with zone1-zone2 Boundary blocks to a zone 
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


#4.3 Extract number of vessels per species range
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


#4.4 Extract number of blocks where shark has been caught within effective area
Tol.blks.whi=length(fn.n.blk.caught(Data.monthly.GN.whiskery,Data.daily.GN.whiskery,17003))
Tol.blks.gum=length(fn.n.blk.caught(Data.monthly.GN.gummy,Data.daily.GN.gummy,17001))
Tol.blks.dus=length(fn.n.blk.caught(Data.monthly.GN.dusky,Data.daily.GN.dusky,c(18003,18001)))
Tol.blks.san=length(fn.n.blk.caught(subset(Data.monthly.GN.sandbar,YEAR.c>=1986),Data.daily.GN.sandbar,18007))

Blks.by.species=cbind(Tol.blks.whi,Tol.blks.gum,Tol.blks.dus,Tol.blks.san)
colnames(Blks.by.species)=c("WH","GM","Dus","San")
write.csv(Blks.by.species,paste(hndl,"All.Blks.by.species.csv",sep=""),row.names=F)


#4.5 Block weights     
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


#4.4 Data fixes
#Remove initial years for sandbar as they were not reported in fishery stats
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,
                               !(FINYEAR%in%paste(1975:1984,substr(1976:1985,3,4),sep="-")))


#4.5 Remove NA effort
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
FINYEAR.ALL=as.character(unique(c(Data.monthly.GN.whiskery$FINYEAR,Data.daily.GN.whiskery$FINYEAR)))
N.yrs.ALL=length(FINYEAR.ALL)

FINYEAR.monthly=as.character(unique(Data.monthly.GN.whiskery$FINYEAR))
N.yrs=length(FINYEAR.monthly)

FINYEAR.daily=as.character(unique(Data.daily.GN.whiskery$FINYEAR))
N.yrs.daily=length(FINYEAR.daily) 

new.q=FINYEAR.monthly[match(Q_change,FINYEAR.monthly):length(FINYEAR.monthly)]  #whiskery q years


#Explore proportion of dusky and copper shark
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


#4.9 Define indicative vessels and blocks 
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


#4.10 Put data into a list
#note: this still has all records, all vessels and all blocks
Species.list=list(whiskery=Data.monthly.GN.whiskery,gummy=Data.monthly.GN.gummy,
                  dusky=Data.monthly.GN.dusky,sandbar=Data.monthly.GN.sandbar)

Species.list.daily=list(whiskery=Data.daily.GN.whiskery,gummy=Data.daily.GN.gummy,
                        dusky=Data.daily.GN.dusky,sandbar=Data.daily.GN.sandbar)

rm(Data.monthly.GN.whiskery,Data.monthly.GN.gummy,Data.monthly.GN.dusky,Data.monthly.GN.sandbar,
   Data.daily.GN.whiskery,Data.daily.GN.gummy,Data.daily.GN.dusky,Data.daily.GN.sandbar)


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


#4.11 Compare nominal all records VS 'good reporters' only
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


#4.12 Construct wide database for analysis 
#steps: 
#   1. select "Good" records (the variable "Reporter" includes good/bad catch and effort reporters)
#   2. Construct a single record for each record (i.e. year-month-vessel-block-gear for monthly
#      returns and year-Session-vessel-block10-gear for daily logbooks), with catch of target
#      and other species as separate columns, giving a 0 catch for column "target" if no catch
#      where a Session is TSNo or SNo.DSNo.TSNo (see below)

#note: this still keeps all blocks and vessels from good reporters

#variables used in standardisation
    #monthly
These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Days.inv","Km.Gillnet.Hours.c","Km.Gillnet.Days.c","Km.Gillnet.Days.c.no.creep",
                "zone","MONTH","BLOCKX","mesh","SHOTS.c","BDAYS.c","HOURS.c","NETLEN.c")

    #daily
These.efforts.daily=c("FINYEAR","date","TSNo","Km.Gillnet.Days.inv","Km.Gillnet.Hours.c","Km.Gillnet.Days.c","Km.Gillnet.Days.c.no.creep",
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

#Add SOI, Freo and moon phase to SNo    
  #Freo lag
Freo$Freo.Lag6=c(rep(NA,6),Freo$Freo[1:(length(Freo$Freo)-6)])
Freo$Freo.Lag12=c(rep(NA,12),Freo$Freo[1:(length(Freo$Freo)-12)])
Freo=subset(Freo,Year>1973)
for(i in 1:N.species)
{
  # #monthly
  # a=SOI
  # a$dummy=with(a,paste(Year.soi, Month.soi))
  # a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c[[i]],paste(YEAR.c, MONTH))))
  # a=a[,-match('dummy',names(a))]
  # DATA.list.LIVEWT.c[[i]]=merge(DATA.list.LIVEWT.c[[i]],a, by.x=c("YEAR.c","MONTH"),
  #                               by.y=c("Year.soi","Month.soi"),all.x=T)
  # 
  # a=Freo
  # a$dummy=with(a,paste(Year, Month))
  # a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c[[i]],paste(YEAR.c, MONTH))))
  # a=a[,-match('dummy',names(a))]
  # DATA.list.LIVEWT.c[[i]]=merge(DATA.list.LIVEWT.c[[i]],a, by.x=c("YEAR.c","MONTH"),
  #                         by.y=c("Year","Month"),all.x=T)
  
  #daily
  a=SOI
  a$dummy=with(a,paste(Year.soi, Month.soi))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a, by.x=c("YEAR.c","MONTH"),
                                by.y=c("Year.soi","Month.soi"),all.x=T)
  
  a=Freo
  a$dummy=with(a,paste(Year, Month))
  a=subset(a,dummy%in%unique(with(DATA.list.LIVEWT.c.daily[[i]],paste(YEAR.c, MONTH))))
  a=a[,-match('dummy',names(a))]
  DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a, by.x=c("YEAR.c","MONTH"),
                                by.y=c("Year","Month"),all.x=T)
  rm(a)
  
  #Add moon phase to daily records
  DATA.list.LIVEWT.c.daily_nfish[[i]]$Moon=lunar.phase(as.Date(DATA.list.LIVEWT.c.daily_nfish[[i]]$date),name=4)
}



#Explore cpue by block
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


#drop first years of sandbar data because vessels don't meet selection 
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


#4.13 Compare wide vs long (raw) data sets
if(Model.run=="First")
{
  #Monthly
  for(i in 1:length(TARGETS))
  {
    HH="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Compare cpue, catch, vessel, effort"
    fn.fig(paste(HH,"/LongVsWide.monthly_",SPECIES.vec[i],sep=""),2400, 2400)
    Wide.vs.long(subset(DATA.list.LIVEWT.c_all.blks[[i]],Catch.Target>0),
                 subset(Species.list_all.blks[[i]],SPECIES%in%TARGETS[[i]] & Reporter=="good"))
    mtext(TARGETS[[i]],3,line=-2)
    dev.off()
  }
  
  #Daily
  for(i in 1:length(TARGETS))
  {
    HH="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Compare cpue, catch, vessel, effort"
    fn.fig(paste(HH,"/LongVsWide.daily_",SPECIES.vec[i],sep=""),2400,2400)
    
    D.wide=subset(DATA.list.LIVEWT.c.daily[[i]])
    D.long=subset(Species.list.daily[[i]],SPECIES%in%TARGETS[[i]] & 
                    Reporter=="good"&TSNo%in%unique(D.wide$TSNo))
    DATA=subset(Species.list.daily[[i]],Reporter=="good")
    Wide.vs.long.daily(D.wide,D.long,DATA)
    mtext(TARGETS[[i]],3,line=-2)
    dev.off()
  }
  
}


#4.14 Corroborate effective area  
if(Model.run=="First")for(s in 1:N.species)
{
  #Daily
  
  fn.check.eff.area(ALL=DATA.list.LIVEWT.c.daily,WHAT="TSNo",
                    Effrt='km.gillnet.days.c',
                    QL_prop_ktch=.9,
                    spname=casefold(unlist(strsplit(SPECIES.vec[s], " "))[1]),
                    TYPE='Daily')
}


#4.15  Identify targeting behaviour   
props="YES"  #response variable as proportion
SQRT="NO"   #transform
CENTRE="YES"
SCALE="YES"

use.what="catch"  #use proportional catch to be independent of abundance
#use.what="cpue"

Ktch.var=c("Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar")

agg_PCA="by_Same.return.SNo"
HndL.Species_targeting="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Species_targeting/"

if(do.cluster=="YES")   
{
  # #Hoyle et al 2015
  # Targt=c(17003,17001,18003,18007)
  # for(i in 1:N.species)  fn.targeting(Species.list.daily[[i]],SP=Targt[i])  
  
  library('cluster')
  library(factoextra) #for plotting
  
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
  
  rm(ALL,df,d.cpue,ALL.shots)
}

if(do.PCA=="YES")     
{
  if(Clusterable=="clusterable")
  {
    PercentExpl=90   #percent of cumulative variance explained for Axis selection
    library(factoextra) #for PCA plotting
    library(plyr)
    library(dplyr)
    library(ggbiplot)
    
    
    #PCA analysis as per Winker et al 2014 on nfish as this has data at Sesssion level
    
    #note:PCA is done on combined data set

    
    #check correlation
    M<-cor(df)
    #M<-cor(d.cpue,use = "pairwise.complete.obs")
    fn.fig(paste(HndL.Species_targeting,"PCA/correlation",sep=""),2400,2400)
    corrplot(M, method="circle")
    dev.off()
    
    #fit PCA
    res.pca <- prcomp(df,center = FALSE, scale. = FALSE)
    #res.pca <- prcomp(d.cpue,center = TRUE, scale. = FALSE)
    
    #Visualize all eigenvalues, from most to least contribution
    fn.fig(paste(HndL.Species_targeting,"PCA/eigenvalues",sep=""),2400,2400)
    fviz_eig(res.pca)
    dev.off()
    
    #Graph of records. Records with a similar profile are grouped together.
    fn.fig(paste(HndL.Species_targeting,"PCA/PCA_records",sep=""),2400,2400)
    g <- ggbiplot(res.pca, obs.scale = 1, var.scale = 1,
                  ellipse = TRUE, circle = TRUE)
    g <- g + scale_color_discrete(name = '')
    g= g 
    g <- g + theme(legend.direction = 'horizontal',legend.position = 'top')
    print(g)
    # fviz_pca_ind(res.pca,
    #                col.ind = "cos2", # Color by the quality of representation
    #                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    #                repel = TRUE)     # Avoid text overlapping
    dev.off()
    
    
    
    #Graph of each species cpue. Positive correlated variables point to the same side of the plot.
    #                     Negative correlated variables point to opposite sides of the grap
    fn.fig(paste(HndL.Species_targeting,"PCA/PCA_species",sep=""),2400,2400)
    fviz_pca_var(res.pca,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE )    # Avoid text overlapping
    dev.off()
    
    #Biplot of records and cpue    #too heavy
    # fviz_pca_biplot(res.pca, repel = TRUE,
    #                 col.var = "#2E9FDF", # cpues color
    #                 col.ind = "#696969")  # records color
    # 
    
    #Acces the results
    Eig.val=get_eigenvalue(res.pca)
    #res.var <- get_pca_var(res.pca)
    #res.var$coord          # Coordinates
    #res.var$contrib        # Contributions to the PCs
    #res.var$cos2           # Quality of representation 
    
    PCA.variance=Eig.val
    keep.axis=which(PCA.variance$cumulative.variance.percent<PercentExpl)
    PCA.scores=get_pca_ind(res.pca)$coord          # PCA Coordinates
    PCA.scores=data.frame(Same.return.SNo=row.names(PCA.scores),PCA.scores[,keep.axis])     #keep only relevant axes     
    row.names(PCA.scores)=NULL                                          
    PCA.axis=as.data.frame(res.pca$rotation)       # PCA axis loading. How species contribute to each axis
    
    #density
   # bivn.kde <- kde2d( PCA.scores$Dim.1, PCA.scores$Dim.2, n = 100)
  #  image(bivn.kde); contour(bivn.kde, add = T)
    
    #Plot variance explained
    fn.fig(paste(HndL.Species_targeting,"PCA/Variance",sep=""),2400,2400)
    plot(PCA.variance$cumulative.variance.percent,ylab="Cumulative  variance explained",xlab="Axis")
    dev.off()
    
    
    #Alocate axis values to each data set
    for(i in 1:N.species)
    {
      a=subset(PCA.scores,Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily_nfish[[i]]$Same.return.SNo))
      DATA.list.LIVEWT.c.daily_nfish[[i]]=merge(DATA.list.LIVEWT.c.daily_nfish[[i]],a,by="Same.return.SNo",all.x=T)
    }

    rm(ALL,d,d.cpue,a)
    detach("package:ggbiplot")  #remove package as it interfers with dplyr
  }
}



#4.16 Table of sensitivity scenarios       #ACA: redo accordingly!!!!
Tab.Sensi=data.frame(Scenario=c("Effective","Stand_Base case","Stand_Block","Stand_Vessel",
                                "Stand_Efficiency","Stand_Area weighting"),
         Records=c("All","Reliable","Reliable","Reliable","Reliable","Reliable"),
         Standardisation=c("No","Yes","Yes","Yes","Yes","Yes"),
         Blocks_used=c("All","Subset","All","Subset","Subset","Subset"),
         Vessels_used=c("All","Indicative","Indicative","All","Indicative","Indicative"),
         Efficiency_increase=c("Yes","Yes","Yes","Yes","No","Yes"),
         Area_weighting=c("No",rep("Yes",4),"No"))

setwd(paste(getwd(),"/Outputs/Paper",sep=""))
fn.word.table(WD=getwd(),TBL=Tab.Sensi,Doc.nm="Sensitivity tests",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


#4.17 Compute foly index for exporting       

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


List.foly.nom=list(whis=list(Foly=Foly.whis),
                   gum=list(Foly=Foly.gum),
                   dus=list(Foly=Foly.dus),
                   san=list(Foly=Foly.san))


#4.18 Compute nominal for exporting 
#note: Ratio = mean(catch)/mean(effort)
#      Mean = mean(cpue)
#     LnMean= exp(mean(log(cpue))+bias corr)
#     DLnMean = exp(log(prop pos)+exp(mean(log(cpue))+bias corr)

QL_expl_ktch_prop=.9   #proportion of explained annual catch for selecting record
PRP.MLCLM=0.1          #proportion of catch of target

Store_nom_cpues_monthly_km.gn.d=Store_nom_cpues_daily_km.gn.d=
Store_nom_cpues_monthly_km.gn.h=Store_nom_cpues_daily_km.gn.h=
Store_nom_cpues_daily_km.gn.d_nfish=Store_nom_cpues_daily_km.gn.h_nfish=DATA.list.LIVEWT.c.daily
Hnd.ains="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Ainsline_different_cpues/"
  
for(s in 1:N.species)
{
  #Monthly
  Efrt='km.gillnet.days.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_monthly_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_monthly_km.gn.d[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "kg/km gillnet days",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used[[s]]),VesL=VES.used[[s]],Show.QL=TRUE)
  dev.off() 
  
  Efrt='km.gillnet.hours.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_monthly_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_monthly_km.gn.h[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "kg/km gillnet hours",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used[[s]]),VesL=VES.used[[s]],Show.QL=FALSE)
  dev.off() 
  
  
  #Daily weight
  Efrt='km.gillnet.days.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_daily_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_daily_km.gn.d[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "kg/km gillnet days",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Show.QL=TRUE)
  dev.off() 
  
  Efrt='km.gillnet.hours.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_daily_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_daily_km.gn.h[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "kg/km gillnet hours",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Show.QL=FALSE)
  dev.off() 
  
  #Daily numbers
  Efrt='km.gillnet.days.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_daily_nfish_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_daily_km.gn.d_nfish[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily_nfish[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "Individuals/km gillnet days",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Show.QL=FALSE)
  dev.off() 
  
  Efrt='km.gillnet.hours.c'
  pdf(paste(Hnd.ains,SPECIES.vec[s],"_daily_nfish_",paste(Efrt),".pdf",sep="")) 
  Store_nom_cpues_daily_km.gn.h_nfish[[s]]=fn.ainslie(dat=DATA.list.LIVEWT.c.daily_nfish[[s]],Ktch.targt='catch.target',
            Effrt=Efrt,explr="NO",QL_prop_ktch=QL_expl_ktch_prop,Prop.Malcolm=PRP.MLCLM,
            cpue.units = "Individuals/km gillnet hours",spname=SPECIES.vec[s],
            BLks=as.numeric(BLKS.used.daily[[s]]),VesL=VES.used.daily[[s]],Show.QL=FALSE)
  dev.off() 
}


#Evaluate balance of data subset based on QL  #ACA
#note: remove blocks with less than Min.Vess.yr (monthly) and Min.Vess.yr.d (daily)
#      for those blocks, keep vessels with more than Min.Vess.yr / Min.Vess.yr.d

BLKS.used.indi=BLKS.used  #keep copy of indicators approach
VES.used.indi=VES.used
BLKS.used.daily.indi=BLKS.used.daily
VES.used.daily.indi=VES.used.daily

HndL="C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/QL_balanced_design/"

for(s in 1:N.species)
{
  fn.check.balanced=function(d,SP,what,MN.YR)
  {
    fn.plt=function(dd)
    {
      a=dd
      a[a>0]=1
      Orderd=rev(sort(rowSums(a)))
      dd=as.matrix(dd[match(names(Orderd),row.names(dd)),])
      Nx=c(1,ncol(dd))
      Ny=c(1,nrow(dd)+3)
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
    
    pdf(paste(HndL,paste(SP,"_",what,sep=""),".pdf",sep="")) 
    
    #First select blocks
    BLK.Yr=with(d,table(blockx,finyear))
    this.blks=fn.plt(BLK.Yr)
    mtext("Block (all)",2,line=1.65,cex=1.25)
    
    #Second, select vessels within selected blocks
    Ves.Yr=with(subset(d,blockx%in%this.blks),table(vessel,finyear))
    this.ves=fn.plt(Ves.Yr)
    mtext("Vessel (for selected blocks)",2,line=1.65,cex=1.25)
    
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
    dev.off()
    
    return(list(this.blks=this.blks,this.ves=this.ves))
    
     
  }
  a=fn.check.balanced(d=Store_nom_cpues_monthly_km.gn.d[[s]]$QL_dat,
                    SP=SPECIES.vec[s],what="monthly",MN.YR=Min.Vess.yr)
  BLKS.used[[s]]=a$this.blks
  VES.used[[s]]=a$this.ves
  
  a=fn.check.balanced(d=Store_nom_cpues_daily_km.gn.d[[s]]$QL_dat,
                    SP=SPECIES.vec[s],what="daily",MN.YR=Min.Vess.yr.d)
  BLKS.used.daily[[s]]=a$this.blks
  VES.used.daily[[s]]=a$this.ves
}



setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")


#Show gummy monthly cpue effect of using km gn d or km g h
fn.fig("Gummy_km.gn.d_VS_km.gn.h",2400,2400)
par(mfrow=c(3,2),mar=c(1,3,2,1),oma=c(2,.5,.5,1.75),mgp=c(1.5,.5,0),cex.lab=1.5)
Get.Mns(d=DATA.list.LIVEWT.c$gummy,grp="FINYEAR",
        Vars=c("HOURS.c","SHOTS.c",'Km.Gillnet.Days.c','Km.Gillnet.Hours.c','Catch.Target'),
        LGND=c("Hours fished per day","Number of shots",
               'km gillnet days','km gillnet hours','Catch (kg)'),
        add.arrow=c(rep(TRUE,4),FALSE))

dev.off()


if(check.interactions=="YES")
{
  for(i in 1:N.species) fn.vessl.cpue.yr(DATA.list.LIVEWT.c[[i]])
  for(i in 1:N.species) fn.vessl.cpue.yr(DATA.list.LIVEWT.c.daily[[i]])
}


#export data to ainslie
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



#4.19 Output data tables
#Output table with number of records available in effective area, numbers used in standardisation
TABle=vector('list',N.species)
names(TABle)=SPECIES.vec
for(s in 1:N.species)
{
  
  #total records
  Tot.m=table(DATA.list.LIVEWT.c_all_reporters[[s]]$FINYEAR)
  Tot.d=table(DATA.list.LIVEWT.c.daily_all_reporters[[s]]$FINYEAR)
  
  Good.m=table(DATA.list.LIVEWT.c[[s]]$FINYEAR)
  Good.d=table(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR)
  
  
  #Used in Standardisation after selection of blocks and vessels
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



#Table of species by zone 
#monthly
Zones.whi=table(DATA.list.LIVEWT.c[[1]]$zone)
Zones.gum=table(DATA.list.LIVEWT.c[[2]]$zone)
Zones.dus=table(DATA.list.LIVEWT.c[[3]]$zone)
Zones.san=table(DATA.list.LIVEWT.c[[4]]$zone)

#Table of years by block for each species within its' effective area
#monthly
Yrs.blk.whis=Yrs.by.blk(DATA.list.LIVEWT.c[[1]])
Yrs.blk.gum=Yrs.by.blk(DATA.list.LIVEWT.c[[2]])
Yrs.blk.dus=Yrs.by.blk(DATA.list.LIVEWT.c[[3]])
Yrs.blk.san=Yrs.by.blk(DATA.list.LIVEWT.c[[4]])

#Blocks-Years with catch
Yr_Blk.ktch.whis=ncol(Yrs.blk.whis$Table1)
Yr_Blk.ktch.gum=ncol(Yrs.blk.gum$Table1)
Yr_Blk.ktch.dus=ncol(Yrs.blk.dus$Table1)
Yr_Blk.ktch.san=ncol(Yrs.blk.san$Table1)

#Vessels-Years with catch
Yr_Ves.ktch.whis=ncol(Yrs.blk.whis$Table2)
Yr_Ves.ktch.gum=ncol(Yrs.blk.gum$Table2)
Yr_Ves.ktch.dus=ncol(Yrs.blk.dus$Table2)
Yr_Ves.ktch.san=ncol(Yrs.blk.san$Table2)

#Year-months with catch
Yr_Mn.ktch.whis=Yrs.blk.whis$Table3
Yr_Mn.ktch.gum=Yrs.blk.gum$Table3
Yr_Mn.ktch.dus=Yrs.blk.dus$Table3
Yr_Mn.ktch.san=Yrs.blk.san$Table3

#Blocks with catch
Blk.ktch.whis=length(Yrs.blk.whis$Table4)
Blk.ktch.gum=length(Yrs.blk.gum$Table4)
Blk.ktch.dus=length(Yrs.blk.dus$Table4)
Blk.ktch.san=length(Yrs.blk.san$Table4)

#Vessels with catch
Ves.ktch.whis=length(Yrs.blk.whis$Table5)
Ves.ktch.gum=length(Yrs.blk.gum$Table5)
Ves.ktch.dus=length(Yrs.blk.dus$Table5)
Ves.ktch.san=length(Yrs.blk.san$Table5)


#Blocks with no catch
Blk.No.ktch.whis=length(Yrs.blk.whis$Table6)
Blk.No.ktch.gum=length(Yrs.blk.gum$Table6)
Blk.No.ktch.dus=length(Yrs.blk.dus$Table6)
Blk.No.ktch.san=length(Yrs.blk.san$Table6)

#Vessels with no catch
Ves.No.ktch.whis=Yrs.blk.whis$Table7
Ves.No.ktch.gum=Yrs.blk.gum$Table7
Ves.No.ktch.dus=Yrs.blk.dus$Table7
Ves.No.ktch.san=Yrs.blk.san$Table7


Ves.no.ktc=list(Ves.No.ktch.whis,Ves.No.ktch.gum,Ves.No.ktch.dus,Ves.No.ktch.san)
Numb.ves.no.ktc=list(whis=length(Ves.No.ktch.whis),gum=length(Ves.No.ktch.gum),
                     dus=length(Ves.No.ktch.dus),san=length(Ves.No.ktch.san))



#4.20 Plot cpue by year and month
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Month")
for(i in 1:N.species)
{
  FILE=paste("Mn_YR.cpue.",SPECIES.vec[i],sep="")
  fn.fig(FILE,2400, 2400)    
  
  if(Separate.monthly.daily=="NO")
  {
    plt.ktch.yr_mn(subset(DATA.list.LIVEWT.c[[i]],select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,MONTH)))
  }
  if(Separate.monthly.daily=="YES")
  {
    plt.ktch.yr_mn(
      rbind(subset(DATA.list.LIVEWT.c[[i]],select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,MONTH)),
            subset(DATA.list.LIVEWT.c.daily[[i]],select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,MONTH))))
  }
  dev.off()
} 
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")


#4.21 Check outliers in catch and effort for removing nonsense values   
  #monthly
Check.outliers=DATA.list.LIVEWT.c
for(i in 1:length(Check.outliers))Check.outliers[[i]]=check.cpue(DATA.list.LIVEWT.c[[i]],names(DATA.list.LIVEWT.c)[i],BYDAY="NO")



#4.22 Construct index of abundance
ZONES=c("West","Zone1","Zone2")
Predictors=c("FINYEAR","BLOCKX","MONTH","VESSEL","Catch.Gummy",
             "Catch.Whiskery","Catch.Dusky","Catch.Sandbar","Catch.Scalefish",
             "SOI","Freo","Freo.Lag6","Freo.Lag12","block10") 
Eff.vars=c("Km.Gillnet.Days.c","Km.Gillnet.Hours.c")
VARIABLES=c("Catch.Target",Predictors)

#Zone.list=list(c("West","Zone1","Zone2"),c("Zone1","Zone2"),c("West","Zone1","Zone2"),c("West","Zone1","Zone2"))

#remove catch column of target species
drop.pred=c(match("Catch.Whiskery",Predictors),match("Catch.Gummy",Predictors),
            match("Catch.Dusky",Predictors),match("Catch.Sandbar",Predictors))


    #4.22.1 Explore data sets
if(Model.run=="First")
{
  #Number of observations per Month-Year-Block
    #Monthly
  Deg.free=DATA.list.LIVEWT.c
  for(i in 1:length(Deg.free)) Deg.free[[i]]=Check.deg.free(DATA.list.LIVEWT.c[[i]])
  
    #Daily
  Deg.free.daily=DATA.list.LIVEWT.c
  for(i in 1:length(Deg.free.daily)) Deg.free.daily[[i]]=Check.deg.free(DATA.list.LIVEWT.c.daily[[i]])
  
    #Exploratory
  RecordThreshold=99
  ZONAS=c("West","Zone1","Zone2")
  Data.Summary=vector('list',length=N.species)
  for ( i in 1:N.species)
  {
    fn.fig(paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Exploratory.monthly",SPECIES.vec[i],
                 sep=""),2400,2400)  
    par(mfrow=c(3,1),mar=c(4,4,1,1),las=1)
    Data.Summary[[i]]=fn.explore(DATA.list.LIVEWT.c[[i]],SPECIES.vec[i])
    mtext(SPECIES.vec[i], side = 3, line = -1, outer = T)
    dev.off()
  }

  #4.22.2 Catch VS CPUE 
  #note: IF catch increases and CPUE decreases then Effective Effort or Q increased!
  for ( i in 1:N.species)
  {
    fn.fig(paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Catch_CPUE.monthly.",SPECIES.vec[i],
                 sep=""),2400, 2400)  
    par(mfrow=c(1,1),mar=c(1,3.6,1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
    fn.catch.cpue(DATA.list.LIVEWT.c[[i]])
    dev.off()
    
    fn.fig(paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Catch_CPUE.daily.",SPECIES.vec[i],
                 sep=""),2400,2400)  
    par(mfrow=c(1,1),mar=c(1,3.6,1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
    fn.catch.cpue(DATA.list.LIVEWT.c.daily[[i]])
    dev.off()
  }  
  
   #Check if splitting peak in cpue by using financial year
  #note: don't split seasonal peak of catch/cpue into different years
  fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Splitting catch peak or not/Monthly",2400,2400)  
  par(mfcol=c(2,2),mar=c(1,.1,1.5,4),oma=c(4,5,.1,.5),las=1,mgp=c(1.9,.6,0))
  for(i in 1:length(SPECIES.vec))fn.see.month(d=DATA.list.LIVEWT.c[[i]])
  mtext("Average cpue (SE)",4,-1,outer=T,las=3,cex=1.35,col=2)
  mtext("Average catch (SE)",2,3,las=3,outer=T,cex=1.35)
  mtext("Month",1,1,outer=T,cex=1.35)
  mtext("Monthly returns",3,-1.5,outer=T,cex=1.35)
  dev.off()
  
  fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Splitting catch peak or not/Daily",2400,2400)  
  par(mfcol=c(2,2),mar=c(1,.1,1.5,4),oma=c(4,5,.1,.5),las=1,mgp=c(1.9,.6,0))
  for(i in 1:length(SPECIES.vec))fn.see.month(d=DATA.list.LIVEWT.c.daily[[i]])
  mtext("Average cpue (SE)",4,-1,outer=T,las=3,cex=1.35,col=2)
  mtext("Average catch (SE)",2,3,las=3,outer=T,cex=1.35)
  mtext("Month",1,1,outer=T,cex=1.35)
  mtext("Daily logbooks",3,-1.5,outer=T,cex=1.35)
  dev.off()
}

    #4.22.2 Combine vessels in categories
if(Combine.ves=="YES")
{
  #Monthly
  New.Ves.Lev=vector('list',N.species)
  names(New.Ves.Lev)=names(DATA.list.LIVEWT.c)
  
  Store.New.Ves=Store.New.Ves.sens=DATA.list.LIVEWT.c
  for(i in 1:N.species)
  {
    Store.New.Ves[[i]]=fn.new.ves(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Per.int.bc)
    Store.New.Ves.sens[[i]]=fn.new.ves(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Per.int.sens)
  }
  
  
  #plot predicted catch by new vessel category
  fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Paper/GLMM.catch.pred.vessel.cat_Monthly",2400,2400)    
  par(mfcol=c(2,2),mar=c(1,3.7,.75,.1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
  for(i in 1:N.species)fn.plot.nw.ves.cat(Store.New.Ves[[i]],names(Store.New.Ves)[i])
  mtext("Vessel category",side=1,line=1,outer=T,cex=1.5)
  mtext("Mean predicted catch (kg) for a standardised effort level (+/- SD)",side=2,
        line=-1.25,outer=T,cex=1.25,las=3)
  dev.off()
  
  
  #Add new vessel categories to data
  #Monthly
  DATA.list.LIVEWT.c.grup.sens=DATA.list.LIVEWT.c.indi.vess=DATA.list.LIVEWT.c
  for ( i in 1:N.species)
  {
    DATA.list.LIVEWT.c.grup.sens[[i]]=fn.add.new.ves(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Store.New.Ves.sens[[i]])
    DATA.list.LIVEWT.c[[i]]=fn.add.new.ves(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Store.New.Ves[[i]])
    DATA.list.LIVEWT.c.indi.vess[[i]]=subset(DATA.list.LIVEWT.c[[i]],VESSEL%in%Indicative.vessels[[i]]$These.vessels)
  }
  
  #Daily
  #note: for daily, set new.level.vessel to the original vessel variable
  for ( i in 1:N.species) DATA.list.LIVEWT.c.daily[[i]]$new.level.vessel=DATA.list.LIVEWT.c.daily[[i]]$VESSEL
  
  #Explore predictor's effect
  if(Model.run=="First")
  {
    setwd("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/Predictor.effect")
    Missing.Values=vector('list',N.species)
    names(Missing.Values)=SPECIES.vec
    
    for(i in 1:N.species)
    {
      #drop vessels that never reported catch of target species
      DATA=subset(DATA.list.LIVEWT.c[[i]],!(VESSEL%in%Ves.no.ktc[[i]])) 
      
      DATA$VESSEL=as.factor(DATA$VESSEL)
      DATA$FINYEAR=as.factor(DATA$FINYEAR)
      DATA$MONTH=as.factor(DATA$MONTH)
      DATA$BLOCKX=as.factor(DATA$BLOCKX)
      DATA$zone=as.factor(DATA$zone)
      DATA$mesh=as.factor(DATA$mesh)
      PREDS=Predictors[-match("Km.Gillnet.Days.c",Predictors)]
      PREDS=PREDS[-drop.pred[i]]
      PREDS=PREDS[-match("block10",PREDS)]
      DATA$cpue=DATA$Catch.Target/DATA$Km.Gillnet.Days.c
      
      fn.fig(paste("Exp.",SPECIES.vec[i],".General",sep=""),2400,2400)
      par(mfcol=c(3,2),mai=c(.6,.65,.3,.1),oma=c(.2,.2,.1,.1))
      hist(log(DATA$cpue),main="hist log(cpue)",xlab="log(cpue)",ylab="Count")
      
      boxcox(cpue+0.00001 ~ log(YEAR.c), data = DATA,lambda = seq(0, 1, length = 10))
      legend("topright","Box Cox (should be small)",bty='n')
      
      # Cook distance to see outliers or overdisperse data (if distance >1)
      M1.1=glm(log(cpue+0.00001)~FINYEAR,family=gaussian,data=DATA)
      plot(M1.1,which=4)
      legend("topright","outliers or overdispersion if distance >1",bty='n',cex=.85, text.col=2)
      
      plot(table(DATA$Catch.Target),type='h',xlab="Catch",ylab="Count",main="Catch zero inflation and right tail")
      
      # Search for missing values per predictor
      Preds=match(Predictors,names(DATA))
      Preds=subset(Preds,!is.na(Preds))
      Missing.Values[[i]]=Predictors_missing_values=colSums(is.na(DATA[,Preds]))/nrow(DATA)
      
      #Outliers response var
      boxplot(DATA$cpue~DATA$FINYEAR,main="Outliers in response var?",ylab="cpue (kg/km.gn.day)")
      dev.off()
      
      
      #boxplot of response var and predictors
      fn.fig(paste("Exp.",SPECIES.vec[i],".Box.plot",sep=""),2400, 2400)
      par(mfcol=c(4,3),mai=c(.45,.5,.2,.1),oma=c(.1,.1,.1,.1),mgp=c(1.1,.35,0),cex.main=1)
      for(d in 1:length(PREDS))
      {
        a=DATA[,match(c("cpue",PREDS[d]),names(DATA))]
        if(!is.factor(a[,2])) a[,2]=cfac(a[,2])
        plot(clog(a[,1])~a[,2],main=PREDS[d],ylab="",xlab="")
      }
      mtext("log(cpue)",side=2,line=-1.5,las=3,outer=T)
      dev.off()
      
      
      #   tiff(paste("Exp.",SPECIES.vec[i],".Design.plot.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      #   plot.design(clog(Catch.Target)~FINYEAR+MONTH+BLOCKX+VESSEL+cfac(SOI)+cfac(Freo)+
      #                 cfac(round(Catch.Gummy))+cfac(round(Catch.other.shk))+cfac(round(Catch.other.shk))
      #               +cfac(round(Catch.other.shk)),data=DATA,cex.lab=1.5,ylab="log(Catch)",
      #               main="Design plot of Predictors",xlab="",las=2,cex.axis=.750)
      #   dev.off()    
      
      #   tiff(paste("Exp.",SPECIES.vec[i],".Interaction.plot.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      #   interaction.plot(DATA$FINYEAR,DATA$zone,clog(DATA$Catch.Target),ylab="log(Catch)",
      #          xlab="FINYEAR",main="Year-zone interaction",fixed=TRUE, col = 1:100, leg.bty = "o")
      #   dev.off() 
      
      
      #Covariate correlations.
      PREDS.cont=PREDS[-match(c("FINYEAR","BLOCKX","VESSEL"),PREDS)]
      fn.fig(paste("Exp.",SPECIES.vec[i],".Correlations.plot",sep=""),2400, 2400)
      Covars=DATA[,match(PREDS.cont,names(DATA))]
      Covars$MONTH=as.numeric(Covars$MONTH)
      pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor)
      dev.off()
    }
    
    Missing.Values  #No 0s?
    rm(M1.1)
    
  }
  
  #4.23.6 Define best model and error structure 
  VARIABLES=c(VARIABLES,"new.level.vessel")
  
}

setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")

    #4.22.3  Define model structure with Yr-Blk interactions    #takes 29 mins (in sequential computation)
if(Def.mod.Str.=="YES" & With.interact=="YES")
{
  #table 1. Model terms
  Tab1.Mon=vector('list',length(SPECIES.vec))
  names(Tab1.Mon)=SPECIES.vec
  Tab1.Day=Tab1.Mon
  
  #Monthly
  for(j in 1:N.species) Tab1.Mon[[j]]=fn.Tabl1(DATA.list.LIVEWT.c[[j]])
  write.csv(do.call(rbind,Tab1.Mon),"Table1_Month.csv")
  
  #Daily
  for(j in 1:N.species) Tab1.Day[[j]]=fn.Tabl1.daily(DATA.list.LIVEWT.c.daily[[j]])
  write.csv(do.call(rbind,Tab1.Day),"Table1_Daily.csv")
  
  #Covariate correlations
  Cor.Ktch=c("Catch.Gummy","Catch.Whiskery", "Catch.Scalefish", "Catch.Dusky", "Catch.Sandbar",
             "SOI","Freo")
  if(Chck.Corr=="YES")
  {
    for(j in 1:N.species)
    {
      d=DATA.list.LIVEWT.c[[j]][,match(Cor.Ktch,names(DATA.list.LIVEWT.c[[j]]))]
      pairs(d, lower.panel=panel.smooth, upper.panel=panel.cor,main=paste(names(DATA.list.LIVEWT.c)[j],"monthly"))
      
      d=DATA.list.LIVEWT.c.daily[[j]][,match(Cor.Ktch,names(DATA.list.LIVEWT.c.daily[[j]]))]
      pairs(d, lower.panel=panel.smooth, upper.panel=panel.cor,main=paste(names(DATA.list.LIVEWT.c.daily)[j],"daily"))
      rm(d)
    }
  }
  
  
  #tested models
  Species.model=vector('list',length=N.species)
  names(Species.model)=SPECIES.vec
  Species.model.sel=Species.model
  N.models=length(Predictors[-match(c("FINYEAR","BLOCKX","VESSEL","Km.Gillnet.Days.c"),Predictors)])
  
  Mon.wd="C:/Matias/Analyses/Catch and effort/Outputs/ModelStructure/monthly"
  Dai.wd="C:/Matias/Analyses/Catch and effort/Outputs/ModelStructure/daily"
  
  for(w in 1:length(Fit.to.what))
  {
    Species.model=fn.define.models(Fit.to.what[w],Interaction="+")    
    All.models=vector('list',length(SPECIES.vec))
    names(All.models)=SPECIES.vec
    All.models.daily=All.models
    DAT=names(DATA.list.LIVEWT.c)  
    LEVELS.ves.blk=All.models
    
    #Monthly
    DIR=paste(Mon.wd,"/",Fit.to.what[w],sep='')
    setwd(DIR)
    
    #run loop in multiple cores
    for ( j in 1:N.species)
    {
      bin.form=Species.model[[j]]$bin.form
      pos.log.form=Species.model[[j]]$pos.log.form
      pos.form=Species.model[[j]]$pos.form
      
      Test.models(DATA.list.LIVEWT.c[[j]],Ves.no.ktc[[j]],SPECIES.vec[j],Min.weight[j])
      Blk.ves.levels=c(length(unique(DATA.list.LIVEWT.c[[j]]$VESSEL)),length(unique(DATA.list.LIVEWT.c[[j]]$BLOCKX)))
      names(Blk.ves.levels)=c("Vessels","Blocks")
      LEVELS.ves.blk[[j]]=Blk.ves.levels
    }
    AIC.bin=Dev.bin=rep(NA,N.models)
    AIC.Pos=Dev.pos=rep(NA,N.models)
    AIC.Pos.Gamma=Dev.pos.Gamma=rep(NA,N.models)
    
    #Best models
    Monthly.best=load.and.select(DIR,1:N.species)
    Selected.pos=vector('list',length=N.species)
    names(Selected.pos)=names(Monthly.best)
    Selected.bin=Selected.pos.Gamma=Selected.pos
    
    #get selected models 
    for ( i in 1:N.species)  
    {
      TTERMS=all.vars(Species.model[[i]]$bin.form[[N.models]])
      
      #bin part
      Dev.expl.bin=round(Monthly.best[[i]]$Dev.explained$bi,3)
      
      Terms.names=rep(NA,length(Dev.expl.bin))
      for(r in 1:length(Terms.names))
      {
        ll=all.vars(Species.model[[i]]$bin.form[[r]])
        Terms.names[r]=ll[length(ll)-1]
      }
      
      names(Dev.expl.bin)=Terms.names
      AIC.mdls.bin=round(Monthly.best[[i]]$AIC.bin)      
      sel.mod=Dev.expl.bin
      sel.mod[1]=0
      for(m in 2:length(Dev.expl.bin))sel.mod[m]=Dev.expl.bin[m]-Dev.expl.bin[m-1]
      sel.mod[sel.mod<0]=0
      
      BESTT=which(sel.mod>(Red.th/100))
      #       if(length(BESTT)==0)
      #       {
      #         sel.mod[2:length(sel.mod)]=round(Dev.expl.bin[2:length(sel.mod)]-Dev.expl.bin[1],2)
      #         BESTT=which(sel.mod>(Red.th/100))
      #         BESTT=which(sel.mod==max(sel.mod[BESTT]))[1]
      #       }
      Selected.bin[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")      
      
      #pos part
      Dev.expl.pos=round(Monthly.best[[i]]$Dev.explained$pos,3)
      names(Dev.expl.pos)=Terms.names
      AIC.mdls.pos=round(Monthly.best[[i]]$AIC.Pos)
      sel.mod=Dev.expl.pos
      sel.mod[1]=0
      for(m in 2:length(Dev.expl.pos))sel.mod[m]=Dev.expl.pos[m]-Dev.expl.pos[m-1]
      sel.mod[sel.mod<0]=0
      BESTT=which(sel.mod>(Red.th/100))
      #       if(length(BESTT)==0)
      #       {
      #         sel.mod[2:length(sel.mod)]=round(Dev.expl.pos[2:length(sel.mod)]-Dev.expl.pos[1],2)
      #         BESTT=which(sel.mod>(Red.th/100))
      #         BESTT=which(sel.mod==max(sel.mod[BESTT]))[1]
      #       }      
      Selected.pos[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")
      
      #pos part Gamma
      Dev.expl.pos=round(Monthly.best[[i]]$Dev.explained$pos.Gamma,3)
      names(Dev.expl.pos)=Terms.names
      AIC.mdls.pos=round(Monthly.best[[i]]$AIC.Pos.Gamma)
      sel.mod=Dev.expl.pos
      sel.mod[1]=0
      for(m in 2:length(Dev.expl.pos))sel.mod[m]=Dev.expl.pos[m]-Dev.expl.pos[m-1]
      sel.mod[sel.mod<0]=0      
      BESTT=which(sel.mod>(Red.th/100))
      Selected.pos.Gamma[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")      
    }
    rm(Monthly.best)
    
    
    #Daily
    if(Separate.monthly.daily=="YES")
    {
      DIR=paste(Dai.wd,"/",Fit.to.what[w],sep='')
      setwd(DIR)
      for ( j in 1:N.species)       
      {
        bin.form=Species.model[[j]]$bin.form
        pos.log.form=Species.model[[j]]$pos.log.form
        pos.form=Species.model[[j]]$pos.form        
        Test.models(DATA.list.LIVEWT.c.daily[[j]],Ves.no.ktc[[j]],SPECIES.vec[j],Min.weight[j])
      }
      
      #Best models
      Daily.best=load.and.select(DIR,1:N.species)
      
      Selected.pos.daily=vector('list',length=N.species)
      names(Selected.pos.daily)=names(Daily.best)
      Selected.bin.daily=Selected.pos.Gamma.daily=Selected.pos.daily
      
      #get selected models
      for ( i in 1:N.species)      
      {
        TTERMS=all.vars(Species.model[[i]]$bin.form[[N.models]])
        
        #bin part
        Dev.expl.bin=round(Daily.best[[i]]$Dev.explained$bi,3)
        Terms.names=rep(NA,length(Dev.expl.bin))
        for(r in 1:length(Terms.names))
        {
          ll=all.vars(Species.model[[i]]$bin.form[[r]])
          Terms.names[r]=ll[length(ll)-1]
        }
        
        names(Dev.expl.bin)=Terms.names  
        AIC.mdls.bin=round(Daily.best[[i]]$AIC.bin)      
        sel.mod=Dev.expl.bin
        sel.mod[1]=0
        for(m in 2:length(Dev.expl.bin))sel.mod[m]=Dev.expl.bin[m]-Dev.expl.bin[m-1]
        sel.mod[sel.mod<0]=0        
        BESTT=which(sel.mod>(Red.th/100))        
        Selected.bin.daily[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")      
        
        #pos part
        Dev.expl.pos=round(Daily.best[[i]]$Dev.explained$pos,3)
        names(Dev.expl.pos)=Terms.names
        AIC.mdls.pos=round(Daily.best[[i]]$AIC.Pos)
        sel.mod=Dev.expl.pos
        sel.mod[1]=0
        for(m in 2:length(Dev.expl.pos))sel.mod[m]=Dev.expl.pos[m]-Dev.expl.pos[m-1]
        sel.mod[sel.mod<0]=0
        BESTT=which(sel.mod>(Red.th/100))
        Selected.pos.daily[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")
        
        #pos part Gamma
        Dev.expl.pos=round(Daily.best[[i]]$Dev.explained$pos.Gamma,3)
        names(Dev.expl.pos)=Terms.names
        AIC.mdls.pos=round(Daily.best[[i]]$AIC.Pos.Gamma)
        sel.mod=Dev.expl.pos
        sel.mod[1]=0
        for(m in 2:length(Dev.expl.pos))sel.mod[m]=Dev.expl.pos[m]-Dev.expl.pos[m-1]
        sel.mod[sel.mod<0]=0
        BESTT=which(sel.mod>(Red.th/100))
        Selected.pos.Gamma.daily[[i]]=c(TTERMS[1:4],TTERMS[BESTT+3],"log.Effort")  
      }
      
      rm(Daily.best)
    }
    
    #Best model by zone (whiskery and dusky)
    Do.by.Zone="NO"
    if(Do.by.Zone=="YES")
    {
      Whi.Dus=match(c("Whiskery shark","Dusky shark"),SPECIES.vec)
      All.models=vector('list',length(Whi.Dus))
      names(All.models)=SPECIES.vec[Whi.Dus]
      
      zn=c("West","Zone1","Zone2")
      Mon.wd.zone=list("C:/Matias/Analyses/Catch and effort/Outputs/ModelStructure/monthly/by zone/West",
                       "C:/Matias/Analyses/Catch and effort/Outputs/ModelStructure/monthly/by zone/Zone1",
                       "C:/Matias/Analyses/Catch and effort/Outputs/ModelStructure/monthly/by zone/Zone2")
      n.zn=length(zn)
      system.time(for(j in Whi.Dus)
      {
        bin.form=Species.model[[j]]$bin.form
        pos.log.form=Species.model[[j]]$pos.log.form
        
        for(z in 1:n.zn)
        {
          setwd(Mon.wd.zone[[z]])
          dat=subset(DATA.list.LIVEWT.c[[j]],zone==zn[z])
          Test.models(dat,Ves.no.ktc[[j]],SPECIES.vec[j],Min.weight[j])
        }
      })
      
      for(z in 1:n.zn)
      {
        Monthly.best=load.and.select(Mon.wd.zone[[z]],Whi.Dus)
        print(paste("___________",zn[z],"___________"))
        for(i in 1:N.species)
        {
          print("______________________________")
          print(names(Monthly.best)[i])
          print(round(Monthly.best[[i]]$Dev.explained$bi,2))
          print(round(Monthly.best[[i]]$AIC.bin))
          print(Species.model[[i]]$bin.form)
          print("___________")
          print(round(Monthly.best[[i]]$Dev.explained$pos,2))
          print(round(Monthly.best[[i]]$AIC.Pos))
          print(Species.model[[i]]$pos.log.form)
          print("______________________________")
        }
      }
      
    }
    
    #re do sandbars for area subsets
    re.do.san="NO"
    if(re.do.san=="YES")
    {
      San.areas=c("North.30.","North.32.","South.32.")
      San.areas.range=list(c(-30,-26),c(-32,-26),c(-36,-32))
      names(San.areas.range)=San.areas
      Monthly.best=San.areas.range
      for(j in 4:4)
      {
        bin.form=Species.model[[j]]$bin.form
        pos.log.form=Species.model[[j]]$pos.log.form
        
        for(w in 1: length(San.areas))
        {
          Q=San.areas[w]
          dummy=subset(DATA.list.LIVEWT.c[[j]],LAT>=(San.areas.range[[w]][1]) & LAT<(San.areas.range[[w]][2]))
          Monthly.best[[w]]=Test.models.sandbar(dummy,Ves.no.ktc[[j]],paste(Q,SPECIES.vec[j],sep=""),Min.weight[j])      
        }
      }
      
      DEV.bin=DEV.pos=Monthly.best
      for(w in 1: length(San.areas))
      {
        Dev.bin=rep(NA,N.models)
        Dev.pos=rep(NA,N.models)
        for(i in 1:N.models)
        {
          Dev.bin[[i]]=fn.Dev.Exp(Monthly.best[[w]]$bi[[i]])
          Dev.pos[[i]]=fn.Dev.Exp(Monthly.best[[w]]$pos[[i]])
        }
        DEV.bin[[w]]=Dev.bin
        DEV.pos[[w]]=Dev.pos
        
      }
      
    }
  }
  
  #Print out
  Selected.bin
  Selected.pos
  Selected.pos.Gamma
  
  if(Separate.monthly.daily=="YES")
  {
    Selected.bin.daily
    Selected.pos.daily
    Selected.pos.Gamma.daily 
  }
  
}

    #4.22.4 Write out best model
Best.Model=vector('list',N.species)
names(Best.Model)=SPECIES.vec

  #Monthly
Best.Model$'Whiskery shark'=list(
  Bi=as.formula(Catch.Target~FINYEAR+VESSEL+MONTH+BLOCKX+offset(log.Effort)), 
  Log=as.formula(log.Catch~FINYEAR+VESSEL+MONTH+BLOCKX+offset(log.Effort)))

Best.Model$'Gummy shark'=Best.Model$'Dusky shark'=Best.Model$'Sandbar shark'=Best.Model$'Whiskery shark'


  #Daily
Best.Model.daily=Best.Model


#Check Proportion of records with mesh data   
if(Model.run=="First")
{
  Mesh.Dat.prop.no.mesh=Species.list
  Mesh.Dat.prop.no.mesh.daily=Species.list.daily
  for(i in 1:length(Mesh.Dat.prop))
  {
    d=table(Mesh.Dat.prop.no.mesh[[i]]$mesh,useNA='ifany')
    Mesh.Dat.prop.no.mesh[[i]]=as.numeric(d[match(NA,names(d))]/sum(d))
    d=table(Mesh.Dat.prop.no.mesh.daily[[i]]$mesh,useNA='ifany')
    Mesh.Dat.prop.no.mesh.daily[[i]]=as.numeric(d[match(NA,names(d))]/sum(d))
    
  }
  Mesh.Dat.prop.no.mesh=round(do.call(cbind,Mesh.Dat.prop.no.mesh),2)
  Mesh.Dat.prop.no.mesh.daily=round(do.call(cbind,Mesh.Dat.prop.no.mesh.daily),2)
  
}

      #4.22.5 Get blocks from each time series
    #Monthly
BLOCKX.zn=DATA.list.LIVEWT.c
for ( i in 1:N.species)BLOCKX.zn[[i]]=fn.blk.zone("BLOCKX",DATA.list.LIVEWT.c[[i]])
BLOCKX.zn=do.call(rbind,BLOCKX.zn)
BLOCKX.zn=as.data.frame.matrix(with(BLOCKX.zn,table(BLOCKX,zone)))
BLKZ.per.zone=list(West=as.numeric(rownames(subset(BLOCKX.zn,West>0))),
                   Zone1=as.numeric(rownames(subset(BLOCKX.zn,Zone1>0))),
                   Zone2=as.numeric(rownames(subset(BLOCKX.zn,Zone2>0))))
    #Daily
BLOCKX.zn.daily=DATA.list.LIVEWT.c.daily
for ( i in 1:N.species)BLOCKX.zn.daily[[i]]=fn.blk.zone("BLOCKX",DATA.list.LIVEWT.c.daily[[i]])
BLOCKX.zn.daily=do.call(rbind,BLOCKX.zn.daily)
BLOCKX.zn.daily=as.data.frame.matrix(with(BLOCKX.zn.daily,table(BLOCKX,zone)))
BLKZ.per.zone.daily=list(West=as.numeric(rownames(subset(BLOCKX.zn.daily,West>0))),
                   Zone1=as.numeric(rownames(subset(BLOCKX.zn.daily,Zone1>0))),
                   Zone2=as.numeric(rownames(subset(BLOCKX.zn.daily,Zone2>0))))

BLKZ.per.zone$West=unique(c(BLKZ.per.zone$West,BLKZ.per.zone.daily$West))
BLKZ.per.zone$Zone1=unique(c(BLKZ.per.zone$Zone1,BLKZ.per.zone.daily$Zone1))
BLKZ.per.zone$Zone2=unique(c(BLKZ.per.zone$Zone2,BLKZ.per.zone.daily$Zone2))

    #Daily block 10
BLOCK_10.zn=DATA.list.LIVEWT.c.daily
for ( i in 1:N.species)BLOCK_10.zn[[i]]=fn.blk.zone("block10",DATA.list.LIVEWT.c.daily[[i]])
BLOCK_10.zn=do.call(rbind,BLOCK_10.zn)
BLOCK_10.zn=as.data.frame.matrix(with(BLOCK_10.zn,table(block10,zone)))
BLK_10.per.zone=list(West=as.numeric(rownames(subset(BLOCK_10.zn,West>0))),
                     Zone1=as.numeric(rownames(subset(BLOCK_10.zn,Zone1>0))),
                     Zone2=as.numeric(rownames(subset(BLOCK_10.zn,Zone2>0))))


      #4.22.6 Show effort creep       
Fish.pow.inc=c(1,(1+Fish.Pow*(1+(1975:1993)-1975)))
names(Fish.pow.inc)=1975:1994
Power.yrs=as.numeric(substr(FINYEAR.ALL,1,4))
Power.yrs=Power.yrs[which(!Power.yrs%in%names(Fish.pow.inc))]
Fish.pow.inc=c(Fish.pow.inc,rep(Fish.pow.inc[length(Fish.pow.inc)],length(Power.yrs)))
names(Fish.pow.inc)=FINYEAR.ALL


#Export default effort creep appliled
fn.fig("Effort_creep/Effort_creep_applied",2400,2400)
par(las=1,mgp=c(2.5,.9,0))
plot(1:length(Fish.pow.inc),Fish.pow.inc,xaxt='n',ylab="Effort creep",
     xlab="Finacial year",type='o',pch=19,cex.lab=2,cex.axis=1.5,cex=1.5)
axis(1,1:length(Fish.pow.inc),F,tck=-0.01)
axis(1,seq(1,length(Fish.pow.inc),10),FINYEAR.ALL[seq(1,length(Fish.pow.inc),10)],tck=-0.04,cex.axis=1.5)
segments(6,1.1,7,1.1,col=2,lwd=2)
segments(7,1.1,7,1.12,col=2,lwd=2)
text(9.5,1.1,"2 %",col=2,cex=2)
dev.off()



#Is there data in metro closures?
Metro.closure.data=vector('list',N.species)
names(Metro.closure.data)=SPECIES.vec
for(i in 1:N.species)Metro.closure.data[[i]]=fn.see.bl10.closure(DATA.list.LIVEWT.c.daily[[i]])


      #4.22.7 Standardise catches

#Check number of records by year for indicative vessels

  #Monthly
n.vesls.per.yer.mn=DATA.list.LIVEWT.c
SCaler=4
fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Vessel_pos_records_by_yr/Records_ind.vessels_monthly",2400,2400)  
par(mfcol=c(2,2),mar=c(1,.75,2.25,2),oma=c(3,2.5,.01,.1),las=1,mgp=c(3,.5,0),xpd=T)
for(n in 1:length(SPECIES.vec))
{
  DD=DATA.list.LIVEWT.c[[n]]
  DD=subset(DD,BLOCKX%in%as.numeric(BLKS.used[[n]]))      
  DD=subset(DD,VESSEL%in%VES.used[[n]])
  ves.yr=aggregate(Catch.Target~FINYEAR+VESSEL,subset(DD,Catch.Target>0),length)
  ves.yr$FIN=as.numeric(substr(ves.yr$FINYEAR,1,4))
  ves.yr$VESSEL=factor(ves.yr$VESSEL)
  ves.yr$VESSEL.level=match(ves.yr$VESSEL,levels(ves.yr$VESSEL))
  plot(ves.yr$FIN,ves.yr$VESSEL.level,pch=19,cex=(ves.yr$Catch.Target/max(ves.yr$Catch.Target))*SCaler,
       yaxt='n',ylab="",xlab="",col="steelblue")
  axis(2,1:length(levels(ves.yr$VESSEL)),levels(ves.yr$VESSEL),las=1,cex.axis=.8)
  mtext(paste(SPECIES.vec[n],"     "),3,0.25,cex=1.2)
  legend("topright",paste(c(5,10,20),"  "),pch=19,pt.cex=(c(5,10,20)/max(ves.yr$Catch.Target))*SCaler,
         bty='n',horiz=T,inset=-0.1,xjust=0,col="steelblue")
  
  ves.yr$single=with(ves.yr,ifelse(Catch.Target>0,1,0))
  n.vesls.per.yer.mn[[n]]=aggregate(single~FINYEAR,ves.yr,sum)
}
mtext(paste("Vessel","    "),2,1.25,las=3,outer=T,cex=1.35)
mtext("Financial year",1,1,outer=T,cex=1.35)
dev.off()

  #Daily
n.vesls.per.yer.daily=DATA.list.LIVEWT.c.daily
fn.fig("C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Vessel_pos_records_by_yr/Records_ind.vessels_daily",2400,2400)  
par(mfcol=c(2,2),mar=c(1,.75,2.25,2),oma=c(3,2.5,.01,.1),las=1,mgp=c(3,.5,0),xpd=T)
for(n in 1:length(SPECIES.vec))
{
  DD=DATA.list.LIVEWT.c.daily[[n]]
  DD=subset(DD,BLOCKX%in%as.numeric(BLKS.used.daily[[n]]))      
  DD=subset(DD,VESSEL%in%VES.used.daily[[n]])
  ves.yr=aggregate(Catch.Target~FINYEAR+VESSEL,subset(DD,Catch.Target>0),length)
  ves.yr$FIN=as.numeric(substr(ves.yr$FINYEAR,1,4))
  ves.yr$VESSEL=factor(ves.yr$VESSEL)
  ves.yr$VESSEL.level=match(ves.yr$VESSEL,levels(ves.yr$VESSEL))
  plot(ves.yr$FIN,ves.yr$VESSEL.level,pch=19,cex=(ves.yr$Catch.Target/max(ves.yr$Catch.Target))*SCaler,
       yaxt='n',ylab="",xlab="",col="steelblue")
  axis(2,1:length(levels(ves.yr$VESSEL)),levels(ves.yr$VESSEL),las=1,cex.axis=.8)
  mtext(paste(SPECIES.vec[n],"       "),3,0.25,cex=1.2)
  legend("topright",paste(c(10,25,50)," "),pch=19,pt.cex=(c(10,25,50)/max(ves.yr$Catch.Target))*SCaler,
         bty='n',horiz=T,inset=-0.1,xjust=0,col="steelblue")
  
  ves.yr$single=with(ves.yr,ifelse(Catch.Target>0,1,0))
  n.vesls.per.yer.daily[[n]]=aggregate(single~FINYEAR,ves.yr,sum)
}
mtext(paste("Vessel","    "),2,1.25,las=3,outer=T,cex=1.35)
mtext("Financial year",1,1,outer=T,cex=1.35)
dev.off()



#ACA
#apply best model to all scenarios   
Stand.out=vector('list',length=N.species)
names(Stand.out)=SPECIES.vec 
Stand.out.daily=Stand.out.zone=Stand.out.daily.zone=Stand.out

Stand.out.hour=Stand.out.daily.hour=Stand.out

system.time(for ( n in 1:N.species)    #<1 min to run with sensitivity tests  
{
  #Create storing lists
  if(do.sensitivity=="YES")
  {
    Zens=Tab.Sensi[-match("Effective",Tab.Sensi$Scenario),]
    Scen.store=vector('list',nrow(Zens))
    names(Scen.store)=Zens$Scenario
  } else 
  {
    Zens=Tab.Sensi[match("Stand_Base case",Tab.Sensi$Scenario),]
    Scen.store=vector('list',nrow(Zens))
    names(Scen.store)=Zens$Scenario
  }
  Scen.store.daily=Scen.store
  
  Scen.store.hour=Scen.store.daily.hour=Scen.store
  
  #zones combined
  for(i in 1:length(Scen.store))
  {
    ZEN=as.character(Zens$Scenario[i])
    Blocks_used=as.character(Zens[i,]$Blocks_used)
    Vessels_used=as.character(Zens[i,]$Vessels_used)
    
      #Monthly  
    DD=DATA.list.LIVEWT.c[[n]]
    if(Blocks_used=="Subset")  DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used[[n]],1,4)))      
    if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used[[n]] 
                     & FINYEAR %in% subset(n.vesls.per.yer.mn[[n]],single>=Threshold.n.vessls.per.yr)$FINYEAR) 
    if(ZEN%in%c("Stand_Efficiency","Stand_Area weighting")) Scen.store[[i]]=Scen.store[[match("Stand_Base case",names(Scen.store))]] else
    Scen.store[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model[[n]]$Bi,Formula.cpue.pos=Best.Model[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")

    Scen.store.hour[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model[[n]]$Bi,Formula.cpue.pos=Best.Model[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Hours")
    
      #Daily
    DD=DATA.list.LIVEWT.c.daily[[n]]
    if(Blocks_used=="Subset")  DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used.daily[[n]],1,4)))
    if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used.daily[[n]]) 
    if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used.daily[[n]] 
        & FINYEAR %in% subset(n.vesls.per.yer.daily[[n]],single>=Threshold.n.vessls.per.yr)$FINYEAR) 
    if(ZEN%in%c("Stand_Efficiency","Stand_Area weighting")) Scen.store.daily[[i]]=Scen.store.daily[[match("Stand_Base case",names(Scen.store.daily))]] else
    Scen.store.daily[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model.daily[[n]]$Bi,Formula.cpue.pos=Best.Model.daily[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")
 
    Scen.store.daily.hour[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model.daily[[n]]$Bi,Formula.cpue.pos=Best.Model.daily[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Hours")
    
    rm(DD)
  }
  Stand.out[[n]]=Scen.store
  Stand.out.daily[[n]]=Scen.store.daily

  Stand.out.hour[[n]]=Scen.store.hour
  Stand.out.daily.hour[[n]]=Scen.store.daily.hour
  
  
  #by zone
  dummy.zn=dummy.zn.daily=vector('list',length(BLKZ.per.zone))
  names(dummy.zn)=names(BLKZ.per.zone)
  if(!SPECIES.vec[n]%in%c("Gummy shark"))
  {
    for(z in 1:length(BLKZ.per.zone))
    {
      
      #Monthly
      DD=DATA.list.LIVEWT.c[[n]]
      DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used[[n]],1,4)))      
      DD=subset(DD,BLOCKX%in%BLKZ.per.zone[[z]])
      DD=subset(DD,VESSEL%in%VES.used[[n]]) 
      
      #remove years with very poor year coefficient estimation when subsetting by zone   
      if(ZONES[z]=="West" & SPECIES.vec[[n]]=="Whiskery shark") DD=subset(DD,!FINYEAR%in%c("1977-78","1978-79"))  
      if(SPECIES.vec[[n]]=="Dusky shark")
      {
        if(ZONES[z]=="West")  DD=subset(DD,!FINYEAR%in%c("1977-78"))
        if(ZONES[z]=="Zone1")  DD=subset(DD,!FINYEAR%in%c("1987-88"))
        if(ZONES[z]=="Zone2")  DD=subset(DD,!FINYEAR%in%c("1975-76"))
      }
      if(SPECIES.vec[[n]]=="Sandbar shark" & ZONES[z]=="Zone2")  DD=NULL
      
      if(!is.null(DD)) if(nrow(DD)>0)dummy.zn[[z]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model[[n]]$Bi,
                            Formula.cpue.pos=Best.Model[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")
      
      
      #Daily
      DD=DATA.list.LIVEWT.c.daily[[n]]
      DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used.daily[[n]],1,4)))
      DD=subset(DD,BLOCKX%in%BLKZ.per.zone[[z]])
      DD=subset(DD,VESSEL%in%VES.used.daily[[n]]) 
      if(SPECIES.vec[[n]]=="Sandbar shark" & ZONES[z]=="Zone2")  DD=NULL

      if(!is.null(DD)) if(nrow(DD)>0)dummy.zn.daily[[z]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model.daily[[n]]$Bi,
                          Formula.cpue.pos=Best.Model.daily[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")
    }
    Stand.out.zone[[n]]=dummy.zn
    Stand.out.daily.zone[[n]]=dummy.zn.daily
  }
})

#Set gummy by zone to Zone2 as only that zone falls within effective area
Stand.out.zone$`Gummy shark`=list(West=NULL,Zone1=NULL,Zone2=Stand.out$`Gummy shark`$`Stand_Base case`)
Stand.out.daily.zone$`Gummy shark`=list(West=NULL,Zone1=NULL,Zone2=Stand.out.daily$`Gummy shark`$`Stand_Base case`)


#Deviance and ANOVA Tables                
if(Extract.Deviance.table=="YES")
{
  #Get anovas
  Anova.Deviance=vector('list',length=N.species)
  names(Anova.Deviance)=SPECIES.vec 
  Anova.Deviance.daily=Anova.Deviance
  for ( i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    Anova.Deviance[[i]]=Anova.and.Dev.exp(Stand.out[[i]][[id.glm]]$GLMbi,Stand.out[[i]][[id.glm]]$GLMlog)
    Anova.Deviance.daily[[i]]=Anova.and.Dev.exp(Stand.out.daily[[i]][[id.glm]]$GLMbi,Stand.out.daily[[i]][[id.glm]]$GLMlog) 
  }
  
  #Output as nice table
  Tab.ano=vector('list',length=N.species)
  names(Tab.ano)=SPECIES.vec 
  Tab.ano.daily=Tab.ano
  for ( i in 1:N.species)
  {
    Tab.ano[[i]]=Comb.anova.and.dev(Anova.Deviance[[i]],SPECIES.vec[i])
    Tab.ano.daily[[i]]=Comb.anova.and.dev(Anova.Deviance.daily[[i]],SPECIES.vec[i])
  }
  
  #Monthly
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
  fn.word.table(WD=getwd(),TBL=do.call(rbind,Tab.ano),Doc.nm="ANOVA_table_monthly",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  #Daily
  fn.word.table(WD=getwd(),TBL=do.call(rbind,Tab.ano.daily),Doc.nm="ANOVA_table_daily",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
}



      #4.22.8 Get CI thru MC simulations   
Niters=1000
options(contrasts=c("contr.treatment","contr.poly")) #default
#set contrasts to sum to 0 rather than 'treatments' so the intercept is only the
# overall mean and doesn't have the first level of all factors
# options(contrasts=c("contr.sum","contr.poly"))

#Monthly
Index.monthly=Index.monthly.hour=Index.monthly.zone=DATA.list.LIVEWT.c
system.time(for ( i in 1:N.species)       #takes 0.2 sec per iteration
{
  #extract base case
  id.glm=match("Stand_Base case",names(Stand.out[[i]]))
  
  #Zones combined
  Index.monthly[[i]]=fn.MC.cpue(MOD=Stand.out[[i]][[id.glm]],
           Formula.cpue=Best.Model[[i]]$Bi,Formula.cpue.pos=Best.Model[[i]]$Log,
           niter=Niters,No.MC.Bi="YES",Area.w=AREA.W[[i]],BLKs="ALL",MESH="NO")   

  Index.monthly.hour[[i]]=fn.MC.cpue(MOD=Stand.out.hour[[i]][[id.glm]],
          Formula.cpue=Best.Model[[i]]$Bi,Formula.cpue.pos=Best.Model[[i]]$Log,
          niter=Niters,No.MC.Bi="YES",Area.w=AREA.W[[i]],BLKs="ALL",MESH="NO")   
  
    
  #By zone
  #Get blocks per zone     
  Store.zn=vector('list',length(ZONES))
  names(Store.zn)=ZONES
  if(fit.to.by.zone=="NO")
  {
    D=table(as.character(DATA.list.LIVEWT.c[[i]]$zone))
    for(g in 1:length(D))
    {
      BLks=unique(subset(DATA.list.LIVEWT.c[[i]],zone==names(D)[g],select="BLOCKX")$BLOCKX)
      Store.zn[[g]]=fn.MC.cpue(MOD=Stand.BaseCase[[i]],
                               Formula.cpue=Best.Model[[i]]$Bi,Formula.cpue.pos=Best.Model[[i]]$Log,
                               niter=Niters,No.MC.Bi="YES",Area.w=AREA.W[[i]],BLKs="byzone",MESH="NO")
    }
  }
  if(fit.to.by.zone=="YES")
  {
    #Run model over each zone
    for(g in 1:length(Store.zn))
    {
      if(!is.null(Stand.out.zone[[i]][[g]]))     #only do for zones with data
      {
        Store.zn[[g]]=fn.MC.cpue(MOD=Stand.out.zone[[i]][[g]],
              Formula.cpue=Best.Model[[i]]$Bi,Formula.cpue.pos=Best.Model[[i]]$Log,
              niter=Niters,No.MC.Bi="YES",Area.w=AREA.W[[i]],BLKs="ALL",MESH="NO")           
      }
    }
  }
  Index.monthly.zone[[i]]=Store.zn
  
  rm(Store.zn)
})

#Daily
Index.daily=Index.daily.hour=Index.daily.zone=DATA.list.LIVEWT.c.daily
system.time(for ( i in 1:N.species)       #takes 0.27 secs per iteration
{
  #extract base case
  id.glm=match("Stand_Base case",names(Stand.out[[i]]))
  
  #Zones combined
  Index.daily[[i]]=fn.MC.cpue(MOD=Stand.out.daily[[i]][[id.glm]],
       Formula.cpue=Best.Model.daily[[i]]$Bi,Formula.cpue.pos=Best.Model.daily[[i]]$Log,
       niter=Niters,No.MC.Bi="YES",Area.w=AREA.W_b10[[i]],BLKs="ALL",MESH="YES")    

  Index.daily.hour[[i]]=fn.MC.cpue(MOD=Stand.out.daily.hour[[i]][[id.glm]],
        Formula.cpue=Best.Model.daily[[i]]$Bi,Formula.cpue.pos=Best.Model.daily[[i]]$Log,
        niter=Niters,No.MC.Bi="YES",Area.w=AREA.W_b10[[i]],BLKs="ALL",MESH="YES")    
  
  #By zone
  #Get blocks per zone
  Store.zn=vector('list',length(ZONES))
  names(Store.zn)=ZONES
  if(fit.to.by.zone=="NO")
  {
    D=table(as.character(DATA.list.LIVEWT.c.daily[[i]]$zone))
    for(g in 1:length(D))
    {
      BLks=unique(subset(DATA.list.LIVEWT.c.daily[[i]],zone==names(D)[g],select="block10")$block10)
      Store.zn[[g]]=fn.MC.cpue(MOD=Stand.BaseCase.daily[[i]],
                               Formula.cpue=Best.Model.daily[[i]]$Bi,Formula.cpue.pos=Best.Model.daily[[i]]$Log,
                               niter=Niters,No.MC.Bi="YES",Area.w=AREA.W_b10[[i]],BLKs="byzone",MESH="YES")
    }
  }
  if(fit.to.by.zone=="YES")
  {
    #Run model over each zone
    for(g in 1:length(Store.zn))
    {
      if(!is.null(Stand.out.daily.zone[[i]][[g]]))    #only do for zones with data
      {
        Store.zn[[g]]=fn.MC.cpue(MOD=Stand.out.daily.zone[[i]][[g]],
             Formula.cpue=Best.Model.daily[[i]]$Bi,Formula.cpue.pos=Best.Model.daily[[i]]$Log,
             niter=Niters,No.MC.Bi="YES",Area.w=AREA.W_b10[[i]],BLKs="ALL",MESH="YES")
      }
    }      
  }
  Index.daily.zone[[i]]=Store.zn
  rm(Store.zn)
  
})


    #4.22.9 Base Case model fit diagnotics

#base case observed VS predicted      
setwd("C:/Matias/Analyses/Catch and effort/Outputs/model.fits")
if(Model.run=="First")
{
  for ( i in 1:N.species)     
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]])) 
    
    #Monthly
    fn.fig(paste(SPECIES.vec[i],".preds.vs.obs",sep=""),2400, 2400)
    see.BC.fit.to.data(GLM=Stand.out[[i]][[id.glm]],RespVar=Fit.to.what)     
    dev.off()
    
    #Daily
    if(Separate.monthly.daily=="YES")
    {
      fn.fig(paste(SPECIES.vec[i],".preds.vs.obs.daily",sep=""),2400, 2400)
      see.BC.fit.to.data(GLM=Stand.out.daily[[i]][[id.glm]],RespVar=Fit.to.what)
      dev.off()
    }
  
    
    #monthly
    fn.fig(paste(SPECIES.vec[i],"Fit.Diag",sep=""),2400,2400)
    fn.plot.diag(Stand.out[[i]][[id.glm]]$GLMlog,"LogNormal",SPECIES.vec[i])
    dev.off()
    
    #daily
    if(Separate.monthly.daily=="YES")
    {
      fn.fig(paste(SPECIES.vec[i],"Fit.Diag.daily",sep=""),2400,2400)
      fn.plot.diag(Stand.out.daily[[i]][[id.glm]]$GLMlog,"LogNormal",SPECIES.vec[i])
      dev.off()      
    }
  }
}

#Paper figure diagnostics
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
if(Model.run=="First")
{
  #Positive catch
  fn.fig("Appendix 6.monthly",2400, 2400)
  par(mfcol=c(3,4),las=1,mar=c(3,4,1.75,0.1),oma=c(1,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
  for ( i in 1:N.species) 
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    Pos.Diag.fn(Stand.out[[i]][[id.glm]]$GLMlog,SPECIES.vec[i]) 
  }
  dev.off()
  
  if(Separate.monthly.daily=="YES")
  {
    fn.fig("Appendix 6.daily",2400, 2400)
    par(mfcol=c(3,4),las=1,mar=c(3,4,1.75,0.1),oma=c(1,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
    for ( i in 1:N.species)
    {
      id.glm=match("Stand_Base case",names(Stand.out[[i]]))
      Pos.Diag.fn(Stand.out.daily[[i]][[id.glm]]$GLMlog,SPECIES.vec[i])
    }
    dev.off()
  }
  
}


    #4.22.10 Show Term Coefficients 
if(Model.run=="First")
{
  
  #--Positive catch
  
  #Paper figures
  
  #a. Month effect            
  WHERE=c("topleft","bottomleft","bottomleft","topright")
  fn.fig("Figure 2_Term_effect_Month",2400, 2400)
  par(mfcol=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    #Monthly
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    fn.show.coef(MOD=Stand.out[[i]][[id.glm]]$GLMlog,SP="",This="MONTH",ShowLevel="YES",WHRE=WHERE[i],CEX=1.1,Pt.cex=2)  
    if(i==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  }
    #Daily
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
    fn.show.coef(MOD=Stand.out.daily[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],
                 This="MONTH",ShowLevel="YES",WHRE=WHERE[i],CEX=1.25,Pt.cex=2)
    if(i==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  }
  mtext("Coefficient value",2,line=-0.8,outer=T,cex=1.75,las=3)
  mtext("Month",1,line=1.4,outer=T,cex=1.75)
  dev.off()
  
  
  #b. Block effect   
  numInt=50
  colfunc <- colorRampPalette(c("grey95","black"))
  couleurs=colfunc(numInt)
  fn.fig("Figure 3_Term_effect_BLOCKX",2000,2400)
  par(mfcol=c(4,2),mar=c(1,2,.5,.1),oma=c(2.5,2,1.25,.85),las=1,mgp=c(1.9,.5,0))
    #Monthly
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    BLk.dummy=as.numeric(substr(BLKS.used[[i]],1,4))
    fn.show.coef.spatial(MOD=Stand.out[[i]][[id.glm]]$GLMlog,SP="",This="BLOCKX",
          dat=subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0 & BLOCKX%in%BLk.dummy,select=c(LONG,LAT,BLOCKX))) 
    if(i%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    axis(2,seq(-36,-26,2),rev(seq(26,36,2)),tck=-0.025,cex.axis=1.35)
    if(i==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    axis(1,seq(113,129,2),F,tck=-0.025,cex.axis=1.35)
    axis(2,seq(-36,-26,2),F,tck=-0.025,cex.axis=1.35)  
  }
  color.legend(128,-35,129,-27,legend=c("Low",rep("",numInt-1),"High"),
               rect.col=couleurs,gradient="y",col=1,cex=0.9)  
    #Daily
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
    BLk.dummy=as.numeric(substr(BLKS.used.daily[[i]],1,4))
    fn.show.coef.spatial(MOD=Stand.out.daily[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],This="BLOCKX",
           dat=subset(DATA.list.LIVEWT.c.daily[[i]],Catch.Target>0 & BLOCKX%in%BLk.dummy,select=c(LONG,LAT,BLOCKX)))
    if(i%in%c(4))axis(1,seq(113,129,2),seq(113,129,2),tck=-0.025,cex.axis=1.35)
    axis(1,seq(113,129,2),F,tck=-0.025,cex.axis=1.35)
    axis(2,seq(-36,-26,2),F,tck=-0.025,cex.axis=1.35)
     if(i==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  }
  mtext("Latitude (ºS)",2,line=0,outer=T,cex=1.7,las=3)
  mtext("Longitude (ºE)",1,line=1.2,outer=T,cex=1.7,las=1)
  dev.off()
  

  #individual figures
  hnDl.term="C:/Matias/Analyses/Catch and effort/Outputs/Term effect/"
  
    #1. Monthly data 
  #Month effect
  WHERE=c("topleft","bottomleft","bottomleft","topright")
  fn.fig(paste(hnDl.term,"Term_effect_Month_Monthly",sep=""),2400, 2400)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(4,4,.1,0.1),las=1,mgp=c(1.9,.6,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    fn.show.coef(MOD=Stand.out[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],This="MONTH",ShowLevel="YES",WHRE=WHERE[i],CEX=1.1,Pt.cex=2)
  }
  mtext("Coefficient value",2,line=1.75,outer=T,cex=2,las=3)
  mtext("Month",1,line=2,outer=T,cex=2)
  dev.off()
  
  #BLOCKx effect
  colfunc <- colorRampPalette(c("yellow","red"))
  couleurs=colfunc(numInt)  
  fn.fig(paste(hnDl.term,"Term_effect_BLOCKX_Spatial",sep=""),2400, 1800)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(3,3,.1,0.1),las=2,mgp=c(1.9,.5,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    BLk.dummy=as.numeric(substr(BLKS.used[[i]],1,4))
    fn.show.coef.spatial(MOD=Stand.out[[i]][[id.glm]]$GLMlog,SP="",This="BLOCKX",
       dat=subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0 & BLOCKX%in%BLk.dummy,select=c(LONG,LAT,BLOCKX))) 
     if(i%in%c(2,4))axis(1,113:129,113:129,tck=-0.015)
    if(i%in%c(1:2))axis(2,-36:-26,rev(26:36),tck=-0.015)
  }
  mtext("Latitude (ºS)",2,line=1.4,outer=T,cex=1.75,las=3)
  mtext("Longitude (ºE)",1,line=1.65,outer=T,cex=1.75,las=1)
  color.legend(128,-35,129,-27,legend=c("Low",rep("",numInt-1),"High"),
               rect.col=couleurs,gradient="y",col=1,cex=0.9)  
  dev.off()
  
  #Vessel effect
  WHERE=c("topleft","topleft","bottomleft","bottomleft")
  fn.fig(paste(hnDl.term,"Term_effect_VESSEL_Monthly",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(3,4,.1,0.1),las=2,mgp=c(1.9,.7,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out[[i]]))
    fn.show.coef(MOD=Stand.out[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],This="VESSEL",ShowLevel="NO",WHRE=WHERE[i],CEX=1,Pt.cex=1)
  }
  mtext("Coefficient value",2,line=1.75,outer=T,cex=2,las=3)
  mtext("Vessel",1,line=1.5,outer=T,cex=2,las=1)
  dev.off()
  
      #2. Daily data 
  #Month effect
  WHERE=c("topleft","bottomleft","bottomleft","topright")
  fn.fig(paste(hnDl.term,"Term_effect_Month_Daily",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(4,4,.1,0.1),las=1,mgp=c(1.9,.7,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
    fn.show.coef(MOD=Stand.out.daily[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],This="MONTH",ShowLevel="YES",WHRE=WHERE[i],CEX=1.25,Pt.cex=2)
  }
  mtext("Coefficient value",2,line=1.75,outer=T,cex=2,las=3)
  mtext("Month",1,line=2,outer=T,cex=2)
  dev.off()
  
  #Block effect
  fn.fig(paste(hnDl.term,"Term_effect_BLOCKX_Spatial_daily",sep=""),2400,1800)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(3,3,.1,0.1),las=2,mgp=c(1.9,.5,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
    BLk.dummy=as.numeric(substr(BLKS.used.daily[[i]],1,4))
    fn.show.coef.spatial(MOD=Stand.out.daily[[i]][[id.glm]]$GLMlog,SP="",This="BLOCKX",
         dat=subset(DATA.list.LIVEWT.c.daily[[i]],Catch.Target>0 & BLOCKX%in%BLk.dummy,select=c(LONG,LAT,BLOCKX))) 
    if(i%in%c(2,4))axis(1,113:129,113:129,tck=-0.015)
    if(i%in%c(1:2))axis(2,-36:-26,rev(26:36),tck=-0.015)
  }
  mtext("Latitude (ºS)",2,line=1.4,outer=T,cex=1.75,las=3)
  mtext("Longitude (ºE)",1,line=1.65,outer=T,cex=1.75,las=1)
  color.legend(128,-35,129,-27,legend=c("Low",rep("",numInt-1),"High"),
               rect.col=couleurs,gradient="y",col=1,cex=0.9)  
  dev.off()
  
  #Vessel effect
  WHERE=c("topleft","topleft","bottomleft","bottomleft")
  fn.fig(paste(hnDl.term,"Term_effect_VESSEL_Daily",sep=""),2400,2400)
  par(mfcol=c(2,2),mar=c(1,1,1,1.1),oma=c(3,4,.1,0.1),las=2,mgp=c(1.9,.7,0))
  for(i in 1:N.species)
  {
    id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
    fn.show.coef(MOD=Stand.out.daily[[i]][[id.glm]]$GLMlog,SP=SPECIES.vec[i],This="VESSEL",ShowLevel="NO",WHRE=WHERE[i],CEX=1,Pt.cex=1)
  }
  mtext("Coefficient value",2,line=1.75,outer=T,cex=2,las=3)
  mtext("Vessel",1,line=1.5,outer=T,cex=2,las=1)
  dev.off()
}

    #4.22.11 Records used
hndl.kept="C:/Matias/Analyses/Catch and effort/Outputs/Kept_blocks_vessels/"
LsT.rec.blk=data.frame(Species=SPECIES.vec,Monthly=NA,Daily=NA)
LsT.blok.blk=LsT.ves.blk=LsT.rec.blk
Y=-36:-27; X=seq(113,129,length.out=length(Y))
for ( i in 1:N.species) 
{
  #Monthly
  id.glm=match("Stand_Base case",names(Stand.out[[i]]))
  DAT=Stand.out[[i]][[id.glm]]$BiData
  BLk.uni=unique(DAT$BLOCKX)
  VLS.uni=unique(DAT$VESSEL)
  LsT.rec.blk$Monthly[i]=1-nrow(DAT)/nrow(DATA.list.LIVEWT.c[[i]])
  LsT.blok.blk$Monthly[i]=1-length(BLk.uni)/length(unique(DATA.list.LIVEWT.c[[i]]$BLOCKX))
  LsT.ves.blk$Monthly[i]=1-length(VLS.uni)/length(unique(DATA.list.LIVEWT.c[[i]]$VESSEL))
  write.csv(BLk.uni,paste(hndl.kept,"blocks_used_",SPECIES.vec[i],"_monthly.csv",sep=""))
  fn.fig(paste(hndl.kept,"block_used_map_",SPECIES.vec[i],"_monthly",sep=""),2400,2400)    
  par(las=1,mai=c(.7,.7,.1,.1),mgp=c(2,.7,0))
  fn.show.blk(BLk.uni,BLk.uni,show.dropped="NO")
  mtext(SPECIES.vec[i],3,-3,cex=2.5)
  dev.off()  
  
    #Daily
  id.glm=match("Stand_Base case",names(Stand.out.daily[[i]]))
  DAT=Stand.out.daily[[i]][[id.glm]]$BiData
  BLk.uni=unique(DAT$BLOCKX)
  VLS.uni=unique(DAT$VESSEL)
  LsT.rec.blk$Daily[i]=1-nrow(DAT)/nrow(DATA.list.LIVEWT.c.daily[[i]])
  LsT.blok.blk$Daily[i]=1-length(BLk.uni)/length(unique(DATA.list.LIVEWT.c.daily[[i]]$BLOCKX))
  LsT.ves.blk$Daily[i]=1-length(VLS.uni)/length(unique(DATA.list.LIVEWT.c.daily[[i]]$VESSEL))
  write.csv(unique(DAT$BLOCKX),paste(hndl.kept,"blocks_used_",SPECIES.vec[i],"_daily.csv",sep=""))
  
  fn.fig(paste(hndl.kept,"block_used_map_",SPECIES.vec[i],"_daily",sep=""),2400,2400)    
  par(las=1,mai=c(.7,.7,.1,.1),mgp=c(2,.7,0))
  fn.show.blk(BLk.uni,BLk.uni,show.dropped="NO")
  mtext(SPECIES.vec[i],3,-3,cex=2.5)
  dev.off()  
}
write.csv(LsT.rec.blk,paste(hndl.kept,"Proportion.lost.records.csv",sep=""))
write.csv(LsT.blok.blk,paste(hndl.kept,"Proportion.lost.blocks.csv",sep=""))
write.csv(LsT.ves.blk,paste(hndl.kept,"Proportion.lost.vessels.csv",sep=""))

Prop.catch.explained=DATA.list.LIVEWT.c
fn.fig("Appendix 5. Vessels and Catch by year",2400,2400) 
par(mfrow=c(4,2),mai=c(.45,.5,.1,.2),oma=c(.4,.1,.1,.95),mgp=c(1.5,.6,0),las=1)
for(i in 1:N.species) 
{
  xx=DATA.list.LIVEWT.c[[i]]
  yy=subset(DATA.list.LIVEWT.c[[i]],VESSEL%in%VES.used[[i]])
  # yy=subset(DATA.list.LIVEWT.c[[i]],BLOCKX%in%as.numeric(substr(BLKS.used[[i]],1,4)) &
  #             VESSEL%in%VES.used[[i]])
  Prop.catch.explained[[i]]=fn.plot.vess.selection(xx,yy)  
  mtext(SPECIES.vec[i],4,las=3,line=0.75,cex=1.25)
}
mtext("Financial year",1,outer=T,line=-.75,cex=1.5)
dev.off()
write.csv(unlist(Prop.catch.explained),"Prop.of.catch.explained.after.dropping.vessels.csv")


    #4.22.12 Plot data gaps 

#Display proportion with and without catch for each species
MAX.list=c(45,100,30,15)
Ymax.vec=c(.24,.24,.24,.24)
Xmax.vec=c(11,11,11,11)
DATA.list.0.catch.2=vector('list',length=N.species)
names(DATA.list.0.catch.2)=SPECIES.vec

  # Daily combined with monthly 
if(Model.run=="First")
{
  if(Separate.monthly.daily=="NO")
  {
    for ( i in 1:N.species)DATA.list.0.catch.2[[i]]=Zero.catch.fun2(DATA.list.LIVEWT.c[[i]])
    
    fn.fig("Appendix 4.Prop of records with 0 catch",2400,2400)    
    par(mfcol=c(2,2),mar=c(1,3.6,1,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
    
    for ( i in 1:N.species) plot0(DATA.list.LIVEWT.c[[i]],MAX.list[i],Ymax.vec[i],Xmax.vec[i],SPECIES.vec[i])
    mtext("Catch per record (tonnes)",side=1,line=1,font=1,las=0,cex=1.25,outer=T)
    mtext("Proportion of monthly records with catch",side=2,line=-1.2,font=1,las=0,cex=1.25,outer=T)
    
    pos.list=list(c(.1,.5,.75,1),c(.1,.5,.25,.5),c(.6,1,.75,1),c(.6,1,.25,.5))
    for ( i in 1:N.species)
    {
      par(fig=pos.list[[i]], new = T,mgp=c(.1,.4,0))
      LABS=DATA.list.0.catch.2[[i]]$FINYEAR
      N.yrs=nrow(DATA.list.0.catch.2[[i]])
      ZeroKTC=DATA.list.0.catch.2[[i]]$proportion
      plot(1:N.yrs,ZeroKTC,ylab="",xlab="",xaxt="n",las=2,pch=19,cex=1,ylim=c(0,1)
           ,cex.axis=1.25,yaxt="n",lwd=1)
      axis(1,at=1:N.yrs,labels=F,tck=-0.02)
      axis(1,at=seq(1,N.yrs,by=5),labels=LABS[seq(1,N.yrs,by=5)],tck=-0.04,cex.axis=.9)
      axis(2,at=seq(0,1,.2),labels=seq(0,1,.2),tck=-0.04,cex.axis=.9)
      legend("bottomright",SPECIES.vec[i],bty='n',cex=1.25)
      mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1,outer=F)
      #mtext("Proportion of records with catch",side=2,line=-1.2,font=1,las=0,cex=.75,outer=F)
    }
    dev.off()
    
  }
  
  #daily separated from monthly
  if(Separate.monthly.daily=="YES")
  {
    for ( i in 1:N.species)
    {
      dummy=rbind(subset(DATA.list.LIVEWT.c[[i]],select=c(FINYEAR,Catch.Target)),
                  subset(DATA.list.LIVEWT.c.daily[[i]],select=c(FINYEAR,Catch.Target)))
      dummy=dummy[order(dummy$FINYEAR),]
      DATA.list.0.catch.2[[i]]=Zero.catch.fun2(dummy)
    }
    
    fn.fig("Appendix 4.Prop of records with 0 catch",2400, 2400)    
    par(mfcol=c(2,2),mar=c(1.5,1,1,1.2),oma=c(2.5,3,.1,.1),las=1,mgp=c(.1,.5,0))
    for ( i in 1:N.species)
    {
      ss=DATA.list.0.catch.2[[i]]
      if(SPECIES.vec[i]=="Sandbar shark") ss=subset(ss,FINYEAR%in%c(San.Yrs,FINYEAR.daily))
      ZeroKTC=1-ss$proportion
      LABS=as.character(ss$FINYEAR)
      NN=length(LABS)
      plot(1:NN,ZeroKTC,ylab="",xlab="",xaxt="n",las=2,pch=19,cex=1,ylim=c(0,1),cex.axis=1.5,yaxt="n")
      polygon(x=c(1,match("2005-06",LABS),match("2005-06",LABS),1),
              y=c(-0.5,-0.5,1.1,1.1),col="grey75",border="grey75")
      polygon(x=c(match("2005-06",LABS),NN,NN,match("2005-06",LABS)),
              y=c(-0.5,-0.5,1.1,1.1),col="grey95",border="grey95")
      points(1:NN,ZeroKTC,las=2,pch=19,cex=1.25)
      box()
      axis(1,at=1:NN,labels=F,tck=-0.015)
      axis(1,at=seq(1,NN,by=5),labels=LABS[seq(1,NN,by=5)],tck=-0.025,cex.axis=1.15)
      axis(2,at=seq(0,1,.2),labels=seq(0,1,.2),tck=-0.025,cex.axis=1.15)
      legend("topleft",SPECIES.vec[i],bty='n',cex=1.5)
      
    }
    mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.5,outer=T)
    mtext("Proportion of records with no catch",side=2,line=1.5,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
  }
  
}

#Plot missing year-blocks for each species
if(Model.run=="First")
{
  Missing.yr.blk=DATA.list.LIVEWT.c
  for (i in 1:N.species)
  {
    fn.fig(paste("Complex/Figure 6",names(Species.list)[i],".Block.byYR"),2400, 2400)    
    dull=DATA.list.LIVEWT.c[[i]]
    if(Separate.monthly.daily=="NO")
    {
      dummy=subset(dull,select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,BLOCKX))
    }
    
    if(Separate.monthly.daily=="YES")
    {
      dummy=rbind(subset(dull,select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,BLOCKX)),
                  subset(DATA.list.LIVEWT.c.daily[[i]],select=c(Catch.Target,Km.Gillnet.Days.c,FINYEAR,BLOCKX)))  
    }
    Missing.yr.blk[[i]]=plt.miss.yr_blk(dummy)
    dev.off()
  }
}


      #4.22.13 Calculate relative standardised catch for MC simulations 
#steps: 1. normalise is applicable

#create objet for storing relative Indices
Index.monthly.rel=Index.monthly    
Index.monthly.zone.rel=Index.monthly.zone
if(Q_change=="1982-83")
{
  Index.monthly_q.change.rel=Index.monthly_q.change
  Index.monthly.zone_q.change.rel=Index.monthly.zone_q.change
}
Index.daily.rel=Index.daily
Index.daily.zone.rel=Index.daily.zone

Index.monthly.hour.rel=Index.monthly.hour 
Index.daily.hour.rel=Index.daily.hour

List.foly.nom.rel=List.foly.nom    

for(i in 1:N.species)
{
  #zones combined
  Index.monthly.rel[[i]]=fn.relative(D=Index.monthly.rel[[i]])
  if(Q_change=="1982-83")if(SPECIES.vec[[i]]=="Whiskery shark")Index.monthly_q.change.rel[[i]]=fn.relative(D=Index.monthly_q.change.rel[[i]])
  Index.daily.rel[[i]]=fn.relative(D=Index.daily.rel[[i]])
  
  Index.monthly.hour.rel[[i]]=fn.relative(D=Index.monthly.hour.rel[[i]])
  Index.daily.hour.rel[[i]]=fn.relative(D=Index.daily.hour.rel[[i]])
  
  #By zone
  Zn=names(Index.monthly.zone[[i]])
  for(z in 1:length(Zn)) 
  {
    if(!is.null(Index.monthly.zone.rel[[i]][[z]])) 
    {
      Index.monthly.zone.rel[[i]][[z]]=fn.relative(D=Index.monthly.zone.rel[[i]][[z]])
      if(Q_change=="1982-83")if(SPECIES.vec[[i]]=="Whiskery shark")Index.monthly.zone_q.change.rel[[i]][[z]]=fn.relative(D=Index.monthly.zone_q.change.rel[[i]][[z]])
    }
  }
  
  Zn=names(Index.daily.zone[[i]])
  for(z in 1:length(Zn))
  {
    if(!is.null(Index.daily.zone.rel[[i]][[z]])) Index.daily.zone.rel[[i]][[z]]=fn.relative(D=Index.daily.zone.rel[[i]][[z]])
  }
  
  
  #Unstandardised   
  #note: Foly already has efficiency creep built in so remove it first, calculate relative, and then apply efficiency
  List.foly.nom.rel[[i]]$Foly=fn.relative.fol.nom(D=List.foly.nom[[i]]$Foly,what="Folly")
  List.foly.nom.rel[[i]]$Nom=fn.relative.fol.nom(D=List.foly.nom[[i]]$Nom,what="Nom") 
  
}


#     4.22.14 Calculate relative standardised catch for sensitivity tests
#steps: 1. predict year-block
#       2. apply block weight is applicable
#       3. aggregate by year
#       4. normalise is applicable
#       5. remove efficiency is applicable

if(do.sensitivity=="YES")       
{
  Sens.list.monthly=Stand.out
  Sens.list.daily=Stand.out.daily
  
  #1. Predict area-weighted annual index (steps 1-3)
  for ( n in 1:N.species)
  {
    dummy=vector('list',length(Scen.store))
    names(dummy)=names(Scen.store)
    dummy.daily=dummy
    for(i in 1:length(Scen.store))
    {
      id=match(names(Stand.out[[n]])[i],Tab.Sensi$Scenario)
      Do.area=as.character(Tab.Sensi$Area_weighting[id])

        #Monthly
      dummy[[i]]=fn.stand.cpue.sens(MOD=Stand.out[[n]][[i]],FORMULA=Best.Model[[n]],Apply.area.w=Do.area,Ar=AREA.W[[n]])
      
        #Daily
      dummy.daily[[i]]=fn.stand.cpue.sens(MOD=Stand.out.daily[[n]][[i]],FORMULA=Best.Model.daily[[n]],Apply.area.w=Do.area,Ar=AREA.W[[n]])
      
    }
    Sens.list.monthly[[n]]=dummy
    Sens.list.daily[[n]]=dummy.daily
  }
  
  #2. Normalise
  Sens.list.monthly.rel=Sens.list.monthly  
  Sens.list.daily.rel=Sens.list.daily
  for ( n in 1:N.species)
  {
    for(i in 1:length(Scen.store))
    {
        #Monthly
      Sens.list.monthly.rel[[n]][[i]]$Index.w=Sens.list.monthly.rel[[n]][[i]]$Index.w/mean(Sens.list.monthly.rel[[n]][[i]]$Index.w)
      Sens.list.monthly.rel[[n]][[i]]$cpue=Sens.list.monthly.rel[[n]][[i]]$cpue/mean(Sens.list.monthly.rel[[n]][[i]]$cpue)
      
        #Daily
      Sens.list.daily.rel[[n]][[i]]$Index.w=Sens.list.daily.rel[[n]][[i]]$Index.w/mean(Sens.list.daily.rel[[n]][[i]]$Index.w)
      Sens.list.daily.rel[[n]][[i]]$cpue=Sens.list.daily.rel[[n]][[i]]$cpue/mean(Sens.list.daily.rel[[n]][[i]]$cpue)
    }
  }
  
  #3. Remove efficiency         #ACA: this is not required as non creep effort is used in glm
  for ( n in 1:N.species)
  {
    for(i in 1:length(Scen.store))
    {
      if(!names(Stand.out[[n]])[i]=="Stand_Efficiency")  
      {
        #absolute
        Sens.list.monthly[[n]][[i]]=fn.remove.eff.creep(d=Sens.list.monthly[[n]][[i]],VRS=c("Index.w","cpue"))
        Sens.list.daily[[n]][[i]]=fn.remove.eff.creep(d=Sens.list.daily[[n]][[i]],VRS=c("Index.w","cpue"))
        
        #relative
        Sens.list.monthly.rel[[n]][[i]]=fn.remove.eff.creep(d=Sens.list.monthly.rel[[n]][[i]],VRS=c("Index.w","cpue"))
        Sens.list.daily.rel[[n]][[i]]=fn.remove.eff.creep(d=Sens.list.daily.rel[[n]][[i]],VRS=c("Index.w","cpue"))
      }
    }  
  }
}

#     4.22.15 Create index without effort creep
Index.monthly_no.creep=Index.monthly
Index.monthly.zone_no.creep=Index.monthly.zone
Index.monthly.rel_no.creep=Index.monthly.rel
Index.monthly.zone.rel_no.creep=Index.monthly.zone.rel

Index.daily_no.creep=Index.daily 
Index.daily.zone_no.creep=Index.daily.zone
Index.daily.rel_no.creep=Index.daily.rel  
Index.daily.zone.rel_no.creep=Index.daily.zone.rel


Index.monthly.hour_no.creep=Index.monthly.hour
Index.monthly.hour.rel_no.creep=Index.monthly.hour.rel
Index.daily.hour_no.creep=Index.daily.hour 
Index.daily.hour.rel_no.creep=Index.daily.hour.rel  

for ( n in 1:N.species)
{
    #Monthly
    Index.monthly_no.creep[[n]]=Index.monthly[[n]]
    Index.monthly.rel_no.creep[[n]]=Index.monthly.rel[[n]]
    Index.monthly.hour_no.creep[[n]]=Index.monthly.hour[[n]]
    Index.monthly.hour.rel_no.creep[[n]]=Index.monthly.hour.rel[[n]]
    
      #remove effort creep   #ACA: this is not required as non.creep effort is used in glm
    Index.monthly[[n]]=fn.remove.eff.creep(d=Index.monthly[[n]],VRS=c("MEAN","SD","LOW","UP"))   
    Index.monthly.rel[[n]]=fn.remove.eff.creep(d=Index.monthly.rel[[n]],VRS=c("MEAN","SD","LOW","UP"))
    Index.monthly.hour[[n]]=fn.remove.eff.creep(d=Index.monthly.hour[[n]],VRS=c("MEAN","SD","LOW","UP"))
    Index.monthly.hour.rel[[n]]=fn.remove.eff.creep(d=Index.monthly.hour.rel[[n]],VRS=c("MEAN","SD","LOW","UP"))
    
    Zn=names(Index.monthly.zone[[n]])
    for(z in 1:length(Zn))
    {
      if(!is.null(Index.monthly.zone[[n]][[z]]))
      {
        Index.monthly.zone_no.creep[[n]][[z]]=Index.monthly.zone[[n]][[z]]
        Index.monthly.zone.rel_no.creep[[n]][[z]]=Index.monthly.zone.rel[[n]][[z]]
        
          #remove effort creep       #ACA: this is not required as non.creep effort is used in glm
        Index.monthly.zone[[n]][[z]]=fn.remove.eff.creep(d=Index.monthly.zone[[n]][[z]],VRS=c("MEAN","SD","LOW","UP"))
        Index.monthly.zone.rel[[n]][[z]]=fn.remove.eff.creep(d=Index.monthly.zone.rel[[n]][[z]],VRS=c("MEAN","SD","LOW","UP"))
        
      }
    }
    
    #Daily
    Index.daily_no.creep[[n]]=Index.daily[[n]]
    Index.daily.rel_no.creep[[n]]=Index.daily.rel[[n]]
    Index.daily.hour_no.creep[[n]]=Index.daily.hour[[n]]
    Index.daily.hour.rel_no.creep[[n]]=Index.daily.hour.rel[[n]]

        #remove effort creep       #ACA: this is not required as non.creep effort is used in glm
    Index.daily[[n]]=fn.remove.eff.creep(d=Index.daily[[n]],VRS=c("MEAN","SD","LOW","UP"))
    Index.daily.rel[[n]]=fn.remove.eff.creep(d=Index.daily.rel[[n]],VRS=c("MEAN","SD","LOW","UP"))
    Index.daily.hour[[n]]=fn.remove.eff.creep(d=Index.daily.hour[[n]],VRS=c("MEAN","SD","LOW","UP"))
    Index.daily.hour.rel[[n]]=fn.remove.eff.creep(d=Index.daily.hour.rel[[n]],VRS=c("MEAN","SD","LOW","UP"))
    
    for(z in 1:length(Zn))
    {
      if(!is.null(Index.monthly.zone[[n]][[z]]))
      {
        Index.daily.zone_no.creep[[n]][[z]]=Index.daily.zone[[n]][[z]]
        Index.daily.zone.rel_no.creep[[n]][[z]]=Index.daily.zone.rel[[n]][[z]]
        
          #remove effort creep     #ACA: this is not required as non.creep effort is used in glm
        Index.daily.zone[[n]][[z]]=fn.remove.eff.creep(d=Index.daily.zone[[n]][[z]],VRS=c("MEAN","SD","LOW","UP"))
        Index.daily.zone.rel[[n]][[z]]=fn.remove.eff.creep(d=Index.daily.zone.rel[[n]][[z]],VRS=c("MEAN","SD","LOW","UP"))
      }
    }
}


#     4.22.16 Plot indices           

#4.22.16.1 Base case and effective
LEGs=c("Effective","Nominal","Standardised")

  #--RElative index
#a. Zones combined     

  #a.1 Species together          
fn.fig("Figure 6.Annual_Index",2000, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for ( i in 1:N.species)
{
  #Monthly
  if(!Q_change=="1982-83")PLOT.Index(Index=Index.monthly.rel[[i]],ADD.his.nom="YES",
             FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=NULL)  
  if(Q_change=="1982-83")PLOT.Index(Index=Index.monthly.rel[[i]],ADD.his.nom="YES",
            FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=Index.monthly_q.change.rel[[i]])  
  if(i==1)
  {
    mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",LEGs,pch=c(21,21,NA),bty='n',cex=1.5,pt.cex=1.5,
           pt.bg=c("grey50","white",NA),lty=c(0,0,1),lwd=c(1,1,2))
  }
  
  #Daily
  PLOT.Index(Index=Index.daily.rel[[i]],ADD.his.nom="YES",
             FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=NULL)
  if(i==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",SPECIES.vec[i],bty='n',cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()

  #a.2 By Species           
for ( i in 1:N.species)
{
  fn.fig(paste("Annual_Index_relative_",SPECIES.vec[i],sep=""),2400,1600)    
  par(mfrow=c(1,2),mar=c(1,3,1.5,.95),oma=c(2.5,.8,.1,.1),las=1,mgp=c(1.9,.7,0))
  
  #Monthly
  if(!Q_change=="1982-83")PLOT.Index(Index=Index.monthly.rel[[i]],ADD.his.nom="YES",
                                     FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=NULL)
  if(Q_change=="1982-83")PLOT.Index(Index=Index.monthly.rel[[i]],ADD.his.nom="YES",
                                    FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=Index.monthly_q.change.rel[[i]]) 
  mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",LEGs,pch=c(21,21,NA),bty='n',cex=1.5,pt.cex=1.5,
         pt.bg=c("grey50","white",NA),lty=c(0,0,1),lwd=c(1,1,2))
  
  
  #Daily
  PLOT.Index(Index=Index.daily.rel[[i]],ADD.his.nom="YES",
             FOLY=List.foly.nom.rel[[i]]$Foly,NOM=List.foly.nom.rel[[i]]$Nom,Index_q.change=NULL)
  mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
}


#b. By zone                                   
CL=c("black","grey60","white")
names(CL)=c("West","Zone1","Zone2")  

  #b.1 Species together
Where.Leg=c("topleft","top","top","topleft")
fn.fig("Figure 6.Annual_Index_by.zone",2400, 2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.1),las=1,mgp=c(1.9,.7,0))
for ( i in 1:N.species)
{
  #Monthly
  YRss=FINYEAR.monthly
  if( SPECIES.vec[i]=="Sandbar shark") YRss=San.Yrs
  if(!Q_change=="1982-83")PLOT.Index.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=NULL,ALL.yrs=YRss)
  if(Q_change=="1982-83")PLOT.Index.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=Index.monthly.zone_q.change.rel[[i]],
                                         ALL.yrs=YRss)
  if(i==1)
  {
    mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",c(" West","Zone 1","Zone 2"),bty='n',pt.cex=2.5,cex=1.5,horiz=T,
           pch=rep(21,3),col=rep(1,3),pt.bg=CL)      
  }
  
  #Daily
  PLOT.Index.zone(LISTA=Index.daily.zone.rel[[i]],LISTA2=NULL,ALL.yrs=FINYEAR.daily)
  legend(Where.Leg[i],SPECIES.vec[i],bty='n',cex=1.8)
  if(i==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()

  #b.2 By species
Where.Leg=c("topleft","top","top","topleft")  
for ( i in 1:N.species)
{
  fn.fig(paste("Annual_Index_relative_zone_",SPECIES.vec[i],sep=""),2400,1400)    
  par(mfrow=c(1,2),mar=c(1,3,1.5,.95),oma=c(2.5,.8,.1,.1),las=1,mgp=c(1.9,.7,0))
  
  #Monthly
  if(!Q_change=="1982-83")PLOT.Index.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=NULL,ALL.yrs=FINYEAR.monthly)
  if(Q_change=="1982-83")PLOT.Index.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=Index.monthly.zone_q.change.rel[[i]],
                                         ALL.yrs=FINYEAR.monthly)
  legend("topright",c(" West","Zone 1","Zone 2"),bty='n',pt.cex=2,cex=1,horiz=T,
         pch=rep(21,3),col=rep(1,3),pt.bg=CL)      
  mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  PLOT.Index.zone(LISTA=Index.daily.zone.rel[[i]],LISTA2=NULL,ALL.yrs=FINYEAR.daily)
  legend(Where.Leg[i],"",bty='n',cex=1.8)    
  mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  
  dev.off()
}


#Temporal trends in CVs by zone (base case by zone)
CL=c("black","grey60","white")
names(CL)=c("West","Zone1","Zone2")  
Where.Leg=c("right","bottomleft","bottomleft","bottomleft")

  #species combined
fn.fig("Appendix 7.CV_by.zone",2400,2400)    
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.1),las=1,mgp=c(1.9,.7,0))
for ( i in 1:N.species)
{
  #Monthly
  YYRS=FINYEAR.monthly
  if(SPECIES.vec[i]=="Sandbar shark")YYRS=San.Yrs
  if(!Q_change=="1982-83")PLOT.CV.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=NULL,ALL.yrs=YYRS)
  if(Q_change=="1982-83")PLOT.CV.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=Index.monthly.zone_q.change.rel[[i]],
                                      ALL.yrs=YYRS)
  if(i==1)
  {
    mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    legend("toprigh",c(" West","Zone 1","Zone 2"),bty='n',pt.cex=2.5,cex=1.5,horiz=T,
           pch=rep(21,3),col=rep(1,3),pt.bg=CL)      
  }
  
  #Daily
  PLOT.CV.zone(LISTA=Index.daily.zone.rel[[i]],LISTA2=NULL,ALL.yrs=FINYEAR.daily)
  legend(Where.Leg[i],SPECIES.vec[i],bty='n',cex=1.8)
  if(i==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("CV",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()

  #by species
for ( i in 1:N.species)
{
  fn.fig(paste("Annual_Index_relative_zone_CV_",SPECIES.vec[i],sep=""),2400,1400)    
  par(mfrow=c(1,2),mar=c(1,3,1.5,.95),oma=c(2.5,.8,.1,.1),las=1,mgp=c(1.9,.7,0))
  
  #Monthly
  YYRS=FINYEAR.monthly
  if(SPECIES.vec[i]=="Sandbar shark")YYRS=San.Yrs
  if(!Q_change=="1982-83")PLOT.CV.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=NULL,ALL.yrs=YYRS)
  if(Q_change=="1982-83")PLOT.CV.zone(LISTA=Index.monthly.zone.rel[[i]],LISTA2=Index.monthly.zone_q.change.rel[[i]],
                                      ALL.yrs=YYRS)
  mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  legend("bottomleft",c(" West","Zone 1","Zone 2"),bty='n',pt.cex=2,cex=1,horiz=T,
         pch=rep(21,3),col=rep(1,3),pt.bg=CL) 
  
  #Daily
  PLOT.CV.zone(LISTA=Index.daily.zone.rel[[i]],LISTA2=NULL,ALL.yrs=FINYEAR.daily)    
  mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("CV",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  
  dev.off()
}


#km.gn.day VS km.gn.hour
fn.fig("km.gn.day_VS_km.gn.hour",2000, 2400)     
par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for ( i in 1:N.species)
{
  #Monthly
  PLOT.km.gn._vs_km.gn.h(Index=Index.monthly.rel[[i]],HOURS=Index.monthly.hour.rel[[i]])  
  if(i==1)
  {
    mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
    legend("topright",c("km.gn.days","km.gn.hours"),bty='n',cex=1.5
           ,lty=1,lwd=2,col=c("black",rgb(red=0.1, green=0.5, blue=0.6)))
  }
  
  #Daily
  PLOT.km.gn._vs_km.gn.h(Index=Index.daily.rel[[i]],HOURS=Index.daily.hour.rel[[i]])  
  if(i==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  legend("topright",SPECIES.vec[i],bty='n',cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()


#     4.22.17 sensitivity tests
if(do.sensitivity=="YES")
{
  SCENARIOS=as.character(Tab.Sensi$Scenario)
  id.eff=match("Effective",SCENARIOS)
  nm.SCENARIOS=c(SCENARIOS[id.eff],substr(SCENARIOS[-id.eff],7,50))
  
    #Plot By species relative cpue
  CL=c("red","black","chartreuse4","orange","forestgreen","blue") 
  names(CL)=SCENARIOS
  LIN=c(rep(1,2),2,4,3,6)
  for ( i in 1:N.species)
  {
    fn.fig(paste("Sensitivity_relative.",SPECIES.vec[i],sep=""),2400, 2200)    
    par(mfcol=c(2,1),mar=c(1,3.75,2,.15),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
    if(i%in%c(1))WHERE="topright"
    if(i%in%c(2,3))WHERE="top"
    if(i%in%c(4))WHERE="topleft"
    if(!Q_change=="1982-83")fn.show.sens(FOLY=List.foly.nom.rel[[i]]$Foly,STAND=Sens.list.monthly.rel[[i]],
                 STAND.d=Sens.list.daily.rel[[i]],COMBO="NO")
     mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.75,outer=T)
    mtext("Relative catch rate",side=2,line=-1.5,font=1,las=0,cex=1.75,outer=T)   
    dev.off()
  }
  

  #Plot Species Combined relative cpue
  CL=c("grey80","black","grey50","grey10","grey40","grey70")
  names(CL)=SCENARIOS
  
  fn.fig("Figure 4_Sensitivity",2400, 2400)    
  par(mfrow=c(4,2),mar=c(1,3.75,1.5,.5),oma=c(2.5,.1,.1,.5),las=1,mgp=c(1.9,.6,0))
  for ( i in 1:N.species)
  {
    WHERE=NA
    if(i%in%c(1))WHERE="topright"
    if(!Q_change=="1982-83")fn.show.sens(FOLY=List.foly.nom.rel[[i]]$Foly,STAND=Sens.list.monthly.rel[[i]],
                                         STAND.d=Sens.list.daily.rel[[i]],COMBO="YES")
     mtext("Financial year",side=1,line=1.20,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative catch rate",side=2,line=-1.6,font=1,las=0,cex=1.5,outer=T)  
  }
  dev.off()
}


#     4.22.18 Influence Plots     
#note: Explain the influence of the different model terms(Bentley et al 2012)
if(do.influence=="YES")
{
  #Postive catch only
  
  #CDI plots
      #MONTHLY  
  Terms=vector('list',N.species)
  names(Terms)=names(DATA.list.LIVEWT.c)
  Term.Type=Names.vec=Terms
  for(i in 1:N.species)
  {
    a=gsub(" ", "", unlist(strsplit(as.character(Best.Model[[i]]$Log)[3],"[+]")), fixed = TRUE)
    IDD=match("offset(log.Effort)",a)
    if(length(IDD)>0) a=a[-IDD]
    a=a[-match("FINYEAR",a)]
    Terms[[i]]=a  
    
    b=a
    b=ifelse(b%in%c("BLOCKX","VESSEL","MONTH"),"CAT","Cont")
    Term.Type[[i]]=b
    
    nm=a
    nm=ifelse(nm=="BLOCKX","Block",ifelse(nm=="MONTH","Mn",ifelse(nm=="VESSEL","Vessel",
              ifelse(nm=="log.Catch.Dusky","Dusky catch",
              ifelse(nm=="log.Catch.Scalefish","Scalefish catch",
              ifelse(nm=="log.Catch.Whiskery","Whiskery catch",NA))))))
    Names.vec[[i]]=nm
  }
  YLIM.vec=list(c(.7,1.3),c(.6,1.25),c(.5,1.65),c(.7,2.55))
  SCALE.vec=list(c(25,30,20),c(25,30,20),c(25,30,20),c(20,30,20))
  Store.Influence=vector('list',N.species)
  names(Store.Influence)=SPECIES.vec
  HnDll="C:/Matias/Analyses/Catch and effort/Outputs/Influence.plot/Monthly."
  for ( i in 1:N.species)
  {
    DD=Stand.out[[i]]$"Stand_Base case" 
    Store.Influence[[i]]=Influence.fn(DD$GLMlog,DD$PosData,Term.Type[[i]],
        Terms[[i]],Names.vec[[i]],YLIM.vec[[i]],SCALE.vec[[i]],add.Influence="NO")
  }
  
      #DAILY
  Terms=vector('list',N.species)
  names(Terms)=names(DATA.list.LIVEWT.c.daily)
  Term.Type=Names.vec=Terms
  for(i in 1:N.species)
  {
    a=gsub(" ", "", unlist(strsplit(as.character(Best.Model.daily[[i]]$Log)[3],"[+]")), fixed = TRUE)
    IDD=match("offset(log.Effort)",a)
    if(length(IDD)>0) a=a[-IDD]
    a=a[-match("FINYEAR",a)]
    Terms[[i]]=a  
    
    b=a
    b=ifelse(b%in%c("block10","BLOCKX","VESSEL","MONTH"),"CAT","Cont")
    Term.Type[[i]]=b
    
    nm=a
    nm=ifelse(nm=="block10","Block",ifelse(nm=="MONTH","Mn",ifelse(nm=="VESSEL","Vessel",
                                                                   ifelse(nm=="log.Catch.Dusky","Dusky catch",
                                                                          ifelse(nm=="log.Catch.Scalefish","Scalefish catch",
                                                                                 ifelse(nm=="log.Catch.Whiskery","Whiskery catch",NA))))))
    Names.vec[[i]]=nm
  }
  YLIM.vec.d=list(c(.9,1.06),c(.9,1.08),c(.85,1.15),c(0.6,1.6))
  SCALE.vec.d=list(c(50,45,50),c(45,40,60),c(35,40,30),c(30,40,20))
  Store.Influence.daily=vector('list',N.species)
  names(Store.Influence.daily)=SPECIES.vec
  HnDll="C:/Matias/Analyses/Catch and effort/Outputs/Influence.plot/Daily."
  for ( i in 1:N.species)
  {
    DD=Stand.out.daily[[i]]$"Stand_Base case"
    Store.Influence.daily[[i]]=Influence.fn(DD$GLMlog,DD$PosData,Term.Type[[i]],
              Terms[[i]],Names.vec[[i]],YLIM.vec.d[[i]],SCALE.vec.d[[i]],add.Influence="NO")
  }
  

  #Figure 5. Compare influence of all terms  
  # LWD=3
  # LTY.col=c("black","grey60","black","grey45","grey55")
  # Whr=c("top","bottomright","topright","topright")
  # Whr.d=c("top","bottomright","bottomright","topright")
  # fn.fig("Figure 5.All.terms.Influence",2400, 2400)
  # par(mfrow=c(4,2),mar=c(1,3,1.5,.9),oma=c(2.5,1,.1,.2),las=1,mgp=c(1.9,.7,0))    
  # for ( i in 1:N.species)
  # {
  #   #Monthly
  #   Fig3.plot.fun(A=Store.Influence[[i]],YLIM=YLIM.vec[[i]],Whr[i])
  #   if(i==1)mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  #   
  #   #Daily
  #   Fig3.plot.fun(A=Store.Influence.daily[[i]],YLIM=YLIM.vec.d[[i]],Whr.d[i])
  #   if(i==1)mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  #   
  #   legend("bottomleft",SPECIES.vec[i],bty='n',cex=1.5)
  # }
  # mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  # mtext("Influence",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  # dev.off()
  
}


#     4.22.19 Compare creep and non-creep       
for ( i in 1:N.species)
{
  fn.fig(paste("Effort_creep/",SPECIES.vec[i],sep=""),2400, 2200)    
  par(mfcol=c(2,1),mar=c(1,3.75,2,.15),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
  #Monthly
  if(!Q_change=="1982-83")Compare.Eff.creep.No.creep(Index=Index.monthly.rel[[i]],
                                                     Index_no.creep=Index.monthly.rel_no.creep[[i]],
                                                     Index_q.change=NULL,
                                                     Index_q.change_no.creep=NULL,
                                                     ADD.CI="NO") 
  
  if(Q_change=="1982-83")Compare.Eff.creep.No.creep(Index=Index.monthly.rel[[i]],
                                                    Index_no.creep=Index.monthly.rel_no.creep[[i]],
                                                    Index_q.change=Index.monthly_q.change.rel[[i]],
                                                    Index_q.change_no.creep=Index.monthly_q.change.rel_no.creep[[i]],
                                                    ADD.CI="NO") 
  #Daily
  Compare.Eff.creep.No.creep(Index=Index.daily.rel[[i]],
                             Index_no.creep=Index.daily.rel_no.creep[[i]],
                             Index_q.change=NULL,
                             Index_q.change_no.creep=NULL,
                             ADD.CI="NO")  
  mtext("Relative catch rate",side=2,line=-1.5,font=1,las=0,cex=1.75,outer=T)   
  mtext("Financial year",side=1,line=1,cex=1.75,outer=T)
  legend("bottomleft",c("Without eff.creep","With eff. creep"),lwd=2,lty=c(2,1),bty='n',cex=1.25)
  dev.off()
}


#     4.22.20 Export relative index 
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Index")

#Relative unstandardised      
for ( i in 1:N.species)
{
  nm=paste(SPECIES.vec[i],".annual.folly.csv",sep="")
  write.csv(List.foly.nom.rel[[i]]$Foly,nm,row.names=F)  
  
  nm=paste(SPECIES.vec[i],".annual.nominal.csv",sep="")
  write.csv(List.foly.nom.rel[[i]]$Nom,nm,row.names=F)    
  
}

#Relative Standardised base case (zones combined)
  # with effort creep
for ( i in 1:N.species)
{  
  #km gn days
    #Monthly
  dat=Index.monthly.rel[[i]] 
  if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") dat=rbind(dat,Index.monthly_q.change.rel[[i]])
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
    #Daily
  dat=Index.daily.rel[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
  #km gn hours
    #Monthly
  dat=Index.monthly.hour.rel[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly_hours.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
    #Daily
  dat=Index.daily.hour.rel[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily_hours.csv",sep="")
  write.csv(dat,nm,row.names=F)
}

  # with no effort creep
for ( i in 1:N.species)
{  
  #km gn days
    #Monthly
  dat=Index.monthly.rel_no.creep[[i]] 
  if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") dat=rbind(dat,Index.monthly_q.change.rel_no.creep[[i]])
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly_no.creep.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
    #Daily
  dat=Index.daily.rel_no.creep[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily_no.creep.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
  #km gn hours
    #Monthly
  dat=Index.monthly.hour.rel_no.creep[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly_no.creep_hours.csv",sep="")
  write.csv(dat,nm,row.names=F)
  
    #Daily
  dat=Index.daily.hour.rel_no.creep[[i]] 
  names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
  nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily_no.creep_hours.csv",sep="")
  write.csv(dat,nm,row.names=F)
}


#Relative Standardised base case by zone
  # with effort creep
for ( i in 1:N.species)
{    
  #Monthly
  D=Index.monthly.zone.rel[[i]]
  if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") D_q.change=Index.monthly.zone_q.change.rel[[i]]
  zn=names(D)
  for(z in 1:length(zn))
  {
    dat=D[[z]] 
    if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") dat=rbind(dat,D_q.change[[z]])
    
    if(!is.null(dat))
    {        
      names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
      nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly.",Zn[z],".csv",sep="")
      write.csv(dat,nm,row.names=F)        
    }
  }
  
  
  #Daily
  D=Index.daily.zone.rel[[i]]     
  zn=names(D)
  for(z in 1:length(zn))
  {
    dat=D[[z]] 
    if(!is.null(dat))
    {        
      names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
      nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily.",Zn[z],".csv",sep="")
      write.csv(dat,nm,row.names=F)        
    }
  }
}

  # with no effort creep
for ( i in 1:N.species)
{    
  #Monthly
  D=Index.monthly.zone.rel_no.creep[[i]]
  if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") D_q.change=Index.monthly.zone_q.change.rel_no.creep[[i]]
  zn=names(D)
  for(z in 1:length(zn))
  {
    dat=D[[z]]
    if(Q_change=="1982-83")if(SPECIES.vec[i]=="Whiskery shark") dat=rbind(dat,D_q.change[[z]])
    if(!is.null(dat))
    {        
      names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
      nm=paste(SPECIES.vec[i],".annual.abundance.basecase.monthly_no.creep.",Zn[z],".csv",sep="")
      write.csv(dat,nm,row.names=F)        
    }
  }
  
  
  #Daily 
  D=Index.daily.zone.rel_no.creep[[i]]     
  zn=names(D)
  for(z in 1:length(zn))
  {
    dat=D[[z]] 
    if(!is.null(dat))
    {        
      names(dat)=c("Finyear","Mean","SD","CV","LOW.CI","UP.CI")  
      nm=paste(SPECIES.vec[i],".annual.abundance.basecase.daily_no.creep.",Zn[z],".csv",sep="")
      write.csv(dat,nm,row.names=F)        
    }
  }
}




#     4.22.21  AMM actions
  #2017 AMM actions
if(do.actions.AMM.2017=="YES")
{
  Best.Model_AMM.2017.action=Best.Model
  Stand.out_AMM.2017.action=vector('list',length=N.species)
  names(Stand.out_AMM.2017.action)=SPECIES.vec
  Stand.out.daily_AMM.2017.action=Stand.out.hour_AMM.2017.action=
    Stand.out.daily.hour_AMM.2017.action=Stand.out_AMM.2017.action
  
  system.time(for ( n in 1:N.species)    #360 seconds to run
  {
    #Create storing lists
    Zens=Tab.Sensi[match("Stand_Base case",Tab.Sensi$Scenario),]
    Scen.store=vector('list',nrow(Zens))
    names(Scen.store)=Zens$Scenario
    Scen.store.daily=Scen.store.hour=Scen.store.daily.hour=Scen.store
    
    #zones combined
    for(i in 1:length(Scen.store))
    {
      ZEN=as.character(Zens$Scenario[i])
      Blocks_used=as.character(Zens[i,]$Blocks_used)
      Vessels_used=as.character(Zens[i,]$Vessels_used)
      
      #Monthly  
      DD=DATA.list.LIVEWT.c[[n]]
      
      Scen.store$Nominal=aggregate((Catch.Target/Km.Gillnet.Days.c)~FINYEAR+BLOCKX,subset(DD,Catch.Target>0),mean,na.rm=T)   #add nominal
      Scen.store.hour$Nominal=aggregate((Catch.Target/Km.Gillnet.Hours.c)~FINYEAR+BLOCKX,subset(DD,Catch.Target>0),mean,na.rm=T)
      
      if(Blocks_used=="Subset")  DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used[[n]],1,4)))      
      if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used[[n]] 
                                                & FINYEAR %in% subset(n.vesls.per.yer.mn[[n]],single>=Threshold.n.vessls.per.yr)$FINYEAR) 
      Scen.store[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model_AMM.2017.action[[n]]$Bi,
                                    Formula.cpue.pos=Best.Model_AMM.2017.action[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")
      
      Scen.store.hour[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model_AMM.2017.action[[n]]$Bi,
                                         Formula.cpue.pos=Best.Model_AMM.2017.action[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Hours")
      
      #Daily
      DD=DATA.list.LIVEWT.c.daily[[n]]
      
      Scen.store.daily$Nominal=aggregate((Catch.Target/Km.Gillnet.Days.c)~FINYEAR+BLOCKX,subset(DD,Catch.Target>0),mean,na.rm=T)   #add nominal
      Scen.store.daily.hour$Nominal=aggregate((Catch.Target/Km.Gillnet.Hours.c)~FINYEAR+BLOCKX,subset(DD,Catch.Target>0),mean,na.rm=T)
      
      
      if(Blocks_used=="Subset")  DD=subset(DD,BLOCKX%in%as.numeric(substr(BLKS.used.daily[[n]],1,4)))
      if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used.daily[[n]]) 
      if(Vessels_used=="Indicative")  DD=subset(DD,VESSEL%in%VES.used.daily[[n]] 
                                                & FINYEAR %in% subset(n.vesls.per.yer.daily[[n]],single>=Threshold.n.vessls.per.yr)$FINYEAR) 
      Scen.store.daily[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model_AMM.2017.action[[n]]$Bi,
                                          Formula.cpue.pos=Best.Model_AMM.2017.action[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Days")
      
      Scen.store.daily.hour[[i]]=fn.stand.cpue(DataFile=DD,Formula.cpue=Best.Model_AMM.2017.action[[n]]$Bi,
                                               Formula.cpue.pos=Best.Model_AMM.2017.action[[n]]$Log,MESH="NO",Effort.calc="Km.Gillnet.Hours")
      
      rm(DD)
    }
    Stand.out_AMM.2017.action[[n]]=Scen.store
    Stand.out.daily_AMM.2017.action[[n]]=Scen.store.daily
    
    Stand.out.hour_AMM.2017.action[[n]]=Scen.store.hour
    Stand.out.daily.hour_AMM.2017.action[[n]]=Scen.store.daily.hour
  })
  
  #Predict year-blk
  fn.pred.cpue.AMM.2017=function(MOD,FORMULA)
  {
    GLMbi=MOD$GLMbi
    GLMlog=MOD$GLMlog
    BiData=MOD$BiData
    PosData=MOD$PosData
    Formula.cpue=FORMULA$Bi
    Formula.cpue.pos=FORMULA$Log
    
    #Get model terms
    bi.terms=all.vars(Formula.cpue)
    pos.terms=all.vars(Formula.cpue.pos)
    
    Find.Blk.10=match("block10",pos.terms)
    
    #Predict year-block
    if(is.na(Find.Blk.10))
    {
      bi.terms=bi.terms[-match(c("Catch.Target", "FINYEAR","BLOCKX"),bi.terms)]
      pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","BLOCKX"),pos.terms)]  
      Bi.dat.pred=fn.match.terms(bi.terms,BiData,GLMbi)
      Pos.dat.pred=fn.match.terms(pos.terms,PosData,GLMlog)
    }
    if(!is.na(Find.Blk.10))
    {
      bi.terms= bi.terms[-match(c("Catch.Target", "FINYEAR","block10"),bi.terms)]
      pos.terms= pos.terms[-match(c("log.Catch", "FINYEAR","block10"),pos.terms)]   
      Bi.dat.pred=fn.match.terms.daily(bi.terms,BiData,GLMbi)
      Pos.dat.pred=fn.match.terms.daily(pos.terms,PosData,GLMlog)
    }
    Bi.dat.new=Bi.dat.pred$Value
    Pos.dat.new=Pos.dat.pred$Value
    Bi.dat.new$prob=predict(GLMbi,Bi.dat.new,type="response")
    Pos.dat.new$pos=predict(GLMlog,Pos.dat.new,type="response")  
    if(is.na(Find.Blk.10))
    {
      Index=merge(Bi.dat.new[,match(c("FINYEAR","BLOCKX","prob"),names(Bi.dat.new))],
                  Pos.dat.new[,match(c("FINYEAR","BLOCKX","pos"),names(Pos.dat.new))],by=c("FINYEAR","BLOCKX"))    
    }
    if(!is.na(Find.Blk.10))
    {
      Index=merge(Bi.dat.new[,match(c("FINYEAR","block10","prob"),names(Bi.dat.new))],
                  Pos.dat.new[,match(c("FINYEAR","block10","pos"),names(Pos.dat.new))],by=c("FINYEAR","block10"))    
    }
    Index$Pos=exp(Index$pos)
    
    #Combine prob of catch and pos catch
    Index$Index=Index$Pos*Index$prob
    
    Index$FINYEAR=as.character(Index$FINYEAR)
    Index=Index[order(Index$FINYEAR),]
    
    return(Index)  
  }
  
  AMM.2017.action.monthly=vector('list',N.species)
  names(AMM.2017.action.monthly)=SPECIES.vec
  AMM.2017.action.daily=AMM.2017.action.monthly
  for ( n in 1:N.species)
  {
    #Monthly
    AMM.2017.action.monthly[[n]]=fn.pred.cpue.AMM.2017(MOD=Stand.out_AMM.2017.action[[n]][[1]],FORMULA=Best.Model_AMM.2017.action[[n]])
    #Daily
    AMM.2017.action.daily[[n]]=fn.pred.cpue.AMM.2017(MOD=Stand.out.daily_AMM.2017.action[[n]][[1]],FORMULA=Best.Model_AMM.2017.action[[n]])
  }
  
  #normalise
  AMM.2017.action.monthly.rel=AMM.2017.action.monthly
  AMM.2017.action.daily.rel=AMM.2017.action.daily
  for ( n in 1:N.species)
  {
    #Monthly
    YrS=sort(unique(AMM.2017.action.monthly[[n]]$FINYEAR))
    Lst=vector('list',length(YrS))
    Lst.nom=Lst
    for(yy in 1:length(YrS))
    {
      dummy=subset(AMM.2017.action.monthly[[n]],FINYEAR==YrS[yy])
      dummy$Pos=dummy$Pos/mean(dummy$Pos)
      dummy$Index=dummy$Index/mean(dummy$Index)
      Lst[[yy]]=dummy
      
      dummy.nom=subset(Stand.out_AMM.2017.action[[n]]$Nominal,FINYEAR==YrS[yy])
      dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`=dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`/mean(dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`)
      Lst.nom[[yy]]=dummy.nom
    }
    AMM.2017.action.monthly.rel[[n]]=do.call(rbind,Lst)
    Stand.out_AMM.2017.action[[n]]$Nominal.rel=do.call(rbind,Lst.nom)
    
    #Daily
    YrS=sort(unique(AMM.2017.action.daily[[n]]$FINYEAR))
    Lst=vector('list',length(YrS))
    Lst.nom=Lst
    for(yy in 1:length(YrS))
    {
      dummy=subset(AMM.2017.action.daily[[n]],FINYEAR==YrS[yy])
      dummy$Pos=dummy$Pos/mean(dummy$Pos)
      dummy$Index=dummy$Index/mean(dummy$Index)
      Lst[[yy]]=dummy
      
      dummy.nom=subset(Stand.out.daily_AMM.2017.action[[n]]$Nominal,FINYEAR==YrS[yy])
      dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`=dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`/mean(dummy.nom$`(Catch.Target/Km.Gillnet.Days.c)`)
      Lst.nom[[yy]]=dummy.nom
    }
    AMM.2017.action.daily.rel[[n]]=do.call(rbind,Lst)
    Stand.out.daily_AMM.2017.action[[n]]$Nominal.rel=do.call(rbind,Lst.nom)
  }
  
  #remove efficience          #ACA: this is not required as non.creep effort is used in glm
  for ( n in 1:N.species)
  {
    AMM.2017.action.monthly.rel[[n]]=fn.remove.eff.creep(d=AMM.2017.action.monthly.rel[[n]],VRS=c("Pos","Index"))
    Stand.out_AMM.2017.action[[n]]$Nominal.rel=fn.remove.eff.creep(d=Stand.out_AMM.2017.action[[n]]$Nominal.rel,VRS=c("(Catch.Target/Km.Gillnet.Days.c)"))
    AMM.2017.action.daily.rel[[n]]=fn.remove.eff.creep(d=AMM.2017.action.daily.rel[[n]],VRS=c("Pos","Index"))
    Stand.out.daily_AMM.2017.action[[n]]$Nominal.rel=fn.remove.eff.creep(d=Stand.out.daily_AMM.2017.action[[n]]$Nominal.rel,VRS=c("(Catch.Target/Km.Gillnet.Days.c)"))
  }
  
  #Year-Block effect   
  hnDl.term="C:/Matias/Analyses/Catch and effort/Outputs/Action_AMM/"
  numInt=50
  colfunc <- colorRampPalette(c("grey95","black"))
  couleurs=colfunc(numInt)
  
  colfunc <- colorRampPalette(c("yellow","red"))
  couleurs=colfunc(numInt)  
  
  fn.show.AMM.2017=function(dat,YRR)
  {
    dat$LAT=-as.numeric(substr(dat$BLOCKX,1,2))
    dat$LONG=100+as.numeric(substr(dat$BLOCKX,3,4))
    All.blks=as.character(unique(dat$BLOCKX))
    all.lat=seq(min(dat$LAT),max(dat$LAT))
    all.lon=seq(min(dat$LONG),max(dat$LONG))
    
    MAT=subset(dat,FINYEAR==YRR)
    Misin.lat=all.lat[which(!all.lat%in%MAT$LAT)]
    Misin.lon=all.lon[which(!all.lon%in%MAT$LONG)]
    if(length(Misin.lat)>0)
    {
      Add=data.frame(FINYEAR=unique(MAT$FINYEAR),BLOCKX=paste(-Misin.lat,"dummy",sep=""),prob=NA,pos=NA,Pos=NA,Index=NA,
                     LONG=115, LAT=Misin.lat)
      MAT=rbind(MAT,Add)
    }
    
    if(length(Misin.lon)>0)
    {
      Add=data.frame(FINYEAR=unique(MAT$FINYEAR),BLOCKX=paste("dummy",Misin.lon-100,sep=""),prob=NA,pos=NA,Pos=NA,Index=NA,
                     LONG=Misin.lon, LAT=-32)
      MAT=rbind(MAT,Add)
    }
    dat=reshape(MAT[match(c("Index","LONG","LAT"),names(MAT))],
                idvar="LONG",timevar="LAT",v.names="Index", direction="wide")
    names(dat)[2:ncol(dat)]=substr(names(dat)[2:ncol(dat)],start=7,stop=15)
    dat=dat[order(dat$LONG),]
    Lon=dat$LONG
    d=dat[,2:ncol(dat)]
    Sort.Lat=as.character(sort(as.numeric(colnames(d))))
    d=d[,match(Sort.Lat,colnames(d))]
    d=as.matrix(d)
    Lat=as.numeric(colnames(d))
    
    Lat=Lat-.5
    Lon=Lon+.5
    
    image(Lon,Lat,d,col =couleurs,breaks=Breaks,xaxt='n',
          yaxt='n',ylab="",xlab="",main="",ylim=c(-36,-26),xlim=c(113,129))
    box()
    mtext(YRR,3,-1.3)
    axis(1,113:129,F,tck=-0.03,cex=.65)
    axis(2,-36:-26,F,tck=-0.03,cex=.65)
    
    axis(1,seq(113,129,2),seq(113,129,2),tck=-0.06,cex=.65)
    axis(2,seq(-36,-26,2),rev(seq(26,36,2)),tck=-0.06,cex=.65)
  }
  
  
  for(n in 1:N.species)
  {
    Breaks=sort(quantile(c(AMM.2017.action.monthly.rel[[n]]$Index,
                           AMM.2017.action.daily.rel[[n]]$Index),probs=seq(0,1,1/numInt),na.rm=T))
    YRs=sort(unique(AMM.2017.action.monthly.rel[[n]]$FINYEAR))
    YRs.daily=sort(unique(AMM.2017.action.daily.rel[[n]]$FINYEAR))
    
    fn.fig(paste(hnDl.term,SPECIES.vec[n],sep=""),2400, 2400)
    smart.par(n.plots=length(c(YRs,YRs.daily)),MAR=c(2,2,.1,.1),OMA=c(2,2,.5,.5),MGP=c(1,.5,0))
    
    #Monthly
    for(y in 1:length(YRs))fn.show.AMM.2017(dat=AMM.2017.action.monthly.rel[[n]],YRR=YRs[y])
    
    #Daily
    for(y in 1:length(YRs.daily))fn.show.AMM.2017(dat=AMM.2017.action.daily.rel[[n]],YRR=YRs.daily[y])
    
    mtext("Latitude (ºS)",2,line=0,outer=T,cex=1.5,las=3)
    mtext("Longitude (ºE)",1,line=0.5,outer=T,cex=1.5,las=1)
    plot(1:10,1:10,ann=F,axes=F,col="transparent")
    color.legend(3,3,9,9,legend=c("Low",rep("",numInt-1),"High"),
                 rect.col=couleurs,gradient="y",col=1,cex=1.1)
    dev.off()
  }
}