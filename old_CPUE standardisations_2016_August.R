#--------- CPUE STANDARDISATIONS OF DIFFERENT DATASETS ---------#

#NOTE:  This script standardises the catch and effort data for the 4 commercial shark species,
#         and compares nominal (folly) and improved cpue (industry involvemnt and scientific rigour)
#       Dusky includes C. obscurus and C. brachyurus

# Error structures: Delta-lognormal and Delta-gamma as response variable is catch (in kg).
#                   Hence, Zero-inflated Poisson or Neg Bin not used


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
#           4.1. Add SOI and Freo sea level data
#           4.2 Deal with zone1-zone2 Boundary blocks to a zone
#           4.3 Extract number of vessels per species range
#           4.4 Data fixes
#           4.5 Remove blocks with few sharks
#           4.6 Remove daily data from monthly records
#           4.7 Remove spatial closures 
#           4.8 Data explorations
#           4.9 Remove NA effort
#           4.10 Put data into a list 
#           4.11 Illustrate folly effect: mean Vs sum
#           4.12 Extract number of records in glms
#           4.13 Cluster analysis of daily data to identify targeting behaviour
#           4.14 Construct wide database for analysis
#           4.15 Keep records within qualification levels
#           4.16 Compute foly index for exporting
#           4.17 Compute nominal for exporting
#           4.18 Plot vessels and catch
#           4.19 Input data tables
#           4.20 Plot cpue by year and month
#           4.21 Check outliers in catch and effort for removing nonsense values
#           4.22 Plot data gaps
#           4.23 Construct index of abundance
#             4.23.1 Explore data sets
#	            4.23.2 Catch VS CPUE
#	            4.23.3 Number of vessels and blocks per zone and species
#	            4.23.4 Combine vessels in categories
#	            4.23.5 Explore predictor's effect
#	            4.23.6 Define best model and error structure
#             4.23.7 Best model
#             4.23.8 Standardise catches
#             4.23.9 Predict standardised catch
#             4.23.10 Sensitivity tests
#             4.23.11 Imputation of missing Year-block
#               4.23.11.1 Construct index
#           		4.23.11.2 Compare indices in sensitivity test
#           		4.23.11.3 Construct indices for sensitivity scenarios
#             4.23.12 Plot sandbar spatio-temporal catches
#             4.23.13 Explain the influence of the different model terms
#             4.23.14  Build confidence intervals and abundance index for stock assessments



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

setwd("C:/Matias/Analyses/SOURCE_SCRIPTS")
source("Delta_lognormal.R")
source("Delta_gamma.R")
#source("Bootstrap_Delta_Methods.R")
source("Compare.error.structure.R")
source("Deviance.explained.R")
source("Sorting.objects.R")




#----1. DATA SECTION-----#

#1.1 Import data

setwd("C:/Matias/Analyses/Catch and effort")
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

#Southern Oscillation Index
SOI=read.csv("C:/Matias/Data/SOI.1975_2013.csv")

#Mean Freo sea level
Freo=read.csv("C:/Matias/Data/Freo_mean_sea_level.csv")  
names(Freo)[3]="Freo"

#note: from 10/2013 is dummy because the data were not #available, however, Freo sea level is not a significant term
Freo.missing=rbind(subset(Freo,Year==2012 & Month%in% 11:12),subset(Freo,Year==2012 & Month%in% 1:12)
                   ,subset(Freo,Year==2012 & Month%in% 1:12))
Freo.missing$Year=c(rep(2013,2),rep(2014,12),rep(2015,12))
Freo=rbind(Freo,Freo.missing)


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
  Whis.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_whiskery_30_70.csv")
  Gum.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_gummy_les_70.csv")
  Dusky.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_dusky_less_60.csv")
  Sand.fishArea=read.csv("C:/Matias/Data/Catch and Effort/FishableArea/BLOCKX_sandbar_30_120.csv")  
}


#1.2 Control what parts of script are activated

  #1.2.1 Data controls

#Control if doing separated analysis of monthly and daily records
Separate.monthly.daily="YES"


#Control how to aggregate daily records
Use.Date="YES"

#Control if removing closures' data or not
Remove.closure="NO"


#control if showing example of folly index
Show.folly.eg="NO"

#control how to remove blocks with low catch
Remove.blk.by="blk_only"   #remove only blocks
#Remove.blk.by="yr_blk"    #remove by yr_block combos

#control if removing or reallocating border blocks
BOUND.BLK="REALLOCATE"
#BOUND.BLK="REMOVE"


  #1.2.2 Procedure controls 

#Control type of model run
Model.run="Standard"  # for running standardisations in subsequent years
#Model.run="First"    # for first time of doing standardisation. This allows selection of 
#  best model and sensitivity tests

#Control if continue exploring fit
improve="NO"

#Control if checking correlation between catch covariates
Chck.Corr="NO"

#Control if doing cluster analysis of daily data
do.cluster="NO"

#Control if comparing lognormal to gamma
Compare.best.lognormal.gamma="NO"

#Control if adding interactions to model
With.interact="YES"
#control criteria for selecting indicative vessels
#Criteria.indi='all'  #criteria for selecting sensitivity of indicative vessel
Criteria.indi='subset'  #use different criteria for indicative vessels

#control criteria for indicative vessels
second.criteria="top.percent"

#Control type of business rule for imputing missing last year of data
Second.rule="combined" #use r if negative slope, use slope if positive
#Second.rule="pop.growth.rate"       #use pop growth rate to project from last year (use average last three years)
#Second.rule="average.last.three"  #use average of last 3 years (used given low number of observations in last year)
#Second.rule="last.ob.carried.forward" 
#Second.rule="linear.extrapolation"  #linear extrapolation

impute.aver=3 #number of years used to impute mean for recent years
impute.YRs.back=10  #number of years back to impute linear extrapolation


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

#Control if sorting factor levels 
#note: this specifies the glm reference level of block and vessel
#Sort.levels="NO"
#Sort.levels="Habitat.area"
#Sort.levels="Most.common"
Sort.levels="Highest.catch"



#Control if extracting SE directly from glm 
Extract.SE="NO"

#Control if extracting glm deviance table
Extract.Deviance.table="NO" 

#Control if showing binomial fit
do.binomial.fit="NO" 

#Control if binomial part is set to same new data as positive part for predictions
Bi.same.Pos="YES"  

#Control new data factors
Sel.lev='weighted.aver'
#Sel.lev='most.common'

#Control if using same new data for all vessels as for base case
All.ves.dat.same.BC.dat="NO"  
#All.ves.dat.same.BC.dat="YES"

#control if using same new data for predicting effort creep
Creep.dat.same.BC.dat="YES"
#Creep.dat.same.BC.dat="NO"

#Control if fitting a yr-block interaction for sandbars
Interaction.sandbar="NO"  

#Control if doing bootstraps
if(Model.run=="Standard") Do.boot="YES"    #calculate CI thru bootstrapping (takes ~1 day)

#Control bootstrapping method
if(Do.boot=="YES")
{
  #Boot.method="DATA"
  Boot.method="RESIDUALS"     
  N.boot=1000
}

#control if further exploring sandbar catch rates
Do.sandbar.explore="NO"


  #1.2.3 Reporting controls

#Control if doing influence plots (Bentley et al 2012)
if(do.sensitivity=="NO") do.influence="NO"
if(do.sensitivity=="YES") do.influence='YES'

#control color of plots
#do.colors="YES"
do.colors="NO"  #grey scale

if(do.colors=="YES") what.color="cols"
if(do.colors=="NO") what.color="black"

#Control if standardised cpue are expressed in relative terms or not
RELATIVE="YES"
#RELATIVE="NO"

#Control standardised cpue relative to what (note that relative to 1st year reduces uncertainty, creates a funnel)
if(RELATIVE=="YES")
{
  relative.to="Mean.1"
  #relative.to="First.yr"
}

#Control if plotting effort creep together
compare.creep.sep="NO"





#----2. PARAMETERS SECTION-----#

Red.th=1  #percentage reduction in deviance to accept term

Per.int.bc=.1  #interval for grouping vessels
Per.int.sens=.05   #sensitivity test

Prop.ktch.exp=.9   #keep vessels that explain 90% of catch

SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")

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


#Min and Max observed FL in catch
Min.FL.w=70;Min.FL.g=60; Min.FL.d=75; Min.FL.s=50

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

Min.weight=c(.2,.2,.5,.25)
names(Min.weight)=SPECIES.vec

#Indicative vessles
Threshold.n.yrs=5 #minimum number of years reporting the species
Threshold.n.yrs.sens=10  #sensitivity

Min.rec.ves=10  #minimum accepted number of records to keep a vessel 
Min.rec.ves.sens=20   #sensitivity



#----3. FUNCTIONS SECTION-----#

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
fn.cum.ktch=function(dat,dat1,SP,Threshold)
{
  Dat=subset(dat,SPECIES==SP & LIVEWT.c>0,select=c(BLOCKX,LIVEWT.c))
  Dat1=subset(dat1,SPECIES==SP & LIVEWT.c>0,select=c(BLOCKX,LIVEWT.c))
  KTCH=aggregate(LIVEWT.c~BLOCKX,Dat,sum,na.rm=T)
  KTCH=KTCH[order(-KTCH$LIVEWT.c),]
  KTCH$Cum=round(cumsum(KTCH$LIVEWT.c)/sum(KTCH$LIVEWT.c),2)
  id=which(KTCH$Cum<=Threshold)
  TopBlok=KTCH$BLOCKX[id]
  return(list(TopBloks=TopBlok,Dropped.blks=KTCH$BLOCKX[-id]))
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
  
  plot(X,Y,xaxt='n',yaxt='n',ylab='LAT',xlab="LONG",col="transparent",
       main=paste("N blocks=",length(dat1)),cex.lab=1.5)
  
  #kept
  for(e in 1:length(LAT.kept))
  {
    dd.y=c(LAT.kept[e]-1,LAT.kept[e]-1,LAT.kept[e],LAT.kept[e])
    dd.x=c(LONG.kept[e],LONG.kept[e]+1,LONG.kept[e]+1,LONG.kept[e])
    polygon(dd.x,dd.y,col=rgb(0, 0, 1,0.25), border=1)
    text(LONG.kept[e]+0.5,LAT.kept[e]-0.5,dat1[e],cex=0.75,col=1)
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
fn.compr.cpue.ktch=function(DATA,VES.No.K,Formula.cpue.pos,Vars,MIN.Wght,Int) #compare cpue and catch as res vars
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
  
  par(mfcol=c(3,1),mai=c(.2,.2,.1,.1))
  Pos.Diag.fn(MODEL,SPECIES.vec[i])
  legend("bottomleft",paste("AIC",round(fn.AICc(MODEL),2)),col=2,bty='n')
  legend("bottomright",paste("Dev. expl.",round(fn.Dev.Exp(MODEL),2)),col=2,bty='n')
  dat$Pred=predict(MODEL,type="response")
  return(dat)    
}

#functions for checking number of records in glm
fun.check.numb.rec=function(DATA)   
{
  DATA$UNIK=with(DATA,paste(MONTH,BLOCKX,FINYEAR,VESSEL))
  return(data.length=length(unique(DATA$UNIK)))
}
fun.check.numb.rec=function(DATA)
{
  DATA$UNIK=with(DATA,paste(Same.return.SNo,MONTH,BLOCKX,FINYEAR,VESSEL))
  return(data.length=length(unique(DATA$UNIK)))
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
  Effort.data=aggregate(cbind(Km.Gillnet.Days.c,SOI,Freo,Freo.Lag6,Freo.Lag12)~zone+
                          FINYEAR+Same.return+MONTH+BLOCKX,Effort.data1,max)
  Effort.data=aggregate(cbind(Km.Gillnet.Days.c,SOI,Freo,Freo.Lag6,Freo.Lag12)~zone+
                          FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,sum)
  
  Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data1,max)
  Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return,Effort.data.inv,sum)
  
  Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return,Effort.data1,max)
  Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return,Effort.data.no.creep,sum)
  
  Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return",all.x=T)
  Effort.data=merge(Effort.data,Effort.data.no.creep,by="Same.return",all.x=T)
  
  #target species catch 
  ID=match(c(ktch),colnames(DATA))
  DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
  
  #catch targeted at other species
  DATA$Catch.Gummy=with(DATA,ifelse(SPECIES==17001,DATA[,ID],0))
  DATA$Catch.Whiskery=with(DATA,ifelse(SPECIES==17003,DATA[,ID],0))
  DATA$Catch.Dusky=with(DATA,ifelse(SPECIES%in%c(18003,18001),DATA[,ID],0))
  DATA$Catch.Sandbar=with(DATA,ifelse(SPECIES==18007,DATA[,ID],0))
  DATA$Catch.Scalefish=with(DATA,ifelse(SPECIES%in%188000:599001,DATA[,ID],0))
  DATA$Catch.Total=with(DATA,ifelse(SPECIES<1000000,DATA[,ID],0))
  
  #reshape catch data
  TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Scalefish,Catch.Dusky,
                        Catch.Sandbar,Catch.Total)~MONTH+FINYEAR+BLOCKX+VESSEL+Same.return+LAT+LONG+YEAR.c,data=DATA,sum,na.rm=T)
  TABLE=TABLE[order(TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  
  #proportion of records with target catch
  prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
  
  #merge catch and effort
  dat=merge(TABLE,Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"),all.x=T)
  
  #create "other shark catch" variable
  #   dat$Catch.other.shk=NA
  #   if(target[1]==17003)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Sandbar
  #   if(target[1]==17001)dat$Catch.other.shk=dat$Catch.Whiskery
  #   if(target[1]%in%c(18003,18001))dat$Catch.other.shk=dat$Catch.Whiskery+dat$Catch.Sandbar
  #   if(target[1]==18007)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Whiskery
  
  
  #recalculate 60 by 60 blocks
  dat$BLOCKX.orignl=dat$BLOCKX
  dat$BLOCKX=as.numeric(substr(dat$BLOCKX,1,4))  
  
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}

  #daily data
Effort.data.fun.daily=function(DATA,target,ktch)
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
  if(Use.Date=="NO")
  {
    Effort.data=aggregate(cbind(Km.Gillnet.Days.c,SOI,Freo,Freo.Lag6,Freo.Lag12)~zone+
                            FINYEAR+Same.return.SNo+MONTH+BLOCKX+block10,Effort.data1,max)
    
    Effort.data.inv=aggregate(Km.Gillnet.Days.inv~Same.return.SNo,Effort.data1,max)
    Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~Same.return.SNo,Effort.data1,max)
    
    Effort.data=merge(Effort.data,Effort.data.inv,by="Same.return.SNo",all.x=T)
    Effort.data=merge(Effort.data,Effort.data.no.creep,by="Same.return.SNo",all.x=T)
  }
  
  if(Use.Date=="YES")
  {
    Effort.data=aggregate(cbind(Km.Gillnet.Days.c,SOI,Freo,Freo.Lag6,Freo.Lag12)~zone+
                            FINYEAR+date+MONTH+BLOCKX+block10+VESSEL,Effort.data1,max)
    Effort.data.inv=aggregate(Km.Gillnet.Days.inv~FINYEAR+date+MONTH+BLOCKX+block10+VESSEL,Effort.data1,max)    
    Effort.data.no.creep=aggregate(Km.Gillnet.Days.c.no.creep~FINYEAR+date+MONTH+BLOCKX+block10+VESSEL,Effort.data1,max)    
    
    Effort.data=merge(Effort.data,Effort.data.inv,by=c("FINYEAR","date","MONTH","BLOCKX","block10","VESSEL"),all.x=T)
    Effort.data=merge(Effort.data,Effort.data.no.creep,by=c("FINYEAR","date","MONTH","BLOCKX","block10","VESSEL"),all.x=T)
  }
  
  
  #target species catch 
  ID=match(c(ktch),colnames(DATA))
  DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
  
  #catch targeted at other species
  DATA$Catch.Gummy=with(DATA,ifelse(SPECIES==17001,DATA[,ID],0))
  DATA$Catch.Whiskery=with(DATA,ifelse(SPECIES==17003,DATA[,ID],0))
  DATA$Catch.Dusky=with(DATA,ifelse(SPECIES%in%c(18003,18001),DATA[,ID],0))
  DATA$Catch.Sandbar=with(DATA,ifelse(SPECIES==18007,DATA[,ID],0))
  DATA$Catch.Scalefish=with(DATA,ifelse(SPECIES%in%188000:599001,DATA[,ID],0))
  DATA$Catch.Total=with(DATA,ifelse(SPECIES<1000000,DATA[,ID],0))
  
  #reshape catch data
  if(Use.Date=="NO")
  {
    TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Scalefish,Catch.Dusky,
                          Catch.Sandbar,Catch.Total)~MONTH+FINYEAR+BLOCKX+block10+VESSEL+Same.return.SNo+LAT+
                      LONG+YEAR.c,data=DATA,sum,na.rm=T)    
  }
  if(Use.Date=="YES")
  {
    TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Scalefish,Catch.Dusky,
                          Catch.Sandbar,Catch.Total)~MONTH+FINYEAR+BLOCKX+block10+VESSEL+date+LAT+LONG+
                      YEAR.c,data=DATA,sum,na.rm=T)    
  }
  
  TABLE=TABLE[order(TABLE$date,TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  
  #proportion of records with target catch
  prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
  
  #merge catch and effort
  if(Use.Date=="NO")dat=merge(TABLE,Effort.data,by=c("Same.return.SNo","FINYEAR","MONTH","BLOCKX","block10"),all.x=T)
  if(Use.Date=="YES")dat=merge(TABLE,Effort.data,by=c("date","FINYEAR","MONTH","BLOCKX","block10","VESSEL"),all.x=T)
  
  #recalculate 60 by 60 blocks
  dat$BLOCKX.orignl=dat$BLOCKX
  dat$BLOCKX=as.numeric(substr(dat$BLOCKX,1,4))
  
  
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
  TABLE6=merge(TABLE6,data.frame(VESSEL=names(TABLE1),Count=TABLE1),by="VESSEL")
  TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
  TABLE6=TABLE6[order(TABLE6$Count),]
  
  TABLE7=aggregate(Catch.Target~BLOCKX,data=DATA,mean)
  names(TABLE7)[2]="Mean Catch"
  TABLE7.1=aggregate(Catch.Target~BLOCKX,data=DATA,sd)
  names(TABLE7.1)[2]="SD Catch"
  TABLE7=merge(TABLE7,data.frame(BLOCKX=names(TABLE2),Count=TABLE2),by="BLOCKX")
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
  TABLE16=data.frame(Records=TABLE16,VESSEL=names(TABLE16))
  TABLE16$CumRecords=cumsum(TABLE16$Records)
  TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
  TopVessels.records=subset(TABLE16,PerCumRecords<=RecordThreshold, select=c(VESSEL))
  TABLE16.b=subset(TABLE16,VESSEL%in%TopVessels.records$VESSEL,select=c("Records","VESSEL"))
  agg=subset(TABLE16,!(VESSEL%in%TopVessels.records$VESSEL))
  Per.grouped.rec.vesl=round(100*(sum(agg$Records)/(sum(TABLE16.b$Records)+sum(agg$Records))),1)
  TABLE16.b=rbind(TABLE16.b,data.frame(Records=sum(agg$Records),VESSEL=as.factor("PlusGroup")))
  
  Blocks.zones=unique(data.frame(BLOCKX=DATA$BLOCKX,zone=DATA$zone))
  TABLE17=rev(sort(table(DATA$BLOCKX)))
  TABLE17=data.frame(Records=TABLE17,BLOCKX=names(TABLE17))
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

fn.define.models=function(res.var)  #function for defining model structures
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
    pos.log.form[[1]]=paste("FINYEAR","*","BLOCKX","+","VESSEL",sep="")
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

#Comparing cpue, catch, random effects, etc functions
fn2=function(NOMINAL,MOD.logcpue,MOD.offset,MOD.main,MOD.logcpue_vessel,MOD.offset_vessel,MOD.main_vessel,
             logcpue,offset,main,logcpue_vessel,offset_vessel,main_vessel)
{
  dev.exp=paste(" (",round(100*c(Dsquared(MOD.logcpue),Dsquared(MOD.offset),Dsquared(MOD.main),
                                 Dsquared(MOD.logcpue_vessel),Dsquared(MOD.offset_vessel),
                                 Dsquared(MOD.main_vessel)),0),"%)",sep="")  
  NOMINAL[,2]=exp(NOMINAL[,2])
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
  
  points(NOMINAL[,2]/mean(NOMINAL[,2]),pch=19,cex=1.5)
  
  legend("bottomleft",
         paste(c("cpue","offset","main term","cpue_vess","offset_vess","main term_vess"),dev.exp),
         lty=c(rep(1,3),rep(3,3)),col=rep(1:3,2),bty='n',cex=1.5)
  
}

Compare.all=function(RAW,RESHAPED,RESHAPED.ALL.VES)   #function for testing models
{
  RAW$log.Effort=log(RAW$Km.Gillnet.Days.c)
  RAW$log.Catch=log(RAW$LIVEWT.c)
  RAW$cpue=RAW$LIVEWT.c/RAW$Km.Gillnet.Days.c
  RAW$log.cpue=log(RAW$cpue)
  RAW$FINYEAR=as.factor(RAW$FINYEAR)
  RAW$BLOCKX=as.factor(RAW$BLOCKX)
  RAW$MONTH=as.factor(RAW$MONTH)
  RAW$VESSEL=as.factor(RAW$VESSEL)
  
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
  Nominal_RAW=aggregate(log.cpue~FINYEAR,RAW,mean)
  Nominal_RESHAPED=aggregate(log.cpue~FINYEAR,RESHAPED,mean)
  Nominal_RESHAPED.ALL.VES=aggregate(log.cpue~FINYEAR,RESHAPED.ALL.VES,mean)
  
  #Lognormal GLM
  #RAW
  MODEL_RAW.logcpue <- glm(log.cpue ~ FINYEAR , data=RAW, family=gaussian, maxit=500)
  MODEL_RAW.offset <- glm(log.Catch ~ FINYEAR + offset(log.Effort), data=RAW, family=gaussian, maxit=500)
  MODEL_RAW.main <- glm(log.Catch ~ FINYEAR + log.Effort, data=RAW, family=gaussian, maxit=500)
  
  MODEL_RAW.logcpue_vessel <- glm(log.cpue ~ FINYEAR+VESSEL, data=RAW, family=gaussian, maxit=500)
  MODEL_RAW.offset_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + offset(log.Effort), data=RAW, family=gaussian, maxit=500)
  MODEL_RAW.main_vessel <- glm(log.Catch ~ FINYEAR+VESSEL + log.Effort, data=RAW, family=gaussian, maxit=500)
  
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
  #RAW
  Preds_RAW.logcpue=lsmeans(MODEL_RAW.logcpue,specs=c("FINYEAR"))
  Preds_RAW.offset=lsmeans(MODEL_RAW.offset,specs=c("FINYEAR"))
  Preds_RAW.main=lsmeans(MODEL_RAW.main,specs=c("FINYEAR"))
  
  Preds_RAW.logcpue_vessel=lsmeans(MODEL_RAW.logcpue_vessel,specs=c("FINYEAR"))
  Preds_RAW.offset_vessel=lsmeans(MODEL_RAW.offset_vessel,specs=c("FINYEAR"))
  Preds_RAW.main_vessel=lsmeans(MODEL_RAW.main_vessel,specs=c("FINYEAR"))
  
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
  
  
  
  #RAW
  logcpue=summary(Preds_RAW.logcpue)
  offset=summary(Preds_RAW.offset)
  main=summary(Preds_RAW.main)  
  logcpue_RAW=as.data.frame(logcpue[c("lsmean")])
  offset_RAW=as.data.frame(offset[c("lsmean")])
  main_RAW=as.data.frame(main[c("lsmean")])
  
  logcpue=summary(Preds_RAW.logcpue_vessel)
  offset=summary(Preds_RAW.offset_vessel)
  main=summary(Preds_RAW.main_vessel)  
  logcpue_RAW_vessel=as.data.frame(logcpue[c("lsmean")])
  offset_RAW_vessel=as.data.frame(offset[c("lsmean")])
  main_RAW_vessel=as.data.frame(main[c("lsmean")])
  
  
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
  
  
  #RAW
  fn2(Nominal_RAW,MODEL_RAW.logcpue,MODEL_RAW.offset,MODEL_RAW.main,
      MODEL_RAW.logcpue_vessel,MODEL_RAW.offset_vessel,MODEL_RAW.main_vessel,
      logcpue_RAW,offset_RAW,main_RAW,logcpue_RAW_vessel,offset_RAW_vessel,main_RAW_vessel)
  mtext(paste("Raw data, all blocks and vessels      ",SPECIES.vec[[i]]),side=3,line=-2)
  
  #_RESHAPED
  fn2(Nominal_RESHAPED,MODEL_RESHAPED.logcpue,MODEL_RESHAPED.offset,MODEL_RESHAPED.main,
      MODEL_RESHAPED.logcpue_vessel,MODEL_RESHAPED.offset_vessel,MODEL_RESHAPED.main_vessel,
      logcpue_RESHAPED,offset_RESHAPED,main_RESHAPED,logcpue_RESHAPED_vessel,offset_RESHAPED_vessel,main_RESHAPED_vessel)
  mtext("Reshaped data, indicative vessels",side=3,line=-2)
  
  #_RESHAPED.ALL.VES
  fn2(Nominal_RESHAPED.ALL.VES,MODEL_RESHAPED.ALL.VES.logcpue,MODEL_RESHAPED.ALL.VES.offset,MODEL_RESHAPED.ALL.VES.main,
      MODEL_RESHAPED.ALL.VES.logcpue_vessel,MODEL_RESHAPED.ALL.VES.offset_vessel,MODEL_RESHAPED.ALL.VES.main_vessel,
      logcpue_RESHAPED.ALL.VES,offset_RESHAPED.ALL.VES,main_RESHAPED.ALL.VES,
      logcpue_RESHAPED.ALL.VES_vessel,offset_RESHAPED.ALL.VES_vessel,main_RESHAPED.ALL.VES_vessel)
  mtext("Reshaped data, all vessels",side=3,line=-2)
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

fn.select.mod=function(DAT,dat.type,do.gamma)  #function for running model selection
{
  Species.model=vector('list',length=N.species)
  names(Species.model)=SPECIES.vec
  Species.model.sel=Species.model
  
  for (j in 1:N.species)
  {
    #1. drop vessels that never reported catch of target species   
    DATA=subset(DAT[[j]],!(VESSEL%in%Ves.no.ktc[[j]]))
    
    # 2. Put data in proper shape
    DataFile=DATA[,match(VARIABLES,names(DATA))]
    
    #convert to factor and sort levels
    DataFile$FINYEAR=as.factor(DataFile$FINYEAR)
    DataFile$BLOCKX=as.factor(DataFile$BLOCKX)
    DataFile$MONTH=as.factor(DataFile$MONTH)
    
    DataFile$FINYEAR=factor(DataFile$FINYEAR,levels=(names(table(DataFile$FINYEAR))))
    DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=(names(table(DataFile$BLOCKX))))
    DataFile$MONTH=factor(DataFile$MONTH,levels=(names(table(DataFile$MONTH))))
    
    #log effort
    DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c+0.00001)
    
    
    # Create binary dataset (0/1 for CPUE)
    BiData <- DataFile
    BiData[,1] <- as.numeric(DataFile[,1]>0)
    
    
    # Create positive (non-zero) dataset
    PosData <- DataFile[DataFile[,1]>0,]  
    PosData$log.Catch=log(PosData$Catch.Target)
    PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c
    PosData$log.cpue=log(PosData$cpue)
    
    #loop over all possible model structures
    for(i in 1:length(bin.form))
    {
      
      # binomial  dat.type="monthly"
      Log.GLMbi=NULL
      if((dat.type=="daily" & !names(DAT[j])=="gummy")|(dat.type=="monthly" & names(DAT[j])=="sandbar"))
      {
        Log.GLMbi <- glm(bin.form[[i]], data=BiData, family="binomial", maxit=500)
      }
      
      
      #Lognormal
      GLMlog <- glm(pos.log.form[[i]], data=PosData, family=gaussian, maxit=500)
      
      #Gamma
      GLMGamma=NULL
      if(do.gamma=="YES")GLMGamma <- glm(pos.form[[i]], data=PosData, family=Gamma(link=log), maxit=500)
      
      
      #Store full model
      form.bin[[i]]=list(Log.GLMbi=Log.GLMbi)
      form.pos.log[[i]]=list(GLMlog=GLMlog)
      form.pos.gam[[i]]=list(GLMGamma=GLMGamma)
      
      
      #free up memory
      rm(GLMlog,GLMGamma,Log.GLMbi)  
      
    }
    
    Species.model[[j]]=list(bin=form.bin,pos.log=form.pos.log,pos.gam=form.pos.gam)
  }
  
  #Select best model through BICs and AICs
  AIC.bin.mod=vector('list',length=N.species)
  names(AIC.bin.mod)=SPECIES.vec
  AIC.pos.log.mod=AIC.pos.gam.mod=BIC.bin.mod=BIC.pos.log.mod=BIC.pos.gam.mod=
    DEV.bin.mod=DEV.pos.log.mod=DEV.pos.gam.mod=DEV2.bin.mod=DEV2.pos.log.mod=DEV2.pos.gam.mod=AIC.bin.mod
  
  
  for (j in 1:N.species)
  {
    Bi=Species.model[[j]]$bin
    LogN=Species.model[[j]]$pos.log
    Gam=Species.model[[j]]$pos.gam
    
    for(i in 1:length(bin.form))
    {
      Log.GLMbi=Bi[[i]]$Log.GLMbi
      GLMlog=LogN[[i]]$GLMlog
      GLMGamma=Gam[[i]]$GLMGamma
      
      #extract BIC
      if(is.null(Log.GLMbi))BIC.form.bin[[i]]=list(BIC=NULL)
      if(!is.null(Log.GLMbi))BIC.form.bin[[i]]=list(BIC=BIC(Log.GLMbi))
      BIC.form.pos.log[[i]]=list(BIC=BIC(GLMlog))
      if(is.null(GLMGamma))BIC.form.pos.gam[[i]]=list(BIC=NULL)
      if(!is.null(GLMGamma))BIC.form.pos.gam[[i]]=list(BIC=BIC(GLMGamma))
      
      #extract AIC.c
      if(is.null(Log.GLMbi))AIC.form.bin[[i]]=list(AIC=NULL)
      if(!is.null(Log.GLMbi))AIC.form.bin[[i]]=list(AIC=fn.AICc(Log.GLMbi))
      AIC.form.pos.log[[i]]=list(AIC=fn.AICc(GLMlog))      
      if(is.null(GLMGamma))AIC.form.pos.gam[[i]]=list(AIC=NULL)
      if(!is.null(GLMGamma))AIC.form.pos.gam[[i]]=list(AIC=fn.AICc(GLMGamma))
      
      #Deviance explained by model
      if(is.null(Log.GLMbi))Dev.bin[[i]]=NULL
      if(!is.null(Log.GLMbi))Dev.bin[[i]]=Log.GLMbi$deviance
      Dev.pos.log[[i]]=GLMlog$deviance
      if(is.null(GLMGamma))Dev.pos.gam[[i]]=NULL
      if(!is.null(GLMGamma))Dev.pos.gam[[i]]=GLMGamma$deviance
      
      
      #D-squared
      if(is.null(Log.GLMbi))Dev2.bin[[i]]=NULL
      if(!is.null(Log.GLMbi))Dev2.bin[[i]]=Dsquared(Log.GLMbi, adjust = T)   
      Dev2.pos.log[[i]]=Dsquared(GLMlog, adjust = T)   
      if(is.null(GLMGamma))Dev2.pos.gam[[i]]=NULL
      if(!is.null(GLMGamma))Dev2.pos.gam[[i]]=Dsquared(GLMGamma, adjust = T)
      
      
    }
    
    #store in data frame
    BIC.bin.mod[[j]]=unlist(as.data.frame(do.call(rbind,BIC.form.bin)))
    BIC.pos.log.mod[[j]]=unlist(as.data.frame(do.call(rbind,BIC.form.pos.log)))
    BIC.pos.gam.mod[[j]]=unlist(as.data.frame(do.call(rbind,BIC.form.pos.gam)))
    
    AIC.bin.mod[[j]]=unlist(as.data.frame(do.call(rbind,AIC.form.bin)))
    AIC.pos.log.mod[[j]]=unlist(as.data.frame(do.call(rbind,AIC.form.pos.log)))
    AIC.pos.gam.mod[[j]]=unlist(as.data.frame(do.call(rbind,AIC.form.pos.gam)))
    
    DEV.bin.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev.bin)))
    DEV.pos.log.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev.pos.log)))
    DEV.pos.gam.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev.pos.gam)))
    
    DEV2.bin.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev2.bin)))
    DEV2.pos.log.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev2.pos.log)))
    DEV2.pos.gam.mod[[j]]=unlist(as.data.frame(do.call(rbind,Dev2.pos.gam)))
  }
  
  
  for (j in 1:N.species)
  {
    BIC.bin=BIC.bin.mod[[j]]
    BIC.pos.log=BIC.pos.log.mod[[j]]
    BIC.pos.gam=BIC.pos.gam.mod[[j]]
    
    AIC.bin=AIC.bin.mod[[j]]
    AIC.pos.log=AIC.pos.log.mod[[j]]
    AIC.pos.gam=AIC.pos.gam.mod[[j]]
    
    Best.BIC.bin=Best.AIC.bin=Best.dev.bin=NULL
    if(!is.null(BIC.bin))Best.BIC.bin=fn.AIC.BIC.ratio(BIC.bin)
    Best.BIC.log=fn.AIC.BIC.ratio(BIC.pos.log)
    if(!is.null(BIC.pos.gam))Best.BIC.gam=fn.AIC.BIC.ratio(BIC.pos.gam)
    
    if(!is.null(AIC.bin))Best.AIC.bin=fn.AIC.BIC.ratio(AIC.bin)
    Best.AIC.log=fn.AIC.BIC.ratio(AIC.pos.log)
    if(!is.null(AIC.pos.gam))Best.AIC.gam=fn.AIC.BIC.ratio(AIC.pos.gam)
    
    if(!is.null(BIC.bin))Best.dev.bin=fn.improved.dev(DEV.bin.mod[[j]],Red.th,names(BIC.bin))
    Best.dev.log=fn.improved.dev(DEV.pos.log.mod[[j]],Red.th,names(BIC.pos.log))
    if(!is.null(BIC.pos.gam))Best.dev.gam=fn.improved.dev(DEV.pos.gam.mod[[j]],Red.th,names(BIC.pos.gam))
    
    Species.model.sel[[j]]=list(
      Best.BIC.bin=Best.BIC.bin,Best.BIC.log=Best.BIC.log,Best.BIC.gam=Best.BIC.gam,
      Best.AIC.bin=Best.AIC.bin,Best.AIC.log=Best.AIC.log,Best.AIC.gam=Best.AIC.gam,
      Best.dev.bin=Best.dev.bin,Best.dev.log=Best.dev.log,Best.dev.gam=Best.dev.gam)
  }
  
  #3.2.4. Export best model structure and define best error structure
  #note: not using Yr-block interaction due to memory issues
  Best.Str=vector('list',N.species)
  names(Best.Str)=SPECIES.vec
  
  
  for (j in 1:N.species)
  {
    BIC.bin=BIC.bin.mod[[j]]
    BIC.pos.gam=BIC.pos.gam.mod[[j]]
    l=vector('list',length=3)
    MOD=Species.model.sel[[j]]
    for(i in 1:length(l))l[[i]]=MOD[[i]]$Best.Mod
    
    if(!is.null(BIC.bin))l=data.frame(part=names(MOD)[1:3],mod.str=names(unlist(l)),list.order=unlist(l))
    if(is.null(BIC.bin))l=data.frame(part=names(MOD)[2:3],mod.str=names(unlist(l)),list.order=unlist(l))
    
    #export best model structure
    write.csv(l,file = paste(SPECIES.vec[j],"best.model.csv",sep=""),row.names=F)
    #dput(l, file = paste("Best.model.",SPECIES.vec[j],".csv",sep="")) 
    
    
    #save best model structure
    ID=l$list.order
    if(!is.null(BIC.pos.gam))
    {
      if(!is.null(BIC.bin))Best.Str[[j]]=list(Bi=bin.form[[ID[1]]], Log=pos.log.form[[ID[2]]], Gam=pos.form[[ID[3]]])
      if(is.null(BIC.bin))Best.Str[[j]]=list(Bi=NULL, Log=pos.log.form[[ID[1]]], Gam=pos.form[[ID[2]]])
      #Define best error structure
      Pos.Log=Species.model[[j]]$pos.log[[ID[(length(ID)-1)]]][[1]]
      Pos.Gam=Species.model[[j]]$pos.gam[[ID[length(ID)]]][[1]]
      
      #plot model fit
      tiff(file=paste(SPECIES.vec[j],"Log.N.fit.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      fn.plot.diag(Pos.Log,"LogN",SPECIES.vec[j])
      dev.off()
      tiff(file=paste(SPECIES.vec[j],"Gamma.fit.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      fn.plot.diag(Pos.Gam,"Gamma",SPECIES.vec[j])
      dev.off()
      
      
      #tabulate model comparison
      AIC.pos.log=AIC.pos.log.mod[[j]][ID[(length(ID)-1)]]
      AIC.pos.gam=AIC.pos.gam.mod[[j]][ID[length(ID)]]
      AIC.compare=c(AIC.pos.log,AIC.pos.gam)
      Best.AIC=fn.AIC.BIC.ratio(AIC.compare)
      
      Dev.exp.log=Dsquared(Pos.Log, adjust = T)
      Dev.exp.gam=Dsquared(Pos.Gam, adjust = T)
      
      rm(MODEL.LIST)
      MODEL.LIST=list(LOGN=Pos.Log,GAMMA=Pos.Gam)
      Best.Err.Tab=Table.fit(MODEL.LIST,"NO")
      Best.Err.Tab$Like=c(Best.AIC$Like)
      Best.Err.Tab$Evidence=c(Best.AIC$Evidence.ratio_how.much.better)
      Best.Err.Tab$Deviance.explained=c(Dev.exp.log,Dev.exp.gam)
      write.table(Best.Err.Tab,paste(SPECIES.vec[j],".best.error.csv",sep=""),sep = ",",row.names=F)
      
    }
    
    if(is.null(BIC.pos.gam))
    {
      if(!is.null(BIC.bin))Best.Str[[j]]=list(Bi=bin.form[[ID[1]]], Log=pos.log.form[[ID[2]]], Gam=NULL)
      if(is.null(BIC.bin))Best.Str[[j]]=list(Bi=NULL, Log=pos.log.form[[ID[1]]], Gam=NULL)
      #Define best error structure
      Pos.Log=Species.model[[j]]$pos.log[[ID[(length(ID)-1)]]][[1]]
      Pos.Gam=NULL
      
      #plot model fit
      tiff(file=paste(SPECIES.vec[j],"Log.N.fit.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      fn.plot.diag(Pos.Log,"LogN",SPECIES.vec[j])
      dev.off()
      
      
      #tabulate model comparison
      AIC.pos.log=AIC.pos.log.mod[[j]][ID[(length(ID)-1)]]
      AIC.pos.gam=1e10
      AIC.compare=c(AIC.pos.log,AIC.pos.gam)
      Best.AIC=fn.AIC.BIC.ratio(AIC.compare)
      
      Dev.exp.log=Dsquared(Pos.Log, adjust = T)
      Dev.exp.gam=NULL
      
      rm(MODEL.LIST)
      MODEL.LIST=list(LOGN=Pos.Log)
      Best.Err.Tab=Table.fit(MODEL.LIST,"NO")
      Best.Err.Tab$Deviance.explained=c(Dev.exp.log)
      write.table(Best.Err.Tab,paste(SPECIES.vec[j],".best.error.csv",sep=""),sep = ",",row.names=F)
      
    }
    
  }
  
  
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

Test.models.sandbar=function(DATA,No.ves,SPEC,MIN.Wght)   #function for testing sandbar models
{
  #1. Create storing objects
  form.bin=vector("list",length=length(bin.form))
  form.pos.log=vector("list",length=length(pos.log.form))
  names(form.bin)=bin.form
  names(form.pos.log)=pos.log.form
  
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
  Store.pos=Store.bin=vector('list',length(bin.form))
  for(i in 1:length(bin.form))
  {
    Bin.mod=NULL
    if(!is.null(bin.form[[i]])) Bin.mod=glm(bin.form[[i]], data=BiData, family="binomial", maxit=500)  
    Store.bin[[i]]=Bin.mod
    rm(Bin.mod)
  }
  #Lognormal
  for(i in 1:length(pos.log.form))
  {
    LogN.model=glm(pos.log.form[[i]], data=PosData, family=gaussian, maxit=500)
    Store.pos[[i]]=LogN.model
    rm(LogN.model)
  }
  return(list(bi=Store.bin,pos=Store.pos))
}

fn.LogN.vs.Gama=function(DATA)  
{
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
  
  
  # create positive (non-zero) dataset
  PosData <- DataFile[DataFile[,1]>0,]  
  PosData$log.Catch=log(PosData$Catch.Target)
  PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c
  PosData$log.cpue=log(PosData$cpue)
  
  
  #4. Loop over all possible model structures    
  #Lognormal
  LogN.model=glm(LogN.form, data=PosData, family=gaussian, maxit=500)
  
  #Gamma
  Gamma.model=glm(Gama.form, data=PosData, family=Gamma(link=log), maxit=500)
  return(list(LOgN=LogN.model,Gama=Gamma.model))
}


#functions for Loading models and selecting best
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

fn.best.error.k_fold=function(DATA,FORMULA.pos,FORMULA.pos.log)  #function for k-fold crossvalidation
{
  # 1. Put data in proper shape
  DataFile=DATA[,match(VARIABLES,names(DATA))]
  
  #drop levels not occurring in data
  for(f in 1:ncol(DataFile))
  {
    if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
  }
  
  #select subsets
  id=round(nrow(DataFile)*training.data)
  
  list.ID=list()
  for(t in 1:K)list.ID[[t]]=sort(sample(1:nrow(DataFile),id))
  
  
  #2. loop for k-fold validation
  Cor.Log=Cor.Gama=Cor.Hur.Poi=Cor.Hur.NB=RMSE.Log=RMSE.Gama=RMSE.Hur.Poi=
    RMSE.Hur.NB=vector('list',length=K)
  
  for(k in 1:K)
  {
    # 2.1. select data subsets
    All.data.train=DataFile[list.ID[[k]],]
    All.data.test=DataFile[-list.ID[[k]],]
    
    # create binary dataset
    Bi.train=All.data.train
    Bi.train[,1]=as.numeric(Bi.train[,1]>0)
    Bi.test=All.data.test
    Bi.test[,1]=as.numeric(Bi.test[,1]>0)
    
    # create positive dataset
    Pos.train=All.data.train[All.data.train[,1]>0,]
    Pos.test=All.data.test[All.data.test[,1]>0,]    
    
    
    #round catch for discrete distributions
    All.data.train$Catch.Target=with(All.data.train,ifelse(Catch.Target>0 &Catch.Target<1,1,
                                                           round(Catch.Target)))
    All.data.test$Catch.Target=with(All.data.test,ifelse(Catch.Target>0 &Catch.Target<1,1,
                                                         round(Catch.Target)))
    
    
    #2.2. Fit different error structures
    
    #2.2.1   Binomial
    Binomial <- glm(FORMULA.pos, data=Bi.train, family="binomial", maxit=100)
    
    #2.2.2 Lognormal
    GLMlog <- glm(FORMULA.pos.log, data=Pos.train, family=gaussian, maxit=100)
    
    #2.2.3 Gamma
    GLMgam<-glm(FORMULA.pos, data=Pos.train, family=Gamma(link=log), maxit=100)
    
    #2.2.4 Zero-truncated Poisson
    #    GLM.hurdle.Pois <- hurdle(FORMULA.pos, data = All.data.train,dist = "poisson",zero.dist = "binomial")
    GLM.hurdle.Pois <- glmmadmb(FORMULA.pos,data=subset(All.data.train,Catch.Target>0),family="truncpoiss",
                                admb.opts=admbControl(maxfn=1000,imaxfn=1000))
    
    #2.2.5 Zero-truncated Negative Binomial
    #    GLM.hurdle.neg.bin <- hurdle(FORMULA.pos,data = All.data.train,dist = "negbin",zero.dist = "binomial")
    GLM.hurdle.neg.bin <- glmmadmb(FORMULA.pos,data=subset(All.data.train,Catch.Target>0),family="truncnbinom1")
    
    
    
    #2.2.6 Zero-inflated Negative Poisson
    #    Zero.inf.pois <- zeroinfl(FORMULA.pos,data = All.data.train, dist = "poisson")
    
    #2.2.7 Zero-inflated Negative Binomial
    #    Zero.inf.negbin <- zeroinfl(FORMULA.pos,data = All.data.train, dist = "negbin")
    
    
    
    #2.3. Predict test.data
    
    IDD=which(All.data.test$Catch.Target>0)
    Pred.this=All.data.test[IDD,]
    
    #Binomial
    Bin.pred=predict(Binomial,newdata=Pred.this[,2:ncol(Bi.test)], type='response')
    # Bin.pred=predict(Binomial,newdata=Bi.test[,2:ncol(Bi.test)], type='response')
    # Bin.pred=Bin.pred[IDD]  #extract only the comparable records
    
    
    #Lognormal
    biasCorr <- GLMlog$deviance/GLMlog$df.residual/2
    Log.pred=exp(predict(GLMlog,newdata=Pred.this[,2:ncol(Pos.test)],type='response')+biasCorr)
    #    Log.pred=exp(predict(GLMlog,newdata=Pos.test[,2:ncol(Pos.test)],type='response')+biasCorr)
    Log.pred=Log.pred*Bin.pred
    
    #Gamma
    Gama.pred=predict(GLMgam,newdata=Pred.this[,2:ncol(Pos.test)],type='response')
    #    Gama.pred=predict(GLMgam,newdata=Pos.test[,2:ncol(Pos.test)],type='response')
    Gama.pred=Gama.pred*Bin.pred
    
    #Zero-truncated Poisson
    Hurdle.Pois.pred=predict(GLM.hurdle.Pois,newdata=Pred.this[,2:ncol(All.data.test)],type='response')
    #    Hurdle.Pois.pred=predict(GLM.hurdle.Pois,newdata=All.data.test[,2:ncol(All.data.test)],type='response')
    #    Hurdle.Pois.pred=Hurdle.Pois.pred[IDD]
    
    #Zero-truncated Negative Binomial
    Hurdle.NB.pred=predict(GLM.hurdle.neg.bin,newdata=Pred.this[,2:ncol(All.data.test)],type='response')
    #    Hurdle.NB.pred=predict(GLM.hurdle.neg.bin,newdata=All.data.test[,2:ncol(All.data.test)],type='response')
    #    Hurdle.NB.pred=Hurdle.NB.pred[IDD]
    
    #Zero-inflated Negative Poisson
    #     ZIN.Pois.pred=predict(Zero.inf.pois,newdata=All.data.test[,2:ncol(All.data.test)],type='response')
    #     ZIN.Pois.pred=ZIN.Pois.pred[IDD]
    
    #Zero-inflated Negative Binomial
    #     ZIN.NB.pred=predict(Zero.inf.negbin,newdata=All.data.test[,2:ncol(All.data.test)],type='response')
    #     ZIN.NB.pred=ZIN.NB.pred[IDD]
    
    
    #2.4. Measure fit
    #Correlation
    Cor.Log[[k]]=cor(Pred.this[,1],Log.pred,method='pearson')
    Cor.Gama[[k]]=cor(Pred.this[,1],Gama.pred,method='pearson')
    Cor.Hur.Poi[[k]]=cor(Pred.this[,1],Hurdle.Pois.pred,method='pearson')
    Cor.Hur.NB[[k]]=cor(Pred.this[,1],Hurdle.NB.pred,method='pearson')
    
    #     cor(All.data.test[IDD,1],ZIN.Pois.pred,method='pearson')
    #     cor(All.data.test[IDD,1],ZIN.NB.pred,method='pearson')
    
    
    #RMSE
    RMSE.Log[[k]]=sqrt(mean((Pred.this[,1]-Log.pred)^2))
    RMSE.Gama[[k]]=sqrt(mean((Pred.this[,1]-Gama.pred)^2))
    RMSE.Hur.Poi[[k]]=sqrt(mean((Pred.this[,1]-Hurdle.Pois.pred)^2))
    RMSE.Hur.NB[[k]]=sqrt(mean((Pred.this[,1]-Hurdle.NB.pred)^2))
    
    #Compare error structures
    # MODELS.LIST=list("Lognormal" = GLMlog, "Gama" = GLMgam, "Hurdle.Pois" = GLM.hurdle.Pois,
    #                  "Hurdle.NB" = GLM.hurdle.neg.bin, "ZI.Pois" =Zero.inf.pois, "ZI.NB" = Zero.inf.negbin)
    #fn.compare.zeroinfl.errors(MODELS.LIST,length(coef(GLMlog)))
    
    
  }
  
  
  return(list(Cor.Log=unlist(Cor.Log),Cor.Gama=unlist(Cor.Gama),Cor.Hur.Poi=unlist(Cor.Hur.Poi),
              Cor.Hur.NB=unlist(Cor.Hur.NB),RMSE.Log=unlist(RMSE.Log),RMSE.Gama=unlist(RMSE.Gama),
              RMSE.Hur.Poi=unlist(RMSE.Hur.Poi),RMSE.Hur.NB=unlist(RMSE.Hur.NB)))
}

#function for applying best model to all data
fn.stand.cpue=function(DATA,VES.No.K,SPEC,Formula.cpue,Formula.cpue.pos,Vars,Eff.Creep)
{
  #1. drop vessels that never reported catch of target species   
  DATA=subset(DATA,!(VESSEL%in%VES.No.K))
  
  #2. put data in proper shape
  DataFile=DATA[,match(Vars,names(DATA))]
  
  #drop levels not occurring in data
  for(f in 1:ncol(DataFile))if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
  
  #convert to factor
  DataFile$FINYEAR=as.factor(DataFile$FINYEAR)
  DataFile$MONTH=as.factor(DataFile$MONTH)  
  DataFile$VESSEL=as.factor(DataFile$VESSEL) 
  
  if(Sort.levels=="Habitat.area")     
  {
    Ar=AREA.W[[i]]
    BL=unique(DataFile$BLOCKX)
    Ar=Ar[which(Ar$BLOCKX%in%BL),]
    Ar=Ar[order(-Ar$Fish.Area),]
    DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=Ar$BLOCKX)
  } 
  
  if(Sort.levels=="Most.common")
  {
    d=with(subset(DataFile,Catch.Target>0),table(BLOCKX))     
    DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=names(rev(sort(d))))       
  }
  
  if(Sort.levels=="Highest.catch")
  {
    agg=aggregate(Catch.Target~BLOCKX,DataFile,sum)
    agg=agg[order(-agg$Catch.Target),]
    DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=agg$BLOCKX)       
  }
  
  if(Sort.levels=="NO")  DataFile$BLOCKX=as.factor(DataFile$BLOCKX)    
  
  
  #DataFile$new.level.vessel=as.factor(DataFile$new.level.vessel)
  
  #log effort
  if(Eff.Creep=="YES") DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c)
  if(Eff.Creep=="NO") DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c.no.creep)
  
  # log catch other
  DataFile$log.Catch.Gummy=log(DataFile$Catch.Gummy+0.000000001)
  DataFile$log.Catch.Whiskery=log(DataFile$Catch.Whiskery+0.000000001)
  DataFile$log.Catch.Dusky=log(DataFile$Catch.Dusky+0.000000001)
  DataFile$log.Catch.Sandbar=log(DataFile$Catch.Sandbar+0.000000001)
  DataFile$log.Catch.Scalefish=log(DataFile$Catch.Scalefish+0.000000001)
  DataFile$log.Freo=log(DataFile$Freo)
  DataFile$log.Freo.Lag6=log(DataFile$Freo.Lag6)
  DataFile$log.Freo.Lag12=log(DataFile$Freo.Lag12)
  
  # Create binary dataset (0/1 for CPUE)
  BiData <- DataFile
  Id.Ktch=match("Catch.Target",names(DataFile))
  BiData[,Id.Ktch] <- as.numeric(DataFile[,Id.Ktch]>0)
  
  # Create positive (non-zero) dataset
  PosData <- DataFile[DataFile[,Id.Ktch]>0,]  
  PosData$log.Catch=log(PosData$Catch.Target)  
  if(Fit.to.what=="cpue")
  {
    if(Eff.Creep=="YES") PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c
    if(Eff.Creep=="NO") PosData$cpue=PosData$Catch.Target/PosData$Km.Gillnet.Days.c.no.creep  
    PosData$log.cpue=log(PosData$cpue)
  }
  
  
  # binomial GLM
  BiData=subset(BiData,VESSEL%in%unique(PosData$VESSEL))#make sure to keep same vessels in bin and pos data
  Log.GLMbi=NULL
  if(!is.null(Formula.cpue))Log.GLMbi <- glm(Formula.cpue, data=BiData, family="binomial", maxit=500)
  
  #Lognormal GLM
  GLMlog <- glm(Formula.cpue.pos, data=PosData, family=gaussian, maxit=500)
  
  return(list(Bi=Log.GLMbi,LogN=GLMlog,PosData=PosData,BiData=BiData))
}

BiasCor.fn=function(Median,SE) biasCorr <- exp(Median+(SE^2)/2) #function for bias corrected mean in normal space

#functions for checking model predictions VS observations
plot.blank=function(LAB) 
{
  plot(1:10,col="transparent",xaxt='n',yaxt='n',cex.lab=1.5,ylab=LAB,xlab="")
  text(5,5,"Not applicable",cex=2)
}
see.BC.fit.to.data=function(GLM,RespVar)
{
  BiData=GLM$BiData
  PosData=GLM$PosData 
  Log.GLMbi=GLM$Bi
  GLMlog=GLM$LogN 
  
  par(mfcol=c(3,1),mai=c(.6,.75,.3,.1),oma=c(1,1,1,1),las=1,mgp=c(2,.5,0))
  
  if(!is.null(Log.GLMbi))
  {
    BiData$Pred=round(predict(Log.GLMbi,type="response"))
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
  
  PosData$Pred=predict(GLMlog,type="response")
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

fn.pred.month=function(MOD,SP)
{
  LogN=MOD$LogN
  This=c("MONTH2","MONTH3","MONTH4","MONTH5","MONTH6","MONTH7","MONTH8","MONTH9",
         "MONTH10","MONTH11","MONTH12") 
  COFS=coef(LogN)[match(This,names(coef(LogN)))]
  plot(2:12,COFS,type='o',cex=2,pch=19,ylab="GLM Coefficient",xlab="Month") 
  legend("top",SP,bty='n',cex=2)
}

fn.get.stuff=function(MOD,SPEC)  #function for extracting SE and coefs
{
  Bi.estimates=NULL
  if(!SPEC=="Dusky shark")
  {
    Bi=MOD$Bi
    Bi.coef=coef(Bi)
    Bi.SE=summary(Bi)$coefficients[, 2]
    Bi.estimates=data.frame(Name=names(Bi.coef),mean=Bi.coef,SE=Bi.SE)    
    rownames(Bi.estimates)=NULL
  }
  
  LogN=MOD$LogN  
  LogN.coef=coef(LogN)
  LogN.coef=data.frame(Mean=LogN.coef,Name=names(LogN.coef))
  LogN.SE=summary(LogN)$coefficients[, 2]
  LogN.SE=data.frame(SE=LogN.SE,Name=names(LogN.SE))
  Log.estimates=merge(LogN.coef,LogN.SE,by="Name",all.x=T)
  
  return(list(Bi=Bi.estimates,LOGN=Log.estimates))
}

fn.get.Significance=function(MOD,SPEC)  #function for extracting term significance
{
  LogN=MOD$LogN
  Anova.bi=NULL
  if(!SPEC=="Dusky shark")
  {
    Bi=MOD$Bi
    Anova.bi=anova(Bi, test = "Chisq")
  }
  Anova.LogN=anova(LogN, test = "Chisq")
  return(list(Bi=Anova.bi,LogN=Anova.LogN))
}

Bi.Diag.fn=function(DAT,VES.No.K,MODEL,Int)  #function for Binomial model diagnostics
{
  if(!is.null(MODEL))
  {
    #1. select variables
    DAT=subset(DAT,!(VESSEL%in%VES.No.K))
    
    #2. put data in proper shape
    DataFile=DAT[,match(VARIABLES,names(DAT))]
    
    #drop levels not occurring in data
    for(f in 1:ncol(DataFile))if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
    
    #convert to factor and sort levels
    DataFile$FINYEAR=as.factor(DataFile$FINYEAR)
    DataFile$BLOCKX=as.factor(DataFile$BLOCKX)
    DataFile$MONTH=as.factor(DataFile$MONTH)
    
    DataFile$FINYEAR=factor(DataFile$FINYEAR,levels=(names(table(DataFile$FINYEAR))))
    DataFile$BLOCKX=factor(DataFile$BLOCKX,levels=(names(table(DataFile$BLOCKX))))
    DataFile$MONTH=factor(DataFile$MONTH,levels=(names(table(DataFile$MONTH))))
    
    #log effort
    DataFile$log.Effort=log(DataFile$Km.Gillnet.Days.c+0.00001)
    
    # Create binary dataset (0/1 for CPUE)
    BiData <- DataFile
    BiData[,1] <- as.numeric(DataFile[,1]>0)
    
    #Prediction
    BiData$Pred=predict(MODEL,type="response")
    
    #Aggregate by bins
    BiData$BINS=cut(BiData$Pred, seq(0,1,Int), include.lowest=TRUE)
    Agg.pred=aggregate(Pred~BINS,BiData,mean)
    Agg.Obs=aggregate(Catch.Target~BINS,BiData,mean)
    
    #plot
    plot(seq(0,1,Int),seq(0,1,Int),type='l',lty=2,lwd=2,col="grey",ylab="",xlab="",xaxt='n',cex.axis=1.5)
    points(Agg.pred$Pred,Agg.Obs$Catch.Target,pch=19,cex=1.75)
    axis(1,seq(0,1,0.1),label=F)
    
  }
  if(is.null(MODEL)) plot.blank("")
}
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

fn.get.coef=function(dat)    #function for defining vessel and month levels from weighted average
{
  #get number of records
  TAB=function(dat2,what) Tabla=table(dat2[,match(what,names(dat2))])
  N.rec.ves=TAB(dat$BiData,"VESSEL")
  N.rec.mnth=TAB(dat$BiData,"MONTH")
  
  #get binominal coefs
  fn.cof=function(mod,End,what)
  {
    a=substr(names(coef(mod)),1,End)
    a=which(a==what)
    return(coef(mod)[a])
  }
  Ves.bi.coef=fn.cof(dat$Bi,6,"VESSEL")
  names(Ves.bi.coef)=substr(names(Ves.bi.coef),7,15)
  Ves.bi.coef=c(0,Ves.bi.coef)  #add 0 as first level is set to 0
  names(Ves.bi.coef)[1]=names(N.rec.ves)[1]
  
  Mnth.bi.coef=fn.cof(dat$Bi,5,"MONTH")  
  if(length(Mnth.bi.coef)>0)
  {
    names(Mnth.bi.coef)=substr(names(Mnth.bi.coef),6,15)
    Mnth.bi.coef=c(0,Mnth.bi.coef)
    names(Mnth.bi.coef)[1]=names(N.rec.mnth)[1]
  }
  
  #get pos catch coefs
  Ves.pos.coef=c(0,fn.cof(dat$LogN,6,"VESSEL"))
  names(Ves.pos.coef)=names(Ves.bi.coef)
  
  Mnth.pos.coef=fn.cof(dat$LogN,5,"MONTH")
  names(Mnth.pos.coef)=substr(names(Mnth.pos.coef),6,15)
  Mnth.pos.coef=c(0,Mnth.pos.coef)
  names(Mnth.pos.coef)[1]=names(N.rec.mnth)[1]
  
  
  #calculate weighted average
  N.rec.ves=subset(N.rec.ves,names(N.rec.ves)%in%names(Ves.bi.coef))
  ves.weight.avg.bi=weighted.mean(Ves.bi.coef,(N.rec.ves/sum(N.rec.ves)))
  if(length(Mnth.bi.coef)>1)mnth.weight.avg.bi=weighted.mean(Mnth.bi.coef,(N.rec.mnth/sum(N.rec.mnth)))
  
  ves.weight.avg.pos=weighted.mean(Ves.pos.coef,(N.rec.ves/sum(N.rec.ves)))
  mnth.weight.avg.pos=weighted.mean(Mnth.pos.coef,(N.rec.mnth/sum(N.rec.mnth)))
  
  #select month and effort that is closed to mean coefficient values
  ID.ves.bi=names(Ves.bi.coef)[which.min(abs(Ves.bi.coef - ves.weight.avg.bi))]
  ID.mnth.bi=NULL
  if(length(Mnth.bi.coef)>1)ID.mnth.bi=names(Mnth.bi.coef)[which.min(abs(Mnth.bi.coef - mnth.weight.avg.bi))]
  
  ID.ves.pos=names(Ves.pos.coef)[which.min(abs(Ves.pos.coef - ves.weight.avg.pos))]
  ID.mnth.pos=names(Mnth.pos.coef)[which.min(abs(Mnth.pos.coef - mnth.weight.avg.pos))]
  return(list(ID.ves.bi=ID.ves.bi,ID.mnth.bi=ID.mnth.bi,ID.ves.pos=ID.ves.pos,ID.mnth.pos=ID.mnth.pos))
}

#functions for predicting standardised cpue
fn.extract=function(what) factor(what,levels=(names(table(what))))  
Most.common=function(DAT,what)
{
  ID=match(what,names(DAT))
  Tabla=sort(table(DAT[,ID]))
  return(names(Tabla[length(Tabla)]))
}
Mn.cov=function(DAT,what) mean(DAT[,match(what,names(DAT))])
Get.Coef=function(MOD)
{
  id=match(c("MONTH2","MONTH12"),names(coef(MOD)))
  MN=sort(coef(MOD)[id[1]:id[2]])
  
  id=match(c("new.level.vessel2","new.level.vessel10"),names(coef(MOD)))
  Ves=sort(coef(MOD)[id[1]:id[2]])
  return(list(Month=MN,VesCat=Ves))
}

predict.ktch=function(MOD,Grid.dat,SPEC)  #function for Predicting new data 
{
  Bi.mod=NULL
  if(!is.null(MOD$Bi)) Bi.mod=MOD$Bi
  Pos.mod=MOD$LogN
  
  Bi.dat=NULL
  if(!is.null(MOD$Bi)) Bi.dat=Grid.dat$bin 
  Pos.dat=Grid.dat$pos
  
  if(!is.null(MOD$Bi)) Bi.dat$Pred.Prob=predict(Bi.mod,Bi.dat,type="response")
  
  #Pos.dat$Pred.Catch=predict(Pos.mod,Pos.dat,type="response")
  Pos.Preds=predict(Pos.mod,Pos.dat,type="response",se.fit=T)
  Pos.dat$SE=Pos.Preds$se.fit
  Pos.dat$Pred.Catch=Pos.Preds$fit
  
  #reset Na coefficients to NA (predict() is replacing with dummy)
  NA.coefs=coef(Pos.mod)[is.na(coef(Pos.mod))]  
  NA.coefs=names(NA.coefs)
  
  Yr.NA=sapply(strsplit(NA.coefs,":"), "[", 1)
  BlK.NA=sapply(strsplit(NA.coefs,":"), "[", 2)
  
  Yr.NA=sapply(Yr.NA, function(x) substr(x, 8, 14))
  BlK.NA=sapply(BlK.NA, function(x) substr(x, 7, 10))
  Yr.NA_BlK.NA=paste(Yr.NA,BlK.NA)
  
  Pos.dat$FINYEAR.BLOCKX=with(Pos.dat,as.character(paste(FINYEAR,BLOCKX)))
  Pos.dat$Pred.Catch=with(Pos.dat,ifelse(FINYEAR.BLOCKX%in%Yr.NA_BlK.NA,NA,Pred.Catch))
  Pos.dat$SE=with(Pos.dat,ifelse(FINYEAR.BLOCKX%in%Yr.NA_BlK.NA,NA,SE))
  return(list(Bi.dat=Bi.dat,Pos.dat=Pos.dat))
}

impute.ktch=function(DAT) #function for imputing missing year-blocks
{
  #order data
  DAT=DAT[order(DAT$BLOCKX,DAT$FINYEAR),]  
  
  #create imputation grid
  DAT$FINYEAR=as.character(DAT$FINYEAR)
  DAT$BLOCKX=as.character(DAT$BLOCKX)
  yr=unique(DAT$FINYEAR)
  blk=unique(DAT$BLOCKX)
  n.yr=length(yr)
  n.blk=length(blk)
  MAT=matrix(DAT$Pred.Catch,nrow=n.yr,ncol=n.blk)
  colnames(MAT)=blk
  rownames(MAT)=yr
  
  MAT.SE=matrix(DAT$SE,nrow=n.yr,ncol=n.blk)
  colnames(MAT.SE)=colnames(MAT)
  rownames(MAT.SE)=rownames(MAT)
  
  
  #First imputation rule (fill in missing early years)
  id.first.NA=which(is.na(MAT[1,]))
  if(length(id.first.NA)>0)
  {
    a=MAT[,id.first.NA]    
    if(is.matrix(a))
    {
      na=ncol(a)
      for(e in 1:na)
      {
        id=which(!is.na(a[,e]))[1]-1
        if(!is.na(id))a[1:id,e]=na.aggregate(a[1:(id+impute.aver),e])[1:id]
      }
    }
    
    if(!is.matrix(a))
    {
      id=which(!is.na(a))[1]-1
      if(!is.na(id))a[1:id]=na.aggregate(a[1:(id+impute.aver)])[1:id]
    }
    
    MAT[,id.first.NA]=a 
  }
  
  #Second imputation rule (fill in missing recent years)  
  id.last.NA=which(is.na(MAT[n.yr,]))
  if(length(id.last.NA)>0)
  {
    a=MAT[,id.last.NA]
    a.SE=MAT.SE[,id.last.NA]    
    
    if(class(a)=="numeric")
    {
      no.Nas=which(!is.na(a))
      id=no.Nas[length(no.Nas)]+1
      if(Second.rule=="last.ob.carried.forward") if(length(id)>0)a[id:n.yr]=a[id-1]
      if(Second.rule=="average.last.three")
      {
        if(length(id)>0)
        {
          idd=no.Nas[(length(no.Nas)-impute.aver+1):length(no.Nas)]
          a[(no.Nas[length(no.Nas)]+1):n.yr]=mean(a[idd],na.rm=T)
        }
      }
      if(Second.rule=="linear.extrapolation")
      {      
        if(length(id)>0)
        {
          idd=no.Nas[(length(no.Nas)-min(length(no.Nas),impute.YRs.back)+1):length(no.Nas)]
          x=idd
          y=a[idd]
          SE=a.SE[idd]           
          a[(no.Nas[length(no.Nas)]+1):n.yr]=fn.linear.imp(x=x,y=y,SE=SE,pred=(x[length(x)]+1):nrow(a),cap=max(a,na.rm=T))
        }
      }
      
      if(Second.rule=="pop.growth.rate")
      {      
        if(length(id)>0)
        {
          
          idd=no.Nas[(length(no.Nas)-impute.aver+1):length(no.Nas)]
          x=idd
          y=a[idd]          
          a[(no.Nas[length(no.Nas)]+1):n.yr]=fn.pop.grw.rte.imp(x=x,y=y,pred=(x[length(x)]+1):length(a),cap=max(a,na.rm=T))
        }
      }  
      
      if(Second.rule=="combined")
      {      
        if(length(id)>0)
        {
          idd=no.Nas[(length(no.Nas)-min(length(no.Nas),impute.YRs.back)+1):length(no.Nas)]
          x=idd
          y=a[idd]
          SE=a.SE[idd]
          SLOPE=coef(lm(y~x))[2]
          
          #if slope is negative, set to r with constant set at last year of data             
          if(SLOPE<0) a[(no.Nas[length(no.Nas)]+1):n.yr]=fn.pop.grw.rte.imp(x=x,y=y[length(y)],pred=(x[length(x)]+1):length(a),cap=2*max(a,na.rm=T))
          
          #if slope is positive, set to linear increase based on data
          if(SLOPE>0) a[(no.Nas[length(no.Nas)]+1):n.yr]=fn.linear.imp(x=x,y=y,SE=SE,pred=(x[length(x)]+1):length(a),cap=2*max(a,na.rm=T))
          
        }
      }
    }  
    
    if(!class(a)=="numeric")
    {
      for(e in 1:ncol(a))
      {
        no.Nas=which(!is.na(a[,e]))
        id=no.Nas[length(no.Nas)]+1
        if(Second.rule=="last.ob.carried.forward")if(length(id)>0)a[id:n.yr,e]=a[id-1,e]
        
        if(Second.rule=="average.last.three")
        {
          if(length(id)>0)
          {
            idd=no.Nas[(length(no.Nas)-impute.aver+1):length(no.Nas)]
            a[(no.Nas[length(no.Nas)]+1):n.yr,e]=mean(a[idd,e],na.rm=T)
          }
        }
        
        if(Second.rule=="linear.extrapolation")
        {      
          if(length(id)>0)
          {
            idd=no.Nas[(length(no.Nas)-min(length(no.Nas),impute.YRs.back)+1):length(no.Nas)]
            x=idd
            y=a[idd,e]
            SE=a.SE[idd,e]  
            a[(no.Nas[length(no.Nas)]+1):n.yr,e]=fn.linear.imp(x=x,y=y,SE=SE,pred=(x[length(x)]+1):nrow(a),cap=max(a[,e],na.rm=T))
          }
        }
        
        if(Second.rule=="pop.growth.rate")
        {      
          if(length(id)>0)
          {
            idd=no.Nas[(length(no.Nas)-impute.aver+1):length(no.Nas)]
            x=idd
            y=a[idd,e]            
            a[(no.Nas[length(no.Nas)]+1):n.yr,e]=fn.pop.grw.rte.imp(x=x,y=y,pred=(x[length(x)]+1):nrow(a),cap=max(a[,e],na.rm=T))
          }
        }
        
        if(Second.rule=="combined")
        {      
          if(length(id)>0)
          {
            idd=no.Nas[(length(no.Nas)-min(length(no.Nas),impute.YRs.back)+1):length(no.Nas)]
            x=idd
            y=a[idd,e]
            SE=a.SE[idd,e]
            SLOPE=coef(lm(y~x))[2]
            
            #if slope is negative, set to r with constant set at last year of data              
            if(SLOPE<0) a[(no.Nas[length(no.Nas)]+1):n.yr,e]=fn.pop.grw.rte.imp(x=x,y=y[length(y)],pred=(x[length(x)]+1):nrow(a),cap=2*max(a[,e],na.rm=T))
            
            #if slope is positive, set to linear increase based on data
            if(SLOPE>0) a[(no.Nas[length(no.Nas)]+1):n.yr,e]=fn.linear.imp(x=x,y=y,SE=SE,pred=(x[length(x)]+1):nrow(a),cap=2*max(a[,e],na.rm=T))
            
          }
        }
        
      }
    }
    
    MAT[,id.last.NA]=a 
  }
  
  #Third imputation rule  (fill in missing in between years)
  id.inbetween.NA=MAT
  id.inbetween.NA[!is.na(id.inbetween.NA)] <- 0
  id.inbetween.NA[is.na(id.inbetween.NA)] <- 1
  id.inbetween.NA=colSums(id.inbetween.NA)[colSums(id.inbetween.NA)>=1]
  id.inbetween.NA=match(names(id.inbetween.NA),colnames(MAT))
  if(length(id.inbetween.NA)>0)
  {
    a=MAT[,id.inbetween.NA]
    if(is.vector(a)) a=na.approx(a)
    if(!is.vector(a)) for(e in 1:ncol(a))if(!sum(is.na(a[,e]))==nrow(a))a[,e]=na.approx(a[,e])
    MAT[,id.inbetween.NA]=a 
  }
  
  #Combine data
  Imputed.DAT=expand.grid(FINYEAR=yr,BLOCKX=blk)
  Imputed.DAT=Imputed.DAT[order(Imputed.DAT$BLOCKX,Imputed.DAT$FINYEAR),]
  Imputed.DAT$Pred.Catch.imp=c(MAT)  
  DAT=merge(DAT,Imputed.DAT,by=c("FINYEAR","BLOCKX"),all.x=T)
  DAT$Pred.Catch.imp=with(DAT,ifelse(is.na(Pred.Catch),Pred.Catch.imp,Pred.Catch))
  return(DAT)
}

fn.linear.imp=function(x,y,SE,pred,cap)  #linear extrapolation
{
  daT=data.frame(x=x,y=y,SE=SE)
  mod=lm(y~x,data=daT,weights=1/(SE)^2)       
  PRDS=predict(mod,newdata=data.frame(x=pred))
  PRDS=ifelse(PRDS>cap,cap,PRDS)   #cap extrapolation at max ever observed
  if(as.numeric(coef(mod)[2])<0)PRDS=rep(y[length(y)],length(pred)) # if extrapolated trend is negative, as there is no fishing, reset to last year of data
  return(PRDS)
}

fn.pop.grw.rte.imp=function(x,y,pred,cap)  #population growth rate
{
  Strt=mean(y)
  PRDS=Strt+Pop.growth*(pred-x[length(x)])
  PRDS=ifelse(PRDS>cap,cap,PRDS)   #cap extrapolation at max ever observed    
  return(PRDS)
}

construct.index=function(LISTA,Area.w,IMPUTE,SPEC,Weight.by.area) #function for constructing index                        
{
  #Get data
  Bi=LISTA$Bi.dat
  Pos=LISTA$Pos.dat
  if(!is.null(Bi))Dat=merge(Pos,Bi[,match(c("FINYEAR","BLOCKX","Pred.Prob"),names(Bi))],by=c("FINYEAR","BLOCKX"),all.x=T)
  if(is.null(Bi))Dat=Pos
  
  #Area weight as a proportion of total habitat area
  Area.w=subset(Area.w,BLOCKX%in%unique(Dat$BLOCKX))
  #Area.w$Fish.Area=Area.w$Fish.Area/sum(Area.w$Fish.Area)  
  
  #Convert to normal space                                  
  #if(IMPUTE=="YES")Dat$KTCH=exp(Dat$Pred.Catch.imp)
  #if(IMPUTE=="NO")Dat$KTCH=exp(Dat$Pred.Catch)
  #bias corrected
  MEAN.SE=mean(Dat$SE,na.rm=T)
  Dat$SE=with(Dat,ifelse(is.na(SE),MEAN.SE,SE))
  if(IMPUTE=="YES") Dat$KTCH=BiasCor.fn(Dat$Pred.Catch.imp,Dat$SE)
  if(IMPUTE=="NO") Dat$KTCH=BiasCor.fn(Dat$Pred.Catch,Dat$SE)
  
  
  #Multiply Positive catch and Prob of catch
  if(!is.null(Bi))Dat$I_y.b=Dat$KTCH*Dat$Pred.Prob
  if(is.null(Bi))Dat$I_y.b=Dat$KTCH
  
  #Add block weight
  names(Area.w)=c("BLOCKX","WEIGHT","MaxDepth")
  Dat=merge(Dat,Area.w[,match(c("BLOCKX","WEIGHT"),names(Area.w))],by="BLOCKX",all.x=T)
  
  #Calculate annual index 
  if(Weight.by.area=="YES")
  {
    #Multiply index by weighted (fishable area)
    Dat$I_y.b.W=Dat$I_y.b*Dat$WEIGHT
    
    #Aggregate index by year
    I_y=aggregate(I_y.b.W~FINYEAR,Dat,sum,na.rm=T)
  }
  
  if(Weight.by.area=="NO")
  {
    Dat$I_y.b.W=Dat$I_y.b
    
    #Aggregate index by year
    I_y=aggregate(I_y.b.W~FINYEAR,Dat,mean,na.rm=T)
  }
  
  #express in selected effort units
  #I_y$I_y.b.W=I_y$I_y.b.W/Eff.unit  #used only if fitting to catch
  
  return(I_y=I_y)
}

construct.index.month=function(LISTA,Area.w,IMPUTE,SPEC,Weight.by.area)    #function for constructing index by month                      
{
  #Get data
  Bi=LISTA$Bi.dat
  Pos=LISTA$Pos.dat
  if(!is.null(Bi))Dat=merge(Pos,Bi[,match(c("FINYEAR","BLOCKX","Pred.Prob"),names(Bi))],by=c("FINYEAR","BLOCKX"),all.x=T)
  if(is.null(Bi))Dat=Pos
  
  #Area weight as a proportion of total habitat area
  Area.w=subset(Area.w,BLOCKX%in%unique(Dat$BLOCKX))
  Area.w$Fish.Area=Area.w$Fish.Area/sum(Area.w$Fish.Area)  
  
  #Convert to normal space                                  
  #if(IMPUTE=="YES")Dat$KTCH=exp(Dat$Pred.Catch.imp)
  #if(IMPUTE=="NO")Dat$KTCH=exp(Dat$Pred.Catch)
  #bias corrected
  MEAN.SE=mean(Dat$SE,na.rm=T)
  Dat$SE=with(Dat,ifelse(is.na(SE),MEAN.SE,SE))  
  if(IMPUTE=="YES") Dat$KTCH=BiasCor.fn(Dat$Pred.Catch.imp,Dat$SE)
  if(IMPUTE=="NO") Dat$KTCH=BiasCor.fn(Dat$Pred.Catch,Dat$SE)
  
  #Multiply Positive catch and Prob of catch
  if(!is.null(Bi))Dat$I_y.b=Dat$KTCH*Dat$Pred.Prob
  if(is.null(Bi))Dat$I_y.b=Dat$KTCH
  
  #Add block weight
  names(Area.w)=c("BLOCKX","WEIGHT","MaxDepth")
  Dat=merge(Dat,Area.w[,match(c("BLOCKX","WEIGHT"),names(Area.w))],by="BLOCKX",all.x=T)
  
  #Calculate monthly index 
  if(Weight.by.area=="YES")
  {
    #Multiply index by weighted (fishable area)
    Dat$I_y.b.W=Dat$I_y.b*Dat$WEIGHT
    
    #Aggregate by month
    I_y=aggregate(I_y.b.W~MONTH+FINYEAR,Dat,sum,na.rm=T)
  }
  
  if(Weight.by.area=="NO")
  {
    Dat$I_y.b.W=Dat$I_y.b
    
    #Aggregate by month
    I_y=aggregate(I_y.b.W~MONTH+FINYEAR,Dat,mean,na.rm=T)
  }
  
  return(I_y=I_y)
}

gm_mean = function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) #geometric mean function

Plot.Comb.cpue=function(BC,Fan,No.creep,Ind.ves.sens,No.wght,SPEC,REL,WHERE,Ves.no.ktch,Min.w)  #function for plotting sensitivity analysis
{
  a=get.foly.and.nominal(List.foly.nom[[i]])
  Nominal=a$Nominal
  names(Nominal)[2]="Nominal"
  Foly=a$Foly
  Unstandardised=cbind(Nominal,Foly[,2])
  ID=which(!Unstandardised$FINYEAR%in%BC$FINYEAR)
  if(length(ID)>0) Unstandardised[ID,2:3]=NA
  
  ID=which(!BC$FINYEAR%in%Ind.ves.sens$FINYEAR)
  if(length(ID)>0)
  {
    ADD=BC[ID,]
    ADD$I_y.b.W=NA
    Ind.ves.sens=rbind(Ind.ves.sens,ADD)
    Ind.ves.sens=Ind.ves.sens[order(Ind.ves.sens$FINYEAR),]
  }
  names(BC)[2]="BaseCase"
  names(Fan)[2]="Fan"   
  names(No.creep)[2]="No.creep"
  names(Ind.ves.sens)[2]="Ind.ves.sens"
  names(No.wght)[2]="No.wght"
  
  Standardised=cbind(Fan,No.creep[,2],Ind.ves.sens[,2],No.wght[,2],BC[,2])
  Series=merge(Unstandardised,Standardised,by="FINYEAR",all.x=T)
  Series=subset(Series,!is.na(Nominal))
  YRs=Series$FINYEAR
  N.YRs=length(YRs)
  
  #normalise  
  if(REL=="YES")if(relative.to=="Mean.1")
  {
    for(t in 2:ncol(Series))
    {
      MN=mean(Series[,t],na.rm=T)
      if(t==5) MN=mean(Series[,8],na.rm=T)
      Series[,t]=Series[,t]/MN 
    }
  }
  
  MAXi=max(Series[,2:ncol(Series)],na.rm=T)
  
  plot(1:N.YRs,Series[,2],ylab="",type='l',xlab="",xaxt='n',col=COLORES[1],
       lty=LINEAS[1],lwd=LWD[1],ylim=c(0,MAXi),cex.axis=1.35)
  for(t in 2:(ncol(Series)-1)) lines(1:N.YRs,Series[,t+1],lwd=LWD[t],lty=LINEAS[t],col=COLORES[t]) 
  axis(1,at=1:N.YRs,labels=F,tck=-0.01)
  #legend(WHERE,SPEC,bty='n',cex=1.75)
  axis(1,at=seq(1,N.YRs,by=5),labels=F,tck=-0.01)
  axis(1,at=seq(1,N.YRs,by=5),labels=YRs[seq(1,N.YRs,by=5)],tck=-0.02,cex.axis=1.15,padj=-.05)
}

Plot.Comb.cpue.daily=function(BC,Fan,BC.daily,Fan.daily,REL)  #function for plotting sensitivity analysis
{
  a=get.foly.and.nominal(List.foly.nom[[i]])
  Nominal=a$Nominal
  names(Nominal)[2]="Nominal"
  Foly=a$Foly
  names(Foly)[2]="Folly"
  YRs=Nominal$FINYEAR
  N.YRs=length(YRs)
  
  names(BC)[2]="BaseCase"
  names(Fan)[2]="Fan"
  names(BC.daily)[2]="BaseCase.daily"
  names(Fan.daily)[2]="Fan.daily"
  
  #normalise  
  if(REL=="YES")if(relative.to=="Mean.1")
  {
    Nominal[,2]=Nominal[,2]/mean(Nominal[,2])
    Foly[,2]=Foly[,2]/mean(Foly[,2])
    BC[,2]=BC[,2]/mean(BC[,2])
    Fan[,2]=Fan[,2]/mean(Fan[,2])
    BC.daily[,2]=BC.daily[,2]/mean(BC.daily[,2])
    Fan.daily[,2]=Fan.daily[,2]/mean(Fan.daily[,2])
  }
  
  #combine series
  Unstandardised=merge(Nominal,Foly,by="FINYEAR")     
  Standardised=merge(Fan,BC,by="FINYEAR")
  Standardised.daily=merge(Fan.daily,BC.daily,by="FINYEAR")
  
  Series=merge(Unstandardised,Standardised,by="FINYEAR",all.x=T)
  Series=merge(Series,Standardised.daily,by="FINYEAR",all.x=T)
  Series=subset(Series,!is.na(Nominal))
  
  MAXi=max(Series[,2:ncol(Series)],na.rm=T)
  
  plot(1:N.YRs,Series[,2],ylab="",type='l',xlab="",xaxt='n',col=COLORES[1],
       lty=LINEAS[1],lwd=LWD[1],ylim=c(0,MAXi),cex.axis=1.35)
  polygon(x=c(-0.5,match("2005-06",YRs),match("2005-06",YRs),-0.5),
          y=c(-0.5,-0.5,MAXi*2,MAXi*2),col="grey75",border="grey75")
  polygon(x=c(match("2005-06",YRs),length(YRs)*2,length(YRs)*2,match("2005-06",YRs)),
          y=c(-0.5,-0.5,MAXi*2,MAXi*2),col="grey95",border="grey95")
  box()
  for(t in 1:(ncol(Series)-1)) lines(1:N.YRs,Series[,t+1],lwd=LWD[t],lty=LINEAS[t],col=COLORES[t]) 
  axis(1,at=1:N.YRs,labels=F,tck=-0.01)
  axis(1,at=seq(1,N.YRs,by=5),labels=F,tck=-0.01)
  axis(1,at=seq(1,N.YRs,by=5),labels=YRs[seq(1,N.YRs,by=5)],tck=-0.02,cex.axis=1.15,padj=-.05)
}

Plot.Comb.cpue.eff.creep=function(BC,No.creep,SPEC,WHERE)  #function for showing effect of creep
{
  names(BC)[2]="BaseCase"
  names(No.creep)[2]="No.creep"
  Series=cbind(No.creep,BC[,2])
  YRs=Series$FINYEAR
  N.YRs=length(YRs)
  MAXi=max(Series[,2:ncol(Series)],na.rm=T)
  
  plot(1:N.YRs,Series[,2],ylab="",type='l',xlab="",xaxt='n',col=COLORES[1],
       lty=LINEAS[1],lwd=LWD,ylim=c(0,MAXi),cex.axis=1.75)
  for(t in 2:(ncol(Series)-1)) lines(1:N.YRs,Series[,t+1],lwd=LWD,lty=LINEAS[t],col=COLORES[t]) 
  axis(1,at=1:N.YRs,labels=F,tck=-0.01)
  legend(WHERE,SPEC,bty='n',cex=2)
  axis(1,at=seq(1,N.YRs,by=1),labels=F,tck=-0.02)
  axis(1,at=seq(1,N.YRs,by=5),labels=YRs[seq(1,N.YRs,by=5)],tck=-0.04,cex.axis=1.5,padj=-.05)
}

get.foly.and.nominal=function(DATA)    #extract foly and nominal cpue
{
  return(list(Nominal=DATA$Nom,Foly=DATA$Foly))
}

fn.prop.imputed=function(DAT,DAT1,SPEC,MOD)   #function for plotting impluted blocks
{
  DAT=DAT[,match(c("FINYEAR","BLOCKX"),names(DAT))]
  DAT$count=NA
  DAT$dummy=with(DAT,paste(FINYEAR,BLOCKX,sep=""))
  
  DAT1$count=1
  DAT1$dummy=with(DAT1,paste(FINYEAR,BLOCKX,sep=""))
  
  DAT=subset(DAT,!dummy%in%DAT1$dummy)
  DAT=rbind(DAT,DAT1)
  
  DAT=DAT[,-match("dummy",names(DAT))]
  
  #add block coefficient
  ID=c(1,grep("[:]", names(coef(MOD))))
  Cofs=coef(MOD)[-ID]
  Store1=as.character(sort(unique(DAT[,match("BLOCKX",names(DAT))]))[-1])
  Store2=paste("BLOCKX",Store1,sep="")
  MatcH=Cofs[match(Store2,names(Cofs))]
  names(MatcH)=substr(names(MatcH),start=7,stop=15)
  ADd=data.frame(BLOCKX=names(MatcH),coef=MatcH)
  
  DAT=merge(DAT,ADd,by="BLOCKX",all.x=T)
  DAT$coef=with(DAT,ifelse(is.na(coef),0,coef))
  DAT$count=with(DAT,ifelse(!is.na(count),coef,count))
  DAT=subset(DAT,!coef==0)
  DAT=DAT[,-match("coef",names(DAT))]
  dat=reshape(DAT,idvar=c("FINYEAR"),timevar="BLOCKX",v.names="count", direction="wide")
  names(dat)[2:ncol(dat)]=substr(names(dat)[2:ncol(dat)],start=7,stop=15)
  dat=dat[order(dat$FINYEAR),]
  Breaks=sort(quantile(unlist(as.matrix(dat[,2:ncol(dat)])),probs=seq(0,1,1/numInt),na.rm=T))
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  image(1:nrow(dat),2:ncol(dat),as.matrix(dat[,2:ncol(dat)]),col =couleurs,breaks=Breaks,xaxt='n',
        yaxt='n',ylab="",xlab="",main=SPEC,cex.main=1.5)
  axis(1,1:nrow(dat),F,tck=-0.015)
  axis(1,seq(1,nrow(dat),5),dat$FINYEAR[seq(1,nrow(dat),5)],tck=-0.03,cex.axis=1.25)
  axis(2,2:ncol(dat),colnames(dat)[2:ncol(dat)],las=1,tck=-0.03,cex.axis=1.05)
  color.legend(nrow(dat)*1.175,ncol(dat),nrow(dat)*1.275,ncol(dat)*.60,round(Breaks,2),
               rect.col=couleurs,gradient="y",col=colLeg,cex=0.95)
  box()
  
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
    points(xo[s,],yo[s,],cex=zo[,s],pch=16,col="grey60")
    points(xo[s,],yo[s,],cex=zo[,s],pch=1,col="black")
  }
}
Influence.fn=function(MOD,DAT,DAT2,Term.type,termS,LABEL,XLIM,SCALE)  #Bentley et al 2012
{
  DAT=DAT[,match(VARIABLES,names(DAT))]
  
  DAT$log.Catch.Gummy=with(DAT,ifelse(Catch.Gummy>0,log(Catch.Gummy),0))
  DAT$log.Catch.Whiskery=with(DAT,ifelse(Catch.Whiskery>0,log(Catch.Whiskery),0))
  DAT$log.Catch.Dusky=with(DAT,ifelse(Catch.Dusky>0,log(Catch.Dusky),0))
  DAT$log.Catch.Sandbar=with(DAT,ifelse(Catch.Sandbar>0,log(Catch.Sandbar),0))
  
  
  #extract main term coefficients for each species
  Cat.term=which(Term.type=="CAT")
  Term.type=Term.type[Cat.term] #only do categorical vars
  termS=termS[Cat.term]
  
  nt=length(termS)
  Store1=Store2=MatcH=COEF.list=COEF.SE.list=vector('list',nt)
  ID=c(1,grep("[:]", names(coef(MOD))))
  Cofs=coef(MOD)[-ID]
  Cofs.SE=summary(MOD)$coefficients[-ID, 2]
  
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      Store1[[p]]=as.character(levels(DAT2[,match(termS[p],names(DAT))]))[-1]
      #Store1[[p]]=as.character(sort(unique(DAT[,match(termS[p],names(DAT))]))[-1])
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
      #A=as.character(sort(unique(DAT[,match(termS[p],names(DAT))]))[1])
      A=as.character(levels(DAT2[,match(termS[p],names(DAT))]))[1]
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
      DAT=cbind(DAT,rep(COEF.list[[p]],nrow(DAT)))
      colnames(DAT)[ncol(DAT)]=paste("Coef.",termS[p],sep="")
    }
  }
  
  
  #ny
  ny=table(DAT$FINYEAR)
  
  #calculate mean of coefficient
  Coef.vec=match(paste("Coef.",termS,sep=""),names(DAT))
  Mean.coef=Annual.Dev=vector('list',nt)
  for(p in 1:nt) Mean.coef[[p]]=mean(DAT[,Coef.vec[p]])
  
  
  
  #calculate annual deviation from mean
  #- Categorical variables
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      dev=rep(NA,length(ny))
      for(t in 1:length(ny))
      {
        a=subset(DAT,FINYEAR==names(ny[t]))
        dev[t]=(sum(a[,Coef.vec[p]]-Mean.coef[[p]]))/ny[t]
      }
      
      #exp because it's multiplicative
      Annual.Dev[[p]]=exp(dev)  
    }
  }
  
  Over.all.influence=rep(NA,nt)
  names(Over.all.influence)=termS
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
      
      tiff(paste(SPECIES.vec[i],".CDI.",termS[p],".tiff",sep=""),width = 2400, height = 1600,units = "px", res = 300, compression = "lzw")
      nf <- layout(matrix(c(1,1,0,2,2,3), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
      par(mar=c(0,0,0,0),oma=c(4,6,1,1),las=1,mgp=c(1,.9,0))
      layout.show(nf)
      # Coefficients
      x=1:length(COEF)
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
      
      #Influence plot
      plot(Annual.Dev[[p]],1:length(ny),type="o",pch=19,xlab="",ylab="",cex=2,cex.axis=1.25,yaxt='n',xlim=XLIM)
      abline(v=1,lty=3,col=1)
      mtext("Influence",side=1,line=2.5,cex=1.5)      
      axis(2,1:length(ny),F,tcl=0.5)
      axis(2,seq(1,length(ny),2),F,tcl=1)
      dev.off()
      
      #Store overall influence of variable
      Over.all.influence[p]=exp(sum(abs(log(Annual.Dev[[p]])))/length(ny))-1
    }
    
  }
  
  
  return(list(Annual.Dev=Annual.Dev,ny=ny,Over.all.influence=Over.all.influence))
}

#functions used in bootstrapping
fn.build.pred.mn.dat.boot=function(dat)  #function for creating standardise data for predictions
{
  BIN=dat$bin
  yr=unique(BIN$FINYEAR)
  blk=unique(BIN$BLOCKX)
  Mn=All.mns
  Ves=unique(BIN$VESSEL)
  Gum=unique(BIN$log.Catch.Gummy)
  Whis=unique(BIN$log.Catch.Whiskery)
  Dus=unique(BIN$log.Catch.Dusky)
  San=unique(BIN$log.Catch.Sandbar)
  
  log.Effort=unique(BIN$log.Effort)
  Freo=unique(BIN$Freo)
  SOI=unique(BIN$SOI)
  BIN=expand.grid(FINYEAR=yr,BLOCKX=blk,MONTH=Mn,VESSEL=Ves,log.Catch.Gummy=Gum,
                  log.Catch.Whiskery=Whis,log.Catch.Dusky=Dus,log.Catch.Sandbar=San,
                  log.Effort=log.Effort,Freo=Freo,SOI=SOI)
  if(Bi.same.Pos=="YES") POS=BIN
  return(list(bin=BIN,pos=POS))
}
Boot.fn.res=function(N.boot)   #function for resampling residuals
{
  LOg=MOD$LogN
  PosData=DAT 
  
  #1. Predict values
  PosData$LOg.pred=predict(LOg,type='response')
  
  #2. Random sample of residuals and add to predicted value
  if(Fit.to.what=="catch")
  {
    if(Res.type=="additive") Residuals=PosData$log.Catch-PosData$LOg.pred   #(Haddon 2000 page 306)
    if(Res.type=="multiplicative") Residuals=PosData$log.Catch/PosData$LOg.pred   
  }
  if(Fit.to.what=="cpue")
  {
    if(Res.type=="additive") Residuals=PosData$log.cpue-PosData$LOg.pred   #(Haddon 2000 page 306)
    if(Res.type=="multiplicative") Residuals=PosData$log.cpue/PosData$LOg.pred   
  }
  
  NEW.dat=vector('list',N.boot)
  for ( b in 1:N.boot)
  {
    NEW=PosData
    LOg.res=sample(Residuals,replace=T) 
    if(Fit.to.what=="catch")
    {
      if(Res.type=="additive") NEW$log.Catch=NEW$LOg.pred+LOg.res   
      if(Res.type=="multiplicative") NEW$log.Catch=NEW$LOg.pred*LOg.res   
    }
    
    if(Fit.to.what=="cpue")
    {
      if(Res.type=="additive") NEW$log.cpue=NEW$LOg.pred+LOg.res   
      if(Res.type=="multiplicative") NEW$log.cpue=NEW$LOg.pred*LOg.res   
    }
    NEW.dat[[b]]=NEW
  }
  return(PosData=NEW.dat)
}
Boot.fn=function(DAT,VES.No.K,vars,SPATIAL)  #function for resampling data
{
  #1. drop vessels that never reported catch of target species   
  DAT=subset(DAT,!(VESSEL%in%VES.No.K))
  
  #2. put data in proper shape
  DataFile=DAT[,match(vars,names(DAT))]
  
  if(SPATIAL=="Zones.combined") DataFile$YrBlk=with(DataFile,paste(FINYEAR,BLOCKX))
  if(SPATIAL=="By.zone") DataFile$YrBlk=with(DataFile,paste(FINYEAR,BLOCKX,MONTH))
  YrBLK=unique(DataFile$YrBlk)
  
  #3. resample with replacement
  Lista=vector('list',length(YrBLK))
  for (x in 1:length(YrBLK))
  {
    dat=subset(DataFile,YrBlk==YrBLK[x])
    ID=sample(1:nrow(dat),replace=T)
    dat=dat[ID,]
    Lista[[x]]=dat
  }
  
  New.dat=do.call(rbind,Lista)
  return(New.dat=New.dat)
}
fn.bootstrap.index=function(DAT,Spc,vrs,best,Grid.d,Grid.d.zone,Ves.no,area.w,zonas,bloques) #run whole boot procedure
{
  #1. Set storing objects
  Store.boot.BaseCase=vector('list',N.boot)
  Store.boot.BaseCase.zone=Store.boot.BaseCase.annual.by.zone=Store.boot.BaseCase
  #2. Run bootstrap
  for (n in 1:N.boot) 
  {
    tryCatch(
{
  #2.1 create bootstrapped data
  if(Boot.method=="DATA") Boot.Dat=Boot.fn(DAT,Ves.no,vrs,"Zones.combined")
  if(Boot.method=="RESIDUALS") Boot.Dat=Bootstrapped.Dat[[i]][[n]]
  
  #2.2. refit GLMs to bootstrapped data
  if(Boot.method=="DATA") Fit.BaseCase=fn.stand.cpue(Boot.Dat,Ves.no,Spc,best$Bi,best$Log,vrs,Eff.Creep="YES")      
  if(Boot.method=="RESIDUALS") Fit.BaseCase=fn.stand.cpue.boot.res(BiData,Boot.Dat,best$Bi,best$Log)      
  
  #2.3. run predictions
  #zones combined
  Pred.BaseCase=predict.ktch(Fit.BaseCase,Grid.d,Spc)
  
  #by zone and year
  Nz=length(zonas)
  dummy=vector('list',Nz)
  names(dummy)=zonas
  Pred.BaseCase.annual.by.zone=Pred.BaseCase.zone=dummy
  
  for(z in 1:Nz)
  {
    #keep only relevant blocks for each zone
    a=Grid.d
    this.zn=zonas[z]
    these.blks=BLKZ.per.zone[[match(this.zn,names(BLKZ.per.zone))]]
    a$bin=subset(a$bin,BLOCKX%in%these.blks)
    a$pos=subset(a$pos,BLOCKX%in%these.blks)
    
    #make prediction
    Pred.BaseCase.annual.by.zone[[z]]=predict.ktch(Fit.BaseCase,a,Spc)
  }
  
  
  #by zone and month
  for(z in 1:Nz)
  {
    #keep only relevant blocks for each zone
    a=Grid.d.zone
    this.zn=zonas[z]
    these.blks=BLKZ.per.zone[[match(this.zn,names(BLKZ.per.zone))]]
    a$bin=subset(a$bin,BLOCKX%in%these.blks)
    a$pos=subset(a$pos,BLOCKX%in%these.blks)
    
    #make prediction
    Pred.BaseCase.zone[[z]]=predict.ktch(Fit.BaseCase,a,Spc)
  }
  
  #2.4. run imputations      
  #zones combined
  Pred.BaseCase$Pos.dat=impute.ktch(Pred.BaseCase$Pos.dat)  
  #by zone and year
  for(z in 1:Nz) Pred.BaseCase.annual.by.zone[[z]]$Pos.dat=impute.ktch(Pred.BaseCase.annual.by.zone[[z]]$Pos.dat) 
  #by zone
  for(z in 1:Nz) Pred.BaseCase.zone[[z]]$Pos.dat=impute.ktch(Pred.BaseCase.zone[[z]]$Pos.dat) 
  
  #2.5. construct abundance indices for each zone or combined
  #zones combined
  Store.boot=construct.index(Pred.BaseCase,area.w,"YES",Spc,Weight.by.area="YES")
  
  #by zone and year
  for(z in 1:Nz)
  {
    BZ=as.numeric(unique(Pred.BaseCase.annual.by.zone[[z]]$Pos.dat$BLOCKX))
    dummy[[z]]=construct.index(Pred.BaseCase.annual.by.zone[[z]],subset(area.w,BLOCKX%in%BZ),"YES",Spc,Weight.by.area="YES")
  }   
  Store.boot.annual.by.zone=dummy
  
  #by zone and month
  for(z in 1:Nz)
  {
    BZ=as.numeric(unique(Pred.BaseCase.zone[[z]]$Pos.dat$BLOCKX))
    dummy[[z]]=construct.index.month(Pred.BaseCase.zone[[z]],subset(area.w,BLOCKX%in%BZ),"YES",Spc,Weight.by.area="YES")
  }
  Store.Zone.boot=dummy
  
  #2.6. store bootstrapped index
  Store.boot.BaseCase[[n]]=Store.boot
  Store.boot.BaseCase.annual.by.zone[[n]]=Store.boot.annual.by.zone
  Store.boot.BaseCase.zone[[n]]=Store.Zone.boot
  
  #2.7. free up memory
  rm(Boot.Dat,Fit.BaseCase,Pred.BaseCase,Pred.BaseCase.zone)
  
}, error=function(e){})
  }

return(list(Store.boot.BaseCase=Store.boot.BaseCase,
            Store.boot.BaseCase.annual.by.zone=Store.boot.BaseCase.annual.by.zone,
            Store.boot.BaseCase.zone=Store.boot.BaseCase.zone))
}
fn.blk.zn=function(DAT)  #function for extracting blocks from zones
{
  zn=sort(unique(DAT$zone))
  LISTA=vector('list',length(zn))
  names(LISTA)=zn
  for(n in 1:length(zn))
  {
    x=subset(DAT,zone==zn[n])
    LISTA[[n]]=sort(unique(x$BLOCKX))
  }
  return(LISTA) 
}
fn.stand.cpue.boot.res=function(Bi,Pos,Formula.cpue,Formula.cpue.pos) #function for running glm of residuals boot
{
  #Binomial GLM
  Log.GLMbi=NULL
  if(!is.null(Formula.cpue))Log.GLMbi <- glm(Formula.cpue, data=Bi, family="binomial", maxit=500)
  
  #Lognormal GLM
  GLMlog <- glm(Formula.cpue.pos, data=Pos, family=gaussian, maxit=500)
  
  return(list(Bi=Log.GLMbi,LogN=GLMlog))
}
fn.get.mean=function(DAT,what)  #function for extracting mean, CV and confidence intervals     
{
  #Annual
  if(what=="ZonesCombined")
  {
    YRs=DAT[[1]]$FINYEAR
    N.YRs=length(YRs)
    N.rep=length(DAT)
    BC=matrix(NA,nrow=N.YRs,ncol=N.rep)    
    for(o in 1:N.rep) if(!is.null(DAT[[o]]$I_y.b.W)) BC[,o]=DAT[[o]]$I_y.b.W
    
    #mean
    MEAN=rowMeans(BC,na.rm=T)
    
    #SD
    SD=apply(BC,1,sd,na.rm=T)
    
    #CV   
    CV=SD/MEAN
    
    #CI
    LOW=apply(BC, 1, function(x) quantile(x, 0.025,na.rm=T))
    UP=apply(BC, 1, function(x) quantile(x, 0.975,na.rm=T))
    
    return(list(MEAN=MEAN,LOW.CI=LOW,UP.CI=UP,FINYEAR=YRs,CV=CV,SD=SD))       
    
  }
  
  #Annual by zone
  if(what=="AnnualByZones")
  {
    YR.Mns=DAT[[1]][[1]]$FINYEAR
    N.YRs=length(YR.Mns)
    N.rep=length(DAT)
    zns=names(DAT[[1]])
    MEAN=vector('list',length(zns))    
    names(MEAN)=zns
    LOW=UP=SD=CV=MEAN
    for(z in 1:length(zns)) 
    {
      BC=matrix(NA,nrow=N.YRs,ncol=N.rep)
      for(o in 1:N.rep)
      {
        dummy=DAT[[o]][[z]]
        if(!is.null(dummy))
        {
          dummy=dummy[order(dummy$FINYEAR),] #order by year-month
          BC[,o]=dummy$I_y.b.W
        }
      }
      
      MEAN[[z]]=rowMeans(BC,na.rm=T) #mean     
      
      #SD
      SD[[z]]=apply(BC,1,sd,na.rm=T)
      
      #CV   
      CV[[z]]=SD[[z]]/MEAN[[z]]
      
      LOW[[z]]=apply(BC, 1, function(x) quantile(x, 0.025,na.rm=T)) #CI
      UP[[z]]=apply(BC, 1, function(x) quantile(x, 0.975,na.rm=T))
    }
    return(list(MEAN=MEAN,LOW.CI=LOW,UP.CI=UP,SD=SD,CV=CV,FINYEAR=YR.Mns))    
  }
  
  
  #Monthly by zone
  if(what=="ByZones")
  {
    YR.Mns=DAT[[1]][[1]]$FINYEAR
    N.YRs=length(YR.Mns)
    N.rep=length(DAT)
    zns=names(DAT[[1]])
    MEAN=vector('list',length(zns))    
    names(MEAN)=zns
    LOW=UP=SD=CV=MEAN
    for(z in 1:length(zns)) 
    {
      BC=matrix(NA,nrow=N.YRs,ncol=N.rep)
      for(o in 1:N.rep)
      {
        dummy=DAT[[o]][[z]]
        if(!is.null(dummy))
        {
          dummy=dummy[order(dummy$FINYEAR,dummy$MONTH),] #order by year-month
          BC[,o]=dummy$I_y.b.W
        }
      }
      
      MEAN[[z]]=rowMeans(BC,na.rm=T) #mean     
      
      #SD
      SD[[z]]=apply(BC,1,sd,na.rm=T)
      
      #CV   
      CV[[z]]=SD[[z]]/MEAN[[z]]
      
      LOW[[z]]=apply(BC, 1, function(x) quantile(x, 0.025,na.rm=T)) #CI
      UP[[z]]=apply(BC, 1, function(x) quantile(x, 0.975,na.rm=T))
    }
    return(list(MEAN=MEAN,LOW.CI=LOW,UP.CI=UP,SD=SD,CV=CV,FINYEAR=YR.Mns))    
  }
}

#function for plotting index
Plot.Index=function(MEAN,LOW1,UP1,FINYEAR,COL,Colr,Colr2,RELATIVE,Time.step,Do.axis)
{
  if(RELATIVE=="YES")
  { 
    Mn=mean(MEAN)
    MEAN=MEAN/Mn
    LOW1=LOW1/Mn
    UP1=UP1/Mn
  }
  if(RELATIVE=="NO")
  {
    MEAN=MEAN
    LOW1=LOW1
    UP1=UP1
  }
  MAX.y=max(UP1)
  N.YRs=length(FINYEAR)
  YR=1:N.YRs
  
  #polygons
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec <- c(LOW1, tail(UP1, 1), rev(UP1), LOW1[1])
  
  #plot  
  plot(YR,MEAN,type='l',lty=1,col=COL,lwd=3,ylim=c(0,MAX.y),
       xaxt='n',xlab="",ylab="",cex.axis=1.25)
  polygon(Year.Vec, Biom.Vec, col = Colr, border = Colr2,lwd=1.25)
  lines(YR,MEAN,type='l',lty=1,col=COL,lwd=3)
  if(Time.step=="YEAR")
  {
    axis(1,YR,labels=F,tck=-0.015)
    if(Do.axis=="NO")axis(1,seq(1,N.YRs,5),labels=F,tck=-0.030,cex.axis=1.25)
    if(Do.axis=="YES")axis(1,seq(1,N.YRs,5),labels=FINYEAR[seq(1,N.YRs,5)],tck=-0.030,cex.axis=1.25)
  }
  if(Time.step=="YR.MN")
  {
    step=seq(0,N.YRs,by=12)
    axis(1,step,labels=F,tck=-0.030)
    step=seq(0,N.YRs,by=12*5)
    if(Do.axis=="NO")axis(1,step,labels=F,tck=-0.060)
    if(Do.axis=="YES")axis(1,step,labels=FINYEAR[1+step],tck=-0.060,cex.axis=1)
  }
  
}
fn.out.boot.zn=function(a,NAME)  #function for exporting all bootstrapped indices
{
  YRs=a[[1]][[1]]$FINYEAR
  N.YRs=length(YRs)
  N.rep=length(a)
  Zns=names(a[[1]])
  for(p in 1:length(Zns))
  {
    BC=matrix(NA,nrow=N.YRs,ncol=N.rep)
    for(o in 1:N.rep) BC[,o]=a[[o]][[p]]$I_y.b.W
    test=data.frame(YRs,BC)
    names(test)=NMs
    nm=paste(SPECIES.vec[i],".",Zns[p],NAME,sep="")
    write.csv(test,nm,row.names=F)
  }
}
export.foly=function(DATA,DATA1)     #function for creating foly indices
{
  #monthly
  Agg.ktc=aggregate(LIVEWT.c~FINYEAR,DATA,sum,na.rm=T)
  Agg.eff=aggregate(Km.Gillnet.Days.c~FINYEAR,DATA,sum,na.rm=T)
  Foly=Agg.ktc$LIVEWT.c/Agg.eff$Km.Gillnet.Days.c
  
  #daily
  Agg.ktc_d=aggregate(LIVEWT.c~FINYEAR,DATA1,sum,na.rm=T)
  Agg.eff_d=aggregate(Km.Gillnet.Days.c~FINYEAR,DATA1,sum,na.rm=T)
  Foly_d=Agg.ktc_d$LIVEWT.c/Agg.eff_d$Km.Gillnet.Days.c
  
  
  Foly=data.frame(FINYEAR=c(Agg.ktc$FINYEAR,Agg.ktc_d$FINYEAR),CPUE=c(Foly,Foly_d))
  return(Foly)
}
export.nominal=function(DATA,DATA1)     #function for creating nominal indices 
{
  DATA$cpue=DATA$LIVEWT.c/DATA$Km.Gillnet.Days.c
  DATA1$cpue=DATA1$LIVEWT.c/DATA1$Km.Gillnet.Days.c
  Nominal=aggregate(cpue~FINYEAR,DATA,mean)
  Nominal_d=aggregate(cpue~FINYEAR,DATA1,mean)
  return(rbind(Nominal,Nominal_d))
}



#----4. PROCEDURE SECTION-----#

#4.1. Add SOI and Freo sea level data

#Freo lag
Freo$Freo.Lag6=c(rep(NA,6),Freo$Freo[1:(length(Freo$Freo)-6)])
Freo$Freo.Lag12=c(rep(NA,12),Freo$Freo[1:(length(Freo$Freo)-12)])

#monthly
Data.monthly.GN.whiskery=merge(Data.monthly.GN.whiskery,SOI, by.x=c("MONTH","YEAR.c"),
                               by.y=c("Month.soi","Year.soi"),all.x=T)
Data.monthly.GN.gummy=merge(Data.monthly.GN.gummy,SOI, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month.soi","Year.soi"),all.x=T)  
Data.monthly.GN.dusky=merge(Data.monthly.GN.dusky,SOI, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month.soi","Year.soi"),all.x=T)
Data.monthly.GN.sandbar=merge(Data.monthly.GN.sandbar,SOI, by.x=c("MONTH","YEAR.c"),
                              by.y=c("Month.soi","Year.soi"),all.x=T)

Data.monthly.GN.whiskery=merge(Data.monthly.GN.whiskery,Freo, by.x=c("MONTH","YEAR.c"),
                               by.y=c("Month","Year"),all.x=T)
Data.monthly.GN.gummy=merge(Data.monthly.GN.gummy,Freo, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month","Year"),all.x=T)  
Data.monthly.GN.dusky=merge(Data.monthly.GN.dusky,Freo, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month","Year"),all.x=T)
Data.monthly.GN.sandbar=merge(Data.monthly.GN.sandbar,Freo, by.x=c("MONTH","YEAR.c"),
                              by.y=c("Month","Year"),all.x=T)

#daily
Data.daily.GN.whiskery=merge(Data.daily.GN.whiskery,SOI, by.x=c("MONTH","YEAR.c"),
                             by.y=c("Month.soi","Year.soi"),all.x=T)
Data.daily.GN.gummy=merge(Data.daily.GN.gummy,SOI, by.x=c("MONTH","YEAR.c"),
                          by.y=c("Month.soi","Year.soi"),all.x=T)  
Data.daily.GN.dusky=merge(Data.daily.GN.dusky,SOI, by.x=c("MONTH","YEAR.c"),
                          by.y=c("Month.soi","Year.soi"),all.x=T)
Data.daily.GN.sandbar=merge(Data.daily.GN.sandbar,SOI, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month.soi","Year.soi"),all.x=T)

Data.daily.GN.whiskery=merge(Data.daily.GN.whiskery,Freo, by.x=c("MONTH","YEAR.c"),
                             by.y=c("Month","Year"),all.x=T)
Data.daily.GN.gummy=merge(Data.daily.GN.gummy,Freo, by.x=c("MONTH","YEAR.c"),
                          by.y=c("Month","Year"),all.x=T)  
Data.daily.GN.dusky=merge(Data.daily.GN.dusky,Freo, by.x=c("MONTH","YEAR.c"),
                          by.y=c("Month","Year"),all.x=T)
Data.daily.GN.sandbar=merge(Data.daily.GN.sandbar,Freo, by.x=c("MONTH","YEAR.c"),
                            by.y=c("Month","Year"),all.x=T)


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

#GIS approach (complied by Dale Smith)
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




#4.6 Data fixes

#Remove initial years for sandbar as they were not reported in fishery stats
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,
         !(FINYEAR%in%c("1975-76","1976-77","1977-78","1978-79",
          "1979-80","1980-81","1981-82","1982-83","1983-84","1984-85")))

#subset sandbar data  to San.areas.range
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,LAT>=(-30) & LAT<=(-26))
Data.daily.GN.sandbar=subset(Data.daily.GN.sandbar,LAT>=(-30) & LAT<=(-26))  



#4.7 Remove NA effort

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



#4.8 Remove blocks with few sharks        
#note: this keeps blocks for which catch is within top 90% at least in one year

#sort blocks by catch (across entire history)
if(Remove.blk.by=="blk_only")
{
  Sort.blks.whi=fn.cum.ktch(Data.monthly.GN.whiskery,Data.daily.GN.whiskery,17003,Threshold=Prop.ktch.exp)
  Sort.blks.gum=fn.cum.ktch(Data.monthly.GN.gummy,Data.daily.GN.gummy,17001,Threshold=Prop.ktch.exp)
  Sort.blks.dus=fn.cum.ktch(Data.monthly.GN.dusky,Data.daily.GN.dusky,18003,Threshold=Prop.ktch.exp)
  Sort.blks.san=fn.cum.ktch(Data.monthly.GN.sandbar,Data.daily.GN.sandbar,18007,Threshold=Prop.ktch.exp)
  
  #these blocks were dropped out
  No.blok.whi=Sort.blks.whi$Dropped.blks
  No.blok.gum=Sort.blks.gum$Dropped.blks
  No.blok.dus=Sort.blks.dus$Dropped.blks
  No.blok.san=Sort.blks.san$Dropped.blks
  
  #these blocks were kept
  blok.whi=Sort.blks.whi$TopBloks
  blok.gum=Sort.blks.gum$TopBloks
  blok.dus=Sort.blks.dus$TopBloks
  blok.san=Sort.blks.san$TopBloks  
}

#number of records before removing blocks and vessels
  #Monthly
N.whi.all=nrow(Data.monthly.GN.whiskery)
N.gum.all=nrow(Data.monthly.GN.gummy )
N.dus.all=nrow(Data.monthly.GN.dusky)
N.san.all=nrow(Data.monthly.GN.sandbar)

  #Daily
N.whi.all_daily=nrow(Data.daily.GN.whiskery)
N.gum.all_daily=nrow(Data.daily.GN.gummy )
N.dus.all_daily=nrow(Data.daily.GN.dusky)
N.san.all_daily=nrow(Data.daily.GN.sandbar)

#keep copy of all blocks
#Monthtly
Data.monthly.GN.whiskery_all.blks=Data.monthly.GN.whiskery
Data.monthly.GN.gummy_all.blks=Data.monthly.GN.gummy
Data.monthly.GN.dusky_all.blks=Data.monthly.GN.dusky
Data.monthly.GN.sandbar_all.blks=Data.monthly.GN.sandbar

#Daily
Data.daily.GN.whiskery_all.blks=Data.daily.GN.whiskery
Data.daily.GN.gummy_all.blks=Data.daily.GN.gummy
Data.daily.GN.dusky_all.blks=Data.daily.GN.dusky
Data.daily.GN.sandbar_all.blks=Data.daily.GN.sandbar


#remove dropped out blocks
  #Monthtly
Data.monthly.GN.whiskery=subset(Data.monthly.GN.whiskery,BLOCKX%in%blok.whi)
Data.monthly.GN.gummy=subset(Data.monthly.GN.gummy,BLOCKX%in%blok.gum)
Data.monthly.GN.dusky=subset(Data.monthly.GN.dusky,BLOCKX%in%blok.dus)
Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,BLOCKX%in%blok.san)

  #Daily
Data.daily.GN.whiskery=subset(Data.daily.GN.whiskery,BLOCKX%in%blok.whi)
Data.daily.GN.gummy=subset(Data.daily.GN.gummy,BLOCKX%in%blok.gum)
Data.daily.GN.dusky=subset(Data.daily.GN.dusky,BLOCKX%in%blok.dus)
Data.daily.GN.sandbar=subset(Data.daily.GN.sandbar,BLOCKX%in%blok.san)

#export number of records lost from using selected blocks
LsT.rec.blk=data.frame(Species=c("Whi","Gum","Dus","San"),
  Monthly=c(100-100*nrow(Data.monthly.GN.whiskery)/N.whi.all,
            100-100*nrow(Data.monthly.GN.gummy)/N.gum.all,
            100-100*nrow(Data.monthly.GN.dusky)/N.dus.all,
            100-100*nrow(Data.monthly.GN.sandbar)/N.san.all),
  Daily=c(100-100*nrow(Data.daily.GN.whiskery)/N.whi.all_daily,
            100-100*nrow(Data.daily.GN.gummy)/N.gum.all_daily,
            100-100*nrow(Data.daily.GN.dusky)/N.dus.all_daily,
            100-100*nrow(Data.daily.GN.sandbar)/N.san.all_daily))
write.csv(LsT.rec.blk,"C:/Matias/Analyses/Catch and effort/Outputs/Blocks_kept/Percent.lost.csv")


#Create useful vars
FINYEAR.ALL=as.character(unique(c(Data.monthly.GN.whiskery$FINYEAR,Data.daily.GN.whiskery$FINYEAR)))
N.yrs.ALL=length(FINYEAR.ALL)

FINYEAR.monthly=as.character(unique(Data.monthly.GN.whiskery$FINYEAR))
N.yrs=length(FINYEAR.monthly)

FINYEAR.daily=as.character(unique(Data.daily.GN.whiskery$FINYEAR))
N.yrs.daily=length(FINYEAR.daily) 


#4.9 Remove spatial closures 
if(Remove.closure=="YES")
{
  #Monthly
    #Remove Metro closure
  Data.monthly.GN.whiskery=subset(Data.monthly.GN.whiskery,LAT> (-31) | LAT < (-33) )
  Data.monthly.GN.gummy=subset(Data.monthly.GN.gummy,LAT> (-31) | LAT < (-33) )
  Data.monthly.GN.dusky=subset(Data.monthly.GN.dusky,LAT> (-31) | LAT < (-33) )
  Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,LAT> (-31) | LAT < (-33) )
  
    #Remove whiskery closure
  Data.monthly.GN.whiskery=subset(Data.monthly.GN.whiskery,!(MONTH %in% 8:10 & LONG <118))
  Data.monthly.GN.gummy=subset(Data.monthly.GN.gummy,!(MONTH %in% 8:10 & LONG <118) )
  Data.monthly.GN.dusky=subset(Data.monthly.GN.dusky,!(MONTH %in% 8:10 & LONG <118) )
  Data.monthly.GN.sandbar=subset(Data.monthly.GN.sandbar,!(MONTH %in% 8:10 & LONG <118) )  
  
  #Daily
    #Remove Metro closure
  Data.daily.GN.whiskery=subset(Data.daily.GN.whiskery,LAT> (-31) | LAT < (-33) )
  Data.daily.GN.gummy=subset(Data.daily.GN.gummy,LAT> (-31) | LAT < (-33) )
  Data.daily.GN.dusky=subset(Data.daily.GN.dusky,LAT> (-31) | LAT < (-33) )
  Data.daily.GN.sandbar=subset(Data.daily.GN.sandbar,LAT> (-31) | LAT < (-33) )
  
    #Remove whiskery closure
  Data.daily.GN.whiskery=subset(Data.daily.GN.whiskery,!(MONTH %in% 8:10 & LONG <118))
  Data.daily.GN.gummy=subset(Data.daily.GN.gummy,!(MONTH %in% 8:10 & LONG <118) )
  Data.daily.GN.dusky=subset(Data.daily.GN.dusky,!(MONTH %in% 8:10 & LONG <118) )
  Data.daily.GN.sandbar=subset(Data.daily.GN.sandbar,!(MONTH %in% 8:10 & LONG <118) )  
}


#4.10 Data explorations

#Explore large catches                       
# Catch.quant.whis=aggregate(LIVEWT.c/1000~SPECIES,data=subset(Data.monthly.GN.whiskery,SPECIES==17003),quantile,probs=seq(0, 1, 0.1))
# Catch.quant.gumm=aggregate(LIVEWT.c/1000~SPECIES,data=subset(Data.monthly.GN.gummy,SPECIES==17001),quantile,probs=seq(0, 1, 0.1))
# Catch.quant.dusk=aggregate(LIVEWT.c/1000~SPECIES,data=subset(Data.monthly.GN.dusky,SPECIES==18003),quantile,probs=seq(0, 1, 0.1))
# Catch.quant.sand=aggregate(LIVEWT.c/1000~SPECIES,data=subset(Data.monthly.GN.sandbar,SPECIES==18007),quantile,probs=seq(0, 1, 0.1))


#Explore proportion of dusky and copper shark
    #monthly
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Data.monthly.GN.dusky,SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Data.monthly.GN.dusky,SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
par(las=1,mai=c(1.25,1.25,.1,.1))
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="copper catch / dusky catch",xlab="Finyear",pch=19,col=2,cex=1.75,cex.lab=1.5)
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)

    #daily
All.dusky=aggregate(LIVEWT.c~FINYEAR,subset(Data.daily.GN.dusky,SPECIES==18003),sum)
All.copper=aggregate(LIVEWT.c~FINYEAR,subset(Data.daily.GN.dusky,SPECIES==18001),sum)
All.dusky=All.dusky[match(All.copper$FINYEAR,All.dusky$FINYEAR),]
Prop.copper_dusky=data.frame(FINYEAR=All.dusky$FINYEAR,proportion=All.copper$LIVEWT.c/All.dusky$LIVEWT.c)
par(las=1,mai=c(1.25,1.25,.1,.1))
plot(1:nrow(Prop.copper_dusky),Prop.copper_dusky$proportion,xaxt='n',
     ylab="copper catch / dusky catch",xlab="Finyear",pch=19,col=2,cex=1.75,cex.lab=1.5)
axis(1,1:nrow(Prop.copper_dusky),Prop.copper_dusky$FINYEAR)



#4.11 Put data into a list
  #Selected blocks
Species.list=list(whiskery=Data.monthly.GN.whiskery,gummy=Data.monthly.GN.gummy,
                  dusky=Data.monthly.GN.dusky,sandbar=Data.monthly.GN.sandbar)

Species.list.daily=list(whiskery=Data.daily.GN.whiskery,gummy=Data.daily.GN.gummy,
                        dusky=Data.daily.GN.dusky,sandbar=Data.daily.GN.sandbar)

#All blocks
Species.list_all.blks=list(whiskery=Data.monthly.GN.whiskery_all.blks,gummy=Data.monthly.GN.gummy_all.blks,
                  dusky=Data.monthly.GN.dusky_all.blks,sandbar=Data.monthly.GN.sandbar_all.blks)

Species.list.daily_all.blks=list(whiskery=Data.daily.GN.whiskery_all.blks,gummy=Data.daily.GN.gummy_all.blks,
                        dusky=Data.daily.GN.dusky_all.blks,sandbar=Data.daily.GN.sandbar_all.blks)



#4.12 Illustrate folly effect: mean Vs sum
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


#4.13 Cluster analysis of daily data to identify targeting behaviour (Hoyle et al 2015)
if(do.cluster=="YES")
{
  Targt=c(17003,17001,18003,18007)
  for(i in 1:N.species)  fn.targeting(Species.list.daily[[i]],SP=Targt[i])  
}


#compare nominal all records VS 'good reporters' only
for(i in 1:N.species)
{
  tiff(file=paste("Outputs/Paper/All records vs Good records/",names(Species.list)[i],".tiff",sep=""),
       width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
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


#4.14 Construct wide database for analysis 
#notes: 
#   1. this constructs a single record for each record (i.e. year-month-vessel-block-gear for monthly
#      returns and year-date-vessel-block10-gear for daily logbooks, with catch of target
#      and other species as separate columns, giving a 0 catch for column "target" if no catch

#   2. only select "Good" records (the variable "Reporter" includes bad catch and bad effort)

#variables used in standardisation
  #monthly
These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Days.inv","Km.Gillnet.Days.c","Km.Gillnet.Days.c.no.creep",
                "SOI","zone","MONTH","BLOCKX","Freo","Freo.Lag6","Freo.Lag12")

  #daily
These.efforts.daily=c("FINYEAR","date","Km.Gillnet.Days.inv","Km.Gillnet.Days.c","Km.Gillnet.Days.c.no.creep",
                      "SOI","zone","MONTH","BLOCKX","block10","VESSEL","Freo","Freo.Lag6","Freo.Lag12")

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
DATA.list.LIVEWT.c.daily=vector('list',length=N.species)
Prop.Catch.daily=vector(length=N.species)
names(DATA.list.LIVEWT.c.daily)=names(Prop.Catch.daily)=names(Species.list.daily)
for ( i in 1:N.species)
{
  #create data sets 
  dummy=Effort.data.fun.daily(subset(Species.list.daily[[i]],Reporter=="good"),TARGETS[[i]],"LIVEWT.c")
  DATA.list.LIVEWT.c.daily[[i]]=dummy$dat
  
  #proportion with catch
  Prop.Catch.daily[i]=dummy$prop.with.catch
}

write.csv(Prop.Catch,paste(hndl,"Prop.records.with.catch.monthly.csv",sep=""),row.names=T)
write.csv(Prop.Catch.daily,paste(hndl,"Prop.records.with.catch.daily.csv",sep=""),row.names=T)


#4.15 Keep records within qualification levels 
QUALIF.LEVL=DATA.list.LIVEWT.c
QUALIF.LEVL.daily=DATA.list.LIVEWT.c.daily
for(i in 1:N.species)
{
  QUALIF.LEVL[[i]]=qual.levl.fun(DATA.list.LIVEWT.c[[i]],Prop.ktch.exp)
  QUALIF.LEVL.daily[[i]]=qual.levl.fun(DATA.list.LIVEWT.c.daily[[i]],Prop.ktch.exp)
}
  

tiff(file="C:/Matias/Analyses/Catch and effort/Outputs/Qualification levels.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(2,2),mar=c(1,3.6,1,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
for(i in 1:N.species) plot.qual.levl.fun(QUALIF.LEVL[[i]],Prop.ktch.exp)
mtext(paste("Proportion explaining",Prop.ktch.exp*100,"% of the catch"),2,outer=T,las=3,line=-1.5,cex=1.5)
mtext("Financial year",1,outer=T,line=1,cex=1.5)
dev.off()

#Keep records from indicative vessels only for monthly returns
Indi.ves=Indi.ves.sens=DATA.list.LIVEWT.c
DATA.list.LIVEWT.c.sens.ind.ves=DATA.list.LIVEWT.c
#Indi.ves.daily=Indi.ves.sens.daily=DATA.list.LIVEWT.c.daily
Store.DATA.list.LIVEWT.c=DATA.list.LIVEWT.c    #keep a copy of all data
#Store.DATA.list.LIVEWT.c.daily=DATA.list.LIVEWT.c.daily
#DATA.list.LIVEWT.c.daily.sens.ind.ves=DATA.list.LIVEWT.c.daily

for(i in 1:N.species)
{
  #monthly
  Indi.ves[[i]]=fn.indicative(Store.DATA.list.LIVEWT.c[[i]],Threshold.n.yrs,Min.rec.ves,second.criteria,Prop.ktch.exp*100)
  Indi.ves.sens[[i]]=fn.indicative(Store.DATA.list.LIVEWT.c[[i]],Threshold.n.yrs.sens,Min.rec.ves.sens,second.criteria,Prop.ktch.exp*100)
  
  if(second.criteria=="top.percent")
  {
    n=Indi.ves[[i]]$vess.store
    DATS=n    
    for(q in 1:length(n))DATS[[q]]=subset(Store.DATA.list.LIVEWT.c[[i]],VESSEL%in%n[[q]] & FINYEAR ==names(n)[q])
    DATA.list.LIVEWT.c[[i]]=do.call(rbind,DATS)
    
    if(!Criteria.indi=='all')
    {
      n=Indi.ves.sens[[i]]$vess.store
      DATS=n
      for(q in 1:length(n))DATS[[q]]=subset(Store.DATA.list.LIVEWT.c[[i]],VESSEL%in%n[[q]] & FINYEAR ==names(n)[q])
      DATA.list.LIVEWT.c.sens.ind.ves[[i]]=do.call(rbind,DATS)
    }
    
    if(Criteria.indi=='all') DATA.list.LIVEWT.c.sens.ind.ves[[i]]=Store.DATA.list.LIVEWT.c[[i]]
  }  
  
#   #daily
#   Indi.ves.daily[[i]]=fn.indicative(Store.DATA.list.LIVEWT.c.daily[[i]],Threshold.n.yrs,Min.rec.ves,second.criteria,Prop.ktch.exp*100)
#   Indi.ves.sens.daily[[i]]=fn.indicative(Store.DATA.list.LIVEWT.c.daily[[i]],Threshold.n.yrs,Min.rec.ves.sens,second.criteria,Prop.ktch.exp*100)
#   
#   if(second.criteria=="top.percent")
#   {
#     n=Indi.ves.daily[[i]]$vess.store
#     DATS=n    
#     for(q in 1:length(n))DATS[[q]]=subset(Store.DATA.list.LIVEWT.c.daily[[i]],VESSEL%in%n[[q]] & FINYEAR ==names(n)[q])
#     DATA.list.LIVEWT.c.daily[[i]]=do.call(rbind,DATS)
#     
#     if(!Criteria.indi=='all')
#     {
#       n=Indi.ves.sens.daily[[i]]$vess.store
#       DATS=n
#       for(q in 1:length(n))DATS[[q]]=subset(Store.DATA.list.LIVEWT.c.daily[[i]],VESSEL%in%n[[q]] & FINYEAR ==names(n)[q])
#       DATA.list.LIVEWT.c.daily.sens.ind.ves[[i]]=do.call(rbind,DATS)
#     }
#     
#     if(Criteria.indi=='all') DATA.list.LIVEWT.c.daily.sens.ind.ves[[i]]=Store.DATA.list.LIVEWT.c[[i]]
#   }
  
}


#Compare qualification levels, indicative vessels and all data

  #compare cpue from all records, qualification levels and indicative vessels
for(i in 1:N.species)
{
  dummy.all=Store.DATA.list.LIVEWT.c[[i]]
  dummy=DATA.list.LIVEWT.c[[i]]
  dummy.qual=QUALIF.LEVL[[i]]
  
  dummy.all$cpue=dummy.all$Catch.Target/dummy.all$Km.Gillnet.Days.c
  dummy$cpue=dummy$Catch.Target/dummy$Km.Gillnet.Days.c
  dummy.qual$cpue=dummy.qual$Catch.Target/dummy.qual$Km.Gillnet.Days.c
  
  all.cpue=aggregate(cpue~FINYEAR,dummy.all,mean)
  indi.cpue=aggregate(cpue~FINYEAR,dummy,mean)
  quali.cpue=aggregate(cpue~FINYEAR,dummy.qual,mean)
  
  par(mfcol=c(2,1),mai=c(.5,1.1,.3,.1),mgp=c(2.5,.5,0))
  plot(all.cpue$cpue,type='l',ylab='',main=names(Store.DATA.list.LIVEWT.c)[i],
       ylim=c(0,max(c(all.cpue$cpue,indi.cpue$cpue,quali.cpue$cpue))))
  lines(indi.cpue$cpue,col=2)
  lines(quali.cpue$cpue,col=3)
  legend('topright',c("all records","quali levls","indi. vess"),bty='n',lty=1,col=1:3)
  
  #normalised
  plot(all.cpue$cpue/mean(all.cpue$cpue),ylab='',type='l',ylim=c(0,1))
  lines(indi.cpue$cpue/mean(indi.cpue$cpue),col=2)
  lines(quali.cpue$cpue/mean(quali.cpue$cpue),col=3)
  mtext('nominal cpue',2,outer=T,line=-2,las=3,cex=2)
  
}


#4.16 Compute folly index for exporting         

#note: this has all records (good and bad reporters) from all blocks within effective area
San.Yrs=unique(DATA.list.LIVEWT.c[[4]]$FINYEAR)

Foly.whis=export.foly(subset(Data.monthly.GN.whiskery,SPECIES==17003),subset(Data.daily.GN.whiskery,SPECIES==17003))
Foly.gum=export.foly(subset(Data.monthly.GN.gummy,SPECIES==17001),subset(Data.daily.GN.gummy,SPECIES==17001))
Foly.dus=export.foly(subset(Data.monthly.GN.dusky,SPECIES%in%c(18003,18001)),subset(Data.daily.GN.dusky,SPECIES%in%c(18003,18001)))
Foly.san=export.foly(subset(Data.monthly.GN.sandbar,SPECIES==18007 & FINYEAR%in%San.Yrs),subset(Data.daily.GN.sandbar,SPECIES==18007))

#4.17 Compute nominal for exporting
#note: this has all records (good and bad reporters) from all blocks within effective area
Nominal.whis=export.nominal(subset(Data.monthly.GN.whiskery,SPECIES==17003),subset(Data.daily.GN.whiskery,SPECIES==17003))
Nominal.gum=export.nominal(subset(Data.monthly.GN.gummy,SPECIES==17001),subset(Data.daily.GN.gummy,SPECIES==17001))
Nominal.dus=export.nominal(subset(Data.monthly.GN.dusky,SPECIES%in%c(18003,18001)),subset(Data.daily.GN.dusky,SPECIES%in%c(18003,18001)))
Nominal.san=export.nominal(subset(Data.monthly.GN.sandbar,SPECIES==18007 & FINYEAR%in%San.Yrs),subset(Data.daily.GN.sandbar,SPECIES==18007))


List.foly.nom=list(whis=list(Foly=Foly.whis,Nom=Nominal.whis),
                   gum=list(Foly=Foly.gum,Nom=Nominal.gum),
                   dus=list(Foly=Foly.dus,Nom=Nominal.dus),
                   san=list(Foly=Foly.san,Nom=Nominal.san))


#4.18 Plot vessels and catch for monthly returns
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")

tiff(file="Appendix 7. Vessels and Catch by year.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfrow=c(4,2),mai=c(.45,.5,.1,.2),oma=c(.4,.1,.1,.95),mgp=c(1.5,.6,0),las=1)
for(i in 1:N.species) 
{
  fn.plot.vess.selection(Store.DATA.list.LIVEWT.c[[i]],DATA.list.LIVEWT.c[[i]])
  mtext(SPECIES.vec[i],4,las=3,line=0.75,cex=1.25)
}
mtext("Financial year",1,outer=T,line=-.75,cex=1.5)
dev.off()

Prop.vessels.dropped=Prop.catch.explained=Indi.ves
for(i in 1:N.species) 
{
  Prop.vessels.dropped[[i]]=Indi.ves[[i]]$Prop.Ves.dropped
  Prop.catch.explained[[i]]=Indi.ves[[i]]$Prop.ktch.exp.
}
write.csv(unlist(Prop.vessels.dropped),"Prop.of.vessels.dropped.csv")
write.csv(unlist(Prop.catch.explained),"Prop.of.catch.explained.after.dropping.vessels.csv")


if(check.interactions=="YES")
{
  for(i in 1:N.species) fn.vessl.cpue.yr(DATA.list.LIVEWT.c[[i]])
  for(i in 1:N.species) fn.vessl.cpue.yr(DATA.list.LIVEWT.c.daily[[i]])
}


#Show blocks, vessels and prop of records used in standardisation
hndl="C:/Matias/Analyses/Catch and effort/Outputs/Blocks_kept/"
Y=-36:-27; X=seq(113,129,length.out=length(Y))
blk.dropped=list(No.blok.whi,No.blok.gum,No.blok.dus,No.blok.san)
USED.VES=rep(NA,N.species)

for(i in 1:N.species)
{
  a=DATA.list.LIVEWT.c[[i]]
  blk.kept=unique(a$BLOCKX)
  tiff(file=paste(hndl,names(DATA.list.LIVEWT.c)[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.show.blk(blk.dropped[[i]],blk.kept,show.dropped="NO")
  dev.off()  
  USED.VES[i]=length(unique(a$VESSEL))
}
  

#Show number of records dropped in standardisation
PROP.recs.used=rep(NA,N.species)
PROP.recs.used[1]=nrow(Data.monthly.GN.whiskery)/N.whi.all
PROP.recs.used[2]=nrow(Data.monthly.GN.gummy )/N.gum.all
PROP.recs.used[3]=nrow(Data.monthly.GN.dusky)/N.dus.all
PROP.recs.used[4]=nrow(Data.monthly.GN.sandbar)/N.san.all
PROP.recs.dropped=round(1-PROP.recs.used,2)



#4.19 Input data tables

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
  FILE=paste("Mn_YR.cpue.",SPECIES.vec[i],".tiff",sep="")
  tiff(file=FILE,width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  
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
check.cpue.whis=check.cpue(DATA.list.LIVEWT.c[[1]],names(DATA.list.LIVEWT.c)[1],BYDAY="NO")
check.cpue.gum=check.cpue(DATA.list.LIVEWT.c[[2]],names(DATA.list.LIVEWT.c)[2],BYDAY="NO")
check.cpue.dus=check.cpue(DATA.list.LIVEWT.c[[3]],names(DATA.list.LIVEWT.c)[3],BYDAY="NO")
check.cpue.san=check.cpue(DATA.list.LIVEWT.c[[4]],names(DATA.list.LIVEWT.c)[4],BYDAY="NO")


#4.22 Plot data gaps

#Display proportion with and without catch for each species
MAX.list=c(45,100,30,15)
Ymax.vec=c(.24,.24,.24,.24)
Xmax.vec=c(11,11,11,11)

DATA.list.0.catch.2=vector('list',length=N.species)
names(DATA.list.0.catch.2)=SPECIES.vec


# Daily combined with monthly 
if(Separate.monthly.daily=="NO")
{
  for ( i in 1:N.species)DATA.list.0.catch.2[[i]]=Zero.catch.fun2(DATA.list.LIVEWT.c[[i]])
  
  tiff(file="Appendix 5.Prop of records with 0 catch.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
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
    plot(1:N.yrs,ZeroKTC,ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1,ylim=c(0,1)
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
  
  tiff(file="Appendix 5.Prop of records with 0 catch.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1.5,1,1,1.2),oma=c(2.5,3,.1,.1),las=1,mgp=c(.1,.5,0))
  for ( i in 1:N.species)
  {
    ZeroKTC=1-DATA.list.0.catch.2[[i]]$proportion
    LABS=as.character(DATA.list.0.catch.2[[i]]$FINYEAR)
    NN=length(LABS)
    plot(1:NN,ZeroKTC,ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1,ylim=c(0,1)
         ,cex.axis=1.5,yaxt="n",lwd=1)
    polygon(x=c(1,match("2005-06",LABS),match("2005-06",LABS),1),
            y=c(-0.5,-0.5,1.1,1.1),col="grey75",border="grey75")
    polygon(x=c(match("2005-06",LABS),NN,NN,match("2005-06",LABS)),
            y=c(-0.5,-0.5,1.1,1.1),col="grey95",border="grey95")
    points(1:NN,ZeroKTC,las=2,type="b",pch=19,cex=1.25,lwd=1.25)
    box()
    axis(1,at=1:NN,labels=F,tck=-0.015)
    axis(1,at=seq(1,NN,by=5),labels=LABS[seq(1,NN,by=5)],tck=-0.025,cex.axis=1.25)
    axis(2,at=seq(0,1,.2),labels=seq(0,1,.2),tck=-0.025,cex.axis=1.25)
    legend("topleft",SPECIES.vec[i],bty='n',cex=1.5)
    
  }
  mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.5,outer=T)
  mtext("Proportion of records with no catch",side=2,line=1.5,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
}


#Plot missing year-blocks for each species
Missing.yr.blk=DATA.list.LIVEWT.c
for (i in 1:N.species)
{
  tiff(file=paste("Figure 6",names(Species.list)[i],".Block.byYR.tiff"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
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




#4.23 Construct index of abundance

ZONES=c("West","Zone1","Zone2")
Predictors=c("FINYEAR","BLOCKX","MONTH","VESSEL","Catch.Gummy",
   "Catch.Whiskery","Catch.Dusky","Catch.Sandbar","Catch.Scalefish",
   "SOI","Freo","Freo.Lag6","Freo.Lag12","Km.Gillnet.Days.c") 

VARIABLES=c("Catch.Target",Predictors)

#Zone.list=list(c("West","Zone1","Zone2"),c("Zone1","Zone2"),c("West","Zone1","Zone2"),c("West","Zone1","Zone2"))

#remove catch column of target species
drop.pred=c(match("Catch.Whiskery",Predictors),match("Catch.Gummy",Predictors),
            match("Catch.Dusky",Predictors),match("Catch.Sandbar",Predictors))


#4.23.1 Explore data sets

#Number of observations per Month-Year-Block
    #Monthly
Whis.deg.free=Check.deg.free(DATA.list.LIVEWT.c[[1]])
Gum.deg.free=Check.deg.free(DATA.list.LIVEWT.c[[2]])
Dus.deg.free=Check.deg.free(DATA.list.LIVEWT.c[[3]])
San.deg.free=Check.deg.free(DATA.list.LIVEWT.c[[4]])

  #Daily
Whis.deg.free.daily=Check.deg.free(DATA.list.LIVEWT.c.daily[[1]])
Gum.deg.free.daily=Check.deg.free(DATA.list.LIVEWT.c.daily[[2]])
Dus.deg.free.daily=Check.deg.free(DATA.list.LIVEWT.c.daily[[3]])
San.deg.free.daily=Check.deg.free(DATA.list.LIVEWT.c.daily[[4]])

#Data summary
RecordThreshold=99 #number of records threshold (not used)
ZONAS=c("West","Zone1","Zone2")
Data.Summary=vector('list',length=N.species)

if(Model.run=="First")
{
  for ( i in 1:N.species)
  {
    tiff(file=paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Exploratory.monthly",SPECIES.vec[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")  
    par(mfrow=c(3,1),mar=c(4,4,1,1),las=1)
    Data.Summary[[i]]=fn.explore(DATA.list.LIVEWT.c[[i]],SPECIES.vec[i])
    mtext(SPECIES.vec[i], side = 3, line = -1, outer = T)
    dev.off()
  }
}


#4.23.2 Catch VS CPUE 
#note: IF catch increases and CPUE decreases then Effective Effort or Q increased!
if(Model.run=="First")
{
  for ( i in 1:N.species)
  {
    tiff(file=paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Catch_CPUE.monthly.",SPECIES.vec[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")  
    par(mfrow=c(1,1),mar=c(1,3.6,1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
    fn.catch.cpue(DATA.list.LIVEWT.c[[i]])
    dev.off()
    
    tiff(file=paste("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory/","Catch_CPUE.daily.",SPECIES.vec[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")  
    par(mfrow=c(1,1),mar=c(1,3.6,1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
    fn.catch.cpue(DATA.list.LIVEWT.c.daily[[i]])
    dev.off()
  }
}



#4.23.3 Number of vessels and blocks per zone and species
  #monthly
N.blk.N.ves.Zone=vector('list',N.species)
names(N.blk.N.ves.Zone)=SPECIES.vec
for ( i in 1:N.species)N.blk.N.ves.Zone[[i]]=fun.ves.blk.zone(DATA.list.LIVEWT.c[[i]])

  #daily
N.blk.N.ves.Zone.daily=N.blk.N.ves.Zone
for ( i in 1:N.species)N.blk.N.ves.Zone.daily[[i]]=fun.ves.blk.zone(DATA.list.LIVEWT.c.daily[[i]])


#4.23.4 Combine vessels in categories
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
  tiff(file="C:/Matias/Analyses/Catch and effort/Outputs/Paper/GLMM.catch.pred.vessel.cat_Monthly.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
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
  
}


#4.23.5 Explore predictor's effect

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
    PREDS=Predictors[-match("Km.Gillnet.Days.c",Predictors)]
    PREDS=PREDS[-drop.pred[i]]
    DATA$cpue=DATA$Catch.Target/DATA$Km.Gillnet.Days.c
    
    tiff(paste("Exp.",SPECIES.vec[i],".General.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
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
    Missing.Values[[i]]=Predictors_missing_values=colSums(is.na(DATA[,Preds]))/nrow(DATA)
    
    #Outliers response var
    boxplot(DATA$cpue~DATA$FINYEAR,main="Outliers in response var?",ylab="cpue (kg/km.gn.day)")
    dev.off()
    
    
    #boxplot of response var and predictors
    tiff(paste("Exp.",SPECIES.vec[i],".Box.plot.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
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
    tiff(paste("Exp.",SPECIES.vec[i],".Correlations.plot.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    Covars=DATA[,match(PREDS.cont,names(DATA))]
    Covars$MONTH=as.numeric(Covars$MONTH)
    pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor)
    dev.off()
  }
  
  Missing.Values  #No 0s?
  rm(M1.1)
  
}



#Compare raw vs reshaped data, cpue vs catch, effort as offset and main term, and vessel effect
if(COMPARE.RAW.etc=="YES")
{
  HNDL.cpre="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Compare cpue, catch, vessel, effort/"
  for(i in 1:N.species)
  {
    #Monthly
    tiff(file=paste(HNDL.cpre,SPECIES.vec[i],"_Monthly.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(3,1),mai=c(.25,.5,.1,.1),mgp=c(1.75,0.6,0),las=1)
    Compare.all(RAW=subset(Species.list[[i]],SPECIES%in%TARGETS[[i]] & LIVEWT.c>0),
                RESHAPED=DATA.list.LIVEWT.c[[i]],
                RESHAPED.ALL.VES= Store.DATA.list.LIVEWT.c[[i]])
    dev.off()
    
    #Daily
    tiff(file=paste(HNDL.cpre,SPECIES.vec[i],"_Daily.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(3,1),mai=c(.25,.5,.1,.1),mgp=c(1.75,0.6,0),las=1)
    Compare.all(RAW=subset(Species.list.daily[[i]],SPECIES==TARGETS[[i]] & LIVEWT.c>0),
                RESHAPED=DATA.list.LIVEWT.c.daily[[i]],
                RESHAPED.ALL.VES= DATA.list.LIVEWT.c.daily[[i]])
    dev.off()
    
    
    #random effects
      #Monthly
    tiff(file=paste(HNDL.cpre,SPECIES.vec[i],"_Monthly_main.VS.random.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(1,1),mai=c(.5,.75,.1,.1),mgp=c(1.75,0.6,0),las=1)
    Compare.vess.random(RESHAPED=DATA.list.LIVEWT.c[[i]])
    dev.off()
    
      #Daily
    tiff(file=paste(HNDL.cpre,SPECIES.vec[i],"_Daily_main.VS.random.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(1,1),mai=c(.5,.75,.1,.1),mgp=c(1.75,0.6,0),las=1)
    Compare.vess.random(RESHAPED=DATA.list.LIVEWT.c.daily[[i]])
    dev.off()    
  }
}


#Deje aca

#4.23.6 Define best model and error structure 
if(Combine.ves=="YES") VARIABLES=c(VARIABLES,"new.level.vessel")

#With Yr-Blk interactions           #this takes 29 mins
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
  for(j in 1:N.species) Tab1.Day[[j]]=fn.Tabl1(DATA.list.LIVEWT.c.daily[[j]])
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
    Species.model=fn.define.models(Fit.to.what[w])    
    All.models=vector('list',length(SPECIES.vec))
    names(All.models)=SPECIES.vec
    All.models.daily=All.models
    DAT=names(DATA.list.LIVEWT.c)  
    LEVELS.ves.blk=All.models
    
    #Monthly
    DIR=paste(Mon.wd,"/",Fit.to.what[w],sep='')
    setwd(DIR)
    for(j in 1:N.species)     
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
    for(i in 1:N.species)
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
      for(j in 1:N.species)     
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
      for(i in 1:N.species)
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


#Write out best model
Best.Model=vector('list',N.species)
names(Best.Model)=SPECIES.vec

  #Monthly
Best.Model$'Whiskery shark'=list(
  Bi=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL+MONTH+log.Catch.Sandbar+offset(log.Effort)), 
  cpue=as.formula(log.cpue~FINYEAR*BLOCKX+VESSEL+MONTH),
  catch=as.formula(log.Catch~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)),
  catch.Gamma=as.formula(Catch.Target~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)))

Best.Model$'Gummy shark'=list(
  Bi=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL+MONTH+log.Catch.Sandbar+offset(log.Effort)), 
  cpue=as.formula(log.cpue~FINYEAR*BLOCKX+VESSEL+MONTH),
  catch=as.formula(log.Catch~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)),
  catch.Gamma=as.formula(Catch.Target~FINYEAR*BLOCKX+VESSEL+offset(log.Effort)))

Best.Model$'Dusky shark'=list(
  Bi=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL+MONTH+log.Catch.Gummy+log.Catch.Whiskery+log.Catch.Sandbar+offset(log.Effort)), 
  cpue=as.formula(log.cpue~FINYEAR*BLOCKX+VESSEL+MONTH+log.Catch.Whiskery),
  catch=as.formula(log.Catch~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)),
  catch.Gamma=as.formula(Catch.Target~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)))

Best.Model$'Sandbar shark'=list(
  Bi=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL+MONTH+log.Catch.Gummy+log.Catch.Whiskery+offset(log.Effort)), 
  cpue=as.formula(log.cpue ~ FINYEAR * BLOCKX + VESSEL + MONTH),
  catch=as.formula(log.Catch~FINYEAR*BLOCKX+VESSEL+MONTH+log.Catch.Dusky+offset(log.Effort)),
  catch.Gamma=as.formula(Catch.Target~FINYEAR*BLOCKX+VESSEL+MONTH+offset(log.Effort)))

  #Daily
Best.Model.daily=Best.Model
Best.Model.daily$'Whiskery shark'$Bi=as.formula(Catch.Target~FINYEAR + BLOCKX + VESSEL+MONTH+log.Catch.Scalefish+offset(log.Effort))
Best.Model.daily$'Whiskery shark'$catch=as.formula(log.Catch ~ FINYEAR * BLOCKX + VESSEL + log.Catch.Scalefish + offset(log.Effort))
  
Best.Model.daily$'Gummy shark'$Bi=as.formula(Catch.Target ~ FINYEAR + BLOCKX + VESSEL+MONTH+log.Catch.Whiskery+offset(log.Effort))
Best.Model.daily$'Gummy shark'$catch=as.formula(log.Catch ~ FINYEAR * BLOCKX + VESSEL + offset(log.Effort))
  
Best.Model.daily$'Dusky shark'$Bi=as.formula(Catch.Target ~ FINYEAR + BLOCKX + VESSEL + MONTH + offset(log.Effort))
  
Best.Model.daily$'Sandbar shark'$Bi=as.formula(Catch.Target ~ FINYEAR + BLOCKX + VESSEL + MONTH + offset(log.Effort))
Best.Model.daily$'Sandbar shark'$catch=as.formula(log.Catch ~ FINYEAR * BLOCKX + VESSEL + MONTH + offset(log.Effort))


#Compare best gamma and best lognormal
if(Compare.best.lognormal.gamma=="YES")
{
  for(i in 1:N.species)
  {
    LogN.form=Best.Model[[i]]$catch
    Gama.form=Best.Model[[i]]$catch.Gamma
    HnDl="C:/Matias/Analyses/Catch and effort/Outputs/ErrorStructure/"
    
    #Monthly
    aa=fn.LogN.vs.Gama(DATA.list.LIVEWT.c[[i]])
    LogN=aa$LOgN
    GaMa=aa$Gama
    DeV.EXP=data.frame(cbind(fn.Dev.Exp(LogN),fn.Dev.Exp(GaMa)))
    

    tiff(file=paste(HnDl,"Monthly.",SPECIES.vec[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mfcol=c(3,2))
    Pos.Diag.fn(LogN,SPECIES.vec[i])
    legend("topright","LOGN",bty='n',cex=1.5)
    Pos.Diag.fn(GaMa,SPECIES.vec[i])
    legend("topright","Gamma",bty='n',cex=1.5)
    dev.off()
    
      
    #Daily
    aa=fn.LogN.vs.Gama(DATA.list.LIVEWT.c.daily[[i]])
    LogN=aa$LOgN
    GaMa=aa$Gama
    
    tiff(file=paste(HnDl,"Daily.",SPECIES.vec[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    
    par(mfcol=c(3,2))
    Pos.Diag.fn(LogN,SPECIES.vec[i])
    legend("topright","LOGN",bty='n',cex=1.5)
    Pos.Diag.fn(GaMa,SPECIES.vec[i])
    legend("topright","Gamma",bty='n',cex=1.5)
    dev.off()
    
    DeV.EXP=rbind(DeV.EXP,data.frame(cbind(fn.Dev.Exp(LogN),fn.Dev.Exp(GaMa))))
    rownames(DeV.EXP)=c("Monthly","Daily")
    colnames(DeV.EXP)=c("LogN","Gamma")
    write.csv(DeV.EXP,paste(HnDl,"Dev.Exp.",SPECIES.vec[i],".csv",sep=""))  
    
    rm(aa,LogN,GaMa,DeV.EXP)
  }  
}

#Compare best cpue and best catch models
if(Compare.best.catch.cpue=="YES")
{
  stores.cpue=DATA.list.LIVEWT.c
  stores.catch=stores.cpue
  for(i in 1:N.species) 
  {
    stores.cpue[[i]]=fn.compr.cpue.ktch(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Best.Model[[i]]$cpue,VARIABLES,Min.weight[i],Int=.1)
    stores.catch[[i]]=fn.compr.cpue.ktch(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Best.Model[[i]]$catch,VARIABLES,Min.weight[i],Int=.1)
  }
  par(mfcol=c(2,1),mai=c(.2,.2,.2,.2))
  for(i in 1:N.species) 
  {
    plot(stores.cpue[[i]]$log.cpue,stores.cpue[[i]]$Pred)
    lines(stores.cpue[[i]]$log.cpue,stores.cpue[[i]]$log.cpue,col=2)
    legend("topleft","cpue",bty='n')
    plot(stores.catch[[i]]$log.Catch,stores.catch[[i]]$Pred)
    lines(stores.catch[[i]]$log.Catch,stores.catch[[i]]$log.Catch,col=2)
    legend("topleft","catch",bty='n')
  }
  
}

#K-fold cross validation approach
if(Do.K_n.fold.test=="YES")
{
  K=5
  training.data=0.9 #training subset (90%)
  
  
  Best.error.catch.k_fold=vector('list',length=N.species)
  names(Best.error.catch.k_fold)=SPECIES.vec
  
  system.time(for ( i in 1:N.species) 
  {
    Best.error.catch.k_fold[[i]]=fn.best.error.k_fold(DATA.list.LIVEWT.c[[i]],Best.pos[[i]],Best.pos.log[[i]])
  })
  
  #extract n-fold values
  MEAN=SD=NULL
  for(k in 1:N.species)
  {
    for (i in 1:4)
    {
      MEAN=c(MEAN,mean(Best.error.catch.k_fold[[j]][[i]]))
      SD=c(SD,sd(Best.error.catch.k_fold[[j]][[i]]))
    }
    DAT=data.frame(x=1:4,mean=MEAN,lower=MEAN-SD,upper=MEAN+SD)
    with(DAT,errbar(x,mean,upper,lower,cex=1.5,lwd=1.5,xaxt='n'))
    grid(nx=NA,ny=NULL)
  }
  
  
  for(i in 1:N.species) write.table(do.call(rbind,Best.error.catch.k_fold[[i]]),paste("N-fold.",SPECIES.vec[i],".csv",sep=""),sep = ",")
  
  
  #Free up memory
  rm(list=c("Best.error.catch"))
  
}


#4.23.7 Best model 
  #Monthly
Best.Model$'Whiskery shark'=list(Bi=Best.Model$'Whiskery shark'$Bi, 
                                 Log=Best.Model$'Whiskery shark'$catch)

Best.Model$'Gummy shark'=list(Bi=Best.Model$'Gummy shark'$Bi, 
                              Log=Best.Model$'Gummy shark'$catch)

Best.Model$'Dusky shark'=list(Bi=Best.Model$'Dusky shark'$Bi, 
                              Log=Best.Model$'Dusky shark'$catch)

Best.Model$'Sandbar shark'=list(Bi=Best.Model$'Sandbar shark'$Bi, 
                                Log=Best.Model$'Sandbar shark'$catch)


  #Daily
Best.Model.daily$'Whiskery shark'=list(Bi=Best.Model.daily$'Whiskery shark'$Bi, 
                                      Log=Best.Model.daily$'Whiskery shark'$catch)

Best.Model.daily$'Gummy shark'=list(Bi=Best.Model.daily$'Gummy shark'$Bi, 
                                    Log=Best.Model.daily$'Gummy shark'$catch)

Best.Model.daily$'Dusky shark'=list(Bi=Best.Model.daily$'Dusky shark'$Bi, 
                                    Log=Best.Model.daily$'Dusky shark'$catch)

Best.Model.daily$'Sandbar shark'=list(Bi=Best.Model.daily$'Sandbar shark'$Bi, 
                                    Log=Best.Model.daily$'Sandbar shark'$catch)


#no interactions
if(With.interact=="NO")
{
  Best.Model.no.int=Best.Model
  Best.Model.no.int$'Whiskery shark'$Log=log.Catch ~ FINYEAR + BLOCKX + VESSEL + MONTH + offset(log.Effort)
  Best.Model.no.int$'Gummy shark'$Log=log.Catch ~ FINYEAR + BLOCKX + VESSEL + MONTH + offset(log.Effort)
  Best.Model.no.int$'Dusky shark'$Log=log.Catch ~ FINYEAR + BLOCKX + VESSEL + MONTH + log.Catch.Whiskery + offset(log.Effort)
  Best.Model.no.int$'Sandbar shark'$Log=log.Catch ~ FINYEAR + BLOCKX + VESSEL + MONTH + offset(log.Effort)  
}
      

#model vars
VARS=c(VARIABLES,"Km.Gillnet.Days.c.no.creep")

#Explore predictions for improving fit
if(Def.mod.Str.=="YES" & improve=="YES")
{
  stores=DATA.list.LIVEWT.c
  for(i in 1:N.species) 
  {
    stores[[i]]=fit.model.to.inspect(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Best.Model[[i]]$Log,VARS,Min.weight[i],Int=.1)
  }
  
  #manual bit:
  #use locator to select points departing from one to one and evaluate distribution of Months, vessels, blocks,
  i=2
  Agg=aggregate(cbind(Pred,log.cpue)~BINS,stores[[i]],mean)
  locator(2)
  Check=subset(Agg,log.cpue<(-1.46) | log.cpue>3.6)
  Unik.bins=unique(Check$BINS)
  dummy=subset(stores[[i]],BINS%in%Unik.bins)
  table(dummy$FINYEAR,dummy$MONTH)
  table(dummy$FINYEAR,dummy$BLOCKX)
  table(dummy$VESSEL)
  dummy=fit.model.to.inspect(subset(stores[[i]],!BINS%in%Unik.bins),Ves.no.ktc[[i]],Best.Model[[i]]$Log,VARS,Min.weight[i],Int=.1)
  rm(stores,dummy) 
}


#Free up memory
Free.up="NO"
if(Free.up=="YES")
{
  rm(list=c("Species.list","Data.monthly.GN.dusky","Data.monthly.GN.whiskery",
            "Data.monthly.GN.sandbar","Data.monthly.GN.gummy"))
  rm(list=c("Species.list.daily","Data.daily.GN.dusky","Data.daily.GN.whiskery",
            "Data.daily.GN.sandbar","Data.daily.GN.gummy"))
}


#4.23.8 Standardise catches

#Apply best model to all scenarios
Stand.BaseCase=vector('list',length=N.species)
names(Stand.BaseCase)=SPECIES.vec
Stand.Fantasy=Stand.GoodRep=Stand.no.creep=Stand.ind.vess=
Stand.BaseCase.daily=Stand.no.interaction=Stand.BaseCase

system.time(for ( i in 1:N.species)     #note: this takes ~ 3 mins
{
  #Monthly
    #base case
  Stand.BaseCase[[i]]=fn.stand.cpue(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],SPECIES.vec[[i]],
                                    Best.Model[[i]]$Bi,Best.Model[[i]]$Log,VARS,Eff.Creep="YES")
    #sensitivities
  if(do.sensitivity=="YES")
  {
    #No effort creep
    Stand.no.creep[[i]]=fn.stand.cpue(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],SPECIES.vec[[i]],
                  Best.Model[[i]]$Bi,Best.Model[[i]]$Log,VARS,Eff.Creep="NO")
    
     #Indicative vessel sensitivity
     Stand.ind.vess[[i]]=fn.stand.cpue(DATA.list.LIVEWT.c.sens.ind.ves[[i]],Ves.no.ktc[[i]],SPECIES.vec[[i]],
                 Best.Model[[i]]$Bi,Best.Model[[i]]$Log,VARS,Eff.Creep="YES")
    
    #No YR-BLOCK interaction
#     Stand.no.interaction[[i]]=fn.stand.cpue(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],SPECIES.vec[[i]],
#                                 Best.Model.no.int[[i]]$Bi,Best.Model.no.int[[i]]$Log,VARS,Eff.Creep="YES")
  }
  

    #Daily
  if(Separate.monthly.daily=="YES")
  {
    #base case
    Stand.BaseCase.daily[[i]]=fn.stand.cpue(DATA.list.LIVEWT.c.daily[[i]],Ves.no.ktc[[i]],SPECIES.vec[[i]],
                                     Best.Model.daily[[i]]$Bi,Best.Model.daily[[i]]$Log,VARS,Eff.Creep="YES")
  
}
  
})


#Plot model diagnostics and extract deviance explained
setwd("C:/Matias/Analyses/Catch and effort/Outputs/model.fits")

#base case observed VS predicted
system.time(for ( i in 1:N.species)     
{
    #Monthly
  tiff(paste(SPECIES.vec[i],".preds.vs.obs.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  see.BC.fit.to.data(GLM=Stand.BaseCase[[i]],RespVar=Fit.to.what)
  dev.off()
  
    #Daily
  if(Separate.monthly.daily=="YES")
  {
    tiff(paste(SPECIES.vec[i],".preds.vs.obs.daily.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    see.BC.fit.to.data(GLM=Stand.BaseCase.daily[[i]],RespVar=Fit.to.what)
     dev.off()
  }
  
})

#Deje aca
#output AIC and deviance explained
Dev.Expl=data.frame(Species=SPECIES.vec,Binomial=NA,LogNormal=NA)
AIC.species=Dev.Expl.daily=AIC.species.daily=Dev.Expl

for ( i in 1:N.species)
{
    #monthly
    tiff(paste(SPECIES.vec[i],"Fit.Diag.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    fn.plot.diag(Stand.BaseCase[[i]]$LogN,"LogNormal",SPECIES.vec[i])
    dev.off()
    
    #daily
    if(Separate.monthly.daily=="YES")
    {
      tiff(paste(SPECIES.vec[i],"Fit.Diag.daily.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      fn.plot.diag(Stand.BaseCase.daily[[i]]$LogN,"LogNormal",SPECIES.vec[i])
      dev.off()
      
    }
    
    #print(Stand.BaseCase[[i]]$LogN, corr = FALSE) #print without the correlation matrix
     
    #monthly
    Dev.Expl$LogNormal[i]=round(fn.Dev.Exp(Stand.BaseCase[[i]]$LogN)*100)
    AIC.species$LogNormal[i]=fn.AICc(Stand.BaseCase[[i]]$LogN)
    if(!is.null(Stand.BaseCase[[i]]$Bi))
      {
        Dev.Expl$Binomial[i]=round(fn.Dev.Exp(Stand.BaseCase[[i]]$Bi)*100)
        AIC.species$Binomial[i]=fn.AICc(Stand.BaseCase[[i]]$Bi)
      }
    
    #daily
    if(Separate.monthly.daily=="YES")
    {
      Dev.Expl.daily$LogNormal[i]=round(fn.Dev.Exp(Stand.BaseCase.daily[[i]]$LogN)*100)
      AIC.species.daily$LogNormal[i]=fn.AICc(Stand.BaseCase.daily[[i]]$LogN)
      if(!is.null(Stand.BaseCase.daily[[i]]$Bi))
      {
        Dev.Expl.daily$Binomial[i]=round(fn.Dev.Exp(Stand.BaseCase.daily[[i]]$Bi)*100)
        AIC.species.daily$Binomial[i]=fn.AICc(Stand.BaseCase.daily[[i]]$Bi)
      }
    }

  }
write.csv(Dev.Expl,"Dev.expalined.csv",row.names=F)
write.csv(AIC.species,"AICc.csv",row.names=F)

if(Separate.monthly.daily=="YES")
{
  write.csv(Dev.Expl.daily,"Dev.expalined.daily.csv",row.names=F)
  write.csv(AIC.species.daily,"AICc.daily.csv",row.names=F)
}

#Extract SE and coefs
if(Extract.SE=="YES")
{
  Store.Coefs.BaseCase=vector('list',length=N.species)
  names(Store.Coefs.BaseCase)=SPECIES.vec
  for ( i in 1:N.species)Store.Coefs.BaseCase[[i]]=fn.get.stuff(Stand.BaseCase[[i]],SPECIES.vec[i])
}

#Deviance Tables
if(Extract.Deviance.table=="YES")
{
  DevianceTable=vector('list',N.species)
  for ( i in 1:N.species)DevianceTable[[i]]=fn.get.Significance(Stand.BaseCase[[i]],SPECIES.vec[i])
}


#Figure diagnostics
#fit diagnostics
setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")

#Binomial
if(do.binomial.fit=="YES")
{
  #monthly
  tiff("Figure Binomial.Fit.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  #par(mfcol=c(3,1),las=1,mar=c(2,6,.1,6),oma=c(2,.1,.1,4),las=1,mgp=c(2,.5,0),cex.axis=1.5,cex.lab=1.1)
  par(mfcol=c(2,2),mar=c(1,3,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
  for ( i in 1:N.species)
  {
    Bi.Diag.fn(DATA.list.LIVEWT.c[[i]],Ves.no.ktc[[i]],Stand.BaseCase[[i]]$Bi,0.1)
    legend("bottomright",SPECIES.vec[i],bty='n',cex=1.75)
  }
  mtext("Predicted fraction of non-zero catches",1,outer=T,line=0.5,cex=1.35)
  mtext("Observed fraction of non-zero catches",2,outer=T,line=-2.5,cex=1.35,las=3)
  axis(1,seq(0,1,0.1),seq(0,1,0.1),cex=1.5)
  dev.off()
  
  
  #daily
  tiff("Figure Binomial.Fit.daily.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,2),mar=c(1,3,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
  for ( i in 1:N.species)
  {
    Bi.Diag.fn(DATA.list.LIVEWT.c.daily[[i]],Ves.no.ktc[[i]],Stand.BaseCase.daily[[i]]$Bi,0.1)
    legend("bottomright",SPECIES.vec[i],bty='n',cex=1.75)
  }
  mtext("Predicted fraction of non-zero catches",1,outer=T,line=1.25,cex=1.35)
  mtext("Observed fraction of non-zero catches",2,outer=T,line=-2.5,cex=1.35,las=3)
  axis(1,seq(0,1,0.1),seq(0,1,0.1),cex=1.5)
  dev.off()
  
}

#Positive catch
tiff("Appendix 8.monthly.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(3,4),las=1,mar=c(3,4,1.75,0.1),oma=c(1,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
for ( i in 1:N.species) Pos.Diag.fn(Stand.BaseCase[[i]]$LogN,SPECIES.vec[i])
dev.off()

if(Separate.monthly.daily=="YES")
{
  tiff("Appendix 8.daily.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,4),las=1,mar=c(3,4,1.75,0.1),oma=c(1,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
  for ( i in 1:N.species) Pos.Diag.fn(Stand.BaseCase.daily[[i]]$LogN,SPECIES.vec[i])
  dev.off()
}


#check the contrasts options used by glm
getOption("contrasts")


#Show month coefficient
tiff("Month coefs_Monthly.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,2),mar=c(1,3,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
for(i in 1:N.species) fn.pred.month(Stand.BaseCase[[i]],SPECIES.vec[i])
dev.off()
tiff("Month coefs_Daily.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,2),mar=c(1,3,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
for(i in 1:N.species) fn.pred.month(Stand.BaseCase.daily[[i]],SPECIES.vec[i])
dev.off()

#4.23.9 Predict standardised catch

#create Base case data 
FINYEAR=fn.extract(sort(unique(DATA.list.LIVEWT.c[[1]]$FINYEAR)))
FINYEAR.san=fn.extract(sort(unique(DATA.list.LIVEWT.c[[4]]$FINYEAR)))

BLOCKX=BLOCKX.good=BLOCKX.daily=vector('list',length=N.species)
for ( i in 1:N.species) BLOCKX[[i]]=unique(Stand.BaseCase[[i]]$BiData$BLOCKX)
#Select levels of factors
new.MN=new.level.vessel=Get.cof=vector('list',length=N.species)

#Approach 1. Most common value of categorical variables
if(Sel.lev=='most.common')
{
  #month
  MONTH=fn.extract(sort(unique(DATA.list.LIVEWT.c[[1]]$MONTH)))
  for ( i in 1:N.species)
  {
    new.MN[[i]]=Most.common(DATA.list.LIVEWT.c[[i]],"MONTH")
    a=MONTH[as.numeric(new.MN[[i]])]
    new.MN[[i]]=list(bi=a,pos=a)
  }
  
  #vessels
  for ( i in 1:N.species)
  {
    new.level.vessel[[i]]=Most.common(DATA.list.LIVEWT.c[[i]],"VESSEL")
    Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c[[i]]$VESSEL)))
    a=Unik.ves[match(new.level.vessel[[i]],Unik.ves)]
    new.level.vessel[[i]]=list(bi=a,pos=a)
  }
}

#Approach 2. Weighted average of month and vessel terms
if(Sel.lev=='weighted.aver')
{
  for ( i in 1:N.species)
  {
    a=fn.get.coef(Stand.BaseCase[[i]])
    if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.bi,pos=a$ID.mnth.pos)
    if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.pos,pos=a$ID.mnth.pos)
    Unik.mn=fn.extract(sort(unique(DATA.list.LIVEWT.c[[i]]$MONTH)))
    new.MN[[i]]=list(bi=Unik.mn[match(dummy$bi,Unik.mn)],pos=Unik.mn[match(dummy$pos,Unik.mn)])
    
    if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.bi,pos=a$ID.ves.pos)
    if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.pos,pos=a$ID.ves.pos)
    Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c[[i]]$VESSEL)))
    new.level.vessel[[i]]=list(bi=Unik.ves[match(dummy$bi,Unik.ves)],pos=Unik.ves[match(dummy$pos,Unik.ves)])
    
  }
}

#Select mean covariates
FREO=soi=KTch.Gum=KTch.Whi=KTch.Dus=KTch.San=KTch.Fish=log.Effort=vector(length=N.species)
for ( i in 1:N.species)
{
  FREO[i]=Mn.cov(DATA.list.LIVEWT.c[[i]],"Freo")
  soi[i]=Mn.cov(DATA.list.LIVEWT.c[[i]],"SOI")
  KTch.Gum[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Catch.Gummy"))
  KTch.Whi[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Catch.Whiskery"))
  KTch.Dus[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Catch.Dusky"))
  KTch.San[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Catch.Sandbar"))
  KTch.Fish[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Catch.Scalefish"))
  log.Effort[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Km.Gillnet.Days.c"))
}
 
#create Indicative vessel sensitivity
if(do.sensitivity=="YES")
{
  BLOCKX.sens.ind.ves=BLOCKX
  for ( i in 1:N.species) BLOCKX.sens.ind.ves[[i]]=fn.extract(sort(unique(DATA.list.LIVEWT.c.sens.ind.ves[[i]]$BLOCKX)))
  
  new.MN.sens.ind.ves=new.level.vessel.sens.ind.ves=new.MN
  
  if(Sel.lev=='most.common')
  {
    #get Most common value of categorical variables
    #month
    for ( i in 1:N.species)
    {
      new.MN.sens.ind.ves[[i]]=Most.common(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"MONTH")
      a=MONTH[as.numeric(new.MN.sens.ind.ves[[i]])]
      new.MN.sens.ind.ves[[i]]=list(bi=a,pos=a)
    }
    
    #vessels
    for ( i in 1:N.species)
    {
      new.level.vessel.sens.ind.ves[[i]]=Most.common(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"VESSEL")
      Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c.sens.ind.ves[[i]]$VESSEL)))
      a=Unik.ves[match(new.level.vessel.sens.ind.ves[[i]],Unik.ves)]   
      new.level.vessel.sens.ind.ves[[i]]=list(bi=a,pos=a)
    }
  }
  
  if(Sel.lev=='weighted.aver')
  {
    #Approach 2. Weighted average of month and vessel terms
    for ( i in 1:N.species)
    {
      if(!All.ves.dat.same.BC.dat=="YES") a=fn.get.coef(Stand.ind.vess[[i]])
      if(All.ves.dat.same.BC.dat=="YES") a=fn.get.coef(Stand.BaseCase[[i]])
      
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.bi,pos=a$ID.mnth.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.pos,pos=a$ID.mnth.pos)
      Unik.mn=fn.extract(sort(unique(DATA.list.LIVEWT.c.sens.ind.ves[[i]]$MONTH)))
      new.MN.sens.ind.ves[[i]]=list(bi=Unik.mn[match(dummy$bi,Unik.mn)],pos=Unik.mn[match(dummy$pos,Unik.mn)])
      
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.bi,pos=a$ID.ves.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.pos,pos=a$ID.ves.pos)
      Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c.sens.ind.ves[[i]]$VESSEL)))
      new.level.vessel.sens.ind.ves[[i]]=list(bi=Unik.ves[match(dummy$bi,Unik.ves)],pos=Unik.ves[match(dummy$pos,Unik.ves)])
    }
  }
  
  
  #mean covariates
  FREO.sens.ind.ves=soi.sens.ind.ves=KTch.Gum.sens.ind.ves=KTch.Whi.sens.ind.ves=
    KTch.Dus.sens.ind.ves=KTch.San.sens.ind.ves=KTch.Fish.sens.ind.ves=log.Effort.sens.ind.ves=FREO
  for ( i in 1:N.species)
  {
    if(All.ves.dat.same.BC.dat=="YES")
    {
      FREO.sens.ind.ves[i]=FREO[i]
      soi.sens.ind.ves[i]=soi[i]
      KTch.Gum.sens.ind.ves[i]=KTch.Gum[i]
      KTch.Whi.sens.ind.ves[i]=KTch.Whi[i]
      KTch.Dus.sens.ind.ves[i]=KTch.Dus[i]
      KTch.San.sens.ind.ves[i]=KTch.San[i]
      KTch.Fish.sens.ind.ves[i]=KTch.Fish[i]
      log.Effort.sens.ind.ves[i]=log.Effort[i]
    }
    if(!All.ves.dat.same.BC.dat=="YES")
    {
      FREO.sens.ind.ves[i]=Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Freo")
      soi.sens.ind.ves[i]=Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"SOI")
      KTch.Gum.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Catch.Gummy"))
      KTch.Whi.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Catch.Whiskery"))
      KTch.Dus.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Catch.Dusky"))
      KTch.San.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Catch.Sandbar"))
      KTch.Fish.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Catch.Scalefish"))
      log.Effort.sens.ind.ves[i]=log(Mn.cov(DATA.list.LIVEWT.c.sens.ind.ves[[i]],"Km.Gillnet.Days.c"))
      
    }
  }
  
}

#create no effort creep sensitivity
if(do.sensitivity=="YES")
{
  new.MN.sens.creep=new.level.vessel.sens.creep=new.MN
  if(Sel.lev=='most.common')
  {
    #get Most common value of categorical variables
    #month
    for ( i in 1:N.species) new.MN.sens.creep[[i]]=new.MN[[i]]
    
    #vessels
    for ( i in 1:N.species)new.level.vessel.sens.creep[[i]]=new.level.vessel[[i]]
  }
  
  if(Sel.lev=='weighted.aver')
  {
    #Approach 2. Weighted average of month and vessel terms
    for ( i in 1:N.species)
    {
      if(!Creep.dat.same.BC.dat=="YES") a=fn.get.coef(Stand.no.creep[[i]])
      if(Creep.dat.same.BC.dat=="YES") a=fn.get.coef(Stand.BaseCase[[i]])
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.bi,pos=a$ID.mnth.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.pos,pos=a$ID.mnth.pos)
      Unik.mn=fn.extract(sort(unique(DATA.list.LIVEWT.c[[i]]$MONTH)))
      new.MN.sens.creep[[i]]=list(bi=Unik.mn[match(dummy$bi,Unik.mn)],pos=Unik.mn[match(dummy$pos,Unik.mn)])
      
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.bi,pos=a$ID.ves.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.pos,pos=a$ID.ves.pos)
      Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c[[i]]$VESSEL)))
      new.level.vessel.sens.creep[[i]]=list(bi=Unik.ves[match(dummy$bi,Unik.ves)],pos=Unik.ves[match(dummy$pos,Unik.ves)])
    }
  }
  
}


#daily
if(Separate.monthly.daily=="YES")
{
  FINYEAR.daily=fn.extract(sort(unique(DATA.list.LIVEWT.c.daily[[1]]$FINYEAR)))
  for ( i in 1:N.species) BLOCKX.daily[[i]]=unique(Stand.BaseCase.daily[[i]]$BiData$BLOCKX)
  new.MN.daily=new.MN
  new.level.vessel.daily=new.level.vessel

  #Approach 1. Most common value of categorical variables
  if(Sel.lev=='most.common')
  {
    #month
    MONTH.daily=MONTH
    MONTH.daily=fn.extract(sort(unique(DATA.list.LIVEWT.c.daily[[1]]$MONTH)))
    for ( i in 1:N.species)
    {
      new.MN.daily[[i]]=Most.common(DATA.list.LIVEWT.c.daily[[i]],"MONTH")
      a=MONTH.daily[as.numeric(new.MN.daily[[i]])]
      new.MN.daily[[i]]=list(bi=a,pos=a)
    }
    
    #vessels    
    for ( i in 1:N.species)
    {
      new.level.vessel.daily[[i]]=Most.common(DATA.list.LIVEWT.c.daily[[i]],"VESSEL")
      Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c.daily[[i]]$VESSEL)))
      a=Unik.ves[match(new.level.vessel[[i]],Unik.ves)]
      new.level.vessel.daily[[i]]=list(bi=a,pos=a)
    }
  }
  
  #Approach 2. Weighted average of month and vessel terms
  if(Sel.lev=='weighted.aver')
  {
    for ( i in 1:N.species)
    {
      a=fn.get.coef(Stand.BaseCase.daily[[i]])
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.bi,pos=a$ID.mnth.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.mnth.pos,pos=a$ID.mnth.pos)
      Unik.mn=fn.extract(sort(unique(DATA.list.LIVEWT.c.daily[[i]]$MONTH)))
      new.MN.daily[[i]]=list(bi=Unik.mn[match(dummy$bi,Unik.mn)],pos=Unik.mn[match(dummy$pos,Unik.mn)])
      
      if(!Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.bi,pos=a$ID.ves.pos)
      if(Bi.same.Pos=="YES") dummy=list(bi=a$ID.ves.pos,pos=a$ID.ves.pos)
      Unik.ves=fn.extract(sort(unique(DATA.list.LIVEWT.c.daily[[i]]$VESSEL)))
      new.level.vessel.daily[[i]]=list(bi=Unik.ves[match(dummy$bi,Unik.ves)],pos=Unik.ves[match(dummy$pos,Unik.ves)])      
    }
  }
  
  
  #covariates
  FREO.daily=soi.daily=KTch.Gum.daily=KTch.Whi.daily=KTch.Dus.daily=
    KTch.San.daily=KTch.Fish.daily=log.Effort.daily=vector(length=N.species)
  for ( i in 1:N.species)
  {
    FREO.daily[i]=Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Freo")
    soi.daily[i]=Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"SOI")
    KTch.Gum.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Catch.Gummy"))
    KTch.Whi.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Catch.Whiskery"))
    KTch.Dus.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Catch.Dusky"))
    KTch.San.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Catch.Sandbar"))
    KTch.Fish.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Catch.Scalefish"))
    log.Effort.daily[i]=log(Mn.cov(DATA.list.LIVEWT.c.daily[[i]],"Km.Gillnet.Days.c"))
  }

}

Grid.data=vector('list',length=N.species)
names(Grid.data)=SPECIES.vec
Grid.data.daily=Grid.data


#Create base case data set for predictions
    #monthly
for ( i in 1:N.species)
{
  if(!i==4) yR.bin=yR.pos=FINYEAR     
  if(i==4)
  {
    yR.bin=FINYEAR.san
    yR.pos=fn.extract(sort(unique(subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0)$FINYEAR)))
  }
    
  if(length(new.MN[[i]]$bi)>0) Mn=new.MN[[i]]$bi
  if(length(new.MN[[i]]$bi)==0) Mn=NA
  dummy.bin=expand.grid(FINYEAR=yR.bin,BLOCKX=BLOCKX[[i]],MONTH=Mn,
      VESSEL=new.level.vessel[[i]]$bi,log.Catch.Gummy=KTch.Gum[i],
      log.Catch.Whiskery=KTch.Whi[i],log.Catch.Dusky=KTch.Dus[i],
      log.Catch.Sandbar=KTch.San[i],log.Catch.Scalefish=KTch.Fish[i],          
      log.Effort=log.Effort[[i]],Freo=FREO[i],SOI=soi[i])

  dummy.pos=expand.grid(FINYEAR=yR.pos,BLOCKX=BLOCKX[[i]],MONTH=new.MN[[i]]$pos,
      VESSEL=new.level.vessel[[i]]$pos,log.Catch.Gummy=KTch.Gum[i],
      log.Catch.Whiskery=KTch.Whi[i],log.Catch.Dusky=KTch.Dus[i],
      log.Catch.Sandbar=KTch.San[i],log.Catch.Scalefish=KTch.Fish[i],          
      log.Effort=log.Effort[[i]],Freo=FREO[i],SOI=soi[i])
  
  Grid.data[[i]]=list(bin=dummy.bin,pos=dummy.pos) 
}
  
    #daily
if(Separate.monthly.daily=="YES")
{
  for ( i in 1:N.species)
  {
     dummy.bin=expand.grid(FINYEAR=FINYEAR.daily,BLOCKX=BLOCKX.daily[[i]],MONTH=new.MN.daily[[i]]$bi,
                          VESSEL=new.level.vessel.daily[[i]]$bi,log.Catch.Gummy=KTch.Gum.daily[i],
                          log.Catch.Whiskery=KTch.Whi.daily[i],log.Catch.Dusky=KTch.Dus.daily[i],
                          log.Catch.Sandbar=KTch.San.daily[i],log.Catch.Scalefish=KTch.Fish.daily[i],          
                          log.Effort=log.Effort.daily[[i]],Freo=FREO.daily[i],SOI=soi.daily[i])
    
    dummy.pos=expand.grid(FINYEAR=FINYEAR.daily,BLOCKX=BLOCKX.daily[[i]],MONTH=new.MN.daily[[i]]$pos,
                          VESSEL=new.level.vessel.daily[[i]]$pos,log.Catch.Gummy=KTch.Gum.daily[i],
                          log.Catch.Whiskery=KTch.Whi.daily[i],log.Catch.Dusky=KTch.Dus.daily[i],
                          log.Catch.Sandbar=KTch.San.daily[i],log.Catch.Scalefish=KTch.Fish.daily[i],          
                          log.Effort=log.Effort.daily[[i]],Freo=FREO.daily[i],SOI=soi.daily[i])

    Grid.data.daily[[i]]=list(bin=dummy.bin,pos=dummy.pos) 
  }
}


#4.23.10 Sensitivity tests
if(do.sensitivity=="YES")
{
  #Indicative vessel sensitivity        
  Grid.data.sens.ind.ves=Grid.data
  for ( i in 1:N.species)
  {
    if(!i==4) yR.bin=yR.pos=FINYEAR
    if(i==4)
    {
      yR.bin=yR.pos=fn.extract(sort(unique(subset(DATA.list.LIVEWT.c.sens.ind.ves[[i]],Catch.Target>0)$FINYEAR)))
    }
        
    if(length(new.MN.sens.ind.ves[[i]]$bi)>0) Mn=new.MN.sens.ind.ves[[i]]$bi
    if(length(new.MN.sens.ind.ves[[i]]$bi)==0) Mn=NA
    dummy.bin=expand.grid(FINYEAR=yR.bin,BLOCKX=BLOCKX.sens.ind.ves[[i]],MONTH=Mn,
      VESSEL=new.level.vessel.sens.ind.ves[[i]]$bi,log.Catch.Gummy=KTch.Gum.sens.ind.ves[i],
      log.Catch.Whiskery=KTch.Whi.sens.ind.ves[i],log.Catch.Dusky=KTch.Dus.sens.ind.ves[i],
      log.Catch.Sandbar=KTch.San.sens.ind.ves[i],log.Catch.Scalefish=KTch.Fish.sens.ind.ves[i],          
      log.Effort=log.Effort.sens.ind.ves[[i]],Freo=FREO.sens.ind.ves[i],SOI=soi.sens.ind.ves[i])
    
    dummy.pos=expand.grid(FINYEAR=yR.pos,BLOCKX=BLOCKX.sens.ind.ves[[i]],MONTH=new.MN.sens.ind.ves[[i]]$pos,
      VESSEL=new.level.vessel.sens.ind.ves[[i]]$pos,log.Catch.Gummy=KTch.Gum.sens.ind.ves[i],
      log.Catch.Whiskery=KTch.Whi.sens.ind.ves[i],log.Catch.Dusky=KTch.Dus.sens.ind.ves[i],
      log.Catch.Sandbar=KTch.San.sens.ind.ves[i],log.Catch.Scalefish=KTch.Fish.sens.ind.ves[i],          
      log.Effort=log.Effort.sens.ind.ves[[i]],Freo=FREO.sens.ind.ves[i],SOI=soi.sens.ind.ves[i])
    
    Grid.data.sens.ind.ves[[i]]=list(bin=dummy.bin,pos=dummy.pos) 
  }
  
  #No effort creep
  Grid.data.no.creep=Grid.data
  log.Effort.no.creep=FREO
  for ( i in 1:N.species)
  {
    log.Effort.no.creep[i]=log(Mn.cov(DATA.list.LIVEWT.c[[i]],"Km.Gillnet.Days.c.no.creep"))
    Grid.data.no.creep[[i]]$bin$log.Effort=log.Effort.no.creep[i]      
    Grid.data.no.creep[[i]]$pos$log.Effort=log.Effort.no.creep[i]  
    
    if(length(new.MN.sens.creep[[i]]$bi)>0) Mn=new.MN.sens.creep[[i]]$bi
    if(length(new.MN.sens.creep[[i]]$bi)==0) Mn=NA
    Grid.data.no.creep[[i]]$bin$MONTH=Mn 
    Grid.data.no.creep[[i]]$bin$VESSEL=new.level.vessel.sens.creep[[i]]$bi
    
    Grid.data.no.creep[[i]]$pos$MONTH=new.MN.sens.creep[[i]]$pos
    Grid.data.no.creep[[i]]$pos$VESSEL=new.level.vessel.sens.creep[[i]]$pos
  }
  
}


#4.23.11 Imputation of missing Year-block 
#rules: (combination of Punt et al 2000 and Carruthers et al 2011)
        #1.   If NAs prior first record, then fill in with average first 3 records with data
        #2.   If NAs after last record, see "Second.rule"
        #3.   If NAs in between records, then use linear interpolation

#note: although whiskery shark shows strong monthly trends, this is already accounted for when adding
#       Month as a model term. As a fixed Month value is used for predicting catch (i.e. the effect of
#       Month is removed by the standardisation), the imputation is not affected by having no observation
#       during the seasonal closure



#4.23.11.1 Construct index 

#4.23.11.2 Compare indices in sensitivity test
Pred.st.BaseCase=vector('list',length=N.species)
names(Pred.st.BaseCase)=SPECIES.vec
Pred.st.Fantasy=Pred.st.no.creep=Pred.st.ind.ves.sens=Pred.st.no.area.weight=
  Pred.st.BaseCase.daily=Pred.st.Fantasy.daily=Pred.st.no.interaction=Pred.st.BaseCase

Yr_Blks.imputed=Prop.yr_blks.imputed=Prop.yr_blks.imputed.daily=Pred.st.BaseCase

#4.23.11.3 Construct indices for sensitivity scenarios
system.time(for (i in 1:N.species)
{
  #Base case
  
    #1. predict catch
  aa=predict.ktch(Stand.BaseCase[[i]],Grid.data[[i]],SPECIES.vec[[i]])        
  Pred.BaseCase=aa
  
    #2. impute missing year-blocks
  Pop.growth=POP.GRW.RATE[i]
  Pred.BaseCase$Pos.dat=impute.ktch(Pred.BaseCase$Pos.dat) 
  ss=subset(Pred.BaseCase$Pos.dat,is.na(Pred.Catch))
  Prop.yr_blks.imputed[[i]]=nrow(ss)/nrow(Pred.BaseCase$Pos.dat)
  Yr_Blks.imputed[[i]]=ss[,match(c("FINYEAR","BLOCKX"),names(ss))]
  
    #3. plot example of how imputation is done
  this.impted.blks=unique(Pred.BaseCase$Pos.dat$BLOCKX)
  DaTT=DATA.list.LIVEWT.c[[i]]  
  DaTT$cpue=DaTT$Catch.Target/DaTT$Km.Gillnet.Days.c 
  COF=coef(Stand.BaseCase[[i]]$LogN)
  Cof.mn=COF[match(paste("MONTH",2:12,sep=""),names(COF))]
  names(COF)=substr(names(COF),7,20)
  Cof.mn=sort(Cof.mn)
  names(Cof.mn)=substr(names(Cof.mn),6,7)  
  for(Q in 1:length(this.impted.blks))
  {
    #show predicted catch
    a=subset(Pred.BaseCase$Pos.dat,BLOCKX==this.impted.blks[Q])
    if(nrow(a)>0)
    {
      tiff(paste("Example of how imputation works/",SPECIES.vec[i],"/",this.impted.blks[Q],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      par(mgp=c(2.5,.7,0),las=1)
      x=1:nrow(a)
      Ylow=min(a$Pred.Catch-1.96*a$SE,na.rm=T)
      Yhigh=max(a$Pred.Catch+1.96*a$SE,na.rm=T)
      plot(x,a$Pred.Catch.imp,pch=1,col=2,cex=2,xlab="Finyear",ylab="Log of catch  +/-95% CI",xaxt='n',
           cex.lab=2,ylim=c(Ylow,Yhigh))
      axis(1,x,a$FINYEAR)
      points(x,a$Pred.Catch,pch=19,col=1,cex=2)
      segments(x,a$Pred.Catch,x,a$Pred.Catch+1.96*a$SE)
      segments(x,a$Pred.Catch,x,a$Pred.Catch-1.96*a$SE)
      legend("topright",c("predicted by model","imputed"),bty="n",col=c(1,2),pch=c(19,1),cex=2)
      Monthly.yrs=x[1:31]
      Ylow=rep(Ylow,length(Monthly.yrs))
      Yhigh=rep(Yhigh,length(Monthly.yrs))
      XX=c(Monthly.yrs,tail(Monthly.yrs, 1),rev(Monthly.yrs),Monthly.yrs[1])
      YY <- c(Ylow, tail(Yhigh, 1), rev(Yhigh), Ylow[1])
      polygon(XX,YY,col=rgb(0,.3,.1,alpha=0.2),border=rgb(0,.3,.1,alpha=0.2))
      
      dev.off()
    }   
    
    #show data
    b=subset(DaTT,BLOCKX==this.impted.blks[Q])      
    if(nrow(b)>0)
    {
      tiff(file=paste("Example of how imputation works/",SPECIES.vec[i],"/","raw.",this.impted.blks[Q],".tiff",sep=""),
           width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")          
      par(mfcol=c(3,1),mai=c(.5,.5,.15,.1),mgp=c(2.5,.5,0))              
      boxplot(Catch.Target~FINYEAR,b,main=this.impted.blks[Q],ylab="Catch",cex.axis=1,cex.lab=1.75,col="grey")
      boxplot(Km.Gillnet.Days.c~FINYEAR,b,ylab="Km.Gillnet.Days.c",cex.axis=1,cex.lab=1.75,col="grey")
      boxplot(cpue~FINYEAR,b,ylab="cpue",cex.axis=1,cex.lab=1.75,col="grey")
      dev.off()   
      
      #vessels for that block
      dis.vs=as.character(unique(b$VESSEL))
      TBl.ves.yr=with(b,table(FINYEAR,VESSEL))
      ID=match(dis.vs,names(COF))
      Ves.CF=sort(COF[ID])        
      TBl.ves.yr=as.data.frame.matrix(TBl.ves.yr)
      TBl.ves.yr=TBl.ves.yr[,match(names(Ves.CF),colnames(TBl.ves.yr))]
      colnames(TBl.ves.yr)=paste(colnames(TBl.ves.yr),".coef.",round(Ves.CF,2),sep="")
      write.csv(TBl.ves.yr,paste("Example of how imputation works/",SPECIES.vec[i],"/",this.impted.blks[Q],"_Vessels.csv",sep=""))
      
      #months for that block
#       dis.mn=as.character(unique(b$MONTH))
#       TBl.mn.yr=with(b,table(FINYEAR,MONTH))        
#       TBl.mn.yr=as.data.frame.matrix(TBl.mn.yr)
#       Cof.mn=
#       xxx=match(names(Cof.mn),colnames(TBl.mn.yr))
#       xxx=xxx[!is.na(xxx)]
#       TBl.mn.yr=TBl.mn.yr[,xxx]
#       colnames(TBl.mn.yr)=paste("Month.",colnames(TBl.mn.yr),".coef.",round(Cof.mn,2),sep="")
#       write.csv(TBl.mn.yr,paste("Example of how imputation works/",SPECIES.vec[i],"/",this.impted.blks[Q],"_Month.csv",sep=""))
      
    }      
  }
  
    #4. construct index   
  Pred.st.BaseCase[[i]]=construct.index(Pred.BaseCase,AREA.W[[i]],"YES",SPECIES.vec[[i]],
                          Weight.by.area="YES")  
  
  #Fantasy
    #1. predict catch
  Pred.Fantasy=aa
  
    #2. construct index
  Pred.st.Fantasy[[i]]=construct.index(Pred.Fantasy,AREA.W[[i]],IMPUTE="NO",SPECIES.vec[[i]],
                                       Weight.by.area="YES")
  
  rm(Pred.BaseCase,Pred.Fantasy)
  
  if(do.sensitivity=="YES") 
  {
    Pred.ind.ves=Stand.ind.vess[[i]]
    Pred.no.creep=Stand.no.creep[[i]]
    Pred.no.weight=Stand.BaseCase[[i]]
    #Pred.no.interaction=Stand.no.interaction[[i]]
   
      #1. predict catch
    Pred.ind.ves=predict.ktch(Pred.ind.ves,Grid.data.sens.ind.ves[[i]],SPECIES.vec[[i]])
    Pred.no.weight=predict.ktch(Pred.no.weight,Grid.data[[i]],SPECIES.vec[[i]])
    Pred.no.creep=predict.ktch(Pred.no.creep,Grid.data.no.creep[[i]],SPECIES.vec[[i]])    
    #Pred.no.interaction=predict.ktch(Pred.no.interaction,Grid.data[[i]],SPECIES.vec[[i]])
    
      #2. impute missing year-blocks
    Pred.ind.ves$Pos.dat=impute.ktch(Pred.ind.ves$Pos.dat)
    Pred.no.creep$Pos.dat=impute.ktch(Pred.no.creep$Pos.dat)
    Pred.no.weight$Pos.dat=impute.ktch(Pred.no.weight$Pos.dat)
    #Pred.no.interaction$Pos.dat=impute.ktch(Pred.no.interaction$Pos.dat)
    
    #3. construct index
    Pred.st.ind.ves.sens[[i]]=construct.index(Pred.ind.ves,AREA.W[[i]],IMPUTE="YES",SPECIES.vec[[i]],
                                  Weight.by.area="YES")  
    Pred.st.no.creep[[i]]=construct.index(Pred.no.creep,AREA.W[[i]],IMPUTE="YES",SPECIES.vec[[i]],
                                 Weight.by.area="YES")    
    Pred.st.no.area.weight[[i]]=construct.index(Pred.no.weight,AREA.W[[i]],IMPUTE="YES",SPECIES.vec[[i]],
                                Weight.by.area="NO")
   # Pred.st.no.interaction[[i]]=construct.index(Pred.no.interaction,AREA.W[[i]],"YES",SPECIES.vec[[i]],
   #                                             Weight.by.area="YES")
      
    rm(Pred.ind.ves,Pred.no.weight,Pred.no.creep)
  }

    #daily
  if(Separate.monthly.daily=="YES")
  {
     
    #1. predict catch
    aa=predict.ktch(Stand.BaseCase.daily[[i]],Grid.data.daily[[i]],SPECIES.vec[[i]])         
    
    Pred.BaseCase=aa
    Pred.Fantasy=aa
    
    
    #2. impute missing year-blocks
    Prop.yr_blks.imputed.daily[[i]]=nrow(subset(Pred.BaseCase$Pos.dat,is.na(Pred.Catch)))/nrow(Pred.BaseCase$Pos.dat)
    Pred.BaseCase$Pos.dat=impute.ktch(Pred.BaseCase$Pos.dat) 
    
    #3. construct index   
    Pred.st.BaseCase.daily[[i]]=construct.index(Pred.BaseCase,AREA.W[[i]],"YES",SPECIES.vec[[i]],Weight.by.area="YES")  
    Pred.st.Fantasy.daily[[i]]=construct.index(Pred.Fantasy,AREA.W[[i]],IMPUTE="NO",SPECIES.vec[[i]],Weight.by.area="YES")
  
    rm(Pred.BaseCase,Pred.Fantasy)    

  }
})

write.csv(do.call(rbind,Prop.yr_blks.imputed),"Proportion.yr-blks.imputed.csv")
write.csv(do.call(rbind,Prop.yr_blks.imputed.daily),"Proportion.yr-blks.imputed.daily.csv")


  #Compare indices

if(do.colors=="YES")
{
  COLORES=c("coral","coral4","chocolate1","red","darkblue","chartreuse4","black") 
  LWD=c(2.5,2.5,3,3,3,3.5,2.75)
  LINEAS=c(3,3,2,2,2,1,1)
}
  
if(!do.colors=="YES")
{
  if(do.sensitivity=="YES")
  {
    COLORES=c("black","grey60","grey30","grey60","grey45","grey65","black")
    LWD=c(2.5,2.5,3,3,3,4,4)
    LINEAS=c(3,3,2,4,4,1,1)
  }

  if(Separate.monthly.daily=="YES")
  {
    COLORES=c("grey10","grey40","grey60","black","grey60","black")
    names(COLORES)=c("Nominal","Folly","Fan","BaseCase","Fan.daily","BaseCase.daily")
    LWD=c(3,3,3,3,3,3)
    LINEAS=c(2,3,1,1,1,1)    
  }
}


setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
if(Separate.monthly.daily=="YES")
{
  #Plot all the indices
  where=c('topright','topright','topright','topright')
  
  #by species
  for ( i in 1:N.species)
  {
    hndl=paste("Figure 10.Cpue.comparison.",SPECIES.vec[i],".tiff",sep="")
    tiff(file=hndl,width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(1,1),mar=c(1.5,3.5,1,.1),oma=c(2.5,.2,.1,.2),las=1,mgp=c(1.9,.55,0))
    Plot.Comb.cpue.daily(BC=Pred.st.BaseCase[[i]],Fan=Pred.st.Fantasy[[i]],
                         BC.daily=Pred.st.BaseCase.daily[[i]],Fan.daily=Pred.st.Fantasy.daily[[i]],
                         REL=RELATIVE)
    legend('topright',c("Nominal","Folly","Fantasy","Base case"),col=COLORES,lty=LINEAS,bty='n',cex=1.75,lwd=LWD)
    mtext("Financial year",side=1,line=2.75,font=1,las=0,cex=2.5)
    if(RELATIVE=="YES") mtext("Relative catch rate",side=2,line=2.05,font=1,las=0,cex=2)
    if(RELATIVE=="NO") mtext("Catch rate (kg/km.gn.day)",side=2,line=2.2,font=1,las=0,cex=2)
    dev.off()
  }
  
  #species combined
  tiff(file="Figure 10.Cpue.comparison.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,3.75,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
  for ( i in 1:N.species)
  {
    Plot.Comb.cpue.daily(BC=Pred.st.BaseCase[[i]],Fan=Pred.st.Fantasy[[i]],
                         BC.daily=Pred.st.BaseCase.daily[[i]],Fan.daily=Pred.st.Fantasy.daily[[i]],
                         REL=RELATIVE)
    if(i==2)legend('top',c("Nominal","Folly","Fantasy","Base case"),col=COLORES,
                   lty=LINEAS,bty='n',cex=1.25,lwd=LWD)
    legend("bottomleft",SPECIES.vec[i],bty='n',cex=1.4)
  }
  mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.75,outer=T)
  if(RELATIVE=="YES") mtext("Relative catch rate",side=2,line=-1.4,font=1,las=0,cex=1.75,outer=T)
  if(RELATIVE=="NO") mtext("Catch rate (kg/km.gn.day)",side=2,line=-1.35,font=1,las=0,cex=1.75,outer=T)
  dev.off() 
  

}


if(do.sensitivity=="YES")
{
  #Plot all the indices
  where=c('topright','topright','topright','topright')
  
    #by species
  for ( i in 1:N.species)
  {
      hndl=paste("Figure 10.Cpue.comparison.",SPECIES.vec[i],".tiff",sep="")
      tiff(file=hndl,width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      par(mfcol=c(1,1),mar=c(1.5,3.5,1,.1),oma=c(2.5,.2,.1,.2),las=1,mgp=c(1.9,.55,0))
      Plot.Comb.cpue(Pred.st.BaseCase[[i]],Pred.st.Fantasy[[i]],Pred.st.no.creep[[i]],
                     Pred.st.ind.ves.sens[[i]],Pred.st.no.area.weight[[i]],
                     SPECIES.vec[i],REL=RELATIVE,where[i],Ves.no.ktc[[i]],Min.weight[i])
      legend('topright',c("Nominal","Folly","Fantasy","No effort creep","Ind. vessels sens.",
                          "No area weighting","Base case"),col=COLORES,lty=LINEAS,bty='n',cex=1.75,lwd=LWD)
      mtext("Financial year",side=1,line=2.75,font=1,las=0,cex=2.5)
      if(RELATIVE=="YES") mtext("Relative catch rate",side=2,line=2.05,font=1,las=0,cex=2)
      if(RELATIVE=="NO") mtext("Catch rate (kg/km.gn.day)",side=2,line=2.2,font=1,las=0,cex=2)
    dev.off()
  }

    #species combined
  tiff(file="Figure 10.Cpue.comparison.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfcol=c(2,2),mar=c(1,3.75,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,.7,0))
    for ( i in 1:N.species)
    {
      Plot.Comb.cpue(Pred.st.BaseCase[[i]],Pred.st.Fantasy[[i]],Pred.st.no.creep[[i]],
                     Pred.st.ind.ves.sens[[i]],Pred.st.no.area.weight[[i]],
                     SPECIES.vec[i],REL=RELATIVE,where[i],Ves.no.ktc[[i]],Min.weight[i])
      if(i==2)legend('top',c("Nominal","Folly","Fantasy","No effort creep","Ind. vessels sens.",
                    "No area weighting","Base case"),col=COLORES,
                     lty=LINEAS,bty='n',cex=1.25,lwd=LWD)
      legend("bottomleft",SPECIES.vec[i],bty='n',cex=1.4)
    }
    mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.75,outer=T)
    if(RELATIVE=="YES") mtext("Relative catch rate",side=2,line=-1.4,font=1,las=0,cex=1.75,outer=T)
    if(RELATIVE=="NO") mtext("Catch rate (kg/km.gn.day)",side=2,line=-1.35,font=1,las=0,cex=1.75,outer=T)
  dev.off() 
  
  #Compare no effort creep with base case

  if(compare.creep.sep=="YES")
  {
    COLORES=c("red","black")
    LINEAS=c(1,1)
    where=c("topright","topright","topright")
    
    tiff(file="Figure 9.Cpue.effort.creep.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mfcol=c(3,1),mar=c(2,5.5,1,.1),oma=c(2.5,.1,.1,.2),las=1,mgp=c(1.9,1,0))
    for ( i in 1:3)
    {
      Plot.Comb.cpue.eff.creep(Pred.st.BaseCase[[i]],Pred.st.no.creep[[i]],SPECIES.vec[i],where[i])
      if(i==1)legend('right',c("No effort creep","Base case"),col=COLORES,lty=LINEAS,bty='n',cex=1.75,lwd=LWD)
    }
    mtext("CPUE (kg/km.gn.day)",side=2,line=-2,font=1,las=0,cex=1.75,outer=T)
    mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.75,outer=T)
    dev.off()
  }
}


#4.23.12 Plot sandbar spatio-temporal catches
#note:  no point in choosing blocks that are permanently closed (Metro closure: 3115, 3215, 3214)
if(Do.sandbar.explore=="YES")
  {
  #Plot spatial sandbar catches 
  i=4
  a=subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0)
  #a$cpue=a$Catch.Target/a$Km.Gillnet.Days.c
  a$LAT=a$LAT-.5
  a$LONG=a$LONG+.5
  yr=sort(unique(a$FINYEAR))
  n=aggregate(Catch.Target~FINYEAR+LONG+LAT,a,sum)
  ss=reshape(n,idvar=c("FINYEAR","LONG"),timevar="LAT",v.names="Catch.Target", direction="wide")  
  LoNgS=sort(unique(a$LONG))
  numInt=10
  colfunc <- colorRampPalette(c("aquamarine","cadetblue"))
  couleurs=colfunc(numInt)
  #couleurs=rev(gray(seq(0.2,0.9,length=numInt)))
  #couleurs=(heat.colors(numInt))
  couleurs=c("darkolivegreen1","darkolivegreen2","darkolivegreen3","chartreuse2","chartreuse3","chartreuse4",
             "forestgreen","darkolivegreen4","darkolivegreen","darkgreen")
  

  tiff(file="Sandbar.Spatial.Catch.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
  par(mfcol=c(5,6),mai=c(.1,.1,.1,.1),oma=c(3,3,.1,.1),las=1,mgp=c(1,.5,0))
  for(p in 1:length(yr))
  {
    aa=subset(a,FINYEAR==yr[p])
    #BRKS=quantile(aa$Catch.Target,probs=seq(0,1,1/numInt),na.rm=T)
    BRKS=seq(0,1,length.out=numInt+1)
    dum=subset(ss,FINYEAR==yr[p])
    if(length(dum$LONG)<length(LoNgS))
    {
      miss=LoNgS[which(!LoNgS%in%dum$LONG)]
      Add=dum[1:length(miss),]
      Add[,3:ncol(Add)]=NA
      Add$LONG=miss
      dum=rbind(dum,Add)
    }
    
    dum=dum[order(dum$LONG),]
    long=dum$LONG  
    dum=as.matrix(dum[,-(1:2)])
    colnames(dum)=substr(colnames(dum),14,20)
    lat=as.numeric(colnames(dum))
    dum=dum/sum(dum,na.rm=T)   #as a proportion of annual catch
    
    image(long,lat,z=dum,xlab="",ylab="",col=couleurs,breaks=BRKS)
    #axis(side = 1, at =113:118, labels = F, tcl = .1)
    #axis(side = 2, at = -36:-27, labels = F,tcl =.1)
    box()
    legend("bottomleft",yr[p],bty='n',cex=1)
  }
  mtext("Longitude",side=1,line=1.5,font=1,las=0,cex=1.5,outer=T)
  mtext("Latitude",side=2,line=1.25,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  sort(table(a$BLOCKX))
  
}


#4.23.13 Explain the influence of the different model terms
#note: implements Bentley et al 2012, "Influence plots"

numInt=10
numberLab=2
if(do.colors=="YES")couleurs <- colorRampPalette(c("darkseagreen1","chartreuse4"))(numInt)
if(do.colors=="NO")couleurs=rev(gray(seq(0.2,0.9,length=numInt)))  

if(do.influence=="YES")
{
    #3.11.1 proportion of imputed blocks by year
  tiff(file="Appendix 6.Imputed.blks_by_yr.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1.25,4,1.5,4.25),oma=c(2.5,.8,.8,.1),las=1,mgp=c(1.9,.5,0))
  for(i in 1:N.species) fn.prop.imputed(Grid.data[[i]]$pos,Yr_Blks.imputed[[i]],SPECIES.vec[i],Stand.BaseCase[[i]]$LogN)
  #par(mfcol=c(3,1),mar=c(2,4,1.5,.5),oma=c(2.5,2,.8,.1),las=1,mgp=c(1.9,.5,0))
  #for(i in 1:3) fn.prop.imputed(Grid.data[[i]]$pos,Yr_Blks.imputed[[i]],SPECIES.vec[i],Stand.BaseCase[[i]]$LogN)
  mtext("Financial year",1,outer=T,line=1.25,cex=2.25)
  mtext("Block",2,outer=T,las=3,line=-1,cex=2.25)
  dev.off()
  
  #numbers of year-blocks imputed
  N.imp=Yr_Blks.imputed
  #for(i in 1:N.species) N.imp[[i]]=nrow(Yr_Blks.imputed[[i]])
  for(i in 1:3) N.imp[[i]]=nrow(Yr_Blks.imputed[[i]])


  #3.11.2 Influence plots
  #       postive catch only
  #terms for each species (from Best model)
  Terms=vector('list',N.species)
  names(Terms)=names(DATA.list.LIVEWT.c)
  Term.Type=Names.vec=Terms
  for(i in 1:N.species)
  {
    a=gsub(" ", "", unlist(strsplit(as.character(Best.Model[[i]]$Log)[3],"[+]")), fixed = TRUE)
    IDD=match("offset(log.Effort)",a)
    if(length(IDD)>0) a=a[-IDD]
    
    if(!(i==4 & Interaction.sandbar=="NO")) a[1]=unlist(strsplit(a[1],"[*]"))[2]
    if((i==4 & Interaction.sandbar=="NO")) a=a[-1]
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
  
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Influence.plot")
  
  XLIM.vec=list(c(.5,2.5),c(.6,1.4),c(.5,1.75),c(.5,1.6))
  SCALE.vec=list(c(25,30,20),c(25,30,20),c(20,30,20),c(20,15,50))
  
  
  Store.Influence=vector('list',N.species)
  names(Store.Influence)=SPECIES.vec
  
  #Figure bubble plots
  for ( i in 1:N.species)
  {
    Store.Influence[[i]]=Influence.fn(Stand.BaseCase[[i]]$LogN,
              subset(DATA.list.LIVEWT.c[[i]],!(VESSEL%in%Ves.no.ktc[[i]])),
              Stand.BaseCase[[i]]$PosData,
              Term.Type[[i]],Terms[[i]],Names.vec[[i]],XLIM.vec[[i]],SCALE.vec[[i]])
  }
  
  #Figure 11
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
  #XLIM.vec=list(c(.65,1.45),c(.65,1.45),c(.65,1.45),c(.65,1.45))
  #Compare influence of all terms
  tiff("Figure 11.All.terms.Influence.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,2),mar=c(1.5,2.5,.5,.1),oma=c(2,2,.1,.6),las=1,mgp=c(1.9,.6,0))
  #par(mfcol=c(3,1),mar=c(.1,.1,2,3),oma=c(5,6,.5,.5),las=1,mgp=c(1.9,.7,0))
  for ( i in 1:N.species)
  {
    Annual.Dev=Store.Influence[[i]]$Annual.Dev
    ny=Store.Influence[[i]]$ny
    nt=length(Annual.Dev)
     
    LWD=3
    if(do.colors=="YES")
    {
      LTY.col=c("red","chartreuse4","darkblue","chocolate1")
      LTY=rep(1,nt)
    }
      
    if(do.colors=="NO")
    {
      LTY.col=c("black","grey60","black","grey45")
      LTY=1:nt
    }
      
    plot(1:length(ny),Annual.Dev[[1]],col=LTY.col[1],type="l",xlab="",ylab="",lwd=LWD,cex.axis=1.35,xaxt='n',ylim=XLIM.vec[[i]])
    abline(h=1,lty=3,col=1)
    axis(1,1:length(ny),F,tck=-0.015)
    axis(1,seq(1,length(ny),5),F,tck=-0.03)
    axis(1,seq(1,length(ny),5),names(ny)[seq(1,length(ny),5)],cex.axis=1.35,tck=-0.03)
    for(p in 2:nt)lines(1:length(ny),Annual.Dev[[p]],lwd=LWD,lty=LTY[p],col=LTY.col[p])
    LEG=paste(Names.vec[[i]]," (",100*round(Store.Influence[[i]]$Over.all.influence,2),"%)",sep="")
    legend("topright",LEG,bty='n',lty=LTY,col=LTY.col,lwd=LWD,cex=1.20)
    legend("bottomright",SPECIES.vec[i],bty='n',cex=1.65)
  }
  mtext("Influence",side=2,line=0,cex=2.25,las=3,outer=T)
  mtext("Financial year",side=1,line=0.9,cex=2.25,las=1,outer=T)
  dev.off()
  
}



#4.23.14  Build confidence intervals and abundance index for stock assessments
#note: it takes 26 hours for 4 species, 1000 boots and annual and monthly indices
system.time(if(Do.boot=="YES")
{
  #note:  Boot.method="DATA" resamples original data in a stratified way to keep all year-blocks (Aires da Silva et al 2008)
  #       Boot.method="RESIDUALS uses Haddon 2001 (page 306) resampling of residuals approach
  
    #get blocks per zone
  A=c(DATA.list.LIVEWT.c[[1]]$BLOCKX,DATA.list.LIVEWT.c[[2]]$BLOCKX,DATA.list.LIVEWT.c[[3]]$BLOCKX,DATA.list.LIVEWT.c[[4]]$BLOCKX)
  B=c(DATA.list.LIVEWT.c[[1]]$zone,DATA.list.LIVEWT.c[[2]]$zone,DATA.list.LIVEWT.c[[3]]$zone,DATA.list.LIVEWT.c[[4]]$zone)
  D=as.data.frame.matrix(table(A,B))
  BLKZ.per.zone=list(West=as.numeric(rownames(subset(D,West>0))),
                     Zone1=as.numeric(rownames(subset(D,Zone1>0))),
                     Zone2=as.numeric(rownames(subset(D,Zone2>0))))
 
  #Extract variables to bootstrap
  common.vars=c("Catch.Target","Km.Gillnet.Days.c","BLOCKX","VESSEL","MONTH","FINYEAR",
                "Catch.Gummy","Catch.Whiskery","Catch.Dusky","Catch.Sandbar")
  VARS=Grid.data
  for(i in 1:N.species)VARS[[i]]=common.vars

  #Set up data sets for predicting index
    #monthly
  Grid.data.zn.mnth=Grid.data
  All.mns=factor(1:12,levels=1:12)
  for (o in 1:N.species) Grid.data.zn.mnth[[o]]=fn.build.pred.mn.dat.boot(Grid.data.zn.mnth[[o]])

    #extract blockx whithin each zone for each species
  ZONES=list(Whis=names(Zones.whi),Gum=names(Zones.gum),Dus=names(Zones.dus),San="West")
  Blks.in.zn=vector('list',N.species)
  names(Blks.in.zn)=SPECIES.vec
  for ( i in 1:N.species) Blks.in.zn[[i]]=fn.blk.zn(DATA.list.LIVEWT.c[[i]])
  
  #get the blocks used in the standardisations
  for(i in 1:N.species)
  {
    BLkS=Blks.in.zn[[i]]
    hndl="C:/Matias/Analyses/Catch and effort/Outputs/Blocks_kept/"
     for(b in 1:length(BLkS))
     {
       x=BLkS[[b]]
#       LAT=-as.numeric(substr(x,1,2))
#       LONG=100+as.numeric(substr(x,3,4))
        write.csv(x,paste(hndl,SPECIES.vec[i],"_",names(BLkS)[b],".csv",sep=""),row.names=F)
     }
  }
  
  
  #3.12.1 create boot data
  
    #3.12.1.1 resampling of residuals approach
  #note: this approach is fast but only boostraps the positive catch
  if(Boot.method=="RESIDUALS")
  {
    Res.type="additive"   
    #Res.type="multiplicative"
    Bootstrapped.Dat=vector('list',length=N.species)
    names(Bootstrapped.Dat)=SPECIES.vec
    for (i in 1:N.species)
    {
      MOD=Stand.BaseCase[[i]]
      DAT=Stand.BaseCase[[i]]$PosData
      vars=VARS[[i]]
      Bootstrapped.Dat[[i]]=Boot.fn.res(N.boot)  
    }
  }
  
  #remove too heavy object
  if(exists(c("Stand.no.creep","Data.Summary","DATA.list.LIVEWT.c.daily",
              "DATA.list.LIVEWT.c.daily.sens","Indi.ves.daily",
              "Store.DATA.list.LIVEWT.c.daily"))) 
  {
    rm(Stand.no.creep,Data.Summary,DATA.list.LIVEWT.c.daily,
       DATA.list.LIVEWT.c.daily.sens,Indi.ves.daily,Store.DATA.list.LIVEWT.c.daily)    
  }

     #3.12.1.2 run bootstraps
    #create objects for storing boots
    Pred.st.BaseCase.CI=vector('list',length=N.species)
    names(Pred.st.BaseCase.CI)=SPECIES.vec
    Pred.st.BaseCase.CI.by.zone=Pred.st.BaseCase.CI.Annual.by.zone=Pred.st.BaseCase.CI
    
    for ( i in 1:N.species)
    {
      BiData=Stand.BaseCase[[i]]$BiData
      
      AA=fn.bootstrap.index(DATA.list.LIVEWT.c[[i]],SPECIES.vec[[i]],
                            VARS[[i]],Best.Model[[i]],Grid.data[[i]],Grid.data.zn.mnth[[i]],Ves.no.ktc[[i]],
                            AREA.W[[i]],ZONES[[i]],Blks.in.zn[[i]])
      
      #Annual (zones combined)
      Pred.st.BaseCase.CI[[i]]=AA$Store.boot.BaseCase 
      
      #Annual (by zone)
      Pred.st.BaseCase.CI.Annual.by.zone[[i]]=AA$Store.boot.BaseCase.annual.by.zone
      
      #Monthly (by zone)
      Pred.st.BaseCase.CI.by.zone[[i]]=AA$Store.boot.BaseCase.zone
    }
  

  #3.12.2 extract statistics and plot indices
  
  #extract mean, CV and confidence intervals  
  Index.BaseCase.CI=vector('list',N.species)
  names(Index.BaseCase.CI)=names(Pred.st.BaseCase.CI)
  Index.BaseCase.CI.by.zone=Index.BaseCase.CI.annual.by.zone=Index.BaseCase.CI
  
    #Annual base case zones combined          
  for(i in 1:N.species)Index.BaseCase.CI[[i]]=fn.get.mean(Pred.st.BaseCase.CI[[i]],"ZonesCombined")

    #Annual base case by zone
  for(i in 1:N.species)Index.BaseCase.CI.annual.by.zone[[i]]=fn.get.mean(Pred.st.BaseCase.CI.Annual.by.zone[[i]],"AnnualByZones")
  
    #Monthly base case by zone
  for(i in 1:N.species)Index.BaseCase.CI.by.zone[[i]]=fn.get.mean(Pred.st.BaseCase.CI.by.zone[[i]],"ByZones")
  
  #plotting 

  if (what.color=="cols")
  {
    COL="dodgerblue4"
    Colr="lightsteelblue1"
    Colr2="lightsteelblue4"
#     COL="forestgreen"
#     Colr="olivedrab2"
#     Colr2="darkgreen"
    
  }
  if (what.color=="black")
  {
    COL="black"
    Colr="grey80"
    Colr2="grey20"
  }


  do.axis=c("YES","YES","YES","YES")
  
  #Plot annual index
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
  tiff(file="Figure 13.Annual.index.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1.5,1.5,1,1),oma=c(2.5,2.5,.1,.1),las=1,mgp=c(.1,.7,0))
  for ( i in 1:N.species)
  {
    a=Index.BaseCase.CI[[i]]
    Plot.Index(a$MEAN,a$LOW.CI,a$UP.CI,a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YEAR",do.axis[i])
    legend('topright',SPECIES.vec[i],bty='n',cex=1.75)
  }
  mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=2.5,outer=T)
  if(RELATIVE=="YES") mtext("Relative catch rate",side=2,line=0,font=1,las=0,cex=2.5,outer=T)
  if(RELATIVE=="NO") mtext("CPUE (kg/km.gn.day)",side=2,line=1,font=1,las=0,cex=2.5,outer=T)
  dev.off()
  
  #Plot annual by zone index
  do.axis.zn=rep("YES",12)
  ZONAS=c("West","Zone1","Zone2")
  Blnk=function()
  {
    plot(1:10,xaxt='n',yaxt='n',ann=F,col="transparent")
    box(col="white")
  }
  tiff(file="Figure 13.Annual.by.zone.index.tiff",width = 2800, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(4,3),mar=c(1,1.25,1.5,1),oma=c(2.5,3,.1,.5),las=1,mgp=c(1.9,.55,0))
  for ( i in 1:N.species)
  {
    a=Index.BaseCase.CI.annual.by.zone[[i]]
    zns=names(a$MEAN)
    names(zns)="OK"
    ID=which(!ZONAS%in%zns)
    missing.zns=ZONAS[ID]
    if(length(missing.zns)>0)
    {
      if(i==2)
      {
        Blnk();Blnk()
        Plot.Index(a$MEAN[[1]],a$LOW.CI[[1]],a$UP.CI[[1]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YEAR",do.axis.zn[i])
      }
      if(i==4)
      {
        Plot.Index(a$MEAN[[1]],a$LOW.CI[[1]],a$UP.CI[[1]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YEAR",do.axis.zn[i])
        Blnk();Blnk()
      }
      
    }
    if(length(missing.zns)==0)
    {
      for(z in 1:length(zns)) 
      {
        Plot.Index(a$MEAN[[z]],a$LOW.CI[[z]],a$UP.CI[[z]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YEAR",do.axis.zn[i])
        if(i==1)mtext(zns[z],side=3,line=.15,font=1,las=0,cex=1.25,outer=F)
      }
    }
    #legend('topright',SPECIES.vec[i],bty='n',cex=1.5)
    mtext(SPECIES.vec[i],4,line=0.5,cex=1.25,las=3)
  }
  mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.5,outer=T)
  if(RELATIVE=="YES") mtext("Relative cpue",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
  if(RELATIVE=="NO") mtext("CPUE (kg/km.gn.day)",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
  dev.off()


  #Plot monthly by zone index
  tiff(file="Figure 13.Monthly.by.zone.index.tiff",width = 2800, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(4,3),mar=c(1,1.25,1.5,1),oma=c(2.5,3,.1,.5),las=1,mgp=c(1.9,.55,0))
  for ( i in 1:N.species)
  {
    a=Index.BaseCase.CI.by.zone[[i]]
    zns=names(a$MEAN)
    names(zns)="OK"
    ID=which(!ZONAS%in%zns)
    missing.zns=ZONAS[ID]
    if(length(missing.zns)>0)
    {
      if(i==2)
      {
        Blnk();Blnk()
        Plot.Index(a$MEAN[[1]],a$LOW.CI[[1]],a$UP.CI[[1]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YR.MN",do.axis.zn[i])
      }
      if(i==4)
      {
        Plot.Index(a$MEAN[[1]],a$LOW.CI[[1]],a$UP.CI[[1]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YR.MN",do.axis.zn[i])
        Blnk();Blnk()
      }
      
    }
    if(length(missing.zns)==0)
    {
      for(z in 1:length(zns)) 
      {
        Plot.Index(a$MEAN[[z]],a$LOW.CI[[z]],a$UP.CI[[z]],a$FINYEAR,COL,Colr,Colr2,RELATIVE,"YR.MN",do.axis.zn[i])
        if(i==1)mtext(zns[z],side=3,line=.15,font=1,las=0,cex=1.25,outer=F)
      }
    }
    #legend('topright',SPECIES.vec[i],bty='n',cex=1.5)
    mtext(SPECIES.vec[i],4,line=0.5,cex=1.25,las=3)
  }
  mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.5,outer=T)
  if(RELATIVE=="YES") mtext("Relative cpue",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
  if(RELATIVE=="NO") mtext("CPUE (kg/km.gn.day)",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
  dev.off()
  
  
  
  #Export index  
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Index")

    #Annual folly index      
  Foly.series=list(Foly.whis,Foly.gum,Foly.dus,Foly.san)
  for ( i in 1:N.species)
  {
    dat=Foly.series[[i]]
    nm=paste(SPECIES.vec[i],".annual.folly.csv",sep="")
    write.csv(dat,nm,row.names=F)
  }
   
    #Annual base case (zones combined)
  for ( i in 1:N.species)
  {    
    a=Index.BaseCase.CI[[i]]
    dat=data.frame(Finyear=a$FINYEAR,Mean=a$MEAN,
          CV=a$CV,LOW.CI=a$LOW.CI,UP.CI=a$UP.CI,SD=a$SD)
    nm=names(Index.BaseCase.CI)[i]
    nm=paste(SPECIES.vec[i],".annual.abundance.basecase.csv",sep="")
    write.csv(dat,nm,row.names=F)
  }
   
    #Annual base case by zone
  CVs=Index.BaseCase.CI.annual.by.zone
  for ( i in 1:N.species)
  {
    a=Index.BaseCase.CI.annual.by.zone[[i]]
    zns=names(a$MEAN)
    cvs=vector('list',length(zns))
    names(cvs)=zns
    for(z in 1:length(zns))
    {
      dat=data.frame(Finyear=a$FINYEAR,Mean=a$MEAN[[z]],
                     CV=a$CV[[z]],LOW.CI=a$LOW.CI[[z]],UP.CI=a$UP.CI[[z]],SD=a$SD[[z]])
      nm=paste(SPECIES.vec[i],".",zns[z],".annual.abundance.basecase.csv",sep="")
      write.csv(dat,nm,row.names=F)
      cvs[[z]]=a$CV[[z]]
    }
    CVs[[i]]=data.frame(Finyear=a$FINYEAR,CV=do.call(cbind,cvs))
  }

    #Monthly base case by zone
  for ( i in 1:N.species)
  {
    a=Index.BaseCase.CI.by.zone[[i]]
    zns=names(a$MEAN)
    for(z in 1:length(zns))
    {
      dat=data.frame(Finyear=a$FINYEAR,Mean=a$MEAN[[z]],
          CV=a$CV[[z]],LOW.CI=a$LOW.CI[[z]],UP.CI=a$UP.CI[[z]],SD=a$SD[[z]])
      nm=paste(SPECIES.vec[i],".",zns[z],".monthly.abundance.basecase.csv",sep="")
      write.csv(dat,nm,row.names=F)
    }
  }
  
  #Plot CVs
  LT=2:4  
  LWD=3
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper")
  tiff(file="CVs.Annual.by.zone.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1.5,1.5,1,1),oma=c(2.5,2.5,.1,.1),las=1,mgp=c(.1,.6,0))
  for ( i in 1:N.species)
  {
    a=CVs[[i]]
    YR=nrow(a)
    zNs=names(a)[-1]
    if(ncol(a)>2)zNs= substr(names(a)[-1], start=4, stop=10)
    plot(a[,2],type='l',ylim=c(0,max(a[,2:ncol(a)])),lwd=LWD,col=LT[1],ylab='',xaxt='n',xlab='',cex.axis=1.25)
    if(length(zNs)>1)for(z in 2:length(zNs)) lines(a[,z+1],lwd=LWD,lty=1,col=LT[z])

    legend('topright',zNs,bty='n',cex=1.25,lty=1,col=LT,lwd=LWD)
    legend('bottomleft',SPECIES.vec[i],bty='n',cex=1.5)
    axis(1,1:YR,labels=F,tck=-0.015)
    axis(1,seq(1,YR,5),a$Finyear[seq(1,YR,5)],tck=-0.03,cex.axis=1.25)
  }
  mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=2.5,outer=T)
  mtext("CV",side=2,line=0.45,font=1,las=0,cex=2.5,outer=T)
  dev.off()

  

  #Export all bootstraps
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/All.bootstraps")

  for ( i in 1:N.species)
  {
    #Annual
    test=do.call(cbind,Pred.st.BaseCase.CI[[i]])
    a=seq(3,ncol(test),by=2)
    test=test[,-a]  #remove copy of years
    nm=paste(SPECIES.vec[i],".annual.abundance.basecase.csv",sep="")
    write.csv(test,nm,row.names=F)
    NMs=names(test)
    
    #Annual by zone
    fn.out.boot.zn(Pred.st.BaseCase.CI.Annual.by.zone[[i]],".annual.abundance.basecase.csv")            
  
    
    #Monthly by zone
    fn.out.boot.zn(Pred.st.BaseCase.CI.by.zone[[i]],".monthly.abundance.basecase.csv")
  }
})










###NOT USED ####

#9.1. Nominal cpue (folly classic: total catch/total effort)
#scenario 1. Folly (classic annual aggregation of catch and effort)
# fn.folly.cpue=function(DATA,ZONAS)
# {
#   DATA=subset(DATA,!is.na(Km.Gillnet.Days.inv))
#   FINYEAR.monthly=unique(DATA$FINYEAR)
#   N.yrs=length(FINYEAR.monthly)
# 
#   if(ZONAS=="Combined")
#   {
#     Aggregated.catch=aggregate(Catch.Target~FINYEAR,data=DATA,FUN=sum,na.rm=T)
#     Aggregated.effort=aggregate(Km.Gillnet.Days.inv~FINYEAR,data=DATA,FUN=sum,na.rm=T)
#     Aggregated.cpue=data.frame(FINYEAR=Aggregated.catch[,1],cpue=Aggregated.catch[,2]/Aggregated.effort[,2],
#                                effort=Aggregated.effort[,2])
#   }
#   
#   if(ZONAS=="ByZone")
#   {
#     Aggregated.catch=Aggregated.effort=matrix(nrow=N.yrs,ncol=length(ZONES))
#     for ( tt in 1:length(ZONES))
#     {
#       DATA1=subset(DATA,zone==ZONES[tt])
#       if(nrow(DATA1)>0)
#       {
#         Ag.Ktch=aggregate(Catch.Target~FINYEAR,data=DATA1,FUN=sum,na.rm=T)  
#         Ag.Eff=aggregate(Km.Gillnet.Days.c~FINYEAR,data=DATA1,FUN=sum,na.rm=T) 
#         if(!(length(Ag.Ktch$FINYEAR)==N.yrs))
#         {
#           missing=FINYEAR.monthly[which(!(FINYEAR.monthly%in%Ag.Ktch$FINYEAR))]
#           add=Ag.Ktch[1:length(missing),]
#           add[,2]=NA
#           add[,1]=missing
#           Ag.Ktch=rbind(Ag.Ktch,add)
#           Ag.Ktch=Ag.Ktch[order(Ag.Ktch$FINYEAR),]
#           
#           add.eff=add
#           names(add.eff)[2]=names(Ag.Eff)[2]
#           Ag.Eff=rbind(Ag.Eff,add.eff)
#           Ag.Eff=Ag.Eff[order(Ag.Eff$FINYEAR),]
#           
# 
#         }
#         Aggregated.catch[,tt]=Ag.Ktch[,2]
#         Aggregated.effort[,tt]=Ag.Eff[,2]
#       }
#     }
#     Aggregated.cpue=as.data.frame(Aggregated.catch/Aggregated.effort)
#     Aggregated.cpue$FINYEAR=FINYEAR.monthly
#     Aggregated.cpue$effort_West=Aggregated.effort[,1]
#     Aggregated.cpue$effort_Zone1=Aggregated.effort[,2]
#     Aggregated.cpue$effort_Zone2=Aggregated.effort[,3]
#     names(Aggregated.cpue)[1:3]=paste("cpue_",ZONES,sep="")
#   }
#   
#   return(cpue=Aggregated.cpue)
# }
# 
# Folly.cpue.LIVEWT.Combined=Folly.cpue.LIVEWT.byZone=vector('list',length=N.species)
# names(Folly.cpue.LIVEWT.Combined)=names(Folly.cpue.LIVEWT.byZone)=SPECIES.vec
# for ( i in 1:N.species)
# {
#   Folly.cpue.LIVEWT.Combined[[i]]=fn.folly.cpue(subset(DATA.list.LIVEWT.c[[i]],!(VESSEL%in%Ves.no.ktc[[i]])),"Combined")
#   Folly.cpue.LIVEWT.byZone[[i]]=fn.folly.cpue(subset(DATA.list.LIVEWT.c[[i]],!(VESSEL%in%Ves.no.ktc[[i]])),"ByZone")
# }


#Habitat areas
#Manual approach
# AREA.W$"Whiskery shark"=data.frame(BLOCKX=as.character(BLOCKX.w),
#                         Fish.Area=c(13,12,1,16,8.5,2.5,1,12.5,2.5,2.5,14.25,1.5,11,17.5,19.5,25,11,9.5,2,3.5,2.5,4,
#                         24,19.5,6,6,10,8.5,12.5,1.75,6,10.5,15,15,15,14.5,5,1,1,1,2,9,5,5.25,1)/25)
# AREA.W$"Gummy shark"=data.frame(BLOCKX=as.character(BLOCKX.g),
#                       Fish.Area=c(3,1.25,11.5,17.75,19.5,25,25,1.5,2.75,2,3.5,23,16.5,3.25,2.75,1,1,5.75,12.25,13
#                       ,12.5,14.75,14.5,5,1,7.5,8.5,5.75,1)/25)
# AREA.W$"Dusky shark"=data.frame(BLOCKX=as.character(BLOCKX.d),Fish.Area=c(13,11,1,18,8,3,1,12,3,16,
#                       11,9.5,9,13,1.5,6.5,10.5,2,8,8.5,6.5,1)/25)
# AREA.W$"Sandbar shark"=data.frame(BLOCKX=as.character(BLOCKX.s),
#           Fish.Area=c(9,17,5,22,2.5,13,11.5,1,18.5,8,3,1,13,2,15,11,10,8.5,13.5,1.5,2,8,8.5)/25)



# #4. Plot CPUE series     
# Leg=c("Folly.West","Folly.Zn1","Folly.Zn2",
#       "Stand.West","Stand.Zn1","Stand.Zn2")
# #COLORES=gray(seq(0.1,0.9,length=length(Leg)))
# COLORES=c(rep("black",3),rep("grey40",3))
# LINEAS=1:length(Leg)
# 
# #By Zone
# tiff(file="Figure X.Cpue.by.zone.comparison.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# par(mfcol=c(2,2),mar=c(1,3.7,1.5,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
# for ( i in 1:N.species)
# {
#   YRs=Folly.cpue.LIVEWT.byZone[[i]]$FINYEAR
#   N.YRs=length(YRs)
#   plot(1:N.YRs,Folly.cpue.LIVEWT.byZone[[i]][,1],ylab="",type='l',xlab="",xaxt='n',col=COLORES[1],
#        lty=LINEAS[1],lwd=1.5,ylim=c(0,max(Folly.cpue.LIVEWT.byZone[[i]][,1:3],na.rm=T)))
#   for(j in 2:3) lines(1:N.YRs,Folly.cpue.LIVEWT.byZone[[i]][,j],col=COLORES[j],lty=LINEAS[j],lwd=1.5)
#   #for(j in 1:3) lines(1:N.yrs,STAND[[i]][,j],col=COLORES[j+3],lty=LINEAS[j+3],lwd=1.5)    #MISSING: standardised cpue
#   axis(1,at=1:N.YRs,labels=F,tck=-0.02)
#   legend('bottomright',SPECIES.vec[i],bty='n',cex=1)
#   axis(1,at=seq(1,N.YRs,by=5),labels=F,tck=-0.04)
#   if(i==1)  legend('topright',Leg,col=COLORES,lty=LINEAS,bty='n',cex=1,lwd=1.5)
#   if (i%in%c(2,4))axis(1,at=seq(1,N.YRs,by=5),labels=YRs[seq(1,N.YRs,by=5)],tck=-0.04,cex.axis=1,padj=-.05)
# }
# mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.15,outer=T)
# mtext("CPUE (kg/km gillnet day)",side=2,line=-1.3,font=1,las=0,cex=1.15,outer=T)
# dev.off()
# 
# 
# #Combined
# tiff(file="Figure X.Cpue.combined.comparison.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# par(mfcol=c(2,2),mar=c(1,3.7,1.5,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
# for ( i in 1:N.species)
# {
#   YRs=Folly.cpue.LIVEWT.byZone[[i]]$FINYEAR
#   N.YRs=length(YRs)
#   
#   plot(1:N.YRs,Folly.cpue.LIVEWT.Combined[[i]]$cpue,ylab="",type='l',xlab="",xaxt='n',col=COLORES[1],
#        lty=LINEAS[1],lwd=1.5)
#   #lines(1:N.yrs,STAND[[i]],col=COLORES[j+3],lty=LINEAS[j+3],lwd=1.5)    #MISSING: standardised cpue
#   axis(1,at=1:N.YRs,labels=F,tck=-0.02)
#   legend('bottomright',SPECIES.vec[i],bty='n',cex=1)
#   axis(1,at=seq(1,N.YRs,by=5),labels=F,tck=-0.04)
#   if(i==1)  legend('topright',c("Folly","Standardised"),col=COLORES,lty=LINEAS,bty='n',cex=1,lwd=1.5)
#   if (i%in%c(2,4))axis(1,at=seq(1,N.YRs,by=5),labels=YRs[seq(1,N.YRs,by=5)],tck=-0.04,cex.axis=1,padj=-.05)
# }
# mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.15,outer=T)
# mtext("CPUE (kg/km gillnet day)",side=2,line=-1.3,font=1,las=0,cex=1.15,outer=T)
# dev.off()
# 
# 
# 
# #Combined effort and cpue for folly
# PLOT.eff.cpue=function(DATA)
# {
#   DATA[,3]=DATA[,3]/1000
#   max1=max(DATA[,2],na.rm=T)*1.05
#   max2=max(DATA[,3],na.rm=T)*1.05
#   FInYEAR=as.character(unique(DATA$FINYEAR))
#   N=length(FInYEAR)
#   plot(1:N,DATA[,2],ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
#        ,cex.axis=1.2,lwd=1.75)
#   axis(1,at=1:N,labels=F,tck=-0.015)
#   axis(1,at=seq(1,N,5),labels=FInYEAR[seq(1,N,5)],tck=-0.0225,cex.axis=1.1)
#   par(new=T)
#   plot(1:N,DATA[,3],col="grey70",type='l',axes=F,ann='F',lwd=2,ylim=c(0,max2))
#   axis(4,at=pretty(DATA[,3]),labels=pretty(DATA[,3]),las=2,cex.axis=1.2)
#   mtext("Year",1,1.25,outer=T,cex=1.35)
#   mtext("CPUE (kg/km gn.d)",2,-1.15,outer=T,las=3,cex=1.35)
#   mtext("Effort ('000 km gn.d)",4,-1.5,outer=T,las=3,cex=1.35)
# }
# 
# tiff(file="FigureX.Combined.eff.cpue.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# 
# par(mfcol=c(2,2),mar=c(1,3.5,1.5,3.5),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
# #par(mfcol=c(2,2),mar=c(1,3.7,1.5,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
# for ( i in 1:N.species)PLOT.eff.cpue(Folly.cpue.LIVEWT.Combined[[i]])
# dev.off()

# do.boot="YES"
# #do.boot="NO"
# 
# #creating boostrapped datasets
# Bootstrapped.Dat=vector('list',length=N.species)
# names(Bootstrapped.Dat)=SPECIES.vec
# Bootstrapped.Dat.good=Bootstrapped.Dat
# 
# if(Boot.method=="DATA")
# {
#   if(do.boot=="YES")
#   {
#     #note: 1000 bootstraps take 3 hours.Run only one and then reload data
#     system.time(for (i in 1:N.species)
#     {
#       a=b=vector('list',N.boot)
#       for ( n in 1:N.boot)
#       {
#         a[[n]]=Boot.fn(DATA.list.LIVEWT.c[[i]])
#         b[[n]]=Boot.fn(DATA.list.LIVEWT.goodRep[[i]])
#       }      
#       Bootstrapped.Dat[[i]]=a
#       Bootstrapped.Dat.good[[i]]=b
#     })
#     setwd("C:/Matias/Analyses/Catch and effort/Outputs/BootData")
#     save(Bootstrapped.Dat,file="Boot.RData")
#     save(Bootstrapped.Dat.good,file="Boot_good.RData")
#     
#   }
#   load("Boot.RData")
#   load("Boot_good.RData")
# }
# 
# 
# if(Boot.method=="RESIDUALS")
# {
#   #note: this approach is fast
#   system.time(for (i in 1:N.species)
#   {
#     a=b=vector('list',N.boot)
#     for ( n in 1:N.boot)
#     {
#       a[[n]]=Boot.fn(Stand.BaseCase[[i]],DATA.list.LIVEWT.c[[i]])
#       b[[n]]=Boot.fn(Stand.GoodRep[[i]],DATA.list.LIVEWT.goodRep[[i]])
#       
#     }
#     Bootstrapped.Dat[[i]]=a
#     Bootstrapped.Dat.good[[i]]=b
#   })
#   
# }
###

#9.2.2. Define new vessel levels
#Create new vessel and block levels               #NOT USED
#note: this fits a mixed model to catch rate and vessel to regroup vessels in 5 categories
# New.Ves.Lev=vector('list',N.species)
# names(New.Ves.Lev)=names(DATA.list.LIVEWT.c)
# STORE.a=New.Ves.Lev
# for(i in 1:N.species)
# {
#   n.zn=length(Zone.list[[i]])
#   VES.coef=vector('list',n.zn)
#   names(VES.coef)=Zone.list[[i]]
#   store.a=VES.coef
#   for(j in 1:n.zn)
#   {
#     #note: this drops the vessels that never reported catch of target species
#     DATA=subset(DATA.list.LIVEWT.c[[i]],zone==Zone.list[[i]][[j]]&!(VESSEL%in%Ves.no.ktc[[i]]))
#     
#     #convert to factors
#     DATA$VESSEL=as.factor(DATA$VESSEL)
#     DATA$BLOCKX=as.factor(DATA$BLOCKX)
#     DATA$MONTH=as.factor(DATA$MONTH)
#     DATA$FINYEAR=as.factor(DATA$FINYEAR)
# 
#     #fit simple mixed model for ranking all vessels 
#     simple.mixed=lme((Catch.Target/Km.Gillnet.Days.c)~1, data = DATA, random = ~ 1 |VESSEL)
#     vessel.coef=data.frame(VESSEL=rownames(ranef(simple.mixed)),coef=ranef(simple.mixed))
#     percentiles=quantile(vessel.coef$X.Intercept.,probs=seq(0,1,.20))
#     
#     a=aggregate((Catch.Target/Km.Gillnet.Days.c)~VESSEL,DATA,mean)
#     b=aggregate((Catch.Target/Km.Gillnet.Days.c)~VESSEL,DATA,sd)
#     a=merge(a,b,by="VESSEL")
#     names(a)[2:3]=c("Mean","SD")
#     a=merge(a,vessel.coef,by="VESSEL",all.x=T)
#     
#     #3. create new vessel levels
#     a$new.level.vessel=with(a,
#             ifelse(X.Intercept.<percentiles[2]&Mean>0,"VL",
#             ifelse(X.Intercept.>=percentiles[2]&X.Intercept.<percentiles[3]&Mean>0,"L",
#             ifelse(X.Intercept.>=percentiles[3]&X.Intercept.<percentiles[4]&Mean>0,"A",
#             ifelse(X.Intercept.>=percentiles[4]&X.Intercept.<percentiles[5]&Mean>0,"H",
#             ifelse(X.Intercept.>=percentiles[5]&Mean>0,"VH",
#             ifelse(Mean==0,"VL",NA)))))))
#     
#     #relevel factors
#     #a$new.level.vessel=relevel(a$new.level.vessel, "VL")#set to specific level
#     a$new.level.vessel <- factor(a$new.level.vessel, levels = c("VL","L","A","H","VH"))#reorder whole thing
#     
#     
#     vessel.coef=a[,c(1,5)]
#     vessel.coef$new.level.vessel=as.factor(vessel.coef$new.level.vessel)
#     VES.coef[[j]]=vessel.coef
#     store.a[[j]]=a
#     
#   }
#   New.Ves.Lev[[i]]=VES.coef
#     STORE.a[[i]]=store.a
# }

# #Extract new levels for blocks
# Block.new.levels.Top=Block.new.levels.Group=vector('list',length=N.species)
# 
# for (i in 1:N.species)
# {
#   Block.new.levels.Top[[i]]=Data.Summary[[i]]$TopBlocks.records
#   Block.new.levels.Group[[i]]=Data.Summary[[i]]$agg
# }


#Create new vessel and block levels grouping few-record levels
#store.a=vector('list',length=N.species)  #store mean catches and reml rankings
# for (i in 1:N.species)
#   {
#   
#   DATA.list.LIVEWT.c[[i]]$BLOCKX.new=ifelse(as.character(DATA.list.LIVEWT.c[[i]]$BLOCKX)%in%
#           as.character(Block.new.levels.Top[[i]]$BLOCKX),as.character(DATA.list.LIVEWT.c[[i]]$BLOCKX),
#     ifelse(as.character(DATA.list.LIVEWT.c[[i]]$BLOCKX)%in%as.character(Block.new.levels.Group[[i]]$BLOCKX) &
#       DATA.list.LIVEWT.c[[i]]$zone=="West","PlusGroup.West",
#     ifelse(as.character(DATA.list.LIVEWT.c[[i]]$BLOCKX)%in%as.character(Block.new.levels.Group[[i]]$BLOCKX) &
#       DATA.list.LIVEWT.c[[i]]$zone=="Zone1","PlusGroup.Zone1",
#     ifelse(as.character(DATA.list.LIVEWT.c[[i]]$BLOCKX)%in%as.character(Block.new.levels.Group[[i]]$BLOCKX) &
#       DATA.list.LIVEWT.c[[i]]$zone=="Zone2","PlusGroup.Zone2",NA))))
#   
#   
#   DATA.list.LIVEWT.c[[i]]$BLOCKX.new=as.factor(DATA.list.LIVEWT.c[[i]]$BLOCKX.new)
#   
#   #1. manipulate data
#   DATA.list.LIVEWT.c[[i]]$VESSEL=DATA.list.LIVEWT.c[[i]]$VESSEL[, drop=TRUE]#remove vessel levels that don't occur in the dataset
#  DATA.list.LIVEWT.c[[i]]$CPUE=DATA.list.LIVEWT.c[[i]]$Catch.Target/DATA.list.LIVEWT.c[[i]]$Km.Gillnet.Days.c
#   DATA.list.LIVEWT.c[[i]]$FINYEAR=as.factor(DATA.list.LIVEWT.c[[i]]$FINYEAR)
#   DATA.list.LIVEWT.c[[i]]$MONTH=as.factor(DATA.list.LIVEWT.c[[i]]$MONTH)
#   DATA.list.LIVEWT.c[[i]]$BLOCKX=as.factor(DATA.list.LIVEWT.c[[i]]$BLOCKX)
#   DATA.list.LIVEWT.c[[i]]$BLOCKX.new=as.factor(DATA.list.LIVEWT.c[[i]]$BLOCKX.new)
#   
#   #2. fit simple mixed model for ranking all vessels 
#   PosDat=subset(DATA.list.LIVEWT.c[[i]],Catch.Target>0)
#   #simple.mixed=lme(CPUE ~ FINYEAR + BLOCKX.new +MONTH, data = PosDat, random = ~ 1 |VESSEL)
#   simple.mixed=lme(Catch.Target~FINYEAR+BLOCKX.new+MONTH, data = PosDat, random = ~ 1 |VESSEL)
#   vessel.coef=data.frame(VESSEL=rownames(ranef(simple.mixed)),coef=ranef(simple.mixed))
#   percentiles=quantile(vessel.coef$X.Intercept.,probs=seq(0,1,.20))
#   
# #    a=aggregate(CPUE~VESSEL,DATA.list.LIVEWT.c[[i]],mean)
# #    b=aggregate(CPUE~VESSEL,DATA.list.LIVEWT.c[[i]],sd)
#   a=aggregate(Catch.Target~VESSEL,DATA.list.LIVEWT.c[[i]],mean)
#   b=aggregate(Catch.Target~VESSEL,DATA.list.LIVEWT.c[[i]],sd)
#   a=merge(a,b,by="VESSEL")
#   names(a)[2:3]=c("Mean","SD")
#   a=merge(a,vessel.coef,by="VESSEL",all.x=T)
#   
#   #3. create new vessel levels
#   a$new.level.vessel=with(a,ifelse(X.Intercept.<percentiles[2]&Mean>0,"VL",
#                   ifelse(X.Intercept.>=percentiles[2]&X.Intercept.<percentiles[3]&Mean>0,"L",
#                   ifelse(X.Intercept.>=percentiles[3]&X.Intercept.<percentiles[4]&Mean>0,"A",
#                   ifelse(X.Intercept.>=percentiles[4]&X.Intercept.<percentiles[5]&Mean>0,"H",
#                   ifelse(X.Intercept.>=percentiles[5]&Mean>0,"VH",
#                   ifelse(Mean==0,"VL",NA)))))))
#   
#       #relevel factors
#   #a$new.level.vessel=relevel(a$new.level.vessel, "VL")#set to specific level
#   a$new.level.vessel <- factor(a$new.level.vessel, levels = c("VL","L","A","H","VH"))#reorder whole thing
#   
#   
#   vessel.coef=a[,c(1,5)]
#   vessel.coef$new.level.vessel=as.factor(vessel.coef$new.level.vessel)
#   DATA.list.LIVEWT.c[[i]]=merge(DATA.list.LIVEWT.c[[i]],vessel.coef,by="VESSEL")
#   store.a[[i]]=a
#   }

# Interviewed.vess=c("E 067","F 517","F 505","B 067","B 142") #vessels interviewed for white shark catch recons.
# Interviewed.vess.rank=vector('list',length=N.species)
# names(Interviewed.vess.rank)=SPECIES.vec
# for (i in 1:N.species)
#   {
#     a=subset(DATA.list.LIVEWT.c[[i]],VESSEL%in%Interviewed.vess)
#     a=a[-which(duplicated(a$VESSEL)),]
#     a=a[,match(c("VESSEL","new.level.vessel"),names(a))]
#     Interviewed.vess.rank[[i]]=a
#   }



#BEST ERROR STRUCTURE
#3.2.4. Define best error structure
#First test without interaction for speeding process and without zone (issues with zero-trunc optim)
# form1=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL.new)
# form2=as.formula(Catch.Target~FINYEAR+BLOCKX+VESSEL.new+MONTH)
# form3=as.formula(log(Catch.Target)~FINYEAR+BLOCKX+VESSEL.new+MONTH)
# 
# Best.bin=list(form1,form1,form1,form2)
# Best.pos=list(form2,form2,form2,form2)
# Best.pos.log=list(form3,form3,form3,form3)
# 
# 
# VARIABLES=c("Catch.Target","FINYEAR","BLOCKX","zone","VESSEL.new","MONTH")
# 
 

#3.2.4.1 Compare models fit
# fn.best.error=function(DATA,FORMULA.pos,FORMULA.pos.log)
# {
#   # 1. Put data in proper shape
#   DataFile=DATA[,match(VARIABLES,names(DATA))]
#   
#   #drop levels not occurring in data
#   for(f in 1:ncol(DataFile))
#   {
#     if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
#   }
#   
#   
#   # Create binary dataset (0/1 for CPUE)
#   BiData <- DataFile
#   BiData[,1] <- as.numeric(DataFile[,1]>0)
#   
#   # Create positive (non-zero) dataset
#   PosData <- DataFile[DataFile[,1]>0,]
#   
#   
#   #round catch for discrete distributions
#   DataFile$Catch.Target=with(DataFile,ifelse(Catch.Target>0 &Catch.Target<1,1,round(Catch.Target)))
#   
#   
#   #Fit different error structures
#   
#   #Binomial
#   Binomial <- glm(FORMULA.pos, data=BiData, family="binomial", maxit=100)
#   
#   #Lognormal
#   GLMlog <- glm(FORMULA.pos.log, data=PosData, family=gaussian, maxit=100)
#   
#   #Gamma
#   GLMgam<-glm(FORMULA.pos, data=PosData, family=Gamma(link=log), maxit=100)
#   
#   #Zero-truncated Poisson
#   #GLM.hurdle.Pois <- hurdle(FORMULA.pos, data = DataFile,dist = "poisson",zero.dist = "binomial")
#   Zero.t.Pois <- glmmadmb(FORMULA.pos,data=subset(DataFile,Catch.Target>0),family="truncpoiss",
#                           admb.opts=admbControl(maxfn=1000,imaxfn=1000))
#   
#   #Zero-truncated Negative Binomial
#   #GLM.hurdle.neg.bin <- hurdle(FORMULA.pos,data = DataFile,dist = "negbin",zero.dist = "binomial")
#   Zero.t.NB <- glmmadmb(FORMULA.pos,data=subset(DataFile,Catch.Target>0),family="truncnbinom1")
#   
#   
#   return(list(GLMlog=GLMlog,GLMgam=GLMgam,Zero.t.Pois=Zero.t.Pois,Zero.t.NB=Zero.t.NB))
# }
# 
# 
# Best.error.catch=Best.error.TABLE=vector('list',length=N.species)
# names(Best.error.catch)=names(Best.error.TABLE)=SPECIES.vec
# 
# system.time(for ( i in 1:N.species) 
# {
#   Best.error.catch[[i]]=fn.best.error(DATA.list.LIVEWT.c[[i]],Best.pos[[i]],Best.pos.log[[i]])
#   
# })

# MOD=names(Best.error.catch[[1]])
#plot model fit

# for ( i in 1:N.species) 
# {
#   for (j in 1:length(MOD))
#   {
#     tiff(file=paste("C:/Matias/Analyses/Catch and effort/Outputs/Best.error/Best.error.",MOD[j],SPECIES.vec[i],".fit.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
#     fn.plot.diag(Best.error.catch[[i]][[j]],MOD[j],SPECIES.vec[i])
#     dev.off()
#   }
#   MODEL.LIST=Best.error.catch[[i]]
#   Best.error.TABLE[[i]]=Table.fit(MODEL.LIST,"NO")
# }  
# 
# for(i in 1:N.species) 
# {
#   write.table(Best.error.TABLE[[i]],
#               paste("C:/Matias/Analyses/Catch and effort/Outputs/Best.error/Best.error.table",SPECIES.vec[i],".csv",sep=""),sep = ",")
# }








# ####
# #STATE OF FISHERIES REPORT
# ####
# 
# setwd("C:/Matias/Analyses/Catch and effort/State of fisheries")
# fun.Fig.5.8.SoFar=function(DAT,DAT1,scaler)
# {
#   FInYEAR=as.character(unique(DAT1$FINYEAR))
#   NN=length(FInYEAR)
#   
#   par(mar=c(4,4.7,2,5),mgp=c(2.5,.65,0),las=1)
#   plot(1:NN,DAT,type="o",col=1,axes=F,ann='F',lwd=2,cex.lab=1.3,pch=21)
#   axis(2,at=pretty(DAT),labels=pretty(DAT),col=1,col.ticks=1,col.axis=1)
#   par(new=T)  
#   plot(1:NN,DAT1[,2]/scaler,type='l',col="grey80",xaxt='n',yaxt='n',ann='F',
#        lwd=2,cex.lab=1.3)
#   axis(4,at=pretty(DAT1[,2]/scaler),labels=pretty(DAT1[,2]/scaler),col=1,col.ticks=1,col.axis=1)
# 
#   axis(1,at=1:NN,labels=F,tck=-0.0125)
#   axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
#   
#   mtext("CPUE (kg/km gn.d)",side=2,line=3.25,las=3,cex=1.5)
#   mtext("Effort (1000 km gn.d)",side=4,line=3.25,las=3,cex=1.5)
#   mtext("Financial year",side=1,line=2,cex=1.5)
# }
# 
# 
# for (i in 1:N.species)
# {
#   annual.effort.days.total=aggregate(Km.Gillnet.Days.c~FINYEAR,data=DATA.list.LIVEWT.c[[i]],sum,na.rm=T)
#   CPUE=Stand.cpue.imp[[i]][,2]
# 
#   jpeg(file=paste("Figure.",i+4,SPECIES.vec[i],".Eff.CPUE.jpeg",sep=""),width = 2400, height = 2400,units = "px", res = 300)
#   fun.Fig.5.8.SoFar(CPUE,annual.effort.days.total,1000)
#   dev.off()
# }
# 




##############





# tiff(file="Figure 4.Data subsetting.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
# par(mfrow=c(2,1),mar=c(3,.1,.1,.1),oma=c(1,4,.1,.1),las=1,mgp=c(1,.6,0))
# 
# a=vector(length=4);for(i in 1:4)a[i]=paste(SPECIES.vec[i]," (",Data.Summary[[i]]$Per.rec.kept.vess,"% records kept)",sep="")
# b=vector(length=4);for(i in 1:4)b[i]=paste(SPECIES.vec[i]," (",Data.Summary[[i]]$Per.rec.kept.block,"% records kept)",sep="")
# 
# line.col=c("black","black","grey65","grey65")
# line.type=1:4
# plot(Data.Summary[[1]]$Ves.Cum.Ca,ylab="", xlab="",xaxt='n',type='l',cex.lab=1.25,lwd=3)
# for (i in 2:N.species)lines(Data.Summary[[i]]$Ves.Cum.Ca,col=line.col[i],lwd=3,lty=line.type[i])
# axis(1,seq(0,length(Data.Summary[[i]]$Ves.Cum.Ca),10),labels=F,tck=-0.015)
# axis(1,seq(100,length(Data.Summary[[i]]$Ves.Cum.Ca),100),labels=seq(100,length(Data.Summary[[i]]$Ves.Cum.Ca),100),tck=-0.03,cex.axis=1)
# #legend("bottomright",a,bty='n',cex=.75,col=line.col,lwd=3,lty=line.type)
# #abline(CatchThreshold,0,col=2,lwd=3)
# mtext("Number of Vessels",side=1,line=1.6,font=1,las=0,cex=1.5,outer=F)
# 
# 
# plot(Data.Summary[[1]]$Block.Cum.Ca,ylab="", xlab="",xaxt='n',type='l',cex.lab=1.25,lwd=3)
# for (i in 2:N.species)lines(Data.Summary[[i]]$Block.Cum.Ca,col=line.col[i],lwd=3,lty=line.type[i])
# axis(1,1:length(Data.Summary[[i]]$Block.Cum.Ca),labels=F,tck=-0.015)
# axis(1,seq(5,length(Data.Summary[[i]]$Block.Cum.Ca),5),labels=seq(5,length(Data.Summary[[i]]$Block.Cum.Ca),5),tck=-0.03,cex.axis=1)
# #abline(CatchThreshold,0,col=2,lwd=2)
# mtext("Number of blocks",side=1,line=1.6,font=1,las=0,cex=1.5,outer=F)
# mtext("Cumulative catch (%)",side=2,line=2,font=1,las=0,cex=1.5,outer=T)
# #legend("bottomright",b,bty='n',cex=.75,col=line.col,lwd=3,lty=line.type)
# legend("topleft",SPECIES.vec,bty='n',cex=.95,col=line.col,lwd=3,lty=line.type)
# 
# pos.list=list(c(.27,.62,.72,.963),c(.27,.62,0.545,.788),c(.63,.98,.72,.963),c(.63,.98,0.545,.788))
# #pos.list=list(c(.27,.62,.72,.963),c(.27,.62,.57,.813),c(.63,.98,.72,.963),c(.63,.98,.57,.813))
# for ( i in 1:N.species)
# {
#   par(fig=pos.list[[i]], new = T,mgp=c(.1,.2,0))
#   
#   plot(Data.Summary[[i]]$TABLE15$Number.obs,Data.Summary[[i]]$TABLE15$Catch.Target/1000,col=Data.Summary[[i]]$TABLE15$Top,pch=19,
#        ylab="",xlab="",xaxt='n',yaxt='n',ylim=c(0,1010),xlim=c(0,720))
#   if(i==4)legend("topleft",c("dropped","selected"),col=unique(Data.Summary[[i]]$TABLE15$Top),pch=19,bty='n',cex=.8)
#   legend("bottomright",SPECIES.vec[i],bty='n',cex=.8)
#   axis(1,at=seq(0,800,200),labels=F,tck=-0.03)
#   axis(2,at=seq(0,1000,250),labels=F,tck=-0.03)
#   if(i%in%c(2,4))axis(1,at=seq(0,800,200),labels=seq(0,800,200),cex.axis=.7,tck=-0.03,padj=-0.5)
#   if(i%in%c(1,2))axis(2,at=seq(0,1000,250),labels=seq(0,1000,250),cex.axis=.7,tck=-0.03)
# }  
# mtext("                                                                                                 Total catch per vessel (tonnes)",
#       side=2,line=-8,font=1,las=0,cex=1,outer=T)
# mtext("                                             Number of records",
#       side=1,line=-23.5,font=1,las=0,cex=1,outer=T)
# 
# 
# pos.list=list(c(.27,.62,.22,.463),c(.27,.62,0.045,0.288),c(.63,.98,.22,.463),c(.63,.98,0.045,0.288))
# #pos.list=list(c(.27,.62,.22,.463),c(.27,.62,.07,0.313),c(.63,.98,.22,.463),c(.63,.98,.07,0.313))
# for ( i in 1:N.species)
# {
#   par(fig=pos.list[[i]], new = T,mgp=c(.1,.2,0))
#   
#   plot(Data.Summary[[i]]$TABLE14$Number.obs,Data.Summary[[i]]$TABLE14$Catch.Target/1000,col=Data.Summary[[i]]$TABLE14$Top,pch=19,
#        ylab="",xlab="",xaxt='n',yaxt='n',ylim=c(0,1010),xlim=c(0,720))
#   axis(1,at=seq(0,800,200),labels=F,tck=-0.03)
#   axis(2,at=seq(0,1000,250),labels=F,tck=-0.03)
#   if(i%in%c(2,4))axis(1,at=seq(0,800,200),labels=seq(0,800,200),cex.axis=.7,tck=-0.03,padj=-0.5)
#   if(i%in%c(1,2))axis(2,at=seq(0,1000,250),labels=seq(0,1000,250),cex.axis=.7,tck=-0.03)
#   if(i==2)mtext("                       Total catch per block (tonnes)",side=2,line=1.65,font=1,las=0,cex=1,outer=F)
# }  
# 
# mtext("                                             Number of records",
#       side=1,line=-4.1,font=1,las=0,cex=1,outer=T)
# dev.off()





# tiff(file=paste("Figure 10.Blocks per year-month fished.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")  
# par(mfrow=c(2,2),mar=c(4,4,1,1),las=1)
# for ( i in 1:N.species)
# {
#   hist(Data.Summary[[i]]$Block.Yr.Mn,breaks=100,xlim=c(0,432),main="",xlab="Number of year-months fished per block",cex.lab=1)
#   legend("topright",SPECIES.vec[i],bty='n',cex=1.25)
#   box()
#   }
# dev.off()

###############

# #OLD ONES
# Boot.fn=function(DAT,VES.No.K,vars)  
# {
#   #1. drop vessels that never reported catch of target species   
#   DAT=subset(DAT,!(VESSEL%in%VES.No.K))
#   
#   #2. put data in proper shape
#   DataFile=DAT[,match(vars,names(DAT))]
#   
#   DataFile$YrBlk=with(DataFile,paste(FINYEAR,BLOCKX))
#   YrBLK=unique(DataFile$YrBlk)
#   
#   #3. resample with replacement
#   Lista=vector('list',length(YrBLK))
#   for (x in 1:length(YrBLK))
#   {
#     dat=subset(DataFile,YrBlk==YrBLK[x])
#     ID=sample(1:nrow(dat),replace=T)
#     dat=dat[ID,]
#     Lista[[x]]=dat
#   }
#   
#   New.dat=do.call(rbind,Lista)
#   return(New.dat=New.dat)
# }
# fn.bootstrap.index=function(DAT,Spc,vrs,best,Grid.d,Ves.no,area.w,zonas,bloques,ZNEs,Min.w) #function for running show
# {
#   #1. Set storing objects
#   Store.boot.BaseCase=vector('list',N.boot)
#   
#   #2. Run bootstrap
#   for (n in 1:N.boot) 
#   {
#     tryCatch(
# {
#   #2.1 create bootstrapped data
#   if(Boot.method=="DATA") Boot.Dat=Boot.fn(DAT,Ves.no,vrs)
#   if(Boot.method=="RESIDUALS") Boot.Dat=Bootstrapped.Dat[[i]][[n]]
#   
#   #2.2. refit GLMs to bootstrapped data
#   Fit.BaseCase=fn.stand.cpue(Boot.Dat,Ves.no,Spc,best$Bi,
#                              best$Log,vrs,Eff.Creep="YES")
#   
#   #2.3. run predictions
#   Pred.BaseCase=predict.ktch(Fit.BaseCase,Grid.d,Spc)
#   
#   
#   #2.4. run imputations      
#   if(ZNEs=="Zones.combined") Pred.BaseCase$Pos.dat=impute.ktch(Pred.BaseCase$Pos.dat)
#   if(ZNEs=="By.zone") Pred.BaseCase$Pos.dat=impute.ktch.month(Pred.BaseCase$Pos.dat)
#   
#   #2.5. construct abundance indices for each zone or combined
#   if(ZNEs=="Zones.combined") Store.Zone.boot=construct.index(Pred.BaseCase,area.w,"YES",Spc,Weight.by.area="YES")
#   
#   if(ZNEs=="By.zone")
#   {
#     ZN=zonas
#     n.zone=length(ZN)
#     BLK.ZN=bloques
#     
#     Store.Zone.boot=vector('list',n.zone)
#     names(Store.Zone.boot)=ZN
#     for (z in 1:n.zone)
#     {
#       BZ=BLK.ZN[[z]]
#       Pred.zn=vector('list',length(Pred.BaseCase))
#       names(Pred.zn)=names(Pred.BaseCase)
#       for(p in 1:length(Pred.BaseCase))
#       {
#         if(length(Pred.BaseCase[[p]])==0)Pred.zn[p]=Pred.BaseCase[p]
#         if(length(Pred.BaseCase[[p]])>0)Pred.zn[[p]]=subset(Pred.BaseCase[[p]],BLOCKX%in%BZ)
#       }
#       
#       Store.Zone.boot[[z]]=construct.index.month(Pred.zn,subset(area.w,BLOCKX%in%BZ),"YES",Spc,Weight.by.area="YES")
#     }
#   }
#   
#   
#   #2.6. store bootstrapped index
#   Store.boot.BaseCase[[n]]=Store.Zone.boot
#   
#   
#   #2.7. free up memory
#   rm(Boot.Dat,Fit.BaseCase,Pred.BaseCase)
#   
# }, error=function(e){})
#   }
# 
# return(Store.boot.BaseCase)
# }
# fn.get.mean=function(DAT,what)  #function for extracting mean, CV and confidence intervals     
# {
#   if(what=="ZonesCombined")
#   {
#     YRs=DAT[[1]]$FINYEAR
#     N.YRs=length(YRs)
#     N.rep=length(DAT)
#     BC=matrix(NA,nrow=N.YRs,ncol=N.rep)    
#     for(o in 1:N.rep) BC[,o]=DAT[[o]]$I_y.b.W
#     
#     #mean
#     MEAN=rowMeans(BC)
#     
#     #SD
#     SD=apply(BC,1,sd)
#     
#     #CV   
#     CV=SD/MED
#     
#     #CI
#     LOW=apply(BC, 1, function(x) quantile(x, 0.025))
#     UP=apply(BC, 1, function(x) quantile(x, 0.975))
#     
#     return(list(MEAN=MEAN,LOW.CI=LOW,UP.CI=UP,FINYEAR=YRs,CV=CV))       
#     
#   }
#   if(what=="ByZones")
#   {
#     YR.Mns=DAT[[1]][[1]]$FINYEAR
#     N.YRs=length(YR.Mns)
#     N.rep=length(DAT)
#     zns=names(DAT[[1]])
#     MEAN=vector('list',length(zns))    
#     names(MEAN)=zns
#     LOW=UP=SD=CV=MEAN
#     for(z in 1:length(zns)) 
#     {
#       BC=matrix(NA,nrow=N.YRs,ncol=N.rep)
#       for(o in 1:N.rep)
#       {
#         dummy=DAT[[o]][[z]]
#         dummy=dummy[order(dummy$FINYEAR,dummy$MONTH),] #order by year-month
#         BC[,o]=dummy$I_y.b.W
#       }
#       
#       MEAN[[z]]=rowMeans(BC) #mean     
#       
#       #SD
#       SD[[z]]=apply(BC,1,sd)
#       
#       #CV   
#       CV[[z]]=SD[[z]]/MEAN[[z]]
#       
#       LOW[[z]]=apply(BC, 1, function(x) quantile(x, 0.025)) #CI
#       UP[[z]]=apply(BC, 1, function(x) quantile(x, 0.975))
#     }
#     return(list(MEAN=MEAN,LOW.CI=LOW,UP.CI=UP,SD=SD,CV=CV,FINYEAR=YR.Mns))    
#   }
# }
# 
