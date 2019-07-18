###--- SHARK GILLNET AND LONGLINE FISHERY CATCH AND EFFORT MANIPULATIONS ----###

#MISSING:  1. charter.mean.w is dummy. Go out with them to calculate properly!
#             another issue with Charter data is species id (e.g. bull sharks for licence FT1L373 may not be)
          #2. For NSF add GN equivalent effort to LL (ask Rory for rule)

#Index for navigation
#SECTION A. ---- MONTHLY RECORDS 
#SECTION B. ---- DAILY LOGBOOKS 
#SECTION C. ---- CATCH MERGING AND CORRECTIONS 
#SECTION D. ---- EFFORT INSPECTIONS 
#SECTION E. ---- EFFORT CORRECTIONS 
#SECTION F. ---- PROCEDURE SECTION 
          # F 1. EXTRACT QUANTITIES 
          # F 2. EXPORT DATA FOR ASSESSMENT AND CPUE STANDARDISATION
#SECTION G. ---- REPORT SECTION 
          # G 1. SPATIO TEMPORAL CATCH AND EFFORT DISTRIBUTION
          # G 2. DATA REQUESTS 
#SECTION H. ----  EXPORT TOTAL CATCH FOR REFERENCE POINT ANALYSIS ---
#SECTION I. ----  EXPLORATORY ANALYSES ---   (including movies)
#SECTION J. ----  DROPPED CODE ---


#NOTES:

#MONTHLY RECORDS:   CAESS data from 1988/89 to 2007/08 (some fishers kept reporting in CAESS 
#                      after introduction of daily logbooks).
#                   Records include shark data from all fisheries, not just TDGDLF (fishery
#                     code: SGL & WCGL) and NSF (code C127 & CL02)
#                   Records include scalefish data from TDGDLF and NSF only
#                   Hence, combine Rory's Table 81d.mdb 1975/76 to 1987/88 and CAESS for full series

#                   In 2002 each unit became 270 m of net (Rory per comm)
#                   After 1988, estuary blocks are true shark shots (Rory per comm)

#FishCUBE index
               #Extract data for FishCUBE
               #G 4.9 FishCUBE 



#Data inspections:
                 # CATCH INSPECTIONS
                 #     Set current.yr; Identify 0 catch shots (are they real?), ID catch records
                 #         where catch is too large of too small compared to species weight range

                 # EFFORT INSPECTIONS
                 #     Visually inspect effort variables, identify suspicious records and raise to data
                 #         entry staff


#Data structure
#     . Data.monthly: CAESS monthly records 1975-1976 to 2005-2006
#               each row is a species' monthly catch per fisher per block, 
#               including the effort exerted by that fisher in that block

#     . Data.daily: daily logbooks 2006-2007 onwards


#Catch correction criteria:
#           1. Reapportion catch composition of main commercial shark species
#           2. All catch data prior to 1990 are increased by 5%


#Effort correction criteria:
#           1. Fix invalid effort records by
#                 - using vessel's annual average parameter values if available
#                 - using the monthly fleet average if average annual values not available
#                   for the individual vessel (e.g. vessel only submitting 1 record for that year)
#           2. All effort data prior to 1990 are increased by 5%


#Effort standardisation:
#           1. Standardise longline effort to gillnet equivalent effort (used for reporting only)


#Effective effort:
#           1. Define area ranges for each species and use only those records for CPUE analysis


#Fishing efficiency:
#           1. Adjust for increase in fishing efficiency for early years

#Estuaries catch:
#           1. Catch from estuaries is kept for total catch calculation but dropped for cpue standardisation


##############--- DATA SECTION ---###################

rm(list=ls(all=TRUE))
library(RODBC)				#include ODBC library for importing Acccess data
library(PBSmapping)    	#needed to obtain maps
library(plotrix)  			#needed for graph legends
library(adehabitatHR)   #for kernel utilization distribution
library(lubridate)
library(zoo)
library(data.table)   #for fast csv exporting  
library(Rcpp)
library(ggplot2)
require(animation) # NB, must install ImageMagick
library(dplyr)
library(tidyr)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
par.default=par()

source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/SoFaR.figs.R")
#source("C:/Matias/Analyses/SOURCE_SCRIPTS/send.emails.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R")
fn.scale=function(x,scaler) ((x/max(x,na.rm=T))^0.5)*scaler

setwd("C:/Matias/Data/Catch and Effort")  # working directory



First.run="NO"    #turn to yes when new year's data is included
#First.run="YES"

Current.yr="2017-18"    #Set current financial year 
Current.yr.dat=paste(substr(Current.yr,1,4),substr(Current.yr,6,7),sep="_")



#Monthly records 

  #CAESS data 
#note: CAESS is dynamic, constantly updated, Table81.d is static so only use Table 81.d for 
#       years prior to 1988-89 as CAESS data starts in 1988-89.
#       Also, don't use CAESS from non TDGDLF of NSF fisheries for catch rate standardisations as
#       0 catch records are missing/incomplete

#Get.CAESS.Logbook="NO"
Get.CAESS.Logbook="YES"

#if(Get.CAESS.Logbook=="YES")Data.monthly.CAESS=read.csv("CAESS/CAES.data.8889.to.1415.csv",stringsAsFactors =F) 		
if(Get.CAESS.Logbook=="YES")Data.monthly.CAESS=read.csv("M:/CAEMaster/Commercial/FishCubeWA/Fisheries/CAES Monthly Fisheries/02. Data/CAES Monthly Data - Shark.csv",stringsAsFactors =F) 		
if(Get.CAESS.Logbook=="YES") Mesh.monthly=subset(Data.monthly.CAESS,METHOD=="GN",select=c(VESSEL,FINYEAR,MONTH,BLOCKX,METHOD,MSLOW,MSHIGH))

  #Table81.d data
channel <- odbcConnectAccess2007("Table 81d.mdb")      
Data.monthly=sqlFetch(channel, "CAESS data", colnames = F) 		#selecting whole table
close(channel)
drop=c("VesselID","BlockAveID","AnnualVesselAveID","SpeciesID","ReturnID")
Data.monthly=Data.monthly[,-match(drop,names(Data.monthly))]

  #combine Table 81d.mdb and CAESS
if(Get.CAESS.Logbook=="YES")
{
  Data.monthly$FINYEAR=as.character(Data.monthly$FINYEAR)
  Data.monthly$VESSEL=as.character(Data.monthly$VESSEL)
  Data.monthly$METHOD=as.character(Data.monthly$METHOD)
  Data.monthly$SNAME=as.character(Data.monthly$SNAME)
  Data.monthly$CONDITN=as.character(Data.monthly$CONDITN)
  keep=c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81","1981-82",
          "1982-83","1983-84","1984-85","1985-86","1986-87","1987-88")
  Data.monthly=subset(Data.monthly,FINYEAR%in%keep)
  Data.monthly$fishery=NA
  Data.monthly$licence=NA
  Data.monthly$PORT=NA
  Data.monthly$rowid=NA
  
  names(Data.monthly.CAESS)[match(c("sname1","FACTOR"),names(Data.monthly.CAESS))]=c("SNAME","Factor")
   # names(Mesh.monthly)[match(c("month"),names(Mesh.monthly))]=c("MONTH")
  
  ID=match(names(Data.monthly),names(Data.monthly.CAESS))
  Data.monthly.CAESS=Data.monthly.CAESS[,ID]
  Data.monthly=rbind(Data.monthly,Data.monthly.CAESS)
}


#Daily logbooks
  #Select previous daily logbooks (static) or Eva's dynamic data dump
#Get.Daily.Logbook="NO"  #use Rory's files
Get.Daily.Logbook="YES"  #use SADA's annual extraction
if(Get.Daily.Logbook=="YES")
{
  #Daily records 
  Daily.source="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/Shark/02. Data/"
  Daily.file="1. Shark Daily Logbook Data.xlsx"
  Daily.dat=paste(Daily.source,Daily.file,sep="")
  channel <- odbcConnectExcel2007(Daily.dat)
  Data.daily<- sqlFetch(channel,"CATCH", colnames = F)
  close(channel)
  

  #teps
  TEPS.current<- read.csv(paste(Current.yr.dat,"/TEPS_PROTECTEDSP.csv",sep=""),stringsAsFactors=F)
  Comments.TEPS.current<- read.csv(paste(Current.yr.dat,"/TEPS_Comments.csv",sep=""),stringsAsFactors=F)
  if(!is.na(match("DSNo",names(Comments.TEPS.current))))
  {
    names(Comments.TEPS.current)[match("DSNo",names(Comments.TEPS.current))]="DailySheetNumber"
  }
  
  TEPS.current=TEPS.current%>%left_join(Comments.TEPS.current,
                                        by=c("DailySheetNumber","fishery","vessel","year","month"))
  if(!is.na(match("financial year",names(TEPS.current))))
  {
    names(TEPS.current)[match("financial year",names(TEPS.current))]="finyear"
    names(TEPS.current)[match("financial.year",names(TEPS.current))]="finyear"
  }
  
  #other fisheries
  #note: can extract from FishCube running this query:
  # http://F01-FIMS-WEBP01/FishCubeWA/Query.aspx?CubeId=CommercialDPIRDOnly&QueryId=8729f44a-3f88-4fb2-80d3-81a80aad9734
  #     however, since 2017 'other.fishery.catch' is provided by Paul F in the updated Data.monthly
  #Other.fishery.catch=read.csv(paste(Current.yr.dat,"/otherSH.csv",sep=""),stringsAsFactors=F)      
  
  #Catch price
  PRICES=read.csv(paste(Current.yr.dat,"/PriceComparison.csv",sep=""),stringsAsFactors=F)  #update 2016-17
}
if(Get.Daily.Logbook=="NO")
{
  channel <- odbcConnectExcel2007("2009_10_(Jan_2011)") 
 Data.daily.2006.07<- sqlFetch(channel,"2006-07")
 Data.daily.2007.08<- sqlFetch(channel,"2007-08")
 Data.daily.2008.09<- sqlFetch(channel,"2008-09")
 Data.daily.2009.10<- sqlFetch(channel,"2009-10")
 close(channel) 
 
 channel <- odbcConnectExcel2007("2012_Sharkfigs")
 Data.daily.2010.11<- sqlFetch(channel,"DATA")
 close(channel)
 
 
 #Daily records 2011-12 (for 2013 report)
 channel <- odbcConnectAccess2007("2011_12/updated/Shark1213.mdb")
 Data.daily.2011.12<- sqlFetch(channel,"shark", colnames = F)
 close(channel)
 Data.daily.2011.12=subset(Data.daily.2011.12,finyear=="2011-12")
 
   #teps 
 channel <- odbcConnectExcel2007("2011_12/comments&PS-1213")
 TEPS.2011_12<- sqlFetch(channel,"protectedSp")
 Comments.TEPS.2011_12<- sqlFetch(channel,"comments")
 close(channel)
 colnames(Comments.TEPS.2011_12)[match("DSNo",names(Comments.TEPS.2011_12))]="DailySheetNumber"
 TEPS.2011_12=merge(TEPS.2011_12,Comments.TEPS.2011_12,by=c("DailySheetNumber","fishery","vessel","year","month"),all.x=T)
 
 
 
 #Daily records 2012-13 (2014 report)
 channel <- odbcConnectAccess2007("2012_13/Shark2013.mdb")
 Data.daily.2012.13<- sqlFetch(channel,"sharklog", colnames = F)
 close(channel)
 Data.daily.2012.13=subset(Data.daily.2012.13,finyear=="2012-13")
 
 #teps
 channel <- odbcConnectExcel2007("2012_13/TEPS")
 TEPS.2012_13<- sqlFetch(channel,"PROTECTEDSP")
 Comments.TEPS.2012_13<- sqlFetch(channel,"comments")
 close(channel)
 colnames(Comments.TEPS.2012_13)[match("DSNo",names(Comments.TEPS.2012_13))]="DailySheetNumber"
 TEPS.2012_13=merge(TEPS.2012_13,Comments.TEPS.2012_13,by=c("DailySheetNumber","fishery","vessel","year","month"),all.x=T)
 names(TEPS.2012_13)[match("financial year",names(TEPS.2012_13))]="finyear"
 
 
 #Other fisheries
 Other.fishery.catch.2011.12=read.csv("2011_12/SOF1213-others.csv")     #2011_12
 Other.fishery.catch.2012.13=read.csv("2012_13/SOF1214-others.csv")      #2012_13      
}


#Historic TEP data (Prior to 2012 report)
Results.pre.2013=read.csv("Historic/Historic.res.csv")
TEPS.pre.current=read.csv("Historic/Historic.TEPS.res.csv",stringsAsFactors =F)
Spec.catch.zone.pre.2013=read.csv("Historic/Spec.catch.zone.csv")  


#Charter operators                            
#note: each year ask Rhonda to do a data extraction from FishCube
#       species identification is an issue
channel <- odbcConnectExcel2007("Charter/Charter.xlsx") 
Charter.fish.catch<- sqlFetch(channel,"Data")
close(channel)
charter.mean.w=read.csv('Charter/Sp_mean_weight.csv')    
charter.blk=read.csv('Charter/Centroid_Grid_Block_5NM_All.csv')


#Rec fisheries
#note:  iSurvey 2011-12 Tables 7-10 (estimated retained number of shark and ray individuals in 2011-12,boat-based
#       phone diary survey (Ryan et al 2013)
Rec.fish.catch=read.csv("Recreational/I.Survey.csv")


#SHAPE FILE PERTH ISLANDS
# PerthIs=read.table("C:/Matias/Data/Mapping/WAislandsPointsNew.txt", header=T)
# Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
# Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))

#bathymetry data
#    bathymetry data downloaded from http://topex.ucsd.edu/cgi-bin/get_data.cgi (Topography option)
Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi")
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)


#Block10 locations
BlOCK_10=read.csv("C:/Matias/Data/Mapping/Blocks_10NM.csv")


#Weight ranges
Wei.range=read.csv("C:/Matias/Data/Length_Weights/Data.Ranges.csv")
Wei.range.names=read.csv("C:/Matias/Data/Length_Weights/Species.names.csv")


#Conditions
Conditions=read.csv("C:/Matias/Data/Catch and Effort/Conditions.csv")


#Rory's manual changes to netlen and nlines
Rory_Alex_net_val=read.csv("C:/Matias/Data/Catch and Effort/Rory_Alex_net/Book2.csv",stringsAsFactors=F)




##################--- PARAMETERS SECTION ---#######################


#control if inspecting the new data 
if(First.run=="NO")Inspect.New.dat="NO"
if(First.run=="YES")Inspect.New.dat="YES"


#control how to aggregate effort
#note: define whether to use DATE or ID (i.e. the combination of SNo and DSNo) to aggregate by MAX?
Use.Date="YES"    #Advice from Rory McAuley (aggregating by DATE)
#Use.Date="NO"     # aggregating by SNo and DSNo 

#Control if doing figures for cpue standardisation paper
plot.cpue.paper.figures="NO"

#Control if doing figures for Mean weight paper
do.mean.weight.figure="NO"

#control if doing State of Fisheries
do.SoFAR="YES"
#Use.Previos.Sofar="YES"   #Select YES if attaching previous Sofar data to current year #     
Use.Previos.Sofar="NO"   #From 2017 SOFAR, set Use.Previos.Sofar="NO"

Report.stand.cpue="YES"    #set to "YES" if reporting standardised cpue 
#Report.stand.cpue="NO"


#Control if daily logbooks are reapportioned 
#note: reapportion rules were designed for monthly returns
Reapportion.daily="NO"

#Control if reapportioning catch from non-TDGDLF
reapportion.ktch.other.method="NO"

#Control if exporting data for Reference Point paper
do.Ref.Points="NO"

#Control if doing movies
do.exploratory="NO"

#Control if exploring data for catch composition sampling
explore.Catch.compo="NO" 

#Control what data requests are done
do.audit="NO"
do.Russels="NO"
do.Carlys="NO"
do.Alexs="NO"
do.Garys="NO"
do.whiskerys="NO"
do.Hammerheads="NO"
do.Jodies.ASL="NO"
do.Jodies.scalies="NO"
Extract.data.FishCUBE="NO"
do.Steves="NO"
do.AMM.2013="NO"
Export.SAFS="NO"
do.Perth.metro.closure="NO"
do.Heather.request="NO"
do.Jeffs="NO"
do.Dave.F=do.Jeff.N=do.Paul.L=do.SoFAR
do.ABARES="NO"
do.Carlie.Telfer="NO"
do.Cetacean.Inter="NO"
do.Jodies.metro.closure="NO"
do.Tim_N.fin="NO"
do.Tim_ASL="NO"
do.Adrian.Gleiss="NO"
do.Alexandra.Hoschke="NO"
ASL.compensation="NO"
do.CITES="NO"
do.Nick_mesh.size.WCDGDLF="NO"
do.ASL.action.2018="NO"
do.financial.ass="NO"
do.Parks.Australia="NO"


#Spatial range TDGDLF
TDGDLF.lat.range=c(-26,-40)

#Fin to total weight ratio
Percent.fin.of.livewt=0.03

#Species definition
Shark.species=5001:24900
Ray.species=25000:31000
Scalefish.species=188000:599001
Indicator.species=c(17001,17003,18001,18003,18007)
names(Indicator.species)=c("Gummy","Whiskery","Bronzy","Dusky","Sandbar")
TARGETS=list(whiskery=17003,gummy=17001,dusky=c(18001,18003),sandbar=18007) 
N.species=length(TARGETS)
Fix.species=c(18003,17001,17003,22999)  #species that need catch reapportioning

#Catch and effort misreporting  as adviced by Colin Simpfendorfer (2001) and Rory McAuley 
Inc.per=1.05  #5% percent increase in catch and effort prior 1990                    

#Assumed increase in fishing power as adviced by Colin Simpfendorfer (2001) and Rory McAuley 
Fish.Pow=.02  #2% annual (i.e. 2%, 4%, 6%, etc) increase in fishing power prior 1994                     

#Maximum net length in TDGDLF
Net.max=8500    #Rory pers comm
#Net.max=12000  #maximum possible net length   (McAuley et al 2005)

#Boundary blocks
Boundary.blk=c(34160,35160,36160)

#Commercial species regions for effective effort calculation (McAuley 2005)
    #Rory's rule 5a-5k
Dusky.range=c(-28,120)
Sandbar.range=c(-26,118)
Whiskery.range=c(-28,129)
Gummy.range=c(116,129)
Dist.range=list(Dusky=Dusky.range,Sandbar=Sandbar.range,Whiskery=Whiskery.range,Gummy=Gummy.range)

#Maximum possible weights for commercial species
Wei.range=Wei.range%>%left_join(Wei.range.names,by=c("Sname"))
#Wei.range=merge(Wei.range,Wei.range.names,by="Sname",all.x=T)
Wei.range=subset(Wei.range,!(is.na(SPECIES)|is.na(TW.min)|is.na(TW.max)))
max.w.whis=Wei.range[Wei.range$SPECIES==17003,]$TW.max   
max.w.gum=Wei.range[Wei.range$SPECIES==17001,]$TW.max
max.w.dus=Wei.range[Wei.range$SPECIES==18003,]$TW.max
max.w.san=Wei.range[Wei.range$SPECIES==18007,]$TW.max


 #how to export figures
Do.tiff="NO"
Do.jpeg="YES"



##################--- DATA MANIPULATION SECTION ---##############

#Non-commercial Fisheries
 #a. recreational
Rec.fish.catch$FinYear=as.character(Rec.fish.catch$Survey)
Rec.fish.catch=subset(Rec.fish.catch,!FinYear=="Grand Total")
names(Rec.fish.catch)[match("Estimated.Kept.Catch..by.numbers.",names(Rec.fish.catch))]="Kept.Number"
names(Rec.fish.catch)[match("Estimated.Released.Catch..by.numbers.",names(Rec.fish.catch))]="Rel.Number"
names(Rec.fish.catch)[match("Species.Common.Name",names(Rec.fish.catch))]="Common.Name"
names(Rec.fish.catch)[match("Species.Scientific.Name",names(Rec.fish.catch))]="Scientific.Name"
if(is.factor(Rec.fish.catch$Kept.Number)) Rec.fish.catch$Kept.Number <- as.numeric(gsub(",","",Rec.fish.catch$Kept.Number))
if(is.factor(Rec.fish.catch$Rel.Number)) Rec.fish.catch$Rel.Number <- as.numeric(gsub(",","",Rec.fish.catch$Rel.Number))

  #b. charter boats
names(Charter.fish.catch)[match("TO Zone",names(Charter.fish.catch))]="Zone"
names(Charter.fish.catch)[match("Calendar Year",names(Charter.fish.catch))]="Year"
names(Charter.fish.catch)[match("05x05NM Block",names(Charter.fish.catch))]="Block5"
names(Charter.fish.catch)[match("Species Common Name",names(Charter.fish.catch))]="SNAME"
names(Charter.fish.catch)[match("Fish Kept Count",names(Charter.fish.catch))]="Kept_N"
names(Charter.fish.catch)[match("Fish Released Count",names(Charter.fish.catch))]="Released_N"
names(Charter.fish.catch)[match("Licence Count",names(Charter.fish.catch))]="Number_boats"
names(Charter.fish.catch)[match("Fishing Trip Count",names(Charter.fish.catch))]="Trips_N"
Charter.fish.catch=subset(Charter.fish.catch,!Zone=='Grand Total')
Charter.fish.catch=Charter.fish.catch%>%left_join(charter.mean.w,by=c("SNAME"))
#Charter.fish.catch=merge(Charter.fish.catch,charter.mean.w,by="SNAME",all.x=T)
Charter.fish.catch$Kept_w=Charter.fish.catch$Kept_N*Charter.fish.catch$Mean_weight_kg
Charter.fish.catch$Released_w=Charter.fish.catch$Released_N*Charter.fish.catch$Mean_weight_kg
Charter.fish.catch$Year=year(Charter.fish.catch$Year)
Charter.fish.catch$Month=month(Charter.fish.catch$Month)
Charter.fish.catch=Charter.fish.catch%>%left_join(charter.blk,by=c("Block5"="BlockNo"))
#Charter.fish.catch=merge(Charter.fish.catch,charter.blk,by.x="Block5",by.y="BlockNo",all.x=T)
Tab.chart.blck.sp=aggregate(cbind(Kept_N,Released_N)~SNAME+Block5,Charter.fish.catch,sum)
Tab.chart.blck.sp_w=aggregate(cbind(Kept_w,Released_w)~SNAME+Block5,Charter.fish.catch,sum)
Tab.chart.yr.sp_w=aggregate(cbind(Kept_w,Released_w)~SNAME+Year,Charter.fish.catch,sum)
par(mfcol=c(2,1),mar=c(2,2,1,1),oma=c(1,1,1,1),las=1,mgp=c(1,.5,0))
with(aggregate(cbind(Kept_N,Released_N)~LAT+LONG,Charter.fish.catch,sum),
{
       plot(LONG,LAT,pch=19,col="steelblue",cex=fn.scale(Kept_N,3),
            main=paste("All years. Kept numbers","(n=",sum(Kept_N),")"),ylab='',xlab='')
       plot(LONG,LAT,pch=19,col="steelblue",cex=fn.scale(Released_N,3),
            main=paste("All years. Released numbers","(n=",sum(Released_N),")"),ylab='',xlab='')
     })
yr=as.numeric(substr(Current.yr,1,4))
with(aggregate(cbind(Kept_N,Released_N)~LAT+LONG,subset(Charter.fish.catch,Year==yr),sum),
{
       plot(LONG,LAT,pch=19,col="steelblue",cex=fn.scale(Kept_N,3),
            main=paste(yr,"Kept numbers","(n=",sum(Kept_N),")"),ylab='',xlab='')
       plot(LONG,LAT,pch=19,col="steelblue",cex=fn.scale(Released_N,3),
            main=paste(yr,"Released numbers","(n=",sum(Released_N),")"),ylab='',xlab='')
     })

#Quick check to determine samples size catch age composition
if(explore.Catch.compo=="YES")
{
  hnd.compo="C:/Matias/Fieldwork and workplans/Catch age composition/"
    #Albany
  fun.sample.catch.age.com=function(d)
  {
    d=d %>%mutate_all( funs(as.character(.)), names( .[,sapply(., is.factor)] ))%>%
      mutate(Name=ifelse(species==17003,"Whiskery",
                         ifelse(species==17001,"Gummy",
                                ifelse(species==18003,"Dusky",
                                       ifelse(species==18007,"Sandbar",NA)))))
    
    Tab1=with(d[!duplicated(d$TSNo),],table(finyear,vessel))
    jpeg(file=paste(hnd.compo,unique(d$port),"_vessel trips.jpg",sep=""),width=2400,height=2400,units="px",res=300)
    par(las=1)
    barplot(Tab1,horiz=T,legend.text=rownames(Tab1),xlab="Number of trips",
            args.legend=list(x='bottomright',bty='n',cex=1.5))
    box()
    dev.off()
    
    d$nfish=as.numeric(d$nfish)
    
    dd=d %>% group_by(Name,finyear) %>%
      summarize(Sum=sum(nfish))%>%
      mutate(Yr=as.numeric(substr(finyear,1,4)))%>%
      as.data.frame
    jpeg(file=paste(hnd.compo,unique(d$port),"_total.fish.jpg",sep=""),width=2400,height=2400,units="px",res=300)
    par(cex.axis=1.35, cex.lab=1.5)
    ggplot(data = dd, aes(x = finyear, y = jitter(Sum,50))) + 
      geom_point(aes(colour  =Name),size=5)+ labs(y = "Total # individuals")+ theme(axis.title = element_text(face="bold", size=22))+
      theme(axis.text=element_text(size=18))+theme(legend.text=element_text(size=16))
    dev.off()
    
    Uni=unique(d$Name)
    for(u in 1:length(Uni))
    {
      jpeg(file=paste(hnd.compo,unique(d$port),"_nfish_",Uni[u],".jpg",sep=""),width=2400,height=2400,units="px",res=300)
      ggplot(data = subset(d,Name==Uni[u]), aes(x = finyear, y = nfish)) + 
        geom_boxplot(aes(fill =vessel ), width = 0.8) + theme_bw()+
        theme(axis.title = element_text(face="bold", size=22))+
        theme(axis.text=element_text(size=18))+theme(legend.text=element_text(size=16))
      dev.off()
    }
    
    # dd=d %>% group_by(Name,vessel,finyear) %>%
    #   summarize(mean=mean(nfish),
    #          sd=sd(nfish))%>%
    #   as.data.frame
    # for(u in 1:length(Uni))
    # {
    #   xx=subset(dd,Name==Uni[u])
    #   Mn=xx %>% select(vessel,finyear,mean)%>% spread(finyear,mean)
    #   SD=xx %>% select(vessel,finyear,sd)%>% spread(finyear,sd)
    # }
  }
  fun.sample.catch.age.com(d=Data.daily %>%
                             filter(port=="ALBANY" & finyear%in%c('2014-15','2015-16','2016-17') &
                                      species%in%c(17003,17001,18003,18007))%>% 
                             select(port,species,nfish,vessel,TSNo,finyear))
}



#SECTION A. ---- MONTHLY RECORDS ----

# A.1. Add year if missing
Data.monthly$Split.Year.1=as.numeric(sapply(strsplit(as.character(Data.monthly$FINYEAR),"-"), "[", 1))
Data.monthly$Split.Year.2=as.numeric(sapply(strsplit(as.character(Data.monthly$FINYEAR),"-"), "[", 2))
Data.monthly$Split.Year.2=with(Data.monthly,ifelse(Split.Year.2<75,Split.Year.2+2000,Split.Year.2+1900))
Data.monthly$YEAR.c=with(Data.monthly,ifelse(is.na(YEAR) & MONTH>=7,Split.Year.1,
                          ifelse(is.na(YEAR) & MONTH<7,Split.Year.2,YEAR)))


# A.2. Remove records post June 2006 if using Table 89.d only
Data.monthly$ID=with(Data.monthly,ifelse(YEAR.c>2006,"drop","keep"))
Data.monthly$ID=with(Data.monthly,ifelse(YEAR.c==2006 & MONTH>6,"drop",ID))
if(Get.CAESS.Logbook=="NO") Data.monthly=subset(Data.monthly,ID=="keep")


# A.3. Remove dummies
Data.monthly=Data.monthly[,-match(c("Split.Year.1","Split.Year.2","ID"),names(Data.monthly))]


# A.4. Sort by year and month and vessel
Data.monthly=Data.monthly[order(Data.monthly$YEAR.c,Data.monthly$MONTH,Data.monthly$VESSEL),]


# A.4.2 Create dummy to group same return (finyear-month-vessel-method-block)             # REVIEW RORY
Data.monthly$Same.return=with(Data.monthly,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
Data.monthly$Same.return.SNo=Data.monthly$Same.return
Data.monthly$TYPE.DATA="monthly"


#remove fins and livers to avoid duplication 
  #note: id records were only fins or livers were reported.
  # for this, check with species were wrongly assigned a liver
  # or fin code and remove records from those where liver/fin and
  # WF or WH weights were reported
A=subset(Data.monthly,CONDITN%in%c("LI","FI") & 
           !SPECIES%in%c(22997,22998))
SPEC.livr.fin=unique(A$SPECIES)
A=subset(Data.monthly,SPECIES%in%SPEC.livr.fin)
Data.monthly=subset(Data.monthly,!SPECIES%in%SPEC.livr.fin)
keep.these=vector('list',length(SPEC.livr.fin))
for(s in 1:length(SPEC.livr.fin))
{
  x=subset(A,SPECIES==SPEC.livr.fin[s])
  Tbl.x=with(x,table(Same.return,CONDITN))
  Tbl.x=as.data.frame.matrix(Tbl.x) 
  Tbl.x$Sums=rowSums(Tbl.x)
  if(!is.null(Tbl.x$LI)) KEEP=subset(Tbl.x,(Sums==1 & FI==1)| (Sums==1 & LI==1)) else
  KEEP=subset(Tbl.x,(Sums==1 & FI==1))
  keep.these[[s]]=row.names(KEEP)
}
keep.these=unlist(keep.these)
keep.these=keep.these[!duplicated(keep.these)]
A$KEEP=with(A,ifelse(Same.return%in%keep.these,"KEEP","DROP"))
A$KEEP=with(A,ifelse(!CONDITN%in%c("FI","LI"),"KEEP",KEEP))
Data.monthly$KEEP=with(Data.monthly,
        ifelse(CONDITN%in%c("FI","LI"),"DROP","KEEP"))
Data.monthly=rbind(Data.monthly,A)
Data.monthly=subset(Data.monthly,!(SPECIES%in%c(22997,22998)))
Data.monthly=subset(Data.monthly,KEEP=="KEEP")
Data.monthly$CONDITN=with(Data.monthly,
      ifelse(CONDITN%in%c("FI","LI"),"WF",CONDITN))
rm(A)
Data.monthly=Data.monthly[,-match("KEEP",names(Data.monthly))]


# A.4.3. Create copy of original file
Data.monthly.original=Data.monthly
Data.monthly.original$LAT=-as.numeric(substr(Data.monthly.original$BLOCKX,1,2))
Data.monthly.original$LONG=100+as.numeric(substr(Data.monthly.original$BLOCKX,3,4))

fn.chk.ktch=function(d1,d2,VAR1,VAR2)
{
  d1sum=round(sum(d1[,match(VAR1,names(d1))],na.rm=T))
  d2sum=round(sum(d2[,match(VAR2,names(d2))],na.rm=T))
  if(d1sum>d2sum)Message="d1 has more catch"
  if(d1sum<d2sum)Message="d2 has more catch"
  if(d1sum==d2sum)Message="d1 and d2 have same catch"
  print(Message)
  
  A=aggregate(d1[,match(VAR1,names(d1))]~FINYEAR+VESSEL,d1,sum)
  names(A)[3]="weight.original"
  B=aggregate(d2[,match(VAR2,names(d2))]~FINYEAR+VESSEL,d2,sum)
  names(B)[3]="weight.changed"
  D=A%>%full_join(B,by=c("FINYEAR","VESSEL"))
  #D=merge(A,B,by=c("FINYEAR","VESSEL"))
  D$delta.w=round(D$weight.original-D$weight.changed,0)
  discrepancy=subset(D,!delta.w==0)
  print(paste("Catch by year-vessel has this many discrepancies=",
              nrow(discrepancy),"records"))
}
fn.chk.ktch(d1=Data.monthly.original,
      d2=subset(Data.monthly,FINYEAR%in%Data.monthly.original$FINYEAR),
      VAR1="LIVEWT",VAR2="LIVEWT")


# A.5. Create variables

  #define estuaries
Estuaries=c(95010,95020,95030,95040,95050,95060,95070,95080,95090,85010,85020,85030,85040,85050,85060,
            85070,85080,85090,85100,85110,85130,85990)
Bays=96000:98000
Shark.Bay=96021
Abrolhos=97011
Geographe.Bay=96010 
Cockburn.sound=96000
King.George.sound=96030 

  #Idenfity estuaries
Data.monthly$Estuary=with(Data.monthly,ifelse(BLOCKX%in%Estuaries,"YES","NO"))
test.Estuaries=subset(Data.monthly,BLOCKX%in%Estuaries)
Table.Estuary=aggregate(LIVEWT~FINYEAR+SPECIES,test.Estuaries,sum)  
#note: there is shark catch in estuaries, must keep records for total catch taken


  # A.5.1 Lat and Long of block (top left corner)

#get lat and long from block but first dummy the bays                                                                      
Data.monthly$blockxFC=Data.monthly$BLOCKX
Data.monthly$BLOCKX=with(Data.monthly,ifelse(BLOCKX%in%c(Shark.Bay),25120,     
                    ifelse(BLOCKX%in%c(96022,96023),26131,
                    ifelse(BLOCKX%in%c(Abrolhos),27132,                        
                    ifelse(BLOCKX%in%c(97012,97013),28132,       
                    ifelse(BLOCKX%in%c(97014,97015),29132,
                    ifelse(BLOCKX%in%c(Geographe.Bay),33151,                        
                    ifelse(BLOCKX%in%c(Cockburn.sound),32150,                        
                    ifelse(BLOCKX%in%c(King.George.sound),35181,BLOCKX)))))))))         

Data.monthly$LAT=-as.numeric(substr(Data.monthly$BLOCKX,1,2))
Data.monthly$LONG=100+as.numeric(substr(Data.monthly$BLOCKX,3,4))

Data.monthly$BLOCKX=Data.monthly$blockxFC   #reset block

# adjust special cells to actual lat and long                                              
#note: use celing for lat and floor for long to allocate lat and long of north-west corner
Data.monthly$LAT=with(Data.monthly,ifelse(BLOCKX%in%c(96021),-25,
                 ifelse(BLOCKX%in%c(96022,96023),-26,
                ifelse(BLOCKX%in%c(97011),-27,
                ifelse(BLOCKX%in%c(97012,97013),-28,       
                ifelse(BLOCKX%in%c(97014,97015),-29,
                ifelse(BLOCKX%in%c(96010),-33,
                ifelse(BLOCKX%in%c(96000),-33,
                ifelse(BLOCKX%in%c(96030,95050,95090,95040),-35,
                ifelse(BLOCKX%in%c(85030),-34.4586,
                ifelse(BLOCKX%in%c(85050),-34.2853,
                ifelse(BLOCKX%in%c(85080),-33.917,
                ifelse(BLOCKX%in%c(95010),-31.9554,
                ifelse(BLOCKX%in%c(95020),-32.5167,
                ifelse(BLOCKX%in%c(95030),-33.3503,
                ifelse(BLOCKX%in%c(95060),-34.9953,
                ifelse(BLOCKX%in%c(95070),-34.9731,
                ifelse(BLOCKX%in%c(95080),-34.9278,
                ifelse(BLOCKX%in%c(85990),NA,
                ifelse(BLOCKX%in%c(85130),-34.9278,LAT
                ))))))))))))))))))))

Data.monthly$LAT=ceiling(Data.monthly$LAT)

Data.monthly$LONG=with(Data.monthly,ifelse(BLOCKX%in%c(96021),113,
                  ifelse(BLOCKX%in%c(96022,96023),113,
                  ifelse(BLOCKX%in%c(97011),113,
                  ifelse(BLOCKX%in%c(97012,97013),113,       
                  ifelse(BLOCKX%in%c(97014,97015),113,
                  ifelse(BLOCKX%in%c(96010),115,
                  ifelse(BLOCKX%in%c(96000),115,
                  ifelse(BLOCKX%in%c(96030,95050,95090,95040),118,
                  ifelse(BLOCKX%in%c(85030),118.8897,
                  ifelse(BLOCKX%in%c(85050),119.4819,
                  ifelse(BLOCKX%in%c(85080),120.050,
                  ifelse(BLOCKX%in%c(95010),115.8585,
                  ifelse(BLOCKX%in%c(95020),115.7167,
                  ifelse(BLOCKX%in%c(95030),115.6494,
                  ifelse(BLOCKX%in%c(95060),117.4117,
                  ifelse(BLOCKX%in%c(95070),116.9744,
                  ifelse(BLOCKX%in%c(95080),116.4489,
                  ifelse(BLOCKX%in%c(85990),NA,
                  ifelse(BLOCKX%in%c(85130),116.4489,LONG
                  ))))))))))))))))))))

Data.monthly$LONG=floor(Data.monthly$LONG)


# A.7. Fix condition
Data.monthly$CONDITN=as.character(Data.monthly$CONDITN)
Data.monthly$CONDITN=with(Data.monthly,ifelse(is.na(SPECIES),"NIL",CONDITN))


# A.8. Add bioregion and zone                                                               
Data.monthly$LONG=with(Data.monthly,ifelse(BLOCKX==85990,NA,
                    ifelse(BLOCKX==31300,130,
                    ifelse(BLOCKX==32300,130,LONG))))
Data.monthly$LAT=with(Data.monthly,ifelse(BLOCKX==85990,NA,
                    ifelse(BLOCKX==31300,-31,
                    ifelse(BLOCKX==32300,-32,LAT))))

Data.monthly$Bioregion=as.character(with(Data.monthly,ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
             ifelse(LONG<115.5 & LAT<=(-27),"WC",
             ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
             ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))))

Data.monthly$Bioregion=with(Data.monthly,
            ifelse(Bioregion=="SC"& LAT>(-34) & LONG <115.91 ,"WC",Bioregion))

Data.monthly$zone=as.character(with(Data.monthly,ifelse(LONG>=116.5 & LAT<=(-26),"Zone2",
                  ifelse(LONG<116.5 & LAT<=(-33),"Zone1",
                  ifelse(LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                  ifelse(LAT>(-26) & LONG<114,"Closed",
                  ifelse(LAT>(-23) & LONG>=114 & LONG<123.75,"North",
                  ifelse(LAT>(-23) & LONG>=123.75,"Joint",NA))))))))


# A.9. Create Monthly effort dataset
Data.monthly$NETLEN=with(Data.monthly,ifelse(METHOD%in%c("LL","HL","DL")&NETLEN>0,NA,NETLEN))
Effort.vars=c("FDAYS","BDAYS","HOURS","HOOKS","SHOTS","NETLEN")
Effort.monthly=Data.monthly[,match(c("FINYEAR","YEAR","MONTH","zone","VESSEL","METHOD","BLOCKX",Effort.vars,
                                      "YEAR.c","Same.return","LAT","LONG"),names(Data.monthly))]
  #add dummy nlines
Effort.monthly$nlines=NA
Effort.vars=c(Effort.vars,"nlines")

    #add variables for FishCUBE
Data.monthly$FisheryZone=Data.monthly$zone
Data.monthly$FisheryCode=Data.monthly$fishery
Data.monthly$Landing.Port=Data.monthly$PORT
Data.monthly$LatMin=NA
Data.monthly$LongMin=NA
Data.monthly$RSCommonName=NA
Data.monthly$RSSpeciesId=NA
Data.monthly$Block=NA
Data.monthly$LatFC=NA
Data.monthly$LongFC=NA

    #add var for merging with daily
Data.monthly$nfish=NA
Data.monthly=Data.monthly[,-match("fishery",names(Data.monthly))]




#SECTION B. ---- DAILY LOGBOOKS ----

#simple financial assessment
if(do.financial.ass=="YES")
{
  b=aggregate(livewt~vessel+species,subset(Data.daily,finyear==Current.yr,select=c(livewt,vessel,species)),sum,na.rm=T)
  bb=merge(b,subset(PRICES,select=c(SPECIES,uv1516)),by.x="species",by.y="SPECIES",all.x=T)
  bb$Revenue_annual=bb$livewt*bb$uv1516
  dd=aggregate(cbind(Revenue_annual,livewt)~vessel,bb,sum,na.rm=T)
  v=subset(Data.daily,vessel%in%dd$vessel & finyear==Current.yr,select=c(TSNo,vessel,BoatName,MastersName,crew,date,fdays,bdays))
  d=v[!duplicated(v$vessel),]
  vv=aggregate(crew~vessel,v,mean,na.rm=T)
  vv=merge(vv,subset(d,select=c(vessel,BoatName,MastersName)),by="vessel")
  
  #correct fishing days
  fDaYs=v
  fDaYs$dummy=with(fDaYs,paste(vessel,TSNo,date))
  fDaYs=fDaYs[!duplicated(fDaYs$dummy),]
  fDaYs$fdays.c=1
  fDaYs=aggregate(fdays.c~vessel,fDaYs,sum)
  
  A=merge(dd,vv,by='vessel')
  A=merge(A,fDaYs,by='vessel')
  
  A$Revenue_per_fishing_day=A$Revenue_annual/A$fdays.c
  A=A[order(A$Revenue_per_fishing_day),]
  A$Costs="?"
  A=A[,match(c("vessel","BoatName","MastersName","crew","fdays.c" ,"livewt","Revenue_annual","Revenue_per_fishing_day","Costs"),names(A))]
  names(A)[match(c("crew","fdays.c","livewt","Revenue_annual","Revenue_per_fishing_day"),names(A))]=
    c("crew_average","fishing_days_annual","livewt_annual(kg)","Revenue_annual (AUD$)","Revenue_per_fishing_day (AUD$)")
  write.csv(A,paste("C:/Matias/Analyses/Catch and effort/Annual.revenue.",Current.yr,".csv",sep=""),row.names=F)
  
}

#Create backup file   
Data.daily.original=Data.daily
names(Data.daily.original)[match(c("finyear","month","vessel","method","bdays","hours","hooks",
               "shots","netlen","LatDeg"),names(Data.daily.original))]=c("FINYEAR","MONTH","VESSEL","METHOD",
                     "BDAYS","HOURS","HOOKS","SHOTS","NETLEN","LAT")
Data.daily.original$LAT=(-Data.daily.original$LAT)
Data.daily.original$Same.return.SNo=with(Data.daily.original,paste(SNo,DSNo,TSNo))


#note: a unique shot (session) is the combination of Sno (shot num), DSNo (daily sheet num) and     
#       TSNo (trip sheet num),i.e. the variable "Session ID"
Session.vars=c("SNo","DSNo","TSNo")
FishCube.vars=c("RSSpeciesId","RSCommonName","Block","LatFC","LongFC","blockxFC","licence","port")
This=c("finyear","year","month","day","date","zone","depthMax","vessel","fdays","method","blockx","block10","bdays",
       "hours","hooks","shots","nlines","netlen","species","sname1","nfish","livewt","conditn","landwt",
       "factor","LatDeg","LongDeg","LatMin","LongMin",Session.vars,"Bioregion","fishery",FishCube.vars)

#create file for checking mesh size
Mesh.size=subset(Data.daily,method=="GN",select=c(DSNo,TSNo,SNo,vessel,date,finyear,method,
                            netlen,hours,shots,mshigh,mslow,LatDeg,LongDeg,zone))

  #Select which file to use
if(Get.Daily.Logbook=="YES")
{
  Data.daily.1=Data.daily
  Data.daily$Bioregion=Data.daily$bioregion
  kk=match("day",names(Data.daily))
  if(is.na(kk))Data.daily$day=day(Data.daily$date)
  kk=match("date",names(Data.daily))
  if(is.na(kk))Data.daily$date=with(Data.daily,as.POSIXct(paste(year,"-",month,"-",day,sep="")))
  Data.daily$finyear=as.character(Data.daily$finyear)
  if(is.na(match("hooks",names(Data.daily)))) Data.daily$hooks=0
  #if(is.na(match("nlines",names(Data.daily)))) Data.daily$nlines=0
  
  #create demersal suite scalefish species
  Suite=subset(Data.daily,type=="demersal")           
  Suite=unique(Suite$species)
  
  #Keep original lat and long for FISHCUBE
  Data.daily$LatFC=as.character(Data.daily$Lat)
  Data.daily$LongFC=as.character(Data.daily$Long)

  #Keep original blockx FISHCUBE
  Data.daily$blockxFC=Data.daily$blockx
  
  Data.daily=Data.daily[,match(c(This,"rowid","flagtype"),names(Data.daily))]
  Data.daily.FC.NA.sp=subset(Data.daily,species==99999 | is.na(species))
  
  Data.daily.FC.2005_06=subset(Data.daily,finyear=="2005-06")
  Data.daily.FC.2005_06=subset(Data.daily.FC.2005_06,!(species==99999 | is.na(species)))
  
  Data.daily=subset(Data.daily,!finyear=="2005-06")
  

  #Fix some records identified as incorrect in CATCH INSPECTIONS and EFFORT INSPECTIONS  
  #note: the ammended values were provided by data entry girls
  
  id=subset(Data.daily,finyear=="2011-12" & TSNo=="TDGLF8000737" & 
              DSNo=="TDGLF8000734" & sname1=="Blue Morwong (queen snapper)")
  if(nrow(id)>=1)Data.daily=subset(Data.daily,!(finyear=="2011-12" & TSNo=="TDGLF8000737" & 
              DSNo=="TDGLF8000734" & sname1=="Blue Morwong (queen snapper)"))
  
  Data.daily$livewt=with(Data.daily,
    ifelse(finyear=="2011-12" & TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & 
             sname1=="Blue Morwong (queen snapper)"& livewt<700,5.04,
    ifelse(finyear=="2011-12" & TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & 
             sname1=="Blue Morwong (queen snapper)"& livewt>700,10.06,
    ifelse(finyear=="2011-12" & TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & 
             sname1=="Blue Morwong (queen snapper)"& livewt<200,13,
    ifelse(finyear=="2011-12" & TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & 
             sname1=="Blue Morwong (queen snapper)"& livewt>200,26,livewt)))))

  #Fixing issues from reporting combined catch but not properly flagged by data entry girls so Eva's
  # correction for combining catch was not applied
  Data.daily$livewt=with(Data.daily,
    ifelse(species%in%c(17001,17003) & SNo%in%c(2,4,6)&DSNo=="TDGLF8009540" & flagtype==0,
           1.59*3.889706*nfish,
    ifelse(species%in%c(17001,17003) & SNo%in%c(1)&DSNo=="TDGLF8007247" & flagtype==0,
           1.59*4.742857*nfish,
    ifelse(species%in%c(17001,17003)& DSNo=="TDGLF8009946" & flagtype==0,
           1.59*3.517375*nfish,
    ifelse(species%in%c(17001,17003)& TSNo=="TDGLF8008604" & flagtype==0,
           1.59*4.98*nfish,
    ifelse(species%in%c(18003,18007)& DSNo=="TDGLF8009501" & flagtype==0,
           1.59*4.661017*nfish,
    ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009505" & flagtype==0,
           1.59*6.943662*nfish,
    ifelse(species%in%c(18003,18007) &TSNo=="TDGLF8009507" & flagtype==0,
           1.59*7.304348*nfish,
    ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009534" & flagtype==0,
           1.59*14.5*nfish,
    ifelse(species%in%c(18003,18007,18001)& TSNo=="TDGLF8006245" & flagtype==0,
           1.59*7.371429*nfish,
    ifelse(species%in%c(18003,18007,18001)& TSNo=="TDGLF8006248" & flagtype==0,
           1.59*5*nfish,
    ifelse(species%in%c(18003,18001)& TSNo=="TDGLF8009402" & flagtype==0,
           1.59*6.653846*nfish,
    ifelse(species%in%c(18023,18001)& TSNo=="TDGLF8008487" & flagtype==0,
           1.59*4.306569*nfish,
    ifelse(species%in%c(17006,18003,19000)& TSNo=="TDGLF8009537" & flagtype==0,
           1.59*7.911765*nfish*807/933.5883,
    ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009503" & flagtype==0,
           1.59*7.090909*nfish,livewt)))))))))))))))
  
  Data.daily$landwt=with(Data.daily,
     ifelse(species%in%c(17001,17003) & SNo%in%c(2,4,6)&DSNo=="TDGLF8009540" & flagtype==0,
            3.889706*nfish,
     ifelse(species%in%c(17001,17003) & SNo%in%c(1)&DSNo=="TDGLF8007247" & flagtype==0,
            4.742857*nfish,
     ifelse(species%in%c(17001,17003)& DSNo=="TDGLF8009946" & flagtype==0,
            3.517375*nfish,
     ifelse(species%in%c(17001,17003)& TSNo=="TDGLF8008604" & flagtype==0,
            4.98*nfish,
     ifelse(species%in%c(18003,18007)& DSNo=="TDGLF8009501" & flagtype==0,
            4.661017*nfish,
     ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009505" & flagtype==0,
            6.943662*nfish,
     ifelse(species%in%c(18003,18007) &TSNo=="TDGLF8009507" & flagtype==0,
            7.304348*nfish,
     ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009534" & flagtype==0,
           14.5*nfish,
     ifelse(species%in%c(18003,18007,18001)& TSNo=="TDGLF8006245" & flagtype==0,
           7.371429*nfish,
     ifelse(species%in%c(18003,18007,18001)& TSNo=="TDGLF8006248" & flagtype==0,
           5*nfish,
     ifelse(species%in%c(18003,18001)& TSNo=="TDGLF8009402" & flagtype==0,
           6.653846*nfish,
     ifelse(species%in%c(18023,18001)& TSNo=="TDGLF8008487" & flagtype==0,
           4.306569*nfish,
     ifelse(species%in%c(17006,18003,19000)& TSNo=="TDGLF8009537" & flagtype==0,
           7.911765*nfish*807/933.5883,
    ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009503" & flagtype==0,
          7.090909*nfish,landwt)))))))))))))))
  
  
  Data.daily=Data.daily[,-match("flagtype",names(Data.daily))]
  
  Data.daily$nfish=with(Data.daily,
    ifelse(finyear=="2011-12" & DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & 
             sname1=="Shark, Whiskery"& nfish>3,52,
    ifelse(finyear=="2011-12" & DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & 
             sname1=="Shark, Whiskery"& nfish<3,12,nfish)))
  
  Data.daily$netlen=with(Data.daily,
    ifelse(finyear=="2011-12" & vessel=="B 067"& netlen==650,6500,
    ifelse(finyear=="2011-12" & vessel=="B 142"& netlen==540,5400,
    ifelse(finyear=="2011-12" & vessel=="E 009"& netlen==420,4320,
    ifelse(finyear=="2011-12" & vessel=="E 009"& netlen==400,4000,
    ifelse(finyear=="2011-12" & vessel=="E 035"& netlen==380,3800,
    ifelse(finyear=="2010-11" & vessel=="G 297"& netlen==7320,4320,
    ifelse(finyear=="2011-12" & vessel=="F 417"& netlen==6000,600,
    ifelse(finyear%in%c("2010-11","2013-14") & vessel=="F 417"& netlen==7500,750,
    ifelse(finyear=="2013-14" & vessel=="F 417"& netlen==6500,650,netlen))))))))))
  
  Data.daily$shots=with(Data.daily,
    ifelse(finyear%in%c("2009-10","2010-11") & vessel=="F 541"& shots==6,1,
    ifelse(finyear=="2011-12" & vessel=="F 541"& shots==5,1,
    ifelse(finyear=="2010-11" & vessel=="E 056"& shots==10,1,shots))))
  
  Data.daily$nfish=with (Data.daily,
    ifelse(finyear=="2012-13" & species==377004 & vessel=="G 297" & landwt == 92.57,7,nfish))
  
  Data.daily$livewt=with (Data.daily,
    ifelse(finyear=="2012-13" & species==377004 & vessel=="E 035" & landwt == 49,6.7,
    ifelse(finyear=="2012-13" & species==377004 & vessel=="G 297" & landwt == 92.57,216,
    ifelse(finyear=="2012-13" & species==17003 & vessel=="E 059" & landwt == 1769,17,livewt))))
  
  
  Data.daily$netlen=with (Data.daily,
    ifelse(finyear=="2012-13" & DSNo=="TDGLF8002239" &  TSNo=="TDGLF8002240" & netlen==165,3500,
    ifelse(finyear=="2012-13" & TSNo=="TDGLF8002261" & netlen==300,3000,
    ifelse(finyear=="2012-13" & TSNo%in%c("TDGLF8008946","TDGLF8008924") & netlen==400,4000,netlen))))
}
if(Get.Daily.Logbook=="NO")
{
  # B.1. Combine data sets and create copy
  
  #2006-2010
  Data.daily=rbind(Data.daily.2006.07,Data.daily.2007.08,Data.daily.2008.09,Data.daily.2009.10)
  Data.daily.1=Data.daily
  
  #note: SNo is the session number, this get's repeated
  #       DSNo is daily sheet number, this is unique and it groups all SNos within a date
  #       TSNo is trip sheet number, which is the DSNo were the total download weigth of the trip is reported
  Data.daily$nfish=Data.daily$orgnfish
  
  Data.daily=Data.daily[,match(This,names(Data.daily))]
  
  
  
  #2010-2011
  Data.daily.2010.11.1=Data.daily.2010.11
  Data.daily.2010.11$date=with(Data.daily.2010.11,as.POSIXct(paste(year,"-",month,"-",day,sep="")))
  Data.daily.2010.11$hooks=0
  Data.daily.2010.11$nlines=0
  Data.daily.2010.11$Bioregion=Data.daily.2010.11$bioregion
  Data.daily.2010.11=Data.daily.2010.11[,match(This,names(Data.daily.2010.11))]
  
  
  
  #2011-2012
  Data.daily.2011.12$hooks=0
  Data.daily.2011.12$nlines=0
  Data.daily.2011.12$Bioregion=Data.daily.2011.12$bioregion
  Data.daily.2011.12$date=with(Data.daily.2011.12,as.POSIXct(paste(year,"-",month,"-",day,sep="")))
  
  #manually fix some catch records identified as incorrect in "#---- CATCH INSPECTIONS -----" 
  #note: the ammended values were provided by data entry girls
  id=subset(Data.daily.2011.12,TSNo=="TDGLF8000737" & DSNo=="TDGLF8000734" & sname1=="Blue Morwong (queen snapper)")
  if(nrow(id)>=1)Data.daily.2011.12=subset(Data.daily.2011.12,!(TSNo=="TDGLF8000737" & DSNo=="TDGLF8000734" & sname1=="Blue Morwong (queen snapper)"))
  
  Data.daily.2011.12$livewt=with(Data.daily.2011.12,
                                 ifelse(TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & sname1=="Blue Morwong (queen snapper)"& livewt<700,5.04,
                                        ifelse(TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & sname1=="Blue Morwong (queen snapper)"& livewt>700,10.06,
                                               ifelse(TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & sname1=="Blue Morwong (queen snapper)"& livewt<200,13,
                                                      ifelse(TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & sname1=="Blue Morwong (queen snapper)"& livewt>200,26,livewt)))))
  
  Data.daily.2011.12$nfish=with(Data.daily.2011.12,
                                ifelse(DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & sname1=="Shark, Whiskery"& nfish>3,52,
                                       ifelse(DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & sname1=="Shark, Whiskery"& nfish<3,12,nfish)))
  
  Data.daily.2011.12$netlen=with(Data.daily.2011.12,ifelse(vessel=="B 067"& netlen==650,6500,
                                                           ifelse(vessel=="B 142"& netlen==540,5400,
                                                                  ifelse(vessel=="E 009"& netlen==420,4320,
                                                                         ifelse(vessel=="E 009"& netlen==400,4000,
                                                                                ifelse(vessel=="E 035"& netlen==380,3800,
                                                                                       ifelse(vessel=="F 417"& netlen==6000,600,netlen)))))))
  
  Data.daily.2011.12$shots=with(Data.daily.2011.12,ifelse(vessel=="F 541"& shots==5,1,shots))
  
  #create historic file
  Data.daily.2011.12.1=Data.daily.2011.12
  
  #select appropriate variables
  Data.daily.2011.12=Data.daily.2011.12[,match(This,names(Data.daily.2011.12))]
  
  #create demersal suite scalefish species
  Suite=subset(Data.daily.2011.12.1,type=="demersal")     #MISSING: NEED PROPER LIST FROM EVA/DAVE
  Suite=unique(Suite$species)
  
  
  
  #2012-2013
  Data.daily.2012.13$hooks=0
  Data.daily.2012.13$nlines=0
  Data.daily.2012.13$Bioregion=Data.daily.2012.13$bioregion
  Data.daily.2012.13$date=with(Data.daily.2012.13,as.POSIXct(paste(year,"-",month,"-",day,sep="")))
  
  #manually fix some catch records identified as incorrect in "#---- CATCH INSPECTIONS -----" 
  #note: the ammended values were provided by data entry girls
  Data.daily.2012.13$nfish=with (Data.daily.2012.13,
                                 ifelse(species==377004 & vessel=="G 297" & landwt == 92.57,7,nfish))
  
  Data.daily.2012.13$livewt=with (Data.daily.2012.13,
                                  ifelse(species==377004 & vessel=="E 035" & landwt == 49,6.7,
                                         ifelse(species==377004 & vessel=="G 297" & landwt == 92.57,216,
                                                ifelse(species==17003 & vessel=="E 059" & landwt == 1769,17,livewt))))
  
  
  Data.daily.2012.13$netlen=with (Data.daily.2012.13,
                                  ifelse(DSNo=="TDGLF8002239" &  TSNo=="TDGLF8002240" & netlen==165,3500,
                                         ifelse(TSNo=="TDGLF8002261" & netlen==300,3000,
                                                ifelse(TSNo%in%c("TDGLF8008946","TDGLF8008924") & netlen==400,4000,netlen))))
  
  #create historic file
  Data.daily.2012.13.1=Data.daily.2012.13
  
  #select appropriate variables
  Data.daily.2012.13=Data.daily.2012.13[,match(This,names(Data.daily.2012.13))]
  
  #Merge daily data sets
  Data.daily=rbind(Data.daily,Data.daily.2010.11)           
  Data.daily=rbind(Data.daily,Data.daily.2011.12)
  Data.daily=rbind(Data.daily,Data.daily.2012.13)
  
  #Set finyear to character
  Data.daily$finyear=as.character(Data.daily$finyear)
  
}

  #Remove incomplete years  
FINYrs=c(sort(as.character(unique(Data.monthly$FINYEAR))),sort(as.character(unique(Data.daily$finyear))))
idss=match("   .- .",FINYrs)
if(!is.na(idss))FINYrs=FINYrs[-idss]
FINYrs=FINYrs[1:match(Current.yr,FINYrs)]

Data.daily.incomplete=subset(Data.daily,!finyear%in%FINYrs)
Data.daily.incomplete=subset(Data.daily.incomplete,!(species==99999 | is.na(species)))


Data.daily=subset(Data.daily,finyear%in%FINYrs)

#remove fins and livers to avoid duplication 
  #no livers reported, only fins. If conditn=="FI" then 
  # species is corrrect, but if species == 22998, then
  # "WF" & "WH" conditns reported. Inspection of these
  # records showed that the fins correspond to shark weights
  # already in the record and range between 1.1 and 6 kg,
  # hence remove directly

#Before removing records, check if any fisher reported just fins
if(Inspect.New.dat=="YES")
{
  A=subset(Data.daily,species==22998 & finyear==Current.yr)
  THIS=unique(with(A,paste(DSNo,TSNo)))
  NN=length(THIS)
  a=Data.daily
  a$THAT=with(a,(paste(DSNo,TSNo)))
  a=subset(a,THAT%in%THIS)
  guarda=vector('list',NN)
  names(guarda)=THIS[1:NN]
  system.time(for(i in 1:NN)
  {
    b=subset(a,THAT%in%THIS[i])
    checK=unique(b$species)
    checK=checK[-match(22998,checK)]
    if(length(which(checK%in%Shark.species))==0)  guarda[[i]]=THIS[i]
  })
  Check=unlist(guarda)
  if(!is.null(Check))
  {
    stop("records with only fins")
    plot.new()
    text(1,1,paste("DSNo", strsplit(Check, "[ ]")[[1]][1]),col=2)
    mtext('WARNING. Reporting Fins only',3,cex=2,col=2)
  }
  rm(A,a,b,guarda)
}

Data.daily=subset(Data.daily,!(species%in%c(22997,22998)))
if(nrow(Data.daily.incomplete))Data.daily.incomplete=subset(Data.daily.incomplete,!(species%in%c(22997,22998)))
if(nrow(Data.daily.FC.2005_06))Data.daily.FC.2005_06=subset(Data.daily.FC.2005_06,!(species%in%c(22997,22998)))


  #fix typos
Data.daily$netlen=with(Data.daily,
      ifelse(vessel=="G 297" & finyear=="2014-15" & netlen==400,4000,
      ifelse(vessel=="E 075" & finyear=="2014-15" & netlen==200,2000,netlen)))

  #Create same return                                                                          
Data.daily$Same.return=with(Data.daily,paste(finyear,month,vessel,method,blockx)) 

  #Create unique session (i.e. shote)
Data.daily$Same.return.SNo=with(Data.daily,paste(SNo,DSNo,TSNo))
Data.daily$ID=with(Data.daily,paste(DSNo,SNo,"_",sep=""))                                 


  #Create data set identifier
Data.daily$TYPE.DATA="daily"

  #Fix estuaries lats and longs                                       
Data.daily$LatDeg=with(Data.daily,ifelse(blockx%in%c(96021),25,
                  ifelse(blockx%in%c(96022,96023),26,
                  ifelse(blockx%in%c(97011),27,
                  ifelse(blockx%in%c(97012,97013),28,       
                  ifelse(blockx%in%c(97014,97015),29,
                  ifelse(blockx%in%c(96010),33,
                  ifelse(blockx%in%c(96000),33,
                  ifelse(blockx%in%c(96030,95050,95090,95040),35,
                  ifelse(blockx%in%c(85030),34.4586,
                  ifelse(blockx%in%c(85050),34.2853,
                  ifelse(blockx%in%c(85080),33.917,
                  ifelse(blockx%in%c(95010),31.9554,
                  ifelse(blockx%in%c(95020),32.5167,
                  ifelse(blockx%in%c(95030),33.3503,
                  ifelse(blockx%in%c(95060),34.9953,
                  ifelse(blockx%in%c(95070),34.9731,
                  ifelse(blockx%in%c(95080),34.9278,
                  ifelse(blockx%in%c(85990),NA,
                  ifelse(blockx%in%c(85130),34.9278,LatDeg
                  ))))))))))))))))))))
Data.daily$LatDeg=floor(Data.daily$LatDeg)

Data.daily$LongDeg=with(Data.daily,ifelse(blockx%in%c(96021),113,
                   ifelse(blockx%in%c(96022,96023),113,
                   ifelse(blockx%in%c(97011),113,
                   ifelse(blockx%in%c(97012,97013),113,       
                   ifelse(blockx%in%c(97014,97015),113,
                   ifelse(blockx%in%c(96010),115,
                   ifelse(blockx%in%c(96000),115,
                   ifelse(blockx%in%c(96030,95050,95090,95040),118,
                   ifelse(blockx%in%c(85030),118.8897,
                   ifelse(blockx%in%c(85050),119.4819,
                   ifelse(blockx%in%c(85080),120.050,
                   ifelse(blockx%in%c(95010),115.8585,
                   ifelse(blockx%in%c(95020),115.7167,
                   ifelse(blockx%in%c(95030),115.6494,
                   ifelse(blockx%in%c(95060),117.4117,
                   ifelse(blockx%in%c(95070),116.9744,
                   ifelse(blockx%in%c(95080),116.4489,
                   ifelse(blockx%in%c(85990),NA,
                   ifelse(blockx%in%c(85130),116.4489,LongDeg
                   ))))))))))))))))))))
Data.daily$LongDeg=floor(Data.daily$LongDeg)    

  #Idenfity estuaries
Data.daily$Estuary=with(Data.daily,ifelse(blockx%in%Estuaries,"YES","NO"))

# fix blocks                                                                        
Data.daily$blockx=with(Data.daily,ifelse(blockx%in%c(96021),25120,     #Shark Bay
                  ifelse(blockx%in%c(96022,96023),26131,
                  ifelse(blockx%in%c(97011),27132,                        #Abrolhos
                  ifelse(blockx%in%c(97012,97013),28132,       
                  ifelse(blockx%in%c(97014,97015),29132,
                  ifelse(blockx%in%c(96010),33151,                        #Geographe Bay
                  ifelse(blockx%in%c(96000),32150,                        #Cockburn sound
                  ifelse(blockx%in%c(96030),35181,blockx)))))))))         # King George sound

#Fix 0 lat degrees or long degrees
Data.daily$LatDeg=with(Data.daily,ifelse(blockx<35600 & LatDeg==0,as.numeric(substr(blockx,1,2)),LatDeg))  
Data.daily$LongDeg=with(Data.daily,ifelse(blockx<35600 & LongDeg==0,100+as.numeric(substr(blockx,3,4)),LongDeg))

Data.daily$LatDeg=with(Data.daily,ifelse(vessel=="B 067" & blockx%in%c(35610,35630,35640,35650) & LatDeg==0,35,LatDeg))  
Data.daily$LongDeg=with(Data.daily,ifelse(vessel=="B 067" & blockx%in%c(35610,35630,35640,35650) & LongDeg==0,116,LongDeg))

#Add lat and long if NA for calculating zone (allocate to top left corner)
Data.daily$LatDeg=with(Data.daily,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
Data.daily$LatMin=with(Data.daily,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
Data.daily$LongDeg=with(Data.daily,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
Data.daily$LongMin=with(Data.daily,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))

Data.daily$LatMin=with(Data.daily,ifelse(is.na(LatMin) & is.na(block10),NA,LatMin))
Data.daily$LongMin=with(Data.daily,ifelse(is.na(LongMin) & is.na(block10),NA,LongMin))

Data.daily$blockx=Data.daily$blockxFC   #reset blockx to original

  #Fix zones                                                                           
Data.daily$zone=as.character(Data.daily$zone)
Data.daily$zone=as.character(with(Data.daily,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",NA))))))))


if(nrow(Data.daily.FC.2005_06)>0)
{
  Data.daily.FC.2005_06$LatDeg=with(Data.daily.FC.2005_06,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
  Data.daily.FC.2005_06$LatMin=with(Data.daily.FC.2005_06,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
  Data.daily.FC.2005_06$LongDeg=with(Data.daily.FC.2005_06,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
  Data.daily.FC.2005_06$LongMin=with(Data.daily.FC.2005_06,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))
  
  Data.daily.FC.2005_06$zone=as.character(Data.daily.FC.2005_06$zone)
  Data.daily.FC.2005_06$zone=as.character(with(Data.daily.FC.2005_06,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",NA))))))))
}
if(nrow(Data.daily.FC.NA.sp)>0)
{
  Data.daily.FC.NA.sp$LatDeg=with(Data.daily.FC.NA.sp,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
  Data.daily.FC.NA.sp$LatMin=with(Data.daily.FC.NA.sp,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
  Data.daily.FC.NA.sp$LongDeg=with(Data.daily.FC.NA.sp,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
  Data.daily.FC.NA.sp$LongMin=with(Data.daily.FC.NA.sp,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))
  
  Data.daily.FC.NA.sp$LatMin=with(Data.daily.FC.NA.sp,ifelse(is.na(LatMin),0,LatMin))
  Data.daily.FC.NA.sp$LongMin=with(Data.daily.FC.NA.sp,ifelse(is.na(LongMin),0,LongMin))
  
  
  Data.daily.FC.NA.sp$zone=as.character(Data.daily.FC.NA.sp$zone)
  Data.daily.FC.NA.sp$zone=as.character(with(Data.daily.FC.NA.sp,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",NA))))))))
}
if(nrow(Data.daily.incomplete))
{
  Data.daily.incomplete$LatDeg=with(Data.daily.incomplete,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
  Data.daily.incomplete$LatMin=with(Data.daily.incomplete,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
  Data.daily.incomplete$LongDeg=with(Data.daily.incomplete,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
  Data.daily.incomplete$LongMin=with(Data.daily.incomplete,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))
  
  Data.daily.incomplete$zone=as.character(Data.daily.incomplete$zone)
  Data.daily.incomplete$zone=as.character(with(Data.daily.incomplete,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",NA))))))))
}

  #set up data set for Alex Hesp
if(do.Alexs=="YES")
{
  Data.daily.Alex=subset(Data.daily,method=="GN")
  Data.daily.Alex=Data.daily.Alex[!(duplicated(Data.daily.Alex$Same.return.SNo)),]
  Data.daily.1$Same.return.SNo=with(Data.daily.1,paste(SNo,DSNo,TSNo))
  sss=subset(Data.daily.1,Same.return.SNo%in%unique(Data.daily.Alex$Same.return.SNo),select=c(Same.return.SNo,Lat,Long))
  sss=sss[!duplicated(sss$Same.return.SNo),]
  Data.daily.Alex=Data.daily.Alex[,match(c("block10","LongDeg","LongMin","LatDeg","LatMin","date",
                                           "year","month","day","depthMax","nlines","Same.return.SNo","fishery"),names(Data.daily.Alex))]
  Data.daily.Alex=Data.daily.Alex%>%left_join(sss,by=c("Same.return.SNo"))
  #Data.daily.Alex=merge(Data.daily.Alex,sss,by="Same.return.SNo",all.x=T)
  
}

  #set missing block10 to the BLOCKX
Data.daily$LatMinFirst=substr(Data.daily$LatMin,1,1)
Data.daily$LongMinFirst=substr(Data.daily$LongMin,1,1)
Data.daily$LongDegSec=Data.daily$LongDeg-100
Data.daily$block101=with(Data.daily,
      ifelse(is.na(block10) & !is.na(LatDeg),paste(LatDeg,LatMinFirst,LongDegSec,LongMinFirst,sep=""),
      ifelse(is.na(block10) & is.na(LatDeg),as.numeric(paste(substr(blockx,1,2),0,substr(blockx,3,5),sep="")),
             block10)))
Data.daily=Data.daily[,-match(c("LatMinFirst","LongMinFirst","LongDegSec"),names(Data.daily))]

 #Create effort dataset
Data.daily$netlen=with(Data.daily,ifelse(method%in%c("LL","HL","DL")&netlen>0,NA,netlen))
Effort.vars.daily=c("fdays","bdays","hours","hooks","shots","netlen","nlines")
Effort.daily=Data.daily[,match(c("SNo","DSNo","TSNo","Same.return.SNo","finyear","year","month","date",
              "zone","vessel","method","blockx","block10",Effort.vars.daily,"Same.return","ID","LatDeg",
              "LongDeg","LatMin","LongMin"),names(Data.daily))]
Effort.daily$year.c=Effort.daily$year
Effort.daily$LAT=Effort.daily$LatDeg+(Effort.daily$LatMin/60)
Effort.daily$LONG=Effort.daily$LongDeg+(Effort.daily$LongMin/60)
Effort.daily=Effort.daily[,-match(c("LatDeg","LongDeg","LatMin","LongMin"),names(Effort.daily))]

#Export reported depth data
write.csv(subset(Data.daily,species%in%c(17003,17001,18003,18007)),"C:/Matias/Analyses/Catch and effort/Data_outs/Daily.log.depth.csv",row.names=F)


# B.1. Catch Inspections             

#note: 1. run this initially, identify errors, get data entry girls to check, 
#      2.  raise issues to Eva for correcting current catch return before analyses 
#      3.  then fix in "#manually fix some catch records identified as incorrect..."

  #Select data for current year assessment
#note: Current.data is a dummy, only used to inspect current catch and effort.
Yr.current=2000+as.numeric(substr(Current.yr,6,7))
Current.data=subset(Data.daily.1,finyear==Current.yr)

#note: for SoFAR, Data.current.Sofar is only used to get the number of crew and licences
Data.current.Sofar=Current.data
Data.current.Sofar$LAT=-(Data.current.Sofar$LatDeg)  
Data.current.Sofar=Data.current.Sofar[,-match(c("LatMin","LongMin"),names(Data.current.Sofar))]  


if(Inspect.New.dat=="YES")
{
  #create file and path for checking new data
  handle=paste("C:/Matias/Data/Catch and Effort/Check_these_vars/Daily/",Current.yr,sep="")
  if(!file.exists(handle)) dir.create(handle)
  
  Top.mon.ktch=mean(c(15000,20000))
  Top.mon.eff=263.5 # 8500 m X 31 days
  # (i.e. monthly catch >Top.mon.ktch and monthly Km.gn.d > Top.mon.eff)  VIP!!!
  
  #1. ID missing catch weight records in current year data (raise livewt NAs to Eva)  
  NAs=subset(Current.data,is.na(livewt))
  NA.table=table(as.character(NAs$sname1),useNA='ifany')
  
  
  #2. Identify cero catch shots (i.e. NA species and NA weight)
  table(Current.data$NilCatch)  
  NA.shots=subset(Current.data,is.na(as.character(sname1)))
  NA.Sheet.shot=unique(paste(as.character(NA.shots$DSNo),as.character(NA.shots$SNo)))
  Current.data$Sheet.shot=with(Current.data,paste(as.character(DSNo),as.character(SNo)))
  Cero.catch.shot=subset(Current.data,Sheet.shot%in%NA.Sheet.shot)
  table(Cero.catch.shot$species,useNA='ifany')
  
  
  #3. Check good weight calculation
  #note:this compares average weights from returns to the range of possible weights by species.
  #     for some species there's no Wei.range info.
  
    #Check if there's numbers data but no weights      
  check.nfish.weight=subset(Current.data, !is.na(nfish) | !nfish==0)
  check.nfish.weight=subset(check.nfish.weight,is.na(livewt) | livewt==0)
  if(nrow(check.nfish.weight)>0)
  {
    par(bg=2)
    plot.new()
    mtext("there are records ",3,cex=3,col="white")
    mtext("with nfish but no weight",3,-2,cex=3,col="white")
    par(bg="white")
    #stop("records with nfish but no weight")
  }
  write.csv(check.nfish.weight,file=paste(handle,"/Check.no_live.weight.csv",sep=""),row.names=F)
  
  
    #calculate average weight for each species
  Current.data$Avg.wt=with(Current.data,livewt/nfish)
  Uniq.sp=Current.data[,match(c("species","sname1"),names(Current.data))]
  Uniq.sp=Uniq.sp[!(duplicated(Uniq.sp$species)),]
  Uniq.sp=subset(Uniq.sp,!(species==22998))  #remove fins
  Uniq.sp.nam=as.character(Uniq.sp$sname1)
  Uniq.sp=Uniq.sp$species
  names(Uniq.sp)=Uniq.sp.nam
  id=match(Wei.range$SPECIES,Uniq.sp)
  id=subset(id,!is.na(id))
  Uniq.sp.with.weight=Uniq.sp[id]
  Uniq.sp.nam.with.weight=Uniq.sp.nam[id]
  tolerance=0.2   #tolerance bounds for acceptable weight
  fn.avg.wt=function(sp,sp.name)
  {
    Data=subset(Current.data,species==sp & !(is.na(Avg.wt)))
    Weit.r=subset(Wei.range, SPECIES==sp)
    
    if(nrow(Data)>1)
    {
      R=range(Data$Avg.wt)
      R.n=range(Data$nfish)
      R.w=range(Data$livewt)
      
      par(mfcol=c(2,2),mai=c(1.1,.85,.2,.1),oma=c(.1,.1,1,.1))
      hist(Data$nfish,xlab="Numbers per shot",main=sp.name,cex.main=.8,col=2)
      hist(Data$livewt,xlab="Weight (kg) per shot",main="",col=2)
      hist(Data$Avg.wt,xlab="Average weight (kg)",main="",col=2)
      
      Q=quantile(Data$Avg.wt,probs = seq(0, 1, 0.025))
      Q.n=quantile(Data$nfish,probs = seq(0, 1, 0.025))
      
      Rev.weight=0
      if(sum(Data$Avg.wt<(Weit.r$TW.min*(1-tolerance)) | (Data$Avg.wt>Weit.r$TW.max*(1+tolerance)))>0)
      {
        id=which(Data$Avg.wt<(Weit.r$TW.min*(1-tolerance)) | Data$Avg.wt>(Weit.r$TW.max*(1+tolerance)))
        Rev.weight=Data[id,]
        plot(Rev.weight$nfish,Rev.weight$livewt,xlab="nfish",ylab="livewt",main=paste("min wt=",round(Weit.r$TW.min,1),
                                                                                      " max wt=",round(Weit.r$TW.max,1),sep=""),pch=19,col=2,cex.main=.8)
      }
      
      return(list(W.range=R.w,Avg.w.range=R,num.range=R.n,Avg.w.quantiles=Q,num.quantiles=Q.n,Review.Record=Rev.weight))
    }
    if(nrow(Data)<=1)return("PROBLEM")
  }
  Avg.wt.list=vector('list',length=length(Uniq.sp.with.weight))
  names(Avg.wt.list)=Uniq.sp.nam.with.weight
  for (i in 1:length(Uniq.sp.with.weight)) Avg.wt.list[[i]]=fn.avg.wt(Uniq.sp.with.weight[i],Uniq.sp.nam.with.weight[i])
  Current.data=Current.data%>%left_join(Wei.range[,6:8],by=c("species"="SPECIES"))
  # Current.data=merge(Current.data,Wei.range[,6:8],
  #                    by.x="species",by.y="SPECIES",all.x=T)
  
  Current.data$Chk.wt=with(Current.data,
            ifelse(Avg.wt<(TW.min*(1-tolerance)) | 
            Avg.wt>(TW.max*(1+tolerance)),"check","ok"))
  check.weights=subset(Current.data,Chk.wt=="check",
     select=c(TSNo,DSNo,SNo,vessel,finyear,month,
              species,sname1,flagtype,conditn,nfish,landwt,livewt,Avg.wt,TW.min,TW.max))
  idd=match(c("landwt","livewt","Avg.wt","TW.min","TW.max"),names(check.weights))
  check.weights[,idd]=round(check.weights[,idd],2)
  check.weights=check.weights[order(check.weights$finyear,check.weights$month,check.weights$TSNo,
                                    check.weights$DSNo),]
  write.csv(check.weights[,-match(c("Avg.wt","TW.min","TW.max"),names(check.weights))],
            file=paste(handle,"/Check.nfish.weight.typo.csv",sep=""),row.names=F)

  Avg.wt.list.ind=vector('list',length=length(Indicator.species))
  names(Avg.wt.list.ind)=names(Indicator.species)
  for (i in 1:length(Indicator.species)) Avg.wt.list.ind[[i]]=fn.avg.wt(Indicator.species[i],names(Indicator.species)[i])  
  
  # Compare average weight indicator species   
  #note: this is not required, already taken care above
  do.indic="NO"
  if(do.indic=="YES")
  {
    fn.Wght.Freq=function(LAT,SPEC,Gear,MAX.w,BRK,XMAX,YMAX,SPEC.name)
    {
    dat=subset(Current.data,LatDeg>LAT & species==SPEC & method==Gear & finyear==Current.yr)
    dat$Aver.w=dat$livewt/dat$nfish
    yrs=unique(dat$finyear)
    dat=subset(dat,!is.na(Aver.w))
    CHeck=subset(dat,Aver.w>MAX.w)
    CHeck=CHeck[order(CHeck$vessel,CHeck$date),]
    return(CHeck[,match(c("TSNo","DSNo","SNo","date","vessel","species",
                          "sname1","flagtype","conditn","nfish","landwt","livewt","Aver.w"),names(CHeck))])
  }
    Check.Whis.Weight=fn.Wght.Freq(26.5,17003,"GN",MAX.w=max.w.whis,BRK=seq(0,max.w.whis,by=5),XMAX=max.w.whis,YMAX=1300,"Whiskery")
    Check.Gum.Weight=fn.Wght.Freq(26.5,17001,"GN",MAX.w=max.w.gum,BRK=seq(0,max.w.gum,by=5),XMAX=25,YMAX=3200,"Gummy")
    Check.Dus.Weight=fn.Wght.Freq(26.5,18003,"GN",MAX.w=max.w.dus,BRK=seq(0,max.w.dus,by=5),XMAX=50,YMAX=1100,"Dusky")
    Check.San.Weight=fn.Wght.Freq(26.5,18007,"GN",MAX.w=max.w.san,BRK=seq(0,max.w.san,by=5),XMAX=25,YMAX=300,"Sandbar")
    if(nrow(Check.Whis.Weight)>0)write.csv(Check.Whis.Weight,paste(handle,"/Check.Weight.Whis.csv",sep=""),row.names=F)
    if(nrow(Check.Gum.Weight)>0)write.csv(Check.Gum.Weight,paste(handle,"/Check.Weight.Gum.csv",sep=""),row.names=F)
    if(nrow(Check.Dus.Weight)>0)write.csv(Check.Dus.Weight,paste(handle,"/Check.Weight.Dus.csv",sep=""),row.names=F)
    if(nrow(Check.San.Weight)>0)write.csv(Check.San.Weight,paste(handle,"/Check.Weight.San.csv",sep=""),row.names=F)
  }
  
  #spatial distribution of species
  Current.data$Lat=-abs(as.numeric(as.character(Current.data$Lat)))
  Current.data$Long=as.numeric(as.character(Current.data$Long))
  Dist.range=list(Gummy=list(c(113,-29),c(129,-36)), Whiskery=list(c(113,-21),c(129,-36)),
                  Bronzy=list(c(113,-29),c(129,-36)),Dusky=list(c(113,-15),c(129,-36)),
                  Sandbar=list(c(113,-15),c(125,-36)))
  smart.par(n.plots=length(Indicator.species),MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(.1, 0.5, 0))
  for (i in 1:length(Indicator.species))
  {
    s=subset(Current.data,species==Indicator.species[i])
    NM=names(Indicator.species)[i]
    plot(s$Long,s$Lat,main=NM,ylim=c(-36,-10),xlim=c(113,129),ylab="",xlab="")
    pol=Dist.range[[match(NM,names(Dist.range))]]
    polygon(x=c(pol[[1]][1],pol[[2]][1],pol[[2]][1],pol[[1]][1]),
            y=c(pol[[1]][2],pol[[1]][2],pol[[2]][2],pol[[2]][2]),border=2)
  }
  rm(Current.data)
}


  #Create dummies to id NIL catch shots and nil effort
Data.daily$sname1=with(Data.daily,ifelse(is.na(species),"no catch",as.character(sname1))) 
Data.daily$landwt=with(Data.daily,ifelse(is.na(species),0,landwt))
Data.daily$livewt=with(Data.daily,ifelse(is.na(species),0,livewt))
Data.daily$conditn=as.character(Data.daily$conditn)
Data.daily$conditn=with(Data.daily,ifelse(is.na(species),"NIL",conditn))
Data.daily$factor=with(Data.daily,ifelse(is.na(species),0,factor))
Data.daily$species=with(Data.daily,ifelse(is.na(species),999999,species))


  #Create data frame for spatial analysis of 10 min blocks
THIS.10=c("finyear","year","month","day","fishery","licence","zone","orgzone","targetSpecies","block10","blockx",
          "LatDeg","LatMin","LongDeg","LongMin","method","species","sname1","vessel","NilCatch",
          "hours","netlen","shots","depthMax","depthMin","landwt","livewt","totlandwt","bdays")     
if(Get.Daily.Logbook=="NO")
{
  Fine.scale.eff=rbind(Data.daily.1[,match(THIS.10,names(Data.daily.1))],
                       Data.daily.2010.11.1[,match(THIS.10,names(Data.daily.2010.11.1))],
                       Data.daily.2011.12.1[,match(THIS.10,names(Data.daily.2011.12.1))])
}
if(Get.Daily.Logbook=="YES") Fine.scale.eff=Data.daily.1
Fine.scale.eff=subset(Fine.scale.eff,method=="GN")


# B.2. Create some variables 
Data.daily$year.c=Data.daily$year


# B.3. Sort by Year, month and vessel
Data.daily=Data.daily[order(Data.daily$year,Data.daily$month,Data.daily$vessel),]



# #Aggregate number of fish
# Data.daily.agg.Numbers=aggregate(nfish~finyear+year+month+vessel+method+blockx+
#               species+sname1+conditn+factor+year.c+LatDeg+LongDeg+Same.return+TYPE.DATA+Bioregion+zone,
#                                  data=Data.daily,sum,na.rm=T)
# names(Data.daily.agg.Numbers)=c("FINYEAR","YEAR","MONTH","VESSEL","METHOD","BLOCKX","SPECIES",
#                                 "SNAME","CONDITN","Factor","YEAR.c","LAT","LONG","Same.return",
#                                 "TYPE.DATA","Bioregion","zone","nfish")


#Combine degrees and minutes
Data.daily$LatDeg=Data.daily$LatDeg+(Data.daily$LatMin/60)
Data.daily$LongDeg=Data.daily$LongDeg+(Data.daily$LongMin/60)
Data.daily$LatDeg=-Data.daily$LatDeg

#Fix some bioregions
Data.daily$Bioregion.old=Data.daily$Bioregion
# Data.daily$Bioregion=as.character(with(Data.daily,
#   ifelse(LongDeg>=115.5 & LongDeg<=129 & LatDeg<=(-26),"SC", 
#   ifelse(LongDeg<115.5 & LatDeg<=(-27),"WC",
#   ifelse(LongDeg<=114.834 & LatDeg>(-27),"Gascoyne",
#   ifelse(LongDeg>114.834 & LatDeg>(-27) & LongDeg<=129,"NC",NA))))))
# Data.daily$Bioregion=with(Data.daily,
#   ifelse(Bioregion=="SC"& LatDeg>(-34) & LongDeg<115.91 ,"WC",Bioregion))


#Fix weights with no or zero data but with nfish
#note: replace by speciess average for the block or if not available, for the species
Data.daily=Data.daily %>% 
            group_by(species,block10) %>% 
            mutate(livewt= ifelse((is.na(livewt)|livewt==0) & nfish>0, 
                                  mean(livewt, na.rm=TRUE)/mean(nfish, na.rm=TRUE), livewt)) %>%
            group_by(species) %>% 
            mutate(livewt= ifelse((is.na(livewt)|livewt==0) & nfish>0, 
                        mean(livewt, na.rm=TRUE)/mean(nfish, na.rm=TRUE), livewt)) %>%
            as.data.frame()


#Export weight and nfish data from TDGDLF for mean weight analysis
#note: for this analysis, remove records with average weight > maximum weight
#       of species as in some cases nfish is under-reported
Mn.wght.dat=subset(Data.daily,method=="GN" & !(blockx%in%Estuaries) & netlen >100)
Mn.wght.dat$Avrg.w=Mn.wght.dat$livewt/Mn.wght.dat$nfish
Mn.wght.dat$Keep=with(Mn.wght.dat,
      ifelse( Avrg.w>max.w.whis & species==17003,"NO",
      ifelse( Avrg.w>max.w.gum & species==17001,"NO",
      ifelse( Avrg.w>max.w.dus & species==18003,"NO",
      ifelse( Avrg.w>max.w.san & species==18007,"NO","YES")))))                     
Mn.wght.dat=subset(Mn.wght.dat,Keep=="YES")               
Mn.wght.dat=Mn.wght.dat[,-match(c("Avrg.w","Keep"),names(Mn.wght.dat))]
write.csv(Mn.wght.dat,file ="C:/Matias/Analyses/Catch and effort/Logbook.data.mean.weight.csv")

Daily.nfish.agg=aggregate(nfish~finyear+year+month+vessel+method+blockx+blockxFC+species+
                    year.c+Same.return+zone+Estuary,data=Data.daily,sum,na.rm=T)
names(Daily.nfish.agg)=c("FINYEAR","YEAR","MONTH","VESSEL","METHOD",
                         "BLOCKX","blockxFC","SPECIES","YEAR.c","Same.return","zone","Estuary","nfish")


# B.5. sort columns and drop some variables
Data.daily$LAT=-as.numeric(substr(Data.daily$blockx,1,2))  
Data.daily$LONG=100+as.numeric(substr(Data.daily$blockx,3,4))
#Data.daily$LatDeg=Data.daily$LAT
#Data.daily$LongDeg=Data.daily$LONG

Data.daily.depth.max=subset(Data.daily,select=c(Same.return.SNo,depthMax))
Data.daily=Data.daily[,-match(c("depthMax"),names(Data.daily))]

#add variables for FishCUBE
Data.daily$FisheryZone=Data.daily$zone
Data.daily$FisheryCode=Data.daily$fishery
Data.daily$Landing.Port=Data.daily$port

this= c("finyear","year","month","vessel","fdays","method","blockx","bdays","hours","hooks","shots",
        "netlen","species","sname1","livewt","conditn","landwt","factor","year.c","LatDeg","LongDeg",
        "Estuary","Same.return","Same.return.SNo","TYPE.DATA","Bioregion","zone","LatMin","LongMin",
        "day","block10","FisheryZone","FisheryCode","Landing.Port","rowid",FishCube.vars,"nfish")
Data.daily=Data.daily[,match(this,names(Data.daily))]
Data.daily=Data.daily[,sort(names(Data.daily))]
names(Data.daily)=sort(c(names(Data.monthly),"day","block10"))   


# B.6. Aggregate Data.daily by Same.Return (i.e. finyear-month-vessel-method-block)    
#note: Monthly records have ~ 3 times less species reported than Daily logbooks so monthly-aggregated daily logbooks
#       have more rows than Monthly records
Data.daily.agg=aggregate(cbind(LANDWT,LIVEWT)~FINYEAR+YEAR+MONTH+VESSEL+METHOD+BLOCKX+
                          blockxFC+SPECIES+SNAME+CONDITN+Factor+YEAR.c+LAT+LONG+Same.return+
                           TYPE.DATA+Bioregion+zone+Estuary+RSCommonName+RSSpeciesId,
                    data=Data.daily[,-match(c("day","block10"),names(Data.daily))],sum,na.rm=T)


#SECTION C. ---- CATCH MERGING AND CORRECTIONS ----

# C.1. Convert factors to characters to avoid levels issues
fn.see.class=function(A)
{
  X=names(A)
  names(X)=X
  for(i in 1:ncol(A))X[i]=class(A[,i])
  print(sort(X))
}
# fn.see.class(Data.monthly)
# fn.see.class(Data.daily.agg)
# fn.see.class(Data.daily)
# fn.see.class(Effort.monthly)
# fn.see.class(Effort.daily)

Factor.to.charact=c("FINYEAR","VESSEL","METHOD","SNAME")
Effort.Factor.to.charact=Factor.to.charact[1:3]
Factor.to.charact.daily=c("finyear","vessel","method")
for (i in 1:length(Factor.to.charact))
{
  id=match(Factor.to.charact[i],names(Data.monthly))
  Data.monthly[,id]=as.character(Data.monthly[,id])
  Data.daily.agg[,id]=as.character(Data.daily.agg[,id])
  Data.daily[,id]=as.character(Data.daily[,id])
}
for (i in 1:length(Effort.Factor.to.charact))
{
  id=match(Effort.Factor.to.charact[i],names(Effort.monthly))
  Effort.monthly[,id]=as.character(Effort.monthly[,id])  
  id=match(Factor.to.charact.daily[i],names(Effort.daily))
  Effort.daily[,id]=as.character(Effort.daily[,id])
}


#Explore monthly data post 2005-06
Monthly.not.in.daily=subset(Data.monthly,FINYEAR%in%unique(Data.daily$FINYEAR))
TAB=table(Monthly.not.in.daily$FINYEAR,Monthly.not.in.daily$METHOD)

# d=subset(Monthly.not.in.daily,METHOD=="GN")
# plot(d$LONG,d$LAT,main="GN")
# d=subset(Monthly.not.in.daily,!METHOD=="GN")
# plot(d$LONG,d$LAT,main="non-GN")

#Get vessels from TDGDLF to check vessel characteristics
TDGDLF.vessels.month=subset(Data.monthly, METHOD=="GN"& LAT<=(-26) & Estuary=="NO",select=VESSEL)
TDGDLF.vessels.daily=subset(Data.daily, METHOD=="GN"& LAT<=(-26) & Estuary=="NO",select=VESSEL)
TDGDLF.vessels=c(as.character(TDGDLF.vessels.month$VESSEL),
                 as.character(TDGDLF.vessels.daily$VESSEL))
TDGDLF.vessels=data.frame(VESSEL=sort(unique(TDGDLF.vessels)))
write.csv(TDGDLF.vessels,"C:/Matias/Data/Fishing power/TDGDLF.vessels.csv",row.names=F)

#Get skippers
SKIPPERS=subset(Data.daily.1, method=="GN"& LatDeg>=(26) ,select=MastersName)
SKIPERS=data.frame(SKIPPERS=unique(as.character(SKIPPERS$MastersName)))
write.csv(SKIPERS,"C:/Matias/Data/Fishing power/TDGDLF.SKIPERS.csv",row.names=F)


# C.2. Merge Monthly and Aggregated Daily data sets

  #remove any overlapping records
X=substr(2007:(2000+as.numeric(substr(Current.yr,6,7))),3,4)
This=paste(2006:as.numeric(substr(Current.yr,1,4)),X,sep="-")
a=subset(Data.monthly,FINYEAR%in%This)
a$dummy=with(a,paste(FINYEAR, MONTH, BLOCKX, VESSEL,METHOD,SPECIES))
b=subset(Data.daily,FINYEAR%in%This)
b$dummy=with(b,paste(FINYEAR, MONTH, BLOCKX, VESSEL,METHOD,SPECIES))
IIdd=which(b$dummy%in%unique(a$dummy))
if(length(IIdd)>0)
{
  Data.monthly$dummy=with(Data.monthly,paste(FINYEAR, MONTH, BLOCKX, VESSEL,METHOD,SPECIES))
  Data.monthly=subset(Data.monthly,dummy%in%unique(a$dummy))
  Data.monthly=Data.monthly[,-match("dummy",names(Data.monthly))]
}
rm(a,b)


Data.daily.agg$FisheryZone=NA
Data.daily.agg$FisheryCode=NA
Data.daily.agg$Landing.Port=NA
Data.daily.agg$BDAYS=NA
Data.daily.agg$licence=NA
Data.daily.agg$Landing.Port=NA
Data.daily.agg$rowid=NA

  #merge monthly and aggregated daily
Data.monthly=rbind(Data.monthly[match(names(Data.daily.agg),
                names(Data.monthly))],Data.daily.agg)


#check catch
fn.chk.ktch(d1=Data.monthly.original,
            d2=subset(Data.monthly,Same.return
                      %in%unique(Data.monthly.original$Same.return)),
            VAR1="LIVEWT",VAR2="LIVEWT")


# C.3. fix the weight conversion factor values            #Rory's rule 1a (i)       
  #monthly
Data.monthly$Factor.c=with(Data.monthly,ifelse(SPECIES%in%c(22997,22998)|CONDITN%in%c("DR","LI","FI"),0,
                    ifelse(CONDITN=="FL",3.33,ifelse(CONDITN=="WH",1,
                    ifelse(SPECIES==17003,1.5,ifelse(SPECIES%in%c(25000:31000),2.86,
                    ifelse(SPECIES%in%c(23900:23903),2.33,1.59)))))))
Data.monthly$Factor.c=with(Data.monthly,ifelse(is.na(CONDITN) & is.na(Factor.c),Factor,Factor.c))
FACTOR=table(Data.monthly$Factor,useNA='ifany')
FACTOR.c=table(Data.monthly$Factor.c,useNA='ifany')

  #daily
Data.daily$Factor.c=with(Data.daily,ifelse(SPECIES%in%c(22997,22998)|CONDITN%in%c("DR","LI","FI"),0,
                     ifelse(CONDITN=="FL",3.33,ifelse(CONDITN=="WH",1,
                     ifelse(SPECIES==17003,1.5,ifelse(SPECIES%in%c(25000:31000),2.86,
                     ifelse(SPECIES%in%c(23900:23903),2.33,1.59)))))))
Data.daily$Factor.c=with(Data.daily,ifelse(is.na(CONDITN) & is.na(Factor.c),Factor,Factor.c))



    #use specific factors for teleosts
Conditions=subset(Conditions,spgroup=="scalefish")

      #monthly
Data.monthly$Factor.c=with(Data.monthly,ifelse(SPECIES>31000 &!(CONDITN=="WH"),NA,Factor.c))
Data.monthly$Factor.c=with(Data.monthly,ifelse(SPECIES%in%c(337000)&CONDITN=="GG",1.1,
                    ifelse(SPECIES%in%c(445001,311100,441900)&CONDITN=="GG",1.15,
                    ifelse(SPECIES%in%c(361004,441002)&(CONDITN=="GG"|CONDITN=="HD"),1.22,
                    ifelse(SPECIES%in%c(377004,363001,369002,599001,599000,361015)&CONDITN=="GG",1.24,
                    ifelse(SPECIES%in%c(384999,367000,296002,264004,337006,465000,354001,
                    258000,337007,337003,353001,334002,337062,320000)&CONDITN=="GG",1.25,
                    ifelse(SPECIES%in%c(384002)&CONDITN=="GG",1.32,
                    ifelse(SPECIES%in%c(337007)&CONDITN=="HG",1.43,
                    ifelse(SPECIES%in%c(599000)&CONDITN=="HG",1.50,
                    ifelse(SPECIES%in%c(465000)&CONDITN=="HG",2.0,
                    ifelse(SPECIES%in%c(441000)&CONDITN=="GG",1.048,       
                    ifelse(SPECIES%in%c(441020)&CONDITN=="GG",1.11,
                    ifelse(SPECIES%in%c(311903)&CONDITN=="GG",1.12,
                    ifelse(SPECIES%in%c(311006,441004)&CONDITN=="GG",1.15,
                    ifelse(SPECIES%in%c(351005,351008,351009)&CONDITN=="GG",1.16,
                    ifelse(SPECIES%in%c(296000,320001,330001,335001,346004,
                          350000,377000,384904,386000,460000)&CONDITN=="GG",1.25,
                    ifelse(SPECIES%in%c(344002)&CONDITN=="GG",1.3,
                    ifelse(SPECIES%in%c(311000,311012)&CONDITN=="GG",1.32,
                    ifelse(SPECIES%in%c(441007)&CONDITN=="HG",1.176,
                    ifelse(SPECIES%in%c(353001)&CONDITN=="HG",1.39,
                    ifelse(SPECIES%in%c(258004,320000,335001,350000,354001,
                                        441002)&CONDITN=="HG",1.43,
                    ifelse(SPECIES%in%c(599001)&CONDITN=="HG",1.5,
                    ifelse(SPECIES%in%c(384999)&CONDITN=="HG",1.52,
                    ifelse(SPECIES%in%c(361015)&CONDITN=="HG",1.54,
                    ifelse(SPECIES%in%c(296000,377004,384002,386000)&CONDITN=="HG",1.59,
                    ifelse(SPECIES%in%c(351005,351008,351009)&CONDITN=="HG",1.75,
                    ifelse(SPECIES%in%c(353998)&CONDITN=="GG",1.16,
                    ifelse(SPECIES%in%c(311005,361003)&CONDITN=="GG",1.24,
                    ifelse(SPECIES%in%c(367003)&CONDITN=="GG",1.25,
                    ifelse(SPECIES%in%c(367000)&CONDITN=="HG",1.5,
                    ifelse(SPECIES%in%c(465000)&CONDITN=="OT",1,
                     Factor.c)))))))))))))))))))))))))))))))
Data.monthly$Factor.c=with(Data.monthly,ifelse(SPECIES>31000 & is.na(Factor.c),Factor,Factor.c))
# NA.condition=subset(Data.monthly, is.na(Factor.c))  
# NA.condition=NA.condition[,match(c("SPECIES","CONDITN","Factor"),names(NA.condition))]
# table(NA.condition$SPECIES,NA.condition$CONDITN) 


#daily
Data.daily$Factor.c=with(Data.daily,ifelse(SPECIES>31000 &!(CONDITN=="WH"),NA,Factor.c))
Data.daily$Factor.c=with(Data.daily,ifelse(SPECIES%in%c(337000)&CONDITN=="GG",1.1,
         ifelse(SPECIES%in%c(445001,311100,441900)&CONDITN=="GG",1.15,
         ifelse(SPECIES%in%c(361004,441002)&(CONDITN=="GG"|CONDITN=="HD"),1.22,
         ifelse(SPECIES%in%c(377004,363001,369002,599001,599000,361015)&CONDITN=="GG",1.24,
         ifelse(SPECIES%in%c(384999,367000,296002,264004,337006,465000,354001,
         258000,337007,337003,353001,334002,337062,320000)&CONDITN=="GG",1.25,
         ifelse(SPECIES%in%c(384002)&CONDITN=="GG",1.32,
         ifelse(SPECIES%in%c(337007)&CONDITN=="HG",1.43,
         ifelse(SPECIES%in%c(599000)&CONDITN=="HG",1.50,
         ifelse(SPECIES%in%c(465000)&CONDITN=="HG",2.0,
         ifelse(SPECIES%in%c(441000)&CONDITN=="GG",1.048,       
         ifelse(SPECIES%in%c(441020)&CONDITN=="GG",1.11,
         ifelse(SPECIES%in%c(311903)&CONDITN=="GG",1.12,
         ifelse(SPECIES%in%c(311006,441004)&CONDITN=="GG",1.15,
         ifelse(SPECIES%in%c(351005,351008,351009)&CONDITN=="GG",1.16,
         ifelse(SPECIES%in%c(296000,320001,330001,335001,346004,
         350000,377000,384904,386000,460000)&CONDITN=="GG",1.25,
         ifelse(SPECIES%in%c(344002)&CONDITN=="GG",1.3,
         ifelse(SPECIES%in%c(311000,311012)&CONDITN=="GG",1.32,
         ifelse(SPECIES%in%c(441007)&CONDITN=="HG",1.176,
         ifelse(SPECIES%in%c(353001)&CONDITN=="HG",1.39,
         ifelse(SPECIES%in%c(258004,320000,335001,350000,354001,
         441002)&CONDITN=="HG",1.43,
         ifelse(SPECIES%in%c(599001)&CONDITN=="HG",1.5,
         ifelse(SPECIES%in%c(384999)&CONDITN=="HG",1.52,
         ifelse(SPECIES%in%c(361015)&CONDITN=="HG",1.54,
         ifelse(SPECIES%in%c(296000,377004,384002,386000)&CONDITN=="HG",1.59,
         ifelse(SPECIES%in%c(351005,351008,351009)&CONDITN=="HG",1.75,
         ifelse(SPECIES%in%c(353998)&CONDITN=="GG",1.16,
         ifelse(SPECIES%in%c(311005,361003)&CONDITN=="GG",1.24,
         ifelse(SPECIES%in%c(367003)&CONDITN=="GG",1.25,
         ifelse(SPECIES%in%c(367000)&CONDITN=="HG",1.5,
         ifelse(SPECIES%in%c(465000)&CONDITN=="OT",1,
         Factor.c)))))))))))))))))))))))))))))))



    #C.4. Update live weights                   # Rory's rule 1a (ii)              
#note:  maintain livewt if monthly as this was already changed in Table81.d        
#       this is only needed if using CAESS and Logbook data from scratch
Data.monthly$LIVEWT.orgnl=Data.monthly$LIVEWT
Data.daily$LIVEWT.orgnl=Data.daily$LIVEWT

FINYEAR.daily=as.character(unique(Data.daily$FINYEAR))

fn.chk.ktch(d1=Data.monthly.original,
            d2=subset(Data.monthly,Same.return
                      %in%unique(Data.monthly.original$Same.return)),
            VAR1="LIVEWT",VAR2="LIVEWT")


# if(Get.CAESS.Logbook=="YES")
# {
#   Data.monthly$LIVEWT=with(Data.monthly,ifelse(TYPE.DATA=="daily" & SPECIES<=31000 & Factor.c>0,
#                               LANDWT*Factor.c,LIVEWT)) 
#   Data.daily$LIVEWT=with(Data.daily,ifelse(TYPE.DATA=="daily" & SPECIES<=31000 & Factor.c>0,
#                               LANDWT*Factor.c,LIVEWT))
# }
   


    #C.5. Update CAESS data                     #Rory's rule 1a (iii)             
      #monthly
Data.monthly$VesselID=with(Data.monthly,paste(FINYEAR,MONTH,BLOCKX,VESSEL))
Data.monthly$BlockAveID=with(Data.monthly,paste(FINYEAR,BLOCKX,VESSEL))
Data.monthly$AnnualVesselAveID=with(Data.monthly,paste(FINYEAR,VESSEL))
Data.monthly$MonthlyID=with(Data.monthly,paste(FINYEAR,MONTH))
Data.monthly$BlockID=with(Data.monthly,paste(FINYEAR,BLOCKX))
Data.monthly$GoodsplitID=with(Data.monthly,paste(FINYEAR,MONTH,BLOCKX))
Data.monthly$ZoneID=with(Data.monthly,paste(FINYEAR,MONTH,zone))
Data.monthly$ZnID=with(Data.monthly,paste(FINYEAR,zone))

      #daily
Data.daily$VesselID=with(Data.daily,paste(FINYEAR,MONTH,BLOCKX,VESSEL))
Data.daily$BlockAveID=with(Data.daily,paste(FINYEAR,BLOCKX,VESSEL))
Data.daily$AnnualVesselAveID=with(Data.daily,paste(FINYEAR,VESSEL))
Data.daily$MonthlyID=with(Data.daily,paste(FINYEAR,MONTH))
Data.daily$BlockID=with(Data.daily,paste(FINYEAR,BLOCKX))
Data.daily$GoodsplitID=with(Data.daily,paste(FINYEAR,MONTH,BLOCKX))
Data.daily$ZoneID=with(Data.daily,paste(FINYEAR,MONTH,zone))
Data.daily$ZnID=with(Data.daily,paste(FINYEAR,zone))


    #C.6. Define useful variables

#Number of years for monthly records
    #monthly
Data.monthly=Data.monthly[order(Data.monthly$YEAR.c,Data.monthly$MONTH),]
FINYEAR.monthly=as.character(unique(Data.monthly$FINYEAR))
FINYEAR.monthly=FINYEAR.monthly[1:match(Current.yr,FINYEAR.monthly)]
YEAR.c.monthly=sort(unique(Data.monthly$YEAR.c))
NN.monthly=length(FINYEAR.monthly)
Yrs.months=unique(paste(Data.monthly$YEAR.c,Data.monthly$MONTH))
names(Yrs.months)=as.numeric(sapply(strsplit(Yrs.months," "), "[", 1))
N.Yrs.months=length(Yrs.months)

    #daily
Data.daily=Data.daily[order(Data.daily$YEAR.c,Data.daily$MONTH),]
YEAR.c.daily=sort(unique(Data.daily$YEAR.c))
NN.daily=length(FINYEAR.daily)
Yrs.days=unique(paste(Data.daily$YEAR.c,Data.daily$MONTH))
names(Yrs.days)=as.numeric(sapply(strsplit(Yrs.days," "), "[", 1))
N.Yrs.days=length(Yrs.days)


    #C.7. Reapportion species catches
Sp.Prop.Catch=c(17003,17001,18003,18007,17008,20000)
names(Sp.Prop.Catch)=c("Whiskery","Gummy","Dusky","Sandbar","School","Dogfish")
All.species=sort(unique(Data.monthly$SPECIES))
Other.species=All.species[-match(Sp.Prop.Catch,All.species)]
Sharks.other=22999
Reported.Shark.species=Shark.species[-which(Sharks.other==Shark.species)]


#Define species categories
Shark.species=Shark.species[which(Shark.species %in% unique(Data.monthly$SPECIES))]
Gummy=17001;Dusky_whaler=c(18003,18001);Whiskery=17003;Sandbar=18007;Hammerheads=19000;
Spinner=18023;Wobbegongs=13000;Common_saw_shark=23002;School=17008
Ray.species=Ray.species[which(Ray.species %in% unique(Data.monthly$SPECIES))]
Elasmo.species=c(Shark.species,Ray.species)


Scalefish.species=Scalefish.species[which(Scalefish.species %in% unique(Data.monthly$SPECIES))]
Redfishes=c(258000,258004,258005,258006)
Blue_morwong=377004;Blue_groper=384002;West_Australian_dhufish=320000;Pink_snapper=353001;
Boarfishes=367000;Samsonfish=337007;Redfishes=Redfishes;
Mulloway=354001;Sweetlips=350000;Baldchin_groper=384999

Mollusc.species=600000:610000
Mollusc.species=Mollusc.species[which(Mollusc.species %in% unique(Data.monthly$SPECIES))]



#C.7.1 calculate proportions                   #Rory's rule 3b                   
fun.prop=function(DAT,SPEC)
{
  #1. calculate proprtion by same.return  
  this.same.returns=unique(DAT$Same.return)
  dat=subset(DAT,Same.return %in% this.same.returns & SPECIES %in% Elasmo.species)
  dat.species=subset(dat,SPECIES==SPEC)
  
  All.1=aggregate(LIVEWT~Same.return,data=dat,sum,na.rm=T)
  names(All.1)[2]="LIVEWT.all"
  
  Target.sp.1=aggregate(LIVEWT~Same.return,data=dat.species,sum,na.rm=T)
  names(Target.sp.1)[2]="LIVEWT.t"
  
  All=All.1%>%left_join(Target.sp.1,by=c("Same.return"))
  #All=merge(All.1,Target.sp.1,by="Same.return",all.x=T)
  All$Proportion=All$LIVEWT.t/All$LIVEWT.all
  All[is.na(All)]=0
  
  #2. Vessel, gear, fin. year, month, block (given by the "Same.return" variable)
  Prop.VesYrMonBlock=All[,match(c("Same.return","Proportion"),names(All))] 
  
  #3. add needed variables for aggregating
  these.vars=c("Same.return","MonthlyID","GoodsplitID","ZoneID","ZnID","zone")
  these.vars=dat.species[!duplicated(dat.species$Same.return),match(these.vars,names(dat.species))]
  All=All%>%full_join(these.vars,by=c("Same.return"))
  #All=merge(All,these.vars, by="Same.return")
  
  #4. Fin. year, month (mean proportion given by the "MonthlyID" variable)
  Prop.FinYrMon=aggregate(Proportion~MonthlyID,data=All,mean,na.rm=T)
  
  #5. Fin. year, month, block (mean proportion given by the "GoodsplitID" variable)
  Prop.GoodsplitID=aggregate(Proportion~GoodsplitID,data=All,mean,na.rm=T)
  
  #6. Fin. year, month, zone (given by the "ZoneID" variable)  
  Prop.FinYrZone=aggregate(Proportion~ZoneID,data=All,mean,na.rm=T)
  
  #7. Yr, zone (given by the "ZnID" variable)  
  Prop.YrZone=aggregate(Proportion~ZnID,data=All,mean,na.rm=T)
  
  #8. zone (given by the "zone" variable)  
  Prop.Zone=aggregate(Proportion~zone,data=All,mean,na.rm=T)
  
  return(list(Prop.VesYrMonBlock=Prop.VesYrMonBlock,Prop.GoodsplitID=Prop.GoodsplitID,
              Prop.FinYrZone=Prop.FinYrZone,Prop.FinYrMon=Prop.FinYrMon,
              Prop.YrZone=Prop.YrZone,Prop.Zone=Prop.Zone))
}

    #monthly
Catch.prop.whiskery=fun.prop(Data.monthly,17003)
Catch.prop.gummy=fun.prop(Data.monthly,17001)
Catch.prop.dusky=fun.prop(Data.monthly,18003)
Catch.prop.sandbar=fun.prop(Data.monthly,18007)
Catch.prop.school=fun.prop(Data.monthly,17008)
Catch.prop.dogfish=fun.prop(Data.monthly,20000)
Catch.prop.other=fun.prop(Data.monthly,Sharks.other)

    #daily
if(Reapportion.daily=="YES")
{
  Catch.prop.whiskery.daily=fun.prop(Data.daily,17003)
  Catch.prop.gummy.daily=fun.prop(Data.daily,17001)
  Catch.prop.dusky.daily=fun.prop(Data.daily,18003)
  Catch.prop.sandbar.daily=fun.prop(Data.daily,18007)
  Catch.prop.other.daily=fun.prop(Data.daily,Sharks.other)  
}
Catch.prop.school.daily=fun.prop(Data.daily,17008)
Catch.prop.dogfish.daily=fun.prop(Data.daily,20000)



  # Define vessel reporting category

    #first add catch proportions of the different species
#note: by merging to all Same.returns, NAs are introduced when the species was not caught in a 
#     given vessel-month-year-block the NAs are set to 0 as no catch of that species was reported

    #monthly
Vessel.Report=data.frame(Same.return=as.character(unique(Data.monthly$Same.return)))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.gummy$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.gummy$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Gum.Ves.ID"
Vessel.Report$Prop.Gum.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Gum.Ves.ID),0,Prop.Gum.Ves.ID))   
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.whiskery$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.whiskery$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Whi.Ves.ID"
Vessel.Report$Prop.Whi.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Whi.Ves.ID),0,Prop.Whi.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.dusky$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.dusky$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Dus.Ves.ID"
Vessel.Report$Prop.Dus.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Dus.Ves.ID),0,Prop.Dus.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.other$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.other$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Other.Ves.ID"
Vessel.Report$Prop.Other.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Other.Ves.ID),0,Prop.Other.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.school$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.school$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Sch.Ves.ID"
Vessel.Report$Prop.Sch.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Sch.Ves.ID),0,Prop.Sch.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.dogfish$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.dogfish$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.DogS.Ves.ID"
Vessel.Report$Prop.DogS.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.DogS.Ves.ID),0,Prop.DogS.Ves.ID))

    #daily
if(Reapportion.daily=="NO")
{
  Vessel.Report.daily=data.frame(Same.return=as.character(unique(Data.daily$Same.return)))
  Vessel.Report.daily=Vessel.Report.daily%>%left_join(Catch.prop.school.daily$Prop.VesYrMonBlock,by=c("Same.return"))
  #Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.school.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Sch.Ves.ID"
  Vessel.Report.daily$Prop.Sch.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Sch.Ves.ID),0,Prop.Sch.Ves.ID))
  Vessel.Report.daily=Vessel.Report.daily%>%left_join(Catch.prop.dogfish.daily$Prop.VesYrMonBlock,by=c("Same.return"))
  #Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.dogfish.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.DogS.Ves.ID"
  Vessel.Report.daily$Prop.DogS.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.DogS.Ves.ID),0,Prop.DogS.Ves.ID))  
}
if(Reapportion.daily=="YES")
{
  Vessel.Report.daily=data.frame(Same.return=as.character(unique(Data.daily$Same.return)))
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.gummy.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Gum.Ves.ID"
  Vessel.Report.daily$Prop.Gum.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Gum.Ves.ID),0,Prop.Gum.Ves.ID))   
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.whiskery.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Whi.Ves.ID"
  Vessel.Report.daily$Prop.Whi.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Whi.Ves.ID),0,Prop.Whi.Ves.ID))
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.dusky.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Dus.Ves.ID"
  Vessel.Report.daily$Prop.Dus.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Dus.Ves.ID),0,Prop.Dus.Ves.ID))
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.other.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Other.Ves.ID"
  Vessel.Report.daily$Prop.Other.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Other.Ves.ID),0,Prop.Other.Ves.ID))
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.school.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Sch.Ves.ID"
  Vessel.Report.daily$Prop.Sch.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Sch.Ves.ID),0,Prop.Sch.Ves.ID))
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.dogfish.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.DogS.Ves.ID"
  Vessel.Report.daily$Prop.DogS.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.DogS.Ves.ID),0,Prop.DogS.Ves.ID))  
}




#_GILLNET REAPPORTIONING (within TDGDLF)_

  # Merge proportions to data
Data.monthly=Data.monthly%>%left_join(Vessel.Report,by=c("Same.return"))
#Data.monthly=merge(Data.monthly,Vessel.Report,by="Same.return",all.x=T)
Data.daily=Data.daily%>%left_join(Vessel.Report.daily,by=c("Same.return"))
#Data.daily=merge(Data.daily,Vessel.Report.daily,by="Same.return",all.x=T)

  # Separate vessels into 'good' and 'bad' reporters
Data.monthly$Reporter=NA
Data.daily$Reporter=NA
if(Reapportion.daily=="NO")  Data.daily$Reporter="good"


#C.7.2 Northern split                               #Rory's rules 3c      
  #monthly
Data.monthly$Reporter=with(Data.monthly,
  ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & TYPE.DATA=="monthly" & METHOD=="GN" & SPECIES%in%c(22999,Indicator.species)) &
           (Prop.Other.Ves.ID==1|(Prop.Whi.Ves.ID==0 | Prop.Dus.Ves.ID==0)),"bad","good"))

  #daily
if(Reapportion.daily=="YES")
{
  Data.daily$Reporter=with(Data.daily,
       ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & METHOD=="GN" & SPECIES%in%c(22999,Indicator.species)) &
       (Prop.Other.Ves.ID==1|(Prop.Whi.Ves.ID==0 | Prop.Dus.Ves.ID==0)),"bad","good"))  
}


#C.7.3 SW split                                     #Rory's rules 3d                
  #monthly
Data.monthly$Reporter=with(Data.monthly,
  ifelse((LAT<=(-32) & LONG<=125 & TYPE.DATA=="monthly" & 
            METHOD=="GN" & 
            SPECIES%in%c(22999,Indicator.species)) & 
    (Prop.Other.Ves.ID==1|
    (Prop.Whi.Ves.ID==0 & Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
    (Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
    (Prop.Dus.Ves.ID==0 & Prop.Whi.Ves.ID==0)|
    (Prop.Gum.Ves.ID==0) & Prop.Whi.Ves.ID==0 |
    (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Gum.Ves.ID)|
    (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Whi.Ves.ID)|
    (Prop.Gum.Ves.ID>0 & Prop.Gum.Ves.ID==Prop.Whi.Ves.ID)),"bad",Reporter))

  #daily
if(Reapportion.daily=="YES")
{
  Data.daily$Reporter=with(Data.daily,
          ifelse((LAT<=(-32) & LONG<=125 & METHOD=="GN" & SPECIES%in%c(22999,Indicator.species)) & 
               (Prop.Other.Ves.ID==1|
               (Prop.Whi.Ves.ID==0 & Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
               (Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
               (Prop.Dus.Ves.ID==0 & Prop.Whi.Ves.ID==0)|
               (Prop.Gum.Ves.ID==0) & Prop.Whi.Ves.ID==0 |
               (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Gum.Ves.ID)|
               (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Whi.Ves.ID)|
               (Prop.Gum.Ves.ID>0 & Prop.Gum.Ves.ID==Prop.Whi.Ves.ID)),"bad",Reporter))  
}


#C.7.4 SE split                                      #Rory's rules 3e               
  #monthly
Data.monthly$Reporter=with(Data.monthly,
          ifelse((LAT<=(-26) & LONG>125 & 
        TYPE.DATA=="monthly" & METHOD=="GN"  & 
        SPECIES%in%c(22999,Indicator.species)) & 
            (Prop.Other.Ves.ID==1|
            (Prop.Gum.Ves.ID==0 | Prop.Sch.Ves.ID==0)),"bad",Reporter))

  #daily
if(Reapportion.daily=="YES")
{
  Data.daily$Reporter=with(Data.daily,
       ifelse((LAT<=(-26) & LONG>125 & METHOD=="GN"  & SPECIES%in%c(22999,Indicator.species)) & 
             (Prop.Other.Ves.ID==1|(Prop.Gum.Ves.ID==0 | Prop.Sch.Ves.ID==0)),"bad",Reporter))  
}


#set to "good" if no catch reported                                                 
  #monthly
Data.monthly$Reporter=with(Data.monthly,ifelse(SPECIES==999999,"good",Reporter))
Data.monthly$Reporter=with(Data.monthly,ifelse(is.na(Reporter),"good",Reporter))

  #daily
Data.daily$Reporter=with(Data.daily,ifelse(SPECIES==999999,"good",Reporter))
Data.daily$Reporter=with(Data.daily,ifelse(is.na(Reporter),"good",Reporter))


#C.7.5 Update school and dogfish boats                #Rory's rules 3f              
  #monthly
Data.monthly$Sch.or.DogS=with(Data.monthly,ifelse(Prop.Sch.Ves.ID>.1 | Prop.DogS.Ves.ID >0,"Yes","No"))
Data.monthly$Reporter=with(Data.monthly,ifelse(Sch.or.DogS=="Yes" & TYPE.DATA=="monthly","good",Reporter))

  #daily
Data.daily$Sch.or.DogS=with(Data.daily,ifelse(Prop.Sch.Ves.ID>.1 | Prop.DogS.Ves.ID >0,"Yes","No"))
Data.daily$Reporter=with(Data.daily,ifelse(Sch.or.DogS=="Yes","good",Reporter))  


Table.Reporter_gum.dus.whis.sch.dog=table(Data.monthly$FINYEAR,Data.monthly$Reporter,useNA='ifany')



#C.7.6.1 Make good split table: Yr-Mn-Block                      #Rory's rules 4b              
#not: set NA to 0 for those that don't occur
  #recalculate proportions for 'good' reporters only

    #monthly
Catch.prop.whiskery=fun.prop(subset(Data.monthly,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2] &
                                      Reporter=="good" & METHOD=="GN"),17003)
Catch.prop.gummy=fun.prop(subset(Data.monthly,LONG >= Gummy.range[1] & LONG <= Gummy.range[2] & 
                                   LAT <=(-30) & Reporter=="good" & METHOD=="GN"),17001)
Catch.prop.dusky=fun.prop(subset(Data.monthly,LAT <= Dusky.range[1] & LONG <= Dusky.range[2] &
                                   Reporter=="good" & METHOD=="GN"),18003)

Good.split=data.frame(GoodsplitID=as.character(unique(Data.monthly$GoodsplitID)))

Good.split=Good.split%>%left_join(Catch.prop.dusky$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.dusky$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Dus.Good.spl"
Good.split$Prop.Dus.Good.spl=with(Good.split,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))

Good.split=Good.split%>%left_join(Catch.prop.gummy$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.gummy$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Gum.Good.spl"
Good.split$Prop.Gum.Good.spl=with(Good.split,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))

Good.split=Good.split%>%left_join(Catch.prop.whiskery$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.whiskery$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Whi.Good.spl"
Good.split$Prop.Whi.Good.spl=with(Good.split,ifelse(is.na(Prop.Whi.Good.spl),0,Prop.Whi.Good.spl))


    #daily
if(Reapportion.daily=="YES")
{
  Catch.prop.whiskery.daily=fun.prop(subset(Data.daily,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2] &
                                              Reporter=="good" & METHOD=="GN"),17003)
  Catch.prop.gummy.daily=fun.prop(subset(Data.daily,LONG >= Gummy.range[1] & LONG <= Gummy.range[2] & 
                                           LAT <=(-30) & Reporter=="good" & METHOD=="GN"),17001)
  Catch.prop.dusky.daily=fun.prop(subset(Data.daily,LAT <= Dusky.range[1] & LONG <= Dusky.range[2] &
                                           Reporter=="good" & METHOD=="GN"),18003)
  
  Good.split.daily=data.frame(GoodsplitID=as.character(unique(Data.daily$GoodsplitID)))
  
  Good.split.daily=merge(Good.split.daily,Catch.prop.dusky.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Dus.Good.spl"
  Good.split.daily$Prop.Dus.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))
  
  Good.split.daily=merge(Good.split.daily,Catch.prop.gummy.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Gum.Good.spl"
  Good.split.daily$Prop.Gum.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))
  
  Good.split.daily=merge(Good.split.daily,Catch.prop.whiskery.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Whi.Good.spl"
  Good.split.daily$Prop.Whi.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Whi.Good.spl),0,Prop.Whi.Good.spl))  
}


#C.7.6.2 Make good split table: Yr-Mn-Zone                                                       
    #monthly
Zone.good.split=data.frame(ZoneID=as.character(unique(Data.monthly$ZoneID)))

Zone.good.split=Zone.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.dusky$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Dus.Zone.Good.spl"
Zone.good.split$Prop.Dus.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                        0,Prop.Dus.Zone.Good.spl))

Zone.good.split=Zone.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.gummy$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Gum.Zone.Good.spl"
Zone.good.split$Prop.Gum.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                        0,Prop.Gum.Zone.Good.spl))

Zone.good.split=Zone.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.whiskery$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Whi.Zone.Good.spl"
Zone.good.split$Prop.Whi.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Whi.Zone.Good.spl),
                                                                        0,Prop.Whi.Zone.Good.spl))


    #daily
if(Reapportion.daily=="YES")
{
  Zone.good.split.daily=data.frame(ZoneID=as.character(unique(Data.daily$ZoneID)))
  
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.dusky.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Dus.Zone.Good.spl"
  Zone.good.split.daily$Prop.Dus.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                                 0,Prop.Dus.Zone.Good.spl))
  
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.gummy.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Gum.Zone.Good.spl"
  Zone.good.split.daily$Prop.Gum.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                                 0,Prop.Gum.Zone.Good.spl))
  
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.whiskery.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Whi.Zone.Good.spl"
  Zone.good.split.daily$Prop.Whi.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Whi.Zone.Good.spl),
                                                                                 0,Prop.Whi.Zone.Good.spl))  
}



#C.7.6.3 Make good split table: Yr-Mn                          #Rory's rules 4c                 
    #monthly
Monthly.good.split=data.frame(MonthlyID=as.character(unique(Data.monthly$MonthlyID)))

Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.dusky$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Dus.Mon.Good.spl"
Monthly.good.split$Prop.Dus.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                        0,Prop.Dus.Mon.Good.spl))

Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.gummy$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Gum.Mon.Good.spl"
Monthly.good.split$Prop.Gum.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                        0,Prop.Gum.Mon.Good.spl))

Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.whiskery$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Whi.Mon.Good.spl"
Monthly.good.split$Prop.Whi.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Whi.Mon.Good.spl),
                                                                        0,Prop.Whi.Mon.Good.spl))


    #daily
if(Reapportion.daily=="YES")
{
  Monthly.good.split.daily=data.frame(MonthlyID=as.character(unique(Data.daily$MonthlyID)))
  
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.dusky.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Dus.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Dus.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                                      0,Prop.Dus.Mon.Good.spl))
  
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.gummy.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Gum.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Gum.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                                      0,Prop.Gum.Mon.Good.spl))
  
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.whiskery.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Whi.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Whi.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Whi.Mon.Good.spl),
                                                                                      0,Prop.Whi.Mon.Good.spl))
  
}


#C.7.6.4 Make good split table: Yr-Zn                                      
    #monthly
YrZn.good.split=data.frame(ZnID=as.character(unique(Data.monthly$ZnID)))

YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.dusky$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.dusky$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Dus.YrZn.Good.spl"
YrZn.good.split$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                        0,Prop.Dus.YrZn.Good.spl))

YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.gummy$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.gummy$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Gum.YrZn.Good.spl"
YrZn.good.split$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                        0,Prop.Gum.YrZn.Good.spl))

YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.whiskery$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.whiskery$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Whi.YrZn.Good.spl"
YrZn.good.split$Prop.Whi.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Whi.YrZn.Good.spl),
                                                                        0,Prop.Whi.YrZn.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  YrZn.good.split.daily=data.frame(ZnID=as.character(unique(Data.daily$ZnID)))
  
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.dusky.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Dus.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                                 0,Prop.Dus.YrZn.Good.spl))
  
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.gummy.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Gum.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                                 0,Prop.Gum.YrZn.Good.spl))
  
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.whiskery.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Whi.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Whi.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Whi.YrZn.Good.spl),
                                                                                 0,Prop.Whi.YrZn.Good.spl))
  
}



#C.7.7  Update good split catches                       #Rory's rules 4e                          
#If valid average month-year-block proportions available, update "bad" records of Dusky, 
# Gummy and whiskery with this

    #-monthly-
Data.monthly$LIVEWT.reap=with(Data.monthly,ifelse(LAT<=(-26),LIVEWT,NA))

    #add total catch of elasmobranchs for each record       
  
Total.Sk.ktch=subset(Data.monthly,SPECIES%in%Fix.species)  #total catch of reap sp
#Total.Sk.ktch=subset(Data.monthly,SPECIES%in%Elasmo.species) #total elasmo catch
Total.Sk.other.ktch=subset(Data.monthly,SPECIES==Sharks.other)
Total.Sk.ktch=aggregate(LIVEWT~Same.return,data=Total.Sk.ktch,sum,na.rm=T)
names(Total.Sk.ktch)[2]="Tot.shk.livewt"
Total.Sk.other.ktch=aggregate(LIVEWT~Same.return,data=Total.Sk.other.ktch,sum,na.rm=T)
names(Total.Sk.other.ktch)[2]="Shark.other.livewt"
Total.Sk.ktch=Total.Sk.ktch%>%left_join(Total.Sk.other.ktch,by=c("Same.return"))
#Total.Sk.ktch=merge(Total.Sk.ktch,Total.Sk.other.ktch,by="Same.return",all.x=T)
Total.Sk.ktch$Shark.other.livewt=with(Total.Sk.ktch,ifelse(is.na(Shark.other.livewt),0,Shark.other.livewt))

Data.monthly=Data.monthly%>%left_join(Total.Sk.ktch,by=c("Same.return"))%>%
                            left_join(Good.split,by=c("GoodsplitID"))%>%
                            left_join(Zone.good.split,by=c("ZoneID"))%>%
                            left_join(Monthly.good.split,by=c("MonthlyID"))%>%
                            left_join(YrZn.good.split,by=c("ZnID"))
#Data.monthly=merge(Data.monthly,Total.Sk.ktch,by="Same.return",all.x=T)
#Data.monthly=merge(Data.monthly,Good.split,by="GoodsplitID",all.x=T)
#Data.monthly=merge(Data.monthly,Zone.good.split,by="ZoneID",all.x=T)
#Data.monthly=merge(Data.monthly,Monthly.good.split,by="MonthlyID",all.x=T)
#Data.monthly=merge(Data.monthly,YrZn.good.split,by="ZnID",all.x=T)


#Normalise proportions if their sum is >1
  #Good.spl
Data.monthly$dummy=with(Data.monthly,Prop.Dus.Good.spl+Prop.Gum.Good.spl+Prop.Whi.Good.spl)
Data.monthly$Prop.Dus.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.Good.spl/dummy,Prop.Dus.Good.spl))
Data.monthly$Prop.Gum.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.Good.spl/dummy,Prop.Gum.Good.spl))
Data.monthly$Prop.Whi.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.Good.spl/dummy,Prop.Whi.Good.spl))

  #Zone.Good.spl
Data.monthly$dummy=with(Data.monthly,Prop.Dus.Zone.Good.spl+Prop.Gum.Zone.Good.spl+Prop.Whi.Zone.Good.spl)
Data.monthly$Prop.Dus.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.Zone.Good.spl/dummy,Prop.Dus.Zone.Good.spl))
Data.monthly$Prop.Gum.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.Zone.Good.spl/dummy,Prop.Gum.Zone.Good.spl))
Data.monthly$Prop.Whi.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.Zone.Good.spl/dummy,Prop.Whi.Zone.Good.spl))

  #YrZn.Good.spl
Data.monthly$dummy=with(Data.monthly,Prop.Dus.YrZn.Good.spl+Prop.Gum.YrZn.Good.spl+Prop.Whi.YrZn.Good.spl)
Data.monthly$Prop.Dus.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.YrZn.Good.spl/dummy,Prop.Dus.YrZn.Good.spl))
Data.monthly$Prop.Gum.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.YrZn.Good.spl/dummy,Prop.Gum.YrZn.Good.spl))
Data.monthly$Prop.Whi.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.YrZn.Good.spl/dummy,Prop.Whi.YrZn.Good.spl))
Data.monthly=Data.monthly[,-match("dummy",names(Data.monthly))]


  #-daily-
if(Reapportion.daily=="NO") Data.daily$LIVEWT.reap=Data.daily$LIVEWT

if(Reapportion.daily=="YES")
{
  Data.daily$LIVEWT.reap=with(Data.daily,ifelse(LAT<=(-26),LIVEWT,NA))
  
  #add total catch for each record                            
  Total.Sk.ktch.daily=subset(Data.daily,SPECIES%in%Elasmo.species)
  Total.Sk.other.ktch.daily=subset(Data.daily,SPECIES==Sharks.other)
  Total.Sk.ktch.daily=aggregate(LIVEWT~Same.return,data=Total.Sk.ktch.daily,sum,na.rm=T)
  names(Total.Sk.ktch.daily)[2]="Tot.shk.livewt"
  Total.Sk.other.ktch.daily=aggregate(LIVEWT~Same.return,data=Total.Sk.other.ktch.daily,sum,na.rm=T)
  names(Total.Sk.other.ktch.daily)[2]="Shark.other.livewt"
  Total.Sk.ktch.daily=merge(Total.Sk.ktch.daily,Total.Sk.other.ktch.daily,by="Same.return",all.x=T)
  Total.Sk.ktch.daily$Shark.other.livewt=with(Total.Sk.ktch.daily,ifelse(is.na(Shark.other.livewt),0,Shark.other.livewt))
  
  
  Data.daily=merge(Data.daily,Total.Sk.ktch.daily,by="Same.return",all.x=T)
  Data.daily=merge(Data.daily,Good.split.daily,by="GoodsplitID",all.x=T)
  Data.daily=merge(Data.daily,Zone.good.split.daily,by="ZoneID",all.x=T)
  Data.daily=merge(Data.daily,Monthly.good.split.daily,by="MonthlyID",all.x=T)
  Data.daily=merge(Data.daily,YrZn.good.split.daily,by="ZnID",all.x=T)
  
  #remove any duplicates
  #ID=anyDuplicated(Data.daily)
  #Data.daily=Data.daily[-ID,]
  
  #Normalise proportions if their sum is >1
  #Good.spl
  Data.daily$dummy=with(Data.daily,Prop.Dus.Good.spl+Prop.Gum.Good.spl+Prop.Whi.Good.spl)
  Data.daily$Prop.Dus.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.Good.spl/dummy,Prop.Dus.Good.spl))
  Data.daily$Prop.Gum.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.Good.spl/dummy,Prop.Gum.Good.spl))
  Data.daily$Prop.Whi.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.Good.spl/dummy,Prop.Whi.Good.spl))
  
  #Zone.Good.spl
  Data.daily$dummy=with(Data.daily,Prop.Dus.Zone.Good.spl+Prop.Gum.Zone.Good.spl+Prop.Whi.Zone.Good.spl)
  Data.daily$Prop.Dus.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.Zone.Good.spl/dummy,Prop.Dus.Zone.Good.spl))
  Data.daily$Prop.Gum.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.Zone.Good.spl/dummy,Prop.Gum.Zone.Good.spl))
  Data.daily$Prop.Whi.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.Zone.Good.spl/dummy,Prop.Whi.Zone.Good.spl))
  
  #YrZn.Good.spl
  Data.daily$dummy=with(Data.daily,Prop.Dus.YrZn.Good.spl+Prop.Gum.YrZn.Good.spl+Prop.Whi.YrZn.Good.spl)
  Data.daily$Prop.Dus.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.YrZn.Good.spl/dummy,Prop.Dus.YrZn.Good.spl))
  Data.daily$Prop.Gum.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.YrZn.Good.spl/dummy,Prop.Gum.YrZn.Good.spl))
  Data.daily$Prop.Whi.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.YrZn.Good.spl/dummy,Prop.Whi.YrZn.Good.spl))
  Data.daily=Data.daily[,-match("dummy",names(Data.daily))]  
}




#C.7.8 Reaportion the catch of gummy, whiskeries and duskies

  #C.7.8.1 create bad reporter files for fixing catches
    #monthly
Bad.Reporters.All=subset(Data.monthly,Reporter=="bad")
Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter))

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.All.daily=subset(Data.daily,Reporter=="bad")
  Data.daily=subset(Data.daily,Reporter=="good"|is.na(Reporter))  
}


  #monthly
Bad.Reporters=subset(Bad.Reporters.All,SPECIES%in%Fix.species)
Bad.Reporters.otherSP=subset(Bad.Reporters.All,!SPECIES%in%Fix.species)
Data.monthly=rbind(Data.monthly,Bad.Reporters.otherSP) #put back other species

  #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily=subset(Bad.Reporters.All.daily,SPECIES%in%Fix.species)
  Bad.Reporters.otherSP.daily=subset(Bad.Reporters.All.daily,!SPECIES%in%Fix.species)
  Data.daily=rbind(Data.daily,Bad.Reporters.otherSP.daily) #put back other species  
}



#C.7.8.2
#note: remove duplicates of Same.return as Tot.shk.livewt is split proportionally
#          among dusky, whiskery and gummy
#note: this is used to just keep ancilliary information, the total catch from this return is reapportioned

  #monthly
Agg.Sp.Ktch=aggregate(LIVEWT~Same.return,Bad.Reporters,sum)
names(Agg.Sp.Ktch)[2]="Total.LIVEWT.reap"

Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),] 
NroW=nrow(Bad.Reporters)

  #daily
if(Reapportion.daily=="YES")
{
  Agg.Sp.Ktch.daily=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Agg.Sp.Ktch.daily)[2]="Total.LIVEWT.reap"
  Bad.Reporters.daily=Bad.Reporters.daily[!duplicated(Bad.Reporters.daily$Same.return),]
  NroW.daily=nrow(Bad.Reporters.daily)  
}


# replicate Bad.Reporters four times to have the three species and shark others as a record
    #monthly
Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters,Bad.Reporters)
Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]

Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
Bad.Reporters$Sname.old=Bad.Reporters$SNAME

Bad.Reporters$SPECIES=rep(Fix.species,NroW)
Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY","SHARK, OTHER"),NroW)

  #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily=rbind(Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily)
  Bad.Reporters.daily=Bad.Reporters.daily[order(Bad.Reporters.daily$Same.return),]
  
  Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
  Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME
  
  Bad.Reporters.daily$SPECIES=rep(Fix.species,NroW.daily)
  Bad.Reporters.daily$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY","SHARK, OTHER"),NroW.daily)  
}


# reapportion catch
Bad.Reporters$LIVEWT.reap=NA
if(Reapportion.daily=="YES") Bad.Reporters.daily$LIVEWT.reap=NA

#first use "Good.spl" (standardise the proportions to sum(split catch)=Tot.shk.livewt)
  #Monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
    ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
    ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
    ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
    LIVEWT.reap))))

  #Daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
       ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
       ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
       ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
       LIVEWT.reap))))
}

#Second, if previous not available, use "Zone.Good.spl"
  #Monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
    ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
    ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
    ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
    LIVEWT.reap))))

  #Daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
      ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
      ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
      LIVEWT.reap))))  
}


#Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
  #Monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
    ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
    ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
    ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
    LIVEWT.reap))))

  #Daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
    ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
    ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
    ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
    LIVEWT.reap))))  
}


#rescale to avoid creating catch
  #Monthly
Bad.Reporters=Bad.Reporters%>%left_join(Agg.Sp.Ktch,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Agg.Sp.Ktch,by="Same.return",all.x=T)

Agg.reap.KtcH=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)
names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
Bad.Reporters=Bad.Reporters%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Agg.reap.KtcH,by="Same.return",all.x=T)


Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,ifelse(Tot.shk.livewt>Total.LIVEWT.reap,
                  LIVEWT.reap*(Total.LIVEWT.reap/Tot.shk.livewt),LIVEWT.reap))
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,ifelse(LIVEWT.reap.scaler>Total.LIVEWT.reap,
                  LIVEWT.reap*(Total.LIVEWT.reap/LIVEWT.reap.scaler),LIVEWT.reap))

  #Daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Agg.Sp.Ktch.daily,by="Same.return",all.x=T)
  
  Agg.reap.KtcH.daily=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Agg.reap.KtcH.daily)[2]="LIVEWT.reap.scaler"
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Agg.reap.KtcH.daily,by="Same.return",all.x=T)
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,ifelse(Tot.shk.livewt>Total.LIVEWT.reap,
                      LIVEWT.reap*(Total.LIVEWT.reap/Tot.shk.livewt),LIVEWT.reap))
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,ifelse(LIVEWT.reap.scaler>Total.LIVEWT.reap,
                      LIVEWT.reap*(Total.LIVEWT.reap/LIVEWT.reap.scaler),LIVEWT.reap))
  
}


#Reset 'shark other' to the remainder
    #monthly
Remainder=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)
names(Remainder)[2]="Remainder"
Bad.Reporters=Bad.Reporters%>%left_join(Remainder,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Remainder,by="Same.return",all.x=T)
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
                  ifelse(SPECIES==22999,Total.LIVEWT.reap-Remainder,LIVEWT.reap))
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,ifelse(LIVEWT.reap<0,0,LIVEWT.reap))


  #check that reapportioned didn't add catch
b=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)                                                   
A=Agg.Sp.Ktch%>%full_join(b,by=c("Same.return"))
#A=merge(Agg.Sp.Ktch,b,by="Same.return")
plot(A[,2],A[,3])
lines(A[,2],A[,2])
rm(A,b)

Bad.Reporters=Bad.Reporters[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler","Remainder"),
                                    names(Bad.Reporters))]

Bad.Reporters$LIVEWT=with(Bad.Reporters,ifelse(SPECIES==22999,LIVEWT.orgnl,0))
Bad.Reporters$LIVEWT.orgnl=with(Bad.Reporters,ifelse(SPECIES==22999,LIVEWT.orgnl,0))


    #daily
if(Reapportion.daily=="YES") 
{
  Remainder=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Remainder)[2]="Remainder"
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Remainder,by="Same.return",all.x=T)
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
                    ifelse(SPECIES==22999,Total.LIVEWT.reap-Remainder,LIVEWT.reap))
  
  Bad.Reporters.daily=Bad.Reporters.daily[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler"),names(Bad.Reporters.daily))]
  Bad.Reporters.daily$LIVEWT=with(Bad.Reporters.daily,ifelse(SPECIES==22999,LIVEWT.orgnl,0))
  Bad.Reporters.daily$LIVEWT.orgnl=with(Bad.Reporters.daily,ifelse(SPECIES==22999,LIVEWT.orgnl,0))  
}


#create new vars
    #monthly
Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
Bad.Reporters$Reporter="good"
Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
Bad.Reporters$Sname.old=Bad.Reporters$SNAME

    #daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily$Reporter.old=Bad.Reporters.daily$Reporter
  Bad.Reporters.daily$Reporter="good"
  Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
  Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME  
}


#add old species column to data
    #monthly
Data.monthly$Reporter.old="good"
Data.monthly$Spec.old=Data.monthly$SPECIES
Data.monthly$Sname.old=Data.monthly$SNAME

    #daily
if(Reapportion.daily=="YES") 
{
  Data.daily$Reporter.old="good"
  Data.daily$Spec.old=Data.daily$SPECIES
  Data.daily$Sname.old=Data.daily$SNAME  
}

#remove reportioned catch equal to 0
Bad.Reporters=subset(Bad.Reporters,!LIVEWT.reap==0)

#merge reapportioned catch
  #monthly
Bad.Reporters=Bad.Reporters[,match(names(Data.monthly),names(Bad.Reporters))]
Data.monthly=rbind(Data.monthly,Bad.Reporters)  

  #daily
if(Reapportion.daily=="YES") 
{
  Bad.Reporters.daily=Bad.Reporters.daily[,match(names(Data.daily),names(Bad.Reporters.daily))]
  Data.daily=rbind(Data.daily,Bad.Reporters.daily)  
}

d2=subset(Data.monthly,LAT<=(-26) & !is.na(LIVEWT.reap) &TYPE.DATA=="monthly")
fn.chk.ktch(d1=subset(Data.monthly.original,Same.return
                      %in%unique(d2$Same.return)),
            d2=d2,
            VAR1="LIVEWT",VAR2="LIVEWT.reap")
rm(d2)


#C.7.9 Make good sandbar reporters and 4h update good sandbard reporters         #Rory's rules 4g      
#note: identify vessels that report sandbars within sandbar area

  #identify if reporting sandbar by same return
    #monthly
names(Catch.prop.sandbar$Prop.VesYrMonBlock)[2]="Prop.sandbar"
Data.monthly=Data.monthly%>%left_join(Catch.prop.sandbar$Prop.VesYrMonBlock,by=c("Same.return"))
#Data.monthly=merge(Data.monthly,Catch.prop.sandbar$Prop.VesYrMonBlock,by="Same.return",all.x=T)
Data.monthly$Prop.sandbar=with(Data.monthly,ifelse(is.na(Prop.sandbar),0,Prop.sandbar))

    #daily
if(Reapportion.daily=="YES")
{
  names(Catch.prop.sandbar.daily$Prop.VesYrMonBlock)[2]="Prop.sandbar"
  Data.daily=merge(Data.daily,Catch.prop.sandbar.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  Data.daily$Prop.sandbar=with(Data.daily,ifelse(is.na(Prop.sandbar),0,Prop.sandbar))  
}


  #id bad reporters
Data.monthly$SanBar.rep=with(Data.monthly,ifelse(Prop.sandbar==0 & YEAR.c>1984 & TYPE.DATA=="monthly" &
        LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2] & METHOD=='GN',"bad",'good'))
if(Reapportion.daily=="YES") Data.daily$SanBar.rep=with(Data.daily,ifelse(Prop.sandbar==0 & YEAR.c>1984 &TYPE.DATA=="monthly" &
        LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2] & METHOD=='GN',"bad",'good'))

Table.Reporter_sandbar=table(Data.monthly$FINYEAR,Data.monthly$SanBar.rep,useNA='ifany')

  #separate bad sandbar reporters and reapportion only Species ==22999
    #monthly
Sanbar.dat=subset(Data.monthly,SanBar.rep=="bad" & SPECIES== 22999)
Data.monthly=subset(Data.monthly,!(SanBar.rep=="bad" & SPECIES== 22999))
Data.monthly$SanBar.rep="good"

    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.daily=subset(Data.daily,SanBar.rep=="bad" & SPECIES== 22999)
  Data.daily=subset(Data.daily,!(SanBar.rep=="bad" & SPECIES== 22999))
  Data.daily$SanBar.rep="good"  
}


  #find first month-year with sandbar reporting for each vessel
fn.sandbar.vessels=function(dat)
{
  SB.ves=unique(as.character(dat$VESSEL))
  DUMMY=vector('list',length(SB.ves))
  
  for(i in 1:length(SB.ves))
  {
    dummy=subset(dat,VESSEL==SB.ves[i])
    dummy=dummy[order(dummy$YEAR.c,dummy$MONTH),]
    id=which(dummy$SPECIES==18007)[1]
    dummy$SanBar.rep="bad"
    if(!is.na(id)) dummy$SanBar.rep[id:nrow(dummy)]="good"
    DUMMY[[i]]=dummy
  }
  
  return(do.call(rbind,DUMMY))
  
}

  #monthly
Sanbar.dat.c=fn.sandbar.vessels(Sanbar.dat)
Sanbar.dat.good=subset(Sanbar.dat.c,SanBar.rep=="good")
Sanbar.dat.bad=subset(Sanbar.dat.c,SanBar.rep=="bad")


Agg.Sp.Ktch=aggregate(LIVEWT.reap~Same.return,Sanbar.dat.bad,sum)
names(Agg.Sp.Ktch)[2]="Total.LIVEWT.reap"


  #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.c.daily=fn.sandbar.vessels(Sanbar.dat.daily)
  Sanbar.dat.good.daily=subset(Sanbar.dat.c.daily,SanBar.rep=="good")
  Sanbar.dat.bad.daily=subset(Sanbar.dat.c.daily,SanBar.rep=="bad")  
}


#merge good sandbar records back
if(nrow(Sanbar.dat.good)>0)Data.monthly=rbind(Data.monthly,Sanbar.dat.good)
if(Reapportion.daily=="YES")if(nrow(Sanbar.dat.good.daily)>0)Data.daily=rbind(Data.daily,Sanbar.dat.good.daily)


# re-calculate sandbar proportions for good reporters
  #monthly
Catch.prop.sandbar=fun.prop(subset(Data.monthly,SanBar.rep=="good" & METHOD=='GN'& 
          LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2]),18007)

  #daily
if(Reapportion.daily=="YES") Catch.prop.sandbar.daily=fun.prop(subset(Data.daily,SanBar.rep=="good" & METHOD=='GN'& 
          LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2]),18007)


#good split table: Yr-Mn-Block 
    #monthly
Good.split.san=data.frame(GoodsplitID=as.character(unique(Data.monthly$GoodsplitID)))
Good.split.san=Good.split.san%>%left_join(Catch.prop.sandbar$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split.san=merge(Good.split.san,Catch.prop.sandbar$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split.san)[match("Proportion",names(Good.split.san))]="Prop.San.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Good.split.san.daily=data.frame(GoodsplitID=as.character(unique(Data.daily$GoodsplitID)))
  Good.split.san.daily=merge(Good.split.san.daily,Catch.prop.sandbar.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.san.daily)[match("Proportion",names(Good.split.san.daily))]="Prop.San.Good.spl"  
}


#good split table: Yr-Mn-Zone                            
    #monthly
Zone.good.split.san=data.frame(ZoneID=as.character(unique(Data.monthly$ZoneID)))
Zone.good.split.san=Zone.good.split.san%>%left_join(Catch.prop.sandbar$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split.san=merge(Zone.good.split.san,Catch.prop.sandbar$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split.san)[match("Proportion",names(Zone.good.split.san))]="Prop.San.Zone.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Zone.good.split.san.daily=data.frame(ZoneID=as.character(unique(Data.daily$ZoneID)))
  Zone.good.split.san.daily=merge(Zone.good.split.san.daily,Catch.prop.sandbar.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.san.daily)[match("Proportion",names(Zone.good.split.san.daily))]="Prop.San.Zone.Good.spl"  
}


# good split table: Yr-Mn                          
    #monthly
Monthly.good.split.san=data.frame(MonthlyID=as.character(unique(Data.monthly$MonthlyID)))
Monthly.good.split.san=Monthly.good.split.san%>%left_join(Catch.prop.sandbar$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split.san=merge(Monthly.good.split.san,Catch.prop.sandbar$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split.san)[match("Proportion",names(Monthly.good.split.san))]="Prop.San.Mon.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Monthly.good.split.san.daily=data.frame(MonthlyID=as.character(unique(Data.daily$MonthlyID)))
  Monthly.good.split.san.daily=merge(Monthly.good.split.san.daily,Catch.prop.sandbar.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.san.daily)[match("Proportion",names(Monthly.good.split.san.daily))]="Prop.San.Mon.Good.spl"  
}


#good split table: Yr-Zone                            
    #monthly
YrZn.good.split.san=data.frame(ZnID=as.character(unique(Data.monthly$ZnID)))
YrZn.good.split.san=YrZn.good.split.san %>%left_join(Catch.prop.sandbar$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split.san=merge(YrZn.good.split.san,Catch.prop.sandbar$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split.san)[match("Proportion",names(YrZn.good.split.san))]="Prop.San.YrZn.Good.spl"
YrZn.good.split.san$Prop.San.YrZn.Good.spl=with(YrZn.good.split.san,ifelse(is.na(Prop.San.YrZn.Good.spl),
                                                                       0,Prop.San.YrZn.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  YrZn.good.split.san.daily=data.frame(ZnID=as.character(unique(Data.daily$ZnID)))
  YrZn.good.split.san.daily=merge(YrZn.good.split.san.daily,Catch.prop.sandbar.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.san.daily)[match("Proportion",names(YrZn.good.split.san.daily))]="Prop.San.YrZn.Good.spl"
  YrZn.good.split.san.daily$Prop.San.YrZn.Good.spl=with(YrZn.good.split.san.daily,ifelse(is.na(Prop.San.YrZn.Good.spl),
                                                                                         0,Prop.San.YrZn.Good.spl))  
}


#good split table: Zone
    #monthly
Zn.good.split.san=data.frame(zone=as.character(unique(Data.monthly$zone)))
Zn.good.split.san=Zn.good.split.san%>%left_join(Catch.prop.sandbar$Prop.Zone,by=c("zone"))
#Zn.good.split.san=merge(Zn.good.split.san,Catch.prop.sandbar$Prop.Zone,by="zone",all.x=T)
names(Zn.good.split.san)[match("Proportion",names(Zn.good.split.san))]="Prop.San.zone.Good.spl"
Zn.good.split.san$Prop.San.zone.Good.spl=with(Zn.good.split.san,ifelse(is.na(Prop.San.zone.Good.spl),
                                                                           0,Prop.San.zone.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  Zn.good.split.san.daily=data.frame(zone=as.character(unique(Data.daily$zone)))
  Zn.good.split.san.daily=merge(Zn.good.split.san.daily,Catch.prop.sandbar.daily$Prop.Zone,by="zone",all.x=T)
  names(Zn.good.split.san.daily)[match("Proportion",names(Zn.good.split.san.daily))]="Prop.San.zone.Good.spl"
  Zn.good.split.san.daily$Prop.San.zone.Good.spl=with(Zn.good.split.san.daily,ifelse(is.na(Prop.San.zone.Good.spl),
                                                                                     0,Prop.San.zone.Good.spl))  
}

  #monthly merging
Sanbar.dat.bad=Sanbar.dat.bad %>% left_join(Good.split.san,by="GoodsplitID")%>%
                                  left_join(Zone.good.split.san,by="ZoneID")%>%
                                  left_join(Monthly.good.split.san,by="MonthlyID")%>%
                                  left_join(YrZn.good.split.san,by="ZnID")%>%
                                  left_join(Zn.good.split.san,by="zone")
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Good.split.san,by="GoodsplitID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Zone.good.split.san,by="ZoneID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Monthly.good.split.san,by="MonthlyID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,YrZn.good.split.san,by="ZnID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Zn.good.split.san,by="zone",all.x=T)

Sanbar.dat.bad$Prop.San.Good.spl=with(Sanbar.dat.bad,
              ifelse(is.na(Prop.San.Good.spl),0,Prop.San.Good.spl))
Sanbar.dat.bad$Prop.San.Zone.Good.spl=with(Sanbar.dat.bad,
              ifelse(is.na(Prop.San.Zone.Good.spl),0,Prop.San.Zone.Good.spl))
Sanbar.dat.bad$Prop.San.Mon.Good.spl=with(Sanbar.dat.bad,
              ifelse(is.na(Prop.San.Mon.Good.spl),0,Prop.San.Mon.Good.spl))

  #daily merging
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Good.split.san.daily,by="GoodsplitID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Zone.good.split.san.daily,by="ZoneID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Monthly.good.split.san.daily,by="MonthlyID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,YrZn.good.split.san.daily,by="ZnID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Zn.good.split.san.daily,by="zone",all.x=T)
  Sanbar.dat.bad.daily$Prop.San.Good.spl=with(Sanbar.dat.bad.daily,
                                              ifelse(is.na(Prop.San.Good.spl),0,Prop.San.Good.spl))
  Sanbar.dat.bad.daily$Prop.San.Zone.Good.spl=with(Sanbar.dat.bad.daily,
                                                   ifelse(is.na(Prop.San.Zone.Good.spl),0,Prop.San.Zone.Good.spl))
  Sanbar.dat.bad.daily$Prop.San.Mon.Good.spl=with(Sanbar.dat.bad.daily,
                                                  ifelse(is.na(Prop.San.Mon.Good.spl),0,Prop.San.Mon.Good.spl))  
}


# replicate Sanbar.dat.bad twice to have sandbar and shark others as a record
#monthly
NroW.san=nrow(Sanbar.dat.bad)
Sanbar.dat.bad=rbind(Sanbar.dat.bad,Sanbar.dat.bad)   
Sanbar.dat.bad=Sanbar.dat.bad[order(Sanbar.dat.bad$Same.return),]

Sanbar.dat.bad$Spec.old=Sanbar.dat.bad$SPECIES
Sanbar.dat.bad$Sname.old=Sanbar.dat.bad$SNAME
Sanbar.dat.bad$SPECIES=rep(c(18007,22999),NroW.san)
Sanbar.dat.bad$SNAME=rep(c("SHARK, THICKSKIN (SANDBAR)","SHARK, OTHER"),NroW.san)

#daily
if(Reapportion.daily=="YES")
{  
  NroW.san.d=nrow(Sanbar.dat.bad.daily)
  Sanbar.dat.bad.daily=rbind(Sanbar.dat.bad.daily,Sanbar.dat.bad.daily)
  Sanbar.dat.bad.daily$Spec.old=Sanbar.dat.bad.daily$SPECIES
  Sanbar.dat.bad.daily$Sname.old=Sanbar.dat.bad.daily$SNAME
  Sanbar.dat.bad.daily$SPECIES=rep(c(18007,22999),NroW.san.d)
  Sanbar.dat.bad.daily$SNAME=rep(c("SHARK, THICKSKIN (SANDBAR)","SHARK, OTHER"),NroW.san.d)  
}


#C.7.10 Reapportion "sharks, other' in "bad" sandbar      #Rory's rules 4i-4o                   
Sanbar.dat.bad$LIVEWT.reap.san=NA
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=NA

#first use "Good.spl" 
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
   ifelse(SPECIES==18007 & Prop.San.Good.spl>0,LIVEWT.reap*Prop.San.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
   ifelse(SPECIES==18007 & Prop.San.Good.spl>0,LIVEWT.reap*Prop.San.Good.spl,LIVEWT.reap.san))


#Second, if previous not available, use "Zone.Good.spl"
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
   ifelse(SPECIES==18007 & Prop.San.Zone.Good.spl>0,LIVEWT.reap*Prop.San.Zone.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
   ifelse(SPECIES==18007 & Prop.San.Zone.Good.spl>0,LIVEWT.reap*Prop.San.Zone.Good.spl,LIVEWT.reap.san))


#third, if previous are not available, use "YrMonth.Good.spl"
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
   ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.Mon.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
    ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.Mon.Good.spl,LIVEWT.reap.san))


# #fourth, if previous are not available, use "YrZn.Good.spl"
#  Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
#    ifelse(SPECIES==18007 & Prop.San.YrZn.Good.spl>0,LIVEWT.reap*Prop.San.YrZn.Good.spl,LIVEWT.reap.san))
# 
# #finally, use just zone if no other available
#  Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
#     ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.zone.Good.spl ,LIVEWT.reap.san))

#put remainder in "shark other'
    #monthly 
InDx=1:(nrow(Sanbar.dat.bad)-1)
Sanbar.dat.bad$Remainder= with(Sanbar.dat.bad,c(NA,LIVEWT.reap[InDx]-LIVEWT.reap.san[InDx]))
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,ifelse(SPECIES==22999,Remainder,LIVEWT.reap.san))
Sanbar.dat.bad=Sanbar.dat.bad[,-match("Remainder",names(Sanbar.dat.bad))]

    #daily
if(Reapportion.daily=="YES")
{
  InDx=1:(nrow(Sanbar.dat.bad.daily)-1)
  Sanbar.dat.bad.daily$Remainder= with(Sanbar.dat.bad.daily,c(NA,LIVEWT.reap[InDx]-LIVEWT.reap.san[InDx]))
  Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,ifelse(SPECIES==22999,Remainder,LIVEWT.reap.san))
  Sanbar.dat.bad.daily=Sanbar.dat.bad.daily[,-match("Remainder",names(Sanbar.dat.bad.daily))]
  
}


#rescale to equal to total catch or less              
  #monthly
Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.Sp.Ktch,by=c("Same.return"))
#Sanbar.dat.bad=merge(Sanbar.dat.bad,Agg.Sp.Ktch,by="Same.return",all.x=T)
Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad,sum,na.rm=T)
names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
#Sanbar.dat.bad=merge(Sanbar.dat.bad,Agg.reap.KtcH,by="Same.return",all.x=T)

Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,ifelse(LIVEWT.reap.scaler>0,
                    LIVEWT.reap.san*Total.LIVEWT.reap/LIVEWT.reap.scaler,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Agg.Sp.Ktch.daily,by="Same.return",all.x=T)
  Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad.daily,sum,na.rm=T)
  names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Agg.reap.KtcH,by="Same.return",all.x=T)
  
  Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,ifelse(LIVEWT.reap.scaler>0,
           LIVEWT.reap.san*Total.LIVEWT.reap/LIVEWT.reap.scaler,LIVEWT.reap.san))
}
  
    #monthly                 
Sanbar.dat.bad$LIVEWT.reap=Sanbar.dat.bad$LIVEWT.reap.san
Sanbar.dat.bad=Sanbar.dat.bad[,-match(c("Prop.San.Good.spl","Prop.San.Zone.Good.spl","Prop.San.Mon.Good.spl",
          "Prop.San.YrZn.Good.spl","Prop.San.zone.Good.spl","LIVEWT.reap.san"),names(Sanbar.dat.bad))]

    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily$LIVEWT.reap=Sanbar.dat.bad.daily$LIVEWT.reap.san
  Sanbar.dat.bad.daily=Sanbar.dat.bad.daily[,-match(c("Prop.San.Good.spl","Prop.San.Zone.Good.spl","Prop.San.Mon.Good.spl",
                  "Prop.San.YrZn.Good.spl","Prop.San.zone.Good.spl","LIVEWT.reap.san"),names(Sanbar.dat.bad.daily))]  
}


#create file for flagging bad reporters
Sanbar.dat.bad$Reporter.old='bad'
if(Reapportion.daily=="YES")Sanbar.dat.bad.daily$Reporter='bad'

#remove reportioned catch equal to 0
Sanbar.dat.bad=subset(Sanbar.dat.bad,!LIVEWT.reap==0)

#merge corrected sandbar records                        
  #monthly
Sanbar.dat.bad=Sanbar.dat.bad[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler"),names(Sanbar.dat.bad))]
Sanbar.dat.bad$Reporter='good'
Data.monthly=rbind(Data.monthly,Sanbar.dat.bad)
Data.monthly$Reporter.old=with(Data.monthly,ifelse(Reporter=="bad","bad",Reporter.old))

  #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily=Sanbar.dat.bad.daily[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler"),names(Sanbar.dat.bad.daily))]
  Sanbar.dat.bad.daily$Reporter.old=Sanbar.dat.bad.daily$Reporter
  Sanbar.dat.bad.daily$Reporter='good'
  Data.daily=rbind(Data.daily,Sanbar.dat.bad.daily)
  Data.daily$Reporter.old=with(Data.daily,ifelse(Reporter=="bad","bad",Reporter.old))  
}


# if NA livewt reapportioned, then set to livewt (deals with LAT>(-26) cases)
Data.monthly$LIVEWT.reap=with(Data.monthly,ifelse(is.na(LIVEWT.reap),LIVEWT,LIVEWT.reap))
if(Reapportion.daily=="YES")  Data.daily$LIVEWT.reap=with(Data.daily,ifelse(is.na(LIVEWT.reap),LIVEWT,LIVEWT.reap))


#check catch
d2=subset(Data.monthly,LAT<=(-26) & !is.na(LIVEWT.reap) &TYPE.DATA=="monthly")
fn.chk.ktch(d1=subset(Data.monthly.original,Same.return
                      %in%unique(d2$Same.return)),
            d2=d2,
            VAR1="LIVEWT",VAR2="LIVEWT.reap")
rm(d2)



# Export percentage of GN records reapportioned by year
Table.Reporter=table(Data.monthly$Reporter.old)
write.csv(round(100*Table.Reporter[1]/sum(Table.Reporter),1),"C:/Matias/Analyses/Catch and effort/Outputs/Paper/Percent.GN.reapportioned.csv")


#_LONGLINE REAPPORTIONING (within TDGDLF)_

# Merge again original proportions to data
    #monthly
Catch.prop.sandbar=fun.prop(subset(Data.monthly,LAT<=Sandbar.range[1] & LONG <= Sandbar.range[2]),18007)
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.sandbar$Prop.VesYrMonBlock,by=c("Same.return"))
#Vessel.Report=merge(Vessel.Report,Catch.prop.sandbar$Prop.VesYrMonBlock,by="Same.return",all.x=T)
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.San.Ves.ID"
Vessel.Report$Prop.San.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.San.Ves.ID),0,Prop.San.Ves.ID))
Drop.these=names(Vessel.Report)[-match(c("Same.return","Prop.San.Ves.ID"),names(Vessel.Report))]
Data.monthly=Data.monthly[,-match(Drop.these,names(Data.monthly))]
Data.monthly=Data.monthly%>%left_join(Vessel.Report,by=c("Same.return"))
#Data.monthly=merge(Data.monthly,Vessel.Report,by="Same.return",all.x=T)

    #daily
if(Reapportion.daily=="YES")
{
  Catch.prop.sandbar.daily=fun.prop(subset(Data.daily,LAT<=Sandbar.range[1] & LONG <= Sandbar.range[2]),18007)
  Vessel.Report.daily=merge(Vessel.Report.daily,Catch.prop.sandbar.daily$Prop.VesYrMonBlock,by="Same.return",all.x=T)
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.San.Ves.ID"
  Vessel.Report.daily$Prop.San.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.San.Ves.ID),0,Prop.San.Ves.ID))
  Drop.these.daily=names(Vessel.Report.daily)[-match(c("Same.return","Prop.San.Ves.ID"),names(Vessel.Report.daily))]
  Data.daily=Data.daily[,-match(Drop.these.daily,names(Data.daily))]
  Data.daily=merge(Data.daily,Vessel.Report.daily,by="Same.return",all.x=T)  
}


#C.7.12 Northern LL split                                     #Rory's rules 6c            
  #monthly
Data.monthly$Reporter=with(Data.monthly,
      ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & YEAR.c<2007 & TYPE.DATA=="monthly" & METHOD=="LL" &
                SPECIES%in%c(22999,Indicator.species)) &
          (Prop.Other.Ves.ID==1| (Prop.Dus.Ves.ID==0 | Prop.San.Ves.ID==0)),"bad","good"))

  #daily
if(Reapportion.daily=="YES") Data.daily$Reporter=with(Data.daily,
      ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & TYPE.DATA=="monthly" & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) &
          (Prop.Other.Ves.ID==1| (Prop.Dus.Ves.ID==0| Prop.San.Ves.ID==0)),"bad","good"))


#C.7.13  SW LL split                                          #Rory's rules 6d           
  #monthly
Data.monthly$Reporter=with(Data.monthly,
      ifelse((LAT<=(-32) & LONG<125& YEAR.c<2007 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) &
           (Prop.Other.Ves.ID==1|
           (Prop.Whi.Ves.ID==0 & Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
           (Prop.Whi.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
           (Prop.Dus.Ves.ID==0 & Prop.Whi.Ves.ID==0)|
           (Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
           (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Gum.Ves.ID)|
           (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Whi.Ves.ID)|
           (Prop.Gum.Ves.ID>0 & Prop.Gum.Ves.ID==Prop.Whi.Ves.ID)),"bad",Reporter))

  #daily
if(Reapportion.daily=="YES") 
{
  Data.daily$Reporter=with(Data.daily,
                           ifelse((LAT<=(-32) & LONG<125 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) &
                                    (Prop.Other.Ves.ID==1|
                                       (Prop.Whi.Ves.ID==0 & Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
                                       (Prop.Whi.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
                                       (Prop.Dus.Ves.ID==0 & Prop.Whi.Ves.ID==0)|
                                       (Prop.Dus.Ves.ID==0 & Prop.Gum.Ves.ID==0)|
                                       (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Gum.Ves.ID)|
                                       (Prop.Dus.Ves.ID>0 & Prop.Dus.Ves.ID==Prop.Whi.Ves.ID)|
                                       (Prop.Gum.Ves.ID>0 & Prop.Gum.Ves.ID==Prop.Whi.Ves.ID)),"bad",Reporter))  
}


#C.7.14 SE split                                              #Rory's rules 6e         
    #monthly
Data.monthly$Reporter=with(Data.monthly,
       ifelse((LAT<=(-26) & LONG>=125 & YEAR.c<2007 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) & 
          (Prop.Other.Ves.ID==1|(Prop.Gum.Ves.ID==0 & Prop.Sch.Ves.ID==0)),"bad",Reporter))
Table.Reporter.LL=table(Data.monthly$Reporter,useNA='ifany')

    #daily
if(Reapportion.daily=="YES") 
{
  Data.daily$Reporter=with(Data.daily,
                           ifelse((LAT<=(-26) & LONG>=125 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) & 
                                    (Prop.Other.Ves.ID==1|(Prop.Gum.Ves.ID==0 & Prop.Sch.Ves.ID==0)),"bad",Reporter))
  Table.Reporter.LL.daily=table(Data.daily$Reporter,useNA='ifany')
  
}


#C.7.15.1 Make good split table: Yr-Mn-Block                      #Apply Rory's rules 4b
#not: set NA to 0 for those that don't occur
rm(Catch.prop.gummy,Catch.prop.whiskery,Catch.prop.dusky,YrZn.good.split)
if(Reapportion.daily=="YES") rm(Catch.prop.gummy.daily,Catch.prop.whiskery.daily,Catch.prop.dusky.daily,YrZn.good.split.daily)


#recalculate proportions for 'good' reporters only
    #monthly
Catch.prop.whiskery=fun.prop(subset(Data.monthly,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2] &
                                      Reporter=="good" & METHOD=="LL"),17003)
Catch.prop.gummy=fun.prop(subset(Data.monthly,LONG >= Gummy.range[1] & LONG <= Gummy.range[2] & 
                                   LAT <=(-30) & Reporter=="good" & METHOD=="LL"),17001)
Catch.prop.dusky=fun.prop(subset(Data.monthly,LAT <= Dusky.range[1] & LONG <= Dusky.range[2] &
                                   Reporter=="good" & METHOD=="LL"),18003)

    #daily
if(Reapportion.daily=="YES") 
{
  Catch.prop.whiskery.daily=fun.prop(subset(Data.daily,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2] &
                                              Reporter=="good" & METHOD=="LL"),17003)
  Catch.prop.gummy.daily=fun.prop(subset(Data.daily,LONG >= Gummy.range[1] & LONG <= Gummy.range[2] & 
                                           LAT <=(-30) & Reporter=="good" & METHOD=="LL"),17001)
  Catch.prop.dusky.daily=fun.prop(subset(Data.daily,LAT <= Dusky.range[1] & LONG <= Dusky.range[2] &
                                           Reporter=="good" & METHOD=="LL"),18003)  
}


  #monthly
Good.split=data.frame(GoodsplitID=as.character(unique(Data.monthly$GoodsplitID)))
Good.split=Good.split%>%left_join(Catch.prop.dusky$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.dusky$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Dus.Good.spl"
Good.split$Prop.Dus.Good.spl=with(Good.split,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))
Good.split=Good.split%>%left_join(Catch.prop.gummy$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.gummy$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Gum.Good.spl"
Good.split$Prop.Gum.Good.spl=with(Good.split,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))
Good.split=Good.split%>%left_join(Catch.prop.whiskery$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split=merge(Good.split,Catch.prop.whiskery$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Whi.Good.spl"
Good.split$Prop.Whi.Good.spl=with(Good.split,ifelse(is.na(Prop.Whi.Good.spl),0,Prop.Whi.Good.spl))


  #daily
if(Reapportion.daily=="YES") 
{
  Good.split.daily=data.frame(GoodsplitID=as.character(unique(Data.daily$GoodsplitID)))
  Good.split.daily=merge(Good.split.daily,Catch.prop.dusky.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Dus.Good.spl"
  Good.split.daily$Prop.Dus.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))
  Good.split.daily=merge(Good.split.daily,Catch.prop.gummy.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Gum.Good.spl"
  Good.split.daily$Prop.Gum.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))
  Good.split.daily=merge(Good.split.daily,Catch.prop.whiskery.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.daily)[match("Proportion",names(Good.split.daily))]="Prop.Whi.Good.spl"
  Good.split.daily$Prop.Whi.Good.spl=with(Good.split.daily,ifelse(is.na(Prop.Whi.Good.spl),0,Prop.Whi.Good.spl))  
}


#C.7.15.2 Make good split table: Yr-Mn-Zone                            

  #monthly
Zone.good.split=data.frame(ZoneID=as.character(unique(Data.monthly$ZoneID)))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.dusky$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Dus.Zone.Good.spl"
Zone.good.split$Prop.Dus.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                   0,Prop.Dus.Zone.Good.spl))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.gummy$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Gum.Zone.Good.spl"
Zone.good.split$Prop.Gum.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                   0,Prop.Gum.Zone.Good.spl))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split=merge(Zone.good.split,Catch.prop.whiskery$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Whi.Zone.Good.spl"
Zone.good.split$Prop.Whi.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Whi.Zone.Good.spl),         
                                                                   0,Prop.Whi.Zone.Good.spl))

#daily
if(Reapportion.daily=="YES") 
{
  Zone.good.split.daily=data.frame(ZoneID=as.character(unique(Data.daily$ZoneID)))
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.dusky.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Dus.Zone.Good.spl"
  Zone.good.split.daily$Prop.Dus.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                                 0,Prop.Dus.Zone.Good.spl))
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.gummy.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Gum.Zone.Good.spl"
  Zone.good.split.daily$Prop.Gum.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                                 0,Prop.Gum.Zone.Good.spl))
  Zone.good.split.daily=merge(Zone.good.split.daily,Catch.prop.whiskery.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.daily)[match("Proportion",names(Zone.good.split.daily))]="Prop.Whi.Zone.Good.spl"
  Zone.good.split.daily$Prop.Whi.Zone.Good.spl=with(Zone.good.split.daily,ifelse(is.na(Prop.Whi.Zone.Good.spl),
                                                                                 0,Prop.Whi.Zone.Good.spl))
  
}


#C.7.15.3 Make good split table: Yr-Mn                          #Rory's rules 4c  
    #monthly
Monthly.good.split=data.frame(MonthlyID=as.character(unique(Data.monthly$MonthlyID)))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.dusky$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Dus.Mon.Good.spl"
Monthly.good.split$Prop.Dus.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                        0,Prop.Dus.Mon.Good.spl))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.gummy$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Gum.Mon.Good.spl"
Monthly.good.split$Prop.Gum.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                        0,Prop.Gum.Mon.Good.spl))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split=merge(Monthly.good.split,Catch.prop.whiskery$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Whi.Mon.Good.spl"
Monthly.good.split$Prop.Whi.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Whi.Mon.Good.spl),
                                                                        0,Prop.Whi.Mon.Good.spl))

#daily
if(Reapportion.daily=="YES") 
{
  Monthly.good.split.daily=data.frame(MonthlyID=as.character(unique(Data.daily$MonthlyID)))
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.dusky.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Dus.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Dus.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                                      0,Prop.Dus.Mon.Good.spl))
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.gummy.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Gum.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Gum.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                                      0,Prop.Gum.Mon.Good.spl))
  Monthly.good.split.daily=merge(Monthly.good.split.daily,Catch.prop.whiskery.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.daily)[match("Proportion",names(Monthly.good.split.daily))]="Prop.Whi.Mon.Good.spl"
  Monthly.good.split.daily$Prop.Whi.Mon.Good.spl=with(Monthly.good.split.daily,ifelse(is.na(Prop.Whi.Mon.Good.spl),
                                                                                      0,Prop.Whi.Mon.Good.spl))
  
}


#C.7.15.4 Make good split table: Yr-Zn                            
    #monthly
YrZn.good.split=data.frame(ZnID=as.character(unique(Data.monthly$ZnID)))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.dusky$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.dusky$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Dus.YrZn.Good.spl"
YrZn.good.split$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                   0,Prop.Dus.YrZn.Good.spl))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.gummy$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.gummy$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Gum.YrZn.Good.spl"
YrZn.good.split$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                   0,Prop.Gum.YrZn.Good.spl))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.whiskery$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split=merge(YrZn.good.split,Catch.prop.whiskery$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Whi.YrZn.Good.spl"
YrZn.good.split$Prop.Whi.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Whi.YrZn.Good.spl),
                                                                   0,Prop.Whi.YrZn.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  YrZn.good.split.daily=data.frame(ZnID=as.character(unique(Data.daily$ZnID)))
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.dusky.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Dus.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                                 0,Prop.Dus.YrZn.Good.spl))
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.gummy.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Gum.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                                 0,Prop.Gum.YrZn.Good.spl))
  YrZn.good.split.daily=merge(YrZn.good.split.daily,Catch.prop.whiskery.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.daily)[match("Proportion",names(YrZn.good.split.daily))]="Prop.Whi.YrZn.Good.spl"
  YrZn.good.split.daily$Prop.Whi.YrZn.Good.spl=with(YrZn.good.split.daily,ifelse(is.na(Prop.Whi.YrZn.Good.spl),
                                                                                 0,Prop.Whi.YrZn.Good.spl))
  
}



#C.7.15.4  Update good split catches                       #Rory's rules 4e                          
    #monthly
Recalc=c(names(Good.split)[-1],names(Zone.good.split)[-1],names(Monthly.good.split)[-1],
         names(YrZn.good.split)[-1])
Data.monthly=Data.monthly[,-match(Recalc,names(Data.monthly))]

Data.monthly=Data.monthly %>% left_join(Good.split,by="GoodsplitID") %>%
                              left_join(Zone.good.split,by="ZoneID")%>%
                              left_join(Monthly.good.split,by="MonthlyID")%>%
                              left_join(YrZn.good.split,by="ZnID")
# Data.monthly=merge(Data.monthly,Good.split,by="GoodsplitID",all.x=T)
# Data.monthly=merge(Data.monthly,Zone.good.split,by="ZoneID",all.x=T)
# Data.monthly=merge(Data.monthly,Monthly.good.split,by="MonthlyID",all.x=T)
# Data.monthly=merge(Data.monthly,YrZn.good.split,by="ZnID",all.x=T)

  #daily
if(Reapportion.daily=="YES")
{
  Recalc.daily=c(names(Good.split.daily)[-1],names(Zone.good.split.daily)[-1],names(Monthly.good.split.daily)[-1],
                 names(YrZn.good.split.daily)[-1])
  Data.daily=Data.daily[,-match(Recalc.daily,names(Data.daily))]
  Data.daily=merge(Data.daily,Good.split.daily,by="GoodsplitID",all.x=T)
  Data.daily=merge(Data.daily,Zone.good.split.daily,by="ZoneID",all.x=T)
  Data.daily=merge(Data.daily,Monthly.good.split.daily,by="MonthlyID",all.x=T)
  Data.daily=merge(Data.daily,YrZn.good.split.daily,by="ZnID",all.x=T)
  
}


#Normalise proportions if their sum is >1
#Good.spl
    #monthly
Data.monthly$dummy=with(Data.monthly,Prop.Dus.Good.spl+Prop.Gum.Good.spl+Prop.Whi.Good.spl)
Data.monthly$Prop.Dus.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.Good.spl/dummy,Prop.Dus.Good.spl))
Data.monthly$Prop.Gum.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.Good.spl/dummy,Prop.Gum.Good.spl))
Data.monthly$Prop.Whi.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.Good.spl/dummy,Prop.Whi.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  Data.daily$dummy=with(Data.daily,Prop.Dus.Good.spl+Prop.Gum.Good.spl+Prop.Whi.Good.spl)
  Data.daily$Prop.Dus.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.Good.spl/dummy,Prop.Dus.Good.spl))
  Data.daily$Prop.Gum.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.Good.spl/dummy,Prop.Gum.Good.spl))
  Data.daily$Prop.Whi.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.Good.spl/dummy,Prop.Whi.Good.spl))
  
}


#Zone.Good.spl
    #monthly
Data.monthly$dummy=with(Data.monthly,Prop.Dus.Zone.Good.spl+Prop.Gum.Zone.Good.spl+Prop.Whi.Zone.Good.spl)
Data.monthly$Prop.Dus.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.Zone.Good.spl/dummy,Prop.Dus.Zone.Good.spl))
Data.monthly$Prop.Gum.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.Zone.Good.spl/dummy,Prop.Gum.Zone.Good.spl))
Data.monthly$Prop.Whi.Zone.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.Zone.Good.spl/dummy,Prop.Whi.Zone.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  Data.daily$dummy=with(Data.daily,Prop.Dus.Zone.Good.spl+Prop.Gum.Zone.Good.spl+Prop.Whi.Zone.Good.spl)
  Data.daily$Prop.Dus.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.Zone.Good.spl/dummy,Prop.Dus.Zone.Good.spl))
  Data.daily$Prop.Gum.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.Zone.Good.spl/dummy,Prop.Gum.Zone.Good.spl))
  Data.daily$Prop.Whi.Zone.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.Zone.Good.spl/dummy,Prop.Whi.Zone.Good.spl))
  
}


#YrZn.Good.spl
    #monthly
Data.monthly$dummy=with(Data.monthly,Prop.Dus.YrZn.Good.spl+Prop.Gum.YrZn.Good.spl+Prop.Whi.YrZn.Good.spl)
Data.monthly$Prop.Dus.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Dus.YrZn.Good.spl/dummy,Prop.Dus.YrZn.Good.spl))
Data.monthly$Prop.Gum.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Gum.YrZn.Good.spl/dummy,Prop.Gum.YrZn.Good.spl))
Data.monthly$Prop.Whi.YrZn.Good.spl=with(Data.monthly,ifelse(dummy>1,Prop.Whi.YrZn.Good.spl/dummy,Prop.Whi.YrZn.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  Data.daily$dummy=with(Data.daily,Prop.Dus.YrZn.Good.spl+Prop.Gum.YrZn.Good.spl+Prop.Whi.YrZn.Good.spl)
  Data.daily$Prop.Dus.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Dus.YrZn.Good.spl/dummy,Prop.Dus.YrZn.Good.spl))
  Data.daily$Prop.Gum.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Gum.YrZn.Good.spl/dummy,Prop.Gum.YrZn.Good.spl))
  Data.daily$Prop.Whi.YrZn.Good.spl=with(Data.daily,ifelse(dummy>1,Prop.Whi.YrZn.Good.spl/dummy,Prop.Whi.YrZn.Good.spl))
  
}


Data.monthly=Data.monthly[,-match("dummy",names(Data.monthly))]
if(Reapportion.daily=="YES")Data.daily=Data.daily[,-match("dummy",names(Data.daily))]


#C.7.16 Reaportion the catch of gummy, whiskeries and duskies

  #C.7.16.1 create bad reporter files for fixing catches
  #monthly
Bad.Reporters.All=subset(Data.monthly,Reporter=="bad")
Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter))
Bad.Reporters=subset(Bad.Reporters.All,SPECIES%in%Fix.species)
Bad.Reporters.otherSP=subset(Bad.Reporters.All,!SPECIES%in%Fix.species)
Data.monthly=rbind(Data.monthly,Bad.Reporters.otherSP) #put back other species
Data.monthly$Reporter='good'

#daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.All.daily=subset(Data.daily,Reporter=="bad")
  Data.daily=subset(Data.daily,Reporter=="good"|is.na(Reporter))
  Bad.Reporters.daily=subset(Bad.Reporters.All.daily,SPECIES%in%Fix.species)
  Bad.Reporters.otherSP.daily=subset(Bad.Reporters.All.daily,!SPECIES%in%Fix.species)
  Data.daily=rbind(Data.daily,Bad.Reporters.otherSP.daily) #put back other species
  Data.daily$Reporter='good'
  
}


#C.7.16.2
#note: remove duplicates of Same.return as Tot.shk.livewt is split proportionally
#          among dusky, whiskery and gummy
Agg.Sp.Ktch=aggregate(LIVEWT~Same.return,Bad.Reporters,sum)
names(Agg.Sp.Ktch)[2]="Total.LIVEWT.reap"
Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),]
NroW=nrow(Bad.Reporters)

if(Reapportion.daily=="YES")
{
  Agg.Sp.Ktch.daily=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Agg.Sp.Ktch.daily)[2]="Total.LIVEWT.reap"  
  Bad.Reporters.daily=Bad.Reporters.daily[!duplicated(Bad.Reporters.daily$Same.return),]
  NroW.daily=nrow(Bad.Reporters.daily)  
}


# replicate Bad.Reporters three times to have the three species and shark others as a record
    #monthly
Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters,Bad.Reporters)
Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]
Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
Bad.Reporters$Sname.old=Bad.Reporters$SNAME
Bad.Reporters$SPECIES=rep(Fix.species,NroW)
Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY","SHARK, OTHER"),NroW)

  #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily=rbind(Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily)
  Bad.Reporters.daily=Bad.Reporters.daily[order(Bad.Reporters.daily$Same.return),]
  Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
  Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME
  Bad.Reporters.daily$SPECIES=rep(Fix.species,NroW.daily)
  Bad.Reporters.daily$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY","SHARK, OTHER"),NroW.daily)
  
}


# reapportion catch
Bad.Reporters$LIVEWT.reap=NA
if(Reapportion.daily=="YES") Bad.Reporters.daily$LIVEWT.reap=NA

#first use "Good.spl" (standardise the proportions to sum(split catch)=Tot.shk.livewt)
    #monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
      ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
      ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
      ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
      LIVEWT.reap))))

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
                                       ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
                                              ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
                                                     ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
                                                            LIVEWT.reap))))
  
}

#Second, if previous not available, use "Zone.Good.spl"
    #monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
  ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
  ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
  ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
  LIVEWT.reap))))

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
                                       ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
                                              ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
                                                     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
                                                            LIVEWT.reap))))
}



#Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
    #monthly
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
        ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
        ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
        ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
        LIVEWT.reap))))

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
                                       ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
                                              ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
                                                     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
                                                            LIVEWT.reap))))
  
}


#rescale to equal to total catch or less
    #monthly
Bad.Reporters=Bad.Reporters%>%left_join(Agg.Sp.Ktch,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Agg.Sp.Ktch,by="Same.return",all.x=T)
Agg.reap.KtcH=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)
names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
Bad.Reporters=Bad.Reporters%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Agg.reap.KtcH,by="Same.return",all.x=T)

Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
    ifelse(Tot.shk.livewt>Total.LIVEWT.reap,LIVEWT.reap*(Total.LIVEWT.reap/Tot.shk.livewt),LIVEWT.reap))
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
    ifelse(LIVEWT.reap.scaler>Total.LIVEWT.reap,LIVEWT.reap*(Total.LIVEWT.reap/LIVEWT.reap.scaler),LIVEWT.reap))


  #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Agg.Sp.Ktch.daily,by="Same.return",all.x=T)
  Agg.reap.KtcH.daily=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Agg.reap.KtcH.daily)[2]="LIVEWT.reap.scaler"
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Agg.reap.KtcH.daily,by="Same.return",all.x=T)
  
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
      ifelse(Tot.shk.livewt>Total.LIVEWT.reap,LIVEWT.reap*(Total.LIVEWT.reap/Tot.shk.livewt),LIVEWT.reap))
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
      ifelse(LIVEWT.reap.scaler>Total.LIVEWT.reap,LIVEWT.reap*(Total.LIVEWT.reap/LIVEWT.reap.scaler),LIVEWT.reap))
  
}


#Reset 'shark other' to the remainder
    #monthly
Remainder=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)
names(Remainder)[2]="Remainder"
Bad.Reporters=Bad.Reporters%>%left_join(Remainder,by=c("Same.return"))
#Bad.Reporters=merge(Bad.Reporters,Remainder,by="Same.return",all.x=T)
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
                               ifelse(SPECIES==22999,Total.LIVEWT.reap-Remainder,LIVEWT.reap))
Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,ifelse(LIVEWT.reap<0,0,LIVEWT.reap))

#check that reapportioned didn't add catch
b=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)                                                   
A=Agg.Sp.Ktch%>%full_join(b,by=c("Same.return"))
#A=merge(Agg.Sp.Ktch,b,by="Same.return")
plot(A[,2],A[,3])
lines(A[,2],A[,2])
sum(A[,2])==sum(A[,3])
rm(A,b)

Bad.Reporters=Bad.Reporters[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler","Remainder"),
                                    names(Bad.Reporters))]
Bad.Reporters$LIVEWT=with(Bad.Reporters,ifelse(SPECIES==22999,LIVEWT.orgnl,0))
Bad.Reporters$LIVEWT.orgnl=with(Bad.Reporters,ifelse(SPECIES==22999,LIVEWT.orgnl,0))


#daily
if(Reapportion.daily=="YES") 
{
  Remainder=aggregate(LIVEWT.reap~Same.return,Bad.Reporters.daily,sum,na.rm=T)
  names(Remainder)[2]="Remainder"
  Bad.Reporters.daily=merge(Bad.Reporters.daily,Remainder,by="Same.return",all.x=T)
  Bad.Reporters.daily$LIVEWT.reap=with(Bad.Reporters.daily,
                                       ifelse(SPECIES==22999,Total.LIVEWT.reap-Remainder,LIVEWT.reap))
  
  Bad.Reporters.daily=Bad.Reporters.daily[,-match(c("Total.LIVEWT.reap","LIVEWT.reap.scaler"),names(Bad.Reporters.daily))]
  Bad.Reporters.daily$LIVEWT=with(Bad.Reporters.daily,ifelse(SPECIES==22999,LIVEWT.orgnl,0))
  Bad.Reporters.daily$LIVEWT.orgnl=with(Bad.Reporters.daily,ifelse(SPECIES==22999,LIVEWT.orgnl,0))  
}


#create new vars
    #monthly
Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
Bad.ReportCatch.prop.dusky=fun.prop(subset(Data.monthly,LAT <= Dusky.range[1] & LONG <= Dusky.range[2]),18003)
Bad.Reporters$Reporter="good"
Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
Bad.Reporters$Sname.old=Bad.Reporters$SNAME

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily$Reporter.old=Bad.Reporters.daily$Reporter
  Bad.ReportCatch.prop.dusky.daily=fun.prop(subset(Data.daily,LAT <= Dusky.range[1] & LONG <= Dusky.range[2]),18003)
  Bad.Reporters.daily$Reporter="good"
  Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
  Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME
  
}

#remove reportioned catch equal to 0
Bad.Reporters=subset(Bad.Reporters,!LIVEWT.reap==0)

#merge reapportioned catch
    #monthly
Bad.Reporters=Bad.Reporters[,match(names(Data.monthly),names(Bad.Reporters))]
Data.monthly=rbind(Data.monthly,Bad.Reporters)

    #daily
if(Reapportion.daily=="YES")
{
  Bad.Reporters.daily=Bad.Reporters.daily[,match(names(Data.daily),names(Bad.Reporters.daily))]
  Data.daily=rbind(Data.daily,Bad.Reporters.daily) 
}



#REAPPORTION SANDBAR (6k-6s)

#C.7.17 Make good sandbar reporters and 4h update good sandbard reporters         #Rory's rules 4g    
#note: identify vessels that report sandbars within sandbar area

#id bad reporters
    #monthly
Data.monthly$SanBar.rep=NA
Data.monthly$SanBar.rep=with(Data.monthly,ifelse(Prop.sandbar==0 & YEAR.c>1984 & TYPE.DATA=="monthly" &
             LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2] & METHOD=='LL',"bad","good"))

    #daily
if(Reapportion.daily=="YES")
{
  Data.daily$SanBar.rep=NA
  Data.daily$SanBar.rep=with(Data.daily,ifelse(Prop.sandbar==0 & YEAR.c>1984 & TYPE.DATA=="monthly" &
     LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2] & METHOD=='LL',"bad","good"))
  
}


#separate bad and  good sandbar reporters and reapportion only Species ==22999
    #monthly
Sanbar.dat=subset(Data.monthly,SanBar.rep=="bad" & SPECIES== 22999)
Data.monthly=subset(Data.monthly,!(SanBar.rep=="bad" & SPECIES== 22999))
Data.monthly$SanBar.rep="good"

    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.daily=subset(Data.daily,SanBar.rep=="bad" & SPECIES== 22999)
  Data.daily=subset(Data.daily,!(SanBar.rep=="bad" & SPECIES== 22999))
  Data.daily$SanBar.rep="good"  
}


#find first month-year with sandbar reporting for each vessel
fn.sandbar.vessels=function(dat)
{
  SB.ves=unique(as.character(dat$VESSEL))
  DUMMY=vector('list',length(SB.ves))
  
  for(i in 1:length(SB.ves))
  {
    dummy=subset(dat,VESSEL==SB.ves[i])
    dummy=dummy[order(dummy$YEAR.c,dummy$MONTH),]
    id=which(dummy$SPECIES==18007)[1]
    dummy$SanBar.rep="bad"
    if(!is.na(id)) dummy$SanBar.rep[id:nrow(dummy)]="good"
    DUMMY[[i]]=dummy
  }
  
  return(do.call(rbind,DUMMY))
  
}

    #monthly
Sanbar.dat.c=fn.sandbar.vessels(Sanbar.dat)
Sanbar.dat.good=subset(Sanbar.dat.c,SanBar.rep=="good")
Sanbar.dat.bad=subset(Sanbar.dat.c,SanBar.rep=="bad")

Agg.Sp.Ktch=aggregate(LIVEWT.reap~Same.return,Sanbar.dat.bad,sum)
names(Agg.Sp.Ktch)[2]="Total.LIVEWT.reap"


    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.c.daily=fn.sandbar.vessels(Sanbar.dat.daily)
  Sanbar.dat.good.daily=subset(Sanbar.dat.c.daily,SanBar.rep=="good")
  Sanbar.dat.bad.daily=subset(Sanbar.dat.c.daily,SanBar.rep=="bad")
}



#merge good sandbar records back
if(nrow(Sanbar.dat.good)>0)Data.monthly=rbind(Data.monthly,Sanbar.dat.good)
if(Reapportion.daily=="YES") if(nrow(Sanbar.dat.good.daily)>0)Data.daily=rbind(Data.daily,Sanbar.dat.good.daily)


# re-calculate sandbar proportions for good reporters
Catch.prop.sandbar=fun.prop(subset(Data.monthly,SanBar.rep=="good" & METHOD=='LL'&
          LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2]),18007)

if(Reapportion.daily=="YES") Catch.prop.sandbar.daily=fun.prop(subset(Data.daily,SanBar.rep=="good" & METHOD=='LL'&
          LAT<=Sandbar.range[1] & LONG<=Sandbar.range[2]),18007)


#good split table: Yr-Mn-Block 
    #monthly
Good.split.san=data.frame(GoodsplitID=as.character(unique(Data.monthly$GoodsplitID)))
Good.split.san=Good.split.san%>%left_join(Catch.prop.sandbar$Prop.GoodsplitID,by=c("GoodsplitID"))
#Good.split.san=merge(Good.split.san,Catch.prop.sandbar$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
names(Good.split.san)[match("Proportion",names(Good.split.san))]="Prop.San.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Good.split.san.daily=data.frame(GoodsplitID=as.character(unique(Data.daily$GoodsplitID)))
  Good.split.san.daily=merge(Good.split.san.daily,Catch.prop.sandbar.daily$Prop.GoodsplitID,by="GoodsplitID",all.x=T)
  names(Good.split.san.daily)[match("Proportion",names(Good.split.san.daily))]="Prop.San.Good.spl"
  
}


#good split table: Yr-Mn-Zone                            
    #monthly
Zone.good.split.san=data.frame(ZoneID=as.character(unique(Data.monthly$ZoneID)))
Zone.good.split.san=Zone.good.split.san%>%left_join(Catch.prop.sandbar$Prop.FinYrZone,by=c("ZoneID"))
#Zone.good.split.san=merge(Zone.good.split.san,Catch.prop.sandbar$Prop.FinYrZone,by="ZoneID",all.x=T)
names(Zone.good.split.san)[match("Proportion",names(Zone.good.split.san))]="Prop.San.Zone.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Zone.good.split.san.daily=data.frame(ZoneID=as.character(unique(Data.daily$ZoneID)))
  Zone.good.split.san.daily=merge(Zone.good.split.san.daily,Catch.prop.sandbar.daily$Prop.FinYrZone,by="ZoneID",all.x=T)
  names(Zone.good.split.san.daily)[match("Proportion",names(Zone.good.split.san.daily))]="Prop.San.Zone.Good.spl"
  
}


# good split table: Yr-Mn                          
    #monthly
Monthly.good.split.san=data.frame(MonthlyID=as.character(unique(Data.monthly$MonthlyID)))
Monthly.good.split.san=Monthly.good.split.san%>%left_join(Catch.prop.sandbar$Prop.FinYrMon,by=c("MonthlyID"))
#Monthly.good.split.san=merge(Monthly.good.split.san,Catch.prop.sandbar$Prop.FinYrMon,by="MonthlyID",all.x=T)
names(Monthly.good.split.san)[match("Proportion",names(Monthly.good.split.san))]="Prop.San.Mon.Good.spl"

    #daily
if(Reapportion.daily=="YES")
{
  Monthly.good.split.san.daily=data.frame(MonthlyID=as.character(unique(Data.daily$MonthlyID)))
  Monthly.good.split.san.daily=merge(Monthly.good.split.san.daily,Catch.prop.sandbar.daily$Prop.FinYrMon,by="MonthlyID",all.x=T)
  names(Monthly.good.split.san.daily)[match("Proportion",names(Monthly.good.split.san.daily))]="Prop.San.Mon.Good.spl"
  
}


#good split table: Yr-Zone                            
    #monthly
YrZn.good.split.san=data.frame(ZnID=as.character(unique(Data.monthly$ZnID)))
YrZn.good.split.san=YrZn.good.split.san%>%left_join(Catch.prop.sandbar$Prop.YrZone,by=c("ZnID"))
#YrZn.good.split.san=merge(YrZn.good.split.san,Catch.prop.sandbar$Prop.YrZone,by="ZnID",all.x=T)
names(YrZn.good.split.san)[match("Proportion",names(YrZn.good.split.san))]="Prop.San.YrZn.Good.spl"
YrZn.good.split.san$Prop.San.YrZn.Good.spl=with(YrZn.good.split.san,ifelse(is.na(Prop.San.YrZn.Good.spl),
                                                                           0,Prop.San.YrZn.Good.spl))
    #daily
if(Reapportion.daily=="YES")
{
  YrZn.good.split.san.daily=data.frame(ZnID=as.character(unique(Data.daily$ZnID)))
  YrZn.good.split.san.daily=merge(YrZn.good.split.san.daily,Catch.prop.sandbar.daily$Prop.YrZone,by="ZnID",all.x=T)
  names(YrZn.good.split.san.daily)[match("Proportion",names(YrZn.good.split.san.daily))]="Prop.San.YrZn.Good.spl"
  YrZn.good.split.san.daily$Prop.San.YrZn.Good.spl=with(YrZn.good.split.san.daily,ifelse(is.na(Prop.San.YrZn.Good.spl),
                                                                                         0,Prop.San.YrZn.Good.spl))
  
}


#good split table: Zone
    #monthly
Zn.good.split.san=data.frame(zone=as.character(unique(Data.monthly$zone)))
Zn.good.split.san=Zn.good.split.san%>%left_join(Catch.prop.sandbar$Prop.Zone,by=c("zone"))
#Zn.good.split.san=merge(Zn.good.split.san,Catch.prop.sandbar$Prop.Zone,by="zone",all.x=T)
names(Zn.good.split.san)[match("Proportion",names(Zn.good.split.san))]="Prop.San.zone.Good.spl"
Zn.good.split.san$Prop.San.zone.Good.spl=with(Zn.good.split.san,ifelse(is.na(Prop.San.zone.Good.spl),
                                                                       0,Prop.San.zone.Good.spl))

    #daily
if(Reapportion.daily=="YES")
{
  Zn.good.split.san.daily=data.frame(zone=as.character(unique(Data.daily$zone)))
  Zn.good.split.san.daily=merge(Zn.good.split.san.daily,Catch.prop.sandbar.daily$Prop.Zone,by="zone",all.x=T)
  names(Zn.good.split.san.daily)[match("Proportion",names(Zn.good.split.san.daily))]="Prop.San.zone.Good.spl"
  Zn.good.split.san.daily$Prop.San.zone.Good.spl=with(Zn.good.split.san.daily,ifelse(is.na(Prop.San.zone.Good.spl),
                                                                                     0,Prop.San.zone.Good.spl))
  
}

    #monthly 
Sanbar.dat.bad=Sanbar.dat.bad%>% left_join(Good.split.san,by="GoodsplitID")%>%
                                 left_join(Zone.good.split.san,by="ZoneID")%>%
                                 left_join(Monthly.good.split.san,by="MonthlyID")%>%
                                 left_join(YrZn.good.split.san,by="ZnID")%>%
                                 left_join(Zn.good.split.san,by="zone")

# Sanbar.dat.bad=merge(Sanbar.dat.bad,Good.split.san,by="GoodsplitID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Zone.good.split.san,by="ZoneID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Monthly.good.split.san,by="MonthlyID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,YrZn.good.split.san,by="ZnID",all.x=T)
# Sanbar.dat.bad=merge(Sanbar.dat.bad,Zn.good.split.san,by="zone",all.x=T)

Sanbar.dat.bad$Prop.San.Good.spl=with(Sanbar.dat.bad,
             ifelse(is.na(Prop.San.Good.spl),0,Prop.San.Good.spl))
Sanbar.dat.bad$Prop.San.Zone.Good.spl=with(Sanbar.dat.bad,
             ifelse(is.na(Prop.San.Zone.Good.spl),0,Prop.San.Zone.Good.spl))
Sanbar.dat.bad$Prop.San.Mon.Good.spl=with(Sanbar.dat.bad,
            ifelse(is.na(Prop.San.Mon.Good.spl),0,Prop.San.Mon.Good.spl))


    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Good.split.san.daily,by="GoodsplitID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Zone.good.split.san.daily,by="ZoneID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Monthly.good.split.san.daily,by="MonthlyID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,YrZn.good.split.san.daily,by="ZnID",all.x=T)
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Zn.good.split.san.daily,by="zone",all.x=T)
  Sanbar.dat.bad.daily$Prop.San.Good.spl=with(Sanbar.dat.bad.daily,
                                              ifelse(is.na(Prop.San.Good.spl),0,Prop.San.Good.spl))
  Sanbar.dat.bad.daily$Prop.San.Zone.Good.spl=with(Sanbar.dat.bad.daily,
                                                   ifelse(is.na(Prop.San.Zone.Good.spl),0,Prop.San.Zone.Good.spl))
  Sanbar.dat.bad.daily$Prop.San.Mon.Good.spl=with(Sanbar.dat.bad.daily,
                                                  ifelse(is.na(Prop.San.Mon.Good.spl),0,Prop.San.Mon.Good.spl))
  
}


# replicate Sanbar.dat.bad twice to have sandbar and shark others as a record
#monthly
NroW.san=nrow(Sanbar.dat.bad)
Sanbar.dat.bad=rbind(Sanbar.dat.bad,Sanbar.dat.bad)
Sanbar.dat.bad=Sanbar.dat.bad[order(Sanbar.dat.bad$Same.return),]

Sanbar.dat.bad$Spec.old=Sanbar.dat.bad$SPECIES
Sanbar.dat.bad$Sname.old=Sanbar.dat.bad$SNAME
Sanbar.dat.bad$SPECIES=rep(c(18007,22999),NroW.san)
Sanbar.dat.bad$SNAME=rep(c("SHARK, THICKSKIN (SANDBAR)","SHARK, OTHER"),NroW.san)

#daily
if(Reapportion.daily=="YES")
{  
  NroW.san.d=nrow(Sanbar.dat.bad.daily)
  Sanbar.dat.bad.daily=rbind(Sanbar.dat.bad.daily,Sanbar.dat.bad.daily)
  Sanbar.dat.bad.daily$Spec.old=Sanbar.dat.bad.daily$SPECIES
  Sanbar.dat.bad.daily$Sname.old=Sanbar.dat.bad.daily$SNAME
  Sanbar.dat.bad.daily$SPECIES=rep(c(18007,22999),NroW.san.d)
  Sanbar.dat.bad.daily$SNAME=rep(c("SHARK, THICKSKIN (SANDBAR)","SHARK, OTHER"),NroW.san.d)  
}


#C.7.17.1 Reapportion "sharks, other' in "bad" sandbar      #Rory's rules 4i-4o                     
Sanbar.dat.bad$LIVEWT.reap.san=NA
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=NA


#first use "Good.spl" 
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
    ifelse(SPECIES==18007 & Prop.San.Good.spl>0,LIVEWT.reap*Prop.San.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
    ifelse(SPECIES==18007 & Prop.San.Good.spl>0,LIVEWT.reap*Prop.San.Good.spl,LIVEWT.reap.san))


#Second, if previous not available, use "Zone.Good.spl"
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
    ifelse(SPECIES==18007 & Prop.San.Zone.Good.spl>0,LIVEWT.reap*Prop.San.Zone.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
    ifelse(SPECIES==18007 & Prop.San.Zone.Good.spl>0,LIVEWT.reap*Prop.San.Zone.Good.spl,LIVEWT.reap.san))


#third, if previous are not available, use "YrMonth.Good.spl"
  #monthly
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
    ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.Mon.Good.spl,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES") Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,
    ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.Mon.Good.spl,LIVEWT.reap.san))


# #fourth, if previous are not available, use "YrZn.Good.spl"

#  Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
#    ifelse(SPECIES==18007 & Prop.San.YrZn.Good.spl>0,LIVEWT.reap*Prop.San.YrZn.Good.spl,LIVEWT.reap.san))
# 
# #finally, use just zone if no other available
#  Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,
#     ifelse(SPECIES==18007,LIVEWT.reap*Prop.San.zone.Good.spl ,LIVEWT.reap.san))

#put remainder in "shark other'
#monthly 
InDx=1:(nrow(Sanbar.dat.bad)-1)
Sanbar.dat.bad$Remainder= with(Sanbar.dat.bad,c(NA,LIVEWT.reap[InDx]-LIVEWT.reap.san[InDx]))
Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,ifelse(SPECIES==22999,Remainder,LIVEWT.reap.san))
Sanbar.dat.bad=Sanbar.dat.bad[,-match("Remainder",names(Sanbar.dat.bad))]

#daily
if(Reapportion.daily=="YES")
{
  InDx=1:(nrow(Sanbar.dat.bad.daily)-1)
  Sanbar.dat.bad.daily$Remainder= with(Sanbar.dat.bad.daily,c(NA,LIVEWT.reap[InDx]-LIVEWT.reap.san[InDx]))
  Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,ifelse(SPECIES==22999,Remainder,LIVEWT.reap.san))
  Sanbar.dat.bad.daily=Sanbar.dat.bad.daily[,-match("Remainder",names(Sanbar.dat.bad.daily))]
  
}


#rescale to equal to total catch or less              
  #monthly
Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.Sp.Ktch,by=c("Same.return"))
#Sanbar.dat.bad=merge(Sanbar.dat.bad,Agg.Sp.Ktch,by="Same.return",all.x=T)
Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad,sum,na.rm=T)
names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
#Sanbar.dat.bad=merge(Sanbar.dat.bad,Agg.reap.KtcH,by="Same.return",all.x=T)

Sanbar.dat.bad$LIVEWT.reap.san=with(Sanbar.dat.bad,ifelse(LIVEWT.reap.scaler>0,
                        LIVEWT.reap.san*Total.LIVEWT.reap/LIVEWT.reap.scaler,LIVEWT.reap.san))

  #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Agg.Sp.Ktch.daily,by="Same.return",all.x=T)
  Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad.daily,sum,na.rm=T)
  names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
  Sanbar.dat.bad.daily=merge(Sanbar.dat.bad.daily,Agg.reap.KtcH,by="Same.return",all.x=T)
  
  Sanbar.dat.bad.daily$LIVEWT.reap.san=with(Sanbar.dat.bad.daily,ifelse(LIVEWT.reap.scaler>0,
                                                                        LIVEWT.reap.san*Total.LIVEWT.reap/LIVEWT.reap.scaler,LIVEWT.reap.san))
}


#monthly
Sanbar.dat.bad$LIVEWT.reap=Sanbar.dat.bad$LIVEWT.reap.san

Sanbar.dat.bad=Sanbar.dat.bad[,-match(c("Prop.San.Good.spl","Prop.San.Zone.Good.spl","Prop.San.Mon.Good.spl",
            "Prop.San.YrZn.Good.spl","Prop.San.zone.Good.spl","LIVEWT.reap.san",
            "Total.LIVEWT.reap" , "LIVEWT.reap.scaler"),names(Sanbar.dat.bad))]

#daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily$LIVEWT.reap=Sanbar.dat.bad.daily$LIVEWT.reap.san
  Sanbar.dat.bad.daily=Sanbar.dat.bad.daily[,-match(c("Prop.San.Good.spl","Prop.San.Zone.Good.spl","Prop.San.Mon.Good.spl",
              "Prop.San.YrZn.Good.spl","Prop.San.zone.Good.spl","LIVEWT.reap.san"),names(Sanbar.dat.bad.daily))]
  
}



#create file for flagging bad reporters
Sanbar.dat.bad$Reporter='bad'
if(Reapportion.daily=="YES")Sanbar.dat.bad.daily$Reporter='bad'

#remove reportioned catch equal to 0
Sanbar.dat.bad=subset(Sanbar.dat.bad,!LIVEWT.reap==0)


#merge corrected sandbar records
    #monthly
Sanbar.dat.bad$Reporter.old=Sanbar.dat.bad$Reporter
Sanbar.dat.bad$Reporter='good'
Data.monthly=rbind(Data.monthly,Sanbar.dat.bad)

    #daily
if(Reapportion.daily=="YES")
{
  Sanbar.dat.bad.daily$Reporter.old=Sanbar.dat.bad.daily$Reporter
  Sanbar.dat.bad.daily$Reporter='good'
  Data.daily=rbind(Data.daily,Sanbar.dat.bad.daily)
  
}

#Check catch
d2=subset(Data.monthly,LAT<=(-26) & !is.na(LIVEWT.reap) &TYPE.DATA=="monthly")
fn.chk.ktch(d1=subset(Data.monthly.original,Same.return
                      %in%unique(d2$Same.return)),d2=d2,
            VAR1="LIVEWT",VAR2="LIVEWT.reap")



# C.7.19 Reclass catches as estuary of combined mullet, whiting and herring if catch is >10%
    #monthly
whiting.herring.mullet=subset(Data.monthly,METHOD=="GN"& SPECIES%in%c(344001,330001))
whiting.herring.mullet=subset(Data.monthly,Same.return%in%whiting.herring.mullet$Same.return)
whit.her.mul=aggregate(LIVEWT.reap~Same.return,subset(whiting.herring.mullet,SPECIES%in%c(344001,330001)),sum)
All.whit.her.mul=aggregate(LIVEWT.reap~Same.return,whiting.herring.mullet,sum)
whit.her.mul=whit.her.mul%>%left_join(All.whit.her.mul,by=c("Same.return"))
#whit.her.mul=merge(whit.her.mul,All.whit.her.mul,by="Same.return")
whit.her.mul$Prop.mul.her.whi=with(whit.her.mul,LIVEWT.reap.x/LIVEWT.reap.y)
whit.her.mul=whit.her.mul[,c(1,4)]
Data.monthly=Data.monthly%>%left_join(whit.her.mul,by=c("Same.return"))
#Data.monthly=merge(Data.monthly,whit.her.mul,by="Same.return",all.x=T)
Data.monthly$Prop.mul.her.whi=with(Data.monthly,ifelse(is.na(Prop.mul.her.whi),0,Prop.mul.her.whi))
Data.monthly$Estuary=with(Data.monthly,ifelse(Prop.mul.her.whi>0.1,"YES",Estuary))
Data.monthly=Data.monthly[,-match("Prop.mul.her.whi",names(Data.monthly))]

    #daily
if(Reapportion.daily=="YES")
{
  whiting.herring.mullet.daily=subset(Data.daily,METHOD=="GN"& SPECIES%in%c(344001,330001))
  whiting.herring.mullet.daily=subset(Data.daily,Same.return%in%whiting.herring.mullet.daily$Same.return)
  whit.her.mul.daily=aggregate(LIVEWT.reap~Same.return,subset(whiting.herring.mullet.daily,SPECIES%in%c(344001,330001)),sum)
  All.whit.her.mul.daily=aggregate(LIVEWT.reap~Same.return,whiting.herring.mullet.daily,sum)
  whit.her.mul.daily=merge(whit.her.mul.daily,All.whit.her.mul.daily,by="Same.return")
  whit.her.mul.daily$Prop.mul.her.whi=with(whit.her.mul.daily,LIVEWT.reap.x/LIVEWT.reap.y)
  whit.her.mul.daily=whit.her.mul.daily[,c(1,4)]
  Data.daily=merge(Data.daily,whit.her.mul.daily,by="Same.return",all.x=T)
  Data.daily$Prop.mul.her.whi=with(Data.daily,ifelse(is.na(Prop.mul.her.whi),0,Prop.mul.her.whi))
  Data.daily$Estuary=with(Data.daily,ifelse(Prop.mul.her.whi>0.1,"YES",Estuary))
  Data.daily=Data.daily[,-match("Prop.mul.her.whi",names(Data.daily))]
  
}



# C.7.20 REAPPORTION CATCH FROM NON-TDGDLF CATCHES
#note: this deals with speciess==22999 for non TDGDFL
if(reapportion.ktch.other.method=="YES")
{
  #monthly
  Data.monthly.other=subset(Data.monthly,SPECIES==22999 & !METHOD%in%c("GN","LL")& Bioregion%in%c("SC","WC"))
  Data.monthly=subset(Data.monthly,!(SPECIES==22999 & !METHOD%in%c("GN","LL") & Bioregion%in%c("SC","WC")))
  
  #daily
  if(Reapportion.daily=="YES")
  {
    Data.daily.other=subset(Data.daily,SPECIES==22999 & !METHOD%in%c("GN","LL")& Bioregion%in%c("SC","WC"))
    Data.daily=subset(Data.daily,!(SPECIES==22999 & !METHOD%in%c("GN","LL") & Bioregion%in%c("SC","WC")))
    
  }
  
  
  Prop.bioregion.yr=function(dat)
  {
    dat1=aggregate(LIVEWT.reap~FINYEAR+Bioregion+SPECIES,subset(dat,SPECIES%in%c(18003,17001,17003,18007)),sum)
    wide=reshape(dat1,v.names="LIVEWT.reap",timevar="SPECIES",idvar=c("FINYEAR","Bioregion"),direction="wide")
    colnames(wide)[3:6]=c("Gummy","Whiskery","Dusky","Sandbar")
    wide[is.na(wide)]=0
    dat2=aggregate(LIVEWT.reap~FINYEAR+Bioregion,dat,sum)
    wide=wide%>%full_join(dat2,by=c("FINYEAR","Bioregion"))
    #wide=merge(wide,dat2,by=c("FINYEAR","Bioregion"))
    wide$Gum.prop=wide$Gummy/wide$LIVEWT.reap
    wide$Whi.prop=wide$Whiskery/wide$LIVEWT.reap
    wide$Dus.prop=wide$Dusky/wide$LIVEWT.reap
    wide$San.prop=wide$Sandbar/wide$LIVEWT.reap
    wide=wide[,c(1:2,8:11)]
    Mean.san.prop=mean(wide$San.prop[wide$San.prop>0])
    wide$San.prop=with(wide,ifelse(San.prop==0,Mean.san.prop,San.prop))
    return(wide)
  }
  
  Prop.bio.yr=Prop.bioregion.yr(subset(Data.monthly,Bioregion%in%c("SC","WC")))
  if(Reapportion.daily=="YES") Prop.bio.yr.daily=Prop.bioregion.yr(subset(Data.daily,Bioregion%in%c("SC","WC")))
  
  #monthly
  Data.monthly.other=merge(Data.monthly.other,Prop.bio.yr,by=c("FINYEAR","Bioregion"),all.x=T)
  Data.monthly.other.gum=Data.monthly.other.whi=Data.monthly.other.dus=
    Data.monthly.other.san=Data.monthly.other
  Data.monthly.other.gum$SPECIES=17001; Data.monthly.other.gum$SNAME="SHARK, GUMMY"
  Data.monthly.other.whi$SPECIES=17003; Data.monthly.other.whi$SNAME="SHARK, WHISKERY"
  Data.monthly.other.dus$SPECIES=18003; Data.monthly.other.dus$SNAME="Shark, Dusky"
  Data.monthly.other.san$SPECIES=18007; Data.monthly.other.san$SNAME="SHARK, THICKSKIN (SANDBAR)"
  Data.monthly.other.gum$LIVEWT.reap=with(Data.monthly.other.gum,LIVEWT.reap*Gum.prop)
  Data.monthly.other.whi$LIVEWT.reap=with(Data.monthly.other.whi,LIVEWT.reap*Whi.prop)
  Data.monthly.other.dus$LIVEWT.reap=with(Data.monthly.other.dus,LIVEWT.reap*Dus.prop)
  Data.monthly.other.san$LIVEWT.reap=with(Data.monthly.other.san,LIVEWT.reap*San.prop)
  Data.monthly.other=rbind(Data.monthly.other.gum,Data.monthly.other.whi,
                           Data.monthly.other.dus,Data.monthly.other.san)
  Data.monthly.other=Data.monthly.other[,-match(c("Gum.prop","Whi.prop","Dus.prop",
                                                  "San.prop"),names(Data.monthly.other))]
  Whaler.ktch.other=sum(Data.monthly.other$LIVEWT.reap[Data.monthly.other$SPECIES%in%c(18003,18007)])
  thss=subset(Data.monthly,SPECIES%in%c(18003,18007) & FINYEAR%in% unique(Data.monthly.other$FINYEAR))
  Whaler.ktch=sum(thss$LIVEWT.reap,na.rm=T)
  Prop.whaler.ktch.other=Whaler.ktch.other/Whaler.ktch
  Data.monthly=rbind(Data.monthly,Data.monthly.other)
  
  
  #daily
  if(Reapportion.daily=="YES")
  {
    if(nrow(Data.daily.other)>0)
    {
      Data.daily.other=merge(Data.daily.other,Prop.bio.yr.daily,by=c("FINYEAR","Bioregion"),all.x=T)
      Data.daily.other.gum=Data.daily.other.whi=Data.daily.other.dus=Data.daily.other.san=Data.daily.other
      Data.daily.other.gum$SPECIES=17001; Data.daily.other.gum$SNAME="SHARK, GUMMY"
      Data.daily.other.whi$SPECIES=17003; Data.daily.other.whi$SNAME="SHARK, WHISKERY"
      Data.daily.other.dus$SPECIES=18003; Data.daily.other.dus$SNAME="Shark, Dusky"
      Data.daily.other.san$SPECIES=18007; Data.daily.other.san$SNAME="SHARK, THICKSKIN (SANDBAR)"
      Data.daily.other.gum$LIVEWT.reap=with(Data.daily.other.gum,LIVEWT.reap*Gum.prop)
      Data.daily.other.whi$LIVEWT.reap=with(Data.daily.other.whi,LIVEWT.reap*Whi.prop)
      Data.daily.other.dus$LIVEWT.reap=with(Data.daily.other.dus,LIVEWT.reap*Dus.prop)
      Data.daily.other.san$LIVEWT.reap=with(Data.daily.other.san,LIVEWT.reap*San.prop)
      Data.daily.other=rbind(Data.daily.other.gum,Data.daily.other.whi,
                             Data.daily.other.dus,Data.daily.other.san)
      Data.daily.other=Data.daily.other[,-match(c("Gum.prop","Whi.prop","Dus.prop",
                                                  "San.prop"),names(Data.daily.other))]
      Whaler.ktch.other.daily=sum(Data.daily.other$LIVEWT.reap[Data.daily.other$SPECIES%in%c(18003,18007)])
      thss.daily=subset(Data.daily,SPECIES%in%c(18003,18007) & FINYEAR%in% unique(Data.daily.other$FINYEAR))
      Whaler.ktch.daily=sum(thss.daily$LIVEWT.reap,na.rm=T)
      Prop.whaler.ktch.other.daily=Whaler.ktch.other.daily/Whaler.ktch.daily
      Data.daily=rbind(Data.daily,Data.daily.other)
      
    }
  }
  
  
}


# C.7.21 REAPPORTION CATCH FROM SHARK BAY
#note: this deals with speciess==22999 from Shark Bay
Data.monthly.other=subset(Data.monthly,SPECIES==22999 & 
                  METHOD%in%c("GN","LL","HL")& LAT>(-26) &
                    LAT<=(-25))
Data.monthly=subset(Data.monthly,!(SPECIES==22999 & 
                METHOD%in%c("GN","LL","HL")& LAT>(-26) & 
                  LAT<=(-25)))

Prop.bioregion.yr=function(dat)
{
  All.yrs=sort(unique(Data.monthly.other$FINYEAR))
  dat1=aggregate(LIVEWT~FINYEAR+SPECIES,subset(dat,SPECIES%in%c(18003,18007)),sum)
  wide=reshape(dat1,v.names="LIVEWT",timevar="SPECIES",idvar=c("FINYEAR"),direction="wide")
  colnames(wide)[2:3]=c("Dusky","Sandbar")
  
  msn.yrs=All.yrs[which(!All.yrs%in%wide$FINYEAR)]
  if(length(msn.yrs)>0)
  {
    add.wide=wide[1:length(msn.yrs),]
    add.wide[,2:3]=NA
    add.wide$FINYEAR=msn.yrs
    wide=rbind(wide,add.wide)
    wide=wide[order(wide$FINYEAR),]
  }
  dat2=aggregate(LIVEWT~FINYEAR,subset(dat,SPECIES%in%Elasmo.species),sum)
  wide=wide%>%left_join(dat2,by=c("FINYEAR"))
  #wide=merge(wide,dat2,by=c("FINYEAR"),all.x=T)
  
  wide$Dus.prop=wide$Dusky/wide$LIVEWT
  wide$San.prop=wide$Sandbar/wide$LIVEWT
  wide[is.na(wide)]=0
  wide=wide[,c(1,5:6)]
  Mean.san.prop=mean(wide$San.prop[wide$San.prop>0])
  Mean.dus.prop=mean(wide$Dus.prop[wide$Dus.prop>0])
  wide$San.prop=with(wide,ifelse(San.prop==0,Mean.san.prop,San.prop))
  wide$Dus.prop=with(wide,ifelse(Dus.prop==0,Mean.dus.prop,Dus.prop))
  TT=rowSums(wide[,2:3])
  wide$Dus.prop=wide$Dus.prop/TT
  wide$San.prop=wide$San.prop/TT
  return(wide)
}
Prop.bio.yr=Prop.bioregion.yr(subset(Data.monthly,LAT>(-26) & LAT<=(-25)))
Data.monthly.other=Data.monthly.other%>%left_join(Prop.bio.yr,by=c("FINYEAR"))
#Data.monthly.other=merge(Data.monthly.other,Prop.bio.yr,by=c("FINYEAR"),all.x=T)
Data.monthly.other.dus=Data.monthly.other.san=Data.monthly.other
Data.monthly.other.dus$SPECIES=18003
Data.monthly.other.dus$SNAME="Shark, Dusky"

Data.monthly.other.san$SPECIES=18007
Data.monthly.other.san$SNAME="SHARK, THICKSKIN (SANDBAR)"

Data.monthly.other.dus$LIVEWT.reap=with(Data.monthly.other.dus,LIVEWT*Dus.prop)
Data.monthly.other.san$LIVEWT.reap=with(Data.monthly.other.san,LIVEWT*San.prop)

Data.monthly.other=rbind(Data.monthly.other.dus,Data.monthly.other.san)
Data.monthly.other=Data.monthly.other[,-match(c("Dus.prop","San.prop"),names(Data.monthly.other))]
Data.monthly=rbind(Data.monthly,Data.monthly.other)


#Compare amended weights, loosing or gaining catch?
d2=subset(Data.monthly,LAT<=(-26) & !is.na(LIVEWT.reap) &TYPE.DATA=="monthly")
fn.chk.ktch(d1=subset(Data.monthly.original,Same.return
                      %in%unique(d2$Same.return)),
            d2=d2,
            VAR1="LIVEWT",VAR2="LIVEWT.reap")
#if "d1 and d2 have same catch" then OK


#Remove records with reported catch set at 0 or NA
Data.monthly=subset(Data.monthly,!LIVEWT.reap==0)
Data.daily=subset(Data.daily,!LIVEWT.reap==0)

Data.monthly=subset(Data.monthly,!is.na(LIVEWT.reap))
Data.daily=subset(Data.daily,!is.na(LIVEWT.reap))

#by zone
Hndl="C:/Matias/Analyses/Catch and effort/Compare_to_previous_approach/"
fun.compare.weight=function(SPEC)
{
  dat=subset(Data.monthly,SPECIES==SPEC & zone%in%c("West","Zone1","Zone2"))
  dat1=subset(dat,Spec.old==SPEC)
  spe=unique(dat$SNAME)[1]
  tab1=aggregate(LIVEWT~FINYEAR+zone,dat1,sum,na.rm=T)
  tab2=aggregate(LIVEWT.reap~FINYEAR+zone,dat,sum,na.rm=T)
  wide1 <- reshape(tab1, v.names ="LIVEWT" , idvar = "FINYEAR",timevar = "zone", direction = "wide")
  wide2 <- reshape(tab2, v.names ="LIVEWT.reap" , idvar = "FINYEAR",timevar = "zone", direction = "wide")
  Yer=1:length(wide1$FINYEAR)
  Yer2=1:length(wide2$FINYEAR)
  if(length(Yer)<length(Yer2))wide2=wide2[match(wide1$FINYEAR,wide2$FINYEAR),]
  Yer2=1:length(wide2$FINYEAR)
  plot(Yer,wide1$LIVEWT.West/1000,col=2,ylim=c(0,max(c(tab1$LIVEWT/1000,tab2$LIVEWT.reap/1000))),
       ylab="",xlab="",xaxt='n',main=spe,pch=19,cex.axis=1.25)
  lines(Yer2,wide2$LIVEWT.reap.West/1000,col=2,pch=19)
  #points(Yer2,wide2$LIVEWT.reap.West,col=2,pch=21)
  points(Yer,wide1$LIVEWT.Zone1/1000,col=3,pch=19)
  lines(Yer2,wide2$LIVEWT.reap.Zone1/1000,col=3)
  #points(Yer2,wide2$LIVEWT.reap.Zone1,col=3,pch=21)
  points(Yer,wide1$LIVEWT.Zone2/1000,col=4,pch=19)
  lines(Yer2,wide2$LIVEWT.reap.Zone2/1000,col=4)
  #points(Yer2,wide2$LIVEWT.reap.Zone2,col=4,pch=21)
  axis(1,Yer,labels=F,tck=-0.015)
  axis(1,seq(1,length(Yer),by=5),labels=wide1$FINYEAR[seq(1,length(Yer),by=5)],tck=-0.03,cex.axis=1.25)
}
tiff(file=paste(Hndl,"Compare.reapKtch_zone.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,2),mar=c(3,4,1,4),oma=c(2.5,.6,.1,.1),las=1)
fun.compare.weight(17001)
legend("topleft",c("LiveWT","LiveWt.reap"),bty='n',pch=c(21,NA),lty=c(0,1),pt.bg=c("black","white"))
fun.compare.weight(17003)
legend("topright",c("West","Zn1","Zn2"),bty='n',lty=1,col=c(2,3,4))
fun.compare.weight(18003)
fun.compare.weight(18007)
mtext("Tonnes",2,outer=T,cex=1.75,line=-1.5,las=3)
mtext("Financial year",1,outer=T,cex=1.75)
dev.off()

  #C.8 Add 5% to catch records prior to 1990  as per instructed by Rory                                              
Data.monthly$LIVEWT.c=with(Data.monthly,ifelse(YEAR.c<1990 & LAT<=(-26),LIVEWT.reap*Inc.per,LIVEWT.reap))
Data.daily$LIVEWT.c=Data.daily$LIVEWT.reap


#Check for duplications in Monthly
  #note: no duplications in daily
# Data.monthly$dupli=with(Data.monthly,paste(Same.return,SPECIES,LIVEWT.c,TYPE.DATA,CONDITN))
# a=Data.monthly[duplicated(Data.monthly$dupli),]
# b=subset(Data.monthly,dupli%in%unique(a$dupli)&TYPE.DATA=="monthly",select=c(dupli))
# a=subset(Data.monthly,dupli%in%unique(b$dupli))
# a=a[order(a$dupli),]
# write.csv(a,"duplicated.records.csv",row.names=F)
# Data.monthly=subset(Data.monthly,!dupli%in%unique(b$dupli))
# a=a[!duplicated(a$dupli),]
# Data.monthly=rbind(Data.monthly,a)
# Data.monthly=Data.monthly[,-match("dupli",names(Data.monthly))]


#Extract data for FishCUBE
  #Monthly
FishCUBE.monthly=subset(Data.monthly,FINYEAR%in%c(unique(Data.monthly.CAESS$FINYEAR),keep),
     select=c(Same.return,FisheryCode,METHOD,Landing.Port,LAT,LONG,
         FINYEAR,YEAR.c,MONTH,FisheryZone,Bioregion,
         blockxFC,LIVEWT.c,VESSEL,SPECIES,SNAME,RSCommonName,
         RSSpeciesId,licence,BDAYS,rowid) )  
FishCUBE.monthly=subset(FishCUBE.monthly,!FINYEAR%in%unique(Data.daily$FINYEAR))          #remove post 2005 records already in daily

Monthly.not.in.daily$LIVEWT.c=Monthly.not.in.daily$LIVEWT    
FishCUBE.monthly=rbind(FishCUBE.monthly,Monthly.not.in.daily[,match(names(FishCUBE.monthly),
                                                    names(Monthly.not.in.daily))])     #Add post 2005-06 records not in Daily
FishCUBE.monthly$TYPE.DATA="monthly"


  #Daily
FishCUBE.daily=Data.daily
FishCUBE.daily=FishCUBE.daily[,-match("BLOCKX",names(FishCUBE.daily))]
Data.monthly=Data.monthly[,-match("rowid",names(Data.monthly))]
Data.daily=Data.daily[,-match("rowid",names(Data.daily))]


  #C.9 Separate fisheries  into North and South                                                                   
#note: no catch reapportioning of effort corrections done on northern data file
  #monthly
Data.monthly.north=subset(Data.monthly,LAT>(-26))   
Data.monthly=subset(Data.monthly,LAT<=(-26))

  #daily
Data.daily.north=subset(Data.daily,LAT>(-26))   
Data.daily=subset(Data.daily,LAT<=(-26))


# C.9.1 Northern and Southern catches
fun.compare.weight.NORTH.tot=function(DAT,SPEC)
{
  dat=subset(DAT,SPECIES==SPEC)
  dat1=subset(dat,Spec.old==SPEC)
  spe=unique(dat$SNAME)[1]
  tab1=aggregate(LIVEWT~FINYEAR,dat1,sum,na.rm=T)
  tab2=aggregate(LIVEWT.reap~FINYEAR,dat,sum,na.rm=T)
  wide1 <- tab1
  wide2 <- tab2
  Yer=1:length(wide1$FINYEAR)
  Yer2=1:length(wide2$FINYEAR)
  if(length(Yer)<length(Yer2))wide2=wide2[match(wide1$FINYEAR,wide2$FINYEAR),]
  Yer2=1:length(wide2$FINYEAR)
  plot(Yer,wide1$LIVEWT/1000,col=2,ylim=c(0,max(c(tab1$LIVEWT/1000,tab2$LIVEWT.reap/1000))),ylab="Tons",xlab="FinYr",
       xaxt='n',main=spe,pch=19)
  lines(Yer2,wide2$LIVEWT.reap/1000,col=3,pch=19)
  axis(1,Yer,labels=F,tck=-0.015)
  axis(1,seq(1,length(Yer),by=5),labels=wide1$FINYEAR[seq(1,length(Yer),by=5)],tck=-0.03)
  legend("topleft",c("LiveWT","LiveWt.reap"),bty='n',pch=c(21,NA),col=c(2,3),lty=c(0,1),pt.bg=c("red","white"))
}

  #North
Show.north="NO"  #change to YES if doing catch reapportions for NSF
if(Show.north=="YES")
{  
  tiff(file=paste(Hndl,"Compare.reapKtch.North.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,2),las=1)
  fun.compare.weight.NORTH.tot(subset(Data.monthly.north,!is.na(LIVEWT.reap)),17001)
  fun.compare.weight.NORTH.tot(subset(Data.monthly.north,!is.na(LIVEWT.reap)),17003)
  fun.compare.weight.NORTH.tot(subset(Data.monthly.north,!is.na(LIVEWT.reap)),18003)
  fun.compare.weight.NORTH.tot(subset(Data.monthly.north,!is.na(LIVEWT.reap)),18007)
  dev.off()
}

  #South
tiff(file=paste(Hndl,"Compare.reapKtch.South.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,2),las=1)
fun.compare.weight.NORTH.tot(Data.monthly,17001)
fun.compare.weight.NORTH.tot(Data.monthly,17003)
fun.compare.weight.NORTH.tot(Data.monthly,18003)
fun.compare.weight.NORTH.tot(Data.monthly,18007)
dev.off()


#Compare with Rory's previous summaries
fn.check=function(dat,WHAT)
{
  FInYEAR=unique(Results.pre.2013$FINYEAR)
  NN=length(FInYEAR)
  if(WHAT=="LIVEWT") this=aggregate(LIVEWT~FINYEAR,dat,sum)
  if(WHAT=="LIVEWT.reap") this=aggregate(LIVEWT.reap~FINYEAR,dat,sum)
  if(WHAT=="LIVEWT.c") this=aggregate(LIVEWT.c~FINYEAR,dat,sum)
  
  plot(Results.pre.2013$TDGDLF.tot.sk.live.wt,type='l',ylim=c(0,2000),ylab="Total shark catch",xaxt='n',main=WHAT)
  lines(this[,2]/1000,col=2)
  axis(1,at=1:NN,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
  legend("topleft",c("Rory's","new"),bty='n',lty=1,col=1:2)
  
}
tiff(file=paste(Hndl,"Compare.Rory.current.ktch.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
fn.check(subset(Data.monthly,SPECIES%in%Elasmo.species & METHOD%in%c("GN","LL") &
                  !(CONDITN %in% c("DR","LI","FI")) & LAT<=(-26)),"LIVEWT.c")
dev.off()


#by zone 
fn.check.zn=function(dat,WHAT)
{
  if(WHAT=="LIVEWT")
  {
    annual.catch.by.zone=aggregate(LIVEWT~FINYEAR+zone,data=dat,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.reap")
  {
    annual.catch.by.zone=aggregate(LIVEWT.reap~FINYEAR+zone,data=dat,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.reap",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.c")
  {
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  
  plot(Results.pre.2013$WC.tot.sk.live.wt,type='l',ylim=c(0,1200),ylab="Total shark catch",main=WHAT)
  lines(Results.pre.2013$Z1.tot.sk.live.wt,lty=2)
  lines(Results.pre.2013$Z2.tot.sk.live.wt,lty=3)
  LTY=1:3;COL=2:4
  for (i in 1:(ncol(wide)-1)) lines(wide[,i+1]/1000,lty=LTY[i],col=COL[i])
  
  legend("topleft",c("WC","Zn1","Zn2"),bty='n',lty=1:3)
}
tiff(file=paste(Hndl,"Compare.Rory.current.ktch.Zone.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
fn.check.zn(subset(Data.monthly,SPECIES%in%Elasmo.species & METHOD%in%c("GN","LL") & LAT<=(-26)),"LIVEWT.reap")
dev.off()


#by zone and species
Read=read.csv("C:/Matias/Data/Catch and Effort/Historic/Spec.catch.zone.csv")
fn.by.species=function(dat,WHAT)
{
  FInYEAR=unique(Read$Fin.yr)
  NN=length(FInYEAR)
  
  #dusky
  d1=subset(dat,SPECIES%in%c(18001,18003))
  if(WHAT=="LIVEWT")
  {
    annual.catch.by.zone=aggregate(LIVEWT~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.reap")
  {
    annual.catch.by.zone=aggregate(LIVEWT.reap~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.reap",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  
  if(WHAT=="LIVEWT.c")
  {
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")    
  }  
  plot(Read$Dusky.catch.WC,type='l',ylim=c(0,350),ylab="Total catch",xaxt='n',main=paste("Dusky",WHAT))
  lines(Read$Dusky.catch.Zn1,lty=2)
  lines(Read$Dusky.catch.Zn2,lty=3)
  LTY=1:3;COL=2:4
  for (i in 1:(ncol(wide)-1)) lines(wide[,i+1]/1000,lty=LTY[i],col=COL[i])
  axis(1,at=1:NN,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
  legend("topleft",c("Rory's","new"),bty='n',lty=1,col=1:2)
  legend("topright",c("WC","Zn1","Zn2"),bty='n',lty=1:3)
  
  #sandbar
  d1=subset(dat,SPECIES%in%c(18007))
  if(WHAT=="LIVEWT")
  {
    annual.catch.by.zone=aggregate(LIVEWT~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.reap")
  {
    annual.catch.by.zone=aggregate(LIVEWT.reap~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.reap",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  
  if(WHAT=="LIVEWT.c")
  {
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")    
    #wide=wide[-1,]
  }  
  Ns=nrow(Read)
  plot(Read$Sandbar.catch.WC[11:Ns],type='l',ylim=c(0,350),ylab="Total catch",xaxt='n',main=paste("Sandbar",WHAT))
  lines(Read$Sandbar.catch.Zn1[11:Ns],lty=2)
  lines(Read$Sandbar.catch.Zn2[11:Ns],lty=3)
  LTY=1:3;COL=2:4
  for (i in 1:(ncol(wide)-1)) lines(wide[(match("1985-86",wide$FINYEAR)):nrow(wide),i+1]/1000,lty=LTY[i],col=COL[i])
  id=match("1985-86",FInYEAR)
  axis(1,at=1:nrow(wide),labels=F,tck=-0.01)
  axis(1,at=seq(1,nrow(wide),5),labels=FInYEAR[id:Ns][seq(1,nrow(wide),5)],tck=-0.02,cex.axis=1.1)
  legend("topleft",c("Rory's","new"),bty='n',lty=1,col=1:2)
  legend("topright",c("WC","Zn1","Zn2"),bty='n',lty=1:3)
  
  
  #gummy
  d1=subset(dat,SPECIES%in%c(17001))
  if(WHAT=="LIVEWT")
  {
    annual.catch.by.zone=aggregate(LIVEWT~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.reap")
  {
    annual.catch.by.zone=aggregate(LIVEWT.reap~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.reap",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.c")
  {
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")    
  }  
  plot(Read$Gummy.catch.WC,type='l',ylim=c(0,700),ylab="Total catch",xaxt='n',main=paste("Gummy",WHAT))
  lines(Read$Gummy.catch.Zn1,lty=2)
  lines(Read$Gummy.catch.Zn2,lty=3)
  LTY=1:3;COL=2:4
  for (i in 1:(ncol(wide)-1)) lines(wide[,i+1]/1000,lty=LTY[i],col=COL[i])
  axis(1,at=1:NN,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
  legend("topleft",c("Rory's","new"),bty='n',lty=1,col=1:2)
  legend("topright",c("WC","Zn1","Zn2"),bty='n',lty=1:3)
  
  
  #Whiskery
  d1=subset(dat,SPECIES%in%c(17003))
  if(WHAT=="LIVEWT")
  {
    annual.catch.by.zone=aggregate(LIVEWT~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.reap")
  {
    annual.catch.by.zone=aggregate(LIVEWT.reap~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.reap",timevar="zone",idvar="FINYEAR",direction="wide")    
  }
  if(WHAT=="LIVEWT.c")
  {
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=d1,sum,na.rm=T)  
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")    
  }  
  plot(Read$Whiskery.catch.WC,type='l',ylim=c(0,350),ylab="Total catch",xaxt='n',main=paste("Whiskery",WHAT))
  lines(Read$Whiskery.catch.Zn1,lty=2)
  lines(Read$Whiskery.catch.Zn2,lty=3)
  LTY=1:3;COL=2:4
  for (i in 1:(ncol(wide)-1)) lines(wide[,i+1]/1000,lty=LTY[i],col=COL[i])
  axis(1,at=1:NN,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
  legend("topleft",c("Rory's","new"),bty='n',lty=1,col=1:2)
  legend("topright",c("WC","Zn1","Zn2"),bty='n',lty=1:3)
}
tiff(file=paste(Hndl,"Compare.Rory.current.ktch.Zone.and.species.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,2),las=1)
fn.by.species(subset(Data.monthly,SPECIES%in%Elasmo.species & METHOD%in%c("GN","LL") & LAT<=(-26)),"LIVEWT.c")
dev.off()


# C.10  Check fishing methods used
Methods=unique(as.character(Data.monthly$METHOD))
Methods.north=unique(as.character(Data.monthly.north$METHOD))


# C.11. Sort out species
All.species=Data.monthly[,match(c("SNAME","SPECIES"),names(Data.monthly))]
All.species$dummy=with(All.species,paste(SNAME,SPECIES))
All.species=All.species[!duplicated(All.species$dummy),]
All.species=All.species[order(All.species$SPECIES),1:2]



#C.13. Proportion of bad catch reporters
  #monthly
Bad.reps.all=table(Data.monthly$Reporter.old,useNA='ifany')
Percent.Bad.reps.all=100*Bad.reps.all[1]/sum(Bad.reps.all)

  #daily
if(Reapportion.daily=="YES")
{
  Bad.reps.all.daily=table(Data.daily$Reporter.old,useNA='ifany')
  Percent.Bad.reps.all.daily=100*Bad.reps.all.daily[1]/sum(Bad.reps.all.daily)  
}


#SECTION D. ---- EFFORT INSPECTIONS -----
#note: this part sets functions for inspecting current year data and raise issues
#       to data-entry girls
#       Visually inspect each plot and raise any inconsistent record to data-entry girls
#       This has already been done for all Daily Records prior to current year.

Dummy.yr=sort(unique(Effort.daily$finyear))
par(mfcol=c(1,1),mai=c(1.3,1.3,.2,.2))

if(Inspect.New.dat=="YES")
{
  Inspect.vars=c("NETLEN","BDAYS","HOURS","SHOTS","nlines")
  
  #id typos and nonsense
  fun.inspect=function(Data,YRS,Type,Min.lim,max.var,VAR) 
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    
    if(Type=="Monthly")
    {
      Data=Data[order(Data$VESSEL,Data$FINYEAR,Data$MONTH),]
      Data=Data[!(duplicated(Data$Same.return)),]
      Unik.ves=as.character(unique(Data$VESSEL))
      these.vars=match(c("VESSEL","Same.return","FINYEAR","MONTH",VAR),names(Data))
      Check=data.frame(Data[1:(nrow(Data)-1),these.vars],Data[2:nrow(Data),match(c("Same.return",VAR),names(Data))])
      Check$diff=Check[,match(VAR,names(Check))]-Check[,match(paste(VAR,".1",sep=""),names(Check))]
      Check$flag=with(Check,ifelse(diff<0 &abs(diff)/Min.lim>(Min.lim*.1),paste("Check ",VAR,".1",sep=""),ifelse(diff<0 &abs(diff)/Min.lim <(Min.lim*.1),
                                                                                                                 paste("Check ",VAR,sep=""),ifelse(diff>0 &diff/Min.lim >(Min.lim*.1),paste("Check ",VAR,sep=""),"OK"))))
      Check$flag=ifelse(abs(Check$diff)<=Check[,match(VAR,names(Check))]*(1+max.var) &abs(Check$diff)>=Min.lim &abs(Check$diff)<(Min.lim*10),"OK",Check$flag)
      Check$flag=ifelse(Check[,match(paste(VAR,".1",sep=""),names(Check))]<=(Check[,match(VAR,names(Check))]/10),paste("Check ",VAR,".1",sep=""),Check$flag)
      Check$Select=ifelse(Check$flag==paste("Check ",VAR,".1",sep=""),Check$Same.return.1,
                          ifelse(Check$flag==paste("Check ",VAR,sep=""),Check$Same.return,"OK"))
      Check=subset(Check,!(Select=="OK"),select=c(VESSEL,Select))
      rownames(Check)=NULL
      names(Check)[2]="Same.return"
      Check.vesl=subset(Data,VESSEL%in%as.character(Check$VESSEL))
      Check.vesl$Check=with(Check.vesl,ifelse(Same.return%in%Check$Same.return,"check","OK"))
      OK.vesl=subset(Check.vesl,Check=="OK")
      OK.vesl=OK.vesl[,match(c(VAR,"VESSEL"),names(OK.vesl))]
      names(OK.vesl)[1]="Dodgy.var"
      
      Mode.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl, function(x) as.numeric(names(which.max(table(x)))))
      names(Mode.net.len)[2]=paste("Mode.",VAR,sep="")
      Mean.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,mean)
      names(Mean.net.len)[2]=paste("Mean.",VAR,sep="")
      Mean.net.len=merge(Mean.net.len,Mode.net.len,by="VESSEL",all.x=T)
      SD.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,sd)
      names(SD.net.len)[2]=paste("SD.",VAR,sep="")
      Mean.net.len=merge(Mean.net.len,SD.net.len,by="VESSEL",all.x=T)
      Min.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,min)
      names(Min.net.len)[2]=paste("Min.",VAR,sep="")
      Max.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,max)
      names(Max.net.len)[2]=paste("Max.",VAR,sep="")
      Range.net.len=merge(Min.net.len,Max.net.len,by="VESSEL",all.x=T)
      Mean.net.len=merge(Mean.net.len,Range.net.len,by="VESSEL",all.x=T)
      Check.vesl=merge(Check.vesl,Mean.net.len,by="VESSEL",all.x=T)
      Check.vesl$Check=ifelse(Check.vesl[match(VAR,names(Check.vesl))]<=(Check.vesl[match(paste("Mode.",VAR,sep=""),names(Check.vesl))]*1.2) &
                                Check.vesl[match(VAR,names(Check.vesl))]>=(Check.vesl[match(paste("Mode.",VAR,sep=""),names(Check.vesl))]*0.8),"OK",Check.vesl$Check)
      Check.vesl$Check=ifelse(nchar(Check.vesl[match(VAR,names(Check.vesl))][,1])<3 |nchar(Check.vesl[match(VAR,names(Check.vesl))][,1])>4,
                              "check",Check.vesl$Check)
      Check.vesl=subset(Check.vesl,Check=="check" & !is.na(Check.vesl[,match(paste("SD.",VAR,sep=""),
                                                                             names(Check.vesl))]))
      Check.vesl=Check.vesl[,match(c("VESSEL","Same.return","FINYEAR","MONTH",VAR,paste("Mode.",VAR,sep=""),paste("Mean.",VAR,sep=""), 
                                     paste("SD.",VAR,sep=""),paste("Min.",VAR,sep=""),paste("Max.",VAR,sep="")),names(Check.vesl))]
    }
    if(Type=="Daily")
    {
      Data=Data[order(Data$VESSEL,Data$year,Data$MONTH),]
      Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))
      
      Data=Data[!(duplicated(Data$Session_ID)),]
      Unik.ves=as.character(unique(Data$VESSEL))
      
      these.vars=match(c("VESSEL","Session_ID","FINYEAR","MONTH",VAR),names(Data))
      Check=data.frame(Data[1:(nrow(Data)-1),these.vars],Data[2:nrow(Data),match(c("Session_ID",VAR),names(Data))])
      Check$diff=Check[,match(VAR,names(Check))]-Check[,match(paste(VAR,".1",sep=""),names(Check))]
      Check$flag=with(Check,ifelse(diff<0 &abs(diff)/Min.lim>(Min.lim*.1),paste("Check ",VAR,".1",sep=""),ifelse(diff<0 &abs(diff)/Min.lim <(Min.lim*.1),
                                                                                                                 paste("Check ",VAR,sep=""),ifelse(diff>0 &diff/Min.lim >(Min.lim*.1),paste("Check ",VAR,sep=""),"OK"))))
      Check$flag=ifelse(abs(Check$diff)<=Check[,match(VAR,names(Check))]*(1+max.var) &abs(Check$diff)>=Min.lim &abs(Check$diff)<(Min.lim*10),"OK",Check$flag)
      
      Check$flag=ifelse(Check[,match(paste(VAR,".1",sep=""),names(Check))]<=(Check[,match(VAR,names(Check))]/10),paste("Check ",VAR,".1",sep=""),Check$flag)
      
      Check$Select=ifelse(Check$flag==paste("Check ",VAR,".1",sep=""),Check$Session_ID.1,
                          ifelse(Check$flag==paste("Check ",VAR,sep=""),Check$Session_ID,"OK"))
      Check=subset(Check,!(Select=="OK"),select=c(VESSEL,Select))
      rownames(Check)=NULL
      names(Check)[2]="Session_ID"
      Check.vesl=subset(Data,VESSEL%in%as.character(Check$VESSEL))
      Check.vesl$Check=with(Check.vesl,ifelse(Session_ID%in%Check$Session_ID,"check","OK"))
      OK.vesl=subset(Check.vesl,Check=="OK")
      OK.vesl=OK.vesl[,match(c(VAR,"VESSEL"),names(OK.vesl))]
      names(OK.vesl)[1]="Dodgy.var"
      
      Mode.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl, function(x) as.numeric(names(which.max(table(x)))))
      names(Mode.net.len)[2]=paste("Mode.",VAR,sep="")
      
      Mean.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,mean)
      names(Mean.net.len)[2]=paste("Mean.",VAR,sep="")
      
      Mean.net.len=merge(Mean.net.len,Mode.net.len,by="VESSEL",all.x=T)
      
      SD.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,sd)
      names(SD.net.len)[2]=paste("SD.",VAR,sep="")
      Mean.net.len=merge(Mean.net.len,SD.net.len,by="VESSEL",all.x=T)
      Min.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,min)
      names(Min.net.len)[2]=paste("Min.",VAR,sep="")
      Max.net.len=aggregate(Dodgy.var~VESSEL,data=OK.vesl,max)
      names(Max.net.len)[2]=paste("Max.",VAR,sep="")
      Range.net.len=merge(Min.net.len,Max.net.len,by="VESSEL",all.x=T)
      Mean.net.len=merge(Mean.net.len,Range.net.len,by="VESSEL",all.x=T)
      Check.vesl=merge(Check.vesl,Mean.net.len,by="VESSEL",all.x=T)
      Check.vesl$Check=ifelse(Check.vesl[match(VAR,names(Check.vesl))]<=(Check.vesl[match(paste("Mode.",VAR,sep=""),names(Check.vesl))]*1.2) &
                                Check.vesl[match(VAR,names(Check.vesl))]>=(Check.vesl[match(paste("Mode.",VAR,sep=""),names(Check.vesl))]*0.8),"OK",Check.vesl$Check)
      Check.vesl$Check=ifelse(nchar(Check.vesl[match(VAR,names(Check.vesl))][,1])<3 |nchar(Check.vesl[match(VAR,names(Check.vesl))][,1])>4,
                              "check",Check.vesl$Check)
      Check.vesl=subset(Check.vesl,Check=="check" & !is.na(Check.vesl[,match(paste("SD.",VAR,sep=""),
                                                                             names(Check.vesl))]))
      Check.vesl=Check.vesl[,match(c("VESSEL","MONTH","TSNo","DSNo","SNo",VAR,paste("Mode.",VAR,sep=""),paste("Mean.",VAR,sep=""),
                                     paste("SD.",VAR,sep=""),paste("Min.",VAR,sep=""),paste("Max.",VAR,sep="")),names(Check.vesl))]
      
    }
    return(Check.vesl)
  }
  
  #visual inspection of variability by vessel
  fun.plot.eff.var=function(Data,YRS,VAR,lim1,lim2,threshold,withCols)     
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(VAR,names(Data))
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Unik.ves=as.character(unique(Data$VESSEL))
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session),]
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      if(length(YRS>1))TIT=paste(c(YRS[1],"to",YRS[length(YRS)],Unik.ves[i]),collapse=" ")else
        TIT=paste(YRS,"",Unik.ves[i])
      
      plot(da$dummy,da[,id],type="b",xlab=expression(paste("session n",degree,sep="")),
           ylab=VAR,main=TIT,cex.main=.8,ylim=c(lim1,lim2))
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
        }      
      }
      
    }
  }
  
  #further checks of inconsistencies and variability
  multi.ggplot <- function(plots, cols=1, layout=NULL)
  {
    library(grid)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout))
    {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) 
    {
      print(plots[[1]])
      
    }else{
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  fun.further.chk=function(Data,YRS,VAR)    
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN" & VESSEL %in%check.vesl)
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(VAR,names(Data))
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session_ID),]
    
    VESL=Data[,match(c(VAR,"VESSEL"),names(Data))]
    names(VESL)[1]="Dodgy.var"
    
    if(VAR=="NETLEN")
    {
      Mode.net.len=aggregate(Dodgy.var~VESSEL,data=VESL, function(x) as.numeric(names(which.max(table(x)))))
      names(Mode.net.len)[2]="Mode"
      VESL=merge(VESL,Mode.net.len,by="VESSEL",all.x=T)
      Data$Check=ifelse(VESL$Dodgy.var>=(VESL$Mode*.8) & VESL$Dodgy.var<=(VESL$Mode*1.2),"OK","check")
      Data=merge(Data,Mode.net.len,by="VESSEL",all.x=T)
      Check.vesl=subset(Data,Check=="check",select=c("VESSEL","TSNo","DSNo","SNo",VAR,"Mode"))
    }
    if(!(VAR=="NETLEN"))
    {
      id=match(VAR,names(Data))
      Data$Ok.var=ifelse(Data[,id]<MIN|Data[,id]>MAX,"Problem","OK")
      Check.vesl=Data[,match(c("VESSEL","MONTH","TSNo","DSNo","SNo",VAR,"Ok.var"),names(Data))]
    }  
    p=vector('list',length(check.vesl))
    for(k in 1:length(check.vesl))
    {
      p[[k]]=ggplot(subset(Check.vesl,VESSEL==check.vesl[k]),aes_string(x=VAR,group='Ok.var',fill='Ok.var'))+
        geom_histogram(position="dodge",binwidth=1)+
        theme_bw()+ggtitle(check.vesl[k])
    }
    multi.ggplot(plots=p)
    return(Check.vesl)
  }
  
  #visual insp combine hours and shots
  fun.plot.combo=function(Data,YRS,VAR,VAR1,lim1,lim2,threshold,withCols)     
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(c(VAR,VAR1),names(Data))
    Data$combo=Data[,id[1]]*Data[,id[2]]
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Unik.ves=as.character(unique(Data$VESSEL))
    
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session),]
    
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da$combo,type="b",xlab="session n?",ylab=paste(VAR,"x",VAR1),main=paste(YRS,"",Unik.ves[i]),ylim=c(lim1,lim2))
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt$combo,col=Cols[j],type='b')
        }      
      }
      
    }
  }
  
  #further checks of inconsistencies and variability
  fun.further.chk.combo=function(Data,YRS,VAR,VAR1)    
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN" & VESSEL %in%check.vesl)
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(c(VAR,VAR1),names(Data))
    Data$combo=Data[,id[1]]*Data[,id[2]]
    
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session_ID),]
    Data$Ok.var=ifelse(Data$combo<MIN|Data$combo>MAX,"Problem","OK")
    Check.vesl=subset(Data,select=c("VESSEL","TSNo","DSNo","SNo",VAR,VAR1,"combo","Ok.var"))
    
    p=vector('list',length(check.vesl))
    for(k in 1:length(check.vesl))
    {
      p[[k]]=ggplot(subset(Check.vesl,VESSEL==check.vesl[k]),aes_string(x=VAR,group='Ok.var',fill='Ok.var'))+
        geom_histogram(position="dodge",binwidth=1)+
        theme_bw()+ggtitle(check.vesl[k])
    }
    multi.ggplot(plots=p)
    
    
    return(Check.vesl)
  }
  
  #visual insp of km.gn.hours
  fun.plot.km.gn.hours=function(Data,YRS,VAR,VAR1,VAR2,withCols)     
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(c(VAR,VAR1,VAR2),names(Data))
    Data$km.gn.hours=Data[,id[1]]*Data[,id[2]]*Data[,id[3]]/1000
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Unik.ves=as.character(unique(Data$VESSEL))
    
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session),]
    
    ddd=vector('list',length(Unik.ves))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da$km.gn.hours,type="b",xlab="session number",ylab="km.gn.hours",main=paste(YRS,"",Unik.ves[i]))
      Med=median(da$km.gn.hours,na.rm=T)
      Up=quantile(da$km.gn.hours,.975,na.rm=T)
      Low=quantile(da$km.gn.hours,.025,na.rm=T)
      abline(h=Med,lwd=2)
      polygon(c(da$dummy,rev(da$dummy)),c(rep(Low,nrow(da)),rep(Up,nrow(da))),
              border='transparent',col=rgb(.1,.3,.1,alpha=.2))
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt$km.gn.hours,col=Cols[j],type='b')
        }      
      }
      da=da%>%
        mutate(VESSEL=as.character(VESSEL),
               TSNo=as.character(TSNo),
               Check=ifelse(km.gn.hours>2*Med,"Above 2 x media","OK"))
      iii=which(da$TSNo%in%unique(subset(da,Check=="Above 2 x media")$TSNo))
      if(length(iii)>0)ddd[[i]]=da[iii,c(match(c("VESSEL","TSNo","DSNo","km.gn.hours","Check"),names(da)),id)]
    }
    return(do.call(rbind,ddd))
  }
  
  #bdays by month
  #note: bday is number of fishing days per block per month
  fun.plot.eff.days=function(Data,yr,what)
  {
    Data=subset(Data,FINYEAR%in%yr)
    Unik.ves=as.character(unique(Data$VESSEL))
    
    Data=Data[!duplicated(Data$TSNo),]
    fDaYs=aggregate(BDAYS~VESSEL+MONTH,Data,sum)
    Data=fDaYs[order(fDaYs$VESSEL,fDaYs$MONTH),]
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      plot(da$MONTH,da$BDAYS,type="b",xlab="Month",ylab="BDAYS",pch=19,main=paste(yr,"",Unik.ves[i]),ylim=c(0,max(Data$BDAYS)))
      abline(30,0,col=2,lwd=2)
      text(5,30*1.15,"30 days",col=2)
    }
    
  }
  
  #fdays by month
  #note: fdays in number of fishing days per month
  fun.plot.fishing.days.month=function(Data,YRS,VAR,lim1,lim2,threshold)     
  {
    Data$LAT=ifelse(is.na(Data$LAT),-as.numeric(substr(Data$blockx,1,2)),Data$LAT)
    Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(VAR,names(Data))
    Data$Session_ID=with(Data,paste(FINYEAR,MONTH,date,VESSEL,sep="_"))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Unik.ves=as.character(unique(Data$VESSEL))
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session),]
    
    #Data$Trip=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data$dummy=1
    #Ag=aggregate(dummy~date+FINYEAR+MONTH+VESSEL+DSNo,Data,sum)
    Ag=aggregate(dummy~FINYEAR+MONTH+VESSEL,Data,sum)
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$MONTH)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      da=da[!duplicated(da$MONTH),]
      if(length(YRS>1))TIT=paste(c(YRS[1],"to",YRS[length(YRS)],Unik.ves[i]),collapse=" ")else
        TIT=paste(YRS,"",Unik.ves[i])
      da=da[order(da$MONTH),]
      plot(da$MONTH,da[,id],type="b",xlab="Month",ylab=VAR,main=TIT,cex.main=.8,ylim=c(lim1,lim2))
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      
    }
  }
  
  ##all years daily logbooks nlines X netlen
  ChEck.nline_netlen="NO"
  if(ChEck.nline_netlen=="NO")
  {
    hndl.net_line="C:\\Matias\\Analyses\\Catch and effort\\Outputs\\Netlen_nlines\\"
    fun.plot.nlines.netlen=function(Data,YRS,lim1,lim2,lim11,lim22,threshold,threshold1)     
    {
      Data=subset(Data, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
      Data=Data[!duplicated(Data$Same.return.SNo),]
      id=match('nlines',names(Data))
      id2=match('NETLEN',names(Data))
      Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
      Data=Data[!(duplicated(Data$Session_ID)),]
      Unik.ves=as.character(unique(Data$VESSEL))
      Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
      Data=Data[order(Data$VESSEL,Data$year,Data$MONTH,Data$Session),]
      
      for (i in 1:length(Unik.ves))
      {
        par(mfcol=c(2,1),mar=c(2,5,1,6),mgp=c(2,.8,0))
        da=subset(Data,VESSEL==Unik.ves[i])
        Unik.sess=unique(da$Session)
        NN=length(Unik.sess)
        da$dummy=1:nrow(da)
        if(length(YRS>1))TIT=paste(c(YRS[1],"to",YRS[length(YRS)],Unik.ves[i]),collapse=" ")else
          TIT=paste(YRS,"",Unik.ves[i])
        
        plot(da$dummy,da[,id],type='o',col=2,xlab="session n?",ylab="",main=TIT,
             cex.main=.8,ylim=c(lim1,lim2),yaxt="n")
        abline(threshold,0,col=2,lwd=2)
        text(1,threshold*0.9,threshold,col=2)
        axis(side = 2,col.axis=2)
        mtext(side = 2, line = 2, 'nlines',col=2)
        
        
        
        #netlen
        par(new = T)
        plot(da$dummy,da[,id2],type='o',col=3,axes=F, xlab=NA, ylab=NA,
             ylim=c(lim11,lim22),col.lab=3)
        abline(threshold1,0,col=3,lwd=2)
        text(2,threshold1*1.15,threshold1,col=3)
        abline(8000,0,col=3,lwd=2)
        text(2,8000*1.15,8000,col=3)
        axis(side = 4,col.axis=3)
        mtext(side = 4, line = 2, 'Netlen',col=3)
        
        plot(da$dummy,da[,id]*da[,id2],type='o',col=1,xlab="session n?",ylab="nlines x netlen",main="",
             ylim=c(lim1,lim22))
        abline(threshold1,0,col=1,lwd=2)
        text(3,threshold1*1.2,threshold1,col=1)
        abline(8000,0,col=1,lwd=2)
        text(2,8000*1.15,8000,col=1)
        
        
      }
    }
    pdf(paste(hndl.net_line,"all_years_nlines_X_netlen.pdf",sep=""))
    fun.plot.nlines.netlen(Data=Data.daily.original,YRS=FINYEAR.daily,lim1=0,lim2=10,
                           lim11=0,lim22=15000,threshold=5,threshold1=Net.max)
    dev.off()
    
    Ves.not.splitting=c("G 297","E 004","F 029","E 056","E 045","E 030","E 009","A 277","B 038","A 153",
                        "B 067","B 042","B 091","F 199")
    
    Ves.splitting_possibly_reporting_total.netlen=c("M 141","G 406","G 375","A 201","A 224","B 142","F 517",
                                                    "E 067","F 190","E 065","E 059","E 034","E 007","E 035",
                                                    "A 271","F 520","F 541","F 244","F 505")
    
    Ves.splitting_possibly_reporting_each.panel=c("A 012","A 035B","F 417","A 212","A 107","E 075","E 055",
                                                  "A 057B")
    
    pdf(paste(hndl.net_line,"not_splitters.pdf",sep=""))
    fun.plot.nlines.netlen(Data=subset(Data.daily.original,VESSEL%in%Ves.not.splitting),
            YRS=FINYEAR.daily,lim1=0,lim2=10,lim11=0,lim22=15000,threshold=5,threshold1=Net.max)
    dev.off()
    
    pdf(paste(hndl.net_line,"splitters_total.len.pdf",sep=""))
    fun.plot.nlines.netlen(Data=subset(Data.daily.original,VESSEL%in%Ves.splitting_possibly_reporting_total.netlen),
           YRS=FINYEAR.daily,lim1=0,lim2=10,lim11=0,lim22=15000,threshold=5,threshold1=Net.max)
    dev.off()
    
    pdf(paste(hndl.net_line,"splitters_partial.len.pdf",sep=""))
    fun.plot.nlines.netlen(Data=subset(Data.daily.original,VESSEL%in%Ves.splitting_possibly_reporting_each.panel),
           YRS=FINYEAR.daily,lim1=0,lim2=10,lim11=0,lim22=15000,threshold=5,threshold1=Net.max)
    dev.off()
    
  }

  #-Net length
  graphics.off()
  fun.plot.eff.var(Data.daily.original,Current.yr,"NETLEN",0,15000,Net.max,"YES")
  #Inspect.net.monthly=fun.inspect(Data.monthly.original,FINYEAR.monthly,"Monthly",100,.2,"NETLEN")
 
    #vessel visually identified as problematic
  Inspect.net.daily=fun.inspect(Data.daily.original,Current.yr,"Daily",100,.2,"NETLEN")
 
  
  #-BDAYS by session
  graphics.off()
  fun.plot.eff.var(Data.daily.original,Current.yr,"BDAYS",0,40,15,"YES")
    #vessell visually identified as problematic
  #check.vesl=c("F 517","E 045","B 142","F 417") 
  #MIN=1; MAX=30
  #Inspect.bday=fun.further.chk(Data.daily.original,Current.yr,"BDAYS")
  
  
  #-Fdays by month (shouldn't be > 30)
  graphics.off()
  fun.plot.fishing.days.month(Data=Data.daily.original,YRS=Current.yr,VAR="fdays",lim1=0,lim2=40,threshold=30)

  
  #-Hours
  graphics.off()
  fun.plot.eff.var(Data.daily.original,Current.yr,"HOURS",0,48,24,"YES")
    #vessell visually identified as problematic
  #check.vesl=c("E 035","B 142","F 517","G 406") 
  #MIN=1; MAX=24
  #Inspect.hours=fun.further.chk(Data.daily.original,Current.yr,"HOURS")
 
  
  #-Shots
  graphics.off()
  Table.Shots.year=table(Effort.daily$finyear,Effort.daily$shots)
  fun.plot.eff.var(Data.daily.original,Current.yr,"SHOTS",0,10,2,"YES")
    #vessell visually identified as problematic
  #check.vesl=c("A 035B") 
  #MIN=1; MAX=2
  #Inspect.shots=fun.further.chk(Data.daily.original,Current.yr,"SHOTS")
  #Check.shot=subset(Data.daily.original,FINYEAR==Current.yr & VESSEL%in%check.vesl) 
  #Check.shot=Check.shot[!duplicated(Check.shot$Same.return.SNo),match(c("date","VESSEL",
  #        "SNo","DSNo","TSNo","SHOTS"),names(Check.shot))]
  #Check.shot=subset(Check.shot, SHOTS==0 | SHOTS>2)
  
  
  #-Shots X Hours
  graphics.off()
  fun.plot.combo(Data.daily.original,Current.yr,"SHOTS","HOURS",0,48,26,"YES")
  # check.vesl=c("E 035","B 142","F 517") 
  # Check.shot.hours=subset(Data.daily.original,FINYEAR==Current.yr & VESSEL%in%check.vesl) 
  # Check.shot.hours=Check.shot.hours[!duplicated(Check.shot.hours$Same.return.SNo),
  #           match(c("date","VESSEL","SNo","DSNo","TSNo","SHOTS","HOURS"),names(Check.shot.hours))]
  # Check.shot.hours$Shots.hours=with(Check.shot.hours,SHOTS*HOURS)
  # Check.shot.hours=subset(Check.shot.hours, Shots.hours>26)
  
    #separate vessell visually identified as problematic
  #MIN=1; MAX=24
  #graphics.off()
  #Inspect.hours.shots=fun.further.chk.combo(Data.daily.original,Current.yr,"SHOTS","HOURS")
  
  
  #-Effort (km.gn.hours)
  graphics.off()
  Explr.ves=fun.plot.km.gn.hours(Data.daily.original,Current.yr,"SHOTS","HOURS","NETLEN","YES")
   
    
  #Export file for checking original returns by data entry girls
  #daily
  if(exists("Check.netlen"))if(nrow(Check.netlen)>0) write.csv(Check.netlen,paste(handle,"/Check.net.csv",sep=""),row.names=F)
  if(exists("Inspect.hours"))if( nrow(Inspect.hours)>0) write.csv(Inspect.hours,paste(handle,"/Check.hours.csv",sep=""),row.names=F)
  if(exists("Check.shot") )if( nrow(Check.shot)>0) write.csv(Check.shot,paste(handle,"/Check.shot.csv",sep=""),row.names=F)
  if(exists("Check.shot.hours"))if( nrow(Check.shot.hours)>0) write.csv(Check.shot.hours,paste(handle,"/Check.shot_hours.csv",sep=""),row.names=F)
  
  #monthly
  #write.csv(Inspect.net.monthly,file="Monthly/Monthly.Check.net.length.csv",row.names=F)
  
}
graphics.off()




#SECTION E. ---- EFFORT CORRECTIONS ----      

#E.1. Make annual vessel gillnet fishing effort param averages        #Rory's rule 1b

    #effort, monthly
Effort.monthly$VesselID=with(Effort.monthly,paste(FINYEAR,MONTH,BLOCKX,VESSEL))
Effort.monthly$BlockAveID=with(Effort.monthly,paste(FINYEAR,BLOCKX,VESSEL))
Effort.monthly$AnnualVesselAveID=with(Effort.monthly,paste(FINYEAR,VESSEL))
Effort.monthly$MonthlyID=with(Effort.monthly,paste(FINYEAR,MONTH))
Effort.monthly$MonthlyZoneID=with(Effort.monthly,paste(FINYEAR,MONTH,zone))
Effort.monthly$BlockID=with(Effort.monthly,paste(FINYEAR,BLOCKX))

    #effort, daily
Effort.daily$VesselID=with(Effort.daily,paste(finyear,month,blockx,vessel))
Effort.daily$BlockAveID=with(Effort.daily,paste(finyear,blockx,vessel))
Effort.daily$AnnualVesselAveID=with(Effort.daily,paste(finyear,vessel))
Effort.daily$MonthlyID=with(Effort.daily,paste(finyear,month))
Effort.daily$MonthlyZoneID=with(Effort.daily,paste(finyear,month,zone))
Effort.daily$BlockID=with(Effort.daily,paste(finyear,blockx))


Ref.Table.data=subset(Effort.daily,method=="GN" & netlen >100 & !(blockx%in%Estuaries),select=c(bdays,
                     hours,shots,netlen,finyear,vessel,VesselID,BlockAveID,AnnualVesselAveID,MonthlyID,MonthlyZoneID))
Ref.Table.data1=subset(Effort.monthly,METHOD=="GN" & NETLEN >100 & !(BLOCKX%in%Estuaries),select=c(BDAYS,
                 HOURS,SHOTS,NETLEN,FINYEAR,VESSEL,VesselID,BlockAveID,AnnualVesselAveID,MonthlyID,MonthlyZoneID))
names(Ref.Table.data)=names(Ref.Table.data1)
Ref.Table.data=rbind(Ref.Table.data1,Ref.Table.data)
rm(Ref.Table.data1)

        #remove nonsense effort variables values for calculating accurate means
Ref.Table.data=subset(Ref.Table.data,SHOTS<=3 & SHOTS>0)
Ref.Table.data=subset(Ref.Table.data,HOURS<=24 & HOURS>0)
Ref.Table.data=subset(Ref.Table.data,NETLEN<=Net.max & NETLEN>0)
Ref.Table.data=subset(Ref.Table.data,BDAYS>0)


        #correct effort table               #Rory's rule 2a                           
Mean.Eff.VesselID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~VesselID, data=Ref.Table.data,mean,na.rm=T)
Mean.Eff.BlockAveID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~BlockAveID, data=Ref.Table.data,mean,na.rm=T)
Mean.Eff.AnnualVesselAveID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~AnnualVesselAveID, data=Ref.Table.data,mean,na.rm=T)
Mean.Eff.MonthlyID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~MonthlyID, data=Ref.Table.data,mean,na.rm=T)
Mean.Eff.VesselAveID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~VESSEL, data=Ref.Table.data,mean,na.rm=T)
Mean.Eff.MonthlyZoneID=aggregate(cbind(BDAYS,HOURS,SHOTS,NETLEN)~MonthlyZoneID, data=Ref.Table.data,mean,na.rm=T)

mean.names=c("_BDAYS.m","_HOURS.m","_SHOTS.m","_NETLEN.m")
names(Mean.Eff.VesselID)=c(names(Mean.Eff.VesselID)[1],paste(names(Mean.Eff.VesselID)[1],mean.names,sep=""))
names(Mean.Eff.BlockAveID)=c(names(Mean.Eff.BlockAveID)[1],paste(names(Mean.Eff.BlockAveID)[1],mean.names,sep=""))
names(Mean.Eff.AnnualVesselAveID)=c(names(Mean.Eff.AnnualVesselAveID)[1],paste(names(Mean.Eff.AnnualVesselAveID)[1],mean.names,sep=""))
names(Mean.Eff.MonthlyID)=c(names(Mean.Eff.MonthlyID)[1],paste(names(Mean.Eff.MonthlyID)[1],mean.names,sep=""))
names(Mean.Eff.VesselAveID)=c(names(Mean.Eff.VesselAveID)[1],paste(names(Mean.Eff.VesselAveID)[1],mean.names,sep=""))
names(Mean.Eff.MonthlyZoneID)=c(names(Mean.Eff.MonthlyZoneID)[1],paste(names(Mean.Eff.MonthlyZoneID)[1],mean.names,sep=""))


#round to integer shot
Mean.Eff.AnnualVesselAveID$AnnualVesselAveID_SHOTS.m=round(Mean.Eff.AnnualVesselAveID$AnnualVesselAveID_SHOTS.m)
Mean.Eff.MonthlyZoneID$MonthlyZoneID_SHOTS.m=round(Mean.Eff.MonthlyZoneID$MonthlyZoneID_SHOTS.m) 

#round to integer bday
Mean.Eff.AnnualVesselAveID$AnnualVesselAveID_BDAYS.m=round(Mean.Eff.AnnualVesselAveID$AnnualVesselAveID_BDAYS.m)
Mean.Eff.MonthlyZoneID$MonthlyZoneID_BDAYS.m=round(Mean.Eff.MonthlyZoneID$MonthlyZoneID_BDAYS.m)




#table vessels by zone
fn.Ves.by.zone=function(YEARS)
{
  DATA=subset(Data.monthly,YEAR.c%in%YEARS)
  Vess.by.Zone=with(DATA,table(VESSEL,zone))
  Vess.by.Zone=as.data.frame(ifelse(Vess.by.Zone>0,1,0))
  Vess.by.Zone$Prop=rowSums(Vess.by.Zone)
  
  Prop.1.zone=sum(Vess.by.Zone$Prop<2)/nrow(Vess.by.Zone)
  Prop.2.zones=sum(Vess.by.Zone$Prop==2)/nrow(Vess.by.Zone) 
  Prop.3.zones=sum(Vess.by.Zone$Prop>2)/nrow(Vess.by.Zone)
  return(data.frame(Prop.1.zone,Prop.2.zones,Prop.3.zones))
}
Prop.Ves.Zon.All.yrs=fn.Ves.by.zone(unique(Data.monthly$YEAR.c))
#Prop.Ves.Zon.All.2006.2012=fn.Ves.by.zone(2006:2012)

#Check if FDAYS match BDAYS for same monthly return
# fun.match.Fdays=function(DATA)
# {
#   DATA=subset(DATA,FINYEAR%in%FINYEAR.monthly[1:30])
#   DATA=DATA[!duplicated(DATA$Same.return),]
#   DATA$Same.rec=with(DATA,paste(FINYEAR,MONTH,VESSEL,BLOCKX))
#   SAME=unique(DATA$Same.rec)
#   test=aggregate(DATA$BDAYS, by=list(DATA$Same.rec),sum)
#   names(test)=c("Same.rec","Sum of BDAYS")
#   test=test[order(test$Same.rec),]
#   id=match(test$Same.rec,DATA$Same.rec)
#   test$FDAYS=DATA$FDAYS[id]
#   test$SAME=ifelse(test$FDAYS==test$"Sum of BDAYS",1,0)
#   id=which(test$SAME==0)
#   a=test[id,]
#   plot(test$FDAYS,test$"Sum of BDAYS",ylab="",xlab="")
#   lines(0:30,0:30,lwd=2,col=2)
#   return(a)
# }
# 
# tiff(file="C:/Matias/Analyses/Catch and effort/Outputs/Figure 1.FDAYS vs SUM BDAYS.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
# par(mfcol=c(2,1),mar=c(3,4,1,4),oma=c(2.5,.1,.1,.1),las=1)
# 
# No.match.data.monthly=fun.match.Fdays(Ref.Table.data)
# legend("topleft",paste("TDGDLF, No matches =",nrow(No.match.data.monthly)),bty='n')
# 
# No.match.data.monthly.north=fun.match.Fdays(Data.monthly.north)
# legend("topleft",paste("NSF, No matches =",nrow(No.match.data.monthly.north)),bty='n')
# 
# mtext("Sum of BDAYS",side=2,line=-1.25,font=1,las=0,cex=1.7,outer=T)
# mtext("FDAYS",side=1,line=1,font=1,las=0,cex=1.7,outer=T)
# dev.off()
#rm(Ref.Table.data)



#E.2. Adjust invalid effort measures                 
#note: this includes 0 values, NAs, and too large values

#E.2.1 correct Bdays (replace 0, NA or >31 with mean value)                      #Rory's rule 2b
#monthly
Effort.monthly=Effort.monthly%>% left_join(Mean.Eff.AnnualVesselAveID,by="AnnualVesselAveID")%>%
                                 left_join(Mean.Eff.MonthlyZoneID,by="MonthlyZoneID")
#Effort.monthly=merge(Effort.monthly,Mean.Eff.AnnualVesselAveID,by="AnnualVesselAveID",all.x=T)
#Effort.monthly=merge(Effort.monthly,Mean.Eff.MonthlyZoneID,by="MonthlyZoneID",all.x=T)

Effort.monthly$BDAYS.c=Effort.monthly$BDAYS
Effort.monthly$BDAYS.c=with(Effort.monthly,
              ifelse((BDAYS==0|is.na(BDAYS)|BDAYS>31) & METHOD=="GN",AnnualVesselAveID_BDAYS.m,BDAYS.c))
Effort.monthly$BDAYS.c=with(Effort.monthly,
              ifelse((BDAYS.c==0|is.na(BDAYS.c)|BDAYS.c>31) & METHOD=="GN",MonthlyZoneID_BDAYS.m,BDAYS.c))
table(Effort.monthly$BDAYS.c,Effort.monthly$METHOD,useNA='ifany') #check correction (look at GN only)

#daily   
Effort.daily=Effort.daily%>%left_join(Mean.Eff.AnnualVesselAveID,by="AnnualVesselAveID")%>%
                            left_join(Mean.Eff.MonthlyZoneID,by="MonthlyZoneID")
#Effort.daily=merge(Effort.daily,Mean.Eff.AnnualVesselAveID,by="AnnualVesselAveID",all.x=T)
#Effort.daily=merge(Effort.daily,Mean.Eff.MonthlyZoneID,by="MonthlyZoneID",all.x=T)
Effort.daily$bdays.c=Effort.daily$bdays
Effort.daily$bdays.c=with(Effort.daily,
              ifelse((bdays==0|is.na(bdays)|bdays>31) & method=="GN",AnnualVesselAveID_BDAYS.m,bdays.c))
Effort.daily$bdays.c=with(Effort.daily,
              ifelse((bdays.c==0|is.na(bdays.c)|bdays.c>31) & method=="GN",MonthlyZoneID_BDAYS.m,bdays.c))
table(Effort.daily$bdays.c,Effort.daily$method,useNA='ifany') #check correction (look at GN only)


#create reporter
Effort.monthly$Eff.Reporter=with(Effort.monthly,ifelse(BDAYS==BDAYS.c,"good","bad")) 
Effort.daily$Eff.Reporter=with(Effort.daily,ifelse(bdays==bdays.c,"good","bad")) #create reporter


#E.2.2 correct Hours                                #Rory's rule 2c                      
Mean.hours=19           #(as per Simpfendorfer et al. 2000)

#monthly
Effort.monthly$HOURS.c=Effort.monthly$HOURS
Effort.monthly$HOURS.c=with(Effort.monthly,
              ifelse((HOURS==0|is.na(HOURS)|HOURS>24) & METHOD=="GN",AnnualVesselAveID_HOURS.m,HOURS.c))
Effort.monthly$HOURS.c=with(Effort.monthly,
              ifelse((HOURS.c==0|is.na(HOURS.c)|HOURS.c>24) & METHOD=="GN",MonthlyZoneID_HOURS.m,HOURS.c))
Effort.monthly$HOURS.c=with(Effort.monthly,
              ifelse((HOURS.c==0|is.na(HOURS.c)|HOURS.c>24) & METHOD=="GN",Mean.hours,HOURS.c))
table(Effort.monthly$HOURS.c,Effort.monthly$METHOD,useNA='ifany') #check correction (look at GN only)


#daily

#Inspect daily records of fishers fishing > 24 hours
#note: random inspection of daily records from fishers fishing >24 hours shows that these records are feasible
#       given that shots prior or post are < 24 hours long. hence, removed statement "|hours.c>24" from hours.c ammendment
insp.more.24.hours="NO"
if(insp.more.24.hours=="YES")
{
  a=subset(Effort.daily,hours>24 & method=="GN")
  ver.TSNo=as.character(unique(a$TSNo))
  ver=vector('list',length(ver.TSNo))
  for(v in 1:length(ver))
  {
    xx=subset(Effort.daily,TSNo==ver.TSNo[v],select=c(Same.return.SNo,SNo,DSNo,TSNo,vessel,date,hours,netlen))
    xx=xx[!duplicated(xx$Same.return.SNo),]
    xx=xx[order(xx$date),]
    xx$delta.day=0
    if(nrow(xx)>1)
    {
      for(x in 2:nrow(xx)) xx$delta.day[x]=xx$date[x]-xx$date[x-1]
      xx$cum.hours=cumsum(xx$hours)
    }
    ver[[v]]=xx
  }
}


#note: given that catch is reported for 0 hour shots, hours need to be imputed
Effort.daily$hours.c=Effort.daily$hours
Effort.daily$hours.c=with(Effort.daily,
              ifelse((hours==0|is.na(hours)) & method=="GN",AnnualVesselAveID_HOURS.m,hours.c))
Effort.daily$hours.c=with(Effort.daily,
              ifelse((hours.c==0|is.na(hours.c)) & method=="GN",MonthlyZoneID_HOURS.m,hours.c))
Effort.daily$hours.c=with(Effort.daily,
              ifelse((hours.c==0|is.na(hours.c)) & method=="GN",Mean.hours,hours.c))
table(Effort.daily$hours.c,Effort.daily$method,useNA='ifany') #check correction (look at GN only)

Effort.monthly$Eff.Reporter=with(Effort.monthly,ifelse(!(HOURS==HOURS.c),"bad",Eff.Reporter))
Effort.daily$Eff.Reporter=with(Effort.daily,ifelse(!(hours==hours.c),"bad",Eff.Reporter)) 


#E.2.3 correct Shots                                     #Rory's rule 2d             
#monthly
Effort.monthly$SHOTS.c=Effort.monthly$SHOTS
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(METHOD=="GN"& (SHOTS==0|is.na(SHOTS))&!(BLOCKX%in%Estuaries),
                          AnnualVesselAveID_SHOTS.m,SHOTS.c))
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(METHOD=="GN"& (SHOTS.c==0|is.na(SHOTS.c))&!(BLOCKX%in%Estuaries),
                          MonthlyZoneID_SHOTS.m,SHOTS.c))

      #correct some specific vessels for typos or splitting panels
Effort.monthly$SHOTS.c=with(Effort.monthly,
      ifelse(METHOD=="GN" & VESSEL=="A 069" & SHOTS%in%c(2:3,7:10) & NETLEN<1000,1,       
      ifelse(METHOD=="GN" & VESSEL=="A 080" & SHOTS%in%c(4,10),1,
      ifelse(METHOD=="GN" & VESSEL=="B 014" & SHOTS%in%c(10),1,
      ifelse(METHOD=="GN" & VESSEL=="D 026" & SHOTS>1,1,
      ifelse(METHOD=="GN" & VESSEL=="E 033" & SHOTS==9,1,
      ifelse(METHOD=="GN" & VESSEL=="F 232" & SHOTS==10 & NETLEN<500,1,
      ifelse(METHOD=="GN" & VESSEL=="F 589" & SHOTS==10 & NETLEN<500,1,
      ifelse(METHOD=="GN" & VESSEL=="M 008" & SHOTS==10,1,
      SHOTS.c)))))))))


        #if still NA default to 1
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(METHOD=="GN"& (SHOTS.c==0|is.na(SHOTS.c)),1,SHOTS.c))
table(Effort.monthly$SHOTS.c,Effort.monthly$METHOD,useNA='ifany')  #check correction (look at GN only)

Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(METHOD=="GN"& VESSEL=="E 033" & is.na(SHOTS.c),1,
                                                  ifelse(METHOD=="GN"& VESSEL=="F 232" & is.na(SHOTS.c),1,            
                                                         SHOTS.c)))


#daily
Effort.daily$shots.c=Effort.daily$shots
Effort.daily$shots.c=with(Effort.daily,ifelse(method=="GN"& (shots==0|is.na(shots))&!(blockx%in%Estuaries),
                          AnnualVesselAveID_SHOTS.m,shots.c))
Effort.daily$shots.c=with(Effort.daily,ifelse(method=="GN"& (shots.c==0|is.na(shots.c))&!(blockx%in%Estuaries),
                          MonthlyZoneID_SHOTS.m,shots.c))
      #if still NA default to 1
Effort.daily$shots.c=with(Effort.daily,ifelse(method=="GN"& (shots.c==0|is.na(shots.c)),1,shots.c))
table(Effort.daily$shots.c,Effort.daily$method,useNA='ifany')  #check correction (look for NAs at GN only)

    #correct some specific vessels for typos or splitting panels
Effort.daily$shots.c=with(Effort.daily,ifelse(method=="GN" & vessel=="E 059" & shots==11| is.na(shots),1,shots.c))


#E.2.3.1 correct high shots                               #Rory's rule 2e            
#monthly
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse((METHOD=="GN"& SHOTS.c %in% 3:6 & NETLEN>1001)|
  (METHOD=="GN"& SHOTS.c %in% 7:8 & NETLEN>601)|(METHOD=="GN"& SHOTS.c >8 & NETLEN>361),1,SHOTS.c))

#daily
Effort.daily$shots.c=with(Effort.daily,ifelse((method=="GN"& shots.c %in% 3:6 & netlen>1001)|
  (method=="GN"& shots.c %in% 7:8 & netlen>601)|(method=="GN"& shots.c >8 & netlen>361),1,shots.c))


#E.2.3.2 correct "B67" multiple shots                               #Rory's rule 2f 
#note: this amends multiple shots from this vessel for some sessions that are not fixed by 2e, based on
#      a discussion with vessel owner/skipper and observations on how he fished
#monthly
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(VESSEL=="B 067" & FINYEAR %in% c("2000-01","2001-02",
                          "2002-03","2003-04","2004-05","2005-06"),1,SHOTS.c))


Effort.monthly$Eff.Reporter=with(Effort.monthly,ifelse(!(SHOTS==SHOTS.c),"bad",Eff.Reporter))
Effort.daily$Eff.Reporter=with(Effort.daily,ifelse(!(shots==shots.c),"bad",Eff.Reporter))


#Correct nlines        
Effort.daily$nlines.c=Effort.daily$nlines
Effort.daily$nlines.c=with(Effort.daily,ifelse(vessel=="F 417" & nlines.c==11 & 
                   Same.return.SNo=="3 TDGLF8001838 TDGLF8001838",1,
                   ifelse(vessel=="E 059" & nlines.c>3,NA,nlines.c)))

#Add Rory's manual changes to netlen and nlines for Alex
#note: Rory manually identified records with discrepancies in the netlen and nlines combinations and
#     changed them accordingly. This bit of code updates those records with Rory's netlen and nlines values
Rory_Alex_net_val=subset(Rory_Alex_net_val,!Change.explanation=="",
          select=c(vessel,SessionID,fishery,finyear,netlen.original,netlen_Matias_adjust,netlen.Rory,
                    nlines.original,nlines_Matias_adjust,nlines.Rory,Change.explanation))
Rory_Alex_net_val$netlen.delta=abs(with(Rory_Alex_net_val,netlen.original-netlen.Rory))
Rory_Alex_net_val$nlines.delta=abs(with(Rory_Alex_net_val,nlines.original-nlines.Rory)) 
Fixed.Rory.Alex=subset(Rory_Alex_net_val,netlen.delta>0 | nlines.delta>0)
for(i in 1:ncol(Effort.daily))if(is.factor(Effort.daily[,i])) Effort.daily[,i]=as.character(Effort.daily[,i])
ddd=subset(Effort.daily,Same.return.SNo%in%unique(Fixed.Rory.Alex$SessionID))
Effort.daily=subset(Effort.daily,!Same.return.SNo%in%unique(Fixed.Rory.Alex$SessionID))
ddd=ddd%>%left_join(subset(Fixed.Rory.Alex,select=c(vessel,SessionID,netlen.Rory,nlines.Rory)),
                    by=c("vessel"="vessel","Same.return.SNo"="SessionID"))
# ddd=merge(ddd,subset(Fixed.Rory.Alex,select=c(vessel,SessionID,netlen.Rory,nlines.Rory)),
#           by.x=c("vessel","Same.return.SNo"),by.y=c("vessel","SessionID"),all.x=T)
ddd$netlen.c=ddd$netlen.Rory
ddd$nlines.c=ddd$nlines.Rory
ddd=ddd[,match(colnames(Effort.daily),colnames(ddd))]
Effort.daily=rbind(Effort.daily,ddd)
rm(ddd,Rory_Alex_net_val,Fixed.Rory.Alex)


#update Reporter
Effort.daily$Eff.Reporter=with(Effort.daily,ifelse(!(nlines==nlines.c),"bad",Eff.Reporter))

#E.2.4 Correct net length                                                           
#monthly
Effort.monthly$NETLEN.c=Effort.monthly$NETLEN
Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(METHOD %in% c("DL","HL","LL"),NA,NETLEN.c))
Effort.monthly$NETLEN.c=with(Effort.monthly,
        ifelse(METHOD=="GN"& (NETLEN==0|is.na(NETLEN))&!(BLOCKX%in%Estuaries),AnnualVesselAveID_NETLEN.m,NETLEN.c))
Effort.monthly$NETLEN.c=with(Effort.monthly,
        ifelse(METHOD=="GN"& (NETLEN.c==0|is.na(NETLEN.c))&!(BLOCKX%in%Estuaries),MonthlyZoneID_NETLEN.m,NETLEN.c))
Tab.mon.net.method=table(Effort.monthly$NETLEN.c,Effort.monthly$METHOD,useNA='ifany')

    #correct splitting panels
Effort.monthly$NETLEN.c=with(Effort.monthly,
       ifelse(METHOD=="GN" & VESSEL=="A 069" & SHOTS%in%c(2:3,7:10) & NETLEN<1000,NETLEN*SHOTS,
       ifelse(METHOD=="GN" & VESSEL=="A 080" & SHOTS%in%c(4),SHOTS*NETLEN,
       ifelse(METHOD=="GN" & VESSEL=="D 026" & SHOTS>1,SHOTS*NETLEN,
       ifelse(METHOD=="GN" & VESSEL=="E 033" & SHOTS==9,SHOTS*NETLEN,
       ifelse(METHOD=="GN" & VESSEL=="F 232" & SHOTS==10 & NETLEN<500,SHOTS*NETLEN,
       ifelse(METHOD=="GN" & VESSEL=="F 589" & SHOTS==10 & NETLEN<500,SHOTS*NETLEN,
       NETLEN.c)))))))  

Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(VESSEL=="E 033" & METHOD=="GN" & NETLEN>7800,7800,NETLEN.c))
Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(VESSEL=="E 033" & METHOD=="GN" & is.na(NETLEN.c),NETLEN,NETLEN.c))
Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(VESSEL=="F 232" & METHOD=="GN" & is.na(NETLEN.c),10*NETLEN,NETLEN.c))


#daily
Effort.daily$netlen.c=Effort.daily$netlen
Effort.daily$netlen.c=with(Effort.daily,ifelse(method %in% c("DL","HL","LL"),NA,netlen.c))
Effort.daily$netlen.c=with(Effort.daily,
            ifelse(method=="GN"& (netlen==0|is.na(netlen))&!(blockx%in%Estuaries),AnnualVesselAveID_NETLEN.m,netlen.c))
Effort.daily$netlen.c=with(Effort.daily,
            ifelse(method=="GN"& (netlen.c==0|is.na(netlen.c))&!(blockx%in%Estuaries),MonthlyZoneID_NETLEN.m,netlen.c))
Tab.daily.net.method=table(Effort.daily$netlen.c,Effort.daily$method,useNA='ifany')


#update Reporter
Effort.monthly$Eff.Reporter=with(Effort.monthly,ifelse(!(NETLEN==NETLEN.c),"bad",Eff.Reporter))
Effort.daily$Eff.Reporter=with(Effort.daily,ifelse(!(netlen==netlen.c),"bad",Eff.Reporter))


#E.2.5 Correct Km.hours.shots > 90,000
Correct.shot.90k="No"
if(Correct.shot.90k=="Yes")
{
  Effort.monthly$Km.hours.shots=with(Effort.monthly,NETLEN.c*HOURS.c*SHOTS.c)
  Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(Km.hours.shots<90000,SHOTS.c,1))
  Effort.monthly=Effort.monthly[,-match("Km.hours.shots",names(Effort.monthly))]
}




#E.2.6 Check ups 

  #Hours-Shots combination                                                 
#Monthly
MoreThan24=subset(Effort.monthly,METHOD=="GN")
MoreThan24=MoreThan24[!duplicated(MoreThan24$Same.return),]
MoreThan24$ShotsHour=MoreThan24$HOURS.c*MoreThan24$SHOTS.c
MoreThan24.monthly=subset(MoreThan24,is.na(ShotsHour)|ShotsHour>24,
        select=c(FINYEAR,MONTH,BLOCKX,VESSEL,HOURS,HOURS.c,SHOTS,SHOTS.c,ShotsHour,NETLEN,NETLEN.c))
MoreThan24.monthly$ShotsHourNet=MoreThan24.monthly$ShotsHour*MoreThan24.monthly$NETLEN.c

#daily
MoreThan24=subset(Effort.daily,method=="GN")
MoreThan24=MoreThan24[!duplicated(MoreThan24$Same.return.SNo),]
MoreThan24$ShotsHour=MoreThan24$hours.c*MoreThan24$shots.c
MoreThan24.daily=subset(MoreThan24,is.na(ShotsHour)|ShotsHour>24,
    select=c(finyear,month,blockx,vessel,hours,hours.c,shots,shots.c,ShotsHour,netlen,netlen.c))
MoreThan24.daily$ShotsHourNet=MoreThan24.daily$ShotsHour*MoreThan24.daily$netlen.c


  #Netlen-Shots combination                                                 
#Monthly
CHECK=subset(Effort.monthly,METHOD=="GN" & LAT<=(-26))
CHECK=CHECK[!duplicated(CHECK$Same.return),]
CHECK$NetlenShots=CHECK$NETLEN.c*CHECK$SHOTS.c
CHECK.monthly=subset(CHECK,is.na(NetlenShots)|NetlenShots>5000,
      select=c(FINYEAR,MONTH,BLOCKX,VESSEL,SHOTS,SHOTS.c,NETLEN,NETLEN.c,NetlenShots,HOURS,HOURS.c))
CHECK.monthly$ShotsHourNet=CHECK.monthly$HOURS.c*CHECK.monthly$NetlenShots

#daily
CHECK=subset(Effort.daily,method=="GN" & LAT>=(26))
CHECK=CHECK[!duplicated(CHECK$Same.return.SNo),]
CHECK$NetlenShots=CHECK$netlen.c*CHECK$shots.c
CHECK.daily=subset(CHECK,is.na(NetlenShots)|NetlenShots>5000,
      select=c(finyear,month,blockx,vessel,shots,shots.c,netlen,netlen.c,NetlenShots,hours,hours.c))
CHECK.daily$ShotsHourNet=CHECK.daily$hours.c*CHECK.daily$NetlenShots




#E.2.6  Inspect effort vars after corrections and final changes
#Inspect.Eff.vars="YES"
Inspect.Eff.vars="NO"
if(Inspect.Eff.vars=="YES")
{
  #Functions. Monthly
  Dummy.month=function(YRS,VAR,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.monthly, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return),]
    
    dropVes=which(table(Data$VESSEL)<4)
    Data=subset(Data,!VESSEL%in%names(dropVes))
    Data=Data[order(Data$VESSEL,Data$YEAR,Data$MONTH,Data$Same.return),]
    id=match(VAR,names(Data))
    Unik.ves=as.character(unique(Data$VESSEL))
    par(mfcol=c(3,3),mai=c(.4,.4,.2,.1))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Same.return)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab=VAR,ylim=c(lim1,lim2))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Same.return==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
        }      
      }
      
    }
  }
  
  NetShot.month=function(YRS,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.monthly, FINYEAR%in% YRS & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return),]
    
    dropVes=which(table(Data$VESSEL)<4)
    Data=subset(Data,!VESSEL%in%names(dropVes))
    Data=Data[order(Data$VESSEL,Data$YEAR,Data$MONTH,Data$Same.return),]
    Data$VAR=Data$NETLEN.c*Data$SHOTS.c
    id=match("VAR",names(Data))
    Unik.ves=as.character(unique(Data$VESSEL))
    
    par(mfcol=c(3,3),mai=c(.4,.4,.2,.1))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Same.return)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab="",ylim=c(lim1,lim2))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Same.return==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          points(datt$dummy,datt$NETLEN.c,col="grey40",type='b',pch=19)
        }      
      }
      legend("bottomleft",c("shots X netlen","netlen"),pch=c(1,19),col=c("red","grey40"),bty='n')
      
    }
    print(Unik.ves[1])
  }
  
  NetShot.month.trouble=function(VES,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.monthly, VESSEL%in% VES & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return),]
    
    dropVes=which(table(Data$VESSEL)<4)
    Data=subset(Data,!VESSEL%in%names(dropVes))
    Data=Data[order(Data$VESSEL,Data$YEAR,Data$MONTH,Data$Same.return),]
    Data$VAR=Data$NETLEN.c*Data$SHOTS.c
    id=match("VAR",names(Data))
    Unik.ves=as.character(unique(Data$VESSEL))
    YRS=unique(Data$FINYEAR)
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Same.return)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab="",ylim=c(lim1,lim2))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Same.return==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          points(datt$dummy,datt$NETLEN.c,col="grey40",type='b',pch=19)
        }      
      }
      legend("bottomleft",c("shots X netlen","netlen"),pch=c(1,19),col=c("red","grey40"),bty='n')
      
    }
  }
  
  
  #Functions. Daily
  Dummy.day=function(YRS,VAR,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.daily, finyear%in% YRS & LAT>=26 & LAT <=40 & method=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    id=match(VAR,names(Data))
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Unik.ves=as.character(unique(Data$vessel))
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$vessel,Data$year,Data$month,Data$Session),]
    par(mfcol=c(3,3),mai=c(.4,.4,.2,.1))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,vessel==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab=VAR,main=paste(YRS,"",Unik.ves[i]),ylim=c(lim1,lim2))
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
        }      
      }
      
    }
  }
  
  NetShot.day=function(YRS,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.daily, finyear%in% YRS & LAT>=26 & LAT <=40 & method=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    
    dropVes=which(table(Data$vessel)<4)
    Data=subset(Data,!vessel%in%names(dropVes))
    
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$vessel,Data$year,Data$month,Data$Session),]
    Data$VAR=Data$netlen.c*Data$shots.c
    id=match("VAR",names(Data))
    Unik.ves=as.character(unique(Data$vessel))
    
    par(mfcol=c(3,3),mai=c(.4,.4,.2,.1))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,vessel==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab="",ylim=c(lim1,lim2))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          points(datt$dummy,datt$netlen.c,col="grey40",type='b',pch=19)
        }      
      }
      legend("bottomleft",c("shots X netlen","netlen"),pch=c(1,19),col=c("red","grey40"),bty='n')
      
    }
    print(Unik.ves[1])
  }
  
  NetShot.day.trouble=function(VES,lim1,lim2,threshold,withCols)     
  {
    Data=subset(Effort.daily, vessel==VES & LAT>=26 & LAT <=40 & method=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    dropVes=which(table(Data$vessel)<4)
    Data=subset(Data,!vessel%in%names(dropVes))
    
    
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$vessel,Data$year,Data$month,Data$Session),]
    Data$VAR=Data$netlen.c*Data$shots.c
    id=match("VAR",names(Data))
    Unik.ves=as.character(unique(Data$vessel))
    YRS=unique(Data$finyear)
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,vessel==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab="",ylim=c(lim1,lim2))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          points(datt$dummy,datt$netlen.c,col="grey40",type='b',pch=19)
        }      
      }
      legend("bottomleft",c("shots X netlen","netlen"),pch=c(1,19),col=c("red","grey40"),bty='n')
      
    }
  }
  
  fn.check=function(ves,MAX,YLIM)
  {
    Data=subset(Effort.daily, vessel==ves & LAT>=26 & LAT <=40 & method=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    Data$hr.sht=Data$hours.c*Data$shots.c
    par(mfcol=c(3,1),mai=c(.8,.8,.01,.01))
    plot(Data$hours.c,xlab="",ylab="Hours",cex.lab=1.5)
    mtext(unique(Data$vessel),3,line=-2,cex=2)
    plot(Data$shots.c,xlab="",ylab="Shots",cex.lab=1.5)
    plot(Data$hr.sht,ylim=c(0,YLIM),xlab="Session",ylab="Hours X shots",cex.lab=1.5)
    abline(h=MAX,col=2,lwd=3)
    a=subset(Data,hr.sht>MAX,select=c(SNo,DSNo,TSNo,date,vessel,fdays,bdays,shots,shots.c,
                                      hours,hours.c,hr.sht,netlen,netlen.c))
    a=a[order(a$date),]
    return(a)
  }
  
  
  #Monthly records
  #general inspection by year
  Dummy.Myr=unique(Effort.monthly$FINYEAR)
  Dummy.month(Dummy.Myr[3],"NETLEN.c",0,15000,Net.max,"YES")
  Dummy.month(Dummy.Myr[1],"HOURS.c",0,48,24,"YES")
  Dummy.month(Dummy.Myr[1],"SHOTS.c",0,10,2,"YES")
  
  #general inspection by year for netlen X shot
  NetShot.month(Dummy.Myr[29],0,15000,Net.max,"YES")
  
  #inspect specific vessels for all years
  Trouble.net.shot=c("F 169","E 010","E 045","F 725","G 361","F 850","F 778","E 018",
                     "E 056","E 054","E 030","E 055","E 062","E 044","E 061","E 004",
                     "E 033","B 103","E 059","E 067","B 151","A 162","A 043","M 141")
  NetShot.month.trouble(Trouble.net.shot[1],0,15000,Net.max,"YES")
  
  
  
  #Daily records
  #general inspection by year
  Dummy.day(Dummy.yr[4],"netlen.c",0,15000,Net.max,"YES")
  Dummy.day(Dummy.yr[7],"hours.c",0,48,24,"YES")
  Dummy.day(Dummy.yr[7],"shots.c",0,10,2,"YES")
  
  #general inspection by year for netlen X shot
  NetShot.day(Dummy.yr[7],0,15000,Net.max,"YES")
  
  
  #inspect specific vessels for all years
  Trouble.net.shot=c("A 224","B 067","B 142","E 009","E 034","E 035",
                     "E 056","E 065","E 067","E 059","F 244","M 141","F 517")
  NetShot.day.trouble(Trouble.net.shot[13],0,15000,Net.max,"YES")
  
  
  
  #shots
  check.shots=c("A 224","E 034","E 067","E 056","E 009","M 141","B 067","E 059","E 035",
                "E 065","F 244","F 517","B 042","E 034","B 142","F 541","A 012")
  
  Store.dodgy=vector('list',length(check.shots))
  names(Store.dodgy)=check.shots
  for (i in 1:length(check.shots))Store.dodgy[[i]]=fn.check(check.shots[i],26,100)
  
}

#fix remaining errors in monthly
Effort.monthly$Shots.net=Effort.monthly$SHOTS.c*Effort.monthly$NETLEN.c
chk=subset(Effort.monthly,Shots.net>Net.max & METHOD=="GN",select=c(Same.return,VESSEL,SHOTS,SHOTS.c,
                                                                    NETLEN,NETLEN.c,BDAYS.c,Shots.net))
chk$Changed=with(chk, ifelse(SHOTS==SHOTS.c & NETLEN==NETLEN.c,"NO","YES"))
fix.this=unique(subset(chk,Changed=="YES" & Shots.net>19000))$VESSEL

a=subset(Effort.monthly,VESSEL%in%fix.this & METHOD=="GN")
Mode <- function(x)
  {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
aggregate(SHOTS~VESSEL, a, Mode)
aggregate(NETLEN~VESSEL, a, Mode)

Effort.monthly$Eff.Reporter=with(Effort.monthly,
              ifelse(VESSEL=="F 428" & SHOTS.c==40 & Shots.net>19000,"bad",
              ifelse(VESSEL=="B 038" & SHOTS.c==4 & Shots.net>19000,"bad",
              ifelse(VESSEL=="E 009" & SHOTS.c==2 & Shots.net>19000,"bad",Eff.Reporter))))

Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(VESSEL=="F 428" & SHOTS.c==40 & Shots.net>19000,4,SHOTS.c))
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(VESSEL=="B 038" & SHOTS.c==4 & Shots.net>19000,1,SHOTS.c))
Effort.monthly$SHOTS.c=with(Effort.monthly,ifelse(VESSEL=="E 009" & SHOTS.c==2 & Shots.net>19000,4,SHOTS.c))

Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(VESSEL=="F 428" & NETLEN.c>1000 & Shots.net>19000,100,NETLEN.c))
Effort.monthly$NETLEN.c=with(Effort.monthly,ifelse(VESSEL=="E 009" & NETLEN.c>1000 & Shots.net>19000,3800,NETLEN.c))

hist(Effort.monthly$SHOTS.c*Effort.monthly$NETLEN.c)
Effort.monthly=Effort.monthly[,-match("Shots.net",names(Effort.monthly))]


#E.2.7 Check the distribution of corrected effort variables   
    
cfac=function(x,breaks=NULL,int) #function for transforming continuous var into factor levels
{
  if(is.null(breaks)) breaks=unique(quantile(x,probs=seq(0,1,int),na.rm=T))
  x=cut(x,breaks,include.lowest=T,right=F)
  levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
            c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
  return(x)
} 
fn.check.eff.vars.c=function(MAIN,DATA,YRS,var1,var2,var3,var4)
{
  if(MAIN=="Monthly returns")DATA=subset(DATA,METHOD=="GN" & !(BLOCKX%in%Estuaries) & FINYEAR%in%YRS)
  if(MAIN=="Daily returns")DATA=subset(DATA,method=="GN" & !(blockx%in%Estuaries) & finyear%in%YRS)
  a=DATA[,match(var1,names(DATA))]
  aa=DATA[,match(var2,names(DATA))]
  aaa=DATA[,match(var3,names(DATA))]
  aaaa=DATA[,match(var4,names(DATA))]
  par(mfcol=c(2,2))
  hist(a,breaks=100,main="",ylab="",xlab="BDAYS.c")
  hist(aa,breaks=100,main="",ylab="",xlab="HOURS.c")
  hist(aaa,breaks=100,main="",ylab="",xlab="SHOTS.c")
  hist(aaaa,breaks=100,main="",ylab="",xlab="NETLEN.c")
  mtext("Numbers",2,-2,outer=T,las=3,cex=1.25)
  mtext(MAIN,3,-2,outer=T,las=1,cex=1.25)
  
  Table.BDAYS.c=table(cfac(a,NULL,.1))
  Table.HOURS.c=table(cfac(aa,NULL,.1))
  Table.SHOTS.c=table(cfac(aaa,NULL,.1))
  Table.NETLEN.c=table(cfac(aaaa,NULL,.1))
  return(list(Table.BDAYS.c=Table.BDAYS.c,Table.HOURS.c=Table.HOURS.c,Table.SHOTS.c=Table.SHOTS.c,
              Table.NETLEN.c=Table.NETLEN.c))
}
par=par.default
Check.effort.changes=fn.check.eff.vars.c("Monthly returns",Effort.monthly,unique(Effort.monthly$FINYEAR),
                      "BDAYS.c","HOURS.c","SHOTS.c","NETLEN.c")
Check.effort.changes=fn.check.eff.vars.c("Daily returns",Effort.daily,unique(Effort.daily$finyear),
                      "bdays.c","hours.c","shots.c","netlen.c")




#E.2.8 Output percentage of corrected effort variables
fn.compare.eff=function(var1,var2,var3,var4)
{
  a=subset(Effort.monthly,METHOD=="GN" & !BLOCKX%in%Estuaries & LAT<=(-26))
  a=a[!duplicated(a$Same.return),]
  b=subset(Effort.daily,method=="GN" & !blockx%in%Estuaries & LAT>=(26))
  b=b[!duplicated(b$ID),]
  v1=a[,match(var1,names(a))]
  v2=a[,match(var2,names(a))]
  v3=b[,match(var3,names(b))]
  v4=b[,match(var4,names(b))]  
  return((sum(v1==v2,na.rm=T)+sum(v3==v4,na.rm=T))/(length(v1)+length(v3)))
}
Per.Eff.wrong.bdays=100*(1-fn.compare.eff("BDAYS","BDAYS.c","bdays","bdays.c"))
Per.Eff.wrong.hours=100*(1-fn.compare.eff("HOURS","HOURS.c","hours","hours.c"))
Per.Eff.wrong.shots=100*(1-fn.compare.eff("SHOTS","SHOTS.c","shots","shots.c"))
Per.Eff.wrong.net=100*(1-fn.compare.eff("NETLEN","NETLEN.c","netlen","netlen.c"))

# Export percentage of GN records reapportioned by year
Table.Eff.Reporter=table(Effort.monthly$Eff.Reporter)
write.csv(round(100*Table.Eff.Reporter[1]/sum(Table.Eff.Reporter),1),"C:/Matias/Analyses/Catch and effort/Outputs/Paper/Percent.GN.effort.fixed.csv")

# Effort for Alex's ASL model
if(do.Alexs=="YES")
{
  Effort.alex=subset(Effort.monthly,YEAR>2005 & !FINYEAR=="2005-06" & METHOD=="GN" & LAT<=(-26))
  Effort.alex.daily= subset(Effort.daily,method=="GN")
  Effort.alex.daily=Effort.alex.daily[,-match(c("LAT","LONG","date","block10","year","year.c","month"),names(Effort.alex.daily))]
}


 

##################--- F. PROCEDURE SECTION ---##############
#SECTION F 1. ---- EXTRACT QUANTITIES ---- 

#1 Check for spatial expansion

#1.1. Annual
Expand.fun=function(DATA)
{
  dummy=1:NN.monthly
  for(i in 1:NN.monthly)
  {  dat=subset(DATA,as.character(FINYEAR)==FINYEAR.monthly[i])
     dummy[i]=length(unique(dat$BLOCKX))
  }
  return(dummy)
}
Effort1.fun=function(DATA)
{
  Num.vess.expan=1:NN.monthly
  Vessels.per.yr=Annual.av.eff.per.ves.GN=Annual.av.eff.per.ves.LL=Annual.records.per.ves.GN=
    Annual.records.per.ves.LL=Annual.av.eff.fleet.GN=Annual.av.eff.fleet.LL=vector('list',length=NN.monthly)
  
  for(i in 1:NN.monthly)
  {         
    dat=subset(DATA,as.character(FINYEAR)==FINYEAR.monthly[i])
    all.vess=as.character(sort(unique(dat$VESSEL)))
    Vessels.per.yr[[i]]=all.vess
    Num.vess.expan[i]=length(all.vess)
    
    #     #remove duplicated effort
    #     dat$Duplicate=with(dat,paste(BLOCKX,MONTH,METHOD,VESSEL))
    #     id=which(duplicated(dat$Duplicate))
    #     if(!length(id)==0)dat=dat[-id,]
    #     
    #     #separate hook from gillnet
    #     dat.GN=subset(dat,METHOD=="GN")
    #     dat.LL=subset(dat,METHOD=="LL")
    #     
    #     if(!nrow(dat.GN)==0)
    #     {
    #       #calculate average annual effort parameters by vessel and method
    #       avg.BDAYS.GN=aggregate(dat.GN['BDAYS'],dat.GN['VESSEL'],FUN=mean,na.rm=T)
    #       avg.HOURS.GN=aggregate(dat.GN['HOURS'],dat.GN['VESSEL'],FUN=mean,na.rm=T)
    #       avg.SHOTS.GN=aggregate(dat.GN['SHOTS'],dat.GN['VESSEL'],FUN=mean,na.rm=T)
    #       avg.NETLEN.GN=aggregate(dat.GN['NETLEN'],dat.GN['VESSEL'],FUN=mean,na.rm=T) 
    #       
    #       avg.GN=merge(merge(merge(avg.BDAYS.GN,avg.HOURS.GN,by="VESSEL"),avg.SHOTS.GN,by="VESSEL"),
    #                    avg.NETLEN.GN,by="VESSEL")  
    #       Annual.av.eff.per.ves.GN[[i]]=avg.GN
    #       
    #       #Determine if only 1 record per year per fisher
    #       Annual.records.per.ves.GN[[i]]=table(as.character(dat.GN$VESSEL))
    #       
    #       #calculate average annual effort parameters by fleet and method
    #       Annual.av.eff.fleet.GN[[i]]=cbind(BDAYS=mean(dat.GN$BDAYS,na.rm=T),HOURS=mean(dat.GN$HOURS,na.rm=T),
    #                                         SHOTS=mean(dat.GN$SHOTS,na.rm=T),NETLEN=mean(dat.GN$NETLEN,na.rm=T))
    #     }
    #     
    #     if(!nrow(dat.LL)==0)
    #     {
    #       avg.BDAYS.LL=aggregate(dat.LL['BDAYS'],dat.LL['VESSEL'],FUN=mean,na.rm=T)
    #       avg.HOURS.LL=aggregate(dat.LL['HOURS'],dat.LL['VESSEL'],FUN=mean,na.rm=T)
    #       avg.NETLEN.LL=aggregate(dat.LL['NETLEN'],dat.LL['VESSEL'],FUN=mean,na.rm=T)   
    #       avg.HOOKS.LL=aggregate(dat.LL['HOOKS'],dat.LL['VESSEL'],FUN=mean,na.rm=T)
    #       
    #       avg.LL=merge(merge(merge(avg.BDAYS.LL,avg.HOURS.LL,by="VESSEL"),avg.HOOKS.LL,by="VESSEL"),
    #                    avg.NETLEN.LL,by="VESSEL")    
    #       Annual.av.eff.per.ves.LL[[i]]=avg.LL 
    #       
    #       #Determine if only 1 record per year per fisher
    #       Annual.records.per.ves.LL[[i]]=table(as.character(dat.LL$VESSEL))
    #       
    #       #calculate average annual effort parameters by fleet and method
    #       Annual.av.eff.fleet.LL[[i]]=cbind(BDAYS=mean(dat.LL$BDAYS,na.rm=T),HOURS=mean(dat.LL$HOURS,na.rm=T),
    #                                         HOOKS=mean(dat.LL$HOOKS,na.rm=T),NETLEN=mean(dat.LL$NETLEN,na.rm=T))
    #     }
    
  }
  
  return(list(Ves.yr=Vessels.per.yr,N.ves.yr=Num.vess.expan,Av.eff.ves.GN=Annual.av.eff.per.ves.GN,
              Av.eff.ves.LL=Annual.av.eff.per.ves.LL,Records.ves.GN=Annual.records.per.ves.GN,
              Records.ves.LL=Annual.records.per.ves.LL,Av.eff.fleet.GN=Annual.av.eff.fleet.GN,
              Av.eff.fleet.LL=Annual.av.eff.fleet.LL))
}

#number of blocks
N.blocks=unique(Data.monthly$BLOCKX)
N.blocks=N.blocks[!N.blocks%in%Estuaries]
N.blocks=length(N.blocks)


#1.2. Proportion of blocks fished per month for each year
# Table.0.catch.block.month=matrix(nrow=NN.monthly,ncol=12)
# rownames(Table.0.catch.block.month)=FINYEAR.monthly
# colnames(Table.0.catch.block.month)=1:12
# 
# for (i in 1:NN.monthly) 
# {
#   datas=subset(Data.monthly,FINYEAR==FINYEAR.monthly[i])
#   a=aggregate(BLOCKX~MONTH,data=datas,unique)
#   store=vector(length=12)
#   for (j in 1:12)store[j]=length(unlist(a[j,2]))
#   Table.0.catch.block.month[i,]=round(store/N.blocks,3)
# }
# 


#2. Calculate gillnet effort 

#note:  RORY DOESN'T USE SHOTS FOR KM.GN.DAYS as this effort measure is only concerned on km of nets used in a day                             
#       Also, shots are not used DAILY for km gn hours because fishers report the total hours fished per session,
#       whether doing 1 or 2 shots per session
#       For monthly km gn hours, shots should not be used because some fishers report combined hours

use.shots="NO"
#use.shots="YES"

  #2.1. Define nominal effort hours and effort days
      #Monthly
        #km gn hours
if(use.shots=="YES")
{
  Effort.monthly$Km.Gillnet.Hours.inv=with(Effort.monthly,
        ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(HOURS*BDAYS*SHOTS*NETLEN)/1000,NA))  #raw
  Effort.monthly$Km.Gillnet.Hours.val=with(Effort.monthly,
        ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(HOURS.c*BDAYS.c*SHOTS.c*NETLEN.c)/1000,NA))  #corrected
}
if(use.shots=="NO")
{
  Effort.monthly$Km.Gillnet.Hours.inv=with(Effort.monthly,
        ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(HOURS*BDAYS*NETLEN)/1000,NA))  #raw
  Effort.monthly$Km.Gillnet.Hours.val=with(Effort.monthly,
        ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(HOURS.c*BDAYS.c*NETLEN.c)/1000,NA))  #corrected
}

        #km gn days
if(use.shots=="YES")
{
  Effort.monthly$Km.Gillnet.Days.inv=with(Effort.monthly,
                    ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(SHOTS*BDAYS*NETLEN)/1000,NA)) 
  Effort.monthly$Km.Gillnet.Days.val=with(Effort.monthly,
                    ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(SHOTS.c*BDAYS.c*NETLEN.c)/1000,NA))   
}
if(use.shots=="NO")
{
  Effort.monthly$Km.Gillnet.Days.inv=with(Effort.monthly,
                    ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(BDAYS*NETLEN)/1000,NA))  
  Effort.monthly$Km.Gillnet.Days.val=with(Effort.monthly,
                    ifelse(METHOD=="GN"& !(BLOCKX%in%Estuaries),(BDAYS.c*NETLEN.c)/1000,NA))   
}


      #Daily
          #km gn hours
Effort.daily$Km.Gillnet.Hours.inv=with(Effort.daily,
              ifelse(method=="GN" & !(blockx%in%Estuaries),(hours*netlen)/1000,NA))  #raw
Effort.daily$Km.Gillnet.Hours.val=with(Effort.daily,
               ifelse(method=="GN"& !(blockx%in%Estuaries),(hours.c*netlen.c)/1000,NA))  #corrected

          #km gn days                                                                
if(use.shots=="YES")
{
  Effort.daily$Km.Gillnet.Days.inv=with(Effort.daily,
                    ifelse(method=="GN"& !(blockx%in%Estuaries),(shots*netlen)/1000,NA))  
  Effort.daily$Km.Gillnet.Days.val=with(Effort.daily,
                    ifelse(method=="GN"& !(blockx%in%Estuaries),(shots.c*netlen.c)/1000,NA))  
}
if(use.shots=="NO")
{
  Effort.daily$Km.Gillnet.Days.inv=with(Effort.daily,
                    ifelse(method=="GN"& !(blockx%in%Estuaries),(netlen)/1000,NA))  
  Effort.daily$Km.Gillnet.Days.val=with(Effort.daily,
                    ifelse(method=="GN"& !(blockx%in%Estuaries),(netlen.c)/1000,NA))  
}


  #2.2. Apply effort manipulations

      #2.2.1 Add 5% to effort records prior to 1989                 
  #Monthly
Effort.monthly$Km.Gillnet.Hours.plus5=with(Effort.monthly,ifelse(YEAR.c<1989,
                                    Km.Gillnet.Hours.val*Inc.per,Km.Gillnet.Hours.val))
Effort.monthly$Km.Gillnet.Days.plus5=with(Effort.monthly,ifelse(YEAR.c<1989,
                                    Km.Gillnet.Days.val*Inc.per,Km.Gillnet.Days.val))

      #2.2.2 Add increase in fishing efficiency                         
#note: this is added to the cpue standardisation data, not the effort data reported
#       in SOFAR, etc, when summarising effort....

#Monthly
Effort.monthly$Km.Gillnet.Hours.c=Effort.monthly$Km.Gillnet.Hours.plus5
Effort.monthly$Km.Gillnet.Days.c=Effort.monthly$Km.Gillnet.Days.plus5
  

#Daily
#note: no need to apply 5% extra to daily records
Effort.daily$Km.Gillnet.Hours.c=Effort.daily$Km.Gillnet.Hours.val
Effort.daily$Km.Gillnet.Days.c=Effort.daily$Km.Gillnet.Days.val



# -- Look at proportional effort by mesh sizes           
# Mesh.Effort.zn=subset(Mesh.Effort,zone%in%c("West","Zone1","Zone2"))
# Mesh.Effort.all=aggregate(Km.Gillnet.Days.c~finyear,Mesh.Effort.zn,sum)

Mesh.size=subset(Mesh.size,mshigh>0)
    #km gn days                                                                
if(use.shots=="YES") Mesh.size$Km.Gillnet.Days.c=with(Mesh.size,(shots*netlen)/1000)
if(use.shots=="NO")  Mesh.size$Km.Gillnet.Days.c=with(Mesh.size,(netlen)/1000)

Mesh.size$mshigh=round(Mesh.size$mshigh)
Mesh.size$LatDeg=-Mesh.size$LatDeg
Mesh.size$zone=as.character(Mesh.size$zone)
Mesh.size$zone=with(Mesh.size,ifelse(zone=="*","West",
                 ifelse(zone=="1","Zone1",
                 ifelse(zone=="2","Zone2",zone))))
Mesh.size$zone=as.character(with(Mesh.size,ifelse(zone=="West" & LatDeg>(-33) & LatDeg<=(-26) & LongDeg<116.5,"West",
                 ifelse(zone=="West" & LatDeg>(-26) & LongDeg<114,"Closed",
                 ifelse(zone=="West" & LatDeg>(-23) & LongDeg>=114 & LongDeg<123.75,"North",
                 ifelse(zone=="West" & LatDeg>(-23) & LongDeg>=123.75,"Joint",zone))))))


Mesh.Effort=aggregate(Km.Gillnet.Days.c~date+vessel+mshigh+zone+finyear,data=Mesh.size,max,na.rm=T) #remove duplicates

Mesh.Effort.zn=aggregate(Km.Gillnet.Days.c~zone+mshigh+finyear,data=Mesh.Effort,sum,na.rm=T)  #sum
Mesh.Effort=aggregate(Km.Gillnet.Days.c~mshigh+finyear,data=Mesh.Effort,sum,na.rm=T)


fn.mesh=function(what,agg)
{
  what=subset(what,finyear%in%FINYEAR.monthly)
  names(what)[1]='mesh'
  if(agg=="YES")
  {
    what$mesh=with(what,ifelse(mesh<165,'156-164',
          ifelse(mesh>165 &mesh<178,'166-177',ifelse(mesh>178,'>178',mesh))))
    what$mesh=factor(what$mesh,levels=c("156-164","165","166-177","178",">178"))
    what=aggregate(Km.Gillnet.Days.c~finyear+mesh,what,sum)
  }
  d=reshape(what,v.names="Km.Gillnet.Days.c",idvar="finyear",timevar="mesh",direction="wide")
  d=d[order(d$finyear),]
  d[is.na(d)]=0
  N=2:ncol(d)
  colnames(d)[N]=substr(colnames(d)[N],19,30)
  d[,N]=d[,N]/rowSums(d[,N],na.rm=T)
  Cls=ClS[match(colnames(d)[2:ncol(d)],names(ClS))]
  barplot(t(as.matrix(d[,N])),names.arg=d$finyear,ylim=c(0,1.195),col=Cls,cex.axis=1.25,
          cex.names=1.5)
  return(d)
}
get.eff.pro.mesh=function(a)
{
  a=a[,match(c("finyear",inch_6.5,inch_7),colnames(a))]
  a[,2:3]=a[,2:3]/rowSums(a[,2:3])
  return(a)
}

hnl.msh="C:/Matias/Analyses/Catch and effort/"

#Overall
AGG="YES"   #aggregate meshes
All.meshes=sort(unique(round(Mesh.size$mshigh)))
inch_6.5=165
inch_7=178
ClS=grey.colors(length(All.meshes),start=0.2,end=0.95)
names(ClS)=All.meshes
ClS[match(c(inch_6.5,inch_7),names(ClS))]=c("brown","darkolivegreen4")
if(AGG=="YES")
{
  ClS=c("#5A5A5A","brown","#979797","darkolivegreen4","#E9E9E9")
  names(ClS)=c("156-164","165","166-177","178",">178")
}

tiff(file=paste(hnl.msh,"mesh.proportional.effort.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mai=c(.8,.8,.01,.01),oma=c(1,1,.2,.1),las=1,xpd=T,mgp=c(.65,1,0))
a=fn.mesh(what=Mesh.Effort,agg=AGG)
mtext("Proportional annual effort",2,las=3,cex=1.75,line=2.7)
mtext("Financial year",1,cex=1.75,line=-1,outer=T)
if(AGG=="YES")
{
  legend(-1.15,1.15,legend = names(a)[-1],bty='n',horiz=T,fill = ClS,cex=1.2)
  
}else
{
  legend(-1.15,1.15,legend = names(ClS)[1:10],bty='n',horiz=T,fill = ClS[1:10],cex=1.2)
  legend(-1.15,1.09,legend = names(ClS)[10:(length(ClS)-1)],bty='n',horiz=T,fill=ClS[10:(length(ClS)-1)],cex=1.2)
}
dev.off()
NMSa=names(a)[-1]
d=get.eff.pro.mesh(a)
write.csv(d,paste(hnl.msh,"mesh.proportional.effort.csv",sep=""),row.names=F)


#by zone
jpeg(file=paste(hnl.msh,"mesh.proportional.effort.zone.jpeg",sep=""),width=2400,height=2400,units="px",res=300)
par(mfcol=c(3,1),mai=c(.1,.5,.15,.01),oma=c(4,2,1.5,.1),las=1,xpd=T,mgp=c(.65,1,0))
ZONE="West"
a=fn.mesh(what=subset(Mesh.Effort.zn,zone==ZONE,select=c(mshigh, finyear, Km.Gillnet.Days.c)),agg=AGG)
mtext("Proportional annual effort",2,las=3,cex=1.75,line=-.85,outer=T)
mtext("Financial year",1,cex=1.75,line=2,outer=T)
nnn=13
if(AGG=="YES")
{
  legend(-1.15,1.15,legend = NMSa,bty='n',horiz=T,fill = ClS,cex=1.2)
}else
{
  legend(-1.15,1.275,legend = names(ClS)[1:nnn],bty='n',horiz=T,fill = ClS[1:nnn],cex=1.25)
  legend(-1.15,1.15,legend = names(ClS)[(nnn+1):(length(ClS)-1)],bty='n',horiz=T,fill=ClS[(nnn+1):(length(ClS)-1)],cex=1.25)
}
par(font=2)
legend("topright",legend =ZONE ,bty='n',cex=1.75)
d=get.eff.pro.mesh(a)
write.csv(d,paste(hnl.msh,"mesh.proportional.effort.",ZONE,".csv",sep=""),row.names=F)

ZONE="Zone1"
a=fn.mesh(what=subset(Mesh.Effort.zn,zone==ZONE,select=c(mshigh, finyear, Km.Gillnet.Days.c)),agg=AGG)
legend("topright",legend =ZONE ,bty='n',cex=1.75)
d=get.eff.pro.mesh(a)
write.csv(d,paste(hnl.msh,"mesh.proportional.effort.",ZONE,".csv",sep=""),row.names=F)

ZONE="Zone2"
a=fn.mesh(what=subset(Mesh.Effort.zn,zone==ZONE,select=c(mshigh, finyear, Km.Gillnet.Days.c)),agg=AGG)
legend("topright",legend =ZONE ,bty='n',cex=1.75)
d=get.eff.pro.mesh(a)
write.csv(d,paste(hnl.msh,"mesh.proportional.effort.",ZONE,".csv",sep=""),row.names=F)
dev.off()


      #2.2.3 Check if effort is calculated for all records and flesh out the NAs
fn.chec=function(what,what2)
{
  if(what=="monthly")dat=subset(Effort.monthly,METHOD=="GN" & LAT<=(-26) & !BLOCKX%in%Estuaries)
  if(what=="daily")dat=subset(Effort.daily,method=="GN" & LAT>=(26) & !blockx%in%Estuaries)
  id=match(what2,names(dat))
  cat("all these are NA =",sum(is.na(dat[,id])))
  return(dat[is.na(dat[,id]),])
}

vec.eff=c("Km.Gillnet.Days.c","Km.Gillnet.Hours.c")
#vec.eff=c("Km.Gillnet.Days.inv","Km.Gillnet.Hours.inv","Km.Gillnet.Days.c","Km.Gillnet.Hours.c")
Missing.monthly=vector('list',length(vec.eff))
names(Missing.monthly)=vec.eff
Missing.daily=Missing.monthly
for(j in 1:length(vec.eff))
{
    Missing.monthly[[j]]=fn.chec("monthly",vec.eff[j])
    Missing.daily[[j]]=fn.chec("daily",vec.eff[j])

  }
#No NAs in corrected effort variables, only in the raw ones.


  #Post effort calculation checks by vessel. Any Inconsistencies?
Do.post.check="NO" 
#Do.post.check="YES"
if(Do.post.check=="YES")
{
  #monthly
  Post.Effort.month=function(VES,VAR,threshold,withCols)     
  {
    Data=subset(Effort.monthly, VESSEL%in% VES & LAT<=(-26) & LAT >=(-40) & METHOD=="GN")
    Data=Data[!duplicated(Data$Same.return),]
    
    Data=Data[order(Data$VESSEL,Data$YEAR,Data$MONTH,Data$Same.return),]
    
    id=match(VAR,names(Data))
    Unik.ves=as.character(unique(Data$VESSEL))
    YRS=unique(Data$FINYEAR)
    
    Data$Km.Gillnet.Days.val=with(Data,ifelse(METHOD=="GN",(SHOTS.c*BDAYS.c*NETLEN.c)/1000,NA))   
    Data$Km.Gillnet.Days.plus5=with(Data,ifelse(YEAR.c<1990,
                      Km.Gillnet.Days.val*Inc.per,Km.Gillnet.Days.val))
    Data$Km.Gillnet.Days.c.shot=with(Data,ifelse(YEAR.c<1995,
          Km.Gillnet.Days.plus5*(1+Fish.Pow*(1+YEAR.c-YEAR.c.monthly[1])),Km.Gillnet.Days.plus5))
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      Unik.sess=unique(da$Same.return)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab=VAR,ylim=c(0,max(da[,id],da$Km.Gillnet.Days.c.shot)))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Same.return==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          if(VAR=="Km.Gillnet.Days.c")points(datt$dummy,datt$Km.Gillnet.Days.c.shot,col="grey40",cex=0.75,pch=19,type='b',lty=3,lwd=.8)
          
        }      
      }
      if(VAR=="Km.Gillnet.Days.c")legend("bottomleft",c("without shot","with shot"),pch=c(1,19),col=c("red","grey40"),bty='n')
      
    }
  }
  Unic.Ves=with(subset(Effort.monthly,LAT<=(-26) & LAT >=(-40) & METHOD=="GN"),table(VESSEL))
  Unic.Ves=which(Unic.Ves>4)
  Unic.Ves=names(Unic.Ves)
  par(mfcol=c(6,1),mai=c(.25,.8,.1,.1))
  NNN=length(Unic.Ves)
  for (p in 300:NNN)Post.Effort.month(Unic.Ves[p],"Km.Gillnet.Days.c",50,"YES")
  for (p in 1:NNN)Post.Effort.month(Unic.Ves[p],"Km.Gillnet.Hours.c",50,"YES")
  
  troubled.ones=c("A 062", "A 043","B 098","B 132","E 004","E 007",
                  "E 010", "F 148","E 067","E 061","E 062","E 054","E 056",
                  "E 055","P 119")
  for (p in 1:length(troubled.ones))Post.Effort.month(troubled.ones[p],"Km.Gillnet.Days.c",50,"YES")
  
  #daily
  Post.Effort.day=function(VES,VAR,threshold,withCols)     
  {
    Data=subset(Effort.daily, vessel==VES & LAT>=26 & LAT <=40 & method=="GN")
    Data=Data[!duplicated(Data$Same.return.SNo),]
    Data$Session_ID=with(Data,paste(TSNo,"_",DSNo,"_",SNo,sep=""))  
    Data=Data[!(duplicated(Data$Session_ID)),]
    Data$Session=with(Data,paste(TSNo,"_",DSNo,sep=""))
    Data=Data[order(Data$vessel,Data$year,Data$month,Data$Session),]
    
    id=match(VAR,names(Data))
    Unik.ves=as.character(unique(Data$vessel))
    YRS=unique(Data$finyear)
    
    Data$Km.Gillnet.Days.val=with(Data,ifelse(method=="GN",(shots.c*netlen.c)/1000,NA))  
    Data$Km.Gillnet.Days.c.shot=Data$Km.Gillnet.Days.val
    
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,vessel==Unik.ves[i])
      Unik.sess=unique(da$Session)
      NN=length(Unik.sess)
      da$dummy=1:nrow(da)
      
      
      plot(da$dummy,da[,id],type="b",xlab="session n?",ylab=VAR,ylim=c(0,max(da[,id],da$Km.Gillnet.Days.c.shot)))
      legend("topright",c(Unik.ves[i],YRS),bty='n')
      abline(threshold,0,col=2,lwd=2)
      text(1,threshold*1.15,threshold,col=2)
      Mean=mean(da[,id],na.rm=T)
      abline(h=Mean,col=3,lwd=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*0.7,col=3,lwd=1,lty=2)
      abline(h=Mean*1.3,col=3,lwd=1,lty=2)
      if(withCols=="YES")
      {
        Cols=sample(rainbow(NN))
        for (j in 1:NN)
        {
          datt=subset(da,Session==Unik.sess[j])
          points(datt$dummy,datt[,id],col=Cols[j],type='b')
          if(VAR=="Km.Gillnet.Days.c")points(datt$dummy,datt$Km.Gillnet.Days.c.shot,col="grey40",cex=0.75,pch=19,type='b',lty=3,lwd=.8)
        }      
      }
      if(VAR=="Km.Gillnet.Days.c")legend("bottomleft",c("without shot","with shot"),pch=c(1,19),col=c("red","grey40"),bty='n')
    }
  }
  Unic.Ves=with(subset(Effort.daily,LAT>=26 & LAT <=40 & method=="GN"),unique(vessel))
  par(mfcol=c(6,1),mai=c(.25,.8,.1,.1))
  for (p in 1:length(Unic.Ves))Post.Effort.day(Unic.Ves[p],"Km.Gillnet.Days.c",50,"YES")
  for (p in 1:length(Unic.Ves))Post.Effort.day(Unic.Ves[p],"Km.Gillnet.Hours.c",50,"YES")
  
  troubled.ones=c("A 224","B 067","F 517","B 142","E 035","M 141","E 034","E 065")
  Post.Effort.day(troubled.ones[8],"Km.Gillnet.Days.c",50,"YES")
  
}


  #2.3. Aggregate effort by year and zone for reporting in SOFAR                         

#2.3.1. Daily

Effort.daily$Fishing_yr=with(Effort.daily,ifelse(month>=6,year,year-1))

        #km gn days 
if(Use.Date=="NO")
{
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~ID+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
  
  Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+Fishing_yr+block10,data=Effort.daily,max,na.rm=T)
}

if(Use.Date=="YES")
{
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~date+vessel+zone+finyear,data=Effort.daily,max,na.rm=T) 
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~date+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
  
  Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~date+vessel+zone+Fishing_yr+block10,data=Effort.daily,max,na.rm=T)
}

Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~zone+finyear,data=Attach.Effort.daily.c,sum,na.rm=T)
Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~zone+finyear,data=Attach.Effort.daily,sum,na.rm=T)
Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~vessel+zone+Fishing_yr+block10,data=Jodie.Effort.daily.c_block10,sum,na.rm=T)


        #km gn hours
Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~finyear+vessel+ID+zone,data=Effort.daily,max,na.rm=T)
Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~finyear+vessel+ID+zone,data=Effort.daily,max,na.rm=T)

Jodie.Effort.daily.hrs.c_block10=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+ID+zone+block10,data=Effort.daily,max,na.rm=T)

# if(Use.Date=="YES")
# {
#   Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
#   Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
# }

Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~zone+finyear,data=Attach.Effort.daily.hrs.c,sum,na.rm=T)
Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~zone+finyear,data=Attach.Effort.daily.hrs,sum,na.rm=T)

Jodie.Effort.daily.hrs.c_block10=aggregate(Km.Gillnet.Hours.c~zone+Fishing_yr+vessel+block10,data=Jodie.Effort.daily.hrs.c_block10,sum,na.rm=T)

        #merge into single file
Attach.Effort.daily.c=Attach.Effort.daily.c %>% left_join(Attach.Effort.daily.hrs.c,by=c("zone","finyear"))
Attach.Effort.daily=Attach.Effort.daily %>% left_join(Attach.Effort.daily.hrs,by=c("zone","finyear"))
Attach.Effort.daily=Attach.Effort.daily.c %>% left_join(Attach.Effort.daily,by=c("zone","finyear"))
#Attach.Effort.daily.c=merge(Attach.Effort.daily.c,Attach.Effort.daily.hrs.c,by=c("zone","finyear"),all=T)
#Attach.Effort.daily=merge(Attach.Effort.daily,Attach.Effort.daily.hrs,by=c("zone","finyear"),all=T)
#Attach.Effort.daily=merge(Attach.Effort.daily.c,Attach.Effort.daily,by=c("zone","finyear"),all.x=T)


Jodie.Effort=Jodie.Effort.daily.c_block10%>%full_join(Jodie.Effort.daily.hrs.c_block10,
                                                by=c("vessel","zone","Fishing_yr","block10"))
#Jodie.Effort=merge(Jodie.Effort.daily.c_block10,Jodie.Effort.daily.hrs.c_block10,
#                   by=c("vessel","zone","Fishing_yr","block10"),all.x=T)



#2.3.2. Monthly 

#Rory's Zone 2 multiplier to account for multiple shots per day
Apply.zn2.multi.sht="No"
if(Use.Date=="YES") Apply.zn2.multi.sht="Yes"

if(Apply.zn2.multi.sht=="Yes")
{
  Effort.monthly$Km.Gillnet.Days.inv=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.inv*1.6,Km.Gillnet.Days.inv))
  Effort.monthly$Km.Gillnet.Days.val=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.val*1.6,Km.Gillnet.Days.val))
  Effort.monthly$Km.Gillnet.Days.plus5=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.plus5*1.6,Km.Gillnet.Days.plus5))
  Effort.monthly$Km.Gillnet.Days.c=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.c*1.6,Km.Gillnet.Days.c))
}


#note: max is used to remove duplicated effort records from same return (i.e. from different species caught in that return)
    #km gn days             
      #max to remove duplicates
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)

      #split boundary blocks
fn.split.boundary=function(DAT,What)
{
  id=match(What,colnames(DAT))
  a=subset(DAT,BLOCKX%in%Boundary.blk)
  b=subset(DAT,!BLOCKX%in%Boundary.blk)  
  x=a
  x$zone="Zone2"
  x[,id]=NA
  x[,id]=ifelse(!x$FINYEAR%in%c("2003-04","2004-05","2005-06"),0.5*a[,id],0.05*a[,id])
  a[,id]=ifelse(!a$FINYEAR%in%c("2003-04","2004-05","2005-06"),0.5*a[,id],0.95*a[,id])
  a=rbind(a,x)
  return(rbind(a,b))
}
Attach.Effort.monthly.c=fn.split.boundary(Attach.Effort.monthly.c,"Km.Gillnet.Days.c")
Attach.Effort.monthly=fn.split.boundary(Attach.Effort.monthly,"Km.Gillnet.Days.inv")

      #sum to aggregage by zone and year
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~zone+FINYEAR,data=Attach.Effort.monthly.c,sum,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~zone+FINYEAR,data=Attach.Effort.monthly,sum,na.rm=T)


    #km gn hours
      #max to remove duplicates
Attach.Effort.monthly.hrs.c=aggregate(Km.Gillnet.Hours.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)
Attach.Effort.monthly.hrs=aggregate(Km.Gillnet.Hours.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)


      #split boundary blocks
Attach.Effort.monthly.hrs.c=fn.split.boundary(Attach.Effort.monthly.hrs.c,"Km.Gillnet.Hours.c")
Attach.Effort.monthly.hrs=fn.split.boundary(Attach.Effort.monthly.hrs,"Km.Gillnet.Hours.inv")

      #sum to aggregage by zone and year
Attach.Effort.monthly.hrs.c=aggregate(Km.Gillnet.Hours.c~zone+FINYEAR,data=Attach.Effort.monthly.hrs.c,sum,na.rm=T)
Attach.Effort.monthly.hrs=aggregate(Km.Gillnet.Hours.inv~zone+FINYEAR,data=Attach.Effort.monthly.hrs,sum,na.rm=T)


  #merge into single file
Attach.Effort.monthly.c=Attach.Effort.monthly.c %>%left_join(Attach.Effort.monthly.hrs.c,by=c("zone","FINYEAR"))
Attach.Effort.monthly=Attach.Effort.monthly %>%left_join(Attach.Effort.monthly.hrs,by=c("zone","FINYEAR"))
Attach.Effort.monthly=Attach.Effort.monthly.c %>%left_join(Attach.Effort.monthly,by=c("zone","FINYEAR"))
#Attach.Effort.monthly.c=merge(Attach.Effort.monthly.c,Attach.Effort.monthly.hrs.c,by=c("zone","FINYEAR"),all=T)
#Attach.Effort.monthly=merge(Attach.Effort.monthly,Attach.Effort.monthly.hrs,by=c("zone","FINYEAR"),all=T)
#Attach.Effort.monthly=merge(Attach.Effort.monthly.c,Attach.Effort.monthly,by=c("zone","FINYEAR"),all.x=T)

rm(Attach.Effort.monthly.c,Attach.Effort.monthly.hrs.c,Attach.Effort.monthly.hrs)    
   

  #2.3. Merge Monthly and daily
names(Attach.Effort.daily)[2]="FINYEAR"

#add daily records reported in CAESS post 2006
Daily.l.years=sort(unique(Data.daily$FINYEAR))

a=subset(Attach.Effort.monthly,FINYEAR%in%Attach.Effort.daily$FINYEAR)
Attach.Effort.monthly=subset(Attach.Effort.monthly,!FINYEAR%in%Daily.l.years)
colnames(a)[3:ncol(a)]=paste(colnames(a)[3:ncol(a)],"1",sep=".")
Attach.Effort.daily=Attach.Effort.daily%>%left_join(a,by=c("zone","FINYEAR"))
# Attach.Effort.daily=merge(Attach.Effort.daily,a,by=c("zone","FINYEAR"),all.x=T)
Attach.Effort.daily[is.na(Attach.Effort.daily)]=0
Attach.Effort.daily$Km.Gillnet.Days.c=with(Attach.Effort.daily,Km.Gillnet.Days.c+Km.Gillnet.Days.c.1)
Attach.Effort.daily$Km.Gillnet.Hours.c=with(Attach.Effort.daily,Km.Gillnet.Hours.c+Km.Gillnet.Hours.c.1)
Attach.Effort.daily$Km.Gillnet.Days.inv=with(Attach.Effort.daily,Km.Gillnet.Days.inv+Km.Gillnet.Days.inv.1)
Attach.Effort.daily$Km.Gillnet.Hours.inv=with(Attach.Effort.daily,Km.Gillnet.Hours.inv+Km.Gillnet.Hours.inv.1)
Attach.Effort.daily=Attach.Effort.daily[,-match(c("Km.Gillnet.Days.c.1","Km.Gillnet.Hours.c.1",
                  "Km.Gillnet.Days.inv.1","Km.Gillnet.Hours.inv.1"),names(Attach.Effort.daily))]

Attach.Effort=rbind(Attach.Effort.monthly,Attach.Effort.daily)

#3. Aggregate gillnet effort for reporting in Sofar           
if(Use.Date=="NO") use.ID="YES"   #select if aggregating effort by SNo and DSNo
if(Use.Date=="YES") use.ID="NO"   #select if aggregating effort by date

  #3.1. Daily records
fn.Eff.Sofar=function(THESE.YRS)
{
  #3.1.1. Aggregate effort by year and zone
  SUBSET=subset(Effort.daily,finyear%in%THESE.YRS)
  
  #km gn days                         
  if(use.shots=="NO")
  {
    if(use.ID=="YES")
    {
      Km.gn.days=aggregate(netlen.c~ID+vessel+zone+finyear,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c~zone+finyear,data=Km.gn.days,sum,na.rm=T)
    }
    if(use.ID=="NO")
    {
      Km.gn.days=aggregate(netlen.c~finyear+vessel+date+zone,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c~finyear+zone,data=Km.gn.days,sum,na.rm=T)
    }
    
    Effort.Km.gn.days <- reshape(Km.gn.days,v.names="netlen.c",timevar="zone",idvar=c("finyear"),direction="wide")
    iid=match(c("netlen.c.West","netlen.c.Zone1","netlen.c.Zone2"),names(Effort.Km.gn.days))
    Effort.Km.gn.days[,iid]=Effort.Km.gn.days[,iid]/1e6
    names(Effort.Km.gn.days)[iid]=c("West","Zone1","Zone2")  
  }
  
  if(use.shots=="YES")
  {
    if(use.ID=="YES")
    {
      Km.gn.days=aggregate(netlen.c*shots.c~ID+vessel+zone+finyear,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c*shots.c~zone+finyear,data=Km.gn.days,sum,na.rm=T)
    }
    if(use.ID=="NO")
    {
      Km.gn.days=aggregate(netlen.c*shots.c~finyear+vessel+date+zone,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c*shots.c~finyear+zone,data=Km.gn.days,sum,na.rm=T)
    }
    
     Effort.Km.gn.days <- reshape(Km.gn.days,v.names="Km.Gillnet.Days.c",timevar="zone",idvar=c("finyear"),direction="wide")
    iid=match(c("Km.Gillnet.Days.c.West","Km.Gillnet.Days.c.Zone1","Km.Gillnet.Days.c.Zone2"),names(Effort.Km.gn.days))
    Effort.Km.gn.days[,iid]=Effort.Km.gn.days[,iid]/1000
    names(Effort.Km.gn.days)[iid]=c("West","Zone1","Zone2")  
  }
  
  
  #km gn hours
  Km.gn.hours=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+ID+zone,data=SUBSET,max,na.rm=T)
  Km.gn.hours=aggregate(Km.Gillnet.Hours.c~finyear+zone,data=Km.gn.hours,sum,na.rm=T)
  Effort.Km.gn.hours <- reshape(Km.gn.hours,v.names="Km.Gillnet.Hours.c",timevar="zone",idvar=c("finyear"),direction="wide")
  iid=match(c("Km.Gillnet.Hours.c.West","Km.Gillnet.Hours.c.Zone1","Km.Gillnet.Hours.c.Zone2"),names(Effort.Km.gn.hours))
  Effort.Km.gn.hours[,iid]=Effort.Km.gn.hours[,iid]/1000
  names(Effort.Km.gn.hours)[iid]=c("West","Zone1","Zone2")
  
  
  #3.1.2. Calculate longline gillnet effort equivalent
  Daily.shark_rays.zone=subset(Data.daily,SPECIES%in%Elasmo.species & FINYEAR%in%THESE.YRS)
  
  #3.1.2.1 Annual catch by zone
  #GN                                
  Daily.Catch.shark_rays.zone=aggregate(LIVEWT~FINYEAR+zone,data=subset(Daily.shark_rays.zone,METHOD=="GN"),sum,na.rm=T)
  Daily.Catch.shark_rays.zone=reshape(Daily.Catch.shark_rays.zone,v.names="LIVEWT",timevar="zone",
                                      idvar=c("FINYEAR"),direction="wide")
  #LL
  LL.dat=subset(Daily.shark_rays.zone,METHOD=="LL")
  if(nrow(LL.dat)>0)
  {
    Daily.Catch.shark_rays.zone.LL=aggregate(LIVEWT~FINYEAR+zone,data=LL.dat,sum,na.rm=T)
    Daily.Catch.shark_rays.zone.LL=reshape(Daily.Catch.shark_rays.zone.LL,v.names="LIVEWT",timevar="zone",
                                           idvar=c("FINYEAR"),direction="wide")
    Daily.Catch.shark_rays.zone.LL=Daily.Catch.shark_rays.zone.LL[order(Daily.Catch.shark_rays.zone.LL$FINYEAR),]
    Daily.Catch.shark_rays.zone.LL[is.na(Daily.Catch.shark_rays.zone.LL)]=0
    
  }
  
  #3.1.2.2 nominal cpue
  FinYR=Daily.Catch.shark_rays.zone$FINYEAR
  Nom.cpue.days=Daily.Catch.shark_rays.zone[,iid]/Effort.Km.gn.days[,iid]
  Nom.cpue.hours=Daily.Catch.shark_rays.zone[,iid]/Effort.Km.gn.hours[,iid]
  
  #3.1.2.3 equiv effort
    #if there's no LL effort in that year
  LL.equiv.Eff.days.zone=Effort.Km.gn.days
  LL.equiv.Eff.hours.zone=Effort.Km.gn.hours
  LL.equiv.Eff.days.zone[,iid]=0
  LL.equiv.Eff.hours.zone[,iid]=0
  
    #if there's LL effort in that year
  if(nrow(LL.dat)>0)
  {
    #add columns to match dimensions, if no longline catch in some zone
    if(!ncol(Daily.Catch.shark_rays.zone.LL)==ncol(Daily.Catch.shark_rays.zone))
    {
      ii=which(!names(Daily.Catch.shark_rays.zone)%in%names(Daily.Catch.shark_rays.zone.LL))
      ADD=names(Daily.Catch.shark_rays.zone)[ii]
      mat=as.data.frame(matrix(ncol=length(ADD),nrow=nrow(Daily.Catch.shark_rays.zone)))
      names(mat)=ADD
      Daily.Catch.shark_rays.zone.LL=cbind(Daily.Catch.shark_rays.zone.LL,mat)
      Daily.Catch.shark_rays.zone.LL=Daily.Catch.shark_rays.zone.LL[,match(names(Daily.Catch.shark_rays.zone),names(Daily.Catch.shark_rays.zone.LL))]
      Daily.Catch.shark_rays.zone.LL[is.na(Daily.Catch.shark_rays.zone.LL)]=0
    }
    LL.equiv.Eff.days.zone=data.frame(FINYEAR=FinYR,Daily.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.days)
    LL.equiv.Eff.hours.zone=data.frame(FINYEAR=FinYR,Daily.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.hours)
  }
  
  #3.1.3. Total effort
  #by year and zone
  Total.effort.zone.days=data.frame(FINYEAR=FinYR,Effort.Km.gn.days[,iid]+LL.equiv.Eff.days.zone[,iid])
  Total.effort.zone.hours=data.frame(FINYEAR=FinYR,Effort.Km.gn.hours[,iid]+LL.equiv.Eff.hours.zone[,iid])
  
  #by year and joint authority
  Total.effort.joint.days=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.days$West,
                                     JASGL=Total.effort.zone.days$Zone1+Total.effort.zone.days$Zone2)
  Total.effort.joint.hours=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.hours$West,
                                      JASGL=Total.effort.zone.hours$Zone1+Total.effort.zone.hours$Zone2)
  #by year and total
  Total.effort.days=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.days$WCGL+Total.effort.joint.days$JASGL)
  Total.effort.hours=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.hours$WCGL+Total.effort.joint.hours$JASGL)
  
  return(list(Total.effort.zone.days=Total.effort.zone.days,Total.effort.zone.hours=Total.effort.zone.hours,
              Total.effort.joint.days=Total.effort.joint.days,Total.effort.joint.hours=Total.effort.joint.hours,
              Total.effort.days=Total.effort.days,Total.effort.hours=Total.effort.hours))
}

Eff.Current.Yr=fn.Eff.Sofar(Current.yr)
Total.effort.zone.days=Eff.Current.Yr$Total.effort.zone.days
Total.effort.zone.hours=Eff.Current.Yr$Total.effort.zone.hours
Total.effort.joint.days=Eff.Current.Yr$Total.effort.joint.days
Total.effort.joint.hours=Eff.Current.Yr$Total.effort.joint.hours


  #3.2. Monthly records                                                       

    #3.2.1. Aggregate effort by year and zone
#km gn days
Effort.Km.gn.days.monthly=subset(Attach.Effort,zone%in%c("West","Zone1","Zone2"))
Effort.Km.gn.days.monthly=Effort.Km.gn.days.monthly[,match(c("zone","FINYEAR","Km.Gillnet.Days.c"),
                                                           names(Effort.Km.gn.days.monthly))]
Effort.Km.gn.days.monthly <- reshape(Effort.Km.gn.days.monthly,v.names="Km.Gillnet.Days.c",timevar="zone",
                                     idvar=c("FINYEAR"),direction="wide")
iid=match(c("Km.Gillnet.Days.c.West","Km.Gillnet.Days.c.Zone1","Km.Gillnet.Days.c.Zone2"),
          names(Effort.Km.gn.days.monthly))
Effort.Km.gn.days.monthly[,iid]=Effort.Km.gn.days.monthly[,iid]/1000
names(Effort.Km.gn.days.monthly)[iid]=c("West","Zone1","Zone2")

#km gn hours
Effort.Km.gn.hours.monthly=subset(Attach.Effort,zone%in%c("West","Zone1","Zone2"))
Effort.Km.gn.hours.monthly=Effort.Km.gn.hours.monthly[,match(c("zone","FINYEAR","Km.Gillnet.Hours.c"),
                                                           names(Effort.Km.gn.hours.monthly))]
Effort.Km.gn.hours.monthly <- reshape(Effort.Km.gn.hours.monthly,v.names="Km.Gillnet.Hours.c",
                          timevar="zone",idvar=c("FINYEAR"),direction="wide")
Effort.Km.gn.hours.monthly[,iid]=Effort.Km.gn.hours.monthly[,iid]/1000
names(Effort.Km.gn.hours.monthly)[iid]=c("West","Zone1","Zone2")


# select financial years
THESE.YRS=FINYEAR.monthly


    #3.2.2. Calculate longline gillnet effort equivalent                           
Monthly.shark_rays.zone=subset(Data.monthly,SPECIES%in%Elasmo.species & FINYEAR%in%THESE.YRS)

      #3.2.2.1 Annual catch by zone
#GN                                
Monthly.Catch.shark_rays.zone=aggregate(LIVEWT~FINYEAR+zone,data=subset(Monthly.shark_rays.zone,METHOD=="GN"),sum,na.rm=T)
Monthly.Catch.shark_rays.zone=reshape(Monthly.Catch.shark_rays.zone,v.names="LIVEWT",timevar="zone",
                                    idvar=c("FINYEAR"),direction="wide")
#LL
Monthly.Catch.shark_rays.zone.LL=aggregate(LIVEWT~FINYEAR+zone,data=subset(Monthly.shark_rays.zone,METHOD=="LL"),
                                           sum,na.rm=T)
Monthly.Catch.shark_rays.zone.LL=reshape(Monthly.Catch.shark_rays.zone.LL,v.names="LIVEWT",timevar="zone",
                                       idvar=c("FINYEAR"),direction="wide")
Monthly.Catch.shark_rays.zone.LL=Monthly.Catch.shark_rays.zone.LL[order(Monthly.Catch.shark_rays.zone.LL$FINYEAR),]
Monthly.Catch.shark_rays.zone.LL[is.na(Monthly.Catch.shark_rays.zone.LL)]=0

        #3.2.2.2 nominal cpue
FinYR=Monthly.Catch.shark_rays.zone$FINYEAR
II=match(THESE.YRS,Effort.Km.gn.days.monthly$FINYEAR)

Monthly.Catch.shark_rays.zone=Monthly.Catch.shark_rays.zone[order(Monthly.Catch.shark_rays.zone$FINYEAR),]
Effort.Km.gn.days.monthly=Effort.Km.gn.days.monthly[order(Effort.Km.gn.days.monthly$FINYEAR),]
Effort.Km.gn.hours.monthly=Effort.Km.gn.hours.monthly[order(Effort.Km.gn.hours.monthly$FINYEAR),]

Nom.cpue.days=Monthly.Catch.shark_rays.zone[,iid]/Effort.Km.gn.days.monthly[II,iid]
Nom.cpue.hours=Monthly.Catch.shark_rays.zone[,iid]/Effort.Km.gn.hours.monthly[II,iid]

Nom.cpue.days$FINYEAR=Monthly.Catch.shark_rays.zone$FINYEAR
Nom.cpue.hours$FINYEAR=Monthly.Catch.shark_rays.zone$FINYEAR

        #3.2.2.3 equiv effort
#add columns to match dimensions, if no longline catch in some zone
if(!ncol(Monthly.Catch.shark_rays.zone.LL)==ncol(Monthly.Catch.shark_rays.zone))
{
  ii=which(!names(Monthly.Catch.shark_rays.zone)%in%names(Monthly.Catch.shark_rays.zone.LL))
  ADD=names(Monthly.Catch.shark_rays.zone)[ii]
  mat=as.data.frame(matrix(ncol=length(ADD),nrow=nrow(Monthly.Catch.shark_rays.zone)))
  names(mat)=ADD
  Monthly.Catch.shark_rays.zone.LL=cbind(Monthly.Catch.shark_rays.zone.LL,mat)
  Monthly.Catch.shark_rays.zone.LL=Monthly.Catch.shark_rays.zone.LL[,match(names(Monthly.Catch.shark_rays.zone),
                                  names(Monthly.Catch.shark_rays.zone.LL))]
  Monthly.Catch.shark_rays.zone.LL[is.na(Monthly.Catch.shark_rays.zone.LL)]=0
}

LL.equiv.Eff.days.zone=merge(Monthly.Catch.shark_rays.zone.LL,Nom.cpue.days,by="FINYEAR",all=T)
LL.equiv.Eff.hours.zone=merge(Monthly.Catch.shark_rays.zone.LL,Nom.cpue.hours,by="FINYEAR",all=T)

LL.equiv.Eff.days.zone$LIVEWT.West=with(LL.equiv.Eff.days.zone,LIVEWT.West.x/LIVEWT.West.y)
LL.equiv.Eff.days.zone$LIVEWT.Zone1=with(LL.equiv.Eff.days.zone,LIVEWT.Zone1.x/LIVEWT.Zone1.y)
LL.equiv.Eff.days.zone$LIVEWT.Zone2=with(LL.equiv.Eff.days.zone,LIVEWT.Zone2.x/LIVEWT.Zone2.y)

LL.equiv.Eff.hours.zone$LIVEWT.West=with(LL.equiv.Eff.hours.zone,LIVEWT.West.x/LIVEWT.West.y)
LL.equiv.Eff.hours.zone$LIVEWT.Zone1=with(LL.equiv.Eff.hours.zone,LIVEWT.Zone1.x/LIVEWT.Zone1.y)
LL.equiv.Eff.hours.zone$LIVEWT.Zone2=with(LL.equiv.Eff.hours.zone,LIVEWT.Zone2.x/LIVEWT.Zone2.y)

LL.equiv.Eff.days.zone=subset(LL.equiv.Eff.days.zone,select=c(FINYEAR,LIVEWT.West,LIVEWT.Zone1,LIVEWT.Zone2))
LL.equiv.Eff.hours.zone=subset(LL.equiv.Eff.hours.zone,select=c(FINYEAR,LIVEWT.West,LIVEWT.Zone1,LIVEWT.Zone2))

LL.equiv.Eff.days.zone[is.na(LL.equiv.Eff.days.zone)]=0
LL.equiv.Eff.hours.zone[is.na(LL.equiv.Eff.hours.zone)]=0


 #3.2.3. Total effort (gillnet  plus longline equivalent. This is what's reported in Sofar)

#by year and zone
Total.effort.zone.days.monthly=data.frame(FINYEAR=FinYR,Effort.Km.gn.days.monthly[II,iid]+LL.equiv.Eff.days.zone[,iid])
Total.effort.zone.hours.monthly=data.frame(FINYEAR=FinYR,Effort.Km.gn.hours.monthly[II,iid]+LL.equiv.Eff.hours.zone[,iid])

#by year and joint authority
Total.effort.joint.days.monthly=data.frame(FINYEAR=FinYR,
            WCGL=Total.effort.zone.days.monthly$West,
            JASGL=Total.effort.zone.days.monthly$Zone1+Total.effort.zone.days.monthly$Zone2)
Total.effort.joint.hours.monthly=data.frame(FINYEAR=FinYR,
            WCGL=Total.effort.zone.hours.monthly$West,
            JASGL=Total.effort.zone.hours.monthly$Zone1+Total.effort.zone.hours.monthly$Zone2)
#by year and total
  #TDGDLF
Total.effort.days.monthly=data.frame(FINYEAR=FinYR,
              Total=Total.effort.joint.days.monthly$WCGL+
                          Total.effort.joint.days.monthly$JASGL)
Total.effort.hours.monthly=data.frame(FINYEAR=FinYR,
              Total=Total.effort.joint.hours.monthly$WCGL+
                          Total.effort.joint.hours.monthly$JASGL)

  #NSF  
#missing:add GN equivalent LL (ask Rory for rule)
Effort.monthly.NSF=subset(Effort.monthly,zone%in%c("Closed","Joint","North"),
                     select=c(FINYEAR,MONTH,BLOCKX,zone,VESSEL,METHOD,Same.return,
                                   FDAYS,BDAYS.c,HOURS.c,HOOKS,SHOTS.c,NETLEN.c,nlines,
                                   Km.Gillnet.Days.c,Km.Gillnet.Hours.c))

Effort.daily.NSF=subset(Effort.daily,zone%in%c("Closed","Joint","North"),
                    Select=c(date,finyear,month,blockx,zone,vessel,method,Same.return,Same.return.SNo,
                                 fdays,bdays.c,hours.c,hooks,shots.c,netlen.c,nlines,
                                 Km.Gillnet.Days.c,Km.Gillnet.Hours.c))

Effort.monthly.NSF$hook.days=with(Effort.monthly.NSF,BDAYS.c*HOOKS)
Effort.monthly.NSF$hook.hours=with(Effort.monthly.NSF,BDAYS.c*HOOKS*HOURS.c)

Effort.daily.NSF$hook.days=with(Effort.daily.NSF,bdays.c*hooks)
Effort.daily.NSF$hook.hours=with(Effort.daily.NSF,bdays.c*hooks*hours.c)


Attach.Effort.monthly.c_NSF=aggregate(cbind(hook.days,hook.hours)~MONTH+BLOCKX+VESSEL+zone+FINYEAR,
                                      data=Effort.monthly.NSF,max,na.rm=T)
Attach.Effort.monthly.c_NSF=aggregate(cbind(hook.days,hook.hours)~FINYEAR,
                                      data=Attach.Effort.monthly.c_NSF,sum,na.rm=T)
Attach.Effort.monthly.c_NSF=subset(Attach.Effort.monthly.c_NSF)

Attach.Effort.daily.c_NSF=aggregate(cbind(hook.days,hook.hours)~date+vessel+zone+finyear,
                                    data=Effort.daily.NSF,max,na.rm=T)
Attach.Effort.daily.c_NSF=aggregate(cbind(hook.days,hook.hours)~finyear,
                                    data=Attach.Effort.daily.c_NSF,sum,na.rm=T)
names(Attach.Effort.monthly.c_NSF)=names(Attach.Effort.daily.c_NSF)

#Total.effort_NFS=rbind(Attach.Effort.monthly.c_NSF,Attach.Effort.daily.c_NSF)  #issue: check daily 2007-08, too high effrot
Total.effort_NFS=Attach.Effort.monthly.c_NSF
Total.effort_NFS$hook.days=Total.effort_NFS$hook.days/1000
Total.effort_NFS$hook.hours=Total.effort_NFS$hook.hours/1000
Total.effort_NFS=aggregate(cbind(hook.days,hook.hours)~finyear,Total.effort_NFS,sum)

#add dummy effort from , sort out issue of 2007-08 daily beeing too high!!
Total.effort_NFS[32:34,2]=300
names(Total.effort_NFS)=c("FINYEAR","Hook days","Hook hours")


#Compare Rory's effort and current script's effort
if(Inspect.New.dat=="YES")
{
  Results.pre.2013=read.csv("C:/Matias/Data/Catch and Effort/Historic/Historic.res.csv")
  B=Total.effort.zone.days.monthly
  A=Total.effort.days.monthly
  COLs=2:4
  Yr=unique(A$FINYEAR)
  nYr=1:length(Yr)
  
  tiff(file=paste(Hndl,"Compare.effort_Rory_new.tiff",sep=""),width = 2400, height = 2400,units = "px",
       res = 300, compression = "lzw")
  par(mfcol=c(2,1))
  
  #Total
  Yr=unique(A$FINYEAR)
  nYr=1:length(Yr)
  plot(nYr,ylim=c(0,max(c(A$Total,Results.pre.2013$TDGDLF.km.gn.days/1000))),col="white",ylab="1000 Km.Gillnet.Days.c",xlab="Year",xaxt="n",main="KM.gn.days")
  axis(1,nYr,F,tck=-0.015)
  axis(1,seq(nYr[5],nYr[length(nYr)],5),Yr[seq(nYr[5],nYr[length(nYr)],5)],tck=-0.035,cex.axis=0.8)
  lines(nYr,A$Total,col=1,lwd=2)
  lines(1:36,Results.pre.2013$TDGDLF.km.gn.days/1000,lwd=2,lty=2)
  legend("topleft",c("Matias","Rory"),bty='n',lty=1:2,col=1)
  
   
  #By zone
  COLs=2:4
  plot(nYr,ylim=c(0,30),col="white",ylab="1000 Km.Gillnet.Days.c",xlab="Year",xaxt="n",main="KM.gn.days")
  axis(1,nYr,F,tck=-0.015)
  axis(1,seq(nYr[5],nYr[length(nYr)],5),Yr[seq(nYr[5],nYr[length(nYr)],5)],tck=-0.035,cex.axis=0.8)
  lines(nYr,B$West[nYr],col=COLs[1])
  lines(nYr,B$Zone1[nYr],col=COLs[2])
  lines(nYr,B$Zone2[nYr],col=COLs[3])
  
  Yr=unique(Results.pre.2013$FINYEAR)
  nYr=1:length(Yr)
  lines(nYr,Results.pre.2013$WC.km.gn.days/1000,col=COLs[1],lty=2)
  lines(nYr,Results.pre.2013$Z1.km.gn.days/1000,col=COLs[2],lty=2)
  lines(nYr,Results.pre.2013$Z2.km.gn.days/1000,col=COLs[3],lty=2)
  
  legend("topright",c("West","Zn1","Zn2"),bty='n',lty=1,col=c(COLs))
  legend("topleft",c("Matias","Rory"),bty='n',lty=1:2,col=1)
  dev.off()
}


#4. Prepare files for cpue standardisation         
 
  #4.1. Some final changes

    #4.1.1 add the variable netlen.c to Data.monthly
    #monthly
Ef.netlen.daily=aggregate(netlen.c~Same.return,Effort.daily,max,na.rm=T)
Ef.netlen.monthly=aggregate(NETLEN.c~Same.return,Effort.monthly,max,na.rm=T)
names(Ef.netlen.daily)=names(Ef.netlen.monthly)
Ef.netlen=rbind(Ef.netlen.monthly,Ef.netlen.daily)
Data.monthly=Data.monthly%>%left_join(Ef.netlen,by=c("Same.return"))
# Data.monthly=merge(Data.monthly,Ef.netlen,by="Same.return",all.x=T)

    #daily
Ef.netlen.daily.1=aggregate(netlen.c~Same.return.SNo,Effort.daily,max,na.rm=T)
Data.daily=Data.daily%>%left_join(Ef.netlen.daily.1,by="Same.return.SNo")
#Data.daily=merge(Data.daily,Ef.netlen.daily.1,by="Same.return.SNo",all.x=T)


    #4.1.2 modify over-reported monthly catches
#Too.high.catch=subset(Data.monthly,LIVEWT.c>Top.mon.ktch)
#Data.monthly$LIVEWT.c=with(Data.monthly,ifelse(YEAR.c < 2006 & LIVEWT.c>Top.mon.ktch,Top.mon.ktch,LIVEWT.c))

  #4.2 Separate gillnet from other methods     
#note: exclude estuaries       
    #monthly
Data.monthly.GN=subset(Data.monthly,METHOD=="GN" & Estuary=="NO")
Data.monthly.LL=subset(Data.monthly,METHOD=="LL" & Estuary=="NO")
Data.monthly.other=subset(Data.monthly,!(METHOD=="LL"|METHOD=="GN")
                    |((METHOD=="LL"|METHOD=="GN") & Estuary=="YES"))

    #daily
Data.daily.GN=subset(Data.daily,METHOD=="GN" & Estuary=="NO")
Data.daily.LL=subset(Data.daily,METHOD=="LL" & Estuary=="NO")
Data.daily.other=subset(Data.daily,!(METHOD=="LL"|METHOD=="GN")
                          |((METHOD=="LL"|METHOD=="GN") & Estuary=="YES"))

#Check use of different gears by year and zone for main species
fn.check.method.contrib=function(D,SPEC)
{
  hndl.gear="C:/Matias/Analyses/Catch and effort/Outputs/"
  a=subset(D,  Estuary=="NO" & LIVEWT.c>0 & LAT<=(-26) & SPECIES==SPEC)
  a$METHOD1=with(a,ifelse(!METHOD%in%c("GN","HL","LL"),"Other",METHOD))
  Rcrds.yr.mthd=with(a,table(FINYEAR,METHOD1))
  Rcrds.yr.mthd=Rcrds.yr.mthd/rowSums(Rcrds.yr.mthd)
  Rcrds.yr.mthd=t(Rcrds.yr.mthd)
  
  
  fn.fig(paste(hndl.gear,"TDGDLF_gear_by_year.",a$SNAME[1],sep=""),2400, 2400)
  par(mfcol=c(1,1),mai=c(1.1,.85,.2,.1),oma=c(.1,.1,1,.1),xpd=TRUE,las=1)
  
  #Zones combined
  barplot(Rcrds.yr.mthd,legend.text=rownames(Rcrds.yr.mthd),
          args.legend=list(c(0,1),cex=1.15,bty='n',horiz=T,
                           border='black',yjust=0))
  mtext("Financial year",1,-2,outer=T,cex=1.5)
  mtext("Proportion of records",2,-2,outer=T,las=3,cex=1.5)
  
  dev.off()
  
  
  Rcrds.GN.LL=with(a,table(FINYEAR,METHOD1,zone))
  
  fn.fig(paste(hndl.gear,"TDGDLF_gear_by_year_zone.",a$SNAME[1],sep=""),2400, 2400)
  par(mfcol=c(3,1),mai=c(.3,.8,.17,.01),oma=c(4,2,1.5,.1),xpd=TRUE,las=1,mgp=c(.65,1,0))
  
  #West
  x=Rcrds.GN.LL[,,1]/rowSums(Rcrds.GN.LL[,,1])
  x=t(x)
  barplot(x,legend.text=rownames(x),
          args.legend=list(x=ncol(x)*1.2,
                           y=max(x)*.975,cex=1.15,bty='n',horiz=T,
                           border='black',yjust=0))
  mtext(dimnames(Rcrds.GN.LL)$zone[1],3)
  
  #Zn1
  x=Rcrds.GN.LL[,,2]/rowSums(Rcrds.GN.LL[,,2])
  x=t(x)
  barplot(x)
  mtext(dimnames(Rcrds.GN.LL)$zone[2],3)
  
  #Zn2
  x=Rcrds.GN.LL[,,3]/rowSums(Rcrds.GN.LL[,,3])
  x=t(x)
  barplot(x)
  mtext(dimnames(Rcrds.GN.LL)$zone[3],3)
  
  mtext("Financial year",1,0.5,outer=T,cex=1.5)
  mtext("Proportion of records",2,-3.5,outer=T,las=3,cex=1.5)
  
  dev.off()
  
}
fn.check.method.contrib(D=Data.monthly,SPEC=17001)
fn.check.method.contrib(D=Data.monthly,SPEC=17003)
fn.check.method.contrib(D=Data.monthly,SPEC=18003)
fn.check.method.contrib(D=Data.monthly,SPEC=18007)


#Data for Hammerhead Listing
if(do.Hammerheads=="YES")
{
  HHEads=c(19000)
  Sel.HH=c("SPECIES","FINYEAR","MONTH","METHOD","VESSEL","BLOCKX","LAT","LONG","zone","LIVEWT.c")
  Hammerheads=subset(Data.monthly,SPECIES%in% HHEads,select=Sel.HH)
}
 

#Data for Jeff Norris
if (do.Jeffs=="YES") Jeff.Norris=subset(Data.monthly,SPECIES%in%Scalefish.species & Bioregion=="SC")  #Data requested by Jeff Norris

#ACA
#add nfish
Data.monthly.GN$BLOCKX=as.numeric(Data.monthly.GN$BLOCKX)
Data.monthly.GN=Data.monthly.GN%>%left_join(Daily.nfish.agg,by=c("Same.return",
           "FINYEAR","YEAR","YEAR.c","MONTH","zone","BLOCKX","VESSEL","METHOD","SPECIES","Estuary"))
# Data.monthly.GN=merge(Data.monthly.GN,Daily.nfish.agg,by=c("Same.return",
#   "FINYEAR","YEAR","YEAR.c","MONTH","zone","BLOCKX","VESSEL","METHOD","SPECIES","Estuary"),all.x=T)

if (do.Jeffs=="YES") Jeff.Norris=merge(Jeff.Norris,Daily.nfish.agg,by=c("Same.return",
        "FINYEAR","YEAR","YEAR.c","MONTH","zone","BLOCKX","VESSEL","METHOD","SPECIES","Estuary"),all.x=T)


  #4.3. check spatial expansion
BLKS.YEAR=Expand.fun(Data.monthly.GN)
Spatial.expan=round(100*BLKS.YEAR/N.blocks,0)
Spatial.expan.north=round(100*Expand.fun(Data.monthly.north)/N.blocks,0)
Effort.expan=Effort1.fun(Data.monthly.GN)
Effort.expan.north=Effort1.fun(Data.monthly.north)


  #4.4. Define effective area effort                                   #Rory's rule 5a-5k
Data.monthly.GN$Boundary.blk=with(Data.monthly.GN,
          ifelse(BLOCKX%in%Boundary.blk & !FINYEAR%in%Daily.l.years,"Y","N"))



#SECTION F 2. ---- EXPORT DATA FOR ASSESSMENT AND CPUE STANDARDISATION ----

setwd("C:/Matias/Analyses/Catch and effort/Data_outs")

#some final amendments
crap=c("GoodsplitID","Prop.sandbar","SanBar.rep",            
       "Prop.Gum.Ves.ID","Prop.Whi.Ves.ID","Prop.Dus.Ves.ID","Prop.Other.Ves.ID",     
       "Prop.Sch.Ves.ID","Prop.DogS.Ves.ID","Prop.San.Ves.ID","Prop.Dus.Good.spl",     
       "Prop.Gum.Good.spl","Prop.Whi.Good.spl","Prop.Dus.Zone.Good.spl",
       "Prop.Gum.Zone.Good.spl","Prop.Whi.Zone.Good.spl","Prop.Dus.Mon.Good.spl",
       "Prop.Gum.Mon.Good.spl","Prop.Whi.Mon.Good.spl","Prop.Dus.YrZn.Good.spl",
       "Prop.Gum.YrZn.Good.spl","Prop.Whi.YrZn.Good.spl")
crap.daily=c("YEAR","GoodsplitID","GoodsplitID","BlockAveID",
             "AnnualVesselAveID","MonthlyID","BlockID","ZoneID",
             "Prop.Sch.Ves.ID","Prop.DogS.Ves.ID")
crap.ef=c("AnnualVesselAveID_BDAYS.m","AnnualVesselAveID_HOURS.m","AnnualVesselAveID_SHOTS.m", 
       "AnnualVesselAveID_NETLEN.m","MonthlyZoneID_BDAYS.m","MonthlyZoneID_HOURS.m",     
       "MonthlyZoneID_SHOTS.m","MonthlyZoneID_NETLEN.m")


#Export all these objects
#note: Total.effort... : annual effort reported in SOFAR (GN plus LL equivalent)
#      Data.monthly has catch from all gears

Exprt.list=list(
  Annual.total.eff.days=Total.effort.days.monthly,
  Annual.total.eff.hours=Total.effort.hours.monthly,
  Annual.zone.eff.hours=Total.effort.zone.hours.monthly,
  Annual.zone.eff.days=Total.effort.zone.days.monthly,
  Annual.total.eff_NSF=Total.effort_NFS,
  Effort.monthly=Effort.monthly[,-match(crap.ef,names(Effort.monthly))],
  Effort.daily=Effort.daily[,-match(crap.ef,names(Effort.daily))],
  Mesh.monthly=Mesh.monthly,
  Mesh.size=Mesh.size,
  TEPS.current=TEPS.current,
  Rec.fish.catch=Rec.fish.catch,
  Data.current.Sofar=Data.current.Sofar,
  PRICES=PRICES,
  Data.monthly=Data.monthly[,-match(crap,names(Data.monthly))],
  Data.monthly.north=Data.monthly.north[,-match(crap,names(Data.monthly.north))],
  Data.monthly.GN=Data.monthly.GN[,-match(crap,names(Data.monthly.GN))],
  Data.monthly.LL=Data.monthly.LL[,-match(crap,names(Data.monthly.LL))],
  Data.monthly.other=Data.monthly.other,
  Data.daily=Data.daily[,-match(crap.daily,names(Data.daily))],
  Data.daily.original=Data.daily.original,
  Data.daily.GN=Data.daily.GN[,-match(crap.daily,names(Data.daily.GN))],
  Data.daily.LL=Data.daily.LL[,-match(crap.daily,names(Data.daily.LL))],
  Data.daily.other=Data.daily.other[,-match(crap.daily,names(Data.daily.other))])
for(i in 1:length(Exprt.list)) fwrite(Exprt.list[[i]],paste(names(Exprt.list)[i],".csv",sep=""),row.names=F)
rm(Exprt.list)



  #Total catch SAFS for EVA
if(Export.SAFS=="YES")
{
  vars.Eva.SAFS=c("FINYEAR","YEAR.c","MONTH","METHOD","SNAME","SPECIES","VESSEL","zone","LIVEWT.c",
                  "FisheryCode","FisheryZone","BLOCKX")
  Eva.SAFS=Data.monthly[,match(vars.Eva.SAFS,names(Data.monthly))]
  Eva.SAFS.north=Data.monthly.north[,match(vars.Eva.SAFS,names(Data.monthly.north))]
  Eva.SAFS=rbind(Eva.SAFS,Eva.SAFS.north)
  write.csv(Eva.SAFS,"M:/Fisheries Research/Common Data Area/Matias/For Eva/SAFS.csv",row.names=F)
  rm(Eva.SAFS.north,Eva.SAFS)  
}


##################--- G. REPORT SECTION ---##############
#------ G 1. SPATIO TEMPORAL CATCH AND EFFORT DISTRIBUTION ------
if(First.run=="YES")
{
  # Plot spatio-temporal catches of 4 main shark species
  Effect.area=function(SPEC,RANGO)
  {
    a=subset(Data.monthly,SPECIES%in%SPEC)
    yr=sort(unique(a$FINYEAR))
    n=length(yr)
    
    agg=aggregate(LIVEWT.c~FINYEAR+BLOCKX+LAT+LONG,a,sum)
    
    for(i in 1:n)
    {
      b=subset(agg,FINYEAR==yr[i])
      z=b$LIVEWT.c/max(b$LIVEWT.c)
      plot(b$LONG,b$LAT,cex=z*2,main=yr[i],ylab="",xlab="",pch=19,ylim=c(-36,-26),
           xlim=c(113,129),cex.axis=0.8,cex.main=.85,col="steelblue4")
      legend("top",paste(round(max(b$LIVEWT.c)/1000),"tons"),pch=19,pt.cex=2,bty='n',col="steelblue4")
      Y1=RANGO[RANGO<0]
      Y2=-36
      if(length(Y1)>0)
      {
        X2=RANGO[RANGO>0]
        X1=112.5
      }
      if(length(Y1)==0)
      {
        Y1=-31
        X1=RANGO[1]
        X2=RANGO[2]
      }
      polygon(c(X1,X2,X2,X1),c(Y2,Y2,Y1,Y1),col=rgb(0.1, .1, .1, 0.1),border="transparent")
    }
  }
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Paper/Catch_spatial_temporal")
  tiff("Spatio.temporal.catch.Whiskery.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  n.graph=length(unique(subset(Data.monthly,SPECIES==Whiskery)$FINYEAR))
  smart.par(n.plots=n.graph,MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  Effect.area(Whiskery,Whiskery.range)
  plot(1:10,col="transparent",xaxt='n',ann=F,yaxt='n',fg="white")
  text(5,7,"Whiskery",cex=1.5)
  text(5,4,"shark",cex=1.75)
  dev.off()
  
  
  tiff("Spatio.temporal.catch.gummy.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  n.graph=length(unique(subset(Data.monthly,SPECIES==Gummy)$FINYEAR))
  smart.par(n.plots=n.graph,MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  Effect.area(Gummy,Gummy.range)
  plot(1:10,col="transparent",xaxt='n',ann=F,yaxt='n',fg="white")
  text(5,7,"Gummy",cex=1.5)
  text(5,4,"shark",cex=1.75)
  dev.off()
  
  
  tiff("Spatio.temporal.catch.dusky.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  n.graph=length(unique(subset(Data.monthly,SPECIES%in%Dusky_whaler)$FINYEAR))
  smart.par(n.plots=n.graph,MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  Effect.area(Dusky_whaler,Dusky.range)
  plot(1:10,col="transparent",xaxt='n',ann=F,yaxt='n',fg="white")
  text(5,7,"Dusky",cex=1.5)
  text(5,4,"shark",cex=1.75)
  dev.off()
  
  
  tiff("Spatio.temporal.catch.sandbar.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  n.graph=length(unique(subset(Data.monthly,SPECIES==Sandbar)$FINYEAR))
  smart.par(n.plots=n.graph,MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  Effect.area(Sandbar,Sandbar.range)
  plot(1:10,col="transparent",xaxt='n',ann=F,yaxt='n',fg="white")
  text(5,7,"Sandbar",cex=1.5)
  text(5,4,"shark",cex=1.75)
  dev.off()
  
  #add effort to Data.monthly.GN
  s=subset(Eff,LAT<=(-26))
  s$METHOD="GN"
  s$Same.return=with(s,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
  s=s[,-match(c("BLOCKX","FINYEAR","MONTH","VESSEL","METHOD","LAT","LONG","Eff.Reporter"),names(s))]
  Data.monthly.GN=Data.monthly.GN%>%left_join(s,by="Same.return")
  #Data.monthly.GN=merge(Data.monthly.GN,s,by="Same.return",all.x=T)
  
  #grouping years
  Last.calendar.Yr=as.numeric(substr(Current.yr,1,4))
  yrs.grupd=4
  Yr.group=seq(1975,(Last.calendar.Yr+1),5)
  Yr.group.plus=Yr.group+yrs.grupd
  Yr.group.plus[length(Yr.group.plus)]=min(Yr.group.plus[length(Yr.group.plus)],Last.calendar.Yr)
  Yr.range=vector('list',length(Yr.group))
  for(e in 1:length(Yr.group))Yr.range[[e]]=paste(Yr.group[e],"-",Yr.group.plus[e])
  Yr.range=unlist(Yr.range)  
  
  #create data list
  DATA.lista=vector('list',length(Yr.range))
  names(DATA.lista)=Yr.range
  for(i in 1:(length(DATA.lista))) DATA.lista[[i]]=subset(Data.monthly.GN,YEAR.c>=Yr.group[i] & YEAR.c<=Yr.group.plus[i])
  
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2))
  
  add.depth="NO"
  if(add.depth=="YES")reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
  
  data(worldLLhigh)
  
  #define coordinates of plots
  South.WA.lat=c(-36,-25); South.WA.long=c(112,129)
  S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
  S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
  
  OZ.lat=c(-44.5,-11);OZ.long=c(113,155)
  
  Perth=c(115.866,-31.95)
  Rotnest=c(115.50,-32.02)
  
  Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
  Lat.seq=c(-26,-28,-30,-32,-34)
  
  numInt=20  #number of intervals for effort 
  
  #couleurs=rev(gray(0:(numInt-1)/(numInt-1)))
  couleurs=rev(gray(seq(0.2,0.9,length=numInt)))
  #couleurs  <- tail(topo.colors(trunc(1.4 * numInt)),numInt)
  numberLab=5
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  
  a=112:129
  b=seq(-37,South.WA.lat[2],length.out=length(a))
  
  
  #max block effort for each year period
  Ef.var="days"   #express effort in km gn days
  #Ef.var="hours"
  
  MAX.EFF=data.frame(Period=Yr.range,Max.eff=NA)
  EFF.bin=NULL
  for (i in 1:length(DATA.lista))
  {
    DATA=DATA.lista[[i]]
    DATA=DATA[-which(duplicated(DATA$Same.return)),]
    if(Ef.var=="days") Max.effort=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX),FUN=sum,na.rm=T))
    if(Ef.var=="hours") Max.effort=with(DATA,aggregate(Km.Gillnet.Hours.c,list(BLOCKX),FUN=sum,na.rm=T))
    EFF.bin=c(EFF.bin,Max.effort[,2])
    MAX.EFF[i,2]=max(Max.effort[,2])
  }
  Max.effort=ceiling(max(MAX.EFF[,2])) 
  # dummy=cbind(128.5,-35.5,Max.effort)
  fn.eff.plot=function(DATA,tcl.1,tcl.2,EffortBreaks)   #function plot effort
  {
    DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
    DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
    DATA$BLOCKX.c=with(DATA,paste(-LAT,LONG-100,sep=""))
    MapEffort=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX.c),FUN=sum,na.rm=T))    
    
    colnames(MapEffort)=c("BLOCKX.c","Total effort")
    id=unique(match(MapEffort$BLOCKX.c,DATA$BLOCKX.c))
    MapEffort$LAT=DATA$LAT[id]
    MapEffort$LONG=DATA$LONG[id]
    lat=sort(unique(MapEffort$LAT))
    ids=match(-27,lat)
    if(is.na(ids))      #add dummy for sorting image out of whack when missing Lat
    {
      adD=MapEffort[1,]
      adD[,2]=NA
      adD$BLOCKX.c="-27 113"
      adD$LAT=-27
      MapEffort=rbind(MapEffort,adD)
    }
    
    MapEffort$LAT.cen=MapEffort$LAT-.5
    MapEffort$LONG.cen=MapEffort$LONG+.5  
    MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
    MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
    long=sort(unique(MapEffort$LONG.cen))
    lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image  
    
    MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Total effort"),names(MapEffort))]  
    #MapEffort=rbind(MapEffort,dummy)#keep in perspective
    
    Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",    #transposed as matrix   
                               timevar="LAT.cen",v.names="Total effort", direction="wide"))  
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]  									
    
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=EffortBreaks,axes = FALSE,add=T)			
    axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
    axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
    
    par(new=T)
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                 nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
  }
  
  EffortBreakSS=quantile(EFF.bin,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
  
  PLATE=c(.01,.9,.075,.9)
  
  LATT=South.WA.lat[2]:South.WA.lat[1]
  LONGG=South.WA.long[1]:South.WA.long[2]
  
  South.WA.lat=c(-37,-25)
  Lat.seq=c(-26,-28,-30,-32,-34,-36)
  
  #Create figures 1 to 5
  if (plot.cpue.paper.figures=="YES")
  {
    tiff(file="Figure 1. Map.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mar=c(2,2,2,2),oma=c(1,1,1,1))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
    polygon(x=c(116.5,116.5,112,112),y=c(-26.5,-33,-33,-26.5),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
    polygon(x=c(116.5,116.5,112,112),y=c(-33,-37,-37,-33),lwd=1.5,col=rgb(.3,.3,.3,alpha=.5))
    polygon(x=c(129,129,116.5,116.5),y=c(-30,-37,-37,-30),lwd=1.5,col=rgb(.7,.7,.7,alpha=.2))
    axis(side = 1, at =LONGG, labels = F, tcl = 34,lty=3,col="grey30")
    axis(side = 4, at = LATT, labels = F,tcl =34,lty=3,col="grey30")
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
    mtext("Latitude (?S)",side=2,line=1.7,las=3,cex=1.75)
    mtext("Longitude (?E)",side=1,line=1.75,cex=1.75)
    
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
  couleurs=rev(heat.colors(numInt))
  
  
  #Spatial Catch by year period
  if(First.run=='YES')
  {
    #catch breaks
    sp=c("whiskery","gummy","dusky","sandbar")
    MAX.CATCH=data.frame(Period=Yr.range,Max.catch=NA)
    CATCH.bin=NULL
    fn.catch.breaks=function(SP)
    {
      
      for (j in 1:length(DATA.lista))
      {
        DATA=subset(DATA.lista[[j]],LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9) & SPECIES%in%SP)
        if(nrow(DATA)>0)
        {
          Max.catch=with(DATA,aggregate(LIVEWT.c,list(BLOCKX),FUN=sum,na.rm=T))
          CATCH.bin=c(CATCH.bin,Max.catch[,2])
          MAX.CATCH[j,2]=max(Max.catch[,2])
        }
      }
      
      Breaks=quantile(CATCH.bin,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
      return(Breaks)
    }
    
    #Plot catch by 5-year period   
    fn.catch.plot=function(DATA,SP,tcl.1,tcl.2,BREAKS)
    {
      DATA=subset(DATA,LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9) & SPECIES%in%SP)
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      
      if(nrow(DATA)<=2)   plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      if(nrow(DATA)>2)
      {
        DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
        
        MapCatch=with(DATA,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))    
        colnames(MapCatch)=c("BLOCKX.c","Total Catch")
        id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
        MapCatch$LAT=DATA$LAT[id]
        MapCatch$LONG=DATA$LONG[id]
        
        MapCatch$LAT.cen=MapCatch$LAT-.5
        MapCatch$LONG.cen=MapCatch$LONG+.5  
        MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
        MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
        long=sort(unique(MapCatch$LONG.cen))
        lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
        
        MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
        
        
        Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                   timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
        Reshaped=Reshaped[order(Reshaped[,1]),]
        Reshaped=Reshaped[,-1]										
        
        
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
        axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
        axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
        
        par(new=T)
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        
      }
      
      if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                   nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
    }
    hndl.sptl.ktch="C:/Matias/Analyses/Catch and effort/Outputs/Paper/Catch_spatial_temporal/"
    Tar=TARGETS
    #Tar[[3]]=18003
    
    # for (ss in 1:length(Tar))   
    # {
    #   fn.fig(paste(hndl.sptl.ktch,"Spatio.temporal.catch.",sp[ss],".Catch_map",sep=""),2400, 2400) 
    #   
    #   opar <- smart.par(n.plots=length(DATA.lista),MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(.1, 0.15, 0))
    #   for (i in 1:length(DATA.lista))
    #   {
    #     Breaks=fn.catch.breaks(Tar[ss])
    #     fn.catch.plot(DATA.lista[[i]],Tar[ss],tcl.1=.5,tcl.2=.5,Breaks)
    #     #fn.catch.plot(DATA.lista[[i]],Tar[ss],tcl.1=16,tcl.2=18.15,Breaks)
    #     mtext(Yr.range[i],side=3,line=-2,cex=.95)
    #     axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .5,las=1,cex.axis=1,padj=-.15)
    #     axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .5,las=2,cex.axis=1,hadj=1.1)
    #     if(i==length(DATA.lista))color.legend(126,-26,129,-30.5,paste(round(Breaks/1000,0),'t'),
    #                                           rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
    #     # if(i==1)
    #     # {
    #     #   par(fig=c(0.2,.35,.8,0.95), new = T,mgp=c(.25,.2,0),las=1)
    #     #   plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
    #     #           col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    #     #   box()
    #     #   polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
    #     #   par(opar) 
    #     # }
    #   }
    #   mtext("Latitude (?S)",side=2,line=0.4,las=3,cex=1.3,outer=T)
    #   mtext("Longitude (?E)",side=1,line=0.6,cex=1.3,outer=T)
    #   dev.off()
    # }
    
    #plot TDGDLF catch by year and grouped years
    Colfunc <- colorRampPalette(c("yellow","red"))
    Couleurs=c("white",Colfunc(numInt-1))
    #Couleurs=Colfunc(numInt)
    fn.ctch.plot.all.yrs=function(DATA,tcl.1,tcl.2,numInt) 
    {
      DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
      DATA$LONG=as.numeric(DATA$LONG)
      DATA$blk=substr(DATA$BLOCKX,1,4)
      A=aggregate(LIVEWT.c~FINYEAR+blk,DATA,sum)
      Ymax=max(A$LIVEWT.c)
      Ymin=min(A$LIVEWT.c)
      #Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
      Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      for(y in 1:length(FINYrS))
      {
        A=subset(DATA,FINYEAR==FINYrS[y])
        MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
        colnames(MapCatch)=c("BLOCKX.c","Total Catch")
        id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
        MapCatch$LAT=DATA$LAT[id]
        MapCatch$LONG=DATA$LONG[id]
        msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
        msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
        if(length(msn.lat)>0)
        {
          dummy=MapCatch[1:length(msn.lat),]
          dummy$`Total Catch`=0
          dummy$LAT=msn.lat
          dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
          MapCatch=rbind(MapCatch,dummy)
        }
        
        if(unique(min(A$YEAR.c))>2007)
        {
          MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
        }
        MapCatch$LAT.cen=MapCatch$LAT-.5
        MapCatch$LONG.cen=MapCatch$LONG+.5  
        MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
        MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
        long=sort(unique(MapCatch$LONG.cen))
        lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
        MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
        Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                   timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
        Reshaped=Reshaped[order(Reshaped[,1]),]
        Reshaped=Reshaped[,-1]	
        numberLab=10
        colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
        
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        image(long,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,axes = FALSE,add=T)			
        axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
        axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
        par(new=T)
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        legend('top',FINYrS[y],bty='n',cex=1.2)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      }
      plot(a,b,ann=F,axes=F,col='transparent')
      color.legend(quantile(a,probs=.25),quantile(b,probs=.95),quantile(a,probs=.6),quantile(b,probs=.25),
                   paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.75)
    }
    fn.ctch.plot.grouped.yrs=function(DATA,tcl.1,tcl.2,numInt,grouping) 
    {
      DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
      DATA$LONG=as.numeric(DATA$LONG)
      DATA$blk=substr(DATA$BLOCKX,1,4)
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      
      DATA=subset(DATA,FINYEAR%in%FINYrS)
      DATA$yyr=as.numeric(substr(DATA$FINYEAR,1,4))
      
      AA=vector('list',length(FINYrS.gped))
      for(y in 1:length(FINYrS.gped))
      {
        A=aggregate(LIVEWT.c~blk,subset(DATA,FINYEAR%in%FINYrS.gped[[y]]),sum,na.rm=T)
        A$y.group=y
        AA[[y]]=A
      }
      A=do.call(rbind,AA)
      
      Ymax=max(A$LIVEWT.c)
      Ymin=min(A$LIVEWT.c)
      Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
      
      couleurs=rev(heat.colors(numInt))
      
      for(y in 1:length(FINYrS.gped))
      {
        A=subset(DATA,FINYEAR%in%FINYrS.gped[[y]])
        MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
        colnames(MapCatch)=c("BLOCKX.c","Total Catch")
        id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
        MapCatch$LAT=DATA$LAT[id]
        MapCatch$LONG=DATA$LONG[id]
        msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
        msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
        if(length(msn.lat)>0)
        {
          dummy=MapCatch[1:length(msn.lat),]
          dummy$`Total Catch`=0
          dummy$LAT=msn.lat
          dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
          MapCatch=rbind(MapCatch,dummy)
        }
        
        if(unique(min(A$YEAR.c))>2007)
        {
          MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
        }
        MapCatch$LAT.cen=MapCatch$LAT-.5
        MapCatch$LONG.cen=MapCatch$LONG+.5  
        MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
        MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
        long=sort(unique(MapCatch$LONG.cen))
        lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
        MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
        Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                   timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
        Reshaped=Reshaped[order(Reshaped[,1]),]
        Reshaped=Reshaped[,-1]	
        numberLab=10
        colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
        
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
        axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
        axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
        par(new=T)
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        text(115.60,-31.96,"Perth",pos=4,cex=1.2)
        legend('top',names(FINYrS.gped)[y],bty='n',cex=1.2)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      }
      plot(a,b,ann=F,axes=F,col='transparent')
      color.legend(quantile(a,probs=.25),quantile(b,probs=.95),quantile(a,probs=.6),quantile(b,probs=.25),
                   paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.75)
    }
    HnD.ctch.exp="C:/Matias/Analyses/Catch and effort/Outputs/Expansion_Catch"    
    for(i in 1:length(Tar))
    {
      ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                   METHOD%in%c("GN","LL")& LAT<=(-26) & SPECIES%in%Tar[[i]])
      FINYrS=unique(ddd$FINYEAR)
      
      #all years
      fn.fig(paste(HnD.ctch.exp,names(Tar)[i],sep="/"),2400, 2400)
      smart.par(n.plots=length(FINYrS)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
      fn.ctch.plot.all.yrs(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=20)
      mtext("Latitude (?S)",side=2,line=0.5,las=3,cex=1.1,outer=T)
      mtext("Longitude (?E)",side=1,line=0.5,cex=1.1,outer=T)
      dev.off()
      
      #grouped years
      grouping=5
      FINYrS.gp=seq(1,length(FINYrS),by=grouping)
      FINYrS.gped=vector('list',length(FINYrS.gp))
      for(f in 1:length(FINYrS.gped))
      {
        if(f==length(FINYrS.gped))
        {
          FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
          if(length(FINYrS.gped[[f]]==1))
          {
            names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
          }else
          {
            names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
          }
          
        }else
        {
          FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
          names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
        }
      }
      fn.fig(paste(HnD.ctch.exp,"/",names(Tar)[i],"_grouped.yrs",sep=""),2000, 2400)
      smart.par(n.plots=length(FINYrS.gped)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
      fn.ctch.plot.grouped.yrs(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50,grouping=5)
      mtext("Latitude (?S)",side=2,line=0.5,las=3,cex=1.1,outer=T)
      mtext("Longitude (?E)",side=1,line=0.5,cex=1.1,outer=T)
      dev.off()
      
      rm(ddd)
    }
    
    #movie
    Frame.speed=.2  #match to talking speed
    ani.options(ani.width=480,ani.height=480)
    fn.ctch.plot.all.yrs.movie=function(DATA,tcl.1,tcl.2,numInt) 
    {
      DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
      DATA$LONG=as.numeric(DATA$LONG)
      DATA$blk=substr(DATA$BLOCKX,1,4)
      A=aggregate(LIVEWT.c~FINYEAR+blk,DATA,sum)
      Ymax=max(A$LIVEWT.c)
      Ymin=min(A$LIVEWT.c)
      Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
      #Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      saveGIF({
        for(y in 1:length(FINYrS))
        {
          A=subset(DATA,FINYEAR==FINYrS[y])
          MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
          colnames(MapCatch)=c("BLOCKX.c","Total Catch")
          id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
          MapCatch$LAT=DATA$LAT[id]
          MapCatch$LONG=DATA$LONG[id]
          msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
          msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
          if(length(msn.lat)>0)
          {
            dummy=MapCatch[1:length(msn.lat),]
            dummy$`Total Catch`=0
            dummy$LAT=msn.lat
            dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
            MapCatch=rbind(MapCatch,dummy)
          }
          
          if(unique(min(A$YEAR.c))>2007)
          {
            MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
          }
          MapCatch$LAT.cen=MapCatch$LAT-.5
          MapCatch$LONG.cen=MapCatch$LONG+.5  
          MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
          MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
          long=sort(unique(MapCatch$LONG.cen))
          lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
          MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
          Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                     timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
          Reshaped=Reshaped[order(Reshaped[,1]),]
          Reshaped=Reshaped[,-1]	
          numberLab=10
          colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
          
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          image(long,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,axes = FALSE,add=T)			
          axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
          axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
          par(new=T)
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          legend('top',FINYrS[y],bty='n',cex=1.2)
          axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
          axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
          mtext("Latitude (?S)",side=2,line=1,las=3,cex=1.75,outer=T)
          mtext("Longitude (?E)",side=1,line=1,cex=1.75,outer=T)
          color.legend(quantile(a,probs=.9),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.5),
                       paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.95)
        }
      },movie.name=paste(HnD.ctch.exp,paste(names(Tar)[i],".movie.gif",sep=''),sep="/"),interval=Frame.speed,loop =1)
    }
    for(i in 1:length(Tar))
    {
      ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                   METHOD%in%c("GN","LL")& LAT<=(-26) & SPECIES==Tar[[i]])
      FINYrS=unique(ddd$FINYEAR)
      fn.ctch.plot.all.yrs.movie(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=20)
      rm(ddd)
    }
    
    
    
    
    #Effort   
    #km gn days 
    #Monthly
    Attach.Effort.monthly.blk.c=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR+LAT+LONG,
                                          data=subset(Effort.monthly,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T)
    Attach.Effort.monthly.blk.c=fn.split.boundary(Attach.Effort.monthly.blk.c,"Km.Gillnet.Days.c")
    Attach.Effort.monthly.blk.c=aggregate(Km.Gillnet.Days.c~BLOCKX+FINYEAR+LAT+LONG,data=Attach.Effort.monthly.blk.c,sum,na.rm=T)
    
    #Daily
    if(Use.Date=="NO")Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~ID+blockx+vessel+finyear+LAT+LONG,
                                                          data=subset(Effort.daily,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T)
    if(Use.Date=="YES")Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~date+blockx+vessel+finyear+LAT+LONG,
                                                           data=subset(Effort.daily,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T) 
    Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~blockx+finyear+LAT+LONG,data=Attach.Effort.daily.blk.c,sum,na.rm=T)
    names(Attach.Effort.daily.blk.c)=names(Attach.Effort.monthly.blk.c)
    
    Effor.monthly.km.gn.d.block=rbind(Attach.Effort.monthly.blk.c,Attach.Effort.daily.blk.c)
    
    fn.effort.breaks=function(DD)
    {
      Max.eff=Yr.rango
      for (j in 1:length(Yr.rango))
      {
        DATA=subset(DD,FINYEAR%in%Yr.rango[[j]])
        dummy=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX),FUN=sum,na.rm=T))
        Max.eff[[j]]=dummy$x
      }
      Max.eff=unlist(Max.eff)
      Breaks=quantile(Max.eff,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
      return(Breaks)
    }
    Yr.rango=vector('list',length(Yr.group))   
    for(e in 1:length(Yr.rango))
    {
      from=seq(Yr.group[e],Yr.group.plus[e])
      to=substr(seq(Yr.group[e]+1,Yr.group.plus[e]+1),3,4)
      Yr.rango[[e]]=paste(from,to,sep="-") 
    }
    
    EffortBreakS=fn.effort.breaks(Effor.monthly.km.gn.d.block)
    
    
    #km gn days
    # fn.fig(paste(hndl.sptl.ktch,"Spatio.temporal.effort.km.gn.d",sep=""),2400, 2400) 
    # smart.par(n.plots=length(Yr.rango),MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(.1, 0.15, 0))
    # 
    # for (i in 1:length(Yr.rango))
    # {
    #   fn.eff.plot(subset(Effor.monthly.km.gn.d.block,FINYEAR%in%Yr.rango[[i]]),
    #               tcl.1=.5,tcl.2=.5,EffortBreakS)
    #   mtext(Yr.range[i],side=3,line=-1.35,cex=.95)
    #   axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .5,las=1,cex.axis=1,padj=-.15)
    #   axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .5,las=2,cex.axis=1,hadj=1.1)
    #   
    #   if(i==length(Yr.rango))color.legend(126,-26,129,-30.5,paste(round(EffortBreakS,0),"km gn d"),
    #                                       rect.col=couleurs,gradient="y",col=colLeg,cex=0.65)
    # }
    # mtext("Latitude (?S)",side=2,line=0.4,las=3,cex=1.3,outer=T)
    # mtext("Longitude (?E)",side=1,line=0.6,cex=1.3,outer=T)
    # dev.off()
    
    
    
    #Each year
    fn.eff.plot.all.yrs=function(DATA,WHAT)   #function plot effort
    {
      MapEffort=subset(DATA,!is.na(Km.Gillnet.Days.c))
      colnames(MapEffort)=c("FINYEAR","BLOCKX.c","LAT","LONG","Total effort")
      lat=sort(unique(MapEffort$LAT))
      # ids=match(-27,lat)
      # if(is.na(ids))      #add dummy for sorting image out of whack when missing Lat
      # {
      #   adD=MapEffort[1,]
      #   adD[,2]=NA
      #   adD$BLOCKX.c="-27 113"
      #   adD$LAT=-27
      #   MapEffort=rbind(MapEffort,adD)
      # }
      
      if(WHAT=="monthly")
      {
        MapEffort$LAT.cen=MapEffort$LAT-.5
        MapEffort$LONG.cen=MapEffort$LONG+.5 
      }else
      {
        MapEffort$LAT.cen=MapEffort$LAT
        MapEffort$LONG.cen=MapEffort$LONG 
        
      }
      
      MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
      MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapEffort$LONG.cen))
      lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image  
      
      MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Total effort"),names(MapEffort))]  
      #MapEffort=rbind(MapEffort,dummy)#keep in perspective
      
      Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",    #transposed as matrix   
                                 timevar="LAT.cen",v.names="Total effort", direction="wide"))  
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]  									
      
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      #plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      #image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=EffortBreaks,axes = FALSE,add=T)			
      #par(new=T)
      #plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,ylim=c(-36,-26),xlim=c(113,129),breaks=EffortBreaks,axes = FALSE)
      box()
      
    }
    HnD.eff.exp="C:/Matias/Analyses/Catch and effort/Outputs/Expansion_Effort"
    
    #Monthly
    s1=subset(Eff.monthly.c,LAT<=(-26) & !FINYEAR%in%This)
    s1$BLOCKX=substr(s1$BLOCKX,1,4)
    s1=aggregate(Km.Gillnet.Days.c~FINYEAR+BLOCKX+LAT+LONG,s1,sum)
    
    tiff(file=paste(HnD.eff.exp,"Effort dynamics.gillnets_all.yrs_monthly.tiff",sep="/"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=length(Mn.yrs)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
    Mn.yrs=sort(unique(s1$FINYEAR))
    EffortBreaks=quantile(s1$Km.Gillnet.Days.c,probs=seq(0,1,1/numInt),na.rm=T)
    for (i in 1:length(Mn.yrs))
    {
      fn.eff.plot.all.yrs(DATA=subset(s1,FINYEAR== Mn.yrs[i]),WHAT="monthly")
      mtext(Mn.yrs[i],side=3,line=0,cex=.95)
      axis(side = 1, at =Long.seq, labels = F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
    }  
    plot(1:1,col="transparent",axes=F,ylab="",xlab="")
    color.legend(xl=0.72,yb=0.39,xr=1.29,yt=1.4,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
    mtext("Latitude (?S)",side=2,line=.5,las=3,cex=1.25,outer=T)
    mtext("Longitude (?E)",side=1,line=.8,cex=1.25,outer=T)
    dev.off()
    
    
    #Daily
    s2=subset(Eff.daily.c.daily,LAT<=(-26))
    s2$LatDeg=with(s2,as.numeric(substr(block10,1,2)))  
    s2$LatMin=with(s2,10*as.numeric(substr(block10,3,3)))  
    s2$LongDeg=with(s2,100+as.numeric(substr(block10,4,5)))
    s2$LongMin=with(s2,10*as.numeric(substr(block10,6,6)))
    s2$LAT=-with(s2,LatDeg+(LatMin/60))
    s2$LONG=with(s2,LongDeg+(LongMin/60))
    s2=aggregate(Km.Gillnet.Days.c~finyear+block10+LAT+LONG,s2,sum)
    
    Dy.yrs=sort(unique(s2$finyear))
    tiff(file=paste(HnD.eff.exp,"Effort dynamics.gillnets_all.yrs_daily.tiff",sep="/"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=length(Dy.yrs)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
    for (i in 1:length(Dy.yrs))
    {
      EffortBreaks=quantile(s2$Km.Gillnet.Days.c,probs=seq(0,1,1/numInt),na.rm=T)
      fn.eff.plot.all.yrs(DATA=subset(s2,finyear== Dy.yrs[i]),WHAT="daily")
      axis(side = 1, at =Long.seq, labels =F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
      mtext(Dy.yrs[i],side=3,line=0,cex=.95)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
    }  
    plot(1:1,col="transparent",axes=F,ylab="",xlab="")
    color.legend(xl=0.72,yb=0.5,xr=1.29,yt=1.4,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
    mtext("Latitude (?S)",side=2,line=.5,las=3,cex=1.25,outer=T)
    mtext("Longitude (?E)",side=1,line=.8,cex=1.25,outer=T)
    dev.off()
    
    
    #Monthly and years in same plot
    fn.eff.plot.all.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt) 
    {
      DATA$blk=with(DATA,paste(LAT,LONG))
      A=aggregate(km.gillnet.days.c~finyear+blk,DATA,sum)
      Ymax=max(A$km.gillnet.days.c)
      Ymin=min(A$km.gillnet.days.c)
      Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
      #Breaks=seq(Ymin,Ymax,length.out=(numInt+1))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      colfunc <- colorRampPalette(c("yellow", "red"))
      couleurs=colfunc(numInt)
      
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      
      for(y in 1:length(FINYrS))
      {
        A=subset(DATA,finyear==FINYrS[y])
        MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
        colnames(MapEffrt)=c("BLOCKX.c","Effort")
        id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
        MapEffrt$LAT=DATA$LAT[id]
        MapEffrt$LONG=DATA$LONG[id]
        MapEffrt$LAT.cen=MapEffrt$LAT-.5
        MapEffrt$LONG.cen=MapEffrt$LONG+.5  
        MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
        MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
        long=sort(unique(MapEffrt$LONG.cen))
        lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
        MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
        Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                   timevar="LAT.cen",v.names="Effort", direction="wide"))	
        Reshaped=Reshaped[order(Reshaped[,1]),]
        Reshaped=Reshaped[,-1]	
        
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
        axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
        axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
        par(new=T)
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        legend('top',FINYrS[y],bty='n',cex=1.2)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      }
      plot(a,b,ann=F,axes=F,col='transparent')
      color.legend(quantile(a,probs=.8),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.25),
                   paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=.7)
    }
    s1=subset(Eff.monthly.c,LAT<=(-26) &FINYEAR%in% Mn.yrs)
    s1$BLOCKX=substr(s1$BLOCKX,1,4)
    names(s1) =  casefold(names(s1))
    s2=subset(Eff.daily.c,LAT<=(-26) &!finyear%in% Mn.yrs)
    names(s2) =  casefold(names(s2))
    ddd=rbind(s1[,match(names(s2),names(s1))],s2)
    FINYrS=sort(unique(ddd$finyear))
    ddd$LAT=as.numeric(substr(ddd$lat,1,3))
    ddd$LONG=as.numeric(substr(ddd$long,1,3))
    fn.fig(paste(HnD.eff.exp,"Effort dynamics.gillnets_all.yrs",sep="/"),2400, 2400)
    smart.par(n.plots=length(FINYrS)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
    fn.eff.plot.all.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50)
    mtext("Latitude (?S)",side=2,line=0.5,las=3,cex=1.1,outer=T)
    mtext("Longitude (?E)",side=1,line=0.5,cex=1.1,outer=T)
    dev.off()
    
    #grouped years
    fn.eff.plot.grouped.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt,grouping) 
    {
      DATA$blk=with(DATA,paste(LAT,LONG))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      AA=vector('list',length(FINYrS.gped))
      for(y in 1:length(FINYrS.gped))
      {
        A=aggregate(km.gillnet.days.c~blk,subset(DATA,finyear%in%FINYrS.gped[[y]]),sum,na.rm=T)
        A$y.group=y
        AA[[y]]=A
      }
      A=do.call(rbind,AA)
      Ymax=max(A$km.gillnet.days.c)
      Ymin=min(A$km.gillnet.days.c)
      Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
      couleurs=rev(heat.colors(numInt))
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      
      for(y in 1:length(FINYrS.gped))
      {
        A=subset(DATA,finyear%in%FINYrS.gped[[y]])
        MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
        colnames(MapEffrt)=c("BLOCKX.c","Effort")
        id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
        MapEffrt$LAT=DATA$LAT[id]
        MapEffrt$LONG=DATA$LONG[id]
        MapEffrt$LAT.cen=MapEffrt$LAT-.5
        MapEffrt$LONG.cen=MapEffrt$LONG+.5  
        MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
        MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
        long=sort(unique(MapEffrt$LONG.cen))
        lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
        MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
        Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                   timevar="LAT.cen",v.names="Effort", direction="wide"))	
        Reshaped=Reshaped[order(Reshaped[,1]),]
        Reshaped=Reshaped[,-1]	
        
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
        axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
        axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
        par(new=T)
        plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
        text(115.60,-31.96,"Perth",pos=4,cex=1.2)
        legend('top',names(FINYrS.gped)[y],bty='n',cex=1.2)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      }
      plot(a,b,ann=F,axes=F,col='transparent')
      color.legend(quantile(a,probs=.8),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.25),
                   paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=.7)
    }
    grouping=5
    FINYrS.gp=seq(1,length(FINYrS),by=grouping)
    FINYrS.gped=vector('list',length(FINYrS.gp))
    for(f in 1:length(FINYrS.gped))
    {
      if(f==length(FINYrS.gped))
      {
        FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
        if(length(FINYrS.gped[[f]]==1))
        {
          names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
        }else
        {
          names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
        }
        
      }else
      {
        FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
        names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
      }
    }
    
    fn.fig(paste(HnD.eff.exp,"Effort dynamics.gillnets_grouped.yrs",sep="/"),2000, 2400)
    smart.par(n.plots=length(FINYrS.gped)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
    fn.eff.plot.grouped.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50,grouping=grouping)
    mtext("Latitude (?S)",side=2,line=0.5,las=3,cex=1.1,outer=T)
    mtext("Longitude (?E)",side=1,line=0.5,cex=1.1,outer=T)
    dev.off()
    
    #MOVIE_Monthly and years in same plot
    Frame.speed=.2  #match to talking speed
    ani.options(ani.width=480,ani.height=480)
    movie.fn.eff.plot.all.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt) 
    {
      DATA$blk=with(DATA,paste(LAT,LONG))
      A=aggregate(km.gillnet.days.c~finyear+blk,DATA,sum)
      Ymax=max(A$km.gillnet.days.c)
      Ymin=min(A$km.gillnet.days.c)
      Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
      #Breaks=seq(Ymin,Ymax,length.out=(numInt+1))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      colfunc <- colorRampPalette(c("yellow", "red"))
      couleurs=colfunc(numInt)
      
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      saveGIF({
        for(y in 1:length(FINYrS))
        {
          par(las=1,mar=c(1,1,.1,1),oma=c(3,4,.1,.1),mgp=c(3.5,.5,0))
          A=subset(DATA,finyear==FINYrS[y])
          MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
          colnames(MapEffrt)=c("BLOCKX.c","Effort")
          id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
          MapEffrt$LAT=DATA$LAT[id]
          MapEffrt$LONG=DATA$LONG[id]
          MapEffrt$LAT.cen=MapEffrt$LAT-.5
          MapEffrt$LONG.cen=MapEffrt$LONG+.5  
          MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
          MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
          long=sort(unique(MapEffrt$LONG.cen))
          lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
          MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
          Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                     timevar="LAT.cen",v.names="Effort", direction="wide"))	
          Reshaped=Reshaped[order(Reshaped[,1]),]
          Reshaped=Reshaped[,-1]	
          
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
          axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
          axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
          par(new=T)
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          legend('top',FINYrS[y],bty='n',cex=1.75)
          axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
          axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
          mtext("Latitude (?S)",side=2,line=1,las=3,cex=1.75,outer=T)
          mtext("Longitude (?E)",side=1,line=1,cex=1.75,outer=T)
          color.legend(quantile(a,probs=.9),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.5),
                       paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=1)
          
        }
      },movie.name=paste(HnD.eff.exp,"Effort.movie.gif",sep="/"),interval=Frame.speed,loop =1)
      
      
    }
    s1=subset(Eff.monthly.c,LAT<=(-26) &FINYEAR%in% Mn.yrs)
    s1$BLOCKX=substr(s1$BLOCKX,1,4)
    names(s1) =  casefold(names(s1))
    s2=subset(Eff.daily.c,LAT<=(-26) &!finyear%in% Mn.yrs)
    names(s2) =  casefold(names(s2))
    ddd=rbind(s1[,match(names(s2),names(s1))],s2)
    FINYrS=sort(unique(ddd$finyear))
    ddd$LAT=as.numeric(substr(ddd$lat,1,3))
    ddd$LONG=as.numeric(substr(ddd$long,1,3))
    movie.fn.eff.plot.all.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50)
    
  }
  
}

#------ G 2. DATA REQUESTS ------
hndl="C:/Matias/Analyses/Catch and effort/Data_Resquests"

#G 4.1 AUDITS
if(do.audit=="YES")
{
  Audit=subset(Data.daily,FINYEAR==Current.yr & METHOD%in%c("GN","LL") & Estuary=="NO" &
                 LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  write.csv(Audit,paste(hndl,"/",Current.yr,".Audit.csv",sep=""),row.names=F)  
}

#G 4.2. Russel Hudson, FishWell consulting
if(do.Russels=="YES")
{
  Russel=subset(Data.monthly,SPECIES%in%c(18032,19000),select=c(LIVEWT,FINYEAR,SPECIES,SNAME))
  Russel$Fishery="TDGDLF"
  Russel2=subset(Data.monthly.north,SPECIES%in%c(18032,19000),select=c(LIVEWT,FINYEAR,SPECIES,SNAME))
  Russel2$Fishery="NSH"
  Russel=rbind(Russel,Russel2)
  Russel.catch=aggregate(LIVEWT~FINYEAR+SPECIES+SNAME+Fishery,Russel,sum)  
  Russel.catch=Russel.catch[order(Russel.catch$SPECIES,Russel.catch$FINYEAR,Russel.catch$Fishery),]
  write.csv(Russel.catch,paste(hndl,"/Russel.catch.csv",sep=""),row.names=F)  
}

#G 4.3. Carly, TEPS
#note: check for Tiger sharks in Comments and inconsistencies in the numbers reported in the
#       columns A, D, and in the comments
Carly.yr=2015   #calendar year
if(do.Carlys=="YES")
{
  These.TEPS=TEPS.current
  TEPS.Carly=subset(These.TEPS,fishery%in%c("SGL","WCGL") & year==Carly.yr)
  TEPS.Carly$LAT=-as.numeric(substr(TEPS.Carly$block10,1,2))
  TEPS.Carly$LONG=100+as.numeric(substr(TEPS.Carly$block10,4,5))
  TEPS.Carly$Bioregion=as.character(with(TEPS.Carly,ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
                                                           ifelse(LONG<115.5 & LAT<=(-27),"WC",
                                                                  ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
                                                                         ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))))
  TEPS.Carly$Bioregion=with(TEPS.Carly,
                            ifelse(Bioregion=="SC"& LAT>(-34) & LONG <115.91 ,"WC",Bioregion))
  
  
  Scan.Carly.teps=TEPS.Carly[,match(c("DailySheetNumber","vessel","month",
                                      "Status","Number","DataEntryName","Comments"),names(TEPS.Carly))]
  write.csv(Scan.Carly.teps,"Scan.Carly.teps.csv",row.names=F)
  Scan.Carly.teps$ID=1:nrow(Scan.Carly.teps)
  #check.comments=TEPS.Carly[c(60:71),]  #these bronzies recorded in methods but not in numbers for 2012
  check.comments$Unic=with(check.comments,paste(DailySheetNumber,SessionNumber))
  check.comments=check.comments[!duplicated(check.comments$Unic),]
  
  if(Carly.yr==2012)
  {
    #Change 'Whale' record entry to 'DUSKY WHALER OVERSIZE' as per comments.
    TEPS.Carly$SpeciesCode=with(TEPS.Carly,ifelse(DataEntryName=="WHALES",18003,SpeciesCode))
    TEPS.Carly$DataEntryName=with(TEPS.Carly,
                                  ifelse(DataEntryName=="WHALES","DUSKY WHALER OVERSIZE",DataEntryName)) 
  }
  
  fun.Tab2.SoFaR=function(Dat)
  {
    Dat$Status=as.character(Dat$Status)
    Dat$Status=with(Dat,ifelse(Status=="a","A",
                               ifelse(Status=="d","D",Status)))
    Dat$Ali.Ded.yr=with(Dat,paste(finyear,Status))
    
    TABLA=aggregate(Number~SpeciesCode+Ali.Ded.yr+Bioregion,data=Dat,sum,na.rm=T)
    Species.id=Dat[,match(c("SpeciesCode","DataEntryName"),names(Dat))]
    Species.id=Species.id[!(duplicated(Species.id$SpeciesCode)),]
    
    wide <- reshape(TABLA,v.names="Number",timevar="Ali.Ded.yr",idvar=c("SpeciesCode","Bioregion"),direction="wide")
    wide=wide[order(wide$SpeciesCode),]
    wide=Species.id%>%left_join(wide,by="SpeciesCode")
    #wide=merge(Species.id,wide,by="SpeciesCode",all.x=T)
    wide=wide[-match("SpeciesCode",names(wide))]
    names(wide)=c("Species","Bioregion","alive","dead")
    wide=wide[order(wide$Bioregion,wide$Species),]
    wide=wide[,match(c("Bioregion","Species","alive","dead"),names(wide))]
    return(wide)
  }
  TEP.Carly=fun.Tab2.SoFaR(TEPS.Carly)
  
  
  if(Carly.yr==2015)
  { 
    #manually add records in comments for 1 tiger DailySheetN TDGLF6000137
    Add.Carly=TEP.Carly[3,]
    Add.Carly$Species="TIGER SHARK OVERSIZE"
    Add.Carly$alive=1
    Add.Carly$dead=NA  
    TEP.Carly=rbind(TEP.Carly,Add.Carly) 
  }
  
  if(Carly.yr==2013)
  { 
    #manually add records in comments for 1 tiger DailySheetN TDGLF6000137
    Add.Carly=TEP.Carly[3,]
    Add.Carly$Species="TIGER SHARK OVERSIZE"
    Add.Carly$alive=1
    Add.Carly$dead=NA  
    TEP.Carly=rbind(TEP.Carly,Add.Carly)   
  }
  
  TEP.Carly=TEP.Carly[order(TEP.Carly$Bioregion,TEP.Carly$Species),]
  
  write.csv(TEP.Carly,paste(hndl,"/Carly/Carly.TEPS.",Carly.yr,".csv",sep=""),row.names=F)
  
}
if(exists("TEPS.2011_12"))
{
  #note: This should be added to catch information for assessment. 
  #     For finyears post 2011-12, must manually scan comments for matching what was reported in columns and in comments
  
  # Also note that for most whalers, (e.g. Tigers), comments don't say if oversized or not. So 
  #assumed they are
  Lista.protected.elasmos=list(
    Tiger=c("tiger","Tiger","TIGER","tigersharks","TG","T/S","Tidgers"),
    Dusky=c("dusky","Bronze","bronze","b/w","Dusky","whalers","BW","Bronzy","b/whaler"),
    Copper=c("copper","Copper"),
    Sandbar=c("Thiskskin","Sandbar","thickskins","TK"),
    Greynurse=c("Nurse","NURSE","nurse","g/n","GN","greynurse","G.N","G/N"),
    Blacktip=c("Blacktip"),
    Whites=c("W/S","WHITE","white","White","Pointer"),
    Mako=c("Mako"),
    Sawfish=c("Sawfish"),
    Manta=c("Manta"))
  
  
  Lista.id.pro.el=Check.These.Sks=vector('list',length=length(Lista.protected.elasmos))
  names(Lista.id.pro.el)=names(Check.These.Sks)=names(Lista.protected.elasmos)
  
  fn.find=function(what,inwhat)
  {
    Find=regexpr(what, inwhat$Comments) > 0
    return(which(Find==T))
  }
  
  for (i in 1:length(Lista.protected.elasmos))
  {
    these=Lista.protected.elasmos[[i]]
    store=NULL
    if(length(these)>1)for(j in 1:length(these))store=c(store,fn.find(these[j],TEPS.2011_12))    
    if(length(these)==1)store=fn.find(these,TEPS.2011_12)
    Lista.id.pro.el[[i]]=store
  }
  
  for (i in 1:length(Lista.id.pro.el))
  {
    id=Lista.id.pro.el[[i]]
    a=TEPS.2011_12[id,match(c("DailySheetNumber","SessionNumber",
                              "DataEntryName","Status","Number","Comments"),names(TEPS.2011_12))]
    a=a[order(a$DailySheetNumber,a$SessionNumber),]
    Check.These.Sks[[i]]=a
  }
  
  #The only manual bit is to go thru each element of Check.These.Sks and compare comments with data
  i=1;print(Check.These.Sks[[i]])
  
  #Table
  All.TEP.table=fun.Tab2.SoFaR(TEPS.2011_12)
  TEPS[is.na(TEPS)] = ""
}

#G 4.4 Alex Hexp. Effort of ASL model
if(do.Alexs=="YES")
{
  Data.daily.Alex=merge(Data.daily.Alex,BlOCK_10,by.x="block10",by.y="BlockNo",all.x=T)
  Effort.alex.daily=Effort.alex.daily[!duplicated(Effort.alex.daily$Same.return.SNo),]
  Effort.alex.daily$Estuary=with(Effort.alex.daily,ifelse(blockx%in%Estuaries,"YES","NO"))
  
  #Effort.alex.daily=subset(Effort.alex.daily,Estuaries=="NO")
  
  KEEp.nms=c(match("Same.return.SNo",names(Data.daily.Alex)),which(!names(Data.daily.Alex)%in%names(Effort.alex.daily)))
  Data.Alex=merge(Data.daily.Alex[,KEEp.nms],Effort.alex.daily,by="Same.return.SNo",all.y=T)
  Data.Alex=Data.Alex[!(duplicated(Data.Alex$Same.return.SNo)),]
  Effort.alex=subset(Effort.alex,YEAR.c>=2005)
  Data.Alex=subset(Data.Alex,year>=2005 & LatDeg>=26)
  Data.Alex=subset(Data.Alex,netlen.c>100 & method=="GN")

  Data.Alex.nlines.check=subset(Data.Alex,select=c(Same.return.SNo,vessel,
                       netlen,netlen.c,hours,hours.c,nlines,shots,shots.c,fishery))
  Data.Alex=subset(Data.Alex,select=c(Same.return.SNo,finyear,year,month,day,date,
                  blockx,block10,netlen,netlen.c,hours.c,nlines,nlines.c,shots.c,Long,
                  LongDeg,LongMin,Lat,
                  LatDeg,LatMin,Latitude_Centroid,Longitude_Centroid,depthMax,fishery,vessel))
   
  change.names=match(c("Same.return.SNo","blockx","date","netlen","netlen.c","hours.c",
                       "nlines","nlines.c","shots.c"),names(Data.Alex))
  names(Data.Alex)[change.names]=c("SessionID","block","Date","netlen.original","netlen","hours",
                                   "nlines.original","nlines","shots")
  This.Alex=c("vessel","SessionID","fishery","finyear","month","year","Date","netlen.original","netlen",
              "hours","nlines.original","nlines","shots","block10","block","Lat","LatDeg","LatMin",
              "Long","LongDeg","LongMin","Latitude_Centroid","Longitude_Centroid","depthMax")
  Data.Alex=Data.Alex[,match(This.Alex,names(Data.Alex))]
  Data.Alex$Date=as.character(Data.Alex$Date)
  
  fwrite(Data.Alex,paste(hndl,"/Alex_ASL/Data.daily.Alex.csv",sep=""),row.names=F)  
  fwrite(Data.Alex.nlines.check,paste(hndl,"/Alex_ASL/Data.Alex.nlines.check.csv",sep=""),row.names=F)
  fwrite(Effort.alex,paste(hndl,"/Alex_ASL/Data.CAESS.Alex.csv",sep=""),row.names=F)
  
  rm(Data.daily.Alex,Data.Alex)
  
}

#G 4.5 SAFS. Number of vessels catching a target species
if(do.Garys=="YES")
{
  fn.gary=function(SPEC,YR)
  {
    dat=subset(Data.monthly,SPECIES==SPEC & FINYEAR==YR)
    dat$zoneVes=with(dat,paste(zone,VESSEL))
    unik.ves=unique(dat$VESSEL)
    dat=dat[!duplicated(dat$zoneVes),]
    unik.ves.zn=table(dat$zone)
    
    return(list(Numb.ves.all=length(unik.ves),Numb.ves.zn=unik.ves.zn))
  }
  
  Pinkies=fn.gary(353001,"2012-13")
  Pinkies=unlist(Pinkies)
  names.Pinkies=names(Pinkies)
  Pinkies=matrix(Pinkies,ncol=4)
  colnames(Pinkies)=names.Pinkies
  
  write.csv(Pinkies,paste(hndl,"/Pinkies.for.gary.csv",sep=""),row.names=F)
  
}

#G 4.6 Whiskery shark catch by month (AMM 2014)
if(do.whiskerys=="YES")
{
  Yr.Whiskery=2005:2013
  Whisk.ktch.West=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="West"),sum)
  Whisk.ktch.Zn1=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="Zone1"),sum)
  Whisk.ktch.Zn2=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="Zone2"),sum)
  
  closure.ktch=data.frame(YEAR.c=2008:2012,MONTH=9,LIVEWT.c=0)
  Whisk.ktch.Zn1=rbind(Whisk.ktch.Zn1,closure.ktch)
  Whisk.ktch.West=rbind(Whisk.ktch.West,closure.ktch)
  
  Whisk.ktch.West=Whisk.ktch.West[order(Whisk.ktch.West$YEAR.c,Whisk.ktch.West$MONTH),]
  Whisk.ktch.Zn1=Whisk.ktch.Zn1[order(Whisk.ktch.Zn1$YEAR.c,Whisk.ktch.Zn1$MONTH),]
  Whisk.ktch.Zn2=Whisk.ktch.Zn2[order(Whisk.ktch.Zn2$YEAR.c,Whisk.ktch.Zn2$MONTH),]
  
  tiff(file=paste(hndl,"/Whiskery.monthly.catch.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  
  par(las=1)
  MAX.y=max(Whisk.ktch.West$LIVEWT.c/1000,Whisk.ktch.Zn1$LIVEWT.c/1000,Whisk.ktch.Zn2$LIVEWT.c/1000)
  plot(Whisk.ktch.West$LIVEWT.c/1000,ylim=c(0,MAX.y),xaxt='n',col=2,lwd=2,
       type='l',xlab="",ylab="",cex.lab=2,cex.axis=1.5)
  Sep=seq(9,nrow(Whisk.ktch.West),12)
  abline(v=Sep,lwd=3,col="grey80")
  lines(Whisk.ktch.West$LIVEWT.c/1000,col=2,lwd=2)
  lines(Whisk.ktch.Zn1$LIVEWT.c/1000,col=3,lwd=2)
  lines(Whisk.ktch.Zn2$LIVEWT.c/1000,col=4,lwd=2)
  legend("topright",c("West coast","Zone 1","Zone 2"),
         lty=1,col=2:4,lwd=2,cex=2,bg="white",box.col="white")
  axis(1,at=1:nrow(Whisk.ktch.West),F,tck=-0.015)
  YERS=seq(1,nrow(Whisk.ktch.West),12)
  axis(1,at=YERS,Yr.Whiskery,tck=-0.03,cex.axis=1.5)
  box()
  mtext("Year",1,3,cex=2)
  mtext("Whiskery catch (tonnes)",2,2.5,cex=2,las=3)
  dev.off()
  
  names(Whisk.ktch.West)[3]="West"
  names(Whisk.ktch.Zn1)[3]="Zone1"
  names(Whisk.ktch.Zn2)[3]="Zone2"
  Whisk.ktch=merge(Whisk.ktch.West,Whisk.ktch.Zn1,by=c("YEAR.c","MONTH"))
  Whisk.ktch=merge(Whisk.ktch,Whisk.ktch.Zn2,by=c("YEAR.c","MONTH"))
  
  write.csv(Whisk.ktch,paste(hndl,"/Whisk.ktch.csv",sep=""),row.names=F)  
}

#G 4.7 Jeff Norris, all catch and effort for South Coast Bioregion
if(do.Jeffs=="YES")
{
  Yrs.of.scalies.data=unique(subset(Data.monthly,SPECIES%in%Scalefish.species,select=FINYEAR))
  Jeff.fn=function(DATA,Eff.d,Eff.m)
  {
    DATA=subset(DATA,select=c(FINYEAR,YEAR,MONTH,VESSEL,METHOD,BLOCKX,SPECIES,SNAME,LAT,LONG,
                              Bioregion,zone,LIVEWT.c))
    SC.blocks=unique(DATA$BLOCKX)
    Eff.m=subset(Eff.m,BLOCKX%in%SC.blocks)  
    Eff.d=subset(Eff.d,blockx%in%SC.blocks)
    
    #Monthly
    Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~VESSEL+FINYEAR+MONTH+BLOCKX,data=Eff.m,max,na.rm=T)
    Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~FINYEAR+MONTH+VESSEL+BLOCKX,data=Eff.monthly.hour.c,sum,na.rm=T) 
    
    #Daily
    Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~ID+vessel+finyear+month+blockx,data=Eff.d,max,na.rm=T)
    Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~finyear+month+vessel+blockx,data=Eff.daily.hour.c,sum,na.rm=T)
    
    names(Eff.daily.hour.c)[1:4]=c("FINYEAR","MONTH","VESSEL","BLOCKX")
    
    return(list(Catch=DATA,Effort=rbind(Eff.monthly.hour.c,Eff.daily.hour.c)))
  }
  Jeff.s=Jeff.fn(subset(Data.monthly,Bioregion=="SC" & FINYEAR%in%Yrs.of.scalies.data$FINYEAR),
                 Effort.daily,subset(Effort.monthly,!FINYEAR%in%Daily.l.years))
  
  write.csv(Jeff.s$Catch,paste(hndl,"/Jeff.Catch.south_coast.csv",sep=""),row.names=F)
  write.csv(Jeff.s$Effort,paste(hndl,"/Jeff.Effort.south_coast.csv",sep=""),row.names=F)
  
}

#G 4.8 Jodie O'Malley   

  #G 4.8.1 ASL exclusion areas catch and effort
hndls=paste(hndl,"/Jodie_OMalley_2016/data/",sep="")
if(do.Jodies.ASL=="YES")
{  
  ASL_exclusions_block10_JASDGDLF=list(
      Twilight_cove=c(321255,321260,321261,322255,322260),
      Esperance=c(335222:335225,335230,335232,333235,333240,333241,
                334235,334240,334241,335235,335240,335241,340215,
                340220:340225,340230:340235,340240,341215,341220:341225,
                341230:341235,341240,342215,342220,342232:342234,
                343215,343220),
      Albany_east=c(343184,344183,344184),
      Hopetoun_west=c(340194,340195,341193,341194,342192:342194),
      Hopetoun_east=c(335203:335205,340201:340205,340210,341202:341205))
  
  ASL_exclusions_block10_WCDGDLF=list(
      Geraldton=c(283133:283135,283140,284133:284135,284140,285134,285135,
                285140,290134,290135,290140),
      JurienBay_north=c(293144,293145,294144,294145,295144,295145,300144,
                300145,301144,301145),
      JurienBay_south=c(303145,303150,304145,304150))
  
  
  Data.daily$Fishing_yr=with(Data.daily,ifelse(MONTH>=6,YEAR.c,YEAR.c-1))
  Data.daily$FISHERY=with(Data.daily,ifelse(zone=="West","WCDGDLF",
                                            ifelse(zone%in%c("Zone1","Zone2"),"JASDGDLF",NA)))
  Jodie.Effort$FISHERY=with(Jodie.Effort,ifelse(zone=="West","WCDGDLF",
                                            ifelse(zone%in%c("Zone1","Zone2"),"JASDGDLF",NA)))
  
  #add zone 3
  Data.daily$zone3=with(Data.daily,ifelse(Data.daily$LONG>=116.5 & Data.daily$LONG<= 116.923 & Data.daily$LAT < (-33),"Zone3",zone))
  Zn3.blks=unique(subset(Data.daily,zone3=="Zone3",select=block10))$block10
  Jodie.Effort$zone3=with(Jodie.Effort,ifelse(block10%in%Zn3.blks,"Zone3",zone))
  
  
  #some functions
  fn.extract.Jodie.catch=function(blk,FSHRY,YRS,SP,Byzone)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP))
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    
    if(Byzone=="NO")
    {
      Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+VESSEL+block10+Spec.tab1,dat.blk,sum)
      Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,dat.fishery,sum)
    }
    
    if(Byzone=="YES")
    {
      Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+VESSEL+block10+Spec.tab1+zone3,dat.blk,sum)
      Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,dat.fishery,sum) 
    }
    
    return(list(BLKs=Aggregated.blk,Fishery=Aggregated.fishery))
  }
  
  fn.extract.Jodie.effort=function(blk,FSHRY,YRS,Byzone)
  {
    Dat=subset(Jodie.Effort)
    
    dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    
    if(Byzone=="NO")
    {
      Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+block10,dat.blk,sum)
      Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    }
    
    if(Byzone=="YES")
    {
      Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+block10+zone3,dat.blk,sum)
      Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+zone3,dat.fishery,sum)     
    }
    
    return(list(BLKs=Aggregated.blk_km.gn.hr,Fishery=Aggregated.fishery_km.gn.hr))
  }
  
  fn.agg.Jodie.efrt=function(dat,Byzone)
  {
    BLK=dat$BLKs
    All=dat$Fishery
    
    if(Byzone=="NO")Blk.ag=aggregate(Km.Gillnet.Hours.c~Fishing_yr,BLK,sum)
    if(Byzone=="YES")Blk.ag=aggregate(Km.Gillnet.Hours.c~Fishing_yr+zone3,BLK,sum)
    
    return(list(BLK=Blk.ag,ALL=All))
  }
  
  fn.tab.eff=function(dat,Byzone)
  {
    if(Byzone=="NO")
    {
      D=merge(dat$ALL,dat$BLK,by="Fishing_yr")
      names(D)[2:3]=c("All","Blocks") 
      D$Prop=D$Blocks/D$All
    }
    
    if(Byzone=="YES")
    {
      D=merge(dat$ALL,dat$BLK,by=c("Fishing_yr","zone3"),all.x=T)
      names(D)[3:4]=c("All","Blocks")    
      D$Prop=D$Blocks/D$All
      
      D=reshape(D, v.names = c("All","Blocks","Prop"), idvar = "Fishing_yr",
              timevar = "zone3", direction = "wide")
      D=subset(D,select=c(Fishing_yr,All.Zone1,All.Zone2,All.Zone3,Blocks.Zone2,Prop.Zone2))
      
    }
     
    
    return(D)
  }
  
  combo.tbl_eff=function(A,B,C,nms)
  {
    D=merge(A,B,by="Fishing_yr")
    D=merge(D,C,by="Fishing_yr")
    D=D[,match(nms,colnames(D))]
    return(D)
  }
  
  fn.agg.Jodie.ktch=function(dat,SORT,Byzone)
  {
    BLK=dat$BLKs
    All=dat$Fishery
    SORT=SORT[which(SORT%in%unique(All$Spec.tab1))]
    
    
    if(Byzone=="NO")
    {
      BLK$Spec.tab1=as.character(BLK$Spec.tab1)
      All$Spec.tab1=as.character(All$Spec.tab1)
      
      Blk.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,BLK,sum)
      Blk.ag.wide<- reshape(Blk.ag,v.names="LIVEWT.c",timevar="Fishing_yr",idvar="Spec.tab1",direction="wide")
      
      All.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,All,sum)
      All.ag.wide<- reshape(All.ag,v.names="LIVEWT.c",timevar="Fishing_yr",idvar="Spec.tab1",direction="wide")
      
      Add.this=All.ag.wide$Spec.tab1[which(!All.ag.wide$Spec.tab1%in%Blk.ag.wide$Spec.tab1)]
      if(length(Add.this)>0)
      {
        addd=Blk.ag.wide[1:length(Add.this),] 
        addd[,]=NA
        addd$Spec.tab1=Add.this
        Blk.ag.wide=rbind(Blk.ag.wide,addd)    
      }
      
      Blk.ag.wide=Blk.ag.wide[match(SORT,Blk.ag.wide$Spec.tab1),]
      All.ag.wide=All.ag.wide[match(SORT,All.ag.wide$Spec.tab1),]
      
      fn.rshp=function(x)
      {
        x=reshape(x,direction="long")
        x=x[order(x$Fishing_yr),]
      }
      Blk.ag.wide=fn.rshp(Blk.ag.wide)
      All.ag.wide=fn.rshp(All.ag.wide)    
    }
    
    if(Byzone=="YES")
    {
      BLK$Spec.tab1=as.character(BLK$Spec.tab1)
      All$Spec.tab1=as.character(All$Spec.tab1)
      
       Blk.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,BLK,sum)  #no need for zone here, only zone2 has closures
       Blk.ag.wide<- reshape(Blk.ag,v.names="LIVEWT.c",timevar="zone3",idvar=c("Fishing_yr","Spec.tab1"),direction="wide")


      All.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,All,sum)
      All.ag.wide<- reshape(All.ag,v.names="LIVEWT.c",timevar=c("zone3"),idvar=c("Fishing_yr","Spec.tab1"),direction="wide")
      
      Add.this=merge(All.ag.wide,Blk.ag.wide,by=c("Fishing_yr","Spec.tab1"),all.x=T)
      Add.this=Add.this[which(!with(Add.this,paste(Spec.tab1,Fishing_yr))%in%with(Blk.ag.wide,paste(Spec.tab1,Fishing_yr))),]
      if(nrow(Add.this)>0)
      {
        addd=Blk.ag.wide[1:nrow(Add.this),] 
        addd[,]=NA
        addd$Fishing_yr=Add.this$Fishing_yr
        addd$Spec.tab1=Add.this$Spec.tab1
        Blk.ag.wide=rbind(Blk.ag.wide,addd)    
      }
      Yrs=unique(Blk.ag.wide$Fishing_yr)
      Yrs.SORT=paste(rep(SORT,each=length(Yrs)),rep(Yrs,times=length(SORT)))
      id.blk=match(Yrs.SORT,paste(Blk.ag.wide$Spec.tab1,Blk.ag.wide$Fishing_yr))
      Blk.ag.wide=Blk.ag.wide[id.blk,]
      id.all=match(Yrs.SORT,paste(All.ag.wide$Spec.tab1,All.ag.wide$Fishing_yr))
      All.ag.wide=All.ag.wide[id.all,]

      if(length(Yrs.SORT[which(is.na(id.blk))])>0)
      {
        Blk.ag.wide$Fishing_yr[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.blk))], " "), "[[", 2)
        Blk.ag.wide$Spec.tab1[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.blk))], " "), "[[", 1)
      }
      if(length(Yrs.SORT[which(is.na(id.all))])>0)
      {
        All.ag.wide$Fishing_yr[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.all))], " "), "[[", 2)
        All.ag.wide$Spec.tab1[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.all))], " "), "[[", 1)
      }

    }
    
    return(list(BLK=Blk.ag.wide,ALL=All.ag.wide))
  }
  
  combo.tbl=function(A,B,C,nms,SORT) 
  {
    if(nrow(A$BLK)==sum(A$BLK$Spec.tab1==A$ALL$Spec.tab1)) A1=cbind(A$ALL,A$BLK$LIVEWT.c)
    colnames(A1)[3:4]=c("WCDGDLF","Closures_WCDGDLF")
    
    if(nrow(B$BLK)==sum(B$BLK$Spec.tab1==B$ALL$Spec.tab1)) B1=cbind(B$ALL,B$BLK$LIVEWT.c)
    colnames(B1)[3:4]=c("JASDGDLF","Closures_JASDGDLF")
    
    if(nrow(C$BLK)==sum(C$BLK$Spec.tab1==C$ALL$Spec.tab1)) C1=merge(C$ALL,C$BLK,by=c("Fishing_yr","Spec.tab1"),all.x=T) 
    colnames(C1)[3:6]=c("JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3","Closures_JASDGDLF")
    
    A1[is.na(A1)]=0
    B1[is.na(B1)]=0
    C1[is.na(C1)]=0
    
    A1$Prop_WCDGDLF=with(A1,Closures_WCDGDLF/WCDGDLF)
    B1$Prop_JASDGDLF=with(B1,Closures_JASDGDLF/JASDGDLF)
    C1$Prop_JASDGDLF_Zone2=with(C1,Closures_JASDGDLF/JASDGDLF_Zone2)
    
    D=merge(A1,B1,by=c("Spec.tab1","Fishing_yr"),all=T)
    D=merge(D,C1,by=c("Spec.tab1","Fishing_yr","Closures_JASDGDLF"),all=T)
    D[is.na(D)]=0
    D$TDGDLF=with(D,WCDGDLF+JASDGDLF)
    D$Closures_TDGDLF=with(D,Closures_WCDGDLF+Closures_JASDGDLF)
    D$Prop_TDGDLF=with(D,Closures_TDGDLF/TDGDLF)                       
    
    D=D[,match(nms,colnames(D))]
    
    SORT=SORT[which(SORT%in%unique(D$Spec.tab1))]
    Yrs=unique(D$Fishing_yr)
    Yrs.SORT=paste(rep(SORT,each=length(Yrs)),rep(Yrs,times=length(SORT)))
    D=D[match(Yrs.SORT,paste(D$Spec.tab1,D$Fishing_yr)),]
    
    
    return(D)
  }
  
  combo.tbl_agg=function(A)
  {
    d=aggregate(cbind(TDGDLF,WCDGDLF,JASDGDLF,
                      JASDGDLF_Zone1,JASDGDLF_Zone2,JASDGDLF_Zone3,
                      Closures_TDGDLF,Closures_WCDGDLF,Closures_JASDGDLF)~Fishing_yr,A,sum)
    d$Prop_TDGDLF=d$Closures_TDGDLF/d$TDGDLF
    d$Prop_WCDGDLF=d$Closures_WCDGDLF/d$WCDGDLF
    d$Prop_JASDGDLF=d$Closures_JASDGDLF/d$JASDGDLF
    d$Prop_JASDGDLF_Zone2=d$Closures_JASDGDLF/d$JASDGDLF_Zone2
    return(d)  
  }
  
  
  #loop over requested periods
  Yrs.Jodie=list(year_2010_14=2010:2014,year_2007_09=2007:2009)  
  
  for(j in 1:length(Yrs.Jodie))
  {
    
    #Extracth catch
    KTCH_WCDGDLF_elasmo=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"NO")
    KTCH_WCDGDLF_scalies=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"NO")
    KTCH_JASDGDLF_elasmo=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"NO")
    KTCH_JASDGDLF_scalies=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"NO")
    
    KTCH_JASDGDLF_elasmo_zone=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"YES")
    KTCH_JASDGDLF_scalies_zone=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"YES")
    
    
    #Extract effort
    EFFORT_WCDGDLF=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],"NO")
    EFFORT_JASDGDLF=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],"NO")
    
    EFFORT_JASDGDLF_zone=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],"YES")
    
    
    #Tabulate effort  
    EFRT.Jodie_WCDGDLF=fn.agg.Jodie.efrt(EFFORT_WCDGDLF,"NO")
    EFRT.Jodie_JASDGDLF=fn.agg.Jodie.efrt(EFFORT_JASDGDLF,"NO")
    EFRT.Jodie_JASDGDLF_zone=fn.agg.Jodie.efrt(EFFORT_JASDGDLF_zone,"YES")
    
    Tab.eff.Jodie_WCDGDLF=fn.tab.eff(EFRT.Jodie_WCDGDLF,"NO")
    Tab.eff.Jodie_JASDGDLF=fn.tab.eff(EFRT.Jodie_JASDGDLF,"NO")
    Tab.eff.Jodie_JASDGDLF_zone=fn.tab.eff(EFRT.Jodie_JASDGDLF_zone,"YES")
    
    names(Tab.eff.Jodie_WCDGDLF)[2:4]=c("WCDGDLF","Closures_WCDGDLF","Prop_WCDGDLF")
    names(Tab.eff.Jodie_JASDGDLF)[2:4]=c("JASDGDLF","Closures_JASDGDLF","Prop_JASDGDLF")
    names(Tab.eff.Jodie_JASDGDLF_zone)[2:6]=c("JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3","Closures_JASDGDLF_Zone2","Prop_JASDGDLF_Zone2")
    Tab.eff.Jodie_JASDGDLF=merge(Tab.eff.Jodie_JASDGDLF,subset(Tab.eff.Jodie_JASDGDLF_zone,
           select=c(Fishing_yr,JASDGDLF_Zone1,JASDGDLF_Zone2,JASDGDLF_Zone3,Prop_JASDGDLF_Zone2)),by="Fishing_yr")
    
    Tab.eff.Jodie=merge(Tab.eff.Jodie_WCDGDLF,Tab.eff.Jodie_JASDGDLF,by="Fishing_yr")  
    Tab.eff.Jodie.TDGDLF=Tab.eff.Jodie[,1:2]
    names(Tab.eff.Jodie.TDGDLF)[2]="TDGDLF"
    Tab.eff.Jodie.TDGDLF$TDGDLF=Tab.eff.Jodie_WCDGDLF$WCDGDLF+Tab.eff.Jodie_JASDGDLF$JASDGDLF
    Tab.eff.Jodie.TDGDLF$Closures_TDGDLF=Tab.eff.Jodie_WCDGDLF$Closures_WCDGDLF+Tab.eff.Jodie_JASDGDLF$Closures_JASDGDLF
    Tab.eff.Jodie.TDGDLF$Prop_TDGDLF=with(Tab.eff.Jodie.TDGDLF,Closures_TDGDLF/TDGDLF)
    Jodie.Eff=combo.tbl_eff(A=Tab.eff.Jodie.TDGDLF,B=Tab.eff.Jodie_WCDGDLF,C=Tab.eff.Jodie_JASDGDLF,
                            nms=c("Fishing_yr","TDGDLF","WCDGDLF","JASDGDLF",
                                  "JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3",
                                  "Closures_TDGDLF","Closures_WCDGDLF","Closures_JASDGDLF",
                                  "Prop_TDGDLF","Prop_WCDGDLF","Prop_JASDGDLF","Prop_JASDGDLF_Zone2"))
    
    Jodie.Eff[,11:14]=100*Jodie.Eff[,11:14]   #percentage
    Jodie.Eff[,2:14]=round(Jodie.Eff[,2:14],1) #round
    
    write.csv(Jodie.Eff,paste(hndls,"Jodie.Eff_km.gn.hr_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    
    
    #tabulate catch
    SORT.s=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads","Spinner","Wobbegongs","Rays","Common_saw_shark","School",
           "Other_elasmobranchs")
    KTCH.Jodie_WCDGDLF_elasmo=fn.agg.Jodie.ktch(KTCH_WCDGDLF_elasmo,SORT.s,"NO")
    KTCH.Jodie_JASDGDLF_elasmo=fn.agg.Jodie.ktch(KTCH_JASDGDLF_elasmo,SORT.s,"NO")
    KTCH.Jodie_JASDGDLF_elasmo_zone=fn.agg.Jodie.ktch(KTCH_JASDGDLF_elasmo_zone,SORT.s,"YES")
    
    SORT.t=c("Blue_morwong","Blue_groper","West_Australian_dhufish","Pink_snapper","Boarfishes","Samsonfish","Redfishes","Mulloway","Sweetlips",
           "Baldchin_groper","Other_scalefish")
    KTCH.Jodie_WCDGDLF_scalies=fn.agg.Jodie.ktch(KTCH_WCDGDLF_scalies,SORT=SORT.t,"NO")
    KTCH.Jodie_JASDGDLF_scalies=fn.agg.Jodie.ktch(KTCH_JASDGDLF_scalies,SORT.t,"NO")
    KTCH.Jodie_JASDGDLF_scalies_zone=fn.agg.Jodie.ktch(KTCH_JASDGDLF_scalies_zone,SORT.t,"YES")
    
    NMS=c("Fishing_yr","Spec.tab1","TDGDLF","WCDGDLF","JASDGDLF",
          "JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3",
          "Closures_TDGDLF","Closures_WCDGDLF","Closures_JASDGDLF",
          "Prop_TDGDLF","Prop_WCDGDLF","Prop_JASDGDLF","Prop_JASDGDLF_Zone2")
    
    
    #elasmos
    Jodie.KTCH.elasmos=combo.tbl(A=KTCH.Jodie_WCDGDLF_elasmo,B=KTCH.Jodie_JASDGDLF_elasmo,
                                 C=KTCH.Jodie_JASDGDLF_elasmo_zone,nms=NMS,SORT=SORT.s)
    Jodie.KTCH.elasmos_agg=combo.tbl_agg(A=Jodie.KTCH.elasmos)
    
    Jodie.KTCH.elasmos[,3:11]=Jodie.KTCH.elasmos[,3:11]/1000        #convert to tonnes
    Jodie.KTCH.elasmos_agg[,2:10]=Jodie.KTCH.elasmos_agg[,2:10]/1000
    
    Jodie.KTCH.elasmos[,12:15]=100*Jodie.KTCH.elasmos[,12:15]       #percentage
    Jodie.KTCH.elasmos_agg[,11:14]=100*Jodie.KTCH.elasmos_agg[,11:14]
    
    Jodie.KTCH.elasmos[,3:15]=round(Jodie.KTCH.elasmos[,3:15],1)  #round
    Jodie.KTCH.elasmos_agg[,2:14]=round(Jodie.KTCH.elasmos_agg[,2:14],1)
    
    write.csv(Jodie.KTCH.elasmos,paste(hndls,"Jodie.KTCH.elasmos_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    write.csv(Jodie.KTCH.elasmos_agg,paste(hndls,"Jodie.KTCH.elasmos_agg_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    
    #teleost
    Jodie.KTCH.scalies=combo.tbl(A=KTCH.Jodie_WCDGDLF_scalies,B=KTCH.Jodie_JASDGDLF_scalies,
                                 C=KTCH.Jodie_JASDGDLF_scalies_zone,nms=NMS,SORT=SORT.t)
    Jodie.KTCH.scalies_agg=combo.tbl_agg(A=Jodie.KTCH.scalies)
    
    Jodie.KTCH.scalies[,3:11]=Jodie.KTCH.scalies[,3:11]/1000        #convert to tonnes
    Jodie.KTCH.scalies_agg[,2:10]=Jodie.KTCH.scalies_agg[,2:10]/1000
    
    Jodie.KTCH.scalies[,12:15]=100*Jodie.KTCH.scalies[,12:15]       #percentage
    Jodie.KTCH.scalies_agg[,11:14]=100*Jodie.KTCH.scalies_agg[,11:14]
    
    Jodie.KTCH.scalies[,3:15]=round(Jodie.KTCH.scalies[,3:15],1)  #round
    Jodie.KTCH.scalies_agg[,2:14]=round(Jodie.KTCH.scalies_agg[,2:14],1)
    
    write.csv(Jodie.KTCH.scalies,paste(hndls,"Jodie.KTCH.scalies_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    write.csv(Jodie.KTCH.scalies_agg,paste(hndls,"Jodie.KTCH.scalies_agg_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
  }

  #Specific Request Block 33240 years 2007-2014
  BLOCK33240_closures=c(333240, 333241, 334240, 334241,334242, 335240, 335241)
  
  fn.extract.Jodie.catch_BLOCK=function(BLK,blk,YRS,SP)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP) & BLOCKX==BLK & Fishing_yr%in%YRS)
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    dat.blk=subset(Dat,block10%in%blk )
    dat.fishery=Dat
    Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr,dat.blk,sum)
    Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr,dat.fishery,sum)
    Agg.dat=merge(Aggregated.fishery,Aggregated.blk,by="Fishing_yr")
    colnames(Agg.dat)[2:3]=c("Total","Closure")
    Agg.dat$Percent=100* Agg.dat$Closure/ Agg.dat$Total
    Agg.dat[2:3]=Agg.dat[2:3]/1000  #convert to tonnes
    Agg.dat[2:4]=round(Agg.dat[2:4],1)
    
    Dat$UNIKS=with(Dat,paste(Fishing_yr,block10,VESSEL))
    dat.Ves=Dat[!duplicated(Dat$UNIKS),]
    
    dat.Ves$VESSEL=as.character(dat.Ves$VESSEL)
    dat.Ves$block10=as.factor(dat.Ves$block10)
    BB=vector('list',length(YRS))
    for(p in 1:length(YRS))
    {
      x=subset(dat.Ves,Fishing_yr==YRS[p])
      x=with(x,table(block10,VESSEL))
      B=data.frame(Fishing_yr=YRS[p])
      BB[[p]]=cbind(B,as.data.frame(t(rowSums(x))))
      
    }
    N.vessels=do.call(rbind,BB)
    
    return(list(Agg.dat=Agg.dat,N.vessels=N.vessels))
  }
  
  fn.extract.Jodie.effort_BLOCK=function(BLK,blk,YRS)
  {
    Dat=subset(Jodie.Effort,block10%in%BLK & Fishing_yr%in%YRS)
    dat.blk=subset(Dat,block10%in%blk)
    dat.fishery=Dat
    Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.blk,sum)
    Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    
    Agg.dat=merge(Aggregated.fishery_km.gn.hr,Aggregated.blk_km.gn.hr,by="Fishing_yr")
    colnames(Agg.dat)[2:3]=c("Total","Closure")
    Agg.dat$Percent=100* Agg.dat$Closure/ Agg.dat$Total
    Agg.dat[2:4]=round(Agg.dat[2:4],1)
    return(Agg.dat)
  }
  
  
  #Extracth catch
  KTCH_33240_elasmo=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.elasmo)$Agg.dat
  KTCH_33240_scalies=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.scalies)$Agg.dat
  Vessels_33240=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.elasmo)$N.vessels
  
  #Extract effort
  blk10in33240=subset(Data.daily,BLOCKX==33240)
  blk10in33240=unique(blk10in33240$block10)
  blk10in33240=blk10in33240[!is.na(blk10in33240)]
  EFFORT_33240=fn.extract.Jodie.effort_BLOCK(BLK=blk10in33240,blk=BLOCK33240_closures,YRS=2007:2014)

  write.csv(EFFORT_33240,paste(hndls,"EFFORT_blk_33240_2007_14.csv",sep=""),row.names=F)  
  write.csv(KTCH_33240_elasmo,paste(hndls,"KTCH_elasmo_blk_33240_2007_14.csv",sep=""),row.names=F)
  write.csv(KTCH_33240_scalies,paste(hndls,"KTCH_scalies_blk_33240_2007_14.csv",sep=""),row.names=F)
  write.csv(Vessels_33240,paste(hndls,"Vessels_blk_33240_2007_14.csv",sep=""),row.names=F)
  
  
  #27 June 2016 Request
  ASL_block_27.6.request=c(33240,32250,33250,32260,33260,32270,33270,32280,33280)
  ASL_block10_27.6.request=c(322255,322260,321260)
  Yrs.Jodie_27.6.request=2007:2014
  
  fn.extract.Jodie.catch.Req.27.6=function(blk,YRS,SP,FSHRY,what.blk,Zn)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP))
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    dat.zn=subset(Dat,zone==Zn & Fishing_yr%in%YRS)
    
    if(what.blk=='block')dat.blk=subset(Dat,BLOCKX%in%blk & Fishing_yr%in%YRS)
    if(what.blk=='block10')dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    
    if(what.blk=='block') Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+BLOCKX,dat.blk,sum)  
    if(what.blk=='block10') Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+block10,dat.blk,sum)
    
    Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr,dat.fishery,sum)
    Aggregated.zone=aggregate(LIVEWT.c~Fishing_yr,dat.zn,sum)
    
    return(list(BLKs=Aggregated.blk,Fishery=Aggregated.fishery,Zone=Aggregated.zone))
  }
  fn.extract.Jodie.effort.Req.27.6=function(blk,FSHRY,YRS,what.blk,Zn)
  {
    Dat=Jodie.Effort
    Dat$BLOCKX=with(Dat, as.numeric(paste(substr(block10,1,2),substr(block10,4,5),0,sep="")))
    
    if(what.blk=='block')dat.blk=subset(Dat,BLOCKX%in%blk & Fishing_yr%in%YRS)
    if(what.blk=='block10')dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    dat.zn=subset(Dat,zone==Zn & Fishing_yr%in%YRS)
    
    if(what.blk=='block')Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+BLOCKX,dat.blk,sum)
    if(what.blk=='block10')Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+block10,dat.blk,sum)
    
    Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    Aggregated.zone_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.zn,sum)
    
    return(list(BLKs=Aggregated.blk_km.gn.hr,Fishery=Aggregated.fishery_km.gn.hr,Zone=Aggregated.zone_km.gn.hr))
  }
  
      #Elasmos catch
  Elasmo.ktch.block_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block_27.6.request,
          FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.elasmo,what.blk='block',Zn="Zone2")
  Elasmo.ktch.block10_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block10_27.6.request,
         FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.elasmo,what.blk='block10',Zn="Zone2")
  
      #Teleosts catch
  Teleost.ktch.block_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block_27.6.request,
          FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.scalies,what.blk='block',Zn="Zone2")
  Teleost.ktch.block10_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block10_27.6.request,
          FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.scalies,what.blk='block10',Zn="Zone2")
  
      #Effort
  Effort.block_27.6.request=fn.extract.Jodie.effort.Req.27.6(blk=ASL_block_27.6.request,
          FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,what.blk='block',Zn="Zone2")
  Effort.block10_27.6.request=fn.extract.Jodie.effort.Req.27.6(blk=ASL_block10_27.6.request,
          FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,what.blk='block10',Zn="Zone2")
  
      #Export
        #Attachment 1 request
  write.csv(Elasmo.ktch.block_27.6.request$BLKs,paste(hndls,"Elasmos.ktch_block_Attachment_1.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block_27.6.request$BLKs,paste(hndls,"Teleosts.ktch_block_Attachment_1.csv",sep=""),row.names=F)  
  write.csv(Effort.block_27.6.request$BLKs,paste(hndls,"Effort.block_Attachment_1.csv",sep=""),row.names=F)  
  
        #Attachment 2 request
  write.csv(Elasmo.ktch.block10_27.6.request$BLKs,paste(hndls,"Elasmos.ktch_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Elasmo.ktch.block10_27.6.request$Fishery,paste(hndls,"Elasmos.ktch_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Elasmo.ktch.block10_27.6.request$Zone,paste(hndls,"Elasmos.ktch_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
  write.csv(Teleost.ktch.block10_27.6.request$BLKs,paste(hndls,"Teleosts.ktch_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block10_27.6.request$Fishery,paste(hndls,"Teleosts.ktch_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block10_27.6.request$Zone,paste(hndls,"Teleosts.ktch_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
  write.csv(Effort.block10_27.6.request$BLKs,paste(hndls,"Effort_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Effort.block10_27.6.request$Fishery,paste(hndls,"Effort_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Effort.block10_27.6.request$Zone,paste(hndls,"Effort_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
}

  #G 4.8.2 Demersal scalefish 2015/16
if(do.Jodies.scalies=="YES")
{

  a=subset(Data.daily,YEAR>=2015 & SPECIES%in%Scalefish.species,select=c(FisheryCode,YEAR,MONTH,METHOD,VESSEL,CONDITN,SPECIES,SNAME,LIVEWT.c))
  b=subset(Data.daily.incomplete,year>=2015 & species%in%Scalefish.species,select=c(fishery,year,month,method,vessel,conditn,species,sname1,livewt))
  names(a)=names(b)  
  a=subset(a,year== 2015 & month ==6)
  a=rbind(a,b)
  write.csv(a,paste(hndls,"Jodie.Demersal.scalies_15_16.csv",sep=""),row.names=F)
  rm(a,b)
}


#G 4.9 FishCUBE                 
if(Extract.data.FishCUBE=="YES")
{
  hndl="C:/Matias/Analyses/Catch and effort/Data_Resquests"
  
  #add variables
    #Daily
  FishCUBE.daily$ExternalDataSourceName="Shark Daily Logbook"
  FishCUBE.daily$DailyorMonthly="D"
  FishCUBE.daily$LogBookPageNumber=sapply(strsplit(FishCUBE.daily$Same.return.SNo, " "), "[", 2)
  FishCUBE.daily$FishingSeason=NA
  FishCUBE.daily$YEAR=FishCUBE.daily$YEAR.c
  FishCUBE.daily$VesselRegistration=FishCUBE.daily$VESSEL
  FishCUBE.daily$GPSLatitude=FishCUBE.daily$LAT
  FishCUBE.daily$GPSLongitude=FishCUBE.daily$LONG
  FishCUBE.daily$Block10by10=FishCUBE.daily$block10
  FishCUBE.daily$Block60by60=FishCUBE.daily$blockxFC
  FishCUBE.daily$FishingMethod=FishCUBE.daily$METHOD
  FishCUBE.daily$LiveWeight=FishCUBE.daily$LIVEWT.c
  FishCUBE.daily$LiveWeight=with(FishCUBE.daily,ifelse(is.na(LiveWeight),LIVEWT.orgnl,LiveWeight))
  FishCUBE.daily$FisheryZone=with(FishCUBE.daily,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
  FishCUBE.daily$MONTH1=with(FishCUBE.daily,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
  FishCUBE.daily$day1=with(FishCUBE.daily,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
  FishCUBE.daily$date=with(FishCUBE.daily,paste(YEAR,"-",MONTH1,"-",day1,sep=""))
  FishCUBE.daily$fishery=FishCUBE.daily$FisheryCode
  
  these.FishCUBE=c("ExternalDataSourceName","DailyorMonthly","LogBookPageNumber","FishingSeason",
                   "FINYEAR","YEAR","MONTH","date","FDAYS","LatFC","LongFC",
                   "VesselRegistration","fishery","FisheryZone","Landing.Port","Block10by10","Block","Block60by60","Bioregion",
                   "FishingMethod","SPECIES","SNAME","RSCommonName","LiveWeight","licence","BDAYS","rowid")
  
  FishCUBE=FishCUBE.daily[,match(c(these.FishCUBE,"RSSpeciesId"),names(FishCUBE.daily))]
  names(FishCUBE)[match(c("Landing.Port","LatFC","LongFC"),names(FishCUBE))]=c("LandingPort","Lat","Long")
  
  #add removed May and June 2006
  Data.daily.FC.2005_06$ExternalDataSourceName="Shark Daily Logbook"
  Data.daily.FC.2005_06$DailyorMonthly="D"
  Data.daily.FC.2005_06$LogBookPageNumber=Data.daily.FC.2005_06$DSNo
  Data.daily.FC.2005_06$FishingSeason=NA
  Data.daily.FC.2005_06$YEAR=Data.daily.FC.2005_06$year
  Data.daily.FC.2005_06$VesselRegistration=Data.daily.FC.2005_06$vessel
  Data.daily.FC.2005_06$GPSLatitude=Data.daily.FC.2005_06$LatDeg
  Data.daily.FC.2005_06$GPSLongitude=Data.daily.FC.2005_06$LongDeg
  Data.daily.FC.2005_06$Block10by10=Data.daily.FC.2005_06$block10
  Data.daily.FC.2005_06$Block60by60=Data.daily.FC.2005_06$blockxFC
  Data.daily.FC.2005_06$FishingMethod=Data.daily.FC.2005_06$method
  Data.daily.FC.2005_06$LiveWeight=Data.daily.FC.2005_06$livewt
  Data.daily.FC.2005_06$FisheryZone=with(Data.daily.FC.2005_06,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
  Data.daily.FC.2005_06$MONTH=Data.daily.FC.2005_06$month
  Data.daily.FC.2005_06$MONTH1=with(Data.daily.FC.2005_06,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
  Data.daily.FC.2005_06$day1=with(Data.daily.FC.2005_06,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
  Data.daily.FC.2005_06$FINYEAR=Data.daily.FC.2005_06$finyear
  Data.daily.FC.2005_06$FDAYS=Data.daily.FC.2005_06$fdays
  Data.daily.FC.2005_06$Lat=Data.daily.FC.2005_06$LatFC  
  Data.daily.FC.2005_06$Long=Data.daily.FC.2005_06$LongFC
  Data.daily.FC.2005_06$LandingPort=Data.daily.FC.2005_06$port
  Data.daily.FC.2005_06$SPECIES=Data.daily.FC.2005_06$species   
  Data.daily.FC.2005_06$SNAME=Data.daily.FC.2005_06$sname1
  Data.daily.FC.2005_06$BDAYS=Data.daily.FC.2005_06$bdays
  FishCUBE.data.daily.2005_06=Data.daily.FC.2005_06[,match(names(FishCUBE),names(Data.daily.FC.2005_06))]
  FishCUBE.data.daily.2005_06$date=as.character(FishCUBE.data.daily.2005_06$date)
  
  COLss=names(FishCUBE)
  for(cl in 1:length(COLss))
  {
    if(is.factor(FishCUBE[,cl])) FishCUBE[,cl]=as.character(FishCUBE[,cl])
    if(is.factor(FishCUBE.data.daily.2005_06[,cl])) FishCUBE.data.daily.2005_06[,cl]=as.character(FishCUBE.data.daily.2005_06[,cl])
  }
  FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.2005_06)
  
  
  #add latest incomplete year to daily
   if(nrow(Data.daily.incomplete)>0)
  {
    Data.daily.incomplete$ExternalDataSourceName="Shark Daily Logbook"
    Data.daily.incomplete$DailyorMonthly="D"
    Data.daily.incomplete$LogBookPageNumber=Data.daily.incomplete$DSNo
    Data.daily.incomplete$FishingSeason=NA
    Data.daily.incomplete$YEAR=Data.daily.incomplete$year
    Data.daily.incomplete$VesselRegistration=Data.daily.incomplete$vessel
    Data.daily.incomplete$GPSLatitude=Data.daily.incomplete$LatDeg
    Data.daily.incomplete$GPSLongitude=Data.daily.incomplete$LongDeg
    Data.daily.incomplete$Block10by10=Data.daily.incomplete$block10
    Data.daily.incomplete$Block60by60=Data.daily.incomplete$blockxFC
    Data.daily.incomplete$FishingMethod=Data.daily.incomplete$method
    Data.daily.incomplete$LiveWeight=Data.daily.incomplete$livewt
    Data.daily.incomplete$FisheryZone=with(Data.daily.incomplete,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
    Data.daily.incomplete$MONTH=Data.daily.incomplete$month
    Data.daily.incomplete$MONTH1=with(Data.daily.incomplete,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
    Data.daily.incomplete$day1=with(Data.daily.incomplete,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
    Data.daily.incomplete$FINYEAR=Data.daily.incomplete$finyear
    Data.daily.incomplete$FDAYS=Data.daily.incomplete$fdays
    Data.daily.incomplete$Lat=Data.daily.incomplete$LatFC  
    Data.daily.incomplete$Long=Data.daily.incomplete$LongFC
    Data.daily.incomplete$LandingPort=Data.daily.incomplete$port
    Data.daily.incomplete$SPECIES=Data.daily.incomplete$species   
    Data.daily.incomplete$SNAME=Data.daily.incomplete$sname1
    Data.daily.incomplete$BDAYS=Data.daily.incomplete$bdays
    #Data.daily.incomplete$rowid=NA
     
    FishCUBE.data.daily.incomplete=Data.daily.incomplete[,match(names(FishCUBE),names(Data.daily.incomplete))]
    FishCUBE.data.daily.incomplete$date=as.character(FishCUBE.data.daily.incomplete$date)
     
    for(cl in 1:length(COLss))if(is.factor(FishCUBE.data.daily.incomplete[,cl])) FishCUBE.data.daily.incomplete[,cl]=as.character(FishCUBE.data.daily.incomplete[,cl])
    FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.incomplete)
  }
  FishCUBE=subset(FishCUBE,!(SPECIES==99999 | is.na(SPECIES)))
  
  #add 0 catch
  if(nrow(Data.daily.FC.NA.sp)>0)
  {
    Data.daily.FC.NA.sp$ExternalDataSourceName="Shark Daily Logbook"
    Data.daily.FC.NA.sp$DailyorMonthly="D"
    Data.daily.FC.NA.sp$LogBookPageNumber=Data.daily.FC.NA.sp$DSNo
    Data.daily.FC.NA.sp$FishingSeason=NA
    Data.daily.FC.NA.sp$YEAR=Data.daily.FC.NA.sp$year
    Data.daily.FC.NA.sp$VesselRegistration=Data.daily.FC.NA.sp$vessel
    Data.daily.FC.NA.sp$GPSLatitude=Data.daily.FC.NA.sp$LatDeg
    Data.daily.FC.NA.sp$GPSLongitude=Data.daily.FC.NA.sp$LongDeg
    Data.daily.FC.NA.sp$Block10by10=Data.daily.FC.NA.sp$block10
    Data.daily.FC.NA.sp$Block60by60=Data.daily.FC.NA.sp$blockxFC
    Data.daily.FC.NA.sp$FishingMethod=Data.daily.FC.NA.sp$method
    Data.daily.FC.NA.sp$LiveWeight=Data.daily.FC.NA.sp$livewt
    Data.daily.FC.NA.sp$FisheryZone=with(Data.daily.FC.NA.sp,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
    Data.daily.FC.NA.sp$MONTH=Data.daily.FC.NA.sp$month
    Data.daily.FC.NA.sp$MONTH1=with(Data.daily.FC.NA.sp,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
    Data.daily.FC.NA.sp$day1=with(Data.daily.FC.NA.sp,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
    Data.daily.FC.NA.sp$FINYEAR=Data.daily.FC.NA.sp$finyear
    Data.daily.FC.NA.sp$FDAYS=Data.daily.FC.NA.sp$fdays
    Data.daily.FC.NA.sp$Lat=Data.daily.FC.NA.sp$LatFC  
    Data.daily.FC.NA.sp$Long=Data.daily.FC.NA.sp$LongFC
    Data.daily.FC.NA.sp$LandingPort=Data.daily.FC.NA.sp$port
    Data.daily.FC.NA.sp$SPECIES=Data.daily.FC.NA.sp$species   
    Data.daily.FC.NA.sp$SNAME=Data.daily.FC.NA.sp$sname1
    Data.daily.FC.NA.sp$BDAYS=Data.daily.FC.NA.sp$bdays
    #Data.daily.FC.NA.sp$rowid=NA
    
    FishCUBE.data.daily.NA.sp=Data.daily.FC.NA.sp[,match(names(FishCUBE),names(Data.daily.FC.NA.sp))]
    FishCUBE.data.daily.NA.sp$date=as.character(FishCUBE.data.daily.NA.sp$date)
    
    for(cl in 1:length(COLss))if(is.factor(FishCUBE.data.daily.NA.sp[,cl])) FishCUBE.data.daily.NA.sp[,cl]=as.character(FishCUBE.data.daily.NA.sp[,cl])
    FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.NA.sp)
  }
  
  
    #Monthly
  FishCUBE.monthly=subset(FishCUBE.monthly,TYPE.DATA=="monthly")    #keep only monthly records
  a=subset(Data.monthly.original,Same.return%in%unique(FishCUBE.monthly$Same.return),select=c(Same.return,FDAYS))
  a=a[!duplicated(a$Same.return),]
  FishCUBE.monthly=merge(FishCUBE.monthly,a,by="Same.return",all.x=T)
  FishCUBE.monthly$FishingMethod=FishCUBE.monthly$Fishing.method
  FishCUBE.monthly$LandingPort=FishCUBE.monthly$Landing.Port
  FishCUBE.monthly$LatFC=FishCUBE.monthly$LAT
  FishCUBE.monthly$LongFC=FishCUBE.monthly$LONG
  FishCUBE.monthly$YEAR=FishCUBE.monthly$YEAR.c
  FishCUBE.monthly$Block=FishCUBE.monthly$blockxFC
  FishCUBE.monthly$Block60by60=FishCUBE.monthly$blockxFC
  FishCUBE.monthly$FishingMethod=FishCUBE.monthly$METHOD
  FishCUBE.monthly$LiveWeight=FishCUBE.monthly$LIVEWT.c
  FishCUBE.monthly$VesselRegistration=FishCUBE.monthly$VESSEL
  FishCUBE.monthly$fishery=FishCUBE.monthly$FisheryCode
  
  FishCUBE.monthly$ExternalDataSourceName="FLAMSCAE"
  FishCUBE.monthly$DailyorMonthly="M"
  
  FishCUBE.monthly$LogBookPageNumber=NA
  FishCUBE.monthly$FishingSeason=NA
  FishCUBE.monthly$day=NA
  FishCUBE.monthly$date=NA
  FishCUBE.monthly$GPSLatitude=NA
  FishCUBE.monthly$GPSLongitude=NA
  FishCUBE.monthly$Block10by10=NA

  
  FishCUBE.monthly=FishCUBE.monthly[,match(these.FishCUBE,names(FishCUBE.monthly))]
  names(FishCUBE.monthly)[match(c("Landing.Port","LatFC","LongFC"),names(FishCUBE.monthly))]=c("LandingPort","Lat","Long")
  
  
  #Fix RSSpeciesId for reapportioned sharks
    #Daily
  FishCUBE$RSSpeciesId=with(FishCUBE,ifelse(SPECIES==17001,63,
              ifelse(SPECIES==18003,70,
              ifelse(SPECIES==17003,64,
              ifelse(SPECIES==18007,74,
              ifelse(SPECIES==22999,943,RSSpeciesId))))))
  
  FishCUBE$RSCommonName=with(FishCUBE,ifelse(SPECIES==17001,"Gummy Shark",
              ifelse(SPECIES==18003,"Dusky Whaler",
              ifelse(SPECIES==17003,"Whiskery Shark",
              ifelse(SPECIES==18007,"Sandbar Shark",
                     ifelse(SPECIES==22999,"Sharks",RSCommonName))))))
  
  FishCUBE$SNAME=with(FishCUBE,ifelse(SPECIES==17001,"SHARK, GUMMY",
              ifelse(SPECIES==18003,"SHARK, BRONZE WHALER",
              ifelse(SPECIES==17003,"SHARK, WHISKERY",
              ifelse(SPECIES==18007,"SHARK, THICKSKIN (SANDBAR)",
              ifelse(SPECIES==22999,"SHARK, OTHER",SNAME))))))
  
  
  #Monthly
  # FishCUBE.monthly$RSSpeciesId=with(FishCUBE.monthly,ifelse(SPECIES==17001,63,
  #              ifelse(SPECIES==18003,70,
  #              ifelse(SPECIES==17003,64,
  #              ifelse(SPECIES==18007,74,
  #              ifelse(SPECIES==22999,943,RSSpeciesId))))))
   
  FishCUBE.monthly$RSCommonName=with(FishCUBE.monthly,ifelse(SPECIES==17001,"Gummy Shark",
               ifelse(SPECIES==18003,"Dusky Whaler",
               ifelse(SPECIES==17003,"Whiskery Shark",
               ifelse(SPECIES==18007,"Sandbar Shark",
               ifelse(SPECIES==22999,"Sharks",RSCommonName))))))

 FishCUBE.monthly$SNAME=with(FishCUBE.monthly,ifelse(SPECIES==17001,"SHARK, GUMMY",
               ifelse(SPECIES==18003,"SHARK, BRONZE WHALER",
               ifelse(SPECIES==17003,"SHARK, WHISKERY",
               ifelse(SPECIES==18007,"SHARK, THICKSKIN (SANDBAR)",
               ifelse(SPECIES==22999,"SHARK, OTHER",SNAME))))))
  
  
  #change NA to blank
    #Daily
  FishCUBE$Lat=as.character(FishCUBE$Lat)
  FishCUBE$Long=as.character(FishCUBE$Long)
  FishCUBE$Lat=with(FishCUBE,ifelse(is.na(Lat),"",Lat))
  FishCUBE$Long=with(FishCUBE,ifelse(is.na(Long),"",Long))
  
    #Monthly
  FishCUBE.monthly$Lat=as.character(FishCUBE.monthly$Lat)
  FishCUBE.monthly$Long=as.character(FishCUBE.monthly$Long)
  FishCUBE.monthly$Lat=with(FishCUBE.monthly,ifelse(is.na(Lat),"",Lat))
  FishCUBE.monthly$Long=with(FishCUBE.monthly,ifelse(is.na(Long),"",Long))
  
  
  #export data
  
    #monthly
  hndl.FishCUBE="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/CAES Monthly Fisheries/02. Data/"
  system.time(fwrite(FishCUBE.monthly,paste(hndl.FishCUBE,"2. CAES Monthly Data - Shark - Corrected with R.csv",sep=""),row.names=F))
  

    #daily
  hndl.FishCUBE.daily="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/Shark/02. Data/"
  system.time(fwrite(FishCUBE,paste(hndl.FishCUBE.daily,"2. Shark Daily Logbook Data - Corrected with R.csv",sep=""),row.names=F))
  #system.time(write.csv(FishCUBE,paste(hndl.FishCUBE,"Shark Data Extract for FishCubeWA.",Sys.Date(),".csv",sep=""),row.names=F))

 # library("openxlsx")
 # options(java.parameters = "-Xmx8000m")
#  system.time(write.xlsx(FishCUBE,paste(hndl.FishCUBE,"Shark Data Extract for FishCubeWA.",Sys.Date(),".xlsx",sep=""), colNames = TRUE))
  
  #notify Vero by email
  # function.send.email(
  #   to ="Paul.Fildes@dpird.wa.gov.au",
  #   subject ="Annual data",
  #   body =paste("Data corrections done and uploaded.",
  #             "Data files stored in",hndl.FishCUBE,"and",
  #             hndl.FishCUBE.daily),                     
  #   Attachment=NULL
  # )
  
}

#G 4.10 Steve's map for interview
if (do.Steves=="YES")
{
  a=South.WA.long[1]:South.WA.long[2]
  b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
  
  GRID="Y"
  if(GRID=="Y")tiff(file="C:/Matias/Analyses/Catch and effort/Outputs/for Steve/Steve,map.grid.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  if(GRID=="N")tiff(file="C:/Matias/Analyses/Catch and effort/Outputs/for Steve/Steve,map.NoGrid.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  
  plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
  if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                               nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
  if(GRID=="Y")
  {
    axis(side = 1, at =LONGG, labels = F, tcl = 30.6)
    axis(side = 4, at = LATT, labels = F,tcl =34)
    par(new=T)
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    
  }
  axis(side = 1, at =LONGG, labels = LONGG, tcl = -.5)
  axis(side = 2, at = LATT, labels = LATT,tcl =-.5,las=1)
  
  text(116,-28.78,("Geraldton"),col="black", cex=1.1)
  points(114.83,-28.86,pch=19)
  text(116.73,Perth[2],("Perth"),col="black", cex=1.1)
  points(115.86,-31.95,pch=19)
  text(116.73,-33.55,("Bunbury"),col="black", cex=1.1)
  points(115.6,-33.55,pch=19)
  text(115.95,-34.085,("Augusta"),col="black", cex=1.1)
  points(115.15,-34.31,pch=19)
  text(117.7,-34.8,("Albany"),col="black", cex=1.1)
  points(117.8,-35,pch=19)
  text(122,-33.66,("Esperance"),col="black", cex=1.1)
  points(121.9,-33.86,pch=19)
  
  mtext("Latitude (?S)",side=2,line=2.75,las=3,cex=1.3)
  mtext("Longitude (?E)",side=1,line=2.75,cex=1.3)
  
  
  dev.off()
  
}

#G 4.11 EXTRA FIGURES FOR AMM 2013
if(do.AMM.2013=="YES")
{
  setwd(paste(hndl,"/AMM_2013",sep=""))
  if(Use.Previos.Sofar=="YES")
  {
    fn.fig.AMM=function(GROUP,LAT1,LAT2,INT,INT2,THIS)
    {
      dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
      
      annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
      annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
      
      wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
      
      #add previous years
      Prev.zn=Spec.catch.zone.pre.2013[,match(THIS,names(Spec.catch.zone.pre.2013))]
      names(Prev.zn)=names(wide)
      Prev.zn[,2:4]=Prev.zn[,2:4]*1000
      wide=rbind(Prev.zn,wide)
      
      #total 
      annual.catch.total=data.frame(FINYEAR=wide$FINYEAR,Total.catch=rowSums(wide[,2:4]))
      
      names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
      names(wide)[match("FINYEAR",names(wide))]="finyear"
      
      #plot
      fun.fig.SoFar(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
      
      return(Tot.Catch=annual.catch.total)
    }
  }
  
  Lista.especies=list(whiskery=c("Whiskery.catch.WC","Whiskery.catch.Zn1","Whiskery.catch.Zn2"),
                      gummy=c("Gummy.catch.WC","Gummy.catch.Zn1","Gummy.catch.Zn2"),
                      dusky=c("Dusky.catch.WC","Dusky.catch.Zn1","Dusky.catch.Zn2"),
                      sandbar=c("Sandbar.catch.WC","Sandbar.catch.Zn1","Sandbar.catch.Zn2"))
  SPec=list(17003,17001,c(18001,18003),18007)
  sp=c("whiskery","gummy","dusky","sandbar")
  Catch.range=list(c(175,225),c(350,450),c(200,300),c(120,120))
  
  STORE.Tot.Catch=vector('list',length=length(TARGETS))
  names(STORE.Tot.Catch)=sp
  for (i in 1:length(TARGETS))
  {
    jpeg(file=paste("Catch.zone.",sp[i],".jpeg",sep=""),width = 2400, height = 2400,units = "px", res = 300)
    STORE.Tot.Catch[[i]]=fn.fig.AMM(SPec[[i]],-27,-40,100,50,c("Fin.yr",Lista.especies[[i]]))
    abline(Catch.range[[i]][1],0,col=2)
    abline(Catch.range[[i]][2],0,col=2)
    dev.off()  
  }
  
}

#G 4.12 Catch and effort Perth metro closures
#Closure
CLOSED=c("31120","31130","31140","31150","32120","32130","32140","32150")
if(do.Perth.metro.closure=="YES")
{
  #1. Aggregate effort by year and blockfor WCDGDLF
  
  #Daily
  a=subset(Effort.daily,zone=="West")
  #km gn days 
  if(Use.Date=="NO")  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~ID+vessel+blockx+finyear,data=a,max,na.rm=T)
  if(Use.Date=="YES") Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~date+vessel+blockx+finyear,data=a,max,na.rm=T) 
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~blockx+finyear,data=Effort.Perth.metro,sum,na.rm=T)
  
  #km gn hours
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~finyear+vessel+ID+blockx,data=a,max,na.rm=T)
  #if(Use.Date=="YES")Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~blockx+finyear,data=Effort.Perth.metro.hrs,sum,na.rm=T)
  Effort.Perth.metro.daily=merge(Effort.Perth.metro.hrs,Effort.Perth.metro,by=c("blockx","finyear"))
  
  #Monthly
  a=subset(Effort.monthly,zone=="West" & !FINYEAR%in%Daily.l.years)
  #km gn days             
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+BLOCKX+FINYEAR,data=a,max,na.rm=T)
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~BLOCKX+FINYEAR,data=Effort.Perth.metro,sum,na.rm=T)
  
  #km gn hours
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~MONTH+BLOCKX+VESSEL+BLOCKX+FINYEAR,data=a,max,na.rm=T)
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~BLOCKX+FINYEAR,data=Effort.Perth.metro.hrs,sum,na.rm=T)
  Effort.Perth.metro.monthly=merge(Effort.Perth.metro.hrs,Effort.Perth.metro,by=c("BLOCKX","FINYEAR"))
  
  Effort.Perth.metro.daily$finyear=as.character(Effort.Perth.metro.daily$finyear)
  Effort.Perth.metro.monthly$FINYEAR=as.character(Effort.Perth.metro.monthly$FINYEAR)
  names(Effort.Perth.metro.monthly)=names(Effort.Perth.metro.daily)
  Effort.Perth.metro=rbind(Effort.Perth.metro.monthly,Effort.Perth.metro.daily)
  
  
  
  #2. Aggregate catch by year and blockfor WCDGDLF
  KTCH.Perth.metro=aggregate(LIVEWT.c~BLOCKX+FINYEAR+SPECIES+SNAME,data=subset(Data.monthly,zone=="West"),sum,na.rm=T)
  
  #3. Extract effort Metro closure
  Effort.Perth.metro.closure=subset(Effort.Perth.metro,blockx%in%as.numeric(CLOSED))
    
  #4. Sharks and Demersal scalefish suite component
  KTCH.Perth.metro.shrk=subset(KTCH.Perth.metro,SPECIES%in%Elasmo.species)
  KTCH.Perth.metro.scale=subset(KTCH.Perth.metro,SPECIES%in%Suite)
  
  #5. Catch from Metro Closure
  KTCH.Perth.metro.closure.shrk=subset(KTCH.Perth.metro,SPECIES%in%Elasmo.species & BLOCKX%in%CLOSED)
  KTCH.Perth.metro.closure.scale=subset(KTCH.Perth.metro,SPECIES%in%Suite  & BLOCKX%in%CLOSED)
  
  #6. Aggregate by year
  Eff=aggregate(Km.Gillnet.Days.c~finyear,Effort.Perth.metro,sum)
  Eff.closure=aggregate(Km.Gillnet.Days.c~finyear,Effort.Perth.metro.closure,sum)
  
  Agg.shark=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.shrk,sum)
  Agg.scale=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.scale,sum)
  
  Agg.shark.closed=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.closure.shrk,sum)
  Agg.scale.closed=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.closure.scale,sum)
  
  
  #Export
  fn.plt=function(a,b,d,e)
  {
    if(length(e)<length(d))
    {
      id=which(!d%in%e)
      test=data.frame(e,b)
      test1=data.frame(d[id],0)
      names(test1)=names(test)
      b1=rbind(test,test1)
      b=b1[,2]
    }
    plot(1:length(a),a,ylab="",xlab="",type='l',lwd=2,ylim=c(0,max(a)),xaxt='n',cex.axis=1.35)    
    lines(1:length(b),b,lwd=2,lty=3)
    axis(1,1:length(a),F,tck=-0.02)
    if(length(a)>20)axis(1,seq(1,length(a),10),d[seq(1,length(a),10)],tck=-0.04,cex.axis=1.5)
    if(length(a)<20)axis(1,seq(1,length(a),2),d[seq(1,length(a),2)],tck=-0.04,cex.axis=1.5)
  }
  tiff(file=paste(hndl,"/Metro_closure_catch_effort.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.5,.6,.1,.1),las=1,mgp=c(2,0.65,0))  
  
  #Effort
  fn.plt(Eff$Km.Gillnet.Days.c/1000,Eff.closure$Km.Gillnet.Days.c/1000,Eff$finyear,Eff.closure$finyear)
  mtext("Effort (1000 km gn d)",2,line=2.7,font=1,las=0,cex=1.45)
  legend("topright",c("WCDGDLF","Metro closure"),bty='n',lty=c(1,3),col=1,lwd=2,cex=2)
  
  #Shark catch
  fn.plt(Agg.shark$LIVEWT.c/1000,Agg.shark.closed$LIVEWT.c/1000,Agg.shark$FINYEAR,Agg.shark.closed$FINYEAR)
  mtext("Catch (tonnes)",2,line=2.7,font=1,las=0,cex=1.5)
  legend("topright",c("Shark"),bty='n',cex=2.5)
  
  #Scalefish suite catch
  fn.plt(Agg.scale$LIVEWT.c/1000,Agg.scale.closed$LIVEWT.c/1000,Agg.scale$FINYEAR,Agg.scale.closed$FINYEAR)
  mtext("Catch (tonnes)",2,line=2.7,font=1,las=0,cex=1.5)
  legend("topright",c("Demersal scalefish suite"),bty='n',cex=2.5)
  
  mtext("Financial year",side=1,line=-1.5,font=1,las=0,cex=1.5,outer=T) 
  dev.off()
  
   
  #write.csv(Effort.Perth.metro,paste(hndl,"/Effort.Perth.metro.csv",sep=""),row.names=F)
  #write.csv(KTCH.Perth.metro,paste(hndl,"/Catch.Perth.metro.csv",sep=""),row.names=F)
}

#G 4.13  Heather's request proportion catch by LL and GN
if(do.Heather.request=="YES")
{
  Fn.prop.ktch.ll.gn=function(datGN,datLL)
  {
    no.scalies=paste(1975:1997,substr(1976:1998,3,4),sep="-")
    datGN=subset(datGN,!FINYEAR%in%no.scalies)   #no scalies reported before 1998-99
    datLL=subset(datLL,!FINYEAR%in%no.scalies)
    datGN$scalie=with(datGN,ifelse(SPECIES%in%Scalefish.species,"YES" ,"NO"))
    datLL$scalie=with(datLL,ifelse(SPECIES%in%Scalefish.species,"YES" ,"NO"))
    TabGN=with(datGN,table(FINYEAR, scalie))
    TabLL=with(datLL,table(FINYEAR, scalie))
    
    TabGN=as.data.frame.matrix(TabGN)
    TabLL=as.data.frame.matrix(TabLL)
    
    TabGN$Tot=TabGN$YES+TabGN$NO
    TabLL$Tot=TabLL$YES+TabLL$NO
    
    TabGN$Prop.scale=with(TabGN, YES/ Tot)
    TabLL$Prop.scale=with(TabLL, YES/ Tot)
    
    TabGN$Prop.scale=with(TabGN,ifelse(Tot<100,NA,Prop.scale))
    TabLL$Prop.scale=with(TabLL,ifelse(Tot<100,NA,Prop.scale))
    
    plot(TabGN$Prop.scale,type='l',col=2,lwd=2,ylab="Proportion of scalefish",xlab="",
         ylim=c(0,1),xaxt='n',cex.lab=1.5)
    lines(TabLL$Prop.scale,lwd=2,col=4)
    legend("topleft",c("Gillnet","Longline"),lty=1,lwd=2,col=c(2,4),bty='n',cex=1.25)
    axis(1,1:nrow(TabGN),F)
    axis(1,seq(1,nrow(TabGN),5),rownames(TabGN)[seq(1,nrow(TabGN),5)])
  }
    
  tiff(file=paste(hndl,"/Heather_prop_LL_GN.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,1),mar=c(3,4,1,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2.5,0.8,0))
  
  Fn.prop.ktch.ll.gn(subset(Data.monthly.GN, zone=="West"),
                     subset(Data.monthly.LL, zone=="West"))
  mtext("West Coast",3,cex=1.5,line=-2)
  
  Metro=Fn.prop.ktch.ll.gn(subset(Data.monthly.GN, BLOCKX%in%as.numeric(CLOSED)), 
                           subset(Data.monthly.LL, BLOCKX%in%as.numeric(CLOSED)))
  mtext("Closure",3,cex=1.5,line=-2)
  mtext("Financial year",1,outer=T,cex=1.5) 
  dev.off()  
}

#G 4.14 Annual scalefish Catch for Dave F.
if(do.Dave.F=="YES")
{
  solo.monthly=subset(FINYEAR.monthly,!FINYEAR.monthly%in%FINYEAR.daily)
  a=subset(Data.monthly,FINYEAR%in%solo.monthly & SPECIES%in%Scalefish.species & Bioregion%in%c("SC","WC"))
  if(nrow(a)>0)
  {
    KTCH.scalies_monthly=aggregate(LIVEWT.c~BLOCKX+Bioregion+FINYEAR+SPECIES+SNAME,data=a,sum,na.rm=T)
    write.csv(KTCH.scalies_monthly,paste(hndl,"/Dave F/Dave_F_scalefish_monthly.csv",sep=""),row.names=F)
    rm(a)
  }
  KTCH.scalies_daily=aggregate(LIVEWT.c~block10+Bioregion+FINYEAR+SPECIES+SNAME,data=subset(Data.daily,SPECIES%in%Scalefish.species),sum,na.rm=T)
  write.csv(KTCH.scalies_daily,paste(hndl,"/Dave F/Dave_F_scalefish_daily.csv",sep=""),row.names=F)

  #email data
  # function.send.email(
  #   to ="David.Fairclough@dpird.wa.gov.au",
  #   subject ="Annual data",
  #   body ="Please find attached the annual monthly records",                     
  #   Attachment=paste(hndl,"/Dave F/Dave_F_scalefish_monthly.csv",sep="")
  # )
  # function.send.email(
  #   to ="David.Fairclough@dpird.wa.gov.au",
  #   subject ="Annual data",
  #   body ="Please find attached the annual daily records",                     
  #   Attachment=paste(hndl,"/Dave F/Dave_F_scalefish_daily.csv",sep="")
  # )
  
  do.table.all.yrs="NO"
  if(do.table.all.yrs=="YES")
  {
    This.fish=match(c(Blue_morwong,Blue_groper,West_Australian_dhufish,Pink_snapper,Boarfishes,Samsonfish,Redfishes,Mulloway,Sweetlips,
                      Baldchin_groper),Scalefish.species)
    This.fish=This.fish[!is.na(This.fish)]
    DAT.D=subset(Data.daily,SPECIES%in%Scalefish.species&Bioregion%in%c("SC","WC"))
    DAT.D$SPECIES.agg=with(DAT.D,
        ifelse(SPECIES%in%Blue_morwong,"Blue_morwong",
        ifelse(SPECIES%in%Blue_groper,"Blue_groper",
        ifelse(SPECIES%in%West_Australian_dhufish,"West_Australian_dhufish",
        ifelse(SPECIES%in%Pink_snapper,"Pink_snapper",
        ifelse(SPECIES%in%Boarfishes,"Boarfishes",
        ifelse(SPECIES%in%Samsonfish,"Samsonfish",
        ifelse(SPECIES%in%Redfishes,"Redfishes",
        ifelse(SPECIES%in%Mulloway,"Mulloway",
        ifelse(SPECIES%in%Sweetlips,"Sweetlips",
        ifelse(SPECIES%in%Baldchin_groper,"Baldchin_groper",
        ifelse(SPECIES%in%Scalefish.species[-This.fish],"Other_scalefish",
        NA))))))))))))
  
    TABLA=aggregate(LIVEWT.c~Bioregion+FINYEAR+SPECIES.agg,data=DAT.D,sum,na.rm=T)
    Tab.suite=aggregate(LIVEWT.c~Bioregion+FINYEAR,data=subset(DAT.D,SPECIES%in%Suite),sum,na.rm=T)
    Tab.suite$SPECIES.agg="Suite"
    TABLA=rbind(TABLA,Tab.suite)
    wide <- reshape(TABLA,v.names="LIVEWT.c",timevar="Bioregion",idvar=c("FINYEAR","SPECIES.agg"),direction="wide")
    write.csv(wide,paste(hndl,"/Dave F/Daily_table_all_yrs_bioregion.csv",sep=""),row.names=F)
  
    TABLA=aggregate(LIVEWT.c~zone+FINYEAR+SPECIES.agg,data=DAT.D,sum,na.rm=T)
    Tab.suite=aggregate(LIVEWT.c~zone+FINYEAR,data=subset(DAT.D,SPECIES%in%Suite),sum,na.rm=T)
    Tab.suite$SPECIES.agg="Suite"
    TABLA=rbind(TABLA,Tab.suite)
    wide <- reshape(TABLA,v.names="LIVEWT.c",timevar="zone",idvar=c("FINYEAR","SPECIES.agg"),direction="wide")
    write.csv(wide,paste(hndl,"/Dave F/Daily_table_all_yrs_zone.csv",sep=""),row.names=F)
  
    rm(DAT.D,wide)
  }
}

#G 4.17 Annual scalefish Catch for Jeff N.
if(do.Jeff.N=="YES")
{
  a=subset(Data.monthly,SPECIES%in%Scalefish.species & Bioregion=="SC")
  
  Demersal.scale=c(224901,228002,258000,258004,258005,258006,264004,287000,288000,311000,
                   311005,311006,311078,311100,311152,311170,320000,320000,330001,346911,
                   346914,353001,361004,367000,369002,377004,384002,384904,384999,386000,
                   439002,445001)
  
  Jeff.list=list(Demersal.scale=Demersal.scale,Redfish=c(258000,258004,258006),
    Pink.snapper=353001,Blue.morwong=377004,Blue.groper=384002,Hapuku=311006) 
  
  b=subset(a,SPECIES%in%Jeff.list$Demersal.scale)
  Dummy1=aggregate(LIVEWT.c~FINYEAR,data=b,sum,na.rm=T)
  names(Dummy1)[2]="Demersal.scalefish"
    
  a$SP.jeff=with(a,
      ifelse(SPECIES%in%Jeff.list$Redfish,"Redfish",
      ifelse(SPECIES%in%Jeff.list$Pink.snapper,"Pink.snapper",
      ifelse(SPECIES%in%Jeff.list$Blue.morwong,"Blue.morwong",
      ifelse(SPECIES%in%Jeff.list$Blue.groper,"Blue.groper",
      ifelse(SPECIES%in%Jeff.list$Hapuku,"Hapuku",NA
        ))))))
  a=subset(a,!is.na(SP.jeff))
  
  Dummy2=aggregate(LIVEWT.c~FINYEAR+SP.jeff,data=a,sum,na.rm=T)
  Dummy2=reshape(Dummy2, v.names = "LIVEWT.c", idvar = "FINYEAR",
                 timevar ="SP.jeff" , direction = "wide")
  names(Dummy2)[2:ncol(Dummy2)]=substr(names(Dummy2)[2:ncol(Dummy2)],10,50)
  Dummy2[is.na(Dummy2)]=0
  Dummy=merge(Dummy1,Dummy2,by="FINYEAR")
  OuTT=paste(hndl,"/Jeff N/","Annual_ktch(kg).csv",sep="")  
  write.csv(Dummy,OuTT,row.names=F)
  
  #email data
  # function.send.email(
  #   to ="Jeffrey.Norriss@dpird.wa.gov.au",
  #   subject ="Annual data",
  #   body ="Please find attached the annual catch records",                     
  #   Attachment=OuTT
  # )
  
  rm(a)

}

#G 4.18 Annual scalefish Catch for Paul L.
if(do.Paul.L=="YES")
{
  Paul.list=list(Samson.fish=337007)  
  a=subset(Data.monthly,SPECIES%in%Scalefish.species)
  
  Dummy=aggregate(LIVEWT.c~FINYEAR+Bioregion+SPECIES,data=subset(a,SPECIES%in%Paul.list[[1]]),sum,na.rm=T)
  Dummy$SPECIES="Samson.fish"
  
  OuTT=paste(hndl,"/Paul L/","Annual_ktch(kg).csv",sep="")
  write.csv(Dummy,OuTT,row.names=F)
  
  #email data
  # function.send.email(
  #   to ="paul.lewis@dpird.wa.gov.au",
  #   subject ="Annual data",
  #   body ="Please find attached the annual daily records",                     
  #   Attachment=OuTT
  # )
  
  rm(a)
}

#G 4.19 Hammerheads listing
if(do.Hammerheads=="YES")
{
  KTCH=aggregate(LIVEWT.c/1000~FINYEAR,Hammerheads,sum)
   names(KTCH)[2]="LIVEWT.c.ton"
  MissYr=FINYrs[which(!FINYrs%in%KTCH$FINYEAR)]
  KTCH=rbind(KTCH,data.frame(FINYEAR=MissYr,LIVEWT.c.ton=NA))
  KTCH=KTCH[order(KTCH$FINYEAR),]
  Ny=length(FINYrs)
  tiff(file=paste(hndl,"/Annual.Hammerhead.Catch.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mgp=c(2.3,.7,0),las=1)
  plot(1:Ny,KTCH$LIVEWT.c.ton,xlab="Financial year",ylab="Total catch (tonnes)",xaxt='n',
       type="b",pch=19,cex=2,cex.lab=2,cex.axis=1.25)
  axis(1,1:Ny,F,tck=-0.01)
  axis(1,seq(1,Ny,5),F,tck=-0.015)
  axis(1,seq(1,Ny,10),FINYrs[seq(1,Ny,10)],tck=-0.02,cex.axis=1.25)
  dev.off()
}
 
#G 4.20 ABARES (James Woodhams)  
if(do.ABARES="YES")
{
  #Catch by species and financial year for 2006-07 to 2014-15
  Abr.yrs=c("2006-07","2007-08","2008-09","2009-10","2010-11",
            "2011-12","2012-13","2013-14","2014-15")
  D.mon=subset(Data.monthly,FINYEAR%in%Abr.yrs & SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c))
  D.mon.nor=subset(Data.monthly.north,FINYEAR%in%Abr.yrs & SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c))
  D.abr=rbind(D.mon,D.mon.nor)
  SP.abare=D.abr[!duplicated(D.abr$SPECIES),match(c("SPECIES","SNAME"),names(D.abr))]
  Ktch.ab=aggregate(LIVEWT.c~FINYEAR+SPECIES,D.abr,sum)
  Ktch.ab=merge(Ktch.ab,SP.abare,by="SPECIES",all.x=T)
  Scien.names.abrs=read.csv("C:/Matias/Data/Comm.sp.Scientific.csv")
  Scien.names.abrs=subset(Scien.names.abrs,select=c(SPECIES,Scientific.name))
  Ktch.ab=merge(Ktch.ab,Scien.names.abrs,by="SPECIES",all.x=T)
  
  
  #Effort by financial year for TDGDLF for 2006-07 to 2014-15
  Eff.ab_1000_km_gn_d.all=subset(Total.effort.days.monthly,FINYEAR%in%Abr.yrs)
  Eff.ab_1000_km_gn_d.zn=subset(Total.effort.zone.days.monthly,FINYEAR%in%Abr.yrs)
  
  #export
  write.csv(Ktch.ab,paste( hndl,"/Ktch.ab.csv",sep=''),row.names=F)
  write.csv(Eff.ab_1000_km_gn_d.all,paste( hndl,"/Eff.ab_1000_km_gn_d.all.csv",sep=""),row.names=F)
  write.csv(Eff.ab_1000_km_gn_d.zn,paste( hndl,"/Eff.ab_1000_km_gn_d.zn.csv",sep=""),row.names=F)
}

#G 4.21 Carlie Telfer
if(do.Carlie.Telfer=="YES")
{  
  hndls=paste(hndl,"/Carlie_Telfer/",sep="")
  Carlie.Telfer.Eff.km.gn.d=aggregate(Km.Gillnet.Days.c~block10+Fishing_yr,Jodie.Effort,sum)
  Carlie.Telfer.Eff.km.gn.h=aggregate(Km.Gillnet.Hours.c~block10+Fishing_yr,Jodie.Effort,sum)
  
  write.csv(Carlie.Telfer.Eff.km.gn.d,paste(hndls,"Eff.km.gn.d.csv",sep=""),row.names=F)
  write.csv(Carlie.Telfer.Eff.km.gn.h,paste(hndls,"Eff.km.gn.h.csv",sep=""),row.names=F)
  
}

#G 4.22 Cetacean Interactions
if(do.Cetacean.Inter=="YES")
{  
  hndls=paste(hndl,"/Cetacean_Inter/",sep="")
  
  #confidential=2
  confidential=0  #provide all records to Rory, so he can decide
  
  #GILLNETS
  KTCH.ly=subset(Data.monthly,METHOD=="GN" & Estuary=="NO" & NETLEN.c>=100 & LAT<=(-26),
                      select=c(Same.return,FINYEAR,MONTH,VESSEL,METHOD,
                               SPECIES,SNAME,LIVEWT.c,BLOCKX,LAT,LONG))
  
  KTCH.ly_daily=subset(Data.daily,METHOD=="GN" & Estuary=="NO" & netlen.c>=100 & LAT<=(-26),
                 select=c(Same.return,FINYEAR,MONTH,VESSEL,METHOD,zone,
                          block10,BLOCKX,SPECIES,SNAME,LIVEWT.c,LAT,LONG))
  
  daily.yrs=unique(Effort.daily$finyear)
  
      #monthly
  SHOTS_c=unique(KTCH.ly$Same.return)
  
  KTCH.ly=subset(KTCH.ly,!FINYEAR%in%daily.yrs)
  KTCH.ly$N=1
  KTCH.ly.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+VESSEL,data=KTCH.ly,max,na.rm=T)
  KTCH.ly.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX,data=KTCH.ly,sum,na.rm=T)
  Remove.these=subset(KTCH.ly.vessels,N<=confidential)
  
  KTCH.ly=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX,data=KTCH.ly,sum,na.rm=T)
  KTCH.ly$Fishery="TDGDLF"
  KTCH.ly$Gear="gillnet"
  names(KTCH.ly)[match("LIVEWT.c",names(KTCH.ly))]="Catch(kg)"
  
  if(nrow(Remove.these)>0)Remove.these$dummy=with(Remove.these,paste(FINYEAR, MONTH, BLOCKX))
  KTCH.ly$dummy=with(KTCH.ly,paste(FINYEAR, MONTH, BLOCKX))
  if(nrow(Remove.these)>0)KTCH.ly=subset(KTCH.ly,!dummy%in%Remove.these$dummy)
  KTCH.ly=KTCH.ly[,-match("dummy",names(KTCH.ly))]
  
  write.csv(KTCH.ly,paste(hndls,"Catch_monthly_records_TDGDLF.csv",sep=""),row.names=F)
  
        #effort
  Eff.ly=subset(Effort.monthly,!FINYEAR%in%daily.yrs & METHOD=="GN")
  Eff.ly=subset(Eff.ly,Same.return%in%SHOTS_c)
  Eff.ly=aggregate(Km.Gillnet.Days.c~Same.return+FINYEAR+MONTH+BLOCKX+VESSEL,data=Eff.ly,max,na.rm=T) #remove duplicates
  Eff.ly=aggregate(Km.Gillnet.Days.c~FINYEAR+MONTH+BLOCKX+VESSEL,data=Eff.ly,sum,na.rm=T)
  Eff.ly$Fishery="TDGDLF"
  Eff.ly$Gear="gillnet"
  names(Eff.ly)[match("Km.Gillnet.Days.c",names(Eff.ly))]="Effort(Km.Gillnet.Days)"
  
  if(nrow(Remove.these)>0)
  {
    Eff.ly$dummy=with(Eff.ly,paste(FINYEAR, MONTH, BLOCKX))
    Eff.ly=subset(Eff.ly,!dummy%in%Remove.these$dummy)
    Eff.ly=Eff.ly[,-match("dummy",names(Eff.ly))]
    Percent.removed.monthly=100*nrow(Remove.these)/nrow(KTCH.ly.vessels)
    write.csv(Percent.removed,paste(hndls,"Percent.removed.csv",sep=""),row.names=F)
  }
  
  write.csv(Eff.ly,paste(hndls,"Efforty_monthly_records_TDGDLF.csv",sep=""),row.names=F)
  
    #daily
  #note: if using block10, 98% of the data had less than 3 vessels
  #       if using BLOCKX, 77% of the data had less than 3 vessels, hence used ZONE
  
  SHOTS_c=unique(KTCH.ly_daily$Same.return)
  
  KTCH.ly_daily$N=1
  KTCH.ly_daily.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+block10+VESSEL,data=KTCH.ly_daily,max,na.rm=T)
  KTCH.ly_daily.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+block10,data=KTCH.ly_daily.vessels,sum,na.rm=T)
  if(nrow(Remove.these)>0)Remove.these_daily=subset(KTCH.ly_daily.vessels,N<=confidential)
  
  KTCH.ly_daily=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+block10+VESSEL,data=KTCH.ly_daily,sum,na.rm=T)
  KTCH.ly_daily$Fishery="TDGDLF"
  KTCH.ly_daily$Gear="gillnet"
  names(KTCH.ly_daily)[match("LIVEWT.c",names(KTCH.ly_daily))]="Catch(kg)"
  
  if(nrow(Remove.these)>0)Remove.these_daily$dummy=with(Remove.these_daily,paste(FINYEAR, MONTH,block10,VESSEL))
  KTCH.ly_daily$dummy=with(KTCH.ly_daily,paste(FINYEAR, MONTH, block10,VESSEL))
  if(nrow(Remove.these)>0)KTCH.ly_daily=subset(KTCH.ly_daily,!dummy%in%Remove.these_daily$dummy)
  KTCH.ly_daily=KTCH.ly_daily[,-match("dummy",names(KTCH.ly_daily))]
  
  write.csv(KTCH.ly_daily,paste(hndls,"Catch_daily_records_TDGDLF.csv",sep=""),row.names=F)
  

    #effort
  Eff.ly_daily=subset(Effort.daily,Same.return%in%SHOTS_c & method=="GN")
  Use.Date="YES"    #Rory's approach (aggregating by DATE)
  #Use.Date="NO"     # aggregating by SNo and DSNo 
  if(Use.Date=="NO")  Eff.ly_daily=aggregate(Km.Gillnet.Days.c~ID+Same.return+finyear+month+blockx+block10+vessel,data=Eff.ly_daily,max,na.rm=T)
  if(Use.Date=="YES") Eff.ly_daily=aggregate(Km.Gillnet.Days.c~date+Same.return+finyear+month+blockx+block10+vessel,data=Eff.ly_daily,max,na.rm=T) 
  Eff.ly_daily=aggregate(Km.Gillnet.Days.c~finyear+month+blockx+block10+vessel,data=Eff.ly_daily,sum,na.rm=T)
  Eff.ly_daily$Fishery="TDGDLF"
  Eff.ly_daily$Gear="gillnet"
  names(Eff.ly_daily)[match("Km.Gillnet.Days.c",names(Eff.ly_daily))]="Effort(Km.Gillnet.Days)"
  
  if(nrow(Remove.these)>0)
  {
    Eff.ly_daily$dummy=with(Eff.ly_daily,paste(finyear, month, block10,vessel))
    Eff.ly_daily=subset(Eff.ly_daily,!dummy%in%Remove.these_daily$dummy)
    Eff.ly_daily=Eff.ly_daily[,-match("dummy",names(Eff.ly_daily))]
    Percent.removed.daily=100*nrow(Remove.these_daily)/nrow(KTCH.ly_daily.vessels)
    write.csv(Percent.removed.daily,paste(hndls,"Percent.removed.daily.csv",sep=""),row.names=F)
  }
  write.csv(Eff.ly_daily,paste(hndls,"Efforty_daily_records_TDGDLF.csv",sep=""),row.names=F)
  
 
}

#G 4.23 Scalefish catch by year for metro closure
if(do.Jodies.metro.closure=="YES")
{
  this.yr.Jodie=c("2002-03","2003-04","2004-05","2005-06","2006-07")
  Day.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,
                    data=subset(Data.daily,SPECIES%in%Scalefish.species & LAT<=(-31)& LAT>=(-33) & 
                                  zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)
  Mn.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,
                   data=subset(Data.monthly,SPECIES%in%Scalefish.species & LAT<=(-31)& LAT>=(-33) & 
                                 zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)  
  
  Day.Jod$METHOD=as.character(Day.Jod$METHOD) 
  Mn.Jod$BLOCKX=as.numeric(Mn.Jod$BLOCKX) 
  Jod=rbind(Day.Jod,Mn.Jod)
  Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,data=Jod,sum,na.rm=T)  
  write.csv(Jod,paste(hndl,"/Scalefish_metro_closure_catch_TDGDLF.csv",sep=""),row.names=F)
  
}

#G 4.24 Shark and ray catch by year for metro closure
if(do.Jodies.metro.closure=="YES")
{
  this.yr.Jodie=c("2002-03","2003-04","2004-05","2005-06","2006-07")
  Day.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+VESSEL,
                    data=subset(Data.daily,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                                  zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)
  Mn.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+VESSEL,
                   data=subset(Data.monthly,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                                 zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)  
  
  SPE.nms.day=subset(Data.daily,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                       zone=="West" & FINYEAR%in%this.yr.Jodie,select=c(SPECIES,SNAME))
  SPE.nms.day=SPE.nms.day[!duplicated(SPE.nms.day$SPECIES),]
  Day.Jod=merge(Day.Jod,SPE.nms.day,by="SPECIES",all.x=T)
  
  SPE.nms.Mn=subset(Data.monthly,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                      zone=="West" & FINYEAR%in%this.yr.Jodie,select=c(SPECIES,SNAME))
  SPE.nms.Mn=SPE.nms.Mn[!duplicated(SPE.nms.Mn$SPECIES),]
  Mn.Jod=merge(Mn.Jod,SPE.nms.Mn,by="SPECIES",all.x=T)
  
  Day.Jod$METHOD=as.character(Day.Jod$METHOD) 
  Mn.Jod$BLOCKX=as.numeric(Mn.Jod$BLOCKX) 
  
  Jod=rbind(Day.Jod,Mn.Jod)
  
  Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,data=Jod,sum,na.rm=T)  
  write.csv(Jod,paste(hndl,"/Elasmobranch_metro_closure_catch_TDGDLF.csv",sep=""),row.names=F)  
}

#G 4.25 Total shark fin by year for Tim Nicholas
if(do.Tim_N.fin=="YES")
{
  DAT=subset(Data.monthly,FINYEAR%in%FINYEAR.daily & METHOD%in%c("GN","LL") & Estuary=="NO" &
               LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  #total fins   
  Fin.Weight=subset(DAT,SPECIES%in%Elasmo.species,select=c(SPECIES,FINYEAR,LIVEWT.c))
  Fin.Weight=subset(Fin.Weight,!(SPECIES%in%c(13000,20000,31000))) #remove no fin species (wobbies, spurdog,rays except shovelnose)
  Fin.Weight$Total.catch.fin=Percent.fin.of.livewt*Fin.Weight$LIVEWT.c
  Total.catch.fin=aggregate(Total.catch.fin~FINYEAR,Fin.Weight,sum,na.rm=T)
  write.csv(Total.catch.fin,paste(hndl,"/Shark_fins_TDGDLF.csv",sep=""),row.names=F)
}

#G 4.26 Total shark and scalefish catch by year block for ASL closure for Tim Nicholas
if(do.Tim_ASL=="YES")
{  
  dummy=subset(Data.daily,FisheryCode%in%c("SGL1","SGL2","WCGL") 
               & FINYEAR%in%c("2006-07","2007-08","2008-09","2009-10","2010-11","2011-12",
                              "2012-13","2013-14","2014-15","2015-16"))
  
  Elas=aggregate(LIVEWT.c~FisheryCode+day+MONTH+FINYEAR+block10,subset(dummy,SPECIES%in%Elasmo.species),sum)
  Scalies=aggregate(LIVEWT.c~FisheryCode+day+MONTH+FINYEAR+block10,subset(dummy,SPECIES%in%Scalefish.species),sum)
  names(Elas)[match("LIVEWT.c",names(Elas))]="Shark_ray_tot_catch_kg"
  names(Scalies)[match("LIVEWT.c",names(Scalies))]="Scalefish_tot_catch_kg"
  Tab=merge(Elas,Scalies,by=c("FisheryCode","day","MONTH","FINYEAR","block10"),all=T)
  rm(dummy)
  write.csv(Tab,paste(hndl,"/ASL_catch_block.csv",sep=""),row.names=F)
}

#G 4.27 Bull shark catches for Adrian Gleiss
if(do.Adrian.Gleiss=="YES")
{
  Bull=subset(Data.monthly,SPECIES==18021)
  Bull.daily=subset(Data.daily,SPECIES==18021)
  Bull=subset(Bull,!Same.return%in%unique(Bull.daily$Same.return))
  Bull=subset(Bull,select=c(MONTH,FINYEAR,METHOD,LAT,LONG,SNAME,LIVEWT.c))
  Bull.daily=subset(Bull.daily,select=c(day,MONTH,FINYEAR,METHOD,LAT,LONG,RSCommonName,LIVEWT.c))
  write.csv(Bull,paste(hndl,"/Bulls.csv",sep=""),row.names=F)
  write.csv(Bull,paste(hndl,"/Bulls_daily.csv",sep=""),row.names=F)
}

#G 4.28 Alexandra.Hoschke data extraction
if(do.Alexandra.Hoschke=="YES")
{
  #reported catch data
  a=subset(Data.monthly,SPECIES==8001,select=c(SPECIES,SNAME,FINYEAR,YEAR.c,MONTH,METHOD,BLOCKX,LAT,LONG,LIVEWT.c))
  #note:no greynurse reported in Daily logbooks
  
  #reported interactions
  TEPS.grey.nurse=subset(TEPS.current,SpeciesCode==8001,select=c(DailySheetNumber,DataEntryName,Status,ScientificName,Comments))
  TEPS.grey.nurse$Status=with(TEPS.grey.nurse,ifelse(Status%in%c("A","a"),"Alive",ifelse(Status%in%c("D","d"),"Dead",NA)))
  names(TEPS.grey.nurse)[match("Status",names(TEPS.grey.nurse))]="Discarding.condition"
  xx=unique(TEPS.grey.nurse$DailySheetNumber)
  # Get=subset(Data.daily.1,TSNo%in%xx,select=c(TSNo,finyear,date,depthMax,method,
  #                                             Block,block10,Lat,LatDeg,LatMin,Long,LongDeg,LongMin))
  # Get=Get[!duplicated(Get$TSNo),]
  # names(Get)[match("LAT",names(Get))]="LatDeg"
  # 
  Get=subset(Data.daily.original,DSNo%in%xx,select=c(DSNo,FINYEAR,date,depthMax,METHOD,
                                                     Block,block10,LAT,LatMin,LongDeg,LongMin))
  names(Get)[match("LAT",names(Get))]="LatDeg"
  Get=Get[!duplicated(Get$DSNo),]
  TEPS.grey.nurse=merge(TEPS.grey.nurse,Get,by.x="DailySheetNumber",by.y="DSNo",all.x=T)
  write.csv(a,paste(hndl,"/Alexandra.Hoschke.csv",sep=""),row.names=F)
  write.csv(TEPS.grey.nurse,paste(hndl,"/Alexandra.Hoschke_TEPS.csv",sep=""),row.names=F)
}

#G 4.29 ASL closures compensation
if(ASL.compensation=="YES")
{
  yr=c("2012-13","2013-14","2014-15","2015-16","2016-17")
  
  
  #catch
  a=subset(Data.daily,FINYEAR%in%yr & zone%in%c('West','Zone1','Zone2'),select=c(FINYEAR,block10,LIVEWT.c))
  Block_ktch=aggregate(LIVEWT.c~block10+FINYEAR,a,sum)
  
  #effort
  if(Use.Date=="NO")
  {
    b=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    b1=aggregate(Km.Gillnet.Hours.c~ID+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
  }
  if(Use.Date=="YES")
  {
    b=aggregate(Km.Gillnet.Days.c~date+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    b1=aggregate(Km.Gillnet.Hours.c~date+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    
  }
  b=aggregate(Km.Gillnet.Days.c~vessel+zone+finyear+block10,data=b,sum,na.rm=T)
  b1=aggregate(Km.Gillnet.Hours.c~vessel+zone+finyear+block10,data=b1,sum,na.rm=T)
  b2=merge(b,b1,by=c("vessel","zone","finyear","block10"),all.x=T)
  b2=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~block10+finyear,data=b2,sum,na.rm=T)
  
  uniqVes <- function(x)  length(unique(x))
  VeslS=aggregate(vessel~ finyear+block10, data = subset(Effort.daily,finyear%in%yr), FUN= uniqVes)
  b2=merge(b2,VeslS,by=c("finyear","block10"))
  
  
  yr_col=paste(yr[1],yr[length(yr)],sep=".to.")
  write.csv(Block_ktch,paste(hndl,"/Block_ktch_",yr_col,".csv",sep=""),row.names=F)
  write.csv(b2,paste(hndl,"/Block_effort_",yr_col,".csv",sep=""),row.names=F)
}


#G 4.30 CITES 2018
if(do.CITES=="YES")
{
  a=subset(Data.monthly,SPECIES%in%c(10001,18001,18003,18004,26999) & FINYEAR %in%FINYrs)
  a$SPECIES=with(a,ifelse(SPECIES==18001,18003,SPECIES))
  
  Tab=aggregate(LIVEWT.c~FINYEAR+SPECIES,a,sum)
  Spnms=a[!duplicated(a$SPECIES),match(c('SPECIES','SNAME'),names(a))]
  Tab=merge(Tab,Spnms,by='SPECIES')
  Tab=subset(Tab,select=c(FINYEAR,LIVEWT.c,SNAME))
  Tab1=reshape(Tab,v.names = "LIVEWT.c", idvar = "FINYEAR",
              timevar = "SNAME", direction = "wide")
  names(Tab1)[2:ncol(Tab1)]=paste(substr(names(Tab1)[2:ncol(Tab1)],10,50),"(kg)")
  Tab1=Tab1[order(Tab1$FINYEAR),]
  write.csv(Tab1,paste(hndl,"/CITES_2018.csv",sep=""),row.names=F)
  PRCES=subset(PRICES,RSCommonName%in%c("Bronze Whaler","Blue Shark","Guitarfishes","Shortfin Mako"))
  colnames(PRCES)[match("uv1617",colnames(PRCES))]="$/kg"
  write.csv(PRCES,paste(hndl,"/CITES_2018_PRICES.csv",sep=""),row.names=F)
}

#G 4.31 Changes in catch composition of WCDGDLF due to changes in mesh size
if(do.Nick_mesh.size.WCDGDLF=="YES")
{
  if("plyr" %in% (.packages())) detach(package:plyr)
  
  TAB=with(subset(Data.daily,zone=='West'),table(RSCommonName))
  fn.get=function(d)
  {
    return(WCDGDLF=d %>%
      filter(zone=='West' & RSCommonName%in% names(TAB[TAB>100])) %>%
      dplyr::select(Same.return.SNo,FINYEAR,RSCommonName,LIVEWT.c) %>% 
      group_by(RSCommonName,FINYEAR)  %>% 
      summarize(Annual.ktch_tonnes=sum(LIVEWT.c/1000)) %>%
      spread(RSCommonName, Annual.ktch_tonnes)    %>%
      as.data.frame)
  }
  
  write.csv(fn.get(d=Data.daily),paste(hndl,"/WCDGDLF_catch_compo_by_year.csv",sep=""),row.names=F)
  write.csv(fn.get(d=subset(Data.daily,METHOD=="GN")),paste(hndl,"/WCDGDLF_catch_compo_by_year_gillnet.csv",sep=""),row.names=F)
  write.csv(fn.get(d=subset(Data.daily,METHOD=="LL")),paste(hndl,"/WCDGDLF_catch_compo_by_year_longline.csv",sep=""),row.names=F)
}

#do.ASL.action.2018
if(do.ASL.action.2018=="YES")
{
  TAB=TEPS.current %>%
    filter(SpeciesCode%in%c(999902)) %>%
    group_by(finyear,month,blockx,DataEntryName) %>%
    summarise(Number = sum(Number)) %>%
    arrange(finyear,month,blockx) %>%
    spread(key = "DataEntryName", value = "Number", fill = 0) %>%
    as.data.frame
  write.csv(TAB,paste(hndl,"/ASL_action2_2018.csv",sep=""),row.names=F)
  
} 

if(do.Parks.Australia=="YES")  
{
  South.WA.lat=c(-36,-25); South.WA.long=c(112,130)
  PLATE=c(.01,.9,.075,.9)
  Yrs=c("2017-18","2018-19")
  aa= Data.daily.original%>%filter(FINYEAR%in%Yrs) %>%
    mutate(LatDeg=as.numeric(substr(block10,1,2)),
           LatMin=ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin),
           Lat=-abs(LatDeg+(LatMin/60)),
           LongDeg=ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg), 
           LongMin=ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin),
           Long=LongDeg+(LongMin/60))%>%
    filter(Lat<=(-26) & Lat>(-36.5)& Long<=(129) & Long >(111.9))
  
  numInt=20
  couleurs=rev(heat.colors(numInt)) 
  tcl.1=.5
  tcl.2=.5
  Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
  Lat.seq=c(-26,-28,-30,-32,-34)
  numberLab=5
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  
  #lodged returns
  TAB1=aa %>% filter(FINYEAR%in%Yrs) %>%
    group_by(FINYEAR)%>%
    summarise(Unique_TSNo=n_distinct(TSNo))%>%
    data.frame
  
  TAB2=aa %>% filter(FINYEAR%in%Yrs) %>%
    group_by(FINYEAR,METHOD)%>%
    summarise(Unique_TSNo=n_distinct(TSNo))%>%
    data.frame
  
  #Who's been using hooks?
  TAB3_trips=aa%>%group_by(METHOD,FINYEAR,VESSEL)%>%
    summarise(Trips=n_distinct(TSNo))%>%
    spread(METHOD, Trips)%>%
    replace(is.na(.), "")%>%
    data.frame

  
  bb=aa%>%filter(METHOD=="LL")%>%
          distinct(Same.return.SNo,.keep_all =T) %>%
          select(VESSEL,BoatName,MastersName,port,block10,FINYEAR,MONTH,bioregion,Lat,Long,depthMax,
                 NilCatch,species,nfish,livewt,
                 HookSize,HookType,HOOKS,HOURS,nlines,SHOTS)
  
  TAB4=bb%>%group_by(VESSEL,BoatName,MastersName,port)%>%
            summarise(mean.hook.n=mean(HOOKS,na.rm=T),
                      mean.hook.size=mean(HookSize,na.rm=T),
                      mean.hook.hours=mean(HOURS,na.rm=T))%>%
              replace(is.na(.), "")%>%
            data.frame
  library(gridExtra)
  library(grid)
  mytheme <- gridExtra::ttheme_default(
    base_size = 10,
    core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .65)),
    colhead = list(fg_params=list(cex = .75)),
    rowhead = list(fg_params=list(cex = .75)))
  
  pdf(file=paste(hndl,"/Parks Australia/Parks_Australia_2018-19.effort_catch.pdf",sep=""))
  
  grid.draw(gridExtra::tableGrob(TAB3_trips, theme = mytheme,rows = NULL))
  grid.newpage()
  
  grid.draw(gridExtra::tableGrob(TAB4, theme = mytheme,rows = NULL))
  
  #effort
  b=aa %>% filter(METHOD=="GN") %>%
    mutate(Km.Gillnet.Hours=HOURS*NETLEN/1000)%>%
    filter(Km.Gillnet.Hours>0)%>%
    group_by(Same.return.SNo,FINYEAR, block10) %>%
    summarize(Km.Gillnet.Hours = max(Km.Gillnet.Hours, na.rm = TRUE))%>%
    group_by(FINYEAR, block10) %>%
    summarize(sum = sum(Km.Gillnet.Hours, na.rm = TRUE))%>%
    mutate(LatDeg=as.numeric(substr(block10,1,2)),
           LatMin=10*as.numeric(substr(block10,3,3)),
           Lat=-abs(LatDeg+(LatMin/60)),
           LongDeg=100+as.numeric(substr(block10,4,5)), 
           LongMin=10*as.numeric(substr(block10,6,6)),
           Long=LongDeg+(LongMin/60))%>%
    data.frame
  b=subset(b,select=c(FINYEAR,sum,Lat,Long))
  BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
  
  par(mfrow=c(2,1),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
  for(y in 1:length(Yrs))
  {
    bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
      arrange(Lat,Long)
    long=sort(unique(bb$Long))
    lat=sort(unique(bb$Lat))      
    Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]	
    
    a=South.WA.long[1]:South.WA.long[2]
    bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
    axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
    axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                         rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
    nnn=with(TAB2%>%filter(FINYEAR==Yrs[y]),paste(paste(Yrs[y]," (gillnet returns= ",Unique_TSNo[1],"; longline returns= ",Unique_TSNo[2],")",sep="")))
    mtext(nnn,3,-2)
  }
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-1,las=3,cex=1.1,outer=T)
  mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=-.5,cex=1.1,outer=T)
  mtext("Effort (Km.gn.hours)",3,-0.75,outer=T)
  
  #catch
  for(s in 1:length(TARGETS))
  {
    b=aa %>% filter(METHOD=="GN" & species%in%TARGETS[[s]]) %>%
      group_by(FINYEAR, block10) %>%
      summarize(sum = sum(livewt, na.rm = TRUE))%>%
      mutate(LatDeg=as.numeric(substr(block10,1,2)),
             LatMin=10*as.numeric(substr(block10,3,3)),
             Lat=-abs(LatDeg+(LatMin/60)),
             LongDeg=100+as.numeric(substr(block10,4,5)), 
             LongMin=10*as.numeric(substr(block10,6,6)),
             Long=LongDeg+(LongMin/60))%>%
      data.frame
    b=subset(b,select=c(FINYEAR,sum,Lat,Long))
    BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
    
    par(mfrow=c(2,1),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
    for(y in 1:length(Yrs))
    {
      bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
        arrange(Lat,Long)
      long=sort(unique(bb$Long))
      lat=sort(unique(bb$Lat))      
      Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      a=South.WA.long[1]:South.WA.long[2]
      bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                           rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
      mtext(Yrs[y],3,-2)
    }
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-1,las=3,cex=1.1,outer=T)
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=-.5,cex=1.1,outer=T)
    mtext(paste(names(TARGETS)[s],"(catch in kg)"),3,-.75,outer=T)
    
    
  }
  
  #all scalefish
  {
    b=aa %>% filter(METHOD=="GN" & species%in%Scalefish.species) %>%
      group_by(FINYEAR, block10) %>%
      summarize(sum = sum(livewt, na.rm = TRUE))%>%
      mutate(LatDeg=as.numeric(substr(block10,1,2)),
             LatMin=10*as.numeric(substr(block10,3,3)),
             Lat=-abs(LatDeg+(LatMin/60)),
             LongDeg=100+as.numeric(substr(block10,4,5)), 
             LongMin=10*as.numeric(substr(block10,6,6)),
             Long=LongDeg+(LongMin/60))%>%
      data.frame
    b=subset(b,select=c(FINYEAR,sum,Lat,Long))
    BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
    
    par(mfrow=c(2,1),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
    for(y in 1:length(Yrs))
    {
      bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
        arrange(Lat,Long)
      long=sort(unique(bb$Long))
      lat=sort(unique(bb$Lat))      
      Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      a=South.WA.long[1]:South.WA.long[2]
      bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                           rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
      mtext(Yrs[y],3,-2)
    }
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-1,las=3,cex=1.1,outer=T)
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=-.5,cex=1.1,outer=T)
    mtext("Scalefish (catch in kg)",3,-.75,outer=T)
    
    
  }
  dev.off()
  
}


########### SECTION H. ----  EXPORT TOTAL CATCH FOR REFERENCE POINT ANALYSIS --- ###########
if(do.Ref.Points=="YES")
{
  for (i in 1:length(TARGETS))write.csv(STORE.Tot.Catch[[i]],
    paste("C:/Matias/Analyses/Reference Points/Tot.c.",sp[i],".csv",sep=""),row.names=F)  
}



###########  SECTION I. ----  EXPLORATORY ANALYSES ---   (including movies) ###########
if(do.exploratory=="YES")
{
  setwd("C:/Matias/Analyses/Catch and effort/Outputs/Exploratory")
  
  #5.1 Annual spatial expansion
  tiff(file="Figure 1.Proportion of blocks fished.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(las=1,mgp=c(2.5,.85,0))
  plot(1:NN.monthly,Spatial.expan,ylim=c(0,100),ylab="Percent of blocks used", xlab="Financial year",
       pch=19,cex=1.75,xaxt='n',cex.axis=1.25,cex.lab=1.75)
  axis(1,at=1:NN.monthly,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN.monthly,5),labels=FINYEAR.monthly[seq(1,NN.monthly,5)],tck=-0.02,cex.axis=1.1)
  dev.off()
  
  
  
  #5.2 Effort expansion and licence holder fishing 
  PLOT.fn=function(DATA1,DATA2,max1,max2,LEG)
  {
    plot(1:NN.monthly,DATA1,ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.5,ylim=c(0,max1)
         ,cex.axis=1.25,lwd=1.75)
    axis(1,at=1:NN.monthly,labels=F,tck=-0.015)
    
    par(new=T)
    plot(1:NN.monthly,DATA2,col="grey60",type='o',pch=19,axes=F,ann='F',cex=1.5,lwd=1.75,ylim=c(0,max2))
    axis(4,at=pretty(DATA2),labels=pretty(DATA2),las=2,cex.axis=1.3)
    legend("topleft",LEG,bty='n',cex=1.5)
    
  }
  
  tiff(file="Figure 2. Effort expansion.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,1),mar=c(1,3.6,.1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.1,.7,0))
  PLOT.fn(Spatial.expan,Effort.expan$N.ves.yr,80,250,"Temperate Gillent and Longline")
  PLOT.fn(Spatial.expan.north,Effort.expan.north$N.ves.yr,80,100,"Northern Shark Fisheries")
  axis(1,at=seq(1,NN.monthly,5),labels=FINYEAR.monthly[seq(1,NN.monthly,5)],tck=-0.03,cex.axis=1.1)
  mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.7,outer=T)
  mtext("Number of blocks fished",side=2,line=-1.2,font=1,las=0,cex=1.7,outer=T)
  mtext("Number of licence holders fishing",side=4,line=-1.2,las=3,cex=1.7,col="grey60",outer=T)
  dev.off()
  
  
  
  
  #5.3 Spatial and temporal expansion
  
  tiff(file="Figure 3.Effort_map.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  
  par(mfrow=c(3,3),mai = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
  plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
          col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
  text(133,-21.5,("Australia"),col="black", cex=2)
  text(118.7,-32,("Perth"),col="black", cex=1.1)
  mtext("Latitude (?S)",side=2,line=0.4,las=3,cex=1.3)
  mtext("Longitude (?E)",side=1,line=0.6,cex=1.3)
  for (i in 1:length(DATA.lista))
  {
    fn.eff.plot(DATA.lista[[i]],tcl.1=16,tcl.2=18,EffortBreaks)
    mtext(Yr.range[i],side=3,line=0,cex=.95)
    if(i%in%6:8) axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
    if(i%in%c(1,3,6)) axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    if(i==8)color.legend(126,-26,129,-30.5,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
  }
  dev.off()
  
 
  #5.5 Species distribution 
  fn.dist=function(SPEC,YR)
  {
    Dat=subset(Data.daily.original,species%in%SPEC & FINYEAR==YR &  LAT>-40 & METHOD=="GN")
    Dat$lati=-with(Dat,-LAT+(LatMin/60))
    Dat$longi=with(Dat,LongDeg+(LongMin/60))
    Dat$Km.gn.hours=Dat$NETLEN*Dat$HOURS*Dat$SHOTS
    quilt.plot(Dat$longi,Dat$lati,Dat$livewt/Dat$Km.gn.hours)
  }
  YR.dist=FINYEAR.monthly[match("2006-07",FINYEAR.monthly):length(FINYEAR.monthly)]
  par(mfcol=c(2,3))
  #for(i in 1:length(YR.dist))fn.dist(17001,YR.dist[i])
  
  
  
  #Depth distributions
  fn.depth.dist=function(DATA,SP)
  {
    DATA=subset(DATA,species%in%SP&depthMax<150&depthMax>0,select=c(depthMax,livewt,species))
    hist(DATA$depthMax,main=sp[ss])
  }
  par(mfcol=c(2,2))
  for (ss in 1:length(Tar))fn.depth.dist(DATA=Data.daily.original,SP=Tar[ss])
  
  
  
  #MOVIES
  # Effort movie
  MAX.EFF.movie=aggregate(Km.Gillnet.Days.c~BLOCKX+YEAR.c,
                          data=Data.monthly.GN[-which(duplicated(Data.monthly.GN$Same.return)),],FUN=sum,na.rm=T)
  EffortBreaks.movie=quantile(MAX.EFF.movie[,3],probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
  
  ani.options(nmax=length(YEAR.c.monthly),outdir = paste(getwd(),"/Effort.movie",sep=""))
  saveGIF({
    for (i in 1:length(YEAR.c.monthly))
    {
      fn.eff.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),tcl.1=.5,tcl.2=.5,EffortBreaks.movie)
      mtext(YEAR.c.monthly[i],side=3,line=0,cex=1.25)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      #  color.legend(126,-26,129,-30.5,round(EffortBreaks.movie,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
      color.legend(126,-26,129,-30.5,"",rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
      mtext("Latitude (?S)",side=2,line=3,las=3,cex=1.25)
      mtext("Longitude (?E)",side=1,line=3,cex=1.25)    
    }
  }, interval=0.6,movie.name="Effort.movie.wmv", ani.width=300,ani.height=300)
  
  tiff(file="Effort.movie/Effort_movie.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfrow=c(7,6),mai = c(0.15, 0.15, 0.15, 0.05))
  for (i in 1:(length(YEAR.c.monthly)-1))
  {
    fn.eff.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),tcl.1=.5,tcl.2=.5,EffortBreaks.movie)
    mtext(YEAR.c.monthly[i],side=3,line=0,cex=.75)
    #  axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
    #  axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    if(i==36)color.legend(126,-26,129,-30.5,round(EffortBreaks.movie,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.35)
  }
  dev.off()
  
  
  
  # Catch movie
  fn.Catch.mov.breaks=function(SP)
  {
    MAX.Catch.movie=aggregate(LIVEWT.c~BLOCKX+YEAR.c,
                              data=subset(Data.monthly.GN,SPECIES%in%SP),FUN=sum,na.rm=T)
    return(quantile(MAX.Catch.movie[,3],probs=seq(0,1,1/numInt),na.rm=T))
  }
  
  
  for (j in 1:length(TARGETS))
  {
    SP=Tar[j]
    Breaks=fn.Catch.mov.breaks(SP) 
    ani.options(nmax=length(YEAR.c.monthly),outdir = paste(getwd(),"/Catch.movie",sep=""))
    saveGIF({
      for (i in 1:length(YEAR.c.monthly))
      {
        fn.catch.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),SP,tcl.1=.5,tcl.2=.5,Breaks)
        mtext(YEAR.c.monthly[i],side=3,line=0,cex=1.25)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
        mtext("Latitude (?S)",side=2,line=3,las=3,cex=1.25)
        mtext("Longitude (?E)",side=1,line=3,cex=1.25)    
      }
    }, interval=0.6,movie.name=paste(SP,".Catch.movie.wmv",sep=""), ani.width=300,ani.height=300)
    
    
    tiff(file=paste("Catch.movie/",SP,".Catch_movie.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfrow=c(6,7),mai = c(0.15, 0.15, 0.15, 0.05))
    for (i in 1:(length(YEAR.c.monthly)-1))
    {
      fn.catch.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),SP,tcl.1=.5,tcl.2=.5,Breaks)
      mtext(YEAR.c.monthly[i],side=3,line=0,cex=.75)
      #  axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      #  axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    dev.off()
    
    
    
  }
  
  
  #CPUE movie
  fn.CPUE.mov.breaks=function(SP)
  {
    MAX.CPUE.movie=aggregate(LIVEWT.c/Km.Gillnet.Hours.c~BLOCKX+YEAR.c,
                             data=subset(Data.monthly.GN,SPECIES==SP),FUN=mean,na.rm=T)
    names(MAX.CPUE.movie)[3]="CPUE"
    MAX.CPUE.movie=subset(MAX.CPUE.movie,CPUE<1000)
    return(quantile(MAX.CPUE.movie[,3],probs=seq(0,1,1/numInt),na.rm=T))
  }
  fn.CPUE.plot=function(DATA,SP,tcl.1,tcl.2,BREAKS)
  {
    DATA=subset(DATA,LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9) & SPECIES==SP)
    
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    
    if(nrow(DATA)<=2)   plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    if(nrow(DATA)>2)
    {
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      
      MapCatch=with(DATA,aggregate(LIVEWT.c/Km.Gillnet.Hours.c,list(BLOCKX.c),FUN=mean,na.rm=T))    
      colnames(MapCatch)=c("BLOCKX.c","CPUE")
      id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
      MapCatch$LAT=DATA$LAT[id]
      MapCatch$LONG=DATA$LONG[id]
      
      MapCatch$LAT.cen=MapCatch$LAT-.5
      MapCatch$LONG.cen=MapCatch$LONG+.5  
      MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
      MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapCatch$LONG.cen))
      lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
      
      MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","CPUE"),names(MapCatch))]  
      
      
      Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="CPUE", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]										
      
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      
    }
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                 nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
  }
  
  for (j in 1:length(TARGETS))
  {
    SP=TARGETS[j]
    Breaks=fn.CPUE.mov.breaks(SP) 
    ani.options(nmax=length(YEAR.c.monthly),outdir = paste(getwd(),"/CPUE.movie",sep=""))
    saveGIF({
      for (i in 1:length(YEAR.c.monthly))
      {
        fn.CPUE.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),SP,tcl.1=.5,tcl.2=.5,Breaks)
        mtext(YEAR.c.monthly[i],side=3,line=0,cex=1.25)
        axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
        axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
        mtext("Latitude (?S)",side=2,line=3,las=3,cex=1.25)
        mtext("Longitude (?E)",side=1,line=3,cex=1.25)    
      }
    }, interval=0.6,movie.name=paste(SP,".CPUE.movie.wmv",sep=""), ani.width=300,ani.height=300)
    
    
    tiff(file=paste("CPUE.movie/",SP,".CPUE_movie.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfrow=c(6,6),mai = c(0.15, 0.15, 0.15, 0.05))
    for (i in 1:(length(YEAR.c.monthly)-1))
    {
      fn.CPUE.plot(subset(Data.monthly.GN,YEAR.c==YEAR.c.monthly[i]),SP,tcl.1=.5,tcl.2=.5,Breaks)
      mtext(YEAR.c.monthly[i],side=3,line=0,cex=.75)
      #  axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      #  axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    dev.off()
    
    
    
  }
  
 
  
  ####Effort and Catch and fine scale (10 by 10 nm blocks)
  
  # 1. Aggregate by finyear-month-vessel-method-block
  Fine.scale.eff$Same.return=with(Fine.scale.eff,paste(finyear,month,vessel,method,block10))
  Fine.scale.eff.agg.catch=aggregate(cbind(livewt,landwt)~finyear+year+month+vessel+method+block10+species+sname1+
                                       LatDeg+LatMin+LongDeg+LongMin+Same.return,data=Fine.scale.eff,sum,na.rm=T)
  
  Fine.scale.depth=Fine.scale.eff[,match(c("Same.return","depthMax","depthMin"),names(Fine.scale.eff))]
  Fine.scale.depth=Fine.scale.depth[!duplicated(Fine.scale.depth$Same.return),]
  
  Fine.scale.eff.agg.effort=aggregate(cbind(hours, shots, netlen)~finyear+year+month+vessel+method+
                                        block10+Same.return,data=Fine.scale.eff,mean,na.rm=T)
  Fine.scale.eff.agg.effort=Fine.scale.eff.agg.effort[,match(c("Same.return","hours","shots","netlen"),
                                                             names(Fine.scale.eff.agg.effort))]
  
  Fine.scale.eff.Bdays=aggregate(cbind(bdays)~finyear+year+month+vessel+method+
                                   block10+Same.return,data=Fine.scale.eff,mean,na.rm=T)
  Fine.scale.eff.Bdays=Fine.scale.eff.Bdays[,match(c("Same.return", "bdays"),names(Fine.scale.eff.Bdays))]
  
  
  Fine.scale.eff.agg=merge(Fine.scale.eff.agg.catch,Fine.scale.eff.agg.effort,by="Same.return",all.x =T)
  Fine.scale.eff.agg=merge(Fine.scale.eff.agg,Fine.scale.eff.Bdays,by="Same.return",all.x =T)
  Fine.scale.eff.agg=merge(Fine.scale.eff.agg,Fine.scale.depth,by="Same.return",all.x =T)
  
  #effort
  Fine.scale.eff.agg$Effort=with(Fine.scale.eff.agg,hours*shots*bdays*netlen/1000)
  Fine.scale.eff.agg$LatDeg=-Fine.scale.eff.agg$LatDeg
  
  Yr.group.fine.scale=seq(2006,2011)
  DATA.lista.Fine.scale.eff=vector('list',length(Yr.group.fine.scale))
  names(DATA.lista.Fine.scale.eff)=Yr.group.fine.scale
  for(i in 1:length(DATA.lista.Fine.scale.eff)) 
  {DATA.lista.Fine.scale.eff[[i]]=subset(Fine.scale.eff.agg,year==Yr.group.fine.scale[i])}
  
  #max block effort for each year period
  MAX.EFF=data.frame(Period=Yr.group.fine.scale,Max.eff=NA)
  EFF.bin=NULL
  for (i in 1:length(DATA.lista.Fine.scale.eff))
  {
    DATA=DATA.lista.Fine.scale.eff[[i]]
    DATA=DATA[-which(duplicated(DATA$Same.return)),]
    Max.effort=with(DATA,aggregate(Effort,list(block10),FUN=sum,na.rm=T))
    EFF.bin=c(EFF.bin,Max.effort[,2])
    MAX.EFF[i,2]=max(Max.effort[,2])
  }
  Max.effort=ceiling(max(MAX.EFF[,2])) 
  EffortBreaks=quantile(EFF.bin,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
  
  
  #ADD centroid lat and long for Block10
  Position.Blk.10=BlOCK_10
  names(Position.Blk.10)=c("block10","LAT.cen","LONG.cen")
  
  
  
  fn.eff.plot.fine=function(DATA,tcl.1,tcl.2,EffortBreaks)
  {
    DATA=subset(DATA,LatDeg<=(-26) & LatDeg>(-36.1)&LongDeg<=(129) & LongDeg >(111.9))
    id=which(duplicated(DATA$Same.return))
    DATA=DATA[-id,]
    
    #add Lat and Long
    DATA=merge(DATA,Position.Blk.10,by="block10")
    
    
    MapEffort=with(DATA,aggregate(Effort,list(block10),FUN=sum,na.rm=T))    
    colnames(MapEffort)=c("block10","Total effort")
    id=unique(match(MapEffort$block10,DATA$block10))
    MapEffort$LAT=DATA$LAT[id]
    MapEffort$LONG=DATA$LONG[id]
    
    MapEffort=merge(MapEffort,Position.Blk.10[,match(c("block10","LAT.cen"),names(Position.Blk.10))],by="block10")
    MapEffort=merge(MapEffort,Position.Blk.10[,match(c("block10","LONG.cen"),names(Position.Blk.10))],by="block10")
    
    MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
    MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
    long=sort(unique(MapEffort$LONG.cen))
    lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image	
    
    MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Total effort"),names(MapEffort))]  
    
    Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",		#transposed as matrix 	
                               timevar="LAT.cen",v.names="Total effort", direction="wide"))	
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]										
    
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=EffortBreaks,axes = FALSE,add=T)			
    axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
    axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
    
    par(new=T)
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                 nlev = 3,labcex=0.5,lty = c(1,2,3),col=c("gray20","gray40","gray60"),add=T)
  }
  
  tiff(file="Figure.Effort_map.fine.scale.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfrow=c(3,2),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
  for (i in 1:length(DATA.lista.Fine.scale.eff))
  {
    fn.eff.plot.fine(DATA.lista.Fine.scale.eff[[i]],tcl.1=16.5,tcl.2=25.5,EffortBreaks)
    mtext(Yr.group.fine.scale[i],side=3,line=0,cex=.95)
    if(i%in%5:6) axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
    if(i%in%c(1,3,5)) axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    if(i==6)color.legend(126,-26,129,-30.5,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
  }
  mtext("Latitude (?S)",side=2,line=-1,las=3,cex=1.1,outer=T)
  mtext("Longitude (?E)",side=1,line=-1,cex=1.1,outer=T)
  dev.off()
  
  
  #CATCH FINER RESOLUTION
  #catch breaks
  MAX.CATCH=data.frame(Period=Yr.group.fine.scale,Max.catch=NA)
  CATCH.bin=NULL
  fn.catch.breaks.fine=function(SP)
  {
    
    for (j in 1:length(DATA.lista.Fine.scale.eff))
    {
      DATA=subset(DATA.lista.Fine.scale.eff[[j]],LatDeg<=(-26) & LatDeg>(-36.1)&LongDeg<=(129) & LongDeg >(111.9) & species==SP)
      if(nrow(DATA)>0)
      {
        Max.catch=with(DATA,aggregate(livewt,list(block10),FUN=sum,na.rm=T))
        CATCH.bin=c(CATCH.bin,Max.catch[,2])
        MAX.CATCH[j,2]=max(Max.catch[,2])
      }
    }
    
    Breaks=quantile(CATCH.bin,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
    return(Breaks)
  }
  
  #function plot catch by period
  fn.catch.plot.fine=function(DATA,SP,tcl.1,tcl.2,BREAKS)
  {
    DATA=subset(DATA,LatDeg<=(-26) & LatDeg>(-36.1)&LongDeg<=(129) & LongDeg >(111.9))
    id=which(duplicated(DATA$Same.return))
    DATA=DATA[-id,]
    
    #add Lat and Long
    DATA=merge(DATA,Position.Blk.10,by="block10")
    
    
    MapEffort=with(DATA,aggregate(livewt,list(block10),FUN=sum,na.rm=T))    
    colnames(MapEffort)=c("block10","Catch")
    id=unique(match(MapEffort$block10,DATA$block10))
    MapEffort$LAT=DATA$LAT[id]
    MapEffort$LONG=DATA$LONG[id]
    
    MapEffort=merge(MapEffort,Position.Blk.10[,match(c("block10","LAT.cen"),names(Position.Blk.10))],by="block10")
    MapEffort=merge(MapEffort,Position.Blk.10[,match(c("block10","LONG.cen"),names(Position.Blk.10))],by="block10")
    
    MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
    MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
    long=sort(unique(MapEffort$LONG.cen))
    lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image  
    
    MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Catch"),names(MapEffort))]  
    
    Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",  	#transposed as matrix 	
                               timevar="LAT.cen",v.names="Catch", direction="wide"))	
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]										
    
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
    axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
    axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
    
    par(new=T)
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                 nlev = 3,labcex=0.5,lty = c(1,2,3),col=c("gray20","gray40","gray60"),add=T)
  }
  
  
  for (ss in 1:length(TARGETS))
  {
    tiff(file=paste("Figure.",TARGETS[ss],".Catch_map.fine.scale.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    
    par(mfrow=c(3,2),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
    
    for (i in 1:length(DATA.lista.Fine.scale.eff))
    {
      Breaks=fn.catch.breaks.fine(TARGETS[ss])
      fn.catch.plot.fine(DATA.lista.Fine.scale.eff[[i]],TARGETS[ss],tcl.1=16.5,tcl.2=25.5,Breaks)
      mtext(Yr.group.fine.scale[i],side=3,line=0,cex=.95)
      if(i%in%5:6) axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      if(i%in%c(1,3,5)) axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      if(i==6)color.legend(126,-26,129,-30.5,round(Breaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
    }
    mtext("Latitude (?S)",side=2,line=-1,las=3,cex=1.1,outer=T)
    mtext("Longitude (?E)",side=1,line=-1,cex=1.1,outer=T)
    dev.off()
    
  }
  
  
  #Catch by Depth by zone by species
  cfac=function(x,breaks=seq(0,200,25))
  {
    #if(is.null(breaks)) breaks=unique(quantile(x,probs = seq(0, 1, 0.2)))
    x=cut(x,breaks,include.lowest=T,right=F)
    levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                                                   c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
    return(x)
  }
  
  ZONES=c("West","Zone1","Zone2")
  fn.Catch.Depth=function(SP)
  {
    DATA=subset(Fine.scale.eff,species==SP & depthMax>0)
    DATA$LatDeg=-DATA$LatDeg
    DATA$zone=with(DATA,ifelse(LongDeg>=116.5 & LatDeg<=(-26),"Zone2",
                               ifelse(LongDeg<116.5&LatDeg<=(-33),"Zone1",
                                      ifelse(LatDeg>(-33) & LatDeg<=(-26) & LongDeg<116.5,"West",NA))))
    
    for (i in 1:length(ZONES))
    {
      DATA1=subset(DATA,zone==ZONES[i] & livewt>0)
      DATA1$Depth.Range=cfac(DATA1$depthMax)
      if(i==2)barplot(t(with(DATA1,table(blockx,Depth.Range))),beside=T,main=ZONES[i],legend.text=T,col=c(1:8))
      if(!(i==2))barplot(t(with(DATA1,table(blockx,Depth.Range))),beside=T,main=ZONES[i],col=c(1:8))
      box()
    }
  }
  
  for (ss in 1:length(TARGETS))
  {
    tiff(file=paste("Figure.",TARGETS[ss],".Catch_depth.fine.scale.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    
    par(mfrow=c(3,1),mai = c(0.5, 0.5, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(2, 0.5, 0),las=1)
    
    fn.Catch.Depth(TARGETS[ss])
    mtext("Block",side=1,line=-1,outer=T)
    mtext("Number of records",side=2,line=-1,outer=T,las=3)
    dev.off()
  }
  
}




########### SECTION J. ----  DROPPED CODE --- ###########

#3.1. Number of blocks and vessels per Yr.Mn
# Expand.fun.Yr.Mn=function(DATA)
# {
#   DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   Tab=table(DATA$Yr.Mn,DATA$BLOCKX)
#   Tab=ifelse(Tab>=1,1,0)
#   Tab=rowSums(Tab)
#   Yr.Mn=names(Tab)
#   Yr=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 1))
#   Mn=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 2))
#   Tab=data.frame(Tab,Yr,Mn)
#   Tab=Tab[order(Tab$Yr,Tab$Mn),1]
# 
#   return(Tab)
# }
# Spatial.expan.Yr.Mn=Expand.fun.Yr.Mn(Data.monthly)
# 
# Effort1.fun.YrMn=function(DATA)
# {
#   DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   Tab=table(DATA$Yr.Mn,DATA$VESSEL)
#   Tab=ifelse(Tab>=1,1,0)
#   Tab=rowSums(Tab)
#   Yr.Mn=names(Tab)
#   Yr=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 1))
#   Mn=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 2))
#   Tab=data.frame(Tab,Yr,Mn)
#   Tab=Tab[order(Tab$Yr,Tab$Mn),1]
#   
#   return(Tab)
# }
# Effort.expan=Effort1.fun.YrMn(Data.monthly)
# 
# tiff(file="Figure 3. Folly and Fantasy.1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# par(mfcol=c(1,1),mar=c(1,3.6,.1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.15,.7,0))
# 
# plot(1:N.Yrs.months,Spatial.expan.Yr.Mn,ylab="",xlab="",xaxt="n",las=2,type="l",pch=19,cex=1.5,ylim=c(0,50)
#      ,cex.axis=1.25,lwd=1.75)
# axis(1,at=seq(1,N.Yrs.months,12),labels=F,tck=-0.01)
# 
# par(new=T)
# plot(1:N.Yrs.months,Effort.expan,col="grey60",type="l",pch=19,axes=F,ann='F',cex=1.5,lwd=1.75,ylim=c(0,105))
# axis(4,at=pretty(Effort.expan),labels=pretty(Effort.expan),las=2,cex.axis=1.3)
# 
# 
# axis(1,at=seq(1,N.Yrs.months,24),labels=names(Yrs.months)[seq(1,N.Yrs.months,24)],tck=-0.02,cex.axis=1.25)
# 
# 
# mtext("Year",side=1,line=1.5,font=1,las=0,cex=1.7,outer=T)
# mtext("Number of blocks fished",side=2,line=-1.2,font=1,las=0,cex=1.7,outer=T)
# mtext("Number of licence holders fishing",side=4,line=-1.2,las=3,cex=1.7,col="grey60",outer=T)
# dev.off()



# Mean.fun.Yr.Mn=function(VAR)
# {
#   DATA=Data.monthly.GN
#   
#   DATA$CPUE=with(DATA,ifelse(BDAYS.c*HOURS.c*SHOTS.c>0,LIVEWT.c/((NETLEN.c/1000)*BDAYS.c*HOURS.c*SHOTS.c*Inc.per),NA))
#   
#   #Yr.Mn
#   #DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   #Tab=aggregate(CPUE~Yr.Mn+BLOCKX+SPECIES,data=DATA,mean,na.rm=T)
#   
#   #Yr only
#   Tab=aggregate(CPUE~YEAR.c+BLOCKX+SPECIES,data=DATA,mean,na.rm=T)
#   
#   
#   Tab1=subset(Tab,SPECIES==VAR)
# #  Tab1$Yr=as.numeric(sapply(strsplit(Tab1$Yr.Mn," "), "[", 1))
# #  Tab1$Mn=as.numeric(sapply(strsplit(Tab1$Yr.Mn," "), "[", 2))
# #  Tab1=Tab1[order(Tab1$Yr,Tab1$Mn),]
# #  Tab1=Tab1[,c(2,4:6)]
# 
# #  Reshaped=as.matrix(reshape(Tab1,idvar=c("Yr","Mn"),  	#transposed as matrix 	
# #                             timevar="BLOCKX",v.names="CPUE", direction="wide"))	
#   
#   Tab1=Tab1[,c(1:2,4)]
#   
#   Reshaped=as.matrix(reshape(Tab1,idvar=c("YEAR.c"),    #transposed as matrix 	
#                              timevar="BLOCKX",v.names="CPUE", direction="wide"))
#   Reshaped=Reshaped[order(Reshaped[,1]),]
#   return(Reshaped)
# }
#  Spatial.expan.Yr.Mn=Mean.fun.Yr.Mn(TARGETS[2])
# # 
# # tiff(file="Figure 3. Folly and Fantasy.3.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
# # 
# # plot(1:N.Yrs.months,Spatial.expan.Yr.Mn[,3],ylim=c(0,10))
# # for(i in 4:ncol(Spatial.expan.Yr.Mn))   lines(1:N.Yrs.months,Spatial.expan.Yr.Mn[,i],col=runif(1,1,100))
# 
#  plot(1:nrow(Spatial.expan.Yr.Mn),Spatial.expan.Yr.Mn[,2],ylim=c(0,20),ann=F,xaxt='n',col='transparent')
#  for(i in 2:ncol(Spatial.expan.Yr.Mn))   lines(1:nrow(Spatial.expan.Yr.Mn),Spatial.expan.Yr.Mn[,i],col=runif(1,1,100))
# 

#dev.off()


# fun.prop=function(DAT,SPEC)
# {
#   #Vessel, gear, fin. year, month, block (given by the "Same.return" variable)
#   ID=which(DAT$SPECIES==SPEC)
#   this.same.returns=unique(DAT[ID,]$Same.return)
#   dat=subset(DAT,Same.return %in% this.same.returns & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~Same.return,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~Same.return,data=dat,sum,na.rm=T)
#   Prop.VesYrMonBlock=data.frame(Same.return=All.1$Same.return,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #   #Vessel, fin. year (mean proportion given by the "AnnualVesselAveID" variable)
#   #   this.same=unique(DAT[ID,]$AnnualVesselAveID)
#   #   dat=subset(DAT,AnnualVesselAveID %in% this.same & SPECIES %in% Shark.species)
#   #   dat.species=subset(dat,SPECIES==SPEC)
#   #   Target.sp.1=aggregate(LIVEWT~AnnualVesselAveID,data=dat.species,sum,na.rm=T)
#   #   All.1=aggregate(LIVEWT~AnnualVesselAveID,data=dat,sum,na.rm=T)
#   #   Prop.VesFinYr=data.frame(AnnualVesselAveID=All.1$AnnualVesselAveID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, month (mean proportion given by the "MonthlyID" variable)
#   this.same=unique(DAT[ID,]$MonthlyID)
#   dat=subset(DAT,MonthlyID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~MonthlyID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~MonthlyID,data=dat,sum,na.rm=T)
#   Prop.FinYrMon=data.frame(MonthlyID=All.1$MonthlyID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #   #Fin. year, block (mean proportion given by the "BlockID" variable) 
#   #   this.same=unique(DAT[ID,]$BlockID)
#   #   dat=subset(DAT,BlockID %in% this.same & SPECIES %in% Shark.species)
#   #   dat.species=subset(dat,SPECIES==SPEC)
#   #   Target.sp.1=aggregate(LIVEWT~BlockID,data=dat.species,sum,na.rm=T)
#   #   All.1=aggregate(LIVEWT~BlockID,data=dat,sum,na.rm=T)
#   #   Prop.FinYrBlok=data.frame(BlockID=All.1$BlockID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, month, block (mean proportion given by the "GoodsplitID" variable)
#   this.same=unique(DAT[ID,]$GoodsplitID)
#   dat=subset(DAT,GoodsplitID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~GoodsplitID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~GoodsplitID,data=dat,sum,na.rm=T)
#   Prop.GoodsplitID=data.frame(GoodsplitID=All.1$GoodsplitID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, zone (given by the "ZoneID" variable)  (mean proportion)
#   this.same=unique(DAT[ID,]$ZoneID)
#   dat=subset(DAT,ZoneID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~ZoneID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~ZoneID,data=dat,sum,na.rm=T)
#   Prop.FinYrZone=data.frame(ZoneID=All.1$ZoneID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #zone (given by the "zone" variable)  (mean proportion)
#   this.same=unique(DAT[ID,]$zone)
#   dat=subset(DAT,zone %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~zone,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~zone,data=dat,sum,na.rm=T)
#   Prop.Zone=data.frame(zone=All.1$zone,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   return(list(Prop.VesYrMonBlock=Prop.VesYrMonBlock,Prop.GoodsplitID=Prop.GoodsplitID,
#               Prop.FinYrZone=Prop.FinYrZone,Prop.FinYrMon=Prop.FinYrMon,Prop.Zone=Prop.Zone))
# }

# Catch.prop.gummy=fun.prop(Data.monthly,17001)
# Catch.prop.whiskery=fun.prop(Data.monthly,17003)
# Catch.prop.dusky=fun.prop(Data.monthly,18003)
# Catch.prop.sandbar=fun.prop(Data.monthly,18007)
# Catch.prop.school=fun.prop(Data.monthly,17008)
# Catch.prop.dogfish=fun.prop(Data.monthly,20000)
# Catch.prop.other=fun.prop(Data.monthly,Sharks.other)


# #create bad reporter files for fixing catches
# Bad.Reporters=subset(Data.monthly,Reporter=="bad")
# 
# Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter))
# 
# Bad.dus.gum.whi=subset(Bad.Reporters,Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0)
# NroW=nrow(Bad.dus.gum.whi)
# 
# Bad.dus.gum.whi.noBMY=subset(Bad.Reporters,!(Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0))
# 
# 
# #Replicate Bad.dus.gum.whi twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi=rbind(Bad.dus.gum.whi,Bad.dus.gum.whi,Bad.dus.gum.whi)
# Bad.dus.gum.whi=Bad.dus.gum.whi[order(Bad.dus.gum.whi$Same.return),]
# 
# Bad.dus.gum.whi$Spec.old=Bad.dus.gum.whi$SPECIES
# Bad.dus.gum.whi$Sname.old=Bad.dus.gum.whi$SNAME
# 
# Bad.dus.gum.whi$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.dus.gum.whi$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
# Bad.dus.gum.whi$LIVEWT.reap=with(Bad.dus.gum.whi,
#                                  ifelse(SPECIES%in%c(18003),Shark.other.livewt*Prop.Dus.Good.spl,
#                                         ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Good.spl,
#                                                ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi=subset(Bad.dus.gum.whi,LIVEWT.reap>0)
# 
# #create new vars
# Bad.dus.gum.whi$Reporter.old=Bad.dus.gum.whi$Reporter
# Bad.dus.gum.whi$Reporter="good"
# 
# 
# #add old species column to data
# Data.monthly$Spec.old=Data.monthly$SPECIES
# Data.monthly$Sname.old=Data.monthly$SNAME
# Data.monthly$Reporter.old=Data.monthly$Reporter
# 
# 
# #update "bad" recorders with reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.8  Update good split ave zone catches
# #If valid month-year-block proportions NOT available, then update "bad" records of Dusky, Gummy and
# #     whiskery with year-month-zone average
# 
# Data.monthly=merge(Data.monthly,Zone.good.split,by="ZoneID",all.x=T)
# Bad.dus.gum.whi.noBMY=merge(Bad.dus.gum.whi.noBMY,Zone.good.split,by="ZoneID",all.x=T)
# 
# Bad.dus.gum.whi.noBMY.month=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999 &
#                                      !(Prop.Dus.Zone.Good.spl>0 & Prop.Gum.Zone.Good.spl>0 & Prop.Whi.Zone.Good.spl>0))
# 
# 
# NroW.noBMY=nrow(Bad.dus.gum.whi.noBMY)
# 
# #Replicated Bad.dus.gum.whi.noBMY twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi.noBMY=rbind(Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY)
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[order(Bad.dus.gum.whi.noBMY$Same.return),]
# 
# Bad.dus.gum.whi.noBMY$Spec.old=Bad.dus.gum.whi.noBMY$SPECIES
# Bad.dus.gum.whi.noBMY$Sname.old=Bad.dus.gum.whi.noBMY$SNAME
# 
# 
# Bad.dus.gum.whi.noBMY$SPECIES=rep(c(18003,17001,17003),NroW.noBMY)
# Bad.dus.gum.whi.noBMY$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY)
# 
# 
# Bad.dus.gum.whi.noBMY$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY,
#                                        ifelse(SPECIES%in%c(18003),Shark.other.livewt*Prop.Dus.Zone.Good.spl,
#                                               ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Zone.Good.spl,
#                                                      ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Zone.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi.noBMY=subset(Bad.dus.gum.whi.noBMY,LIVEWT.reap>0)
# 
# Bad.dus.gum.whi.noBMY$Reporter.old=Bad.dus.gum.whi.noBMY$Reporter
# Bad.dus.gum.whi.noBMY$Reporter="good"
# 
# 
# #update "bad" recorders with reapportioned catch
# ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY))
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[,ID.names]
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.9  Update good split ave monthly catches           #Rory's rules 4f                          
# #If valid month-year-block proportions NOT available or month-year-block proportions NOT available,
# #       then update "bad" records of Dusky, Gummy and whiskery with month-year average
# 
# Data.monthly=merge(Data.monthly,Monthly.good.split,by="MonthlyID",all.x=T)
# 
# if(nrow(Bad.dus.gum.whi.noBMY.month)>0)
# {
#   Bad.dus.gum.whi.noBMY.month=merge(Bad.dus.gum.whi.noBMY.month,Monthly.good.split,by="MonthlyID",all.x=T)
#   
#   NroW.noBMY.month=nrow(Bad.dus.gum.whi.noBMY.month)
#   
#   #Replicated Bad.dus.gum.whi.noBMY.month twice to add catch of dusky, gummy and whiskery
#   Bad.dus.gum.whi.noBMY.month=rbind(Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[order(Bad.dus.gum.whi.noBMY.month$Same.return),]
#   
#   Bad.dus.gum.whi.noBMY.month$Spec.old=Bad.dus.gum.whi.noBMY.month$SPECIES
#   Bad.dus.gum.whi.noBMY.month$Sname.old=Bad.dus.gum.whi.noBMY.month$SNAME
#   
#   
#   Bad.dus.gum.whi.noBMY.month$SPECIES=rep(c(18003,17001,17003),NroW.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY.month)
#   
#   Bad.dus.gum.whi.noBMY.month$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY.month,
#                                                ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Mon.Good.spl,
#                                                       ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Mon.Good.spl,
#                                                              ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Mon.Good.spl,LIVEWT.reap))))
#   
#   #remove artificially created 0 catches
#   Bad.dus.gum.whi.noBMY.month=subset(Bad.dus.gum.whi.noBMY.month,LIVEWT.reap>0)
#   
#   Bad.dus.gum.whi.noBMY.month$Reporter.old=Bad.dus.gum.whi.noBMY.month$Reporter
#   Bad.dus.gum.whi.noBMY.month$Reporter="good"
#   
#   
#   #update "bad" recorders with reapportioned catch
#   ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY.month))
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[,ID.names]
#   Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY.month)
#   
#   #remove duplicates
#   Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
#   Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
#   
# }
# 
# # 
# #C.7.15 Reapportion catch                                     #Rory's rules 6k-6s            
# #note: uses same rules as for southern catch (#C7.7- #C7.9)
# 
# #create bad reporter files for fixing catches
# Bad.Reporters=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999)
# Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter)| !(Reporter=="bad" & SPECIES== 22999))
# 
# 
# Bad.dus.gum.whi=subset(Bad.Reporters,Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0)
# NroW=nrow(Bad.dus.gum.whi)
# 
# Bad.dus.gum.whi.noBMY=subset(Bad.Reporters,!(Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0))
# 
# #C.7.15.1 Good.spl criteria
# #Replicated Bad.dus.gum.whi twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi=rbind(Bad.dus.gum.whi,Bad.dus.gum.whi,Bad.dus.gum.whi)
# Bad.dus.gum.whi=Bad.dus.gum.whi[order(Bad.dus.gum.whi$Same.return),]
# 
# Bad.dus.gum.whi$Spec.old=Bad.dus.gum.whi$SPECIES
# Bad.dus.gum.whi$Sname.old=Bad.dus.gum.whi$SNAME
# 
# Bad.dus.gum.whi$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.dus.gum.whi$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
# 
# Bad.dus.gum.whi$LIVEWT.reap=with(Bad.dus.gum.whi,
#                                  ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Good.spl,
#                                         ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Good.spl,
#                                                ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi=subset(Bad.dus.gum.whi,LIVEWT.reap>0)
# 
# #create new vars
# Bad.dus.gum.whi$Reporter.old=Bad.dus.gum.whi$Reporter
# Bad.dus.gum.whi$Reporter="good"
# 
# #update "bad" recorders with reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.15.2 ZoneID criteria
# NroW.noBMY=nrow(Bad.dus.gum.whi.noBMY)
# 
# #Replicated Bad.dus.gum.whi.noBMY twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi.noBMY=rbind(Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY)
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[order(Bad.dus.gum.whi.noBMY$Same.return),]
# 
# 
# Bad.dus.gum.whi.noBMY$SPECIES=rep(c(18003,17001,17003),NroW.noBMY)
# Bad.dus.gum.whi.noBMY$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY)
# 
# 
# Bad.dus.gum.whi.noBMY$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY,
#                                        ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Zone.Good.spl,
#                                               ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Zone.Good.spl,
#                                                      ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Zone.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi.noBMY=subset(Bad.dus.gum.whi.noBMY,LIVEWT.reap>0)
# 
# Bad.dus.gum.whi.noBMY$Reporter.old=Bad.dus.gum.whi.noBMY$Reporter
# Bad.dus.gum.whi.noBMY$Reporter="good"
# 
# 
# #update "bad" recorders with reapportioned catch
# ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY))
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[,ID.names]
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.15.3 Ave monthly catches criteria                               
# 
# Bad.dus.gum.whi.noBMY.month=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999)
# 
# if(nrow(Bad.dus.gum.whi.noBMY.month)>0)
# {
#   NroW.noBMY.month=nrow(Bad.dus.gum.whi.noBMY.month)
#   
#   #Replicated Bad.dus.gum.whi.noBMY.month twice to add catch of dusky, gummy and whiskery
#   Bad.dus.gum.whi.noBMY.month=rbind(Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[order(Bad.dus.gum.whi.noBMY.month$Same.return),]
#   
#   
#   Bad.dus.gum.whi.noBMY.month$SPECIES=rep(c(18003,17001,17003),NroW.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY.month)
#   
#   Bad.dus.gum.whi.noBMY.month$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY.month,
#                                                ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Mon.Good.spl,
#                                                       ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Mon.Good.spl,
#                                                              ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Mon.Good.spl,LIVEWT.reap))))
#   
#   #remove artificially created 0 catches
#   Bad.dus.gum.whi.noBMY.month=subset(Bad.dus.gum.whi.noBMY.month,LIVEWT.reap>0)
#   
#   Bad.dus.gum.whi.noBMY.month$Reporter.old=Bad.dus.gum.whi.noBMY.month$Reporter
#   Bad.dus.gum.whi.noBMY.month$Reporter="good"
#   
#   
#   #update "bad" recorders with reapportioned catch
#   ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY.month))
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[,ID.names]
#   Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY.month)
#   
#   #remove duplicates
#   Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
#   Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
#   
# }
# 



# #Compare mine and Rory's catch and effort
# setwd("C:/Users/myb/Desktop/New folder")
# 
# Use.Previos.Sofar="YES"   #Select YES if attaching previous Sofar data to current year
# #Use.Previos.Sofar="NO"
# 
# 
# Ind.spe.list=list(Gummy=17001,Whiskery=17003,Bronzy.Dusky=c(18001,18003),sandbar=18007)
# 
# #Fishing effort limits
# FishEffLims=data.frame(zone=c("West","Zone1","Zone2"),Km.Gillnet.Hours.c=c(67692,84075,144102),
#                        Km.Gillnet.Days.c=c(2832,3503,7205))
# 
# #Current year data set
# DAT=subset(Data.monthly,FINYEAR==Current.yr)
# 
# #Main Feature table (total catch, indicator species catch, teleost catch, etc)
# Other=subset(Other.fishery.catch,financial.year==Current.yr)        
# 
# 
# #Add Fishing effort                                       #REVIEW RORY
# 
# C.yr=match(Current.yr,Total.effort.days.monthly$FINYEAR)
# 
# #annual
# 
# Curr.annual.1000km.gn.hours=Total.effort.hours.monthly[C.yr,]
# Curr.annual.km.gn.days=Total.effort.days.monthly[C.yr,]
# Curr.annual.km.gn.days[1,2]=Curr.annual.km.gn.days[1,2]*1000
# 
# #annual by zone
# Curr.annual.1000km.gn.hours.zone=Total.effort.zone.hours.monthly[C.yr,]
# Curr.annual.km.gn.days.zone=Total.effort.zone.days.monthly[C.yr,]
# Curr.annual.km.gn.days.zone[,2:4]=Curr.annual.km.gn.days.zone[,2:4]*1000
# 
# #percentages of effort limits
# #annual
# Per.annual.lim.1000km.gn.hours=100*Curr.annual.1000km.gn.hours$Total/
#   (sum(FishEffLims$Km.Gillnet.Hours.c)/1000)
# Per.annual.lim.km.gn.days=100*Curr.annual.km.gn.days$Total/sum(FishEffLims$Km.Gillnet.Days.c)
# 
# #annual by zone
# Per.annual.lim.1000km.gn.hours.zone=100*Curr.annual.1000km.gn.hours.zone[,2:4]/
#   (FishEffLims$Km.Gillnet.Hours.c/1000)
# Per.annual.lim.km.gn.days.zone=100*Curr.annual.km.gn.days.zone[,2:4]/FishEffLims$Km.Gillnet.Days.c
# 
# 
# 
# 
# #Figures 2-3. 
# par(mfcol=c(1,1),mar=c(3.5,3.6,.1,1),oma=c(1,.5,.1,.1))
# LINE=c(5,1,1)
# TYPE=c("l","o","o")
# PCH=21
# COL=1
# BG=c("black","black","white")
# 
# #Redefine DAT for "Use.Previos.Sofar=="YES"
# DAT=subset(Data.monthly,FINYEAR%in%c("2011-12",Current.yr) & METHOD%in%c("GN","LL"))
# 
# 
# 
# if(Use.Previos.Sofar=="YES")
# {
#   
#   fun.fig.all=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT[,2],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     plot(1:NN,DAT[id:N,2]/scaler,type='l',col="grey80",ylim=c(0,MAX),xaxt='n',yaxt='n',
#          ylab=TITLE1, xlab=TITLE2,las=1,lwd=2,cex.lab=1.3)
#     axis(1,at=1:NN,labels=F,tck=-0.01)
#     axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
#     
#     axis(2,at=seq(0,MAX,INT),labels=F,tck=-0.01)
#     axis(2,at=seq(0,MAX,INT2),labels=seq(0,MAX,INT2),tck=-0.02,cex.axis=1.1)
#     
#     # for(i in 1:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COL,lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   fun.fig.zn=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT1[,2:4],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     plot(1:NN,DAT1[id:N,2]/scaler,type='l',col=COL,ylim=c(0,MAX),xaxt='n',yaxt='n',
#          ylab=TITLE1, xlab=TITLE2,las=1,lwd=1,cex.lab=1.3,lty=LINE[1])
#     axis(1,at=1:NN,labels=F,tck=-0.01)
#     axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
#     
#     axis(2,at=seq(0,MAX,INT),labels=F,tck=-0.01)
#     axis(2,at=seq(0,MAX,INT2),labels=seq(0,MAX,INT2),tck=-0.02,cex.axis=1.1)
#     
#     for(i in 2:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COL,lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   fn.figs2.3.SoFaR.all=function(GROUP,LAT1,LAT2,INT,INT2,what)
#   {
#     dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     #add previous years
#     if(what=="Elasmos")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.sk.live.wt","Z1.tot.sk.live.wt","Z2.tot.sk.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.sk.live.wt"),names(Results.pre.2013))]
#     }
#     
#     if(what=="Teleosts")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.tel.live.wt","Z1.tot.tel.live.wt","Z2.tot.tel.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.tel.live.wt"),names(Results.pre.2013))]
#     }
#     
#     names(Prev.zn)=names(wide)
#     Prev.zn[,2:4]=Prev.zn[,2:4]*1000
#     wide=rbind(Prev.zn,wide)
#     
#     
#     Prev[,2]=Prev[,2]*1000
#     names(Prev)=names(annual.catch.total)
#     annual.catch.total=rbind(Prev,annual.catch.total)
#     
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.all(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
#   
#   fn.figs2.3.SoFaR.zn=function(GROUP,LAT1,LAT2,INT,INT2,what)
#   {
#     dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     #add previous years
#     if(what=="Elasmos")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.sk.live.wt","Z1.tot.sk.live.wt","Z2.tot.sk.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.sk.live.wt"),names(Results.pre.2013))]
#     }
#     
#     if(what=="Teleosts")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.tel.live.wt","Z1.tot.tel.live.wt","Z2.tot.tel.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.tel.live.wt"),names(Results.pre.2013))]
#     }
#     
#     names(Prev.zn)=names(wide)
#     Prev.zn[,2:4]=Prev.zn[,2:4]*1000
#     wide=rbind(Prev.zn,wide)
#     
#     
#     Prev[,2]=Prev[,2]*1000
#     names(Prev)=names(annual.catch.total)
#     annual.catch.total=rbind(Prev,annual.catch.total)
#     
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.zn(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
# }
# 
# if(Use.Previos.Sofar=="NO")
# {
#   fun.fig.all=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT[,2],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     lines(1:NN,DAT[id:N,2]/scaler,type='l',col="red",lwd=2)
#   }
#   
#   fun.fig.zn=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT1[,2:4],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     BG=c("red","green","blue") 
#     COLs=BG
#     for(i in 1:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COLs[i],lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   
#   fn.figs2.3.SoFaR.all=function(GROUP,LAT1,LAT2,INT,INT2)
#   {
#     dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2 & METHOD%in%c("GN","LL"))
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     names(wide)[match("FINYEAR",names(wide))]="finyear"
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.all(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
#   
#   fn.figs2.3.SoFaR.zn=function(GROUP,LAT1,LAT2,INT,INT2)
#   {
#     dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     names(wide)[match("FINYEAR",names(wide))]="finyear"
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.zn(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
# }
# 
# #total
# jpeg(file="Figure 2.TotalElasmoCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fn.figs2.3.SoFaR.all(Elasmo.species,-27,-40,100,500,"Elasmos")
# if(Use.Previos.Sofar=="NO")fn.figs2.3.SoFaR.all(Elasmo.species,-27,-40,100,500)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("grey80","red"),lwd=2)
# dev.off()
# 
# #by zone
# jpeg(file="Figure 2.ZoneElasmoCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fn.figs2.3.SoFaR.zn(Elasmo.species,-27,-40,100,500,"Elasmos")
# if(Use.Previos.Sofar=="NO")fn.figs2.3.SoFaR.zn(Elasmo.species,-27,-40,100,500)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("black","red"),lwd=2)
# legend("topleft",c("West","Zn1","Zn2"),bty='n',lty=LINE,col=c("red","green","blue"),lwd=2)
# dev.off()
# 
# 
# 
# 
# #Figure 4
# 
# #add 2011-12 to Current
# Current.yr=c("2011-12","2012-13")
# Eff.Current.Yr=fn.Eff.Sofar(Current.yr)
# Total.effort.zone.days=Eff.Current.Yr$Total.effort.zone.days
# Total.effort.zone.hours=Eff.Current.Yr$Total.effort.zone.hours
# Total.effort.joint.days=Eff.Current.Yr$Total.effort.joint.days
# Total.effort.joint.hours=Eff.Current.Yr$Total.effort.joint.hours
# Total.effort.days=Eff.Current.Yr$Total.effort.days
# Total.effort.hours=Eff.Current.Yr$Total.effort.hours
# 
# if(Use.Previos.Sofar=="YES")
# {
#   Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.km.gn.days","Z1.km.gn.days","Z2.km.gn.days"),
#                                   names(Results.pre.2013))]
#   Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.km.gn.days"),names(Results.pre.2013))]
#   
#   idi=match(Current.yr,Total.effort.zone.days$FINYEAR)
#   wide=data.frame(finyear=Current.yr,WC.km.gn.days=Total.effort.zone.days[idi,2]*1000,
#                   Z1.km.gn.days=Total.effort.zone.days[idi,3]*1000,
#                   Z2.km.gn.days=Total.effort.zone.days[idi,4]*1000)
#   
#   names(Prev.zn)=names(wide)
#   wide=rbind(Prev.zn,wide)
#   
#   annual.effort.days.total=data.frame(finyear=Current.yr,TDGDLF.km.gn.days=Total.effort.days[idi,2]*1000)
#   names(Prev)=names(annual.effort.days.total)
#   annual.effort.days.total=rbind(Prev,annual.effort.days.total)
#   
#   
# }
# 
# if(Use.Previos.Sofar=="NO")
# {
#   annual.effort.days.total=Total.effort.days.monthly
#   wide=Total.effort.zone.days.monthly
#   names(wide)[match("FINYEAR",names(wide))]="finyear"
#   names(annual.effort.days.total)[match("FINYEAR",names(annual.effort.days.total))]="finyear"
# }
# 
# #all
# jpeg(file="Figure 4.StandardisedEffort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fun.fig.all(annual.effort.days.total,wide,1000,"Effort (1000km gn.d)","Financial year",2,10)
# if(Use.Previos.Sofar=="NO")fun.fig.all(annual.effort.days.total,wide,1,"Effort (1000km gn.d)","Financial year",2,10)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("grey80","red"),lwd=2)
# dev.off()
# 
# #by zone
# jpeg(file="Figure 4.ZoneStandardisedEffort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fun.fig.zn(annual.effort.days.total,wide,1000,"Effort (1000km gn.d)","Financial year",2,10)
# if(Use.Previos.Sofar=="NO")fun.fig.zn(annual.effort.days.total,wide,1,"Effort (1000km gn.d)","Financial year",2,10)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("black","red"),lwd=2)
# legend("topleft",c("West","Zn1","Zn2"),bty='n',lty=LINE,col=c("red","green","blue"),lwd=2)
# 
# dev.off()




#REMOVED REAPPORTIONING CODE
# 
#   #C.7.8.2 First fix Bad.dusky, Bad.gummy, Bad.whiskery                         #REVIEW RORY
# Bad.Reporters.Dus.Gum.Whi=subset(Bad.Reporters,SPECIES%in%c(18003,17001,17003))
# 
#     #reapportion catch
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=NA
# 
#         #first use "Good.spl" 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#         ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#         ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#         ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#         LIVEWT.reap))))
# 
#         #second, if previous not available, use "Zone.Good.spl" (i.e. Yr-Mn-Zone)
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#     ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #USE YR-MONTH AGAIN!!!!
#       #third, if previous are not available, use "Mon.Good.spl"  (i.e. Yr-Mn)            #Rory's rules 4f 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
#         (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
#          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
#         (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
#     LIVEWT.reap))))
# 
# #         #third, if previous are not available, use "YrZn.Good.spl" (i.e. Yr-Zone)         #REVIEW RORY 
# # Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
# #      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
# #     LIVEWT.reap))))
# 
# #create new vars
# Bad.Reporters.Dus.Gum.Whi$Reporter.old=Bad.Reporters.Dus.Gum.Whi$Reporter
# Bad.Reporters.Dus.Gum.Whi$Reporter="good"
# Bad.Reporters.Dus.Gum.Whi$Spec.old=Bad.Reporters.Dus.Gum.Whi$SPECIES
# Bad.Reporters.Dus.Gum.Whi$Sname.old=Bad.Reporters.Dus.Gum.Whi$SNAME
# 
# 
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters.Dus.Gum.Whi)
# 
# 
#   #C.7.8.3 Then fix 22999                               #REVIEW RORY
#     #note: remove duplicates of Same.return as Tot.shk.livewt is split proportionally
# #          among dusky, whiskery and gummy
# Bad.Reporters=subset(Bad.Reporters,SPECIES==22999)
# Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),]
# NroW=nrow(Bad.Reporters)
# 
#     # replicate Bad.Reporters twice to have the three species as a record
# Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters)
# Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]
# 
# Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
# Bad.Reporters$Sname.old=Bad.Reporters$SNAME
# 
# Bad.Reporters$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
#   # reapportion catch
# Bad.Reporters$LIVEWT.reap=NA
# 
# #first use "Good.spl" (standardise the proportions to sum(split catch)=Tot.shk.livewt)
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#      ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#      ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#      LIVEWT.reap))))
# 
# #Second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#     LIVEWT.reap))))
# 
# # #Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# # Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
# #     ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #         (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #     LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #create file for flagging bad reporters
# #Flag.bad.rep1=Bad.Reporters[,match(c("Same.return","Spec.old","Reporter"),names(Bad.Reporters))]
# 
# #create new vars
# Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
# Bad.Reporters$Reporter="good"
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters)

#   #C.7.16.2 First fix Bad.dusky, Bad.gummy, Bad.whiskery 
# Bad.Reporters.Dus.Gum.Whi=subset(Bad.Reporters,SPECIES%in%c(18003,17001,17003))
# 
# #reapportion catch
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=NA
# 
# #first use "Good.spl" 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#            ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,
#            (Tot.shk.livewt*Prop.Dus.Good.spl),
#            ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,
#            (Tot.shk.livewt*Prop.Gum.Good.spl),
#            ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,
#            (Tot.shk.livewt*Prop.Whi.Good.spl),
#             LIVEWT.reap))))
# 
# #second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#             ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,
#             (Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#             ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,
#            (Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#            ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,
#           (Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#           LIVEWT.reap))))
# 
# # #third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# # Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
# #                  ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #                 (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #                 ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #                (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #                ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #               (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #             LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# 
# #create new vars
# Bad.Reporters.Dus.Gum.Whi$Reporter.old=Bad.Reporters.Dus.Gum.Whi$Reporter
# Bad.Reporters.Dus.Gum.Whi$Reporter="good"
# Bad.Reporters.Dus.Gum.Whi$Spec.old=Bad.Reporters.Dus.Gum.Whi$SPECIES
# Bad.Reporters.Dus.Gum.Whi$Sname.old=Bad.Reporters.Dus.Gum.Whi$SNAME
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters.Dus.Gum.Whi)
# 
# 
# 
#   #C.7.16.3 Then fix 22999
#     # remove duplicates of Same.return as Tot.shk.livewt is split proportionally
# Bad.Reporters=subset(Bad.Reporters,SPECIES==22999)
# Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),]
# NroW=nrow(Bad.Reporters)
# 
#     # replicate Bad.Reporters twice to have the three species as a record
# Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters)
# Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]
# 
# Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
# Bad.Reporters$Sname.old=Bad.Reporters$SNAME
# 
# Bad.Reporters$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
#     # reapportion catch
# Bad.Reporters$LIVEWT.reap=NA
# 
# #first use "Good.spl" 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#      ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#      ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#      LIVEWT.reap))))
# 
# 
# #Second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,
#         (Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#      ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,
#          (Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#      ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,
#          (Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#      LIVEWT.reap))))
# 
# #Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# #note: this is nonsense because it aggregates across all zones when gummy and whiskery
# #       do not occur in the north so it's not applied
# # Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
# #      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #           (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #         (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #     LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #remove artificially created 0 catches
# #Bad.Reporters=subset(Bad.Reporters,!is.na(LIVEWT.reap))
# 
# #create file for flagging bad reporters
# #Flag.bad.rep3=Bad.Reporters[,match(c("Same.return","Spec.old","Reporter"),names(Bad.Reporters))]
# 
# #create new vars
# Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
# Bad.Reporters$Reporter="good"
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters)



# #4.9 Create data for tracking mean weight
# ALL.THIS=c(names(Data.daily.agg.Numbers)[-18],"Km.Gillnet.Days.c","Km.Gillnet.Hours.c","LIVEWT.c","NETLEN.c")
# Mean.w.whiskery=Data.monthly.GN.whiskery[,match(ALL.THIS,names(Data.monthly.GN.whiskery))]
# Mean.w.gummy=Data.monthly.GN.gummy[,match(ALL.THIS,names(Data.monthly.GN.gummy))]
# Mean.w.dusky=Data.monthly.GN.dusky[,match(ALL.THIS,names(Data.monthly.GN.dusky))]
# Mean.w.sandbar=Data.monthly.GN.sandbar[,match(ALL.THIS,names(Data.monthly.GN.sandbar))]
# 
# This.yrs.weight=sort(unique(Data.daily.agg.Numbers$FINYEAR))
# 
# Mean.w.whiskery=subset(Mean.w.whiskery,FINYEAR%in%This.yrs.weight)
# Mean.w.gummy=subset(Mean.w.gummy,FINYEAR%in%This.yrs.weight)
# Mean.w.dusky=subset(Mean.w.dusky,FINYEAR%in%This.yrs.weight)
# Mean.w.sandbar=subset(Mean.w.sandbar,FINYEAR%in%This.yrs.weight)
# 
#   
# THIS.N=subset(Data.daily.agg.Numbers,METHOD=="GN",select=c(Same.return,nfish,SPECIES))
# Mean.w.whiskery=merge(Mean.w.whiskery,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.gummy=merge(Mean.w.gummy,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.dusky=merge(Mean.w.dusky,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.sandbar=merge(Mean.w.sandbar,THIS.N,by=c("Same.return","SPECIES"),all.x=T)


#   #Mean weigth data
# write.csv(Mean.w.whiskery,file ="Mean.w.whiskery.GN.csv")
# write.csv(Mean.w.gummy,file ="Mean.w.gummy.GN.csv")
# write.csv(Mean.w.dusky,file ="Mean.w.dusky.GN.csv")
# write.csv(Mean.w.sandbar,file ="Mean.w.sandbar.GN.csv")


# #Rescale catches again
# #note: recalculate the reapportioned catches of dusky, whiskery, gummy and 'other'
# #     considering the real catches of other shark species
# a=unique(Bad.Reporters$Same.return)
# b=subset(Data.monthly,Same.return%in%a & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26))
# 
# if(Reapportion.daily=="YES") 
# {
#   a.daily=unique(Bad.Reporters.daily$Same.return)
#   b.daily=subset(Data.daily,Same.return%in%a.daily & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26))  
# }
# 
# xx=function(dat)
# {
#   x=NA
#   if(round(sum(dat$LIVEWT.reap,na.rm=T))>round(unique(dat$Tot.shk.livewt))) x=unique(dat$Same.return)
#   return(x)
# }
# 
# VECT=rep(NA,length(a))
# for ( i in 1:length(a))VECT[i]=xx(subset(b,Same.return==a[i]))
# vect=VECT[!is.na(VECT)]   #all returns where reapportioned weight > total shark weight
# 
# if(Reapportion.daily=="YES")
# {
#   VECT.daily=rep(NA,length(a.daily))
#   for ( i in 1:length(a.daily))VECT.daily[i]=xx(subset(b.daily,Same.return==a.daily[i]))
#   vect.daily=VECT.daily[!is.na(VECT.daily)]   #all returns where reapportioned weight > total shark weight  
# }
# 
#   #monthly
# STOREss=vector("list",length(vect))
# for (i in 1:length(vect))
# {
#   s=subset(b,Same.return==vect[i])
#   Non.Fixed.shk.wgt=unique(s$Tot.shk.livewt)-sum(s$LIVEWT[which(!s$SPECIES%in%Fix.species)])
#   s$Non.Fixed.shk.wgt=Non.Fixed.shk.wgt  
#   Tot.reap=sum(s$LIVEWT.reap[which(s$SPECIES%in%Fix.species)])  
#   s$LIVEWT.reap=with(s,ifelse(SPECIES%in%Fix.species,
#                               Non.Fixed.shk.wgt*(LIVEWT.reap/Tot.reap),LIVEWT.reap))   
#   STOREss[[i]]=s
# }
# bb=do.call(rbind,STOREss)
# bb=bb[,-match("Non.Fixed.shk.wgt",names(bb))]
# Data.monthly=subset(Data.monthly,!(Same.return%in%vect & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26)))
# Data.monthly=rbind(Data.monthly,bb)
# 
# 
# #daily
# if(Reapportion.daily=="YES")
# {
#   STOREss.daily=vector("list",length(vect.daily))
#   for (i in 1:length(vect.daily))
#   {
#     s=subset(b.daily,Same.return==vect.daily[i])
#     Non.Fixed.shk.wgt=unique(s$Tot.shk.livewt)-sum(s$LIVEWT[which(!s$SPECIES%in%Fix.species)])
#     s$Non.Fixed.shk.wgt=Non.Fixed.shk.wgt  
#     Tot.reap=sum(s$LIVEWT.reap[which(s$SPECIES%in%Fix.species)])  
#     s$LIVEWT.reap=with(s,ifelse(SPECIES%in%Fix.species,
#                                 Non.Fixed.shk.wgt*(LIVEWT.reap/Tot.reap),LIVEWT.reap))   
#     STOREss.daily[[i]]=s
#   }
#   bb.daily=do.call(rbind,STOREss.daily)
#   bb.daily=bb.daily[,-match("Non.Fixed.shk.wgt",names(bb.daily))]
#   Data.daily=subset(Data.daily,!(Same.return%in%vect.daily & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26)))
#   Data.daily=rbind(Data.daily,bb.daily) 
# }

