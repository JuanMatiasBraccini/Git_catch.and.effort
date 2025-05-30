# ---------- SHARK GILLNET AND LONGLINE FISHERY CATCH AND EFFORT MANIPULATIONS ------------------------------------------------------

#NOTES:
  #NEW YEAR OF DATA:
#                     Set First.run to 'YES' and update Current.yr
#                     Bring in Monthly and Daily data thru SQL queries
#                     Also bring in TEPS and PRICE data thru SQL queries

  #Before SQL:   CAESS data from 1988/89 to 2007/08 (some fishers kept reporting in CAESS 
#                      after introduction of daily logbooks).
#                   Records include shark data from all fisheries, not just TDGDLF (fishery
#                     code: SGL & WCGL) and NSF (code C127 & CL02)
#                   Records include scalefish data from TDGDLF and NSF only
#                   Hence, combine Rory's Table 81d.mdb (1975/76 to 1987/88) with CAESS for full series

#               This is no longer relevant as SQL supersedes all these....

#                   In 2002 each unit became 270 m of net (Rory per comm)
#                   After 1988, estuary blocks are true shark shots (Rory per comm)

  #FishCUBE index
               #Extract data for FishCUBE
               #G 4.9 FishCUBE 

  #Data inspections:
                 # CATCH INSPECTIONS
                 #     Set current.yr; Identify 0 catch shots (are they real?), ID catch records
                 #         where catch is too large or too small compared to species weight range,
                 #         identify suspicious records and raise to data entry staff

                 # EFFORT INSPECTIONS
                 #     Visually inspect effort variables, identify suspicious records and raise to data
                 #         entry staff

  #Data structure
#     . Data.monthly: CAESS monthly records 1975-1976 to 2005-2006
#               each row is a species' monthly catch per vessel per block, 
#               including the effort exerted by that vessel in that block

#     . Data.daily: daily logbooks 2006-2007 onwards

  #Catch correction criteria for monthly:
#           1. Reapportion catch composition of main commercial shark species
#           2. All catch data prior to 1990 are increased by 5%

  #Effort correction criteria for monthly:
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


#check this email from VERO:
# 	There are 4 new fields in fcTourOperatorCatch, fcTourOperatorEffort, fcTourOperatorSpeciesLW:
# Your existing code should still run but depending on how your sql queries are written, 
# the output might have an additional 4 fields.  They are [FinancialYear],[Year],[Month],[YyyyMm].
# 
# 	5 additional species stats groupings have also been added:
#   [dbo].[fcTourOperatorSpeciesStatsPerFinancialYear]
# [dbo].[fcTourOperatorSpeciesStatsPerFinancialYearZone]
# [dbo].[fcTourOperatorSpeciesStatsPerYear]
# [dbo].[fcTourOperatorSpeciesStatsPerYearZone]
# [dbo].[fcTourOperatorSpeciesStatsPerZone]
# 
# They could easily be derived from the base table called [dbo].[fcTourOperatorSpeciesLW].
# But now it is even easier to get different flavours of species standard weights. More tips at the bottom of this email.
# 
# Once again, the one stop-shop for catch and effort data is:
#   Server:                 reports.dpird.wa.gov.au    # previous one was: CP-SDBS0001P-19\RESP01
# Database:                ResearchDataWarehouseQuery
# Tables:                 [dbo].[fcTourOperatorCatch]
# [dbo].[fcTourOperatorEffort]
# [dbo].[fcTourOperatorSpeciesLW], etc
# 
# FishCubeWA production links:
#   Public Cube:
#   http://f01-fims-webp01/FishCubeWA/Query.aspx?CubeId=TourOperatorPublic
# DPIRD Cube:
#   http://f01-fims-webp01/FishCubeWA/Query.aspx?CubeId=TourOperatorDPIRDOnly

# ---------- DATA SECTION ------------------------------------------------------
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
library(rgdal)
library(purrr)
library(janitor)
library(rnaturalearth)
library(ggpubr)  

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240,dplyr.summarise.inform = FALSE) 
par.default=par()
fn.user=function(x1,x2)paste(x1,Sys.getenv("USERNAME"),x2,sep='/')
if(!exists('handl_OneDrive')) source(fn.user(x1='C:/Users',
                                             x2='OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R'))

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/SoFaR.figs.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/send.emails.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))
fn.scale=function(x,scaler) ((x/max(x,na.rm=T))^0.5)*scaler

setwd(handl_OneDrive("Data/Catch and Effort"))  # working directory

First.run="NO"    
#First.run="YES"

Current.yr="2023-24"    #Set current financial year 
Current.yr.dat=paste(substr(Current.yr,1,4),substr(Current.yr,6,7),sep="_")
xx=paste(getwd(),Current.yr.dat,sep="/")
if(!file.exists(xx)) dir.create(xx) 
rm(xx)

#Email addresses
Email.data.checks="anja.giltay@dpird.wa.gov.au"
Email.FishCube="veronique.vanderklift@dpird.wa.gov.au"
Email.data.checks2="sarah.vanRyssen@dpird.wa.gov.au"

#1. Extract catch and effort data
do.sql.extraction=TRUE
  #1.1 SQL Server database
if(do.sql.extraction)
{
  #SQL connections
  use.new.server=TRUE #new server. To access previous Monthly data see '#Access Monthly previous version' 
  #use.new.server=FALSE
  
  if(use.new.server) Server <-"reports.dpird.wa.gov.au"  
  if(!use.new.server) Server <-"CP-SDBS0001P-19\\RESP01"  #original, retired 
  #conn=odbcDriverConnect(connection=paste("Driver={SQL Server};server=",Server,";database=ResearchShared;trusted_connection=yes;",sep=""))  #original, retired
  conn1=odbcDriverConnect(connection=paste("Driver={SQL Server};server=",Server,";database=ResearchDataWarehouseQuery;trusted_connection=yes;",sep=""))

  if(use.new.server) Monthly.log <- 'ceMonthlyLog'  #Use ceMonthlyLog if mismatch with monthlyrawlog is resolved. See Mark/Vero
  if(!use.new.server) Monthly.log <- 'monthlyrawlog' #original
  #notes: 1. 'monthlyrawlog' is the raw monthly data always used; 
  #         'monthlylog' is the one after the reapportioning of indicator species, etc that goes to FishCube.
  #           So I need to keep using 'monthlyrawlog' for the purpose of cpue standardisation, catch reconstructions, etc
  #           The variable SCRref identifies records that were modified, e.g. hammerheads: %>%filter(!is.na(SCRref) & grepl('Hammerhead',RSSpeciesCommonName))
  #       2. The variable VESSEL was recorded in SQL (spaces and zeros dropped) so re-map it to allow my scripts to run
  if(Monthly.log=='ceMonthlyLog')
  {
    Server.monthly <-"reports.dpird.wa.gov.au"
    Database.monthly <- "ResearchDataWarehouseQuery"
  }
  if(Monthly.log=='monthlyrawlog')
  {
    Server.monthly <-"CP-SDBS0001P-19\\RESP01"
    Database.monthly <- "ResearchDataWarehouseCE"
  }
  conn.monthly <- odbcDriverConnect(connection=paste("Driver={SQL Server};server=",Server.monthly,
                                                     ";database=",Database.monthly,
                                                     ";trusted_connection=yes;",sep=""))
  Query.monthly=paste("SELECT * FROM",Monthly.log,
                      "WHERE RSSpeciesEcologicalSuite = 'elasmobranchs' or RSInitiallyReportedFisheryCode in ('SGL','WCGL','C127','WANCS','CL02','JANS')")
  
  #Extract fishery codes table
  FisheryCodeTable=sqlFetch(channel=odbcConnectExcel2007('FisheryCodeTable.xlsx'),'CAEStoFISHCUBE')
  FisheryCodeTable <- remove_empty(FisheryCodeTable,which="cols")
  
  #Extract variable mapping
  channel <- odbcConnectExcel2007(handl_OneDrive("Data/Catch and Effort/SQL Shark Data Mapping.xlsx"))
  SQL.data.mapping_daily=sqlFetch(channel,"DailyMappings")
  SQL.data.mapping_monthly=sqlFetch(channel,"MonthlyMappings")
  close(channel)
  
  #Extract fishing methods
  FishingMethod <- sqlQuery(conn1, query="SELECT * FROM [dbo].[rsFishingMethod]")
  
  #Extract fish conditions
  fish.conditions <- sqlQuery(conn1, query="SELECT * FROM [dbo].[rsSpeciesCondition]")
  #fish.conditions <- sqlQuery(conn, query="SELECT * FROM Research.[vwSpeciesCondition]")
  
  #Extract fish conversion factors
  fish.conversion <- sqlQuery(conn1, query="SELECT * FROM [dbo].[rsSpeciesConditionFactor]")
  #fish.conversion <- sqlQuery(conn, query="SELECT * FROM Research.[vwSpeciesConditionFactor]")
  
  #Extract species mapping
  Species.mapping <- sqlQuery(conn1, query="SELECT * FROM [dbo].[rsSASSpecies]")
  
  #Extract vessel mapping
  Vessel.mapping <- sqlQuery(conn1, query="SELECT * FROM [dbo].[rsBoatLicence]")  
  
  #Extract Daily returns from TDGDLF and NSF
  system.time({Data.daily<- sqlQuery(conn1, query="SELECT * FROM ceSharkLog")})
  
    #check duplication
  Data.daily=Data.daily%>%
    mutate(Dupli=paste(DailySheetNumber,TripSheetNumber,SessionIndex,LiveWeight,TripLandedCondition,RSSpeciesCommonName,RSSpeciesCode))
  Dups=table(Data.daily$Dupli)  
  Dups=names(Dups[Dups>1])
  if(length(Dups)>0)
  {
    if(First.run=='YES')
    {
      handle.duplis=handl_OneDrive("Data/Catch and Effort")
      Duplis=Data.daily%>%
        filter(Dupli%in%Dups)%>%
        dplyr::select(VesselName,VesselRegistration,Skipper,FinancialYear,
                      Month,DailySheetNumber,TripSheetNumber,SessionIndex,RSSpeciesCommonName,LiveWeight)
      write.csv(Duplis,file=paste(handle.duplis,"/Duplis.csv",sep=""),row.names=F)
      send.email(TO=Email.data.checks,
                 CC=Email.data.checks2,
                 BCC=Email.FishCube,
                 Subject=paste("Duplicated records in daily logbooks",Sys.time(),sep=' _ '),
                 Body= "Hi,
              I've attached a spredsheet with duplicated records
              (i.e. same DailySheetNumber, SessionIndex, LiveWeight & RSSpeciesCommonName).
              Could you please have a look?
              Cheers
              Matias",  
                 Attachment=paste(handle.duplis,"/Duplis.csv",sep="")) 
    }
      
    Data.daily=Data.daily%>%
                distinct(Dupli,.keep_all = T)
  }
  rm(Dups)
  Data.daily=Data.daily%>%
    dplyr::select(-Dupli)
  
  rename.daily=SQL.data.mapping_daily$OLD
  names(rename.daily)=SQL.data.mapping_daily$NEW
  rename.daily=subset(rename.daily,!is.na(names(rename.daily)))
  rename.daily=subset(rename.daily,!is.na(rename.daily))  
  for(r in 1:length(rename.daily))
  {
    id=match(names(rename.daily)[r],colnames(Data.daily))
    colnames(Data.daily)[id]=rename.daily[r]
    rm(id)
  }
  Data.daily=Data.daily%>%
    mutate(fishery=ifelse(fishery=='SDGDL','JASDGDL',fishery),  #Changed from JASDGDL to SDGDL in 2018 but need JASDGDL to run scripts
           Lat=ifelse(is.na(Lat) & !is.na(block10),
                      as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6),
                      Lat),
           Long=ifelse(is.na(Long) & !is.na(block10),
                      100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6),
                      Long),
           Lat=abs(Lat),
           LatDeg=trunc(Lat),
           LatMin=60*(Lat%%1),
           LongDeg=trunc(Long),
           LongMin=60*(Long%%1),
           Lat=-Lat,
           sname1=RSCommonName,
           species=case_when(type=='elasmobranchs'~substr(RSSpeciesCode,4,10),
                             TRUE~substr(RSSpeciesCode,3,10)),
           species=as.numeric(species),
           species=case_when(species== 5902 ~ 5002,   #reset to previous code for consistency thru all scripts
                             species== 5901 ~ 5001,
                             species== 12001 ~ 12000,
                             species== 13900 ~ 13000,
                             species== 18901 ~ 18014,
                             species== 27909 ~ 26999,
                             species== 287003 ~ 288000,
                             species== 311010 ~ 311058,
                             species== 311912 ~ 311199,
                             species== 311991 ~ 311903,
                             species== 320004 ~ 320000,
                             species== 344004 ~ 344002,
                             species== 346014 ~ 346914,
                             species== 353013 ~ 353998,
                             species== 381010 ~ 381006,
                             species== 384039 ~ 384999,
                             species== 441911 ~ 441000,
                             species== 441006 ~ 441020,
                             species== 441912 ~ 441900,
                             species== 990009 ~ 460000,
                             species== 615000 ~ 600000,
                             species== 607901 ~ 610000,
                             species== 911005 ~ 702003,
                             species== 346015 ~ 346912,
                             species== 353000 ~ 353901,
                             species== 85018 ~ 358998,
                             species== 90001 ~ 990001,
                             species== 999997 ~ 599000,
                             species== 916901 ~ 702009,
                             species== 99993 ~ 22997,
                             species== 99994 ~ 22998,
                             species== 99992 ~ 22999,
                             species== 346916 ~ 346911,
                             species== 384902 ~ 384901,
                             species== 999997 ~ 599001,
                             species== 351000 ~ 351910,
                             TRUE~species),
           spgroup=case_when(species %in% 5001:24900 ~ 'sharks',
                             species %in% c(25000:31000,39001,990001) ~ 'rays',
                             type=='inverts'~'inverts',
                             TRUE~'scalefish'),
           zone=ifelse(fishery%in%c('NCS','WCDGDL'),'*',zone),
           bioregion=case_when(RSBioregionName=='North Coast' ~ 'NC',
                               RSBioregionName=='Gascoyne Coast' ~ 'GC',
                               RSBioregionName=='West Coast' ~ 'WC',
                               RSBioregionName=='South Coast' ~ 'SC',
                               RSBioregionName=='Outside WA State Boundaries'& Lat>(-20) ~ 'NC',
                               RSBioregionName=='Northern Inland'& Lat<(-20) & Lat>(-27)~ 'GC',
                               RSBioregionName=='Northern Inland'& Lat<(-27) & Lat>(-30)~ 'WC',
                               RSBioregionName=='Southern Inland'& Long<=115.5~ 'WC',
                               RSBioregionName=='Southern Inland'& Long>115.5~ 'SC'))  
  Keep.these.columns_daily=c("bioregion","spgroup","finyear","targetSpecies",
                             "date","DSNo","TSNo","SNo","year","month",
                             "block10","blockx","species","LatDeg","LongDeg","method","fishery","zone",
                             "sname1","vessel","port","BoatName","MastersName","crew","Block",
                             "Lat","Long","hooks","netlen","hours","mslow","mshigh","shots",
                             "nlines","nNets","depthMin","depthMax","NilCatch","HookSize",        #NEW check nNets
                             "HookType","nfish","LatMin","LongMin","conditn",
                             "flagtype","factor","livewt","b10days",
                             "bdays","fdays","type","RSCommonName","RSSpeciesCode",
                             "RSSpeciesId","RSBioregionName",
                             "TripLandedWeight")
  Data.daily=Data.daily[,Keep.these.columns_daily]%>%
                mutate(nlines=case_when(method=='GN' & is.na(nlines) ~nNets,
                                        TRUE~nlines))%>%
                dplyr::select(-nNets)

    #re-map vessel codes to previous codes
  Data.daily=Data.daily%>%
    rename(FLAMS=vessel)%>%
    left_join(Vessel.mapping%>%
                rename(FLAMS='LFB In FLAMS Format Without LFB String',
                       CAES='LFB In CAES Format')%>%
                dplyr::select(FLAMS,CAES),
              by='FLAMS')%>%
    rename(vessel=CAES)%>%
    dplyr::select(-FLAMS)
  
  #Extract Monthly returns from TDGDLF and NDS & shark catch in other fisheries with monthly returns 
  system.time({Data.monthly <- sqlQuery(channel=conn.monthly, query=Query.monthly)})
  odbcClose(conn.monthly)
  
    #check duplication
  Data.monthly=Data.monthly%>%
    mutate(Dupli=paste(FinancialYear, Year, Month,VesselRegistration,FishingMethod,CAESBlock,LiveWeight,LandedCondition,RSSpeciesCommonName,RSSpeciesCode)) 
  Dups=table(Data.monthly$Dupli)  
  Dups=names(Dups[Dups>1])
  if(length(Dups)>0)
  {
    if(First.run=='YES')
    {
      handle.duplis=handl_OneDrive("Data/Catch and Effort")
      Duplis=Data.monthly%>%
        filter(Dupli%in%Dups)%>%
        arrange(Dupli)%>%
        dplyr::select(FinancialYear, Year, Month,VesselRegistration,FishingMethod,CAESBlock,
                      LiveWeight,LandedCondition,RSSpeciesCommonName,RSSpeciesCode)
      write.csv(Duplis,file=paste(handle.duplis,"/Duplis_monlthly.csv",sep=""),row.names=F)
      send.email(TO=Email.data.checks,
                 CC=Email.data.checks2,
                 BCC=Email.FishCube,
                 Subject=paste("Duplicated records in monthly returns",Sys.time(),sep=' _ '),
                 Body= "Hi,
              I've attached a spredsheet with duplicated records
              (i.e. same FinancialYear, Month, VesselRegistration, FishingMethod, CAESBlock, LiveWeight, LandedCondition, RSSpeciesCommonName).
              Could you please have a look?
              Cheers
              Matias",  
                 Attachment=paste(handle.duplis,"/Duplis_monlthly.csv",sep="")) 
    }
    Data.monthly=Data.monthly%>%
              distinct(Dupli,.keep_all = T)
  }
  rm(Dups)
  Data.monthly=Data.monthly%>%
        dplyr::select(-Dupli)
  
  rename.monthly=SQL.data.mapping_monthly$OLD
  names(rename.monthly)=SQL.data.mapping_monthly$NEW
  rename.monthly=subset(rename.monthly,!is.na(names(rename.monthly)))
  for(r in 1:length(rename.monthly))
  {
    id=match(names(rename.monthly)[r],colnames(Data.monthly))
    colnames(Data.monthly)[id]=rename.monthly[r]
    rm(id)
  }
  Data.monthly=Data.monthly%>%
                    rename(zone=FisheryZone,
                           SNAME=RSSpeciesCommonName,
                           VESSEL=VesselRegistration,
                           type=RSSpeciesEcologicalSuite)%>%            
                    mutate(SPECIES=case_when(type=='elasmobranchs'~substr(RSSpeciesCode,4,10),
                                             TRUE~substr(RSSpeciesCode,3,10)),
                           SPECIES=as.numeric(SPECIES),
                           SPECIES=case_when(SPECIES== 5902 ~ 5002,   #reset to previous code for consistency thru all scripts
                                             SPECIES== 5901 ~ 5001,
                                             SPECIES== 12001 ~ 12000,
                                             SPECIES== 13900 ~ 13000,
                                             SPECIES== 18901 ~ 18014,
                                             SPECIES== 27909 ~ 26999,
                                             SPECIES== 287003 ~ 288000,
                                             SPECIES== 311010 ~ 311058,
                                             SPECIES== 311912 ~ 311199,
                                             SPECIES== 311991 ~ 311903,
                                             SPECIES== 320004 ~ 320000,
                                             SPECIES== 344004 ~ 344002,
                                             SPECIES== 346014 ~ 346914,
                                             SPECIES== 353013 ~ 353998,
                                             SPECIES== 381010 ~ 381006,
                                             SPECIES== 384039 ~ 384999,
                                             SPECIES== 441911 ~ 441000,
                                             SPECIES== 441006 ~ 441020,
                                             SPECIES== 441912 ~ 441900,
                                             SPECIES== 990009 ~ 460000,
                                             SPECIES== 615000 ~ 600000,
                                             SPECIES== 607901 ~ 610000,
                                             SPECIES== 911005 ~ 702003,
                                             SPECIES== 346015 ~ 346912,
                                             SPECIES== 353000 ~ 353901,
                                             SPECIES== 85018 ~ 358998,
                                             SPECIES== 90001 ~ 990001,
                                             SPECIES== 999997 ~ 599000,
                                             SPECIES== 916901 ~ 702009,
                                             SPECIES== 99993 ~ 22997,
                                             SPECIES== 99994 ~ 22998,
                                             SPECIES== 99992 ~ 22999,
                                             SPECIES== 346916 ~ 346911,
                                             SPECIES== 384902 ~ 384901,
                                             SPECIES== 999997 ~ 599001,
                                             SPECIES== 351000 ~ 351910,
                                             TRUE~SPECIES))
  Mesh.monthly=subset(Data.monthly,METHOD=="GN",select=c(VESSEL,FINYEAR,MONTH,BLOCKX,METHOD,MeshSizeHigh,MeshSizeLow))
 
  Keep.these.columns_monthly=c("FINYEAR","YEAR","MONTH","VESSEL","FDAYS","METHOD","BLOCKX","BDAYS","HOURS","HOOKS","SHOTS",
                               "NETLEN","SPECIES","SNAME","LIVEWT","CONDITN","LANDWT","Factor","fishery","PORT",
                               "RSSpeciesCode","type","zone")
  if('lat'%in%names(Data.monthly))
  {
    Keep.these.columns_monthly=c(Keep.these.columns_monthly,"lat","long")
    Data.monthly=Data.monthly%>%
                    mutate(long=ifelse(nchar(long)<3,100+long,long))
  }
    
  Data.monthly=Data.monthly[,Keep.these.columns_monthly]
  
    #re-map vessel codes to previous codes
  Data.monthly=Data.monthly%>%
                rename(FLAMS=VESSEL)%>%
                left_join(Vessel.mapping%>%
                            rename(FLAMS='LFB In FLAMS Format Without LFB String',
                                   CAES='LFB In CAES Format')%>%
                            dplyr::select(FLAMS,CAES),
                          by='FLAMS')%>%
                rename(VESSEL=CAES)%>%
                dplyr::select(-FLAMS)
  
  #set to bad reporters if catch reapportioned by FishCube 
  Bad.monthly.reporters=read.csv('Bad.monthly.reporters_catch_more_80per_22999.csv')
  Data.monthly=Data.monthly%>%
    mutate(Same.return=paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX),
           Reporter=ifelse(Same.return%in%Bad.monthly.reporters$Same.return,'bad','good'))%>%
    dplyr::select(-Same.return)
  
  
  #Access Monthly previous version
  Acces.Monthly.previous.version=FALSE  #March 2025
  if(Acces.Monthly.previous.version)
  {
    library(stringr)
    Server.monthly.dummy='reports.dpird.wa.gov.au'
    Database.monthly.dummy='ResearchDataWarehouseCE'
    
    conn.dummy <- odbcDriverConnect(connection=paste("Driver={SQL Server};server=",Server.monthly.dummy,
                                                     ";database=",Database.monthly.dummy,
                                                     ";trusted_connection=yes;",sep=""))
    Monthly.log.dummy= 'monthlyrawlog'
    
    Query.dummy=paste("SELECT * FROM",Monthly.log.dummy,
                      "WHERE RSSpeciesEcologicalSuite = 'elasmobranchs' or RSInitiallyReportedFisheryCode in ('SGL','WCGL','C127','WANCS','CL02','JANS')")
    
    
    Data.monthly.old <- sqlQuery(channel=conn.dummy, query=Query.dummy)
    
    Data.monthly.old=Data.monthly.old%>%
      mutate(FDAYS=NA,BDAYS=NA,Reporter=NA,
             RSSpeciesCode=substr(RSSpeciesCode,4,10),
             SPECIES=RSSpeciesCode,)%>%
      rename(FINYEAR=FinancialYear,
             YEAR=Year,
             MONTH=Month,
             METHOD=FishingMethod,
             BLOCKX=CAESBlock,
             HOURS=nHoursFishing,
             HOOKS=nHooks,
             SHOTS=nShots,
             NETLEN=NetLength,
             SNAME=RSSpeciesCommonName,
             LIVEWT=LiveWeight,
             CONDITN=LandedCondition,
             LANDWT=LandedWeight,
             RSSpeciesCode=RSSpeciesCode,
             Factor=ConversionFactor,
             fishery=RSInitiallyReportedFisheryCode,
             PORT=DepartPort,
             type=RSSpeciesEcologicalSuite,
             VESSEL=VesselRegistration)
    Data.monthly.old=Data.monthly.old[,match(names(Data.monthly),names(Data.monthly.old))]
    Data.monthly.old=Data.monthly.old%>%
      mutate(SPECIES=ifelse(SNAME=='Sharks',22999,SPECIES))

    Data.monthly.old=Data.monthly.old%>%
      mutate(VESSEL=case_when(VESSEL=="A100A"~"A 100A",      
                              VESSEL=="A197A"~"A 197A",
                              VESSEL=="A238A"~"A 238A",
                              VESSEL=="F188A"~"F 188A",
                              VESSEL=="F200A"~"F 200A",
                              nchar(VESSEL)==2~gsub("(\\d*)(\\D*)\\s*(\\d*)","\\2 00\\3",VESSEL),
                              nchar(VESSEL)==3~gsub("(\\d*)(\\D*)\\s*(\\d*)","\\2 0\\3",VESSEL),
                              nchar(VESSEL)==4~gsub("(\\d*)(\\D*)\\s*(\\d*)","\\2 \\3",VESSEL),
                              TRUE~VESSEL),
             VESSEL=case_when(VESSEL=="A 07A 0"~"A 007A", 
                              VESSEL%in%paste0("BR 0",1:7)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("BR ",10:20)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("DY 0",1:9)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("DY ",10:20)~str_replace(VESSEL," ","0"),
                              VESSEL=="DY 08"~"DY008",
                              VESSEL=="DY 09"~"DY009",
                              VESSEL=="F 010"~"F10",
                              VESSEL=="F 105"~"F105",
                              VESSEL=="F 221"~"F221",
                              VESSEL=="F 248"~"F248",
                              VESSEL=="F 540"~"F540",
                              VESSEL=="F 550"~"F550",
                              VESSEL=="F 600"~"F600",
                              VESSEL=="F 711"~"F711",
                              VESSEL=="G 296"~"G296",
                              VESSEL=="PH 5A "~"PH005A",
                              VESSEL%in%paste0("PH 0",1:7)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("PS 0",1:9)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("PS ",10:28)~str_replace(VESSEL," ","0"),
                              VESSEL=="SB 1A "~"SB001A",
                              VESSEL%in%paste0("SB 0",1:9)~str_replace(VESSEL," ","0"),
                              VESSEL%in%paste0("SB ",10:74)~str_replace(VESSEL," ","0"),
                              TRUE~VESSEL))
  }
  
  #3. Daily returns from other fisheries that have reported shark
  Data.mart.Logs=tolower(c('ceGDSLog','ceMACLog','ceNDSLog','cePFTLog','ceCWCDSLog'))
  Data.daily_other=vector('list',length(Data.mart.Logs))
  names(Data.daily_other)=Data.mart.Logs
  Data.monthly_other=Data.daily_other
  for(l in 1:length(Data.daily_other))
  {
    dummy <- sqlQuery(conn1, query=paste("SELECT * FROM",Data.mart.Logs[l],"WHERE RSSpeciesEcologicalSuite = 'elasmobranchs'"))
    if(is.data.frame(dummy))
    {
      if(nrow(dummy)>0)
      {
        if(!'RSCommonName'%in%names(dummy)) dummy$RSCommonName=dummy$RSSpeciesCommonName
        if('RSCAABCode'%in%names(dummy)) dummy$RSSpeciesCode=dummy$RSCAABCode
        dummy=dummy%>%
          mutate(sname1=RSCommonName,
                 species=case_when(RSSpeciesEcologicalSuite=='elasmobranchs'~substr(RSSpeciesCode,4,10),
                                   TRUE~substr(RSSpeciesCode,3,10)),
                 species=as.numeric(species),
                 species=case_when(species== 5902 ~ 5002,   #reset to previous code for consistency thru all scripts
                                   species== 5901 ~ 5001,
                                   species== 12001 ~ 12000,
                                   species== 13900 ~ 13000,
                                   species== 18901 ~ 18014,
                                   species== 27909 ~ 26999,
                                   species== 287003 ~ 288000,
                                   species== 311010 ~ 311058,
                                   species== 311912 ~ 311199,
                                   species== 311991 ~ 311903,
                                   species== 320004 ~ 320000,
                                   species== 344004 ~ 344002,
                                   species== 346014 ~ 346914,
                                   species== 353013 ~ 353998,
                                   species== 381010 ~ 381006,
                                   species== 384039 ~ 384999,
                                   species== 441911 ~ 441000,
                                   species== 441006 ~ 441020,
                                   species== 441912 ~ 441900,
                                   species== 990009 ~ 460000,
                                   species== 615000 ~ 600000,
                                   species== 607901 ~ 610000,
                                   species== 911005 ~ 702003,
                                   species== 346015 ~ 346912,
                                   species== 353000 ~ 353901,
                                   species== 85018 ~ 358998,
                                   species== 90001 ~ 990001,
                                   species== 999997 ~ 599000,
                                   species== 916901 ~ 702009,
                                   species== 99993 ~ 22997,
                                   species== 99994 ~ 22998,
                                   species== 99992 ~ 22999,
                                   species== 346916 ~ 346911,
                                   species== 384902 ~ 384901,
                                   species== 999997 ~ 599001,
                                   species== 351000 ~ 351910,
                                   TRUE~species),
                 spgroup=case_when(species %in% 5001:24900 ~ 'sharks',
                                   species %in% c(25000:31000,39001,990001) ~ 'rays',
                                   RSSpeciesEcologicalSuite=='inverts'~'inverts',
                                   TRUE~'scalefish'))  
        if(!"LongitudeDecimalDegrees"%in%colnames(dummy) & "SessionMiddleLongitudeDD"%in%colnames(dummy))
        {
          dummy$LongitudeDecimalDegrees=dummy$SessionMiddleLongitudeDD
          dummy$LatitudeDecimalDegrees=dummy$SessionMiddleLatitudeDD
        }
        if(!"Block10"%in%colnames(dummy))
        {
          dummy=dummy%>%
            mutate(Block10=Block,
                   CAESBlock=SessionMiddleCAESBlock,
                   CAESBlock=ifelse(is.na(CAESBlock) & Fishery=="PFT",19118,CAESBlock))
        }
        
        dummy=dummy%>%
          mutate(LatitudeDecimalDegrees=ifelse(is.na(LatitudeDecimalDegrees),
                                               -abs(as.numeric(substr(Block10,1,2))+(10*as.numeric(substr(Block10,3,3))/60)),  
                                               LatitudeDecimalDegrees),
                 LongitudeDecimalDegrees=ifelse(is.na(LongitudeDecimalDegrees),
                                                100+as.numeric(substr(Block10,4,5))+(10*as.numeric(substr(Block10,6,6))/60),
                                                LongitudeDecimalDegrees),
                 LongitudeDecimalDegrees=ifelse(is.na(LongitudeDecimalDegrees) & Fishery=="PFT",118,LongitudeDecimalDegrees),
                 LatitudeDecimalDegrees=ifelse(is.na(LatitudeDecimalDegrees) & Fishery=="PFT",-19.5,LatitudeDecimalDegrees))
        
        #split daily from monthly
        dummy.daily=dummy%>%filter(FinancialYear%in%unique(Data.daily$finyear))
        dummy.monthly=dummy%>%filter(!FinancialYear%in%unique(Data.daily$finyear))
        
        #daily
        for(r in 1:length(rename.daily))
        {
          kk=rename.daily[r]
          if(kk=="fishery" & 'Fishery'%in%colnames(dummy.daily)) names(kk)='Fishery'
          id=match(names(kk),colnames(dummy.daily))
          colnames(dummy.daily)[id]=kk
          rm(id)
        }
        IId=which(!Keep.these.columns_daily%in%names(dummy.daily))
        if(length(IId)>0)
        {
          dummy2=as.data.frame(matrix(nrow=nrow(dummy.daily),ncol=length(Keep.these.columns_daily[IId])))
          colnames(dummy2)=Keep.these.columns_daily[IId]
          dummy.daily=cbind(dummy.daily,dummy2)
          rm(dummy2)
        }
        dummy.daily=dummy.daily[,Keep.these.columns_daily]
        Data.daily_other[[l]]=dummy.daily
        rm(dummy.daily)
        
        #monthly
        if(nrow(dummy.monthly)>0)
        {
          dummy.monthly=dummy.monthly%>%
            rename(SNAME=sname1,
                   SPECIES=species,
                   fishery=Fishery,
                   VESSEL=LFB,
                   LIVEWT=LiveWeight,
                   PORT=DeparturePort)%>%
            mutate(METHOD=ifelse(fishery=='PFT','TW',NA))
          for(r in 1:length(rename.monthly))
          {
            id=match(names(rename.monthly)[r],colnames(dummy.monthly))
            colnames(dummy.monthly)[id]=rename.monthly[r]
            rm(id)
          }
          IId=which(!Keep.these.columns_monthly%in%names(dummy.monthly))
          if(length(IId)>0)
          {
            dummy2=as.data.frame(matrix(nrow=nrow(dummy.monthly),ncol=length(Keep.these.columns_monthly[IId])))
            colnames(dummy2)=Keep.these.columns_monthly[IId]
            dummy.monthly=cbind(dummy.monthly,dummy2)
            rm(dummy2)
          }
          dummy.monthly=dummy.monthly[,Keep.these.columns_monthly]
          Data.monthly_other[[l]]=dummy.monthly
          rm(dummy.monthly)
        }
      }
    }
    rm(dummy)
  }
  Data.daily_other=compact(Data.daily_other)
  Data.daily_other=do.call(rbind,Data.daily_other)
  Data.monthly_other=compact(Data.monthly_other)
  Data.monthly_other=do.call(rbind,Data.monthly_other)
  if('lat'%in%names(Data.monthly_other))
  {
    Data.monthly_other=Data.monthly_other%>%  
      mutate(lat=ifelse(is.na(lat),-as.numeric(substr(BLOCKX,1,2)),lat),
             long=ifelse(is.na(long),as.numeric(substr(BLOCKX,3,5)),long))
  }
  #add other' to 'main databases
  Data.monthly_other$Reporter=NA
  Data.daily=rbind(Data.daily,Data.daily_other[,match(names(Data.daily),names(Data.daily_other))]) 
  Data.monthly=rbind(Data.monthly,Data.monthly_other[,match(names(Data.monthly),names(Data.monthly_other))])
  
  rm(Data.monthly_other,Data.daily_other)
}
  #1.2 Extract catch and effort data from Excel and Access databases (superseded by SQL)
if(!do.sql.extraction)
{
  #1. Monthly records 
  
  #CAESS data 
  #note: CAESS is dynamic, constantly updated, Table81.d is static so only use Table 81.d for 
  #       years prior to 1988-89 as CAESS data starts in 1988-89.
  #       Also, don't use CAESS from non-TDGDLF or NSF fisheries for catch rate standardisations as
  #       0 catch records are missing/incomplete
  #       Data.monthly.CAESS has all WA fisheries reporting shark
  Get.CAESS.Logbook="YES"
  
  #CAES extract done by Paul Fildes from FishCUBE. This has all fisheries reporting shark
  if(Get.CAESS.Logbook=="YES")Data.monthly.CAESS=read.csv("M:/CAEMaster/Commercial/FishCubeWA/Fisheries/CAES Monthly Fisheries/02. Data/CAES Monthly Data - Shark.csv") 		
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
    keep=c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81","1981-82",
           "1982-83","1983-84","1984-85","1985-86","1986-87","1987-88")
    Data.monthly=subset(Data.monthly,FINYEAR%in%keep)%>%
      mutate(fishery=NA,
             licence=NA,
             PORT=NA,
             rowid=NA)
    Data.monthly.CAESS=Data.monthly.CAESS%>%
      rename(LIVEWT=livewt,
             Factor=FACTOR)
    if(!"SNAME"%in%names(Data.monthly.CAESS)) Data.monthly.CAESS=Data.monthly.CAESS%>%rename(SNAME=sname1)
    ID=match(names(Data.monthly),names(Data.monthly.CAESS))
    Data.monthly.CAESS=Data.monthly.CAESS[,ID]
    Data.monthly=rbind(Data.monthly,Data.monthly.CAESS)
  }
  
  
  #2. Daily logbooks
  #Select previous daily logbooks (static) or Eva's dynamic data dump
  #Get.Daily.Logbook="NO"  #use Rory's files
  Get.Daily.Logbook="YES"  #use SADA's annual extraction
  #Daily records from shark fisheries
  Daily.source="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/Shark/02. Data/"
  Daily.file="1. Shark Daily Logbook Data.xlsx"
  Daily.dat=paste(Daily.source,Daily.file,sep="")
  channel <- odbcConnectExcel2007(Daily.dat)
  Data.daily<- sqlFetch(channel,"CATCH", colnames = F)
  close(channel)
  

  #3. Daily records from other fisheries reporting shark catch in daily returns
  #note: run this FishCube query (copy to browser and update file in folder) 
  #      http://F01-FIMS-WEBP01/FishCubeWA/Query.aspx?CubeId=CommercialDPIRDOnly&QueryId=7c2d0c88-fe4a-4818-8d3a-0405eab6199b
  #and update the excel file that the query creates
  channel <- odbcConnectExcel2007("Fish Cube WA_daily.other.xlsx") 
  Data.daily.other.fisheries<- sqlFetch(channel,"Data")
  close(channel)
  channel <- odbcConnectExcel2007("Fish Cube WA_daily.other_names.xlsx") 
  Data.daily.other.fisheries.names<- sqlFetch(channel,"Sheet1")
  close(channel)
  Data.daily.other.fisheries.names=Data.daily.other.fisheries.names%>%
    dplyr::select(RSSpeciesId,species,CAESname)%>%
    filter(RSSpeciesId%in%unique(Data.daily.other.fisheries$`Species Id`)&
             species<31000)
  
  Data.daily.other.fisheries=Data.daily.other.fisheries%>%
    rename(method='Fishing Method Code',
           fishery='Fishery Code',
           fishery.name='Fishery Name',
           livewt='Weight (kg) Total',
           date=Date,
           block10='10x10NM Block',
           finyear='Financial Year',
           year='Calendar Year',
           RSCommonName='Species Common Name')%>%
    left_join(Data.daily.other.fisheries.names,by=c("Species Id" = "RSSpeciesId"))%>%
    dplyr::select(Bioregion,date,Vessel,method,fishery,fishery.name,block10,finyear,year,
                  species,RSCommonName,CAESname,livewt)
  
}


#2. TEPS                                
if(do.sql.extraction)
{
  #Prior to 2019 (not available in SQL)
  TEPS.prior2019=sqlFetch(channel=odbcConnectExcel2007('M:/Fisheries Research/Stock Assessment & Data Analysis/CAEMaster/Commercial/SAS Programs/Shark/Output/SHARK Protected Species.xlsx'),
                          'PROTECTEDSP')
  TEPS.prior2019_comments=sqlFetch(channel=odbcConnectExcel2007('M:/Fisheries Research/Stock Assessment & Data Analysis/CAEMaster/Commercial/SAS Programs/Shark/Output/SHARK Comments.xlsx'),
                                   'COMMDATA')
  TEPS.prior2019=TEPS.prior2019%>%
                left_join(TEPS.prior2019_comments%>%
                            filter(DSNo%in%unique(TEPS.prior2019$DailySheetNumber))%>%
                            distinct(Comments,DSNo),
                          by=c('DailySheetNumber'='DSNo'))
  
  #2019 onwards   
   TEPS.current <- sqlQuery(conn1, query="SELECT  etp.*, s.StartBlock AS SessionStartBlock
                                     FROM [ResearchDataWarehouseQuery].[dbo].[flTripEtp] etp
                                     LEFT JOIN [ResearchDataWarehouseQuery].[dbo].[flTripSession] s ON etp.TripSessionId = s.TripSessionId
                                     WHERE etp.[FisheryCode] in ('WCDGDL','OASC,OT','JASDGDL','SDGDL')")
  TEPS.current=TEPS.current%>%
    rename(DailySheetNumber=TripHeaderId,
         Comments=Notes,
         fishery=FisheryCode,
         Status=CatchCode,
         RSCommonName=SpeciesName)%>%
    mutate(CommonName=RSCommonName,
           finyear=ifelse(TripMonth>6,paste(TripYear,substr(TripYear+1,3,4),sep='-'),
                          paste(TripYear-1,substr(TripYear,3,4),sep='-')))%>%
    dplyr::select(DailySheetNumber,finyear,fishery,Status,SpeciesCode,
                  SheetStartDate,VesselName,RSCommonName,CommonName,
                  Number,SessionStartBlock,Comments)%>%
    rename(date=SheetStartDate,
                     vessel=VesselName,
           blockx=SessionStartBlock)%>%
    mutate(blockx=as.numeric(paste(substr(blockx,1,2),substr(blockx,4,5),0,sep='')))

  
  #combine prior to 2019 and current
  TEPS.current=rbind(TEPS.current,
                     TEPS.prior2019%>%dplyr::select(names(TEPS.current))%>%
                       mutate(CommonName=RSCommonName)%>%
                       filter(!finyear%in%unique(TEPS.current$finyear)))%>%
                mutate(fishery=case_when(fishery=="SGL"~'JASDGDL',
                                         fishery=="*"~'OASC,OT',
                                         fishery=="WCGL"~'WCDGDL',
                                         TRUE~fishery))
  

}
if(!do.sql.extraction)
{
  #2.1. Current
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
  
  #2.2. Historic TEP data (Prior to 2012 report)
  Results.pre.2013=read.csv("Historic/Historic.res.csv")
  TEPS.pre.current=read.csv("Historic/Historic.TEPS.res.csv",stringsAsFactors =F)
  Spec.catch.zone.pre.2013=read.csv("Historic/Spec.catch.zone.csv")  
}


#3. Catch price 
if(do.sql.extraction)
{
  PRICES <- sqlQuery(conn1, "select * from dbo.rsSpeciesPrice") 
  Current.yr.price=Current.yr
  if(!Current.yr.price%in%unique(PRICES$FinancialYear)) Current.yr.price=max(unique(PRICES$FinancialYear))
  PRICES=PRICES%>%
          filter(FinancialYear==Current.yr.price)%>%
          rename(RSSpeciesId=SpeciesId,
                 CommonName=SpeciesCommonName,
                 'Beach Price (Adjusted)'=UnitValue)%>%
          dplyr::select(RSSpeciesId,CommonName,'Beach Price (Adjusted)')
}
if(!do.sql.extraction)
{
  PRICES=sqlFetch(channel=odbcConnectExcel2007(paste(Current.yr.dat,"/PriceComparison.xlsx",sep="")),
                  'Sheet1')
}
PRICES=PRICES%>%
          left_join(Species.mapping%>%dplyr::select(RSSpeciesId,species),by='RSSpeciesId')%>%
          rename(ASA.Species.Code=species)


#4. SHAPE FILE PERTH ISLANDS
# PerthIs=read.table("C:/Matias/Data/Mapping/WAislandsPointsNew.txt", header=T)
# Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
# Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))


#5. bathymetry data
#    bathymetry data downloaded from http://topex.ucsd.edu/cgi-bin/get_data.cgi (Topography option)
Bathymetry_120=read.table(handl_OneDrive("Data/Mapping/get_data112_120.cgi"))
Bathymetry_138=read.table(handl_OneDrive("Data/Mapping/get_data120.05_138.cgi"))
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)


#6. Block10 locations
BlOCK_10=read.csv(handl_OneDrive("Data/Mapping/Blocks_10NM.csv"))


#7. All block 60 lat and long, including estuaries
BLOCK_60=read.csv(handl_OneDrive("Data/Mapping/Block60s.csv"))


#8. Weight ranges
Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"))
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"))


#9. Conditions
#Conditions=read.csv(handl_OneDrive("Data/Catch and Effort/Conditions.csv"))  #superseeded by fish.conversion


#10. Rory's manual changes to netlen and nlines
Rory_Alex_net_val=read.csv(handl_OneDrive("Data/Catch and Effort/Rory_Alex_net/Book2.csv"),stringsAsFactors=F)


#11. ASL block10 closures
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


#12. define estuaries
Estuaries=c(95010,95020,95030,95040,95050,95060,95070,95080,95090,85010,85020,85030,85040,85050,85060,
            85070,85080,85090,85100,85110,85120,85130,85210,85220,85990)
Bays=96000:98000
Shark.Bay=96021
Abrolhos=97011
Geographe.Bay=96010 
Cockburn.sound=96000
King.George.sound=96030 
85010-> Hardy.inlet  
85020-> Parrys.inlet
85030-> Beaufort.inlet
85050-> Gordon.inlet
85080-> Culham.inlet


95010-> Swan.inlet
95020-> Peel.inlet
95030-> Leschenault.inlet

95050-> Oyster.harbour
95060-> Wilson.inlet
95070-> Irwin.inlet
95080-> Broke.inlet


# ---------- PARAMETERS SECTION ----------------------------------------------
#control if inspecting the new data 
if(First.run=="NO")Inspect.New.dat="NO"
if(First.run=="YES")Inspect.New.dat="YES"

#control if plotting spatio-temporal trends in catch and effort for TDGDLF
do.spatio.temporal.ktch.effort=First.run
do.tdgdlf.effort.density="NO"

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
do.data.requests="YES"
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
do.Jeffs="NO"   #old request
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
do.annual.TEPS.extraction="NO"
do.Paul.Rogers_ASL="NO"
do.ASL.closure_effort_overlap="NO"
Check.school.shark.targeting="NO"
do.Peter.Rogers=FALSE
do.DBCA_2022=FALSE
get.Mats.data.2023=FALSE
do.Ngari.2023=FALSE
do.SARDI.2024=FALSE

#Spatial range TDGDLF
TDGDLF.lat.range=c(-26,-40)

#Fin to total weight ratio
Percent.fin.of.livewt=0.03

#Shark fisheries
FishCubeCode_Shark.fisheries=c('C070','JASDGDL','WCDGDL','OANCGCWC','WANCS','JANS','NCS','OASC')

#Species definition    
Shark.species=5001:24900
Ray.species=c(25000:31000,39001,90030)
Scalefish.species=117001:599001   
Indicator.species=c(17001,17003,18001,18003,18007)
names(Indicator.species)=c("Gummy","Whiskery","Bronzy","Dusky","Sandbar")
TARGETS=list(whiskery=17003,gummy=17001,dusky=c(18001,18003),sandbar=18007) 
N.species=length(TARGETS)
Fix.species=c(18003,17001,17003,22999)  #species that need catch reapportioning
HammerheadSpecies=c(19001,19002,19004)
reset.hammerhead.to.reported=FALSE
if(use.new.server) reset.hammerhead.to.reported=TRUE    #if using new server, reset hammerhead to originally reported record
                                                        # because FishCube does an automatic reapportioning

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
Wei.range=subset(Wei.range,!(is.na(SPECIES)|is.na(TW.min)|is.na(TW.max)))
max.w.whis=Wei.range[Wei.range$SPECIES==17003,]$TW.max   
max.w.gum=Wei.range[Wei.range$SPECIES==17001,]$TW.max
max.w.dus=Wei.range[Wei.range$SPECIES==18003,]$TW.max
max.w.san=Wei.range[Wei.range$SPECIES==18007,]$TW.max

min.w.whis=Wei.range[Wei.range$SPECIES==17003,]$TW.min  
min.w.gum=Wei.range[Wei.range$SPECIES==17001,]$TW.min
min.w.dus=Wei.range[Wei.range$SPECIES==18003,]$TW.min
min.w.san=Wei.range[Wei.range$SPECIES==18007,]$TW.min


 #how to export figures
Do.tiff="NO"
Do.jpeg="YES"



#--- Quick check to determine samples size catch age composition ---------

yr=as.numeric(substr(Current.yr,1,4))
if(explore.Catch.compo=="YES")
{
  hnd.compo=handl_OneDrive("Fieldwork and workplans/Catch age composition/")
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
                             dplyr::select(port,species,nfish,vessel,TSNo,finyear))
}



#SECTION A. ---- DATA MANIPULATION - MONTHLY RECORDS ----

#Set hammerhead species to the original report to overwrite FishCube's reapportioning  
if(reset.hammerhead.to.reported)
{
  all.hammers=Data.monthly%>%
    filter(SPECIES%in%c(19000,HammerheadSpecies))%>%
    mutate(SNAME=ifelse(SPECIES%in%HammerheadSpecies,"Hammerhead Sharks",SNAME),
           RSSpeciesCode=ifelse(SPECIES%in%HammerheadSpecies,37019000,RSSpeciesCode),
           SPECIES=ifelse(SPECIES%in%HammerheadSpecies,19000,SPECIES))
  Nms.hmrs=names(all.hammers)
  Nms.hmrs=Nms.hmrs[-match(c('LIVEWT','LANDWT'),Nms.hmrs)]
  all.hammers=all.hammers%>%
    group_by_at(Nms.hmrs)%>%
    summarise(LANDWT=sum(LANDWT,na.rm=T),
              LIVEWT=sum(LIVEWT,na.rm = T))
  Data.monthly=Data.monthly%>%
    filter(!SPECIES%in%c(19000,HammerheadSpecies))
  Data.monthly=rbind(Data.monthly,all.hammers%>%relocate(names(Data.monthly)))
}


# A.1. Add year if missing
Data.monthly$Split.Year.1=as.numeric(sapply(strsplit(as.character(Data.monthly$FINYEAR),"-"), "[", 1))
Data.monthly$Split.Year.2=as.numeric(sapply(strsplit(as.character(Data.monthly$FINYEAR),"-"), "[", 2))
Data.monthly$Split.Year.2=with(Data.monthly,ifelse(Split.Year.2<75,Split.Year.2+2000,Split.Year.2+1900))
Data.monthly$YEAR.c=with(Data.monthly,ifelse(is.na(YEAR) & MONTH>=7,Split.Year.1,
                          ifelse(is.na(YEAR) & MONTH<7,Split.Year.2,YEAR)))


# A.2. Remove records post June 2006 if using Table 89.d only
if(!do.sql.extraction)
{
  Data.monthly$ID=with(Data.monthly,ifelse(YEAR.c>2006,"drop","keep"))
  Data.monthly$ID=with(Data.monthly,ifelse(YEAR.c==2006 & MONTH>6,"drop",ID))
  if(Get.CAESS.Logbook=="NO") Data.monthly=subset(Data.monthly,ID=="keep")
  Data.monthly=Data.monthly[,-match(c("ID"),names(Data.monthly))]
  
}


# A.3. Remove dummies
Data.monthly=Data.monthly[,-match(c("Split.Year.1","Split.Year.2"),names(Data.monthly))]


# A.4. Sort by year and month and vessel
Data.monthly=Data.monthly[order(Data.monthly$YEAR.c,Data.monthly$MONTH,Data.monthly$VESSEL),]


# A.4.2 Create dummy to group same return (finyear-month-vessel-method-block)             
Data.monthly$Same.return=with(Data.monthly,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
Data.monthly$Same.return.SNo=Data.monthly$Same.return
Data.monthly$TYPE.DATA="monthly"


#remove fins and livers to avoid duplication when calculating live weight
  #note: id records were only fins or livers were reported.
  # for this, check which species were wrongly assigned a liver
  # or fin code and remove records from those where liver/fin and
  # WF or WH weights were reported
A=Data.monthly%>%
        filter(CONDITN%in%c("LI","FI") &
               !SPECIES%in%c(22997,22998))
if(nrow(A)>0)
{
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
    if(!is.null(Tbl.x$LI)) KEEP=subset(Tbl.x,(Sums==1 & FI==1) | (Sums==1 & LI==1)) else
      KEEP=subset(Tbl.x,(Sums==1 & FI==1))
    keep.these[[s]]=row.names(KEEP)
  }
  keep.these=unlist(keep.these)
  keep.these=keep.these[!duplicated(keep.these)]
  A$KEEP=with(A,ifelse(Same.return%in%keep.these,"KEEP","DROP"))
  A$KEEP=with(A,ifelse(!CONDITN%in%c("FI","LI"),"KEEP",KEEP))
}
if(nrow(A)>0) Data.monthly$KEEP=with(Data.monthly,ifelse(CONDITN%in%c("FI","LI"),"DROP","KEEP"))
if(nrow(A)>0) Data.monthly=rbind(Data.monthly,A)

  #a. reconstruct landed and livewt records that only reported 'liver' or 'fin'
only.liver.fin=Data.monthly%>%
                  filter(SPECIES%in%c(22997,22998))%>%
                  distinct(Same.return)%>%
                  pull(Same.return)
if(!use.new.server & length(only.liver.fin)>0) 
{
  only.liver.fin=Data.monthly%>%
    filter(Same.return%in%only.liver.fin & SPECIES<50000)%>%
    mutate(SNAME=tolower(SNAME))%>%
    group_by(Same.return,SNAME)%>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))%>%
    dplyr::select(-n)%>%
    spread(SNAME,freq,fill=0)%>%
    data.frame%>%
    mutate(Prop=shark.fins+shark.liver)
  Prop.fin=only.liver.fin%>%
    filter(shark.fins>0)%>%
    pull(Same.return)
  Prop.fin=Data.monthly%>%
    filter(Same.return%in%Prop.fin & SPECIES<50000)%>%
    mutate(fin.no.fin=ifelse(SPECIES==22998,"Fin","No.Fin"))%>%
    group_by(fin.no.fin)%>%
    summarise(livewt=sum(LIVEWT),
              landwt=sum(LANDWT))
  Prop.fin.livewt=Prop.fin$livewt[1]/Prop.fin$livewt[2]
  Prop.fin.landwt=Prop.fin$landwt[1]/Prop.fin$landwt[2]
  
  Prop.livr=only.liver.fin%>%
    filter(shark.liver>0)%>%
    pull(Same.return)
  Prop.livr=Data.monthly%>%
    filter(Same.return%in%Prop.livr & SPECIES<50000)%>%
    mutate(livr.no.livr=ifelse(SPECIES==22997,"Livr","No.Livr"))%>%
    group_by(livr.no.livr)%>%
    summarise(livewt=sum(LIVEWT),
              landwt=sum(LANDWT))
  Prop.livr.livewt=Prop.livr$livewt[1]/Prop.livr$livewt[2]
  Prop.livr.landwt=Prop.livr$landwt[1]/Prop.livr$landwt[2]
  
  add.only.liver.fin=only.liver.fin%>%filter(Prop==1)%>%pull(Same.return)
  add.only.liver.fin=Data.monthly%>%
    filter(Same.return%in%add.only.liver.fin & 
             SPECIES%in%c(22997,22998))%>%
    mutate(SPECIES=22999,
           SNAME='shark, other',
           LIVEWT=ifelse(SPECIES==22998,LIVEWT/Prop.fin.livewt,
                         LIVEWT/Prop.livr.livewt),
           LANDWT=ifelse(SPECIES==22998,LANDWT/Prop.fin.landwt,
                         LANDWT/Prop.livr.livewt),
           CONDITN='WF',
           KEEP='KEEP')
  
  Data.monthly=Data.monthly%>%     
    filter(!(Same.return%in%add.only.liver.fin$Same.return & 
               SPECIES%in%c(22997,22998)))
  Data.monthly=rbind(Data.monthly,add.only.liver.fin)
}

if(!do.sql.extraction)
{
  Data.daily.other.fisheries=Data.daily.other.fisheries%>%
    filter(!species==22998)    #remove fins as already part of 'shark other'
}

#b. now remove fins and livers 
Data.monthly=subset(Data.monthly,!(SPECIES%in%c(22997,22998)))
if(nrow(A)>0) Data.monthly=subset(Data.monthly,KEEP=="KEEP")
Data.monthly$CONDITN=with(Data.monthly,
      ifelse(CONDITN%in%c("FI","LI"),"WF",CONDITN))
if(nrow(A)>0) Data.monthly=Data.monthly[,-match("KEEP",names(Data.monthly))]
if(!use.new.server & length(only.liver.fin)>0) rm(A,add.only.liver.fin,only.liver.fin)

# A.4.3. Create copy of original file
Data.monthly.original=Data.monthly
if(!'lat'%in%colnames(Data.monthly.original)) Data.monthly.original=Data.monthly.original%>%mutate(lat=NA,long=NA)
Data.monthly.original=Data.monthly.original%>%
                        rename(LAT=lat,
                               LONG=long)%>%
                        mutate(LAT=-abs(LAT),
                               LAT=ifelse(is.na(LAT),-as.numeric(substr(BLOCKX,1,2)),LAT),
                               LONG=ifelse(is.na(LONG),100+as.numeric(substr(BLOCKX,3,4)),LONG))

fn.chk.ktch=function(d1,d2,VAR1,VAR2)
{
  d1sum=round(sum(d1[,match(VAR1,names(d1))],na.rm=T))
  d2sum=round(sum(d2[,match(VAR2,names(d2))],na.rm=T))
  if(d1sum>d2sum)Message="d1 has more catch"
  if(d1sum<d2sum)Message="d2 has more catch"
  if(d1sum==d2sum)Message="d1 and d2 have same catch"
  print(Message)
  
  A=aggregate(d1[,match(VAR1,names(d1))]~FINYEAR+VESSEL+Same.return,d1,sum)
  names(A)[4]="weight.original"
  B=aggregate(d2[,match(VAR2,names(d2))]~FINYEAR+VESSEL+Same.return,d2,sum)
  names(B)[4]="weight.changed"
  D=A%>%full_join(B,by=c("FINYEAR","VESSEL","Same.return"))
  #D=merge(A,B,by=c("FINYEAR","VESSEL"))
  D$delta.w=round(D$weight.original-D$weight.changed,0)
  discrepancy=subset(D,!delta.w==0)
  print(paste("Catch by year-vessel has this many discrepancies=",
              nrow(discrepancy),"records"))
  if(nrow(discrepancy)>0) return(discrepancy)
}
fn.chk.ktch(d1=Data.monthly.original,
      d2=subset(Data.monthly,FINYEAR%in%Data.monthly.original$FINYEAR),
      VAR1="LIVEWT",VAR2="LIVEWT")


# A.5. Create variables

  #Identify estuaries
Data.monthly$Estuary=with(Data.monthly,ifelse(BLOCKX%in%Estuaries,"YES","NO"))
test.Estuaries=subset(Data.monthly,BLOCKX%in%Estuaries)
Table.Estuary=aggregate(LIVEWT~FINYEAR+SPECIES,test.Estuaries,sum)  
#note: there is shark catch in estuaries, must keep records for total catch taken


  # A.5.1 Lat and Long of block (top left corner)

#get lat and long from block but first dummy the bays                                                                      
Data.monthly=Data.monthly%>%
  mutate(BLOCKX=case_when(BLOCKX=='3518s'~'35181',
                          BLOCKX=='2613s'~'26131 ',
                          BLOCKX=='2513s'~'96021',
                          BLOCKX=='2413s'~'24131',
                          TRUE~BLOCKX),
         BLOCKX=as.numeric(BLOCKX),
         blockxFC=BLOCKX,
         BLOCKX=case_when(BLOCKX%in% Shark.Bay ~ 25120,
                          BLOCKX%in% c(96022,96023) ~ 26131,
                          BLOCKX%in% Abrolhos ~ 27132,
                          BLOCKX%in% c(97012,97013) ~ 28132,
                          BLOCKX%in% c(97014,97015) ~ 29132,
                          BLOCKX%in% Geographe.Bay ~ 33151,
                          BLOCKX%in% Cockburn.sound ~ 32150,
                          BLOCKX%in% King.George.sound ~ 35181,
                          BLOCKX%in% Swan.inlet ~ 31150,
                          BLOCKX%in% Peel.inlet ~ 32150,
                          BLOCKX%in% Leschenault.inlet ~ 33150,
                          BLOCKX%in% Oyster.harbour ~ 35181,
                          BLOCKX%in% Parrys.inlet ~ 35170,
                          BLOCKX%in% Wilson.inlet ~ 34170,
                          BLOCKX%in% Irwin.inlet ~ 34160,
                          BLOCKX%in% Broke.inlet ~ 34160,
                          TRUE~BLOCKX))        

# Dumi.blk=BLOCK_60%>%
#           mutate(BlockCode=ifelse(grepl("s",BlockCode),paste(substr(BlockCode,1,4),0,sep=''),BlockCode),
#                  BLOCKX=as.numeric(BlockCode),
#                  LAT=ifelse(BLOCKX==85990,NA,NorthWestPointGPSLatitude),
#                  LONG=ifelse(BLOCKX==85990,NA,NorthWestPointGPSLongitude))%>%
#           select(BLOCKX,LAT,LONG)%>%
#           filter(BLOCKX%in%unique(Data.monthly$BLOCKX))

if(!'lat'%in%colnames(Data.monthly)) Data.monthly=Data.monthly%>%mutate(lat=NA,long=NA)
Data.monthly=Data.monthly%>%
              rename(LAT=lat,
                     LONG=long)%>%
              mutate(LAT=-abs(LAT),
                     LAT=ifelse(is.na(LAT),-as.numeric(substr(BLOCKX,1,2)),LAT),
                     LONG=ifelse(is.na(LONG),100+as.numeric(substr(BLOCKX,3,4)),LONG))
Data.monthly$BLOCKX=Data.monthly$blockxFC   #reset block

# adjust special cells to actual lat and long                                              
#note: use ceiling for lat and floor for long to allocate lat and long of north-west corner
Data.monthly$LAT=with(Data.monthly,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96021),-25,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96022,96023),-26,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97011),-27,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97012,97013),-28,       
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97014,97015),-29,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96010),-33,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96000),-33,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96030,95050,95090,95040),-35,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85030),-34.4586,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85050),-34.2853,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85080),-33.917,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85110),-34.28,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95010),-31.9554,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95020),-32.5167,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95030),-33.3503,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95060),-34.9953,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95070),-34.9731,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95080),-34.9278,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85990),NA,
                ifelse(!is.na(BLOCKX) & BLOCKX==85010,-34.06976,
                ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85130),-34.9278,LAT
                ))))))))))))))))))))))
Data.monthly=Data.monthly%>%
                mutate(LAT=ifelse(is.na(LAT)& PORT=="ALBANY" & Estuary=='YES',-35.1,LAT))
Data.monthly$LAT=ceiling(Data.monthly$LAT)

if(min(Data.monthly$LAT,na.rm=T)<(-40))
{
  plot(1,1,col='transparent',ann=F,axes=F)
  text(1,1,"ERROR IN LATITUDE",col=2,cex=3)
  stop("CHECK LATITUDE")
}
  

Data.monthly$LONG=with(Data.monthly,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96021),113,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96022,96023),113,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97011),113,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97012,97013),113,       
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(97014,97015),113,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96010),115,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96000),115,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(96030,95050,95090,95040),118,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85030),118.8897,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85050),119.4819,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85080),120.050,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85110),115.18,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95010),115.8585,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95020),115.7167,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95030),115.6494,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95060),117.4117,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95070),116.9744,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(95080),116.4489,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85990),NA,
                  ifelse(!is.na(BLOCKX) & BLOCKX==85010,115.1438 ,
                  ifelse(!is.na(BLOCKX) & BLOCKX%in%c(85130),116.4489,LONG
                  ))))))))))))))))))))))
Data.monthly=Data.monthly%>%
        mutate(LONG=ifelse(is.na(LONG)& PORT=="ALBANY" & Estuary=='YES',117.8,LONG))

Data.monthly$LONG=floor(Data.monthly$LONG)

#Reset Same.return after fixing blockx
Data.monthly$Same.return=with(Data.monthly,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
Data.monthly$Same.return.SNo=Data.monthly$Same.return

#add Fishery code if missing
Data.monthly=Data.monthly%>%
  mutate(fishery=ifelse(is.na(fishery) & METHOD%in%c("GN","LL")& LAT <(-26) & 
                                   LAT >=(-33) & LONG<116.5,'WCGL',
                 ifelse(is.na(fishery) & METHOD%in%c("GN","LL")& LAT <(-33) &
                                          LONG<116.5,'SGL1',
                 ifelse(is.na(fishery) & METHOD%in%c("GN","LL")& LAT <(-26) &
                                                 LONG>=116.5,'SGL2',fishery))),
         fishery=case_when(fishery=='*'~NA_character_,
                           fishery=='WCGL'~'WCDGDL',
                           fishery%in%c('C070','OANCGCWC') & LAT <=(-26) & LAT>=(-32)~'WCDGDL',
                           fishery%in%c('WCDGDL','OANCGCWC') & LAT <=(-33)~'JASDGDL',
                           fishery%in%c('WCDGDL') & LAT <=(-31) & LONG>118~'JASDGDL',
                           fishery=='CL02'~'JANS',
                           fishery=='C127'~'WANCS', 
                           fishery=='C051'~'JANS',  
                           fishery%in%c('SGL2','SGL1','SGL')~'JASDGDL',
                           TRUE~fishery))

# A.7. Fix condition
Data.monthly$CONDITN=as.character(Data.monthly$CONDITN)
Data.monthly$CONDITN=with(Data.monthly,ifelse(is.na(SPECIES),"NIL",CONDITN))

#fix some shark records with no fishery
Data.monthly=Data.monthly%>%
              mutate(fishery=case_when(is.na(fishery) & type=='elasmobranchs'& METHOD%in%c('GN','LL') & 
                                         zone%in%c('Closed','Joint','North','West') ~'OANCGCWC',
                                       TRUE~fishery))

# A.8. Add bioregion and zone                                                               
Data.monthly$Bioregion=as.character(with(Data.monthly,ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
             ifelse(LONG<115.5 & LAT<=(-27),"WC",
             ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
             ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))))

Data.monthly$Bioregion=with(Data.monthly,
            ifelse(Bioregion=="SC"& LAT>(-34) & LONG <115.91 ,"WC",Bioregion))

Data.monthly=Data.monthly%>%
                mutate(zone=case_when(
                         fishery%in%(FishCubeCode_Shark.fisheries) & LONG>=116.5 & LAT<=(-26)~"Zone2",
                         fishery%in%(FishCubeCode_Shark.fisheries) & LONG<116.5 & LAT<=(-33)~"Zone1",
                         fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-33) & LAT<=(-26) & LONG<116.5~"West",
                         fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-26) & LONG<114~"Closed",
                         fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-23) & LONG>=114 & LONG<123.75~"North",
                         fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-23) & LONG>=123.75~"Joint",
                                      TRUE~zone))
Data.monthly$zone=as.character(with(Data.monthly,
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LONG>=116.5 & LAT<=(-26),"Zone2",
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LONG<116.5 & LAT<=(-33),"Zone1",
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-26) & LONG<114,"Closed",
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-23) & LONG>=114 & LONG<123.75,"North",
                  ifelse(fishery%in%(FishCubeCode_Shark.fisheries) & LAT>(-23) & LONG>=123.75,"Joint",
                         NA))))))))

#Check spatial blocks and catch assigned to right zone
world <- ne_countries(scale = "medium", returnclass = "sf")
fn.ck.spatial=function(d,BLK.data,Depth.data,BLK.size,add.yr=FALSE,pt.size=3,pt.alpha=1)
{
  Limx=range(d$ln,na.rm=T); Limy=range(d$la,na.rm=T)
  p=ggplot(data = world) +
    geom_sf(color = "black", fill = "darkorange4") +
    coord_sf(xlim =Limx , ylim = Limy, expand = T) +
    xlab("") + ylab("")+
    geom_contour(data = Depth.data%>%filter(V1>=Limx[1] & V1<=Limx[2] & V2>=Limy[1] & V2<=Limy[2]), 
                 aes(x=V1, y=V2, z=V3),
                 breaks=c(-50,-100,-200,-500),linetype="solid",colour="grey70")+
    geom_point(data=d,aes(ln,la,color=zn),size=pt.size,alpha=pt.alpha)+
    geom_text(data=BLK.data%>%mutate(zn='')%>%filter(x>=Limx[1] & x<=Limx[2] & y>=Limy[1] & y<=Limy[2]),
              aes(x,y,label=BlockCode),angle=45,size=BLK.size,alpha=0.5)+
    theme(legend.position = 'top',legend.text = element_text(size=14))+
    scale_x_continuous(breaks=seq(round(Limx)[1],round(Limx)[2]))+
    scale_y_continuous(breaks=seq(round(Limy)[1],round(Limy)[2]))
  if(add.yr) p=p+ggtitle(unique(d$yr))
  print(p)
  
}
do.this=FALSE
if(do.this)
{
  fn.ck.spatial(d=Data.monthly%>%
                  filter(fishery%in%c('SGL','WCGL','SGL1','SGL2','JASDGDL','WCDGDL'))%>%
                  mutate(yr=FINYEAR,
                         zn=zone,
                         ln=LONG,
                         la=LAT,
                         blk=BLOCKX,
                         Unik=Same.return)%>%
                  distinct(Unik,.keep_all = T),
                BLK.data=BLOCK_60%>%
                  mutate(dumi=as.numeric(BlockCode))%>%
                  filter(dumi<85010)%>%
                  rowwise()%>% 
                  mutate(x = mean(c(NorthWestPointGPSLongitude, SouthEastPointGPSLongitude)),
                         y = mean(c(NorthWestPointGPSLatitude, SouthEastPointGPSLatitude))),
                Depth.data=Bathymetry,
                BLK.size=5)
  ggsave(handl_OneDrive("Data/Catch and Effort/Check_these_vars/Monthly/Map_blocks_monthly.tiff"),
         width = 14,height = 10,compression = "lzw")
  
}

# A.9. Create Monthly effort data set   
Data.monthly$NETLEN=with(Data.monthly,ifelse(METHOD%in%c("LL","HL","DL")&NETLEN>0,NA,NETLEN))
Effort.vars=c("FDAYS","BDAYS","HOURS","HOOKS","SHOTS","NETLEN")
Effort.monthly=Data.monthly[,match(c("FINYEAR","YEAR","MONTH","zone","VESSEL","METHOD","BLOCKX",Effort.vars,
                                      "YEAR.c","Same.return","LAT","LONG","fishery"),names(Data.monthly))]
Effort.monthly=Effort.monthly%>%rename(FisheryCode=fishery)

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

#Remove rows will all NA records
Data.monthly=Data.monthly%>%filter(!(is.na(YEAR) & is.na(MONTH) & FINYEAR=='' & is.na(BLOCKX)))

#Get data for Matt Koopman Fishery Review 2023    
if(get.Mats.data.2023)
{
  Mat.monthly=Data.monthly%>%
    dplyr::select(FINYEAR,YEAR,MONTH,FDAYS,METHOD,BLOCKX,BDAYS,HOURS,HOOKS,SHOTS,NETLEN,
                  SPECIES,SNAME,LIVEWT,fishery,PORT,type,VESSEL,LAT,LONG,zone)
  Vesl.id=rbind(Mat.monthly%>%
                  distinct(VESSEL),
                Data.daily%>%
                  distinct(vessel)%>%
                  rename_with(toupper))%>%
    distinct(VESSEL)%>%
    mutate(VESSEL.d=paste('Vessel',row_number()))
  Mat.monthly=Mat.monthly%>%
    left_join(Vesl.id,by='VESSEL')%>%
    dplyr::select(-VESSEL)%>%
    filter(LAT<=(-26))%>%
    mutate(Shark.fishery=case_when(zone=='Joint' & fishery%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                   zone=='North' & fishery%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                   zone=='Closed' & fishery%in%c('OT')~'WANCS',
                                   zone=='West' & fishery%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                   is.na(zone) & fishery=='WCGL'~'WCDGDL',
                                   zone=='Zone1' & fishery=='WCGL'~'JASDGDL',
                                   fishery%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                   zone%in%c('Zone1','Zone2') & fishery=='*'~'OASC',
                                   zone%in%c('Zone1','Zone2') & fishery=="WL"~'OASC',
                                   zone%in%c('West') & fishery=='*'~'OANCGCWC',
                                   zone%in%c('West','Zone1') & fishery=='OT'~'OANCGCWC',
                                   METHOD=='HL' & is.na(fishery) ~ 'Handline shark fishery',
                                   METHOD=='HL' & fishery=='SBS' ~ 'Shark Bay Pink Snapper fishery',
                                   METHOD=='DL' & is.na(fishery) ~ 'Dropline shark fishery', 
                                   METHOD=='LL' & fishery=='CSLP' ~ 'estuary fishery',
                                   METHOD=='GN' & fishery%in%c('C019','C066','CSFN','MBC','SCE','WCE') ~ 'estuary fishery',
                                   TRUE~"other fishery"))%>%
    filter(Shark.fishery%in%c('Dropline shark fishery','Handline shark fishery','JASDGDL','OANCGCWC','OASC','WCDGDL'))%>%
    rename_with(toupper)
}

# Add Factor if not available but landwt and conditn available
Data.monthly=Data.monthly%>%
                mutate(Factor=ifelse(grepl('Skates',SNAME) & is.na(Factor) & CONDITN%in%c('HD','WD','WF'),1.59,Factor),
                       LIVEWT=ifelse(grepl('Skates',SNAME) & is.na(LIVEWT) & CONDITN%in%c('HD','WD','WF'),LANDWT*Factor,
                                     LIVEWT))
#Check for 0 or NA livewt records
a=Data.monthly%>%filter(LANDWT==0 | LIVEWT==0 | is.na(LANDWT) | is.na(LIVEWT))%>%
  dplyr::select(CONDITN,Factor,SNAME,fishery,LIVEWT,LANDWT)%>%
  filter(is.na(LIVEWT) | LIVEWT==0)%>%
  distinct(SNAME,CONDITN,Factor,LANDWT,LIVEWT)
if(nrow(a)>0)
{
  par(bg=2)
  plot.new()
  mtext(paste("there are ",nrow(a),"records"),3,cex=3,col="white")
  mtext("with no weight",3,-2,cex=3,col="white")
  par(bg="white")
  
}

#SECTION B. ---- DATA MANIPULATION - DAILY LOGBOOKS ----

#Check catch composition deeper than 100 m
  #note: only using Daily logbooks, Monthly return have very little data (mostly 'shark other')
do.ktch.deep=FALSE
if(do.ktch.deep)
{
  Deep.ktch=Data.daily%>%
    filter(method%in%c('GN','LL'))%>%
    mutate(LatDeg=ifelse(blockx<35600 & LatDeg==0,as.numeric(substr(blockx,1,2)),LatDeg),
           LatMin=ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin),
           LongDeg=ifelse(blockx<35600 & LongDeg==0,100+as.numeric(substr(blockx,3,4)),LongDeg),
           LongMin=ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin),
           LAT=LatDeg+(LatMin/60),
           LON=LongDeg+(LongMin/60),
           LAT=-abs(LAT))%>%
    filter(depthMax>100 & LAT<=(-26) & !is.na(sname1))%>%
    mutate(zone=ifelse(zone=='*','west',paste('zone',zone)),
           Depth.interval=10*round(depthMax/10),
           Depth.interval=ifelse(Depth.interval>150,'>150',Depth.interval),
           Depth.interval=paste(Depth.interval,'m'),
           Depth.interval=factor(Depth.interval,levels=c('100 m','110 m','120 m',
                                                         '130 m','140 m','150 m',
                                                         '>150 m')))%>%
    dplyr::select(species,RSCommonName,livewt,method,zone,depthMax,Depth.interval)%>%
    filter(!is.na(livewt))
  
  dummy=Deep.ktch%>%
    group_by(RSCommonName)%>%
    summarise(livewt=sum(livewt)/1000)%>%
    arrange(-livewt)%>%
    mutate(CumKtch=100*cumsum(livewt)/sum(livewt))%>%
    filter(CumKtch>95)%>%pull(RSCommonName)
  
  Deep.ktch1=Deep.ktch%>%
    mutate(RSCommonName1=ifelse(RSCommonName%in%dummy,"Other",RSCommonName),
           RSCommonName1=ifelse(RSCommonName1=="Gulper sharks, Sleeper Sharks & Dogfishes",'Dogfishes',RSCommonName1),
           Zone.method=paste(zone,method,sep='-'))%>%
    group_by(Zone.method,method,zone,Depth.interval,RSCommonName1)%>%
    summarise(livewt=sum(livewt)/1000)
  
  tiff(file=handl_OneDrive('Analyses/Catch and effort/Data_Resquests/Catch_comp_deeper100_Daily.tiff') ,
       width=2400,height=2000,units="px",res=300,compression="lzw")
  Deep.ktch1%>%
    ggplot(aes(x = RSCommonName1, y = livewt))+
    geom_col(aes(fill = Zone.method), width = 0.7)+
    facet_wrap(~Depth.interval,scales='free_x')+
    coord_flip()+
    ylab("Total catch (tonnes)")+xlab('')+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text=element_text(size=12),
          axis.text=element_text(size=7),
          axis.title=element_text(size=14))+
    scale_fill_manual("legend", 
                      values = c("west-GN" = "brown1", "west-LL" = "brown4",
                                 "zone 1-GN" = "dodgerblue3",
                                 "zone 2-GN" = "green3","zone 2-LL" = "darkgreen"))
  dev.off()
  
}

#Set hammerhead species to the original report to overwrite FishCube's reapportioning     
if(reset.hammerhead.to.reported)
{
  species.spec.reporting.year=2023
  all.hammers=Data.daily%>%
              filter(species%in%c(19000,HammerheadSpecies))%>%
              mutate(sname1=ifelse(species%in%HammerheadSpecies & year<species.spec.reporting.year,'Hammerhead Sharks',sname1),
                     RSCommonName=ifelse(species%in%HammerheadSpecies & year<species.spec.reporting.year,'Hammerhead Sharks',RSCommonName),
                     RSSpeciesCode=ifelse(species%in%HammerheadSpecies & year<species.spec.reporting.year,37019000,RSSpeciesCode),
                     RSSpeciesId=ifelse(species%in%HammerheadSpecies & year<species.spec.reporting.year,90,RSSpeciesId),
                     species=ifelse(species%in%HammerheadSpecies & year<species.spec.reporting.year,19000,species))
  Nms.hmrs=names(all.hammers)
  Nms.hmrs=Nms.hmrs[-match(c('nfish','livewt','TripLandedWeight'),Nms.hmrs)]
  all.hammers=all.hammers%>%
            group_by_at(Nms.hmrs)%>%
            summarise(nfish=sum(nfish,na.rm=T),
                      livewt=sum(livewt,na.rm = T),
                      TripLandedWeight=sum(TripLandedWeight,na.rm = T))
  Data.daily=Data.daily%>%
          filter(!species%in%c(19000,HammerheadSpecies))
  
  Data.daily=rbind(Data.daily,all.hammers%>%relocate(names(Data.daily)))
}

#simple financial assessment
if(do.financial.ass=="YES")
{
  b=aggregate(landwt~vessel+species,subset(Data.daily,finyear==Current.yr,select=c(landwt,vessel,species)),sum,na.rm=T)
  bb=merge(b,subset(PRICES,select=c(SPECIES,uv1516)),by.x="species",by.y="SPECIES",all.x=T)
  bb$Revenue_annual=bb$livewt*bb$uv1516
  dd=aggregate(cbind(Revenue_annual,landwt)~vessel,bb,sum,na.rm=T)
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
  A=A[,match(c("vessel","BoatName","MastersName","crew","fdays.c" ,"landwt","Revenue_annual","Revenue_per_fishing_day","Costs"),names(A))]
  names(A)[match(c("crew","fdays.c","landwt","Revenue_annual","Revenue_per_fishing_day"),names(A))]=
    c("crew_average","fishing_days_annual","livewt_annual(kg)","Revenue_annual (AUD$)","Revenue_per_fishing_day (AUD$)")
  write.csv(A,paste(handl_OneDrive("Analyses/Catch and effort/Annual.revenue."),Current.yr,".csv",sep=""),row.names=F)
  
}

#Extract individual landwt not just trip landed weight
Data.daily=Data.daily%>%
            mutate(landwt=livewt/(factor))

#Create backup file   
Data.daily.original=Data.daily
names(Data.daily.original)[match(c("finyear","month","vessel","method","bdays","hours","hooks",
               "shots","netlen","LatDeg"),names(Data.daily.original))]=c("FINYEAR","MONTH","VESSEL","METHOD",
                     "BDAYS","HOURS","HOOKS","SHOTS","NETLEN","LAT")
Data.daily.original$LAT=-abs(Data.daily.original$LAT)
Data.daily.original$Same.return.SNo=with(Data.daily.original,paste(SNo,DSNo,TSNo))


#Fix boat and masters names
Data.daily=Data.daily%>%
  mutate(BoatName=tolower(BoatName),
         BoatName=ifelse(BoatName=="kabralee 11","kabralee ii",
                  ifelse(BoatName%in%c("st gerard m","st. gerard m","st. gerard m"),"st gerard m",
                  ifelse(BoatName=="san marco","san margo",
                  ifelse(BoatName=="tara-marie","tara marie",
                  ifelse(BoatName=="barbarosa ll","barbarosa ii",
                  ifelse(BoatName%in%c("chivers regal 11","chivers regal 2"),"chivers regal ii",
                  ifelse(BoatName%in%c("discovery  111","discovery 111"),"discovery iii",
                  ifelse(BoatName%in%c("dorreen","dorren"),"doreen",
                  ifelse(BoatName%in%c("elizabeth maria 11","elizabith maria 11","elizabeht maria 11"),"elizabeth maria ii",
                  ifelse(BoatName%in%c("falcon 11","falcon 2","falcon ll"),"falcon ii",
                  ifelse(BoatName=="fish tales","fishtales",
                  ifelse(BoatName=="giuliano 2","giuliano ii",
                  ifelse(BoatName=="lone hand","lonehand",
                  ifelse(BoatName=="maniki 2","maniki ii",
                  ifelse(BoatName%in%c("planjak 11","planjak 2"),"planjak ii",
                  ifelse(BoatName=="sea venture 11","sea venture ii",
                  ifelse(BoatName=="steve-mayree d","steve mayree d",
                  ifelse(BoatName=="sveti-nikola","sveti nikola",
                  ifelse(BoatName=="tracey-lea","tracey lea",
                  ifelse(BoatName=="carado","corado",
                         ifelse(BoatName=="st gerard","st gerard m",
                  ifelse(BoatName=="catch fillet release","catch fillet & release",
                         BoatName)))))))))))))))))))))),
         MastersName=tolower(MastersName),
         MastersName=case_when(MastersName%in% c("a. joy","andrew joy","andrew f joy")~"joy, andrew francis",
                               MastersName%in% c("j.e.robb","j.e. robb","j e robb","j. e. robb","j. robb","james robb",
                                                "j.robb","james e robb","james e. robb","james edward robb")~  "robb,james",
                               MastersName%in%c("marcus branderhorst","marcus allen branderhorst") ~ "branderhorst, marcus allen",
                               MastersName%in% c("steve buckeridge","stephen buckeridge",
                                                "buckeridge, steve","stephen grant buckeridge")~ "buckeridge, stephen grant",
                               MastersName%in% c("ryan bradley","ryan james bradley",'ryan, bradley')~ "bradley, ryan",
                               MastersName%in% c("chris bradley","christopher anthony bradley")~ "bradley, chris",
                               MastersName== "graeme edward sell"~ "sell, graeme edward",
                               MastersName%in% c("john patrick richardson","john richardson","j.richardson")~ "richardson,john",
                               MastersName%in% c("paul douglas murch","paul murch","p. murch")~ "murch, paul douglas",
                               MastersName== "stephen charles mcwhirter"~ "mcwhirter, stephen charles",
                               MastersName== "storm mansted"~ "mansted, storm",
                               MastersName%in% c("anthony  mansted","anthony david mansted","anthony mansted")~ "mansted, anthony",
                               MastersName== "mason thomas"~ "thomas, mason",
                               MastersName%in% c("james tindall","james stuart tindall")~ "tindall, james stuart",
                               MastersName%in% c("m.tonkin","m. tonkin","michael tonkin")~ "tonkin, michael",
                               MastersName%in% c("jayson lindsay & scott farrant","jayson lindsay / scott farrant",
                                                "jayson lyndsay  / scott farrant")~ "lindsay, jayson / scott farrant",
                               MastersName%in% c("scott  farrant")~ "scott farrant",
                               MastersName%in% c("anthony james cooke","anthony cooke") ~ "cooke, anthony",
                               MastersName%in% c("jeffrey frank cooke","jeffrey cooke","jeff cooke")~ "cooke, jeffrey",
                               MastersName%in% c("c. black","chris black","christopher black","christopher barry black")~  "black, christopher barry",
                               MastersName%in% c("tim  goodall","tim goodall")~  "goodall, tim",
                               MastersName%in% c("jamie walter thornton","j.thornton",
                                                 "jw thornton","j.w.thornton")~  "thornton, j.w.",
                               MastersName%in% c("n. triantafyllou","neoclis triantafyllou")~ "triantafyllou, neoclis",
                               MastersName%in% c("d. hawkins","david joseph hawkins")~ "hawkins, david joseph" ,
                               MastersName%in% c("darcy madgen","darcy james madgen")~  "madgen, darcy james",
                               MastersName%in% c("n e soulos","n.e. soulos","emanuel soulos","emanuel nicholas soulos","n.e soulos",
                                                "n.e.soulas","n.e.soulos","ne soulos")~ "soulos, emanuel nicholas" ,
                               MastersName%in% c("g.m. sharp","greg sharp","gregory mark sharp")~  "sharp, gregory mark",
                               MastersName%in% c("glen brodi","glen brodie","glen nicholas brodie")~ "brodie, glen",
                               MastersName%in% c("geoffrey foster campbell","geoff campbell","geoffrey campbell")~ "campbell,geoffrey",
                               MastersName%in% c("g. whetstone","greg whetstone","whetstone, gregory neil",
                                                "g whetstone","gregory neil whetstone")~ "whetstone,greg" ,
                               MastersName%in% c("brian scimone","brian gregory scimone",
                                                "scimone, brian gregory")~ "scimone, brian" ,
                               MastersName%in% c("jason & brian scimone","brian & jason  scimone","brian & jason scimone",
                                                "brian / jason scimone")~ "brian / jason scimone" ,
                               MastersName%in% c("c. henderson","chris henderson","chris  henderson",
                                                 "chris henderson")~ "henderson, christian william",
                               MastersName%in% c("andrew  joy","andrew francis joy","andrew f joy")~ "joy, andrew f",
                               MastersName%in% c("smythe, john lawrence","j smythe","john smythe")~ "smythe, john",
                               MastersName%in% c("m.bub","mark robert bubb","mark robert bubb(m.tonkin)")~ "bubb, mark robert",
                               MastersName%in% c("manual karaterpos","manuel karaterpos","m.karaterpos")~ "karaterpos, manuel",
                               MastersName%in% c("warrilow, peter charles","peter charles warrilow",
                                                "p. warrilow","peter  warrilow","peter warrilow")~ "warrilow, peter",
                               MastersName== "phil  toumazos"~ "phil toumazos",
                               MastersName%in% c("parker, roger joseph","roger joseph parker","roger parker")~ "parker, roger",
                               TRUE~MastersName))
#export list of boats and skippers
write.csv(Data.daily%>%distinct(vessel,BoatName,MastersName)%>%arrange(MastersName),
          'List_boats_skippers.csv',row.names = F)
List.of.daily.vessel.boat=Data.daily%>%
  mutate(Same.return.SNo=paste(SNo,DSNo,TSNo))%>%
  distinct(vessel,BoatName,MastersName,Same.return.SNo) 

#Fishers operating in recent years  
recent.yrs=as.numeric(substr(Current.yr,1,4))
recent.yrs=paste(seq((recent.yrs-2),recent.yrs),substr(seq((recent.yrs-1),recent.yrs+1),3,4),sep='-')
Recent.fishers.n.shots=Data.daily%>%
                          filter(finyear%in%recent.yrs)%>%
                          distinct(date,DSNo,TSNo,SNo,.keep_all = T)%>%
                          group_by(MastersName,BoatName,vessel,finyear)%>%
                          tally()%>%
                          spread(finyear,n,fill=0)%>%
                          arrange(BoatName,MastersName)

#note: a unique shot (session) is the combination of Sno (shot num), DSNo (daily sheet num) and     
#       TSNo (trip sheet num),i.e. the variable "Session ID"
Session.vars=c("SNo","DSNo","TSNo")
FishCube.vars=c("RSSpeciesId","RSCommonName","Block","LatFC","LongFC","blockxFC","port")
This=c("finyear","year","month","day","date","zone","depthMax","vessel","fdays","method","blockx","block10","bdays",
       "hours","hooks","shots","nlines","netlen","species","sname1","nfish","livewt","conditn","landwt",
       "factor","LatDeg","LongDeg","LatMin","LongMin",Session.vars,"Bioregion","fishery",FishCube.vars,
       'type','RSSpeciesCode')

#create file for checking mesh size
Mesh.size=subset(Data.daily,method=="GN",select=c(DSNo,TSNo,SNo,vessel,date,finyear,method,
                            netlen,hours,shots,mshigh,mslow,LatDeg,LongDeg,zone))

#Set RSCommonName== 'Nil Fish Caught' from NA livewt to 0 livewt
Data.daily=Data.daily%>%
  mutate(nfish=ifelse(is.na(nfish) & RSCommonName=='Nil Fish Caught',0,nfish),
         landwt=ifelse(is.na(landwt) & RSCommonName=='Nil Fish Caught',0,landwt),
         livewt=ifelse(is.na(livewt) & RSCommonName=='Nil Fish Caught',0,livewt))

  #Select which file to use
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

Data.daily=Data.daily[,match(c(This,"flagtype"),names(Data.daily))]
Data.daily.FC.NA.sp=subset(Data.daily,species==99999 | is.na(species))

Data.daily.FC.2005_06=subset(Data.daily,(finyear=="2005-06" & fishery%in%c('*','JASDGDL','WCDGDL')))
Data.daily.FC.2005_06=subset(Data.daily.FC.2005_06,!(species==99999 | is.na(species)))

Data.daily=subset(Data.daily,!(finyear=="2005-06" & fishery%in%c('*','JASDGDL','WCDGDL')))

#Daily data fixes
#Fix some records identified as incorrect in CATCH INSPECTIONS and EFFORT INSPECTIONS  
#note: the amended values were provided by data entry girls

id=subset(Data.daily,finyear=="2011-12" & TSNo=="TDGLF8000737" & 
            DSNo=="TDGLF8000734" & sname1=="Blue Morwong")
if(nrow(id)>=1)Data.daily=subset(Data.daily,!(finyear=="2011-12" & TSNo=="TDGLF8000737" & 
                                                DSNo=="TDGLF8000734" & sname1=="Blue Morwong"))

Data.daily$livewt=with(Data.daily,
                       ifelse(finyear=="2011-12" & TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & 
                                sname1=="Blue Morwong"& livewt<700,5.04,
                       ifelse(finyear=="2011-12" & TSNo=="TDGLF6007702" & DSNo=="TDGLF6007702" & 
                                sname1=="Blue Morwong"& livewt>700,10.06,
                       ifelse(finyear=="2011-12" & TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & 
                                sname1=="Blue Morwong"& livewt<200,13,
                       ifelse(finyear=="2011-12" & TSNo=="TDGLF8004653" & DSNo=="TDGLF8004653" & 
                                sname1=="Blue Morwong"& livewt>200,26,
                              livewt)))))

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
                        ifelse(species%in%c(17006,18003,19000,HammerheadSpecies)& TSNo=="TDGLF8009537" & flagtype==0,
                                 1.59*7.911765*nfish*807/933.5883,
                        ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009503" & flagtype==0,
                                 1.59*7.090909*nfish,
                        livewt)))))))))))))))

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
                        ifelse(species%in%c(17006,18003,19000,HammerheadSpecies)& TSNo=="TDGLF8009537" & flagtype==0,
                               7.911765*nfish*807/933.5883,
                        ifelse(species%in%c(18003,18007)& TSNo=="TDGLF8009503" & flagtype==0,
                               7.090909*nfish,
                        landwt)))))))))))))))


Data.daily=Data.daily[,-match("flagtype",names(Data.daily))]

Data.daily$nfish=with(Data.daily,
                      ifelse(finyear=="2011-12" & DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & 
                               sname1=="Whiskery Shark"& nfish>3,52,
                      ifelse(finyear=="2011-12" & DSNo=="TDGLF8003608" & TSNo=="TDGLF8003611" & 
                              sname1=="Whiskery Shark"& nfish<3,12,
                      nfish)))

Data.daily$netlen=with(Data.daily,
                       ifelse(finyear=="2011-12" & vessel=="B 067"& netlen==650,6500,
                       ifelse(finyear=="2011-12" & vessel=="B 142"& netlen==540,5400,
                       ifelse(finyear=="2011-12" & vessel=="E 009"& netlen==420,4320,
                       ifelse(finyear=="2011-12" & vessel=="E 009"& netlen==400,4000,
                       ifelse(finyear=="2011-12" & vessel=="E 035"& netlen==380,3800,
                       ifelse(finyear=="2010-11" & vessel=="G 297"& netlen==7320,4320,
                       ifelse(finyear=="2011-12" & vessel=="F 417"& netlen==6000,600,
                       ifelse(finyear%in%c("2010-11","2013-14") & vessel=="F 417"& netlen==7500,750,
                       ifelse(finyear=="2013-14" & vessel=="F 417"& netlen==6500,650,
                       netlen))))))))))

Data.daily$shots=with(Data.daily,
                      ifelse(finyear%in%c("2009-10","2010-11") & vessel=="F 541"& shots==6,1,
                      ifelse(finyear=="2011-12" & vessel=="F 541"& shots==5,1,
                      ifelse(finyear=="2010-11" & vessel=="E 056"& shots==10,1,
                      ifelse(finyear=='2020-21' & vessel=="E 007"& (shots==0| is.na(shots)),1,
                             shots)))))

Data.daily$nfish=with(Data.daily,
                       ifelse(finyear=="2012-13" & species==377004 & vessel=="G 297" & landwt == 92.57,7,nfish))

Data.daily$livewt=with(Data.daily,
                        ifelse(finyear=="2012-13" & species==377004 & vessel=="E 035" & landwt == 49,6.7,
                        ifelse(finyear=="2012-13" & species==377004 & vessel=="G 297" & landwt == 92.57,216,
                        ifelse(finyear=="2012-13" & species==17003 & vessel=="E 059" & landwt == 1769,17,
                        livewt))))


Data.daily$netlen=with (Data.daily,
                        ifelse(finyear=="2012-13" & DSNo=="TDGLF8002239" &  TSNo=="TDGLF8002240" & netlen==165,3500,
                        ifelse(finyear=="2012-13" & TSNo=="TDGLF8002261" & netlen==300,3000,
                        ifelse(finyear=="2012-13" & TSNo%in%c("TDGLF8008946","TDGLF8008924") & netlen==400,4000,
                        netlen))))


  #Remove incomplete years 
#FINYrs=c(sort(as.character(unique(Data.monthly$FINYEAR))),sort(as.character(unique(Data.daily$finyear))))
FINYrs=sort(unique(c(as.character(unique(Data.monthly$FINYEAR))),sort(as.character(unique(Data.daily$finyear)))))
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

#Check if any fisher reported just fins 
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
Data.daily=subset(Data.daily,!(species%in%c(22997,22998))) #remove fins & livers to avoid double counting when scaling carcass to whole 
if(nrow(Data.daily.incomplete)>0)Data.daily.incomplete=subset(Data.daily.incomplete,!(species%in%c(22997,22998)))
if(nrow(Data.daily.FC.2005_06)>0)Data.daily.FC.2005_06=subset(Data.daily.FC.2005_06,!(species%in%c(22997,22998)))


  #fix typos
Data.daily$netlen=with(Data.daily,
      ifelse(vessel=="G 297" & finyear=="2014-15" & netlen==400,4000,
      ifelse(vessel=="E 075" & finyear=="2014-15" & netlen==200,2000,netlen)))

  #Create same return                                                                          
Data.daily$Same.return=with(Data.daily,paste(finyear,month,vessel,method,blockx)) 

  #Create unique session 
Data.daily$Same.return.SNo=with(Data.daily,paste(SNo,DSNo,TSNo))
Data.daily$ID=with(Data.daily,paste(DSNo,SNo,"_",sep=""))                                 


  #Create data set identifier
Data.daily$TYPE.DATA="daily"

  #Fix NA lat and long degrees but availabe FC
Data.daily=Data.daily%>%
            mutate(LatFC=as.numeric(LatFC),
                   LongFC=as.numeric(LongFC),
                   LatDeg=ifelse(is.na(LatDeg) & !is.na(LatFC),-floor(abs(LatFC)),LatDeg),
                   LatMin=ifelse(is.na(LatMin) & !is.na(LatFC),(LatFC%%1)*60,LatMin),
                   LongDeg=ifelse(is.na(LongDeg) & !is.na(LongFC),floor(LongFC),LongDeg),
                   LongMin=ifelse(is.na(LongMin) & !is.na(LongFC),(LongFC%%1)*60,LongMin))

  #Fix estuaries lats and longs 
Data.daily$LongDeg=with(Data.daily,
                        ifelse(!is.na(blockx) & blockx%in%c(96021) & LatDeg>60,113,
                        ifelse(!is.na(blockx) & blockx%in%c(96022,96023) & LatDeg,113,
                        ifelse(!is.na(blockx) & blockx%in%c(97011) & LatDeg,113,
                        ifelse(!is.na(blockx) & blockx%in%c(97012,97013) & LatDeg,113,       
                        ifelse(!is.na(blockx) & blockx%in%c(97014,97015) & LatDeg,113,
                        ifelse(!is.na(blockx) & blockx%in%c(96010) & LatDeg,115,
                        ifelse(!is.na(blockx) & blockx%in%c(96000) & LatDeg,115,
                        ifelse(!is.na(blockx) & blockx%in%c(96030,95050,95090,95040) & LatDeg,118,
                        ifelse(!is.na(blockx) & blockx%in%c(85030) & LatDeg,118.8897,
                        ifelse(!is.na(blockx) & blockx%in%c(85050) & LatDeg,119.4819,
                        ifelse(!is.na(blockx) & blockx%in%c(85080) & LatDeg,120.050,
                        ifelse(!is.na(blockx) & blockx%in%c(95010) & LatDeg,115.8585,
                        ifelse(!is.na(blockx) & blockx%in%c(95020) & LatDeg,115.7167,
                        ifelse(!is.na(blockx) & blockx%in%c(95030) & LatDeg,115.6494,
                        ifelse(!is.na(blockx) & blockx%in%c(95060) & LatDeg,117.4117,
                        ifelse(!is.na(blockx) & blockx%in%c(95070) & LatDeg,116.9744,
                        ifelse(!is.na(blockx) & blockx%in%c(95080) & LatDeg,116.4489,
                        ifelse(!is.na(blockx) & blockx%in%c(85990) & LatDeg,NA,
                        ifelse(!is.na(blockx) & blockx==85010 & LatDeg,115.1438,
                        ifelse(!is.na(blockx) & blockx%in%c(85130) & LatDeg,116.4489,
                               LongDeg)))))))))))))))))))))
Data.daily$LongDeg=floor(Data.daily$LongDeg) 

Data.daily$LatDeg=with(Data.daily,  
                  ifelse(!is.na(blockx) & blockx%in%c(96021) & LatDeg>60,25,
                  ifelse(!is.na(blockx) & blockx%in%c(96022,96023) & LatDeg>60,26,
                  ifelse(!is.na(blockx) & blockx%in%c(97011) & LatDeg>60,27,
                  ifelse(!is.na(blockx) & blockx%in%c(97012,97013) & LatDeg>60,28,       
                  ifelse(!is.na(blockx) & blockx%in%c(97014,97015) & LatDeg>60,29,
                  ifelse(!is.na(blockx) & blockx%in%c(96010) & LatDeg>60,33,
                  ifelse(!is.na(blockx) & blockx%in%c(96000) & LatDeg>60,33,
                  ifelse(!is.na(blockx) & blockx%in%c(96030,95050,95090,95040) & LatDeg>60,35,
                  ifelse(!is.na(blockx) & blockx%in%c(85030) & LatDeg>60,34.4586,
                  ifelse(!is.na(blockx) & blockx%in%c(85050) & LatDeg>60,34.2853,
                  ifelse(!is.na(blockx) & blockx%in%c(85080) & LatDeg>60,33.917,
                  ifelse(!is.na(blockx) & blockx%in%c(95010) & LatDeg>60,31.9554,
                  ifelse(!is.na(blockx) & blockx%in%c(95020) & LatDeg>60,32.5167,
                  ifelse(!is.na(blockx) & blockx%in%c(95030) & LatDeg>60,33.3503,
                  ifelse(!is.na(blockx) & blockx%in%c(95060) & LatDeg>60,34.9953,
                  ifelse(!is.na(blockx) & blockx%in%c(95070) & LatDeg>60,34.9731,
                  ifelse(!is.na(blockx) & blockx%in%c(95080) & LatDeg>60,34.9278,
                  ifelse(!is.na(blockx) & blockx%in%c(85990) & LatDeg>60,NA,
                  ifelse(!is.na(blockx) & blockx==85010 & LatDeg>60,34.06976,
                  ifelse(!is.na(blockx) & blockx%in%c(85130) & LatDeg>60,34.9278,
                         LatDeg)))))))))))))))))))))
Data.daily$LatDeg=-floor(abs(Data.daily$LatDeg))

  #Identify estuaries
Data.daily$Estuary=with(Data.daily,ifelse(blockx%in%Estuaries,"YES","NO"))

# fix blocks 
#note: this creates a dummy blockx just for extracting lat and long, it gets reset below in '#reset blockx to original'
Data.daily$blockx=with(Data.daily,ifelse(blockx%in%c(96021),25120,     #Shark Bay
                  ifelse(blockx%in%c(96022,96023),26131,
                  ifelse(blockx%in%c(97011),27132,                        #Abrolhos
                  ifelse(blockx%in%c(97012,97013),28132,       
                  ifelse(blockx%in%c(97014,97015),29132,
                  ifelse(blockx%in%c(96010),33151,                        #Geographe Bay
                  ifelse(blockx%in%c(96000),32150,                        #Cockburn sound
                  ifelse(blockx%in%c(96030),35181,blockx)))))))))         # King George sound

#Extract reported ray catches in West and South coast estuaries
do.ray.ktch.estuary=FALSE
if(do.ray.ktch.estuary)
{
  Data.daily%>%
    filter(blockx%in%c(96010,96000,96030) & species%in%Ray.species)%>%
    mutate(Est=case_when(blockx==96000~'Cockburn sound',
                         blockx==96010~'Geographe Bay',
                         blockx==96030~'King George Sound'))%>%
    group_by(finyear,RSCommonName,Est)%>%
    summarise(Tot=sum(livewt)/1000)%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    filter(!RSCommonName=='Guitarfishes')%>%
    ggplot(aes(year,Tot,color=RSCommonName))+
    geom_line(size=1.5)+
    facet_wrap(~Est)+
    theme(strip.text.x = element_text(size = 14),
          axis.text=element_text(size=12),
          legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          legend.text = element_text(size = 14),
          title=element_text(size=16))+
    xlab("Financial year") +ylab("Catch (tonnes)")
  ggsave("C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/Ray catches.tiff", width = 14,height = 10, dpi = 300, compression = "lzw")
  
}

#Fix 0 lat degrees or long degrees
Data.daily$LatDeg=with(Data.daily,ifelse(blockx<35600 & LatDeg==0,as.numeric(substr(blockx,1,2)),LatDeg))  
Data.daily$LongDeg=with(Data.daily,ifelse(blockx<35600 & LongDeg==0,100+as.numeric(substr(blockx,3,4)),LongDeg))

Data.daily$LatDeg=with(Data.daily,ifelse(vessel=="B 067" & blockx%in%c(35610,35630,35640,35650) & LatDeg==0,-35,LatDeg))  
Data.daily$LongDeg=with(Data.daily,ifelse(vessel=="B 067" & blockx%in%c(35610,35630,35640,35650) & LongDeg==0,116,LongDeg))

Data.daily=Data.daily%>%
  mutate(LatDeg=ifelse(is.na(LatDeg) & fishery=="PFT",-19.5,LatDeg),
         LongDeg=ifelse(is.na(LongDeg) & fishery=="PFT",118,LongDeg))

#Add lat and long if NA for calculating zone (allocate to top left corner)
Data.daily=Data.daily%>%
  mutate(LatDeg=ifelse(is.na(LatDeg) & !is.na(LatFC), as.integer(LatFC), LatDeg),
         LatMin=ifelse(is.na(LatMin) & !is.na(LatFC), 60*(as.numeric(LatFC)-floor(as.numeric(LatFC))), LatMin),
         LongDeg=ifelse(is.na(LongDeg) & !is.na(LongFC), as.integer(LongFC) , LongDeg),
         LongMin=ifelse(is.na(LongMin) & !is.na(LongFC), 60*(as.numeric(LongFC)-floor(as.numeric(LongFC))), LongMin))

Data.daily$LatDeg=with(Data.daily,ifelse(is.na(LatDeg)& !is.na(block10),as.numeric(substr(block10,1,2)),LatDeg))  
Data.daily$LatMin=with(Data.daily,ifelse(is.na(LatMin)& !is.na(block10),10*as.numeric(substr(block10,3,3)),LatMin))  
Data.daily$LongDeg=with(Data.daily,ifelse(is.na(LongDeg)& !is.na(block10),100+as.numeric(substr(block10,4,5)),LongDeg))
Data.daily$LongMin=with(Data.daily,ifelse(is.na(LongMin)& !is.na(block10),10*as.numeric(substr(block10,6,6)),LongMin))

Data.daily$LatMin=with(Data.daily,ifelse(is.na(LatMin) & is.na(block10),NA,LatMin))
Data.daily$LongMin=with(Data.daily,ifelse(is.na(LongMin) & is.na(block10),NA,LongMin))

Data.daily$blockx=Data.daily$blockxFC   #reset blockx to original
Data.daily$LatDeg=-abs(Data.daily$LatDeg)

#Reset same return after fixing blocks                                                                          
Data.daily$Same.return=with(Data.daily,paste(finyear,month,vessel,method,blockx)) 

  #fix some fishery codes 
Data.daily=Data.daily%>%
            mutate(fishery=case_when(fishery=='JASDGDL' & LatDeg>(-33) & LongDeg<118~'WCDGDL',
                                     TRUE~fishery))
  #Fix zones                                                                           
Data.daily=Data.daily%>%
  mutate(LatDeg=abs(LatDeg),
         zone=as.character(zone),
         zone=ifelse(is.na(LatDeg) | is.na(LatMin) | is.na(LongDeg) | is.na(LongMin),zone,
             ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
              ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
              ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
              ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
              ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
              ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",
              NA))))))),
         LatDeg=-abs(LatDeg))
Data.daily$zone=with(Data.daily,
                ifelse(zone=='2',"Zone2",
                ifelse(zone=='1',"Zone1",
                ifelse(zone=='*' & as.numeric(LatFC)>(-33),"West",
                ifelse(zone=='*' & as.numeric(LatFC)<=(-33),"Zone1",
                ifelse(zone=='*' & as.numeric(LatFC)<=(-33) & as.numeric(LongFC)>=116.5,"Zone2",
                zone))))))

#only allocate zone to shark fishery codes
Data.daily=Data.daily%>%
              mutate(zone=ifelse(!fishery%in%FishCubeCode_Shark.fisheries,NA,zone))
                     
  #Reallocate zone 3 to 1 or 2 depending on shot proximity (requested by fishers at AMM 2019)
Zn3.lim.eas=116+55.4/60
Zn3.lim.wes=116+30/60
Zn3.mid=(Zn3.lim.wes+Zn3.lim.eas)/2
Data.daily=Data.daily%>%
            mutate(zone=ifelse(zone=="Zone2"& !is.na(LongMin) & (LongDeg+LongMin/60)<=Zn3.mid,"Zone1",zone))


if(nrow(Data.daily.FC.2005_06)>0)
{
  Data.daily.FC.2005_06$LatDeg=with(Data.daily.FC.2005_06,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
  Data.daily.FC.2005_06$LatMin=with(Data.daily.FC.2005_06,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
  Data.daily.FC.2005_06$LongDeg=with(Data.daily.FC.2005_06,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
  Data.daily.FC.2005_06$LongMin=with(Data.daily.FC.2005_06,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))
  
  Data.daily.FC.2005_06$zone=as.character(Data.daily.FC.2005_06$zone)
  Data.daily.FC.2005_06$zone=as.character(with(Data.daily.FC.2005_06,
                ifelse(is.na(LatDeg) | is.na(LatMin) | is.na(LongDeg) | is.na(LongMin),zone,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",
                NA)))))))))
  Data.daily.FC.2005_06$zone=with(Data.daily.FC.2005_06,
                ifelse(zone=='2',"Zone2",
                ifelse(zone=='1',"Zone1",
                ifelse(zone=='*',"West",
                zone))))
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
                ifelse(is.na(LatDeg) | is.na(LatMin) | is.na(LongDeg) | is.na(LongMin),zone,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",
                       NA)))))))))
  Data.daily.FC.NA.sp$zone=with(Data.daily.FC.NA.sp,
                ifelse(zone=='2',"Zone2",
                ifelse(zone=='1',"Zone1",
                ifelse(zone=='*',"West",
                zone))))
}
if(nrow(Data.daily.incomplete)>0)
{
  Data.daily.incomplete$LatDeg=with(Data.daily.incomplete,ifelse(is.na(LatDeg),as.numeric(substr(block10,1,2)),LatDeg))  
  Data.daily.incomplete$LatMin=with(Data.daily.incomplete,ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin))  
  Data.daily.incomplete$LongDeg=with(Data.daily.incomplete,ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg))
  Data.daily.incomplete$LongMin=with(Data.daily.incomplete,ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin))
  
  Data.daily.incomplete$zone=as.character(Data.daily.incomplete$zone)
  Data.daily.incomplete$zone=as.character(with(Data.daily.incomplete,
                ifelse(is.na(LatDeg) | is.na(LatMin) | is.na(LongDeg) | is.na(LongMin),zone,
                ifelse((LongDeg+LongMin/60)>=116.5 & (LatDeg+LatMin/60)>=(26),"Zone2",
                ifelse((LongDeg+LongMin/60)<116.5 & (LatDeg+LatMin/60)>=(33),"Zone1",
                ifelse((LatDeg+LatMin/60)<(33) & (LatDeg+LatMin/60)>=(26) & (LongDeg+LongMin/60)<116.5,"West",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)<114,"Closed",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=114 & (LongDeg+LongMin/60)<123.75,"North",
                ifelse((LatDeg+LatMin/60)<(26) & (LongDeg+LongMin/60)>=123.75,"Joint",
                NA)))))))))
  Data.daily.incomplete$zone=with(Data.daily.incomplete,
                ifelse(zone=='2',"Zone2",
                ifelse(zone=='1',"Zone1",
                ifelse(zone=='*',"West",
                zone))))
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
              "LongDeg","LatMin","LongMin","fishery"),names(Data.daily))]
Effort.daily=Effort.daily%>%rename(FisheryCode=fishery)
Effort.daily$year.c=Effort.daily$year
Effort.daily$LatDeg=abs(Effort.daily$LatDeg)
Effort.daily$LAT=Effort.daily$LatDeg+(Effort.daily$LatMin/60)
Effort.daily$LAT=-abs(Effort.daily$LAT)
Effort.daily$LONG=Effort.daily$LongDeg+(Effort.daily$LongMin/60)
Effort.daily=Effort.daily[,-match(c("LatDeg","LongDeg","LatMin","LongMin"),names(Effort.daily))]

#Export reported depth data for indicator species
write.csv(subset(Data.daily,species%in%c(17003,17001,18003,18007)),
          handl_OneDrive("Analyses/Data_outs/Daily.log.depth.csv"),row.names=F)

#Get data for Matt Koopman Fishery Review 2023 
if(get.Mats.data.2023)
{
  Mat.daily=Data.daily%>%
    rename_with(toupper)%>%
    rename(SNAME=SNAME1,
           LAT=LATFC,
           LONG=LONGFC)%>%
    dplyr::select(FINYEAR,YEAR,MONTH,FDAYS,METHOD,BLOCKX,BDAYS,HOURS,HOOKS,SHOTS,NETLEN,
                  SPECIES,SNAME,LIVEWT,FISHERY,PORT,TYPE,VESSEL,LAT,LONG,ZONE,DEPTHMAX)
  
  Mat.daily=Mat.daily%>%
    left_join(Vesl.id,by='VESSEL')%>%
    dplyr::select(-VESSEL)%>%
    mutate(LAT=as.numeric(LAT),
           LONG=as.numeric(LONG))%>%
    filter(LAT<=(-26))%>%
    mutate(Shark.fishery=case_when(ZONE=='Joint' & FISHERY%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                   ZONE=='North' & FISHERY%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                   ZONE=='Closed' & FISHERY%in%c('OT')~'WANCS',
                                   ZONE=='West' & FISHERY%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                   is.na(ZONE) & FISHERY=='WCGL'~'WCDGDL',
                                   ZONE=='Zone1' & FISHERY=='WCGL'~'JASDGDL',
                                   FISHERY%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                   ZONE%in%c('Zone1','Zone2') & FISHERY=='*'~'OASC',
                                   ZONE%in%c('Zone1','Zone2') & FISHERY=="WL"~'OASC',
                                   ZONE%in%c('West') & FISHERY=='*'~'OANCGCWC',
                                   ZONE%in%c('West','Zone1') & FISHERY=='OT'~'OANCGCWC',
                                   METHOD=='HL' & is.na(FISHERY) ~ 'Handline shark fishery',
                                   METHOD=='HL' & FISHERY=='SBS' ~ 'Shark Bay Pink Snapper fishery',
                                   METHOD=='DL' & is.na(FISHERY) ~ 'Dropline shark fishery', 
                                   METHOD=='LL' & FISHERY=='CSLP' ~ 'estuary fishery',
                                   METHOD=='GN' & FISHERY%in%c('C019','C066','CSFN','MBC','SCE','WCE') ~ 'estuary fishery',
                                   TRUE~"other fishery"))%>%
    filter(Shark.fishery%in%c('Dropline shark fishery','Handline shark fishery','JASDGDL','OANCGCWC','OASC','WCDGDL'))%>%
    rename_with(toupper)
  
  
}


# B.1. Catch Inspections             

#note: 1. run this initially, identify errors, get data entry girls to check, 
#      2.  raise issues to Eva for correcting current catch return before analyses 
#      3.  then fix in "Daily data fixes"

  #Select data for current year assessment
#note: Current.data is a dummy, only used to inspect current catch and effort.
Yr.current=2000+as.numeric(substr(Current.yr,6,7))
Current.data=subset(Data.daily.1,finyear==Current.yr)

#note: for SoFAR, Data.current.Sofar is only used to get the number of crew and licences
Data.current.Sofar=Current.data
Data.current.Sofar$LAT=-abs(Data.current.Sofar$LatDeg)  
Data.current.Sofar=Data.current.Sofar[,-match(c("LatMin","LongMin"),names(Data.current.Sofar))]  

Current.data=Current.data%>%
  rename(SessionStartDate=date,
         VesselName=BoatName,
         TripSheetNumber=TSNo,
         DailySheetNumber=DSNo,
         SessionIndex=SNo,
         RSSpeciesCommonName=RSCommonName,
         Nfish=nfish,
         Landwt=landwt,
         TripLandedCondition=conditn)

if(Inspect.New.dat=="YES")  
{
  #create file and path for checking new data
  handle=paste(handl_OneDrive("Data/Catch and Effort/Check_these_vars/Daily/"),Current.yr,sep="")
  if(!file.exists(handle)) dir.create(handle)
  
  
  #vars needed by Vero to inspect data
  Vero.vars=c('SessionStartDate','VesselName','TripSheetNumber','DailySheetNumber',
              'SessionIndex','RSSpeciesCommonName','Nfish','Landwt','livewt','RSSpeciesId','RSSpeciesCode')

  Top.mon.ktch=mean(c(15000,20000))
  Top.mon.eff=263.5 # 8500 m X 31 days
  # (i.e. monthly catch >Top.mon.ktch and monthly Km.gn.d > Top.mon.eff)  VIP!!!
  
  
  #Check spatial blocks and catch assigned to right zone
  fn.ck.spatial(d=Current.data%>%
                  filter(fishery%in%c('JASDGDL','WCDGDL'))%>%
                  rename(date=SessionStartDate,
                         BoatName=VesselName,
                         TSNo=TripSheetNumber,
                         DSNo=DailySheetNumber,
                         SNo=SessionIndex)%>%
                  mutate(yr=finyear,
                         zn=zone,
                         ln=Long,
                         la=Lat,
                         blk=block10,
                         Unik=paste(SNo,DSNo,TSNo))%>%
                  distinct(Unik,.keep_all = T),
                BLK.data=BLOCK_60%>%
                  mutate(dumi=as.numeric(BlockCode))%>%
                  filter(dumi<85010)%>%
                  rowwise()%>% 
                  mutate(x = mean(c(NorthWestPointGPSLongitude, SouthEastPointGPSLongitude)),
                         y = mean(c(NorthWestPointGPSLatitude, SouthEastPointGPSLatitude))),
                Depth.data=Bathymetry,
                BLK.size=4.5,
                add.yr=TRUE,
                pt.size=1.5,
                pt.alpha=0.35)
  ggsave(paste(handle,"Map_blocks_daily.tiff",sep='/'),width = 14,height = 10,compression = "lzw")
  
  ODD.blocks=c(33260,35190) #update after visual inspection of map 
  out.catch.block=Current.data%>%
    filter(fishery%in%c('JASDGDL','WCDGDL'))%>%
    rename(date=SessionStartDate,
           BoatName=VesselName,
           TSNo=TripSheetNumber,
           DSNo=DailySheetNumber,
           SNo=SessionIndex)%>%
    mutate(Unik=paste(SNo,DSNo,TSNo))%>%
    distinct(Unik,.keep_all = T)%>%
    filter(blockx%in%ODD.blocks)%>%
    filter(!(blockx==33260 & Lat>(-33.5)))%>%
    dplyr::select(fishery,BoatName,vessel,date,SNo,DSNo,TSNo,blockx,Lat,Long)
  if(nrow(out.catch.block)>0)
  {
    write.csv(out.catch.block,file=paste(handle,"/Check.catch outside normal fishing grounds.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'NA_nfish.weight.csv',sep=' _ '),
               Body= "Hi,
              I have started the data validation process for the new year data stream.
              I will be sending a series of emails with queries.
              In the attached, I extracted records where fishing occurs outside normal fishing grounds
              Cheers
              Matias",  
               Attachment=paste(handle,"/Check.catch outside normal fishing grounds.csv",sep="")) 
  }
  
    #NA in zone
  out.catch.no.zone=Current.data%>%
    filter(fishery%in%c('JASDGDL','WCDGDL'))%>%
    rename(date=SessionStartDate,
           BoatName=VesselName,
           TSNo=TripSheetNumber,
           DSNo=DailySheetNumber,
           SNo=SessionIndex)%>%
    mutate(Unik=paste(SNo,DSNo,TSNo))%>%
    distinct(Unik,.keep_all = T)%>%
    filter(is.na(zone))%>%
    dplyr::select(fishery,BoatName,vessel,date,SNo,DSNo,TSNo,zone)
  if(nrow(out.catch.no.zone)>0)
  {
    write.csv(out.catch.no.zone,file=paste(handle,"/Check.NA zone.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'NA_nfish.weight.csv',sep=' _ '),
               Body= "Hi,
              I have started the data validation process for the new year data stream.
              I will be sending a series of emails with queries
              In the attached, I extracted records with NA in zone
              Cheers
              Matias",  
               Attachment=paste(handle,"/Check.NA zone.csv",sep="")) 
  }
  
    #wrong zone
  out.check.zone=Current.data%>%
    filter(fishery%in%c('JASDGDL','WCDGDL'))%>%
    rename(date=SessionStartDate,
           BoatName=VesselName,
           TSNo=TripSheetNumber,
           DSNo=DailySheetNumber,
           SNo=SessionIndex)%>%
    mutate(Unik=paste(SNo,DSNo,TSNo))%>%
    distinct(Unik,.keep_all = T)%>%
    mutate(zone.true=zone,
           zone.true=case_when(Lat<(-33) & Long>116.5 & Long<(116+55.40/60)~'3',
                               Lat<(-30) & Long>=(116+55.40/60)~'2',
                               Lat<=(-33) & Long<=116.5~'1',
                               TRUE~zone.true))%>%
    filter(!zone.true==zone)%>%
    dplyr::select(fishery,vessel,date,SNo,DSNo,TSNo,blockx,Lat,Long,zone,zone.true)
  if(nrow(out.check.zone)>0)
  {
    write.csv(out.check.zone,file=paste(handle,"/Check.wrong zone.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'NA_nfish.weight.csv',sep=' _ '),
               Body= "Hi,
              I have started the data validation process for the new year data stream.
              I will be sending a series of emails with queries
              In the attached, I extracted records for which the zone is incorrect
              Cheers
              Matias",  
               Attachment=paste(handle,"/Check.wrong zone.csv",sep="")) 
  }
  
  
  #1. ID missing catch weight records in current year data when nfish reported  
  NA.livewt=Current.data%>%
              filter(is.na(livewt) & !TripLandedCondition%in%c('SC','NR'))%>%
              distinct(SessionStartDate,vessel,TripSheetNumber,DailySheetNumber,SessionIndex,
                       RSSpeciesId,RSSpeciesCommonName,RSSpeciesCode,Nfish,Landwt,livewt,
                       TripLandedWeight,TripLandedCondition,factor)%>%
              arrange(vessel)%>%
              rename(Sp.Conv.Factor=factor)
  if(nrow(NA.livewt)>0)
  {
    write.csv(NA.livewt,file=paste(handle,"/Check.NALivewt_missing  conversion factor for condition.species.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'NA_nfish.weight.csv',sep=' _ '),
               Body= "Hi,
              I have started the data validation process for the new year data stream.
              I will be sending a series of emails with queries.
              In the attached, I extracted records with species and ‘nfish’ information but 
              NA in ‘livewt’ and ‘landwt’. Is this a typo or system error?
              Cheers
              Matias",  
               Attachment=paste(handle,"/Check.NALivewt_missing  conversion factor for condition.species.csv",sep="")) 
  }
  #NAs=subset(Current.data,is.na(livewt))
  #NA.table=table(as.character(NAs$sname1),useNA='ifany')
  
  
  #2. Identify zero catch shots (i.e. NA species and NA weight)
  table(Current.data$NilCatch)  
  NA.shots=subset(Current.data,is.na(as.character(sname1)))
  NA.Sheet.shot=unique(paste(as.character(NA.shots$DailySheetNumber),as.character(NA.shots$SessionIndex)))
  Current.data$Sheet.shot=with(Current.data,paste(as.character(DailySheetNumber),as.character(SessionIndex)))
  Cero.catch.shot=subset(Current.data,Sheet.shot%in%NA.Sheet.shot)
  table(Cero.catch.shot$species,useNA='ifany')
  
  Cero.catch.shot=Current.data%>%   #not needed, this is taken care below 'if(nrow(check.nfish.weight)>0)'
    filter((Landwt==0 | livewt==0) & !RSSpeciesCommonName=="Nil Fish Caught")%>%
    filter(!TripLandedCondition%in%c('SC','NR') )  #remove self consumption

  #3. Check good weight calculation
  #note:this compares average weights from returns to the range of possible weights by species.
  #     for some species there's no Wei.range info.
  
    #Check if there's numbers data and landwt but no livewt  (this is amended in 'Fix weights with NA weight or zero weight but with nfish')  
  check.nfish.weight=Current.data%>%
                    filter(!is.na(Nfish) | !Nfish==0)%>%
                    filter(!RSSpeciesCommonName=="Nil Fish Caught")%>%
                    filter(!TripLandedCondition%in%c('SC','NR'))%>%
                    filter(is.na(livewt) | livewt==0)%>%
                    dplyr::select(vessel,SessionStartDate,DailySheetNumber,TripSheetNumber,
                                  species,RSSpeciesCode,RSSpeciesCommonName,Nfish,Landwt,
                                  livewt,TripLandedCondition)%>%
                    arrange(vessel,SessionStartDate)
  if(nrow(check.nfish.weight)>0)
  {
    par(bg=2)
    plot.new()
    mtext("there are records ",3,cex=3,col="white")
    mtext("with nfish but no weight",3,-2,cex=3,col="white")
    par(bg="white")
    #stop("records with nfish but no weight")
    if(nrow(NA.livewt)>0)
    {
      check.nfish.weight=check.nfish.weight%>%
        mutate(dummy=paste(DailySheetNumber,TripSheetNumber,RSSpeciesCommonName))%>%
        filter(!dummy%in%paste(NA.livewt$DailySheetNumber,NA.livewt$TripSheetNumber,NA.livewt$RSSpeciesCommonName))%>%
        dplyr::select(-dummy)
    }
    if(nrow(check.nfish.weight)>0)
    {
      write.csv(check.nfish.weight%>%mutate(livewt=ifelse(is.na(livewt),'',livewt)),
                file=paste(handle,"/Check.no_live.weight.csv",sep=""),row.names=F)
      send.email(TO=Email.data.checks,
                 CC=Email.data.checks2,
                 BCC=Email.FishCube,
                 Subject=paste("Shark validation",Sys.time(),'no_live.weight.csv',sep=' _ '),
                 Body= "Hi,
              In the attached, I extracted records where ‘landwt’ and/or ‘livewt’ is NA or 0 but there are 
              values for ‘nfish’.
              Is this a typo or there’s a legit reason for not having a ‘livewt’ value?
              Cheers
              Matias",  
                 Attachment=paste(handle,"/Check.no_live.weight.csv",sep=""))
      
    }
  }

  #4.1 Check average weight for each species
  Current.data$Avg.wt=with(Current.data,livewt/Nfish)
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
  tolerance=0.5   #tolerance bounds for acceptable weight
  fn.avg.wt=function(sp,sp.name)
  {
    Data=subset(Current.data,species==sp & !(is.na(Avg.wt)))
    Weit.r=subset(Wei.range, SPECIES==sp)
    
    if(nrow(Data)>1)
    {
      R=range(Data$Avg.wt)
      R.n=range(Data$Nfish)
      R.w=range(Data$livewt)
      
      par(mfcol=c(2,2),mai=c(1.1,.85,.2,.1),oma=c(.1,.1,1,.1))
      hist(Data$Nfish,xlab="Numbers per shot",main=sp.name,cex.main=1.2,col=2)
      hist(Data$livewt,xlab="Weight (kg) per shot",main="",col=2)
      hist(Data$Avg.wt,xlab="Average weight (kg)",main=paste0("Biological w8t range: ",round(Weit.r$TW.min,1),
                                                             "-",round(Weit.r$TW.max,1)," kg"),col=2)
      
      Q=quantile(Data$Avg.wt,probs = seq(0, 1, 0.025))
      Q.n=quantile(Data$Nfish,probs = seq(0, 1, 0.025))
      
      Rev.weight=0
      if(sum(Data$Avg.wt<(Weit.r$TW.min*(1-tolerance)) | (Data$Avg.wt>Weit.r$TW.max*(1+tolerance)))>0)
      {
        id=which(Data$Avg.wt<(Weit.r$TW.min*(1-tolerance)) | Data$Avg.wt>(Weit.r$TW.max*(1+tolerance)))
        Rev.weight=Data[id,]
        plot(Rev.weight$Nfish,Rev.weight$livewt,xlab="nfish",ylab="livewt",
             main='Records outside biologcial w8t range',pch=19,col=2,cex.main=.8,xlim=c(0,max(Rev.weight$Nfish)),ylim=c(0,max(Rev.weight$livewt)))
      }else
      {
        plot.new()
        mtext("no records outside w8t range",3,cex=1,col="forestgreen")
      }
      
      return(list(W.range=R.w,Avg.w.range=R,num.range=R.n,Avg.w.quantiles=Q,num.quantiles=Q.n,Review.Record=Rev.weight))
    }
    if(nrow(Data)<=1)return("PROBLEM")
  }
  Avg.wt.list=vector('list',length=length(Uniq.sp.with.weight))
  names(Avg.wt.list)=Uniq.sp.nam.with.weight
  for (i in 1:length(Uniq.sp.with.weight)) Avg.wt.list[[i]]=fn.avg.wt(sp=Uniq.sp.with.weight[i],
                                                                      sp.name=Uniq.sp.nam.with.weight[i])
  Current.data=Current.data%>%left_join(Wei.range%>%dplyr::select(TW.min,TW.max,SPECIES),
                                        by=c("species"="SPECIES"))
  Current.data$Chk.wt=with(Current.data,
            ifelse(Avg.wt<(TW.min*(1-tolerance)) | 
            Avg.wt>(TW.max*(1+tolerance)),"check","ok"))
  check.weights=subset(Current.data,Chk.wt=="check",
     select=c(SessionStartDate,vessel,TripSheetNumber,DailySheetNumber,SessionIndex,vessel,
              RSSpeciesId,RSSpeciesCommonName,flagtype,factor,TripLandedCondition,Nfish,Landwt,
              livewt,Avg.wt,TW.min,TW.max))%>%
            rename(Sp.Conv.Factor=factor)
  idd=match(c("Landwt","livewt","Avg.wt","TW.min","TW.max"),names(check.weights))
  if(nrow(check.weights)>0)
  {
    check.weights[,idd]=round(check.weights[,idd],2)
    check.weights=check.weights%>%
      arrange(vessel,TripSheetNumber,DailySheetNumber,RSSpeciesCommonName)%>%
      rename(Derived_Avg.wt=Avg.wt,
             Min.wt=TW.min,
             Max.wt=TW.max)
    write.csv(check.weights,file=paste(handle,"/Check.nfish.weight.typo.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'nfish.weight.typo.csv',sep=' _ '),
               Body= "For the attached file, could you please check that there are no typos in the nfish, landed weight and condition
                      columns as the average weight resulting from these 3 columns are outside the range for the species (TW.min & TW.max). 
                    Cheers
                    Matias",
               Attachment=paste(handle,"/Check.nfish.weight.typo.csv",sep=""))
  }

  
  #4.2 Compare average weight indicator species   
  #note: this is not required, already taken care above
  do.indic="NO"
  if(do.indic=="YES")
  {
    Avg.wt.list.ind=vector('list',length=length(Indicator.species))
    names(Avg.wt.list.ind)=names(Indicator.species)
    for (i in 1:length(Indicator.species)) Avg.wt.list.ind[[i]]=fn.avg.wt(Indicator.species[i],names(Indicator.species)[i])  
    
    
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
  
  #5.flag species outside spatial distribution 
  library(sp)
  Current.data$Lat=-abs(as.numeric(as.character(Current.data$Lat)))
  Current.data$Long=as.numeric(as.character(Current.data$Long))
  Dist.range=list('17001'=list(c(113,-29),c(129,-36)),  #Gummy Shark
                  '17003'=list(c(113,-21),c(129,-36)),  #Whiskery Shark
                  '18001'=list(c(113,-29),c(129,-36)),  #Bronze Whaler
                  '18003'=list(c(113,-15),c(129,-36)),  #Dusky Whaler
                  '18007'=list(c(113,-15),c(125,-36)),  #Sandbar Shark
                  '23002'=list(c(113,-30),c(129,-36)),  #"Common Sawshark"
                  '17006'=list(c(113,-20),c(129,-36)),  #"Pencil Shark"   
                  '17008'=list(c(113,-30),c(129,-36)),  #"School Shark"
                  '5001'=list(c(115,-30),c(129,-36)),   #"Sevengill Sharks"
                  '10001'=list(c(113,-10),c(129,-36)),  #"Shortfin Mako"
                  '39001'=list(c(113,-30),c(129,-36)),  #"Southern Eagle Ray" 
                  '18023'=list(c(113,-16),c(126,-36)),  #"Spinner Shark"
                  '12000'=list(c(113,-10),c(129,-36)),  #"Thresher Shark"
                  '18022' =list(c(113,-14),c(129,-36)),   #"Tiger Shark"
                  '24900' =list(c(113,-29),c(129,-36)),   #"Angel Shark"
                  '13000' =list(c(113,-26),c(129,-36)),   #"Wobbegong"
                  '19001' =list(c(113,-10),c(118,-35)),   #"Scalloped Hammerhead"
                  '19002' =list(c(113,-10),c(115,-30)),   #"Great Hammerhead"
                  '19004' =list(c(113,-20),c(129,-36))   #"Smooth Hammerhead"
              )
  dis.sp=sort(unique(Data.daily%>%
                filter(finyear==Current.yr &species%in%c(Shark.species,Ray.species))%>%
                  filter(!species%in%c(31000,20000,22999,9998,90030,27000))%>%
                  pull(species)))
  not.considered=dis.sp[which(!dis.sp%in%as.numeric(names(Dist.range)))]
  if(length(not.considered)>0)
  {
    par(bg=2)
    plot.new()
    mtext("add distribution for these species",3,cex=1,col="white")
    mtext(paste(not.considered,collapse = '-'),3,-2,cex=1,col="white")
    par(bg="white")
  }
  
  Outer=vector('list',length(Dist.range))
  smart.par(n.plots=length(Dist.range),MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(.1, 0.5, 0))
  pdf(paste(handle,"/Check.lat.and.long.pdf",sep=""))
  for (i in 1:length(Dist.range))
  {
    s=subset(Current.data,species==as.numeric(names(Dist.range)[i]))
    if(nrow(s)>0)
    {
      NM=unique(s$RSSpeciesCommonName)
      plot(s$Long,s$Lat,main=NM,ylim=c(-36,-10),xlim=c(113,129),ylab="",xlab="")
      pol=Dist.range[[i]]
      polygon(x=c(pol[[1]][1],pol[[2]][1],pol[[2]][1],pol[[1]][1]),
              y=c(pol[[1]][2],pol[[1]][2],pol[[2]][2],pol[[2]][2]),border=2)
      
      
      Outside=point.in.polygon(point.x=s$Long,
                               point.y=s$Lat,
                               pol.x=c(pol[[1]][1],pol[[2]][1],pol[[2]][1],pol[[1]][1]),
                               pol.y=c(pol[[1]][2],pol[[1]][2],pol[[2]][2],pol[[2]][2]))
      Outside=s[which(!Outside==1),]
      if(nrow(Outside)>0) Outer[[i]]=Outside%>%dplyr::select(c(Vero.vars[1:length(Vero.vars)],vessel,RSSpeciesCode,Lat,Long))
      
    }
  }
  dev.off()
  Outer=do.call(rbind,Outer)
  Outer=Outer%>%filter(!grepl('Hammerhead',RSSpeciesCommonName)) #don't consider hammerheads due to reapportioning
  if(nrow(Outer)>0)
  {
    write.csv(Outer,file=paste(handle,"/Check.lat.and.long.csv",sep=""),row.names=F)
    send.email(TO=Email.data.checks,
               CC=Email.data.checks2,
               BCC=Email.FishCube,
               Subject=paste("Shark validation",Sys.time(),'lat.and.long.typo.csv',sep=' _ '),
               Body= "For the attached file, could you please check the latitude and longitude of the shot?
                      Currently they fall outside the species distribution. 
                    Cheers
                    Matias",
               Attachment=paste(handle,"/Check.lat.and.long.csv",sep=""))
  }
  
  Current.data%>%
    mutate(dummy=paste(Long,Lat),
           fishery.zone=paste(fishery,zone))%>%
    distinct(dummy,.keep_all = T)%>%
    ggplot(aes(Long,Lat,color=fishery.zone))+
    geom_point()
  ggsave(paste(handle,"/Lat.and.long.of.records.tiff",sep=""),width = 6,height = 6,compression = "lzw")
  
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
THIS.10=c("finyear","year","month","day","fishery","zone","orgzone","targetSpecies","block10","blockx",
          "LatDeg","LatMin","LongDeg","LongMin","method","species","sname1","vessel","NilCatch",
          "hours","netlen","shots","depthMax","depthMin","landwt","livewt","bdays")     
Fine.scale.eff=Data.daily.1
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
Data.daily$LatDeg=as.numeric(substr(abs(Data.daily$LatDeg),1,2))
Data.daily$LatDeg=as.numeric(substr(Data.daily$LatDeg,1,3))+(Data.daily$LatMin/60)
Data.daily$LongDeg=Data.daily$LongDeg+(Data.daily$LongMin/60)
Data.daily$LatDeg=-abs(Data.daily$LatDeg)

#Fix some bioregions
Data.daily$Bioregion.old=Data.daily$Bioregion
# Data.daily$Bioregion=as.character(with(Data.daily,
#   ifelse(LongDeg>=115.5 & LongDeg<=129 & LatDeg<=(-26),"SC", 
#   ifelse(LongDeg<115.5 & LatDeg<=(-27),"WC",
#   ifelse(LongDeg<=114.834 & LatDeg>(-27),"Gascoyne",
#   ifelse(LongDeg>114.834 & LatDeg>(-27) & LongDeg<=129,"NC",NA))))))
# Data.daily$Bioregion=with(Data.daily,
#   ifelse(Bioregion=="SC"& LatDeg>(-34) & LongDeg<115.91 ,"WC",Bioregion))


#Export weight and nfish data from TDGDLF for mean weight analysis
#note: for this analysis, remove records outside weight range
if(First.run=='YES')
{
  Mn.wght.dat=subset(Data.daily,method=="GN" & !(blockx%in%Estuaries) & netlen >100)
  Mn.wght.dat$Avrg.w=Mn.wght.dat$livewt/Mn.wght.dat$nfish
  Mn.wght.dat$Keep=with(Mn.wght.dat,
                        ifelse( Avrg.w>max.w.whis & Avrg.w<min.w.whis & species==17003,"NO",
                                ifelse( Avrg.w>max.w.gum & Avrg.w<min.w.gum & species==17001,"NO",
                                        ifelse( Avrg.w>max.w.dus & Avrg.w<min.w.dus & species==18003,"NO",
                                                ifelse( Avrg.w>max.w.san & Avrg.w<min.w.san & species==18007,"NO",
                                                        ifelse(is.na(livewt) | livewt<=0,"NO","YES"))))))                     
  Mn.wght.dat=subset(Mn.wght.dat,Keep=="YES")               
  Mn.wght.dat=Mn.wght.dat[,-match(c("Avrg.w","Keep"),names(Mn.wght.dat))]
  write.csv(Mn.wght.dat%>%
              filter(!is.na(livewt))%>%
              left_join(Mesh.size%>%
                          mutate(mshigh=round(mshigh))%>%
                          dplyr::distinct(DSNo,TSNo,SNo,mshigh),
                        by = c("SNo", "DSNo", "TSNo")),
            file =handl_OneDrive("Analyses/Catch and effort/Logbook.data.mean.weight.csv"))
  rm(Mn.wght.dat)
}


#Amend weights with NA weight or zero weight but with nfish
#note: This is needed to account for population extractions but not used for cpue stand. 
#      Criteria:  Use species average for the block; if not available, use the overall average for the species
Avrg.w8ts_daily_spcies.blk10=Data.daily%>%
              left_join(Wei.range%>%dplyr::select(SPECIES,TW.min,TW.max),by=c('species'='SPECIES'))%>%
              group_by(species,block10) %>% 
              mutate(n=n()) %>%
              mutate(avrg_livewt=ifelse(n>10 & !is.na(livewt) & !is.na(nfish),livewt/nfish,NA),
                     within.range=ifelse(!is.na(TW.min) & !is.na(avrg_livewt) &
                                           (avrg_livewt>TW.max | avrg_livewt<TW.min),'NO','YES'))%>%
              filter(avrg_livewt>0.1 & avrg_livewt<100)%>%
              filter(within.range=='YES')%>%
              group_by(species,block10) %>% 
              summarise(avrg_livewt_sp_blk10=mean(avrg_livewt, na.rm=TRUE))%>%
              ungroup()%>%
              data.frame()

Avrg.w8ts_daily_spcies=Data.daily%>%
              left_join(Wei.range%>%dplyr::select(SPECIES,TW.min,TW.max),by=c('species'='SPECIES'))%>%
              group_by(species) %>% 
              mutate(n=n()) %>%
              mutate(avrg_livewt=ifelse(n>10 & !is.na(livewt) & !is.na(nfish),livewt/nfish,NA),
                     within.range=ifelse(!is.na(TW.min) & !is.na(avrg_livewt) &
                                           (avrg_livewt>TW.max | avrg_livewt<TW.min),'NO','YES'))%>%
              filter(avrg_livewt>0.1 & avrg_livewt<100)%>%
              filter(within.range=='YES')%>%
              group_by(species) %>% 
              summarise(avrg_livewt_sp=mean(avrg_livewt, na.rm=TRUE))%>%
              ungroup()%>%
              data.frame()

Data.daily=Data.daily %>% 
            left_join(Avrg.w8ts_daily_spcies.blk10,by=c('species','block10'))%>% 
            left_join(Avrg.w8ts_daily_spcies,by=c('species'))%>%
            mutate(Reporter= ifelse((is.na(livewt)|livewt==0) & nfish>0 & !conditn=='FI',"bad","good"),    
                   livewt= ifelse((is.na(livewt)|livewt==0) & nfish>0, 
                                  nfish*avrg_livewt_sp_blk10, livewt)) %>%
            mutate(livewt= ifelse((is.na(livewt)|livewt==0) & nfish>0, 
                                  nfish*avrg_livewt_sp, livewt)) %>%
            mutate(livewt=ifelse(is.na(livewt) & !is.na(landwt) & !is.na(nfish),
                                 landwt,livewt))%>%
            dplyr::select(-c(avrg_livewt_sp_blk10,avrg_livewt_sp))

#self consumption and no weight
Average.self.con=Data.daily%>%
                  filter(conditn%in%c("OT","SC"))%>%
                  group_by(species)%>%
                  summarise(Average.weight.self.consumed=sum(livewt,na.rm=T)/sum(nfish,na.rm=T))
Data.daily=Data.daily %>% 
            left_join(Average.self.con,by=c('species')) %>%
  mutate(Reporter= ifelse((is.na(livewt)|livewt==0) & is.na(nfish) & conditn%in%c("OT","SC"),"bad",Reporter),    
         livewt= ifelse((is.na(livewt)|livewt==0) & is.na(nfish) & conditn%in%c("OT","SC"), 
                        Average.weight.self.consumed, livewt)) %>%
  dplyr::select(-Average.weight.self.consumed)
  

#Aggregate daily nfish 
Daily.nfish.agg=aggregate(nfish~finyear+year+month+vessel+method+blockx+blockxFC+species+
                    year.c+Same.return+zone+Estuary,data=Data.daily,sum,na.rm=T)
names(Daily.nfish.agg)=c("FINYEAR","YEAR","MONTH","VESSEL","METHOD",
                         "BLOCKX","blockxFC","SPECIES","YEAR.c","Same.return","zone","Estuary","nfish")


# B.5. sort columns and drop some variables
Data.daily$LAT=-as.numeric(substr(abs(Data.daily$LatDeg),1,2))  
Data.daily$LONG=as.numeric(substr(Data.daily$LongDeg,1,3))

Data.daily.depth.max=subset(Data.daily,select=c(Same.return.SNo,depthMax))
Data.daily=Data.daily[,-match(c("depthMax"),names(Data.daily))]

#add variables for FishCUBE
Data.daily$FisheryZone=Data.daily$zone
Data.daily$FisheryCode=Data.daily$fishery
Data.daily$Landing.Port=Data.daily$port

this= c("finyear","year","month","vessel","fdays","method","blockx","bdays","hours","hooks","shots",
        "netlen","species","sname1","livewt","conditn","landwt","factor","year.c","LAT","LONG",
        "Estuary","Same.return","Same.return.SNo","TYPE.DATA","Bioregion","zone","LatMin","LongMin",
        "day","block10","FisheryZone","fishery","FisheryCode","Landing.Port",FishCube.vars,"nfish",
        "type","RSSpeciesCode","Reporter")
Data.daily=Data.daily[,match(this,names(Data.daily))]
Data.daily=Data.daily[,sort(names(Data.daily))]
names(Data.daily)=sort(c(names(Data.monthly),"day","block10"))   


# B.6. Aggregate Data.daily by Same.Return (i.e. finyear-month-vessel-method-block)    
#note: Monthly records have ~ 3 times less species reported than Daily logbooks so monthly-aggregated daily logbooks
#       have more rows than Monthly records
Data.daily.agg=aggregate(cbind(LANDWT,LIVEWT)~FINYEAR+YEAR+MONTH+VESSEL+METHOD+BLOCKX+blockxFC+
                           SPECIES+SNAME+CONDITN+Factor+YEAR.c+LAT+LONG+Same.return+
                           TYPE.DATA+Bioregion+zone+Estuary+RSCommonName+RSSpeciesId+
                           FisheryZone+FisheryCode+PORT+Reporter,
                         data=Data.daily[,-match(c("day","block10"),names(Data.daily))]%>%
                           mutate(LANDWT=ifelse(is.na(LANDWT),1e9,LANDWT),
                                  CONDITN=ifelse(is.na(CONDITN),'unknown',CONDITN),
                                  Factor=ifelse(is.na(Factor),1e9,Factor),
                                  METHOD=case_when(is.na(METHOD) & FisheryCode=='PFT'~'TW',
                                                   is.na(METHOD) & FisheryCode=='NCS'~'LL',
                                                   is.na(METHOD) & FisheryCode=='NDS'~'FT',
                                                   TRUE~METHOD),
                                  LAT=ifelse(is.na(LAT),as.numeric(LatFC),LAT),
                                  LONG=ifelse(is.na(LONG),as.numeric(LongFC),LongFC),
                                  TYPE.DATA=ifelse(is.na(TYPE.DATA),'daily',TYPE.DATA),
                                  Bioregion=case_when(is.na(Bioregion) & LONG>=115.5 & LONG<=129 & LAT<=(-26)~"SC", 
                                                      is.na(Bioregion) & LONG<115.5 & LAT<=(-27)~"WC",
                                                      is.na(Bioregion) & LONG<=114.834 & LAT>(-27)~"Gascoyne",
                                                      is.na(Bioregion) & LONG>114.834 & LAT>=(-27) & LONG<=129~"NC",
                                                      TRUE~Bioregion),
                                  zone=ifelse(is.na(zone),'',zone),
                                  FisheryZone=ifelse(is.na(FisheryZone),'',FisheryZone)),
                         sum,na.rm=T)%>%
                          mutate(LANDWT=ifelse(LANDWT>1e8,NA,LANDWT),
                                 CONDITN=ifelse(CONDITN=='unknown',NA,CONDITN),
                                 Factor=ifelse(Factor>1e8,NA,Factor))
Data.daily.agg=Data.daily.agg%>%rename(Landing.Port=PORT)

#Check for 0 or NA livewt records
a=Data.daily%>%filter(LANDWT==0 | LIVEWT==0 | is.na(LANDWT) | is.na(LIVEWT))%>%
  filter(!SNAME=='Nil Fish Caught')%>%
  filter(is.na(LIVEWT) | LIVEWT==0)
if(nrow(a)>0)
{
  par(bg=2)
  plot.new()
  mtext(paste("there are",nrow(a),"records"),3,cex=3,col="white")
  mtext("with no weight",3,-2,cex=3,col="white")
  mtext("& SNAME is not Nil Fish Caught",3,-4,cex=2,col="white")
  mtext(paste("they are",paste(unique(a$RSCommonName),collapse=',')),3,-6,cex=2,col="white")
  par(bg="white")
  write.csv(a%>%
              mutate(SNo=substr(Same.return.SNo,1,2),
                     DSNo=substr(Same.return.SNo,3,14),
                     TSNo=substr(Same.return.SNo,16,27))%>%
              dplyr::select(VESSEL,FINYEAR,SNo,DSNo,TSNo,block10,fishery,Landing.Port,LANDWT,LIVEWT,nfish,
                            RSCommonName,RSSpeciesCode,SPECIES,CONDITN)%>%
              arrange(VESSEL,TSNo,FINYEAR),
            paste0(handl_OneDrive("Data/Catch and Effort/Check_these_vars/Daily/"),Current.yr,
                   '/catch with no nfish, landwt or livewt but species is not nil catch.csv'),row.names = F)
  
}

#check one off spatial dist (don't repeat, only do it in Inspections of new year's data) 
if(Inspect.New.dat=="YES")
{
  handle1=handl_OneDrive("Data/Catch and Effort/Check_these_vars/Daily/")
  #yrs=sort(unique(Data.daily.original$FINYEAR))
  yrs=Current.yr
  for(i in 1:length(yrs))
  {
    fn.ck.spatial(d=Data.daily.original%>%
                    filter(FINYEAR==yrs[i] & fishery%in%c('JASDGDL','WCDGDL'))%>%
                    mutate(yr=FINYEAR,
                           zn=zone,
                           ln=Long,
                           la=Lat,
                           blk=block10,
                           Unik=Same.return.SNo)%>%
                    distinct(Unik,.keep_all = T),
                  BLK.data=BLOCK_60%>%
                    mutate(dumi=as.numeric(BlockCode))%>%
                    filter(dumi<85010)%>%
                    rowwise()%>% 
                    mutate(x = mean(c(NorthWestPointGPSLongitude, SouthEastPointGPSLongitude)),
                           y = mean(c(NorthWestPointGPSLatitude, SouthEastPointGPSLatitude))),
                  Depth.data=Bathymetry,
                  BLK.size=4.5,
                  add.yr=TRUE,
                  pt.size=1.5,
                  pt.alpha=0.35)
    ggsave(paste(handle1,yrs[i],"/Map_blocks_daily.tiff",sep='/'),width = 14,height = 10,compression = "lzw")
    
  }
  
}

#SECTION C. ---- DATA MANIPULATION - CATCH MERGING AND CORRECTIONS ----

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
  
  id=match(Factor.to.charact[i],names(Data.daily.agg))
  Data.daily.agg[,id]=as.character(Data.daily.agg[,id])
  
  id=match(Factor.to.charact[i],names(Data.daily))
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
write.csv(TDGDLF.vessels,handl_OneDrive("Data/Fishing power/TDGDLF.vessels.csv"),row.names=F)

#Get skippers
SKIPPERS=subset(Data.daily.1, method=="GN"& LatDeg>=(26) ,select=MastersName)
SKIPERS=data.frame(SKIPPERS=unique(as.character(SKIPPERS$MastersName)))
write.csv(SKIPERS,handl_OneDrive("Data/Fishing power/TDGDLF.SKIPERS.csv"),row.names=F)


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

Data.daily.agg$LONG=as.numeric(Data.daily.agg$LONG)
Data.daily.agg$LAT=as.numeric(Data.daily.agg$LAT)
Data.daily.agg$BDAYS=NA

  #merge monthly and aggregated daily 
Data.monthly=rbind(Data.monthly[match(names(Data.daily.agg),names(Data.monthly))],
                   Data.daily.agg)


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
#Conditions=subset(Conditions,spgroup=="scalefish")

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
FINYEAR.daily=subset(FINYEAR.daily,!FINYEAR.daily=="2005-06")
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
   
#Export conversion factors for fishers to check
out.ratios=FALSE
if(out.ratios)
{
  Condition.codes=data.frame(CONDITN=c("HG","WH","WF","WD","SC","FL","NU","OT","GG","HD"),
                             Description=c("Headed and gutted","Whole","Without head, without tail and fins",
                                           "Without head, with tail and fins","Self consumption","Filleted",
                                           "Number of fish (no weight reported)","Other",
                                           "Gutted and gilled","Headed"))
  
  out1=Data.daily%>%mutate(dummy=paste(SPECIES,CONDITN,Factor.c))%>%
    distinct(dummy,.keep_all = T)%>%
    filter(!SPECIES==999999 & !is.na(Factor.c) & !CONDITN%in%c("SC","NU"))%>%
    left_join(Condition.codes,by='CONDITN')%>%
    arrange(SPECIES)%>%
    rename(condition=CONDITN,Species=SNAME,Conversion.factor=Factor.c)%>%
    dplyr::select(c(Species,condition,Description,Conversion.factor))
  write.csv(out1,handl_OneDrive("Presentations\\DoF\\AMM\\2019\\Conversion.ratios.csv"),row.names = F)
}



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
Gummy=17001;Dusky_whaler=c(18003,18001);Whiskery=17003;Sandbar=18007;Hammerheads=c(19000,HammerheadSpecies);
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
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Gum.Ves.ID"
Vessel.Report$Prop.Gum.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Gum.Ves.ID),0,Prop.Gum.Ves.ID))   
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.whiskery$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Whi.Ves.ID"
Vessel.Report$Prop.Whi.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Whi.Ves.ID),0,Prop.Whi.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.dusky$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Dus.Ves.ID"
Vessel.Report$Prop.Dus.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Dus.Ves.ID),0,Prop.Dus.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.other$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Other.Ves.ID"
Vessel.Report$Prop.Other.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Other.Ves.ID),0,Prop.Other.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.school$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.Sch.Ves.ID"
Vessel.Report$Prop.Sch.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.Sch.Ves.ID),0,Prop.Sch.Ves.ID))
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.dogfish$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.DogS.Ves.ID"
Vessel.Report$Prop.DogS.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.DogS.Ves.ID),0,Prop.DogS.Ves.ID))

    #daily
if(Reapportion.daily=="NO")
{
  Vessel.Report.daily=data.frame(Same.return=as.character(unique(Data.daily$Same.return)))
  Vessel.Report.daily=Vessel.Report.daily%>%left_join(Catch.prop.school.daily$Prop.VesYrMonBlock,by=c("Same.return"))
  names(Vessel.Report.daily)[match("Proportion",names(Vessel.Report.daily))]="Prop.Sch.Ves.ID"
  Vessel.Report.daily$Prop.Sch.Ves.ID=with(Vessel.Report.daily,ifelse(is.na(Prop.Sch.Ves.ID),0,Prop.Sch.Ves.ID))
  Vessel.Report.daily=Vessel.Report.daily%>%left_join(Catch.prop.dogfish.daily$Prop.VesYrMonBlock,by=c("Same.return"))
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
Data.daily=Data.daily%>%left_join(Vessel.Report.daily,by=c("Same.return"))

  # Separate vessels into 'good' and 'bad' reporters
#C.7.2 Northern split                               #Rory's rules 3c      
  #monthly
if(Monthly.log=='monthlyrawlog')
{
  Data.monthly$Reporter=with(Data.monthly,
                             ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & TYPE.DATA=="monthly" & METHOD=="GN" & SPECIES%in%c(22999,Indicator.species)) &
                                      (Prop.Other.Ves.ID==1|(Prop.Whi.Ves.ID==0 | Prop.Dus.Ves.ID==0)),"bad","good"))
  
}

  #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog')
{
  Data.daily$Reporter=with(Data.daily,
       ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & METHOD=="GN" & SPECIES%in%c(22999,Indicator.species)) &
       (Prop.Other.Ves.ID==1|(Prop.Whi.Ves.ID==0 | Prop.Dus.Ves.ID==0)),"bad","good"))  
}


#C.7.3 SW split                                     #Rory's rules 3d                
  #monthly
if(Monthly.log=='monthlyrawlog')
{
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
  
}

  #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog')
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
if(Monthly.log=='monthlyrawlog')
{
  Data.monthly$Reporter=with(Data.monthly,
                             ifelse((LAT<=(-26) & LONG>125 & 
                                       TYPE.DATA=="monthly" & METHOD=="GN"  & 
                                       SPECIES%in%c(22999,Indicator.species)) & 
                                      (Prop.Other.Ves.ID==1|
                                         (Prop.Gum.Ves.ID==0 | Prop.Sch.Ves.ID==0)),"bad",Reporter))
  
}

  #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog')
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
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Dus.Good.spl"
Good.split$Prop.Dus.Good.spl=with(Good.split,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))

Good.split=Good.split%>%left_join(Catch.prop.gummy$Prop.GoodsplitID,by=c("GoodsplitID"))
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Gum.Good.spl"
Good.split$Prop.Gum.Good.spl=with(Good.split,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))

Good.split=Good.split%>%left_join(Catch.prop.whiskery$Prop.GoodsplitID,by=c("GoodsplitID"))
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
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Dus.Zone.Good.spl"
Zone.good.split$Prop.Dus.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                        0,Prop.Dus.Zone.Good.spl))

Zone.good.split=Zone.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrZone,by=c("ZoneID"))
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Gum.Zone.Good.spl"
Zone.good.split$Prop.Gum.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                        0,Prop.Gum.Zone.Good.spl))

Zone.good.split=Zone.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrZone,by=c("ZoneID"))
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
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Dus.Mon.Good.spl"
Monthly.good.split$Prop.Dus.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                        0,Prop.Dus.Mon.Good.spl))

Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrMon,by=c("MonthlyID"))
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Gum.Mon.Good.spl"
Monthly.good.split$Prop.Gum.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                        0,Prop.Gum.Mon.Good.spl))

Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrMon,by=c("MonthlyID"))
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
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Dus.YrZn.Good.spl"
YrZn.good.split$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                        0,Prop.Dus.YrZn.Good.spl))

YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.gummy$Prop.YrZone,by=c("ZnID"))
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Gum.YrZn.Good.spl"
YrZn.good.split$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                        0,Prop.Gum.YrZn.Good.spl))

YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.whiskery$Prop.YrZone,by=c("ZnID"))
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
Total.Sk.ktch$Shark.other.livewt=with(Total.Sk.ktch,ifelse(is.na(Shark.other.livewt),0,Shark.other.livewt))

Data.monthly=Data.monthly%>%left_join(Total.Sk.ktch,by=c("Same.return"))%>%
                            left_join(Good.split,by=c("GoodsplitID"))%>%
                            left_join(Zone.good.split,by=c("ZoneID"))%>%
                            left_join(Monthly.good.split,by=c("MonthlyID"))%>%
                            left_join(YrZn.good.split,by=c("ZnID"))

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


#C.7.8 Reapportion the catch of gummy, whiskeries and duskies   
if(Monthly.log=='monthlyrawlog') 
{
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
  Bad.Reporters$SNAME=rep(c("Bronze Whaler","Gummy Shark","Whiskery Shark","shark, other"),NroW)
  
  #daily
  if(Reapportion.daily=="YES")
  {
    Bad.Reporters.daily=rbind(Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily)
    Bad.Reporters.daily=Bad.Reporters.daily[order(Bad.Reporters.daily$Same.return),]
    
    Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
    Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME
    
    Bad.Reporters.daily$SPECIES=rep(Fix.species,NroW.daily)
    Bad.Reporters.daily$SNAME=rep(c("Bronze Whaler","Gummy Shark","Whiskery Shark","shark, other"),NroW.daily)  
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
}


#C.7.9 Make good sandbar reporters and 4h update good sandbard reporters         #Rory's rules 4g      
#note: identify vessels that report sandbars within sandbar area
if(Monthly.log=='monthlyrawlog')
{
  #identify if reporting sandbar by same return
  #monthly
  names(Catch.prop.sandbar$Prop.VesYrMonBlock)[2]="Prop.sandbar"
  Data.monthly=Data.monthly%>%left_join(Catch.prop.sandbar$Prop.VesYrMonBlock,by=c("Same.return"))
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
  Sanbar.dat.bad$SNAME=rep(c("Sandbar Shark","shark, other"),NroW.san)
  
  #daily
  if(Reapportion.daily=="YES")
  {  
    NroW.san.d=nrow(Sanbar.dat.bad.daily)
    Sanbar.dat.bad.daily=rbind(Sanbar.dat.bad.daily,Sanbar.dat.bad.daily)
    Sanbar.dat.bad.daily$Spec.old=Sanbar.dat.bad.daily$SPECIES
    Sanbar.dat.bad.daily$Sname.old=Sanbar.dat.bad.daily$SNAME
    Sanbar.dat.bad.daily$SPECIES=rep(c(18007,22999),NroW.san.d)
    Sanbar.dat.bad.daily$SNAME=rep(c("Sandbar Shark","shark, other"),NroW.san.d)  
  }
  
  
  # Reapportion "sharks, other' in "bad" sandbar      #Rory's rules 4i-4o                   
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
  Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad,sum,na.rm=T)
  names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
  Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
  
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
  write.csv(round(100*Table.Reporter[1]/sum(Table.Reporter),1),
            handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/Percent.GN.reapportioned.csv"))
}



#_LONGLINE REAPPORTIONING (within TDGDLF)_

# Merge again original proportions to data
    #monthly
Catch.prop.sandbar=fun.prop(subset(Data.monthly,LAT<=Sandbar.range[1] & LONG <= Sandbar.range[2]),18007)
Vessel.Report=Vessel.Report%>%left_join(Catch.prop.sandbar$Prop.VesYrMonBlock,by=c("Same.return"))
names(Vessel.Report)[match("Proportion",names(Vessel.Report))]="Prop.San.Ves.ID"
Vessel.Report$Prop.San.Ves.ID=with(Vessel.Report,ifelse(is.na(Prop.San.Ves.ID),0,Prop.San.Ves.ID))
Drop.these=names(Vessel.Report)[-match(c("Same.return","Prop.San.Ves.ID"),names(Vessel.Report))]
Data.monthly=Data.monthly[,-match(Drop.these,names(Data.monthly))]
Data.monthly=Data.monthly%>%left_join(Vessel.Report,by=c("Same.return"))

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


#C.8.1 Northern LL split                                     #Rory's rules 6c            
  #monthly
if(Monthly.log=='monthlyrawlog')
{
  Data.monthly$Reporter=with(Data.monthly,
                             ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & YEAR.c<2007 & TYPE.DATA=="monthly" & METHOD=="LL" &
                                       SPECIES%in%c(22999,Indicator.species)) &
                                      (Prop.Other.Ves.ID==1| (Prop.Dus.Ves.ID==0 | Prop.San.Ves.ID==0)),"bad","good"))
  
}

  #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog')
{
  Data.daily$Reporter=with(Data.daily,
                           ifelse((LAT<=(-26) & LAT>(-32) & LONG<125 & TYPE.DATA=="monthly" & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) &
                                    (Prop.Other.Ves.ID==1| (Prop.Dus.Ves.ID==0| Prop.San.Ves.ID==0)),"bad","good"))
}



#C.8.2  SW LL split                                          #Rory's rules 6d           
  #monthly
if(Monthly.log=='monthlyrawlog')
{
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
}


  #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog') 
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


#C.8.3 SE split                                              #Rory's rules 6e         
    #monthly
if(Monthly.log=='monthlyrawlog')
{
  Data.monthly$Reporter=with(Data.monthly,
                             ifelse((LAT<=(-26) & LONG>=125 & YEAR.c<2007 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) & 
                                      (Prop.Other.Ves.ID==1|(Prop.Gum.Ves.ID==0 & Prop.Sch.Ves.ID==0)),"bad",Reporter))
  Table.Reporter.LL=table(Data.monthly$Reporter,useNA='ifany')
  
}

    #daily
if(Reapportion.daily=="YES" & Monthly.log=='monthlyrawlog') 
{
  Data.daily$Reporter=with(Data.daily,
                           ifelse((LAT<=(-26) & LONG>=125 & METHOD=="LL" & SPECIES%in%c(22999,Indicator.species)) & 
                                    (Prop.Other.Ves.ID==1|(Prop.Gum.Ves.ID==0 & Prop.Sch.Ves.ID==0)),"bad",Reporter))
  Table.Reporter.LL.daily=table(Data.daily$Reporter,useNA='ifany')
  
}


#C.8.3.1 Make good split table: Yr-Mn-Block                      #Apply Rory's rules 4b
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
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Dus.Good.spl"
Good.split$Prop.Dus.Good.spl=with(Good.split,ifelse(is.na(Prop.Dus.Good.spl),0,Prop.Dus.Good.spl))
Good.split=Good.split%>%left_join(Catch.prop.gummy$Prop.GoodsplitID,by=c("GoodsplitID"))
names(Good.split)[match("Proportion",names(Good.split))]="Prop.Gum.Good.spl"
Good.split$Prop.Gum.Good.spl=with(Good.split,ifelse(is.na(Prop.Gum.Good.spl),0,Prop.Gum.Good.spl))
Good.split=Good.split%>%left_join(Catch.prop.whiskery$Prop.GoodsplitID,by=c("GoodsplitID"))
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


#C.8.3.2 Make good split table: Yr-Mn-Zone                            

  #monthly
Zone.good.split=data.frame(ZoneID=as.character(unique(Data.monthly$ZoneID)))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrZone,by=c("ZoneID"))
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Dus.Zone.Good.spl"
Zone.good.split$Prop.Dus.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Dus.Zone.Good.spl),
                                                                   0,Prop.Dus.Zone.Good.spl))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrZone,by=c("ZoneID"))
names(Zone.good.split)[match("Proportion",names(Zone.good.split))]="Prop.Gum.Zone.Good.spl"
Zone.good.split$Prop.Gum.Zone.Good.spl=with(Zone.good.split,ifelse(is.na(Prop.Gum.Zone.Good.spl),
                                                                   0,Prop.Gum.Zone.Good.spl))
Zone.good.split=Zone.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrZone,by=c("ZoneID"))
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


#C.8.3.3 Make good split table: Yr-Mn                          #Rory's rules 4c  
    #monthly
Monthly.good.split=data.frame(MonthlyID=as.character(unique(Data.monthly$MonthlyID)))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.dusky$Prop.FinYrMon,by=c("MonthlyID"))
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Dus.Mon.Good.spl"
Monthly.good.split$Prop.Dus.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Dus.Mon.Good.spl),
                                                                        0,Prop.Dus.Mon.Good.spl))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.gummy$Prop.FinYrMon,by=c("MonthlyID"))
names(Monthly.good.split)[match("Proportion",names(Monthly.good.split))]="Prop.Gum.Mon.Good.spl"
Monthly.good.split$Prop.Gum.Mon.Good.spl=with(Monthly.good.split,ifelse(is.na(Prop.Gum.Mon.Good.spl),
                                                                        0,Prop.Gum.Mon.Good.spl))
Monthly.good.split=Monthly.good.split%>%left_join(Catch.prop.whiskery$Prop.FinYrMon,by=c("MonthlyID"))
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


#C.8.4.4 Make good split table: Yr-Zn                            
    #monthly
YrZn.good.split=data.frame(ZnID=as.character(unique(Data.monthly$ZnID)))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.dusky$Prop.YrZone,by=c("ZnID"))
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Dus.YrZn.Good.spl"
YrZn.good.split$Prop.Dus.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Dus.YrZn.Good.spl),
                                                                   0,Prop.Dus.YrZn.Good.spl))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.gummy$Prop.YrZone,by=c("ZnID"))
names(YrZn.good.split)[match("Proportion",names(YrZn.good.split))]="Prop.Gum.YrZn.Good.spl"
YrZn.good.split$Prop.Gum.YrZn.Good.spl=with(YrZn.good.split,ifelse(is.na(Prop.Gum.YrZn.Good.spl),
                                                                   0,Prop.Gum.YrZn.Good.spl))
YrZn.good.split=YrZn.good.split%>%left_join(Catch.prop.whiskery$Prop.YrZone,by=c("ZnID"))
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



#C.8.4.5  Update good split catches                       #Rory's rules 4e                          
    #monthly
Recalc=c(names(Good.split)[-1],names(Zone.good.split)[-1],names(Monthly.good.split)[-1],
         names(YrZn.good.split)[-1])
Data.monthly=Data.monthly[,-match(Recalc,names(Data.monthly))]

Data.monthly=Data.monthly %>% left_join(Good.split,by="GoodsplitID") %>%
                              left_join(Zone.good.split,by="ZoneID")%>%
                              left_join(Monthly.good.split,by="MonthlyID")%>%
                              left_join(YrZn.good.split,by="ZnID")

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


#C.8.5 Reapportion the catch of gummy, whiskeries and duskies
if(Monthly.log=='monthlyrawlog')
{
  #C.8.5.1 create bad reporter files for fixing catches
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
  
  
  #C.8.5.2
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
  Bad.Reporters$SNAME=rep(c("Bronze Whaler","Gummy Shark","Whiskery Shark","shark, other"),NroW)
  
  #daily
  if(Reapportion.daily=="YES")
  {
    Bad.Reporters.daily=rbind(Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily,Bad.Reporters.daily)
    Bad.Reporters.daily=Bad.Reporters.daily[order(Bad.Reporters.daily$Same.return),]
    Bad.Reporters.daily$Spec.old=Bad.Reporters.daily$SPECIES
    Bad.Reporters.daily$Sname.old=Bad.Reporters.daily$SNAME
    Bad.Reporters.daily$SPECIES=rep(Fix.species,NroW.daily)
    Bad.Reporters.daily$SNAME=rep(c("Bronze Whaler","Gummy Shark","Whiskery Shark","shark, other"),NroW.daily)
    
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
  Agg.reap.KtcH=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)
  names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
  Bad.Reporters=Bad.Reporters%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
  
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
  Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
                                 ifelse(SPECIES==22999,Total.LIVEWT.reap-Remainder,LIVEWT.reap))
  Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,ifelse(LIVEWT.reap<0,0,LIVEWT.reap))
  
  #check that reapportioned didn't add catch
  b=aggregate(LIVEWT.reap~Same.return,Bad.Reporters,sum,na.rm=T)                                                   
  A=Agg.Sp.Ktch%>%full_join(b,by=c("Same.return"))
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
}


#REAPPORTION SANDBAR (6k-6s)
if(Monthly.log=='monthlyrawlog')
{
  #C.8.6 Make good sandbar reporters and 4h update good sandbard reporters         #Rory's rules 4g    
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
  Sanbar.dat.bad$SNAME=rep(c("Sandbar Shark","shark, other"),NroW.san)
  
  #daily
  if(Reapportion.daily=="YES")
  {  
    NroW.san.d=nrow(Sanbar.dat.bad.daily)
    Sanbar.dat.bad.daily=rbind(Sanbar.dat.bad.daily,Sanbar.dat.bad.daily)
    Sanbar.dat.bad.daily$Spec.old=Sanbar.dat.bad.daily$SPECIES
    Sanbar.dat.bad.daily$Sname.old=Sanbar.dat.bad.daily$SNAME
    Sanbar.dat.bad.daily$SPECIES=rep(c(18007,22999),NroW.san.d)
    Sanbar.dat.bad.daily$SNAME=rep(c("Sandbar Shark","shark, other"),NroW.san.d)  
  }
  
  
  #C.8.6.1 Reapportion "sharks, other' in "bad" sandbar      #Rory's rules 4i-4o                     
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
  Agg.reap.KtcH=aggregate(LIVEWT.reap.san~Same.return,Sanbar.dat.bad,sum,na.rm=T)
  names(Agg.reap.KtcH)[2]="LIVEWT.reap.scaler"
  Sanbar.dat.bad=Sanbar.dat.bad%>%left_join(Agg.reap.KtcH,by=c("Same.return"))
  
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
}


# C.8.6 Reclass catches as estuary of combined mullet, whiting and herring if catch is >10%
    #monthly
whiting.herring.mullet=subset(Data.monthly,METHOD=="GN"& SPECIES%in%c(344001,330001))
whiting.herring.mullet=subset(Data.monthly,Same.return%in%whiting.herring.mullet$Same.return)
whit.her.mul=aggregate(LIVEWT.reap~Same.return,subset(whiting.herring.mullet,SPECIES%in%c(344001,330001)),sum)
All.whit.her.mul=aggregate(LIVEWT.reap~Same.return,whiting.herring.mullet,sum)
whit.her.mul=whit.her.mul%>%left_join(All.whit.her.mul,by=c("Same.return"))
whit.her.mul$Prop.mul.her.whi=with(whit.her.mul,LIVEWT.reap.x/LIVEWT.reap.y)
whit.her.mul=whit.her.mul[,c(1,4)]
Data.monthly=Data.monthly%>%left_join(whit.her.mul,by=c("Same.return"))
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



# C.8.7 REAPPORTION CATCH FROM NON-TDGDLF CATCHES
#note: this deals with species==22999 for non TDGDFL
if(Monthly.log=='monthlyrawlog'& reapportion.ktch.other.method=="YES")
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
  Data.monthly.other.gum$SPECIES=17001; Data.monthly.other.gum$SNAME="Gummy Shark"
  Data.monthly.other.whi$SPECIES=17003; Data.monthly.other.whi$SNAME="Whiskery Shark"
  Data.monthly.other.dus$SPECIES=18003; Data.monthly.other.dus$SNAME="Bronze Whaler"
  Data.monthly.other.san$SPECIES=18007; Data.monthly.other.san$SNAME="Sandbar Shark"
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
      Data.daily.other.gum$SPECIES=17001; Data.daily.other.gum$SNAME="Gummy Shark"
      Data.daily.other.whi$SPECIES=17003; Data.daily.other.whi$SNAME="Whiskery Shark"
      Data.daily.other.dus$SPECIES=18003; Data.daily.other.dus$SNAME="Bronze Whaler"
      Data.daily.other.san$SPECIES=18007; Data.daily.other.san$SNAME="Sandbar Shark"
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


# C.8.8 REAPPORTION CATCH FROM SHARK BAY
#note: this deals with species==22999 from Shark Bay
if(Monthly.log=='monthlyrawlog')
{
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
  Data.monthly.other.dus=Data.monthly.other.san=Data.monthly.other
  Data.monthly.other.dus$SPECIES=18003
  Data.monthly.other.dus$SNAME="Bronze Whaler"
  
  Data.monthly.other.san$SPECIES=18007
  Data.monthly.other.san$SNAME="Sandbar Shark"
  
  Data.monthly.other.dus$LIVEWT.reap=with(Data.monthly.other.dus,LIVEWT*Dus.prop)
  Data.monthly.other.san$LIVEWT.reap=with(Data.monthly.other.san,LIVEWT*San.prop)
  
  Data.monthly.other=rbind(Data.monthly.other.dus,Data.monthly.other.san)
  Data.monthly.other=Data.monthly.other[,-match(c("Dus.prop","San.prop"),names(Data.monthly.other))]
  Data.monthly=rbind(Data.monthly,Data.monthly.other)
  
}


#Compare amended weights, loosing or gaining catch?
d2=subset(Data.monthly,LAT<=(-26) & !is.na(LIVEWT.reap) &TYPE.DATA=="monthly")
fn.chk.ktch(d1=subset(Data.monthly.original,Same.return
                      %in%unique(d2$Same.return)),
            d2=d2,
            VAR1="LIVEWT",VAR2="LIVEWT.reap")
#if "d1 and d2 have same catch" then OK


#Remove records with reported catch set at 0 or NA
if(Monthly.log=='monthlyrawlog')
{
  Data.monthly=subset(Data.monthly,!LIVEWT.reap==0)
  Data.daily=subset(Data.daily,!LIVEWT.reap==0)
  
  Data.monthly=subset(Data.monthly,!is.na(LIVEWT.reap))
  Data.daily=subset(Data.daily,!is.na(LIVEWT.reap))
}


#by zone
Hndl=handl_OneDrive("Analyses/Catch and effort/Compare_to_previous_approach/")
if(Monthly.log=='monthlyrawlog')
{
  if(First.run=="YES")
  {
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
  }
}



  #C.9 Add 5% to catch records prior to 1990  as per instructed by Rory                                              
Data.monthly=Data.monthly%>%
              mutate(LIVEWT.c=case_when(!is.na(LIVEWT.reap) & YEAR.c<1990 & LAT<=(-26)~LIVEWT.reap*Inc.per,
                                        !is.na(LIVEWT.reap) & YEAR.c>=1990 & LAT<=(-26)~LIVEWT.reap,
                                        is.na(LIVEWT.reap) ~LIVEWT))
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
if(!do.sql.extraction)
{
  #Monthly
  fishkiuv.monthly.yrs=sort(c(keep,unique(Data.monthly.CAESS$FINYEAR)))
  FishCUBE.monthly=subset(Data.monthly,FINYEAR%in%fishkiuv.monthly.yrs,
                          select=c(Same.return,FisheryCode,METHOD,Landing.Port,LAT,LONG,
                                   FINYEAR,YEAR.c,MONTH,FisheryZone,Bioregion,
                                   blockxFC,LIVEWT.c,VESSEL,SPECIES,SNAME,RSCommonName,
                                   RSSpeciesId,BDAYS) )  
  FishCUBE.monthly=subset(FishCUBE.monthly,!FINYEAR%in%unique(Data.daily$FINYEAR))          #remove post 2005 records already in daily
  
  Monthly.not.in.daily$LIVEWT.c=Monthly.not.in.daily$LIVEWT    
  FishCUBE.monthly=rbind(FishCUBE.monthly,Monthly.not.in.daily[,match(names(FishCUBE.monthly),
                                                                      names(Monthly.not.in.daily))])     #Add post 2005-06 records not in Daily
  FishCUBE.monthly$TYPE.DATA="monthly"
  
  
  #Daily
  FishCUBE.daily=Data.daily
  FishCUBE.daily=FishCUBE.daily[,-match("BLOCKX",names(FishCUBE.daily))]
}

#Data for Agustin Parks 2025
do.this=FALSE
if(do.this)
{
  Agus.monthly=Data.monthly%>%filter(SPECIES<=90030)%>%filter(!SPECIES%in%c(31000,22999,90030))%>%
    mutate(SNAME=ifelse(SPECIES==18003,'Dusky Whaler',SNAME))%>%
    group_by(SPECIES,SNAME,LAT,LONG,YEAR,MONTH,METHOD)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    ungroup()%>%
    rename(Latitude=LAT,
           Longitude=LONG)
  write.csv(Agus.monthly,handl_OneDrive('Parks Australia/2025_project/Data/Data sets/WA/Commercial_monthly.csv'),row.names = F)
  Agus.daily=Data.daily%>%filter(SPECIES<=90030)%>%filter(!SPECIES%in%c(31000,22999,90030))%>%
    group_by(SPECIES,SNAME,LatFC,LongFC,YEAR,MONTH,METHOD)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    ungroup()%>%
    rename(Latitude=LatFC,
           Longitude=LongFC)
  write.csv(Agus.daily,handl_OneDrive('Parks Australia/2025_project/Data/Data sets/WA/Commercial_daily.csv'),row.names = F)
  
}

  #C.10 Separate fisheries  into North and South                                                                   
#note: no catch reapportioning of effort corrections done on northern data file
  #monthly
n.row.Data.monthly=nrow(Data.monthly)
Data.monthly.north=subset(Data.monthly,LAT>(-26))   
Data.monthly=subset(Data.monthly,LAT<=(-26))
if(!n.row.Data.monthly==(nrow(Data.monthly.north)+nrow(Data.monthly)))
{
  par(bg=2)
  plot.new()
  mtext('Split of Data.monthly into North ',3,cex=2,col="white")
  mtext(' and South did not work.',3,-2,cex=2,col="white")
  mtext("See if NA's in Latitude",3,-4,cex=2,col="white")
  par(bg="white")
}
  
  #daily
n.row.Data.daily=nrow(Data.daily)
Data.daily.north=subset(Data.daily,LAT>(-26))   
Data.daily=subset(Data.daily,LAT<=(-26))
if(!n.row.Data.daily==(nrow(Data.daily.north)+nrow(Data.daily)))
{
  par(bg=2)
  plot.new()
  mtext('Split of Data.daily into North ',3,cex=2,col="white")
  mtext(' and South did not work.',3,-2,cex=2,col="white")
  mtext("See if NA's in Latitude",3,-4,cex=2,col="white")
  par(bg="white")
}


#Some Fishery code fixes
Data.monthly=Data.monthly%>%
              mutate(FisheryCode=case_when(FisheryCode%in%c('CL02')~'JANS',
                                           FisheryCode%in%c('SGL','SGL1','SGL2')~'JASDGDL',
                                           FisheryCode%in%c('WL')~'OASC',
                                           FisheryCode%in%c('C127')~'WANCS',
                                           FisheryCode%in%c('WCGL')~'WCDGDL',
                                           FisheryCode%in%c('WCDGDL') & LAT<=(-33)~'JASDGDL',
                                           TRUE~FisheryCode),
                     FisheryCode=case_when(zone=='West' & FisheryCode=='JASDGDL'~'WCDGDL',
                                           TRUE~FisheryCode))
Data.daily=Data.daily%>%
              mutate(FisheryCode=case_when(zone=='West' & FisheryCode=='JASDGDL'~'WCDGDL',
                                           TRUE~FisheryCode))

# C.9.1 Northern and Southern catches
compare.monthly.northern.southern=FALSE
if(compare.monthly.northern.southern)
{
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
  
}


#Compare with Rory's previous summaries
if(!do.sql.extraction)
{
  if(First.run=="YES")
  {
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
  }
  
  #by zone and species
  if(First.run=="YES")
  {
    Read=read.csv(handl_OneDrive("Data/Catch and Effort/Historic/Spec.catch.zone.csv"))
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
  }
}




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


#SECTION D. ---- DATA MANIPULATION - EFFORT INSPECTIONS -----
#note: this part sets functions for inspecting current year data and raise issues
#       to data-entry girls
#       Visually inspect each plot and raise any inconsistent record to data-entry girls
#       This has already been done for all Daily Records prior to current year.

Dummy.yr=sort(unique(Effort.daily$finyear))
par(mfcol=c(1,1),mai=c(1.3,1.3,.2,.2))
Max.list=list(NETLEN=Net.max, BDAYS=30, HOURS=10,
              SHOTS=3, nlines=5, HOOKS=2500,fdays=30)  #max realistic effort var values for daily effort
Inspect.vars.combo=list(GN=c("NETLEN","SHOTS","HOURS"),
                        LL=c("HOOKS","SHOTS","HOURS"))
Max.list.combo=list(GN=Max.list$NETLEN*Max.list$SHOTS*Max.list$HOURS, 
                    LL=Max.list$HOOKS*Max.list$SHOTS*Max.list$HOURS)

if(Inspect.New.dat=="YES")
{
  Inspect.vars=c("NETLEN","BDAYS","HOURS","SHOTS","nlines","HOOKS","fdays")
  
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
    Data=subset(Data,!is.na(km.gn.hours))
    ddd=vector('list',length(Unik.ves))
    for (i in 1:length(Unik.ves))
    {
      da=subset(Data,VESSEL==Unik.ves[i])
      
      if(nrow(da)>0)
      {
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
    hndl.net_line=handl_OneDrive("Analyses\\Catch and effort\\Outputs\\Netlen_nlines\\")
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
    
    pdf(paste(hndl.net_line,"non_splitters.pdf",sep=""))
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
  
  #-- Single effort variable exploration
  handle=paste(handl_OneDrive("Data/Catch and Effort/Check_these_vars/Daily/"),Current.yr,sep="")
  
  fn.inspect.singl.var=function(d,var,ylim,valid.max)
  {
    d=d%>%
      arrange(VESSEL,date)%>%
      group_by(VESSEL)%>%
      mutate(Session=row_number())%>%
      filter(!is.na(!!rlang::sym(var)))
    if(nrow(d)>0)
    {
      d$valid.max=valid.max
      if(var=="HOURS")
      {
        d=d%>%
          group_by(VESSEL)%>%
          mutate(valid.max=mean(!!rlang::sym(var)))
      }
      rm(valid.max)
      ylim[2]=max(ylim[2],max(d[,var]))
      p=d%>%
        ggplot(aes_string(x='Session',y=var))+
        geom_point()+
        geom_hline(aes(yintercept = valid.max),color=2)+
        facet_wrap(~VESSEL,scales='free')+
        ylim(ylim)
      print(p)
      
      var1=var
      if(var1=="HOURS")  var1=c(var,'SHOTS')
      out.dodgy=d%>%filter(!!rlang::sym(var)>2*!!rlang::sym('valid.max'))%>%
        dplyr::select(c("VESSEL","BoatName","date","TSNo","SNo","DSNo",'valid.max',all_of(var1)))%>%
        rename(TripSheetNumber=TSNo,
               DailySheetNumber=DSNo,
               SessionIndex=SNo)%>%
        mutate(Max=valid.max)
      if(nrow(out.dodgy)>0) return(out.dodgy)
    }
    
  }

  for(v in 1:length(Inspect.vars))
  {
    Out=fn.inspect.singl.var(d=Data.daily.original%>%
                               distinct(Same.return.SNo,.keep_all=T)%>%
                               filter(FINYEAR==Current.yr),
                             var=Inspect.vars[v],
                             ylim=c(0,Max.list[[v]]),
                             valid.max=Max.list[[v]])
    if(!is.null(Out))
    {
      write.csv(Out,
                file=paste(handle,"/Check.effort variable beyond range_",Inspect.vars[v],".csv",sep=""),
                row.names=F)
      send.email(TO=Email.data.checks,
                 CC=Email.data.checks2,
                 BCC=Email.FishCube,
                 Subject=paste("Shark validation",Sys.time(),Inspect.vars[v] ,"beyond range",sep=' _ '),
                 Body= "Hi,
                      This effort variable is beyond the maximum (see column Max). Is this a typo or system error?
                      Cheers
                      Matias",
                 Attachment=paste(handle,"/Check.effort variable beyond range_",Inspect.vars[v],".csv",sep=""))
      
    }
  }
  
  #-- Combined effort variables exploration
  fn.inspect.combo.var=function(d,var,ylim,valid.max)
  {
    d=d%>%
      arrange(VESSEL,date)%>%
      group_by(VESSEL)%>%
      mutate(Session=row_number())%>%
      filter(!is.na(!!rlang::sym(var[1])))%>%
      data.frame
    if(nrow(d)>0)
    {
      d$Effort=d[,var[1]]*d[,var[2]]*d[,var[3]]
      
      ylim[2]=max(ylim[2],max(d$Effort,na.rm=T))
      p=d%>%
        ggplot(aes_string(x='Session',y='Effort'))+
        geom_point()+
        geom_hline(yintercept = valid.max,color=2)+
        facet_wrap(~VESSEL,scales='free')+
        ylim(ylim)
      print(p)
      
      out.dodgy=d%>%filter(Effort>valid.max)%>%
        dplyr::select(c("VESSEL","BoatName","date","TSNo","SNo","DSNo",all_of(var)))%>%
        rename(TripSheetNumber=TSNo,
               DailySheetNumber=DSNo,
               SessionIndex=SNo)%>%
        mutate(Max=valid.max)
      if(nrow(out.dodgy)>0) return(out.dodgy)
    }
    
  }
  for(v in 1:length(Inspect.vars.combo))
  {
    Out=fn.inspect.combo.var(d=Data.daily.original%>%
                               distinct(Same.return.SNo,.keep_all=T)%>%
                               filter(FINYEAR==Current.yr),
                             var=Inspect.vars.combo[[v]],
                             ylim=c(0,Max.list.combo[[v]]),
                             valid.max=Max.list.combo[[v]])
    if(!is.null(Out))
    {
      write.csv(Out,
                file=paste(handle,"/Check.effort variable beyond range_",paste(Inspect.vars.combo[[v]],collapse='_'),".csv",sep=""),
                row.names=F)
      send.email(TO=Email.data.checks,
                 CC=Email.data.checks2,
                 BCC=Email.FishCube,
                 Subject=paste("Shark validation",Sys.time(),paste(Inspect.vars.combo[[v]],collapse='_') ,"beyond range",sep=' _ '),
                 Body= "Hi,
                      These effort variables are beyond the maximum (see column Max). Is this a typo or system error?
                      Cheers
                      Matias",
                 Attachment=paste(handle,"/Check.effort variable beyond range_",paste(Inspect.vars.combo[[v]],collapse='_'),".csv",sep=""))
      
    }
  }
  
  #-Net length
  #graphics.off()
  #fun.plot.eff.var(Data.daily.original,Current.yr,"NETLEN",0,15000,Net.max,"YES")
  #Inspect.net.monthly=fun.inspect(Data.monthly.original,FINYEAR.monthly,"Monthly",100,.2,"NETLEN")
 
    #vessel visually identified as problematic
  #Inspect.net.daily=fun.inspect(Data.daily.original,Current.yr,"Daily",100,.2,"NETLEN")
 
  
  #-BDAYS by session
  #graphics.off()
  #fun.plot.eff.var(Data.daily.original,Current.yr,"BDAYS",0,40,15,"YES")

  #-Fdays by month (shouldn't be > 30)
  #graphics.off()
  #fun.plot.fishing.days.month(Data=Data.daily.original,YRS=Current.yr,VAR="fdays",lim1=0,lim2=40,threshold=30)

  
  #-Hours
  #graphics.off()
  #fun.plot.eff.var(Data.daily.original,Current.yr,"HOURS",0,48,24,"YES")
  #check.vesl=c("B 091","E 009") #vessell visually identified as problematic
  #MIN=1; MAX=24
  #Inspect.hours=fun.further.chk(Data.daily.original,Current.yr,"HOURS")%>%
  #                   filter(Ok.var=='Problem')

  
  #-Shots   
 # graphics.off()
  #Table.Shots.year=table(Effort.daily$finyear,Effort.daily$shots)
  #fun.plot.eff.var(Data.daily.original,Current.yr,"SHOTS",0,10,2,"YES")
  # check.vesl=c("G 297")  #vessell visually identified as problematic
  # MIN=1; MAX=2
  # Inspect.shots=fun.further.chk(Data.daily.original,Current.yr,"SHOTS")%>%
  #   filter(Ok.var=='Problem')%>%
  #   arrange(VESSEL,date)


  
  #-Shots X Hours     
  #graphics.off()
  # fun.plot.combo(Data.daily.original,Current.yr,"SHOTS","HOURS",0,48,26,"YES")
  #  check.vesl=c("E 035","E 009","G 297") 
  #  Check.shot.hours=subset(Data.daily.original,FINYEAR==Current.yr & VESSEL%in%check.vesl)%>%
  #    distinct(Same.return.SNo,date,VESSEL,SNo,DSNo,TSNo,SHOTS,HOURS)%>%
  #    mutate(Shots.hours=SHOTS*HOURS)%>%
  #    filter(Shots.hours>26)%>%
  #    arrange(VESSEL,date)
  #  if(nrow(Check.shot.hours)>0)
  #  {
  #    write.csv(Check.shot.hours,paste(handle,"/Inspect.shot.hours.csv",sep=""),row.names = F)
  #    send.email(TO=Email.data.checks,
  #               CC=Email.FishCube,
  #               Subject="TDGDLF data validation check SHOTS & HOURS",
  #               Body= "For the attached file, could you please check that there are no typos in the ‘SHOTS’ &  'HOURS' columns. 
  #                    The values are way outside the range reported by these fishers in other shots. 
  #                   Cheers
  #                   Matias",
  #               Attachment=paste(handle,"/Inspect.shots.csv",sep=""))
  #  }
  # 
  # MIN=1; MAX=24
  # graphics.off()
  # Inspect.hours.shots=fun.further.chk.combo(Data.daily.original,Current.yr,"SHOTS","HOURS")
  
  
  #-Effort (km.gn.hours)
  #graphics.off()
  #Explr.ves=fun.plot.km.gn.hours(Data.daily.original,Current.yr,"SHOTS","HOURS","NETLEN","YES")
 
    

}

#inspect overall daily effort
if(Inspect.New.dat=='YES')
{
  dodgy.effort=Inspect.vars.combo
  fn.inspect.dodgy.daily=function(d,var,valid.max)
  {
    names(d)=tolower(names(d))
    var=tolower(var)
    d=d%>%
      group_by(vessel,same.return.sno)%>%
      mutate(Session=row_number())%>%
      filter(!is.na(!!rlang::sym(var[1])))%>%
      data.frame
    if(nrow(d)>0)
    {
      d$Effort=d[,var[1]]*d[,var[2]]*d[,var[3]]
      out.dodgy=d%>%filter(Effort>valid.max)%>%
        dplyr::select(c("vessel","year","month","same.return.sno",all_of(var)))%>%
        mutate(Max=valid.max)%>%
        arrange(vessel,year,month)
      if(nrow(out.dodgy)>0) return(out.dodgy)
    }
  }
  for(v in 1:length(Inspect.vars.combo))
  {
    dodgy.effort[[v]]=fn.inspect.dodgy.daily(d=Data.daily%>%
                                               distinct(Same.return.SNo,.keep_all=T),
                                             var=Inspect.vars.combo[[v]],
                                             valid.max=Max.list.combo[[v]])
  }
}
#note: most issues due to fishers confusing number of 'shots' with number of 'net split'
#      Hence, not using 'shots' in effort calculations. See '#2. Calculate gillnet effort'


#SECTION E. ---- DATA MANIPULATION - EFFORT CORRECTIONS ----      

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
      ifelse(METHOD=="GN" & VESSEL=="F 471" & SHOTS==100,1,
      ifelse(METHOD=="GN" & VESSEL=="F 520" & SHOTS==200,1,
      ifelse(METHOD=="GN" & VESSEL=="M 141" & SHOTS==2100,1,
      ifelse(METHOD=="GN" & VESSEL=="F 428" & SHOTS==40,1,
      SHOTS.c)))))))))))))



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
Effort.daily$shots.c=with(Effort.daily,
                ifelse(method=="GN" & vessel%in%c("E 059",'E 067','F 517') & shots==11| is.na(shots),1,shots.c))


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
#update Reporter
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
if(nrow(a)>0)
{
  aggregate(SHOTS~VESSEL, a, Mode)
  aggregate(NETLEN~VESSEL, a, Mode)
}


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
write.csv(round(100*Table.Eff.Reporter[1]/sum(Table.Eff.Reporter),1),
          handl_OneDrive("Analyses/Catch and effort/Outputs/Paper/Percent.GN.effort.fixed.csv"))

# Effort for Alex's ASL model
if(do.Alexs=="YES")
{
  Effort.alex=subset(Effort.monthly,YEAR>2005 & !FINYEAR=="2005-06" & METHOD=="GN" & LAT<=(-26))
  Effort.alex.daily= subset(Effort.daily,method=="GN")
  Effort.alex.daily=Effort.alex.daily[,-match(c("LAT","LONG","date","block10","year","year.c","month"),names(Effort.alex.daily))]
}


#E.3 Add effort variables to FishCUBE monthly as requested by Vero
if(!do.sql.extraction)
{
  FishCUBE.monthly.effort=Effort.monthly%>%
    filter(Same.return%in%FishCUBE.monthly$Same.return)%>%
    distinct(Same.return,.keep_all=T)%>%
    dplyr::select(Same.return,HOURS.c,HOOKS,SHOTS.c,NETLEN.c)
  
  FishCUBE.monthly=FishCUBE.monthly%>%
    left_join(FishCUBE.monthly.effort,by='Same.return')
  
}


#SECTION F.1. ---- EXTRACT QUANTITIES ---- 

#Check for spatial expansion

  #Annual
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


  #Proportion of blocks fished per month for each year
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

#1. Define shark fisheries in the north and south       (see FishCubeCode_Shark.fisheries)
Effort.monthly=Effort.monthly%>%
       mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                      zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                      zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                      zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                      is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                      zone=='Zone1' & FisheryCode=='WCGL'~'JASDGDL',
                                      FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                      zone%in%c('Zone1','Zone2') & METHOD%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                      zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                      zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                      zone%in%c('Closed','Joint','North') & METHOD%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                      zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                      zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                      TRUE~"non.shark.fishery"))
Effort.daily=Effort.daily%>% 
  mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                 zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                 zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                 zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                 is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                 zone=='Zone1' & FisheryCode%in%c('WCDGDL','WCGL')~'JASDGDL',
                                 FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                 zone%in%c('Zone1','Zone2') & method%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                 zone%in%c('Closed','Joint','North') & method%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                 zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                 zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                 TRUE~"non.shark.fishery"))

Data.monthly=Data.monthly%>%
  mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                 zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                 zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                 zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                 is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                 zone=='Zone1' & FisheryCode%in%c('WCDGDL','WCGL')~'JASDGDL',
                                 FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                 zone%in%c('Zone1','Zone2') & METHOD%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                 zone%in%c('Closed','Joint','North') & METHOD%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                 zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                 zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                 TRUE~"non.shark.fishery"))

Data.monthly.north=Data.monthly.north%>%
  mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                 zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                 zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                 zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                 is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                 zone=='Zone1' & FisheryCode%in%c('WCDGDL','WCGL')~'JASDGDL',
                                 FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                 zone%in%c('Zone1','Zone2') & METHOD%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                 zone%in%c('Closed','Joint','North') & METHOD%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                 zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                 zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                 TRUE~"non.shark.fishery"))

Data.daily=Data.daily%>%
  mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                 zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                 zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                 zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                 is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                 zone=='Zone1' & FisheryCode%in%c('WCDGDL','WCGL')~'JASDGDL',
                                 FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                 zone%in%c('Zone1','Zone2') & METHOD%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                 zone%in%c('Closed','Joint','North') & METHOD%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                 zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                 zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                 TRUE~"non.shark.fishery"))

Data.daily.north=Data.daily.north%>%
  mutate(Shark.fishery=case_when(zone=='Joint' & FisheryCode%in%c('CL02','OT','C051','JANS','NCS','WANCS')~'JANS',
                                 zone=='North' & FisheryCode%in%c('C127','OT','WANCS','NCS','C051','JANS')~'WANCS',
                                 zone=='Closed' & FisheryCode%in%c('OT')~'WANCS',
                                 zone=='West' & FisheryCode%in%c('WCGL','C070','WCDGDL')~'WCDGDL',
                                 is.na(zone) & FisheryCode=='WCGL'~'WCDGDL',
                                 zone=='Zone1' & FisheryCode%in%c('WCDGDL','WCGL')~'JASDGDL',
                                 FisheryCode%in%c("SGL1","SGL2","SGL",'JASDGDL')~'JASDGDL',
                                 zone%in%c('Zone1','Zone2') & METHOD%in%c('GN','LL') & FisheryCode=='OASC'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=='*'~'OASC',
                                 zone%in%c('Zone1','Zone2') & FisheryCode=="WL"~'OASC',
                                 zone%in%c('Closed','Joint','North') & METHOD%in%c('GN','LL') & FisheryCode=='OANCGCWC'~'OANCGCWC',
                                 zone%in%c('West') & FisheryCode=='*'~'OANCGCWC',
                                 zone%in%c('West','Zone1') & FisheryCode=='OT'~'OANCGCWC',
                                 TRUE~"non.shark.fishery"))

Exprt.shark.fisheries.code=rbind(Data.daily.north%>%distinct(zone,FisheryCode,Shark.fishery)%>%mutate(DATA='daily.north'),
                                 Data.daily%>%distinct(zone,FisheryCode,Shark.fishery)%>%mutate(DATA='daily'),
                                 Data.monthly.north%>%distinct(zone,FisheryCode,Shark.fishery)%>%mutate(DATA='monthly.north'),
                                 Data.monthly%>%distinct(zone,FisheryCode,Shark.fishery)%>%mutate(DATA='monthly'))%>%
                            distinct(zone,FisheryCode,Shark.fishery)

write.csv(Exprt.shark.fisheries.code%>%left_join(FisheryCodeTable%>%
                                                   distinct(FishCubeCode,.keep_all = T)%>%
                                                   dplyr::select(-SASCode, CurrentAt2017),
                                                 by=c("FisheryCode"="FishCubeCode")),
          "Shark.fisheries_and_non.shark_fisheries.csv",row.names = F)

#2. Calculate gillnet effort 

#note:  . RORY DOESN'T USE SHOTS FOR KM.GN.DAYS as this effort measure is only concerned on km of nets used in a day                             
#       . For daily km gn hours, shots are not used because fishers report the total hours fished per session,
#           whether doing 1, 2 or more shots per session
#       . For monthly km gn hours, shots should not be used because some fishers report combined hours

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

hnl.msh=handl_OneDrive("Analyses/Catch and effort/")

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
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+finyear,
                                  data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~ID+vessel+zone+finyear,
                                  data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
  
  Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+Fishing_yr+block10,
                                  data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
}

if(Use.Date=="YES")
{
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~date+vessel+zone+finyear,
                                  data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T) 
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~date+vessel+zone+finyear,
                                data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
  
  Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~date+vessel+zone+Fishing_yr+block10,
                                data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
}

Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~zone+finyear,data=Attach.Effort.daily.c,sum,na.rm=T)
Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~zone+finyear,data=Attach.Effort.daily,sum,na.rm=T)
Jodie.Effort.daily.c_block10=aggregate(Km.Gillnet.Days.c~vessel+zone+Fishing_yr+block10,data=Jodie.Effort.daily.c_block10,sum,na.rm=T)


        #km gn hours
Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~finyear+vessel+ID+zone,
                                    data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~finyear+vessel+ID+zone,
                                  data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)

Jodie.Effort.daily.hrs.c_block10=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+ID+zone+block10,
                                           data=subset(Effort.daily,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)

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


Jodie.Effort=Jodie.Effort.daily.c_block10%>%full_join(Jodie.Effort.daily.hrs.c_block10,
                                                by=c("vessel","zone","Fishing_yr","block10"))




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
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,
                                  data=subset(Effort.monthly,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,
                                  data=subset(Effort.monthly,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)

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

      #sum to aggregate by zone and year
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~zone+FINYEAR,data=Attach.Effort.monthly.c,sum,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~zone+FINYEAR,data=Attach.Effort.monthly,sum,na.rm=T)


    #km gn hours
      #max to remove duplicates
Attach.Effort.monthly.hrs.c=aggregate(Km.Gillnet.Hours.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,
                                      data=subset(Effort.monthly,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)
Attach.Effort.monthly.hrs=aggregate(Km.Gillnet.Hours.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,
                                    data=subset(Effort.monthly,!Shark.fishery=='non.shark.fishery'),max,na.rm=T)


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
  SUBSET=Effort.daily%>%
      filter(finyear%in%THESE.YRS & !Shark.fishery=='non.shark.fishery')
  
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
  
  
  #3.1.2. Calculate longline gillnet effort equivalent (McAuley & Simpfendorfer 2003)
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


    #3.2.2. Calculate longline gillnet effort equivalent  (McAuley et al 2003)                         
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
Effort.monthly.NSF=subset(Effort.monthly,zone%in%c("Closed","Joint","North") & 
                          Shark.fishery%in%c('JANS','WANCS'),
                     select=c(FINYEAR,MONTH,BLOCKX,zone,VESSEL,METHOD,Same.return,
                              FDAYS,BDAYS.c,HOURS.c,HOOKS,SHOTS.c,NETLEN.c,nlines,
                              Km.Gillnet.Days.c,Km.Gillnet.Hours.c,
                              FisheryCode))

Effort.daily.NSF=subset(Effort.daily,zone%in%c("Closed","Joint","North") & 
                          Shark.fishery%in%c('JANS','WANCS'),
                    Select=c(date,finyear,month,blockx,zone,vessel,method,Same.return,Same.return.SNo,
                             fdays,bdays.c,hours.c,hooks,shots.c,netlen.c,nlines,
                             Km.Gillnet.Days.c,Km.Gillnet.Hours.c,
                             FisheryCode))

Effort.monthly.NSF$hook.days=with(Effort.monthly.NSF,BDAYS.c*HOOKS)
Effort.monthly.NSF$hook.hours=with(Effort.monthly.NSF,BDAYS.c*HOOKS*HOURS.c)

Effort.daily.NSF$hook.days=with(Effort.daily.NSF,hooks)
Effort.daily.NSF$hook.hours=with(Effort.daily.NSF,hooks*hours.c)


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

#issue: check daily 2007-08, too high effort
Total.effort_NFS=rbind(Attach.Effort.monthly.c_NSF,Attach.Effort.daily.c_NSF)%>% 
      group_by(finyear)%>%
      summarise(hook.days=sum(hook.days),
                hook.hours=sum(hook.hours))%>%
      mutate(hook.days=hook.days/1000,
             hook.hours=hook.hours/1000)
names(Total.effort_NFS)=c("FINYEAR","Hook days","Hook hours")


#Compare Rory's effort and current script's effort
if(Inspect.New.dat=="YES")
{
  Results.pre.2013=read.csv(handl_OneDrive("Data/Catch and Effort/Historic/Historic.res.csv"))
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


#Add corresponding effort equivalents
fn.quiv.eff=function(ktch,efF,SP,FORM.pred)
{
  efF=efF%>%
    dplyr::select(-c(block10))
  if(is.na(match('hook.hours',names(efF))))
  {
    efF$hook.days=with(efF,hooks)
    efF$hook.hours=with(efF,hooks*hours.c)
  }
  
  ktch=ktch%>%filter(SPECIES%in%SP)
  d=left_join(ktch,efF,by='Same.return.SNo')%>%
    dplyr::select(Same.return.SNo,METHOD,day,YEAR.c,MONTH,BLOCKX,SPECIES,LIVEWT.c,
           Km.Gillnet.Hours.c,hook.hours,Km.Gillnet.Days.c,hook.days)%>%
    filter(METHOD%in%c('GN','LL'))%>%
    mutate(m.Gillnet.Hours.c=1000*Km.Gillnet.Hours.c,
           m.Gillnet.Days.c=1000*Km.Gillnet.Days.c,
           cpue.hour=ifelse(METHOD=="GN",LIVEWT.c/m.Gillnet.Hours.c,   #have conversion m to hooks
                     ifelse(METHOD=="LL",LIVEWT.c/hook.hours,NA)),
           cpue.day=ifelse(METHOD=="GN",LIVEWT.c/m.Gillnet.Days.c,
                     ifelse(METHOD=="LL",LIVEWT.c/hook.days,NA)))%>%
           filter(!cpue.hour=="Inf")
  Rango.LL=unique(subset(d,METHOD=='LL')$BLOCKX)
  Rango.GN=unique(subset(d,METHOD=='GN')$BLOCKX)
  
  d=d%>%filter(BLOCKX%in%intersect(Rango.LL,Rango.GN))%>%
        mutate(METHOD=as.factor(METHOD),
               MONTH=as.factor(MONTH),
               YEAR.c=as.factor(YEAR.c),
               BLOCKX=as.factor(BLOCKX))
  
  mod.day=glm(formula(paste('log(cpue.day)',FORM.pred,sep="~")),data=d,family=gaussian)
  mod.hour=glm(formula(paste('log(cpue.hour)',FORM.pred,sep="~")),data=d,family=gaussian)

  NEWD=data.frame(METHOD=c("GN","LL"),
                  YEAR.c=factor(names(sort(-table(d$YEAR.c))[1]),levels(d$YEAR.c)),
                  MONTH=factor(names(sort(-table(d$MONTH))[1]),levels(d$MONTH)),
                  BLOCKX=factor(names(sort(-table(d$BLOCKX))[1]),levels(d$BLOCKX)))
  a=predict(mod.day,newdata=NEWD, type="response",se.fit=T)
  NEWD$cpue.day=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
  a=predict(mod.hour,newdata=NEWD, type="response",se.fit=T)
  NEWD$cpue.hour=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
  
  #ratios
  LL.to.GN=data.frame(t(with(NEWD,c(day=cpue.day[2]/cpue.day[1],hour=cpue.hour[2]/cpue.hour[1]))))
    
  return(LL.to.GN)

}
  #North
top.sp=aggregate(LIVEWT.c~SPECIES,subset(Data.daily.north,!Shark.fishery=='non.shark.fishery'),sum)%>%
          arrange(-LIVEWT.c)
top.spID=cumsum(top.sp$LIVEWT.c)/sum(top.sp$LIVEWT.c)
top.spID=1:which.min(abs(top.spID - .95))
LL.to.GN.North=fn.quiv.eff(ktch=subset(Data.daily.north,!Shark.fishery=='non.shark.fishery'),
                           efF=subset(Effort.daily.NSF,!Shark.fishery=='non.shark.fishery')%>%distinct(Same.return.SNo,.keep_all = T),
                            SP=top.sp$SPECIES[top.spID],FORM.pred=paste('METHOD','YEAR.c','BLOCKX',sep="+"))

  #South
LL.to.GN.South=fn.quiv.eff(ktch=subset(Data.daily,!Shark.fishery=='non.shark.fishery'),
                           efF=subset(Effort.daily,!Shark.fishery=='non.shark.fishery')%>%distinct(Same.return.SNo,.keep_all = T),
                            SP=unlist(TARGETS),FORM.pred=paste('METHOD','YEAR.c','MONTH','BLOCKX',sep="+"))




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
    #monthly TDGDLF
Data.monthly.GN=subset(Data.monthly,METHOD=="GN" & Estuary=="NO" & !Shark.fishery=='non.shark.fishery')
Data.monthly.LL=subset(Data.monthly,METHOD=="LL" & Estuary=="NO"& !Shark.fishery=='non.shark.fishery')
Data.monthly.other=subset(Data.monthly,!(METHOD=="LL"|METHOD=="GN")
                    |((METHOD=="LL"|METHOD=="GN") & Estuary=="YES"))

    #daily TDGDLF
Data.daily.GN=subset(Data.daily,METHOD=="GN" & Estuary=="NO" & !Shark.fishery=='non.shark.fishery')
Data.daily.LL=subset(Data.daily,METHOD=="LL" & Estuary=="NO" & !Shark.fishery=='non.shark.fishery')
Data.daily.other=subset(Data.daily,!(METHOD=="LL"|METHOD=="GN")
                          |((METHOD=="LL"|METHOD=="GN") & Estuary=="YES"))

#Check use of different gears by year and zone for main species
fn.check.method.contrib=function(D,SPEC)
{
  hndl.gear=handl_OneDrive("Analyses/Catch and effort/Outputs/")
  a=subset(D,  Estuary=="NO" & LIVEWT.c>0 & LAT<=(-26) & SPECIES==SPEC)
  a$METHOD1=with(a,ifelse(!METHOD%in%c("GN","HL","LL"),"Other",METHOD))
  Rcrds.yr.mthd=with(a,table(FINYEAR,METHOD1))
  Rcrds.yr.mthd=Rcrds.yr.mthd/rowSums(Rcrds.yr.mthd)
  Rcrds.yr.mthd=t(Rcrds.yr.mthd)
  
  COLS=c("darkolivegreen4","black","brown","grey60")  
  fn.fig(paste(hndl.gear,"TDGDLF_gear_by_year.",a$SNAME[1],sep=""),2400, 2400)
  par(mfcol=c(1,1),mai=c(1.1,.85,.2,.1),oma=c(.1,.1,1,.1),xpd=TRUE,las=1)
  
  #Zones combined
  barplot(Rcrds.yr.mthd,legend.text=rownames(Rcrds.yr.mthd),col=COLS,
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
  barplot(x,legend.text=rownames(x),col=COLS,
          args.legend=list(x=ncol(x)*1.2,
                           y=max(x)*.975,cex=1.15,bty='n',horiz=T,
                           border='black',yjust=0))
  mtext(dimnames(Rcrds.GN.LL)$zone[1],3)
  
  #Zn1
  x=Rcrds.GN.LL[,,2]/rowSums(Rcrds.GN.LL[,,2])
  x=t(x)
  barplot(x,col=COLS)
  mtext(dimnames(Rcrds.GN.LL)$zone[2],3)
  
  #Zn2
  x=Rcrds.GN.LL[,,3]/rowSums(Rcrds.GN.LL[,,3])
  x=t(x)
  barplot(x,col=COLS)
  mtext(dimnames(Rcrds.GN.LL)$zone[3],3)
  
  mtext("Financial year",1,0.5,outer=T,cex=1.5)
  mtext("Proportion of records",2,-3.5,outer=T,las=3,cex=1.5)
  
  dev.off()
  
}
fn.check.method.contrib(D=Data.monthly%>%filter(!(is.na(zone)|zone=='')),SPEC=17001)
fn.check.method.contrib(D=Data.monthly%>%filter(!(is.na(zone)|zone=='')),SPEC=17003)
fn.check.method.contrib(D=Data.monthly%>%filter(!(is.na(zone)|zone=='')),SPEC=18003)
fn.check.method.contrib(D=Data.monthly%>%filter(!(is.na(zone)|zone=='')),SPEC=18007)


#Data for Hammerhead Listing 
if(do.Hammerheads=="YES")
{
  HHEads=c(19000,HammerheadSpecies)
  Sel.HH=c("SPECIES","FINYEAR","MONTH","METHOD","VESSEL","BLOCKX","LAT","LONG","zone","LIVEWT.c")
  Hammerheads=subset(Data.monthly,SPECIES%in% HHEads,select=Sel.HH)
}
 

#Data for Jeff Norris
if (do.Jeffs=="YES") Jeff.Norris=subset(Data.monthly,SPECIES%in%Scalefish.species & Bioregion=="SC")  #Data requested by Jeff Norris


#add nfish
Data.monthly.GN$BLOCKX=as.numeric(Data.monthly.GN$BLOCKX)
Data.monthly.GN=Data.monthly.GN%>%left_join(Daily.nfish.agg,by=c("Same.return",
           "FINYEAR","YEAR","YEAR.c","MONTH","zone","BLOCKX","VESSEL","METHOD","SPECIES","Estuary"))
# Data.monthly.GN=merge(Data.monthly.GN,Daily.nfish.agg,by=c("Same.return",
#   "FINYEAR","YEAR","YEAR.c","MONTH","zone","BLOCKX","VESSEL","METHOD","SPECIES","Estuary"),all.x=T)

if (do.Jeffs=="YES")
{
  Jeff.Norris=merge(Jeff.Norris,Daily.nfish.agg,by=c("Same.return","FINYEAR","YEAR","YEAR.c","MONTH","zone",
                                                     "BLOCKX","VESSEL","METHOD","SPECIES","Estuary"),all.x=T)
}


  #4.3. check spatial expansion
BLKS.YEAR=Expand.fun(Data.monthly.GN)
Spatial.expan=round(100*BLKS.YEAR/N.blocks,0)
Spatial.expan.north=round(100*Expand.fun(Data.monthly.north)/N.blocks,0)
Effort.expan=Effort1.fun(Data.monthly.GN)
Effort.expan.north=Effort1.fun(Data.monthly.north)


  #4.4. Define effective area effort                                   #Rory's rule 5a-5k
Data.monthly.GN$Boundary.blk=with(Data.monthly.GN,
          ifelse(BLOCKX%in%Boundary.blk & !FINYEAR%in%Daily.l.years,"Y","N"))



#SECTION F.2. ---- EXPORT DATA FOR ASSESSMENT AND CPUE STANDARDISATION ----

setwd(handl_OneDrive("Analyses/Data_outs"))

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
fn.crap=function(Krap,neim)
{
  neim.id=match(Krap,neim)
  neim.id=subset(neim.id,!is.na(neim.id))
  return(neim.id)
}
Exprt.list=list(
      Annual.total.eff.days=Total.effort.days.monthly,
      Annual.total.eff.hours=Total.effort.hours.monthly,
      Annual.zone.eff.hours=Total.effort.zone.hours.monthly,
      Annual.zone.eff.days=Total.effort.zone.days.monthly,
      Annual.total.eff_NSF=Total.effort_NFS,  
      Effort.monthly=Effort.monthly%>%
                        filter(!zone%in%c("Closed","Joint","North"))%>%
                        dplyr::select(-crap.ef[1:length(crap.ef)]),
      Effort.daily=Effort.daily%>%
                        filter(!zone%in%c("Closed","Joint","North"))%>%
                        dplyr::select(-crap.ef[1:length(crap.ef)]),
      Effort.monthly.north=Effort.monthly.NSF,
      Effort.daily.north=Effort.daily.NSF%>%
                          dplyr::select(-crap.ef[1:length(crap.ef)]),
      Mesh.monthly=Mesh.monthly,
      Mesh.size=Mesh.size,
      Data.current.Sofar=Data.current.Sofar, 
      Suite=data.frame(Suite=Suite),
      LL.equiv.Eff.days.zone=LL.equiv.Eff.days.zone,
      Equivalent.LL_to_GN_South=LL.to.GN.South,
      Equivalent.LL_to_GN_North=LL.to.GN.North,
      Data.monthly=Data.monthly[,-fn.crap(crap,names(Data.monthly))],    
      Data.monthly.north=Data.monthly.north[,-fn.crap(crap,names(Data.monthly.north))],
      Data.monthly.GN=Data.monthly.GN[,-fn.crap(crap,names(Data.monthly.GN))],
      Data.monthly.LL=Data.monthly.LL[,-fn.crap(crap,names(Data.monthly.LL))],
      Data.daily=Data.daily[,-fn.crap(crap.daily,names(Data.daily))]%>%
                  left_join(List.of.daily.vessel.boat,by=c('VESSEL'='vessel','Same.return.SNo')),
      Data.daily.north=Data.daily.north[,-fn.crap(crap.daily,names(Data.daily.north))]%>%
                  left_join(List.of.daily.vessel.boat,by=c('VESSEL'='vessel','Same.return.SNo')),
      Data.daily.original=Data.daily.original,
      Data.daily.GN=Data.daily.GN[,-fn.crap(crap.daily,names(Data.daily.GN))]%>%
                  left_join(List.of.daily.vessel.boat,by=c('VESSEL'='vessel','Same.return.SNo')),
      Data.daily.LL=Data.daily.LL[,-fn.crap(crap.daily,names(Data.daily.LL))]%>%
                  left_join(List.of.daily.vessel.boat,by=c('VESSEL'='vessel','Same.return.SNo'))
      )

if(exists('TEPS.current'))Exprt.list$TEPS.current=TEPS.current
if(exists('PRICES'))Exprt.list$PRICES=PRICES
do.dis=FALSE
if(Monthly.log=='ceMonthlyLog' & do.dis)  
{
  dont.export=grep(paste(c('Data.monthly','Effort.monthly'),collapse='|'),names(Exprt.list))
  Exprt.list=Exprt.list[-dont.export]
}
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


#SECTION G.1. ---- SPATIO TEMPORAL CATCH AND EFFORT DISTRIBUTION ------
if(do.spatio.temporal.ktch.effort=="YES")
{
  source(handl_OneDrive("Analyses/Catch and effort/Git_catch.and.effort/do.spatio.temporal.ktch.effort.R"))
}

if(do.tdgdlf.effort.density=="YES")
{
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  p=list(Gillnet=Data.daily%>%
           filter(zone %in%c('West','Zone1','Zone2') & 
                    fishery%in%c('JASDGDL','WCDGDL','*') & 
                    METHOD%in%c('GN'))%>%
           mutate(LAT=as.numeric(LatFC),
                  LONG=as.numeric(LongFC))%>%
           distinct(Same.return.SNo,FINYEAR,LAT,LONG,METHOD)%>%
           filter(!is.na(LAT)),
         Longline=Data.daily%>%
           filter(zone %in%c('West','Zone1','Zone2') & 
                    fishery%in%c('JASDGDL','WCDGDL','*') & 
                    METHOD%in%c('LL'))%>%
           mutate(LAT=as.numeric(LatFC),
                  LONG=as.numeric(LongFC))%>%
           distinct(Same.return.SNo,FINYEAR,LAT,LONG,METHOD)%>%
           filter(!is.na(LAT)))
  
  fun.mup=function(d,LAB)
  {
    plt=ggplot(data = world) +
      geom_sf(color = "black", fill = "darkorange4") +
      coord_sf(xlim =c(112.5,130) , ylim = c(-36,-25), expand = T) +
      xlab("") + ylab("")+
      geom_point(data=d,aes(x=LONG, y=LAT),size=1.5,fill='black',shape=21,alpha=.4)+
      stat_density_2d(data=d, geom = "polygon", contour = TRUE,
                      aes(x=LONG, y=LAT, fill = after_stat(level)), colour = "black",
                      bins = 10,alpha = 0.7)+
      scale_fill_distiller(palette = "Blues", direction = 1) +
      ggtitle(LAB)+
      theme(legend.position = 'none',
            plot.title = element_text(size=22),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 18),
            plot.margin=unit(c(0,0.1,0,0),"cm"))+
      geom_contour(data = Bathymetry, 
                   aes(x=V1, y=V2, z=V3),
                   breaks=c(-100,-200),linetype="dashed",colour="grey10")
    
    return(plt)
  }
  
  
  for(i in 1:length(p)) p[[i]]=fun.mup(d=p[[i]],LAB=names(p)[[i]])
  
  fig=ggarrange(plotlist=p,
                nrow=2,ncol=1,
                common.legend=FALSE)
  annotate_figure(fig,
                  bottom = text_grob("Longitude",size = 20),
                  left = text_grob("Latitude",size = 20,rot = 90))
  ggsave(paste(handl_OneDrive("Analyses/Catch and effort/Outputs/Spatio.temporal_Effort"),
               "Map_TDGDLF_effort density by method.tiff",sep="/"),
         width = 10,height = 14,compression = "lzw")
}

#SECTION G.2. ---- DATA REQUESTS ------
if(do.data.requests=="YES")
{
  source(handl_OneDrive("Analyses/Catch and effort/Git_catch.and.effort/do.data.requests.R"))
}
  

#SECTION H. ---- EXPORT TOTAL CATCH FOR REFERENCE POINT ANALYSIS ------
if(do.Ref.Points=="YES")
{
  for (i in 1:length(TARGETS))write.csv(STORE.Tot.Catch[[i]],
    paste(handl_OneDrive("Analyses/Reference Points/Tot.c."),sp[i],".csv",sep=""),row.names=F)  
}

#SECTION I. ---- DROPPED CODE ------
if(do.exploratory=="YES")
{
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Exploratory"))
  
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
  data(worldLLhigh)
  OZ.long=c(112,155)
  OZ.lat=c(-46,-10)
  tiff(file="Figure 3.Effort_map.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  
  par(mfrow=c(3,3),mai = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
  plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
          col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
  text(133,-21.5,("Australia"),col="black", cex=2)
  text(118.7,-32,("Perth"),col="black", cex=1.1)
  mtext(Lat.exp,side=2,line=0.4,las=3,cex=1.3)
  mtext(Lon.exp,side=1,line=0.6,cex=1.3)
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
      mtext(Lat.exp,side=2,line=3,las=3,cex=1.25)
      mtext(Lon.exp,side=1,line=3,cex=1.25)    
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
        mtext(Lat.exp,side=2,line=3,las=3,cex=1.25)
        mtext(Lon.exp,side=1,line=3,cex=1.25)    
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
        mtext(Lat.exp,side=2,line=3,las=3,cex=1.25)
        mtext(Lon.exp,side=1,line=3,cex=1.25)    
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
  mtext(Lat.exp,side=2,line=-1,las=3,cex=1.1,outer=T)
  mtext(Lon.exp,side=1,line=-1,cex=1.1,outer=T)
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
    mtext(Lat.exp,side=2,line=-1,las=3,cex=1.1,outer=T)
    mtext(Lon.exp,side=1,line=-1,cex=1.1,outer=T)
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
#source(handl_OneDrive("Analyses/Catch and effort/Git_catch.and.effort/do.dropped.code.R"))