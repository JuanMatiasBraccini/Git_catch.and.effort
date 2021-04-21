#-------SCRIPT FOR ANALYSIS VESSEL CHARACTERISTIC SURVEY FOR TDGDLF ---#

library('stringr')
library('lubridate')

# DATA SECTION
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

setwd(handl_OneDrive("Data/Fishing power"))

TDGDLF.vessels=read.csv("TDGDLF.vessels.csv",stringsAsFactors=F)
TDGDLF.skippers=read.csv("TDGDLF.SKIPERS.csv",stringsAsFactors=F)
Vessel.survey=read.csv("VesselGearSurveyData.csv",stringsAsFactors=F)



# PROCEDURE SECTION
TDGDLF.survey=subset(Vessel.survey,BOATREGO%in%TDGDLF.vessels$VESSEL)
TDGDLF.survey=subset(TDGDLF.survey,!DATE=="")
TDGDLF.survey$DATE=as.Date(TDGDLF.survey$DATE,"%d/%m/%Y")
TDGDLF.survey$year=year(TDGDLF.survey$DATE)
TDGDLF.survey$month=month(TDGDLF.survey$DATE)

TDGDLF.survey$year=ifelse(TDGDLF.survey$year==1903,1993,TDGDLF.survey$year)
Table.year=table(TDGDLF.survey$year)
Table.vessels=table(TDGDLF.survey$BOATREGO)
Table.skippers=table(TDGDLF.survey$SKIPPER)


#Find certain fishers
Zn2.fishers=c("MANSTEAD","TONKIN","ISARAYEPPOS","HENDERSON",
              "SOUMELIDIS","SHARP","GILBERT","KENNEDY","PETERS","BROWN",
              "BRADLEY","TINDALL","WHETSONE","WHETSTONE","STEEL","STEELE","BLACK",
              "BRANDERHORST","TOUMAZOS","TRIANTAFYLLOU","GOODALL")
Zn1.fishers=c("COOKE","SCIMONE","SOULOS","MILES","COCKMAN","ROBB","WARRILOW",
              "REAY","SOULIS","ADAMS","SELL")
WC.fishers=c("PARKER","AITCHISON","CARR","MCWHIRTER","BUCKERIDGE")
Fisher.lookup=c(Zn2.fishers,Zn1.fishers,WC.fishers)

include <- function (theList, toMatch) unique (grep(paste(toMatch,collapse="|"),theList, value=TRUE))
FISHERS=include(names(Table.skippers),Fisher.lookup)


TDGDLF.survey$GROSSTON=with(TDGDLF.survey,ifelse(GROSSTON=="750KG",.75,ifelse(GROSSTON=="A",NA,GROSSTON)))
TDGDLF.survey$GROSSTON=as.numeric(TDGDLF.survey$GROSSTON)

#drop irrelevant variables
irrelevant.vars=c("VIDEOM","VIDEOT","TRBSL2","TRBSB2","TRNT2","TRBSL4","TRNT4","TRHRL4","TRBSL5",
               "TRBSB5","TRNT5","TRHRL5",
               "TRBSL6","TRBSB6","TRNT6","TRHRL6","ROXM","ROXT","ROXYR",
               "TrapOtherTypeCount","TrapOtherTypeSpecifications",
               "TrapOtherTypeLength","TrapOtherTypeWidth","TrapOtherTypeHeight",
               "POTBAT3","POTBAT4","BEEHIVE","OTHERTYP","OTHERNUM","TRRIG1","TRRIG2",
               "TRRIGOTH","TRBSL1","TRBSB1","TRNT1","TRHRL1","TRBSL2","TRBSB2","TRNT2","TRHRL2","TRBSL3","TRBSB3",
               "TRNT3","TRHRL3","TRBSL4","TRBSB4","TRNT4","TRHRL4","TRBSL5","TRBSB5","TRNT5","TRHRL5",
               "TRBSL6","TRBSB6","TRNT6","TRHRL6",
               "TROPNO","LNHRNO","LNHRAV",
               "LNHRLN","LNRDNO","LNRDAV","LNRDLN","LNGRNO","LNGRAV","LNGRLN","LNHYNO","LNHYAV","LNHYLN",
               "LNHYANO","LNHYAAV","LNHYALN","LNELNO","LNELAV","LNELLN","LNELANO","LNELAAV","LNELALN",
               "TrapOtherTypeCount","TrapOtherTypeSpecifications","TrapOtherTypeLength","TrapOtherTypeWidth",
               "TrapOtherTypeHeight")

TDGDLF.survey=TDGDLF.survey[,-match(irrelevant.vars,names(TDGDLF.survey))]

#Select records for identified fishers
TDGDLF.survey.not.selected=subset(TDGDLF.survey,!SKIPPER%in%FISHERS)
TDGDLF.survey.selected=subset(TDGDLF.survey,SKIPPER%in%FISHERS)

Table.year.SKIPPER=table(TDGDLF.survey$SKIPPER,TDGDLF.survey$year)
Table.year.vessel=table(TDGDLF.survey$BOATREGO,TDGDLF.survey$year)

