#-------SCRIPT FOR ANALYSIS VESSEL CHARACTERISTIC SURVEY FOR TDGDLF ---#

library('stringr')
library('lubridate')
library(readxl)
library(tidyverse)


# --------------------Import data-----------------------------
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Data/Fishing power"))

TDGDLF.vessels=read.csv("TDGDLF.vessels.csv",stringsAsFactors=F)
List_boats_skippers=read.csv(handl_OneDrive("Data/Catch and Effort/List_boats_skippers.csv"),stringsAsFactors=F)
TDGDLF.skippers=read.csv("TDGDLF.SKIPERS.csv",stringsAsFactors=F)

#Fishing licence registration survey
Vessel.survey=read.csv("VesselGearSurveyData/VesselGearSurveyData.csv",stringsAsFactors=F)
Vessel.survey_meta=read_excel("VesselGearSurveyData/VesselGearSurveyData.xlsx", sheet = "MetaData")

#use questionnaire to validate survey answers
Matias.questionnaire=read_excel("Questionnaire responses/Questionaire.xlsx", sheet = "Questionaire")
Steve.T.questionnaire=read_excel("Effort Creep.Inteviews_Steve_Taylor/Fishing Power extract for Matias.xlsx", sheet = "Fishing Power") 

#Mackerel vessels
Mackerel.vessels=read_excel("Mackerel vessels for history.xlsx", sheet = "Sheet1")


# --------------------Manipulate questionnaire data-----------------------------
colnames(Steve.T.questionnaire)=gsub(" ", ".", colnames(Steve.T.questionnaire))
Steve.T.questionnaire=Steve.T.questionnaire%>%
                        mutate(Last.name=gsub("^\\S+ ", "", Fisher.name))

# --------------------Manipulate survey data-----------------------------
# Select TDGDLF only vessels
TDGDLF.vessels=unique(c(TDGDLF.vessels$VESSEL,List_boats_skippers$vessel))
TDGDLF.skippers=unique(c(TDGDLF.skippers$SKIPPERS,List_boats_skippers$MastersName))

TDGDLF.Vessel.survey_meta=Vessel.survey_meta%>%
                            filter(Fishery%in%c('All','Net'))%>%
                      arrange(Type,Type2,Variable)
relevant.variables=TDGDLF.Vessel.survey_meta$Variable
key.vars=c("BOATREGO","LICYEAR","BOATNAME","SKIPPER","PFL","DATE")

Mackerel.survey=Vessel.survey%>%
                  filter(BOATREGO%in%Mackerel.vessels$LFB1 | 
                           tolower(Vessel.survey$BOATNAME)%in%tolower(Mackerel.vessels$Name))

TDGDLF.survey=Vessel.survey%>%
                filter(BOATREGO%in%TDGDLF.vessels)%>%
                filter(!DATE=="")%>%
                mutate(across(where(is.character), toupper))%>%
                mutate(DATE=as.Date(DATE,"%d/%m/%Y"),
                       year=year(DATE),
                       year=ifelse(year==1903,1993,year),
                       month=month(DATE),
                       BOATNAME=case_when(BOATREGO=="E 059"~"FALCON 2",
                                          BOATNAME=="ALBA-MARINA II"~"ALBA MARINA II",
                                          BOATNAME=="AVESPUCCI II"~"A VESPUCCI II",
                                          BOATNAME=="KARIN-E"~"KARINE-E",
                                          BOATNAME=="CO SAN II"~"COSAN II",
                                          BOATNAME=="PLANJAK II"~"PLANJAK",
                                          BOATNAME=="SEA VIEW II"~"SEAVIEW II",
                                          BOATNAME%in%c("ST. GERARD M","ST GERARD M")~"ST. GERARD",
                                          BOATNAME=="SVETI-NIKOLA"~"SVET-NIKOLA",
                                          BOATNAME=="WAVE DANCE I"~"WAVE DANCER I",
                                          TRUE~BOATNAME),
                       DATEBUILT=as.integer(DATEBUILT),
                       DATEBUILT=case_when(DATEBUILT==14 & BOATNAME=="MOBYDICK II"~1994,
                                           TRUE~DATEBUILT),
                       REG_LENGTH=as.numeric(REG_LENGTH),
                       MAXBEAM=as.numeric(MAXBEAM),
                       MAXDRAU=as.numeric(gsub('[A-Za]','',MAXDRAU)),
                       GROSSTON=as.numeric(gsub('[A-Za]','',GROSSTON)),
                       ENGMAKE=case_when(ENGMAKE=="VOLVO - PENTA"~"VOLVO PENTA",
                                         ENGMAKE=="CATT"~"CAT",
                                         ENGMAKE=="DETRIOT"~"DETROIT",
                                         ENGMAKE=="FORD LEES (X 2)"~"FORD LEES (X2)",
                                         ENGMAKE%in%c("GARDENER","GARDNER")~"GARDINER",
                                         ENGMAKE%in%c("MAN","M A N")~"M.A.N.",
                                         ENGMAKE=="M T U"~"M.T.U.",
                                         ENGMAKE=="MARINARE"~"MARINER",
                                         ENGMAKE%in%c("MERCRUSER","MERCRUISERS")~"MERCRUISER",
                                         ENGMAKE%in%c("SCANIS","SCARIA")~"SCANIA",
                                         ENGMAKE%in%c("YAMAH","YAMMAHA")~"YAMAHA",
                                         TRUE~ENGMAKE),
                       HULLCONS=case_when(HULLCONS%in%c('','1')~NA_character_,
                                          TRUE~HULLCONS),
                       HULLNUMB=case_when(HULLNUMB%in%c('')~NA_character_,
                                          TRUE~HULLNUMB),
                       HULLTYPE=case_when(HULLTYPE%in%c('')~NA_character_,
                                          TRUE~HULLTYPE),
                       WHEELHOU=case_when(WHEELHOU%in%c('')~NA_character_,
                                          TRUE~WHEELHOU),
                       FLYBRIDGE=case_when(FLYBRIDGE%in%c('')~NA_character_,
                                          TRUE~FLYBRIDGE),
                       ENGMAKE=case_when(ENGMAKE%in%c('')~NA_character_,
                                           TRUE~ENGMAKE),
                       ENGMODL=case_when(ENGMODL%in%c('')~NA_character_,
                                         TRUE~ENGMODL),
                       ENGNUM=case_when(ENGNUM%in%c('')~NA_character_,
                                         TRUE~ENGNUM),
                       ENGPOWR=case_when(ENGPOWR%in%c('')~NA_character_,
                                         ENGPOWR%in%c('336X2')~"672",
                                        TRUE~ENGPOWR),
                       ENGPOWR=as.numeric(ENGPOWR),
                       ENGDERAT=case_when(ENGDERAT%in%c('')~NA_character_,
                                         TRUE~ENGDERAT),
                       ENGSPD=as.numeric(ENGSPD),
                       ENGSPD=case_when(ENGSPD%in%c('')~NA_real_,
                                         TRUE~ENGSPD),
                       LHTABV=as.numeric(LHTABV),
                       LHTABV=case_when(LHTABV%in%c('')~NA_real_,
                                        TRUE~LHTABV),
                       LHTBLW=as.numeric(LHTBLW),
                       LHTBLW=case_when(LHTABV%in%c('')~NA_real_,
                                        TRUE~LHTBLW),
                       FRZCAP=as.numeric(FRZCAP),
                       FRZCAP=case_when(FRZCAP%in%c('')~NA_real_,
                                        TRUE~FRZCAP),
                       BRNTNK=as.numeric(BRNTNK),
                       BRNTNK=case_when(BRNTNK%in%c('')~NA_real_,
                                        BRNTNK%in%c(15000)~15,
                                        TRUE~BRNTNK),
                       ICEBOX=case_when(ICEBOX%in%c('')~NA_character_,
                                        ICEBOX%in%c('500')~'5',
                                        ICEBOX%in%c('3 X 2')~'6',
                                        TRUE~ICEBOX),
                       BWECHOM=case_when(BWECHOM%in%c('')~NA_character_,
                                         TRUE~BWECHOM),
                       BWECHOT=case_when(BWECHOT%in%c('')~NA_character_,
                                         TRUE~BWECHOT),
                       COECHOM=case_when(COECHOM%in%c('')~NA_character_,
                                         COECHOM%in%c("FURNUNO","FURUNO X 2","FURUNO  / FURUNO")~"FURUNO",
                                         COECHOM%in%c("JRC + FURUNO","FURUNO  / JRC")~"JRC / FURUNO",
                                         COECHOM%in%c("JRC 120 DUAL FREQ.")~"JRC",
                                         TRUE~COECHOM),
                       COECHOT=case_when(COECHOT%in%c('')~NA_character_,
                                         COECHOT%in%c('130 + 1100L')~"130  / 1100L",
                                         TRUE~COECHOT),
                       RADARM=case_when(RADARM%in%c('')~NA_character_,
                                         TRUE~RADARM),
                       RADART=case_when(RADART%in%c('')~NA_character_,
                                        TRUE~RADART),
                       SONARM=case_when(SONARM%in%c('')~NA_character_,
                                        TRUE~SONARM),
                       SONART=case_when(SONART%in%c('')~NA_character_,
                                        TRUE~SONART),
                       SNAVTRM=case_when(SNAVTRM%in%c('')~NA_character_,
                                        TRUE~SNAVTRM),
                       SNAVTRT=case_when(SNAVTRT%in%c('')~NA_character_,
                                         TRUE~SNAVTRT),
                       SNAVGPM=case_when(SNAVGPM%in%c('')~NA_character_,
                                         TRUE~SNAVGPM),
                       SNAVGPT=case_when(SNAVGPT%in%c('')~NA_character_,
                                         TRUE~SNAVGPT),
                       SKIPYEAR=as.numeric(SKIPYEAR),
                       SKIPYEAR=case_when(SKIPYEAR%in%c('')~NA_real_,
                                        TRUE~SKIPYEAR),
                       DECKYEAR=as.numeric(DECKYEAR),
                       DECKYEAR=case_when(DECKYEAR%in%c('')~NA_real_,
                                          TRUE~DECKYEAR),
                       PLOTM=case_when(PLOTM%in%c('')~NA_character_,
                                       PLOTM%in%c("C--PLOT","C PLOT")~"C-PLOT",
                                         TRUE~PLOTM),
                       PLOTT=case_when(PLOTT%in%c('')~NA_character_,
                                       PLOTT%in%c("C--PLOT","C PLOT")~"C-PLOT",
                                       TRUE~PLOTT),
                       LHTABV_T=case_when(LHTABV_T%in%c('')~NA_character_,
                                         TRUE~LHTABV_T),
                       LHTBLW_T=case_when(LHTBLW_T%in%c('')~NA_character_,
                                          LHTBLW_T%in%c('.200')~'0.2',
                                          LHTBLW_T%in%c('800')~'8',
                                          TRUE~LHTBLW_T),
                       LHTBLW_T=as.numeric(LHTBLW_T),
                       FRZCAP_T=as.numeric(FRZCAP_T),
                       FRZCAP_T=case_when(FRZCAP_T%in%c('')~NA_real_,
                                          TRUE~FRZCAP_T),
                       BRNTNK_T=as.numeric(BRNTNK_T),
                       BRNTNK_T=case_when(BRNTNK_T%in%c('')~NA_real_,
                                          TRUE~BRNTNK_T),
                       ICEBOX_T=case_when(ICEBOX_T%in%c('')~NA_character_,
                                         TRUE~ICEBOX_T),
                       PLOTYR=case_when(PLOTYR%in%c('')~NA_character_,
                                         TRUE~PLOTYR),
                       GPSYR=case_when(GPSYR%in%c('')~NA_character_,
                                       GPSYR%in%c('5600')~"2005",
                                        TRUE~GPSYR),
                       COECHOYR=case_when(COECHOYR%in%c('','JRC')~NA_character_,
                                       TRUE~COECHOYR),
                       BOATNAME.BOATREGO=paste(BOATNAME,BOATREGO,sep='-'),
                       BRNTNK=ifelse(is.na(BRNTNK),BRNTNK_T,BRNTNK),
                       FRZCAP=ifelse(is.na(FRZCAP),FRZCAP_T,FRZCAP),
                       ICEBOX=ifelse(is.na(ICEBOX),ICEBOX_T,ICEBOX),
                       ICEBOX=as.numeric(ICEBOX),
                       LHTABV=ifelse(is.na(LHTABV),LHTABV_T,LHTABV),
                       LHTABV=as.numeric(LHTABV),
                       LHTBLW=ifelse(is.na(LHTBLW),LHTBLW_T,LHTBLW),
                       LHTABV=ifelse(LHTABV<1,NA,LHTABV),
                       LHTBLW=ifelse(LHTBLW<1,NA,LHTBLW))

relevant.variables=c(relevant.variables,"year","month","BOATNAME.BOATREGO")
relevant.variables=subset(relevant.variables,!relevant.variables%in%c("COMMENT","VesselGearSurveyDataId",
                                                                      "VIDEOM","VIDEOT",
                                                                      "NETDINGY","NETJETB","NETFISH",
                                                                      "NETHAND","ROXM","ROXT",
                                                                      "DGPSM","DGPST","NOBOAT",
                                                                      "ROXYR","DGPSYR",
                                                                      "ENGFUEL","ENGKORT",
                                                                      "BRNTNK_T","FRZCAP_T",
                                                                      "ICEBOX_T","LHTABV_T",
                                                                      "LHTBLW_T",paste0("NETTP",1:14)))
TDGDLF.survey=TDGDLF.survey%>%
                dplyr::select(all_of(relevant.variables))

vessel.vars=colnames(TDGDLF.survey)
vessel.vars=vessel.vars[-match(c("BOATREGO","LICYEAR","BOATNAME","SKIPPER","PFL","DATE",
                                 "year","month","BOATNAME.BOATREGO"),vessel.vars)]
names(vessel.vars)=TDGDLF.Vessel.survey_meta%>%filter(Variable%in%vessel.vars)%>%pull(Description)
  
Table.year=table(TDGDLF.survey$year)
Table.vessels=table(TDGDLF.survey$BOATREGO)
Table.skippers=table(TDGDLF.survey$SKIPPER)
Table.year.SKIPPER=table(TDGDLF.survey$SKIPPER,TDGDLF.survey$year)
Table.year.vessel=table(TDGDLF.survey$BOATREGO,TDGDLF.survey$year)

#Top TDGDLF fishers (90% of catch)
Zn2.fishers=c("MANSTEAD","MANSTED","TONKIN","ISARAYEPPOS","HENDERSON",
              "SOUMELIDIS","SHARP","GILBERT","KENNEDY","PETERS","BROWN",
              "BRADLEY","TINDALL","WHETSONE","WHETSTONE","STEEL","STEELE","BLACK",
              "BRANDERHORST","TOUMAZOS","TRIANTAFYLLOU","GOODALL","SANDERS","HAMILTON",
              "DIMOFF","KARATERPOS","MADGEN")
Zn1.fishers=c("COOKE","SCIMONE","SOULOS","MILES","COCKMAN","ROBB","WARRILOW",
              "REAY","SOULIS","ADAMS","SELL","THOMAS","THORNTON","ENGLISH")
WC.fishers=c("PARKER","AITCHISON","CARR","MCWHIRTER","BUCKERIDGE")
Fisher.lookup=c(Zn2.fishers,Zn1.fishers,WC.fishers)

include <- function (theList, toMatch) unique (grep(paste(toMatch,collapse="|"),theList, value=TRUE))
FISHERS=include(theList=names(Table.skippers),toMatch=Fisher.lookup)

#Top Vessel regos (90% of catch)
VESSELs=c("B 067","B 142","E 009","E 004","E 034","E 067","G 297","E 056","E 045","F 517","E 035","B 022",
          "F 520","E 007","A 201","E 010","B 038","F 768","E 055","B 048","F 417","F 244","B 091","E 061",
          "F 505","F 756","B 012","E 030","B 098","E 059")

#Select records for identified fishers
TDGDLF.survey.selected=TDGDLF.survey%>%
                        filter(SKIPPER%in%FISHERS)%>%
                        filter(BOATREGO%in%VESSELs)




# Visualize records-------------------------------------------------------------------------
setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Efficiency creep/All_Vessel_characteristics"))

unik.BOATNAME.BOATREGO=unique(TDGDLF.survey.selected$BOATNAME.BOATREGO)
unik.BOATNAME=unique(TDGDLF.survey.selected$BOATNAME)
out.this=FALSE
if(out.this)
{
  fn.viz=function(dis,sortby)
  {
    d=TDGDLF.survey.selected%>%
      mutate(BOATREGO.BOATNAME=paste(BOATREGO,BOATNAME,sep='.'))
    
    d=d[,c(key.vars,sortby,dis)]%>%
      mutate(x=!!sym(dis))
    if(is.numeric(d$x)|is.integer(d$x))
    {
      LvLs=sort(unique(d$x))
      d$x=factor(d$x,levels=LvLs)
    }
    NKl=1
    if(length(unique(d$x))>20) NKl=2
    d$var=d[,sortby]
    p=d%>%
      ggplot(aes(x=DATE,y=var,color=x))+
      geom_line(linewidth=5)+
      theme(panel.background = element_rect(fill = "cadetblue3",
                                            colour = "black",
                                            size = 0.5, linetype = "solid"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.title=element_blank(),
            legend.text = element_text(size=6))+
      ggtitle(names(dis))+scale_color_grey(start=0,end=1)+guides(color=guide_legend(ncol =NKl))
    print(p)
  }
  pdf("All vessel variables_sorted_boatname.pdf")
  for(i in 1:length(vessel.vars))
  {
    fn.viz(dis=vessel.vars[i],sortby='BOATNAME.BOATREGO')
    print(vessel.vars[i])
  }
  dev.off()
  
  pdf("All vessel variables_sorted_boatrego.pdf")
  for(i in 1:length(vessel.vars))
  {
    fn.viz(dis=vessel.vars[i],sortby='BOATREGO.BOATNAME')
    print(vessel.vars[i])
  }
  dev.off()
}


# Validate records-------------------------------------------------------------------------
a=TDGDLF.survey.selected%>%
  filter(grepl(paste(toupper(Steve.T.questionnaire$Last.name),collapse='|'),SKIPPER))%>%
  distinct(SKIPPER,SONARM,GPSYR,COECHOYR,PLOTYR)
b=Steve.T.questionnaire%>%distinct(Fisher.name,Sonar,GPS,Colour.echo.sounder,Plotter,
                                   Efficiency.1,Efficiency.2,Efficiency.3)


#No boat in Matias' questionnaire found in TDGDLF.survey
# a=TDGDLF.survey.selected%>%filter(grepl("Warrilow",tolower(SKIPPER)))
# b=Matias.questionnaire%>%filter(Surname=="Warrilow")
# a$ENGPOWR
# b$'EnginePower(hp)'


# Extract useful vessel data for cpue standard.-------------------------------------------
useful.vessel.var=vessel.vars[-c(grep(paste(c('NETDP','NETLN','NETMSL','NETMSS'),collapse = '|'),vessel.vars),
                                  match(c("SNAVTRM","SNAVTRT"),vessel.vars))]
ID.irrelevant.variables=c(grep(paste(c('NETDP','NETLN','NETMSL','NETMSS'),collapse = '|'),relevant.variables),
                          match(c("SNAVTRM","SNAVTRT"),relevant.variables))

#Get year of when technology was incorporated
#Note:  If technology incorporated before first year of responding survey and there is no specific question asking
# the 'year' of purchase of the technology (BWECHO, SONAR, RADAR), then year of incorporation will be biased high
TDGDLF.cpue.stand=TDGDLF.survey.selected%>%
                  dplyr::select(all_of(relevant.variables[-ID.irrelevant.variables]))%>%
                  mutate(GPSYR=as.numeric(GPSYR),
                         GPSYR=case_when(is.na(GPSYR) & !is.na(SNAVGPM)~year,
                                         is.na(GPSYR) & is.na(SNAVGPM) & !is.na(SNAVGPT)~year,
                                         TRUE~GPSYR),
                         GPSYR=ifelse(is.na(GPSYR),0,GPSYR),
                         GPS=GPSYR,
                         PLOTYR=as.numeric(PLOTYR),
                         PLOTYR=case_when(is.na(PLOTYR) & !is.na(PLOTM)~year,
                                          is.na(PLOTYR) & is.na(PLOTM) & !is.na(PLOTT)~year,
                                          TRUE~PLOTYR),
                         PLOTYR=ifelse(is.na(PLOTYR),0,PLOTYR),
                         PLOT=PLOTYR,
                         COECHOYR=as.numeric(COECHOYR),
                         COECHOYR=case_when(is.na(COECHOYR) & !is.na(COECHOM)~year,
                                            is.na(COECHOYR) & is.na(COECHOM) & !is.na(COECHOT)~year,
                                            TRUE~COECHOYR),
                         COECHOYR=ifelse(is.na(COECHOYR),0,COECHOYR),
                         COECHO=COECHOYR,
                         RADARYR=case_when(!is.na(RADARM)~year,
                                           is.na(RADARM) & !is.na(RADART)~year,
                                           TRUE~NA_real_),
                         RADARYR=ifelse(is.na(RADARYR),0,RADARYR),
                         RADAR=RADARYR,
                         SONARYR=case_when(!is.na(SONARM)~year,
                                           is.na(SONARM) & !is.na(SONART)~year,
                                           TRUE~NA_real_),
                         SONARYR=ifelse(is.na(SONARYR),0,SONARYR),
                         SONAR=SONARYR,
                         BWECHOYR=case_when(!is.na(BWECHOM)~year,
                                            is.na(BWECHOM) & !is.na(BWECHOT)~year,
                                            TRUE~NA_real_),
                         BWECHOYR=ifelse(is.na(BWECHOYR),0,BWECHOYR),
                         BWECHO=BWECHOYR)%>%
  dplyr::select(-c(GPSYR,SNAVGPM,SNAVGPT,
                   PLOTYR,PLOTM,PLOTT,
                   COECHOYR,COECHOM,COECHOT,
                   RADARYR,RADARM,RADART,
                   SONARYR,SONARM,SONART,
                   BWECHOYR,BWECHOM,BWECHOT,
                   year,DATE,month,BOATNAME.BOATREGO,
                   ENGMAKE,ENGMODL,
                   DECKYEAR, SKIPYEAR,  #meaningless, not updated each year of survey year and NA for some key players
                   NETHWINC,NETPOWER,NETVEHI,NETVWINC))%>% #hauling vars are incomplete, e.g. no data for top vessel ('SOUTH WESTERN')
  group_by(BOATNAME,BOATREGO)%>%
  mutate(GPS=min(GPS),
         PLOT=min(PLOT),
         COECHO=min(COECHO),
         RADAR=min(RADAR),
         SONAR=min(SONAR),
         BWECHO=min(BWECHO))%>%
  ungroup()
#Fix some typos
TDGDLF.cpue.stand=TDGDLF.cpue.stand%>%
  mutate(SONAR=ifelse(SONAR==0,NA,SONAR),
         PLOT=ifelse(PLOT==0,NA,PLOT),
         GPS=ifelse(GPS==0,NA,GPS),
         RADAR=ifelse(RADAR==0,NA,RADAR),
         BWECHO=ifelse(BWECHO==0,NA,BWECHO),
         COECHO=ifelse(COECHO==0,NA,COECHO),
         ENGDERAT=case_when(BOATNAME=='DESTINY' & is.na(ENGDERAT)~'N',
                            BOATNAME=='SANTA BARBARA II' & is.na(ENGDERAT)~'N',
                            TRUE~ENGDERAT),
         ENGNUM=case_when(BOATNAME=='STORMRAKER' & is.na(ENGNUM)~'S',
                          BOATNAME=='SANTA BARBARA II' & is.na(ENGNUM)~'S',
                          TRUE~ENGNUM),
         ENGPOWR=case_when(BOATNAME=='FATAL ATTRACTION' & ENGPOWR==574~380,
                           BOATNAME=='STORMRAKER' & is.na(ENGPOWR)~270,
                           BOATNAME=='SANTA BARBARA II' & is.na(ENGPOWR)~283,
                             TRUE~ENGPOWR),
         ENGSPD=case_when(BOATNAME=='STORMRAKER' & is.na(ENGSPD)~13,
                          BOATNAME=='SANTA BARBARA II' & is.na(ENGSPD)~16,
                           TRUE~ENGSPD),
         BRNTNK=case_when(BOATNAME=='OCEAN MISTRESS' & BRNTNK==6~8,
                          BOATNAME=='STORMRAKER' & is.na(BRNTNK)~20,
                          TRUE~BRNTNK),
         ICEBOX=case_when(BOATNAME=='SANTA BARBARA II' & is.na(ICEBOX)~2,
                          TRUE~ICEBOX),
         DATEBUILT=case_when(BOATNAME=='STORMRAKER' & is.na(DATEBUILT)~1988,
                             BOATNAME=='SANTA BARBARA II' & is.na(DATEBUILT)~1979,
                             BOATNAME=='ELIZABETH MARIA II' & is.na(DATEBUILT)~1978,
                             BOATNAME=='FALCON 2' ~1978,
                             BOATNAME=='WARNBRO LADY' ~1974,
                             BOATNAME=='SAN MARGO' ~1977,
                          TRUE~DATEBUILT),
         FLYBRIDGE=case_when(BOATNAME=='STORMRAKER' & is.na(FLYBRIDGE)~'Y',
                             BOATNAME=='SANTA BARBARA II' & is.na(FLYBRIDGE)~'N',
                             TRUE~FLYBRIDGE),
         GROSSTON=case_when(BOATNAME=='SANTA BARBARA II' & is.na(GROSSTON)~11,
                            BOATNAME=='DOREEN' & is.na(GROSSTON)~80,
                          TRUE~GROSSTON),
         HULLCONS=case_when(BOATNAME=='DESTINY' & is.na(HULLCONS)~'F',
                            BOATNAME=='STORMRAKER' & is.na(HULLCONS)~'A',
                            BOATNAME=='SANTA BARBARA II' & is.na(HULLCONS)~'P',
                            TRUE~HULLCONS),
         HULLNUMB=case_when(BOATNAME=='STORMRAKER' & is.na(HULLNUMB)~'S',
                            BOATNAME=='SANTA BARBARA II' & is.na(HULLNUMB)~'S',
                            TRUE~HULLNUMB),
         HULLTYPE=case_when(BOATNAME=='DESTINY' & is.na(HULLTYPE)~'H',
                            BOATNAME=='STORMRAKER' & is.na(HULLTYPE)~'H',
                            BOATNAME=='FALCON 2' & is.na(HULLTYPE)~'H',
                            TRUE~HULLTYPE),
         MAXBEAM=case_when(BOATNAME=='STORMRAKER' & is.na(MAXBEAM)~5,
                           BOATNAME=='SANTA BARBARA II' & is.na(MAXBEAM)~4,
                           TRUE~MAXBEAM),
         MAXDRAU=case_when(BOATNAME=='STORMRAKER' & is.na(MAXDRAU)~1.9,
                           BOATNAME=='SANTA BARBARA II' & is.na(MAXDRAU)~1.3,
                           TRUE~MAXDRAU),
         REG_LENGTH=case_when(BOATNAME=='STORMRAKER' & is.na(REG_LENGTH)~17.6,
                              BOATNAME=='SANTA BARBARA II' & is.na(REG_LENGTH)~11.6,
                              BOATNAME=='DOREEN' & is.na(REG_LENGTH)~18.3,
                              TRUE~REG_LENGTH),
         WHEELHOU=case_when(BOATNAME=='STORMRAKER' & is.na(WHEELHOU)~'F',
                            BOATNAME=='SANTA BARBARA II' & is.na(WHEELHOU)~'F',
                            TRUE~WHEELHOU))

#Export
setwd(handl_OneDrive("Data/Fishing power"))
write.csv(TDGDLF.cpue.stand,'Vessel.charac.for_TDGDLF.cpue.stand.csv',row.names = F)
write.csv(Mackerel.survey,'Vessel.charac.Mackerel.csv',row.names = F)
