#------ G 3. STATE OF FISHERIES REPORT ------

library(dplyr)
library(ggplot2)
library(ggrepel)
library(zoo)
library(lubridate)
library(gganimate)
library(gapminder)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(stringr)

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/SoFaR.figs.R"))

# Must run:
#       catch rate standardisation script (2.CPUE standardisations.R) before so it can bring in the index with the latest year
#      'recons_recreational.csv' to bring in current year of reconstructed rec catch

#note:  TDGDLF is defined as METHOD= "GN" or "LL" and non-estuaries and south of 26 S. 
#         'Other fisheries' are extracted from 'Other.fishery.catch' (this includes 'Data.monthly'  
#         for which gear is not 'GN' or 'LL' and 'Daily logbooks' from other fisheries that report 
#         sharks (Pilbara trawl, WCDF, mackerel, Kimberley barramundi, Shark Bay prawn, etc) 
#       table 81.d for monthly records has no teleost catches 
#     SOFAR effort is reported as km.gn.d. calculated with Use.Date to "YES" (as instructed by Rory). 
#     If Use.Date = "NO" then 'ID' is used in the aggregation rather than 'date', then effort
#       in recent years is larger due to the multiple shot per day practice of
#       vessels from zone 2. Rory instructed this is biased. See '4.Script for testing effect of DATE or ID 
#       in effort calculation'.  This is not an issue for km.gn.hr

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240)

Current.yr="2023-24"   #financial year for latest catch

Percent.fin.of.livewt=0.03  #used as default if species not in Shark.Fin.Price.List
Shark.Fin.Price.List=read.csv(handl_OneDrive("Data/Catch and Effort/Shark Fin Price List.csv")) #provided by Rory to Eva Lai


TDGDLF.lat.range=c(-26,-40)

#bring in data from 1.Manipulate data.R
setwd(handl_OneDrive("Analyses/Data_outs"))

Data.monthly=read.csv("Data.monthly.csv")
Data.daily=read.csv("Data.daily.csv")
Rec.fish.catch=read.csv('recons_recreational.csv')   
Data.current.Sofar=read.csv("Data.current.Sofar.csv")
PRICES=read.csv("PRICES.csv")
Total.effort.days.monthly=read.csv("Annual.total.eff.days.csv")
Total.effort.hours.monthly=read.csv("Annual.total.eff.hours.csv")
Total.effort.zone.hours.monthly=read.csv("Annual.zone.eff.hours.csv")
Total.effort.zone.days.monthly=read.csv("Annual.zone.eff.days.csv")
TEPS.current=read.csv("TEPS.current.csv")
Suite=read.csv("suite.csv")$Suite
Results.pre.2013=read.csv("Results.pre.2013.csv")
Effort.monthly=read.csv("Effort.monthly.csv")   
Effort.daily=read.csv("Effort.daily.csv")
LL.equiv.Eff.days.zone=read.csv("LL.equiv.Eff.days.zone.csv")


All.species.names=read.csv(handl_OneDrive("Data/Species.code.csv"))
Shark.species=5001:24999
Ray.species=c(25000:31000,38000,39001,90030)
Elasmo.species=c(Shark.species,Ray.species)
Gummy=17001;Dusky_whaler=c(18003,18001);Whiskery=17003;Sandbar=18007;Hammerheads=19000;
Spinner=18023;Wobbegongs=13000;Common_saw_shark=23002;School=17008
Order.Elas.Sp.SoFAR=data.frame(Species=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads",
                                         "Spinner","Wobbegongs","Rays","Common_saw_shark","School","Other_elasmobranchs",
                                         "Total Elasmobranchs"),
                               Species_or_taxon=c("Mustelus antarcticus","Carcharhinus obscurus","Furgaleus macki",
                                                  "Carcharhinus plumbeus","F.Sphyrnidae","Carcharhinus brevipinna",
                                                  "F. Orectolobidae","Batoidea","Pristiophorus cirratus","Galeorhinus galeus","",""))


Scalefish.species=c(117001,180000:599001)
Scalefish.species=Scalefish.species[which(Scalefish.species %in% unique(Data.monthly$SPECIES))]
Redfishes=c(258000,258004,258005,258006)
Blue_morwong=377004;Blue_groper=384002;West_Australian_dhufish=320000;Pink_snapper=353001;
Boarfishes=367000;Samsonfish=337007;Redfishes=Redfishes;
Mulloway=354001;Sweetlips=350000;Baldchin_groper=384999
Order.Scale.Sp.SoFAR=data.frame(Species=c("Blue_morwong","Blue_groper","Pink_snapper","West_Australian_dhufish",
                                          "Boarfishes","Samsonfish","Redfishes","Mulloway","Sweetlips",
                                          "Baldchin_groper","Other_scalefish","Total Scalefish"),
                                Species_or_taxon=c("Nemadactylus valenciennesi","Achoerodus gouldii",
                                                   "Chrysophrys auratus","Glaucosoma hebraicum","F. Pentacerotidae","Seriola hippos",
                                                   "Centroberyx spp.","Argyrosomus japonicus","F. Haemulidae",
                                                   "Choerodon rubescens","",""))

Mollusc.species=600000:610000
Mollusc.species=Mollusc.species[which(Mollusc.species %in% unique(Data.monthly$SPECIES))]



#Fishing effort limits (2001-02)
FishEffLims=data.frame(zone=c("West","Zone1","Zone2"),Km.Gillnet.Hours.c=c(67692,84075,144102),
                       Km.Gillnet.Days.c=c(2832,3503,7205))


#species status for SOFAR
Gummy.status="Adequate"
Dusky.status="Recovering"
Sandbar.status="Adequate"
Whiskery.status="Adequate"


#Catch range of indicator species (Gummy, whiskery, dusky, sandbar)
Catch.range.key.species=data.frame(SPECIES=c("Gummy","Whiskery","Bronzy.Dusky","Sandbar"),
                                   Min.catch=c(350,175,200,0),
                                   Max.catch=c(450,225,300,120))%>%
                        mutate(SPECIES=as.character(SPECIES))

handle.Sofar=paste(handl_OneDrive("Analyses/Catch and effort/State of fisheries/"),Current.yr,sep="")
if(!file.exists(handle.Sofar))dir.create(handle.Sofar)

Display.current.yr=paste(substr(Current.yr,1,4),"/",substr(Current.yr,6,7),sep="")
Ind.spe.list=list(Gummy=17001,Whiskery=17003,Bronzy.Dusky=c(18001,18003),sandbar=18007)


#Create current year data set for TDGDLF
DAT=subset(Data.monthly,FINYEAR==Current.yr & METHOD%in%c("GN","LL") & Estuary=="NO" &
             LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])

#add small catch to WC (Gascoyne appears because of upper boundary set to 26S)
DAT$Bioregion=with(DAT,ifelse(Bioregion=="Gascoyne","WC",Bioregion))  


#Create current year data set for other fisheries in temperate region (i.e non TDGDLF or NSF)  
#select south and west coast fisheries
Other.fishery.catch=subset(Data.monthly,SPECIES<=31000 & Shark.fishery=='non.shark.fishery')  
Other.fishery.catch=subset(Other.fishery.catch, !(METHOD=="LL"&Estuary=="NO"))
Other.fishery.catch$LAT=-as.numeric(substr(Other.fishery.catch$BLOCKX,1,2))
Other.fishery.catch=subset(Other.fishery.catch, !(METHOD=="GN"&Estuary=="NO" & LAT<(-26)))

Other=subset(Other.fishery.catch,FINYEAR==Current.yr & LAT<(-26))  
#  Other=subset(Other.fishery.catch,FINYEAR==Current.yr & fishery%in%c("WCRL","WCE"," SWT","SCE","OASC","CSC","CSLP"))  
names(Other)[match('FINYEAR',names(Other))]="financial.year"

#--- Main Feature table (total catch, indicator species catch, teleost catch, etc) ---#

#Shark species indicators
Tot.wt=Ind.spe.list
for(i in 1:length(Tot.wt))
{
  Data=subset(DAT,SPECIES%in%Ind.spe.list[[i]])
  Tot.wt[[i]]=round(sum(Data$LIVEWT.c,na.rm=T)/1000)                                   
}
Tot.wt=as.data.frame(unlist(Tot.wt))
Tot.wt$SPECIES=rownames(Tot.wt)
rownames(Tot.wt)=NULL
colnames(Tot.wt)[1]="Catch.tons"
Tot.wt=rbind(Tot.wt,data.frame(Catch.tons=round(sum(Tot.wt[,1])),SPECIES="Total indicators"))

#Total shark catch other fisheries
Tot.wt.other=data.frame(Catch.tons=round(sum(as.numeric(Other$LIVEWT))/1000,1),
                        SPECIES="Sharks and rays by other commercial fisheries") 

Tot.wt=rbind(Tot.wt,Tot.wt.other)     

#Total catch of elasmos and teleosts  
a=subset(DAT,SPECIES%in%Elasmo.species)
Tot.elasmos=round(sum(a$LIVEWT.c)/1000,0)
Tot.elasmos=data.frame(Tot.elasmos)
names(Tot.elasmos)="Catch.tons"
Tot.elasmos=data.frame(Tot.elasmos,SPECIES="Total sharks and rays")

a=subset(DAT,SPECIES%in%Scalefish.species)
Tot.teleosts=round(sum(a$LIVEWT.c)/1000,0)
Tot.teleosts=data.frame(Tot.teleosts)
names(Tot.teleosts)="Catch.tons"
Tot.teleosts=data.frame(Tot.teleosts,SPECIES="Total scalefish")

MainFeatures=rbind(Tot.elasmos,Tot.teleosts,Tot.wt)
MainFeatures=MainFeatures[,match(c("SPECIES","Catch.tons"),names(MainFeatures))]


#Add latest data on Recreational fishing catch
Rec.fi.Tot.catch=Rec.fish.catch%>%
    group_by(FINYEAR)%>%
    summarise(Catch.tons=sum(LIVEWT.c)/1000)%>%
    mutate(yr=as.numeric(substr(FINYEAR,1,4)),
           yr.min=abs(yr-as.numeric(substr(Current.yr,1,4))))%>%
    filter(yr.min==min(yr.min))%>%
    mutate(SPECIES=paste("Reconstructed statewide recreational catch ","(",FINYEAR,")",sep=""),
           Percnt.Tot.Com= Catch.tons*100/MainFeatures$Catch.tons[which(MainFeatures$SPECIES=="Total sharks and rays")],
           Catch.tons=paste0(ceiling(Percnt.Tot.Com),'% of commercial catch'))%>%
    data.frame
REC.yr=Rec.fi.Tot.catch$FINYEAR
Rec.fi.Tot.catch=Rec.fi.Tot.catch%>%
                  dplyr::select(SPECIES,Catch.tons)

Ktch.indic=subset(MainFeatures,SPECIES=="Total indicators",select=Catch.tons)$Catch.tons

MainFeatures$Catch.tons=paste(round(MainFeatures$Catch.tons),"t")
MainFeatures=rbind(MainFeatures,Rec.fi.Tot.catch)


#number of licenses operating in current year
dummy=subset(Data.current.Sofar,finyear==Current.yr,select=c(finyear,month,blockx,vessel,crew,MastersName))
dummy$blockx=with(dummy,ifelse(blockx%in%c(96021),25120,     #Shark Bay
                        ifelse(blockx%in%c(96022,96023),26131,
                        ifelse(blockx%in%c(97011),27132,                        #Abrolhos
                        ifelse(blockx%in%c(97012,97013),28132,       
                        ifelse(blockx%in%c(97014,97015),29132,
                        ifelse(blockx%in%c(96010),33151,                        #Geographe Bay
                        ifelse(blockx%in%c(96000),32150,                        #Cockburn sound
                        ifelse(blockx%in%c(96030),35181,blockx)))))))))         # King George sound
DAT$FINYEAR=as.character(DAT$FINYEAR)
DAT$VESSEL=as.character(DAT$VESSEL)
dummy$finyear=as.character(dummy$finyear)
dummy$vessel=as.character(dummy$vessel)
crap=c("FINYEAR","MONTH","BLOCKX","VESSEL")
colnames(dummy)[1:4]=crap
dummy$dupl=with(dummy,paste(VESSEL,crew,MastersName,FINYEAR,MONTH,BLOCKX))
dummy=merge(dummy,subset(DAT,select=c(FINYEAR,MONTH,BLOCKX,VESSEL,zone)),by=crap,all.x=T)
dummy=dummy[!duplicated(dummy$dupl),]

Licenses=dummy[,match(c("VESSEL","zone"),names(dummy))]
Licenses$Unic=with(Licenses,paste(VESSEL,zone))
Licenses=Licenses[!duplicated(Licenses$Unic),]
Licenses=table(Licenses$zone)


#number of crew in current year
Crew=dummy[,match(c("VESSEL","crew","MastersName","zone"),names(dummy))]
N.Crew.min=aggregate(crew~zone+VESSEL,Crew,min,na.rm=T)
N.Crew.min=N.Crew.min[!duplicated(N.Crew.min$VESSEL),]
N.Crew.max=aggregate(crew~zone+VESSEL,Crew,max,na.rm=T)
N.Crew.max=N.Crew.max[!duplicated(N.Crew.max$VESSEL),]

N.Crew.min=aggregate(crew~zone,N.Crew.min,sum)
N.Crew.max=aggregate(crew~zone,N.Crew.max,sum)

N.Crew=merge(N.Crew.min,N.Crew.max,by="zone")
names(N.Crew)=c("zone","min","max")


#fins value  
#note: this applies an average fin price to all sharks with valuable fins using
#       Prices provided by Eva
#dolar.per.kg.fin=35 
PRICES=PRICES%>%mutate(dolar.per.kg=Beach.Price..Adjusted.,
            #          dolar.per.kg=PRICES[,match(paste("uv",substr(Current.yr,3,4),
            #                            substr(Current.yr,6,7),sep=""),names(PRICES))],
                       dolar.per.kg=as.numeric(gsub("\\$", "", dolar.per.kg)),
                       SPECIES=ASA.Species.Code)
            #%>%select(-uv1718)
names(Shark.Fin.Price.List)[grep('Percent',names(Shark.Fin.Price.List))]='Percent'
Shark.Fin.Price.List=Shark.Fin.Price.List%>%
  mutate(Prop=Percent/100)%>%
  dplyr::select(Percent,Prop,species)%>%
  rbind(data.frame(Percent=0, Prop=0, species=27000)) 

Fin.Weight=subset(DAT,SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,RSCommonName,zone,LIVEWT.c))%>%
              filter(!SPECIES%in%c(22998))%>%
              left_join(Shark.Fin.Price.List,by=c('SPECIES'='species'))%>%
              mutate(Per.fin.livewt=Prop,
                     Per.fin.livewt=ifelse(is.na(Per.fin.livewt),Percent.fin.of.livewt,Per.fin.livewt))%>%
              filter(!Per.fin.livewt==0)
Fin.Weight$dolar.per.kg.fin=PRICES[match(22998,PRICES$SPECIES),match('dolar.per.kg',names(PRICES))]
Fin.Weight$fin.weight=Fin.Weight$Per.fin.livewt*Fin.Weight$LIVEWT.c  #using species-specific fin-livewt proportions
#Fin.Weight$fin.weight=Percent.fin.of.livewt*Fin.Weight$LIVEWT.c #previously using blanket 3% (upto 2022 SOFAR)
Fin.Weight$Total.price=Fin.Weight$dolar.per.kg.fin*Fin.Weight$fin.weight
FINS.value=aggregate(Total.price~zone,Fin.Weight,sum,na.rm=T)
FINS.value.species.zone=aggregate(cbind(Total.price,fin.weight)~SNAME+zone,Fin.Weight,sum,na.rm=T)

# Fin.species.composition=Fin.Weight%>%
#   group_by(RSCommonName)%>%
#   summarise(fin.weight=sum(fin.weight))%>%
#   ungroup()%>%
#   mutate(Prop=fin.weight/sum(fin.weight))%>%
#   arrange(-Prop)

#catch value (GVP)
get.yr.price="dolar.per.kg"
PRICES1=PRICES[,match(c("SPECIES",get.yr.price),names(PRICES))]
names(PRICES1)=c("SPECIES","ABARE.dolar.per.kg")

PRICES1=DAT[,match(c("SPECIES","SNAME","zone","LANDWT","LIVEWT.c","CONDITN"),names(DAT))]%>%
          left_join(PRICES1%>%filter(SPECIES%in%unique(DAT$SPECIES)),by="SPECIES")
              

# PRICES1=merge(DAT[,match(c("SPECIES","SNAME","zone","LANDWT","LIVEWT.c","CONDITN"),names(DAT))],
#               PRICES1,by="SPECIES",all.x=T)

# PRICES1=PRICES[,match(c("SPECIES","Weighted.Average.Price.kg","Abare.unit.."),names(PRICES))]
# names(PRICES1)=c("SPECIES","DOF.dolar.per.kg","ABARE.dolar.per.kg")
# PRICES1=merge(DAT[,match(c("SPECIES","SNAME","zone","LIVEWT.c"),names(DAT))],PRICES1,by="SPECIES",all.x=T)

PRICES1=subset(PRICES1,!is.na(LIVEWT.c))

No.price=PRICES1[is.na(PRICES1$ABARE.dolar.per.kg),]
No.price=No.price[!duplicated(No.price$SPECIES),]

Default.price=1
PRICES1$ABARE.dolar.per.kg=with(PRICES1,ifelse(is.na(ABARE.dolar.per.kg),Default.price,ABARE.dolar.per.kg))  #add default if no price data

  #calculate catch price (live weight X price) #As per Eva's calculations of GVP: "I multiple the live weight by 'Final.beach.price'"
PRICES1$Total.price=PRICES1$LIVEWT.c*PRICES1$ABARE.dolar.per.kg
#PRICES1$Total.price=PRICES1$LANDWT*PRICES1$ABARE.dolar.per.kg

CATCH.value=aggregate(Total.price~zone,PRICES1,sum,na.rm=T)

Total.value=merge(CATCH.value,FINS.value,by="zone")
Total.value$Annual.Value.millions=Total.value$Total.price.x+Total.value$Total.price.y
Total.value$Annual.Value.millions=Total.value$Annual.Value.millions/1e6
Total.value=Total.value[,c(1,4)]


#Table 1. catch summary of selected species for current year
Spec.tab.1.elasmo=list(Gummy=Gummy,Dusky_whaler=Dusky_whaler,Whiskery=Whiskery,Sandbar=Sandbar,
                       Hammerheads=Hammerheads,Spinner=Spinner,Wobbegongs=Wobbegongs,Rays=Ray.species,
                       Common_saw_shark=Common_saw_shark,School=School,
                       Other_elasmobranchs=Elasmo.species[-match(c(Gummy,Dusky_whaler,Whiskery,Sandbar,
                                                                   Hammerheads,Spinner,Wobbegongs,Common_saw_shark,School,
                                                                   Ray.species),Elasmo.species)])

This.fish=match(c(Blue_morwong,Blue_groper,West_Australian_dhufish,Pink_snapper,Boarfishes,Samsonfish,Redfishes,Mulloway,Sweetlips,
                  Baldchin_groper),Scalefish.species)
This.fish=This.fish[!is.na(This.fish)]
Spec.tab.1.scalies=list(Blue_morwong=Blue_morwong,Blue_groper=Blue_groper,West_Australian_dhufish=West_Australian_dhufish,
                        Pink_snapper=Pink_snapper,
                        Boarfishes=Boarfishes,Samsonfish=Samsonfish,Redfishes=Redfishes,
                        Mulloway=Mulloway,Sweetlips=Sweetlips,Baldchin_groper=Baldchin_groper,
                        Other_scalefish=Scalefish.species[-This.fish])

fun.Tab1.SoFaR=function(SP)
{
  Dat=DAT  
  Dat$Spec.tab1=NA
  SPECIES=names(SP)
  for(p in 1:length(SPECIES))
  {
    Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
  }
  Dat=subset(Dat,!(is.na(Spec.tab1)))
  TABLA=aggregate(LIVEWT.c~Spec.tab1+Bioregion,data=Dat,sum,na.rm=T)
  TABLA2=aggregate(LIVEWT.c~Spec.tab1+zone,data=Dat,sum,na.rm=T)
  wide <- reshape(TABLA,v.names="LIVEWT.c",timevar="Bioregion",idvar="Spec.tab1",direction="wide")
  colnames(wide)=c("Species","South","West")
  wide=wide[,match(c("Species","West","South"),names(wide))]
  wide2 <- reshape(TABLA2,v.names="LIVEWT.c",timevar="zone",idvar="Spec.tab1",direction="wide")
  colnames(wide2)=c("Species","West.dem","Zone1","Zone2")
  wide=merge(wide,wide2,by="Species")
  wide[,2:ncol(wide)]=round(wide[,2:ncol(wide)])
  wide$Total=rowSums(wide[,2:3],na.rm=T)
  
  return(wide)
}

Table1.SoFar.Elasmos=fun.Tab1.SoFaR(Spec.tab.1.elasmo)
ID.nm=match(c("West","South","West.dem","Zone1","Zone2","Total"),names(Table1.SoFar.Elasmos))
Table1.SoFar.Elasmos[,ID.nm]= Table1.SoFar.Elasmos[,ID.nm]/1000    #convert to tonnes
Table1.SoFar.Scalies=fun.Tab1.SoFaR(Spec.tab.1.scalies)
Table1.SoFar.Scalies[,ID.nm]= Table1.SoFar.Scalies[,ID.nm]/1000


#order and add scientific name
id.elasmo=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads","Spinner","Wobbegongs","Rays","Common_saw_shark","School",
     "Other_elasmobranchs")
Matched=match(id.elasmo,Table1.SoFar.Elasmos$Species)
id.0=id.elasmo[which(is.na(Matched))]
Matched=Matched[!is.na(Matched)]
Table1.SoFar.Elasmos=Table1.SoFar.Elasmos[Matched,]
if(length(id.0)>0)
{
  add.0s=as.data.frame(matrix(id.0,nrow=length(id.0),ncol=ncol(Table1.SoFar.Elasmos)))
  add.0s[,2:ncol(add.0s)]=0
  names(add.0s)=names(Table1.SoFar.Elasmos)
  Table1.SoFar.Elasmos=rbind(Table1.SoFar.Elasmos,add.0s)
}

Group.total=colSums(Table1.SoFar.Elasmos[,2:ncol(Table1.SoFar.Elasmos)],na.rm=T)
Group.total=data.frame(Species="Total Elasmobranchs",rbind(Group.total))
Table1.SoFar.Elasmos=rbind(Table1.SoFar.Elasmos,Group.total)
rownames(Table1.SoFar.Elasmos)=NULL
Table1.SoFar.Elasmos=merge(Table1.SoFar.Elasmos,Order.Elas.Sp.SoFAR,by="Species",all.x=T,sort=F)
THISS=c("Species","Species_or_taxon","Zone1","Zone2","West.dem","South","West","Total")
Table1.SoFar.Elasmos=Table1.SoFar.Elasmos[,match(THISS,names(Table1.SoFar.Elasmos))]

id.scale=c("Blue_morwong","Blue_groper","West_Australian_dhufish","Pink_snapper","Boarfishes","Samsonfish","Redfishes","Mulloway","Sweetlips",
     "Baldchin_groper","Other_scalefish")
Table1.SoFar.Scalies=Table1.SoFar.Scalies[match(id.scale,Table1.SoFar.Scalies$Species),]

Group.total=colSums(Table1.SoFar.Scalies[,2:ncol(Table1.SoFar.Scalies)],na.rm=T)
Group.total=data.frame(Species="Total Scalefish",rbind(Group.total))
Table1.SoFar.Scalies=rbind(Table1.SoFar.Scalies,Group.total)
rownames(Table1.SoFar.Scalies)=NULL
Table1.SoFar.Scalies=merge(Table1.SoFar.Scalies,Order.Scale.Sp.SoFAR,by="Species",all.x=T,sort=F)
Table1.SoFar.Scalies=Table1.SoFar.Scalies[,match(THISS,names(Table1.SoFar.Scalies))]


# add Demersal scalefish suite component (needed for West Coast demersal scalefish Assessment)
WC.Tel=subset(DAT,Bioregion=="WC" & SPECIES%in%Suite)
Tot.Catch.WC.Tel=sum(WC.Tel$LIVEWT.c)/1000

Dem.scale.suit.bio=aggregate(LIVEWT.c~Bioregion,data=subset(DAT,SPECIES%in%Suite),sum)
Dem.scale.suit.bio[,2]=Dem.scale.suit.bio[,2]/1000
Dem.scale.suit.zone=aggregate(LIVEWT.c~zone,data=subset(DAT,SPECIES%in%Suite),sum)
Dem.scale.suit.zone[,2]=Dem.scale.suit.zone[,2]/1000

Table1.SoFar.Scalies=rbind(Table1.SoFar.Scalies,data.frame(Species="Demersal scalefish suite component",Species_or_taxon="",
                                                           Zone1=Dem.scale.suit.zone[2,2],Zone2=Dem.scale.suit.zone[3,2],West.dem=Dem.scale.suit.zone[1,2],
                                                           South=Dem.scale.suit.bio[1,2],West=Dem.scale.suit.bio[2,2],Total=sum(Dem.scale.suit.zone[,2])))


#Add Fishing effort   
C.yr=match(Current.yr,Total.effort.days.monthly$FINYEAR)

#annual
Curr.annual.1000km.gn.hours=Total.effort.hours.monthly[C.yr,]
Curr.annual.km.gn.days=Total.effort.days.monthly[C.yr,]
Curr.annual.km.gn.days[1,2]=Curr.annual.km.gn.days[1,2]*1000

#annual by zone
Curr.annual.1000km.gn.hours.zone=Total.effort.zone.hours.monthly[C.yr,]
Curr.annual.km.gn.days.zone=Total.effort.zone.days.monthly[C.yr,]
Curr.annual.km.gn.days.zone[,2:4]=Curr.annual.km.gn.days.zone[,2:4]*1000

#percentages of effort limits
#annual
Per.annual.lim.1000km.gn.hours=100*Curr.annual.1000km.gn.hours$Total/
  (sum(FishEffLims$Km.Gillnet.Hours.c)/1000)
Per.annual.lim.km.gn.days=100*Curr.annual.km.gn.days$Total/sum(FishEffLims$Km.Gillnet.Days.c)

#annual by zone
Per.annual.lim.1000km.gn.hours.zone=100*Curr.annual.1000km.gn.hours.zone[,2:4]/
  (FishEffLims$Km.Gillnet.Hours.c/1000)
Per.annual.lim.km.gn.days.zone=100*Curr.annual.km.gn.days.zone[,2:4]/FishEffLims$Km.Gillnet.Days.c

#combine in table
add.KM.GN.days="YES"
#add.KM.GN.days="NO"
if(add.KM.GN.days=="NO")
{
  Table1.SoFar.Effort=data.frame(Name=c("Fishing effort (1000 km gn hr)","(percent of reference level)"),
                                 West.dem=NA,Zone1=NA,Zone2=NA,South=NA,West=NA,Total=NA)
  Table1.SoFar.Effort[1,2:7]=c(unlist(Curr.annual.1000km.gn.hours.zone[1,2:4]),NA,NA,Curr.annual.1000km.gn.hours[1,2])
  Table1.SoFar.Effort[2,2:7]=c(Per.annual.lim.1000km.gn.hours.zone,NA,NA,Per.annual.lim.1000km.gn.hours)
  Table1.SoFar.Effort=Table1.SoFar.Effort[,match(c("Name","Zone1","Zone2","West.dem","South","West","Total"),
                                                 names(Table1.SoFar.Effort))]
}

if(add.KM.GN.days=="YES")
{
  Table1.SoFar.Effort=data.frame(Name=c("Fishing effort (km gn d)","","Fishing effort (1000 km gn hr)",""),
                                 West.dem=NA,Zone1=NA,Zone2=NA,South=NA,West=NA,Total=NA)
  Table1.SoFar.Effort[1,2:7]=c(unlist(Curr.annual.km.gn.days.zone[1,2:4]),NA,NA,Curr.annual.km.gn.days[1,2])
  Table1.SoFar.Effort[2,2:7]=c(Per.annual.lim.km.gn.days.zone,NA,NA,Per.annual.lim.km.gn.days)
  Table1.SoFar.Effort[3,2:7]=c(unlist(Curr.annual.1000km.gn.hours.zone[1,2:4]),NA,NA,Curr.annual.1000km.gn.hours[1,2])
  Table1.SoFar.Effort[4,2:7]=c(Per.annual.lim.1000km.gn.hours.zone,NA,NA,Per.annual.lim.1000km.gn.hours)
  Table1.SoFar.Effort=Table1.SoFar.Effort[,match(c("Name","Zone1","Zone2","West.dem","South","West","Total"),
                                                 names(Table1.SoFar.Effort))]
}

Titles=c("Name","Species or taxon","JASDGLF.Zn1","JASDGLF.Zn2","WCDGDLF","Bio.SC","Bio.WC","Total")
names(Table1.SoFar.Elasmos)=names(Table1.SoFar.Scalies)=Titles
names(Table1.SoFar.Effort)=Titles[-2]

Table1.SoFar.Elasmos[,3:ncol(Table1.SoFar.Elasmos)]=round(Table1.SoFar.Elasmos[,3:ncol(Table1.SoFar.Elasmos)],1)
Table1.SoFar.Scalies[,3:ncol(Table1.SoFar.Scalies)]=round(Table1.SoFar.Scalies[,3:ncol(Table1.SoFar.Scalies)],1)
Table1.SoFar.Effort[,2:ncol(Table1.SoFar.Effort)]=round(Table1.SoFar.Effort[,2:ncol(Table1.SoFar.Effort)],1)

Table1.SoFar.Elasmos[is.na(Table1.SoFar.Elasmos)] = ""
Table1.SoFar.Scalies[is.na(Table1.SoFar.Scalies)] = ""
Table1.SoFar.Effort[is.na(Table1.SoFar.Effort)] = ""

Table1.SoFar.Elasmos[Table1.SoFar.Elasmos==0] = "<0.1"
Table1.SoFar.Scalies[Table1.SoFar.Scalies==0] = "<0.1"



#Table 2. catch summary of TEPS
#note: this combines the previous Sofar data with the current data. For current year, must
#      manually scan comments for matching what was reported in columns and in comments
setwd(handle.Sofar)

these.yrs=sort(unique(Data.monthly$FINYEAR))
these.yrs=these.yrs[1:match(Current.yr,these.yrs)]
fun.Tab2.SoFaR=function(Dat)
{
  Dat=Dat%>%
    mutate(Status=ifelse(Status=="a","A",
                  ifelse(Status=="d","D",
                  ifelse(Status=="NR","D",Status))),
           Ali.Ded.yr=paste(finyear,Status),
           CommonName=ifelse(CommonName=='Muttonbird','Shearwater',CommonName))%>%
    filter(finyear%in%these.yrs) 
  
  TABLA=Dat%>%
        group_by(SpeciesCode,Ali.Ded.yr)%>%
        summarise(Number=sum(Number,na.rm=T))%>%
        spread(Ali.Ded.yr,Number,fill='')%>%
        left_join(Dat%>%distinct(SpeciesCode,CommonName),by='SpeciesCode')%>%
        data.frame%>%
        arrange(SpeciesCode)%>% 
        relocate(CommonName, .after = SpeciesCode)
  oldnames=colnames(TABLA)[-(1:2)]
  newnames=gsub("X","",oldnames)
  newnames=gsub(x = newnames, pattern = "A", replacement = "Alive")
  newnames=gsub(x = newnames, pattern = "D", replacement = "Dead")
  TABLA=TABLA%>%
    rename(Species=CommonName)%>%
    mutate(Species=capitalize(tolower(Species)))

  colnames(TABLA)[-(1:2)]=newnames
  
  return(TABLA)
}

TEPS.WC.SC=subset(TEPS.current,fishery%in%c('JASDGDL','OASC,OT','WCDGDL','SGL','WCGL','*'))  


#review comments in case reported in comments and not in logbook columns
if(!file.exists("scan.teps.csv"))   
{
  # look at comments and see if comments match what is recorded in number; 
  #     if not, maybe it is a new interaction record.
  # Remember that comments are duplicated so only use unique record. 
  Scan.current=TEPS.WC.SC%>%
                filter(finyear==Current.yr & !is.na(Comments))%>%
    dplyr::select(DailySheetNumber,Status,Number,CommonName,RSCommonName,Comments)%>%
    arrange(DailySheetNumber)
  write.csv(Scan.current,"scan.teps.csv",row.names=F) 
  cat(paste("scan.teps.csv FILE LOCATED IN",getwd()))
}

#remove sightings from interactions
Sighting.only=which(grepl("sighting", tolower(TEPS.WC.SC$Comments)))
if(length(Sighting.only)>0) TEPS.WC.SC=TEPS.WC.SC[-Sighting.only,]

#tabulate interactions
TEPS=fun.Tab2.SoFaR(Dat=TEPS.WC.SC)
  
#reconcile species
TEPS=TEPS%>%
  mutate(Species=case_when(SpeciesCode==37008001 ~ "Greynurse shark",
                           SpeciesCode==37010003 ~ "White shark",
                           SpeciesCode==37018001 ~ "Bronze whaler",
                           SpeciesCode==37018003 ~ "Dusky shark",
                           SpeciesCode==37018022 ~ "Tiger shark",
                           SpeciesCode==37025000 ~ "Sawfishes",
                           SpeciesCode==37041004 ~ "Giant manta ray",
                           SpeciesCode==37990030 ~ "Stingrays",
                           SpeciesCode==40041050 ~ "Shearwaters",
                           SpeciesCode==41131001 ~ "New zealand fur seal",
                           SpeciesCode==41131005 ~ "Australian Sea lion",
                           SpeciesCode==39125000 ~ "Sea snakes",
                           TRUE~Species))
                                           
#remove whalers
TEPS.whalers=TEPS%>%filter(SpeciesCode%in%37018001:37018025)
TEPS=TEPS%>%filter(!SpeciesCode%in%37018001:37018025)


#Export tables

  #Main features
FIN.yr.slash=paste(substr(Current.yr,1,4),"/",substr(Current.yr,6,7),sep="")
MainFeatures$SPECIES=as.character(MainFeatures$SPECIES)
MainFeatures.doc=rbind(MainFeatures[1:3,],MainFeatures)
MainFeatures.doc[1,]=c(paste("Current Landings (",FIN.yr.slash,")",sep=""),"")
MainFeatures.doc[2,]=c("Demersal Gillnet and Demersal Longline Fishery","")
MainFeatures.doc[3,]=c("Indicator species","")
OrdeR=c(paste("Current Landings (",FIN.yr.slash,")",sep=""),
        "Demersal Gillnet and Demersal Longline Fishery","Total sharks and rays",
        "Total scalefish","Indicator species","Gummy","Bronzy.Dusky",                                        
        "sandbar","Whiskery","Sharks and rays by other commercial fisheries",
        paste("Reconstructed statewide recreational catch ","(",REC.yr,")",sep=""))
MainFeatures.doc=MainFeatures.doc[match(OrdeR,MainFeatures.doc$SPECIES),]
MAIN.F.sp=c("Gummy shark","Dusky shark","Sandbar shark","Whiskery shark")
MainFeatures.doc$SPECIES[match(c("Gummy","Bronzy.Dusky","sandbar","Whiskery"),
                               MainFeatures.doc$SPECIES)]=MAIN.F.sp

MainFeatures.doc.status=MainFeatures.doc
MainFeatures.doc.status$SPECIES=c("Status","Stock level",MAIN.F.sp,"Fishing Level",
                                  "JASDGDLF Zone 1","JASDGDLF Zone 2","WCDGDLF","")            

fn.stat.fishing=function(zn)
{
  iddz=match(zn,colnames(Table1.SoFar.Effort))
  return(ifelse(max(c(Table1.SoFar.Effort[4,iddz],Table1.SoFar.Effort[2,iddz]))>100,"Unacceptable","Acceptable")) 
}
Fishing.Zone1.status=fn.stat.fishing("JASDGLF.Zn1")
Fishing.Zone2.status=fn.stat.fishing("JASDGLF.Zn2")
Fishing.WC.status=fn.stat.fishing("WCDGDLF")
MainFeatures.doc.status$Catch.tons=c("","",
                                     Gummy.status,Dusky.status,Sandbar.status,Whiskery.status,"",
                                     Fishing.Zone1.status,Fishing.Zone2.status,Fishing.WC.status,"")

MainFeatures.doc=cbind(MainFeatures.doc.status,MainFeatures.doc)
write.csv(MainFeatures.doc,"1.1.MainFeatures.csv",row.names=F)


#Retained Species
OrdeR=c("Total sharks and rays","Total indicators","Gummy","Bronzy.Dusky","Whiskery","sandbar")
Retained.species=MainFeatures[MainFeatures$SPECIES%in%OrdeR,]
Retained.species=Retained.species[match(OrdeR,Retained.species$SPECIES),]
Retained.species$Catch.tons=paste(Retained.species$Catch.tons,"onnes",sep="")

Change=c("All sharks (and rays):","Indicator shark species:","Gummy:","Dusky2:","Whiskery:","Sandbar:")
Retained.species$SPECIES[match(OrdeR,Retained.species$SPECIES)]=Change
write.csv(Retained.species,"Retained.species.csv",row.names=F)


#Export if catch levels are Acceptable or Not 
Ktch_comm_accept=Tot.wt%>%mutate(SPECIES=ifelse(SPECIES=="sandbar","Sandbar",SPECIES))%>%
                          filter(SPECIES%in%c("Gummy","Whiskery","Bronzy.Dusky","Sandbar"))%>%
                          dplyr::select(SPECIES,Catch.tons)%>%
                          left_join(Catch.range.key.species,by="SPECIES")%>%
                          mutate(Acceptable=ifelse(Catch.tons>=Min.catch*.9 & Catch.tons<=Max.catch*1.1,"Acceptable",
                                                   "Unacceptable"))
write.csv(Ktch_comm_accept,"Ktch_comm_accept.csv",row.names=F)


#Export tables as .csv
Yr.Rec.Ktch=max(Rec.fish.catch$FINYEAR)
Rec.Catch.stats=Rec.fish.catch%>%filter(FINYEAR==Yr.Rec.Ktch)%>%
  group_by(Common.Name)%>%summarise(Total.catch.kg=sum(LIVEWT.c))
write.csv(Rec.Catch.stats,paste("1.2.Rec.Catch.stats.",Yr.Rec.Ktch,".csv",sep=""),row.names=F)

write.csv(Licenses,"1.4.Number.of.vessels.csv",row.names=F)
write.csv(N.Crew,"1.5.Number.of.crew.csv",row.names=F)
write.csv(Total.value,"1.6.EconomicEffects.csv",row.names=F)
write.csv(FINS.value.species.zone,"1.7.FINS.value_only.in.dollars.csv",row.names=F)


write.table(Table1.SoFar.Elasmos,"2.Table1.SoFar.Elasmos.csv",sep = ",",row.names = F)
write.table(Table1.SoFar.Scalies,"2.Table1.SoFar.Scalies.csv",sep = ",",row.names = F)
write.table(Table1.SoFar.Effort,"2.Table1.SoFar.Effort.csv",sep = ",",row.names = F)

if(file.exists("scan.teps.csv"))
{
  write.table(TEPS,"3.Table2.TEPS.csv",sep = ",",row.names = F)
  write.table(TEPS.whalers,"3.TEPS_whalers.csv",sep = ",",row.names = F)
 #  write.csv(TEPS,handl_OneDrive("Data/Catch and Effort/Historic/Historic.TEPS.res.csv"),row.names = F)  
}

#Prices of main species
write.table(PRICES%>%filter(ASA.Species.Code%in%c(17001,17003,18003,18007)),"4.Table.Beach Prices_indicator species.csv",sep = ",",row.names = F)

#Figures 2-3. 
Boundary.blk=c(34160,35160,36160) #Boundary blocks
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
par(mfcol=c(1,1),mar=c(3.5,3.6,.1,1),oma=c(1,.5,.1,.1))
LINE=c(5,1,1)
TYPE=c("l","o","o")
PCH=21
COL=1
BG=c("black","black","white")


fn.figs2.3.SoFaR=function(GROUP,LAT1,LAT2,INT,INT2)
{
    dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2 & METHOD%in%c("GN","LL") &
                 Estuary=="NO")%>%
      filter(!Shark.fishery=='non.shark.fishery')
    
    #Split boundary zone1-zone2 catch
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone+BLOCKX,data=dat,sum,na.rm=T)
    annual.catch.by.zone=fn.split.boundary(annual.catch.by.zone,"LIVEWT.c")
    
    #aggregate by zone
    annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=annual.catch.by.zone,sum,na.rm=T)
    
    #aggregate by total
    annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
    
    wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
    
    names(wide)[match("FINYEAR",names(wide))]="finyear"
    names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
    
    fun.fig.SoFar(DAT=annual.catch.total,DAT1=wide,scaler=1000,
                  TITLE1="Catch (tonnes live wt.)",TITLE2="Financial year",INT,INT2)
  }


#Figure 2   
LEG=function()legend("topright",c("TDGDLF","WCDGDLF","JASDGDLF Zone1","JASDGDLF Zone2"),lty=c(1,LINE),lwd=1.75,
                     col=c("grey80",rep("black",3)),pch=c(NA,NA,PCH,PCH),pt.bg=c(NA,NA,BG[2:3]),bty='n',cex=1.5)
PaR=function()par(mai=c(.9,1,.1,.15),oma=c(.01,.5,.01,.1),mgp=c(3.3,.7,0),las=1)
CeX=3
WHere=3
DO="YES"

jpeg(file="Figure 2.TotalElasmoCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
PaR()
fn.figs2.3.SoFaR(GROUP=Elasmo.species,LAT1=TDGDLF.lat.range[1],LAT2=TDGDLF.lat.range[2],INT=100,INT2=500)
LEG()
dev.off()


#Figure 3
#note: Monthly data provided by Rory has no scalefish so use a specifc function that uses
#       the historic excel file
Yr.current=2000+as.numeric(substr(Current.yr,6,7))
fn.figs2.3.SoFaR.scalies=function(GROUP,LAT1,LAT2,INT,INT2)
{
  a=2011:(Yr.current-1);b=2012:Yr.current
  Curr.yr=paste(a,"-",substr(b,3,4),sep="")
  DAT=subset(Data.monthly,FINYEAR%in%Curr.yr & METHOD%in%c("GN","LL") & Estuary=="NO" &
               LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])%>%
    filter(!Shark.fishery=="non.shark.fishery")
  
  dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
  
  annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
  annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
  
  wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
  
  Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.tel.live.wt","Z1.tot.tel.live.wt","Z2.tot.tel.live.wt"),
                                  names(Results.pre.2013))]
  Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.tel.live.wt"),names(Results.pre.2013))]
  
  names(Prev.zn)=names(wide)
  Prev.zn[,2:4]=Prev.zn[,2:4]*1000
  wide=rbind(Prev.zn,wide)
  
  
  Prev[,2]=Prev[,2]*1000
  names(Prev)=names(annual.catch.total)
  annual.catch.total=rbind(Prev,annual.catch.total)
  
  names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
  
  fun.fig.SoFar(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
}
jpeg(file="Figure 3.TotalScalefishCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
PaR()
fn.figs2.3.SoFaR.scalies(GROUP=Scalefish.species,LAT1=TDGDLF.lat.range[1],LAT2=TDGDLF.lat.range[2],INT=20,INT2=100)
LEG()
dev.off()


#Figure 4
fn.eff=function(DAT,DAT1,TITLE1,TITLE2,INT,INT2)
{
  MAX=max(DAT[,2],na.rm=T)
  FInYEAR=as.character(unique(DAT$FINYEAR))
  N=length(FInYEAR)
  
  #id=match(start.yr,FInYEAR)
  id=which.min(is.na(DAT[,2]))
  FInYEAR=FInYEAR[id:length(FInYEAR)]
  NN=length(FInYEAR)
  
  
  plot(1:NN,DAT[id:N,2],type='l',col="grey80",ylim=c(0,MAX),xaxt='n',yaxt='n',
       ylab="", xlab="",las=1,lwd=2,cex.lab=1.3)
  axis(1,at=1:NN,labels=F,tck=-0.01)
  axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.35)
  
  axis(2,at=seq(0,MAX,INT),labels=F,tck=-0.01)
  axis(2,at=seq(0,MAX,INT2),labels=seq(0,MAX,INT2),tck=-0.02,cex.axis=1.35)
  mtext(TITLE1,2,las=3,cex=CeX,line=3.2)
  mtext(TITLE2,1,cex=CeX,line=WHere)
  
  
  for(i in 1:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1],type=TYPE[i],lty=LINE[i],col=COL,lwd=1.5,pch=PCH,bg=BG[i])
}
jpeg(file="Figure 4.StandardisedEffort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
PaR()
fn.eff(Total.effort.days.monthly,Total.effort.zone.days.monthly,"Effort (1000 km gn d)","Financial year",2,10)
LEG()
dev.off()

#Figure Catch with Effort
CeX=2
WHere=2
jpeg(file="Figure TotalElasmoCatch_Effort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
par(mfcol=c(2,1),mar=c(1,4.5,1,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2,.6,0))
DO="NO"
fn.figs2.3.SoFaR(Elasmo.species,TDGDLF.lat.range[1],TDGDLF.lat.range[2],100,500)
DO="YES"
fn.eff(Total.effort.days.monthly,Total.effort.zone.days.monthly,"Effort (1000 km gn d)","Financial year",2,10)
LEG()
dev.off()



#Figures 5-8. Standardised cpue (zones combined)
  #read in data
HNDL=handl_OneDrive("Analyses/Data_outs/")
#Whiskery
whis.mon=read.csv(paste(HNDL,"Whiskery shark/Whiskery Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
whis.daily=read.csv(paste(HNDL,"Whiskery shark/Whiskery Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#gummy
gum.mon=read.csv(paste(HNDL,"Gummy Shark/Gummy Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
gum.daily=read.csv(paste(HNDL,"Gummy Shark/Gummy Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#dusky
dus.mon=read.csv(paste(HNDL,"Dusky Shark/Dusky Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
dus.daily=read.csv(paste(HNDL,"Dusky Shark/Dusky Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#sandbar
san.mon=read.csv(paste(HNDL,"Sandbar Shark/Sandbar Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
san.daily=read.csv(paste(HNDL,"Sandbar Shark/Sandbar Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#plot cpues
Plot.cpue.delta=function(cpuedata,cpuedata.daily,CL,CxS,
                         Yvar,add.lines,firstyear,ADD.nomnl)    #plot cpues
{
  if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
  if(length(cpuedata)<=3)tc=seq(-.5*0.25,.5*0.25,length.out=length(cpuedata))
  
  Yrs=c(as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4)),
        as.numeric(substr(cpuedata.daily[[1]][,match(Yvar,names(cpuedata.daily[[1]]))],1,4)))
  Tops=c(unlist(lapply(cpuedata, `[`, "upper.CL")),
         unlist(lapply(cpuedata.daily, `[`, "upper.CL")),
         ADD.nomnl$response)
  ymax=max(Tops)
  Quant=quantile(Tops,probs=c(.9,1))
  if(diff(Quant)>3) ymax=quantile(Tops,probs=.99)
  
  
  plot(Yrs,Yrs,ylim=c(0,ymax),xlim=c(firstyear,max(Yrs)),ylab="",xlab="",
       col="transparent",cex.axis=1.25)
  for(l in 1:length(cpuedata))
  {
    aaa=cpuedata[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
    aaa.daily=cpuedata.daily[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
    
    msn=Yrs[which(!Yrs%in%c(aaa$finyear,aaa.daily$finyear))]
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
           if(add.lines=="NO") points(finyear+tc[l], response, pch=19, lty=2, col=CL[l],cex=CxS)
           if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=19, lty=2, col=CL[l],cex=CxS)
           arrows(x0=finyear+tc[l], y0=lower.CL, 
                  x1=finyear+tc[l], y1=upper.CL, 
                  code=3, angle=90, length=0.05, col=CL[l])
         })
    with(aaa.daily,
         {
           if(l==1)polygon(x=c(finyear[1]-.5,finyear[length(finyear)]+.5,finyear[length(finyear)]+.5,finyear[1]-.5),
                           y=c(0,0,ymax*.99,ymax*.99),col='grey92',border="transparent")
           arrows(x0=finyear+tc[l], y0=lower.CL, 
                  x1=finyear+tc[l], y1=upper.CL, 
                  code=3, angle=90, length=0.05, col=CL[l])
           if(add.lines=="NO") points(finyear+tc[l], response, pch=21, lty=2, col=CL[l],cex=CxS)
           if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=21,bg="white", lty=2, col=CL[l],cex=CxS)
         })
  }
  if(!is.null(ADD.nomnl))
  {
    with(ADD.nomnl,points(as.numeric(substr(finyear,1,4))+.1,response,pch=19,col=rgb(.1,.1,.1,alpha=.2),cex=CxS)) 
  }
  
}

Pred.creep=list(whis.mon,gum.mon,dus.mon,san.mon)
Pred.daily.creep=list(whis.daily,gum.daily,dus.daily,san.daily)
SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")

jpeg(file="Figure All.cpues.jpeg",width = 1800, height = 2400,units = "px", res = 300)
par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
for(s in 1:length(Pred.creep))
{
  dummy=Pred.creep[[s]]
  colnames(dummy)=tolower(colnames(dummy))
  dummy=dummy%>%rename(response=mean,upper.CL=up.ci,lower.CL=low.ci)
  dummy.daily=Pred.daily.creep[[s]]
  colnames(dummy.daily)=tolower(colnames(dummy.daily))
  dummy.daily=dummy.daily%>%rename(response=mean,upper.CL=up.ci,lower.CL=low.ci)
  
  Plot.cpue.delta(cpuedata=list(Standardised=dummy),
                  cpuedata.daily=list(Standardised=dummy.daily),
                  CL=1,CxS=1.5,Yvar="finyear",add.lines="YES",firstyear=1975,
                  ADD.nomnl=NULL)
  legend("topright",SPECIES.vec[s],bty='n',cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("RELATIVE CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
dev.off()

#Number of trips using longline or gillnet by year
a=Effort.monthly%>%
  filter(METHOD%in%c("LL","GN"))%>%
  distinct(FINYEAR,Same.return,METHOD)%>%
  mutate(N=1,
         Yr=as.numeric(substr(FINYEAR,1,4)))
b=Effort.daily%>%
  filter(method%in%c("LL","GN"))%>%
  distinct(finyear,Same.return,method)%>%
  mutate(N=1,
         Yr=as.numeric(substr(finyear,1,4)))%>%
  rename(METHOD=method,
         FINYEAR=finyear)

jpeg(file="Proportion.trips.GN.LL.jpeg",width=2400,height=2200,units="px",res=300)
rbind(a,b)%>% 
  group_by(Yr,METHOD)%>%
  summarise(N=sum(N))%>%
  mutate(METHOD=case_when(METHOD=="GN"~"Gillnet",
                          METHOD=="LL"~"Longline"))%>%
  ggplot(aes(fill=METHOD, y=N, x=Yr)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Number of trips")+
  xlab("Financial year")+
  theme(axis.text= element_text( size = 14),
        axis.title=element_text( size = 20),
        legend.text = element_text(size = 13),
        legend.title=element_blank(),
        legend.position="top",
        plot.margin=unit(c(.5,.5,.5,.5),"cm"))
dev.off()


# Is reduction in catch due to reduction in number of shots per vessel
fn.check.drop.in.ktch=function(years,methodS)
{
  DaT=Data.daily%>%
    filter(FINYEAR%in%years & METHOD%in%methodS & Estuary=="NO" &
             FisheryCode%in%c('JASDGDL','WCDGDL') &
             LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  top10=sort(table(DaT$SNAME))
  DaT=DaT%>%
    mutate(SP=ifelse(!SNAME%in%names(tail(top10,10)),'Other',SNAME))
  
  Annual.N.shot.by.vessel=DaT%>%
    distinct(Same.return.SNo,.keep_all = T)%>%
    group_by(Same.return.SNo,VESSEL,LongFC,LatFC,METHOD,FINYEAR,FisheryZone)%>%
    tally()%>%
    ungroup()%>%
    group_by(METHOD,FINYEAR,FisheryZone,VESSEL)%>%
    summarise(n=sum(n))
  
  Annual.Ktch.by.vessel=DaT%>%
    group_by(METHOD,FINYEAR,FisheryZone,VESSEL)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T)/1000)%>%
    ungroup()
  
  p1=Annual.N.shot.by.vessel%>%
    left_join(Annual.Ktch.by.vessel,by=c('METHOD','FINYEAR','FisheryZone','VESSEL'))%>%
    mutate(Zone.method=paste(FisheryZone,METHOD,sep='-'),
           year=as.numeric(substr(FINYEAR,1,4)))%>%
    ggplot(aes(year,n,color=Zone.method))+
    geom_point(aes(size=Tonnes),alpha=0.5)+
    geom_line(alpha=0.5)+
    facet_wrap(~VESSEL)+
    ylab('Total number of shots')+xlab('Financial year')+
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.75))+
    scale_color_manual(values = c("West-GN"="red","West-LL"="brown",
                                  "Zone1-GN"="forestgreen","Zone1-LL"="green",
                                  "Zone2-GN"="blue","Zone2-LL"="steelblue"))+
    labs(title="Number of shots by fisher per financial year",
         subtitle = "Bubble size is proportional to total catch")
  
  Annual.Ktch.per.sp=DaT%>%
    group_by(SP,METHOD,FINYEAR,FisheryZone,VESSEL)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T)/1000)%>%
    ungroup()
  p2=Annual.Ktch.per.sp%>%
    mutate(Zone.method=paste(FisheryZone,METHOD,sep='-'),
           year=as.numeric(substr(FINYEAR,1,4)))%>%
    ggplot(aes(year,Tonnes,color=SP))+
    geom_point(size=2)+
    geom_line()+
    facet_wrap(~VESSEL)+
    ylab('Total catch (tonnes)')+xlab('Financial year')+
    theme(legend.position = 'top',
          legend.title=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.75))+
    scale_color_manual(values = c("Blue Morwong"="cadetblue4","Boarfishes"="cadetblue2","Pink Snapper"="aquamarine3", 
                                  "West Australian Dhufish"="chartreuse3", "Western Blue Groper"="chartreuse4",
                                  "Bronze Whaler"="chocolate4", "Dusky Whaler"="chocolate","Gummy Shark"="red",
                                  "Hammerhead Sharks"="darkorange", "Whiskery Shark"="brown4",
                                  "Other"="darkorchid2"))+
    labs(title="Total catch by fisher per financial year and species")
  print(p1)
  print(p2)
  
}
nn1=as.numeric(substr(Current.yr,1,4))

pdf(file="Shots and catch per vessel for last 5 years.pdf")
fn.check.drop.in.ktch(years=paste((nn1-4):nn1,substr((nn1-3):(nn1+1),3,4),sep='-'),
                      methodS=c("GN","LL"))
dev.off()


# For ERA -----------------------------------------------------------------
#Create data for TDGDLF Ecological Risk Assessment 
do.ERA=TRUE   #needed for Reconstructions of TDGDLF discards
if(do.ERA)
{
  era.yrs=as.numeric(substr(Current.yr,1,4))
  era.yrs1=seq(era.yrs-3,era.yrs+1)
  era.yrs=seq(era.yrs-4,era.yrs)
  era.yrs=paste(era.yrs,substr(era.yrs1,3,4),sep='-')
  ERA=subset(Data.monthly,METHOD%in%c("GN","LL") & Estuary=="NO" &
               LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])%>%
    filter(!Shark.fishery=="non.shark.fishery")%>%
    filter(FINYEAR%in%era.yrs)%>%
    dplyr::select(FINYEAR,SPECIES,LIVEWT.c)%>%
    group_by(FINYEAR,SPECIES)%>%
    summarise(Tons=sum(LIVEWT.c)/1000)%>%
    group_by(SPECIES)%>%
    mutate(Average.catch=mean(Tons))%>%
    spread(FINYEAR,Tons,fill = 0)%>%
    arrange(-Average.catch)%>%
    ungroup()%>%
    mutate(Percent.of.total.retained=round(100*Average.catch/sum(Average.catch),1))%>%
    left_join(All.species.names%>%
                distinct(CAES_Code,.keep_all=T)%>%
                dplyr::select(CAES_Code,COMMON_NAME,SCIENTIFIC_NAME), 
              by=c('SPECIES'='CAES_Code'))%>%
    mutate(Percent.of.total.retained=ifelse(Percent.of.total.retained<0.1,
                                            '<0.1',Percent.of.total.retained))%>%
    filter(!is.na(COMMON_NAME))%>%
    dplyr::select(COMMON_NAME,SCIENTIFIC_NAME,all_of(era.yrs),Average.catch,Percent.of.total.retained)%>%
    mutate(across(where(is.numeric), round, 2),
           COMMON_NAME=ifelse(COMMON_NAME=="Skates and rays, other","Other skates and rays",COMMON_NAME))%>%
    filter(Average.catch>0)
  
  
  #1. Catch
  #Table of retained species
  write.csv(ERA,"ERA_table.retained.species.csv",row.names=F)
  
  #Annual catches main species
  ERA.main.SP=data.frame(SPECIES=c(17001,18003,17003,18001,
                                   19000,18007,18023,13000,
                                   384002,377004,353001,320000),
                         Names=c('Gummy shark', 'Dusky shark', 'Whiskery shark', 'Bronze whaler',
                                 'Hammerheads','Sandbar shark','Spinner shark','Wobbegongs',
                                 'Blue groper','Queen snapper','Pink snapper','West Australian dhufish'))
  ERA.main.SP$Names=factor(ERA.main.SP$Names,levels=ERA.main.SP$Names)
  
  #myColors=c(brewer.pal(8,"YlOrRd"),brewer.pal(4,"Blues"))
  myColors=c("brown","brown1","chocolate","deeppink2","pink","coral2","goldenrod","darkorange",
             "blue","cyan3","deepskyblue","cornflowerblue")
  names(myColors)=ERA.main.SP$Names
  jpeg(file="ERA_Catch_main.species.jpeg",width=2400,height=2200,units="px",res=300)
  Data.monthly%>%
    filter(METHOD%in%c("GN","LL") & Estuary=="NO" & SPECIES%in%ERA.main.SP$SPECIES  &
             LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])%>%
    filter(!Shark.fishery=="non.shark.fishery")%>%
    left_join(ERA.main.SP,by='SPECIES')%>%
    group_by(FINYEAR,Names)%>%
    summarise(Catch.tons=sum(LIVEWT.c)/1000)%>%
    mutate(Yr=as.numeric(substr(FINYEAR,1,4)))%>%
    ggplot(aes(Yr,Catch.tons))+
    geom_line(aes(colour=Names),linewidth=1.5)+
    ylab("Catch (tonnes live wt.)")+
    xlab("Financial year")+
    theme(axis.text= element_text( size = 14),
          axis.title=element_text( size = 20),
          legend.text = element_text(size = 13),
          legend.title=element_blank(),
          legend.position="top",
          plot.margin=unit(c(.5,.5,.5,.5),"cm"))+
    scale_colour_manual(name = 'Names',values = myColors)
  dev.off()
  
  #2. Effort
  # Gillnet and longline equivalent effort
  LL.equiv.Eff.days.total=LL.equiv.Eff.days.zone%>%
    mutate(Total=LIVEWT.West + LIVEWT.Zone1 + LIVEWT.Zone2,
           Data="Longline equivalent")%>%
    dplyr::select(FINYEAR,Total,Data)
  d=Total.effort.days.monthly%>%
    mutate(Total=Total-LL.equiv.Eff.days.total$Total,
           Data="Gillnet")
  jpeg(file="ERA_Effort_Gillnet and longline equivalent.jpeg",width=2400,height=2400,units="px",res=300)
  rbind(d,LL.equiv.Eff.days.total)%>%
    mutate(Year=as.numeric(substr(FINYEAR,1,4)))%>%
    ggplot(aes(Year,Total))+
    geom_line(aes(colour=Data),linewidth=3)+
    ylab("Effort (1000 km gn d)")+
    xlab("Finacial year")+
    theme(axis.text= element_text( size = 14),
          axis.title=element_text( size = 20),
          legend.text = element_text(size = 16),
          legend.title=element_blank(),
          legend.position="top",
          plot.margin=unit(c(.75,1.5,1,.75),"cm"))
  dev.off()
  
  write.csv(100*LL.equiv.Eff.days.total$Total/Total.effort.days.monthly$Total,
            "ERA_Effort_Longline.equiv.eff_percent.total.csv",row.names = Total.effort.days.monthly$FINYEAR)
  
  
  # Gillnet and longline effort
  LL.effort.monthly=Effort.monthly%>%
    filter(METHOD=='LL' & LAT<=(-26) & LAT>(-36.1) & LONG<=(129) & LONG>(111))%>%
    mutate(hook.days=BDAYS.c*HOOKS)%>%
    group_by(VESSEL,BLOCKX,FINYEAR,MONTH,YEAR.c)%>%
    summarise(hook.days=max(hook.days,na.rm=T))%>%
    group_by(FINYEAR)%>%
    summarise(hook.days=sum(hook.days,na.rm=T))
  
  LL.effort.daily=Effort.daily%>%
    mutate(LAT=-abs(LAT))%>%
    filter(method=='LL' & LAT<=(-26) & LAT>(-36.1) & LONG<=(129) & LONG>(111))%>%
    mutate(hook.days=hooks)
  Use.Date="YES"
  if(Use.Date=="NO")
  {
    LL.effort.daily=LL.effort.daily%>%
      group_by(ID,vessel,finyear,month,blockx)%>%
      summarise(hook.days=max(hook.days,na.rm=T))
  }
  if(Use.Date=="YES")
  {
    LL.effort.daily=LL.effort.daily%>%
      group_by(date,vessel,finyear,month,blockx)%>%
      summarise(hook.days=max(hook.days,na.rm=T))
  }
  LL.effort.daily=LL.effort.daily%>%
    group_by(finyear,month,vessel,blockx)%>%
    summarise(hook.days=sum(hook.days,na.rm=T))%>%
    group_by(finyear)%>%
    summarise(hook.days=sum(hook.days,na.rm=T))
  colnames(LL.effort.daily)=colnames(LL.effort.monthly)
  
  LL.effort=rbind(LL.effort.monthly,LL.effort.daily)%>%
    group_by(FINYEAR)%>%
    summarise(thousand.hook.days=sum(hook.days)/1000)%>%
    mutate(Yr=as.numeric(substr(FINYEAR,1,4)))
  
  coeff=30
  jpeg(file="ERA_Effort_Gillnet and longline.jpeg",width=2400,height=2400,units="px",res=300)
  full_join(d%>%mutate(Yr=as.numeric(substr(FINYEAR,1,4)))%>%dplyr::select(-Data),
            LL.effort,by=c('FINYEAR','Yr'))%>%
    mutate(thousand.hook.days=ifelse(is.na(thousand.hook.days),0,thousand.hook.days))%>%
    ggplot(aes(x=Yr))+
    geom_line(aes(y=Total), size=3, color="#F8766D")+
    geom_line(aes(y=thousand.hook.days/coeff),linewidth=3,col="#00BFC4")+
    scale_y_continuous(name = "Gillnet effort (1000 km gn d)",
                       sec.axis = sec_axis(~.*coeff, name="Longline effort (1000 hook d)"))+
    xlab("Finacial year")+
    theme(axis.text= element_text( size = 14),
          axis.title=element_text( size = 20),
          legend.text = element_text(size = 16),
          axis.title.y = element_text(color = "#F8766D", size=18),
          axis.title.y.right = element_text(color = "#00BFC4", size=18)
    )
  dev.off()
  
}

# For AMM -----------------------------------------------------------------
do.AMM=TRUE
if(do.AMM)
{
  do.mngmnt.timeline=FALSE
  if(do.mngmnt.timeline)
  {
    library(stringr)
    library(ggpubr)
    Management=read.csv(handl_OneDrive('Management/Sharks/Timeline management measures/Management_timeline.csv'))
    Management=Management%>%
      mutate(date=as.POSIXct(StartDate,format="%d/%m/%Y"),
             finyear=decimal_date(date))
    Management.north=Management%>%
      filter(Relevant.to%in%c('NSF',"TDGDLF & NSF") & 
               round(finyear)<=as.numeric(substr(Current.yr,1,4)))
    
    Management=Management%>%
      filter(Relevant.to%in%c('TDGDLF',"TDGDLF & NSF") & 
               round(finyear)<=as.numeric(substr(Current.yr,1,4)))
    
    
    #1. Catch and management measures timeline
    AMM.ktch.eff.managmtnt=function(GROUP,LAT1,LAT2,MGMT,SEP,Lab.font,Lab.back,
                                    nn,frms,n.frms,do.animation,speed)
    {
      dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2 & METHOD%in%c("GN","LL") &
                   Estuary=="NO")%>%
        filter(!Shark.fishery=="non.shark.fishery")
      
      #aggregate by total
      annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)%>%
        mutate(finyear=1+as.numeric(substr(FINYEAR,1,4)),
               LIVEWT.c=LIVEWT.c/1000)%>%
        dplyr::select(-FINYEAR)
      
      dummy=annual.catch.total[1:nrow(MGMT),]%>%
        mutate(LIVEWT.c=NA,
               finyear=MGMT$finyear)%>%
        filter(!finyear%in%annual.catch.total$finyear)
      annual.catch.total=rbind(annual.catch.total,dummy)
      
      Man=MGMT%>%dplyr::select(Event,finyear,Category)%>%mutate(id=1:nrow(MGMT))
      
      d=annual.catch.total%>%
        left_join(Man,by='finyear')%>%
        arrange(finyear)%>%
        filter(finyear>=1975)%>%
        mutate(LIVEWT.c=ifelse(is.na(LIVEWT.c),na.approx(LIVEWT.c),LIVEWT.c),
               Category=ifelse(is.na(Category),"",Category),
               Event=ifelse(is.na(Event),"",Event))
      
      set.seed(666)
      colors=c(Data='red',Closure='steelblue',Other='forestgreen','NA'='transparent')
      p=ggplot(d,
               aes(finyear, LIVEWT.c,label=Event,color = factor(Category))) +
        geom_line(colour="orange",size=1.25) + geom_point() +
        geom_label_repel(box.padding=SEP,hjust = 0,size = Lab.font,fill=Lab.back) + 
        geom_line(colour="orange",size=1.25, alpha = 0.3) + geom_point(alpha = 0.4)+
        ylab("Cach (tonnes live wt)") + xlab("Financial year")+
        theme(legend.position = "none",
              axis.text=element_text(size=12),
              axis.title=element_text(size=14)
              # panel.background = element_rect(fill = "black", color  =  NA),
              # panel.border = element_rect(fill = NA, color = "white"),  
              # panel.grid.major = element_line(color = "grey20"),  
              # panel.grid.minor = element_line(color = "grey20"),  
              # panel.spacing = unit(0.5, "lines")
        )+
        scale_color_manual(values = colors)
      p
      ggsave("AMM_catch.management.tiff", width = 8,height = 5, dpi = 300, compression = "lzw")
      
      
      
      #animation
      if(do.animation)
      {
        Min.yr=min(dat$YEAR.c)
        
        a=dat%>%
          mutate(finyear=1+as.numeric(substr(FINYEAR,1,4)))%>%
          group_by(finyear)%>%
          summarise(LIVEWT.c=sum(LIVEWT.c)/1000)
        
        b=a[rep(1:nrow(a),nn),]%>%
          arrange(finyear)
        SEQ=seq(1,nrow(b),by=nn)
        b$finyear[-SEQ]=NA
        b$LIVEWT.c[-SEQ]=NA
        b$finyear=with(b,ifelse(is.na(finyear),na.approx(finyear),finyear))
        b$LIVEWT.c=with(b,ifelse(is.na(LIVEWT.c),na.approx(LIVEWT.c),LIVEWT.c))
        
        Man1=b%>%
          full_join(Man%>%filter(finyear>=Min.yr),by='finyear')%>%
          filter(finyear>=Min.yr)%>%
          mutate(Event=ifelse(is.na(Event),'',Event),
                 Category=ifelse(is.na(Category),'',Category))%>%
          group_by(Event)%>%
          mutate(n=n())%>%
          mutate(temp = case_when((n >= n.frms) ~ (1),
                                  (n<n.frms) ~ (n.frms))) %>%
          uncount(temp)%>%
          dplyr::select(-c(id,n))%>%
          arrange(finyear)%>%
          data.frame
        #dummy=Man1[rep(1,5),]
        #dummy$Event=dummy$Category=''
        #dummy$finyear=seq(dummy$finyear[1]-.5,dummy$finyear[1]-.99,length.out=nrow(dummy))
        #Man1=rbind(dummy,Man1)
        Man1$ID=1:nrow(Man1)
        Man1$LIVEWT.c=with(Man1,ifelse(is.na(LIVEWT.c),na.approx(LIVEWT.c, rule=2),LIVEWT.c))
        Man1$Event.yr=with(Man1,ifelse(Event!='',paste(round(finyear),Event,sep=': '),''))
        anim <- ggplot(Man1, aes(finyear, LIVEWT.c,label=Event.yr)) +
          geom_line(colour="orange",size=3) +
          geom_label_repel(box.padding=SEP,hjust = 0,size = 9,fill=Lab.back)+
          theme(legend.position = "none",
                axis.text=element_text(size=20),
                axis.title=element_text(size=24))+
          scale_color_manual(values = colors)+
          ylab("Cach (tonnes live wt)") + xlab("Financial year")+
          transition_reveal(as.numeric(ID))
        p.animate=animate(anim,nframes=frms,duration = speed, height = 700, width = 1000)
        anim_save("AMM_catch.management.gif", p.animate)
      }
      
    }
    AMM.ktch.eff.managmtnt(GROUP=Elasmo.species,
                           LAT1=TDGDLF.lat.range[1],
                           LAT2=TDGDLF.lat.range[2],
                           MGMT=Management,
                           SEP=.1,
                           Lab.font=3.5,
                           Lab.back="white",
                           nn=50,
                           frms=300,
                           n.frms=600,
                           do.animation=TRUE,
                           speed=100)  #lower speed value is faster
    
    
    #2. Effort and management measures timeline
    AMM.eff.managmtnt=function(dat,MGMT,SEP,Lab.font,Lab.back,Effort.lab)
    {
      d=dat%>%
        mutate(finyear=as.numeric(substr(FINYEAR,1,4)))
      Man=MGMT%>%
              dplyr::select(Event,finyear,Category)%>%
              mutate(id=1:nrow(MGMT))%>%
              filter(finyear>=min(d$finyear))
      d=d%>%
        full_join(Man,by='finyear')%>%
        arrange(finyear)%>%
        filter(finyear>=1975)%>%
        mutate(Total=ifelse(is.na(Total),na.approx(Total),Total),
               Category=ifelse(is.na(Category),"",Category),
               Event=ifelse(is.na(Event),"",Event))
      
      this=which(!is.na(d$FINYEAR))
      this=this[length(this)]
      
      if(is.na(d$FINYEAR[nrow(d)])) 
      {
        d$Total[(this+1):nrow(d)]=d$Total[this]
      }
      
      colors=c(Data='red',Closure='steelblue',Other='forestgreen','None'='orange')
      p=d%>%
        mutate(Category=ifelse(Category=='','None',Category),
               Category=factor(Category))%>%
        ggplot(aes(finyear, Total,label=Event,color = Category)) +
        geom_line(colour="orange",size=1.25) + 
        #geom_point() +
        geom_label_repel(box.padding=SEP,min.segment.length=0,hjust = 0,size = Lab.font,fill=Lab.back,max.overlaps=30) + 
        geom_line(colour="orange",size=1.25, alpha = 0.3) + 
        #geom_point(alpha = 0.4)+
        ylab(Effort.lab) + xlab("Financial year")+
        theme(legend.position = "none",
              axis.text=element_text(size=12),
              axis.title=element_text(size=14))+
        scale_color_manual(values = colors)
      return(print(p))
    }
    
    #TDGDFL
    set.seed(11)
    p=AMM.eff.managmtnt(dat=Total.effort.days.monthly,
                        MGMT=Management,
                        SEP=.15,
                        Lab.font=3,
                        Lab.back="white",
                        Effort.lab="Effort (1000 km gn d)")
    
    ggsave(handl_OneDrive('Management/Sharks/Timeline management measures/Effort.management_TDGDLF.tiff'), 
           width = 8,height = 5, dpi = 300, compression = "lzw")
    
    
    #NSF
    #note: NSF only have a fishery code since 1988-89
    Total.effort.NSF=read.csv(paste(HNDL,"Annual.total.eff_NSF.csv",sep=""),stringsAsFactors=F)
    NSF_add.0.effort=seq(2009,as.numeric(substr(Current.yr,1,4)))
    NSF_add.0.effort=paste(NSF_add.0.effort,substr(NSF_add.0.effort+1,3,4),sep='-')
    NSF_add.0.effort=data.frame(FINYEAR=NSF_add.0.effort,
                                Hook.days=0,
                                Hook.hours=0)%>%
                      filter(!FINYEAR%in%unique(Total.effort.NSF$FINYEAR))
    
    Total.effort.NSF=rbind(Total.effort.NSF,NSF_add.0.effort)%>%
            arrange(FINYEAR)%>%
            mutate(Total=Hook.days)
    set.seed(666)
    p.N=AMM.eff.managmtnt(dat=Total.effort.NSF,
                          MGMT=Management.north,
                          SEP=.15,
                          Lab.font=3,
                          Lab.back="white",
                          Effort.lab="Effort (1000 hook d)")
    ggsave(handl_OneDrive('Management/Sharks/Timeline management measures/Effort.management_NSF.tiff'),
           width = 8,height = 5, dpi = 300, compression = "lzw")
    
    
    #Two fisheries combined
    ggarrange(p+scale_y_continuous(trans='log2'), p.N+scale_y_continuous(trans='log2'), ncol = 1, nrow = 2)
    ggsave(handl_OneDrive('Management/Sharks/Timeline management measures/Effort.management_TDGDLF & NSF.tiff'),
           width = 8,height = 8, dpi = 300, compression = "lzw")
    
  }
   
  
  #3. Summary landed weight
  Summary.landwt=DAT%>%
        mutate(Taxa=case_when(
                      SPECIES%in%Elasmo.species ~ "Elasmo",
                      SPECIES%in%Scalefish.species ~ "Scalefish"))%>%
        group_by(Taxa)%>%
        summarise(Catch.tons=round(sum(LANDWT)/1000,0))
  fun.Tab1.SoFaR.LANDWT=function(SP,id,Tab.livewt,Add)
  {
    Dat=DAT  
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))
    {
      Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    }
    Dat=subset(Dat,!(is.na(Spec.tab1)))
    Dat$LANDWT=Dat$LANDWT/1000
    TABLA=aggregate(LANDWT~Spec.tab1+Bioregion,data=Dat,sum,na.rm=T)
    TABLA2=aggregate(LANDWT~Spec.tab1+zone,data=Dat,sum,na.rm=T)
    wide <- reshape(TABLA,v.names="LANDWT",timevar="Bioregion",idvar="Spec.tab1",direction="wide")
    colnames(wide)=c("Species","South","West")
    wide=wide[,match(c("Species","West","South"),names(wide))]
    wide2 <- reshape(TABLA2,v.names="LANDWT",timevar="zone",idvar="Spec.tab1",direction="wide")
    colnames(wide2)=c("Species","West.dem","Zone1","Zone2")
    wide=merge(wide,wide2,by="Species")
    wide[,2:ncol(wide)]=round(wide[,2:ncol(wide)])
    wide$Total=rowSums(wide[,2:3],na.rm=T)
    
    Matched=match(id,wide$Species)
    id.0=id[which(is.na(Matched))]
    Matched=Matched[!is.na(Matched)]
    wide=wide[Matched,]
    
    all=colSums(wide[,-1],na.rm=T)
    nms=names(all)
    all=matrix(all,nrow=1)
    colnames(all)=nms
    all=cbind(Species=Add,as.data.frame(all))
    wide=rbind(wide,all)
    wide[,3:ncol(wide)]=round(wide[,3:ncol(wide)],1)
    wide[is.na(wide)] = ""
    wide[wide==0] = "<0.1"
    wide=wide%>%rename(Total.land=Total)
    
    Tab.livewt=Tab.livewt%>%
      mutate_all(funs(ifelse(.=='<0.1', 0, .)))%>%
      mutate_at(vars(JASDGLF.Zn1,JASDGLF.Zn2,WCDGDLF,Total),
                as.numeric)%>%
      mutate_if(is.numeric, round, 0)
    Tab.livewt[is.na(Tab.livewt)] = ""
    Tab.livewt[Tab.livewt==0] = "<0.1"
    
    wide=Tab.livewt%>%
      left_join(wide,by=c('Name' = 'Species'))%>%
      mutate(JASDGLF.Zn1=paste(JASDGLF.Zn1," (",Zone1,")",sep=''),
             JASDGLF.Zn2=paste(JASDGLF.Zn2," (",Zone2,")",sep=''),
             WCDGDLF=paste(WCDGDLF," (",West.dem,")",sep=''),
             Total=paste(Total," (",Total.land,")",sep=''))%>%
      dplyr::select(Name,JASDGLF.Zn1,JASDGLF.Zn2,WCDGDLF,Total)
    
    wide=wide%>%mutate_all(function(x) ifelse(x==' ()','',x))
    return(wide)
  }
  AAM.table.Elasmos=fun.Tab1.SoFaR.LANDWT(SP=Spec.tab.1.elasmo,
                                          id=id.elasmo,
                                          Tab.livewt=Table1.SoFar.Elasmos,
                                          Add='Total Elasmobranchs')
  AAM.table.Scalies=fun.Tab1.SoFaR.LANDWT(SP=Spec.tab.1.scalies,
                                          id=id.scale,
                                          Tab.livewt=Table1.SoFar.Scalies%>%
                                            filter(Name!='Demersal scalefish suite component'),
                                          Add='Total Scalefish')
  dummy=MainFeatures.doc
  colnames(dummy)[1:2]=c("First.col","Second.col")
  dummy=dummy%>%
    mutate(Catch.tons=str_remove(Catch.tons, " t"))
  dummy=dummy%>%
    mutate(Catch.tons=
             case_when(SPECIES=="Total sharks and rays"~
                         paste(Catch.tons," (",
                               Summary.landwt[match('Elasmo',Summary.landwt$Taxa),'Catch.tons'],
                               ")",sep=''),
                       SPECIES=="Total scalefish"~
                         paste(Catch.tons," (",
                               Summary.landwt[match('Scalefish',Summary.landwt$Taxa),'Catch.tons'],
                               ")",sep=''),
                       SPECIES=="Gummy shark"~AAM.table.Elasmos[match('Gummy',AAM.table.Elasmos$Name),'Total'],
                       SPECIES=="Dusky shark"~AAM.table.Elasmos[match('Dusky_whaler',AAM.table.Elasmos$Name),'Total'],
                       SPECIES=="Sandbar shark"~AAM.table.Elasmos[match('Sandbar',AAM.table.Elasmos$Name),'Total'],
                       SPECIES=="Whiskery shark"~AAM.table.Elasmos[match('Whiskery',AAM.table.Elasmos$Name),'Total'],
                       TRUE ~ Catch.tons))
  write.table(dummy,"AMM_1.1.MainFeatures.LANDWT.csv",sep = ",",row.names = F)
  write.table(AAM.table.Elasmos,"AMM_Table1.SoFar.Elasmos.csv",sep = ",",row.names = F)
  write.table(AAM.table.Scalies,"AMM_Table1.SoFar.Scalies.csv",sep = ",",row.names = F)
}
