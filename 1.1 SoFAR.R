#------ G 3. STATE OF FISHERIES REPORT ------

library(dplyr)
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/SoFaR.figs.R")

#note:  TDGDLF is defined as METHOD= "GN" or "LL" and non-estuaries and south of 26 S. 
#       'Other fisheries' are extracted from 'Other.fishery.catch' (this includes 'Data.monthly'  
#         for which gear is not 'GN' or 'LL' and 'Daily logbooks' from other fisheries that report 
#         sharks (Pilbara trawl, WCDF, mackerel, Kimberley barramundi, Shark Bay prawn, etc) 

#       table 81.d for monthly records has no teleost catches 

#       Must run catch rate standardisation script before so it can bring in the index with the latest year

#     SOFAR effort is reported as km.gn.d. calculated with Use.Date to "YES" (as instructed by Rory). 
#     If Use.Date = "NO" then 'ID' is used in the aggregation rather than 'date', then effort
#      in recent years is larger due to the multiple shot per day practice of
#     vessels from zone 2. Rory instructed this is biased. See '4.Script for testing effect of DATE or ID 
#     in effort calculation'.  This is not an issue for km.gn.hr

#     Before executing SOFAR bit, run "2.CPUE standardisations.R" to construct latest standardised cpue

Current.yr="2017-18"

Percent.fin.of.livewt=0.03
TDGDLF.lat.range=c(-26,-40)

Ray.species=25000:31000
Elasmo.species=5001:31000
Gummy=17001;Dusky_whaler=c(18003,18001);Whiskery=17003;Sandbar=18007;Hammerheads=19000;
Spinner=18023;Wobbegongs=13000;Common_saw_shark=23002;School=17008
Order.Elas.Sp.SoFAR=data.frame(Species=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads",
                                         "Spinner","Wobbegongs","Rays","Common_saw_shark","School","Other_elasmobranchs",
                                         "Total Elasmobranchs"),
                               Species_or_taxon=c("Mustelus antarcticus","Carcharhinus obscurus","Furgaleus macki",
                                                  "Carcharhinus plumbeus","F.Sphyrnidae","Carcharhinus brevipinna",
                                                  "F. Orectolobidae","Batoidea","Pristiophorus cirratus","Galeorhinus galeus","",""))


Scalefish.species=180000:599001
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
Sandbar.status="Recovering"
Whiskery.status="Adequate"


#Catch range of indicator species (Gummy, whiskery, dusky, sandbar)
Catch.range.key.species=data.frame(SPECIES=c("Gummy","Whiskery","Bronzy.Dusky","Sandbar"),
                                   Min.catch=c(350,175,200,0),
                                   Max.catch=c(450,225,300,120))%>%
                        mutate(SPECIES=as.character(SPECIES))

#bring in data from 1.Manipulate data.R
setwd("C:/Matias/Analyses/Catch and effort/Data_outs")

Data.monthly=read.csv("Data.monthly.csv",stringsAsFactors = F)
Rec.fish.catch=read.csv("Rec.fish.catch.csv",stringsAsFactors = F)
Data.current.Sofar=read.csv("Data.current.Sofar.csv",stringsAsFactors = F)
PRICES=read.csv("PRICES.csv",stringsAsFactors = F)
Total.effort.days.monthly=read.csv("Annual.total.eff.days.csv",stringsAsFactors = F)
Total.effort.hours.monthly=read.csv("Annual.total.eff.hours.csv",stringsAsFactors = F)
Total.effort.zone.hours.monthly=read.csv("Annual.zone.eff.hours.csv",stringsAsFactors = F)
Total.effort.zone.days.monthly=read.csv("Annual.zone.eff.days.csv",stringsAsFactors = F)
TEPS.current=read.csv("TEPS.current.csv",stringsAsFactors = F)
Suite=read.csv("suite.csv",stringsAsFactors = F)$x
TEPS.pre.current=read.csv("TEPS.pre.current.csv",stringsAsFactors = F)
Results.pre.2013=read.csv("Results.pre.2013.csv",stringsAsFactors = F)

  
handle.Sofar=paste("C:/Matias/Analyses/Catch and effort/State of fisheries/",Current.yr,sep="")
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
CAESS.shark.fishery.codes=c("SGL1","SGL2","SGL","WCGL","C127","CL02")
Other.fishery.catch=subset(Data.monthly,SPECIES<=31000 & !FisheryCode%in%CAESS.shark.fishery.codes)
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
Rec.Groups=Rec.Catch.stats=NULL
REC=subset(Rec.fish.catch,FinYear==Current.yr)  
#use latest estimate if no estimates available for current year
if(nrow(REC)==0)
{
  REC.yrs=sort(as.character(unique(Rec.fish.catch$FinYear)))
  REC=subset(Rec.fish.catch,FinYear==REC.yrs[length(REC.yrs)])  
}
REC.yr=as.character(unique(REC$FinYear))
REC.yr=paste(substr(REC.yr,1,4),"/",substr(REC.yr,6,7),sep="")
Tot.Ktch=sum(REC$Kept.Number)+sum(REC$Rel.Number)
Yr.Rec.Ktch=as.character(unique(REC$FinYear))
Tot.Kept=sum(REC$Kept.Number)
Rec.Catch.stats=data.frame(Group="All sharks and rays",Total.Catch=Tot.Ktch,Total.Retained=Tot.Kept)                           
REC=subset(REC,Bioregion%in%c("South Coast","West Coast"))
REC$Group=with(REC,
               ifelse(Common.Name%in%c("Bronze Whaler","Sandbar Shark","Tiger Shark","Whaler Sharks"),"Whalers",
                      ifelse(Common.Name%in%c("Gummy Sharks","School Shark","Whiskery Shark"),"Hounds",
                             ifelse(Common.Name%in%c("Wobbegong"),"wobbegongs",
                                    ifelse(Common.Name%in%c("Hammerhead Sharks"),"hammerheads ","other")))))
Rec.Groups=aggregate(Kept.Number~Group,REC,sum)
Rec.Groups$Percent=round(Rec.Groups$Kept.Number*100/sum(Rec.Groups$Kept.Number))
Rec.Groups$Bioregion="WC and SC"
Rec.Catch.stats=rbind(Rec.Catch.stats,data.frame(Group="All sharks and rays.SC&WC",
                                                 Total.Catch=sum(REC$Caught.Number),Total.Retained=sum(REC$Kept.Number)))
Rec.fish.catch.WC=subset(REC,Bioregion=="West Coast")
Rec.Catch.stats=rbind(Rec.Catch.stats,data.frame(Group="All sharks and rays.WC",
                                                 Total.Catch=sum(Rec.fish.catch.WC$Caught.Number),Total.Retained=sum(Rec.fish.catch.WC$Kept.Number)))
Rec.Groups.WC=aggregate(Kept.Number~Group,Rec.fish.catch.WC,sum)
Rec.Groups.WC$Percent=round(Rec.Groups.WC$Kept.Number*100/sum(Rec.Groups.WC$Kept.Number))
Rec.Groups.WC$Bioregion="WC"
Rec.Groups=rbind(Rec.Groups,Rec.Groups.WC)
REC=sum(REC$Kept.Number)
Avg.wt=5        #average weight of a shark                                                   
Rec.fi.Tot.catch=data.frame(SPECIES=paste("Recreational catch ","(",REC.yr,")",sep=""), Catch.tons=REC*Avg.wt/1000)
Re.Percnt.Tot.Com= Rec.fi.Tot.catch$Catch.tons*100/MainFeatures$Catch.tons[which(MainFeatures$SPECIES=="Total sharks and rays")]  
Rec.fi.Tot.catch$Catch.tons=ifelse(Re.Percnt.Tot.Com<5,"< 5% of commercial catch",
                                   ifelse(Re.Percnt.Tot.Com>=5 &Re.Percnt.Tot.Com<10,"< 10% of commercial catch",
                                          Re.Percnt.Tot.Com))

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
PRICES=PRICES%>%mutate(dolar.per.kg=PRICES[,match(paste("uv",substr(Current.yr,3,4),
                                        substr(Current.yr,6,7),sep=""),names(PRICES))],
                       dolar.per.kg=as.numeric(gsub("\\$", "", dolar.per.kg)))%>%
                select(-uv1718)

Fin.Weight=subset(DAT,SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,RSCommonName,zone,LIVEWT.c))
Fin.Weight=subset(Fin.Weight,!(SPECIES%in%c(13000,20000,31000,26999))) #remove no fin species (wobbies, spurdog,rays except shovelnose)
Fin.Weight$dolar.per.kg.fin=PRICES[match(22998,PRICES$SPECIES),match('dolar.per.kg',names(PRICES))]
Fin.Weight$fin.weight=Percent.fin.of.livewt*Fin.Weight$LIVEWT.c 
Fin.Weight$Total.price=Fin.Weight$dolar.per.kg.fin*Fin.Weight$fin.weight
FINS.value=aggregate(Total.price~zone,Fin.Weight,sum,na.rm=T)
FINS.value.species.zone=aggregate(cbind(Total.price,fin.weight)~SNAME+zone,Fin.Weight,sum,na.rm=T)



#catch value 
get.yr.price="dolar.per.kg"
PRICES1=PRICES[,match(c("SPECIES",get.yr.price),names(PRICES))]
names(PRICES1)=c("SPECIES","ABARE.dolar.per.kg")
PRICES1=merge(DAT[,match(c("SPECIES","SNAME","zone","LIVEWT.c"),names(DAT))],PRICES1,by="SPECIES",all.x=T)

# PRICES1=PRICES[,match(c("SPECIES","Weighted.Average.Price.kg","Abare.unit.."),names(PRICES))]
# names(PRICES1)=c("SPECIES","DOF.dolar.per.kg","ABARE.dolar.per.kg")
# PRICES1=merge(DAT[,match(c("SPECIES","SNAME","zone","LIVEWT.c"),names(DAT))],PRICES1,by="SPECIES",all.x=T)

PRICES1=subset(PRICES1,!is.na(LIVEWT.c))

No.price=PRICES1[is.na(PRICES1$ABARE.dolar.per.kg),]
No.price=No.price[!duplicated(No.price$SPECIES),]

Default.price=1

PRICES1$ABARE.dolar.per.kg=with(PRICES1,ifelse(is.na(ABARE.dolar.per.kg),Default.price,ABARE.dolar.per.kg))  #add default if no price data
PRICES1$Total.price=PRICES1$LIVEWT.c*PRICES1$ABARE.dolar.per.kg

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
id=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads","Spinner","Wobbegongs","Rays","Common_saw_shark","School",
     "Other_elasmobranchs")
Matched=match(id,Table1.SoFar.Elasmos$Species)
id.0=id[which(is.na(Matched))]
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

id=c("Blue_morwong","Blue_groper","West_Australian_dhufish","Pink_snapper","Boarfishes","Samsonfish","Redfishes","Mulloway","Sweetlips",
     "Baldchin_groper","Other_scalefish")
Table1.SoFar.Scalies=Table1.SoFar.Scalies[match(id,Table1.SoFar.Scalies$Species),]

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

fun.Tab2.SoFaR=function(Dat)
{
  Dat$finyear=as.character(Dat$finyear)
  Dat$Status=as.character(Dat$Status)
  Dat$Status=with(Dat,ifelse(Status=="a","A",
                             ifelse(Status=="d","D",Status)))
  Dat$Ali.Ded.yr=with(Dat,paste(finyear,Status))
  
  TABLA=aggregate(Number~SpeciesCode+Ali.Ded.yr,data=Dat,sum,na.rm=T)
  Species.id=Dat[,match(c("SpeciesCode","DataEntryName"),names(Dat))]
  Species.id=Species.id[!(duplicated(Species.id$SpeciesCode)),]
  
  wide <- reshape(TABLA,v.names="Number",timevar="Ali.Ded.yr",idvar="SpeciesCode",direction="wide")
  wide=wide[order(wide$SpeciesCode),]
  wide=merge(Species.id,wide,by="SpeciesCode",all.x=T)
  wide=wide[-match("SpeciesCode",names(wide))]
  names(wide)=c("Species",paste("alive.",unique(Dat$finyear),sep=""),paste("dead.",unique(Dat$finyear),sep=""))
  return(wide)
}

TEPS.WC.SC=subset(TEPS.current,fishery%in%c("SGL","WCGL") & finyear==Current.yr)

#Fix typo. change whales to whaler sharks
TEPS.WC.SC$SpeciesCode=with(TEPS.WC.SC,ifelse(DataEntryName=="DUSKY WHALER OVERSIZE",18003,SpeciesCode))

#review comments in case reported in comments and not in logbook columns
if(!file.exists("scan.teps.csv"))   
{
  # look at comments and see if comments match what is recorded in number; 
  #     if not, maybe it is a new interaction record.
  # Remember that comments are duplicated so only use unique record. 
  Scan.current=TEPS.WC.SC[,match(c("DailySheetNumber","Status","Number","DataEntryName","Comments"),names(TEPS.WC.SC))]
  Scan.current=subset(Scan.current,!is.na(Comments))                                                                                                 
  Scan.current=Scan.current[order(Scan.current$DailySheetNumber),]
  write.csv(Scan.current,"scan.teps.csv",row.names=F) 
  cat(paste("scan.teps.csv FILE LOCATED IN",getwd()))
}

Sighting.only=which(grepl("sighting", tolower(TEPS.WC.SC$Comments)))
TEPS.WC.SC.without.sightings=TEPS.WC.SC
if(length(Sighting.only)>0) TEPS.WC.SC.without.sightings=TEPS.WC.SC[-Sighting.only,]
TEP.table=fun.Tab2.SoFaR(TEPS.WC.SC.without.sightings)
TEP.table$Species=with(TEP.table,ifelse(Species=="WHITE POINTER","WHITE SHARK",Species))  

#remove current in case overwrote historic
II=match(paste(c("alive.","dead."),paste(substr(Current.yr,1,4),substr(Current.yr,6,7),sep="."),sep=""),colnames(TEPS.pre.current))
if(sum(II,na.rm=T)>0) TEPS.pre.current=TEPS.pre.current[,-II]

id=which(!(TEPS.pre.current$Species%in%TEP.table$Species))
Add.this=as.character(TEPS.pre.current$Species[id])
Add.this=data.frame(Species=Add.this,Alive=NA,Dead=NA)
names(Add.this)[2:3]=c(names(TEP.table)[2],names(TEP.table)[3])
TEP.table=rbind(TEP.table,Add.this)
TEP.table=TEP.table[order(TEP.table$Species),]
TEP.table=subset(TEP.table,!(Species=="DUSKY WHALER OVERSIZE"))  #not reporting oversized whalers
TEPS=merge(TEPS.pre.current,TEP.table,by="Species")
TEPS[is.na(TEPS)] = ""


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
        paste("Recreational catch ","(",REC.yr,")",sep=""))
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
                          select(SPECIES,Catch.tons)%>%
                          left_join(Catch.range.key.species,by="SPECIES")%>%
                          mutate(Acceptable=ifelse(Catch.tons>=Min.catch & Catch.tons<=Max.catch,"Acceptable",
                                                   "Unacceptable"))
 
write.csv(Ktch_comm_accept,"Ktch_comm_accept.csv",row.names=F)



#Export tables as .csv
if(!is.null(Rec.Catch.stats))write.csv(Rec.Catch.stats,paste("1.2.Rec.Catch.stats.",Yr.Rec.Ktch,".csv",sep=""),row.names=F)
if(!is.null(Rec.Catch.stats))write.csv(Rec.Groups,paste("1.3.Rec.Groups.SC&WC.",Yr.Rec.Ktch,".csv",sep=""),row.names=F)
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
  write.csv(TEPS,"C:/Matias/Data/Catch and Effort/Historic/Historic.TEPS.res.csv",row.names = F)  
}


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
                 Estuary=="NO")
    
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
    
    fun.fig.SoFar(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
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
fn.figs2.3.SoFaR(Elasmo.species,TDGDLF.lat.range[1],TDGDLF.lat.range[2],100,500)
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
               LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  
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
fn.figs2.3.SoFaR.scalies(Scalefish.species,TDGDLF.lat.range[1],TDGDLF.lat.range[2],20,100)
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
HNDL="C:/Matias/Analyses/Catch and effort/Data_outs/"
#Whiskery
whis.mon=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
whis.daily=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#gummy
gum.mon=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
gum.daily=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#dusky
dus.mon=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
dus.daily=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#sandbar
san.mon=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
san.daily=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

#plot cpues
Plot.cpue=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar)    
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
    ymax = max(cpuedata$UP.CI)
    Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
    plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
    CL="black"
    points(Yrs, cpuedata$Mean, "o", pch=16, lty=2, col=CL,cex=CxS)
    arrows(x0=Yrs, y0=cpuedata$LOW.CI, 
           x1=Yrs, y1=cpuedata$UP.CI, 
           code=3, angle=90, length=0.05, col=CL)
  }
}

jpeg(file="Figure All.cpues.jpeg",width = 2200, height = 2400,units = "px", res = 300)
Pred.creep=list(whis.mon,gum.mon,dus.mon,san.mon)
Pred.daily.creep=list(whis.daily,gum.daily,dus.daily,san.daily)
SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
for(s in 1:length(Pred.creep))
{
  #Monthly
  Mon.dat=Pred.creep[[s]]
  LgND="NO"
  if(s==1)LgND="YES"
  Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="Finyear")
  if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
  
  #Daily
  Daily.dat=Pred.daily.creep[[s]]
  LgND="NO"
  Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="Finyear")
  if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
  mtext(SPECIES.vec[s],4,line=1,las=3,cex=1.5)
}
mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
mtext("CPUE (kg/ km gillnet hour)",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
dev.off()







