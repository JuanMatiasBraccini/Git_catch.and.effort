#                 SCRIPT FOR DERIVING STANDARDISED CPUE FROM OBSERVER DATA ON TDGDLF    #

#notes:

rm(list=ls(all=TRUE))

options(stringsAsFactors = FALSE)
library(rlang)
library(tidyverse)
library(doParallel)
library(zoo)
library(abind)
library(stringr)
library(Hmisc)
library(mgcv)
library(zigam)
library(pscl)
library(MASS)
library(stringr)
library(yarrr)

# 1. Data ---------------------------------------------------------
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


# Observers data
User="Matias"
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))

# Species names
All.species.names=read.csv(handl_OneDrive("Data/Species_names_shark.only.csv")) #for catch

setwd(handl_OneDrive('Analyses/Catch and effort/Observer_TDGDLF'))

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Nominal_cpue_functions.R"))

# Parameters section ---------------------------------------------------------

# Criteria to define species to analyse
Min.yrs1=3
Min.obs.per.yr=10

#Standardisation criteria
Min.yr.blks=2             #at least this number of observations per year-block of positive catches
core.per=90               #core area defined as 90% of catch
use.core.area=FALSE
Min.boat= 5                 #use boats with at least these number of records
Min.yrs=Min.obs.per.yr      #at least 10 positive records per year
Min.yrs.data=Min.yrs1            #do standardisation for species with > 3 years of data
this.var=c('SHEET_NO','Effort','date','Month','Finyear','BOAT','zone','BLOCK','Mid.Lat','Mid.Long','MESH_SIZE')


# Manipulate species names ---------------------------------------------------------
All.species.names=All.species.names%>%
  mutate(Name=capitalize(tolower(Name)),
         Name=case_when(Name=='Port jackson shark'~'Port Jackson shark',
                        TRUE~Name))


# Manipulate observer data  ---------------------------------------------------------
Res.ves=c("HAM","HOU","NAT","FLIN","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")
Dat_obs=DATA %>%
  filter(Mid.Lat<=(-26) & !BOAT%in%Res.ves & Taxa%in%c('Elasmobranch') & !COMMON_NAME=='WHALE')  %>% 
  dplyr::select(c(SHEET_NO,date,Month,year,BOAT,SKIPPER,zone,BLOCK,SOAK.TIME,MESH_SIZE,
                  MESH_DROP,NET_LENGTH,Mid.Lat,Mid.Long,Method,SPECIES,Taxa,
                  COMMON_NAME,SCIENTIFIC_NAME,CAES_Code,TL,FL,Disc.width,Number))%>%
  filter(Method=="GN" & !is.na(NET_LENGTH) & !is.na(SOAK.TIME))%>%
  mutate(COMMON_NAME=ifelse(SPECIES=='PD','Spurdogs',COMMON_NAME),
         SCIENTIFIC_NAME=ifelse(SPECIES=='PD','Squalus spp.',SCIENTIFIC_NAME),
         CAES_Code=ifelse(SPECIES=='PD',20000,CAES_Code),
         SPECIES=ifelse(SPECIES=='PD','SD',SPECIES))%>%
  left_join(All.species.names%>%rename(Species=SPECIES),by=c('SPECIES'='SP'))%>%
  mutate(Name=ifelse(is.na(Name),COMMON_NAME,Name),
         Scien.nm=ifelse(is.na(Scien.nm),SCIENTIFIC_NAME,Scien.nm))%>%
  filter(NET_LENGTH>=0.2 & !is.na(MESH_SIZE))%>%
  mutate(BOAT=ifelse(BOAT=="TRACEY LEA","E35",BOAT),
         Finyear=ifelse(Month>6,paste(year,substr(year+1,3,4),sep='-'),
                 ifelse(Month<=6,paste(year-1,substr(year,3,4),sep='-'),
                        NA)))


  # Group Angel sharks
Dat_obs=Dat_obs%>%
  mutate(COMMON_NAME=case_when(SPECIES%in%c("AO","AU")~"Angel Sharks (general)", TRUE~COMMON_NAME),
         SCIENTIFIC_NAME=case_when(SPECIES%in%c("AO","AU")~"Family Squatinidae", TRUE~SCIENTIFIC_NAME),
         CAES_Code=case_when(SPECIES%in%c("AO","AU")~24900, TRUE~CAES_Code),
         Species=as.double(Species),
         Species=case_when(SPECIES%in%c("AO","AU")~24900, TRUE~Species),
         Name=case_when(SPECIES%in%c("AO","AU")~"Angel sharks", TRUE~Name),
         Scien.nm=case_when(SPECIES%in%c("AO","AU")~"Squatinidae", TRUE~Scien.nm),
         SPECIES=case_when(SPECIES%in%c("AO","AU")~"AA",TRUE~SPECIES))



# Set FL to disc width for stingrays for computational purposes----------------------------------------
Stingrays=35000:40000
Dat_obs=Dat_obs%>%
  mutate(FL=ifelse(CAES_Code%in%Stingrays,Disc.width,FL),
         TL=ifelse(CAES_Code%in%Stingrays,NA,TL))

# Fill in missing FL info  ---------------------------------------------------------
#1. Set to NA FL or TL records less than size at birth
size.birth=30 #FL, in cm
Dat_obs$FL=with(Dat_obs,ifelse(FL<size.birth,NA,FL))
Dat_obs$TL=with(Dat_obs,ifelse(TL<(size.birth/.85),NA,TL))


#2. Derive FL as a proportion of TL
Dat_obs$FL=with(Dat_obs,ifelse(is.na(FL),TL*.875,FL))

#3. Define species to analyse
This.sp=Dat_obs%>%
  group_by(SPECIES,Finyear)%>%
  tally()%>%
  mutate(n=ifelse(n<Min.obs.per.yr,0,1))%>%
  group_by(SPECIES)%>%
  summarise(n=sum(n))%>%
  filter(n>Min.yrs1)%>%
  pull(SPECIES)




# Size frequency of analysed species  ---------------------------------------------------------
Dat_obs%>%
  filter(SPECIES%in%This.sp)%>%
  filter(!Name=='Stingrays')%>%  #Stringrays set to average size as not measured, no point displaying
  mutate(Name=ifelse(Name=="Eagle ray","Southern eagle ray",Name))%>%
  ggplot( aes(x=FL, color=Name, fill=Name)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  theme(legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size =9),
        axis.title=element_text(size=16)) +
  xlab("Size (cm)") +
  ylab("Frequency") +
  facet_wrap(~Name, scales = "free")
ggsave('Figure_Size.frequency.tiff',width = 10,height = 10,compression = "lzw")

# Arrange data to one row per shot  ---------------------------------------------------------
Dat=Dat_obs%>%
  filter(SPECIES%in%This.sp)%>%
  group_by(SHEET_NO,date,Month,Finyear,BOAT,zone,BLOCK,Mid.Lat,Mid.Long,MESH_SIZE,SPECIES)%>%
  summarise(N=sum(Number),
            SOAK.TIME=max(SOAK.TIME),
            NET_LENGTH=max(NET_LENGTH))%>%
  data.frame

Ktch=Dat%>%
  dplyr::select(-c(SOAK.TIME,NET_LENGTH))%>%
  spread(SPECIES,N,fill = 0)
Effort=Dat%>%
        dplyr::select(SHEET_NO,date,BLOCK,SOAK.TIME,NET_LENGTH)%>%
        mutate(Effort=SOAK.TIME*NET_LENGTH,
               dupli=paste(SHEET_NO,date,BLOCK))%>%
  distinct(dupli,.keep_all = T)%>%
  dplyr::select(SHEET_NO,date,BLOCK,Effort)

Ktch=left_join(Ktch,Effort,by=c('SHEET_NO','date','BLOCK'))


# Nominal cpue  ---------------------------------------------------------
#note: this has all records
Ktch$season=as.numeric(substr(Ktch$Finyear,1,4))

Nominal.cpue=vector('list',length(This.sp))
names(Nominal.cpue)=This.sp
pdf("Figure_Nominal cpues.pdf")
for(s in 1:length(This.sp))
{
  x=This.sp[s]
  nm=unique(Dat_obs%>%filter(SPECIES==x)%>%pull(Name))
  Nominal.cpue[[s]]=CalcMeanCPUE(cpuedata = Ktch, catch.column=x, effort.column="Effort",
                                 plot.title = nm, cpue.units = "Number per km gn hours", 
                                 draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
}
dev.off()


# Standardised cpue  ---------------------------------------------------------
#Check which species have enough records
this.sp.enough=This.sp
Store.dat=vector('list',length(This.sp))
names(Store.dat)=This.sp
pdf("Figure_Core areas.pdf")
for(s in 1:length(This.sp))
{
  d=Ktch[,c(this.var,This.sp[s])]%>%
    mutate(yr.blk=paste(Finyear,BLOCK))%>%
    rename(Catch=!!sym(This.sp[s]))%>%
    filter(!is.na(BOAT))
  
  plot(d$Mid.Long,d$Mid.Lat,ylim=c(-36,-26),xlim=c(113,129),
       pch=21,cex=.8,col="grey80",bg="white",xlab="Longitude",ylab='Latitude',
       main=unique(Dat_obs%>%filter(SPECIES==This.sp[s])%>%pull(Name)))
  
  #Core area
  if(use.core.area)
  {
    Ag.blk=d%>%
      group_by(BLOCK)%>%
      summarise(tot=sum(Catch))%>%
      arrange(-tot)%>%
      mutate(cumsum=cumsum(tot),
             percentile=cumsum/sum(tot))%>%
      filter(percentile<=core.per/100)%>%
      pull(BLOCK)
    d=d%>%filter(BLOCK%in%Ag.blk)
  }

  
  #Min.yr.blks
  Tab=d%>%filter(Catch>0)%>%
    group_by(yr.blk)%>%
    tally()%>%
    filter(n>=Min.yr.blks)
  d=d%>%filter(yr.blk%in%Tab$yr.blk)
  
  #boat records
  Tab=d%>%
    group_by(BOAT)%>%
    tally()%>%
    filter(n>=Min.boat)
  d=d%>%filter(BOAT%in%Tab$BOAT)
  
  Tab=table(d$BOAT,d$Finyear)
  Tab[Tab<=2]=0
  Tab[Tab>2]=1
  Tab=rowSums(Tab)>2
  d=d%>%filter(BOAT%in%names(which(Tab==TRUE)))
  
  #Min.pos.years
  Tab=d%>%filter(Catch>0)%>%
    group_by(Finyear)%>%
    tally()%>%
    filter(n>=Min.yrs)
  d=d%>%filter(Finyear%in%Tab$Finyear)
  
  if(nrow(d)>0)
  {
    points(d$Mid.Long,d$Mid.Lat,pch=21,bg="steelblue",cex=3*d$Catch/max(d$Catch))
    legend("topright",c("All shots","Positive shot (prop to catch)"),
           pch=c(21,21),cex=1.25,col=c("grey70","black"),pt.bg=c("white","steelblue"),bty='n')
  }
    
  if(nrow(d)==0)
  {
    this.sp.enough[s]=NA
    legend("topright","All shots",pch=21,cex=1.25,col="grey70",bg="white",bty='n')
  }
    
  
  #Remove vessels for which coefficients cannot be estimated
  if(This.sp[s]%in%c("BW","WH","HZ","PN")) d=subset(d,!BOAT=='E67')
  if(This.sp[s]%in%c("GM")) d=subset(d,!BOAT=='E7')
  if(This.sp[s]%in%c("WW","WD")) d=subset(d,!BOAT=='F517')
  
  
  #Manipulations for standardisation
  if(nrow(d)>0)
  {
    d=d%>%
      mutate(BOAT=as.factor(BOAT),
             BLOCK=as.factor(BLOCK),
             Finyear=as.factor(Finyear),
             MESH=ifelse(!MESH_SIZE%in%c(6,6.5,7),'other',MESH_SIZE),
             MESH=as.factor(MESH),
             Mn=factor(Month,levels=1:12),
             log.Effort=log(Effort))
    Store.dat[[s]]=d
  }
  
  if(nrow(d)>0 & length(table(d$Finyear))<=Min.yrs.data) this.sp.enough[s]=NA
  rm(d)
}
dev.off()

#select suitable species
this.sp.enough=subset(this.sp.enough,!is.na(this.sp.enough))
this.sp.enough=subset(this.sp.enough,!this.sp.enough%in%c("AA","SR","WB","SH"))  #also remove species complex
Store.dat=Store.dat[match(this.sp.enough,names(Store.dat))]

#Select  best model
Stand.cpue=vector('list',length(this.sp.enough))
names(Stand.cpue)=this.sp.enough

formula.gam=formula(Catch ~ Finyear + BLOCK + s(BOAT,bs="re") + s(Month, bs = "cc") + offset(log.Effort))
formula.gam=replicate(length(this.sp.enough),formula.gam,FALSE)
names(formula.gam)=this.sp.enough
formula.gam$GN=formula.gam$SD=formula(Catch ~ Finyear + BLOCK + s(BOAT,bs="re") + offset(log.Effort))
#formula.gam$SD=formula(Catch ~ Finyear + BLOCK + offset(log.Effort))
#formula.gam$WC=formula(Catch ~ Finyear + s(Month, bs = "cc") + offset(log.Effort))

fn.stand=function(d,FORMULA)
{
  pois=gam(FORMULA,family=poisson(),data=d,method="REML")
  NB=gam(FORMULA,family=nb(),data=d,method="REML")
  return(list(data=d, pois=pois,NB=NB))
}
system.time({for(s in 1:length(this.sp.enough))
{
  Stand.cpue[[s]]=fn.stand(d=Store.dat[[s]],FORMULA=formula.gam[[s]])
}})

do.glm=FALSE
if(do.glm)
{
  formula.glm=formula("Catch~ Finyear + BLOCK + BOAT + Mn + offset(log.Effort)")
  formula.glm=replicate(length(this.sp.enough),formula.glm,FALSE)
  names(formula.glm)=this.sp.enough
  formula.glm$PJ=formula("Catch~ Finyear + BLOCK + BOAT + Mn + MESH + offset(log.Effort)")
  formula.glm$PN=formula.glm$WD=formula("Catch~ Finyear + BLOCK + BOAT  + offset(log.Effort)")
  formula.glm$WC=formula("Catch ~ Finyear + Mn + offset(log.Effort)")
  formula.glm$WH=formula("Catch ~ Finyear + BLOCK +  Mn + offset(log.Effort)")
  formula.glm$LG=formula("Catch ~ Finyear + Mn + offset(log.Effort)")
  formula.glm$ER=formula("Catch ~ Finyear + BOAT + Mn + offset(log.Effort)")
  formula.glm$TK=formula("Catch ~ Finyear +  BOAT + Mn + offset(log.Effort)")
  formula.glm$GM=formula("Catch ~ Finyear + BLOCK+   Mn + offset(log.Effort)")
  formula.glm$WW=formula("Catch ~ Finyear +  BOAT + Mn + offset(log.Effort)")
  formula.glm$WD=formula("Catch ~ Finyear +  BOAT + Mn + offset(log.Effort)")
  formula.glm.zero.zip=formula("Catch~ Finyear + BLOCK + BOAT + Mn + offset(log.Effort) |
                                       Finyear + BLOCK + BOAT + Mn + offset(log.Effort)")
  formula.glm.zero.zip=replicate(length(this.sp.enough),formula.glm.zero.zip,FALSE)
  names(formula.glm.zero.zip)=this.sp.enough
  formula.glm.zero.znb=formula.glm.zero.zip
  formula.glm.zero.zip$BW=formula("Catch ~ Finyear + BLOCK + Mn + MESH + offset(log.Effort) |
                      Finyear + BOAT + Mn + MESH + offset(log.Effort)")   #count part then zero part
  formula.glm.zero.zip$TK=formula("Catch ~ Finyear + BLOCK+ BOAT+ Mn+ offset(log.Effort) |
                      Finyear + BLOCK+ BOAT + Mn + MESH + offset(log.Effort)")
  formula.glm.zero.zip$GM=formula("Catch ~ Finyear  + BOAT + Mn + MESH + offset(log.Effort) |
                      Finyear  + BOAT + offset(log.Effort)")
  formula.glm.zero.zip$WH=formula("Catch ~ Finyear +  BOAT + Mn + MESH + offset(log.Effort) |
                       Finyear + BLOCK + Mn + offset(log.Effort)")
  formula.glm.zero.zip$HZ=formula("Catch ~ Finyear + BOAT + Mn + MESH + offset(log.Effort) |
                      Finyear + BOAT + Mn  + MESH + offset(log.Effort)")
  formula.glm.zero.zip$LG=formula("Catch ~ Finyear + BLOCK + BOAT + Mn + MESH + offset(log.Effort) |
                      Finyear + BLOCK + BOAT  + Mn + offset(log.Effort)")
  formula.glm.zero.zip$WW=formula("Catch ~ Finyear + Mn + offset(log.Effort) |
                      Finyear + BOAT  + offset(log.Effort)")
  formula.glm.zero.zip$PN=formula("Catch ~ Finyear + BOAT + Mn + offset(log.Effort) |
                      Finyear + BOAT +  offset(log.Effort)")
  formula.glm.zero.zip$WD=formula("Catch ~ Finyear +  BOAT + Mn + offset(log.Effort) |
                      Finyear + BLOCK + offset(log.Effort)")
  formula.glm.zero.zip$PJ=formula("Catch ~ Finyear + BLOCK + BOAT + Mn + MESH + offset(log.Effort) |
                      Finyear + BLOCK + BOAT + Mn  + MESH + offset(log.Effort)")
  formula.glm.zero.zip$ER=formula("Catch ~ Finyear + BLOCK + BOAT + Mn  + MESH + offset(log.Effort) |
                      Finyear + BLOCK + BOAT  + MESH + offset(log.Effort)")
  formula.glm.zero.zip$WC=formula("Catch ~ Finyear + Mn + offset(log.Effort) |
                      Finyear  + offset(log.Effort)")
  
  formula.glm.zero.znb$BW=formula("Catch ~ Finyear +  offset(log.Effort) | 
                     Finyear +  offset(log.Effort)")
  formula.glm.zero.znb$TK=formula("Catch ~ Finyear + BLOCK+ BOAT+ Mn+ offset(log.Effort) | 
                     Finyear + BLOCK+ BOAT + Mn + offset(log.Effort)")
  formula.glm.zero.znb$GM=formula("Catch ~ Finyear  + BOAT + Mn + offset(log.Effort) |
                      Finyear  + BOAT + offset(log.Effort)")
  formula.glm.zero.znb$WH=formula("Catch ~ Finyear + BOAT  + Mn + MESH + offset(log.Effort) |
                      Finyear + offset(log.Effort)")
  formula.glm.zero.znb$HZ=formula("Catch ~ Finyear + BOAT  + Mn + MESH + offset(log.Effort) |
                      BLOCK  + MESH + offset(log.Effort)")
  formula.glm.zero.znb$LG=formula("Catch ~ Finyear + offset(log.Effort) |
                      Finyear + BOAT + offset(log.Effort)")
  formula.glm.zero.znb$WW=formula("Catch ~ Finyear  + offset(log.Effort) |
                      Finyear + BOAT  + offset(log.Effort)")
  formula.glm.zero.znb$PN=formula("Catch ~ Finyear + BOAT + Mn + offset(log.Effort) |
                      Finyear + BOAT +  offset(log.Effort)")
  formula.glm.zero.znb$WD=formula("Catch ~ Finyear +  BOAT + Mn + offset(log.Effort) |
                      Finyear + BLOCK + offset(log.Effort)")
  formula.glm.zero.znb$PJ=formula("Catch ~ Finyear + BLOCK + BOAT + Mn + MESH + offset(log.Effort) |
                      Finyear + BLOCK + BOAT + offset(log.Effort)")
  formula.glm.zero.znb$ER=formula("Catch ~ Finyear + BLOCK  + Mn  + offset(log.Effort) |
                      Finyear + BOAT  + offset(log.Effort)")
  formula.glm.zero.znb$WC=formula("Catch ~ Finyear + Mn + offset(log.Effort) |
                      Finyear  + offset(log.Effort)")
  
  fn.stand=function(d,FORMULA,FORMULA.zero.zip,FORMULA.zero.znb)
  {
    pois=glm(FORMULA,data=d,family=poisson)
    NB=glm.nb(FORMULA, data=d)
    ZIP=zeroinfl(FORMULA.zero.zip, data=d)
    ZINB=zeroinfl(FORMULA.zero.znb, dist = "negbin", data=d)
    
    #Fit.pois=gam(FORMULA,data=d,method = "REML",family=poisson)
    #Fit.NB=gam(FORMULA, data=d,method = "REML",family = nb)
    #Fit.ZIP=zipgam(lambda.formula=FORMULA,pi.formula=FORMULA,data=d)
    #Fit.ZINB=zinbgam(mu.formula=FORMULA,pi.formula=FORMULA, data=d)
    
    return(list(data=d, pois=pois,NB=NB,ZIP=ZIP,ZINB=ZINB))
    
  }
  system.time({for(s in 1:length(this.sp.enough))
  {
    
    Stand.cpue[[s]]=fn.stand(d=Store.dat[[s]],
                             FORMULA=formula.glm[[s]],
                             FORMULA.zero.zip=formula.glm.zero.zip[[match(this.sp.enough[s],names(formula.glm.zero.zip))]],
                             FORMULA.zero.znb=formula.glm.zero.znb[[match(this.sp.enough[s],names(formula.glm.zero.znb))]])
  }})
  
  Best.model=list(
    BW=list(formula=formula.glm$BW,
            error='NB'),
    TK=list(formula=formula.glm$TK,
            error='NB'),
    GM=list(formula=formula.glm$GM,
            error='NB'),
    WH=list(formula=formula.glm$WH,
            error='NB'),
    HZ=list(formula=formula.glm$HZ,
            error='NB'),
    LG=list(formula=formula.glm$LG,
            error='NB'),
    WW=list(formula=formula.glm$WW,
            error='NB'),
    PN=list(formula=formula.glm$PN,
            error='NB'),
    WD=list(formula=formula.glm$WD,
            error='NB'),
    PJ=list(formula=formula.glm$PJ,
            error='NB'),
    ER=list(formula=formula.glm$ER,
            error='NB'),
    WC=list(formula=formula.glm$WC,
            error='NB')
  )
  
}

AIC.tab=vector('list',length(this.sp.enough))
for(s in 1:length(this.sp.enough))
{
  pois=Stand.cpue[[s]]$pois
  NB=Stand.cpue[[s]]$NB
  AIC.tab[[s]]=data.frame(SP=names(Stand.cpue)[s],
                          Pois=pois$aic,
                          NB=NB$aic)
  if(do.glm)
  {
    ZIP=Stand.cpue[[s]]$ZIP
    ZINB=Stand.cpue[[s]]$ZINB
    AIC.tab[[s]]=data.frame(SP=names(Stand.cpue)[s],
                            Pois=pois$aic,
                            NB=NB$aic,
                            ZIP=AIC(ZIP),
                            ZNB=AIC(ZINB))
    #logLik(pois)  
    print(names(Stand.cpue)[s])
    vuong(ZIP,pois)
    vuong(ZIP,NB)
    vuong(ZIP,ZINB)
    vuong(ZINB,pois)
    vuong(ZINB,NB) 
    print("______________________________________________")
  }
  
}
AIC.tab=do.call(rbind,AIC.tab)
AIC.tab$Best.mod=colnames(AIC.tab)[apply(AIC.tab[,-1],1,which.min)+1]


# Standardise using best model
Best.model=vector('list',length=nrow(AIC.tab))
names(Best.model) <- AIC.tab$SP
for(s in 1:nrow(AIC.tab)) Best.model[[s]]=list(formula=formula.gam[[s]],error=AIC.tab$Best.mod[s])

fn.stand=function(d,FORMULA,error)
{
  if(error=="Pois") mod=gam(FORMULA,family=poisson(),data=d,method="REML")
  if(error=='NB')   mod=gam(FORMULA,family=nb(),data=d,method="REML")
  return(list(data=d, mod=mod))
}
if(do.glm)
{
  fn.stand=function(d,FORMULA,error)
  {
    if(error=="Pois") mod=glm(FORMULA,data=d,family=poisson)
    if(error=='NB') mod=glm.nb(FORMULA, data=d)
    if(error=='ZIP') mod=zeroinfl(FORMULA, data=d)
    if(error=='ZNB') mod=zeroinfl(FORMULA, dist = "negbin", data=d)
    
    return(list(data=d, mod=mod))
    
  }
}

system.time({for(s in 1:length(this.sp.enough))
{
  Stand.cpue[[s]]=fn.stand(d=Store.dat[[s]],
                           FORMULA=Best.model[[s]]$formula,
                           error=Best.model[[s]]$error)
}})

#Predict year effect
fn.factor=function(x)
{
  Tab=sort(table(x))
  return(factor(names(Tab[length(Tab)]),levels=levels(x)))
}
pred.fn.boot=function(MOD,Dta,response,Formula)
{
  this=c(labels(terms(Formula)),'log.Effort')
  this=unique(unlist(str_split(this, " +")))
  this=subset(this,this%in%colnames(Dta))
  this=this[-match(response,this)]
  new.dat=Dta%>%
          dplyr::select(all_of(this))%>%
          mutate_if(is.factor, fn.factor)%>%
          mutate_if(is.double,mean)%>%
          distinct(!!sym(this[1]),.keep_all = T)
  Pred.res=Dta%>%
            dplyr::select(all_of(response))%>%
            distinct(!!sym(response))%>%
            arrange(!!sym(response))
  new.dat=cbind(Pred.res,new.dat)

  PREDS=predict(MOD,new.dat, type='response')
  out=cbind(new.dat[response],Mean=PREDS)
  out=out%>%
    mutate(year=as.numeric(substr(Finyear,1,4)))
  
  return(out)
}
pred.fn=function(MOD,Dta,response,Formula)
{
  #this=c(labels(terms(Formula)),'log.Effort')
  this=all.vars(Formula)
  this=this[-match(response,this)][-1]
  new.dat=Dta%>%
    dplyr::select(all_of(this))%>%
    mutate_if(is.factor, fn.factor)%>%
    mutate_if(is.double,mean)%>%
    distinct(!!sym(this[1]),.keep_all = T)
  Pred.res=Dta%>%
    dplyr::select(all_of(response))%>%
    distinct(!!sym(response))%>%
    arrange(!!sym(response))
  new.dat=cbind(Pred.res,new.dat)
  
  PREDS=predict(MOD,new.dat, type='response',se.fit = T)
  out=cbind(new.dat[response],Mean=PREDS$fit,se=PREDS$se.fit)
  out=out%>%
    mutate(year=as.numeric(as.numeric(substr(Finyear,1,4))),
           CV=se/Mean*100,
           LOW.CI=Mean-1.96*se,
           UP.CI=Mean+1.96*se,
           LOW.CI=LOW.CI/mean(Mean),
           UP.CI=UP.CI/mean(Mean),
           Mean=Mean/mean(Mean))
  return(out)
}
PREDS=vector('list',length(this.sp.enough))
names(PREDS)=this.sp.enough
n.boot=1000
library(boot)
system.time({ for(s in 1:length(this.sp.enough))
{
  if(Best.model[[s]]$error%in%c("Pois","NB"))
  {
    PREDS[[s]]=pred.fn(MOD=Stand.cpue[[s]]$mod,
                            Dta=Stand.cpue[[s]]$data,
                            response="Finyear",
                            Formula=Best.model[[s]]$formula)
  }
  if(Best.model[[s]]$error%in%c("ZIP","ZNB"))
  {
    set.seed(10)
    prEds=vector('list',n.boot)
    for(n in 1:n.boot)
    {
      #re-fit model
      d=Store.dat[[s]]
      ii=sample(1:nrow(d),nrow(d),replace = T)
      
      tryCatch({
        mod=zeroinfl(Best.model[[s]]$formula, dist = "negbin",
                     data=d[ii,],
                     start = list(count=round(coef(Stand.cpue[[s]]$mod, "count"), 4),
                                  zero=round(coef(Stand.cpue[[s]]$mod, "zero"), 4)))
      }, error = function(e) 
      {
      })
      
      #preds
      if(exists("mod"))
      {
        prEds[[n]]=pred.fn.boot(MOD=mod,
                           Dta=d[ii,],
                           response="Finyear",
                           Formula=Best.model[[s]]$formula)

        remove(d,mod)  
      }
      
    }
    PREDS[[s]]=prEds
  }
  print(s)  

}})    #takes 0.025 per iteration


# Plot --------------------------------------------------------------------
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
AXIS1=seq(10*round(min(Dat_obs$year)/10),10*round(max(Dat_obs$year)/10),by=5)
AXIS2=seq(10*round(min(Dat_obs$year)/10),10*round(max(Dat_obs$year)/10),by=10)
         
fn.plt=function(d,nm,dat,size=.8)
{
  d$year=as.numeric(substr(d$Finyear,1,4))
  nm=unique(Dat_obs%>%filter(SPECIES==nm)%>%pull(Name))
  plot(d$year,d$Mean,ylim=c(0,max(d$UP.CI)),pch=19,cex=size,
       xlim=c(min(AXIS2),max(AXIS2)),ylab='',xlab='',xaxt='n')
  lines(d$year,d$Mean,lty=2,col=transparent("black",.5))
  segments(d$year,d$LOW.CI,d$year,d$UP.CI,lwd=size)
  mtext(nm,3,cex=.9)
  axis(1,AXIS1,F)
  axis(1,AXIS2,AXIS2)
  
  
  # Nominal CPUE
  dat=dat%>%
    mutate(year=as.numeric(substr(Finyear,1,4)))
  out1 = dat %>%
  group_by(year) %>%
    summarise(My = mean(Catch),
              Mx = mean(Effort),
              Sy = sd(Catch),
              Sx = sd(Effort),
              r = cor(Catch, Effort),
              n = length(Catch)) %>%
    as.data.frame
  out1$r[is.na(out1$r)] = 0
  out1 = out1 %>%
    mutate(mean=My/Mx,
           se =  sqrt(1/n*(My^2*Sx^2/(Mx^4) + Sy^2/(Mx^2) - 2*My*r*Sx*Sy/(Mx^3))),
           lowCL = mean - 1.96*se,
           uppCL = mean + 1.96*se) %>%
    as.data.frame%>%
    mutate(lowCL=lowCL/mean(mean),
           uppCL=uppCL/mean(mean),
           mean=mean/mean(mean),
           year=as.numeric(as.character(year)))
  with(out1,points(year-.25,mean,col="grey70",pch=19,cex=size))
  with(out1,lines(year,mean,lty=2,col=transparent("grey70",.6)))
  with(out1,segments(year-.25,lowCL,year-.25,uppCL,lwd=size,col="grey70"))

  
  
  #Arithmetic Mean CPUE
  out2 = dat %>%
    mutate(cpue=Catch/Effort)%>%
    group_by(year) %>%
    summarise(mean = mean(cpue),
              n = length(cpue),
              sd = sd(cpue)) %>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame %>%
    mutate(lowCL=lowCL/mean(mean),
           uppCL=uppCL/mean(mean),
           mean=mean/mean(mean),
           year=as.numeric(as.character(year)))
  with(out2,points(year+.25,mean,col="steelblue",pch=19,cex=size))
  with(out2,lines(year,mean,lty=2,col=transparent("steelblue",.6)))
  with(out2,segments(year+.25,lowCL,year+.25,uppCL,lwd=size,col="steelblue"))
  
}
tiff(file="Figure_Stand.CPUE.tiff",width = 2400, height = 2000,
     units = "px", res = 300, compression = "lzw")  
par(cex.axis=.85)
smart.par(n.plots=length(this.sp.enough),MAR=c(1,1.5,1.5,.5),OMA=c(2,2,.15,.1),MGP=c(.1, 0.5, 0))
for(s in 1:length(this.sp.enough)) fn.plt(d=PREDS[[s]],nm=names(PREDS)[s],dat=Stand.cpue[[s]]$data)
plot.new()
legend("bottomright",c("Standardised","Nominal","Arithmetic mean"),
       col=c("black","grey70","steelblue"),bty='n',pch=19,cex=1.25)
mtext("Relative cpue",2,outer=T,las=3,line=0.5)
mtext("Financial year",1,outer=T,line=0.5)
dev.off()

# b=Store.dat[[s]]%>%
#   mutate(year=factor(year))
# ggplot(data=b,aes(Mid.Long,Mid.Lat))+
#   geom_point() +
#   geom_point(data=subset(b,Catch>0),aes(Mid.Long,Mid.Lat,color=year))

# Export --------------------------------------------------------------------
hndl=handl_OneDrive("Analyses/Data_outs")
for(s in 1:length(this.sp.enough))
{
  nm=unique(Dat_obs%>%filter(SPECIES==names(PREDS)[s])%>%pull(Name))
  nm.original=nm
  if(nm%in%c("Cobbler wobbegong","Western wobbegong","Banded wobbegong")) nm='Wobbegongs'
  out=PREDS[[s]]%>%
    dplyr::select(-se)
  write.csv(out,paste(hndl,'/',nm,'/',nm.original,".CPUE_Observer_TDGDLF.csv",sep=""),row.names=F)
  rm(nm,out)
}

