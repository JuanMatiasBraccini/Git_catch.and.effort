# Script for calculating temporal changes in spatial distribution index of catch and effort

library(tidyverse)
library(data.table)
library(ggpmisc)

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

hndl.in=handl_OneDrive("Analyses/Data_outs")
hndl.out=handl_OneDrive("Analyses/Catch and effort/Spatial distribution index")


# Data section ------------------------------------------------------------
setwd(hndl.in)
Data.monthly.GN=fread("Data.monthly.GN.csv",data.table=FALSE)

Effort.daily=fread("Effort.daily.csv",data.table=FALSE)   #MISSING: need to aggregate and combine the efforts
Effort.monthly=fread("Effort.monthly.csv",data.table=FALSE)

All.species.names=read.csv(handl_OneDrive("Data/Species.code.csv"))%>%
  filter(!COMMON_NAME=='Common blacktip shark')


# Parameters section ------------------------------------------------------------
Scenario=data.frame(Scenario=c('S1','S2'),
                    Percentile=c(.95,.5))

Protected.sp=c(8001,10003) #protected species (greynurse and white shark) 


# Data manipulations ------------------------------------------------------------
Data.monthly.GN$BLOCKX=as.integer(substr(Data.monthly.GN$BLOCKX,1,4))
Effort.monthly$BLOCKX=as.integer(substr(Effort.monthly$BLOCKX,1,4))
Effort.daily$blockx=as.integer(substr(Effort.daily$blockx,1,4))

Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~VESSEL+BLOCKX+FINYEAR+MONTH,
                             data=subset(Effort.monthly,NETLEN.c>100 & METHOD=="GN"),max,na.rm=T)
names(Eff.monthly.hour.c)=tolower(names(Eff.monthly.hour.c))

Eff.daily.hour.c.daily=aggregate(Km.Gillnet.Hours.c~Same.return.SNo+vessel+blockx+finyear+month,
                            data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)%>%
                      group_by(vessel,blockx,finyear,month)%>%
                      summarise(Km.Gillnet.Hours.c=sum(Km.Gillnet.Hours.c))%>%
                      data.frame
names(Eff.daily.hour.c.daily)=tolower(names(Eff.daily.hour.c.daily))
Eff=rbind(Eff.monthly.hour.c,Eff.daily.hour.c.daily)


# Calculate Distribution index for sharks and rays ------------------------------------------------------------
First.occurrence=setDT(Data.monthly.GN)%>%
                    mutate(Yr=as.numeric(substr(FINYEAR,1,4)))%>%
                    filter(!SPECIES%in%c(22999,31000))%>%
                    filter(SPECIES<5e4)

First.occurrence=First.occurrence[, .SD[which.min(Yr)], by = SPECIES]%>%
                  dplyr::select(FINYEAR,Yr,SPECIES)%>%
                  arrange(SPECIES)%>%
                  rename(FINYEAR.first=FINYEAR,
                         Yr.first=Yr)

Distrib.inx=Data.monthly.GN%>%
              filter(!SPECIES%in%c(22999,31000))%>%
              filter(SPECIES<5e4)%>%
              mutate(Yr=as.numeric(substr(FINYEAR,1,4)))%>%
              group_by(SPECIES,Yr,BLOCKX)%>%
              summarise(Catch=sum(LIVEWT.c))%>%
              left_join(First.occurrence,by='SPECIES')%>%
              mutate(N.blk.fished=length(unique(Data.monthly.GN$BLOCKX)))


N.blks.fished=length(unique(Distrib.inx$BLOCKX))

function.dist.in=function(d,Percentil)
{
  d1=d%>%
    ungroup()%>%
    data.frame%>%
    arrange(Yr,-Catch)%>%
    group_by(Yr)%>%
    mutate(cs = cumsum(Catch),
           Percentl=quantile(cs,Percentil),
           delta=Percentl-cs)%>%
    filter(delta>=0)%>%
    group_by(Yr)%>%
    tally()%>%
    mutate(Dist.indx=n/N.blks.fished,
           Species=unique(d$SPECIES))
  return(d1)
}

sp.list=Distrib.inx%>%   #species explaining 95% of catch
          group_by(SPECIES)%>%
          summarise(Tot=sum(Catch)/1000)%>%
          arrange(-Tot)%>%
          mutate(cs = cumsum(Tot),
                 TOT=sum(Tot),
                 Cum.per=cs/TOT,
                 delta=0.95-Cum.per)%>%
          filter(delta>-0)%>%
          pull(SPECIES)
input.scen=vector('list',nrow(Scenario))
names(input.scen)=Scenario$Scenario
input.list=vector('list',length(sp.list))
for(i in 1:nrow(Scenario))
{
  for(s in 1:length(sp.list))
  {
    input.list[[s]]=function.dist.in(d=Distrib.inx%>%filter(SPECIES==sp.list[s]),
                                     Percentil=Scenario$Percentile[i])
  }
  dummy=do.call(rbind,input.list)%>%
    mutate(Scenario=Scenario$Percentile[i])
  input.scen[[i]]=dummy
}



#Effort percentiles
N.blks.fished=length(unique(Eff$blockx))
function.dist.in=function(d,Percentil)    
{
  d1=d%>%
    group_by(blockx,finyear)%>%
    summarise(km.gillnet.hours.c=sum(km.gillnet.hours.c))%>%
    ungroup()%>%
    mutate(Yr=as.numeric(substr(finyear,1,4)))%>%
    data.frame%>%
    arrange(Yr,-km.gillnet.hours.c)%>%
    group_by(Yr)%>%
    mutate(cs = cumsum(km.gillnet.hours.c),
           Percentl=quantile(cs,Percentil),
           delta=Percentl-cs)%>%
    filter(delta>=0)%>%
    group_by(Yr)%>%
    tally()%>%
    mutate(Dist.indx=n/N.blks.fished)
  return(d1)
}
input.scen.effort=vector('list',nrow(Scenario))
names(input.scen.effort)=Scenario$Scenario

for(i in 1:nrow(Scenario))
{
  dummy=function.dist.in(d=Eff,
                   Percentil=Scenario$Percentile[i])%>%
           mutate(Scenario=Scenario$Percentile[i])
  input.scen.effort[[i]]=dummy
}


#Plotting
my.formula <- y ~ x   #fit linear model for trend

fn.plt=function(d,d1)
{
  d%>%
    left_join(All.species.names%>%distinct(CAAB_code,COMMON_NAME),by=c("Species"="CAAB_code"))%>%
    ggplot(aes(x=Yr,y=Dist.indx))+
    geom_point(colour="steelblue",size=2.5)+
    geom_smooth(method = "lm",se=T, color="red", formula = y ~ x)+
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE,size = 5)+
    facet_wrap(~COMMON_NAME)+
    xlab("Financial year")+ylab('Distribution index')+
    theme_bw() +
    theme(strip.text = element_text(size = 14),
          legend.text = element_text(size = 16),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16))+
    ylim(0,1)+
    geom_line(data=d1,aes(Yr,Dist.indx),colour=adjustcolor( "darkgreen", alpha.f = 0.6),size=1.5)
}

#S1
fn.plt(d=input.scen$S1,d1=input.scen.effort$S1)
ggsave(paste(hndl.out,"S1.tiff",sep='/'),width = 12,height = 8,compression = "lzw")

#S2
fn.plt(d=input.scen$S2,d1=input.scen.effort$S2)
ggsave(paste(hndl.out,"S2.tiff",sep='/'),width = 12,height = 8,compression = "lzw")
