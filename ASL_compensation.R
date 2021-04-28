# --- Script for calculating distance and time travel in relation to closures

library(data.table)
library(tidyverse)
library(dplyr)
library(geosphere)
library(ggpmisc)
library(mgcv)
library(emmeans) 

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


# Data section ------------------------------------------------------------
setwd(handl_OneDrive('Analyses/Data_outs'))
Data.daily.GN=fread("Data.daily.GN.csv",data.table=FALSE)
BlOCK_10=read.csv(handl_OneDrive("Data/Mapping/Blocks_10NM.csv"))

#ports
Ports=read.csv(handl_OneDrive('Data/Ports.csv'))

n.trips=1  #minimum number of trips per year
N.min=1  #minimum number of months per year
N.yrs=1  #minimum number of years with data

# Procedure section ------------------------------------------------------------

Data=Data.daily.GN%>%
    group_by(Same.return.SNo,PORT,VESSEL,zone,block10,day,MONTH,YEAR.c,FINYEAR)%>%
    summarise(Catch=sum(LIVEWT.c))%>%
    left_join(BlOCK_10,by=c('block10'='BlockNo'))%>%
    left_join(Ports,by="PORT")

Data=Data%>%
      separate(Same.return.SNo,c("SNo", "DSNo","TSNo"))%>%
      filter(!TSNo=='A2012019072700')%>%
      arrange(YEAR.c,TSNo,MONTH,day)

Tab=Data%>%
    distinct(TSNo,.keep_all = T)%>%
    group_by(VESSEL,YEAR.c)%>%
    summarise(N=n())%>%
    spread(VESSEL,N,fill=0)%>%
    data.frame%>%
  dplyr::select(-YEAR.c)

Tab[Tab<n.trips]=0    
Tab[Tab>=n.trips]=1   
vesls=names(which(colSums(Tab)>N.yrs))
vesls=gsub(".", " ", vesls, fixed=TRUE)

Distance.travld=vector('list',length(vesls))
for(v in 1:length(vesls))
{
  dummy=Data%>%
    filter(VESSEL==vesls[v])
  dummy$port.arrive=lead(dummy$PORT,1)
  dummy$Port_Longitude.arrive=lead(dummy$Port_Longitude,1)
  dummy$Port_Latitude.arrive=lead(dummy$Port_Latitude,1)
   
  trips=sort(unique(dummy$TSNo))
  Store=vector('list',length(trips))
  for(i in 1:length(trips))
  {
    a=dummy%>%
      filter(TSNo==trips[i])%>%
      arrange(MONTH,day)
    n=nrow(a)
    if(is.na(a$Port_Longitude.arrive[n]))
    {
      if(n>1)
      {
        a$Port_Longitude.arrive[n]=a$Port_Longitude.arrive[n-1]
        a$Port_Latitude.arrive[n]=a$Port_Latitude.arrive[n-1]
      }
      if(n==1)
      {
        a$Port_Longitude.arrive[n]=a$Port_Longitude[n]
        a$Port_Latitude.arrive[n]=a$Port_Latitude[n]
      }
    }
    Mat=with(a,cbind(c(Port_Longitude[1],Longitude_Centroid,Port_Longitude.arrive[n]),
                     c(Port_Latitude[1],Latitude_Centroid,Port_Latitude.arrive[n])))
    Distance=sum(distHaversine(Mat[-1,],Mat[-nrow(Mat),]))/1000 #in km
    out=a[1,c('TSNo','PORT','VESSEL','zone','MONTH','YEAR.c','FINYEAR')]
    out$km.travelled=Distance
    Store[[i]]=out
    rm(out)
  }
  Distance.travld[[v]]=do.call(rbind,Store)
}

Distance.travld=do.call(rbind,Distance.travld)


Dist.trvl.sumery=Distance.travld%>%
                group_by(PORT,MONTH,YEAR.c)%>%
                summarise(mean=mean(km.travelled),
                          sd=sd(km.travelled),
                          n=length(km.travelled),
                          se=sd/sqrt(n))%>%
                mutate(Yr.mnz=YEAR.c+MONTH/13)

#Output mean distance and error for full data set
setwd(handl_OneDrive('Analyses\\Catch and effort\\ASL.closure_compensation'))

#export data
write.csv(Distance.travld,'Distance.travld.csv',row.names = F)

Tab=Dist.trvl.sumery%>%
  group_by(PORT,YEAR.c)%>%
  summarise(n=sum(n))%>%
  spread(YEAR.c,n,fill=0)
nms=Tab$PORT
Tab=as.matrix(Tab[,-1])
rownames(Tab)=nms
Tab[Tab<N.min]=0
Tab[Tab>=N.min]=1
Sum=rowSums(Tab)

this.port=names(which(Sum>N.yrs))

my.formula <- y ~ x

Dist.trvl.sumery%>%
  filter(PORT%in%this.port)%>%
  ggplot(aes(Yr.mnz,mean,color=as.factor(MONTH)))+
    expand_limits(y=0)+ 
  # geom_smooth(method = "lm", formula = my.formula)+
  #  stat_poly_eq(formula = my.formula, 
  #               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #               parse = TRUE)+
  geom_point(size=.9) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(vars(PORT), scales = "free")+  labs(color = "Month")+
  theme(axis.text.x=element_text(size=8))+
  ylab("Mean distance travelled (km)")+
  xlab("Year")
ggsave('Mean distance travelled.tiff', width = 10,height = 6, dpi = 300, compression = "lzw")


#GAM
do.gam=FALSE
if(do.gam)
{
  d=Distance.travld%>%
    filter(PORT%in%this.port)%>%
    mutate(ln.km.travelled=log(km.travelled),
           PORT=as.factor(PORT),
           VESSEL=as.factor(VESSEL))
  d$YEAR.c=factor(d$YEAR.c,levels=sort(unique(d$YEAR.c)))
  
  mod<-gam(ln.km.travelled~PORT+YEAR.c+s(MONTH,k=12,bs='cc')+s(VESSEL,bs='re'),data=d,method="REML")
  gam.check(mod)
  plot(mod)
  #predict port and year
  lsm=summary(emmeans(mod, c('PORT','YEAR.c'), type="link"))%>%
    mutate(response=exp(emmean)*exp(SE^2/2),
           lower.CL=exp(lower.CL)*exp(SE^2/2),
           upper.CL=exp(upper.CL)*exp(SE^2/2))
  lsm%>%
    ggplot(aes(YEAR.c,response))+
    geom_point(size=.9) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                  position=position_dodge(0.05))+
    expand_limits(y=0)+ 
    facet_wrap(vars(PORT), scales = "free")+
    theme(axis.text.x=element_text(size=8))+
    ylab("Mean distance travelled (km)")+
    xlab("Year")
  ggsave('GAM_Mean distance travelled.tiff', width = 8,height = 6, dpi = 300, compression = "lzw")
  
  
}
