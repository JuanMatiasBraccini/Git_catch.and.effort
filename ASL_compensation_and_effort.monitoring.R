# --- Script for calculating distance and time traveled in relation to closures and spatial overlap

library(data.table)
library(tidyverse)
library(dplyr)
library(geosphere)
library(ggpmisc)
library(mgcv)
library(emmeans) 
library(ggforce)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


# Data section ------------------------------------------------------------
setwd(handl_OneDrive('Analyses/Data_outs'))
Data.daily.GN=fread("Data.daily.GN.csv",data.table=FALSE)
Data.daily.LL=fread("Data.daily.LL.csv",data.table=FALSE)
BlOCK_10=read.csv(handl_OneDrive("Data/Mapping/Blocks_10NM.csv"))

#ports
Ports=read.csv(handl_OneDrive('Data/Ports.csv'))

n.trips=1  #minimum number of trips per year
N.min=1  #minimum number of months per year
N.yrs=1  #minimum number of years with data

# Distance traveled ------------------------------------------------------------

these.finyrs=sort(unique(Data.daily.GN$FINYEAR))
these.finyrs=these.finyrs[(grep('18',these.finyrs)[1]-1):length(these.finyrs)]

Data=Data.daily.GN%>%
    filter(FINYEAR%in%these.finyrs)%>%
    group_by(Same.return.SNo,PORT,VESSEL,zone,block10,LatFC,LongFC,day,MONTH,YEAR.c,FINYEAR)%>%
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

consider.port.arrive=FALSE   #do not include dummy Port.arrive as this is not reported so it's assumed to be the depature Port

system.time({for(v in 1:length(vesls))    #takes 8 secs per year
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
    if(consider.port.arrive)
    {
      Mat=with(a,cbind(c(Port_Longitude[1],Longitude_Centroid,Port_Longitude.arrive[n]),
                       c(Port_Latitude[1],Latitude_Centroid,Port_Latitude.arrive[n])))
    }
    if(!consider.port.arrive)
    {
      Mat=with(a,cbind(c(Port_Longitude[1],Longitude_Centroid),
                       c(Port_Latitude[1],Latitude_Centroid)))
    } 
    
    Distance=sum(distHaversine(Mat[-1,],Mat[-nrow(Mat),]))/1000 #in km
    out=a[1,c('TSNo','PORT','VESSEL','zone','MONTH','YEAR.c','FINYEAR')]
    out$km.travelled=Distance
    Store[[i]]=out
    rm(out)
  }
  Distance.travld[[v]]=do.call(rbind,Store)
}})

Distance.travld=do.call(rbind,Distance.travld)%>%
            mutate(YEAR.c=as.numeric(substr(FINYEAR,1,4)))

N.records=table(Distance.travld$PORT)
dis.ports=names(N.records[N.records>20])

Distance.travld=Distance.travld%>%
  filter(PORT%in%dis.ports)

#Output mean distance and error for full data set
setwd(handl_OneDrive('Analyses\\Catch and effort\\ASL.closure_compensation'))

Distance.travld%>%
  mutate(FINYEAR=substr(FINYEAR,1,4))%>%
  ggplot(aes(x=FINYEAR,y=km.travelled,color=PORT))+
  geom_violin(alpha=0.2)+
  geom_sina(alpha=0.5)+
  stat_summary(fun='mean', geom='point', size=3, col='black')+
  stat_summary(fun='median', geom='point', size=3, col='grey60')+
  facet_wrap(~PORT,scales='free')+
  theme(legend.position = 'none',
        axis.text.x=element_text(size=7))+
  ylab("Mean distance travelled (km)")+
  xlab("Financial year")+expand_limits(y=0)
ggsave('Distribution distance traveled.tiff', width = 10,height = 6, dpi = 300, compression = "lzw")



# Dist.trvl.sumery=Distance.travld%>%
#                 group_by(PORT,MONTH,YEAR.c)%>%
#                 summarise(mean=mean(km.travelled),
#                           sd=sd(km.travelled),
#                           n=length(km.travelled),
#                           se=sd/sqrt(n))%>%
#                 mutate(Yr.mnz=YEAR.c+MONTH/13)

Dist.trvl.sumery=Distance.travld%>%
  group_by(PORT,YEAR.c)%>%
  summarise(mean=mean(km.travelled),
            sd=sd(km.travelled),
            n=length(km.travelled),
            se=sd/sqrt(n))


#export data
write.csv(Distance.travld,'Distance.travld.csv',row.names = F)

# Tab=Dist.trvl.sumery%>%
#   group_by(PORT,YEAR.c)%>%
#   summarise(n=sum(n))%>%
#   spread(YEAR.c,n,fill=0)
# nms=Tab$PORT
# Tab=as.matrix(Tab[,-1])
# rownames(Tab)=nms
# Tab[Tab<N.min]=0
# Tab[Tab>=N.min]=1
# Sum=rowSums(Tab)
# 
# this.port=names(which(Sum>N.yrs))
# 
# my.formula <- y ~ x
# 
# Dist.trvl.sumery%>%
#   filter(PORT%in%this.port)%>%
#   ggplot(aes(YEAR.c,mean,color=PORT))+
#     expand_limits(y=0)+ 
#    # geom_smooth(method = "lm", formula = my.formula)+
#    #  stat_poly_eq(formula = my.formula, 
#    #               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#    #               parse = TRUE)+
#   geom_point(size=.9) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
#                 position=position_dodge(0.05))+
#   facet_wrap(vars(PORT), scales = "free")+
#   #labs(color = "Month")+
#   theme(legend.position = 'none',
#               axis.text.x=element_text(size=7))+
#   ylab("Mean distance travelled (km)")+
#   xlab("Financial year")
# ggsave('Mean distance traveled.tiff', width = 10,height = 6, dpi = 300, compression = "lzw")


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



# Effort overlap map ------------------------------------------------------
#notes:   Ideally, shot location should be derived from VMS as many shots have lat and long derived from block10
#          see 2. CPUE standardisations_spatial.protection

#1. Map showing shot location and ASL closure
library(rgdal)
paz=handl_OneDrive("Data/Mapping/Closures/")
ASL_Closures_SC=readOGR(paste(paz,"ASL_Closures/ASL_Closures.shp",sep=''),
                     layer="ASL_Closures") 
ASL_Closures_WC=readOGR(paste(paz,"ASL_Closures/Prohibition_on_Fishing_WCDGDLIMF_Order_2018_DPIRD_088.shp",sep=''),
                     layer="Prohibition_on_Fishing_WCDGDLIMF_Order_2018_DPIRD_088") 
WA_Commonwealth_Marine_Parks=readOGR(paste(paz,"WA_Commonwealth_Marine_Parks/WA_Commonwealth_Marine_Parks.shp",sep=''),
                                     layer="WA_Commonwealth_Marine_Parks") 
Shark_Fishery_Closures=readOGR(paste(paz,"Shark_Fishery_Closures/Shark_Fishery_Closures.shp",sep=''),
                               layer="Shark_Fishery_Closures") 

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))

South.WA.long=c(109,129)
South.WA.lat=c(-40,-9)
Xlim=South.WA.long
Ylim=South.WA.lat
col.shark.closure=rgb(.1,.7,.3,alpha=.6)
col.Marine.park.closure=rgb(1,.1,.1,alpha=.6)
col.ASL.closure=rgb(.1,.3,.6,alpha=.6)


Overlap.GN=Data.daily.GN%>%
  distinct(Same.return.SNo,FINYEAR,fishery,LatFC,LongFC,FINYEAR,YEAR.c)%>%
  filter(FINYEAR%in%these.finyrs)

Overlap.LL=Data.daily.LL%>%
  distinct(Same.return.SNo,FINYEAR,fishery,LatFC,LongFC,FINYEAR,YEAR.c)%>%
  filter(FINYEAR%in%these.finyrs)

Do.non.ASL=FALSE
fn.plt.map=function(Xlim,Ylim,d,non.ASL.closure.on=Do.non.ASL)
{
  plot(1,xlim=c(Xlim[1]*0.9995,Xlim[2]*0.99),ylim=Ylim,xlab="",ylab="",axes=F,main="")

  #closures
  if(non.ASL.closure.on) plot(WA_Commonwealth_Marine_Parks,add=T,col=col.Marine.park.closure)
  if(non.ASL.closure.on) plot(Shark_Fishery_Closures,add=T,col=col.shark.closure)
  plot(ASL_Closures_SC,add=T,col=col.ASL.closure)
  plot(ASL_Closures_WC,add=T,col=col.ASL.closure)
  
  #shots
  points(d$LongFC,d$LatFC,bg=rgb(.1,.1,.1,alpha=.1),col=rgb(.1,.1,.1,alpha=.1),pch=21,cex=.8)
  
  #coast
  polygon(WAcoast$Longitude,WAcoast$Latitude, col="burlywood3")

  box()
}

tiff("Map_overlap_GN.tiff",width=2400,height=2400,units="px",res=300,compression="lzw")
smart.par(n.plots=length(these.finyrs),MAR=c(1,1,1,1),OMA=c(3,3,1,.3),MGP=c(.04,.6,0))
par(las=1)
for(f in 1:length(these.finyrs))
{
  fn.plt.map(Xlim=c(113,129),Ylim=c(-35,-27),d=Overlap.GN%>%filter(FINYEAR==these.finyrs[f]))
  mtext(these.finyrs[f],3,line=-1)
  #axes
  axis(side = 1, seq(South.WA.long[1],South.WA.long[2],2), labels =F, tck = -.02)
  axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],2), labels = F, tck = -.02)
  
  axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
  axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),seq(round(Xlim[1]),round(Xlim[2]),2),cex.axis=1.25)
  
}
if(Do.non.ASL)
{
  legend('center',c("Fishery closures","Marine parks","ASL closures"),
         fill =c(col.shark.closure,col.Marine.park.closure,col.ASL.closure),
         bty='n',cex=1.2)
}else
{
  legend('center',c("ASL closures"),fill =c(col.ASL.closure),bty='n',cex=1.2)
}
mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1,font=1,las=0,cex=1.35,outer=T)
mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
dev.off()


tiff("Map_overlap_LL.tiff",width=2400,height=2400,units="px",res=300,compression="lzw")
smart.par(n.plots=length(these.finyrs),MAR=c(1,1,1,1),OMA=c(3,3,1,.3),MGP=c(.04,.6,0))
par(las=1)
for(f in 1:length(these.finyrs))
{
  fn.plt.map(Xlim=c(113,129),Ylim=c(-35,-27),d=Overlap.LL%>%filter(FINYEAR==these.finyrs[f]))
  mtext(these.finyrs[f],3,line=-1)
  #axes
  axis(side = 1, seq(South.WA.long[1],South.WA.long[2],2), labels =F, tck = -.02)
  axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],2), labels = F, tck = -.02)
  
  axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
  axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),seq(round(Xlim[1]),round(Xlim[2]),2),cex.axis=1.25)
  
}
if(Do.non.ASL)
{
  legend('center',c("Fishery closures","Marine parks","ASL closures"),
         fill =c(col.shark.closure,col.Marine.park.closure,col.ASL.closure),
         bty='n',cex=1.2)
}else
{
  legend('center',c("ASL closures"),fill =c(col.ASL.closure),bty='n',cex=1.2)
}
mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1,font=1,las=0,cex=1.35,outer=T)
mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
dev.off()



#2. Calculate intersection with ASL colonies
library(sf)
Intersection_SC = read_sf(paste(paz,"ASL_Closures/ASL_Closures.shp",sep=''))
Intersection_WC = read_sf(paste(paz,"ASL_Closures/Prohibition_on_Fishing_WCDGDLIMF_Order_2018_DPIRD_088.shp",sep=''))


Method_list=list(GN=Overlap.GN,LL=Overlap.LL)
ASL_list=list(SC=Intersection_SC,WC=Intersection_WC)

Intersect.fn=function(pnts,shp)
{
  pnts_sf <- st_as_sf(pnts, coords = c('y', 'x'), crs = st_crs(shp))
  pnts <- pnts_sf%>%
    mutate(intersection = as.integer(st_intersects(geometry, shp)),
           intersection=ifelse(is.na(intersection),'NO',"YES"))%>%
    data.frame
  return(pnts)
}

Store.intersec=vector('list',length(Method_list))
names(Store.intersec)=names(Method_list)
for(m in 1:length(Method_list))
{
  dummy=vector('list',length(ASL_list))
  names(dummy)=names(ASL_list)
  for(a in 1:length(ASL_list))
  {
    xx=Intersect.fn(pnts=Method_list[[m]]%>%rename(x=LatFC,y=LongFC),
                    shp=ASL_list[[a]])
    dummy[[a]]=xx%>%
      group_by(FINYEAR,intersection)%>%
      tally()%>%
      spread(intersection,n,fill=0)%>%
      ungroup()%>%
      mutate(Intersec.prop=YES/(NO+YES),
             Bioregion=names(ASL_list)[a])
  }
  Store.intersec[[m]]=do.call(rbind,dummy)%>%
    mutate(Method=names(Method_list)[m])
}

Store.intersec=do.call(rbind,Store.intersec)


Store.intersec%>%
  ggplot(aes(x=FINYEAR,y=Intersec.prop,fill=Method))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  facet_wrap(~Bioregion,ncol=1)+
  ylab("Proportion of fishing shots inside ASL closures")+
  ylim(0,1)+
  theme(strip.text.x = element_text(size = 14),
        axis.text=element_text(size=12),
        legend.position="top",
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 14),
        title=element_text(size=16))

ggsave("Proportion of shots inside ASL closures.tiff", width = 8,height = 10, dpi = 300, compression = "lzw")

