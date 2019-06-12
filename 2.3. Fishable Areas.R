
#Fishable Areas for the main commercial species


library(RODBC)  			#include ODBC library for importing Acccess data
library(PBSmapping)    	#needed to obtain maps
par.default=par()



setwd("C:/Matias/Data/Catch and Effort")  # working directory




###################################################
#--- DATA SECTION ---#
###################################################



#Daily records( 2006-2011)
channel <- odbcConnectExcel("2009_10_(Jan_2011)") 
Data.daily.2006.07<- sqlFetch(channel,"2006-07")
Data.daily.2007.08<- sqlFetch(channel,"2007-08")
Data.daily.2008.09<- sqlFetch(channel,"2008-09")
Data.daily.2009.10<- sqlFetch(channel,"2009-10")
close(channel) 

channel <- odbcConnectExcel("2012_Sharkfigs")
Data.daily.2010.11<- sqlFetch(channel,"DATA")
close(channel)


#2013 report
channel <- odbcConnectAccess("2011_12/updated/Shark1213.mdb")
Data.daily.2011.12<- sqlFetch(channel,"shark", colnames = F)
close(channel)
Data.daily.2011.12=subset(Data.daily.2011.12,finyear=="2011-12")

#2014 report
channel <- odbcConnectAccess("2012_13/Shark2013.mdb")
Data.daily.2012.13<- sqlFetch(channel,"sharklog", colnames = F)
close(channel)
Data.daily.2012.13=subset(Data.daily.2012.13,finyear=="2012-13")




#bathymetry data
#    bathymetry data downloaded from http://topex.ucsd.edu/cgi-bin/get_data.cgi (Topography option)
Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi")
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)




#_____Procedure section____

#2006-2010
Data.daily=rbind(Data.daily.2006.07,Data.daily.2007.08,Data.daily.2008.09,Data.daily.2009.10)
Data.daily.1=Data.daily
Data.daily$nfish=Data.daily$orgnfish

This=c("finyear","LatDeg","LongDeg","year","depthMax","method","blockx","block10","species","sname1","nfish","livewt")

Data.daily=Data.daily[,match(This,names(Data.daily))]

#2010-2011
Data.daily.2010.11=Data.daily.2010.11[,match(This,names(Data.daily.2010.11))]

#2011-2012
Data.daily.2011.12=Data.daily.2011.12[,match(This,names(Data.daily.2011.12))]

#2012-2013
Data.daily.2012.13=Data.daily.2012.13[,match(This,names(Data.daily.2012.13))]

#Merge daily data sets
Data.daily=rbind(Data.daily,Data.daily.2010.11)           
Data.daily=rbind(Data.daily,Data.daily.2011.12)
Data.daily=rbind(Data.daily,Data.daily.2012.13)

Data.daily$LatDeg=-Data.daily$LatDeg

#Select Species of interest
Species=c(17003,17001,18003,18007)
SPECIES=c("Whiskery","Gummy","Dusky","Sandbar")
Data.daily=subset(Data.daily,species%in%Species)


#Block of interest (from GLM model)
Block.list=list(whis=c("2813","2814","2913","2914","3014","3015","3114","3115","3128","3214",
                       "3215","3224","3225","3226","3227","3228","3314","3315","3320","3321",
                       "3322","3323","3324","3325","3326","3327","3328","3414","3415","3416",
                       "3418","3419","3420","3421","3422","3423","3424","3425","3426","3514",
                       "3515","3516","3517","3518","3519"),
                gum=c("3128","3224","3225","3226","3227","3228","3229","3320","3321","3322",
                      "3323","3324","3325","3326","3327","3328","3417","3418","3419","3420",
                      "3421","3422","3423","3424","3428","3516","3517","3518","3519"),
                dus=c("2813","2814","2913","2914","3014","3015","3114","3115","3214","3215",
                      "3314","3315","3320","3414","3415","3416","3418","3419","3420","3515",
                      "3516","3517","3518","3519"),
                san=c("2612","2613","2712","2713","2714","2813","2814","2913","2914","3014",
                      "3015","3114","3115","3214","3215","3314","3315","3414","3415","3416",
                      "3418","3515","3516","3517","3518"))


Data.daily$blockx=with(Data.daily,ifelse(blockx%in%c(96021),25120,     #Shark Bay
                                         ifelse(blockx%in%c(96022,96023),26131,
                                                ifelse(blockx%in%c(97011),27132,                        #Abrolhos
                                                       ifelse(blockx%in%c(97012,97013),28132,       
                                                              ifelse(blockx%in%c(97014,97015),29132,
                                                                     ifelse(blockx%in%c(96010),33151,                        #Geographe Bay
                                                                            ifelse(blockx%in%c(96000),32150,                        #Cockburn sound
                                                                                   ifelse(blockx%in%c(96030),35181,blockx)))))))))         # King George sound

Data.daily$blockx=as.numeric(substr(Data.daily$blockx,1,4))


setwd("C:/Matias/Analyses/Catch and effort/Outputs/Fishable Area")

#Histogram by species
fun.hist=function(Sp,blk,Spname)
{
  dat=subset(Data.daily,species==Sp & blockx%in%blk & depthMax<200 & LatDeg>(-38))
  hist(dat$depthMax,main=Spname,ylab="",xlab="",col="grey")
  box()
  LAT=range(dat$LatDeg)
  LONG=range(dat$LongDeg)
  return(list(LAT=LAT,LONG=LONG))
}
LIMS=vector('list',4)
tiff("Depth_distribution.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(2,2),mar=c(1.5,3.7,1.5,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(1.9,.7,0))
for ( i in 1:length(Species)) LIMS[[i]]=fun.hist(Species[i],Block.list[[i]],SPECIES[i])
mtext("Depth (m)",side=1,line=1,font=1,las=0,cex=1.5,outer=T)
mtext("Frequency",side=2,line=-1.3,font=1,las=0,cex=1.5,outer=T)
dev.off()


Max.depth=c(200,100,200,200)   #max depth set for all species

#Calculate proportion of area within max depth

Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
xbat=sort(unique(Bathymetry$V1))
ybat=sort(unique(Bathymetry$V2))

#reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
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
couleurs=rev(gray(seq(0.2,0.9,length=numInt)))
numberLab=5
colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))

a=112:129
b=seq(-37,South.WA.lat[2],length.out=length(a))



fn.map=function(XLIM,YLIM,BLK,MAXZ)
{
  dat=subset(Data.daily,species==Sp & blockx%in%BLK & depthMax<200 & LatDeg>(-38))
  dat=dat[!duplicated(dat$blockx),]
  YLIM1=YLIM
  YLIM[1]=YLIM[1]-1
  plotMap(worldLLhigh, xlim=XLIM,ylim=YLIM,plt = c(.1, 1, 0.075, 1),
          col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)

  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-MAXZ),
          levels =-MAXZ,labcex=1.5,lty = 1,col=1:4,add=T)
#          nlevels = 3,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
#1 degree grid
axis(side = 1, at =XLIM[1]:XLIM[2], labels = XLIM[1]:XLIM[2], tcl = 50,lty=2,col="grey20")
axis(side = 2, at = YLIM1[1]:YLIM1[2], labels = YLIM1[1]:YLIM1[2],tcl =50,las=2,lty=2,col="grey20")

#10 mins grid
axis(side = 1, at =seq(XLIM[1],XLIM[2],0.2), labels = F, tcl = 50,lty=3,col="grey60")
axis(side = 2, at = seq(YLIM[1],YLIM[2],0.2), labels =F,tcl =50,lty=3,col="grey60")

dat$LongDeg=with(dat,ifelse(blockx%in%c("2813","2913"),113,LongDeg))
text(dat$LongDeg+0.5,dat$LatDeg-0.5,dat$blockx,cex=0.95)
box()
}

for ( i in 1:length(Species))
{
  tiff(paste(SPECIES[i],"fishable_area.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.map(LIMS[[i]]$LONG,LIMS[[i]]$LAT,Block.list[[i]],Max.depth[i])
  dev.off()
}
