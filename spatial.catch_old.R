
#Script for evaluating spatial catch

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


setwd(handl_OneDrive("Analyses/Catch and effort"))
library(PBSmapping)   #for polygon

#DATA SECTION

#Total catch data
Data.monthly=read.csv("Data.monthly.csv")
Data.monthly.north=read.csv("Data.monthly.north.csv")

  #Depth
Depth.catch=read.csv("Daily.log.depth.csv")


#Onboard observing
    #Sharks data base
Observer=read.csv(handl_OneDrive("Analyses/Size and sex patterns/DATA.csv"))


#SST
Temp=read.csv(handl_OneDrive("Data/SST.nice.format.csv"))

PerthIs=read.table(handl_OneDrive("Data/Mapping/WAislandsPointsNew.txt"), header=T) #function for reading txt file
Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))

data(worldLLhigh)

#PARAMETERS SECTION

Blocks.north=c(21130,22130,23130)
names(Blocks.north)=c("Aatams.North","Aatams.Mid","Aatams.South")

Blocks=c(31150,32150,34140,35160,35180)
names(Blocks)=c("SMN.Perth1","SMN.Perth2","SMN.Hamelin","SMN.Bald","SMN.Albany")


#PROCEDURE SECTION

#1. Eastern-most distribution of dusky and sandbars with month

Tem1=aggregate(Temperature~Month+Year,subset(Temp,Year%in%2006:2013 & Lat <(-26) & Long>113),mean)
fn.plot.Temp=function()
{
  plot(Tem1$Temperature,ylab="",xaxt="n",col=col1,type='o',lwd=2,pch=19,cex=2)
  axis(1,seq(1,nrow(Tem1),12),unique(Tem1$Year),tck=-0.035)
  axis(1,1:nrow(Tem1),F,tck=-0.01)
  axis(1,seq(1,nrow(Tem1),6),F,tck=-0.02)
  mtext("Mean temperature (?C)",2,line=3,las=3,cex=1.25)
}


Axis1=function(D) axis(1,1:nrow(D),F,tck=-0.01)
Axis2=function(D) axis(1,seq(1,nrow(D),6),F,tck=-0.02)
Axis3=function(D) axis(1,seq(1,nrow(D),12),F,tck=-0.035)
Axis4=function(D) axis(1,seq(1,nrow(D),12),unique(D$YEAR.c),tck=-0.035)

fn.polygon=function(Var,COL1,COL2,LGND)
{
  Unic=unique(dat$Yr.Mn)  
  Low.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (0+Nper)/100))
  High.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (100-Nper)/100))
  
  lista=vector('list',length(Unic))
  for(i in 1:length(lista))
  {
    dummy=subset(dat,Yr.Mn==Unic[i])
    dummy=dummy[,match(Var,names(dummy))]
    lista[[i]]=matrix(dummy)
  }
  
  
  LOW.50=UP.50=LOW.75=UP.75=LOW.95=UP.95=LOW.100=UP.100=Med=rep(NA,length(lista))
  for(i in 1:length(lista))
  {
    DAT=lista[[i]]
    
    Med[i]=quantile(unlist(DAT),.5)
    
    Nper=(100-50)/2
    LOW.50[i]=Low.percentile(Nper)
    UP.50[i]=High.percentile(Nper)
    
    #75% of data
    Nper=(100-75)/2
    LOW.75[i]=Low.percentile(Nper)
    UP.75[i]=High.percentile(Nper)
    
    #95% of data
    Nper=(100-95)/2
    LOW.95[i]=Low.percentile(Nper)
    UP.95[i]=High.percentile(Nper)
    
    #100% of data
    Nper=(100-100)/2
    LOW.100[i]=Low.percentile(Nper)
    UP.100[i]=High.percentile(Nper)
  }  
  
  YR=1:length(lista)  
  
  #construct polygons
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec.50 <- c(LOW.50, tail(UP.50, 1), rev(UP.50), LOW.50[1])
  Biom.Vec.75 <- c(LOW.75, tail(UP.75, 1), rev(UP.75), LOW.75[1])
  Biom.Vec.95 <- c(LOW.95, tail(UP.95, 1), rev(UP.95), LOW.95[1])
  Biom.Vec.100 <- c(LOW.100, tail(UP.100, 1), rev(UP.100), LOW.100[1])
  
  YLIM=range(dat[,match(Var,names(dat))])
  #plot
  plot(YR,UP.95,ylim=YLIM,type="l",ylab="",xlab="",cex.axis=1.25,
       cex.lab=1.75,xaxt='n',col='transparent')
  
  colfunc <- colorRampPalette(c(COL1,COL2))
  COLS=colfunc(3)
  #polygon(Year.Vec, Biom.Vec.100, col = COLS[4], border = "grey20")
  polygon(Year.Vec, Biom.Vec.95, col = COLS[1], border = "grey20")
  polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
  polygon(Year.Vec, Biom.Vec.50, col = COLS[3], border = "grey20")
  #lines(YR,Med,lwd=2,col="black")
  if(LGND=="Yes")legend("top",c("95%","75%","50%"),horiz=T,bty='n',fill=COLS,title="Percentile",cex=1.5)
}


fun.long.dist=function(dat,WHAT,Var,COL1,COL2,LGND)
{
  East=aggregate(LONG~MONTH+YEAR.c,dat,max)
  East.max.yr=aggregate(LONG~YEAR.c,East,max)
  names(East.max.yr)[2]="Long.max"
  East=merge(East,East.max.yr,by="YEAR.c")
  East$stand=East$LONG/East$Long.max
  D=East
  
  if(WHAT=='raw')
  {
    plot(1:nrow(East),East$LONG,xaxt='n',ylab='',xlab='',col=COL1,type='o',lwd=2,pch=19,cex=2)
    Axis1(D);Axis2(D);Axis4(D)
  }
    
  if(WHAT=='standardised')
  {
    plot(1:nrow(East),East$stand,xaxt='n',ylab='',xlab='',col=COL1,type='o',lwd=2,pch=19,cex=2)
    Axis1(D);Axis2(D);Axis4(D)
  }

  
  if(WHAT=='range')
  {
    #ranges    
    fn.polygon(Var,col1,col2,LGND)
    Axis1(D);Axis2(D);Axis3(D)
  }
}

fun.lat.dist=function(dat,WHAT,Var,COL1,COL2,LGND)
{
  South=aggregate(LAT~MONTH+YEAR.c,dat,min)
  South.max.yr=aggregate(LAT~YEAR.c,South,min)
  names(South.max.yr)[2]="LAT.max"
  South=merge(South,South.max.yr,by="YEAR.c")
  South$stand=South$LAT/South$LAT.max
  D=South
  
  if(WHAT=='raw')
  {
    plot(1:nrow(South),South$LAT,xaxt='n',ylab='',xlab='',col=COL1,type='o',lwd=2,pch=19,cex=2)
    Axis1(D);Axis2(D);Axis4(D)
  }
    
  if(WHAT=='standardised')
  {
    plot(1:nrow(South),South$stand,xaxt='n',ylab='',xlab='',col=COL1,type='o',lwd=2,pch=19,cex=2)
    Axis1(D);Axis2(D);Axis4(D)
  }
    
  
  
  if(WHAT=='range')
  {
    #ranges
    fn.polygon(Var,col1,col2,LGND)
    Axis1(D);Axis2(D);Axis4(D)
  }
}

col1="cadetblue"
col2="lightblue"
#col1="olivedrab1"
#col2="olivedrab4"
  
setwd(handl_OneDrive("Analyses/Spatial catch"))

#Sandbar

  #range
dat=subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005)
dat$Yr.Mn=with(dat,paste(YEAR.c,MONTH))

tiff(file="Sandbar.range.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mai=c(.2,.6,.3,.05),oma=c(3,3,.1,.1),las=1)
fun.long.dist(dat,'range',"LONG",col1,col2,"Yes")
mtext("Longitude (?E)",2,line=3,las=3,cex=2)
fun.lat.dist(dat,'range',"LAT",col1,col2,"No")
mtext("Latitude (?S)",2,line=3,las=3,cex=2)
mtext("Year",1,outer=T,line=2,cex=2)
dev.off()

  #Raw data
tiff(file="Sandbar.raw.max.lat.long.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,1),mai=c(.2,.5,.3,.05),oma=c(3,3,.1,.1),las=1)
fun.long.dist(subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005),
              'raw',"LONG",col1,col2,"Yes")
mtext("Longitude (?E)",2,line=3,las=3,cex=1.25)
fun.lat.dist(subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005),
             'raw',"LAT",col1,col2,"No")
mtext("Latitude (?S)",2,line=3,las=3,cex=1.25)
fn.plot.Temp()
mtext("Year",1,outer=T,line=2,cex=2)
dev.off()

  #Standardised data
par(mfcol=c(3,1),mai=c(.2,.6,.3,.05),oma=c(3,3,.1,.1),las=1)
fun.long.dist(subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005),
              'standardised',"LONG",col1,col2,"Yes")
mtext("Longitude (?E)",2,line=3,las=3,cex=2)
fun.lat.dist(subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005),
             'standardised',"LAT",col1,col2,"No")
mtext("Latitude (?S)",2,line=3,las=3,cex=2)
fn.plot.Temp()
mtext("Year",1,outer=T,line=2,cex=2)


#Dusky
  
  #range
dat=subset(Data.monthly,SPECIES==18003 & METHOD=="GN" & YEAR.c>2005)
dat$Yr.Mn=with(dat,paste(YEAR.c,MONTH))

tiff(file="Dusky.range.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mai=c(.2,.6,.3,.05),oma=c(3,3,.1,.1),las=1)
fun.long.dist(dat,'range',"LONG",col1,col2,"Yes")
mtext("Longitude (?E)",2,line=3,las=3,cex=2)
fun.lat.dist(dat,'range',"LAT",col1,col2,"No")
mtext("Latitude (?S)",2,line=3,las=3,cex=2)
mtext("Year",1,outer=T,line=2,cex=2)
dev.off()


  #Raw data
tiff(file="Dusky.raw.max.lat.long.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,1),mai=c(.2,.5,.3,.05),oma=c(3,3,.1,.1),las=1)
fun.long.dist(subset(Data.monthly,SPECIES==18003 & METHOD=="GN" & YEAR.c>2005),
              'raw',"LONG",col1,col2,"Yes")
mtext("Longitude (?E)",2,line=3,las=3,cex=1.25)
fun.lat.dist(subset(Data.monthly,SPECIES==18003 & METHOD=="GN" & YEAR.c>2005),
             'raw',"LAT",col1,col2,"No")
mtext("Latitude (?S)",2,line=3,las=3,cex=1.25)
fn.plot.Temp()
mtext("Year",1,outer=T,line=2,cex=2)
dev.off()



#2. Monthly proportion of catch for acoustic array blocks
fn.prop.ktc=function(dat)
{
  Prop=aggregate(LIVEWT.c~BLOCKX+MONTH,dat,sum)
  Prop.max.yr=aggregate(LIVEWT.c~BLOCKX,Prop,sum)
  names(Prop.max.yr)[2]="LIVEWT.c.T"
  Prop=merge(Prop,Prop.max.yr,by=c("BLOCKX"))
  Prop=Prop[order(Prop$BLOCKX,Prop$MONTH),]
  Prop$stand=Prop$LIVEWT.c/Prop$LIVEWT.c.T
  
  BLKS=unique(Prop$BLOCKX)
  
  for(i in 1:length(BLKS))
  {
    a=subset(Prop,BLOCKX==BLKS[i])
    #plot(c(1,12),c(0,max(a$stand)),col="transparent",cex.main=1.5,
    plot(c(1,12),c(0,1),col="transparent",cex.main=1.5,
         ylab="",xlab="",main=paste("Block",BLKS[i]),cex.axis=1.3)
    
    lines(a$MONTH,a$stand,lwd=2.5,col=col1)
  }

}

fn.map=function(PLY,Blks)
{

  plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="grey75", axes=F, xlab="", ylab="",
          border="black",bg="grey99",plt = NULL)
  box()
  axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = seq(range.long[1],range.long[2],2),
       tck=-0.035,las=1,cex.axis=1.2)

  axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = seq(range.lat[1],range.lat[2],2),
       tck=-0.035,las=1,cex.axis=1.2)

  box()
  for(x in 1:length(PLY))
  {
    l=PLY[[x]]
    polygon(c(l[3],l[4],l[4],l[3]),c(l[2],l[2],l[1],l[1]),col=3)
    text(mean(l[3:4]),mean(l[1:2]),Blks[x])
  }
  mtext("Latitude (?S)",side=2,line=3,font=1,las=0,cex=1.5)
  mtext("Longitude (?E)",side=1,line=3,font=1,las=0,cex=1.5)
    
}


  #TDGDLF
range.long=c(114,119.5)
range.lat=c(-36,-31)

POLYS=list(c(-31,-32,115,116),c(-32,-33,115,116),c(-34,-35,114,115),c(-35,-36,116,117)
           ,c(-35,-36,118,119))

    #Dusky
tiff(file="Dusky.prop.ktch.receiver.blk.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,2),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc(subset(Data.monthly,SPECIES==18003 & METHOD=="GN" & YEAR.c>2005 & BLOCKX%in%Blocks))
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual dusky catch",2,outer=T,line=1,cex=2,las=3)
fn.map(POLYS,Blocks)
dev.off()


    #Sandbar
tiff(file="Sandbar.prop.ktch.receiver.blk.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,2),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc(subset(Data.monthly,SPECIES==18007 & METHOD=="GN" & YEAR.c>2005 &  BLOCKX%in%Blocks))
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual sandbar catch",2,outer=T,line=1,cex=2,las=3)
fn.map(POLYS,Blocks)
dev.off()



#3. Monthly proportion of catch for all blocks

fn.prop.ktc.all=function(dat,do3D,N)
{
  a=sort(table(dat$BLOCKX))
  this=as.numeric(names(a[a>N]))
  dat=subset(dat,BLOCKX%in%this)
  dat$Long=floor(dat$LONG)
  Prop=aggregate(LIVEWT.c~Long+MONTH,dat,sum)
  Prop.max.yr=aggregate(LIVEWT.c~Long,Prop,sum)
  names(Prop.max.yr)[2]="LIVEWT.c.T"
  Prop=merge(Prop,Prop.max.yr,by=c("Long"))
  Prop=Prop[order(Prop$Long,Prop$MONTH),]
  Prop$stand=Prop$LIVEWT.c/Prop$LIVEWT.c.T
  
  LGs=sort(unique(Prop$Long))
  
  #3D
  if(do3D=="Yes")
  {
    a=Prop[,match(c("Long","MONTH","stand"),names(Prop))]
    a=a[order(a$Long,a$MONTH),]
    Matrix.r=reshape(a,v.names ="stand",timevar="Long",idvar="MONTH",direction='wide')
    Matrix.r=Matrix.r[order(Matrix.r$MONTH),]
    z=as.matrix(Matrix.r[,c(2:ncol(Matrix.r))])
    persp(x=Matrix.r$MONTH,y=LGs,z, theta = 40,
          phi = 30, expand = 0.8, col = "cadetblue",ylab="Longitude (?E)",
          xlab="Month",zlab="Proportion of annual catch", ticktype = "detailed",
          cex.axis=1,cex.lab=1.55,border="cadetblue2",zlim=c(0,1))
  }

  #2D
  if(do3D=="No")
  {
    for(i in 1:length(LGs))
    {
      a=subset(Prop,Long==LGs[i])
#      plot(c(1,12),c(0,max(a$stand)),col="transparent",cex.main=1.5,
      plot(c(1,12),c(0,1),col="transparent",cex.main=1.5,
           ylab="",xlab="",main=paste(LGs[i],"-",(LGs[i]+1),"? East",sep=""),cex.axis=1.3)    
      lines(a$MONTH,a$stand,lwd=3,col=col1)
    }
  }

}


#Dusky
tiff(file="Dusky.prop.ktch.all.blk.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,4),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc.all(subset(Data.monthly,SPECIES==18003 & METHOD=="GN"),do3D="No",100)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual dusky catch",2,outer=T,line=1,cex=2,las=3)
dev.off()

tiff(file="Dusky.prop.ktch.all.blk.3D.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
fn.prop.ktc.all(subset(Data.monthly,SPECIES==18003 & METHOD=="GN"),do3D="Yes",100)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual dusky catch",2,outer=T,line=1,cex=2,las=3)
dev.off()


#Sandbar
tiff(file="Sandbar.prop.ktch.all.blk.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,2),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc.all(subset(Data.monthly,SPECIES==18007 & METHOD=="GN"),do3D="No",100)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual sandbar catch",2,outer=T,line=1,cex=2,las=3)
dev.off()

tiff(file="Sandbar.prop.ktch.all.blk.3D.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
fn.prop.ktc.all(subset(Data.monthly,SPECIES==18007 & METHOD=="GN"),do3D="Yes",100)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual sandbar catch",2,outer=T,line=1,cex=2,las=3)
dev.off()



    #NSF

      #Sandbar
tiff(file="Sandbar.prop.ktch.all.blk_NSF.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,4),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc.all(subset(Data.monthly.north,SPECIES==18007 & YEAR.c>2000),do3D="No",10)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual sandbar catch",2,outer=T,line=1,cex=2,las=3)
dev.off()

tiff(file="Sandbar.prop.ktch.all.blk.3D_NSF.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
fn.prop.ktc.all(subset(Data.monthly.north,SPECIES==18007 & YEAR.c>2000),do3D="Yes",10)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual sandbar catch",2,outer=T,line=1,cex=2,las=3)
dev.off()


      #Dusky
tiff(file="Dusky.prop.ktch.all.blk_NSF.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mai=c(.2,.4,.45,.05),oma=c(3,3,.1,.1),mgp=c(1,.7,0),las=1)
fn.prop.ktc.all(subset(Data.monthly.north,SPECIES==18003 & YEAR.c>2000),do3D="No",10)
mtext("Month",1,outer=T,line=1.5,cex=2)
mtext("Proportion of annual dusky catch",2,outer=T,line=1,cex=2,las=3)
dev.off()


#4. Depth
fn.depth=function(dat)
{
  dat=subset(dat,depthMax<250)

  par(mfcol=c(2,2))
  boxplot(depthMax~LongDeg,dat,main="Long")
  boxplot(depthMax~month,dat,main="Month")
  boxplot(depthMax~as.factor(finyear),dat,main="Finyear")
  d1=aggregate(depthMax~LongDeg+month,dat,mean,na.rm=T)
  a=d1[order(d1$LongDeg,d1$month),]
  Matrix.r=reshape(a,v.names ="depthMax",timevar="LongDeg",idvar="month",direction='wide')
  Matrix.r=Matrix.r[order(Matrix.r$month),]
  z=as.matrix(Matrix.r[,c(2:ncol(Matrix.r))])
  persp(x=Matrix.r$month,y=sort(unique(d1$LongDeg)),z, theta = 40,
        phi = 30, expand = 0.8, col = "cadetblue",ylab="Longitude (?E)",
        xlab="Month",zlab="Depth (m)", ticktype = "detailed",
        cex.axis=1,cex.lab=1.55,border="cadetblue2")
  mtext("Depth (m)",2,outer=T,las=3,cex=2,line=-2.5)
  
}

tiff(file="Sandbar.depth.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
fn.depth(subset(Depth.catch,species==18007 & LatDeg>26))
dev.off()

tiff(file="Dusky.depth.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
fn.depth(subset(Depth.catch,species==18003 & LatDeg>26))
dev.off()



#5. Observer data

#5.1 Mean size
fn.obs.long=function(dat,METHOD)
{
  boxplot(FL~Long.round,dat,main=METHOD,xlab="Long")
  boxplot(FL~Month,dat,xlab="Month")
  boxplot(FL~year,dat,xlab="Year")
  
  
#   d1=aggregate( FL~Long.round+Month,dat,mean,na.rm=T)
#   a=d1[order(d1$Long.round,d1$Month),]
#   Matrix.r=reshape(a,v.names ="FL",timevar="Long.round",idvar="Month",direction='wide')
#   Matrix.r=Matrix.r[order(Matrix.r$Month),]
#   z=as.matrix(Matrix.r[,c(2:ncol(Matrix.r))])
#   persp(x=Matrix.r$Month,y=sort(unique(d1$Long.round)),z, theta = 40,
#         phi = 30, expand = 0.8, col = "cadetblue",ylab="Longitude (?E)",
#         xlab="Month",zlab="Depth (m)", ticktype = "detailed",
#         cex.axis=1,cex.lab=1.55,border="cadetblue2")
   mtext("FL (cm)",2,outer=T,las=3,cex=2,line=-3)
  
}

fn.obs.lat=function(dat,METHOD)
{
  boxplot(FL~Lat.round,dat,main=METHOD,xlab="Lat")
  boxplot(FL~Month,dat,xlab="Month")
  boxplot(FL~year,dat,xlab="Year")
  mtext("FL (cm)",2,outer=T,las=3,cex=2,line=-3)
}
a=function() par(mfcol=c(3,2),mai=c(.6,.7,.2,.1))


  #Dusky shark

    #Long pattern
tiff(file="Dusky.FL.obs.long.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
a()
fn.obs.long(subset(Observer,SPECIES=="BW" & Lat.round<(-33) & Long.round>112 & Method=="GN"),"Gillnet")
fn.obs.long(subset(Observer,SPECIES=="BW" & Lat.round<(-33) & Long.round>112 & Method=="LL"),"Longline")
dev.off()

    #Lat pattern
tiff(file="Dusky.FL.obs.lat.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
a()
fn.obs.lat(subset(Observer,SPECIES=="BW" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="GN"),"Gillnet")
fn.obs.lat(subset(Observer,SPECIES=="BW" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="LL"),"Longline")
dev.off()



  #Sandbar shark

    #Long pattern
tiff(file="Sandbar.FL.obs.long.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
a()
fn.obs.long(subset(Observer,SPECIES=="TK" & Lat.round<(-33) & Long.round>112 & Method=="GN"),"Gillnet")
fn.obs.long(subset(Observer,SPECIES=="TK" & Lat.round<(-33) & Long.round>112 & Method=="LL"),"Longline")
dev.off()

    #Lat pattern
tiff(file="Sandbar.FL.obs.lat.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
a()
fn.obs.lat(subset(Observer,SPECIES=="TK" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="GN"),"Gillnet")
fn.obs.lat(subset(Observer,SPECIES=="TK" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="LL"),"Longline")
dev.off()



#5.2 Histogram
fn.obs.long=function(dat,METHOD)
{
  y=sort(unique(dat$Long.round))
  for(i in 1:length(y))
  {
    da1=subset(dat,Long.round==y[i])
    hist(da1$FL,xlim=c(40,max(dat$FL,na.rm=T)),main=y[i],xlab="FL",ylab="")
  }
}


fn.obs.lat=function(dat,METHOD)
{
  y=sort(unique(dat$Lat.round))
  for(i in 1:length(y))
  {
    da1=subset(dat,Lat.round==y[i])
    if(nrow(da1)>1)hist(da1$FL,xlim=c(40,max(dat$FL,na.rm=T)),main=y[i],xlab="FL",ylab="")
  }

}



b=function() par(mfcol=c(2,4),mai=c(.6,.7,.2,.1))


#Dusky shark

#Long pattern
tiff(file="Dusky.FL.obs.long.hist.GN.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,3),mai=c(.6,.7,.2,.1))
fn.obs.long(subset(Observer,SPECIES=="BW" & Lat.round<(-33) & Long.round>112 & Method=="GN"),"Gillnet")
dev.off()

tiff(file="Dusky.FL.obs.long.hist.LL.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.6,.7,.2,.1))
fn.obs.long(subset(Observer,SPECIES=="BW" & Lat.round<(-33) & Long.round>112 & Method=="LL"),"Longline")
dev.off()

#Lat pattern
tiff(file="Dusky.FL.obs.lat.hist.GN.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mai=c(.6,.7,.2,.1))
fn.obs.lat(subset(Observer,SPECIES=="BW" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="GN"),"Gillnet")
dev.off()

tiff(file="Dusky.FL.obs.lat.hist.LL.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,3),mai=c(.6,.7,.2,.1))
fn.obs.lat(subset(Observer,SPECIES=="BW" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="LL"),"Longline")
dev.off()



#Sandbar shark
#Long pattern
tiff(file="Sandbar.FL.obs.long.hist.GN.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mai=c(.6,.7,.2,.1))
fn.obs.long(subset(Observer,SPECIES=="TK" & Lat.round<(-33) & Long.round>112 & Method=="GN"),"Gillnet")
dev.off()

tiff(file="Sandbar.FL.obs.long.hist.LL.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(1,1),mai=c(.6,.7,.2,.1))
fn.obs.long(subset(Observer,SPECIES=="TK" & Lat.round<(-33) & Long.round>112 & Method=="LL"),"Longline")
dev.off()

#Lat pattern
tiff(file="Sandbar.FL.obs.lat.hist.GN.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mai=c(.6,.7,.2,.1))
fn.obs.lat(subset(Observer,SPECIES=="TK" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="GN"),"Gillnet")
dev.off()

tiff(file="Sandbar.FL.obs.lat.hist.LL.tiff",width = 2200, height = 2300,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,3),mai=c(.6,.7,.2,.1))
fn.obs.lat(subset(Observer,SPECIES=="TK" & Lat.round<(-10) & Long.round<115 & Long.round>112 & Method=="LL"),"Longline")
dev.off()




