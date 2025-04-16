#Some functions

PLOT.fn=function(DATA1,DATA2,max1,max2,LEG,CL2)   #function plot vessels and blocks
{
  plot(1:NN.monthly,DATA1,ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1,ylim=c(0,max1)
       ,cex.axis=.75,lwd=1.1)
  axis(1,at=1:NN.monthly,labels=F,tck=-0.015)
  
  par(new=T)
  plot(1:NN.monthly,DATA2,col=CL2,type='o',pch=19,axes=F,ann='F',cex=1,lwd=1.1,ylim=c(0,max2))
  #  axis(4,at=pretty(DATA2),labels=pretty(DATA2),las=2,cex.axis=.75)
  axis(4,at=seq(0,max2,40),labels=seq(0,max2,40),las=2,cex.axis=.75,col.axis=CL2)
  legend("topleft",LEG,bty='n',cex=.75)
}

fn.eff.probs=function(gear,eff.var,eff.var1,letra,CLs)
{
  fn.wrong=function()
  {
    ID1=match(corr,names(dummy))
    dummy$VAR=ifelse(!(dummy[,ID]==dummy[,ID1]),"Wrong","OK")
    dummy$VAR=with(dummy,ifelse(is.na(VAR),"Wrong",VAR))
    Agg=table(dummy[,ID2],dummy$VAR)
    Agg=as.matrix(round(100*Agg/rowSums(Agg),4))
    if(ncol(Agg)==2)Tab=cbind(Agg[,1],Agg[,2])
    if(ncol(Agg)==1)Tab=cbind(Agg[,1],rep(0,length(Agg[,1])))
    colnames(Tab)=c("ok","wrong")
    return(Tab)
  }
  
  #monthly
  dummy=subset(Effort.monthly,!FINYEAR%in%Daily.l.years & METHOD==gear & !(BLOCKX%in%Estuaries))
  corr=paste(eff.var,".c",sep="")
  dummy=dummy[,match(c("FINYEAR","METHOD","Same.return",eff.var,corr),names(dummy))]
  dummy=dummy[!duplicated(dummy$Same.return),]
  THIS="FINYEAR"
  eff.var.1=eff.var
  ID=match(eff.var,names(dummy))
  ID2=match(THIS,names(dummy))
  
  
  Tab.m=fn.wrong()
  
  #daily
  dummy=subset(Effort.daily,method==gear & !(blockx%in%Estuaries))
  corr=paste(eff.var1,".c",sep="")
  dummy=dummy[,match(c("finyear","method","date","ID",eff.var1,corr),names(dummy))]
  dummy$ID=with(dummy,paste(date,ID))
  dummy=dummy[!duplicated(dummy$ID),]
  ID=match(eff.var1,names(dummy))
  ID2=match("finyear",names(dummy))
  
  Tab.d=fn.wrong()
  
  TAB=rbind(Tab.m,Tab.d)
  rownames(TAB)=NULL
  barplot(rbind(TAB[,1],TAB[,2]),beside=F,col=CLs,cex.axis=1.25,cex.lab=1.25,
          ylim=c(0,110),yaxs="i",xaxs="i",legend.text=paste(letra),          
          args.legend=list(x = "topleft",cex=1,bty='n',
                           fill="transparent",border='transparent'),ylab="")
  box()
  AXIS1()
  AXIS2()
  return(data.frame(OK=mean(TAB[,1]),wrong=mean(TAB[,2])))
}

fn.explore=function(DATA,NAMES)
{
  #remove initial years when sandbar was not reported
  if(NAMES=="Sandbar shark")
  {
    DATA=subset(DATA,!(FINYEAR%in%c("1975-76","1976-77","1977-78",
                                    "1978-79","1979-80","1980-81","1981-82","1982-83","1983-84","1984-85")))
  }
  
  #Remove vessel levels that don't occur in data
  DATA$VESSEL=DATA$VESSEL[, drop=TRUE]
  
  DATA$CPUE=with(DATA,Catch.Target/Km.Gillnet.Days.c)
  DATA$SP.Target=ifelse(DATA$CPUE>0,1,0)
  
  #Tables
  TABLE1=sort(with(DATA,table(VESSEL)))
  TABLE1=data.frame(VESSEL=names(TABLE1),Count=as.numeric(TABLE1))
  
  TABLE6=aggregate(Catch.Target~VESSEL,data=DATA,mean)
  names(TABLE6)[2]="Mean Catch"
  TABLE6.1=aggregate(Catch.Target~VESSEL,data=DATA,sd)
  names(TABLE6.1)[2]="SD Catch"
  TABLE6=merge(TABLE6,TABLE1,by="VESSEL")
  TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
  TABLE6=TABLE6[order(TABLE6$Count),]
  
  
  #cumulative catch
  TABLE12=aggregate(Catch.Target~VESSEL,data=DATA,sum)
  TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
  TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
  TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
  
  TABLE13=aggregate(Catch.Target~BLOCKX,data=DATA,sum)
  TABLE13=TABLE13[order(-TABLE13$Catch.Target),]
  TABLE13$CumCatch=cumsum(TABLE13$Catch.Target)
  TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$Catch.Target),2)
  
  #cumulative records
  TABLE16=rev(sort(table(DATA$VESSEL)))
  TABLE16=data.frame(Records=as.numeric(TABLE16),VESSEL=names(TABLE16))
  TABLE16$CumRecords=cumsum(TABLE16$Records)
  TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
  
  TABLE17=rev(sort(table(DATA$BLOCKX)))
  TABLE17=data.frame(Records=as.numeric(TABLE17),BLOCKX=names(TABLE17))
  TABLE17$CumRecords=cumsum(TABLE17$Records)
  TABLE17$PerCumRecords=round(TABLE17$CumRecords*100/sum(TABLE17$Records),2)
  
  
  
  return(list(Ves.Cum.Ca=TABLE12$PerCumCatch,Block.Cum.Ca=TABLE13$PerCumCatch,
              TABLE17=TABLE17,TABLE16=TABLE16,Mean_catch_Vessel=TABLE6))
}

fun.rec.per.ves=function(DATA)
{
  Count.ves.records=table(DATA$Count)
  Count.ves.records=c(subset(Count.ves.records,as.numeric(names(Count.ves.records))<Min.rec.ves),sum(
    subset(Count.ves.records,as.numeric(names(Count.ves.records))<Min.rec.ves)))
  names(Count.ves.records)[length(Count.ves.records)]=paste(Min.rec.ves-1,"+",sep="")
  b=barplot(Count.ves.records,xaxt='n',yaxt='n')
  axis(1,at=b,labels=F,tck=-0.016)
  axis(2,at=seq(0,300,100),labels=seq(0,300,100),tck=-0.016,cex=1.25,las=1)
  axis(1,at=b[seq(2,(length(b)-2),by=2)],labels=names(Count.ves.records)[seq(2,(length(b)-2),by=2)],
       tck=-0.032,cex=1.25)
  mtext("Number of records per vessel",side=1,line=1.2,font=1,las=0,cex=1.125,outer=F)
  mtext("Count",side=2,line=1.75,font=1,las=0,cex=1.25,outer=F)
  axis(1,at=b[length(b)],labels=paste(">",19,sep=""),tck=-0.032,cex.axis=1.15)
  
}

fn.Figure5=function(DATA,NAMES,SP)
{
  #remove initial years when sandbar was not reported
  if(NAMES=="Sandbar shark")
  {
    DATA=subset(DATA,!(FINYEAR%in%c("1975-76","1976-77","1977-78",
                                    "1978-79","1979-80","1980-81","1981-82","1982-83","1983-84","1984-85")))
  }
  
  #Remove vessel levels that don't occur in data
  DATA$VESSEL=DATA$VESSEL[, drop=TRUE]
  
  DATA=subset(DATA,SPECIES==SP)
  DATA$CPUE=with(DATA,LIVEWT.c/Km.Gillnet.Hours.c)
  DATA$SP.Target=ifelse(DATA$CPUE>0,1,0)
  DATA$Catch.Target=DATA$LIVEWT.c
  
  #Tables
  TABLE1=sort(with(DATA,table(VESSEL)))
  
  TABLE6=aggregate(Catch.Target~VESSEL,data=DATA,mean)
  names(TABLE6)[2]="Mean Catch"
  TABLE6.1=aggregate(Catch.Target~VESSEL,data=DATA,sd)
  names(TABLE6.1)[2]="SD Catch"
  TABLE6=merge(TABLE6,data.frame(VESSEL=names(TABLE1),Count=as.numeric(TABLE1)),by="VESSEL")
  TABLE6=merge(TABLE6,TABLE6.1,by="VESSEL")
  TABLE6=TABLE6[order(TABLE6$Count),]
  
  
  #cumulative catch
  TABLE12=aggregate(Catch.Target~VESSEL,data=DATA,sum)
  TABLE12=TABLE12[order(-TABLE12$Catch.Target),]
  TABLE12$CumCatch=cumsum(TABLE12$Catch.Target)
  TABLE12$PerCumCatch=round(TABLE12$CumCatch*100/sum(TABLE12$Catch.Target),2)
  
  TABLE13=aggregate(Catch.Target~BLOCKX,data=DATA,sum)
  TABLE13=TABLE13[order(-TABLE13$Catch.Target),]
  TABLE13$CumCatch=cumsum(TABLE13$Catch.Target)
  TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$Catch.Target),2)
  
  #cumulative records
  TABLE16=rev(sort(table(DATA$VESSEL)))
  TABLE16=data.frame(Records=as.numeric(TABLE16),VESSEL=names(TABLE16))
  TABLE16$CumRecords=cumsum(TABLE16$Records)
  TABLE16$PerCumRecords=round(TABLE16$CumRecords*100/sum(TABLE16$Records),2)
  
  TABLE17=rev(sort(table(DATA$BLOCKX)))
  TABLE17=data.frame(Records=as.numeric(TABLE17),BLOCKX=names(TABLE17))
  TABLE17$CumRecords=cumsum(TABLE17$Records)
  TABLE17$PerCumRecords=round(TABLE17$CumRecords*100/sum(TABLE17$Records),2)
  
  
  
  return(list(Ves.Cum.Ca=TABLE12$PerCumCatch,Block.Cum.Ca=TABLE13$PerCumCatch,
              TABLE17=TABLE17,TABLE16=TABLE16,Mean_catch_Vessel=TABLE6))
}

#Extract monthly records by Year-Month-Vessel-Block

These.efforts=c("FINYEAR","Same.return","Km.Gillnet.Days.inv",
                "Km.Gillnet.Days.c","zone","MONTH","BLOCKX")

Effort.data.fun=function(DATA,target,ktch)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  # remove nonsense lat
  DATA=subset(DATA,LAT>=(-36))
  
  #calculate effort
  Match.these.eff=match(These.efforts,names(DATA))
  Effort.data=DATA[,Match.these.eff]
  Effort.data=aggregate(cbind(Km.Gillnet.Days.inv,Km.Gillnet.Days.c)~zone+
                          FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,max)
  Effort.data=aggregate(cbind(Km.Gillnet.Days.inv,Km.Gillnet.Days.c)~zone+
                          FINYEAR+Same.return+MONTH+BLOCKX,Effort.data,sum)
  
  #target species catch 
  ID=match(c(ktch),colnames(DATA))
  DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
  
  
  #catch targeted at other species
  DATA$Catch.Gummy=with(DATA,ifelse(SPECIES==17001,DATA[,ID],0))
  DATA$Catch.Whiskery=with(DATA,ifelse(SPECIES==17003,DATA[,ID],0))
  DATA$Catch.Dusky=with(DATA,ifelse(SPECIES%in%c(18003),DATA[,ID],0))
  DATA$Catch.Sandbar=with(DATA,ifelse(SPECIES==18007,DATA[,ID],0))
  DATA$Catch.Scalefish=with(DATA,ifelse(SPECIES%in%188000:599001,DATA[,ID],0))
  DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
  
  #reshape catch data
  TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Scalefish,
                        Catch.Dusky,Catch.Sandbar,Catch.Total)~MONTH+
                    FINYEAR+BLOCKX+VESSEL+Same.return+LAT+LONG+
                    YEAR.c,data=DATA,sum,na.rm=T)
  TABLE=TABLE[order(TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
  
  #proportion of records with target catch
  prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
  
  #merge catch and effort
  dat=merge(TABLE,Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"),all.x=T)
  
  #create "other shark catch" variable
  dat$Catch.other.shk=NA
  if(target[1]==17003)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Sandbar
  if(target[1]==17001)dat$Catch.other.shk=dat$Catch.Whiskery
  if(target[1]%in%c(18003))dat$Catch.other.shk=dat$Catch.Whiskery+dat$Catch.Sandbar
  if(target[1]==18007)dat$Catch.other.shk=dat$Catch.Dusky+dat$Catch.Whiskery
  
  #recalculate 60 by 60 blocks
  dat$BLOCKX.orignl=dat$BLOCKX
  dat$BLOCKX=as.numeric(substr(dat$BLOCKX,1,4))
  
  
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}

DATA.list.LIVEWT.c=vector('list',length=N.species)
names(DATA.list.LIVEWT.c)=Nms.sp[Tar.sp]   


#Create data sets for plotting cpue paper figures
for ( i in Tar.sp)DATA.list.LIVEWT.c[[i]]=Effort.data.fun(DATA=Species.list[[i]],
                                                          target=TARGETS[[i]],
                                                          ktch="LIVEWT.c")$dat

#Create figures 1 to 5
tiff(file="Figure 1. Map.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mar=c(2,2,2,2),oma=c(1,1,1,1))
plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
polygon(x=c(116.5,116.5,112,112),y=c(-26.5,-33,-33,-26.5),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
polygon(x=c(116.5,116.5,112,112),y=c(-33,-37,-37,-33),lwd=1.5,col=rgb(.3,.3,.3,alpha=.5))
polygon(x=c(129,129,116.5,116.5),y=c(-30,-37,-37,-30),lwd=1.5,col=rgb(.7,.7,.7,alpha=.2))

axis(side = 1, at =seq(LONGG[1],LONGG[length(LONGG)],length.out = 7+6*(length(LONGG)-2)), labels = F, tcl = 34,lty=3,col="grey60")
axis(side = 4, at = seq(LATT[1],LATT[length(LATT)],length.out = 7+6*(length(LATT)-2)), labels = F,tcl =34,lty=3,col="grey30")
axis(side = 1, at =LONGG, labels = F, tcl = 34,lty=1,col="grey30")
axis(side = 4, at = LATT, labels = F,tcl =34,lty=1,col="grey30")

if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=c(-37,-25),xlim=South.WA.long, zlim=c(-1,-300),
                             nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)

par(new=T,mar=c(2,2,2,2),oma=c(1,1,1,1))
plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
axis(side = 1, at =seq(112,129,2), labels = seq(112,129,2), tcl = .35,las=1,cex.axis=1.25,padj=-1.25)
axis(side = 2, at = seq(-36,-25,2), labels = -seq(-36,-25,2),tcl = .35,las=2,cex.axis=1.25,hadj=.3)
text(116.73,Perth[2],("Perth"),col="black", cex=1.1)
points(115.86,-31.95,pch=19,cex=1.5)
text(116.73,-33.55,("Bunbury"),col="black", cex=1.1)
points(115.6,-33.55,pch=19,cex=1.5)
text(117.7,-34.75,("Albany"),col="black", cex=1.1)
points(117.8,-35,pch=19,cex=1.5)
text(122,-33.62,("Esperance"),col="black", cex=1.1)
points(121.9,-33.86,pch=19,cex=1.5)
mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1.2,las=3,cex=1.75)
mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.75,cex=1.75)

text(113.5,-30.5,("WCDGDLF"),col="black", cex=1.4)
text(114,-34.75,("JASDGDLF"),col="black", cex=1.4) 
text(114,-35.55,("(Zone 1)"),col="black", cex=1.4)
text(122,-34.75,("JASDGDLF"),col="black", cex=1.4)
text(122,-35.55,("(Zone 2)"),col="black", cex=1.4)

par(fig=c(.5,.92,.5,.92), new = T,mgp=c(.1,.4,0))
plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
        col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
box()
polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
text(134,-22.5,("Australia"),col="black", cex=2)
dev.off()


#Appendix 1.                

#Stack-up catch plot
these.ones=c(22999,17001,17003,18007,18003)
names(these.ones)=c("SHARK, OTHER","SHARK, GUMMY","SHARK, WHISKERY",
                    "SHARK, THICKSKIN (SANDBAR)","SHARK, BRONZE WHALER")
#Monthly
FINYEAR.monthly=sort(unique(Data.monthly.GN$FINYEAR))
Mn.yr=subset(FINYEAR.monthly,!FINYEAR.monthly %in% FINYEAR.daily)
NN.monthly=length(Mn.yr)
STORE=matrix(nrow=NN.monthly,ncol=length(these.ones))
colnames(STORE)=sort(these.ones)
for(i in 1:NN.monthly)
{
  Dat=subset(Data.monthly.GN, FINYEAR==FINYEAR.monthly[i])    
  test=aggregate(LIVEWT~FINYEAR+Spec.old,data=Dat,sum,na.rm=T)
  test$SP=ifelse(test$Spec.old%in%these.ones,test$Spec.old,"OTHER")
  test=aggregate(LIVEWT~FINYEAR+ SP,data=test,sum,na.rm=T)
  prop=100*test[,3]/sum(test[,3])
  names(prop)=test[,2]
  if(!(names(prop)[4]=="18007"))prop=c(prop[1:3],0,prop[4])
  if(names(prop)[4]=="18007")prop=prop[1:5]
  STORE[i,]=prop
}

#Daily
NN.daily=length(FINYEAR.daily)
STORE.daily=matrix(nrow=NN.daily,ncol=length(these.ones))
colnames(STORE.daily)=sort(these.ones)
for(i in 1:NN.daily)
{
  Dat=subset(Data.daily.GN, FINYEAR==FINYEAR.daily[i])    
  Dat$Spec.old=Dat$SPECIES
  test=aggregate(LIVEWT~FINYEAR+Spec.old,data=Dat,sum,na.rm=T)
  test$SP=ifelse(test$Spec.old%in%these.ones,test$Spec.old,"OTHER")
  test=aggregate(LIVEWT~FINYEAR+ SP,data=test,sum,na.rm=T)
  prop=100*test[,3]/sum(test[,3])
  names(prop)=test[,2]
  if(!(names(prop)[4]=="18007"))prop=c(prop[1:3],0,prop[4])
  if(names(prop)[4]=="18007")prop=prop[1:5]
  STORE.daily[i,]=prop
}
STORE=rbind(STORE,STORE.daily)


#Stack-up barplot of effort problems
COL.BAR=c("white","grey35","grey55","grey75","black")
NN=NN.monthly+NN.daily
AXIS1=function()axis(1,at=b,labels=F,tck=-0.016)
AXIS2=function()axis(1,at=b[seq(1,NN,by=5)],labels=F,tck=-0.03)
AXIS3=function()axis(1,at=b[seq(1,NN,by=5)],labels=FINYEAR.monthly[seq(1,NN,by=5)],tck=-0.035,cex.axis=1.25)

tiff(file="Appendix 1. Data problems_All.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
layout(matrix(c(1,1,2:5), 3, 2, byrow = TRUE))
par(mar=c(2,.75,.25,2),oma=c(2,4,.25,.01),las=1,mgp=c(1,.75,0))

#Catch (species id problem)
b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
          ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
          args.legend=list(x = "topleft",horiz=T,cex=1.7,pt.cex=1.7,bty='n',inset=c(0, -0.025)),ylab="")
box()
AXIS1()
AXIS3()

#BDAYS                                        
print(fn.eff.probs("GN","BDAYS","bdays","",c("grey85","black"))  )
mtext("Number of days fished per month",3,line=-1.375,cex=1.25)

#HOURS
print(fn.eff.probs("GN","HOURS","hours","",c("grey85","black")))
mtext("Number of hours fished per day",3,line=-1.375,cex=1.25)

#SHOTS
print(fn.eff.probs("GN","SHOTS","shots","",c("grey85","black")))
AXIS3()
mtext("Number of shots per day",3,line=-1.375,cex=1.25)

#NETLEN     
print(fn.eff.probs("GN","NETLEN","netlen","",c("grey85","black")))
AXIS3()
mtext("Net length per shot",3,line=-1.375,cex=1.25)

mtext("Financial year",side=1,line=0.8,font=1,las=0,cex=1.75,outer=T)
mtext("Percentage",side=2,line=2,font=1,las=0,cex=1.75,outer=T)
dev.off()


tiff(file="Appendix 1. Data problems.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(2,1),mar=c(1,3,.1,.1),oma=c(2.5,.5,.1,.5),las=1,mgp=c(1.8,.6,0))

#Catch (species id problem)
b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
          ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
          args.legend=list(x = "topleft",horiz=T,cex=1.17,pt.cex=1.25,bty='n',inset=c(0, -0.025)),ylab="")
box()
AXIS1()
AXIS2()

# #NETLEN     
fn.eff.probs("GN","NETLEN","netlen","",c("grey85","black"))
AXIS3()

mtext("Financial year",side=1,line=1,font=1,las=0,cex=1.75,outer=T)
mtext("Percentage",side=2,line=-0.75,font=1,las=0,cex=1.75,outer=T)
dev.off()


tiff(file="Appendix 1. Data problems_catch_only.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(1,1),mar=c(1,3,.15,.1),oma=c(3,.5,.15,.65),las=1,mgp=c(1.8,.9,0))

#Catch (species id problem)
b=barplot(rbind(STORE[,1],STORE[,2],STORE[,3],STORE[,4],STORE[,5]),beside=F,col=COL.BAR,cex.axis=1.25,cex.lab=1.25,
          ylim=c(0,110),yaxs="i",xaxs="i",legend.text=c("Gummy","Whiskery","Dusky","Sandbar","Unid. shark"),
          args.legend=list(x = "topleft",horiz=T,cex=1.17,pt.cex=1.25,bty='n',inset=c(0, -0.015)),ylab="")
box()
AXIS1()
AXIS2()
AXIS3()
mtext("Financial year",side=1,line=1.5,font=1,las=0,cex=1.75,outer=T)
mtext("Percentage",side=2,line=-0.75,font=1,las=0,cex=1.75,outer=T)
dev.off()


tiff(file="Appendix 1. Data problems_effort_only.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(2,2),mar=c(1,3.2,.1,.1),oma=c(2.5,.5,.1,.65),las=1,mgp=c(1.8,.6,0))

#BDAYS                                        
print(fn.eff.probs("GN","BDAYS","bdays","(a)",c("grey85","black"))  )

#HOURS
print(fn.eff.probs("GN","HOURS","hours","(b)",c("grey85","black")))
AXIS3()
#SHOTS
print(fn.eff.probs("GN","SHOTS","shots","(c)",c("grey85","black")))

#NETLEN     
print(fn.eff.probs("GN","NETLEN","netlen","(d)",c("grey85","black")))
AXIS3()

mtext("Financial year",side=1,line=1.25,font=1,las=0,cex=1.75,outer=T)
mtext("Percentage",side=2,line=-1,font=1,las=0,cex=1.75,outer=T)
dev.off()



#Appendix 2. Flowchart diagram 
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/flow_chart.R"))
#text
LABELS <- list("raw data",
               c("Exclude records from estuaries, non-gillnet gear","and school shark and dogfish targeting"),
               c("Valid catch composition"),
               c("Correct catch but do not use record"," in catch and effort standardisation"),
               "Valid effort",
               c("Correct effort but do not use record","in catch and effort standardisation"),
               c("Financial year > 1989-90"),
               c("Adjust incomplete","catch and effort"),
               c("Record within effective area"),
               c("Do not use record in","catch and effort standardisation"),
               c("Standardise catch and effort","Construct abundance index"),
               c("Adjust for increase in","fishing efficiency"))

#type of shape
#Note:
#oval: start and terminal points
#square or round: process
#diammond: decision (yes/no)

SHAPES=c("oval","round","diamond","round","diamond","round","diamond",
         "round","diamond","round","round","oval")

#Shape coordinates
MaInX=0.725
X.COOR=c(rep(MaInX,3),MaInX*0.35,MaInX,MaInX*0.35,MaInX,MaInX*0.35,
         MaInX,MaInX*0.35,MaInX,MaInX)
N.labl=length(LABELS)
Y.COOR=rep(NA,N.labl)
Y.COOR[1]=0.975
delta=c(rep(0.115,3),rep(0.07,8),0.115)
for(q in 2:N.labl) Y.COOR[q]=Y.COOR[q-1]-delta[q]

#Shape size
X.size=c(.08,.19,.13,.165,.09,.165,.125,.125,.125,.125,.13,.125)

#arrows
ArROW=c(rep("Straight",2),rep(c("Side.left","Side.right"),4),"Straight","Side.back")

tiff(file="Appendix 2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mar=c(0.1,0.1,0.1,0.1))
fn.flow.chart(lab=LABELS,SHPE=SHAPES,X.coor=X.COOR,Y.coor=Y.COOR,SX=X.size,ARRw=ArROW,CEX=.9,n=8,n1=9)
dev.off()  


#Appendix 3. Effort dynamics (expansion and contraction)
NN.monthly=NN
Lat.seq=c(-26,-28,-30,-32,-34,-36)

tiff(file="Appendix 3. Effort dynamics.gillnets.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
smart.par(n.plots=(1+length(DATA.lista)),MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(.1, 0.5, 0))
#number of blocks and vessels per year
PLOT.fn(BLKS.YEAR,Effort.expan$N.ves.yr,60,180,"",CL2="grey65")
axis(1,at=seq(1,NN.monthly,5),labels=FINYEAR.monthly[seq(1,NN.monthly,5)],tck=-0.03,cex.axis=.9)
mtext("Financial year",side=1,line=1.5,font=1,las=0,cex=1,outer=F)
mtext("Number of blocks fished",side=2,line=1.35,font=1,las=0,cex=.85,outer=F)
mtext("Number of licence holders fishing",side=4,line=1.3,las=3,cex=.75,outer=F,col="grey65")

#effort by block per calendar year groups
for (i in 1:length(DATA.lista))
{
  DATA=DATA.lista[[i]][-which(duplicated(DATA.lista[[i]]$Same.return)),]
  
  fn.eff.plot(DATA,tcl.1=0,tcl.2=0,EffortBreakSS)
  mtext(Yr.range[i],side=3,line=-1.25,cex=.95)
  axis(side = 1, at =Long.seq, labels = F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
  axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
  if(i%in%c(10,7,8)) axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
  if(i%in%c(3,6,9)) axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
  if(i==8) color.legend(126,-26,129,-30.5,round(EffortBreakSS,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
}  
mtext("Latitude (?S)",side=2,line=0.4,las=3,cex=1.25,outer=T)
mtext("Longitude (?E)",side=1,line=0.75,cex=1.25,outer=T)  
dev.off()



#Appendix 4. Cumulative catch and vessels
SPECIES.vec=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
id=70
Min.rec.ves=20  #minimum accepted number of records to keep a vessel


Data.Summary=vector('list',length=N.species)

for ( i in 1:N.species)Data.Summary[[i]]=fn.explore(DATA.list.LIVEWT.c[[i]],SPECIES.vec[i])
Data.Fig5=vector('list',length=N.species)
spe=c(17003,17001,18003,18007)
for ( i in 1:N.species)Data.Fig5[[i]]=fn.Figure5(Species.list[[i]],SPECIES.vec[i],spe[i])

tiff(file="Appendix 4.Cummulative.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
layout(matrix(c(1,2,1,2,3,4,3,4), 4, 2, byrow = TRUE))
par(mar=c(3,.75,.1,3),oma=c(1,4,.1,.01),las=1,mgp=c(1,.9,0))
line.col=c("black","grey20","grey55","grey80")
line.type=c(1,2,1,1)
plot(Data.Fig5[[1]]$Block.Cum.Ca[1:id],ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
for (i in 2:N.species)lines(Data.Fig5[[i]]$Block.Cum.Ca[1:id],col=line.col[i],lwd=3,lty=line.type[i])
axis(1,1:length(Data.Fig5[[1]]$Block.Cum.Ca[1:id]),labels=F,tck=-0.015)
axis(1,seq(5,length(Data.Fig5[[1]]$Block.Cum.Ca[1:id]),5),
     labels=F,tck=-0.03)
mtext("Cumulative catch (%)",side=2,line=2.75,font=1,las=0,cex=1.5,outer=F)

plot(Data.Fig5[[1]]$Ves.Cum.Ca,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
for (i in 2:N.species)lines(Data.Fig5[[i]]$Ves.Cum.Ca,col=line.col[i],lwd=3,lty=line.type[i])
axis(1,seq(0,length(Data.Fig5[[1]]$Ves.Cum.Ca),10),labels=F,tck=-0.015)
axis(1,seq(100,length(Data.Fig5[[1]]$Ves.Cum.Ca),100),labels=F,tck=-0.03)
legend("bottomright",SPECIES.vec,bty='n',cex=1.5,col=line.col,lwd=3,lty=line.type)

plot(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]*.99,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
for (i in 2:N.species)lines(Data.Fig5[[i]]$TABLE17$PerCumRecords[1:id],col=line.col[i],lwd=3,lty=line.type[i])
axis(1,1:length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),labels=F,tck=-0.015)
axis(1,seq(5,length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),5),
     labels=seq(5,length(Data.Fig5[[1]]$TABLE17$PerCumRecords[1:id]),5),tck=-0.03,cex.axis=1.65)
mtext("Number of blocks",side=1,line=2.5,font=1,las=0,cex=1.5,outer=F)
mtext("Cumulative records (%)",side=2,line=2.75,font=1,las=0,cex=1.5,outer=F)


plot(Data.Fig5[[1]]$TABLE16$PerCumRecords*.99,ylab="", xlab="",xaxt='n',type='l',cex.axis=1.65,lwd=3)
for (i in 2:N.species)lines(Data.Fig5[[i]]$TABLE16$PerCumRecords,col=line.col[i],lwd=3,lty=line.type[i])
axis(1,seq(0,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),10),labels=F,tck=-0.015)
axis(1,seq(100,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),100),
     labels=seq(100,length(Data.Fig5[[1]]$TABLE16$PerCumRecords),100),tck=-0.03,cex.axis=1.65)
mtext("Number of vessels",side=1,line=2.5,font=1,las=0,cex=1.5,outer=F)


#Plot of records per vessel
par(fig=c(0.60,1.00,.1,0.35), new = T,mgp=c(.25,.2,0),las=1)
fun.rec.per.ves(Data.Fig5[[1]]$Mean_catch_Vessel)
dev.off()


#Plot effective area
Dusky=c(X1=South.WA.long[1],X2=Dusky.range[2],Y1=South.WA.lat[1],Y2=Dusky.range[1])
Sandbar=c(X1=South.WA.long[1],X2=Sandbar.range[2],Y1=South.WA.lat[1],Y2=Sandbar.range[1])
Whiskery=c(X1=South.WA.long[1],X2=Whiskery.range[2],Y1=South.WA.lat[1],Y2=Whiskery.range[1])
Gummy=c(X1=Gummy.range[1],X2=Gummy.range[2],Y1=South.WA.lat[1],Y2=-31.6)
LISta=list(Whiskery=Whiskery,Gummy=Gummy,Dusky=Dusky,Sandbar=Sandbar)
fn.show=function(X1,X2,Y1,Y2)
{
  plotmap(a,b,PLATE,"grey85",South.WA.long,South.WA.lat)
  polygon(c(X1,X2,X2,X1),c(Y2,Y2,Y1,Y1),col='grey35',border="transparent")
  par(new=T)
  plotmap(a,b,PLATE,"grey85",South.WA.long,South.WA.lat)
}
jpeg(file="Effective_area.jpeg",width = 2400, height = 2400,units = "px", res = 300)
par(mfcol=c(2,2),mai=c(.5,.5,.2,.1),oma=c(.1,.1,.65,.1))
for(x in 1:length(LISta))
{
  fn.show(X1=LISta[[x]][[1]],X2=LISta[[x]][[2]],Y1=LISta[[x]][[3]],Y2=LISta[[x]][[4]]) 
  mtext(names(LISta)[x],3,cex=2)
  At=c(LISta[[x]][[1]],LISta[[x]][[2]])
  axis(side = 1, at =At, labels = At, tcl = .35,las=1,cex.axis=1,padj=-0.9)
  At=c(LISta[[x]][[3]],LISta[[x]][[4]])
  axis(side = 2, at =At , labels = -At,tcl = .35,las=2,cex.axis=1,hadj=.65)
}
dev.off()  



Check.num.vs.weight=FALSE
if(Check.num.vs.weight)
{
  D=Data.daily.GN%>%filter(!is.na(nfish) & SPECIES<35000)%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(SPECIES,year)%>%
    summarise(N.tot=sum(nfish),
              W.tot=sum(LIVEWT.c/1000),
              N.mean=mean(nfish),
              W.mean=mean(LIVEWT.c),
              Mean.weight=mean(LIVEWT.c/nfish),
              Mean.weight.sd=sd(LIVEWT.c/nfish),
              n = length(LIVEWT.c))%>%
    mutate(se=Mean.weight.sd/sqrt(n))
  A1=table(D$SPECIES)
  SP=names(A1[A1>10])
  D=subset(D,SPECIES%in%SP)
  
  Nms=Data.daily.GN%>%distinct(SPECIES,.keep_all = TRUE)%>%dplyr::select(c(SPECIES,SNAME))%>%
    filter(SPECIES%in%SP)%>%arrange(SPECIES)
  smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
  
  tiff(file="total_number.vs.weight.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(length(SP),MAR=c(1.5,4,1.5,2),OMA=c(3,3,.5,3),MGP=c(1,.5,0))
  for(s in 1:length(SP))
  {
    with(subset(D,SPECIES==SP[s]),{
      plot(year,N.tot,type='l',lwd=2,col=2,ylab='',xlab='',main=tolower(Nms$SNAME[s]))
      par(new=T)
      plot(year,W.tot,type='l',lwd=2,col=3,ylab='',xlab='',yaxt='n',xaxt='n')
      axis(4,seq(round(min(W.tot)),round(max(W.tot)),length.out = 5))
    })
  }
  mtext("Financial year",1,outer=T,line=.5,cex=1.2)
  mtext("total numbers",2,outer=T,line=1.25,cex=1.2,las=3)
  mtext("total weight (tonnes)",4,outer=T,line=1.25,cex=1.2,las=3)
  plot.new()
  legend("top",c("Number","Weight"),col=2:3,lty=1,bty='n',lwd=2,cex=1.25)
  dev.off()
  
  
  tiff(file="mean_number.vs.weight.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(length(SP),MAR=c(1.5,4,1.5,2),OMA=c(3,3,.5,3),MGP=c(1,.5,0))
  for(s in 1:length(SP))
  {
    with(subset(D,SPECIES==SP[s]),{
      plot(year,N.mean,type='l',lwd=2,col=2,ylab='',xlab='',main=tolower(Nms$SNAME[s]))
      par(new=T)
      plot(year,W.mean,type='l',lwd=2,col=3,ylab='',xlab='',yaxt='n',xaxt='n')
      axis(4,seq(round(min(W.mean)),round(max(W.mean)),length.out = 5))
    })
  }
  mtext("Financial year",1,outer=T,line=.5,cex=1.2)
  mtext("Mean numbers",2,outer=T,line=1.25,cex=1.2,las=3)
  mtext("Mean weight (kg)",4,outer=T,line=1.25,cex=1.2,las=3)
  plot.new()
  legend("top",c("Number","Weight"),col=2:3,lty=1,bty='n',lwd=2,cex=1.25)
  dev.off()
  
  
  tiff(file="mean_size.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(length(SP),MAR=c(1.5,4,1.5,2),OMA=c(3,3,.5,3),MGP=c(1,.5,0))
  for(s in 1:length(SP))
  {
    with(subset(D,SPECIES==SP[s]),{
      plot(year,Mean.weight,type='l',lwd=2,ylab='',xlab='',ylim=c(min(Mean.weight-1.96*se),max(Mean.weight+1.96*se)),main=tolower(Nms$SNAME[s]))
      segments(year,Mean.weight+1.96*se,year,Mean.weight-1.96*se)
    })
  }
  mtext("Financial year",1,outer=T,line=.5,cex=1.2)
  mtext("Mean size (kg)",2,outer=T,line=1.25,cex=1.2,las=3)
  dev.off()
}
