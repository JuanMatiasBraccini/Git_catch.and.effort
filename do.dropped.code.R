#3.1. Number of blocks and vessels per Yr.Mn
# Expand.fun.Yr.Mn=function(DATA)
# {
#   DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   Tab=table(DATA$Yr.Mn,DATA$BLOCKX)
#   Tab=ifelse(Tab>=1,1,0)
#   Tab=rowSums(Tab)
#   Yr.Mn=names(Tab)
#   Yr=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 1))
#   Mn=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 2))
#   Tab=data.frame(Tab,Yr,Mn)
#   Tab=Tab[order(Tab$Yr,Tab$Mn),1]
# 
#   return(Tab)
# }
# Spatial.expan.Yr.Mn=Expand.fun.Yr.Mn(Data.monthly)
# 
# Effort1.fun.YrMn=function(DATA)
# {
#   DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   Tab=table(DATA$Yr.Mn,DATA$VESSEL)
#   Tab=ifelse(Tab>=1,1,0)
#   Tab=rowSums(Tab)
#   Yr.Mn=names(Tab)
#   Yr=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 1))
#   Mn=as.numeric(sapply(strsplit(Yr.Mn," "), "[", 2))
#   Tab=data.frame(Tab,Yr,Mn)
#   Tab=Tab[order(Tab$Yr,Tab$Mn),1]
#   
#   return(Tab)
# }
# Effort.expan=Effort1.fun.YrMn(Data.monthly)
# 
# tiff(file="Figure 3. Folly and Fantasy.1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# par(mfcol=c(1,1),mar=c(1,3.6,.1,4),oma=c(2.5,.1,.1,.1),las=1,mgp=c(.15,.7,0))
# 
# plot(1:N.Yrs.months,Spatial.expan.Yr.Mn,ylab="",xlab="",xaxt="n",las=2,type="l",pch=19,cex=1.5,ylim=c(0,50)
#      ,cex.axis=1.25,lwd=1.75)
# axis(1,at=seq(1,N.Yrs.months,12),labels=F,tck=-0.01)
# 
# par(new=T)
# plot(1:N.Yrs.months,Effort.expan,col="grey60",type="l",pch=19,axes=F,ann='F',cex=1.5,lwd=1.75,ylim=c(0,105))
# axis(4,at=pretty(Effort.expan),labels=pretty(Effort.expan),las=2,cex.axis=1.3)
# 
# 
# axis(1,at=seq(1,N.Yrs.months,24),labels=names(Yrs.months)[seq(1,N.Yrs.months,24)],tck=-0.02,cex.axis=1.25)
# 
# 
# mtext("Year",side=1,line=1.5,font=1,las=0,cex=1.7,outer=T)
# mtext("Number of blocks fished",side=2,line=-1.2,font=1,las=0,cex=1.7,outer=T)
# mtext("Number of licence holders fishing",side=4,line=-1.2,las=3,cex=1.7,col="grey60",outer=T)
# dev.off()



# Mean.fun.Yr.Mn=function(VAR)
# {
#   DATA=Data.monthly.GN
#   
#   DATA$CPUE=with(DATA,ifelse(BDAYS.c*HOURS.c*SHOTS.c>0,LIVEWT.c/((NETLEN.c/1000)*BDAYS.c*HOURS.c*SHOTS.c*Inc.per),NA))
#   
#   #Yr.Mn
#   #DATA$Yr.Mn=paste(DATA$YEAR.c,DATA$MONTH)
#   #Tab=aggregate(CPUE~Yr.Mn+BLOCKX+SPECIES,data=DATA,mean,na.rm=T)
#   
#   #Yr only
#   Tab=aggregate(CPUE~YEAR.c+BLOCKX+SPECIES,data=DATA,mean,na.rm=T)
#   
#   
#   Tab1=subset(Tab,SPECIES==VAR)
# #  Tab1$Yr=as.numeric(sapply(strsplit(Tab1$Yr.Mn," "), "[", 1))
# #  Tab1$Mn=as.numeric(sapply(strsplit(Tab1$Yr.Mn," "), "[", 2))
# #  Tab1=Tab1[order(Tab1$Yr,Tab1$Mn),]
# #  Tab1=Tab1[,c(2,4:6)]
# 
# #  Reshaped=as.matrix(reshape(Tab1,idvar=c("Yr","Mn"),  	#transposed as matrix 	
# #                             timevar="BLOCKX",v.names="CPUE", direction="wide"))	
#   
#   Tab1=Tab1[,c(1:2,4)]
#   
#   Reshaped=as.matrix(reshape(Tab1,idvar=c("YEAR.c"),    #transposed as matrix 	
#                              timevar="BLOCKX",v.names="CPUE", direction="wide"))
#   Reshaped=Reshaped[order(Reshaped[,1]),]
#   return(Reshaped)
# }
#  Spatial.expan.Yr.Mn=Mean.fun.Yr.Mn(TARGETS[2])
# # 
# # tiff(file="Figure 3. Folly and Fantasy.3.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
# # 
# # plot(1:N.Yrs.months,Spatial.expan.Yr.Mn[,3],ylim=c(0,10))
# # for(i in 4:ncol(Spatial.expan.Yr.Mn))   lines(1:N.Yrs.months,Spatial.expan.Yr.Mn[,i],col=runif(1,1,100))
# 
#  plot(1:nrow(Spatial.expan.Yr.Mn),Spatial.expan.Yr.Mn[,2],ylim=c(0,20),ann=F,xaxt='n',col='transparent')
#  for(i in 2:ncol(Spatial.expan.Yr.Mn))   lines(1:nrow(Spatial.expan.Yr.Mn),Spatial.expan.Yr.Mn[,i],col=runif(1,1,100))
# 

#dev.off()


# fun.prop=function(DAT,SPEC)
# {
#   #Vessel, gear, fin. year, month, block (given by the "Same.return" variable)
#   ID=which(DAT$SPECIES==SPEC)
#   this.same.returns=unique(DAT[ID,]$Same.return)
#   dat=subset(DAT,Same.return %in% this.same.returns & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~Same.return,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~Same.return,data=dat,sum,na.rm=T)
#   Prop.VesYrMonBlock=data.frame(Same.return=All.1$Same.return,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #   #Vessel, fin. year (mean proportion given by the "AnnualVesselAveID" variable)
#   #   this.same=unique(DAT[ID,]$AnnualVesselAveID)
#   #   dat=subset(DAT,AnnualVesselAveID %in% this.same & SPECIES %in% Shark.species)
#   #   dat.species=subset(dat,SPECIES==SPEC)
#   #   Target.sp.1=aggregate(LIVEWT~AnnualVesselAveID,data=dat.species,sum,na.rm=T)
#   #   All.1=aggregate(LIVEWT~AnnualVesselAveID,data=dat,sum,na.rm=T)
#   #   Prop.VesFinYr=data.frame(AnnualVesselAveID=All.1$AnnualVesselAveID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, month (mean proportion given by the "MonthlyID" variable)
#   this.same=unique(DAT[ID,]$MonthlyID)
#   dat=subset(DAT,MonthlyID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~MonthlyID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~MonthlyID,data=dat,sum,na.rm=T)
#   Prop.FinYrMon=data.frame(MonthlyID=All.1$MonthlyID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #   #Fin. year, block (mean proportion given by the "BlockID" variable) 
#   #   this.same=unique(DAT[ID,]$BlockID)
#   #   dat=subset(DAT,BlockID %in% this.same & SPECIES %in% Shark.species)
#   #   dat.species=subset(dat,SPECIES==SPEC)
#   #   Target.sp.1=aggregate(LIVEWT~BlockID,data=dat.species,sum,na.rm=T)
#   #   All.1=aggregate(LIVEWT~BlockID,data=dat,sum,na.rm=T)
#   #   Prop.FinYrBlok=data.frame(BlockID=All.1$BlockID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, month, block (mean proportion given by the "GoodsplitID" variable)
#   this.same=unique(DAT[ID,]$GoodsplitID)
#   dat=subset(DAT,GoodsplitID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~GoodsplitID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~GoodsplitID,data=dat,sum,na.rm=T)
#   Prop.GoodsplitID=data.frame(GoodsplitID=All.1$GoodsplitID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #Fin. year, zone (given by the "ZoneID" variable)  (mean proportion)
#   this.same=unique(DAT[ID,]$ZoneID)
#   dat=subset(DAT,ZoneID %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~ZoneID,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~ZoneID,data=dat,sum,na.rm=T)
#   Prop.FinYrZone=data.frame(ZoneID=All.1$ZoneID,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   #zone (given by the "zone" variable)  (mean proportion)
#   this.same=unique(DAT[ID,]$zone)
#   dat=subset(DAT,zone %in% this.same & SPECIES %in% Shark.species)
#   dat.species=subset(dat,SPECIES==SPEC)
#   Target.sp.1=aggregate(LIVEWT~zone,data=dat.species,sum,na.rm=T)
#   All.1=aggregate(LIVEWT~zone,data=dat,sum,na.rm=T)
#   Prop.Zone=data.frame(zone=All.1$zone,Proportion=Target.sp.1$LIVEWT/All.1$LIVEWT)
#   
#   return(list(Prop.VesYrMonBlock=Prop.VesYrMonBlock,Prop.GoodsplitID=Prop.GoodsplitID,
#               Prop.FinYrZone=Prop.FinYrZone,Prop.FinYrMon=Prop.FinYrMon,Prop.Zone=Prop.Zone))
# }

# Catch.prop.gummy=fun.prop(Data.monthly,17001)
# Catch.prop.whiskery=fun.prop(Data.monthly,17003)
# Catch.prop.dusky=fun.prop(Data.monthly,18003)
# Catch.prop.sandbar=fun.prop(Data.monthly,18007)
# Catch.prop.school=fun.prop(Data.monthly,17008)
# Catch.prop.dogfish=fun.prop(Data.monthly,20000)
# Catch.prop.other=fun.prop(Data.monthly,Sharks.other)


# #create bad reporter files for fixing catches
# Bad.Reporters=subset(Data.monthly,Reporter=="bad")
# 
# Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter))
# 
# Bad.dus.gum.whi=subset(Bad.Reporters,Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0)
# NroW=nrow(Bad.dus.gum.whi)
# 
# Bad.dus.gum.whi.noBMY=subset(Bad.Reporters,!(Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0))
# 
# 
# #Replicate Bad.dus.gum.whi twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi=rbind(Bad.dus.gum.whi,Bad.dus.gum.whi,Bad.dus.gum.whi)
# Bad.dus.gum.whi=Bad.dus.gum.whi[order(Bad.dus.gum.whi$Same.return),]
# 
# Bad.dus.gum.whi$Spec.old=Bad.dus.gum.whi$SPECIES
# Bad.dus.gum.whi$Sname.old=Bad.dus.gum.whi$SNAME
# 
# Bad.dus.gum.whi$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.dus.gum.whi$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
# Bad.dus.gum.whi$LIVEWT.reap=with(Bad.dus.gum.whi,
#                                  ifelse(SPECIES%in%c(18003),Shark.other.livewt*Prop.Dus.Good.spl,
#                                         ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Good.spl,
#                                                ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi=subset(Bad.dus.gum.whi,LIVEWT.reap>0)
# 
# #create new vars
# Bad.dus.gum.whi$Reporter.old=Bad.dus.gum.whi$Reporter
# Bad.dus.gum.whi$Reporter="good"
# 
# 
# #add old species column to data
# Data.monthly$Spec.old=Data.monthly$SPECIES
# Data.monthly$Sname.old=Data.monthly$SNAME
# Data.monthly$Reporter.old=Data.monthly$Reporter
# 
# 
# #update "bad" recorders with reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.8  Update good split ave zone catches
# #If valid month-year-block proportions NOT available, then update "bad" records of Dusky, Gummy and
# #     whiskery with year-month-zone average
# 
# Data.monthly=merge(Data.monthly,Zone.good.split,by="ZoneID",all.x=T)
# Bad.dus.gum.whi.noBMY=merge(Bad.dus.gum.whi.noBMY,Zone.good.split,by="ZoneID",all.x=T)
# 
# Bad.dus.gum.whi.noBMY.month=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999 &
#                                      !(Prop.Dus.Zone.Good.spl>0 & Prop.Gum.Zone.Good.spl>0 & Prop.Whi.Zone.Good.spl>0))
# 
# 
# NroW.noBMY=nrow(Bad.dus.gum.whi.noBMY)
# 
# #Replicated Bad.dus.gum.whi.noBMY twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi.noBMY=rbind(Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY)
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[order(Bad.dus.gum.whi.noBMY$Same.return),]
# 
# Bad.dus.gum.whi.noBMY$Spec.old=Bad.dus.gum.whi.noBMY$SPECIES
# Bad.dus.gum.whi.noBMY$Sname.old=Bad.dus.gum.whi.noBMY$SNAME
# 
# 
# Bad.dus.gum.whi.noBMY$SPECIES=rep(c(18003,17001,17003),NroW.noBMY)
# Bad.dus.gum.whi.noBMY$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY)
# 
# 
# Bad.dus.gum.whi.noBMY$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY,
#                                        ifelse(SPECIES%in%c(18003),Shark.other.livewt*Prop.Dus.Zone.Good.spl,
#                                               ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Zone.Good.spl,
#                                                      ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Zone.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi.noBMY=subset(Bad.dus.gum.whi.noBMY,LIVEWT.reap>0)
# 
# Bad.dus.gum.whi.noBMY$Reporter.old=Bad.dus.gum.whi.noBMY$Reporter
# Bad.dus.gum.whi.noBMY$Reporter="good"
# 
# 
# #update "bad" recorders with reapportioned catch
# ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY))
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[,ID.names]
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.9  Update good split ave monthly catches           #Rory's rules 4f                          
# #If valid month-year-block proportions NOT available or month-year-block proportions NOT available,
# #       then update "bad" records of Dusky, Gummy and whiskery with month-year average
# 
# Data.monthly=merge(Data.monthly,Monthly.good.split,by="MonthlyID",all.x=T)
# 
# if(nrow(Bad.dus.gum.whi.noBMY.month)>0)
# {
#   Bad.dus.gum.whi.noBMY.month=merge(Bad.dus.gum.whi.noBMY.month,Monthly.good.split,by="MonthlyID",all.x=T)
#   
#   NroW.noBMY.month=nrow(Bad.dus.gum.whi.noBMY.month)
#   
#   #Replicated Bad.dus.gum.whi.noBMY.month twice to add catch of dusky, gummy and whiskery
#   Bad.dus.gum.whi.noBMY.month=rbind(Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[order(Bad.dus.gum.whi.noBMY.month$Same.return),]
#   
#   Bad.dus.gum.whi.noBMY.month$Spec.old=Bad.dus.gum.whi.noBMY.month$SPECIES
#   Bad.dus.gum.whi.noBMY.month$Sname.old=Bad.dus.gum.whi.noBMY.month$SNAME
#   
#   
#   Bad.dus.gum.whi.noBMY.month$SPECIES=rep(c(18003,17001,17003),NroW.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY.month)
#   
#   Bad.dus.gum.whi.noBMY.month$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY.month,
#                                                ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Mon.Good.spl,
#                                                       ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Mon.Good.spl,
#                                                              ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Mon.Good.spl,LIVEWT.reap))))
#   
#   #remove artificially created 0 catches
#   Bad.dus.gum.whi.noBMY.month=subset(Bad.dus.gum.whi.noBMY.month,LIVEWT.reap>0)
#   
#   Bad.dus.gum.whi.noBMY.month$Reporter.old=Bad.dus.gum.whi.noBMY.month$Reporter
#   Bad.dus.gum.whi.noBMY.month$Reporter="good"
#   
#   
#   #update "bad" recorders with reapportioned catch
#   ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY.month))
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[,ID.names]
#   Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY.month)
#   
#   #remove duplicates
#   Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
#   Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
#   
# }
# 
# # 
# #C.7.15 Reapportion catch                                     #Rory's rules 6k-6s            
# #note: uses same rules as for southern catch (#C7.7- #C7.9)
# 
# #create bad reporter files for fixing catches
# Bad.Reporters=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999)
# Data.monthly=subset(Data.monthly,Reporter=="good"|is.na(Reporter)| !(Reporter=="bad" & SPECIES== 22999))
# 
# 
# Bad.dus.gum.whi=subset(Bad.Reporters,Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0)
# NroW=nrow(Bad.dus.gum.whi)
# 
# Bad.dus.gum.whi.noBMY=subset(Bad.Reporters,!(Prop.Dus.Good.spl>0 | Prop.Gum.Good.spl>0 | Prop.Whi.Good.spl>0))
# 
# #C.7.15.1 Good.spl criteria
# #Replicated Bad.dus.gum.whi twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi=rbind(Bad.dus.gum.whi,Bad.dus.gum.whi,Bad.dus.gum.whi)
# Bad.dus.gum.whi=Bad.dus.gum.whi[order(Bad.dus.gum.whi$Same.return),]
# 
# Bad.dus.gum.whi$Spec.old=Bad.dus.gum.whi$SPECIES
# Bad.dus.gum.whi$Sname.old=Bad.dus.gum.whi$SNAME
# 
# Bad.dus.gum.whi$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.dus.gum.whi$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
# 
# Bad.dus.gum.whi$LIVEWT.reap=with(Bad.dus.gum.whi,
#                                  ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Good.spl,
#                                         ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Good.spl,
#                                                ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi=subset(Bad.dus.gum.whi,LIVEWT.reap>0)
# 
# #create new vars
# Bad.dus.gum.whi$Reporter.old=Bad.dus.gum.whi$Reporter
# Bad.dus.gum.whi$Reporter="good"
# 
# #update "bad" recorders with reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.15.2 ZoneID criteria
# NroW.noBMY=nrow(Bad.dus.gum.whi.noBMY)
# 
# #Replicated Bad.dus.gum.whi.noBMY twice to add catch of dusky, gummy and whiskery
# Bad.dus.gum.whi.noBMY=rbind(Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY,Bad.dus.gum.whi.noBMY)
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[order(Bad.dus.gum.whi.noBMY$Same.return),]
# 
# 
# Bad.dus.gum.whi.noBMY$SPECIES=rep(c(18003,17001,17003),NroW.noBMY)
# Bad.dus.gum.whi.noBMY$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY)
# 
# 
# Bad.dus.gum.whi.noBMY$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY,
#                                        ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Zone.Good.spl,
#                                               ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Zone.Good.spl,
#                                                      ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Zone.Good.spl,LIVEWT.reap))))
# 
# #remove artificially created 0 catches
# Bad.dus.gum.whi.noBMY=subset(Bad.dus.gum.whi.noBMY,LIVEWT.reap>0)
# 
# Bad.dus.gum.whi.noBMY$Reporter.old=Bad.dus.gum.whi.noBMY$Reporter
# Bad.dus.gum.whi.noBMY$Reporter="good"
# 
# 
# #update "bad" recorders with reapportioned catch
# ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY))
# Bad.dus.gum.whi.noBMY=Bad.dus.gum.whi.noBMY[,ID.names]
# Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY)
# 
# #remove duplicates
# Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
# Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
# 
# 
# 
# #C.7.15.3 Ave monthly catches criteria                               
# 
# Bad.dus.gum.whi.noBMY.month=subset(Data.monthly,Reporter=="bad" & SPECIES== 22999)
# 
# if(nrow(Bad.dus.gum.whi.noBMY.month)>0)
# {
#   NroW.noBMY.month=nrow(Bad.dus.gum.whi.noBMY.month)
#   
#   #Replicated Bad.dus.gum.whi.noBMY.month twice to add catch of dusky, gummy and whiskery
#   Bad.dus.gum.whi.noBMY.month=rbind(Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month,Bad.dus.gum.whi.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[order(Bad.dus.gum.whi.noBMY.month$Same.return),]
#   
#   
#   Bad.dus.gum.whi.noBMY.month$SPECIES=rep(c(18003,17001,17003),NroW.noBMY.month)
#   Bad.dus.gum.whi.noBMY.month$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW.noBMY.month)
#   
#   Bad.dus.gum.whi.noBMY.month$LIVEWT.reap=with(Bad.dus.gum.whi.noBMY.month,
#                                                ifelse(SPECIES==18003,Shark.other.livewt*Prop.Dus.Mon.Good.spl,
#                                                       ifelse(SPECIES==17001,Shark.other.livewt*Prop.Gum.Mon.Good.spl,
#                                                              ifelse(SPECIES==17003,Shark.other.livewt*Prop.Whi.Mon.Good.spl,LIVEWT.reap))))
#   
#   #remove artificially created 0 catches
#   Bad.dus.gum.whi.noBMY.month=subset(Bad.dus.gum.whi.noBMY.month,LIVEWT.reap>0)
#   
#   Bad.dus.gum.whi.noBMY.month$Reporter.old=Bad.dus.gum.whi.noBMY.month$Reporter
#   Bad.dus.gum.whi.noBMY.month$Reporter="good"
#   
#   
#   #update "bad" recorders with reapportioned catch
#   ID.names=match(names(Data.monthly),names(Bad.dus.gum.whi.noBMY.month))
#   Bad.dus.gum.whi.noBMY.month=Bad.dus.gum.whi.noBMY.month[,ID.names]
#   Data.monthly=rbind(Data.monthly,Bad.dus.gum.whi.noBMY.month)
#   
#   #remove duplicates
#   Data.monthly$Dupli=with(Data.monthly,paste(VesselID,METHOD,SPECIES,LIVEWT.reap))
#   Data.monthly=Data.monthly[!duplicated(Data.monthly$Dupli),]
#   
# }
# 



# #Compare mine and Rory's catch and effort
# setwd("C:/Users/myb/Desktop/New folder")
# 
# Use.Previos.Sofar="YES"   #Select YES if attaching previous Sofar data to current year
# #Use.Previos.Sofar="NO"
# 
# 
# Ind.spe.list=list(Gummy=17001,Whiskery=17003,Bronzy.Dusky=c(18001,18003),sandbar=18007)
# 
# #Fishing effort limits
# FishEffLims=data.frame(zone=c("West","Zone1","Zone2"),Km.Gillnet.Hours.c=c(67692,84075,144102),
#                        Km.Gillnet.Days.c=c(2832,3503,7205))
# 
# #Current year data set
# DAT=subset(Data.monthly,FINYEAR==Current.yr)
# 
# #Main Feature table (total catch, indicator species catch, teleost catch, etc)
# Other=subset(Other.fishery.catch,financial.year==Current.yr)        
# 
# 
# #Add Fishing effort                                       #REVIEW RORY
# 
# C.yr=match(Current.yr,Total.effort.days.monthly$FINYEAR)
# 
# #annual
# 
# Curr.annual.1000km.gn.hours=Total.effort.hours.monthly[C.yr,]
# Curr.annual.km.gn.days=Total.effort.days.monthly[C.yr,]
# Curr.annual.km.gn.days[1,2]=Curr.annual.km.gn.days[1,2]*1000
# 
# #annual by zone
# Curr.annual.1000km.gn.hours.zone=Total.effort.zone.hours.monthly[C.yr,]
# Curr.annual.km.gn.days.zone=Total.effort.zone.days.monthly[C.yr,]
# Curr.annual.km.gn.days.zone[,2:4]=Curr.annual.km.gn.days.zone[,2:4]*1000
# 
# #percentages of effort limits
# #annual
# Per.annual.lim.1000km.gn.hours=100*Curr.annual.1000km.gn.hours$Total/
#   (sum(FishEffLims$Km.Gillnet.Hours.c)/1000)
# Per.annual.lim.km.gn.days=100*Curr.annual.km.gn.days$Total/sum(FishEffLims$Km.Gillnet.Days.c)
# 
# #annual by zone
# Per.annual.lim.1000km.gn.hours.zone=100*Curr.annual.1000km.gn.hours.zone[,2:4]/
#   (FishEffLims$Km.Gillnet.Hours.c/1000)
# Per.annual.lim.km.gn.days.zone=100*Curr.annual.km.gn.days.zone[,2:4]/FishEffLims$Km.Gillnet.Days.c
# 
# 
# 
# 
# #Figures 2-3. 
# par(mfcol=c(1,1),mar=c(3.5,3.6,.1,1),oma=c(1,.5,.1,.1))
# LINE=c(5,1,1)
# TYPE=c("l","o","o")
# PCH=21
# COL=1
# BG=c("black","black","white")
# 
# #Redefine DAT for "Use.Previos.Sofar=="YES"
# DAT=subset(Data.monthly,FINYEAR%in%c("2011-12",Current.yr) & METHOD%in%c("GN","LL"))
# 
# 
# 
# if(Use.Previos.Sofar=="YES")
# {
#   
#   fun.fig.all=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT[,2],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     plot(1:NN,DAT[id:N,2]/scaler,type='l',col="grey80",ylim=c(0,MAX),xaxt='n',yaxt='n',
#          ylab=TITLE1, xlab=TITLE2,las=1,lwd=2,cex.lab=1.3)
#     axis(1,at=1:NN,labels=F,tck=-0.01)
#     axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
#     
#     axis(2,at=seq(0,MAX,INT),labels=F,tck=-0.01)
#     axis(2,at=seq(0,MAX,INT2),labels=seq(0,MAX,INT2),tck=-0.02,cex.axis=1.1)
#     
#     # for(i in 1:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COL,lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   fun.fig.zn=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT1[,2:4],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     plot(1:NN,DAT1[id:N,2]/scaler,type='l',col=COL,ylim=c(0,MAX),xaxt='n',yaxt='n',
#          ylab=TITLE1, xlab=TITLE2,las=1,lwd=1,cex.lab=1.3,lty=LINE[1])
#     axis(1,at=1:NN,labels=F,tck=-0.01)
#     axis(1,at=seq(1,NN,5),labels=FInYEAR[seq(1,NN,5)],tck=-0.02,cex.axis=1.1)
#     
#     axis(2,at=seq(0,MAX,INT),labels=F,tck=-0.01)
#     axis(2,at=seq(0,MAX,INT2),labels=seq(0,MAX,INT2),tck=-0.02,cex.axis=1.1)
#     
#     for(i in 2:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COL,lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   fn.figs2.3.SoFaR.all=function(GROUP,LAT1,LAT2,INT,INT2,what)
#   {
#     dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     #add previous years
#     if(what=="Elasmos")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.sk.live.wt","Z1.tot.sk.live.wt","Z2.tot.sk.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.sk.live.wt"),names(Results.pre.2013))]
#     }
#     
#     if(what=="Teleosts")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.tel.live.wt","Z1.tot.tel.live.wt","Z2.tot.tel.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.tel.live.wt"),names(Results.pre.2013))]
#     }
#     
#     names(Prev.zn)=names(wide)
#     Prev.zn[,2:4]=Prev.zn[,2:4]*1000
#     wide=rbind(Prev.zn,wide)
#     
#     
#     Prev[,2]=Prev[,2]*1000
#     names(Prev)=names(annual.catch.total)
#     annual.catch.total=rbind(Prev,annual.catch.total)
#     
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.all(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
#   
#   fn.figs2.3.SoFaR.zn=function(GROUP,LAT1,LAT2,INT,INT2,what)
#   {
#     dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     #add previous years
#     if(what=="Elasmos")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.sk.live.wt","Z1.tot.sk.live.wt","Z2.tot.sk.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.sk.live.wt"),names(Results.pre.2013))]
#     }
#     
#     if(what=="Teleosts")
#     {
#       Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.tot.tel.live.wt","Z1.tot.tel.live.wt","Z2.tot.tel.live.wt"),
#                                       names(Results.pre.2013))]
#       Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.tot.tel.live.wt"),names(Results.pre.2013))]
#     }
#     
#     names(Prev.zn)=names(wide)
#     Prev.zn[,2:4]=Prev.zn[,2:4]*1000
#     wide=rbind(Prev.zn,wide)
#     
#     
#     Prev[,2]=Prev[,2]*1000
#     names(Prev)=names(annual.catch.total)
#     annual.catch.total=rbind(Prev,annual.catch.total)
#     
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.zn(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
# }
# 
# if(Use.Previos.Sofar=="NO")
# {
#   fun.fig.all=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT[,2],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     par(mgp=c(2.5,.65,0),las=1)
#     lines(1:NN,DAT[id:N,2]/scaler,type='l',col="red",lwd=2)
#   }
#   
#   fun.fig.zn=function(DAT,DAT1,scaler,TITLE1,TITLE2,INT,INT2)
#   {
#     MAX=max(DAT1[,2:4],na.rm=T)/scaler
#     FInYEAR=as.character(unique(DAT$finyear))
#     N=length(FInYEAR)
#     
#     #id=match(start.yr,FInYEAR)
#     id=which.min(is.na(DAT[,2]))
#     FInYEAR=FInYEAR[id:length(FInYEAR)]
#     NN=length(FInYEAR)
#     
#     BG=c("red","green","blue") 
#     COLs=BG
#     for(i in 1:(ncol(DAT1)-1))points(1:NN,DAT1[id:N,i+1]/scaler,type=TYPE[i],lty=LINE[i],col=COLs[i],lwd=1.5,pch=PCH,bg=BG[i])
#   }
#   
#   
#   fn.figs2.3.SoFaR.all=function(GROUP,LAT1,LAT2,INT,INT2)
#   {
#     dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2 & METHOD%in%c("GN","LL"))
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     names(wide)[match("FINYEAR",names(wide))]="finyear"
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.all(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
#   
#   fn.figs2.3.SoFaR.zn=function(GROUP,LAT1,LAT2,INT,INT2)
#   {
#     dat=subset(Data.monthly,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
#     
#     annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
#     annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
#     
#     wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
#     
#     names(wide)[match("FINYEAR",names(wide))]="finyear"
#     names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
#     
#     fun.fig.zn(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
#   }
# }
# 
# #total
# jpeg(file="Figure 2.TotalElasmoCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fn.figs2.3.SoFaR.all(Elasmo.species,-27,-40,100,500,"Elasmos")
# if(Use.Previos.Sofar=="NO")fn.figs2.3.SoFaR.all(Elasmo.species,-27,-40,100,500)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("grey80","red"),lwd=2)
# dev.off()
# 
# #by zone
# jpeg(file="Figure 2.ZoneElasmoCatch.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fn.figs2.3.SoFaR.zn(Elasmo.species,-27,-40,100,500,"Elasmos")
# if(Use.Previos.Sofar=="NO")fn.figs2.3.SoFaR.zn(Elasmo.species,-27,-40,100,500)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("black","red"),lwd=2)
# legend("topleft",c("West","Zn1","Zn2"),bty='n',lty=LINE,col=c("red","green","blue"),lwd=2)
# dev.off()
# 
# 
# 
# 
# #Figure 4
# 
# #add 2011-12 to Current
# Current.yr=c("2011-12","2012-13")
# Eff.Current.Yr=fn.Eff.Sofar(Current.yr)
# Total.effort.zone.days=Eff.Current.Yr$Total.effort.zone.days
# Total.effort.zone.hours=Eff.Current.Yr$Total.effort.zone.hours
# Total.effort.joint.days=Eff.Current.Yr$Total.effort.joint.days
# Total.effort.joint.hours=Eff.Current.Yr$Total.effort.joint.hours
# Total.effort.days=Eff.Current.Yr$Total.effort.days
# Total.effort.hours=Eff.Current.Yr$Total.effort.hours
# 
# if(Use.Previos.Sofar=="YES")
# {
#   Prev.zn=Results.pre.2013[,match(c("FINYEAR","WC.km.gn.days","Z1.km.gn.days","Z2.km.gn.days"),
#                                   names(Results.pre.2013))]
#   Prev=Results.pre.2013[,match(c("FINYEAR","TDGDLF.km.gn.days"),names(Results.pre.2013))]
#   
#   idi=match(Current.yr,Total.effort.zone.days$FINYEAR)
#   wide=data.frame(finyear=Current.yr,WC.km.gn.days=Total.effort.zone.days[idi,2]*1000,
#                   Z1.km.gn.days=Total.effort.zone.days[idi,3]*1000,
#                   Z2.km.gn.days=Total.effort.zone.days[idi,4]*1000)
#   
#   names(Prev.zn)=names(wide)
#   wide=rbind(Prev.zn,wide)
#   
#   annual.effort.days.total=data.frame(finyear=Current.yr,TDGDLF.km.gn.days=Total.effort.days[idi,2]*1000)
#   names(Prev)=names(annual.effort.days.total)
#   annual.effort.days.total=rbind(Prev,annual.effort.days.total)
#   
#   
# }
# 
# if(Use.Previos.Sofar=="NO")
# {
#   annual.effort.days.total=Total.effort.days.monthly
#   wide=Total.effort.zone.days.monthly
#   names(wide)[match("FINYEAR",names(wide))]="finyear"
#   names(annual.effort.days.total)[match("FINYEAR",names(annual.effort.days.total))]="finyear"
# }
# 
# #all
# jpeg(file="Figure 4.StandardisedEffort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fun.fig.all(annual.effort.days.total,wide,1000,"Effort (1000km gn.d)","Financial year",2,10)
# if(Use.Previos.Sofar=="NO")fun.fig.all(annual.effort.days.total,wide,1,"Effort (1000km gn.d)","Financial year",2,10)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("grey80","red"),lwd=2)
# dev.off()
# 
# #by zone
# jpeg(file="Figure 4.ZoneStandardisedEffort.jpeg",width = 2400, height = 2400,units = "px", res = 300)
# if(Use.Previos.Sofar=="YES")fun.fig.zn(annual.effort.days.total,wide,1000,"Effort (1000km gn.d)","Financial year",2,10)
# if(Use.Previos.Sofar=="NO")fun.fig.zn(annual.effort.days.total,wide,1,"Effort (1000km gn.d)","Financial year",2,10)
# legend("bottomright",c("Previous","Mine"),bty='n',lty=1,col=c("black","red"),lwd=2)
# legend("topleft",c("West","Zn1","Zn2"),bty='n',lty=LINE,col=c("red","green","blue"),lwd=2)
# 
# dev.off()




#REMOVED REAPPORTIONING CODE
# 
#   #C.7.8.2 First fix Bad.dusky, Bad.gummy, Bad.whiskery                         #REVIEW RORY
# Bad.Reporters.Dus.Gum.Whi=subset(Bad.Reporters,SPECIES%in%c(18003,17001,17003))
# 
#     #reapportion catch
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=NA
# 
#         #first use "Good.spl" 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#         ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#         ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#         ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#         LIVEWT.reap))))
# 
#         #second, if previous not available, use "Zone.Good.spl" (i.e. Yr-Mn-Zone)
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#     ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #USE YR-MONTH AGAIN!!!!
#       #third, if previous are not available, use "Mon.Good.spl"  (i.e. Yr-Mn)            #Rory's rules 4f 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
#         (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
#          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
#         (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
#     LIVEWT.reap))))
# 
# #         #third, if previous are not available, use "YrZn.Good.spl" (i.e. Yr-Zone)         #REVIEW RORY 
# # Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
# #      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
# #     LIVEWT.reap))))
# 
# #create new vars
# Bad.Reporters.Dus.Gum.Whi$Reporter.old=Bad.Reporters.Dus.Gum.Whi$Reporter
# Bad.Reporters.Dus.Gum.Whi$Reporter="good"
# Bad.Reporters.Dus.Gum.Whi$Spec.old=Bad.Reporters.Dus.Gum.Whi$SPECIES
# Bad.Reporters.Dus.Gum.Whi$Sname.old=Bad.Reporters.Dus.Gum.Whi$SNAME
# 
# 
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters.Dus.Gum.Whi)
# 
# 
#   #C.7.8.3 Then fix 22999                               #REVIEW RORY
#     #note: remove duplicates of Same.return as Tot.shk.livewt is split proportionally
# #          among dusky, whiskery and gummy
# Bad.Reporters=subset(Bad.Reporters,SPECIES==22999)
# Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),]
# NroW=nrow(Bad.Reporters)
# 
#     # replicate Bad.Reporters twice to have the three species as a record
# Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters)
# Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]
# 
# Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
# Bad.Reporters$Sname.old=Bad.Reporters$SNAME
# 
# Bad.Reporters$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
#   # reapportion catch
# Bad.Reporters$LIVEWT.reap=NA
# 
# #first use "Good.spl" (standardise the proportions to sum(split catch)=Tot.shk.livewt)
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#      ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#      ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#      LIVEWT.reap))))
# 
# #Second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#     LIVEWT.reap))))
# 
# # #Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# # Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
# #     ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #         (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #     LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #create file for flagging bad reporters
# #Flag.bad.rep1=Bad.Reporters[,match(c("Same.return","Spec.old","Reporter"),names(Bad.Reporters))]
# 
# #create new vars
# Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
# Bad.Reporters$Reporter="good"
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters)

#   #C.7.16.2 First fix Bad.dusky, Bad.gummy, Bad.whiskery 
# Bad.Reporters.Dus.Gum.Whi=subset(Bad.Reporters,SPECIES%in%c(18003,17001,17003))
# 
# #reapportion catch
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=NA
# 
# #first use "Good.spl" 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#            ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,
#            (Tot.shk.livewt*Prop.Dus.Good.spl),
#            ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,
#            (Tot.shk.livewt*Prop.Gum.Good.spl),
#            ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,
#            (Tot.shk.livewt*Prop.Whi.Good.spl),
#             LIVEWT.reap))))
# 
# #second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#             ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,
#             (Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#             ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,
#            (Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#            ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,
#           (Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#           LIVEWT.reap))))
# 
# # #third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# # Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
# #                  ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #                 (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #                 ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #                (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #                ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #               (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #             LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters.Dus.Gum.Whi$LIVEWT.reap=with(Bad.Reporters.Dus.Gum.Whi,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# 
# #create new vars
# Bad.Reporters.Dus.Gum.Whi$Reporter.old=Bad.Reporters.Dus.Gum.Whi$Reporter
# Bad.Reporters.Dus.Gum.Whi$Reporter="good"
# Bad.Reporters.Dus.Gum.Whi$Spec.old=Bad.Reporters.Dus.Gum.Whi$SPECIES
# Bad.Reporters.Dus.Gum.Whi$Sname.old=Bad.Reporters.Dus.Gum.Whi$SNAME
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters.Dus.Gum.Whi)
# 
# 
# 
#   #C.7.16.3 Then fix 22999
#     # remove duplicates of Same.return as Tot.shk.livewt is split proportionally
# Bad.Reporters=subset(Bad.Reporters,SPECIES==22999)
# Bad.Reporters=Bad.Reporters[!duplicated(Bad.Reporters$Same.return),]
# NroW=nrow(Bad.Reporters)
# 
#     # replicate Bad.Reporters twice to have the three species as a record
# Bad.Reporters=rbind(Bad.Reporters,Bad.Reporters,Bad.Reporters)
# Bad.Reporters=Bad.Reporters[order(Bad.Reporters$Same.return),]
# 
# Bad.Reporters$Spec.old=Bad.Reporters$SPECIES
# Bad.Reporters$Sname.old=Bad.Reporters$SNAME
# 
# Bad.Reporters$SPECIES=rep(c(18003,17001,17003),NroW)
# Bad.Reporters$SNAME=rep(c("SHARK, BRONZE WHALER","SHARK, GUMMY","SHARK, WHISKERY"),NroW)
# 
#     # reapportion catch
# Bad.Reporters$LIVEWT.reap=NA
# 
# #first use "Good.spl" 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & Prop.Dus.Good.spl>0,(Tot.shk.livewt*Prop.Dus.Good.spl),
#      ifelse(SPECIES==17001 & Prop.Gum.Good.spl>0,(Tot.shk.livewt*Prop.Gum.Good.spl),
#      ifelse(SPECIES==17003 & Prop.Whi.Good.spl>0,(Tot.shk.livewt*Prop.Whi.Good.spl),
#      LIVEWT.reap))))
# 
# 
# #Second, if previous not available, use "Zone.Good.spl"
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Zone.Good.spl>0,
#         (Tot.shk.livewt*Prop.Dus.Zone.Good.spl),
#      ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Zone.Good.spl>0,
#          (Tot.shk.livewt*Prop.Gum.Zone.Good.spl),
#      ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Zone.Good.spl>0,
#          (Tot.shk.livewt*Prop.Whi.Zone.Good.spl),
#      LIVEWT.reap))))
# 
# #Third, if previous are not available, use "Mon.Good.spl"             #Rory's rules 4f 
# #note: this is nonsense because it aggregates across all zones when gummy and whiskery
# #       do not occur in the north so it's not applied
# # Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
# #      ifelse(SPECIES==18003 & is.na(LIVEWT.reap) & Prop.Dus.Mon.Good.spl>0,
# #           (Tot.shk.livewt*Prop.Dus.Mon.Good.spl),
# #     ifelse(SPECIES==17001 & is.na(LIVEWT.reap) & Prop.Gum.Mon.Good.spl>0,
# #          (Tot.shk.livewt*Prop.Gum.Mon.Good.spl),
# #     ifelse(SPECIES==17003 & is.na(LIVEWT.reap) & Prop.Whi.Mon.Good.spl>0,
# #         (Tot.shk.livewt*Prop.Whi.Mon.Good.spl),
# #     LIVEWT.reap))))
# 
#     #third, if previous are not available, use "YrZn.Good.spl"             #NEW 
# Bad.Reporters$LIVEWT.reap=with(Bad.Reporters,
#      ifelse(SPECIES==18003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Dus.YrZn.Good.spl),
#     ifelse(SPECIES==17001 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Gum.YrZn.Good.spl),
#     ifelse(SPECIES==17003 & is.na(LIVEWT.reap),(Tot.shk.livewt*Prop.Whi.YrZn.Good.spl),
#     LIVEWT.reap))))
# 
# 
# #remove artificially created 0 catches
# #Bad.Reporters=subset(Bad.Reporters,!is.na(LIVEWT.reap))
# 
# #create file for flagging bad reporters
# #Flag.bad.rep3=Bad.Reporters[,match(c("Same.return","Spec.old","Reporter"),names(Bad.Reporters))]
# 
# #create new vars
# Bad.Reporters$Reporter.old=Bad.Reporters$Reporter
# Bad.Reporters$Reporter="good"
# 
# 
# #merge reapportioned catch
# Data.monthly=rbind(Data.monthly,Bad.Reporters)



# #4.9 Create data for tracking mean weight
# ALL.THIS=c(names(Data.daily.agg.Numbers)[-18],"Km.Gillnet.Days.c","Km.Gillnet.Hours.c","LIVEWT.c","NETLEN.c")
# Mean.w.whiskery=Data.monthly.GN.whiskery[,match(ALL.THIS,names(Data.monthly.GN.whiskery))]
# Mean.w.gummy=Data.monthly.GN.gummy[,match(ALL.THIS,names(Data.monthly.GN.gummy))]
# Mean.w.dusky=Data.monthly.GN.dusky[,match(ALL.THIS,names(Data.monthly.GN.dusky))]
# Mean.w.sandbar=Data.monthly.GN.sandbar[,match(ALL.THIS,names(Data.monthly.GN.sandbar))]
# 
# This.yrs.weight=sort(unique(Data.daily.agg.Numbers$FINYEAR))
# 
# Mean.w.whiskery=subset(Mean.w.whiskery,FINYEAR%in%This.yrs.weight)
# Mean.w.gummy=subset(Mean.w.gummy,FINYEAR%in%This.yrs.weight)
# Mean.w.dusky=subset(Mean.w.dusky,FINYEAR%in%This.yrs.weight)
# Mean.w.sandbar=subset(Mean.w.sandbar,FINYEAR%in%This.yrs.weight)
# 
#   
# THIS.N=subset(Data.daily.agg.Numbers,METHOD=="GN",select=c(Same.return,nfish,SPECIES))
# Mean.w.whiskery=merge(Mean.w.whiskery,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.gummy=merge(Mean.w.gummy,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.dusky=merge(Mean.w.dusky,THIS.N,by=c("Same.return","SPECIES"),all.x=T)
# Mean.w.sandbar=merge(Mean.w.sandbar,THIS.N,by=c("Same.return","SPECIES"),all.x=T)


#   #Mean weigth data
# write.csv(Mean.w.whiskery,file ="Mean.w.whiskery.GN.csv")
# write.csv(Mean.w.gummy,file ="Mean.w.gummy.GN.csv")
# write.csv(Mean.w.dusky,file ="Mean.w.dusky.GN.csv")
# write.csv(Mean.w.sandbar,file ="Mean.w.sandbar.GN.csv")


# #Rescale catches again
# #note: recalculate the reapportioned catches of dusky, whiskery, gummy and 'other'
# #     considering the real catches of other shark species
# a=unique(Bad.Reporters$Same.return)
# b=subset(Data.monthly,Same.return%in%a & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26))
# 
# if(Reapportion.daily=="YES") 
# {
#   a.daily=unique(Bad.Reporters.daily$Same.return)
#   b.daily=subset(Data.daily,Same.return%in%a.daily & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26))  
# }
# 
# xx=function(dat)
# {
#   x=NA
#   if(round(sum(dat$LIVEWT.reap,na.rm=T))>round(unique(dat$Tot.shk.livewt))) x=unique(dat$Same.return)
#   return(x)
# }
# 
# VECT=rep(NA,length(a))
# for ( i in 1:length(a))VECT[i]=xx(subset(b,Same.return==a[i]))
# vect=VECT[!is.na(VECT)]   #all returns where reapportioned weight > total shark weight
# 
# if(Reapportion.daily=="YES")
# {
#   VECT.daily=rep(NA,length(a.daily))
#   for ( i in 1:length(a.daily))VECT.daily[i]=xx(subset(b.daily,Same.return==a.daily[i]))
#   vect.daily=VECT.daily[!is.na(VECT.daily)]   #all returns where reapportioned weight > total shark weight  
# }
# 
#   #monthly
# STOREss=vector("list",length(vect))
# for (i in 1:length(vect))
# {
#   s=subset(b,Same.return==vect[i])
#   Non.Fixed.shk.wgt=unique(s$Tot.shk.livewt)-sum(s$LIVEWT[which(!s$SPECIES%in%Fix.species)])
#   s$Non.Fixed.shk.wgt=Non.Fixed.shk.wgt  
#   Tot.reap=sum(s$LIVEWT.reap[which(s$SPECIES%in%Fix.species)])  
#   s$LIVEWT.reap=with(s,ifelse(SPECIES%in%Fix.species,
#                               Non.Fixed.shk.wgt*(LIVEWT.reap/Tot.reap),LIVEWT.reap))   
#   STOREss[[i]]=s
# }
# bb=do.call(rbind,STOREss)
# bb=bb[,-match("Non.Fixed.shk.wgt",names(bb))]
# Data.monthly=subset(Data.monthly,!(Same.return%in%vect & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26)))
# Data.monthly=rbind(Data.monthly,bb)
# 
# 
# #daily
# if(Reapportion.daily=="YES")
# {
#   STOREss.daily=vector("list",length(vect.daily))
#   for (i in 1:length(vect.daily))
#   {
#     s=subset(b.daily,Same.return==vect.daily[i])
#     Non.Fixed.shk.wgt=unique(s$Tot.shk.livewt)-sum(s$LIVEWT[which(!s$SPECIES%in%Fix.species)])
#     s$Non.Fixed.shk.wgt=Non.Fixed.shk.wgt  
#     Tot.reap=sum(s$LIVEWT.reap[which(s$SPECIES%in%Fix.species)])  
#     s$LIVEWT.reap=with(s,ifelse(SPECIES%in%Fix.species,
#                                 Non.Fixed.shk.wgt*(LIVEWT.reap/Tot.reap),LIVEWT.reap))   
#     STOREss.daily[[i]]=s
#   }
#   bb.daily=do.call(rbind,STOREss.daily)
#   bb.daily=bb.daily[,-match("Non.Fixed.shk.wgt",names(bb.daily))]
#   Data.daily=subset(Data.daily,!(Same.return%in%vect.daily & SPECIES%in%Elasmo.species & METHOD=="GN" & LAT<=(-26)))
#   Data.daily=rbind(Data.daily,bb.daily) 
# }

