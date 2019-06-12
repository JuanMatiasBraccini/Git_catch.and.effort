new.level.vessel=sort(unique(DATA.list.LIVEWT.c[[1]]$new.level.vessel))[1]  #1= worst vessel
#whiksery 
Grid.data[[1]]=list(bin=expand.grid(FINYEAR=FINYEAR,BLOCKX=BLOCKX.w,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[1]]],Freo=FREO[1],log.Effort=log.Effort),
                    pos=expand.grid(FINYEAR=FINYEAR,BLOCKX=BLOCKX.w,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[1]]],log.Effort=log.Effort))

#gummy
Grid.data[[2]]=list(bin=expand.grid(FINYEAR=FINYEAR,BLOCKX=BLOCKX.g,new.level.vessel=new.level.vessel,log.Effort=log.Effort),
                    pos=expand.grid(FINYEAR=FINYEAR,BLOCKX=BLOCKX.g,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[2]]],log.Effort=log.Effort))

#dusky
Grid.data[[3]]=list(bin=NULL,
                    pos=expand.grid(FINYEAR=FINYEAR,BLOCKX=BLOCKX.d,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[3]]],Freo=FREO[3],log.Effort=log.Effort))

#sandbar
Grid.data[[4]]=list(bin=expand.grid(FINYEAR=FINYEAR.san,BLOCKX=BLOCKX.s,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[4]]],Freo=FREO[4],log.Effort=log.Effort),
                    pos=expand.grid(FINYEAR=FINYEAR.san,BLOCKX=BLOCKX.s,new.level.vessel=new.level.vessel,MONTH=MONTH[Top.MN[[4]]],log.Effort=log.Effort))




i=4


Pred.BaseCase=Stand.BaseCase[[i]]
Pred.GoodRep=Stand.GoodRep[[i]]

#predict catch
Pred.BaseCase=predict.ktch(Pred.BaseCase,Grid.data[[i]],SPECIES.vec[[i]])
Pred.GoodRep=predict.ktch(Pred.GoodRep,Grid.data[[i]],SPECIES.vec[[i]])
Pred.Fantasy=Pred.BaseCase

#impute missing year-blocks      
Pred.BaseCase$Pos.dat=impute.ktch(Pred.BaseCase$Pos.dat)
Pred.GoodRep$Pos.dat=impute.ktch(Pred.GoodRep$Pos.dat)


LISTA=Pred.BaseCase
#BiasCorr=BiasCor[[i]]
BiasCorr=0                    #NEW
  Area.w=AREA.W[[i]]
  IMPUTE="YES"
SPEC=SPECIES.vec[[i]]

Bi=LISTA$Bi.dat
Pos=LISTA$Pos.dat

#area weight as a proportion
Area.w$WEIGHT=Area.w$Fish.Area/sum(Area.w$Fish.Area)

if(!SPEC=="Dusky shark")Dat=merge(Pos,Bi[,match(c("FINYEAR","BLOCKX","Pred.Prob"),names(Bi))],by=c("FINYEAR","BLOCKX"),all.x=T)
if(SPEC=="Dusky shark")Dat=Pos

#correct log bias  
if(IMPUTE=="YES")Dat$KTCH=exp(Dat$Pred.Catch.imp+BiasCorr) 
if(IMPUTE=="NO")Dat$KTCH=exp(Dat$Pred.Catch+BiasCorr) 

#Multiply Positive catch and Prob of catch
if(!SPEC=="Dusky shark")Dat$I_y.b=Dat$KTCH*Dat$Pred.Prob
if(SPEC=="Dusky shark")Dat$I_y.b=Dat$KTCH

#Add block weight
Dat=merge(Dat,Area.w[,match(c("BLOCKX","WEIGHT"),names(Area.w))],by="BLOCKX",all.x=T)

#Multiply index by weighted (fishable area)
Dat$I_y.b.W=Dat$I_y.b*Dat$WEIGHT

#Aggregate index by year
I_y=aggregate(I_y.b.W~FINYEAR,Dat,sum,na.rm=T)

#express in 1 km gn day
I_y$I_y.b.W=I_y$I_y.b.W/Eff.unit




#Nominal
a=DATA.list.LIVEWT.c[[i]]
a$CPUE=a$Catch.Target/a$Km.Gillnet.Days.c
CPUE=aggregate(CPUE~FINYEAR,a,mean,na.rm=T)


#Folly
Folly.catch=aggregate(Catch.Target~FINYEAR,a,sum,na.rm=T)
Folly.eff=aggregate(Km.Gillnet.Days.c~FINYEAR,a,sum,na.rm=T)
Folly.cpue=Folly.catch$Catch.Target/Folly.eff$Km.Gillnet.Days.c


setwd("C:/Matias/Analyses/Catch and effort/Outputs/CompareCPUES")
tiff(file=paste(SPEC,".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(2,1),mai=c(1,1.25,.1,.1),oma=c(1,1,.1,.1),las=1)

MAXI=max(c(I_y[,2],CPUE$CPUE,Folly.cpue))
plot(I_y[,2],type="l",lwd=2,ylim=c(0,MAXI),ylab="Kg/km gn day",xaxt='n',xlab="",cex.lab=1.5)
lines(CPUE$CPUE,col=2,lwd=2)
lines(Folly.cpue,col=3,lwd=2)
legend("topleft",c("Standardised","Nominal","Folly"),col=1:3,lty=1,bty='n',lwd=2,cex=1.5)
axis(1,1:length(CPUE$FINYEAR),CPUE$FINYEAR)
mtext(paste(SPEC,", vessel cat.=",as.character(new.level.vessel),", top month",sep=""),3,line=-1,outer=F,cex=1.1)



#relative
MAXI=max(c(I_y[,2]/I_y[1,2],CPUE$CPUE/CPUE$CPUE[1],Folly.cpue/Folly.cpue[1]))
plot(I_y[,2]/I_y[1,2],type="l",lwd=2,ylim=c(0,MAXI),ylab="Relative cpue",
     cex.lab=1.5,xaxt='n',xlab="")
lines(CPUE$CPUE/CPUE$CPUE[1],col=2,lwd=2)
lines(Folly.cpue/Folly.cpue[1],col=3,lwd=2)
axis(1,1:length(CPUE$FINYEAR),CPUE$FINYEAR)
mtext("Financial year",1,line=-1,outer=T,cex=1.5)
dev.off()


