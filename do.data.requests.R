
hndl=handl_OneDrive("Analyses/Catch and effort/Data_Resquests")


#G 4.1 AUDITS -------------------------------------------------------------------------
if(do.audit=="YES")
{
  Audit=subset(Data.daily,FINYEAR==Current.yr & METHOD%in%c("GN","LL") & Estuary=="NO" &
                 LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  write.csv(Audit,paste(hndl,"/",Current.yr,".Audit.csv",sep=""),row.names=F)  
}

#G 4.2. Russel Hudson, FishWell consulting -------------------------------------------------------------------------
if(do.Russels=="YES")
{
  Russel=subset(Data.monthly,SPECIES%in%c(18032,19000),select=c(LIVEWT,FINYEAR,SPECIES,SNAME))
  Russel$Fishery="TDGDLF"
  Russel2=subset(Data.monthly.north,SPECIES%in%c(18032,19000),select=c(LIVEWT,FINYEAR,SPECIES,SNAME))
  Russel2$Fishery="NSH"
  Russel=rbind(Russel,Russel2)
  Russel.catch=aggregate(LIVEWT~FINYEAR+SPECIES+SNAME+Fishery,Russel,sum)  
  Russel.catch=Russel.catch[order(Russel.catch$SPECIES,Russel.catch$FINYEAR,Russel.catch$Fishery),]
  write.csv(Russel.catch,paste(hndl,"/Russel.catch.csv",sep=""),row.names=F)  
}

#G 4.3. Carly, TEPS -------------------------------------------------------------------------
#note: check for Tiger sharks in Comments and inconsistencies in the numbers reported in the
#       columns A, D, and in the comments
Carly.yr=2015   #calendar year
if(do.Carlys=="YES")
{
  These.TEPS=TEPS.current
  TEPS.Carly=subset(These.TEPS,fishery%in%c("SGL","WCGL") & year==Carly.yr)
  TEPS.Carly$LAT=-as.numeric(substr(TEPS.Carly$block10,1,2))
  TEPS.Carly$LONG=100+as.numeric(substr(TEPS.Carly$block10,4,5))
  TEPS.Carly$Bioregion=as.character(with(TEPS.Carly,ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
                                                           ifelse(LONG<115.5 & LAT<=(-27),"WC",
                                                                  ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
                                                                         ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))))
  TEPS.Carly$Bioregion=with(TEPS.Carly,
                            ifelse(Bioregion=="SC"& LAT>(-34) & LONG <115.91 ,"WC",Bioregion))
  
  
  Scan.Carly.teps=TEPS.Carly[,match(c("DailySheetNumber","vessel","month",
                                      "Status","Number","DataEntryName","Comments"),names(TEPS.Carly))]
  write.csv(Scan.Carly.teps,"Scan.Carly.teps.csv",row.names=F)
  Scan.Carly.teps$ID=1:nrow(Scan.Carly.teps)
  #check.comments=TEPS.Carly[c(60:71),]  #these bronzies recorded in methods but not in numbers for 2012
  check.comments$Unic=with(check.comments,paste(DailySheetNumber,SessionNumber))
  check.comments=check.comments[!duplicated(check.comments$Unic),]
  
  if(Carly.yr==2012)
  {
    #Change 'Whale' record entry to 'DUSKY WHALER OVERSIZE' as per comments.
    TEPS.Carly$SpeciesCode=with(TEPS.Carly,ifelse(DataEntryName=="WHALES",18003,SpeciesCode))
    TEPS.Carly$DataEntryName=with(TEPS.Carly,
                                  ifelse(DataEntryName=="WHALES","DUSKY WHALER OVERSIZE",DataEntryName)) 
  }
  
  fun.Tab2.SoFaR=function(Dat)
  {
    Dat$Status=as.character(Dat$Status)
    Dat$Status=with(Dat,ifelse(Status=="a","A",
                               ifelse(Status=="d","D",Status)))
    Dat$Ali.Ded.yr=with(Dat,paste(finyear,Status))
    
    TABLA=aggregate(Number~SpeciesCode+Ali.Ded.yr+Bioregion,data=Dat,sum,na.rm=T)
    Species.id=Dat[,match(c("SpeciesCode","DataEntryName"),names(Dat))]
    Species.id=Species.id[!(duplicated(Species.id$SpeciesCode)),]
    
    wide <- reshape(TABLA,v.names="Number",timevar="Ali.Ded.yr",idvar=c("SpeciesCode","Bioregion"),direction="wide")
    wide=wide[order(wide$SpeciesCode),]
    wide=Species.id%>%left_join(wide,by="SpeciesCode")
    #wide=merge(Species.id,wide,by="SpeciesCode",all.x=T)
    wide=wide[-match("SpeciesCode",names(wide))]
    names(wide)=c("Species","Bioregion","alive","dead")
    wide=wide[order(wide$Bioregion,wide$Species),]
    wide=wide[,match(c("Bioregion","Species","alive","dead"),names(wide))]
    return(wide)
  }
  TEP.Carly=fun.Tab2.SoFaR(TEPS.Carly)
  
  
  if(Carly.yr==2015)
  { 
    #manually add records in comments for 1 tiger DailySheetN TDGLF6000137
    Add.Carly=TEP.Carly[3,]
    Add.Carly$Species="TIGER SHARK OVERSIZE"
    Add.Carly$alive=1
    Add.Carly$dead=NA  
    TEP.Carly=rbind(TEP.Carly,Add.Carly) 
  }
  
  if(Carly.yr==2013)
  { 
    #manually add records in comments for 1 tiger DailySheetN TDGLF6000137
    Add.Carly=TEP.Carly[3,]
    Add.Carly$Species="TIGER SHARK OVERSIZE"
    Add.Carly$alive=1
    Add.Carly$dead=NA  
    TEP.Carly=rbind(TEP.Carly,Add.Carly)   
  }
  
  TEP.Carly=TEP.Carly[order(TEP.Carly$Bioregion,TEP.Carly$Species),]
  
  write.csv(TEP.Carly,paste(hndl,"/Carly/Carly.TEPS.",Carly.yr,".csv",sep=""),row.names=F)
  
}
if(exists("TEPS.2011_12"))
{
  #note: This should be added to catch information for assessment. 
  #     For finyears post 2011-12, must manually scan comments for matching what was reported in columns and in comments
  
  # Also note that for most whalers, (e.g. Tigers), comments don't say if oversized or not. So 
  #assumed they are
  Lista.protected.elasmos=list(
    Tiger=c("tiger","Tiger","TIGER","tigersharks","TG","T/S","Tidgers"),
    Dusky=c("dusky","Bronze","bronze","b/w","Dusky","whalers","BW","Bronzy","b/whaler"),
    Copper=c("copper","Copper"),
    Sandbar=c("Thiskskin","Sandbar","thickskins","TK"),
    Greynurse=c("Nurse","NURSE","nurse","g/n","GN","greynurse","G.N","G/N"),
    Blacktip=c("Blacktip"),
    Whites=c("W/S","WHITE","white","White","Pointer"),
    Mako=c("Mako"),
    Sawfish=c("Sawfish"),
    Manta=c("Manta"))
  
  
  Lista.id.pro.el=Check.These.Sks=vector('list',length=length(Lista.protected.elasmos))
  names(Lista.id.pro.el)=names(Check.These.Sks)=names(Lista.protected.elasmos)
  
  fn.find=function(what,inwhat)
  {
    Find=regexpr(what, inwhat$Comments) > 0
    return(which(Find==T))
  }
  
  for (i in 1:length(Lista.protected.elasmos))
  {
    these=Lista.protected.elasmos[[i]]
    store=NULL
    if(length(these)>1)for(j in 1:length(these))store=c(store,fn.find(these[j],TEPS.2011_12))    
    if(length(these)==1)store=fn.find(these,TEPS.2011_12)
    Lista.id.pro.el[[i]]=store
  }
  
  for (i in 1:length(Lista.id.pro.el))
  {
    id=Lista.id.pro.el[[i]]
    a=TEPS.2011_12[id,match(c("DailySheetNumber","SessionNumber",
                              "DataEntryName","Status","Number","Comments"),names(TEPS.2011_12))]
    a=a[order(a$DailySheetNumber,a$SessionNumber),]
    Check.These.Sks[[i]]=a
  }
  
  #The only manual bit is to go thru each element of Check.These.Sks and compare comments with data
  i=1;print(Check.These.Sks[[i]])
  
  #Table
  All.TEP.table=fun.Tab2.SoFaR(TEPS.2011_12)
  TEPS[is.na(TEPS)] = ""
}

#G 4.4 Alex Hexp. Effort of ASL model -------------------------------------------------------------------------
if(do.Alexs=="YES")
{
  Data.daily.Alex=merge(Data.daily.Alex,BlOCK_10,by.x="block10",by.y="BlockNo",all.x=T)
  Effort.alex.daily=Effort.alex.daily[!duplicated(Effort.alex.daily$Same.return.SNo),]
  Effort.alex.daily$Estuary=with(Effort.alex.daily,ifelse(blockx%in%Estuaries,"YES","NO"))
  
  #Effort.alex.daily=subset(Effort.alex.daily,Estuaries=="NO")
  
  KEEp.nms=c(match("Same.return.SNo",names(Data.daily.Alex)),which(!names(Data.daily.Alex)%in%names(Effort.alex.daily)))
  Data.Alex=merge(Data.daily.Alex[,KEEp.nms],Effort.alex.daily,by="Same.return.SNo",all.y=T)
  Data.Alex=Data.Alex[!(duplicated(Data.Alex$Same.return.SNo)),]
  Effort.alex=subset(Effort.alex,YEAR.c>=2005)
  Data.Alex=subset(Data.Alex,year>=2005 & LatDeg>=26)
  Data.Alex=subset(Data.Alex,netlen.c>100 & method=="GN")
  
  Data.Alex.nlines.check=subset(Data.Alex,select=c(Same.return.SNo,vessel,
                                                   netlen,netlen.c,hours,hours.c,nlines,shots,shots.c,fishery))
  Data.Alex=subset(Data.Alex,select=c(Same.return.SNo,finyear,year,month,day,date,
                                      blockx,block10,netlen,netlen.c,hours.c,nlines,nlines.c,shots.c,Long,
                                      LongDeg,LongMin,Lat,
                                      LatDeg,LatMin,Latitude_Centroid,Longitude_Centroid,depthMax,fishery,vessel))
  
  change.names=match(c("Same.return.SNo","blockx","date","netlen","netlen.c","hours.c",
                       "nlines","nlines.c","shots.c"),names(Data.Alex))
  names(Data.Alex)[change.names]=c("SessionID","block","Date","netlen.original","netlen","hours",
                                   "nlines.original","nlines","shots")
  This.Alex=c("vessel","SessionID","fishery","finyear","month","year","Date","netlen.original","netlen",
              "hours","nlines.original","nlines","shots","block10","block","Lat","LatDeg","LatMin",
              "Long","LongDeg","LongMin","Latitude_Centroid","Longitude_Centroid","depthMax")
  Data.Alex=Data.Alex[,match(This.Alex,names(Data.Alex))]
  Data.Alex$Date=as.character(Data.Alex$Date)
  
  fwrite(Data.Alex,paste(hndl,"/Alex_ASL/Data.daily.Alex.csv",sep=""),row.names=F)  
  fwrite(Data.Alex.nlines.check,paste(hndl,"/Alex_ASL/Data.Alex.nlines.check.csv",sep=""),row.names=F)
  fwrite(Effort.alex,paste(hndl,"/Alex_ASL/Data.CAESS.Alex.csv",sep=""),row.names=F)
  
  rm(Data.daily.Alex,Data.Alex)
  
}

#G 4.5 SAFS. Number of vessels catching a target species -------------------------------------------------------------------------
if(do.Garys=="YES")
{
  fn.gary=function(SPEC,YR)
  {
    dat=subset(Data.monthly,SPECIES==SPEC & FINYEAR==YR)
    dat$zoneVes=with(dat,paste(zone,VESSEL))
    unik.ves=unique(dat$VESSEL)
    dat=dat[!duplicated(dat$zoneVes),]
    unik.ves.zn=table(dat$zone)
    
    return(list(Numb.ves.all=length(unik.ves),Numb.ves.zn=unik.ves.zn))
  }
  
  Pinkies=fn.gary(353001,"2012-13")
  Pinkies=unlist(Pinkies)
  names.Pinkies=names(Pinkies)
  Pinkies=matrix(Pinkies,ncol=4)
  colnames(Pinkies)=names.Pinkies
  
  write.csv(Pinkies,paste(hndl,"/Pinkies.for.gary.csv",sep=""),row.names=F)
  
}

#G 4.6 Whiskery shark catch by month (AMM 2014) -------------------------------------------------------------------------
if(do.whiskerys=="YES")
{
  Yr.Whiskery=2005:2013
  Whisk.ktch.West=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="West"),sum)
  Whisk.ktch.Zn1=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="Zone1"),sum)
  Whisk.ktch.Zn2=aggregate(LIVEWT.c~YEAR.c+MONTH,subset(Data.monthly,SPECIES==17003 & YEAR.c%in%Yr.Whiskery &zone=="Zone2"),sum)
  
  closure.ktch=data.frame(YEAR.c=2008:2012,MONTH=9,LIVEWT.c=0)
  Whisk.ktch.Zn1=rbind(Whisk.ktch.Zn1,closure.ktch)
  Whisk.ktch.West=rbind(Whisk.ktch.West,closure.ktch)
  
  Whisk.ktch.West=Whisk.ktch.West[order(Whisk.ktch.West$YEAR.c,Whisk.ktch.West$MONTH),]
  Whisk.ktch.Zn1=Whisk.ktch.Zn1[order(Whisk.ktch.Zn1$YEAR.c,Whisk.ktch.Zn1$MONTH),]
  Whisk.ktch.Zn2=Whisk.ktch.Zn2[order(Whisk.ktch.Zn2$YEAR.c,Whisk.ktch.Zn2$MONTH),]
  
  tiff(file=paste(hndl,"/Whiskery.monthly.catch.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  
  par(las=1)
  MAX.y=max(Whisk.ktch.West$LIVEWT.c/1000,Whisk.ktch.Zn1$LIVEWT.c/1000,Whisk.ktch.Zn2$LIVEWT.c/1000)
  plot(Whisk.ktch.West$LIVEWT.c/1000,ylim=c(0,MAX.y),xaxt='n',col=2,lwd=2,
       type='l',xlab="",ylab="",cex.lab=2,cex.axis=1.5)
  Sep=seq(9,nrow(Whisk.ktch.West),12)
  abline(v=Sep,lwd=3,col="grey80")
  lines(Whisk.ktch.West$LIVEWT.c/1000,col=2,lwd=2)
  lines(Whisk.ktch.Zn1$LIVEWT.c/1000,col=3,lwd=2)
  lines(Whisk.ktch.Zn2$LIVEWT.c/1000,col=4,lwd=2)
  legend("topright",c("West coast","Zone 1","Zone 2"),
         lty=1,col=2:4,lwd=2,cex=2,bg="white",box.col="white")
  axis(1,at=1:nrow(Whisk.ktch.West),F,tck=-0.015)
  YERS=seq(1,nrow(Whisk.ktch.West),12)
  axis(1,at=YERS,Yr.Whiskery,tck=-0.03,cex.axis=1.5)
  box()
  mtext("Year",1,3,cex=2)
  mtext("Whiskery catch (tonnes)",2,2.5,cex=2,las=3)
  dev.off()
  
  names(Whisk.ktch.West)[3]="West"
  names(Whisk.ktch.Zn1)[3]="Zone1"
  names(Whisk.ktch.Zn2)[3]="Zone2"
  Whisk.ktch=merge(Whisk.ktch.West,Whisk.ktch.Zn1,by=c("YEAR.c","MONTH"))
  Whisk.ktch=merge(Whisk.ktch,Whisk.ktch.Zn2,by=c("YEAR.c","MONTH"))
  
  write.csv(Whisk.ktch,paste(hndl,"/Whisk.ktch.csv",sep=""),row.names=F)  
}

#G 4.7 Jeff Norris, all catch and effort for South Coast Bioregion (old request) -------------------------------------------------------------------------
if(do.Jeffs=="YES")
{
  Yrs.of.scalies.data=unique(subset(Data.monthly,SPECIES%in%Scalefish.species,select=FINYEAR))
  Jeff.fn=function(DATA,Eff.d,Eff.m)
  {
    DATA=subset(DATA,select=c(FINYEAR,YEAR,MONTH,VESSEL,METHOD,BLOCKX,SPECIES,SNAME,LAT,LONG,
                              Bioregion,zone,LIVEWT.c))
    SC.blocks=unique(DATA$BLOCKX)
    Eff.m=subset(Eff.m,BLOCKX%in%SC.blocks)  
    Eff.d=subset(Eff.d,blockx%in%SC.blocks)
    
    #Monthly
    Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~VESSEL+FINYEAR+MONTH+BLOCKX,data=Eff.m,max,na.rm=T)
    Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~FINYEAR+MONTH+VESSEL+BLOCKX,data=Eff.monthly.hour.c,sum,na.rm=T) 
    
    #Daily
    Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~ID+vessel+finyear+month+blockx,data=Eff.d,max,na.rm=T)
    Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~finyear+month+vessel+blockx,data=Eff.daily.hour.c,sum,na.rm=T)
    
    names(Eff.daily.hour.c)[1:4]=c("FINYEAR","MONTH","VESSEL","BLOCKX")
    
    return(list(Catch=DATA,Effort=rbind(Eff.monthly.hour.c,Eff.daily.hour.c)))
  }
  Jeff.s=Jeff.fn(subset(Data.monthly,Bioregion=="SC" & FINYEAR%in%Yrs.of.scalies.data$FINYEAR),
                 Effort.daily,subset(Effort.monthly,!FINYEAR%in%Daily.l.years))
  
  write.csv(Jeff.s$Catch,paste(hndl,"/Jeff.Catch.south_coast.csv",sep=""),row.names=F)
  write.csv(Jeff.s$Effort,paste(hndl,"/Jeff.Effort.south_coast.csv",sep=""),row.names=F)
  
}

#G 4.8 Jodie O'Malley -------------------------------------------------------------------------
#G 4.8.1 ASL exclusion areas catch and effort
hndls=paste(hndl,"/Jodie_OMalley_2016/data/",sep="")
if(do.Jodies.ASL=="YES")
{  
  Data.daily$Fishing_yr=with(Data.daily,ifelse(MONTH>=6,YEAR.c,YEAR.c-1))
  Data.daily$FISHERY=with(Data.daily,ifelse(zone=="West","WCDGDLF",
                                            ifelse(zone%in%c("Zone1","Zone2"),"JASDGDLF",NA)))
  Jodie.Effort$FISHERY=with(Jodie.Effort,ifelse(zone=="West","WCDGDLF",
                                                ifelse(zone%in%c("Zone1","Zone2"),"JASDGDLF",NA)))
  
  #add zone 3
  Data.daily$zone3=with(Data.daily,ifelse(Data.daily$LONG>=116.5 & Data.daily$LONG<= 116.923 & Data.daily$LAT < (-33),"Zone3",zone))
  Zn3.blks=unique(subset(Data.daily,zone3=="Zone3",select=block10))$block10
  Jodie.Effort$zone3=with(Jodie.Effort,ifelse(block10%in%Zn3.blks,"Zone3",zone))
  
  
  #some functions
  fn.extract.Jodie.catch=function(blk,FSHRY,YRS,SP,Byzone)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP))
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    
    if(Byzone=="NO")
    {
      Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+VESSEL+block10+Spec.tab1,dat.blk,sum)
      Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,dat.fishery,sum)
    }
    
    if(Byzone=="YES")
    {
      Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+VESSEL+block10+Spec.tab1+zone3,dat.blk,sum)
      Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,dat.fishery,sum) 
    }
    
    return(list(BLKs=Aggregated.blk,Fishery=Aggregated.fishery))
  }
  
  fn.extract.Jodie.effort=function(blk,FSHRY,YRS,Byzone)
  {
    Dat=subset(Jodie.Effort)
    
    dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    
    if(Byzone=="NO")
    {
      Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+block10,dat.blk,sum)
      Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    }
    
    if(Byzone=="YES")
    {
      Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+vessel+block10+zone3,dat.blk,sum)
      Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+zone3,dat.fishery,sum)     
    }
    
    return(list(BLKs=Aggregated.blk_km.gn.hr,Fishery=Aggregated.fishery_km.gn.hr))
  }
  
  fn.agg.Jodie.efrt=function(dat,Byzone)
  {
    BLK=dat$BLKs
    All=dat$Fishery
    
    if(Byzone=="NO")Blk.ag=aggregate(Km.Gillnet.Hours.c~Fishing_yr,BLK,sum)
    if(Byzone=="YES")Blk.ag=aggregate(Km.Gillnet.Hours.c~Fishing_yr+zone3,BLK,sum)
    
    return(list(BLK=Blk.ag,ALL=All))
  }
  
  fn.tab.eff=function(dat,Byzone)
  {
    if(Byzone=="NO")
    {
      D=merge(dat$ALL,dat$BLK,by="Fishing_yr")
      names(D)[2:3]=c("All","Blocks") 
      D$Prop=D$Blocks/D$All
    }
    
    if(Byzone=="YES")
    {
      D=merge(dat$ALL,dat$BLK,by=c("Fishing_yr","zone3"),all.x=T)
      names(D)[3:4]=c("All","Blocks")    
      D$Prop=D$Blocks/D$All
      
      D=reshape(D, v.names = c("All","Blocks","Prop"), idvar = "Fishing_yr",
                timevar = "zone3", direction = "wide")
      D=subset(D,select=c(Fishing_yr,All.Zone1,All.Zone2,All.Zone3,Blocks.Zone2,Prop.Zone2))
      
    }
    
    
    return(D)
  }
  
  combo.tbl_eff=function(A,B,C,nms)
  {
    D=merge(A,B,by="Fishing_yr")
    D=merge(D,C,by="Fishing_yr")
    D=D[,match(nms,colnames(D))]
    return(D)
  }
  
  fn.agg.Jodie.ktch=function(dat,SORT,Byzone)
  {
    BLK=dat$BLKs
    All=dat$Fishery
    SORT=SORT[which(SORT%in%unique(All$Spec.tab1))]
    
    
    if(Byzone=="NO")
    {
      BLK$Spec.tab1=as.character(BLK$Spec.tab1)
      All$Spec.tab1=as.character(All$Spec.tab1)
      
      Blk.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,BLK,sum)
      Blk.ag.wide<- reshape(Blk.ag,v.names="LIVEWT.c",timevar="Fishing_yr",idvar="Spec.tab1",direction="wide")
      
      All.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1,All,sum)
      All.ag.wide<- reshape(All.ag,v.names="LIVEWT.c",timevar="Fishing_yr",idvar="Spec.tab1",direction="wide")
      
      Add.this=All.ag.wide$Spec.tab1[which(!All.ag.wide$Spec.tab1%in%Blk.ag.wide$Spec.tab1)]
      if(length(Add.this)>0)
      {
        addd=Blk.ag.wide[1:length(Add.this),] 
        addd[,]=NA
        addd$Spec.tab1=Add.this
        Blk.ag.wide=rbind(Blk.ag.wide,addd)    
      }
      
      Blk.ag.wide=Blk.ag.wide[match(SORT,Blk.ag.wide$Spec.tab1),]
      All.ag.wide=All.ag.wide[match(SORT,All.ag.wide$Spec.tab1),]
      
      fn.rshp=function(x)
      {
        x=reshape(x,direction="long")
        x=x[order(x$Fishing_yr),]
      }
      Blk.ag.wide=fn.rshp(Blk.ag.wide)
      All.ag.wide=fn.rshp(All.ag.wide)    
    }
    
    if(Byzone=="YES")
    {
      BLK$Spec.tab1=as.character(BLK$Spec.tab1)
      All$Spec.tab1=as.character(All$Spec.tab1)
      
      Blk.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,BLK,sum)  #no need for zone here, only zone2 has closures
      Blk.ag.wide<- reshape(Blk.ag,v.names="LIVEWT.c",timevar="zone3",idvar=c("Fishing_yr","Spec.tab1"),direction="wide")
      
      
      All.ag=aggregate(LIVEWT.c~Fishing_yr+Spec.tab1+zone3,All,sum)
      All.ag.wide<- reshape(All.ag,v.names="LIVEWT.c",timevar=c("zone3"),idvar=c("Fishing_yr","Spec.tab1"),direction="wide")
      
      Add.this=merge(All.ag.wide,Blk.ag.wide,by=c("Fishing_yr","Spec.tab1"),all.x=T)
      Add.this=Add.this[which(!with(Add.this,paste(Spec.tab1,Fishing_yr))%in%with(Blk.ag.wide,paste(Spec.tab1,Fishing_yr))),]
      if(nrow(Add.this)>0)
      {
        addd=Blk.ag.wide[1:nrow(Add.this),] 
        addd[,]=NA
        addd$Fishing_yr=Add.this$Fishing_yr
        addd$Spec.tab1=Add.this$Spec.tab1
        Blk.ag.wide=rbind(Blk.ag.wide,addd)    
      }
      Yrs=unique(Blk.ag.wide$Fishing_yr)
      Yrs.SORT=paste(rep(SORT,each=length(Yrs)),rep(Yrs,times=length(SORT)))
      id.blk=match(Yrs.SORT,paste(Blk.ag.wide$Spec.tab1,Blk.ag.wide$Fishing_yr))
      Blk.ag.wide=Blk.ag.wide[id.blk,]
      id.all=match(Yrs.SORT,paste(All.ag.wide$Spec.tab1,All.ag.wide$Fishing_yr))
      All.ag.wide=All.ag.wide[id.all,]
      
      if(length(Yrs.SORT[which(is.na(id.blk))])>0)
      {
        Blk.ag.wide$Fishing_yr[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.blk))], " "), "[[", 2)
        Blk.ag.wide$Spec.tab1[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.blk))], " "), "[[", 1)
      }
      if(length(Yrs.SORT[which(is.na(id.all))])>0)
      {
        All.ag.wide$Fishing_yr[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.all))], " "), "[[", 2)
        All.ag.wide$Spec.tab1[which(is.na(id.blk))]=sapply(strsplit(Yrs.SORT[which(is.na(id.all))], " "), "[[", 1)
      }
      
    }
    
    return(list(BLK=Blk.ag.wide,ALL=All.ag.wide))
  }
  
  combo.tbl=function(A,B,C,nms,SORT) 
  {
    if(nrow(A$BLK)==sum(A$BLK$Spec.tab1==A$ALL$Spec.tab1)) A1=cbind(A$ALL,A$BLK$LIVEWT.c)
    colnames(A1)[3:4]=c("WCDGDLF","Closures_WCDGDLF")
    
    if(nrow(B$BLK)==sum(B$BLK$Spec.tab1==B$ALL$Spec.tab1)) B1=cbind(B$ALL,B$BLK$LIVEWT.c)
    colnames(B1)[3:4]=c("JASDGDLF","Closures_JASDGDLF")
    
    if(nrow(C$BLK)==sum(C$BLK$Spec.tab1==C$ALL$Spec.tab1)) C1=merge(C$ALL,C$BLK,by=c("Fishing_yr","Spec.tab1"),all.x=T) 
    colnames(C1)[3:6]=c("JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3","Closures_JASDGDLF")
    
    A1[is.na(A1)]=0
    B1[is.na(B1)]=0
    C1[is.na(C1)]=0
    
    A1$Prop_WCDGDLF=with(A1,Closures_WCDGDLF/WCDGDLF)
    B1$Prop_JASDGDLF=with(B1,Closures_JASDGDLF/JASDGDLF)
    C1$Prop_JASDGDLF_Zone2=with(C1,Closures_JASDGDLF/JASDGDLF_Zone2)
    
    D=merge(A1,B1,by=c("Spec.tab1","Fishing_yr"),all=T)
    D=merge(D,C1,by=c("Spec.tab1","Fishing_yr","Closures_JASDGDLF"),all=T)
    D[is.na(D)]=0
    D$TDGDLF=with(D,WCDGDLF+JASDGDLF)
    D$Closures_TDGDLF=with(D,Closures_WCDGDLF+Closures_JASDGDLF)
    D$Prop_TDGDLF=with(D,Closures_TDGDLF/TDGDLF)                       
    
    D=D[,match(nms,colnames(D))]
    
    SORT=SORT[which(SORT%in%unique(D$Spec.tab1))]
    Yrs=unique(D$Fishing_yr)
    Yrs.SORT=paste(rep(SORT,each=length(Yrs)),rep(Yrs,times=length(SORT)))
    D=D[match(Yrs.SORT,paste(D$Spec.tab1,D$Fishing_yr)),]
    
    
    return(D)
  }
  
  combo.tbl_agg=function(A)
  {
    d=aggregate(cbind(TDGDLF,WCDGDLF,JASDGDLF,
                      JASDGDLF_Zone1,JASDGDLF_Zone2,JASDGDLF_Zone3,
                      Closures_TDGDLF,Closures_WCDGDLF,Closures_JASDGDLF)~Fishing_yr,A,sum)
    d$Prop_TDGDLF=d$Closures_TDGDLF/d$TDGDLF
    d$Prop_WCDGDLF=d$Closures_WCDGDLF/d$WCDGDLF
    d$Prop_JASDGDLF=d$Closures_JASDGDLF/d$JASDGDLF
    d$Prop_JASDGDLF_Zone2=d$Closures_JASDGDLF/d$JASDGDLF_Zone2
    return(d)  
  }
  
  
  #loop over requested periods
  Yrs.Jodie=list(year_2010_14=2010:2014,year_2007_09=2007:2009)  
  
  for(j in 1:length(Yrs.Jodie))
  {
    
    #Extracth catch
    KTCH_WCDGDLF_elasmo=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"NO")
    KTCH_WCDGDLF_scalies=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"NO")
    KTCH_JASDGDLF_elasmo=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"NO")
    KTCH_JASDGDLF_scalies=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"NO")
    
    KTCH_JASDGDLF_elasmo_zone=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.elasmo,"YES")
    KTCH_JASDGDLF_scalies_zone=fn.extract.Jodie.catch(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],Spec.tab.1.scalies,"YES")
    
    
    #Extract effort
    EFFORT_WCDGDLF=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_WCDGDLF),"WCDGDLF",Yrs.Jodie[[j]],"NO")
    EFFORT_JASDGDLF=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],"NO")
    
    EFFORT_JASDGDLF_zone=fn.extract.Jodie.effort(unlist(ASL_exclusions_block10_JASDGDLF),"JASDGDLF",Yrs.Jodie[[j]],"YES")
    
    
    #Tabulate effort  
    EFRT.Jodie_WCDGDLF=fn.agg.Jodie.efrt(EFFORT_WCDGDLF,"NO")
    EFRT.Jodie_JASDGDLF=fn.agg.Jodie.efrt(EFFORT_JASDGDLF,"NO")
    EFRT.Jodie_JASDGDLF_zone=fn.agg.Jodie.efrt(EFFORT_JASDGDLF_zone,"YES")
    
    Tab.eff.Jodie_WCDGDLF=fn.tab.eff(EFRT.Jodie_WCDGDLF,"NO")
    Tab.eff.Jodie_JASDGDLF=fn.tab.eff(EFRT.Jodie_JASDGDLF,"NO")
    Tab.eff.Jodie_JASDGDLF_zone=fn.tab.eff(EFRT.Jodie_JASDGDLF_zone,"YES")
    
    names(Tab.eff.Jodie_WCDGDLF)[2:4]=c("WCDGDLF","Closures_WCDGDLF","Prop_WCDGDLF")
    names(Tab.eff.Jodie_JASDGDLF)[2:4]=c("JASDGDLF","Closures_JASDGDLF","Prop_JASDGDLF")
    names(Tab.eff.Jodie_JASDGDLF_zone)[2:6]=c("JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3","Closures_JASDGDLF_Zone2","Prop_JASDGDLF_Zone2")
    Tab.eff.Jodie_JASDGDLF=merge(Tab.eff.Jodie_JASDGDLF,subset(Tab.eff.Jodie_JASDGDLF_zone,
                                                               select=c(Fishing_yr,JASDGDLF_Zone1,JASDGDLF_Zone2,JASDGDLF_Zone3,Prop_JASDGDLF_Zone2)),by="Fishing_yr")
    
    Tab.eff.Jodie=merge(Tab.eff.Jodie_WCDGDLF,Tab.eff.Jodie_JASDGDLF,by="Fishing_yr")  
    Tab.eff.Jodie.TDGDLF=Tab.eff.Jodie[,1:2]
    names(Tab.eff.Jodie.TDGDLF)[2]="TDGDLF"
    Tab.eff.Jodie.TDGDLF$TDGDLF=Tab.eff.Jodie_WCDGDLF$WCDGDLF+Tab.eff.Jodie_JASDGDLF$JASDGDLF
    Tab.eff.Jodie.TDGDLF$Closures_TDGDLF=Tab.eff.Jodie_WCDGDLF$Closures_WCDGDLF+Tab.eff.Jodie_JASDGDLF$Closures_JASDGDLF
    Tab.eff.Jodie.TDGDLF$Prop_TDGDLF=with(Tab.eff.Jodie.TDGDLF,Closures_TDGDLF/TDGDLF)
    Jodie.Eff=combo.tbl_eff(A=Tab.eff.Jodie.TDGDLF,B=Tab.eff.Jodie_WCDGDLF,C=Tab.eff.Jodie_JASDGDLF,
                            nms=c("Fishing_yr","TDGDLF","WCDGDLF","JASDGDLF",
                                  "JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3",
                                  "Closures_TDGDLF","Closures_WCDGDLF","Closures_JASDGDLF",
                                  "Prop_TDGDLF","Prop_WCDGDLF","Prop_JASDGDLF","Prop_JASDGDLF_Zone2"))
    
    Jodie.Eff[,11:14]=100*Jodie.Eff[,11:14]   #percentage
    Jodie.Eff[,2:14]=round(Jodie.Eff[,2:14],1) #round
    
    write.csv(Jodie.Eff,paste(hndls,"Jodie.Eff_km.gn.hr_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    
    
    #tabulate catch
    SORT.s=c("Gummy","Dusky_whaler","Whiskery","Sandbar","Hammerheads","Spinner","Wobbegongs","Rays","Common_saw_shark","School",
             "Other_elasmobranchs")
    KTCH.Jodie_WCDGDLF_elasmo=fn.agg.Jodie.ktch(KTCH_WCDGDLF_elasmo,SORT.s,"NO")
    KTCH.Jodie_JASDGDLF_elasmo=fn.agg.Jodie.ktch(KTCH_JASDGDLF_elasmo,SORT.s,"NO")
    KTCH.Jodie_JASDGDLF_elasmo_zone=fn.agg.Jodie.ktch(KTCH_JASDGDLF_elasmo_zone,SORT.s,"YES")
    
    SORT.t=c("Blue_morwong","Blue_groper","West_Australian_dhufish","Pink_snapper","Boarfishes","Samsonfish","Redfishes","Mulloway","Sweetlips",
             "Baldchin_groper","Other_scalefish")
    KTCH.Jodie_WCDGDLF_scalies=fn.agg.Jodie.ktch(KTCH_WCDGDLF_scalies,SORT=SORT.t,"NO")
    KTCH.Jodie_JASDGDLF_scalies=fn.agg.Jodie.ktch(KTCH_JASDGDLF_scalies,SORT.t,"NO")
    KTCH.Jodie_JASDGDLF_scalies_zone=fn.agg.Jodie.ktch(KTCH_JASDGDLF_scalies_zone,SORT.t,"YES")
    
    NMS=c("Fishing_yr","Spec.tab1","TDGDLF","WCDGDLF","JASDGDLF",
          "JASDGDLF_Zone1","JASDGDLF_Zone2","JASDGDLF_Zone3",
          "Closures_TDGDLF","Closures_WCDGDLF","Closures_JASDGDLF",
          "Prop_TDGDLF","Prop_WCDGDLF","Prop_JASDGDLF","Prop_JASDGDLF_Zone2")
    
    
    #elasmos
    Jodie.KTCH.elasmos=combo.tbl(A=KTCH.Jodie_WCDGDLF_elasmo,B=KTCH.Jodie_JASDGDLF_elasmo,
                                 C=KTCH.Jodie_JASDGDLF_elasmo_zone,nms=NMS,SORT=SORT.s)
    Jodie.KTCH.elasmos_agg=combo.tbl_agg(A=Jodie.KTCH.elasmos)
    
    Jodie.KTCH.elasmos[,3:11]=Jodie.KTCH.elasmos[,3:11]/1000        #convert to tonnes
    Jodie.KTCH.elasmos_agg[,2:10]=Jodie.KTCH.elasmos_agg[,2:10]/1000
    
    Jodie.KTCH.elasmos[,12:15]=100*Jodie.KTCH.elasmos[,12:15]       #percentage
    Jodie.KTCH.elasmos_agg[,11:14]=100*Jodie.KTCH.elasmos_agg[,11:14]
    
    Jodie.KTCH.elasmos[,3:15]=round(Jodie.KTCH.elasmos[,3:15],1)  #round
    Jodie.KTCH.elasmos_agg[,2:14]=round(Jodie.KTCH.elasmos_agg[,2:14],1)
    
    write.csv(Jodie.KTCH.elasmos,paste(hndls,"Jodie.KTCH.elasmos_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    write.csv(Jodie.KTCH.elasmos_agg,paste(hndls,"Jodie.KTCH.elasmos_agg_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    
    #teleost
    Jodie.KTCH.scalies=combo.tbl(A=KTCH.Jodie_WCDGDLF_scalies,B=KTCH.Jodie_JASDGDLF_scalies,
                                 C=KTCH.Jodie_JASDGDLF_scalies_zone,nms=NMS,SORT=SORT.t)
    Jodie.KTCH.scalies_agg=combo.tbl_agg(A=Jodie.KTCH.scalies)
    
    Jodie.KTCH.scalies[,3:11]=Jodie.KTCH.scalies[,3:11]/1000        #convert to tonnes
    Jodie.KTCH.scalies_agg[,2:10]=Jodie.KTCH.scalies_agg[,2:10]/1000
    
    Jodie.KTCH.scalies[,12:15]=100*Jodie.KTCH.scalies[,12:15]       #percentage
    Jodie.KTCH.scalies_agg[,11:14]=100*Jodie.KTCH.scalies_agg[,11:14]
    
    Jodie.KTCH.scalies[,3:15]=round(Jodie.KTCH.scalies[,3:15],1)  #round
    Jodie.KTCH.scalies_agg[,2:14]=round(Jodie.KTCH.scalies_agg[,2:14],1)
    
    write.csv(Jodie.KTCH.scalies,paste(hndls,"Jodie.KTCH.scalies_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
    write.csv(Jodie.KTCH.scalies_agg,paste(hndls,"Jodie.KTCH.scalies_agg_",names(Yrs.Jodie)[j],".csv",sep=""),row.names=F)
  }
  
  #Specific Request Block 33240 years 2007-2014
  BLOCK33240_closures=c(333240, 333241, 334240, 334241,334242, 335240, 335241)
  
  fn.extract.Jodie.catch_BLOCK=function(BLK,blk,YRS,SP)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP) & BLOCKX==BLK & Fishing_yr%in%YRS)
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    dat.blk=subset(Dat,block10%in%blk )
    dat.fishery=Dat
    Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr,dat.blk,sum)
    Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr,dat.fishery,sum)
    Agg.dat=merge(Aggregated.fishery,Aggregated.blk,by="Fishing_yr")
    colnames(Agg.dat)[2:3]=c("Total","Closure")
    Agg.dat$Percent=100* Agg.dat$Closure/ Agg.dat$Total
    Agg.dat[2:3]=Agg.dat[2:3]/1000  #convert to tonnes
    Agg.dat[2:4]=round(Agg.dat[2:4],1)
    
    Dat$UNIKS=with(Dat,paste(Fishing_yr,block10,VESSEL))
    dat.Ves=Dat[!duplicated(Dat$UNIKS),]
    
    dat.Ves$VESSEL=as.character(dat.Ves$VESSEL)
    dat.Ves$block10=as.factor(dat.Ves$block10)
    BB=vector('list',length(YRS))
    for(p in 1:length(YRS))
    {
      x=subset(dat.Ves,Fishing_yr==YRS[p])
      x=with(x,table(block10,VESSEL))
      B=data.frame(Fishing_yr=YRS[p])
      BB[[p]]=cbind(B,as.data.frame(t(rowSums(x))))
      
    }
    N.vessels=do.call(rbind,BB)
    
    return(list(Agg.dat=Agg.dat,N.vessels=N.vessels))
  }
  
  fn.extract.Jodie.effort_BLOCK=function(BLK,blk,YRS)
  {
    Dat=subset(Jodie.Effort,block10%in%BLK & Fishing_yr%in%YRS)
    dat.blk=subset(Dat,block10%in%blk)
    dat.fishery=Dat
    Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.blk,sum)
    Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    
    Agg.dat=merge(Aggregated.fishery_km.gn.hr,Aggregated.blk_km.gn.hr,by="Fishing_yr")
    colnames(Agg.dat)[2:3]=c("Total","Closure")
    Agg.dat$Percent=100* Agg.dat$Closure/ Agg.dat$Total
    Agg.dat[2:4]=round(Agg.dat[2:4],1)
    return(Agg.dat)
  }
  
  
  #Extracth catch
  KTCH_33240_elasmo=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.elasmo)$Agg.dat
  KTCH_33240_scalies=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.scalies)$Agg.dat
  Vessels_33240=fn.extract.Jodie.catch_BLOCK(BLK=33240,blk=BLOCK33240_closures,YRS=2007:2014,SP=Spec.tab.1.elasmo)$N.vessels
  
  #Extract effort
  blk10in33240=subset(Data.daily,BLOCKX==33240)
  blk10in33240=unique(blk10in33240$block10)
  blk10in33240=blk10in33240[!is.na(blk10in33240)]
  EFFORT_33240=fn.extract.Jodie.effort_BLOCK(BLK=blk10in33240,blk=BLOCK33240_closures,YRS=2007:2014)
  
  write.csv(EFFORT_33240,paste(hndls,"EFFORT_blk_33240_2007_14.csv",sep=""),row.names=F)  
  write.csv(KTCH_33240_elasmo,paste(hndls,"KTCH_elasmo_blk_33240_2007_14.csv",sep=""),row.names=F)
  write.csv(KTCH_33240_scalies,paste(hndls,"KTCH_scalies_blk_33240_2007_14.csv",sep=""),row.names=F)
  write.csv(Vessels_33240,paste(hndls,"Vessels_blk_33240_2007_14.csv",sep=""),row.names=F)
  
  
  #27 June 2016 Request
  ASL_block_27.6.request=c(33240,32250,33250,32260,33260,32270,33270,32280,33280)
  ASL_block10_27.6.request=c(322255,322260,321260)
  Yrs.Jodie_27.6.request=2007:2014
  
  fn.extract.Jodie.catch.Req.27.6=function(blk,YRS,SP,FSHRY,what.blk,Zn)
  {
    Dat=subset(Data.daily,SPECIES%in%unlist(SP))
    Dat$Spec.tab1=NA
    SPECIES=names(SP)
    for(p in 1:length(SPECIES))Dat$Spec.tab1=with(Dat,ifelse(SPECIES%in%SP[[p]],names(SP[p]),Spec.tab1))
    
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    dat.zn=subset(Dat,zone==Zn & Fishing_yr%in%YRS)
    
    if(what.blk=='block')dat.blk=subset(Dat,BLOCKX%in%blk & Fishing_yr%in%YRS)
    if(what.blk=='block10')dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    
    if(what.blk=='block') Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+BLOCKX,dat.blk,sum)  
    if(what.blk=='block10') Aggregated.blk=aggregate(LIVEWT.c~Fishing_yr+block10,dat.blk,sum)
    
    Aggregated.fishery=aggregate(LIVEWT.c~Fishing_yr,dat.fishery,sum)
    Aggregated.zone=aggregate(LIVEWT.c~Fishing_yr,dat.zn,sum)
    
    return(list(BLKs=Aggregated.blk,Fishery=Aggregated.fishery,Zone=Aggregated.zone))
  }
  fn.extract.Jodie.effort.Req.27.6=function(blk,FSHRY,YRS,what.blk,Zn)
  {
    Dat=Jodie.Effort
    Dat$BLOCKX=with(Dat, as.numeric(paste(substr(block10,1,2),substr(block10,4,5),0,sep="")))
    
    if(what.blk=='block')dat.blk=subset(Dat,BLOCKX%in%blk & Fishing_yr%in%YRS)
    if(what.blk=='block10')dat.blk=subset(Dat,block10%in%blk & Fishing_yr%in%YRS)
    
    dat.fishery=subset(Dat,FISHERY==FSHRY & Fishing_yr%in%YRS)
    dat.zn=subset(Dat,zone==Zn & Fishing_yr%in%YRS)
    
    if(what.blk=='block')Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+BLOCKX,dat.blk,sum)
    if(what.blk=='block10')Aggregated.blk_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr+block10,dat.blk,sum)
    
    Aggregated.fishery_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.fishery,sum)
    Aggregated.zone_km.gn.hr=aggregate(Km.Gillnet.Hours.c~Fishing_yr,dat.zn,sum)
    
    return(list(BLKs=Aggregated.blk_km.gn.hr,Fishery=Aggregated.fishery_km.gn.hr,Zone=Aggregated.zone_km.gn.hr))
  }
  
  #Elasmos catch
  Elasmo.ktch.block_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block_27.6.request,
                                                                 FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.elasmo,what.blk='block',Zn="Zone2")
  Elasmo.ktch.block10_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block10_27.6.request,
                                                                   FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.elasmo,what.blk='block10',Zn="Zone2")
  
  #Teleosts catch
  Teleost.ktch.block_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block_27.6.request,
                                                                  FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.scalies,what.blk='block',Zn="Zone2")
  Teleost.ktch.block10_27.6.request=fn.extract.Jodie.catch.Req.27.6(blk=ASL_block10_27.6.request,
                                                                    FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,SP=Spec.tab.1.scalies,what.blk='block10',Zn="Zone2")
  
  #Effort
  Effort.block_27.6.request=fn.extract.Jodie.effort.Req.27.6(blk=ASL_block_27.6.request,
                                                             FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,what.blk='block',Zn="Zone2")
  Effort.block10_27.6.request=fn.extract.Jodie.effort.Req.27.6(blk=ASL_block10_27.6.request,
                                                               FSHRY="JASDGDLF",YRS=Yrs.Jodie_27.6.request,what.blk='block10',Zn="Zone2")
  
  #Export
  #Attachment 1 request
  write.csv(Elasmo.ktch.block_27.6.request$BLKs,paste(hndls,"Elasmos.ktch_block_Attachment_1.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block_27.6.request$BLKs,paste(hndls,"Teleosts.ktch_block_Attachment_1.csv",sep=""),row.names=F)  
  write.csv(Effort.block_27.6.request$BLKs,paste(hndls,"Effort.block_Attachment_1.csv",sep=""),row.names=F)  
  
  #Attachment 2 request
  write.csv(Elasmo.ktch.block10_27.6.request$BLKs,paste(hndls,"Elasmos.ktch_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Elasmo.ktch.block10_27.6.request$Fishery,paste(hndls,"Elasmos.ktch_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Elasmo.ktch.block10_27.6.request$Zone,paste(hndls,"Elasmos.ktch_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
  write.csv(Teleost.ktch.block10_27.6.request$BLKs,paste(hndls,"Teleosts.ktch_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block10_27.6.request$Fishery,paste(hndls,"Teleosts.ktch_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Teleost.ktch.block10_27.6.request$Zone,paste(hndls,"Teleosts.ktch_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
  write.csv(Effort.block10_27.6.request$BLKs,paste(hndls,"Effort_block10_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Effort.block10_27.6.request$Fishery,paste(hndls,"Effort_JASDGDLF_Attachment_2.csv",sep=""),row.names=F)  
  write.csv(Effort.block10_27.6.request$Zone,paste(hndls,"Effort_Zone2_Attachment_2.csv",sep=""),row.names=F)  
  
}

#G 4.8.2 Demersal scalefish 2015/16
if(do.Jodies.scalies=="YES")
{
  
  a=subset(Data.daily,YEAR>=2015 & SPECIES%in%Scalefish.species,select=c(FisheryCode,YEAR,MONTH,METHOD,VESSEL,CONDITN,SPECIES,SNAME,LIVEWT.c))
  b=subset(Data.daily.incomplete,year>=2015 & species%in%Scalefish.species,select=c(fishery,year,month,method,vessel,conditn,species,sname1,livewt))
  names(a)=names(b)  
  a=subset(a,year== 2015 & month ==6)
  a=rbind(a,b)
  write.csv(a,paste(hndls,"Jodie.Demersal.scalies_15_16.csv",sep=""),row.names=F)
  rm(a,b)
}

#G 4.9 FishCUBE-------------------------------------------------------------------------
if(Extract.data.FishCUBE=="YES")
{
  hndl=handl_OneDrive("Analyses/Catch and effort/Data_Resquests")
  
  #add variables
  #Daily
  FishCUBE.daily$ExternalDataSourceName="Shark Daily Logbook"
  FishCUBE.daily$DailyorMonthly="D"
  FishCUBE.daily$LogBookPageNumber=sapply(strsplit(FishCUBE.daily$Same.return.SNo, " "), "[", 2)
  FishCUBE.daily$FishingSeason=NA
  FishCUBE.daily$YEAR=FishCUBE.daily$YEAR.c
  FishCUBE.daily$VesselRegistration=FishCUBE.daily$VESSEL
  FishCUBE.daily$GPSLatitude=FishCUBE.daily$LAT
  FishCUBE.daily$GPSLongitude=FishCUBE.daily$LONG
  FishCUBE.daily$Block10by10=FishCUBE.daily$block10
  FishCUBE.daily$Block60by60=FishCUBE.daily$blockxFC
  FishCUBE.daily$FishingMethod=FishCUBE.daily$METHOD
  FishCUBE.daily$LiveWeight=FishCUBE.daily$LIVEWT.c
  FishCUBE.daily$LiveWeight=with(FishCUBE.daily,ifelse(is.na(LiveWeight),LIVEWT.orgnl,LiveWeight))
  FishCUBE.daily$FisheryZone=with(FishCUBE.daily,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
  FishCUBE.daily$MONTH1=with(FishCUBE.daily,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
  FishCUBE.daily$day1=with(FishCUBE.daily,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
  FishCUBE.daily$date=with(FishCUBE.daily,paste(YEAR,"-",MONTH1,"-",day1,sep=""))
  FishCUBE.daily$fishery=FishCUBE.daily$FisheryCode
  
  these.FishCUBE=c("ExternalDataSourceName","DailyorMonthly","LogBookPageNumber","FishingSeason",
                   "FINYEAR","YEAR","MONTH","date","FDAYS","LatFC","LongFC",
                   "VesselRegistration","fishery","FisheryZone","Landing.Port","Block10by10","Block","Block60by60","Bioregion",
                   "FishingMethod","SPECIES","SNAME","RSCommonName","LiveWeight","BDAYS")
  
  FishCUBE=FishCUBE.daily[,match(c(these.FishCUBE,"RSSpeciesId"),names(FishCUBE.daily))]
  names(FishCUBE)[match(c("Landing.Port","LatFC","LongFC"),names(FishCUBE))]=c("LandingPort","Lat","Long")
  
  #add removed May and June 2006
  Data.daily.FC.2005_06$ExternalDataSourceName="Shark Daily Logbook"
  Data.daily.FC.2005_06$DailyorMonthly="D"
  Data.daily.FC.2005_06$LogBookPageNumber=Data.daily.FC.2005_06$DSNo
  Data.daily.FC.2005_06$FishingSeason=NA
  Data.daily.FC.2005_06$YEAR=Data.daily.FC.2005_06$year
  Data.daily.FC.2005_06$VesselRegistration=Data.daily.FC.2005_06$vessel
  Data.daily.FC.2005_06$GPSLatitude=Data.daily.FC.2005_06$LatDeg
  Data.daily.FC.2005_06$GPSLongitude=Data.daily.FC.2005_06$LongDeg
  Data.daily.FC.2005_06$Block10by10=Data.daily.FC.2005_06$block10
  Data.daily.FC.2005_06$Block60by60=Data.daily.FC.2005_06$blockxFC
  Data.daily.FC.2005_06$FishingMethod=Data.daily.FC.2005_06$method
  Data.daily.FC.2005_06$LiveWeight=Data.daily.FC.2005_06$livewt
  Data.daily.FC.2005_06$FisheryZone=with(Data.daily.FC.2005_06,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
  Data.daily.FC.2005_06$MONTH=Data.daily.FC.2005_06$month
  Data.daily.FC.2005_06$MONTH1=with(Data.daily.FC.2005_06,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
  Data.daily.FC.2005_06$day1=with(Data.daily.FC.2005_06,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
  Data.daily.FC.2005_06$FINYEAR=Data.daily.FC.2005_06$finyear
  Data.daily.FC.2005_06$FDAYS=Data.daily.FC.2005_06$fdays
  Data.daily.FC.2005_06$Lat=Data.daily.FC.2005_06$LatFC  
  Data.daily.FC.2005_06$Long=Data.daily.FC.2005_06$LongFC
  Data.daily.FC.2005_06$LandingPort=Data.daily.FC.2005_06$port
  Data.daily.FC.2005_06$SPECIES=Data.daily.FC.2005_06$species   
  Data.daily.FC.2005_06$SNAME=Data.daily.FC.2005_06$sname1
  Data.daily.FC.2005_06$BDAYS=Data.daily.FC.2005_06$bdays
  FishCUBE.data.daily.2005_06=Data.daily.FC.2005_06[,match(names(FishCUBE),names(Data.daily.FC.2005_06))]
  FishCUBE.data.daily.2005_06$date=as.character(FishCUBE.data.daily.2005_06$date)
  
  COLss=names(FishCUBE)
  for(cl in 1:length(COLss))
  {
    if(is.factor(FishCUBE[,cl])) FishCUBE[,cl]=as.character(FishCUBE[,cl])
    if(is.factor(FishCUBE.data.daily.2005_06[,cl])) FishCUBE.data.daily.2005_06[,cl]=as.character(FishCUBE.data.daily.2005_06[,cl])
  }
  FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.2005_06)
  
  
  #add latest incomplete year to daily
  if(nrow(Data.daily.incomplete)>0)
  {
    Data.daily.incomplete$ExternalDataSourceName="Shark Daily Logbook"
    Data.daily.incomplete$DailyorMonthly="D"
    Data.daily.incomplete$LogBookPageNumber=Data.daily.incomplete$DSNo
    Data.daily.incomplete$FishingSeason=NA
    Data.daily.incomplete$YEAR=Data.daily.incomplete$year
    Data.daily.incomplete$VesselRegistration=Data.daily.incomplete$vessel
    Data.daily.incomplete$GPSLatitude=Data.daily.incomplete$LatDeg
    Data.daily.incomplete$GPSLongitude=Data.daily.incomplete$LongDeg
    Data.daily.incomplete$Block10by10=Data.daily.incomplete$block10
    Data.daily.incomplete$Block60by60=Data.daily.incomplete$blockxFC
    Data.daily.incomplete$FishingMethod=Data.daily.incomplete$method
    Data.daily.incomplete$LiveWeight=Data.daily.incomplete$livewt
    Data.daily.incomplete$FisheryZone=with(Data.daily.incomplete,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
    Data.daily.incomplete$MONTH=Data.daily.incomplete$month
    Data.daily.incomplete$MONTH1=with(Data.daily.incomplete,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
    Data.daily.incomplete$day1=with(Data.daily.incomplete,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
    Data.daily.incomplete$FINYEAR=Data.daily.incomplete$finyear
    Data.daily.incomplete$FDAYS=Data.daily.incomplete$fdays
    Data.daily.incomplete$Lat=Data.daily.incomplete$LatFC  
    Data.daily.incomplete$Long=Data.daily.incomplete$LongFC
    Data.daily.incomplete$LandingPort=Data.daily.incomplete$port
    Data.daily.incomplete$SPECIES=Data.daily.incomplete$species   
    Data.daily.incomplete$SNAME=Data.daily.incomplete$sname1
    Data.daily.incomplete$BDAYS=Data.daily.incomplete$bdays
    #Data.daily.incomplete$rowid=NA
    
    FishCUBE.data.daily.incomplete=Data.daily.incomplete[,match(names(FishCUBE),names(Data.daily.incomplete))]
    FishCUBE.data.daily.incomplete$date=as.character(FishCUBE.data.daily.incomplete$date)
    
    for(cl in 1:length(COLss))if(is.factor(FishCUBE.data.daily.incomplete[,cl])) FishCUBE.data.daily.incomplete[,cl]=as.character(FishCUBE.data.daily.incomplete[,cl])
    FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.incomplete)
  }
  FishCUBE=subset(FishCUBE,!(SPECIES==99999 | is.na(SPECIES)))
  
  #add 0 catch
  if(nrow(Data.daily.FC.NA.sp)>0)
  {
    Data.daily.FC.NA.sp$ExternalDataSourceName="Shark Daily Logbook"
    Data.daily.FC.NA.sp$DailyorMonthly="D"
    Data.daily.FC.NA.sp$LogBookPageNumber=Data.daily.FC.NA.sp$DSNo
    Data.daily.FC.NA.sp$FishingSeason=NA
    Data.daily.FC.NA.sp$YEAR=Data.daily.FC.NA.sp$year
    Data.daily.FC.NA.sp$VesselRegistration=Data.daily.FC.NA.sp$vessel
    Data.daily.FC.NA.sp$GPSLatitude=Data.daily.FC.NA.sp$LatDeg
    Data.daily.FC.NA.sp$GPSLongitude=Data.daily.FC.NA.sp$LongDeg
    Data.daily.FC.NA.sp$Block10by10=Data.daily.FC.NA.sp$block10
    Data.daily.FC.NA.sp$Block60by60=Data.daily.FC.NA.sp$blockxFC
    Data.daily.FC.NA.sp$FishingMethod=Data.daily.FC.NA.sp$method
    Data.daily.FC.NA.sp$LiveWeight=Data.daily.FC.NA.sp$livewt
    Data.daily.FC.NA.sp$FisheryZone=with(Data.daily.FC.NA.sp,ifelse(GPSLongitude>=116.5 & GPSLongitude<= 116.923 & GPSLatitude < (-33),"Zone3",zone))
    Data.daily.FC.NA.sp$MONTH=Data.daily.FC.NA.sp$month
    Data.daily.FC.NA.sp$MONTH1=with(Data.daily.FC.NA.sp,ifelse(MONTH<10,paste("0",MONTH,sep=""),as.character(MONTH)))
    Data.daily.FC.NA.sp$day1=with(Data.daily.FC.NA.sp,ifelse(day<10,paste("0",day,sep=""),as.character(day)))
    Data.daily.FC.NA.sp$FINYEAR=Data.daily.FC.NA.sp$finyear
    Data.daily.FC.NA.sp$FDAYS=Data.daily.FC.NA.sp$fdays
    Data.daily.FC.NA.sp$Lat=Data.daily.FC.NA.sp$LatFC  
    Data.daily.FC.NA.sp$Long=Data.daily.FC.NA.sp$LongFC
    Data.daily.FC.NA.sp$LandingPort=Data.daily.FC.NA.sp$port
    Data.daily.FC.NA.sp$SPECIES=Data.daily.FC.NA.sp$species   
    Data.daily.FC.NA.sp$SNAME=Data.daily.FC.NA.sp$sname1
    Data.daily.FC.NA.sp$BDAYS=Data.daily.FC.NA.sp$bdays
    #Data.daily.FC.NA.sp$rowid=NA
    
    FishCUBE.data.daily.NA.sp=Data.daily.FC.NA.sp[,match(names(FishCUBE),names(Data.daily.FC.NA.sp))]
    FishCUBE.data.daily.NA.sp$date=as.character(FishCUBE.data.daily.NA.sp$date)
    
    for(cl in 1:length(COLss))if(is.factor(FishCUBE.data.daily.NA.sp[,cl])) FishCUBE.data.daily.NA.sp[,cl]=as.character(FishCUBE.data.daily.NA.sp[,cl])
    FishCUBE=rbind(FishCUBE,FishCUBE.data.daily.NA.sp)
  }
  
  
  #Monthly
  FishCUBE.monthly=subset(FishCUBE.monthly,TYPE.DATA=="monthly")    #keep only monthly records
  a=subset(Data.monthly.original,Same.return%in%unique(FishCUBE.monthly$Same.return),select=c(Same.return,FDAYS))
  a=a[!duplicated(a$Same.return),]
  FishCUBE.monthly=merge(FishCUBE.monthly,a,by="Same.return",all.x=T)
  FishCUBE.monthly$FishingMethod=FishCUBE.monthly$Fishing.method
  FishCUBE.monthly$LandingPort=FishCUBE.monthly$Landing.Port
  FishCUBE.monthly$LatFC=FishCUBE.monthly$LAT
  FishCUBE.monthly$LongFC=FishCUBE.monthly$LONG
  FishCUBE.monthly$YEAR=FishCUBE.monthly$YEAR.c
  FishCUBE.monthly$Block=FishCUBE.monthly$blockxFC
  FishCUBE.monthly$Block60by60=FishCUBE.monthly$blockxFC
  FishCUBE.monthly$FishingMethod=FishCUBE.monthly$METHOD
  FishCUBE.monthly$LiveWeight=FishCUBE.monthly$LIVEWT.c
  FishCUBE.monthly$VesselRegistration=FishCUBE.monthly$VESSEL
  FishCUBE.monthly$fishery=FishCUBE.monthly$FisheryCode
  
  FishCUBE.monthly$ExternalDataSourceName="FLAMSCAE"
  FishCUBE.monthly$DailyorMonthly="M"
  
  FishCUBE.monthly$LogBookPageNumber=NA
  FishCUBE.monthly$FishingSeason=NA
  FishCUBE.monthly$day=NA
  FishCUBE.monthly$date=NA
  FishCUBE.monthly$GPSLatitude=NA
  FishCUBE.monthly$GPSLongitude=NA
  FishCUBE.monthly$Block10by10=NA
  
  
  FishCUBE.monthly=FishCUBE.monthly[,match(c(these.FishCUBE,'HOURS.c','HOOKS','SHOTS.c','NETLEN.c'),names(FishCUBE.monthly))]
  names(FishCUBE.monthly)[match(c("Landing.Port","LatFC","LongFC",'HOURS.c','HOOKS','SHOTS.c','NETLEN.c'),
                                names(FishCUBE.monthly))]=c("LandingPort","Lat","Long",'Hours','Hooks','Shots','Netlen')
  
  
  #Fix RSSpeciesId for reapportioned sharks
  #Daily
  FishCUBE$RSSpeciesId=with(FishCUBE,ifelse(SPECIES==17001,63,
                                            ifelse(SPECIES==18003,70,
                                                   ifelse(SPECIES==17003,64,
                                                          ifelse(SPECIES==18007,74,
                                                                 ifelse(SPECIES==22999,943,RSSpeciesId))))))
  
  FishCUBE$RSCommonName=with(FishCUBE,ifelse(SPECIES==17001,"Gummy Shark",
                                             ifelse(SPECIES==18003,"Dusky Whaler",
                                                    ifelse(SPECIES==17003,"Whiskery Shark",
                                                           ifelse(SPECIES==18007,"Sandbar Shark",
                                                                  ifelse(SPECIES==22999,"Sharks",RSCommonName))))))
  
  FishCUBE$SNAME=with(FishCUBE,ifelse(SPECIES==17001,"Gummy Shark",
                                      ifelse(SPECIES==18003,"Bronze Whaler",
                                             ifelse(SPECIES==17003,"Whiskery Shark",
                                                    ifelse(SPECIES==18007,"Sandbar Shark",
                                                           ifelse(SPECIES==22999,"shark, other",SNAME))))))
  
  
  #Monthly
  # FishCUBE.monthly$RSSpeciesId=with(FishCUBE.monthly,ifelse(SPECIES==17001,63,
  #              ifelse(SPECIES==18003,70,
  #              ifelse(SPECIES==17003,64,
  #              ifelse(SPECIES==18007,74,
  #              ifelse(SPECIES==22999,943,RSSpeciesId))))))
  
  FishCUBE.monthly$RSCommonName=with(FishCUBE.monthly,ifelse(SPECIES==17001,"Gummy Shark",
                                                             ifelse(SPECIES==18003,"Dusky Whaler",
                                                                    ifelse(SPECIES==17003,"Whiskery Shark",
                                                                           ifelse(SPECIES==18007,"Sandbar Shark",
                                                                                  ifelse(SPECIES==22999,"Sharks",RSCommonName))))))
  
  FishCUBE.monthly$SNAME=with(FishCUBE.monthly,ifelse(SPECIES==17001,"Gummy Shark",
                                                      ifelse(SPECIES==18003,"Bronze Whaler",
                                                             ifelse(SPECIES==17003,"Whiskery Shark",
                                                                    ifelse(SPECIES==18007,"Sandbar Shark",
                                                                           ifelse(SPECIES==22999,"shark, other",SNAME))))))
  
  
  #change NA to blank
  #Daily
  FishCUBE$Lat=as.character(FishCUBE$Lat)
  FishCUBE$Long=as.character(FishCUBE$Long)
  FishCUBE$Lat=with(FishCUBE,ifelse(is.na(Lat),"",Lat))
  FishCUBE$Long=with(FishCUBE,ifelse(is.na(Long),"",Long))
  
  #Monthly
  FishCUBE.monthly$Lat=as.character(FishCUBE.monthly$Lat)
  FishCUBE.monthly$Long=as.character(FishCUBE.monthly$Long)
  FishCUBE.monthly$Lat=with(FishCUBE.monthly,ifelse(is.na(Lat),"",Lat))
  FishCUBE.monthly$Long=with(FishCUBE.monthly,ifelse(is.na(Long),"",Long))
  
  
  #export data
  
  #monthly
  hndl.FishCUBE="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/CAES Monthly Fisheries/02. Data/"
  system.time(fwrite(FishCUBE.monthly,paste(hndl.FishCUBE,"2. CAES Monthly Data - Shark - Corrected with R.csv",sep=""),row.names=F))
  
  
  #daily
  hndl.FishCUBE.daily="M:/CAEMaster/Commercial/FishCubeWA/Fisheries/Shark/02. Data/"
  system.time(fwrite(FishCUBE,paste(hndl.FishCUBE.daily,"2. Shark Daily Logbook Data - Corrected with R.csv",sep=""),row.names=F))
  #system.time(write.csv(FishCUBE,paste(hndl.FishCUBE,"Shark Data Extract for FishCubeWA.",Sys.Date(),".csv",sep=""),row.names=F))
  
  # library("openxlsx")
  # options(java.parameters = "-Xmx8000m")
  #  system.time(write.xlsx(FishCUBE,paste(hndl.FishCUBE,"Shark Data Extract for FishCubeWA.",Sys.Date(),".xlsx",sep=""), colNames = TRUE))
  
  #notify FishCube by email 
  send.email(TO=Email.FishCube,
             Subject="Annual data",
             Body= paste("Data corrections done and uploaded.",
                         "Data files stored in",hndl.FishCUBE,"and",
                         hndl.FishCUBE.daily),  
             Attachment=NULL)
}

#G 4.10 Steve's map for interview-------------------------------------------------------------------------
if (do.Steves=="YES")
{
  a=South.WA.long[1]:South.WA.long[2]
  b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
  
  GRID="Y"
  if(GRID=="Y")tiff(file=handl_OneDrive("Analyses/Catch and effort/Outputs/for Steve/Steve,map.grid.tiff"),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  if(GRID=="N")tiff(file=handl_OneDrive("Analyses/Catch and effort/Outputs/for Steve/Steve,map.NoGrid.tiff"),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  
  plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
  if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                               nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
  if(GRID=="Y")
  {
    axis(side = 1, at =LONGG, labels = F, tcl = 30.6)
    axis(side = 4, at = LATT, labels = F,tcl =34)
    par(new=T)
    plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    
  }
  axis(side = 1, at =LONGG, labels = LONGG, tcl = -.5)
  axis(side = 2, at = LATT, labels = LATT,tcl =-.5,las=1)
  
  text(116,-28.78,("Geraldton"),col="black", cex=1.1)
  points(114.83,-28.86,pch=19)
  text(116.73,Perth[2],("Perth"),col="black", cex=1.1)
  points(115.86,-31.95,pch=19)
  text(116.73,-33.55,("Bunbury"),col="black", cex=1.1)
  points(115.6,-33.55,pch=19)
  text(115.95,-34.085,("Augusta"),col="black", cex=1.1)
  points(115.15,-34.31,pch=19)
  text(117.7,-34.8,("Albany"),col="black", cex=1.1)
  points(117.8,-35,pch=19)
  text(122,-33.66,("Esperance"),col="black", cex=1.1)
  points(121.9,-33.86,pch=19)
  
  mtext(Lat.exp,side=2,line=2.75,las=3,cex=1.3)
  mtext(Lon.exp,side=1,line=2.75,cex=1.3)
  
  
  dev.off()
  
}

#G 4.11 EXTRA FIGURES FOR AMM 2013-------------------------------------------------------------------------
if(do.AMM.2013=="YES")
{
  setwd(paste(hndl,"/AMM_2013",sep=""))
  if(Use.Previos.Sofar=="YES")
  {
    fn.fig.AMM=function(GROUP,LAT1,LAT2,INT,INT2,THIS)
    {
      dat=subset(DAT,SPECIES%in%GROUP & LAT<=LAT1 & LAT >=LAT2)
      
      annual.catch.by.zone=aggregate(LIVEWT.c~FINYEAR+zone,data=dat,sum,na.rm=T)
      annual.catch.total=aggregate(LIVEWT.c~FINYEAR,data=dat,sum,na.rm=T)
      
      wide=reshape(annual.catch.by.zone,v.names="LIVEWT.c",timevar="zone",idvar="FINYEAR",direction="wide")
      
      #add previous years
      Prev.zn=Spec.catch.zone.pre.2013[,match(THIS,names(Spec.catch.zone.pre.2013))]
      names(Prev.zn)=names(wide)
      Prev.zn[,2:4]=Prev.zn[,2:4]*1000
      wide=rbind(Prev.zn,wide)
      
      #total 
      annual.catch.total=data.frame(FINYEAR=wide$FINYEAR,Total.catch=rowSums(wide[,2:4]))
      
      names(annual.catch.total)[match("FINYEAR",names(annual.catch.total))]="finyear"
      names(wide)[match("FINYEAR",names(wide))]="finyear"
      
      #plot
      fun.fig.SoFar(annual.catch.total,wide,1000,"Catch (tonnes live wt.)","Financial year",INT,INT2)
      
      return(Tot.Catch=annual.catch.total)
    }
  }
  
  Lista.especies=list(whiskery=c("Whiskery.catch.WC","Whiskery.catch.Zn1","Whiskery.catch.Zn2"),
                      gummy=c("Gummy.catch.WC","Gummy.catch.Zn1","Gummy.catch.Zn2"),
                      dusky=c("Dusky.catch.WC","Dusky.catch.Zn1","Dusky.catch.Zn2"),
                      sandbar=c("Sandbar.catch.WC","Sandbar.catch.Zn1","Sandbar.catch.Zn2"))
  SPec=list(17003,17001,c(18001,18003),18007)
  sp=c("whiskery","gummy","dusky","sandbar")
  Catch.range=list(c(175,225),c(350,450),c(200,300),c(120,120))
  
  STORE.Tot.Catch=vector('list',length=length(TARGETS))
  names(STORE.Tot.Catch)=sp
  for (i in 1:length(TARGETS))
  {
    jpeg(file=paste("Catch.zone.",sp[i],".jpeg",sep=""),width = 2400, height = 2400,units = "px", res = 300)
    STORE.Tot.Catch[[i]]=fn.fig.AMM(SPec[[i]],-27,-40,100,50,c("Fin.yr",Lista.especies[[i]]))
    abline(Catch.range[[i]][1],0,col=2)
    abline(Catch.range[[i]][2],0,col=2)
    dev.off()  
  }
  
}

#G 4.12 Catch and effort Perth metro closures-------------------------------------------------------------------------
#Closure
CLOSED=c("31120","31130","31140","31150","32120","32130","32140","32150")
if(do.Perth.metro.closure=="YES")
{
  #1. Aggregate effort by year and blockfor WCDGDLF
  
  #Daily
  a=subset(Effort.daily,zone=="West")
  #km gn days 
  if(Use.Date=="NO")  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~ID+vessel+blockx+finyear,data=a,max,na.rm=T)
  if(Use.Date=="YES") Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~date+vessel+blockx+finyear,data=a,max,na.rm=T) 
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~blockx+finyear,data=Effort.Perth.metro,sum,na.rm=T)
  
  #km gn hours
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~finyear+vessel+ID+blockx,data=a,max,na.rm=T)
  #if(Use.Date=="YES")Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~blockx+finyear,data=Effort.Perth.metro.hrs,sum,na.rm=T)
  Effort.Perth.metro.daily=merge(Effort.Perth.metro.hrs,Effort.Perth.metro,by=c("blockx","finyear"))
  
  #Monthly
  a=subset(Effort.monthly,zone=="West" & !FINYEAR%in%Daily.l.years)
  #km gn days             
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+BLOCKX+FINYEAR,data=a,max,na.rm=T)
  Effort.Perth.metro=aggregate(Km.Gillnet.Days.c~BLOCKX+FINYEAR,data=Effort.Perth.metro,sum,na.rm=T)
  
  #km gn hours
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~MONTH+BLOCKX+VESSEL+BLOCKX+FINYEAR,data=a,max,na.rm=T)
  Effort.Perth.metro.hrs=aggregate(Km.Gillnet.Hours.c~BLOCKX+FINYEAR,data=Effort.Perth.metro.hrs,sum,na.rm=T)
  Effort.Perth.metro.monthly=merge(Effort.Perth.metro.hrs,Effort.Perth.metro,by=c("BLOCKX","FINYEAR"))
  
  Effort.Perth.metro.daily$finyear=as.character(Effort.Perth.metro.daily$finyear)
  Effort.Perth.metro.monthly$FINYEAR=as.character(Effort.Perth.metro.monthly$FINYEAR)
  names(Effort.Perth.metro.monthly)=names(Effort.Perth.metro.daily)
  Effort.Perth.metro=rbind(Effort.Perth.metro.monthly,Effort.Perth.metro.daily)
  
  
  
  #2. Aggregate catch by year and blockfor WCDGDLF
  KTCH.Perth.metro=aggregate(LIVEWT.c~BLOCKX+FINYEAR+SPECIES+SNAME,data=subset(Data.monthly,zone=="West"),sum,na.rm=T)
  
  #3. Extract effort Metro closure
  Effort.Perth.metro.closure=subset(Effort.Perth.metro,blockx%in%as.numeric(CLOSED))
  
  #4. Sharks and Demersal scalefish suite component
  KTCH.Perth.metro.shrk=subset(KTCH.Perth.metro,SPECIES%in%Elasmo.species)
  KTCH.Perth.metro.scale=subset(KTCH.Perth.metro,SPECIES%in%Suite)
  
  #5. Catch from Metro Closure
  KTCH.Perth.metro.closure.shrk=subset(KTCH.Perth.metro,SPECIES%in%Elasmo.species & BLOCKX%in%CLOSED)
  KTCH.Perth.metro.closure.scale=subset(KTCH.Perth.metro,SPECIES%in%Suite  & BLOCKX%in%CLOSED)
  
  #6. Aggregate by year
  Eff=aggregate(Km.Gillnet.Days.c~finyear,Effort.Perth.metro,sum)
  Eff.closure=aggregate(Km.Gillnet.Days.c~finyear,Effort.Perth.metro.closure,sum)
  
  Agg.shark=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.shrk,sum)
  Agg.scale=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.scale,sum)
  
  Agg.shark.closed=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.closure.shrk,sum)
  Agg.scale.closed=aggregate(LIVEWT.c~FINYEAR,KTCH.Perth.metro.closure.scale,sum)
  
  
  #Export
  fn.plt=function(a,b,d,e)
  {
    if(length(e)<length(d))
    {
      id=which(!d%in%e)
      test=data.frame(e,b)
      test1=data.frame(d[id],0)
      names(test1)=names(test)
      b1=rbind(test,test1)
      b=b1[,2]
    }
    plot(1:length(a),a,ylab="",xlab="",type='l',lwd=2,ylim=c(0,max(a)),xaxt='n',cex.axis=1.35)    
    lines(1:length(b),b,lwd=2,lty=3)
    axis(1,1:length(a),F,tck=-0.02)
    if(length(a)>20)axis(1,seq(1,length(a),10),d[seq(1,length(a),10)],tck=-0.04,cex.axis=1.5)
    if(length(a)<20)axis(1,seq(1,length(a),2),d[seq(1,length(a),2)],tck=-0.04,cex.axis=1.5)
  }
  tiff(file=paste(hndl,"/Metro_closure_catch_effort.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.5,.6,.1,.1),las=1,mgp=c(2,0.65,0))  
  
  #Effort
  fn.plt(Eff$Km.Gillnet.Days.c/1000,Eff.closure$Km.Gillnet.Days.c/1000,Eff$finyear,Eff.closure$finyear)
  mtext("Effort (1000 km gn d)",2,line=2.7,font=1,las=0,cex=1.45)
  legend("topright",c("WCDGDLF","Metro closure"),bty='n',lty=c(1,3),col=1,lwd=2,cex=2)
  
  #Shark catch
  fn.plt(Agg.shark$LIVEWT.c/1000,Agg.shark.closed$LIVEWT.c/1000,Agg.shark$FINYEAR,Agg.shark.closed$FINYEAR)
  mtext("Catch (tonnes)",2,line=2.7,font=1,las=0,cex=1.5)
  legend("topright",c("Shark"),bty='n',cex=2.5)
  
  #Scalefish suite catch
  fn.plt(Agg.scale$LIVEWT.c/1000,Agg.scale.closed$LIVEWT.c/1000,Agg.scale$FINYEAR,Agg.scale.closed$FINYEAR)
  mtext("Catch (tonnes)",2,line=2.7,font=1,las=0,cex=1.5)
  legend("topright",c("Demersal scalefish suite"),bty='n',cex=2.5)
  
  mtext("Financial year",side=1,line=-1.5,font=1,las=0,cex=1.5,outer=T) 
  dev.off()
  
  
  #write.csv(Effort.Perth.metro,paste(hndl,"/Effort.Perth.metro.csv",sep=""),row.names=F)
  #write.csv(KTCH.Perth.metro,paste(hndl,"/Catch.Perth.metro.csv",sep=""),row.names=F)
}

#G 4.13  Heather's request proportion catch by LL and GN-------------------------------------------------------------------------
if(do.Heather.request=="YES")
{
  Fn.prop.ktch.ll.gn=function(datGN,datLL)
  {
    no.scalies=paste(1975:1997,substr(1976:1998,3,4),sep="-")
    datGN=subset(datGN,!FINYEAR%in%no.scalies)   #no scalies reported before 1998-99
    datLL=subset(datLL,!FINYEAR%in%no.scalies)
    datGN$scalie=with(datGN,ifelse(SPECIES%in%Scalefish.species,"YES" ,"NO"))
    datLL$scalie=with(datLL,ifelse(SPECIES%in%Scalefish.species,"YES" ,"NO"))
    TabGN=with(datGN,table(FINYEAR, scalie))
    TabLL=with(datLL,table(FINYEAR, scalie))
    
    TabGN=as.data.frame.matrix(TabGN)
    TabLL=as.data.frame.matrix(TabLL)
    
    TabGN$Tot=TabGN$YES+TabGN$NO
    TabLL$Tot=TabLL$YES+TabLL$NO
    
    TabGN$Prop.scale=with(TabGN, YES/ Tot)
    TabLL$Prop.scale=with(TabLL, YES/ Tot)
    
    TabGN$Prop.scale=with(TabGN,ifelse(Tot<100,NA,Prop.scale))
    TabLL$Prop.scale=with(TabLL,ifelse(Tot<100,NA,Prop.scale))
    
    plot(TabGN$Prop.scale,type='l',col=2,lwd=2,ylab="Proportion of scalefish",xlab="",
         ylim=c(0,1),xaxt='n',cex.lab=1.5)
    lines(TabLL$Prop.scale,lwd=2,col=4)
    legend("topleft",c("Gillnet","Longline"),lty=1,lwd=2,col=c(2,4),bty='n',cex=1.25)
    axis(1,1:nrow(TabGN),F)
    axis(1,seq(1,nrow(TabGN),5),rownames(TabGN)[seq(1,nrow(TabGN),5)])
  }
  
  tiff(file=paste(hndl,"/Heather_prop_LL_GN.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,1),mar=c(3,4,1,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2.5,0.8,0))
  
  Fn.prop.ktch.ll.gn(subset(Data.monthly.GN, zone=="West"),
                     subset(Data.monthly.LL, zone=="West"))
  mtext("West Coast",3,cex=1.5,line=-2)
  
  Metro=Fn.prop.ktch.ll.gn(subset(Data.monthly.GN, BLOCKX%in%as.numeric(CLOSED)), 
                           subset(Data.monthly.LL, BLOCKX%in%as.numeric(CLOSED)))
  mtext("Closure",3,cex=1.5,line=-2)
  mtext("Financial year",1,outer=T,cex=1.5) 
  dev.off()  
}

#G 4.14 Annual scalefish Catch for Dave F.-------------------------------------------------------------------------
if(do.Dave.F=="YES")
{
  solo.monthly=subset(FINYEAR.monthly,!FINYEAR.monthly%in%FINYEAR.daily)
  a=subset(Data.monthly,FINYEAR%in%solo.monthly & SPECIES%in%Scalefish.species & Bioregion%in%c("SC","WC"))
  if(nrow(a)>0)
  {
    KTCH.scalies_monthly=aggregate(LIVEWT.c~BLOCKX+Bioregion+FINYEAR+SPECIES+SNAME,data=a,sum,na.rm=T)
    write.csv(KTCH.scalies_monthly,paste(hndl,"/Dave F/Dave_F_scalefish_monthly.csv",sep=""),row.names=F)
    rm(a)
  }
  KTCH.scalies_daily=aggregate(LIVEWT.c~block10+Bioregion+FINYEAR+SPECIES+SNAME,data=subset(Data.daily,SPECIES%in%Scalefish.species),sum,na.rm=T)
  write.csv(KTCH.scalies_daily,paste(hndl,"/Dave F/Dave_F_scalefish_daily.csv",sep=""),row.names=F)
  
  #email data
  send.email(TO="David.Fairclough@dpird.wa.gov.au",
             Subject="Annual data",
             Body= "Please find attached the annual monthly records",  
             Attachment=paste(hndl,"/Dave F/Dave_F_scalefish_monthly.csv",sep=""))
  send.email(TO="David.Fairclough@dpird.wa.gov.au",
             Subject="Annual data",
             Body= "Please find attached the annual daily records",  
             Attachment=paste(hndl,"/Dave F/Dave_F_scalefish_daily.csv",sep=""))
  
  do.table.all.yrs="NO"
  if(do.table.all.yrs=="YES")
  {
    This.fish=match(c(Blue_morwong,Blue_groper,West_Australian_dhufish,Pink_snapper,Boarfishes,Samsonfish,Redfishes,Mulloway,Sweetlips,
                      Baldchin_groper),Scalefish.species)
    This.fish=This.fish[!is.na(This.fish)]
    DAT.D=subset(Data.daily,SPECIES%in%Scalefish.species&Bioregion%in%c("SC","WC"))
    DAT.D$SPECIES.agg=with(DAT.D,
                           ifelse(SPECIES%in%Blue_morwong,"Blue_morwong",
                                  ifelse(SPECIES%in%Blue_groper,"Blue_groper",
                                         ifelse(SPECIES%in%West_Australian_dhufish,"West_Australian_dhufish",
                                                ifelse(SPECIES%in%Pink_snapper,"Pink_snapper",
                                                       ifelse(SPECIES%in%Boarfishes,"Boarfishes",
                                                              ifelse(SPECIES%in%Samsonfish,"Samsonfish",
                                                                     ifelse(SPECIES%in%Redfishes,"Redfishes",
                                                                            ifelse(SPECIES%in%Mulloway,"Mulloway",
                                                                                   ifelse(SPECIES%in%Sweetlips,"Sweetlips",
                                                                                          ifelse(SPECIES%in%Baldchin_groper,"Baldchin_groper",
                                                                                                 ifelse(SPECIES%in%Scalefish.species[-This.fish],"Other_scalefish",
                                                                                                        NA))))))))))))
    
    TABLA=aggregate(LIVEWT.c~Bioregion+FINYEAR+SPECIES.agg,data=DAT.D,sum,na.rm=T)
    Tab.suite=aggregate(LIVEWT.c~Bioregion+FINYEAR,data=subset(DAT.D,SPECIES%in%Suite),sum,na.rm=T)
    Tab.suite$SPECIES.agg="Suite"
    TABLA=rbind(TABLA,Tab.suite)
    wide <- reshape(TABLA,v.names="LIVEWT.c",timevar="Bioregion",idvar=c("FINYEAR","SPECIES.agg"),direction="wide")
    write.csv(wide,paste(hndl,"/Dave F/Daily_table_all_yrs_bioregion.csv",sep=""),row.names=F)
    
    TABLA=aggregate(LIVEWT.c~zone+FINYEAR+SPECIES.agg,data=DAT.D,sum,na.rm=T)
    Tab.suite=aggregate(LIVEWT.c~zone+FINYEAR,data=subset(DAT.D,SPECIES%in%Suite),sum,na.rm=T)
    Tab.suite$SPECIES.agg="Suite"
    TABLA=rbind(TABLA,Tab.suite)
    wide <- reshape(TABLA,v.names="LIVEWT.c",timevar="zone",idvar=c("FINYEAR","SPECIES.agg"),direction="wide")
    write.csv(wide,paste(hndl,"/Dave F/Daily_table_all_yrs_zone.csv",sep=""),row.names=F)
    
    rm(DAT.D,wide)
  }
}

#G 4.15 Annual scalefish Catch for Jeff N.-------------------------------------------------------------------------
if(do.Jeff.N=="YES")
{
  a=subset(Data.monthly,SPECIES%in%Scalefish.species & Bioregion=="SC")
  
  Demersal.scale=c(224901,228002,258000,258004,258005,258006,264004,287000,288000,311000,
                   311005,311006,311078,311100,311152,311170,320000,320000,330001,346911,
                   346914,353001,361004,367000,369002,377004,384002,384904,384999,386000,
                   439002,445001)
  
  Jeff.list=list(Demersal.scale=Demersal.scale,Redfish=c(258000,258004,258006),
                 Pink.snapper=353001,Blue.morwong=377004,Blue.groper=384002,Hapuku=311006) 
  
  b=subset(a,SPECIES%in%Jeff.list$Demersal.scale)
  Dummy1=aggregate(LIVEWT.c~FINYEAR,data=b,sum,na.rm=T)
  names(Dummy1)[2]="Demersal.scalefish"
  
  a$SP.jeff=with(a,
                 ifelse(SPECIES%in%Jeff.list$Redfish,"Redfish",
                        ifelse(SPECIES%in%Jeff.list$Pink.snapper,"Pink.snapper",
                               ifelse(SPECIES%in%Jeff.list$Blue.morwong,"Blue.morwong",
                                      ifelse(SPECIES%in%Jeff.list$Blue.groper,"Blue.groper",
                                             ifelse(SPECIES%in%Jeff.list$Hapuku,"Hapuku",NA
                                             ))))))
  a=subset(a,!is.na(SP.jeff))
  
  Dummy2=aggregate(LIVEWT.c~FINYEAR+SP.jeff,data=a,sum,na.rm=T)
  Dummy2=reshape(Dummy2, v.names = "LIVEWT.c", idvar = "FINYEAR",
                 timevar ="SP.jeff" , direction = "wide")
  names(Dummy2)[2:ncol(Dummy2)]=substr(names(Dummy2)[2:ncol(Dummy2)],10,50)
  Dummy2[is.na(Dummy2)]=0
  Dummy=merge(Dummy1,Dummy2,by="FINYEAR")
  OuTT=paste(hndl,"/Jeff N/","Annual_ktch(kg).csv",sep="")  
  write.csv(Dummy,OuTT,row.names=F)
  
  #email data
  send.email(TO="Jeffrey.Norriss@dpird.wa.gov.au",
             Subject="Annual data",
             Body= "Please find attached the annual catch records",  
             Attachment=OuTT)
  
  rm(a)
  
}

#G 4.16 Annual scalefish Catch for Paul L.-------------------------------------------------------------------------
if(do.Paul.L=="YES")
{
  Paul.list=list(Samson.fish=337007)  
  a=subset(Data.monthly,SPECIES%in%Scalefish.species)
  
  Dummy=aggregate(LIVEWT.c~FINYEAR+Bioregion+SPECIES,data=subset(a,SPECIES%in%Paul.list[[1]]),sum,na.rm=T)
  Dummy$SPECIES="Samson.fish"
  
  OuTT=paste(hndl,"/Paul L/","Annual_ktch(kg).csv",sep="")
  write.csv(Dummy,OuTT,row.names=F)
  
  #email data
  send.email(TO="paul.lewis@dpird.wa.gov.au",
             Subject="Annual data",
             Body= "Please find attached the annual daily records",  
             Attachment=OuTT)
  rm(a)
}

#G 4.17 Hammerheads listing-------------------------------------------------------------------------
if(do.Hammerheads=="YES")
{
  KTCH=aggregate(LIVEWT.c/1000~FINYEAR,Hammerheads,sum)
  names(KTCH)[2]="LIVEWT.c.ton"
  MissYr=FINYrs[which(!FINYrs%in%KTCH$FINYEAR)]
  KTCH=rbind(KTCH,data.frame(FINYEAR=MissYr,LIVEWT.c.ton=NA))
  KTCH=KTCH[order(KTCH$FINYEAR),]
  Ny=length(FINYrs)
  tiff(file=paste(hndl,"/Annual.Hammerhead.Catch.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mgp=c(2.3,.7,0),las=1)
  plot(1:Ny,KTCH$LIVEWT.c.ton,xlab="Financial year",ylab="Total catch (tonnes)",xaxt='n',
       type="b",pch=19,cex=2,cex.lab=2,cex.axis=1.25)
  axis(1,1:Ny,F,tck=-0.01)
  axis(1,seq(1,Ny,5),F,tck=-0.015)
  axis(1,seq(1,Ny,10),FINYrs[seq(1,Ny,10)],tck=-0.02,cex.axis=1.25)
  dev.off()
}

#G 4.18 ABARES (James Woodhams) -------------------------------------------------------------------------
if(do.ABARES="YES")
{
  #Catch by species and financial year for 2006-07 to 2014-15
  Abr.yrs=c("2006-07","2007-08","2008-09","2009-10","2010-11",
            "2011-12","2012-13","2013-14","2014-15")
  D.mon=subset(Data.monthly,FINYEAR%in%Abr.yrs & SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c))
  D.mon.nor=subset(Data.monthly.north,FINYEAR%in%Abr.yrs & SPECIES%in%Elasmo.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c))
  D.abr=rbind(D.mon,D.mon.nor)
  SP.abare=D.abr[!duplicated(D.abr$SPECIES),match(c("SPECIES","SNAME"),names(D.abr))]
  Ktch.ab=aggregate(LIVEWT.c~FINYEAR+SPECIES,D.abr,sum)
  Ktch.ab=merge(Ktch.ab,SP.abare,by="SPECIES",all.x=T)
  Scien.names.abrs=read.csv(handl_OneDrive("Data/Comm.sp.Scientific.csv"))
  Scien.names.abrs=subset(Scien.names.abrs,select=c(SPECIES,Scientific.name))
  Ktch.ab=merge(Ktch.ab,Scien.names.abrs,by="SPECIES",all.x=T)
  
  
  #Effort by financial year for TDGDLF for 2006-07 to 2014-15
  Eff.ab_1000_km_gn_d.all=subset(Total.effort.days.monthly,FINYEAR%in%Abr.yrs)
  Eff.ab_1000_km_gn_d.zn=subset(Total.effort.zone.days.monthly,FINYEAR%in%Abr.yrs)
  
  #export
  write.csv(Ktch.ab,paste( hndl,"/Ktch.ab.csv",sep=''),row.names=F)
  write.csv(Eff.ab_1000_km_gn_d.all,paste( hndl,"/Eff.ab_1000_km_gn_d.all.csv",sep=""),row.names=F)
  write.csv(Eff.ab_1000_km_gn_d.zn,paste( hndl,"/Eff.ab_1000_km_gn_d.zn.csv",sep=""),row.names=F)
}

#G 4.19 Carlie Telfer -------------------------------------------------------------------------
if(do.Carlie.Telfer=="YES")
{  
  hndls=paste(hndl,"/Carlie_Telfer/",sep="")
  Carlie.Telfer.Eff.km.gn.d=aggregate(Km.Gillnet.Days.c~block10+Fishing_yr,Jodie.Effort,sum)
  Carlie.Telfer.Eff.km.gn.h=aggregate(Km.Gillnet.Hours.c~block10+Fishing_yr,Jodie.Effort,sum)
  
  write.csv(Carlie.Telfer.Eff.km.gn.d,paste(hndls,"Eff.km.gn.d.csv",sep=""),row.names=F)
  write.csv(Carlie.Telfer.Eff.km.gn.h,paste(hndls,"Eff.km.gn.h.csv",sep=""),row.names=F)
  
}

#G 4.20 Cetacean Interactions-------------------------------------------------------------------------
if(do.Cetacean.Inter=="YES")
{  
  hndls=paste(hndl,"/Cetacean_Inter/",sep="")
  
  #confidential=2
  confidential=0  #provide all records to Rory, so he can decide
  
  #GILLNETS
  KTCH.ly=subset(Data.monthly,METHOD=="GN" & Estuary=="NO" & NETLEN.c>=100 & LAT<=(-26),
                 select=c(Same.return,FINYEAR,MONTH,VESSEL,METHOD,
                          SPECIES,SNAME,LIVEWT.c,BLOCKX,LAT,LONG))
  
  KTCH.ly_daily=subset(Data.daily,METHOD=="GN" & Estuary=="NO" & netlen.c>=100 & LAT<=(-26),
                       select=c(Same.return,FINYEAR,MONTH,VESSEL,METHOD,zone,
                                block10,BLOCKX,SPECIES,SNAME,LIVEWT.c,LAT,LONG))
  
  daily.yrs=unique(Effort.daily$finyear)
  
  #monthly
  SHOTS_c=unique(KTCH.ly$Same.return)
  
  KTCH.ly=subset(KTCH.ly,!FINYEAR%in%daily.yrs)
  KTCH.ly$N=1
  KTCH.ly.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+VESSEL,data=KTCH.ly,max,na.rm=T)
  KTCH.ly.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX,data=KTCH.ly,sum,na.rm=T)
  Remove.these=subset(KTCH.ly.vessels,N<=confidential)
  
  KTCH.ly=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX,data=KTCH.ly,sum,na.rm=T)
  KTCH.ly$Fishery="TDGDLF"
  KTCH.ly$Gear="gillnet"
  names(KTCH.ly)[match("LIVEWT.c",names(KTCH.ly))]="Catch(kg)"
  
  if(nrow(Remove.these)>0)Remove.these$dummy=with(Remove.these,paste(FINYEAR, MONTH, BLOCKX))
  KTCH.ly$dummy=with(KTCH.ly,paste(FINYEAR, MONTH, BLOCKX))
  if(nrow(Remove.these)>0)KTCH.ly=subset(KTCH.ly,!dummy%in%Remove.these$dummy)
  KTCH.ly=KTCH.ly[,-match("dummy",names(KTCH.ly))]
  
  write.csv(KTCH.ly,paste(hndls,"Catch_monthly_records_TDGDLF.csv",sep=""),row.names=F)
  
  #effort
  Eff.ly=subset(Effort.monthly,!FINYEAR%in%daily.yrs & METHOD=="GN")
  Eff.ly=subset(Eff.ly,Same.return%in%SHOTS_c)
  Eff.ly=aggregate(Km.Gillnet.Days.c~Same.return+FINYEAR+MONTH+BLOCKX+VESSEL,data=Eff.ly,max,na.rm=T) #remove duplicates
  Eff.ly=aggregate(Km.Gillnet.Days.c~FINYEAR+MONTH+BLOCKX+VESSEL,data=Eff.ly,sum,na.rm=T)
  Eff.ly$Fishery="TDGDLF"
  Eff.ly$Gear="gillnet"
  names(Eff.ly)[match("Km.Gillnet.Days.c",names(Eff.ly))]="Effort(Km.Gillnet.Days)"
  
  if(nrow(Remove.these)>0)
  {
    Eff.ly$dummy=with(Eff.ly,paste(FINYEAR, MONTH, BLOCKX))
    Eff.ly=subset(Eff.ly,!dummy%in%Remove.these$dummy)
    Eff.ly=Eff.ly[,-match("dummy",names(Eff.ly))]
    Percent.removed.monthly=100*nrow(Remove.these)/nrow(KTCH.ly.vessels)
    write.csv(Percent.removed,paste(hndls,"Percent.removed.csv",sep=""),row.names=F)
  }
  
  write.csv(Eff.ly,paste(hndls,"Efforty_monthly_records_TDGDLF.csv",sep=""),row.names=F)
  
  #daily
  #note: if using block10, 98% of the data had less than 3 vessels
  #       if using BLOCKX, 77% of the data had less than 3 vessels, hence used ZONE
  
  SHOTS_c=unique(KTCH.ly_daily$Same.return)
  
  KTCH.ly_daily$N=1
  KTCH.ly_daily.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+block10+VESSEL,data=KTCH.ly_daily,max,na.rm=T)
  KTCH.ly_daily.vessels=aggregate(N~FINYEAR+MONTH+BLOCKX+block10,data=KTCH.ly_daily.vessels,sum,na.rm=T)
  if(nrow(Remove.these)>0)Remove.these_daily=subset(KTCH.ly_daily.vessels,N<=confidential)
  
  KTCH.ly_daily=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+block10+VESSEL,data=KTCH.ly_daily,sum,na.rm=T)
  KTCH.ly_daily$Fishery="TDGDLF"
  KTCH.ly_daily$Gear="gillnet"
  names(KTCH.ly_daily)[match("LIVEWT.c",names(KTCH.ly_daily))]="Catch(kg)"
  
  if(nrow(Remove.these)>0)Remove.these_daily$dummy=with(Remove.these_daily,paste(FINYEAR, MONTH,block10,VESSEL))
  KTCH.ly_daily$dummy=with(KTCH.ly_daily,paste(FINYEAR, MONTH, block10,VESSEL))
  if(nrow(Remove.these)>0)KTCH.ly_daily=subset(KTCH.ly_daily,!dummy%in%Remove.these_daily$dummy)
  KTCH.ly_daily=KTCH.ly_daily[,-match("dummy",names(KTCH.ly_daily))]
  
  write.csv(KTCH.ly_daily,paste(hndls,"Catch_daily_records_TDGDLF.csv",sep=""),row.names=F)
  
  
  #effort
  Eff.ly_daily=subset(Effort.daily,Same.return%in%SHOTS_c & method=="GN")
  Use.Date="YES"    #Rory's approach (aggregating by DATE)
  #Use.Date="NO"     # aggregating by SNo and DSNo 
  if(Use.Date=="NO")  Eff.ly_daily=aggregate(Km.Gillnet.Days.c~ID+Same.return+finyear+month+blockx+block10+vessel,data=Eff.ly_daily,max,na.rm=T)
  if(Use.Date=="YES") Eff.ly_daily=aggregate(Km.Gillnet.Days.c~date+Same.return+finyear+month+blockx+block10+vessel,data=Eff.ly_daily,max,na.rm=T) 
  Eff.ly_daily=aggregate(Km.Gillnet.Days.c~finyear+month+blockx+block10+vessel,data=Eff.ly_daily,sum,na.rm=T)
  Eff.ly_daily$Fishery="TDGDLF"
  Eff.ly_daily$Gear="gillnet"
  names(Eff.ly_daily)[match("Km.Gillnet.Days.c",names(Eff.ly_daily))]="Effort(Km.Gillnet.Days)"
  
  if(nrow(Remove.these)>0)
  {
    Eff.ly_daily$dummy=with(Eff.ly_daily,paste(finyear, month, block10,vessel))
    Eff.ly_daily=subset(Eff.ly_daily,!dummy%in%Remove.these_daily$dummy)
    Eff.ly_daily=Eff.ly_daily[,-match("dummy",names(Eff.ly_daily))]
    Percent.removed.daily=100*nrow(Remove.these_daily)/nrow(KTCH.ly_daily.vessels)
    write.csv(Percent.removed.daily,paste(hndls,"Percent.removed.daily.csv",sep=""),row.names=F)
  }
  write.csv(Eff.ly_daily,paste(hndls,"Efforty_daily_records_TDGDLF.csv",sep=""),row.names=F)
  
  
}

#G 4.21 Scalefish catch by year for metro closure-------------------------------------------------------------------------
if(do.Jodies.metro.closure=="YES")
{
  this.yr.Jodie=c("2002-03","2003-04","2004-05","2005-06","2006-07")
  Day.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,
                    data=subset(Data.daily,SPECIES%in%Scalefish.species & LAT<=(-31)& LAT>=(-33) & 
                                  zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)
  Mn.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,
                   data=subset(Data.monthly,SPECIES%in%Scalefish.species & LAT<=(-31)& LAT>=(-33) & 
                                 zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)  
  
  Day.Jod$METHOD=as.character(Day.Jod$METHOD) 
  Mn.Jod$BLOCKX=as.numeric(Mn.Jod$BLOCKX) 
  Jod=rbind(Day.Jod,Mn.Jod)
  Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,data=Jod,sum,na.rm=T)  
  write.csv(Jod,paste(hndl,"/Scalefish_metro_closure_catch_TDGDLF.csv",sep=""),row.names=F)
  
}

#G 4.22 Shark and ray catch by year for metro closure-------------------------------------------------------------------------
if(do.Jodies.metro.closure=="YES")
{
  this.yr.Jodie=c("2002-03","2003-04","2004-05","2005-06","2006-07")
  Day.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+VESSEL,
                    data=subset(Data.daily,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                                  zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)
  Mn.Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+VESSEL,
                   data=subset(Data.monthly,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                                 zone=="West" & FINYEAR%in%this.yr.Jodie),sum,na.rm=T)  
  
  SPE.nms.day=subset(Data.daily,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                       zone=="West" & FINYEAR%in%this.yr.Jodie,select=c(SPECIES,SNAME))
  SPE.nms.day=SPE.nms.day[!duplicated(SPE.nms.day$SPECIES),]
  Day.Jod=merge(Day.Jod,SPE.nms.day,by="SPECIES",all.x=T)
  
  SPE.nms.Mn=subset(Data.monthly,SPECIES%in%Elasmo.species & LAT<=(-31)& LAT>=(-33) & 
                      zone=="West" & FINYEAR%in%this.yr.Jodie,select=c(SPECIES,SNAME))
  SPE.nms.Mn=SPE.nms.Mn[!duplicated(SPE.nms.Mn$SPECIES),]
  Mn.Jod=merge(Mn.Jod,SPE.nms.Mn,by="SPECIES",all.x=T)
  
  Day.Jod$METHOD=as.character(Day.Jod$METHOD) 
  Mn.Jod$BLOCKX=as.numeric(Mn.Jod$BLOCKX) 
  
  Jod=rbind(Day.Jod,Mn.Jod)
  
  Jod=aggregate(LIVEWT.c~FINYEAR+MONTH+BLOCKX+METHOD+SPECIES+SNAME+VESSEL,data=Jod,sum,na.rm=T)  
  write.csv(Jod,paste(hndl,"/Elasmobranch_metro_closure_catch_TDGDLF.csv",sep=""),row.names=F)  
}

#G 4.23 Total shark fin by year for Tim Nicholas-------------------------------------------------------------------------
if(do.Tim_N.fin=="YES")
{
  DAT=subset(Data.monthly,FINYEAR%in%FINYEAR.daily & METHOD%in%c("GN","LL") & Estuary=="NO" &
               LAT<=TDGDLF.lat.range[1] & LAT >=TDGDLF.lat.range[2])
  #total fins   
  Fin.Weight=subset(DAT,SPECIES%in%Elasmo.species,select=c(SPECIES,FINYEAR,LIVEWT.c))
  Fin.Weight=subset(Fin.Weight,!(SPECIES%in%c(13000,20000,31000))) #remove no fin species (wobbies, spurdog,rays except shovelnose)
  Fin.Weight$Total.catch.fin=Percent.fin.of.livewt*Fin.Weight$LIVEWT.c
  Total.catch.fin=aggregate(Total.catch.fin~FINYEAR,Fin.Weight,sum,na.rm=T)
  write.csv(Total.catch.fin,paste(hndl,"/Shark_fins_TDGDLF.csv",sep=""),row.names=F)
}

#G 4.24 Total shark and scalefish catch by year block for ASL closure for Tim Nicholas-------------------------------------------------------------------------
if(do.Tim_ASL=="YES")
{  
  dummy=subset(Data.daily,FisheryCode%in%c("SGL1","SGL2","WCGL") 
               & FINYEAR%in%c("2006-07","2007-08","2008-09","2009-10","2010-11","2011-12",
                              "2012-13","2013-14","2014-15","2015-16"))
  
  Elas=aggregate(LIVEWT.c~FisheryCode+day+MONTH+FINYEAR+block10,subset(dummy,SPECIES%in%Elasmo.species),sum)
  Scalies=aggregate(LIVEWT.c~FisheryCode+day+MONTH+FINYEAR+block10,subset(dummy,SPECIES%in%Scalefish.species),sum)
  names(Elas)[match("LIVEWT.c",names(Elas))]="Shark_ray_tot_catch_kg"
  names(Scalies)[match("LIVEWT.c",names(Scalies))]="Scalefish_tot_catch_kg"
  Tab=merge(Elas,Scalies,by=c("FisheryCode","day","MONTH","FINYEAR","block10"),all=T)
  rm(dummy)
  write.csv(Tab,paste(hndl,"/ASL/ASL_catch_block.csv",sep=""),row.names=F)
}

#G 4.25 Bull shark catches for Adrian Gleiss-------------------------------------------------------------------------
if(do.Adrian.Gleiss=="YES")
{
  Bull=subset(Data.monthly,SPECIES==18021)
  Bull.daily=subset(Data.daily,SPECIES==18021)
  Bull=subset(Bull,!Same.return%in%unique(Bull.daily$Same.return))
  Bull=subset(Bull,select=c(MONTH,FINYEAR,METHOD,LAT,LONG,SNAME,LIVEWT.c))
  Bull.daily=subset(Bull.daily,select=c(day,MONTH,FINYEAR,METHOD,LAT,LONG,RSCommonName,LIVEWT.c))
  write.csv(Bull,paste(hndl,"/Bulls.csv",sep=""),row.names=F)
  write.csv(Bull,paste(hndl,"/Bulls_daily.csv",sep=""),row.names=F)
}

#G 4.26 Alexandra.Hoschke data extraction-------------------------------------------------------------------------
if(do.Alexandra.Hoschke=="YES")
{
  #reported catch data
  a=subset(Data.monthly,SPECIES==8001,select=c(SPECIES,SNAME,FINYEAR,YEAR.c,MONTH,METHOD,BLOCKX,LAT,LONG,LIVEWT.c))
  #note:no greynurse reported in Daily logbooks
  
  #reported interactions
  TEPS.grey.nurse=subset(TEPS.current,SpeciesCode==8001,select=c(DailySheetNumber,DataEntryName,Status,ScientificName,Comments))
  TEPS.grey.nurse$Status=with(TEPS.grey.nurse,ifelse(Status%in%c("A","a"),"Alive",ifelse(Status%in%c("D","d"),"Dead",NA)))
  names(TEPS.grey.nurse)[match("Status",names(TEPS.grey.nurse))]="Discarding.condition"
  xx=unique(TEPS.grey.nurse$DailySheetNumber)
  # Get=subset(Data.daily.1,TSNo%in%xx,select=c(TSNo,finyear,date,depthMax,method,
  #                                             Block,block10,Lat,LatDeg,LatMin,Long,LongDeg,LongMin))
  # Get=Get[!duplicated(Get$TSNo),]
  # names(Get)[match("LAT",names(Get))]="LatDeg"
  # 
  Get=subset(Data.daily.original,DSNo%in%xx,select=c(DSNo,FINYEAR,date,depthMax,METHOD,
                                                     Block,block10,LAT,LatMin,LongDeg,LongMin))
  names(Get)[match("LAT",names(Get))]="LatDeg"
  Get=Get[!duplicated(Get$DSNo),]
  TEPS.grey.nurse=merge(TEPS.grey.nurse,Get,by.x="DailySheetNumber",by.y="DSNo",all.x=T)
  write.csv(a,paste(hndl,"/Alexandra.Hoschke.csv",sep=""),row.names=F)
  write.csv(TEPS.grey.nurse,paste(hndl,"/Alexandra.Hoschke_TEPS.csv",sep=""),row.names=F)
}

#G 4.27 ASL closures compensation-------------------------------------------------------------------------
if(ASL.compensation=="YES")
{
  yr=c("2012-13","2013-14","2014-15","2015-16","2016-17")
  
  
  #catch
  a=subset(Data.daily,FINYEAR%in%yr & zone%in%c('West','Zone1','Zone2'),select=c(FINYEAR,block10,LIVEWT.c))
  Block_ktch=aggregate(LIVEWT.c~block10+FINYEAR,a,sum)
  
  #effort
  if(Use.Date=="NO")
  {
    b=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    b1=aggregate(Km.Gillnet.Hours.c~ID+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
  }
  if(Use.Date=="YES")
  {
    b=aggregate(Km.Gillnet.Days.c~date+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    b1=aggregate(Km.Gillnet.Hours.c~date+vessel+zone+finyear+block10,data=subset(Effort.daily,finyear%in%yr),max,na.rm=T)
    
  }
  b=aggregate(Km.Gillnet.Days.c~vessel+zone+finyear+block10,data=b,sum,na.rm=T)
  b1=aggregate(Km.Gillnet.Hours.c~vessel+zone+finyear+block10,data=b1,sum,na.rm=T)
  b2=merge(b,b1,by=c("vessel","zone","finyear","block10"),all.x=T)
  b2=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~block10+finyear,data=b2,sum,na.rm=T)
  
  uniqVes <- function(x)  length(unique(x))
  VeslS=aggregate(vessel~ finyear+block10, data = subset(Effort.daily,finyear%in%yr), FUN= uniqVes)
  b2=merge(b2,VeslS,by=c("finyear","block10"))
  
  
  yr_col=paste(yr[1],yr[length(yr)],sep=".to.")
  write.csv(Block_ktch,paste(hndl,"/Block_ktch_",yr_col,".csv",sep=""),row.names=F)
  write.csv(b2,paste(hndl,"/Block_effort_",yr_col,".csv",sep=""),row.names=F)
}

#G 4.28 CITES 2018-------------------------------------------------------------------------
if(do.CITES=="YES")
{
  a=subset(Data.monthly,SPECIES%in%c(10001,18001,18003,18004,26999) & FINYEAR %in%FINYrs)
  a$SPECIES=with(a,ifelse(SPECIES==18001,18003,SPECIES))
  
  Tab=aggregate(LIVEWT.c~FINYEAR+SPECIES,a,sum)
  Spnms=a[!duplicated(a$SPECIES),match(c('SPECIES','SNAME'),names(a))]
  Tab=merge(Tab,Spnms,by='SPECIES')
  Tab=subset(Tab,select=c(FINYEAR,LIVEWT.c,SNAME))
  Tab1=reshape(Tab,v.names = "LIVEWT.c", idvar = "FINYEAR",
               timevar = "SNAME", direction = "wide")
  names(Tab1)[2:ncol(Tab1)]=paste(substr(names(Tab1)[2:ncol(Tab1)],10,50),"(kg)")
  Tab1=Tab1[order(Tab1$FINYEAR),]
  write.csv(Tab1,paste(hndl,"/CITES_2018.csv",sep=""),row.names=F)
  PRCES=subset(PRICES,RSCommonName%in%c("Bronze Whaler","Blue Shark","Guitarfishes","Shortfin Mako"))
  colnames(PRCES)[match("uv1617",colnames(PRCES))]="$/kg"
  write.csv(PRCES,paste(hndl,"/CITES_2018_PRICES.csv",sep=""),row.names=F)
}

#G 4.29 Changes in catch composition of WCDGDLF due to changes in mesh size-------------------------------------------------------------------------
if(do.Nick_mesh.size.WCDGDLF=="YES")
{
  if("plyr" %in% (.packages())) detach(package:plyr)
  
  TAB=with(subset(Data.daily,zone=='West'),table(RSCommonName))
  fn.get=function(d)
  {
    return(WCDGDLF=d %>%
             filter(zone=='West' & RSCommonName%in% names(TAB[TAB>100])) %>%
             dplyr::select(Same.return.SNo,FINYEAR,RSCommonName,LIVEWT.c) %>% 
             group_by(RSCommonName,FINYEAR)  %>% 
             summarize(Annual.ktch_tonnes=sum(LIVEWT.c/1000)) %>%
             spread(RSCommonName, Annual.ktch_tonnes)    %>%
             as.data.frame)
  }
  
  write.csv(fn.get(d=Data.daily),paste(hndl,"/WCDGDLF_catch_compo_by_year.csv",sep=""),row.names=F)
  write.csv(fn.get(d=subset(Data.daily,METHOD=="GN")),paste(hndl,"/WCDGDLF_catch_compo_by_year_gillnet.csv",sep=""),row.names=F)
  write.csv(fn.get(d=subset(Data.daily,METHOD=="LL")),paste(hndl,"/WCDGDLF_catch_compo_by_year_longline.csv",sep=""),row.names=F)
}

#G 4.30 do.ASL.action.2018-------------------------------------------------------------------------
if(do.ASL.action.2018=="YES")
{
  TAB=TEPS.current %>%
    filter(SpeciesCode%in%c(999902)) %>%
    group_by(finyear,month,blockx,DataEntryName) %>%
    summarise(Number = sum(Number)) %>%
    arrange(finyear,month,blockx) %>%
    spread(key = "DataEntryName", value = "Number", fill = 0) %>%
    as.data.frame
  write.csv(TAB,paste(hndl,"/ASL/ASL_action2_2018.csv",sep=""),row.names=F)
  
} 

#G 4.31 Do annual extraction of TEPS by calendar year-------------------------------------------------------------------------
TEP.yr=2018
if(do.annual.TEPS.extraction=="YES")
{
  These.TEPS=read.csv(handl_OneDrive("Analyses/Data_outs/TEPS.current.csv"),stringsAsFactors = F)
  TEPs=These.TEPS%>%filter(fishery%in%c("SGL","WCGL") & year==TEP.yr)%>%
    mutate(LAT=-(as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6)),
           LONG=100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6),
           Bioregion=as.character(ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
                                         ifelse(LONG<115.5 & LAT<=(-27),"WC",
                                                ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
                                                       ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))),
           Bioregion=ifelse(Bioregion=="SC"& LAT>(-34) & LONG <115.91 ,"WC",Bioregion))
  
  #scan comments to check if not reported as a record
  write.csv(TEPs%>%dplyr::select(DailySheetNumber,SessionNumber,SpeciesCode,Status,
                                 Number,DataEntryName,ScientificName,Comments),paste(hndl,"/SCAN.TEPS.csv",sep=""))
  
  
  #export table
  TEPs=TEPs%>%group_by(DataEntryName,Bioregion,block10,Status)%>%
    summarise(N=sum(Number))%>%
    data.frame
  write.csv(TEPs,paste(hndl,"/TEPS_annual_reporting/TEPS.",TEP.yr,".csv",sep=""),row.names=F)
  
}

#G 4.32 Total shark and scalefish catch by year block for ASL closure for Paul Rogers-------------------------------------------------------------------------
if(do.Paul.Rogers_ASL=="YES")  
{
  Hndl.Rogr=paste(hndl,"/ASL/Peter Rogers/",sep='')
  Yrs=c("2014-15","2015-16","2016-17","2017-18","2018-19")  
  
  #Catch by year-block-gear of sharks and scalefish
  dummy=Data.daily%>%
    filter(FINYEAR%in%Yrs & FisheryCode%in%c("SGL1","SGL2","WCGL") &
             METHOD%in%c("GN",'LL'))
  Elas=dummy%>%filter(SPECIES%in%Elasmo.species)%>%
    group_by(FINYEAR,block10,METHOD)%>%
    summarise(Shark_ray_tot_catch_kg=sum(LIVEWT.c))
  
  Scalies=dummy%>%filter(SPECIES%in%Scalefish.species)%>%
    group_by(FINYEAR,block10,METHOD)%>%
    summarise(Scalefish_tot_catch_kg=sum(LIVEWT.c))
  
  Tab=full_join(Elas,Scalies,by=c("FINYEAR","block10","METHOD"))
  
  write.csv(Tab,paste(Hndl.Rogr,"Table_catch_by.year_block_gear.csv",sep=""),row.names=F)
  
  
  #Maps
  South.WA.lat=c(-36,-25); South.WA.long=c(112,130)
  PLATE=c(.01,.9,.075,.9)
  aa= Data.daily.original%>%filter(FINYEAR%in%Yrs) %>%
    mutate(LatDeg=as.numeric(substr(block10,1,2)),
           LatMin=ifelse(is.na(LatMin),10*as.numeric(substr(block10,3,3)),LatMin),
           Lat=-abs(LatDeg+(LatMin/60)),
           LongDeg=ifelse(is.na(LongDeg),100+as.numeric(substr(block10,4,5)),LongDeg), 
           LongMin=ifelse(is.na(LongMin),10*as.numeric(substr(block10,6,6)),LongMin),
           Long=LongDeg+(LongMin/60))%>%
    filter(Lat<=(-26) & Lat>(-36.5)& Long<=(129) & Long >(111.9))%>%
    filter(METHOD%in%c('GN',"LL"))%>%
    mutate(METHOD=ifelse(METHOD=="GN","Gillnet",ifelse(METHOD=="LL","Long line",NA)))
  
  numInt=20
  couleurs=rev(heat.colors(numInt)) 
  tcl.1=.5
  tcl.2=.5
  Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
  Lat.seq=c(-26,-28,-30,-32,-34)
  numberLab=5
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  
  #lodged returns
  TAB2=aa %>% filter(FINYEAR%in%Yrs) %>%
    group_by(FINYEAR,METHOD)%>%
    summarise(Unique_TSNo=n_distinct(TSNo))%>%
    data.frame%>%
    rename(Lodged.returns=Unique_TSNo)
  write.csv(TAB2,paste(Hndl.Rogr,"Table_lodged.returns.by.year.csv",sep=""),row.names = F)
  
  TAB3=aa%>%
    group_by(METHOD,FINYEAR,VESSEL)%>%
    summarise(Trips=n_distinct(TSNo))%>%
    spread(METHOD, Trips)%>%
    replace(is.na(.), "")%>%
    data.frame%>%
    left_join(aa%>%distinct(VESSEL,.keep_all = T)%>%dplyr::select(VESSEL,BoatName),by="VESSEL")%>%
    dplyr::select(FINYEAR,VESSEL,BoatName,Gillnet,Long.line)%>%
    rename(Gillnet.trips=Gillnet,
           Long.line.trips=Long.line)
  write.csv(TAB3,paste(Hndl.Rogr,"Table_trips_by_gear_fisher.csv",sep=""),row.names = F)
  
  
  #effort
  b=aa %>% filter(METHOD=="Gillnet") %>%
    mutate(Km.Gillnet.Hours=HOURS*NETLEN/1000)%>%
    filter(Km.Gillnet.Hours>0)%>%
    group_by(Same.return.SNo,FINYEAR, block10) %>%
    summarize(Km.Gillnet.Hours = max(Km.Gillnet.Hours, na.rm = TRUE))%>%
    group_by(FINYEAR, block10) %>%
    summarize(sum = sum(Km.Gillnet.Hours, na.rm = TRUE))%>%
    mutate(LatDeg=as.numeric(substr(block10,1,2)),
           LatMin=10*as.numeric(substr(block10,3,3)),
           Lat=-abs(LatDeg+(LatMin/60)),
           LongDeg=100+as.numeric(substr(block10,4,5)), 
           LongMin=10*as.numeric(substr(block10,6,6)),
           Long=LongDeg+(LongMin/60))%>%
    data.frame
  b=subset(b,select=c(FINYEAR,sum,Lat,Long))
  BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
  
  for(y in 1:length(Yrs))
  {
    tiff(file=paste(Hndl.Rogr,"Gillnet effort_",Yrs[y],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.75, 0.75, 0.2, 0.1),mgp=c(.1, 0.15, 0))
    
    bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
      arrange(Lat,Long)
    long=sort(unique(bb$Long))
    lat=sort(unique(bb$Lat))      
    Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]	
    
    a=South.WA.long[1]:South.WA.long[2]
    bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
    axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
    axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                         rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
    #nnn=with(TAB2%>%filter(FINYEAR==Yrs[y]),Yrs[y])
    nnn=with(TAB2%>%filter(FINYEAR==Yrs[y]),paste(paste(Yrs[y]," (gillnet returns= ",Lodged.returns[1],
                                                        "; longline returns= ",Lodged.returns[2],")",sep="")))
    mtext(nnn,3,-2,cex=1.25)
    mtext(Lat.exp,side=2,line=-1,las=3,cex=1.25,outer=T)
    mtext(Lon.exp,side=1,line=0,cex=1.25,outer=T)
    mtext("Effort (Km.gn.hours)",3,-0.75,outer=T)
    dev.off()
  }
  
  
  # shark & ray catch
  {
    b=aa %>% filter(METHOD=="Gillnet" & species%in%Elasmo.species) %>%
      group_by(FINYEAR, block10) %>%
      summarize(sum = sum(livewt, na.rm = TRUE))%>%
      mutate(LatDeg=as.numeric(substr(block10,1,2)),
             LatMin=10*as.numeric(substr(block10,3,3)),
             Lat=-abs(LatDeg+(LatMin/60)),
             LongDeg=100+as.numeric(substr(block10,4,5)), 
             LongMin=10*as.numeric(substr(block10,6,6)),
             Long=LongDeg+(LongMin/60))%>%
      data.frame
    b=subset(b,select=c(FINYEAR,sum,Lat,Long))
    BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
    
    for(y in 1:length(Yrs))
    {
      tiff(file=paste(Hndl.Rogr,"Catch_shark_ray_",Yrs[y],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      par(mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.75, 0.75, 0.2, 0.1),mgp=c(.1, 0.15, 0))
      
      bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
        arrange(Lat,Long)
      long=sort(unique(bb$Long))
      lat=sort(unique(bb$Lat))      
      Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      a=South.WA.long[1]:South.WA.long[2]
      bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                           rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
      mtext(Yrs[y],3,-2,cex=1.25)
      mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-1,las=3,cex=1.25,outer=T)
      mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=0,cex=1.25,outer=T)
      mtext("Shark & ray catch (kg)",3,-0.75,outer=T)
      dev.off()
    }
  }
  
  # scalefish catch
  {
    b=aa %>% filter(METHOD=="Gillnet" & species%in%Scalefish.species) %>%
      group_by(FINYEAR, block10) %>%
      summarize(sum = sum(livewt, na.rm = TRUE))%>%
      mutate(LatDeg=as.numeric(substr(block10,1,2)),
             LatMin=10*as.numeric(substr(block10,3,3)),
             Lat=-abs(LatDeg+(LatMin/60)),
             LongDeg=100+as.numeric(substr(block10,4,5)), 
             LongMin=10*as.numeric(substr(block10,6,6)),
             Long=LongDeg+(LongMin/60))%>%
      data.frame
    b=subset(b,select=c(FINYEAR,sum,Lat,Long))
    BREAKS=quantile(b$sum,probs=seq(0,1,1/numInt),na.rm=T)
    
    par(mfrow=c(2,1),mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.5, 0.4, 0.2, 0.1),mgp=c(.1, 0.15, 0))
    for(y in 1:length(Yrs))
    {
      tiff(file=paste(Hndl.Rogr,"Catch_scalefish_",Yrs[y],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      par(mai = c(0.3, 0.4, 0.15, 0.2),oma = c(0.75, 0.75, 0.2, 0.1),mgp=c(.1, 0.15, 0))
      
      bb=subset(b,FINYEAR==Yrs[y],select=-FINYEAR)%>%
        arrange(Lat,Long)
      long=sort(unique(bb$Long))
      lat=sort(unique(bb$Lat))      
      Reshaped=as.matrix(reshape(bb,idvar="Long",timevar="Lat",v.names="sum", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      a=South.WA.long[1]:South.WA.long[2]
      bx=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      plotmap(a,bx,PLATE,"transparent",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.1)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
      if(y==2)color.legend(129.5,South.WA.lat[2],South.WA.long[2],-33,round(BREAKS,0),
                           rect.col=couleurs,gradient="y",col=colLeg,cex=0.85)
      mtext(Yrs[y],3,-2,cex=1.25)
      mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-1,las=3,cex=1.25,outer=T)
      mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=0,cex=1.25,outer=T)
      mtext("Scalefish catch (kg)",3,-0.75,outer=T)
      dev.off()
    }
  }
  
}

#G 4.33 Overlap between ASL closures and effort-------------------------------------------------------------------------
if(do.ASL.closure_effort_overlap)
{
  hndl.ASL.closures.overlap=paste(hndl,'ASL/ASL_overlap.closure.effort',sep='/')
  paz=handl_OneDrive("Data/Mapping/Closures/")
  ASL_Closures=readOGR(paste(paz,"ASL_Closures/ASL_Closures.shp",sep=''),
                       layer="ASL_Closures") 
  ASL_Closures_west.coast=data.frame(
    Point=c(paste("A",1:18,sep="."),
            paste("B",1:12,sep="."),
            paste("C",1:8,sep=".")),
    Lon=c(113.75660,113.89002,114.00565,114.0727,
          114.17,114.24888,114.29735,114.27432,
          114.1833,114.04712,113.88545,113.7455,
          113.6421,113.59458,113.61688,113.58822,
          113.58155,113.64345,
          
          114.97348,114.853,114.73920,114.6717,
          114.6753,114.74808,114.80042,114.72990,
          114.7305,114.80023,114.91950,115.04328,
          
          115.07697,115.00548,114.94362,114.93220,
          114.96693,115.0388,115.13683,115.19965),
    
    
    Lat=c(-28.45407,-28.45165,-28.50912,-28.61142,
          -28.66167,-28.74223,-28.87392,-29.01473,
          -29.13078,-29.19832,-29.20673,-29.1564,
          -29.05877,-28.92743,-28.78988,-28.73745,
          -28.62003,-28.51593,
          
          -29.64328,-29.6233,-29.66638,-29.7584,
          -29.8679,-29.95597,-29.98138,-30.07547,
          -30.1861,-30.27803,-30.32247,-30.30342,
          
          -30.50578,-30.52662,-30.59515,-30.67653,
          -30.75267,-30.80168,-30.81597,-30.79598))
  paz=handl_OneDrive("Data/Mapping/Closures/")
  Closed.ASL=read.csv(paste(paz,'ASL_Closures_Block_Intersection.csv',sep=""))
  source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
  
  if(!Use.Date=="YES") Effrt=aggregate(Km.Gillnet.Days.c~Same.return.SNo+vessel+finyear+month+
                                         LAT+LONG+block10,data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)    
  if(Use.Date=="YES") Effrt=aggregate(Km.Gillnet.Days.c~date+Same.return.SNo+vessel+finyear+month+
                                        LAT+LONG+block10,data=subset(Effort.daily,netlen.c>100 & method=="GN"),max,na.rm=T)    
  
  fn.plt.map=function(Xlim,Ylim,col.ASL.closure,Efrt)
  {
    plot(1,xlim=c(Xlim[1]*0.9995,Xlim[2]*0.99),ylim=Ylim,xlab="",ylab="",axes=F,main="")
    
    #ASL closure
    plot(ASL_Closures,add=T,col=col.ASL.closure)
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("A",1:18,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("B",1:12,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("C",1:8,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    
    #Effort
    with(Efrt,points(LONG,LAT,pch=21,cex=.5,col="transparent",bg=rgb(.8,.1,.1,alpha=0.3)))
    
    #Coast
    polygon(WAcoast$Longitude,WAcoast$Latitude, col="grey80")
    box()
  }
  Inset=function(Xlim,Ylim,col.ASL.closure)
  {
    plot(1,xlim=c(Xlim[1]*0.9995,Xlim[2]*0.99),ylim=Ylim,xlab="",ylab="",axes=F,main="")
    
    #ASL closure
    plot(ASL_Closures,add=T,col=col.ASL.closure)
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("A",1:18,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("B",1:12,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    with(ASL_Closures_west.coast%>%filter(Point%in%paste("C",1:8,sep=".")),
         {
           polygon(Lon,Lat, col=col.ASL.closure)
         })
    #Coast
    polygon(WAcoast$Longitude,WAcoast$Latitude, col="grey70")
    mtext("ASL exclusion zones",3,-2,cex=1.5)
    
    text(114.58,-28.77,"Geraldton",pos=4)
    text(123,-33.6,"Esperance",pos=3,srt=25)
  }
  ASL.yrs=paste(2016:as.numeric(substr(Current.yr,1,4)),
                substr(2017:(1+as.numeric(substr(Current.yr,1,4))),3,4),sep="-")
  tiff(file=paste(hndl.ASL.closures.overlap,"ASL.closures_daily.effort_overlap.tiff",sep="/"),
       width = 2400, height = 2200,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(ASL.yrs)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  Inset(Xlim=c(113,129.5),Ylim=c(-36,-26),"steelblue") #ASL closures
  for (i in 1:length(ASL.yrs))
  {
    fn.plt.map(Xlim=c(113,129.5),
               Ylim=c(-36,-26),
               col.ASL.closure='steelblue',
               Efrt=subset(Effrt,finyear==ASL.yrs[i]))
    axis(side = 1, at =Long.seq, labels =F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
    mtext(ASL.yrs[i],side=3,line=-1.75,cex=1.5)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
  }  
  mtext(Lat.exp,side=2,line=.1,las=3,cex=1.25,outer=T)
  mtext(Lon.exp,side=1,line=.8,cex=1.25,outer=T)
  dev.off()
  
  
  #Percentage change
  ASL_closure_block10=data.frame(Block10=c(unlist(ASL_exclusions_block10_JASDGDLF),
                                           unlist(ASL_exclusions_block10_WCDGDLF)))
  a=Effrt%>%
    filter(finyear%in%ASL.yrs & block10%in%ASL_closure_block10$Block10)%>%
    group_by(finyear,block10)%>%
    summarise(Total=sum(Km.Gillnet.Days.c))%>%
    arrange(block10,finyear)%>%
    spread(finyear,Total,fill=0)%>%
    mutate('2016-17 to 2017-18'=ifelse(`2016-17`>0,round(100*`2017-18`/`2016-17`,2),0),
           '2017-18 to 2018-19'=ifelse(`2017-18`>0,round(100*`2018-19`/`2017-18`,2),0),
           '2018-19 to 2019-20'=ifelse(`2018-19`>0,round(100*`2019-20`/`2018-19`,2),0))%>%
    data.frame%>%
    dplyr::select(-X2016.17, -X2017.18, -X2018.19, -X2019.20)
  
  write.csv(a,paste(hndl.ASL.closures.overlap,"ASL.closures_percent.change.csv",sep="/"))
}

#G 4.34 Explore school shark targeting-------------------------------------------------------------------------
if(Check.school.shark.targeting=="YES")
{
  library(ozmaps)
  library(sf)
  library(ggrepel)
  library(ggpubr)
  library(gridExtra)
  
  #Targeted trip
  d=Data.daily.original%>%filter(TSNo%in%c('TDGLF8010741','TDGLF6003155'))%>%
    group_by(DSNo,TSNo,SNo,RSCommonName,date,VESSEL,port,BoatName,MastersName,Lat,Long,block10)%>%
    summarise(Tonnes=round(sum(livewt)/1000,2))%>%
    data.frame%>%
    mutate(Legend=paste(VESSEL,TSNo,port,MastersName,sep=', '),
           Latitude=-(as.numeric(substr(block10,1,2))+10*as.numeric(substr(block10,3,3))/60),
           Longitude=100+as.numeric(substr(block10,4,5))+10*as.numeric(substr(block10,6,6))/60)
  
  p1=d%>%
    ggplot(aes(x=SNo,y=Tonnes, fill=RSCommonName))+
    geom_bar(stat="identity")+ 
    facet_wrap(vars(Legend),scales="free_y")+
    xlab("Shot number")+
    theme(legend.position = 'top',
          legend.title=element_blank(),
          legend.text=element_text(size=rel(1.3)),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          strip.text = element_text(size = 14))
  
  #map
  oz_states <- ozmaps::ozmap_states
  d1=d%>%distinct(block10,.keep_all = T)
  p2=ggplot(oz_states) + 
    geom_sf() + 
    coord_sf(xlim = c(112,129),ylim = c(-36,-26))+
    geom_point(data=d1,aes(Longitude,Latitude),size=2)+
    geom_label_repel(data=d1,aes(Longitude,Latitude,label=block10))+
    theme(legend.position = 'top',
          legend.title=element_blank(),
          legend.text=element_text(size=rel(1.3)),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12))
  
  
  infographic=grid.arrange(p1,p2, nrow = 2,ncol=1,heights=c(3,3))
  annotate_figure(infographic,
                  bottom = text_grob("",size = 20),
                  left = text_grob("",rot = 90,size = 20))
  
  ggsave(paste(hndl,'School.shark.targeting.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
}

#G 4.35 Zone 1 catch east and west Black Point-------------------------------------------------------------------------
if(do.Peter.Rogers)
{
  Suite=read.csv(handl_OneDrive("Analyses/Data_outs/suite.csv"))$Suite
  
  Black.Poin.lat=-34.4174
  Black.Poin.lon=115.5432
  
  Rogers=Data.daily%>%filter(zone=='Zone1' & fishery=='JASDGDL')%>%
    rename(species=SPECIES,
           sname1=SNAME,
           method=METHOD,
           finyear=FINYEAR)%>%
    mutate(latitude=-abs(as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6)),
           longitude=100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6),
           zone1_Black.point=ifelse(longitude<115.5,'West.BP','East.BP'),
           zone1_Black.point=ifelse(latitude>(-34),'West.BP',zone1_Black.point),
           Suite=ifelse(species%in%Suite,'demersal suite',NA))
  
  SPP=Rogers%>%distinct(species,.keep_all=T)%>%dplyr::select(species,sname1)%>%arrange(species)
  Rogers=Rogers%>%
    dplyr::select(-sname1)
  
  Rogers=Rogers%>%
    left_join(SPP,by='species')
  
  Rogers%>%
    distinct(longitude,latitude,zone1_Black.point)%>%
    ggplot(aes(longitude,latitude,color=zone1_Black.point))+
    geom_point(size=2)+
    theme(legend.title = element_blank())
  
  
  A=Rogers%>%
    mutate(sname1=tolower(sname1))%>%
    filter(!is.na(sname1))%>%
    filter(!finyear=='2020-21')%>%
    group_by(Suite,sname1,finyear,zone1_Black.point,method)%>%
    summarise(Tons=sum(LIVEWT)/1000)
  
  
  write.csv(A,paste(hndl,"/Peter Rogers/PeterRogers_2022_data.csv",sep=""),row.names=F)
  
  #Check snapper and dhufish hotspots
  a_ktch=Rogers%>%
    filter(method=="GN")%>%
    mutate(Target=case_when(species==353001~'Snapper',
                            species==320000~'Dhufish',
                            TRUE~'other'))%>%
    group_by(Target,MONTH,latitude,longitude)%>%
    summarise(Total.catch.kg=sum(LIVEWT.c,na.rm=T))%>%
    ungroup()
  
  
  a_ktch%>%
    filter(Target=='Dhufish')%>%
    ggplot(aes(longitude,latitude,color=Total.catch.kg))+
    geom_point(aes(size=Total.catch.kg))+
    facet_wrap(~MONTH)+
    ggtitle('Dhufish')
  ggsave(paste(hndl,'/Peter Rogers/Hotspots_Dhufish.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  a_ktch%>%
    filter(Target=='Snapper')%>%
    ggplot(aes(longitude,latitude,color=Total.catch.kg))+
    geom_point(aes(size=Total.catch.kg))+
    facet_wrap(~MONTH)+
    ggtitle('Snapper')
  ggsave(paste(hndl,'/Peter Rogers/Hotspots_Snapper.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  #Check if removing targeted shots makes a difference
  #no major effect
  do.this=FALSE
  if(do.this)
  {
    a=Rogers%>%
      filter(method=="GN")%>%
      mutate(Targeting=case_when(species==353001~'Snapper',
                                 species==320000~'Dhufish',
                                 TRUE~'other'))%>%
      group_by(Same.return.SNo,finyear,MONTH,VESSEL,Targeting)%>%
      summarise(Targeting.total=sum(LIVEWT.c,na.rm=T))%>%
      ungroup()%>%
      spread(Targeting,Targeting.total,fill=0)%>%
      mutate(Total=Dhufish +other + Snapper,
             Snapper.prop=Snapper/Total,
             Dhufish.prop=Dhufish/Total,
             Snapper.prop.bin=cut(Snapper.prop,breaks=10),
             Dhufish.prop.bin=cut(Dhufish.prop,breaks=10))
    
    p1=a%>%
      mutate(Group=ifelse(Dhufish.prop>0.25,'>0.25','<0.25'))%>%
      ggplot(aes(Dhufish.prop,fill=Group))+
      geom_histogram()+xlab("Proportion of catch")+
      ggtitle('WA Dhufish')
    
    p2=a%>%
      mutate(Group=ifelse(Snapper.prop>0.25,'>0.25','<0.25'))%>%
      ggplot(aes(Snapper.prop,fill=Group))+
      geom_histogram()+xlab("Proportion of catch")+
      ggtitle('Pink Snapper')
    
    
    Tot=a%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Snapper)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='All')
    Tot.minus.0.25.plus.Snapper.prop=a%>%
      filter(Snapper.prop<=0.25)%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Snapper)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='<0.25')
    
    Tot.plus.0.25.plus.Snapper.prop=a%>%
      filter(Snapper.prop>0.25)%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Snapper)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='>0.25')
    
    p4=rbind(Tot,Tot.minus.0.25.plus.Snapper.prop,Tot.plus.0.25.plus.Snapper.prop)%>%
      ggplot(aes(year,Ktch,color=Group))+
      geom_line(size=1.5)+ylab("Total catch (tonnes)")+
      scale_color_manual(values=c("#F8766D", "#00BFC4", "#00BA38"))
    
    
    Tot=a%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Dhufish)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='All')
    Tot.minus.0.25.plus.Dhufish.prop=a%>%
      filter(Dhufish.prop<=0.25)%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Dhufish)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='<0.25')
    Tot.plus.0.25.plus.Dhufish.prop=a%>%
      filter(Dhufish.prop>0.25)%>%
      group_by(finyear)%>%
      summarise(Ktch=sum(Dhufish)/1000)%>%
      ungroup()%>%
      mutate(year=as.numeric(substr(finyear,1,4)),
             Group='>0.25')
    
    
    p3=rbind(Tot,Tot.minus.0.25.plus.Dhufish.prop,Tot.plus.0.25.plus.Dhufish.prop)%>%
      ggplot(aes(year,Ktch,color=Group))+
      geom_line(size=1.5)+ylab("Total catch (tonnes)")+
      scale_color_manual(values=c("#F8766D", "#00BFC4", "#00BA38"))
    
    ggarrange(plotlist=list(p1,p2,p3,p4),nrow=2,ncol=2)
    ggsave(paste(hndl,'/Peter Rogers/Species composition and catch.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
    
  }
  
  #Effort Zn 1 east and west Black point
  a=Effort.daily%>%
    filter(FisheryCode=='JASDGDL' & zone=='Zone1' & !finyear=='2020-21' & method=='GN')%>%
    mutate(latitude1=-abs(as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6)),
           longitude1=100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6),
           zone1_Black.point=ifelse(longitude1<115.5,'West.BP','East.BP'),
           zone1_Black.point=ifelse(latitude1>(-34),'West.BP',zone1_Black.point))
  
  # a%>%distinct(latitude1,longitude1,zone1_Black.point)%>%
  #   ggplot(aes(longitude1,latitude1))+
  #   geom_point(aes(color=zone1_Black.point))
  
  #temporal effort
  a%>%
    group_by(date,vessel,zone,finyear,zone1_Black.point)%>%
    summarise(Km.Gillnet.Days.c=max(Km.Gillnet.Days.c,na.rm=T))%>%
    ungroup()%>%
    group_by(finyear,zone1_Black.point)%>%
    summarise(Km.Gillnet.Days.c=sum(Km.Gillnet.Days.c,na.rm=T))%>%
    ungroup()%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    ggplot(aes(year,Km.Gillnet.Days.c,color=zone1_Black.point))+
    geom_point(size=2)+
    geom_line(linetype = "dashed")+expand_limits(y = 0)+
    xlab('Financial year')+
    theme(legend.position = 'top',
          legend.title = element_blank(),
          text = element_text(size = 18))+
    scale_x_continuous(breaks = seq(min(a$year),max(a$year),2))
  ggsave(paste(hndl,'/Peter Rogers/Effort_zone1.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  
  #spatial effort group_by(month,latitude1,longitude1)
  b=a%>%
    group_by(date,vessel,zone,finyear,zone1_Black.point,month,latitude1,longitude1)%>%
    summarise(Km.Gillnet.Days.c=max(Km.Gillnet.Days.c,na.rm=T))%>%
    ungroup()%>%
    group_by(zone1_Black.point,month,latitude1,longitude1)%>%
    summarise(Km.Gillnet.Days.c=sum(Km.Gillnet.Days.c,na.rm=T))%>%
    ungroup()
  b%>%
    ggplot(aes(longitude1,latitude1,color=Km.Gillnet.Days.c))+
    geom_point(aes(size=Km.Gillnet.Days.c))+
    facet_wrap(~month)+
    ggtitle('Effort (Km.Gillnet.Days)')
  ggsave(paste(hndl,'/Peter Rogers/Effort_zone1_spatial.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  
  #nominal cpue
  a_ktch%>%
    filter(Target=='Dhufish')%>% 
    left_join(b,by=c('MONTH'='month', 'latitude'='latitude1', 'longitude'='longitude1'))%>%
    mutate(cpue=Total.catch.kg/Km.Gillnet.Days.c)%>%
    ggplot(aes(longitude,latitude,color=cpue))+
    geom_point(aes(size=cpue))+
    facet_wrap(~MONTH)+
    ggtitle('CPUE (kg/Km.Gillnet.Days)')
  ggsave(paste(hndl,'/Peter Rogers/Hotspots_Dhufish_cpue.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  a_ktch%>%
    filter(Target=='Snapper')%>% 
    left_join(b,by=c('MONTH'='month', 'latitude'='latitude1', 'longitude'='longitude1'))%>%
    mutate(cpue=Total.catch.kg/Km.Gillnet.Days.c)%>%
    ggplot(aes(longitude,latitude,color=cpue))+
    geom_point(aes(size=cpue))+
    facet_wrap(~MONTH)+
    ggtitle('CPUE (kg/Km.Gillnet.Days)')
  ggsave(paste(hndl,'/Peter Rogers/Hotspots_Snapper_cpue.tiff',sep='/'), width = 12, height = 8,compression = "lzw")
  
  
  
}

#G 4.36 DBCA Marine Planning 2022-------------------------------------------------------------------------
if(do.DBCA_2022)
{
  Park.boundaries="Data/Mapping/park boundaries/"
  DBCA_SouthCoast_hexagons=readOGR(handl_OneDrive(paste0(Park.boundaries,"DBCA_SCMP_hexagons/south-coast_shapefile.shp")), layer="south-coast_shapefile") 
  DBCA_SouthCoast_sanctuary=readOGR(handl_OneDrive(paste0(Park.boundaries,"DBCA_SCMP/pmcr-scmp-sanctuary-hwm-proposed-sco-20240111.shp")), layer="pmcr-scmp-sanctuary-hwm-proposed-sco-20240111") 
  Closures="Data/Mapping/Closures/"
  ASL_closures=readOGR(handl_OneDrive(paste0(Closures,"ASL_closures/ASL_Closures.shp")), layer="ASL_Closures") 
  
  Boat.skippers=read.csv("C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Data/Catch and Effort/List_boats_skippers.csv")
  Boat.skippers=Boat.skippers%>%distinct(VESSEL,.keep_all = T)
    
    
  library(ggspatial) 
  library(sf)
  library("rnaturalearth")
  library("rnaturalearthdata")
  library(ozmaps)
  dis.yrs=c('2019-20','2020-21','2021-22')
  
  DBCA=Data.daily%>%filter(zone=='Zone2' & fishery=='JASDGDL')%>%
    rename(species=SPECIES,
           sname1=SNAME,
           method=METHOD,
           finyear=FINYEAR)%>%
    mutate(latitude=as.numeric(LatFC),
           longitude=as.numeric(LongFC),
           latitude1=-abs(as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6)),
           longitude1=100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6))%>%
    filter(longitude1>119)%>%
    filter(finyear%in%dis.yrs)
  
  DBCA.effort.GN=Effort.daily%>%filter(zone=='Zone2' & method=='GN')%>%
    mutate(latitude=as.numeric(LAT),
           longitude=as.numeric(LONG),
           latitude1=-abs(as.numeric(substr(block10,1,2))+(as.numeric(substr(block10,3,3))/6)),
           longitude1=100+as.numeric(substr(block10,4,5))+(as.numeric(substr(block10,6,6))/6),
           latitude=ifelse(is.na(latitude),latitude1,latitude),
           longitude=ifelse(is.na(longitude),longitude1,longitude))%>%
    filter(longitude>119)%>%
    filter(finyear%in%dis.yrs)%>%
    dplyr::select(Same.return.SNo,latitude,longitude,Km.Gillnet.Days.c,Km.Gillnet.Hours.c,
                  finyear,vessel)%>%
  distinct(Same.return.SNo,.keep_all = T)
  
  oz_states <- ozmaps::ozmap_states
  Chunks=list(one=c(119,120,-33.5,-34.5),
              two=c(120,122,-33.5,-34.5),
              three=c(122,124,-33.5,-34.5),
              four=c(124,126,-32.5,-34),
              five=c(126,128,-32,-32.5),
              six=c(128,129,-31,-32.25))

  #Overlap of catch with hexagons
  fn.plt.DBCA=function(ktch,XLIM,YLIM,ptsize=5,ParkShapeFile)
  {
    p=ggplot(data = oz_states) +
      geom_sf(color = "black", fill = "sienna4") +
      xlab("Longitude") + ylab("Latitude") +
      geom_sf(data=st_as_sf(ParkShapeFile), color="grey", inherit.aes = FALSE)+
      geom_point(data=ktch,aes(x,y,color=Total.catch.kg),size=ptsize)+
      ylab("Latitude")+xlab('Longitude')+
      labs(title = paste('Total catch', paste(range(DBCA$finyear),collapse = ' to ')),
           caption = 'note: honeycomb overlaid')+
      theme_bw()+
      theme(legend.position = "top",
            legend.justification="right",
            legend.text = element_text(size=18),
            legend.title = element_text(size=20),
            axis.title = element_text(size=20),
            axis.text = element_text(size=16),
            plot.title = element_text(hjust = 0.5,size = 20))+
      scale_colour_gradient(low = "yellow2", high = "red2", na.value = NA)+
      guides(color=guide_legend(title="Tonnes"))+
      coord_sf(xlim = XLIM, ylim = YLIM)
    return(p)
  }
  fn.plt.DBCA(ktch=DBCA%>%
                rename(y=latitude1,
                       x=longitude1)%>%
                group_by(x,y)%>%
                summarise(Total.catch.kg=sum(LIVEWT.c,na.rm=T)/1000)%>%
                ungroup(),
              XLIM=c(119, 129),
              YLIM=c(-30, -35),
              ptsize=3,
              ParkShapeFile=DBCA_SouthCoast_hexagons)
  ggsave(paste(hndl,'DBCA/2022/Spatial catch.tiff',sep='/'), 
         width = 20, height = 10,compression = "lzw",dpi=300)
  
  for(k in 1:length(Chunks))
  {
    fn.plt.DBCA(ktch=DBCA%>%
                  rename(y=latitude1,
                         x=longitude1)%>%
                  group_by(x,y)%>%
                  summarise(Total.catch.kg=sum(LIVEWT.c,na.rm=T)/1000)%>%
                  ungroup(),
                XLIM=c(Chunks[[k]][1], Chunks[[k]][2]),
                YLIM=c(Chunks[[k]][3], Chunks[[k]][4]),
                ParkShapeFile=DBCA_SouthCoast_hexagons)
    ggsave(paste(hndl,'/DBCA/2022/',paste('Spatial catch_chunk_',
                                          paste(Chunks[[k]][1],Chunks[[k]][2],sep='-'),'.tiff',sep=''),sep='/'), 
           width = 20, height = 10,compression = "lzw",dpi=300)
  }
  
  
  #Overlap with sanctuary zones  
  fn.plt.DBCA2=function(Input,XLIM,YLIM,ParkShapeFile,ASLShapeFile,TITL,Leyen)
  {
    p=ggplot(data = oz_states) +
      geom_sf(color = "black", fill = "sienna4") +
      xlab("Longitude") + ylab("Latitude") +
      geom_sf(data=st_as_sf(ASLShapeFile), fill="red", alpha=0.1,inherit.aes = FALSE)+
      geom_sf(data=st_as_sf(ParkShapeFile), fill="green", inherit.aes = FALSE)+
      ylab("Latitude")+xlab('Longitude')+
      labs(title = TITL,caption = 'green: Sanctuary zone\n pink: ASL closure')+
      theme_bw()+
      theme(legend.position = "top",
            legend.justification="right",
            legend.text = element_text(size=18),
            legend.title = element_text(size=20),
            axis.title = element_text(size=20),
            axis.text = element_text(size=16),
            plot.title = element_text(hjust = 0.5,size = 20))+
      scale_colour_gradient(low = "yellow2", high = "red2", na.value = NA)+
      guides(color=guide_legend(title="Tonnes"))+
      coord_sf(xlim = XLIM, ylim = YLIM)+
      geom_point(data=Input,aes(x,y,size=z),alpha=.2)+labs(size=Leyen)
    #geom_contour_filled(data=Input,aes(x=x,y=y, z = z),bins=20)+
    #scale_fill_brewer(palette = "Spectral")
    
    return(p)
  }
  
  dumi='Spatial catch overlap with sanctuary zone'
   #Catch whole area
  fn.plt.DBCA2(Input=DBCA%>%
                 mutate(y=round(latitude,1),
                        x=round(longitude,1))%>%
                 group_by(x,y)%>%
                 summarise(z=sum(LIVEWT.c,na.rm=T)/1000)%>%
                 ungroup(),
               XLIM=c(119, 129),
               YLIM=c(-30, -35),
               ParkShapeFile=DBCA_SouthCoast_sanctuary,
               ASLShapeFile=ASL_closures,
               TITL=paste('Total catch', paste(range(DBCA$finyear),collapse = ' to ')),
               Leyen='Tonnes')
  ggsave(paste0(hndl,'/DBCA/2022/',paste0(dumi,'.tiff')), 
         width = 20, height = 10,compression = "lzw",dpi=300)
  
  #Catch by longitude chunk
  for(k in 1:length(Chunks))
  {
    print(fn.plt.DBCA2(Input=DBCA%>%
                   mutate(y=round(latitude,1),
                          x=round(longitude,1))%>%
                   group_by(x,y)%>%
                   summarise(z=sum(LIVEWT.c,na.rm=T)/1000)%>%
                   ungroup(),
                 XLIM=c(Chunks[[k]][1], Chunks[[k]][2]),
                 YLIM=c(Chunks[[k]][3], Chunks[[k]][4]),
                 ParkShapeFile=DBCA_SouthCoast_sanctuary,
                 ASLShapeFile=ASL_closures,
                 TITL=paste0('Total catch between ',paste(Chunks[[k]][1],Chunks[[k]][2],sep=' & ')),
                 Leyen='Tonnes'))
    ggsave(paste0(hndl,'/DBCA/2022/',paste0(dumi,' ',paste(Chunks[[k]][1],Chunks[[k]][2],sep='_'),'.tiff')), 
           width = 20, height = 10,compression = "lzw",dpi=300)
  }
  
  #Catch by vessel
  dumi='Spatial catch overlap with sanctuary zone vessel'
  Skipper.list=sort(unique(DBCA$VESSEL))
  for(k in 1:length(Skipper.list))
  {
    Input=DBCA%>%
      mutate(y=round(latitude,1),
             x=round(longitude,1))%>%
      filter(VESSEL==Skipper.list[k])
    Skipper=Boat.skippers%>%filter(VESSEL==Skipper.list[k])
    print(fn.plt.DBCA2(Input=Input%>%
                         group_by(x,y)%>%
                         summarise(z=sum(LIVEWT.c,na.rm=T)/1000)%>%
                         ungroup(),
                       XLIM=range(Input$x),
                       YLIM=-range(abs(Input$y)),
                       ParkShapeFile=DBCA_SouthCoast_sanctuary,
                       ASLShapeFile=ASL_closures,
                       TITL=paste0('Total catch. Vessel: ',Skipper$BoatName,', Skipper: ',Skipper$MastersName),
                       Leyen='Tonnes'))
    ggsave(paste0(hndl,'/DBCA/2022/',paste0(dumi,' ',Skipper.list[[k]],'.tiff')), 
           width = 20, height = 10,compression = "lzw",dpi=300)
  }
  
  #Effort by vessel
  Skipper.list=sort(unique(DBCA.effort.GN$vessel))
  dumi='Spatial Gillnet effort overlap with sanctuary zone vessel'
  for(k in 1:length(Skipper.list))
  {
    Input=DBCA.effort.GN%>%
      mutate(y=round(latitude,1),
             x=round(longitude,1))%>%
      filter(vessel==Skipper.list[k])
    Skipper=Boat.skippers%>%filter(VESSEL==Skipper.list[k])
    print(fn.plt.DBCA2(Input=Input%>%
                         group_by(x,y)%>%
                         summarise(z=sum(Km.Gillnet.Hours.c,na.rm=T)/1000)%>%
                         ungroup(),
                       XLIM=range(Input$x),
                       YLIM=-range(abs(Input$y)),
                       ParkShapeFile=DBCA_SouthCoast_sanctuary,
                       ASLShapeFile=ASL_closures,
                       TITL=paste0('Total effort. Vessel: ',Skipper$BoatName,', Skipper: ',Skipper$MastersName),
                       Leyen='Km.Gillnet.Hours.c'))
    ggsave(paste0(hndl,'/DBCA/2022/',paste0(dumi,' ',Skipper.list[[k]],'.tiff')), 
           width = 20, height = 10,compression = "lzw",dpi=300)
  }
  
}

#G 4.37 Data extraction for Matt Koopman 2023 Fishery Review-------------------------------------------------------------------------
if(get.Mats.data.2023)
{
  write.csv(Mat.monthly%>%dplyr::select(-PORT),paste(hndl,'/Matt Koopman/monthly catch and effort.csv',sep=''),row.names = F)
  write.csv(Mat.daily%>%dplyr::select(-PORT),paste(hndl,'/Matt Koopman/daily catch and effort.csv',sep=''),row.names = F)
  
}

#G 4.38 Ngari Capes Compensation 2023 - Tim Nicholas-------------------------------------------------------------------------
if(do.Ngari.2023)
{
  #Ngari Capes Marine Park extends from Busselton (33.6516° S, 115.3473° E) to Augusta (34.3110° S, 115.1573° E)
  
  LatLower=-35.6
  LatUpper=-33.2
  LongRight=115.35
  LongLeft=114.6
  
  hndl.out=paste(hndl,'Compensations/2023_Ngari Capes',sep='/')
  
  a=Data.daily%>%
    distinct(Same.return.SNo,Block,block10,BLOCKX,LAT,LatMin,LONG,LongMin)%>%
    mutate(Lon=LONG+(as.numeric(LongMin)/60),
           Lat=-(abs(LAT)+(as.numeric(LatMin)/60)))%>%
    filter(Lon<LongRight & Lat<=(LatUpper))%>%
    filter(Lon>=LongLeft & Lat>(LatLower))
  a%>%
    ggplot(aes())+
    geom_point(aes(Lon,Lat))+
    geom_polygon(data=WAcoast,aes(Longitude,Latitude),fill='brown4')+
    ylim(LatLower,LatUpper)+xlim(LongLeft,LongRight)
  ggsave(paste(hndl.out,'Map of fishing shots.tiff',sep='/'),width = 6, height = 8,compression = "lzw")

  
  b=Effort.daily%>%
    filter(method=='GN')%>%
    filter(Same.return.SNo%in%unique(a$Same.return.SNo))%>%
    mutate(Hook.hours=hooks*hours.c)%>%
    group_by(Same.return.SNo)%>%
    summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c,na.rm=T))
  
  a=Data.daily%>%
    mutate(Lon=LONG+(as.numeric(LongMin)/60),
           Lat=-(abs(LAT)+(as.numeric(LatMin)/60)))%>%
    filter(Lon<LongRight & Lat<=(LatUpper))%>%
    filter(Lon>=LongLeft & Lat>(LatLower))%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(Same.return.SNo,VESSEL,finyear,METHOD,FisheryCode,Lon,Lat,SNAME,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    left_join(b,by='Same.return.SNo')

  Tot.ktch=a%>%
    group_by(finyear)%>%
    summarise(Tons=sum(LIVEWT.c,na.rm=T)/1000)
  
  Tot.ktch_vessel=a%>%
                  group_by(finyear,VESSEL)%>%
                  summarise(Tons=sum(LIVEWT.c,na.rm=T)/1000)
  
  Tot.ktch_method=a%>%
    group_by(finyear,METHOD)%>%
    summarise(Tons=sum(LIVEWT.c,na.rm=T)/1000)
  
  Tot.ktch_vessel_method=a%>%
    group_by(finyear,VESSEL,METHOD)%>%
    summarise(Tons=sum(LIVEWT.c,na.rm=T)/1000)
  
  Tot.Km_gn_hour_vessel_method=a%>%
    group_by(finyear,VESSEL)%>%
    summarise(Km.Gillnet.Hours.c=sum(Km.Gillnet.Hours.c,na.rm=T))%>%
    ungroup()%>%
    mutate(Km.Gillnet.Hours.c=Km.Gillnet.Hours.c/max(Km.Gillnet.Hours.c))%>%
    filter(!is.infinite(Km.Gillnet.Hours.c))
  

    Tot.ktch_vessel_method%>%
    ggplot(aes(finyear,Tons,color=METHOD))+
    geom_point()+
    facet_wrap(~VESSEL,scales='free_y')+
    ylab("Total reported landings (tonnes)")+
    geom_vline(xintercept=2019,color="orange",size=3,alpha=.4)+
    geom_line(data=Tot.Km_gn_hour_vessel_method, 
              aes(x=finyear,y=Km.Gillnet.Hours.c*1e2),size=1.25,color='black',alpha=0.25)
  
    ggsave(paste(hndl.out,'Landings and GN effort.tiff',sep='/'),width = 6, height = 6,compression = "lzw")
    

  
}

#G 4.39 Roger Kirkwood SARDI. effort for ASL ERA-------------------------------------------------------------------------
if(do.SARDI.2024)
{
  dis.yrs=c(paste(2012:2021,substr(2013:2022,3,4),sep='-'))
  a=Effort.daily%>%
    distinct(Same.return.SNo,method,.keep_all = T)%>%
    mutate(hook.hours=ifelse(method=='LL',hooks*hours.c,NA))%>%
    dplyr::select(vessel,block10,finyear,method,hook.hours,Km.Gillnet.Hours.c)%>%
    filter(method%in%c('GN',"LL") & finyear%in%dis.yrs)
  
  N.vessls.per.blk=a%>%
    distinct(vessel,block10,finyear,method)%>%  
    group_by(block10,finyear,method)%>%
    tally()%>%
    rename(n.ves.per.blk=n)
  
  a1=a%>%
    group_by(block10,finyear,method)%>%
    summarise(hook.hours=sum(hook.hours,na.rm=T), 
              Km.Gillnet.Hours.c=sum(Km.Gillnet.Hours.c,na.rm=T))
  
  a1=a1%>%
    left_join(N.vessls.per.blk,by=c('block10','finyear','method'))%>%
    mutate(hook.hours=ifelse(n.ves.per.blk<3 & method=='LL',NA,hook.hours),
           Km.Gillnet.Hours.c=ifelse(n.ves.per.blk<3 & method=='GN',NA,Km.Gillnet.Hours.c))
  
  a1=a1%>%
    mutate(fishery.Name='Southern Demersal Gillnet and Demersal Longline Managed Fishery',
           n.ves.per.blk=ifelse(n.ves.per.blk<3,"<3",">=3"))%>%
    relocate(fishery.Name,block10,finyear,n.ves.per.blk)%>%
    rename('Fishery Name'=fishery.Name,
           '10x10NM Block'=block10,
           'Financial Year'=finyear,
           'Vessel Count'=n.ves.per.blk,
           Km.Gillnet.Hours=Km.Gillnet.Hours.c)
  
  write.csv(a1,paste(hndl,'/SARDI_2024_ASL_ERA.csv',sep=''),row.names = F)
  
}