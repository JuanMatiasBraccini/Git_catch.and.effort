#Script used for testing the effect of using DATE or ID in the calculation of effort

Use.Date="YES"
#Use.Date="NO"

#km gn days 
if(Use.Date=="NO")
{
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~ID+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~ID+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
}

if(Use.Date=="YES")
{
  Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~date+vessel+zone+finyear,data=Effort.daily,max,na.rm=T) 
  Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~date+vessel+zone+finyear,data=Effort.daily,max,na.rm=T)
}

Attach.Effort.daily.c=aggregate(Km.Gillnet.Days.c~zone+finyear,data=Attach.Effort.daily.c,sum,na.rm=T)
Attach.Effort.daily=aggregate(Km.Gillnet.Days.inv~zone+finyear,data=Attach.Effort.daily,sum,na.rm=T)


#km gn hours
  Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~finyear+vessel+ID+zone,data=Effort.daily,max,na.rm=T)
  Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~finyear+vessel+ID+zone,data=Effort.daily,max,na.rm=T)

# if(Use.Date=="YES")
# {
#   Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
#   Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~finyear+vessel+date+zone,data=Effort.daily,max,na.rm=T)
# }

Attach.Effort.daily.hrs.c=aggregate(Km.Gillnet.Hours.c~zone+finyear,data=Attach.Effort.daily.hrs.c,sum,na.rm=T)
Attach.Effort.daily.hrs=aggregate(Km.Gillnet.Hours.inv~zone+finyear,data=Attach.Effort.daily.hrs,sum,na.rm=T)


#merge into single file
Attach.Effort.daily.c=merge(Attach.Effort.daily.c,Attach.Effort.daily.hrs.c,by=c("zone","finyear"),all=T)
Attach.Effort.daily=merge(Attach.Effort.daily,Attach.Effort.daily.hrs,by=c("zone","finyear"),all=T)
Attach.Effort.daily=merge(Attach.Effort.daily.c,Attach.Effort.daily,by=c("zone","finyear"),all.x=T)


#2.3.2. Monthly 

#Rory's Zone 2 multiplier to account for multiple shots per day
Apply.zn2.multi.sht="No"
if(Use.Date=="YES") Apply.zn2.multi.sht="Yes"

if(Apply.zn2.multi.sht=="Yes")
{
  Effort.monthly$Km.Gillnet.Days.inv=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.inv*1.6,Km.Gillnet.Days.inv))
  Effort.monthly$Km.Gillnet.Days.val=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.val*1.6,Km.Gillnet.Days.val))
  Effort.monthly$Km.Gillnet.Days.plus5=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.plus5*1.6,Km.Gillnet.Days.plus5))
  Effort.monthly$Km.Gillnet.Days.c=with(Effort.monthly,ifelse(zone=="Zone2",Km.Gillnet.Days.c*1.6,Km.Gillnet.Days.c))
}


#note: max is used to remove duplicated effort records from same return (i.e. from different species caught in that return)
#km gn days             
#max to remove duplicates
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)

#split boundary blocks
fn.split.boundary=function(DAT,What)
{
  id=match(What,colnames(DAT))
  a=subset(DAT,BLOCKX%in%Boundary.blk)
  b=subset(DAT,!BLOCKX%in%Boundary.blk)  
  x=a
  x$zone="Zone2"
  x[,id]=NA
  x[,id]=ifelse(!x$FINYEAR%in%c("2003-04","2004-05","2005-06"),0.5*a[,id],0.05*a[,id])
  a[,id]=ifelse(!a$FINYEAR%in%c("2003-04","2004-05","2005-06"),0.5*a[,id],0.95*a[,id])
  a=rbind(a,x)
  return(rbind(a,b))
}
Attach.Effort.monthly.c=fn.split.boundary(Attach.Effort.monthly.c,"Km.Gillnet.Days.c")
Attach.Effort.monthly=fn.split.boundary(Attach.Effort.monthly,"Km.Gillnet.Days.inv")

#sum to aggregage by zone and year
Attach.Effort.monthly.c=aggregate(Km.Gillnet.Days.c~zone+FINYEAR,data=Attach.Effort.monthly.c,sum,na.rm=T)
Attach.Effort.monthly=aggregate(Km.Gillnet.Days.inv~zone+FINYEAR,data=Attach.Effort.monthly,sum,na.rm=T)


#km gn hours
#max to remove duplicates
Attach.Effort.monthly.hrs.c=aggregate(Km.Gillnet.Hours.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)
Attach.Effort.monthly.hrs=aggregate(Km.Gillnet.Hours.inv~MONTH+BLOCKX+VESSEL+zone+FINYEAR,data=Effort.monthly,max,na.rm=T)


#split boundary blocks
Attach.Effort.monthly.hrs.c=fn.split.boundary(Attach.Effort.monthly.hrs.c,"Km.Gillnet.Hours.c")
Attach.Effort.monthly.hrs=fn.split.boundary(Attach.Effort.monthly.hrs,"Km.Gillnet.Hours.inv")

#sum to aggregage by zone and year
Attach.Effort.monthly.hrs.c=aggregate(Km.Gillnet.Hours.c~zone+FINYEAR,data=Attach.Effort.monthly.hrs.c,sum,na.rm=T)
Attach.Effort.monthly.hrs=aggregate(Km.Gillnet.Hours.inv~zone+FINYEAR,data=Attach.Effort.monthly.hrs,sum,na.rm=T)


#merge into single file
Attach.Effort.monthly.c=merge(Attach.Effort.monthly.c,Attach.Effort.monthly.hrs.c,by=c("zone","FINYEAR"),all=T)
Attach.Effort.monthly=merge(Attach.Effort.monthly,Attach.Effort.monthly.hrs,by=c("zone","FINYEAR"),all=T)
Attach.Effort.monthly=merge(Attach.Effort.monthly.c,Attach.Effort.monthly,by=c("zone","FINYEAR"),all.x=T)

rm(Attach.Effort.monthly.c,Attach.Effort.monthly.hrs.c,Attach.Effort.monthly.hrs)    


#2.3. Merge Monthly and daily
names(Attach.Effort.daily)[2]="FINYEAR"
Attach.Effort=rbind(Attach.Effort.monthly,Attach.Effort.daily)



#3. Aggregate gillnet effort for reporting in Sofar           
if(Use.Date=="YES") use.ID="NO"
if(Use.Date=="NO") use.ID="YES" #select if aggregating effort by SNo and DSNo


#3.1. Daily records
fn.Eff.Sofar=function(THESE.YRS)
{
  #3.1.1. Aggregate effort by year and zone
  SUBSET=subset(Effort.daily,finyear%in%THESE.YRS)
  
  #km gn days                         
  if(use.shots=="NO")
  {
    if(use.ID=="YES")
    {
      Km.gn.days=aggregate(netlen.c~ID+vessel+zone+finyear,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c~zone+finyear,data=Km.gn.days,sum,na.rm=T)
    }
    if(use.ID=="NO")
    {
      Km.gn.days=aggregate(netlen.c~finyear+vessel+date+zone,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c~finyear+zone,data=Km.gn.days,sum,na.rm=T)
    }
    
    Effort.Km.gn.days <- reshape(Km.gn.days,v.names="netlen.c",timevar="zone",idvar=c("finyear"),direction="wide")
    iid=match(c("netlen.c.West","netlen.c.Zone1","netlen.c.Zone2"),names(Effort.Km.gn.days))
    Effort.Km.gn.days[,iid]=Effort.Km.gn.days[,iid]/1e6
    names(Effort.Km.gn.days)[iid]=c("West","Zone1","Zone2")  
  }
  
  if(use.shots=="YES")
  {
    if(use.ID=="YES")
    {
      Km.gn.days=aggregate(netlen.c*shots.c~ID+vessel+zone+finyear,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c*shots.c~zone+finyear,data=Km.gn.days,sum,na.rm=T)
    }
    if(use.ID=="NO")
    {
      Km.gn.days=aggregate(netlen.c*shots.c~finyear+vessel+date+zone,data=SUBSET,max,na.rm=T)
      Km.gn.days=aggregate(netlen.c*shots.c~finyear+zone,data=Km.gn.days,sum,na.rm=T)
    }
    
    Effort.Km.gn.days <- reshape(Km.gn.days,v.names="Km.Gillnet.Days.c",timevar="zone",idvar=c("finyear"),direction="wide")
    iid=match(c("Km.Gillnet.Days.c.West","Km.Gillnet.Days.c.Zone1","Km.Gillnet.Days.c.Zone2"),names(Effort.Km.gn.days))
    Effort.Km.gn.days[,iid]=Effort.Km.gn.days[,iid]/1000
    names(Effort.Km.gn.days)[iid]=c("West","Zone1","Zone2")  
  }
  
  
  #km gn hours
  Km.gn.hours=aggregate(Km.Gillnet.Hours.c~finyear+vessel+date+ID+zone,data=SUBSET,max,na.rm=T)
  Km.gn.hours=aggregate(Km.Gillnet.Hours.c~finyear+zone,data=Km.gn.hours,sum,na.rm=T)
  Effort.Km.gn.hours <- reshape(Km.gn.hours,v.names="Km.Gillnet.Hours.c",timevar="zone",idvar=c("finyear"),direction="wide")
  iid=match(c("Km.Gillnet.Hours.c.West","Km.Gillnet.Hours.c.Zone1","Km.Gillnet.Hours.c.Zone2"),names(Effort.Km.gn.hours))
  Effort.Km.gn.hours[,iid]=Effort.Km.gn.hours[,iid]/1000
  names(Effort.Km.gn.hours)[iid]=c("West","Zone1","Zone2")
  
  
  #3.1.2. Calculate longline gillnet effort equivalent
  Daily.shark_rays.zone=subset(Data.daily,SPECIES%in%Elasmo.species & FINYEAR%in%THESE.YRS)
  
  #3.1.2.1 Annual catch by zone
  #GN                                
  Daily.Catch.shark_rays.zone=aggregate(LIVEWT~FINYEAR+zone,data=subset(Daily.shark_rays.zone,METHOD=="GN"),sum,na.rm=T)
  Daily.Catch.shark_rays.zone=reshape(Daily.Catch.shark_rays.zone,v.names="LIVEWT",timevar="zone",
                                      idvar=c("FINYEAR"),direction="wide")
  #LL
  Daily.Catch.shark_rays.zone.LL=aggregate(LIVEWT~FINYEAR+zone,data=subset(Daily.shark_rays.zone,METHOD=="LL"),sum,na.rm=T)
  Daily.Catch.shark_rays.zone.LL=reshape(Daily.Catch.shark_rays.zone.LL,v.names="LIVEWT",timevar="zone",
                                         idvar=c("FINYEAR"),direction="wide")
  Daily.Catch.shark_rays.zone.LL=Daily.Catch.shark_rays.zone.LL[order(Daily.Catch.shark_rays.zone.LL$FINYEAR),]
  Daily.Catch.shark_rays.zone.LL[is.na(Daily.Catch.shark_rays.zone.LL)]=0
  
  #3.1.2.2 nominal cpue
  FinYR=Daily.Catch.shark_rays.zone$FINYEAR
  Nom.cpue.days=Daily.Catch.shark_rays.zone[,iid]/Effort.Km.gn.days[,iid]
  Nom.cpue.hours=Daily.Catch.shark_rays.zone[,iid]/Effort.Km.gn.hours[,iid]
  
  #3.1.2.3 equiv effort
  #add columns to match dimensions, if no longline catch in some zone
  if(!ncol(Daily.Catch.shark_rays.zone.LL)==ncol(Daily.Catch.shark_rays.zone))
  {
    ii=which(!names(Daily.Catch.shark_rays.zone)%in%names(Daily.Catch.shark_rays.zone.LL))
    ADD=names(Daily.Catch.shark_rays.zone)[ii]
    mat=as.data.frame(matrix(ncol=length(ADD),nrow=nrow(Daily.Catch.shark_rays.zone)))
    names(mat)=ADD
    Daily.Catch.shark_rays.zone.LL=cbind(Daily.Catch.shark_rays.zone.LL,mat)
    Daily.Catch.shark_rays.zone.LL=Daily.Catch.shark_rays.zone.LL[,match(names(Daily.Catch.shark_rays.zone),names(Daily.Catch.shark_rays.zone.LL))]
    Daily.Catch.shark_rays.zone.LL[is.na(Daily.Catch.shark_rays.zone.LL)]=0
  }
  
  LL.equiv.Eff.days.zone=data.frame(FINYEAR=FinYR,Daily.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.days)
  LL.equiv.Eff.hours.zone=data.frame(FINYEAR=FinYR,Daily.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.hours)
  
  
  #3.1.3. Total effort
  #by year and zone
  Total.effort.zone.days=data.frame(FINYEAR=FinYR,Effort.Km.gn.days[,iid]+LL.equiv.Eff.days.zone[,iid])
  Total.effort.zone.hours=data.frame(FINYEAR=FinYR,Effort.Km.gn.hours[,iid]+LL.equiv.Eff.hours.zone[,iid])
  
  #by year and joint authority
  Total.effort.joint.days=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.days$West,
                                     JASGL=Total.effort.zone.days$Zone1+Total.effort.zone.days$Zone2)
  Total.effort.joint.hours=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.hours$West,
                                      JASGL=Total.effort.zone.hours$Zone1+Total.effort.zone.hours$Zone2)
  #by year and total
  Total.effort.days=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.days$WCGL+Total.effort.joint.days$JASGL)
  Total.effort.hours=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.hours$WCGL+Total.effort.joint.hours$JASGL)
  
  return(list(Total.effort.zone.days=Total.effort.zone.days,Total.effort.zone.hours=Total.effort.zone.hours,
              Total.effort.joint.days=Total.effort.joint.days,Total.effort.joint.hours=Total.effort.joint.hours,
              Total.effort.days=Total.effort.days,Total.effort.hours=Total.effort.hours))
}

Eff.Current.Yr=fn.Eff.Sofar(Current.yr)
Total.effort.zone.days=Eff.Current.Yr$Total.effort.zone.days
Total.effort.zone.hours=Eff.Current.Yr$Total.effort.zone.hours
Total.effort.joint.days=Eff.Current.Yr$Total.effort.joint.days
Total.effort.joint.hours=Eff.Current.Yr$Total.effort.joint.hours


#3.2. Monthly records                                                       

#3.2.1. Aggregate effort by year and zone
#km gn days
Effort.Km.gn.days.monthly=subset(Attach.Effort,zone%in%c("West","Zone1","Zone2"))
Effort.Km.gn.days.monthly=Effort.Km.gn.days.monthly[,match(c("zone","FINYEAR","Km.Gillnet.Days.c"),
                                                           names(Effort.Km.gn.days.monthly))]
Effort.Km.gn.days.monthly <- reshape(Effort.Km.gn.days.monthly,v.names="Km.Gillnet.Days.c",timevar="zone",
                                     idvar=c("FINYEAR"),direction="wide")
iid=match(c("Km.Gillnet.Days.c.West","Km.Gillnet.Days.c.Zone1","Km.Gillnet.Days.c.Zone2"),
          names(Effort.Km.gn.days.monthly))
Effort.Km.gn.days.monthly[,iid]=Effort.Km.gn.days.monthly[,iid]/1000
names(Effort.Km.gn.days.monthly)[iid]=c("West","Zone1","Zone2")

#km gn hours
Effort.Km.gn.hours.monthly=subset(Attach.Effort,zone%in%c("West","Zone1","Zone2"))
Effort.Km.gn.hours.monthly=Effort.Km.gn.hours.monthly[,match(c("zone","FINYEAR","Km.Gillnet.Hours.c"),
                                                             names(Effort.Km.gn.hours.monthly))]
Effort.Km.gn.hours.monthly <- reshape(Effort.Km.gn.hours.monthly,v.names="Km.Gillnet.Hours.c",
                                      timevar="zone",idvar=c("FINYEAR"),direction="wide")
Effort.Km.gn.hours.monthly[,iid]=Effort.Km.gn.hours.monthly[,iid]/1000
names(Effort.Km.gn.hours.monthly)[iid]=c("West","Zone1","Zone2")


# select financial years
#THESE.YRS=FINYEAR.monthly[1:31]
THESE.YRS=FINYEAR.monthly


#3.2.2. Calculate longline gillnet effort equivalent                              #REVIEW RORY
Monthly.shark_rays.zone=subset(Data.monthly,SPECIES%in%Elasmo.species & FINYEAR%in%THESE.YRS)

#3.2.2.1 Annual catch by zone
#GN                                
Monthly.Catch.shark_rays.zone=aggregate(LIVEWT~FINYEAR+zone,data=subset(Monthly.shark_rays.zone,METHOD=="GN"),sum,na.rm=T)
Monthly.Catch.shark_rays.zone=reshape(Monthly.Catch.shark_rays.zone,v.names="LIVEWT",timevar="zone",
                                      idvar=c("FINYEAR"),direction="wide")
#LL
Monthly.Catch.shark_rays.zone.LL=aggregate(LIVEWT~FINYEAR+zone,data=subset(Monthly.shark_rays.zone,METHOD=="LL"),
                                           sum,na.rm=T)
Monthly.Catch.shark_rays.zone.LL=reshape(Monthly.Catch.shark_rays.zone.LL,v.names="LIVEWT",timevar="zone",
                                         idvar=c("FINYEAR"),direction="wide")
Monthly.Catch.shark_rays.zone.LL=Monthly.Catch.shark_rays.zone.LL[order(Monthly.Catch.shark_rays.zone.LL$FINYEAR),]
Monthly.Catch.shark_rays.zone.LL[is.na(Monthly.Catch.shark_rays.zone.LL)]=0

#3.2.2.2 nominal cpue
FinYR=Monthly.Catch.shark_rays.zone$FINYEAR
II=match(THESE.YRS,Effort.Km.gn.days.monthly$FINYEAR)
Nom.cpue.days=Monthly.Catch.shark_rays.zone[,iid]/Effort.Km.gn.days.monthly[II,iid]
Nom.cpue.hours=Monthly.Catch.shark_rays.zone[,iid]/Effort.Km.gn.hours.monthly[II,iid]

#3.2.2.3 equiv effort
#add columns to match dimensions, if no longline catch in some zone
if(!ncol(Monthly.Catch.shark_rays.zone.LL)==ncol(Monthly.Catch.shark_rays.zone))
{
  ii=which(!names(Monthly.Catch.shark_rays.zone)%in%names(Monthly.Catch.shark_rays.zone.LL))
  ADD=names(Monthly.Catch.shark_rays.zone)[ii]
  mat=as.data.frame(matrix(ncol=length(ADD),nrow=nrow(Monthly.Catch.shark_rays.zone)))
  names(mat)=ADD
  Monthly.Catch.shark_rays.zone.LL=cbind(Monthly.Catch.shark_rays.zone.LL,mat)
  Monthly.Catch.shark_rays.zone.LL=Monthly.Catch.shark_rays.zone.LL[,match(names(Monthly.Catch.shark_rays.zone),
                                                                           names(Monthly.Catch.shark_rays.zone.LL))]
  Monthly.Catch.shark_rays.zone.LL[is.na(Monthly.Catch.shark_rays.zone.LL)]=0
}

LL.equiv.Eff.days.zone=data.frame(FINYEAR=FinYR,Monthly.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.days)
LL.equiv.Eff.hours.zone=data.frame(FINYEAR=FinYR,Monthly.Catch.shark_rays.zone.LL[,iid]/Nom.cpue.hours)


#3.2.3. Total effort (gillnet  plus longline equivalent. This is what's reported in Sofar)
#by year and zone
Total.effort.zone.days.monthly=data.frame(FINYEAR=FinYR,Effort.Km.gn.days.monthly[II,iid]+LL.equiv.Eff.days.zone[,iid])
Total.effort.zone.hours.monthly=data.frame(FINYEAR=FinYR,Effort.Km.gn.hours.monthly[II,iid]+LL.equiv.Eff.hours.zone[,iid])

#by year and joint authority
Total.effort.joint.days.monthly=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.days.monthly$West,
                                           JASGL=Total.effort.zone.days.monthly$Zone1+Total.effort.zone.days.monthly$Zone2)
Total.effort.joint.hours.monthly=data.frame(FINYEAR=FinYR,WCGL=Total.effort.zone.hours.monthly$West,
                                            JASGL=Total.effort.zone.hours.monthly$Zone1+Total.effort.zone.hours.monthly$Zone2)
#by year and total
Total.effort.days.monthly=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.days.monthly$WCGL+
                                       Total.effort.joint.days.monthly$JASGL)
Total.effort.hours.monthly=data.frame(FINYEAR=FinYR,Total=Total.effort.joint.hours.monthly$WCGL+
                                        Total.effort.joint.hours.monthly$JASGL)


Total.effort.days.monthly

Total.effort.hours.monthly