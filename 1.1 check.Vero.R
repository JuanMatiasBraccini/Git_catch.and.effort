fn.check.Vero=function(D,what)
{
  if(what=='fishery')
  {
    Tab=aggregate(LIVEWT~fishery+FINYEAR,D,sum)
    names(Tab)[3]="WT"
    Tab <- reshape(Tab, v.names = "WT", idvar = "fishery",
                   timevar = "FINYEAR", direction = "wide")
  }
  if(what=='FisheryCode')
  {
    Tab=aggregate(LIVEWT~FisheryCode+FINYEAR,D,sum)
    names(Tab)[3]="WT"
    Tab <- reshape(Tab, v.names = "WT", idvar = "FisheryCode",
                   timevar = "FINYEAR", direction = "wide")
  }
  
  Tab[is.na(Tab)]=0
  return(Tab)
  
}
fn.check.Vero.reap=function(D,what)
{
  if(what=='fishery')
  {
    Tab=aggregate(LIVEWT.reap~fishery+FINYEAR,D,sum)
    names(Tab)[3]="WT"
    Tab <- reshape(Tab, v.names = "WT", idvar = "fishery",
                   timevar = "FINYEAR", direction = "wide")
  }
  if(what=='FisheryCode')
  {
    Tab=aggregate(LIVEWT.reap~FisheryCode+FINYEAR,D,sum)
    names(Tab)[3]="WT"
    Tab <- reshape(Tab, v.names = "WT", idvar = "FisheryCode",
                   timevar = "FINYEAR", direction = "wide")
  }
  
  Tab[is.na(Tab)]=0
  return(Tab)
  
}

fn.plt.Vero=function(A,B)
{
  for(i in 1:nrow(A))
  {
    plot(unlist(A[i,2:ncol(A)]),ylab="catch",xlab="year",main=A$fishery[i])
    points(unlist(B[i,2:ncol(B)]),pch=19,col=2,cex=.7)
    legend("topright",c("CAESS","Matias"),
           pch=21,bty="n",col=c("black","white"),
           pt.bg=c("white","red"))
  } 
}

#Data.monthly has a few more kg because it resets
#dodgy livers and fins to Species 22999 and sets condition to "WF"
# but when removing fins and livers from Data.monthly.CAESS, those
# changes are not taken into accont

#Same logic applies to Data.daily

#upto lines 2248, both match
Data.monthly.CAESS$Same.return=with(Data.monthly.CAESS,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))

D.CAESS=subset(Data.monthly.CAESS,
               !SPECIES%in%c(22997,22998,22999))

D.caess=fn.check.Vero(D=D.CAESS,what='fishery')
D.mon=fn.check.Vero(D=subset(Data.monthly,
                  !SPECIES==22999&FINYEAR%in%unique(D.CAESS$FINYEAR)),what='FisheryCode')
ComPr=D.caess==D.mon
SAME=sum(ComPr)==nrow(D.caess)*ncol(D.caess)
SAME
if(SAME=="FALSE") fn.plt.Vero(A=D.caess,B=D.mon)

A=D.mon
A[,]=0
for(i in 1:nrow(A))A[i,]=D.mon[i,2:ncol(D.mon)]-D.caess[i,2:ncol(D.mon)]
which(A>0)   #if any record >0, then D.caess < D.mon
which(A<0)   #if any record <0, then D.caess > D.mon


 
a=subset(Data.monthly.CAESS,fishery=="SBS" & FINYEAR=="2002-03")
b=subset(Data.monthly,FisheryCode=="SBS" & FINYEAR=="2002-03")
sum(a$LIVEWT)
sum(b$LIVEWT.reap)



#upto lines 3836, both match
D.CAESS=subset(Data.monthly.CAESS,
               !SPECIES%in%c(22997,22998))

D.caess=fn.check.Vero(D=D.CAESS,what='fishery')
D.mon=fn.check.Vero.reap(D=subset(Data.monthly,
                FINYEAR%in%unique(D.CAESS$FINYEAR)),what='FisheryCode')
ComPr=D.caess==D.mon
SAME=sum(ComPr)==nrow(D.caess)*ncol(D.caess)
SAME
if(SAME=="FALSE") fn.plt.Vero(A=D.caess,B=D.mon)



a=subset(Data.monthly.CAESS,
         !SPECIES%in%c(22997,22998)&fishery=="CSFN",
         select=c(Same.return,FINYEAR,SPECIES,LIVEWT))
b=subset(Data.monthly,
    FINYEAR%in%unique(a$FINYEAR)&FisheryCode=="CSFN",
    select=c(Same.return,FINYEAR,SPECIES,LIVEWT.reap))
nn=unique(a$Same.return)
for(n in 1:length(nn))
{
  a1=subset(a,Same.return==nn[n])
  b1=subset(b,Same.return==nn[n])
  round(sum(a1$LIVEWT))
  round(sum(b1$LIVEWT.reap))
}


ay=aggregate(LIVEWT~FINYEAR,a,sum)
by=aggregate(LIVEWT.reap~FINYEAR,b,sum)
plot(ay[,2])
points(by[,2],pch=19,col=2,cex=.75)
