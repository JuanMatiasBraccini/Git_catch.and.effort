#CATCH COMPOSITION FROM SCIENTIFIC OBSERVATIONS #

#note:  script for analysing the observed catch composition of observed shots
#       this includes scientific surveys and trips on commercial gillnet and longline vessels

rm(list=ls(all=TRUE))
library(RODBC)				#include ODBC library for importing Acccess data



if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#DATA
setwd("M:/Fisheries Research/Production Databases/Shark")  

  
channel <- odbcConnectAccess("Sharks.mdb")			
Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F)   
Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)
Scalefish=sqlFetch(channel, "Scalefish", colnames = F)
close(channel)

Boat_hdr$year=as.numeric(substr(Boat_hdr$DATE,1,4))


setwd(handl_OneDrive("Data"))
Species_Code=read.csv("Species_names.csv")
Len.Wt.scale=read.csv("Length_Weights/SCALEY LENGTH-WEIGHT.csv")
Len.Wt.shark=read.csv("Length_Weights/SHARK LENGTH-WEIGHT.csv")





#PARAMETERS
#Growth pars (female)
Linf.w=120.7
Linf.d=374.4
Linf.s=244.2
K.w=0.369
K.d=.0367
K.s=.040
to.w=-0.544
to.d=-3.3
to.s=-4.8

growth.pars=data.frame(Linf=c(Linf.w,Linf.d,Linf.s),K=c(K.w,K.d,K.s),to=c(to.w,to.d,to.s))

#MAx age
max.age=c(20,55,39)



#PROCEDURE

#1. Merge data sets
useful.boat=c("SHEET_NO","BOAT","NET_LENGTH","NO HOOKS","Method","BLOCK","Month","year")
Boat=Boat_hdr[,match(useful.boat,names(Boat_hdr))]

useful.bio=c("SHEET_NO","SPECIES","TL","FL","PL","SEX")
Bio=Boat_bio[,match(useful.bio,names(Boat_bio))]

DATA=merge(Bio,Boat, by="SHEET_NO",all.x=T)
DATA$Scaley="NO"        #to separate duplicated codes

useful.scl=c("SHEET_NO","SPECIES","TL","SEX")
Bio.scl=Scalefish[,match(useful.scl,names(Scalefish))]
Bio.scl$Scaley="YES"
Bio.scl$PL=NA
Bio.scl$FL=NA

DATA.scl=merge(Bio.scl,Boat, by="SHEET_NO",all.x=T)
DATA=rbind(DATA,DATA.scl)

DATA$SPECIES=with(DATA,ifelse(SPECIES=="bt","BT",ifelse(SPECIES=="bw","BW",
                  ifelse(SPECIES=="tk","TK",ifelse(SPECIES=="wd","WD",
                  ifelse(SPECIES=="ww","WW",as.character(DATA$SPECIES)))))))

DATA=merge(DATA,Species_Code,all.x=T)

Len.wt=rbind(Len.Wt.shark,Len.Wt.scale)

DATA=subset(DATA,!(is.na(SPECIES)))
DATA=subset(DATA,!(SPECIES=="XX"))

  #correct typos in TL
id=which((DATA$SPECIES=="QS" & DATA$TL>100)| (DATA$Species_name=="West Australian Dhufish" & DATA$TL>200))
DATA=DATA[-id,]

  #scientific names list
List.sci.nam=DATA[!duplicated(DATA$Scientific_name),c(15:16)]

  #Data requests
    #Russel Hudson
These.sp=c("HH","HS","HZ","HW","HG","OW")
Russel=subset(DATA,SPECIES%in%These.sp,select=c(SPECIES, TL,FL))
Size.Freq=function(SP)
{
  dat=subset(Russel,SPECIES==SP)
  aa=data.frame(FL.bin.cm=NA,Number=NA,SPECIES=SP)
  if(nrow(dat)>2)
    {
      a=hist(dat$FL,plot=F,breaks=seq(0,max(dat$FL,na.rm=T)+10,by=10))
      aa=data.frame(FL.bin.cm=a$breaks[-length(a$breaks)],Number=a$counts,SPECIES=SP)
    }

  return(aa)
}
Russel.size=NULL
for (i in 1:length(These.sp)) Russel.size=rbind(Russel.size,Size.Freq(These.sp[i]))
id=subset(Species_Code,SPECIES%in%These.sp,select=c(SPECIES,Species_name))
Russel.size=merge(Russel.size,id,by="SPECIES")

write.csv(Russel.size,handl_OneDrive("Analyses/Catch and effort/Data_Resquests/Russel.observer.size.freq.csv"),row.names=F)



#2. Catch composition for commercial gillnet vessels
  #2.1. Subset data
Research.ves=c("FLIN","HAM","HOU","NAT")
Commercial=subset(DATA,!(BOAT%in%Research.ves)) #separate commercial

Commercial.gn=subset(Commercial,NET_LENGTH>0)   #separate gillnet

Commercial.gn=subset(Commercial.gn,BLOCK>=2600)  #south of 26


  #2.2 Size range by species
Uniq.sp=as.character(unique(Commercial.gn$SPECIES))
Uniq.sp=Commercial.gn[,match(c("SPECIES","Species_name"),names(Commercial.gn))]
Uniq.sp=Uniq.sp[!(duplicated(Uniq.sp$SPECIES)),]
Uniq.sp=subset(Uniq.sp,!(is.na(Species_name)))

Uniq.sp.nam=as.character(Uniq.sp[,2])
Uniq.sp=Uniq.sp[,1]

Problems=NULL
fn.size.by.sp=function(sp,sp.name)
{
  Data=subset(Commercial.gn,SPECIES==sp)
  Data$TL=ifelse(Data$TL<5,NA,Data$TL)  #remove nonsense small values
  Data$FL=ifelse(Data$FL<5,NA,Data$FL)  
  
  if(sp=="PS")Data$FL=  0.8460*Data$TL+0.3    #convert TL to FL for calculating weight (Derived from FishBase,except PS and RS)
  if(sp=="RS")Data$FL=  0.8670*Data$TL-20.9831
  if(sp=="NW")Data$FL=  0.9747*Data$TL+2e-14
  if(sp=="PA")Data$FL=  0.9149*Data$TL+4e-14
  if(sp=="SB")Data$FL=  0.9259*Data$TL
  if(sp=="SP")Data$FL=  0.9434*Data$TL-4.0321
  if(sp=="ST")Data$FL=  0.768*Data$TL-3e-14
  if(sp=="YT")Data$FL=  0.9025*Data$TL
  if(sp=="HH")Data$TL=  1.25*Data$FL
 # if(sp=="TG")Data$FL=  0.8761*Data$TL-13.3535


  TL.range=c(NA,NA)
  if(sum(is.na(Data$TL))<length(Data$TL)) TL.range=range(Data$TL,na.rm=T)
  
  

  #Ammend teleost records based on Dave F. observations
    #change max TL observed (out of range)
  if(sp.name=="Knife Jaw")TL.range[2]=48*1.1
  if(sp.name=="Knightfish")TL.range[2]=28*1.1
  if(sp.name=="Moonlighter")TL.range[2]=30*1.1
  if(sp.name=="Sea Sweep")TL.range[2]=30*1.1
  
    #set min length to minimum legal length to detect landing of undersized fish
  if(sp.name=="Breaksea Cod")TL.range[1]=30
  if(sp.name=="Baldchin Groper")TL.range[1]=40
  if(sp.name=="Blue Groper")TL.range[1]=50
  if(sp.name=="Coral Trout") TL.range[1]= 45
  if(sp.name=="Emperor (general)")TL.range[1]=28
  if(sp.name=="Estuary  Cod")TL.range[1]=40
  if(sp.name=="Pink Snapper")TL.range[1]=41
  if(sp.name=="Queen Snapper")TL.range[1]=41
  if(sp.name=="Red Emperor")TL.range[1]=41
  if(sp.name=="Red Snapper, Redfish, Bight Redfish, Nannygai")TL.range[1]=30
  if(sp.name=="Spangled Emperor")TL.range[1]=41
  if(sp.name=="Swallow Tail")TL.range[1]=30
  if(sp.name=="West Australian Dhufish")TL.range[1]=50


    
  FL.range=c(NA,NA)
  if(sum(is.na(Data$FL))<length(Data$FL)) FL.range=range(Data$FL,na.rm=T)
  
  Data$TW=NA
  TW.range=range(Data$TW)
  
  if(!(is.na(match(sp,Len.wt$CODE))))
  {
    len.wt=subset(Len.wt,CODE==sp)
    length=unique(as.character(len.wt$Length_type))
    if(sp=="TK")length="TL"
    id=match(length,colnames(Data))
    Data$TW=len.wt$a[1]*Data[,id]^len.wt$b[1]
    if(sp=="SA")Data$TW=len.wt$a[1]*(Data[,id]*10)^len.wt$b[1]
    Data$TW=ifelse(Data$TW<0.1,NA,Data$TW)  #remove nonsense small values
  # if(sum(is.na(Data$TW))<length(Data$TW)) TW.range=range(Data$TW,na.rm=T)
    
    if(length=="TL")TW.range=len.wt$a*TL.range^len.wt$b
    if(length=="TL" & sp=="SA")TW.range=len.wt$a*(TL.range*10)^len.wt$b
    if(length=="FL")TW.range=len.wt$a*FL.range^len.wt$b
    if(sp=="WW") TW.range[2]=75   #cap at 75 because length-weight pars give meaningless weights
    if(sp=="WR") TW.range[2]=5   #cap at 5 because length-weight pars give meaningless weights
    
    
    if(sum(is.na(Data[,id]))==length(Data[,id]))Problems<<-c(Problems,sp.name)
    
    if(sum(is.na(Data[,id]))<length(Data[,id]))
    {
      par(mfcol=c(3,1),mai=c(1.1,1.1,.2,.1),oma=c(.1,1,1,1))
      plot(Data[,id],Data$TW,xlab=paste(length,"(cm)"),ylab="Total weight (kg)")
      hist(Data[,id],xlab=paste(length,"(cm)"),main=sp.name)
      hist(Data$TW,xlab="Total weight (kg)",main=sp.name)
    }
  }
  
  return(list(TL.range=TL.range,FL.range=FL.range,TW.range=TW.range))
}

Size.list=vector('list',length=length(Uniq.sp))
names(Size.list)=Uniq.sp.nam
for (i in 1:length(Uniq.sp)) Size.list[[i]]=fn.size.by.sp(Uniq.sp[i],Uniq.sp.nam[i])

Size.matrix=matrix(unlist(Size.list),ncol=6,byrow=T)
colnames(Size.matrix)=c("TL.min","TL.max","FL.min","FL.max","TW.min","TW.max")
Size.matrix=as.data.frame(Size.matrix)
Size.matrix$Sname=Uniq.sp.nam





write.csv(Size.matrix,handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"),row.names=F)



#Size composition trends
Sel.species=c("WH","BW","TK")
id=match(Sel.species,Uniq.sp)
names(Sel.species)=Uniq.sp.nam[id]
size.at.birth=c(22,43,43) #from von B pars

fn.size.comp=function(sp,min.size,YR)
{
  Dat=subset(Commercial.gn,SPECIES==sp & FL>=min.size)
  Rango=range(Dat$FL,na.rm=T)
  Rango=c(floor(Rango[1]/10),ceiling(Rango[2]/10))*10
  
  Data=subset(Commercial.gn,SPECIES==sp & FL>=min.size & year%in%YR)
  if(nrow(Data)>2)
  {
    Seq=seq(Rango[1],Rango[2],10)
    a=hist(Data$FL,plot=F,breaks=Seq)
    probs=a$counts/sum(a$counts)
    return(data.frame(size=a$mids,prob=probs))
  }
  if(nrow(Data)<=2) return(0)
}

Size.props=vector('list',length(Sel.species))
names(Size.props)=names(Sel.species)
YRS=1993:2013

YRS.list=vector('list',length(YRS)+1)
names(YRS.list)=c("ALL",YRS)
YRS.list[[1]]=YRS
for (i in 2:length(YRS.list))YRS.list[[i]]=YRS[i-1]

# YRS.store=vector('list',length(YRS.list))
# names(YRS.store)=names(YRS.list)
YRS.store=YRS.list
for(i in 1:length(Sel.species))
{
  for(j in 1:length(YRS.store))  YRS.store[[j]]=fn.size.comp(Sel.species[i],size.at.birth[i],YRS.list[[j]])
  Size.props[[i]]=YRS.store
}


#selectivity (length to age conversion)
cv=0.2  #assumed CV of vonB function
Sel.fn=function(A,Linf,k,to,high.length,low.length,prop.sel)
{
  #Lower bound of age interval
  age.low=0:(A-1)
  
  #Upper bound of age interval
  age.up=1:A
  
  #von B function
  mean.length=Linf*(1-exp(-k*(age.low-to)))
  sd=mean.length*cv
  
  #prob of age given length
  prob=function(data)
  {
    pnorm(data,mean.length,sd)
  }
  high.bound=sapply(high.length,prob)
  low.bound=sapply(low.length,prob)
  mat.prob=high.bound-low.bound
  
  #normalised matrix
  st.mat.prob=t(t(mat.prob)/rowSums(t(mat.prob)))
  
  #calculate normalised selectivity
  Sel=rowSums(st.mat.prob*matrix(rep(prop.sel,nrow(st.mat.prob)),ncol=(ncol(st.mat.prob)),byrow=T))
  Sel=Sel/max(Sel)
  
  return(data.frame(age=age.low,pred.FL=mean.length,Sel=Sel))
}

Selectivity=vector('list',length(Sel.species))
names(Selectivity)=names(Sel.species)
for(i in 1:length(Sel.species))
  {
    Selectivity[[i]]=Sel.fn(max.age[i],growth.pars[i,1],growth.pars[i,2],growth.pars[i,3],
                  Size.props[[i]]$ALL$size*1.1,Size.props[[i]]$ALL$size*.9,Size.props[[i]]$ALL$prob)
  }

setwd(handl_OneDrive("Analyses/Catch and effort/Catch composition"))

SPECIES=c("whiskery","dusky","sandbar")

for(i in 1:length(Sel.species))
{
  tiff(file=paste(SPECIES[i],"Sel.tiff",sel=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.6,.5,.2,.1),oma=c(.1,.1,.1,.1))
  plot(Selectivity[[i]]$age, Selectivity[[i]]$pred.FL,xlab="age",ylab="Pred FL",main=SPECIES[i])
  plot(Selectivity[[i]]$age, Selectivity[[i]]$Sel,xlab="age",ylab="Selectivity")
  plot(Selectivity[[i]]$pred.FL, Selectivity[[i]]$Sel,xlab="PRed FL",ylab="Selectivity")
  dev.off()
  
  write.csv(Selectivity[[i]],paste(handl_OneDrive("Data/Selectivity/"),SPECIES[i],"Sel.csv",sel=""),row.names=F)
}