#---DATA SECTION-----
setwd("C:/Matias/Analyses/Catch and effort")

#OBSERVERS DATA
Survey=read.csv("Survey.weight.csv")


#FL to TL pars             
b.w=13.171   
b.g=4.6424
b.d=2.9835
b.s=0.2133

a.w=1.0044
a.g=1.0837
a.d=1.1849
a.s=1.2185

Survey$FL=with(Survey,ifelse(is.na(FL)| FL==0,CALCULATED.FL,FL))

Survey$FL=with(Survey,ifelse(is.na(FL)& !is.na(TL) & SPECIES=="WH",(TL-b.w)/a.w,
                             ifelse(is.na(FL)& !is.na(TL) & SPECIES=="GM",(TL-b.g)/a.g,
                                    ifelse(is.na(FL)& !is.na(TL) & SPECIES=="BW",(TL-b.d)/a.d,
                                           ifelse(is.na(FL)& !is.na(TL) & SPECIES=="TK",(TL-b.s)/a.s,FL)))))


Survey1=subset(Survey,SPECIES%in%c("GM","WH","BW","TK") & Mid.Lat<=(-26),select=c(SPECIES,Month,year,SEX,
          TL,FL,CALCULATED.FL,ZONE,Lat.round,Long.round,BOTDEPTH,Mid.Lat,Mid.Long))

Max.depth=aggregate(BOTDEPTH~ZONE+SPECIES,Survey1,max)
Min.depth=aggregate(BOTDEPTH~ZONE+SPECIES,Survey1,min)


fn.hist=function(dat,COL,RANGO)
{
  hist(dat$BOTDEPTH,col=COL,main="",ylab="",xlab="",xlim=c(0,max(Survey1$BOTDEPTH,na.rm=T)))
  legend("topright",paste('(',RANGO[1],'-',RANGO[2],' m)',sep=''),bty='n',cex=1.25)
  box()
}
 
SpecieS=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
Spec=c("WH","GM","BW","TK")
Zona=c("West coast","Zone 1","Zone 2")
Zns=c("WC","1","2")
COls=2:4

setwd("C:/Matias/Analyses/Depth distribution")
tiff(file="Depth distribution.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfrow=c(4,3),mar=c(1,1.25,1.5,1.25),oma=c(2.5,3,.1,.5),las=1,mgp=c(1.9,.55,0))
for(i in 1:length(Spec))
{
  for(j in 1:length(Zns))
  {
    Rang= c(subset(Min.depth,SPECIES==Spec[i] &ZONE==Zns[j])$BOTDEPTH,
            subset(Max.depth,SPECIES==Spec[i] &ZONE==Zns[j])$BOTDEPTH)
    fn.hist(subset(Survey1,SPECIES==Spec[i] & ZONE==Zns[j]),COls[j],Rang)
    if(i==1) mtext(Zona[j],3)
  }
    
  mtext(SpecieS[i],4,las=3,cex=1.25,line=0.5)
}
mtext("Depth (m)",side=1,line=1.25,font=1,las=0,cex=1.5,outer=T)
mtext("Frequency",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
dev.off()


fn.size.depth=function(dat,COL) plot(dat$BOTDEPTH,dat$FL,col=COL)
par(mfrow=c(4,3),mar=c(1,1.25,1.5,1.25),oma=c(2.5,3,.1,.5),las=1,mgp=c(1.9,.55,0))

for(i in 1:length(Spec))
{
  for(j in 1:length(Zns))
  {
    fn.size.depth(subset(Survey1,SPECIES==Spec[i] & ZONE==Zns[j]),COls[j])
    if(i==1) mtext(Zona[j],3)
  }
  
  mtext(SpecieS[i],4,las=3,cex=1.25,line=0.5)
}
mtext("Depth (m)",side=1,line=1.25,font=1,las=0,cex=1.5,outer=T)
mtext("Fork length (cm)",side=2,line=1,font=1,las=0,cex=1.5,outer=T)
