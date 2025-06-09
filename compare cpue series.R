# compare 2022 index used in SA and 2025 version with sensitivites -----------------------------------------------------------
Family='tw'
hndl.cmpr=handl_OneDrive('Analyses/Catch and effort/Outputs/Compare all index types/Compare 2022 with current/')

#Daily
for(s in Tar.sp)
{
  #Index used in 2022 assessment
  nm=capitalize(tolower(names(Best.Model.daily)[s]))
  if(nm=="Dusky whaler") nm="Dusky shark"
  prev.indx=read.csv(handl_OneDrive(paste0('Analyses/Data_outs/',nm,"/",nm,".annual.abundance.basecase.daily.csv")))%>%
    mutate(finyear=Finyear,
           response= Mean,
           Model='fit to SA',
           lower.CL=LOW.CI,
           upper.CL=UP.CI)%>%
    dplyr::select(Model,finyear,response,lower.CL,upper.CL)
  
  #current index
  Fmula=list(Best=Best.Model.daily[[s]]) 
  if(names(Best.Model.daily)[s]=='Gummy Shark')
  {
    Fmula$Old.best=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) + step.mcal_target_group')
    Fmula$No.targeting=formula("cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+s(long10.corner,lat10.corner)")
    Fmula$Year.continuous=formula('cpue ~ s(year) + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) + step.mcal_target_group + mesh')
  }
  if(names(Best.Model.daily)[s]%in%c("Whiskery Shark","Dusky Whaler"))
  {
    Fmula$Old.best=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) +s(mean.depth)+ step.mcal_target_group')
    Fmula$No.targeting=formula("cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+s(long10.corner,lat10.corner)+s(mean.depth)")
    Fmula$Year.continuous=formula('cpue ~ s(year) + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) + step.mcal_target_group + shots.c + mesh')
  }
  if(names(Best.Model.daily)[s]%in%c("Sandbar Shark"))
  {
    Fmula$Old.best=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) +s(mean.depth)+ step.mcal_target_group')
    Fmula$No.targeting=formula("cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+s(long10.corner,lat10.corner)+s(mean.depth)")
    Fmula$Year.continuous=formula('cpue ~ s(year) + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) + step.mcal_target_group')
  }
  d=DATA.list.LIVEWT.c.daily[[s]]%>%
    filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])%>%
    mutate(shots.c=ifelse(shots.c>2,2,shots.c),
           nlines.c=ifelse(nlines.c>3,3,nlines.c),
           mesh=ifelse(!mesh%in%c(165,178),'other',mesh),
           year=as.numeric(substr(FINYEAR,1,4)))
  
  d.indi=DATA.list.LIVEWT.c.daily[[s]]%>%
    filter(VESSEL%in%VES.used.daily.indi[[s]] & BLOCKX%in%BLKS.used.daily.indi[[s]])%>%
    mutate(shots.c=ifelse(shots.c>2,2,shots.c),
           nlines.c=ifelse(nlines.c>3,3,nlines.c),
           mesh=ifelse(!mesh%in%c(165,178),'other',mesh),
           year=as.numeric(substr(FINYEAR,1,4)))
  
  Terms=Predictors_daily
  Continuous=Covariates.daily
  colnames(d)=tolower(colnames(d))
  colnames(d.indi)=tolower(colnames(d.indi))
  Terms=tolower(Terms)
  Continuous=tolower(Continuous)
  Factors=Terms[!Terms%in%Continuous]
  d=d%>%
    dplyr::select(c(catch.target,km.gillnet.hours.c,any_of(Terms),year))%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  d <- makecategorical(Factors,d)
  
  d.indi=d.indi%>%
    dplyr::select(c(catch.target,km.gillnet.hours.c,any_of(Terms),year))%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  d.indi <- makecategorical(Factors,d.indi)
  
  
  MODS=foreach(i=1:length(Fmula),.packages=c('dplyr','mgcv')) %dopar%
    {
      mod<-bam(formula(Fmula[[i]]),data=d,family=Family,select=TRUE,method="fREML",discrete=TRUE)
      return(mod)
    }
  names(MODS)=names(Fmula)
  
  Fmula$ind.ves_blk=Fmula$Best
  MODS$ind.ves_blk=bam(formula(Fmula$ind.ves_blk),data=d.indi,family=Family,select=TRUE,method="fREML",discrete=TRUE)
  
  Preds=vector('list',length(MODS))
  for(x in 1:length(MODS))
  {
    plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,
         pages=1,rug=FALSE)
    mtext(Fmula[[x]],side=3, line = 2, adj = 0.65,col='brown',cex=.8)
    
    if(names(MODS)[x]=="Year.continuous")
    {
      Preds[[x]]=pred.fun.continuous(d=d,mod=MODS[[x]],PRED='year',Formula=Fmula[[x]])%>%
                  mutate(finyear=as.integer(year),
                         Model=names(MODS)[x],
                         response=Pred,
                         lower.CL=Pred-1.96*Pred.SE,upper.CL=Pred+1.96*Pred.SE)%>%
        group_by(Model,finyear)%>%
        summarise(response=mean(response),
                  lower.CL=mean(lower.CL),
                  upper.CL=mean(upper.CL))%>%
        mutate(finyear=factor(finyear))
    }else
    {
      Preds[[x]]=pred.fun(mod=MODS[[x]],biascor="YES",PRED="finyear")%>%mutate(Model=names(MODS)[x])%>%
        dplyr::select(Model,finyear,response,lower.CL,upper.CL)
    }
      
    
  }
  names(Preds)=names(MODS)
  Preds$prev.indx=prev.indx
  
  Preds$No.bias.corr=pred.fun(mod=MODS$Best,biascor="NO",PRED="finyear")%>%mutate(Model='No.bias.corr')%>%
                        dplyr::select(Model,finyear,response,lower.CL,upper.CL)
  
  p1=do.call(rbind,Preds)%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    ggplot(aes(year,response,color=Model))+
    geom_line(alpha=0.35,linetype='dashed')+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
    geom_point(size=2)+
    theme_PA()+theme(legend.position = 'top')+ylim(0,NA)+
    ylab('CPUE (+/- 95% CI)')+xlab('Financial year')+ggtitle("Normal scale")
  
  
  p2=do.call(rbind,Preds)%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    group_by(Model)%>%
    mutate(MEAN=mean(response),response=response/MEAN, lower.CL=lower.CL/MEAN, upper.CL=upper.CL/MEAN)%>%
    ggplot(aes(year,response,color=Model))+
    geom_line(alpha=0.35,linetype='dashed')+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
    geom_point(size=2)+
    theme_PA()+theme(legend.position = 'top')+ylim(0,NA)+
    ylab('CPUE (+/- 95% CI)')+xlab('Financial year')+ggtitle("Relative scale")
  ggarrange(p1,p2,ncol=1,common.legend = T)
  ggsave(paste0(hndl.cmpr,nm,'_daily.tiff'), width = 6,height = 6, dpi = 300, compression = "lzw")
}

#Monthly
for(s in Tar.sp)
{
  #Index used in 2022 assessment
  nm=capitalize(tolower(names(Best.Model)[s]))
  if(nm=="Dusky whaler") nm="Dusky shark"
  prev.indx=read.csv(handl_OneDrive(paste0('Analyses/Data_outs/',nm,"/",nm,".annual.abundance.basecase.monthly.csv")))%>%
    mutate(finyear=Finyear,
           response= Mean,
           Model='fit to SA',
           lower.CL=LOW.CI,
           upper.CL=UP.CI)%>%
    dplyr::select(Model,finyear,response,lower.CL,upper.CL)
  
  #current index
  Fmula=list(Best=Best.Model[[s]]) 
  d=DATA.list.LIVEWT.c[[s]]%>%
    filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])%>%
    mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
  d.indi=DATA.list.LIVEWT.c[[s]]%>%
    filter(VESSEL%in%VES.used.indi[[s]] & BLOCKX%in%BLKS.used.indi[[s]])%>%
    mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
  
  Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
  Continuous=Covariates.monthly
  
  colnames(d)=tolower(colnames(d))
  colnames(d.indi)=tolower(colnames(d.indi))
  Terms=tolower(Terms)
  Continuous=tolower(Continuous)
  Factors=Terms[!Terms%in%Continuous]
  d=d%>%
    dplyr::select(c(catch.target,km.gillnet.hours.c,any_of(Terms)))%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  d <- makecategorical(Factors,d)
  d.indi=d.indi%>%
    dplyr::select(c(catch.target,km.gillnet.hours.c,any_of(Terms)))%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  d.indi <- makecategorical(Factors,d.indi)
  
  MODS=foreach(i=1:length(Fmula),.packages=c('dplyr','mgcv')) %dopar%
    {
      mod<-bam(formula(Fmula[[i]]),data=d,family=Family,select=TRUE,method="fREML",discrete=TRUE)
      return(mod)
    }
  names(MODS)=names(Fmula)
  
  Fmula$ind.ves_blk=Fmula$Best
  MODS$ind.ves_blk=bam(formula(Fmula$ind.ves_blk),data=d.indi,family=Family,select=TRUE,method="fREML",discrete=TRUE)
  
  if(names(DATA.list.LIVEWT.c)[s]%in%c("Sandbar Shark","Gummy Shark"))
  {
    Fmula$drop.year=Fmula$Best
    if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")dd=d%>%filter(!finyear%in%c('1975-76'))
    if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")dd=d%>%filter(!finyear%in%c('1986-87','1987-88','1988-89'))
    MODS$drop.year<-bam(formula(Fmula$drop.year),data=droplevels(dd),family=Family,select=TRUE,method="fREML",discrete=TRUE)
  }

  
  
  Preds=vector('list',length(MODS))
  for(x in 1:length(MODS))
  {
    plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,
         pages=1,rug=FALSE)
    mtext(Fmula[[x]],side=3, line = 2, adj = 0.65,col='brown',cex=.8)
    
    Preds[[x]]=pred.fun(mod=MODS[[x]],biascor="YES",PRED="finyear")%>%mutate(Model=names(MODS)[x])%>%
      dplyr::select(Model,finyear,response,lower.CL,upper.CL)
    
  }
  names(Preds)=names(MODS)
  Preds$prev.indx=prev.indx
  Preds$No.bias.corr=pred.fun(mod=MODS$Best,biascor="NO",PRED="finyear")%>%mutate(Model='No.bias.corr')%>%
    dplyr::select(Model,finyear,response,lower.CL,upper.CL)
  
  p1=do.call(rbind,Preds)%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    ggplot(aes(year,response,color=Model))+
    geom_line(alpha=0.35,linetype='dashed')+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
    geom_point(size=2)+
    theme_PA()+theme(legend.position = 'top')+ylim(0,NA)+
    ylab('CPUE (+/- 95% CI)')+xlab('Financial year')+ggtitle("Normal scale")
  
  
  p2=do.call(rbind,Preds)%>%
    mutate(year=as.numeric(substr(finyear,1,4)))%>%
    group_by(Model)%>%
    mutate(MEAN=mean(response),response=response/MEAN, lower.CL=lower.CL/MEAN, upper.CL=upper.CL/MEAN)%>%
    ggplot(aes(year,response,color=Model))+
    geom_line(alpha=0.35,linetype='dashed')+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
    geom_point(size=2)+
    theme_PA()+theme(legend.position = 'top')+ylim(0,NA)+
    ylab('CPUE (+/- 95% CI)')+xlab('Financial year')+ggtitle("Relative scale")
  ggarrange(p1,p2,ncol=1,common.legend = T)
  ggsave(paste0(hndl.cmpr,nm,'.tiff'), width = 6,height = 6, dpi = 300, compression = "lzw")
}

# Explore wigglyness of month and depth -----------------------------------------------------------
    #Daily
for(s in Tar.sp)  
{

  #current index
  Fmula=list(Best=Best.Model.daily[[s]]) 
  Fmula$free_month=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
    s(long10.corner, lat10.corner) + step.mcal_target_group+  mesh')
  
  if(names(Best.Model.daily)[s]%in%c("Whiskery Shark","Dusky Whaler"))
  {
    Fmula$free_month=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
    s(long10.corner, lat10.corner) +s(mean.depth)+ step.mcal_target_group+shots.c + mesh')
    Fmula$k_3_depth=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, k = 12, bs = "cc") + 
    s(long10.corner, lat10.corner) +s(mean.depth,k=3)+ step.mcal_target_group+shots.c + mesh')
    
  }
  if(names(Best.Model.daily)[s]%in%c("Sandbar Shark"))
  {
    Fmula$free_month=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
    s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group')
    Fmula$k_3_depth=formula('cpue ~ finyear + s(vessel, bs = "re") + s(month, k=12, bs = "cc") + 
    s(long10.corner, lat10.corner) + s(mean.depth, k=3) + step.mcal_target_group')
  }
  
  
  d=DATA.list.LIVEWT.c.daily[[s]]%>%
    filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])%>%
    mutate(shots.c=ifelse(shots.c>2,2,shots.c),
           nlines.c=ifelse(nlines.c>3,3,nlines.c),
           mesh=ifelse(!mesh%in%c(165,178),'other',mesh),
           year=as.numeric(substr(FINYEAR,1,4)))
  
  Terms=Predictors_daily
  Continuous=Covariates.daily
  colnames(d)=tolower(colnames(d))
  Terms=tolower(Terms)
  Continuous=tolower(Continuous)
  Factors=Terms[!Terms%in%Continuous]
  d=d%>%
    dplyr::select(c(catch.target,km.gillnet.hours.c,any_of(Terms),year))%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  d <- makecategorical(Factors,d)
  MODS=foreach(i=1:length(Fmula),.packages=c('dplyr','mgcv')) %dopar%
    {
      mod<-bam(formula(Fmula[[i]]),data=d,family=Family,select=TRUE,method="fREML",discrete=TRUE)
      return(mod)
    }
  names(MODS)=names(Fmula)
  
  
  Preds=vector('list',length(MODS))
  for(x in 1:length(MODS))
  {
    plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,
         pages=1,rug=FALSE)
    mtext(Fmula[[x]],side=3, line = 2, adj = 0.65,col='brown',cex=.8)

  }
}

# Some old crap -----------------------------------------------------------
do.old.crap=FALSE
if(do.old.crap)
{
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
  
  if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
  
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Compare all index types"))
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
  
  
  
}

