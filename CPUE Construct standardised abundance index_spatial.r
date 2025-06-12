#-- Select model structure ----------------------------------------------
add.CITES.term=TRUE  #CITES listing term for smooth HH
#remove predictors identified as highly correlated
cNSTNT=Categorical[!Categorical=="block10"] 
cNSTNT.daily=Categorical[!Categorical=="blockx"]
if(any(c(do_cluster,do_pca,do_Stephens_McCall)=='YES')) cNSTNT.daily=c(cNSTNT.daily,Targeting.vars)

#define best model   
Best.Model.zone=vector('list',length(SP.list)) 
names(Best.Model.zone)=names(SP.list)
Best.Model.zone.daily=Best.Model.zone

if(Def.mod.Str=="YES")     #takes 25 mins 
{
  hndl.modl.sel=handl_OneDrive("Analyses/Catch and effort/Outputs/Model Selection/")
  if(Use.Tweedie)       
  {
    tic()
    #1. Indicator species   
    foreach(s=Tar.sp,.packages=c('dplyr','mgcv','doParallel')) %do%
    {
      zns=sort(unique(DATA.list.LIVEWT.c[[s]]$zone))
      if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #outside core area so not enough observations, issues with deg. freedom 
      if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")  
        
      for(z in 1:length(zns))
      {
            #1.1. Monthly
            if(do.monthly.def.str)
            {
              d=DATA.list.LIVEWT.c[[s]]%>%  
                filter(VESSEL%in%VES.used.zone[[s]][[z]] & BLOCKX%in%BLKS.used.zone[[s]][[z]] & zone==zns[z])%>%
                mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
              if(Exclude.yr.gummy) if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
              if(Exclude.yr.sandbar) if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
              Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
              Continuous=Covariates.monthly
              Fmula=list(Year.block="cpue~finyear + blockx",
                         Month=     "cpue~finyear + blockx + s(month,k=12,bs='cc')",
                         Shots=     "cpue~finyear + blockx + s(month,k=12,bs='cc') + shots.c",
                         Vessel=    "cpue~finyear + blockx + s(month,k=12,bs='cc') + shots.c + s(vessel,bs='re')") 
              Family='tw'
              FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],'_',zns[z],"_GAM_monthly",sep="")
              colnames(d)=tolower(colnames(d))
              Terms=tolower(Terms)
              Continuous=tolower(Continuous)
              Factors=Terms[!Terms%in%Continuous]
              d=d%>%
                dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
                mutate(cpue=catch.target/km.gillnet.hours.c)
              d <- makecategorical(Factors,d)
              MODS=foreach(i=1:length(Fmula),.packages=c('dplyr','mgcv','doParallel')) %dopar%
                {
                  mod<-gam(formula(Fmula[[i]]),data=d,family=Family,select=TRUE,method="REML")
                  return(mod)
                }
              names(MODS)=names(Fmula)
              pdf(paste(hndl.modl.sel,FILE,".pdf",sep=''))
              mod.sumery=vector('list',length(MODS))
              Explained.dev=Preds=mod.sumery
              for(x in 1:length(MODS))
              {
                plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,
                     pages=1,rug=FALSE)
                mtext(Fmula[[x]],side=3, line = 2, adj = 0.65,col='brown',cex=.8)
                mod.sumery[[x]]=tidy(MODS[[x]])
                Preds[[x]]=pred.fun(mod=MODS[[x]],biascor="YES",PRED="finyear")%>%mutate(Model=names(MODS)[x])
                if(nrow(mod.sumery[[x]])>0)
                {
                  grid.newpage()
                  class(mod.sumery[[x]]) <- "data.frame"
                  grid.table(mod.sumery[[x]])
                }
                Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100  
              }
              #Model prediction by each term
              print(do.call(rbind,Preds)%>%
                      mutate(year=as.numeric(substr(finyear,1,4)))%>%
                      ggplot(aes(year,response,color=Model))+
                      geom_line(alpha=0.35,linetype='dashed')+
                      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
                      geom_point(size=2)+
                      theme_PA()+theme(legend.position = 'top')+
                      ylab('CPUE (+/- 95% CI)')+xlab('Financial year'))
              
              #explained deviance by each term
              Explained.dev=data.frame(Model=names(Fmula),
                                       Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
                mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
              grid.newpage()
              class(Explained.dev) <- "data.frame"
              grid.table(Explained.dev)
              dev.off()
            }
            
            #1.2. Daily
            d=DATA.list.LIVEWT.c.daily[[s]]%>%
              filter(VESSEL%in%VES.used.zone.daily[[s]][[z]] & BLOCKX%in%BLKS.used.zone.daily[[s]][[z]] & zone==zns[z])%>%
              mutate(shots.c=ifelse(shots.c>2,2,shots.c),
                     nlines.c=ifelse(nlines.c>3,3,nlines.c),
                     mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
            Terms=Predictors_daily
            Continuous=Covariates.daily
            Fmula=list(Year=      "cpue~finyear",
                       Shots=     "cpue~finyear+shots.c",
                       Mesh=      "cpue~finyear+shots.c+mesh",
                       Nlines=    "cpue~finyear+shots.c+mesh+nlines.c",
                       Targeting= "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group",
                       Month=     "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group+s(month,k=12,bs='cc')",
                       Depth=     "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group+s(month,k=12,bs='cc')+s(mean.depth)",
                       Lunar=     "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group+s(month,k=12,bs='cc')+s(mean.depth)+s(lunar,k=5,bs='cc')",
                       Vessel=    "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group+s(month,k=12,bs='cc')+s(mean.depth)+s(lunar,k=5,bs='cc')+s(vessel,bs='re')",
                       Lat.long=  "cpue~finyear+shots.c+mesh+nlines.c+step.mcal_target_group+s(month,k=12,bs='cc')+s(mean.depth)+s(lunar,k=5,bs='cc')+s(vessel,bs='re')+s(long10.corner,lat10.corner)")
            Family='tw'
            FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],'_',zns[z],"_GAM_daily",sep="")
            colnames(d)=tolower(colnames(d))
            Terms=tolower(Terms)
            Continuous=tolower(Continuous)
            Factors=Terms[!Terms%in%Continuous]
            d=d%>%
              dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
              mutate(cpue=catch.target/km.gillnet.hours.c)
            d <- makecategorical(Factors,d)
            MODS=foreach(i=1:length(Fmula),.packages=c('dplyr','mgcv')) %dopar%
              {
                mod<-gam(formula(Fmula[[i]]),data=d,family=Family,select=TRUE,method="REML")
                return(mod)
              }
            names(MODS)=names(Fmula)
            pdf(paste(hndl.modl.sel,FILE,".pdf",sep=''))
            mod.sumery=vector('list',length(MODS))
            Explained.dev=Preds=mod.sumery
            for(x in 1:length(MODS))
            {
              plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,
                   pages=1,rug=FALSE)
              mtext(Fmula[[x]],side=3, line = 2, adj = 0.65,col='brown',cex=.8)
              mod.sumery[[x]]=tidy(MODS[[x]])
              Preds[[x]]=pred.fun(mod=MODS[[x]],biascor="YES",PRED="finyear")%>%mutate(Model=names(MODS)[x])
              if(nrow(mod.sumery[[x]])>0)
              {
                grid.newpage()
                class(mod.sumery[[x]]) <- "data.frame"
                grid.table(mod.sumery[[x]])
              }
              Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100  
            }
            
            #Model prediction by each term
            print(do.call(rbind,Preds)%>%
                    mutate(year=as.numeric(substr(finyear,1,4)))%>%
                    ggplot(aes(year,response,color=Model))+
                    geom_line(alpha=0.35,linetype='dashed')+
                    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
                    geom_point(size=2)+
                    theme_PA()+theme(legend.position = 'top')+
                    ylab('CPUE (+/- 95% CI)')+xlab('Financial year'))
            
            #explained deviance by each term
            Explained.dev=data.frame(Model=names(Fmula),
                                     Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
              mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
            grid.newpage()
            class(Explained.dev) <- "data.frame"
            grid.table(Explained.dev)
            dev.off()
            
          }
    }
    toc()
  }
}
if(Def.mod.Str=="NO")   
{
  if(Use.Tweedie)
  {
    for(s in Tar.sp)
    {
      #Monthly
      zns=sort(unique(DATA.list.LIVEWT.c[[s]]$zone))
      if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #outside core area so not enough observations, issues with deg. freedom 
      if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")  
      
      dummy=vector('list',length(zns))
      names(dummy)=zns
      for(z in 1:length(zns))
      {
        dummy[[z]]=formula(cpue~finyear + blockx + s(vessel,bs='re') + s(month,k=12,bs="cc"))
        if(names(SP.list)[s]=="Gummy Shark" & zns[z]=='Zone2')  dummy[[z]]=formula(cpue~finyear + blockx + s(vessel,bs='re'))
      }
      Best.Model.zone[[s]]=dummy
      
      #Daily
      for(z in 1:length(zns))
      {
        dummy[[z]]=formula(cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
                             s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group)
        if(names(SP.list)[s]=="Gummy Shark" & zns[z]=='Zone1') dummy[[z]]=formula(cpue ~ finyear + s(month, bs = "cc") + 
                                                                                    s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group+shots.c+mesh+nlines.c)
        if(names(SP.list)[s]=="Gummy Shark" & zns[z]=='Zone2')  dummy[[z]]=formula(cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                                                     s(long10.corner, lat10.corner)+mesh)
        if(names(SP.list)[s]=="Dusky Whaler" & zns[z]=='Zone1')  dummy[[z]]=formula(cpue ~ finyear + s(month, bs = "cc") + 
                                                                                     s(long10.corner, lat10.corner) + s(mean.depth) +shots.c+mesh)
        if(names(SP.list)[s]=="Dusky Whaler" & zns[z]=='West')  dummy[[z]]=formula(cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                                                     s(long10.corner, lat10.corner) + s(mean.depth)+ step.mcal_target_group +mesh+nlines.c)
        if(names(SP.list)[s]=="Whiskery Shark" & zns[z]=='West')  dummy[[z]]=formula(cpue ~ finyear + s(month, bs = "cc") + 
                                                                                       s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group + mesh)
        if(names(SP.list)[s]=="Whiskery Shark" & zns[z]=='Zone1')  dummy[[z]]=formula(cpue ~ finyear  + s(month, bs = "cc") + 
                                                                                        s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group+shots.c+mesh)
        if(names(SP.list)[s]=="Sandbar Shark" & zns[z]=='West')  dummy[[z]]=formula(cpue ~ finyear + s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                                                       s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group+mesh)
        if(names(SP.list)[s]=="Sandbar Shark" & zns[z]=='Zone1')  dummy[[z]]=formula(cpue ~ finyear  + s(month, bs = "cc") + 
                                                                                        s(long10.corner, lat10.corner) + s(mean.depth) + step.mcal_target_group+shots.c+mesh)
      }
      Best.Model.zone.daily[[s]]=dummy
    }
  }
}

#-- Run standardisation on selected model and plot indices  ----------------------------------------------
if(Use.Delta)  
{
  system.time({for(s in match(TARGETS[-match(17001,TARGETS)],SP.list))  
  {
  if(Use.Qualif.level)
  {
    #Monthly
    DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used[[s]])
    ZonEs=unique(DAT$zone)
    if(names(SP.list)[s]=="Sandbar Shark") ZonEs=subset(ZonEs,!ZonEs=="Zone2")   #no enough observations in zone2
    pred.temp=vector('list',length(ZonEs))
    names(pred.temp)=ZonEs
    pred.temp.crip=pred.temp.nrm=pred.temp
    
    for(z in 1:length(ZonEs))
    {
      #1. Fit model to zone data
      model=fn.stand(d=subset(DAT,zone==ZonEs[z]),Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                     efrt="km.gillnet.hours.c",Formula=Best.Model[[s]],Formula.gam=NULL)
      
      #2. Predict for selected blocks
      d=model$DATA
      Pred.1=pred.fun(MOD=model$res,biascor="YES",PRED="finyear",Pred.type="link")
      
      #3. Apply creep
      Pred.1.c=Pred.1
      add.crp=Eff.creep$effort.creep[match(Pred.1.c$finyear,Eff.creep$finyear)]
      Pred.1.c$response=Pred.1.c$response*(1-add.crp)
      Pred.1.c$lower.CL=Pred.1.c$lower.CL*(1-add.crp)
      Pred.1.c$upper.CL=Pred.1.c$upper.CL*(1-add.crp)
      
      #4.Normalize
      Pred.1.c.norm=Pred.1.c
      Mn=mean(Pred.1.c.norm$response)
      Pred.1.c.norm$response=Pred.1.c.norm$response/Mn
      Pred.1.c.norm$lower.CL=Pred.1.c.norm$lower.CL/Mn
      Pred.1.c.norm$upper.CL=Pred.1.c.norm$upper.CL/Mn
      
      #5.Store
      pred.temp[[z]]=Pred.1
      pred.temp.crip[[z]]=Pred.1.c
      pred.temp.nrm[[z]]=Pred.1.c.norm
      rm(d,Pred.1,Pred.1.c,Pred.1.c.norm)
    }
    Pred.zone[[s]]=pred.temp
    Pred.zone.creep[[s]]=pred.temp.crip
    Pred.zone.nrmlzd[[s]]=pred.temp.nrm
    
    
    #Daily
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    ZonEs=unique(DAT$zone)
    if(names(SP.list)[s]=="Sandbar Shark") ZonEs=subset(ZonEs,!ZonEs=="Zone2")   #no enough observations in zone2
    pred.temp=vector('list',length(ZonEs))
    names(pred.temp)=ZonEs
    pred.temp.crip=pred.temp.nrm=pred.temp
    
    for(z in 1:length(ZonEs))
    {
      #1. Fit model to zone data
      model=fn.stand(d=subset(DAT,zone==ZonEs[z]),Response="catch.target",RESPNS="LNcpue",
                     PREDS=Predictors_daily,efrt="km.gillnet.hours.c",
                     Formula=NULL,Formula.gam=Best.Model.daily.gam[[s]])
      
      #2. Predict for selected blocks
      d=model$DATA
      Pred.1=pred.fun(MOD=model$res.gam,biascor="YES",PRED="finyear",Pred.type="link")
      
      #3. Apply creep
      Pred.1.c=Pred.1
      add.crp=Eff.creep$effort.creep[match(Pred.1.c$finyear,Eff.creep$finyear)]
      Pred.1.c$response=Pred.1.c$response*(1-add.crp)
      Pred.1.c$lower.CL=Pred.1.c$lower.CL*(1-add.crp)
      Pred.1.c$upper.CL=Pred.1.c$upper.CL*(1-add.crp)
      
      #4.Normalize
      Pred.1.c.norm=Pred.1.c
      Mn=mean(Pred.1.c.norm$response)
      Pred.1.c.norm$response=Pred.1.c.norm$response/Mn
      Pred.1.c.norm$lower.CL=Pred.1.c.norm$lower.CL/Mn
      Pred.1.c.norm$upper.CL=Pred.1.c.norm$upper.CL/Mn
      
      #5.Store
      pred.temp[[z]]=Pred.1
      pred.temp.crip[[z]]=Pred.1.c
      pred.temp.nrm[[z]]=Pred.1.c.norm
      rm(d,Pred.1,Pred.1.c,Pred.1.c.norm)
    }
    Pred.daily.zone[[s]]=pred.temp
    Pred.daily.zone.creep[[s]]=pred.temp.crip
    Pred.daily.zone.nrmlzd[[s]]=pred.temp.nrm
    
  }
  
  #monthly
  system.time({Zone_preds.monthly=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %dopar%
    {
      zns=sort(unique(DATA.list.LIVEWT.c[[s]]$zone))
      if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #no enough observations 
      if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="Zone1")   #no enough observations 
      out=vector('list',length(zns))
      names(out)=zns
      for(z in 1:length(zns))
      {
        #1. Fit models
        DAT=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & zone==zns[z] &
                                               BLOCKX%in%BLKS.used[[s]])
        colnames(DAT)=tolower(colnames(DAT))
        #select years with a min number of vessels
        DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) 
        #Binomial
        DAT.bi=DAT%>%mutate(catch.target=ifelse(catch.target>0,1,0))
        Bi=fn.stand.delta(d=DAT.bi,Response="catch.target",PREDS=Predictors_monthly,
                          efrt="km.gillnet.hours.c",Formula=Best.Model_delta[[s]]$bi,
                          Formula.gam=NULL,Family="binomial")
        #Positive catch
        DAT=DAT%>%filter(catch.target>0)
        Pos=fn.stand.delta(d=DAT,Response="catch.target",PREDS=Predictors_monthly,
                           efrt="km.gillnet.hours.c",Formula=Best.Model_delta[[s]]$pos,
                           Formula.gam=NULL,Family="gaussian")
        
        #2.Make predictions
        Preds=fn.MC.delta.cpue(BiMOD=Bi$res,
                               MOD=Pos$res,
                               BiData=Bi$DATA%>%mutate(LNeffort=LN.effort),
                               PosData=Pos$DATA,
                               niter=Niter,
                               pred.term='finyear',
                               ALL.terms=Predictors_monthly)
        #3.Calculate CV
        Preds=Preds%>%mutate(CV=SD/response)
        
        #4.Apply creep
        Preds.creep=Preds
        yrs=unique(DAT$finyear)
        Preds.creep=subset(Preds.creep,finyear%in%yrs)
        add.crp=Eff.creep$effort.creep[match(Preds.creep$finyear,Eff.creep$finyear)]
        Preds.creep=Preds.creep%>%
          mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
                 upper.CL=upper.CL-(response-response*(1-add.crp)),
                 response=response*(1-add.crp))
        
        #5.Normalised
        Preds.nrmlzd=Preds.creep
        Mn=mean(Preds.nrmlzd$response)
        Preds.nrmlzd=Preds.nrmlzd%>%
          mutate(response=response/Mn,
                 lower.CL=lower.CL/Mn,
                 upper.CL=upper.CL/Mn) 
        out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
        rm(DAT)
      }
      return(out)
    }
  })   
  names(Zone_preds.monthly)=names(SP.list)[Tar.sp]
  
  #daily
  system.time({Zone_preds.daily=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm','mgcv')) %dopar%
    {
      zns=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$zone))
      if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #no enough observations 
      if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="Zone1")   #no enough observations 
      out=vector('list',length(zns))
      names(out)=zns
      for(z in 1:length(zns))
      {
        #1. Fit models
        DAT=DATA.list.LIVEWT.c.daily[[s]]%>%filter(VESSEL%in%VES.used.daily[[s]] & zone==zns[z]
                                                   & BLOCKX%in%BLKS.used.daily[[s]])
        colnames(DAT)=tolower(colnames(DAT))
        #select years with a min number of vessels
        DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT))
        DAT=DAT%>% mutate(nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                          mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
        
        ##
        #Binomial
        DAT.bi=DAT%>%mutate(catch.target=ifelse(catch.target>0,1,0))
        Bi=fn.stand.delta(d=DAT.bi,Response="catch.target",PREDS=Predictors_daily,
                          efrt="km.gillnet.hours.c",Formula=NULL,
                          Formula.gam=Best.Model.daily.gam_delta[[s]]$bi,Family="binomial")
        #Positive catch
        DAT=DAT%>%filter(catch.target>0)
        Pos=fn.stand.delta(d=DAT,Response="catch.target",PREDS=Predictors_daily,
                           efrt="km.gillnet.hours.c",Formula=NULL,
                           Formula.gam=Best.Model.daily.gam_delta[[s]]$pos,Family="gaussian")
        
        #2.Make predictions
        Preds=fn.MC.delta.cpue(BiMOD=Bi$res.gam,
                               MOD=Pos$res.gam,
                               BiData=Bi$DATA%>%mutate(LNeffort=LN.effort),
                               PosData=Pos$DATA,
                               niter=Niter,
                               pred.term='finyear',
                               ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner'))
        
        #3.Calculate CV
        Preds=Preds%>%mutate(CV=SD/response)
        
        #4.Apply creep
        Preds.creep=Preds
        yrs=unique(DAT$finyear)
        Preds.creep=subset(Preds.creep,finyear%in%yrs)
        add.crp=Eff.creep$effort.creep[match(Preds.creep$finyear,Eff.creep$finyear)]
        Preds.creep=Preds.creep%>%
          mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
                 upper.CL=upper.CL-(response-response*(1-add.crp)),
                 response=response*(1-add.crp))
        
        #5.Normalised
        Preds.nrmlzd=Preds.creep
        Mn=mean(Preds.nrmlzd$response)
        Preds.nrmlzd=Preds.nrmlzd%>%
          mutate(response=response/Mn,
                 lower.CL=lower.CL/Mn,
                 upper.CL=upper.CL/Mn) 
        out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
        rm(DAT)
      }
      return(out)
    }
  })   
  names(Zone_preds.daily)=names(SP.list)[Tar.sp]
  stopCluster(cl)
  
}})
}
if(Use.Tweedie)  
{
  #1. Get index
    #1.1 monthly  3 secs per species
  tic()
  Zone_preds.monthly=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
    zns=sort(unique(DATA.list.LIVEWT.c[[s]]$zone))
    if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #outside core area so not enough observations, issues with deg. freedom 
    if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")    
    out=vector('list',length(zns))
    names(out)=zns
    for(z in 1:length(zns))
    {
      d=DATA.list.LIVEWT.c[[s]]%>%
        filter(VESSEL%in%VES.used.zone[[s]][[z]] & BLOCKX%in%BLKS.used.zone[[s]][[z]] & zone==zns[z])%>%
        mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
      
      #remove first years with very few positive catch observation
      if(Exclude.yr.gummy) if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
      if(Exclude.yr.sandbar) if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1985-86','1986-87','1987-88','1988-89'))
      if(names(DATA.list.LIVEWT.c)[s]=="Whiskery Shark" & zns[z]=="West")d=d%>%filter(!FINYEAR%in%c('1975-76'))
      
      
      Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
      Continuous=Covariates.monthly
      colnames(d)=tolower(colnames(d))
      Terms=tolower(Terms)
      Continuous=tolower(Continuous)
      Factors=Terms[!Terms%in%Continuous]
      Terms=all.vars(Best.Model.zone[[s]][[z]])[-1]
      d <- d%>%
        dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
        mutate(cpue=catch.target/km.gillnet.hours.c)
      d <- makecategorical(Factors[Factors%in%Terms],d)
      
      #1. Fit model
      mod<-bam(Best.Model.zone[[s]][[z]],data=d,family='tw',method="fREML",discrete=TRUE)
      
      
      #2.Predict year effect (considering log bias corr if required)
      Preds=pred.fun(mod=mod,biascor="NO",PRED="finyear")
      
      #3.Calculate CV
      Preds=Preds%>%mutate(CV=SD/response)
      
      #4.Apply creep
      Preds.creep=Preds
      yrs=unique(d$finyear)
      Preds.creep=subset(Preds.creep,finyear%in%yrs)
      add.crp=Eff.creep$effort.creep[match(Preds.creep$finyear,Eff.creep$finyear)]
      Preds.creep=Preds.creep%>%
        mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
               upper.CL=upper.CL-(response-response*(1-add.crp)),
               response=response*(1-add.crp))
      
      #5.Normalised
      Preds.nrmlzd=Preds
      Mn=mean(Preds.nrmlzd$response)
      Preds.nrmlzd=Preds.nrmlzd%>%
        mutate(response=response/Mn,
               lower.CL=lower.CL/Mn,
               upper.CL=upper.CL/Mn)%>%
        mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
               upper.CL=upper.CL-(response-response*(1-add.crp)),
               response=response*(1-add.crp)) 
      
      out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
      rm(d,mod)
    }
    return(out)  
  }
  toc()   
  names(Zone_preds.monthly)=names(SP.list)[Tar.sp]
  
    #1.2 daily  4 secs per species
  tic()
  Zone_preds.daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
    zns=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$zone))
    if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #not enough observations, issues with deg. freedom 
    if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")    
    
    out=vector('list',length(zns))
    names(out)=zns
    
    for(z in 1:length(zns))
    {
      d=DATA.list.LIVEWT.c.daily[[s]]%>%
        filter(VESSEL%in%VES.used.zone.daily[[s]][[z]] & BLOCKX%in%BLKS.used.zone.daily[[s]][[z]] & zone==zns[z])%>%
        mutate(shots.c=ifelse(shots.c>2,2,shots.c),
               nlines.c=ifelse(nlines.c>3,3,nlines.c),
               mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
      
      #remove first year of transition from Monthly to Daily returns due to effort misreporting
      d=d%>%filter(!FINYEAR=="2006-07")
      
      Terms=Predictors_daily
      Continuous=Covariates.daily
      colnames(d)=tolower(colnames(d))
      Terms=tolower(Terms)
      Continuous=tolower(Continuous)
      Factors=Terms[!Terms%in%Continuous]
      Terms=all.vars(Best.Model.zone.daily[[s]][[z]])[-1]
      d <- d%>%
        dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
        mutate(cpue=catch.target/km.gillnet.hours.c)
      d <- makecategorical(Factors[Factors%in%Terms],d)
      
      #1. Fit model
      mod<-bam(Best.Model.zone.daily[[s]][[z]],data=d,family='tw',method="fREML",discrete=TRUE) 
      
      #2.Predict year effect (considering log bias corr if required)
      Preds=pred.fun(mod=mod,biascor="NO",PRED="finyear")
      
      #3.Calculate CV
      Preds=Preds%>%mutate(CV=SD/response)
      
      #4.Apply creep
      Preds.creep=Preds
      yrs=unique(d$finyear)
      Preds.creep=subset(Preds.creep,finyear%in%yrs)
      add.crp=Eff.creep$effort.creep[match(Preds.creep$finyear,Eff.creep$finyear)]
      Preds.creep=Preds.creep%>%
        mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
               upper.CL=upper.CL-(response-response*(1-add.crp)),
               response=response*(1-add.crp))
      
      #5.Normalised
      Preds.nrmlzd=Preds
      Mn=mean(Preds.nrmlzd$response)
      Preds.nrmlzd=Preds.nrmlzd%>%
        mutate(response=response/Mn,
               lower.CL=lower.CL/Mn,
               upper.CL=upper.CL/Mn)%>%
        mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
               upper.CL=upper.CL-(response-response*(1-add.crp)),
               response=response*(1-add.crp)) 

      out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
      
      rm(d,mod)
    }
    
    return(out)
  }
  toc()
  names(Zone_preds.daily)=names(SP.list)[Tar.sp]
  
  
  #2. Plot index
  fn.plot.spatial.indices=function(sp,Pred1,Pred2,Pred.spatial,Pred1.daily,Pred2.daily,Pred.spatial.daily) 
  {
    #Monthly
    d5=Pred1[[match(sp, names(Pred1))]]%>%
      mutate(method='Single zone',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='1.monthly')
    
    mn=Pred1[[match(sp, names(Pred1))]]$response
    d6=Pred2[[match(sp, names(Pred2))]]%>%
      mutate(method='Single zone with creep',
             lower.CL=lower.CL/mean(mn),
             upper.CL=upper.CL/mean(mn),
             mean=response/mean(mn),
             year=0.05+as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='1.monthly')
    
    dd=Pred.spatial[[match(sp, names(Pred.spatial))]]
    d7=list.map(dd,Preds)
    for(x in 1:length(d7))
    {
      d7[[x]]=d7[[x]]%>%
              mutate(method=names(d7)[x],
                     lower.CL=lower.CL/mean(response),
                     upper.CL=upper.CL/mean(response),
                     mean=response/mean(response),
                     year=-0.05+as.numeric(substr(finyear,1,4)))%>%
              dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
              mutate(period='1.monthly')
    }
    d7=do.call(rbind,d7)
    
    d8=list.map(dd,Preds.nrmlzd)
    for(x in 1:length(d8))
    {
      d8[[x]]=d8[[x]]%>%
        mutate(method=paste(names(d8)[x],'with creep'),
               mean=response,
               year=-0.15+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='1.monthly')
    }
    d8=do.call(rbind,d8)
    
    
    #Daily
    d5.daily=Pred1.daily[[match(sp, names(Pred1.daily))]]%>%
      mutate(method='Single zone',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='2.daily')
    
    d6.daily=Pred2.daily[[match(sp, names(Pred2.daily))]]%>%
      mutate(method='Single zone with creep',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=0.05+as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='2.daily')
    
    
    dd=Pred.spatial.daily[[match(sp, names(Pred.spatial.daily))]]
    d7.daily=list.map(dd,Preds)
    for(x in 1:length(d7.daily))
    {
      d7.daily[[x]]=d7.daily[[x]]%>%
        mutate(method=names(d7.daily)[x],
               lower.CL=lower.CL/mean(response),
               upper.CL=upper.CL/mean(response),
               mean=response/mean(response),
               year=-0.05+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='2.daily')
    }
    d7.daily=do.call(rbind,d7.daily)
    
    d8.daily=list.map(dd,Preds.nrmlzd)
    for(x in 1:length(d8.daily))
    {
      d8.daily[[x]]=d8.daily[[x]]%>%
        mutate(method=paste(names(d8.daily)[x],'with creep'),
               mean=response,
               year=-0.15+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='2.daily')
    }
    d8.daily=do.call(rbind,d8.daily)
    
    
    
    
    #plots
    LVLS=c("Single zone","Single zone with creep",
           "West","West with creep",
           "Zone1","Zone1 with creep",
           "Zone2", "Zone2 with creep")
    Kls=c('grey50','grey10','chartreuse3','forestgreen','firebrick4','brown1','cyan3','blue4')
    
    
    p1=rbind(d5,d6,d7,d8)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Monthly')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())+
      scale_color_manual(values=Kls)
    
    p2=rbind(d5.daily,d6.daily,d7.daily,d8.daily)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Daily')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())+
      scale_color_manual(values=Kls)
    
    p=ggarrange(p1,p2,ncol=1,common.legend = T)
    annotate_figure(p,
                    left = textGrob("Relative CPUE (+/- 95% CI)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                    bottom = textGrob("Financial year", gp = gpar(cex = 1.3)))
    return(p)
    
  }
  for(d in 1:length(Zone_preds.monthly))
  {
    NM=names(Zone_preds.monthly)[d]
    fn.plot.spatial.indices(sp=NM,
                            Pred1=Pred,
                            Pred2=Pred.creep,
                            Pred.spatial=Zone_preds.monthly,
                            Pred1.daily=Pred.daily,
                            Pred2.daily=Pred.daily.creep,
                            Pred.spatial.daily=Zone_preds.daily)
    ggsave(handl_OneDrive(paste('Analyses/Catch and effort/Outputs/Compare all index types/',
                                NM,'_Spatial.jpeg',sep='')),width=10,height= 10)
  }
}

#-- Compare approaches for deriving spatial index  ----------------------------------------------
compare.spatial.approach=FALSE  
if(compare.spatial.approach)
{
  #Single model with zone interaction
    #Monthly
  Best.Model.single.model=Best.Model
  Best.Model.single.model$`Gummy Shark`=formula(cpue ~ blockx + s(vessel, bs = "re") + znyr)
  Best.Model.single.model$`Whiskery Shark`=formula(cpue ~ blockx + s(vessel, bs = "re") + znyr)
  Best.Model.single.model$`Dusky Whaler`=formula(cpue ~ blockx + s(vessel, bs = "re") + znyr)
  Best.Model.single.model$`Sandbar Shark`=formula(cpue ~  blockx + s(vessel, bs = "re") + s(month, k = 12, bs = "cc")+znyr)
  Zone_preds.monthly_single_model=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
       d=DATA.list.LIVEWT.c[[s]]%>%
        filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])%>%
        mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
      
       if(names(SP.list)[s]=="Sandbar Shark") d=subset(d,!zone=="Zone2")   #not enough observations, issues with deg. freedom 
       if(names(SP.list)[s]=="Gummy Shark") d=subset(d,!zone=="West")    
       
       
      #remove first years with very few positive catch observation
      if(Exclude.yr.gummy) if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
      if(Exclude.yr.sandbar) if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1985-86','1986-87','1987-88','1988-89'))
      
      d=d%>%mutate(znyr = interaction(zone, FINYEAR, drop = TRUE))
      
      Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
      Continuous=Covariates.monthly
      colnames(d)=tolower(colnames(d))
      Terms=tolower(Terms)
      Continuous=tolower(Continuous)
      Factors=c("zone",Terms[!Terms%in%Continuous])
      Terms=all.vars(Best.Model.single.model[[s]])[-1]
      d <- d%>%
        mutate(year=as.integer(substr(finyear,1,4)))%>%
        dplyr::select(c(year,catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
        mutate(cpue=catch.target/km.gillnet.hours.c)
      d <- makecategorical(Factors[Factors%in%Terms],d)
      
      #1. Fit model
      mod<-bam(Best.Model.single.model[[s]],data=d,family='tw',method="fREML",discrete=TRUE)

      #2.Predict year effect (considering log bias corr if required)
      Preds=summary(emmeans(mod, 'znyr', type="response", rg.limit = 2e5,data=d))
      #Preds=summary(emmeans(mod, ~ finyear | zone, type="response", rg.limit = 2e5,data=d))
      Preds=Preds%>%
        mutate(zone=sub("\\..*", "", znyr),
               finyear=sub(".*\\.", "", znyr))
   
      return(Preds)
    }
  names(Zone_preds.monthly_single_model)=names(SP.list)[Tar.sp]
  
    #Daily
  Best.Model.daily.single.model=Best.Model.daily
  Best.Model.daily.single.model$`Gummy Shark`=formula(cpue ~  s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                        block10 + step.mcal_target_group + 
                                                        mesh + znyr)
  Best.Model.daily.single.model$`Whiskery Shark`=formula(cpue ~  s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                           block10 + s(mean.depth) + step.mcal_target_group + 
                                                           shots.c + mesh + znyr)
  Best.Model.daily.single.model$`Dusky Whaler`=formula(cpue ~ s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                         block10 + s(mean.depth) + step.mcal_target_group + 
                                                          shots.c + mesh+ znyr)
  Best.Model.daily.single.model$`Sandbar Shark`=formula(cpue ~ s(vessel, bs = "re") + s(month, bs = "cc") + 
                                                          block10 + s(mean.depth) + step.mcal_target_group+ znyr)
  
  Zone_preds.daily_single_model=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
    d=DATA.list.LIVEWT.c.daily[[s]]%>%
      filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])%>%
      mutate(shots.c=ifelse(shots.c>2,2,shots.c),
             nlines.c=ifelse(nlines.c>3,3,nlines.c),
             mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    
    if(names(SP.list)[s]=="Sandbar Shark") d=subset(d,!zone=="Zone2")   #not enough observations, issues with deg. freedom 
    if(names(SP.list)[s]=="Gummy Shark") d=subset(d,!zone=="West")    
    
    #remove first year of transition from Monthly to Daily returns due to effort misreporting
    d=d%>%filter(!FINYEAR=="2006-07")
    d=d%>%mutate(znyr = interaction(zone, FINYEAR, drop = TRUE))
    Terms=Predictors_daily
    Continuous=Covariates.daily
    colnames(d)=tolower(colnames(d))
    Terms=tolower(Terms)
    Continuous=tolower(Continuous)
    Factors=c(Terms[!Terms%in%Continuous],'zone')
    Terms=c(all.vars(Best.Model.daily.single.model[[s]])[-1],'zone')
    d <- d%>%
      dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
      mutate(cpue=catch.target/km.gillnet.hours.c)
    d <- makecategorical(Factors[Factors%in%Terms],d)
    
    #1. Fit model
    mod<-bam(Best.Model.daily.single.model[[s]],data=d,family='tw',method="fREML",discrete=TRUE) 
    
    #2.Predict year effect (considering log bias corr if required)
    Preds=summary(emmeans(mod, 'znyr', type="response", rg.limit = 2e5,data=d))
    Preds=Preds%>%
      mutate(zone=sub("\\..*", "", znyr),
             finyear=sub(".*\\.", "", znyr))
    return(Preds)
    }
  names(Zone_preds.daily_single_model)=names(SP.list)[Tar.sp]
  
  #Compare 2 approaches for spatial cpue
  fn.plot.spatial.indices2=function(sp,Pred,Pred.spatial,Pred.spatial_single,
                                    Pred.daily,Pred.spatial.daily,Pred.spatial.daily_single) 
  {
    #Monthly
    d5=Pred[[match(sp, names(Pred))]]%>%
      mutate(method='Single zone',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='1.monthly')
    
    dd=Pred.spatial[[match(sp, names(Pred.spatial))]]
    d7=list.map(dd,Preds)
    for(x in 1:length(d7))
    {
      d7[[x]]=d7[[x]]%>%
        mutate(method=names(d7)[x],
               lower.CL=lower.CL/mean(response),
               upper.CL=upper.CL/mean(response),
               mean=response/mean(response),
               year=-0.1+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='1.monthly')
    }
    d7=do.call(rbind,d7)
    
    d8=Pred.spatial_single[[match(sp, names(Pred.spatial_single))]]%>%
      mutate(method=paste0(zone,'_single modl'),
             period='1.monthly',
             year=0.1+as.numeric(substr(finyear,1,4)))%>%
      group_by(method,period)%>%
      mutate(Mn=mean(response))%>%
      mutate(lower.CL=lower.CL/Mn,
             upper.CL=upper.CL/Mn,
             mean=response/Mn)%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL,period)
    
    
    #Daily
    d5.daily=Pred.daily[[match(sp, names(Pred.daily))]]%>%
      mutate(method='Single zone',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='2.daily')
    
    dd=Pred.spatial.daily[[match(sp, names(Pred.spatial.daily))]]
    d7.daily=list.map(dd,Preds)
    for(x in 1:length(d7.daily))
    {
      d7.daily[[x]]=d7.daily[[x]]%>%
        mutate(method=names(d7.daily)[x],
               lower.CL=lower.CL/mean(response),
               upper.CL=upper.CL/mean(response),
               mean=response/mean(response),
               year=-0.1+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='2.daily')
    }
    d7.daily=do.call(rbind,d7.daily)
    
    d8.daily=Pred.spatial.daily_single[[match(sp, names(Pred.spatial.daily_single))]]%>%
      mutate(method=paste0(zone,'_single modl'),
             period='2.daily',
             year=0.1+as.numeric(substr(finyear,1,4)))%>%
      group_by(method,period)%>%
      mutate(Mn=mean(response))%>%
      mutate(lower.CL=lower.CL/Mn,
             upper.CL=upper.CL/Mn,
             mean=response/Mn)%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL,period)
    
    
    #plots
    LVLS=c("Single zone",
           "West","West_single modl","Zone1","Zone1_single modl","Zone2",
           "Zone2_single modl")
    Kls=c('grey50','chartreuse3','forestgreen','firebrick4','brown1','cyan3','blue4')
    names(Kls)=LVLS  
    p1=rbind(d5,d7,d8)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Monthly')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())+
      scale_color_manual(values=Kls)
    
    p2=rbind(d5.daily,d7.daily,d8.daily)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Daily')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())+
      scale_color_manual(values=Kls)
    
    p=ggarrange(p1,p2,ncol=1,common.legend = T)
    annotate_figure(p,
                    left = textGrob("Relative CPUE (+/- 95% CI)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                    bottom = textGrob("Financial year", gp = gpar(cex = 1.3)))
    return(p)
    
  }
  for(d in 1:length(Zone_preds.monthly))
  {
    NM=names(Zone_preds.monthly)[d]
    fn.plot.spatial.indices2(sp=NM,
                             Pred=Pred,
                             Pred.spatial=Zone_preds.monthly,
                             Pred.spatial_single=Zone_preds.monthly_single_model,
                             Pred.daily=Pred.daily,
                             Pred.spatial.daily=Zone_preds.daily,
                             Pred.spatial.daily_single=Zone_preds.daily_single_model)
    ggsave(handl_OneDrive(paste('Analyses/Catch and effort/Outputs/Compare all index types/',
                                NM,'_Spatial_two approaches.jpeg',sep='')),width=10,height= 10)
  }
  
}
 
