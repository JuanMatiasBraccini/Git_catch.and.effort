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
  #Get index
    #monthly  3 secs per species
  tic()
  Zone_preds.monthly=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
      if(!is.null(BLKS.used[[s]]))
      {
        zns=sort(unique(DATA.list.LIVEWT.c[[s]]$zone))
        if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #not enough observations, issues with deg. freedom 
        if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")    
        out=vector('list',length(zns))
        names(out)=zns
        for(z in 1:length(zns))
        {
          d=DATA.list.LIVEWT.c[[s]]%>%
            filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]] & zone==zns[z])%>%
            mutate(SHOTS.c=ifelse(SHOTS.c>2,2,SHOTS.c))
          
          #remove first years with very few positive catch observation
          if(Exclude.yr.gummy) if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
          if(Exclude.yr.sandbar) if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
          
          Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
          Continuous=Covariates.monthly
          colnames(d)=tolower(colnames(d))
          Terms=tolower(Terms)
          Continuous=tolower(Continuous)
          Factors=Terms[!Terms%in%Continuous]
          Terms=all.vars(Best.Model[[s]])[-1]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          d <- makecategorical(Factors[Factors%in%Terms],d)
          
          #1. Fit model
          mod<-bam(Best.Model[[s]],data=d,family='tw',method="fREML",discrete=TRUE)
          
          
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
          Preds.nrmlzd=Preds.creep
          Mn=mean(Preds.nrmlzd$response)
          Preds.nrmlzd=Preds.nrmlzd%>%
            mutate(response=response/Mn,
                   lower.CL=lower.CL/Mn,
                   upper.CL=upper.CL/Mn) 
          
          out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
          rm(d,mod)
        }
        return(out)
      }
    }
  toc()   
  names(Zone_preds.monthly)=names(SP.list)[Tar.sp]
  
    #daily  4 secs per species
  tic()
  Zone_preds.daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','emmeans')) %dopar%
  {
      if(!is.null(BLKS.used.daily[[s]]))
      {
        zns=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$zone))
        if(names(SP.list)[s]=="Sandbar Shark") zns=subset(zns,!zns=="Zone2")   #not enough observations, issues with deg. freedom 
        if(names(SP.list)[s]=="Gummy Shark") zns=subset(zns,!zns=="West")    
        
        out=vector('list',length(zns))
        names(out)=zns
        
        for(z in 1:length(zns))
        {
          d=DATA.list.LIVEWT.c.daily[[s]]%>%
            filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]] & zone==zns[z])%>%
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
          Terms=all.vars(Best.Model.daily[[s]])[-1]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,all_of(Terms)))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          d <- makecategorical(Factors[Factors%in%Terms],d)
          
          #1. Fit model
          mod<-bam(Best.Model.daily[[s]],data=d,family='tw',method="fREML",discrete=TRUE) 
          
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
          Preds.nrmlzd=Preds.creep
          Mn=mean(Preds.nrmlzd$response)
          Preds.nrmlzd=Preds.nrmlzd%>%
            mutate(response=response/Mn,
                   lower.CL=lower.CL/Mn,
                   upper.CL=upper.CL/Mn) 
          
          out[[z]]=list(Preds=Preds,Preds.creep=Preds.creep,Preds.nrmlzd=Preds.nrmlzd)
          
          rm(d,mod)
        }
        
        return(out)
      }
    }
  toc()
  names(Zone_preds.daily)=names(SP.list)[Tar.sp]
  
  
  
  
  #Plot index
  fn.plot.spatial.indices=function(sp,Pred1,Pred2,Pred.spatial,Pred1.daily,Pred2.daily,Pred.spatial.daily) 
  {
    #Monthly
    d5=Pred1[[match(sp, names(Pred1))]]%>%
      mutate(method='Standardised',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=0.15+as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='1.monthly')
    
    mn=Pred1[[match(sp, names(Pred1))]]$response
    d6=Pred2[[match(sp, names(Pred2))]]%>%
      mutate(method='Standardised with creep',
             lower.CL=lower.CL/mean(mn),
             upper.CL=upper.CL/mean(mn),
             mean=response/mean(mn),
             year=0.1+as.numeric(substr(finyear,1,4)))%>%
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
    
    d8=list.map(dd,Preds.nrmlzd)
    for(x in 1:length(d8))
    {
      d8[[x]]=d8[[x]]%>%
        mutate(method=paste(names(d8)[x],'with creep'),
               lower.CL=lower.CL/mean(response),
               upper.CL=upper.CL/mean(response),
               mean=response/mean(response),
               year=-0.15+as.numeric(substr(finyear,1,4)))%>%
        dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
        mutate(period='1.monthly')
    }
    d8=do.call(rbind,d8)
    
    
    #Daily
    d5.daily=Pred1.daily[[match(sp, names(Pred1.daily))]]%>%
      mutate(method='Standardised',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=0.15+as.numeric(substr(finyear,1,4)))%>%
      dplyr::select(year,method,mean,lower.CL,upper.CL)%>%
      mutate(period='2.daily')
    
    d6.daily=Pred2.daily[[match(sp, names(Pred2.daily))]]%>%
      mutate(method='Standardised with creep',
             lower.CL=lower.CL/mean(response),
             upper.CL=upper.CL/mean(response),
             mean=response/mean(response),
             year=0.1+as.numeric(substr(finyear,1,4)))%>%
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
    LVLS=c("Standardised","Standardised with creep",
           "West","West with creep",
           "Zone1","Zone1 with creep",
           "Zone2","Zone2 with creep")
    p1=rbind(d5,d6,d7,d8)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Monthly')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())
    
    p2=rbind(d5.daily,d6.daily,d7.daily,d8.daily)%>%
      mutate(method=factor(method,levels=LVLS))%>%
      ggplot(aes(year,mean,color=method))+
      geom_line(alpha=0.35,linetype='dashed')+
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,alpha=0.35)+
      geom_point(size=2.5)+ylim(0,NA)+
      ggtitle('Daily')+theme_PA()+ylab('')+xlab('')+
      theme(legend.position = 'top',legend.title = element_blank())
    
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
      

    ggsave(handl_OneDrive(paste('Analyses/Catch and effort/Outputs/Compare all index types/Spatial/',
                                NM,'.jpeg',sep='')),width=10,height= 10)  
    
  }
}



 
