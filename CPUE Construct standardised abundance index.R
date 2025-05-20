ZONES=c("West","Zone1","Zone2")
Response="catch.target"    #cpue and positive catch are calculated inside functions
Eff.vars=c("km.gillnet.hours.c")
Categorical=c("finyear","vessel","blockx","shots.c")
Covariates.monthly=c("MONTH","LONG10.corner","LAT10.corner","Freo","Temp.res","SOI","Yrs.of.experience")
if(do_Stephens_McCall=="YES") Targeting.vars=c("Step.MCal_target_group") #"Step.MCal_target_prob
if(do_cluster=="YES") Targeting.vars=c("cluster_clara")
if(do_pca=="YES") Targeting.vars=c("Dim.1","Dim.2","Dim.3")
Covariates.daily=c(Covariates.monthly,"Mean.depth","Lunar","mesh","nlines.c")
if(any(do_pca=="YES"|do_cluster=="YES"|do_Stephens_McCall=="YES"))Covariates.daily=c(Covariates.daily,Targeting.vars)

Predictors_monthly=c(Categorical,Covariates.monthly)
Predictors_daily=c(Categorical,Covariates.daily)


#-- check data properties and degrees of freedom ----------------------------------------------
if(Model.run=="First")
{
  Prop.deg.free.m=Prop.deg.free.d=Prop.deg.free.d_n=Store_nom_cpues_monthly
  for(s in Tar.sp)
  {
    dummy=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used[[s]])
    Prop.m=properties(dummy)
    d.f=Prop.m[match(Predictors_monthly,rownames(Prop.m)),]
    d.f$dummy=as.numeric(with(d.f,ifelse(Class%in%c('numeric'),1,Unique)))
    d.f.m=data.frame(Deg.F=sum(d.f$dummy),Obser=nrow(dummy))
    Prop.deg.free.m[[s]]=list(dat=d.f,deg.f=d.f.m)
    
    dummy=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    dummy=subset(dummy,blockx%in%BLKS.used.daily[[s]])
    Prop.d=properties(dummy)
    d.f=Prop.d[match(Predictors_daily,rownames(Prop.d)),]
    d.f$dummy=as.numeric(with(d.f,ifelse(Class%in%c('numeric'),1,Unique)))
    d.f.d=data.frame(Deg.F=sum(d.f$dummy),Obser=nrow(dummy))
    Prop.deg.free.d[[s]]=list(dat=d.f,deg.f=d.f.d)
  }
}

#-- Show applied effort creep ----------------------------------------------
Fish.pow.inc=FINYEAR.ALL
names(Fish.pow.inc)=FINYEAR.ALL
Fish.pow.inc=as.numeric(substr(Fish.pow.inc,1,4))
Fish.pow.inc=c(seq(0,0.4,by=Fish.Pow),
               rep(0.4,length(1996:Fish.pow.inc[length(Fish.pow.inc)])))
if(Model.run=="First")
{
  fn.fig("Effort_creep/Effort_creep_applied",2400,2400)
  par(las=1,mgp=c(2.5,.9,0))
  plot(1:length(Fish.pow.inc),Fish.pow.inc*100,xaxt='n',ylab="Effort creep (% increase)",
       xlab="Finacial year",type='o',pch=19,cex.lab=1.6,cex.axis=1.35,cex=1.5)
  axis(1,1:length(Fish.pow.inc),F,tck=-0.01)
  axis(1,seq(1,length(Fish.pow.inc),10),FINYEAR.ALL[seq(1,length(Fish.pow.inc),10)],tck=-0.04,cex.axis=1.35)
  segments(6,Fish.pow.inc[6]*100,6,Fish.pow.inc[5]*100,col=2,lwd=2)
  segments(5,Fish.pow.inc[5]*100,6,Fish.pow.inc[5]*100,col=2,lwd=2)
  text(9,Fish.pow.inc[5]*100,"2 %",col=2,cex=2)
  dev.off()
}
Eff.creep=data.frame(finyear=FINYEAR.ALL,effort.creep=Fish.pow.inc)


#-- Select model structure ----------------------------------------------
add.CITES.term=TRUE  #CITES listing term for smooth HH
#remove predictors identified as highly correlated
cNSTNT=Categorical[!Categorical=="block10"] 
cNSTNT.daily=Categorical[!Categorical=="blockx"]
if(do_cluster=="YES") cNSTNT.daily=c(Categorical[!Categorical=="blockx"],Targeting.vars)

#extract best model   
Best.Model=vector('list',length(SP.list)) 
names(Best.Model)=names(SP.list)
Best.Model.daily=Best.Model.daily.gam=Best.Model_delta=Best.Model.daily.gam_delta=Best.Model.daily_delta=
  Store.Best.Model=Store.Best.Model.daily=Best.Model
if(Def.mod.Str=="YES")     
{
  hndl.modl.sel=handl_OneDrive("Analyses/Catch and effort/Outputs/Model Selection/")
  if(Use.Tweedie)       #takes 3.5 hours
  {
    #target species   #takes 45 mins per species
    system.time({foreach(s=Tar.sp,.packages=c('dplyr','mgcv','doParallel')) %do%
        {
          #Monthly
          do.monthly=FALSE
          if(do.monthly)
          {
            d=DATA.list.LIVEWT.c[[s]]%>%
              filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])
            if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
            if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
            Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
            Continuous=Covariates.monthly
            Fmula=list(Year.block="cpue~finyear+blockx",
                       Shots="cpue~finyear+blockx+shots.c",
                       Month="cpue~finyear+blockx+shots.c+s(month,k=12,bs='cc')",
                       Temp.res="cpue~finyear+blockx+shots.c+s(month,k=12,bs='cc')+s(temp.res)",
                       Vessel="cpue~finyear+blockx+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(vessel,bs='re')")
            Family='tw'
            FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],"_GAM_monthly",sep="")
            colnames(d)=tolower(colnames(d))
            Terms=tolower(Terms)
            Continuous=tolower(Continuous)
            Factors=Terms[!Terms%in%Continuous]
            d=d%>%
              dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
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
            Explained.dev=mod.sumery
            for(x in 1:length(MODS))
            {
              plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,pages=1)
              mod.sumery[[x]]=tidy(MODS[[x]])
              if(length(mod.sumery[[x]])>0)
              {
                grid.newpage()
                class(mod.sumery[[x]]) <- "data.frame"
                grid.table(mod.sumery[[x]])
              }
              Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100  
            }
            Explained.dev=data.frame(Model=names(Fmula),
                                     Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
              mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
            grid.newpage()
            class(Explained.dev) <- "data.frame"
            grid.table(Explained.dev)
            dev.off()
            
          }
          
          
          #Daily
          d=DATA.list.LIVEWT.c.daily[[s]]%>%
            filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
          Terms=Predictors_daily
          Continuous=Covariates.daily
          Fmula=list(Year=    "cpue~finyear",
                     Shots=   "cpue~finyear+shots.c",
                     Month=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')",
                     Temp.res="cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)",
                     Depth=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)",
                     Lunar=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)",
                     Lat.long="cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(long10.corner,lat10.corner)",
                     Vessel=  "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(long10.corner,lat10.corner)+s(vessel,bs='re')"
                     #PCA.1=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(dim.1)",
                     #PCA.2=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(dim.1)+s(dim.2)",
                     #PCA.3=   "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(dim.1)+s(dim.2)+s(dim.3)",
                     #Lat.long="cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(dim.1)+s(dim.2)+s(dim.3)+s(long10.corner,lat10.corner)",
                     #Vessel=  "cpue~finyear+shots.c+s(month,k=12,bs='cc')+s(temp.res)+s(mean.depth)+s(lunar)+s(dim.1)+s(dim.2)+s(dim.3)+s(long10.corner,lat10.corner)+s(vessel,bs='re')"
          )
          Family='tw'
          FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],"_GAM_daily",sep="")
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
          Explained.dev=mod.sumery
          for(x in 1:length(MODS))
          {
            plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,pages=1)
            mod.sumery[[x]]=tidy(MODS[[x]])
            if(length(mod.sumery[[x]])>0)
            {
              grid.newpage()
              class(mod.sumery[[x]]) <- "data.frame"
              grid.table(mod.sumery[[x]])
            }
            Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100  
          }
          Explained.dev=data.frame(Model=names(Fmula),
                                   Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
            mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
          grid.newpage()
          class(Explained.dev) <- "data.frame"
          grid.table(Explained.dev)
          dev.off()
        }
    })
    
    #other species    #takes 2 mins per species
    system.time({foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mgcv','doParallel')) %do%
        {
          #Monthly
          if(!is.null(BLKS.used[[s]]))
          {
            d=DATA.list.LIVEWT.c[[s]]%>%
              filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])
            if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
            if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
            Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
            Continuous=Covariates.monthly
            Fmula=list(Year.block="cpue~finyear+blockx",
                       Month="cpue~finyear+blockx+s(month,k=12,bs='cc')",
                       Vessel="cpue~finyear+blockx+s(month,k=12,bs='cc')+s(vessel,bs='re')")
            Family='tw'
            FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],"_GAM_monthly",sep="")
            colnames(d)=tolower(colnames(d))
            Terms=tolower(Terms)
            Continuous=tolower(Continuous)
            Factors=Terms[!Terms%in%Continuous]
            d=d%>%
              dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
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
            Explained.dev=mod.sumery
            for(x in 1:length(MODS))
            {
              plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,pages=1)
              mod.sumery[[x]]=tidy(MODS[[x]])
              if(length(mod.sumery[[x]])>0)
              {
                grid.newpage()
                class(mod.sumery[[x]]) <- "data.frame"
                grid.table(mod.sumery[[x]])
              }
              Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100
            }
            Explained.dev=data.frame(Model=names(Fmula),
                                     Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
              mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
            grid.newpage()
            class(Explained.dev) <- "data.frame"
            grid.table(Explained.dev)
            dev.off()
          }
          
          #Daily
          if(!is.null(BLKS.used.daily[[s]]))
          {
            d=DATA.list.LIVEWT.c.daily[[s]]%>%
              filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
            Terms=Predictors_daily[!Predictors_daily%in%c("Dim.1","Dim.2","Dim.3")]
            Continuous=Covariates.daily[!Covariates.daily%in%c("Dim.1","Dim.2","Dim.3")]
            Fmula=list(Year=    "cpue~finyear",
                       Month=   "cpue~finyear+s(month,k=12,bs='cc')",
                       Depth=   "cpue~finyear+s(month,k=12,bs='cc')+s(mean.depth)",
                       Lat.long="cpue~finyear+s(month,k=12,bs='cc')+s(mean.depth)+s(long10.corner,lat10.corner)",
                       Vessel=  "cpue~finyear+s(month,k=12,bs='cc')+s(mean.depth)+s(long10.corner,lat10.corner)+s(vessel,bs='re')")
            Family='tw'
            FILE=paste(names(DATA.list.LIVEWT.c.daily)[s],"_GAM_daily",sep="")
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
            Explained.dev=mod.sumery
            for(x in 1:length(MODS))
            {
              plot(MODS[[x]],all.terms = TRUE, shade = TRUE,shade.col = "orange",scale=0,pages=1)
              mod.sumery[[x]]=tidy(MODS[[x]])
              if(length(mod.sumery[[x]])>0)
              {
                grid.newpage()
                class(mod.sumery[[x]]) <- "data.frame"
                grid.table(mod.sumery[[x]])
              }
              Explained.dev[[x]]=summary(MODS[[x]])$dev.expl*100
            }
            Explained.dev=data.frame(Model=names(Fmula),
                                     Percent.deviance.explained.by.model=do.call(rbind,Explained.dev))%>%
              mutate(Percent.explained.by.term=Percent.deviance.explained.by.model-lag(Percent.deviance.explained.by.model,1))
            grid.newpage()
            class(Explained.dev) <- "data.frame"
            grid.table(Explained.dev)
            dev.off()
          }
        }
    })
    
    
  }
  
  if(Use.Delta)
  {
    select.best.method="dev.exp"
    if(select.best.method=="dev.exp")     #takes  8 mins
    {
      fn.modl.sel.delta=function(K,RESPNS,dat,Family,MOD.type)   #function for delta model structure selection
      {
        PREDS.test=c(1,PREDS[-match(K,PREDS)])
        mod=vector('list',length(PREDS.test))
        names(mod)=PREDS.test
        if(MOD.type=="GAM") K=c('finyear+vessel+s(long10.corner,lat10.corner)+month')
        for(p in 1:length(PREDS.test))
        {
          if(MOD.type=="GLM")
          {
            if(RESPNS=="LNcpue")Formula=formula(paste("LNcpue",paste(c(K,PREDS.test[1:p]),collapse="+"),sep="~"))
            if(RESPNS=="catch.target")Formula=formula(paste(Response,paste(paste(c(K,PREDS.test[1:p]),collapse="+"),
                                                                           paste("offset(","LN.effort)",sep=""),sep="+"),sep="~"))
            mod[[p]]=glm(Formula,data=dat,family = Family, maxit=500)
          }
          if(MOD.type=="GAM")
          {
            add.this=PREDS.test[1:p]
            if(add.this[p]%in%Covariates) add.this[p]=paste("s(",add.this[p],",k=6)",sep='')
            if(RESPNS=="LNcpue")Formula=formula(paste("LNcpue",paste(c(K,add.this),collapse="+"),sep="~"))
            if(RESPNS=="catch.target")Formula=formula(paste(Response,paste(paste(c(K,add.this),collapse="+"),
                                                                           paste("offset(","LN.effort)",sep=""),sep="+"),sep="~"))
            mod[[p]]=gam(Formula,data=dat,family = Family,method="REML")
          }
        }
        return(mod)
      }
      efrt="km.gillnet.hours.c"
      system.time({Model.structure=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %dopar%
        {
          #need for only daily 
          d=DATA.list.LIVEWT.c.daily[[s]]%>%filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
          colnames(d)=tolower(colnames(d))
          PREDS=Predictors_daily[which(Predictors_daily%in%names(d))]
          d=d%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%BLKS.used.daily[[s]])
          d=d %>% mutate(nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                         mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
          id.cov=which(PREDS%in%Covariates)
          #binomial part
          d.bi=d%>%mutate(catch.target=ifelse(catch.target>0,1,0))
          d.bi=makecategorical(PREDS[which(PREDS%in%Categorical)],d.bi)
          d.bi$LN.effort=log(d.bi[,match(efrt,names(d.bi))])
          daily.mods.bi=fn.modl.sel.delta(K=cNSTNT.daily,RESPNS='catch.target',dat=d.bi,
                                          Family="binomial",MOD.type="GAM")
          #pos part
          d=d%>%filter(catch.target>0)
          d=makecategorical(PREDS[which(PREDS%in%Categorical)],d)
          d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
          daily.mods.pos=fn.modl.sel.delta(K=cNSTNT.daily,RESPNS='LNcpue',dat=d,
                                           Family="gaussian",MOD.type="GAM")
          
          return(list(daily.mods.bi=daily.mods.bi,daily.mods.pos=daily.mods.pos))
          
        }
      })
      stopCluster(cl)
      names(Model.structure)=names(SP.list)[Tar.sp]
      
      #calculate AIC and deviance explained per term
      Tab.monthly=vector('list',length(Tar.sp))
      Tab.daily=Tab.monthly
      for(i in 1:length(Tar.sp))
      {
        #Daily bi
        n=length(Model.structure[[i]]$daily.mods.bi)
        Aic.daily.bi=rep(NA,n)
        names(Aic.daily.bi)=names(Model.structure[[i]]$daily.mods.bi)
        Res.DEv.daily.bi=Aic.daily.bi
        for( l in 1:n)
        {
          Aic.daily.bi[l]=aicc(Model.structure[[i]]$daily.mods.bi[[l]])
          Res.DEv.daily.bi[l]=Model.structure[[i]]$daily.mods.bi[[l]]$deviance
        }
        
        #Daily pos
        n=length(Model.structure[[i]]$daily.mods.pos)
        Aic.daily.pos=rep(NA,n)
        names(Aic.daily.pos)=names(Model.structure[[i]]$daily.mods.pos)
        Res.DEv.daily.pos=Aic.daily.pos
        for( l in 1:n)
        {
          Aic.daily.pos[l]=aicc(Model.structure[[i]]$daily.mods.pos[[l]])
          Res.DEv.daily.pos[l]=Model.structure[[i]]$daily.mods.pos[[l]]$deviance
        }
        PerExp.daily.bi=round(100*(1-Res.DEv.daily.bi[-1]/Res.DEv.daily.bi[i]),1)
        PerExp.daily.pos=round(100*(1-Res.DEv.daily.pos[-1]/Res.DEv.daily.pos[i]),1)
        Tab.daily[[i]]=cbind(data.frame(Species=names(Model.structure)[i],
                                        Model=c('Daily_bi','Daily_pos')),
                             rbind(PerExp.daily.bi,PerExp.daily.pos))
        write.csv(Tab.daily[[i]],paste(hndl.modl.sel,"Dev.expl_",names(SP.list)[Tar.sp[i]],"_daily.csv",sep=""),row.names = F)
      }
    }
    if(select.best.method=="glmmulti")
    {
      fn.show.mod.sel=function(MODS,outs)
      {
        plot(MODS, type="s")
        plot.new()
        tmp <- weightable(MODS)
        mytheme <- gridExtra::ttheme_default(
          core = list(fg_params=list(cex = .65)),
          colhead = list(fg_params=list(cex = .8)),
          rowhead = list(fg_params=list(cex = .8)))
        myt <- gridExtra::tableGrob(tmp[1:outs,], theme = mytheme)
        grid.draw(myt)
      }
      fitFun= function(formula, data,always="", ...) 
      {
        glm(as.formula(paste(deparse(formula), always)), data=data, ...)
      }
      fn.modl.sel=function(RESPNS)   #function for model structure selection
      {
        PREDS[id.cov]=paste("LN",PREDS[id.cov],sep="")
        PREDS=PREDS[-match(always,PREDS)]
        if(RESPNS=="LNcpue")Formula=formula(paste("LNcpue",paste(PREDS,collapse="+"),sep="~"))
        if(RESPNS=="catch")Formula=formula(paste(Response,paste(paste(PREDS,collapse="+"),
                                                                paste("offset(","LN",efrt,")",sep=""),sep="+"),sep="~"))
        
        res <- glmulti(Formula,data=d,level=ifelse(Inter=="MainTerm",1,ifelse(Inter=="2way",2,"3way")),
                       method="h",fitfunction=fitFun,
                       always=paste('+',paste(always,collapse="+"),sep=""),
                       crit="aicc",confsetsize=2^length(PREDS),plotty=F,report=T)
        
        return(list(res=res,BEST=res@formulas[[1]]))
      }
      efrt="km.gillnet.hours.c"
      Inter="MainTerm"
      system.time({
        for(s in Tar.sp)
        {
          #monthly
          PREDS=Predictors_monthly
          d=Store_nom_cpues_monthly[[s]]$QL_dat%>%filter(vessel%in%VES.used[[s]] & blockx%in%BLKS.used[[s]])
          d=makecategorical(PREDS[which(PREDS%in%Categorical)],d)
          d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
          id.cov=which(PREDS%in%Covariates)
          d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
          always=cNSTNT
          Store.Best.Model[[s]]=fn.modl.sel(RESPNS="LNcpue")
          best.fm=Store.Best.Model[[s]]$BEST    
          best.fm=all.vars(best.fm)
          best.fm=as.formula(paste(best.fm[1],'~',paste(c(always,best.fm[-1]),collapse='+')))
          Best.Model[[s]]=best.fm
          rm(d,id.cov,PREDS)
          
          #daily
          d=Store_nom_cpues_daily[[s]]$QL_dat%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%BLKS.used.daily[[s]])
          PREDS=Predictors_daily
          d=d %>% mutate(mean.depth=10*round(mean.depth/10),
                         nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                         mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
          d=makecategorical(PREDS[which(PREDS%in%Categorical)],d)
          d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
          id.cov=which(PREDS%in%Covariates)
          d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
          always=cNSTNT.daily
          Store.Best.Model.daily[[s]]=fn.modl.sel(RESPNS="LNcpue")
          best.fm=Store.Best.Model.daily[[s]]$BEST    
          best.fm=all.vars(best.fm)
          best.fm=as.formula(paste(best.fm[1],'~',paste(c(always,best.fm[-1]),collapse='+')))
          Best.Model.daily[[s]]=best.fm
          rm(d,id.cov,PREDS)
        }
        
        #4.22.3.3. show selection outcomes
        hndl.modl.sel=handl_OneDrive("Analyses/Catch and effort/Outputs/Model Selection/")
        for(s in Tar.sp)
        {
          pdf(paste(hndl.modl.sel,names(SP.list)[s],"_monthly.pdf",sep=""))
          fn.show.mod.sel(MODS=Store.Best.Model[[s]]$res,outs=length(Store.Best.Model[[s]]$res@objects))
          dev.off()
          
          pdf(paste(hndl.modl.sel,names(SP.list)[s],"_daily.pdf",sep=""))
          fn.show.mod.sel(MODS=Store.Best.Model.daily[[s]]$res,outs=length(Store.Best.Model.daily[[s]]$res@objects))
          dev.off()
        }})
      
    }
  }
  
}   
if(Def.mod.Str=="NO")   
{
  if(Use.Tweedie)
  {
    #Monthly        
    for(s in nnn) Best.Model[[s]]=formula(cpue~finyear + blockx + s(vessel,bs='re') + s(month,k=12,bs='cc'))
    Best.Model['Bronze Whaler']=list(NULL)  #no species code 
    Best.Model$`Tiger Shark`=formula(cpue~finyear + blockx + vessel + s(month,k=12,bs='cc'))
    Best.Model$`Dusky Whaler`=Best.Model$`Whiskery Shark`=formula(cpue~finyear + blockx + s(vessel,bs='re'))
    
    #Daily
    for(s in nnn) Best.Model.daily[[s]]=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+
                                                  s(long10.corner,lat10.corner)+s(mean.depth))
    if("Greynurse Shark"%in%names(Best.Model.daily)) Best.Model.daily['Greynurse Shark']=list(NULL)  #protected
    if("Smooth Hammerhead Shark"%in%names(Best.Model.daily))if(add.CITES.term) Best.Model.daily$"Smooth Hammerhead Shark"=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+s(long10.corner,lat10.corner)+s(mean.depth)+cites)
    #add targeting to target species
    Best.Model.daily$`Sandbar Shark`=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+
                                               s(long10.corner,lat10.corner)+s(mean.depth)+cluster_clara)
    Best.Model.daily$`Whiskery Shark`=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+
                                                s(long10.corner,lat10.corner)+s(mean.depth)+cluster_clara)
    Best.Model.daily$`Gummy Shark`=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+
                                             s(long10.corner,lat10.corner)+cluster_clara)
    Best.Model.daily$`Dusky Whaler`=formula(cpue~finyear+s(vessel,bs='re')+s(month,k=12,bs='cc')+
                                              s(long10.corner,lat10.corner)+s(mean.depth)+cluster_clara)
    
    
  }
  
  if(Use.Delta)
  {
    #Other species
    for(s in nnn[-sort(Tar.sp)])
    {
      Best.Model[[s]]=formula(LNcpue ~ finyear + vessel + blockx + month)
      Best.Model.daily[[s]]=Best.Model[[s]]
      Best.Model.daily.gam[[s]]=formula(LNcpue ~ finyear + vessel + s(long10.corner,lat10.corner) + month)
    }
    
    #only 1 vessel level for 7gill sharks after data selection
    Best.Model$`Sevengill Sharks`=Best.Model.daily$`Sevengill Sharks`=
      formula(LNcpue ~ finyear + blockx + month)
    Best.Model.daily.gam$`Sevengill Sharks`=formula(LNcpue ~ finyear + 
                                                      s(long10.corner,lat10.corner) + month)
    
    #Target species
    #Monthly
    Best.Model$`Gummy Shark`=formula(LNcpue ~ finyear + vessel + blockx  + month)
    Best.Model$`Whiskery Shark`=  Best.Model$`Dusky Whaler Bronze Whaler`=
      Best.Model$`Sandbar Shark`=Best.Model$`Gummy Shark`
    
    #Daily
    Best.Model.daily.gam$`Gummy Shark`=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+
                                                 shots.c+lunar+s(mean.depth)+mesh)
    Best.Model.daily.gam$`Whiskery Shark`=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+
                                                    shots.c+lunar+mesh)
    Best.Model.daily.gam$`Dusky Whaler Bronze Whaler`=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+
                                                                shots.c+lunar+s(mean.depth))
    Best.Model.daily.gam$`Sandbar Shark`=formula(LNcpue ~finyear+vessel+s(long10.corner,lat10.corner)+month+
                                                   shots.c+lunar+s(mean.depth))
    
    #Delta monthly
    Best.Model_delta$`Whiskery Shark`=list(bi=formula(catch.target~finyear+vessel+blockx+month+offset(LN.effort)),
                                           pos=formula(LNcpue~finyear+vessel+blockx+month))
    Best.Model_delta$`Gummy Shark`=  Best.Model_delta$`Dusky Whaler Bronze Whaler`=
      Best.Model_delta$`Sandbar Shark`=Best.Model_delta$`Whiskery Shark`
    
    
    #Delta daily
    #GLM
    Best.Model.daily_delta$`Whiskery Shark`=list(
      bi=formula(catch.target~finyear+vessel+block10+month+mean.depth+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+block10+month+mesh+mean.depth+
                    dim.1+dim.2+dim.3))
    Best.Model.daily_delta$`Gummy Shark`=list(
      bi=formula(catch.target~finyear+vessel+block10+month+mean.depth+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+block10+month+mesh+mean.depth+
                    dim.1+dim.2))
    Best.Model.daily_delta$`Dusky Whaler Bronze Whaler`=list(
      bi=formula(catch.target~finyear+vessel+block10+month+mean.depth+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+block10+month+mean.depth+
                    dim.1+dim.2+dim.3))
    
    Best.Model.daily_delta$`Sandbar Shark`=Best.Model.daily_delta$`Dusky Whaler Bronze Whaler`
    
    
    #GAM   #note: targeting predictors excluded from binomial
    Best.Model.daily.gam_delta$`Whiskery Shark`=list(
      bi=formula(catch.target~finyear+vessel+s(long10.corner,lat10.corner)+
                   month+s(mean.depth,k=6)+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+
                    s(dim.1,k=6)+s(dim.2,k=6)+s(dim.3,k=6)))
    Best.Model.daily.gam_delta$`Gummy Shark`=list(
      bi=formula(catch.target~finyear+vessel+s(long10.corner,lat10.corner)+
                   month+s(mean.depth,k=6)+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+mesh+s(dim.1,k=6)))
    Best.Model.daily.gam_delta$`Dusky Whaler Bronze Whaler`=list(
      bi=formula(catch.target~finyear+vessel+s(long10.corner,lat10.corner)+
                   month+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+
                    s(dim.1,k=6)+s(dim.2,k=6)))
    Best.Model.daily.gam_delta$`Sandbar Shark`=list(
      bi=formula(catch.target~finyear+vessel+s(long10.corner,lat10.corner)+
                   month+s(mean.depth,k=6)+offset(LN.effort)),
      pos=formula(LNcpue~finyear+vessel+s(long10.corner,lat10.corner)+month+s(mean.depth,k=6)+
                    s(dim.1,k=6)+s(dim.2,k=6)))
    
  }
}


# Export table of term levels ----------------------------------------------
if(Model.run=="First")
{
  fn.table.terms=function(d,PREDS)
  {
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d)
    
    store=vector('list',length(PREDS))
    names(store)=PREDS
    for(p in 1:length(store))
    {
      xx=d[,match(PREDS[p],names(d))]
      if(is.factor(xx))
      {
        levels=length(levels(xx))
        type="Categorical"
      }
      if(!is.factor(xx))
      {
        Unik=1
        type="Continuous"
        levels=""
      }
      a=data.frame(Term=PREDS[p],Type=type,Levels=levels)
      store[[p]]=a
    }
    return(do.call(rbind,store))
  }
  terms.table=vector('list',length(Tar.sp))
  for(s in Tar.sp)   
  {
    #monthly
    DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used[[s]])
    Mon=fn.table.terms(d=DAT,PREDS=tolower(Predictors_monthly))
    Mon=cbind(Data="Monthly",Mon)
    
    #daily
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    DAT=DAT %>% mutate(mean.depth=10*round(mean.depth/10),
                       nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                       mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    Day=fn.table.terms(d=DAT,PREDS=tolower(Predictors_daily))
    Day=cbind(Data="Daily",Day)
    DAT=rbind(Mon,Day)
    terms.table[[s]]=cbind(Species=names(SP.list)[s],DAT)
  }
  Tab.Terms=do.call(rbind,terms.table)
  row.names(Tab.Terms)=NULL
  fn.word.table(WD=getwd(),TBL=Tab.Terms,Doc.nm="Terms_table",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}

# Show proportion of 0 catch records
fn.sel.yrs.used.glm=function(DD)
{
  a=with(DD,table(finyear,vessel))
  a[a>0]=1
  a=rowSums(a)
  a[a<Threshold.n.vessls.per.yr]=NA
  return(names(a[which(!is.na(a))]))
}
if(Model.run=="First")
{
  fn.fig("Appendix 4.Prop of records with 0 catch",2000, 2400)    
  par(mfrow=c(2,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
  for(s in Tar.sp)
  {
    DAT=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])
    colnames(DAT)=tolower(colnames(DAT))
    DAT=DAT%>%filter(finyear%in%fn.sel.yrs.used.glm(DAT))%>%
      mutate(catch.target=ifelse(catch.target>0,1,0))%>%
      count(catch.target,finyear)%>%
      group_by(finyear)%>%
      mutate(prop = prop.table(n))%>%
      filter(catch.target==1)%>%data.frame%>%
      mutate(year=as.numeric(substr(finyear,1,4)))
    
    DAT1=DATA.list.LIVEWT.c.daily[[s]]%>%filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
    colnames(DAT1)=tolower(colnames(DAT1))
    DAT1=DAT1%>%filter(finyear%in%fn.sel.yrs.used.glm(DAT1))%>%
      mutate(catch.target=ifelse(catch.target>0,1,0))%>%
      count(catch.target,finyear)%>%
      group_by(finyear)%>%
      mutate(prop = prop.table(n))%>%
      filter(catch.target==1)%>%data.frame%>%
      mutate(year=as.numeric(substr(finyear,1,4)))
    
    
    with(DAT,plot(year,1-prop,type='o',pch=19,ylab="",xlab="",cex=1.25,ylim=c(0,1),xlim=c(min(year),max(DAT1$year))))
    with(DAT1,polygon(x=c(year[1],year[length(year)],year[length(year)],year[1]),
                      y=c(0,0,1,1),col='grey90',border="transparent"))
    with(DAT1,points(year,1-prop,type='o',pch=21,bg="white",cex=1.25))
    legend("top",Nms.sp[s],bty='n',cex=1.25)
    
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
  mtext("Proportion of records with no catch",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
  dev.off()
}


#-- Run standardisation on selected model  ----------------------------------------------
if(Use.Qualif.level)
{
  pred.fun=function(MOD,biascor,PRED,Pred.type)             
  {
    lsm=summary(emmeans(MOD, PRED, type=Pred.type))
    id.low=which(names(lsm)%in%c("lower.CL","asymp.LCL"))
    id.up=which(names(lsm)%in%c("upper.CL","asymp.UCL"))
    if(biascor=="YES")
    {
      sigma.glm=sqrt(summary(MOD)$dispersion) # residuals standard error
      lsm$response=exp(lsm$emmean)*exp(sigma.glm^2/2)
      lsm$lower.CL=exp(lsm[,id.low])*exp(sigma.glm^2/2)
      lsm$upper.CL=exp(lsm[,id.up])*exp(sigma.glm^2/2)
    }
    if(biascor=="NO")
    {
      if(is.na(match("response",names(lsm))))lsm$response=lsm$emmean
      lsm$lower.CL=lsm[,id.low]
      lsm$upper.CL=lsm[,id.up]
    }
    return(lsm)
  }
  fn.stand=function(d,Response,RESPNS,PREDS,efrt,Formula,Formula.gam)   #function for standardisation
  {
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    res=NULL
    if(!is.null(Formula)) res <- glm(Formula,data=d)
    res.gam=NULL
    if(!is.null(Formula.gam)) res.gam <-gam(Formula.gam,data=d,method="REML")
    return(list(res=res,res.gam=res.gam,DATA=d))
  }
  
  #monthly
  system.time({Stand.out=foreach(s=Tar.sp,.packages=c('dplyr')) %dopar%
    {
      DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%VES.used[[s]])  #selet blocks and vessels
      DAT=subset(DAT,blockx%in%BLKS.used[[s]])
      DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select years with a min number of vessels
      return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                      efrt="km.gillnet.hours.c",Formula=Best.Model[[s]],Formula.gam=NULL))
      rm(DAT)
    }
  })
  
  #daily
  system.time({Stand.out.daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %dopar%
    {
      if(Model.run=="First") this.form=Best.Model.daily[[s]] else
        this.form=NULL
      DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
      DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
      DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
      DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                        nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                        mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
      return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",
                      PREDS=Predictors_daily,efrt="km.gillnet.hours.c",
                      Formula=this.form,Formula.gam=Best.Model.daily.gam[[s]]))
      rm(DAT)
    }
  }) 
  
  names(Stand.out)=names(Stand.out.daily)=names(SP.list)[Tar.sp]
}
if(Use.Delta)  
{
  fn.stand.delta=function(d,Response,PREDS,efrt,Formula,Formula.gam,Family)   #standardisation using Delta method
  {
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    d$LN.effort=log(d[,match(efrt,names(d))])
    res=NULL
    if(!is.null(Formula)) res <- glm(Formula,data=d,family=Family)
    res.gam=NULL
    if(!is.null(Formula.gam)) res.gam <-gam(Formula.gam,data=d,family=Family,method="REML")
    return(list(res=res,res.gam=res.gam,DATA=d))
  }
  fn.delta=function(d,Response,PREDS,efrt,Formula,Formula.gam)   #used for other species
  {
    ALLvars=all.vars(Formula)[-1]
    Formula.bi=as.formula(paste('catch.pos',"~",paste(paste(ALLvars,collapse="+"),"offset(LNeffort)",sep="+")))
    id.fctr=which(PREDS%in%Categorical)
    if(names(SP.list[s])=="Tiger Shark") Formula=formula(LNcpue ~ finyear + vessel + month)
    Bi <- d %>%mutate(catch.pos=as.numeric(catch.target>0))
    TAB=table(Bi$catch.pos,Bi$finyear)
    TAB[TAB>0]=1
    drop.yrs=names(which(TAB[2,]==0))
    if(length(drop.yrs)>0)
    {
      Bi=Bi%>%filter(!finyear%in%drop.yrs)
      d=d%>%filter(!finyear%in%drop.yrs)
    }
    Bi=makecategorical(PREDS[id.fctr],Bi)
    Bi$LNeffort=log(Bi[,match(efrt,names(Bi))])
    
    d=d%>%filter(catch.target>0)
    d=makecategorical(PREDS[id.fctr],d)
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    
    res=res_bi=NULL
    if(!is.null(Formula))
    {
      res_bi <- glm(Formula.bi, data=Bi, family="binomial", maxit=100)
      res <- glm(Formula,data=d)
      while(any(is.na(coef(res)))|any(is.na(coef(res_bi))))
      {
        na.coef=names(which(is.na(coef(res_bi))))
        if(isTRUE(any(grepl('blockx', na.coef))))
        {
          drp=substr(na.coef[grepl('blockx', na.coef)],nchar('blockx')+1,15)
          Bi=subset(Bi,!blockx%in%drp)%>%
            mutate(blockx=droplevels(blockx))
        }
        if(isTRUE(any(grepl('vessel', na.coef))))
        {
          drp=substr(na.coef[grepl('vessel', na.coef)],nchar('vessel')+1,15)
          Bi=subset(Bi,!vessel%in%drp)%>%
            mutate(vessel=droplevels(vessel))
        }
        if(isTRUE(any(grepl('month', na.coef))))
        {
          drp=substr(na.coef[grepl('vessel', na.coef)],nchar('month')+1,15)
          Bi=subset(Bi,!month%in%drp)%>%
            mutate(month=droplevels(month))
        }
        res_bi <- glm(Formula.bi, data=Bi, family="binomial", maxit=100)
        
        
        na.coef=names(which(is.na(coef(res))))
        if(isTRUE(any(grepl('blockx', na.coef))))
        {
          drp=substr(na.coef[grepl('blockx', na.coef)],nchar('blockx')+1,15)
          d=subset(d,!blockx%in%drp)%>%
            mutate(blockx=droplevels(blockx))
        }
        if(isTRUE(any(grepl('vessel', na.coef))))
        {
          drp=substr(na.coef[grepl('vessel', na.coef)],nchar('vessel')+1,15)
          d=subset(d,!vessel%in%drp)%>%
            mutate(vessel=droplevels(vessel))
        }
        if(isTRUE(any(grepl('month', na.coef))))
        {
          drp=substr(na.coef[grepl('vessel', na.coef)],nchar('month')+1,15)
          d=subset(d,!month%in%drp)%>%
            mutate(month=droplevels(month))
        }
        res <- glm(Formula,data=d)
      }
    }
    
    res.gam=res.gam_bi=NULL
    if(!is.null(Formula.gam))
    {
      ALLvars.gam=all.vars(Formula.gam)[-1]
      ALLvars.gam=ALLvars.gam[-match(c('long10.corner','lat10.corner'),ALLvars.gam)]
      Formula.bi.gam=as.formula(paste('catch.pos',"~",paste(paste(paste(ALLvars.gam,collapse="+"),
                                                                  "s(long10.corner,lat10.corner)",sep='+'),"offset(LNeffort)",sep="+")))
      res.gam_bi <-gam(Formula.bi.gam,data=Bi, family="binomial",method="REML")
      res.gam <-gam(Formula.gam,data=d,method="REML")
    }
    
    return(list(res=res,res_bi=res_bi,res.gam=res.gam,res.gam_bi=res.gam_bi,DATA=d,DATA_bi=Bi))
  }
  fn.MC.delta.cpue=function(BiMOD,MOD,BiData,PosData,niter,pred.term,ALL.terms)
  {
    #don't use binomial part if no contrast (applicable to monthly records where catch is mostly positive)
    if(class(BiMOD)[1]=="gam") Bi.cont=quantile(summary(BiMOD)$se,probs=.6)
    if(class(BiMOD)[1]=="glm") Bi.cont=quantile(summary(BiMOD)$coefficients[,2],probs=.6)  
    
    #get terms
    ALL.terms=ALL.terms[which(ALL.terms%in%colnames(BiData))]
    Covar.bi=as.matrix(vcov(BiMOD))
    Covar.pos=as.matrix(vcov(MOD))
    dummy.BiMOD=BiMOD
    dummy.MOD=MOD
    knstnt.terms=ALL.terms[-match(pred.term,ALL.terms)]
    id.fctr=knstnt.terms[which(knstnt.terms%in%Categorical)]
    id.cont=knstnt.terms[which(!knstnt.terms%in%Categorical)]
    newdata.pos=matrix(nrow=1,ncol=length(knstnt.terms))
    colnames(newdata.pos)=c(id.fctr,id.cont) 
    newdata.pos=as.data.frame(newdata.pos)
    newdata.bi=newdata.pos
    for(ii in 1:ncol(newdata.pos))
    {
      if(colnames(newdata.pos)[ii]%in%id.fctr)
      {
        id=match(colnames(newdata.bi)[ii],names(BiData))
        if(!is.na(id))
        {
          dummy=sort(table(BiData[,id]))
          newdata.bi[,ii]= factor(names(dummy[length(dummy)]),levels(BiData[,id]))
        }
        id=match(colnames(newdata.pos)[ii],names(PosData))
        if(!is.na(id))
        {
          dummy=sort(table(PosData[,id]))
          newdata.pos[,ii]= factor(names(dummy[length(dummy)]),levels(PosData[,id]))
        }
      }
      if(colnames(newdata.pos)[ii]%in%id.cont)
      {
        id=match(colnames(newdata.bi)[ii],names(BiData))
        newdata.bi[,ii]= mean(BiData[,id])
        
        id=match(colnames(newdata.pos)[ii],names(PosData))
        newdata.pos[,ii]= mean(PosData[,id])
      }
    }
    nms.coef=names(unlist(dummy.coef(MOD)))
    pred.var=sapply(strsplit(nms.coef[grepl(pred.term, nms.coef)], paste(pred.term,".",sep="")), "[", 2)
    
    pred.dat.pos=data.frame(factor(pred.var,levels=pred.var))
    colnames(pred.dat.pos)=pred.term
    newdata.pos=cbind(pred.dat.pos,newdata.pos)
    
    newdata.bi=cbind(pred.dat.pos,newdata.bi)
    newdata.bi=cbind(newdata.bi,LNeffort=mean(BiData$LNeffort),LN.effort=mean(BiData$LNeffort))
    
    set.seed(999)
    
    Bi.pars.rand=rmvnorm(niter,mean=coef(BiMOD),sigma=Covar.bi)
    Pos.pars.rand=rmvnorm(niter,mean=coef(MOD),sigma=Covar.pos)
    
    MC.preds=matrix(nrow=niter,ncol=nrow(newdata.bi))
    
    
    for(n in 1:niter)
    {
      #Binomial part
      if(Bi.cont<5) dummy.BiMOD$coefficients=Bi.pars.rand[n,] else
        dummy.BiMOD=BiMOD
      newdata.bi$Pred.bi=predict(dummy.BiMOD,newdata=newdata.bi, type="response")
      
      #Positive part
      dummy.MOD$coefficients=Pos.pars.rand[n,]
      a=predict(dummy.MOD,newdata=newdata.pos, type="response",se.fit=T)
      newdata.pos$Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
      
      dummy=left_join(newdata.bi,newdata.pos,by=pred.term)%>%
        #mutate(Index=Pred)%>%
        mutate(Index=Pred.bi*Pred)%>%
        select(Index)%>%
        unlist
      
      MC.preds[n,]=dummy
    }
    
    #Get summary stats
    MEAN=colMeans(MC.preds,na.rm=T)
    SD=apply(MC.preds,2,sd,na.rm=T)
    LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T))
    UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T))
    
    Stats=cbind(pred.dat.pos,data.frame(MEAN=MEAN,SD=SD,LOW=LOW,UP=UP))
    Stats=Stats[order(Stats[,1]),]
    Stats=Stats%>%rename(response=MEAN, lower.CL=LOW, upper.CL=UP)
    rownames(Stats)=NULL
    return(Stats)
  }
  fn.MC.delta.cpue_spatial=function(BiMOD,MOD,BiData,PosData,niter,pred.term,ALL.terms,Spatial.grid)
  {
    #don't use binomial part if no contrast (applicable to monthly records where catch is mostly positive)
    if(class(BiMOD)[1]=="gam") Bi.cont=quantile(summary(BiMOD)$se,probs=.6)
    if(class(BiMOD)[1]=="glm") Bi.cont=quantile(summary(BiMOD)$coefficients[,2],probs=.6)  
    
    #get terms
    ALL.terms=ALL.terms[which(ALL.terms%in%colnames(BiData))]
    Covar.bi=as.matrix(vcov(BiMOD))
    Covar.pos=as.matrix(vcov(MOD))
    dummy.BiMOD=BiMOD
    dummy.MOD=MOD
    knstnt.terms=ALL.terms[-match(pred.term,ALL.terms)]
    id.fctr=knstnt.terms[which(knstnt.terms%in%Categorical)]
    id.cont=knstnt.terms[which(!knstnt.terms%in%Categorical)]
    newdata.pos=matrix(nrow=1,ncol=length(knstnt.terms))
    colnames(newdata.pos)=c(id.fctr,id.cont) 
    newdata.pos=as.data.frame(newdata.pos)
    newdata.bi=newdata.pos
    for(ii in 1:ncol(newdata.pos))
    {
      if(colnames(newdata.pos)[ii]%in%id.fctr)
      {
        id=match(colnames(newdata.bi)[ii],names(BiData))
        if(!is.na(id))
        {
          dummy=sort(table(BiData[,id]))
          newdata.bi[,ii]= factor(names(dummy[length(dummy)]),levels(BiData[,id]))
        }
        id=match(colnames(newdata.pos)[ii],names(PosData))
        if(!is.na(id))
        {
          dummy=sort(table(PosData[,id]))
          newdata.pos[,ii]= factor(names(dummy[length(dummy)]),levels(PosData[,id]))
        }
      }
      if(colnames(newdata.pos)[ii]%in%id.cont)
      {
        id=match(colnames(newdata.bi)[ii],names(BiData))
        newdata.bi[,ii]= mean(BiData[,id])
        
        id=match(colnames(newdata.pos)[ii],names(PosData))
        newdata.pos[,ii]= mean(PosData[,id])
      }
    }
    nms.coef=names(unlist(dummy.coef(MOD)))
    newdata.pos=cbind(Spatial.grid,newdata.pos)
    newdata.bi=cbind(Spatial.grid,newdata.bi)
    newdata.bi=cbind(newdata.bi,LNeffort=mean(BiData$LNeffort),LN.effort=mean(BiData$LNeffort))
    
    set.seed(999)
    
    Bi.pars.rand=rmvnorm(niter,mean=coef(BiMOD),sigma=Covar.bi)
    Pos.pars.rand=rmvnorm(niter,mean=coef(MOD),sigma=Covar.pos)
    MC.preds=matrix(nrow=niter,ncol=nrow(newdata.bi))
    for(n in 1:niter)
    {
      #Binomial part
      if(Bi.cont<5) dummy.BiMOD$coefficients=Bi.pars.rand[n,] else
        dummy.BiMOD=BiMOD
      newdata.bi$Pred.bi=predict(dummy.BiMOD,newdata=newdata.bi, type="response")
      
      #Positive part
      dummy.MOD$coefficients=Pos.pars.rand[n,]
      a=predict(dummy.MOD,newdata=newdata.pos, type="response",se.fit=T)
      newdata.pos$Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
      
      dummy=left_join(newdata.bi,newdata.pos,by=pred.term)%>%
        #mutate(Index=Pred)%>%
        mutate(Index=Pred.bi*Pred)%>%
        select(Index)%>%
        unlist
      
      MC.preds[n,]=dummy
    }
    
    #Get summary stats
    MEAN=colMeans(MC.preds,na.rm=T)
    SD=apply(MC.preds,2,sd,na.rm=T)
    LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T))
    UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T))
    
    Stats=cbind(Spatial.grid,data.frame(MEAN=MEAN,SD=SD,LOW=LOW,UP=UP))
    Stats=Stats[order(Stats[,1]),]
    Stats=Stats%>%rename(response=MEAN, lower.CL=LOW, upper.CL=UP)
    rownames(Stats)=NULL
    return(Stats)
  }
  Plot.cpue.delta=function(cpuedata,cpuedata.daily,CL,CxS,
                           Yvar,add.lines,firstyear,ADD.nomnl)    
  {
    if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
    if(length(cpuedata)<=3)tc=seq(-.5*0.25,.5*0.25,length.out=length(cpuedata))
    
    Yrs=c(as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4)),
          as.numeric(substr(cpuedata.daily[[1]][,match(Yvar,names(cpuedata.daily[[1]]))],1,4)))
    Tops=c(unlist(lapply(cpuedata, `[`, "upper.CL")),
           unlist(lapply(cpuedata.daily, `[`, "upper.CL")),
           ADD.nomnl$response)
    Tops=Tops[!is.infinite(Tops)]
    ymax=max(Tops)
    Quant=quantile(Tops,probs=c(.9,1))
    if(diff(Quant)>3) ymax=quantile(Tops,probs=.99)
    
    
    plot(Yrs,Yrs,ylim=c(0,ymax),xlim=c(firstyear,max(Yrs)),ylab="",xlab="",
         col="transparent",cex.axis=1.25)
    for(l in 1:length(cpuedata))
    {
      aaa=cpuedata[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
      aaa.daily=cpuedata.daily[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
      
      msn=Yrs[which(!Yrs%in%c(aaa$finyear,aaa.daily$finyear))]
      if(length(msn)>0)
      {
        ad=aaa[length(msn),]
        ad[,]=NA
        ad$finyear=msn
        aaa=rbind(aaa,ad)
        aaa=aaa[order(aaa$finyear),]
      }
      with(aaa,
           {
             if(add.lines=="NO") points(finyear+tc[l], response, pch=19, lty=2, col=CL[l],cex=CxS)
             if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=19, lty=2, col=CL[l],cex=CxS)
             arrows(x0=finyear+tc[l], y0=lower.CL, 
                    x1=finyear+tc[l], y1=upper.CL, 
                    code=3, angle=90, length=0.05, col=CL[l])
           })
      with(aaa.daily,
           {
             if(l==1)polygon(x=c(finyear[1]-.5,finyear[length(finyear)]+.5,finyear[length(finyear)]+.5,finyear[1]-.5),
                             y=c(0,0,ymax*.99,ymax*.99),col='grey92',border="transparent")
             arrows(x0=finyear+tc[l], y0=lower.CL, 
                    x1=finyear+tc[l], y1=upper.CL, 
                    code=3, angle=90, length=0.05, col=CL[l])
             if(add.lines=="NO") points(finyear+tc[l], response, pch=21, lty=2, col=CL[l],cex=CxS)
             if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=21,bg="white", lty=2, col=CL[l],cex=CxS)
           })
    }
    if(!is.null(ADD.nomnl))
    {
      with(ADD.nomnl,points(as.numeric(substr(finyear,1,4))+.1,response,pch=19,col=rgb(.1,.1,.1,alpha=.2),cex=CxS)) 
    }
    
  }
  
  #Target species
  #monthly
  system.time({Stand.out=foreach(s=Tar.sp,.packages=c('dplyr')) %dopar%
    {
      DAT=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])
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
      return(list(Bi=Bi,Pos=Pos))
      rm(DAT)
    }
  })   
  
  #daily
  system.time({Stand.out.daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %dopar%
    {
      DAT=DATA.list.LIVEWT.c.daily[[s]]%>%filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
      colnames(DAT)=tolower(colnames(DAT))
      #select years with a min number of vessels
      DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT))
      DAT=DAT%>% mutate(nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                        mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
      
      #Binomial
      DAT.bi=DAT%>%mutate(catch.target=ifelse(catch.target>0,1,0))
      if(Model.run=="First") this.form=Best.Model.daily_delta[[s]]$bi else
        this.form=NULL
      Bi=fn.stand.delta(d=DAT.bi,Response="catch.target",PREDS=Predictors_daily,
                        efrt="km.gillnet.hours.c",Formula=this.form,
                        Formula.gam=Best.Model.daily.gam_delta[[s]]$bi,Family="binomial")
      #Positive catch
      DAT=DAT%>%filter(catch.target>0)
      if(Model.run=="First") this.form=Best.Model.daily_delta[[s]]$pos else
        this.form=NULL
      Pos=fn.stand.delta(d=DAT,Response="catch.target",PREDS=Predictors_daily,
                         efrt="km.gillnet.hours.c",Formula=this.form,
                         Formula.gam=Best.Model.daily.gam_delta[[s]]$pos,Family="gaussian")
      return(list(Bi=Bi,Pos=Pos))
      rm(DAT)
    }})
  
  names(Stand.out)=names(Stand.out.daily)=names(SP.list)[Tar.sp]
  
  
  # Other species (based on delta method due to excess zeros)
  #monthly          takes 2 sec
  system.time({Stand.out.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr')) %dopar%
    {
      if(!is.null(DATA.list.LIVEWT.c[[s]]) & !is.null(BLKS.used[[s]]))
      {
        DAT=DATA.list.LIVEWT.c[[s]]
        colnames(DAT)=tolower(colnames(DAT)) 
        DAT=DAT%>%filter(vessel%in%VES.used[[s]] & blockx%in%as.numeric(BLKS.used[[s]]))  #selet blocks and vessels
        return(fn.delta(d=DAT,Response="catch.target",PREDS=Predictors_monthly,
                        efrt="km.gillnet.hours.c",Formula=Best.Model[[s]],Formula.gam=NULL))
        rm(DAT)
      }
    }
  })
  #daily            takes 3.5 sec
  system.time({Stand.out.daily.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mgcv')) %dopar%
    {
      if(!is.null(DATA.list.LIVEWT.c.daily[[s]])& !is.null(BLKS.used.daily[[s]]))
      {
        DAT=DATA.list.LIVEWT.c.daily[[s]]
        colnames(DAT)=tolower(colnames(DAT)) 
        DAT=DAT%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%as.numeric(BLKS.used.daily[[s]]))  #selet blocks and vessels
        return(fn.delta(d=DAT,Response="catch.target",PREDS=Predictors_monthly,
                        efrt="km.gillnet.hours.c",
                        Formula=NULL,Formula.gam=Best.Model.daily.gam[[s]]))
        rm(DAT)
      }
    }
  })
  names(Stand.out.other)=names(Stand.out.daily.other)=names(SP.list)[nnn[-sort(Tar.sp)]]
  
  
  #Combine target and other species lists
  Stand.out=c(Stand.out,Stand.out.other)
  Stand.out=Stand.out[names(SP.list)]
  Stand.out.daily=c(Stand.out.daily,Stand.out.daily.other)
  Stand.out.daily=Stand.out.daily[names(SP.list)]
  
  rm(Stand.out.other,Stand.out.daily.other)
  
}
if(Use.Tweedie)    #takes 5 mins with gam(),  0.63 mins with bam()
{
  #predictions
  #   use marginal means which accounts for unbalanced data
  #   see https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
  pred.fun=function(mod,biascor,PRED)             
  {
    if(biascor=="YES")  #apply bias correction for log transf
    {
      lsm=summary(emmeans(mod, PRED, type="link", rg.limit = 2e5))%>%
        mutate(response=exp(emmean)*exp(SE^2/2),
               lower.CL=exp(lower.CL)*exp(SE^2/2),
               upper.CL=exp(upper.CL)*exp(SE^2/2))
    }
    if(biascor=="NO") lsm=summary(emmeans(mod, PRED, type="response", rg.limit = 2e5))
    
    lsm$SD=lsm$SE
    
    return(lsm)
  }
  
  #continuous preds
  pred.fun.continuous=function(d,mod,PRED,Formula)
  {
    NewD=NULL
    VARS=all.vars(Formula)[-1]
    if(PRED[1]%in%VARS)
    {
      #create new data
      VARS=VARS[!VARS%in%PRED]
      NewD=d[,PRED]
      if(length(PRED)==1)
      {
        NewD=seq(min(NewD),max(NewD))
        if(length(NewD)<50) NewD=seq(min(NewD),max(NewD),length.out = 100)
        NewD=as.data.frame(NewD)
        names(NewD)=PRED
      }else
      {
        NewD$dummy=paste(NewD[,1],NewD[,2])
        NewD=NewD%>%
          distinct(dummy,.keep_all = TRUE)%>%
          dplyr::select(-dummy)
      }
      
      FixedVar=as.data.frame(matrix(NA,nrow=1,ncol=length(VARS)))
      names(FixedVar)=VARS
      for(l in 1:length(VARS))
      {
        vv=d[,VARS[l]]
        
        if(is.factor(vv))
        {
          FF=names(sort(-table(vv)))[1]
          FixedVar[,l]=factor(FF,levels(vv))
        }else
          FixedVar[,l]=mean(vv,na.rm=T)
        
        rm(vv)
      }
      NewD=cbind(NewD,FixedVar)
      
      #predict new data
      a=predict(mod,newdata=NewD,type="response",se.fit=T)
      NewD$Pred=a$fit
      NewD$Pred.SE=a$se.fit
    }
    return(NewD)  
  }
  
  
  #monthly  
  system.time({Stand.out=foreach(s=nnn,.packages=c('dplyr','mgcv')) %dopar%
    {
      if(!is.null(BLKS.used[[s]]))
      {
        iid=match(SpiSis[match(names(DATA.list.LIVEWT.c)[s],names(SpiSis))],names(First.year.catch))
        Firs.yr.ktch=names(First.year.catch[[iid]])
        theseyears=sort(unique(DATA.list.LIVEWT.c[[s]]$FINYEAR))
        theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
        
        d=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]] & FINYEAR%in%theseyears )
        
        #remove first years with very few positive catch observation
        if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
        if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
        
        Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
        Continuous=Covariates.monthly
        colnames(d)=tolower(colnames(d))
        Terms=tolower(Terms)
        Continuous=tolower(Continuous)
        Factors=Terms[!Terms%in%Continuous]
        Terms=all.vars(Best.Model[[s]])[-1]
        d <- d[,c('catch.target','km.gillnet.hours.c',Terms)]%>%
          mutate(cpue=catch.target/km.gillnet.hours.c)
        d <- makecategorical(Factors[Factors%in%Terms],d)
        mod<-bam(Best.Model[[s]],data=d,family='tw',method="fREML",discrete=TRUE)
        
        return(list(DATA=d,res.gam=mod))
        
        rm(d,mod)
      }
    }
  })   
  
  #daily 
  if("Smooth Hammerhead Shark"%in%names(DATA.list.LIVEWT.c.daily) & add.CITES.term)
  {
    DATA.list.LIVEWT.c.daily$"Smooth Hammerhead Shark"=DATA.list.LIVEWT.c.daily$"Smooth Hammerhead Shark"%>%
      mutate(cites=ifelse(date<as.Date("2015-01-01"),"No","Yes"))
  }
  system.time({Stand.out.daily=foreach(s=nnn,.packages=c('dplyr','mgcv')) %dopar%
    {
      if(!is.null(BLKS.used.daily[[s]]))
      {
        iid=match(SpiSis[match(names(DATA.list.LIVEWT.c.daily)[s],names(SpiSis))],names(First.year.catch.daily))
        Firs.yr.ktch=names(First.year.catch.daily[[iid]])
        theseyears=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR))
        theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
        
        #remove first year of transition from Monthly to Daily returns due to effort misreporting
        drop.this=match("2006-07",theseyears)
        if(!is.na(drop.this))theseyears=theseyears[-drop.this]
        
        d=DATA.list.LIVEWT.c.daily[[s]]%>%
          filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]] & FINYEAR%in%theseyears)
        Terms=Predictors_daily
        if(add.CITES.term & names(DATA.list.LIVEWT.c.daily)[s]=="Smooth Hammerhead Shark") Terms=c(Terms,'cites')
        Continuous=Covariates.daily   
        colnames(d)=tolower(colnames(d))
        Terms=tolower(Terms)
        Continuous=tolower(Continuous)
        Factors=Terms[!Terms%in%Continuous]
        Terms=all.vars(Best.Model.daily[[s]])[-1]
        d <- d[,c('catch.target','km.gillnet.hours.c',Terms)]%>%
          mutate(cpue=catch.target/km.gillnet.hours.c)
        d <- makecategorical(Factors[Factors%in%Terms],d)
        mod<-bam(Best.Model.daily[[s]],data=d,family='tw',method="fREML",discrete=TRUE)
        
        return(list(DATA=d,res.gam=mod))
        rm(d,mod)
      }
    }
  })
  
  names(Stand.out)=names(Stand.out.daily)=names(SP.list)
}

# rm(Species.list.daily,Species.list,Data.daily.GN,
#    Data.monthly.GN,Effort.monthly,Effort.daily,DATA.list.LIVEWT.c.daily_all_reporters,
#    DATA.list.LIVEWT.c_all_reporters)


#-- Run sensitivity tests  ----------------------------------------------    
if(Model.run=="First")      
{
  sens=Tab.Sensi%>%filter(!Scenario=='Nominal')
  sens=sens%>%mutate(Efrt.used="km.gillnet.hours.c",
                     run.model=ifelse(!Scenario%in%c('Base case','No efficiency'),"YES","NO"))
  Sens.pred=vector('list',length(SP.list))
  names(Sens.pred)=names(SP.list)
  Sens.pred.daily=Sens.pred.normlzd=Sens.pred.daily.normlzd=Sens.pred
  
  
  if(Use.Qualif.level)
  {
    system.time({for(s in Tar.sp)   #takes 55 secs
    {
      #1. Fit models
      sens_monthly=foreach(o=1:nrow(sens),.packages=c('dplyr')) %dopar%   
        {
          MiN.YR=as.numeric(substr(sens$Vessels_used[o],1,1))
          EFrT=sens$Efrt.used[o]
          a=fn.check.balanced(d=Store_nom_cpues_monthly[[s]]$QL_dat,SP=names(SP.list)[s],
                              what="monthly",MN.YR=MiN.YR,pLot=F)
          if(names(SP.list)[s]=="Sandbar shark") a$this.blks=subset(a$this.blks,!a$this.blks=="3517") #cannot estimate this coef for 0==1
          DAT=subset(Store_nom_cpues_monthly[[s]]$QL_dat,vessel%in%a$this.ves)
          DAT=subset(DAT,blockx%in%a$this.blks)
          DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
          return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,efrt=EFrT,
                          Formula=Best.Model[[s]],Formula.gam=NULL))
          rm(DAT)
        }
      sens_daily=foreach(o=1:nrow(sens),.packages=c('dplyr')) %dopar%
        {
          MiN.YR=as.numeric(substr(sens$Vessels_used[o],1,1))
          EFrT=sens$Efrt.used[o]
          a=fn.check.balanced(d=Store_nom_cpues_daily[[s]]$QL_dat,SP=names(SP.list)[s],
                              what="daily",MN.YR=MiN.YR,pLot=F)
          
          DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%a$this.ves)
          DAT=subset(DAT,blockx%in%a$this.blks)
          DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
          DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                            nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                            mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
          return(fn.stand(d=DAT,Response="catch.target",RESPNS="LNcpue",PREDS=Predictors_monthly,
                          efrt=EFrT,Formula=Best.Model[[s]],Formula.gam=NULL))
          rm(DAT)
        }
      names(sens_monthly)=names(sens_daily)=sens$Scenario
      
      #2. Predict years based on emmeans (formerly lsmeans) considering log bias corr if required
      dummy=vector('list',length=nrow(sens))
      names(dummy)=sens$Scenario
      dummy.daily=dummy
      for(o in 1:nrow(sens))
      {
        d=sens_monthly[[o]]$DATA   #note: need data as global for ref_grid
        dummy[[o]]=pred.fun(MOD=sens_monthly[[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
        
        d=sens_daily[[o]]$DATA
        dummy.daily[[o]]=pred.fun(MOD=sens_daily[[o]]$res,biascor="YES",PRED="finyear",Pred.type="link")
        rm(d)
      }
      
      #3. Apply efficiency creep where required     
      for(o in 1:nrow(sens))
      {
        if(sens$Efficiency_increase[o]=="Yes")
        {
          #monthly
          add.crp=Eff.creep$effort.creep[match(dummy[[o]]$finyear,Eff.creep$finyear)]
          dummy[[o]]$response=dummy[[o]]$response*(1-add.crp)
          dummy[[o]]$lower.CL=dummy[[o]]$lower.CL*(1-add.crp)
          dummy[[o]]$upper.CL=dummy[[o]]$upper.CL*(1-add.crp)
          
          #daily
          add.crp=Eff.creep$effort.creep[match(dummy.daily[[o]]$finyear,Eff.creep$finyear)]
          dummy.daily[[o]]$response=dummy.daily[[o]]$response*(1-add.crp)
          dummy.daily[[o]]$lower.CL=dummy.daily[[o]]$lower.CL*(1-add.crp)
          dummy.daily[[o]]$upper.CL=dummy.daily[[o]]$upper.CL*(1-add.crp)
        }
      }
      
      #4. Normalise
      dummy.normlzd=dummy
      dummy.daily.normlzd=dummy.daily
      for(o in 1:nrow(sens))
      {
        #monthly
        Mn=mean(dummy[[o]]$response)
        dummy.normlzd[[o]]$response=dummy[[o]]$response/Mn
        dummy.normlzd[[o]]$lower.CL=dummy[[o]]$lower.CL/Mn
        dummy.normlzd[[o]]$upper.CL=dummy[[o]]$upper.CL/Mn
        
        #daily
        Mn=mean(dummy.daily[[o]]$response)
        dummy.daily.normlzd[[o]]$response=dummy.daily[[o]]$response/Mn
        dummy.daily.normlzd[[o]]$lower.CL=dummy.daily[[o]]$lower.CL/Mn
        dummy.daily.normlzd[[o]]$upper.CL=dummy.daily[[o]]$upper.CL/Mn
      }
      
      Sens.pred[[s]]=dummy
      Sens.pred.daily[[s]]=dummy.daily
      Sens.pred.normlzd[[s]]=dummy.normlzd
      Sens.pred.daily.normlzd[[s]]=dummy.daily.normlzd
    }})
    #Plot   
    fn.fig("Appendix_Sensitivity",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      LgND="NO"
      if(s==5)LgND="YES"
      Plot.cpue(cpuedata=Sens.pred[[s]],ADD.LGND=LgND,whereLGND='topright',
                COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
      if(s==5) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      LgND="NO"
      Plot.cpue(cpuedata=Sens.pred.daily[[s]],ADD.LGND=LgND,whereLGND='topright',
                COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
      if(s==5) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
    fn.fig("Appendix_Sensitivity_nomalised",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      LgND="NO"
      if(s==5)LgND="YES"
      Plot.cpue(cpuedata=Sens.pred.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
                COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
      if(s==5) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      LgND="NO"
      Plot.cpue(cpuedata=Sens.pred.daily.normlzd[[s]],ADD.LGND=LgND,whereLGND='topright',
                COL="color",CxS=1.15,Yvar="finyear",add.lines="YES")
      if(s==5) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    dev.off()
  }
  if(Use.Delta)  #takes 7 mins
  {
    #1. Fit models
    #monthly
    system.time({sens_monthly=foreach(s=Tar.sp,.packages=c('dplyr')) %dopar%
      {
        DAT=DATA.list.LIVEWT.c[[s]]
        if(Nms.sp[s]=="Sandbar Shark") DAT=DAT%>%filter(!BLOCKX==9600 & 
                                                          !FINYEAR%in%c('1986-87')) #cannot estimate this coefficient
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
        return(list(Bi=Bi,Pos=Pos))
        rm(DAT)
      }
    })
    #daily
    system.time({sens_daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %dopar%
      {
        DAT=DATA.list.LIVEWT.c.daily[[s]]
        colnames(DAT)=tolower(colnames(DAT))
        #select years with a min number of vessels
        DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT))
        DAT=DAT%>% mutate(nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                          mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
        
        #Binomial
        #note: had to remove spline bit from bin form because with all records, there's not enough info for gam
        Bi.form=formula(catch.target ~ finyear + vessel + month + offset(LN.effort))
        DAT.bi=DAT%>%mutate(catch.target=ifelse(catch.target>0,1,0))
        Bi=fn.stand.delta(d=DAT.bi,Response="catch.target",PREDS=Predictors_daily,
                          efrt="km.gillnet.hours.c",Formula=NULL,
                          Formula.gam=Bi.form,Family="binomial")
        #Positive catch
        DAT=DAT%>%filter(catch.target>0)
        Pos=fn.stand.delta(d=DAT,Response="catch.target",PREDS=Predictors_daily,
                           efrt="km.gillnet.hours.c",Formula=NULL,
                           Formula.gam=Best.Model.daily.gam_delta[[s]]$pos,Family="gaussian")
        return(list(Bi=Bi,Pos=Pos))
        rm(DAT)
      }})
    
    names(sens_monthly)=names(sens_daily)=names(SP.list)[Tar.sp]
    Store.sen=list(Stand.out,sens_monthly)
    Store.sen.daily=list(Stand.out.daily,sens_daily)
    names(Store.sen)=names(Store.sen.daily)=sens$Scenario[-match("No efficiency",sens$Scenario)]
    
    #2. Predict years considering log bias corr if required  
    All.preds=vector('list',length(Store.sen))
    names(All.preds)=names(Store.sen)
    system.time({for(i in 1:length(Store.sen))
    {
      #Monthly
      Pred.sens=foreach(s=1:length(Tar.sp),.packages=c('dplyr','mvtnorm')) %dopar%
        {
          return(fn.MC.delta.cpue(BiMOD=Store.sen[[i]][[s]]$Bi$res,
                                  MOD=Store.sen[[i]][[s]]$Pos$res,
                                  BiData=Store.sen[[i]][[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                  PosData=Store.sen[[i]][[s]]$Pos$DATA,
                                  niter=100,
                                  pred.term='finyear',
                                  ALL.terms=Predictors_monthly))
        }
      names(Pred.sens)=Nms.sp[Tar.sp]
      
      #Daily
      Pred.sens.daily=foreach(s=1:length(Tar.sp),.packages=c('dplyr','mvtnorm')) %dopar%
        {
          return(fn.MC.delta.cpue(BiMOD=Store.sen.daily[[i]][[s]]$Bi$res.gam,
                                  MOD=Store.sen.daily[[i]][[s]]$Pos$res.gam,
                                  BiData=Store.sen.daily[[i]][[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                  PosData=Store.sen.daily[[i]][[s]]$Pos$DATA,
                                  niter=100,
                                  pred.term='finyear',
                                  ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
          
        }
      names(Pred.sens)=names(Pred.sens.daily)=Nms.sp[Tar.sp]
      All.preds[[i]]=list(monthly=Pred.sens,daily=Pred.sens.daily)
    }
    })
    All.preds$`No efficiency`=All.preds$`Base case`
    
    
    #3. Apply efficiency creep where required
    All.preds.creep=All.preds
    for(o in 1:nrow(sens))
    {
      if(sens$Efficiency_increase[o]=="Yes")
      {
        for(i in 1:length(All.preds.creep[[o]]$monthly))
        {
          #monthly
          All.preds.creep[[o]]$monthly[[i]]=All.preds.creep[[o]]$monthly[[i]]%>%
            mutate(add.crp=Eff.creep$effort.creep[match(finyear,Eff.creep$finyear)],
                   lower.CL=lower.CL-(response-response*(1-add.crp)),
                   upper.CL=upper.CL-(response-response*(1-add.crp)),
                   response=response*(1-add.crp))
          #daily
          All.preds.creep[[o]]$daily[[i]]=All.preds.creep[[o]]$daily[[i]]%>%
            mutate(add.crp=Eff.creep$effort.creep[match(finyear,Eff.creep$finyear)],
                   lower.CL=lower.CL-(response-response*(1-add.crp)),
                   upper.CL=upper.CL-(response-response*(1-add.crp)),
                   response=response*(1-add.crp))
        }
      }
    }
    
    #Plot
    fn.fig("Appendix_Sensitivity",2000, 2400)    
    par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in 1:length(Tar.sp))
    {
      Da=list(All.preds.creep$`Base case`$monthly[[s]],
              All.preds.creep$`All vessels & blocks`$monthly[[s]],
              All.preds.creep$`No efficiency`$monthly[[s]])
      Da.daily=list(All.preds.creep$`Base case`$daily[[s]],
                    All.preds.creep$`All vessels & blocks`$daily[[s]],
                    All.preds.creep$`No efficiency`$daily[[s]])
      names(Da)=names(Da.daily)=sens$Scenario
      
      Plot.cpue.delta(cpuedata=Da,cpuedata.daily=Da.daily,
                      CL=c("black","forestgreen","red"),CxS=1.15,Yvar="finyear",
                      add.lines="YES",firstyear=1975,
                      ADD.nomnl=NULL)
      
      legend("top",Nms.sp[Tar.sp[s]],bty='n',cex=1.75)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("CPUE (kg/km gillnet hour)",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    legend('topright',sens$Scenario,bty='n',col=c("black","forestgreen","red"),
           pch=19,cex=1.25)
    dev.off()
    
  }
  
  if(Use.Tweedie)  #takes 10 mins
  {
    #1. Fit models
    #monthly
    system.time({sens_monthly=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %dopar%
      {
        d=DATA.list.LIVEWT.c[[s]]
        if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
        if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
        Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
        Continuous=Covariates.monthly
        colnames(d)=tolower(colnames(d))
        Terms=tolower(Terms)
        Continuous=tolower(Continuous)
        Factors=Terms[!Terms%in%Continuous]
        Terms=all.vars(Best.Model[[s]])[-1]
        d <- d%>%
          dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
          mutate(cpue=catch.target/km.gillnet.hours.c)
        d <- makecategorical(Factors[Factors%in%Terms],d)
        mod<-bam(Best.Model[[s]],data=d,family='tw',method="fREML",discrete=TRUE)
        
        return(list(DATA=d,res.gam=mod))
        rm(d,mod)
      }
    })
    
    #daily
    system.time({sens_daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv')) %do%
      {
        theseyears=as.character(unique(Stand.out.daily[[s]]$DATA$finyear))
        d=DATA.list.LIVEWT.c.daily[[s]]%>%
          filter(FINYEAR%in%theseyears)
        Terms=Predictors_daily
        Continuous=Covariates.daily
        colnames(d)=tolower(colnames(d))
        Terms=tolower(Terms)
        Continuous=tolower(Continuous)
        Factors=Terms[!Terms%in%Continuous]
        Terms=all.vars(Best.Model.daily[[s]])[-1]
        d <- d%>%
          dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
          mutate(cpue=catch.target/km.gillnet.hours.c)
        d <- makecategorical(Factors[Factors%in%Terms],d)
        mod<-bam(Best.Model.daily[[s]],data=d,family='tw',method="fREML",discrete=TRUE)
        
        return(list(DATA=d,res.gam=mod))
        rm(d)
      }})
    
    names(sens_monthly)=names(sens_daily)=names(SP.list)[Tar.sp]
    Store.sen=list(Stand.out[Tar.sp],sens_monthly)
    Store.sen.daily=list(Stand.out.daily[Tar.sp],sens_daily)
    names(Store.sen)=names(Store.sen.daily)=sens$Scenario[-match("No efficiency",sens$Scenario)]
    
    #2. Predict years considering log bias corr if required  
    All.preds=vector('list',length(Store.sen))
    names(All.preds)=names(Store.sen)
    for(i in 1:length(Store.sen))
    {
      #Monthly
      Pred.sens=foreach(s=1:length(Tar.sp),.packages=c('dplyr','emmeans')) %do%
        {
          d=Store.sen[[i]][[s]]$DATA   #note: need data as global for ref_grid
          return(pred.fun(mod=Store.sen[[i]][[s]]$res.gam,
                          biascor="NO",
                          PRED="finyear"))
          rm(d)
        }
      names(Pred.sens)=Nms.sp[Tar.sp]
      
      #Daily
      Pred.sens.daily=foreach(s=1:length(Tar.sp),.packages=c('dplyr','emmeans')) %do%
        {
          d=Store.sen.daily[[i]][[s]]$DATA
          return(pred.fun(mod=Store.sen.daily[[i]][[s]]$res.gam,
                          biascor="NO",
                          PRED="finyear"))
          rm(d)
        }
      names(Pred.sens)=names(Pred.sens.daily)=Nms.sp[Tar.sp]
      All.preds[[i]]=list(monthly=Pred.sens,daily=Pred.sens.daily)
    }
    
    All.preds$`No efficiency`=All.preds$`Base case`
    
    
    #3. Apply efficiency creep where required
    All.preds.creep=All.preds
    for(o in 1:nrow(sens))
    {
      if(sens$Efficiency_increase[o]=="Yes")
      {
        for(i in 1:length(All.preds.creep[[o]]$monthly))
        {
          #monthly
          All.preds.creep[[o]]$monthly[[i]]=All.preds.creep[[o]]$monthly[[i]]%>%
            mutate(add.crp=Eff.creep$effort.creep[match(finyear,Eff.creep$finyear)],
                   lower.CL=lower.CL-(response-response*(1-add.crp)),
                   upper.CL=upper.CL-(response-response*(1-add.crp)),
                   response=response*(1-add.crp))
          #daily
          All.preds.creep[[o]]$daily[[i]]=All.preds.creep[[o]]$daily[[i]]%>%
            mutate(add.crp=Eff.creep$effort.creep[match(finyear,Eff.creep$finyear)],
                   lower.CL=lower.CL-(response-response*(1-add.crp)),
                   upper.CL=upper.CL-(response-response*(1-add.crp)),
                   response=response*(1-add.crp))
        }
      }
    }
    
    #4. Compare annual indices for different scenarios   
    Plot.cpue.tweedie=function(cpuedata,CL,CxS,Yvar,add.lines,firstyear,ADD.nomnl)    #plot cpues
    {
      if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
      if(length(cpuedata)<=3)tc=seq(-.5*0.25,.5*0.25,length.out=length(cpuedata))
      
      Yrs=as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4))
      Tops=c(unlist(lapply(cpuedata, `[`, "upper.CL")),ADD.nomnl$response)
      Tops=Tops[!is.infinite(Tops)]
      ymax=max(Tops)
      Quant=quantile(Tops,probs=c(.9,1))
      if(diff(Quant)>3) ymax=quantile(Tops,probs=.99)
      plot(Yrs,Yrs,ylim=c(0,ymax),xlim=c(firstyear,max(Yrs)),ylab="",xlab="",
           col="transparent",cex.axis=1.25)
      for(l in 1:length(cpuedata))
      {
        if(!is.null(cpuedata[[l]]))
        {
          aaa=cpuedata[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
          msn=Yrs[which(!Yrs%in%c(aaa$finyear,aaa.daily$finyear))]    
          if(length(msn)>0)
          {
            ad=aaa[length(msn),]
            ad[,]=NA
            ad$finyear=msn
            aaa=rbind(aaa,ad)
            aaa=aaa[order(aaa$finyear),]
          }
          with(aaa,
               {
                 if(add.lines=="NO") points(finyear+tc[l], response, pch=19, lty=2, col=CL[l],cex=CxS)
                 if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=19, lty=2, col=CL[l],cex=CxS)
                 arrows(x0=finyear+tc[l], y0=lower.CL, 
                        x1=finyear+tc[l], y1=upper.CL, 
                        code=3, angle=90, length=0.05, col=CL[l])
               })
        }
      }
      if(!is.null(ADD.nomnl))
      {
        with(ADD.nomnl,points(as.numeric(substr(finyear,1,4))+.1,response,pch=19,col=rgb(.1,.1,.1,alpha=.2),cex=CxS)) 
      }
    }
    All.preds.creep$`All vessels & blocks`$monthly['Sandbar shark']=list(NULL) #cannot estimate all vessels
    fn.fig("Appendix_Sensitivity",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in 1:length(Tar.sp))
    {
      aaa.daily=Stand.out.daily[[Tar.sp[s]]]$DATA%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
      
      #Monthly
      Da=list(All.preds.creep$`Base case`$monthly[[s]],
              All.preds.creep$`All vessels & blocks`$monthly[[s]],
              All.preds.creep$`No efficiency`$monthly[[s]])
      names(Da)=sens$Scenario
      Plot.cpue.tweedie(cpuedata=Da,CL=c("black","forestgreen","red"),
                        CxS=1.15,Yvar="finyear",
                        add.lines="YES",firstyear=1975,
                        ADD.nomnl=NULL)
      if(s==1) mtext("Monthly",3,cex=1.5)
      if(s==4)     legend('topleft',sens$Scenario,bty='n',col=c("black","forestgreen","red"),
                          pch=19,cex=1.25)
      
      #Daily
      Da.daily=list(All.preds.creep$`Base case`$daily[[s]],
                    All.preds.creep$`All vessels & blocks`$daily[[s]],
                    All.preds.creep$`No efficiency`$daily[[s]])
      names(Da.daily)=sens$Scenario
      Plot.cpue.tweedie(cpuedata=Da.daily,CL=c("black","forestgreen","red"),
                        CxS=1.15,Yvar="finyear",
                        add.lines="YES",firstyear=2006,
                        ADD.nomnl=NULL)
      legend("topright",Nms.sp[Tar.sp[s]],bty='n',cex=1.5)
      if(s==1) mtext("Daily",3,cex=1.5)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.35,outer=T)
    mtext("CPUE (kg/km gillnet hour)",side=2,line=-0.75,font=1,las=0,cex=1.35,outer=T)
    dev.off()
    
  }
  
  stopCluster(cl)
}  


#-- Deviance explained  ----------------------------------------------
if(Model.run=="First")  #takes 4 mins
{
  if(Use.Delta)  
  {
    Anova.and.Dev.exp=function(MOD,SP,type,gam.extra)   #function for extracting term significance and deviance explained
    {
      #Anovas
      if(class(MOD)[1]=='glm')
      {
        Anova.tab=anova(MOD, test = "Chisq")
        n=2:length(Anova.tab$Deviance)
        Term.dev.exp=100*(Anova.tab$Deviance[n]/MOD$null.deviance)
        names(Term.dev.exp)=rownames(Anova.tab)[n]
        Dev.exp=sum(Term.dev.exp)
        ANOVA=as.data.frame.matrix(Anova.tab)
        Term=data.frame(Percent.dev.exp=Term.dev.exp)
        Table=ANOVA[-1,match(c("Df","Pr(>Chi)"),names(ANOVA))]
        Term=Term[match(rownames(Term), rownames(Table)),]
        Table=cbind(Table,Term)
        names(Table)[match("Term",names(Table))]="Percent.dev.exp"
        Table$term=rownames(ANOVA)[2:nrow(ANOVA)]
        Table=Table%>%select('term','Df','Pr(>Chi)','Percent.dev.exp')
        All=Table[1,]
        All[,1:ncol(All)]=NA
        All$term="model"
        All$Percent.dev.exp=round(Dev.exp,3)
        Table=rbind(Table,All)
        vars <- c(df = "Df", 'p-value' ="Pr(>Chi)")
        Table= Table %>% mutate_at(c("Pr(>Chi)","Percent.dev.exp"), round, 3) %>%
          rename(!!vars)
      }
      
      if(class(MOD)[1]=="gam")
      {
        Anova.tab=anova(MOD)
        ANOVA=as.data.frame(Anova.tab$pTerms.table[,-2])
        s.mat=as.data.frame(Anova.tab$s.table)
        s.mat=s.mat[,-(2:3)]
        names(s.mat)=names(ANOVA)
        ANOVA=rbind(ANOVA,s.mat)
        ANOVA$term=rownames(ANOVA)
        gamo=gam.extra
        for(l in 2:length(gam.extra)) gamo[l]=gam.extra[l]-gam.extra[l-1]
        gamo=100*gamo
        gam.dev.exp=data.frame(Percent.dev.exp=gamo,term=names(gam.extra))
        gam.dev.exp$term=str_remove(gam.dev.exp$term, ", k = 6")
        gam.dev.exp$term=str_remove(gam.dev.exp$term, " ")
        
        ANOVA=ANOVA%>%left_join(gam.dev.exp,by="term")
        
        ANOVA = ANOVA %>% select(term, df, 'p-value', Percent.dev.exp)
        model=ANOVA[1,]
        model[,]=NA
        model$Percent.dev.exp=sum(gamo)
        model$term='model'
        ANOVA=rbind(ANOVA,model)
        Table= ANOVA %>% mutate_at(c("p-value","Percent.dev.exp"), round, 3)  %>%
          mutate_at(c("df"), round, 0)
        
      }
      Table$"p-value"=ifelse(Table$"p-value"<0.001,"<0.001",Table$"p-value")
      dummy=Table[1:2,]
      dummy[,]=NA
      dummy$term=c(type,SP)
      Table=rbind(dummy,Table)
      Table[is.na(Table)] <- ""
      rownames(Table)=NULL
      return(Table)
    }
    
    #calculate gam deviance explained by each term
    gam.list=Stand.out.daily
    
    system.time({dummy1=foreach(s=nnn,.packages=c('mgcv')) %dopar%
      {
        if(!is.null(gam.list[[s]]))
        {
          if(s %in% Tar.sp)
          {
            ALLvars.gam.bi=c(1,labels(terms(Best.Model.daily.gam_delta[[s]]$bi)))
            ALLvars.gam.pos=c(1,labels(terms(Best.Model.daily.gam_delta[[s]]$pos)))
            
            dev.exp.bi=rep(NA,length(ALLvars.gam.bi))
            names(dev.exp.bi)=ALLvars.gam.bi
            dev.exp.pos=rep(NA,length(ALLvars.gam.pos))
            names(dev.exp.pos)=ALLvars.gam.pos
            
            for(g in 1:length(ALLvars.gam.bi))
            {
              if(g==1) added.bit=paste(ALLvars.gam.bi[g],collapse="+") else
                added.bit=paste(ALLvars.gam.bi[1:g],collapse="+")
              Formula.gam=as.formula(paste('catch.target',"~",paste(added.bit,'offset(LN.effort)',sep="+")))
              res.gam <-gam(Formula.gam,data=gam.list[[s]]$Bi$DATA,family='binomial',method="REML")
              dev.exp.bi[g]=summary(res.gam)$dev.expl
            }
            dev.exp.bi=dev.exp.bi[-1]  
            
            for(g in 1:length(ALLvars.gam.pos))
            {
              if(g==1) added.bit=paste(ALLvars.gam.pos[g],collapse="+") else
                added.bit=paste(ALLvars.gam.pos[1:g],collapse="+")
              Formula.gam=as.formula(paste('LNcpue',"~",added.bit))
              res.gam <-gam(Formula.gam,data=gam.list[[s]]$Pos$DATA,method="REML")
              dev.exp.pos[g]=summary(res.gam)$dev.expl
            }
            dev.exp.pos=dev.exp.pos[-1]  
            return(list(dev.exp.bi=dev.exp.bi,dev.exp.pos=dev.exp.pos))
          }
          if(!s %in% Tar.sp)
          {
            ALLvars.gam=all.vars(Best.Model.daily.gam[[s]])[-1]
            ALLvars.gam=c(1,"s(long10.corner,lat10.corner)",ALLvars.gam[-match(c('long10.corner','lat10.corner'),ALLvars.gam)])
            dev.exp=rep(NA,length(ALLvars.gam))
            names(dev.exp)=ALLvars.gam
            dev.exp.BI=dev.exp
            for(g in 1:length(ALLvars.gam))
            {
              if(g==1) added.bit=paste(ALLvars.gam[g],collapse="+") else
                added.bit=paste(ALLvars.gam[1:g],collapse="+")
              
              Formula.bi.gam=as.formula(paste('catch.pos',"~",paste(added.bit,"offset(LNeffort)",sep="+")))
              res.gam_bi <-gam(Formula.bi.gam,data=gam.list[[s]]$DATA_bi, family="binomial",method="REML")
              
              Formula.gam=as.formula(paste('LNcpue',"~",added.bit))
              res.gam <-gam(Formula.gam,data=gam.list[[s]]$DATA,method="REML")
              
              dev.exp.BI[g]=summary(res.gam_bi)$dev.expl
              dev.exp[g]=summary(res.gam)$dev.expl
            }
            dev.exp.BI=dev.exp.BI[-1]
            dev.exp=dev.exp[-1] 
            return(list(Bi=dev.exp.BI,Pos=dev.exp))
          }
        }
      }
    })
    stopCluster(cl)
    names(dummy1)=names(gam.list)
    gam.list=dummy1
    
    Dev.exp=vector('list',length(SP.list))  
    names(Dev.exp)=names(SP.list)
    Dev.exp.daily=Dev.exp.bi=Dev.exp.daily.bi=Dev.exp
    system.time({for(s in nnn)
    {
      #monthly
      if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out[[s]]$res)))
      {
        if(s%in%Tar.sp)
        {
          MOD.pos=Stand.out[[s]]$Pos$res
          MOD.bi=Stand.out[[s]]$Bi$res
        }else
        {
          MOD.pos=Stand.out[[s]]$res
          MOD.bi=Stand.out[[s]]$res_bi
        }
        Dev.exp[[s]]=Anova.and.Dev.exp(MOD=MOD.pos,SP=Nms.sp[s],type="Monthly.pos",gam.extra=NULL)
        Dev.exp.bi[[s]]=Anova.and.Dev.exp(MOD=MOD.bi,SP=Nms.sp[s],type="Monthly.bi",gam.extra=NULL)
      }
      
      #daily
      if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out.daily[[s]]$res.gam)))
      {
        if(s%in%Tar.sp)
        {
          MOD.pos=Stand.out.daily[[s]]$Pos$res.gam
          MOD.bi=Stand.out.daily[[s]]$Bi$res.gam
          Gm.ext=gam.list[[s]]$dev.exp.pos
          Gm.ext.bi=gam.list[[s]]$dev.exp.bi
        }else
        {
          MOD.pos=Stand.out.daily[[s]]$res.gam
          MOD.bi=Stand.out.daily[[s]]$res.gam_bi
          Gm.ext=gam.list[[s]]$Pos
          Gm.ext.bi=gam.list[[s]]$Bi
        }
        Dev.exp.daily[[s]]=Anova.and.Dev.exp(MOD=MOD.pos,SP=Nms.sp[s],type="Daily.pos",gam.extra=Gm.ext)
        Dev.exp.daily.bi[[s]]=Anova.and.Dev.exp(MOD=MOD.bi,SP=Nms.sp[s],type="Daily.bi",gam.extra=Gm.ext.bi)
        
      }
      
    }})  
    
    Tab.Dev.Exp=rbind(do.call(rbind,Dev.exp),do.call(rbind,Dev.exp.bi),
                      do.call(rbind,Dev.exp.daily),do.call(rbind,Dev.exp.daily.bi))
    rownames(Tab.Dev.Exp)=NULL
  }
  
  if(Use.Tweedie)     #takes 7 mins
  {
    Anova.and.Dev.exp=function(MOD,SP,Dev.Exp)   #function for extracting term significance and deviance explained
    {
      Anova.tab=anova(MOD)
      ANOVA=as.data.frame(Anova.tab$pTerms.table)
      n.an=1:nrow(ANOVA)
      s.mat=as.data.frame(Anova.tab$s.table)
      Para.terms=rownames(ANOVA)
      Smooth.terms=rownames(s.mat)
      
      misn.col=ncol(s.mat)-ncol(ANOVA)
      if(misn.col>0)
      {
        dummy=as.data.frame(ANOVA[,1:misn.col])
        names(dummy)=NULL
        dummy[,]=NA
        ANOVA=cbind(ANOVA,dummy)
        names(ANOVA)=names(s.mat)
      }
      
      ANOVA=rbind(ANOVA,rep(NA,ncol(ANOVA)),rep(NA,ncol(ANOVA)))
      ANOVA=rbind(ANOVA,s.mat) %>% mutate_all(round, 3)
      ANOVA$term=c(Para.terms,NA,NA,Smooth.terms)
      
      Dev.exp.mod=round(Dev.Exp$Dev.exp[nrow(Dev.Exp)],3)
      Dev.Exp$Term <- str_remove_all(Dev.Exp$Term, paste(c(", bs = \"re\"", ", k = 12, bs = \"cc\""), collapse = "|"))
      Dev.Exp=Dev.Exp%>%
        mutate(Lag=lag(Dev.exp,1,default = 0),
               Percent.deviance.explained=round(Dev.exp-Lag,3))%>%
        dplyr::select(Term,Percent.deviance.explained)
      
      Dev.Exp$Term=gsub(" ", "", Dev.Exp$Term, fixed = TRUE)
      
      ANOVA=ANOVA%>%
        left_join(Dev.Exp,by=c("term"="Term"))%>%
        dplyr::select(term,edf,Ref.df,'F','p-value',Percent.deviance.explained)
      ANOVA$"p-value"=ifelse(ANOVA$"p-value"<0.001,"<0.001",ANOVA$"p-value")
      ANOVA$F[n.an]=ifelse(ANOVA$F[n.an]<0.001,"<0.001",ANOVA$F[n.an])
      Which.NA=which(is.na(ANOVA$term))
      ANOVA[Which.NA[1],1]="Smooth terms"
      ANOVA[Which.NA[2],]=colnames(ANOVA)
      colnames(ANOVA)=c('term','df','F','p-value',NA,'Percent.deviance.explained')
      
      model=ANOVA[1,]
      model[,]=NA
      model$Percent.deviance.explained=Dev.exp.mod
      model$term='model'
      Table=rbind(ANOVA,model)
      
      dummy=Table[1:3,]
      dummy[,]=NA
      dummy$term[1:2]=c(SP,"Parametric terms")
      dummy[3,]=colnames(Table)
      Table=rbind(dummy,Table)
      Table[is.na(Table)] <- ""
      rownames(Table)=NULL
      colnames(Table)=rep('dummy',ncol(Table))
      return(Table)
    }
    
    #monthly
    system.time({Dev.monthly=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','stringr')) %dopar%
      {
        if(!is.null(BLKS.used[[s]]))
        {
          d=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]])
          if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
          if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
          Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
          Continuous=Covariates.monthly
          colnames(d)=tolower(colnames(d))
          Terms=tolower(Terms)
          Continuous=tolower(Continuous)
          Factors=Terms[!Terms%in%Continuous]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          d <- makecategorical(Factors[Factors%in%Terms],d)
          
          Terms=c(1,labels(terms(Best.Model[[s]])))
          dev.exp=data.frame(Term=Terms,Dev.exp=NA)
          for(g in 1:length(Terms))
          {
            if(g==1) added.bit=paste(Terms[g],collapse="+") else
              added.bit=paste(Terms[1:g],collapse="+")
            Formula.gam=as.formula(paste('cpue',"~",added.bit))
            mod<-gam(Formula.gam,data=d,family='tw',method="REML")
            dev.exp$Dev.exp[g]=summary(mod)$dev.expl*100
          }
          
          AOV=Anova.and.Dev.exp(MOD=mod,SP=Nms.sp[s],Dev.Exp=dev.exp)
          
          return(AOV)
          
          rm(d,mod)
        }
      }
    })   
    
    #daily
    system.time({Dev.daily=foreach(s=Tar.sp,.packages=c('dplyr','mgcv','stringr')) %dopar%
      {
        if(!is.null(BLKS.used.daily[[s]]))
        {
          d=DATA.list.LIVEWT.c.daily[[s]]%>%
            filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]])
          Terms=Predictors_daily
          Continuous=Covariates.daily
          colnames(d)=tolower(colnames(d))
          Terms=tolower(Terms)
          Continuous=tolower(Continuous)
          Factors=Terms[!Terms%in%Continuous]
          Terms=all.vars(Best.Model.daily[[s]])[-1]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          d <- makecategorical(Factors[Factors%in%Terms],d)
          
          Terms=c(1,labels(terms(Best.Model.daily[[s]])))
          
          dev.exp=data.frame(Term=Terms,Dev.exp=NA)
          for(g in 1:length(Terms))
          {
            if(g==1) added.bit=paste(Terms[g],collapse="+") else
              added.bit=paste(Terms[1:g],collapse="+")
            Formula.gam=as.formula(paste('cpue',"~",added.bit))
            mod<-gam(Formula.gam,data=d,family='tw',method="REML")
            dev.exp$Dev.exp[g]=summary(mod)$dev.expl*100
          }
          
          AOV=Anova.and.Dev.exp(MOD=mod,SP=Nms.sp[s],Dev.Exp=dev.exp)
          
          return(AOV)
          
          rm(d,mod)
        }
      }
    })
    
    
    #add power parameter
    for(s in 1:length(Tar.sp))
    {
      #Monthly
      Dev.monthly[[s]]$power.p=summary(Stand.out[[Tar.sp[s]]]$res.gam)[11]$family$family
      
      #Daily
      Dev.daily[[s]]$power.p=summary(Stand.out.daily[[Tar.sp[s]]]$res.gam)[11]$family$family
      
    }
    
    Tab.Dev.Exp=rbind(do.call(rbind,Dev.monthly),do.call(rbind,Dev.daily))
    colnames(Tab.Dev.Exp)=NULL
    
    stopCluster(cl)
  }
  
  fn.word.table(WD=getwd(),TBL=as.matrix(Tab.Dev.Exp),Doc.nm="ANOVA_table",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
}


#-- Fit diagnostics  ---------------------------------------------- 
#par(mfrow = c(2,2))
#gam.check(gam_y)  #to see gam fit
if(Model.run=="First") 
{
  if(Use.Delta)  #positive part only
  {
    Pos.Diag.fn=function(MODEL,SPECIES,M)   #function for positive catch diagnostics
    {
      RES=MODEL$residuals   #residuals
      Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
      PRED=predict(MODEL)
      
      qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="",xlab="")
      qqline(RES, col = 'grey40',lwd=1.5,lty=2)
      mtext(SPECIES,3,outer=F,line=0.25,cex=1.3)
      if(s==1) mtext("Residuals",2,outer=F,line=2,las=3,cex=M)
      if(s==2) mtext("                        Theoretical quantiles",1,outer=F,line=1.5,cex=M)
      
      hist(Std.RES,xlim=c(-5,5),ylab="",xlab="",main="",col="grey",breaks=50)
      box()
      if(s==1) mtext("Frequency",2,outer=F,line=2.5,las=3,cex=M)
      if(s==2) mtext("                      Standardised residuals",1,outer=F,line=1.5,cex=M)
      
      plot(PRED,Std.RES,ylab="",xlab="",ylim=c(-5,5))
      abline(0,0,lwd=1.5,lty=2,col='grey40')
      if(s==1) mtext("Standardised residuals",2,outer=F,line=2,las=3,cex=M)
      if(s==2) mtext("                         Fitted values",1,outer=F,line=1.5,cex=M)
      
    }
    
    fn.fig("Appendix 6",2000, 2400)
    par(mfcol=c(2*3,4),las=1,mar=c(2,2,1.75,1),oma=c(1,2,.1,2),las=1,mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.1)
    for(s in Tar.sp) 
    {
      Pos.Diag.fn(MODEL=Stand.out[[s]]$Pos$res,SPECIES=Nms.sp[s],M=.9)
      
      #Daily
      Pos.Diag.fn(MODEL=Stand.out.daily[[s]]$Pos$res.gam,SPECIES="",M=.9)
    }
    mtext(c("Daily logbooks                                          Monthly returns     "),4,
          outer=T,las=3,line=0,cex=1.3)
    dev.off()
  }
  
  if(Use.Tweedie)
  {
    
    #Fit
    #notes: 
    # small p-values indicate that residuals are not randomly distributed. 
    # This often means there are not enough basis functions.
    
    #Q-Q plot compares the model residuals to a normal distribution. 
    # A well-fit model's residuals will be close to a straight line. 
    
    # Histogram of residuals should be a symmetrical bell shape.
    
    # Residual values should be evenly distributed around zero.
    
    # Response against fitted values; the pattern should cluster around the 1-to-1 line.
    Fit.Diag.fn=function(mod,SP)   
    {
      out=capture.output(gam.check(mod,cex.main=.9))  
      mtext(SP,4,las=3,cex=.75,line=1)
      return(out)
    }
    
    QQ.fn=function(mod) qq(mod,cex=2)
    
    fn.fig("Appendix 6",2400, 1800)
    par(mfrow=c(2,4),las=1,mar=c(1.8,1,1.25,1.2),oma=c(1.5,1.75,.1,1),las=1,
        mgp=c(2,.5,0),cex.axis=1,cex.lab=1.1)
    #monthly
    for(s in Tar.sp)
    {
      QQ.fn(mod=Stand.out[[s]]$res.gam)
      mtext(Nms.sp[s],3,cex=.95,line=0.25)
      if(s==Tar.sp[4]) mtext("Monthly returns",4,las=3,line=.5,cex=1.15)
    }
    
    #daily
    for(s in Tar.sp)
    {
      QQ.fn(mod=Stand.out.daily[[s]]$res.gam)
      if(s==Tar.sp[4])mtext("Daily logbooks",4,las=3,line=.5,cex=1.15)
    }
    mtext("Deviance residuals",2,las=3,line=0.25,cex=1.2,outer=T)
    mtext("Theoretical quantiles",1,line=0.1,cex=1.2,outer=T)
    
    dev.off()
    
    Store.fit=vector('list',length(nnn))
    names(Store.fit)=Nms.sp
    Store.fit.daily=Store.fit
    
    #monthly
    for(s in Tar.sp) Store.fit[[s]]=Fit.Diag.fn(mod=Stand.out[[s]]$res.gam,SP=Nms.sp[s])
    
    #daily
    for(s in Tar.sp) Store.fit.daily[[s]]=Fit.Diag.fn(mod=Stand.out.daily[[s]]$res.gam,SP=Nms.sp[s])
    
    
    
    Store.fit <- Store.fit[which(!sapply(Store.fit, is.null))]
    Store.fit.daily <- Store.fit.daily[which(!sapply(Store.fit.daily, is.null))]
    
    capture.output(Store.fit, file = "Appendix 6_monthly.txt")
    capture.output(Store.fit.daily, file = "Appendix 6_daily.txt")
    
    #colinearity
    #concurvity(mod,full=TRUE)
    #concurvity(mod,full=FALSE)
    
    #autocorrelation
    #layout(matrix(1:2, ncol = 2))
    #acf(resid(mod), lag.max = 36, main = "ACF")
    #pacf(resid(mod), lag.max = 36, main = "pACF")
    
    
  }
  
}


#-- Plot base case, unstandardised and nominal  ----------------------------------------------  

#Extract comparable unstandardised cpue
Unstandardised=vector('list',length(SP.list)) 
names(Unstandardised)=names(SP.list)
Unstandardised.daily=Unstandardised

if(Use.Delta)
{
  for(s in nnn)   
  {
    if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out[[s]]$DATA)))
    {
      DAT=DATA.list.LIVEWT.c[[s]]
      colnames(DAT)=tolower(colnames(DAT)) 
      DAT=DAT%>%filter(vessel%in%VES.used[[s]] & blockx%in%as.numeric(BLKS.used[[s]]))
      if(s%in%Tar.sp) DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select years with a min number of vessels
      Unstandardised[[s]]=fn.out.nominal(d=DAT,method="DLnMean")
    }
    
    if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out.daily[[s]]$DATA)))
    {
      DAT=DATA.list.LIVEWT.c.daily[[s]]
      colnames(DAT)=tolower(colnames(DAT)) 
      DAT=DAT%>%filter(vessel%in%VES.used.daily[[s]] & blockx%in%as.numeric(BLKS.used.daily[[s]]))
      if(s%in%Tar.sp) DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select years with a min number of vessels
      Unstandardised.daily[[s]]=fn.out.nominal(d=DAT,method="DLnMean")
    }
    
  }
}

if(Use.Tweedie)   
{
  for(s in nnn)   
  {
    #Monthly
    if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out[[s]]$DATA)))
    {
      DAT=Stand.out[[s]]$DATA%>%
        mutate(effort=km.gillnet.hours.c,
               catch=catch.target)
      Unstandardised[[s]]=fn.out.nominal(d=DAT,method="Mean")
    }
    
    #Daily
    if(s%in%Tar.sp | (!s%in%Tar.sp & !is.null(Stand.out.daily[[s]]$DATA)))
    {
      DAT=Stand.out.daily[[s]]$DATA%>%
        mutate(effort=km.gillnet.hours.c,
               catch=catch.target)
      Unstandardised.daily[[s]]=fn.out.nominal(d=DAT,method="Mean")
    }
    
  }
}

#Rename unstandardised for plotting purposes ----------------------------------------------
for(s in nnn)
{
  if(!is.null(Unstandardised[[s]]))
  {
    Unstandardised[[s]]=Unstandardised[[s]]%>%
      rename(response=mean,
             lower.CL=lowCL,
             upper.CL= uppCL) 
  }
  
  if(!is.null(Unstandardised.daily[[s]]))
  {
    Unstandardised.daily[[s]]=Unstandardised.daily[[s]]%>% 
      rename(response=mean,
             lower.CL=lowCL,
             upper.CL=uppCL)
    
  }
}

#Compare glm and gam spatial predictions ----------------------------------------------
if(compare.glm.gam_spatial=="YES")
{
  AOV.tab=function(Anova.tab,GLM)
  {
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/GLM$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    Dev.exp=sum(Term.dev.exp)
    
    ANOVA=as.data.frame.matrix(Anova.tab)
    Term=data.frame(Percent.dev.exp=Term.dev.exp)
    Table=ANOVA[-1,match(c("Deviance","Pr(>Chi)"),names(ANOVA))]
    Term=Term[match(rownames(Term), rownames(Table)),]
    Table=cbind(Table,Term)
    names(Table)[match("Term",names(Table))]="Percent.dev.exp"
    Table= Table %>% mutate_at(c("Deviance","Percent.dev.exp"), round, 3) 
    All=Table[1,]
    rownames(Table)=rownames(ANOVA)[2:nrow(ANOVA)]
    rownames(All)="Model"
    All[,1:ncol(All)]=NA
    All$Percent.dev.exp=round(Dev.exp,3)
    return(rbind(Table,All))
  }
  fn.compare.glm.gam=function(d)
  {
    Formula=formula(LNcpue ~ finyear + vessel + blockx )
    Formula_10=formula(LNcpue ~ finyear + vessel + block10 )
    Formula.gam=formula(LNcpue ~ finyear + vessel + s(long10.corner, lat10.corner) )
    
    id.fctr=which(PREDS%in%Categorical)
    d=makecategorical(PREDS[id.fctr],d) 
    d$LNcpue=log(d[,match(Response,names(d))]/d[,match(efrt,names(d))])
    id.cov=which(PREDS%in%Covariates)
    d=mutate_at(d, setNames(c(PREDS[id.cov],efrt), paste0("LN", c(PREDS[id.cov],efrt),sep="")), log)
    
    res <- glm(Formula,data=d)
    res_10 <- glm(Formula_10,data=d)
    res.gam <-gam(Formula.gam,data=d,method="REML")
    
    Tab.aov=AOV.tab(anova(res, test = "Chisq"),res)
    row.names(Tab.aov)=paste("glm",row.names(Tab.aov),sep='...')
    Tab.aov_10=AOV.tab(anova(res_10, test = "Chisq"),res_10)
    row.names(Tab.aov_10)=paste("glm_10",row.names(Tab.aov_10),sep='...')
    
    aa=summary(res.gam)
    Tab.aov.gam=data.frame(Deviance=NA,'Pr(>Chi)'=NA,Percent.dev.exp=round(100*aa$dev.expl,2))
    row.names(Tab.aov.gam)=paste("gam","Model",sep='...')
    colnames( Tab.aov.gam)=colnames(Tab.aov)
    
    grid.table(rbind(Tab.aov,Tab.aov_10,Tab.aov.gam))
    
    
    #finyear
    BLKr=sort(table(d$blockx))
    BLKr=names(BLKr[length(BLKr)])
    BLKr10=sort(table(d$block10))
    BLKr10=names(BLKr10[length(BLKr10)])
    long10.corner=as.numeric(d[d$block10==BLKr10,match('long10.corner',names(d))][1])
    lat10.corner=as.numeric(d[d$block10==BLKr10,match('lat10.corner',names(d))][1])
    
    VSl= sort(table(d$vessel))
    VSl=names(VSl[length(VSl)])
    
    new.block=data.frame(finyear=factor(levels(d$finyear)),
                         vessel=factor(VSl,levels(d$vessel)),
                         blockx=factor(BLKr,levels(d$blockx)))
    new.block$cpue=predict(res,newdata=new.block,type='response')
    
    new.block10=data.frame(finyear=factor(levels(d$finyear)),
                           vessel=factor(VSl,levels(d$vessel)),
                           block10=factor(BLKr10,levels(d$block10)))
    new.block10$cpue=predict(res_10,newdata=new.block10,type='response')
    
    new.gam=data.frame(finyear=factor(levels(d$finyear)),
                       vessel=factor(VSl,levels(d$vessel)),
                       long10.corner=long10.corner,
                       lat10.corner=lat10.corner)
    new.gam$cpue=predict(res.gam,newdata=new.gam,type='response')
    
    new.block$cpue=new.block$cpue/mean(new.block$cpue)
    new.block10$cpue=new.block10$cpue/mean(new.block10$cpue)
    new.gam$cpue=new.gam$cpue/mean(new.gam$cpue)
    
    yrs=new.gam$finyear
    plot(1:length(yrs),new.block$cpue,xaxt='n',xlab="financial year",ylab="normalised cpue",
         cex.lab=1.5,pch=19,cex=1.5,type='o',
         ylim=c(min(c(new.block$cpue,new.block10$cpue,new.gam$cpue)),max(c(new.block$cpue,new.block10$cpue,new.gam$cpue))))
    axis(1,1:length(yrs),yrs)
    points((1:length(yrs))+.15,new.block10$cpue,pch=19,col='cyan3',cex=1.5)
    lines((1:length(yrs))+.15,new.block10$cpue,pch=19,col='cyan3')
    points((1:length(yrs))-.15,new.gam$cpue,pch=19,col='grey60',cex=1.5)
    lines((1:length(yrs))-.15,new.gam$cpue,pch=19,col='grey60')
    legend("topright",c("block","block10","gam"),bty='n',pch=19,col=c("black","cyan3","grey60"),cex=1.5)
    #space
    par(mfcol=c(3,1),mar=c(2,2,2,2),oma=c(1,1,.1,.1),mgp=c(1.5,.6,0))
    
    FINYr=sort(table(d$finyear))
    FINYr=names(FINYr[length(FINYr)])
    VSl= sort(table(d$vessel))
    VSl=names(VSl[length(VSl)])
    
    new.block=data.frame(finyear=factor(FINYr,levels(d$finyear)),
                         vessel=factor(VSl,levels(d$vessel)),
                         blockx=factor(levels(d$blockx)))
    new.block$cpue=predict(res,newdata=new.block,type='response')
    new.block=new.block%>%mutate( Lat=-as.numeric(substr(blockx,1,2)),
                                  Long=100+as.numeric(substr(blockx,3,4)))%>%
      select(-c(blockx,finyear,vessel))
    YLIM=c(min(new.block$Lat),max(new.block$Lat))
    XLIM=c(min(new.block$Long),max(new.block$Long))
    seq.lat=seq(YLIM[1],YLIM[2])
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.block$Lat))]
    seq.lon=seq(XLIM[1],XLIM[2])
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.block$Long))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(Lat=seq.lat, Long=seq.lon)
      new.block=combo%>%left_join(new.block,by=c("Lat","Long"))
    }
    new.block=new.block%>%spread(Lat,cpue)
    Lon=as.numeric(new.block$Long)
    new.block=as.matrix(new.block[,-1]) 
    LaT=as.numeric(colnames(new.block))
    brk<- quantile( c(new.block),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.block, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Blockx",bty='n',cex=1.5)
    
    
    
    new.block10=data.frame(finyear=factor(FINYr,levels(d$finyear)),
                           vessel=factor(VSl,levels(d$vessel)),
                           block10=factor(levels(d$block10)))
    new.block10$cpue=predict(res_10,newdata=new.block10,type='response')
    dummy=subset(d,select=c(block10,lat10.corner,long10.corner)) %>%
      mutate(block10=as.character(block10)) %>%
      distinct(block10,.keep_all =T)
    new.block10=new.block10%>%mutate(block10=as.character(block10))%>%
      left_join(dummy,by="block10") %>%
      select(-c(block10,finyear,vessel)) %>%
      mutate(lat10.corner=round(lat10.corner,2),
             long10.corner=round(long10.corner,2))
    YLIM=c(min(new.block10$lat10.corner),max(new.block10$lat10.corner))
    XLIM=c(min(new.block10$long10.corner),max(new.block10$long10.corner))
    #seq(112,129,length.out = 7+6*16)
    seq.lat=c(-26.83, -26.67, -26.50, -26.33, -26.17, -26.00,
              -27.83, -27.67, -27.50, -27.33, -27.17, -27.00,
              -28.83, -28.67, -28.50, -28.33, -28.17, -28.00,
              -29.83, -29.67, -29.50, -29.33, -29.17, -29.00,
              -30.83, -30.67, -30.50, -30.33, -30.17, -30.00,
              -31.83, -31.67, -31.50, -31.33, -31.17, -31.00,
              -32.83, -32.67, -32.50, -32.33, -32.17, -32.00,
              -33.83, -33.67, -33.50, -33.33, -33.17, -33.00,
              -34.83, -34.67, -34.50, -34.33, -34.17, -34.00,
              -35.83, -35.67, -35.50, -35.33, -35.17, -35.00)
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.block10$lat10.corner))]
    seq.lon=c(113.83, 113.67, 113.50, 113.33, 113.17, 113.00,
              114.83, 114.67, 114.50, 114.33, 114.17, 114.00,
              115.83, 115.67, 115.50, 115.33, 115.17, 115.00,
              116.83, 116.67, 116.50, 116.33, 116.17, 116.00,
              117.83, 117.67, 117.50, 117.33, 117.17, 117.00,
              118.83, 118.67, 118.50, 118.33, 118.17, 118.00,
              119.83, 119.67, 119.50, 119.33, 119.17, 119.00,
              120.83, 120.67, 120.50, 120.33, 120.17, 120.00,
              121.83, 121.67, 121.50, 121.33, 121.17, 121.00,
              122.83, 122.67, 122.50, 122.33, 122.17, 122.00,
              123.83, 123.67, 123.50, 123.33, 123.17, 123.00,
              124.83, 124.67, 124.50, 124.33, 124.17, 124.00,
              125.83, 125.67, 125.50, 125.33, 125.17, 125.00,
              126.83, 126.67, 126.50, 126.33, 126.17, 126.00,
              127.83, 127.67, 127.50, 127.33, 127.17, 127.00,
              128.83, 128.67, 128.50, 128.33, 128.17, 128.00,
              129.83, 129.67, 129.50, 129.33, 129.17, 129.00)
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.block10$long10.corner))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(lat10.corner=seq.lat, long10.corner=seq.lon)
      new.block10=combo%>%left_join(new.block10,by=c("lat10.corner","long10.corner"))
    }
    new.block10=new.block10%>%spread(lat10.corner,cpue)
    Lon=as.numeric(new.block10$long10.corner)
    new.block10=as.matrix(new.block10[,-1]) 
    LaT=as.numeric(colnames(new.block10))
    brk<- quantile( c(new.block10),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.block10, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Block10",bty='n',cex=1.5)
    
    
    
    ##option using lsmeans
    # dummy=subset(d,select=c(block10,lat10.corner,long10.corner)) %>%
    #   mutate(block10=as.character(block10)) %>%
    #   distinct(block10,.keep_all =T)%>%
    #   arrange(block10)
    # Long.seq=dummy$long10.corner
    # Lat.seq=dummy$lat10.corner
    # a=summary(ref_grid(res.gam, at = list(finyear=new.finyr,lat10.corner =Lat.seq ,long10.corner = Long.seq)))
    
    new.gam=d%>%select(block10,lat10.corner,long10.corner)%>%
      distinct(block10,.keep_all=T)%>%
      mutate(finyear=factor(FINYr,levels(d$finyear)),
             vessel=factor(VSl,levels(d$vessel)))
    new.gam$cpue=predict(res.gam,newdata=new.gam,type='response')
    new.gam=new.gam%>%select(-c(block10,finyear,vessel))%>%
      mutate(lat10.corner=round(lat10.corner,2),
             long10.corner=round(long10.corner,2))
    
    misn.lat=seq.lat[which(!seq.lat%in%unique(new.gam$lat10.corner))]
    misn.lon=seq.lon[which(!seq.lon%in%unique(new.gam$long10.corner))]
    if(length(misn.lat)>0 | length(misn.lon)>0)
    {
      combo=expand.grid(lat10.corner=seq.lat, long10.corner=seq.lon)
      new.gam=combo%>%left_join(new.gam,by=c("lat10.corner","long10.corner"))
    }
    new.gam=new.gam%>%spread(lat10.corner,cpue)
    Lon=as.numeric(new.gam$long10.corner)
    new.gam=as.matrix(new.gam[,-1]) 
    LaT=as.numeric(colnames(new.gam))
    brk<- quantile( c(new.gam),probs=seq(0,1,.1),na.rm=T)
    YLIM[1]=YLIM[1]-0.5
    YLIM[2]=YLIM[2]+0.5
    XLIM[1]=XLIM[1]-0.5
    XLIM[2]=XLIM[2]+0.5
    image.plot(Lon,LaT,new.gam, breaks=brk, col=rev(heat.colors(length(brk)-1)), 
               lab.breaks=names(brk),ylim=YLIM,xlim=XLIM,ylab="",xlab="")
    legend('topright',"Gam",bty='n',cex=1.5)
    
    mtext("Longitude",1,outer=T)
    mtext("Latitude",2,outer=T)
    
    
    #explore gam
    xlpr.gam="NO"
    if(xlpr.gam=="YES")
    {
      par(mfrow = c(2,2))
      gam.check(res.gam)
      #small p-values indicate that residuals are not randomly distributed. 
      # This often means there are not enough basis functions
      
      
      par(mfrow = c(1,2))
      plot(res.gam, residuals = TRUE, pch = 1)
      plot(res.gam, residuals = TRUE, pch = 1, scheme = 2)
      #In this plot the axes represent values of our predictor variables, x1 and x2. 
      #The interior is a topographic map of predicted values. 
      #The contour lines represent points of equal predicted values, and they are labeled. 
      #The dotted lines show uncertainty in prediction; they represent how contour
      # lines would move if predictions were one standard error higher or lower.
      
      vis.gam(x = res.gam,                # GAM object
              view = c("long10.corner", "lat10.corner"),   # variables
              plot.type = "contour", too.far = 0.05)    # kind of plot 
      
    }
    
    
  }
  system.time({for(s in Tar.sp)
  {
    DAT=subset(Store_nom_cpues_daily[[s]]$QL_dat,vessel%in%VES.used.daily[[s]])
    DAT=subset(DAT,blockx%in%BLKS.used.daily[[s]])
    DAT=subset(DAT,finyear%in%fn.sel.yrs.used.glm(DAT)) #select min vessels per year
    DAT=DAT%>% mutate(mean.depth=10*round(mean.depth/10),
                      nlines.c=ifelse(nlines.c>2,'>2',nlines.c),
                      mesh=ifelse(!mesh%in%c(165,178),'other',mesh))
    
    Response="catch.target"
    RESPNS="LNcpue"
    PREDS=c("finyear","vessel","block10","blockx")
    efrt="km.gillnet.hours.c"
    
    pdf(paste(getwd(),"/Compare glm and gam/Compare.glm.vs.gam_",Nms.sp[s],".pdf",sep=""))
    fn.compare.glm.gam(d=DAT)
    dev.off()
    rm(DAT)
  }
  })
}

#Predict year effect (considering log bias corr if required) ----------------------------------------------
if(Use.Qualif.level)
{
  #Target species  
  #monthly          takes 80 sec
  system.time({Pred.tar=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
      return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="finyear",Pred.type="link"))
      rm(d)
    }
  })
  #Daily              takes 70 sec
  system.time({Pred.daily.tar=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
    {
      d=Stand.out.daily[[s]]$DATA
      return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="finyear",Pred.type="link"))
      rm(d)
    }
  })
  
}
if(Use.Delta)
{
  #Other species
  Niter=1000   #MC interations
  #monthly          takes 0.3 sec per iteration
  system.time({Pred.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mvtnorm')) %dopar%
    {
      if(!is.null(Stand.out[[s]]))
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$res_bi,
                                MOD=Stand.out[[s]]$res,
                                BiData=Stand.out[[s]]$DATA_bi,
                                PosData=Stand.out[[s]]$DATA,
                                niter=Niter,
                                pred.term='finyear',
                                ALL.terms=Predictors_monthly))
      } 
    }
  })
  #Daily              takes 0.2 sec per iteration
  system.time({Pred.daily.other=foreach(s=nnn[-sort(Tar.sp)],.packages=c('dplyr','mvtnorm')) %do%
    {
      if(!is.null(Stand.out.daily[[s]]))
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$res.gam_bi,
                                MOD=Stand.out.daily[[s]]$res.gam,
                                BiData=Stand.out.daily[[s]]$DATA_bi,
                                PosData=Stand.out.daily[[s]]$DATA,
                                niter=Niter,
                                pred.term='finyear',
                                ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
      } 
    }
  })
  names(Pred.other)=names(Pred.daily.other)=names(SP.list)[nnn[-sort(Tar.sp)]]
  
  
  #Target species
  #monthly          takes 80 sec
  system.time({Pred.tar=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
    {
      return(fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$Bi$res,
                              MOD=Stand.out[[s]]$Pos$res,
                              BiData=Stand.out[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                              PosData=Stand.out[[s]]$Pos$DATA,
                              niter=Niter,
                              pred.term='finyear',
                              ALL.terms=Predictors_monthly))
    }
  })
  #Daily              takes 70 sec
  system.time({Pred.daily.tar=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
    {
      return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$Bi$res.gam,
                              MOD=Stand.out.daily[[s]]$Pos$res.gam,
                              BiData=Stand.out.daily[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                              PosData=Stand.out.daily[[s]]$Pos$DATA,
                              niter=Niter,
                              pred.term='finyear',
                              ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
    }
  })
  
  names(Pred.tar)=names(Pred.daily.tar)=names(SP.list)[Tar.sp]
  
  
  Pred=c(Pred.tar,Pred.other)
  Pred=Pred[names(SP.list)]
  Pred.daily=c(Pred.daily.tar,Pred.daily.other)
  Pred.daily=Pred.daily[names(SP.list)]
}
if(Use.Tweedie)     #takes <1 min  
{
  #Monthly
  system.time({
    Pred=foreach(s=nnn,.packages=c('dplyr','emmeans')) %do%
      {
        if(!is.null(Stand.out[[s]]))
        {
          d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
          return(pred.fun(mod=Stand.out[[s]]$res.gam,
                          biascor="NO",
                          PRED="finyear"))
          rm(d)
        }
      }
  })
  
  #Daily
  system.time({
    Pred.daily=foreach(s=nnn,.packages=c('dplyr','emmeans')) %do%
      {
        if(!is.null(Stand.out.daily[[s]]))
        {
          d=Stand.out.daily[[s]]$DATA
          return(pred.fun(mod=Stand.out.daily[[s]]$res.gam,
                          biascor="NO",
                          PRED="finyear"))
          rm(d)
        }
      }
  })
  
  names(Pred)=names(Pred.daily)=names(SP.list)
}


#Remove early effective years for sandbar and gummy shark ----------------------------------------------
Effective$`Sandbar Shark`=Effective$`Sandbar Shark`%>%
  filter(finyear%in%Pred$`Sandbar Shark`$finyear)
Effective$`Gummy Shark`=Effective$`Gummy Shark`%>%
  filter(finyear%in%Pred$`Gummy Shark`$finyear)

#Calculate CV ----------------------------------------------
for(s in nnn)
{
  #monthly
  if(!is.null(Pred[[s]]))
  {
    Pred[[s]]=Pred[[s]]%>%mutate(CV=SD/response)
    if(!is.null(Effective[[s]])) Effective[[s]]=Effective[[s]]%>%mutate(CV=se/mean)
  }
  
  #daily
  if(!is.null(Pred.daily[[s]]))
  {
    Pred.daily[[s]]=Pred.daily[[s]]%>%mutate(CV=SD/response)
    
    if(!is.null(Effective_daily[[s]])) Effective_daily[[s]]=Effective_daily[[s]]%>%mutate(CV=se/mean)
  }
}

#Apply efficiency creep  ----------------------------------------------     
Pred.creep=Pred
Pred.daily.creep=Pred.daily
Unstandardised.creep=Unstandardised               
Unstandardised.daily.creep=Unstandardised.daily
Effective.creep=Effective
Effective.daily.creep=Effective_daily
for(s in nnn)
{
  #monthly
  if(!is.null(Pred.creep[[s]]))
  {
    yrs=Unstandardised.creep[[s]]$finyear
    Pred.creep[[s]]=subset(Pred.creep[[s]],finyear%in%yrs)
    add.crp=Eff.creep$effort.creep[match(Pred.creep[[s]]$finyear,Eff.creep$finyear)]
    
    Pred.creep[[s]]=Pred.creep[[s]]%>%
      mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
             upper.CL=upper.CL-(response-response*(1-add.crp)),
             response=response*(1-add.crp))
    
    yrs=as.character(Pred.creep[[s]]$finyear)
    Unstandardised.creep[[s]]=subset(Unstandardised.creep[[s]],finyear%in%yrs)%>%
      mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
             upper.CL=upper.CL-(response-response*(1-add.crp)),
             response=response*(1-add.crp))
    if(!is.null(Effective.creep[[s]]))
    {
      add.crp=Eff.creep$effort.creep[match(Effective.creep[[s]]$finyear,Eff.creep$finyear)]
      Effective.creep[[s]]=Effective.creep[[s]]%>%
        mutate(lower.CL=lowCL-(mean-mean*(1-add.crp)),
               upper.CL=uppCL-(mean-mean*(1-add.crp)),
               response=mean*(1-add.crp))
    }
  }
  #daily
  if(!is.null(Pred.daily.creep[[s]]))
  {
    yrs=Unstandardised.daily.creep[[s]]$finyear
    Pred.daily.creep[[s]]=subset(Pred.daily.creep[[s]],finyear%in%yrs)
    add.crp=Eff.creep$effort.creep[match(Pred.daily.creep[[s]]$finyear,Eff.creep$finyear)]
    
    Pred.daily.creep[[s]]=Pred.daily.creep[[s]]%>%
      mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
             upper.CL=upper.CL-(response-response*(1-add.crp)),
             response=response*(1-add.crp))
    
    yrs=as.character(Pred.daily.creep[[s]]$finyear)
    Unstandardised.daily.creep[[s]]=subset(Unstandardised.daily.creep[[s]],finyear%in%yrs)%>%
      mutate(lower.CL=lower.CL-(response-response*(1-add.crp)),
             upper.CL=upper.CL-(response-response*(1-add.crp)),
             response=response*(1-add.crp))
    if(!is.null(Effective.daily.creep[[s]]))
    {
      Effective.daily.creep[[s]]=Effective.daily.creep[[s]]%>%
        filter(finyear%in%yrs)%>%
        mutate(lower.CL=lowCL-(mean-mean*(1-add.crp)),
               upper.CL=uppCL-(mean-mean*(1-add.crp)),
               response=mean*(1-add.crp))
    }
  }
}

#Keep only meaningful years ----------------------------------------------
for(s in nnn)
{
  if(!is.null(Unstandardised[[s]]))
  {
    Unstandardised[[s]]=subset(Unstandardised[[s]],finyear%in%Unstandardised.creep[[s]]$finyear)
  }
  
  if(!is.null(Unstandardised.daily[[s]]))
  {
    Unstandardised.daily[[s]]=subset(Unstandardised.daily[[s]],finyear%in%
                                       Unstandardised.daily.creep[[s]]$finyear)
  }
}

#Plot Target standardised and nominal (i.e. effective)  ---------------------------------------------- 
show.how='together'
CxS=1.35
CLs=c('black',"grey70")

if(Use.Delta)
{
  Plot.cpue=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar,add.lines)    #plot cpues
  {
    if(inherits(cpuedata, "list")) 
    {
      if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
      if(length(cpuedata)<=3)tc=seq(-.5*0.15,.5*0.15,length.out=length(cpuedata))
      ymax = max(unlist(lapply(cpuedata, `[`, "upper.CL")),na.rm=T)
      Yrs=as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4))
      plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
      if(COL=='color')CL=c("black","forestgreen", "red","bisque3","blue2","dodgerblue")
      if(COL=='grey') CL=gray.colors(length(cpuedata),start=0.2,end=0.65)
      for(l in 1:length(cpuedata))
      {
        aaa=cpuedata[[l]]
        aaa$finyear=as.character(aaa$finyear)
        msn=Yrs[which(!Yrs%in%as.numeric(substr(aaa$finyear,1,4)))]
        if(length(msn)>0)
        {
          ad=aaa[length(msn),]
          ad[,]=NA
          ad$finyear=msn
          aaa=rbind(aaa,ad)
          aaa=aaa[order(aaa$finyear),]
        }
        
        with(aaa,
             {
               if(add.lines=="NO") points(Yrs+tc[l], response, pch=16, lty=2, col=CL[l],cex=CxS)
               if(add.lines=="YES") points(Yrs+tc[l], response, "o", pch=16, lty=2, col=CL[l],cex=CxS)
               arrows(x0=Yrs+tc[l], y0=lower.CL, 
                      x1=Yrs+tc[l], y1=upper.CL, 
                      code=3, angle=90, length=0.05, col=CL[l])
             })
        if(ADD.LGND=="YES") legend(whereLGND,names(cpuedata),bty='n',pch=16,col=CL,cex=1.45)
      }
    }
    
    if(inherits(cpuedata, "data.frame"))
    {
      ymax = max(cpuedata$upper.CL,na.rm=T)
      Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
      plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent",cex.axis=1.25)
      CL="black"
      if(add.lines=="NO")points(Yrs, cpuedata$response, pch=16, lty=2, col=CL,cex=CxS)
      if(add.lines=="YES")points(Yrs, cpuedata$response, "o", pch=16, lty=2, col=CL,cex=CxS)
      arrows(x0=Yrs, y0=cpuedata$lower.CL, 
             x1=Yrs, y1=cpuedata$upper.CL, 
             code=3, angle=90, length=0.05, col=CL)
    }
  }
  
  if(show.how=='together')
  {
    fn.fig("Annual_Index",1800, 2400)  
    par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      Plot.cpue.delta(cpuedata=list(Standardised=Pred.creep[[s]]),
                      cpuedata.daily=list(Standardised=Pred.daily.creep[[s]]),
                      CL=CLs,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975,
                      ADD.nomnl=rbind(Effective.creep[[s]],Effective.daily.creep[[s]]))
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    legend("bottomleft",c("Standardised","Nominal"),bty='n',pch=c(16,16),
           col=c(CLs[1],rgb(.1,.1,.1,alpha=.2)),cex=1.45)
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("CPUE (kg/km gillnet hour)",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    dev.off()
  }
  if(show.how=='separate')
  {
    fn.fig("Annual_Index",2000, 2400)   
    par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      #Monthly
      Mon.dat=list(Standardised=Pred.creep[[s]],Unstandardised=Unstandardised.creep[[s]])
      LgND="NO"
      if(s==Tar.sp[1])LgND="YES"
      Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
      if(s==Tar.sp[1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      #Daily
      Daily.dat=list(Standardised=Pred.daily.creep[[s]],Unstandardised=Unstandardised.daily.creep[[s]])
      LgND="NO"
      Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
      if(s==Tar.sp[1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      mtext(Nms.sp[s],4,line=1,las=3,cex=1.5)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("CPUE (kg/ km gillnet hour)",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
    dev.off()  
  }
}

if(Use.Tweedie)
{
  Plot.cpue=function(cpuedata,cpuedata.daily,CL,CxS,Yvar,add.lines,firstyear,ADD.nomnl,Add.CI.nominal)    
  {
    if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
    if(length(cpuedata)<=3)tc=seq(-.5*0.25,.5*0.25,length.out=length(cpuedata))
    if(length(cpuedata)==1)tc=seq(-.5*0.5,.5*0.5,length.out=length(cpuedata))
    
    Yrs=c(as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4)),
          as.numeric(substr(cpuedata.daily[[1]][,match(Yvar,names(cpuedata.daily[[1]]))],1,4)))
    Tops=c(unlist(lapply(cpuedata, `[`, "upper.CL")),
           unlist(lapply(cpuedata.daily, `[`, "upper.CL")),
           unlist(lapply(ADD.nomnl, `[`, "response")))
    Tops=Tops[!is.infinite(Tops)]
    ymax=max(Tops)
    Quant=quantile(Tops,probs=c(.95,1))
    if(diff(Quant)>3) ymax=quantile(Tops,probs=.995)
    
    plot(Yrs,Yrs,ylim=c(0,ymax),xlim=c(firstyear,max(Yrs)),ylab="",xlab="",
         col="transparent",cex.axis=1.25)
    for(l in 1:length(cpuedata))
    {
      if(!is.null(cpuedata[[l]]))aaa=cpuedata[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
      if(!is.null(cpuedata.daily[[l]]))aaa.daily=cpuedata.daily[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
      
      # msn=Yrs[which(!Yrs%in%c(aaa$finyear,aaa.daily$finyear))]
      # if(length(msn)>0)
      # {
      #   ad=aaa[length(msn),]
      #   ad[,]=NA
      #   ad$finyear=msn
      #   aaa=rbind(aaa,ad)
      #   aaa=aaa[order(aaa$finyear),]
      # }
      if(!is.null(cpuedata[[l]]))with(aaa,
                                      {
                                        if(add.lines=="NO") points(finyear+tc[l], response, pch=19, lty=2, col=CL[l],cex=CxS)
                                        if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=19, lty=2, col=CL[l],cex=CxS)
                                        arrows(x0=finyear+tc[l], y0=lower.CL, 
                                               x1=finyear+tc[l], y1=upper.CL, 
                                               code=3, angle=90, length=0.05, col=CL[l])
                                      })
      if(!is.null(cpuedata.daily[[l]]))with(aaa.daily,
                                            {
                                              if(l==1)polygon(x=c(finyear[1]-.5,finyear[length(finyear)]+.5,finyear[length(finyear)]+.5,finyear[1]-.5),
                                                              y=c(0,0,ymax*.99,ymax*.99),col='grey92',border="transparent")
                                              arrows(x0=finyear+tc[l], y0=lower.CL, 
                                                     x1=finyear+tc[l], y1=upper.CL, 
                                                     code=3, angle=90, length=0.05, col=CL[l])
                                              if(add.lines=="NO") points(finyear+tc[l], response, pch=21, lty=2, col=CL[l],cex=CxS)
                                              if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=19, lty=2, col=CL[l],cex=CxS)
                                            })
    }
    if(!is.null(ADD.nomnl))
    {
      CL.b=c("white","grey70")
      tc2=c(.1,-.1)
      for(l in 1:length(ADD.nomnl))
      {
        ADD.nomnl[[l]]=ADD.nomnl[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4))) 
        with(ADD.nomnl[[l]],
             {
               points(finyear+tc2[l],response,pch=21,bg=CL.b[l],cex=CxS)
               if(Add.CI.nominal) arrows(x0=finyear+tc2[l], y0=lower.CL, 
                                         x1=finyear+tc2[l], y1=upper.CL, 
                                         code=3, angle=90, length=0.05)
             })
      }
    }
  }
  
  if(show.how=='together')
  {
    fn.fig("Figure 3.Annual_Index",1800, 2400)  
    par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      Plot.cpue(cpuedata=list(Standardised=Pred.creep[[s]]),
                cpuedata.daily=list(Standardised=Pred.daily.creep[[s]]),
                CL=CLs,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975,
                ADD.nomnl=list(nominal=rbind(Unstandardised.creep[[s]],Unstandardised.daily.creep[[s]]),
                               effective=rbind(Effective.creep[[s]],Effective.daily.creep[[s]])),
                Add.CI.nominal=FALSE)
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    legend("bottomleft",c("Standardised","Nominal","Effective"),bty='n',pch=21,
           pt.bg=c(CLs[1],"grey80","white"),cex=1.45)
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("CPUE (kg/km gillnet hour)",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    dev.off()
  }
}


#Plot Target and nominal (and effective) normalised   ---------------------------------------------- 
show.nominal.effective=FALSE
show.raw=TRUE

Pred.normlzd=Pred.creep
Pred.daily.normlzd=Pred.daily.creep
Unstandardised.normlzd=Unstandardised.creep
Unstandardised.daily.normlzd=Unstandardised.daily.creep
Effective.normlzd=Effective.creep
Effective.daily.normlzd=Effective.daily.creep
for(s in nnn)
{
  #monthly
  if(!is.null(Pred.normlzd[[s]]))
  {
    Mn=mean(Pred.normlzd[[s]]$response)
    Pred.normlzd[[s]]$response=Pred.normlzd[[s]]$response/Mn
    Pred.normlzd[[s]]$lower.CL=Pred.normlzd[[s]]$lower.CL/Mn
    Pred.normlzd[[s]]$upper.CL=Pred.normlzd[[s]]$upper.CL/Mn
    
    Mn=mean(Unstandardised.normlzd[[s]]$response)
    Unstandardised.normlzd[[s]]=Unstandardised.normlzd[[s]]%>%
      mutate(response=response/Mn,
             lower.CL=lower.CL/Mn,
             upper.CL=upper.CL/Mn)
    
    if(!is.null(Effective.normlzd[[s]]))
    {
      Mn=mean(Effective.normlzd[[s]]$response)
      Effective.normlzd[[s]]=Effective.normlzd[[s]]%>%
        mutate(response=response/Mn,
               lower.CL=lower.CL/Mn,
               upper.CL=upper.CL/Mn)
    }
    
  }
  
  #daily
  if(!is.null(Pred.daily.normlzd[[s]]))
  {
    Mn=mean(Pred.daily.normlzd[[s]]$response)
    Pred.daily.normlzd[[s]]$response=Pred.daily.normlzd[[s]]$response/Mn
    Pred.daily.normlzd[[s]]$lower.CL=Pred.daily.normlzd[[s]]$lower.CL/Mn
    Pred.daily.normlzd[[s]]$upper.CL=Pred.daily.normlzd[[s]]$upper.CL/Mn
    
    Mn=mean(Unstandardised.daily.normlzd[[s]]$response)
    Unstandardised.daily.normlzd[[s]]=Unstandardised.daily.normlzd[[s]]%>%
      mutate(response=response/Mn,
             lower.CL=lower.CL/Mn,
             upper.CL=upper.CL/Mn)
    if(!is.null(Effective.daily.normlzd[[s]]))
    {
      Mn=mean(Effective.daily.normlzd[[s]]$response)
      
      Effective.daily.normlzd[[s]]=Effective.daily.normlzd[[s]]%>%
        mutate(response=response/Mn,
               lower.CL=lower.CL/Mn,
               upper.CL=upper.CL/Mn)
    }
  }
}

if(Use.Delta)
{
  if(show.how=='together')
  {
    fn.fig("Figure 3.Annual_Index_normalised",1800, 2400)
    par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      Plot.cpue.delta(cpuedata=list(Standardised=Pred.normlzd[[s]]),
                      cpuedata.daily=list(Standardised=Pred.daily.normlzd[[s]]),
                      CL=CLs,CxS=CxS,Yvar="finyear",
                      add.lines="YES",firstyear=1975,
                      ADD.nomnl=rbind(Effective.normlzd[[s]],Effective.daily.normlzd[[s]]))
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    legend("bottomleft",c("Standardised","Nominal"),bty='n',pch=c(16,16),
           col=c(CLs[1],rgb(.1,.1,.1,alpha=.2)),cex=1.45)
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    
    dev.off()
  }   
  if(show.how=='separate')
  {
    fn.fig("Figure 3.Annual_Index_normalised",2000, 2400) 
    par(mfrow=c(4,2),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      #Monthly
      Mon.dat=list(Standardised=Pred.normlzd[[s]],Unstandardised=Unstandardised.normlzd[[s]])
      LgND="NO"
      if(s==Tar.sp[1])LgND="YES"
      Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
      if(s==Tar.sp[1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      #Daily
      Daily.dat=list(Standardised=Pred.daily.normlzd[[s]],Unstandardised=Unstandardised.daily.normlzd[[s]])
      LgND="NO"
      Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
      if(s==Tar.sp[1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      mtext(Nms.sp[s],4,line=1,las=3,cex=1.5)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
    dev.off()
  }
}

if(Use.Tweedie)     
{
  if(show.how=='together') 
  {
    if(show.nominal.effective)
    {
      fn.fig("Figure 3.Annual_Index_normalised",1800, 2400)  
      par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
      for(s in Tar.sp)
      {
        Plot.cpue(cpuedata=list(Standardised=Pred.normlzd[[s]]),
                  cpuedata.daily=list(Standardised=Pred.daily.normlzd[[s]]),
                  CL=CLs,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975,
                  ADD.nomnl=list(nominal=rbind(Unstandardised.normlzd[[s]],Unstandardised.daily.normlzd[[s]]),
                                 effective=rbind(Effective.normlzd[[s]],Effective.daily.normlzd[[s]])),
                  Add.CI.nominal=FALSE)
        legend("topright",Nms.sp[s],bty='n',cex=1.75)
      }
      legend("bottomleft",c("Standardised","Nominal","Effective"),bty='n',pch=21,
             pt.bg=c(CLs[1],"white","grey70"),cex=1.45)
      mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
      mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
      dev.off()
    }
    if(show.raw)
    {
      fn.fig("Figure 3.Annual_Index_normalised",1800, 2400)  
      par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
      for(s in Tar.sp)
      {
        nm=match(names(Pred.normlzd)[s],names(Raw.index))
        dummy.raw=Raw.index[[nm]]%>%
          mutate(lower.CL=lowCL/mean(mean),
                 upper.CL=uppCL/mean(mean),
                 response=mean/mean(mean))%>%
          rename(finyear=FINYEAR)%>%
          filter(finyear%in%Pred.normlzd[[s]]$finyear)
        dummy.raw.daily=Raw.index.daily[[nm]]%>%
          mutate(lower.CL=lowCL/mean(mean),
                 upper.CL=uppCL/mean(mean),
                 response=mean/mean(mean))%>%
          rename(finyear=FINYEAR)%>%
          filter(finyear%in%Pred.daily.normlzd[[s]]$finyear)
        
        Plot.cpue(cpuedata=list(Standardised=Pred.normlzd[[s]]),
                  cpuedata.daily=list(Standardised=Pred.daily.normlzd[[s]]),
                  CL=CLs,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975,
                  ADD.nomnl=list(raw=rbind(dummy.raw,dummy.raw.daily)),
                  Add.CI.nominal=FALSE)
        legend("topright",Nms.sp[s],bty='n',cex=1.75)
      }
      legend("bottomleft",c("Standardised","Nominal"),bty='n',pch=21,
             pt.bg=c(CLs[1],"white"),cex=1.45)
      mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
      mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
      dev.off()
    }
  }
}


# -- Sensitivity tests comparing Tweedie, Lognormal and Delta-lognormal ----------------------------------------------
if(Use.Tweedie)
{
  if(Model.run=="First") 
  {
    #1. Run standardisation
    #monthly
    system.time({Stand.out.Tweedie.sens=foreach(s=nnn,.packages=c('dplyr','mgcv')) %dopar%
      {
        if(s %in% Tar.sp)
        {
          iid=match(SpiSis[match(names(DATA.list.LIVEWT.c)[s],names(SpiSis))],names(First.year.catch))
          Firs.yr.ktch=names(First.year.catch[[iid]])
          theseyears=sort(unique(DATA.list.LIVEWT.c[[s]]$FINYEAR))
          theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
          
          d=DATA.list.LIVEWT.c[[s]]%>%filter(VESSEL%in%VES.used[[s]] & BLOCKX%in%BLKS.used[[s]] & FINYEAR%in%theseyears )
          
          #remove first years with very few positive catch observation
          if(names(DATA.list.LIVEWT.c)[s]=="Gummy Shark")d=d%>%filter(!FINYEAR%in%c('1975-76'))
          if(names(DATA.list.LIVEWT.c)[s]=="Sandbar Shark")d=d%>%filter(!FINYEAR%in%c('1986-87','1987-88','1988-89'))
          
          Terms=Predictors_monthly[!Predictors_monthly%in%c("block10")]
          Continuous=Covariates.monthly
          colnames(d)=tolower(colnames(d))
          Terms=tolower(Terms)
          Continuous=tolower(Continuous)
          Factors=Terms[!Terms%in%Continuous]
          Terms=all.vars(Best.Model[[s]])[-1]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          
          
          #NEW#
          #Fit lognormal + constant (set at 10% mean cpue, Campbell 2004)
          Mn.cpue=mean(d%>%filter(cpue>0)%>%pull(cpue))*.1
          d.logN=d%>%mutate(cpue=ifelse(cpue==0,Mn.cpue,cpue),
                            cpue=log(cpue))
          d.logN <- makecategorical(Factors[Factors%in%Terms],d.logN)
          mod.lognormal=bam(Best.Model[[s]],data=d.logN,family="gaussian",method="fREML",discrete=TRUE)
          
          #Fit Delta
          # binomial part
          Bi <- d %>%mutate(cpue=as.numeric(catch.target>0))
          Bi <- makecategorical(Factors[Factors%in%Terms],Bi)
          mod.bi<-bam(Best.Model[[s]],data=Bi,family="binomial",method="fREML",discrete=TRUE)
          
          #lognormal part
          d.pos<- d%>%
            filter(cpue>0)%>%
            mutate(cpue=log(cpue))
          d.pos <- makecategorical(Factors[Factors%in%Terms],d.pos)
          mod.pos<-bam(Best.Model[[s]],data=d.pos,family="gaussian",method="fREML",discrete=TRUE)
          
          
          return(list(d.logN=d.logN,mod.lognormal=mod.lognormal,
                      Bi=Bi,mod.bi=mod.bi,
                      d.pos=d.pos,mod.pos=mod.pos))
          
          rm(d,d.logN,Bi,d.pos,mod.pos,mod.bi,mod.lognormal)
        }
      }
    })   
    
    #daily
    system.time({Stand.out.daily.Tweedie.sens=foreach(s=nnn,.packages=c('dplyr','mgcv')) %dopar%
      {
        if(s %in% Tar.sp)
        {
          iid=match(SpiSis[match(names(DATA.list.LIVEWT.c.daily)[s],names(SpiSis))],names(First.year.catch.daily))
          Firs.yr.ktch=names(First.year.catch.daily[[iid]])
          theseyears=sort(unique(DATA.list.LIVEWT.c.daily[[s]]$FINYEAR))
          theseyears=theseyears[match(Firs.yr.ktch,theseyears):length(theseyears)]
          
          
          d=DATA.list.LIVEWT.c.daily[[s]]%>%
            filter(VESSEL%in%VES.used.daily[[s]] & BLOCKX%in%BLKS.used.daily[[s]] & FINYEAR%in%theseyears)
          Terms=Predictors_daily
          Continuous=Covariates.daily
          colnames(d)=tolower(colnames(d))
          Terms=tolower(Terms)
          Continuous=tolower(Continuous)
          Factors=Terms[!Terms%in%Continuous]
          Terms=all.vars(Best.Model.daily[[s]])[-1]
          d <- d%>%
            dplyr::select(c(catch.target,km.gillnet.hours.c,Terms))%>%
            mutate(cpue=catch.target/km.gillnet.hours.c)
          
          
          #Fit lognormal + constant (set at 10% mean cpue, Campbell 2004)
          Mn.cpue=mean(d%>%filter(cpue>0)%>%pull(cpue))*.1
          d.logN=d%>%mutate(cpue=ifelse(cpue==0,Mn.cpue,cpue),
                            cpue=log(cpue))
          d.logN <- makecategorical(Factors[Factors%in%Terms],d.logN)
          mod.lognormal=bam(Best.Model.daily[[s]],data=d.logN,family="gaussian",method="fREML",discrete=TRUE)
          
          #Fit Delta
          # binomial part
          Bi <- d %>%mutate(cpue=as.numeric(catch.target>0))
          Bi <- makecategorical(Factors[Factors%in%Terms],Bi)
          mod.bi<-bam(Best.Model.daily[[s]],data=Bi,family="binomial",method="fREML",discrete=TRUE)
          
          #lognormal part
          d.pos<- d%>%
            filter(cpue>0)%>%
            mutate(cpue=log(cpue))
          d.pos <- makecategorical(Factors[Factors%in%Terms],d.pos)
          mod.pos<-bam(Best.Model.daily[[s]],data=d.pos,family="gaussian",method="fREML",discrete=TRUE)
          
          
          return(list(d.logN=d.logN,mod.lognormal=mod.lognormal,
                      Bi=Bi,mod.bi=mod.bi,
                      d.pos=d.pos,mod.pos=mod.pos))
          
          rm(d,d.logN,Bi,d.pos,mod.pos,mod.bi,mod.lognormal)
        }
      }
    })
    
    names(Stand.out.Tweedie.sens)=names(Stand.out.daily.Tweedie.sens)=names(SP.list)
    
    
    #2. Predict annual effect
    #Monthly
    system.time({
      Pred.Tweedie.sens=foreach(s=nnn,.packages=c('dplyr','emmeans')) %do%
        {
          if(!is.null(Stand.out.Tweedie.sens[[s]]))
          {
            #Lognormal
            d.logN=Stand.out.Tweedie.sens[[s]]$d.logN   #note: need data as global for ref_grid
            Lognormal=pred.fun(mod=Stand.out.Tweedie.sens[[s]]$mod.lognormal,biascor="YES",PRED="finyear")%>%
              dplyr::select(finyear,response)
            
            #Delta
            Bi=Stand.out.Tweedie.sens[[s]]$Bi   
            Bi.pred=pred.fun(mod=Stand.out.Tweedie.sens[[s]]$mod.bi,biascor="NO",PRED="finyear")
            
            d.pos=Stand.out.Tweedie.sens[[s]]$d.pos   
            Pos.pred=pred.fun(mod=Stand.out.Tweedie.sens[[s]]$mod.pos,biascor="YES",PRED="finyear")
            
            Delta=full_join(Bi.pred%>%dplyr::select(finyear,prob),
                            Pos.pred%>%dplyr::select(finyear,response),by='finyear')%>%
              mutate(response=response*prob)%>%
              dplyr::select(-prob)
            
            
            return(list(Lognormal=Lognormal,Delta=Delta))
            rm(d.logN,Bi,d.pos)
          }
        }
    })
    
    #Daily   
    system.time({
      Pred.daily.Tweedie.sens=foreach(s=nnn,.packages=c('dplyr','emmeans')) %do%
        {
          if(!is.null(Stand.out.daily.Tweedie.sens[[s]]))
          {
            #Lognormal
            d.logN=Stand.out.daily.Tweedie.sens[[s]]$d.logN   #note: need data as global for ref_grid
            Lognormal=pred.fun(mod=Stand.out.daily.Tweedie.sens[[s]]$mod.lognormal,biascor="YES",PRED="finyear")%>%
              dplyr::select(finyear,response)
            
            #Delta
            Bi=Stand.out.daily.Tweedie.sens[[s]]$Bi   
            Bi.pred=pred.fun(mod=Stand.out.daily.Tweedie.sens[[s]]$mod.bi,biascor="NO",PRED="finyear")
            
            d.pos=Stand.out.daily.Tweedie.sens[[s]]$d.pos   
            Pos.pred=pred.fun(mod=Stand.out.daily.Tweedie.sens[[s]]$mod.pos,biascor="YES",PRED="finyear")
            
            Delta=full_join(Bi.pred%>%dplyr::select(finyear,prob),
                            Pos.pred%>%dplyr::select(finyear,response),by='finyear')%>%
              mutate(response=response*prob)%>%
              dplyr::select(-prob)
            
            
            return(list(Lognormal=Lognormal,Delta=Delta))
            rm(d.logN,Bi,d.pos)
          }
        }
    })
    
    names(Pred.Tweedie.sens)=names(Pred.daily.Tweedie.sens)=names(SP.list)
    
    stopCluster(cl)
    
    
    #3. Display annual indices  
    Plot.cpue.sens=function(cpuedata,cpuedata.daily,CL,CxS,Yvar,add.lines,firstyear)    
    {
      if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
      if(length(cpuedata)<=3)tc=seq(-.5*0.25,.5*0.25,length.out=length(cpuedata))
      if(length(cpuedata)==1)tc=seq(-.5*0.5,.5*0.5,length.out=length(cpuedata))
      
      #normalize
      for(l in 1:length(cpuedata)) cpuedata[[l]]$response=cpuedata[[l]]$response/mean(cpuedata[[l]]$response)
      for(l in 1:length(cpuedata.daily)) cpuedata.daily[[l]]$response=cpuedata.daily[[l]]$response/mean(cpuedata.daily[[l]]$response)
      
      
      Yrs=c(as.numeric(substr(cpuedata[[1]][,match(Yvar,names(cpuedata[[1]]))],1,4)),
            as.numeric(substr(cpuedata.daily[[1]][,match(Yvar,names(cpuedata.daily[[1]]))],1,4)))
      Tops=c(unlist(lapply(cpuedata, `[`, "response")),
             unlist(lapply(cpuedata.daily, `[`, "response")))
      Tops=Tops[!is.infinite(Tops)]
      ymax=max(Tops)
      Quant=quantile(Tops,probs=c(.95,1))
      if(diff(Quant)>3) ymax=quantile(Tops,probs=.995)
      
      plot(Yrs,Yrs,ylim=c(0,ymax),xlim=c(firstyear,max(Yrs)),ylab="",xlab="",
           col="transparent",cex.axis=1.25)
      for(l in 1:length(cpuedata))
      {
        if(!is.null(cpuedata[[l]]))aaa=cpuedata[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
        if(!is.null(cpuedata.daily[[l]]))aaa.daily=cpuedata.daily[[l]]%>%mutate(finyear=as.numeric(substr(finyear,1,4)))
        
        if(!is.null(cpuedata[[l]]))with(aaa,
                                        {
                                          if(add.lines=="NO") points(finyear+tc[l], response, pch=19, lty=2, col=CL[l],cex=CxS)
                                          if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=21, lty=2, bg=CL[l],cex=CxS)
                                          
                                        })
        if(!is.null(cpuedata.daily[[l]]))with(aaa.daily,
                                              {
                                                if(l==1)polygon(x=c(finyear[1]-.5,finyear[length(finyear)]+.5,finyear[length(finyear)]+.5,finyear[1]-.5),
                                                                y=c(0,0,ymax*.99,ymax*.99),col='grey92',border="transparent")
                                                if(add.lines=="NO") points(finyear+tc[l], response, pch=21, lty=2, col=CL[l],cex=CxS)
                                                if(add.lines=="YES") points(finyear+tc[l], response, "o", pch=21, lty=2, bg=CL[l],cex=CxS)
                                              })
      }
      
    }
    CL.sens=c("black","grey50","white")
    fn.fig("Appendix 5_Annual_Index_normalised_Tweedie_Delta_lognormal",1800, 2400)  
    par(mfrow=c(4,1),mar=c(1,3,1.5,.6),oma=c(2.5,1,.1,.3),las=1,mgp=c(1.9,.7,0))
    for(s in Tar.sp)
    {
      Plot.cpue.sens(cpuedata=list(Tweedie=Pred[[s]],
                                   Lognormal=Pred.Tweedie.sens[[s]]$Lognormal,
                                   Delta=Pred.Tweedie.sens[[s]]$Delta),
                     cpuedata.daily=list(Tweedie=Pred.daily[[s]],
                                         Lognormal=Pred.daily.Tweedie.sens[[s]]$Lognormal,
                                         Delta=Pred.daily.Tweedie.sens[[s]]$Delta),
                     CL=CL.sens,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975)
      legend("topright",Nms.sp[s],bty='n',cex=1.75)
    }
    legend("bottomleft",c("Tweedie","Lognormal","Delta-lognormal"),bty='n',pch=21,
           pt.bg=CL.sens,cex=1.45)
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative catch rate",side=2,line=-0.75,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
  }
}

#-- Sensitivity test for Delta. Get glm predictions for daily and compare to gam ----------------------------------------------
if(Use.Delta)
{
  if(Model.run=="First")
  {
    system.time({Pred.daily.tar.glm=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$Bi$res,
                                MOD=Stand.out.daily[[s]]$Pos$res,
                                BiData=Stand.out.daily[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                PosData=Stand.out.daily[[s]]$Pos$DATA,
                                niter=Niter,
                                pred.term='finyear',
                                ALL.terms=c(Predictors_daily)))
      }
    })
    names(Pred.daily.tar.glm)=names(SP.list)[Tar.sp]
    stopCluster(cl)
    
    #relative glm and gam standardised
    fn.fig("Figure.Annual_Index_glm_gam.stand",1200, 2400)    
    par(mfrow=c(4,1),mar=c(1,1,1.5,2),oma=c(2.5,3,.1,.2),las=1,mgp=c(1.9,.7,0))
    for(s in 1:length(Tar.sp))
    {
      #Daily
      GAM=Pred.daily[[Tar.sp[s]]]
      Mn=mean(GAM$response)
      GAM=GAM%>%mutate(response=response/Mn,
                       lower.CL=lower.CL/Mn,
                       upper.CL=upper.CL/Mn)
      GLM=Pred.daily.tar.glm[[s]]
      Mn=mean(GLM$response)
      GLM=GLM%>%mutate(response=response/Mn,
                       lower.CL=lower.CL/Mn,
                       upper.CL=upper.CL/Mn)
      Mon.dat=list(GAM=GAM,GLM=GLM)
      LgND="NO"
      if(s==1) LgND="YES"
      Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="finyear",add.lines="YES")
      if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      mtext(Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
    }
    mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
  }
}

#-- Plot Other species normalised ----------------------------------------------
if(Use.Delta)
{
  Plot.cpue.other=function(cpuedata,ADD.LGND,whereLGND,COL,CxS,Yvar,All.yrs)    #plot cpues other species
  {
    if(is.null(cpuedata[[1]])) plot.new()
    if(!is.null(cpuedata[[1]]))
    {
      if(inherits(cpuedata, "list")) 
      {
        if(length(cpuedata)>3)tc=seq(-1.5*0.15,1.5*0.15,length.out=length(cpuedata))
        if(length(cpuedata)<=3)tc=seq(-.5*0.15,.5*0.15,length.out=length(cpuedata))
        ymax = max(unlist(lapply(cpuedata, `[`, "upper.CL")),na.rm=T)
        Yrs=as.numeric(substr(All.yrs,1,4))
        plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent")
        if(COL=='color')CL=c("black","forestgreen", "red","bisque3","blue2","dodgerblue")
        if(COL=='grey') CL=gray.colors(length(cpuedata),start=0.2,end=0.65)
        for(l in 1:length(cpuedata))
        {
          aaa=cpuedata[[l]]%>%mutate(finyear=as.character(finyear))
          aaa=data.frame(finyear=All.yrs)%>%left_join(aaa,by="finyear")
          with(aaa,
               {
                 points(Yrs+tc[l], response, "o", pch=16, lty=2, col=CL[l],cex=CxS)
                 arrows(x0=Yrs+tc[l], y0=lower.CL, 
                        x1=Yrs+tc[l], y1=upper.CL, 
                        code=3, angle=90, length=0.05, col=CL[l])
               })
          if(ADD.LGND=="YES") legend(whereLGND,names(cpuedata),bty='n',pch=16,col=CL,cex=1.45)
        }
      }
      if(inherits(cpuedata, "data.frame"))
      {
        ymax = max(cpuedata$upper.CL,na.rm=T)
        Yrs=as.numeric(substr(cpuedata[,match(Yvar,names(cpuedata))],1,4))
        plot(Yrs,Yrs,ylim=c(0,ymax),ylab="",xlab="",col="transparent")
        CL="black"
        points(Yrs, cpuedata$response, "o", pch=16, lty=2, col=CL,cex=CxS)
        arrows(x0=Yrs, y0=cpuedata$lower.CL, 
               x1=Yrs, y1=cpuedata$upper.CL, 
               code=3, angle=90, length=0.05, col=CL)
      }
    }
  }
  fn.fig("Figure 3.Annual_Index_normalised_other species",1200, 2400)    
  par(mfrow=c(length(nnn[-sort(Tar.sp)]),2),mar=c(1,1,.75,.95),oma=c(2.5,3,1,.25),las=1,mgp=c(1.9,.5,0),cex.axis=.8)
  for(s in nnn[-sort(Tar.sp)])
  {
    #Monthly
    Mon.dat=list(Standardised=Pred.normlzd[[s]])
    LgND="NO"
    suppressWarnings(Plot.cpue.other(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1,Yvar="finyear",All.yrs=FINYEAR.monthly))
    if(s==nnn[-sort(Tar.sp)][1]) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.25)
    
    #Daily
    Daily.dat=list(Standardised=Pred.daily.normlzd[[s]])
    LgND="NO"
    suppressWarnings(Plot.cpue.other(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1,Yvar="finyear",All.yrs=FINYEAR.daily))
    if(s==nnn[-sort(Tar.sp)][1]) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.25)
    mtext(gsub(paste("Shark", collapse="|"), "", Nms.sp[s]),4,0.1,cex=.7,las=3)
  }
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.25,outer=T)
  mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.25,outer=T)
  dev.off()
}

if(Use.Tweedie)   
{
  fn.fig("Figure 3.Annual_Index_normalised_other species",1200, 2400)    
  par(mfrow=c(length(nnn[-sort(Tar.sp)]),1),mar=c(1,1,.75,.95),oma=c(2.5,3,1,.25),las=1,mgp=c(1.9,.5,0),cex.axis=.8)
  for(s in nnn[-sort(Tar.sp)])
  {
    Plot.cpue(cpuedata=list(Standardised=Pred.normlzd[[s]]),
              cpuedata.daily=list(Standardised=Pred.daily.normlzd[[s]]),
              CL=CLs,CxS=CxS,Yvar="finyear",add.lines="YES",firstyear=1975,
              ADD.nomnl=list(rbind(Unstandardised.normlzd[[s]],Unstandardised.daily.normlzd[[s]])),
              Add.CI.nominal=FALSE)
    legend("topleft",Nms.sp[s],bty='n',cex=1.25)
  }
  legend("bottomleft",c("Standardised","Nominal"),bty='n',pch=21,
         pt.bg=c(CLs[1],"white"),cex=1)
  mtext("Financial year",side=1,line=1.2,font=1,las=0,cex=1.25,outer=T)
  mtext("Relative CPUE",side=2,line=1.15,font=1,las=0,cex=1.25,outer=T)
  dev.off()
}


#-- Boxplots of factors ----------------------------------------------
if(Use.Delta)
{
  if(Model.run=="First")
  {
    for(s in Tar.sp)
    {
      fn.fig(paste("boxplot.daily.preds_",Nms.sp[s],sep=""),2000, 2400)  
      par(mfcol=c(3,2))  
      with(Stand.out.daily[[s]]$DATA,
           {
             Ylim=c(0,quantile(catch.target/km.gillnet.hours.c,probs=0.99))
             boxplot((catch.target/km.gillnet.hours.c)~mean.depth,col="grey70",ylim=Ylim)
             boxplot((catch.target/km.gillnet.hours.c)~lunar,col="grey70",ylim=Ylim)
             boxplot((catch.target/km.gillnet.hours.c)~month,col="grey70",ylim=Ylim)
             boxplot((catch.target/km.gillnet.hours.c)~mesh,col="grey70",ylim=Ylim)
             boxplot((catch.target/km.gillnet.hours.c)~shots.c,col="grey70",ylim=Ylim)
           })
      dev.off()
    }
  }
}


# -- Show effects for other terms  ----------------------------------------------
if(Model.run=="First")
{
  #Vessel effect
  if(Use.Qualif.level)
  {
    #emmeans (formerly lsmeans) considering log bias corr if required
    #monthly          takes 80 sec
    system.time({Pred.vess=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
      {
        d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
        return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="vessel",Pred.type="link"))
        rm(d)
      }
    })
    #Daily              takes 70 sec
    system.time({Pred.vess.daily=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
      {
        d=Stand.out.daily[[s]]$DATA
        return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="vessel",Pred.type="link"))
        rm(d)
      }
    })
    
  }
  if(Use.Delta)  
  {
    #monthly          takes 80 sec
    system.time({Pred.vess=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$Bi$res,
                                MOD=Stand.out[[s]]$Pos$res,
                                BiData=Stand.out[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                PosData=Stand.out[[s]]$Pos$DATA,
                                niter=Niter,
                                pred.term='vessel',
                                ALL.terms=Predictors_monthly))
      }
    })
    #Daily              takes 70 sec
    system.time({Pred.vess.daily=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$Bi$res.gam,
                                MOD=Stand.out.daily[[s]]$Pos$res.gam,
                                BiData=Stand.out.daily[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                PosData=Stand.out.daily[[s]]$Pos$DATA,
                                niter=Niter,
                                pred.term='vessel',
                                ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
      }
    })
    
    names(Pred.vess)=names(Pred.vess.daily)=names(SP.list)[Tar.sp]
    fn.fig("Figure 2.Vessel effect",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,2.5,.1,.2),las=1,mgp=c(1.9,.7,0))
    for(s in 1:length(Pred.vess))
    {
      #Monthly
      Mon.dat=Pred.vess[[s]]
      LgND="NO"
      Mon.dat$vessel=1:nrow(Mon.dat)
      Mn=mean(Mon.dat$response)
      Mon.dat=Mon.dat%>%mutate(response=response/Mn,
                               lower.CL=lower.CL/Mn,
                               upper.CL=upper.CL/Mn)
      if(s==1)LgND="YES"
      Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="vessel",add.lines="NO")
      if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      #Daily
      Daily.dat=Pred.vess.daily[[s]]
      LgND="NO"
      Daily.dat$vessel=1:nrow(Daily.dat)
      Mn=mean(Daily.dat$response)
      Daily.dat=Daily.dat%>%mutate(response=response/Mn,
                                   lower.CL=lower.CL/Mn,
                                   upper.CL=upper.CL/Mn)
      Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="vessel",add.lines="NO")
      if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      mtext( Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
    }
    mtext("Vessel",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=0.5,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
  }
  if(Use.Tweedie)  
  {
    Store.vessel=vector('list',length(Tar.sp))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      d=Stand.out.daily[[s]]$DATA
      if(length(grep("vessel",attr(terms(Best.Model.daily[[s]]), "term.labels")))>0)
      {
        a=pred.fun(Stand.out.daily[[s]]$res.gam,biascor="NO",PRED='vessel')%>%
          arrange(response)
        Mn=mean(a$response)
        Store.vessel[[i]]=a%>%
          mutate(response=response/Mn,
                 lower.CL=lower.CL/Mn,
                 upper.CL=upper.CL/Mn)
      }
      
    }
  }
  
  #Month effect    
  if(Use.Qualif.level)
  {
    #emmeans (formerly lsmeans) considering log bias corr if required
    #monthly          takes 80 sec
    system.time({Pred.month=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
      {
        d=Stand.out[[s]]$DATA   #note: need data as global for ref_grid
        return(pred.fun(MOD=Stand.out[[s]]$res,biascor="YES",PRED="month",Pred.type="link"))
        rm(d)
      }
    })
    #Daily              takes 70 sec
    system.time({Pred.month.daily=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %do%
      {
        d=Stand.out.daily[[s]]$DATA
        return(pred.fun(MOD=Stand.out.daily[[s]]$res.gam,biascor="YES",PRED="month",Pred.type="link"))
        rm(d)
      }
    })
  }
  if(Use.Delta)  
  {
    #monthly          takes 80 sec
    system.time({Pred.month=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out[[s]]$Bi$res,
                                MOD=Stand.out[[s]]$Pos$res,
                                BiData=Stand.out[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                PosData=Stand.out[[s]]$Pos$DATA,
                                niter=Niter,
                                pred.term='month',
                                ALL.terms=Predictors_monthly))
      }
    })
    #Daily              takes 70 sec
    system.time({Pred.month.daily=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue(BiMOD=Stand.out.daily[[s]]$Bi$res.gam,
                                MOD=Stand.out.daily[[s]]$Pos$res.gam,
                                BiData=Stand.out.daily[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                PosData=Stand.out.daily[[s]]$Pos$DATA,
                                niter=Niter,
                                pred.term='month',
                                ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner')))
      }
    })
    
    names(Pred.month)=names(Pred.month.daily)=names(SP.list)[Tar.sp]
    stopCluster(cl)
    
    fn.fig("Figure.Monthly effects",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,2.5,.1,.2),las=1,mgp=c(1.9,.7,0))
    for(s in 1:length(Pred.vess))
    {
      #Monthly
      Mon.dat=Pred.month[[s]]
      LgND="NO"
      if(s==1)LgND="YES"
      Mn=mean(Mon.dat$response)
      Mon.dat=Mon.dat%>%mutate(response=response/Mn,
                               lower.CL=lower.CL/Mn,
                               upper.CL=upper.CL/Mn)
      Plot.cpue(cpuedata=Mon.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month",add.lines="YES")
      if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.5)
      
      #Daily
      Daily.dat=Pred.month.daily[[s]]
      LgND="NO"
      Mn=mean(Daily.dat$response)
      Daily.dat=Daily.dat%>%mutate(response=response/Mn,
                                   lower.CL=lower.CL/Mn,
                                   upper.CL=upper.CL/Mn)
      Plot.cpue(cpuedata=Daily.dat,ADD.LGND=LgND,whereLGND='topright',COL="grey",CxS=1.35,Yvar="month",add.lines="YES")
      if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.5)
      mtext( Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
    }
    mtext("Month",side=1,line=1.2,font=1,las=0,cex=1.5,outer=T)
    mtext("Relative CPUE",side=2,line=0,font=1,las=0,cex=1.5,outer=T)
    dev.off()
    
  }
  if(Use.Tweedie)  
  {
    Store.month=vector('list',length(Tar.sp))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      Store.month[[i]]=pred.fun.continuous(d=Stand.out.daily[[s]]$DATA,
                                           mod=Stand.out.daily[[s]]$res.gam,
                                           PRED='month',
                                           Formula=Best.Model.daily[[s]])
    }
  }
  
  #Depth effect
  if(Use.Tweedie)  
  {
    Store.depth=vector('list',length(Tar.sp))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      Store.depth[[i]]=pred.fun.continuous(d=Stand.out.daily[[s]]$DATA,
                                           mod=Stand.out.daily[[s]]$res.gam,
                                           PRED='mean.depth',
                                           Formula=Best.Model.daily[[s]])
    }
  }
  
  #Spatial cpue (monthly blocks and daily lat/long)  
  if(Use.Qualif.level|Use.Delta)
  {
    pred.fun.spatial=function(DAT,MOD,PRED,FORM,Spatial.grid)
    {
      TermS=all.vars(FORM)
      TermS=TermS[-match(c("LNcpue",PRED),TermS)]
      id.cat=TermS[which(TermS%in%Categorical)]
      NewDat=as.data.frame(matrix(nrow=1,ncol=length(id.cat)))
      names(NewDat)=id.cat
      for(ii in 1:length(id.cat))
      {
        dummy=sort(table(DAT[,match(id.cat[ii],names(DAT))]))
        NewDat[,ii]=factor(names(dummy[length(dummy)]),levels(DAT[,match(id.cat[ii],names(DAT))]))
      }
      id.cont=TermS[which(!TermS%in%Categorical)]
      if(length(id.cont)>0)
      {
        NewDat.cont=as.data.frame(matrix(nrow=1,ncol=length(id.cont)))
        names(NewDat.cont)=id.cont
        for(ii in 1:length(id.cont)) NewDat.cont[,ii]=mean(DAT[,match(id.cont[ii],names(DAT))])
        NewDat=cbind(NewDat,NewDat.cont)
      }
      NewDat=cbind(Spatial.grid,NewDat)
      PRD=predict(MOD,newdata=NewDat,type='link',se.fit=T)
      NewDat$cpue=exp(PRD$fit+(PRD$se.fit^2)/2)
      return(NewDat)
    }
  }
  if(Use.Qualif.level)
  {
    #monthly          takes 4 sec
    system.time({Pred.spatial.monthly=foreach(s=Tar.sp,.packages=c('dplyr','emmeans')) %dopar%
      {
        
        return(pred.fun.spatial(DAT=Stand.out[[s]]$DATA,MOD=Stand.out[[s]]$res,
                                PRED='blockx',FORM=Best.Model[[s]],
                                Spatial.grid=data.frame(blockx=factor(levels(Stand.out[[s]]$DATA$blockx))))
        )
        
      }
    })
    #Daily              takes 0.1 sec
    system.time({Pred.spatial.daily=foreach(s=Tar.sp,.packages=c('dplyr')) %do%
      {
        return(pred.fun.spatial(DAT=Stand.out.daily[[s]]$DATA,MOD=Stand.out.daily[[s]]$res.gam,
                                PRED=c('long10.corner','lat10.corner'),FORM=Best.Model.daily.gam[[s]],
                                Spatial.grid=Stand.out.daily[[s]]$DATA%>%select(block10,lat10.corner,long10.corner)%>%
                                  distinct(block10,.keep_all=T)))
      }
    })
  }
  if(Use.Delta)  
  {
    Plot.cpue.spatial=function(cpuedata,var,Pol.x,Pol.y,add.pol,delta,show)
    {
      if(var[1]=='blockx')cpuedata=cpuedata%>%
          mutate( Lat=-round(as.numeric(substr(get(var),1,2)),2),
                  Long=round(100+as.numeric(substr(get(var),3,4)),2))else
                    cpuedata=cpuedata%>%mutate( Lat=round(get(var[2]),2),
                                                Long=round(get(var[1]),2))
                  
                  YLIM=floor(range(Full.lat))    
                  XLIM=floor(range(Full.long)) 
                  
                  misn.lat=sort(Full.lat[which(!Full.lat%in%unique(cpuedata$Lat))])
                  misn.lon=sort(Full.long[which(!Full.long%in%unique(cpuedata$Long))])
                  if(length(misn.lat)>0 | length(misn.lon)>0)
                  {
                    if(var[1]=='blockx')
                    {
                      combo=expand.grid(Lon=Full.long,Lat=Full.lat)%>%
                        mutate(blockx=paste(abs(Lat),Lon-100,sep=''))%>%
                        select(blockx)
                      cpuedata=cpuedata%>%mutate(blockx=as.character(blockx))
                      cpuedata=combo%>%left_join(cpuedata,by=var)%>%
                        mutate( Lat=-round(as.numeric(substr(get(var),1,2)),2),
                                Long=round(100+as.numeric(substr(get(var),3,4)),2))
                    }else
                    {
                      combo=expand.grid(Long=Full.long,Lat=Full.lat)
                      cpuedata=combo%>%left_join(cpuedata,by=c('Long','Lat'))
                    }
                  }
                  cpuedata=cpuedata%>%select(c(cpue,Lat,Long)) 
                  cpuedata.spread=cpuedata%>%spread(Lat,cpue)
                  Lon=as.numeric(cpuedata.spread$Long)
                  cpuedata.spread=as.matrix(cpuedata.spread[,-1]) 
                  LaT=as.numeric(colnames(cpuedata.spread))
                  brk<- quantile( c(cpuedata.spread),probs=seq(0,1,.1),na.rm=T)
                  YLIM[1]=YLIM[1]-0.5
                  YLIM[2]=YLIM[2]+0.5
                  XLIM[1]=XLIM[1]-0.5
                  XLIM[2]=XLIM[2]+0.5
                  plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
                  if(add.pol) polygon(Pol.x,Pol.y,col=rgb(.1,.2,.1,alpha=.2),border="transparent")
                  image(Lon+delta,LaT-delta,cpuedata.spread, breaks=brk, col=rev(heat.colors(length(brk)-1)),add=T)
                  par(new=T)
                  plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
                  if(show)suppressWarnings(image.plot(Lon+delta,LaT-delta,cpuedata.spread, breaks=brk, 
                                                      col=rev(heat.colors(length(brk)-1)),lab.breaks=names(brk),add=T,
                                                      legend.only=T,legend.shrink=.5,legend.mar=c(25,5)))
                  axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = 0.15)
                  axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =0.15)
                  n=seq(South.WA.long[1],South.WA.long[2],2)
                  axis(side = 1, at =n, labels = n, tcl = 0.3)
                  n=seq(South.WA.lat[1],South.WA.lat[2],2)
                  axis(side = 2, at = n, labels = abs(n),las=2,tcl =0.3)
    }
    
    #monthly          takes 4 sec
    system.time({Pred.spatial.monthly=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %dopar%
      {
        
        return(fn.MC.delta.cpue_spatial(BiMOD=Stand.out[[s]]$Bi$res,
                                        MOD=Stand.out[[s]]$Pos$res,
                                        BiData=Stand.out[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                        PosData=Stand.out[[s]]$Pos$DATA,
                                        niter=Niter,
                                        pred.term='blockx',
                                        ALL.terms=Predictors_monthly,
                                        Spatial.grid=data.frame(blockx=factor(levels(Stand.out[[s]]$Pos$DATA$blockx)))))
      }
    })
    #Daily              takes 0.1 sec
    system.time({Pred.spatial.daily=foreach(s=Tar.sp,.packages=c('dplyr','mvtnorm')) %do%
      {
        return(fn.MC.delta.cpue_spatial(BiMOD=Stand.out.daily[[s]]$Bi$res.gam,
                                        MOD=Stand.out.daily[[s]]$Pos$res.gam,
                                        BiData=Stand.out.daily[[s]]$Bi$DATA%>%mutate(LNeffort=LN.effort),
                                        PosData=Stand.out.daily[[s]]$Pos$DATA,
                                        niter=Niter,
                                        pred.term=c('long10.corner','lat10.corner'),
                                        ALL.terms=c(Predictors_daily,'long10.corner','lat10.corner'),
                                        Spatial.grid=Stand.out.daily[[s]]$Pos$DATA%>%select(block10,lat10.corner,long10.corner)%>%
                                          distinct(block10,.keep_all=T)%>%select(-block10)))
      }
    })
    
    
    names(Pred.spatial.monthly)=names(Pred.spatial.daily)=names(SP.list)[Tar.sp]
    stopCluster(cl)
    
    South.WA.lat=c(-36,-25); South.WA.long=c(112,129)
    data(worldLLhigh)
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    PLATE=c(.01,.9,.075,.9)
    Dusky=c(X1=South.WA.long[1],X2=Dusky.range[2],Y1=South.WA.lat[1],Y2=Dusky.range[1])
    Sandbar=c(X1=South.WA.long[1],X2=Sandbar.range[2],Y1=South.WA.lat[1],Y2=Sandbar.range[1])
    Whiskery=c(X1=South.WA.long[1],X2=Whiskery.range[2],Y1=South.WA.lat[1],Y2=Whiskery.range[1])
    Gummy=c(X1=Gummy.range[1],X2=Gummy.range[2],Y1=South.WA.lat[1],Y2=-31.6)
    LISta=list(Whiskery=Whiskery,Gummy=Gummy,Dusky=Dusky,Sandbar=Sandbar)
    
    fn.fig("Figure 3.Spatial effect",2000, 2400)    
    par(mfrow=c(4,2),mar=c(1,2,1.5,2),oma=c(2.5,1,.1,1),las=1,mgp=c(1.9,.5,0))
    for(s in 1:length(Pred.vess))
    {
      sHow=F
      #Monthly
      Full.long=seq(113,129)
      Full.lat=seq(-36,-26) 
      if(Use.Delta) Pred.spatial.monthly[[s]]$cpue=Pred.spatial.monthly[[s]]$response
      Plot.cpue.spatial(cpuedata=Pred.spatial.monthly[[s]],var='blockx',
                        Pol.x=c(LISta[[s]][1],LISta[[s]][2],LISta[[s]][2],LISta[[s]][1]),
                        Pol.y=c(LISta[[s]][4],LISta[[s]][4],LISta[[s]][3],LISta[[s]][3]),
                        add.pol=T,delta=.5,show=sHow)
      if(s==1) mtext("Monthly returns",side=3,line=0,font=1,las=0,cex=1.35)
      
      
      #Daily
      if(s==1)sHow=T
      if(Use.Delta) Pred.spatial.daily[[s]]$cpue=Pred.spatial.daily[[s]]$response
      Full.long=apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(113:129)),rep(113:129,each=6)),1,sum)
      Full.lat=-apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(35:26)),rep(35:26,each=6)),1,sum)
      Plot.cpue.spatial(cpuedata=Pred.spatial.daily[[s]],var=c('long10.corner','lat10.corner'),
                        Pol.x=c(LISta[[s]][1],LISta[[s]][2],LISta[[s]][2],LISta[[s]][1]),
                        Pol.y=c(LISta[[s]][4],LISta[[s]][4],LISta[[s]][3],LISta[[s]][3]),
                        add.pol=T,delta=.16,show=sHow)
      
      
      if(s==1) mtext("Daily logbooks",side=3,line=0,font=1,las=0,cex=1.35)
      mtext( Nms.sp[Tar.sp[s]],4,line=1,las=3,cex=1.5)
    }
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.4,font=1,las=0,cex=1.35,outer=T)
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-0.85,font=1,las=0,cex=1.35,outer=T)
    dev.off()
  }
  if(Use.Tweedie)  
  {
    Plot.cpue.spatial=function(cpuedata,var,Pol.x,Pol.y,add.pol,delta,show)
    {
      cpuedata=cpuedata%>%mutate( Lat=round(get(var[2]),2),
                                  Long=round(get(var[1]),2))
      YLIM=floor(range(Full.lat))    
      XLIM=floor(range(Full.long)) 
      misn.lat=sort(Full.lat[which(!Full.lat%in%unique(cpuedata$Lat))])
      misn.lon=sort(Full.long[which(!Full.long%in%unique(cpuedata$Long))])
      if(length(misn.lat)>0 | length(misn.lon)>0)
      {
        combo=expand.grid(Long=Full.long,Lat=Full.lat)
        cpuedata=combo%>%left_join(cpuedata,by=c('Long','Lat'))
      }
      cpuedata=cpuedata%>%
        mutate(cpue=Pred)%>%
        dplyr::select(c(cpue,Lat,Long)) 
      cpuedata.spread=cpuedata%>%
        spread(Lat,cpue)
      Lon=as.numeric(cpuedata.spread$Long)
      cpuedata.spread=as.matrix(cpuedata.spread[,-1]) 
      LaT=as.numeric(colnames(cpuedata.spread))
      brk<- quantile( c(cpuedata.spread),probs=seq(0,1,.1),na.rm=T)
      YLIM[1]=YLIM[1]-0.5
      YLIM[2]=YLIM[2]+0.5
      XLIM[1]=XLIM[1]-0.5
      XLIM[2]=XLIM[2]+0.5
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      if(add.pol) polygon(Pol.x,Pol.y,col=rgb(.1,.2,.1,alpha=.2),border="transparent")
      image(Lon+delta,LaT-delta,cpuedata.spread, breaks=brk, col=rev(heat.colors(length(brk)-1)),add=T)
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      
      if(show)suppressWarnings(image.plot(Lon+delta,LaT-delta,cpuedata.spread, breaks=brk, 
                                          col=rev(heat.colors(length(brk)-1)),lab.breaks=names(brk),add=T,
                                          legend.only=T,legend.shrink=.5,legend.mar=c(25,5)))
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = -0.15)
      axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =-0.15)
      n=seq(South.WA.long[1],South.WA.long[2],2)
      axis(side = 1, at =n, labels = n, tcl = -0.3)
      n=seq(South.WA.lat[1],South.WA.lat[2],2)
      axis(side = 2, at = n, labels = abs(n),las=2,tcl =-0.3)
    }
    South.WA.lat=c(-36,-25)
    #South.WA.lat=c(-35,-25)
    South.WA.long=c(112,129)
    data(worldLLhigh)
    PLATE=c(.01,.9,.075,.9)
    LISta=list(Whiskery=Core$`Whiskery Shark`,
               Gummy=Core$`Gummy Shark`,
               Dusky=Core$`Dusky Whaler`,
               Sandbar=Core$`Sandbar Shark`)
    for(l in 1:length(LISta)) LISta[[l]]$Lat[1]=-35.25
    Full.long=apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(113:129)),rep(113:129,each=6)),1,sum)
    Full.lat=-apply(cbind(rep(c(.83,.67,.5,.33,.17,0),length(abs(South.WA.lat[1]):26)),
                          rep(abs(South.WA.lat[1]):26,each=6)),1,sum)
    
    Store.spatial=vector('list',length(Tar.sp))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      Store.spatial[[i]]=pred.fun.continuous(d=Stand.out.daily[[s]]$DATA,
                                             mod=Stand.out.daily[[s]]$res.gam,
                                             PRED=c('long10.corner', 'lat10.corner'),
                                             Formula=Best.Model.daily[[s]])
    }
  }
  
  #Cluster
  if(Use.Tweedie)  
  {
    Store.cluster=vector('list',length(Tar.sp))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      d=Stand.out.daily[[s]]$DATA
      if(length(grep("cluster_clara",attr(terms(Best.Model.daily[[s]]), "term.labels")))>0)
      {
        a=pred.fun(Stand.out.daily[[s]]$res.gam,biascor="NO",PRED='cluster_clara')%>%
          arrange(response)
        Mn=mean(a$response)
        Store.cluster[[i]]=a%>%
          mutate(response=response/Mn,
                 lower.CL=lower.CL/Mn,
                 upper.CL=upper.CL/Mn)
      }
      
    }
  } 
  
  #Plot effects  
  if(Use.Tweedie)
  {
    CX.t=1.25
    show.vessel=TRUE
    
    fn.fig("Figure 2.Other terms effect",2400, 2400)   
    if(show.vessel) 
    {
      par(mfrow=c(5,4),mar=c(2.5,2,1.2,1),oma=c(1,2.5,.5,.2),las=1,mgp=c(1.9,.65,0),cex.axis=1.25)
    }else
    {
      par(mfrow=c(4,4),mar=c(2.5,2,1.2,1),oma=c(1,2.5,.5,.2),las=1,mgp=c(1.9,.65,0),cex.axis=1.25)
    }
    
    #Vessel effect
    if(show.vessel)
    {
      for(i in 1:length(Tar.sp))
      {
        s=Tar.sp[i]
        a=Store.vessel[[i]]
        plot(1:nrow(a),a$response,pch=19,cex=1.5,ylim=c(0,max(a$upper.CL)),
             #plot(1:nrow(a),a$response,pch=19,cex=1.5,ylim=c(min(a$lower.CL),max(a$upper.CL)),
             ylab='',xlab='')
        arrows(1:nrow(a), a$lower.CL,1:nrow(a), a$upper.CL,code=3, angle=90, length=0.05)
        mtext(Nms.sp[s],3,0.15,cex=1.25)
        if(i==2) mtext("                             Vessel",
                       side=1,line=2,font=1,las=0,cex=CX.t)
      }
    }
    
    #Cluster effect
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      a=Store.cluster[[i]]
      plot(1:nrow(a),a$response,pch=19,cex=1.5,ylim=c(0,max(a$upper.CL)),
           #plot(1:nrow(a),a$response,pch=19,cex=1.5,ylim=c(min(a$lower.CL),max(a$upper.CL)),
           ylab='',xlab='',xaxt='n',xlim=c(0.5,2.5))
      arrows(1:nrow(a), a$lower.CL,1:nrow(a), a$upper.CL,code=3, angle=90, length=0.05)
      axis(1,1:2,c("Group 1","Group 2"))
      if(i==2) mtext("                             Fishing tactic",
                     side=1,line=2,font=1,las=0,cex=CX.t)
    }
    
    #Month effect    
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      PRED=Store.month[[i]]
      MN=mean(PRED$Pred)
      PRED=PRED%>%
        mutate(Mean=Pred/MN,
               Upper=(Pred+1.96*Pred.SE)/MN,
               Lower=(Pred-1.96*Pred.SE)/MN)
      with(PRED,
           {
             plot(month,Mean,type='l',lwd=2,ylim=c(0,max(Upper)),xlab='')
             #plot(month,Mean,type='l',lwd=2,ylim=c(min(Lower),max(Upper)),xlab='')
             polygon(c(month,rev(month)),
                     c(Upper,rev(Lower)),
                     col=rgb(.1,.1,.1,alpha=.2),border=rgb(.1,.1,.1,alpha=.4))
           })
      if(!show.vessel) mtext(Nms.sp[s],3,0.15,cex=1.25)
      if(i==2) mtext("                             Month",
                     side=1,line=2,font=1,las=0,cex=CX.t)
      if(i==1) mtext("Relative catch rate",side=2,line=2.5,font=1,las=0,cex=CX.t)
    }
    
    #Depth effect    
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      PRED=Store.depth[[i]]
      if(!is.null(PRED))
      {
        MN=mean(PRED$Pred)
        PRED=PRED%>%
          mutate(Mean=Pred/MN,
                 Upper=(Pred+1.96*Pred.SE)/MN,
                 Lower=(Pred-1.96*Pred.SE)/MN)
        with(PRED,
             {
               plot(mean.depth,Mean,type='l',lwd=2,ylim=c(0,max(Upper)),xlab='')
               #plot(mean.depth,Mean,type='l',lwd=2,ylim=c(min(Lower),max(Upper)),xlab='')
               polygon(c(mean.depth,rev(mean.depth)),
                       c(Upper,rev(Lower)),
                       col=rgb(.1,.1,.1,alpha=.2),border=rgb(.1,.1,.1,alpha=.4))
             })
      }else
      {
        plot.new()
      }
      if(i==2) mtext("                             Depth (m)",
                     side=1,line=2,font=1,las=0,cex=CX.t)
    }
    
    #Lat long effect
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    for(i in 1:length(Tar.sp))
    {
      s=Tar.sp[i]
      PRED=Store.spatial[[i]]
      sHow.legend=FALSE
      if(i==1)sHow.legend=TRUE
      Plot.cpue.spatial(cpuedata=PRED,
                        var=c('long10.corner','lat10.corner'),
                        Pol.x=c(LISta[[i]]$Long[1],LISta[[i]]$Long[2],LISta[[i]]$Long[2],LISta[[i]]$Long[1]),
                        Pol.y=c(LISta[[i]]$Lat[2],LISta[[i]]$Lat[2],LISta[[i]]$Lat[1],LISta[[i]]$Lat[1]),
                        add.pol=T,
                        delta=.16,
                        show=sHow.legend)
      if(i==2) mtext(expression(paste("                             Longitude (",degree,"E)",sep="")),
                     side=1,line=2.4,font=1,las=0,cex=CX.t)
      if(i==1) mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=2.5,font=1,las=0,cex=CX.t)
    }
    
    dev.off() 
  }
  
  
}
stopCluster(cl) 
