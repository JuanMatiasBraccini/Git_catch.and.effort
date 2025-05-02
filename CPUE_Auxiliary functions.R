
# Exploratory analysis -------------------------------------------------------
fn.choose.sp=function(d,crit,Titl)
{
  d=d%>%filter(!SPECIES==22999)
  Tab=table(d$SPECIES)
  Tab=Tab[Tab>50]
  this.sp=names(Tab)
  d=d%>%filter(SPECIES%in%this.sp)
  
  
  s=d%>%
    mutate(N=1,
           YR=as.numeric(substr(FINYEAR,1,4)),
           key=!!as.name(crit),
           dummy=paste(key,SPECIES,FINYEAR))%>%
    distinct(dummy,.keep_all = T)
  s1=s%>%
    group_by(SPECIES,YR)%>%
    summarise(n=sum(N))
  
  s2=s%>%distinct(key,.keep_all = T)%>%group_by(YR)%>%tally()%>%rename(Shots=n)
  
  s1=left_join(s1,s2,by='YR')%>%
    mutate(Presence=100*n/Shots,
           Absence=100*(Shots-n)/Shots)%>%
    dplyr::select(-c(n,Shots))
  
  s3=s1%>%gather("pos","n",-SPECIES,-YR)%>%arrange(SPECIES,YR)
  p=s3%>%
    ggplot(aes(x = YR,y = n)) + 
    geom_bar(aes(fill = pos), position = "stack", stat="identity")+
    facet_wrap(~SPECIES)+ ylab("Percentage") + labs(fill = "")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    ggtitle(Titl)
  print(p)
}
cfac=function(x,breaks=NULL)    #function for converting continuous var to factor
{
  if(is.null(breaks)) breaks=unique(quantile(x,probs = seq(0, 1, 0.1),na.rm = T))
  x=cut(x,breaks,include.lowest=T,right=F)
  levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                                                 c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
  return(x)
}
clog=function(x) log(x+0.05)        #function for applying log

makecategorical=function (labelModel, indat) 
{
  Interact <- grep(":", labelModel)
  nInteract <- length(Interact)
  numvars <- length(labelModel) - nInteract
  for (fac in 1:numvars) {
    if (length(indat[, labelModel[fac]]) > 0) {
      indat[, labelModel[fac]] <- factor(indat[, labelModel[fac]])
    }
    else {
      warning(paste0("Factor name ", labelModel[fac], 
                     "does not appear in data.frame"))
    }
  }
  return(indat)
}
fn.box.plt.year=function(d)
{
  d$cpue=d$Catch.Target/d$Km.Gillnet.Days.c
  par(mfcol=c(3,1),mar=c(4,4,1,1),mgp=c(2,.7,0))
  boxplot(Catch.Target~FINYEAR,d,ylab="Catch (kg)")
  boxplot(Km.Gillnet.Days.c~FINYEAR,d,ylab="Km gn d")
  boxplot(cpue~FINYEAR,d,ylab="cpue")
  
  FINYRS=sort(unique(d$FINYEAR))
  for(xx in 1:length(FINYRS))    boxplot(cpue~BLOCKX,subset(d,FINYEAR==FINYRS[xx]),ylab="cpue",main=FINYRS[xx])
  
  Yrs=as.numeric(substr(sort(unique(d$FINYEAR)),1,4))
  
  Ag.blk=aggregate(cpue~FINYEAR+BLOCKX,d,mean)
  Ag.blk= reshape(Ag.blk, v.names = "cpue", idvar="BLOCKX",timevar ="FINYEAR", direction = "wide")
  blks=Ag.blk$BLOCKX
  Ag.blk=Ag.blk[,-match('BLOCKX',names(Ag.blk))]
  CL=rainbow(ncol(Ag.blk))
  plot(Yrs,Ag.blk[1,],col=CL[1],pch=19,ylim=c(0,max(Ag.blk,na.rm=T)))
  for(qq in 2:ncol(Ag.blk)) points(Yrs,Ag.blk[qq,],pch=19,col=CL[qq])
  legend("topright",paste(blks),bty='n',pch=19,col=CL)
  
  Ag.vsl=aggregate(cpue~FINYEAR+VESSEL,d,mean)
  Ag.vsl= reshape(Ag.vsl, v.names = "cpue", idvar="VESSEL",timevar ="FINYEAR", direction = "wide")
  vsls=Ag.vsl$VESSEL
  Ag.vsl=Ag.vsl[,-match('VESSEL',names(Ag.vsl))]
  CL=rainbow(ncol(Ag.vsl))
  plot(Yrs,Ag.vsl[1,],col=CL[1],pch=19,ylim=c(0,max(Ag.vsl,na.rm=T)))
  for(qq in 2:ncol(Ag.vsl)) points(Yrs,Ag.vsl[qq,],pch=19,col=CL[qq])
  legend("topright",paste(vsls),bty='n',pch=19,col=CL)
}
fn.plt=function(a,y,TITL)
{
  plot(1:nrow(a),ylim=c(0,max(a,na.rm=T)),col="transparent",ann=F,axes=F)
  CL=rainbow(ncol(a))
  for(pp in 1:ncol(a)) lines(1:nrow(a),a[,pp],col=CL[pp],lwd=4)
  axis(1,1:nrow(a),rownames(a))
  nn=seq(0,max(a,na.rm=T),length.out=5)
  axis(2,nn,round(nn))
  mtext(y,2,3,las=3,cex=1.5)
  legend("topright",colnames(a),text.col=CL,bty='n',title=TITL)
}
fn.expl.cede=function(d,PREDS,kg,Do.ggplts)    #function for exploratory analysis
{
  PREDS=PREDS[which(PREDS%in%colnames(d))]
  
  Yrs=length(unique(d$year.c))
  div=1
  if(kg) div=1000  #in tonnes
  
  fn.plt(tapsum(d,"catch.target","finyear","month",div=div),"Catch","month")
  fn.plt(tapsum(d,"catch.target","finyear","zone",div=div),"Catch","zone")
  fn.plt(tapsum(d,"km.gillnet.hours.c","finyear","zone",div=1.0),"Effort","km.gillnet.hours")
  
  #explore turn over of vessel per year
  cbv <- tapsum(d,"catch.target","vessel","year.c",div=div) # often more vessels than years
  total <- rowSums(cbv,na.rm=TRUE)
  cbv1 <- cbv[order(total),] 
  to <- turnover(cbv1)    
  yearBubble(cbv1,ylabel="sqrt(catch-per-vessel)",diam=0.125,txt=c(2,3,4,5),hline=TRUE)
  
  plot.new()
  grid.table(to)
  
  #depth bin selection
  if(!is.na(match("mean.depth",PREDS)))
  {
    par(mfcol=c(2,2),mar=c(2,2,2,.1))
    barplot(table(trunc(d$mean.depth/2) * 2),main="2 m bin")
    barplot(table(trunc(d$mean.depth/5) * 5),main="5 m bin")
    barplot(table(trunc(d$mean.depth/10) * 10),main="10 m bin")
    barplot(table(trunc(d$mean.depth/25) * 25),main="25 m bin")
    mtext("Depth categories",3,-2,outer=T,col=2)
    
    cc <- histyear(d,Lbound=0,Rbound=max(d$mean.depth)*1.1,inc=10,pickvar="mean.depth",
                   years="year.c",varlabel="Depth (m)",plots=n2mfrow(Yrs),vline=120)
    
    d$DepCat=trunc(d$mean.depth/10) * 10
    
  }
  
  if(!is.na(match("nlines.c",PREDS)))
  {
    par(mfcol=c(1,1),mar=c(1,1,2,1))
    barplot(table(d$nlines.c),main="n lines")
  }
  if(!is.na(match("mesh",PREDS)))
  {
    barplot(table(d$mesh),main="mesh")
  }
  
  #effort distribution
  outH <- histyear(d,Lbound=0,Rbound=max(d$km.gillnet.hours.c)*1.1,inc=10,pickvar="km.gillnet.hours.c",
                   years="year.c",varlabel="km.gillnet.hours",plots=n2mfrow(Yrs),vline=NA)
  
  #catch vs effort
  par(mfrow=c(1,1),mai=c(0.45,0.45,0.05,0.05),cex=0.85, 
      mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)  
  plot(d$km.gillnet.hours.c,d$catch.target,type="p",pch=16,col=rgb(1,0,0,1/5),
       ylim=c(0,max(d$catch.target)),xlab="km.gillnet.hours",ylab="Catch")
  abline(h=0.0,col="grey")
  
  #Exploration of Spatial distribution of data
  leftlong <- 113;  rightlong <- 129
  uplat <- -26;  downlat <- -36
  plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
  dd=subset(d,catch.target>0,select=c(lat,long,catch.target))
  names(dd)=c("Lat","Long","catch.target")
  addpoints(dd,intitle="Location of Positive catches")
  
  plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
  plotpolys(dd,leftlong,rightlong,uplat,downlat,gridon=1,leg="left",
            intitle="1 degree squares",mincount=2,namecatch="catch.target",textout = F)
  
  
  #cpue distribution by year and by month
  d$cpue=d$catch.target/d$km.gillnet.hours.c
  d$LnCE=log(d$cpue)
  cc=histyear(d,Lbound=min(d$LnCE)*1.25,Rbound=max(d$LnCE)*1.25,inc=0.2,pickvar="LnCE",
              years="year.c",varlabel="log(CPUE)",plots=n2mfrow(Yrs))
  cc=histyear(d,Lbound=min(d$LnCE)*1.25,Rbound=max(d$LnCE)*1.25,inc=0.2,pickvar="LnCE",
              years="month",varlabel="log(CPUE)",plots=n2mfrow(12))
  
  #cpue boxplots
  Cat=PREDS[which(PREDS%in%Categorical)]
  dd <- makecategorical(Cat,d)
  
  smart.par(length(PREDS),c(2.5,2.5,.1,.1),c(1,1,1,1),c(1.5,.5,0))
  for(pp in 1:length(PREDS))
  {
    x=dd[,match(c("cpue",PREDS[pp],"finyear"),names(dd))]
    if(!(is.factor(x[,2])|is.integer(x[,2]))) x[,2]=cut(x[,2],breaks=quantile(x[,2]))
    boxplot(x$cpue~x[,2],ylab="cpue",xlab=PREDS[pp],notch=F,varwidth=T)
    #"varwidth=T": box widths proportional to the square roots of the sample sizes
  }  
  
  if(Do.ggplts)
  {
    for(pp in 1:length(PREDS))
    {
      x=dd[,match(c("cpue",PREDS[pp],"finyear"),names(dd))]
      if(!(is.factor(x[,2])|is.integer(x[,2])))
      {
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue,color = finyear))+geom_point()
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue))+geom_density_2d()+
          xlab(PREDS[pp]) +geom_point(alpha=0.2,color="brown",size=1.5)
        ggplot(data=x,mapping=aes(x=x[,2],y=cpue))+geom_point()+geom_smooth(aes(color = finyear)) +facet_wrap( ~ finyear)
      }
    } 
    ggplot(dd, aes(x = vessel, y = year.c)) +  geom_jitter()
    ggplot(dd, aes(x = year.c, y = cpue)) + geom_jitter(alpha = 0.6) + facet_wrap( ~ vessel) + coord_flip()
    ggplot(dd, aes(x = cpue, fill = vessel)) +geom_histogram(bins = 25)
    if(do_pca=="YES") ggplot(dd, aes(x = dim.1,y=dim.2, color = log(cpue)))+geom_point()
  }
  
}
check.cpue=function(DATA,NAME,cl)   #function for checking cpue outliers
{
  Ktc.q=quantile(DATA$catch.target,probs=seq(0,1,.1))
  Eff.q=quantile(DATA$km.gillnet.hours.c,probs=seq(0,1,.1))
  CPUE.q=quantile(DATA$cpue.target,probs=seq(0,1,.1))
  par(mfcol=c(3,1),mai=c(.2,.5,.3,.1),mgp=c(2.5,.5,0))
  boxplot(catch.target~finyear,DATA,ylab="KG",main=NAME,col="grey80")
  abline(h=Ktc.q[6],col=cl,lwd=2)
  text(1,Ktc.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=Ktc.q[10],col=cl,lwd=2)
  text(1,Ktc.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(Ktc.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
  
  boxplot(km.gillnet.hours.c~finyear,DATA,ylab="km gn hr",col="grey80")
  abline(h=Eff.q[6],col=cl,lwd=2)
  text(1,Eff.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=Eff.q[10],col=cl,lwd=2)
  text(1,Eff.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(Eff.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
  
  boxplot(cpue.target~finyear,DATA,ylab="KG / km gn hr",col="grey80")
  abline(h=CPUE.q[6],col=cl,lwd=2)
  text(1,CPUE.q[6],"50%",pos=3,col=cl,font=2)
  abline(h=CPUE.q[10],col=cl,lwd=2)
  text(1,CPUE.q[10],"90%",pos=3,col=cl,font=2)
  legend("topright",paste("median=",round(CPUE.q[6],3),"       "),text.col=cl,cex=1.15,bty='n')
} 
fn.pred.effect <- function(DATA,PREDS) 
{
  colnames(DATA)=tolower(colnames(DATA))
  PREDS=tolower(PREDS)
  Categorical=tolower(Categorical)
  PREDS=PREDS[which(PREDS%in%colnames(DATA))]
  DATA=DATA %>% mutate_each_(funs(factor(.)),PREDS[which(PREDS%in%Categorical)])%>%
    mutate(cpue=catch.target/km.gillnet.hours.c)
  
  par(mfcol=c(3,2),mai=c(.6,.65,.3,.1),oma=c(.2,.2,.1,.1))
  hist(log(DATA$cpue),main="hist log(cpue)",xlab="log(cpue)",ylab="Count")
  
  boxcox(cpue+0.00001 ~ log(year.c), data = DATA,lambda = seq(0, 1, length = 10))
  legend("topright","Box Cox (should be small)",bty='n')
  
  # Cook distance to see outliers or overdisperse data (if distance >1)
  M1.1=glm(log(cpue+0.00001)~finyear,family=gaussian,data=DATA)
  plot(M1.1,which=4)
  legend("topright","outliers or overdispersion if distance >1",bty='n',cex=.85, text.col=2)
  
  plot(table(DATA$catch.target),type='h',xlab="Catch",ylab="Count",main="Catch zero inflation and right tail")
  
  #Outliers response var
  boxplot(DATA$cpue~DATA$finyear,main="Outliers in response var?",ylab="cpue (kg/km.gn.day)")
  
  #boxplot of response var and predictors
  smart.par(length(PREDS),c(2,2,2,.1),c(.1,.3,.1,.1),c(1.1,.35,0))
  for(d in 1:length(PREDS))
  {
    a=DATA[,match(c("cpue",PREDS[d]),names(DATA))]
    if(!is.factor(a[,2])) a[,2]=cfac(a[,2])
    plot(clog(a[,1])~a[,2],main=PREDS[d],ylab="",xlab="")
  }
  mtext("log(cpue)",side=2,line=-1,las=3,outer=T)
  
  #Covariate correlations.
  Covars=DATA%>%mutate(month=as.numeric(as.character(month)))%>%
    select(month,PREDS[which(!PREDS%in%Categorical)])
  if(ncol(Covars)>1)
  {
    pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor)
    par(mfcol=c(1,1))
    corrplot(cor(Covars), order = "AOE", addCoef.col = "grey70")
  }
  
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)   
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
fn.pred.plt=function(mod,newD,Y)
{
  Mod.pred=predict(mod,newdata=newD,type='response',se.fit=T)
  Mn=mean(Mod.pred$fit)
  plot(Y,Mod.pred$fit/Mn,pch=19,type='o',main='prediction',ylab='relative cpue')
  segments(Y,(Mod.pred$fit+1.96*Mod.pred$se.fit)/Mn,
           Y,(Mod.pred$fit-1.96*Mod.pred$se.fit)/Mn)
}
fn.mod.plt=function(MOD,MAIN,YLAB,ALL.T,DD)
{
  if(DD=='bi')
  {
    plot(MOD,main=MAIN, trans = plogis,all.terms = ALL.T, shade = TRUE,shade.col = "pink",
         shift = coef(MOD)[1],ylab=YLAB)
  }else
  {
    plot(MOD,main=MAIN, all.terms = ALL.T, shade = TRUE,shade.col = "pink",
         shift = coef(MOD)[1],ylab=YLAB)  #to see data in absolute scale (shift) and 95 CI
  }
}
fn.explr.gam.rel=function(d,Fktrs,Covars,OUT)
{
  colnames(d)=tolower(colnames(d))
  Fktrs=tolower(Fktrs)
  Covars=tolower(Covars)
  
  d=d%>%mutate(ln.effort=log(km.gillnet.hours.c),
               month.as.factor=factor(month,ordered = T))
  
  d.bi=d%>%
    mutate(catch.target=ifelse(catch.target>0,1,0),
           finyear=factor(finyear,ordered = T),
           vessel=factor(vessel),
           shots.c=factor(shots.c,ordered = T))
  
  d.pos=d%>%
    filter(catch.target>0)%>%
    mutate(cpue=catch.target/km.gillnet.hours.c,
           ln.cpue=log(cpue))%>%
    mutate(finyear=factor(finyear,ordered = T),
           vessel=factor(vessel),
           shots.c=factor(shots.c,ordered = T))
  #mutate_if(is.character,as.factor)
  
  #Exploratory graphs
  NN=1:(length(Covars)+length(Fktrs))
  nnn.i=4*NN
  plot_list=vector('list',length(NN)*4)
  for(x in 1:length(Covars))
  {
    plot_list[[nnn.i[x]-3]]=ggplot(d.pos) + aes_string(x = Covars[x], y = 'cpue')+
      geom_bin2d(bins = 100,binwidth = c(0.25, .5))+
      labs(y="Positive cpue",fill ="density") +geom_smooth(method = "lm")
    plot_list[[nnn.i[x]-2]]=ggplot(d.bi) + aes_string(x = Covars[x], y = 'catch.target')+
      geom_bin2d(bins = 100,binwidth = c(0.1, .1))+
      labs(y='Presence /absence',fill ="density")
  }
  for(x in 1:length(Fktrs))
  {
    plot_list[[nnn.i[x]-1]]=ggplot(d.pos) + aes_string(x = Fktrs[x], y = 'cpue') + 
      geom_boxplot() + labs(y="Positive cpue")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    TAB=as.data.frame(table(d.bi[,Fktrs[x]],d.bi$catch.target))
    plot_list[[nnn.i[x]]]=ggplot(TAB) + aes(x = Var1, y = Freq,fill=factor(Var2)) +
      geom_bar(stat='identity',position="stack")+
      labs( x = Fktrs[x],fill ="")+ 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  #Interaction
  # qplot(x = year.c, y = cpue, data = d.pos, color = factor(blockx)) +
  #    geom_smooth(method = "lm")
  
  plot_list <- plot_list[which(!sapply(plot_list, is.null))]
  multi.page <-ggarrange(plotlist=plot_list, nrow = 2, ncol = 1)
  ggexport(multi.page, filename = paste(hndl.expl,OUT,"_GAM.pdf",sep=''))
  
  
  pdf(paste(hndl.expl,OUT,"_GAM_explore_year.month.vessel.pdf",sep=''))
  
  #Tweedie vs normal
  a=d%>%mutate(finyear=factor(finyear,ordered = T,levels=sort(unique(finyear))),
               cpue=catch.target/km.gillnet.hours.c)
  Mod1=gam(cpue~finyear,dat=a,family=tw,method="REML")  
  Mod2=gam(cpue~finyear,dat=a,method="REML")
  
  par(mfcol=c(2,3),mar=c(3,4,1,1),mgp=c(2,.9,0))
  
  fn.mod.plt(Mod1,"","cpue",TRUE,'pos');legend('top',"Mod1 Tweedie",bty='n',cex=1.5)
  fn.mod.plt(Mod2,"","cpue",TRUE,'pos');legend('top',"Mod2 Normal",bty='n',cex=1.5)
  
  New.d=data.frame(finyear=factor(levels(a$finyear)))%>%
    mutate(year.c=as.numeric(substr(finyear,1,4)))
  fn.pred.plt(mod=Mod1,newD=New.d,Y=New.d$year.c)
  fn.pred.plt(mod=Mod2,newD=New.d,Y=New.d$year.c)
  
  plot(coef(Mod1)[-1],coef(Mod2)[-1],pch=19,cex=2,xlab="coef mod1",ylab="coef mod2")
  lines(coef(Mod1)[-1],coef(Mod1)[-1],lwd=2,col=2)
  
  
  #Model year as an autocorrelated continuous variable
  #pos part
  Mod1=gam(ln.cpue~s(year.c,bs='gp'),dat=d.pos,method="REML")  #specify a gaussian process
  Mod2=gam(ln.cpue~s(year.c),dat=d.pos,method="REML")
  Mod3=gam(ln.cpue~finyear,dat=d.pos,method="REML")
  
  #bi part
  Mod1.bi=gam(catch.target~s(year.c,bs='gp')+offset(ln.effort),dat=d.bi,family=binomial,method="REML")  #specify a gaussian process
  Mod2.bi=gam(catch.target~s(year.c)+offset(ln.effort),dat=d.bi,family=binomial,method="REML")
  Mod3.bi=gam(catch.target~finyear+offset(ln.effort),dat=d.bi,family=binomial,method="REML")
  
  
  par(mfcol=c(3,2),mar=c(3,4,1,1))
  #pos part
  fn.mod.plt(Mod1,"Mod1 s(year,bs='gp')","cpue",FALSE,'pos')
  fn.mod.plt(Mod2,"Mod2 s(year)","cpue",FALSE,'pos')
  fn.mod.plt(Mod3,"Mod2 s(year)","cpue",TRUE,'pos')
  legend('top',"Mod3 finyear",bty='n',cex=1.5)
  
  New.d=data.frame(finyear=factor(levels(d.pos$finyear)))%>%
    mutate(year.c=as.numeric(substr(finyear,1,4)))
  fn.pred.plt(mod=Mod1,newD=New.d,Y=New.d$year.c)
  fn.pred.plt(mod=Mod2,newD=New.d,Y=New.d$year.c)
  fn.pred.plt(mod=Mod3,newD=New.d,Y=New.d$year.c)
  
  grid.newpage()
  mod.sumery=rbind(tidy(Mod1),tidy(Mod2),tidy(Mod3))
  class(mod.sumery) <- "data.frame"
  rownames(mod.sumery)=paste('Mod',1:nrow(mod.sumery),sep='')
  grid.table(mod.sumery)
  
  AIC.tab=AICtab(Mod1,Mod2,Mod3)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  #bi part
  par(mfcol=c(3,2),mar=c(3,4,1,1))
  fn.mod.plt(Mod1.bi,"Mod1.bi s(year,bs='gp')","catch",FALSE,'bi')
  fn.mod.plt(Mod2.bi,"Mod2.bi s(year)","catch",FALSE,'bi')
  fn.mod.plt(Mod3.bi,"Mod2.bi s(year)","catch",TRUE,'bi')
  legend('top',"Mod3.bi finyear",bty='n',cex=1.5)
  
  New.d=data.frame(finyear=factor(levels(d.bi$finyear)))%>%
    mutate(year.c=as.numeric(substr(finyear,1,4)),
           ln.effort=mean(d.bi$ln.effort))
  fn.pred.plt(mod=Mod1.bi,newD=New.d,Y=New.d$year.c)
  fn.pred.plt(mod=Mod2.bi,newD=New.d,Y=New.d$year.c)
  fn.pred.plt(mod=Mod3.bi,newD=New.d,Y=New.d$year.c)
  
  grid.newpage()
  mod.sumery=rbind(tidy(Mod1.bi),tidy(Mod2.bi),tidy(Mod3.bi))
  class(mod.sumery) <- "data.frame"
  rownames(mod.sumery)=paste('Mod',1:nrow(mod.sumery),sep='')
  grid.table(mod.sumery)
  
  AIC.tab=AICtab(Mod1.bi,Mod2.bi,Mod3.bi)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  
  #Model month as smoothing term
  #pos part
  Mod1=gam(ln.cpue~s(month,k=12,bs='cc'),dat=d.pos,method="REML")  
  Mod2=gam(ln.cpue~s(month,k=12),dat=d.pos,method="REML")
  Mod3=gam(ln.cpue~month.as.factor,dat=d.pos,method="REML")
  
  #bi part
  Mod1.bi=gam(catch.target~s(month,k=12,bs='cc')+offset(ln.effort),dat=d.bi,family=binomial,method="REML")
  Mod2.bi=gam(catch.target~s(month,k=12)+offset(ln.effort),dat=d.bi,family=binomial,method="REML")
  Mod3.bi=gam(catch.target~month.as.factor+offset(ln.effort),dat=d.bi,family=binomial,method="REML")
  
  
  par(mfcol=c(3,2),mar=c(3,4,1,1))
  #pos part
  fn.mod.plt(Mod1,"Mod1 s(month,k=12,bs='cc')","cpue",FALSE,'pos')
  fn.mod.plt(Mod2,"Mod2 s(month,k=12)","cpue",FALSE,'pos')
  fn.mod.plt(Mod3,"month.as.factor","cpue",TRUE,'pos')
  legend('top',"Mod3_month.as.factor",bty='n',cex=1.5)
  
  New.d=data.frame(month.as.factor=factor(levels(d.pos$month.as.factor)))%>%
    mutate(month=as.numeric(as.character(month.as.factor)))
  fn.pred.plt(mod=Mod1,newD=New.d,Y=New.d$month)
  fn.pred.plt(mod=Mod2,newD=New.d,Y=New.d$month)
  fn.pred.plt(mod=Mod3,newD=New.d,Y=New.d$month)
  
  grid.newpage()
  mod.sumery=rbind(tidy(Mod1),tidy(Mod2),tidy(Mod3))
  class(mod.sumery) <- "data.frame"
  rownames(mod.sumery)=paste('Mod',1:nrow(mod.sumery),sep='')
  grid.table(mod.sumery)
  
  AIC.tab=AICtab(Mod1,Mod2,Mod3)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  #bi part
  par(mfcol=c(3,2),mar=c(3,4,1,1))
  fn.mod.plt(Mod1.bi,"Mod1.bi s(month,k=12,bs='cc')","catch",FALSE,'bi')
  fn.mod.plt(Mod2.bi,"Mod2.bi s(month,k=12)","catch",FALSE,'bi')
  fn.mod.plt(Mod3.bi,"month.as.factor","catch",TRUE,'bi')
  legend('top',"Mod3.bi month.as.factor",bty='n',cex=1.5)
  
  New.d=data.frame(month.as.factor=factor(levels(d.bi$month.as.factor)))%>%
    mutate(month=as.numeric(as.character(month.as.factor)),
           ln.effort=mean(d.bi$ln.effort))
  fn.pred.plt(mod=Mod1.bi,newD=New.d,Y=New.d$month)
  fn.pred.plt(mod=Mod2.bi,newD=New.d,Y=New.d$month)
  fn.pred.plt(mod=Mod3.bi,newD=New.d,Y=New.d$month)
  
  grid.newpage()
  mod.sumery=rbind(tidy(Mod1.bi),tidy(Mod2.bi),tidy(Mod3.bi))
  class(mod.sumery) <- "data.frame"
  rownames(mod.sumery)=paste('Mod',1:nrow(mod.sumery),sep='')
  grid.table(mod.sumery)
  
  AIC.tab=AICtab(Mod1.bi,Mod2.bi,Mod3.bi)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  
  #Model vessel as random effect
  
  #pos part
  Mod1=gam(ln.cpue~finyear+s(vessel,bs='re'),dat=d.pos,method="REML")  
  Mod2=gam(ln.cpue~finyear+vessel,dat=d.pos,method="REML")
  
  #bi part
  Mod1.bi=bam(catch.target~finyear+s(vessel,bs='re')+offset(ln.effort),dat=d.bi,family=binomial,method="fREML")
  Mod2.bi=bam(catch.target~finyear+vessel+offset(ln.effort),dat=d.bi,family=binomial,method="fREML")
  
  par(mfcol=c(2,1),mar=c(3,4,1,1))
  New.d=data.frame(finyear=factor(levels(d.pos$finyear)))%>%
    mutate(vessel=factor(names(sort(-table(d.pos$vessel)))[1],levels=levels(d.pos$vessel)))
  fn.pred.plt(mod=Mod1,newD=New.d,Y=as.numeric(substr(New.d$finyear,1,4)))
  fn.pred.plt(mod=Mod2,newD=New.d,Y=as.numeric(substr(New.d$finyear,1,4)))
  
  AIC.tab=AICtab(Mod1,Mod2)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  #bi part
  par(mfcol=c(2,1),mar=c(3,4,1,1))
  New.d=data.frame(finyear=factor(levels(d.bi$finyear)))%>%
    mutate(vessel=factor(names(sort(-table(d.bi$vessel)))[1],levels=levels(d.bi$vessel)),
           ln.effort=mean(d.bi$ln.effort))
  fn.pred.plt(mod=Mod1.bi,newD=New.d,Y=as.numeric(substr(New.d$finyear,1,4)))
  fn.pred.plt(mod=Mod2.bi,newD=New.d,Y=as.numeric(substr(New.d$finyear,1,4)))
  
  AIC.tab=AICtab(Mod1.bi,Mod2.bi)
  class(AIC.tab) <- "data.frame"
  grid.newpage()
  grid.table(AIC.tab)
  
  dev.off()
  
  
  #plotting stuff
  #plot(Mod, rug = TRUE, residuals = TRUE, pch = 1, cex = 1) # see residuals
  #plot(Mod, seWithMean = TRUE, shift = coef(Mod)[1])  #shows the partial effect and intercept uncertainty (better approach) in an interpretable scale
  #plot(Mod,main="Mod", page=1, all.terms = TRUE, rug = TRUE, shade = TRUE,shade.col = "pink")
  
  #better summary tables
  #tidy(Mod)
  #glance(Mod)
  #augment(Mod)
}

# CREATE SPECIES DATA SETS FOR STANDARDISATIONS ----------------------------------------------
fn.cpue.data=function(Dat,EffrrT,sp)  
{
  TAB=with(subset(Dat,SPECIES%in%sp),unique(YEAR.c))
  if(length(TAB)>=N.keep)
  {
    #aggregate records by Same return (drop issues with Condition and Bioregion...)
    Dat=Dat%>%group_by(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX,Boundary.blk,SPECIES,SNAME,YEAR.c,
                       LAT,LONG,Same.return,TYPE.DATA,zone,Reporter,Sch.or.DogS,
                       Temperature,Temp.res,SOI,Freo,Freo_lag6,Freo_lag12,Yrs.of.experience)%>%
      summarise(LIVEWT = sum(LIVEWT),
                LIVEWT.c = sum(LIVEWT.c))%>%
      data.frame()
    
    #add effort
    Ids=match(c("LAT","LONG"),names(EffrrT))
    Dat=Dat%>%left_join(EffrrT[,-Ids],by=c("BLOCKX","FINYEAR","MONTH","VESSEL"))
    
    #consider effort reporter
    Dat$Reporter=with(Dat,ifelse(Eff.Reporter=="bad","bad",Reporter))  
    
    #remove records with NA effort
    Dat=subset(Dat,!(is.na(Km.Gillnet.Days.c) | is.na(Km.Gillnet.Hours.c)))
    
    #remove school shark or dogfish shots if not the target species
    idd=which(sp%in%c(17008,20000))
    if(length(idd)==0) Dat=subset(Dat,!(Sch.or.DogS=="Yes"))
    
    #keep records from first year with data as some species (e.g. Sandbars) 
    # didn't have a code for reporting early on
    TAB=with(subset(Dat,SPECIES%in%sp),table(YEAR.c))
    Dat=subset(Dat,YEAR.c>=as.numeric(names(TAB[1])))
    
    #for greynurse, drop years post protection
    idd=which(sp%in%8001)
    if(length(idd)==1 & length(sp)<3)
    {
      Greyn.yrs=sort(unique(Dat$FINYEAR))
      Greyn.yrs=Greyn.yrs[1:match(Greynurse.protection,Greyn.yrs)]
      Dat=subset(Dat,FINYEAR%in%Greyn.yrs)
    }
  }else
  {
    Dat=NULL
  }
  return(Dat)
}
fn.cpue.data.daily=function(Dat,EffrrT,sp)
{
  TAB=with(subset(Dat,SPECIES%in%sp),unique(YEAR.c))
  if(length(TAB)>=N.keep)
  {
    Dat=Dat%>%group_by(date,Same.return.SNo,TSNo,Same.return,FINYEAR,MONTH,
                       VESSEL,METHOD,BLOCKX,block10,SPECIES,SNAME,YEAR.c,LAT,LONG,
                       TYPE.DATA,zone,Reporter,Sch.or.DogS,ZnID,
                       Temperature,Temp.res,SOI,Freo,Freo_lag6,Freo_lag12,Lunar,
                       BoatName,MastersName,Yrs.of.experience)%>%
      summarise(LIVEWT = sum(LIVEWT),
                LIVEWT.c = sum(LIVEWT.c),
                nfish = sum(nfish))%>%
      data.frame()
    
    #add effort
    Ids=match(c("BLOCKX","FINYEAR","MONTH","LAT","LONG","VESSEL","block10"),names(EffrrT))
    Dat=Dat%>%left_join(EffrrT[,-Ids],by=c("Same.return.SNo"))
    
    #consider effort reporter
    Dat$Reporter=with(Dat,ifelse(Eff.Reporter=="bad","bad",Reporter))  
    
    #remove records with NA effort
    Dat=subset(Dat,!(is.na(Km.Gillnet.Days.c) | is.na(Km.Gillnet.Hours.c)))
    
    #remove school shark or dogfish shots if not the target species
    idd=which(sp%in%c(17008,20000))
    if(length(idd)==0) Dat=subset(Dat,!(Sch.or.DogS=="Yes"))
  }else
  {
    Dat=NULL
  }
  return(Dat)
}

# VESSELS & BLOCKS----------------------------------------------
fn.see.all.yrs.ves.blks=function(a,SP,NM,what,Ves.sel.BC,Ves.sel.sens,BLK.sel.BC,BLK.sel.sens,Min.ktch)
{
  All.ves=unique(as.character(a$VESSEL))
  All.blk=unique(as.character(a$BLOCKX))
  a=subset(a,Reporter=="good")
  dddd=subset(a,SPECIES%in%SP)
  if(nrow(dddd)>10)
  {
    CATCH.sp=dddd %>% group_by(YEAR.c,VESSEL)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
      spread(VESSEL, LIVEWT.c) %>%
      arrange(YEAR.c) %>%
      data.frame()
    Vess=names(CATCH.sp)[2:ncol(CATCH.sp)]
    Vess=chartr(".", " ", Vess)
    Yrs=CATCH.sp$YEAR.c
    Z=as.matrix(CATCH.sp[,-1])
    
    #Step 1. Select vessels with > X years of records of a minimum catch
    ZZ=Z
    ZZ[ZZ<Min.ktch]=NA
    ZZ[ZZ>=Min.ktch]=1
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    
    pdf(paste(handl_OneDrive("Analyses/Catch and effort/Outputs/Kept_blocks_vessels/Vessel_pos_records_by_yr/"),paste(NM,what,sep=""),".pdf",sep="")) 
    
    #Ves.sel.BC
    par(mar=c(3,3.5,.8,.8))
    WHICh=which(Yrs.with.ktch>Ves.sel.BC)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.BC=Vess[WHICh]
    if(length(Ves.BC)==1)
    {
      plot.new()
      legend("center","only 1 vessel selected",bty='n',cex=2)
    }
    
    if(length(Ves.BC)>1)
    {
      ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
      Z.this=Z.this[,ID.sort]
      if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
      Ves.BC=Ves.BC[ID.sort]
      image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
      axis(1,1:length(Yrs),Yrs)
      axis(2,1:ncol(Z.this),Ves.BC,las=1,cex.axis=.9)
      legend("top",paste("vessels with >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
    }
    Drop.ves=All.ves[which(!All.ves%in%Ves.BC)]
    
    #Number of records for selected vessels
    Tab.sel.ves=dddd%>%
      filter(VESSEL%in%Ves.BC)%>%
      group_by(Same.return,VESSEL,FINYEAR)%>%
      summarise(tot=sum(LIVEWT.c))%>%
      mutate(FINYEAR=substr(FINYEAR,1,4))%>%
      group_by(VESSEL,FINYEAR)%>%
      summarise(n=n())%>%
      spread(FINYEAR,n,fill='')%>%
      data.frame
    colnames(Tab.sel.ves)[2:ncol(Tab.sel.ves)]=substr(colnames(Tab.sel.ves)[2:ncol(Tab.sel.ves)],2,10)
    grid.newpage()
    grid.draw(tableGrob(Tab.sel.ves,rows = NULL,theme=ttheme_minimal(base_size = 7) ))
    
    
    #Ves.sel.sens
    WHICh=which(Yrs.with.ktch>Ves.sel.sens)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.Sens=Vess[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Ves.Sens=Ves.Sens[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Ves.Sens,las=1,cex.axis=.6)
    legend("top",paste("vessels with >=",Ves.sel.sens, "years of records"),bty='n')
    
    #plot CPUEs
    a$CPUE.km.gn.day=a$LIVEWT.c/a$Km.Gillnet.Days.c
    a$CPUE.km.gn.h=a$LIVEWT.c/a$Km.Gillnet.Hours.c
    a.mean.cpue.km.day_all=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,SPECIES%in%SP),mean)
    a.mean.cpue.km.h_all=aggregate(CPUE.km.gn.h~YEAR.c,subset(a, SPECIES%in%SP),mean)
    a.mean.cpue.km.day_Sens=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES%in%SP),mean)
    a.mean.cpue.km.h_Sens=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.Sens& SPECIES%in%SP),mean)
    a.mean.cpue.km.day_BC=aggregate(CPUE.km.gn.day~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES%in%SP),mean)
    a.mean.cpue.km.h_BC=aggregate(CPUE.km.gn.h~YEAR.c,subset(a,VESSEL%in%Ves.BC& SPECIES%in%SP),mean)
    
    #kmgday
    par(mar=c(3,3,.5,4),mgp=c(2,.7,0))
    Yrs=a.mean.cpue.km.day_Sens$YEAR.c
    plot(a.mean.cpue.km.day_all$YEAR.c,a.mean.cpue.km.day_all$CPUE.km.gn.day,ylab="",xlab="")
    points(a.mean.cpue.km.day_Sens$YEAR.c,a.mean.cpue.km.day_Sens$CPUE.km.gn.day,pch=19,col=2)
    if(nrow(a.mean.cpue.km.day_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.day_BC$CPUE.km.gn.day,pch=19,col=3)
    
    #kmgday  
    par(new = T)
    plot(a.mean.cpue.km.h_all$YEAR.c,a.mean.cpue.km.h_all$CPUE.km.gn.h,ylab=NA, axes=F,xlab=NA,pch=0,cex=2)
    points(a.mean.cpue.km.h_Sens$YEAR.c,a.mean.cpue.km.h_Sens$CPUE.km.gn.h,pch=15,col=2,cex=2)
    if(nrow(a.mean.cpue.km.h_BC)==length(Yrs))points(Yrs,a.mean.cpue.km.h_BC$CPUE.km.gn.h,pch=15,col=3,cex=2)
    axis(side = 4)
    mtext("Financial year",1,line=2)
    mtext("Nominal CPUE (Kg/km.gn.day)",2,line=2)
    mtext("Nominal CPUE (Kg/km.gn.hour)",4,line=2)
    legend("top",c("Kg/km.gn.day","Kg/km.gn.hour"),bty='n',pch=c(0,19))
    legend("bottomleft",c("all",paste(Ves.sel.sens,"y"),paste(Ves.sel.BC,"y")),bty='n',pch=c(0,19,19),col=c(1,2,3))
    
    
    # step 2. For selected vessels, plot number of blocks by year
    #Ves.BC
    #all blocks
    AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.BC) %>%
      group_by(YEAR.c,BLOCKX)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
      spread(BLOCKX, LIVEWT.c) %>%
      arrange(YEAR.c) %>%
      data.frame()
    BLOCs=substr(names(AA)[2:ncol(AA)],2,6)
    Yrs=AA$YEAR.c
    Z=as.matrix(AA[,-1])
    ZZ=Z
    ZZ[ZZ>0]=1
    ZZZ=ZZ
    ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
    ZZZ=ZZZ[,ID.sort]
    if(!is.matrix(ZZZ)) ZZZ=t(as.matrix(ZZZ))
    par(mar=c(3,3.5,.8,.8))
    image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
    legend("top",paste("all blocks for vessels >=",Ves.sel.BC, "years of records and >",Min.ktch,"kg per year"),bty='n')
    
    #blocks with > bc records
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    WHICh=which(Yrs.with.ktch>BLK.sel.BC)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.BC=BLOCs[WHICh]
    if(length(Blks.BC)==0)
    {
      plot.new()
      legend("center","no blocks selected",bty='n',cex=2)
    }
    
    if(length(Blks.BC)>0)
    {
      ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
      Z.this=Z.this[,ID.sort]
      if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
      Blks.BC=Blks.BC[ID.sort]
      image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
      axis(1,1:length(Yrs),Yrs)
      axis(2,1:ncol(Z.this),Blks.BC,las=1,cex.axis=.9)
      legend("top",paste("blocks with >=",BLK.sel.BC, "years of records for vessels >=",Ves.sel.BC,"years of records and >",Min.ktch,"kg per year"),
             cex=0.75,bty='n')
      
    }
    Drop.blks=All.blk[which(!All.blk%in%Blks.BC)]  
    
    #Ves.Sens
    #all blocks
    AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.Sens) %>%
      group_by(YEAR.c,BLOCKX)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
      spread(BLOCKX, LIVEWT.c) %>%
      arrange(YEAR.c) %>%
      data.frame()
    BLOCs=substr(names(AA)[2:ncol(AA)],2,6)
    Yrs=AA$YEAR.c
    Z=as.matrix(AA[,-1])
    ZZ=Z
    ZZ[ZZ>0]=1
    ZZZ=ZZ
    ID.sort=match(names(sort(colSums(ZZZ,na.rm=TRUE))),colnames(ZZZ))
    ZZZ=ZZZ[,ID.sort]
    if(!is.matrix(ZZZ)) ZZZ=t(as.matrix(ZZZ))
    par(mar=c(3,3.5,.8,.8))
    image(x=1:length(Yrs),y=1:ncol(ZZZ),ZZZ,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:length(BLOCs),BLOCs[ID.sort],las=1,cex.axis=.9)
    legend("top",paste("all blocks for vessels >=",Ves.sel.sens, "years of records"),bty='n')
    
    #blocks with > Sens records
    Yrs.with.ktch=colSums(ZZ,na.rm=T)
    WHICh=which(Yrs.with.ktch>BLK.sel.sens)
    Z.this=ZZ[,WHICh]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.Sens=BLOCs[WHICh]
    ID.sort=match(names(sort(colSums(Z.this,na.rm=TRUE))),colnames(Z.this))
    Z.this=Z.this[,ID.sort]
    if(!is.matrix(Z.this)) Z.this=t(as.matrix(Z.this))
    Blks.Sens=Blks.Sens[ID.sort]
    image(x=1:length(Yrs),y=1:ncol(Z.this),Z.this,xaxt='n',yaxt='n',ann=F)
    axis(1,1:length(Yrs),Yrs)
    axis(2,1:ncol(Z.this),Blks.Sens,las=1,cex.axis=.9)
    legend("top",paste("blocks with >=",BLK.sel.sens, "years of records for vessels >=",Ves.sel.sens,"years of records of",SP),bty='n')
    dev.off()
    
    Drop.blks_10=Blks.BC_10=NULL
    if(what==".daily")
    {
      AA=subset(a,SPECIES%in%SP & VESSEL%in%Ves.BC) %>%
        group_by(YEAR.c,block10)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c)) %>%
        spread(block10, LIVEWT.c) %>%
        arrange(YEAR.c) %>%
        data.frame()
      AA=as.matrix(AA[,-1])
      AA[AA>0]=1
      Yrs.with.ktch=colSums(AA,na.rm=T)
      Blks.BC_10=substr(names(which(Yrs.with.ktch>BLK.sel.BC)),2,50)
      Drop.blks_10=unique(a$block10)[which(!unique(a$block10)%in%as.numeric(Blks.BC_10))]
    }
    return(list(Ves.BC=Ves.BC, Ves.Sens=Ves.Sens, Blks.BC=Blks.BC,Blks.BC_10=Blks.BC_10, Blks.Sens=Blks.Sens,
                Drop.ves=Drop.ves, Drop.blks=Drop.blks,Drop.blks_10=Drop.blks_10))
    
  }
}
fn.comp.used_not.used=function(blk.used,blk.not.used,vsl.used,vsl.not.used,NMS.arg)
{
  barplot(prop.table(matrix(c(length(blk.used),length(blk.not.used),
                              length(vsl.used),length(vsl.not.used)),ncol=2,nrow=2),2),
          col=c("chartreuse3","brown1"),names.arg=NMS.arg)
}

# CONSTRUCT WIDE DATABASE FOR STANDARDISATIONS ----------------------------------------------
Effort.data.fun=function(DATA,target,ktch)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  if(nrow(DATA)>0)
  {
    #remove nonsense lat
    DATA=subset(DATA,LAT>=(-36))
    
    #calculate effort (only need max effort)
    Match.these.eff=match(These.efforts,names(DATA))
    Effort.data1=DATA[,Match.these.eff]
    Effort.data=Effort.data1%>%
      group_by(zone,FINYEAR,Same.return,MONTH,BLOCKX,SHOTS.c,HOURS.c)%>%
      summarise(Km.Gillnet.Days.c = max(Km.Gillnet.Days.c),
                Km.Gillnet.Hours.c = max(Km.Gillnet.Hours.c))%>%
      data.frame()
    
    #target species catch 
    #note: catch targeted at other species: pointless as these are multiple trips combined in one month
    ID=match(c(ktch),colnames(DATA))
    DATA$Catch.Target=with(DATA,ifelse(SPECIES%in%target,DATA[,ID],0))
    DATA$Catch.Total=with(DATA,ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),DATA[,ID],0))
    
    #reshape catch data
    TABLE=DATA%>%group_by(MONTH,FINYEAR,BLOCKX,VESSEL,Same.return,LAT,LONG,YEAR.c)%>%
      summarise(Catch.Target = sum(Catch.Target,na.rm=T),
                Catch.Total = sum(Catch.Total,na.rm=T))
    Enviro=DATA%>%group_by(MONTH,FINYEAR,BLOCKX)%>%
      summarise(Temperature=mean(Temperature),
                Temp.res=mean(Temp.res),
                Freo=mean(Freo,na.rm=T),
                SOI=mean(SOI,na.rm=T))
    TABLE=TABLE%>%left_join(Enviro,by=c("FINYEAR","MONTH","BLOCKX"))%>%
      arrange(FINYEAR,MONTH,BLOCKX) %>%
      data.frame()
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    
    #merge catch and effort
    dat=TABLE%>%left_join(Effort.data,by=c("Same.return","FINYEAR","MONTH","BLOCKX"))
    
    #Add mesh size                                
    # d=subset(DATA,select=c(Same.return,mesh))
    # d=d[!duplicated(d$Same.return),]
    # dat=dat%>%left_join(d,by=c("Same.return"))
    
  }else
  {
    dat=DATA
    prop.with.catch=0
  }
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}
Effort.data.fun.daily=function(DATA,target,ktch,Aggregtn)
{
  #remove record if no effort data
  ID=which(DATA$Km.Gillnet.Days.c==0) #no effort records
  ID=c(ID,which(is.na(DATA$Km.Gillnet.Days.c))) #NA records
  if(length(ID)>0)DATA=DATA[-ID,]
  
  if(nrow(DATA)>0)
  {
    #remove nonsense lat
    DATA=subset(DATA,LAT>=(-36))
    
    #unique lat and long per Same.return.SNo
    DATA=DATA%>%
      group_by(Same.return.SNo)%>%
      mutate(LAT=mean(LAT,na.rm=T),
             LONG=mean(LONG,na.rm=T))%>%
      ungroup()
    
    #calculate effort
    Match.these.eff=match(These.efforts.daily,names(DATA))
    Effort.data1=DATA[,Match.these.eff]
    
    #aggregate at shot level
    if(Use.Date=="NO")
    {
      #max effort by Sno, DsNo & TSNo
      Effort.data=Effort.data1%>%      
        group_by(zone,FINYEAR,Same.return.SNo,MONTH,BLOCKX,block10,shots.c,BoatName,MastersName)%>%
        summarise(Km.Gillnet.Days.c = max(Km.Gillnet.Days.c),
                  Km.Gillnet.Hours.c = max(Km.Gillnet.Hours.c))%>%
        data.frame()
      
      #aggregate at TSNo if required
      if(Aggregtn=="TSNo")
      {
        Effort.data$TSNo=word(Effort.data$Same.return.SNo,3)
        Effort.data=aggregate(cbind(Km.Gillnet.Days.c,Km.Gillnet.Hours.c)~zone+
                                FINYEAR+TSNo+MONTH+BLOCKX+BoatName+MastersName,Effort.data,sum)
      }
    }
    
    #target species and other main species catch 
    DATA=DATA%>%
      mutate(Catch.Target=ifelse(SPECIES%in%target,!!sym(ktch),0),
             Catch.Gummy=ifelse(SPECIES==17001,!!sym(ktch),0),
             Catch.Whiskery=ifelse(SPECIES==17003,!!sym(ktch),0),
             Catch.Dusky=ifelse(SPECIES%in%c(18003),!!sym(ktch),0),
             Catch.Sandbar=ifelse(SPECIES==18007,!!sym(ktch),0),
             Catch.Groper=ifelse(SPECIES%in%c(384002),!!sym(ktch),0),
             Catch.Snapper=ifelse(SPECIES%in%c(353001),!!sym(ktch),0),
             Catch.Blue_mor=ifelse(SPECIES%in%c(377004),!!sym(ktch),0),
             Catch.Total=ifelse(SPECIES%in%c(5001:24900,25000:31000,188000:599001),!!sym(ktch),0))

    #reshape catch data
    if(Use.Date=="NO")
    {
      if(Aggregtn=="SNo") 
      {
        TABLE=DATA%>%group_by(MONTH,FINYEAR,BLOCKX,block10,VESSEL,Same.return.SNo,date,LAT,LONG,YEAR.c,Lunar)%>%
          summarise(Catch.Target = sum(Catch.Target,na.rm=T),
                    Catch.Gummy=sum(Catch.Gummy,na.rm=T),
                    Catch.Whiskery=sum(Catch.Whiskery,na.rm=T),
                    Catch.Dusky=sum(Catch.Dusky,na.rm=T),
                    Catch.Sandbar=sum(Catch.Sandbar,na.rm=T),
                    Catch.Groper=sum(Catch.Groper,na.rm=T),
                    Catch.Snapper=sum(Catch.Snapper,na.rm=T),
                    Catch.Blue_mor=sum(Catch.Blue_mor,na.rm=T),
                    Catch.Total=sum(Catch.Total,na.rm=T))
        Enviro=DATA%>%group_by(MONTH,FINYEAR,BLOCKX)%>%
          summarise(Temperature=mean(Temperature),
                    Temp.res=mean(Temp.res),
                    Freo=mean(Freo),
                    SOI=mean(SOI,na.rm=T))
        TABLE=TABLE%>%left_join(Enviro,by=c("FINYEAR","MONTH","BLOCKX"))    %>%
          arrange(Same.return.SNo,FINYEAR,MONTH,BLOCKX) %>%
          data.frame()
      }
      #aggregating by trip
      if(Aggregtn=="TSNo")   
      {
        TABLE=aggregate(cbind(Catch.Target,Catch.Gummy,Catch.Whiskery,Catch.Dusky,Catch.Sandbar,
                              Catch.Groper,Catch.Snapper,Catch.Blue_mor,Catch.Dhufish,
                              Catch.Other.shrk,Catch.Other.scalefish,
                              Catch.non_indicators,Catch.Total)~MONTH+FINYEAR+BLOCKX+VESSEL+
                          TSNo+YEAR.c,data=DATA,sum,na.rm=T)
        xx=subset(DATA,select=c(BLOCKX,LAT,LONG))
        xx=xx[!duplicated(xx$BLOCKX),]
        xx$LAT=round(xx$LAT)
        xx$LONG=round(xx$LONG)
        TABLE=merge(TABLE,xx,by="BLOCKX",all.x=T)
        TABLE=TABLE[order(TABLE$TSNo,TABLE$FINYEAR,TABLE$MONTH,TABLE$BLOCKX),]
      }
    }
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    
    #merge catch and effort
    if(Aggregtn=="SNo") dat=TABLE%>%left_join(Effort.data,by=c("Same.return.SNo","FINYEAR","MONTH","BLOCKX","block10"))
    if(Aggregtn=="TSNo") dat=TABLE%>%left_join(Effort.data,by=c("TSNo","FINYEAR","MONTH","BLOCKX","block10"))
    
    
    #Add mesh size, shots, depth and nlines for each session
    if(Aggregtn=="SNo")
    {
      d=subset(DATA,select=c(Same.return.SNo,VESSEL,mesh,nlines.c,Mean.depth,hours.c))
      d=d[!duplicated(paste(d$Same.return.SNo,d$VESSEL)),]
      dat=dat%>%left_join(d,by=c("Same.return.SNo","VESSEL"))
    }
    
  }else
  {
    dat=DATA
    prop.with.catch=0
  }
  
  return(list(dat=dat,prop.with.catch=prop.with.catch))
}


# GET NUMBER OF SPECIES BY DAILY RECORD ----------------------------------------------
fn.species.per.session=function(d,target)
{
  d=d%>%
    group_by(Same.return.SNo,SNAME)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  a=d%>%
    group_by(Same.return.SNo)%>%
    tally()
  p1=a%>%
    ggplot(aes(n))+
    geom_bar()+xlab('Number of species')+ggtitle(target)+ylab('Number of sessions')
  
  x=a%>%filter(n==1)%>%pull(Same.return.SNo)
  p_one.species=d%>%
    filter(Same.return.SNo%in%x)%>%
    group_by(SNAME)%>%
    tally()%>%
    arrange(n)
  p_one.species=p_one.species%>%
    mutate(SNAME=factor(SNAME,levels=unique(p_one.species$SNAME)))%>%
    ggplot(aes(x=SNAME,y=n))+
    geom_bar(stat = 'identity')+
    ylab('Number of sessions')+xlab('')+
    coord_flip()+
    labs(subtitle='records with only 1 species reported')
  
  return(ggarrange(p1,p_one.species,ncol=1))  
}
# PROPORTION NO CATCH THRU TIME BY FISHER TO ID FISHING EFFICIENCY CREEP ----------------------------------------------
fn.prop.0.catch.by.fisher=function(d,explained.ktch,NM,series)
{
  CUM.ktch=d%>%
    group_by(Ves.var)%>%
    summarise(Ktch=sum(Catch.Target))%>%
    ungroup()%>%
    arrange(-Ktch)%>%
    mutate(CumKtch=cumsum(Ktch),
           Percent=CumKtch/sum(Ktch))%>%
    filter(Percent<=explained.ktch)
  
  d1=d%>%
    filter(Ves.var%in%CUM.ktch$Ves.var)%>%
    group_by(Ves.var,time.var)%>%
    summarise(cpue=mean(cpue))%>%
    mutate(yr=floor(time.var))
  
  d1.bin=d1%>%mutate(Target=ifelse(cpue>0,1,0))
  N.row=ceiling(length(unique(d1$Ves.var))/10)
  TitlE=paste("Vessels explaining",100*explained.ktch,"% of the catch")
  p=d1%>%
    group_by(Ves.var,yr)%>%
    summarise(Target=mean(cpue))%>%
    ggplot(aes(yr,Target, group = Ves.var, color = Ves.var))+ 
    geom_point() + 
    geom_line() +
    ylab('Mean kg/km gn h') + ggtitle(TitlE)+
    theme(legend.title = element_blank(),
          legend.position = 'top')+
    guides(color=guide_legend(nrow=N.row,byrow=TRUE))
  
  p.bin=d1.bin%>%
    group_by(Ves.var,yr)%>%
    summarise(Target=mean(Target))%>%
    ggplot(aes(yr,Target, group = Ves.var, color = Ves.var))+ 
    geom_point() + 
    geom_line() +
    ylab('Proportion positive catch') +
    theme(legend.title = element_blank(),
          legend.position = 'top')+
    guides(color=guide_legend(nrow=N.row,byrow=TRUE))
  
  print(ggarrange(p,p.bin, ncol=1,common.legend = T))
  Out.name=handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/Proportion No catch thru time by fisher/',NM,'_',series))
  ggsave(paste0(Out.name,'.tiff'), width = 12,height = 10, dpi = 300, compression = "lzw")
  
  print(d1%>%
          group_by(Ves.var,yr)%>%
          summarise(Target=mean(cpue))%>%
          ggplot(aes(yr,Target))+ 
          geom_point() + 
          geom_line() +facet_wrap(~Ves.var)+ ggtitle(TitlE)+
          ylab('Mean kg/km gn h') )
  ggsave(paste0(Out.name,'_cpue only.tiff'), width = 8,height = 8, dpi = 300, compression = "lzw")
  
  
  print(d1.bin%>%
          group_by(Ves.var,yr)%>%
          summarise(Target=mean(Target))%>%
          ggplot(aes(yr,Target))+ 
          geom_point() + 
          geom_line() +
          facet_wrap(~Ves.var)+ ggtitle(TitlE)+
          ylab('Proportion positive catch')) 
  ggsave(paste0(Out.name,'_prop only.tiff'), width = 8,height = 8, dpi = 300, compression = "lzw")
  
  print(d%>%
          mutate(yr=floor(time.var))%>%
          group_by(yr)%>%
          summarise(Target=mean(cpue),
                    SD=sd(cpue))%>%
          ggplot(aes(yr,Target))+
          geom_line(data=d1%>%
                      group_by(Ves.var,yr)%>%
                      summarise(Target=mean(cpue)),
                    aes(yr,Target,group = Ves.var),color='grey50')+ 
          geom_point(color='darkred',size=5) + 
          geom_line(color='darkred') +
          geom_errorbar(aes(ymin=Target-SD, ymax=Target+SD), width=.2,
                        position=position_dodge(.9),color='darkred')+
          ylab('Mean kg/km gn h (+/- SD)') + ggtitle('All vessels'))
  ggsave(paste0(Out.name,'_all vessels.tiff'), width = 8,height = 8, dpi = 300, compression = "lzw")
  
  
}

# IDENTIFY FISHING ON DIFFERENT HABITATS (~TARGETING BEHAVIOUR, only applicable to Daily logbooks) ----------------------------------------------

  #Stephens & McCall
library(forcats)
Species.catch.ranking=function(d,TITL,min.avrg.catch=50,minyears=5,minlocs=1,target,drop.noncollocated)
{
  #Create species matrix
  options(scipen=999)
  d_multi=d%>%
    dplyr::select(Same.return.SNo,zone,LIVEWT.c,SPECIES,SNAME)%>%
    dplyr::select(-c(SNAME))%>%spread(SPECIES,LIVEWT.c,fill=0)
  
  #Define co-located species and non-co-located species
  id.sp=match(sort(unique(d$SPECIES)),names(d_multi))
  bindata=d_multi[,id.sp]/d_multi[,id.sp]
  bindata[is.na(bindata)] = 0
  species_bincolnum = match(target, names(bindata))
  sp_colocation = colSums(bindata[, species_bincolnum] * bindata) ## check how many events contain catches of target + species
  out_table = data.frame(RSSpeciesId = as.numeric(names(sp_colocation)),
                         nlocs = unname(sp_colocation)) %>%
    arrange(desc(nlocs))
  noncollocated = out_table$RSSpeciesId[which(out_table$nlocs<minlocs)]

  SpeciesCatch = d %>%
    group_by(SNAME, SPECIES) %>%
    summarise(Total = sum(LIVEWT.c, na.rm=TRUE)/1000,
              nYears = length(unique(FINYEAR)),
              .groups="drop") %>%
    arrange(desc(Total)) %>%
    mutate(Select = case_when(SPECIES%in%noncollocated ~ "noncollocated",
                              .default = "included"),
           AverageCatch = Total / nYears,
           Select = ifelse(AverageCatch<min.avrg.catch | nYears<minyears, "omit", Select)) %>%
    as.data.frame
  
  SpeciesAnnualCatch = d %>%
    mutate(FishingSeason=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(SNAME, SPECIES, FishingSeason) %>%
    summarise(Catch = sum(LIVEWT.c, na.rm=TRUE)/1000, .groups="drop") %>%
    left_join(SpeciesCatch,by = c('SNAME', 'SPECIES')) %>%
    arrange(desc(Total)) %>%
    as.data.frame

  cc <- scales::seq_gradient_pal("darkturquoise", "firebrick4", "Lab")(seq(0,1,length.out=length(unique(d$FINYEAR))))
  p=ggplot(SpeciesAnnualCatch, aes(x = fct_reorder(SNAME, Total), y = Catch, fill=factor(FishingSeason))) +
    geom_col() +
    coord_flip() + scale_fill_manual(name = "Financial year",values=cc)+
    labs(x = "Species", y = "Catch (tonnes)", title = TITL)+
    geom_hline(yintercept=min.avrg.catch,color='orange')+
    theme_PA()+
    theme(legend.position = 'right',
          axis.text.y = element_text(size=8))
  x_cols = as.matrix(get_guide_data(p,"y"))[,2] 
  x_cols=data.frame(SNAME=x_cols)%>%
          left_join(SpeciesCatch%>%dplyr::select(SNAME,Select),by="SNAME")%>%
          mutate(Axis.col=ifelse(Select=='included','chartreuse4','grey50'))
  
  p=p+theme(axis.text.y = element_text(colour= x_cols$Axis.col))
  
  if(drop.noncollocated)
  {
    SpeciesCatch=SpeciesCatch%>%filter(!SPECIES%in%noncollocated)
    d_multi=d_multi[,-match(noncollocated,names(d_multi))]
    bindata=bindata[,-match(noncollocated,names(bindata))] 
    noncollocated=NA
  }
  
  xx = SpeciesCatch %>%
          filter(Select=="included" & !SPECIES==target) 
  indspecies_ids=xx%>%pull(SPECIES)
  indspecies_names=xx%>%pull(SNAME)
  
  return(list(p=p,dat=SpeciesCatch,noncollocated=noncollocated,d_multi=d_multi,bindata=bindata,
              indspecies_ids=indspecies_ids, indspecies_names=indspecies_names))
  
}
fn.prop.by.shot=function(NMS,dd,tar)
{
  dd1=dd%>%
    dplyr::select(where(is.numeric))%>%
    mutate(Tot=rowSums(.))%>%
    filter(Tot>0)
  dd1_prop <- dd1 %>%
    mutate(across(everything()), . / Tot)%>%
    dplyr::select(-Tot) %>%
    gather(Species,Proportion) %>%
    group_by(Species) %>%
    summarise(Mean.prop=mean(Proportion),
              SD.prop=sd(Proportion))%>%
    ungroup()%>%
    mutate(Species=as.numeric(Species))%>%
    arrange(-Mean.prop)%>%
    mutate(CL=ifelse(Species==tar,'brown','grey60'))%>%
    left_join(NMS,by=c('Species'='SPECIES'))
  dd1_prop$SNAME=factor(dd1_prop$SNAME,levels=rev(dd1_prop$SNAME))
  p=dd1_prop%>%
    ggplot(aes(SNAME,Mean.prop))+
    geom_col(aes(fill = CL))+
    geom_errorbar(aes(ymin = Mean.prop, ymax = Mean.prop+SD.prop,color= CL))+
    coord_flip() +ylim(0,NA)+
    theme_PA()+theme(legend.position = 'none')+
    ylab('Shot mean proportion')+xlab('')
  return(p)
}
FitStephensMacCallModel = function(smdat, species_ids, indspecies_ids,indspecies_names, use_model=MODL, refit_sig=FALSE, sigval=0.05)
{
  smdat$zone=factor(smdat$zone)
  zones = levels(smdat$zone)
  nzones = length(zones)
  indspecies_ids=paste0("ind", indspecies_ids)
  ## Fit glm (binomial) model to indicator variables
  if (use_model==1)
  {
    print("Fit Model ...")
    formula1_text = paste0("ind", species_ids, " ~ ", paste0(indspecies_ids, " ", collapse="+"))
    formula1 = as.formula(formula1_text)
    mod1 = glm(formula1, family = binomial(link = "logit"), data=smdat)
  }
  if (use_model==2)
  {
    print("Fit Model (including zone) ...")
    formula1_text = paste0("ind", species_ids, " ~ ", paste0(indspecies_ids, " ", collapse="+"), "+ zone + ", 
                           paste0(indspecies_ids, ":zone ", collapse="+"))
    formula1 = as.formula(formula1_text)
    mod1 = glm(formula1, family = binomial(link = "logit"),  data=smdat)
    
    # Refit model with only significant variables
    if (refit_sig)
    {
      # If your goal is causal inference, keep key control variables even if they're not significant.
      significant_vars <- names(which(summary(mod1)$coefficients[, 4] < 0.05))
      significant_vars <- setdiff(significant_vars, "(Intercept)")
      significant_spp = as.numeric(unique(str_remove(significant_vars, 'ind')))
      significant_spp = significant_spp[!is.na(significant_spp)]
      retained_ids = intersect(str_remove(indspecies_ids, 'ind'), significant_spp)
      
      formula1_text = paste0("ind", species_ids, " ~ ", paste0("ind", retained_ids, " ", collapse="+"), "+ zone + ",
                             paste0("ind", retained_ids, ":zone ", collapse="+"))
      formula1 = as.formula(formula1_text)
      mod1 = glm(formula1, family = binomial(link = "logit"), data=smdat)
    }

  }
  # summary(mod1)$coefficients  # Look at p-values
  lookuptable = data.frame(species=str_remove(indspecies_ids, 'ind'), 
                           RSSpeciesCommonName=indspecies_names)%>%
    mutate(RSSpeciesCommonName=ifelse(species=='ind99','noncollocated',RSSpeciesCommonName))
  lookuptable$RSSpeciesId = as.numeric(lookuptable$species)

  ## extract model coefficients
  modcoefs = data.frame(summary(mod1)$coefficients)
  ## meed to calc coeffs for each species in each zone
  modcoefs$Name = rownames(modcoefs)
  
  ## extract intercept and zone coefs (and stderrors)
  modfit_intercept = modcoefs %>%
    filter(Name=="(Intercept)")
  if (use_model==2)
  {
    modfit_zone = modcoefs %>%
      filter(str_detect(Name, "^zone"))
    modfit_zone_ests = data.frame(t(modfit_zone$Estimate))
    names(modfit_zone_ests) = paste0("Const_", rownames(modfit_zone))
    modfit_zone_stderr = data.frame(t(modfit_zone$Std..Error))
    names(modfit_zone_stderr) = paste0("Const_", rownames(modfit_zone))
  }
  ## now species coefs and species:zone coefs (and stderrors)
  modfit_species = modcoefs %>%
    filter(Name!="(Intercept)") %>%
    filter(!str_detect(Name, "^zone")) %>%
    separate(Name, into = c("RSSpeciesId", "zone"), sep = ":", fill = "right", remove = FALSE) %>%
    mutate(zone = ifelse(is.na(zone), "Species", zone))
  rownames(modfit_species) = NULL

  ## now overall estimates
  modfit_ests = modfit_species %>%
    dplyr::select(RSSpeciesId, Estimate, zone) %>%
    pivot_wider(names_from=zone, values_from=Estimate) %>%
    as.data.frame
  head(modfit_ests)
  if (use_model==2)
  {
    for (izone in zones[-1])
    {
      dcol = grep(izone, names(modfit_ests))
      dcolzone = grep(izone, names(modfit_zone_ests))
      modfit_ests[, dcol] = modfit_ests$Species + modfit_ests[, dcol] + modfit_zone_ests[, dcolzone]
    }
    names(modfit_ests)[-1] = zones
  }
  ## and now overall stderrors
  modfit_stderrs = modfit_species %>%
    dplyr::select(RSSpeciesId, Std..Error, zone) %>%
    pivot_wider(names_from=zone, values_from=Std..Error) %>%
    as.data.frame
  if (use_model==2)
  {
    for (izone in zones[-1])
    {
      dcol = grep(izone, names(modfit_stderrs))
      dcolzone = grep(izone, names(modfit_zone_stderr))
      modfit_stderrs[, dcol] = (0 + 
                                  modfit_stderrs$Species^2 + modfit_stderrs[, dcol]^2 + 
                                  modfit_zone_stderr[, dcolzone]^2)^0.5
    }
    names(modfit_stderrs)[-1] = zones

  }
  # modfit_stderrs$Species = (modfit_intercept$Std..Error^2 + modfit_stderrs$Species^2)^0.5
  modfit_ests_rev = modfit_ests %>%
    pivot_longer(cols=-RSSpeciesId, names_to="zone", values_to="Estimate") %>%
    as.data.frame
  modfit_stderrs_rev = modfit_stderrs %>%
    pivot_longer(cols=-RSSpeciesId, names_to="zone", values_to="StdErr") %>%
    as.data.frame
  modfit = modfit_ests_rev %>%
    left_join(modfit_stderrs_rev, by = join_by(RSSpeciesId, zone)) %>%
    mutate(zvalue = Estimate / StdErr,
           pvalue = 2 * (1 - pnorm(abs(zvalue))),
           RSSpeciesId = as.numeric(substring(RSSpeciesId, 4, nchar(RSSpeciesId)))) %>%
    as.data.frame
  modfit = modfit %>%
    left_join(lookuptable, by = join_by(RSSpeciesId)) %>%
   # mutate(zone = factor(zone, levels=zones)) %>%
    arrange(zone, Estimate)%>%
    mutate(zone=ifelse(is.na(zone),"Zones combined",zone),
           RSSpeciesCommonName=ifelse(is.na(RSSpeciesCommonName),'noncollocated',RSSpeciesCommonName))
  order_names = unique(modfit$RSSpeciesCommonName)
  order_ids = unique(modfit$RSSpeciesId)
  modfit$RSSpeciesId = factor(modfit$RSSpeciesId, levels=order_ids)
  modfit$RSSpeciesCommonName = factor(modfit$RSSpeciesCommonName, levels=order_names)
  modfit$sig = ifelse(modfit$pvalue<sigval, "sig", "nonsig")
  
  
  ## now compute predictions and find critical value
  smdat$Pred = predict(mod1, type="response")
  minprob = floor(quantile(smdat$Pred, 0.005) * 10) / 10
  maxprob = ceiling(quantile(smdat$Pred, 0.995) * 10) / 10
  output = data.frame(critval = seq(minprob, maxprob, 0.01),
                      diff = NA)
  spidcolnum = which(names(smdat)==paste0("ind", species_ids))
  for (critval in output$critval){
    smdat$pred_presence = ifelse(smdat$Pred >= critval, 1, 0)
    output$diff[match(critval, output$critval)] = abs(sum(smdat[, spidcolnum]) - sum(smdat$pred_presence))
  }
  smdat$pred_presence = NULL
  critical_value = output$critval[which.min(output$diff)]
  critical_diff = output$diff[which.min(output$diff)]
  smdat$critval = critical_value
  smdat$TargetSM = ifelse(smdat$Pred >= critical_value, 1, 0)

  
  return(list(use_model = use_model,
              modfit = modfit, 
              summary = summary(mod1),
              critical = c(critical_value=critical_value, critical_diff=critical_diff),
              difftable = output,
              data = smdat))
}
PlotStephensMacCallModel = function(smfit)
{
  modfit = smfit$modfit
  ## Extract significant variables -----
  sigfit = modfit %>%
    filter(sig=="sig")
  id_levels = unique(sigfit$RSSpeciesId)
  name_levels = unique(sigfit$RSSpeciesCommonName)
  sigfit$RSSpeciesCommonName = factor(sigfit$RSSpeciesCommonName, levels=name_levels)
  sigfit$RSSpeciesId = factor(sigfit$RSSpeciesId, levels=id_levels)
  if (smfit$use_model==2) sigfit$zone = factor(sigfit$zone, levels=unique(modfit$zone))
  
  p_all_RSSpeciesId = ggplot(modfit, aes(Estimate, RSSpeciesId, fill=sig)) +
                      geom_col(width = 0.6) +
                      scale_fill_manual(values=c(4,2)) +
                      geom_vline(xintercept = 0, lty=2) +
                      facet_grid(.~zone, scale="free_x") +
                      labs(fill="Sig.")
  p_all_RSSpeciesCommonName = ggplot(modfit, aes(Estimate, RSSpeciesCommonName, fill=sig)) +
                              geom_col(width = 0.6) +
                              scale_fill_manual(values=c(4,2)) +
                              geom_vline(xintercept = 0, lty=2) +
                              facet_grid(.~zone, scale="free_x") +
                              labs(fill="Sig.")
  p_sig_RSSpeciesId = ggplot(sigfit, aes(Estimate, RSSpeciesId)) +
                      geom_col(fill = 2, width = 0.6) +
                      geom_vline(xintercept = 0, lty=2) +
                      facet_grid(.~zone)
  p_sig_RSSpeciesCommonName = ggplot(sigfit, aes(Estimate, RSSpeciesCommonName)) +
                              geom_col(fill = 2, width = 0.6) +
                              geom_vline(xintercept = 0, lty=2) +
                              facet_grid(.~zone)

  return(list(p_all_RSSpeciesId=p_all_RSSpeciesId,
              p_all_RSSpeciesCommonName=p_all_RSSpeciesCommonName,
              p_sig_RSSpeciesId=p_sig_RSSpeciesId,
              p_sig_RSSpeciesCommonName=p_sig_RSSpeciesCommonName))
}
PlotStephensMacCallCriticalValues = function(smfit, species_ids, species_names)
{
  smdat = smfit$data

  ymax1 = max(smfit$difftable["diff"], na.rm=TRUE)
  gg1 = ggplot(smfit$difftable, aes(critval, diff)) +
    geom_line() +
    geom_point(x=smfit$critical["critical_value"], 
               y=smfit$critical["critical_diff"], col=2, size=3) +
    geom_label(x=smfit$critical["critical_value"], 
               y=smfit$critical["critical_diff"] + 0.1*ymax1, 
               label=smfit$critical["critical_value"], col=2) +
    theme_bw() +
    labs(title = paste0("Critical value (", species_names,")"),
         x="Probability", y="Events: abs(obs-pred)",
         subtitle = paste0("The critical probability at which the difference between the\nobserved and expected number of events encountering the\nspecies is minimised."))
  

  gg2 = ggplot(smdat, aes(Pred, fill=factor(TargetSM))) +
    geom_histogram(breaks=seq(0, 1, 0.05), col=1, show.legend=FALSE) +
    theme_bw() +
    labs(title = paste0("Predicted probabilities for presence of ", species_names,""),
         subtitle = "Blue shading indicates selected records for calculating CPUE index.",
         x="Probability", y="Frequency")
  
  gg3 = ggplot(smdat, aes(FishingSeason, Proportion, group=interaction(FishingSeason, TargetSM), fill=factor(TargetSM))) +
    geom_boxplot(show.legend = FALSE) +
    theme_bw() +
    scale_x_continuous(breaks=seq(2006, 2024, 2)) +
    labs(title = paste0("Proportion of ", species_names),
         x="Fishing Season", y="Proportion",
         subtitle = paste0("Critical value of predicted probability (Stephens-MacCall) is ", smfit$critical["critical_value"] ,".\nBlue shading indicates selected records for calculating CPUE index."),
         fill="Target")
  

  gg4 = ggplot(smdat, aes(FishingSeason, CPUE, group=interaction(FishingSeason, TargetSM), fill=factor(TargetSM))) +
    geom_boxplot(show.legend = FALSE) +
    theme_bw() +
    scale_x_continuous(breaks=seq(2006, 2024, 2)) +
    labs(title = paste0("CPUE of ", species_names),
         x="Fishing Season", y="CPUE",
         subtitle = paste0("Critical value of predicted probability (Stephens-MacCall) is ", smfit$critical["critical_value"] ,".\nBlue shading indicates selected records for calculating CPUE index."),
         fill="Target")
  
  ylim1 = boxplot.stats(smdat$CPUE)$stats[c(1, 5)] # compute lower and upper whiskers
  gg5 = gg4 + coord_cartesian(ylim = ylim1*1.05)# scale y limits based on ylim1
  
  
  
  return(list(gg1=gg1,gg2=gg2,gg3=gg3,gg4=gg4,gg5=gg5))
}
PlotStephensMacCallTarget.vs.All_cpue=function(smfit, species_names)
{
  smdata = smfit$data
  
  nomcpueALL = smdata %>%
    group_by(FishingSeason) %>%
    summarise(n = length(CPUE),
              m = length(CPUE[CPUE>0]),
              mean.lognz = mean(log(CPUE[CPUE>0])),
              sd.lognz = sd(log(CPUE[CPUE>0]))) %>%
    mutate(p.nz = m/n,
           theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
           c = (1-p.nz)^(n-1),
           d = 1+(n-1)*p.nz,
           vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
             sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
           mean = exp(theta),
           lowCL = exp(theta - 1.96*sqrt(vartheta)),
           uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
    mutate(Group = "All") %>%
    as.data.frame
  
  nomcpueTARG = smdata %>%
    filter(TargetSM==1) %>%
    group_by(FishingSeason) %>%
    summarise(n = length(CPUE),
              m = length(CPUE[CPUE>0]),
              mean.lognz = mean(log(CPUE[CPUE>0])),
              sd.lognz = sd(log(CPUE[CPUE>0]))) %>%
    mutate(p.nz = m/n,
           theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
           c = (1-p.nz)^(n-1),
           d = 1+(n-1)*p.nz,
           vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
             sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
           mean = exp(theta),
           lowCL = exp(theta - 1.96*sqrt(vartheta)),
           uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
    mutate(Group = "Target") %>%
    as.data.frame
  
  nomcpue = suppressMessages(nomcpueALL %>%
                               full_join(nomcpueTARG))
  
  p=ggplot(nomcpue, aes(FishingSeason, mean, ymin=lowCL, ymax=uppCL, col=Group)) +
    geom_errorbar(width=0, linewidth=1, position=position_dodge(width=0.3)) +
    geom_line(lty=2, position=position_dodge(width=0.3)) +
    geom_point(shape=21, fill="white", position=position_dodge(width=0.3), size=2) +
    theme_bw() +
    theme(legend.position.inside=c(0.1, 0.85),
          legend.title = element_blank()) +
    scale_color_manual(values=c(8, 2)) +
    scale_x_continuous(breaks=seq(2008, 2024, 2), expand = c(0.01, 0.01)) +
    expand_limits(x=c(2008, 2024)) +
    labs(title = paste0("Nominal CPUE of ", species_names),
         x="Fishing Season", y="CPUE",
         subtitle = paste0("Deltalognormal distribution"))
  return(p)
}
fn.untangle=function(dat,axs,tar.sp)
{
  dat.high.cpue=dat%>%
    filter(CPUE>=quantile(dat$CPUE,probs=.90))%>%
    filter(Pred>=0.975 | Pred<=0.5)%>%
    select(starts_with("a_"),CPUE,Pred)%>%
    gather(Species,Catch,-c(CPUE,Pred))%>%
    mutate(Species=as.numeric(str_remove(Species,'a_')),
           Pred=paste('Prob=',round(Pred,1)))%>%
    filter(CPUE<quantile(dat$CPUE,probs=.95))
  dat.low.cpue=dat%>%
    filter(CPUE<=quantile(dat$CPUE,probs=0.15))%>%
    filter(Pred>=0.975 | Pred<=0.1)%>%
    select(starts_with("a_"),CPUE,Pred)%>%
    gather(Species,Catch,-c(CPUE,Pred))%>%
    mutate(Species=as.numeric(str_remove(Species,'a_')),
           Pred=paste('Prob=',round(Pred,1)))
  
  # a=dat%>%
  #   filter(CPUE>10 & Pred<0.5)%>%
  #   select(starts_with("a_"),CPUE,Pred)%>%
  #   gather(Species,Catch,-c(CPUE,Pred))%>%
  #   mutate(Species=as.numeric(str_remove(Species,'a_')),
  #          Pred=paste('Prob=',round(Pred,1)))
  
  Lvls=sort(unique(dat.high.cpue$Species))
  names(Lvls)=rep('black',length(Lvls))
  names(Lvls)[match(tar.sp,Lvls)]='brown'
  p1=dat.high.cpue%>%
    mutate(Species=factor(Species,levels=Lvls))%>%
    ggplot(aes(Species,Catch,color=CPUE))+
    geom_point()+
    facet_wrap(~Pred,nrow=1)+coord_flip()+ylab('Catch of co-occurring species')+xlab('')+
    theme_PA(lgT.siz=8,leg.siz=6,strx.siz=8,axs.t.siz=6,axs.T.siz=10)+
    theme(legend.position = 'top',axis.text.y = element_text(size = axs))+
    ylim(0,quantile(dat.high.cpue$Catch,probs=0.99))+
    theme(axis.text.y = element_text(colour = names(Lvls)))
  
  p2=NULL
  if(nrow(dat.low.cpue)>0)
  {
    Lvls=sort(unique(dat.low.cpue$Species))
    #Lvls=sort(unique(a$Species))
    names(Lvls)=rep('black',length(Lvls))
    names(Lvls)[match(tar.sp,Lvls)]='brown'
    p2=dat.low.cpue%>%
      mutate(Species=factor(Species,levels=Lvls))%>%
      ggplot(aes(Species,Catch,color=CPUE))+
      geom_point()+
      facet_wrap(~Pred,nrow=1)+coord_flip()+ylab('Catch of co-occurring species')+xlab('')+
      theme_PA(lgT.siz=8,leg.siz=6,strx.siz=8,axs.t.siz=6,axs.T.siz=10)+
      theme(legend.position = 'top',axis.text.y = element_text(size = axs))+
      ylim(0,quantile(dat.low.cpue$Catch,probs=0.99))+
      theme(axis.text.y = element_text(colour = names(Lvls)))
  }
 
  return(list(p1=p1,p2=p2))
}

  #cluster
fn.cluster.Stephen.MacCall=function(a,selected.species,n.clus,target,check.clust.num,out.clara,
                                    apply.scale,do.proportion)
{
  a=a[,match(selected.species,names(a))]
  
  if(apply.scale) a=scale(a)
  if(do.proportion)  a=a/rowSums(a)
  

  #step 2. Determine optimum number of clusters
  if(check.clust.num)
  {
    ran.samp=sample(1:nrow(a),round(nrow(a)*.5),replace=F) #random sample to reduce computation time
    b=fviz_nbclust(a[ran.samp,], clara, method = "silhouette",print.summary=T,k.max=6,nboot=10)
    print(b+theme_classic())
    ggsave(paste(HndL.Species_targeting,"Cluster/CLARA_optimal_numbers_",target,".tiff",sep=""),
           width = 8,height = 8, dpi = 300, compression = "lzw")
    num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
  }
  
  #step 3. fit clara
  if(!exists("num.clus")) num.clus=n.clus
  clara.res <- clara(a, num.clus, samples = 50, pamLike = TRUE)
  
  #step 4. visualize CLARA clusters in data scattergram
  if(out.clara)
  {
    print(fviz_cluster(clara.res, 
                 palette = rainbow(num.clus), # color palette
                 ellipse.type = "t", # Concentration ellipse
                 geom = "point", pointsize = 1,
                 ggtheme = theme_classic()))
    ggsave(paste(HndL.Species_targeting,"Cluster/CLARA_cluster_",target,".tiff",sep=""),
           width = 8,height = 8, dpi = 300, compression = "lzw") 
  }
  
  #step 5. export cluster to add to input data
  dd.clara <- cbind(as.data.frame(a), cluster_clara = clara.res$cluster)
  dd.clara=dd.clara%>%
    rownames_to_column(var = "Same.return.SNo")%>%
    dplyr::select(cluster_clara,Same.return.SNo)%>%remove_rownames()
  return(dd.clara)
}
fn.cluster=function(data,TarSp,n.clus,target,check.clustrbl)
{
  a=data[[TarSp]]%>%
    filter(Same.return.SNo%in%unique(DATA.list.LIVEWT.c.daily[[TarSp]]$Same.return.SNo))%>%
    filter(!SPECIES%in%31000)%>%
    mutate(SPECIES=ifelse(SPECIES==19004,19000,SPECIES),
           Same.return.SNo.block10=paste(Same.return.SNo,block10))
  TOP.ktch=a%>%group_by(SPECIES)%>%
    summarise(Catch=sum(LIVEWT.c))%>%
    arrange(-Catch)%>%
    mutate(Cum.ktch=cumsum(Catch),
           Quantiles=Cum.ktch/sum(Catch))
  top.sp=TOP.ktch$SPECIES[1:which.min(abs(TOP.ktch$Quantiles-95/100))]
  
  a=a%>%filter(SPECIES%in%top.sp)%>%
    group_by(Same.return.SNo.block10,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T),
              Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c))%>%
    mutate(cpue=LIVEWT.c/Km.Gillnet.Hours.c)%>%
    dplyr::select(SPECIES,cpue,Same.return.SNo.block10)%>% 
    spread(SPECIES,cpue)%>%
    column_to_rownames(var = "Same.return.SNo.block10")
  
  #proportion
  a[is.na(a)]=0
  a=a/rowSums(a,na.rm = T)
  
  
  #step 1. Define if data are clusterable
  if(check.clustrbl=="YES")
  {
    #random sample to reduce computation time
    ran.samp=sample(1:nrow(a),5000,replace=F)
    
    res <- get_clust_tendency(a[ran.samp,], n = nrow(a[ran.samp,])-1, graph = FALSE)
    if(1-res$hopkins_stat>0.75) clusterable="YES"else  clusterable="NO"
    print(clusterable)
  }
  
  #step 2. Determine optimum number of clusters
  if(check.clustrbl=="YES")
  {
    fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_optimal_numbers_",target,sep=""),2400,2400)
    b=fviz_nbclust(a, clara, method = "silhouette",print.summary=T)
    b+theme_classic()
    dev.off()
    num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
  }
  
  #step 3. fit clara
  if(!exists("num.clus")) num.clus=n.clus
  clara.res <- clara(a, num.clus, samples = 50, pamLike = TRUE)
  
  #step 4. visualize CLARA clusters in data scattergram
  if(Model.run=="First")
  {
    fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster_",target,sep=""),2400,2400)
    fviz_cluster(clara.res, 
                 palette = rainbow(num.clus), # color palette
                 ellipse.type = "t", # Concentration ellipse
                 geom = "point", pointsize = 1,
                 ggtheme = theme_classic())
    dev.off()
  }
  
  #step 5. export cluster to add to input data
  dd.clara <- cbind(as.data.frame(a), cluster_clara = clara.res$cluster)
  dd.clara=dd.clara%>%rownames_to_column(var = "Same.return.SNo.block10")%>%
    dplyr::select(cluster_clara,Same.return.SNo.block10)%>%remove_rownames()
  return(dd.clara)
}
fn.cluster.catch=function(data,TarSp,target,varS,scaling,check.clustrbl,n.clus)
{
  a=data[[TarSp]]%>%
    mutate(Same.return.SNo.block10=paste(Same.return.SNo,block10))%>%
    column_to_rownames(var = "Same.return.SNo.block10")%>%
    dplyr::select(varS[-match(target,varS)])
  if(scaling=="YES")a=scale(a)
  
  #step 1. Define if data are clusterable
  if(check.clustrbl=="YES")
  {
    #random sample to reduce computation time
    ran.samp=sample(1:nrow(a),15000,replace=F)
    
    res <- get_clust_tendency(a[ran.samp,], n = nrow(a[ran.samp,])-1, graph = FALSE)
    if(1-res$hopkins_stat>0.75) clusterable="YES"else  clusterable="NO"
    print(clusterable)
  }
  
  #step 2. Determine optimum number of clusters
  if(check.clustrbl=="YES")
  {
    fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_optimal_numbers_",target,sep=""),2400,2400)
    b=fviz_nbclust(a, clara, method = "silhouette",print.summary=T)
    b+theme_classic()
    dev.off()
    num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
  }
  
  #step 3. fit clara
  if(!exists("num.clus")) num.clus=n.clus
  clara.res <- clara(a, num.clus, samples = 50, pamLike = TRUE)
  
  #step 4. visualize CLARA clusters in data scattergram
  fn.fig(paste(HndL.Species_targeting,"Cluster/CLARA_cluster_",target,sep=""),2400,2400)
  fviz_cluster(clara.res, 
               palette = rainbow(num.clus), # color palette
               ellipse.type = "t", # Concentration ellipse
               geom = "point", pointsize = 1,
               ggtheme = theme_classic())
  dev.off()
  
  #step 5. add cluster to input data
  dd.clara <- cbind(as.data.frame(a), cluster_clara = clara.res$cluster)
  dd.clara=dd.clara%>%rownames_to_column(var = "Same.return.SNo")%>%
    select(cluster_clara,Same.return.SNo)%>%remove_rownames()
  return(dd.clara)
}
fn.compare.targeting=function(DAT,Drop.var,Title)
{
  d=DAT%>% dplyr::select(c(cluster_clara,Clus.vars[-match(Drop.var,Clus.vars)]))
  d[-1]=d[-1]/rowSums(d[-1])
  colnames(d)[-1]=gsub("Catch.", "", colnames(d)[-1], fixed = TRUE)
  p=d%>%
    gather('ID','value',-cluster_clara)%>%           
    ggplot(aes(x=ID, y=value, fill=factor(cluster_clara))) +
    geom_boxplot()+
    theme(legend.title = element_blank())+
    labs(title = Title,x = "",y="Proportion of catch")
  print(p)
}


#Simulated multivariate data 
do.sim.multivar=FALSE
if(do.sim.multivar)
{
  n.samps=1e3
  dumy.dat=data.frame(Species=rep(c("Gummy","Dusky","Whiskery","Sandbar","Snapper","Blue morwong"),each=n.samps))%>%
    mutate(habitat=case_when(Species%in%c("Gummy","Blue morwong","Snapper")~sample(c('sand','reef','deep'),n(),replace=T,prob=c(.7,.2,.1)),
                             Species%in%c("Dusky","Whiskery")~sample(c('sand','reef','deep'),n(),replace=T,prob=c(.2,.7,.1)),
                             Species%in%c("Sandbar")~sample(c('sand','reef','deep'),n(),replace=T,prob=c(.1,.2,.7))),
           Catch=case_when(Species%in%c("Gummy","Blue morwong","Snapper") & habitat=='sand'~runif(n(),min=100,max=500),
                           Species%in%c("Gummy","Blue morwong","Snapper") & habitat=='reef'~round(runif(n(),min=0,max=10)),
                           Species%in%c("Gummy","Blue morwong","Snapper") & habitat=='deep'~round(runif(n(),min=0,max=1)),
                           Species%in%c("Dusky","Whiskery") & habitat=='sand'~round(runif(n(),min=0,max=10)),
                           Species%in%c("Dusky","Whiskery") & habitat=='reef'~runif(n(),min=100,max=500),
                           Species%in%c("Dusky","Whiskery") & habitat=='deep'~round(runif(n(),min=0,max=1)),
                           Species%in%c("Sandbar") & habitat=='sand'~round(runif(n(),min=0,max=1)),
                           Species%in%c("Sandbar") & habitat=='reef'~round(runif(n(),min=0,max=10)),
                           Species%in%c("Sandbar") & habitat=='deep'~runif(n(),min=100,max=500)))%>%
    group_by(habitat,Species)%>%
    mutate(Count = 1:n())%>%
    ungroup()%>%
    mutate(Sample.number=paste(habitat,Count))
}

#Test performance of Stephen & MacCall on simulated data
test.sim.StephMac=FALSE
if(test.sim.StephMac)
{
  output_dir=paste0(HndL.Species_targeting,'Stephens_McCall/z_test on simulated data/')
  minlocs.vec=1
  MODL=1
  dumy.dat=dumy.dat%>%
    rename(SNAME=Species)%>%
    mutate(Species=case_when(SNAME=='Gummy'  ~  17001,
                             SNAME=='Whiskery'  ~  17003,
                             SNAME=='Dusky' ~   18003,
                             SNAME=='Sandbar'  ~  18007,
                             SNAME=='Snapper' ~ 353001,
                             SNAME=='Blue morwong' ~ 377004))
  
  p1=dumy.dat%>%
    ggplot(aes(SNAME,Catch,fill=habitat))+
    geom_boxplot()+coord_flip()+theme(legend.position = 'top')
  UNIK=dumy.dat%>%distinct(Species,SNAME)
  UNIK.names=UNIK$SNAME
  UNIK=UNIK$Species
  
  for(i in 1:length(UNIK))
  {
    target=UNIK[i]
    names(target)=UNIK.names[i]
    
    indspecies_ids=dumy.dat%>%distinct(SNAME,Species)
    indspecies_names=indspecies_ids$SNAME
    indspecies_ids=indspecies_ids$Species
    indspecies_ids=subset(indspecies_ids,!indspecies_ids==target)
    indspecies_names=subset(indspecies_names,!indspecies_names==names(target))
    
    dumy.dat.spread=dumy.dat%>%
      dplyr::select(-SNAME)%>%
      spread(Species,Catch,fill=0)%>%
      data.frame
    colnames(dumy.dat.spread)=str_remove(colnames(dumy.dat.spread),'X')
    id.sp=match(c(target,indspecies_ids),names(dumy.dat.spread))
    bindata=dumy.dat.spread[,id.sp]
    bindata[bindata>0]=1
    names(bindata)=paste0('ind',names(bindata))
    
    dumy.dat.spread$Catch=dumy.dat.spread[,match(target,names(dumy.dat.spread))] 
    dumy.dat.spread$Proportion=dumy.dat.spread$Catch/rowSums(dumy.dat.spread[,id.sp])
    dumy.dat.spread=dumy.dat.spread%>%
      mutate(Effort=1,
             zone=1,
             FishingSeason=1,
             CPUE=Catch/Effort)
    
    dat.comb=cbind(dumy.dat.spread,bindata)
    
    
    smfit=FitStephensMacCallModel(smdat=dat.comb, species_ids=target, indspecies_ids,indspecies_names,use_model=MODL)
    
    #plot model coefs
    Plots=PlotStephensMacCallModel(smfit)
    names(Plots)=c('coefs.All_id','coefs.All_name','coefs.Sig_id','coefs.Sig_name')
    

    #plot cpue vs pred prob
    p3=smfit$data%>%
      mutate(TargetSM=as.character(TargetSM))%>%
      ggplot(aes(CPUE,Pred,color = TargetSM))+
      geom_point()+facet_wrap(~habitat)+
      coord_cartesian(xlim = c(0,quantile(smfit$data$CPUE,probs=0.95)))+
      theme_PA()
    
    ggarrange(plotlist = list(p1+ggtitle(paste('Simulated data. Target=',names(target))),Plots[[2]],p3),ncol=1,nrow=3)
    plot_name=paste0(names(target),'.tiff')
    ggsave(paste0(output_dir, plot_name),width = 5, height = 10,compression="lzw")
    
  }
  
  
}

#Test performance of cluster analysis on simulated data
test.sim.cluster=FALSE
if(test.sim.cluster)
{
  p1=dumy.dat%>%
    ggplot(aes(Species,Catch,fill=habitat))+
    geom_boxplot()+coord_flip()+theme(legend.position = 'top')
  
  fn.sim.cluster=function(a,selected.species,n.clus,check.clustrbl,apply.scale,do.bin,do.proportion)
  {
    a.ori=a
    row.names(a)=a$Sample.number
    a=a[,match(selected.species,names(a))]
    if(apply.scale) a=scale(a)
    if(do.bin)  a[a>0]=1
    if(do.proportion)  a=a/rowSums(a)
    
    #step 2. Determine optimum number of clusters
    if(check.clustrbl)
    {
      ran.samp=sample(1:nrow(a),round(nrow(a)*.5),replace=F) #random sample to reduce computation time
      b=fviz_nbclust(a[ran.samp,], clara, method = "silhouette",print.summary=T,k.max=6,nboot=10)+
        theme_classic()
      num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
    }
    
    #step 3. fit clara
    if(!exists("num.clus")) num.clus=n.clus
    clara.res <- clara(a, num.clus, samples = 50, pamLike = TRUE)
    
    #step 4. visualize CLARA clusters in data scattergram
    p=fviz_cluster(clara.res, 
                   palette = rainbow(num.clus), # color palette
                   ellipse.type = "t", # Concentration ellipse
                   geom = "point", pointsize = 1,
                   ggtheme = theme_classic())
    
    
    #step 5. export cluster to add to input data
    dd.clara <- cbind(as.data.frame(a), cluster_clara = clara.res$cluster)
    dd.clara=dd.clara%>%
      rownames_to_column(var = "Sample.number")%>%
      dplyr::select(cluster_clara,Sample.number)%>%remove_rownames()
    
    a.ori=a.ori%>%left_join(dd.clara,by='Sample.number')
    
    
    p2=a.ori%>%
      dplyr::select(-c(Count,Sample.number))%>%
      gather(Species,Catch,-c(habitat,cluster_clara))%>%
      mutate(cluster_clara=paste('cluster',cluster_clara))%>%
      ggplot(aes(Species,Catch,fill=cluster_clara))+
      geom_boxplot()+coord_flip()+theme(legend.position = 'top')
    
    return(list(b=b,p=p,p2=p2))
    
  }
  
  dumy.dat.spread=dumy.dat%>%
    spread(Species,Catch,fill=0)%>%data.frame
  
  dd_scale=fn.sim.cluster(a=dumy.dat.spread,
                          selected.species=c("Gummy","Dusky","Whiskery","Sandbar","Snapper","Blue.morwong"),
                          n.clus=NA,
                          check.clustrbl=TRUE,
                          apply.scale=TRUE,
                          do.bin=FALSE,
                          do.proportion=FALSE)
  
  dd_proportion=fn.sim.cluster(a=dumy.dat.spread,
                               selected.species=c("Gummy","Dusky","Whiskery","Sandbar","Snapper","Blue.morwong"),
                               n.clus=NA,
                               check.clustrbl=TRUE,
                               apply.scale=FALSE,
                               do.bin=FALSE,
                               do.proportion=TRUE)
  
  dd_bin=fn.sim.cluster(a=dumy.dat.spread,
                        selected.species=c("Gummy","Dusky","Whiskery","Sandbar","Snapper","Blue.morwong"),
                        n.clus=NA,
                        check.clustrbl=TRUE,
                        apply.scale=FALSE,
                        do.bin=TRUE,
                        do.proportion=FALSE)
  
  dd_bin.scaled=fn.sim.cluster(a=dumy.dat.spread,
                               selected.species=c("Gummy","Dusky","Whiskery","Sandbar","Snapper","Blue.morwong"),
                               n.clus=NA,
                               check.clustrbl=TRUE,
                               apply.scale=TRUE,
                               do.bin=TRUE,
                               do.proportion=FALSE)
  
  ggarrange(plotlist = list(dd_scale$b,dd_scale$p,p1,dd_scale$p2),ncol=2,nrow=2)
  ggarrange(plotlist = list(dd_proportion$b,dd_proportion$p,p1,dd_proportion$p2),ncol=2,nrow=2)
  ggarrange(plotlist = list(dd_bin$b,dd_bin$p,p1,dd_bin$p2),ncol=2,nrow=2)
  ggarrange(plotlist = list(dd_bin.scaled$b,dd_bin.scaled$p,p1,dd_bin.scaled$p2),ncol=2,nrow=2)
  
}


# COMPUTE FOLLY AND NOMINAL INDICES FOR EXPORTING ----------------------------------------------
fn.out.effective=function(a)
{
  colnames(a)=tolower(colnames(a))
  out = a %>%
    group_by(finyear) %>%
    summarise(Sumy = sum(livewt.c),
              Sumx = sum(km.gillnet.hours.c),
              My = mean(livewt.c),
              Mx = mean(km.gillnet.hours.c),
              Sy = sd(livewt.c),
              Sx = sd(km.gillnet.hours.c),
              r = cor(livewt.c, km.gillnet.hours.c),
              n = length(livewt.c)) %>%
    as.data.frame
  out$r[is.na(out$r)] = 0
  out = out %>%
    mutate(mean=Sumy/Sumx,
           se =  sqrt(1/n*(My^2*Sx^2/(Mx^4) + Sy^2/(Mx^2) - 2*My*r*Sx*Sy/(Mx^3))),
           lowCL = mean - 1.96*se,
           uppCL = mean + 1.96*se) %>%
    as.data.frame
  return(out)
}
fn.out.nominal=function(d,method)
{
  if(is.na(match('cpue.target',names(d))))d=d%>%mutate(cpue.target=catch.target/km.gillnet.hours.c)
  
  # Nominal CPUE
  if(method == "Nominal")
  {
    out = d %>%
      group_by(finyear) %>%
      summarise(My = mean(catch.target),
                Mx = mean(effort),
                Sy = sd(catch.target),
                Sx = sd(effort),
                r = cor(catch.target, effort),
                n = length(catch.target)) %>%
      as.data.frame
    out$r[is.na(out$r)] = 0
    out = out %>%
      mutate(mean=My/Mx,
             se =  sqrt(1/n*(My^2*Sx^2/(Mx^4) + Sy^2/(Mx^2) - 2*My*r*Sx*Sy/(Mx^3))),
             lowCL = mean - 1.96*se,
             uppCL = mean + 1.96*se) %>%
      as.data.frame
  }
  
  # Arithmetic Mean CPUE
  if(method == "Mean")
  {
    out = d %>%
      group_by(finyear) %>%
      summarise(mean = mean(cpue.target),
                n = length(cpue.target),
                sd = sd(cpue.target)) %>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      as.data.frame
  }
  
  # Lognormal CPUE
  if(method == "LnMean")
  {
    out = d %>%
      filter(cpue.target>0) %>%
      group_by(finyear) %>%
      summarise(ymean = mean(LNcpue),
                n = length(LNcpue),
                ysigma = sd(LNcpue)) %>%
      mutate(mean=exp(ymean + ysigma^2/2) ,
             ySE = sqrt(ysigma^2/n+ysigma^4/(2*(n-1))),
             lowCL = exp(ymean + ysigma^2/2 - 1.96*ySE),
             uppCL = exp(ymean + ysigma^2/2 + 1.96*ySE)) %>%
      as.data.frame
    
  }
  
  # Delta - Lognormal CPUE
  if(method == "DLnMean")
  {
    out = d %>%
      group_by(finyear) %>%
      summarise(n = length(cpue.target),
                m = length(cpue.target[cpue.target>0]),
                mean.lognz = mean(log(cpue.target[cpue.target>0])),
                sd.lognz = sd(log(cpue.target[cpue.target>0]))) %>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
      as.data.frame
  }
  
  return(out)
}
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

fn.ainslie=function(dat,Ktch.targt,Effrt,explr,QL_prop_ktch,Prop.Malcolm,cpue.units,spname,BLks,VesL,Type)
{
  names(dat) =  casefold(names(dat))
  
  ## some variable names
  dat$catch = dat[,match(Ktch.targt,names(dat))]
  dat$year = dat$year.c
  dat$fymonth = factor(dat$month, levels=c(7:12, 1:6))
  dat$season = as.numeric(substring(dat$finyear, 1, 4))
  dat$smonth = factor(dat$month, levels=c(7:12, 1:6))
  
  CPUE.All=vector('list',length(Effrt))
  names(CPUE.All)=Effrt
  CPUE.blk_vess=CPUE.non_zero=CPUE.QL_target=CPUE.Malcolm=CPUE.All
  
  for(ef in 1:length(Effrt))
  {
    pdf(paste(Hnd.ains,spname,Type,paste(Effrt[ef]),".pdf",sep=""))
    
    dat$effort = dat[,match(Effrt[ef],names(dat))]
    
    if(explr=="YES")
    {
      table(dat$year, dat$month)
      table(dat$season, dat$smonth)
      
      all.years = sort(unique(dat$year))
      par(mfrow=c(2,2))
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.gummy, list(year, month), sum)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.whiskery, list(year, month), sum)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.dusky, list(year, month), sum)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.sandbar, list(year, month), sum)), xlab="", ylab="", las=1, main="Sandbar")
      
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.gummy/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.whiskery/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.dusky/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.years, y=1:12, z=with(dat, tapply(catch.sandbar/km.gillnet.days.c, list(year, month), mean)), xlab="", ylab="", las=1, main="Sandbar")
      
      
      all.seasons = sort(unique(dat$season))
      par(mfrow=c(2,2))
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.gummy, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.whiskery, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.dusky, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.sandbar, list(season, smonth), sum)), xlab="", ylab="", las=1, main="Sandbar")
      
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.gummy/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Gummy")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.whiskery/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Whiskery")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.dusky/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Dusky")
      image(x=all.seasons, y=1:12, z=with(dat, tapply(catch.sandbar/km.gillnet.days.c, list(season, smonth), mean)), xlab="", ylab="", las=1, main="Sandbar")
      
      
      par(mfrow=c(2,1))
      boxplot(catch.target/km.gillnet.days.c ~ year, dat)
      boxplot(catch.target/km.gillnet.days.c ~ season, dat)
      
      boxplot(catch.target ~ month, dat, main="")
      boxplot(catch.target/km.gillnet.days.c ~ month, dat)
      
    }
    
    
    ## calculate year-specific QL and add column to data set
    
    dat$prop = dat$catch.target / dat$catch.total
    
    if(ef==1)
    {
      par(mfrow=c(1,1))
      boxplot(prop ~ season, dat)
      mtext("Proportion of target species catch out of total catch",3,-1,col="red")
    } 
    qldat = CalcQL(dat, prop.catch=QL_prop_ktch)
    dat = merge(dat, qldat, all=TRUE)
    dat$target = ifelse(dat$prop > dat$ql, 1, 0)
    if(ef==1)
    {
      par(mfrow=c(2,1))
      boxplot(catch.target/km.gillnet.days.c ~ season, dat)
      mtext("Qualification levels_all",3,-1,col="red")
      boxplot(catch.target/km.gillnet.days.c ~ season, subset(dat, target==1))
      mtext("Qualification levels_target only",3,-1,col="red")
    }
    
    
    ## malcolm's targeting - proportion below which no variation exists in cpue
    
    dat$cpue.target = dat$catch / dat$effort
    all.seasons = unique(dat$season)
    if(ef==1)
    {
      smart.par(n.plots=length(all.seasons+1),MAR=c(2,2,.1,.1),OMA=c(2,2,.5,.5),MGP=c(1,.5,0))
      with(dat, plot(prop, cpue.target,  pch=16, col=rgb(1,0,0,0.1), ylab='',xlab=''))
      legend('top',"All Seasons",bty='n')
      ## hmmm this is hard to see - probably very low <0.1
      for (i in all.seasons)
      {
        with(subset(dat, season==i), plot(prop, cpue.target, pch=16, ylab='',xlab='', col=rgb(0,0,1,0.1)))
        legend('top',paste(i),bty='n')
      }
      ## still say <0.1
      mtext(paste("Cpue (",cpue.units[ef],")",sep=""),2,outer=T,las=3)
      mtext("Proportion",1,outer=T)
    }
    
    ## plot raw mean cpues and CIs using 4 different data sets
    
    #compare different subsets of the data
    par(mfrow=c(3,2),mar=c(3,3,1,1),oma=c(1,1,1,1),mgp=c(2,.5,0))
    
    CPUE.All[[ef]] = CalcMeanCPUE(cpuedata = dat, catch.column="catch", effort.column="effort",
                                  plot.title = paste(spname, "_All rec"), cpue.units = cpue.units[ef], 
                                  draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
    
    CPUE.blk_vess[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, blockx%in%BLks & vessel%in%VesL), catch.column="catch", 
                                       effort.column="effort",plot.title = paste(spname, "_indicative_blk_ves"),
                                       cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="NO")
    
    CPUE.non_zero[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, catch>0), catch.column="catch", 
                                       effort.column="effort",plot.title = paste(spname, "_Nonzero"),
                                       cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
    
    CPUE.QL_target[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, target==1), catch.column="catch",
                                        effort.column="effort",plot.title = paste(spname, "_QL Target"),
                                        cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=F,PaR="NO",showLNMean="YES")
    
    if(nrow(subset(dat, prop>Prop.Malcolm))>100)CPUE.Malcolm[[ef]] = CalcMeanCPUE(cpuedata = subset(dat, prop>Prop.Malcolm), catch.column="catch",
                                                                                  effort.column="effort",plot.title = paste(spname, "_Malcolm_Prop>",Prop.Malcolm),
                                                                                  cpue.units = cpue.units[ef],draw.plot=TRUE, show.legend=TRUE,PaR="NO",showLNMean="YES")
    dev.off()
  }
  
  #show proportion of records selected by year
  #note: show that prob of catching doesn't change thru time
  TAB=subset(dat,select=c(finyear,target)) %>%
    group_by(finyear,target) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n))%>%
    as.data.frame
  TAB$yr=substr(TAB$finyear,1,4)
  
  return(list(CPUE.All=CPUE.All,CPUE.blk_vess=CPUE.blk_vess,CPUE.non_zero=CPUE.non_zero,
              CPUE.QL_target=CPUE.QL_target,CPUE.Malcolm=CPUE.Malcolm,
              QL_dat=subset(dat, target==1),Prop.selected=TAB))
}

# EVALUATE BALANCE OF DATA FOR QL ----------------------------------------------
fn.check.balanced=function(d,SP,what,MN.YR,pLot)
{
  fn.plt=function(dd)
  {
    a=dd
    a[a>0]=1
    Orderd=rev(sort(rowSums(a)))
    dd=as.matrix(dd[match(names(Orderd),row.names(dd)),])
    Nx=c(1,ncol(dd))
    Ny=c(1,nrow(dd)+2)
    Mx=max(dd)
    par(mar=c(2,3,.8,.8),mgp=c(2,.5,0))
    plot(1,1,xlim=Nx,ylim=Ny,col="transparent",ann=F,xaxt='n',yaxt='n')
    Selected=names(Orderd[Orderd>=MN.YR])
    show.pol=match(Selected,row.names(dd))
    Nx.p=c(Nx[1]-1,Nx[2]+1)
    if(length(Selected)>0)polygon(x=c(Nx.p,rev(Nx.p)),y=c(rep(show.pol[1]-1,2),rep(show.pol[length(show.pol)],2)),col=rgb(.1,.1,.1,.25),border='transparent')
    for(i in 1:nrow(dd)) points(Nx[1]:Nx[2],rep(i,Nx[2]),pch=21,cex=fn.scale(dd[i,],2.5),bg=rgb(.1,.1,.1,.4))
    axis(1,1:ncol(dd),colnames(dd))
    axis(2,1:nrow(dd),rownames(dd),las=1,cex.axis=.5)
    mtext("number of records per year",3,1,cex=1.25)
    Lab=round(quantile(dd,probs=c(.95,.995,1)))
    legend('topright',paste(Lab),pch=21,pt.bg="grey70",
           pt.cex=fn.scale(Lab,2.5),horiz=T,title="# of records")
    
    return(Selected)
  }
  
  if(pLot)pdf(paste(HndL,paste(SP,"_",what,sep=""),".pdf",sep="")) 
  
  #First, select vessels 
  Ves.Yr=with(d,table(vessel,finyear))
  this.ves=fn.plt(Ves.Yr)
  mtext("Vessel (all)",2,line=1.65,cex=1.25)
  
  #Second, select blocks for selected vessel
  BLK.Yr=with(subset(d,vessel%in%this.ves),table(blockx,finyear))
  this.blks=fn.plt(BLK.Yr)
  mtext("Block (for selected Vessels)",2,line=1.65,cex=1.25)
  
  
  #Third, keep only selected blocks and vessels
  d=subset(d,blockx%in%this.blks)
  d=subset(d,vessel%in%this.ves)
  
  BLK.Yr=with(d,table(blockx,finyear))
  this.blks=fn.plt(BLK.Yr)
  mtext("Block (for selected block and vessel)",2,line=1.65,cex=1.25)
  
  Ves.Yr=with(d,table(vessel,finyear))
  this.ves=fn.plt(Ves.Yr)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  
  MN.Yr=with(d,table(month,finyear))
  kk=fn.plt(MN.Yr)
  mtext("Month (for selected block and vessel)",2,line=1.65,cex=1.25)
  
  
  MN.Ves=with(d,table(vessel,month))
  kk=fn.plt(MN.Ves)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Month",1,line=1,cex=1.1)
  
  MN.blk=with(d,table(blockx,month))
  kk=fn.plt(MN.blk)
  mtext("Block (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Month",1,line=1,cex=1.1)
  
  ves.blk=with(d,table(vessel,blockx))
  kk=fn.plt(ves.blk)
  mtext("Vessel (for selected block and vessel)",2,line=1.65,cex=1.25)
  mtext("Block",1,line=1,cex=1.1)
  
  if(pLot)dev.off()
  
  
  return(list(this.blks=this.blks,this.ves=this.ves))
}
fn.show.blk=function(dat,CEX,SRt,dat.all) 
{
  dat=sort(dat)
  LAT.kept=sapply(dat, function(x) -as.numeric(substr(x, 1, 2)))
  LONG.kept=sapply(dat, function(x) 100+as.numeric(substr(x, 3, 4)))
  
  Y=-36:-26; X=seq(113,129,length.out=length(Y))
  plot(X,Y,ylab='',xlab="",col="transparent",cex.lab=1.5,cex.axis=1.25)
  for(e in 1:length(LAT.kept))
  {
    dd.y=c(LAT.kept[e]-1,LAT.kept[e]-1,LAT.kept[e],LAT.kept[e])
    dd.x=c(LONG.kept[e],LONG.kept[e]+1,LONG.kept[e]+1,LONG.kept[e])
    polygon(dd.x,dd.y,col=rgb(0, 0, 1,0.25), border=rgb(0, 0, 1,0.5))
    text(LONG.kept[e]+0.5,LAT.kept[e]-0.5,dat[e],cex=CEX,col=1,srt=SRt,font=2)
  }
  
  #add not used
  dat.all=sort(dat.all)
  id.not=which(!dat.all%in%dat)
  if(length(id.not)>0)
  {
    dat.not.used=dat.all[id.not]
    LAT.kept=sapply(dat.not.used, function(x) -as.numeric(substr(x, 1, 2)))
    LONG.kept=sapply(dat.not.used, function(x) 100+as.numeric(substr(x, 3, 4)))
    for(e in 1:length(LAT.kept))
    {
      dd.y=c(LAT.kept[e]-1,LAT.kept[e]-1,LAT.kept[e],LAT.kept[e])
      dd.x=c(LONG.kept[e],LONG.kept[e]+1,LONG.kept[e]+1,LONG.kept[e])
      polygon(dd.x,dd.y,col=rgb(1, 0, 0,0.01), border=rgb(1, 0, 0,0.1))
      text(LONG.kept[e]+0.5,LAT.kept[e]-0.5,dat[e],cex=CEX,col=scales::alpha(2,.5),srt=SRt,font=2)
    }
  }
    

}

# SHOW EFFECT OF USING km gn d VS km g h for gummy ----------------------------------------------
Get.Mns=function(d,grp,Vars,LGND,add.arrow)
{
  do.efforts=FALSE
  
  d=d[,match(c(grp,Vars),names(d))]
  d$cpue.d=d$Catch.Target/d$Km.Gillnet.Days.c
  d$cpue.h=d$Catch.Target/d$Km.Gillnet.Hours.c
  #d$cpue.h_shot=d$Catch.Target/(d$Km.Gillnet.Hours_shot.c)
  for(v in 1:length(Vars))
  {
    if(Vars[v]=="HOURS.c")
    {
      ddd=subset(d,SHOTS.c==1)
      B= ddd[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
        mutate(lowCL = mean - 1.96*sd/sqrt(n),
               uppCL = mean + 1.96*sd/sqrt(n)) %>%
        as.data.frame
      B$yr=as.numeric(substr(B$FINYEAR,1,4))
      
      ddd=subset(d,SHOTS.c==2)
      B2= ddd[,match(c(grp,Vars[v]),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
        mutate(lowCL = mean - 1.96*sd/sqrt(n),
               uppCL = mean + 1.96*sd/sqrt(n)) %>%
        as.data.frame
      B2$yr=as.numeric(substr(B2$FINYEAR,1,4)) 
      
      plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col=1,ylim=c(min(c(B$lowCL,B2$lowCL)),max(c(B$uppCL,B2$uppCL))))
      arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col=1)
      #if(add.arrow[v])
      #{
      #   Is=(nrow(B)-5):nrow(B)
      #    arrows(B$yr[1],B$mean[1],mean(B$yr[Is]),mean(B$mean[Is]),col=1,lwd=2)
      #  legend("bottomright",paste(round(mean(B$mean[Is])/B$mean[1],1),"fold",sep="-"),bty='n',cex=1)
      #  with(B[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.1,.8,.15)))
      #   }
      
      points(B2$yr,B2$mean,pch=19,col='grey65')
      B2$lowCL=ifelse(B2$lowCL==B2$mean,B2$lowCL*.999,B2$lowCL)
      B2$uppCL=ifelse(B2$uppCL==B2$mean,B2$uppCL*1.0001,B2$uppCL)
      arrows(x0=B2$yr, y0=B2$lowCL,x1=B2$yr, y1=B2$uppCL,code=3, angle=90, length=0.05, col='grey65')
      #  if(add.arrow[v])
      #  {
      #    Is=(nrow(B2)-5):nrow(B2)
      #   arrows(B2$yr[1],B2$mean[1],mean(B2$yr[Is]),mean(B2$mean[Is]),col='forestgreen',lwd=2)
      #   legend("topright",paste(round(mean(B2$mean[Is])/B2$mean[1],1),"fold",sep="-"),bty='n',cex=1)
      #    with(B2[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.8,.1,.15)))
      #  }
      legend("topleft",c("1 shot","2 shots"),text.col=c("black","grey65"),bty='n',cex=1.5)
      
    }
    if(Vars[v]=="SHOTS.c")
    {
      ddd=subset(d,SHOTS.c%in%c(1,2))
      
      B=ddd[,match(c(grp,Vars[v]),names(d))] %>%
        group_by(FINYEAR,SHOTS.c) %>%
        summarise (n = n()) %>%
        mutate(freq = n / sum(n))%>%
        as.data.frame
      B$yr=as.numeric(substr(B$FINYEAR,1,4))
      B1=reshape(subset(B,select=c(SHOTS.c,freq,yr)),
                 v.names = "freq", idvar = "SHOTS.c",
                 timevar = "yr", direction = "wide")
      barplot(as.matrix(B1[,-1]),names.arg=unique(B$yr),legend.text=c("1 shot","2 shots"))
      box()
      
    }
    if(do.efforts)if(!Vars[v]%in%c("HOURS.c","SHOTS.c")){
      # B= d[,match(c(grp,Vars[v]),names(d))] %>%
      #   na.omit()%>%
      #   group_by(FINYEAR) %>%
      #   summarise_all(funs(mean=mean,sd=sd,n=length)) %>%
      #   mutate(lowCL = mean - 1.96*sd/sqrt(n),
      #          uppCL = mean + 1.96*sd/sqrt(n)) %>%
      #   as.data.frame
      B= d[,match(c(grp,Vars[v]),names(d))] %>%
        na.omit()%>%
        group_by(FINYEAR) %>%
        summarise_all(funs(mean=sum)) %>%
        as.data.frame
      
      B$yr=as.numeric(substr(B$FINYEAR,1,4))
      plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col=1)
      #plot(B$yr,B$mean,xlab="",ylab=LGND[v],pch=19,col=1,ylim=c(min(B$lowCL),max(B$uppCL)))
      #arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col=1)
      #    if(add.arrow[v])
      #    {
      #     Is=(nrow(B)-5):nrow(B)
      #    arrows(B$yr[1],B$mean[1],mean(B$yr[Is]),mean(B$mean[Is]),col='black',lwd=2)
      #     legend("bottomright",paste(round(mean(B$mean[Is])/B$mean[1],1),"fold",sep="-"),bty='n',cex=1.5)
      #    with(B[Is,],polygon(x=c(yr,rev(yr)),y=c(rep(min(lowCL),length(Is)),rep(max(uppCL),length(Is))),col=rgb(.1,.1,.8,.15)))
      # }
    }
    
  }
  
  #cpues
  B= d[,match(c(grp,"cpue.d"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr,B$mean,xlab="",ylab="CPUE (kg gn d)",pch=19,col='black',ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr, y0=B$lowCL,x1=B$yr, y1=B$uppCL,code=3, angle=90, length=0.05, col='black')
  axis(2,col="black",col.axis = "black")
  
  par(new=T)
  B= d[,match(c(grp,"cpue.h"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  B$yr=as.numeric(substr(B$FINYEAR,1,4))
  plot(B$yr+.5,B$mean,xlab="",ylab="",pch=19,col='grey50',axes=F,ylim=c(min(B$lowCL),max(B$uppCL)))
  arrows(x0=B$yr+.5, y0=B$lowCL,x1=B$yr+.5, y1=B$uppCL,code=3, angle=90, length=0.05, col='grey50')
  
  
  # B= d[,match(c(grp,"cpue.h_shot"),names(d))]%>%group_by(FINYEAR) %>% summarise_all(funs(mean=mean,sd=sd,n=length))%>%
  #   mutate(lowCL = mean - 1.96*sd/sqrt(n),
  #          uppCL = mean + 1.96*sd/sqrt(n)) %>%
  #   as.data.frame
  # B$yr=as.numeric(substr(B$FINYEAR,1,4))
  # arrows(x0=B$yr+.5, y0=B$lowCL,x1=B$yr+.5, y1=B$uppCL,code=3, angle=90, length=0.05, col='black')
  # points(B$yr+.5,B$mean,xlab="",ylab="",pch=21,col='black',bg="white")
  axis(side = 4,col="grey50",col.axis = "grey50")
  mtext("CPUE (kg gn h)",4,line=1.5,col="grey50",cex=1,las=3)
  legend("top",c("kg/km gillnet days","kg/km gillnet hours"),
         bty='n',col=c("black","grey50"),cex=1.35,
         pt.bg=c("black","grey50"),pch=21)
}

# OUTPUT DATA TABLES ----------------------------------------------
Export.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                    body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                    HDR.names,HDR.span,HDR.2nd)
{
  mydoc = docx(Doc.nm)  #create r object
  mydoc = addSection( mydoc, landscape = T )   #landscape table
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                     header.cell.props = cellProperties(background.color=HdR.bg), 
                     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                        font.weight="bold",font.family =Fnt.hdr), 
                     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
  #Add header
  MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
  
  #Add second header
  MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
  
  
  # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
                                 inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                 outer.vertical = borderNone(),
                                 outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  mydoc = addFlexTable( mydoc, MyFTable)   
  mydoc = addSection( mydoc, landscape = F ) 
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}



# EFFICIENCY_CALCULATE EFFORT CREEP THRU TIME----------------------------------------------
fun.check.ves.char.on.cpue=function(d,NM)
{
  d=d%>%
    mutate(CPUE=Catch.Target/Km.Gillnet.Hours.c,
           finyear=substr(FINYEAR,1,4))%>%
    group_by(finyear)%>%
    mutate(Annual.rel.cpue=CPUE/max(CPUE))%>%
    ungroup()%>%
    data.frame()
  
  #temporal dynamics of vess char and cpue
  pdf(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/CPUE vs Vessel_char_thru time_',NM,'.pdf')),
      width=12)
  for(v in 1:length(ves.vars))
  {
    d$var=d[,ves.vars[v]]
    if(is.integer(d$var) | is.numeric(d$var)) d$var=factor(d$var,levels=sort(unique(d$var)))
    nN=unique(d%>%filter(!is.na(var))%>%pull(VESSEL))
    p=d%>%
      filter(VESSEL%in%nN)%>%
      filter(!is.na(var))%>%
      ggplot(aes(finyear,Annual.rel.cpue,fill=var))+
      geom_boxplot()+
      ylim(0,quantile(d$Annual.rel.cpue,probs=0.99))+
      scale_fill_manual(name = paste0(ves.vars[v],' (',length(nN),' vessels)'),
                        values=colorRampPalette(c("cadetblue1", "dodgerblue4"))(length(unique(d$var))))+
      theme_PA()+
      theme(legend.position = 'top',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p)
  }
  dev.off()
  
  #vess char correlation
  x=d[,ves.vars]
  x=x[!rowSums(!is.na(x)) == 0,]
  
  pdf(handl_OneDrive(paste0('Analyses/Catch and effort/Outputs/Efficiency creep/CPUE vs Vessel_char_correlation_',NM,'.pdf')))
  ggcorrplot(x%>%correlation(), lab = T,type='upper', show.diag = F)
  dev.off()
  
}

# INFLUENCE PLOTS ---------------------------------------------------------
Influence.fn=function(MOD,DAT,Term.type,termS,add.Influence,SCALER)  
{
  termS=subset(termS,termS%in%names(Term.type))
  #extract main term coefficients for each species
  nt=length(termS)
  Store1=Store2=MatcH=COEF.list=COEF.SE.list=vector('list',nt)
  ID=c(1,grep("[:]", names(coef(MOD))))
  Cofs=coef(MOD)[-ID]
  Cofs.intercept=coef(MOD)[ID]
  if(class(MOD)[1]=="glm")
  {
    Cofs.SE=summary(MOD)$coefficients[-ID, 2]
    Cofs.SE.intercept=summary(MOD)$coefficients[ID, 2]
  }
  
  if(class(MOD)[1]=="gam")
  {
    Cofs.SE=summary(MOD)$se[-ID]
    Cofs.SE.intercept=summary(MOD)$se[ID]
  }
  
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      Store1[[p]]=as.character(levels(DAT[,match(termS[p],names(DAT))]))[-1]
      Store2[[p]]=paste(termS[p],Store1[[p]],sep="")
    }
  }
  for(p in 1:nt)MatcH[[p]]=if (Term.type[p]=="CAT") match(Store2[[p]],names(Cofs))
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") COEF.list[[p]]=Cofs[MatcH[[p]]]
    if (Term.type[p]=="Cont") 
    {
      if(class(MOD)[1]=="glm") COEF.list[[p]]=Cofs[match(termS[p],names(Cofs))]
      if(class(MOD)[1]=="gam")
      {
        aa=dummy.coef(MOD)
        COEF.list[[p]]=aa[[match(termS[p],names(aa))]]
      }
    }
  }
  for(p in 1:nt) 
  {
    if (Term.type[p]=="CAT")
    {
      COEF.list[[p]]=data.frame(Dummy=Store1[[p]],coef=COEF.list[[p]])
      COEF.list[[p]]$Dummy=as.character(COEF.list[[p]]$Dummy)
    }
  }
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      A=as.character(levels(DAT[,match(termS[p],names(DAT))]))[1]
      COEF.list[[p]]=rbind(COEF.list[[p]],data.frame(Dummy=A,coef=0))
    }
  }
  for(p in 1:nt) if (Term.type[p]=="CAT")colnames(COEF.list[[p]])=c(termS[p],paste("Coef.",termS[p],sep=""))
  
  #attach coefficients to data
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT") DAT=merge(DAT,COEF.list[[p]],by=termS[p],all.x=T)
    if (Term.type[p]=="Cont")
    {
      DAT=cbind(DAT,COEF.list[[p]]*DAT[,match(names(COEF.list[[p]]),names(DAT))]) #coef X value
      colnames(DAT)[ncol(DAT)]=paste("Coef.",termS[p],sep="")      
    }
  }
  
  #ny
  ny=table(DAT$finyear)
  
  #calculate rho
  Coef.vec=match(paste("Coef.",termS,sep=""),names(DAT))
  Mean.coef=Annual.Dev=vector('list',nt)
  names(Annual.Dev)=termS
  Over.all.influence=rep(NA,nt)
  names(Over.all.influence)=termS
  for(p in 1:nt) Mean.coef[[p]]=mean(DAT[,Coef.vec[p]],na.rm=T)
  
  #calculate overall and annual deviation from mean (i.e. influence)
  for(p in 1:nt)
  {
    #Calculate lambda y
    dev=rep(NA,length(ny))
    for(t in 1:length(ny))
    {
      a=subset(DAT,finyear==names(ny[t]))
      dev[t]=(sum(a[,Coef.vec[p]]-Mean.coef[[p]]))/ny[t]
    }  
    
    #Store Annual deviance
    #note: exp because it's multiplicative
    Annual.Dev[[p]]=exp(dev)  
    
    #Store Overall influence of variable
    Over.all.influence[p]=exp(sum(abs(dev))/length(ny))-1    
  }
  
  #plot CDI (categorical vars only)
  STORe=vector('list',nt)
  names(STORe)=termS
  for(p in 1:nt)
  {
    if (Term.type[p]=="CAT")
    {
      COEF=COEF.list[[p]][,2]
      names(COEF)=COEF.list[[p]][,1]
      #add intercept
      COEF=Cofs.intercept+COEF      
      COEF.SE=c(Cofs.SE.intercept,sqrt(Cofs.SE.intercept^2+Cofs.SE[MatcH[[p]]]^2))
      names(COEF.SE)=names(COEF)
      
      COEF=sort(COEF)
      COEF.SE=COEF.SE[match(names(COEF),names(COEF.SE))]
      x=1:length(COEF)
      
      if(length(x)>1)
      {
        nf <- layout(matrix(c(1,1,1,2,2,2), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
        if(add.Influence=="YES")nf <- layout(matrix(c(1,1,0,2,2,3), 2, 3, byrow = TRUE), widths=c(1.5,1),heights=c(1,1.5))
        par(mar=c(0,0,0,0),oma=c(4,6,1,1),las=1,mgp=c(1,.9,0))
        #layout.show(nf)
        
        # Coefficients
        minSE=COEF-COEF.SE
        maxSE=COEF+COEF.SE
        plot(x,COEF,xlab="",xaxt="n",ylim=c(min(minSE),max(maxSE)),cex.axis=1.25,pch=19,cex=2)
        arrows(x, minSE, x, maxSE, code=3, angle=90, length=0.1)
        axis(1,1:length(COEF),F,tcl=0.5)
        axis(1,seq(1,length(COEF),2),F,cex.axis=1.15,tcl=1)      
        mtext("Coefficient",side=2,line=4,cex=1.5,las=3)
        
        # Bubble plot of records
        TAb=table(DAT$finyear,DAT[,match(termS[p],names(DAT))])
        Prop.Rec=TAb/rowSums(TAb)
        Nombres=gsub("[^[:digit:]]", "", names(COEF))
        if(termS[p]=="vessel") Nombres=1:length(Nombres)   #change vessel name for dummy
        Prop.Rec=Prop.Rec[,match(names(COEF),colnames(Prop.Rec))]
        bubble.plot(x,1:length(ny),Prop.Rec,scaler=SCALER,termS[p],"Financial year")
        axis(1,1:length(COEF),F,tck=-0.015)
        axis(1,seq(1,length(COEF),1),Nombres[seq(1,length(COEF),1)],cex.axis=1.15,tck=-0.025)
        axis(2,1:length(ny),F,tck=-0.015)
        axis(2,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1,tck=-0.025)   
        mtext("Financial year",side=2,line=4.35,cex=1.5,las=3)
        mtext(termS[p],side=1,line=2.5,cex=1.5)
        
        if(add.Influence=="YES")
        {
          #Influence plot
          plot(Annual.Dev[[p]],1:length(ny),type="o",pch=19,xlab="",ylab="",cex=2,cex.axis=1.25,yaxt='n')
          abline(v=1,lty=3,col=1)
          mtext("Influence",side=1,line=2.5,cex=1.5)      
          axis(2,1:length(ny),F,tcl=0.5)
          axis(2,seq(1,length(ny),2),F,tcl=1)
        }
        
        STORe[[p]]=list(x=x,COEF=COEF,COEF.SE=COEF.SE,Prop.Rec=Prop.Rec)
      }
    }
  }
  return(list(Annual.Dev=Annual.Dev,ny=ny,Over.all.influence=Over.all.influence,store=STORe,
              termS=termS,ny=ny,SCALER=SCALER))
}
Compare.term.infl.fun=function(A,WHERE,WHERE2,YLIM,spliT)
{
  Annual.Dev=A$Annual.Dev
  ny=A$ny
  if(is.null(YLIM))YLIM=c(min(unlist(lapply(Annual.Dev,min))),max(unlist(lapply(Annual.Dev,max))))
  NamE=names(A$Over.all.influence)
  
  nt=length(Annual.Dev)  
  LTY=c(1,4,3,1,3,2)
  plot(1:length(ny),Annual.Dev[[1]],col=LTY.col[1],type="l",xlab="",ylab="",lwd=LWD,
       cex.axis=1.35,xaxt='n',ylim=YLIM)
  abline(h=1,lty=3,col=1)
  axis(1,1:length(ny),F,tck=-0.02)
  axis(1,seq(1,length(ny),2),F,tck=-0.04)
  axis(1,seq(1,length(ny),2),names(ny)[seq(1,length(ny),2)],cex.axis=1.35,tck=-0.04)
  for(p in 2:nt)lines(1:length(ny),Annual.Dev[[p]],lwd=LWD,lty=LTY[p],col=LTY.col[p])
  LEG=paste(NamE," (",round(100*A$Over.all.influence,1),"%)",sep="")
  if(spliT=="YES")
  {
    nn=1:(length(LEG)/2)
    legend(WHERE,LEG[nn],bty='n',lty=LTY[nn],col=LTY.col[nn],lwd=LWD,cex=1.25,pt.cex=1.5)
    nn=1+(length(LEG)/2):length(LEG)
    legend(WHERE2,LEG[nn],bty='n',lty=LTY[nn],col=LTY.col[nn],lwd=LWD,cex=1.25,pt.cex=1.5)
  }else
    legend(WHERE,LEG,bty='n',lty=LTY,col=LTY.col,lwd=LWD,cex=1.25,pt.cex=1.5)
}
Fig.CDI.paper.fn=function(store,SCALER,termS)
{
  ny=rownames(store[[1]]$Prop.Rec)
  nt=length(termS)
  YLABs=termS
  YLABs=ifelse(YLABs=="blockx","Block",ifelse(YLABs=="vessel","Vessel",
                                              ifelse(YLABs=="month","Month",NA)))
  for(p in 1:nt)
  {
    x=store[[p]]$x
    COEF=store[[p]]$COEF
    COEF.SE=store[[p]]$COEF.SE
    Prop.Rec=store[[p]]$Prop.Rec
    
    # Coefficients
    minSE=COEF-COEF.SE
    maxSE=COEF+COEF.SE
    plot(x,COEF,xlab="",ylab="",xaxt="n",ylim=c(min(minSE),max(maxSE)),cex.axis=1.25,pch=19,cex=1.75)
    arrows(x, minSE, x, maxSE, code=3, angle=90, length=0.1)
    axis(1,1:length(COEF),F,tcl=0.5)
    axis(1,seq(1,length(COEF),2),F,cex.axis=1.15,tcl=1)      
    
    #Bubbleplot
    Nombres=gsub("[^[:digit:]]", "", names(COEF))
    if(termS[p]=="vessel") Nombres=1:length(Nombres)   #change vessel name for dummy
    Prop.Rec=Prop.Rec[,match(names(COEF),colnames(Prop.Rec))]
    bubble.plot(x,1:length(ny),Prop.Rec,scaler=SCALER,"","")
    axis(1,1:length(COEF),F,tck=-0.015)
    axis(1,seq(1,length(COEF),1),Nombres[seq(1,length(COEF),1)],cex.axis=1.15,tck=-0.025)
    axis(2,1:length(ny),F,tck=-0.015)
    axis(2,seq(1,length(ny),2),F,cex.axis=1,tck=-0.025) 
    if(p==1)axis(2,seq(1,length(ny),2),ny[seq(1,length(ny),2)],cex.axis=1,tck=-0.025)   
    mtext(YLABs[p],side=1,line=1.75,cex=1)
  }
}
bubble.plot=function(x,y,z,scaler,Xlab,Ylab)  
{
  xo=outer(x,rep(1,length=length(y)))
  yo=t(outer(y,rep(1,length=length(x))))
  zo=z
  for(zz in 1:nrow(zo))zo[zz,]=((zo[zz,]/max(zo[zz,]))^0.5)*scaler
  matplot(xo,yo,type="n",xlab=Xlab,ylab=Ylab,xaxt='n',yaxt='n')
  for(s in 1:length(x))
  {
    points(xo[s,],yo[s,],cex=zo[,s],pch=16,col="grey80")
    points(xo[s,],yo[s,],cex=zo[,s],pch=1,col="black")
  }
}

# Compare all different ways of calculating cpues -----------------------------------------------------------------------
fn.plot.all.indices=function(sp)
{
  d1=Unstandardised.creep[[match(sp, names(Unstandardised.creep))]]%>%
    mutate(method='Unstandardised with creep',
           mean=response/mean(response),
           year=as.numeric(substr(finyear,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='1.monthly')
  d2=Raw.index[[match(sp, names(Raw.index))]]%>%
    mutate(method='raw',
           mean=mean/mean(mean),
           year=as.numeric(substr(FINYEAR,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='1.monthly')
  d3=Store_nom_cpues_monthly[[match(sp, names(Store_nom_cpues_monthly))]]$CPUE.All$km.gillnet.hours.c%>%
    rename(year=season)%>%
    group_by(method)%>%
    mutate(mean1=mean(mean))%>%
    ungroup()%>%
    mutate(mean=mean/mean1)%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='1.monthly')
  
  d4=Pred[[match(sp, names(Pred))]]%>%
    mutate(method='standardised with creep',
           mean=response/mean(response),
           year=as.numeric(substr(finyear,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='1.monthly')
  
  d1.daily=Unstandardised.daily.creep[[match(sp, names(Unstandardised.daily.creep))]]%>%
    mutate(method='Unstandardised with creep',
           mean=response/mean(response),
           year=as.numeric(substr(finyear,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='2.daily')
  d2.daily=Raw.index.daily[[match(sp, names(Raw.index.daily))]]%>%
    mutate(method='raw',
           mean=mean/mean(mean),
           year=as.numeric(substr(FINYEAR,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='2.daily')
  d3.daily=Store_nom_cpues_daily[[match(sp, names(Store_nom_cpues_daily))]]$CPUE.All$km.gillnet.hours.c%>%
    rename(year=season)%>%
    group_by(method)%>%
    mutate(mean1=mean(mean))%>%
    ungroup()%>%
    mutate(mean=mean/mean1)%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='2.daily')
  
  d4.daily=Pred.daily[[match(sp, names(Pred.daily))]]%>%
    mutate(method='standardised with creep',
           mean=response/mean(response),
           year=as.numeric(substr(finyear,1,4)))%>%
    dplyr::select(year,method,mean)%>%
    mutate(period='2.daily')
  
  
  p1=rbind(d1,d2,d3,d4)%>%
    ggplot(aes(year,mean,color=method))+
    geom_point(size=2.5)+geom_line()+
    facet_wrap(~period)+
    theme(legend.position = 'top')+ylab('Normalised cpue')
  p2=rbind(d1.daily,d2.daily,d3.daily,d4.daily)%>%
    ggplot(aes(year,mean,color=method))+
    geom_point(size=2.5)+geom_line()+
    facet_wrap(~period)+
    theme(legend.position = 'top')+ylab('Normalised cpue')
  p=ggarrange(p1,p2,ncol=1)
  print(p)
}