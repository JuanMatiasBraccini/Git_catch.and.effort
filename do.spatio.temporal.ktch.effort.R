
  setwd(handl_OneDrive("Analyses/Catch and effort/Outputs/Spatio.temporal_Catch"))
  
  # Plot spatio-temporal catches of 4 main shark species
  do.bubble.plots=TRUE
  if(do.bubble.plots)
  {
    Indicator.ranges=list(Gummy.range,Whiskery.range,Dusky.range,Dusky.range,Sandbar.range)
    Effect.area=function(SPEC,RANGO)
    {
      a=subset(Data.monthly,SPECIES%in%SPEC)
      yr=sort(unique(a$FINYEAR))
      n=length(yr)
      
      agg=aggregate(LIVEWT.c~FINYEAR+BLOCKX+LAT+LONG,a,sum)
      
      for(i in 1:n)
      {
        b=subset(agg,FINYEAR==yr[i])
        z=b$LIVEWT.c/max(b$LIVEWT.c)
        plot(b$LONG,b$LAT,cex=z*2,main=yr[i],ylab="",xlab="",pch=19,ylim=c(-36,-26),
             xlim=c(113,129),cex.axis=0.8,cex.main=.85,col="steelblue4")
        legend("top",paste(round(max(b$LIVEWT.c)/1000),"tons"),pch=19,pt.cex=2,bty='n',col="steelblue4")
        Y1=RANGO[RANGO<0]
        Y2=-36
        if(length(Y1)>0)
        {
          X2=RANGO[RANGO>0]
          X1=112.5
        }
        if(length(Y1)==0)
        {
          Y1=-31
          X1=RANGO[1]
          X2=RANGO[2]
        }
        polygon(c(X1,X2,X2,X1),c(Y2,Y2,Y1,Y1),col=rgb(0.1, .1, .1, 0.1),border="transparent")
      }
    }
    for(i in 1:length(Indicator.species))
    {
      print(paste('Catch bubble plots for ----',names(Indicator.species)[i]))
      NM=names(Indicator.species)[i]
      tiff(paste0("Bubble.plot.catch_",NM,".tiff"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      n.graph=length(unique(subset(Data.monthly,SPECIES==Indicator.species[i])$FINYEAR))+1
      smart.par(n.plots=n.graph,MAR=c(1,1.5,1.5,1.5),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
      Effect.area(Indicator.species[i],Indicator.ranges[[i]])
      plot(1:10,col="transparent",xaxt='n',ann=F,yaxt='n',fg="white")
      text(5,7,NM,cex=1.25)
      text(5,4,"shark",cex=1.25)
      dev.off()
    }
  }
  
  #Catch densities
  do.density=FALSE
  if(do.density)
  {
    for(i in 1:length(Indicator.species))
    {
      print(paste('Catch density plots for ----',names(Indicator.species)[i]))
      tiff(paste0("Catch.density_",names(Indicator.species)[i],".tiff"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
      fn1=function(x) as.numeric(as.character(cut((x-trunc(x)),breaks=seq(0,60,10)/60,labels=seq(0,50,10)/60)))
      a1=Data.daily%>%filter(SPECIES==Indicator.species[i])%>%
        mutate(LAT=abs(LAT),
               LAT10=-round(trunc(LAT)+fn1(LAT+(LatMin/60)),2),
               LONG10=round(trunc(LONG)+fn1(LONG+(LongMin/60)),2))%>%
        filter(!is.na(LONG10)&!is.na(LAT10))%>%
        group_by(LAT10,LONG10)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c))%>%
        spread(LAT10,LIVEWT.c,fill=0)
      
      a2=as.matrix(a1[,-1])
      x=a1$LONG10
      y=as.numeric(colnames(a1)[-1])
      filled.contour(x,y,a2,
                     plot.axes = { axis(1,seq(trunc(min(x)),trunc(max(x)),by=1) )
                       axis(2, seq(trunc(min(y)),trunc(max(y)),by=1))},
                     key.axes = axis(4, seq(round(min(a2)), round(max(a2)), by = 5)))
      dev.off()
    }
    
  }
  
  #define what effort to plot
  Ef.var="days"   #express effort in km gn days
  #Ef.var="hours"
  
  #define spatial variables for plotting
  data(worldLLhigh)
  South.WA.lat=c(-36,-25); South.WA.long=c(112,129)
  S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
  S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
  OZ.lat=c(-44.5,-11);OZ.long=c(113,155)
  Perth=c(115.866,-31.95)
  Rotnest=c(115.50,-32.02)
  Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
  Lat.seq=c(-26,-28,-30,-32,-34)
  numInt=20  #number of breaks 
  in.color="YES"
  if(!in.color=="YES")couleurs=rev(gray(seq(0.2,0.9,length=(numInt-1))))
  if(in.color=="YES")
  {
    couleurs=rev(heat.colors((numInt-1)))
    Colfunc <- colorRampPalette(c("yellow","red"))
    Couleurs=Colfunc(numInt-1)
  }
  couleurs=c("white",couleurs)
  numberLab=5
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  a=112:129
  b=seq(-37,South.WA.lat[2],length.out=length(a))
  PLATE=c(.01,.9,.075,.9)
  LATT=South.WA.lat[2]:South.WA.lat[1]
  LONGG=South.WA.long[1]:South.WA.long[2]
  #South.WA.lat=c(-37,-25)
  Lat.seq=c(-26,-28,-30,-32,-34,-36)
  OZ.lat=c(-44.5,-11);OZ.long=c(113,155)
  Perth=c(115.866,-31.95)
  Rotnest=c(115.50,-32.02)
  Lat.exp=expression(paste("Latitude (",degree,"S)",sep=""))
  Lon.exp=expression(paste("Longitude (",degree,"E)",sep=""))
  
  #grouping years
  Last.calendar.Yr=as.numeric(substr(Current.yr,1,4))
  yrs.grupd=4
  Yr.group=seq(1975,(Last.calendar.Yr+1),5)
  Yr.group.plus=Yr.group+yrs.grupd
  Yr.group.plus[length(Yr.group.plus)]=min(Yr.group.plus[length(Yr.group.plus)],Last.calendar.Yr)
  Yr.range=vector('list',length(Yr.group))
  for(e in 1:length(Yr.group))
  {
    if(Yr.group[e]<Yr.group.plus[e]) Yr.range[[e]]=paste(Yr.group[e],"-",Yr.group.plus[e])
  }
  
  Yr.range=unlist(Yr.range)  
  
  #create data list
  DATA.lista=vector('list',length(Yr.range))
  names(DATA.lista)=Yr.range
  for(i in 1:(length(DATA.lista)))
  {
    DATA.lista[[i]]=  subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                               METHOD%in%c("GN","LL")& LAT<=(-26) & 
                               YEAR.c>=Yr.group[i] & YEAR.c<=Yr.group.plus[i])
  }  
  
  add.depth="NO"
  
  #This bit is irrelevant
  do.this=FALSE
  if(do.this)
  {
    #add effort to Data.monthly.GN 
    s=subset(Eff,LAT<=(-26))
    s$METHOD="GN"
    s$Same.return=with(s,paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))
    s=s[,-match(c("BLOCKX","FINYEAR","MONTH","VESSEL","METHOD","LAT","LONG","Eff.Reporter"),names(s))]
    Data.monthly.GN=Data.monthly.GN%>%left_join(s,by="Same.return")
    #Data.monthly.GN=merge(Data.monthly.GN,s,by="Same.return",all.x=T)
    
    Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
    xbat=sort(unique(Bathymetry$V1))
    ybat=sort(unique(Bathymetry$V2))
    if(add.depth=="YES")reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
    data(worldLLhigh)
    #define coordinates of plots
    South.WA.lat=c(-36,-25); South.WA.long=c(112,129)
    S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
    S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
    Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
    Lat.seq=c(-26,-28,-30,-32,-34)
    numInt=20  #number of intervals for effort 
    #couleurs=rev(gray(0:(numInt-1)/(numInt-1)))
    couleurs=rev(gray(seq(0.2,0.9,length=numInt)))
    #couleurs  <- tail(topo.colors(trunc(1.4 * numInt)),numInt)
    numberLab=5
    colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
    a=112:129
    b=seq(-37,South.WA.lat[2],length.out=length(a))
    #max block effort for each year period
    Ef.var="days"   #express effort in km gn days
    #Ef.var="hours"
    MAX.EFF=data.frame(Period=Yr.range,Max.eff=NA)
    EFF.bin=NULL
    for (i in 1:length(DATA.lista))
    {
      DATA=DATA.lista[[i]]
      DATA=DATA[-which(duplicated(DATA$Same.return)),]
      if(Ef.var=="days") Max.effort=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX),FUN=sum,na.rm=T))
      if(Ef.var=="hours") Max.effort=with(DATA,aggregate(Km.Gillnet.Hours.c,list(BLOCKX),FUN=sum,na.rm=T))
      EFF.bin=c(EFF.bin,Max.effort[,2])
      MAX.EFF[i,2]=max(Max.effort[,2])
    }
    Max.effort=ceiling(max(MAX.EFF[,2])) 
    # dummy=cbind(128.5,-35.5,Max.effort)
    fn.eff.plot=function(DATA,tcl.1,tcl.2,EffortBreaks)   #function plot effort
    {
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      DATA$BLOCKX.c=with(DATA,paste(-LAT,LONG-100,sep=""))
      MapEffort=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX.c),FUN=sum,na.rm=T))    
      
      colnames(MapEffort)=c("BLOCKX.c","Total effort")
      id=unique(match(MapEffort$BLOCKX.c,DATA$BLOCKX.c))
      MapEffort$LAT=DATA$LAT[id]
      MapEffort$LONG=DATA$LONG[id]
      lat=sort(unique(MapEffort$LAT))
      ids=match(-27,lat)
      if(is.na(ids))      #add dummy for sorting image out of whack when missing Lat
      {
        adD=MapEffort[1,]
        adD[,2]=NA
        adD$BLOCKX.c="-27 113"
        adD$LAT=-27
        MapEffort=rbind(MapEffort,adD)
      }
      
      MapEffort$LAT.cen=MapEffort$LAT-.5
      MapEffort$LONG.cen=MapEffort$LONG+.5  
      MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
      MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapEffort$LONG.cen))
      lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image  
      
      MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Total effort"),names(MapEffort))]  
      #MapEffort=rbind(MapEffort,dummy)#keep in perspective
      
      Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",    #transposed as matrix   
                                 timevar="LAT.cen",v.names="Total effort", direction="wide"))  
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]  									
      
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=EffortBreaks,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      
      if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                   nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
    }
    EffortBreakSS=quantile(EFF.bin,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
  }
  
  
  #--Catch
  
  #catch breaks
  sp=c("whiskery","gummy","dusky","sandbar")
  MAX.CATCH=data.frame(Period=Yr.range,Max.catch=NA)
  CATCH.bin=NULL
  fn.catch.breaks=function(SP)
  {
    for (j in 1:length(DATA.lista))
    {
      if(!is.null(SP))DATA=subset(DATA.lista[[j]],LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9) & SPECIES%in%SP)
      if(is.null(SP)) DATA=subset(DATA.lista[[j]],LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9))
      if(nrow(DATA)>0)
      {
        Max.catch=with(DATA,aggregate(LIVEWT.c,list(BLOCKX),FUN=sum,na.rm=T))
        CATCH.bin=c(CATCH.bin,Max.catch[,2])
        MAX.CATCH[j,2]=max(Max.catch[,2])
      }
    }
    Breaks=quantile(CATCH.bin,probs=seq(0,1,1/numInt),na.rm=T)  
    return(Breaks)
  }
  
  # 1. Plot TDGDLF catch by 5-year period and 1 nm block  
  fn.catch.plot=function(DATA,SP,tcl.1,tcl.2,BREAKS)
  {
    if(!is.null(SP))DATA=subset(DATA,LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9) & SPECIES%in%SP)
    if(is.null(SP)) DATA=subset(DATA,LAT<=(-26) & LAT>(-36.1)&LONG<=(129) & LONG >(111.9))
    DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
    DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    
    if(nrow(DATA)<=2)   plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    if(nrow(DATA)>2)
    {
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      
      MapCatch=with(DATA,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))    
      colnames(MapCatch)=c("BLOCKX.c","Total Catch")
      id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
      MapCatch$LAT=DATA$LAT[id]
      MapCatch$LONG=DATA$LONG[id]
      
      MapCatch$LAT.cen=MapCatch$LAT-.5
      MapCatch$LONG.cen=MapCatch$LONG+.5  
      MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
      MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapCatch$LONG.cen))
      lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
      
      MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
      
      
      Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]										
      
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=BREAKS,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      
    }
    
    if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-300),
                                 nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
  }
  hndl.sptl.ktch=handl_OneDrive("Analyses/Catch and effort/Outputs/Spatio.temporal_Catch/")
  Tar=TARGETS
  Tar[[3]]=18003
  #1.1 by indicator species
  for (ss in 1:length(Tar))
  {
    print(paste('Heat_map_5.yr.group. for ----',Tar[ss]))
    fn.fig(paste(hndl.sptl.ktch,"Heat_map_5.yr.group.",sp[ss],sep=""),2400, 2400)
    opar <- smart.par(n.plots=length(DATA.lista),MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(.1, 0.15, 0))
    Breaks=fn.catch.breaks(Tar[ss])
    for (i in 1:length(DATA.lista))
    {
      fn.catch.plot(DATA.lista[[i]],Tar[ss],tcl.1=.5,tcl.2=.5,Breaks)
      #fn.catch.plot(DATA.lista[[i]],Tar[ss],tcl.1=16,tcl.2=18.15,Breaks)
      mtext(Yr.range[i],side=3,line=-2,cex=.95)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .5,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .5,las=2,cex.axis=1,hadj=1.1)
      if(i==length(DATA.lista))color.legend(126,-26,129,-30.5,paste(round(Breaks/1000,0),'t'),
                                            rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
      # if(i==1)
      # {
      #   par(fig=c(0.2,.35,.8,0.95), new = T,mgp=c(.25,.2,0),las=1)
      #   plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
      #           col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
      #   box()
      #   polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
      #   par(opar)
      # }
    }
    mtext(Lat.exp,side=2,line=0.4,las=3,cex=1.3,outer=T) 
    mtext(Lon.exp,side=1,line=0.6,cex=1.3,outer=T)
    dev.off()
    rm(Breaks)
  }
  
  #1.2 total catch
  fn.fig(paste(hndl.sptl.ktch,"Heat_map_5.yr.group.","Total",sep=""),2400, 2400)
  opar <- smart.par(n.plots=length(DATA.lista),MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(.1, 0.15, 0))
  Breaks=fn.catch.breaks(NULL)
  for (i in 1:length(DATA.lista))
  {
    fn.catch.plot(DATA.lista[[i]],NULL,tcl.1=.5,tcl.2=.5,Breaks)
    mtext(Yr.range[i],side=3,line=-2,cex=.95)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .5,las=1,cex.axis=1,padj=-.15)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .5,las=2,cex.axis=1,hadj=1.1)
    if(i==length(DATA.lista))color.legend(126,-26,129,-30.5,paste(round(Breaks/1000,0),'t'),
                                          rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
  }
  mtext(Lat.exp,side=2,line=0.4,las=3,cex=1.3,outer=T) 
  mtext(Lon.exp,side=1,line=0.6,cex=1.3,outer=T)
  dev.off()
  rm(Breaks)
  
  # 2. Plot TDGDLF catch by year and 1 nm block
  Couleurs=couleurs
  fn.ctch.plot.all.yrs=function(DATA,tcl.1,tcl.2,numInt) 
  {
    DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
    DATA$LONG=as.numeric(DATA$LONG)
    DATA$blk=substr(DATA$BLOCKX,1,4)
    A=aggregate(LIVEWT.c~FINYEAR+blk,DATA,sum)
    Ymax=max(A$LIVEWT.c)
    Ymin=min(A$LIVEWT.c)
    #Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
    Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
    DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
    DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
    for(y in 1:length(FINYrS))
    {
      A=subset(DATA,FINYEAR==FINYrS[y])
      MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
      colnames(MapCatch)=c("BLOCKX.c","Total Catch")
      id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
      MapCatch$LAT=DATA$LAT[id]
      MapCatch$LONG=DATA$LONG[id]
      msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
      msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
      if(length(msn.lat)>0)
      {
        dummy=MapCatch[1:length(msn.lat),]
        dummy$`Total Catch`=0
        dummy$LAT=msn.lat
        dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
        MapCatch=rbind(MapCatch,dummy)
      }
      
      if(unique(min(A$YEAR.c))>2007)
      {
        MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
      }
      MapCatch$LAT.cen=MapCatch$LAT-.5
      MapCatch$LONG.cen=MapCatch$LONG+.5  
      MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
      MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapCatch$LONG.cen))
      lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
      MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
      Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      legend('top',FINYrS[y],bty='n',cex=1.2)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    plot(a,b,ann=F,axes=F,col='transparent')
    color.legend(quantile(a,probs=.25),quantile(b,probs=.95),quantile(a,probs=.6),quantile(b,probs=.25),
                 paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.75)
  }
  fn.ctch.plot.grouped.yrs=function(DATA,tcl.1,tcl.2,numInt,grouping) 
  {
    DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
    DATA$LONG=as.numeric(DATA$LONG)
    DATA$blk=substr(DATA$BLOCKX,1,4)
    DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
    DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
    
    DATA=subset(DATA,FINYEAR%in%FINYrS)
    DATA$yyr=as.numeric(substr(DATA$FINYEAR,1,4))
    
    AA=vector('list',length(FINYrS.gped))
    for(y in 1:length(FINYrS.gped))
    {
      A=aggregate(LIVEWT.c~blk,subset(DATA,FINYEAR%in%FINYrS.gped[[y]]),sum,na.rm=T)
      A$y.group=y
      AA[[y]]=A
    }
    A=do.call(rbind,AA)
    
    Ymax=max(A$LIVEWT.c)
    Ymin=min(A$LIVEWT.c)
    Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
    
    couleurs=rev(heat.colors(numInt))
    
    for(y in 1:length(FINYrS.gped))
    {
      A=subset(DATA,FINYEAR%in%FINYrS.gped[[y]])
      MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
      colnames(MapCatch)=c("BLOCKX.c","Total Catch")
      id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
      MapCatch$LAT=DATA$LAT[id]
      MapCatch$LONG=DATA$LONG[id]
      msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
      msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
      if(length(msn.lat)>0)
      {
        dummy=MapCatch[1:length(msn.lat),]
        dummy$`Total Catch`=0
        dummy$LAT=msn.lat
        dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
        MapCatch=rbind(MapCatch,dummy)
      }
      
      if(unique(min(A$YEAR.c))>2007)
      {
        MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
      }
      MapCatch$LAT.cen=MapCatch$LAT-.5
      MapCatch$LONG.cen=MapCatch$LONG+.5  
      MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
      MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapCatch$LONG.cen))
      lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
      MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
      Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      text(115.60,-31.96,"Perth",pos=4,cex=1.2)
      legend('top',names(FINYrS.gped)[y],bty='n',cex=1.2)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    plot(a,b,ann=F,axes=F,col='transparent')
    color.legend(quantile(a,probs=.25),quantile(b,probs=.95),quantile(a,probs=.6),quantile(b,probs=.25),
                 paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.75)
  }
  HnD.ctch.exp=handl_OneDrive("Analyses/Catch and effort/Outputs/Spatio.temporal_Catch")
  
  #2.1 by indicator species
  for(i in 1:length(Tar))
  {
    print(paste('Heat_map for ----',Tar[i]))
    
    ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                 METHOD%in%c("GN","LL")& LAT<=(-26) & SPECIES%in%Tar[[i]])
    FINYrS=unique(ddd$FINYEAR)
    
    #all years
    fn.fig(paste(HnD.ctch.exp,paste('Heat_map',names(Tar)[i],sep='_'),sep="/"),2400, 2400)
    smart.par(n.plots=length(FINYrS)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
    fn.ctch.plot.all.yrs(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=numInt)
    mtext(Lat.exp,side=2,line=0.5,las=3,cex=1.1,outer=T)
    mtext(Lon.exp,side=1,line=0.5,cex=1.1,outer=T)
    dev.off()
    
    # #grouped years
    # grouping=5
    # FINYrS.gp=seq(1,length(FINYrS),by=grouping)
    # FINYrS.gped=vector('list',length(FINYrS.gp))
    # for(f in 1:length(FINYrS.gped))
    # {
    #   if(f==length(FINYrS.gped))
    #   {
    #     FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
    #     if(length(FINYrS.gped[[f]]==1))
    #     {
    #       names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
    #     }else
    #     {
    #       names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
    #     }
    #     
    #   }else
    #   {
    #     FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
    #     names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
    #   }
    # }
    # fn.fig(paste(HnD.ctch.exp,"/",names(Tar)[i],"_grouped.yrs",sep=""),2000, 2400)
    # smart.par(n.plots=length(FINYrS.gped)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
    # fn.ctch.plot.grouped.yrs(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50,grouping=5)
    # mtext(Lat.exp,side=2,line=0.5,las=3,cex=1.1,outer=T)
    # mtext(Lon.exp,side=1,line=0.5,cex=1.1,outer=T)
    # dev.off()
    
    rm(ddd)
  }
  
  #2.2 total catch
  ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
               METHOD%in%c("GN","LL")& LAT<=(-26))
  FINYrS=unique(ddd$FINYEAR)
  fn.fig(paste(HnD.ctch.exp,paste('Heat_map','Total',sep='_'),sep="/"),2400, 2400)
  smart.par(n.plots=length(FINYrS)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
  fn.ctch.plot.all.yrs(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=numInt)
  mtext(Lat.exp,side=2,line=0.5,las=3,cex=1.1,outer=T)
  mtext(Lon.exp,side=1,line=0.5,cex=1.1,outer=T)
  dev.off()
  rm(ddd)
  
  
  #3. Movie
  do.movi=FALSE
  if(do.movi)
  {
    Frame.speed=.4  #match to talking speed
    ani.options(ani.width=480,ani.height=480)
    fn.ctch.plot.all.yrs.movie=function(DATA,tcl.1,tcl.2,numInt,SP) 
    {
      DATA$LIVEWT.c=DATA$LIVEWT.c/1000   #in tonnes
      DATA$LONG=as.numeric(DATA$LONG)
      DATA$blk=substr(DATA$BLOCKX,1,4)
      A=aggregate(LIVEWT.c~FINYEAR+blk,DATA,sum)
      Ymax=max(A$LIVEWT.c)
      Ymin=min(A$LIVEWT.c)
      Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
      #Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
      DATA$LAT=as.numeric(substr(DATA$LAT,1,3))
      DATA$LONG=as.numeric(substr(DATA$LONG,1,3))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      saveGIF({
        for(y in 1:length(FINYrS))
        {
          A=subset(DATA,FINYEAR==FINYrS[y])
          MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
          colnames(MapCatch)=c("BLOCKX.c","Total Catch")
          id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
          MapCatch$LAT=DATA$LAT[id]
          MapCatch$LONG=DATA$LONG[id]
          msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
          msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
          if(length(msn.lat)>0)
          {
            dummy=MapCatch[1:length(msn.lat),]
            dummy$`Total Catch`=0
            dummy$LAT=msn.lat
            dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
            MapCatch=rbind(MapCatch,dummy)
          }
          
          if(unique(min(A$YEAR.c))>2007)
          {
            MapCatch$`Total Catch`=ifelse(MapCatch$LAT>=(-32) & MapCatch$LAT<=(-31) & MapCatch$LONG<116,0,MapCatch$`Total Catch`)
          }
          MapCatch$LAT.cen=MapCatch$LAT-.5
          MapCatch$LONG.cen=MapCatch$LONG+.5  
          MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
          MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
          long=sort(unique(MapCatch$LONG.cen))
          lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
          MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
          Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
                                     timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
          Reshaped=Reshaped[order(Reshaped[,1]),]
          Reshaped=Reshaped[,-1]	
          numberLab=10
          colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
          
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          image(long,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,axes = FALSE,add=T)			
          axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
          axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
          par(new=T)
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          legend('top',FINYrS[y],bty='n',cex=1.2)
          axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
          axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
          mtext(Lat.exp,side=2,line=1,las=3,cex=1.75,outer=T)
          mtext(Lon.exp,side=1,line=1,cex=1.75,outer=T)
          color.legend(quantile(a,probs=.9),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.5),
                       paste(round(Breaks,0),"t"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.95)
        }
      },movie.name=paste(HnD.ctch.exp,paste('movie.',SP,".gif",sep=''),sep="/"),interval=Frame.speed,loop =1)
    }
    
    #3.1 by indicator species
    for(i in 1:length(Tar))
    {
      ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                   METHOD%in%c("GN","LL")& LAT<=(-26) & SPECIES==Tar[[i]])
      FINYrS=unique(ddd$FINYEAR)
      fn.ctch.plot.all.yrs.movie(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=20,SP=names(Tar)[i])
      rm(ddd)
    }
    
    #3.2 total catch
    ddd=subset(Data.monthly,Estuary=="NO" & !is.na(LIVEWT.c) &
                 METHOD%in%c("GN","LL")& LAT<=(-26))
    FINYrS=unique(ddd$FINYEAR)
    fn.ctch.plot.all.yrs.movie(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=20,SP="Total")
    rm(ddd)
    
  }
  
  
  #--Effort 
  
  #1. km gn days 
  
  #Monthly
  Attach.Effort.monthly.blk.c=aggregate(Km.Gillnet.Days.c~MONTH+BLOCKX+VESSEL+zone+FINYEAR+LAT+LONG,
                                        data=subset(Effort.monthly,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T)
  Attach.Effort.monthly.blk.c=fn.split.boundary(Attach.Effort.monthly.blk.c,"Km.Gillnet.Days.c")
  Attach.Effort.monthly.blk.c=aggregate(Km.Gillnet.Days.c~BLOCKX+FINYEAR+LAT+LONG,data=Attach.Effort.monthly.blk.c,sum,na.rm=T)
  
  #Daily
  if(Use.Date=="NO")
  {
    Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~ID+blockx+vessel+finyear+LAT+LONG,
                                        data=subset(Effort.daily,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T)
    
  }
  if(Use.Date=="YES")
  {
    Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~date+blockx+vessel+finyear+LAT+LONG,
                                        data=subset(Effort.daily,LAT<=(-26)&LAT>(-36.1)&LONG<=(129)&LONG>(111)),max,na.rm=T) 
    
  }
  Attach.Effort.daily.blk.c=aggregate(Km.Gillnet.Days.c~blockx+finyear+LAT+LONG,data=Attach.Effort.daily.blk.c,sum,na.rm=T)
  names(Attach.Effort.daily.blk.c)=names(Attach.Effort.monthly.blk.c)
  
  Effor.monthly.km.gn.d.block=rbind(Attach.Effort.monthly.blk.c,Attach.Effort.daily.blk.c)
  
  fn.effort.breaks=function(DD)
  {
    Max.eff=Yr.rango
    for (j in 1:length(Yr.rango))
    {
      DATA=subset(DD,FINYEAR%in%Yr.rango[[j]])
      dummy=with(DATA,aggregate(Km.Gillnet.Days.c,list(BLOCKX),FUN=sum,na.rm=T))
      Max.eff[[j]]=dummy$x
    }
    Max.eff=unlist(Max.eff)
    Breaks=quantile(Max.eff,probs=seq(0,1,1/numInt),na.rm=T)  #breaks for effort
    return(Breaks)
  }
  Yr.rango=vector('list',length(Yr.group))    
  for(e in 1:length(Yr.rango))
  {
    from=seq(Yr.group[e],Yr.group.plus[e])
    to=substr(seq(Yr.group[e]+1,Yr.group.plus[e]+1),3,4)
    Yr.rango[[e]]=paste(from,to,sep="-") 
  }
  Yr.rango=Yr.rango[which(unlist(lapply(Yr.rango,function(x) length(x)==5)))]
  EffortBreakS=fn.effort.breaks(Effor.monthly.km.gn.d.block)
  
  #get monthly effort
  Monthly=subset(Effort.monthly,NETLEN.c>100 & METHOD=="GN")
  Block.lat.long=Effort.monthly[!duplicated(Effort.monthly$BLOCKX),match(c("BLOCKX","LAT","LONG"),names(Effort.monthly))]
  Eff.monthly.c=aggregate(Km.Gillnet.Days.c~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)
  Eff.monthly.hour.c=aggregate(Km.Gillnet.Hours.c~VESSEL+BLOCKX+FINYEAR+MONTH+YEAR.c,data=Monthly,max,na.rm=T)
  THESE=c("BLOCKX","FINYEAR","MONTH","VESSEL","YEAR.c")
  Eff.monthly.c=merge(Eff.monthly.c,Eff.monthly.hour.c,by=THESE,all=T)
  Eff.monthly.c=merge(Eff.monthly.c,Block.lat.long,by=c("BLOCKX"),all.x=T)
  Eff.monthly.c$LAT=-abs(Eff.monthly.c$LAT)
  
  #get daily effort
  Daily=subset(Effort.daily,netlen.c>100 & method=="GN")   
  Block.lat.long=Effort.daily[!duplicated(Effort.daily$blockx),match(c("blockx","LAT","LONG"),names(Effort.daily))]
  if(!Use.Date=="YES") Eff.daily.c.daily=aggregate(Km.Gillnet.Days.c~Same.return.SNo+vessel+finyear+month+
                                                     blockx+block10,data=Daily,max,na.rm=T)    
  if(Use.Date=="YES") Eff.daily.c.daily=aggregate(Km.Gillnet.Days.c~date+Same.return.SNo+vessel+finyear+month+
                                                    blockx+block10,data=Daily,max,na.rm=T)    
  if(!Use.Date=="YES") Eff.daily.hour.c.daily=aggregate(Km.Gillnet.Hours.c~Same.return.SNo+vessel+finyear+month+
                                                          blockx+block10,data=Daily,max,na.rm=T)    
  if(Use.Date=="YES") Eff.daily.hour.c.daily=aggregate(Km.Gillnet.Hours.c~date+Same.return.SNo+vessel+finyear+month+
                                                         blockx+block10,data=Daily,max,na.rm=T)    
  Eff.daily.c.daily=merge(Eff.daily.c.daily,Eff.daily.hour.c.daily,
                          by=c("blockx","block10","finyear","month","vessel","Same.return.SNo"),all=T)
  Eff.daily.c.daily=merge(Eff.daily.c.daily,Block.lat.long,by=c("blockx"),all.x=T)  
  
  #get daily effort aggregated at monthly level
  if(Use.Date=="NO")  Eff.daily.c=aggregate(Km.Gillnet.Days.c~ID+vessel+finyear+month+blockx,data=Daily,max,na.rm=T)    
  if(Use.Date=="YES") Eff.daily.c=aggregate(Km.Gillnet.Days.c~date+vessel+finyear+month+blockx,data=Daily,max,na.rm=T) 
  Eff.daily.c=aggregate(Km.Gillnet.Days.c~finyear+month+vessel+blockx,data=Eff.daily.c,sum,na.rm=T)
  if(Use.Date=="NO") Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~ID+vessel+finyear+month+blockx,data=Daily,max,na.rm=T)
  if(Use.Date=="YES") Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~date+vessel+finyear+month+blockx,data=Daily,max,na.rm=T)
  Eff.daily.hour.c=aggregate(Km.Gillnet.Hours.c~finyear+month+vessel+blockx,data=Eff.daily.hour.c,sum,na.rm=T)
  Eff.daily.c=merge(Eff.daily.c,Eff.daily.hour.c,by=c("blockx","finyear","month","vessel"),all=T)
  Eff.daily.c=merge(Eff.daily.c,Block.lat.long,by=c("blockx"),all.x=T)
  
  
  #1.1 Each year 
  fn.eff.plot.all.yrs=function(DATA,WHAT)   #function plot effort
  {
    MapEffort=subset(DATA,!is.na(Km.Gillnet.Days.c))
    colnames(MapEffort)=c("FINYEAR","BLOCKX.c","LAT","LONG","Total effort")
    lat=sort(unique(MapEffort$LAT))
    # ids=match(-27,lat)
    # if(is.na(ids))      #add dummy for sorting image out of whack when missing Lat
    # {
    #   adD=MapEffort[1,]
    #   adD[,2]=NA
    #   adD$BLOCKX.c="-27 113"
    #   adD$LAT=-27
    #   MapEffort=rbind(MapEffort,adD)
    # }
    
    if(WHAT=="monthly")
    {
      MapEffort$LAT.cen=MapEffort$LAT-.5
      MapEffort$LONG.cen=MapEffort$LONG+.5 
    }else
    {
      MapEffort$LAT.cen=MapEffort$LAT
      MapEffort$LONG.cen=MapEffort$LONG 
      
    }
    
    MapEffort=MapEffort[order(MapEffort$LAT.cen,MapEffort$LONG.cen),]
    MapEffort=subset(MapEffort,LONG.cen<=South.WA.long[2])
    long=sort(unique(MapEffort$LONG.cen))
    lat=sort(unique(MapEffort$LAT.cen))      #latitude vector for image  
    
    MapEffort=MapEffort[,match(c("LONG.cen","LAT.cen","Total effort"),names(MapEffort))]  
    #MapEffort=rbind(MapEffort,dummy)#keep in perspective
    
    Reshaped=as.matrix(reshape(MapEffort,idvar="LONG.cen",    #transposed as matrix   
                               timevar="LAT.cen",v.names="Total effort", direction="wide"))  
    Reshaped=Reshaped[order(Reshaped[,1]),]
    Reshaped=Reshaped[,-1]  									
    
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    #plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    #image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=EffortBreaks,axes = FALSE,add=T)			
    #par(new=T)
    #plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
    image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,ylim=c(-36,-26),xlim=c(113,129),breaks=EffortBreaks,axes = FALSE)
    box()
    
  }
  HnD.eff.exp=handl_OneDrive("Analyses/Catch and effort/Outputs/Spatio.temporal_Effort")
  
  #Monthly
  s1=subset(Eff.monthly.c,LAT<=(-26) & !FINYEAR%in%FINYEAR.daily)  
  s1$BLOCKX=substr(s1$BLOCKX,1,4)
  s1=aggregate(Km.Gillnet.Days.c~FINYEAR+BLOCKX+LAT+LONG,s1,sum)
  Mn.yrs=sort(unique(s1$FINYEAR))  
  tiff(file=paste(HnD.eff.exp,"Gillnets_all.yrs_monthly.tiff",sep="/"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(Mn.yrs)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  EffortBreaks=quantile(s1$Km.Gillnet.Days.c,probs=seq(0,1,1/numInt),na.rm=T)
  for (i in 1:length(Mn.yrs))
  {
    print(paste('Gillnets_all.yrs_monthly for ----',Mn.yrs[i]))
    fn.eff.plot.all.yrs(DATA=subset(s1,FINYEAR== Mn.yrs[i]),WHAT="monthly")
    mtext(Mn.yrs[i],side=3,line=0,cex=.95)
    axis(side = 1, at =Long.seq, labels = F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
  }  
  plot(1:1,col="transparent",axes=F,ylab="",xlab="")
  color.legend(xl=0.92,yb=0.39,xr=1.29,yt=1.4,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
  mtext(Lat.exp,side=2,line=.5,las=3,cex=1.25,outer=T)
  mtext(Lon.exp,side=1,line=.8,cex=1.25,outer=T)
  dev.off()
  
  #Daily      
  s2=subset(Eff.daily.c.daily,LAT<=(-26))
  s2$LatDeg=with(s2,as.numeric(substr(block10,1,2)))  
  s2$LatMin=with(s2,10*as.numeric(substr(block10,3,3)))  
  s2$LongDeg=with(s2,100+as.numeric(substr(block10,4,5)))
  s2$LongMin=with(s2,10*as.numeric(substr(block10,6,6)))
  s2$LAT=-abs(with(s2,LatDeg+(LatMin/60)))
  s2$LONG=with(s2,LongDeg+(LongMin/60))
  s2=aggregate(Km.Gillnet.Days.c~finyear+block10+LAT+LONG,s2,sum)
  Dy.yrs=sort(unique(s2$finyear))
  tiff(file=paste(HnD.eff.exp,"Gillnets_all.yrs_daily.tiff",sep="/"),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(Dy.yrs)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.5, 0))
  for (i in 1:length(Dy.yrs))
  {
    print(paste('Gillnets_all.yrs_daily for ----',Dy.yrs[i]))
    
    EffortBreaks=quantile(s2$Km.Gillnet.Days.c,probs=seq(0,1,1/numInt),na.rm=T)
    fn.eff.plot.all.yrs(DATA=subset(s2,finyear== Dy.yrs[i]),WHAT="daily")
    axis(side = 1, at =Long.seq, labels =F, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = F,tcl = .35,las=2,cex.axis=1,hadj=.65)
    mtext(Dy.yrs[i],side=3,line=0,cex=.95)
    axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-0.9)
    axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=.65)
  }  
  plot(1:1,col="transparent",axes=F,ylab="",xlab="")
  color.legend(xl=0.92,yb=0.5,xr=1.29,yt=1.4,round(EffortBreaks,0),rect.col=couleurs,gradient="y",col=colLeg,cex=0.75)
  mtext(Lat.exp,side=2,line=.25,las=3,cex=1.25,outer=T)
  mtext(Lon.exp,side=1,line=.8,cex=1.25,outer=T)
  dev.off()
  
  #Monthly and daily in same plot
  fn.eff.plot.all.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt) 
  {
    DATA$blk=with(DATA,paste(LAT,LONG))
    A=aggregate(km.gillnet.days.c~finyear+blk,DATA,sum)
    Ymax=max(A$km.gillnet.days.c)
    Ymin=min(A$km.gillnet.days.c)
    Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
    #Breaks=seq(Ymin,Ymax,length.out=(numInt+1))
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
    colfunc <- colorRampPalette(c("yellow", "red"))
    couleurs=colfunc(numInt)
    
    numberLab=10
    colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
    
    for(y in 1:length(FINYrS))
    {
      A=subset(DATA,finyear==FINYrS[y])
      MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
      colnames(MapEffrt)=c("BLOCKX.c","Effort")
      id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
      MapEffrt$LAT=DATA$LAT[id]
      MapEffrt$LONG=DATA$LONG[id]
      MapEffrt$LAT.cen=MapEffrt$LAT-.5
      MapEffrt$LONG.cen=MapEffrt$LONG+.5  
      MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
      MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapEffrt$LONG.cen))
      lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
      MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
      Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="Effort", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      legend('top',FINYrS[y],bty='n',cex=1.2)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    plot(a,b,ann=F,axes=F,col='transparent')
    color.legend(quantile(a,probs=.8),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.25),
                 paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=.7)
  }
  s1=subset(Eff.monthly.c,LAT<=(-26) &FINYEAR%in% Mn.yrs)
  s1$BLOCKX=substr(s1$BLOCKX,1,4)
  names(s1) =  casefold(names(s1))
  s2=subset(Eff.daily.c,LAT<=(-26) &!finyear%in% Mn.yrs)
  names(s2) =  casefold(names(s2))
  ddd=rbind(s1[,match(names(s2),names(s1))],s2)
  FINYrS=sort(unique(ddd$finyear))
  ddd$LAT=as.numeric(substr(ddd$lat,1,3))
  ddd$LONG=as.numeric(substr(ddd$long,1,3))
  fn.fig(paste(HnD.eff.exp,"Gillnets_all.yrs",sep="/"),2400, 2400)
  smart.par(n.plots=length(FINYrS)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
  fn.eff.plot.all.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50)
  mtext(Lat.exp,side=2,line=0.5,las=3,cex=1.1,outer=T)
  mtext(Lon.exp,side=1,line=0.5,cex=1.1,outer=T)
  dev.off()
  
  
  #1.2 Grouped years 
  fn.eff.plot.grouped.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt,grouping) 
  {
    DATA$blk=with(DATA,paste(LAT,LONG))
    a=South.WA.long[1]:South.WA.long[2]
    b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
    DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
    AA=vector('list',length(FINYrS.gped))
    for(y in 1:length(FINYrS.gped))
    {
      A=aggregate(km.gillnet.days.c~blk,subset(DATA,finyear%in%FINYrS.gped[[y]]),sum,na.rm=T)
      A$y.group=y
      AA[[y]]=A
    }
    A=do.call(rbind,AA)
    Ymax=max(A$km.gillnet.days.c)
    Ymin=min(A$km.gillnet.days.c)
    Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
    couleurs=rev(heat.colors(numInt))
    numberLab=10
    colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
    
    for(y in 1:length(FINYrS.gped))
    {
      A=subset(DATA,finyear%in%FINYrS.gped[[y]])
      MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
      colnames(MapEffrt)=c("BLOCKX.c","Effort")
      id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
      MapEffrt$LAT=DATA$LAT[id]
      MapEffrt$LONG=DATA$LONG[id]
      MapEffrt$LAT.cen=MapEffrt$LAT-.5
      MapEffrt$LONG.cen=MapEffrt$LONG+.5  
      MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
      MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
      long=sort(unique(MapEffrt$LONG.cen))
      lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
      MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
      Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                 timevar="LAT.cen",v.names="Effort", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      Reshaped=Reshaped[,-1]	
      
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
      axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
      axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
      par(new=T)
      plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
      text(115.60,-31.96,"Perth",pos=4,cex=1.2)
      legend('top',names(FINYrS.gped)[y],bty='n',cex=1.2)
      axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
      axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
    }
    plot(a,b,ann=F,axes=F,col='transparent')
    color.legend(quantile(a,probs=.8),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.25),
                 paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=.7)
  }
  grouping=5
  FINYrS.gp=seq(1,length(FINYrS),by=grouping)
  FINYrS.gped=vector('list',length(FINYrS.gp))
  for(f in 1:length(FINYrS.gped))
  {
    if(f==length(FINYrS.gped))
    {
      FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
      if(length(FINYrS.gped[[f]])==1)
      {
        names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
      }else
      {
        names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
      }
      
    }else
    {
      FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
      names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
    }
  }
  FINYrS.gped=FINYrS.gped[which(unlist(lapply(FINYrS.gped,function(x) length(x)==5)))]
  fn.fig(paste(HnD.eff.exp,"Gillnets_grouped.yrs",sep="/"),2000, 2400)
  smart.par(n.plots=length(FINYrS.gped)+1,MAR=c(1,1,1,1),OMA=c(2,2,.1,.1),MGP=c(1, 0.35, 0))
  fn.eff.plot.grouped.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50,grouping=grouping)
  mtext(Lat.exp,side=2,line=0.5,las=3,cex=1.1,outer=T)
  mtext(Lon.exp,side=1,line=0.5,cex=1.1,outer=T)
  dev.off()
  
  
  #1.3 Movie (Monthly and Daily in same plot)
  if(do.movi)
  {
    ani.options(ani.width=480,ani.height=480)
    movie.fn.eff.plot.all.yrs.mon.and.daily=function(DATA,tcl.1,tcl.2,numInt) 
    {
      DATA$blk=with(DATA,paste(LAT,LONG))
      A=aggregate(km.gillnet.days.c~finyear+blk,DATA,sum)
      Ymax=max(A$km.gillnet.days.c)
      Ymin=min(A$km.gillnet.days.c)
      Breaks=quantile(A$km.gillnet.days.c,probs=seq(0,1,1/numInt),na.rm=T)
      #Breaks=seq(Ymin,Ymax,length.out=(numInt+1))
      a=South.WA.long[1]:South.WA.long[2]
      b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
      DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
      colfunc <- colorRampPalette(c("yellow", "red"))
      couleurs=colfunc(numInt)
      
      numberLab=10
      colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
      saveGIF({
        for(y in 1:length(FINYrS))
        {
          par(las=1,mar=c(1,1,.1,1),oma=c(3,4,.1,.1),mgp=c(3.5,.5,0))
          A=subset(DATA,finyear==FINYrS[y])
          MapEffrt=with(A,aggregate(km.gillnet.days.c,list(BLOCKX.c),FUN=sum,na.rm=T))
          colnames(MapEffrt)=c("BLOCKX.c","Effort")
          id=unique(match(MapEffrt$BLOCKX.c,DATA$BLOCKX.c))
          MapEffrt$LAT=DATA$LAT[id]
          MapEffrt$LONG=DATA$LONG[id]
          MapEffrt$LAT.cen=MapEffrt$LAT-.5
          MapEffrt$LONG.cen=MapEffrt$LONG+.5  
          MapEffrt=MapEffrt[order(MapEffrt$LAT.cen,MapEffrt$LONG.cen),]
          MapEffrt=subset(MapEffrt,LONG.cen<=South.WA.long[2])
          long=sort(unique(MapEffrt$LONG.cen))
          lat=sort(unique(MapEffrt$LAT.cen))      #latitude vector for image  
          MapEffrt=MapEffrt[,match(c("LONG.cen","LAT.cen","Effort"),names(MapEffrt))]  
          Reshaped=as.matrix(reshape(MapEffrt,idvar="LONG.cen",  	#transposed as matrix 	
                                     timevar="LAT.cen",v.names="Effort", direction="wide"))	
          Reshaped=Reshaped[order(Reshaped[,1]),]
          Reshaped=Reshaped[,-1]	
          
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          image(long,lat,z=Reshaped,xlab="",ylab="",col =couleurs,breaks=Breaks,axes = FALSE,add=T)			
          axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
          axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
          par(new=T)
          plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
          legend('top',FINYrS[y],bty='n',cex=1.75)
          axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
          axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
          mtext(Lat.exp,side=2,line=1,las=3,cex=1.75,outer=T)
          mtext(Lon.exp,side=1,line=1,cex=1.75,outer=T)
          color.legend(quantile(a,probs=.9),quantile(b,probs=.91),quantile(a,probs=.95),quantile(b,probs=.5),
                       paste(round(Breaks,0),"km gn d"),rect.col=couleurs,gradient="y",col=colLeg,cex=1)
          
        }
      },movie.name=paste(HnD.eff.exp,"movie.Effort.gif",sep="/"),interval=Frame.speed,loop =1)
      
      
    }
    s1=subset(Eff.monthly.c,LAT<=(-26) &FINYEAR%in% Mn.yrs)
    s1$BLOCKX=substr(s1$BLOCKX,1,4)
    names(s1) =  casefold(names(s1))
    s2=subset(Eff.daily.c,LAT<=(-26) &!finyear%in% Mn.yrs)
    names(s2) =  casefold(names(s2))
    ddd=rbind(s1[,match(names(s2),names(s1))],s2)
    FINYrS=sort(unique(ddd$finyear))
    ddd$LAT=as.numeric(substr(ddd$lat,1,3))
    ddd$LONG=as.numeric(substr(ddd$long,1,3))
    movie.fn.eff.plot.all.yrs.mon.and.daily(DATA=ddd,tcl.1=.1,tcl.2=.1,numInt=50)
    
  }
  
  
