#Compare Rory's effort with mine

test=aggregate(Km.Gillnet.Days.inv~FINYEAR+Same.return,Data.monthly,max,na.rm=T)
test1=aggregate(Km.Gillnet.Days.inv~FINYEAR,test,sum,na.rm=T)
plot(1:nrow(test1),test1[,2]/1e3,col=2)

lines(Results.pre.2013$FINYEAR,Results.pre.2013$TDGDLF.km.gn.days/1000)

test=aggregate(Km.Gillnet.Days.c~FINYEAR+Same.return,Data.monthly,max,na.rm=T)
test1=aggregate(Km.Gillnet.Days.c~FINYEAR,test,sum,na.rm=T)
lines(1:nrow(test1),test1[,2]/1e3,col=4)


#Possible differences in the way monthly effort is aggregated? My case: by ~Same.return, max function?

