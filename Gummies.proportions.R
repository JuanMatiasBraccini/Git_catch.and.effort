
library(RODBC)
require(lubridate)  #for dates manipulation

#Species names
SPECIES.names=read.csv("C:/Matias/Data/Species.code.csv")

#Sharks data base
setwd("M:/Fisheries Research/Production Databases/Shark")  # working directory
channel <- odbcConnectAccess("Sharks.mdb")      
Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F)   
Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)   
close(channel)

DATA=merge(Boat_bio,Boat_hdr,by="SHEET_NO",all.x=T)
drop=c("UMBIL_SCAR","NO DISCARDS","CLASPLENTH","CLASP_CALC","GON_STAGE","RUN_SPERM","MAXOVRYDIA",
       "NO_YOLKOVA","UTERINESTG","NO_EMBRYOS","NO_UNDEVELOPED","EMBLEN_1","EMBLEN_2","EMBLEN_3","EMBLEN_4",
       "EMBLEN_5","EMBLEN_6","EMBLEN_7","EMBLEN_8","EMBLEN_9","EMBLEN_10","EMBLEN_11","EMBLEN_12",
       "STMCH_FULL","STMCH_CONT","BAG NO","VERT_SAMPL","Count","BOTTYPE","CTCHABLITY","RECORDER","SUBBLOCK","Buffer Zone",
       "MESH_DROP")
DATA=DATA[,-match(drop,names(DATA))]

names(DATA)[match(c("MID LAT","MID LONG"),names(DATA))]=c('Mid.Lat','Mid.Long')
DATA$Mid.Lat=-DATA$Mid.Lat
DATA$END1LATD=-DATA$END1LATD
DATA$END2LATD=-DATA$END2LATD


#Extract month,year, day, hour
DATA$date=as.Date(DATA$DATE,format="%Y-%m-%d")
DATA$Day=mday(DATA$date)
DATA$Month=month(DATA$date)
DATA$year=year(DATA$date)


#fix dodgy longitude                                      
DATA$END1LNGD=with(DATA,ifelse(SHEET_NO=="R00890",113.96,END1LNGD))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="R00890",113.9632,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00099",113.4272,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00597",113.2404,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00558",113.2055,Mid.Long))
#fix dodgy latitude
DATA$Mid.Lat=with(DATA,ifelse(SHEET_NO=="N00401",-20.90,Mid.Lat))


DATA$zone=as.character(with(DATA,ifelse(Mid.Long>=116.5 & Mid.Lat<=(-26),"Zone2",
                                        ifelse(Mid.Long<116.5 & Mid.Lat<=(-33),"Zone1",
                                               ifelse(Mid.Lat>(-33) & Mid.Lat<=(-26) & Mid.Long<116.5,"West",
                                                      ifelse(Mid.Lat>(-26) & Mid.Long<114,"Closed",
                                                             ifelse(Mid.Lat>(-26) & Mid.Long>=114 & Mid.Long<123.75,"North",
                                                                    ifelse(Mid.Lat>(-26) & Mid.Long>=123.75,"Joint",NA))))))))


DATA=merge(DATA,SPECIES.names,by.x="SPECIES",by.y="Species",all.x=T)


#Look at gummy sharks proportions
DATA=subset(DATA,SPECIES%in%c("GM","WG","GG"))



DATA$Bioregion=as.character(with(DATA,ifelse(Mid.Long>=115.5 & Mid.Long<=129 & Mid.Lat<=(-26),"SC", 
                     ifelse(Mid.Long<115.5 & Mid.Lat<=(-27),"WC",
                    ifelse(Mid.Long<=114.834 & Mid.Lat>(-27),"Gascoyne",
                   ifelse(Mid.Long>114.834 & Mid.Lat>=(-27) & Mid.Long<=129,"NC",NA))))))

DATA$Number=1
ag=aggregate(Number~SPECIES+Bioregion,DATA,sum)
ag1=aggregate(Number~Bioregion,DATA,sum)

Prop=merge(ag,ag1,by="Bioregion")
Prop$Prop=round(Prop$Number.x/Prop$Number.y,2)


write.csv(Prop,"C:/Matias/Analyses/Catch and effort/Gummies.prop.csv",row.names=F)
