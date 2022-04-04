library(tidyverse)
library(mgcv)
library(emmeans)  #for model predictions
library(flextable)
# Data section  ---------------------------------------------------------
#1. Sharks data base
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
User="Matias"
if(User=="Matias") source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R'))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Nominal_cpue_functions.R"))

# Manipulation  ---------------------------------------------------------
Res.ves=c("HAM","HOU","NAT","FLIN","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")

Dusky.nursery=DATA%>%
              filter(SPECIES=='BW' &Mid.Lat<(-27))%>%
              left_join(DATA.bio%>%
                          filter(SPECIES=='BW' & UMBIL_SCAR%in%c('P','Y'))%>%
                          dplyr::select(UNIQUE_ID,SPECIES,UMBIL_SCAR),
                        by=c('UNIQUE_ID','SPECIES'))%>%
              dplyr::select(SHEET_NO,Mid.Lat,Mid.Long,BLOCK,Month,year,BOAT,SKIPPER,
                            SPECIES,Numbers,TL,FL,Method,SOAK.TIME,NET_LENGTH,BOTDEPTH,MESH_SIZE,
                            Set.time,Haul.time,COMMON_NAME,CAES_Code,CAAB_code,UMBIL_SCAR)%>%
  mutate(FL=ifelse(is.na(FL) & SPECIES=='BW',(TL-1.9133)/1.2062,FL),
         Neonate=ifelse(UMBIL_SCAR%in%c('P','Y'),'Y',
                 ifelse(is.na(UMBIL_SCAR) & FL<85,'Y',
                 'N')))%>%
  filter(!is.na(FL))%>%
  mutate(Neonate=ifelse(FL>100,'N',Neonate))


Dusky=DATA%>%
      filter(Method=='GN' & !is.na(NET_LENGTH) & !is.na(SOAK.TIME) & Mid.Lat<(-27) & !BOAT%in%Res.ves)%>%
      left_join(DATA.bio%>%
                  filter(SPECIES=='BW' & UMBIL_SCAR%in%c('P','Y'))%>%
                  dplyr::select(UNIQUE_ID,SPECIES,UMBIL_SCAR),
                by=c('UNIQUE_ID','SPECIES'))%>%
      dplyr::select(SHEET_NO,Mid.Lat,Mid.Long,BLOCK,Month,year,BOAT,SKIPPER,
                    SPECIES,Numbers,TL,FL,Method,SOAK.TIME,NET_LENGTH,BOTDEPTH,MESH_SIZE,
                    Set.time,Haul.time,COMMON_NAME,CAES_Code,CAAB_code,UMBIL_SCAR)
  
Dusky=Dusky%>%
  mutate(FL=ifelse(is.na(FL) & SPECIES=='BW',(TL-1.9133)/1.2062,FL),
         Neonate=ifelse(UMBIL_SCAR%in%c('P','Y'),'Y',
                  ifelse(is.na(UMBIL_SCAR) & FL<85,'Y',
                         'N')),
         SP=ifelse(SPECIES=="BW" & Neonate=='Y',"BW","Other"))%>%
  filter(!is.na(FL))

Dusky%>%
  filter(SPECIES=='BW')%>%
  ggplot(aes(FL,colour=Neonate))+
  geom_density()

Dusky%>%
  filter(Mid.Lat>-40 & Neonate=='Y')%>%
  ggplot(aes(Mid.Long,Mid.Lat))+geom_point(alpha = 0.2,color="red") +
  geom_density_2d(colour="black",bins=50)

Dusky.h=Dusky%>%
          group_by(SHEET_NO,Mid.Lat,Mid.Long,BLOCK,Month,year,BOAT,SKIPPER,
                   SP,SOAK.TIME,NET_LENGTH,BOTDEPTH,MESH_SIZE,Set.time,Haul.time)%>%
          summarise(N=sum(Numbers))%>%
          ungroup()%>%
          spread(SP,N,fill=0)%>%
          mutate(Effort=SOAK.TIME*NET_LENGTH,  #effort in km hm hour
                 log.Effort=log(Effort),
                 YEAR=factor(year),
                 BLOCK=factor(BLOCK),
                 BOAT=factor(BOAT),
                 cpue=BW/Effort)

Dusky.h%>%
  ggplot(aes(YEAR,cpue))+
  geom_boxplot()+ylim(0,5)




# Standardisation  ---------------------------------------------------------
formula.gam=formula(BW ~ YEAR + BLOCK + s(BOAT,bs="re") + s(Month, bs = "cc") +s(BOTDEPTH)+ offset(log.Effort))

fn.stand=function(d,FORMULA)
{
  pois=gam(FORMULA,family=poisson(),data=d,method="REML")
  NB=gam(FORMULA,family=nb(),data=d,method="REML")
  return(list(data=d, pois=pois,NB=NB))
}
system.time({ MOD=fn.stand(d=Dusky.h%>%
                             filter(Mid.Long<119 & Mid.Lat<(-30)),
                           FORMULA=formula.gam)})
AIC.tab=data.frame(Pois=MOD$pois$aic,
                   NB=MOD$NB$aic)
AIC.tab$Best.mod=colnames(AIC.tab)[apply(AIC.tab,1,which.min)]

Best=MOD[[AIC.tab$Best.mod]]
as_flextable(Best)

#Predict annual index
d=MOD$data
Index=summary(emmeans(Best, 'YEAR', type="response"))  

  #relative index
Relative.index=Index%>%
                mutate(yr=YEAR,
                       CV=SE/response,
                       UppCI=upper.CL/mean(response),
                       LowCI=lower.CL/mean(response),
                       SE=SE/mean(response),
                       MeAn=response/mean(response))%>%
                dplyr::select(yr,MeAn,UppCI,LowCI,CV,SE)

#Plot index
Nominal.index=Dusky.h%>%
              group_by(year)%>%
              summarise(Mean=mean(cpue),
                        SD=sd(cpue),
                        n=n(),
                        SE=SD/sqrt(n))%>%
              ungroup()%>%
              mutate(Rel.Mean=Mean/mean(Mean),
                     Rel.SE=SE/mean(Mean),
                     UpperCI=Rel.Mean+1.96*Rel.SE,
                     LowerCI=Rel.Mean-1.96*Rel.SE)

Relative.index%>%
  mutate(year=as.numeric(as.character(yr)))%>%
  ggplot(aes(year,MeAn))+
  geom_point()+
  geom_errorbar(aes(ymin=LowCI, ymax=UppCI), width=.2,
                position=position_dodge(.9))+
  geom_point(data=Nominal.index,aes(year+.25,Rel.Mean), inherit.aes = FALSE,color='orange',pch=19)+
  geom_errorbar(data=Nominal.index,aes(x=year+.25,y=Rel.Mean,ymin=LowerCI, ymax=UpperCI), width=.2,
                position=position_dodge(.9),color='orange',pch=19)


#Export index
write.csv(Relative.index,
          handl_OneDrive("Analyses/Data_outs/Dusky shark/Dusky shark.neonate.index.csv"),
          row.names = F)
# Explore 2006  ---------------------------------------------------------
Dusky.h%>%
  mutate(Dummy=ifelse(BW>0,1,0))%>%
  group_by(year)%>%
  summarise(PresAbs=sum(Dummy),
            Pos=sum(BW))

Dusky.h%>%
  filter(BW>0)%>%
  group_by(year,BOAT)%>%
  tally()%>%
  spread(year,n)%>%
  data.frame

Dusky.h%>%
  filter(BW>0 & BOAT=='E35')%>%
  group_by(year,Month)%>%
  tally()%>%
  spread(year,n)%>%
  data.frame

Dusky.h%>%
  filter(Mid.Lat>-37)%>%
  ggplot(aes(Mid.Long, Mid.Lat,colour=BOAT))+
  geom_point()+
  facet_wrap(~year)


# Nurseries  ---------------------------------------------------------
#Notes:
#Do All species. The Point of analysis is location of 0+.
# Class 0+ as 'YES umb scar' or within birth range. Class non 0+. Add to Dani's thesis


# Do zero inflated gam with 0 catches. Response var is number of neonates and predictors are method, month, lat-long interaction. Combine years. Produce heat maps. Could try with other species?
#   Nursery criteria:
#   
#   1) sharks are more commonly encountered in the area than other areas; 
#   (2) sharks have a tendency to remain or return for extended periods; 
#   (3) the area or habitat is repeatedly used across years.
# 
# For Criteria 3, i need the 1, 2 year olds, as done for snapper


Dusky.nursery%>%
  filter(Neonate=='Y')%>%
  mutate(FL.bin=10*round(FL/10),
         Mid.Long.jit=jitter(Mid.Long),
         Mid.Lat.jit=jitter(Mid.Lat))%>%
  ggplot(aes(Mid.Long.jit,Mid.Lat.jit))+
  geom_point(aes(fill=FL.bin),size=1.5,shape=21,alpha=.2)+
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(x=Mid.Long.jit, y=Mid.Lat.jit, fill = after_stat(level)), colour = "black",
                  bins = 20,alpha = 0.7)+
  scale_fill_distiller(palette = "Blues", direction = 1)+
  #facet_wrap(~Month)+
  theme(legend.position = 'none',
        plot.title = element_text(size=22),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        plot.margin=unit(c(0,0.1,0,0),"cm"))