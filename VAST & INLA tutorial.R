#devtools::install_github("James-Thorson-NOAA/VAST")

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#inla.upgrade()

library(VAST)
library(INLA)
library(tictoc)
library(simglm)
library(tidyverse)
library(mgcv)

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
setwd(handl_OneDrive('Analyses/Catch and effort/VAST & INLA'))

# Index standardization-------------------------------------------------------------------------
#Load data
example = load_example( data_set="EBS_pollock" )

#Configure model settings
settings = make_settings( n_x = 100, Region = example$Region, purpose = "index2" )
#settings = make_settings( n_x = 100, Region = 'Other', purpose = "index2" )

#Run the model time: 492.71 sec elapsed
tic()
fit = fit_model( settings = settings, 
                 Lat_i = example$sampling_data[,'Lat'], 
                 Lon_i = example$sampling_data[,'Lon'], 
                 t_i = example$sampling_data[,'Year'], 
                 b_i = example$sampling_data[,'Catch_KG'], 
                 a_i = example$sampling_data[,'AreaSwept_km2'] )
toc()

#Summarise model
summary(fit)

#Plot the model
plot( fit )

# Simulated data ----------------------------------------------------------
sim_arguments <- list(
  formula = y ~ 1 + weight + age + sex,
  fixed = list(weight = list(var_type = 'continuous', mean = 180, sd = 30),
               age = list(var_type = 'ordinal', levels = 30:60),
               sex = list(var_type = 'factor', levels = c('male', 'female'))),
  error = list(variance = 1),
  sample_size = 100,
  reg_weights = c(0, 0.5, -1, 0.5)
)

Dat=simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments)

Dat%>%
  ggplot(aes(weight,y,color=age))+
  geom_point()+
  facet_wrap(~sex)

mod=glm(y ~ weight + age + sex,data = Dat)
summary(mod)


#random effect
sim_arguments <- list(
  formula = y ~ 1 + weight + age + sex + (1 | neighborhood),
  reg_weights = c(4, -0.03, 0.2, 0.33),
  fixed = list(weight = list(var_type = 'continuous', mean = 180, sd = 30),
               age = list(var_type = 'ordinal', levels = 30:60),
               sex = list(var_type = 'factor', levels = c('male', 'female'))),
  error = list(variance = 1),
  randomeffect = list(int_neighborhood = list(variance = 8, var_level = 2)),
  sample_size = list(level1 = 1000, level2 = 20)
)

nested_data <- sim_arguments %>%
  simulate_fixed(data = NULL, .) %>%
  simulate_randomeffect(sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments)

b <- gam(y ~ weight + age + sex +s(neighborhood,bs="re"),data=nested_data,method="REML")
summary(b)         


#
sim_arguments <- list(
  formula = y ~ 1 + Lat+ Long + age + sex + (1 | neighborhood),
  reg_weights = c(4, -0.03, 0.2, 0.33),
  fixed = list(weight = list(var_type = 'continuous', mean = 180, sd = 30),
               age = list(var_type = 'ordinal', levels = 30:60),
               sex = list(var_type = 'factor', levels = c('male', 'female'))),
  error = list(variance = 1),
  randomeffect = list(int_neighborhood = list(variance = 8, var_level = 2)),
  sample_size = list(level1 = 1000, level2 = 20)
)

dat <- sim_arguments %>%
  simulate_fixed(data = NULL, .) %>%
  simulate_randomeffect(sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments)

Dat=expand.grid(Year=rep(1982:2014,each=10),
                AreaSwept_km2=0.01,
                Lat=seq(54,62,by=5),
                Lon=-seq(158,178,by=10))
Dat=Dat%>%
        mutate(Lat=Lat*runif(nrow(Dat),0.99,1.01),
               Lon=Lon*runif(nrow(Dat),0.99,1.01),
               Catch_KG=(25-0.01*Year+1*AreaSwept_km2+0.05*Lat+-0.01*Lon)*runif(nrow(Dat),0.99,1.01))

Dat=Dat[-sample(1:nrow(Dat),round(0.5*nrow(Dat)),replace=F),]

Dat%>%
  ggplot(aes(Year,Catch_KG/AreaSwept_km2))+
  geom_point(aes(color=Lon))

mod=glm(Catch_KG ~ Year + AreaSwept_km2 + Lat + Lon,data = Dat)
summary(mod)



#Configure model settings
settings = make_settings( n_x = 100, Region = example$Region, purpose = "index2" )

#Run the model time: 
tic()
fit = fit_model( settings = settings, 
                 Lat_i = Dat$Lat, 
                 Lon_i = Dat$Lon, 
                 t_i = Dat$Year, 
                 b_i = Dat$Catch_KG, 
                 a_i = Dat$AreaSwept_km2 )
toc()

#Summarise model
summary(fit)

