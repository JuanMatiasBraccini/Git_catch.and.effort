
options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 

setwd('C:/Matias/Analyses/Catch and effort/Data_outs')

#vip: use WHOLE catch and effort data set (south and north) to show entire distribution across WA (though may have issues
#of spatial resolution if north doesn't have block10....)

#SPIN: "Varying degrees of spatial overlap between shark abundance hotspots 
#      and fishing closures in south-western Australia (use package VAST)"
# need to consider overlap (same approach as recent Nature movement paper and longline fishing)

#-Daily catch and effort for all methods, north of 26 S
Data.daily.NSF=read.csv('Data.daily.NSF.csv')
Effort.daily.NSF=read.csv('Effort.daily.NSF.csv')

#-Daily catch and effort for all methods, south of 26 S 
Data.daily=read.csv('Data.daily.csv')
Effort.daily=read.csv('Effort.daily.csv')

#Use function fn.compare.glm.gam() from 2.cpue standardisation.R
#predict only for block10 with data
#compare predictions from glm with block10 and gam as gam smoothes over the spatial heterogeneity of 
# the block10s, e.g blocks with similar lat and long may have different habitats
#use the predict function rather than lsmeans as I only care for the mean value, not the CI
# consider that seasonal patterns affect abundance (e.g. sandbars not around south in winter) so show
#different maps for different times of year
#maps are in relative terms (by setting year, month, etc to fix levels)
#use Delta method for probabilty and for positive catches, no need to combine as I'm not presenting index of abundance,
#just spatial hotspots
#need to convert longline to gillnet equivalent to have same effort measures
#Only focus on DAILY data, ditch monthly because it doesn't have the spatial resolution
#Formulas: use same formulas used for target and other species: Best.Model.daily.gam

#hotspots defined here as those areas with
#≥75th percentile of weighted daily location density;
