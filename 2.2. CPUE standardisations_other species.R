#Pseudo code to catch and effort standardisation other species

#Get species used in Martel & Froese method

#CPUE.stand=function(RANGO (use Last and Stevens), SP)
#{
  # data=subset(Data.monthly.GN/Data,daily.GN, SPECIES=SP
  #             & LAT and LONG within RANGO )
  # aggregate (catch, effort by record)
  # reshape(data) to have 1 row per record
  # apply delta method (Catch~year+month+block, etc offset(effort))
  # GET CI
  # plot
#}