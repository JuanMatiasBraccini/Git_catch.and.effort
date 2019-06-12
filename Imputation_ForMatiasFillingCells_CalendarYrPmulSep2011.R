# Following are manipulations of output from Mick ONeil Genstat CPUE
# standardisations 23/May/2011.

# The function 'Impute.Block' imputes catches for blocks for missing years
# in the series, following C. Walters algorithm in Carruthers et al. 2011.

setwd(
'C:/Matias/Analyses/Catch and effort/other scripts/Ross_Marriot_Imputation')

Blocks <- as.vector(read.csv('BlockNames.csv', header = F))
Years <- as.vector(c(seq(1993,2008,by=1)))
m1 <- as.matrix(read.csv('FullPrediction.csv', header = F))
m2 <- as.matrix(read.csv('PresentPrediction.csv', header = F))
nums <- as.matrix(read.csv('nPerCell.csv', header = F))

colnames(m1) <- t(Blocks); rownames(m1) <- Years
colnames(m2) <- t(Blocks); rownames(m2) <- Years
colnames(nums) <- t(Blocks); rownames(nums) <- Years

# Note: Upon revisiting this 16-Sep-11, m2 (with the zeros, from Mick O'Neil
#       Genstat program) is the one we want to impute.

imputed <- m2

# Have a look at this for other blocks: Need to revise code to account for 
# multiple imputation periods:

colnames(imputed)
# [1] "11230" "11240" "11250" "11260" "12230" "12240" "12250" "12260" "12280" "12290" "13220"
#[12] "13230" "13240" "13250" "13260" "14210" "14220" "14230" "14240" "15200" "15210" "15220"
#[23] "15230" "16200" "16210" "16220" "17200" "17210" "17220" "18200" "18210" "18220" "19200"
#[34] "19210"

# Matias: this example looks only at cases where there are at least 3 records of cpue (for
# simplicity).  You should get the idea from this.

# I can also provide the code for n<=3 but it is more convoluted.

N.Yrs.Per.Block <- numeric(ncol(imputed))
Pres.Absence <- matrix(,nrow=nrow(imputed),ncol=ncol(imputed))
colnames(Pres.Absence) <- colnames(imputed)

for (i in 1:nrow(imputed)){
  for (j in 1:ncol(imputed)){
    Pres.Absence[i,j] <- ifelse(imputed[i,j]>0,1,0)
    }
  }
N.Yrs.Per.Block <- colSums(Pres.Absence)
To.Include <- N.Yrs.Per.Block[N.Yrs.Per.Block>2]
imputed <- imputed[,colnames(imputed) %in% names(To.Include)==TRUE]

##**********************************************************
##**********************************************************

Impute.Block <- function(dataset,blockname){
#This labels the different year types for imputation: "Before","Interim","After"
Imputing.Block <- function(dataset,blockname){ # will need to input blockname in quotes
dataset <- as.matrix(dataset)
Subset <- as.vector(imputed[,colnames(dataset)==blockname])
Subset <- data.frame("Year"=rownames(dataset),"Catch"=Subset)
Subset$To.impute <- ifelse(Subset$Catch==0,"TRUE","FALSE")
nonzeroC <- Subset[Subset$Catch>0,]
Subset$LUMatcher <- charmatch(Subset$Year,nonzeroC$Year)
Subset$Year <- as.numeric(as.character(Subset$Year))
nonzeroC$Year <- as.numeric(as.character(nonzeroC$Year))
Subset$LUMatcher <- ifelse(
              Subset$Year<min(nonzeroC$Year),"Before",
              (ifelse(Subset$Year>max(nonzeroC$Year),"After",
              (ifelse(Subset$To.impute=="TRUE","Interim",
              Subset$LUMatcher)))))
Subset$Check.Interim <- as.numeric(as.character(ifelse(Subset$LUMatcher=="Interim",1,0)))
#Subset$Check.ns <- as.numeric(as.character(ifelse(Subset$Catch>0,1,0)))
return(Subset)
}

# This separates different periods of Interim years for imputation as groups:                      
Impute.Diffnt.Interims <- function(Subset){
if(sum(Subset$Check.Interim)==1){
  Subset$Interim <- ifelse(Subset$LUMatcher != "Interim", 0, "Interim")
} else {
  Interim <- Subset[Subset$LUMatcher=="Interim",]
  for(i in 1:(length(Interim$Year)-1)){Interim$Temp[1]<-1;Interim$Temp[i+1]<-Interim$Year[i+1]-Interim$Year[i]}
  Interim$NewGroup <- ifelse(Interim$Temp!=1,"Yes","No");Interim$NewGroup[1]<-"Yes"
  for(i in 1:length(Interim$Year)){
  Interim$Group[i] <- ifelse(Interim$NewGroup[i]=="Yes",
                              Interim$Year[i],
                              Interim$Group[i-1])
                              }
  Interim$Interim<-c(paste(Interim$LUMatcher,Interim$Group,sep=""))
  Subset$InterimRelevant <- Subset$Year %in% Interim$Year
  Interim$Year <- as.character(Interim$Year)
  Subset$Year <- as.character(Subset$Year)
  Subset$MatchInterim <- charmatch(Subset$Year,Interim$Year)
  Subset$MatchInterim[is.element(Subset$MatchInterim,c(NA))]<- 0
  for (i in (1:length(Subset$Year))){
  Subset$Interim[i] <- ifelse(Subset$InterimRelevant[i]==TRUE,Interim$Interim[Subset$MatchInterim[i]],0)}
}

return(Subset)
}

Impute.Catch.No.Interim <- function(Subset){
Subset$Catch <- ifelse(Subset$LUMatcher=="Before",
                          mean(Subset[Subset$Catch>0,]$Catch[1:3]),
                (ifelse(Subset$LUMatcher=="After",
                          Subset[Subset$Catch>0,]$Catch[length(Subset[Subset$Catch>0,]$Catch)],
                          Subset$Catch))) 
return(Subset)
}                          
                          
Impute.Catch <- function(Subset){
Subset$Catch <- ifelse(Subset$LUMatcher=="Before",
                          mean(Subset[Subset$Catch>0,]$Catch[1:3]),
                      (ifelse(Subset$LUMatcher=="After",
                          Subset[Subset$Catch>0,]$Catch[length(Subset[Subset$Catch>0,]$Catch)],
                          Subset$Catch)))                          
Subset.2 <- as.matrix(aggregate(as.numeric(as.character(Subset$Year)), list(Subset$Interim), summary))[,c(1:2,7)]
colnames(Subset.2) <- c("Group","Start","End"); Subset.2 <- as.data.frame(Subset.2)
Subset.3 <- Subset.2[Subset.2$Group!="0",]
No.Interim.Imputes <- length(Subset.3$Start)
Subset.3$Start <- as.numeric(as.character(Subset.3$Start))
Subset.3$End <- as.numeric(as.character(Subset.3$End))
for (i in 1:No.Interim.Imputes){
Subset.3$Imputed.Catch[i] <- mean(c(Subset[Subset$Year==(Subset.3$Start[i]-1),]$Catch,
                              Subset[Subset$Year==(Subset.3$End[i]+1),]$Catch))}
Subset$InterimGpRelevant <- Subset$Interim %in% Subset.3$Group
Subset.3$Group <- as.character(Subset.3$Group)
Subset$Interim <- as.character(Subset$Interim)
Subset$MatchInterimGp <- charmatch(Subset$Interim,Subset.3$Group)
Subset$MatchInterimGp[is.element(Subset$MatchInterimGp,c(NA))]<- 0
for (i in (1:length(Subset$Interim))){
Subset$Catch[i] <- ifelse(Subset$InterimGpRelevant[i]==TRUE,
                    Subset.3$Imputed.Catch[Subset$MatchInterimGp[i]],
                    Subset$Catch[i])}
return(Subset)
} 

Imputed <- matrix(,nrow=nrow(dataset),ncol=2)
Interim <- ifelse((sum(Imputing.Block(dataset,blockname)$Check.Interim)>0),1,0)
if(Interim==1){
  Imputed <- Impute.Catch(Impute.Diffnt.Interims(Imputing.Block(dataset,blockname)))[,1:2]
} else {
  Imputed <- Impute.Catch.No.Interim(Imputing.Block(dataset,blockname))[,1:2] 
}

return(Imputed)                
}

attributes(imputed)
#$dim
#[1] 16 24

#$dimnames
#$dimnames[[1]]
# [1] "1993" "1994" "1995" "1996" "1997" "1998" "1999" "2000" "2001" "2002" "2003"
#[12] "2004" "2005" "2006" "2007" "2008"

#$dimnames[[2]]
# [1] "11240" "11250" "11260" "12230" "12240" "12250" "12260" "13220" "13230" "13240"
#[11] "14210" "14220" "14230" "14240" "15200" "15210" "15220" "16200" "16210" "17200"
#[21] "17210" "18200" "18210" "19200"

copy.imputed <- imputed

for (i in 1:(ncol(imputed))){
 copy.imputed[,i] <- Impute.Block(imputed,colnames(imputed)[i])[,2]
 }
 
# Results.  
 
imputed # This is output from Genstat for model:
# pmul ~ boat+logdays+logother+c12+cs12+c6+cs6+calendaryr*block
# DISTRIBUTION=poisson; LINK=logarithm
# using predict [print=desc,pred; combinations=present;pred=a4_pred_cal;] calendaryr,block
 
#     11240 11250 11260 12230 12240 12250 12260 13220 13230 13240 14210 14220 14230
#1993     0     0     0  1200  1077    44     0     0  1739  1192   914  1271  3183
#1994     0     0     0  1042  1250     0     0     0  1924  1781  1162   861  1264
#1995     0     0     0   958   960     0     0     0   960   264   115  1100   896
#1996     0     0     0  1146  1474  1083     0     0  1146  1696     0   601   757
#1997     0     0     0   695  1167     0     0     0  1170   859     0  1036   835
#1998     0     0     0     0     0     0     0     0     0     0     0  2298  3272
#1999     0     0     0  3322  2806  1572  2810  1433  2295     0  2579  3043  4421
#2000     0  2975     0   742  1310  2295  2059     0  1506  2096     0  1441  1556
#2001     0     0     0  1199  1969  1488   908  1120     0  2181  2559  2820  1437
#2002  3682  4261     0  2409  2830  2373   374  4094     0  2381  1880  2627  1795
#2003     0  2604     0  1662  2350  2482     0     0  2941  4189  4353  3041  4521
#2004     0     0     0  5487  3462  1722     0     0  1787  3141   302  2276  3271
#2005  3760  3931  5031  5144  4196  3614     0     0  5494  4623     0  4491  4159
#2006  7427     0  5161  7570  6372  3788  5161     0  5389  3219   626  2715  3646
#2007 11552  4962  5139  3305  6388  4126     0     0  3723  5382     0  2836  4015
#2008     0  5599     0     0 14181     0     0     0  5254  3719  1972  3518  3234
#     14240 15200 15210 15220 16200 16210 17200 17210 18200 18210 19200
#1993     0     0  1713  1370  3355   827  1059   800     0     0     0
#1994     0     0   621   743  1047   502  1251   982    30     0     0
#1995  1236     0  1615  1387   717   492   960   544   233   196     0
#1996     0     0   738  1923   607   852  1417  1399   402     0     0
#1997     0     0   974   605  1441   684   693   909   607     0     0
#1998     0  2046   596     0  1940  1542  2238  1514   381     0     0
#1999     0     0     0     0  2294  3309  2298  2572  1043     0     0
#2000     0     0  1223     0   797  2456  1480  1186   608     0     0
#2001     0  2912  1712     0  1437  2793  1331  2657   309     0     0
#2002     0   984     0     0     0  1560  2015  1983     0     0     0
#2003     0     0  1744     0  2102  2000  4720  2389   299     0     0
#2004  2093   851  2241   672  2581  2698  4099  3111   946     0   207
#2005  5579     0  3085     0  8156  3919  5704  1374   117  1052    78
#2006  1752     0  4603   570     0  3752  4974  1910    35    29     0
#2007     0     0  2763  1580     0  1857  4968  2945   592  1791   713
#2008     0     0  1477  2099  2098  2167  6253  2518  1139  1583   114

copy.imputed # This is the imputed dataset using above R code:
#         11240  11250    11260  12230   12240  12250    12260    13220  13230  13240
#1993  4956.333 3280.0 5110.333 1200.0  1077.0   44.0 1925.667 2215.667 1739.0 1192.0
#1994  4956.333 3280.0 5110.333 1042.0  1250.0  563.5 1925.667 2215.667 1924.0 1781.0
#1995  4956.333 3280.0 5110.333  958.0   960.0  563.5 1925.667 2215.667  960.0  264.0
#1996  4956.333 3280.0 5110.333 1146.0  1474.0 1083.0 1925.667 2215.667 1146.0 1696.0
#1997  4956.333 3280.0 5110.333  695.0  1167.0 1327.5 1925.667 2215.667 1170.0  859.0
#1998  4956.333 3280.0 5110.333 2008.5  1986.5 1327.5 1925.667 2215.667 1732.5 1477.5
#1999  4956.333 3280.0 5110.333 3322.0  2806.0 1572.0 2810.000 1433.000 2295.0 1477.5
#2000  4956.333 2975.0 5110.333  742.0  1310.0 2295.0 2059.000 1276.500 1506.0 2096.0
#2001  4956.333 3618.0 5110.333 1199.0  1969.0 1488.0  908.000 1120.000 2223.5 2181.0
#2002  3682.000 4261.0 5110.333 2409.0  2830.0 2373.0  374.000 4094.000 2223.5 2381.0
#2003  3721.000 2604.0 5110.333 1662.0  2350.0 2482.0 2767.500 4094.000 2941.0 4189.0
#2004  3721.000 3267.5 5110.333 5487.0  3462.0 1722.0 2767.500 4094.000 1787.0 3141.0
#2005  3760.000 3931.0 5031.000 5144.0  4196.0 3614.0 2767.500 4094.000 5494.0 4623.0
#2006  7427.000 4446.5 5161.000 7570.0  6372.0 3788.0 5161.000 4094.000 5389.0 3219.0
#2007 11552.000 4962.0 5139.000 3305.0  6388.0 4126.0 5161.000 4094.000 3723.0 5382.0
#2008 11552.000 5599.0 5139.000 3305.0 14181.0 4126.0 5161.000 4094.000 5254.0 3719.0
#     14210 14220 14230    14240    15200  15210  15220  16200 16210 17200 17210
#1993   914  1271  3183 2969.333 1980.667 1713.0 1370.0 3355.0   827  1059   800
#1994  1162   861  1264 2969.333 1980.667  621.0  743.0 1047.0   502  1251   982
#1995   115  1100   896 1236.000 1980.667 1615.0 1387.0  717.0   492   960   544
#1996  1347   601   757 1664.500 1980.667  738.0 1923.0  607.0   852  1417  1399
#1997  1347  1036   835 1664.500 1980.667  974.0  605.0 1441.0   684   693   909
#1998  1347  2298  3272 1664.500 2046.000  596.0  638.5 1940.0  1542  2238  1514
#1999  2579  3043  4421 1664.500 2479.000  909.5  638.5 2294.0  3309  2298  2572
#2000  2569  1441  1556 1664.500 2479.000 1223.0  638.5  797.0  2456  1480  1186
#2001  2559  2820  1437 1664.500 2912.000 1712.0  638.5 1437.0  2793  1331  2657
#2002  1880  2627  1795 1664.500  984.000 1728.0  638.5 1769.5  1560  2015  1983
#2003  4353  3041  4521 1664.500  917.500 1744.0  638.5 2102.0  2000  4720  2389
#2004   302  2276  3271 2093.000  851.000 2241.0  672.0 2581.0  2698  4099  3111
#2005   464  4491  4159 5579.000  851.000 3085.0  621.0 8156.0  3919  5704  1374
#2006   626  2715  3646 1752.000  851.000 4603.0  570.0 5127.0  3752  4974  1910
#2007  1299  2836  4015 1752.000  851.000 2763.0 1580.0 5127.0  1857  4968  2945
#2008  1972  3518  3234 1752.000  851.000 1477.0 2099.0 2098.0  2167  6253  2518
#         18200     18210    19200
#1993  221.6667  425.6667 332.6667
#1994   30.0000  425.6667 332.6667
#1995  233.0000  196.0000 332.6667
#1996  402.0000  624.0000 332.6667
#1997  607.0000  624.0000 332.6667
#1998  381.0000  624.0000 332.6667
#1999 1043.0000  624.0000 332.6667
#2000  608.0000  624.0000 332.6667
#2001  309.0000  624.0000 332.6667
#2002  304.0000  624.0000 332.6667
#2003  299.0000  624.0000 332.6667
#2004  946.0000  624.0000 207.0000
#2005  117.0000 1052.0000  78.0000
#2006   35.0000   29.0000 395.5000
#2007  592.0000 1791.0000 713.0000
#2008 1139.0000 1583.0000 114.0000

# Export results to Excel Spreadsheet for Mike. 

write.table(copy.imputed, file = "GoldbandHistoricalCPUEImputedBlockxYr.dat")


 

    







                 
                        




                        













