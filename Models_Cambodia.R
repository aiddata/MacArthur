
library(sp)
library(jsonlite)
library(xtable)
library(stargazer)
library(RCurl)
library(maptools)
library(geosphere)
library(ncf)
library(gridExtra)
require(grid)
library(multiwayvcov)
library(lmtest)


Panel_Data <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia.csv")

Panel_Data_test <- Panel_Data
Panel_Data_test=Panel_Data_test[Panel_Data_test$MinYr!=0,]
min <- aggregate(MinYr ~ ID, Panel_Data_test, function(x) min(x))
colnames(min)[2] <- "MinDist"
Panel_Data_min <- merge(Panel_Data, min, by.x="ID", by.y="ID")



#--------------------------------------------------#
#Modeling
#--------------------------------------------------#
initModel <- lm(Forest_Loss ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + factor(District) + factor(Year), data=Panel_Data)
cluster <- cluster.vcov(initModel, cbind(Panel_Data$Year, Panel_Data$ID_2), force_posdef=TRUE)
CMREG <- coeftest(initModel, cluster)

initModel1 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + Year + factor(District), data=Panel_Data)
cluster1 <- cluster.vcov(initModel1, cbind(Panel_Data$Year, Panel_Data$ID_2), force_posdef=TRUE)
CMREG1 <- coeftest(initModel1, cluster1)
