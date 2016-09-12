
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

##Load Uganda panel datasets
Panel_Data_infra <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/uganda_infra.csv")
Panel_Data_soc <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/uganda_soc.csv")

## Subset social sector panel data and rename columns to prepare for merge
Panel_Data_soc=Panel_Data_soc[,c("ID","Year","MinYr","MinYr_additive","DecayYr","DecayYr_additive","DecayAddControl")]

names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr"] = "MinYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr_additive"] = "MinYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr"] = "DecayYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr_additive"] = "DecayYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayAddControl"] = "DecayAddControl_soc"

##Merge social and infrastructure datasets to create 1 panel dataset with cell distance to infra and social project locations
Panel_Data_Uganda<- merge(Panel_Data_infra,Panel_Data_soc, by=c("ID","Year"))

##Load Kenya panel datasets
Panel_Data_infra <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/kenya_infra.csv")
Panel_Data_soc <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/kenya_soc.csv")

## Subset social sector panel data and rename columns to prepare for merge
Panel_Data_soc=Panel_Data_soc[,c("ID","Year","MinYr","MinYr_additive","DecayYr","DecayYr_additive","DecayAddControl")]

names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr"] = "MinYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr_additive"] = "MinYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr"] = "DecayYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr_additive"] = "DecayYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayAddControl"] = "DecayAddControl_soc"

##Merge social and infrastructure datasets to create 1 panel dataset with cell distance to infra and social project locations
Panel_Data_Kenya<- merge(Panel_Data_infra,Panel_Data_soc, by=c("ID","Year"))

##Load Tanzania panel datasets
Panel_Data_infra <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra.csv")
Panel_Data_soc <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_social_22502cells.csv")

## Subset social sector panel data and rename columns to prepare for merge
Panel_Data_soc=Panel_Data_soc[,c("ID","Year","MinYr","MinYr_additive","DecayYr","DecayYr_additive","DecayAddControl")]

names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr"] = "MinYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr_additive"] = "MinYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr"] = "DecayYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr_additive"] = "DecayYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayAddControl"] = "DecayAddControl_soc"

##Merge social and infrastructure datasets to create 1 panel dataset with cell distance to infra and social project locations
Panel_Data_Tanz<- merge(Panel_Data_infra,Panel_Data_soc, by=c("ID","Year"))

#----------------------------------#
#Create Master Panel Dataset for Tanzania, Kenya, Uganda
#---------------------------------#
Panel_Data<-rbind(Panel_Data_Tanz,Panel_Data_Kenya)
Panel_Data1<-rbind(Panel_Data_Uganda,Panel_Data)
Panel_Data<-Panel_Data1


#--------------------------------------------------#
#Modeling
#--------------------------------------------------#

#infra + soc 

# initModel <- lm(Forest_Loss ~ DecayYr_additive + DecayAddControl + 
#                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
#                   Elevation + Slope + factor(District) + factor(Year), data=Panel_Data)
# cluster <- cluster.vcov(initModel, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
# CMREG <- coeftest(initModel, cluster)
# CMREG

initModel1 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayYr_additive_soc + DecayAddControl + DecayAddControl_soc + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + factor(District) + factor(Year), data=Panel_Data)
cluster1 <- cluster.vcov(initModel1, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1 <- coeftest(initModel1, cluster1)
CMREG1

initModel2 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayYr_additive_soc + DecayAddControl + DecayAddControl_soc + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + Year + factor(District), data=Panel_Data)
cluster2 <- cluster.vcov(initModel2, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2 <- coeftest(initModel2, cluster2)
CMREG2

#infra only

initModel3 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + factor(District) + factor(Year), data=Panel_Data)
cluster3 <- cluster.vcov(initModel3, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3 <- coeftest(initModel3, cluster3)
CMREG3

initModel4 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + Year + factor(District), data=Panel_Data)
cluster4 <- cluster.vcov(initModel4, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4 <- coeftest(initModel4, cluster4)
CMREG4

#soc only

initModel5 <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + factor(Year) + factor(District), data=Panel_Data)
cluster5 <- cluster.vcov(initModel5, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5 <- coeftest(initModel5, cluster5)
CMREG5

initModel6 <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + Year + factor(District), data=Panel_Data)
cluster6 <- cluster.vcov(initModel6, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG6 <- coeftest(initModel6, cluster6)
CMREG6




stargazer(CMREG3, CMREG4, CMREG5, CMREG6, CMREG1, CMREG2, type="html", keep=c("additive","Control","Min","Max","Elevation","Slope","Year"))


