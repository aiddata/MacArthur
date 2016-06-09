
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

##Panel Data, Thresh=10
Panel_Data <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia.csv")

#Add Post-2003 indicator
Panel_Data$post_2003<-0
Panel_Data$post_2003[Panel_Data$Year>2003]<-1


Panel_Data_test <- Panel_Data
Panel_Data_test=Panel_Data_test[Panel_Data_test$MinYr!=0,]
min <- aggregate(MinYr ~ ID, Panel_Data_test, function(x) min(x))
colnames(min)[2] <- "MinDist"
Panel_Data_min <- merge(Panel_Data, min, by.x="ID", by.y="ID")

mean<- aggregate(DecayYr_additive~Year, Panel_Data, function(x) mean(x))
mean

##Panel Data, Thresh=5

Panel_Data_5 <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_Thresh5.csv")

#Add Post-2003 indicator
Panel_Data_5$post_2003<-0
Panel_Data_5$post_2003[Panel_Data_5$Year>2003]<-1

##Panel Data, Thresh=15

Panel_Data_15<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_Thresh15.csv")

#Add Post-2003 indicator
Panel_Data_15$post_2003<-0
Panel_Data_15$post_2003[Panel_Data_15$Year>2003]<-1

#--------------------------------------------------#
#Modeling
#--------------------------------------------------#
# initModel <- lm(Forest_Loss ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
#                   MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + 
#                   Slope + factor(District) + factor(Year), data=Panel_Data)
# cluster <- cluster.vcov(initModel, cbind(Panel_Data$Year, Panel_Data$ID_2), force_posdef=TRUE)
# CMREG <- coeftest(initModel, cluster)

## Cumulative Forest Loss, Thresh=10
Model1<- lm(Forest_Loss_additive ~ DecayYr_additive, data=Panel_Data)
cluster1 <- cluster.vcov(Model1, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1 <- coeftest(Model1, cluster1)

Model1.1<- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl, data=Panel_Data)
cluster1.1 <- cluster.vcov(Model1.1, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1.1 <- coeftest(Model1.1, cluster1.1)

Model2 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                   MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                   UrbTravTime, data=Panel_Data)
cluster2 <- cluster.vcov(Model2, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2 <- coeftest(Model2, cluster2)

Model3 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(District), data=Panel_Data)
cluster3 <- cluster.vcov(Model3, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3 <- coeftest(Model3, cluster3)

Model4 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster4 <- cluster.vcov(Model4, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4 <- coeftest(Model4, cluster4)

Model5 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
                   PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                    MaxPrecip + MeanPrecip + MinPrecip + 
                   Elevation + Slope + UrbTravTime + Year + Year*post_2003 + factor(District), data=Panel_Data)
cluster5 <- cluster.vcov(Model5, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5 <- coeftest(Model5, cluster5)

##Cumulative Forest Loss, treatment is distance decay within 100km

Model400d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400d <- cluster.vcov(Model400d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400d <- coeftest(Model400d, cluster400d)


##Cumulative Forest Loss, treatment is project count within 100km,Thresh=10

Model400 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400 <- cluster.vcov(Model400, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400 <- coeftest(Model400, cluster400)

Model500 <- lm(Forest_Loss_additive ~ ProjCnt100 + 
               PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
               MaxPrecip + MeanPrecip + MinPrecip + 
               Elevation + Slope + UrbTravTime + Year + Year*post_2003 + factor(District), data=Panel_Data)
cluster500 <- cluster.vcov(Model500, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG500 <- coeftest(Model500, cluster500)

##Cumulative Forest Loss, Thresh=5

Model51 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(Year) + factor(District), data=Panel_Data_5)
cluster51 <- cluster.vcov(Model51, cbind(Panel_Data_5$Year, Panel_Data_5$District), force_posdef=TRUE)
CMREG51 <- coeftest(Model51, cluster51)

Model52 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
               PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
               MaxPrecip + MeanPrecip + MinPrecip + 
               Elevation + Slope + UrbTravTime + Year + Year*post_2003 + factor(District), data=Panel_Data_5)
cluster52 <- cluster.vcov(Model52, cbind(Panel_Data_5$Year, Panel_Data_5$District), force_posdef=TRUE)
CMREG52 <- coeftest(Model52, cluster52)

##Cumulative Forest Loss, Thresh=15

Model151 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + factor(Year) + factor(District), data=Panel_Data_15)
cluster151 <- cluster.vcov(Model151, cbind(Panel_Data_15$Year, Panel_Data_15$District), force_posdef=TRUE)
CMREG151 <- coeftest(Model151, cluster151)

Model152 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
                PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                MaxPrecip + MeanPrecip + MinPrecip + 
                Elevation + Slope + UrbTravTime + Year + Year*post_2003 + factor(District), data=Panel_Data_15)
cluster152 <- cluster.vcov(Model152, cbind(Panel_Data_15$Year, Panel_Data_15$District), force_posdef=TRUE)
CMREG152 <- coeftest(Model152, cluster152)

##Cumulative Forest Loss, treatment is project count within 100km,Thresh=10



#-------------------------#
#Stargazer#
#-------------------------#

stargazer(CMREG1, CMREG4, type="html", keep=c("Forest_Loss","additive","Control","Min","Max","Elevation","Slope","Year"))

stargazer(CMREG1, CMREG1.1, CMREG2, CMREG3, CMREG4, CMREG5,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("District Fixed Effects?","No","No","No","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","No","Yes","No")),
          title="Cambodia Infra Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG400, CMREG500, CMREG51, CMREG52, CMREG151, CMREG152,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post","Proj"),
          omit.stat=c("f","ser"),
          add.lines=list(c("District Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","Yes","No","Yes","No","Yes","No"),
                         c("Threshold=5?","No","No","Yes","Yes","No","No"),
                         c("Threshold=10?","Yes","Yes","No","No","No","No"),
                         c("Threshold=15","No","No","No","No","Yes","Yes")),
          title="Cambodia Infra Regression Results: Alternate Treatment",
          dep.var.labels=c("Cumulative Forest Loss"))

