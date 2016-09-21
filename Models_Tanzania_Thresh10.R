
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

##Load panel datasets
Panel_Data_infra <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra.csv")
Panel_Data_infra_add <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add.csv")
Panel_Data_infra_add_aug<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add_AUG.csv") 
Panel_Data_infra <- Panel_Data_infra_add_aug
Panel_Data_soc <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_soc_AUG.csv")

## Subset social sector panel data and rename columns to prepare for merge
Panel_Data_soc=Panel_Data_soc[,c("ID","Year","MinYr","DecayYr","DecayYr_additive",
                                 "DecayAddControl","ProjCnt100","ProjCnt100_additive","DecayYr100",
                                 "DecayYr100_additive","DecayAddControl100")]

names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr"] = "MinYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "MinYr_additive"] = "MinYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr"] = "DecayYr_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr_additive"] = "DecayYr_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayAddControl"] = "DecayAddControl_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "ProjCnt100"] = "ProjCnt100_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "ProjCnt100_additive"] = "ProjCnt100_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr100"] = "DecayYr100_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayYr100_additive"] = "DecayYr100_additive_soc"
names(Panel_Data_soc)[names(Panel_Data_soc) == "DecayAddControl100"] = "DecayAddControl100_soc"


##Merge social and infrastructure datasets to create 1 panel dataset with cell distance to infra and social project locations
Panel_Data<- merge(Panel_Data_infra,Panel_Data_soc, by=c("ID","Year"))

#Add Post-2007 indicator (because no projects until 2008)
Panel_Data$post_2007<-0
Panel_Data$post_2007[Panel_Data$Year>2007]<-1

#Create lagged treatment variables for both infra and social sectors

#create lagged treatment variables
lagpad <- function(x, k=1) {
  i<-is.vector(x)
  if(is.vector(x)) x<-matrix(x) else x<-matrix(x,nrow(x))
  if(k>0) {
    x <- rbind(matrix(rep(NA, k*ncol(x)),ncol=ncol(x)), matrix(x[1:(nrow(x)-k),], ncol=ncol(x)))
  }
  else {
    x <- rbind(matrix(x[(-k+1):(nrow(x)),], ncol=ncol(x)),matrix(rep(NA, -k*ncol(x)),ncol=ncol(x)))
  }
  if(i) x[1:length(x)] else x
}
#infra sector treatment
Panel_Data$treat_minus1<-NA
Panel_Data$treat_minus1<-lagpad(Panel_Data$DecayYr_additive,-1)
Panel_Data$treat_minus1[Panel_Data$Year==2014]<-NA

Panel_Data$treat_minus2<-NA
Panel_Data$treat_minus2<-lagpad(Panel_Data$DecayYr_additive,-2)
Panel_Data$treat_minus2[Panel_Data$Year>=2013]<-NA

Panel_Data$treat_minus3<-NA
Panel_Data$treat_minus3<-lagpad(Panel_Data$DecayYr_additive,-3)
Panel_Data$treat_minus3[Panel_Data$Year>=2012]<-NA

Panel_Data$treat_minus4<-NA
Panel_Data$treat_minus4<-lagpad(Panel_Data$DecayYr_additive,-4)
Panel_Data$treat_minus4[Panel_Data$Year>=2011]<-NA

Panel_Data$treat_minus5<-NA
Panel_Data$treat_minus5<-lagpad(Panel_Data$DecayYr_additive,-5)
Panel_Data$treat_minus5[Panel_Data$Year>=2010]<-NA

#social sector treatment
Panel_Data$treat_minus1_soc<-NA
Panel_Data$treat_minus1_soc<-lagpad(Panel_Data$DecayYr_additive_soc,-1)
Panel_Data$treat_minus1_soc[Panel_Data$Year==2014]<-NA

Panel_Data$treat_minus2_soc<-NA
Panel_Data$treat_minus2_soc<-lagpad(Panel_Data$DecayYr_additive_soc,-2)
Panel_Data$treat_minus2_soc[Panel_Data$Year>=2013]<-NA

Panel_Data$treat_minus3_soc<-NA
Panel_Data$treat_minus3_soc<-lagpad(Panel_Data$DecayYr_additive_soc,-3)
Panel_Data$treat_minus3_soc[Panel_Data$Year>=2012]<-NA

Panel_Data$treat_minus4_soc<-NA
Panel_Data$treat_minus4_soc<-lagpad(Panel_Data$DecayYr_additive_soc,-4)
Panel_Data$treat_minus4_soc[Panel_Data$Year>=2011]<-NA

Panel_Data$treat_minus5_soc<-NA
Panel_Data$treat_minus5_soc<-lagpad(Panel_Data$DecayYr_additive_soc,-5)
Panel_Data$treat_minus5_soc[Panel_Data$Year>=2010]<-NA

Panel_Data_trtlag <- Panel_Data[Panel_Data$Year<=2009,]

#-------------------------
#Identifying Cells for Visualizations
#-------------------------

Panel_Data2014<-Panel_Data[Panel_Data$Year==2014,]
Panel_Data2014$Forest_Loss_end<-0
Panel_Data2014$Forest_Loss_end=Panel_Data2014$Forest_Loss_additive
Panel_Data2014<-Panel_Data2014[,c("ID","Forest_Loss_end")]

Panel_Data07<-Panel_Data[Panel_Data$Year==2007,]
Panel_Data07<-merge(Panel_Data07, Panel_Data2014, by.x="ID", by.y="ID")
Panel_Data07<-Panel_Data07[Panel_Data07$Forest_Loss_additive==0,]
Panel_Data07<-Panel_Data07[Panel_Data07$Forest_Loss_additive==Panel_Data07$Forest_Loss_end,]



# ##Correlogram
# #Load and identify x intercept
# load("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_correl.RData")
# correlogram_data$x.intercept
# plot(correlogram_data)
# 
# #Trim data by x intercept threshold
# Panel_Data_test <- Panel_Data
# Panel_Data_test=Panel_Data_test[Panel_Data_test$MinYr!=0,]
# min <- aggregate(MinYr ~ ID, Panel_Data_test, function(x) min(x))
# #the max value when you run the below should be less than or equal to correlogram_data$x.intercept
# summary(min$MinYr)
# 
# colnames(min)[2] <- "MinDist"
# Panel_Data_min <- merge(Panel_Data, min, by.x="ID", by.y="ID")

#---------------
#Reading in panel data set at different standing forest thresholds
#---------------

##Panel Data, Thresh=5

#Panel_Data_5 <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_Thresh5.csv")



##Panel Data, Thresh=15

#Panel_Data_15<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_Thresh15.csv")



#--------------------------------------------------#
#Modeling
#--------------------------------------------------#
# initModel <- lm(Forest_Loss ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
#                   MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + 
#                   Slope + factor(District) + factor(Year), data=Panel_Data)
# cluster <- cluster.vcov(initModel, cbind(Panel_Data$Year, Panel_Data$ID_2), force_posdef=TRUE)
# CMREG <- coeftest(initModel, cluster)

## Cumulative Forest Loss, Thresh=10, Infra
Model1<- lm(Forest_Loss_additive ~ DecayYr_additive, data=Panel_Data)
cluster1 <- cluster.vcov(Model1, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1 <- coeftest(Model1, cluster1)

Model1.1<- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl, data=Panel_Data)
cluster1.1 <- cluster.vcov(Model1.1, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1.1 <- coeftest(Model1.1, cluster1.1)

Model2 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl +
             MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope +
             UrbTravTime + 
             ntl_pretrend + 
             Pop +
             factor(District),   
             data=Panel_Data)
cluster2 <- cluster.vcov(Model2, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2 <- coeftest(Model2, cluster2)

Model3 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl +
              MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope +
              UrbTravTime + 
              ntl_pretrend + 
              Pop +
              factor(District) + factor(Year),   
              data=Panel_Data)
cluster3 <- cluster.vcov(Model3, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3 <- coeftest(Model3, cluster3)

Model4 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
               Pop + DecayYr_additive*Pop+
               factor(Year) + factor(District), data=Panel_Data)
cluster4 <- cluster.vcov(Model4, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4 <- coeftest(Model4, cluster4)

Model5 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + 
               MaxPrecip + MeanPrecip + MinPrecip + 
               Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
               Pop+
               wdpapct_2007 + DecayYr_additive*wdpapct_2007+
               factor(Year) + factor(District), data=Panel_Data)
cluster5 <- cluster.vcov(Model5, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5 <- coeftest(Model5, cluster5)

Model6 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend +
               Pop + DecayYr_additive*Pop +
               wdpapct_2007 + DecayYr_additive*wdpapct_2007+ 
               Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster6<-cluster.vcov(Model6,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG6 <- coeftest(Model6,cluster6)


Model9 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend  +
               Pop + DecayYr_additive*Pop +
               wdpapct_2007 + DecayYr_additive*wdpapct_2007+
               factor(Year) + factor(District), data=Panel_Data)
cluster9<-cluster.vcov(Model9,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG9 <- coeftest(Model9,cluster9)

Model10 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend + 
               Pop + DecayYr_additive*Pop +
               wdpapct_2007 + DecayYr_additive*wdpapct_2007+
               NTL+
               factor(Year) + factor(District), data=Panel_Data)
cluster10<-cluster.vcov(Model10,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG10 <- coeftest(Model10,cluster10)


Model11 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime  + 
                ntl_pretrend  +
                Pop + DecayYr_additive*Pop +
                wdpapct_2007 + DecayYr_additive*wdpapct_2007+
                treat_minus1 + treat_minus2+treat_minus3+treat_minus4+treat_minus5+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster11<-cluster.vcov(Model11,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG11 <- coeftest(Model11,cluster11)

Model12 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + DecayYr_additive*Pop +
                wdpapct_2007 + DecayYr_additive*wdpapct_2007+
                treat_minus1 + treat_minus2+treat_minus3+treat_minus4+treat_minus5+
                NTL+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster12<-cluster.vcov(Model12,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG12 <- coeftest(Model12,cluster12)


##Cumulative Forest Loss, thresh=10, social

Model1s<- lm(Forest_Loss_additive ~ DecayYr_additive_soc, data=Panel_Data)
cluster1s <- cluster.vcov(Model1s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1s <- coeftest(Model1s, cluster1s)

Model1.1s<- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc, data=Panel_Data)
cluster1.1s <- cluster.vcov(Model1.1s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1.1s <- coeftest(Model1.1s, cluster1.1s)

Model2s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime +
                ntl_pretrend+
                Pop+
                factor(District),
              data=Panel_Data)
cluster2s <- cluster.vcov(Model2s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2s <- coeftest(Model2s, cluster2s)


Model3s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime +
                ntl_pretrend+
                Pop+
                factor(Year) + factor(District), data=Panel_Data)
cluster3s <- cluster.vcov(Model3s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3s <- coeftest(Model3s, cluster3s)

Model4s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  +
               ntl_pretrend +
               Pop + DecayYr_additive_soc*Pop +
               factor(Year) + factor(District), data=Panel_Data)
cluster4s <- cluster.vcov(Model4s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4s <- coeftest(Model4s, cluster4s)

Model5s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + 
                ntl_pretrend+
                Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data)
cluster5s <- cluster.vcov(Model5s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5s <- coeftest(Model5s, cluster5s)

Model6s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 +
               ntl_pretrend + NTL_2007+
               wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
               factor(Year) + factor(District), data=Panel_Data)
cluster6s<-cluster.vcov(Model6s,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG6s <- coeftest(Model6s,cluster6s)

Model9s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 +
               ntl_pretrend + NTL_2007+
               Pop + DecayYr_additive_soc*Pop+
               wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
               factor(Year) + factor(District), data=Panel_Data)
cluster9s<-cluster.vcov(Model9s,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG9s <- coeftest(Model9s,cluster9s)

Model10s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + Pop_2000 +
                ntl_pretrend + NTL_2007+NTL+
                Pop + DecayYr_additive_soc*Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                factor(Year) + factor(District), data=Panel_Data)
cluster10s<-cluster.vcov(Model10s,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG10s <- coeftest(Model10s,cluster10s)

Model11s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + Pop_2000 +
                 ntl_pretrend + NTL_2007+
                 Pop + DecayYr_additive_soc*Pop+
                 wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                 treat_minus1_soc + treat_minus2_soc + treat_minus3_soc+treat_minus4_soc+treat_minus5_soc+
                 factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster11s<-cluster.vcov(Model11s,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG11s <- coeftest(Model11s,cluster11s)

Model12s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + Pop_2000 +
                ntl_pretrend + NTL_2007+NTL+
                Pop + DecayYr_additive_soc*Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                treat_minus1_soc + treat_minus2_soc + treat_minus3_soc+treat_minus4_soc+treat_minus5_soc+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster12s<-cluster.vcov(Model12s,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG12s <- coeftest(Model12s,cluster12s)



##Cumulative Forest Loss, Infra + Soc, Thresh=10

Model1si<- lm(Forest_Loss_additive ~ DecayYr_additive + DecayYr_additive_soc, data=Panel_Data)
cluster1si <- cluster.vcov(Model1si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1si <- coeftest(Model1si, cluster1si)

Model1.1si<- lm(Forest_Loss_additive ~ DecayYr_additive + DecayYr_additive_soc + DecayAddControl +DecayAddControl_soc, 
               data=Panel_Data)
cluster1.1si <- cluster.vcov(Model1.1si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG1.1si <- coeftest(Model1.1si, cluster1.1si)

Model2si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime, data=Panel_Data)
cluster2si <- cluster.vcov(Model2si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2si <- coeftest(Model2si, cluster2si)

Model3si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + factor(District), data=Panel_Data)
cluster3si <- cluster.vcov(Model3si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3si <- coeftest(Model3si, cluster3si)

Model4si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + UrbTravTime*DecayYr_additive+UrbTravTime*DecayYr_additive_soc+
                 Pop_2000 + Pop + Pop*DecayYr_additive+Pop*DecayYr_additive_soc+
                 factor(Year) + factor(District), data=Panel_Data)
cluster4si <- cluster.vcov(Model4si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4si <- coeftest(Model4si, cluster4si)

Model5si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + 
                 PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                 MaxPrecip + MeanPrecip + MinPrecip + 
                 Elevation + Slope + 
                 UrbTravTime +UrbTravTime*DecayYr_additive+UrbTravTime*DecayYr_additive_soc+
                 Pop_2000 + Pop + Pop*DecayYr_additive+Pop*DecayYr_additive_soc+
                 Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster5si <- cluster.vcov(Model5si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5si <- coeftest(Model5si, cluster5si)

Model6si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + UrbTravTime*DecayYr_additive+UrbTravTime*DecayYr_additive_soc+
                 Pop_2000 + Pop + Pop*DecayYr_additive+Pop*DecayYr_additive_soc+
                 wdpapct_2007 + wdpapct_2007*DecayYr_additive + wdpapct_2007*DecayYr_additive_soc+
                 NTL + ntl_pretrend + NTL_2007+
                 treat_minus1 + treat_minus2 + treat_minus3 + treat_minus4 + treat_minus5+
                 treat_minus1_soc + treat_minus2_soc + treat_minus3_soc + treat_minus4_soc + treat_minus5_soc+
                 factor(Year) + factor(District), data=Panel_Data)
cluster6si <- cluster.vcov(Model6si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG6si <- coeftest(Model6si, cluster6si)

Mode75si <- lm(Forest_Loss_additive ~ DecayYr_additive+DecayYr_additive_soc + DecayAddControl+DecayAddControl_soc + 
                 PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                 MaxPrecip + MeanPrecip + MinPrecip + 
                 Elevation + Slope + 
                 UrbTravTime +UrbTravTime*DecayYr_additive+UrbTravTime*DecayYr_additive_soc+
                 Pop_2000 + Pop + Pop*DecayYr_additive+Pop*DecayYr_additive_soc+
                 wdpapct_2007 + wdpapct_2007*DecayYr_additive + wdpapct_2007*DecayYr_additive_soc+
                 NTL + ntl_pretrend + NTL_2007+
                 treat_minus1 + treat_minus2 + treat_minus3 + treat_minus4 + treat_minus5+
                 treat_minus1_soc + treat_minus2_soc + treat_minus3_soc + treat_minus4_soc + treat_minus5_soc+
                 Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster7si <- cluster.vcov(Model7si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG7si <- coeftest(Model7si, cluster7si)


##Cumulative Forest Loss, Infra + Soc, Thresh=10, Distance Decay applied only for projects within 100km

Model104d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster104d <- cluster.vcov(Model104d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG104d <- coeftest(Model104d, cluster4)

Model105d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + 
               PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
               MaxPrecip + MeanPrecip + MinPrecip + 
               Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster105d <- cluster.vcov(Model105d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG105d <- coeftest(Model105d, cluster105d)

Model104ds <- lm(Forest_Loss_additive ~ DecayYr100_additive_soc + DecayAddControl100_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster104ds <- cluster.vcov(Model104ds, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG104ds <- coeftest(Model104ds, cluster104ds)

Model105ds <- lm(Forest_Loss_additive ~ DecayYr100_additive_soc + DecayAddControl100_soc + 
                PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                MaxPrecip + MeanPrecip + MinPrecip + 
                Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster105ds <- cluster.vcov(Model105ds, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG105ds <- coeftest(Model105ds, cluster105ds)

Model104dsi <- lm(Forest_Loss_additive ~ DecayYr100_additive+DecayYr100_additive_soc + DecayAddControl100+
                    DecayAddControl100_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster104dsi <- cluster.vcov(Model104dsi, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG104dsi <- coeftest(Model104dsi, cluster104dsi)

Model105dsi <- lm(Forest_Loss_additive ~ DecayYr100_additive+DecayYr100_additive_soc + DecayAddControl100+
                    DecayAddControl100_soc + 
                 PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                 MaxPrecip + MeanPrecip + MinPrecip + 
                 Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster105dsi <- cluster.vcov(Model105dsi, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG105dsi <- coeftest(Model105dsi, cluster105dsi)

##Cumulative Forest Loss, treatment is project count within 100km,Thresh=10

Model400 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400 <- cluster.vcov(Model400, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400 <- coeftest(Model400, cluster400)

Model500 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + 
                 PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                 MaxPrecip + MeanPrecip + MinPrecip + 
                 Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster500 <- cluster.vcov(Model500, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG500 <- coeftest(Model500, cluster500)

Model400s <- lm(Forest_Loss_additive ~ ProjCnt100_additive_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400s <- cluster.vcov(Model400s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400s <- coeftest(Model400s, cluster400s)

Model500s <- lm(Forest_Loss_additive ~ ProjCnt100_additive_soc + 
                 PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                 MaxPrecip + MeanPrecip + MinPrecip + 
                 Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster500s <- cluster.vcov(Model500s, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG500s <- coeftest(Model500s, cluster500s)

Model400si <- lm(Forest_Loss_additive ~ ProjCnt100_additive+ ProjCnt100_additive_soc + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400si <- cluster.vcov(Model400si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400si <- coeftest(Model400si, cluster400si)

Model500si <- lm(Forest_Loss_additive ~ ProjCnt100_additive + ProjCnt100_additive_soc + 
                  PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                  MaxPrecip + MeanPrecip + MinPrecip + 
                  Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data)
cluster500si <- cluster.vcov(Model500si, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG500si <- coeftest(Model500si, cluster500si)

##Cumulative Forest Loss, Thresh=5

Model51 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + factor(Year) + factor(District), data=Panel_Data_5)
cluster51 <- cluster.vcov(Model51, cbind(Panel_Data_5$Year, Panel_Data_5$District), force_posdef=TRUE)
CMREG51 <- coeftest(Model51, cluster51)

Model52 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
                PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
                MaxPrecip + MeanPrecip + MinPrecip + 
                Elevation + Slope + UrbTravTime + Year + Year*post_2007 + factor(District), data=Panel_Data_5)
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

#WORKING PAPER SUMMARY STATS#
#rename vars directly in html file
#treatment
stargazer(Panel_Data, type = "html", nobs = FALSE, mean.sd = TRUE, median = TRUE,
          iqr = TRUE,
          keep=c("Forest_Loss_additive","DecayYr_additive","DecayYr100_additive","ProjCnt100_additive",
                 "DecayAddControl","DecayAddControl100"))
#covars
stargazer(Panel_Data,type="html", 
          keep=c("MinTemp","MinPrecip","Max","Mean","Elevation","Slope","UrbTravTime","ntl_pretrend","NTL",
                 "Pop","Pct","PreLevelControl","PreTrendControl"))


stargazer(CMREG1, CMREG4, type="html", 
          keep=c("Forest_Loss","additive","Control","Min","Max","Elevation","Slope","Year"))

stargazer(CMREG1, CMREG1.1, CMREG5,CMREG4,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes")),
          title="Tanzania Infra Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

#WORKING PAPER RESULTS TABLE#
#have to take out automatically inserted "factor" lines directly in html
stargazer(CMREG1, CMREG1.1, CMREG2, CMREG3,CMREG9,CMREG10,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),omit.labels=c("factor","Temp","Precip"),
          omit.stat=c("f","ser"),
          covariate.labels=c("Treatment (Proximity)","Proximity Control","Baseline NDVI","NDVI Pre-Trend",
                             "Elevation","Slope",
                             "Urban Travel Time","Nighttime Lights Pre-Trend","Population",
                             "Baseline Protected Areas",
                             "Nighttime Lights",
                             "Population*Treatment","Protected Area*Treatment"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028","315,028","315,028","315,028"),
                         c("Climate Controls?","No","No","Yes","Yes","Yes","Yes"),
                         c("District Fixed Effects?","No","No","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes","Yes")),
          title="Tanzania Infrastructure Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG4s,CMREG6s,CMREG9s,CMREG10s,CMREG11s,CMREG12s,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),omit.labels=c("factor","Temp","Precip"),
          omit.stat=c("f","ser"),
          # add.lines=list(c("Observations","58,996","58,996","58,996","58,996"),
          #                c("District Fixed Effects?","No","No","Yes","Yes"),
          #                c("Year Fixed Effects?","No","No","No","Yes")),
          title="Tanzania Soc Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG1s, CMREG1.1s, CMREG5s, CMREG4s,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes")),
          title="Tanzania Soc Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG1si, CMREG1.1si, CMREG5si, CMREG4si,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes")),
          title="Tanzania Infra+Soc Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG500, CMREG400, CMREG500s, CMREG400s, CMREG400si,CMREG500si,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post","Proj"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","Yes","No","Yes","No","Yes")),
          title="Tanzania Infra+Soc Regression Results: 100km Project Count",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG105d, CMREG104d, CMREG105ds, CMREG104ds, CMREG105dsi, CMREG104dsi,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post","Proj"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","Yes","No","Yes","No","Yes")),
          title="Tanzania Infra+Soc Regression Results: 100km Distance Decay",
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
          title="Tanzania Regression Results: Alternate Treatment",
          dep.var.labels=c("Cumulative Forest Loss"))


