
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
#Panel_Data_infra <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra.csv")
#Panel_Data_infra_add <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add.csv")
#Panel_Data_infra_add_aug<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add_AUG.csv") 
Panel_Data<-read.csv("/Users/rbtrichler/Box Sync/MacArthur/modelData/tanzania_infra_add_oct2017.csv")

Panel_Data <- Panel_Data_infra_add_aug
#Panel_Data_soc <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_soc_ONLY_AUG.csv")
Panel_Data_soc<-read.csv("/Users/rbtrichler/Box Sync/MacArthur/tanzania_soc_oct2017.csv")


## Subset social sector panel data and rename columns to prepare for merge with infrastructure panel dataset 
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

##Merge in outcome, treatment, covariate data from Panel_Data_infra
#Panel_Data_soc will now have both infra and social treatment values for overlapping cells, which should number 300,314 total in dataset
Panel_Data_soc<-merge(Panel_Data_soc,Panel_Data, by=c("ID","Year"))

#Add Post-2007 indicator (because no projects until 2008) for both infra and social sector datasets
Panel_Data$post_2007<-0
Panel_Data$post_2007[Panel_Data$Year>2007]<-1

Panel_Data_soc$post_2007<-0
Panel_Data_soc$post_2007[Panel_Data_soc$Year>2007]<-1

#Create lead treatment variables for both infra and social sectors

#function to create lead treatment variables
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
Panel_Data_soc$treat_minus1_soc<-NA
Panel_Data_soc$treat_minus1_soc<-lagpad(Panel_Data_soc$DecayYr_additive_soc,-1)
Panel_Data_soc$treat_minus1_soc[Panel_Data_soc$Year==2014]<-NA

Panel_Data_soc$treat_minus2_soc<-NA
Panel_Data_soc$treat_minus2_soc<-lagpad(Panel_Data_soc$DecayYr_additive_soc,-2)
Panel_Data_soc$treat_minus2_soc[Panel_Data_soc$Year>=2013]<-NA

Panel_Data_soc$treat_minus3_soc<-NA
Panel_Data_soc$treat_minus3_soc<-lagpad(Panel_Data_soc$DecayYr_additive_soc,-3)
Panel_Data_soc$treat_minus3_soc[Panel_Data_soc$Year>=2012]<-NA

Panel_Data_soc$treat_minus4_soc<-NA
Panel_Data_soc$treat_minus4_soc<-lagpad(Panel_Data_soc$DecayYr_additive_soc,-4)
Panel_Data_soc$treat_minus4_soc[Panel_Data_soc$Year>=2011]<-NA

Panel_Data_soc$treat_minus5_soc<-NA
Panel_Data_soc$treat_minus5_soc<-lagpad(Panel_Data_soc$DecayYr_additive_soc,-5)
Panel_Data_soc$treat_minus5_soc[Panel_Data_soc$Year>=2010]<-NA

#Create data subset for 2001 to 2009 so no years with NA values for treatment leads
Panel_Data_trtlag <- Panel_Data[Panel_Data$Year<=2009,]
Panel_Data_trtlag_soc <- Panel_Data_soc[Panel_Data_soc$Year<=2009,]

#Create data subset for each year with treatment lead (i.e. take out one year of dataset at a time) and rename variables to make stargazer easier
Panel_Data_trtlag13 <- Panel_Data[Panel_Data$Year<=2013,]
names(Panel_Data_trtlag13)[names(Panel_Data_trtlag13) == "treat_minus1"] = "trtlead"

Panel_Data_trtlag12 <- Panel_Data[Panel_Data$Year<=2012,]
names(Panel_Data_trtlag12)[names(Panel_Data_trtlag12) == "treat_minus2"] = "trtlead"

Panel_Data_trtlag11 <- Panel_Data[Panel_Data$Year<=2011,]
names(Panel_Data_trtlag11)[names(Panel_Data_trtlag11) == "treat_minus3"] = "trtlead"

Panel_Data_trtlag10 <- Panel_Data[Panel_Data$Year<=2010,]
names(Panel_Data_trtlag10)[names(Panel_Data_trtlag10) == "treat_minus4"] = "trtlead"

names(Panel_Data_trtlag)[names(Panel_Data_trtlag) == "treat_minus5"] = "trtlead"


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

##Panel Data, Thresh=15, Infrastructure Treatment only

Panel_Data_15<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_Thresh15.csv")

#merge in additional covariate data from thresh=10 panel dataset
ntl<-Panel_Data[,c("ID","Year","wdpapct_2000","wdpapct_2007",
                    "Pop_2000","gpw_v4_density.2005.mean","gpw_v4_density.2010.mean","gpw_v4_density.2015.mean",
                    "Pop","ntltrend_0913","ntl_pretrend","NTL","NTL_2007")]
Panel_Data_15_ntl<-merge(Panel_Data_15,ntl,by=c("ID","Year"))
Panel_Data_15<-Panel_Data_15_ntl

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

Model13 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_trtlag13)
cluster13<-cluster.vcov(Model13,cbind(Panel_Data_trtlag13$Year,Panel_Data_trtlag13$District), force_posdef=TRUE)
CMREG13 <- coeftest(Model13,cluster13)

Model14 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl+PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_trtlag12)
cluster14<-cluster.vcov(Model14,cbind(Panel_Data_trtlag12$Year,Panel_Data_trtlag12$District), force_posdef=TRUE)
CMREG14 <- coeftest(Model14,cluster14)

Model15 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl+PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_trtlag11)
cluster15<-cluster.vcov(Model15,cbind(Panel_Data_trtlag11$Year,Panel_Data_trtlag11$District), force_posdef=TRUE)
CMREG15 <- coeftest(Model15,cluster15)

Model16 <- lm(Forest_Loss_additive ~ trtlead+DecayAddControl+ PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_trtlag10)
cluster16<-cluster.vcov(Model16,cbind(Panel_Data_trtlag10$Year,Panel_Data_trtlag10$District), force_posdef=TRUE)
CMREG16 <- coeftest(Model16,cluster16)

Model17 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl +PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster17<-cluster.vcov(Model17,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG17 <- coeftest(Model17,cluster17)




##Cumulative Forest Loss, thresh=10, social

Model1s<- lm(Forest_Loss_additive ~ DecayYr_additive_soc, data=Panel_Data_soc)
cluster1s <- cluster.vcov(Model1s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG1s <- coeftest(Model1s, cluster1s)

Model1.1s<- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc, data=Panel_Data_soc)
cluster1.1s <- cluster.vcov(Model1.1s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG1.1s <- coeftest(Model1.1s, cluster1.1s)

Model2s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime +
                ntl_pretrend+
                Pop+
                factor(District),
              data=Panel_Data_soc)
cluster2s <- cluster.vcov(Model2s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG2s <- coeftest(Model2s, cluster2s)


Model3s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime +
                ntl_pretrend+
                Pop+
                factor(Year) + factor(District), data=Panel_Data_soc)
cluster3s <- cluster.vcov(Model3s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG3s <- coeftest(Model3s, cluster3s)

Model4s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  +
               ntl_pretrend +
               Pop + DecayYr_additive_soc*Pop +
               factor(Year) + factor(District), data=Panel_Data_soc)
cluster4s <- cluster.vcov(Model4s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG4s <- coeftest(Model4s, cluster4s)

Model5s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + 
                ntl_pretrend+
                Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007+
                factor(Year) + factor(District), data=Panel_Data_soc)
cluster5s <- cluster.vcov(Model5s, cbind(Panel_Data_soc$Year, Panel_Data_soc$District), force_posdef=TRUE)
CMREG5s <- coeftest(Model5s, cluster5s)

Model6s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  +
               ntl_pretrend +
                Pop + DecayYr_additive_soc*Pop+
               wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
               Year + Year*post_2007, data=Panel_Data_soc)
cluster6s<-cluster.vcov(Model6s,cbind(Panel_Data_soc$Year,Panel_Data_soc$District), force_posdef=TRUE)
CMREG6s <- coeftest(Model6s,cluster6s)

Model9s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  +
               ntl_pretrend +
               Pop + DecayYr_additive_soc*Pop+
               wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
               factor(Year) + factor(District), data=Panel_Data_soc)
cluster9s<-cluster.vcov(Model9s,cbind(Panel_Data_soc$Year,Panel_Data_soc$District), force_posdef=TRUE)
CMREG9s <- coeftest(Model9s,cluster9s)

Model10s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime  +
                ntl_pretrend +
                Pop + DecayYr_additive_soc*Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                 NTL+
                factor(Year) + factor(District), data=Panel_Data_soc)
cluster10s<-cluster.vcov(Model10s,cbind(Panel_Data_soc$Year,Panel_Data_soc$District), force_posdef=TRUE)
CMREG10s <- coeftest(Model10s,cluster10s)

Model11s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime  +
                 ntl_pretrend +
                 Pop + DecayYr_additive_soc*Pop+
                 wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                 treat_minus1_soc + treat_minus2_soc + treat_minus3_soc+treat_minus4_soc+treat_minus5_soc+
                 factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster11s<-cluster.vcov(Model11s,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG11s <- coeftest(Model11s,cluster11s)

Model12s <- lm(Forest_Loss_additive ~ DecayYr_additive_soc + DecayAddControl_soc + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime  +
                ntl_pretrend+NTL+
                Pop + DecayYr_additive_soc*Pop+
                wdpapct_2007 + DecayYr_additive_soc*wdpapct_2007 +
                treat_minus1_soc + treat_minus2_soc + treat_minus3_soc+treat_minus4_soc+treat_minus5_soc+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster12s<-cluster.vcov(Model12s,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG12s <- coeftest(Model12s,cluster12s)

##Cumulative Forest Loss, Infra, Thresh=10, Distance Decay applied only for projects within 100km

Model109d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime  + 
                  ntl_pretrend  +
                  Pop + DecayYr100_additive*Pop +
                  wdpapct_2007 + DecayYr100_additive*wdpapct_2007+
                  factor(Year) + factor(District), data=Panel_Data)
cluster109d <- cluster.vcov(Model109d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG109d <- coeftest(Model109d, cluster109d)


##Cumulative Forest Loss, Infra, treatment is project count within 100km, Thresh=10

Model900 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime  + 
                 ntl_pretrend  +
                 Pop + ProjCnt100_additive*Pop +
                 wdpapct_2007 + ProjCnt100_additive*wdpapct_2007+
                 factor(Year) + factor(District), data=Panel_Data)
cluster900 <- cluster.vcov(Model900, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG900 <- coeftest(Model900, cluster900)

##Cumulative Forest Loss, Infra, Thresh=10, treatment proximity includes only projects in implementation or completion from TUFF dataset, not pipeline:commitment
#first load dataset with treatment values that exclude pipeline:commitment projects
Panel_Data_status<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add.csv")

Model9st <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend  +
               Pop + DecayYr_additive*Pop +
               wdpapct_2007 + DecayYr_additive*wdpapct_2007+
               factor(Year) + factor(District), data=Panel_Data_status)
cluster9st<-cluster.vcov(Model9st,cbind(Panel_Data_status$Year,Panel_Data_status$District), force_posdef=TRUE)
CMREG9st <- coeftest(Model9st,cluster9st)


##Cumulative Forest Loss, Infra, Thresh=15

Model159 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime  + 
                 ntl_pretrend  +
                 Pop + DecayYr_additive*Pop +
                 wdpapct_2007 + DecayYr_additive*wdpapct_2007+
                 factor(Year) + factor(District), data=Panel_Data_15)
cluster159 <- cluster.vcov(Model159, cbind(Panel_Data_15$Year, Panel_Data_15$District), force_posdef=TRUE)
CMREG159 <- coeftest(Model159, cluster159)


#-------------------------#
#Stargazer#
#-------------------------#

#WORKING PAPER SUMMARY STATS#
#rename vars directly in html file
#infra treatment stats table
stargazer(Panel_Data_infra, type = "html", nobs = TRUE, mean.sd = TRUE, median = TRUE,
          iqr = TRUE,
          keep=c("Forest_Loss_additive","DecayYr_additive","DecayYr100_additive","ProjCnt100_additive",
                 "DecayAddControl","DecayAddControl100"))
#infra covars stats table
stargazer(Panel_Data,type="html", 
          keep=c("MinTemp","MinPrecip","Max","Mean","Elevation","Slope","UrbTravTime","ntl_pretrend","NTL",
                 "Pop","Pct","PreLevelControl","PreTrendControl"))

#soc treatment stats table
stargazer(Panel_Data_soc, type = "html", nobs = TRUE, mean.sd = TRUE, median = TRUE,
          iqr = TRUE,
          keep=c("Forest_Loss_additive","DecayYr_additive_soc","DecayAddControl100_soc"))
#soc covars stats table
stargazer(Panel_Data_soc,type="html", 
          keep=c("MinTemp","MinPrecip","Max","Mean","Elevation","Slope","UrbTravTime","ntl_pretrend","NTL",
                 "Pop","Pct","PreLevelControl","PreTrendControl"))

#RESULTS TABLES#
stargazer(CMREG1, CMREG1.1, CMREG5,CMREG4,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","315,028","315,028","315,028","315,028"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes")),
          title="Tanzania Infra Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

#WORKING PAPER INFRA RESULTS TABLE#
stargazer(CMREG1, CMREG1.1, CMREG2, CMREG3,CMREG9,CMREG10,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),
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

#Working Paper Social Sector Results Table
stargazer(CMREG1s, CMREG1.1s, CMREG2s, CMREG3s,CMREG9s,CMREG10s,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),
          omit.stat=c("f","ser"),
          covariate.labels=c("Treatment (Proximity)","Proximity Control","Baseline NDVI","NDVI Pre-Trend",
                             "Elevation","Slope",
                             "Urban Travel Time","Nighttime Lights Pre-Trend","Population",
                             "Baseline Protected Areas",
                             "Nighttime Lights",
                             "Population*Treatment","Protected Area*Treatment"),
          add.lines=list(c("Observations","300,314","300,314","300,314","300,314","300,314","300,314","300,314"),
                         c("Climate Controls?","No","No","Yes","Yes","Yes","Yes"),
                         c("District Fixed Effects?","No","No","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes","Yes")),
          title="Tanzania Social Sector Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

#Working Paper Infra Robustness Checks
#remove pop, wdpa manually in html code (but interaction effects should remain)
stargazer(CMREG109d, CMREG900,CMREG9st, CMREG159, 
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip","Elevation","Slope","UrbTravTime","ntl",
                 "PreLevelControl","PreTrendControl","Constant"),
          omit.stat=c("f","ser"),
          covariate.labels=c("Treatment (Proximity, 100km)","100km Proximity Control",
                             "Treatment (100km Project Count)",
                             "Treatment (Proximity)","Proximity Control",
                             "Pop","wdpa",
                             "Population*Treatment (100km Proximity)","Protected Area*Treatment (100km Proximity)",
                             "Population*Treatment (100km Count)","Protected Area*Treatment (100km Count)",
                             "Population*Treatment (Proximity)","Protected Area*Treatment (Proximity)"),
          add.lines=list(c("Observations","315,028","315,028","315,028","261,870"),
                         c("Standing Forest Threshold","10%","10%","10%","15%")),
          title="Tanzania Infrastructure Regression Results: Robustness Checks",
          dep.var.labels=c("Cumulative Forest Loss"))

#Working Paper Infra Treatment Leads Table
stargazer(CMREG13,CMREG14,CMREG15,CMREG16,CMREG17,
          type="html",align=TRUE,
          keep=c("trtlead"),
          omit.stat=c("f","ser"),
          covariate.labels=c("Treatment Lead"),
          column.labels=c("1 Year","2 Year","3 Year","4 Year","5 Year"),
          add.lines=list(c("Observations","292,526","270,024","247,522","225,020","202,518")),
          omit.table.layout="d",
          title="Tanzania Infrastructure Treatment Leads")



