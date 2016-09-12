
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
library(pastecs)

##Panel Data, Thresh=10
Panel_Data <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia.csv")
Panel_Data_add <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_panel_data_add.csv")
Panel_Data<-Panel_Data_add

#Add Post-2003 indicator
Panel_Data$post_2003<-0
Panel_Data$post_2003[Panel_Data$Year>2003]<-1

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

Panel_Data_trtlag <- Panel_Data[Panel_Data$Year<=2009,]

# Panel_Data_test <- Panel_Data
# Panel_Data_test=Panel_Data_test[Panel_Data_test$MinYr!=0,]
# min <- aggregate(MinYr ~ ID, Panel_Data_test, function(x) min(x))
# colnames(min)[2] <- "MinDist"
# Panel_Data_min <- merge(Panel_Data, min, by.x="ID", by.y="ID")
# 
# mean<- aggregate(DecayYr_additive~Year, Panel_Data, function(x) mean(x))
# mean

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

# Model2 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
#                    MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
#                    UrbTravTime, data=Panel_Data)
# cluster2 <- cluster.vcov(Model2, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
# CMREG2 <- coeftest(Model2, cluster2)
# 
# Model3 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
#                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
#                UrbTravTime + factor(District), data=Panel_Data)
# cluster3 <- cluster.vcov(Model3, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
# CMREG3 <- coeftest(Model3, cluster3)

Model4 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               Pop + DecayYr_additive*Pop +
               factor(Year) + factor(District), data=Panel_Data)
cluster4 <- cluster.vcov(Model4, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4 <- coeftest(Model4, cluster4)

Model5 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + 
               PreLevelControl + PreTrendControl + MinTemp + MaxTemp + MeanTemp + 
               MaxPrecip + MeanPrecip + MinPrecip + 
               Elevation + Slope + UrbTravTime + Pop_2000 + Pop + DecayYr_additive*Pop +
               Year + Year*post_2003 + factor(District), data=Panel_Data)
cluster5 <- cluster.vcov(Model5, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG5 <- coeftest(Model5, cluster5)

Model6 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+ 
               factor(Year) + factor(District), data=Panel_Data)
cluster6<-cluster.vcov(Model6,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG6 <- coeftest(Model6,cluster6)

Model7 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+ 
               factor(Year) + factor(District), data=Panel_Data)
cluster7<-cluster.vcov(Model7,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG7 <- coeftest(Model7,cluster7)

Model8 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               plantation_pct + DecayYr_additive*plantation_pct + 
               factor(Year) + factor(District), data=Panel_Data)
cluster8<-cluster.vcov(Model8,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG8 <- coeftest(Model8,cluster8)

Model9 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               Pop + DecayYr_additive*Pop +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+
               plantation_pct + DecayYr_additive*plantation_pct + 
               factor(Year) + factor(District), data=Panel_Data)
cluster9<-cluster.vcov(Model9,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG9 <- coeftest(Model9,cluster9)

Model10 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + Pop_2000 + 
                ntl_pretrend + NTL_2003 + NTL +
                Pop + DecayYr_additive*Pop +
                wdpapct_2003 + DecayYr_additive*wdpapct_2003+
                concessionpct_2003 + DecayYr_additive*concessionpct_2003+
                plantation_pct + DecayYr_additive*plantation_pct + 
                factor(Year) + factor(District), data=Panel_Data)
cluster10<-cluster.vcov(Model10,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG10 <- coeftest(Model10,cluster10)

Model11 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + Pop_2000 + 
               ntl_pretrend + NTL_2003 +
               Pop + DecayYr_additive*Pop +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+
               plantation_pct + DecayYr_additive*plantation_pct + 
               treat_minus1 + treat_minus2+treat_minus3+treat_minus4+treat_minus5+
               factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster11<-cluster.vcov(Model11,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG11 <- coeftest(Model11,cluster11)

Model12 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + Pop_2000 + 
                ntl_pretrend + NTL_2003 + NTL+
                Pop + DecayYr_additive*Pop +
                wdpapct_2003 + DecayYr_additive*wdpapct_2003+
                concessionpct_2003 + DecayYr_additive*concessionpct_2003+
                plantation_pct + DecayYr_additive*plantation_pct + 
                treat_minus1 + treat_minus2+treat_minus3+treat_minus4+treat_minus5+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster12<-cluster.vcov(Model12,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG12 <- coeftest(Model12,cluster12)


##Cumulative Forest Loss, treatment is distance decay within 100km

Model104d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster104d <- cluster.vcov(Model104d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG104d <- coeftest(Model104d, cluster104d)

Model105d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime + Year +Year*post_2003 + factor(District), data=Panel_Data)
cluster105d <- cluster.vcov(Model105d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG105d <- coeftest(Model105d, cluster105d)

##Cumulative Forest Loss, treatment is distance decay within 25km

Model2504d <- lm(Forest_Loss_additive ~ DecayYr25_additive + DecayAddControl25 + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster2504d <- cluster.vcov(Model2504d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2504d <- coeftest(Model2504d, cluster2504d)

Model2505d <- lm(Forest_Loss_additive ~ DecayYr25_additive + DecayAddControl25 + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime + Year +Year*post_2003 + factor(District), data=Panel_Data)
cluster2505d <- cluster.vcov(Model2505d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2505d <- coeftest(Model2505d, cluster2505d)


##Cumulative Forest Loss, treatment is cumulative project count within 100km,Thresh=10

Model400 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + factor(Year) + factor(District), data=Panel_Data)
cluster400 <- cluster.vcov(Model400, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG400 <- coeftest(Model400, cluster400)

Model500 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + 
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

stargazer(Panel_Data,type="html", 
          keep=c("Forest_Loss_additive","DecayYr_additive","Control","MinTemp","MinPrecip",
                 "Max","Mean","Elevation","Slope","UrbTravTime","RivDist","RoadDist","ProjCnt100",
                 "DecayYr100_additive","DecayYr25_additive"))


stargazer(CMREG1, CMREG4, type="html", keep=c("Forest_Loss","additive","Control","Min","Max","Elevation","Slope","Year"))

stargazer(CMREG1, CMREG1.1, CMREG5, CMREG4,CMREG7,CMREG6,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),omit.labels=c("factor","Temp","Precip"),
          #keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","58,996","58,996","58,996","58,996"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes")),
          title="Cambodia Infra Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG4,CMREG6,CMREG7,CMREG8,CMREG9,CMREG10,CMREG11,CMREG12,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),omit.labels=c("factor","Temp","Precip"),
          omit.stat=c("f","ser"),
          # add.lines=list(c("Observations","58,996","58,996","58,996","58,996"),
          #                c("District Fixed Effects?","No","No","Yes","Yes"),
          #                c("Year Fixed Effects?","No","No","No","Yes")),
          title="Cambodia Infra Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG1, CMREG1.1, CMREG5, CMREG4,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Year","Elevation","Slope","UrbTravTime","Pop","Post"),
          covariate.labels=c("Treatment (Proximity)","Proximity Control","Baseline NDVI","NDVI Pre-Trend",
                             "Elevation","Slope",
                             "Urban Travel Time","Baseline Population","Population","Year"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","58,996","58,996","58,996","58,996"),
                         c("District Fixed Effects?","No","No","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes"),
                         c("Climate Controls?","No","No","Yes","Yes")),
          title="Cambodia Infrastructure Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

stargazer(CMREG500, CMREG400, CMREG105d, CMREG104d,CMREG2505d, CMREG2504d,
          type="html", align=TRUE,
          keep=c("Forest_Loss","additive","Control","Min","Max","Mean","Year","Elevation","Slope","UrbTravTime","Post"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","58,996","58,996","58,996","58,996","58,996","58,996"),
                         c("District Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","Yes","No","Yes","No","Yes")),
          title="Cambodia Infra Regression Results: Alternate Treatments",
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

