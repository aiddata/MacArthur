
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
#Panel_Data <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia.csv")
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



# Panel_Data_test <- Panel_Data
# Panel_Data_test=Panel_Data_test[Panel_Data_test$MinYr!=0,]
# min <- aggregate(MinYr ~ ID, Panel_Data_test, function(x) min(x))
# colnames(min)[2] <- "MinDist"
# Panel_Data_min <- merge(Panel_Data, min, by.x="ID", by.y="ID")
# 
# mean<- aggregate(DecayYr_additive~Year, Panel_Data, function(x) mean(x))
# mean

##Panel Data, Thresh=5

Panel_Data_5 <- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_panel_data_add_thresh5.csv")

#Add Post-2003 indicator
Panel_Data_5$post_2003<-0
Panel_Data_5$post_2003[Panel_Data_5$Year>2003]<-1

##Panel Data, Thresh=15

Panel_Data_15<- read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/cambodia_panel_data_add_thresh15.csv")
# ntl<-Panel_Data[,c("ID","Year","wdpapct_2000","wdpapct_2003","concessionpct_all","concessionpct_2003",
#                    "plantation_pct","Pop_2000","gpw_v4_density.2005.mean","gpw_v4_density.2010.mean","gpw_v4_density.2015.mean",
#                    "Pop","ntltrend_0913","ntl_pretrend","NTL","NTL_2003")]

#Add Post-2003 indicator
Panel_Data_15$post_2003<-0
Panel_Data_15$post_2003[Panel_Data_15$Year>2003]<-1

##Panel Data, Thresh=10, only cells with Chinese activities within 100km
Panel_Data_100<-Panel_Data[Panel_Data$DecayAddControl100!=0,]

## Panel_Data, Thresh=10, only cells with Chinese activities within 25km
Panel_Data_25<-Panel_Data[Panel_Data$DecayAddControl25!=0,]


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
               UrbTravTime + 
               ntl_pretrend+
               Pop+
               factor(District), data=Panel_Data)
cluster2 <- cluster.vcov(Model2, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2 <- coeftest(Model2, cluster2)

Model3 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + 
                ntl_pretrend+
                Pop+
                factor(Year) + factor(District), data=Panel_Data)
cluster3 <- cluster.vcov(Model3, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG3 <- coeftest(Model3, cluster3)

Model4 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend  +
               Pop + DecayYr_additive*Pop +
               factor(Year) + factor(District), data=Panel_Data)
cluster4 <- cluster.vcov(Model4, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG4 <- coeftest(Model4, cluster4)

Model6 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+ 
               factor(Year) + factor(District), data=Panel_Data)
cluster6<-cluster.vcov(Model6,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG6 <- coeftest(Model6,cluster6)

Model7 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend  +
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+ 
               factor(Year) + factor(District), data=Panel_Data)
cluster7<-cluster.vcov(Model7,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG7 <- coeftest(Model7,cluster7)

Model8 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
               plantation_pct + DecayYr_additive*plantation_pct + 
               factor(Year) + factor(District), data=Panel_Data)
cluster8<-cluster.vcov(Model8,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG8 <- coeftest(Model8,cluster8)

#this is the main model to use for robustness checks!
Model9 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
               Pop + DecayYr_additive*Pop +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+
               plantation_pct + DecayYr_additive*plantation_pct + 
               factor(Year) + factor(District), data=Panel_Data)
cluster9<-cluster.vcov(Model9,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG9 <- coeftest(Model9,cluster9)

Model10 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
               Pop + DecayYr_additive*Pop +
               wdpapct_2003 + DecayYr_additive*wdpapct_2003+
               concessionpct_2003 + DecayYr_additive*concessionpct_2003+
               plantation_pct + DecayYr_additive*plantation_pct +
               NTL+
               factor(Year) + factor(District), data=Panel_Data)
cluster10<-cluster.vcov(Model10,cbind(Panel_Data$Year,Panel_Data$District), force_posdef=TRUE)
CMREG10 <- coeftest(Model10,cluster10)

#Cumulative Forest Loss, add in treatment lead

Model11 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime  + 
               ntl_pretrend +
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
                ntl_pretrend + NTL+
                Pop + DecayYr_additive*Pop +
                wdpapct_2003 + DecayYr_additive*wdpapct_2003+
                concessionpct_2003 + DecayYr_additive*concessionpct_2003+
                plantation_pct + DecayYr_additive*plantation_pct + 
                treat_minus1 + treat_minus2+treat_minus3+treat_minus4+treat_minus5+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster12<-cluster.vcov(Model12,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG12 <- coeftest(Model12,cluster12)

Model13 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2003 + concessionpct_2003 + plantation_pct+
                factor(Year) + factor(District), data=Panel_Data_trtlag13)
cluster13<-cluster.vcov(Model13,cbind(Panel_Data_trtlag13$Year,Panel_Data_trtlag13$District), force_posdef=TRUE)
CMREG13 <- coeftest(Model13,cluster13)

Model14 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl+PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2003 + concessionpct_2003 + plantation_pct+
                factor(Year) + factor(District), data=Panel_Data_trtlag12)
cluster14<-cluster.vcov(Model14,cbind(Panel_Data_trtlag12$Year,Panel_Data_trtlag12$District), force_posdef=TRUE)
CMREG14 <- coeftest(Model14,cluster14)

Model15 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl+PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2003 + concessionpct_2003 + plantation_pct+
                factor(Year) + factor(District), data=Panel_Data_trtlag11)
cluster15<-cluster.vcov(Model15,cbind(Panel_Data_trtlag11$Year,Panel_Data_trtlag11$District), force_posdef=TRUE)
CMREG15 <- coeftest(Model15,cluster15)

Model16 <- lm(Forest_Loss_additive ~ trtlead+DecayAddControl+ PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2003 + concessionpct_2003 + plantation_pct+
                factor(Year) + factor(District), data=Panel_Data_trtlag10)
cluster16<-cluster.vcov(Model16,cbind(Panel_Data_trtlag10$Year,Panel_Data_trtlag10$District), force_posdef=TRUE)
CMREG16 <- coeftest(Model16,cluster16)

Model17 <- lm(Forest_Loss_additive ~ trtlead+ DecayAddControl +PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime +  
                ntl_pretrend +
                Pop + wdpapct_2003 + concessionpct_2003 + plantation_pct+
                factor(Year) + factor(District), data=Panel_Data_trtlag)
cluster17<-cluster.vcov(Model17,cbind(Panel_Data_trtlag$Year,Panel_Data_trtlag$District), force_posdef=TRUE)
CMREG17 <- coeftest(Model17,cluster17)


##Cumulative Forest Loss, treatment is distance decay within 100km

Model109d <- lm(Forest_Loss_additive ~ DecayYr100_additive + DecayAddControl100 + PreLevelControl + PreTrendControl + 
                  MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                  UrbTravTime + 
                  ntl_pretrend +
                  Pop + DecayYr100_additive*Pop +
                  wdpapct_2003 + DecayYr100_additive*wdpapct_2003+
                  concessionpct_2003 + DecayYr100_additive*concessionpct_2003+
                  plantation_pct + DecayYr100_additive*plantation_pct + 
                  factor(Year) + factor(District), data=Panel_Data)
cluster109d <- cluster.vcov(Model109d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG109d <- coeftest(Model109d, cluster109d)

##Cumulative Forest Loss, treatment is distance decay within 25km

Model2509d <- lm(Forest_Loss_additive ~ DecayYr25_additive + DecayAddControl25 + PreLevelControl + PreTrendControl + 
                   MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                   UrbTravTime + 
                   ntl_pretrend +
                   Pop + DecayYr25_additive*Pop +
                   wdpapct_2003 + DecayYr25_additive*wdpapct_2003+
                   concessionpct_2003 + DecayYr25_additive*concessionpct_2003+
                   plantation_pct + DecayYr25_additive*plantation_pct + 
                   factor(Year) + factor(District), data=Panel_Data)
cluster2509d <- cluster.vcov(Model2509d, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG2509d <- coeftest(Model2509d, cluster2509d)


##Cumulative Forest Loss, treatment is cumulative project count within 100km,Thresh=10

Model900 <- lm(Forest_Loss_additive ~ ProjCnt100_additive + PreLevelControl + PreTrendControl + 
               MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
               UrbTravTime + 
               ntl_pretrend +
                 Pop + ProjCnt100_additive*Pop +
                 wdpapct_2003 + ProjCnt100_additive*wdpapct_2003+
                 concessionpct_2003 + ProjCnt100_additive*concessionpct_2003+
                 plantation_pct + ProjCnt100_additive*plantation_pct + 
               factor(Year) + factor(District), data=Panel_Data)
cluster900 <- cluster.vcov(Model900, cbind(Panel_Data$Year, Panel_Data$District), force_posdef=TRUE)
CMREG900 <- coeftest(Model900, cluster900)

##Cumulative Forest Loss, Thresh=5

Model59 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                UrbTravTime + 
                ntl_pretrend +
                Pop + DecayYr_additive*Pop +
                wdpapct_2003 + DecayYr_additive*wdpapct_2003+
                concessionpct_2003 + DecayYr_additive*concessionpct_2003+
                plantation_pct + DecayYr_additive*plantation_pct + 
                factor(Year) + factor(District), data=Panel_Data_5)
cluster59 <- cluster.vcov(Model59, cbind(Panel_Data_5$Year, Panel_Data_5$District), force_posdef=TRUE)
CMREG59 <- coeftest(Model59, cluster59)

##Cumulative Forest Loss, Thresh=15

Model159 <- lm(Forest_Loss_additive ~ DecayYr_additive + DecayAddControl + PreLevelControl + PreTrendControl + 
                 MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + Elevation + Slope + 
                 UrbTravTime + 
                 ntl_pretrend +
                 Pop + DecayYr_additive*Pop +
                 wdpapct_2003 + DecayYr_additive*wdpapct_2003+
                 concessionpct_2003 + DecayYr_additive*concessionpct_2003+
                 plantation_pct + DecayYr_additive*plantation_pct + 
                 factor(Year) + factor(District), data=Panel_Data_15)
cluster159 <- cluster.vcov(Model159, cbind(Panel_Data_15$Year, Panel_Data_15$District), force_posdef=TRUE)
CMREG159 <- coeftest(Model159, cluster159)



#-------------------------#
#Stargazer#
#-------------------------#

#WORKING PAPER SUMMARY STATS#
#rename vars directly in html file
stargazer(Panel_Data,type="html", 
          keep=c("MinTemp","MinPrecip","Max","Mean","Elevation","Slope","UrbTravTime","ntl_pretrend","NTL",
                 "Pop","Pct","PreLevelControl","PreTrendControl"))


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

#WORKING PAPER MAIN TABLE#
stargazer(CMREG1, CMREG1.1, CMREG2,CMREG3, CMREG9, CMREG10,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip"),
          #omit.labels=c("factor","Temp","Precip"),
          covariate.labels=c("Treatment (Proximity)","Proximity Control","Baseline NDVI","NDVI Pre-Trend",
                             "Elevation","Slope",
                             "Urban Travel Time","Nighttime Lights Pre-Trend","Population",
                             "Baseline Protected Areas","Baseline Concessions","Plantations",
                             "Nighttime Lights",
                             "Population*Treatment","Protected Area*Treatment",
                             "Concession*Treatment","Plantation*Treatment"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","58,996","58,996","58,996","58,996","58,996","58,996"),
                         c("Climate Controls?","No","No","Yes","Yes","Yes","Yes"),
                         c("District Fixed Effects?","No","No","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes","Yes","Yes")),
          title="Cambodia Infrastructure Regression Results",
          dep.var.labels=c("Cumulative Forest Loss"))

#WORKING PAPER ROBUSTNESS CHECKS#
stargazer(CMREG109d,CMREG2509d,CMREG900,CMREG59, CMREG159,
          type="html", align=TRUE,
          omit=c("factor","Temp","Precip","Elevation","Slope","UrbTravTime","ntl",
                 "PreLevelControl","PreTrendControl","Constant"),
          #order=c(5,6,1,2,3,4,7),
          covariate.labels=c("Treatment (Proximity, 100km)","100km Proximity Control",
                             "Treatment (Proximity, 25km)","25km Proximity Control",
                             "Treatment (100km Project Count)",
                             "Treatment (Proximity)","Proximity Control",
                             "Pop","wdpa","conc","plant",
                             "Population*Treatment (100km Proximity)","Protected Area*Treatment (100km Proximity)","Concession*Treatment (100km Proximity)","Plantation*Treatment (100km Proximity)",
                             "Population*Treatment (25km Proximity)","Protected Area*Treatment (25km Proximity)","Concession*Treatment (25km Proximity)","Plantation*Treatment (25km Proximity)",
                             "Population*Treatment (100km Count)","Protected Area*Treatment (100km Count)","Concession*Treatment (100km Count)","Plantation*Treatment (100km Count)",
                             "Population*Treatment (Proximity)","Protected Area*Treatment (Proximity)","Concession*Treatment (Proximity)","Plantation*Treatment (Proximity)"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","58,996","58,996","58,996","64,106","54,740"),
                         c("Standing Forest Threshold","10%","10%","10%","5%","15%")),
          title="Cambodia Infrastructure Regression Results: Robustness Checks",
          dep.var.labels=c("Cumulative Forest Loss"))

#Working Paper Infra Treatment Leads Table
stargazer(CMREG13,CMREG14,CMREG15,CMREG16,CMREG17,
          type="html",align=TRUE,
          keep=c("trtlead"),
          omit.stat=c("f","ser"),
          covariate.labels=c("Treatment Lead"),
          column.labels=c("1 Year","2 Year","3 Year","4 Year","5 Year"),
          add.lines=list(c("Observations","54,782","50,568","46,354","42,140","37,926")),
          title="Cambodia Infrastructure Treatment Leads")



#Scratch

panel_sub<-Panel_Data[,c("ID","Year","DecayYr_additive")]
panel_merge<-merge(Panel_Data_5,panel_sub,by=c("ID","Year"),all=T)
panel_merge$trtdiff<-NA
panel_merge$trtdiff<-panel_merge$DecayYr_additive.x-panel_merge$DecayYr_additive.y

panel_merge_117204 <- panel_merge[panel_merge$ID=="117204",]
