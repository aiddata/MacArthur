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

#---------------------------------------------------#
#Settings
#---------------------------------------------------#

forest_thresh = 15
restrict_analysis = FALSE

#---------------------------------------------------#
#Download Data for Analysis - This step may take a while (up to hours)
#---------------------------------------------------#
source("RDownload.R")
mDir = getwd()
active_dir_path  <- downlad_data(mDir)

#---------------------------------------------------#
#Load the dataframes in for analysis after download.
#---------------------------------------------------#
csv <- paste(active_dir_path, "/extracts/africa.csv", sep="")
json <- paste(active_dir_path, "/extracts/africa.json", sep="")
dta <- read.csv(csv)
vars <- fromJSON(txt=json)

#--------------------------------------------------#
#Subset the Cell Dataframe 
#--------------------------------------------------#
dta2 <- dta[dta$NAME_0 == "Tanzania",]

#--------------------------------------------------#
#Load and Subset the Cell Spatial Datframe
#--------------------------------------------------#
spdf_cells <- paste(active_dir_path, "/grids/africa_grid.shp", sep="")
cells <- readShapePoly(spdf_cells)
#Keep only relevant cells
AOI_cells <- sp::merge(cells, dta2, by="ID", all.x=FALSE)
if(restrict_analysis != FALSE)
{
  AOI_cells<- AOI_cells[1:restrict_analysis,]
}
proj4string(AOI_cells) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


#--------------------------------------------------#
#Load and Subset the Spatial ADM zone data
#--------------------------------------------------#
spdf_adm_path <- paste(active_dir_path, "/ADM2/GADM_MacEcohotspotSubset_ADM2.shp", sep="")
spdf_adm <- readShapePoly(spdf_adm_path)
#Keep only relevant ADM data
spdf_adm <- spdf_adm[spdf_adm@data$NAME_0 == "Tanzania",]

#--------------------------------------------------#
#Load and Subset the MacArthur Aid Data
#--------------------------------------------------#
location_csv <- paste(active_dir_path, "/MacArthur_Geocoded_data/locations.csv", sep="")
locations <- read.csv(location_csv)
locations2 <- locations[grep("Tanzania", locations$gazetteer_adm_name),]
coords = cbind(locations2$longitude, locations2$latitude)

#--------------------------------------------------#
#Create a Spatial Dataframe of the MacArthur Aid Data
#--------------------------------------------------#
Mac_spdf <- SpatialPointsDataFrame(coords, locations2)
proj4string(Mac_spdf) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
#Join in the year of aid allocation from projects_ancillary
Mac_years <- paste(active_dir_path, "/MacArthur_Geocoded_data/projects_ancillary.csv", sep="")
Mac_yr_dta <- read.csv(Mac_years)
Mac_spdf <- merge(Mac_spdf, Mac_yr_dta, by="project_id")
#Subset by precision code = 1 or 2
Mac_prec <- Mac_spdf[Mac_spdf@data$precision_code<=2,]
Mac_spdf <- Mac_prec
#Subset by sector code (combined infrastructure projects, no sector code "160")
Mac_sector <- Mac_spdf[Mac_spdf@data$crs_sector_code%in%c("210","220","230","320"),]
Mac_spdf <- Mac_sector
#Subset by status = implementation or completion (not pipeline)
Mac_status <- Mac_spdf[Mac_spdf@data$status_code%in%c("2","3"),]
Mac_spdf <- Mac_status
# #Take out sector 160 projects that aren't building infrastructure
# Mac_160<- Mac_spdf[Mac_spdf@data$project_id!="33053",]
# Mac_160<- Mac_160[Mac_160@data$project_id!="1919",]
# Mac_spdf <- Mac_160

#write.csv(Mac_spdf@data,"/home/aiddata/Desktop/Github/MacArthur/modelData/Mac_spdf_TanzaniaInfra_Thresh15.csv")
#writePointsShape(Mac_spdf, "/home/aiddata/Desktop/Github/MacArthur/modelData/Mac_spdf_TanzaniaInfra_Thresh15.shp")

#--------------------------------------------------#
#Create threshdolded forest datasets
#--------------------------------------------------#
ndviDTA <- AOI_cells#dta2[c("tc00_e", "lnyx_2000e")]
ndviDTA_for <- ndviDTA[ndviDTA$tc00_e >= forest_thresh,]
ndviDTA_notfor <- ndviDTA[ndviDTA$tc00_e <= forest_thresh,]
AOI_cells = ndviDTA_for

#--------------------------------------------------#
#Calculate the distances between cells and MacArthur points
#--------------------------------------------------#
source("RDist.R")
dMatrix <- RDist(AOI_cells, Mac_spdf)
#In dMatrix, every column is a cell (referenced in order to AOI_cells)
#Every row is a MacArthur project (referenced in order to Mac_spdf)

#Average distance in KM:
avgDistKm <- mean(dMatrix) / 1000

#Minimum
col_mins <- do.call(pmin, lapply(1:nrow(dMatrix), function(i)dMatrix[i,]))
minDistKm <- mean(col_mins) / 1000

#--------------------------------------------------#
#Calculate the correlogram
#--------------------------------------------------#
#correlogram_data <- correlog(x = coordinates(AOI_cells)[,1], y = coordinates(AOI_cells)[,2], z=AOI_cells$lnyx_1999e, increment=5, latlon=TRUE, na.rm=TRUE, resamp=5)

#save (correlogram_data, file="/home/aiddata/Desktop/Github/MacArthur/modelData/tanz_correl_Thresh15.RData")

load("/home/aiddata/Desktop/Github/MacArthur/modelData/tanz_correl_Thresh15.RData")

#save data into a function to calculate the distance-decay penalty later.
#Chinese projects are "weighted" according to their distance.
#The absolulute correlation for a given distance is used as a weight.
#Projects at distances with a higher positive or inverse correlation are given the highest
#weights.

dVals <- abs(correlogram_data$mean.of.class)
cVals <- abs(correlogram_data$correlation)

#--------------------------------------------------#
#Calculate the average spatial decays for visualization
#--------------------------------------------------#
#no weights, all before the first time X = 0 are counted.
thresh_dMatrix <- dMatrix
thresh_dMatrix[thresh_dMatrix > (as.numeric(correlogram_data$x.intercept)*1000)] <- NA
total_distance_km <- colSums(thresh_dMatrix, na.rm=TRUE) / 1000

AOI_cells$thresh_tot_proj <- apply(thresh_dMatrix, 2, function(x) length(which(!is.na(x))))
AOI_cells$thresh_totDist <- total_distance_km
AOI_cells$thresh_avgDist <- AOI_cells$thresh_totDist / AOI_cells$thresh_tot_proj 


#distance decay
decay_dMatrix <- dMatrix
decay_dMatrix_adj <- apply(decay_dMatrix, 1:2, function(x){(cVals[which.min(abs(dVals - x))])[[1]]})
AOI_cells$thresh_weightedDist <- colSums(decay_dMatrix_adj) / 1000

#Drop cells which have no projects within the threshold
AOI_cells2 <- AOI_cells[!(is.na(AOI_cells@data$thresh_avgDist)),]
AOI_cellsBack <- AOI_cells
AOI_cells <- AOI_cells2

#writePolyShape(AOI_cells, "/home/aiddata/Desktop/Github/MacArthur/modelData/AOI_cells_Tanzania_Thresh15.shp")

#--------------------------------------------------#
#Calculate the over-time treatment effects
#--------------------------------------------------#
all_years <- unique(Mac_spdf$year)
all_years <- all_years[!is.na(all_years)]
record_length <- c(2001:2014)

dYears <- list()
#Drop all MacArthur projects that have no start or end date.
Mac_spdf <- Mac_spdf[!is.na(Mac_spdf$year),]

for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    year <- record_length[[years]]
    ThisYearMac <- Mac_spdf[Mac_spdf@data$year == year,]
    dYears[[years]] <- RDist(AOI_cells, ThisYearMac)
  }
  else
  {
    dYears[[years]] <- 0
  }
}

#In dMatrix, every column is a cell (referenced in order to AOI_cells)
#Every row is a MacArthur project for that year (referenced in order to Mac_spdf)

AvgYears <- vector()
for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    AvgYears[[years]] <- mean(dYears[[years]]) / 1000
  }
  else
  {
    AvgYears[[years]] <- 0
  }
}

Avg_MinYears <- vector()
for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    col_mins_year <- do.call(pmin, lapply(1:nrow(dYears[[years]]), function(i)dYears[[years]][i,]))
    nameRef <- paste("MinYr_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- col_mins_year / 1000
    Avg_MinYears[[years]] <- mean(col_mins_year) / 1000
  }
  else
  {
    Avg_MinYears[[years]] <- 0
    nameRef <- paste("MinYr_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- 0
  }
}

Dist_Decay_Yrs <- vector()
dvz <- cVals
dvz[dVals > correlogram_data$x.intercept] <- 0
for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    t_dyears <- dYears[[years]] / 1000
    decay_dMatrix_adj <- apply(t_dyears, 1:2, function(x){(dvz[which.min(abs(dVals - x))])[[1]]})
    nameRef <- paste("DecayYr_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- colSums(decay_dMatrix_adj) 
    Dist_Decay_Yrs[[years]] <- sum(colSums(decay_dMatrix_adj))
  }
  else
  {
    Dist_Decay_Yrs[[years]] <- 0
    nameRef <- paste("DecayYr_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- 0
  }
}


#---------------------------
#PROJECT YEARS START COUNT (count of all projects within 100km)
#---------------------------
Proj_Thresh_Count_Yrs <- vector()
#in KM
thresh <- 100
cthreshVals <- cVals

cthreshVals[dVals <= thresh] <- 1
cthreshVals[dVals > thresh] <- 0

for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    t_dyears <- dYears[[years]] / 1000
    decay_dMatrix_adj <- apply(t_dyears, 1:2, function(x){(cthreshVals[which.min(abs(dVals - x))])[[1]]})
    nameRef <- paste("Proj_Thresh_Count_Yrs_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- colSums(decay_dMatrix_adj)
    Proj_Thresh_Count_Yrs[[years]] <- sum(colSums(decay_dMatrix_adj))
  }
  else
  {
    Proj_Thresh_Count_Yrs[[years]] <- 0
    nameRef <- paste("Proj_Thresh_Count_Yrs_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- 0
  }
}
#---------------------------
#Limit distance decay threshold to 100km (rather than the x-intercept)
#---------------------------

DistDecay100 <- vector()
#in KM
thresh <- 100
cthreshVals_decay100 <- cVals
cthreshVals_decay100[dVals > thresh] <- 0

for(years in 1:length(record_length))
{
  if(record_length[[years]] %in% all_years)
  {
    t_dyears <- dYears[[years]] / 1000
    decay_dMatrix_adj <- apply(t_dyears, 1:2, function(x){(cthreshVals_decay100[which.min(abs(dVals - x))])[[1]]})
    nameRef <- paste("DistDecay100_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- colSums(decay_dMatrix_adj)
    DistDecay100[[years]] <- sum(colSums(decay_dMatrix_adj))
  }
  else
  {
    DistDecay100[[years]] <- 0
    nameRef <- paste("DistDecay100_",record_length[[years]], sep="")
    AOI_cells@data[nameRef] <- 0
  }
}


# CountProj_Years <- vector()
# for(years in 1:length(record_length))
# {
#   if(record_length[[years]] %in% all_years)
#   {
#   CountProj_Years[[years]] <- nrow(dYears[[years]])
#   }
#   else
#   {
#     CountProj_Years[[years]] <- 0
#   }
# }
# 
# #Build a quick temporal dataframe for plotting and ordering
# TempDF <- cbind.data.frame(record_length, AvgYears, Avg_MinYears, CountProj_Years, Dist_Decay_Yrs, Proj_Thresh_Count_Yrs)
# TempDF <- TempDF[with(TempDF, order(TempDF[,1])),]




#--------------------------------------------------#
#Pre-processing for analysis
#--------------------------------------------------#
DFa <- AOI_cells@data
#Drop irrelevant variables:
dropvars <- c("XMIN","XMAX","YMIN","YMAX","OBJECTID","ID_0","ISO","NAME_0","HASC_2","ID_1","NAME_1","NAME_2",
              "CCN_2","CCA_2","TYPE_2","ENGTYPE_2","NL_NAME_2","VARNAME_2","Shape_Leng","Shape_Area", 
              "thresh_tot_proj","thresh_totDist","thresh_avgDist","thresh_weightedDist")

DFa <- DFa[,!(names(DFa) %in% dropvars)]
DFa_hist <- DFa
DFa <- DFa[, -grep("(19)", names(DFa))]

#Prep for wide to long translation
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

rename_header <- function(x,sub)
{
  t <- paste(substr(x, 1, 0), sub, substr(x, 1, nchar(x)), sep = "")
  substr(t, 1, nchar(t)-1)
}

#names(DFa)[25:143]<- sapply(names(DFa)[25:143], function(x) {rename_header(x, substrRight(x,1))})

#names(DFa)[6:19] <- sapply(names(DFa)[6:19], function(x){substr(x, 5, nchar(x))})

#Reorder variables

for (i in 2:length(DFa)) {
  
  if (substr(colnames(DFa)[i], 1, 4) == "at41"){
    
    name = substr(colnames(DFa)[i],1, 4)
    year = substr(colnames(DFa)[i], 6, 9)
    letter = substr(colnames(DFa)[i], 10,10)
    dt = paste(letter,name,"_",year,sep="")
    colnames(DFa)[i] <- dt
  }
}

for (i in 2:length(DFa)) {
  
  if (substr(colnames(DFa)[i], 1, 4) == "pc41"){
    
    name = substr(colnames(DFa)[i],1, 4)
    year = substr(colnames(DFa)[i], 6, 9)
    letter = substr(colnames(DFa)[i], 10,10)
    dt = paste(letter,name,"_",year,sep="")
    colnames(DFa)[i] <- dt
  }
}

for (i in 2:length(DFa)) {
  
  if (substr(colnames(DFa)[i], 1, 4) == "ncc4"){
    
    name = substr(colnames(DFa)[i],1, 4)
    year = substr(colnames(DFa)[i], 6, 9)
    letter = substr(colnames(DFa)[i], 10,10)
    dt = paste(letter,name,"_",year,sep="")
    colnames(DFa)[i] <- dt
  }
}

for (i in 2:length(DFa)) {
  
  if (substr(colnames(DFa)[i], 1, 4) == "lnyx"){
    
    name = substr(colnames(DFa)[i],1, 4)
    year = substr(colnames(DFa)[i], 6, 9)
    letter = substr(colnames(DFa)[i], 10,10)
    dt = paste(letter,name,"_",year,sep="")
    colnames(DFa)[i] <- dt
  }
}

for (i in 2:length(DFa))
{
  colnames(DFa)[i] <- sub("per_loss_","loss_",colnames(DFa)[i])
}

for (i in 2:length(DFa))
{
  colnames(DFa)[i] <- sub("Proj_Thresh_Count_Yrs","ProjCnt100",colnames(DFa)[i])
}


#drop data for the year 2000, 2013-2015 (NTL ends in 2012)
DFa2 <- DFa[, -grep("(2000)", names(DFa))]
# DFa3 <- DFa2[, -grep("(2013)", names(DFa2))]
#DFa <- DFa[, -grep("(2014)", names(DFa3))]

#DFa <- DFa4

DFa3 <- DFa2[, -grep("^encc4", names(DFa2))]

PCloss <- grep("^loss", names(DFa3))
mean_ln <- grep("^elnyx", names(DFa3))
NTL <- grep("^encc4", names(DFa3))
minairTemp <- grep("^mat41", names(DFa3))
maxairTemp <- grep("^xat41", names(DFa3))
meanairTemp <- grep("^eat41", names(DFa3))
minPre <- grep("^mpc41", names(DFa3))
maxPre <- grep("^xpc41", names(DFa3))
meanPre <- grep("^epc41", names(DFa3))
MinDist <- grep("^MinYr", names(DFa3))
DecayDist <- grep("^DecayYr", names(DFa3))
ProjCount <- grep("^ProjCnt100", names(DFa3))
DecayDist100<- grep("^DistDecay100", names(DFa3))

#--------------------------------------------------#
#Selection of temporally-varying variables and shift from wide- to long-form
#--------------------------------------------------#

all_reshape <- c(PCloss, mean_ln, minairTemp, maxairTemp, meanairTemp, minPre, maxPre, meanPre, MinDist, DecayDist, ProjCount, DecayDist100)
DFa4 <- reshape(DFa3, varying=all_reshape,direction="long", idvar="ID", sep="_", timevar="Year")

DFa <- DFa4

#Rename names to something interpretable...
names(DFa)[names(DFa) == "ID_2"] = "District"
names(DFa)[names(DFa) == "loss"] = "Forest_Loss"
names(DFa)[names(DFa) == "encc4"] = "NighttimeLights"
names(DFa)[names(DFa) == "mat41"] = "MinTemp"
names(DFa)[names(DFa) == "xat41"] = "MaxTemp"
names(DFa)[names(DFa) == "eat41"] = "MeanTemp"
names(DFa)[names(DFa) == "mpc41"] = "MinPrecip"
names(DFa)[names(DFa) == "xpc41"] = "MaxPrecip"
names(DFa)[names(DFa) == "epc41"] = "MeanPrecip"
names(DFa)[names(DFa) == "selv_e"] = "Elevation"
names(DFa)[names(DFa) == "sslp_e"] = "Slope"
names(DFa)[names(DFa) == "dari_e"] = "RivDist"
names(DFa)[names(DFa) == "droa_e"] = "RoadDist"
names(DFa)[names(DFa) == "am50_e"] = "UrbTravTime"
names(DFa)[names(DFa) == "DistDecay100"] = "DecayYr100"


#--------------------------------------------------#
#Additive Year-on-Year 
#--------------------------------------------------#
Panel_Data <- DFa
#Panel_Data$MinYr_additive <- NA
Panel_Data$DecayYr_additive <- NA
Panel_Data$Forest_Loss_additive <- NA
Panel_Data$DecayYr100_additive <- NA
Panel_Data$ProjCnt100_additive <- NA

Panel_Data <- Panel_Data[with(Panel_Data, order(Panel_Data["ID"], Panel_Data["Year"])),]

calc_add <- function(Fdta, year, ID, var)
{
  a.dta <- Fdta[Fdta$ID == ID,]
  b.dta <- a.dta[a.dta$Year <= year,]
  
  exec_st <- paste("sum(b.dta$",var,")",sep="")
  return(eval(parse(text=exec_st)))
  
}

for(i in 1:length(Panel_Data[[1]]))
{
  Panel_Data["DecayYr_additive"][i,] <- calc_add(Panel_Data, Panel_Data[i,]["Year"][[1]], Panel_Data[i,]["ID"][[1]], "DecayYr")
  Panel_Data["DecayYr100_additive"][i,] <- calc_add(Panel_Data, Panel_Data[i,]["Year"][[1]], Panel_Data[i,]["ID"][[1]], "DecayYr100")
  #Panel_Data["MinYr_additive"][i,] <- calc_add(Panel_Data, Panel_Data[i,]["Year"][[1]], Panel_Data[i,]["ID"][[1]], "MinYr")
  Panel_Data["Forest_Loss_additive"][i,] <- calc_add(Panel_Data, Panel_Data[i,]["Year"][[1]], Panel_Data[i,]["ID"][[1]], "Forest_Loss")
  Panel_Data["ProjCnt100_additive"][i,] <- calc_add(Panel_Data, Panel_Data[i,]["Year"][[1]], Panel_Data[i,]["ID"][[1]], "ProjCnt100")
}

#Control Variables

val_lookup <- function(dta, var, id)
{
  ret_find_exec <- paste("dta$",var)
}

pre_trend_func <- function(dta, id)
{
  #build a dataframe
  NDVI_reshape <- c("lnyx_1990e", "lnyx_1991e", "lnyx_1992e", "lnyx_1993e", "lnyx_1994e", "lnyx_1995e", "lnyx_1996e", "lnyx_1997e", "lnyx_1998e", "lnyx_1999e")
  DFa_mdl <- DFa_hist[DFa_hist$ID == id,]
  mdl_dta <- reshape(DFa_mdl, varying=NDVI_reshape,direction="long", idvar="ID", sep="_", timevar="Year")
  return(lm(lnyx ~ Year, data=mdl_dta)$coefficients["Year"][[1]])
}

Panel_Data["DecayAddControl"] <- NA
Panel_Data["DecayAddControl100"]<- NA
Panel_Data["PreLevelControl"] <- NA
Panel_Data["PreTrendControl"] <- NA
for(i in 1:length(Panel_Data[[1]]))
{
  Panel_Data["DecayAddControl"][i,] <- calc_add(Panel_Data, max(Panel_Data["Year"][[1]]), Panel_Data[i,]["ID"][[1]], "DecayYr")
  Panel_Data["DecayAddControl100"][i,] <- calc_add(Panel_Data, max(Panel_Data["Year"][[1]]), Panel_Data[i,]["ID"][[1]], "DecayYr100")
  Panel_Data["PreLevelControl"][i,] <- AOI_cells@data$lnyx_1999e[AOI_cells@data$ID == Panel_Data[i,]["ID"][[1]]]
  Panel_Data["PreTrendControl"][i,] <- pre_trend_func(AOI_cells@data, Panel_Data[i,]["ID"][[1]])
}

write.csv(Panel_Data,"/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_Thresh15_JUNE.csv")
#-------------------------------------------------
#Add in additional data directly into panel dataset from AFRcells#
#Population, Baseline Protected Areas, Nighttime Lights#
#-------------------------------------------------

##Add in baseline protected area data, for pre-2001 and for pre-2008 since no projects in Tanzania until 2008
pa_2000 <- read.csv("/home/aiddata/Desktop/Github/MacArthur/ProtectedAreas_Data/merge_africa_grid_pre2001.csv")
pa_2007 <- read.csv("/home/aiddata/Desktop/Github/MacArthur/ProtectedAreas_Data/merge_africa_grid_wdpa_pre2008.csv")
#Create new column with percentage of cell covered by protected area
pa_2000$wdpapct_2000 <- NA
pa_2000$wdpapct_2000 <- pa_2000$wdpa_pre2001_africa.na.sum/pa_2000$wdpa_pre2001_africa.na.count

pa_2007$wdpapct_2007 <- NA
pa_2007$wdpapct_2007 <- pa_2007$wdpa_pre2008_africa.na.sum/pa_2007$wdpa_pre2008_africa.na.count 
#Merge percentage of cell covered by protected area into AFR cell dataset for pre2001 and pre2008, drop out the sum and count columns
AFRcells<-merge(pa_2000,pa_2007,by.x="ID",by.y="ID")
AFRcells <- AFRcells[,-grep("(africa)", names(AFRcells))]

#merge AFRcells into Panel_Data
Panel_Data_add <- Panel_Data
Panel_Data_add<-merge(Panel_Data_add,AFRcells, by.x="ID",by.y="ID")

## Add in GPW4 Pop Density Data, updated data for 2000,2005,2010,2015
pop <- read.csv("/home/aiddata/Desktop/Github/MacArthur/GPW4_Extracts/merge_africa_grid.csv")
Panel_Data_add<-merge(Panel_Data_add,pop,by.x="ID",by.y="ID")
#Apply 2000 values to years 2001-2004 of Panel_Data, 2005 values to years 2005-2009, 2010 values to years 2010-2014
Panel_Data_add$Pop<-NA
Panel_Data_add$Pop[Panel_Data_add$Year<=2014]<-Panel_Data_add$gpw_v4_density.2010.mean[Panel_Data_add$Year<=2014]
Panel_Data_add$Pop[Panel_Data_add$Year<=2009]<-Panel_Data_add$gpw_v4_density.2005.mean[Panel_Data_add$Year<=2009]
Panel_Data_add$Pop[Panel_Data_add$Year<=2004]<-Panel_Data_add$gpw_v4_density.2000.mean[Panel_Data_add$Year<=2004]
#Maintain Pop_2000 for baseline value
names(Panel_Data_add)[names(Panel_Data_add)=="gpw_v4_density.2000.mean"]="Pop_2000"

##create ntl panel dataset that matches 2001-2014 years of Panel_Data to merge into main dataset
source("SciClone_functions.R")
ntl2013<-read.csv("/home/aiddata/Desktop/Github/MacArthur/ntl_extracts/merge_africa_grid_2013.csv")
# get pre-trends for creation of imputed NTL_2014 values (ntl still in ncc4_1992e name format)
#requires shape file with ntl 1992 to 2012 data, then need to merge in 2013 ntl data
AOI_cells_ntl<-AOI_cells[c(1:102)]
AOI_cells_ntl<-merge(AOI_cells_ntl,ntl2013,by.x="ID",by.y="ID")
names(AOI_cells_ntl)[names(AOI_cells_ntl)=="v4composites_calibrated.2013.mean"]="ncc4_2013e"
#create five year trend to impute 2014 value
AOI_cells_ntl$ntltrend_0913<-timeRangeTrend(AOI_cells_ntl,"ncc4_[0-9][0-9][0-9][0-9]e",2009,2013,"ID","y")
AOI_cells_ntl@data$ncc4_2014e<-NA
AOI_cells_ntl@data$ncc4_2014e=AOI_cells_ntl@data$ncc4_2013e+AOI_cells_ntl@data$ntltrend_0913
AOI_cells_ntl@data$neg2014[AOI_cells_ntl@data$ncc4_2014e<0]<-1
AOI_cells_ntl@data$ncc4_2014e[AOI_cells_ntl@data$neg2014==1]<-0
#create ntl pre-trend for 1992-2003
AOI_cells_ntl$ntl_pretrend<-timeRangeTrend(AOI_cells_ntl,"ncc4_[0-9][0-9][0-9][0-9]e",1992,2007,"ID","y")
#create non-shape file and rename to something obvious
ntl<-AOI_cells_ntl@data
for (i in 2:length(ntl)) {
  
  if (substr(colnames(ntl)[i], 1, 4) == "ncc4"){
    
    name = "NTL"
    year = substr(colnames(ntl)[i], 6, 9)
    dt = paste(name,"_",year,sep="")
    colnames(ntl)[i] <- dt
  }
}

#convert it from wide to long form and reshape into panel
ntl_late<-ntl[c(1,82:105,107)]
ntl_long <- grep("NTL", names(ntl_late))
ntl_reshape <- c(ntl_long)
ntl_panel <- reshape(ntl_late, varying=ntl_reshape,direction="long", idvar="ID", sep="_", timevar="Year")
ntl_panel$NTL_2007<-NA
ntl_panel$NTL_2007<-ntl_panel$NTL[ntl_panel$Year==2007]
#merge into Panel_Data_add
Panel_Data_add1<-merge(Panel_Data_add,ntl_panel,by=c("ID","Year"))

Panel_Data_add<-Panel_Data_add1



#write.csv(Panel_Data, "/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra.csv")
#Panel_Data<-read.csv("/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra.csv")
write.csv(Panel_Data_add,"/home/aiddata/Desktop/Github/MacArthur/modelData/tanzania_infra_panel_data_add_Thresh15_JUNE.csv")

