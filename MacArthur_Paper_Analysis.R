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

forest_thresh = 10
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
csv <- paste(active_dir_path, "/extracts/sea.csv", sep="")
json <- paste(active_dir_path, "/extracts/sea.json", sep="")
dta <- read.csv(csv)
vars <- fromJSON(txt=json)

#--------------------------------------------------#
#Subset the Cell Dataframe 
#--------------------------------------------------#
dta2 <- dta[dta$NAME_0 == "Cambodia",]

#--------------------------------------------------#
#Load and Subset the Cell Spatial Datframe
#--------------------------------------------------#
spdf_cells <- paste(active_dir_path, "/grids/sea_grid.shp", sep="")
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
spdf_adm <- spdf_adm[spdf_adm@data$NAME_0 == "Cambodia",]

#--------------------------------------------------#
#Load and Subset the MacArthur Aid Data
#--------------------------------------------------#
location_csv <- paste(active_dir_path, "/MacArthur_Geocoded_data/locations.csv", sep="")
locations <- read.csv(location_csv)
locations2 <- locations[grep("Cambodia", locations$gazetteer_adm_name),]
coords = cbind(locations2$longitude, locations2$latitude)

#--------------------------------------------------#
#Create a Spatial Dataframe of the MacArthur Aid Data
#--------------------------------------------------#
Mac_spdf <- SpatialPointsDataFrame(coords, locations2)
proj4string(Mac_spdf) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
#Join in the year of aid allocation
Mac_years <- paste(active_dir_path, "/MacArthur_Geocoded_data/projects.csv", sep="")
Mac_yr_dta <- read.csv(Mac_years)
Mac_spdf <- merge(Mac_spdf, Mac_yr_dta, by="project_id")

#--------------------------------------------------#
#Create threshdolded forest datasets
#--------------------------------------------------#
ndviDTA <- dta2[c("tc00_e", "lnyx_2000e")]
ndviDTA_for <- ndviDTA[ndviDTA$tc00_e >= forest_thresh,]
ndviDTA_notfor <- ndviDTA[ndviDTA$tc00_e <= forest_thresh,]


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
correlogram_data <- correlog(x = coordinates(AOI_cells)[,1], y = coordinates(AOI_cells)[,2], z=AOI_cells$lnyx_1999e, increment=5, latlon=TRUE, na.rm=TRUE, resamp=50)

#save data into a function to calculate the distance-decay penalty later.
#Chinese projects are "weighted" according to their distance.
#The absolulute correlation for a given distance is used as a weight.
#Projects at distances with a higher positive or inverse correlation are given the highest
#weights.

dVals <<- abs(correlogram_data$mean.of.class)
cVals <<- abs(correlogram_data$correlation)

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
decay_dMatrix_adj <- apply(decay_dMatrix, 1:2, function(x){(cVals[which.min(abs(dVals - x))] * x)[[1]]})
AOI_cells$thresh_weightedDist <- colSums(decay_dMatrix_adj) / 1000


#--------------------------------------------------#
#Calculate the over-time treatment effects
#--------------------------------------------------#
all_years <- unique(Mac_spdf$transactions_start_year)
all_years <- all_years[!is.na(all_years)]

dYears <- list()
#Drop all MacArthur projects that have no start or end date.
Mac_spdf <- Mac_spdf[!is.na(Mac_spdf$transactions_start_year),]

for(years in 1:length(all_years))
{
  year <- all_years[years]
  ThisYearMac <- Mac_spdf[Mac_spdf@data$transactions_start_year == year,]
  dYears[[years]] <- RDist(AOI_cells, ThisYearMac)
}

#In dMatrix, every column is a cell (referenced in order to AOI_cells)
#Every row is a MacArthur project for that year (referenced in order to Mac_spdf)

AvgYears <- vector()
for(years in 1:length(all_years))
{
  AvgYears[[years]] <- mean(dYears[[years]]) / 1000
}

Avg_MinYears <- vector()
for(years in 1:length(all_years))
{
  col_mins_year <- do.call(pmin, lapply(1:nrow(dYears[[years]]), function(i)dYears[[years]][i,]))
  nameRef <- paste("MinYr_",all_years[years], sep="")
  AOI_cells@data[nameRef] <- col_mins_year / 1000
  Avg_MinYears[[years]] <- mean(col_mins_year) / 1000
}

Dist_Decay_Yrs <- vector()
for(years in 1:length(all_years))
{
  t_dyears <- dYears[[years]] / 1000
  decay_dMatrix_adj <- apply(t_dyears, 1:2, function(x){(cVals[which.min(abs(dVals - x))] * x)[[1]]})
  nameRef <- paste("DecayYr_",all_years[years], sep="")
  AOI_cells@data[nameRef] <- colMeans(decay_dMatrix_adj) / 1000
  Dist_Decay_Yrs[[years]] <- mean(colMeans(decay_dMatrix_adj) / 1000)
}

CountProj_Years <- vector()
for(years in 1:length(all_years))
{
  CountProj_Years[[years]] <- nrow(dYears[[years]])
}

#Build a quick temporal dataframe for plotting and ordering
TempDF <- cbind.data.frame(all_years, AvgYears, Avg_MinYears, CountProj_Years, Dist_Decay_Yrs)
TempDF <- TempDF[with(TempDF, order(TempDF[,1])),]




#--------------------------------------------------#
#Pre-processing for analysis
#--------------------------------------------------#
DFa <- AOI_cells@data
#Drop irrelevant variables:
dropvars <- c("XMIN","XMAX","YMIN","YMAX","OBJECTID","ID_0","ISO","NAME_0","HASC_2","ID_1","NAME_1","NAME_2","CCN_2","CCA_2","TYPE_2","ENGTYPE_2","NL_NAME_2","VARNAME_2","Shape_Leng","Shape_Area", "thresh_tot_proj","thresh_totDist","thresh_avgDist","thresh_weightedDist")

DFa <- DFa[,!(names(DFa) %in% dropvars)]

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

names(DFa)[25:143]<- sapply(names(DFa)[25:143], function(x) {rename_header(x, substrRight(x,1))})

names(DFa)[6:19] <- sapply(names(DFa)[6:19], function(x){substr(x, 5, nchar(x))})

#drop data for the year 2000, 2013-2015 (NTL ends in 2012)
DFa <- DFa[, -grep("(2000)", names(DFa))]
DFa <- DFa[, -grep("(2015)", names(DFa))]
DFa <- DFa[, -grep("(2013)", names(DFa))]
DFa <- DFa[, -grep("(2014)", names(DFa))]

PCloss <- grep("^loss", names(DFa))
mean_ln <- grep("^elnyx", names(DFa))
NTL <- grep("^encc4", names(DFa))
minairTemp <- grep("^mat41", names(DFa))
maxairTemp <- grep("^xat41", names(DFa))
meanairTemp <- grep("^eat41", names(DFa))
minPre <- grep("^mpc41", names(DFa))
maxPre <- grep("^xpc41", names(DFa))
meanPre <- grep("^epc41", names(DFa))
MinDist <- grep("^MinYr", names(DFa))
DecayDist <- grep("^DecayYr", names(DFa))


#--------------------------------------------------#
#Selection of temporally-varying variables and shift from wide- to long-form
#--------------------------------------------------#
all_reshape <- c(PCloss, mean_ln, NTL, minairTemp, maxairTemp, meanairTemp, minPre, maxPre, meanPre, MinDist, DecayDist)
DFa <- reshape(DFa, varying=all_reshape,direction="long", idvar="ID", sep="_", timevar="Year")

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


#--------------------------------------------------#
#Modeling
#--------------------------------------------------#
initModel <- lm(Forest_Loss ~ DecayYr + NighttimeLights + MinTemp + MaxTemp + MeanTemp + MaxPrecip + MeanPrecip + MinPrecip + factor(District), data=DFa)
cluster <- cluster.vcov(initModel, cbind(DFa$Year, DFa$ID_2), force_posdef=TRUE)
CMREG <- coeftest(initModel, cluster)