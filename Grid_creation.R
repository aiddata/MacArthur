##############
# MacArthur Grid #
##############library(rgdal)
library(maptools)
library(gstat)
library(sp)

#import the border shp file
border<-readShapePoly("/home/aid_data/Desktop/GitRepo/MacArthur/Regions_Ind/SEA.shp")
proj4string(border) <- CRS("+proj=longlat +datum=WGS84")
borderProj <- spTransform(border, CRS("+init=epsg:28992")) 
vals<-borderProj@bbox
deltaLong <- as.integer((vals[1,2] - vals[1,1]) + 1.5)
deltaLat <- as.integer((vals[2,2] - vals[2,1]) + 1.5)
gridRes <-5000   #change this value to change the grid size (in meters)
gridSizeX <- deltaLong / gridRes
gridSizeY <- deltaLat / gridRes
grd <- GridTopology(vals[,1],c(gridRes,gridRes),c(gridSizeX,gridSizeY))
pts <- SpatialPoints(coordinates(grd))
pts1 <- SpatialPointsDataFrame(as.data.frame(pts), data=as.data.frame(rep(1,nrow(as.data.frame(pts)))))
Overlay=overlay(pts1,border)
pts1$border=Overlay
nona<-na.exclude(as.data.frame(pts1))
coordinates(nona)=~x+y
gridded(nona) <- TRUE
proj4string(nona)=CRS("+init=epsg:28992")  

SPGRD <- as.SpatialPolygons.GridTopology(grd)

IDs <- sapply(slot(SPGRD, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(rep(0, length(IDs)), row.names=IDs)

SPDF_x <- SpatialPolygonsDataFrame(SPGRD, df)

writePolyShape(SPDF_x, "/home/aid_data/Desktop/GitRepo/MacArthur/Regions_Grids/SEA_grid.shp")
writePolyShape(borderProj, "/home/aid_data/Desktop/GitRepo/MacArthur/Regions_Grids/SEA_Proj.shp")

