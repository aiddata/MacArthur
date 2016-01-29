##############
# MacArthur Grid #
##############

library(rgdal)
library(maptools)
library(gstat)
library(sp)

#import the border shp file
border<-readShapePoly("/home/dan/Desktop/GitRepo/MacArthur/Regions_Ind/SEA.shp")
proj4string(border) <- CRS("+proj=longlat +datum=WGS84")
borderProj <- spTransform(border, CRS("+init=epsg:28992"))



vals <-borderProj@bbox
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
writeAsciiGrid(nona,"/home/dan/Desktop/GitRepo/MacArthur/Regions_Ind/SEA_grid.asc")

writePolyShape(borderProj,"/home/dan/Desktop/GitRepo/MacArthur/Region_Grids/SEA_reProj.shp" )
SPGRD <- as.SpatialPolygons.GridTopology(grd)
writePolyShape(as.SpatialPolygons.GridTopology(grd), "/home/dan/Desktop/GitRepo/MacArthur/Region_Grids/SEA_GRID.shp")












gridWGS <- spTransform(border, CRS("+proj=longlat +datum=WGS84"))






#For the Co-kriging we need to obtain the value of the covariate for each observation
#for doing that we can use the function overlay
over=overlay(zinc,data)
data$zinc=over$zinc.asc
str(as.data.frame(data))

#also the prediction grid need to be overlayed with the covariate
over=overlay(zinc,nona)
nona$zinc=over$zinc.asc


#for the cokriging, the first thing to do is create an object with the
#function gstat() that contains both the variable and the covariate
g<-gstat(id="lead",formula=lead~1,data=data)
g<-gstat(g,id="zinc",formula=zinc~1,data=data)


#Fitting the variogram
#first, plot the residual variogram
vario<-variogram(g)
plot(vario)

#now we can fit the linear model of coregionalization
g<-gstat(g,id=c("lead","zinc"),model=vgm(psill=cov(data$lead,data$zinc),model="Sph",range=sqrt(areaSpatialGrid(zinc))/4,nugget=0))
g<-fit.lmc(vario,g,model=vgm(psill=cov(data$lead,data$zinc),model="Sph",range=sqrt(areaSpatialGrid(zinc))/4,nugget=0))

plot(vario,g$model)

k<-predict.gstat(g,nona)

#Validation
#create training and test
#repredict the LMC
#perfomr the prediction

i<-sample(nrow(data),round(nrow(data)*10/100)) #exclude 10% of the data
training<-data[!data$ID%in%i,]
test<-data[data$ID%in%i,]

gv<-gstat(id="lead",formula=lead~1,data=training)
gv<-gstat(gv,id="zinc",formula=zinc~1,data=training)

gv<-gstat(gv,id=c("lead","zinc"),model=vgm(psill=cov(training$lead,training$zinc),model="Sph",range=sqrt(areaSpatialGrid(zinc))/4,nugget=0))
gv<-fit.lmc(variogram(gv),gv,model=vgm(psill=cov(training$lead,training$zinc),model="Sph",range=sqrt(areaSpatialGrid(zinc))/4,nugget=0))

plot(variogram(gv),gv$model)

krige_cross<-predict.gstat(gv,test)
str(krige_cross)

#Goodness of fit indexes
RSQR<-as.numeric(cor.test(test$lead,krige_cross$lead.pred)$estimate)^2      					#Pearson's R Squared
RMSD<-sqrt(sum((test$lead-krige_cross$lead.pred)^2)/length(test$lead))      #Root Mean Square Deviation

#Print Variograms and Map
jpeg("Variogram.jpg",1200,1200,res=300)
plot(vario,g$model)
dev.off()

jpeg("Prediction_Map.jpg",1200,1200,res=300)
spplot(k,"lead.pred",col.regions=terrain.colors(50),main="Prediction Map",scales=list(draw=T))
dev.off()

jpeg("Error_Map.jpg",1200,1200,res=300)
spplot(k,"lead.var",col.regions=heat.colors(50),main="Error Map",scales=list(draw=T))
dev.off()


#References:
#- Applied Spatial Data Analysis with R. Bivand,Pebesma,G?mez-Rubio (2008)
#- cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf
#- www.itc.nl/~rossiter/teach/R/R_ck.pdf
#- www.ic.arizona.edu/ic/math574/class.../cokriging%20in%20gstat.pdf
