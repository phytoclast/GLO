library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#dem
# solar1 <- rast('gis/dem/tifs/mi_dem_30_solar.tif')#SAGA solar radiation 14 steps 
# solar <- mean(solar1) 
# names(solar) <- 'solar'
# writeRaster(solar, 'gis/dem/tifs/solar.tif')
solar <- rast('gis/dem/tifs/solar.tif') 
openess <- rast('gis/dem/tifs/mi_dem_30_layers.tif') #SAGA topographic openess to 10000 m
names(openess) <- c('shade','popen','nopen')
slope <- rast('gis/dem/tifs/mi_dem_30_slope.tif') #slope SAGA radians
names(slope) <- c('slope','aspect')
dem <- rast('gis/dem/tifs/mi_dem_30.tif')#DEM from USGS Feb 2024
slope500 <- rast('gis/dem/tifs/slope500.tif');names(slope500) <- 'slope500' #mean slope 500 m radius
names(dem) <- c('elev')
dem <- c(dem, slope, slope500, openess, solar)

#soil
misoil <- rast(paste0('gis/soilcov/','MI','soil.tif'))
wisoil <- rast(paste0('gis/soilcov/','WI','soil.tif'))
ohsoil <- rast(paste0('gis/soilcov/','OH','soil.tif'))
insoil <- rast(paste0('gis/soilcov/','IN','soil.tif'))
ilsoil <- rast(paste0('gis/soilcov/','IL','soil.tif'))

if(file.exists('gis/soilcov/misoilcrop.tif')){
misoil <- rast('gis/soilcov/misoilcrop.tif')}else{misoil <- project(misoil, dem, method = 'near', filename='gis/soilcov/misoilcrop.tif')}
if(file.exists('gis/soilcov/wisoilcrop.tif')){
wisoil <- rast('gis/soilcov/wisoilcrop.tif')}else{wisoil <- project(wisoil, dem, method = 'near', filename='gis/soilcov/wisoilcrop.tif')}
if(file.exists('gis/soilcov/ohsoilcrop.tif')){
ohsoil <- rast('gis/soilcov/ohsoilcrop.tif')}else{ohsoil <- project(ohsoil, dem, method = 'near', filename='gis/soilcov/ohsoilcrop.tif')}
if(file.exists('gis/soilcov/insoilcrop.tif')){
insoil <- rast('gis/soilcov/insoilcrop.tif')}else{insoil <- project(insoil, dem, method = 'near', filename='gis/soilcov/insoilcrop.tif')}
if(file.exists('gis/soilcov/ilsoilcrop.tif')){
ilsoil <- rast('gis/soilcov/ilsoilcrop.tif')}else{ilsoil <- project(ilsoil, dem, method = 'near', filename='gis/soilcov/ilsoilcrop.tif')}


if(file.exists('gis/soilcov/allsoilcrop.tif')){
  allsoil <- rast('gis/soilcov/allsoilcrop.tif')}else{allsoil <- merge(misoil,insoil,ohsoil,ilsoil,wisoil, na.rm=TRUE, filename='gis/soilcov/allsoilcrop.tif')}

plot(allsoil$spodic)

#climate 
# remotes::install_github("phytoclast/climatools", dependencies = FALSE)
library(climatools)
prismTpath <- 'C:/a/Ecological_Sites/GIS/Climate/PRISM2010/CorrectedT/'
prismPpath <- 'C:/a/Ecological_Sites/GIS/Climate/PRISM2010/P/'

month <- c('01','02','03','04','05','06','07','08','09','10','11','12')
#load basic files ----
for (i in 1:12){#i=1
  x = rast(paste0(prismPpath,'p',month[i],'/w001001.adf'))
  names(x) <- paste0('p',month[i])
  assign(paste0('p',month[i]), x)
}
for (i in 1:12){
  x = rast(paste0(prismTpath,'tr',month[i],'/w001001.adf'))
  names(x) <- paste0('tr',month[i])
  assign(paste0('tr',month[i]),x)
}
for (i in 1:12){
  x = rast(paste0(prismTpath,'t',month[i],'/w001001.adf'))
  names(x) <- paste0('tr',month[i])
  assign(paste0('t',month[i]),x)
}
for (i in 1:12){
  xh <- get(paste0('t',month[i])) + get(paste0('tr',month[i]))/2
  xl <- get(paste0('t',month[i])) - get(paste0('tr',month[i]))/2
  names(xh) <- paste0('th',month[i])
  names(xl) <- paste0('tl',month[i])
  assign(paste0('th',month[i]), xh)
  assign(paste0('tl',month[i]), xl)
}
for (i in 1:12){assign(paste0('tr',month[i]), NULL)};rm(x,xh,xl)
latlon <- as.data.frame(t01, xy=T); colnames(latlon) <- c('x','y','z')
latlon[,3] <- latlon[,2]; 
lat = rast(x=latlon, type="xyz", crs=crs(t01)); names(lat) <- 'lat'
latlon[,3] <- latlon[,1]; 
lon = rast(x=latlon, type="xyz", crs=crs(t01));  names(lon) <- 'lon'
usadem <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/PRISM_us_dem_800m_bil/PRISM_us_dem_800m_bil.bil');names(usadem) <- 'elev'
lat <- project(lat, t01);lon <- project(lon, t01)
prism <- c(lat,lon, usadem)
prism <- c(prism, rast(mget(paste0('th',month))))
prism <- c(prism, rast(mget(paste0('tl',month))))
prism <- c(prism, rast(mget(paste0('p',month))))
writeRaster(prism, 'gis/climate/prism.tif', overwrite=T)


###load ----
library(climatools)
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
prism <- rast('gis/climate/prism.tif')
timeA <-  Sys.time()
e01 <- GetPET.block(1,prism)
e02 <- GetPET.block(2,prism)
e03 <- GetPET.block(3,prism)

e04 <- GetPET.block(4,prism)
e05 <- GetPET.block(5,prism)
e06 <- GetPET.block(6,prism)

e07 <- GetPET.block(7,prism)
e08 <- GetPET.block(8,prism)
e09 <- GetPET.block(9,prism)

e10 <- GetPET.block(10,prism)
e11 <- GetPET.block(11,prism)
e12 <- GetPET.block(12,prism)
Sys.time() - timeA
timeA <-  Sys.time()
pet <- sum(e01,e02,e03,e04,e05,e06,e07,e08,e09,e10,e11,e12); names(pet) <- 'e'
Sys.time() - timeA
timeA <-  Sys.time()
writeRaster(pet, 'gis/climate/pet.tif', overwrite=T)
Sys.time() - timeA
#----
pet <- rast('gis/climate/pet.tif')
names(pet) <- 'e'

p <- ApplyClim.rast(prism, name='p')
p1 <- ApplyClim.rast(prism, jan='p01', mons=c(12,1,2), fun='sum', name='pq1')
p2 <- ApplyClim.rast(prism, jan='p01', mons=c(3,4,5), fun='sum', name='pq2')
p3 <- ApplyClim.rast(prism, jan='p01', mons=c(6,7,8), fun='sum', name='pq3')
p4 <- ApplyClim.rast(prism, jan='p01', mons=c(9,10,11), fun='sum', name='pq4')

tblock <- climatools::meanT.rast(prism)

Twh <- ApplyClim.rast(prism, jan='th01', fun='max', name='Twh')
Tw <- ApplyClim.rast(tblock, jan='t01', fun='max', name='Tw')
Tcl <- ApplyClim.rast(prism, jan='tl01', fun='min', name='Tcl')
Tc <- ApplyClim.rast(tblock, jan='t01', fun='min', name='Tc')

btblock <- ifel(tblock > 0, tblock, 0)

Tg1 <- ApplyClim.rast(btblock, jan='t01', mons=c(12,1,2,3,4), fun='mean')
Tg2 <- ApplyClim.rast(btblock, jan='t01', mons=c(5,6,7,8,9,10), fun='mean')
Tg <- max(Tg1, Tg2);names(Tg) <- 'Tg'
m <- p/(p+pet+0.0001); names(m) <- 'm'
prism2 <- c(p,p1,p2,p3,p4,Twh,Tw,Tc,Tcl,Tg,pet,m)
writeRaster(prism2, 'gis/climate/prism2.tif', overwrite=T)

#project to project area ---- 
library(climatools)
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
prism <- rast('gis/climate/prism.tif')
prism2 <- rast('gis/climate/prism2.tif')
solar <- rast('gis/dem/tifs/solar.tif')
openess <- rast('gis/dem/tifs/mi_dem_30_layers.tif')
names(openess) <- c('shade','popen','nopen')
slope <- rast('gis/dem/tifs/mi_dem_30_slope.tif')
names(slope) <- c('slope','aspect')
dem <- rast('gis/dem/tifs/mi_dem_30.tif')
names(dem) <- c('elev')
dem <- c(dem, slope, openess, solar)

prism.repro <- project(prism2, dem, method='bilinear')
writeRaster(prism.repro, 'gis/climate/prism.repro.tif', overwrite=T)
prism.repro <- rast('gis/climate/prism.repro.tif')
Tg30 <- enhanceRast(prism2$Tg, prism$elev, dem$elev); names(Tg30) <- 'Tg30'
writeRaster(Tg30, 'gis/climate/Tg30.tif', overwrite=T)
# Tg30 <- rast('gis/climate/Tg30.tif')
# Tg500 <- climatools::focalmed(Tg30, 500)
# plot(Tg500)
# writeRaster(Tg500, 'gis/climate/Tg500.tif', overwrite=T)

