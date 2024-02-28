library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#dem
# solar1 <- rast('gis/dem/tifs/mi_dem_30_solar.tif')
# solar <- mean(solar1) 
# names(solar) <- 'solar'
# writeRaster(solar, 'gis/dem/tifs/solar.tif')
solar <- rast('gis/dem/tifs/solar.tif')
openess <- rast('gis/dem/tifs/mi_dem_30_layers.tif')
names(openess) <- c('shade','popen','nopen')
slope <- rast('gis/dem/tifs/mi_dem_30_slope.tif')
names(slope) <- c('slope','aspect')
dem <- rast('gis/dem/tifs/mi_dem_30.tif')
names(dem) <- c('elev')
dem <- c(dem, slope, openess, solar)

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
for (i in 1:12){
  prism <- c(prism, get(paste0('th',month[i])))
}
for (i in 1:12){
  prism <- c(prism, get(paste0('tl',month[i])))
}
for (i in 1:12){
  prism <- c(prism, get(paste0('p',month[i])))
}

