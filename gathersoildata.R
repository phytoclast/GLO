library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
state='MI'
drlist <- list.files('gis/soil')
drlist <- drlist[grepl('.tif$', drlist)]
for(i in 1:length(drlist)){#i=1
  grid = rast(paste0('gis/soil/', drlist[i]))
  if(i==1){
    soil <- grid
  }else{
    soil <- c(soil, grid)
  }
}
plot(soil$floodfrq)
writeRaster(soil, paste0('gis/soilcov/',state,'soil.tif'), overwrite=T)


