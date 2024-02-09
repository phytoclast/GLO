#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

solar <- rast('D:/GIS/DEM/hmnfsolarmean.tif'); names(solar) = 'solar'
toi <- rast('D:/GIS/DEM/hmnf-openess.tif'); names(toi) = 'toi'
toip <- rast('D:/GIS/DEM/hmnf-pos-openess.tif'); names(toip) = 'toip'
toin <- rast('D:/GIS/DEM/hmnf-neg-openess.tif'); names(toin) = 'toin'
toip100 <- rast('D:/GIS/DEM/hmnf-pos-openess100.tif'); names(toip100) = 'toip100'
toin100 <- rast('D:/GIS/DEM/hmnf-neg-openess100.tif'); names(toin100) = 'toin100'
twi <- rast('D:/GIS/DEM/hmnftwi.tif'); names(twi) = 'twi'
tpi <- rast('D:/GIS/DEM/hmnftopographicpositionindex.tif'); names(tpi) = 'tpi'
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'
slope <- rast('D:/GIS/DEM/hmnfslope.tif'); names(slope) = 'slope'#radians
aspect <- rast('D:/GIS/DEM/hmnfaspect.tif'); names(aspect) = 'aspect'#radians
bt <- rast('D:/scripts/snow/output/bt.90.alt.tif'); names(bt) = 'bt'
tgs <- rast('D:/scripts/snow/output/tgs.90.alt.tif'); names(tgs) = 'tgs'
ppt <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/P/p0112/w001001.adf'); names(ppt) = 'ppt'
p1 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p1/w001001.adf'); names(p1) = 'p1'
p2 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p2/w001001.adf'); names(p2) = 'p2'
p3 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p3/w001001.adf'); names(p3) = 'p3'
p4 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p4/w001001.adf'); names(p4) = 'p4'
slope.500 <- rast('D:/GIS/DEM/slope.500.tif'); names(slope.500) = 'slope.500'

bt <- project(bt, dem)
tgs <- project(tgs, dem)
ppt <- project(ppt, dem)
p1 <- project(p1, dem, method='bilinear')
p2 <- project(p2, dem, method='bilinear')
p3 <- project(p3, dem, method='bilinear')
p4 <- project(p4, dem, method='bilinear')

brk <- c(solar, toi, toip, toin, toip100, toin100, twi, tpi, dem, slope, aspect, slope.500, bt, tgs, ppt,p1,p2,p3,p4)
soil <- rast('D:/scripts/mattnasis/milev30.tif')
soil <- project(soil, dem, method='bilinear')
plot(tgs)