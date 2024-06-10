library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

vars90 <- rast('gis/vars90.tif')

sandy <- ifel(!(vars90$sand150 >= 80| vars90$sand150 > 70), 3, ifel(vars90$pH50 >= 6 | vars90$carbdepth < 100,2, ifel(vars90$spodic >= 0.8, 1,0)))

plot(sandy)
catena <- ifel(is.na(vars90$hydric),NA, ifel(vars90$hydric > 0.5, 1, ifel(vars90$watertable < 120, 2, ifel(vars90$slope >= tan(15/100), 4,ifel(vars90$rockdepth < 100, 5,3)))))
plot(catena)

soilcatena <- catena * 10 + sandy
writeRaster(soilcatena, 'gis/soilcatena.tif', overwrite=T)