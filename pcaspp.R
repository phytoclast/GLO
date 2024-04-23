library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dirlst <- list.files('gis/models90')
dirlstnms <- stringr::str_split_fixed(dirlst, pattern='\\.',  n=2)[,1]
for(i in 1:length(dirlst)){#i=1
x <- rast(paste0('gis/models90/',dirlst[i]))
assign(dirlstnms[i],x)
}

spp <- rast(mget(dirlstnms))

pca <- terra::prcomp(spp, maxcell=ncell(spp)*0.1)
pca.rast <- predict(spp,pca)

plot(pca.rast$PC2)
writeRaster(pca.rast, 'gis/pcaspp.tif', overwrite=T)
pca.norm <- terra::prcomp(spp, maxcell=ncell(spp)*0.1, scale.=TRUE)
pca.rast.norm <- predict(spp,pca.norm)

plot(pca.rast.norm$PC2)
writeRaster(pca.rast.norm, 'gis/pcasppnorm.tif', overwrite=T)


pca$rotation