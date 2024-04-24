library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
mnfi <- read_sf('C:/a/Ecological_Sites/GIS/Vegetation/Veg1800/MichLP_1800veg.shp')
dirlst <- list.files('gis/models90')
dirlstnms <- stringr::str_split_fixed(dirlst, pattern='\\.',  n=2)[,1]
for(i in 1:length(dirlst)){#i=1
x <- rast(paste0('gis/models90/',dirlst[i]))
assign(dirlstnms[i],x)
}
all = Abies+Acer+Betula+Carya+Fagus+Fraxinus+Larix+Picea+Pinus+PopLirio+Quercus+Thuja+Tilia+Tsuga+Ulmus
for(i in 1:length(dirlst)){#i=1
   assign(dirlstnms[i], rast(mget(dirlstnms[i]))/all)
}
basalarea <-  rast('gis/models90/basalarea.tif')

spp <- rast(mget(dirlstnms))
mnfi <- st_transform(mnfi, crs(spp))
mnfi.vect <- vect(mnfi)
spam <- spatSample(mnfi.vect, 10000)
spam.mnfi <- extract(mnfi.vect, spam)
spam.spp <- extract(spp, spam)

spam.all <- spam.mnfi |> left_join(spam.spp, join_by(id.y==ID))

library(rpart)
library(rpart.plot)

rp <- rpart(COVERTYPE ~ basalarea+Abies+Acer+Betula+Carya+Fagus+Fraxinus+Larix+Picea+Pinus+PopLirio+Quercus+Thuja+Tilia+Tsuga+Ulmus, data = spam.all,
            method="class", control = list(maxdepth = 4, cp=0.002, minsplit=100))
png(filename="rpart.png",width = 10, height = 3, units = 'in', res = 600)
rpart.plot(rp, extra=108,legend.cex=0.5, digits=2)
dev.off()