library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
mnfi <- read_sf('C:/a/Ecological_Sites/GIS/Vegetation/Veg1800/MichLP_1800veg.shp')
dirlst <- list.files('gis/models90')
dirlstnms <- stringr::str_split_fixed(dirlst, pattern='\\.',  n=2)[,1] |> unique()
for(i in 1:length(dirlstnms)){#i=1
x <- rast(paste0('gis/models90/',dirlstnms[i],'.tif'))
assign(dirlstnms[i],x)
}
all = Abies+Acer+Betula+Carya+Fagus+Fraxinus+Larix+Picea+Pinus+PopLirio+Quercus+Thuja+Tilia+Tsuga+Ulmus
for(i in 1:length(dirlstnms)){#i=1
  assign(dirlstnms[i], rast(mget(dirlstnms[i]))/(all+0.001))
}
for(i in 1:length(dirlstnms)){#i=1
  assign(dirlstnms[i], rast(mget(dirlstnms[i]))/(minmax(rast(mget(dirlstnms[i])))[2]))
}

basalarea <-  rast('gis/models90/basalarea.tif')

spp <- rast(mget(dirlstnms))
mnfi <- st_transform(mnfi, crs(spp))
mnfi.vect <- vect(mnfi)
spam <- spatSample(mnfi.vect, 10000)
spam.mnfi <- extract(mnfi.vect, spam)
spam.spp <- extract(spp, spam)

spam.all <- spam.mnfi |> left_join(spam.spp, join_by(id.y==ID))
spam.all <- spam.all |> subset(!is.na(COVERTYPE))
spam.all.summary <- spam.all |> group_by(COVERTYPE, VEGCODE) |> 
  summarise(across(.cols=c(basalarea,Abies,Acer,Betula,Carya,Fagus,Fraxinus,Larix,Picea,Pinus,PopLirio,Quercus,Thuja,Tilia,Tsuga,Ulmus), ~ round(mean(.x, na.rm = TRUE),3)))
spam.all.summary <- spam.all |> group_by(COVERTYPE) |> 
  summarise(across(.cols=c(basalarea,Abies,Acer,Betula,Carya,Fagus,Fraxinus,Larix,Picea,Pinus,PopLirio,Quercus,Thuja,Tilia,Tsuga,Ulmus), ~ round(mean(.x, na.rm = TRUE),3)))
spam.all.summary.t <- spam.all.summary |> subset(select = -COVERTYPE) |> t()
colnames(spam.all.summary.t) <- spam.all.summary$COVERTYPE
spam.all <- spam.all |> group_by(COVERTYPE) |> mutate(wts = 1/(sum(length(COVERTYPE))+5)) |> ungroup()
library(rpart)
library(rpart.plot)

rp <- rpart(COVERTYPE ~ basalarea+Abies+Acer+Betula+Carya+Fagus+Fraxinus+Larix+Picea+Pinus+PopLirio+Quercus+Thuja+Tilia+Tsuga+Ulmus, data = spam.all,
            method="class", control = list(maxdepth = 5, cp=0.002, minsplit=100), weights = spam.all$wts)
png(filename="rpart.png",width = 10, height = 3, units = 'in', res = 600)
rpart.plot(rp, extra=108,legend.cex=0.5, digits=2)
dev.off()









