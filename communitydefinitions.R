library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- readRDS('points/pts.geo1.RDS')


pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
vars270 <- rast('gis/vars270.tif')
# pca90 <- princomp(vars270, cor = T, maxcell=ncell(vars270)*0.1)
# pca90grid <- predict(pca90, pca90)
# writeRaster(pca90grid,'gis/pca90grid.tif', overwrite=T)
pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2
cors <- readRDS('output/cors.RDS')

for(i in 1:length(ttt)){
x <- rast(paste0('gis/models270/',ttt[i],'.tif'))
assign(ttt[i], x)         };rm(x)

fire <- ifel(Acer+Fagus+Tsuga+Tilia+Abies+Picea+Thuja >= Pinus+Quercus, 0,1)
fire <- (Pinus+Quercus)/(Acer+Fagus+Tsuga+Tilia+Abies+Picea+Thuja+Pinus+Quercus+0.001)
plot(fire)
wetness <- Abies*cors[cors$column %in% 'Abies','hydric']+
  Acer*cors[cors$column %in% 'Acer','hydric']+
  Betula*cors[cors$column %in% 'Betula','hydric']+
  Carya*cors[cors$column %in% 'Carya','hydric']+
  Fagus*cors[cors$column %in% 'Fagus','hydric']+
  Fraxinus*cors[cors$column %in% 'Fraxinus','hydric']+
  Larix*cors[cors$column %in% 'Larix','hydric']+
  Picea*cors[cors$column %in% 'Picea','hydric']+
  Pinus*cors[cors$column %in% 'Pinus','hydric']+
  PopLirio*cors[cors$column %in% 'PopLirio','hydric']+
  Quercus*cors[cors$column %in% 'Quercus','hydric']+
  Thuja*cors[cors$column %in% 'Thuja','hydric']+
  Tilia*cors[cors$column %in% 'Tilia','hydric']+
  Tsuga*cors[cors$column %in% 'Tsuga','hydric']+
  Ulmus*cors[cors$column %in% 'Ulmus','hydric']

plot(wetness)
writeRaster(wetness,'gis/wetness.tif', overwrite=T)
warmth <- Abies*cors[cors$column %in% 'Abies','Tg']+
  Acer*cors[cors$column %in% 'Acer','Tg']+
  Betula*cors[cors$column %in% 'Betula','Tg']+
  Carya*cors[cors$column %in% 'Carya','Tg']+
  Fagus*cors[cors$column %in% 'Fagus','Tg']+
  Fraxinus*cors[cors$column %in% 'Fraxinus','Tg']+
  Larix*cors[cors$column %in% 'Larix','Tg']+
  Picea*cors[cors$column %in% 'Picea','Tg']+
  Pinus*cors[cors$column %in% 'Pinus','Tg']+
  PopLirio*cors[cors$column %in% 'PopLirio','Tg']+
  Quercus*cors[cors$column %in% 'Quercus','Tg']+
  Thuja*cors[cors$column %in% 'Thuja','Tg']+
  Tilia*cors[cors$column %in% 'Tilia','Tg']+
  Tsuga*cors[cors$column %in% 'Tsuga','Tg']+
  Ulmus*cors[cors$column %in% 'Ulmus','Tg']
writeRaster(warmth,'gis/warmth.tif', overwrite=T)
chem <- Abies*cors[cors$column %in% 'Abies','pH50']+
  Acer*cors[cors$column %in% 'Acer','pH50']+
  Betula*cors[cors$column %in% 'Betula','pH50']+
  Carya*cors[cors$column %in% 'Carya','pH50']+
  Fagus*cors[cors$column %in% 'Fagus','pH50']+
  Fraxinus*cors[cors$column %in% 'Fraxinus','pH50']+
  Larix*cors[cors$column %in% 'Larix','pH50']+
  Picea*cors[cors$column %in% 'Picea','pH50']+
  Pinus*cors[cors$column %in% 'Pinus','pH50']+
  PopLirio*cors[cors$column %in% 'PopLirio','pH50']+
  Quercus*cors[cors$column %in% 'Quercus','pH50']+
  Thuja*cors[cors$column %in% 'Thuja','pH50']+
  Tilia*cors[cors$column %in% 'Tilia','pH50']+
  Tsuga*cors[cors$column %in% 'Tsuga','pH50']+
  Ulmus*cors[cors$column %in% 'Ulmus','pH50']
writeRaster(chem,'gis/chem.tif', overwrite=T)


all = Abies+Acer+Betula+Carya+Fagus+Fraxinus+Larix+Picea+Pinus+PopLirio+Quercus+Thuja+Tilia+Tsuga+Ulmus
plot(all)

pineoak <- (Pinus+Quercus)/all
plot(pineoak)
oakvpine <- Quercus/(Pinus+Quercus+0.001)
plot(oakvpine)
oakhickory <- (Carya+Quercus)/all
plot(oakhickory)
hickoyvoak <- Carya/(Carya+Quercus+0.001)
plot((hickoyvoak))

northvsouth <- (Tsuga+Pinus+Abies+Thuja)-(Carya+Tilia+Quercus)
plot(northvsouth>0)
hardwoods <- (wetness >= 0)*100 + (fire >= 0.50)*10 + (northvsouth >= 0)
plot(hardwoods)
writeRaster(hardwoods, 'gis/hardwoods.tif',overwrite=T)
pine <- Pinus/all
oak <- Quercus/all
beech <- Fagus/all
pineoakmix <- (pine < 2*oak)*(oak < 2*pine)*(pineoak >0.5)
plot(pineoakmix)

pinebeech <- (Pinus+Fagus)/all
plot((fire > 0.333)+(fire > 0.667))

boreal <- (vars270$hydric < 0.333)*1+((Abies+Picea+Thuja)/all >= 0.333)*10
boreal <- ((vars270$watertable >= 50) & (vars270$hydric < 0.333))*1+((Abies+Picea+Thuja)/all >= 0.333)*10

writeRaster(boreal, 'gis/boreal2.tif', overwrite=T)