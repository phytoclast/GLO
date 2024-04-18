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
pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2

tts <- pts.vars90 |> group_by(Level2,Species) |> summarize(npts = length(Level2))
tts <- subset(tts, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10)

sandypts <- pts.vars90 |> subset((sand50 > 70| sand150 > 80) & Tg >= 16)

for(i in 1:length(ttt)){
  pts.vars90 <- pts.vars90 |> mutate(x = ifelse(Level2 %in% ttt[i], 1,0))
  colnames(pts.vars90)[colnames(pts.vars90) %in% 'x'] <- ttt[i]
}

corvars <- subset(pts.vars90, !is.na(BA)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect), select = -c(shade,aspect))

vars  <- colnames(corvars)[grep("^p$",colnames(corvars)):grep("^solar$",colnames(corvars))]
corvars <- corvars[,c(ttt,'BA',vars)] |> sf::st_drop_geometry()

cors <- as.data.frame(cor(corvars))


c('Quercus', 'Acer', 'Fagus', 'Tsuga', 'Thuja', 'Pinus', 'Abies', 'Picea', 'Carya')
ggplot(subset(pts.vars90, Level2 %in% c('Quercus', 'Acer', 'Fagus', 'Tsuga', 'Thuja', 'Pinus', 'Abies', 'Picea', 'Carya')))+
  geom_density(aes(x=Twh, fill=Level2), alpha=0.2)
  
  