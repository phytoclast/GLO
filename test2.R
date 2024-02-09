#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

topoclim <- rast('gis/topoclim.tif')
soil <- rast('gis/soil.tif')

pts <- read.csv('msb-paleon.31.0/PLS_northernMichigan_trees_Level0_v1.0.csv')
pts.geo <- st_as_sf(pts, coords=c(x='x_alb', y='y_alb'), crs=3175)
brick=c(topoclim,soil)
brick90 <- aggregate(brick, 3)
pts$SPP1 |> unique()|> sort()

pts.geo <- pts.geo |> mutate(pos= ifelse(SPP1 %in% c('bch','BCH','BE') | SPP2 %in% c('bch','BCH','BE') |SPP3 %in% c('bch','BCH','BE') |SPP4 %in% c('bch','BCH','BE') ,1,0))

pos <- subset(pts.geo,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.geo,pos %in% 0) |> vect() |> project(soil)


pos <- extract(brick,pos)
neg <- extract(brick,neg)
pos$pos = 1
neg$pos = 0
train <- rbind(pos,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 

beech <- predict(brick90, rf, na.rm=T)
plot(beech)
writeRaster(beech,'gis/beech.tif')





pts.geo <- pts.geo |> mutate(pos= ifelse(SPP1 %in% c('hem','HE','Hem') | SPP2 %in% c('hem','HE','Hem')|SPP3 %in% c('hem','HE','Hem')|SPP4 %in% c('hem','HE','Hem'),1,0))

pos <- subset(pts.geo,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.geo,pos %in% 0) |> vect() |> project(soil)

brick=c(topoclim,soil)

pos <- extract(brick,pos)
neg <- extract(brick,neg)
pos$pos = 1
neg$pos = 0
train <- rbind(pos,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

hemlock <- predict(brick90, rf, na.rm=T)
plot(hemlock)
writeRaster(hemlock,'gis/hemlock.tif')





