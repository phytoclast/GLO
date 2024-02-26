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

pts <- readRDS('points/pts.geo1.RDS')
brick=c(topoclim,soil)
brick90 <- aggregate(brick, 3)



#Tsuga ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Tsuga', 1,0))


poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Tsuga <- predict(brick90, rf, na.rm=T)
plot(Tsuga)
writeRaster(Tsuga,'gis/Tsuga.tif', overwrite=T)

#ba
BA <- subset(pts, !is.na(BA)  & BA < 600) |> mutate(BA = ifelse(BA > 100,100,BA))
train <- extract(brick,BA)
train <- cbind(train, BA = st_drop_geometry(BA$BA))
train <- subset(train, !is.na(BA)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(BA ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


basalarea <- predict(brick90, rf, na.rm=T)
plot(basalarea)
writeRaster(basalarea,'gis/basalarea.tif', overwrite=T)

#Forest
f <- subset(pts, !is.na(BA)  & BA < 600) |> mutate(f = ifelse(BA >= 5,1,0))
train <- extract(brick,f)
train <- cbind(train, f = st_drop_geometry(f$f))
train <- subset(train, !is.na(f)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(f ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


forest <- predict(brick90, rf, na.rm=T)
plot(forest)
writeRaster(forest,'gis/forest.tif', overwrite=T)
#Wooded
f <- subset(pts, !is.na(BA)  & BA < 600) |> mutate(f = ifelse(DT >= 10,1,0))
train <- extract(brick,f)
train <- cbind(train, f = st_drop_geometry(f$f))
train <- subset(train, !is.na(f)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(f ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


wooded <- predict(brick90, rf, na.rm=T)
plot(wooded)
writeRaster(wooded,'gis/wooded.tif', overwrite=T)



library(ggplot2)
BAsum <- pts |> subset(!is.na(BA), select=c(id, BA)) |> st_drop_geometry() |> unique()
BAsum <- BAsum |> summarise(B001 = quantile(BA, 0.001),
                            B01 = quantile(BA, 0.01),
                            B25 = quantile(BA, 0.25),
                            B50 = quantile(BA, 0.50),
                            B75 = quantile(BA, 0.75),
                            B99 = quantile(BA, 0.99),
                            B999 = quantile(BA, 0.999))

BAsum <- pts |> subset(!is.na(BA) & BA < 200 & DT < 200, select=c(id, BA, DT)) |> st_drop_geometry() |> unique()
BAsum <- BAsum[sample(1:nrow(BAsum), size=nrow(BAsum)/10),]
ggplot(BAsum, aes(x=DT, y=BA))+
  geom_point(alpha=0.05)+
  geom_smooth()
#Fagus ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Fagus', 1,0))


poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Fagus <- predict(brick90, rf, na.rm=T)
plot(Fagus)
writeRaster(Fagus,'gis/Fagus.tif', overwrite=T)
#Thuja ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Thuja', 1,0))


poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Thuja <- predict(brick90, rf, na.rm=T)
plot(Thuja)
writeRaster(Thuja,'gis/Thuja.tif', overwrite=T)
#Abies ----
pts.pos <- pts |> mutate(pos= ifelse(Level2 %in% 'Abies', 1,0))


poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Abies <- predict(brick90, rf, na.rm=T)
plot(Abies)
writeRaster(Abies,'gis/Abies.tif', overwrite=T)
#Picea ----
pts.pos <- pts |> mutate(pos= ifelse(Level2 %in% 'Picea', 1,0))


poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Picea <- predict(brick90, rf, na.rm=T)
plot(Picea)
writeRaster(Picea,'gis/Picea.tif', overwrite=T)
#Pinus ----
pts.pos <- pts |> mutate(pos= ifelse(Level2 %in% 'Pinus', 1,0))


pos <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0) |> vect() |> project(soil)

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


Pinus <- predict(brick90, rf, na.rm=T)
plot(Pinus)
writeRaster(Pinus,'gis/Pinus.tif', overwrite=T)

#Pinus strobus ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Pinus strobus', 1,0))
pts.pos <- pts.pos |> mutate(pos= ifelse(Species %in% 'Pinus', 
                                     NA,pos))

pos <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0) |> vect() |> project(soil)

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


Pinus_strobus <- predict(brick90, rf, na.rm=T)
plot(Pinus_strobus)
writeRaster(Pinus_strobus,'gis/Pinus_strobus.tif', overwrite=T)
#Pinus banksiana ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Pinus banksiana', 1,0))
pts.gen <- pts.pos |> subset(Species %in% 'Pinus')

poss <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id & !id %in% pts.gen$id) |> vect() |> project(soil)

poss <- extract(brick,poss)
neg <- extract(brick,neg)
poss$pos = 1
neg$pos = 0
train <- rbind(poss,neg)
train <- subset(train, !is.na(pos)&!is.na(solar)&!is.na(toi)&!is.na(tgs)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope.500)&!is.na(aspect))

rf <- ranger(pos ~ solar+toi+toip+toin+toip100+toin100+twi+tpi+dem+
               slope+aspect+slope.500+bt+tgs+ppt+p1+p2+p3+p4+
               sand50+sand150+pH50+watertable+hydric+OM150+rockdepth+carbdepth+floodfrq, 
             data=train, sample.fraction = 0.333, num.trees=100,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 


Pinus_banksiana <- predict(brick90, rf, na.rm=T)
plot(Pinus_banksiana)
writeRaster(Pinus_banksiana,'gis/Pinus_banksiana.tif', overwrite=T)
#Pinus resinosa ----
pts.pos <- pts |> mutate(pos= ifelse(Species %in% 'Pinus resinosa', 1,0))
pts.pos <- pts.pos |> mutate(pos= ifelse(Species %in% 'Pinus', 
                                     NA,pos))

pos <- subset(pts.pos,pos %in% 1) |> vect() |> project(soil)
neg <- subset(pts.pos,pos %in% 0) |> vect() |> project(soil)

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


Pinus_resinosa <- predict(brick90, rf, na.rm=T)
plot(Pinus_resinosa)
writeRaster(Pinus_resinosa,'gis/Pinus_resinosa.tif', overwrite=T)
