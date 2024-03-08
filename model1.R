#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- readRDS('points/pts.geo1.RDS')

solar <- rast('gis/dem/tifs/solar.tif')
openess <- rast('gis/dem/tifs/mi_dem_30_layers.tif')
names(openess) <- c('shade','popen','nopen')
slope <- rast('gis/dem/tifs/mi_dem_30_slope.tif')
names(slope) <- c('slope','aspect')
dem <- rast('gis/dem/tifs/mi_dem_30.tif')
slope500 <- rast('gis/dem/tifs/slope500.tif');names(slope500) <- 'slope500' #mean slope 500 m radius
names(dem) <- c('elev')
dem <- c(dem, slope, slope500, openess, solar)

prism <- c(rast('gis/climate/prism.repro.tif'), rast('gis/climate/Tg30.tif'))
soil <- rast('gis/soilcov/allsoilcrop.tif')

vars <- c(prism, soil, dem)
names(vars)

# pts.vars <- pts |> vect() |> project(vars)
# pts.vars <- extract(vars,pts.vars)
# pts.vars <- pts |> cbind(pts.vars)
# saveRDS(pts.vars, 'pts.vars.RDS')
# vars90 <- aggregate(vars, 3)
# writeRaster(vars90,'gis/vars90.tif', overwrite=T)

#Load data ----
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
ttt <- pts.vars |> group_by(Level2) |> summarize(npts = length(Level2))

#Taxa ----
taxon = 'Thuja'
pts.pos <- pts.vars |> mutate(pos= ifelse(Species %in% taxon, 1,0)) |> st_drop_geometry()

poss <- subset(pts.pos,pos %in% 1) 
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) 

poss$pos = 1
neg$pos = 0
train0 <- rbind(poss,neg)
train0 <- subset(train0, !is.na(pos)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))
takeout <- sample(1:nrow(train0), nrow(train0)*0.05)
train <- train0[-takeout,]
test <- train0[takeout,]
#random forest

rf <- ranger(pos ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
               Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
               hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
               slope+slope500+popen+nopen+solar, 
             data=train0, sample.fraction = 0.333, num.trees=200,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 

test <- test |> mutate(pred = predictions(predict(rf, test)))
corr <- cor(test$pos, test$pred)
corr

prediction <-  predict(vars90, rf, na.rm=T);  names(prediction) <- taxon
plot(prediction)
writeRaster(prediction, paste0('gis/models/',taxon,'.tif'), overwrite=T)



Acer <- rast('gis/models/Acer.tif')
Fagus <- rast('gis/models/Fagus.tif')
Tsuga <- rast('gis/models/Tsuga.tif')
Betula <- rast('gis/models/Betula.tif')
Tilia <- rast('gis/models/Tilia.tif')
Quercus <- rast('gis/models/Quercus.tif')
Carya <- rast('gis/models/Carya.tif')
Fraxinus <- rast('gis/models/Fraxinus.tif')
Ulmus <- rast('gis/models/Ulmus.tif')
Pinus <- rast('gis/models/Pinus.tif')
Abies <- rast('gis/models/Abies.tif')
Picea <- rast('gis/models/Picea.tif')
Thuja <- rast('gis/models/Thuja.tif')
Larix <- rast('gis/models/Larix.tif')
BA <- rast('gis/models/basalarea.tif')

Trees <- c(Acer, Fagus, Tilia, Tsuga, Betula, Fraxinus, Ulmus, Carya, Quercus, Pinus, Abies, Picea, Thuja, Larix, BA)
writeRaster(Trees, 'gis/models/Trees.tif', overwrite=T)

#gam ----
gm.full <- gam(pos ~ s(p)+
                 s(pq1)+
                 s(pq2)+
                 s(pq3)+
                 s(pq4)+
                 s(Twh)+
                 s(Tw)+
                 s(Tc)+
                 s(Tcl)+
                 s(Tg)+
                 s(e)+
                 s(m)+
                 s(Tg30)+
                 Bhs+
                 carbdepth+
                 clay150+
                 floodfrq+
                 histic+
                 humic+
                 humicdepth+
                 hydric+
                 ksatdepth+
                 OM150+
                 pH50+
                 rockdepth+
                 sand150+
                 sand50+
                 spodic+
                 watertable+
                 slope+
                 slope500+
                 popen+
                 nopen+
                 solar, 
               data=train)
gm.null <- gam(pos ~ 1, 
               data=train)

step.Gam(object = gm.full, scope= list(gm.null))
summary(gm.full)
summary(gm)
prediction <-  predict(vars90, gm.full, na.rm=T);  names(prediction) <- taxon
plot(prediction)
writeRaster(prediction, paste0('gis/',taxon,'.gam.tif'), overwrite=T)

#glm ----
gm <- glm(pos ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
            Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
            hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
            slope+slope500+popen+nopen+solar+
            p:e+Tg30:watertable+spodic:watertable+I(Tg30^2)+I(sand150^2)+I(m^2)
          , 
          data=train0)


mod <- step(gm, scope=list(~1))
summary(mod)
prediction <-  predict(vars90, mod, na.rm=T);  names(prediction) <- taxon
plot(prediction)
writeRaster(prediction, paste0('gis/',taxon,'.glm.tif'), overwrite=T)

# library(Cubist)
# covs <- train |> select(c(p, pq1, pq2, pq3, pq4, Twh, Tw, Tc, Tcl, Tg, e, m, Tg30, 
#                           Bhs, carbdepth, clay150, floodfrq, histic, humic, humicdepth, 
#                           hydric, ksatdepth, OM150, pH50, rockdepth, sand150, sand50, spodic, watertable, 
#                           slope, slope500, popen, nopen, solar))
# y <- train$pos
# 
# cmod <- cubist(covs, y)
# 
# summary(cmod)
# prediction <-  predict(vars90, cmod, na.rm=T);  names(prediction) <- taxon
# plot(prediction)
# writeRaster(prediction, paste0('gis/',taxon,'.cub.tif'), overwrite=T)

#ba ----
train <- subset(pts.pos, !is.na(BA)  & BA < 600) |> mutate(BA = ifelse(BA > 20,20,BA)) |> st_drop_geometry()

train <- subset(train, !is.na(BA)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))

rf <- ranger(BA ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
               Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
               hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
               slope+slope500+popen+nopen+solar, 
             data=train, sample.fraction = 0.333, num.trees=200,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 

basalarea <- predict(vars90, rf, na.rm=T);  names(basalarea) <- 'basalarea'
plot(basalarea)
writeRaster(basalarea,'gis/model/basalarea.tif', overwrite=T)

# gm <- glm(BA ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
#             Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
#             hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
#             slope+slope500+popen+nopen+solar, 
#           data=train)
# 
# 
# summary(gm)
# basalarea <-  predict(vars90, gm, na.rm=T);  names(basalarea) <- 'basalarea'
# plot(basalarea)
# writeRaster(basalarea, paste0('gis/model/basalarea.glm.tif'), overwrite=T)

