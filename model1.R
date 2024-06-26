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
vars90 <- aggregate(vars, 3, na.rm=T)
writeRaster(vars90,'gis/vars90.tif', overwrite=T)

#Load data ----
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)
library(biomod2)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
ttt <- pts.vars |> group_by(Level2) |> summarize(npts = length(Level2))
# 
# ggplot(subset(pts.vars, diam > 0), aes(x=diam))+
#   geom_density(fill='red', alpha=0.5)+
#   scale_x_sqrt(breaks=c((0:9)*5,(5:30)*10), name='diameter (cm)')+
#   scale_y_continuous(name='relative frequency')+
#   coord_cartesian(xlim=c(0,150))
# ggplot(subset(pts.vars, dists > 0), aes(x=dists))+
#   geom_density(fill='red', alpha=0.5)+
#   scale_x_sqrt(breaks=c(0:9,c(1:20)*10), name='distance (m)')+
#   scale_y_continuous(name='relative frequency')+
#   coord_cartesian(xlim=c(0,100))
# 
# median(subset(pts.vars, diam > 0)$diam)
# median(subset(pts.vars, dists > 0)$dists)

#Taxa ----
taxon = 'Carya'
pts.pos <- pts.vars |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()

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
             data=train0, sample.fraction = 0.333, num.trees=200, max.depth = NULL, importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 
# library(randomForest)
# rf <- randomForest(pos ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
#                Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
#                hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
#                slope+slope500+popen+nopen+solar, 
#              data=train0, sampsize = 0.333*nrow(train0), ntree=200,  importance = TRUE)
# 
# varImpPlot(rf)

# gm <- gam(pos ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+s(m)+s(Tg30)+
#                Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
#                hydric+ksatdepth+OM150+s(pH50)+rockdepth+s(sand150)+sand50+spodic+s(watertable)+
#                slope+slope500+popen+nopen+solar+spodic:watertable, 
#              data=train0)
# 
# summary.Gam(gm)
# 
# test <- test |> mutate(pred = (predict.Gam(gm, test)))
# corr <- cor(test$pos, test$pred)
# corr

test <- test |> mutate(pred = predictions(predict(rf, test)))
corr <- cor(test$pos, test$pred)
corr

prediction <-  predict(vars90, rf, na.rm=T);  names(prediction) <- taxon
plot(prediction)
writeRaster(prediction, paste0('gis/models/',taxon,'.tif'), overwrite=T)
#rpart ----
library(rpart)
library(rpart.plot)
rp <- rpart(pos ~ Twh+Tw+Tc+Tcl+Tg+m+Tg30+
               Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
               hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
               slope+slope500+popen+nopen+solar, 
             data=train0, control =  rpart.control(maxdepth = 5, cp=0.001, minsplit=10))
rpart.plot(rp, digits = 3)


#bring together ----
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


Fagus2 <- (Fagus/(minmax(Fagus)[2]))
Tsuga2 <- (Tsuga/(minmax(Tsuga)[2]))^0.5
Acer2 <- (Acer/(minmax(Acer)[2]))^0.5
plot(Tsuga2 > 0.5)

northmesicforest <- c(Acer2, Fagus2, Tsuga2)
writeRaster(northmesicforest,'gis/models/northmesicforest.tif', overwrite=T)

Picea2 <- (Picea/(minmax(Picea)[2]))^0.333
Abies2 <- (Abies/(minmax(Abies)[2]))^0.333
Thuja2 <- (Thuja/(minmax(Thuja)[2]))^0.333
plot(Abies2 > 0.5)

borealforest <- c(Abies2, Picea2, Thuja2)
writeRaster(borealforest,'gis/models/borealforest.tif', overwrite=T)

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
#rpart ----
library(rpart)
library(rpart.plot)
rp <- rpart(BA ~ Twh+Tw+Tc+Tcl+Tg+m+Tg30+
              Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
              hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
              slope+slope500+popen+nopen+solar, 
            data=train0, control =  rpart.control(maxdepth = 5, cp=0.001, minsplit=10))
rpart.plot(rp, digits = 3)

#sample multi models ----
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)
library(biomod2)
library(vegan)

vars90 <- rast('gis/vars90.tif')
trees <- rast('gis/models/Trees.tif')
mlra <- st_read('C:/a/Ecological_Sites/GIS/Ecoregion/MLRA_2017.shp')
mlra <- st_transform(mlra, crs(vars90))
selectedmlra <- mlra |> subset(MLRA %in% c('97','98','99')) |> vect()
sampts <- terra::spatSample(selectedmlra, size=10000, "random")

env0 <- extract(vars90, sampts)
spp0 <- extract(trees, sampts)

#move basal area to environment
env <- env0 |> cbind(subset(spp0, select = c(basalarea)))
spp <- subset(spp0, select = -c(basalarea))
comb <- cbind(env, spp)

#remove na rows
keep <- !is.na(comb$Tg) & !is.na(comb$Tsuga)
env <- env[keep,] |> subset(select=-c(ID)) 
spp <- spp[keep,] |> subset(select=-c(ID))
keep <- spp |> apply(MARGIN = 1, FUN='sum')
keep <- keep > 0
spp <- spp[keep,] 
env <- env[keep,] 
#cca 
mod <- vegan::cca(spp, env)

plot(mod, type="n", scaling="sites")
text(mod, dis="cn", scaling="sites")
#points(mod, pch=21, col="red", bg="yellow", cex=1.2, scaling="sites")
text(mod, "species", col="blue", cex=0.8, scaling="sites")