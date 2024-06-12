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


pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
# vars270 <- vars90 |> aggregate(fact=3, fun="mean", na.rm=T)
# writeRaster(vars270,'gis/vars270.tif', overwrite=T)
vars270 <- rast('gis/vars270.tif')
#supplement MNFI openings and extract variables
# mnfi <- read_sf('C:/a/Ecological_Sites/GIS/Vegetation/Veg1800/MichLP_1800veg.shp')
# mnfi <- mnfi |> subset(COVERTYPE %in% c('EXPOSED BEDROCK', 'SHRUB SWAMP/EMERGENT MARSH', 'WET PRAIRIE', 'GRASSLAND','SAND DUNE'))
# mnfi <- st_transform(mnfi, crs(pts))
# openarea <- st_area(mnfi) |> sum() |> as.numeric()
# nopenpoints <- 3*openarea/1600^2
# mnfi.vect <- vect(mnfi)
# openpoints <- spatSample(mnfi.vect, nopenpoints)
# openpoints <- st_as_sf(openpoints)
# openpoints <- openpoints[,'geometry']
# openpoints <- openpoints |> mutate(Level2 = 'open', Species = 'open', BA = 0, DT= 0)
# pts <- pts |> bind_rows(openpoints)
# 
# pts.vars90 <- pts |> vect() |> project(vars90)
# pts.vars90 <- extract(vars90,pts.vars90)
# pts.vars90 <- pts |> cbind(pts.vars90)
# saveRDS(pts.vars90, 'pts.vars90.RDS')
# pts.vars270 <- pts |> vect() |> project(vars270)
# pts.vars270 <- extract(vars270,pts.vars270)
# pts.vars270 <- pts |> cbind(pts.vars270)
# saveRDS(pts.vars270, 'pts.vars270.RDS')

# pca90 <- princomp(vars270, cor = T, maxcell=ncell(vars270)*0.1)
# pca90grid <- predict(vars270, pca90)
# writeRaster(pca90grid,'gis/pca90grid.tif', overwrite=T)

pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
#extract a list of taxon names
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2

tts <- pts.vars90 |> group_by(Level2,Species) |> summarize(npts = length(Level2))
tts <- subset(tts, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10)


ttt[14]


#Taxa ----
#loop through each taxon producing a model for each
for (i in 1:length(ttt)){ #i=14
  set.seed(4345)# for reproducability
  taxon = ttt[i]
  pts.pos <- pts.vars90 |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()
  # taxon = "Liriodendron"
  # pts.pos <- pts.vars90 |> mutate(pos= ifelse(Species %in% taxon, 1,0)) |> st_drop_geometry()
  
  poss <- subset(pts.pos,pos %in% 1) 
  neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) 
  
  poss$pos = 1
  neg$pos = 0
  train0 <- rbind(poss,neg)
  train0 <- subset(train0, !is.na(pos)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))
  #weights ----
  train0 <- train0 |> group_by(pos) |> mutate(wts = 100000/length(pos)) |> ungroup()
  
  positives <- subset(train0, pos %in% 1)
  negatives <- subset(train0, pos %in% 0)
  ntest = 0.1
  takeout.p <- sample(1:nrow(positives), nrow(positives)*ntest)
  takeout.n <- sample(1:nrow(negatives), nrow(negatives)*ntest)
  train.p <- positives[-takeout.p,]
  test.p <- positives[takeout.p,]
  train.n <- negatives[-takeout.n,]
  test.n <- negatives[takeout.n,]
  train <- rbind(train.p,train.n)
  test <- rbind(test.p,test.n)
  saveRDS(test, paste0('test/',ttt[i],'.RDS'))
  #generalize additive model

  
  gm <- gam(pos ~ p+e+s+d+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
              Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
              hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
              slope+slope500+popen+nopen+solar, 
               data=train)

  summary.Gam(gm)

  
  
  
  #generate prediction raster ----
  prediction <-  predict(vars90, gm, na.rm=T);  names(prediction) <- taxon
  plot(prediction)
  writeRaster(prediction, paste0('gis/gam_models90/',taxon,'.tif'), overwrite=T)
}
