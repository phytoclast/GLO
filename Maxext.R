#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(maxnet)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- readRDS('points/pts.geo1.RDS')


pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')

pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
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
  ntest = 0.8
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

  
  mxnt <- maxnet(p=train$pos, data= train[c('p','e','s','d','Twh','Tw','Tc','Tcl','Tg','e','m','Tg30',
                             'Bhs','carbdepth','clay150','floodfrq','histic','humic','humicdepth',
                             'hydric','ksatdepth','OM150','pH50','rockdepth','sand150','sand50',
                             'spodic','watertable','slope','slope500','popen','nopen','solar')])
  
  
  glmnet

  
  #generate prediction raster ----
  prediction <-  predict(vars90, mxnt, na.rm=T);  names(prediction) <- taxon
  map <- ((prediction - minmax(prediction)[1])/(minmax(prediction)[2]-minmax(prediction)[1]))^2
  plot(map)
  writeRaster(map, paste0('gis/maxent_models90/',taxon,'.tif'), overwrite=T)
}
