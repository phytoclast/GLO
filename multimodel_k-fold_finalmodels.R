#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)
library(Metrics)
library(maxnet)
library(gbm)
# library(ubl)
maxKappa <- function(actual, predicted){ for(i in 1:99){
  k <- i/100
  Kappa0 <- ModelMetrics::kappa(actual=actual, predicted=predicted, cutoff = k)
  if(i == 1){ maxkappa = k}else{maxkappa = max(Kappa0, maxkappa)}
}
  return(maxkappa)}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- readRDS('points/pts.geo1.RDS')


pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
vars270 <- rast('gis/vars270.tif')


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

n <- 25000 #total number of positive samples
maxkeep <-  5000 #total number of positive samples for maxent that runs more slowly
ntest = 0.2
i=1
#Should Random Oversampling (ROF) be employed?
ROF = T
#Should linear models be weighted back to original prevalence>
WT = F
#Should Testing data be class balanced?
ROFTEST = T

#K-Fold replicate seeds
seed <- c(2037, 4392, 6190, 6960, 8848) #high points East, CONUS, North America, Americas, World




#Taxa ----
#loop through each taxon producing a model for each
for (i in 1:length(ttt)){ 
  for (k in 1:length(seed)){
    
    
    
    #i=14;k=1
    set.seed(seed[k])# for predictability
    taxon = ttt[i]
    pts.pos <- pts.vars90 |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()
    # taxon = "Liriodendron"; pts.pos <- pts.vars90 |> mutate(pos= ifelse(Species %in% taxon, 1,0)) |> st_drop_geometry()
    
    poss <- subset(pts.pos,pos %in% 1)
    neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id)
    
    poss$pos = 1
    neg$pos = 0
    train0 <- rbind(poss,neg)
    train0 <- subset(train0, !is.na(pos)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))
    #weights ----
    train0 <- train0 |> group_by(pos) |> mutate(wts = 100000/length(pos)) |> ungroup()
    positivity <- mean(train0$pos)
    
    # train0 <- train0 |> group_by(pos) |> mutate(wts = 1) |> ungroup()#use for unweighted
    positives <- subset(train0, pos %in% 1)
    negatives <- subset(train0, pos %in% 0)
    takeout.p <- sample(1:nrow(positives), nrow(positives)*ntest)
    takeout.n <- sample(1:nrow(negatives), nrow(negatives)*ntest)
    train.p <- positives[-takeout.p,]
    test.p <- positives[takeout.p,]
    train.n <- negatives[-takeout.n,]
    test.n <- negatives[takeout.n,]
    if(ROF){#Random Over Sample
      train.p <- train.p[sample(1:nrow(train.p), n, replace = TRUE),]
      train.n <- train.n[sample(1:nrow(train.n), n, replace = TRUE),]
    }
    if(ROFTEST){#Random Over Sample
      test.p <- test.p[sample(1:nrow(test.p), n, replace = TRUE),]
      test.n <- test.n[sample(1:nrow(test.n), n, replace = TRUE),]
    }
    train <- rbind(train.p,train.n)
    test <- rbind(test.p,test.n)
    #reduce data while avoiding fewer than 50 positives if model is not already resampled
    if(nrow(train) > n*3){
      train.p <- subset(train, pos %in% 1)
      train.n <- subset(train, pos %in% 0)
      keep.p <- sample(1:nrow(train.p), pmax(n*3*(positivity), 50),replace = TRUE)
      keep.n <- sample(1:nrow(train.n), pmax(n*3*(1-positivity),50), replace = TRUE)
      train.p <- train.p[keep.p,]
      train.n <- train.n[keep.n,]
      train <- rbind(train.p, train.n)
    }
    if(WT){
      train <- train |> group_by(pos) |> mutate(wts = ifelse(pos %in% 1, positivity/length(pos), (1-positivity)/length(pos))) |> ungroup()
    }
    #reduce data further for Maxnet models
    train2 <- train
    if(nrow(train2) > maxkeep*2){
      train2.p <- subset(train2, pos %in% 1)
      train2.n <- subset(train2, pos %in% 0)
      keep.p <- sample(1:nrow(train2.p), pmax(maxkeep*2*(positivity), 50),replace = TRUE)
      keep.n <- sample(1:nrow(train2.n), pmax(maxkeep*2*(1-positivity),50), replace = TRUE)
      train2.p <- train2.p[keep.p,]
      train2.n <- train2.n[keep.n,]
      train2 <- rbind(train2.p, train2.n)
    }
    
    
    #models
    timeA <- Sys.time()
    
    smoothvars <- c('p','e','s','d','Twh','Tw','Tc','Tcl','Tg','m',
                    'Bhs','carbdepth','clay150','floodfrq','histic','humic','humicdepth',
                    'hydric','ksatdepth','OM150','pH50','rockdepth','sand150','sand50',
                    'spodic','watertable','slope','slope500','popen','nopen','solar')
    

    
    formular.rf <- as.formula(paste(paste("pos",paste(paste(smoothvars, collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    
    
    rf <- ranger(formular.rf,
                 
                 data=train)
    
#generate prediction raster ----
rf.prediction <-  predict(vars90, rf, na.rm=T);  names(rf.prediction) <- taxon

writeRaster(rf.prediction, paste0('gis/finalrun/',taxon,'.',seed[k],'_rf.tif'), overwrite=T)
  }}

Sys.time() - timeA
