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
    
    
    
    
    formular.gam <- as.formula(paste(paste("pos",paste(paste("s(",smoothvars,")", collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    
    formular.glm <- as.formula(paste(paste("pos",paste(paste("I(",smoothvars,"^2)", collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    
    formular.rf <- as.formula(paste(paste("pos",paste(paste(smoothvars, collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    
    
    rf <- ranger(formular.rf,
                 
                 data=train)
    if(WT){
      gm <- gam(formular.gam,
                family='binomial',
                data=train, weights = train$wts)
    }else{
      gm <- gam(formular.gam,
                family='binomial',
                data=train)
    }
    if(WT){
      gl <- glm(formular.glm,
                family='binomial',
                data=train, weights = train$wts)
    }else{
      gl <- glm(formular.glm,
                family='binomial',
                data=train)
    }
    
    
    
    summary(gl)
    
    mxnt <- maxnet(p=train2$pos, data= train2[,smoothvars])
    
    
    gb <- gbm(formular.rf,
              distribution = "bernoulli",
              n.trees = 200,
              bag.fraction = 0.5,
              interaction.depth=7,
              shrinkage = 0.1,
              data=train)
    
    
    
    sum(train$pos)/nrow(train)
    train.gam <- train |> mutate(prediction = gam::predict.Gam(gm, train, na.rm=T, type = "response"))
    test.gam <- test |> mutate(prediction = gam::predict.Gam(gm, test, na.rm=T, type = "response"))
    train.rf <- train |> mutate(prediction = predictions(predict(rf, train, na.rm=T)))
    test.rf <- test |> mutate(prediction = predictions(predict(rf, test, na.rm=T)))
    train.gb <- train |> mutate(prediction = gbm::predict.gbm(gb, train, na.rm=T, type = "response"))
    test.gb <- test |> mutate(prediction = gbm::predict.gbm(gb, test, na.rm=T, type = "response"))
    
    train.mxnt <- train |> mutate(prediction = predict(mxnt, as.data.frame(train), na.rm=T, type='logistic'))
    test.mxnt <- test |> mutate(prediction = predict(mxnt, as.data.frame(test), na.rm=T, type='logistic'))
    train.gl <- train |> mutate(prediction = predict.glm(gl, train, na.rm=T, type = "response"))
    test.gl <- test |> mutate(prediction = predict.glm(gl, test, na.rm=T, type = "response"))
    
    modmets <- data.frame(taxon = ttt[i],
                          model = c("GAM", "RF","GBM", "MAXNET","GLM"),
                          seed = seed[k],
                          AUCtrain = c(Metrics::auc(actual=train.gam$pos, predicted=train.gam$prediction),
                                       Metrics::auc(actual=train.rf$pos, predicted=train.rf$prediction),
                                       Metrics::auc(actual=train.gb$pos, predicted=train.gb$prediction),
                                       Metrics::auc(actual=train.mxnt$pos, predicted=train.mxnt$prediction),
                                       Metrics::auc(actual=train.gl$pos, predicted=train.gl$prediction)),
                          AUCtest = c(Metrics::auc(actual=test.gam$pos, predicted=test.gam$prediction),
                                      Metrics::auc(actual=test.rf$pos, predicted=test.rf$prediction),
                                      Metrics::auc(actual=test.gb$pos, predicted=test.gb$prediction),
                                      Metrics::auc(actual=test.mxnt$pos, predicted=test.mxnt$prediction),
                                      Metrics::auc(actual=test.gl$pos, predicted=test.gl$prediction)),
                          maxKappatrain = c(maxKappa(actual=train.gam$pos, predicted=train.gam$prediction),
                                            maxKappa(actual=train.rf$pos, predicted=train.rf$prediction),
                                            maxKappa(actual=train.gb$pos, predicted=train.gb$prediction),
                                            maxKappa(actual=train.mxnt$pos, predicted=train.mxnt$prediction),
                                            maxKappa(actual=train.gl$pos, predicted=train.gl$prediction)),
                          maxKappatest = c(maxKappa(actual=test.gam$pos, predicted=test.gam$prediction),
                                           maxKappa(actual=test.rf$pos, predicted=test.rf$prediction),
                                           maxKappa(actual=test.gb$pos, predicted=test.gb$prediction),
                                           maxKappa(actual=test.mxnt$pos, predicted=test.mxnt$prediction),
                                           maxKappa(actual=test.gl$pos, predicted=test.gl$prediction)))
    
    
    

    if(i==1 & k==1){Kfold = modmets}else{Kfold = rbind(Kfold, modmets)}
    
  }}
write.csv(Kfold, 'kfold.csv', row.names = FALSE)



correction = log(positivity)/log(0.5)#exponent to raise the predicted outcome to adjust the overall abundance to match unweighted model

#generate prediction raster ----
rf.prediction <-  predict(vars90, rf, na.rm=T);  names(rf.prediction) <- taxon
gm.prediction <-  predict(vars90, gm, na.rm=T, type = "response");  names(gm.prediction) <- taxon
mxnt.prediction <-  predict(vars90, mxnt, na.rm=T, type='logistic');  names(mxnt.prediction) <- taxon
gb.prediction <-  predict(vars90, gb, na.rm=T, type='response');  names(gb.prediction) <- taxon
gl.prediction <-  predict(vars90, gl, na.rm=T, type='response');  names(gl.prediction) <- taxon

writeRaster(rf.prediction, paste0('gis/modelcompare/',taxon,'_rf.tif'), overwrite=T)
writeRaster(gm.prediction, paste0('gis/modelcompare/',taxon,'_gam.tif'), overwrite=T)
writeRaster(gl.prediction, paste0('gis/modelcompare/',taxon,'_glm.tif'), overwrite=T)
writeRaster(gb.prediction, paste0('gis/modelcompare/',taxon,'_gb.tif'), overwrite=T)
writeRaster(mxnt.prediction, paste0('gis/modelcompare/',taxon,'_mxnt.tif'), overwrite=T)

plot(gm.prediction)
Sys.time() - timeA

