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

n <- 2000 #total number of positive samples
ntest = 0.2
#Taxa ----
#loop through each taxon producing a model for each
for (i in 1:length(ttt)){ #i=14
  set.seed(4345)# for predictability
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
  # train0 <- train0 |> group_by(pos) |> mutate(wts = 1) |> ungroup()#use for unweighted
  positives <- subset(train0, pos %in% 1)
  negatives <- subset(train0, pos %in% 0)
  start.p <- positives[sample(1:nrow(positives), n, replace = TRUE),]
  start.n <- negatives[sample(1:nrow(negatives), n, replace = TRUE),]
  takeout.p <- sample(1:nrow(start.p), nrow(start.p)*ntest)
  takeout.n <- sample(1:nrow(start.n), nrow(start.n)*ntest)
  train.p <- start.p[-takeout.p,]
  test.p <- start.p[takeout.p,]
  train.n <- start.n[-takeout.n,]
  test.n <- start.n[takeout.n,]
  train <- rbind(train.p,train.n)
  test <- rbind(test.p,test.n)
  # takeout.p2 <- sample(1:nrow(train.p), nrow(train.p)*0.7)
  # takeout.n2 <- sample(1:nrow(train.n), nrow(train.n)*0.7)
  # train.p2 <- train.p[-takeout.p2,]
  # train.n2 <- train.n[-takeout.n2,]
  # train2 <- rbind(train.p2,train.n2)
  # test <- test |> group_by(pos) |> mutate(wts = 1/length(pos)) |> ungroup()
  # overtest <- test[sample(x = 1:nrow(test), size = nrow(test)*10, replace=TRUE, prob = test$wts),]
  # overtrain <- train[sample(x = 1:nrow(train), size = nrow(train)*5, replace=TRUE, prob = train$wts),]
  # overtrain2 <- train[sample(x = 1:nrow(train), size = 25000, replace=TRUE, prob = train$wts),]
  
  # saveRDS(test, paste0('test/',ttt[i],'.RDS'))
  #generalize additive model
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

  gm <- gam(formular.gam,
            family='binomial',
            data=train)
  summary(gm)
  
  gl <- glm(formular.glm,
            family='binomial',
            data=train)
  summary(gl)

  mxnt <- maxnet(p=train$pos, data= train[,smoothvars])
  

  gb <- gbm(formular.rf,
            distribution = "bernoulli",
            n.trees = 200,
            bag.fraction = 0.5,
            interaction.depth=7,
            shrinkage = 0.1,
            data=train)
  
  
  # rf <- ranger(pos ~ p+e+s+d+Twh+Tw+Tc+Tcl+Tg+m+
  #                Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
  #                hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
  #                slope+slope500+popen+nopen+solar,
  #              data=train, sample.fraction = 1, num.trees=200, max.depth = NULL, importance = 'impurity',
  #              classification=FALSE,  write.forest = TRUE)#, case.weights = train$wts)
  # 
  # mxnt <- maxnet(p=overtrain2$pos, data= overtrain2[c('p','e','s','d','Twh','Tw','Tc','Tcl','Tg','m',
  #                                             'Bhs','carbdepth','clay150','floodfrq','histic','humic','humicdepth',
  #                                             'hydric','ksatdepth','OM150','pH50','rockdepth','sand150','sand50',
  #                                             'spodic','watertable','slope','slope500','popen','nopen','solar')])
  
  # gl <- glm(pos ~ p+e+s+d+Twh+Tw+Tc+Tcl+Tg+m+
  #             Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
  #             hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
  #             slope+slope500+popen+nopen+solar,
  #           data=train)#, weights = train$wts)
  

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
  
  modmets <- data.frame(model = c("GAM", "RF","GBM", "MAXNET","GLM"),
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
  
  
  
  ppp <- pdp::partial(rf, pred.var = 'spodic')
  pdp::plotPartial(ppp)
  
  caret::confusionMatrix(reference=train.rf$pos, data=train.rf$prediction)
  
  

  
  
  #generate prediction raster ----
  rf.prediction <-  predict(vars90, rf, na.rm=T);  names(rf.prediction) <- taxon
  gm.prediction <-  predict(vars90, gm, na.rm=T, type = "response");  names(gm.prediction) <- taxon
  mxnt.prediction <-  predict(vars90, mxnt, na.rm=T, type='logistic');  names(mxnt.prediction) <- taxon
  gb.prediction <-  predict(vars90, gb, na.rm=T, type='response');  names(gb.prediction) <- taxon
  gl.prediction <-  predict(vars90, gl, na.rm=T, type='response');  names(gl.prediction) <- taxon
  # prediction <- (prediction - minmax(prediction)[1])/(minmax(prediction)[2]-minmax(prediction)[1]+0.00001)

  writeRaster(rf.prediction, paste0('gis/models90/',taxon,'.tif'), overwrite=T)
  writeRaster(gm.prediction, paste0('gis/gam_models90/',taxon,'splines2.tif'), overwrite=T)
  writeRaster(mxnt.prediction, paste0('gis/maxent_models90/',taxon,'.tif'), overwrite=T)
  writeRaster(gb.prediction, paste0('gis/gb_models90/',taxon,'.tif'), overwrite=T)
  plot(gm.prediction)
  Sys.time() - timeA
}
