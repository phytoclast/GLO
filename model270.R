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
# pts.vars90 <- pts |> vect() |> project(vars)
# pts.vars90 <- extract(vars90,pts.vars90)
# pts.vars90 <- pts |> cbind(pts.vars90)
# saveRDS(pts.vars90, 'pts.vars90.RDS')
# pts.vars270 <- pts |> vect() |> project(vars)
# pts.vars270 <- extract(vars270,pts.vars270)
# pts.vars270 <- pts |> cbind(pts.vars270)
# saveRDS(pts.vars270, 'pts.vars270.RDS')
vars270 <- rast('gis/vars270.tif')
# pca90 <- princomp(vars270, cor = T, maxcell=ncell(vars270)*0.1)
# pca90grid <- predict(pca90, pca90)
# writeRaster(pca90grid,'gis/pca90grid.tif', overwrite=T)
pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2

tts <- pts.vars90 |> group_by(Level2,Species) |> summarize(npts = length(Level2))
tts <- subset(tts, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10)



# pca <- princomp(vars270)

ttt[14]


#Taxa ----
#
for (i in 1:length(ttt)){ #i=14
taxon = ttt[i]
pts.pos <- pts.vars270 |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()

poss <- subset(pts.pos,pos %in% 1) 
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) 

poss$pos = 1
neg$pos = 0
train0 <- rbind(poss,neg)
train0 <- subset(train0, !is.na(pos)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))

positives <- subset(train0, pos %in% 1)
negatives <- subset(train0, pos %in% 0)
ntest = 0.5
takeout.p <- sample(1:nrow(positives), nrow(positives)*ntest)
takeout.n <- sample(1:nrow(negatives), nrow(negatives)*ntest)
train.p <- positives[-takeout.p,]
test.p <- positives[takeout.p,]
train.n <- negatives[-takeout.n,]
test.n <- negatives[takeout.n,]
train <- rbind(train.p,train.n)
test <- rbind(test.p,test.n)

#random forest

rf <- ranger(pos ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
               Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
               hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
               slope+slope500+popen+nopen+solar, 
             data=train, sample.fraction = 0.333, num.trees=200, max.depth = NULL, importance = 'impurity',
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
roc <- pROC::roc(test$pos,test$pred)
roc$auc#area under curve
1-2*sum((test$pos-(test$pred))^2)/nrow(test)#reverse brier score
# corr <- cor(test$pos, test$pred)
# corr

#loop to find probability with maximum kappa ----
df=data.frame(ii=0,kappa=0)
for(i in 1:99){#i=1
  ii=i/100
test2 <- test[,c('pos','pred')]
test2 <- test2 |> mutate(obs =  factor(ifelse(pos %in% 1, 'present','absent')), pred1 =  factor(ifelse(pred >= ii, 'present','absent')))
if(length(unique(test2$pred1))>1){
xtab <- table(test2$obs, test2$pred1)
agree <- xtab[1,1]/sum(xtab)+xtab[2,2]/sum(xtab)
chance <- sum(xtab[1,])/sum(xtab[,])*sum(xtab[,1])/sum(xtab[,])+
  sum(xtab[2,])/sum(xtab[,])*sum(xtab[,2])/sum(xtab[,])
kappa = (agree-chance)/(1-chance)
df0=data.frame(ii=ii,kappa=kappa)
df=rbind(df,df0)}
}
maxkappa <- subset(df, kappa == max(kappa))
maxkappa$ii; maxkappa$kappa

# #brier score ----
# 2*sum((test$pos-(test$pred>=maxkappa$ii))^2)/nrow(test)#multiplied by 2 as it is both presence and absence are the same
# 2*sum((test$pos-(test$pred>=0.5))^2)/nrow(test)




#generate prediction raster ----
prediction <-  predict(vars270, rf, na.rm=T);  names(prediction) <- taxon
plot(prediction)
# p2 <- ifel(prediction >= 0.26,1,0)
# plot(p2)
# writeRaster(p2, paste0('gis/models270/p2.tif'), overwrite=T)
writeRaster(prediction, paste0('gis/models270/',taxon,'.tif'), overwrite=T)
}


#ba ----
#
pts.pos <- pts.vars270 |> st_drop_geometry()

train <- subset(pts.pos, !is.na(BA)  & BA < 600 & !dataset %in% 'OH') |> mutate(BA = ifelse(BA > 20,20,BA)) 

train <- subset(train, !is.na(BA)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))

rf <- ranger(BA ~ p+pq1+pq2+pq3+pq4+Twh+Tw+Tc+Tcl+Tg+e+m+Tg30+
               Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
               hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
               slope+slope500+popen+nopen+solar, 
             data=train, sample.fraction = 0.333, num.trees=200,  importance = 'impurity',
             classification=FALSE,  write.forest = TRUE)

vimp <- data.frame(imp = rf$variable.importance) |> mutate(var = names(rf$variable.importance))  |> arrange(by=imp) 

basalarea <- predict(vars270, rf, na.rm=T);  names(basalarea) <- 'basalarea'
plot(basalarea)
writeRaster(basalarea,'gis/models270/basalarea.tif', overwrite=T)