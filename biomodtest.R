# devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
library(biomod2)
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
pcagrid <- rast('gis/pca.tif')
pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2
i=14
taxon = ttt[i]
# taxon = 'Liriodendron'

pts.pos <- pts.vars270 |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()
# pts.pos <- pts.vars270 |> mutate(pos= ifelse(Species %in% taxon, 1,0)) |> st_drop_geometry()

poss <- subset(pts.pos,pos %in% 1) 
neg <- subset(pts.pos,pos %in% 0 & !id %in% poss$id) 

poss$pos = 1
neg$pos = 0
train0 <- rbind(poss,neg)
train0 <- subset(train0, !is.na(pos)&!is.na(solar)&!is.na(popen)&!is.na(Tg30)&!is.na(watertable)&!is.na(rockdepth)&!is.na(OM150)&!is.na(sand50)&!is.na(slope)&!is.na(slope500)&!is.na(aspect))

positives <- subset(train0, pos %in% 1)
negatives <- subset(train0, pos %in% 0)
ntest = 0.8
takeout.p <- sample(1:nrow(positives), nrow(positives)*ntest)
takeout.n <- sample(1:nrow(negatives), nrow(negatives)*ntest)

# takeout.n <- sample(1:nrow(negatives), nrow(negatives)*0.90)
train.p <- positives[-takeout.p,]
# train.p <- positives[,]
test.p <- positives[takeout.p,]
train.n <- negatives[-takeout.n,]
test.n <- negatives[takeout.n,]
train <- rbind(train.p,train.n)
test <- rbind(test.p,test.n)

spts <- st_as_sf(train, coords	= c(x='lon', y='lat'), crs='EPSG: 4326')
spts_trans <- st_transform(spts, crs=crs(vars270))
train <- train |> mutate(X = st_coordinates(spts_trans)[,1], Y = st_coordinates(spts_trans)[,2])
# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = train$pos,
                                     expl.var = pcagrid,
                                     resp.xy = train[,c('X','Y')],
                                     resp.name = taxon)



myBiomodData
plot(myBiomodData)






# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    models = c("MAXNET", "RF","GAM"),#"GBM", "GLM", "SRE", "XGBOOST", "ANN", "CTA", "FDA","MAXENT", "MARS"),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 3,
                                    metric.eval = c('KAPPA','ROC'))
# seed.val = 123)
# nb.cpu = 8)
myBiomodModelOut@models.evaluation



myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = pcagrid,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)



x1 <- load('D:/scripts/GLO/Tsuga/Tsuga.1712953122.models.out')
x2 <- load('D:/scripts/GLO/Tsuga/Tsuga.1713214456.models.out')
load('D:/scripts/GLO/Tsuga/Tsuga.1712953122.models.out')
load('D:/scripts/GLO/Tsuga/Tsuga.1713214456.models.out')
x1 <- get(x1)
x2 <- get(x2)


get_evaluations(x1)
get_evaluations(x2)

myBiomodDataRAW <- BIOMOD_FormatingData(resp.var = train$pos,
                                     expl.var = vars270,
                                     resp.xy = train[,c('X','Y')],
                                     resp.name = taxon)
myBiomodDataPCA <- BIOMOD_FormatingData(resp.var = train$pos,
                                     expl.var = pcagrid,
                                     resp.xy = train[,c('X','Y')],
                                     resp.name = taxon)

myBiomodModelRaw <- BIOMOD_Modeling(bm.format = myBiomodDataRAW,
                                    models = c("MAXNET", "RF","GAM"),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 1,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 1,
                                    metric.eval = c('KAPPA','ROC'))
myBiomodModelPCA <- BIOMOD_Modeling(bm.format = myBiomodDataPCA,
                                    models = c("MAXNET", "RF","GAM"),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 1,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 1,
                                    metric.eval = c('KAPPA','ROC'))

get_evaluations(myBiomodModelRaw)
get_evaluations(myBiomodModelPCA)


