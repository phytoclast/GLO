#############################
#Western Great Lakes RF models
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(Metrics)
library(climatools)
####
library(rgbif)
library(vegnasis)
# bonapfips <- readRDS("C:/scripts/bonap2020/data/final/b.fips.RDS")
# isleroyale <- c(26083)
# beaverisland <- c(26029)
# neohio <- c(18074, 18084, 18088, 18119, 18156, 18175, 18177, 18216, 18255, 18276, 18310, 18342)
# allfips <- list(isleroyale,beaverisland,neohio)
# area1 <- c(-89.5,-88.3, 47.7, 48.2)#Isle Royale
# area2 <- c(-85.7,-85.2, 45.5, 45.89)#Beaver Island
# area3 <- c(-82,-80, 40, 42)#NE Ohio
# allareas <- list(area1,area2,area3)
# targetcols <- c("key","scientificName",
#   "decimalLatitude","decimalLongitude",
#   "continent","stateProvince","year","month","day","basisOfRecord",
#   "occurrenceStatus","taxonKey","name","institutionCode","recordedBy","identifiedBy")
# 
# 
# alld <- NULL
# 
# 
# 
# for(k in 1:3){#k=1
# thesefips <- allfips[[k]]
# thesespp <- bonapfips |> subset(FIPS %in% thesefips)
# thesespp <- thesespp |> mutate(genus = extractTaxon(ac.binomial, 'genus')) |> subset(Status %in% "N", select=c(genus, ac.binomial)) |> unique()
# thesegenera <- c("Abies","Acer","Aesculus","Betula","Carya","Celtis","Cornus",
#                  "Fagus","Fraxinus","Juglans","Larix","Liquidambar","Liriodendron","Nyssa",      
#                  "Ostrya","Carpinus","Picea","Pinus","Platanus","Populus","Prunus","Quercus","Salix",
#                  "Thuja","Tilia","Tsuga","Ulmus")
# thesespp <- thesespp |> subset(genus %in% thesegenera) |> mutate(habit.code= get.habit.code(ac.binomial)) |> subset((grepl("^T2",habit.code) | genus %in% c('Carpinus','Ostrya')) & !ac.binomial %in% 'Acer negundo')
# 
# 
# thisarea <- allareas[[k]]
# lat <- paste0(as.character(thisarea[3]), ",",as.character(thisarea[4]))
# lon <- paste0(as.character(thisarea[1]), ",",as.character(thisarea[2]))
# for (i in 1:nrow(thesespp)){#i=1
#   taxon <- thesespp$ac.binomial[i]
#   try({   #suppress server error and move on to next taxon
#     biogeopts2 <- occ_search(limit=1000,
#                              # phylumKey = 7707728, 
#                              scientificName = taxon,
#                              # gadmGid="USA.23_1",
#                              decimalLatitude=lat, decimalLongitude =lon)
#   })
#   d <- biogeopts2
#   if(!is.null(d$data)){
#   d <- d$data 
#   missingcols <- targetcols[!targetcols %in% colnames(d)]
#   if(length(missingcols) >0){d[missingcols] <- NA}
#   d <- d |> select("key","scientificName",
#                    "decimalLatitude","decimalLongitude",
#                    "continent","stateProvince","year","month","day","basisOfRecord",
#                    "occurrenceStatus","taxonKey","name","institutionCode","recordedBy","identifiedBy")
#   if(is.null(alld)){alld <- d}else{alld <- rbind(alld,d)}
#   }
# }}

#saveRDS(alld, 'C:/scripts/GLO/data/gbifdata.RDS')
# vct <- d |> vect( geom=c("decimalLongitude", "decimalLatitude"),crs=crs('epsg:4326')) 
# library(mapview)
# mapview(vct)




#Begin Models ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pts <- readRDS('points/pts.geo1.RDS')
covr <- rast("newcovariates/covr.tif")
# pts.vars <- pts |> cbind(extract(covr, project(vect(pts), covr)))
# saveRDS(pts.vars,'pts.vars.new.RDS')
pts.vars90 <- readRDS('pts.vars.new.RDS')
# xyrast <- covr$elev |> aggregate(3) |> climatools::makexyrast(6)
# xyrast <- xyrast |> project(covr)
# writeRaster(xyrast, "newcovariates/xyrast.tif")
xyrast <- rast("newcovariates/xyrast.tif")
names(xyrast[[1]]) <- 'lat0';names(xyrast[[2]]) <- 'lon0'
covs <- c(names(covr), names(xyrast))
covrast <- c(covr, xyrast)
pxy <- 0.1 #relative strength of rotatedxy covs
covwts <- c(c(1:length(names(covr)))*0+1, (c(1:length(names(xyrast)))*0+1)*length(names(xyrast))/length(names(covr))*pxy)
covwts <- covwts/max(covwts)#ensure it doesn't exceed 1

# covrast270 <- aggregate(covrast, 3)
covrast270 <- rast("newcovariates/covrast270.tif")
# writeRaster(covrast270, "newcovariates/covrast270.tif")
# 
# pts.vars90xy <- extract(xyrast, vect(pts.vars90)) 
# saveRDS(pts.vars90xy,'pts.vars.xyrast.RDS')
pts.vars90xy <- readRDS('pts.vars.xyrast.RDS')
pts.vars90 <- pts.vars90 |> cbind(pts.vars90xy)
#gbifdata ---- 
gbifdata <- readRDS('data/gbifdata.RDS')
# newcovr <- crop(covrast270, c(1100000, 1403971, 1900000, 2300000))
# newcovr <- crop(covrast270, c(500000, 600000, 2750000, 2850000))
gbifdata$Level2 <- extractTaxon(gbifdata$name, 'genus')
gbifdata$Level2 <- ifelse(gbifdata$Level2 %in% c('Ostrya', 'Carpinus'),'OstryaCarpinus',gbifdata$Level2)
gbifdata <- gbifdata |> st_as_sf(coords=c('decimalLongitude', 'decimalLatitude'), crs=crs('epsg:4326')) |> st_transform(crs(covrast))
gbifcovs <- extract(covrast, vect(gbifdata)) 
gbifdata <- gbifdata |> cbind(gbifcovs)
gbifdata <- gbifdata |> subset(!is.na(sand150) & !is.na(elev) & !is.na(popen) & !is.na(lat1))
# n=25000
# taxon = 'Picea'
# train <- gbifdata |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()
# # plot(crop(covrast270$elev, c(500000, 600000, 2750000, 2850000)))
# train <- train |> group_by(pos) |> mutate(wts = 1000/(length(pos)+1)) |> ungroup()
# sample.fraction <- pmin(1,2*n/nrow(train))
# #i=18;k=1
# set.seed(1)# for predictability
# 
# 
# 
# 
# formular.rf <- as.formula(paste(paste("pos",paste(paste(covs, collapse = " + ", sep = ""),""), sep = " ~ ")
# ))
# 
# rf <- ranger(formular.rf,
#              case.weights = train$wts,
#              split.select.weights = covwts,
#              sample.fraction=sample.fraction,
#              data=train)
# rf.prediction <-  predict(newcovr, rf, na.rm=T);  names(rf.prediction) <- taxon
# plot(rf.prediction)


#----
mnfi <- read_sf('C:/GIS/Vegetation/Veg1800/MichLP_1800veg.shp')
mnfi <- read_sf("C:/GIS/Vegetation/presettlement2026/Land_Cover_Circa_1800.shp")
mnfi <- mnfi |> subset(COVERTYPE %in% c('EXPOSED BEDROCK', 'GRASSLAND', 'WET PRAIRIE', 'SHRUB SWAMP/EMERGENT MARSH','SAND DUNE'), select='COVERTYPE')
mnfi <- st_transform(mnfi, crs = st_crs(pts.vars90))
mnfi <- extract(vect(mnfi),vect(pts.vars90)) |> unique()
pts.vars90 <- pts.vars90 |> cbind(mnfi)

pts.vars90 <- pts.vars90 |> mutate(Level2 = case_when(COVERTYPE %in% c('EXPOSED BEDROCK','SAND DUNE')  ~ 'zBare',
                                                      (((Species %in% 'no tree' &  is.na(COVERTYPE))|(dists >= 250 &  is.na(COVERTYPE))) & (watertable >=25)) | COVERTYPE %in% c('GRASSLAND','WET PRAIRIE') ~ 'zPrairie',
                                                      (((Species %in% 'no tree' &  is.na(COVERTYPE))|(dists >= 250 &  is.na(COVERTYPE))) & (watertable < 25) ) | COVERTYPE %in% c('SHRUB SWAMP/EMERGENT MARSH') ~ 'zMArsh',
                                                      Level2 %in% 'Platanus' & lat >= 44 ~ 'unk',
                                                      TRUE ~ Level2))
                                                      
pts.vars90 <- pts.vars90 |> mutate(id.y=NULL, ID=NULL)

# st_write(pts.vars90, 'temppts1.shp', append = F)

#solve Liro
poplirio <- subset(pts.vars90, Level2 %in% 'PopLirio')|> st_drop_geometry()
unk <- subset(poplirio, Species %in% 'PopLirio')
known <-  subset(poplirio, !Species %in% 'PopLirio')
known <- known |> mutate(spp = as.factor(ifelse(Species %in% 'Liriodendron', "Liriodendron","Populus")))


formular.rf <- as.formula(paste(paste("spp",paste(paste(covs, collapse = " + ", sep = ""),""), sep = " ~ ")
))


known <- known |> group_by(spp) |> mutate(wts = 100000/length(spp)) |> ungroup()


rf <- ranger(formular.rf,
             case.weights = known$wts,
             data=known)

unk <- unk |> mutate(spp = predictions(predict(rf, unk)))
unk <- unk[,c("dataset","id","treeseq","tname","spp")]

pts.vars90a <-  pts.vars90 |> left_join(unk)

pts.vars90a <- pts.vars90a |> mutate(Level2 = case_when(Species %in% "PopLirio" & !is.na(spp) ~ spp,
                                                        Level2 %in% "PopLirio" & Species %in% "Liriodendron" ~ "Liriodendron",
                                                        Level2 %in% "PopLirio" ~ "Populus",
                                                        Species %in% "Thuja" ~ "Thuja",
                                                        TRUE ~ Level2))

#solve gums
gums <- subset(pts.vars90, Level2 %in% 'NyssaLiquidambar')|> st_drop_geometry()
unk <- subset(gums, Species %in% 'NyssaLiquidambar')
known <-  subset(gums, !Species %in% 'NyssaLiquidambar')
known <- known |> mutate(spp = as.factor(ifelse(Species %in% 'Liquidambar', "Liquidambar","Nyssa")))
formular.rf <- as.formula(paste(paste("spp",paste(paste(covs, collapse = " + ", sep = ""),""), sep = " ~ ")
))


known <- known |> group_by(spp) |> mutate(wts = 100000/length(spp)) |> ungroup()


rf <- ranger(formular.rf,
             case.weights = known$wts,
             data=known)

unk <- unk |> mutate(spp = predictions(predict(rf, unk)))
unk <- unk[,c("dataset","id","treeseq","tname","spp")]

pts.vars90a <-  pts.vars90a |> mutate(spp=NULL) |> left_join(unk)

pts.vars90a <- pts.vars90a |> mutate(Level2 = case_when(Species %in% "NyssaLiquidambar" & !is.na(spp) ~ spp,
                                                        Level2 %in% "NyssaLiquidambar" & Species %in% "Liquidambar" ~ "Liquidambar",
                                                        Level2 %in% "NyssaLiquidambar" ~ "Nyssa",
                                                        TRUE ~ Level2))
pts.vars90a <-  pts.vars90a |> mutate(spp=NULL)


#merge with gbif data
pts.vars90a$origin <- 'glo'
gbifdata$origin <- 'gbif'
newcols <- intersect(colnames(pts.vars90a), colnames(gbifdata))

#reweigh gbif species due to increase in abundance in disturbance species
genwts1 <- pts.vars90a|>st_drop_geometry() |> mutate(tots = length(Level2)) |> group_by(Level2, tots) |> summarise(gentots = length(Level2)) |> ungroup() |> mutate(gentots = round(100*gentots/tots,3), tots=NULL) 
genwts2 <- gbifdata |>st_drop_geometry() |> mutate(tots = length(Level2)) |> group_by(Level2, tots) |> summarise(gentots2 = length(Level2)) |> ungroup() |> mutate(gentots2 = round(100*gentots2/tots,3), tots=NULL) 
genwts1 <- genwts1 |> left_join(genwts2) |> mutate(gratio = gentots/gentots2)

pts.vars90b <- st_drop_geometry(pts.vars90a[,newcols]) |> rbind(st_drop_geometry(gbifdata[,newcols]))
pts.vars90b <- pts.vars90b |> left_join(genwts1[,c('Level2', 'gratio')]) 

# st_write(pts.vars90a, 'temppts2.shp', append = F)
ttt <- pts.vars90b |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('rock','unk','no tree','no data') & npts >= 1000 | grepl('^z',Level2))
ttt <- ttt$Level2






startwith <- 1
n <- 25000 #total number of positive samples
timeA <- Sys.time()
#Taxa ----
#loop through each taxon producing a model for each
for (i in startwith:length(ttt)){ #i=17
  taxon = ttt[i]#'Castanea'
  train <- pts.vars90b |> mutate(pos= ifelse(Level2 %in% taxon, 1,0)) |> st_drop_geometry()
    
  train <- train |> group_by(pos) |> mutate(wts = 1/(length(pos)+1)) |> ungroup() |> mutate(wts = ifelse(origin %in% 'gbif', wts*gratio, wts)) 
  sample.fraction <- 2*n/nrow(train)
    #i=18;k=1
    set.seed(1)# for predictability

        

    
    formular.rf <- as.formula(paste(paste("pos",paste(paste(covs, collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    
    
    rf <- ranger(formular.rf,
                 case.weights = train$wts,
                 split.select.weights = covwts,
                 sample.fraction=sample.fraction,
                 data=train)
    
#generate prediction raster ----
rf.prediction <-  predict(covrast, rf, na.rm=T);  names(rf.prediction) <- taxon
plot(rf.prediction)
# pop <- rast(paste0('gis/NWT_Models_GBIF/','Pinus','.tif'))
# plot(pop)
writeRaster(rf.prediction, paste0('gis/NWT_Models_GBIF/',taxon,'.tif'), overwrite=T)
  }

Sys.time() - timeA


#load back in
fils <- list.files('gis/NWT_Models_GBIF')
fils <- fils[grepl('.tif$', fils)]
alltrees <-  rast(paste0('gis/NWT_Models_GBIF/',fils))
prob1 = max(alltrees)
DomSpp = which.lyr(alltrees == prob1)
plot(DomSpp)
writeRaster(DomSpp, paste0('gis/NWT_Comm/DomSpp_gbif2.tif'), overwrite=T)
data.frame(id = 1:length(names(alltrees)), name = names(alltrees))



# ohio <- read.delim("C:/GIS/Vegetation/doi_10_5061_dryad_wwpzgmsv8__v20241129/Flinn_et_al._JVS_forest_type.txt")
# ohio <- ohio |> vect( geom=c("POINT_X", "POINT_Y"),crs=crs('epsg:4326')) 
# st_write(st_as_sf(ohio), 'C:/GIS/Vegetation/doi_10_5061_dryad_wwpzgmsv8__v20241129/ohiotrees.shp')