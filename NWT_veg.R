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
#Core area
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

#Western New York ----
spref <- crs('+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs')
data1 <- read.csv('data/3569064/lot_line_section_with_predictors.csv') |> mutate(fil = 'WNY_lot')
data2 <- read.csv('data/3569064/township_bt_post_with_predictors.csv') |> mutate(fil = 'WNY_post')
data3 <- read.csv('data/3569064/township_line_section_with_predictors.csv')|> mutate(fil = 'WNY_line')
data.all <- rbind(data1[,colnames(data2)],data2,data3)
taxa <- c("fag.gra","ace.sac","til.ame","tsu.can","ulm.ame","fra.ame","mag.acu","cas.den","pin.str","fra.nig","que.vel","que.alb","jug.cin","car.spp","lir.tul","bet.all","pru.ser","aln.inc","ace.rub","pla.occ","abi.bal","ost.vir","jug.nig")
genera <- c("Fagus","Acer","Tilia","Tsuga","Ulmus","Fraxinus","Magnolia","Castanea","Pinus","Fraxinus","Quercus","Quercus","Juglans","Carya","Liriodendron","Betula","Prunus","Alnus","Acer","Platanus","Abies","Ostrya","Juglans")
spp <- c("Fagus","Acer","Tilia","Tsuga","Ulmus","Fraxinus americana","Magnolia","Castanea","Pinus strobus","Fraxinus nigra","Quercus velutina","Quercus alba","Juglans cinerea","Carya","Liriodendron","Betula alleghaniensis","Prunus","Alnus","Acer","Platanus","Abies","Ostrya","Juglans nigra")
for(i in 1:length(taxa)){#i=1
  taxon=taxa[i]
  thesecols <- c('fil','x_coord','y_coord', taxon)
  df0 <- data.all[,thesecols]
  colnames(df0)<- c("dataset","x","y","tname")
  df0 <- subset(df0, tname %in% 1)
  df0$tname <- taxon
  df0$Level2 <- genera[i]
  df0$Species <- spp[i]
  if(i==1){pts.wny <- df0}else{pts.wny <- rbind(pts.wny,df0)}
}
covswny <-  rast("newcovariates/WNYcovr.tif")
pts.wny <- st_as_sf(pts.wny, coords=c('x', 'y'), crs=spref) |> st_transform(crs=st_crs(covswny))
covswny0 <- extract(covswny, vect(pts.wny))
pts.wny <- cbind(pts.wny,covswny0)

# data.all <- st_as_sf(data.all, coords=c('x_coord', 'y_coord'), crs=spref)
# st_write(data.all, 'points/wny.shp', append = FALSE)

#Monongahela ----
covsmon <-  rast("newcovariates/Monongahela_covr.tif")
dat1 <- st_read('data/RDS-2014-0022_Data/Data/corners_1752_1899_revised_052012.shp')|> st_transform(crs=st_crs(covsmon))
spp1 <- read.csv('data/RDS-2014-0022_Data/Data/species_codes_witness_trees_Monongahela.csv')
spp1$Genus <- vegnasis::extractTaxon(spp1$scientific.name, 'genus')
spp1$Genus <- ifelse(grepl("Maclura", spp1$scientific.name), 'Maclura', spp1$Genus) 
spp1$Genus <- ifelse(grepl("spruce-pine", spp1$common.name), 'TsugaPicea', spp1$Genus)#to resolve later
for(i in 1:6){
  dat0 <- dat1[c(1:3,3+i)]
  colnames(dat0) <- c("id","DATE_","treeseq", "species.code", "geometry")
  if(i==1){pts.mon <- dat0}else{pts.mon <- rbind(pts.mon,dat0)}
}
pts.mon <- pts.mon |> left_join(spp1) |> subset(!is.na(Genus))
pts.mon$dataset <- 'Monongahela' 
covsmon0 <- extract(covsmon, vect(pts.mon))
pts.mon <- cbind(pts.mon,covsmon0)



#continued ....----
# covrast270 <- aggregate(covrast, 3)
covrast270 <- rast("newcovariates/covrast270.tif")
# writeRaster(covrast270, "newcovariates/covrast270.tif")
# 
# pts.vars90xy <- extract(xyrast, vect(pts.vars90)) 
# saveRDS(pts.vars90xy,'pts.vars.xyrast.RDS')
pts.vars90xy <- readRDS('pts.vars.xyrast.RDS')
pts.vars90 <- pts.vars90 |> cbind(pts.vars90xy)
# #gbifdata ---- 2026-11-10 stop using gbif data
# gbifdata <- readRDS('data/gbifdata.RDS')
# # newcovr <- crop(covrast270, c(1100000, 1403971, 1900000, 2300000))
# # newcovr <- crop(covrast270, c(500000, 600000, 2750000, 2850000))
# gbifdata$Level2 <- extractTaxon(gbifdata$name, 'genus')
# gbifdata$Level2 <- ifelse(gbifdata$Level2 %in% c('Ostrya', 'Carpinus'),'OstryaCarpinus',gbifdata$Level2)
# gbifdata <- gbifdata |> st_as_sf(coords=c('decimalLongitude', 'decimalLatitude'), crs=crs('epsg:4326')) |> st_transform(crs(covrast))
# gbifcovs <- extract(covrast, vect(gbifdata)) 
# gbifdata <- gbifdata |> cbind(gbifcovs)
# gbifdata <- gbifdata |> subset(!is.na(sand150) & !is.na(elev) & !is.na(popen) & !is.na(lat1))
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


# #merge with gbif data; 2026-11-10: GBIF was found to create too much Populus and Prunus, even when adjusting the weights of the data sets to match the species ratios of the original dataset.
# pts.vars90a$origin <- 'glo'
# gbifdata$origin <- 'gbif'
# newcols <- intersect(colnames(pts.vars90a), colnames(gbifdata))
# 
# #reweigh gbif species due to increase in abundance in disturbance species
# genwts1 <- pts.vars90a|>st_drop_geometry() |> mutate(tots = length(Level2)) |> group_by(Level2, tots) |> summarise(gentots = length(Level2)) |> ungroup() |> mutate(gentots = round(100*gentots/tots,3), tots=NULL) 
# genwts2 <- gbifdata |>st_drop_geometry() |> mutate(tots = length(Level2)) |> group_by(Level2, tots) |> summarise(gentots2 = length(Level2)) |> ungroup() |> mutate(gentots2 = round(100*gentots2/tots,3), tots=NULL) 
# genwts1 <- genwts1 |> left_join(genwts2) |> mutate(gratio = gentots/gentots2)
# 
# pts.vars90b <- st_drop_geometry(pts.vars90a[,newcols]) |> rbind(st_drop_geometry(gbifdata[,newcols]))
# pts.vars90b <- pts.vars90b |> left_join(genwts1[,c('Level2', 'gratio')]) 

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