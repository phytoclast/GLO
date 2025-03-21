library(vegan)

library(climatools)
library(sf)
library(terra)
library(dplyr)

library(ggplot2)
library(ranger)
library(vegan)
library(cluster)
library(goeveg)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
MLRA <- st_read('gis/MLRA_2017.shp')
pts <- readRDS('points/pts.geo1.RDS')
seed1 = 4392 #Mount Rainier for resampling of Map
seed2 = 2010 #Mount LeConte for NMDS
seed3 = 4207 #Mauna Kea for clustering
seed4 = c(1345, 1917, 4809, 5636, 5895) #Ben Nevis, Mount Washington, Mont Blanc, Citlaltépetl, Kilimanjaro for KFold Cross Validation and Random Forest

#Load Data ----
set.seed(seed1)
pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
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
cors <- readRDS('output/cors.RDS')

#make species raster stack ----
for(i in 1:length(ttt)){
  x <- rast(paste0('gis/finalrun_mean/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Species <- rast(mget(ttt))

#make and extract sample points
MLRA <- st_transform(MLRA, crs(Species))
MLRA <- MLRA |> mutate(LRU = case_when(LRU %in% c('98A1','98A3','98A4','98A5') ~ '98A',
                                       LRU %in% c('97A1','97A2') ~ '97A',
                                       LRU %in% c('97B1','97B2') ~ '97B',
                                       TRUE ~ LRU))

MLRAall <- subset(MLRA, MLRA %in% c("94A","94C","96","97","98","99","111C") | LRU %in% "111BA", select="LRU")
MLRA98 <- subset(MLRA, MLRA %in% c("97","98","99","111C")| LRU %in% "111BA", select="LRU")
MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C","96"), select="LRU")
MLRAselect <- MLRAall
plot(MLRAselect)
MLRA.vect <- vect(MLRAselect)
spts <- terra::spatSample(MLRA.vect, size = 5000)
spts.mlra <- st_as_sf(spts) |> st_drop_geometry()
spts.spp <- terra::extract(Species, spts)
spts.env <- terra::extract(vars90, spts)

files <- data.frame(filename = list.files('gis/clusters'))
files <- files |> mutate(taxon = str_split_fixed(filename, '\\.', 2)[,1], model = str_split_fixed(filename, '\\.', 3)[,2])
ttt = unique(files$taxon)

#reassemble mean models
for(i in 1:length(ttt)){
  x <- rast(paste0('gis/clusters_mean/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Clusters <- rast(mget(ttt))


spts.clusters <- terra::extract(Clusters, spts)

spts.es <- cbind(spts.spp, spts.clusters[,-1], spts.env[,-1], spts.mlra) 
spts.es <- spts.es |> subset(!is.na(pH50) & !is.na(Tsuga))


spts.es <- spts.es |> mutate(spptotal = apply(spts.es[,ttt], MARGIN=1, FUN='sum')) |> subset(spptotal > 0)


#ES site key ---

spts.es <- spts.es |> mutate(ES = case_when(rockdepth < 150 ~ 'Limestone',
                                            floodfrq > 0 ~ 'Floodplain',
                                            histic > 0.5 &  hydric >= 0.5 ~ 'Mucks',
                                            sand150 >= 80 & sand50 >= 70 | sand50 >80 ~ 'Sandy',
                                            TRUE ~'Loamy'))

spts.es <- spts.es |> mutate(ES = case_when(ES %in%  'Floodplain' & hydric >= 0.5 ~ 'Wet Floodplain',
                                            ES %in%  'Floodplain' ~ 'Moist Floodplain',
                                            ES %in% 'Mucks' & pH50 >= 5 ~ "Euic Muck",
                                            ES %in% 'Mucks'  ~ "Acidic Peat",
                                            ES %in% 'Sandy' &  watertable > 100 ~ 'Dry Sand',
                                            ES %in% 'Sandy' &  hydric < 0.5 ~ 'Moist Sand',
                                            ES %in% 'Sandy' ~ 'Wet Sand',
                                            ES %in% 'Loamy' &  watertable > 100 ~ 'Dry Loam',
                                            ES %in% 'Loamy' &  hydric < 0.5 ~ 'Moist Loam',
                                            ES %in% 'Loamy' ~ 'Wet Loam',
                                            TRUE ~ ES))

spts.es <- spts.es |> mutate(ES = case_when(ES %in%  'Dry Sand' & LRU %in% '94AA' & Bhs >= 0.5 & spodic >= 0.5 ~ 'Rich Sand',
                                            ES %in%  'Dry Sand' & LRU %in% '94AA'  ~ 'Dry Sand',
                                            ES %in%  'Dry Sand' & (spodic >= 0.5 | pH50 >= 5 | carbdepth < 100) ~ 'Rich Sand',
                                            ES %in%  'Dry Sand'  ~ 'Dry Sand',
                                            ES %in%  'Moist Sand'  & (spodic >= 0.5 | pH50 <= 5) ~ 'Acidic Moist Sand',
                                            ES %in%  'Moist Sand'   ~ 'Euic Moist Sand',
                                            ES %in%  'Wet Sand'  & (spodic >= 0.5 | pH50 <= 5.5) ~ 'Acidic Wet Sand',
                                            ES %in%  'Wet Sand'   ~ 'Euic Wet Sand',
                                            TRUE ~ ES))
spts.es <- spts.es |> mutate(ES = ifelse(tan(slope)*100 >= 15 &  watertable > 100 & floodfrq <= 0 &  hydric < 0.5, paste(ES, "Slopes"), paste(ES, "Plains")))

spts.es <- spts.es |> mutate(ESD = paste(LRU, ES))

spts.es <- spts.es |> mutate(ESD = ifelse(rockdepth < 150, 'MLRA 94C Limestone Plains', ESD))


essummary <- spts.es |> group_by(ESD) |> summarise(across(.fns=mean))
write.csv(essummary, 'essummaryclusters.csv', row.names = F)
